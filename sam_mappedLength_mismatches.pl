#' This Source Code Form is subject to the terms of the Mozilla Public
#' License, v. 2.0. If a copy of the MPL was not distributed with this
#' file, You can obtain one at http://mozilla.org/MPL/2.0/.
#'
#' Youtao Lu@Kim Lab, 2016-2020

use warnings;
use strict;
use Getopt::Long;
use Log::Log4perl qw(get_logger :levels);
use Log::Log4perl::Layout::PatternLayout;

our $VERSION = '0.3';
our $LOGGER  = get_logger(__PACKAGE__);
my $debug      = "info";
my $nameSorted = 0;
my $help       = 0;
my $version    = 0;
my $ncores     = 1;
my $endType    = "PE";
my ( $inFile,  $inFh );
my ( $outFile, $outFh );

sub usage {
    print STDERR <<DOC;
Summary:
    Calculate the nonoverlapping mapped length and the number of mismatches per read pair and per mate (if tag NM available)

Usage:
    perl sam_mappedLength_mismatches.pl --inFile <input.bam|sam> [--outFile <output.txt>] [--endType PE] [--debug info] [--version] [--help] [--ncores 1]

Options: 
    --inFile, -i    BAM or SAM input. If multiple mapping, only computes for the primary alignment;
    --outFile, -o   Optional, if omitted, it will output to STDOUT. If the filename ends with ".gz", the output will be GZ compressed. The output table consists of eight columns for each read (SE) or read pair (PE):
        1) read ID
        2) raw mapped length per read pair
        3) nonoverlapping mapped length per read pair
        4) #. of mismatches per read pair
        5) mapped length of R1
        6) mapped length of R2
        7) #. of mismatches in R1
        8) #. of mismatches in R2
    --endType       PE or SE (default: PE)
    --nameSorted    whether the input is sorted by query name (default no). If not, the input will be name sorted and the temporary file has a '.nameSorted.bam' suffix;
    --ncores        how many cores for sorting;
    --debug         debug level, choose from fatal, error, warn, info, debug, trace

Notes:
    Column 7) and 8) will be available only if tag NM (not nM) is set;
    If SE, column 2) will be identical to 3) and 5), column 6) and 8) will be blank, and --nameSorted, --ncores will be dismissed because we don't need to pair two mate before computing the stats.
    If PE but R1 and R2 are mapped to different chromosome (namely discordant mapping), 2) will be reported as the sum of each mate but 3) will be -1.
DOC
}

GetOptions(
    "inFile|i=s"  => \$inFile,
    "outFile|o:s" => \$outFile,
    "nameSorted"  => \$nameSorted,
    "endType:s"   => \$endType,
    "ncores:i"    => \$ncores,
    "version|v"   => \$version,
    "debug=s"     => \$debug,
    "help|h"      => \$help,
) or print STDERR "Wrong arguments!\n" && &usage() && exit(-1);

&usage() && exit(-1) if $help or !defined($inFile);
if ( $endType ne "PE" && $endType ne "SE" ) {
    die("--endType can only be PE or SE!\n");
}
if ( !defined($inFile) || !-e $inFile ) {
    die("--inFile doesn't exist!\n");
}
if ( $inFile !~ /\.(sam|bam)$/g ) {
    die("--inFile doesn't look like SAM or BAM!\n");
}
if ( !defined($outFile) ) {
    $outFile = '';
}

if ( $debug eq "fatal" ) {
    $LOGGER->level($FATAL);
}
elsif ( $debug eq "error" ) {
    $LOGGER->level($ERROR);
}
elsif ( $debug eq "warn" ) {
    $LOGGER->level($WARN);
}
elsif ( $debug eq "info" ) {
    $LOGGER->level($INFO);
}
elsif ( $debug eq "debug" ) {
    $LOGGER->level($DEBUG);
}
elsif ( $debug eq "trace" ) {
    $LOGGER->level($TRACE);
}
my $appender = Log::Log4perl::Appender->new("Log::Log4perl::Appender::Screen");
my $layout   = Log::Log4perl::Layout::PatternLayout->new(
    "[%d{yyyy-MM-dd HH:mm:ss.SSS Z} %p] %m");
$appender->layout($layout);
$LOGGER->add_appender($appender);

$LOGGER->info(
"{ VERSION = $VERSION; inFile = $inFile; outFile = $outFile; endType = $endType; nameSorted = $nameSorted; ncores = $ncores; help = $help; debug = $debug; version = $version }\n"
);

{

    package Read;
    my $endType = "PE";

    sub import {
        my $class = shift;
        $endType = shift;
    }

    sub new {
        my $class = shift;
        my ( $readID, $flag, $chr, $pos, $cigar, @tags ) = @_;
        $LOGGER->trace(
            "(Read->new) $readID, $flag, $chr, $pos, $cigar, @tags\n");
        my $self = bless {
            readID       => $readID,
            flag         => $flag,
            chr          => $chr,
            pos          => $pos,
            whichInPair  => undef,
            isRev        => undef,
            cigar        => $cigar,
            tags         => \@tags,
            mappedLength => undef,
            mismatches   => undef,
        }, $class;
        return $self;
    }

    sub parse_flag {
        my $self = shift;
        $self->{whichInPair} =
          $endType eq "PE" ? ( $self->{flag} & 0x40 ? "R1" : "R2" ) : "R1";
        $self->{isRev} = $self->{flag} & 0x10 ? 1 : 0;
    }

    sub parse_cigar {
        my $self         = shift;
        my $mappedLength = 0;
        while ( $self->{cigar} =~ s/(\d+)[MX=]// ) {
            $mappedLength += $1;
        }
        $self->{mappedLength} = $mappedLength;
    }

    sub parse_mismatches {
        my $self = shift;
        my $NM   = "";
        my @t    = grep { /^NM/ } @{ $self->{tags} };
        $NM = substr( $t[0], 5 ) if ( $#t == 0 );
        $self->{mismatches} = $NM;
    }
}

{

    package ReadPair;
    my $endType = "PE";

    sub import {
        my $class = shift;
        $endType = shift;
    }

    sub new {
        my $class = shift;
        my $pair  = shift;
        my $self  = bless {
            R1               => $pair->{R1},
            R2               => $pair->{R2},
            rawLength        => undef,
            nonoverlapLength => undef,
            mismatches       => undef,
        }, $class;
    }

    sub is_discordant {
        my $self = shift;
        if ( $self->{R1}->{chr} ne $self->{R2}->{chr} ) {
            return 1;
        }
        else {
            return 0;
        }
    }

    sub parse_mappedLength {
        my $self = shift;
        if ( $endType eq "PE" ) {
            if ( defined( $self->{R1} ) && defined( $self->{R2} ) ) {
                $self->{R1}->parse_cigar();
                $self->{R2}->parse_cigar();
                $self->{rawLength} =
                  $self->{R1}->{mappedLength} + $self->{R2}->{mappedLength};
                if ( $self->is_discordant() ) {
                    $self->{nonoverlapLength} = -1;
                    return;
                }
                my $overlap =
                    $self->{R1}->{isRev}
                  ? $self->{R2}->{pos} + $self->{R2}->{mappedLength} -
                  $self->{R1}->{pos}
                  : $self->{R1}->{pos} +
                  $self->{R1}->{mappedLength} -
                  $self->{R2}->{pos};
                $overlap = $overlap > 0 ? $overlap : 0;
                $self->{nonoverlapLength} = $self->{rawLength} - $overlap;
                return;
            }
            elsif ( defined( $self->{R1} ) ) {
                $self->{R1}->parse_cigar();
                $self->{rawLength} = $self->{nonoverlapLength} =
                  $self->{R1}->{mappedLength};
                return;
            }
            else {
                $self->{R2}->parse_cigar();
                $self->{rawLength} = $self->{nonoverlapLength} =
                  $self->{R2}->{mappedLength};
                return;
            }
        }
        else {
            $self->{rawLength} = $self->{nonoverlapLength} =
              $self->{R1}->{mappedLength};
            return;
        }
    }

    sub parse_mismatches {
        my $self = shift;
        my @t;
        if ( $endType eq "PE" ) {
            if ( defined( $self->{R1} ) ) {
                $self->{R1}->parse_mismatches();
                @t = grep { /^nM/ } @{ $self->{R1}->{tags} };
            }
            if ( defined( $self->{R2} ) ) {
                $self->{R2}->parse_mismatches();
                if ( !defined( $self->{R1} ) ) {
                    @t = grep { /^nM/ } @{ $self->{R2}->{tags} };
                }
            }
            $self->{mismatches} = $#t == 0 ? substr( $t[0], 5 ) : "";
            $LOGGER->trace(
"(ReadPair->parse_mismatches) \$self->{mismatches} = $self->{mismatches}\n"
            );
            return;
        }
        else {
            $self->{R1}->parse_mismatches();
            @t = grep { /^nM/ } @{ $self->{R1}->{"tags"} };
            $self->{mismatches} = $#t == 0 ? substr( $t[0], 5 ) : "";
            $LOGGER->trace(
"(ReadPair->parse_mismatches) \$self->{mismatches} = $self->{mismatches}\n"
            );
            return;
        }
    }

    sub output {
        my $self  = shift;
        my $outFh = shift;
        my $line  = (
            defined( $self->{R1} )
            ? $self->{R1}->{readID}
            : $self->{R2}->{readID} )
          . "\t"
          . $self->{rawLength} . "\t"
          . $self->{nonoverlapLength} . "\t"
          . $self->{mismatches} . "\t"
          . ( defined( $self->{R1} ) ? $self->{R1}->{mappedLength} : "" ) . "\t"
          . ( defined( $self->{R2} ) ? $self->{R2}->{mappedLength} : "" ) . "\t"
          . ( defined( $self->{R1} ) ? $self->{R1}->{mismatches}   : "" ) . "\t"
          . ( defined( $self->{R2} ) ? $self->{R2}->{mismatches} : "" ) . "\n";
        print $outFh $line;
    }
}

{

    package SamReader;
    use IO::Zlib;

    sub new {
        my $class = shift;
        my ( $inFile, $outFile, $endType, $nameSorted, $ncores ) = @_;
        Read->import($endType);        # Initialize class static
        ReadPair->import($endType);    # Initialize class static
        my $self = bless {
            inFile     => $inFile,
            outFile    => $outFile,
            inFh       => undef,
            outFh      => undef,
            endType    => $endType,
            nameSorted => $nameSorted,
            ncores     => $ncores,
        }, $class;
        return $self;
    }

    sub init_fh {
        my $self = shift;
        if ( $self->{endType} eq "PE" && !$self->{nameSorted} ) {
            my $tmpFile = $self->{inFile};
            $tmpFile =~ s/\.(bam|sam)$//i;
            my $suffix = $1;
            $tmpFile .= ".nameSorted" . ".$suffix";
            $LOGGER->warn(
"$self->{inFile} is not name sorted. Name sorting it and saving to $tmpFile...\n"
            );
            $LOGGER->warn("$tmpFile exists already! It will be overwritten.\n")
              if -e $tmpFile;
            my $exit = system(
"samtools sort -n -\@ $self->{ncores} -o $tmpFile $self->{inFile}"
            );
            $LOGGER->fatal("Error happended when sorting $self->{inFile}!")
              && die()
              if $exit;
            $self->{inFile} = $tmpFile;
        }
        $LOGGER->info("Opening file handle for $self->{inFile}...\n");
        open( $self->{inFh}, "samtools view -F 0x100 $self->{inFile} |" )
          or $LOGGER->fatal("Cannot open $self->{inFile} for read!\n") && die();
        if ( $self->{outFile} eq '' ) {
            $LOGGER->info("Opening file handle for STDOUT...\n");
            $self->{outFh} = *STDOUT;
        }
        else {
            $LOGGER->info("Opening file handle for $self->{outFile}...\n");
            if ( $self->{outFile} =~ /\.gz$/i ) {
                ( $self->{outFh} = IO::Zlib->new( $self->{outFile}, 'w' ) )
                  or $LOGGER->fatal("Cannot open $self->{outFile} for write!\n")
                  && die();
            }
            else {
                open( $self->{outFh}, ">", $self->{outFile} )
                  or $LOGGER->fatal("Cannot open $self->{outFile} for write!\n")
                  && die();
            }
        }
    }

    sub fin_fh {
        my $self = shift;
        $LOGGER->info("Closing file handle for $self->{outFile}...\n");
        close $self->{outFh}
          or $LOGGER->warn("Cannot close $self->{outFile}!\n");
        $LOGGER->info("Closing file handle for $self->{inFile}...\n");
        close $self->{inFh} or $LOGGER->warn("Cannot close $self->{inFile}!\n");
    }

    sub iterate {
        my $self = shift;
        my ( $inFh, $outFh ) = ( $self->{inFh}, $self->{outFh} );
        my ( $readID, $flag, $chr, $pos, $cigar, @tags );
        my @F;
        print { $self->{outFh} }
"ReadID\tRawLength\tNonoverlapLength\tMismatches\tR1_MappedLength\tR2_MappedLength\tR1_Mismatches\tR2_Mismatches\n";
        if ( $self->{endType} eq "PE" ) {
            my ( $read, $prev_read, $readPair );
            my $whichInPair;
            while ( readline($inFh) ) {
                chomp;
                @F      = split( "\t", $_ );
                $readID = $F[0];
                $flag   = $F[1];
                $chr    = $F[2];
                $pos    = $F[3];
                $cigar  = $F[5];
                @tags   = @F[ 11 .. $#F ];
                $read = Read->new( $readID, $flag, $chr, $pos, $cigar, @tags );

                if ( !defined( $prev_read->{readID} ) ) {
                    $prev_read = $read;
                    next;
                }
                if ( $prev_read->{readID} eq $read->{readID} ) {
                    $prev_read->parse_flag();
                    $read->parse_flag();
                    $LOGGER->debug(
"(SamReader->iterate) \$prev_read->{whichInPair} = $prev_read->{whichInPair}\n"
                    );
                    $LOGGER->debug(
"(SamReader->iterate) \$read->{whichInPair} = $read->{whichInPair}\n"
                    );
                    $readPair = ReadPair->new(
                        {
                            $prev_read->{whichInPair} => $prev_read,
                            $read->{whichInPair}      => $read,
                        }
                    );
                    $readPair->parse_mappedLength();
                    $readPair->parse_mismatches();
                    $readPair->output( $self->{outFh} );
                    $prev_read = undef;
                }
                else {
                    $prev_read->parse_flag();
                    $LOGGER->debug(
"(SamReader->iterate) \$prev_read->{whichInPair} = $prev_read->{whichInPair}\n"
                    );
                    $readPair = ReadPair->new(
                        { $prev_read->{whichInPair} => $prev_read } );
                    $readPair->parse_mappedLength();
                    $readPair->parse_mismatches();
                    $readPair->output( $self->{outFh} );
                    $prev_read = $read;
                }
            }
            if ( defined( $prev_read->{readID} ) ) {
                $prev_read->parse_flag();
                $LOGGER->debug(
"(SamReader->iterate) \$prev_read->{whichInPair} = $prev_read->{whichInPair}\n"
                );
                $readPair =
                  ReadPair->new( { $prev_read->{whichInPair} => $prev_read } );
                $readPair->parse_mappedLength();
                $readPair->parse_mismatches();
                $readPair->output( $self->{outFh} );
            }
        }
        else {
            my ( $inFh, $outFh ) = ( $self->{inFh}, $self->{outFh} );
            my ( $readID, $flag, $chr, $pos, $cigar, @tags );
            my ( $read, $readPair );
            while ( readline($inFh) ) {
                chomp;
                @F      = split( "\t", $_ );
                $readID = $F[0];
                $flag   = $F[1];
                $chr    = $F[2];
                $pos    = $F[3];
                $cigar  = $F[5];
                @tags   = @F[ 11 .. $#F ];
                $read = Read->new( $readID, $flag, $chr, $pos, $cigar, @tags );
                $readPair = ReadPair->new( { R1 => $read } );
                $readPair->parse_mappedLength();
                $readPair->parse_mismatches();
                $readPair->output( $self->{outFh} );
            }
        }
    }
}

my $samReader =
  SamReader->new( $inFile, $outFile, $endType, $nameSorted, $ncores );
$samReader->init_fh();
$samReader->iterate();
$samReader->fin_fh();
$LOGGER->info("All done.\n");
