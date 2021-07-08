## Summary:
Calculate the nonoverlapping mapped length and the number of mismatches per read pair and per mate (if tag NM available)

## Usage:
    perl sam_mappedLength_mismatches.pl --inFile <input.bam|sam> [--outFile <output.txt>] [--endType PE] [--debug info] [--version] [--help] [--ncores 1]

## Options: 

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

## Notes:
* Column 7) and 8) will be available only if tag NM (not nM) is set;
* If SE, column 2) will be identical to 3) and 5), column 6) and 8) will be blank, and --nameSorted, --ncores will be dismissed because we don't need to pair two mate before computing the stats;
* If PE but R1 and R2 are mapped to different chromosome (namely discordant mapping), 2) will be reported as the sum of each mate but 3) will be -1.

## Examples:

    $ perl sam_mappedLength_mismatches.pl -i test/star_input.bam 2>/dev/null
    
    ReadID  RawLength       NonoverlapLength        Mismatches      R1_MappedLength R2_MappedLength R1_Mismatches   R2_Mismatches
    NB501328:197:HMK3KBGX7:1:11101:1498:13562       56      56      3  56   
    NB501328:197:HMK3KBGX7:1:11101:1744:6598        48      24      4  24   24
    NB501328:197:HMK3KBGX7:1:11101:2468:19747       107     107     0  54   53
    NB501328:197:HMK3KBGX7:1:11101:2835:13039       50      25      6  25   25
    

    $ perl sam_mappedLength_mismatches.pl -i test/bt2_input.bam 2>/dev/null
    
    ReadID  RawLength       NonoverlapLength        Mismatches      R1_MappedLength R2_MappedLength R1_Mismatches   R2_Mismatches
    NB501328:251:HTM5TBGXC:4:23612:26830:1974       112     112        55   57      7       2
    NB501328:251:HTM5TBGXC:4:23612:26831:7200       119     119        53   66      0       0
    NB501328:251:HTM5TBGXC:4:23612:26837:13055      111     111        53   58      1       0
    NB501328:251:HTM5TBGXC:4:23612:26842:11567      52      26         26   26      2       0
    NB501328:251:HTM5TBGXC:4:23612:26848:9630       112     112        55   57      0       1
