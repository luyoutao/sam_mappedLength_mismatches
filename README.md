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
