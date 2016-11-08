#SUPER-FOCUS plots
Written by Daniel Cuevas (https://github.com/dacuevas)

_R_ scripts to create bar charts and tables of SUPER-FOCUS output.

###DEPENDENCIES
Required _R_ packages.

1. getopt
2. ggplot2
3. reshape2
4. plyr
5. gridExtra

###SAMPLE PLOTS
![alldata](https://github.com/Adrian-Cantu/cf_pipeline/blob/master/scripts/superfocus/sample/all_top_functions.png "All Top Functions")

![samplebar](https://github.com/Adrian-Cantu/cf_pipeline/blob/master/scripts/superfocus/sample/sample1.fasta_top_functions.png "Sample Top Functions Barchart")

![sampletab](https://github.com/Adrian-Cantu/cf_pipeline/blob/master/scripts/superfocus/sample/sample1.fasta_top_functions_table.png "Sample Top Functions Table")

###USAGE
```
usage: superfocus_plots.sh -d SF_dir [Options]

Required
   -d [SF_dir]             : SUPER-FOCUS directory of files

Optional
   -vir                    : Create virulence-specific plots
   -h, -?, --help          : This help message
   -v                      : Verbose output
```
