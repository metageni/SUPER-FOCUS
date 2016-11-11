#!/bin/bash
# superfocus_plots.sh
# Create figures for SUPER-FOCUS output
# GitHub: https://github.com/metageni/SUPER-FOCUS/
# Ref: Silva GGZ, Green K., B. E. Dutilh, and R. A. Edwards: SUPER-FOCUS: A tool for agile functional analysis of shotgun metagenomic data.
#      (Bioinformatics. 2015 Oct 9. pii: btv584.)
#
# Author: Daniel A Cuevas (dcuevas08.at.gmail.com)
# Created on 08 Nov 2016
# Updated on 11 Nov 2016

VERSION="0.1"

usage() {
    echo "$1
usage: $scriptname -d SF_dir [Options]

Required
   -d [SF_dir]             : SUPER-FOCUS directory of files

Optional
   --skip [NUM]            : Number of non-blank lines to skip before columns [Default: 2]
   --vir                   : Create virulence-specific plots
   -h, -?, --help          : This help message
   -v                      : Verbose output

Notes
   - This program specifically uses the SUPER-FOCUS output file:
      *_____results__all_levels_and_function.xls

   - R scripts specifically look for the columns:
      1) Subsystem Level 3 (case specific)
      2) SEED Function (case specific)
      3) Columns with the word 'abundance' (not case specific)

" >&2
}


getTime() {
    currtime=$(date "+[%F %H:%M:%S]")
}


timeStamp() {
    timestamp=$(date "+%Y%m%dT%H%M%S")
}


####################################################
#ARGUMENT PARSING
####################################################
scriptname=$(echo $0 | perl -ne '/\/?.*\/(.+)/; print $1;')
sfdir=""
vir=0
skip=2
verbose=0

# Set pipefail for catching errors in piped commands
set -o pipefail

while [[ $# != 0 ]]; do
    case $1 in
    -h|-\?|--help)
        usage $frmmax
        exit 2
        ;;
    -d)
        shift
        [[ ! $1 || $(printf "%s" "$1" | perl -ne 'm/(^-.$)/; print $1;') ]] && echo "Missing -d value" >&2 && usage && exit 2
        sfdir=$1
        ;;
    -v)
        verbose=1
        ;;
    --skip)
        shift
        [[ ! $1 || $(printf "%s" "$1" | perl -ne 'm/(^-.$)/; print $1;') ]] && echo "Missing --skip value" >&2 && usage && exit 2
        skip=$1
        ;;
    --vir)
        vir=1
        ;;
    *)
        echo "Unknown option $1" >&2
        usage
        exit 2
    esac
    shift
done

# Check if required variables are set
if [[ ! $sfdir ]]; then
    usage "Missing one or more required arguments."
    exit 2
fi


# Plotting functions
getTime && echo "${currtime}    *****Starting plotting scripts!*****"  >&1
if (( !$verbose )); then
    getTime && echo "${currtime}    Note: verbose flag was not set."  >&1
fi
cmd="Rscript superfocus_functions.R -d ${sfdir}/ -s ${skip}"
(( $verbose )) && getTime && echo "${currtime}    Executing $cmd"  >&1
eval $cmd  2>&1 | tee -a $log
[[ $? -ne 0 ]] && getTime && error "${currtime}    Fail on command: $cmd"

if (( $vir )); then
    # Plotting virulence functions
    cmd="Rscript superfocus_virulence.R -d ${sfdir}/ -s ${skip}"
    (( $verbose )) && getTime && echo "${currtime}    Executing $cmd"  >&1
    eval $cmd  2>&1 | tee -a $log
    [[ $? -ne 0 ]] && getTime && error "${currtime}    Fail on command: $cmd"
fi
getTime && echo "${currtime}    *****Completed!*****"  >&1
