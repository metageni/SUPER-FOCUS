#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import os
import sys
import logging

from pathlib import Path
from shutil import which

from superfocus_app import version

LOGGER_FORMAT = '[%(asctime)s - %(levelname)s] %(message)s'
logging.basicConfig(format=LOGGER_FORMAT, level=logging.INFO)
LOGGER = logging.getLogger(__name__)

WORK_DIRECTORY = str(Path(__file__).parents[0])


def check_aligners():
    """Check if aligners are installed on the system path.

    Returns:
        list: Paths for aligners if exists.

    """
    return [
        info
        for info in [which(x) for x in ['prerapsearch', 'diamond', 'makeblastdb']]
        if info != None
    ]


def format_database(aligners, target_files):
    """Format database for target aligner(s) and folder with files.

    Args:
        aligners (list): Aligners to have the database formatted.
        target_files (str): Download database.


    """
    LOGGER.info('  Joining files')
    [
        os.system(
            'cat {}/{}/*.faa > {}.fasta'.format(target_files, cluster, cluster)
        )
        for cluster in os.listdir(target_files) if cluster.endswith("_clusters")
    ]

    # Format database
    LOGGER.info('  Formatting Database')
    cluster_identities = ['100', '98', '95', '90']
    for dbname in cluster_identities:
        for aligner in aligners:
            if aligner == 'prerapsearch':
                LOGGER.info('  RAPSearch2: DB_{}'.format(dbname))
                os.system(
                    'prerapsearch -d {}_clusters.fasta '
                    ' -n {}/db/static/rapsearch2/{}_clusters.db'.format(dbname, WORK_DIRECTORY, dbname)
                )
            elif aligner == 'diamond':
                LOGGER.info('  DIAMOND: DB_{}'.format(dbname))
                os.system(
                    'diamond makedb --in  {}_clusters.fasta '
                    '--db {}/db/static/diamond/{}_clusters.db'.
                        format(dbname, WORK_DIRECTORY, dbname)
                )
            elif aligner == 'makeblastdb':
                LOGGER.info('  BLAST: DB_{}'.format(dbname))
                os.system(
                    'makeblastdb -in {}_clusters.fasta '
                    '-out {}/db/static/blast/{}.db -title {}.db -dbtype prot'.
                        format(dbname, WORK_DIRECTORY, dbname, dbname)
                )
        # remove joint file
        os.remove("{}_clusters.fasta".format(dbname))


def parse_args():
    """Parse args entered by the user.

    Returns:
        argparse.Namespace: Parsed arguments.

    """
    parser = argparse.ArgumentParser(description="SUPER-FOCUS: A tool for agile functional analysis of shotgun "
                                                 "metagenomic data",
                                     epilog="superfocus_downloadDB -a diamond,rapsearch,blast -i clusters/")
    # basic parameters
    parser.add_argument("-a", "--aligner", help="Aligner name separed by ',' if more than one", required=True)
    parser.add_argument("-i", "--input", help="Target input files to be formatted into the database", required=True)
    parser.add_argument('-v', '--version', action='version', version='superfocus_downloadDB version {}'.format(
        version))

    return parser.parse_args()


def main():
    args = parse_args()

    valid_aligners = {'rapsearch', 'diamond', 'blast'}
    aligner_db_creators = {os.path.basename(x) for x in check_aligners() if x != 'None'}
    if not aligner_db_creators:
        LOGGER.critical('  None of the required aligners are installed {}'.format(list(valid_aligners)))
        sys.exit(1)

    # Parse which aligner(s) the user wants to format the database to
    requested_aligners = {aligner.lower() for aligner in args.aligner.split(",")}
    for aligner in requested_aligners:
        if aligner not in valid_aligners:
            LOGGER.critical('  None of the required aligners are installed {}'.format(list(valid_aligners)))
            sys.exit(1)

    aligners = []
    if 'rapsearch' in requested_aligners and 'prerapsearch' in aligner_db_creators:
        aligners.append("prerapsearch")
    if 'diamond' in requested_aligners and 'diamond' in aligner_db_creators:
        aligners.append("diamond")
    if 'blast' in requested_aligners and 'makeblastdb' in aligner_db_creators:
        aligners.append("makeblastdb")

    os.system('rm {}/db/clusters/ -r 2> /dev/null'.format(WORK_DIRECTORY))  # delete folder if exists

    if aligners:
        format_database(aligners, args.input)
        LOGGER.info('  Done! Now you can run superfocus')
    else:
        LOGGER.critical('  No valid aligner. We cannot move on!')


if __name__ == "__main__":
    main()
