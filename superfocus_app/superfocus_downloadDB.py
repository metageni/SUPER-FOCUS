#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SUPER-FOCUS: A tool for agile functional analysis of shotgun metagenomic data"""

import argparse
import os
import sys
import logging

from pathlib import Path
from shutil import which

from superfocus_app import version
import tempfile

LOGGER_FORMAT = '[%(asctime)s - %(levelname)s] %(message)s'
logging.basicConfig(format=LOGGER_FORMAT, level=logging.INFO)
LOGGER = logging.getLogger(__name__)


def check_aligners():
    """Check if aligners are installed on the system path.

    Returns:
        list: Paths for aligners if exists.

    """
    return [
        info
        for info in [which(x) for x in ['prerapsearch', 'diamond', 'makeblastdb', 'mmseqs']]
        if info != None
    ]


def format_database(aligners, target_files, cluster_identities, db_dir):
    """Format database(s) for target aligner(s) and folder with files.

    Args:
        aligners (list): Aligners to have the database formatted.
        target_files (str): Path to clusters database (fasta).
        cluster_identities (list): Cluster identities used when formatting database.
        db_dir (str): Optional path to user-specified working directory to build databases in.

    """
    if db_dir:
        workdir = Path(db_dir).resolve()
    else:
        workdir = Path(__file__).parents[0]
    db_static_dir = workdir / 'db/static/'

    LOGGER.info('Preparing database(s) in workdir: {}'.format(workdir))
    for dbname in cluster_identities:
        fasta_file = '{}_clusters.fasta'.format(dbname)
        cluster_dir = Path(target_files) / '{}_clusters'.format(dbname)
        if not cluster_dir.exists():
            LOGGER.error('Cluster directory "{}" does not exist'.format(cluster_dir))
            LOGGER.error('Did you specify the input directory correctly?')
            sys.exit(1)

        LOGGER.info('Concatenating {}/*.faa > {}'.format(cluster_dir, fasta_file))
        cat = os.system(
            'cat {}/*.faa > {}'.format(cluster_dir, fasta_file)
        )
        if cat != 0:
            LOGGER.error('Could not concatenate clusters at level {}'.format(dbname))
            sys.exit(1)
        
        LOGGER.info('Formatting database(s) ...')
        return_codes = []
        if 'prerapsearch' in aligners:
            LOGGER.info('RAPSearch2: DB_{}'.format(dbname))
            outdir = Path('{}/{}'.format(db_static_dir, "rapsearch2"))
            outdir.mkdir(parents=True, exist_ok=True)
            return_codes.append((
                "prerapsearch", 
                os.system('prerapsearch -d {} -n {}/{}_clusters.db'.
                    format(fasta_file, outdir, dbname)
            )))
        if 'diamond' in aligners:
            LOGGER.info('DIAMOND: DB_{}'.format(dbname))
            outdir = Path('{}/{}'.format(db_static_dir, "diamond"))
            outdir.mkdir(parents=True, exist_ok=True)
            return_codes.append((
                "diamond", 
                os.system('diamond makedb --in {} --db {}/{}_clusters.db'.
                    format(fasta_file, outdir, dbname)
            )))
        if 'makeblastdb' in aligners:
            LOGGER.info('BLAST: DB_{}'.format(dbname))
            outdir = Path('{}/{}'.format(db_static_dir, "blast"))
            outdir.mkdir(parents=True, exist_ok=True)
            return_codes.append((
                "blast", 
                os.system(
                    'makeblastdb -in {} -out {}/{}_clusters.db '
                    '-title {}_clusters.db -dbtype prot'.
                    format(fasta_file, outdir, dbname, dbname)
            )))
        if 'mmseqs' in aligners:
            LOGGER.info('MMSEQS: DB_{}'.format(dbname))
            outdir = Path('{}/{}'.format(db_static_dir, "mmseqs2"))
            outdir.mkdir(parents=True, exist_ok=True)
            return_codes.append((
                "mmseqs",
                os.system(
                    'mmseqs createdb {} {}/{}_clusters.db {}'.format(fasta_file, outdir, dbname)
            )))


        for aligner, return_code in return_codes:
            if not return_code == 0:
                LOGGER.error('Something went wrong with {}'.format(aligner))
                sys.exit(1)

        LOGGER.info('Removing temporary concatenated file "{}"'.format(fasta_file))
        os.remove(fasta_file)


def parse_args():
    """Parse args entered by the user.

    Returns:
        argparse.Namespace: Parsed arguments.

    """
    parser = argparse.ArgumentParser(
            description=__doc__, 
            epilog="superfocus_downloadDB -a diamond,rapsearch,blast,mmseqs -i clusters/")

    parser.add_argument("-a", "--aligner", required=True,
            help="Aligner name separed by ',' if more than one.")
    parser.add_argument("-c", "--clusters", default="90",
            help="DB types separed by ',' if more than one (e.g 90,95,98,100). Default: 90.")
    parser.add_argument("-i", "--input", required=True, 
            help="Target input files to be formatted into the database.")
    parser.add_argument("-d", "--db-dir", 
            help="Alternate database directory to store DB files in.")
    parser.add_argument('-v', '--version', 
            action='version', 
            version='superfocus_downloadDB version {}'.format(version))

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args()


def main():
    args = parse_args()

    valid_aligners = {'rapsearch', 'diamond', 'blast', 'mmseqs'}
    aligner_db_creators = {os.path.basename(x) for x in check_aligners() if x != 'None'}
    if not aligner_db_creators:
        LOGGER.critical('None of the required aligners are installed {}'.format(list(valid_aligners)))
        sys.exit(1)

    db_dir_exists = Path(args.db_dir).exists() and Path(args.db_dir).is_dir()
    if args.db_dir and not db_dir_exists:
        LOGGER.critical('Provided --db-dir "{}" does not exist'.format(args.db_dir))
        sys.exit(1)

    # Parse which aligner(s) the user wants to format the database to
    requested_aligners = {aligner.lower() for aligner in args.aligner.split(",")}
    for aligner in requested_aligners:
        if aligner not in valid_aligners:
            LOGGER.critical('None of the required aligners are installed {}'.format(list(valid_aligners)))
            sys.exit(1)

    aligners = []
    if 'rapsearch' in requested_aligners and 'prerapsearch' in aligner_db_creators:
        aligners.append("prerapsearch")
    if 'diamond' in requested_aligners and 'diamond' in aligner_db_creators:
        aligners.append("diamond")
    if 'blast' in requested_aligners and 'makeblastdb' in aligner_db_creators:
        aligners.append("makeblastdb")
    if 'mmseqs' in requested_aligners and 'mmseqs' in aligner_db_creators:
        aligners.append("mmseqs")
    clusters_target = args.clusters.split(",")

    for cluster in clusters_target:
        if cluster not in ["90","95", "98", "100"]:
            raise ValueError(
                '{} is not sa valid cluster. Please enter 90, 95, 98, or 100'.format(args.clusters))

    if aligners:
        LOGGER.info('Formatting databases for aligner(s): {} and cluster(s): {}'.format(args.aligner, args.clusters))
        format_database(aligners, args.input, clusters_target, args.db_dir)
        LOGGER.info('Done :)')
    else:
        LOGGER.critical('No valid aligner. We cannot move on!')
        sys.exit(1)


if __name__ == "__main__":
    main()
