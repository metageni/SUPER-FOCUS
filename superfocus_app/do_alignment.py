# !/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import csv

from pathlib import Path
from collections import defaultdict


def normalise_counts(data):
    """Normalise query counts based on number of hits.

    Args:
        data (collections.defaultdict): Data to be normalised

    Returns:
        dict: Normalised data

    """
    sum_hits = sum(data.values())

    return {level: data[level] / sum_hits for level in data}


def update_results(results, sample_index, data, normalise, number_samples):
    """Update final result dict.

    Args:
        results (collections.defaultdict): Results that is updated for each sample
        sample_index (int): Sample index in the result
        data (dict): Data to be added into the final results
        normalise (int): 0/1 Normalise or not the data
        number_samples (int): Number of samples in the analyses

    Returns:
        collections.defaultdict: Updated results

    """
    # if any hit was recorded
    if data:
        # if len(data) not > 1, no need to normalise - it is normalised already -- ha!
        if normalise and len(data) > 1:
            data = normalise_counts(data)

        # update final result
        for level in data:
            if level not in results:
                results[level] = [0] * number_samples
            results[level][sample_index] += data[level]

    return results


def align_reads(query, output_dir, aligner, database, evalue, threads, fast_mode, WORK_DIRECTORY, amino_acid):
    """Align FAST(A/Q) file to database.

    Args:
        query (PosixPath): Path to query in FAST(A/Q) file
        output_dir (PosixPath): Path to alignment output
        aligner (str): Aligner name
        database (str): Database name
        evalue (str): E-value
        threads (str): Number of threads (default = all)
        fast_mode (str): Fast or sensitive mode (default = 1)
        WORK_DIRECTORY (str): Path to directory where works happens
        amino_acid (str): 0 input nucleotides, 1 amino acid

    Returns:
        str: Path to alignment that was written

    """
    # prepare variables
    output_name = "{}/{}_alignments".format(output_dir, query.parts[-1])

    if aligner == "diamond":
        temp_folder = Path("{}/db/tmp/".format(WORK_DIRECTORY))

        if not temp_folder.exists():
            temp_folder.mkdir(parents=True, mode=511)

        mode_diamond = "" if fast_mode == "1" else "--sensitive"
        database_diamond = "{}/db/static/diamond/{}.db".format(WORK_DIRECTORY, database)

        # align
        os.system("diamond blastx -t {} -d {} -q {} -a {} -p {} -e {} {}".format(temp_folder, database_diamond, query,
                                                                                 output_name, threads, evalue,
                                                                                 mode_diamond))
        # dump
        os.system("diamond view -a {}.daa -o {}.m8".format(output_name, output_name))
        # delete binary file
        os.system("rm {}/*.daa".format(output_dir))
        # add aligner extension to output
        output_name = "{}.m8".format(output_name)

    elif aligner == "rapsearch":
        mode_rapsearch = "T" if fast_mode == "1" else "F"
        database_rapsearch = "{}/db/static/rapsearch2/{}.db".format(WORK_DIRECTORY, database)

        os.system('rapsearch -a {} -q {} -d {} -o {} -v 250 -z {} -e {} -b 0 -s f'.format(mode_rapsearch, query,
                                                                                          database_rapsearch,
                                                                                          output_name, threads, evalue))
    elif aligner == "blast":
        database_blast = "{}/db/static/blast/{}.db".format(WORK_DIRECTORY, database)
        blast_mode = 'blastp' if amino_acid == '1' else 'blastx'

        os.system('{} -db {} -query {} -out {} -outfmt 6 -evalue {} -max_target_seqs 250 -num_threads {}'.format(
            blast_mode, database_blast, query, output_name, evalue, threads))

    return '{}.m8'.format(output_name) if aligner == 'rapsearch' else output_name


def parse_alignments(alignment, results, normalise, number_samples, sample_index,
                     minimum_identity, minimum_alignment, subsystems_translation, aligner, binning_reads, query_name):
    """Parses alignment.

    Args:
        alignment (PosixPath): Path to sample alignment
        results (collections.defaultdict): Results that is updated for each sample
        normalise (int): 0/1 Normalise or not the data
        number_samples (int): Number of samples in the analyses
        sample_index (int): Sample index in the result
        minimum_identity (int): Minimum identity to consider a hit
        minimum_alignment (int) Minimum alignment (bp) to be consider a hit
        subsystems_translation (dict): subsystems translation lookup table
        aligner (str): aligner name
        binning_reads (collections.defaultdict): reads binning results
        query_name (str): fasta/q file name used in the alignment

    Returns:
        collections.defaultdict: Updated results

    """
    previous_read_name = None
    best_evalue = None

    # test for an empty file
    if os.stat(alignment).st_size == 0:
        return defaultdict(int)

    with open(alignment) as alignment_file:
        alignment_reader = csv.reader(alignment_file, delimiter='\t')
        if aligner == "rapsearch":
            # skip header
            [next(alignment_reader, None) for _ in range(5)]

        temp_results = defaultdict(int)
        for row in alignment_reader:
            # extract need info from hit
            current_read_name = row[0]
            temp_mi = float(row[2])
            temp_ml = float(row[3])
            current_evalue = row[10]

            # create needed info
            current_hit = row[1].split("__")
            current_subsystem_id = current_hit[1]  # Subsystem PK on SF database
            current_function_name = current_hit[-1].replace("\n", "").replace("\r", "")
            aggregate_levels = subsystems_translation[current_subsystem_id] + "\t" + current_function_name

            # found a different read
            if current_read_name != previous_read_name:
                if previous_read_name:
                    update_results(results, sample_index, temp_results, normalise, number_samples)
                temp_results = defaultdict(int)
                best_evalue = current_evalue

            # check if alignment wanted
            if temp_mi >= minimum_identity and temp_ml >= minimum_alignment and current_evalue == best_evalue:
                temp_results[aggregate_levels] = 1
                binning_reads[query_name][current_read_name].append([temp_mi, temp_ml, current_evalue,
                                                                     aggregate_levels])

            previous_read_name = current_read_name

        # last group of reads
        update_results(results, sample_index, temp_results, normalise, number_samples)

    return results, binning_reads
