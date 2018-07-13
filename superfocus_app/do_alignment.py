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

    return {level: data[level]/sum_hits for level in data}


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


def align_reads(query, output_dir, aligner, database, evalue, threads, fast_mode, WORK_DIRECTORY):
    """Align FAST(A/Q) file to database.

    Args:
        query_file (PosixPath): Path to query in FAST(A/Q) file
        output_dir (PosixPath): Path to alignment output
        aligner (str): Aligner choice
        evalue (str): E-value
        threads (str): Number of threads (default = all)
        fast_mode (str): Fast or sensitive mode (default = 1)
        WORK_DIRECTORY (str): Path to directory where works happens

    Returns:
        str: Path to alignment that was written

    """
    if aligner == "diamond":
        temp_folder = Path("{}/db/tmp/".format(WORK_DIRECTORY))

        if not temp_folder.exists():
            temp_folder.mkdir(parents=True, mode=511)

        database = "{}/db/static/diamond/{}.db".format(WORK_DIRECTORY, database)
        output_name = "{}/{}_alignments".format(output_dir, query.parts[-1])

        # prepare variables
        threads = "T" if threads == "all" else threads
        mode = "" if fast_mode == "1" else "--sensitive"

        # align
        os.system("diamond blastx -t {} -d {} -q {} -a {} -p {} -e {} {}".format(temp_folder, database, query,
                                                                                 output_name, threads, evalue, mode))
        # dump
        os.system("diamond view -a {}.daa -o {}.m8".format(output_name, output_name))
        # delete binary file
        os.system("rm {}/*.daa".format(output_dir))
        # add aligner extension to output
        output_name = "{}.m8".format(output_name)


    return output_name


def parse_alignments(alignment, results, normalise, number_samples, sample_index,
                     minimum_identity, minimum_alignment, subsystems_translation):
    """Parses alignment.

    Args:
        alignment (PosixPath): Path to sample alignment
        results (collections.defaultdict): Results that is updated for each sample
        normalise (int): 0/1 Normalise or not the data
        number_samples (int): Number of samples in the analyses
        sample_index (int): Sample index in the result
        minimum_identity (int): Minimum identity to consider a hit
        minimum_alignment (int) Minimum alignment (bp) to be consider a hit
        subsystems_translation (dict): Subsystems translation lookup table

    Returns:
        collections.defaultdict: Updated results

    """
    previous_read_name = None
    best_evalue = None

    with open(alignment) as alignment_file:
        alignment_reader = csv.reader(alignment_file, delimiter='\t')
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
                temp_results[aggregate_levels] += 1 # <<<<<<<<<<<<<<<<<< = 1 after show that results are the same

            previous_read_name = current_read_name

        # last group of reads
        update_results(results, sample_index, temp_results, normalise, number_samples)

    return results
