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


def align_reads(query, output_dir, aligner, database, WORK_DIRECTORY):
    """Align FAST(A/Q) file to database.

    Args:
        query_file (PosixPath): Path to query in FAST(A/Q) file
        output_dir (PosixPath): Path to alignment output
        aligner (str): aligner name
        WORK_DIRECTORY (str): Path to directory where works happens

    Returns:
        str: Path to alignment that was written

    """
    if aligner == "diamond":
        temp_folder = Path("{}/db/tmp/".format(WORK_DIRECTORY))

        if not temp_folder.exists():
            temp_folder.mkdir(parents=True, mode=511)

        database = "{}/db/static/diamond/{}.db".format(WORK_DIRECTORY, database)
        evalue = "0.00001"
        output_name = "{}/{}_alignments".format(output_dir, query.parts[-1])

        # align
        os.system("diamond blastx -t {} -d {} -q {} -a {} -p T -e {}".format(temp_folder, database, query, output_name, evalue))
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
        subsystems_translation (dict): subsystems translation lookup table

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
            aggregate_levels = subsystems_translation[current_subsystem_id] + " " + current_function_name # <<<<<<<<<<<<<<<< CHECK JOIN

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





######


def write_results(results, header, output_name, query_path, database):
    """Write results in tabular format.

    Args:
        results (dict): dict with results to be written
        header (list): header to be wrritten
        output_name (str): Path to output
        query_path (str): Path to query
        database (str): Database used

    """
    with open(output_name, 'w') as outfile:
        writer = csv.writer(outfile, delimiter='\t', lineterminator='\n')

        # run info
        writer.writerow(["Query: {}".format(query_path)])
        writer.writerow(["Database used: {}".format(database)])
        writer.writerow([""])

        # subsystem and files header
        writer.writerow(header)
        for row in results:
            if sum(results[row]) > 0:
                writer.writerow(row.split("\t") + results[row])




#results = defaultdict(list)

#alignment = "/Users/geni.silva/Desktop/superfocus/sf2/simShort_single_sub.fasta_alignments.m8"
#number_samples = 3
#normalise = 1
#minimum_identity = 60
#minimum_alignment = 15
#sample_position = 1
#WORK_DIRECTORY = 'superfocus_app'


#subsystems_translation =  get_subsystems("/Users/geni.silva/Desktop/superfocus/superfocus_app/db/database_PKs.txt")

#results = parse_alignments(alignment, results, normalise, number_samples, sample_position, minimum_identity, minimum_alignment, subsystems_translation)


#for i in results:
#    print(i, results[i])