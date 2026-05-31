#!/usr/bin/env python3
"""Click CLI entry point for SUPER-FOCUS."""

import logging

import click

from superfocus_app import version
from superfocus_app.pipeline import run

LOGGER_FORMAT = '[%(asctime)s - %(levelname)s] %(message)s'


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.version_option(version, "-v", "--version", prog_name="SUPER-FOCUS")
@click.option("-q", "--query", required=True, multiple=True,
              help="Path to FAST(A/Q) file or directory. Repeatable.")
@click.option("-dir", "--output_directory", required=True,
              help="Path to output directory.")
@click.option("-o", "--output_prefix", default="output_", show_default=True,
              help="Output file prefix.")
@click.option("-a", "--aligner", default="diamond", show_default=True,
              help="Aligner: diamond or mmseqs2.")
@click.option("-db", "--database", default="DB_90", show_default=True,
              help="Database: DB_90, DB_95, DB_98, or DB_100.")
@click.option("-e", "--evalue", default="0.00001", show_default=True,
              help="E-value threshold.")
@click.option("-t", "--threads", default="4", show_default=True,
              help="Number of threads.")
@click.option("-mi", "--minimum_identity", default=60.0, show_default=True,
              help="Minimum percent identity.")
@click.option("-ml", "--minimum_alignment", default=15, show_default=True,
              help="Minimum alignment length (amino acids).")
@click.option("-f", "--fast", "fast_mode", default="1", show_default=True,
              help="Fast mode: 1 (fast) or 0 (sensitive).")
@click.option("-p", "--amino_acid", default="0", show_default=True,
              help="Input type: 0 nucleotides, 1 amino acids.")
@click.option("-n", "--normalise_output", default=1, show_default=True,
              help="Normalise counts: 1 yes, 0 no.")
@click.option("-m", "--focus", "run_focus", default="0", show_default=True,
              help="Run FOCUS reduction: 1 yes, 0 no.")
@click.option("-b", "--alternate_directory", default="",
              help="Alternate directory for databases.")
@click.option("-d", "--delete_alignments", is_flag=True, default=False,
              help="Delete alignment files after parsing.")
@click.option("-s", "--subsample", default=None, type=int,
              help="Subsample reads to this count.")
@click.option("-w", "--latency_wait", default=0, show_default=True,
              help="Seconds to wait after alignment (NFS latency).")
@click.option("-tmp", "--temp_directory", default=None,
              help="Alternate temporary directory.")
@click.option("-l", "--log", default=None,
              help="Path to log file (default: STDOUT).")
def main(query, output_directory, output_prefix, aligner, database, evalue, threads,
         minimum_identity, minimum_alignment, fast_mode, amino_acid, normalise_output,
         run_focus, alternate_directory, delete_alignments, subsample, latency_wait,
         temp_directory, log):
    """SUPER-FOCUS: agile functional analysis of shotgun metagenomic data."""
    logging.basicConfig(format=LOGGER_FORMAT, level=logging.INFO, filename=log)

    run(
        query_paths=list(query),
        output_directory=output_directory,
        prefix=output_prefix,
        aligner=aligner,
        database=database,
        evalue=evalue,
        threads=threads,
        fast_mode=fast_mode,
        amino_acid=amino_acid,
        minimum_identity=minimum_identity,
        minimum_alignment=minimum_alignment,
        normalise_output=normalise_output,
        alternate_directory=alternate_directory,
        delete_alignments=delete_alignments,
        subsample=subsample,
        latency_wait=latency_wait,
        run_focus=run_focus,
        temp_directory=temp_directory,
    )
