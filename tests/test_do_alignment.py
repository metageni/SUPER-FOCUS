#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from collections import defaultdict
from pathlib import Path
from unittest.mock import patch

import pytest

from superfocus_app.do_alignment import normalise_counts, update_results, parse_alignments


# ---------------------------------------------------------------------------
# normalise_counts
# ---------------------------------------------------------------------------

def test_normalise_counts_empty():
    assert normalise_counts({}) == {}


def test_normalise_counts_single():
    assert normalise_counts({"f1": 4}) == {"f1": 1.0}


def test_normalise_counts_equal():
    result = normalise_counts({"f1": 1, "f2": 1})
    assert result == {"f1": 0.5, "f2": 0.5}


def test_normalise_counts_unequal():
    result = normalise_counts({"f1": 3, "f2": 1})
    assert result["f1"] == pytest.approx(0.75)
    assert result["f2"] == pytest.approx(0.25)


# ---------------------------------------------------------------------------
# update_results
# ---------------------------------------------------------------------------

def test_update_results_empty_data():
    results = {}
    assert update_results(results, 0, {}, 0, 1) == {}


def test_update_results_new_function():
    results = {}
    assert update_results(results, 0, {"f1": 1}, 0, 1) == {"f1": [1]}


def test_update_results_accumulates():
    results = {"f1": [1]}
    assert update_results(results, 0, {"f1": 1}, 0, 1) == {"f1": [2]}


def test_update_results_multiple_samples():
    results = {"f1": [0, 0]}
    update_results(results, 1, {"f1": 5}, 0, 2)
    assert results["f1"] == [0, 5]


def test_update_results_normalised():
    results = {}
    out = update_results(results, 0, {"f1": 1, "f2": 1}, 1, 1)
    assert out["f1"] == [0.5]
    assert out["f2"] == [0.5]


def test_update_results_no_normalise_single_function():
    # normalise=1 but only 1 function → no normalisation applied
    results = {}
    out = update_results(results, 0, {"f1": 3}, 1, 1)
    assert out["f1"] == [3]


# ---------------------------------------------------------------------------
# parse_alignments  (integration — fake alignment file, no real DB)
# ---------------------------------------------------------------------------

# Fake subsystems lookup: subsystem_id → "L1\tL2\tL3"
FAKE_SUBSYSTEMS = {
    "001": "Metabolism\tCarbon\tGlycolysis",
    "002": "Metabolism\tNitrogen\tFixation",
}

# blast/diamond tabular format (outfmt 6):
# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
def _make_m8_row(read, sub_id, func, pident=80.0, length=20, evalue="1e-10"):
    # sseqid format: <prefix>__<subsystem_id>__<function>  (split("__")[1] == sub_id)
    sseqid = f"db__{sub_id}__{func}"
    return f"{read}\t{sseqid}\t{pident}\t{length}\t0\t0\t1\t{length}\t1\t{length}\t{evalue}\t100\n"


def test_parse_alignments_basic(tmp_path):
    aln = tmp_path / "sample.m8"
    aln.write_text(
        _make_m8_row("read1", "001", "GlycolysisFunc") +
        _make_m8_row("read2", "002", "FixationFunc")
    )
    results = defaultdict(list)
    binning = defaultdict(lambda: defaultdict(list))
    results, binning = parse_alignments(
        str(aln), results, normalise=0, number_samples=1, sample_index=0,
        minimum_identity=60.0, minimum_alignment=15,
        subsystems_translation=FAKE_SUBSYSTEMS,
        aligner="diamond", binning_reads=binning,
        query_name="sample.fasta", delete_alignments=False,
    )
    assert len(results) == 2


def test_parse_alignments_identity_filter(tmp_path):
    aln = tmp_path / "sample.m8"
    # pident=50 is below minimum_identity=60 → should be excluded
    aln.write_text(_make_m8_row("read1", "001", "GlycolysisFunc", pident=50.0))
    results = defaultdict(list)
    binning = defaultdict(lambda: defaultdict(list))
    results, binning = parse_alignments(
        str(aln), results, normalise=0, number_samples=1, sample_index=0,
        minimum_identity=60.0, minimum_alignment=15,
        subsystems_translation=FAKE_SUBSYSTEMS,
        aligner="diamond", binning_reads=binning,
        query_name="sample.fasta", delete_alignments=False,
    )
    assert len(results) == 0


def test_parse_alignments_alignment_length_filter(tmp_path):
    aln = tmp_path / "sample.m8"
    # length=5 is below minimum_alignment=15 → excluded
    aln.write_text(_make_m8_row("read1", "001", "GlycolysisFunc", length=5))
    results = defaultdict(list)
    binning = defaultdict(lambda: defaultdict(list))
    results, binning = parse_alignments(
        str(aln), results, normalise=0, number_samples=1, sample_index=0,
        minimum_identity=60.0, minimum_alignment=15,
        subsystems_translation=FAKE_SUBSYSTEMS,
        aligner="diamond", binning_reads=binning,
        query_name="sample.fasta", delete_alignments=False,
    )
    assert len(results) == 0


def test_parse_alignments_empty_file(tmp_path):
    aln = tmp_path / "empty.m8"
    aln.write_text("")
    # Pre-populate results from a previous sample
    results = defaultdict(list)
    results["Metabolism\tCarbon\tGlycolysis\tGlycolysisFunc"] = [1]
    binning = defaultdict(lambda: defaultdict(list))
    out_results, out_binning = parse_alignments(
        str(aln), results, normalise=0, number_samples=1, sample_index=0,
        minimum_identity=60.0, minimum_alignment=15,
        subsystems_translation=FAKE_SUBSYSTEMS,
        aligner="diamond", binning_reads=binning,
        query_name="sample.fasta", delete_alignments=False,
    )
    # Bug fix (PR #90): prior results must be preserved, not discarded
    assert "Metabolism\tCarbon\tGlycolysis\tGlycolysisFunc" in out_results


def test_parse_alignments_missing_file(tmp_path):
    # Pre-populate results from a previous sample
    results = defaultdict(list)
    results["Metabolism\tCarbon\tGlycolysis\tGlycolysisFunc"] = [1]
    binning = defaultdict(lambda: defaultdict(list))
    out_results, out_binning = parse_alignments(
        str(tmp_path / "nonexistent.m8"),
        results, normalise=0, number_samples=1, sample_index=0,
        minimum_identity=60.0, minimum_alignment=15,
        subsystems_translation=FAKE_SUBSYSTEMS,
        aligner="diamond", binning_reads=binning,
        query_name="sample.fasta", delete_alignments=False,
    )
    # Bug fix (PR #90): prior results must be preserved, not discarded
    assert "Metabolism\tCarbon\tGlycolysis\tGlycolysisFunc" in out_results


def test_parse_alignments_empty_file_preserves_binning(tmp_path):
    """Empty alignment file must not wipe binning_reads from prior samples."""
    aln = tmp_path / "empty.m8"
    aln.write_text("")
    results = defaultdict(list)
    binning = defaultdict(lambda: defaultdict(list))
    binning["prior_sample.fasta"]["read1"].append([80.0, 20, "1e-5", "some\tfunction"])
    out_results, out_binning = parse_alignments(
        str(aln), results, normalise=0, number_samples=2, sample_index=1,
        minimum_identity=60.0, minimum_alignment=15,
        subsystems_translation=FAKE_SUBSYSTEMS,
        aligner="diamond", binning_reads=binning,
        query_name="sample.fasta", delete_alignments=False,
    )
    assert "read1" in out_binning["prior_sample.fasta"]


def test_parse_alignments_delete_file(tmp_path):
    aln = tmp_path / "sample.m8"
    aln.write_text(_make_m8_row("read1", "001", "GlycolysisFunc"))
    results = defaultdict(list)
    binning = defaultdict(lambda: defaultdict(list))
    parse_alignments(
        str(aln), results, normalise=0, number_samples=1, sample_index=0,
        minimum_identity=60.0, minimum_alignment=15,
        subsystems_translation=FAKE_SUBSYSTEMS,
        aligner="diamond", binning_reads=binning,
        query_name="sample.fasta", delete_alignments=True,
    )
    assert not aln.exists()


def test_parse_alignments_binning_populated(tmp_path):
    aln = tmp_path / "sample.m8"
    aln.write_text(_make_m8_row("read1", "001", "GlycolysisFunc"))
    results = defaultdict(list)
    binning = defaultdict(lambda: defaultdict(list))
    results, binning = parse_alignments(
        str(aln), results, normalise=0, number_samples=1, sample_index=0,
        minimum_identity=60.0, minimum_alignment=15,
        subsystems_translation=FAKE_SUBSYSTEMS,
        aligner="diamond", binning_reads=binning,
        query_name="sample.fasta", delete_alignments=False,
    )
    assert "read1" in binning["sample.fasta"]


def test_parse_alignments_mmseqs_fraction_identity(tmp_path):
    """mmseqs2 reports identity as fraction (0-1); parser should multiply by 100."""
    aln = tmp_path / "sample.m8"
    # pident=0.80 → 80% after conversion, passes minimum_identity=60
    aln.write_text(_make_m8_row("read1", "001", "GlycolysisFunc", pident=0.80))
    results = defaultdict(list)
    binning = defaultdict(lambda: defaultdict(list))
    results, binning = parse_alignments(
        str(aln), results, normalise=0, number_samples=1, sample_index=0,
        minimum_identity=60.0, minimum_alignment=15,
        subsystems_translation=FAKE_SUBSYSTEMS,
        aligner="mmseqs", binning_reads=binning,
        query_name="sample.fasta", delete_alignments=False,
    )
    assert len(results) == 1


def test_parse_alignments_multiple_samples(tmp_path):
    aln1 = tmp_path / "s1.m8"
    aln2 = tmp_path / "s2.m8"
    aln1.write_text(_make_m8_row("read1", "001", "GlycolysisFunc"))
    aln2.write_text(_make_m8_row("read1", "002", "FixationFunc"))

    results = defaultdict(list)
    binning = defaultdict(lambda: defaultdict(list))

    results, binning = parse_alignments(
        str(aln1), results, normalise=0, number_samples=2, sample_index=0,
        minimum_identity=60.0, minimum_alignment=15,
        subsystems_translation=FAKE_SUBSYSTEMS,
        aligner="diamond", binning_reads=binning,
        query_name="s1.fasta", delete_alignments=False,
    )
    results, binning = parse_alignments(
        str(aln2), results, normalise=0, number_samples=2, sample_index=1,
        minimum_identity=60.0, minimum_alignment=15,
        subsystems_translation=FAKE_SUBSYSTEMS,
        aligner="diamond", binning_reads=binning,
        query_name="s2.fasta", delete_alignments=False,
    )

    # each result vector has length 2 (one per sample)
    for v in results.values():
        assert len(v) == 2
