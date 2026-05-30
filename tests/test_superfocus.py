#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
import os
from pathlib import Path
from collections import defaultdict

import numpy as np
import pytest

from superfocus_app.superfocus import (
    is_wanted_file,
    is_valid_number,
    get_denominators,
    add_relative_abundance,
    aggregate_level,
    get_subsystems,
    write_results,
    write_binning,
    subsample_reads,
)


# ---------------------------------------------------------------------------
# is_valid_number
# ---------------------------------------------------------------------------

def test_is_valid_number_valid():
    for v in ("1", "1.0", "0.00001", "0", "100"):
        assert is_valid_number(v) is True


def test_is_valid_number_invalid():
    for v in ("-0.00001", "-1", "-1.0", "aa", "1a", "-1a", ""):
        assert is_valid_number(v) is False


# ---------------------------------------------------------------------------
# is_wanted_file
# ---------------------------------------------------------------------------

def test_is_wanted_file_keeps_valid(tmp_path):
    files = [tmp_path / n for n in ("a.fasta", "b.fastq", "c.fna")]
    for f in files:
        f.touch()
    result = is_wanted_file(files)
    assert [p.name for p in result] == ["a.fasta", "b.fastq", "c.fna"]


def test_is_wanted_file_filters_invalid(tmp_path):
    files = [tmp_path / n for n in ("a.fasta", "b.png", "c.txt")]
    for f in files:
        f.touch()
    result = is_wanted_file(files)
    assert len(result) == 1
    assert result[0].name == "a.fasta"


def test_is_wanted_file_case_insensitive(tmp_path):
    files = [tmp_path / n for n in ("A.FASTA", "B.FASTQ", "C.FNA")]
    for f in files:
        f.touch()
    assert len(is_wanted_file(files)) == 3


def test_is_wanted_file_empty():
    assert is_wanted_file([]) == []


# ---------------------------------------------------------------------------
# get_denominators
# ---------------------------------------------------------------------------

def test_get_denominators_basic():
    assert list(get_denominators({"a": [1, 2], "b": [3, 4]})) == [4, 6]


def test_get_denominators_zeros():
    assert list(get_denominators({"a": [0, 0], "b": [0, 0]})) == [0, 0]


def test_get_denominators_three_samples():
    r = {"a": [1, 2, 3], "b": [4, 8, 10]}
    assert list(get_denominators(r)) == [5, 10, 13]


# ---------------------------------------------------------------------------
# add_relative_abundance
# ---------------------------------------------------------------------------

def test_add_relative_abundance_basic():
    results = {"f1": [1, 3]}
    normalizer = np.array([2, 6])
    out = add_relative_abundance(results, normalizer)
    assert out["f1"] == [1, 3, 50.0, 50.0]


def test_add_relative_abundance_zero_normalizer():
    results = {"f1": [0, 0]}
    normalizer = np.array([0, 0])
    # Should not raise; values where normalizer==0 are undefined but the call completes
    out = add_relative_abundance(results, normalizer)
    assert len(out["f1"]) == 4  # raw counts + relative abundance appended


# ---------------------------------------------------------------------------
# aggregate_level
# ---------------------------------------------------------------------------

def test_aggregate_level_position_0():
    results = {
        "a1\ta2\ta3\tfn": [1, 2, 3],
        "b1\tb2\tb3\tfn": [4, 8, 10],
        "c1\tc2\tc3\tfn": [1, 2, 3],
    }
    normalizer = np.array([6, 12, 16])
    out = aggregate_level(results, 0, normalizer)
    assert list(out["a1"][:3]) == [1, 2, 3]
    assert list(out["b1"][:3]) == [4, 8, 10]


def test_aggregate_level_merges_same_level():
    # two entries share the same level-0 value → should be summed
    results = {
        "shared\tl2a\tl3a\tfn": [1, 0],
        "shared\tl2b\tl3b\tfn": [2, 0],
    }
    normalizer = np.array([3, 0])
    out = aggregate_level(results, 0, normalizer)
    assert list(out["shared"][:2]) == [3, 0]


# ---------------------------------------------------------------------------
# get_subsystems
# ---------------------------------------------------------------------------

def test_get_subsystems(tmp_path):
    db = tmp_path / "database_PKs.txt"
    db.write_text("PK\tL1\tL2\tL3\n001\tMeta\tSub\tFunc\n002\tOther\tX\tY\n")
    result = get_subsystems(db)
    assert result["001"] == "Meta\tSub\tFunc"
    assert result["002"] == "Other\tX\tY"


# ---------------------------------------------------------------------------
# write_results
# ---------------------------------------------------------------------------

def test_write_results(tmp_path):
    results = {"L1\tL2\tL3\tFn": [5, 3]}
    header = ["Subsystem 1", "s1", "s2"]
    out = str(tmp_path / "out.xls")
    write_results(results, header, out, "query.fasta", "90", "diamond")
    lines = Path(out).read_text().splitlines()
    assert any("diamond" in l for l in lines)
    assert any("L1" in l for l in lines)


def test_write_results_skips_zero_rows(tmp_path):
    results = {"L1\tL2\tL3\tFn": [0, 0], "L1\tL2\tL3\tFn2": [1, 0]}
    out = str(tmp_path / "out.xls")
    write_results(results, ["H", "s1", "s2"], out, "q", "90", "diamond")
    content = Path(out).read_text()
    assert "Fn2" in content
    assert "Fn\t" not in content  # zero row excluded


# ---------------------------------------------------------------------------
# write_binning
# ---------------------------------------------------------------------------

def test_write_binning(tmp_path):
    binning = defaultdict(lambda: defaultdict(list))
    binning["sample.fasta"]["read1"].append([75.0, 20.0, "1e-5", "L1\tL2\tL3\tFn"])
    out = str(tmp_path / "binning.xls")
    write_binning(binning, out, "sample.fasta", "90", "diamond")
    content = Path(out).read_text()
    assert "read1" in content
    assert "sample.fasta" in content


# ---------------------------------------------------------------------------
# subsample_reads  (fasta)
# ---------------------------------------------------------------------------

def test_subsample_reads_fasta(tmp_path):
    fasta = tmp_path / "seqs.fasta"
    fasta.write_text(">r1\nACGT\n>r2\nTTTT\n>r3\nGGGG\n")
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    result = subsample_reads(fasta, str(out_dir), 2)
    content = result.read_text()
    assert content.count(">") == 2


def test_subsample_reads_fastq(tmp_path):
    fastq = tmp_path / "seqs.fastq"
    fastq.write_text(
        "@r1\nACGT\n+\nIIII\n"
        "@r2\nTTTT\n+\nIIII\n"
        "@r3\nGGGG\n+\nIIII\n"
    )
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    result = subsample_reads(fastq, str(out_dir), 2)
    lines = result.read_text().splitlines()
    assert len(lines) == 8  # 2 reads × 4 lines
