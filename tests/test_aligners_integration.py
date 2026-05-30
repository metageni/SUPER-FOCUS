#!/usr/bin/env python3
"""End-to-end CLI integration tests for SUPER-FOCUS.

Fixtures (tests/fixtures/):
  db.fasta        — 20 protein sequences across 10 subsystems (2 per subsystem)
  query.fasta     — 10 reads, one per subsystem, identical to DB entries
  work_dir/       — fake WORK_DIRECTORY with indexed DBs and database_PKs.txt

Subsystems covered (PK → L1 / L2 / L3):
  55  Amino Acids / Alanine, serine, and glycine / Sarcosine oxidases
  62  Amino Acids / Alanine, serine, and glycine / Serine Biosynthesis
  190 Amino Acids / Lysine... / Threonine and Homoserine Biosynthesis
  191 Amino Acids / Lysine... / Threonine degradation
  233 Amino Acids / Aromatic... / Tryptophan synthesis
  276 Amino Acids / Arginine... / Urea decomposition
  277 Amino Acids / Arginine... / Urease subunits
  283 Amino Acids / Branched-chain / Valine degradation
  424 Amino Acids / Branched-chain / BCAA Biosynthesis
  572 Amino Acids / Aromatic... / Chorismate Synthesis
"""

import subprocess
import sys
import csv
import pytest
from pathlib import Path

FIXTURES = Path(__file__).parent / "fixtures"
QUERY = str(FIXTURES / "query.fasta")
WORK_DIR = str(FIXTURES / "work_dir")
OLD_PKG = Path("/Users/gsilva/Desktop/sf/old/SUPER-FOCUS")

EXPECTED_OUTPUTS = [
    "output_all_levels_and_function.xls",
    "output_subsystem_level_1.xls",
    "output_subsystem_level_2.xls",
    "output_subsystem_level_3.xls",
    "output_binning.xls",
]

# Exact 20 functions expected in all_levels (2 per subsystem × 10 subsystems)
EXPECTED_FUNCTIONS = {
    "Sarcosine_oxidase_alpha", "Sarcosine_oxidase_beta",
    "Serine_hydroxymethyltransferase", "Serine_acetyltransferase",
    "Threonine_synthase", "Homoserine_kinase",
    "Threonine_dehydratase", "Threonine_aldolase",
    "Tryptophan_synthase_alpha", "Tryptophan_synthase_beta",
    "Urease_accessory_protein_UreD", "Urease_accessory_protein_UreF",
    "Urease_alpha_subunit", "Urease_beta_subunit",
    "Valine_dehydrogenase", "Branched_chain_aminotransferase",
    "BCAA_biosynthesis_enzyme", "Acetolactate_synthase",
    "Chorismate_synthase", "DAHP_synthase",
}

# Expected L2 subsystems
EXPECTED_L2 = {
    "Alanine, serine, and glycine",
    "Arginine; urea cycle, polyamines",
    "Aromatic amino acids and derivatives",
    "Branched-chain amino acids",
    "Lysine, threonine, methionine, and cysteine",
}

# All 10 reads in query.fasta
EXPECTED_READS = {
    "read_sarcosine_1", "read_sarcosine_2",
    "read_serine_1", "read_threonine_1",
    "read_tryptophan_1", "read_urease_1",
    "read_valine_1", "read_bcaa_1",
    "read_chorismate_1", "read_urea_decomp_1",
}


def _run(module, aligner, out_dir):
    """Invoke a superfocus CLI module as a subprocess.

    Args:
        module (str): Python module path (e.g. 'superfocus_app.superfocus').
        aligner (str): Aligner name.
        out_dir (Path): Output directory.

    Returns:
        subprocess.CompletedProcess: Result of the invocation.
    """
    return subprocess.run(
        [
            sys.executable, "-m", module,
            "-q", QUERY,
            "-dir", str(out_dir),
            "-a", aligner,
            "-db", "DB_90",
            "-b", WORK_DIR,
            "-p", "1",
            "-e", "0.001",
            "-t", "1",
            "-n", "0",
        ],
        capture_output=True, text=True,
        cwd=str(OLD_PKG) if "old" in module else "/Users/gsilva/Desktop/sf/SUPER-FOCUS",
    )


def _run_new(aligner, out_dir):
    """Run the current (new) superfocus version.

    Args:
        aligner (str): Aligner name.
        out_dir (Path): Output directory.

    Returns:
        subprocess.CompletedProcess: Result of the invocation.
    """
    return _run("superfocus_app.superfocus", aligner, out_dir)


def _data_rows(path, skip=5):
    """Return non-empty data rows from a TSV output file, skipping the header block.

    Args:
        path (Path): Path to the .xls output file.
        skip (int): Number of header lines to skip.

    Returns:
        list: List of tab-split row lists.
    """
    lines = path.read_text().splitlines()[skip:]
    return [l.split("\t") for l in lines if l.strip()]


def _assert_outputs(out_dir):
    """Assert all expected output files exist, are non-empty, and contain correct hits.

    Args:
        out_dir (Path): Directory where superfocus wrote its output.

    Returns:
        None
    """
    for fname in EXPECTED_OUTPUTS:
        p = out_dir / fname
        assert p.exists(), f"Missing: {fname}"
        assert p.stat().st_size > 0, f"Empty: {fname}"

    all_levels = out_dir / "output_all_levels_and_function.xls"
    rows = _data_rows(all_levels)
    assert len(rows) == 20, f"Expected 20 function rows, got {len(rows)}"

    content = all_levels.read_text()
    for fn in EXPECTED_FUNCTIONS:
        assert fn in content, f"Function missing: {fn}"

    # L1: only one category, 100% relative abundance
    lvl1_rows = _data_rows(out_dir / "output_subsystem_level_1.xls")
    assert len(lvl1_rows) == 1
    assert lvl1_rows[0][0] == "Amino Acids and Derivatives"
    assert float(lvl1_rows[0][2]) == pytest.approx(100.0)

    # L2: exactly 5 categories
    lvl2_rows = _data_rows(out_dir / "output_subsystem_level_2.xls")
    assert len(lvl2_rows) == 5
    assert {r[0] for r in lvl2_rows} == EXPECTED_L2

    # L2 relative abundances sum to 100%
    l2_pct = sum(float(r[2]) for r in lvl2_rows)
    assert l2_pct == pytest.approx(100.0, abs=0.01)

    # L3: exactly 10 subsystems
    lvl3_rows = _data_rows(out_dir / "output_subsystem_level_3.xls")
    assert len(lvl3_rows) == 10

    # binning: all 10 reads present, identity = 100%
    binning = (out_dir / "output_binning.xls").read_text()
    for read in EXPECTED_READS:
        assert read in binning, f"Read missing from binning: {read}"
    assert "100.0" in binning, "Expected 100% identity hits in binning"


# ---------------------------------------------------------------------------
# helpers shared across per-aligner parametrized tests
# ---------------------------------------------------------------------------

ALIGNERS = pytest.mark.parametrize("aligner", ["diamond", "mmseqs2", "blast"])


# ---------------------------------------------------------------------------
# exit code
# ---------------------------------------------------------------------------

@ALIGNERS
@pytest.mark.integration
def test_cli_exits_zero(aligner, tmp_path):
    """CLI exits with code 0 for each aligner.

    Args:
        aligner (str): Aligner under test.
        tmp_path (Path): pytest temporary directory.

    Returns:
        None
    """
    result = _run_new(aligner, tmp_path)
    assert result.returncode == 0, f"{aligner} failed:\n{result.stderr}"


# ---------------------------------------------------------------------------
# output correctness
# ---------------------------------------------------------------------------

@ALIGNERS
@pytest.mark.integration
def test_cli_all_outputs_correct(aligner, tmp_path):
    """All 5 output files are produced with correct content for each aligner.

    Args:
        aligner (str): Aligner under test.
        tmp_path (Path): pytest temporary directory.

    Returns:
        None
    """
    _run_new(aligner, tmp_path)
    _assert_outputs(tmp_path)


@ALIGNERS
@pytest.mark.integration
def test_cli_all_20_functions_present(aligner, tmp_path):
    """All 20 expected functions appear in all_levels output.

    Args:
        aligner (str): Aligner under test.
        tmp_path (Path): pytest temporary directory.

    Returns:
        None
    """
    _run_new(aligner, tmp_path)
    content = (tmp_path / "output_all_levels_and_function.xls").read_text()
    missing = [fn for fn in EXPECTED_FUNCTIONS if fn not in content]
    assert not missing, f"Missing functions: {missing}"


@ALIGNERS
@pytest.mark.integration
def test_cli_relative_abundance_sums_to_100(aligner, tmp_path):
    """Relative abundance at every subsystem level sums to 100%.

    Args:
        aligner (str): Aligner under test.
        tmp_path (Path): pytest temporary directory.

    Returns:
        None
    """
    _run_new(aligner, tmp_path)
    for level in (1, 2, 3):
        rows = _data_rows(tmp_path / f"output_subsystem_level_{level}.xls")
        total = sum(float(r[2]) for r in rows)
        assert total == pytest.approx(100.0, abs=0.01), \
            f"L{level} relative abundance sums to {total}, not 100"


@ALIGNERS
@pytest.mark.integration
def test_cli_binning_identity_and_reads(aligner, tmp_path):
    """Binning output contains all reads with 100% identity hits.

    Args:
        aligner (str): Aligner under test.
        tmp_path (Path): pytest temporary directory.

    Returns:
        None
    """
    _run_new(aligner, tmp_path)
    rows = _data_rows(tmp_path / "output_binning.xls", skip=5)
    read_names = {r[1] for r in rows if len(r) > 1}
    for read in EXPECTED_READS:
        assert read in read_names, f"Read missing from binning: {read}"
    identities = [float(r[6]) for r in rows if len(r) > 6]
    assert all(i >= 60.0 for i in identities), \
        f"Binning hit below minimum_identity=60: {[i for i in identities if i < 60]}"


@ALIGNERS
@pytest.mark.integration
def test_cli_custom_prefix(aligner, tmp_path):
    """Custom output prefix is applied to all output files.

    Args:
        aligner (str): Aligner under test.
        tmp_path (Path): pytest temporary directory.

    Returns:
        None
    """
    subprocess.run(
        [sys.executable, "-m", "superfocus_app.superfocus",
         "-q", QUERY, "-dir", str(tmp_path),
         "-a", aligner, "-db", "DB_90", "-b", WORK_DIR,
         "-p", "1", "-e", "0.001", "-t", "1", "-n", "0",
         "-o", "myrun_"],
        capture_output=True,
    )
    for fname in EXPECTED_OUTPUTS:
        assert (tmp_path / fname.replace("output_", "myrun_")).exists(), \
            f"Custom-prefix file missing: {fname}"


@ALIGNERS
@pytest.mark.integration
def test_cli_delete_alignments_flag(aligner, tmp_path):
    """Alignment files are deleted when -d flag is set.

    Args:
        aligner (str): Aligner under test.
        tmp_path (Path): pytest temporary directory.

    Returns:
        None
    """
    subprocess.run(
        [sys.executable, "-m", "superfocus_app.superfocus",
         "-q", QUERY, "-dir", str(tmp_path),
         "-a", aligner, "-db", "DB_90", "-b", WORK_DIR,
         "-p", "1", "-e", "0.001", "-t", "1", "-n", "0", "-d"],
        capture_output=True,
    )
    m8_files = list(tmp_path.glob("*.m8"))
    assert not m8_files, f"Alignment files not deleted: {m8_files}"


# ---------------------------------------------------------------------------
# cross-aligner consistency
# ---------------------------------------------------------------------------

@pytest.mark.integration
def test_all_aligners_agree_on_subsystem_level1(tmp_path):
    """diamond, mmseqs2, and blast all report identical L1 subsystems.

    Args:
        tmp_path (Path): pytest temporary directory.

    Returns:
        None
    """
    results = {}
    for aligner in ("diamond", "mmseqs2", "blast"):
        out = tmp_path / aligner
        out.mkdir()
        _run_new(aligner, out)
        rows = _data_rows(out / "output_subsystem_level_1.xls")
        results[aligner] = {r[0] for r in rows}

    assert results["diamond"] == results["mmseqs2"] == results["blast"], \
        f"L1 subsystems differ: {results}"


@pytest.mark.integration
def test_all_aligners_agree_on_subsystem_level2(tmp_path):
    """diamond, mmseqs2, and blast all report identical L2 subsystems.

    Args:
        tmp_path (Path): pytest temporary directory.

    Returns:
        None
    """
    results = {}
    for aligner in ("diamond", "mmseqs2", "blast"):
        out = tmp_path / aligner
        out.mkdir()
        _run_new(aligner, out)
        rows = _data_rows(out / "output_subsystem_level_2.xls")
        results[aligner] = {r[0] for r in rows}

    assert results["diamond"] == results["mmseqs2"] == results["blast"], \
        f"L2 subsystems differ: {results}"


@pytest.mark.integration
def test_all_aligners_agree_on_function_set(tmp_path):
    """diamond, mmseqs2, and blast all report the same set of functions.

    Args:
        tmp_path (Path): pytest temporary directory.

    Returns:
        None
    """
    results = {}
    for aligner in ("diamond", "mmseqs2", "blast"):
        out = tmp_path / aligner
        out.mkdir()
        _run_new(aligner, out)
        rows = _data_rows(out / "output_all_levels_and_function.xls")
        results[aligner] = {r[3] for r in rows if len(r) > 3}

    assert results["diamond"] == results["mmseqs2"] == results["blast"], \
        f"Function sets differ: {results}"


# ---------------------------------------------------------------------------
# new vs old version parity
# ---------------------------------------------------------------------------

@ALIGNERS
@pytest.mark.integration
def test_new_vs_old_all_levels_identical(aligner, tmp_path):
    """New and old versions produce byte-identical all_levels_and_function output.

    Args:
        aligner (str): Aligner under test.
        tmp_path (Path): pytest temporary directory.

    Returns:
        None
    """
    new_out = tmp_path / "new"
    old_out = tmp_path / "old"
    new_out.mkdir()
    old_out.mkdir()

    _run_new(aligner, new_out)
    _run("superfocus_app.superfocus", aligner, old_out)

    new_content = (new_out / "output_all_levels_and_function.xls").read_text()
    old_content = (old_out / "output_all_levels_and_function.xls").read_text()
    assert new_content == old_content, \
        f"all_levels_and_function differs between new and old for {aligner}"


@ALIGNERS
@pytest.mark.integration
def test_new_vs_old_all_files_identical(aligner, tmp_path):
    """New and old versions produce byte-identical output for all 5 output files.

    Args:
        aligner (str): Aligner under test.
        tmp_path (Path): pytest temporary directory.

    Returns:
        None
    """
    new_out = tmp_path / "new"
    old_out = tmp_path / "old"
    new_out.mkdir()
    old_out.mkdir()

    _run_new(aligner, new_out)
    _run("superfocus_app.superfocus", aligner, old_out)

    for fname in EXPECTED_OUTPUTS:
        new_c = (new_out / fname).read_text()
        old_c = (old_out / fname).read_text()
        assert new_c == old_c, f"{fname} differs between new and old for {aligner}"
