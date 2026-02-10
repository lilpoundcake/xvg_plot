#!/usr/bin/env python3
"""
Comprehensive test runner for XVG Analysis.
Tests single file plotting, two-file comparisons, and multi-file comparisons
using test data from prot1, prot2, prot3, and other_xvg.

Results layout:
  tests/results/
  ├── prot1/          # single-file plots from prot1
  ├── prot2/          # single-file plots from prot2
  ├── prot3/          # single-file plots from prot3
  ├── other_xvg/      # single-file plots from other_xvg
  ├── *_compare.html  # two-file comparison plots
  └── *_multi_compare.html  # multi-file comparison plots
"""

import shutil
import sys
from pathlib import Path

# Add parent directory to path to import xvg_analysis
sys.path.insert(0, str(Path(__file__).parent.parent))

from xvg_analysis.xvg_analysis import (
    read_xvg, create_plot, compare_two_plots, compare_multi_plots
)


def setup_results_dir() -> Path:
    """Create clean results directory structure."""
    results_dir = Path(__file__).parent / "results"
    # Clean previous results
    if results_dir.exists():
        shutil.rmtree(results_dir)
    results_dir.mkdir()
    return results_dir


def test_single_file(file_path: Path, output_dir: Path) -> tuple[bool, str]:
    """Test single file plotting into a subdirectory."""
    try:
        df, metadata = read_xvg(str(file_path))
        fig = create_plot(df, metadata)

        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / f"{file_path.stem}_plot.html"
        fig.write_html(str(output_path))

        rel = output_path.relative_to(output_dir.parent)
        return True, f"✓ {file_path.name:25} -> {rel}"
    except Exception as e:
        return False, f"✗ {file_path.name:25} ERROR: {str(e)}"


def test_comparison(
    file1: Path, file2: Path, results_dir: Path,
    prot1: str = "prot1", prot2: str = "prot2"
) -> tuple[bool, str]:
    """Test comparison plotting into results root."""
    try:
        df1, meta1 = read_xvg(str(file1))
        df2, meta2 = read_xvg(str(file2))

        fig = compare_two_plots(df1, meta1, df2, meta2, prot1=prot1, prot2=prot2)

        output_path = results_dir / f"{file1.stem}_compare.html"
        fig.write_html(str(output_path))

        return True, f"✓ {file1.stem:15} vs {file2.stem:15} -> {output_path.name}"
    except Exception as e:
        return False, f"✗ {file1.stem:15} vs {file2.stem:15} ERROR: {str(e)}"


def test_multi_compare(
    pattern: str, output_name: str, results_dir: Path, roll_avg: int = 50
) -> tuple[bool, str]:
    """Test multi-file comparison into results root."""
    try:
        matched = sorted(Path(".").glob(pattern))
        if len(matched) < 2:
            return False, f"✗ {output_name:30} ERROR: <2 files for '{pattern}'"

        data_list = []
        for f in matched:
            df, meta = read_xvg(str(f))
            label = f"{f.parent.name}/{f.stem}"
            data_list.append((df, meta, label))

        fig = compare_multi_plots(data_list, roll_avg=roll_avg)
        output_path = results_dir / output_name
        fig.write_html(str(output_path))

        return True, f"✓ {output_name:30} ({len(matched)} files from {pattern})"
    except Exception as e:
        return False, f"✗ {output_name:30} ERROR: {str(e)}"


def main() -> int:
    """Run comprehensive test suite."""
    results_dir = setup_results_dir()
    test_dir = Path(__file__).parent

    # Input directories
    input_dirs = ["prot1", "prot2", "prot3", "other_xvg"]
    source_dirs = {name: test_dir / name for name in input_dirs}

    print("\n" + "=" * 80)
    print("XVG Analysis Test Suite")
    print("=" * 80)

    # ── Single file tests ──────────────────────────────────────────────
    print("\n📊 SINGLE FILE TESTS")
    print("-" * 80)

    single_results = []

    for dir_name, src_dir in source_dirs.items():
        if not src_dir.exists():
            continue
        xvg_files = sorted(src_dir.glob("*.xvg"))
        if not xvg_files:
            continue

        out_subdir = results_dir / dir_name
        print(f"\n  [{dir_name}]")
        for xvg_file in xvg_files:
            success, message = test_single_file(xvg_file, out_subdir)
            single_results.append((success, message))
            print(f"  {message}")

    single_success = sum(1 for s, _ in single_results if s)
    single_total = len(single_results)

    # ── Two-file comparison tests ──────────────────────────────────────
    print("\n📊 COMPARISON TESTS (prot1 vs prot2)")
    print("-" * 80)

    comparison_results = []
    prot1_dir = source_dirs["prot1"]
    prot2_dir = source_dirs["prot2"]

    if prot1_dir.exists() and prot2_dir.exists():
        for prot1_file in sorted(prot1_dir.glob("*.xvg")):
            prot2_file = prot2_dir / prot1_file.name
            if prot2_file.exists():
                success, message = test_comparison(
                    prot1_file, prot2_file, results_dir, "prot1", "prot2"
                )
                comparison_results.append((success, message))
                print(message)

    comp_success = sum(1 for s, _ in comparison_results if s)
    comp_total = len(comparison_results)

    # ── Multi-file comparison tests ────────────────────────────────────
    print("\n📊 MULTI-COMPARE TESTS (pastel color scheme)")
    print("-" * 80)

    multi_results = []

    # prot3 chain-specific comparisons
    prot3_dir = source_dirs["prot3"]
    if prot3_dir.exists():
        chain_types = [
            ("tests/prot3/rmsd_ch*.xvg", "rmsd_multi_compare.html"),
            ("tests/prot3/rmsf_ch*.xvg", "rmsf_multi_compare.html"),
            ("tests/prot3/gyrate_ch*.xvg", "gyrate_multi_compare.html"),
        ]
        for pattern, output_name in chain_types:
            success, message = test_multi_compare(pattern, output_name, results_dir)
            multi_results.append((success, message))
            print(message)

    # Cross-protein comparisons (prot1 vs prot2 via multi-compare)
    cross_types = [
        ("tests/prot*/energy.xvg", "energy_multi_compare.html"),
        ("tests/prot*/2dproj.xvg", "2dproj_multi_compare.html"),
    ]
    for pattern, output_name in cross_types:
        success, message = test_multi_compare(pattern, output_name, results_dir)
        multi_results.append((success, message))
        print(message)

    # other_xvg multi-compare (cphmd files)
    other_xvg_dir = source_dirs["other_xvg"]
    if other_xvg_dir.exists():
        success, message = test_multi_compare(
            "tests/other_xvg/cphmd*.xvg", "cphmd_multi_compare.html", results_dir
        )
        multi_results.append((success, message))
        print(message)

    multi_success = sum(1 for s, _ in multi_results if s)
    multi_total = len(multi_results)

    # ── Summary ────────────────────────────────────────────────────────
    total_success = single_success + comp_success + multi_success
    total_tests = single_total + comp_total + multi_total

    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Single file plots:     {single_success:2d}/{single_total:2d} ✓")
    print(f"Comparison plots:      {comp_success:2d}/{comp_total:2d} ✓")
    print(f"Multi-compare plots:   {multi_success:2d}/{multi_total:2d} ✓")
    print(f"Total:                 {total_success:2d}/{total_tests:2d} ✓")

    # Show directory layout
    print(f"\nResults saved to: {results_dir.resolve()}")
    for sub in sorted(results_dir.iterdir()):
        if sub.is_dir():
            count = len(list(sub.glob("*.html")))
            print(f"  {sub.name}/ ({count} files)")
        else:
            pass  # root-level files shown below
    root_files = [f for f in results_dir.iterdir() if f.is_file()]
    if root_files:
        print(f"  ./ ({len(root_files)} compare/multi-compare files)")

    total_files = sum(1 for _ in results_dir.rglob("*.html"))
    print(f"Total files generated: {total_files}")
    print("=" * 80 + "\n")

    return 0 if total_success == total_tests else 1


if __name__ == "__main__":
    sys.exit(main())
