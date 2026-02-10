#!/usr/bin/env python3
"""
XVG Plot - Automatic plot type detection for GROMACS XVG files.

A console utility for visualizing GROMACS simulation data with interactive Plotly plots.

Usage examples:
    # Plot any file (auto-detects type)
    xvg_plot plot energy.xvg

    # Compare two files (auto-detects type)
    xvg_plot compare prot1/energy.xvg prot2/energy.xvg

    # Save plot to output directory
    xvg_plot plot energy.xvg -o output.html

    # With rolling average
    xvg_plot plot energy.xvg --roll-avg 20

    # Detect only
    xvg_plot detect energy.xvg

    # Batch process
    xvg_plot batch prot1/ -o results/
"""

import argparse
import glob
import os
import sys
from pathlib import Path
from typing import Optional, List

from . import xvg_analysis as xp


def add_output_arg(parser: argparse.ArgumentParser) -> None:
    """Add output file argument to parser."""
    parser.add_argument(
        '-o', '--output',
        type=str,
        nargs='?',
        const=None,
        help='Output file path (HTML for interactive, PNG/SVG for static)'
    )


def add_roll_avg_arg(parser: argparse.ArgumentParser, default: int = 10) -> None:
    """Add rolling average argument to parser."""
    parser.add_argument(
        '--roll-avg',
        type=int,
        default=default,
        help=f'Rolling average window size (default: {default})'
    )


def add_compare_args(parser: argparse.ArgumentParser) -> None:
    """Add common arguments for compare commands."""
    parser.add_argument(
        'file1',
        type=str,
        help='Path to first XVG file'
    )
    parser.add_argument(
        'file2',
        type=str,
        help='Path to second XVG file'
    )
    parser.add_argument(
        '--prot1',
        type=str,
        default='prot_1',
        help='Name for first protein (default: prot_1)'
    )
    parser.add_argument(
        '--prot2',
        type=str,
        default='prot_2',
        help='Name for second protein (default: prot_2)'
    )
    add_roll_avg_arg(parser, default=50)
    add_output_arg(parser)


def open_in_browser(fig) -> None:
    """Open plotly figure in browser."""
    try:
        fig.show()
    except Exception as e:
        print(f"Warning: Could not open in browser: {e}", file=sys.stderr)


def save_figure(fig, output_path: str) -> None:
    """Save figure to HTML or image file."""
    path = Path(output_path)
    suffix = path.suffix.lower()

    if suffix == '.html':
        fig.write_html(output_path, include_plotlyjs='cdn')
        print(f"Saved HTML to: {output_path}")
    elif suffix in ('.png', '.jpg', '.jpeg'):
        fig.write_image(output_path, scale=2)
        print(f"Saved image to: {output_path}")
    elif suffix in ('.svg', '.pdf'):
        fig.write_image(output_path)
        print(f"Saved vector graphic to: {output_path}")
    else:
        # Default to HTML
        fig.write_html(output_path, include_plotlyjs='cdn')
        print(f"Saved HTML to: {output_path}")


def cmd_plot(args: argparse.Namespace) -> int:
    """Handle plot command with automatic type detection."""
    file_path = args.file

    if not os.path.exists(file_path):
        print(f"Error: File not found - {file_path}", file=sys.stderr)
        return 1

    df, metadata = xp.read_xvg(file_path)

    if df.empty:
        print(f"Error: Could not parse data from {file_path}", file=sys.stderr)
        return 1

    print(f"File: {file_path}")
    print(f"  Plot Type: {metadata.plot_type}")
    print(f"  Title: {metadata.title}")
    print(f"  X-axis: {metadata.x_label}")
    print(f"  Y-axis: {metadata.y_label}")
    print(f"  Columns: {list(df.columns)}")
    print(f"  Data points: {len(df)}")

    # Use calculated roll_avg from data if not explicitly set
    if args.roll_avg is None:
        roll_avg = xp.calculate_roll_avg(df)
        print(f"  Calculated roll-avg: {roll_avg}")
    else:
        roll_avg = args.roll_avg

    fig = xp.create_plot(df, metadata, roll_avg=roll_avg)

    open_in_browser(fig)

    # Determine output file path
    if args.output:
        output_path = args.output
    else:
        # Use basename of input file with .html extension in current directory
        stem = Path(file_path).stem
        output_path = f"{stem}_plot.html"
        print(f"  Saving to: {output_path}")

    save_figure(fig, output_path)

    return 0


def cmd_compare(args: argparse.Namespace) -> int:
    """Handle compare command with automatic type detection."""
    file1 = args.file1
    file2 = args.file2

    if not os.path.exists(file1):
        print(f"Error: File not found - {file1}", file=sys.stderr)
        return 1
    if not os.path.exists(file2):
        print(f"Error: File not found - {file2}", file=sys.stderr)
        return 1

    df1, meta1 = xp.read_xvg(file1)
    df2, meta2 = xp.read_xvg(file2)

    if df1.empty:
        print(f"Error: Could not parse data from {file1}", file=sys.stderr)
        return 1
    if df2.empty:
        print(f"Error: Could not parse data from {file2}", file=sys.stderr)
        return 1

    print(f"File 1: {file1}")
    print(f"  Plot Type: {meta1.plot_type}")
    print(f"  Title: {meta1.title}")
    print(f"  Data points: {len(df1)}")
    print(f"File 2: {file2}")
    print(f"  Plot Type: {meta2.plot_type}")
    print(f"  Title: {meta2.title}")
    print(f"  Data points: {len(df2)}")

    # Validate both files have same plot type for comparison
    if meta1.plot_type != meta2.plot_type:
        print(f"Warning: Files have different types ({meta1.plot_type} vs {meta2.plot_type})", file=sys.stderr)

    fig = xp.compare_two_plots(df1, meta1, df2, meta2,
                                prot1=args.prot1, prot2=args.prot2,
                                roll_avg=args.roll_avg)

    open_in_browser(fig)

    if args.output:
        save_figure(fig, args.output)

    return 0


def cmd_pca_compare(args: argparse.Namespace) -> int:
    """Handle PCA comparison."""
    file1 = args.file1
    file2 = args.file2

    if not os.path.exists(file1):
        print(f"Error: File not found - {file1}", file=sys.stderr)
        return 1
    if not os.path.exists(file2):
        print(f"Error: File not found - {file2}", file=sys.stderr)
        return 1

    df1, meta1 = xp.read_xvg(file1)
    df2, meta2 = xp.read_xvg(file2)

    if df1.empty:
        print(f"Error: Could not parse data from {file1}", file=sys.stderr)
        return 1
    if df2.empty:
        print(f"Error: Could not parse data from {file2}", file=sys.stderr)
        return 1

    print(f"File 1: {file1}")
    print(f"  Plot Type: {meta1.plot_type}")
    print(f"  Title: {meta1.title}")
    print(f"  Data points: {len(df1)}")
    print(f"File 2: {file2}")
    print(f"  Plot Type: {meta2.plot_type}")
    print(f"  Title: {meta2.title}")
    print(f"  Data points: {len(df2)}")

    fig = xp.compare_pca_plots(df1, meta1, df2, meta2,
                               prot1=args.prot1, prot2=args.prot2)

    open_in_browser(fig)

    if args.output:
        save_figure(fig, args.output)

    return 0


def cmd_rmsf_compare(args: argparse.Namespace) -> int:
    """Handle RMSF comparison."""
    file1 = args.file1
    file2 = args.file2

    if not os.path.exists(file1):
        print(f"Error: File not found - {file1}", file=sys.stderr)
        return 1
    if not os.path.exists(file2):
        print(f"Error: File not found - {file2}", file=sys.stderr)
        return 1

    df1, meta1 = xp.read_xvg(file1)
    df2, meta2 = xp.read_xvg(file2)

    if df1.empty:
        print(f"Error: Could not parse data from {file1}", file=sys.stderr)
        return 1
    if df2.empty:
        print(f"Error: Could not parse data from {file2}", file=sys.stderr)
        return 1

    print(f"File 1: {file1}")
    print(f"  Plot Type: {meta1.plot_type}")
    print(f"  Title: {meta1.title}")
    print(f"  Data points: {len(df1)}")
    print(f"File 2: {file2}")
    print(f"  Plot Type: {meta2.plot_type}")
    print(f"  Title: {meta2.title}")
    print(f"  Data points: {len(df2)}")

    fig = xp.compare_rmsf_plots(df1, meta1, df2, meta2,
                                prot1=args.prot1, prot2=args.prot2)

    open_in_browser(fig)

    if args.output:
        save_figure(fig, args.output)

    return 0


def cmd_sasa_compare(args: argparse.Namespace) -> int:
    """Handle SASA comparison."""
    file1 = args.file1
    file2 = args.file2

    if not os.path.exists(file1):
        print(f"Error: File not found - {file1}", file=sys.stderr)
        return 1
    if not os.path.exists(file2):
        print(f"Error: File not found - {file2}", file=sys.stderr)
        return 1

    df1, meta1 = xp.read_xvg(file1)
    df2, meta2 = xp.read_xvg(file2)

    if df1.empty:
        print(f"Error: Could not parse data from {file1}", file=sys.stderr)
        return 1
    if df2.empty:
        print(f"Error: Could not parse data from {file2}", file=sys.stderr)
        return 1

    print(f"File 1: {file1}")
    print(f"  Plot Type: {meta1.plot_type}")
    print(f"  Title: {meta1.title}")
    print(f"  Data points: {len(df1)}")
    print(f"File 2: {file2}")
    print(f"  Plot Type: {meta2.plot_type}")
    print(f"  Title: {meta2.title}")
    print(f"  Data points: {len(df2)}")

    fig = xp.compare_sasa_plots(df1, meta1, df2, meta2,
                                prot1=args.prot1, prot2=args.prot2,
                                roll_avg=args.roll_avg)

    open_in_browser(fig)

    if args.output:
        save_figure(fig, args.output)

    return 0


def cmd_multi_compare(args: argparse.Namespace) -> int:
    """Handle multi-compare command - compare N XVG files matching a glob pattern."""
    pattern = args.pattern

    # Expand glob pattern to find matching files
    matched_files = sorted(glob.glob(pattern))

    # Filter to .xvg files only
    matched_files = [f for f in matched_files if f.endswith('.xvg')]

    if len(matched_files) < 2:
        print(f"Error: Need at least 2 .xvg files matching '{pattern}', "
              f"found {len(matched_files)}", file=sys.stderr)
        if matched_files:
            print(f"  Found: {matched_files}", file=sys.stderr)
        return 1

    print(f"Pattern: {pattern}")
    print(f"Found {len(matched_files)} XVG file(s):")

    # Parse all files
    data_list: List[tuple] = []
    detected_types: List[str] = []

    for file_path in matched_files:
        if not os.path.exists(file_path):
            print(f"  Warning: Skipping {file_path} (not found)", file=sys.stderr)
            continue

        df, metadata = xp.read_xvg(file_path)

        if df.empty:
            print(f"  Warning: Skipping {file_path} (could not parse)", file=sys.stderr)
            continue

        # Use filename stem as label, or parent_dir/stem for disambiguation
        label = Path(file_path).stem
        # If multiple files have same stem, use parent/stem
        if sum(1 for f in matched_files if Path(f).stem == label) > 1:
            label = f"{Path(file_path).parent.name}/{label}"

        print(f"  {file_path}: type={metadata.plot_type}, "
              f"points={len(df)}, label={label}")

        data_list.append((df, metadata, label))
        detected_types.append(metadata.plot_type)

    if len(data_list) < 2:
        print("Error: Need at least 2 valid files for comparison", file=sys.stderr)
        return 1

    # Validate all files have the same type
    unique_types = set(detected_types)
    if len(unique_types) > 1:
        print(f"\nError: Files have different plot types: {unique_types}",
              file=sys.stderr)
        print("Multi-compare only works with files of the same type.",
              file=sys.stderr)
        return 1

    print(f"\nPlot type: {detected_types[0]}")
    print(f"Creating multi-comparison plot for {len(data_list)} files...")

    fig = xp.compare_multi_plots(data_list, roll_avg=args.roll_avg)

    open_in_browser(fig)

    # Determine output file path
    if args.output:
        output_path = args.output
    else:
        # Auto-generate name from pattern
        base = Path(pattern).stem.replace('*', 'multi').replace('?', '_')
        output_path = f"{base}_multi_compare.html"
        print(f"  Saving to: {output_path}")

    save_figure(fig, output_path)

    return 0


def cmd_detect(args: argparse.Namespace) -> int:
    """Handle detect command - only detect and print plot type."""
    file_path = args.file

    if not os.path.exists(file_path):
        print(f"Error: File not found - {file_path}", file=sys.stderr)
        return 1

    df, metadata = xp.read_xvg(file_path)

    print(f"File: {file_path}")
    print(f"  Plot Type: {metadata.plot_type}")
    print(f"  Title: {metadata.title}")
    print(f"  X-axis Label: {metadata.x_label}")
    print(f"  Y-axis Label: {metadata.y_label}")
    print(f"  Number of Columns: {metadata.n_columns}")
    print(f"  Column Names: {list(df.columns)}")

    return 0


def cmd_batch(args: argparse.Namespace) -> int:
    """Handle batch command - process all XVG files in directory."""
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)

    if not input_dir.exists():
        print(f"Error: Input directory not found - {input_dir}", file=sys.stderr)
        return 1

    # Create output directory if needed
    output_dir.mkdir(parents=True, exist_ok=True)

    # Find all XVG files
    xvg_files = sorted(input_dir.glob('*.xvg'))

    if not xvg_files:
        print(f"Error: No XVG files found in {input_dir}", file=sys.stderr)
        return 1

    print(f"Found {len(xvg_files)} XVG file(s) to process")
    print(f"Output directory: {output_dir}")

    results = []
    for xvg_file in xvg_files:
        print(f"\nProcessing: {xvg_file.name}")

        df, metadata = xp.read_xvg(str(xvg_file))

        if df.empty:
            print(f"  Error: Could not parse data")
            results.append((xvg_file.name, 'error', 'Could not parse data'))
            continue

        # Create plot
        fig = xp.create_plot(df, metadata, roll_avg=args.roll_avg)

        # Save to output
        output_file = output_dir / f"{xvg_file.stem}_plot.html"
        save_figure(fig, str(output_file))

        results.append((xvg_file.name, 'success', metadata.plot_type))

    # Print summary
    print("\n" + "=" * 50)
    print("BATCH PROCESSING SUMMARY")
    print("=" * 50)
    for filename, status, info in results:
        print(f"  {filename}: {status} ({info})")

    return 0


def main() -> int:
    """Main entry point for CLI."""
    parser = argparse.ArgumentParser(
        prog='xvg_plot',
        description='XVG Plot - Automatic plot type detection for GROMACS XVG files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  xvg_plot plot energy.xvg                          # Plot energy components
  xvg_plot plot energy.xvg -o energy.html          # Save to HTML
  xvg_plot plot energy.xvg --roll-avg 20           # With rolling average
  xvg_plot compare ref.xvg sample.xvg --prot1 Ref --prot2 Sample
  xvg_plot multi-compare "*/energy.xvg"            # Compare N files by glob
  xvg_plot multi-compare "run_*/rmsd.xvg" -o out.html
  xvg_plot detect energy.xvg                        # Only detect type
  xvg_plot batch ./prot1/ -o ./output/             # Batch process directory
        '''
    )

    subparsers = parser.add_subparsers(dest='command', help='Available commands')

    # Plot command
    plot_parser = subparsers.add_parser('plot', help='Plot a single XVG file (auto-detects type)')
    plot_parser.add_argument('file', type=str, help='Path to XVG file')
    add_roll_avg_arg(plot_parser, default=10)
    add_output_arg(plot_parser)
    plot_parser.set_defaults(func=cmd_plot)

    # Compare command
    cmp_parser = subparsers.add_parser('compare', help='Compare two XVG files (auto-detects type)')
    cmp_parser.add_argument('file1', type=str, help='Path to first XVG file')
    cmp_parser.add_argument('file2', type=str, help='Path to second XVG file')
    cmp_parser.add_argument('--prot1', type=str, default='prot_1', help='Name for protein 1')
    cmp_parser.add_argument('--prot2', type=str, default='prot_2', help='Name for protein 2')
    add_roll_avg_arg(cmp_parser, default=50)
    add_output_arg(cmp_parser)
    cmp_parser.set_defaults(func=cmd_compare)

    # PCA compare command
    pca_cmp_parser = subparsers.add_parser('pca-compare', help='Compare PCA projection for two proteins')
    pca_cmp_parser.add_argument('file1', type=str, help='Path to first XVG file')
    pca_cmp_parser.add_argument('file2', type=str, help='Path to second XVG file')
    pca_cmp_parser.add_argument('--prot1', type=str, default='prot_1', help='Name for protein 1')
    pca_cmp_parser.add_argument('--prot2', type=str, default='prot_2', help='Name for protein 2')
    add_output_arg(pca_cmp_parser)
    pca_cmp_parser.set_defaults(func=cmd_pca_compare)

    # RMSF compare command
    rmsf_cmp_parser = subparsers.add_parser('rmsf-compare', help='Compare RMSF for two proteins')
    rmsf_cmp_parser.add_argument('file1', type=str, help='Path to first XVG file')
    rmsf_cmp_parser.add_argument('file2', type=str, help='Path to second XVG file')
    rmsf_cmp_parser.add_argument('--prot1', type=str, default='prot_1', help='Name for protein 1')
    rmsf_cmp_parser.add_argument('--prot2', type=str, default='prot_2', help='Name for protein 2')
    add_output_arg(rmsf_cmp_parser)
    rmsf_cmp_parser.set_defaults(func=cmd_rmsf_compare)

    # SASA compare command
    sasa_cmp_parser = subparsers.add_parser('sasa-compare', help='Compare SASA for two proteins')
    sasa_cmp_parser.add_argument('file1', type=str, help='Path to first XVG file')
    sasa_cmp_parser.add_argument('file2', type=str, help='Path to second XVG file')
    sasa_cmp_parser.add_argument('--prot1', type=str, default='prot_1', help='Name for protein 1')
    sasa_cmp_parser.add_argument('--prot2', type=str, default='prot_2', help='Name for protein 2')
    add_roll_avg_arg(sasa_cmp_parser, default=50)
    add_output_arg(sasa_cmp_parser)
    sasa_cmp_parser.set_defaults(func=cmd_sasa_compare)

    # Multi-compare command
    multi_cmp_parser = subparsers.add_parser(
        'multi-compare',
        help='Compare multiple XVG files matching a glob pattern (same type only)')
    multi_cmp_parser.add_argument(
        'pattern', type=str,
        help='Glob pattern for XVG files (e.g. "energy_*.xvg", "*/rmsd.xvg")')
    add_roll_avg_arg(multi_cmp_parser, default=50)
    add_output_arg(multi_cmp_parser)
    multi_cmp_parser.set_defaults(func=cmd_multi_compare)

    # Detect command
    detect_parser = subparsers.add_parser('detect', help='Detect and print plot type only')
    detect_parser.add_argument('file', type=str, help='Path to XVG file')
    detect_parser.set_defaults(func=cmd_detect)

    # Batch command
    batch_parser = subparsers.add_parser('batch', help='Batch process all XVG files in directory')
    batch_parser.add_argument('input_dir', type=str, nargs='?', default='.',
                             help='Input directory with XVG files (default: current directory)')
    batch_parser.add_argument('-o', '--output-dir', type=str, nargs='?', default='.',
                             help='Output directory for plots (default: current directory)')
    add_roll_avg_arg(batch_parser, default=10)
    batch_parser.set_defaults(func=cmd_batch)

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        return 1

    try:
        return args.func(args)
    except FileNotFoundError as e:
        print(f"Error: File not found - {e}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
