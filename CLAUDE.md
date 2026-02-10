# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

XVG Analysis is an interactive plotting tool for GROMACS XVG output files. It automatically detects simulation data types and generates publication-ready visualizations using Plotly. Originally converted from R/ggplot2 to Python/Plotly.

## Setup & Installation

### Environment Setup
```bash
# Using Poetry (recommended for development)
poetry install
poetry shell

# Using pip with venv
python -m venv venv
source venv/bin/activate
pip install -e .

# Using micromamba/conda
micromamba env create -f environment.yml
micromamba activate xvg-analysis
```

## Development Commands

### Run the CLI tool
```bash
# After installation, available globally as:
xvg_plot [command] [options]

# Or run directly:
python -m xvg_analysis.xvg_plot [command] [options]
```

### Code Quality
```bash
# Format code (black)
black xvg_analysis/

# Lint code (flake8)
flake8 xvg_analysis/

# Type checking (mypy)
mypy xvg_analysis/
```

## Architecture

### Core Modules

**`xvg_analysis/xvg_analysis.py`** (main analysis engine):
- `XVGMetadata`: Dataclass storing parsed metadata (title, labels, units, plot_type, legend_labels)
- `parse_xvg_metadata()`: Reads XVG file, extracts headers and metadata
- `parse_data_to_dataframe()`: Converts numeric data lines to pandas DataFrame with intelligent column naming
- `classify_plot_type()`: Detects plot type from metadata (energy, gyration, rmsd, rmsf, sasa, pca, xy)
- `read_xvg()`: Main entry point - returns DataFrame and XVGMetadata
- `calculate_roll_avg()`: Computes optimal rolling window (5% of data range)
- `create_plot()`: Generates Plotly figure based on detected type
- `compare_two_plots()`: Creates side-by-side comparison plots
- Type-specific comparison functions: `compare_pca_plots()`, `compare_rmsf_plots()`, `compare_sasa_plots()`

**`xvg_analysis/xvg_plot.py`** (CLI interface):
- Argparse-based command dispatcher with multiple subcommands
- Command handlers: `cmd_plot()`, `cmd_compare()`, `cmd_detect()`, `cmd_batch()`, `cmd_pca_compare()`, `cmd_rmsf_compare()`, `cmd_sasa_compare()`
- Helper functions: `save_figure()`, `open_in_browser()`, shared argument builders (`add_output_arg()`, `add_roll_avg_arg()`)

### Data Flow

1. **Parsing**: XVG file → metadata extraction (headers) + data lines
2. **Type Detection**: Metadata analysis → plot type classification
3. **Naming**: Smart column naming based on title/labels (Time, Residue, PC1, PC2, etc.)
4. **Plotting**: DataFrame + metadata → Plotly figure with type-specific layout
5. **Output**: HTML (interactive via Plotly), PNG/SVG/PDF (static via Plotly)

### Plot Types & Detection

| Type | Detection Keys | Typical Columns | Special Features |
|------|---|---|---|
| `energy` | "energy"/"potential" + 4+ cols + time | Time, Potential, Temperature, Pressure, Density | Dynamic grid layout, 2-column subplots |
| `gyration` | "gyration"/"rg" | Time, Rg, [RgX, RgY, RgZ] | Single or multi-dimensional Rg |
| `rmsd` | "rmsd" + time | Time, RMSD | Standard deviation handling |
| `rmsf` | "rmsf"/"fluctuation" + residue | Residue, RMSF | Per-residue flexibility analysis |
| `sasa` | "sasa"/"solvent"/"area" | Time/Residue, Area, [StdDev] | Time-series or per-residue modes |
| `pca` | "projection" + 2 cols | PC1, PC2 | 2D density contour with viridis colormap |
| `xy` | Unknown types | Auto-detected from legend | Generic fallback with legend-based names |

### Key Design Patterns

- **Auto-detection**: No type specification needed - metadata analysis determines plot strategy
- **Legend-based naming**: For unknown types, uses `@ s0 legend`, `@ s1 legend` entries as column headers
- **Dynamic layout**: Energy plots adapt subgrid to column count; width scales with data range
- **Rolling average**: Optional smoothing with auto-calculation (5% of data range when not specified)
- **Flexible output**: Interactive HTML (default) or static PNG/SVG/PDF (via Plotly renderers)

## CLI Commands

### Single File Plotting
- `plot [file] [--roll-avg N] [-o output]`: Plot single XVG with optional smoothing
- `detect [file]`: Show detected type and metadata without plotting

### Comparisons
- `compare [file1] [file2] [--prot1 name] [--prot2 name] [--roll-avg N] [-o output]`: Auto-detect type comparison
- `pca-compare [file1] [file2] [--prot1 name] [--prot2 name] [-o output]`: 2D density overlay
- `rmsf-compare [file1] [file2] [--prot1 name] [--prot2 name] [-o output]`: Per-residue comparison
- `sasa-compare [file1] [file2] [--prot1 name] [--prot2 name] [--roll-avg N] [-o output]`: SASA comparison

### Batch Processing
- `batch [input_dir] [-o output_dir] [--roll-avg N]`: Process all XVG files in directory

## Common Development Tasks

### Adding Support for a New Plot Type

1. **Update `classify_plot_type()`** in xvg_analysis.py with detection keywords
2. **Add column naming logic** in `parse_data_to_dataframe()` (look for title/label patterns)
3. **Create plot generation function** (e.g., `create_energy_plot()`) using plotly subplots/go
4. **Update `create_plot()`** to dispatch to new function based on `metadata.plot_type`
5. **Add comparison function** if needed (e.g., `compare_xy_plots()`)
6. **Test** with sample XVG files in `tests/` directory

### Debugging Type Detection Issues

- Use `xvg_plot detect [file]` to inspect: detected type, parsed title, axis labels, column count, column names
- Check raw file headers: `head -20 [file.xvg]` for metadata patterns
- Verify detection logic in `classify_plot_type()` - order matters (energy checked before gyration, SASA before RMSD)

### Testing New Features

Test data is in `tests/` with sample files:
- `tests/prot1/` and `tests/prot2/`: Complete protein simulation outputs (energy, gyration, rmsd, rmsf, sasa, 2dproj)
- `tests/other_xvg/`: Additional unfamiliar formats for edge case testing
- `tests/results/`: Generated plots for visual inspection

## Dependencies

- **Core**: numpy, pandas, plotly
- **Dev**: black (formatting), flake8 (linting), mypy (type checking)
- **Python**: 3.10+

## Important Notes

- Type annotations are strict (mypy config enforces disallow_untyped_defs) - all functions must have full signatures
- Rolling average calculation: 5% of (max - min) of data, with minimum of 1
- Plotly figures use CDN-hosted js for HTML to keep file sizes smaller
- Legend labels from XVG `@ s0 legend` directives are preferred over auto-generated names
- Energy plots use dynamic 2-column grid for optimal visualization regardless of column count
