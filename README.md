# XVG Analysis

Interactive plotting tool for GROMACS XVG output files with automatic plot type detection, responsive layouts, and publication-ready visualizations. Originally converted from R/ggplot2 to Python/Plotly.

## Installation

### From wheel (recommended)

```bash
pip install dist/xvg_analysis-1.3.0-py3-none-any.whl
```

### Using micromamba/conda

```bash
micromamba env create -f environment.yml
micromamba activate xvg-analysis
```

### Using pip (editable / development)

```bash
python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate
pip install -e .
```

After installation, the `xvg_plot` console utility will be available system-wide.

## Usage

### CLI (Console Utility)

The `xvg_plot` command can be used from anywhere after installation.

```bash
# Show all commands
xvg_plot --help

# Plot a single file (auto-detects type)
xvg_plot plot energy.xvg

# With custom rolling average
xvg_plot plot energy.xvg --roll-avg 50

# Save to custom location
xvg_plot plot energy.xvg -o output.html

# Save to image
xvg_plot plot rmsd.xvg -o rmsd.png

# Detect file type only
xvg_plot detect energy.xvg

# Compare two proteins
xvg_plot compare ref.xvg sample.xvg --prot1 Reference --prot2 Sample

# Compare specific plot types
xvg_plot pca-compare prot1/2dproj.xvg prot2/2dproj.xvg
xvg_plot rmsf-compare prot1/rmsf.xvg prot2/rmsf.xvg
xvg_plot sasa-compare prot1/sasa.xvg prot2/sasa.xvg

# Multi-file comparison (N files by glob pattern, pastel color scheme)
xvg_plot multi-compare "*/energy.xvg"
xvg_plot multi-compare "run_*/rmsd.xvg" -o comparison.html
xvg_plot multi-compare "prot3/rmsf_ch*.xvg" --roll-avg 30

# Batch process directory
xvg_plot batch ./prot1/ -o ./output/
```

### Available Commands

| Command | Description |
|---------|-------------|
| `plot` | Plot a single XVG file (auto-detects type) |
| `compare` | Compare two XVG files |
| `pca-compare` | Compare PCA projection for two proteins |
| `rmsf-compare` | Compare RMSF for two proteins |
| `sasa-compare` | Compare SASA for two proteins |
| `multi-compare` | Compare N XVG files matching a glob pattern (same type only, pastel colors) |
| `detect` | Detect and print plot type only |
| `batch` | Batch process all XVG files in directory |

### Plot Type Detection

Plot type is detected automatically from the `@ title "..."` header line in each XVG file.
Falls back to axis-label / column-count heuristics only when the title does not match any known pattern.

| Plot Type | `@ title` pattern | Description |
|-----------|-------------------|-------------|
| `energy` | "GROMACS Energies", "Energy", ... | Energy components (Potential, Temperature, Pressure, Density) in dynamic grid |
| `gyration` | "Radius of gyration ..." | Radius of gyration (Rg only) |
| `rmsd` | "RMSD" | RMSD of backbone |
| `rmsf` | "RMS fluctuation", "RMSF" | RMSF of Ca atoms |
| `sasa` | "Solvent Accessible Surface", "Area per residue ..." | SASA — time-series or per-residue mode |
| `pca` | "2D projection of trajectory" | PCA projection (2D density contour with viridis) |
| `xy` | Any unrecognized title | Generic XY plot with legend-based column names |

### Features

- **Title-driven plot type detection**: Uses `@ title "..."` from XVG headers to identify energy, gyration, rmsd, rmsf, sasa, pca, and unknown types
- **Dynamic rolling average**: Calculates optimal rolling window size based on 5% of data range, with calculated window size displayed in legend
- **Legend-based column names**: For unfamiliar XVG types, uses `@ s0 legend`, `@ s1 legend`, etc. as column headers
- **Flexible grid layout**: Energy and unknown type plots adapt to number of available columns
- **Responsive plot sizing**: Plots fill browser window - 100% responsive width with autosize enabled
- **Multi-file comparison**: Compare N files matching a glob pattern with pastel color palette (12 distinct pastel colors)
- **Interactive plots**: Uses Plotly for browser-based interactive visualization with zoom, pan, and export capabilities

## Project Structure

```
.
├── xvg_analysis/           # Main package
│   ├── __init__.py
│   ├── xvg_analysis.py    # Core analysis module with plot functions
│   └── xvg_plot.py        # CLI interface
├── dist/                   # Built packages
│   ├── xvg_analysis-*.whl # Python wheel (pip install)
│   └── xvg_analysis-*.tar.gz # Source distribution
├── tests/                  # Test data and automated tests
│   ├── prot1/             # Example XVG files (protein 1)
│   ├── prot2/             # Example XVG files (protein 2)
│   ├── prot3/             # Chain-specific XVG files (chA/chB/chC) for multi-compare
│   ├── other_xvg/         # Additional XVG types for testing
│   ├── results/           # Generated plots (see layout below)
│   ├── test_runner.py     # Automated integration test suite
│   └── test_xvg_analysis.py  # Unit tests
├── environment.yml         # Micromamba/conda environment
├── pyproject.toml         # Package configuration
├── requirements.txt       # Python dependencies
├── CLAUDE.md              # Development guidelines for Claude Code
└── README.md              # This file
```

### Test Results Layout

```
tests/results/
├── prot1/                 # Single-file plots from prot1
├── prot2/                 # Single-file plots from prot2
├── prot3/                 # Single-file plots from prot3
├── other_xvg/             # Single-file plots from other_xvg
├── *_compare.html         # Two-file comparison plots
└── *_multi_compare.html   # Multi-file comparison plots (pastel)
```

## As a Python Module

```python
from xvg_analysis import xvg_analysis as xp

# Read XVG file
df, metadata = xp.read_xvg("energy.xvg")
print(f"Plot type: {metadata.plot_type}")
print(f"Columns: {list(df.columns)}")

# Create plot with auto-calculated rolling average
fig = xp.create_plot(df, metadata)  # roll_avg window size calculated automatically
fig.show()  # Opens in browser

# Or save to file
fig.write_html("output.html")
fig.write_image("output.png")

# Custom rolling average window size
fig = xp.create_plot(df, metadata, roll_avg=25)

# Compare two files
df1, meta1 = xp.read_xvg("prot1/energy.xvg")
df2, meta2 = xp.read_xvg("prot2/energy.xvg")
fig = xp.compare_two_plots(df1, meta1, df2, meta2, prot1="Protein1", prot2="Protein2")
fig.write_html("comparison.html")

# Multi-file comparison (N files, pastel color scheme)
import glob
data_list = []
for f in sorted(glob.glob("prot3/rmsd_ch*.xvg")):
    df, meta = xp.read_xvg(f)
    data_list.append((df, meta, Path(f).stem))
fig = xp.compare_multi_plots(data_list, roll_avg=50)
fig.write_html("multi_comparison.html")
```

## Building

```bash
# Install build tool
pip install build

# Build wheel and sdist
python -m build

# Output in dist/:
#   xvg_analysis-1.3.0-py3-none-any.whl   (wheel)
#   xvg_analysis-1.3.0.tar.gz             (source)
```

The wheel is platform-independent (`py3-none-any`) and can be installed on any system with Python 3.10+.

## Running Tests

### Test Runner (Automated Integration Tests)

```bash
# Run comprehensive test suite (40 tests)
python tests/test_runner.py

# Results organized in tests/results/:
#   prot1/, prot2/, prot3/, other_xvg/  — single-file plots
#   *_compare.html                       — two-file comparisons
#   *_multi_compare.html                 — multi-file comparisons
```

### Unit Tests (pytest)

```bash
# Run all tests
pytest tests/

# Run with verbose output
pytest tests/ -v

# Run specific test class
pytest tests/test_xvg_analysis.py::TestXVGParsing -v

# Run with coverage report
pytest tests/ --cov=xvg_analysis
```
