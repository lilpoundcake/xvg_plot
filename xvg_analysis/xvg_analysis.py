"""
XVG analysis module for GROMACS output files.
Provides automatic plot type detection and parsing of XVG files.
"""

import re
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from pathlib import Path
from typing import Optional, Dict, List, Tuple, Any
from dataclasses import dataclass


@dataclass
class XVGMetadata:
    """Metadata extracted from XVG file headers."""
    title: str = ""
    x_label: str = ""
    y_label: str = ""
    plot_type: str = ""  # energy, gyration, rmsd, rmsf, sasa, pca
    x_unit: str = ""
    y_unit: str = ""
    legend_labels: List[str] = None
    n_columns: int = 0
    file_path: str = ""
    command: str = ""  # Original GROMACS command
    gmx_version: str = ""

    def __post_init__(self):
        if self.legend_labels is None:
            self.legend_labels = []


def parse_xvg_metadata(file_path: str) -> XVGMetadata:
    """
    Parse XVG file and extract metadata including plot type, axis labels, and units.
    Also reads all numeric data into a DataFrame.
    """
    metadata = XVGMetadata(file_path=file_path)
    data_lines = []
    legend_labels = {}
    gmx_version = ""
    command = ""

    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()

            # Skip empty lines
            if not line:
                continue

            # Parse GROMACS header comments
            if line.startswith('#'):
                if 'GROMACS' in line:
                    match = re.search(r'GROMACS - gmx (\w+), (\d+)', line)
                    if match:
                        command = match.group(1)
                        gmx_version = match.group(2)
                continue

            # Parse @ directives (normalize spaces after @)
            if line.startswith('@'):
                # Normalize: remove leading @ and extra spaces
                directive = line[1:].strip()

                # Title
                if directive.startswith('title'):
                    match = re.search(r'"([^"]*)"', line)
                    if match:
                        metadata.title = match.group(1)

                # X axis label (handle variable spacing: xaxis  label or xaxis label)
                elif re.match(r'^xaxis\s+label\s*', directive):
                    match = re.search(r'"([^"]*)"', line)
                    if match:
                        metadata.x_label = match.group(1)

                # Y axis label (handle variable spacing: yaxis  label or yaxis label)
                elif re.match(r'^yaxis\s+label\s*', directive):
                    match = re.search(r'"([^"]*)"', line)
                    if match:
                        metadata.y_label = match.group(1)

                # Legend labels (s0, s1, etc.)
                elif directive.startswith('s') and 'legend' in directive:
                    match = re.search(r's(\d+)\s+legend\s+"([^"]+)"', line)
                    if match:
                        idx = int(match.group(1))
                        label = match.group(2)
                        legend_labels[idx] = label

                # Command line info
                elif 'command' in directive.lower():
                    if 'gmx' in line:
                        match = re.search(r'gmx (\w+)', line)
                        if match:
                            command = match.group(1)

                # Store version if found
                elif 'GROMACS' in line and '20' in line:
                    match = re.search(r'GROMACS - gmx (\w+), (\d+)', line)
                    if match:
                        command = match.group(1)
                        gmx_version = match.group(2)

            # Parse numeric data lines
            elif not line.startswith(('@', '#', '&')):
                # Check if line contains numeric data
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        # Validate first two columns are numeric
                        float(parts[0])
                        float(parts[1])
                        data_lines.append(line)
                    except ValueError:
                        continue

    metadata.legend_labels = [legend_labels.get(i, f"Col{i}") for i in range(len(legend_labels))]
    metadata.command = command
    metadata.gmx_version = gmx_version

    # Parse numeric data into DataFrame
    if data_lines:
        metadata = parse_data_to_dataframe(data_lines, metadata)

    # Classify plot type based on title and labels
    metadata.plot_type = classify_plot_type(metadata)

    # Ensure df attribute exists
    if not hasattr(metadata, 'df'):
        metadata.df = pd.DataFrame()

    return metadata


def parse_data_to_dataframe(data_lines: List[str], metadata: XVGMetadata) -> XVGMetadata:
    """Parse data lines into DataFrame and set column names."""
    parsed_data = []

    for line in data_lines:
        # Split by whitespace, comma, or tab
        parts = re.split(r'[,\s\t]+', line)
        parsed_data.append([p for p in parts if p])

    if not parsed_data:
        metadata.df = pd.DataFrame()
        return metadata

    # Determine number of columns from first row
    num_cols = len(parsed_data[0])
    metadata.n_columns = num_cols

    # Try to detect specific column types from title or labels
    title_lower = metadata.title.lower()
    x_label_lower = metadata.x_label.lower()
    y_label = metadata.y_label

    # Create initial column names
    if num_cols >= 2:
        columns = ['X'] + [f'Y{i}' for i in range(1, num_cols)]

    # Try to extract meaningful names from Y label (only if not energy file)
    if y_label and 'energy' not in title_lower:
        # Check for multiple units in parentheses (like energy files)
        if '),' in y_label:
            unit_matches = re.findall(r'\(([^)]+)\)', y_label)
            if len(unit_matches) >= num_cols - 1:
                columns = ['X'] + [f"{m.strip()}" for m in unit_matches[:num_cols-1]]

    # Set X column name based on context
    # Check time first (most common for simulation data)
    if 'time' in x_label_lower or 'time' in title_lower:
        columns[0] = 'Time'
    elif 'step' in x_label_lower:
        columns[0] = 'Step'
    elif 'residue' in x_label_lower or re.search(r'\bresidue\b', title_lower):
        columns[0] = 'Residue'
    elif 'atom' in x_label_lower:
        columns[0] = 'Atom'
    elif 'eigen' in x_label_lower or 'eigenvector' in x_label_lower:
        columns[0] = 'Eigenvalue'
    elif 'projection' in x_label_lower or 'projection' in title_lower:
        if 'eigenvector' in x_label_lower:
            columns[0] = 'PC1'
            columns[1] = 'PC2'

    # Set Y column names based on context
    # Energy check - handle "energy", "energies", "GROMACS Energies"
    if 'potential' in title_lower or 'energ' in title_lower:
        energy_cols = ['Potential', 'Temperature', 'Pressure', 'Density']
        columns = ['Time'] + energy_cols[:num_cols-1]
    elif re.search(r'\bgyration\b', title_lower) or re.search(r'\brg\b', title_lower):
        if num_cols >= 4:
            columns = ['Time', 'Rg', 'RgX', 'RgY', 'RgZ']
        else:
            columns = ['Time', 'Rg']
    # Check SASA before RMSD since both have similar units
    elif 'area' in y_label.lower() or re.search(r'\bsasa\b', title_lower) or 'solvent' in title_lower:
        # SASA files can have Area + Standard deviation columns
        # Use Time for SASA (most common), Residue for surface per residue analysis
        if 'time' in x_label_lower:
            columns = ['Time', 'Area', 'StdDev'][:num_cols]
        else:
            columns = ['Residue', 'Area', 'StdDev'][:num_cols]
    elif re.search(r'\brmsd\b', title_lower) or 'rmsd' in y_label.lower():
        columns = ['Time', 'RMSD']
    elif re.search(r'\brmsf\b', title_lower) or 'fluctuation' in title_lower:
        if 'residue' in x_label_lower:
            columns = ['Residue', 'RMSF']
        else:
            columns = ['Residue', 'RMSF']
    elif 'pca' in title_lower or 'projection' in title_lower:
        columns = ['PC1', 'PC2']

    # Convert to DataFrame
    df = pd.DataFrame(parsed_data, columns=columns[:num_cols])

    # Convert all columns to numeric
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    metadata.df = df
    return metadata


def classify_plot_type(metadata: XVGMetadata) -> str:
    """
    Classify the plot type based on the @ title directive from the XVG file.

    Title-driven detection using the '@ title "..."' header line.
    Falls back to axis label / column count heuristics only when the title
    does not match any known pattern.
    """
    title = metadata.title.lower()

    # ── Title-based detection (primary) ────────────────────────────────
    # Title patterns ordered from most specific to least specific.
    # Each entry: (title substring, plot type)
    title_rules: List[Tuple[str, str]] = [
        # PCA / 2D projection
        ('2d projection', 'pca'),
        ('projection of trajectory', 'pca'),
        # Energy
        ('energies', 'energy'),
        ('energy', 'energy'),
        # SASA (check before RMSD — "area" could appear in RMSD titles)
        ('solvent accessible surface', 'sasa'),
        ('area per residue', 'sasa'),
        ('sasa', 'sasa'),
        # RMSF (check before RMSD — "rms fluctuation" contains "rms")
        ('rms fluctuation', 'rmsf'),
        ('rmsf', 'rmsf'),
        # RMSD
        ('rmsd', 'rmsd'),
        # Gyration
        ('radius of gyration', 'gyration'),
        ('gyration', 'gyration'),
    ]

    for pattern, plot_type in title_rules:
        if pattern in title:
            return plot_type

    # ── Fallback: axis-label / column-count heuristics ─────────────────
    x_label = metadata.x_label.lower()
    y_label = metadata.y_label.lower()
    n_cols = metadata.n_columns

    # Energy by units (no recognizable title)
    if n_cols >= 4 and 'time' in x_label:
        if any(kw in y_label for kw in ['kj/mol', 'k), (bar', 'kg/m']):
            return 'energy'

    # PCA by axis labels
    if n_cols == 2 and 'eigenvector' in x_label:
        return 'pca'

    # Default: xy for any file with 2+ columns
    if n_cols >= 2:
        return 'xy'

    return 'unknown'


def calculate_roll_avg(df: pd.DataFrame, x_col: str = 'Time') -> int:
    """
    Calculate rolling average window size as 5% of data range.

    Args:
        df: DataFrame with data
        x_col: X-axis column name (default: 'Time')

    Returns:
        Rolling average window size (minimum 3)
    """
    if x_col not in df.columns or len(df) < 10:
        return 10  # Default fallback

    x_range = df[x_col].max() - df[x_col].min()
    # 5% of range, minimum 3, maximum 100
    roll_avg = max(3, min(100, max(3, int(x_range * 0.05))))
    return roll_avg


def read_xvg(file_path: str) -> Tuple[pd.DataFrame, XVGMetadata]:
    """
    Read an XVG file and return DataFrame with metadata.

    Args:
        file_path: Path to XVG file

    Returns:
        Tuple of (DataFrame, XVGMetadata)
    """
    metadata = parse_xvg_metadata(file_path)
    return metadata.df, metadata


# Plotting functions
def plot_energy(df: pd.DataFrame, metadata: XVGMetadata, roll_avg: Optional[int] = None) -> go.Figure:
    """
    Plot energy components with rolling average.
    Dynamically adjusts grid based on available columns.
    """
    # Calculate roll_avg from data range if not provided (default)
    if roll_avg is None:
        roll_avg = calculate_roll_avg(df, 'Time')

    # Standard energy columns in preferred order
    energy_cols = ['Potential', 'Temperature', 'Pressure', 'Density']

    # Filter to only columns that exist in dataframe
    available_cols = [col for col in energy_cols if col in df.columns]

    # Add any additional columns (non-energy) that weren't matched
    additional_cols = [col for col in df.columns if col not in ['Time'] + energy_cols]
    all_y_cols = available_cols + additional_cols

    n_cols = len(all_y_cols)
    if n_cols == 0:
        # Fallback if no columns found
        return plot_xy(df, metadata, roll_avg)

    # Calculate grid dimensions
    n_rows = (n_cols + 1) // 2 if n_cols > 2 else 1
    if n_cols > 4:
        n_rows = (n_cols + 2) // 3  # 3 columns for many items
        n_cols = 3
    else:
        n_cols = min(2, n_cols) if n_cols > 1 else 1

    # Create subplots with dynamic grid
    fig = make_subplots(rows=n_rows, cols=n_cols,
                        subplot_titles=all_y_cols)

    # Map y_columns to titles
    y_titles_map = {
        'Pressure': 'Pressure (bar)',
        'Temperature': 'Temperature (K)',
        'Potential': 'Potential (kJ/mol)',
        'Density': 'Density (kg/m^3)'
    }

    for idx, y_col in enumerate(all_y_cols):
        row = (idx // n_cols) + 1
        col = (idx % n_cols) + 1

        if y_col in df.columns:
            y_data = df[y_col]
            y_ra = y_data.rolling(window=roll_avg, min_periods=1, center=True).mean()

            fig.add_trace(
                go.Scatter(x=df['Time'], y=y_data, mode='lines', name=y_col,
                          line=dict(color='lightblue'), opacity=0.7),
                row=row, col=col
            )
            fig.add_trace(
                go.Scatter(x=df['Time'], y=y_ra, mode='lines', name=f'RA (n={roll_avg})',
                          line=dict(color='navy')),
                row=row, col=col
            )
            y_title = y_titles_map.get(y_col, y_col)
            fig.update_yaxes(title_text=y_title, row=row, col=col)

    fig.update_xaxes(title_text='Time (ps)', row=n_rows, col=1)
    if n_cols > 1 and n_rows > 1:
        fig.update_xaxes(title_text='Time (ps)', row=n_rows, col=n_cols)

    fig.update_layout(
        title_text='Energy Components',
        height=400 + (n_rows - 1) * 500,
        autosize=True
    )
    return fig


def plot_gyration(df: pd.DataFrame, metadata: XVGMetadata, roll_avg: Optional[int] = None) -> go.Figure:
    """Plot radius of gyration with rolling average (Rg only, no RgX/Y/Z)."""
    # Calculate roll_avg from data range if not provided
    if roll_avg is None:
        roll_avg = calculate_roll_avg(df, 'Time')

    fig = go.Figure()

    # Plot Rg (first Y column) with RA only
    if 'Rg' in df.columns:
        y_data = df['Rg']
        y_ra = y_data.rolling(window=roll_avg, min_periods=1, center=True).mean()
        fig.add_trace(go.Scatter(x=df['Time'], y=y_data, mode='lines',
                                name='Rg', line=dict(color='lightblue'), opacity=0.7))
        fig.add_trace(go.Scatter(x=df['Time'], y=y_ra, mode='lines',
                                name=f'RA (n={roll_avg})', line=dict(color='navy')))

    fig.update_layout(
        title=f'Radius of gyration - {metadata.title}',
        xaxis_title='Time (ps)',
        yaxis_title='Rg (nm)',
        template='plotly_white',
        autosize=True
    )
    return fig


def plot_rmsd(df: pd.DataFrame, metadata: XVGMetadata, roll_avg: Optional[int] = None) -> go.Figure:
    """Plot RMSD with rolling average."""
    # Calculate roll_avg from data range if not provided
    if roll_avg is None:
        roll_avg = calculate_roll_avg(df, 'Time')

    fig = go.Figure()

    if 'RMSD' in df.columns:
        y_data = df['RMSD']
        y_ra = y_data.rolling(window=roll_avg, min_periods=1, center=True).mean()
        fig.add_trace(go.Scatter(x=df['Time'], y=y_data, mode='lines',
                                name='RMSD', line=dict(color='lightblue'), opacity=0.7))
        fig.add_trace(go.Scatter(x=df['Time'], y=y_ra, mode='lines',
                                name=f'RA (n={roll_avg})', line=dict(color='navy')))

    fig.update_layout(
        title=f'RMSD of backbone - {metadata.title}',
        xaxis_title='Time (ns)',
        yaxis_title='RMSD (nm)',
        template='plotly_white',
        autosize=True
    )
    return fig


def plot_rmsf(df: pd.DataFrame, metadata: XVGMetadata) -> go.Figure:
    """Plot RMSF of Ca atoms."""
    fig = go.Figure()

    if 'RMSF' in df.columns:
        fig.add_trace(go.Scatter(x=df['Residue'], y=df['RMSF'], mode='lines',
                                name='RMSF', line=dict(color='lightblue')))

    fig.update_layout(
        title=f'RMSF of Ca atoms - {metadata.title}',
        xaxis_title='Residue',
        yaxis_title='RMSF (nm)',
        template='plotly_white',
        autosize=True
    )
    return fig


def plot_sasa(df: pd.DataFrame, metadata: XVGMetadata, roll_avg: Optional[int] = None) -> go.Figure:
    """Plot SASA — handles both time-series and per-residue modes."""
    fig = go.Figure()

    # Determine X column: Time (time-series) or Residue (per-residue)
    if 'Time' in df.columns:
        x_col = 'Time'
        x_title = 'Time (ns)'
    elif 'Residue' in df.columns:
        x_col = 'Residue'
        x_title = 'Residue'
    else:
        x_col = df.columns[0]
        x_title = x_col

    if 'Area' in df.columns:
        y_data = df['Area']

        if x_col == 'Residue':
            # Per-residue: no rolling average, show as line
            fig.add_trace(go.Scatter(x=df[x_col], y=y_data, mode='lines',
                                    name='Area', line=dict(color='lightblue')))
        else:
            # Time-series: with rolling average
            if roll_avg is None:
                roll_avg = calculate_roll_avg(df, x_col)
            y_ra = y_data.rolling(window=roll_avg, min_periods=1, center=True).mean()
            fig.add_trace(go.Scatter(x=df[x_col], y=y_data, mode='lines',
                                    name='Area', line=dict(color='lightblue'), opacity=0.7))
            fig.add_trace(go.Scatter(x=df[x_col], y=y_ra, mode='lines',
                                    name=f'RA (n={roll_avg})', line=dict(color='navy')))

    fig.update_layout(
        title=f'Solvent Accessible Surface Area - {metadata.title}',
        xaxis_title=x_title,
        yaxis_title='Area (nm^2)',
        template='plotly_white',
        autosize=True
    )
    return fig


def plot_pca(df: pd.DataFrame, metadata: XVGMetadata) -> go.Figure:
    """Plot PCA projection with 2D density contour using viridis coloring."""
    fig = go.Figure()

    if 'PC1' in df.columns and 'PC2' in df.columns:
        fig.add_trace(go.Histogram2dContour(
            x=df['PC1'], y=df['PC2'],
            name='PCA',
            colorscale='Viridis',
            showscale=False
        ))

    fig.update_layout(
        title=f'PCA Projection - {metadata.title}',
        xaxis_title='PC1',
        yaxis_title='PC2',
        template='plotly_white',
        autosize=True
    )
    return fig


def plot_eigenval(df: pd.DataFrame, metadata: XVGMetadata) -> go.Figure:
    """Plot eigenvalues from PCA/covariance analysis."""
    fig = go.Figure()

    x_col = 'Eigenvalue' if 'Eigenvalue' in df.columns else df.columns[0]
    y_col = 'Value' if 'Value' in df.columns else df.columns[1]

    fig.add_trace(go.Scatter(x=df[x_col], y=df[y_col], mode='lines+markers',
                            name='Eigenvalue', line=dict(color='navy')))

    fig.update_layout(
        title=f'Eigenvalues - {metadata.title}',
        xaxis_title='Eigenvalue Index',
        yaxis_title='Value (nm^2)',
        template='plotly_white',
        autosize=True
    )
    return fig


def plot_eigrmsf(df: pd.DataFrame, metadata: XVGMetadata) -> go.Figure:
    """Plot eigrmsf data (RMSF with eigenvalue contribution)."""
    fig = go.Figure()

    if 'Residue' in df.columns and 'RMSF' in df.columns:
        fig.add_trace(go.Scatter(x=df['Residue'], y=df['RMSF'], mode='lines',
                                name='RMSF', line=dict(color='lightblue')))

    fig.update_layout(
        title=f'Eigenvector RMSF - {metadata.title}',
        xaxis_title='Residue',
        yaxis_title='RMSF (nm)',
        template='plotly_white',
        autosize=True
    )
    return fig


def plot_rmsf_compare(df1: pd.DataFrame, df2: pd.DataFrame,
                      prot1: str = 'prot_1', prot2: str = 'prot_2') -> go.Figure:
    """Create comparison plot for RMSF (RMS fluctuation)."""
    fig = go.Figure()

    if 'RMSF' in df1.columns:
        fig.add_trace(go.Scatter(x=df1['Residue'], y=df1['RMSF'], mode='lines',
                                name=prot1, line=dict(color='lightblue')))
        fig.add_trace(go.Scatter(x=df2['Residue'], y=df2['RMSF'], mode='lines',
                                name=prot2, line=dict(color='#e84f82')))

    fig.update_layout(
        title='RMSF of Ca atoms',
        xaxis_title='Residue',
        yaxis_title='RMSF (nm)',
        template='plotly_white',
        legend_title='Protein',
        autosize=True
    )
    return fig


def plot_sasa_compare(df1: pd.DataFrame, df2: pd.DataFrame,
                      prot1: str = 'prot_1', prot2: str = 'prot_2',
                      roll_avg: int = 50) -> go.Figure:
    """Create comparison plot for SASA (Solvent Accessible Surface)."""
    fig = go.Figure()

    if 'Area' in df1.columns:
        y_ra1 = df1['Area'].rolling(window=roll_avg, min_periods=1, center=True).mean()
        y_ra2 = df2['Area'].rolling(window=roll_avg, min_periods=1, center=True).mean()

        fig.add_trace(go.Scatter(x=df1['Time'], y=y_ra1, mode='lines',
                                name=f'{prot1} SASA', line=dict(color='lightblue')))
        fig.add_trace(go.Scatter(x=df2['Time'], y=y_ra2, mode='lines',
                                name=f'{prot2} SASA', line=dict(color='#e84f82')))

    fig.update_layout(
        title='Solvent Accessible Surface Area',
        xaxis_title='Time (ns)',
        yaxis_title='Area (nm^2)',
        template='plotly_white',
        legend_title='Protein',
        autosize=True
    )
    return fig


def plot_xy(df: pd.DataFrame, metadata: XVGMetadata, roll_avg: Optional[int] = None) -> go.Figure:
    """Generic XY plot."""
    # Calculate roll_avg from data range if not provided
    if roll_avg is None:
        roll_avg = calculate_roll_avg(df, df.columns[0])

    fig = go.Figure()

    x_col = df.columns[0]
    y_col = df.columns[1] if len(df.columns) > 1 else 'Y'

    if x_col in df.columns and y_col in df.columns:
        y_data = df[y_col]
        y_ra = y_data.rolling(window=roll_avg, min_periods=1, center=True).mean()
        fig.add_trace(go.Scatter(x=df[x_col], y=y_data, mode='lines',
                                name=y_col, line=dict(color='lightblue'), opacity=0.7))
        fig.add_trace(go.Scatter(x=df[x_col], y=y_ra, mode='lines',
                                name=f'RA (n={roll_avg})', line=dict(color='navy')))

    fig.update_layout(
        title=f'Plot - {metadata.title}',
        xaxis_title=metadata.x_label or x_col,
        yaxis_title=metadata.y_label or y_col,
        template='plotly_white',
        autosize=True
    )
    return fig


def create_flexible_plot(df: pd.DataFrame, metadata: XVGMetadata, roll_avg: Optional[int] = None) -> go.Figure:
    """
    Create a flexible plot for unfamiliar XVG types.
    Automatically detects columns and creates appropriate subplots based on metadata.
    Uses legend labels from XVG file as column names.
    """
    # Calculate roll_avg from data range if not provided
    if roll_avg is None:
        roll_avg = calculate_roll_avg(df, df.columns[0])

    fig = make_subplots(rows=1, cols=1)

    title = metadata.title or "XVG Plot"
    x_label = metadata.x_label if metadata.x_label else (df.columns[0] if len(df.columns) > 0 else "X")
    y_label = metadata.y_label if metadata.y_label else "Y"

    # Get legend labels for column names (priority)
    y_column_names = []
    for idx, label in enumerate(metadata.legend_labels):
        # Clean up legend label - extract meaningful name from "Atom Collection ASPP ASPP, index (1)"
        clean_label = label.strip()
        # Remove common prefixes/suffixes
        if clean_label.startswith('"') and clean_label.endswith('"'):
            clean_label = clean_label[1:-1]
        if clean_label:
            y_column_names.append(clean_label)

    # Get all columns except the first one (assumed to be X)
    y_columns = list(df.columns[1:]) if len(df.columns) > 1 else []

    # Use legend labels if we have matching number of columns
    if y_column_names and len(y_column_names) >= len(y_columns):
        y_column_names = y_column_names[:len(y_columns)]
    else:
        y_column_names = y_columns

    # Default height for single plot (consistent with subplot formula base)
    height = 600

    if not y_columns:
        # Single XY line
        x_col = df.columns[0]
        fig.add_trace(go.Scatter(x=df[x_col], y=df[x_col], mode='lines', name=x_col))
    else:
        # Multiple Y columns - create separate subplots in a grid
        n_y_cols = len(y_columns)
        n_cols = min(3, max(1, n_y_cols))
        n_rows = (n_y_cols + n_cols - 1) // n_cols

        fig = make_subplots(rows=n_rows, cols=n_cols, subplot_titles=y_column_names)

        for idx, (y_col, col_name) in enumerate(zip(y_columns, y_column_names)):
            row = (idx // n_cols) + 1
            col = (idx % n_cols) + 1

            if y_col in df.columns:
                y_data = df[y_col]
                y_ra = y_data.rolling(window=roll_avg, min_periods=1, center=True).mean()

                fig.add_trace(
                    go.Scatter(x=df[df.columns[0]], y=y_data, mode='lines',
                              name=col_name, line=dict(color='lightblue'), opacity=0.7),
                    row=row, col=col
                )
                fig.add_trace(
                    go.Scatter(x=df[df.columns[0]], y=y_ra, mode='lines',
                              name=f'RA (n={roll_avg})', line=dict(color='navy')),
                    row=row, col=col
                )

        fig.update_xaxes(title_text=df.columns[0], row=n_rows, col=1)
        if n_cols > 1:
            fig.update_xaxes(title_text=df.columns[0], row=n_rows, col=n_cols)

        # Calculate height based on number of rows (minimum 600px per row)
        # Add extra space for title and labels
        height = 600 + (n_rows - 1) * 250

    fig.update_layout(autosize=True, height=height)
    return fig


def create_plot(df: pd.DataFrame, metadata: XVGMetadata, roll_avg: Optional[int] = None) -> go.Figure:
    """
    Create appropriate plot based on detected plot type.

    Args:
        df: DataFrame with data
        metadata: XVGMetadata with plot information
        roll_avg: Rolling average window size

    Returns:
        Plotly figure object
    """
    plot_type = metadata.plot_type

    if plot_type == 'energy':
        return plot_energy(df, metadata, roll_avg)
    elif plot_type == 'gyration':
        return plot_gyration(df, metadata, roll_avg)
    elif plot_type == 'rmsd':
        return plot_rmsd(df, metadata, roll_avg)
    elif plot_type == 'rmsf':
        return plot_rmsf(df, metadata)
    elif plot_type == 'sasa':
        return plot_sasa(df, metadata, roll_avg)
    elif plot_type == 'pca':
        return plot_pca(df, metadata)
    elif plot_type == 'unknown' or plot_type == 'xy':
        return create_flexible_plot(df, metadata, roll_avg)
    else:
        return plot_xy(df, metadata, roll_avg)


def compare_energy_plots(df1: pd.DataFrame, df2: pd.DataFrame,
                        meta1: XVGMetadata, meta2: XVGMetadata,
                        prot1: str = 'prot_1', prot2: str = 'prot_2',
                        roll_avg: int = 50) -> go.Figure:
    """Create energy comparison plot with dynamic grid layout.
    Each subplot shows both proteins with different colors."""
    # Standard energy columns in preferred order
    energy_cols = ['Potential', 'Temperature', 'Pressure', 'Density']

    # Filter to only columns that exist in both dataframes
    available_cols = [col for col in energy_cols if col in df1.columns and col in df2.columns]

    # Add any additional columns that exist in both
    additional_cols = [col for col in df1.columns
                       if col not in ['Time'] + energy_cols and col in df2.columns]
    all_y_cols = available_cols + additional_cols

    n_cols = len(all_y_cols)
    if n_cols == 0:
        return plot_xy(df1, meta1, roll_avg)

    # Calculate grid dimensions
    n_rows = (n_cols + 1) // 2 if n_cols > 2 else 1
    if n_cols > 4:
        n_rows = (n_cols + 2) // 3
        n_cols = 3
    else:
        n_cols = min(2, n_cols) if n_cols > 1 else 1

    # Create subplots with dynamic grid
    fig = make_subplots(rows=n_rows, cols=n_cols,
                        subplot_titles=all_y_cols)

    y_titles_map = {
        'Pressure': 'Pressure (bar)',
        'Temperature': 'Temperature (K)',
        'Potential': 'Potential (kJ/mol)',
        'Density': 'Density (kg/m^3)'
    }

    for idx, y_col in enumerate(all_y_cols):
        row = (idx // n_cols) + 1
        col = (idx % n_cols) + 1

        if y_col in df1.columns and y_col in df2.columns:
            y_ra1 = df1[y_col].rolling(window=roll_avg, min_periods=1, center=True).mean()
            y_ra2 = df2[y_col].rolling(window=roll_avg, min_periods=1, center=True).mean()

            fig.add_trace(
                go.Scatter(x=df1['Time'], y=y_ra1, mode='lines',
                          name=prot1, line=dict(color='lightblue')),
                row=row, col=col
            )
            fig.add_trace(
                go.Scatter(x=df2['Time'], y=y_ra2, mode='lines',
                          name=prot2, line=dict(color='#e84f82')),
                row=row, col=col
            )
            y_title = y_titles_map.get(y_col, y_col)
            fig.update_yaxes(title_text=y_title, row=row, col=col)

    fig.update_xaxes(title_text='Time (ps)', row=n_rows, col=1)
    if n_cols > 1 and n_rows > 1:
        fig.update_xaxes(title_text='Time (ps)', row=n_rows, col=n_cols)

    fig.update_layout(
        title_text='Energy Components Comparison',
        height=400 + (n_rows - 1) * 500,
        autosize=True
    )
    return fig


def compare_two_plots(df1: pd.DataFrame, meta1: XVGMetadata,
                      df2: pd.DataFrame, meta2: XVGMetadata,
                      prot1: str = 'prot_1', prot2: str = 'prot_2',
                      roll_avg: int = 50) -> go.Figure:
    """
    Create comparison plot for two XVG files.
    """
    plot_type = meta1.plot_type

    if plot_type == 'energy':
        return compare_energy_plots(df1, df2, meta1, meta2, prot1, prot2, roll_avg)

    fig = go.Figure()

    if plot_type == 'gyration':
        for col in ['Rg', 'RgX', 'RgY', 'RgZ']:
            if col in df1.columns:
                y_ra1 = df1[col].rolling(window=roll_avg, min_periods=1, center=True).mean()
                y_ra2 = df2[col].rolling(window=roll_avg, min_periods=1, center=True).mean()

                fig.add_trace(go.Scatter(x=df1['Time'], y=y_ra1, mode='lines',
                                        name=f'{prot1} {col} RA', line=dict(color='lightblue')))
                fig.add_trace(go.Scatter(x=df2['Time'], y=y_ra2, mode='lines',
                                        name=f'{prot2} {col} RA', line=dict(color='#e84f82')))

    elif plot_type == 'rmsd':
        if 'RMSD' in df1.columns:
            y_ra1 = df1['RMSD'].rolling(window=roll_avg, min_periods=1, center=True).mean()
            y_ra2 = df2['RMSD'].rolling(window=roll_avg, min_periods=1, center=True).mean()

            fig.add_trace(go.Scatter(x=df1['Time'], y=y_ra1, mode='lines',
                                    name=f'{prot1} RMSD RA', line=dict(color='lightblue')))
            fig.add_trace(go.Scatter(x=df2['Time'], y=y_ra2, mode='lines',
                                    name=f'{prot2} RMSD RA', line=dict(color='#e84f82')))

    elif plot_type == 'rmsf':
        if 'RMSF' in df1.columns:
            fig.add_trace(go.Scatter(x=df1['Residue'], y=df1['RMSF'], mode='lines',
                                    name=prot1, line=dict(color='lightblue')))
            fig.add_trace(go.Scatter(x=df2['Residue'], y=df2['RMSF'], mode='lines',
                                    name=prot2, line=dict(color='#e84f82')))

    elif plot_type == 'sasa':
        if 'Area' in df1.columns:
            y_ra1 = df1['Area'].rolling(window=roll_avg, min_periods=1, center=True).mean()
            y_ra2 = df2['Area'].rolling(window=roll_avg, min_periods=1, center=True).mean()

            fig.add_trace(go.Scatter(x=df1['Time'], y=y_ra1, mode='lines',
                                    name=f'{prot1} SASA', line=dict(color='lightblue')))
            fig.add_trace(go.Scatter(x=df2['Time'], y=y_ra2, mode='lines',
                                    name=f'{prot2} SASA', line=dict(color='#e84f82')))

    elif plot_type == 'pca':
        return compare_pca_plots(df1, meta1, df2, meta2, prot1, prot2)


    else:
        # Generic XY comparison
        x_col = df1.columns[0]
        y_col = df1.columns[1] if len(df1.columns) > 1 else 'Y'

        if x_col in df1.columns and y_col in df1.columns:
            y_ra1 = df1[y_col].rolling(window=roll_avg, min_periods=1, center=True).mean()
            y_ra2 = df2[y_col].rolling(window=roll_avg, min_periods=1, center=True).mean()

            fig.add_trace(go.Scatter(x=df1[x_col], y=y_ra1, mode='lines',
                                    name=f'{prot1} {y_col} RA', line=dict(color='lightblue')))
            fig.add_trace(go.Scatter(x=df2[x_col], y=y_ra2, mode='lines',
                                    name=f'{prot2} {y_col} RA', line=dict(color='#e84f82')))

    fig.update_layout(
        title=f'Comparison - {meta1.title}',
        xaxis_title=meta1.x_label or 'X',
        yaxis_title=meta1.y_label or 'Y',
        template='plotly_white',
        legend_title='Protein',
        autosize=True
    )
    return fig


def compare_pca_plots(df1: pd.DataFrame, meta1: XVGMetadata,
                      df2: pd.DataFrame, meta2: XVGMetadata,
                      prot1: str = 'prot_1', prot2: str = 'prot_2') -> go.Figure:
    """Create PCA comparison plot with 2D density contour using viridis.
    Two plots side-by-side with same color scale."""
    # Calculate global min/max for consistent color scale
    all_x = list(df1['PC1']) + list(df2['PC1'])
    all_y = list(df1['PC2']) + list(df2['PC2'])

    # Use make_subplots for side-by-side layout
    fig = make_subplots(rows=1, cols=2, subplot_titles=[prot1, prot2])

    # Add PCA 1 (left)
    fig.add_trace(go.Histogram2dContour(
        x=df1['PC1'], y=df1['PC2'],
        name=prot1,
        colorscale='Viridis',
        showscale=False,
        ncontours=30
    ), row=1, col=1)

    # Add PCA 2 (right)
    fig.add_trace(go.Histogram2dContour(
        x=df2['PC1'], y=df2['PC2'],
        name=prot2,
        colorscale='Viridis',  # Same viridis scale for comparison
        showscale=False,
        ncontours=30
    ), row=1, col=2)

    fig.update_xaxes(title_text='PC1', row=1, col=1)
    fig.update_yaxes(title_text='PC2', row=1, col=1)
    fig.update_xaxes(title_text='PC1', row=1, col=2)
    fig.update_yaxes(title_text='PC2', row=1, col=2)

    fig.update_layout(
        title=f'PCA Comparison - {meta1.title}',
        height=900,
        autosize=True,
        template='plotly_white'
    )
    return fig


def compare_rmsf_plots(df1: pd.DataFrame, meta1: XVGMetadata,
                       df2: pd.DataFrame, meta2: XVGMetadata,
                       prot1: str = 'prot_1', prot2: str = 'prot_2') -> go.Figure:
    """Create RMSF comparison plot."""
    fig = go.Figure()

    if 'RMSF' in df1.columns:
        fig.add_trace(go.Scatter(x=df1['Residue'], y=df1['RMSF'], mode='lines',
                                name=prot1, line=dict(color='lightblue')))
        fig.add_trace(go.Scatter(x=df2['Residue'], y=df2['RMSF'], mode='lines',
                                name=prot2, line=dict(color='#e84f82')))

    fig.update_layout(
        title='RMSF of Ca atoms',
        xaxis_title='Residue',
        yaxis_title='RMSF (nm)',
        template='plotly_white',
        legend_title='Protein',
        autosize=True
    )
    return fig


def compare_sasa_plots(df1: pd.DataFrame, meta1: XVGMetadata,
                       df2: pd.DataFrame, meta2: XVGMetadata,
                       prot1: str = 'prot_1', prot2: str = 'prot_2',
                       roll_avg: int = 50) -> go.Figure:
    """Create SASA comparison plot."""
    fig = go.Figure()

    if 'Area' in df1.columns:
        y_ra1 = df1['Area'].rolling(window=roll_avg, min_periods=1, center=True).mean()
        y_ra2 = df2['Area'].rolling(window=roll_avg, min_periods=1, center=True).mean()

        fig.add_trace(go.Scatter(x=df1['Time'], y=y_ra1, mode='lines',
                                name=f'{prot1} SASA', line=dict(color='lightblue')))
        fig.add_trace(go.Scatter(x=df2['Time'], y=y_ra2, mode='lines',
                                name=f'{prot2} SASA', line=dict(color='#e84f82')))

    fig.update_layout(
        title='Solvent Accessible Surface Area',
        xaxis_title='Time (ns)',
        yaxis_title='Area (nm^2)',
        template='plotly_white',
        legend_title='Protein',
        autosize=True
    )
    return fig


# Pastel color palette for multi-file comparison
PASTEL_COLORS: List[str] = [
    '#AEC6CF',  # pastel blue
    '#FFB7B2',  # pastel salmon
    '#B5EAD7',  # pastel mint
    '#E2B6CF',  # pastel lavender
    '#FFDAC1',  # pastel peach
    '#C7CEEA',  # pastel periwinkle
    '#FFD8B1',  # pastel orange
    '#B4F8C8',  # pastel green
    '#A0E7E5',  # pastel teal
    '#FBE7C6',  # pastel yellow
    '#D4A5A5',  # pastel rose
    '#9ED2C6',  # pastel sage
]


def get_pastel_color(index: int) -> str:
    """Return a pastel color for the given index, cycling if needed."""
    return PASTEL_COLORS[index % len(PASTEL_COLORS)]


def compare_multi_plots(
    data_list: List[Tuple[pd.DataFrame, XVGMetadata, str]],
    roll_avg: int = 50
) -> go.Figure:
    """
    Create comparison plot for multiple XVG files of the same type.
    Dispatches to type-specific multi-compare functions.

    Args:
        data_list: List of (DataFrame, XVGMetadata, label) tuples
        roll_avg: Rolling average window size

    Returns:
        Plotly figure object
    """
    if not data_list:
        raise ValueError("No data provided for comparison")

    plot_type = data_list[0][1].plot_type

    if plot_type == 'energy':
        return compare_multi_energy(data_list, roll_avg)
    elif plot_type == 'pca':
        return compare_multi_pca(data_list)
    elif plot_type == 'rmsf':
        return compare_multi_rmsf(data_list)
    elif plot_type == 'sasa':
        return compare_multi_sasa(data_list, roll_avg)
    elif plot_type == 'rmsd':
        return compare_multi_timeseries(data_list, 'RMSD', 'RMSD (nm)',
                                         'RMSD of backbone', roll_avg)
    elif plot_type == 'gyration':
        return compare_multi_timeseries(data_list, 'Rg', 'Rg (nm)',
                                         'Radius of gyration', roll_avg)
    else:
        return compare_multi_generic(data_list, roll_avg)


def compare_multi_energy(
    data_list: List[Tuple[pd.DataFrame, XVGMetadata, str]],
    roll_avg: int = 50
) -> go.Figure:
    """Multi-file energy comparison with pastel colors."""
    energy_cols = ['Potential', 'Temperature', 'Pressure', 'Density']

    # Find columns available in ALL dataframes
    available_cols = [col for col in energy_cols
                      if all(col in d[0].columns for d in data_list)]
    additional_cols = [col for col in data_list[0][0].columns
                       if col not in ['Time'] + energy_cols
                       and all(col in d[0].columns for d in data_list)]
    all_y_cols = available_cols + additional_cols

    n_y = len(all_y_cols)
    if n_y == 0:
        return compare_multi_generic(data_list, roll_avg)

    n_grid_cols = min(2, n_y) if n_y > 1 else 1
    n_rows = (n_y + n_grid_cols - 1) // n_grid_cols

    fig = make_subplots(rows=n_rows, cols=n_grid_cols,
                        subplot_titles=all_y_cols)

    y_titles_map = {
        'Pressure': 'Pressure (bar)',
        'Temperature': 'Temperature (K)',
        'Potential': 'Potential (kJ/mol)',
        'Density': 'Density (kg/m^3)'
    }

    for idx, y_col in enumerate(all_y_cols):
        row = (idx // n_grid_cols) + 1
        col = (idx % n_grid_cols) + 1

        for file_idx, (df, meta, label) in enumerate(data_list):
            if y_col in df.columns:
                y_ra = df[y_col].rolling(window=roll_avg, min_periods=1,
                                          center=True).mean()
                show_legend = (idx == 0)
                fig.add_trace(
                    go.Scatter(x=df['Time'], y=y_ra, mode='lines',
                              name=label, legendgroup=label,
                              showlegend=show_legend,
                              line=dict(color=get_pastel_color(file_idx))),
                    row=row, col=col
                )

        y_title = y_titles_map.get(y_col, y_col)
        fig.update_yaxes(title_text=y_title, row=row, col=col)

    fig.update_xaxes(title_text='Time (ps)', row=n_rows, col=1)
    if n_grid_cols > 1 and n_rows > 1:
        fig.update_xaxes(title_text='Time (ps)', row=n_rows, col=n_grid_cols)

    fig.update_layout(
        title_text=f'Energy Components Comparison ({len(data_list)} files)',
        height=400 + (n_rows - 1) * 500,
        autosize=True
    )
    return fig


def compare_multi_pca(
    data_list: List[Tuple[pd.DataFrame, XVGMetadata, str]]
) -> go.Figure:
    """Multi-file PCA comparison with side-by-side contour plots."""
    n_files = len(data_list)
    n_cols = min(3, n_files)
    n_rows = (n_files + n_cols - 1) // n_cols

    titles = [d[2] for d in data_list]
    fig = make_subplots(rows=n_rows, cols=n_cols, subplot_titles=titles)

    for idx, (df, meta, label) in enumerate(data_list):
        row = (idx // n_cols) + 1
        col = (idx % n_cols) + 1

        if 'PC1' in df.columns and 'PC2' in df.columns:
            fig.add_trace(go.Histogram2dContour(
                x=df['PC1'], y=df['PC2'],
                name=label,
                colorscale='Viridis',
                showscale=False,
                ncontours=30
            ), row=row, col=col)

        fig.update_xaxes(title_text='PC1', row=row, col=col)
        fig.update_yaxes(title_text='PC2', row=row, col=col)

    fig.update_layout(
        title=f'PCA Comparison ({n_files} files)',
        height=max(900, 400 + (n_rows - 1) * 500),
        autosize=True,
        template='plotly_white'
    )
    return fig


def compare_multi_rmsf(
    data_list: List[Tuple[pd.DataFrame, XVGMetadata, str]]
) -> go.Figure:
    """Multi-file RMSF comparison with pastel colors."""
    fig = go.Figure()

    for idx, (df, meta, label) in enumerate(data_list):
        if 'RMSF' in df.columns and 'Residue' in df.columns:
            fig.add_trace(go.Scatter(
                x=df['Residue'], y=df['RMSF'], mode='lines',
                name=label, line=dict(color=get_pastel_color(idx))
            ))

    fig.update_layout(
        title=f'RMSF of Ca atoms ({len(data_list)} files)',
        xaxis_title='Residue',
        yaxis_title='RMSF (nm)',
        template='plotly_white',
        legend_title='Source',
        autosize=True
    )
    return fig


def compare_multi_sasa(
    data_list: List[Tuple[pd.DataFrame, XVGMetadata, str]],
    roll_avg: int = 50
) -> go.Figure:
    """Multi-file SASA comparison with pastel colors."""
    fig = go.Figure()

    for idx, (df, meta, label) in enumerate(data_list):
        if 'Area' in df.columns and 'Time' in df.columns:
            y_ra = df['Area'].rolling(window=roll_avg, min_periods=1,
                                       center=True).mean()
            fig.add_trace(go.Scatter(
                x=df['Time'], y=y_ra, mode='lines',
                name=f'{label} SASA',
                line=dict(color=get_pastel_color(idx))
            ))

    fig.update_layout(
        title=f'Solvent Accessible Surface Area ({len(data_list)} files)',
        xaxis_title='Time (ns)',
        yaxis_title='Area (nm^2)',
        template='plotly_white',
        legend_title='Source',
        autosize=True
    )
    return fig


def compare_multi_timeseries(
    data_list: List[Tuple[pd.DataFrame, XVGMetadata, str]],
    y_col: str,
    y_label: str,
    title: str,
    roll_avg: int = 50
) -> go.Figure:
    """Multi-file time-series comparison (RMSD, gyration, etc.) with pastel colors."""
    fig = go.Figure()

    for idx, (df, meta, label) in enumerate(data_list):
        if y_col in df.columns and 'Time' in df.columns:
            y_ra = df[y_col].rolling(window=roll_avg, min_periods=1,
                                      center=True).mean()
            fig.add_trace(go.Scatter(
                x=df['Time'], y=y_ra, mode='lines',
                name=label, line=dict(color=get_pastel_color(idx))
            ))

    fig.update_layout(
        title=f'{title} ({len(data_list)} files)',
        xaxis_title='Time (ns)',
        yaxis_title=y_label,
        template='plotly_white',
        legend_title='Source',
        autosize=True
    )
    return fig


def compare_multi_generic(
    data_list: List[Tuple[pd.DataFrame, XVGMetadata, str]],
    roll_avg: int = 50
) -> go.Figure:
    """Multi-file generic/xy comparison with pastel colors."""
    fig = go.Figure()

    for idx, (df, meta, label) in enumerate(data_list):
        x_col = df.columns[0]
        y_col = df.columns[1] if len(df.columns) > 1 else df.columns[0]

        if x_col in df.columns and y_col in df.columns:
            y_ra = df[y_col].rolling(window=roll_avg, min_periods=1,
                                      center=True).mean()
            fig.add_trace(go.Scatter(
                x=df[x_col], y=y_ra, mode='lines',
                name=label, line=dict(color=get_pastel_color(idx))
            ))

    meta0 = data_list[0][1]
    fig.update_layout(
        title=f'Comparison ({len(data_list)} files) - {meta0.title}',
        xaxis_title=meta0.x_label or data_list[0][0].columns[0],
        yaxis_title=meta0.y_label or 'Y',
        template='plotly_white',
        legend_title='Source',
        autosize=True
    )
    return fig
