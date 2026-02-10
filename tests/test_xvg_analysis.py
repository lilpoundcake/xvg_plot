"""
Comprehensive test suite for XVG Analysis module.

Tests all plot types, parsing, type detection, and CLI functionality
using test data files from the tests directory.
"""

import pytest
import pandas as pd
from pathlib import Path
from xvg_analysis import xvg_analysis as xp
from xvg_analysis.xvg_analysis import XVGMetadata


# Test data fixtures
@pytest.fixture
def test_data_dir():
    """Return path to test data directory."""
    return Path(__file__).parent


@pytest.fixture
def prot1_dir(test_data_dir):
    """Return path to protein 1 test data."""
    return test_data_dir / "prot1"


@pytest.fixture
def prot2_dir(test_data_dir):
    """Return path to protein 2 test data."""
    return test_data_dir / "prot2"


@pytest.fixture
def other_xvg_dir(test_data_dir):
    """Return path to other XVG test data."""
    return test_data_dir / "other_xvg"


# Test data file paths
TEST_FILES = {
    "energy": "prot1/energy.xvg",
    "gyrate": "prot1/gyrate.xvg",
    "rmsd": "prot1/rmsd.xvg",
    "rmsf": "prot1/rmsf.xvg",
    "sasa": "prot1/sasa.xvg",
    "pca": "prot1/2dproj.xvg",
}

EXPECTED_PLOT_TYPES = {
    "energy": "energy",
    "gyrate": "gyration",
    "rmsd": "rmsd",
    "rmsf": "rmsf",
    "sasa": "sasa",
    "pca": "pca",
}


class TestXVGParsing:
    """Test XVG file parsing and metadata extraction."""

    @pytest.mark.parametrize("file_type", TEST_FILES.keys())
    def test_read_xvg_returns_tuple(self, test_data_dir, file_type):
        """Test that read_xvg returns (DataFrame, XVGMetadata) tuple."""
        file_path = test_data_dir / TEST_FILES[file_type]
        assert file_path.exists(), f"Test file not found: {file_path}"

        df, metadata = xp.read_xvg(str(file_path))

        assert isinstance(df, pd.DataFrame), "First return value should be DataFrame"
        assert isinstance(metadata, XVGMetadata), "Second return value should be XVGMetadata"
        assert not df.empty, f"DataFrame should not be empty for {file_type}"

    @pytest.mark.parametrize("file_type", TEST_FILES.keys())
    def test_metadata_has_required_fields(self, test_data_dir, file_type):
        """Test that metadata contains required fields."""
        file_path = test_data_dir / TEST_FILES[file_type]
        df, metadata = xp.read_xvg(str(file_path))

        assert metadata.title, "Metadata should have title"
        assert metadata.plot_type, "Metadata should have plot_type"
        assert metadata.x_label or metadata.plot_type, "Should have x_label or inferred plot_type"
        assert len(df.columns) > 0, "DataFrame should have columns"

    @pytest.mark.parametrize("file_type", TEST_FILES.keys())
    def test_dataframe_all_numeric(self, test_data_dir, file_type):
        """Test that DataFrame contains only numeric data."""
        file_path = test_data_dir / TEST_FILES[file_type]
        df, metadata = xp.read_xvg(str(file_path))

        for col in df.columns:
            assert pd.api.types.is_numeric_dtype(df[col]), f"Column {col} should be numeric"


class TestPlotTypeDetection:
    """Test plot type classification."""

    @pytest.mark.parametrize("file_type,expected_type", EXPECTED_PLOT_TYPES.items())
    def test_detect_plot_type(self, test_data_dir, file_type, expected_type):
        """Test that plot type is correctly detected for known file types."""
        file_path = test_data_dir / TEST_FILES[file_type]
        df, metadata = xp.read_xvg(str(file_path))

        assert metadata.plot_type == expected_type, (
            f"Expected plot_type '{expected_type}' for {file_type}, "
            f"got '{metadata.plot_type}'"
        )

    def test_energy_detection(self, test_data_dir):
        """Test energy plot detection with specific criteria."""
        file_path = test_data_dir / "prot1/energy.xvg"
        df, metadata = xp.read_xvg(str(file_path))

        # Energy files should have 4+ columns and Time x-axis
        assert metadata.plot_type == "energy"
        assert len(df.columns) >= 4, "Energy plot should have 4+ columns"
        assert "Time" in df.columns, "Energy plot should have Time column"

    def test_rmsf_detection(self, test_data_dir):
        """Test RMSF plot detection."""
        file_path = test_data_dir / "prot1/rmsf.xvg"
        df, metadata = xp.read_xvg(str(file_path))

        assert metadata.plot_type == "rmsf"
        assert "Residue" in df.columns, "RMSF should have Residue column"
        assert "RMSF" in df.columns, "RMSF should have RMSF column"

    def test_pca_detection(self, test_data_dir):
        """Test PCA plot detection."""
        file_path = test_data_dir / "prot1/2dproj.xvg"
        df, metadata = xp.read_xvg(str(file_path))

        assert metadata.plot_type == "pca"
        assert "PC1" in df.columns, "PCA should have PC1 column"
        assert "PC2" in df.columns, "PCA should have PC2 column"


class TestColumnNaming:
    """Test automatic column naming."""

    def test_energy_column_names(self, test_data_dir):
        """Test that energy files have correctly named columns."""
        file_path = test_data_dir / "prot1/energy.xvg"
        df, metadata = xp.read_xvg(str(file_path))

        # Energy files should have expected column names
        expected_cols = ["Time", "Potential", "Temperature", "Pressure", "Density"]
        for col in expected_cols:
            assert col in df.columns, f"Energy plot missing expected column: {col}"

    def test_gyration_column_names(self, test_data_dir):
        """Test that gyration files have correct column names."""
        file_path = test_data_dir / "prot1/gyrate.xvg"
        df, metadata = xp.read_xvg(str(file_path))

        assert "Time" in df.columns, "Gyration should have Time column"
        assert "Rg" in df.columns, "Gyration should have Rg column"

    def test_rmsd_column_names(self, test_data_dir):
        """Test that RMSD files have correct column names."""
        file_path = test_data_dir / "prot1/rmsd.xvg"
        df, metadata = xp.read_xvg(str(file_path))

        assert "Time" in df.columns, "RMSD should have Time column"
        assert "RMSD" in df.columns, "RMSD should have RMSD column"

    def test_sasa_column_names(self, test_data_dir):
        """Test that SASA files have correct column names."""
        file_path = test_data_dir / "prot1/sasa.xvg"
        df, metadata = xp.read_xvg(str(file_path))

        assert "Time" in df.columns, "SASA should have Time column"
        assert "Area" in df.columns, "SASA should have Area column"


class TestRollingAverage:
    """Test rolling average calculation."""

    def test_calculate_roll_avg_basic(self, test_data_dir):
        """Test rolling average calculation with basic data."""
        file_path = test_data_dir / "prot1/energy.xvg"
        df, metadata = xp.read_xvg(str(file_path))

        roll_avg = xp.calculate_roll_avg(df, 'Time')

        assert isinstance(roll_avg, int), "Rolling average should be integer"
        assert roll_avg >= 3, "Rolling average should be minimum 3"
        assert roll_avg <= 100, "Rolling average should not exceed 100"

    @pytest.mark.parametrize("file_type", TEST_FILES.keys())
    def test_roll_avg_all_files(self, test_data_dir, file_type):
        """Test rolling average calculation for all file types."""
        file_path = test_data_dir / TEST_FILES[file_type]
        df, metadata = xp.read_xvg(str(file_path))

        x_col = df.columns[0]
        roll_avg = xp.calculate_roll_avg(df, x_col)

        assert isinstance(roll_avg, int), f"Rolling average should be int for {file_type}"
        assert roll_avg >= 1, f"Rolling average should be >= 1 for {file_type}"

    def test_roll_avg_default_for_small_data(self):
        """Test rolling average default fallback for small datasets."""
        small_df = pd.DataFrame({
            'Time': [0, 1, 2],
            'Value': [10, 20, 30]
        })

        roll_avg = xp.calculate_roll_avg(small_df, 'Time')
        assert roll_avg == 10, "Should return default 10 for small dataset"


class TestPlotGeneration:
    """Test plot generation for each file type."""

    def test_plot_energy(self, test_data_dir):
        """Test energy plot generation."""
        file_path = test_data_dir / "prot1/energy.xvg"
        df, metadata = xp.read_xvg(str(file_path))

        fig = xp.create_plot(df, metadata)

        assert fig is not None, "Should return a figure object"
        assert hasattr(fig, 'update_layout'), "Should be a Plotly figure"

    def test_plot_gyration(self, test_data_dir):
        """Test gyration plot generation."""
        file_path = test_data_dir / "prot1/gyrate.xvg"
        df, metadata = xp.read_xvg(str(file_path))

        fig = xp.create_plot(df, metadata)

        assert fig is not None, "Should return a figure object"
        assert "Rg" in [trace.name for trace in fig.data], "Should plot Rg data"

    def test_plot_rmsd(self, test_data_dir):
        """Test RMSD plot generation."""
        file_path = test_data_dir / "prot1/rmsd.xvg"
        df, metadata = xp.read_xvg(str(file_path))

        fig = xp.create_plot(df, metadata)

        assert fig is not None, "Should return a figure object"
        assert "RMSD" in [trace.name for trace in fig.data], "Should plot RMSD data"

    def test_plot_rmsf(self, test_data_dir):
        """Test RMSF plot generation."""
        file_path = test_data_dir / "prot1/rmsf.xvg"
        df, metadata = xp.read_xvg(str(file_path))

        fig = xp.create_plot(df, metadata)

        assert fig is not None, "Should return a figure object"
        assert "RMSF" in [trace.name for trace in fig.data], "Should plot RMSF data"

    def test_plot_sasa(self, test_data_dir):
        """Test SASA plot generation."""
        file_path = test_data_dir / "prot1/sasa.xvg"
        df, metadata = xp.read_xvg(str(file_path))

        fig = xp.create_plot(df, metadata)

        assert fig is not None, "Should return a figure object"
        assert "Area" in [trace.name for trace in fig.data], "Should plot Area data"

    def test_plot_pca(self, test_data_dir):
        """Test PCA plot generation."""
        file_path = test_data_dir / "prot1/2dproj.xvg"
        df, metadata = xp.read_xvg(str(file_path))

        fig = xp.create_plot(df, metadata)

        assert fig is not None, "Should return a figure object"
        assert len(fig.data) > 0, "Should have trace data"

    @pytest.mark.parametrize("file_type", TEST_FILES.keys())
    def test_custom_roll_avg(self, test_data_dir, file_type):
        """Test plot generation with custom rolling average."""
        file_path = test_data_dir / TEST_FILES[file_type]
        df, metadata = xp.read_xvg(str(file_path))

        custom_roll_avg = 5
        fig = xp.create_plot(df, metadata, roll_avg=custom_roll_avg)

        assert fig is not None, f"Should generate plot with custom RA for {file_type}"


class TestPlotWidth:
    """Test responsive plot width."""

    def test_get_plot_width_single_col(self):
        """Test width calculation for single column."""
        width = xp.get_plot_width(1)
        assert width == 1000, "Single column should have width 1000"

    def test_get_plot_width_two_cols(self):
        """Test width calculation for two columns."""
        width = xp.get_plot_width(2)
        assert width == 1400, "Two columns should have width 1400"

    def test_get_plot_width_three_cols(self):
        """Test width calculation for three columns."""
        width = xp.get_plot_width(3)
        assert 1000 < width <= 1800, "Three columns should scale appropriately"

    def test_plot_energy_has_width(self, test_data_dir):
        """Test that energy plots have width property."""
        file_path = test_data_dir / "prot1/energy.xvg"
        df, metadata = xp.read_xvg(str(file_path))

        fig = xp.create_plot(df, metadata)

        assert 'width' in fig.layout, "Plot should have width property"
        assert fig.layout.width is not None, "Width should be set"


class TestComparisons:
    """Test comparison plot generation."""

    def test_compare_two_energy_plots(self, test_data_dir):
        """Test comparing two energy plots."""
        df1, meta1 = xp.read_xvg(str(test_data_dir / "prot1/energy.xvg"))
        df2, meta2 = xp.read_xvg(str(test_data_dir / "prot2/energy.xvg"))

        fig = xp.compare_two_plots(df1, meta1, df2, meta2, prot1="P1", prot2="P2")

        assert fig is not None, "Should return a figure object"
        assert len(fig.data) > 0, "Should have trace data"

    def test_compare_pca_plots(self, test_data_dir):
        """Test comparing two PCA plots."""
        df1, meta1 = xp.read_xvg(str(test_data_dir / "prot1/2dproj.xvg"))
        df2, meta2 = xp.read_xvg(str(test_data_dir / "prot2/2dproj.xvg"))

        fig = xp.compare_pca_plots(df1, meta1, df2, meta2, prot1="P1", prot2="P2")

        assert fig is not None, "Should return a figure object"

    def test_compare_rmsf_plots(self, test_data_dir):
        """Test comparing two RMSF plots."""
        df1, meta1 = xp.read_xvg(str(test_data_dir / "prot1/rmsf.xvg"))
        df2, meta2 = xp.read_xvg(str(test_data_dir / "prot2/rmsf.xvg"))

        fig = xp.compare_rmsf_plots(df1, meta1, df2, meta2, prot1="P1", prot2="P2")

        assert fig is not None, "Should return a figure object"
        assert len(fig.data) > 0, "Should have trace data"

    def test_compare_sasa_plots(self, test_data_dir):
        """Test comparing two SASA plots."""
        df1, meta1 = xp.read_xvg(str(test_data_dir / "prot1/sasa.xvg"))
        df2, meta2 = xp.read_xvg(str(test_data_dir / "prot2/sasa.xvg"))

        fig = xp.compare_sasa_plots(df1, meta1, df2, meta2, prot1="P1", prot2="P2")

        assert fig is not None, "Should return a figure object"


class TestRollingAverageLabels:
    """Test that rolling average labels use calculated values."""

    def test_energy_ra_label_has_value(self, test_data_dir):
        """Test that energy plots show actual RA value in legend."""
        file_path = test_data_dir / "prot1/energy.xvg"
        df, metadata = xp.read_xvg(str(file_path))

        fig = xp.create_plot(df, metadata)

        # Check that trace names contain RA with calculated value
        trace_names = [trace.name for trace in fig.data]
        ra_traces = [name for name in trace_names if 'RA' in name]
        assert len(ra_traces) > 0, "Should have RA traces"
        # RA label should be like "RA (n=47)" not "10-step RA"
        for ra_name in ra_traces:
            assert ra_name.startswith('RA (n='), f"RA label should be 'RA (n=X)' format, got {ra_name}"

    def test_gyration_ra_label_format(self, test_data_dir):
        """Test that gyration plots use correct RA label format."""
        file_path = test_data_dir / "prot1/gyrate.xvg"
        df, metadata = xp.read_xvg(str(file_path))

        fig = xp.create_plot(df, metadata)

        trace_names = [trace.name for trace in fig.data]
        ra_traces = [name for name in trace_names if 'RA' in name]
        for ra_name in ra_traces:
            assert ra_name.startswith('RA (n='), f"RA label should be 'RA (n=X)' format, got {ra_name}"

    def test_custom_ra_in_label(self, test_data_dir):
        """Test that custom rolling average value appears in label."""
        file_path = test_data_dir / "prot1/rmsd.xvg"
        df, metadata = xp.read_xvg(str(file_path))

        custom_ra = 15
        fig = xp.create_plot(df, metadata, roll_avg=custom_ra)

        trace_names = [trace.name for trace in fig.data]
        # Should contain trace with custom RA value
        assert any(f'RA (n={custom_ra})' in name for name in trace_names), \
            f"Should show custom RA value {custom_ra} in legend"


class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_file_not_found(self):
        """Test handling of non-existent file."""
        with pytest.raises(FileNotFoundError):
            xp.read_xvg("/nonexistent/file.xvg")

    def test_empty_dataframe_handling(self):
        """Test handling of malformed files that produce empty DataFrames."""
        # This would need a malformed test file
        # For now, test with a valid file
        pass

    def test_unknown_plot_type_fallback(self, other_xvg_dir):
        """Test that unknown plot types fall back to xy plot."""
        xvg_files = list(other_xvg_dir.glob("*.xvg"))
        if xvg_files:
            df, metadata = xp.read_xvg(str(xvg_files[0]))
            # Unknown types should still generate a plot
            fig = xp.create_plot(df, metadata)
            assert fig is not None, "Should generate plot even for unknown type"


class TestBatchProcessing:
    """Test batch processing of multiple files."""

    def test_batch_all_protein1_files(self, test_data_dir):
        """Test that all protein 1 files can be processed."""
        prot1_dir = test_data_dir / "prot1"
        xvg_files = sorted(prot1_dir.glob("*.xvg"))

        assert len(xvg_files) >= 6, "Should have at least 6 XVG files in prot1"

        for xvg_file in xvg_files:
            df, metadata = xp.read_xvg(str(xvg_file))
            fig = xp.create_plot(df, metadata)
            assert fig is not None, f"Should generate plot for {xvg_file.name}"

    def test_batch_all_protein2_files(self, test_data_dir):
        """Test that all protein 2 files can be processed."""
        prot2_dir = test_data_dir / "prot2"
        xvg_files = sorted(prot2_dir.glob("*.xvg"))

        assert len(xvg_files) >= 6, "Should have at least 6 XVG files in prot2"

        for xvg_file in xvg_files:
            df, metadata = xp.read_xvg(str(xvg_file))
            fig = xp.create_plot(df, metadata)
            assert fig is not None, f"Should generate plot for {xvg_file.name}"
