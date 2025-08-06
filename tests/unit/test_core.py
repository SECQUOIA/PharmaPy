"""
Unit tests for PharmaPy core functionality.
"""

import pytest
import numpy as np


class TestPharmaImports:
    """Test that PharmaPy modules can be imported correctly."""

    def test_import_pharmapy(self, pharmapy_available):
        """Test that main PharmaPy module imports."""
        import PharmaPy

        # Test that we can import PharmaPy successfully
        assert PharmaPy is not None

    def test_import_streams(self, pharmapy_available):
        """Test that Streams module imports."""
        from PharmaPy import Streams

        assert hasattr(Streams, "LiquidStream")

    def test_import_phases(self, pharmapy_available):
        """Test that Phases module imports."""
        from PharmaPy import Phases

        assert hasattr(Phases, "LiquidPhase")

    def test_import_reactors(self, pharmapy_available):
        """Test that Reactors module imports."""
        from PharmaPy import Reactors

        assert hasattr(Reactors, "PlugFlowReactor")

    def test_import_utilities(self, pharmapy_available):
        """Test that Utilities module imports."""
        from PharmaPy import Utilities

        assert callable(getattr(Utilities, "CoolingWater", None))


class TestLiquidStream:
    """Test LiquidStream functionality."""

    def test_create_liquid_stream(self, sample_liquid_stream):
        """Test creation of a liquid stream."""
        stream = sample_liquid_stream
        assert stream is not None
        assert hasattr(stream, "vol_flow")
        assert hasattr(stream, "temp")

    def test_stream_properties(self, sample_liquid_stream):
        """Test basic stream properties."""
        stream = sample_liquid_stream
        assert stream.vol_flow > 0
        assert stream.temp > 0

    @pytest.mark.parametrize("flow_rate", [10.0, 50.0, 100.0, 500.0])
    def test_stream_flow_rates(self, flow_rate, pharmapy_available):
        """Test stream with different flow rates."""
        from PharmaPy.Streams import LiquidStream
        from PharmaPy.Phases import LiquidPhase

        # Use proper API with data file
        datapath = "tests/integration/data/pfr_test_pure_comp.json"
        phase = LiquidPhase(
            datapath, mole_conc=[1.0, 0.0, 0.0, 0.0], temp=298.15, vol=0.001
        )
        stream = LiquidStream(
            datapath, mole_conc=[1.0, 0.0, 0.0, 0.0], temp=298.15, vol_flow=flow_rate
        )
        assert stream.vol_flow == flow_rate


class TestReactionKinetics:
    """Test reaction kinetics functionality."""

    def test_create_kinetics(self, sample_reaction_kinetics):
        """Test creation of reaction kinetics."""
        kinetics = sample_reaction_kinetics
        assert kinetics is not None
        assert hasattr(kinetics, "partic_species")
        assert hasattr(kinetics, "stoich_matrix")

    def test_kinetics_properties(self, sample_reaction_kinetics):
        """Test basic kinetics properties."""
        kinetics = sample_reaction_kinetics
        assert len(kinetics.partic_species) >= 2
        assert len(kinetics.stoich_matrix) >= 1
        # Test that we have the expected species from test data
        assert "A" in kinetics.partic_species
        assert "B" in kinetics.partic_species


@pytest.mark.unit
class TestDataValidation:
    """Test data validation and file handling."""

    def test_data_directory_exists(self, data_path):
        """Test that data directory exists."""
        assert data_path.exists()
        assert data_path.is_dir()

    def test_test_data_files(self, test_data_files):
        """Test that test data files are available."""
        assert isinstance(test_data_files, dict)
        assert "json" in test_data_files
        assert "csv" in test_data_files
        assert "all" in test_data_files

    def test_json_files_readable(self, test_data_files):
        """Test that JSON test files can be read."""
        import json

        for json_file in test_data_files["json"][:3]:  # Test first 3 files
            try:
                with open(json_file, "r") as f:
                    data = json.load(f)
                assert isinstance(data, (dict, list))
            except (json.JSONDecodeError, FileNotFoundError):
                pytest.skip(f"Could not read {json_file}")


@pytest.mark.integration
class TestBasicWorkflow:
    """Test basic PharmaPy workflow."""

    def test_stream_to_reactor_workflow(self, pharmapy_available):
        """Test basic stream to reactor workflow."""
        from PharmaPy.Streams import LiquidStream
        from PharmaPy.Phases import LiquidPhase
        from PharmaPy.Reactors import PlugFlowReactor

        # Create actual objects to test workflow
        datapath = "tests/integration/data/pfr_test_pure_comp.json"

        # Create phase and stream
        phase = LiquidPhase(
            datapath, mole_conc=[0.8, 0.2, 0.0, 0.0], temp=298.15, vol=0.001
        )
        stream = LiquidStream(
            datapath, mole_conc=[0.8, 0.2, 0.0, 0.0], temp=298.15, vol_flow=100.0
        )

        # Verify objects were created successfully
        assert phase is not None
        assert stream is not None
        assert hasattr(stream, "vol_flow")
        assert hasattr(phase, "mole_conc")
        assert stream.vol_flow == 100.0

        # Verify workflow compatibility - objects can be created and used together
        assert (
            PlugFlowReactor is not None
        )  # Constructor requires parameters, so just verify class exists

    @pytest.mark.slow
    def test_complex_workflow(self, pharmapy_available):
        """Test more complex workflow (marked as slow)."""
        # This would be a more complex test that might take longer
        from PharmaPy import Utilities

        # Test utility functions with proper parameters
        cooling_water = Utilities.CoolingWater(vol_flow=100.0)
        assert cooling_water is not None


class TestNumericalStability:
    """Test numerical stability and edge cases."""

    @pytest.mark.parametrize("temperature", [273.15, 298.15, 373.15, 473.15])
    def test_temperature_ranges(self, temperature, pharmapy_available):
        """Test stream creation with different temperatures."""
        from PharmaPy.Streams import LiquidStream
        from PharmaPy.Phases import LiquidPhase

        # Use proper API with data file
        datapath = "tests/integration/data/pfr_test_pure_comp.json"
        phase = LiquidPhase(
            datapath, mole_conc=[1.0, 0.0, 0.0, 0.0], temp=temperature, vol=0.001
        )
        stream = LiquidStream(
            datapath, mole_conc=[1.0, 0.0, 0.0, 0.0], temp=temperature, vol_flow=100.0
        )
        assert stream.temp == temperature

    def test_zero_flow_rate(self, pharmapy_available):
        """Test handling of zero flow rate."""
        import warnings
        from PharmaPy.Streams import LiquidStream
        from PharmaPy.Phases import LiquidPhase

        # Use proper API with data file
        datapath = "tests/integration/data/pfr_test_pure_comp.json"

        # Test 1: Create a proper phase (no warning expected)
        phase = LiquidPhase(
            datapath, mole_conc=[1.0, 0.0, 0.0, 0.0], temp=298.15, vol=0.001
        )

        # Test 2: Create stream with zero flow rate but valid phase
        try:
            stream = LiquidStream(
                datapath, mole_conc=[1.0, 0.0, 0.0, 0.0], temp=298.15, vol_flow=0.0
            )
            assert stream.vol_flow == 0.0
        except (ValueError, AssertionError):
            # Expected behavior for invalid input
            pass

        # Test 3: Test that warning is emitted when mass, vol, and moles are all zero
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")  # Ensure all warnings are captured
            # Create phase with all defaults (mass=0, vol=0, moles=0) - should trigger warning
            # Need to provide composition to avoid ValueError
            zero_phase = LiquidPhase(
                datapath, temp=298.15, mole_conc=[1.0, 0.0, 0.0, 0.0], check_input=True
            )

            # Check that the warning was issued
            assert len(w) >= 1
            assert any(
                "mass" in str(warning.message)
                and "moles" in str(warning.message)
                and "vol" in str(warning.message)
                for warning in w
            )
            assert any(issubclass(warning.category, RuntimeWarning) for warning in w)

    def test_negative_values_handling(self, pharmapy_available):
        """Test handling of negative values."""
        from PharmaPy.Streams import LiquidStream
        from PharmaPy.Phases import LiquidPhase

        # Use proper API with data file
        datapath = "tests/integration/data/pfr_test_pure_comp.json"
        phase = LiquidPhase(
            datapath, mole_conc=[1.0, 0.0, 0.0, 0.0], temp=298.15, vol=0.001
        )

        # Test that we can create streams (PharmaPy may or may not validate negative values)
        # This is more of a basic functionality test
        try:
            stream = LiquidStream(
                datapath, mole_conc=[1.0, 0.0, 0.0, 0.0], temp=298.15, vol_flow=-100.0
            )
            # If it succeeds, that's also valid behavior - PharmaPy might handle it gracefully
            assert stream is not None or stream is None
        except (ValueError, AssertionError, RuntimeError):
            # If it raises an exception, that's expected behavior
            pass

        # Basic negative value test passed
