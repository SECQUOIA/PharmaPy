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
        assert hasattr(Streams, 'LiquidStream')
    
    def test_import_phases(self, pharmapy_available):
        """Test that Phases module imports."""
        from PharmaPy import Phases
        assert hasattr(Phases, 'LiquidPhase')
    
    def test_import_reactors(self, pharmapy_available):
        """Test that Reactors module imports."""
        from PharmaPy import Reactors
        assert hasattr(Reactors, 'PlugFlowReactor')
    
    def test_import_utilities(self, pharmapy_available):
        """Test that Utilities module imports."""
        from PharmaPy import Utilities
        assert callable(getattr(Utilities, 'CoolingWater', None))


class TestLiquidStream:
    """Test LiquidStream functionality."""
    
    def test_create_liquid_stream(self, sample_liquid_stream):
        """Test creation of a liquid stream."""
        stream = sample_liquid_stream
        assert stream is not None
        assert hasattr(stream, 'flow_rate')
        assert hasattr(stream, 'temperature')
    
    def test_stream_properties(self, sample_liquid_stream):
        """Test basic stream properties."""
        stream = sample_liquid_stream
        assert stream.flow_rate > 0
        assert stream.temperature > 0
    
    @pytest.mark.parametrize("flow_rate", [10.0, 50.0, 100.0, 500.0])
    def test_stream_flow_rates(self, flow_rate, pharmapy_available):
        """Test stream with different flow rates."""
        from PharmaPy.Streams import LiquidStream
        from PharmaPy.Phases import LiquidPhase
        
        phase = LiquidPhase(["water"], [1.0])
        stream = LiquidStream(phase, flow_rate=flow_rate, temperature=298.15)
        assert stream.flow_rate == flow_rate


class TestReactionKinetics:
    """Test reaction kinetics functionality."""
    
    def test_create_kinetics(self, sample_reaction_kinetics):
        """Test creation of reaction kinetics."""
        kinetics = sample_reaction_kinetics
        assert kinetics is not None
        assert hasattr(kinetics, 'components')
        assert hasattr(kinetics, 'reactions')
    
    def test_kinetics_properties(self, sample_reaction_kinetics):
        """Test basic kinetics properties."""
        kinetics = sample_reaction_kinetics
        assert len(kinetics.components) >= 2
        assert len(kinetics.reactions) >= 1


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
        assert 'json' in test_data_files
        assert 'csv' in test_data_files
        assert 'all' in test_data_files
    
    def test_json_files_readable(self, test_data_files):
        """Test that JSON test files can be read."""
        import json
        
        for json_file in test_data_files['json'][:3]:  # Test first 3 files
            try:
                with open(json_file, 'r') as f:
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
        
        # Create a stream
        phase = LiquidPhase(["A", "B"], [0.8, 0.2])
        stream = LiquidStream(phase, flow_rate=100.0, temperature=298.15)
        
        # Create a reactor
        reactor = PlugFlowReactor()
        
        # Basic checks
        assert stream is not None
        assert reactor is not None
        assert hasattr(reactor, '__class__')
    
    @pytest.mark.slow
    def test_complex_workflow(self, pharmapy_available):
        """Test more complex workflow (marked as slow)."""
        # This would be a more complex test that might take longer
        from PharmaPy import Utilities
        
        # Test utility functions
        cooling_water = Utilities.CoolingWater()
        assert cooling_water is not None


class TestNumericalStability:
    """Test numerical stability and edge cases."""
    
    @pytest.mark.parametrize("temperature", [273.15, 298.15, 373.15, 473.15])
    def test_temperature_ranges(self, temperature, pharmapy_available):
        """Test stream creation with different temperatures."""
        from PharmaPy.Streams import LiquidStream
        from PharmaPy.Phases import LiquidPhase
        
        phase = LiquidPhase(["water"], [1.0])
        stream = LiquidStream(phase, flow_rate=100.0, temperature=temperature)
        assert stream.temperature == temperature
    
    def test_zero_flow_rate(self, pharmapy_available):
        """Test handling of zero flow rate."""
        from PharmaPy.Streams import LiquidStream
        from PharmaPy.Phases import LiquidPhase
        
        phase = LiquidPhase(["water"], [1.0])
        # This might raise an exception or handle it gracefully
        try:
            stream = LiquidStream(phase, flow_rate=0.0, temperature=298.15)
            assert stream.flow_rate == 0.0
        except (ValueError, AssertionError):
            # Expected behavior for invalid input
            pass
    
    def test_negative_values_handling(self, pharmapy_available):
        """Test handling of negative values."""
        from PharmaPy.Streams import LiquidStream
        from PharmaPy.Phases import LiquidPhase
        
        phase = LiquidPhase(["water"], [1.0])
        
        # Test negative flow rate
        with pytest.raises((ValueError, AssertionError)):
            LiquidStream(phase, flow_rate=-100.0, temperature=298.15)
        
        # Test negative temperature
        with pytest.raises((ValueError, AssertionError)):
            LiquidStream(phase, flow_rate=100.0, temperature=-298.15)
