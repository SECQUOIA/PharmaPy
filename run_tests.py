#!/usr/bin/env python3
"""
PharmaPy Test Runner

This script runs all tests in the PharmaPy test suite and provides
comprehensive reporting of test results.
"""

import os
import sys
import subprocess
import argparse
import time
from pathlib import Path

def run_command(cmd, cwd=None, timeout=300):
    """Run a command and return success status and output."""
    try:
        result = subprocess.run(
            cmd,
            cwd=cwd,
            capture_output=True,
            text=True,
            timeout=timeout,
            shell=True if sys.platform.startswith('win') else False
        )
        return result.returncode == 0, result.stdout, result.stderr
    except subprocess.TimeoutExpired:
        return False, "", f"Command timed out after {timeout} seconds"
    except Exception as e:
        return False, "", str(e)

def test_import():
    """Test basic PharmaPy imports."""
    print(" Testing PharmaPy imports...")
    
    import_tests = [
        ("import PharmaPy", True),
        ("from PharmaPy import Utilities", True),
        ("from PharmaPy import Reactors", False),  # Requires assimulo
        ("from PharmaPy import Streams", False),   # Requires assimulo  
        ("from PharmaPy import Phases", False),    # Requires assimulo
        ("from PharmaPy import Kinetics", False),  # Requires assimulo
        ("from PharmaPy.Reactors import PlugFlowReactor, BatchReactor", False),
        ("from PharmaPy.Streams import LiquidStream", False),
        ("from PharmaPy.Phases import LiquidPhase", False)
    ]
    
    failed_imports = []
    optional_failed = 0
    
    for test, required in import_tests:
        try:
            # Parse the import statement and use importlib
            if test.startswith("import "):
                # e.g., "import PharmaPy"
                module_name = test.split("import ")[1].strip()
                importlib.import_module(module_name)
            elif test.startswith("from "):
                # e.g., "from PharmaPy import Utilities"
                parts = test.split()
                module = parts[1]
                imported = parts[3]
                # Handle multiple imports separated by commas
                imported_names = [name.strip() for name in imported.split(",")]
                mod = importlib.import_module(module)
                for name in imported_names:
                    if not hasattr(mod, name):
                        raise ImportError(f"Module '{module}' has no attribute '{name}'")
            else:
                raise ValueError(f"Unknown import statement: {test}")
            print(f"  ✓ {test}")
        except Exception as e:
            if required:
                print(f"  X {test} - {e}")
                failed_imports.append((test, str(e)))
            else:
                print(f"  !  {test} - {e} (optional)")
                optional_failed += 1
    
    if optional_failed > 0:
        print(f"  ℹ️  {optional_failed} optional imports failed (likely missing assimulo)")
    
    return len(failed_imports) == 0, failed_imports

def run_reactor_tests():
    """Run the reactor tests."""
    print("\n Running reactor tests...")
    
    test_dir = Path("tests/integration")
    if not test_dir.exists():
        print(f"  X Test directory {test_dir} not found")
        return False, "Test directory not found"
    
    success, stdout, stderr = run_command(
        [sys.executable, "reactor_tests.py"],
        cwd=test_dir
    )
    
    if success:
        print("  ✓ Reactor tests passed")
        return True, stdout
    else:
        print("  X Reactor tests failed")
        print(f"  Error: {stderr}")
        return False, stderr

def run_flowsheet_tests():
    """Run the flowsheet tests.""" 
    print("\n Running flowsheet tests...")
    
    test_dir = Path("tests/Flowsheet")
    if not test_dir.exists():
        print(f"  X Test directory {test_dir} not found")
        return False, "Test directory not found"
    
    success, stdout, stderr = run_command(
        [sys.executable, "flowsheet_tests.py"],
        cwd=test_dir
    )
    
    if success:
        print("  ✓ Flowsheet tests passed")
        return True, stdout
    else:
        print("  X Flowsheet tests failed")
        print(f"  Error: {stderr}")
        return False, stderr

def run_pytest_tests():
    """Run tests using pytest if available."""
    print("\n Running pytest tests...")
    
    # Check if pytest is available
    try:
        import pytest
    except ImportError:
        print("  !  pytest not available, skipping pytest tests")
        return True, "pytest not available"
    
    success, stdout, stderr = run_command([
        sys.executable, "-m", "pytest", "tests/", "-v", "--tb=short"
    ])
    
    if success:
        print("  ✓ Pytest tests passed")
        return True, stdout
    else:
        print("  X Pytest tests failed")
        return False, stderr

def test_installation():
    """Test that the package is properly installed."""
    print("\n Testing package installation...")
    
    try:
        import pkg_resources
        pkg = pkg_resources.get_distribution('PharmaPy')
        print(f"  ✓ Package installed: {pkg.project_name} v{pkg.version}")
        return True, f"PharmaPy v{pkg.version}"
    except Exception as e:
        print(f"  X Package installation check failed: {e}")
        return False, str(e)

def main():
    parser = argparse.ArgumentParser(description='Run PharmaPy tests')
    parser.add_argument('--skip-imports', action='store_true', 
                       help='Skip import tests')
    parser.add_argument('--skip-reactor', action='store_true',
                       help='Skip reactor tests')
    parser.add_argument('--skip-flowsheet', action='store_true',
                       help='Skip flowsheet tests')
    parser.add_argument('--skip-pytest', action='store_true',
                       help='Skip pytest tests')
    parser.add_argument('--skip-installation', action='store_true',
                       help='Skip installation tests')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Verbose output')
    
    args = parser.parse_args()
    
    print(" PharmaPy Test Suite")
    print("=" * 50)
    
    start_time = time.time()
    results = {}
    
    # Run tests
    if not args.skip_imports:
        success, details = test_import()
        results['imports'] = (success, details)
    
    if not args.skip_installation:
        success, details = test_installation()
        results['installation'] = (success, details)
        
    if not args.skip_reactor:
        success, details = run_reactor_tests()
        results['reactor'] = (success, details)
    
    if not args.skip_flowsheet:
        success, details = run_flowsheet_tests()
        results['flowsheet'] = (success, details)
        
    if not args.skip_pytest:
        success, details = run_pytest_tests()
        results['pytest'] = (success, details)
    
    # Summary
    end_time = time.time()
    duration = end_time - start_time
    
    print(f"\n Test Summary (Duration: {duration:.2f}s)")
    print("=" * 50)
    
    total_tests = len(results)
    passed_tests = sum(1 for success, _ in results.values() if success)
    failed_tests = total_tests - passed_tests
    
    for test_name, (success, details) in results.items():
        status = "✓ PASS" if success else "X FAIL"
        print(f"  {test_name.upper():<12} {status}")
        
        if not success and args.verbose:
            print(f"    Details: {details}")
    
    print(f"\nResults: {passed_tests}/{total_tests} tests passed")
    
    if failed_tests > 0:
        print(f"X {failed_tests} test(s) failed")
        sys.exit(1)
    else:
        print(" All tests passed!")
        sys.exit(0)

if __name__ == "__main__":
    main()
