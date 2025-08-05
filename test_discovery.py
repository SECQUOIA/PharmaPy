#!/usr/bin/env python3
"""
Test discovery and execution wrapper for PharmaPy.
This script finds and runs all available tests in a structured way.
"""

import os
import sys
import unittest
import importlib
import traceback
from pathlib import Path

def discover_test_modules():
    """Discover all test modules in the tests directory."""
    test_modules = []
    tests_dir = Path("tests")
    
    if not tests_dir.exists():
        print("X Tests directory not found")
        return test_modules
    
    # Find all Python files that look like tests
    for test_file in tests_dir.rglob("*test*.py"):
        if test_file.name.startswith("__"):
            continue
            
        # Convert file path to module path
        relative_path = test_file.relative_to(Path("."))
        module_path = str(relative_path.with_suffix("")).replace(os.sep, ".")
        
        test_modules.append((module_path, test_file))
    
    return test_modules

def run_module_tests(module_path, test_file):
    """Run tests from a specific module."""
    print(f"\n Testing module: {module_path}")
    print(f"   File: {test_file}")
    
    try:
        # Add the directory containing the test to Python path
        test_dir = test_file.parent
        if str(test_dir) not in sys.path:
            sys.path.insert(0, str(test_dir))
        
        # Import the module
        module = importlib.import_module(module_path)
        
        # Create a test suite
        loader = unittest.TestLoader()
        suite = loader.loadTestsFromModule(module)
        
        # Run the tests
        runner = unittest.TextTestRunner(verbosity=2, stream=sys.stdout)
        result = runner.run(suite)
        
        # Check if tests passed
        if result.wasSuccessful():
            print(f"✓ {module_path}: All tests passed")
            return True, result
        else:
            print(f"X {module_path}: {len(result.failures)} failures, {len(result.errors)} errors")
            return False, result
            
    except Exception as e:
        print(f"X {module_path}: Failed to run - {e}")
        traceback.print_exc()
        return False, None

def run_direct_execution_tests():
    """Run tests that are designed for direct execution."""
    print("\n Running direct execution tests...")
    
    direct_test_files = [
        "tests/integration/reactor_tests.py",
        "tests/Flowsheet/flowsheet_tests.py"
    ]
    
    results = []
    for test_file in direct_test_files:
        if not Path(test_file).exists():
            print(f"!  {test_file} not found, skipping")
            continue
            
        print(f"\n▶️  Executing {test_file}")
        try:
            # Change to the test directory
            old_cwd = os.getcwd()
            test_dir = Path(test_file).parent
            os.chdir(test_dir)
            
            # Execute the test file using subprocess
            result = subprocess.run([sys.executable, Path(test_file).name], capture_output=True, text=True)
            if result.returncode == 0:
                print(f"✓ {test_file} executed successfully")
                print(result.stdout)
                results.append((test_file, True, None))
            else:
                print(f"X {test_file} failed with return code {result.returncode}")
                print(result.stdout)
                print(result.stderr)
                results.append((test_file, False, result.stderr))
            
        except Exception as e:
            print(f"X {test_file} failed: {e}")
            traceback.print_exc()
            results.append((test_file, False, e))
        finally:
            os.chdir(old_cwd)
    
    return results

def main():
    """Main test discovery and execution function."""
    print(" PharmaPy Test Discovery and Execution")
    print("=" * 60)
    
    # Test basic imports first
    print("\n Testing basic imports...")
    try:
        import PharmaPy
        print("✓ PharmaPy imported successfully")
    except Exception as e:
        print(f"X Failed to import PharmaPy: {e}")
        return 1
    
    # Discover and run unittest-based tests
    print("\n Discovering test modules...")
    test_modules = discover_test_modules()
    
    if not test_modules:
        print("!  No test modules discovered")
    else:
        print(f"Found {len(test_modules)} test modules:")
        for module_path, test_file in test_modules:
            print(f"  - {module_path}")
    
    unittest_results = []
    for module_path, test_file in test_modules:
        success, result = run_module_tests(module_path, test_file)
        unittest_results.append((module_path, success, result))
    
    # Run direct execution tests
    direct_results = run_direct_execution_tests()
    
    # Summary
    print("\n Test Summary")
    print("=" * 60)
    
    total_tests = 0
    passed_tests = 0
    
    print("\nUnittest-based tests:")
    for module_path, success, result in unittest_results:
        status = "✓ PASS" if success else "X FAIL"
        print(f"  {module_path:<40} {status}")
        total_tests += 1
        if success:
            passed_tests += 1
    
    print("\nDirect execution tests:")
    for test_file, success, error in direct_results:
        status = "✓ PASS" if success else "X FAIL"
        print(f"  {Path(test_file).name:<40} {status}")
        total_tests += 1
        if success:
            passed_tests += 1
    
    print(f"\nOverall: {passed_tests}/{total_tests} tests passed")
    
    if passed_tests == total_tests:
        print(" All tests passed!")
        return 0
    else:
        print(f"X {total_tests - passed_tests} test(s) failed")
        return 1

if __name__ == "__main__":
    sys.exit(main())
