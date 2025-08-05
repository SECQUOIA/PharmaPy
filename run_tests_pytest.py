"""
Test runner script for PharmaPy with different test configurations.
"""
import sys
import subprocess
import argparse
from pathlib import Path


def run_tests(test_type="all", verbose=False, coverage=False, parallel=False, markers=None):
    """
    Run PharmaPy tests with different configurations.
    
    Args:
        test_type (str): Type of tests to run ('unit', 'integration', 'all')
        verbose (bool): Enable verbose output
        coverage (bool): Enable coverage reporting
        parallel (bool): Enable parallel test execution
        markers (str): Pytest markers to filter tests
    """
    
    # Base pytest command
    cmd = ["pytest"]
    
    # Add test paths based on test_type
    if test_type == "unit":
        cmd.append("tests/unit/")
    elif test_type == "integration":
        cmd.append("tests/integration/")
    elif test_type == "all":
        cmd.append("tests/")
    else:
        cmd.append(test_type)  # Custom path
    
    # Add options
    if verbose:
        cmd.append("-v")
    
    if coverage:
        cmd.extend(["--cov=PharmaPy", "--cov-report=html", "--cov-report=term-missing"])
    
    if parallel:
        cmd.extend(["-n", "auto"])
    
    if markers:
        cmd.extend(["-m", markers])
    
    # Run the tests
    try:
        result = subprocess.run(cmd, check=False)
        return result.returncode
    except FileNotFoundError:
        print("Error: pytest not found. Please install pytest:")
        print("conda install pytest pytest-cov pytest-xdist")
        return 1


def main():
    """Main entry point for test runner."""
    parser = argparse.ArgumentParser(description="PharmaPy Test Runner")
    
    parser.add_argument(
        "--type", 
        choices=["unit", "integration", "all"],
        default="all",
        help="Type of tests to run"
    )
    
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable verbose output"
    )
    
    parser.add_argument(
        "--coverage", "-c",
        action="store_true",
        help="Enable coverage reporting"
    )
    
    parser.add_argument(
        "--parallel", "-p",
        action="store_true",
        help="Enable parallel test execution"
    )
    
    parser.add_argument(
        "--markers", "-m",
        help="Pytest markers to filter tests (e.g., 'not slow', 'unit', 'assimulo')"
    )
    
    parser.add_argument(
        "--fast",
        action="store_true",
        help="Run only fast tests (excludes slow and integration tests)"
    )
    
    parser.add_argument(
        "--assimulo-only",
        action="store_true",
        help="Run only assimulo-related tests"
    )
    
    parser.add_argument(
        "--no-assimulo",
        action="store_true",
        help="Skip assimulo-related tests"
    )
    
    args = parser.parse_args()
    
    # Set up markers based on arguments
    markers = args.markers
    
    if args.fast:
        markers = "not slow and not integration"
    elif args.assimulo_only:
        markers = "assimulo"
    elif args.no_assimulo:
        markers = "not assimulo"
    
    # Run tests
    exit_code = run_tests(
        test_type=args.type,
        verbose=args.verbose,
        coverage=args.coverage,
        parallel=args.parallel,
        markers=markers
    )
    
    sys.exit(exit_code)


if __name__ == "__main__":
    main()
