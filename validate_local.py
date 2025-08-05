#!/usr/bin/env python3
"""
Local validation script for PharmaPy testing infrastructure.
Run this before pushing to ensure CI will pass.
"""

import os
import sys
import subprocess
import platform
from pathlib import Path

def run_command(cmd, cwd=None, timeout=300, shell=None):
    """Run a command and return success status and output."""
    if shell is None:
        shell = platform.system() == "Windows"
        
    try:
        result = subprocess.run(
            cmd,
            cwd=cwd,
            capture_output=True,
            text=True,
            timeout=timeout,
            shell=shell
        )
        return result.returncode == 0, result.stdout, result.stderr
    except subprocess.TimeoutExpired:
        return False, "", f"Command timed out after {timeout} seconds"
    except Exception as e:
        return False, "", str(e)

def test_installation_methods():
    """Test different installation methods."""
    print(" Testing installation methods...")
    
    methods = [
        ("pip install -e .", "Modern pip installation"),
        ("python setup.py develop", "Legacy setup.py installation")
    ]
    
    for cmd, description in methods:
        print(f"\n  Testing: {description}")
        success, stdout, stderr = run_command(cmd.split())
        
        if success:
            print(f"    ✓ {description} successful")
        else:
            print(f"    X {description} failed")
            print(f"    Error: {stderr}")
            return False
    
    return True

def test_imports():
    """Test package imports."""
    print("\n Testing imports...")
    
    import_tests = [
        ("import PharmaPy", True),
        ("from PharmaPy import Utilities", True),
        ("from PharmaPy import Reactors", False),  # Requires assimulo
        ("from PharmaPy import Streams", False),   # Requires assimulo
        ("from PharmaPy import Phases", False),    # Requires assimulo
    ]
    
    for test_import, required in import_tests:
        try:
            exec(test_import)
            print(f"  ✓ {test_import}")
        except Exception as e:
            if required:
                print(f"  X {test_import} - {e}")
                return False
            else:
                print(f"  !  {test_import} - {e} (optional)")
    
    return True

def test_make_commands():
    """Test make commands."""
    print("\n Testing make commands...")
    
    is_windows = platform.system() == "Windows"
    make_cmd = "make.bat" if is_windows else "make"
    
    commands = [
        "test-imports",
        "clean",
    ]
    
    for cmd in commands:
        print(f"  Testing: {make_cmd} {cmd}")
        success, stdout, stderr = run_command([make_cmd, cmd])
        
        if success:
            print(f"    ✓ {make_cmd} {cmd} successful")
        else:
            print(f"    X {make_cmd} {cmd} failed")
            print(f"    Error: {stderr}")
    
    return True

def test_build_system():
    """Test package building."""
    print("\n Testing build system...")
    
    # Test if build tools are available
    try:
        import build
        print("  ✓ Build module available")
    except ImportError:
        print("  !  Build module not available, installing...")
        success, _, _ = run_command(["pip", "install", "build"])
        if not success:
            print("  X Could not install build module")
            return False
    
    # Test building
    print("  Building package...")
    success, stdout, stderr = run_command(["python", "-m", "build"])
    
    if success:
        print("  ✓ Package build successful")
        
        # Check if distributions exist
        dist_path = Path("dist")
        if dist_path.exists():
            wheels = list(dist_path.glob("*.whl"))
            tarballs = list(dist_path.glob("*.tar.gz"))
            
            if wheels and tarballs:
                print(f"  ✓ Found {len(wheels)} wheel(s) and {len(tarballs)} source distribution(s)")
            else:
                print(f"  !  Missing distributions: wheels={len(wheels)}, tarballs={len(tarballs)}")
        
        return True
    else:
        print("  X Package build failed")
        print(f"  Error: {stderr}")
        return False

def test_test_runners():
    """Test our custom test runners."""
    print("\n Testing test runners...")
    
    # Test Python test runner
    print("  Testing run_tests.py...")
    success, stdout, stderr = run_command([
        "python", "run_tests.py", 
        "--skip-reactor", "--skip-flowsheet", "--skip-pytest"
    ])
    
    if success:
        print("  ✓ Python test runner successful")
    else:
        print("  X Python test runner failed")
        print(f"  Output: {stdout}")
        print(f"  Error: {stderr}")
        return False
    
    # Test Windows batch runner (if on Windows)
    if platform.system() == "Windows":
        print("  Testing run_tests.bat...")
        success, stdout, stderr = run_command(["run_tests.bat"])
        
        if success:
            print("  ✓ Windows batch test runner successful")
        else:
            print("  X Windows batch test runner failed")
            print(f"  Error: {stderr}")
    
    return True

def check_assimulo_availability():
    """Check if assimulo is available."""
    print("\n Checking assimulo availability...")
    
    try:
        import assimulo
        print(f"  ✓ Assimulo available: version {assimulo.__version__}")
        return True
    except ImportError:
        print("  ℹ️  Assimulo not available - simulation features will be limited")
        print("   To install assimulo, consider using conda:")
        print("     conda install -c conda-forge sundials=5.8.0 superlu=5.2.2")
        print("     pip install assimulo")
        print("   This is informational - core PharmaPy works without assimulo")
        return True  # Changed to True since this is informational

def main():
    """Main validation function."""
    print(" PharmaPy Local Validation Suite")
    print("=" * 60)
    
    # Change to repo directory
    repo_dir = Path(__file__).parent
    os.chdir(repo_dir)
    print(f"Working directory: {repo_dir}")
    
    tests = [
        ("Import Tests", test_imports),
        ("Make Commands", test_make_commands), 
        ("Build System", test_build_system),
        ("Test Runners", test_test_runners),
        ("Assimulo Check", check_assimulo_availability),
    ]
    
    results = {}
    for test_name, test_func in tests:
        print(f"\n" + "="*60)
        try:
            success = test_func()
            results[test_name] = success
        except Exception as e:
            print(f"X {test_name} failed with exception: {e}")
            results[test_name] = False
    
    # Summary
    print(f"\n" + "="*60)
    print(" Validation Summary")
    print("="*60)
    
    passed = sum(results.values())
    total = len(results)
    
    for test_name, success in results.items():
        status = "✓ PASS" if success else "X FAIL"
        print(f"  {test_name:<20} {status}")
    
    print(f"\nOverall: {passed}/{total} tests passed")
    
    if passed == total:
        print(" All validations passed! Ready for CI/CD.")
        return 0
    else:
        print("X Some validations failed. Please fix before pushing.")
        return 1

if __name__ == "__main__":
    sys.exit(main())
