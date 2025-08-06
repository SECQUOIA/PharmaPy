#!/usr/bin/env python3
"""
Local validation script for PharmaPy testing infrastructure.
Run this before pushing to ensure CI will pass.

This script extracts and runs the same validation steps as GitHub Actions CI,
ensuring perfect sync between local and remote validation.
"""

import os
import sys
import subprocess
import yaml
import platform
from pathlib import Path


def run_command(
    cmd, shell: bool = None, cwd: str = None, timeout: int = 300
) -> tuple[bool, str, str]:
    """
    Run a shell command (string or list) and return (success, stdout, stderr).
    Automatically detects shell usage based on platform and command type.
    """
    if shell is None:
        # Use shell=True for string commands, False for list commands
        shell = isinstance(cmd, str)
    try:
        result = subprocess.run(
            cmd,
            shell=shell,
            capture_output=True,
            text=True,
            timeout=timeout,
            cwd=cwd,
        )
        return result.returncode == 0, result.stdout, result.stderr
    except subprocess.TimeoutExpired:
        return False, "", f"Command timed out after {timeout} seconds"
    except Exception as e:
        return False, "", str(e)


def parse_ci_workflow(workflow_path: Path) -> dict:
    """Parse GitHub Actions workflow YAML file."""
    try:
        with open(workflow_path, "r") as f:
            return yaml.safe_load(f)
    except Exception as e:
        print(f"Error parsing CI workflow: {e}")
        return {}


def extract_quality_steps(workflow: dict) -> list[dict]:
    """Extract the quality and documentation steps from CI workflow."""
    try:
        quality_job = workflow["jobs"]["quality-and-docs"]
        steps = quality_job["steps"]

        # Filter to the validation steps we want to run locally
        validation_steps = []
        for step in steps:
            name = step.get("name", "")
            if any(
                keyword in name.lower()
                for keyword in [
                    "linting",
                    "black",
                    "coverage",
                    "documentation",
                    "pandoc",
                ]
            ):
                validation_steps.append(step)

        return validation_steps
    except Exception as e:
        print(f"Error extracting quality steps: {e}")
        return []


def adapt_command_for_local(run_commands: str) -> list[str]:
    """
    Adapt CI run commands for local execution.
    Remove CI-specific parts and adapt for local environment.
    """
    commands = []
    for line in run_commands.strip().split("\n"):
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        # Skip CI-specific commands but keep informational echoes
        if any(skip in line for skip in ["pip install -e", "pip install -r"]):
            continue

        # For conda install commands, convert to availability checks
        if "conda install" in line and "pandoc" in line:
            line = "pandoc --version"  # Check if pandoc is available
        elif "conda install" in line:
            continue  # Skip other conda installs

        # Adapt paths and commands for local execution
        if "cd doc" in line:
            # Keep track that we need to change directory for subsequent commands
            continue
        elif "make html" in line:
            # Use the fixed Makefile approach
            line = "cd doc && make html"
        elif line.startswith("echo "):
            # Keep echo commands for context
            pass

        commands.append(line)

    return commands


def run_validation_step(step: dict) -> bool:
    """Run a single validation step locally."""
    name = step.get("name", "Unknown step")
    run_commands = step.get("run", "")

    print(f"\n[RUNNING] {name}")
    print("-" * 50)

    if not run_commands:
        print("  [SKIP] No run commands found")
        return True

    # Extract and adapt commands
    commands = adapt_command_for_local(run_commands)

    if not commands:
        print("  [SKIP] No applicable commands for local execution")
        return True

    # Run each command
    overall_success = True
    for cmd in commands:
        print(f"  -> {cmd}")
        success, stdout, stderr = run_command(cmd)

        if success:
            print("    [OK] Success")
            if stdout.strip():
                # Show condensed output for successful commands
                lines = stdout.strip().split("\n")
                if len(lines) <= 3:
                    for line in lines:
                        if line.strip():
                            print(f"    OUTPUT: {line}")
                else:
                    print(f"    OUTPUT: {lines[0]}")
                    if len(lines) > 1:
                        print(f"    OUTPUT: ... ({len(lines)-1} more lines)")
        else:
            print("    [FAIL] Failed")
            if stderr:
                print(f"    ERROR: {stderr[:200]}...")
            if stdout:
                print(f"    OUTPUT: {stdout[:200]}...")

            # For linting and formatting, failure is informational, not blocking
            if any(keyword in name.lower() for keyword in ["linting", "black"]):
                print("    [INFO] Non-blocking: continuing...")
            else:
                overall_success = False

    return overall_success


def main():
    """Main function to run CI validation locally."""
    print("PharmaPy CI-Based Local Validation")
    print("=" * 60)
    print("This script runs the same validation steps as GitHub Actions CI")
    print("No more maintaining separate validation logic!")

    # Change to repo directory
    repo_dir = Path(__file__).parent
    os.chdir(repo_dir)
    print(f"Working directory: {repo_dir}")

    # Parse CI workflow
    ci_path = repo_dir / ".github" / "workflows" / "ci.yml"
    if not ci_path.exists():
        print(f"[ERROR] CI workflow not found at {ci_path}")
        return 1

    print(f"Loading CI workflow from {ci_path}")
    workflow = parse_ci_workflow(ci_path)

    if not workflow:
        print("[ERROR] Failed to parse CI workflow")
        return 1

    # Extract quality steps
    quality_steps = extract_quality_steps(workflow)
    if not quality_steps:
        print("[ERROR] No quality steps found in CI workflow")
        return 1

    print(f"Found {len(quality_steps)} validation steps to run")

    # Prerequisites check
    print("\nChecking prerequisites...")
    prereq_commands = [
        "python --version",
        "pip --version",
        "flake8 --version",
        "black --version",
        "pytest --version",
    ]

    missing_tools = []
    for cmd in prereq_commands:
        success, stdout, stderr = run_command(cmd)
        tool = cmd.split()[0]
        if success:
            version = stdout.strip().split("\n")[0]
            print(f"  [OK] {tool}: {version}")
        else:
            print(f"  [MISSING] {tool}: Not found")
            missing_tools.append(tool)

    if missing_tools:
        print(f"\nWarning: Missing tools: {', '.join(missing_tools)}")
        print("   Install with: pip install flake8 black pytest pytest-cov")
        print("   Or use: pip install -r requirements-dev.txt")

    # Run validation steps
    print(f"\nRunning {len(quality_steps)} validation steps...")
    results = {}

    for step in quality_steps:
        step_name = step.get("name", "Unknown")
        try:
            success = run_validation_step(step)
            results[step_name] = success
        except Exception as e:
            print(f"[EXCEPTION] Exception in {step_name}: {e}")
            results[step_name] = False

    # Summary
    print(f"\nValidation Summary")
    print("=" * 60)

    passed = sum(results.values())
    total = len(results)

    for step_name, success in results.items():
        status = "[PASS]" if success else "[FAIL]"
        print(f"  {status} {step_name}")

    print(f"\nOverall: {passed}/{total} steps passed")

    if passed == total:
        print("All validations passed! Ready for CI/CD.")
        return 0
    else:
        print("Some validations failed. Please fix before pushing.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
