#!/usr/bin/env python3
"""
Consistency test runner.

Runs the complete consistency test suite and generates reports.
"""

import sys
import subprocess
from pathlib import Path
import argparse


def run_consistency_tests(test_pattern="*", verbose=False, generate_report=True):
    """Run consistency tests with optional filtering."""
    
    # Get the project root directory
    project_root = Path(__file__).parent
    
    # Ensure we're in the right directory
    import os
    os.chdir(project_root)
    
    # Build pytest command
    cmd = [
        sys.executable, "-m", "pytest",
        "tests/consistency/",
        "-v" if verbose else "-q",
        "--tb=short",
        "-x",  # Stop on first failure for faster feedback
    ]
    
    # Add test pattern if specified
    if test_pattern != "*":
        cmd.extend(["-k", test_pattern])
    
    # Add HTML report generation
    if generate_report:
        cmd.extend(["--html=consistency_test_report.html", "--self-contained-html"])
    
    print(f"Running consistency tests...")
    print(f"Command: {' '.join(cmd)}")
    print("-" * 60)
    
    # Run the tests
    try:
        result = subprocess.run(cmd, capture_output=False, text=True)
        return result.returncode == 0
    except KeyboardInterrupt:
        print("\nTest execution interrupted by user")
        return False
    except Exception as e:
        print(f"Error running tests: {e}")
        return False


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description="Run SNP Primer Pipeline consistency tests")
    
    parser.add_argument(
        "-k", "--pattern",
        default="*",
        help="Test pattern to match (default: run all tests)"
    )
    
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Verbose output"
    )
    
    parser.add_argument(
        "--no-report",
        action="store_true",
        help="Skip HTML report generation"
    )
    
    parser.add_argument(
        "--parser-only",
        action="store_true",
        help="Run only parser consistency tests"
    )
    
    parser.add_argument(
        "--blast-only",
        action="store_true",
        help="Run only BLAST consistency tests"
    )
    
    args = parser.parse_args()
    
    # Determine test pattern
    test_pattern = args.pattern
    if args.parser_only:
        test_pattern = "test_parser_consistency"
    elif args.blast_only:
        test_pattern = "test_blast_consistency"
    
    # Run tests
    success = run_consistency_tests(
        test_pattern=test_pattern,
        verbose=args.verbose,
        generate_report=not args.no_report
    )
    
    if success:
        print("\n✅ All consistency tests passed!")
        if not args.no_report:
            print("📊 HTML report generated: consistency_test_report.html")
    else:
        print("\n❌ Some consistency tests failed!")
        print("Review the output above and fix V3 code as needed.")
        sys.exit(1)


if __name__ == "__main__":
    main()