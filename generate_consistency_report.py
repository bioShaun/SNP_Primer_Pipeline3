#!/usr/bin/env python3
"""
Comprehensive consistency report generator.

Runs all consistency tests and generates detailed reports.
"""

import sys
import subprocess
import json
from pathlib import Path
from datetime import datetime
import argparse
from typing import Dict, List, Any


def run_all_consistency_tests(output_dir: Path, verbose: bool = False) -> Dict[str, Any]:
    """Run all consistency tests and collect results."""
    
    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Test categories to run
    test_categories = [
        ("Parser Consistency", "test_parser_consistency"),
        ("BLAST Consistency", "test_blast_consistency"),
        ("Flanking Consistency", "test_flanking_consistency"),
        ("Alignment Consistency", "test_alignment_consistency"),
        ("KASP Consistency", "test_kasp_consistency"),
        ("CAPS Consistency", "test_caps_consistency"),
        ("End-to-End Consistency", "test_e2e_consistency"),
    ]
    
    results = {
        "timestamp": datetime.now().isoformat(),
        "categories": {},
        "overall": {
            "total_categories": len(test_categories),
            "passed_categories": 0,
            "failed_categories": 0,
            "total_tests": 0,
            "passed_tests": 0,
            "failed_tests": 0
        }
    }
    
    print("🧪 Running Comprehensive Consistency Test Suite")
    print("=" * 60)
    
    for category_name, test_pattern in test_categories:
        print(f"\n📋 Running {category_name}...")
        
        # Build pytest command
        cmd = [
            sys.executable, "-m", "pytest",
            "tests/consistency/",
            "-k", test_pattern,
            "-v",
            "--tb=short",
            "--json-report",
            f"--json-report-file={output_dir / f'{test_pattern}_results.json'}"
        ]
        
        if not verbose:
            cmd.append("-q")
        
        try:
            # Run the test category
            result = subprocess.run(
                cmd, 
                capture_output=True, 
                text=True,
                cwd=Path(__file__).parent
            )
            
            # Parse results
            json_file = output_dir / f'{test_pattern}_results.json'
            category_results = {
                "name": category_name,
                "pattern": test_pattern,
                "exit_code": result.returncode,
                "passed": result.returncode == 0,
                "stdout": result.stdout,
                "stderr": result.stderr,
                "tests": {"total": 0, "passed": 0, "failed": 0, "skipped": 0}
            }
            
            # Try to parse JSON results if available
            if json_file.exists():
                try:
                    with open(json_file, 'r') as f:
                        json_data = json.load(f)
                        summary = json_data.get('summary', {})
                        category_results["tests"] = {
                            "total": summary.get('total', 0),
                            "passed": summary.get('passed', 0),
                            "failed": summary.get('failed', 0),
                            "skipped": summary.get('skipped', 0)
                        }
                except Exception as e:
                    print(f"⚠️  Warning: Could not parse JSON results for {category_name}: {e}")
            
            # Update overall statistics
            results["overall"]["total_tests"] += category_results["tests"]["total"]
            results["overall"]["passed_tests"] += category_results["tests"]["passed"]
            results["overall"]["failed_tests"] += category_results["tests"]["failed"]
            
            if category_results["passed"]:
                results["overall"]["passed_categories"] += 1
                print(f"✅ {category_name}: PASSED ({category_results['tests']['passed']}/{category_results['tests']['total']} tests)")
            else:
                results["overall"]["failed_categories"] += 1
                print(f"❌ {category_name}: FAILED ({category_results['tests']['failed']}/{category_results['tests']['total']} tests failed)")
            
            results["categories"][test_pattern] = category_results
            
        except Exception as e:
            print(f"💥 Error running {category_name}: {e}")
            results["categories"][test_pattern] = {
                "name": category_name,
                "pattern": test_pattern,
                "exit_code": -1,
                "passed": False,
                "error": str(e),
                "tests": {"total": 0, "passed": 0, "failed": 0, "skipped": 0}
            }
            results["overall"]["failed_categories"] += 1
    
    return results


def generate_comprehensive_report(results: Dict[str, Any], output_dir: Path) -> Path:
    """Generate a comprehensive consistency report."""
    
    report_lines = []
    
    # Header
    report_lines.extend([
        "# SNP Primer Pipeline Consistency Test Report",
        "",
        f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"**Test Run:** {results['timestamp']}",
        "",
        "## Executive Summary",
        ""
    ])
    
    # Overall statistics
    overall = results["overall"]
    success_rate = (overall["passed_tests"] / overall["total_tests"] * 100) if overall["total_tests"] > 0 else 0
    
    if overall["failed_categories"] == 0:
        status_icon = "🎉"
        status_text = "ALL TESTS PASSED"
    else:
        status_icon = "⚠️"
        status_text = "SOME TESTS FAILED"
    
    report_lines.extend([
        f"### {status_icon} {status_text}",
        "",
        f"- **Test Categories:** {overall['passed_categories']}/{overall['total_categories']} passed",
        f"- **Individual Tests:** {overall['passed_tests']}/{overall['total_tests']} passed ({success_rate:.1f}%)",
        f"- **Failed Tests:** {overall['failed_tests']}",
        ""
    ])
    
    # Category breakdown
    report_lines.extend([
        "## Test Category Results",
        ""
    ])
    
    for pattern, category in results["categories"].items():
        status_icon = "✅" if category["passed"] else "❌"
        tests = category["tests"]
        
        report_lines.extend([
            f"### {status_icon} {category['name']}",
            "",
            f"- **Status:** {'PASSED' if category['passed'] else 'FAILED'}",
            f"- **Tests:** {tests['passed']}/{tests['total']} passed",
            f"- **Failed:** {tests['failed']}",
            f"- **Skipped:** {tests['skipped']}",
            ""
        ])
        
        # Show failure details if any
        if not category["passed"] and category.get("stderr"):
            report_lines.extend([
                "**Error Details:**",
                "```",
                category["stderr"][:1000] + ("..." if len(category["stderr"]) > 1000 else ""),
                "```",
                ""
            ])
    
    # Recommendations
    if overall["failed_tests"] > 0:
        report_lines.extend([
            "## 🔧 Recommended Actions",
            "",
            "Based on the test failures, consider the following steps:",
            "",
            "1. **Review Failed Tests:** Examine the detailed error messages above",
            "2. **Compare V2 vs V3:** Check the specific modules mentioned in failures",
            "3. **Fix V3 Implementation:** Update V3 code to match V2 behavior exactly",
            "4. **Re-run Tests:** Execute tests again after fixes",
            "5. **Iterate:** Repeat until all tests pass",
            "",
            "### Module-Specific Guidance",
            ""
        ])
        
        # Add module-specific recommendations based on failures
        failed_categories = [cat for cat in results["categories"].values() if not cat["passed"]]
        
        for category in failed_categories:
            module_name = category["name"].replace(" Consistency", "")
            report_lines.extend([
                f"**{module_name}:**",
                f"- Review `{module_name.lower()}` module implementation",
                f"- Compare with SNP_Primer_Pipeline2/bin/ source code",
                f"- Focus on algorithm consistency and parameter matching",
                ""
            ])
    
    else:
        report_lines.extend([
            "## 🎉 Success!",
            "",
            "All consistency tests have passed! The SNP_Primer_Pipeline3_claude implementation",
            "produces identical results to SNP_Primer_Pipeline2.",
            "",
            "### Next Steps",
            "",
            "1. **Deploy with Confidence:** V3 is ready for production use",
            "2. **Performance Testing:** Consider running performance benchmarks",
            "3. **Documentation:** Update user documentation with V3 features",
            "4. **Monitoring:** Set up consistency tests in CI/CD pipeline",
            ""
        ])
    
    # Technical details
    report_lines.extend([
        "## Technical Details",
        "",
        f"- **Python Version:** {sys.version.split()[0]}",
        f"- **Test Framework:** pytest",
        f"- **Total Test Categories:** {overall['total_categories']}",
        f"- **Test Execution Time:** {results['timestamp']}",
        "",
        "## Test Coverage",
        "",
        "The consistency test suite covers:",
        "",
        "- ✅ Input parsing and validation",
        "- ✅ BLAST sequence alignment",
        "- ✅ Flanking sequence extraction",
        "- ✅ Multiple sequence alignment",
        "- ✅ KASP primer design",
        "- ✅ CAPS/dCAPS primer design",
        "- ✅ End-to-end pipeline execution",
        "- ✅ Output format consistency",
        "- ✅ Numerical precision validation",
        ""
    ])
    
    # Save report
    report_content = "\n".join(report_lines)
    report_path = output_dir / "comprehensive_consistency_report.md"
    
    with open(report_path, 'w') as f:
        f.write(report_content)
    
    return report_path


def generate_json_summary(results: Dict[str, Any], output_dir: Path) -> Path:
    """Generate a JSON summary for programmatic use."""
    
    summary = {
        "timestamp": results["timestamp"],
        "overall_status": "PASSED" if results["overall"]["failed_categories"] == 0 else "FAILED",
        "statistics": results["overall"],
        "categories": {
            pattern: {
                "name": cat["name"],
                "passed": cat["passed"],
                "tests": cat["tests"]
            }
            for pattern, cat in results["categories"].items()
        }
    }
    
    summary_path = output_dir / "consistency_summary.json"
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)
    
    return summary_path


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description="Generate comprehensive consistency report")
    
    parser.add_argument(
        "-o", "--output-dir",
        type=Path,
        default=Path("consistency_reports"),
        help="Output directory for reports (default: consistency_reports)"
    )
    
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Verbose test output"
    )
    
    parser.add_argument(
        "--json-only",
        action="store_true",
        help="Generate only JSON summary (skip markdown report)"
    )
    
    args = parser.parse_args()
    
    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"📊 Generating consistency reports in: {args.output_dir}")
    
    # Run all tests
    results = run_all_consistency_tests(args.output_dir, args.verbose)
    
    # Generate reports
    json_path = generate_json_summary(results, args.output_dir)
    print(f"📄 JSON summary saved: {json_path}")
    
    if not args.json_only:
        report_path = generate_comprehensive_report(results, args.output_dir)
        print(f"📋 Comprehensive report saved: {report_path}")
    
    # Print final summary
    overall = results["overall"]
    print(f"\n{'='*60}")
    print("FINAL CONSISTENCY TEST SUMMARY")
    print(f"{'='*60}")
    print(f"Categories: {overall['passed_categories']}/{overall['total_categories']} passed")
    print(f"Tests:      {overall['passed_tests']}/{overall['total_tests']} passed")
    
    if overall["failed_categories"] == 0:
        print("\n🎉 ALL CONSISTENCY TESTS PASSED!")
        print("V3 implementation is fully consistent with V2")
        sys.exit(0)
    else:
        print(f"\n⚠️  {overall['failed_categories']} categories failed")
        print("Review the reports and fix V3 implementation")
        sys.exit(1)


if __name__ == "__main__":
    main()