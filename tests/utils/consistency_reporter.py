#!/usr/bin/env python3
"""
Consistency reporter for generating test reports.

Generates detailed reports of consistency test results.
"""

from pathlib import Path
from typing import List, Dict, Any
from datetime import datetime
from .output_comparator import ComparisonResult


class ConsistencyReporter:
    """Generates consistency test reports."""
    
    def __init__(self, output_dir: Path):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.results: List[ComparisonResult] = []
        self.test_sections: Dict[str, List[ComparisonResult]] = {}
    
    def add_result(self, result: ComparisonResult) -> None:
        """Add a single comparison result."""
        self.results.append(result)
    
    def add_results(self, results: List[ComparisonResult]) -> None:
        """Add multiple comparison results."""
        self.results.extend(results)
    
    def add_section_results(self, section_name: str, results: List[ComparisonResult]) -> None:
        """Add results for a specific test section."""
        if section_name not in self.test_sections:
            self.test_sections[section_name] = []
        self.test_sections[section_name].extend(results)
        self.add_results(results)
    
    def generate_report(self) -> str:
        """Generate a comprehensive consistency report."""
        report_lines = []
        
        # Header
        report_lines.extend([
            "# SNP Primer Pipeline Consistency Test Report",
            "",
            f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            "",
            "## Summary",
            ""
        ])
        
        # Overall summary
        summary = self.get_summary()
        report_lines.extend([
            f"- **Total Tests:** {summary['total_tests']}",
            f"- **Passed:** {summary['passed_tests']} ({summary['success_rate']:.1f}%)",
            f"- **Failed:** {summary['failed_tests']} ({summary['failure_rate']:.1f}%)",
            ""
        ])
        
        # Section summaries
        if self.test_sections:
            report_lines.extend([
                "## Test Section Results",
                ""
            ])
            
            for section_name, section_results in self.test_sections.items():
                section_summary = self._get_section_summary(section_results)
                status_icon = "✅" if section_summary['failed_tests'] == 0 else "❌"
                
                report_lines.extend([
                    f"### {status_icon} {section_name}",
                    "",
                    f"- **Tests:** {section_summary['total_tests']}",
                    f"- **Passed:** {section_summary['passed_tests']} ({section_summary['success_rate']:.1f}%)",
                    f"- **Failed:** {section_summary['failed_tests']} ({section_summary['failure_rate']:.1f}%)",
                    ""
                ])
                
                # Show failed tests for this section
                failed_results = [r for r in section_results if not r.is_match]
                if failed_results:
                    report_lines.extend([
                        "**Failed Tests:**",
                        ""
                    ])
                    
                    for result in failed_results[:10]:  # Limit to first 10 failures
                        report_lines.append(f"- {result.message}")
                    
                    if len(failed_results) > 10:
                        report_lines.append(f"- ... and {len(failed_results) - 10} more failures")
                    
                    report_lines.append("")
        
        # Detailed failures
        failed_results = [r for r in self.results if not r.is_match]
        if failed_results:
            report_lines.extend([
                "## Detailed Failure Analysis",
                "",
                "The following tests failed and require attention:",
                ""
            ])
            
            # Group failures by category
            failure_categories = self._categorize_failures(failed_results)
            
            for category, failures in failure_categories.items():
                report_lines.extend([
                    f"### {category} ({len(failures)} failures)",
                    ""
                ])
                
                for failure in failures[:20]:  # Limit to first 20 per category
                    report_lines.extend([
                        f"**{failure.field_name}**",
                        f"- Expected: `{failure.expected_value}`",
                        f"- Actual: `{failure.actual_value}`"
                    ])
                    
                    if failure.tolerance is not None:
                        report_lines.append(f"- Tolerance: {failure.tolerance}")
                    
                    report_lines.append("")
                
                if len(failures) > 20:
                    report_lines.append(f"*... and {len(failures) - 20} more {category.lower()} failures*")
                    report_lines.append("")
        
        # Recommendations
        if failed_results:
            report_lines.extend([
                "## Recommendations",
                "",
                "Based on the test failures, consider the following actions:",
                ""
            ])
            
            recommendations = self._generate_recommendations(failed_results)
            for rec in recommendations:
                report_lines.append(f"- {rec}")
            
            report_lines.append("")
        
        # Success message
        if not failed_results:
            report_lines.extend([
                "## 🎉 All Tests Passed!",
                "",
                "The SNP_Primer_Pipeline3 implementation is fully consistent with SNP_Primer_Pipeline2.",
                "No code changes are required.",
                ""
            ])
        
        return "\n".join(report_lines)
    
    def save_report(self, filename: str = "consistency_report.md") -> Path:
        """Save the report to a file."""
        report_content = self.generate_report()
        report_path = self.output_dir / filename
        
        with open(report_path, 'w') as f:
            f.write(report_content)
        
        return report_path
    
    def save_detailed_failures(self, filename: str = "detailed_failures.txt") -> Path:
        """Save detailed failure information to a separate file."""
        failed_results = [r for r in self.results if not r.is_match]
        
        if not failed_results:
            return None
        
        failure_path = self.output_dir / filename
        
        with open(failure_path, 'w') as f:
            f.write(f"Detailed Failure Report - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write("=" * 80 + "\n\n")
            
            for i, result in enumerate(failed_results, 1):
                f.write(f"Failure #{i}: {result.field_name}\n")
                f.write(f"Expected: {result.expected_value}\n")
                f.write(f"Actual:   {result.actual_value}\n")
                
                if result.tolerance is not None:
                    f.write(f"Tolerance: {result.tolerance}\n")
                
                f.write(f"Message: {result.message}\n")
                f.write("-" * 40 + "\n\n")
        
        return failure_path
    
    def get_summary(self) -> Dict[str, Any]:
        """Get summary statistics."""
        total_tests = len(self.results)
        passed_tests = sum(1 for r in self.results if r.is_match)
        failed_tests = total_tests - passed_tests
        
        return {
            'total_tests': total_tests,
            'passed_tests': passed_tests,
            'failed_tests': failed_tests,
            'success_rate': (passed_tests / total_tests * 100) if total_tests > 0 else 100.0,
            'failure_rate': (failed_tests / total_tests * 100) if total_tests > 0 else 0.0
        }
    
    def _get_section_summary(self, results: List[ComparisonResult]) -> Dict[str, Any]:
        """Get summary statistics for a section."""
        total_tests = len(results)
        passed_tests = sum(1 for r in results if r.is_match)
        failed_tests = total_tests - passed_tests
        
        return {
            'total_tests': total_tests,
            'passed_tests': passed_tests,
            'failed_tests': failed_tests,
            'success_rate': (passed_tests / total_tests * 100) if total_tests > 0 else 100.0,
            'failure_rate': (failed_tests / total_tests * 100) if total_tests > 0 else 0.0
        }
    
    def _categorize_failures(self, failures: List[ComparisonResult]) -> Dict[str, List[ComparisonResult]]:
        """Categorize failures by type."""
        categories = {
            'Sequence Mismatches': [],
            'Numerical Differences': [],
            'Position/Coordinate Errors': [],
            'Count Mismatches': [],
            'Other Failures': []
        }
        
        for failure in failures:
            field_name = failure.field_name.lower()
            
            if 'seq' in field_name or 'sequence' in field_name:
                categories['Sequence Mismatches'].append(failure)
            elif failure.tolerance is not None or any(x in field_name for x in ['tm', 'gc', 'score', 'penalty', 'evalue']):
                categories['Numerical Differences'].append(failure)
            elif any(x in field_name for x in ['start', 'end', 'position', 'coordinate']):
                categories['Position/Coordinate Errors'].append(failure)
            elif 'count' in field_name:
                categories['Count Mismatches'].append(failure)
            else:
                categories['Other Failures'].append(failure)
        
        # Remove empty categories
        return {k: v for k, v in categories.items() if v}
    
    def _generate_recommendations(self, failures: List[ComparisonResult]) -> List[str]:
        """Generate recommendations based on failure patterns."""
        recommendations = []
        
        # Analyze failure patterns
        failure_categories = self._categorize_failures(failures)
        
        if 'Sequence Mismatches' in failure_categories:
            recommendations.append(
                "Review sequence handling in core/parser.py and core/alignment.py for FASTA parsing and sequence processing"
            )
        
        if 'Numerical Differences' in failure_categories:
            recommendations.append(
                "Check numerical calculations in primers/kasp.py and primers/caps.py, especially Tm and GC content calculations"
            )
        
        if 'Position/Coordinate Errors' in failure_categories:
            recommendations.append(
                "Verify coordinate calculations in core/blast.py, particularly flanking region extraction and position mapping"
            )
        
        if 'Count Mismatches' in failure_categories:
            recommendations.append(
                "Check filtering logic in BLAST processing and primer design modules"
            )
        
        # Module-specific recommendations
        module_failures = {}
        for failure in failures:
            field_name = failure.field_name.lower()
            
            if 'kasp' in field_name:
                module_failures.setdefault('KASP', []).append(failure)
            elif 'caps' in field_name:
                module_failures.setdefault('CAPS', []).append(failure)
            elif 'blast' in field_name:
                module_failures.setdefault('BLAST', []).append(failure)
            elif 'alignment' in field_name:
                module_failures.setdefault('Alignment', []).append(failure)
        
        for module, module_failure_list in module_failures.items():
            if len(module_failure_list) > 5:  # Significant number of failures
                recommendations.append(
                    f"Focus on {module} module - {len(module_failure_list)} failures detected"
                )
        
        if not recommendations:
            recommendations.append("Review all modules systematically to identify the source of inconsistencies")
        
        return recommendations
    
    def print_summary(self) -> None:
        """Print a brief summary to console."""
        summary = self.get_summary()
        
        print(f"\n{'='*60}")
        print("CONSISTENCY TEST SUMMARY")
        print(f"{'='*60}")
        print(f"Total Tests: {summary['total_tests']}")
        print(f"Passed:      {summary['passed_tests']} ({summary['success_rate']:.1f}%)")
        print(f"Failed:      {summary['failed_tests']} ({summary['failure_rate']:.1f}%)")
        
        if summary['failed_tests'] > 0:
            print(f"\n❌ {summary['failed_tests']} tests failed - review required")
        else:
            print(f"\n✅ All tests passed - implementation is consistent!")
        
        print(f"{'='*60}\n")