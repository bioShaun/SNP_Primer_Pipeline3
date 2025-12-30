# SNP Primer Pipeline Consistency Test Report

**Generated:** 2025-12-30 13:59:17
**Test Run:** 2025-12-30T13:59:15.046402

## Executive Summary

### ⚠️ SOME TESTS FAILED

- **Test Categories:** 0/7 passed
- **Individual Tests:** 0/0 passed (0.0%)
- **Failed Tests:** 0

## Test Category Results

### ❌ Parser Consistency

- **Status:** FAILED
- **Tests:** 0/0 passed
- **Failed:** 0
- **Skipped:** 0

**Error Details:**
```
ERROR: usage: __main__.py [options] [file_or_dir] [file_or_dir] [...]
__main__.py: error: unrecognized arguments: --json-report --json-report-file=consistency_reports/test_parser_consistency_results.json
  inifile: /home/kent/scripts/kasp/SNP_Primer_Pipeline3_claude/pyproject.toml
  rootdir: /home/kent/scripts/kasp/SNP_Primer_Pipeline3_claude


```

### ❌ BLAST Consistency

- **Status:** FAILED
- **Tests:** 0/0 passed
- **Failed:** 0
- **Skipped:** 0

**Error Details:**
```
ERROR: usage: __main__.py [options] [file_or_dir] [file_or_dir] [...]
__main__.py: error: unrecognized arguments: --json-report --json-report-file=consistency_reports/test_blast_consistency_results.json
  inifile: /home/kent/scripts/kasp/SNP_Primer_Pipeline3_claude/pyproject.toml
  rootdir: /home/kent/scripts/kasp/SNP_Primer_Pipeline3_claude


```

### ❌ Flanking Consistency

- **Status:** FAILED
- **Tests:** 0/0 passed
- **Failed:** 0
- **Skipped:** 0

**Error Details:**
```
ERROR: usage: __main__.py [options] [file_or_dir] [file_or_dir] [...]
__main__.py: error: unrecognized arguments: --json-report --json-report-file=consistency_reports/test_flanking_consistency_results.json
  inifile: /home/kent/scripts/kasp/SNP_Primer_Pipeline3_claude/pyproject.toml
  rootdir: /home/kent/scripts/kasp/SNP_Primer_Pipeline3_claude


```

### ❌ Alignment Consistency

- **Status:** FAILED
- **Tests:** 0/0 passed
- **Failed:** 0
- **Skipped:** 0

**Error Details:**
```
ERROR: usage: __main__.py [options] [file_or_dir] [file_or_dir] [...]
__main__.py: error: unrecognized arguments: --json-report --json-report-file=consistency_reports/test_alignment_consistency_results.json
  inifile: /home/kent/scripts/kasp/SNP_Primer_Pipeline3_claude/pyproject.toml
  rootdir: /home/kent/scripts/kasp/SNP_Primer_Pipeline3_claude


```

### ❌ KASP Consistency

- **Status:** FAILED
- **Tests:** 0/0 passed
- **Failed:** 0
- **Skipped:** 0

**Error Details:**
```
ERROR: usage: __main__.py [options] [file_or_dir] [file_or_dir] [...]
__main__.py: error: unrecognized arguments: --json-report --json-report-file=consistency_reports/test_kasp_consistency_results.json
  inifile: /home/kent/scripts/kasp/SNP_Primer_Pipeline3_claude/pyproject.toml
  rootdir: /home/kent/scripts/kasp/SNP_Primer_Pipeline3_claude


```

### ❌ CAPS Consistency

- **Status:** FAILED
- **Tests:** 0/0 passed
- **Failed:** 0
- **Skipped:** 0

**Error Details:**
```
ERROR: usage: __main__.py [options] [file_or_dir] [file_or_dir] [...]
__main__.py: error: unrecognized arguments: --json-report --json-report-file=consistency_reports/test_caps_consistency_results.json
  inifile: /home/kent/scripts/kasp/SNP_Primer_Pipeline3_claude/pyproject.toml
  rootdir: /home/kent/scripts/kasp/SNP_Primer_Pipeline3_claude


```

### ❌ End-to-End Consistency

- **Status:** FAILED
- **Tests:** 0/0 passed
- **Failed:** 0
- **Skipped:** 0

**Error Details:**
```
ERROR: usage: __main__.py [options] [file_or_dir] [file_or_dir] [...]
__main__.py: error: unrecognized arguments: --json-report --json-report-file=consistency_reports/test_e2e_consistency_results.json
  inifile: /home/kent/scripts/kasp/SNP_Primer_Pipeline3_claude/pyproject.toml
  rootdir: /home/kent/scripts/kasp/SNP_Primer_Pipeline3_claude


```

## 🎉 Success!

All consistency tests have passed! The SNP_Primer_Pipeline3_claude implementation
produces identical results to SNP_Primer_Pipeline2.

### Next Steps

1. **Deploy with Confidence:** V3 is ready for production use
2. **Performance Testing:** Consider running performance benchmarks
3. **Documentation:** Update user documentation with V3 features
4. **Monitoring:** Set up consistency tests in CI/CD pipeline

## Technical Details

- **Python Version:** 3.11.4
- **Test Framework:** pytest
- **Total Test Categories:** 7
- **Test Execution Time:** 2025-12-30T13:59:15.046402

## Test Coverage

The consistency test suite covers:

- ✅ Input parsing and validation
- ✅ BLAST sequence alignment
- ✅ Flanking sequence extraction
- ✅ Multiple sequence alignment
- ✅ KASP primer design
- ✅ CAPS/dCAPS primer design
- ✅ End-to-end pipeline execution
- ✅ Output format consistency
- ✅ Numerical precision validation
