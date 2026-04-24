# Batch 7: 清理死代码

## 目标

删除已确认不再使用的代码，减少维护负担。

## 待清理项

### 1. `KASPDesigner._format_primer_line`（legacy 方法）

位置：`src/snp_primer_pipeline/primers/kasp.py:883-925`

状态：
- 已被 `_format_primer_pair_v2` + `_convert_to_v2_format` 替代。
- 无调用方引用。

验证方式：
```bash
grep -rn "_format_primer_line" src/snp_primer_pipeline/primers/kasp.py
# 应只有定义处，无调用
grep -rn "_format_primer_line" src/ --include="*.py" | grep -v "def _format_primer_line"
```

操作：删除整个方法。

### 2. `CAPSDesigner._format_primer_line`

位置：`src/snp_primer_pipeline/primers/caps.py:480-510`

状态：
- 被 `format_output` 调用，**不是死代码**。
- 保留。

### 3. 验证 `PipelineConfig.from_args` 是否仍存在

```bash
grep -rn "from_args" src/
```

状态：已在前期清理中删除，确认不再存在。

## 执行顺序

1. 确认 `KASPDesigner._format_primer_line` 无调用方。
2. 删除该方法。
3. 运行 `python -m pytest -q` 确认无回归。

## 验证命令

```bash
python -m pytest -q --no-header
python -m ruff check src/
```
