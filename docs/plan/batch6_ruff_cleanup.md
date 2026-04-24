# Batch 6: ruff 报错清零

## 状态

当前 101 个报错，分布如下：

| 文件 | 数量 | 主要类型 |
|------|------|---------|
| `core/alignment.py` | 24 | TRY003(15), PLR2004(4), PLR0912(2), PLR0915(1), SIM105(2) |
| `core/parser.py` | 16 | TRY003(13), RUF012(1), PLC0415(1), PLR2004(0) |
| `config.py` | 12 | TRY003(9), PLR2004(2) |
| `primers/kasp.py` | 11 | PLR0913(3), TRY003(2), RUF012(1), ARG002(2), PLR0915(1) |
| `primers/caps.py` | 11 | PLR0913(4), TRY003(2), RUF012(2), PLR0912(1) |
| `core/blast.py` | 10 | TRY003(5), PLR2004(2), PLR0913(1), PLR0912(1) |
| `main.py` | 8 | PLR0913(3), PLR0912(1), PLR0915(1), TRY400(2) |
| `core/primer3_parser.py` | 6 | TRY003(6) |
| `core/specificity.py` | 3 | TRY003(2), PLR0913(1) |

## 处理策略

### 批量忽略（加到 `pyproject.toml` 的 ignore 列表）

- **TRY003**（54 个）：项目使用自定义异常类 + f-string 消息是有意设计，直接 ignore。

### 逐文件忽略（加到 per-file-ignores）

- **PLR0913**（Too many arguments，12 个）：核心设计函数参数较多是领域需要。在 `pyproject.toml` 中对 `src/` 加 per-file-ignores。
- **PLR0912/PLR0915**（Too many branches/statements，6 个）：这些属于 #3（拆分大函数），留给后续重构解决，先 ignore。

### 逐个修复

- **PLR2004**（Magic values，8 个）：提取为常量。
- **RUF012**（Mutable class attrs，4 个）：加 `ClassVar` 注解。
- **TRY400**（2 个）：`logger.error` → `logger.exception`。
- **SIM105**（2 个）：`try/except OSError: pass` → `contextlib.suppress(OSError)`。
- **PLC0415**（2 个）：延迟 import，已有注释说明原因，加 `# noqa: PLC0415`。
- **ARG002**（2 个）：未使用参数，加 `# noqa: ARG002` 或实际移除。

## 执行顺序

1. 在 `pyproject.toml` 的 `[tool.ruff.lint]` ignore 列表加入 `"TRY003"`。
2. 在 per-file-ignores 中对 `src/**/*.py` 加 `"PLR0913"`, `"PLR0912"`, `"PLR0915"`。
3. 逐个修复剩余的 PLR2004、RUF012、TRY400、SIM105、PLC0415、ARG002。
4. 运行 `python -m ruff check src/` 确认 0 errors。
5. 运行 `python -m pytest -q` 确认全量测试通过。

## 验证命令

```bash
python -m ruff check src/
python -m pytest -q --no-header
```
