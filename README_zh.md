# SNP Primer Pipeline 3

[English Documentation](README.md) | **中文文档**

这是一个现代化的、模块化的 Python 流程，用于设计任何具有参考基因组物种的 SNP 基因分型所需的 KASP 和 CAPS/dCAPS 引物。这是 SNP_Primer_Pipeline2 的重构版本，具有改进的架构、可测试性和可维护性。

## 功能特性

- **KASP 引物设计**：设计用于 KASP 基因分型的等位基因特异性引物
- **CAPS/dCAPS 引物设计**：设计用于基于限制性内切酶基因分型的引物
- **多物种支持**：适用于任何具有参考基因组的物种
- **模块化架构**：关注点分离清晰，组件可测试
- **现代 Python**：使用 Python 3.10+，包含类型提示和 dataclasses
- **全面测试**：使用 hypothesis 进行基于属性的测试
- **CLI 接口**：易于使用的命令行界面

## 安装指南

### 前置要求

- Python 3.10 或更高版本
- BLAST+ (blastn, blastdbcmd)

### 源码安装

```bash
git clone <repository-url>
cd SNP_Primer_Pipeline3
pip install -e .
```

## 快速开始

### 基本用法

```bash
# 同时设计 KASP 和 CAPS 引物
snp-primer input.csv reference_db -o output_dir

# 仅设计 KASP 引物
snp-primer input.csv reference_db --no-caps

# 仅设计 CAPS 引物，并设置自定义酶价格限制
snp-primer input.csv reference_db --no-kasp --max-price 150
```

### 输入格式

输入文件应为 CSV 格式，包含以下内容：

```
SNP_name,chromosome,flanking_sequence
SNP1,chr1A,ATCGATCGATCG[A/G]TCGATCGATCG
SNP2,chr2B,GCTAGCTAGCTA[C/T]AGCTAGCTAGT
```

其中：
- `SNP_name`：SNP 的唯一标识符
- `chromosome`：目标染色体名称
- `flanking_sequence`：包含 IUPAC 括号标记的 SNP 的 DNA 序列

### Python API 用法

```python
from snp_primer_pipeline import (
    PipelineConfig, 
    run_pipeline,
    PolymarkerParser,
    KASPDesigner,
    CAPSDesigner
)

# 创建配置
config = PipelineConfig(
    input_file="input.csv",
    reference_file="reference_db",
    output_dir="output",
    design_kasp=True,
    design_caps=True
)

# 运行完整流程
run_pipeline(config)

# 或者使用单独的组件
parser = PolymarkerParser("input.csv")
snps = parser.parse()

kasp_designer = KASPDesigner()
primers = kasp_designer.design_primers(
    template_sequence="ATCGATCG...",
    snp_position=100,
    snp_alleles=("A", "G")
)
```

## 架构设计

本流程分为几个模块：

### 核心模块 (`snp_primer_pipeline.core`)

- **`parser.py`**：解析 polymarker 输入文件并转换 IUPAC 代码
- **`blast.py`**：BLAST 执行、结果解析和侧翼区域提取
- **`alignment.py`**：多序列比对和变异位点识别
- **`primer3_parser.py`**：用于引物设计的 Primer3 接口

### 引物设计模块 (`snp_primer_pipeline.primers`)

- **`kasp.py`**：设计带有等位基因特异性引物的 KASP 引物
- **`caps.py`**：通过限制性内切酶分析设计 CAPS/dCAPS 引物

### 数据模型 (`snp_primer_pipeline.models`)

- **`SNP`**：包含等位基因和位置的 SNP 信息
- **`BlastHit`**：BLAST 比对结果
- **`FlankingRegion`**：SNP 周围的基因组区域
- **`Primer`**：单个引物属性
- **`PrimerPair`**：带有评分的引物对
- **`RestrictionEnzyme`**：用于 CAPS 设计的酶属性

### 配置 (`snp_primer_pipeline.config`)

- **`PipelineConfig`**：主流水线配置
- **`SoftwarePaths`**：外部软件路径管理

## 输出文件

流程生成以下输出：

### KASP 引物
- `KASP_primers_{SNP_name}.txt`：设计的 KASP 引物及其评分
- 每个 SNP 等位基因的特异性引物
- 用于扩增的通用引物

### CAPS 引物
- `CAPS_primers_{SNP_name}.txt`：设计的 CAPS/dCAPS 引物
- 酶信息和切割位点
- 用于限制性内切酶分析的引物对

### 中间文件
- `for_blast.fa`：用于 BLAST 搜索的 FASTA 文件
- `blast_out.txt`：BLAST 结果
- `flanking_sequences.fa`：提取的侧翼序列
- `alignment_{SNP_name}.fa`：多序列比对文件

## 配置选项

### 命令行选项

```
--max-tm FLOAT          最大引物 Tm 值 (默认: 63.0)
--max-size INT          最大引物长度 (默认: 25)
--max-price INT         CAPS 酶的最大价格 (默认: 200)
--pick-anyway           即使违反约束条件也强制挑选引物
--threads INT           BLAST 线程数 (默认: 1)
--log-level LEVEL       日志级别 (DEBUG/INFO/WARNING/ERROR)
```

### 配置文件

创建 `config.yaml` 文件以进行高级配置：

```yaml
# 流程设置
flanking_size: 500
max_hits: 6
primer_product_size_range: [50, 250]

# 引物设计设置
max_tm: 63.0
max_primer_size: 25
pick_anyway: false

# CAPS 设置
max_price: 200

# 软件路径 (如未指定则自动检测)
primer3_path: "/usr/local/bin/primer3_core"
muscle_path: "/usr/local/bin/muscle"
```

## 测试指南

本流程包含全面的测试：

```bash
# 运行所有测试
pytest tests/

# 带覆盖率报告运行
pytest tests/ --cov=snp_primer_pipeline --cov-report=html

# 运行基于属性的测试
pytest tests/ -k "property"
```

## 开发指南

### 项目结构

```
SNP_Primer_Pipeline3/
├── src/snp_primer_pipeline/     # 主包
│   ├── core/                    # 核心处理模块
│   ├── primers/                 # 引物设计模块
│   ├── config.py               # 配置管理
│   ├── models.py               # 数据模型
│   ├── exceptions.py           # 自定义异常
│   └── main.py                 # CLI 接口
├── tests/                      # 测试套件
├── resources/                  # 资源文件
│   └── NEB_parsed_REs.txt     # 酶数据库
├── pyproject.toml             # 项目配置
└── README.md                  # 说明文件 (英文)
```

### 添加新功能

1. **新引物类型**：扩展 `primers` 模块
2. **新输入格式**：扩展 `parser` 模块
3. **新比对工具**：扩展 `alignment` 模块
4. **新数据模型**：添加到 `models.py`

### 代码质量

项目遵循现代 Python 最佳实践：

- 全程使用类型提示
- 使用 Dataclasses 定义数据模型
- 使用 hypothesis 进行基于属性的测试
- 全面的错误处理
- 模块化、可测试的架构

## 与 SNP_Primer_Pipeline2 的对比

### 改进点

1. **现代 Python**：使用 Python 3.10+ 特性
2. **模块化设计**：关注点分离清晰
3. **类型安全**：完整的类型提示，更好的 IDE 支持
4. **可测试性**：全面的测试套件，包含基于属性的测试
5. **错误处理**：使用自定义异常进行稳健的错误处理
6. **文档**：全面的文档和示例
7. **可维护性**：遵循最佳实践的清晰代码结构

### 迁移指南

对于 SNP_Primer_Pipeline2 的用户：

1. **输入格式**：相同的 polymarker CSV 格式
2. **输出格式**：相似，但组织更好
3. **依赖**：需要 Python 3.10+，不再支持 Python 2/3
4. **安装**：使用现代 pip/pyproject.toml，无需手动设置
5. **配置**：使用 YAML 配置代替硬编码值

## 许可证

本项目采用 GNU General Public License v2.0 - 详情见原版 SNP_Primer_Pipeline2。

## 贡献

欢迎贡献！请：

1. Fork 本仓库
2. 创建特性分支
3. 为新功能添加测试
4. 确保所有测试通过
5. 提交 Pull Request

## 支持

如有问题和支持需求：

1. 查看文档
2. 查看现有 issue
3. 创建包含详细信息的新 issue
4. 报告 bug 时请包含输入文件和错误信息
