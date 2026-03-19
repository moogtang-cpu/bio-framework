# Bio-Framework v1.0.0 安装部署与使用说明

**产品名称**：Bio-Framework — AI 驱动的生物信息学分析框架
**版本**：v1.0.0
**类型**：Claude Code 技能插件
**作者**：Hanqing Tang \<moogtang@gmail.com\>
**许可证**：MIT
**发布日期**：2026-02-25

---

## 目录

1. [产品概述](#一产品概述)
2. [数据安全与隐私保障](#二数据安全与隐私保障)
3. [技术指标](#三技术指标)
4. [系统要求](#四系统要求)
5. [安装步骤](#五安装步骤)
6. [快速开始](#六快速开始)
7. [工作流概览](#七工作流概览)
8. [常见问题](#八常见问题faq)
9. [获取帮助](#九获取帮助)

---

## 一、产品概述

Bio-Framework 是一个基于 Claude Code 的生物信息学分析技能插件，提供从原始数据到 SCI 论文初稿的端到端全流程分析能力。

框架的核心设计理念是：**AI 负责分析决策与代码生成，本地环境负责数据处理与执行**。用户只需描述科学问题，框架即可自动规划分析流程、生成可复现的代码、执行质量控制，并最终生成符合 SCI 发表标准的图表与论文初稿。

### 主要特性

- **端到端工作流**：Step 0 初始化 → 动态多阶段分析 → Step 6 SCI 稿件，全流程自动化
- **多组学支持**：涵盖 scRNA-seq、空间转录组、Bulk RNA-seq、蛋白质组、代谢组、脂质组、ATAC-seq 等 7 种以上组学类型
- **知识库驱动**：内置 350KB+ 专业知识库，包含错误模式库、标准图形定义、SCI 设计标准
- **可调教的 AI 助手**：三层配置系统记忆你的参数偏好，经验库自动沉淀问题解决方案，越用越了解你的分析习惯
- **三语言支持**：支持英文、中文、日文命令交互
- **三层架构**：Project（项目级）> User Global（用户级）> Plugin（插件级），灵活的参数覆盖机制
- **自适应思维深度**：四档思维模式（Quick / Standard / Deep / Ultrathink），按分析复杂度自动切换

---

## 二、数据安全与隐私保障

**这是选择 Bio-Framework 的重要理由之一。**

### 完全本地运行

所有分析计算在用户的本地服务器或个人电脑上执行。数据文件（基因组数据、测序结果、蛋白质组数据等）始终保留在本地文件系统，**永不离开本地环境**。

### Claude AI 不接触用户实验数据

Bio-Framework 是 Claude Code 的技能插件，其工作方式如下：

| 角色 | 处理内容 |
|------|----------|
| Claude AI | 分析策略设计、R/Python 代码生成、结果解读与质量判断 |
| 本地 R/Python 环境 | 实际读取数据文件、执行计算、输出分析结果 |

Claude AI 只接收代码执行后的**摘要输出**（如统计指标、图表描述），不读取、不传输、不存储任何原始实验数据文件。

> **说明**：Claude Code 本身作为 AI 服务，会将代码文本和输出摘要发送至 Anthropic API 进行模型处理。但原始数据矩阵（如表达矩阵、测序文件）不在传输范围内。有关 Anthropic 的数据处理政策，请参阅 Anthropic 官方文档。

### 数据流向

```
用户数据文件 (fastq/h5ad/csv/...)
        |
        v
本地 R/Python 环境（Seurat/Scanpy/...）
        |
        v
本地输出文件（图表/报告/稿件草稿）
        |
        v (仅摘要信息)
Claude AI（分析指导与代码生成）
```

### 框架本身无网络传输

框架本身不含任何网络传输模块或云端存储功能。所有数据处理均通过本地代码完成，框架提供的是**分析指导能力**，而非数据处理服务。

### 适合敏感数据场景

- 临床患者样本数据
- 医院科研数据
- 企业内部保密数据
- 需要通过 IRB 审批的研究数据

Bio-Framework 的架构设计完全符合上述场景对数据隐私保护的要求。

---

## 三、技术指标

| 指标 | 数值 |
|------|------|
| 核心技能代码 | 10,318+ 行 |
| 知识库规模 | 350KB+（含 20KB 错误模式库、68KB 标准图形定义、36KB SCI 设计标准）|
| 顶级命令数 | 25+ |
| 子命令数 | 50+ |
| 多语言命令别名 | 110+（英文 / 中文 / 日文）|
| 支持组学类型 | 6 种标准流水线（scRNA-seq、空间转录组、Bulk RNA-seq、蛋白质组、代谢组、脂质组）+ ATAC-seq 等组学的知识库覆盖 |
| 期刊适配 | 6 种（Nature、Cell Reports、JCI、Nature Methods、Science、PLOS ONE）|
| 语言支持 | 3 种（英文、中文、日文）|
| 测试用例 | 898 个（自动化测试用例），通过率 100% |
| 累计审计会话 | 76+ 会话深度审计，400+ Bug 修复 |

---

## 四、系统要求

### 必要依赖

| 组件 | 最低版本 | 说明 |
|------|----------|------|
| Claude Code CLI | >= 2024.01 | 框架运行的基础环境，必须安装 |
| 操作系统 | Linux / macOS / WSL2 | Windows 用户请使用 WSL2 |

### 分析依赖（自动安装）

Bio-Framework 在分析过程中会根据课题需求**自动检测并安装**所需的全部依赖，包括 R、Python 基础环境以及具体的分析包（如 Seurat、DESeq2、scanpy 等）。你无需提前手动安装任何软件。

框架会在 Step 0 初始化阶段根据组学类型自动判断需要的依赖，并在首次使用时生成安装代码供执行。

### 可选组件

**claude-router 插件**：智能模型路由插件，可降低 API 使用成本（节省幅度取决于具体工作负载）。安装脚本运行时会自动询问是否安装，推荐安装。

---

## 五、安装步骤

### 步骤 1：确认 Claude Code CLI 已安装

在终端运行以下命令，确认 Claude Code CLI 已正确安装：

```bash
claude --version
```

如果提示命令不存在，请先参考 [Claude Code 官方文档](https://docs.anthropic.com/claude-code) 完成 CLI 安装。

### 步骤 2：获取 Bio-Framework

购买后，你将获得 `bio-framework-v1.0.0.tar.gz` 压缩包文件。

```bash
# 将压缩包复制到你的工作目录
cp ~/Downloads/bio-framework-v1.0.0.tar.gz ./

# 解压
tar -xzf bio-framework-v1.0.0.tar.gz
cd bio-framework-v1.0.0
```

### 步骤 3：运行安装脚本

```bash
# 添加执行权限
chmod +x install.sh

# 全局安装（推荐）
./install.sh --global
```

安装脚本将自动完成以下操作：

1. 将框架文件复制到 `~/.claude/bio-framework/`
2. 将技能文件安装到 `~/.claude/skills/bio/`
3. 询问是否安装 claude-router 插件
4. 验证安装完整性

**安装脚本参数说明**

| 参数 | 说明 |
|------|------|
| `--global` / `-g` | 全局安装（推荐，安装到 `~/.claude/bio-framework/`）|
| `--project` / `-p` | 项目级安装（仅当前目录生效）|
| `--upgrade` / `-u` | 升级安装（自动清理旧版文件后重新安装）|
| `--skip-router` | 跳过 claude-router 安装（降级模式，无成本优化）|

安装过程示例输出：

```
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Installing Bio-Framework v1.0.0
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

[INFO] Installing to global layer (~/.claude/bio-framework)
[INFO] Copying framework files...
[INFO] Installing skill files...
[INFO] Verifying installation...
[INFO] Installation complete!

Install claude-router now? [Y/n]:
```

### 步骤 4：初始化项目

切换到你的分析项目目录，运行初始化命令：

```bash
cd /path/to/your/analysis/project
bio-init
```

此命令会在当前目录创建 `.claude/` 配置目录和 `TOPIC.yml` 课题文件模板。

### 步骤 5：验证安装

在 Claude Code 中运行：

```
/bio version
```

如果输出类似以下内容，说明安装成功：

```
Bio-Framework v1.0.0
Framework: ~/.claude/bio-framework/
Status: Active
```

---

## 六、快速开始

Bio-Framework 支持两种启动方式，适应不同的准备阶段。

---

### 方式 A：从 TOPIC.yml 开始（标准流程）

适用于：你有明确的科学问题，希望框架自动规划完整分析方案。

**第一步：编写 TOPIC.yml 课题文件**

`TOPIC.yml` 是描述你的科学问题的配置文件，Bio-Framework 基于此文件自动规划分析流程。

在项目目录创建 `TOPIC.yml`：

```yaml
metadata:
  project_name: "肝癌单细胞转录组分析"
  created_date: "2026-03-10"
  author: "张三"
  research_type: "basic_research"
  quality_requirements: "standard"

scientific_questions:
  - Q1: "肝癌肿瘤微环境中存在哪些主要细胞类型？"
  - Q2: "不同细胞类型之间的细胞通讯模式如何？"
  - Q3: "肿瘤相关巨噬细胞的亚群特征是什么？"

hypotheses:
  - H1: "肝癌微环境中存在免疫耗竭的 T 细胞亚群"
  - H2: "M2 型巨噬细胞极化与肿瘤预后相关"

biological_context:
  organism: "Homo sapiens"
  tissue: "肝脏肿瘤组织"
  condition: "肝细胞癌 vs 癌旁正常组织"
  key_cell_types: ["T cells", "Macrophages", "NK cells", "Hepatocytes"]
```

TOPIC.yml 完整模板位于 `framework/templates/TOPIC_template.yml`，包含所有可配置字段的说明。

**第二步：启动分析**

打开 Claude Code，在项目目录中运行：

```
/bio start
```

或使用中文命令：

```
/生信 开始
```

框架将自动读取 TOPIC.yml，进行深度理解（ultrathink 模式），并规划适合你课题的完整分析方案。

---

### 方式 B：从预备材料启动（带数据描述和分析计划）

适用于：你已经准备好了数据集描述文件和/或分析方案，希望 AI 基于这些材料直接开始分析。

这种方式常见于以下场景：
- 已经从 GEO 等数据库下载好数据，整理了数据描述文件
- 已有明确的生信分析方案（Phase 设计、科学问题规划）
- 希望 AI 直接理解你的分析设计并执行

**第一步：准备材料文件**

在项目目录下准备以下文件（命名不限）：

**数据描述文件**（例如 `data/data_info.txt`）— 描述你使用的数据集：

```
GSE165816：

Organism    Homo sapiens
Experiment type    Expression profiling by high throughput sequencing
Summary    Diabetic foot ulcers (DFUs) are a devastating complication
           of diabetes. In order to identify systemic and local factors
           associated with DFU healing, we examined the cellular landscape
           of DFUs by single-cell RNA-seq analysis...

Overall design    Single-cell RNA sequencing of foot and forearm skin cells...

分组信息：
GSM5050521    G1: Forearm skin of Healthy non-diabetic subject
GSM5050522    G1A: Foot skin of Diabetic subject without DFU
GSM5050523    G2: Foot skin of subject with healing DFU
...

文件信息：
GSE165816_RAW.tar
```

**分析计划文件**（例如 `docs/bio_analyze.md`）— 描述你的分析设计：

```markdown
# HIF-1α与NF-κB缺氧-炎症耦联课题

## 生信可直接回答的科学问题

Q1: 糖尿病创面中哪些细胞类型是HIF1A/NFKB1共表达的主要载体？
Q2: HIF1A/NFKB1的共表达模式与创面愈合结局有何关联？
Q3: HIF-1α与NF-κB共同调控的靶基因有哪些？
...

## 生信分析Phase设计

Phase 1: 单细胞核心分析（GSE165816）
  - 1.1 数据预处理与质控
  - 1.2 细胞类型注释
  - 1.3 HIF1A/NFKB1表达分析
  ...

Phase 2: Bulk RNA-seq多队列验证
Phase 3: 时间动力学分析
Phase 4: 巨噬细胞HIF/NF-κB机制深度分析
```

**第二步：使用 `@` 引用材料启动**

打开 Claude Code，通过 `@` 符号引用你的材料文件：

```
请帮我基于以下材料开始生信分析：
@data/data_info.txt
@docs/bio_analyze.md
```

或者直接运行：

```
/bio start
```

然后在对话中用 `@` 引用文件提供上下文。AI 会自动读取你的数据描述和分析计划，理解课题设计，并生成对应的 TOPIC.yml 和分析方案。

> **提示**：方式 B 的核心优势在于你可以用自然语言和自由格式描述课题，无需严格遵循 TOPIC.yml 的 YAML 语法。AI 会从你的材料中提取科学问题、数据信息和分析需求，自动转化为结构化的工作流。

---

### 两种方式的对比

| 对比维度 | 方式 A：TOPIC.yml | 方式 B：预备材料 |
|---------|-----------------|----------------|
| 适合人群 | 熟悉 YAML 格式的用户 | 所有用户 |
| 准备工作 | 编写结构化 TOPIC.yml | 准备自由格式的数据描述和分析计划 |
| AI 理解深度 | 直接解析结构化字段 | 通过 AI 理解自然语言描述 |
| 灵活性 | 字段固定，覆盖全面 | 自由度高，可包含任意信息 |
| 适用场景 | 新建课题、标准化流程 | 已有详细分析设计、多数据集课题 |

---

### 启动后的自动执行流程

无论采用哪种方式启动，Bio-Framework 均采用结构化的自动执行策略：

- **自动执行**：代码生成、脚本运行、文件读写、质量指标计算、R/Python 包安装
- **暂停等待用户决策**：分析方案确认、生物学解读判断、图表视觉审查
- **遇到异常时**：自动查询知识库尝试修复，无法自动处理时暂停并报告

你可以随时运行 `/bio status` 查看当前进度。

### 常用命令速查

| 命令 | 中文别名 | 功能 |
|------|----------|------|
| `/bio start` | `/生信 开始` | 启动新的分析工作流 |
| `/bio continue` | `/生信 继续` | 从断点恢复分析 |
| `/bio stop` | `/生信 停止` | 暂停并保存当前进度 |
| `/bio status` | `/生信 状态` | 查看分析进度 |
| `/bio prime` | `/生信 预热` | 生成会话上下文快照（新会话开始时使用）|
| `/bio help` | `/生信 帮助` | 显示完整帮助信息 |
| `/bio phase list` | - | 列出所有分析阶段 |
| `/bio finding add` | - | 记录分析发现或异常 |
| `/bio pipeline show` | - | 查看流水线合规状态 |
| `/bio manuscript check` | - | 检查稿件质量 |
| `/bio config show` | - | 查看当前参数配置 |
| `/bio version` | - | 显示框架版本信息 |

---

## 七、工作流概览

Bio-Framework 采用结构化的七步工作流，从课题理解到 SCI 稿件一气呵成。

```
Step 0: 初始化
  |-- 读取 TOPIC.yml，深度理解科学问题
  |-- 探测数据文件，识别组学类型
  |-- 规划 Phase 1-N 分析阶段
  |-- 确认流水线合规参数
  v
Phase 1-N: 动态分析阶段
  |-- Phase 1: 数据质量控制与预处理
  |-- Phase 2: 降维聚类与细胞类型注释
  |-- Phase N: 根据科学问题动态扩展（细胞通讯 / 轨迹分析 / 差异表达...）
  v
Step 1.5: 质量审计
  |-- 三层质量控制验证
  |-- 数据质量综合评估
  |-- 异常记录与修复建议
  v
Step 2: 反思与探索
  |-- 回顾分析结果，深度反思生物学意义
  |-- 发现潜在规律，决定是否补充分析
  v
Step 3: 研究总结与图表规划
  |-- 提炼核心发现（ultrathink）
  |-- 规划发表图表（figure plan）
  |-- 验证必需图表完整性
  v
Step 4: 图表生成
  |-- 按 SCI 标准生成发表级图表
  |-- 视觉审查与用户确认
  v
Step 5: 技术报告
  |-- 生成技术基础文档
  |-- 方法学描述、参数记录、可复现性验证
  v
Step 6: SCI 稿件草稿
  |-- 基于全流程结果生成论文初稿（ultrathink）
  |-- 符合目标期刊格式要求
  |-- 图表引用完整，参考文献规范
```

### 典型执行时间参考

| 工作流阶段 | 预计时间（参考）|
|------------|-----------------|
| Step 0 初始化 + Phase 规划 | 10-20 分钟 |
| Phase 1-N 动态分析 | 数小时至数天（取决于数据量和分析复杂度）|
| Step 1.5-3 质量审计与总结 | 30-60 分钟 |
| Step 4 图表生成 | 1-3 小时 |
| Step 5-6 报告与稿件 | 1-2 小时 |

---

## 八、常见问题（FAQ）

### Q1：安装失败，提示"Permission denied"

**原因**：安装脚本没有执行权限，或目标目录权限不足。

**解决方法**：

```bash
# 添加执行权限
chmod +x install.sh

# 确认 ~/.claude/ 目录存在且有写权限
ls -la ~/.claude/
```

如果 `~/.claude/` 目录不存在，Claude Code CLI 可能未正确安装，请重新检查 Claude Code 的安装步骤。

---

### Q2：运行 `/bio start` 后提示找不到 TOPIC.yml

**原因**：TOPIC.yml 未创建，或 Claude Code 未在正确的项目目录下运行。

**解决方法**：

```bash
# 确认在项目目录下
pwd

# 确认 TOPIC.yml 存在
ls TOPIC.yml

# 如果不存在，运行 bio-init 创建模板
bio-init
```

---

### Q3：分析中途中断，如何恢复？

Bio-Framework 在每个分析节点自动保存进度到 `project_skills/runtime/workflow_state.yml`。

**恢复方法**：

```bash
# 在新会话中运行
/bio prime    # 恢复上下文快照
/bio continue # 从断点继续
```

系统会自动检测上次中断位置并从该点恢复，无需重复已完成的分析步骤。

---

### Q4：如何更新到新版本？

**方法一：使用 --upgrade 参数（推荐）**

```bash
cd /path/to/bio-framework-new-version
./install.sh --upgrade
```

`--upgrade` 模式会自动清理旧版框架文件，再安装新版本。用户数据（`rules/`、`user_knowledge/` 等）不受影响。

**方法二：手动升级**

```bash
# 备份用户数据（可选，通常不受影响）
cp -r ~/.claude/bio-framework/rules ~/.claude/bio-framework/rules.backup

# 运行安装脚本
./install.sh --global
```

---

### Q5：如何卸载 Bio-Framework？

```bash
# 移除框架文件
rm -rf ~/.claude/bio-framework
rm -rf ~/.claude/skills/bio

# 如果安装了 claude-router 且不再需要，同样可以移除
rm -rf ~/.claude/plugins/claude-router
```

卸载后，项目目录中的 `.claude/` 文件夹（包含分析配置和结果）不受影响，你的分析数据和产出物均保留在本地。

---

### Q6：是否支持 HPC 集群环境？

支持。Bio-Framework 对运行环境没有特殊要求，只要以下条件满足即可：

- Claude Code CLI 已安装
- R 和/或 Python 环境可用
- 用户对项目目录有读写权限

在集群环境下，通常在登录节点运行 Claude Code（负责代码生成和分析决策），实际的 R/Python 脚本可以提交到计算节点作为作业运行。

---

### Q7：如何配置针对特定项目的参数偏好？

在项目目录的 `.claude/rules/project_rules.yml` 中定义项目级参数：

```yaml
# .claude/rules/project_rules.yml
clustering:
  resolution:
    value: 1.2
    reason: "本课题数据异质性高，需要更细的聚类分辨率"

quality_control:
  mt_percent:
    value: 10
    reason: "心肌细胞样本，线粒体比例天然较高，使用宽松阈值"
```

参数优先级：项目级 > 用户全局 > 框架默认值。

---

## 九、获取帮助

### 框架内置帮助

```
/bio help                   # 显示完整命令帮助
/bio help start             # 查看特定命令的详细说明
/bio debug why [行为描述]   # 调试 AI 决策原因
```

### GitHub Issues

在以下地址提交 Bug 报告或功能建议：

```
https://github.com/moogtang-cpu/bio-framework/issues
```

提交 Issue 时，请附上：
- 操作系统和 Claude Code 版本
- 复现步骤
- 错误信息（如有）
- `TOPIC.yml` 相关内容（请移除敏感信息）

### 联系作者

- **作者**：Hanqing Tang
- **邮件**：moogtang@gmail.com
- **项目主页**：https://bio-framework.github.io

---

*Bio-Framework v1.0.0 — AI 驱动的生物信息学分析框架*
*从原始数据到 SCI 论文初稿的端到端解决方案*
