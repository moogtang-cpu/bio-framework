# Bio-Framework

> **不是又一个 AI 聊天机器人。是 350 KB+ 结构化知识驱动的分析导航系统。**

---

## 先用人话说

你有一批组学数据（单细胞、空间、蛋白质组……），想发 SCI 论文。

Bio-Framework 让 AI 按照生物信息学领域的最佳实践，帮你完成**从数据质控到出版级图表再到论文初稿**的全部工作。在关键节点（数据质量、统计正确性、期刊格式）设置了硬性门控，防止 AI 犯错或编造数据。

**你不需要自己写 R/Python 代码**——AI 会自动生成、执行、调试。但你需要理解自己的科学问题，并在关键判断点（细胞类型注释、分析方向选择）做出决策。

所有分析**在你的本地环境运行**。原始数据不会上传到任何第三方服务器。

<div align="center">
<img src="docs/images/screenshot_bio_status.png" width="600"/>
<br/>
<em><code>/bio status</code> — 真实项目生命周期：86,139 细胞、4 个 Phase，从 GSE165816 质控到论文初稿，全程自动管理。</em>
</div>

---

## 引擎生成的图表展示

以下每一张图都是 **Bio-Framework 自动生成**的——从原始数据到出版级输出，零手动绘图代码。

<table>
<tr>
<td align="center"><b>UMAP 细胞类型着色</b><br>(scRNA-seq，DFU 研究)</td>
<td align="center"><b>QC 小提琴图</b><br>(质控指标)</td>
</tr>
<tr>
<td><img src="docs/images/showcase_umap_celltype.png" width="400"/></td>
<td><img src="docs/images/showcase_qc_violin.png" width="400"/></td>
</tr>
<tr>
<td align="center"><b>出版级组合图</b><br>(多面板组装 + 统一样式)</td>
<td align="center"><b>空间转录组 + 轨迹分析</b><br>(子宫内膜干细胞图谱)</td>
</tr>
<tr>
<td><img src="docs/images/showcase_figure1_dfu.png" width="400"/></td>
<td><img src="docs/images/showcase_spatial_trajectory.png" width="400"/></td>
</tr>
<tr>
<td align="center"><b>图谱总览</b><br>(314,805 细胞，5 数据集整合)</td>
<td align="center"><b>通路与耦合分析</b><br>(DFU 中 HIF-1α/NF-κB)</td>
</tr>
<tr>
<td><img src="docs/images/showcase_atlas_overview.png" width="400"/></td>
<td><img src="docs/images/showcase_figure3_dfu.png" width="400"/></td>
</tr>
</table>

> 所有图表 300 DPI、色盲友好配色、英文标签自动强制。完整输出见 [DFU 示例](examples/DFU-HIF-NFkB-scRNA-analysis/) 和 [子宫内膜图谱示例](examples/endometrial-stem-cell-atlas/)。

---

## 为什么不能直接用 Claude + Ultrathink？

**一句话回答**：Claude 是发动机，Bio-Framework 是整辆车——包括方向盘、导航、安全气囊和行车记录仪。

| 维度 | 单独用 Claude | + Bio-Framework |
|:---|:---|:---|
| **跨会话记忆** | 每次从零开始 | `/bio continue` 精确恢复到子步骤 |
| **质量门控** | 没有强制检查点 | 五维度审计 + 9 种强制 Ultrathink 场景 |
| **错误处理** | 每次从头推理 | 16 种已知错误模式 + 跨项目经验积累 |
| **SCI 图表** | 不验证格式 | SCI 标准 + AI 视觉检查 + 用户审查门控 |
| **论文** | 能写段落 | 21 步端到端稿件 + 交叉验证 |
| **可重复性** | 不自动记录参数 | 种子、版本、参数全部自动记录 |

> 完整对比：[为什么要使用 Bio-Framework](docs/why-bio-framework-zh.md)

---

## 我们解决的"坑"

Bio-Framework 的知识库编码了 **16 种错误模式**，来自真实生信分析的惨痛教训。这里是每个分析师都遇到过的三个：

### 坑 1：批次效应伪装成生物学差异

**场景**：你整合了两批 scRNA-seq 样本。聚类看起来很棒——15 个干净的 cluster。但其中 6 个实际上是批次效应伪影。你的 DEG 比较的是测序化学差异，不是细胞类型差异。

**框架如何拦截**：Step 1.5（强制质量审计）在任何生物学解读开始之前量化残留批次效应。批次校正不充分？**HARD_STOP**——工作流停下等你决策。你在分析阶段修复批次校正，不是在审稿人 #2 的 rebuttal 中才发现。

### 坑 2：过聚类 → 幻影细胞类型

**场景**：分辨率设太高。一个 T 细胞群体分裂成三个 marker 几乎相同的 cluster。你报告了"新型 T 细胞亚型"——审稿人要求验证，你什么都没有。

**框架如何拦截**：聚类分辨率是**强制 Ultrathink** 场景。AI 必须评估多个分辨率、检查过聚类信号（>80% marker 重叠的 cluster）、显式论证选定的分辨率。不允许"用 0.8 因为它是默认值"。

### 坑 3：静默数据造假

**场景**：GEO 下载失败。AI 为了"帮忙"，悄悄生成模拟表达数据继续工作流。你的图表现在包含捏造数据，而你永远不会知道。

**框架如何拦截**：错误 F003（`severity: critical`，HARD_STOP）。真实数据缺失或解析失败时，AI **必须停下来等待**。模拟数据生成被全面禁止。稿件中每个统计数字在 Step 6e 交叉验证中追溯到实际分析输出文件。

> 这不是假设场景。这些是浪费数月工作、导致撤稿的真实失败模式。

---

## 核心亮点

| 能力 | 一句话说明 |
|:---|:---|
| **全流程自动化** | 从 Step 0 初始化到 Step 6 稿件，7 个步骤自动推进，仅在需要科学判断时暂停 |
| **四级思考调度** | 创建文件夹用 Quick，细胞注释强制 Ultrathink——不是所有任务都值得同样的认真程度 |
| **五维度质量审计** | Step 1.5 硬性门控，发现样本混淆或严重批次效应直接停下等你决策 |
| **6 大期刊合规** | Nature/Science/Cell Reports 等出版标准自动检查，图表直接对齐 DPI 和配色要求 |
| **数据真实性保护** | 禁止生成模拟数据，稿件中每个统计数字追溯到源文件，五步提交前交叉验证 |
| **断点续传** | 分析中断后 `/bio continue` 精确恢复到子步骤级别，不丢失任何进度 |

---

## 架构一览

```
┌─────────────────────────────────────────────────────────┐
│                    你的 Claude Code                       │
│                                                          │
│  ┌──────────────────────────────────────────────────┐   │
│  │              Bio-Framework（Skill 插件）           │   │
│  │                                                    │   │
│  │  ┌────────────┐  ┌────────────┐  ┌────────────┐  │   │
│  │  │ 编排引擎    │  │  知识库     │  │  工作流     │  │   │
│  │  │ ORCHESTRATOR│  │ KNOWLEDGE  │  │ CONTROLLER  │  │   │
│  │  │  6 个子模块 │  │  350 KB+   │  │             │  │   │
│  │  │ • 自动执行  │  │ • 流水线   │  │ • /bio start│  │   │
│  │  │ • 思考调度  │  │ • 错误模式 │  │ • /bio cont │  │   │
│  │  │ • 断点恢复  │  │ • 期刊规则 │  │ • /bio stat │  │   │
│  │  │ • 步骤管理  │  │ • SCI 标准 │  │ • 60+ 命令  │  │   │
│  │  └────────────┘  └────────────┘  └────────────┘  │   │
│  └──────────────────────────────────────────────────┘   │
│                          │                               │
│              生成并执行代码                                │
│                          ▼                               │
│  ┌──────────────────────────────────────────────────┐   │
│  │         你的本地 R/Python 环境                      │   │
│  │         (Seurat, Scanpy, DESeq2 等)                │   │
│  └──────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────┘

三层配置优先级：
  优先级 1：Project .claude/rules/     ← 项目级覆盖（最高）
  优先级 2：User ~/.claude/rules/      ← 你的全局偏好
  优先级 3：Plugin framework/rules/    ← 框架默认值
```

> 完整架构文档：[Architecture Documentation](docs/architecture-zh.md)

---

## 让"湿实验"博士无法拒绝的 9 大卖点

### 1. 拒绝"AI 幻觉"：四级思考深度调度

**痛点：** AI 给建议很随性。创建文件夹和鉴定细胞类型用同一种逻辑深度。

**Bio-Framework 的方案：** 四级思考深度，根据任务关键性自动调度：

| 思考级别 | 触发场景 | AI 行为 |
|:---|:---|:---|
| **Quick** | 文件操作、格式调整 | 快速执行 |
| **Standard** | 常规分析步骤 | 标准推理 |
| **Deep** | 参数选择、方法决策 | 多方案比较 |
| **Ultrathink** | 9 种关键场景 | 完整推理链：观察→假设→验证→结论 |

**9 种强制 Ultrathink 场景**：细胞类型注释、数据质量审查（Step 1.5）、差异表达解读、生物学结论推导、稿件关键章节、数据验证（6e）、聚类分辨率选择、批次效应评估、异常结果判断。

---

### 2. 拒绝"带毒数据"：Step 1.5 五维度深度审计

在 Step 1（计算分析）完成后、Step 2（解读）开始前，框架**强制插入**质量审计。五个维度，一个都不能跳：

| 审计维度 | 检查内容 |
|:---|:---|
| **样本身份** | 样本标签是否正确、交叉污染信号 |
| **批次效应** | 批次校正是否充分、残留批次效应量化 |
| **技术质量** | 测序深度、检测基因数、线粒体比例 |
| **生物学合理性** | 已知 marker 是否在预期细胞群中 |
| **统计假设** | 样本量是否支撑统计结论、多重检验校正 |

结果写入 `quality_report.yml`。**只有全部通过，工作流才继续。** HARD_STOP 门控。

---

### 3. 拒绝"格式退稿"：内置 6 大顶刊合规引擎

框架内置 **6 种顶刊**的精确合规配置（机器可读 YAML）：

| 期刊 | 关键约束（示例） |
|:---|:---|
| **Nature** | 摘要 ≤150 词，正文 ≤3000 词，引文 ≤50 条 |
| **Science** | 独立稿件格式规范 |
| **Cell Reports** | 图表规范、STAR Methods 必填 |
| **JCI** | 临床研究附加要求、数据共享声明 |
| **Nature Methods** | 方法验证基准、可重复性声明 |
| **PLOS ONE** | 开放获取格式、Data Availability 必填 |

AI 自动检查图表分辨率（300 DPI）、字体大小、色盲友好配色、字数限制，并生成 `submission_readiness.md` 投稿前自查报告。

---

### 4. 拒绝"数据造假"：数据真实性硬性保护

**三道防线，零容忍：**

1. **禁止生成模拟数据（HARD_STOP）** — 错误 F003，`severity: critical`。数据缺失→停下等待，永不生成。
2. **稿件绝对禁止清单** — 不得捏造统计值、引用或结果。信息不足→`[Ref]` 占位符。
3. **五步提交前交叉验证（6e–6i）** — 每个统计数字追溯到源文件，每个结论有图表支撑，跨文件数字一致性验证。

---

### 5. 六条标准化组学流水线

| 组学类型 | 支持的技术平台 |
|:---|:---|
| **单细胞 RNA-seq** | 10x Genomics / Seurat / Scanpy |
| **空间转录组** | Visium / MERFISH / Slide-seq / STARmap / CODEX |
| **批量 RNA-seq** | DESeq2 / edgeR / limma |
| **蛋白质组** | MaxQuant / DIA-NN / TMT / iTRAQ |
| **代谢组** | XCMS / MetaboAnalyst / MZmine |
| **脂质组** | 继承代谢组 + 脂质特异性扩展 |

每条流水线定义 **mandatory steps**，实时合规检查。跳过 doublet removal？Warning。没做 batch correction？Warning。过程中拦截，不是事后补救。

**模块化使用**：已有分析结果？可直接从 Step 4（图表）或 Step 6（稿件）开始。

---

### 6. 自学习经验库：用得越多越聪明

```
遇到问题 → 查框架知识库（16种常见错误）
         → 查用户历史经验（learned_solutions.yml）
         → AI 自主推理解决
         → 解决过程 > 2 轮 → 自动保存到用户知识库
```

Seurat v5 API 变更、CellChat 内存溢出、R 包编译报错——全部自动记录。下个项目遇到同样问题，直接调用方案，`times_reused` 追踪复用。

---

### 7. 断点续传：跑到一半断了不用从头来

```bash
/bio continue
```

系统扫描已完成输出、与计划对比、**从第一个缺失项开始续做**。Step 6 有 **21 个独立子步骤**，每个写独立文件——只重试失败步骤。

---

### 8. SCI 图表知识库：审稿人看什么，它就检查什么

内置 732 行 SCI 设计标准（`sci_design_standards.yml`）：

- **图表决策树** — n<10 点图+配对线，n=10-30 箱线图+抖动点，n>30 小提琴图
- **错误图表防御** — 8 种常见错误及替代方案（条形图隐藏分布→箱线图）
- **AI 写作真实性检测** — 标记 7 种 AI 生成文字特征
- **审稿人必查清单** — 15+ 个高频质疑项，每项含修复方案

---

### 9. 三层可定制参数：不是黑箱

```
Project .claude/rules/     ← 项目级覆盖（最高优先级）
User ~/.claude/rules/      ← 你的全局偏好
Plugin framework/rules/    ← 框架默认值
```

12 个参数域全部可定制。你的偏好（"MT% 阈值 15%"、"聚类分辨率 0.6"）保存为全局默认，所有新项目自动继承。

---

## 完整工作流一览

```
Step 0  项目初始化 ─── 自动检测组学类型、创建标准目录结构
  │
Step 1  动态分析 ───── 多阶段计算分析（Phase 1-N），每阶段强制检查 mandatory steps
  │
Step 1.5 数据质量审计 ─ 五维度深度审查，HARD_STOP 门控
  │
Step 2  反思与探索 ─── 科学问题完成度打分（80%/70%），定向探索建议
  │
Step 3  研究总结 ───── 核心发现提炼、图表计划制定
  │
Step 4  出版级图表 ─── SCI 标准输出 + 视觉审查闭环
  │
Step 5  报告生成 ───── 结构化分析报告
  │
Step 6  稿件撰写 ───── 期刊合规 + 五步交叉验证（6e-6i）
```

---

## 真实项目示例

Bio-Framework 附带**两个完整示例项目**——真实研究工作流及全部输出，不是占位符。

### 示例 1：糖尿病足溃疡中的 HIF-1α/NF-κB（scRNA-seq）

- **规模**：50 个 scRNA-seq 样本 + 3 个 Bulk RNA-seq 队列 + 4 个芯片数据集
- **分析阶段**：4 个 Phase（单细胞→批量验证→时间动态→巨噬细胞深挖）
- **输出**：7 张主图 + 4 张附图 + 论文初稿 + 技术报告

```
examples/DFU-HIF-NFkB-scRNA-analysis/
├── project_skills/TOPIC.yml              # 研究问题定义
├── scripts/                              # 11 个 R 脚本（自动生成）
│   ├── phase1_1_qc.R → phase1_6_8_subpop.R
│   ├── phase2_bulk.R
│   ├── phase3_temporal.R
│   └── phase4_macrophage.R
└── Phase_output/
    ├── Phase1/                           # QC、聚类、表达、通路...
    ├── Phase2-4/                         # 批量验证、时间动态、巨噬细胞
    ├── publication_figures/              # Figure1-7 + FigureS1-S4（PDF + PNG）
    ├── manuscript/manuscript_draft.md    # 完整 SCI 论文初稿
    └── final_report/technical_report.md  # 结构化分析报告
```

> [浏览完整 DFU 示例 →](examples/DFU-HIF-NFkB-scRNA-analysis/)

### 示例 2：子宫内膜干细胞图谱（多数据集整合 + 空间转录组）

- **规模**：5 个 scRNA-seq 数据集 314,805 细胞 + 8 个 Visium 空间样本
- **分析阶段**：8 个 Phase（QC→整合→空间→疾病→轨迹→通讯→药物靶点）
- **输出**：6 张主图 + 论文初稿 + 药物靶点优先级列表

```
examples/endometrial-stem-cell-atlas/
├── project_skills/TOPIC.yml              # 研究问题定义
├── scripts/                              # Python 分析脚本
├── env/environment.yml                   # Conda 环境规范
└── Phase_output/
    ├── phase3_integration/               # 图谱构建、UMAP、注释
    ├── phase4_spatial/                   # Visium 反卷积、空间图
    ├── phase5_disease/                   # DEG 表、GSEA 结果（40+ 文件）
    ├── phase6_trajectory/                # DPT、PAGA、TF 活性
    ├── phase7_communication/             # LIANA 细胞间通讯
    ├── phase8_drug_targets/              # 463 候选靶点、优先级评分
    ├── publication_figures/              # Fig1-6（PDF + PNG）
    └── manuscript/manuscript_draft.md    # 完整 SCI 论文初稿
```

> [浏览完整子宫内膜图谱示例 →](examples/endometrial-stem-cell-atlas/)

---

## 数据安全与隐私

Bio-Framework 完全运行在 Claude Code 本地环境中：

- **数据不离机**：所有原始数据（h5、rds、csv 等）仅由你本地的 R/Python 环境读取和处理
- **本地代码执行**：AI 生成的分析代码在你的机器上运行，不在云端
- **API 通信范围**：Claude Code 与 Anthropic 服务器通信，发送的是代码逻辑和输出摘要——原始数据矩阵不在传输范围内

> **注意**：对于数据安全要求极高的场景（如受 IRB 约束的临床数据），用户应自行评估是否符合其机构的数据管理政策。详情请参阅 [Anthropic 的数据处理政策](https://www.anthropic.com/policies/privacy)。

---

## 使用前提与成本

| 前提 | 说明 |
|:---|:---|
| **Claude Code** | 需要 Anthropic 官方的 [Claude Code](https://docs.anthropic.com/en/docs/claude-code/overview) CLI 工具（需 Max 计划或 API 额度） |
| **本地环境** | R >= 4.3 和/或 Python >= 3.9，推荐通过 Conda 管理 |
| **操作系统** | macOS / Linux（WSL2 亦可） |
| **硬件建议** | 16GB+ 内存；大型单细胞数据（>50k cells）建议 64GB+ |
| **Bio-Framework** | 一次性购买，包含所有后续更新 |

**总成本 = Bio-Framework 一次性购买 + Claude Code 使用费用（按 Anthropic 定价）。**

---

## 快速上手

<div align="center">
<img src="docs/images/screenshot_bio_help.png" width="600"/>
<br/>
<em><code>/bio help</code> — 不是黑箱。工作流控制、主题管理、Phase 管理——一个完全可控的系统。</em>
</div>

<br/>

**第一步：安装框架**

```bash
./install.sh --global
```

**第二步：在你的数据目录初始化项目**

```bash
cd your_project && bio-init
```

填写生成的 `TOPIC.yml`——用自然语言描述你的研究问题、数据路径、目标期刊。

**第三步：开始分析**

```
/bio start --topic TOPIC.yml
```

框架自动识别组学类型、加载对应流水线、推进全流程。中途可以随时：

```
/bio status      # 查看进度
/bio continue    # 断点续传
/bio help        # 帮助
```

---

## 预期效率（基于设计预期）

| 阶段 | 传统方式（估计） | 使用框架（估计） | 说明 |
|:---|:---|:---|:---|
| QC + 标准化 + 聚类 | 3-7 天 | 数小时 | 标准流程自动化 |
| 细胞注释 | 2-5 天 | 1-2 天 | AI 辅助，但仍需人工验证 |
| 出版级图表 | 3-7 天 | 1-2 天 | 自动生成 + 视觉审查闭环 |
| 论文初稿 | 2-4 周 | 1-3 天 | Step 6 自动生成，需人工审阅 |

> **实际效率取决于数据复杂度、研究问题难度和用户的参与程度。** 以上数据来自框架的设计预期，不是经过大规模用户验证的统计结论。

---

## Bio-Framework 不适合什么

| 场景 | 原因 | 替代方案 |
|:---|:---|:---|
| 大规模生产级流水线 | 框架是交互式的，不适合每天处理数百样本的测序中心 | Nextflow / Snakemake |
| 底层算法开发 | 框架使用现有工具（Seurat、DESeq2 等），不适合方法学创新 | 直接编程 |
| 纯上游分析 | 基因组组装、比对（alignment）等计算密集型任务 | BWA / GATK / SPAdes |
| 不使用 Claude Code | 框架是 Claude Code 的 Skill 插件，依赖该平台 | — |

---

## 了解更多

| 你想了解… | 去这里 |
|:---|:---|
| 与"直接用 Claude"的深度对比 | [为什么要使用 Bio-Framework](docs/why-bio-framework-zh.md) |
| 详细功能展示与使用场景 | [产品展示](docs/product-showcase-zh.md) |
| 安装与配置 | [安装指南](docs/installation-guide-zh.md) |
| 架构与技术实现 | [架构文档](docs/architecture-zh.md) |
| 所有 60+ 命令参考 | [命令系统](framework/manifest.yml) |

---

## 许可证

本仓库采用 [Bio-Framework 公开评估许可证（非商业）](LICENSE.md) 发布。

你可以克隆和查阅本仓库用于个人学习和非商业学术研究。**商业使用、二次分发、衍生作品及 Prompt 框架提取均被严格禁止。** 完整条款请查看 [LICENSE.md](LICENSE.md)。

---

## 获取 Bio-Framework

完整商业版本（350KB+ 知识库、16 种错误模式、完整 Pipeline 自动化）可通过以下渠道购买：

**[在 Gumroad 获取 Bio-Framework →](https://howler26873.gumroad.com/l/gspdf)**

一次购买，所有更新。无论你是刚入门生信的湿实验室博士，还是需要标准化分析流程的课题组 PI——Bio-Framework 让 AI 按照你的领域标准工作。
