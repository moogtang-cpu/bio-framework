# Bio-Framework v1.0.0 系统架构与功能说明

**版本**：v1.0.0
**文档日期**：2026-03-10
**定位**：Claude Code 技能插件，面向生物信息学端到端分析

---

## 目录

1. [系统架构概览](#一系统架构概览)
   - 1.1 [三层资源架构](#11-三层资源架构)
   - 1.2 [核心模块架构](#12-核心模块架构)
   - 1.3 [知识库系统](#13-知识库系统)
2. [功能说明](#二功能说明)
   - 2.1 [端到端工作流引擎](#21-端到端工作流引擎)
   - 2.2 [四级思维深度系统](#22-四级思维深度系统)
   - 2.3 [命令系统](#23-命令系统)
   - 2.4 [知识积累与错误处理](#24-知识积累与错误处理)
   - 2.5 [检查点恢复与会话管理](#25-检查点恢复与会话管理)
3. [系统特点](#三系统特点)
   - 3.1 [自适应分析引擎](#31-自适应分析引擎)
   - 3.2 [生产级质量保证](#32-生产级质量保证)
   - 3.3 [完全本地化运行](#33-完全本地化运行)
   - 3.4 [模块化可扩展](#34-模块化可扩展)
   - 3.5 [SCI 论文全流程支持](#35-sci-论文全流程支持)
   - 3.6 [多语言国际化](#36-多语言国际化)
4. [架构图](#四架构图)
   - 4.1 [数据流向图](#41-数据流向图)
   - 4.2 [命令路由图](#42-命令路由图)
   - 4.3 [资源查找图](#43-资源查找图)

---

## 一、系统架构概览

### 1.1 三层资源架构

Bio-Framework 采用三层优先级架构，高优先级层的配置覆盖或合并低优先级层的配置，确保灵活的定制能力同时保持框架默认行为。这种设计使框架能够在不同规模的使用场景下保持一致性：插件层提供经过验证的领域默认值，用户全局层沉淀个人最佳实践，项目层支持针对特定数据集的精细调整，三层互不干扰，修改任一层都不影响其他层的设置。

```
Project Layer (.claude/)                  <- 优先级 1（最高，项目特定配置）
    |
    | 未找到时向下查找
    v
User Global (~/.claude/bio-framework/)    <- 优先级 2（用户全局偏好）
    |
    | 未找到时向下查找
    v
Plugin Layer (framework/)                 <- 优先级 3（框架默认值）
```

**各层职责说明：**

| 层级 | 路径 | 典型用途 |
|------|------|---------|
| Project Layer | `<project>/.claude/` | 项目特定的规则、知识库、参考代码 |
| User Global Layer | `~/.claude/bio-framework/` | 用户个人偏好参数（如 QC 阈值、聚类分辨率） |
| Plugin Layer | `framework/` | 框架内置默认值、标准流水线、通用知识库 |

**合并策略：**

| 资源类型 | 合并策略 | 说明 |
|---------|---------|------|
| 技能文件（SKILL.md） | override | 高优先级层文件完全覆盖低优先级层 |
| 规则参数（rules/） | override | 参数值取高优先级层定义，未定义项回退到低层 |
| 参考代码（ref_code/） | override | 函数库定义以高优先级层为准 |
| 知识库（knowledge/） | union | 所有层的内容合并，Project 层优先处理冲突 |
| 触发词（triggers/） | union | 所有层的触发词合并，扩展识别范围 |
| 错误模式（COMMON_ERRORS） | union | 所有层的错误模式合并，丰富诊断能力 |

---

### 1.2 核心模块架构

框架核心由技能入口、主编排器及其子模块、以及多个专职控制器组成：

```
framework/
├── templates/skills/bio/
│   └── SKILL.md                    <- 技能入口（551 行，自动加载）
│
└── skills/
    ├── ORCHESTRATOR.md             <- 主编排器（398 行）
    │   └── orchestrator/           <- 子模块目录（13 个，共 10,318 行）
    │       ├── execution_steps.md          执行步骤规范（443 行）
    │       ├── checkpoint_recovery.md      检查点恢复（271 行）
    │       ├── data_quality_review.md      Step 1.5 质量审计（1,840 行）
    │       ├── reflection_exploration.md   Step 2 反思探索（1,128 行）
    │       ├── research_summary.md         Step 3 研究总结（1,658 行）
    │       ├── publication_figures.md      Step 4 出版图表（1,852 行）
    │       ├── manuscript_drafting.md      Step 6 SCI 稿件（1,574 行）
    │       ├── final_report.md             Step 5 技术报告（221 行）
    │       ├── thinking_dispatch.md        思维模式分配（343 行）
    │       ├── pipeline_compliance.md      管道合规检查（489 行）
    │       ├── NEW_TOPIC_CLEANUP.md        新主题清理（含于 skills/）
    │       ├── topic_change_rules.md       主题变更规则（302 行）
    │       └── claude_router.md            模型路由（95 行）
    │
    ├── WORKFLOW_CONTROLLER.md      <- 工作流状态控制（命令执行调度）
    ├── COMMAND_PARSER.md           <- 命令解析 + 多语言路由
    ├── HELP_SYSTEM.md              <- 帮助信息系统（三语言）
    ├── THINKING_MODES.md           <- 四级思维深度调度
    ├── TOPIC_LOADING.md            <- 课题解析和理解
    ├── CONFIG_MANAGER.md           <- 参数配置管理
    ├── AUTO_EXECUTE.md             <- 自动执行策略
    ├── decision_nodes.md           <- 决策节点定义
    └── NEW_TOPIC_CLEANUP.md        <- 新主题初始化清理
```

**工作流状态枚举（`current_step`）：**

| 枚举值 | 含义 |
|--------|------|
| `init` | 项目初始化（Step 0） |
| `phase_running` | 动态分析阶段进行中（Phase 1-N） |
| `analysis_complete` | 所有分析阶段完成 |
| `step1_5_quality_review` | Step 1.5 数据质量审计 |
| `step2_reflection` | Step 2 反思与探索 |
| `step3_summary` | Step 3 研究总结与图表规划 |
| `step4_figures` | Step 4 出版质量图表生成 |
| `step5_report` | Step 5 技术报告 |
| `step6_manuscript` | Step 6 SCI 稿件草稿 |
| `completed` | 全流程完成 |

**`workflow_state` 状态枚举：**`pending | running | paused | completed | failed`

---

### 1.3 知识库系统

```
framework/knowledge/
├── COMMON_ERRORS.yml           <- 20KB，索引化错误模式（16 种已知模式）+ 验证解决方案
├── figure_types.yml            <- 68KB，7 种组学的标准图形类型定义
├── sci_design_standards.yml    <- 36KB，SCI 论文视觉与设计标准
├── task_quality_signals.yml    <- 质量信号和评估阈值
├── pipelines/                  <- 6 个标准分析流水线
│   ├── scrna_seq.yml               单细胞 RNA-seq
│   ├── bulk_rnaseq.yml             Bulk RNA-seq
│   ├── spatial_transcriptomics.yml 空间转录组
│   ├── proteomics.yml              蛋白质组学
│   ├── metabolomics.yml            代谢组学
│   ├── lipidomics.yml              脂质组学
│   └── _schema.yml                 流水线 Schema 定义
├── journals/                   <- 6 个期刊合规配置
│   ├── nature.yml                  Nature
│   ├── science.yml                 Science
│   ├── cell_reports.yml            Cell Reports
│   ├── nature_methods.yml          Nature Methods
│   ├── jci.yml                     JCI
│   ├── plos_one.yml                PLOS ONE
│   └── _schema.yml                 期刊配置 Schema
└── methods/                    <- 分析方法论文档
```

**知识库核心功能：**

- `figure_types.yml`：定义每种组学分析类型的标准图形（`mandatory: true` 图形在 Step 3 验证门中强制检查）
- `COMMON_ERRORS.yml`：结构化错误模式索引，支持三层错误解决机制
- `sci_design_standards.yml`：内置 SCI 出版标准，驱动 Step 4 图表质量判断
- `journals/`：各期刊的数据要求、格式规范、字数限制，支持 Step 6 合规检查

---

## 二、功能说明

### 2.1 端到端工作流引擎

Bio-Framework 实现从原始数据到 SCI 稿件的完整分析流程，共 7 个阶段：

#### Step 0：项目初始化

- 读取 `TOPIC.yml`（课题定义文件）
- 调用 `TOPIC_LOADING.md` 进行课题解析
- 动态生成分析计划（Phase 数量和内容由 AI 根据组学类型自动决定）
- 创建运行时目录结构，初始化 `workflow_state.yml` 和 `findings.yml`

#### Phase 1-N：动态分析阶段

- 阶段数量和内容完全由 AI 根据课题动态设计，不预设固定数量
- 每个 Phase 完成后写入检查点，支持跨会话恢复
- 异常发现实时追踪到 `findings.yml`
- 执行策略：遵循 `execution_policies.md` 定义的自动执行规则

#### Step 1.5：数据质量深度审计

五维度质量审计（由 `data_quality_review.md` 驱动）：

| 维度 | 审计内容 |
|------|---------|
| 样本身份验证 | 样本标签、元数据完整性、性别/组织标记一致性 |
| 批次效应检测 | 批次变量识别、批次效应量化评估 |
| 生物学合理性 | 已知标记基因表达、细胞类型分布合理性 |
| 离群值识别 | 统计离群样本/细胞识别与处置建议 |
| 跨阶段一致性 | 各 Phase 结论是否相互支持，无矛盾 |

#### Step 2：反思与探索循环

- 间隙分析：识别现有结果的证据空白
- 用户决策门：展示候选探索方向，由用户选择
- 补充执行：针对选中方向执行补充分析
- 重评估：更新研究叙事，决定是否再次循环

#### Step 3：研究总结与图表规划

- 故事合成：将所有 Phase 结论整合为连贯的研究故事
- 强制图形验证：读取 `figure_types.yml` 中 `mandatory: true` 的图形类型，确保关键标准图不遗漏
- 图表需求分析：生成 `figure_plan.yml`，定义每张图的数据来源、设计方案、叙事价值
- 完成度评估：加权公式计算分析完成度，低于阈值触发 HARD_STOP

#### Step 4：出版质量图表

四阶段图表生产流程：

1. **4a-4c 代码生成与执行**：按 `figure_plan.yml` 生成 R/Python 图表代码并执行
2. **4d AI 视觉检查**：AI 通过 Read 工具实际查看生成图像，对照 SCI 设计标准评审
3. **4e 修正迭代**：根据视觉检查结果自动修正代码并重新生成
4. **4f 用户视觉审查门**：展示最终图表给用户确认，记录视觉审查日志（`visual_review_log.yml`）

#### Step 5：技术报告

- 完整方法学文档（参数设置、版本信息）
- 可重现性包（代码 + 环境配置）
- `quality_report.md` 生成

#### Step 6：SCI 稿件草稿

六步精细化流程：

1. 初稿生成（基于研究总结和图表说明自动起草）
2. 期刊合规检查（对照目标期刊配置文件验证）
3. 语言精炼（学术写作规范优化）
4. 结构优化（逻辑流、过渡段调整）
5. 数据引用核查（确保统计数字与原始结果一致）
6. 五阶段提交前验证（`submission_readiness.md` 生成）

---

### 2.2 四级思维深度系统

Bio-Framework 根据任务复杂度自动调度四种思维深度，由 `THINKING_MODES.md` 和 `thinking_dispatch.md` 协同控制：

| 级别 | 标签 | Token 预算 | 典型应用场景 |
|------|------|-----------|------------|
| Quick | `[QUICK]` | 低 | 文件操作、状态更新、格式转换 |
| Standard | `[STANDARD]` | 中 | 已知模式的代码生成、参数配置 |
| Deep | `[DEEP]` | 高 | 质量评估、错误诊断、结果解读 |
| Ultrathink | `[ULTRATHINK]` | 最高 | 主题理解、细胞注释、综合结论 |

**9 项强制性 Ultrathink 场景（必须使用最高深度思维）：**

1. 课题初始化与分析计划生成（Step 0）
2. 细胞类型注释与亚群定义
3. 批次效应评估与处置决策
4. 跨 Phase 结论综合（Step 1.5）
5. 研究故事合成（Step 3）
6. 图表叙事规划（Step 3 figure_plan）
7. 探索方向决策（Step 2）
8. 稿件核心论点提炼（Step 6）
9. 提交前最终质量评估（Step 6 submission_readiness）

---

### 2.3 命令系统

命令系统由 `COMMAND_PARSER.md` 解析，`WORKFLOW_CONTROLLER.md` 调度执行：

**规模：**
- 25+ 顶级命令
- 50+ 子命令
- 110+ 多语言别名（英文 / 中文 / 日文）

**命令前缀（等价）：**

```
/bio          <- 标准英文前缀
/bioinformatics
/生信          <- 中文前缀
/バイオ         <- 日文前缀
```

**核心命令分类：**

| 类别 | 代表命令 | 功能 |
|------|---------|------|
| 工作流控制 | `/bio start`, `/bio continue`, `/bio pause` | 启动、继续、暂停分析 |
| 状态查询 | `/bio status`, `/bio progress` | 查看当前进度和阶段状态 |
| 步骤跳转 | `/bio goto step3`, `/bio step4` | 直接进入指定步骤 |
| 稿件操作 | `/bio manuscript`, `/bio journal` | 稿件生成和期刊设置 |
| 配置管理 | `/bio config`, `/bio rules` | 查看和修改参数配置 |
| 帮助系统 | `/bio help`, `/帮助` | 多语言帮助信息 |

**自然语言回退：** 当命令不匹配任何已知命令时，COMMAND_PARSER 调用自然语言识别逻辑，尝试从用户输入中推断意图，回退到最接近的命令。

---

### 2.4 知识积累与错误处理

**三层错误解决机制：**

```
第 1 层：框架知识库（COMMON_ERRORS.yml）
    |- 索引化错误模式（按错误类型、组学类型、工具分类）
    |- 验证过的解决方案
    |- 覆盖常见 R/Python 环境错误、数据格式错误、分析逻辑错误

第 2 层：用户经验库（User Global / Project Layer）
    |- 用户在历史分析中积累的特定错误解决方案
    |- 项目特定的数据格式规范

第 3 层：AI 自主推理
    |- 前两层均未命中时，AI 基于上下文推理解决方案
    |- 复杂问题自动保存解决方案到用户经验库（自动学习）
```

**自动学习机制：** 当 AI 通过自主推理解决了一个复杂且可复用的问题时，框架自动将问题模式和解决方案写入用户经验库，供后续分析复用。

---

### 2.5 检查点恢复与会话管理

框架设计为支持跨会话长期运行，单次会话中断不丢失进度：

**核心机制（由 `checkpoint_recovery.md` 定义）：**

- 每个 Phase 完成后自动写入检查点（`workflow_state.yml` 更新）
- 新会话启动时，ORCHESTRATOR 自动检测 `workflow_state.yml` 并提示恢复
- 恢复策略：根据 `current_step` 枚举值定位到具体的恢复入口

**关键运行时文件：**

| 文件 | 路径 | 用途 |
|------|------|------|
| `workflow_state.yml` | `<project>/.claude/` | 工作流状态持久化 |
| `findings.yml` | `<project>/.claude/` | 跨 Phase 异常发现追踪 |
| `figure_plan.yml` | `<project>/.claude/` | 图表规划（Step 3 输出） |
| `visual_review_log.yml` | `<project>/.claude/` | 视觉审查记录（Step 4 输出） |
| `submission_readiness.md` | `<project>/output/` | 提交就绪检查单（Step 6 输出） |

**暂停与恢复条件：** `AUTO_EXECUTE.md` 定义了触发暂停的条件（如数据质量异常、关键决策点、用户确认门），确保自动化执行与人工干预的平衡。

---

## 三、系统特点

### 3.1 自适应分析引擎

Bio-Framework 不预设固定的分析阶段数，而是由 AI 在 Step 0 根据课题内容动态设计分析计划：

- Phase 数量：AI 根据组学类型、研究问题复杂度、数据规模自动决定（通常 3-8 个 Phase）
- 分析内容：每个 Phase 的具体分析步骤由 AI 动态规划，不限于框架内置流水线
- 参数调整：AI 根据数据特征（如细胞数量、基因数量、批次数）自动建议调整 QC 阈值和分析参数
- 组学类型：支持任意组学类型分析，包括框架内置流水线未覆盖的新兴分析类型

### 3.2 生产级质量保证

框架本身经过多轮深度审计和测试：

| 指标 | 数值 |
|------|------|
| 自动化测试用例 | 898+ 个，100% 通过率 |
| 深度审计会话 | 76+（S12 至 S76） |
| 历史 Bug 修复 | 400+ 个 |
| 一致性自动检查 | 13 项（consistency_checker.py） |
| 核心文件覆盖率 | 100%（S72 全仓库 77 文件覆盖） |

一致性检查器（`framework/scripts/consistency_checker.py`）在每次审计后自动运行 13 项跨文件一致性检查，确保枚举值、字段引用、Schema 定义在所有模块间保持同步。

### 3.3 完全本地化运行

- 所有数据分析在本地 R/Python 环境中执行，代码由 Claude Code 工具链调用
- 原始数据零传输：原始数据、中间结果、最终产出物均不离开本地环境（Claude Code 仅将代码和输出摘要发送至 Anthropic API 进行模型处理）
- 零云端存储：适合处理临床数据、患者数据、保密研究数据
- 环境隔离：通过 `environment_isolation.md` 定义的策略管理 R/Python 包安装，避免污染全局环境

### 3.4 模块化可扩展

**三层定制能力：**

- **项目级定制（Project Layer）**：为特定项目添加专用知识库、调整分析参数、定义项目特定触发词
- **用户级定制（User Global Layer）**：设置个人偏好参数（如 `~/.claude/rules/user_rules.yml`），对所有项目生效
- **扩展知识库**：在任意层添加 `knowledge/` 目录，框架自动合并；新增期刊配置文件即可扩展合规检查
- **新分析类型**：在 `framework/knowledge/pipelines/` 添加新的流水线 YAML 即可支持新组学类型

### 3.5 SCI 论文全流程支持

Bio-Framework 提供从数据分析到论文初稿的一站式支持：

- **标准图强制验证**：`figure_types.yml` 中标记为 `mandatory: true` 的图形（如单细胞 UMAP 细胞类型大图、空间散点图）在 Step 3 验证门中强制检查，缺失则触发 HARD_STOP
- **期刊合规**：内置 6 个主要期刊（Nature、Science、Cell Reports、Nature Methods、JCI、PLOS ONE）的数据要求、格式规范配置
- **视觉审查闭环**：Step 4 要求 AI 实际查看生成图像（非仅检查文件元数据），并通过用户视觉审查门确认
- **SCI 设计标准**：`sci_design_standards.yml` 内置 36KB 的 SCI 出版视觉标准，驱动图表质量判断
- **稿件精细化**：Step 6 包含 6 轮精细化和 5 阶段提交前验证

### 3.6 多语言国际化

框架支持英文、中文、日文三语言操作：

| 语言 | 命令别名数量 | 帮助系统 | 触发词组 |
|------|-----------|---------|---------|
| 英文（en） | 标准命令 | 完整支持 | 标准英文触发词 |
| 中文（zh） | 87 个别名 | 完整支持 | 10 组中文触发词 |
| 日文（ja） | 已声明 | 完整支持 | 日文触发词 |

- **自动语言检测**：COMMAND_PARSER 自动识别输入语言，无需用户指定
- **语言回退**：未知命令按检测到的语言进行自然语言识别
- **本地化帮助**：`HELP_SYSTEM.md` 提供三语言完整帮助内容

---

## 四、架构图

### 4.1 数据流向图

```
用户输入（TOPIC.yml + 原始数据）
    |
    v
+------------------+
| SKILL.md         |  技能入口，自动加载配置，触发 ORCHESTRATOR
| (551 行)         |
+------------------+
    |
    v
+------------------+
| ORCHESTRATOR.md  |  主编排器，协调全局工作流
| + 13 子模块      |
+------------------+
    |
    +----------+----------+----------+----------+----------+
    |          |          |          |          |          |
    v          v          v          v          v          v
[Step 0]   [Phase     [Step 1.5]  [Step 2]  [Step 3]  [Step 4]
 初始化     1-N]       质量审计    反思探索   研究总结   出版图表
           动态分析
    |          |          |          |          |          |
    v          v          v          v          v          v
[workflow   [findings  [audit     [gap       [figure    [visual
 _state.yml] .yml]      _report]   _report]   _plan.yml] _review
                                              table      _log.yml]
                                              _plan.yml]
    |
    +----------+----------+
    |          |          |
    v          v          v
[Step 5]   [Step 6]  [completed]
 技术报告   SCI 稿件
    |          |
    v          v
[quality   [submission
 _report   _readiness
 .md]      .md]
    |          |
    v          v
    最终产出物（本地）
```

### 4.2 命令路由图

```
用户输入（任意语言）
    |
    v
+------------------+
| COMMAND_PARSER   |  词法分析 + 多语言别名映射
| .md              |  前缀识别：/bio /bioinformatics /生信 /バイオ
+------------------+
    |
    +------------------+------------------+
    |                  |                  |
    v                  v                  v
[精确命令匹配]    [多语言别名匹配]    [自然语言回退]
    |                  |                  |
    v                  v                  v
    +------------------+------------------+
                       |
                       v
          +------------------+
          | WORKFLOW_        |  状态检查，确定执行路径
          | CONTROLLER.md    |  读取 workflow_state.yml
          +------------------+
               |        |
               v        v
         [工作流命令]  [非工作流命令]
          start/        config/
          continue/     help/
          pause 等      status 等
               |        |
               v        v
          +------------------+
          | ORCHESTRATOR.md  |  调度对应子模块执行
          +------------------+
               |
               v
          [对应子模块]
          execution_steps /
          data_quality_review /
          publication_figures 等
               |
               v
          [本地 R/Python 执行环境]
               |
               v
          [输出文件写入本地目录]
```

### 4.3 资源查找图

```
需要加载资源（如 clustering.resolution 参数）
    |
    v
+--------------------------------------+
| Step 1：查找 Project Layer           |
| <project>/.claude/rules/             |
| <project>/.claude/knowledge/         |
+--------------------------------------+
    |
    | 找到？--> YES --> 使用该值，停止查找
    |
    | NO
    v
+--------------------------------------+
| Step 2：查找 User Global Layer       |
| ~/.claude/bio-framework/rules/        |
| ~/.claude/bio-framework/knowledge/    |
| ~/.claude/rules/user_rules.yml        |
+--------------------------------------+
    |
    | 找到？--> YES --> 使用该值，停止查找
    |
    | NO
    v
+--------------------------------------+
| Step 3：查找 Plugin Layer            |
| framework/rules/default_rules.yml    |
| framework/knowledge/                 |
| framework/defaults/                  |
+--------------------------------------+
    |
    | 找到？--> YES --> 使用框架默认值
    |
    | NO（参数未定义）
    v
+--------------------------------------+
| AI 内置默认值 / 触发用户确认         |
+--------------------------------------+

合并规则（override 类型资源）：
    Project > User Global > Plugin（高层覆盖低层）

合并规则（union 类型资源）：
    Project + User Global + Plugin（所有层内容合并）
    冲突时以高优先级层条目为准
```

---

*文档结束*

*本文档描述 Bio-Framework v1.0.0 内部架构，面向框架开发者和高级用户。*
*普通用户使用指南请参考 `docs/user_guides/` 目录。*
