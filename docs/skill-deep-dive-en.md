# Bio-Framework v1.0.0 --- Skill Deep Dive

**Document Type**: Technical Deep Dive
**Target Audience**: Developers, advanced users, and contributors who want to understand internal mechanisms
**Framework Version**: 1.0.0
**Last Updated**: 2026-03-10

---

## Table of Contents

1. [Architecture](#1-architecture)
2. [Workflow](#2-workflow)
3. [Workflow and Thinking Mode Matrix](#3-workflow--thinking-mode-matrix)
4. [Core Mechanism Deep Dive](#4-core-mechanism-deep-dive)
5. [ORCHESTRATOR Deep Dive](#5-orchestrator-deep-dive)
6. [claude-router Integration](#6-claude-router-integration)
7. [Command System](#7-command-system)
8. [State Management](#8-state-management)
9. [Data Security Architecture](#9-data-security-architecture)
10. [v1.0.0 Major Evolution from Earlier Versions](#10-v100-major-evolution-from-earlier-versions)

---

## 1. Architecture

### 1.1 Core Architecture --- Layered + Cyclic + Adaptive

Bio-Framework is organized as a five-layer architecture. Static layers provide immutable
logic and knowledge; dynamic layers are generated and mutated by the AI at runtime.
The key insight is that the **Framework Layer** and **Knowledge Layer** are never modified
during execution --- they serve as the invariant foundation. All runtime state flows
through the **Configuration Layer** and **Execution Layer**, with final deliverables
landing in the **Output Layer**.

```
┌────────────────────────────────────────────────────────────────────────────┐
│                    BIO-FRAMEWORK v1.0.0 ARCHITECTURE                     │
├────────────────────────────────────────────────────────────────────────────┤
│                                                                          │
│  ┌────────────────────────────────────────────────────────────────────┐   │
│  │  Layer 1: Framework Layer (framework/skills/)        [STATIC]     │   │
│  │                                                                    │   │
│  │  SKILL.md ─────────────── Entry point + triggers (551 lines)      │   │
│  │       │                                                            │   │
│  │       ├── ORCHESTRATOR.md ─── Main orchestrator (398 lines)       │   │
│  │       │       └── 13 sub-modules (10,318 lines total)             │   │
│  │       │           ├── execution_steps.md       (443 lines)        │   │
│  │       │           ├── data_quality_review.md   (1,840 lines)      │   │
│  │       │           ├── reflection_exploration.md(1,128 lines)      │   │
│  │       │           ├── research_summary.md      (1,658 lines)      │   │
│  │       │           ├── publication_figures.md    (1,852 lines)      │   │
│  │       │           ├── manuscript_drafting.md    (1,574 lines)      │   │
│  │       │           ├── pipeline_compliance.md   (489 lines)        │   │
│  │       │           ├── thinking_dispatch.md     (343 lines)        │   │
│  │       │           ├── checkpoint_recovery.md   (271 lines)        │   │
│  │       │           ├── final_report.md          (221 lines)        │   │
│  │       │           ├── topic_change_rules.md    (302 lines)        │   │
│  │       │           ├── integration_guide.md     (102 lines)        │   │
│  │       │           └── claude_router.md         (95 lines)         │   │
│  │       │                                                            │   │
│  │       ├── WORKFLOW_CONTROLLER.md (1,102 lines)                    │   │
│  │       ├── COMMAND_PARSER.md      (556 lines)                      │   │
│  │       ├── HELP_SYSTEM.md         (769 lines)                      │   │
│  │       ├── THINKING_MODES.md      (456 lines)                      │   │
│  │       ├── TOPIC_LOADING.md       (720 lines)                      │   │
│  │       └── CONFIG_MANAGER.md      (802 lines)                      │   │
│  └────────────────────────────────────────────────────────────────────┘   │
│       │                                                                   │
│       ▼                                                                   │
│  ┌────────────────────────────────────────────────────────────────────┐   │
│  │  Layer 2: Knowledge Layer (framework/knowledge/) [STATIC+EXTEND]  │   │
│  │                                                                    │   │
│  │  COMMON_ERRORS.yml ────── Error pattern library (20 KB)           │   │
│  │  figure_types.yml ─────── Standard figure definitions (68 KB)     │   │
│  │  sci_design_standards.yml  SCI design standards (36 KB)           │   │
│  │  pipelines/ ──────────── 6 standard analysis pipelines            │   │
│  │    ├── scrna_seq.yml        Single-Cell RNA-seq                   │   │
│  │    ├── bulk_rnaseq.yml      Bulk RNA-seq                          │   │
│  │    ├── spatial_transcriptomics.yml  Spatial Transcriptomics       │   │
│  │    ├── proteomics.yml       Proteomics                            │   │
│  │    ├── metabolomics.yml     Metabolomics                          │   │
│  │    ├── lipidomics.yml       Lipidomics                            │   │
│  │    └── _schema.yml          Pipeline schema definition            │   │
│  │  (Note: ATAC-seq covered in figure_types.yml, no pipeline file)   │   │
│  │  journals/ ──────────── 6 journal compliance configs              │   │
│  │    ├── nature.yml / cell_reports.yml / jci.yml                    │   │
│  │    ├── nature_methods.yml / plos_one.yml / science.yml            │   │
│  │    └── _schema.yml                                                │   │
│  └────────────────────────────────────────────────────────────────────┘   │
│       │                                                                   │
│       ▼                                                                   │
│  ┌────────────────────────────────────────────────────────────────────┐   │
│  │  Layer 3: Configuration Layer (project_skills/)  [AI-GENERATED]   │   │
│  │                                                                    │   │
│  │  TOPIC.yml ──────────── Research topic definition                 │   │
│  │  phases/                                                           │   │
│  │    ├── phase_index.yml ── Phase index (generated by Step 0)       │   │
│  │    └── phaseNN_xxx.yml ── Individual phase configs                │   │
│  │  runtime/                                                          │   │
│  │    ├── workflow_state.yml  Current execution state                │   │
│  │    ├── phase_summaries.yml  Phase summary accumulator             │   │
│  │    └── findings.yml ────── Findings tracker                       │   │
│  └────────────────────────────────────────────────────────────────────┘   │
│       │                                                                   │
│       ▼                                                                   │
│  ┌────────────────────────────────────────────────────────────────────┐   │
│  │  Layer 4: Execution Layer (7-step workflow + thinking dispatch)    │   │
│  │                                                                    │   │
│  │  Step 0  ──► Phase 1-N  ──► Step 1.5  ──► Step 2  ──►            │   │
│  │  Step 3  ──► Step 4  ──► Step 5  ──► Step 6  ──► COMPLETED       │   │
│  │                                                                    │   │
│  │  Thinking Dispatch: Quick / Standard / Deep / Ultrathink          │   │
│  │  Auto-Escalation: Quick ──► Standard ──► Deep ──► Ultrathink     │   │
│  └────────────────────────────────────────────────────────────────────┘   │
│       │                                                                   │
│       ▼                                                                   │
│  ┌────────────────────────────────────────────────────────────────────┐   │
│  │  Layer 5: Output Layer (Phase_output/)                            │   │
│  │                                                                    │   │
│  │  Phase_output/                                                     │   │
│  │    ├── phase{NN}_{name}/current/   Phase analysis outputs         │   │
│  │    ├── data_quality_review/current/ Step 1.5 quality report       │   │
│  │    ├── exploration_synthesis/current/ Step 2 synthesis            │   │
│  │    ├── research_summary/current/   Step 3 story + figure plan     │   │
│  │    ├── publication_figures/        Step 4 figures + tables         │   │
│  │    ├── final_report/              Step 5 methods + reproducibility │   │
│  │    └── manuscript/                Step 6 SCI manuscript            │   │
│  └────────────────────────────────────────────────────────────────────┘   │
│                                                                          │
└────────────────────────────────────────────────────────────────────────────┘
```

**Design rationale**: The separation between static layers (Framework + Knowledge) and
dynamic layers (Configuration + Execution + Output) means the framework can be updated
independently of any running project. The Knowledge Layer is marked `STATIC+EXTENSIBLE`
because users can add custom pipelines or journal configs without modifying core logic.

---

### 1.2 Three-Layer Resource Priority

All configuration values, rules, templates, and knowledge files are resolved through a
three-layer priority system. Higher layers override lower layers for the same key.

```
┌──────────────────────────────────────────────────────────┐
│              THREE-LAYER RESOLUTION ORDER                │
├──────────────────────────────────────────────────────────┤
│                                                          │
│   [HIGHEST]  Project Layer                               │
│              .claude/rules/project_rules.yml             │
│              Purpose: Project-specific overrides          │
│                        │                                 │
│                        ▼ (if not found)                  │
│   [MEDIUM]   User Global Layer                           │
│              ~/.claude/bio-framework/rules/user_rules.yml│
│              Purpose: Personal preferences               │
│                        │                                 │
│                        ▼ (if not found)                  │
│   [LOWEST]   Plugin Layer                                │
│              framework/rules/default_rules.yml           │
│              Purpose: Framework defaults                 │
│                                                          │
└──────────────────────────────────────────────────────────┘
```

**Merge strategy by field type**:

| Field Type | Strategy | Example |
|------------|----------|---------|
| Scalar values (numbers, strings) | Override (higher wins) | `mt_percent: 15` overrides default `20` |
| Lists (arrays) | Override (entire list replaced) | Custom `skip_phases: [3, 5]` replaces default `[]` |
| Dictionaries (maps) | Shallow merge (key-level override) | Project `clustering.resolution: 1.2` overrides User Global `0.6` |
| Boolean flags | Override | `quality_review.enabled: false` overrides default `auto` |

**Resolution example** (clustering resolution):

```yaml
# Plugin Layer (framework/rules/default_rules.yml)
clustering:
  resolution:
    value: 0.8       # Framework default

# User Global Layer (~/.claude/bio-framework/rules/user_rules.yml)
clustering:
  resolution:
    value: 0.6       # User preference: finer clusters

# Project Layer (.claude/rules/project_rules.yml)
clustering:
  resolution:
    value: 1.2       # This dataset needs higher resolution

# Resolution result: 1.2 (Project layer wins)
```

---

### 1.3 ORCHESTRATOR Sub-Module Architecture

The ORCHESTRATOR (398 lines) serves as the central dispatch hub, delegating detailed
logic to 13 specialized sub-modules totaling 10,318 lines. Sub-modules are loaded
**on demand** --- only the relevant module is read when a workflow step is entered.

```
ORCHESTRATOR.md (398 lines) ─── Central Dispatch
│
├── execution_steps.md ──────────── (443 lines)  Workflow step definitions
│                                                  7 core steps + Step 1.5
│
├── thinking_dispatch.md ────────── (343 lines)  Per-step thinking mode routing
│                                                  Maps tasks to Quick/Std/Deep/Ultra
│
├── data_quality_review.md ──────── (1,840 lines) [NEW] Step 1.5 deep audit
│                                                  5-dimension quality assessment
│
├── reflection_exploration.md ───── (1,128 lines) Step 2 gap analysis + exploration
│                                                  Max 2 rounds, 3 suggestions/round
│
├── research_summary.md ─────────── (1,658 lines) Step 3 story + figure planning
│                                                  Completeness validation gate
│
├── publication_figures.md ──────── (1,852 lines) Step 4 figure generation
│                                                  AI visual inspection + user review
│
├── manuscript_drafting.md ──────── (1,574 lines) Step 6 SCI manuscript
│                                                  21 independent sub-steps
│
├── pipeline_compliance.md ──────── (489 lines)   Soft pipeline validation
│                                                  Detection signals + compliance score
│
├── checkpoint_recovery.md ──────── (271 lines)   Session recovery
│                                                  4 recovery types (A/B/C/D)
│
├── final_report.md ─────────────── (221 lines)   Step 5 methods + reproducibility
│                                                  Parameter recording, env export
│
├── topic_change_rules.md ───────── (302 lines)   Topic switch governance
│                                                  Version protection rules
│
├── integration_guide.md ────────── (102 lines)   Third-party integration guide
│
└── claude_router.md ────────────── (95 lines)    Model selection integration
                                                   thinking_mode to model mapping
```

**Total framework skill system**: SKILL.md (551) + ORCHESTRATOR (398) + 7 top-level
skills (4,803) + 13 sub-modules (10,318) = **16,070 lines** of structured instructions.

---

## 2. Workflow

The execution workflow consists of 7 core steps plus a conditional Step 1.5.
Steps 0-1 handle initialization and dynamic analysis. Steps 1.5-6 handle
post-analysis processing. The entire workflow is designed to be **auto-executing** ---
the AI proceeds without asking for confirmation unless a defined HARD_STOP condition
is met.

```
┌──────────────────────────────────────────────────────────────────────────────────┐
│                         EXECUTION WORKFLOW (DETAILED)                           │
├──────────────────────────────────────────────────────────────────────────────────┤
│                                                                                │
│  ┌──────────────────────────────────────────────────────────────────────┐       │
│  │  STEP 0: PROJECT INITIALIZATION                    [ULTRATHINK]    │       │
│  │                                                                      │       │
│  │  0a. Read TOPIC.yml ─────── Understand scientific questions          │       │
│  │  0b. Read phase_index.yml ─ Understand planned workflow              │       │
│  │  0c. Read data_info.txt ─── Understand data types and locations      │       │
│  │  0d. Create virtual env ─── conda/venv (Conda-First policy)         │       │
│  │  0e. Init runtime/ ──────── workflow_state.yml, findings.yml         │       │
│  │  0f. Check data integrity ─ Verify files exist and are readable     │       │
│  │  0g. Pipeline compliance ── Match to standard pipeline (if exists)   │       │
│  └──────────────┬───────────────────────────────────────────────────────┘       │
│                 ▼                                                               │
│  ┌──────────────────────────────────────────────────────────────────────┐       │
│  │  STEP 1: ANALYSIS PHASE LOOP (Phase 1-N)          [MIXED MODES]    │       │
│  │                                                                      │       │
│  │  FOR each phase in phase_index.yml:                                  │       │
│  │    ├── Load phase_summaries.yml (context from previous phases)       │       │
│  │    ├── Load current phaseNN_xxx.yml config                           │       │
│  │    ├── Reasoning ──► Act ──► Observe ──► Record                     │       │
│  │    │   (thinking mode selected per task type --- see Section 3)      │       │
│  │    ├── Update phase_summaries.yml with results                      │       │
│  │    └── Update workflow_state.yml                                     │       │
│  │                                                                      │       │
│  │  After last phase:                                                   │       │
│  │    └── Evaluate Step 1.5 trigger conditions                         │       │
│  └──────────────┬───────────────────────────────────────────────────────┘       │
│                 ▼                                                               │
│  ┌──────────────────────────────────────────────────────────────────────┐       │
│  │  STEP 1.5: DATA QUALITY DEEP REVIEW (CONDITIONAL) [ULTRATHINK]     │       │
│  │  [NEW in v1.0.0]                                                     │       │
│  │                                                                      │       │
│  │  Trigger check:                                                      │       │
│  │    ├── TOPIC.yml quality_review.enabled == true? ──► Force trigger   │       │
│  │    ├── TOPIC.yml quality_review.enabled == false? ─► Force skip     │       │
│  │    └── Auto mode:                                                    │       │
│  │        ├── Clinical/method_comparison research? ──► Trigger          │       │
│  │        ├── quality_requirements >= high? ─────────► Trigger          │       │
│  │        └── Accumulated warnings >= 3? ────────────► Trigger          │       │
│  │                                                                      │       │
│  │  If triggered (5 audit dimensions):                                  │       │
│  │    D1. Sample Identity Verification                                  │       │
│  │    D2. Batch Effect Assessment                                       │       │
│  │    D3. Biological Plausibility                                       │       │
│  │    D4. Outlier Sample Detection                                      │       │
│  │    D5. Cross-Phase Consistency                                       │       │
│  │                                                                      │       │
│  │  Decision: PASS ──► proceed                                          │       │
│  │            WARNING ──► proceed with caution                          │       │
│  │            CRITICAL ──► HARD_STOP (user decision required)           │       │
│  │                                                                      │       │
│  │  If not triggered: skip directly to Step 2 (zero overhead)           │       │
│  └──────────────┬───────────────────────────────────────────────────────┘       │
│                 ▼                                                               │
│  ┌──────────────────────────────────────────────────────────────────────┐       │
│  │  STEP 2: REFLECTION & EXPLORATION LOOP            [ULTRATHINK]     │       │
│  │                                                                      │       │
│  │  2a. Gap Analysis ──────── Evaluate question completion levels       │       │
│  │      ├── Weighted completion scoring (per question importance)       │       │
│  │      ├── Unexpected finding detection                                │       │
│  │      ├── Quality signal detection (AI auto-adjust)                  │       │
│  │      └── Publication data readiness inventory                       │       │
│  │                                                                      │       │
│  │  2b. Exploration Gate ──── [PAUSE] User selects directions           │       │
│  │      ├── 3 exploration suggestions generated                        │       │
│  │      └── User: execute selected / skip all                          │       │
│  │                                                                      │       │
│  │  2c. Supplementary Execution ── Run selected explorations            │       │
│  │                                                                      │       │
│  │  2d. Re-evaluation ────── Integrate results, update completion       │       │
│  │      └── Still major gaps? ──► Return to 2b (max 2 rounds)          │       │
│  │                                                                      │       │
│  │  Threshold: Entry=80%, Continue-loop=70% (diminishing returns)       │       │
│  └──────────────┬───────────────────────────────────────────────────────┘       │
│                 ▼                                                               │
│  ┌──────────────────────────────────────────────────────────────────────┐       │
│  │  STEP 3: RESEARCH SUMMARY & FIGURE PLANNING       [ULTRATHINK]     │       │
│  │                                                                      │       │
│  │  3a. Research Story Synthesis ── Construct scientific narrative      │       │
│  │  3b. Figure Requirements ────── Literature-driven figure specs      │       │
│  │      └── Canonical figure WebSearch (for niche omics types)         │       │
│  │  3c. Figure Specification ────── figure_plan.yml output             │       │
│  │      └── Read research type figure framework (Module blueprint)     │       │
│  │  3d. Table Specification ─────── table_plan.yml output              │       │
│  │  3e. Completeness Validation Gate                                    │       │
│  │      └── mandatory_type_id_check (HARD_STOP if mandatory missing)   │       │
│  └──────────────┬───────────────────────────────────────────────────────┘       │
│                 ▼                                                               │
│  ┌──────────────────────────────────────────────────────────────────────┐       │
│  │  STEP 4: PUBLICATION FIGURES & TABLES              [DEEP]          │       │
│  │                                                                      │       │
│  │  4a. Figure Selection & Enhancement ── Apply color, font, annotate  │       │
│  │      └── WebSearch for visualization best practices                 │       │
│  │  4b. Supplementary Figure Generation ── Create missing figures      │       │
│  │  4c. Multi-panel Assembly ──────────── Composite A/B/C/D panels    │       │
│  │  4d. Figure Quality Check ──────────── Checklist + AI visual read   │       │
│  │      └── [NEW] AI reads PNG previews via Read tool                  │       │
│  │  4e. Supplementary Table Generation ── Per table_plan.yml           │       │
│  │  4f. User Visual Review ────────────── [PAUSE] Display for approval │       │
│  │      └── [NEW] User confirms before proceeding to Step 5            │       │
│  └──────────────┬───────────────────────────────────────────────────────┘       │
│                 ▼                                                               │
│  ┌──────────────────────────────────────────────────────────────────────┐       │
│  │  STEP 5: TECHNICAL FOUNDATION                      [DEEP]          │       │
│  │                                                                      │       │
│  │  5a. Methods Record ── parameters + tool versions ─► methods_record │       │
│  │  5b. Reproducibility Package ── env + scripts + data provenance     │       │
│  └──────────────┬───────────────────────────────────────────────────────┘       │
│                 ▼                                                               │
│  ┌──────────────────────────────────────────────────────────────────────┐       │
│  │  STEP 6: SCI MANUSCRIPT DRAFTING                   [ULTRATHINK]    │       │
│  │                                                                      │       │
│  │  Pre: Readiness Gate ───── Context reset + parameter validation     │       │
│  │  6a. Preparation & Outline                                           │       │
│  │  6b. Core Writing (6 sections):                                      │       │
│  │      6b-1 Title & Abstract    [ULTRATHINK]                          │       │
│  │      6b-2 Introduction        [ULTRATHINK]                          │       │
│  │      6b-3 Methods             [DEEP]                                │       │
│  │      6b-4 Results             [ULTRATHINK]                          │       │
│  │      6b-5 Discussion          [ULTRATHINK]                          │       │
│  │      6b-6 Figure Legends      [STANDARD]                           │       │
│  │  6c. Refinement (6 polishing passes)              [DEEP]            │       │
│  │  6d. Assembly & Quality Check (4 sub-steps)       [STANDARD]        │       │
│  │  6e-6i. Submission Readiness Validation:                             │       │
│  │      6e Data Validation       [ULTRATHINK]                          │       │
│  │      6f Claim Review          [STANDARD]                            │       │
│  │      6g Structural Check      [STANDARD]                            │       │
│  │      6h Cross-File Sync       [STANDARD]                            │       │
│  │      6i Readiness Report      [STANDARD]                            │       │
│  │                                                                      │       │
│  │  Total: 21 independent sub-steps (timeout-resilient)                 │       │
│  └──────────────┬───────────────────────────────────────────────────────┘       │
│                 ▼                                                               │
│           [COMPLETED]                                                           │
│                                                                                │
└──────────────────────────────────────────────────────────────────────────────────┘
```

**HARD_STOP pause points** (workflow auto-executes except at these gates):

| Pause Point | Step | Condition |
|-------------|------|-----------|
| Data unavailable | Step 0 | Data files do not exist, user needs to provide them |
| Data unsuitable | Step 0 | Format/sample size/quality issues require user decision |
| Manual review needed | Step 1 | Cell type annotation or similar biological judgment |
| Unresolvable error | Step 1 | Error not in COMMON_ERRORS.yml and auto-fix failed |
| Severely abnormal results | Step 1 | Any metric reaches HARD_STOP threshold per decision_nodes |
| Critical quality issues | Step 1.5 | Quality audit detects CRITICAL-level problems |
| Exploration Gate | Step 2b | User selects which explorations to run (or skip all) |
| Missing mandatory figures | Step 3e | mandatory_type_id_check fails (e.g., no UMAP for scRNA) |
| User Visual Review | Step 4f | Figures displayed; user must approve before Step 5 |
| Version change request | Any | Version field modifications require explicit user approval |

---

## 3. Workflow & Thinking Mode Matrix

### Per-Step Thinking Mode Assignment

Each workflow step and sub-task is assigned a thinking mode. The system automatically
selects the appropriate depth based on task type and can auto-escalate when anomalies
are detected.

| Step | Sub-Task | Mode | Rationale |
|------|----------|------|-----------|
| **Step 0** | Topic understanding | ULTRATHINK | Scientific question design requires deepest reasoning |
| Step 0 | Read phase_index | Standard | Structural understanding, no deep reasoning |
| Step 0 | Read data_info | Standard | Information extraction (upgrade to Deep if abnormal) |
| Step 0 | Check data / setup env | Quick | File existence checks, environment operations |
| Step 0 | Init runtime | Quick | File creation only |
| **Step 1** | Load context | Standard | Read and understand previous phase context |
| Step 1 | Code generation | Standard | Template-based (upgrade to Deep/Ultra for biology) |
| Step 1 | Code execution | Quick | Execution and monitoring only |
| Step 1 | Result interpretation | Deep | Multi-metric quality assessment |
| Step 1 | Cell annotation | ULTRATHINK | Core biological decision |
| Step 1 | Error diagnosis | Standard/Deep | Known pattern vs. unknown error |
| Step 1 | Phase summary | Standard | Structured recording |
| **Step 1.5** | Quality audit | ULTRATHINK | Cross-phase holistic data assessment |
| **Step 2** | Gap analysis | ULTRATHINK | Evaluate scientific question completeness |
| Step 2 | Exploration execution | Standard/Deep | Depends on task complexity |
| Step 2 | Re-evaluation | ULTRATHINK | Integrate exploration results holistically |
| **Step 3** | Story synthesis | ULTRATHINK | Construct coherent scientific narrative |
| Step 3 | Figure requirements | ULTRATHINK | Literature-driven specification |
| Step 3 | Figure specification | Deep | Detailed plan generation |
| Step 3 | Table specification | Deep | Data mapping and column design |
| **Step 4** | Figure enhancement | Deep | Visual and scientific quality judgment |
| Step 4 | Multi-panel assembly | Standard | Layout mechanics |
| Step 4 | Quality check + visual | Deep | AI visual inspection of generated figures |
| Step 4 | Table generation | Standard | Data compilation |
| **Step 5** | Methods record | Deep | Accurate parameter and version tracking |
| Step 5 | Reproducibility pkg | Deep | Environment and provenance capture |
| **Step 6** | Outline (6a) | ULTRATHINK | Manuscript structure design |
| Step 6 | Abstract/Intro/Results/Discussion | ULTRATHINK | Core scientific writing |
| Step 6 | Methods (6b-3) | Deep | Technical accuracy focus |
| Step 6 | Figure Legends (6b-6) | Standard | Structured description |
| Step 6 | Refinement (6c) | Deep | Language and consistency polishing |
| Step 6 | Assembly (6d) | Standard | File compilation |
| Step 6 | Data validation (6e) | ULTRATHINK | Statistical claim verification |
| Step 6 | Other validations (6f-6i) | Standard | Checklist-based verification |

### 9 Mandatory Ultrathink Scenarios

These scenarios **always** trigger Ultrathink mode, regardless of what thinking level
was previously active. This is enforced by the thinking_dispatch system.

```
┌──────────────────────────────────────────────────────────────────────┐
│                  9 MANDATORY ULTRATHINK SCENARIOS                   │
├──────────────────────────────────────────────────────────────────────┤
│                                                                      │
│  [1] New topic understanding and analysis plan design (ST001)        │
│      Trigger: Step 0 reads TOPIC.yml for the first time              │
│      Why: Entire analysis strategy depends on correct understanding  │
│                                                                      │
│  [2] Cell type annotation (ST002)                                    │
│      Trigger: Phase involves cell type annotation                    │
│      Why: Core biological decision requiring domain expertise        │
│                                                                      │
│  [3] Subpopulation definition and naming (ST003)                     │
│      Trigger: Defining or naming cell subpopulations                 │
│      Why: Requires literature knowledge and biological reasoning     │
│                                                                      │
│  [4] Trajectory analysis interpretation                              │
│      Trigger: Interpreting developmental trajectories or pseudotime  │
│      Why: Complex biological inference with multiple valid paths     │
│                                                                      │
│  [5] Cell communication result interpretation                        │
│      Trigger: CellChat/NicheNet/LIANA result analysis                │
│      Why: Multi-dimensional interaction data requires deep synthesis │
│                                                                      │
│  [6] Major analysis plan adjustment                                  │
│      Trigger: Significant deviation requiring plan revision          │
│      Why: Must re-evaluate entire analysis strategy holistically     │
│                                                                      │
│  [7] Project comprehensive restructuring (ST004)                     │
│      Trigger: Cross-phase integration of all findings                │
│      Why: Requires global perspective across all analysis phases     │
│                                                                      │
│  [8] Final conclusions and paper abstract (ST005)                    │
│      Trigger: Step 6b-1 (Title & Abstract) and Step 6b-5 (Discussion)│
│      Why: Paper-level rigor and scientific accuracy                  │
│                                                                      │
│  [9] Deep anomalous result analysis                                  │
│      Trigger: Results seriously inconsistent with literature (AT002) │
│      Why: Could be important discovery or serious error              │
│                                                                      │
└──────────────────────────────────────────────────────────────────────┘
```

### Auto-Escalation Rules

The thinking system supports automatic upgrade based on runtime conditions:

```yaml
# From THINKING_MODES.md - Auto-upgrade triggers
auto_upgrade_rules:
  scenario_triggers:
    - trigger: "New project framework generation"
      from: any -> to: ultrathink

    - trigger: "Cell type annotation decision"
      from: any -> to: ultrathink

    - trigger: "Subpopulation definition and naming"
      from: any -> to: ultrathink

  anomaly_triggers:
    - trigger: "Results deviate from expectation > 30%"
      from: standard -> to: deep

    - trigger: "Results seriously inconsistent with literature"
      from: deep -> to: ultrathink

    - trigger: "Abnormal cell type proportions (>3x known)"
      from: standard -> to: deep

    - trigger: "Unexpected cell subpopulation discovered"
      from: deep -> to: ultrathink

  failure_triggers:
    - trigger: "2 consecutive decision failures"
      from: any -> to: next_higher_level
```

---

## 4. Core Mechanism Deep Dive

### 4.1 Thinking Mode System

The four-level thinking system is defined in `THINKING_MODES.md` (456 lines) and
dispatched by `thinking_dispatch.md` (343 lines). Each level defines a token budget,
use cases, output style, and upgrade triggers.

```
┌──────────────────────────────────────────────────────────────────┐
│                  FOUR-LEVEL THINKING SYSTEM                     │
├──────────────────────────────────────────────────────────────────┤
│                                                                  │
│  QUICK ──────── Token budget: minimal                           │
│  │               Use: File I/O, state updates, env checks        │
│  │               Output: Concise, results only                   │
│  │               Never: Biological judgment, parameter selection  │
│  │                                                                │
│  ▼ (upgrade if: biological judgment needed)                      │
│  STANDARD ───── Token budget: normal                             │
│  │               Use: Known patterns, formula-based decisions     │
│  │               Output: Decision + brief rationale               │
│  │               Upgrade if: >30% deviation from expected         │
│  │                                                                │
│  ▼ (upgrade if: complex multi-metric judgment)                   │
│  DEEP ───────── Token budget: extended                           │
│  │               Use: Quality assessment, batch effects, errors   │
│  │               Output: Reasoning process, factors listed        │
│  │               Upgrade if: core biological interpretation       │
│  │                                                                │
│  ▼ (upgrade if: cell annotation, core biology)                   │
│  ULTRATHINK ─── Token budget: maximum                            │
│                  Use: Cell annotation, trajectory, paper writing  │
│                  Output: Full 6-step reasoning structure:         │
│                    1. Information gathering                       │
│                    2. Hypothesis generation                       │
│                    3. Validation and refutation                   │
│                    4. Conclusion and confidence                   │
│                    5. Uncertainty statement                       │
│                    6. Alternative considerations                 │
│                                                                  │
└──────────────────────────────────────────────────────────────────┘
```

**Ultrathink mandatory output fields**:

```yaml
# Every Ultrathink decision must produce:
mandatory_output:
  reasoning_process: "Complete reasoning chain"
  conclusion: "Clear actionable conclusion"
  confidence: 0.85          # Float 0-1
  uncertainty: "Main uncertainties and their impact"
  alternatives: "Other possibilities that were considered and why rejected"
```

---

### 4.2 Standard Pipeline Materialization

The framework includes 6 standard analysis pipelines that serve as **soft guidance** ---
they define expected steps and detection signals but do not rigidly enforce execution order.

```
┌────────────────────────────────────────────────────────────────────┐
│                   6 STANDARD PIPELINES                            │
├────────────────────────────────────────────────────────────────────┤
│                                                                    │
│  scrna_seq.yml ───────── Single-Cell RNA-seq                      │
│  bulk_rnaseq.yml ─────── Bulk RNA-seq                             │
│  spatial_transcriptomics.yml ── Spatial Transcriptomics           │
│  proteomics.yml ──────── Proteomics (Mass Spectrometry)           │
│  metabolomics.yml ────── Metabolomics                             │
│  lipidomics.yml ──────── Lipidomics                               │
│  _schema.yml ─────────── Pipeline definition schema               │
│                                                                    │
└────────────────────────────────────────────────────────────────────┘
```

**scRNA-seq pipeline structure** (as an example):

```
┌────────────────────────────────────────────────────────────────────┐
│           scRNA-seq STANDARD PIPELINE                             │
│  Minimal: QC -> Normalize -> HVG -> Scale -> PCA -> Cluster      │
│           -> Annotate -> DEG                                      │
├────────────────────────────────────────────────────────────────────┤
│                                                                    │
│  [1] Quality Control ──────────────────── [REQUIRED]              │
│      ├── Cell-level filtering (nFeature, nCount, percent.mt)      │
│      ├── Doublet removal (DoubletFinder / scDblFinder)            │
│      └── Gene-level filtering                                     │
│      Detection: nFeature_RNA, percent.mt, subset(, pp.filter_*    │
│                                                                    │
│  [2] Normalization ────────────────────── [REQUIRED]              │
│      ├── SCTransform or NormalizeData + LogNormalize               │
│      └── scran::computeSumFactors (alternative)                   │
│      Detection: SCTransform, NormalizeData, pp.normalize_total     │
│                                                                    │
│  [3] HVG Selection ───────────────────── [REQUIRED]               │
│      └── FindVariableFeatures / pp.highly_variable_genes           │
│      Detection: FindVariableFeatures, VariableFeatures, vst       │
│                                                                    │
│  [4] Data Scaling ─────────────────────── [REQUIRED]              │
│      └── ScaleData / pp.scale / pp.regress_out                    │
│      Detection: ScaleData, pp.scale                                │
│                                                                    │
│  [5] Dimensionality Reduction ─────────── [REQUIRED]              │
│      ├── PCA on variable features                                  │
│      ├── Determine significant PCs (ElbowPlot)                    │
│      └── UMAP or tSNE embedding                                   │
│      Detection: RunPCA, RunUMAP, RunTSNE, pp.pca, tl.umap        │
│                                                                    │
│  [6] Cell Clustering ─────────────────── [REQUIRED]               │
│      ├── Build nearest-neighbor graph                              │
│      ├── Clustering (Leiden or Louvain)                            │
│      └── Evaluate cluster resolution                               │
│      Detection: FindNeighbors, FindClusters, tl.leiden             │
│                                                                    │
│  [7] Cell Type Annotation ────────────── [REQUIRED]               │
│      ├── Identify marker genes per cluster                         │
│      └── Annotate (manual or automated: SingleR, scType, etc.)    │
│      Detection: FindAllMarkers, SingleR, scType, CellTypist       │
│                                                                    │
│  [8] Differential Expression ─────────── [REQUIRED]               │
│      └── DEG analysis between conditions/cell types                │
│      Detection: FindMarkers, DESeq2, MAST                         │
│                                                                    │
└────────────────────────────────────────────────────────────────────┘
```

**Pipeline compliance** is tracked via `pipeline_compliance.md`. It uses **detection
signals** (code patterns, metric patterns, output patterns) to automatically determine
which pipeline steps have been completed, producing a `compliance_score` (0-1 float).
This is a **soft** validation --- it warns but does not block if optional steps are skipped.

**Journal compliance** is configured per journal. Each journal file defines formatting
requirements, figure specifications, word limits, and data availability requirements:

| Journal | Key Requirements |
|---------|-----------------|
| Nature | Figure limit: 6-8, Extended Data permitted, strict word limits |
| Science | Brief main text, most data in Supplementary Materials |
| Cell Reports | STAR Methods format, Lead Contact required |
| JCI | Clinical focus, CONSORT/STROBE compliance for clinical studies |
| Nature Methods | Methods-focused, detailed protocol requirements |
| PLOS ONE | Open data mandate, less restrictive formatting |

---

### 4.3 Reproducibility System

Reproducibility is enforced through three mechanisms that work together to ensure any
analysis can be exactly reproduced.

```
┌────────────────────────────────────────────────────────────────────┐
│                  REPRODUCIBILITY SYSTEM                           │
├────────────────────────────────────────────────────────────────────┤
│                                                                    │
│  [1] Seed Management                                              │
│      ├── default_rules.yml defines: default_seed: 42              │
│      ├── Seed is recorded in every phase_summary                  │
│      └── Config resolution: Project > User > Plugin seed          │
│                                                                    │
│  [2] Parameter Recording                                          │
│      ├── Every phase_summary records parameters_used:             │
│      │   param_name:                                               │
│      │     value: <actual value used>                              │
│      │     default: <framework default>                            │
│      │     reason: "why this value was chosen"                     │
│      └── Step 5a compiles methods_record.md from all summaries    │
│                                                                    │
│  [3] Environment Capture (Step 5b)                                │
│      ├── Conda environment export (conda env export)               │
│      ├── R sessionInfo() capture                                   │
│      ├── Python package versions (pip freeze)                     │
│      └── Data provenance chain (input files + checksums)          │
│                                                                    │
└────────────────────────────────────────────────────────────────────┘
```

---

### 4.4 Error Memory & Auto-Recovery

The error handling system operates in three layers, with priority resolution
similar to the configuration system.

```
┌────────────────────────────────────────────────────────────────────┐
│              THREE-LAYER ERROR RESOLUTION                         │
├────────────────────────────────────────────────────────────────────┤
│                                                                    │
│  [LAYER 1] Framework Knowledge Base (COMMON_ERRORS.yml, 20 KB)   │
│  │  ├── Pre-defined error patterns with solutions                │
│  │  ├── Categories: R, Python, Seurat, Scanpy, Bioconductor, ... │
│  │  ├── Each entry: pattern + cause + solution + severity         │
│  │  └── Searched first on every error encounter                  │
│  │                                                                │
│  ▼ (if not found in KB)                                          │
│  [LAYER 2] User Experience Library (learned_solutions.yml)           │
│  │  ├── Project-local error history                               │
│  │  ├── Auto-populated when AI resolves novel errors              │
│  │  └── Format: pattern + context + resolution + timestamp        │
│  │                                                                │
│  ▼ (if not found in experience)                                  │
│  [LAYER 3] AI Reasoning                                          │
│     ├── Thinking mode auto-escalates to Deep or Ultrathink       │
│     ├── Uses general programming + biology knowledge             │
│     └── If resolved: auto-saves to experience library            │
│                                                                    │
└────────────────────────────────────────────────────────────────────┘
```

**Experience library auto-population example**:

```yaml
# learned_solutions.yml - auto-populated entry
- id: "ERR-20260301-001"
  timestamp: "2026-03-01T14:23:00Z"
  error_pattern: "Error in CreateSeuratObject: No cell names"
  context:
    phase: "phase01_scrna"
    tool: "Seurat::CreateSeuratObject"
    r_version: "4.4.0"
  resolution: |
    Issue: Input matrix had no column names (cell barcodes).
    Fix: Added colnames(matrix) <- paste0("Cell_", seq_len(ncol(matrix)))
    before calling CreateSeuratObject.
  thinking_mode_used: "deep"
  reusable: true
```

---

### 4.5 Decision Node System

Critical analysis decisions are governed by a graduated response system defined
in `decision_nodes.md`. Each decision node defines thresholds for different
response levels.

```yaml
# decision_nodes structure example
decision_nodes:
  overclustering_check:
    description: "Evaluate whether clustering resolution is too high"
    metrics:
      - name: "min_cluster_size"
        thresholds:
          ok: ">= 50 cells"
          soft_stop: "10-50 cells"
          hard_stop: "< 10 cells"
      - name: "silhouette_score"
        thresholds:
          ok: "> 0.2"
          soft_stop: "0.05-0.2"
          hard_stop: "< 0.05"
    responses:
      ok: "Proceed normally"
      soft_stop: "Log warning, consider resolution adjustment"
      hard_stop: "PAUSE execution, report to user with evidence"

  differential_expression_validity:
    description: "Check DEG result quality"
    metrics:
      - name: "n_significant_genes"
        thresholds:
          ok: "> 20"
          soft_stop: "5-20"
          hard_stop: "< 5 or > 5000"
    responses:
      hard_stop: "Possible technical artifact or wrong comparison"
```

**Graduated response levels**:

| Level | Symbol | Behavior |
|-------|--------|----------|
| OK | -- | Proceed normally, no logging |
| SOFT_STOP | [WARN] | Log warning in findings.yml, continue with caution |
| HARD_STOP | [STOP] | Halt execution, display evidence to user, wait for decision |

---

### 4.6 Adaptive Capability Technical Foundation

The framework's adaptive capabilities span three dimensions:

```
┌────────────────────────────────────────────────────────────────────┐
│                 ADAPTIVE CAPABILITIES                             │
├────────────────────────────────────────────────────────────────────┤
│                                                                    │
│  [1] Dynamic Phase Adjustment                                     │
│      ├── Step 0 generates phase_index.yml based on TOPIC.yml      │
│      │   (phase count and content are topic-dependent)             │
│      ├── /bio phase add: Add phases during execution               │
│      ├── /bio phase skip: Skip phases that become unnecessary     │
│      ├── /bio phase modify: Adjust phase parameters mid-run       │
│      └── Step 2 can trigger supplementary analysis phases         │
│                                                                    │
│  [2] Parameter Adaptation                                         │
│      ├── Three-layer config resolution adapts to project context  │
│      ├── AI adjusts parameters based on data characteristics      │
│      │   (e.g., resolution based on cell count)                    │
│      ├── upgrade_triggers auto-escalate thinking depth            │
│      └── Anomaly triggers adjust analysis strategy                │
│                                                                    │
│  [3] Experience Accumulation                                      │
│      ├── learned_solutions.yml: Error solutions persist across phases │
│      ├── findings.yml: Anomalies tracked with resolution status   │
│      ├── phase_summaries.yml: Decisions + rationale for future ref│
│      └── /bio experience search: Query accumulated knowledge      │
│                                                                    │
└────────────────────────────────────────────────────────────────────┘
```

**Key design point**: Phase count is **never hardcoded**. A simple scRNA-seq project
might have 4 phases; a multi-omics integration study might have 12. The framework
dynamically generates `phase_index.yml` at Step 0 based on the research questions
and data types defined in TOPIC.yml.

---

### 4.7 Visual Review Closed Loop [New in v1.0.0]

Previous versions had no mechanism to verify that generated figures were visually
correct. The AI would check metadata (file size, DPI, resolution) but never actually
look at the figures. v1.0.0 introduces a two-stage visual review system.

```
┌────────────────────────────────────────────────────────────────────┐
│              VISUAL REVIEW CLOSED LOOP (v1.0.0)                  │
├────────────────────────────────────────────────────────────────────┤
│                                                                    │
│  Stage 1: AI Visual Inspection (Step 4d)                          │
│  ├── For each generated figure:                                    │
│  │   1. Read PNG preview via Read tool (multimodal)                │
│  │   2. Check: labels readable? colors distinguishable?            │
│  │   3. Check: axes correct? legend present?                       │
│  │   4. Check: scientific content matches figure_plan.yml spec?    │
│  │   5. If issues found: auto-regenerate with corrections          │
│  └── Output: figure_checklist.md with per-figure pass/fail        │
│                                                                    │
│  Stage 2: User Visual Review (Step 4f)                            │
│  ├── Display all figures to user                                   │
│  ├── HARD_STOP: Wait for user approval                            │
│  ├── User options:                                                 │
│  │   ├── Approve all ──► Proceed to Step 5                        │
│  │   ├── Request revision for specific figures                    │
│  │   │   └── AI regenerates ──► Re-display for approval           │
│  │   └── Major concerns ──► Return to Step 4a for rework          │
│  └── visual_review_status tracks: pending/in_review/approved/     │
│      revision_in_progress                                          │
│                                                                    │
│  WHY THIS MATTERS:                                                │
│  - "File exists and is 300 DPI" does not mean the figure is good  │
│  - A UMAP might have overlapping labels or wrong color scheme     │
│  - A heatmap might have unreadable gene names                     │
│  - Only visual inspection can catch these issues                  │
│                                                                    │
└────────────────────────────────────────────────────────────────────┘
```

---

### 4.8 Data Quality Deep Audit [New in v1.0.0]

Step 1.5 fills a critical gap between "analysis execution" and "result interpretation."
Phase-level QC catches per-phase issues, but cross-phase problems (sample swaps, severe
batch effects affecting downstream interpretation) require a holistic review.

```
┌────────────────────────────────────────────────────────────────────┐
│          STEP 1.5: DATA QUALITY DEEP AUDIT                       │
│          (1,840 lines of specification)                           │
├────────────────────────────────────────────────────────────────────┤
│                                                                    │
│  TRIGGER CONDITIONS:                                              │
│  ├── User explicit: TOPIC.yml quality_review.enabled = true/false │
│  └── Auto mode (default):                                         │
│      ├── Research type: clinical / method_comparison ──► trigger  │
│      ├── Quality requirements: high ──► trigger                   │
│      └── Accumulated warnings >= 3 ──► trigger                    │
│      (Routine exploratory analyses: skip, zero overhead)          │
│                                                                    │
│  5 AUDIT DIMENSIONS:                                              │
│                                                                    │
│  D1. Sample Identity Verification                                 │
│      ├── Cross-reference sample metadata with analysis results    │
│      ├── Check for potential sample swaps                         │
│      └── Verify sample grouping matches experimental design       │
│                                                                    │
│  D2. Batch Effect Assessment                                      │
│      ├── Evaluate residual batch effects post-correction          │
│      ├── Check if batch dominates biological variation            │
│      └── LISI/kBET score evaluation                               │
│                                                                    │
│  D3. Biological Plausibility                                      │
│      ├── Are cell type proportions within expected ranges?        │
│      ├── Do marker gene patterns match literature?                │
│      └── Are DE results biologically coherent?                    │
│                                                                    │
│  D4. Outlier Sample Detection                                     │
│      ├── PCA/UMAP-based outlier detection                         │
│      ├── Library complexity comparison across samples             │
│      └── Gene detection rate consistency                          │
│                                                                    │
│  D5. Cross-Phase Consistency                                      │
│      ├── Do results from different phases agree?                  │
│      ├── Are there contradictory findings?                        │
│      └── Data provenance chain integrity                          │
│                                                                    │
│  OUTPUT: quality_report.yml                                       │
│  ├── per_dimension: { status: PASS|WARNING|CRITICAL, details }   │
│  ├── overall_decision: PASS|WARNING|CRITICAL                     │
│  └── recommendations: [ list of suggested actions ]              │
│                                                                    │
│  DECISION MATRIX:                                                 │
│  ├── All PASS ──────────────────► Auto-proceed to Step 2         │
│  ├── Any WARNING, no CRITICAL ──► Proceed with caution (5s delay)│
│  └── Any CRITICAL ──────────────► HARD_STOP (user must decide)   │
│                                                                    │
└────────────────────────────────────────────────────────────────────┘
```

---

## 5. ORCHESTRATOR Deep Dive

The ORCHESTRATOR is the central command hub of Bio-Framework. It receives all input
(commands and natural language), routes it to the appropriate handler, and manages
the execution lifecycle.

```
┌────────────────────────────────────────────────────────────────────┐
│                 ORCHESTRATOR COMMAND CENTER                       │
├────────────────────────────────────────────────────────────────────┤
│                                                                    │
│  ┌──────────────┐                                                 │
│  │  User Input  │                                                 │
│  └──────┬───────┘                                                 │
│         ▼                                                          │
│  ┌──────────────────────────────────────────┐                     │
│  │        ORCHESTRATOR (398 lines)          │                     │
│  │                                          │                     │
│  │  1. Input Classification                 │                     │
│  │     ├── /bio command? ─► COMMAND_PARSER  │                     │
│  │     ├── Natural language? ─► Parse intent│                     │
│  │     └── Continuation? ─► Resume workflow │                     │
│  │                                          │                     │
│  │  2. Layer Resolution                     │                     │
│  │     Project > User Global > Plugin       │                     │
│  │                                          │                     │
│  │  3. Context Management                   │                     │
│  │     ├── Init budget: ~2700 tokens        │                     │
│  │     ├── Per-phase: ~5800 tokens          │                     │
│  │     └── Auto /compact every 3 phases     │                     │
│  └────────────┬─────────────────────────────┘                     │
│               │                                                    │
│       ┌───────┼───────────────────────────────────┐               │
│       ▼       ▼       ▼       ▼       ▼           ▼               │
│  ┌────────┐ ┌──────┐ ┌──────┐ ┌──────┐ ┌──────┐ ┌──────┐        │
│  │workflow│ │topic │ │config│ │phase │ │help  │ │debug │        │
│  │control │ │load  │ │mgr   │ │mgmt  │ │system│ │tools │        │
│  └───┬────┘ └──────┘ └──────┘ └──────┘ └──────┘ └──────┘        │
│      │                                                             │
│      ▼                                                             │
│  ┌───────────────────────────────────────────────────┐            │
│  │              EXECUTION ENGINE                     │            │
│  │                                                    │            │
│  │  Step 0 ──► Phases ──► Step 1.5 ──► Step 2 ──►   │            │
│  │  Step 3 ──► Step 4 ──► Step 5 ──► Step 6          │            │
│  │                                                    │            │
│  │  Sub-modules loaded on demand:                    │            │
│  │  execution_steps / data_quality_review /           │            │
│  │  reflection_exploration / research_summary /       │            │
│  │  publication_figures / manuscript_drafting /        │            │
│  │  final_report / thinking_dispatch /                │            │
│  │  pipeline_compliance / checkpoint_recovery         │            │
│  └───────────────────────────────────────────────────┘            │
│                                                                    │
└────────────────────────────────────────────────────────────────────┘
```

### Core Responsibilities

| Responsibility | Description | Importance |
|----------------|-------------|------------|
| Input routing | Classify user input and route to correct handler | [HIGH] Must never misroute |
| Auto-execution | Execute without confirmation unless HARD_STOP | [HIGH] Core UX principle |
| Layer resolution | Resolve resources through 3-layer priority | [HIGH] Ensures correct configs |
| Context management | Budget tokens, auto-compact, selective loading | [HIGH] Prevents context overflow |
| Phase lifecycle | Trigger phase completion handler, update state | [HIGH] State consistency |
| Sub-module dispatch | Load sub-modules on demand for each step | [MED] Performance optimization |
| Error escalation | Route errors to KB, then experience, then AI | [MED] Resilience |

### SKILL.md -> ORCHESTRATOR -> Sub-Module Dispatch

```
SKILL.md (551 lines)
│
├── Section 1: Role Definition
│   └── Capabilities, supported domains, core principles
│
├── Section 2: Command System
│   ├── 2.1 Command Recognition (prefix aliases)
│   ├── 2.2 Command Routing Table (25+ commands)
│   └── 2.3 Command Parsing Algorithm
│
├── Section 3-7: Behavioral Rules
│   ├── Three-layer resolution
│   ├── Thinking mode integration
│   ├── Forbidden behaviors
│   └── Session lifecycle (auto-detect, auto-save)
│
├── Section 8: Sub-Module Index
│   └── Points to ORCHESTRATOR.md and all sub-modules
│
└── On any /bio command ──► ORCHESTRATOR.md
    │
    ├── Workflow commands ──► WORKFLOW_CONTROLLER.md
    │   └── start/continue/stop/status/debug
    │
    ├── Topic commands ──► TOPIC_LOADING.md
    │
    ├── Config commands ──► CONFIG_MANAGER.md
    │
    ├── Help commands ──► HELP_SYSTEM.md
    │
    └── Workflow execution ──► Sub-modules (on demand)
        ├── Step 0-1 ──► execution_steps.md + thinking_dispatch.md
        ├── Step 1.5 ──► data_quality_review.md
        ├── Step 2   ──► reflection_exploration.md
        ├── Step 3   ──► research_summary.md
        ├── Step 4   ──► publication_figures.md
        ├── Step 5   ──► final_report.md
        ├── Step 6   ──► manuscript_drafting.md
        ├── Recovery ──► checkpoint_recovery.md
        └── Pipeline ──► pipeline_compliance.md
```

---

## 6. claude-router Integration

Bio-Framework is designed to work with `claude-router`, a model-selection layer that
routes requests to the optimal Claude model (Opus/Sonnet/Haiku) based on task
complexity. Integration is achieved through the `thinking_mode` field in
`workflow_state.yml`.

```
┌────────────────────────────────────────────────────────────────────┐
│              CLAUDE-ROUTER INTEGRATION                            │
├────────────────────────────────────────────────────────────────────┤
│                                                                    │
│  Bio-Framework                        claude-router                │
│  ┌──────────────────┐                ┌──────────────────┐         │
│  │ thinking_dispatch │                │ Model Selector   │         │
│  │                   │    writes      │                  │         │
│  │ Determine mode    │──────────────►│ Read             │         │
│  │ per task type     │  thinking_mode │ workflow_state   │         │
│  │                   │  to state file │                  │         │
│  └──────────────────┘                │ Map to model:    │         │
│                                       │  ultra → Opus    │         │
│                                       │  deep  → Sonnet  │         │
│                                       │  std   → Sonnet  │         │
│                                       │  quick → Haiku   │         │
│                                       └──────────────────┘         │
│                                                                    │
│  WORKFLOW:                                                         │
│  1. Bio-Framework skill determines thinking_mode for current task  │
│  2. Updates workflow_state.yml with thinking_mode field             │
│  3. claude-router reads workflow_state.yml before each API call    │
│  4. Selects model based on thinking_mode mapping                   │
│  5. Executes request with selected model                           │
│                                                                    │
└────────────────────────────────────────────────────────────────────┘
```

### Thinking Mode to Model Mapping

| Thinking Mode | Model | Typical Tasks |
|---------------|-------|---------------|
| `ultrathink` | Opus | Cell annotation, manuscript writing, topic understanding |
| `deep` | Sonnet | Quality assessment, batch evaluation, figure quality check |
| `standard` | Sonnet | Known patterns, routine code generation, parameter selection |
| `quick` | Haiku | File I/O, state updates, environment checks |

### Cost Optimization Impact

| Scenario | Without Router | With Router | Savings |
|----------|---------------|-------------|---------|
| Simple file operations | Sonnet | Haiku | ~80% |
| Standard analysis code | Opus | Sonnet | ~40% |
| Cell type annotation | Sonnet | Opus | Quality improvement |
| Manuscript writing | Random | Opus | Quality guaranteed |
| **Estimated overall** | -- | -- | **Significant cost reduction (varies by workload)** |

The key insight is that most tokens in a bioinformatics session are spent on file
operations and standard code execution (Quick/Standard modes), while only a small
fraction requires Opus-level reasoning. By matching model capability to task
complexity, significant cost savings are achieved without sacrificing quality where
it matters most.

---

## 7. Command System

### Command Routing Architecture

```
┌────────────────────────────────────────────────────────────────────┐
│                  COMMAND ROUTING ARCHITECTURE                     │
├────────────────────────────────────────────────────────────────────┤
│                                                                    │
│  User Input                                                        │
│  ┌────────────────────────────────────┐                           │
│  │ "/bio start --quick"               │                           │
│  │ "/生信 开始"                        │                           │
│  │ "/バイオ start"                     │                           │
│  │ "start the analysis"               │                           │
│  └──────────────┬─────────────────────┘                           │
│                 ▼                                                   │
│  ┌────────────────────────────────────┐                           │
│  │ COMMAND_PARSER.md (556 lines)      │                           │
│  │                                     │                           │
│  │ 1. Extract: command + options       │                           │
│  │ 2. Alias resolution (manifest.yml)  │                           │
│  │ 3. Subcommand parsing               │                           │
│  │ 4. Parameter validation             │                           │
│  │ 5. Natural language fallback        │                           │
│  └──────────────┬─────────────────────┘                           │
│                 ▼                                                   │
│  ┌────────────────────────────────────┐                           │
│  │ SKILL.md Routing Table             │                           │
│  │                                     │                           │
│  │ start/continue/stop/status/debug    │                           │
│  │   └──► WORKFLOW_CONTROLLER.md       │                           │
│  │ help                                │                           │
│  │   └──► HELP_SYSTEM.md               │                           │
│  │ topic new/load/show/edit            │                           │
│  │   └──► TOPIC_LOADING.md             │                           │
│  │ config show/set/reset               │                           │
│  │   └──► CONFIG_MANAGER.md            │                           │
│  │ phase list/show/rerun/add/modify    │                           │
│  │   └──► Phase Management             │                           │
│  │ pipeline show/adjust/add-step       │                           │
│  │   └──► Pipeline Management          │                           │
│  │ finding add/list/show/resolve       │                           │
│  │   └──► Finding Tracker              │                           │
│  │ experience list/search              │                           │
│  │   └──► Experience Store             │                           │
│  │ manuscript check/sync               │                           │
│  │   └──► manuscript_checker.md        │                           │
│  │ prime                               │                           │
│  │   └──► Context Snapshot             │                           │
│  │ version                             │                           │
│  │   └──► Version Info                 │                           │
│  └────────────────────────────────────┘                           │
│                                                                    │
└────────────────────────────────────────────────────────────────────┘
```

### 25+ Command Quick Reference

| Command | Description | Key Options |
|---------|-------------|-------------|
| `start` | Begin new analysis workflow | `--topic`, `--from`, `--skip`, `--quick` |
| `continue` | Resume from checkpoint | `--from` |
| `stop` | Pause and save state | -- |
| `status` | Display progress and state | -- |
| `help` | Show help information | `<topic>` |
| `topic new` | Create new research topic | -- |
| `topic load` | Load existing TOPIC.yml | `<file>` |
| `topic show` | Display current topic | -- |
| `topic edit` | Edit topic definition | -- |
| `phase list` | List all phases | -- |
| `phase show` | Show phase details | `<N>` |
| `phase rerun` | Re-execute a phase | `<N>` |
| `phase add` | Add new phase | -- |
| `phase modify` | Modify phase config | -- |
| `phase skip` | Skip a phase | `<N>` |
| `config show` | Show configuration | `<key>` |
| `config set` | Set configuration value | `<key> <value>` |
| `config reset` | Reset to defaults | `<key>` |
| `pipeline show` | Show pipeline compliance | -- |
| `pipeline adjust` | Adjust pipeline | -- |
| `finding add` | Record anomaly/finding | `<description>` |
| `finding list` | List findings | -- |
| `experience list` | List experience entries | -- |
| `experience search` | Search experience library | `<query>` |
| `prime` | Generate context snapshot | (~800 tokens) |
| `manuscript check` | Check manuscript quality | -- |
| `manuscript sync` | Sync manuscript with data | -- |
| `debug why` | Explain tool choice | `<tool_name>` |
| `debug why-not` | Explain tool exclusion | `<tool_name>` |
| `debug on/off` | Toggle debug logging | -- |
| `debug trace` | Show execution trace | -- |
| `debug report` | Generate debug report | -- |
| `version` | Show framework version | -- |

### Multilingual Alias System

The command system supports three languages with 110+ total aliases:

```yaml
# Prefix aliases (all equivalent to /bio):
prefixes:
  en: ["/bio", "/bioinformatics"]
  zh: ["/生信"]
  ja: ["/バイオ"]

# Command alias examples:
start:
  en: ["run", "begin", "execute"]
  zh: ["开始", "启动", "运行"]

continue:
  en: ["resume"]
  zh: ["继续", "恢复"]

status:
  en: ["info", "progress"]
  zh: ["状态", "进度"]

help:
  en: ["h"]
  zh: ["帮助"]

finding:
  en: ["find", "discovery", "note"]
  zh: ["发现", "记录", "异常"]

# Natural language fallback (no prefix required):
# "start analysis" -> /bio start
# "show me the status" -> /bio status
# "what is the progress" -> /bio status
```

---

## 8. State Management

### current_step Enum (10 values)

The `current_step` field in `workflow_state.yml` tracks exactly where the workflow
is positioned. This enum is used by checkpoint recovery, the status display, and
the thinking dispatch system.

```
┌────────────────────────────────────────────────────────────────────┐
│                  current_step STATE MACHINE                      │
├────────────────────────────────────────────────────────────────────┤
│                                                                    │
│  init ──────────────────────► phase_running                       │
│  (Step 0 complete)             (Step 1 active)                    │
│                                     │                             │
│                                     ▼                             │
│                              analysis_complete                    │
│                              (all phases done)                    │
│                                     │                             │
│                          ┌──────────┴──────────┐                 │
│                          ▼                      ▼                 │
│                 step1_5_quality_review    (skip to step2)         │
│                 (Step 1.5 active)               │                 │
│                          │                      │                 │
│                          └──────────┬───────────┘                 │
│                                     ▼                             │
│                              step2_reflection                     │
│                              (Step 2 active)                      │
│                                     │                             │
│                                     ▼                             │
│                              step3_summary                        │
│                              (Step 3 active)                      │
│                                     │                             │
│                                     ▼                             │
│                              step4_figures                        │
│                              (Step 4 active)                      │
│                                     │                             │
│                                     ▼                             │
│                              step5_report                         │
│                              (Step 5 active)                      │
│                                     │                             │
│                                     ▼                             │
│                              step6_manuscript                     │
│                              (Step 6 active)                      │
│                                     │                             │
│                                     ▼                             │
│                                completed                          │
│                                                                    │
└────────────────────────────────────────────────────────────────────┘
```

### workflow_state.yml Complete Structure

```yaml
# project_skills/runtime/workflow_state.yml
workflow_state:
  # Core state
  current_phase: "phase01_scrna"              # Active or last phase ID
  current_step: "phase_running"               # Enum (10 values, see above)
  completed_phases: ["phase01_scrna"]          # List of completed phase IDs
  interrupted: false                           # true if session ended mid-step

  # Step 1.5 fields [NEW v1.0.0]
  step_1_5_triggered: false                    # Was Step 1.5 executed?
  quality_review_status: null                  # null|in_progress|passed|
                                               #   passed_with_warnings|failed
  critical_issues: []                          # CRITICAL issue IDs from audit

  # Step 2 exploration state
  exploration_round: 0                         # 0=not started, 1=first, 2=second
  exploration_status: null                     # null|gap_analysis|gate|executing|
                                               #   re_evaluation|completed|skipped
  gate_decision: null                          # null|explore|skip_all

  # Step 4 visual review [NEW v1.0.0]
  visual_review_status: null                   # null|pending|in_review|approved|
                                               #   revision_in_progress

  # Pipeline compliance
  pipeline_compliance:
    mandatory_pending: []                      # Steps not yet completed
    compliance_score: null                     # Float 0-1

  # Interruption tracking
  interruption_tracking:
    last_checkpoint: "2026-03-01T14:00:00Z"   # ISO timestamp
    stop_reason: null                          # null or HARD_STOP reason
    resume_instructions: null                  # What user should do

  # Resources
  resources:
    virtual_env: "~/conda/envs/bioanalysis"
    r_version: "4.4.0"
    python_version: "3.10.0"

  # Lifecycle
  status: "running"                            # pending|running|paused|completed|failed
  topic_file: "project_skills/TOPIC.yml"
  started_at: "2026-03-01T10:00:00Z"
  updated_at: "2026-03-01T14:00:00Z"
```

### Checkpoint Recovery Strategy

The framework supports four types of recovery when a workflow is interrupted:

| Type | Trigger | Thinking Mode | Action |
|------|---------|---------------|--------|
| **Type A: Continue** | User interruption (session end) | Standard | Restore to exact interruption point, continue |
| **Type B: Rerun** | Error during execution | Deep | Backup output, reset state, re-execute phase |
| **Type C: Start From** | User wants to skip phases | Deep | Verify predecessor outputs exist, start from specified phase |
| **Type D: Restart** | User wants fresh start | -- | Requires confirmation, resets everything |

Recovery flow:

```
/bio continue
    │
    ├── Read workflow_state.yml
    ├── Read phase_summaries.yml
    ├── Check partial output files
    │
    ├── Diagnose interruption type:
    │   ├── Error? ──► Type B (ULTRATHINK for root cause analysis)
    │   ├── User? ──► Type A (STANDARD, continue directly)
    │   └── Timeout? ──► Deep analysis (complexity/data volume)
    │
    └── Execute recovery:
        ├── Activate virtual environment
        ├── Load completed phase objects
        ├── Restore to interruption point
        └── Continue execution
```

---

## 9. Data Security Architecture

Bio-Framework is designed for **fully local execution**. No experimental data ever
leaves the user's machine. The AI operates on data through code execution, not data
upload.

```
┌────────────────────────────────────────────────────────────────────┐
│                  DATA SECURITY ARCHITECTURE                       │
├────────────────────────────────────────────────────────────────────┤
│                                                                    │
│  ┌──────────────────────────────────────────┐                     │
│  │         USER'S LOCAL MACHINE              │                     │
│  │                                           │                     │
│  │  ┌─────────────────┐   ┌──────────────┐  │                     │
│  │  │ Experimental    │   │ Bio-Framework│  │                     │
│  │  │ Data            │   │ (Skills)     │  │                     │
│  │  │                 │   │              │  │                     │
│  │  │ - .h5ad files   │   │ - SKILL.md   │  │                     │
│  │  │ - .qs/.rds      │◄──│ - ORCH.md    │  │                     │
│  │  │ - .csv/.tsv     │   │ - Knowledge  │  │                     │
│  │  │ - FASTQ         │   │   base       │  │                     │
│  │  │                 │   │              │  │                     │
│  │  └────────┬────────┘   └──────┬───────┘  │                     │
│  │           │                    │          │                     │
│  │           ▼                    ▼          │                     │
│  │  ┌───────────────────────────────────┐   │                     │
│  │  │       Claude Code Runtime         │   │                     │
│  │  │                                    │   │                     │
│  │  │  AI generates R/Python code       │   │                     │
│  │  │  Code executes LOCALLY            │   │                     │
│  │  │  Results stay on machine          │   │                     │
│  │  │  Only code + text go to Claude API│   │                     │
│  │  └───────────────────────────────────┘   │                     │
│  │                                           │                     │
│  └──────────────────────────────────────────┘                     │
│                      │                                             │
│                      │ Code + prompts only                        │
│                      │ (NO experimental data)                     │
│                      ▼                                             │
│  ┌──────────────────────────────────────────┐                     │
│  │           Claude API (Cloud)             │                     │
│  │                                           │                     │
│  │  Receives: Code text, error messages,    │                     │
│  │            summary statistics             │                     │
│  │  Never receives: Raw expression matrices,│                     │
│  │            patient data, FASTQ files      │                     │
│  └──────────────────────────────────────────┘                     │
│                                                                    │
│  WHAT GOES TO CLOUD:           WHAT STAYS LOCAL:                  │
│  ├── R/Python code text        ├── Expression matrices (.h5ad)    │
│  ├── Error messages            ├── Count matrices (.csv/.tsv)     │
│  ├── Summary statistics        ├── Seurat objects (.qs/.rds)      │
│  ├── Plot descriptions         ├── Raw sequencing data (FASTQ)    │
│  ├── Phase summaries           ├── Generated figures (PNG/PDF)    │
│  └── Workflow state            ├── Patient/sample metadata        │
│                                └── All intermediate files          │
│                                                                    │
└────────────────────────────────────────────────────────────────────┘
```

**Key security guarantees**:

1. **No data exfiltration**: The framework generates code that runs locally. Raw data
   files are read by R/Python on the local machine, not uploaded to any API.

2. **Metadata minimization**: Phase summaries contain only summary statistics
   (e.g., "15,234 cells after QC"), never individual cell/sample data.

3. **User control**: All generated code is visible to the user before and during
   execution. The user can inspect, modify, or reject any operation.

4. **No persistent cloud storage**: Claude API calls are stateless. No analysis data
   persists on Anthropic servers beyond the conversation context.

---

## 10. v1.0.0 Major Evolution from Earlier Versions

Bio-Framework v1.0.0 represents a significant architectural evolution from earlier
internal development versions. The following table summarizes key differences.

| Dimension | Pre-v1.0 (Internal) | v1.0.0 | Impact |
|-----------|---------------------|--------|--------|
| **Quality gate** | No systematic data quality audit | Step 1.5: 5-dimension deep audit (conditional) | Catches sample swaps, batch effects before interpretation |
| **Visual review** | Metadata-only figure validation (DPI, size) | AI visual inspection + user approval gate (Step 4d/4f) | Eliminates unreadable labels, wrong colors, misaligned panels |
| **Figure planning** | Figures generated ad-hoc during analysis | Step 3: Dedicated figure planning with figure_types.yml cross-reference | Ensures mandatory figures are never missed, literature-aligned specs |
| **Mandatory figures** | No enforcement of required figure types | mandatory_type_id_check with HARD_STOP | scRNA must have UMAP, spatial must have spatial scatter, etc. |
| **Manuscript structure** | Monolithic manuscript generation | 21 independent sub-steps with per-step timeout resilience | Any sub-step can fail and retry without losing other sections |
| **Thinking system** | 3 levels (no Ultrathink) | 4 levels with 9 mandatory Ultrathink scenarios | Critical biological decisions get maximum reasoning depth |
| **Exploration loop** | Single-pass gap analysis | 2-round exploration with threshold decay (80% then 70%) | Diminishing returns respected, critical gaps force exploration |
| **Pipeline coverage** | 3 pipelines (scRNA, bulk, proteomics) | 6 standard pipeline files (+ ATAC-seq knowledge coverage via figure_types.yml) | Broader omics coverage with detection signals |
| **Journal compliance** | 2 journals | 6 journals (Nature, Science, Cell Reports, JCI, Nature Methods, PLOS ONE) | Manuscript formatting pre-validated per target journal |
| **Sub-module architecture** | Monolithic orchestrator | 13 specialized sub-modules (10,318 lines) loaded on demand | Better context management, reduced token usage |
| **Command system** | 10 commands, English only | 25+ commands, 110+ multilingual aliases (EN/ZH/JA) | Accessibility for non-English-speaking researchers |
| **State machine** | 5 states | 10 current_step values + 5 status values | Precise checkpoint recovery, accurate progress reporting |
| **Decision nodes** | Binary pass/fail | Graduated response (OK / SOFT_STOP / HARD_STOP) with quantified thresholds | Nuanced handling of borderline results |
| **Error recovery** | KB lookup only | Three-layer: KB -> experience library -> AI reasoning (auto-save) | Learning from project-specific errors |
| **Reproducibility** | Manual parameter recording | Automated seed management + methods_record.md + environment capture | Every analysis is reproducible by default |
| **Context budget** | Unbounded (context overflow risk) | Budgeted (~2700 init, ~5800/phase, auto-compact every 3 phases) | Prevents context window overflow in long analyses |

### Architecture Line Count Comparison

```
┌──────────────────────────────────────────────────────────┐
│           FRAMEWORK COMPLEXITY (v1.0.0)                  │
├──────────────────────────────────────────────────────────┤
│                                                          │
│  SKILL.md (entry point):               551 lines        │
│  ORCHESTRATOR.md (main):               398 lines        │
│  Top-level skills (7 files):         4,803 lines        │
│  Sub-modules (13 files):            10,318 lines        │
│  ─────────────────────────────────────────────           │
│  Total framework skills:            16,070 lines        │
│                                                          │
│  Knowledge base:                                         │
│  - COMMON_ERRORS.yml:                  20 KB             │
│  - figure_types.yml:                   68 KB             │
│  - sci_design_standards.yml:           36 KB             │
│  - Pipelines (6 files):            ~1,200 lines         │
│  - Journals (6 files):             ~1,000 lines         │
│  ─────────────────────────────────────────────           │
│  Total knowledge base:              ~350 KB              │
│                                                          │
└──────────────────────────────────────────────────────────┘
```

---

## Appendix: Key File Paths

| File | Path | Purpose |
|------|------|---------|
| SKILL.md | `framework/templates/skills/bio/SKILL.md` | Entry point, command routing |
| ORCHESTRATOR.md | `framework/skills/ORCHESTRATOR.md` | Main workflow orchestrator |
| WORKFLOW_CONTROLLER | `framework/skills/WORKFLOW_CONTROLLER.md` | Lifecycle management |
| COMMAND_PARSER | `framework/skills/COMMAND_PARSER.md` | Command parsing |
| HELP_SYSTEM | `framework/skills/HELP_SYSTEM.md` | Help generation |
| THINKING_MODES | `framework/skills/THINKING_MODES.md` | Thinking depth definitions |
| TOPIC_LOADING | `framework/skills/TOPIC_LOADING.md` | Topic YAML parsing |
| CONFIG_MANAGER | `framework/skills/CONFIG_MANAGER.md` | Three-layer config management |
| Sub-modules | `framework/skills/orchestrator/*.md` | 13 specialized sub-modules |
| Default rules | `framework/rules/default_rules.yml` | Framework default parameters |
| Manifest | `framework/manifest.yml` | Plugin metadata + command definitions |
| Pipelines | `framework/knowledge/pipelines/*.yml` | Standard analysis pipelines |
| Journals | `framework/knowledge/journals/*.yml` | Journal compliance configs |
| Error KB | `framework/knowledge/COMMON_ERRORS.yml` | Error pattern library |
| Figure types | `framework/knowledge/figure_types.yml` | Standard figure definitions |
| SCI standards | `framework/knowledge/sci_design_standards.yml` | SCI design standards |

---

*Bio-Framework v1.0.0 -- Skill Deep Dive*
*Generated: 2026-03-10*
