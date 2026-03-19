# Bio-Framework v1.0.0 — System Architecture and Feature Documentation

**Type:** Claude Code Skill Plugin
**Version:** 1.0.0 (Stable)
**Build Date:** 2026-02-25
**License:** [Public Evaluation License (Non-Commercial)](../LICENSE.md)

---

## Table of Contents

1. [System Architecture Overview](#1-system-architecture-overview)
   - 1.1 [Three-Layer Resource Architecture](#11-three-layer-resource-architecture)
   - 1.2 [Core Module Architecture](#12-core-module-architecture)
   - 1.3 [Knowledge Base System](#13-knowledge-base-system)
2. [Feature Description](#2-feature-description)
   - 2.1 [End-to-End Workflow Engine](#21-end-to-end-workflow-engine)
   - 2.2 [Four-Level Thinking Depth System](#22-four-level-thinking-depth-system)
   - 2.3 [Command System](#23-command-system)
   - 2.4 [Knowledge Accumulation and Error Handling](#24-knowledge-accumulation-and-error-handling)
   - 2.5 [Checkpoint Recovery and Session Management](#25-checkpoint-recovery-and-session-management)
3. [System Characteristics](#3-system-characteristics)
   - 3.1 [Adaptive Analysis Engine](#31-adaptive-analysis-engine)
   - 3.2 [Production-Grade Quality Assurance](#32-production-grade-quality-assurance)
   - 3.3 [Fully Local Execution](#33-fully-local-execution)
   - 3.4 [Modular and Extensible](#34-modular-and-extensible)
   - 3.5 [Full SCI Paper Workflow Support](#35-full-sci-paper-workflow-support)
   - 3.6 [Multilingual Internationalization](#36-multilingual-internationalization)
4. [Architecture Diagrams](#4-architecture-diagrams)
   - 4.1 [Data Flow Diagram](#41-data-flow-diagram)
   - 4.2 [Command Routing Diagram](#42-command-routing-diagram)
   - 4.3 [Resource Lookup Diagram](#43-resource-lookup-diagram)

---

## 1. System Architecture Overview

### 1.1 Three-Layer Resource Architecture

Bio-Framework uses a three-layer priority architecture that allows project-level and user-level customization to override framework defaults without modifying any plugin files.

```
Priority 1 (Highest) — Project Layer
    .claude/
    Project-specific rules, reference code, knowledge, and triggers.
    Applied to the current project only.
          |
          | not found at this layer
          v
Priority 2 — User Global Layer
    ~/.claude/bio-framework/
    Personal default preferences applied across all projects.
    Overrides framework defaults for the individual user.
          |
          | not found at this layer
          v
Priority 3 (Lowest) — Plugin Layer
    framework/
    Shipped framework defaults: rules, knowledge, pipelines, reference code.
    Active when no higher-priority override exists.
```

**Merge strategy by resource type:**

| Resource Type    | Strategy | Description                                                  |
|------------------|----------|--------------------------------------------------------------|
| Skill files      | override | Higher-priority layer fully replaces the lower-priority file |
| Rule parameters  | override | Individual parameter values override lower-priority values   |
| Reference code   | override | Function libraries at higher layers replace lower-layer ones |
| Knowledge base   | union    | All layers are merged; higher-priority entries take precedence on key collision |
| Trigger words    | union    | All layers are merged so custom triggers extend, not replace, framework triggers |
| Error patterns   | union    | All layers are merged to build the broadest possible error pattern index |

The resolution order is enforced at runtime by `framework/scripts/layer_resolver.py`. Helper utilities are provided in `framework/scripts/layer_helpers.sh` and `framework/scripts/bio_extract.py`.

---

### 1.2 Core Module Architecture

```
framework/
├── manifest.yml                        Plugin manifest and metadata
├── rules/
│   └── default_rules.yml               Framework default parameter rules
├── triggers/
│   └── default_triggers.yml            Framework default trigger word definitions
├── ref_code/
│   ├── R/                              Curated R reference function library
│   └── Python/                         Curated Python reference function library
├── scripts/
│   ├── layer_resolver.py               Three-layer priority resolution engine
│   ├── layer_helpers.sh                Shell helpers for layer operations
│   ├── bio_extract.py                  Knowledge base extraction utilities
│   ├── consistency_checker.py          Automated consistency checker (13 checks)
│   └── update_version.py               Version management utility
├── knowledge/                          Knowledge base (see Section 1.3)
├── tests/                              Test specifications and reports
│   ├── specs/                          80 test specification files
│   └── reports/                        Test reports (TR-001 through TR-110)
└── skills/                             Skill runtime modules
    ├── SKILL.md                        Skill entry point (551 lines)
    │                                   Auto-loaded by Claude Code; routes all
    │                                   /bio commands into the framework.
    ├── ORCHESTRATOR.md                 Main orchestrator (398 lines)
    │   └── orchestrator/               13 sub-modules (10,318 lines total)
    │       ├── execution_steps.md      Step definitions and gate conditions (443 lines)
    │       ├── checkpoint_recovery.md  Cross-session checkpoint recovery (271 lines)
    │       ├── data_quality_review.md  Step 1.5 quality audit engine (1,840 lines)
    │       ├── reflection_exploration.md  Step 2 reflection and exploration loop (1,128 lines)
    │       ├── research_summary.md     Step 3 research summary and figure planning (1,658 lines)
    │       ├── publication_figures.md  Step 4 publication-quality figures (1,852 lines)
    │       ├── manuscript_drafting.md  Step 6 SCI manuscript drafting (1,574 lines)
    │       ├── final_report.md         Step 5 technical report (221 lines)
    │       ├── thinking_dispatch.md    Thinking mode dispatch table (343 lines)
    │       ├── pipeline_compliance.md  Pipeline compliance checks (489 lines)
    │       ├── NEW_TOPIC_CLEANUP.md    New topic state reset (included in orchestrator)
    │       ├── topic_change_rules.md   Topic change rules and governance (302 lines)
    │       └── claude_router.md        Model routing logic (95 lines)
    ├── WORKFLOW_CONTROLLER.md          Workflow state machine (1,102 lines)
    │                                   Handles /bio start, /bio continue, /bio status,
    │                                   /bio pause, /bio reset and state transitions.
    ├── COMMAND_PARSER.md               Command parsing and multilingual routing (556 lines)
    ├── HELP_SYSTEM.md                  Help system with EN/ZH/JA output (769 lines)
    ├── THINKING_MODES.md               Four-level thinking depth dispatch (456 lines)
    ├── TOPIC_LOADING.md                Topic parsing and project understanding (720 lines)
    ├── CONFIG_MANAGER.md               Parameter configuration management (802 lines)
    ├── AUTO_EXECUTE.md                 Autonomous execution policies and pause conditions
    ├── decision_nodes.md               Decision node definitions and thresholds
    ├── execution_policies.md           Execution policy rules
    ├── environment_isolation.md        Environment isolation specifications
    ├── package_installation_policy.md  R/Python package installation policy
    ├── code_executor.md                Code execution orchestration
    ├── code_generation_strategy.md     Code generation strategy and conventions
    ├── manuscript_checker.md           Manuscript quality checker
    └── reference_assistant.md         Reference and citation assistant
```

**workflow_state.yml — runtime state schema:**

The framework tracks execution state in a `workflow_state.yml` file written to the project's `.claude/` directory. The `current_step` field uses a fixed enumeration:

| Value | Meaning |
|-------|---------|
| `init` | Project initialized, analysis not started |
| `phase_running` | A dynamic analysis phase is executing |
| `analysis_complete` | All phases completed |
| `step1_5_quality_review` | Step 1.5 data quality audit in progress |
| `step2_reflection` | Step 2 reflection and exploration |
| `step3_summary` | Step 3 research summary and figure planning |
| `step4_figures` | Step 4 publication figure generation |
| `step5_report` | Step 5 technical report |
| `step6_manuscript` | Step 6 SCI manuscript drafting |
| `completed` | All steps completed |

The `status` field uses a separate enumeration: `pending | running | paused | completed | failed`.

---

### 1.3 Knowledge Base System

```
framework/knowledge/
├── COMMON_ERRORS.yml       (20 KB) Indexed error patterns with verified solutions.
│                                   Covers R/Python runtime errors, YAML parsing
│                                   failures, and bioinformatics pipeline errors.
├── figure_types.yml        (68 KB) Standard figure type definitions for 7 omics
│                                   categories. Each entry includes type_id, name,
│                                   description, mandatory flag, and validation signals.
├── sci_design_standards.yml (36 KB) SCI paper design standards covering figure
│                                   aesthetics, color palettes, font sizes, panel
│                                   layout, and journal submission requirements.
├── task_quality_signals.yml (20 KB) Quality signals and evaluation thresholds used
│                                   by Step 1.5 and the reflection loop.
├── pipelines/              6 standard analysis pipeline files + ATAC-seq knowledge coverage
│   ├── scrna_seq.yml           Single-cell RNA-seq (scRNA-seq)
│   ├── bulk_rnaseq.yml         Bulk RNA-seq
│   ├── spatial_transcriptomics.yml  Spatial transcriptomics
│   ├── proteomics.yml          Proteomics
│   ├── metabolomics.yml        Metabolomics
│   ├── lipidomics.yml          Lipidomics
│   └── _schema.yml             Pipeline schema definition
│   (Note: ATAC-seq has entries in figure_types.yml but no standalone pipeline file)
├── journals/               6 journal compliance configurations
│   ├── nature.yml              Nature
│   ├── science.yml             Science
│   ├── cell_reports.yml        Cell Reports
│   ├── jci.yml                 Journal of Clinical Investigation
│   ├── nature_methods.yml      Nature Methods
│   ├── plos_one.yml            PLOS ONE
│   └── _schema.yml             Journal config schema definition
├── methods/                Analysis methodology documentation
└── troubleshooting/        Troubleshooting guides
```

**figure_types.yml coverage by omics type:**

| Omics Type | Mandatory Figure Types | Notes |
|------------|------------------------|-------|
| scRNA-seq | UMAP (cell type labeled), violin, dot plot, heatmap | `mandatory: true` enforced at Step 3 gate |
| Bulk RNA-seq | PCA score, volcano, heatmap, enrichment | Standard DEG workflow |
| Spatial transcriptomics | Spatial scatter (tissue layout), spatial expression | Explicit entry; does not fall through to `default` |
| Proteomics | PCA score, volcano, protein heatmap | — |
| Metabolomics | PCA score, pathway enrichment, metabolite heatmap | — |
| Lipidomics | PCA score (lipid-specific ID), lipid class composition | Unique `type_id` values to avoid collision |
| ATAC-seq | Fragment size distribution, peak annotation | Auto-detected via `detection_signals` |

---

## 2. Feature Description

### 2.1 End-to-End Workflow Engine

The framework provides a structured multi-step workflow from raw data to SCI manuscript draft. Steps 0 and 1 are dynamic; Steps 1.5 through 6 are fixed post-analysis stages.

```
Step 0   Project Initialization
         Read TOPIC.yml → Classify research type → Dynamically generate
         Phase plan → Write workflow_state.yml

Phase 1-N  Dynamic Analysis (count and content designed by AI)
         AI designs the number and content of phases based on research topic.
         Each phase: plan → execute → checkpoint → next phase.
         No hard limit on analysis types or phase count.

Step 1.5  Data Quality Deep Audit
         Five audit dimensions:
           (1) Sample identity verification
           (2) Batch effect detection and assessment
           (3) Biological plausibility checks
           (4) Outlier detection and evaluation
           (5) Cross-stage data consistency
         Output: quality_report.md + findings.yml

Step 2   Reflection and Exploration Loop
         Gap analysis → Present choices to user → Execute supplementary
         analyses → Re-evaluate → Repeat until user satisfied.
         Driven by reflection_exploration.md sub-module.

Step 3   Research Summary and Figure Planning
         Story synthesis across all phases → Figure requirement analysis
         → Generate figure_plan.yml → Validate mandatory type_ids
         against figure_types.yml (HARD_STOP if missing mandatory figures).

Step 4   Publication-Quality Figures
         4a-4c: Figure code generation and execution
         4d: AI visual inspection (Read tool — actual visual review)
         4e: table_index.yml and figure index verification
         4f: User visual review gate (HARD_STOP until user confirms)

Step 5   Technical Report
         Methods documentation → Reproducibility package →
         Environment specification → quality_report.md integration.

Step 6   SCI Manuscript Draft
         Auto-drafting → 6 polishing rounds → 5-stage submission
         verification → submission_readiness.md

Completed  All artifacts present and verified.
```

**Key governance rules:**
- `step2_to_step3` transition: `skip_all` still requires minimum file output; no OR-bypass allowed.
- `step6_to_completed` condition: must verify `submission_readiness.md` (final artifact), not only intermediate outputs.
- Any `HARD_STOP` condition must define an AI fallback path or graceful degradation to prevent autonomous execution deadlock.

---

### 2.2 Four-Level Thinking Depth System

The framework dispatches one of four thinking modes for every AI action, controlled by `THINKING_MODES.md` and routed via `thinking_dispatch.md`.

| Level | Tag | Typical Use Cases |
|-------|-----|-------------------|
| Quick | `[QUICK]` | File reads, status updates, directory listings, simple lookups |
| Standard | `[STANDARD]` | Known-pattern code generation, parameter substitution, template filling |
| Deep | `[DEEP]` | Quality assessment, error diagnosis, cross-file consistency checks |
| Ultrathink | `[ULTRATHINK]` | Research topic understanding, cell type annotation, comprehensive conclusions, architectural decisions |

**9 mandatory Ultrathink scenarios** (cannot be downgraded):

1. Initial topic understanding and research type classification (Step 0)
2. Dynamic Phase plan design
3. Cell type annotation (scRNA-seq)
4. Step 1.5 data quality comprehensive conclusion
5. Step 2 reflection gap analysis
6. Step 3 story synthesis and figure requirement analysis
7. Step 4 AI visual inspection
8. Step 6 manuscript completeness verification
9. Any scenario where findings contradict initial hypotheses

---

### 2.3 Command System

**Top-level command summary:**

| Category | Command Count | Examples |
|----------|---------------|---------|
| Workflow control | 6 | `/bio start`, `/bio continue`, `/bio pause`, `/bio reset`, `/bio status`, `/bio skip` |
| Analysis | 5 | `/bio analyze`, `/bio phase`, `/bio quality`, `/bio explore`, `/bio summary` |
| Output generation | 4 | `/bio figures`, `/bio report`, `/bio manuscript`, `/bio export` |
| Knowledge and config | 5 | `/bio config`, `/bio rules`, `/bio knowledge`, `/bio errors`, `/bio journal` |
| Help and navigation | 5+ | `/bio help`, `/bio topics`, `/bio checkpoint`, `/bio history`, `/bio version` |

**Multilingual alias coverage:**

| Language | Prefix Aliases | Command Aliases |
|----------|---------------|-----------------|
| English | `/bio`, `/bioinformatics` | Primary command names |
| Chinese (Simplified) | `/生信` | 87 aliases |
| Japanese | `/バイオ` | 10 trigger word groups |

**Command routing behavior:**
- `COMMAND_PARSER.md` normalizes all input (including natural language) to canonical command tokens.
- Unrecognized input triggers natural language fallback recognition before returning an error.
- All four prefix forms (`/bio`, `/bioinformatics`, `/生信`, `/バイオ`) are equivalent and route to the same handlers.
- `WORKFLOW_CONTROLLER.md` enforces state-machine guards: commands that are invalid in the current `current_step` state are rejected with a descriptive error.

---

### 2.4 Knowledge Accumulation and Error Handling

**Three-layer error resolution:**

```
Layer 1 — Framework Knowledge Base
    COMMON_ERRORS.yml: indexed patterns, confidence scores, verified solutions.
    Fastest lookup; highest reliability for known error classes.
          |
          | no match
          v
Layer 2 — User Experience Library
    ~/.claude/bio-framework/knowledge/ (user global layer)
    Accumulated from the user's previous analysis sessions.
    Covers project-specific or domain-specific errors.
          |
          | no match
          v
Layer 3 — AI Reasoning
    Claude reasons from first principles using available context,
    documentation, and installed package documentation.
    Solution is a candidate for auto-saving to Layer 2.
```

**Auto-learning policy:**
- When a novel error requires Layer 3 reasoning and the solution is verified to work, it is automatically saved to the user experience library.
- Saves include: error pattern, context, verified solution, and session timestamp.
- This incrementally enriches Layer 2 across sessions without modifying the plugin layer.

---

### 2.5 Checkpoint Recovery and Session Management

**Checkpoint mechanics:**

| Event | Checkpoint Action |
|-------|------------------|
| Phase completion | Write `workflow_state.yml` with `current_step`, phase ID, and output file manifest |
| Step transition (1.5, 2, 3, 4, 5, 6) | Write `workflow_state.yml` with updated `current_step` and `status` |
| Anomaly finding detected | Append to `findings.yml`; continue or pause per `execution_policies.md` |
| Session interrupted | State preserved in `workflow_state.yml`; recovery available on next session start |

**Cross-session recovery flow (`checkpoint_recovery.md`):**

1. On `/bio continue`: read `workflow_state.yml` to determine `current_step`.
2. Validate that required files for the current step exist (step-specific `check_files` list).
3. If files are present: resume from the identified step.
4. If files are missing: determine the last valid completed step and resume from there.
5. Step 6 recovery verifies: `quality_report.md` plus 5 submission verification files (6e through 6i).

**Key runtime files:**

| File | Location | Purpose |
|------|----------|---------|
| `workflow_state.yml` | `.claude/` | Primary state tracking; `current_step`, `status`, phase progress |
| `findings.yml` | `.claude/` | Anomaly and quality finding log across all phases |
| `figure_plan.yml` | `.claude/` | Figure requirement plan generated at Step 3 |
| `table_index.yml` | `.claude/` | Table index validated at Step 4e |
| `quality_report.md` | `.claude/` | Step 1.5 comprehensive quality report |
| `submission_readiness.md` | `.claude/` | Step 6 final submission checklist |
| `visual_review_log.yml` | `.claude/` | Step 4 visual review records |

---

## 3. System Characteristics

### 3.1 Adaptive Analysis Engine

The framework imposes no fixed analysis template. At Step 0, the AI reads `TOPIC.yml` and dynamically designs the full Phase plan:

- **Research type auto-classification:** scRNA-seq, bulk RNA-seq, spatial transcriptomics, proteomics, metabolomics, lipidomics, ATAC-seq, multi-omics, clinical, and others are detected from natural language topic descriptions and file contents.
- **Phase count:** Determined by research complexity; no minimum or maximum.
- **Parameter adaptation:** Clustering resolution, QC thresholds, and normalization strategies are adapted to detected data characteristics.
- **New omics type fallback:** For emerging or uncommon analysis types not explicitly covered by the 6 standard pipeline definitions, `research_summary.md` triggers a WebSearch lookup to identify canonical figure types before planning.

---

### 3.2 Production-Grade Quality Assurance

The framework is validated through a systematic multi-session audit program:

| Metric | Value |
|--------|-------|
| Total test specifications | 80 files |
| Latest audit session | Session 76 (2026-03-07) |
| Total audit sessions (S12-S76) | 65 sessions |
| Cumulative bugs fixed | 400+ |
| Consistency checker checks | 13 automated checks |
| Consistency checker pass rate | 13/13 (100%) across all audit sessions |
| Automated test pass rate | 100% on all test reports (TR-001 through TR-110) |

**Consistency checker (`framework/scripts/consistency_checker.py`)** validates 13 cross-module consistency rules on every audit session, including:
- `current_step` enumeration completeness across all modules
- `workflow_state.status` enumeration alignment
- Mandatory Ultrathink scenario count consistency (SKILL.md vs ORCHESTRATOR vs THINKING_MODES)
- Command prefix alias parity (SKILL.md vs COMMAND_PARSER.md)
- `type_id` global uniqueness in `figure_types.yml`
- Checkpoint recovery file list completeness
- Journal config schema compliance

---

### 3.3 Fully Local Execution

All data processing occurs in the user's local R and Python environment:

- **R minimum version:** 4.3.0
- **Python minimum version:** 3.9.0
- **Network usage:** Raw data never leaves the local environment. Claude Code sends code and result summaries to the Anthropic API for model processing, but raw data files are not transmitted.
- **Storage:** All outputs written to the project directory. No cloud storage.
- **Suitability:** Clinical data, patient genomics, proprietary institutional datasets, and any data subject to privacy or confidentiality requirements.

The framework generates R and Python code that runs in the user's local environment. Code generation uses reference libraries from `framework/ref_code/` as quality anchors, with package installation governed by `package_installation_policy.md`.

---

### 3.4 Modular and Extensible

**Extension points at each layer:**

| Extension Point | Project Layer | User Global Layer |
|----------------|---------------|-------------------|
| Analysis rules | `.claude/rules/project_rules.yml` | `~/.claude/bio-framework/rules/user_rules.yml` |
| Trigger words | `.claude/triggers/` | `~/.claude/bio-framework/triggers/` |
| Reference code | `.claude/ref_code/` | `~/.claude/bio-framework/ref_code/` |
| Knowledge base | `.claude/knowledge/` | `~/.claude/bio-framework/knowledge/` |
| Skill overrides | `.claude/skills/` | `~/.claude/bio-framework/skills/` |

No plugin files need to be modified for customization. All overrides are isolated to the project or user global layer and survive framework upgrades.

---

### 3.5 Full SCI Paper Workflow Support

The framework covers the complete path from raw data to a structured first draft:

**Step 4 — Publication-quality figures:**
- Figures generated to journal-specific dimension and DPI standards (from `sci_design_standards.yml`).
- AI performs actual visual inspection (Step 4d) using the Read tool to view generated images.
- User visual review gate (Step 4f) is a hard stop: the workflow does not advance until the user confirms figure quality.
- `mandatory_type_id_check` validates that all `mandatory: true` figures from `figure_types.yml` are present before Step 4 completes.

**Step 6 — SCI manuscript:**
- Auto-drafting of all standard sections (Introduction, Methods, Results, Discussion, Abstract).
- 6 polishing rounds targeting grammar, flow, logical consistency, and journal style.
- 5-stage submission verification producing `submission_readiness.md`.
- Journal compliance checks available for 6 journals: Nature, Science, Cell Reports, JCI, Nature Methods, PLOS ONE.

---

### 3.6 Multilingual Internationalization

| Language | Status | Command Aliases | Help Output |
|----------|--------|-----------------|-------------|
| English (EN) | Primary | All commands | Full |
| Chinese Simplified (ZH) | Full | 87 aliases | Full |
| Japanese (JA) | Full | 10 trigger groups | Full |

**Language detection:** `COMMAND_PARSER.md` detects input language automatically and routes to the matching alias table. When detection is ambiguous, English is used as the fallback.

**Help system:** `HELP_SYSTEM.md` outputs help text in the language matching the current topic's `language` setting (configured in `TOPIC.yml`). Runtime translation is performed by the AI at output time, not pre-compiled.

**Declared language support in manifest:** `languages: ["en", "zh", "ja"]`

---

## 4. Architecture Diagrams

### 4.1 Data Flow Diagram

```
User Input (natural language or command)
    |
    v
SKILL.md  (skill entry point, auto-loaded by Claude Code)
    |
    | Normalize input
    v
COMMAND_PARSER.md  (multilingual alias resolution, command tokenization)
    |
    | Canonical command token
    v
WORKFLOW_CONTROLLER.md  (state machine guard)
    |
    | Valid command for current state
    v
ORCHESTRATOR.md  (main execution orchestrator)
    |
    +---> THINKING_MODES.md  (select thinking depth for this action)
    |         |
    |         | [QUICK / STANDARD / DEEP / ULTRATHINK]
    |         v
    |     Claude AI reasoning at selected depth
    |
    +---> TOPIC_LOADING.md  (on new topic: parse TOPIC.yml, classify research type)
    |
    +---> CONFIG_MANAGER.md  (resolve parameters via three-layer lookup)
    |         |
    |         | Layer 1: .claude/rules/project_rules.yml
    |         | Layer 2: ~/.claude/bio-framework/rules/user_rules.yml
    |         | Layer 3: framework/rules/default_rules.yml
    |         v
    |     Resolved parameter set
    |
    +---> Sub-module dispatch (based on current_step)
    |         |
    |         +-- execution_steps.md      (Phase 1-N dynamic analysis)
    |         +-- data_quality_review.md  (Step 1.5)
    |         +-- reflection_exploration.md (Step 2)
    |         +-- research_summary.md     (Step 3)
    |         +-- publication_figures.md  (Step 4)
    |         +-- final_report.md         (Step 5)
    |         +-- manuscript_drafting.md  (Step 6)
    |
    v
Knowledge Base Lookup
    |
    | COMMON_ERRORS.yml / figure_types.yml / pipelines/ / journals/
    v
Code Generation (R / Python)
    |
    | Executed in local R/Python environment
    v
Output Artifacts
    |
    +-- Figures (.pdf, .tiff, .png) in project figures/
    +-- workflow_state.yml  (state update)
    +-- findings.yml        (quality findings)
    +-- Reports and manuscript (.md, .docx) in project output/
    |
    v
User
```

---

### 4.2 Command Routing Diagram

```
User Command (e.g., "/生信 继续", "/bio continue", "/バイオ 続ける")
    |
    v
COMMAND_PARSER.md
    |
    +-- Detect language (EN / ZH / JA)
    +-- Expand prefix alias (/生信 -> /bio, /バイオ -> /bio)
    +-- Resolve command alias (继续 -> continue, 続ける -> continue)
    +-- Tokenize arguments and flags (e.g., --quick, --deep, phase ID)
    |
    v
Canonical token: { command: "continue", args: {...}, flags: {...} }
    |
    v
WORKFLOW_CONTROLLER.md
    |
    +-- Read workflow_state.yml
    +-- Check: is this command valid for current_step?
    |       |
    |       | No --> Return error with valid commands for current state
    |       | Yes --> Continue
    |
    +-- /bio start   --> Check: existing workflow? Confirm or abort.
    |                    Initialize workflow_state.yml, call TOPIC_LOADING.md.
    |
    +-- /bio continue --> Determine current_step from workflow_state.yml.
    |                     Dispatch to appropriate sub-module.
    |                     Handle phase ID format (phase01_scrna -> phase0 key lookup).
    |
    +-- /bio status  --> Read workflow_state.yml, compute progress, format output.
    |
    +-- /bio pause   --> Set status = paused, write checkpoint.
    |
    +-- /bio reset   --> Confirm with user, clear workflow_state.yml.
    |
    +-- /bio skip    --> Write minimum required files for skipped step,
    |                    advance current_step (no OR-bypass allowed).
    |
    v
ORCHESTRATOR.md sub-module execution
    |
    v
Output to user + state update
```

---

### 4.3 Resource Lookup Diagram

```
Request: "What is the clustering resolution parameter?"

    Step 1: Check Project Layer
    +-----------------------------------------+
    |  .claude/rules/project_rules.yml         |
    |                                           |
    |  clustering:                              |
    |    resolution:                            |
    |      value: 1.2      <-- FOUND            |
    |      reason: "high heterogeneity"         |
    +-----------------------------------------+
            |
            | Found at Project Layer
            v
        Use value: 1.2
        (Stop lookup. Do not check lower layers.)


    Step 1: Check Project Layer
    +-----------------------------------------+
    |  .claude/rules/project_rules.yml         |
    |  (clustering.resolution not present)      |
    +-----------------------------------------+
            |
            | Not found. Descend to User Global Layer.
            v

    Step 2: Check User Global Layer
    +--------------------------------------------+
    |  ~/.claude/bio-framework/rules/             |
    |  user_rules.yml                             |
    |                                             |
    |  clustering:                                |
    |    resolution:                              |
    |      value: 0.6      <-- FOUND              |
    |      reason: "personal preference"          |
    +--------------------------------------------+
            |
            | Found at User Global Layer
            v
        Use value: 0.6
        (Stop lookup. Do not check Plugin Layer.)


    Step 1: Check Project Layer   --> Not found
    Step 2: Check User Global     --> Not found
            |
            | Descend to Plugin Layer.
            v

    Step 3: Check Plugin Layer
    +--------------------------------------------+
    |  framework/rules/default_rules.yml         |
    |                                             |
    |  clustering:                                |
    |    resolution:                              |
    |      value: 0.8      <-- FOUND              |
    |      reason: "framework default"            |
    +--------------------------------------------+
            |
            | Found at Plugin Layer (lowest)
            v
        Use value: 0.8
        (Framework default applies.)


Union-type resources (knowledge base, trigger words, error patterns):

    Project Layer entries
            +
    User Global Layer entries
            +
    Plugin Layer entries
            |
            v
    Merged collection
    (Higher-priority layer wins on key collision)
    (All unique entries from all layers included)
```

---

## Appendix: File Count Summary

| Component | Files | Total Lines (approx.) |
|-----------|-------|-----------------------|
| Skill entry (SKILL.md) | 1 | 551 |
| ORCHESTRATOR.md + 13 sub-modules | 14 | 10,318 |
| WORKFLOW_CONTROLLER.md | 1 | 1,102 |
| COMMAND_PARSER, HELP, THINKING, TOPIC, CONFIG | 5 | 3,303 |
| Other skill modules (AUTO_EXECUTE, decision_nodes, etc.) | 8 | ~1,500 |
| Knowledge base (YAML) | 4 core + 6 pipelines + 6 journals | ~2,500 |
| Test specifications | 80 | — |
| Test reports | 110 (TR-001 through TR-110) | — |

---

*Documentation generated: 2026-03-10*
*Bio-Framework v1.0.0 — Public Evaluation License (Non-Commercial)*
