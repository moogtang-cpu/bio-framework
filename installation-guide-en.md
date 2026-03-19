# Bio-Framework v1.0.0 — Installation and Usage Guide

**Product**: Bio-Framework — AI-Powered Bioinformatics Analysis Framework
**Version**: v1.0.0
**Type**: Claude Code Skill Plugin
**Author**: Hanqing Tang `<moogtang@gmail.com>`
**License**: MIT
**Release Date**: 2026-02-25

---

## Table of Contents

1. [Product Overview](#1-product-overview)
2. [Technical Specifications](#2-technical-specifications)
3. [Data Security and Privacy](#3-data-security-and-privacy)
4. [System Requirements](#4-system-requirements)
5. [Installation](#5-installation)
6. [Quick Start Guide](#6-quick-start-guide)
7. [Workflow Overview](#7-workflow-overview)
8. [FAQ](#8-faq)
9. [Getting Help](#9-getting-help)

---

## 1. Product Overview

Bio-Framework is an end-to-end bioinformatics analysis framework built as a Claude Code skill plugin. It takes researchers from raw experimental data all the way to a publication-ready SCI manuscript draft, orchestrating every analysis step through AI-guided workflows.

The framework is designed around four principles:

- **Systematic**: A structured 7-step workflow ensures no critical analysis stage is skipped, from quality control through manuscript writing.
- **Flexible**: A three-layer configuration system (Project > User Global > Plugin) lets individual researchers or teams override defaults at any level without touching the core framework.
- **Trainable**: Your parameter preferences, problem-solving experience, and custom knowledge accumulate over time — the more you use it, the better it understands your analysis style.
- **Rigorous**: 76+ sessions of deep auditing, 400+ bugs fixed, and a 898-case automated test suite running at 100% pass rate guarantee production-quality stability.

### What Bio-Framework Is Not

Bio-Framework is a **skill plugin for Claude Code**. It does not perform computation itself — it provides AI-guided analysis orchestration, code generation, and interpretation logic. Actual computational work (R scripts, Python scripts) runs in your local environment.

---

## 2. Technical Specifications

| Attribute | Value |
|---|---|
| Core skill code | 10,318+ lines |
| Knowledge base size | 350 KB+ total |
| Error pattern library | 20 KB |
| Standard figure definitions | 68 KB |
| SCI design standards | 36 KB |
| Top-level commands | 25+ |
| Sub-commands | 50+ |
| Multilingual aliases | 110+ |
| Supported languages | English, Chinese (中文), Japanese (日本語) |
| Test cases | 898 automated test cases (100% pass rate) |
| Audit rounds completed | 76+ |
| Bugs fixed | 400+ |

### Supported Omics Types

| Omics Type | Status |
|---|---|
| Single-cell RNA-seq (scRNA-seq) | Fully supported |
| Spatial Transcriptomics | Fully supported |
| Bulk RNA-seq | Fully supported |
| Proteomics | Fully supported |
| Metabolomics | Fully supported |
| Lipidomics | Fully supported |
| ATAC-seq | Fully supported |

### Supported Target Journals

| Journal | Impact Factor Range |
|---|---|
| Nature | Top-tier |
| Science | Top-tier |
| Cell Reports | High |
| Journal of Clinical Investigation (JCI) | High |
| Nature Methods | High |
| PLOS ONE | Standard |

---

## 3. Data Security and Privacy

**This section is critical for institutional and clinical research users. Please read carefully.**

### Fully Local Execution

All analyses run on your local server or computer. Your experimental data never leaves your local environment under any circumstances. Bio-Framework contains no network transmission capabilities and no cloud storage components.

### Claude Does Not Access Your Data

Bio-Framework is a skill plugin that runs within Claude Code. It provides analysis guidance, generates R and Python code, and interprets results. It does **not** upload, transmit, or store any user experimental data.

The data flow is strictly one-directional and local:

```
Your Raw Data (local)
       |
       v
Local R / Python Processing
       |
       v
Local Output Files (plots, tables, reports)
       |
       v
Claude AI reads: code logic, file paths, and result summaries
                 (never raw experimental data)
```

### What Claude AI Sees vs. What It Does Not See

| Claude AI Sees | Claude AI Does NOT See |
|---|---|
| The R/Python code it generates | Raw sequencing reads (.fastq, .bam) |
| File path structures | Expression matrices (raw counts) |
| Summary statistics you paste | Patient identifiers or clinical metadata |
| Analysis results you share | Proteomic or metabolomic raw tables |
| Configuration files (TOPIC.yml) | Any file unless you explicitly paste it |

### Zero Data Leakage Risk

All data processing is performed exclusively by local R and Python environments. The framework contains no mechanisms for data exfiltration. This architecture makes Bio-Framework suitable for:

- Clinical trial data
- Patient sample data with privacy requirements
- Proprietary institutional datasets
- Data subject to HIPAA, GDPR, or equivalent regulations

### Important Caveat: Claude Code API Communication

Bio-Framework runs inside Claude Code, which communicates with Anthropic's API to power the AI. This means that the **code Claude generates, result summaries you share, and conversational context** are sent to Anthropic's servers as part of normal Claude Code operation. This is inherent to how Claude Code works and is not specific to Bio-Framework.

**What is sent to Anthropic**: Code snippets, file path summaries, result summaries, and any text you paste into the chat.

**What is NOT sent**: Raw data files (FASTQ, BAM, count matrices, raw tables). Your data files remain entirely on your local machine and are never transmitted.

If your institution has strict data governance requirements that prohibit any analysis-related text from leaving on-premises servers, consult your data security officer before using any Claude Code-based tool.

### Verification

You can verify the framework contains no network calls by inspecting the source:

```bash
# Confirm no curl/wget/http calls exist in framework logic files
grep -r "curl\|wget\|http\|upload\|transmit" ~/.claude/bio-framework/framework/skills/
# Expected output: empty (no matches)
```

---

## 4. System Requirements

### Required Software

| Component | Minimum Version | Purpose |
|---|---|---|
| Claude Code CLI | >= 2024.01 | Runtime environment for the plugin |
| Operating System | Linux / macOS / WSL2 | Windows users should use WSL2 |

### Analysis Dependencies (Automatically Installed)

Bio-Framework **automatically detects and installs** all required dependencies during the analysis process, including R, Python base environments and specific analysis packages (such as Seurat, DESeq2, scanpy, etc.). You do not need to pre-install any software manually.

During Step 0 initialization, the framework identifies the dependencies needed for your specific analysis based on the omics type and generates installation code for execution when a package is first needed.

### Optional Component

**claude-router plugin**: An intelligent model routing plugin that can reduce API usage costs (savings vary by workload). The installation script will prompt you to install it. Recommended.

### Hardware Recommendations

| Omics Type | Minimum RAM | Recommended RAM |
|---|---|---|
| Bulk RNA-seq (standard) | 8 GB | 16 GB |
| scRNA-seq (< 10k cells) | 16 GB | 32 GB |
| scRNA-seq (> 50k cells) | 32 GB | 64 GB+ |
| Spatial Transcriptomics | 32 GB | 64 GB+ |

---

## 5. Installation

### Step 1: Install Claude Code CLI

If you have not already installed Claude Code CLI, follow the official installation instructions at [claude.ai/code](https://claude.ai/code). Verify the installation:

```bash
claude --version
# Expected: claude 2024.xx.xx or later
```

### Step 2: Obtain Bio-Framework

After purchase, you will receive a `bio-framework-v1.0.0.tar.gz` archive file.

```bash
# Copy the archive to your working directory
cp ~/Downloads/bio-framework-v1.0.0.tar.gz ./

# Extract the archive
tar -xzf bio-framework-v1.0.0.tar.gz
cd bio-framework-v1.0.0
```

### Step 3: Run the Installation Script

```bash
# Make the script executable
chmod +x install.sh

# Install to your global user layer (recommended)
./install.sh --global
```

The installer will:

1. Copy framework files to `~/.claude/bio-framework/`
2. Register the `/bio` skill command in Claude Code
3. Prompt you to install `claude-router` for cost optimization
4. Verify the installation is complete

**Installation Options**

| Option | Description |
|---|---|
| `--global` / `-g` | Install to user global layer `~/.claude/bio-framework/` (recommended) |
| `--project` / `-p` | Install to the current project directory only |
| `--upgrade` / `-u` | Clean reinstall — removes old `framework/` before copying fresh files |
| `--skip-router` | Skip `claude-router` installation (runs in degraded mode without cost optimization) |
| `--help` / `-h` | Show help message |

### Step 4: Initialize a Project

After global installation, navigate to your research project directory and run:

```bash
cd /path/to/your/project
bio-init
```

This creates the necessary project-level configuration structure under `.claude/` in your project directory.

Alternatively, specify the path directly:

```bash
bio-init --path /path/to/your/project
```

### Step 5: Verify Installation

Open Claude Code CLI in your project directory and run:

```
/bio version
```

Expected output: Bio-Framework version information confirming `v1.0.0`.

```
/bio help
```

Expected output: Full command reference with all available `/bio` commands.

### Step 6 (Optional): Install claude-router for Cost Optimization

`claude-router` is a companion plugin that automatically routes requests to the most appropriate Claude model (Haiku for quick operations, Sonnet for standard analysis, Opus for critical scientific decisions). The `install.sh` script will prompt you to install it.

**Benefits of claude-router**:

| Benefit | Detail |
|---|---|
| Cost reduction | Significant cost reduction (varies by workload) |
| Speed improvement | ~3x faster for simple operations (Haiku routing) |
| Quality improvement | Critical decisions use Opus automatically |

**Model routing used by Bio-Framework**:

| Thinking Mode | Routed Model | Used For |
|---|---|---|
| `ultrathink` | Opus | Critical scientific decisions, manuscript drafting |
| `deep` | Sonnet | Complex analysis, figure generation |
| `standard` | Sonnet | Standard workflow steps |
| `quick` | Haiku | Simple status checks and lookups |

If you skip `claude-router`, Bio-Framework runs in degraded mode — all thinking modes still function via prompt-based guidance, but automatic model switching is disabled.

To install `claude-router` later:

```
/plugin marketplace add 0xrdan/claude-plugins
/plugin install claude-router
```

---

## 6. Quick Start Guide

Bio-Framework supports two ways to start a new project, depending on your preparation stage.

---

### Option A: Start from TOPIC.yml (Standard Flow)

Best for: You have clear scientific questions and want the framework to plan the complete analysis.

**Step 1: Create Your Research Topic File**

Create a `TOPIC.yml` file in your project directory. This file defines your scientific questions, hypotheses, and analysis configuration.

A minimal example:

```yaml
metadata:
  project_name: "Wound Healing scRNA-seq"
  author: "Your Name"
  research_type: "basic_research"
  quality_requirements: "standard"

scientific_questions:
  - Q1: "What cell populations are present in wound healing tissue?"
  - Q2: "How do fibroblast subpopulations change across healing stages?"
  - Q3: "What are the key signaling pathways driving macrophage polarization?"

hypotheses:
  - H1: "Myofibroblasts expand during the proliferative phase of wound healing"

biological_context:
  organism: "Mus musculus"
  tissue: "Skin"
  condition: "Wound healing (Day 0, 3, 7, 14)"
  key_cell_types: ["fibroblasts", "macrophages", "keratinocytes", "endothelial cells"]
```

A full template with all available options is available at:
`~/.claude/bio-framework/framework/templates/TOPIC_template.yml`

**Step 2: Start the Analysis**

Open Claude Code CLI in your project directory, then run:

```
/bio start
```

Bio-Framework will initialize the workflow, load your research configuration, and begin Phase 0 (project initialization).

---

### Option B: Start from Pre-Prepared Materials (Data Description + Analysis Plan)

Best for: You already have dataset descriptions and/or a detailed analysis plan, and want the AI to work directly from those materials.

This approach is common when:
- You have downloaded datasets from GEO or similar repositories and compiled data descriptions
- You already have a concrete bioinformatics analysis design (Phase structure, scientific questions)
- You want the AI to understand your analysis plan and begin execution directly

**Step 1: Prepare Your Material Files**

Place the following files in your project directory (naming is flexible):

**Data description file** (e.g., `data/data_info.txt`) — describes your datasets:

```
GSE165816:

Organism    Homo sapiens
Experiment type    Expression profiling by high throughput sequencing
Summary    Diabetic foot ulcers (DFUs) are a devastating complication
           of diabetes. We examined the cellular landscape of DFUs
           by single-cell RNA-seq analysis of foot and forearm skin
           specimens...

Overall design    Single-cell RNA sequencing of foot and forearm skin cells...

Sample groups:
GSM5050521    G1: Forearm skin of Healthy non-diabetic subject
GSM5050522    G1A: Foot skin of Diabetic subject without DFU
GSM5050523    G2: Foot skin of subject with healing DFU
...

Data files:
GSE165816_RAW.tar
```

**Analysis plan file** (e.g., `docs/bio_analyze.md`) — describes your analysis design:

```markdown
# HIF-1α and NF-κB Hypoxia-Inflammation Coupling Study

## Scientific Questions Addressable by Bioinformatics

Q1: Which cell types are the main carriers of HIF1A/NFKB1 co-expression
    in diabetic wounds?
Q2: How does the HIF1A/NFKB1 co-expression pattern relate to wound
    healing outcome?
Q3: What are the shared target genes regulated by both HIF-1α and NF-κB?
...

## Bioinformatics Analysis Phase Design

Phase 1: Single-cell core analysis (GSE165816)
  - 1.1 Data preprocessing and QC
  - 1.2 Cell type annotation
  - 1.3 HIF1A/NFKB1 expression analysis
  ...

Phase 2: Bulk RNA-seq multi-cohort validation
Phase 3: Temporal dynamics analysis
Phase 4: Macrophage HIF/NF-κB mechanism deep analysis
```

**Step 2: Start with `@` File References**

Open Claude Code and reference your material files using the `@` symbol:

```
Please start the bioinformatics analysis based on these materials:
@data/data_info.txt
@docs/bio_analyze.md
```

Or start with the command first:

```
/bio start
```

Then provide context via `@` file references in the conversation. The AI will read your data descriptions and analysis plan, understand the project design, and automatically generate the corresponding TOPIC.yml and analysis workflow.

> **Tip**: The key advantage of Option B is that you can describe your project in natural language and free-form text, without needing to follow TOPIC.yml's strict YAML syntax. The AI extracts scientific questions, dataset information, and analysis requirements from your materials and converts them into a structured workflow automatically.

---

### Comparison of the Two Options

| Dimension | Option A: TOPIC.yml | Option B: Pre-Prepared Materials |
|-----------|--------------------|---------------------------------|
| Best for | Users comfortable with YAML | All users |
| Preparation | Write a structured TOPIC.yml | Prepare free-form data descriptions and analysis plans |
| AI comprehension | Parses structured fields directly | Understands natural language descriptions |
| Flexibility | Fixed fields, comprehensive coverage | High freedom, can include any information |
| Use cases | New projects, standardized workflows | Detailed existing analysis designs, multi-dataset projects |

---

### Automatic Execution After Starting

Regardless of which option you choose, the framework executes the workflow automatically once started:

- **Automatic execution**: Code generation, script execution, file I/O, quality metric calculation, R/Python package installation
- **Pauses for user decisions**: Analysis plan confirmation, biological interpretation review, figure visual inspection
- **On exceptions**: Automatically queries the knowledge base for fixes; pauses and reports when unable to self-resolve

You can check progress at any time with `/bio status`.

### Common Commands Quick Reference

| Command | Description |
|---|---|
| `/bio help` | Show full command reference |
| `/bio start` | Start a new analysis workflow |
| `/bio start --quick` | Start in quick mode (faster, skips non-critical validations) |
| `/bio continue` | Resume an interrupted analysis |
| `/bio status` | Show current workflow progress |
| `/bio stop` | Pause the current analysis |
| `/bio config show` | Display current configuration |
| `/bio version` | Show Bio-Framework version |
| `/bio prime` | Re-initialize AI context at the start of a new session |
| `/bio topic show` | Display the loaded research topic |
| `/bio phase list` | List all analysis phases and their status |
| `/bio phase jump <n>` | Jump to a specific phase number |
| `/bio step jump <step>` | Jump to a specific post-analysis step |
| `/bio figure list` | List planned figures |
| `/bio manuscript status` | Show manuscript drafting progress |
| `/bio export` | Export analysis results |

**Chinese aliases are fully supported.** For example:

```
/bio 开始
/bio 继续
/bio 状态
```

**Japanese aliases are fully supported.** For example:

```
/バイオ 開始
/バイオ 継続
```

---

## 7. Workflow Overview

Bio-Framework implements a structured 7-step workflow that covers the complete research lifecycle.

```
/bio start
    |
    v
Step 0: Initialization
  - Load TOPIC.yml
  - Detect omics type
  - Set up project structure
  - Initialize runtime files
    |
    v
Phase 1 to N: Dynamic Analysis  [ULTRATHINK / DEEP]
  - AI determines required analysis phases based on your data
  - Each phase runs specific bioinformatics pipelines
  - Examples: QC, normalization, clustering, cell type annotation,
              trajectory analysis, differential expression, etc.
    |
    v
Step 1.5: Data Quality Audit  [DEEP]
  - Systematic review of analysis results
  - Pipeline compliance validation
  - Decision: proceed / re-analyze / flag concerns
    |
    v
Step 2: Reflection and Exploration  [DEEP]
  - AI reviews all findings across phases
  - Identifies key biological insights
  - Explores unexpected patterns
    |
    v
Step 3: Research Summary and Figure Planning  [ULTRATHINK]
  - Synthesizes findings into coherent narrative
  - Plans the complete figure set for publication
  - Validates mandatory figures for your omics type
    |
    v
Step 4: Publication Figure Generation  [DEEP]
  - Generates all planned figures
  - Applies SCI design standards (journal-specific formatting)
  - Visual inspection and user review
    |
    v
Step 5: Technical Foundation  [DEEP]
  - Technical report covering all methods
  - Reproducibility documentation
  - Supplementary materials
    |
    v
Step 6: SCI Manuscript Drafting  [ULTRATHINK]
  - Full manuscript draft (Introduction through Discussion)
  - Journal-specific formatting applied
  - Submission readiness check
    |
    v
Completed
```

### Thinking Mode by Step

The framework uses different levels of AI reasoning depth at each stage:

| Step | Thinking Mode | Rationale |
|---|---|---|
| Dynamic Phases | ultrathink / deep | Scientific analysis requires deep reasoning |
| Step 1.5 | deep | Systematic quality assessment |
| Step 2 | deep | Comprehensive reflection |
| Step 3 | ultrathink | Critical narrative synthesis and figure planning |
| Step 4 | deep | Accurate figure generation |
| Step 5 | deep | Technical documentation |
| Step 6 | ultrathink | Manuscript requires highest quality reasoning |

---

## 8. FAQ

### Installation fails with "Must run install.sh from Bio-Framework root directory"

You must run the installer from the directory containing the `VERSION` file (the Bio-Framework root). Check your current directory:

```bash
ls VERSION   # Should show: VERSION
./install.sh --global
```

### The `/bio` command is not recognized after installation

1. Restart your Claude Code CLI session — skill plugins are loaded at session start.
2. Verify installation files exist:

```bash
ls ~/.claude/bio-framework/framework/
ls ~/.claude/skills/bio/SKILL.md
```

3. If files are missing, re-run the installer:

```bash
./install.sh --upgrade
```

### bio-init command not found

The `bio-init` script is installed to `~/.local/bin/` or a similar user-local bin directory. Ensure this directory is in your `PATH`:

```bash
echo $PATH | grep -o "$HOME/.local/bin"
# If not present, add to your shell config:
echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

### Analysis stops unexpectedly — how do I resume?

Use `/bio continue` to resume from where the analysis was interrupted. Bio-Framework maintains checkpoint files that allow resumption from any completed step.

```
/bio continue
```

If the automatic checkpoint is not picked up correctly, you can specify a phase:

```
/bio continue --from 3
```

### How do I update Bio-Framework?

Download the latest release archive, then run:

```bash
cd bio-framework-v1.x.x
./install.sh --upgrade
```

The `--upgrade` flag performs a clean reinstall: it removes the old `framework/` directory and copies fresh files. Your project data, configuration rules, and user knowledge files are never affected by an upgrade.

### How do I uninstall Bio-Framework?

```bash
# Remove the global installation
rm -rf ~/.claude/bio-framework
rm -rf ~/.claude/skills/bio

# Remove the bio-init command
rm -f ~/.local/bin/bio-init

# Remove project-level files (if installed per-project)
rm -rf /your/project/.claude/bio-framework
```

### Can I use Bio-Framework without claude-router?

Yes. Bio-Framework runs in degraded mode without `claude-router`. All thinking modes (ultrathink, deep, standard, quick) continue to function through prompt-based guidance. The only difference is that automatic model switching is disabled, meaning all requests go to whichever model you have configured in Claude Code.

To install `claude-router` at any time after the initial setup:

```
/plugin marketplace add 0xrdan/claude-plugins
/plugin install claude-router
```

Then restart your Claude Code session.

### How do I add a custom analysis phase?

Advanced customization is possible through the three-layer configuration system. Create a project-level rules file at `.claude/rules/project_rules.yml` to override framework defaults for your project without modifying core files. See `framework/rules/default_rules.yml` for all configurable parameters.

---

## 9. Getting Help

### In-App Help

The fastest way to get help is through the built-in help system:

```
/bio help
```

For help on a specific topic:

```
/bio help start
/bio help config
/bio help manuscript
```

### GitHub Issues

For bug reports, feature requests, and technical questions:

**GitHub Issues**: [https://github.com/bio-framework/claude-bio/issues](https://github.com/bio-framework/claude-bio/issues)

When filing a bug report, include:

- Bio-Framework version (`/bio version`)
- Operating system and version
- Claude Code CLI version (`claude --version`)
- Steps to reproduce the issue
- The contents of your `TOPIC.yml` (omit any sensitive data)
- Relevant output or error messages

### Author Contact

**Hanqing Tang**
Email: [moogtang@gmail.com](mailto:moogtang@gmail.com)

### Additional Resources

| Resource | Location |
|---|---|
| Project homepage | [https://bio-framework.github.io](https://bio-framework.github.io) |
| Changelog | `CHANGELOG.md` in the project root |
| Framework architecture | `framework/README.md` |
| Migration guide (from v5.4) | `docs/migration/MIGRATION_GUIDE_EN.md` |
| TOPIC.yml full template | `~/.claude/bio-framework/framework/templates/TOPIC_template.yml` |
| Default configuration rules | `~/.claude/bio-framework/framework/rules/default_rules.yml` |

---

*Bio-Framework v1.0.0 — AI-Powered Bioinformatics Analysis Framework*
*Author: Hanqing Tang | License: MIT | Release: 2026-02-25*
