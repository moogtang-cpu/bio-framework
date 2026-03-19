# Bio-Framework

> **Not another AI chatbot. A 350 KB+ structured knowledge-driven analysis navigation system.**

---

## The 30-Second Version

You have omics data (single-cell, spatial, proteomics...) and want to publish in a peer-reviewed journal.

Bio-Framework guides AI through bioinformatics best practices to take you **from raw data QC to publication-ready figures to a manuscript first draft**. Hard gates at critical checkpoints (data quality, statistical integrity, journal formatting) prevent the AI from making mistakes or fabricating data.

**You don't need to write R/Python code yourself** — the AI generates, executes, and debugs it automatically. But you do need to understand your scientific question and make decisions at key judgment points (cell type annotation, analysis direction).

All analysis **runs in your local environment**. Raw data never leaves your machine.

---

## What Is It, Exactly?

Let's start with what it's **not**:

- **Not a compute engine** — it doesn't replace R/Python, Seurat, or DESeq2
- **Not a GPT prompt** — not a "you are a bioinformatics expert, please analyze my scRNA-seq data" system prompt
- **Not a pipeline competitor** — it doesn't replace Nextflow/Snakemake for batch computing; it provides navigation at the analysis decision layer

**Bio-Framework is the Guidance Layer.**

It's a Skill plugin for [Claude Code](https://docs.anthropic.com/en/docs/claude-code/overview). Your R/Python environment still does the computation, but Bio-Framework tells the AI at every step:

- What **must be done** at this step (mandatory step enforcement)
- What **commonly goes wrong** (16 error patterns with code-level fixes)
- What **standards the output must meet** (SCI publication standards, journal compliance)

Its core asset isn't code — it's the **350 KB+ structured knowledge base**: error patterns, SCI design standards, workflow logic, journal compliance rules, and methodology documentation, all stored in machine-readable YAML and queried by the AI in real time.

---

## Key Highlights

| Capability | One-line summary |
|:---|:---|
| **End-to-end automation** | From Step 0 initialization to Step 6 manuscript — 7 steps auto-advancing, pausing only when scientific judgment is needed |
| **4-tier thinking dispatch** | Folder creation gets Quick mode; cell annotation forces Ultrathink — not every task deserves the same level of rigor |
| **5-dimension quality audit** | Step 1.5 hard gate — sample mix-ups or severe batch effects trigger an immediate stop for your decision |
| **6-journal compliance** | Nature/Science/Cell Reports and more — automated publication standard checks, figures aligned to DPI and color requirements |
| **Data integrity protection** | Simulated data generation is banned; every statistic in the manuscript traces back to source files; 5-step pre-submission cross-validation |
| **Checkpoint recovery** | After interruption, `/bio continue` resumes to the exact sub-step — no progress lost |

> Detailed feature descriptions below. For an in-depth comparison with "just using Claude + Ultrathink," see [Why Bio-Framework](why-bio-framework-en.md).

---

## 9 Features That Wet-Lab Scientists Can't Refuse

### 1. No More AI Hallucinations: 4-Tier Thinking Depth

**The problem:** Current AI treats all tasks with the same level of thought. Creating a folder and identifying a cell type get the same reasoning depth. Ask "what cell type is this cluster?" and it glances at a few markers and gives you an answer — no cross-validation, no negative evidence exclusion, no literature comparison.

**Bio-Framework's approach:**

Four built-in thinking levels, automatically dispatched based on task criticality:

| Thinking Level | Trigger | AI Behavior |
|:---|:---|:---|
| **Quick** | File operations, formatting | Immediate execution |
| **Standard** | Routine analysis steps | Standard reasoning |
| **Deep** | Parameter selection, method decisions | Multi-option comparison |
| **Ultrathink** | 9 critical scenarios | Full reasoning chain: Observe → Hypothesize → Validate → Conclude |

**9 mandatory Ultrathink scenarios:**

- Cell type identification and annotation
- Data quality review (Step 1.5)
- Differential expression interpretation
- Biological conclusion derivation
- Key manuscript sections (Results / Discussion)
- Data validation (Step 6e: every statistic traced to source)
- Clustering resolution selection
- Batch effect assessment
- Anomalous result evaluation

In Ultrathink mode, the AI **must** complete the full reasoning chain — it doesn't just give you an answer, it shows you **how it got there**.

---

### 2. No More Poisoned Data: Step 1.5 Five-Dimension Deep Audit

**The problem:** You've been running analysis for two weeks, figures are done, and then you discover the sample labels were swapped. Or a reviewer points out inadequate batch correction, and all downstream analysis is invalidated.

**Bio-Framework's approach:**

After Step 1 (computational analysis) and before Step 2 (interpretation), the framework **mandates** a Step 1.5 data quality audit. Five dimensions, none skippable:

| Audit Dimension | What's Checked |
|:---|:---|
| **Sample identity** | Label correctness, cross-contamination signals |
| **Batch effects** | Correction adequacy, residual batch effect quantification |
| **Technical quality** | Sequencing depth, gene detection counts, mitochondrial ratios |
| **Biological plausibility** | Known markers appearing in expected cell populations |
| **Statistical assumptions** | Sample size supporting conclusions, multiple testing correction |

Results are written to `quality_report.yml`. **The workflow only continues if everything passes.** This isn't a suggestion — it's a HARD_STOP gate.

---

### 3. No More Format Rejections: Built-in 6-Journal Compliance Engine

**The problem:** The science is fine, but you get desk-rejected because of insufficient resolution, wrong figure caption format, or abstract word count exceeded.

**Bio-Framework's approach:**

Built-in compliance configurations for **6 major journals** (machine-readable YAML, not natural language suggestions):

| Journal | Key Constraints (examples) |
|:---|:---|
| **Nature** | Abstract ≤150 words, main text ≤3000 words, references ≤50, avoid "reveals/shows/uncovers" in titles |
| **Science** | Research article/report format constraints |
| **Cell Reports** | Figure specifications, STAR Methods required |
| **JCI** | Clinical research requirements, data sharing statement |
| **Nature Methods** | Method validation benchmarks, reproducibility statement |
| **PLOS ONE** | Open access format, Data Availability required |

During Step 6 (manuscript drafting), the AI automatically:

- Verifies figure resolution reaches **300 DPI**
- Checks font sizes meet journal minimums
- Validates **colorblind-friendly** color schemes
- Checks abstract/body word counts against limits
- Verifies Data Availability / Code Availability statements
- Generates a `submission_readiness.md` pre-submission checklist

---

### 4. No More Fabricated Data: Hard Data Integrity Protection

**The problem:** The most unsettling question about AI for researchers — will it fabricate data? If a download fails, will it quietly generate simulated data and keep going? Will the P-values cited in the text actually match the analysis results?

**Bio-Framework's approach:**

**Three lines of defense, zero tolerance:**

**Defense 1: Simulated data generation is banned (HARD_STOP)**

Error F003 in the framework's error library is `severity: critical`. When real data is missing, a download fails, or parsing errors occur, the AI **must stop and wait for user input** — never allowed to "creatively" generate simulated data to continue the workflow.

**Defense 2: Absolute manuscript prohibitions**

- No fabricating experimental data, statistics, P-values, or fold changes
- No inventing literature citations or misrepresenting published research
- No exaggerating results or concealing key limitations
- When information is insufficient, `[Ref]` or `[XX]` placeholders are required — no fabrication

**Defense 3: Five-step pre-submission cross-validation (6e–6i)**

| Stage | What's Validated |
|:---|:---|
| **6e Data Validation** | Every statistic in the text traced to analysis output files (Ultrathink level) |
| **6f Claim Review** | Every conclusion backed by corresponding figures/data |
| **6g Structural Completeness** | All required sections present |
| **6h Cross-File Sync** | Numbers in text cross-checked against figures/supplementary materials ("text says 1,234 DEGs, Table S1 must have exactly 1,234 rows") |
| **6i Submission Readiness** | Final submission readiness report generated |

---

### 5. Six Standardized Omics Pipelines: Not Just "Knowing" — Enforcing

**Coverage:**

| Omics Type | Supported Platforms |
|:---|:---|
| **Single-cell RNA-seq** | 10x Genomics / Seurat / Scanpy |
| **Spatial transcriptomics** | Visium / MERFISH / Slide-seq / STARmap / CODEX |
| **Bulk RNA-seq** | DESeq2 / edgeR / limma |
| **Proteomics** | MaxQuant / DIA-NN / TMT / iTRAQ |
| **Metabolomics** | XCMS / MetaboAnalyst / MZmine |
| **Lipidomics** | Inherits metabolomics + lipid-specific extensions |

Each pipeline defines **mandatory steps**. The AI checks in real time whether your analysis covers every required step — skipped doublet removal? Warning. No batch correction? Warning. Not a post-hoc fix; it's in-process interception.

**Modular use**: You don't have to run Step 0 through Step 6. Already have analysis results? Start directly from Step 4 (figure generation) or Step 6 (manuscript drafting).

---

### 6. Self-Learning Experience Store: Gets Smarter the More You Use It

**How it works:**

```
Problem encountered → Check framework knowledge base (16 known errors)
                   → Check user history (learned_solutions.yml)
                   → AI autonomous reasoning
                   → Resolution took > 2 rounds → Auto-save to user knowledge base
```

Seurat v5 API breaking changes from your first project, CellChat memory overflow fixes, R package compilation errors — all automatically recorded. When the same problem appears in your next project, the stored solution is applied immediately, with a `times_reused` counter tracking usage.

**You're accumulating not just data, but analysis expertise.**

---

### 7. Checkpoint Recovery: Pick Up Where You Left Off

Analysis reached Step 4 (figure generation) and the server crashed?

```bash
/bio continue
```

The system will:
- Scan completed `Figure_*.pdf` files
- Compare against `figure_plan.yml`
- **Resume from the first missing figure, without regenerating completed ones**

Step 6 (manuscript writing) has **21 independent sub-steps**, each producing a separate file. If any step times out or crashes, only that step is retried — completed sections are untouched.

---

### 8. SCI Figure Knowledge Base: Checks What Reviewers Check

Built-in 732-line SCI design standards (`sci_design_standards.yml`) — not teaching you to plot, but helping you **avoid pitfalls**:

- **Figure decision tree** — n<10: dot plot + paired lines; n=10-30: box plot + jittered points; n>30: violin plot
- **Bad figure defense** — 8 common figure mistakes with alternatives (bar charts hiding distributions → box plots, dual Y-axes misleading → faceted plots)
- **AI writing authenticity detection** — Identifies 7 signature patterns of AI-generated text (overuse of "Furthermore/Moreover," uniform paragraph lengths, etc.) to help you write in an authentic expert voice
- **Reviewer checklist** — 15+ frequently raised reviewer concerns, each with specific remediation

---

### 9. Three-Layer Customizable Parameters: Not a Black Box

```
Project .claude/rules/     ← Project-level overrides (highest priority)
User ~/.claude/rules/      ← Your global preferences
Plugin framework/rules/    ← Framework defaults
```

All 12 parameter domains are customizable: QC thresholds, normalization methods, dimensionality reduction parameters, clustering resolution, DE methods, batch correction, cell annotation, visualization style, performance settings, output formats...

Your accumulated experience ("my datasets usually have low MT%, set QC threshold to 15%", "I prefer higher clustering resolution at 0.6") can be saved as global preferences. All new projects inherit them automatically; individual projects can override at any time.

---

## Full Workflow Overview

```
Step 0   Initialization ──── Auto-detect omics type, create standard directory structure
  │
Step 1   Dynamic Analysis ── Multi-phase computation (Phase 1-N), mandatory step checks per phase
  │
Step 1.5 Quality Audit ───── Five-dimension deep review, HARD_STOP gate
  │
Step 2   Reflection ──────── Scientific question completion scoring (80%/70%), directed exploration
  │
Step 3   Research Summary ── Core findings, figure planning
  │
Step 4   Publication Figures  SCI-standard output + visual review loop
  │
Step 5   Report Generation ─ Structured analysis report
  │
Step 6   Manuscript Draft ── Journal compliance + five-step cross-validation (6e-6i)
```

---

## Data Security & Privacy

Bio-Framework runs entirely within the Claude Code local environment:

- **Data stays on your machine**: All raw data (h5, rds, csv, etc.) is read and processed only by your local R/Python environment
- **Local code execution**: AI-generated analysis code runs on your machine, not in the cloud
- **API communication scope**: Claude Code communicates with Anthropic's servers, sending code logic and output summaries — raw data matrices are not included in the transmission

> **Note**: For scenarios with strict data security requirements (e.g., IRB-regulated clinical data), users should assess whether this meets their institution's data governance policies. See [Anthropic's Privacy Policy](https://www.anthropic.com/policies/privacy) for details.

---

## Prerequisites & Cost

Before getting started:

| Requirement | Details |
|:---|:---|
| **Claude Code** | Anthropic's official [Claude Code](https://docs.anthropic.com/en/docs/claude-code/overview) CLI tool (requires Max plan or API credits) |
| **Local environment** | R >= 4.3 and/or Python >= 3.9, Conda recommended |
| **Operating system** | macOS / Linux (WSL2 supported) |
| **Hardware** | 16GB+ RAM recommended; for large single-cell datasets (>50k cells), 64GB+ recommended |
| **Bio-Framework** | One-time purchase, all future updates included |

**Total cost = Bio-Framework one-time purchase + Claude Code usage fees (per Anthropic's pricing).**

---

## Quick Start

**Step 1: Install the framework**

```bash
./install.sh --global
```

**Step 2: Initialize a project in your data directory**

```bash
cd your_project && bio-init
```

Fill in the generated `TOPIC.yml` — describe your research question, data path, and target journal in natural language.

**Step 3: Start your analysis**

```
/bio start --topic TOPIC.yml
```

The framework auto-detects your omics type, loads the corresponding pipeline, and drives the full workflow. You can check in at any time:

```
/bio status      # Check progress
/bio continue    # Resume from checkpoint
/bio help        # Help
```

---

## Expected Efficiency (Based on Design Estimates)

| Phase | Traditional (estimate) | With Framework (estimate) | Notes |
|:---|:---|:---|:---|
| QC + normalization + clustering | 3-7 days | Hours | Standard workflow automation |
| Cell annotation | 2-5 days | 1-2 days | AI-assisted, still requires manual verification |
| Publication figures | 3-7 days | 1-2 days | Auto-generated + visual review loop |
| Manuscript first draft | 2-4 weeks | 1-3 days | Step 6 auto-generated, requires manual review |

> **Actual efficiency depends on data complexity, research question difficulty, and user involvement.** Complex multi-sample, multi-batch datasets may require more time. These estimates reflect design expectations, not statistically validated benchmarks from large-scale user studies.

---

## What Bio-Framework Is NOT For

Staying honest:

| Scenario | Why | Alternative |
|:---|:---|:---|
| Production-scale pipelines | The framework is interactive; not suited for sequencing centers processing hundreds of samples daily | Nextflow / Snakemake |
| Algorithm development | The framework uses existing tools (Seurat, DESeq2, etc.); not suited for methods innovation | Direct programming |
| Pure upstream analysis | Genome assembly, alignment, and other compute-intensive tasks | BWA / GATK / SPAdes |
| Non-Claude Code environments | The framework is a Claude Code Skill plugin and depends on that platform | — |

> Full limitations in [Why Bio-Framework — Section 6](why-bio-framework-en.md#6-limitations-and-considerations).

---

## Learn More

| You want to... | Go here |
|:---|:---|
| Compare with "just using Claude" in depth | [Why Bio-Framework](why-bio-framework-en.md) |
| See detailed features and use cases | [Product Showcase](product-showcase-en.md) |
| Install and configure | [Installation Guide](installation-guide-en.md) |
| Understand the architecture | [Architecture Docs](architecture-en.md) |
| Browse all 60+ commands | [Command Reference](../framework/manifest.yml) |

---

## Get Bio-Framework

Bio-Framework is available now on Gumroad.

**Get it here:** [https://howler26873.gumroad.com/l/gspdf](https://howler26873.gumroad.com/l/gspdf)

Whether you're a wet-lab PhD student entering bioinformatics or a PI standardizing your team's analysis workflow — Bio-Framework makes the AI work to your field's standards, not the other way around.

**One purchase, all updates.** Your analysis guidance layer, starting today.
