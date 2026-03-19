# Bio-Framework v1.0 - Quick Start Guide

## Installation

### Step 1: Global Installation
```bash
tar -xzf bio-framework-v1.0.0.tar.gz
cd bio-framework-v1.0.0
./install.sh --global
```

### Step 2: Initialize Your Project
```bash
cd /your/project/directory
bio-init
```

Or specify path directly:
```bash
bio-init --path /your/project/directory
```

## First Steps

1. Open Claude Code CLI in your project directory
2. Try: `/bio help`
3. Check version: `/bio version`
4. Show configuration: `/bio config show`

## Dependencies

### Required
- Claude CLI >= 2024.01
- claude-router >= 1.0.0 (prompted for installation by install.sh)

### Optional (for analysis)
- R >= 4.3.0
- Python >= 3.9.0
- Project-specific packages (defined in TOPIC.yml)

