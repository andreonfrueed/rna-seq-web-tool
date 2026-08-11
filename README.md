# 🧬 RNA 分析小助手 · RNA-seq Web Tool

![Python](https://img.shields.io/badge/Python-3.10-blue)
![Streamlit](https://img.shields.io/badge/Streamlit-1.61-FF4B4B)
![Platform](https://img.shields.io/badge/Platform-Windows%2010%2F11%20%2B%20WSL2-lightgrey)
![Tests](https://img.shields.io/badge/tests-51%20passed-brightgreen)

**中文** | [English](#-english)

---

## 中文

一个给**湿实验室**用的本地 RNA 测序分析网页工具：不用写代码，在浏览器里上传测序数据，剩下的比对、数基因、找差异、做富集全自动完成。**数据全程在你自己的电脑上，不上传任何服务器。**

### ✨ 功能特性

- **一键安装**：双击 `一键安装.bat`，自动装好 WSL 环境、Miniconda 和全部分析工具
- **参考文件自动准备**：人（GRCh38）/ 小鼠（GRCm39）基因组与注释自动从 Ensembl 下载并校验，下一次永久离线可用
- **智能样本识别**：自动配对 `样本_R1/R2`（或 `_1/_2`）文件，双端缺对、命名不规范都会明确提示
- **开跑前自检**：磁盘空间（能穿透 WSL 虚拟盘检测到真实的 Windows 盘）、FASTQ 完整性、分析环境，有问题当场拦下，而不是几小时后才发现
- **实时进度与断线重连**：浏览器关了重开也能接上进度；中途想停，一键安全停止
- **结果一目了然**：火山图/热图/PCA 直接预览，表格和图单个下载或打包 ZIP；GO/KEGG 富集按「比较 × 上调/下调」分别出结果，目录直接标明方向（如 `LPS高于C`）
- **放心清理**：一键删除 BAM、索引等中间大文件，结果表格与图片完整保留

### 🔬 分析流程

```mermaid
flowchart LR
    A[上传 FASTQ] --> B[质控/修剪<br/>FastQC · Trim Galore]
    B --> C[比对<br/>HISAT2 / STAR]
    C --> D[基因定量<br/>Genomic Overlaps]
    D --> E[标准化 RPKM]
    E --> F[差异分析<br/>pydiffexpress]
    F --> G[聚类与可视化<br/>PCA · t-SNE · 热图]
    G --> H[GO / KEGG 富集]
    H --> I[HTML 报告]
```

### 🚀 快速开始

**环境要求**：Windows 10/11，内存 ≥ 16GB（人源数据建议 32GB），硬盘剩余 ≥ 50GB。

1. 把整个文件夹拷到电脑上，双击 **`一键安装.bat`**（约 30–60 分钟，全自动）
2. 双击 **`打开分析网页.bat`**，浏览器自动打开 http://localhost:8501
3. 网页共 5 个页面，顺着点就行：
   - 🔍 **环境体检**：装没装好，一眼可见
   - 📤 **数据与参考文件**：选物种 → 下载参考文件 → 传数据（大文件可拖进上传文件夹再扫描）
   - 👥 **分组与参数**：每个样本勾选分组（第一个为对照组），阈值默认即可
   - ▶️ **运行分析**：起个名字 → 开始 → 看实时进度
   - 📥 **结果下载**：预览图片、下载表格、打包 ZIP、清理中间文件

详细说明见 [使用说明.md](使用说明.md)。

### 🗂 项目结构

```
├── 一键安装.bat / 打开分析网页.bat / 打开上传文件夹.bat / 更新网页.bat
├── setup_env.sh              # WSL 内的一键安装脚本
├── app/
│   ├── app.py                # Streamlit 页面与交互
│   ├── lib/                  # 后端：运行器、预检、参考文件、结果、富集等
│   └── tests/                # pytest 测试（51 项）
└── archive/                  # 归档的一次性数据修复脚本
```

### 🧪 测试

```bash
cd app && python -m pytest tests/    # 51 passed
```

### 🙏 致谢

分析引擎基于 [PySeqRNA](https://github.com/navduhan/pyseqrna)（GPL-3.0），并使用了 HISAT2、STAR、FastQC、Trim Galore、gseapy / Enrichr、Streamlit 等优秀开源项目。

---

## 🇬🇧 English

A **local-first RNA-seq analysis web app for wet labs**: upload sequencing data in your browser and get alignment, quantification, differential expression and enrichment analysis done automatically — no coding required. **Your data never leaves your computer.**

### ✨ Features

- **One-click install**: double-click `一键安装.bat` to set up WSL, Miniconda and all bioinformatics tools automatically
- **Automatic references**: human (GRCh38) / mouse (GRCm39) genome and annotation downloaded from Ensembl with content validation, then available offline forever
- **Smart sample detection**: auto-pairs `sample_R1/R2` (or `_1/_2`) FASTQ files, with clear messages for missing mates or bad names
- **Pre-flight checks**: disk space (sees through the WSL virtual disk to the real Windows drive), FASTQ integrity and environment readiness — problems are caught before a 10-hour run, not after
- **Live progress & reconnect**: close the browser and come back later, the progress page reconnects; a safe stop button is always available
- **Friendly results**: preview volcano/heatmap/PCA plots inline, download tables individually or as one ZIP; GO/KEGG enrichment is run separately per comparison × up/down, with self-explanatory folder names (e.g. `LPS高于C`)
- **Safe cleanup**: delete BAM/index intermediates in one click while keeping every table and figure

### 🔬 Pipeline

FASTQ → QC/trimming (FastQC, Trim Galore) → alignment (HISAT2 default, STAR optional) → quantification (genomic overlaps) → RPKM normalization → differential expression (pydiffexpress) → clustering & plots (PCA, t-SNE, heatmaps) → GO/KEGG enrichment → HTML report.

### 🚀 Quick start

**Requirements**: Windows 10/11, ≥ 16GB RAM (32GB recommended for human data), ≥ 50GB free disk.

1. Copy this folder to your PC and double-click **`一键安装.bat`** (30–60 min, fully automatic)
2. Double-click **`打开分析网页.bat`** — your browser opens http://localhost:8501
3. Walk through the 5 pages: environment check → species & data → grouping & parameters → run → results

See [使用说明.md](使用说明.md) for the full guide (Chinese).

### 🗂 Project structure

```
├── *.bat                     # one-click install / launch / update helpers (Windows)
├── setup_env.sh              # installer script inside WSL
├── app/
│   ├── app.py                # Streamlit UI
│   ├── lib/                  # backend: runner, preflight, reference, results, enrichment
│   └── tests/                # pytest suite (51 tests)
└── archive/                  # archived one-off data-repair scripts
```

### 🧪 Tests

```bash
cd app && python -m pytest tests/    # 51 passed
```

### 🙏 Acknowledgements

Built on top of [PySeqRNA](https://github.com/navduhan/pyseqrna) (GPL-3.0), HISAT2, STAR, FastQC, Trim Galore, gseapy / Enrichr and Streamlit.
