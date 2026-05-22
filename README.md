
<div align="center">

# 🧬 Chroma-3D

### Hi-C → 3D Chromatin Structure Pipeline & Interactive Web Viewer

[![Python](https://img.shields.io/badge/Python-3.8%2B-3776AB?style=for-the-badge&logo=python&logoColor=white)](https://python.org)
[![NumPy](https://img.shields.io/badge/NumPy-1.24%2B-013243?style=for-the-badge&logo=numpy&logoColor=white)](https://numpy.org)
[![scikit-learn](https://img.shields.io/badge/scikit--learn-MDS-F7931E?style=for-the-badge&logo=scikit-learn&logoColor=white)](https://scikit-learn.org)
[![License](https://img.shields.io/badge/License-MIT-green?style=for-the-badge)](LICENSE)
[![WebGL](https://img.shields.io/badge/Renderer-Canvas2D-E34F26?style=for-the-badge&logo=html5&logoColor=white)](index.html)

<br/>

> **Chroma-3D** is an end-to-end bioinformatics pipeline that transforms raw Hi-C contact matrices (`.npy`) into biologically accurate 3D chromatin fiber models, then renders them in an interactive browser-based viewer — no server required.

<br/>

---

</div>

## 📖 Table of Contents

- [Overview](#-overview)
- [Pipeline Architecture](#-pipeline-architecture)
- [Features](#-features)
- [Installation](#-installation)
- [Usage](#-usage)
  - [Step 1 — Python: .npy → .json](#step-1--python-npy--json)
  - [Step 2 — Browser: .json → 3D Viewer](#step-2--browser-json--3d-viewer)
- [Input / Output Format](#-input--output-format)
- [Algorithm Deep Dive](#-algorithm-deep-dive)
- [Viewer Controls](#-viewer-controls)
- [File Structure](#-file-structure)
- [Gene Annotation Support](#-gene-annotation-support)
- [Batch PNG Export](#-batch-png-export)
- [Dependencies](#-dependencies)
- [Contributing](#-contributing)

---

## 🔬 Overview

Chroma-3D bridges the gap between raw Hi-C sequencing data and scientific visualization. Given a Hi-C contact matrix for any genomic region, the pipeline:

1. **Detects TAD boundaries** using a signal-processing Insulation Score algorithm (no hardcoding)
2. **Folds the chromatin fiber** in 3D space using constrained Multidimensional Scaling (MDS) with polymer physics
3. **Exports a compact JSON** containing 3D coordinates, TAD assignments, and genomic metadata
4. **Renders an interactive 3D scene** in the browser using a custom Canvas2D engine with perspective projection, Catmull-Rom splines, and per-TAD coloring

The viewer supports hovering over individual genomic bins to inspect TAD statistics and overlapping gene annotations loaded from a Parquet file.

---

## 🏗 Pipeline Architecture

```
┌─────────────────────────────────────────────────────────────────────────┐
│                        CHROMA-3D PIPELINE                               │
└─────────────────────────────────────────────────────────────────────────┘

  Hi-C Contact Matrix            Python Backend              Browser Viewer
  ┌──────────────────┐     ┌────────────────────────┐     ┌──────────────────┐
  │                  │     │                        │     │                  │
  │  chrN_start_end  │     │  1. Log-normalization  │     │  index.html      │
  │     .npy file    │────▶│  2. Insulation Score   │     │                  │
  │                  │     │     TAD Detection      │     │  • 3D Canvas     │
  │  [N × N matrix]  │     │  3. Distance Matrix    │     │    Renderer      │
  │                  │     │     (power-law)        │────▶│  • TAD Coloring  │
  └──────────────────┘     │  4. MDS 3D Embedding   │     │  • Hover Popups  │
                           │  5. Physics Refinement │     │  • Gene Lookup   │
  Gene Annotations         │  6. JSON Export        │     │  • PNG Export    │
  ┌──────────────────┐     └────────────────────────┘     │                  │
  │  mouse.parquet   │──────────────────────────────────▶│  mouse.parquet   │
  │  (gene_id, name, │                                     │  (WebAssembly)   │
  │   chr, start,end)│                                     └──────────────────┘
  └──────────────────┘
```

---

## ✨ Features

### Python Pipeline (`Batchn2j.py`)

| Feature | Description |
|---|---|
| **Biological TAD Detection** | Insulation Score + Gaussian smoothing + valley-finding; no hardcoded boundaries |
| **Log-Normalization** | `log1p` transform to handle Hi-C's naturally log-normal contact distribution |
| **Power-Law Distance Matrix** | Converts interactions to physical distances using a `0.7` polymer exponent |
| **Backbone Stiffness** | Adjacent bins are constrained to `0.8` units to model a continuous chromatin chain |
| **Pinch-Point Contacts** | Top 1% off-diagonal interactions are identified and pulled into contact distance |
| **MDS 3D Folding** | scikit-learn `MDS` with `n_init=10`, `max_iter=2000` for global minimum |
| **Iterative Physics Refinement** | Double-pass Gaussian smoothing + 3× contact attraction loop |
| **Batch Processing** | Processes entire folders of `.npy` files automatically |
| **Graceful Fallback** | Generates a simulated structure if the input file is missing |

### Web Viewer (`index.html`)

| Feature | Description |
|---|---|
| **Custom 3D Engine** | Zero-dependency Canvas2D renderer with perspective projection and depth-based fog |
| **Catmull-Rom Splines** | Smooth backbone interpolation with 6 sub-steps per bin segment |
| **TAD Coloring** | Per-TAD HSL color assignment using the golden-angle (137°) for maximum contrast |
| **Auto-Orbit** | Smooth auto-rotation with damped camera interpolation |
| **Interactive Hover** | Hover any bin to see its TAD ID, genomic range, size, and overlapping genes |
| **Parquet Gene Annotations** | Loads gene data via WebAssembly (`parquet-wasm`) + Apache Arrow |
| **Batch PNG Export** | Renders a per-bin animation and saves each frame as a PNG to a local folder |
| **High-DPI Support** | Proper `devicePixelRatio` handling for retina displays |
| **Academic Aesthetic** | Parchment/sepia UI theme with serif typography |

---

## 🚀 Installation

### Prerequisites

- Python 3.8 or higher
- A modern browser (Chrome/Edge recommended for File System Access API)

### Clone the Repository

```bash
git clone https://github.com/your-username/chroma-3d.git
cd chroma-3d
```

### Install Python Dependencies

```bash
pip install numpy scipy scikit-learn
```

Or use the requirements file:

```bash
pip install -r requirements.txt
```

**`requirements.txt`**
```
numpy>=1.24.0
scipy>=1.10.0
scikit-learn>=1.2.0
```

---

## 📋 Usage

### Step 1 — Python: `.npy` → `.json`

Place your Hi-C contact matrix files in a folder. Files **must** follow this naming convention so the pipeline can extract genomic coordinates:

```
chrN_START_END.npy

Examples:
  chr9_10000000_20000000.npy
  chrX_5000000_15000000.npy
```

Then run the batch processor:

```bash
python Batchn2j.py
```

By default this reads from `./chr9/` and writes JSON to `./json_output/`. Edit the `process_batch()` call at the bottom of the script to change paths:

```python
if __name__ == "__main__":
    process_batch(
        input_dir="./my_hic_data",   # Folder with .npy files
        output_dir="./json_output"   # Output folder for .json files
    )
```

**Expected output:**

```
Found 3 files.

Folding polymer with 47 biological contacts...
Saved: json_output/chr9_10000000_20000000.json
Realism Success! 12 TADs detected. Exported 200 bins.

Folding polymer with 61 biological contacts...
Saved: json_output/chr9_20000000_30000000.json
Realism Success! 15 TADs detected. Exported 200 bins.
```

---

### Step 2 — Browser: `.json` → 3D Viewer

1. Copy your generated JSON file and rename or symlink it to `chromatin_data.json` in the same directory as `index.html`:

```bash
cp json_output/chr9_10000000_20000000.json ./chromatin_data.json
```

2. Optionally, place a `mouse.parquet` gene annotation file in the same directory.

3. Serve the directory with any static file server (required for `fetch()` and Parquet loading):

```bash
# Python built-in server
python -m http.server 8080

# Node.js (npx)
npx serve .

# VS Code Live Server extension also works
```

4. Open `http://localhost:8080` in your browser.

> **Note:** If no `chromatin_data.json` is present, the viewer automatically generates a procedural fallback structure so the UI is always functional.

---

## 📂 Input / Output Format

### Input — Hi-C Contact Matrix (`.npy`)

```
Shape: [N, N]  — symmetric square matrix
dtype: float32 or float64
Values: Raw or KR/ICE normalized Hi-C contact counts
Filename: chrN_START_END.npy  (genomic coordinates in bp)
```

### Output — 3D Structure JSON

```json
{
  "x": [0.12, -3.45, ...],        // N x-coordinates (units: viewer space)
  "y": [1.23, 0.67, ...],         // N y-coordinates
  "z": [-0.89, 2.11, ...],        // N z-coordinates
  "Chromosom": "chr9",            // Chromosome name from filename
  "start": 10000000,              // Genomic start (bp)
  "end": 20000000,                // Genomic end (bp)
  "colors": ["#3b82f6", ...],     // Per-bin hex color (by TAD)
  "tad_boundaries": [0, 18, 45, 99, ...]  // Bin indices of TAD boundaries
}
```

---

## 🧪 Algorithm Deep Dive

### 1. Log-Normalization

Raw Hi-C counts span many orders of magnitude. A `log1p` transform compresses this range so that moderate interactions contribute meaningfully to the 3D fold:

```python
heatmap_log = np.log1p(heatmap)
```

### 2. TAD Detection via Insulation Score

Rather than hardcoding domain boundaries, Chroma-3D computes the **Insulation Score** — the mean local interaction density in a sliding window around each bin. Valleys in this signal correspond to CTCF binding sites and TAD boundaries.

```
Insulation Score (sliding window = 5)
  ▲
  │  ╭──╮      ╭───╮   ╭──╮
  │ ╭╯  ╰╮    ╭╯   ╰╮ ╭╯  ╰──
  │╭╯    ╰────╯     ╰─╯
  └─────────────────────────────▶ Bin
         ▼     ▼         ▼
     TAD boundary (valley)
```

Gaussian smoothing (`σ = 1.5`) removes noise before peak-finding.

### 3. Constrained Distance Matrix

Interactions are converted to physical distances using a polymer power law:

```
distance(i, j) = 1 / (log_interaction + ε)^0.7
```

Two biological constraints are then enforced:
- **Backbone stiffness:** `dist[i, i+1] = 0.8` for all adjacent bins
- **Contact pinching:** Top 1% off-diagonal interactions → `dist[i, j] = 1.5`

### 4. MDS 3D Embedding

scikit-learn's `MDS` is run on the constrained distance matrix to find a 3D configuration that best satisfies all pairwise distance constraints simultaneously:

```python
mds = MDS(n_components=3, dissimilarity="precomputed",
          random_state=42, n_init=10, max_iter=2000, eps=1e-9)
coords = mds.fit_transform(dist_matrix)
```

### 5. Polymer Physics Refinement

A two-pass smoothing + iterative contact attraction loop refines the structure:

```
Pass 1: Global smoothing (σ=1.5) — removes MDS artifacts
For 3 iterations:
  For each contact pair (i, j):
    If dist(i, j) > 2.0: pull toward each other (step = 0.1)
  Pass 2: Local smoothing (σ=0.5) — preserves pinch geometry
```

---

## 🖱 Viewer Controls

| Control | Action |
|---|---|
| **Left-click drag** | Rotate the structure |
| **Scroll wheel** | Zoom in / out |
| **Hover a bin** | Show TAD info and gene annotations |
| **Click a bin** | Pin the popup |
| **Auto Orbit** button | Toggle auto-rotation |
| **Bins** button | Toggle bead rendering |
| **Reset** button | Return camera to default view |
| **Start Pre-Render** button | Animate through bins, optionally export PNGs |
| **Batch Process Folder** button | Load a folder of JSONs and export PNG sequences |

---

## 📁 File Structure

```
chroma-3d/
│
├── Batchn2j.py            # Python pipeline: .npy → .json
├── index.html             # Browser viewer (no build step needed)
│
├── chromatin_data.json    # Active structure file (place your output here)
├── mouse.parquet          # Optional: gene annotation table
│
├── chr9/                  # Example input folder
│   ├── chr9_10000000_20000000.npy
│   └── chr9_20000000_30000000.npy
│
├── json_output/           # Generated by Batchn2j.py
│   ├── chr9_10000000_20000000.json
│   └── chr9_20000000_30000000.json
│
└── requirements.txt
```

---

## 🧬 Gene Annotation Support

The viewer optionally loads a **Parquet file** (`mouse.parquet`) containing gene annotations. When a user hovers over a bin, the viewer queries the gene table for any genes whose genomic range overlaps with the hovered TAD.

**Expected Parquet schema:**

| Column | Type | Description |
|---|---|---|
| `chr` | string | Chromosome (e.g. `"chr9"`) |
| `start` | int64 | Gene start position (bp) |
| `end` | int64 | Gene end position (bp) |
| `gene_id` | string | Ensembl gene ID |
| `gene_name` | string | HGNC gene symbol |

The file is parsed entirely in the browser using [`parquet-wasm`](https://github.com/kylebarron/parquet-wasm) and Apache Arrow — no server needed.

---

## 🎬 Batch PNG Export

The **Pre-Render** feature animates through every genomic bin sequentially, highlighting it in the 3D view. Each frame can be saved as a PNG directly to a user-selected local folder using the [File System Access API](https://developer.mozilla.org/en-US/docs/Web/API/File_System_Access_API).

This produces a frame sequence suitable for:
- Creating GIFs or videos with `ffmpeg`
- Scientific figures showing the full chromatin fold
- Comparative visualizations across conditions

```bash
# Assemble exported PNGs into a video with ffmpeg
ffmpeg -framerate 10 -i chromatin_data_bin_%04d.png \
       -c:v libx264 -pix_fmt yuv420p chromatin.mp4
```

> Batch export requires **Chrome or Edge** (Firefox does not support `showDirectoryPicker`).

---

## 📦 Dependencies

### Python

| Package | Purpose |
|---|---|
| `numpy` | Matrix operations and array math |
| `scipy` | Gaussian smoothing (`ndimage`), peak-finding (`signal`) |
| `scikit-learn` | Multidimensional Scaling (`MDS`) |
| `pathlib`, `re`, `json`, `os` | Standard library I/O |

### Browser (CDN — no install needed)

| Library | Version | Purpose |
|---|---|---|
| [Tailwind CSS](https://tailwindcss.com) | latest | UI utility classes |
| [parquet-wasm](https://github.com/kylebarron/parquet-wasm) | 0.6.1 | WebAssembly Parquet reader |
| [Apache Arrow](https://arrow.apache.org) | 16.0.0 | In-memory columnar data format |

---

## 🤝 Contributing

Contributions are welcome! Here are some areas where help would be valuable:

- **Alternative normalizations** — KR, ICE, VC normalization pre-processing
- **Additional TAD callers** — directionality index, HiCseg integration
- **WebGL renderer** — replace Canvas2D with a Three.js backend for larger structures
- **Comparative view** — side-by-side rendering of two conditions
- **Export formats** — `.pdb` or `.mmCIF` export for molecular visualization software

Please open an issue before submitting a large pull request.

---

## 📄 License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

---

<div align="center">

Made with 🧬 for computational genomics

</div>
