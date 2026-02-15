# MIHC Breast Cancer Spatial Network Analysis

This repository contains data, scripts, notebooks, and figure outputs for graph-based analysis of MIBI (Multiplexed Ion Beam Imaging) breast cancer samples. The core idea is to convert spatially localized single-cell detections into proximity graphs and compute tumor-immune mixing and stromal barrier properties.

## Project goals

- Build cell-cell spatial graphs from per-patient MIBI detections.
- Standardize immune/tumor/stromal cell type labels.
- Compute network-derived phenotypes (Cold, Compartmentalized, Mixed).
- Analyze neighborhood affinity around tumor nodes.
- Generate publication-style visualizations for whole networks, subnetworks, degree distributions, and phenotype-level summaries.

## Repository contents (comprehensive)

### Root-level analysis code

- `net-viz.py`
  - Cleans per-cell detection tables.
  - Normalizes class labels using biologically meaningful mappings (e.g., CD3+CD4 → Helper T cell; CD68 → Macrophage).
  - Splits overlapping marker labels into multiple rows where needed.
  - Builds undirected cell-cell graphs by connecting cells whose Euclidean distance is `< 35` µm.
  - Writes graphs to `.gml`.

- `patient_network_properties.py`
  - Computes **mixing score** = immune-to-tumor edges / immune-to-immune edges.
  - Classifies samples into:
    - `Cold` if immune cell count < 250,
    - `Compartmentalized` if mixing score < 0.22,
    - `Mixed` otherwise.
  - Computes stromal clustering coefficient for `Other` (stroma-like) nodes.
  - Computes stromal barrier metrics as shortest-path-based counts of `Other` nodes between immune types and tumor cells.
  - Exports a patient-level CSV of metrics.

### Root-level notebooks

- `network_viz.ipynb`
  - Interactive workflow for turning raw MIBI detections into a spatial graph.
  - Includes class cleanup / overlap handling and graph construction logic.

- `imaging_subset.ipynb`
  - Exploratory analysis of selected high/low tumor regions.
  - Computes affinity matrices and produces boxplots/network visual diagnostics.

- `patient_centrality.ipynb`
  - Centrality and phenotype-oriented analyses (mixing score, clustering, stromal barrier, degree-related exploration).

- `patient_net_analysis.ipynb`
  - Additional patient-level graph analysis and subtype-specific tumor neighborhood exploration.

### Data directories

- `MIBI Image Data/` (41 CSVs)
  - Per-patient detection tables with columns:
    - `Object ID, Name, Class, Parent, Centroid X µm, Centroid Y µm`
  - Appears to be the full cohort input set.

- `data/Detections.csv`
  - Rich, feature-heavy detections file containing full nucleus/cell/cytoplasm marker measurements plus centroid coordinates.

- `patients_data/` (17 CSVs)
  - Processed or selected per-patient detection CSVs.

- `rem_patients_data/` (24 CSVs)
  - Additional/remainder per-patient detection CSVs used in network generation.

### Graph/network directories

- `patients_gml/` (17 GML files)
  - Patient-level graphs.

- `rem_patients_gml/` (24 GML files)
  - Additional generated patient graphs (includes files with `PP*` naming in this snapshot).

- `whole-network.gml`
  - Combined/global graph representation.

- `whole-network.dot`
  - GraphViz DOT representation of the whole network.

- `subgraphs/`
  - Phenotype-stratified subgraphs:
    - `subgraphs/Cold/` (7 GMLs)
    - `subgraphs/Comp/` (7 GMLs)
    - `subgraphs/Mixed/` (7 GMLs)
  - Each phenotype folder includes per-cell-type subnetworks (e.g., B cell, T cell, Macrophage, Regulatory T cell).

- `patient_fig/tumor_subgraphs/`
  - Tumor-focused subgraphs as `.gml` files for phenotype exemplars.

### Visualization/artifact directories

- `MIBI_Network_Viz/` (41 PNGs)
  - Per-patient rendered network figures.

- `figures/` (17 files)
  - Summary and layout figures (spring/kamada-kawai layouts, louvain clustering, degree distributions, tumor subnetwork panels, etc.).

- `degree_dist/` (3 PNGs)
  - Degree distribution visualizations by phenotype grouping.

- `subnetgraph/` (7 PNGs)
  - Rendered subnetwork images by cell type.

- `patient_fig/`
  - Organizes phenotype-specific patient visuals:
    - `cold_patient/`
    - `mixed_patient/`
    - `compartmentalized_patient/`
  - Contains spatial maps, neighborhood distributions, and tumor-focused plots.

- `gif/` (55 PNG frames) + `image.gif`
  - Frame sequence and compiled animation.

- `tumor_network.png`, `degree-dist.png`
  - Additional key standalone summary plots.

### Miscellaneous

- `coords`
  - Serialized binary artifact containing coordinate-like data (not plaintext; likely intermediate/cache output).

## Data schema notes

### Detection CSVs

Most patient-level detection files use this compact schema:

```text
Object ID,Name,Class,Parent,Centroid X µm,Centroid Y µm
```

The scripts rely on:

- `Class` (cell marker/classification label),
- fallback to `Parent` when `Class` is missing,
- centroid coordinates for graph construction.

### GML graph schema

Typical node attributes include:

- `pos` (x/y coordinates)
- `type` (cell type, e.g., Tumor, B cell, T cell, Macrophage, Other, Regulatory T cell)

Edges represent spatial adjacency under the chosen distance threshold.

## Reproducing the core pipeline

> The project is research-code style and does not currently include a locked dependency file.

### 1) Set up Python environment

Recommended packages (inferred from scripts and notebooks):

- `pandas`
- `numpy`
- `networkx`
- `matplotlib`
- `seaborn`

Example:

```bash
python -m venv .venv
source .venv/bin/activate
pip install pandas numpy networkx matplotlib seaborn jupyter
```

### 2) Generate patient graph(s)

- Edit `net-viz.py` input folder/selection logic as needed.
- Run:

```bash
python net-viz.py
```

This writes `.gml` output(s) in the configured destination (currently `rem_patients_gml/` in script defaults).

### 3) Compute patient metrics

For a specific graph:

```bash
python patient_network_properties.py patients_gml/P06.gml
```

Expected output: a per-patient CSV containing mixing score, phenotype label, stromal clustering, and stromal barrier metrics.

### 4) Explore and reproduce figures

Open notebooks for detailed exploratory and plotting workflows:

```bash
jupyter notebook
```

Then run:
- `network_viz.ipynb`
- `imaging_subset.ipynb`
- `patient_centrality.ipynb`
- `patient_net_analysis.ipynb`

## Method summary

1. Ingest per-cell detections with centroid coordinates.
2. Normalize and expand marker-defined classes to canonical immune/tumor/stromal labels.
3. Create a spatial proximity graph with distance-threshold edges.
4. Derive phenotype and centrality/network metrics.
5. Compare phenotypes and cell-type subnetworks through visualizations.

## Known caveats and implementation details

- Some scripts are configured for one-off execution paths (e.g., `net-viz.py` currently processes only `Detections_OP_P5.csv` in its main block).
- No formal packaging/test harness is included.
- The `coords` file is binary and undocumented in-repo.
- A few graph filenames in `rem_patients_gml/` use `PP*` naming; this may be intentional legacy naming or an inconsistency.

## Suggested next improvements

- Add `requirements.txt` or `environment.yml`.
- Convert scripts into a parameterized CLI (input dir, output dir, threshold, patient filters).
- Add unit tests for label mapping, distance-threshold graph construction, and phenotype classification.
- Add provenance metadata (data version/source and preprocessing history).

## Quick directory map

```text
.
├── MIBI Image Data/         # full per-patient detections (CSV)
├── patients_data/           # selected/processed detections (CSV)
├── rem_patients_data/       # remainder detections (CSV)
├── patients_gml/            # generated patient graphs
├── rem_patients_gml/        # additional generated patient graphs
├── subgraphs/               # phenotype + cell-type subgraphs
├── MIBI_Network_Viz/        # per-patient rendered network PNGs
├── figures/                 # aggregate manuscript-style figures
├── degree_dist/             # degree distribution plots
├── subnetgraph/             # subnetwork plots by cell type
├── patient_fig/             # phenotype-specific patient visual outputs
├── gif/                     # animation frames
├── image.gif                # compiled animation
├── net-viz.py               # graph generation script
├── patient_network_properties.py  # phenotype/metric extraction script
├── *.ipynb                  # exploratory notebooks
└── whole-network.{gml,dot}  # whole-network exports
```

---

If you want, this README can be further tailored into a publication supplement format (Methods, Results, Figure legend mapping), or split into `README` + `docs/` pages for reproducibility.
