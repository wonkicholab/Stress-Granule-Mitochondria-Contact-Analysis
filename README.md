# Stress-Granule-Mitochondria-Contact-Analysis
**Won-Ki Cho Lab, KAIST**

---

## Abstract  
This repository contains codes used for quantifying interactions between stress granules (SGs) and mitochondria in live-cell imaging data.
By extending Fiji’s TrackMate with a custom merge-cost threshold and combining Python-based distance calculations with MATLAB-based contact analysis,
our workflow yields detailed metrics of SG–mitochondria proximity before merging events.

Here, we introduce a modified TrackMate pipeline that incorporates spot radii into the linking cost, followed by automated distance and contact analysis.

---

## 2. Steps

### 2.1 Spot Tracking in Fiji/TrackMate  
- **Image data**: Time-lapse stacks of SG-labeled cells.  
- **Plugin**: TrackMate vX.X (Tinevez _et al._, 2017).  
- **Modification**:  
  ```java
  // Only link if cost ≤ combined radii + sqrt(threshold)
  final double cost = mCostFunction.linkingCost(source, target);
  final double sr   = source.getFeature("RADIUS");
  final double tr   = target.getFeature("RADIUS");
  if (Math.sqrt(cost) > (sr + tr + Math.sqrt(mCostThreshold))) {
      continue;
  }
  ```  
- **Outputs**:  
  - `data/merging_G.xml`  (graph of nodes and edges)  
  - `data/branch.csv`     (flags for predecessor/successor branches)  
  - `data/rois.csv`       (spot ID ↔ ROI pixel lists)  
  - `data/mitomasks/*.bmp` (binary mitochondria masks per frame)

### 2.2 Distance Calculation (Python)  
**Script**: `scripts/spot_distance_calculation.py`  
- **Dependencies**: Python 3.8+, pandas, numpy, scikit-learn, imageio, OpenCV.  
- **Procedure**:  
  1. Load and clean `spots.csv` (TrackMate output).  
  2. For each spot, read its ROI pixel coordinates from `rois.csv`.  
  3. Load the corresponding frame’s mitochondrial mask (fallback to previous frame if empty).  
  4. Construct a KD-Tree of mitochondria pixel coordinates; query each spot’s pixels to obtain the minimum Euclidean distance.  
- **Output**: `results/spot_distances.xlsx` (columns: Spot ID, Distance)

### 2.3 Contact Analysis (MATLAB)  
**Script**: `scripts/contact_analysis_using_branch_for_merging_track.m`  
- **Dependencies**: MATLAB R2020b+, Image Processing Toolbox.  
- **Procedure**:  
  1. Import `spot_distances.xlsx`, `merging_G.xml`, and `branch.csv`.  
  2. Use `shortestpath` on the graph to extract:  
     - pre-merge branch nodes  
     - post-merge branch nodes  
     - full branch nodes  
  3. For each branch segment, compute:  
     - **Mean distance**  
     - **Contact frequency** (fraction of frames where distance = 0)  
     - **Average contact duration** (mean length of consecutive zero-distance runs)  
  4. Export summaries as Excel files in `results/`:  
     - `before_merge_summary.xlsx`  
  5. Print overall averages to the MATLAB console.

---

## References  
1. Tinevez J.-Y. _et al._ (2017) _TrackMate: An open and extensible platform for single-particle tracking_. **Methods**.  
2. Pedregosa F. _et al._ (2011) _Scikit-learn: Machine Learning in Python_. **JMLR**.  
3. MATLAB® (2020b), The MathWorks, Inc., Natick, MA, USA.  
