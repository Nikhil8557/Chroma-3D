import numpy as np
import json
import os
from pathlib import Path
from sklearn.manifold import MDS
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
import re

def generate_chromatin(input_file, output_dir="json_output"):
    # ==========================================
    # 0. DATA LOADING & LOG-NORMALIZATION
    # ==========================================
    filename = Path(input_file).name
    chromosome, start, end = re.search(
    r'(chr\w+)_(\d+)_(\d+)',
    filename).groups()
    start1 = int(start)
    end1 = int(end)

    if os.path.exists(input_file):
        heatmap = np.load(input_file)
        N = heatmap.shape[0]
    else:
        print("File not found")
        N = 200
        heatmap = np.eye(N) * 10.0
        # Simulated backbone
        for i in range(N-1): heatmap[i, i+1] = heatmap[i+1, i] = 15
        # Simulated biological loops (Enhancer-Promoter)
        heatmap[30, 70] = heatmap[70, 30] = 40
        heatmap[110, 160] = heatmap[160, 110] = 45

    # Log-transform interactions (Hi-C data is naturally log-normal)
    # This prevents extreme interaction peaks from "crushing" the model
    heatmap_log = np.log1p(heatmap)

    # ==========================================
    # 1. BIOLOGICAL TAD DETECTION (Insulation Score)
    # ==========================================
    # Realism: We don't hardcode TADs. We find them using the Insulation Score.
    window = 5
    insulation = []
    for i in range(N):
        s, e = max(0, i-window), min(N, i+window)
        # Sum local interactions (TAD density)
        insulation.append(np.mean(heatmap_log[s:e, s:e]))
    
    # Valleys in insulation score represent TAD boundaries (CTCF sites)
    valleys, _ = find_peaks(-gaussian_filter1d(insulation, 1.5), distance=10)
    boundaries = np.unique(np.sort(np.concatenate(([0], valleys, [N-1])))).astype(int)

    # ==========================================
    # 2. CONSTRAINED DISTANCE MATRIX
    # ==========================================
    # Convert interactions to physical distances using a 1/3 power law
    # (Typical for polymer volume-to-distance relationships)
    dist_matrix = 1.0 / (heatmap_log + 1e-6)**0.7

    # FORCE BACKBONE STIFFNESS: 
    # Bins i and i+1 must be exactly 1 unit apart to represent a continuous chain.
    for i in range(N-1):
        dist_matrix[i, i+1] = dist_matrix[i+1, i] = 0.8 

    # IDENTIFY PINCH POINTS:
    # Find significant off-diagonal interactions
    contacts = []
    threshold = np.percentile(heatmap_log, 99.0)
    for i in range(N):
        for j in range(i + 10, N):
            if heatmap_log[i, j] > threshold:
                dist_matrix[i, j] = dist_matrix[j, i] = 1.5 # Contact distance
                contacts.append((i, j))

    # ==========================================
    # 3. GLOBAL FOLDING (MDS)
    # ==========================================
    print(f"Folding polymer with {len(contacts)} biological contacts...")
    mds = MDS(n_components=3, dissimilarity="precomputed", random_state=42, 
              n_init=10, max_iter=2000, eps=1e-9)
    coords = mds.fit_transform(dist_matrix)

    # ==========================================
    # 4. POLYMER PHYSICS REFINEMENT
    # ==========================================
    # Improvement: Use a double-smoothing pass.
    # Pass 1: Global stiffness (high sigma)
    coords = gaussian_filter1d(coords, sigma=1.5, axis=0)
    
    # Pass 2: Local refinement (low sigma) to preserve the 'pinches'
    for _ in range(3): # Iterative refinement
        for i, j in contacts:
            # Physically pull contact points toward each other
            vec = coords[j] - coords[i]
            dist = np.linalg.norm(vec)
            if dist > 2.0: # If too far, pull them
                pull = vec * 0.1
                coords[i] += pull
                coords[j] -= pull
        coords = gaussian_filter1d(coords, sigma=0.5, axis=0)

    # ==========================================
    # 5. SCALE & EXPORT
    # ==========================================
    coords -= np.mean(coords, axis=0)
    # Fit into a 300-unit box for the web viewer
    scale = 280 / (np.max(np.abs(coords)) or 1)
    coords *= scale
    

    # Color by TAD IDs
    palette = ["#3b82f6", "#8b5cf6", "#ec4899", "#f59e0b", "#10b981", "#ef4444", "#06b6d4"]
    colors = []
    for i in range(N):
        tad_idx = np.searchsorted(boundaries, i) % len(palette)
        colors.append(palette[tad_idx])

    output = {
        "x": coords[:, 0].tolist(),
        "y": coords[:, 1].tolist(),
        "z": coords[:, 2].tolist(),
        "Chromosom": chromosome,
        "start": start1,
        "end": end1,
        "colors": colors,
        "tad_boundaries": boundaries.tolist(),
    }

    os.makedirs(output_dir, exist_ok=True)

        # Output filename
    output_name = Path(filename).stem + ".json"
    output_path = os.path.join(output_dir, output_name)

    with open(output_path, "w") as f:
        json.dump(output, f)

    print(f"Saved: {output_path}")
    print(f"Realism Success! {len(boundaries)-1} TADs detected. Exported {N} bins.")

def process_batch(input_dir="./", output_dir="json_output"):
    npy_files = sorted(Path(input_dir).glob("*.npy"))

    if not npy_files:
        print("No .npy files found.")
        return

    print(f"Found {len(npy_files)} files.\n")

    for file_path in npy_files:
        try:
             generate_chromatin(str(file_path), output_dir)
        except Exception as e:
            print(f"Failed processing {file_path.name}: {e}")


if __name__ == "__main__":
    process_batch(
        input_dir="./chr9",
        output_dir="./json_output"
    )








