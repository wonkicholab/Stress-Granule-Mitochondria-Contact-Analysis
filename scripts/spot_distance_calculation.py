#!/usr/bin/env python3
# Compute the minimum distance between each stress granule spot and mitochondria.

import pandas as pd
import numpy as np
import re
import imageio.v2 as ii
import cv2 as cv
from sklearn.neighbors import KDTree

# Prompt user for data folder containing input files.
folder = input("Enter the folder: ")

# Load TrackMate spots.csv and drop header rows
spots = pd.read_csv(folder+'/spots.csv', encoding='unicode_escape', low_memory=False)
spots = spots.drop([0,1,2], axis=0)
# Convert FRAME and TRACK_ID columns to numeric
dtype_cols = ['FRAME', 'TRACK_ID']
for col in dtype_cols:
    spots[col] = spots[col].apply(pd.to_numeric)

# Load ROI mappings for stress granule spots
rois_file = pd.read_csv(folder+"/rois.csv")
rois_file['Index'] = rois_file['Index'].astype(str)

# Prepare lists to store results
spot_id = []
spot_d  = []

# Iterate over each spot entry
for idx in range(len(spots)):
    # Skip entries without a valid TRACK_ID
    if np.isnan(spots['TRACK_ID'].iloc[idx]):
        continue
    # Get the frame number for this spot
    time = spots['FRAME'].iloc[idx]

    # Identify ROI CSV file for this spot
    target_label = spots['LABEL'].iloc[idx]
    sel = rois_file['Name'] == target_label
    roi_csv = f"{folder}/{rois_file.loc[sel, 'Index'].values[0]}.csv"
    pixels = pd.read_csv(roi_csv)[['X','Y']].to_numpy()

    # Construct the BMP filename for the mitochondrial mask
    frame_str = f"{int(time):04d}.bmp"
    mask_path = f"{folder}/{frame_str}"

    # Read the mitochondrial mask image
    MT_image = ii.imread(mask_path)
    y_coords, x_coords = np.where(MT_image > 0)

    # If no mitochondria pixels found, search previous frames
    while x_coords.size == 0:
        time -= 1
        frame_str = f"{int(time):04d}.bmp"
        MT_image = ii.imread(f"{folder}/{frame_str}")
        y_coords, x_coords = np.where(MT_image > 0)

    # Build KDTree for mitochondria pixels and query distances
tree = KDTree(np.column_stack((x_coords, y_coords)))
dists, _ = tree.query(pixels, k=1)
min_dist = dists.min()

    # Clean spot label and store results
    spot_id.append(int(re.sub(r'[^0-9]', '', str(target_label))))
    spot_d.append(min_dist)

# Save results to Excel
result_df = pd.DataFrame({'SpotID': spot_id, 'Distance': spot_d})
result_df.to_excel(folder+'/result_whole_spot_distance.xlsx', index=False)
