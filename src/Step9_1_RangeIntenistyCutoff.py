# ====================================== XL Clemson ======================================
# First, range of intensity cutoffs applied to the testing pool
# ====================================== XL Clemson ======================================

import numpy as np
import os
import pandas as pd
from scipy.optimize import nnls
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'Arial'

current_path = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_path)
result_path = os.path.join(project_root, 'output/')

optimal_thresholds = np.load(os.path.join(result_path, '6.DetermineUnmixingThresholdsOnTrainingCells_FigS4/average_optimal_threshold.npy'))
thresholds = np.array(optimal_thresholds)

loading_cell_path = [os.path.join(result_path, '4_2.SeparateTrainingAndTestingData/', f'test{i + 1}/') for i in range(3)]

saving_path = [os.path.join(result_path, '9_1.RangeIntensityCutoff/', f'test{i + 1}/') for i in range(3)]
for path in saving_path:
    os.makedirs(path, exist_ok=True)

# Load required data files：
RF = np.load(result_path + '2.MFI_RF.npy')
RF_peak = np.load(result_path + '2.RF_peak.npy', allow_pickle=True).item()

fp_name = ['01.EBFP2', '02.mTagBFP2', '03.mT_Sapphire', '04.mAmetrine', '05.mCerulean3', '06.LSSmOrange', '07.mBeRFP',
           '08.mTFP1', '09.EGFP', '10.CyOFP1', '11.mClover3', '12.mVenus', '13.mPapaya', '14.mOrange2', '15.mRuby3',
           '16.mKate2', '17.mCardinal', '18.miRFP670']

# Transpose RF once before looping, or it may cause errors
RF = np.array(RF).T
x_label = 'Unmixed Relative Abundance'
y_label = 'Frequency'

# Set the intensity cutoff thresholds
cutoff_thresholds1 = np.arange(1000, 25000, 1000)
cutoff_thresholds2 = np.arange(25000, 50001, 2500)
cutoff_thresholds = np.concatenate((cutoff_thresholds1, cutoff_thresholds2))

np.save(os.path.join(result_path, '9_1.cutoff_thresholds.npy'), cutoff_thresholds)


for i in range(3):
    print('Test' + str(i + 1))
    # Unmix the testing_pool using NNLS
    testing_pool = np.load(loading_cell_path[i] + 'testing_pool.npy')  # testing_pool is a list
    scale_x = []
    unmixed_scale_x = []
    original_cells = []
    for idx, each_cell in enumerate(testing_pool):
        x, residuals = nnls(RF, each_cell)
        x = np.round(x, decimals=4)
        scale_x.append([idx, x[1:]]) # Keep only the 18 reference spectra coefficients (excluding x[0] which is the
        # reference for autofluorescence)
        unmixed_scale_x.append(x[1:])
        original_cells.append(each_cell)

    # Create a dataframe to store unmixed results
    df = pd.DataFrame({
        'idx': range(len(unmixed_scale_x)),
        'original_array': original_cells,
        'unmixed_array': unmixed_scale_x,
    })

    # Assign the true FP labels to each cell in the testing pool
    testing_dict = np.load(loading_cell_path[i] + 'testing_dict.npy', allow_pickle=True).item()
    fp_label = []
    true_fp_label = {}
    for fp, cells in testing_dict.items():
        for idx, test_cell in enumerate(testing_pool):
            if np.any(np.all(np.isclose(cells, test_cell, atol = 1e-6), axis=1)):
                fp_label.append([idx, fp, test_cell])
                true_fp_label[idx] = fp
    df['actual_fp'] = df['idx'].map(true_fp_label)

    # Predict FP label for each cell by normalizing and selecting the maximum score
    predicted_fp = []
    scale_x_array = np.array(unmixed_scale_x)
    df.insert(df.columns.get_loc('actual_fp'), 'normalized_array', (scale_x_array / thresholds).tolist())

    for score in df['normalized_array']:
        max_val = max(score)  # Get the maximum normalized value
        max_idx = np.argmax(score)  # Identify the index of the most-likely FP
        most_likely_fp = fp_name[max_idx]
        predicted_fp.append(most_likely_fp)
    df['predicted_fp'] = np.array(predicted_fp)

    # Retrieve fluorescence intensity for each FP-related peak channel
    peak_FI_list = []
    for idx, row in df.iterrows():
        for key in RF_peak.keys(): # key: 01.EBFP2
            if df['actual_fp'][idx][:2] == key[-2:]:
                peak_channel = RF_peak[key][0]
                peak_FI = df['original_array'][idx][peak_channel]
                peak_FI_list.append(peak_FI)
    df['fluorescence_intensity_peak_channel'] = peak_FI_list

    col_index = df.columns.get_loc('original_array')
    df.insert(col_index + 1, 'fluorescence_intensity_peak_channel', df.pop('fluorescence_intensity_peak_channel'))

    df.to_excel(os.path.join(saving_path[i], 'identification_result_original.xlsx'), index=False)

    # Apply intensity cutoff to filter cells based on their fluorescence intensity
    cutoff_cell_id_dict = {}
    for each_cutoff in cutoff_thresholds:
        filtered_df = df[df['fluorescence_intensity_peak_channel'] >= each_cutoff]
        if not filtered_df.empty:
            file_path = os.path.join(saving_path[i], f'cutoff_{each_cutoff}.xlsx')
            filtered_df.to_excel(file_path, index=False)
