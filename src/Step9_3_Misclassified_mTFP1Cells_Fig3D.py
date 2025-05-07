import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import math
from matplotlib.ticker import FuncFormatter
from my_module import scientific_formatter

plt.rcParams['font.family'] = 'Arial'
current_path = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_path)
result_path = os.path.join(project_root, 'output/')

loading_cell_path = [os.path.join(result_path, '4_2.SeparateTrainingAndTestingData/', f'test{i + 1}/') for i in range(3)]
loading_excel_path = [os.path.join(result_path, '9_1.RangeIntensityCutoff/', f'test{i + 1}/') ### JRH Changed from 'CutOff' to 'Cutoff'
                      for i in range(3)]

negative_thresholds = np.load(os.path.join(result_path, '4_1.negative_control_thresholds.npy'))
saving_path = os.path.join(result_path, '9_3.Misclassified_mTFP1Cells_Fig3D')
os.makedirs(saving_path, exist_ok=True)

total_misclassified_mTFP1_FI, total_misclassified_other_FI = [], []
for i in range(3):
    testing_pool = np.load(loading_cell_path[i] + 'testing_pool.npy')
    total_cells_original = len(testing_pool)
    df = pd.read_excel(os.path.join(loading_excel_path[i], 'cutoff_1000.xlsx'))
    total_cells_cutoff1000 = df.shape[0]
    if total_cells_original == total_cells_cutoff1000:
        print('original data == cutoff_1000 data')

    intensity_mTFP1_val, intensity_other_val = [], []
    for idx, row in df.iterrows():
        if row['actual_fp'] == '08.mTFP1':
            if row['predicted_fp'] != '08.mTFP1':
                intensity_mTFP1_val.append(float(row['fluorescence_intensity_peak_channel']))
                total_misclassified_mTFP1_FI.append(float(row['fluorescence_intensity_peak_channel']))
        else:
            if row['predicted_fp'] != row['actual_fp']:
                intensity_other_val.append(float(row['fluorescence_intensity_peak_channel']))
                total_misclassified_other_FI.append(float(row['fluorescence_intensity_peak_channel']))
    print(len(intensity_mTFP1_val), len(total_misclassified_mTFP1_FI))
    print(len(intensity_other_val), len(total_misclassified_other_FI))

# to create False Positive(FP) FP_cell histogram
log_value_fp = []
for each in total_misclassified_other_FI:
    log_value = math.log10(each)
    log_value_fp.append(round(log_value, 2))

log_value_fp = np.array(log_value_fp)

bins_num = round((max(log_value_fp) - min(log_value_fp)) * 20)

x_ticks = np.logspace(np.log10(min(total_misclassified_other_FI)), np.log10(max(total_misclassified_other_FI)), num=20)

# bin_width is an array containing each bin's right edge value
bin_width = np.logspace(np.log10(min(total_misclassified_other_FI)), np.log10(max(total_misclassified_other_FI)), num=bins_num)

negative_threshold = min(negative_thresholds)
print('negative_threshold', negative_threshold)
print('min value of misclassified mTFP1 cells: ', min(total_misclassified_mTFP1_FI))
print('min value of misclassified other cells: ', min(total_misclassified_other_FI))

fig, ax = plt.subplots(figsize=(6, 4))
ax.hist(total_misclassified_mTFP1_FI, bins=bin_width, edgecolor='navy', facecolor='none', linewidth=1.5, alpha=0.6,
        label='mTFP1', zorder=5)
ax.hist(total_misclassified_other_FI, bins=bin_width, edgecolor='orange', facecolor='none', linewidth=1.5, alpha=0.6,
        label='All Other Fluorescent Proteins')

plt.axvline(x=negative_threshold, ymin=0, ymax=1, color='darkgray', linestyle='--', label='Transfected Cell Cutoff')
ax.set_xlabel('Maximum Fluorescence Intensity of Testing Cells', fontsize=14)
ax.set_ylabel('Misclassified Cell Count', fontsize=14)
ax.set_xlim(10 ** 2, 10 ** 7)
ax.set_ylim(0, 300)
plt.xscale('log')
plt.gca().xaxis.set_major_formatter(FuncFormatter(scientific_formatter))
ax.tick_params(axis='x', labelsize=12)
ax.tick_params(axis='y', labelsize=12)
plt.tight_layout()
plt.subplots_adjust(bottom=0.2)
plt.legend(loc='upper left', bbox_to_anchor=(0.42, 1), fontsize=12, framealpha=0, frameon=False,  handlelength=1.5, handleheight=1, handletextpad=0.3, labelspacing=0.2)

filename = 'histogram of misclassified mTFP1 cells (for fig3D).png'
plt.savefig(os.path.join(saving_path, filename), transparent=True)

filename2 = 'Fig.3D(Step9_3)misclassified mTFP1 cells.png'
plt.savefig(os.path.join(project_root, 'figures_for_paper/', filename2), dpi=300, transparent=True)
plt.close()
