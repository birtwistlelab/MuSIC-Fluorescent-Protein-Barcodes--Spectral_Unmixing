# From the Step7, we know that 8.pR-mTFP1 has the lowest AUC among the 18 pR_fps. Here to confirm that 8.pR_mTFP1 is
# not a proper reference, we plot its FPR heatmap. Below is the condition1: original.

import numpy as np
import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'DejaVu Sans'
project_root = os.path.dirname(os.path.abspath(__file__))
result_path = os.path.join(project_root, 'output/')

optimal_threshold = np.load(result_path + '7.optimal_threshold_dict.npy', allow_pickle=True).item()
pos_pR_fp_scale_x = np.load(result_path + '5.pos_pR_fp_scale_x.npy', allow_pickle=True).item()
pos_cell_dict = np.load(result_path + '4.pos_cell_dict_pR_fp_singlets.npy', allow_pickle=True).item()
RF_peak = np.load(result_path + '2.RF_peak.npy', allow_pickle=True).item()

fp_name = ["EBFP2", "mTagBFP2", "mT-Sapphire", "mAmetrine", "mCerulean3", "LSSmOrange", "mBeRFP", "mTFP1", "EGFP",
           "CyOFP1", "mClover3", "mVenus", "mPapaya", "mOrange2", "mRuby3", "mKate2", "mCardinal", "miRFP670"]

threshold_list = []
for key, val in optimal_threshold.items():
    threshold_list.append(val)
threshold = np.array(threshold_list)

result_list = []
fpr_list = []
tpr_list = []

for key, val in sorted(pos_pR_fp_scale_x.items()):
    data = val[:, 1:]
    total_cells = len(data)
    print("total cells of " + key, total_cells)

    result = []
    cell_idx = []

    for each_key in RF_peak.keys():
        if each_key[-2:] in key:
            peak_channel = RF_peak[each_key][0]
            print('peak_channel', peak_channel)

    for each_cell in data:
        score = each_cell / threshold

        max_val = max(score)
        max_idx = np.argmax(score)

        most_likely_fp = fp_name[max_idx]
        score_list = [0] * 18
        score_list[max_idx] = 1
        result.append(score_list)
        cell_idx.append(max_idx)

    result = np.array(result)
    column_sums = np.sum(result, axis=0)
    percentage = np.round(column_sums / total_cells, decimals=4)
    result_list.append(percentage)

    # find the TPR
    tpr = np.max(percentage)
    fpr = 1 - tpr
    tpr_list.append(tpr)
    fpr_list.append(fpr)

# to create the 18x18 FPR heatmap
result_array = np.array(result_list)
np.fill_diagonal(result_array, np.nan)
df = pd.DataFrame(result_array, columns=[f'Reference {i + 1}' for i in range(18)],
                  index=[f'Reference {i + 1}' for i in range(18)])

fig, ax = plt.subplots(figsize=(11, 11))
sns.heatmap(df, annot=False, cmap='Blues', mask=df.isnull(), cbar=True, linecolor='white', linewidths=1,
            xticklabels=fp_name, yticklabels=fp_name, vmax=0.2)

ax.set_aspect('equal')
for i in range(len(df)):
    ax.add_patch(
        plt.Rectangle((i + 0.05, i + 0.05), 0.9, 0.9, fill=True, facecolor='gray', alpha=0.5, edgecolor='white', lw=1))

colorbar = ax.collections[0].colorbar
colorbar.ax.tick_params(labelsize=18)

plt.ylabel('Actual', fontsize=22)
plt.xlabel('Inferred', fontsize=22)
plt.yticks(fontsize=18)
plt.xticks(fontsize=18)
plt.tight_layout()
plt.subplots_adjust(bottom=0.4)

file_path = result_path + '8.original_condition1/8_1_FPR/'
os.makedirs(file_path, exist_ok=True)
filename = '8_1.heatmap of FPR(condition1:original).png'
plt.savefig(file_path + filename, transparent=True)

fig_path = os.path.join(project_root, 'paper_fig/')
filename1 = 'Fig3C(Step8).png'
plt.savefig(fig_path + filename1, transparent=True)
plt.close()

th = pd.DataFrame(result_array)
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
print(th)

th.to_excel(file_path + '8.condition1_original_fpr_result_18x18.xlsx', index=False)

np.save(result_path + '8_1.condition1_original_FPR.npy', fpr_list)
np.save(result_path + '8_1.condition1_original_TPR.npy', tpr_list)
