# After removing mTFP1, find the optimal threshold for each pR_fp by identifying the threshold point on the ROC curve
# that balances higher TPR and lower FPR which is the one closest to the ideal point(0,1) as (FPR, TPR).

import numpy as np
from sklearn.metrics import auc
import matplotlib.pyplot as plt
import os

plt.rcParams['font.family'] = 'DejaVu Sans'
project_root = os.path.dirname(os.path.abspath(__file__))
result_path = os.path.join(project_root, 'output/')

RF_name = np.load(result_path + '2.RF_name.npy')
RF_name_17 = np.load(result_path + '10_1.RF_name_17.npy')
x_axis = np.load(result_path + '1.x_axis_channel.npy')
channel_num = len(x_axis)
pos_pR_fp_scale_x_17 = np.load(result_path + '10_1.pos_pR_fp_scale_x_17.npy', allow_pickle=True).item()

thresholds_1 = np.arange(0, 1, 0.01)
thresholds_2 = np.arange(1, 10, 0.1)
thresholds = np.concatenate((thresholds_1, thresholds_2))

count = -1
result_dict_17 = {}
auc_dict_17 = {}

# to create the 17 ROC plots
for key, val in pos_pR_fp_scale_x_17.items():
    print(key)
    count += 1
    # to skip the RF0
    new_val = val[:, 1:]
    result_list = []
    TPR = []
    FPR = []
    for each_threshold in thresholds:
        true_posCell_list = []
        false_negCell_list = []
        true_negCell_list = []
        false_posCell_list = []
        for col in range(new_val.shape[1]):
            RF_array = new_val[:, col]
            for idx in range(len(RF_array)):
                if col == count:
                    if RF_array[idx] < each_threshold:
                        false_negCell_list.append(idx)
                    else:
                        true_posCell_list.append(idx)
                else:
                    if RF_array[idx] < each_threshold:
                        true_negCell_list.append(idx)
                    else:
                        false_posCell_list.append(idx)

        result = [len(true_posCell_list), len(false_negCell_list), len(true_negCell_list), len(false_posCell_list)]

        tpr = result[0] / (result[0] + result[1])
        fpr = result[3] / (result[2] + result[3])
        TPR.append(tpr)
        FPR.append(fpr)
        result_list.append([fpr, tpr])
    roc_auc = auc(np.array(FPR), np.array(TPR))
    auc_dict_17[key] = roc_auc

    plt.figure(figsize=(6, 4))
    plt.plot(FPR, TPR, label=f'{key} (AUC = {roc_auc:.3f})')
    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.xlim([-0.01, 1.05])
    plt.ylim([-0.01, 1.05])
    plt.legend(loc='lower right')

    file_path = result_path + '10.condition3_remove_mTFP1/10_2.ROC_folder/'
    os.makedirs(file_path, exist_ok=True)
    filename = key + '_ROC.png'
    plt.savefig(file_path + filename)
    plt.close()

    result_dict_17[key] = np.array(result_list)

# to create a summary plot with 17 ROC curves, each ROC curve for one pR-fp
plt.figure(figsize=(12, 8))
for key, result_array in result_dict_17.items():
    FPR = result_array[:, 0]
    TPR = result_array[:, 1]
    roc_auc = auc(FPR, TPR)
    plt.plot(FPR, TPR, label=f'{key} (AUC = {roc_auc:.3f})')

plt.plot([0, 1], [0, 1], 'k--')
plt.xlabel('False Positive Rate', fontsize=18)
plt.ylabel('True Positive Rate', fontsize=18)
plt.xlim([-0.01, 1.05])
plt.ylim([-0.01, 1.05])
plt.legend(loc='lower right')

file_path2 = result_path + '10.condition3_remove_mTFP1/10_2.ROC_summary/'
os.makedirs(file_path2, exist_ok=True)
filename2 = 'All_17Plasmids_ROC.png'
plt.savefig(file_path2 + filename2, transparent=True)
plt.close()

# to create 3x6 subplots with 3 decimals of each AUC

fig, axs = plt.subplots(3, 6, figsize=(20, 11))

data_list = []
for key, result_array in result_dict_17.items():
    FPR = result_array[:, 0]
    TPR = result_array[:, 1]
    roc_auc = auc(FPR, TPR)
    data_list.append([key[3:], FPR, TPR, roc_auc])

for i in range(3):
    for j in range(6):
        ax = axs[i, j]
        ax_index = i * 6 + j
        if ax_index < len(data_list):
            data = data_list[ax_index]
            ax.plot(data[1], data[2], color='blue')
            ax.plot([], [], ' ', label=f'{data[0]}\n(AUC = {data[3]:.3f})')
            ax.plot([0, 1], [0, 1], 'k--')
            ax.set_xlim([-0.01, 1.05])
            ax.set_ylim([-0.01, 1.05])
            ax.legend(frameon=False, loc='lower right', fontsize=12)
            ax.set_aspect('equal', adjustable='box')

        if i == 1 and j == 0:
            ax.set_ylabel('True Positive Rate', fontsize=16)

        x_label_horizontal_pos = 0.5
        x_label_vertical_pos = 0.01
        fig.text(x_label_horizontal_pos, x_label_vertical_pos, 'False Positive Rate', ha='center', fontsize=16)

plt.subplots_adjust(hspace=0.02, wspace=0.02)
plt.tight_layout()

filename3 = '10_2_17Plasmids_ROC.png'
plt.savefig(file_path2 + filename3)
plt.close()

# to find the optimal thresholds for these 17 pR-fps:

optimal_threshold = {}
for key, val in result_dict_17.items():
    print(key)
    np.set_printoptions(suppress=True)
    FPR = np.array(val[:, 0])
    TPR = np.array(val[:, 1])
    roc_auc = auc(FPR, TPR)

    # ideal_point would be (0,1)
    distances = np.sqrt((FPR - 0) ** 2 + (TPR - 1) ** 2)

    optimal_idx = np.argmin(distances)
    optimal_FPR = FPR[optimal_idx]
    optimal_TPR = TPR[optimal_idx]

    print('the optimal_idx of pR_' + key + ' is ' + str(optimal_idx))

    print('the best threshold is ' + str(round(thresholds[optimal_idx], 2)))

    plt.plot(FPR, TPR, label=f'{key} (AUC = {roc_auc:.3f})')
    plt.scatter(optimal_FPR, optimal_TPR, color='red', label='Optimal Point', zorder=2)

    offset = 0.02
    plt.text(optimal_FPR + offset, optimal_TPR, f'({optimal_FPR:.2f}, {optimal_TPR:.3f})',
             fontsize=9, verticalalignment='top', horizontalalignment='left')

    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.legend(loc='lower right')

    file_path3 = result_path + '10.condition3_remove_mTFP1/10_2.ROC_threshold_folder/'
    os.makedirs(file_path3, exist_ok=True)
    filename4 = key + '_ROC.png'
    plt.savefig(file_path3 + filename4)
    plt.close()

    optimal_threshold[key] = round(thresholds[optimal_idx], 2)
    print(f'Optimal Point: FPR = {round(optimal_FPR, 2)}, TPR = {round(optimal_TPR, 2)}')

# to create the 3x6 subplots with 3 decimals of each AUC, add the best threshold

fig, axs = plt.subplots(3, 6, figsize=(20, 12))

data_list = []
pR_fp_lack_08mTFP1 = {}
for key, result_array in result_dict_17.items():
    FPR = result_array[:, 0]
    TPR = result_array[:, 1]
    roc_auc = auc(FPR, TPR)

    # ideal_point would be (0,1)
    distances = np.sqrt((FPR - 0) ** 2 + (TPR - 1) ** 2)

    optimal_idx = np.argmin(distances)
    optimal_FPR = FPR[optimal_idx]
    optimal_TPR = TPR[optimal_idx]

    data_list.append([key[3:], FPR, TPR, roc_auc, optimal_FPR, optimal_TPR])
    pR_fp_lack_08mTFP1[key] = [FPR, TPR, roc_auc]

for i in range(3):
    for j in range(6):
        ax = axs[i, j]
        ax_index = i * 6 + j
        if ax_index < len(data_list):
            data = data_list[ax_index]
            ax.plot(data[1], data[2], color='blue')
            ax.plot([], [], ' ', label=f'{data[0]}\n(AUC = {data[3]:.3f})')
            ax.plot([0, 1], [0, 1], 'k--')
            ax.scatter(data[4], data[5], color='red', zorder=2)
            offset = 0.02
            ax.text(data[4] + offset, data[5], f'({data[4]:.3f}, {data[5]:.3f})',
                    fontsize=10, verticalalignment='top', horizontalalignment='left')
            ax.set_xlim([-0.01, 1.05])
            ax.set_ylim([-0.01, 1.05])
            ax.legend(frameon=False, loc='lower right', fontsize=12)
            ax.set_aspect('equal', adjustable='box')
        else:
            ax.axis('off')

        if i == 1 and j == 0:
            ax.set_ylabel('True Positive Rate', fontsize=16)

        x_label_horizontal_pos = 0.5
        x_label_vertical_pos = 0.01
        fig.text(x_label_horizontal_pos, x_label_vertical_pos, 'False Positive Rate', ha='center', fontsize=16)

plt.subplots_adjust(hspace=0.02, wspace=0.02)
plt.tight_layout()
file_path4 = result_path + '10.condition3_remove_mTFP1/10_2.ROC_threshold_summary/'
os.makedirs(file_path4, exist_ok=True)
filename5 = '17_Plasmids_ROC.png'
plt.savefig(file_path4 + filename5, dpi=300, transparent=True)
plt.close()

np.save(result_path + '10_2.full_info_17fp.npy', pR_fp_lack_08mTFP1, allow_pickle=True)
np.save(result_path + '10_2.optimal_threshold_17dict.npy', optimal_threshold, allow_pickle=True)
np.save(result_path + '10_2.result_dict_for_17ROC.npy', result_dict_17, allow_pickle=True)
np.save(result_path + '10_2.auc_dict_for_17ROC.npy', auc_dict_17, allow_pickle=True)
