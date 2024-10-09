# from the 18x18 plots, we can get 1 True_Positive(TP), 1 False_Negative(FN), 17 True_negatives(TNs), and 17
# False_positives(FPs) to generate 18 ROCs (one for each pR_fp), and 18 in 1 ROC for better visualization.

import numpy as np
from sklearn.metrics import auc
import matplotlib.pyplot as plt
import os

plt.rcParams['font.family'] = 'DejaVu Sans'

current_path = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_path)
result_path = os.path.join(project_root, 'output/')

RF_name = np.load(result_path + '2.RF_name.npy')
x_axis = np.load(result_path + '1.x_axis_channel.npy')
channel_num = len(x_axis)
print(RF_name)

pos_pR_fp_scale_x = np.load(result_path + '5.pos_pR_fp_scale_x.npy', allow_pickle=True).item()

thresholds_1 = np.arange(0, 1, 0.01)
thresholds_2 = np.arange(1, 10, 0.1)
thresholds = np.concatenate((thresholds_1, thresholds_2))

count = -1
result_dict = {}
auc_dict = {}
for key, val in pos_pR_fp_scale_x.items():
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
    auc_dict[key] = roc_auc

    plt.figure(figsize=(6, 4))

    plt.plot(FPR, TPR, label=f'{key} (AUC = {roc_auc:.3f})')
    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.xlim([-0.01, 1.05])
    plt.ylim([-0.01, 1.05])
    plt.legend(loc='lower right')

    file_path = result_path + '6.unmixing_pos_ROC/'
    os.makedirs(file_path, exist_ok=True)
    filename = key + '_ROC.png'
    plt.savefig(file_path + filename)
    plt.close()

    result_dict[key] = np.array(result_list)

plt.figure(figsize=(12, 8))

for key, result_array in result_dict.items():
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

file_path2 = result_path + '6.unmixing_pos_ROC_summary/'
os.makedirs(file_path2, exist_ok=True)
filename2 = 'All_Plasmids_ROCs.png'
plt.savefig(file_path2 + filename2, transparent=True)
plt.close()

# to create the 3x6 subplot with 3 decimals of each AUC

fig, axs = plt.subplots(3, 6, figsize=(20, 11))
data_list = []
for key, result_array in result_dict.items():
    FPR = result_array[:, 0]
    TPR = result_array[:, 1]
    roc_auc = auc(FPR, TPR)
    data_list.append([key[3:], FPR, TPR, roc_auc])

for i in range(3):
    for j in range(6):
        ax = axs[i, j]
        data = data_list[i * 6 + j]

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

filename3 = 'Summary_of_Plasmids_ROC(3x6).png'
plt.savefig(file_path2 + filename3)
plt.close()

np.save(result_path + '6.result_dict_for_ROC.npy', result_dict, allow_pickle=True)  # key: [fpr, tpr]
np.save(result_path + '6.auc_dict_for_ROC.npy', auc_dict, allow_pickle=True)
