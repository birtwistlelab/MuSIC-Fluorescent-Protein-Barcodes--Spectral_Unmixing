# to find the optimal threshold for each pR_fp by identifying the threshold point on the ROC curve that balances higher
# TPR and lower FPR which is the one closest to the ideal point(0,1) as (FPR, TPR).

import numpy as np
from sklearn.metrics import auc
import matplotlib.pyplot as plt
import os
plt.rcParams['font.family'] = 'DejaVu Sans'
project_root = os.path.dirname(os.path.abspath(__file__))
result_path = os.path.join(project_root, 'output/')

x_axis = np.load(result_path + '1.x_axis_channel.npy')
channel_num = len(x_axis)
pos_pR_fp_scale_x = np.load(result_path + '5.pos_pR_fp_scale_x.npy', allow_pickle=True).item()

thresholds_1 = np.arange(0, 1, 0.01)
thresholds_2 = np.arange(1, 10, 0.1)
thresholds = np.concatenate((thresholds_1, thresholds_2), axis=0)

result_dict = np.load(result_path + '6.result_dict_for_ROC.npy', allow_pickle=True).item()

optimal_threshold = {}
for key, val in result_dict.items():
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

    file_path = result_path + '7.find_the_best_threshold/'
    os.makedirs(file_path, exist_ok=True)
    filename5 = key + '_ROC.png'
    plt.savefig(file_path + filename5)
    plt.close()

    optimal_threshold[key] = round(thresholds[optimal_idx], 2)
    print(f'Optimal Point: FPR = {round(optimal_FPR, 2)}, TPR = {round(optimal_TPR, 2)}')

# to create the Fig.3B for the paper, 6x3 subplot with 3 decimals of each AUC,
# add the best threshold

fig, axs = plt.subplots(6, 3, figsize=(9.6, 16))

data_list = []
for key, result_array in result_dict.items():
    FPR = result_array[:, 0]
    TPR = result_array[:, 1]
    roc_auc = auc(FPR, TPR)

    # ideal_point would be (0,1)
    distances = np.sqrt((FPR - 0) ** 2 + (TPR - 1) ** 2)

    optimal_idx = np.argmin(distances)
    optimal_FPR = FPR[optimal_idx]
    optimal_TPR = TPR[optimal_idx]

    data_list.append([key[3:], FPR, TPR, roc_auc, optimal_FPR, optimal_TPR])

for i in range(6):
    for j in range(3):
        ax = axs[i, j]
        data = data_list[i * 3 + j]

        ax.plot(data[1], data[2], color='blue')
        ax.plot([], [], ' ', label=f'{data[0]}\nAUC = {data[3]:.3f}')
        ax.plot([0, 1], [0, 1], color='gray', linestyle='--')

        ax.scatter(data[4], data[5], color='red', zorder=2)
        offset = 0.02
        ax.text(data[4] + offset, data[5], f'({data[4]:.3f}, {data[5]:.3f})',
                fontsize=16, va='top', ha='left')

        ax.set_xlim([-0.01, 1.05])
        ax.set_ylim([-0.01, 1.05])
        ax.tick_params(axis='x', labelsize=14)
        ax.tick_params(axis='y', labelsize=14)
        ax.legend(frameon=False, loc='lower right', fontsize=16, bbox_to_anchor=(1.05, -0.05))
        ax.set_aspect('equal', adjustable='box')

        if i == 5 and j == 1:
            ax.set_xlabel('False Positive Rate', fontsize=20, labelpad=12)

        y_label_horizontal_pos = 0.015
        y_label_vertical_pos = 0.5
        fig.text(y_label_horizontal_pos, y_label_vertical_pos, 'True Positive Rate',
                 rotation=90, ha='center', fontsize=20)

plt.subplots_adjust(left=0.2, right=0.8, hspace=0.02, wspace=0.005, bottom=0.4)
plt.tight_layout()

file_path2 = result_path + '7.find_the_best_threshold(3x6)/'
os.makedirs(file_path2, exist_ok=True)
filename2 = 'Summary_of_Plasmids_ROC_with_Optimal_Thresholds.png'
plt.savefig(file_path2 + filename2, dpi=300, transparent=True)

fig_path = os.path.join(project_root, 'paper_fig/')
filename2 = 'Fig3B(Step7).png'
plt.savefig(fig_path + filename2, dpi=300, transparent=True)
plt.close()

np.save(result_path + '7.optimal_threshold_dict.npy', optimal_threshold, allow_pickle=True)