# ================================================== XL Clemson ==================================================
# After removing mTFP1, find the optimal threshold for each pR_fp by identifying the threshold point on the ROC curve
# that balances higher TPR and lower FPR which is the one closest to the ideal point(0,1) as (FPR, TPR).
# ================================================== XL Clemson ==================================================

import numpy as np
from sklearn.metrics import auc
import matplotlib.pyplot as plt
import os
import pandas as pd
plt.rcParams['font.family'] = 'Arial'

current_path = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_path)
result_path = os.path.join(project_root, 'output/')

RF_name = np.load(result_path + '2.RF_name.npy')
x_axis = np.load(result_path + '1.x_axis_channel.npy')
channel_num = len(x_axis)

loading_paths = [os.path.join(result_path, '13_1.UnmixingTransfectedTrainingCells_Remove_mTFP1', f'test{i + 1}/') for i in range(3)]

saving_paths = [os.path.join(result_path, '13_2.DetermineUnmixingThresholdsOnTrainingCells_Remove_mTFP1', f'test{i + 1}/') for i in range(3)]
for path in saving_paths:
    os.makedirs(path, exist_ok=True)

saving_paths2 = os.path.join(result_path, '13_2.DetermineUnmixingThresholdsOnTrainingCells_Remove_mTFP1', 'average_data')
os.makedirs(saving_paths2, exist_ok=True)


fp_name_17 = ['01.EBFP2', '02.mTagBFP2', '03.mT_Sapphire', '04.mAmetrine', '05.mCerulean3', '06.LSSmOrange', '07.mBeRFP',
              '09.EGFP', '10.CyOFP1', '11.mClover3', '12.mVenus', '13.mPapaya', '14.mOrange2', '15.mRuby3', '16.mKate2',
              '17.mCardinal', '18.miRFP670']

thresholds_1 = np.arange(0, 1, 0.01)
thresholds_2 = np.arange(1, 10, 0.1)
thresholds = np.concatenate((thresholds_1, thresholds_2))

scale_x_threshold_test_list = []
result_dicts = []

for i in range(3):
    transfected_pR_fp_scale_x_17 = np.load(os.path.join(loading_paths[i],
                            '13_1.transfected_Training_pR_fp_scale_x_Remove_mTFP1.npy'), allow_pickle=True).item()
    count = -1
    result_dict = {}
    auc_dict = {}

    optimal_threshold = {}
    data_list = []

    # Iterate through each reference factor
    for key, val in transfected_pR_fp_scale_x_17.items():
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

        FPR = np.array(FPR)
        TPR = np.array(TPR)

        # Compute AUC (Area Under Curve):
        roc_auc = auc(FPR, TPR)
        auc_dict[key] = roc_auc

        # Plot individual ROC curve:
        plt.figure(figsize=(6, 4))

        # The ideal_point would be (0,1)
        distances = np.sqrt((FPR - 0) ** 2 + (TPR - 1) ** 2)

        # Find the index of the closest point to (0,1)
        optimal_idx = np.argmin(distances)
        optimal_FPR = FPR[optimal_idx]
        optimal_TPR = TPR[optimal_idx]
        data_list.append([key[3:], FPR, TPR, roc_auc, optimal_FPR, optimal_TPR])

        print('the optimal_idx of pR_' + key + ' is ' + str(optimal_idx))
        print('the best threshold is ' + str(round(thresholds[optimal_idx], 2)))

        plt.plot(FPR, TPR, label=f'{key} (AUC = {roc_auc:.3f})')
        plt.scatter(optimal_FPR, optimal_TPR, color='red', label='Optimal Point', zorder=2)

        offset = 0.02
        plt.text(optimal_FPR + offset, optimal_TPR, f'({optimal_FPR:.3f}, {optimal_TPR:.3f})',
                 fontsize=9, verticalalignment='top', horizontalalignment='left')

        plt.plot([0, 1], [0, 1], 'k--')
        plt.xlabel('False Positive Rate of Training Cells')
        plt.ylabel('True Positive Rate of Training Cells')
        plt.xlim([-0.01, 1.05])
        plt.ylim([-0.01, 1.05])
        plt.legend(loc='lower right')

        # Save ROC plots:
        file_path = os.path.join(saving_paths[i], '13_2.unmixing_pos_ROC/')
        os.makedirs(file_path, exist_ok=True)
        filename = key + '_ROC.png'
        plt.savefig(os.path.join(file_path, filename))
        plt.close()

        result_dict[key] = np.array(result_list)

        # Store optimal threshold for the current plasmid
        optimal_threshold[key] = round(thresholds[optimal_idx], 4)
        print(f'Optimal Point: FPR = {round(optimal_FPR, 3)}, TPR = {round(optimal_TPR, 3)}')
    result_dicts.append(result_dict)

    scale_x_threshold_test_list.append(optimal_threshold)
    np.save(os.path.join(saving_paths[i], '13_2.result_dict_for_ROC.npy'), result_dict,
            allow_pickle=True)  # key: [fpr, tpr]
    np.save(os.path.join(saving_paths[i], '13_2.auc_dict_for_ROC.npy'), auc_dict, allow_pickle=True)
    np.save(os.path.join(saving_paths[i], '13_2.optimal_threshold_dict.npy'), optimal_threshold, allow_pickle=True)

average_dict = {}
ave_data_list = []
for key in result_dicts[0].keys():
    optimal_idx = None
    arrays = [result[key] for result in result_dicts]
    stacked_arrays = np.stack(arrays, axis=0)

    ave_FPR = np.mean(stacked_arrays[:, :, 0], axis=0)
    ave_TPR = np.mean(stacked_arrays[:, :, 1], axis=0)

    average_dict[key] = np.column_stack((ave_FPR, ave_TPR))

    FPR_list = average_dict[key][:, 0]
    TPR_list = average_dict[key][:, 1]
    roc_auc = auc(FPR_list, TPR_list)

    plt.plot(FPR_list, TPR_list, label=f'{key} (AUC = {roc_auc:.3f})')

    distances = np.sqrt((FPR_list - 0) ** 2 + (TPR_list - 1) ** 2)

    # Find the index of the closest point to (0,1) to show how to get the optimal points of this ROC
    optimal_idx = np.argmin(distances)
    example_fpr = FPR_list[optimal_idx]
    example_tpr = TPR_list[optimal_idx]

    plt.xlabel('False Positive Rate of Training Cells')
    plt.ylabel('True Positive Rate of Training Cells')
    plt.legend(loc='lower right')

    filename2 = key + '_ROC.png'
    file_path2 = os.path.join(saving_paths2, '13_2.unmixing_average_18_ROCs/')
    os.makedirs(file_path2, exist_ok=True)
    plt.savefig(os.path.join(file_path2, filename2))
    plt.close()

    ave_data_list.append([key[3:], FPR_list, TPR_list, roc_auc, example_fpr, example_tpr])


# Create summary figure with 6x3 subplots showing all ROC curves on training cells for the paper Fig.S9
fig, axs = plt.subplots(6, 3, figsize=(11, 16))
for m in range(6):
    for n in range(3):
        optimal_idx = None
        ax = axs[m, n]
        idx = m * 3 + n
        if idx >= len(ave_data_list):
            ax.axis('off')
            continue

        data = ave_data_list[idx]

        ax.plot(data[1], data[2], color='navy', alpha=0.75)
        ax.plot([], [], ' ', label=f'{data[0]}\nAUC = {data[3]:.3f}')
        ax.plot([0, 1], [0, 1], color='gray', linestyle='--')

        ax.scatter(data[4], data[5], color='orange', zorder=2)
        ax.scatter(0, 1, color='black', s=10, zorder=2)

        ax.set_xlim([-0.01, 1.05])
        ax.set_ylim([-0.01, 1.05])
        ax.tick_params(axis='x', labelsize=14)
        ax.tick_params(axis='y', labelsize=14)
        ax.set_xticks([n * 0.25 for n in range(5)])
        ax.legend(frameon=False, loc='lower right', fontsize=12, bbox_to_anchor=(1.05, -0.05))
        ax.set_aspect('equal', adjustable='box')

        if m == 5 and n == 1:
            ax.set_xlabel('False Positive Rate of Training Cells', fontsize=20, labelpad=12)

        y_label_horizontal_pos = 0.015
        y_label_vertical_pos = 0.4
        fig.text(y_label_horizontal_pos, y_label_vertical_pos, 'True Positive Rate of Training Cells',
                 rotation=90, ha='center', fontsize=20)
plt.subplots_adjust(left=0.28, right=0.85, hspace=0.3, wspace=0.3, bottom=0.05)
plt.tight_layout()

filename3 = '(Step13_2)Summary_of_Plasmids_ROC_with_Optimal_Thresholds_Remove_mTFP1.png'
plt.savefig(os.path.join(saving_paths2, filename3), dpi=300, transparent=True)
plt.close()

# Calculate the average thresholds
df = pd.DataFrame({
    'FP_reference': fp_name_17,
    'Test1_Thresholds': list(scale_x_threshold_test_list[0].values()),
    'Test2_Thresholds': list(scale_x_threshold_test_list[1].values()),
    'Test3_Thresholds': list(scale_x_threshold_test_list[2].values())
})

df['Average_thresholds'] = round(df[['Test1_Thresholds', 'Test2_Thresholds', 'Test3_Thresholds']].mean(axis=1), 3)
df.to_excel(
    os.path.join(result_path, '13_2.DetermineUnmixingThresholdsOnTrainingCells_Remove_mTFP1/', 'optimal_threshold_scale_x_Remove_mTFP1'
                                                                                               '.xlsx'),
    index=False)

average_thresholds_dict = dict(zip(df['FP_reference'], df['Average_thresholds']))
np.save(os.path.join(result_path, '13_2.DetermineUnmixingThresholdsOnTrainingCells_Remove_mTFP1/', 'average_optimal_threshold_Remove_mTFP1.npy'),
        df['Average_thresholds'])
np.save(os.path.join(result_path, '13_2.DetermineUnmixingThresholdsOnTrainingCells_Remove_mTFP1/', 'average_optimal_thresholds_dict_Remove_mTFP1.npy'),
    average_thresholds_dict, allow_pickle=True)
