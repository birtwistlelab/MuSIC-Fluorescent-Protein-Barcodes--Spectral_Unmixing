# ===================================================== XL Clemson =====================================================
# From the 18x18 plots from Step5, we systematically increase the threshold of scale_x. For each histogram along each
# column, we obtain 1 True Positive (TP), 1 False Negative (FN), 17 True Negatives (TNs), and 17 False Positives (FPs).
# These results allow us to generate 18 individual ROC curves—one for each pR-fp, as well as a combined ROC curve that
# consolidates all 18 for better visualization. By analyzing the ROC curves, we can determine the optimal threshold that
# best balances the True Positive Rate (TPR) and False Positive Rate (FPR). Thus, for each fp, we will have one optimal
# thresholds of unmixed scale_x. Since we have 3 tests, the final optimal thresholds would be the average of all tests.
# ===================================================== XL Clemson =====================================================

import numpy as np
from sklearn.metrics import auc
import matplotlib.pyplot as plt
import os
import pandas as pd
plt.rcParams['font.family'] = 'Arial'

current_path = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_path)
result_path = os.path.join(project_root, 'output/')

loading_paths = [os.path.join(result_path, '5.UnmixingTransfectedTrainingCells/', f'test{i + 1}/') for i in range(3)]

saving_paths = [os.path.join(result_path, '6.DetermineUnmixingThresholdsOnTrainingCells_FigS4/', f'test{i + 1}/') for i in range(3)]
for path in saving_paths:
    os.makedirs(path, exist_ok=True)

saving_paths2 = os.path.join(result_path, '6.DetermineUnmixingThresholdsOnTrainingCells_FigS4', 'average_data')
os.makedirs(saving_paths2, exist_ok=True)

saving_paths3 = os.path.join(project_root, 'figures_for_paper/')
os.makedirs(saving_paths3, exist_ok=True)

# Load required data files：
RF_name = np.load(result_path + '2.RF_name.npy')
x_axis = np.load(result_path + '1.x_axis_channel.npy')
channel_num = len(x_axis)

fp_name = ['01.EBFP2', '02.mTagBFP2', '03.mT_Sapphire', '04.mAmetrine', '05.mCerulean3', '06.LSSmOrange', '07.mBeRFP',
           '08.mTFP1', '09.EGFP', '10.CyOFP1', '11.mClover3', '12.mVenus', '13.mPapaya', '14.mOrange2', '15.mRuby3',
           '16.mKate2', '17.mCardinal', '18.miRFP670']

# Define threshold values for ROC computation:
thresholds_1 = np.arange(0, 1, 0.01)  # Fine-grained thresholds from 0 to 1
thresholds_2 = np.arange(1, 10, 0.1)  # Coarser thresholds from 1 to 10
thresholds = np.concatenate((thresholds_1, thresholds_2))

scale_x_threshold_test_list = []  # Store optimal thresholds for each test
result_dicts = []  # Store ROC results for each test

# Loop through three test sets：
for i in range(3):
    transfected_pR_fp_scale_x = np.load(os.path.join(loading_paths[i], '5.transfected_Training_pR_fp_scale_x.npy'), allow_pickle=True).item()
    count = -1
    result_dict = {}  # Store ROC results for each FP
    auc_dict = {}  # Store AUC values for each FP
    optimal_threshold = {}  # Store optimal thresholds for each FP
    data_list = []  # Store FP data for combined visualization

    # Iterate through each reference factor
    for key, val in transfected_pR_fp_scale_x.items():
        count += 1
        # Skip RF0 column and process other RF channels
        new_val = val[:, 1:]
        result_list = []
        TPR = []
        FPR = []

        # Iterate through all thresholds to compute TPR and FPR
        for each_threshold in thresholds:
            true_posCell_list = []
            false_negCell_list = []
            true_negCell_list = []
            false_posCell_list = []

            # Loop through all columns (FPs) to classify cells
            for col in range(new_val.shape[1]):
                RF_array = new_val[:, col]
                for idx in range(len(RF_array)):
                    if col == count:  # True Positives and False Negatives
                        if RF_array[idx] < each_threshold:
                            false_negCell_list.append(idx)
                        else:
                            true_posCell_list.append(idx)
                    else:  # True Negatives and False Positives
                        if RF_array[idx] < each_threshold:
                            true_negCell_list.append(idx)
                        else:
                            false_posCell_list.append(idx)

            # Calculate TPR and FPR
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

        # ideal_point would be (0,1)
        distances = np.sqrt((FPR - 0) ** 2 + (TPR - 1) ** 2)

        # Find the index of the closest point to (0,1)
        optimal_idx = np.argmin(distances)
        optimal_FPR = FPR[optimal_idx]
        optimal_TPR = TPR[optimal_idx]
        data_list.append([key[3:], FPR, TPR, roc_auc, optimal_FPR, optimal_TPR])

        print('the optimal_idx of pR_' + key + ' is ' + str(optimal_idx))
        print('the best threshold is ' + str(round(thresholds[optimal_idx], 2)))

        # Plot the ROC curve with the optimal point
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

        # Save individual ROC plot
        file_path = os.path.join(saving_paths[i], '6.unmixing_pos_ROC/')
        os.makedirs(file_path, exist_ok=True)
        filename = key + '_ROC.png'
        plt.savefig(os.path.join(file_path, filename))
        plt.close()

        # Store ROC and optimal threshold data
        result_dict[key] = np.array(result_list)
        optimal_threshold[key] = round(thresholds[optimal_idx], 4)
        print(f'Optimal Point: FPR = {round(optimal_FPR, 3)}, TPR = {round(optimal_TPR, 3)}')

    result_dicts.append(result_dict)

    scale_x_threshold_test_list.append(optimal_threshold)
    np.save(os.path.join(saving_paths[i], '6.result_dict_for_ROC.npy'), result_dict, allow_pickle=True)  # key: [fpr, tpr]
    np.save(os.path.join(saving_paths[i], '6.auc_dict_for_ROC.npy'), auc_dict, allow_pickle=True)
    np.save(os.path.join(saving_paths[i], '6.optimal_threshold_dict.npy'), optimal_threshold, allow_pickle=True)

# =================== XL Clemson ===================
# Averaging Optimal Thresholds Across Tests
# =================== XL Clemson ===================

average_dict = {}
ave_data_list = []

# Compute average ROC and thresholds across 3 tests
for key in result_dicts[0].keys():
    optimal_idx = None
    arrays = [result[key] for result in result_dicts]
    stacked_arrays = np.stack(arrays, axis=0)

    ave_FPR = np.mean(stacked_arrays[:, :, 0], axis=0)
    ave_TPR = np.mean(stacked_arrays[:, :, 1], axis=0)

    average_dict[key] = np.column_stack((ave_FPR, ave_TPR))

    # Compute average ROC AUC
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

    # Save averaged ROC plot
    filename2 = key + '_ROC.png'
    file_path2 = os.path.join(saving_paths2, '6.unmixing_average_18_ROCs/')
    os.makedirs(file_path2, exist_ok=True)
    plt.savefig(os.path.join(file_path2, filename2))
    plt.close()

    ave_data_list.append([key[3:], FPR_list, TPR_list, roc_auc, example_fpr, example_tpr])


# =================== XL Clemson ===================
# Generate Combined Figure of All ROC Curves (Fig.S9)
# =================== XL Clemson ===================

fig, axs = plt.subplots(6, 3, figsize=(11, 16))
for m in range(6):
    for n in range(3):
        optimal_idx = None
        ax = axs[m, n]
        data = ave_data_list[m * 3 + n]

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

filename3 = 'Fig.S4(Step6)Summary_of_Plasmids_ROC_with_Optimal_Thresholds.png'
plt.savefig(os.path.join(saving_paths2, filename3), dpi=300, transparent=True)
plt.savefig(os.path.join(saving_paths3, filename3), dpi=300, transparent=True)
plt.close()

# =================== XL Clemson ===================
# Save Final Average Thresholds to Excel and .npy
# =================== XL Clemson ===================

# Calculate the average thresholds
df = pd.DataFrame({
    'FP_reference': fp_name,
    'Test1_Thresholds': list(scale_x_threshold_test_list[0].values()),
    'Test2_Thresholds': list(scale_x_threshold_test_list[1].values()),
    'Test3_Thresholds': list(scale_x_threshold_test_list[2].values())
})

df['Average_thresholds'] = round(df[['Test1_Thresholds', 'Test2_Thresholds', 'Test3_Thresholds']].mean(axis=1), 3)
df.to_excel(os.path.join(result_path, '6.DetermineUnmixingThresholdsOnTrainingCells_FigS4/', 'optimal_threshold_scale_x.xlsx'),
            index=False)

average_thresholds_dict = dict(zip(df['FP_reference'], df['Average_thresholds']))
np.save(os.path.join(result_path, '6.DetermineUnmixingThresholdsOnTrainingCells_FigS4/', 'average_optimal_threshold.npy'), df['Average_thresholds'])
np.save(os.path.join(result_path, '6.DetermineUnmixingThresholdsOnTrainingCells_FigS4/', 'average_optimal_thresholds_dict.npy'), average_thresholds_dict, allow_pickle=True)
