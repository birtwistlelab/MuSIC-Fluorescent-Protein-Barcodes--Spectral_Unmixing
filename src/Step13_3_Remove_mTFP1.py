# ================================================ XL Clemson ================================================
# After the optimal thresholds and the cells were ready, we could start to unmix those pMuSICs without mTFP1
# ================================================ XL Clemson ================================================

import numpy as np
import pandas as pd
from scipy.optimize import nnls
import matplotlib.pyplot as plt
import os
import re
from my_module import print_lower_triangular_list, find_combination_index
import seaborn as sns

plt.rcParams['font.family'] = 'Arial'

current_path = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_path)
result_path = os.path.join(project_root, 'output/')

saving_path = os.path.join(result_path, '13_3.Remove_mTFP1')
os.makedirs(saving_path, exist_ok=True)

pos_pMuSIC = np.load(os.path.join(result_path, '10.all_pos_pMuSICs.npy'), allow_pickle=True).item()
RF_original = np.load(result_path + '2.MFI_RF.npy')

# Remove the reference of mTFP1 (index 8) from the original RF data
RF_17 = np.delete(RF_original, 8, axis=0)

pos_pMuSIC_filtered = np.load(os.path.join(result_path, '11.all_pos_pMuSICs_dict_Cutoff.npy'), allow_pickle=True).item()

# Extract and filter pMuSICs based on sequencing results (from PlasmidSaurus)
excel_path = os.path.join(project_root, 'csv_files', 'Sequencing_Results_PlasmidSaurus',
                          'summary of sequencing results of pMuSICs.xlsx')
df_actual = pd.read_excel(excel_path, dtype={'fp_m': str, 'fp_n': str})

# Filter out pMuSICs containing mTFP1 based on the sequencing data
pos_pMuSIC_remove_mTFP1 = []
for i in range(len(df_actual['Name'])):
    if '08' not in df_actual['fp_m'][i] and '08' not in df_actual['fp_n'][i]:
        pos_pMuSIC_remove_mTFP1.append(df_actual['Name'][i])


# Create a DataFrame with pMuSICs that do not contain mTFP1
df_actual_remove_mTFP1 = pd.DataFrame({
    'Name': pos_pMuSIC_remove_mTFP1
})

df_actual_remove_mTFP1['fp_m'] = df_actual_remove_mTFP1['Name'].map(df_actual.set_index('Name')['fp_m'])
df_actual_remove_mTFP1['fp_n'] = df_actual_remove_mTFP1['Name'].map(df_actual.set_index('Name')['fp_n'])
df_actual_remove_mTFP1['short_name'] = df_actual_remove_mTFP1['Name'].str.extract(r'pMuSIC_(.*)')[0]

pd.set_option('display.max_rows', None)
print(df_actual_remove_mTFP1)

df_actual_remove_mTFP1.to_excel(os.path.join(saving_path, '13_3.pMuSICs after removing mTFP1 (actual).xlsx'), index=False)

# Remove pMuSICs that contain mTFP1 from the filtered pMuSICs dictionary
short_names_list = df_actual_remove_mTFP1['short_name'].tolist()
pos_pMuSIC_remove_mTFP1_dict = {}
for key in pos_pMuSIC_filtered.keys():
    match = re.match(r'^S(\d+)', key)
    if match:
        key_short_name = 'S' + match.group(1)
        print(key_short_name)
        if key_short_name in short_names_list:
            pos_pMuSIC_remove_mTFP1_dict[key] = pos_pMuSIC_filtered[key]


# Get the optimal thresholds of all 17 FPs
optimal_thresholds_17 = np.load(os.path.join(result_path, '13_2.DetermineUnmixingThresholdsOnTrainingCells_Remove_mTFP1/average_optimal_threshold_Remove_mTFP1.npy'))

# Define the control library for 136 combinations of FPs
fp_name_17 = ['EBFP2', 'mTagBFP2', 'mT_Sapphire', 'mAmetrine', 'mCerulean3', 'LSSmOrange', 'mBeRFP', 'EGFP',
           'CyOFP1', 'mClover3', 'mVenus', 'mPapaya', 'mOrange2', 'mRuby3', 'mKate2', 'mCardinal', 'miRFP670']

fp_num = len(fp_name_17)
triangle_matrix = np.zeros((fp_num, fp_num), dtype=object)

# Create all combinations of 17 FPs
combinations = []
for i in range(fp_num):
    for j in range(i+1, fp_num):
        combination = (fp_name_17[i], fp_name_17[j])
        triangle_matrix[i, j] = combination
        combinations.append(combination)
len_combinations = len(combinations)
print('the number of total fp combinations with 17x17 fp is:', len_combinations)
print_lower_triangular_list(combinations)
np.save(os.path.join(saving_path, '13_3.pMuSIC_combination_list_remove_mTFP1.npy'), combinations)

# Get the control list for 136 combos
ctrl_list = []
for idx, combo in enumerate(combinations):
    ctrl_list.append([idx, combo])
    print(f'Index {idx}: {combo}')
ctrl_array = np.array([(x, y) for x, y in ctrl_list], dtype=object)
np.save(os.path.join(saving_path, '13_3.pMuSIC_control_array_remove_mTFP1.npy'), ctrl_array)


# Function to create a heatmap for the combination scores
def create_triangle_heatmap(data_list, name):
    n = len(fp_name_17)
    triangle_matrix = np.full((n, n), np.nan)  # Initialize the matrix with NaNs

    # Fill the lower triangle of the matrix (excluding the diagonal)
    index = 0
    for j in range(n):
        for i in range(j + 1, n):
            if index < len(data_list):
                triangle_matrix[i, j] = data_list[index]
                index += 1

    # Create a DataFrame
    df = pd.DataFrame(triangle_matrix)

    # Create a mask for the upper triangle and the diagonal
    mask = np.triu(np.ones_like(triangle_matrix, dtype=bool))
    np.fill_diagonal(mask, True)

    # Plot the heatmap
    plt.figure(figsize=(10, 8))
    sns.heatmap(df, annot=True, fmt=".1f", cmap="Blues", mask=mask, cbar=True, vmax=1,
                linewidths=0.5, linecolor='white', xticklabels=fp_name_17, yticklabels=fp_name_17)

    # Color the diagonal gray
    for i in range(n):
        plt.gca().add_patch(plt.Rectangle((i, i), 1, 1, fill=True, facecolor='lightgray', edgecolor='gray'))

    plt.subplots_adjust(bottom=0.18)
    file_path1 = result_path + '13_3.Remove_mTFP1/triangular_matrix_for_each_fp/'
    os.makedirs(file_path1, exist_ok=True)
    filename1 = name + '.png'
    plt.savefig(file_path1 + filename1)
    plt.close()

# Start the unmixing process
RF = np.array(RF_17).T
pMuSIC_name = []
combination_index = []
potential_combination = []
combination_score = []
Score_lists = []

pMuSIC_scale_x = {}
pMuSIC_score = {}
for key, value in pos_pMuSIC_remove_mTFP1_dict.items():
    pMuSIC_name.append(key)
    scale_x = []
    for each_cell in value:
        x, residuals = nnls(RF, each_cell)
        x = np.round(x, decimals=4)
        scale_x.append(x)
    scale_x = np.array(scale_x)
    pMuSIC_scale_x[key] = scale_x

    total_cells = len(scale_x)

    result_list = []
    for each_cell in scale_x:
        score = each_cell[1:] / optimal_thresholds_17  # Normalize scores based on optimal thresholds

        max_val = max(score)
        max_idx = np.argmax(score)

        second_max_val = -float('inf')
        second_max_idx = -1

        for idx, val in enumerate(score):
            if val > second_max_val and val != max_val:
                second_max_val = val
                second_max_idx = idx

        # Determine the most likely combination based on the highest scores
        most_likely_combination = (fp_name_17[max_idx], fp_name_17[second_max_idx])
        score_list = [0] * len(combinations)
        index = find_combination_index(most_likely_combination, ctrl_array)
        score_list[index] = 1
        result_list.append(score_list)

    result = np.array(result_list)
    column_sums = np.sum(result, axis=0)
    percentage = np.round(column_sums / total_cells, decimals=4)
    Score_lists.append(percentage)

    # Find the combination with the highest score
    possible_index = np.argmax(percentage)
    positive_rate = np.max(percentage)

    print("the index of the most likely pMuSIC is " + str(possible_index))
    combination_index.append(possible_index)

    print("the positive rate is " + str(positive_rate))
    combination_score.append(positive_rate)

    # ctrl_list =[idx, combo]
    pMuSIC_combination = ctrl_list[possible_index][1]
    print("the most likely combination of pMuSIC is ", pMuSIC_combination)
    potential_combination.append(pMuSIC_combination)

    percentage_list = percentage.tolist()
    create_triangle_heatmap(percentage_list, key)

    print_lower_triangular_list(percentage_list)
    pMuSIC_score[key] = percentage

# Create a DataFrame with the results
df2 = pd.DataFrame({
    'Name': pMuSIC_name,
    'ID': combination_index,  # inferred ID from the unmixing results
    'Combination': potential_combination,  # inferred combination from the unmixing results
    'Score': combination_score,  # highest scores from the unmixing results
    'Score_list': Score_lists  # score list (1x136)
})
file_path2 = os.path.join(result_path, '13_3.Remove_mTFP1/excel_summary/')
os.makedirs(file_path2, exist_ok=True)

df2.to_excel(os.path.join(file_path2, '13_3.all_pos_pMuSIC_list(mTFP1_removed).xlsx'))

np.save(os.path.join(saving_path, '13_3.pos_pMuSICs_remove_mTFP1_dict.npy'), pos_pMuSIC_remove_mTFP1_dict, allow_pickle=True)
np.save(os.path.join(saving_path, '13_3.pos_pMuSIC_score_remove_mTFP1.npy'), pMuSIC_score, allow_pickle=True)
np.save(os.path.join(saving_path, '13_3.pos_pMuSIC_scale_x_remove_mTFP1.npy'), pMuSIC_scale_x, allow_pickle=True)
