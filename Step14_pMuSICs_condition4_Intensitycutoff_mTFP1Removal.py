# As in Step 10, we have already removed mTFP1 from the pR-FP reference and created its own optimal threshold
# dictionary. Now, we will use the unmixed results and these thresholds to determine the most likely FP combination for
# each pMuSIC.
# a. Remove all pMuSICs containing mTFP1 based on the sequencing results.
# b. Perform unmixing using the 17-FP reference plus autofluorescence, and retain only the last 17 scaling factors for
# each reference.
# c. Calculate the score for each cell by dividing the unmixing scaling factors [length: 17] by the optimal thresholds
# [length: 17].
# d. Identify the most likely FP combination for each cell by finding the index of the highest score, referencing the
# control library.
# e. Assign a value of 1 to the index of the most likely combination in the control library, and 0 to others. This will
# create a score list for each cell. For the entire pMuSIC-positive population, this will result in a score array with
# dimensions (136, cell number).


import numpy as np
import pandas as pd
import re
from scipy.optimize import nnls
import matplotlib.pyplot as plt
import os
from my_module import print_lower_triangular_list, find_combination_index
import seaborn as sns

plt.rcParams['font.family'] = 'DejaVu Sans'
project_root = os.path.dirname(os.path.abspath(__file__))
result_path = os.path.join(project_root, 'output/')

RF_peak_remove_mTFP1 = np.load(result_path + '10_1.RF_peak_17.npy', allow_pickle=True).item()
RF_name_remove_mTFP1 = np.load(result_path + '10_1.RF_name_17.npy')
optimal_threshold_remove_mTFP1 = np.load(result_path + '10_2.optimal_threshold_17dict.npy', allow_pickle=True).item()
all_pos_pMuSICs = np.load(result_path + '13.all_pos_pMuSICs.npy', allow_pickle=True).item()
RF_remove_mTFP1 = np.load(result_path + '10_1.RF_17.npy')

# a. remove all pMuSICs containing mTFP1 according to the sequencing result.
# to extract the data from the summary list of sequencing result from PlasmidSaurus and set it as df_actual
excel_path = 'summary_of_sequencing_results_PlasmidSaurus/summary of sequencing results of pMuSICs.xlsx'
df_actual = pd.read_excel(excel_path, dtype={'fp_m': str, 'fp_n': str})

# to filter out any pMuSICs containing mTFP1 in the summary of the sequencing result (actual pMuSIC)
pos_pMuSIC_remove_mTFP1 =  []
for i in range(len(df_actual['Name'])):
    if '08' not in df_actual['fp_m'][i] and '08' not in df_actual['fp_n'][i]:
        pos_pMuSIC_remove_mTFP1.append(df_actual['Name'][i])

# to create the dataframe for df_actual_remove_mTFP1
df_actual_remove_mTFP1 = pd.DataFrame({
    'Name': pos_pMuSIC_remove_mTFP1
})

df_actual_remove_mTFP1['fp_m'] = df_actual_remove_mTFP1['Name'].map(df_actual.set_index('Name')['fp_m'])
df_actual_remove_mTFP1['fp_n'] = df_actual_remove_mTFP1['Name'].map(df_actual.set_index('Name')['fp_n'])
df_actual_remove_mTFP1['short_name'] = df_actual_remove_mTFP1['Name'].str.extract(r'pMuSIC_(.*)')[0]

pd.set_option('display.max_rows', None)
print(df_actual_remove_mTFP1)

file_path = result_path + '14.pMuSICs_condition4_unmixing/remove_mTFP1_on_sequencing_result/'
os.makedirs(file_path, exist_ok=True)
df_actual_remove_mTFP1.to_excel(file_path + '14.pMuSICs after removing mTFP1 (actual).xlsx', index=False)

# to remove all pMuSICs containing mTFP1 according to the 'df_actual_remove_mTFP1'
short_names_list = df_actual_remove_mTFP1['short_name'].tolist()
pos_pMuSIC_remove_mTFP1_dict = {}
for key in all_pos_pMuSICs.keys():
    match = re.match(r'^S(\d+)', key)
    if match:
        key_short_name = 'S' + match.group(1)
        print(key_short_name)
        if key_short_name in short_names_list:
            pos_pMuSIC_remove_mTFP1_dict[key] = all_pos_pMuSICs[key]


# now we have the full pMuSIC dict that doesn't contain any mTFP1
# b. unmixing using the 17 fp reference + auto_fluorescence and only take the last 17 scale_x for each reference
fp_name_remove_mTFP1 = ["01.EBFP2", "02.mTagBFP2", "03.mT-Sapphire", "04.mAmetrine", "05.mCerulean3", "06.LSSmOrange",
                        "07.mBeRFP","09.EGFP", "10.CyOFP1", "11.mClover3", "12.mVenus", "13.mPapaya", "14.mOrange2",
                        "15.mRuby3", "16.mKate2", "17.mCardinal", "18.miRFP670"]

# to match the name of each fp in RF_peak_remove_mTFP1 to the fp_name_remove_mTFP1
# thus, we could get the peak channel of each reference if using the same key in positive pMuSIC dictionary
updated_RF_peak = {}
for key, val in RF_peak_remove_mTFP1.items():
    # print(key) # output: RF_fp01, RF_fp02,...RF_fp18
    for name in fp_name_remove_mTFP1:
        if key[-2:] in name:
            updated_RF_peak[name] = int(val[0])
# output:
# {'01.EBFP2': 2, '02.mTagBFP2': 2, '03.mT-Sapphire': 4, '04.mAmetrine': 6, '05.mCerulean3': 4, '06.LSSmOrange': 7,
# '07.mBeRFP': 9, '09.EGFP': 16, '10.CyOFP1': 20, '11.mClover3': 17, '12.mVenus': 18, '13.mPapaya': 18,
# '14.mOrange2': 30, '15.mRuby3': 30, '16.mKate2': 33, '17.mCardinal': 33, '18.miRFP670': 41}

# to filter out all low FI cells
pos_pMuSIC_remove_mTFP1_FIcutoff = {}
for key, val in pos_pMuSIC_remove_mTFP1_dict.items():
    high_FI_cell = []
    for each_cell in val:
        peak_channel = np.argmax(each_cell)
        if each_cell[peak_channel] > 10**4:
            high_FI_cell.append(each_cell)
        else:
            continue

    pos_pMuSIC_remove_mTFP1_FIcutoff[key] = np.array(high_FI_cell)

# to get the optimal threshold list
threshold_list = []
for key, val in optimal_threshold_remove_mTFP1.items():
    threshold_list.append(val)
threshold = np.array(threshold_list)


# define the control library of 136 combination
fp_num = len(threshold_list)
triangle_matrix = np.zeros((fp_num, fp_num), dtype=object)

# to get 136 combo of fps
combinations = []
for i in range(fp_num):
    for j in range(i+1, fp_num):
        combination = (fp_name_remove_mTFP1[i], fp_name_remove_mTFP1[j])
        triangle_matrix[i, j] = combination
        combinations.append(combination)
len_combinations = len(combinations)
print('the number of total fp combinations with 17x17 fp is:' ,len_combinations)
print_lower_triangular_list(combinations)

# to get the control list for 136 combos
ctrl_list = []
for idx, combo in enumerate(combinations):
    ctrl_list.append([idx, combo])
    print(f'Index {idx}: {combo}')
ctrl_array = np.array([(x, y) for x, y in ctrl_list], dtype=object)


def create_triangle_heatmap(data_list, name):
    fp_name = ["EBFP2", "mTagBFP2", "mT-Sapphire", "mAmetrine", "mCerulean3", "LSSmOrange", "mBeRFP", "EGFP", "CyOFP1",
               "mClover3", "mVenus", "mPapaya", "mOrange2", "mRuby3", "mKate2", "mCardinal", "miRFP670"]
    n = 17  # Matrix size
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
                linewidths=0.5, linecolor='white', xticklabels=fp_name, yticklabels=fp_name)

    # Color the diagonal gray
    for i in range(n):
        plt.gca().add_patch(plt.Rectangle((i, i), 1, 1, fill=True, facecolor='lightgray', edgecolor='gray'))

    plt.subplots_adjust(bottom=0.18)
    file_path1 = result_path + '14.pMuSICs_condition4_unmixing/triangular_matrix_for_each_fp/'
    os.makedirs(file_path1, exist_ok=True)
    filename1 = name + '.png'
    plt.savefig(file_path1 + filename1)
    plt.close()


# start to do the unmixing
RF = np.array(RF_remove_mTFP1).T
pMuSIC_name = []
combination_index = []
potential_combination = []
combination_score = []
Score_lists = []

pMuSIC_scale_x = {}
pMuSIC_score = {}
for key, value in pos_pMuSIC_remove_mTFP1_FIcutoff.items():
    print(key)
    pMuSIC_name.append(key)
    scale_x = []
    for each_cell in value:
        x, residuals = nnls(RF, each_cell)
        x = np.round(x, decimals=4)
        scale_x.append(x)
    scale_x = np.array(scale_x)
    # print('scale_x', len(scale_x), scale_x)
    pMuSIC_scale_x[key] = scale_x

    total_cells = len(scale_x)
    print("total_cells", total_cells)

    result_list = []
    for each_cell in scale_x:
        score = each_cell[1:] / threshold

        max_val = max(score)
        max_idx = np.argmax(score)

        second_max_val = -float('inf')
        second_max_idx = -1

        for idx, val in enumerate(score):
            if val > second_max_val and val != max_val:
                second_max_val = val
                second_max_idx = idx

        most_likely_combination = (fp_name_remove_mTFP1[max_idx], fp_name_remove_mTFP1[second_max_idx])
        score_list = [0] * len(combinations)
        index = find_combination_index(most_likely_combination, ctrl_array)
        score_list[index] = 1
        result_list.append(score_list)

    result = np.array(result_list)
    column_sums = np.sum(result, axis=0)
    percentage = np.round(column_sums / total_cells, decimals=4)
    Score_lists.append(percentage)

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

df = pd.DataFrame({
    'Name': pMuSIC_name,
    'ID_136combos': combination_index,
    'Combination_136combos': potential_combination,
    'Score': combination_score,
    'Score_list_136combos': Score_lists
})
file_path2 = result_path + '14.pMuSICs_condition4_unmixing/excel_summary/'
os.makedirs(file_path2, exist_ok=True)

df.to_excel(file_path2 + '14.all_pos_pMuSIC_list_136combos.xlsx', index=False)

np.save(result_path + '14.condition4_remove_mTFP1_pMuSIC_136_control_array.npy', ctrl_array)
np.save(result_path + '14.condition4_remove_mTFP1_pMuSIC_136_combination_list.npy', combinations)
np.save(result_path + '14.condition4_remove_mTFP1_all_pos_pMuSICs_dict.npy', pos_pMuSIC_remove_mTFP1_dict, allow_pickle=True)
np.save(result_path + '14.condition4_remove_mTFP1_FI_cutoff_all_pos_pMuSICs_dict.npy', pos_pMuSIC_remove_mTFP1_FIcutoff, allow_pickle=True)
np.save(result_path + '14.condition4_remove_mTFP1_FI_cutoff_all_pos_pMuSIC_score_136combos.npy', pMuSIC_score, allow_pickle=True)
np.save(result_path + '14.condition4_remove_mTFP1_FI_cutoff_all_pos_pMuSIC_scale_x.npy', pMuSIC_scale_x, allow_pickle=True)
