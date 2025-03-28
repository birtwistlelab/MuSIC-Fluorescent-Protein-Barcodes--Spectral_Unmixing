import numpy as np
import pandas as pd
from scipy.optimize import nnls
import matplotlib.pyplot as plt
import os
from my_module import print_lower_triangular_list, find_combination_index
import seaborn as sns

plt.rcParams['font.family'] = 'DejaVu Sans'

current_path = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_path)
result_path = os.path.join(project_root, 'output/')

pos_pMuSIC = np.load(os.path.join(result_path, '10.all_pos_pMuSICs.npy'), allow_pickle=True).item()
RF = np.load(result_path + '2.MFI_RF.npy')

# To filter out all low FI cells in cells expressing pMuSICs:
pos_pMuSIC_filtered = {}
for key, val in pos_pMuSIC.items():
    high_FI_cell = []
    for each_cell in val:
        peak_channel = np.argmax(each_cell)
        if each_cell[peak_channel] >= 14000:
            high_FI_cell.append(each_cell)
        else:
            continue
    pos_pMuSIC_filtered[key] = np.array(high_FI_cell)

# Get the optimal thresholds of all FPs
optimal_thresholds = np.load(os.path.join(result_path, '6.DetermineUnmixingThresholdsOnTrainingCells_FigS4/average_optimal_threshold.npy'))
thresholds = np.array(optimal_thresholds)

# define the control library of 153 combination
fp_name = ['EBFP2', 'mTagBFP2', 'mT_Sapphire', 'mAmetrine', 'mCerulean3', 'LSSmOrange', 'mBeRFP', 'mTFP1', 'EGFP',
           'CyOFP1', 'mClover3', 'mVenus', 'mPapaya', 'mOrange2', 'mRuby3', 'mKate2', 'mCardinal', 'miRFP670']

fp_num = len(fp_name)
triangle_matrix = np.zeros((fp_num, fp_num), dtype=object)

combinations = []
for i in range(fp_num):
    for j in range(i+1, fp_num):
        combination = (fp_name[i], fp_name[j])
        triangle_matrix[i, j] = combination
        combinations.append(combination)
len_combinations = len(combinations)
print('the number of total fp combinations with 18x18 fp is:', len_combinations)
print_lower_triangular_list(combinations)

# Get the control list for 153 combos
ctrl_list = []
for idx, combo in enumerate(combinations):
    ctrl_list.append([idx, combo])
    print(f'Index {idx}: {combo}')
ctrl_array = np.array([(x, y) for x, y in ctrl_list], dtype=object)


def create_triangle_heatmap(data_list, name):
    n = len(fp_name)
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
    file_path1 = result_path + '11.pMuSICsUnmixing_IntensityCutoff/triangular_matrix_for_each_fp/'
    os.makedirs(file_path1, exist_ok=True)
    filename1 = name + '.png'
    plt.savefig(file_path1 + filename1)
    plt.close()


# Start to do the unmixing
RF = np.array(RF).T
pMuSIC_name = []
combination_index = []
potential_combination = []
combination_score = []
Score_lists = []

pMuSIC_scale_x = {}
pMuSIC_score = {}
for key, value in pos_pMuSIC_filtered.items():
    print(key)
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
        score = each_cell[1:] / optimal_thresholds

        max_val = max(score)
        max_idx = np.argmax(score)

        second_max_val = -float('inf')
        second_max_idx = -1

        for idx, val in enumerate(score):
            if val > second_max_val and val != max_val:
                second_max_val = val
                second_max_idx = idx

        most_likely_combination = (fp_name[max_idx], fp_name[second_max_idx])
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

df2 = pd.DataFrame({
    'Name': pMuSIC_name,
    'ID': combination_index,  # inferred ID from the unmixing results
    'Combination': potential_combination,  # inferred combination from the unmixing results
    'Score': combination_score,  # highest scores from the unmixing results
    'Score_list': Score_lists  # score list (1x153)
})
file_path2 = result_path + '11.pMuSICsUnmixing_IntensityCutoff/excel_summary/'
os.makedirs(file_path2, exist_ok=True)

df2.to_excel(file_path2 + '11.all_pos_pMuSIC_list.xlsx', index=False)

np.save(os.path.join(result_path, '11.pMuSIC_control_array.npy'), ctrl_array)
np.save(os.path.join(result_path, '11.pMuSIC_combination_list.npy'), combinations)  # combinations from ctrl_array
np.save(os.path.join(result_path, '11.all_pos_pMuSICs_dict_Cutoff.npy'), pos_pMuSIC_filtered, allow_pickle=True)
np.save(os.path.join(result_path, '11.all_pos_pMuSIC_Cutoff_score.npy'), pMuSIC_score, allow_pickle=True)
np.save(os.path.join(result_path, '11.all_pos_pMuSIC_Cutoff_scale_x.npy'), pMuSIC_scale_x, allow_pickle=True)
