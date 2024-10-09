import numpy as np
import matplotlib.pyplot as plt
import os
import re
import statistics
import pandas as pd
from my_module import print_lower_triangular_list, find_combination_index

plt.rcParams['font.family'] = 'DejaVu Sans'
project_root = os.path.dirname(os.path.abspath(__file__))
result_path = os.path.join(project_root, 'output/')

condition4_pos_pMuSIC = np.load(result_path + '14.condition4_remove_mTFP1_FI_cutoff_all_pos_pMuSICs_dict.npy',
                                           allow_pickle=True).item()

condition4_all_pos_pMuSIC_scale_x = np.load(result_path +
                                            '14.condition4_remove_mTFP1_FI_cutoff_all_pos_pMuSIC_scale_x.npy',
                                            allow_pickle=True).item()

condition4_optimal_threshold = np.load(result_path + '10_2.optimal_threshold_17dict.npy', allow_pickle=True).item()

RF = np.load(result_path + '2.MFI_RF.npy')
RF_autoFI = RF[0]
original_x_axis = np.load(result_path + '1.x_axis_channel.npy')
x_axis = [item.replace('-A', '') for item in original_x_axis]
channel_num = len(x_axis)
folder_path = result_path + '16.troubleshooting/'
os.makedirs(folder_path, exist_ok=True)

threshold_list = []
for key, val in condition4_optimal_threshold.items():
    threshold_list.append(val)
threshold = np.array(threshold_list)

fp_name_17 = ["01.EBFP2", "02.mTagBFP2", "03.mT-Sapphire", "04.mAmetrine", "05.mCerulean3", "06.LSSmOrange",
              "07.mBeRFP", "09.EGFP", "10.CyOFP1", "11.mClover3", "12.mVenus", "13.mPapaya", "14.mOrange2", "15.mRuby3",
              "16.mKate2", "17.mCardinal", "18.miRFP670"]

# define the control library
fp_num = len(threshold_list)
triangle_matrix = np.zeros((fp_num, fp_num), dtype=object)

combinations = []
for i in range(fp_num):
    for j in range(i+1, fp_num):
        combination = (fp_name_17[i], fp_name_17[j])
        triangle_matrix[i, j] = combination
        combinations.append(combination)
len_combinations = len(combinations)
print_lower_triangular_list(combinations)

ctrl_list = []
for idx, combo in enumerate(combinations):
    ctrl_list.append([idx, combo])
    print(f'Index {idx}: {combo}')
ctrl_array = np.array([(x, y) for x, y in ctrl_list], dtype=object)

# determine those problematic pMuSICs from Step22
problematic_pMuSICs_list = ['S7', 'S91', 'S92', 'S88', 'S65']

problematic_pMuSICs_x_scale_dict = {}
for key in condition4_all_pos_pMuSIC_scale_x.keys():
    for each in problematic_pMuSICs_list:
        match = re.match(r"(S\d+)_", key)
        if match:
            name = match.group(1)
        if each == name:
            problematic_pMuSICs_x_scale_dict[key] = condition4_all_pos_pMuSIC_scale_x[key]

problematic_pMuSICs_dict = {}
for key in condition4_pos_pMuSIC.keys():
    for each in problematic_pMuSICs_list:
        match = re.match(r"(S\d+)_", key)
        if match:
            name = match.group(1)
        if each == name:
            problematic_pMuSICs_dict[key] = condition4_pos_pMuSIC[key]
print(problematic_pMuSICs_dict.keys())

pMuSIC_combination_dict = {}
pMuSIC_combination_mfi_dict ={}
pMuSIC_combination_perc_dict = {}
for key, value in problematic_pMuSICs_x_scale_dict.items():
    print(key)
    total_cells = len(value)
    print("total_cells", total_cells)

    # to score each cell in each pMuSIC group
    result_list = []
    index_list = []
    for each_cell in value:
        score = each_cell[1:] / threshold

        max_val = max(score)
        max_idx = np.argmax(score)

        second_max_val = -float('inf')
        second_max_idx = -1

        for idx, val in enumerate(score):
            if val > second_max_val and val != max_val:
                second_max_val = val
                second_max_idx = idx

        most_likely_combination = (fp_name_17[max_idx], fp_name_17[second_max_idx])
        score_list = [0] * len_combinations
        index = find_combination_index(most_likely_combination, ctrl_array)
        index_list.append(index)
        score_list[index] = 1
        result_list.append(score_list)

    result = np.array(result_list)
    column_sums = np.sum(result, axis=0)
    percentage = np.round(column_sums / total_cells, decimals=4)
    possible_index = np.argmax(percentage)
    positive_rate = np.max(percentage)

    print("the index of the most likely pMuSIC is " + str(possible_index))
    print("the positive rate is " + str(positive_rate))
    pMuSIC_combination = ctrl_list[possible_index][1]
    print("the most likely combination of pR is ", pMuSIC_combination)

    # List categories with a percentage higher than 0.05
    categories = [(i, ctrl_list[i], perc) for i, perc in enumerate(percentage) if perc >= 0.03]
    print(categories)

    # Each listed category stands for a high-percentage combination in each pMuSIC.
    # get the cell index and their single_cell FI value across 48 channels
    pMuSIC_category_dict = {}
    pMuSIC_category_mfi_dict = {}
    pMuSIC_category_perc_dict = {}
    for category in categories:
        category_name = tuple(category[1])
        print(category_name, category[2])
        cell_index = [idx for idx, ind in enumerate(index_list) if ind == category[0]]

        # the cell_lindex stands for those cells share a same high-percentage combination such that we can get their FI
        # values across 48 channels
        cell_FI_list = []
        for each_ind in cell_index:
            cell_FI = problematic_pMuSICs_dict[key][each_ind]
            cell_FI_list.append(cell_FI)

        # here we got a cell array for each category with a shape (cell_number, 48)
        cell_FI_array = np.array(cell_FI_list)
        pMuSIC_category_dict[category_name] = cell_FI_array
        print(len(cell_FI_array))

        sample_MFI = []
        # this is the MFI of each category of each pMuSIC
        for i in range(channel_num):
            pure_FI = cell_FI_array - RF_autoFI
            col_data = pure_FI[:, i]
            pure_MFI = statistics.median(sorted(col_data))
            sample_MFI.append(np.round(pure_MFI, decimals=1))
        pMuSIC_category_mfi_dict[category_name] = sample_MFI
        pMuSIC_category_perc_dict[category_name] = category[2]
    pMuSIC_combination_mfi_dict[key] = pMuSIC_category_mfi_dict
    pMuSIC_combination_perc_dict[key] = pMuSIC_category_perc_dict

np.save(result_path + '16.troubleshooting_5pMuSICs_dict.npy', pMuSIC_combination_dict, allow_pickle=True)
np.save(result_path + '16.troubleshooting_5pMuSICs_MFI_dict.npy', pMuSIC_combination_mfi_dict, allow_pickle=True)
np.save(result_path + '16.troubleshooting_5pMuSICs_perc_dict.npy', pMuSIC_combination_perc_dict, allow_pickle=True)

rows = []
for key, subdict in sorted(pMuSIC_combination_perc_dict.items()):
    sorted_items = sorted(subdict.items(), key=lambda x: x[1], reverse=True)
    for idx, ((index, combination), perc) in enumerate(sorted_items):
        rows.append([key, index, combination, perc])
    rows.append([None, None, None, None])


df = pd.DataFrame(rows, columns=['Experiment', '136 Combo Indices', '136 Combinations', 'Percentage'])
print(df)

df.to_excel(folder_path + '16.summary of categories of each pMuSIC.xlsx', index=False)

ideal_pMuSICs_dict = {
    'S7 (mPapaya-mClover3)': RF[11] + RF[13],
    'S91 (CyOFP1-mRuby3)': RF[10] + RF[15],
    'S92 (mTagBFP2-mKate2)': RF[2] + RF[16],
    'S88 (mCardinal-mAmetrine)': RF[4] + RF[17],
    'S65 (CyOFP1-mKate2)': RF[10] + RF[16]
}

folder_path = result_path + '16.troubleshooting/problematic_pMuSICs_spectrum/'
os.makedirs(folder_path, exist_ok=True)


def extract_common_identifier(key):
    return key.split()[0].split('_')[0]


# to visualize the spectrum of each pMuSIC and its ideal version
for key1, val1 in ideal_pMuSICs_dict.items():
    sample1 = extract_common_identifier(key1)
    max_val1 = max(val1)
    normalized_sample1_MFI = val1/max_val1

    sample2_replicate = {}
    for key2, val2 in problematic_pMuSICs_dict.items():
        sample2 = extract_common_identifier(key2)
        if sample1 == sample2:
            pure_FI = val2 - RF_autoFI
            sample2_MFI = []
            for i in range(channel_num):
                col_data = pure_FI[:, i]
                pure_MFI = statistics.median(sorted(col_data))
                sample2_MFI.append(pure_MFI)
            max_val2 = max(sample2_MFI)
            normalized_sample2_mfi = sample2_MFI/max_val2
            sample2_replicate[key2] = normalized_sample2_mfi

    plt.figure(figsize=(9, 4))
    plt.plot(x_axis, normalized_sample1_MFI, label='Ideal_' + key1)

    if sample1 == 'S65':
        plt.tick_params('x', labelsize=11, rotation=90)
        plt.xlabel('Channel (Wavelength)', fontsize=18)
        plt.subplots_adjust(bottom=0.2)
    else:
        plt.xticks([])
        plt.xlabel('')
        plt.subplots_adjust(bottom=0.16)

    for sub_key, sub_val in sample2_replicate.items():
        plt.plot(x_axis, sub_val, label=sub_key)

    plt.ylabel('Normalized Intensity', fontsize=18)
    if sample1 == 'S88':
        plt.legend(loc='upper center', bbox_to_anchor=(0.42, 1), fontsize=12)
    else:
        plt.legend(loc='best',fontsize=12)

    plt.grid(False)

    file_path = folder_path + 'spectrum of problematic_pos_pMuSICs/'
    os.makedirs(file_path, exist_ok=True)
    filename = key1 + '.png'
    plt.savefig(file_path + filename, transparent=True)
    plt.close()


# to visualize the spectrum of each reference in a pMuSIC

reference_ideal_pMuSICs_dict = {
    'S7 (mPapaya-mClover3)': [RF[11], RF[13]],
    'S91 (CyOFP1-mRuby3)': [RF[10], RF[15]],
    'S92 (mTagBFP2-mKate2)': [RF[2], RF[16]],
    'S88 (mCardinal-mAmetrine)': [RF[4], RF[17]],
    'S65 (CyOFP1-mKate2)': [RF[10], RF[16]]
}


def extract_common_labels(key):
    labels = key.split('(')[1].split(')')[0].split('-')
    return labels[0], labels[1]


for key, val in reference_ideal_pMuSICs_dict.items():
    plt.figure(figsize=(9, 4))

    labels = extract_common_labels(key)
    plt.plot(x_axis, val[0], label=labels[0])
    plt.plot(x_axis, val[1], label=labels[1])

    plt.ylabel('Emission Intensity', fontsize=18)
    if 'S65' in key:
        plt.tick_params('x', labelsize=11, rotation=90)
        plt.xlabel('Channel (Wavelength)', fontsize=18)
        plt.subplots_adjust(bottom=0.2)
    else:
        plt.xticks([])
        plt.xlabel('')
        plt.subplots_adjust(bottom=0.16)

    plt.legend(loc='best', fontsize=12)
    plt.grid(False)

    file_path = folder_path + 'reference spectrum of problematic_pos_pMuSICs/'
    os.makedirs(file_path, exist_ok=True)
    filename = key + '.png'
    plt.savefig(file_path + filename, transparent=True)
    plt.close()



