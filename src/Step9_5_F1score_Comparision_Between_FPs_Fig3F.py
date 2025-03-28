# ====================================== XL Clemson ======================================
# Calculate F1_scores under both original and IntensityCutoff at 10000
# FPR = FP / (FP =TN)
# ====================================== XL Clemson ======================================

import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
from my_module import scientific_formatter
plt.rcParams['font.family'] = 'Arial'

current_path = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_path)
result_path = os.path.join(project_root, 'output/')

loading_path_original = [os.path.join(result_path, '9_1.RangeIntensityCutoff', f'test{i + 1}/') for i in range(3)]
loading_path_cutoff = [os.path.join(result_path, '9_2.CalculationAtEachCutoff',  f'test{i + 1}/') for i in range(3)]

saving_path = [os.path.join(result_path, '9_5.F1score_Comparision_Between_FPs_Fig.3F', f'test{i + 1}/') for i in range(3)]
for i in range(3):
    os.makedirs(saving_path[i], exist_ok=True)

# Initialize dictionaries to store F1 scores for both original and cutoff data
original_F1_dict = {}
cutoff_F1_dict = {}

for i in range(3):
    # Load the original cell number of each fp
    df_original = pd.read_excel(os.path.join(loading_path_original[i], 'identification_result_original.xlsx'))
    original_cell_number = df_original.shape[0]

    # Load total cell numbers for cutoff = 1000 from the cutoff data
    total_cells_dict1 = np.load(os.path.join(loading_path_cutoff[i], 'total_cells_dict1.npy'), allow_pickle=True).item()
    cutoff1000_cell_number = total_cells_dict1[1000]

    # Check if the original and cutoff_1000 data have the same number of cells
    if original_cell_number == cutoff1000_cell_number:
        print(f'For test {i + 1}, original data == cutoff_1000 data')

    # use the results from cutoff_1000 data in loading_path_Cutoff as the original
    df_original_result = pd.read_excel(os.path.join(loading_path_cutoff[i], 'cutoff_1000.xlsx'))
    df_cutoff_result = pd.read_excel(os.path.join(loading_path_cutoff[i], 'cutoff_10000.xlsx'))

    # Extract FP names and F1 scores for both original and cutoff conditions
    fp_name = df_original_result['ref_fp'].str[3:].tolist()
    F1_score_original = df_original_result['F1_score']
    F1_score_cutoff = df_cutoff_result['F1_score']

    original_F1_dict[f'test{i + 1}'] = F1_score_original
    cutoff_F1_dict[f'test{i + 1}'] = F1_score_cutoff

    df = pd.DataFrame({
        'fp_Name': fp_name,
        'F1_Score_Original': F1_score_original,
        'F1_Score_Cutoff': F1_score_cutoff
    })
    df.to_excel(saving_path[i] + 'F1_scores.xlsx', index=False)

# Stack F1 scores for all tests and calculate average and standard deviation
original_f1 = np.vstack(list(original_F1_dict.values()))
average_original_f1 = np.mean(original_f1, axis=0)
std_original_f1 = np.std(original_f1, axis=0)

cutoff_f1 = np.vstack(list(cutoff_F1_dict.values()))
average_cutoff_f1 = np.mean(cutoff_f1, axis=0)
std_cutoff_f1 = np. std(cutoff_f1, axis=0)

# Create a DataFrame to store F1 scores and statistics for all tests
df2 = pd.DataFrame({
    'fp_name': fp_name,
    'original_test1_f1': original_f1[0, :],
    'original_test2_f1': original_f1[1, :],
    'original_test3_f1': original_f1[2, :],
    'cutoff_test1_f1': cutoff_f1[0, :],
    'cutoff_test2_f1': cutoff_f1[1, :],
    'cutoff_test3_f1': cutoff_f1[2, :],
    'average_original_f1': average_original_f1,
    'std_original_f1': std_original_f1,
    'average_cutoff_f1': average_cutoff_f1,
    'std_cutoff_f1': std_cutoff_f1
})


# Create the F1 score comparison bar plot
fp_name = df2['fp_name']
original_data = df2['average_original_f1']
cutoff_data = df2['average_cutoff_f1']
std_original = df2['std_original_f1']
std_cutoff = df2['std_cutoff_f1']

# ====================== XL Clemson =====================
# Plot the bar chart
# ====================== XL Clemson =====================
fig, ax = plt.subplots(figsize=(10, 4))
bar_width = 4
group_gap = 3
index = np.arange(len(fp_name)) * (2 * bar_width + group_gap)


# Plot the bars for original and cutoff data
bars1 = plt.bar(index, original_data, bar_width, label='No Intensity Cutoff', color='navy', alpha=0.8, edgecolor ='black')
bars2 = plt.bar(index + bar_width, cutoff_data, bar_width, label=f'Intensity Cutoff at {scientific_formatter(10000, None)}', color='white', alpha=0.8, edgecolor ='black')

# Add error bars to the plot
plt.errorbar(index, original_data, yerr=[np.zeros_like(std_original)], fmt='none', ecolor='black', capsize=2, capthick=1)
plt.errorbar(index + bar_width, cutoff_data, yerr=[np.zeros_like(std_cutoff)], fmt='none', ecolor='black', capsize=2, capthick=1)

plt.ylabel('F1 Score', fontsize=14)
ax.set_ylim([0.8, 1.1])
plt.yticks(np.arange(0.8, 1.1, 0.05))

ax.set_xticks(index + bar_width / 2)
ax.set_xticklabels(fp_name, rotation=45, ha='right', fontsize=12)
ax.tick_params(axis='y', labelsize=12)


plt.legend(framealpha=0, frameon=False, fontsize=14)
plt.subplots_adjust(bottom=0.3)
plt.savefig(os.path.join(result_path, '9_5.F1score_Comparision_Between_FPs_Fig.3F', 'F1 Score Comparison.png'), dpi=300, transparent=True)

plt.savefig(os.path.join(project_root, 'figures_for_paper/', 'Fig.3F(Step9_5)_F1Score_Comparison.png'), dpi=300, transparent=True)
