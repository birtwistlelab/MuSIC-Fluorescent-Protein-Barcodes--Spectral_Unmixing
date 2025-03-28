# ================================================ XL Clemson ================================================
# The following code is used to process pMuSIC unmixing results, comparing the inferred combinations with the actual
# combinations from sequencing data. The results are then visualized using a triangular heatmap for a better
# understanding of the relationships between the inferred and expected FP combinations.
# ================================================ XL Clemson ================================================

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams['font.family'] = 'Arial'

current_path = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_path)
result_path = os.path.join(project_root, 'output/')
loading_path = os.path.join(result_path, '13_3.Remove_mTFP1')

saving_path = os.path.join(result_path, '14.pMuSICsCalculateCorrectUnmixingPercentage_Remove_mTFP1')
os.makedirs(saving_path, exist_ok=True)

excel_path = os.path.join(loading_path, 'excel_summary/13_3.all_pos_pMuSIC_list(mTFP1_removed).xlsx')
df1 = pd.read_excel(excel_path)

ctrl_array = np.load(os.path.join(loading_path, '13_3.pMuSIC_control_array_remove_mTFP1.npy'), allow_pickle=True)
ctrl_lib = ctrl_array[:, 1].tolist()

fp_name = ['01.EBFP2', '02.mTagBFP2', '03.mT_Sapphire', '04.mAmetrine', '05.mCerulean3', '06.LSSmOrange', '07.mBeRFP',
           '09.EGFP', '10.CyOFP1', '11.mClover3', '12.mVenus', '13.mPapaya', '14.mOrange2', '15.mRuby3',
           '16.mKate2', '17.mCardinal', '18.miRFP670']

# Sort samples by 'Sample_Name', 'Dose', and 'Date'
df1['Sample_Name'] = df1['Name'].str.extract(r'S(\d+)_')[0].astype(int)
df1['Date'] = df1['Name'].str.extract(r'ng_(\d+)')[0]
df1['Dose'] = df1['Name'].str.extract(r'_(\d+)ng_')[0].astype(int)
df1_sorted = df1.sort_values(by=['Sample_Name', 'Dose', 'Date']).reset_index(drop=True)
df1_sorted['Name2'] = 'S' + df1_sorted['Sample_Name'].astype(str)
pd.set_option('display.max_columns', None)
print(df1_sorted)

# Load sequencing results to get the actual FP combinations
loading_excel_path2 = os.path.join(project_root, 'csv_files/Sequencing_Results_PlasmidSaurus/',
                                   'summary of sequencing results of pMuSICs.xlsx')
df2 = pd.read_excel(loading_excel_path2, dtype={'fp_m': str, 'fp_n': str})

# Extract and match the actual FP combinations from sequencing results
new_combinations = []
new_index = []

for j in range(len(df2['Name'])):
    if pd.notna(df2['fp_m'][j]):
        for n in range(len(fp_name)):
            if str(df2['fp_m'][j]) in fp_name[n]:
                new_fp_m = fp_name[n][3:]  # Get the FP name without the index number

            if str(df2['fp_n'][j]) in fp_name[n]:
                new_fp_n = fp_name[n][3:]

        if new_fp_m != new_fp_n:  # Only consider combinations with two different FPs
            for l in ctrl_lib:
                if new_fp_m in l and new_fp_n in l:
                    combo_index = ctrl_lib.index(l)
                    combo = l
            new_index.append(combo_index)
            new_combinations.append(combo)
        else:
            new_index.append(' ')
            new_combinations.append(' ')

# Update df2 with the actual sequencing information
df2['Name'] = 'S' + df2.iloc[:, 0].str.extract(r'pMuSIC_S(\d+)')[0]
df2.insert(3, 'ID_Sequencing', new_index)
df2.insert(4, 'Combination_Sequencing', new_combinations)
df2.to_excel(os.path.join(saving_path, '14.full_info_of_pMuSICs(actual)_remove_mTFP1.xlsx'), index=False)

# Check for duplicates in the samples tested and sequenced
print(df1_sorted['Sample_Name'].unique())
print('the number of MuSIC barcodes that we tested: ', len(df1_sorted['Sample_Name'].unique()))
# to provide the names of all pMuSICs listed in sequencing result
print(df2['Name'].unique())
print('the number of pMuSICs with no mTFP1 that we sent for sequencing', len(df2['Name']))


# Map the expected FP combinations (from sequencing) to the inferred results (from unmixing)
df1_sorted['Expected_ID'] = df1_sorted['Name2'].map(df2.set_index('Name')['ID_Sequencing'])
df1_sorted['Expected_Combination'] = df1_sorted['Name2'].map(df2.set_index('Name')['Combination_Sequencing'])

# Extract the inferred score for each pMuSIC based on the actual FP combination
Inferred_scores = []
for i in range(len(df1_sorted['Score_list'])):
  score_list = df1_sorted['Score_list'][i]
  score_list = score_list.strip('[]').strip()  # Remove the [] and the ' 's before and after the '[' and ']' to clean the score list string
  elements = score_list.split()  # to split each score using the ' ' inbetween scores
  new_row_list = [float(elem) for elem in elements]  # convert each score from the str type to float type
  expected_index = int(df1_sorted['Expected_ID'][i])
  inferred_score = new_row_list[expected_index]
  Inferred_scores.append(inferred_score)

df1_sorted['Inferred_Score'] = Inferred_scores
df1_sorted.to_excel(os.path.join(saving_path, '14.full_info_of_unmixed_pos_pMuSICs_remove_mTFP1(inferred).xlsx'), index=False)

# Prepare a summary dataframe that includes both actual and inferred results
df3 = pd.DataFrame({
    'Full_Name': df1_sorted['Name'],  # full name of each sample
    'Name':df1_sorted['Name2'],   # short name as 'S1'
    'Inferred_ID':df1_sorted['ID'],  # inferred ID from the unmixing results from Step13
    'Inferred_Combination':df1_sorted['Combination'],  # inferred combination from the unmixing results from Step13
    'Expected_ID':df1_sorted['Expected_ID'],  # actual ID from the sequencing results from PlasmidSaurus
    'Expected_Combination':df1_sorted['Expected_Combination'],  # actual combination from the sequencing results from PlasmidSaurus
    'Inferred_Score':df1_sorted['Inferred_Score']  # the inferred score of each pMuSIC under the actual combo index
})


# Average the inferred scores for replicates
df3['Sample_Number'] = df3['Name'].str.extract(r'S(\d+)')[0].astype(int)
df3_sorted = df3.sort_values(by='Sample_Number').reset_index(drop=True)
df3_avg = df3_sorted.groupby('Sample_Number')['Inferred_Score'].mean().reset_index()

# Add actual sequencing info to the averaged results
df3_avg['Name']= 'S' + df3_avg['Sample_Number'].astype(str)
df3_avg['Actual_ID'] = df3_avg['Name'].map(df2.set_index('Name')['ID_Sequencing'])

df3_avg = df3_avg.drop(columns=['Sample_Number'])
df3_avg = df3_avg[['Name','Actual_ID', 'Inferred_Score']]

combinations = np.load(os.path.join(loading_path, '13_3.pMuSIC_combination_list_remove_mTFP1.npy'))
num_samples_3 = len(df3_avg)
num_scores_3 = len(ctrl_lib)

score_list_AVG = []
actual_combo_AVG = []

# Process the average scores and assign to each sample
for i in range(num_samples_3):
  score_list_avg = [0] * num_scores_3  # initialize the average score list to 0 with a shape of (1 × N), where N = 136.
  avg_index = int(df3_avg['Actual_ID'][i])  # get index for actual sequencing ID in the ctrl_lib

  # Assign the inferred average score to its corresponding index in the score list.
  score_list_avg[avg_index] = df3_avg['Inferred_Score'][i]
  actual_combo = combinations[avg_index] # Get the actual combination from the sequencing data

  score_list_AVG.append(score_list_avg)
  actual_combo_AVG.append(actual_combo)
df3_avg['Actual_Combination'] = actual_combo_AVG
df3_avg['Score_list'] = score_list_AVG

df3_avg = df3_avg[['Name','Actual_ID', 'Actual_Combination', 'Inferred_Score',
                   'Score_list']]

print(df3_avg)
df3_avg.to_excel(os.path.join(saving_path, '14.Average of inferred scores of pMuSICs (Actual and inferred)_remove mTFP1.xlsx'), index=False)


# Group the pMuSIC sample by their expected combo index(actual combo136)
df3_fraction = df3_avg.sort_values(by=['Actual_ID']).reset_index(drop=True)
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
print(df3_fraction)

# Check if there is more than one samples share a same actual_combo_Index
# If yes, print out the index and put the average fraction under that actual_combo_Index
column_to_check = 'Actual_ID'
duplicates = df3_fraction[column_to_check][df3_fraction[column_to_check].duplicated()]

if not duplicates.empty:
    print(f"Column '{column_to_check}' contains duplicates:\n{duplicates}")
    # it will show as " the index.       the repeated value"

else:
    print(f"Column '{column_to_check}' has no duplicates.")

duplicate_values = duplicates.unique()
print('duplicate_values', duplicate_values)
# Output: duplicate_values [64 67 71 90 92 100 102 113 121 128]

duplicate_rows = df3_fraction[df3_fraction[column_to_check].isin(duplicate_values)]

print(f"Rows with duplicate values in column '{column_to_check}':\n{duplicate_rows}")
duplicate_indices = duplicate_rows.index

df3_fraction_avg = df3_fraction.groupby('Actual_ID')['Inferred_Score'].mean().reset_index()
df3_fraction_names = df3_fraction.groupby('Actual_ID')['Name'].apply(
    lambda x: sorted(x, key=lambda y: int(y[1:]))
).reset_index()
df3_fraction_avg['Name_group'] = df3_fraction_names['Name']

print('df3_fraction_avg')
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
print(df3_fraction_avg)

print(duplicate_indices)
print(df3_fraction.loc[duplicate_indices, 'Actual_ID'])

# Create the triangle matrix
fraction_list = [-1] * 136
for i in range(len(df3_fraction_avg['Actual_ID'])):
    fraction_index = int(df3_fraction_avg['Actual_ID'][i])  # to get the index of the position to insert the fp
    # to get the inferred score under the actual combo index that calculated above
    fraction_list[fraction_index] = df3_fraction_avg['Inferred_Score'][i]
    print(fraction_index, fraction_list[fraction_index])

n = 17  # Matrix size
triangle_matrix = np.full((n, n), np.nan)  # Initialize the matrix with NaNs

# Fill the lower triangle of the matrix (excluding the diagonal)
index = 0
for j in range(n):
    for i in range(j + 1, n):
        if index < len(fraction_list):
            if fraction_list[index] != -1:
                if fraction_list[index] < 0.999:
                    triangle_matrix[i, j] = round(fraction_list[index], 3)
                elif 0.999 <= fraction_list[index] < 1:
                    triangle_matrix[i, j] = 0.999
                else:
                    triangle_matrix[i, j] = 1

            else:
                triangle_matrix[i, j] = np.nan
            index += 1

# Create a DataFrame
df = pd.DataFrame(triangle_matrix)
fp_names = [fp_name[i][3:] for i in range(17)]  # remove the number before each FP

# Create a mask for the upper triangle and the diagonal
mask = np.triu(np.ones_like(triangle_matrix, dtype=bool))
np.fill_diagonal(mask, True)

# Plot the heatmap
plt.figure(figsize=(16, 16))
ax = sns.heatmap(df, annot=True, fmt=".3f", cmap="Blues", mask=mask, cbar=True, vmax=1,
            linewidths=0.5, linecolor='white', xticklabels=fp_names, yticklabels=fp_names, annot_kws={"size": 12},
            cbar_kws={'label': 'Fraction Abundance', 'shrink': 0.85, 'orientation': 'vertical', 'format': '%.2f',
                      'pad': 0.05}, square=True)
# Increase the colorbar label font size
cbar = plt.gca().collections[0].colorbar
cbar.ax.tick_params(labelsize=26)
cbar.set_label('Fractional Abundance', size=30, labelpad=15)

# Color the diagonal gray
for i in range(n):
    plt.gca().add_patch(plt.Rectangle((i, i), 1, 1, fill=True, facecolor='lightgray', edgecolor='gray'))

plt.yticks(fontsize=26)
plt.xticks(fontsize=26)

plt.subplots_adjust(left=0.2, bottom=0.2)
filename1 = '14.fraction_triangle_matrix_remove_mTFP1.png'
plt.savefig(os.path.join(saving_path, filename1), dpi=300, transparent=True)

fig_path = os.path.join(project_root, 'figures_for_paper/')
filename2 = 'Fig4D(Step14)_fraction_triangle_matrix_remove_mTFP1.png'
plt.savefig(os.path.join(fig_path, filename2), dpi=300, transparent=True)

plt.close()

print("================= Clemson XL =================")
print(fraction_list)
print(len(fraction_list))
np.save(os.path.join(saving_path,'14.fraction_list_heatmap.npy'), fraction_list)

df3_fraction_avg.to_excel(os.path.join(saving_path, '14.for triangle heatmap_remove_mTFP1.xlsx'), index=False)
