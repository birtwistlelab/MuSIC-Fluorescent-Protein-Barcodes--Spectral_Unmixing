# To get the sequencing result list (actual) and make the new list with the expected index from the sequenced result
# list and the scores from the unmixing result list (inferred);
# To obtain the mean value of the scores if there are multiple replicates from one pMuSIC.

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams['font.family'] = 'Arial'

current_path = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_path)
result_path = os.path.join(project_root, 'output/')
loading_excel_path =  os.path.join(result_path,
                                   '11.pMuSICsUnmixing_IntensityCutoff/excel_summary/11.all_pos_pMuSIC_list.xlsx')
df1 = pd.read_excel(loading_excel_path)

ctrl_array = np.load(os.path.join(result_path, '11.pMuSIC_control_array.npy'), allow_pickle=True)
ctrl_lib = ctrl_array[:, 1].tolist()

fp_name = ['01.EBFP2', '02.mTagBFP2', '03.mT_Sapphire', '04.mAmetrine', '05.mCerulean3', '06.LSSmOrange', '07.mBeRFP',
           '08.mTFP1', '09.EGFP', '10.CyOFP1', '11.mClover3', '12.mVenus', '13.mPapaya', '14.mOrange2', '15.mRuby3',
           '16.mKate2', '17.mCardinal', '18.miRFP670']

# Get the sample name and sort them by ascending numerical order and the experiment date
df1['Sample_Name'] = df1['Name'].str.extract(r'S(\d+)_')[0].astype(int)
df1['Date'] = df1['Name'].str.extract(r'ng_(\d+)')[0]
df1['Dose'] = df1['Name'].str.extract(r'_(\d+)ng_')[0].astype(int)
df1_sorted = df1.sort_values(by=['Sample_Name', 'Dose', 'Date']).reset_index(drop=True)
df1_sorted['Name2'] = 'S' + df1_sorted['Sample_Name'].astype(str)
pd.set_option('display.max_columns', None)
print(df1_sorted)

# Extract the data from the summary list of sequencing result
loading_excel_path2 = os.path.join(project_root, 'csv_files/Sequencing_Results_PlasmidSaurus/',
                                   'summary of sequencing results of pMuSICs.xlsx')
df2 = pd.read_excel(loading_excel_path2, dtype={'fp_m': str, 'fp_n': str})

# Get those tested pMuSICs the actual combo info from those sequenced pMuSIC list:
new_combinations = []
new_index = []

for j in range(len(df2['Name'])):
    if pd.notna(df2['fp_m'][j]):
        for n in range(len(fp_name)):
            if str(df2['fp_m'][j]) in fp_name[n]:
                new_fp_m = fp_name[n][3:]

            if str(df2['fp_n'][j]) in fp_name[n]:
                new_fp_n = fp_name[n][3:]

        if new_fp_m != new_fp_n:  # only when fp_m and fp_n are different
            for l in ctrl_lib:
                if new_fp_m in l and new_fp_n in l:
                    combo_index = ctrl_lib.index(l)
                    combo = l
            new_index.append(combo_index)
            new_combinations.append(combo)
        else:
            new_index.append(' ')
            new_combinations.append(' ')

df2['Name'] = 'S' + df2.iloc[:, 0].str.extract(r'pMuSIC_S(\d+)')[0]
df2.insert(3, 'ID_Sequencing', new_index)
df2.insert(4, 'Combination_Sequencing', new_combinations)

saving_path = os.path.join(result_path, '12.pMuSICsCalculateCorrectUnmixingPercentage')
os.makedirs(saving_path, exist_ok=True)

df2.to_excel(os.path.join(saving_path, '12.full_info_of_pMuSICs(actual).xlsx'), index=False)

# Removes any duplicates of one pMuSIC and provides a list of all distinct pMuSICs that we tested
print(df1_sorted['Sample_Name'].unique())
print('the number of MuSIC barcodes that we tested: ', len(df1_sorted['Sample_Name'].unique()))

# Provide the names of all pMuSICs listed in sequencing result
print(df2['Name'].unique())
print('the number of pMuSICs that we sent for sequencing', len(df2['Name']))

# Add the actual combo index of pMuSICs and their combinations to the inferred unmixed result to verify if the
# unmixing of the tested pMuSICs aligns with the expected results from sequencing.

df1_sorted['Expected_ID'] = df1_sorted['Name2'].map(df2.set_index('Name')['ID_Sequencing'])
df1_sorted['Expected_Combination'] = df1_sorted['Name2'].map(df2.set_index('Name')['Combination_Sequencing'])

# Get the inferred score of each pMuSIC under the actual combo index

Inferred_scores = []
for i in range(len(df1_sorted['Score_list'])):
  score_list = df1_sorted['Score_list'][i]
  score_list = score_list.strip('[]').strip()  # Remove the [] and the ' 's before and after the '[' and ']'
  elements = score_list.split()  # Split each score using the ' ' inbetween scores
  new_row_list = [float(elem) for elem in elements]  # Convert each score from the str type to float type
  expected_index = int(df1_sorted['Expected_ID'][i])
  inferred_score = new_row_list[expected_index]
  Inferred_scores.append(inferred_score)

df1_sorted['Inferred_Score'] = Inferred_scores
print(df1_sorted)

df1_sorted.to_excel(os.path.join(saving_path, '12.full_info_of_unmixed_pos_pMuSICs(inferred).xlsx'), index=False)

# Get the brief info that only contains actual vs inferred
df3 = pd.DataFrame({
    'Full_Name': df1_sorted['Name'],  # Full name of each sample
    'Name':df1_sorted['Name2'],   # Short name as 'S1'
    'Inferred_ID':df1_sorted['ID'],  # Inferred ID from the unmixing results from Step11
    'Inferred_Combination':df1_sorted['Combination'],  # Inferred combination from the unmixing results from Step11
    'Expected_ID':df1_sorted['Expected_ID'],  # Actual ID from the sequencing results from PlasmidSaurus
    'Expected_Combination':df1_sorted['Expected_Combination'],  # Actual combination from the sequencing results from PlasmidSaurus
    'Inferred_Score':df1_sorted['Inferred_Score']  # The inferred score of each pMuSIC under the actual combo index
})


# Take average of the inferred score of the pMuSIC replicates and make a new list
df3['Sample_Number'] = df3['Name'].str.extract(r'S(\d+)')[0].astype(int)
df3_sorted = df3.sort_values(by='Sample_Number').reset_index(drop=True)
df3_avg = df3_sorted.groupby('Sample_Number')['Inferred_Score'].mean().reset_index()

df3_avg['Name']= 'S' + df3_avg['Sample_Number'].astype(str)
df3_avg['Actual_ID'] = df3_avg['Name'].map(df2.set_index('Name')['ID_Sequencing'])

print('Test1', df3_avg)
df3_avg = df3_avg.drop(columns=['Sample_Number'])
print('Test2', df3_avg)
df3_avg = df3_avg[['Name','Actual_ID', 'Inferred_Score']]
print('Test3', df3_avg)

combinations = np.load(result_path + '11.pMuSIC_combination_list.npy')
num_samples_3 = len(df3_avg['Name'])
num_scores_3 = len(ctrl_lib)

score_list_AVG = []
actual_combo_AVG = []
for i in range(num_samples_3):
  score_list_avg = [0] * num_scores_3  # Initialize the average score list to 0 with a shape of (1 × N), where N = 153.
  avg_index = int(df3_avg['Actual_ID'][i])  # Get index for actual sequencing ID in the ctrl_lib

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
df3_avg.to_excel(os.path.join(saving_path, '12.Average of inferred scores of pMuSICs (Actual and inferred).xlsx'), index=False)

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
    # It will show as " the index.       the repeated value"

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

# ================= XL Clemson =================
# Create the triangle matrix
# ================= XL Clemson =================

fraction_list = [-1] * 153
for i in range(len(df3_fraction_avg['Actual_ID'])):
    fraction_index = int(df3_fraction_avg['Actual_ID'][i])  # to get the index of the position to insert the fp
    # to get the inferred score under the actual combo index that calculated above
    fraction_list[fraction_index] = df3_fraction_avg['Inferred_Score'][i]
    print(fraction_index, fraction_list[fraction_index])

n = 18  # Matrix size
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
                    triangle_matrix[i, j] = 1.0

            else:
                triangle_matrix[i, j] = np.nan
            index += 1

# Create a DataFrame
df = pd.DataFrame(triangle_matrix)
fp_names = [fp_name[i][3:] for i in range(18)]  # remove the number before each FP

# Create a mask for the upper triangle and the diagonal
mask = np.triu(np.ones_like(triangle_matrix, dtype=bool))
np.fill_diagonal(mask, True)

# ================= XL Clemson =================
# Plot the heatmap
# ================= XL Clemson =================

plt.figure(figsize=(16, 16))
ax = sns.heatmap(df, annot=True, fmt=".3f", cmap="Blues", mask=mask, cbar=True, vmax=1,
            linewidths=0.5, linecolor='white', xticklabels=fp_names, yticklabels=fp_names, annot_kws={"size": 12},
            cbar_kws={"label": "Fraction Abundance", "shrink": 0.85, "orientation": "vertical", "format": "%.2f",
                      "pad": 0.05}, square=True)
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
filename1 = '12.fraction_triangle_matrix.png'
plt.savefig(os.path.join(saving_path, filename1),dpi=300, transparent=True)

fig_path = os.path.join(project_root, 'figures_for_paper/')
filename2 = 'Fig.4C(Step12)_fraction_triangle_martix.png'
plt.savefig(os.path.join(fig_path, filename2), dpi=300, transparent=True)

plt.close()

np.save(os.path.join(saving_path, '12.fraction_list_heatmap.npy'), fraction_list)

df3_fraction_avg.to_excel(os.path.join(saving_path, '12.for triangle heatmap.xlsx'), index=False)
