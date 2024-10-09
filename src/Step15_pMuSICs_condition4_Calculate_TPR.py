# to get the sequencing result list (actual) and make the new list with the expected index from the sequenced result
# list and the scores from the unmixing result list (inferred);
# do average of the scores if there are multiple replicates from one pMuSIC.

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams['font.family'] = 'DejaVu Sans'

current_path = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_path)
result_path = os.path.join(project_root, 'output/')

excel_path = result_path + '14.pMuSICs_condition4_unmixing/excel_summary/14.all_pos_pMuSIC_list_136combos.xlsx'
df1 = pd.read_excel(excel_path)

ctrl_array = np.load(result_path + '14.condition4_remove_mTFP1_pMuSIC_136_control_array.npy', allow_pickle=True)
ctrl_lib = ctrl_array[:, 1].tolist()

fp_name_remove_mTFP1 = ["01.EBFP2", "02.mTagBFP2", "03.mT-Sapphire", "04.mAmetrine", "05.mCerulean3", "06.LSSmOrange",
                        "07.mBeRFP", "09.EGFP", "10.CyOFP1", "11.mClover3", "12.mVenus", "13.mPapaya", "14.mOrange2",
                        "15.mRuby3","16.mKate2", "17.mCardinal", "18.miRFP670"]

# get the sample name and sort them by ascending numerical order and the experiment date
df1['Sample_Name'] = df1['Name'].str.extract(r'S(\d+)_')[0].astype(int)
df1['Date'] = df1['Name'].str.extract(r'ng_(\d+)')[0]
df1['Dose'] = df1['Name'].str.extract(r'_(\d+)ng_')[0].astype(int)
df1_sorted = df1.sort_values(by=['Sample_Name', 'Dose', 'Date']).reset_index(drop=True)
df1_sorted['Name2'] = 'S' + df1_sorted['Sample_Name'].astype(str)
pd.set_option('display.max_columns', None)
print(df1_sorted)

# to extract the data from the summary list of sequencing result which removing all mTPF1
excel_path2 = result_path + '14.pMuSICs_condition4_unmixing/remove_mTFP1_on_sequencing_result/14.pMuSICs after removing mTFP1 (actual).\
xlsx'
df2 = pd.read_excel(excel_path2, dtype={'fp_m': str, 'fp_n': str})


# to get those tested pMuSICs the actual 136combo info from those sequenced pMuSIC list
new_combinations = []
new_index = []
for j in range(len(df2['Name'])):
    if pd.notna(df2['fp_m'][j]):
        for n in range(len(fp_name_remove_mTFP1)):
            if str(df2['fp_m'][j]) in fp_name_remove_mTFP1[n]:
                new_fp_m = fp_name_remove_mTFP1[n]
            if str(df2['fp_n'][j]) in fp_name_remove_mTFP1[n]:
                new_fp_n = fp_name_remove_mTFP1[n]
        for l in ctrl_lib:
            if new_fp_m != new_fp_n: # only when fp_m and fp_n are different
                if new_fp_m in l and new_fp_n in l:
                    combo_index = ctrl_lib.index(l)
                    combo = l
        new_index.append(combo_index)
        new_combinations.append(combo)
    else:
        new_index.append(' ')
        new_combinations.append(' ')

df2['Name'] = 'S' + df2.iloc[:, 0].str.extract(r'pMuSIC_S(\d+)')[0]
df2.insert(3, 'ID_136combos', new_index)
df2.insert(4, 'Combination_136combos', new_combinations)

file_path = result_path + '15.pMuSICs_Condition4_TPR_plotting/'
os.makedirs(file_path, exist_ok=True)
df2.to_excel(file_path + '15.full_info_of_pMuSICs(actual)_136combos.xlsx', index=False)

# to removes any duplicates of one pMuSIC and provides a list of all distinct pMuSICs that we tested
print(df1_sorted['Sample_Name'].unique())
print('the number of MuSIC barcodes that we tested: ', len(df1_sorted['Sample_Name'].unique()))
# to provide the names of all pMuSICs listed in sequencing result
print(df2['Name'].unique())
print('the number of pMuSICs that we sent for sequencing',len(df2['Name']))

# to add the actual 136combo index of pMuSICs and their combinations to the inferred unmixed result to verify if the
# unmixing of the tested pMuSICs aligns with the expected results from sequencing.

df1_sorted['Expected_ID_136combos'] = df1_sorted['Name2'].map(df2.set_index('Name')['ID_136combos'])
df1_sorted['Expected_Combination_136combos'] = df1_sorted['Name2'].map(df2.set_index('Name')['Combination_136combos'])

# to get the inferred score of each pMuSIC under the actual 136combo index
Inferred_scores = []
for i in range(len(df1_sorted['Score_list_136combos'])):
  score_list = df1_sorted['Score_list_136combos'][i]
  score_list = score_list.strip('[]').strip()  # to remove the [] and the ' 's before and after the '[' and ']'
  elements = score_list.split()  # to split each score using the ' ' inbetween scores
  new_row_list = [float(elem) for elem in elements]  # convert each score from the str type to float type
  expected_index = int(df1_sorted['Expected_ID_136combos'][i])
  inferred_score = new_row_list[expected_index]
  Inferred_scores.append(inferred_score)

df1_sorted['Inferred_Score_136combos'] = Inferred_scores
print(df1_sorted)

file_path = result_path + '15.pMuSICs_Condition4_TPR_plotting/'
os.makedirs(file_path, exist_ok=True)

# now we could have the full info of the tested pMuSICs containing experiment date, plasmid dose for transfection,
# unmixing result(inferred) and sequencing result(actual)
df1_sorted.to_excel(file_path + '15.full_info_of_unmixed_pos_pMuSICs(inferred)_136combos.xlsx', index=False)

# to get the brief info that only contains actual vs inferred
df3 = pd.DataFrame({
    'Full_Name': df1_sorted['Name'],
    'Name':df1_sorted['Name2'],
    'Inferred_ID_136combos':df1_sorted['ID_136combos'],
    'Inferred_Combination_136combos':df1_sorted['Combination_136combos'],
    'Expected_ID_136combos':df1_sorted['Expected_ID_136combos'],
    'Expected_Combination_136combos':df1_sorted['Expected_Combination_136combos'],
    'Inferred_Score_136combos':df1_sorted['Inferred_Score_136combos']
})

# take average of the inferred score of the pMuSIC replicates and make a new list
df3['Sample_Number'] = df3['Name'].str.extract(r'S(\d+)')[0].astype(int)
df3_sorted = df3.sort_values(by='Sample_Number').reset_index(drop=True)
df3_avg = df3_sorted.groupby('Sample_Number')['Inferred_Score_136combos'].mean().reset_index()

df3_avg['Name']= 'S' + df3_avg['Sample_Number'].astype(str)
df3_avg['Actual_ID_136combos'] = df3_avg['Name'].map(df2.set_index('Name')['ID_136combos'])

print('Test1', df3_avg)
df3_avg = df3_avg.drop(columns=['Sample_Number'])
print('Test2', df3_avg)
df3_avg = df3_avg[['Name','Actual_ID_136combos', 'Inferred_Score_136combos']]
print('Test3', df3_avg)

combinations = np.load(result_path + '14.condition4_remove_mTFP1_pMuSIC_136_combination_list.npy')
num_samples_3 = len(df3_avg['Name'])
num_scores_3 = len(ctrl_lib)

score_list_AVG = []
actual_combo_AVG = []
for i in range(num_samples_3):
  score_list_avg = [0] * num_scores_3
  avg_index = int(df3_avg['Actual_ID_136combos'][i])  # actual index for 136 combo
  score_list_avg[avg_index] = df3_avg['Inferred_Score_136combos'][i]  # to get the average score under the actual index
  actual_combo = combinations[avg_index] # to get the actual combination from the sequencing data

  score_list_AVG.append(score_list_avg)
  actual_combo_AVG.append(actual_combo)
df3_avg['Actual_combination_136combos'] = actual_combo_AVG
df3_avg['Score_list_136combos'] = score_list_AVG

df3_avg = df3_avg[['Name','Actual_ID_136combos', 'Actual_combination_136combos', 'Inferred_Score_136combos',
                   'Score_list_136combos']]

print(df3_avg)
df3_avg.to_excel(file_path + '15.Average of inferred scores of pMuSICs (Actual and inferred) for 136 combos.xlsx',
                 index=False)


# group the pMuSIC sample by their expected combo index(actual combo136)
df3_tpr = df3_avg.sort_values(by=['Actual_ID_136combos']).reset_index(drop=True)
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
print(df3_tpr)

# check if there is more than one samples share a same actual_combo_Index
# If yes, print out the index and put the average TPR under that actual_combo_Index
column_to_check = 'Actual_ID_136combos'
duplicates = df3_tpr[column_to_check][df3_tpr[column_to_check].duplicated()]

if not duplicates.empty:
    print(f"Column '{column_to_check}' contains duplicates:\n{duplicates}")
    # it will show as " the index.       the repeated value"

else:
    print(f"Column '{column_to_check}' has no duplicates.")

duplicate_values = duplicates.unique()
print('duplicate_values', duplicate_values)
# output: duplicate_values [130 262]

duplicate_rows = df3_tpr[df3_tpr[column_to_check].isin(duplicate_values)]

print(f"Rows with duplicate values in column '{column_to_check}':\n{duplicate_rows}")
duplicate_indices = duplicate_rows.index

print('the expected index in the 136 combo is:')
print(df3_tpr['Actual_ID_136combos'][duplicate_indices])

df3_tpr_avg = df3_tpr.groupby('Actual_ID_136combos')['Inferred_Score_136combos'].mean().reset_index()
df3_tpr_names = df3_tpr.groupby('Actual_ID_136combos')['Name'].apply(
    lambda x: sorted(x, key=lambda y: int(y[1:]))
).reset_index()
df3_tpr_avg['Name_group'] = df3_tpr_names['Name']

print('df3_tpr_avg')
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
print(df3_tpr_avg)

print(duplicate_indices)
print(df3_tpr.loc[duplicate_indices, 'Actual_ID_136combos'])

# to create the triangle matrix
TPR_list_136 = [-1] * 136
for i in range(len(df3_tpr_avg['Actual_ID_136combos'])):
    TPR_index = int(df3_tpr_avg['Actual_ID_136combos'][i])  # to get the index of the position to insert the fp
    # to get the inferred score under the actual 136combo index that calculated from step15
    TPR_list_136[TPR_index] = df3_tpr_avg['Inferred_Score_136combos'][i]
    print(TPR_index, TPR_list_136[TPR_index])

fp_name = ["EBFP2", "mTagBFP2", "mT-Sapphire", "mAmetrine", "mCerulean3", "LSSmOrange", "mBeRFP", "EGFP", "CyOFP1",
           "mClover3", "mVenus", "mPapaya", "mOrange2", "mRuby3", "mKate2", "mCardinal", "miRFP670"]
n = 17  # Matrix size
triangle_matrix = np.full((n, n), np.nan)  # Initialize the matrix with NaNs

# Fill the lower triangle of the matrix (excluding the diagonal)
index = 0
for j in range(n):
    for i in range(j + 1, n):
        if index < len(TPR_list_136):
            if TPR_list_136[index] != -1:
                if TPR_list_136[index] < 0.999:
                    triangle_matrix[i, j] = round(TPR_list_136[index], 3)
                elif TPR_list_136[index] >= 0.999:
                    triangle_matrix[i, j] = 0.999

            else:
                triangle_matrix[i, j] = np.nan
            index += 1
# Create a DataFrame
df = pd.DataFrame(triangle_matrix)

# Create a mask for the upper triangle and the diagonal
mask = np.triu(np.ones_like(triangle_matrix, dtype=bool))
np.fill_diagonal(mask, True)

# Plot the heatmap
plt.figure(figsize=(16, 16))
ax = sns.heatmap(df, annot=True, fmt=".3f", cmap="Blues", mask=mask, cbar=True, vmax=1,
            linewidths=0.5, linecolor='white', xticklabels=fp_name, yticklabels=fp_name, annot_kws={"size": 12},
            cbar_kws={"label": "True Positive Rate", "shrink": 0.85, "orientation": "vertical", "format": "%.2f",
                      "pad": 0.05}, square=True)
# Increase the colorbar label font size
cbar = plt.gca().collections[0].colorbar
cbar.ax.tick_params(labelsize=26)
cbar.set_label('True Positive Rate', size=30, labelpad=15)

# Color the diagonal gray
for i in range(n):
    plt.gca().add_patch(plt.Rectangle((i, i), 1, 1, fill=True, facecolor='lightgray', edgecolor='gray'))

plt.yticks(fontsize=26)
plt.xticks(fontsize=26)
# plt.title('Intensity Cutoff and Remove mTFP1', fontsize=20)

plt.subplots_adjust(left=0.2, bottom=0.2)
filename1 = 'triangle matrix of condition4(Intensity cutoff and remove mTFP1).png'
plt.savefig(file_path + filename1, transparent=True)

fig_path = os.path.join(project_root, 'paper_fig/')
filename2 = 'Fig4C(Step15).png'
plt.savefig(fig_path + filename2, transparent=True)

plt.close()

print("================= Clemson XL =================")
print(TPR_list_136)
print(len(TPR_list_136))
np.save(result_path + '15.TPR_list_136(condition4)_of_17x17_heatmap.npy', TPR_list_136)

df3_tpr_avg.to_excel(file_path + '15.for triangle heatmap.xlsx', index=False)
