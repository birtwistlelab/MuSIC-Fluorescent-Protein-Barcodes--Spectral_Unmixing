# ====================================== XL Clemson ======================================
# Identify Each Cell in Testing Dictionary and Generate Fig.3C
# ====================================== XL Clemson ======================================
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams['font.family'] = 'Arial'

current_path = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_path)
result_path = os.path.join(project_root, 'output/')

loading_cell_path = [os.path.join(result_path, '4_2.SeparateTrainingAndTestingData/', f'test{i + 1}/') for i in
                     range(3)]
loading_scale_x_paths = [os.path.join(result_path, '7.UnmixingTransfectedTestingCells', f'test{i + 1}')
                         for i in range(3)]

saving_paths = [os.path.join(result_path, '8.IdentifyEachCellInTestingDict_Fig3C/', f'test{i + 1}/') for i in range(3)]
for path in saving_paths:
    os.makedirs(path, exist_ok=True)

saving_path2 = os.path.join(project_root, 'figures_for_paper/')

# Load optimal thresholds for classification
optimal_thresholds = np.load(os.path.join(result_path, '6.DetermineUnmixingThresholdsOnTrainingCells_FigS4/average_optimal_threshold.npy'))
thresholds = np.array(optimal_thresholds)

fp_name = ['01.EBFP2', '02.mTagBFP2', '03.mT_Sapphire', '04.mAmetrine', '05.mCerulean3', '06.LSSmOrange', '07.mBeRFP',
           '08.mTFP1', '09.EGFP', '10.CyOFP1', '11.mClover3', '12.mVenus', '13.mPapaya', '14.mOrange2', '15.mRuby3',
           '16.mKate2', '17.mCardinal', '18.miRFP670']

# =========================== XL Clemson =========================
# Main Loop to Analyze and Classify Cells in 3 Test Sets
# =========================== XL Clemson =========================

percentage_tests = []
for i in range(3):
    print(f'Test{i + 1}:')
    transfected_pR_fp_scale_x = np.load(os.path.join(loading_scale_x_paths[i], '7.transfected_Testing_pR_fp_scale_x.npy'),
        allow_pickle=True).item()

    testing_dict = np.load(os.path.join(loading_cell_path[i], 'testing_dict.npy'), allow_pickle=True).item()

    result_list = []  # Store final classification results

    output_path = os.path.join(saving_paths[i], 'unmixing_results.xlsx')
    with pd.ExcelWriter(output_path, engine='xlsxwriter') as writer:

        # Loop through each pR_fp and classify the unmixed results
        for key, val in transfected_pR_fp_scale_x.items():
            original_cells = testing_dict[key]
            total_cells = val.shape[0]

            result = []  # List to store the one-hot encoded classification results
            predicted_fp = []  # List to store the predicted FP labels for each cell

            # Classify each cell by comparing its unmixed fluorescence with thresholds
            for each_cell in val:
                score = each_cell / thresholds  # Normalize fluorescence intensities using optimal thresholds
                max_val = max(score)  # Identify the highest score
                max_idx = np.argmax(score)  # Identify the index of the most-likely FP
                most_likely_fp = fp_name[max_idx]  # Get the corresponding FP name
                predicted_fp.append(most_likely_fp)

                # Generate one-hot encoded classification result
                score_list = [0] * 18
                score_list[max_idx] = 1
                result.append(score_list)

            result = np.array(result)
            column_sums = np.sum(result, axis=0)  # Sum of correctly classified cells for each FP
            percentage = np.round(column_sums / total_cells, decimals=4)  # Calculate classification percentage
            result_list.append(percentage)

            df_unmixing = pd.DataFrame({
                'idx': range(total_cells),
                'original_array': list(original_cells),
                'unmixed_array': list(val),  # Ensure arrays are stored correctly
                'true_fp': key,
                'predicted_fp': predicted_fp
            })

            df_unmixing.to_excel(writer, sheet_name=key, index=False)  # Excel sheet names must be ≤31 chars

    result_array = np.array(result_list)
    percentage_tests.append(result_array)
    th = pd.DataFrame(result_array)
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    print('The percentage table is:')
    print(th)

    # Save classification percentages as Excel file
    th.to_excel(os.path.join(saving_paths[i], '8.original_percentage_18x18.xlsx'), index=False)

# ====================================== XL Clemson ======================================
# Calculate and Save Average Classification Percentage Across Test Sets
# ====================================== XL Clemson ======================================

array1 = percentage_tests[0]
array2 = percentage_tests[1]
array3 = percentage_tests[2]

# Calculate the average classification percentage across 3 test sets
average_array = np.mean([array1, array2, array3], axis=0)

df_average = pd.DataFrame(average_array).round(4)
df_average.to_excel(os.path.join(result_path, '8.IdentifyEachCellInTestingDict_Fig3C', '8.Average_Fraction.xlsx'),
                    index=False)

df = pd.DataFrame(average_array, columns=[f'Reference {i + 1}' for i in range(18)],
                  index=[f'Reference {i + 1}' for i in range(18)]).round(4)

# ========================= XL Clemson ========================
# Generate and Save the 18x18 Heatmap (Fig.3C)
# ========================= XL Clemson ========================

fig, ax = plt.subplots(figsize=(11, 11))

sns.heatmap(df, annot=True, fmt=".2f", cmap='Blues', mask=df.isnull(), cbar=True,
            linecolor='white', linewidths=1, xticklabels=fp_name, yticklabels=fp_name, vmax=1,
            annot_kws={"size": 10})

ax.set_aspect('equal')
colorbar = ax.collections[0].colorbar
colorbar.ax.tick_params(labelsize=18)
colorbar.set_label('Fraction of Testing Cells ', size=22, labelpad=15)

plt.ylabel('Actual', fontsize=22)
plt.xlabel('Inferred', fontsize=22)
plt.yticks(fontsize=18)
plt.xticks(fontsize=18)
plt.tight_layout()
plt.subplots_adjust(bottom=0.4)

filename = '8.FractionHeatmap(original).png'
plt.savefig(os.path.join(result_path, '8.IdentifyEachCellInTestingDict_Fig3C', filename), transparent=True)

filename2 = 'Fig.3C(Step8)_AverageFraction.png'
plt.savefig(os.path.join(saving_path2, filename2), dpi=300, transparent=True)
