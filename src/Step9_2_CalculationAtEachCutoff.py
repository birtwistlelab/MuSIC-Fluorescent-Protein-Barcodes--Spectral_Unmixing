# ====================================== XL Clemson ======================================
# Second: Calculate OverallF1Score at each intensity cutoff
# ====================================== XL Clemson ======================================
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import time

plt.rcParams['font.family'] = 'Arial'
current_path = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_path)
result_path = os.path.join(project_root, 'output/')

# Load the average optimal thresholds
optimal_thresholds = np.load(
    os.path.join(result_path, '6.DetermineUnmixingThresholdsOnTrainingCells_FigS4/average_optimal_threshold.npy'))
thresholds = np.array(optimal_thresholds)

loading_excel_path = [os.path.join(result_path, '9_1.RangeIntensityCutoff/', f'test{i + 1}/') ### JRH changed from 'CutOff' to 'Cutoff'.
                      for i in range(3)]

saving_path = [os.path.join(result_path, '9_2.CalculationAtEachCutoff/', f'test{i + 1}/')
               for i in range(3)]
for path in saving_path:
    os.makedirs(path, exist_ok=True)

fp_name = ['01.EBFP2', '02.mTagBFP2', '03.mT_Sapphire', '04.mAmetrine', '05.mCerulean3', '06.LSSmOrange', '07.mBeRFP',
           '08.mTFP1', '09.EGFP', '10.CyOFP1', '11.mClover3', '12.mVenus', '13.mPapaya', '14.mOrange2', '15.mRuby3',
           '16.mKate2', '17.mCardinal', '18.miRFP670']

# Load cutoff threshold values for analysis
cutoff_thresholds = np.load(os.path.join(result_path, '9_1.cutoff_thresholds.npy'))
x_label = 'Unmixed Relative Abundance'
y_label = 'Frequency'

# Initialize dictionaries to store metrics for each cutoff threshold
precision_dict = {cutoff: [] for cutoff in cutoff_thresholds}
recall_dict = {cutoff: [] for cutoff in cutoff_thresholds}
f1_score_dict = {cutoff: [] for cutoff in cutoff_thresholds}
mcc_dict = {cutoff: [] for cutoff in cutoff_thresholds}


# Function to compute weighted F1 score
def compute_weighted_f1_score(f1_scores, total_cells_dict, cutoff):
    total_cells = sum(total_cells_dict[cutoff].values())
    weighted_f1 = 0  # Initialize weighted F1 score
    total_weight = 0  # Initialize total weight

    for i, f1_score in enumerate(f1_scores):
        ref_fp = fp_name[i]  # Get the fluorescent protein name
        weight = total_cells_dict[cutoff].get(ref_fp, 0) / total_cells  # Calculate weight based on cell count
        weighted_f1 += f1_score * weight  # Add weighted F1 score
        total_weight += weight  # Update total weight

    return weighted_f1  # Return the weighted F1 score


# Function to process each test case
def process_test(i):
    print(f'Processing Test {i + 1}...')
    start_time = time.time()  # Record start time for processing
    readable_start_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start_time))
    print(f"start_time of this test: {readable_start_time}")

    # Initialize dictionaries to store results for each test
    total_cells_dict1 = {}  # key = cutoff, val = total cells number of the testing pool
    total_cells_dict2 = {}  # key = cutoff, val = sub_dict {sub_key: fp name, sub_val: cell number}
    summary_dict = {}  # key = cutoff, val = [ref_fp, TP, FP, FN, TN, recall, precision, f1_score, mcc]
    overall_f1_dict = {}  # key = cutoff, val = overall_f1

    for j, cutoff in enumerate(cutoff_thresholds):
        df = pd.read_excel(os.path.join(loading_excel_path[i], f'cutoff_{cutoff}.xlsx'))

        # Extract the unmixed abundance array and calculate the total number of cells
        scale_x_array = np.vstack(df['unmixed_array'].apply(lambda x: np.fromstring(x[1:-1], sep=' ', dtype=float)))
        total_cells = scale_x_array.shape[0]
        print(f'total cells at {cutoff}:', total_cells)

        result = []  # Store results for each fluorescent protein
        f1_scores = []  # Store F1 scores for each fluorescent protein
        protein_cell_count_dict = {}  # key = fp_name, val = cell number for that fp
        for col, ref_fp in enumerate(fp_name):
            threshold = optimal_thresholds[col]  # Get the optimal threshold for the current fluorescent protein
            ref_values = scale_x_array[:, col]  # Get the intensity values for the current protein

            # Classify cells based on whether their intensity exceeds the threshold
            high_scale_x = ref_values >= threshold
            low_scale_x = ~high_scale_x

            # Get the actual and predicted labels for the current protein
            actual_is_FP = df['actual_fp'].str.contains(ref_fp)
            predicted_is_FP = df['predicted_fp'] == ref_fp

            # Calculate True Positives (TP), False Positives (FP), False Negatives (FN), True Negatives (TN)
            TP = np.sum(high_scale_x & actual_is_FP & predicted_is_FP)
            FP = np.sum(high_scale_x & ~actual_is_FP & predicted_is_FP)
            FN = np.sum(low_scale_x & actual_is_FP & ~predicted_is_FP)
            TN = np.sum(low_scale_x & ~actual_is_FP & ~predicted_is_FP)

            # Calculate precision, recall, F1 score, and MCC
            precision = TP / (TP + FP) if (TP + FP) > 0 else 0
            recall = TP / (TP + FN) if (TP + FN) > 0 else 0
            f1_score = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
            f1_scores.append(f1_score)

            # Calculate MCC
            mcc_numerator = (TP * TN) - (FP * FN)
            mcc_denominator = np.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
            mcc = mcc_numerator / mcc_denominator if mcc_denominator != 0 else 0

            result.append([ref_fp, TP, FP, FN, TN, recall, precision, f1_score, mcc])
            print(f'{ref_fp}: recall={recall}, precision={precision}, F1 Score={f1_score}, mcc = {mcc}')

            # Count the number of cells for the current fluorescent protein
            count = df['actual_fp'].str.contains(ref_fp).sum()
            protein_cell_count_dict[ref_fp] = count

        total_cells_dict1[cutoff] = total_cells
        total_cells_dict2[cutoff] = protein_cell_count_dict
        summary_dict[f'cutoff_{cutoff}'] = np.array(result, dtype=object)

        # Save the results to an Excel file
        df = pd.DataFrame(result, columns=['ref_fp', 'True Positive(TP)', 'False Positive(FP)',
                                           'False Negative(FN)', 'True Negative(TN)',
                                           'Recall', 'Precision', 'F1_score', 'MCC'])
        saving_path2 = os.path.join(saving_path[i], f'cutoff_{cutoff}.xlsx')
        df.to_excel(saving_path2, index=False)

        # Calculate weighted F1 score for each cutoff
        weighted_f1 = compute_weighted_f1_score(f1_scores, total_cells_dict2, cutoff)
        overall_f1_dict[cutoff] = weighted_f1
        print(f"Weighted F1 score at cutoff {cutoff}: {weighted_f1}")

    np.save(os.path.join(saving_path[i], 'overall_f1_dict.npy'), overall_f1_dict, allow_pickle=True)
    np.save(os.path.join(saving_path[i], 'result_summary_dict.npy'), summary_dict, allow_pickle=True)
    np.save(os.path.join(saving_path[i], 'total_cells_dict1.npy'), total_cells_dict1, allow_pickle=True)
    np.save(os.path.join(saving_path[i], 'total_cells_dict2.npy'), total_cells_dict2, allow_pickle=True)

    end_time = time.time()  # Record end time for processing
    readable_end_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end_time))
    print(f"end_time of this test: {readable_end_time}")
    print(f'Test {i + 1} Done!')


for i in range(3):
    process_test(i)

# Load and summarize results after all tests
for i in range(3):
    summary_path = os.path.join(result_path, f'9_2.CalculationAtEachCutoff/test{i + 1}/result_summary_dict.npy')
    summary_dict = np.load(summary_path, allow_pickle=True).item()

    for cutoff in cutoff_thresholds:
        results = summary_dict[f'cutoff_{cutoff}']
        recall_dict[cutoff].append([row[5] for row in results])
        precision_dict[cutoff].append([row[6] for row in results])
        f1_score_dict[cutoff].append([row[7] for row in results])
        mcc_dict[cutoff].append([row[8] for row in results])


# Compute the average values for each metric
def compute_average(metric_dict):
    avg_values = {cutoff: np.mean(metric_dict[cutoff], axis=0) for cutoff in cutoff_thresholds}
    df = pd.DataFrame.from_dict(avg_values, orient='index', columns=fp_name)
    df.insert(0, 'Cutoff Threshold', cutoff_thresholds)
    return df


# Calculate the average precision, recall, F1 score, and MCC for each cutoff threshold
precision_df = compute_average(precision_dict)
recall_df = compute_average(recall_dict)
f1_score_df = compute_average(f1_score_dict)
mcc_df = compute_average(mcc_dict)


# Check if all metrics exceed 0.95 and track the cutoff if they do
def check_greater_than(check_df, col_start, col_end, threshold_val, keyword):
    optimal_cutoff = None
    for index, row in check_df.iterrows():
        check_val = row.iloc[col_start: col_end]  # Get the subset of the row from col_start to col_end
        if (check_val > threshold_val).all():  # Check if all values are greater than threshold_val
            optimal_cutoff = check_df['Cutoff Threshold'][index]  # Get the cutoff value for this row
            print(f'For {keyword}, all fluorescent proteins at intensity cutoff {optimal_cutoff} have {keyword} > 0.95')
            break  # Exit the loop after finding the first matching row
    if optimal_cutoff is None:  # If no row met the condition
        print(f"No row meets the condition of all values > {threshold_val}")
        return None  # Return None to indicate no match
    return float(optimal_cutoff)  # Return the cutoff value for the row that met the condition


precision_cutoff = check_greater_than(precision_df, 1, 19, 0.95, 'precision')
recall_cutoff = check_greater_than(recall_df, 1, 19, 0.95, 'recall')
f1_cutoff = check_greater_than(f1_score_df, 1, 10, 0.95, 'f1 score')
mcc_cutoff = check_greater_than(mcc_df, 1, 19, 0.95, 'mcc')

# Determine the best cutoff threshold based on the metrics
best_cutoff = max(precision_cutoff, recall_cutoff, f1_cutoff, mcc_cutoff)
print(f'So, the best intensity cutoff should be above {best_cutoff}')

# Save metrics to an Excel file
saving_excel_path = os.path.join(result_path, '9_2.CalculationAtEachCutoff/average_metrics.xlsx')
with pd.ExcelWriter(saving_excel_path) as writer:
    precision_df.to_excel(writer, sheet_name='Precision', index=False)
    recall_df.to_excel(writer, sheet_name='Recall', index=False)
    f1_score_df.to_excel(writer, sheet_name='F1 Score', index=False)
    mcc_df.to_excel(writer, sheet_name='MCC', index=False)

