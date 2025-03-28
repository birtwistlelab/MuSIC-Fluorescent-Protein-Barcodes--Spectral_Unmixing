import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import statistics
from collections import defaultdict

plt.rcParams['font.family'] = 'Arial'
current_path = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_path)
result_path = os.path.join(project_root, 'output/')

loading_excel_path = os.path.join(result_path, '14.pMuSICsCalculateCorrectUnmixingPercentage_Remove_mTFP1')

# Load average scores and pMuSICs data from Excel files
average_score_file = pd.read_excel(os.path.join(loading_excel_path,
                                        '14.Average of inferred scores of pMuSICs (Actual and inferred)_remove mTFP1.xlsx'))
score_file = pd.read_excel(os.path.join(loading_excel_path, '14.full_info_of_unmixed_pos_pMuSICs_remove_mTFP1(inferred).xlsx'))
actual_barcode_file = pd.read_excel(os.path.join(loading_excel_path, '14.full_info_of_pMuSICs(actual)_remove_mTFP1.xlsx'))

saving_path = os.path.join(result_path, '16.Troubleshooting_FigS6')
os.makedirs(saving_path, exist_ok=True)

saving_path2 = os.path.join(project_root, 'figures_for_paper')
original_x_axis = np.load(os.path.join(result_path, '1.x_axis_channel.npy'))
x_axis = [item.replace('-A', '') for item in original_x_axis]

# Determine pMuSICs with perfect unmixing fraction (Inferred_Score = 1) and with bad unmixing fraction (Inferred_Score < 0.9)
bad_pMuSICs = average_score_file[average_score_file['Inferred_Score'] < 0.9][
    ['Name', 'Actual_Combination', 'Inferred_Score']].to_numpy()
good_pMuSICs = average_score_file[average_score_file['Inferred_Score'] == 1][
    ['Name', 'Actual_Combination', 'Inferred_Score']].to_numpy()

# Create DataFrame for bad and good pMuSICs
df_bad = pd.DataFrame({
    'Name': bad_pMuSICs[:, 0],
    'Actual_Combination': bad_pMuSICs[:, 1],
    'Inferred_Score': bad_pMuSICs[:, 2]
})
print('Difficult barcodes are:')
print(df_bad)

df_good = pd.DataFrame({
    'Name': good_pMuSICs[:, 0],
    'Actual_Combination': good_pMuSICs[:, 1],
    'Inferred_Score': good_pMuSICs[:, 2]
})
print('Good barcodes are:')
print(df_good)

# Check if there are additional barcodes that share the same combination as these difficult barcodes.
barcode_counts = average_score_file['Actual_Combination'].value_counts()
df_filtered = df_bad[df_bad['Actual_Combination'].map(barcode_counts).fillna(0) >= 2]
if not df_filtered.empty:
    print('Duplicate combinations found in the following pMuSIC samples:')
    print(df_filtered)
else:
    print('No duplicates found.')

# Find the barcode that shares a same combination with the filtered difficult barcode
matched_rows = average_score_file[average_score_file['Actual_Combination'].isin(df_filtered['Actual_Combination'])]
duplicate_pMuSIC = []
for each in matched_rows['Name'].tolist():
    if each not in df_filtered['Name'].values:
        print(f'The duplicate pMuSIC is {each}!')
        for i, row in score_file.iterrows():
            if each in row['Name2']:
                score = row['Inferred_Score']
        print(f'The score of {each} is {score}.')
        duplicate_pMuSIC.append([each, score])
print(duplicate_pMuSIC)

# Combine good and bad pMuSICs into a dictionary with their scores
bad_dict = {name: score for name, score in zip(df_bad['Name'], df_bad['Inferred_Score'])}
good_dict = {name: score for name, score in zip(df_good['Name'], df_good['Inferred_Score'])}
inferred_score_dict = {**bad_dict, **good_dict}
print(inferred_score_dict)

fp_name = ['01.EBFP2', '02.mTagBFP2', '03.mT_Sapphire', '04.mAmetrine', '05.mCerulean3', '06.LSSmOrange', '07.mBeRFP',
           '08.mTFP1', '09.EGFP', '10.CyOFP1', '11.mClover3', '12.mVenus', '13.mPapaya', '14.mOrange2', '15.mRuby3',
           '16.mKate2', '17.mCardinal', '18.miRFP670']

# Prepare actual barcode dictionary with corresponding combinations of fluorescent proteins
all_barcodes = list(inferred_score_dict.keys())
all_barcodes.append(duplicate_pMuSIC[0][0])
actual_barcode_dict = {}
for i, row in actual_barcode_file.iterrows():
    for barcode in all_barcodes:
        if row['Name'] == barcode:
            actual_1st_fp = []
            actual_2nd_fp = []
            for each_fp in fp_name:
                if str(row['fp_m']) in each_fp:
                    actual_1st_fp.append(each_fp)
                    print(actual_1st_fp)
            for each_fp in fp_name:
                if str(row['fp_n']) in each_fp:
                    actual_2nd_fp.append(each_fp)
            actual_combo = [actual_1st_fp[0], actual_2nd_fp[0]]
            actual_barcode_dict[barcode] = actual_combo

# Get the MFI of reference
RF = np.load(os.path.join(result_path, '2.MFI_RF.npy'))
auto_FI = RF[0]

# Create a dictionary of ideal pMuSICs for each combination of fluorophores
reference_ideal_pMuSICs_dict = {
    'S7_(mClover3-mPapaya)': [RF[11], RF[13]],
    'S91_(CyOFP1-mRuby3)': [RF[10], RF[15]],
    'S92_(mTagBFP2-mKate2)': [RF[2], RF[16]],
    'S26_(mCardinal-mKate2)': [RF[16], RF[17]],
    'S57_(mAmetrine-mVenus)': [RF[4], RF[12]],
    'S85_(EGFP-mAmetrine)': [RF[4], RF[9]]
}

# Create a dictionary of ideal spectra for each combination
ideal_pMuSICs_dict = {
    'S7_(mPapaya-mClover3)': RF[13] + RF[11],
    'S91_(CyOFP1-mRuby3)': RF[10] + RF[15],
    'S92_(mTagBFP2-mKate2)': RF[2] + RF[16],
    'S26_(mCardinal-mKate2)': RF[17] + RF[16],
    'S57_(mAmetrine-mVenus)': RF[4] + RF[12],
    'S85_(EGFP-mAmetrine)': RF[9] + RF[4]
}

# Normalize the ideal spectra
normalized_ideal_dict = {}
for key, val in ideal_pMuSICs_dict.items():
    max_val = max(val)
    normalized_sample_MFI = val / max_val
    normalized_ideal_dict[key] = normalized_sample_MFI

# Load pMuSICs data after removing mTFP1
pos_pMuSIC_remove_mTFP1_dict = np.load(os.path.join(result_path, '13_3.Remove_mTFP1/13_3.pos_pMuSICs_remove_mTFP1_dict.npy'), allow_pickle=True).item()


# Helper functions to extract plasmid and barcode names from pMuSICs
def extract_plasmid_name(name):
    return name.split('_')[0]


def extract_barcode_name(name):
    return name.split('(')[1].split(')')[0]


# Find and store problematic pMuSICs in the experiment dictionary
experiment_dict = {}
for key in pos_pMuSIC_remove_mTFP1_dict.keys():
    for k in ideal_pMuSICs_dict.keys():
        if extract_plasmid_name(k) == extract_plasmid_name(key):
            experiment_dict[key] = pos_pMuSIC_remove_mTFP1_dict[key]
        if extract_plasmid_name(key) == 'S59':
            experiment_dict[key] = pos_pMuSIC_remove_mTFP1_dict[key]

# Prepare for visualization
normalized_experiment_dict = {}
correct_rate = {}

sum_spectra = defaultdict(lambda: np.zeros(len(sample_MFI)))
key_counts = defaultdict(int)
for key, val in experiment_dict.items():
    pure_FI = val - auto_FI
    sample_MFI = []
    for i in range(val.shape[1]):
        col_data = pure_FI[:, i]
        pure_MFI = statistics.median(sorted(col_data))
        sample_MFI.append(pure_MFI)
    max_val = max(sample_MFI)
    normalized_sample_MFI = sample_MFI / max_val

    base_key = extract_plasmid_name(key)

    sum_spectra[base_key] += normalized_sample_MFI
    key_counts[base_key] += 1

    normalized_experiment_dict = {k: sum_spectra[k] / key_counts[k] for k in sum_spectra}


def extract_common_labels(name):
    labels = name.split('(')[1].split(')')[0].split('-')
    return list(labels)

# Generate plots comparing actual and ideal spectra
samples = list(normalized_ideal_dict.keys())

num_samples = len(samples)
fig, axes = plt.subplots(num_samples, 2, figsize=(14, 20))


for i, sample in enumerate(samples):
    ax_left = axes[i, 0]
    ax_right = axes[i, 1]
    name_ideal = extract_plasmid_name(sample)  # example: S7

    for key in reference_ideal_pMuSICs_dict.keys():
        if name_ideal in key:
            for j, spectra in enumerate(reference_ideal_pMuSICs_dict[key]):
                ax_left.plot(spectra, label=f"{extract_common_labels(key)[j]}")

    ax_left.set_ylabel('Emission Intensity', fontsize=16)

    if i == 2 or i == 3:
        ax_left.legend(loc='upper left', frameon=False, fontsize=12)
    else:
        ax_left.legend(loc='upper right', frameon=False, fontsize=12)

    # plot the ideal spectra
    ax_right.plot(normalized_ideal_dict[sample], label=f'Expected {extract_barcode_name(sample)}', linestyle="--", color='orangered')
    # plot the actual experimental data
    ax_right.plot(normalized_experiment_dict[name_ideal], label=f'Measured {extract_barcode_name(sample)}'
                                                                f'\nFractional Abundance: {inferred_score_dict[name_ideal]}',
                  linestyle="-", color='navy', alpha=0.6)
    if name_ideal == 'S91':
        ax_right.plot(normalized_experiment_dict['S59'], label=f'Measured mRuby3-CyOFP1 \nFractional Abundance: {duplicate_pMuSIC[0][1]}',
                      linestyle=":", color='green')

    ax_right.set_ylabel("Normalized Intensity", fontsize=16)
    ax_right.legend(loc='best', frameon=False, fontsize=12)

    if i == 5:
        ax_left.set_xlabel('Channel (Wavelength)', fontsize=16)
        ax_right.set_xlabel('Channel (Wavelength)', fontsize=16)

    ax_left.set_xticks([])
    ax_right.set_xticks([])

plt.savefig(os.path.join(saving_path, 'Fig.S6 troubleshooting.png'), dpi=300, transparent=True)
plt.savefig(os.path.join(saving_path2, 'Fig.S6(Step16)pMuSIC illustration.png'), dpi=300, transparent=True)
plt.close()
