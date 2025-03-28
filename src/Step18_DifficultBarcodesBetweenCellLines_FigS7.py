import os
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'Arial'

current_path = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_path)
result_path = os.path.join(project_root, 'output/')

loading_excel_path = os.path.join(project_root, 'csv_files/Spectra_data')
saving_path = os.path.join(project_root, 'figures_for_paper/',
                           'Fig.S7(Step18)Spectra of Difficult barcodes on HEK293T and MDA-MB-231 cells.png')

# Load the data for difficult barcodes from the Excel file
df_difficult_barcodes= pd.read_excel(os.path.join(loading_excel_path, 'Spectrum of difficult barcodes on both cell lines.xlsx'), sheet_name='difficult_barcodes')
original_columns1 = ['Barcode'] + [col for col in df_difficult_barcodes.columns if col not in ['Sample_Name', 'Barcode']]
df_difficult_barcodes['Barcode'] = df_difficult_barcodes['Sample_Name'].str.extract(r'([^\_]+)')
print(df_difficult_barcodes['Barcode'])

# Convert numeric columns to numeric type (coerce errors to NaN)
numeric_cols_barcodes1 = df_difficult_barcodes.columns.difference(['Sample_Name', 'Barcode'])
df_difficult_barcodes[numeric_cols_barcodes1] = df_difficult_barcodes[numeric_cols_barcodes1].apply(pd.to_numeric, errors='coerce')

# Average the data for difficult barcodes, grouped by Barcode
averaged_difficult_barcode = df_difficult_barcodes.groupby('Barcode')[numeric_cols_barcodes1].mean()
averaged_difficult_barcode = averaged_difficult_barcode[original_columns1[1:]]
pd.set_option('display.max_columns', 50)
print(averaged_difficult_barcode)
averaged_difficult_barcode.insert(0, 'Barcode', averaged_difficult_barcode.index)

# Load the data for good barcodes from the Excel file
df_good_barcodes= pd.read_excel(os.path.join(loading_excel_path, 'Spectrum of difficult barcodes on both cell lines.xlsx'), sheet_name='good_barcodes')
original_columns2 = ['Barcode'] + [col for col in df_good_barcodes.columns if col not in ['Sample_Name', 'Barcode']]
df_good_barcodes['Barcode'] = df_good_barcodes['Sample_Name'].str.extract(r'([^\_]+)')
print(df_good_barcodes['Barcode'])

# Convert numeric columns to numeric type (coerce errors to NaN)
numeric_cols_barcodes2 = df_good_barcodes.columns.difference(['Sample_Name', 'Barcode'])
df_good_barcodes[numeric_cols_barcodes2] = df_good_barcodes[numeric_cols_barcodes2].apply(pd.to_numeric, errors='coerce')

# Average the data for good barcodes, grouped by Barcode
averaged_good_barcode = df_good_barcodes.groupby('Barcode')[numeric_cols_barcodes2].mean()
averaged_good_barcode = averaged_good_barcode[original_columns2[1:]]
pd.set_option('display.max_columns', 50)
print(averaged_good_barcode)
averaged_good_barcode.insert(0, 'Barcode', averaged_good_barcode.index)

# Define lists of good and difficult barcode names for plotting
difficult_barcode_names = ['S7_(mClover3-mPapaya)', 'S91_(CyOFP1-mRuby3)', 'S92_(mTagBFP2-mKate2)']
good_barcode_names = ['S26_(mCardinal-mKate2)', 'S57_(mAmetrine-mVenus)', 'S85_(EGFP-mAmetrine)']

# ==================== XL Clemson ====================
# Create the subplots for good and difficult barcodes
# ==================== XL Clemson ====================

fig, axes = plt.subplots(3, 2, figsize=(8, 6))
for i, barcode in enumerate(good_barcode_names):
    ax = axes[i, 0]
    good_barcode_name = barcode.split('_')[0]  # Extract the barcode name before the underscore

    # Check if the barcode exists in the averaged data for good barcodes (for both HEK293T and MDA-MB-231 cells)
    if f'{good_barcode_name}-HEK293T' in averaged_good_barcode['Barcode'].tolist():
        barcode_HEK293T = averaged_good_barcode[averaged_good_barcode['Barcode'] == f'{good_barcode_name}-HEK293T'].iloc[:, 1:].values.flatten()
        if f'{good_barcode_name}-MDA-MB-231' in averaged_good_barcode['Barcode'].tolist():
            barcode_MDA231 = averaged_good_barcode[averaged_good_barcode['Barcode'] == f'{good_barcode_name}-MDA-MB-231'].iloc[:, 1:].values.flatten()

        # Plot the spectra for both cell lines
        ax.plot(barcode_HEK293T, label='HEK293T', color='orange')
        ax.plot(barcode_MDA231, label='MDA-MB-231', color='navy', linestyle=':')

        ax.set_title(barcode.split('(')[-1].split(')')[0])
        ax.set_xticks([])
        ax.legend(loc='best', framealpha=0, frameon=False)
    if i == 2:
        ax.set_xlabel('Channel (Wavelength) of Good Barcodes', fontsize=12)
    if i ==1:
        ax.set_ylabel('Emission Intensity', fontsize=14)

for i, barcode in enumerate(difficult_barcode_names):
    ax = axes[i, 1]
    difficult_barcode_name = barcode.split('_')[0]  # Extract the barcode name before the underscore

    # Check if the barcode exists in the averaged data for difficult barcodes (for both HEK293T and MDA-MB-231 cells)
    if f'{difficult_barcode_name}-HEK293T' in df_difficult_barcodes['Barcode'].tolist():
        barcode_HEK293T = averaged_difficult_barcode[averaged_difficult_barcode['Barcode'] == f'{difficult_barcode_name}-HEK293T'].iloc[:, 1:].values.flatten()
        if f'{difficult_barcode_name}-MDA-MB-231' in df_difficult_barcodes['Barcode'].tolist():
            barcode_MDA231 = averaged_difficult_barcode[averaged_difficult_barcode['Barcode'] == f'{difficult_barcode_name}-MDA-MB-231'].iloc[:, 1:].values.flatten()

        # Plot the spectra for both cell lines
        ax.plot(barcode_HEK293T, label='HEK293T', color='orange')
        ax.plot(barcode_MDA231, label='MDA-MB-231', color='navy', linestyle=':')

        ax.set_title(barcode.split('(')[-1].split(')')[0])
        ax.set_xticks([])
        ax.legend(loc='best', framealpha=0, frameon=False)

    if i == 2:
        ax.set_xlabel('Channel (Wavelength) of Difficult Barcodes', fontsize=12)

plt.tight_layout()
plt.savefig(saving_path, dpi=300, transparent=True)

print('Step18 Completed.')
