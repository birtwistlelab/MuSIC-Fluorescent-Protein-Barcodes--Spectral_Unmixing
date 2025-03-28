import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt


plt.rcParams['font.family'] = 'Arial'

current_path = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_path)
result_path = os.path.join(project_root, 'output/')

loading_excel_path = os.path.join(project_root, 'csv_files/Spectra_data')
saving_path = os.path.join(project_root, 'figures_for_paper/', 'Fig.S3(Step17)Spectra on HEK293T and MDA-MB-231 cells.png')

# Load the data for HEK293T and MDA-MB-231 cells from Excel files
df_HEK293T = pd.read_excel(os.path.join(loading_excel_path, 'Spectrum on HEK293T Cells_FigS2.xlsx'))
df_MDA231 = pd.read_excel(os.path.join(loading_excel_path, 'Spectrum on MDA-MB-231 Cells.xlsx'))

# Define the original columns, excluding 'Sample_Name' and 'Protein'
original_columns = ['Protein'] + [col for col in df_HEK293T.columns if col not in ['Sample_Name', 'Protein']]

# Extract 'Protein' names from 'Sample_Name' in HEK293T data and convert numerical columns to numeric types
df_HEK293T['Protein'] = df_HEK293T['Sample_Name'].str.extract(r'([^\(]+)')
numeric_cols_HEK293T = df_HEK293T.columns.difference(['Sample_Name', 'Protein'])
df_HEK293T[numeric_cols_HEK293T] = df_HEK293T[numeric_cols_HEK293T].apply(pd.to_numeric, errors='coerce')

# Average the data for each protein in HEK293T
averaged_HEK293T = df_HEK293T.groupby('Protein')[numeric_cols_HEK293T].mean()
averaged_HEK293T = averaged_HEK293T[original_columns[1:]]
averaged_HEK293T.insert(0, 'Protein', averaged_HEK293T.index)

# Extract 'Protein' names from 'Sample_Name' in MDA-MB-231 data and convert numerical columns to numeric types
df_MDA231['Protein'] = df_MDA231['Sample_Name'].str.extract(r'([^\(]+)')
numeric_cols_MDA231 = df_MDA231.columns.difference(['Sample_Name', 'Protein'])
df_MDA231[numeric_cols_MDA231] = df_MDA231[numeric_cols_MDA231].apply(pd.to_numeric, errors='coerce')

# Average the data for each protein in MDA-MB-231
averaged_MDA231 = df_MDA231.groupby('Protein')[numeric_cols_MDA231].mean()
averaged_MDA231 = averaged_MDA231[original_columns[1:]]
averaged_MDA231.insert(0, 'Protein', averaged_MDA231.index)

# Find common proteins between the two datasets (HEK293T and MDA-MB-231)
proteins = sorted(set(averaged_HEK293T.index) & set(averaged_MDA231.index))

# Create a 6x3 grid of subplots to visualize the spectra of the common proteins
fig, axes = plt.subplots(6, 3, figsize=(12, 15))
axes = axes.flatten()

for i, protein in enumerate(proteins):
    ax = axes[i]
    fp_HEK293T = averaged_HEK293T[averaged_HEK293T['Protein'] == protein].iloc[:, 1:].values.flatten()
    fp_MDA231 = averaged_MDA231[averaged_MDA231['Protein'] == protein].iloc[:, 1:].values.flatten()

    ax.plot(fp_HEK293T, label='HEK293T', color='orange')
    ax.plot(fp_MDA231, label='MDA-MB-231', color='navy', linestyle=':')

    axes[i].set_title(protein[6:])
    axes[i].set_xticks([])
    axes[i].legend(loc='best', framealpha=0, frameon=False)

    if i == 16:
        axes[i].set_xlabel('Channel (Wavelength)', fontsize=20, labelpad=10)

    ax.legend(framealpha=0, frameon=False)
fig.supylabel('Normalized Median Fluorescence Intensity', fontsize=20)
plt.tight_layout()

plt.subplots_adjust(bottom=0.05)
plt.savefig(saving_path, dpi=300, transparent=True)
