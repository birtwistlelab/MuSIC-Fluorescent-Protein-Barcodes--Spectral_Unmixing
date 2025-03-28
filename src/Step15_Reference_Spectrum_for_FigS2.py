import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

plt.rcParams['font.family'] = 'Arial'
current_path = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_path)
result_path = os.path.join(project_root, 'output/')

saving_path = os.path.join(result_path, '15.Reference_Spectrum_for_FigS2')
os.makedirs(saving_path, exist_ok=True)

excel_path = os.path.join(project_root, 'csv_files/Spectra_data/Spectrum on HEK293T Cells_FigS2.xlsx')
data = pd.read_excel(excel_path)

# Extract 'Sample_Name' column as a string and assign it to a new column 'Name'
data['Name'] = data['Sample_Name'].astype(str)

original_x_axis = np.load(os.path.join(result_path, '1.x_axis_channel.npy'))
x_axis = [item.replace('-A', '') for item in original_x_axis]
channel_num = len(x_axis)


# Function to extract plasmid name from 'Name' column
def extract_plasmid_name(name):
    return name.split('pR-')[1].split('(')[0]

# Function to extract plasmid date from 'Name' column
def extract_plasmid_date(name):
    return name.split('(')[1].split(')')[0]


# Apply the plasmid extraction function to create a new column
data['Plasmid'] = data['Name'].apply(extract_plasmid_name)

unique_plasmids = data['Plasmid'].unique()


# Create a 6x3 grid for subplots (6 rows, 3 columns)
fig, axes = plt.subplots(6, 3, figsize=(12, 15))
axes = axes.flatten()

# Group the data by plasmid and get the number of samples for each plasmid
group_sizes = data.groupby('Plasmid').size()

for i, plasmid in enumerate(unique_plasmids):
    subset = data[data['Plasmid'] == plasmid]
    num_curves = group_sizes[plasmid]
    print('==============')
    print(plasmid)

    if num_curves == 3:
        colors = ['yellow', 'red', 'black']
        linewidths = [3, 1, 2]
        linestyles = ['-', '--', ':']

    else:
        linewidths = [1.5, 1.5]
        colors = ['yellow', 'black']
        linestyles = ['-', ':']
    count = 0
    for index, row in subset.iterrows():
        # Get the normalized MFI values from the current row (1st to 48th columns)
        normalized_mfi_values = row.iloc[1:49].values
        print(index)
        print(count)

        # Plot the spectra for the current plasmid
        axes[i].plot(x_axis, normalized_mfi_values, color=colors[count], linewidth=linewidths[count],
                     linestyle=linestyles[count], label=extract_plasmid_date(row['Name']))
        count += 1
    axes[i].set_title(plasmid)
    axes[i].set_xticks([])

    axes[i].legend(loc='best', frameon=False, framealpha=0)
    if i == 16:
        axes[i].set_xlabel('Channel (Wavelength)', fontsize=20, labelpad=10)

fig.supylabel('Normalized Median Fluorescence Intensity', fontsize=20)
plt.tight_layout()

fig_path = os.path.join(project_root, 'figures_for_paper/')
filename = 'Fig.S2(Step15)_Reference_Spectra.png'
plt.savefig(os.path.join(fig_path, filename), dpi=300, transparent=True)

plt.savefig(os.path.join(saving_path, filename), dpi=300, transparent=True)
plt.close()
