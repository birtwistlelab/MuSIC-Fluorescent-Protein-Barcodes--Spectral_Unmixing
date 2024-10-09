import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

plt.rcParams['font.family'] = 'DejaVu Sans'
project_root = os.path.dirname(os.path.abspath(__file__))
result_path = os.path.join(project_root, 'output/')

excel_path = 'Normalized_spectral_data_of_each_reference_CYTEK/Spectrum for Fig.2S.xlsx'
data = pd.read_excel(excel_path)
data['Name'] = data['Name'].astype(str)

original_x_axis = np.load(result_path + '1.x_axis_channel.npy')
x_axis = [item.replace('-A', '') for item in original_x_axis]
channel_num = len(x_axis)


def extract_plasmid_name(name):
    return name.split('pR-')[1].split('(')[0]


def extract_plasmid_date(name):
    return name.split('(')[1].split(')')[0]

data['Plasmid'] = data['Name'].apply(extract_plasmid_name)

unique_plasmids = data['Plasmid'].unique()

fig, axes = plt.subplots(6, 3, figsize=(12, 15))
axes = axes.flatten()


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
        normalized_mfi_values = row.iloc[1:49].values
        print(index)
        print(count)
        axes[i].plot(x_axis, normalized_mfi_values, color=colors[count], linewidth=linewidths[count],
                     linestyle=linestyles[count], label=extract_plasmid_date(row['Name']))
        count += 1
    axes[i].set_title(plasmid)
    axes[i].set_xticks([])
    # axes[i].set_ylabel('Normalized Median Fluorescence Intensity')
    axes[i].legend(loc='best')
    if i == 16:
        axes[i].set_xlabel('Channel (Wavelength)', fontsize=20, labelpad=10)

fig.supylabel('Normalized Median fluorescence Intensity', fontsize=20)
plt.tight_layout()

fig_path = os.path.join(project_root, 'paper_fig/')
filename = 'Fig.S2(Step17).png'
plt.savefig(fig_path + filename, transparent=True)
plt.close()

folder_path1 = result_path + '17.Reference_spectrum_for_fig.S2/'
os.makedirs(folder_path1, exist_ok=True)
filename1 = 'Fig.S2(Step17).png'
plt.savefig(folder_path1 + filename1, transparent=True)
plt.close()