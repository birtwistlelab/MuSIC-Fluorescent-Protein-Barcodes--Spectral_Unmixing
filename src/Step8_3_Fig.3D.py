# from step8_2, we knew that 08.mTFP1 is the most naughty one among all 18 fluorescent proteins.
# here, we will create a histogram with the false positive cell frequency as y-axis and the FI as the x-axis
# to see how different 08.mTFP1 from the others

import numpy as np
import math
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from my_module import scientific_formatter

plt.rcParams['font.family'] = 'DejaVu Sans'

current_path = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_path)
result_path = os.path.join(project_root, 'output/')
fig_path = os.path.join(project_root, 'paper_fig/')

pos_cell_dict = np.load(result_path + '4.pos_cell_dict_pR_fp_singlets.npy', allow_pickle=True).item()
FP_fp_dict = np.load(result_path + '8_2.FalsePositiveCells(pos_pR_fp_singlets).npy', allow_pickle=True).item()
RF_peak = np.load(result_path + '2.RF_peak.npy', allow_pickle=True).item()
fp_MFI_RF = np.load(result_path + '2.fp_MFI_RF.npy', allow_pickle=True).item()
RF_autoFI = fp_MFI_RF['00.unstained']

# to get the fluorescence intensity value at the peak channel of all FP cells
fp_dict = {}
for key, val in FP_fp_dict.items():
    print(key)
    print(len(val))

    for each_key in RF_peak.keys():
        if each_key[-2:] in key:
            peak_channel = RF_peak[each_key][0]
            print('peak_channel', peak_channel)

    fp_val = []
    for each_cell in val:
        pure_fl = each_cell[peak_channel] - RF_autoFI[peak_channel]
        fp_val.append(pure_fl)

    fp_dict[key] = np.array(fp_val)

print(fp_dict.keys())

# get the FI vals of the total cells
totalfp_cellList = [value for sublist in fp_dict.values() for value in sublist]
totalfp_cellArr = np.array(totalfp_cellList)
print('total cells', len(totalfp_cellArr))

# to discriminate mTFP1 from the other fps
mTFP1_cellArr = fp_dict.pop('08.mTFP1')
cellnumber_mTFP1 = len(mTFP1_cellArr)
print('false positive cells of mTFP1', cellnumber_mTFP1)

# to get all other false positive cells in one group
otherfp_cellList = [value for sublist in fp_dict.values() for value in sublist]
otherfp_cellArr = np.array(otherfp_cellList)
cellnumber_otherfp = len(otherfp_cellArr)
print('total false positive cells of other fps', cellnumber_otherfp)

# to create False Positive(FP) FP_cell histogram
log_value_fp = []
for each in totalfp_cellArr:
    log_value = math.log10(each)
    log_value_fp.append(round(log_value, 2))

log_value_fp = np.array(log_value_fp)

bins_num = round((max(log_value_fp) - min(log_value_fp)) * 20)
# print('bins_num', bins_num)
x_ticks = np.logspace(np.log10(min(totalfp_cellArr)), np.log10(max(totalfp_cellArr)), num=20)
# bin_width is an array containing each bin's right edge value
bin_width = np.logspace(np.log10(min(totalfp_cellArr)), np.log10(max(totalfp_cellArr)), num=bins_num)

fig, ax = plt.subplots(figsize=(6, 6))
ax.hist(mTFP1_cellArr, bins=bin_width, color='orange', alpha=0.6, label='mTFP1')
ax.hist(otherfp_cellArr, bins=bin_width, color='gray', alpha=0.6, label='all other fluorescent proteins')
ax.set_xlabel('Fluorescence Intensity', fontsize=16)
ax.set_ylabel('False positive cell frequency', fontsize=16)
ax.set_xlim(10 ** 2, 10 ** 7)
ax.set_ylim(0, 250)
plt.xscale('log')
plt.gca().xaxis.set_major_formatter(FuncFormatter(scientific_formatter))
ax.legend(fontsize=12, frameon=False)
ax.tick_params(axis='x', labelsize=12)
ax.tick_params(axis='y', labelsize=12)
plt.tight_layout()
plt.subplots_adjust(bottom=0.2)

file_path = result_path + '8.original_condition1/8_3_Fig.3D/'
os.makedirs(file_path, exist_ok=True)

filename = 'histogram of False positive cells(for fig3D).png'
plt.savefig(file_path + filename, transparent=True)

filename2 = 'Fig3D(Step8_3).png'
plt.savefig(fig_path + filename2, transparent=True)
