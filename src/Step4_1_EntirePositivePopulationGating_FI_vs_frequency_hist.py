# ================================================= XL Clemson =================================================
# To plot the histograms of fluorescence intensity vs frequency and gate the transfected population of pR_fp_singlets
# ================================================= XL Clemson =================================================

import numpy as np
import os
import matplotlib.pyplot as plt
import math
from matplotlib.ticker import FuncFormatter
from my_module import scientific_formatter
plt.rcParams['font.family'] = 'Arial'

current_path = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_path)
result_path = os.path.join(project_root, 'output/')

pR_fp_singlets = np.load(os.path.join(result_path, '3_pR_fp_singlets.npy'), allow_pickle=True).item()
fp_reference = np.load(os.path.join(result_path, '1.fp_reference_2023.npy'), allow_pickle=True).item()
RF_peak = np.load(os.path.join(result_path, '2.RF_peak.npy'), allow_pickle=True).item()
RF = np.load(os.path.join(result_path, '2.MFI_RF.npy'))
RF_name = np.load(os.path.join(result_path, '2.RF_name.npy'))
x_axis = np.load(os.path.join(result_path, '1.x_axis_channel.npy'))
channel_num = len(x_axis)

# Get a new dict with the same key of pR_fp_singlets
pR_peak = {}
for k1 in RF_peak.keys():
    for k2 in pR_fp_singlets.keys():
        if k1[-2:] == k2[:2]:
            pR_peak[k2] = RF_peak[k1][0]

# Get the autofluorescence
auto_FI = fp_reference['00.unstained']

# Gate the transfected cells using unstained cells expressing pReceiver:
transfected_cell_dict = {}
transfected_cell_ind_dict = {}
negative_control_thresholds =[]

for key, val in sorted(pR_fp_singlets.items()):
    total_cells = val.shape[0]

    peak_channel = pR_peak[key]
    values = val[:, peak_channel]
    unstained = auto_FI[:, peak_channel]

    transfected_cells = []
    transfected_cell_list = []
    cell_list = values.tolist()

    neg_threshold = np.percentile(unstained, 100)
    negative_control_thresholds.append(neg_threshold)

    for each_val in cell_list:
        each_index = cell_list.index(each_val)
        if each_val > neg_threshold:
            transfected_cell_list.append(each_index)
            transfected_cells.append(val[each_index])
    transfected_cell_num = len(transfected_cell_list)
    transfected_cell_rate = transfected_cell_num/total_cells * 100

    transfected_cell_dict[key] = np.array(transfected_cells)
    transfected_cell_ind_dict[key] = transfected_cell_list

np.save(os.path.join(result_path, '4_1.transfected_cell_dict_pR_fp_singlets.npy'), transfected_cell_dict, allow_pickle=True)
np.save(os.path.join(result_path, '4_1.transfected_cell_ind_dict_pR_fp_singlets.npy'), transfected_cell_ind_dict, allow_pickle=True)
np.save(os.path.join(result_path, '4_1.negative_control_thresholds.npy'), negative_control_thresholds)

# Plot the histogram of fluorescence intensity (FI) versus frequency to visually demonstrate the gating process for
# transfected cells. Since some FI values are negative in the original data exported from FlowJo, we adjust them to
# positive values for logarithmic calculations. This adjustment is solely for visualization purposes—only the raw,
# unmodified data (above) is used for gating transfected cells.

for key, val in sorted(pR_fp_singlets.items()):
    total_cells = val.shape[0]

    peak_channel = pR_peak[key]
    values = val[:, peak_channel]

    # Ensure all values are greater than 0 for subsequent logarithmic operations.
    pure_value = values - np.min(values) + 0.1
    unstained = auto_FI[:, peak_channel] - np.min(auto_FI[:, peak_channel]) + 0.1

    log_value_list = []
    for each in pure_value:
        log_value = math.log10(each)
        log_value_list.append(round(log_value, 2))
    log_value_array = np.array(log_value_list)

    bins_num = round((max(log_value_array) - min(log_value_array)) * 30)
    # print('bins_num', bins_num)
    x_ticks = np.logspace(np.log10(min(pure_value)), np.log10(max(pure_value)), num=20)

    # bin_width here is an array containing each bin's right edge value
    bin_width = np.logspace(np.log10(min(pure_value)), np.log10(max(pure_value)), num=bins_num)

    fig = plt.figure(figsize=(6, 4))
    ax = plt.subplot2grid((3, 3), (0, 0), colspan=3, rowspan=3)

    hist, bins, _ = plt.hist(pure_value, bins=bin_width, color='green', alpha=0.5, label=key)
    # Plot the second histogram without bars (only calculate)
    hist2, bins2 = np.histogram(unstained, bins=bin_width)

    # Make sure the main neg peak of unstained cells close to the main neg peak of sample cells
    max_hist_sample = np.argmax(hist)
    max_hist2_unstained = np.argmax(hist2)
    max_hist_x = (bins[max_hist_sample] + bins[max_hist_sample + 1]) / 2

    max_hist2_x = (bins2[max_hist2_unstained] + bins2[max_hist2_unstained + 1]) / 2
    print(max_hist_x, max_hist2_x)
    distance = max_hist_x - max_hist2_x

    correct_unstained = unstained + distance

    neg_threshold = np.percentile(correct_unstained, 100)

    plt.hist(correct_unstained, bins=bin_width, color='gray', alpha=0.4, label='pR cells')

    plt.xscale('log')

    plt.gca().xaxis.set_major_formatter(FuncFormatter(scientific_formatter))
    plt.xlabel('Fluorescence intensity', fontsize=9)
    plt.ylabel('Frequency', fontsize=9)
    plt.xticks(fontsize=9)
    plt.yticks(fontsize=9)
    plt.axvline(x=neg_threshold, color='black', linestyle='--', label='negative threshold')

    plt.xlim(100)
    plt.legend(loc='upper right', fontsize=9)
    plt.tight_layout()

    file_path = result_path + '4.FI_histogram/'
    os.makedirs(file_path, exist_ok=True)
    filename = key + '.png'
    plt.savefig(file_path + filename)
