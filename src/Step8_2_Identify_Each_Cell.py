# After unmixing, 18 x-scale values for each cell were determined. Using the 18 thresholds from the ROC curves, we can
# calculate a score for each cell by dividing each x-scale value by its corresponding threshold, generating a score
# list of length 18. The highest score indicates the most likely fluorescent protein (FP) for that pR-FP-positive cell.
# Based on the index of the highest score, we can categorize cells into true positive (TP) and false positive (FP)
# groups. Within the FP group, we could further classify cells into four fluorescence intensity (FI) subclasses:
# low (10³–10⁴), medium (10⁴–10⁵), high (10⁵–10⁶), and top (>10⁶).
# By plotting a histogram with FI on the x-axis and false positive cell frequency on the y-axis, we can observe the
# distribution of positive cells across each FI subclass.

import numpy as np
import math
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from my_module import scientific_formatter

plt.rcParams['font.family'] = 'DejaVu Sans'
project_root = os.path.dirname(os.path.abspath(__file__))
result_path = os.path.join(project_root, 'output/')

optimal_threshold = np.load(result_path + '7.optimal_threshold_dict.npy', allow_pickle=True).item()
pos_pR_fp_scale_x = np.load(result_path + '5.pos_pR_fp_scale_x.npy', allow_pickle=True).item()
pos_cell_dict = np.load(result_path + '4.pos_cell_dict_pR_fp_singlets.npy', allow_pickle=True).item()
RF_peak = np.load(result_path + '2.RF_peak.npy', allow_pickle=True).item()

fp_MFI_RF = np.load(result_path + '2.fp_MFI_RF.npy', allow_pickle=True).item()
RF_autoFI = fp_MFI_RF['00.unstained']

fp_name = ["01.EBFP2", "02.mTagBFP2", "03.mT-Sapphire", "04.mAmetrine", "05.mCerulean3", "06.LSSmOrange", "07.mBeRFP",
           "08.mTFP1", "09.EGFP", "10.CyOFP1", "11.mClover3", "12.mVenus", "13.mPapaya", "14.mOrange2", "15.mRuby3",
           "16.mKate2", "17.mCardinal", "18.miRFP670"]

classes = [10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7]
labels = ['Low FI (10^3-10^4)', 'Mid FI (10^4-10^5)', 'High FI (10^5-10^6)', 'Top FI (>10^6)']
colors = ['navy', 'orange', 'gray', 'green']

threshold_list = []
for key, val in optimal_threshold.items():
    threshold_list.append(val)
threshold = np.array(threshold_list)

Score_lists = []
TP_fp_dict = {}
FP_fp_dict = {}

fp_val_dict = {}

fp_category = []  # add up all false positive cells into one category and plotting
for key, val in sorted(pos_pR_fp_scale_x.items()):
    data = val[:, 1:]
    total_cells = len(data)
    print("total cells of " + key, total_cells)

    result_list = []
    cell_idx = []

    for each_key in RF_peak.keys():
        if each_key[-2:] in key:
            peak_channel = RF_peak[each_key][0]
            print('peak_channel', peak_channel)

    for each_cell in data:
        score = each_cell / threshold

        max_val = max(score)
        max_idx = np.argmax(score)

        most_likely_fp = fp_name[max_idx]
        score_list = [0] * 18
        score_list[max_idx] = 1
        result_list.append(score_list)
        cell_idx.append(max_idx)

    result = np.array(result_list)
    column_sums = np.sum(result, axis=0)
    percentage = np.round(column_sums / total_cells, decimals=4)
    Score_lists.append(percentage)

    max_idx_group = np.argmax(percentage)

    # to separate the true positive pR-fp cells and the false positive pR-fp cells
    tp_pos_fp = []
    fp_pos_fp = []
    for i in range(len(cell_idx)):
        if cell_idx[i] == max_idx_group:
            tp_pos_fp.append(pos_cell_dict[key][i])
        else:
            fp_pos_fp.append(pos_cell_dict[key][i])
    print('false positive cell number of ', len(fp_pos_fp))

    tp_pos_fp = np.array(tp_pos_fp)
    fp_pos_fp = np.array(fp_pos_fp)

    TP_fp_dict[key] = tp_pos_fp
    FP_fp_dict[key] = fp_pos_fp

    fp_val = []
    for each_cell in fp_pos_fp:
        pure_fl = each_cell[peak_channel] - RF_autoFI[peak_channel]
        fp_val.append(pure_fl)
    fp_val = np.array(fp_val)
    fp_val_dict[key] = fp_val

    # to create FP_cell histogram
    log_value_fp = []
    for each in fp_val:
        log_value = math.log10(each)
        log_value_fp.append(round(log_value, 2))
    log_value_fp = np.array(log_value_fp)
    bins_num = round((max(log_value_fp) - min(log_value_fp)) * 20)
    # print('bins_num', bins_num)
    x_ticks = np.logspace(np.log10(min(fp_val)), np.log10(max(fp_val)), num=20)

    # bin_width is an array containing each bin's right edge value
    bin_width = np.logspace(np.log10(min(fp_val)), np.log10(max(fp_val)), num=bins_num)

    fig, ax = plt.subplots()
    for i in range(len(classes) - 1):
        mask = (fp_val >= classes[i]) & (fp_val < classes[i + 1])
        ax.hist(fp_val[mask], bins=bin_width, color=colors[i], alpha=0.6, label=labels[i])

    ax.legend()
    ax.set_xlabel('Fluorescence Intensity of ' + key)
    ax.set_ylabel('False positive cell frequency')
    ax.set_xlim(10 ** 2, 10 ** 7)
    ax.set_ylim(0, 120)
    plt.xscale('log')
    plt.gca().xaxis.set_major_formatter(FuncFormatter(scientific_formatter))

    file_path = result_path + '8.original_condition1/8_2_False_Positive_cells/'
    os.makedirs(file_path, exist_ok=True)
    filename = 'False positive cells of ' + key + '.png'
    plt.savefig(file_path + filename)
    plt.close()

np.save(result_path + '8_2.FalsePositiveCells(pos_pR_fp_singlets).npy', FP_fp_dict, allow_pickle=True)
