# to profile the 18 references to identify their spectrums and peak channels from the references
import os
import numpy as np
import statistics
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
plt.rcParams['font.family'] = 'DejaVu Sans'

current_path = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_path)
result_path = os.path.join(project_root, 'output/')

fp_reference = np.load(result_path + '1.fp_reference_2023.npy', allow_pickle=True).item()
x_axis = np.load(result_path + '1.x_axis_channel.npy')
channel_num = len(x_axis)

auto_FI = fp_reference['00.unstained']
noise = []
for i in range(channel_num):
    col_data = auto_FI[:, i]  # now we get the col_data which stands for all cell's FI per channel
    MFI = statistics.median(sorted(col_data))
    noise.append(MFI)
RF_auto_FI = np.array(noise)

fp_MFI_RF = {}
for key, val in fp_reference.items():
    pure_FI = val - RF_auto_FI
    RF_MFI = []
    for i in range(channel_num):
        col_data = pure_FI[:, i]
        pure_MFI = statistics.median(sorted(col_data))
        RF_MFI.append(np.round(pure_MFI, decimals=1))
    fp_MFI_RF[key] = np.array(RF_MFI)

fp_MFI_RF['00.unstained'] = RF_auto_FI

RF_name = []
for i in range(19):
    if i < 10:
        name = 'RF_fp0' + str(i)
    else:
        name = 'RF_fp' + str(i)
    RF_name.append(name)
print(RF_name)
# output of RF_name is ['RF_fp00', 'RF_fp01', 'RF_fp02', 'RF_fp03', 'RF_fp04', 'RF_fp05', 'RF_fp06', 'RF_fp07',
# 'RF_fp08', 'RF_fp09', 'RF_fp10', 'RF_fp11', 'RF_fp12', 'RF_fp13', 'RF_fp14', 'RF_fp15', 'RF_fp16', 'RF_fp17',
# 'RF_fp18']

RF = []
for each_name in RF_name:
    num = each_name[-2:]
    for key in fp_MFI_RF.keys():
        if num in key:
            print(key)
            rf_mfi = fp_MFI_RF[key]
            RF.append(rf_mfi)

# to visualize the spectrum of each Reference
RF_peak = {}
for i in range(1, 19):
    label = RF_name[i]
    each_ref = RF[i]
    height_threshold = round(max(each_ref) * 0.9, 1)
    peaks, _ = find_peaks(RF[i], height=height_threshold)
    RF_peak[label] = peaks

    plt.figure(figsize=(10, 6))
    plt.xlabel('Channels', fontsize=12)
    plt.plot(x_axis, each_ref, label=label)
    plt.tick_params('x', labelsize=8, rotation=90)
    plt.xlabel('Channels', fontsize=12)
    plt.ylabel('Emission Fluorescence Intensity', fontsize=12)
    plt.legend()
    plt.axhline(y=height_threshold, color='orange', linestyle='--')
    plt.subplots_adjust(bottom=0.14)
    plt.grid(True)

    file_path = result_path + '2.RF_spectrum/'
    os.makedirs(file_path, exist_ok=True)
    filename = label + '.png'
    plt.savefig(file_path + filename, transparent=True)
    plt.close()

print(RF_peak)
# output of RF_peak is : {'RF_fp01': array([2]), 'RF_fp02': array([2]), 'RF_fp03': array([4]), 'RF_fp04': array([6]),
# 'RF_fp05': array([4]), 'RF_fp06': array([7]), 'RF_fp07': array([9]), 'RF_fp08': array([16]), 'RF_fp09': array([16]),
# 'RF_fp10': array([20]), 'RF_fp11': array([17]), 'RF_fp12': array([18]), 'RF_fp13': array([18]),
# 'RF_fp14': array([30]), 'RF_fp15': array([30]), 'RF_fp16': array([33]), 'RF_fp17': array([33]),
# 'RF_fp18': array([41])}

np.save(result_path + '2.RF_peak.npy', RF_peak, allow_pickle=True)  # dictionary: 'RF_fp01': array([2])
np.save(result_path + '2.fp_MFI_RF.npy', fp_MFI_RF, allow_pickle=True)  # dictionary, pure_MFI of each reference
np.save(result_path + '2.MFI_RF.npy', RF)  # list, pure_MFI of each reference
np.save(result_path + '2.RF_name.npy', RF_name)  # list: 'RF_fp00'

# to visualize the peak range of each Reference with a threshold of 0.5
RF_peak = {}
for i in range(1, 19):
    label = RF_name[i]
    each_ref = RF[i]
    height_threshold = round(max(each_ref) * 0.3, 1)
    plt.figure(figsize=(10, 6))
    plt.xlabel('Channels', fontsize=12)
    plt.plot(x_axis, each_ref, label=label)
    plt.tick_params('x', labelsize=8, rotation=90)
    plt.xlabel('Channels', fontsize=12)
    plt.ylabel('Emission Fluorescence Intensity', fontsize=12)
    plt.legend()
    plt.axhline(y=height_threshold, color='orange', linestyle='--')
    plt.subplots_adjust(bottom=0.14)
    plt.grid(True)

    file_path = result_path + '2.RF_spectrum/main_peak_range/'
    os.makedirs(file_path, exist_ok=True)
    filename = label + '.png'
    plt.savefig(file_path + filename, transparent=True)
    plt.close()
