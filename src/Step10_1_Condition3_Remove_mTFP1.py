# To remove mTFP1, we need to remove mTFP1 from reference, re-perform unmixing, and re-create ROC curves for thresholds.

import numpy as np
from scipy.optimize import nnls
import matplotlib.pyplot as plt
import os

plt.rcParams['font.family'] = 'DejaVu Sans'

current_path = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_path)
result_path = os.path.join(project_root, 'output/')

pos_cell_dict = np.load(result_path + '4.pos_cell_dict_pR_fp_singlets.npy', allow_pickle=True).item()
pR_fp_singlets = np.load(result_path + '3.pR_fp_singlets.npy', allow_pickle=True).item()
RF_peak = np.load(result_path + '2.RF_peak.npy', allow_pickle=True).item()
RF = np.load(result_path + '2.MFI_RF.npy')
RF_name = np.load(result_path + '2.RF_name.npy')

# delete the positive-population of 08.mTFP1 in the dictionary
pos_cell_dict_17 = pos_cell_dict
del pos_cell_dict_17['08.mTFP1']

# delete the reference of 08.mTFP1 in the RF
RF_17 = RF
RF_17 = np.delete(RF_17, 8, axis=0)

# delete the 08.mTFP1 in RF_name
RF_name_17 = RF_name
RF_name_17 = np.delete(RF_name_17, 8, axis=0)

# delete the 08.mTFP1 in RF_peak
RF_peak_17 = RF_peak
del RF_peak_17['RF_fp08']

x_axis = np.load(result_path + '1.x_axis_channel.npy')
channel_num = len(x_axis)

row_labels = RF_name_17[1:]
col_labels = ["01.EBFP2", "02.mTagBFP2", "03.mT-Sapphire", "04.mAmetrine", "05.mCerulean3", "06.LSSmOrange",
              "07.mBeRFP", "09.EGFP", "10.CyOFP1", "11.mClover3", "12.mVenus", "13.mPapaya", "14.mOrange2",
              "15.mRuby3", "16.mKate2", "17.mCardinal", "18.miRFP670"]
x_label = 'Unmixed Relative Abundance'
y_label = 'Frequency'

RF = np.array(RF_17).T
pos_pR_fp_scale_x_17 = {}
for key, value in sorted(pos_cell_dict.items()):
    print(key)

    scale_x = []
    for each_cell in value:
        x, residuals = nnls(RF, each_cell)
        x = np.round(x, decimals=4)
        scale_x.append(x)
    scale_x = np.array(scale_x)
    total_cells = len(scale_x)
    print("total_cells", total_cells)

    pos_pR_fp_scale_x_17[key] = scale_x

    for col in range(scale_x.shape[1]):
        if col > 0:
            RF_array = scale_x[:, col]
            RF_list = []
            for each_val in RF_array:
                if each_val > 0:
                    RF_list.append(each_val)

            RF_list = np.array(RF_list)

            bins_num = round((max(RF_list) - min(RF_list)) * 20)
            print('RF', col)

            plt.figure(figsize=(6, 4))
            hist, bins, _ = plt.hist(RF_list, bins=bins_num, color='navy', alpha=0.6, label='Reference' + str(col))
            # find the bins width
            bin_width = np.histogram_bin_edges(RF_list, bins=bins_num)[1] - np.histogram_bin_edges(RF_list,
                                                                                                   bins=bins_num)[0]

            plt.xlabel('scale_x of Reference' + str(col))

            plt.ylabel('Frequency')
            plt.title('Unmixing of positive population of pR_fp_' + key)
            plt.legend(loc='upper right')
            plt.xlim(-0.1, 10)

            plt.tight_layout()

            file_path = result_path + '10.condition3_remove_mTFP1/10_1.unmixing_without_mTFP1/'
            os.makedirs(file_path, exist_ok=True)
            filename = 'pR_fp_' + key + '_RF' + str(col) + '.png'
            plt.savefig(file_path + filename)
            plt.close()

np.save(result_path + '10_1.pos_pR_fp_scale_x_17.npy', pos_pR_fp_scale_x_17, allow_pickle=True)
np.save(result_path + '10_1.pos_cell_dict_17.npy', pos_cell_dict_17, allow_pickle=True)
np.save(result_path + '10_1.RF_peak_17.npy', RF_peak_17, allow_pickle=True)
np.save(result_path + '10_1.RF_17.npy', RF_17)
np.save(result_path + '10_1.RF_name_17.npy', RF_name_17)
