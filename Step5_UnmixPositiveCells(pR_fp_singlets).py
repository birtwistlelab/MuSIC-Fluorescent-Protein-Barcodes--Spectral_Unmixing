import numpy as np
from scipy.optimize import nnls
import matplotlib.pyplot as plt
import os
plt.rcParams['font.family'] = 'DejaVu Sans'
project_root = os.path.dirname(os.path.abspath(__file__))
result_path = os.path.join(project_root, 'output/')

pos_cell_dict = np.load(result_path + '4.pos_cell_dict_pR_fp_singlets.npy', allow_pickle=True).item()
RF_peak = np.load(result_path + '2.RF_peak.npy', allow_pickle=True).item()
RF = np.load(result_path + '2.MFI_RF.npy')
RF_name = np.load(result_path + '2.RF_name.npy')

x_axis = np.load(result_path + '1.x_axis_channel.npy')
channel_num = len(x_axis)

row_labels = RF_name[1:]
col_labels = ["01.EBFP2", "02.mTagBFP2", "03.mT-Sapphire", "04.mAmetrine", "05.mCerulean3", "06.LSSmOrange",
              "07.mBeRFP", "08.mTFP1", "09.EGFP", "10.CyOFP1", "11.mClover3", "12.mVenus", "13.mPapaya", "14.mOrange2",
              "15.mRuby3", "16.mKate2", "17.mCardinal", "18.miRFP670"]
x_label = 'Unmixed Relative Abundance'
y_label = 'Frequency'

RF = np.array(RF).T
pos_pR_fp_scale_x = {}
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

    pos_pR_fp_scale_x[key] = scale_x

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

            file_path = result_path + '5.unmixed_pos_pR_singlets/'
            os.makedirs(file_path, exist_ok=True)
            filename = 'pR_fp_' + key + '_RF' + str(col) + '.png'
            plt.savefig(file_path + filename)
            plt.close()

np.save(result_path + '5.pos_pR_fp_scale_x.npy', pos_pR_fp_scale_x, allow_pickle=True)