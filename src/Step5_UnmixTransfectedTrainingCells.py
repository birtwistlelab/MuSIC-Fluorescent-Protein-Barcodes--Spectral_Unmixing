# ====================================== XL Clemson ======================================
# Now we will unmix the transfected training cells and create their 18x18 plots
# ====================================== XL Clemson ======================================

import numpy as np
from scipy.optimize import nnls
import matplotlib.pyplot as plt
import os
plt.rcParams['font.family'] = 'Arial'

current_path = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_path)
result_path = os.path.join(project_root, 'output/')
loading_paths = [os.path.join(result_path, '4_2.SeparateTrainingAndTestingData/', f'test{i + 1}/') for i in range(3)]

saving_paths = [os.path.join(result_path, '5.UnmixingTransfectedTrainingCells/', f'test{i + 1}/') for i in range(3)]
for path in saving_paths:
    os.makedirs(path, exist_ok=True)

# Load required data files：
RF_peak = np.load(os.path.join(result_path, '2.RF_peak.npy'), allow_pickle=True).item()
RF = np.load(os.path.join(result_path, '2.MFI_RF.npy'))
RF_name = np.load(os.path.join(result_path,'2.RF_name.npy'))
x_axis = np.load(os.path.join(result_path,'1.x_axis_channel.npy'))
channel_num = len(x_axis)

# Transpose RF once before looping, or it may cause errors
RF = np.array(RF).T

# Loop through three test sets：
for i in range(3):
    # Define the loading path for the current test set：
    loading_path = loading_paths[i]
    # Load the training dictionary containing transfected cells for each pR-fp：
    transfected_cell_dict = np.load(os.path.join(loading_path, 'training_dict.npy'), allow_pickle=True).item()

    # Define labels for visualization：
    row_labels = RF_name[1:]
    col_labels = ['01.EBFP2', '02.mTagBFP2', '03.mT_Sapphire', '04.mAmetrine', '05.mCerulean3', '06.LSSmOrange',
                  '07.mBeRFP', '08.mTFP1', '09.EGFP', '10.CyOFP1', '11.mClover3', '12.mVenus', '13.mPapaya',
                  '14.mOrange2', '15.mRuby3', '16.mKate2', '17.mCardinal', '18.miRFP670']
    x_label = 'Unmixed Relative Abundance'
    y_label = 'Frequency'

    # Transpose RF to align with cell fluorescence values：
    transfected_pR_fp_scale_x = {}
    for key, value in sorted(transfected_cell_dict.items()):
        scale_x = []  # List to store unmixed fluorescence intensity values

        # Apply non-negative least squares (NNLS) unmixing for each cell：
        for each_cell in value:
            x, residuals = nnls(RF, each_cell)
            x = np.round(x, decimals=4)
            scale_x.append(x)
        scale_x = np.array(scale_x)
        total_cells = len(scale_x)

        transfected_pR_fp_scale_x[key] = scale_x

        # Generate histograms for each reference fluorescence channel：
        for col in range(scale_x.shape[1]):
            if col > 0:  # Skip first column which is the background (untransfected/autofluorescence)
                RF_array = scale_x[:, col]
                RF_list = []
                for each_val in RF_array:
                    if each_val > 0:
                        RF_list.append(each_val)

                RF_list = np.array(RF_list)

                # Determine the number of histogram bins:
                if len(RF_list) > 0:
                    bins_num = round((max(RF_list) - min(RF_list)) * 20)
                else:
                    bins_num = 10  # Default to 10 bins if list is empty

                # Plot histogram of fluorescence intensity values:
                plt.figure(figsize=(6, 4))
                hist, bins, _ = plt.hist(RF_list, bins=bins_num, color='navy', alpha=0.6, label='Reference' + str(col))
                # Calculate bin width:
                bin_width = np.histogram_bin_edges(RF_list, bins=bins_num)[1] - np.histogram_bin_edges(RF_list,
                                                                                                       bins=bins_num)[0]

                plt.xlabel('scale_x of Reference' + str(col), fontsize=18)

                plt.ylabel('Frequency', fontsize=18)
                # plt.title('Unmixing of transfected population of pR_fp_' + key)
                plt.legend(loc='upper right', fontsize=18, frameon=False)
                plt.xlim(-0.1, 10)
                plt.tick_params(labelsize=14)

                plt.tight_layout()

                file_path = os.path.join(saving_paths[i], '18x18_plots', f'{key[:2]}.pR_{key[3:]}_TrainingCells')
                os.makedirs(file_path, exist_ok=True)
                filename = key + '_RF' + str(col) + '.png'
                plt.savefig(os.path.join(file_path, filename), dpi=300, transparent=True)
                plt.close()
    np.save(os.path.join(saving_paths[i], '5.transfected_Training_pR_fp_scale_x.npy'), transfected_pR_fp_scale_x, allow_pickle=True)
