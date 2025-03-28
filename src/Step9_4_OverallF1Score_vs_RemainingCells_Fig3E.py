import numpy as np
import os
import pandas as pd
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline

plt.rcParams['font.family'] = 'Arial'

current_path = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_path)
result_path = os.path.join(project_root, 'output/')

# Load the autofluorescence reference (MFI of pR cells) data
RF = np.load(result_path + '2.MFI_RF.npy')
print(f'ref_autofluorescence (MFI of pR cells): \n {RF[0]}')
print(f'highest FI of autofluorescence: {RF[0].max()}')

# Load summary data paths from multiple test folders
loading_summary_path = [os.path.join(result_path, '9_2.CalculationAtEachCutoff/', f'test{i + 1}/')
                        for i in range(3)]
cutoff_thresholds = np.load(os.path.join(result_path, '9_1.cutoff_thresholds.npy'))
print(cutoff_thresholds)
saving_path = os.path.join(result_path, '9_4.OverallF1Score_vs_RemainingCells_Fig3E/')
os.makedirs(saving_path, exist_ok=True)

# Load summary dictionaries from all three tests
summary_mcc_dicts = {}
summary_cell_number_dicts = {}
for i in range(3):
    # Load the F1 score summary and cell number summary for each test
    summary_dict = np.load(os.path.join(loading_summary_path[i], 'overall_f1_dict.npy'), allow_pickle=True).item()
    cell_number_dict = np.load(os.path.join(loading_summary_path[i], 'total_cells_dict1.npy'), allow_pickle=True).item()

    # Store the MCC values for each test in a list
    mcc_list = []
    for key, val in summary_dict.items():
        mcc_list.append(val)
    summary_mcc_dicts[f'test{i + 1}'] = mcc_list

    # Store the cell numbers for each cutoff in a list
    cell_number_list =[]
    for key2, val2 in cell_number_dict.items():
        cell_number_list.append(val2)
    summary_cell_number_dicts[f'test{i + 1}'] = cell_number_list
print(summary_mcc_dicts['test1'])

# Convert the summary data into DataFrames for easier manipulation and saving
df = pd.DataFrame(summary_mcc_dicts, columns=['test1', 'test2', 'test3'])
df2 = pd.DataFrame(summary_cell_number_dicts, columns=['test1', 'test2', 'test3'])

# Add 'Intensity_Cutoff' column to df and calculate the average and standard deviation for F1 score
df.insert(0, 'Intensity_Cutoff', cutoff_thresholds)
df['Average_F1'] = df[['test1', 'test2', 'test3']].mean(axis=1).round(3)
df['Std_F1'] = df[['test1', 'test2', 'test3']].std(axis=1).round(3)

df.to_excel(saving_path + 'summary of Overall F1 Score of tests.xlsx', index=False)

# Calculate the average and standard deviation for cell numbers across the tests
df2['Average_cell_number'] = df2[['test1', 'test2', 'test3']].mean(axis=1).round(3)
df2['Std_cell_number'] = df2[['test1', 'test2', 'test3']].std(axis=1).round(3)
df2.to_excel(saving_path + 'summary of cell number of tests.xlsx', index=False)


# Function to format numbers in scientific notation for plot labels
def scientific_formatter(x, pos):
    exponent = int(np.log10(x))
    mantissa = x / (10 ** exponent)
    if mantissa > 1:
        return fr"${mantissa:.0f} \times 10^{exponent}$"
    if mantissa == 1:
        return fr"$10^{{{int(np.log10(x))}}}$"

# Extract the average cell numbers and F1 scores
x = np.array(df2['Average_cell_number'])
y = np.array(df['Average_F1'])

# Calculate the standardized error bars for the x and y values
std_x = np.array(df2['Std_cell_number'])/max(x)
std_y = np.array(df['Std_F1'])

# Normalize x values
x_fraction = x/max(x)
print(x_fraction)
print(y)

# Identify the optimal cutoff (with a signal-to-noise ratio of 10)
best_cutoff = 10000
optimal_idx = list(cutoff_thresholds).index(10000)
optimal_x = x_fraction[optimal_idx]
optimal_y = y[optimal_idx]

# Ensure that the x values are sorted for smooth plotting
sorted_indices = np.argsort(x_fraction)
x_fraction_sorted = x_fraction[sorted_indices]
y_sorted = y[sorted_indices]

# Smooth the curve using cubic spline interpolation
x_smooth = np.linspace(x_fraction_sorted.min(), x_fraction_sorted.max(), 525)

spl = make_interp_spline(x_fraction_sorted, y_sorted, k=3)  # k=3 for cubic splines
y_smooth = spl(x_smooth)

# Interpolate cutoff values to match the 500 smooth points
cutoffs_sorted = np.array(cutoff_thresholds)[sorted_indices]

# Since the highest MFI of autofluorescence is around 1000, setting the cutoff at 1000 ensures an SNR of 1.
# We could normalize all cutoffs to 1000 accordingly.
cutoffs_sorted_normalized = cutoffs_sorted / 1000
cutoff_interpolator = interp1d(x_fraction_sorted, cutoffs_sorted_normalized, kind='linear', fill_value="extrapolate")
cutoffs_smooth = cutoff_interpolator(x_smooth)

# Interpolate cutoff values to match the 525 smooth points
norm = plt.Normalize(cutoffs_sorted_normalized.min(), cutoffs_sorted_normalized.max())
cmap = plt.get_cmap('viridis')

# ====================================== XL Clemson ======================================
# Plot the smoothed curve with color representing the SNR
# ====================================== XL Clemson ======================================

plt.figure(figsize=(6, 4))
for i in range(len(x_smooth) - 1):
    color = cmap(norm(cutoffs_smooth[i]))
    plt.plot(x_smooth[i:i+2], y_smooth[i:i+2], color=color, lw=2)

plt.gca().invert_xaxis()
plt.scatter(optimal_x, optimal_y, color='orange', label=f'Intensity cutoff at {scientific_formatter(best_cutoff, None)} \n SNR = 10', zorder=5, s=50)
offset = 0.02
offset2 = 0.001
plt.text(optimal_x - offset, optimal_y - offset2, f'({optimal_x:.3f}, {optimal_y:.3f})', fontsize=12, verticalalignment='top', horizontalalignment='left')

plt.tick_params(labelsize=12)
plt.xlabel("Fraction of Cells Remaining ", fontsize=14)
plt.ylabel("Overall F1 Score", fontsize=14)
plt.legend(loc='lower right', fontsize=14, framealpha=0, frameon=False,  handlelength=1.5, handleheight=1, handletextpad=0.3, labelspacing=0.2)
plt.grid(False)

# Add color bar
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])  # Empty array to create color bar
cbar_ax = plt.gcf().add_axes([0.85, 0.15, 0.03, 0.7])
cbar = plt.colorbar(sm, cax=cbar_ax)
cbar.set_label("Signal-to-Noise Ratio (SNR)", fontsize=14)
# cbar.set_ticks([1000, 10000, 20000, 30000, 40000, 50000])
plt.subplots_adjust(bottom=0.15, right=0.8, left=0.15, top=0.85)
# cbar.ax.yaxis.set_major_formatter(FuncFormatter(scientific_formatter))

filename = 'OverallF1Score.png'
plt.savefig(os.path.join(saving_path, filename), dpi=300, transparent=True)

filename2 = 'Fig.3E(Step9_4)_OverallF1ScorevsIntensityCutoff.png'
plt.savefig(os.path.join(project_root, 'figures_for_paper/', filename2), dpi=300, transparent=True)
