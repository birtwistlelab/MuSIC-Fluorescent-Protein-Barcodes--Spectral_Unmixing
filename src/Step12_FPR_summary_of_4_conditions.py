import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from my_module import scientific_formatter

plt.rcParams['font.family'] = 'DejaVu Sans'

current_path = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_path)
result_path = os.path.join(project_root, 'output/')
fig_path = os.path.join(project_root, 'paper_fig/')

condition1_fpr = np.load(result_path + '8_1.condition1_original_FPR.npy')
condition1_fpr = condition1_fpr.tolist()
condition2_fpr = np.load(result_path + '9.condition2_Intensity_cutoff_FPR.npy')
condition2_fpr = condition2_fpr.tolist()
condition3_fpr = np.load(result_path + '10_3.condition3_remove_mTFP1_FPR.npy')
condition3_fpr = condition3_fpr.tolist()
condition4_fpr = np.load(result_path + '11.condition4_remove_and_cutoff_FPR.npy')
condition4_fpr = condition4_fpr.tolist()

fp_name = ["EBFP2", "mTagBFP2", "mT-Sapphire", "mAmetrine", "mCerulean3", "LSSmOrange", "mBeRFP", "mTFP1", "EGFP",
           "CyOFP1", "mClover3", "mVenus", "mPapaya", "mOrange2", "mRuby3", "mKate2", "mCardinal", "miRFP670"]

# to make the FPR summary bar plot

# to make sure the position of 08.mTFP1 was N/A
condition3_fpr.insert(7, np.nan)
condition4_fpr.insert(7, np.nan)

# to make sure the 4 conditions have the same length
assert len(condition1_fpr) == len(condition2_fpr) == len(condition3_fpr) == len(condition4_fpr)

data = [condition1_fpr, condition2_fpr, condition3_fpr, condition4_fpr]
fig, ax = plt.subplots(figsize=(8, 4))
bar_width = 0.5
group_gap = 2
index = np.arange(len(fp_name)) * (4 * bar_width + group_gap)

plt.bar(index, condition1_fpr, bar_width, label='Original', color='navy', alpha=0.6)
plt.bar(index + bar_width, condition2_fpr, bar_width, label='Intensity Cutoff', color='orange', alpha=0.6)
plt.bar(index + 2 * bar_width, condition3_fpr, bar_width, label='Remove mTFP1', color='gray', alpha=0.6)
plt.bar(index + 3 * bar_width, condition4_fpr, bar_width, label='Intensity Cutoff and Remove mTFP1', color='purple',alpha=0.6)

plt.ylabel('False Possitive Rate', fontsize=16)

ax.set_xticks(index + 1.5 * bar_width)
ax.set_xticklabels(fp_name, rotation=45, ha='right', fontsize=12)
ax.tick_params(axis='y', labelsize=12)
ax.set_ylim([0, 0.3])

plt.legend()
plt.subplots_adjust(bottom=0.3)

file_path = result_path + '12.FPR_Summary_of_4_conditions/'
os.makedirs(file_path, exist_ok=True)
filename = 'Summary of FPR on 4 condtions.png'
plt.savefig(file_path + filename, transparent=True)
plt.close()

# to make the log y-axis of FPR summary plots for fig.3E
data = [condition1_fpr, condition2_fpr, condition3_fpr, condition4_fpr]
colors = ['navy', 'orange', 'gray', 'purple']
conditions = ['Original', 'Intensity Cutoff', 'Remove mTFP1', 'Intensity Cutoff and Remove mTFP1']

fig, ax = plt.subplots(figsize=(8, 4))
bar_width = 0.5
group_gap = 2
index = np.arange(len(fp_name)) * (4 * bar_width + group_gap)

for i, fpr in enumerate(data):
    ax.bar(index + i * bar_width, fpr, bar_width, label=conditions[i], color=colors[i], alpha=0.6)

ax.set_ylabel('False Possitive Rate', fontsize=12)
ax.set_yscale('log')
ax.set_ylim(10 ** -4, 10 ** 0)
plt.gca().yaxis.set_major_formatter(FuncFormatter(scientific_formatter))

ax.set_xticks(index + 1.5 * bar_width)
ax.set_xticklabels(fp_name, rotation=45, ha='right', fontsize=10)
ax.tick_params(axis='y', labelsize=10)
ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1), ncol=4, fontsize=8, frameon=False)
plt.subplots_adjust(bottom=0.22)

file_path2 = result_path + '12.FPR_Summary_of_4_conditions/log_y_axis_for_Fig3E/'
os.makedirs(file_path2, exist_ok=True)
filename3 = 'Summary of log_FPR on 4 condtions.png'
plt.savefig(file_path2 + filename3, transparent=True)

filename4 = 'Fig3E(Step12).png'
plt.savefig(fig_path+ filename4, transparent=True)
plt.close()
