# to extract the data from the csv files of positive pR-fps by manually gating to create the reference for
# further unmixing
import numpy as np
import pandas as pd
import os
import glob
from my_module import get_str

project_root = os.path.dirname(os.path.abspath(__file__))
print(project_root)

result_path = os.path.join(project_root, 'output/')
os.makedirs(result_path, exist_ok=True)

fig_path = os.path.join(project_root, 'paper_fig/')
os.makedirs(fig_path, exist_ok=True)

path = os.path.join(project_root, 'csv_files/pos_pR_fp/')
print(path)

files = glob.glob(os.path.join(path, '*.csv'))
print(len(files))
new_order = ['V1-A', 'V2-A', 'V3-A', 'V4-A', 'V5-A', 'V6-A', 'V7-A', 'V8-A', 'V9-A', 'V10-A', 'V11-A',
             'V12-A', 'V13-A', 'V14-A', 'V15-A', 'V16-A', 'B1-A', 'B2-A', 'B3-A', 'B4-A', 'B5-A', 'B6-A',
             'B7-A', 'B8-A', 'B9-A', 'B10-A', 'B11-A', 'B12-A', 'B13-A', 'B14-A', 'YG1-A', 'YG2-A',
             'YG3-A', 'YG4-A', 'YG5-A', 'YG6-A', 'YG7-A', 'YG8-A', 'YG9-A', 'YG10-A', 'R1-A', 'R2-A',
             'R3-A', 'R4-A', 'R5-A', 'R6-A', 'R7-A', 'R8-A']

fp_reference = {}
for file in files:
    data = pd.read_csv(file)
    col_names = data.columns.tolist()

    updated_col = []
    for i in range(len(new_order)):
        if new_order[i] in col_names:
            col_index = col_names.index(new_order[i])
            updated_col.append(col_index)

    new_file = []
    for each_row in data.values:
        row = each_row.tolist()
        new_row_list = []
        for j in range(len(updated_col)):
            new_row_list.append(row[updated_col[j]])
        new_file.append(new_row_list)

    front_str = 'csv_files/pos_pR_fp/'
    back_str = '.csv'
    name = get_str(file, front_str, back_str)

    pos_reference = np.array(new_file)
    fp_reference.update({name: pos_reference})

np.set_printoptions(suppress=True)
print(sorted(fp_reference))

print(len(fp_reference['00.unstained']))

np.save(result_path + '1.fp_reference_2023.npy', fp_reference, allow_pickle=True)
np.save(result_path + '1.x_axis_channel.npy', new_order)


