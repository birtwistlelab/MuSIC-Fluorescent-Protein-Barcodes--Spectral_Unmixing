# from here, we will start to do data analysis of singly transfected pMuSICs

import numpy as np
import pandas as pd
import os
import glob
from my_module import get_str

current_path = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_path)
result_path = os.path.join(project_root, 'output/')

path = os.path.join(project_root, 'csv_files/pos_pMuSICs/')
files = glob.glob(os.path.join(path, '*.csv'))
print(len(files))

new_order = ['V1-A', 'V2-A', 'V3-A', 'V4-A', 'V5-A', 'V6-A', 'V7-A', 'V8-A', 'V9-A', 'V10-A', 'V11-A',
             'V12-A', 'V13-A', 'V14-A', 'V15-A', 'V16-A', 'B1-A', 'B2-A', 'B3-A', 'B4-A', 'B5-A', 'B6-A',
             'B7-A', 'B8-A', 'B9-A', 'B10-A', 'B11-A', 'B12-A', 'B13-A', 'B14-A', 'YG1-A', 'YG2-A',
             'YG3-A', 'YG4-A', 'YG5-A', 'YG6-A', 'YG7-A', 'YG8-A', 'YG9-A', 'YG10-A', 'R1-A', 'R2-A',
             'R3-A', 'R4-A', 'R5-A', 'R6-A', 'R7-A', 'R8-A']

pos_pMuSICs = {}

for file in files:
    # get the data from each file
    data = pd.read_csv(file)
    # now we have the data for each singlet, it's dataFrame containing index as data.columns and the fluorescence
    # intensity as the data.values
    col_names = data.columns.tolist()
    # print(data.values)

    updated_col = []  # this is the index list for col_names in new_order
    for i in range(len(new_order)):
        if new_order[i] in col_names:
            col_index = col_names.index(new_order[i])
            updated_col.append(col_index)

    new_file = []
    for each_row in data.values:
        row = each_row.tolist()
        new_row_list = []  # get each row with updated col order
        for i in range(len(updated_col)):
            new_row_list.append(row[updated_col[i]])  # get each value with updated col order in each row
        new_file.append(new_row_list)

    front_str = 'csv_files/pos_pMuSICs/export_'
    back_str = '2024 P'
    name = get_str(file, front_str, back_str)

    singlets = np.array(new_file)
    pos_pMuSICs.update({name: singlets})

print(pos_pMuSICs)
print(sorted(pos_pMuSICs))

np.save(result_path + '13.all_pos_pMuSICs.npy', pos_pMuSICs, allow_pickle=True)
