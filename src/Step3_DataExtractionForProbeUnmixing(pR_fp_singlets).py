import numpy as np
import pandas as pd
import os
import glob
from my_module import get_str

project_root = os.path.dirname(os.path.abspath(__file__))
result_path = os.path.join(project_root, 'output/')

path = os.path.join(project_root, 'csv_files/pR_fp_singlets/')
files = glob.glob(os.path.join(path, '*.csv'))
new_order = ['V1-A', 'V2-A', 'V3-A', 'V4-A', 'V5-A', 'V6-A', 'V7-A', 'V8-A', 'V9-A', 'V10-A', 'V11-A',
             'V12-A', 'V13-A', 'V14-A', 'V15-A', 'V16-A', 'B1-A', 'B2-A', 'B3-A', 'B4-A', 'B5-A', 'B6-A',
             'B7-A', 'B8-A', 'B9-A', 'B10-A', 'B11-A', 'B12-A', 'B13-A', 'B14-A', 'YG1-A', 'YG2-A',
             'YG3-A', 'YG4-A', 'YG5-A', 'YG6-A', 'YG7-A', 'YG8-A', 'YG9-A', 'YG10-A', 'R1-A', 'R2-A',
             'R3-A', 'R4-A', 'R5-A', 'R6-A', 'R7-A', 'R8-A']

pR_fp_singlets = {}
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
        for i in range(len(updated_col)):
            new_row_list.append(row[updated_col[i]])
        new_file.append(new_row_list)

    front_str = 'csv_files/pR_fp_singlets/'
    back_str = '_singlets.csv'
    name = get_str(file, front_str, back_str)

    singlets = np.array(new_file)
    pR_fp_singlets.update({name: singlets})

print(pR_fp_singlets)
print(sorted(pR_fp_singlets))

np.save(result_path + '3.pR_fp_singlets.npy', pR_fp_singlets, allow_pickle=True)
