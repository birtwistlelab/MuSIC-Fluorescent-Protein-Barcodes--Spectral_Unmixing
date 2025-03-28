# ================================================= XL Clemson =================================================
# To randomly allocate 70% of the transfected cells from each pR-fp into the training group, while assigning the
# remaining 30% to the testing pool. Repeat 3 times.
# ================================================= XL Clemson =================================================

import os
import numpy as np
from sklearn.model_selection import train_test_split
import csv

current_path = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_path)
result_path = os.path.join(project_root, 'output/')

# Create a list of paths for 3 separate test directories to store training and testing data
results_path = [os.path.join(result_path, '4_2.SeparateTrainingAndTestingData/'f'test{i + 1}/') for i in range(3)]
for path in results_path:
    os.makedirs(path, exist_ok=True)

# Create paths for storing training and testing CSVs in each test directory
training_result_paths = [os.path.join(path, 'training_csv/') for path in results_path]
testing_result_paths = [os.path.join(path, 'testing_csv/') for path in results_path]

for path in training_result_paths:
    os.makedirs(path, exist_ok=True)

for path in testing_result_paths:
    os.makedirs(path, exist_ok=True)

# Load the dictionary containing transfected cell data
# Each key corresponds to a pR-fp name, and the value is an array of transfected cells (Nx48)
data = np.load(result_path + '4_1.transfected_cell_dict_pR_fp_singlets.npy', allow_pickle=True).item()

# Define random states to ensure reproducibility for the 3 splits
random_states = [42, 47, 888]
for i in range(3):
    training_dict = {}
    testing_dict = {}

    for key, val in sorted(data.items()):
        # Get the total number of transfected cells for the current pR-fp
        total_cells = val.shape[0]

        # Split the data: 70% training, 30% testing, using the current random state
        training, testing = train_test_split(val, test_size=0.3, random_state=random_states[i])

        print('training group', training.shape)
        print('testing group', testing.shape)

        training_dict.update({key: training})
        testing_dict.update({key: testing})

        # Create a list to store all cells assigned to the testing pool
        testing_pool = []
        for key, test_data in testing_dict.items():
            # Add all testing data from each pR-fp to the testing pool
            testing_pool.extend(test_data)
        print('cells in testing_pool', len(testing_pool))

        np.save(os.path.join(results_path[i], 'training_dict.npy'), training_dict, allow_pickle=True)
        np.save(os.path.join(results_path[i], 'testing_dict.npy'), testing_dict, allow_pickle=True)
        np.save(os.path.join(results_path[i], 'testing_pool.npy'), testing_pool)

        channel = ['V1-A', 'V2-A', 'V3-A', 'V4-A', 'V5-A', 'V6-A', 'V7-A', 'V8-A', 'V9-A', 'V10-A', 'V11-A', 'V12-A',
                   'V13-A', 'V14-A', 'V15-A', 'V16-A', 'B1-A', 'B2-A', 'B3-A', 'B4-A', 'B5-A', 'B6-A', 'B7-A', 'B8-A',
                   'B9-A', 'B10-A', 'B11-A', 'B12-A', 'B13-A', 'B14-A', 'YG1-A', 'YG2-A', 'YG3-A', 'YG4-A', 'YG5-A',
                   'YG6-A', 'YG7-A', 'YG8-A', 'YG9-A', 'YG10-A', 'R1-A', 'R2-A', 'R3-A', 'R4-A', 'R5-A', 'R6-A', 'R7-A',
                   'R8-A']

        # Save training data as CSV files
        for key, val in training_dict.items():
            training_csv = os.path.join(training_result_paths[i], f'{key}.csv')

            with open(training_csv, mode='w', newline='') as file:
                writer = csv.writer(file)
                writer.writerow(channel)
                writer.writerows(val)

        # Save the testing pool data as a CSV file
        testing_pool_csv = testing_result_paths[i] + 'pR_fp_singlets_testing_pool.csv'
        with open(testing_pool_csv, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(channel)
            writer.writerows(testing_pool)
