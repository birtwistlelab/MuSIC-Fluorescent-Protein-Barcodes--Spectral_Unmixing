import numpy as np


# to get the short name from each sample file
def get_str(s, f, b):
    par = s.partition(f)
    return(par[2].partition(b))[0][:]


# to get the scientific format of each number
def scientific_formatter(x, pos):
    return fr"$10^{{{int(np.log10(x))}}}$"


# to get the lower triangular combination list of 2 fps in a pMuSIC
def print_lower_triangular_list(lst):
    # Determine the size of the lower triangular matrix
    n = int(np.ceil(np.sqrt(2 * len(lst))))

    # Initialize a matrix with all elements as 0
    triangle_matrix = np.zeros((n, n), dtype=object)

    # Fill the lower triangle of the matrix with elements from the list
    index = 0
    for j in range(n):
        for i in range(j, n):
            if index < len(lst):
                triangle_matrix[i, j] = lst[index]
                index += 1

    # Print the lower triangular matrix
    for i in range(n):
        row = ""
        for j in range(n):
            if j <= i:
                row += f"{triangle_matrix[i, j]} "
            else:
                row += "  "  # Add spaces for upper triangle
        print(row)


# to determine the combination index among 153 possible combinations (0-152), regardless of the orientation of the
# fluorescent protein (fp)
def find_combination_index(combination, array):
    for idx, combo in array:
        if combo == combination or combo == combination[::-1]:
            return idx
    return -1

