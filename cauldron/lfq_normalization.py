import pandas as pd
import numpy as np
from scipy import linalg
from typing import List


def spread_group_label(column_index: int, data_matrix: np.ndarray, group_labels: List[int], current_group_count: int):
    """
    This function assigns the current group count as the group label for the current column and any connected columns.
    It does this by iterating over each row in the data matrix and checking if the current element is not NaN.
    If it's not NaN, it iterates over each column in the data matrix and checks if the current element is not NaN and its group label has not been assigned.
    If these conditions are met, it recursively calls itself for the current column.
    """
    # Assign the current group count as the group label for the current column
    group_labels[column_index] = current_group_count

    # Iterate over each row in the data matrix
    for row_index, row in enumerate(data_matrix):
        # If the current element is not NaN
        if not np.isnan(row[column_index]):
            # Iterate over each column in the data matrix
            for inner_column_index, element in enumerate(row):
                # If the current element is not NaN and its group label has not been assigned
                if not np.isnan(element) and np.isnan(group_labels[inner_column_index]):
                    # Recursively call the spread_group_label function for the current column
                    spread_group_label(inner_column_index, data_matrix, group_labels, current_group_count)


def compute_maxLFQ_normalization(data_subset: np.ndarray):
    """
    This function computes the MaxLFQ normalization for a subset of the data matrix.

    Parameters:
    data_subset (np.ndarray): A 2D numpy array where each row represents a different data point and each column represents a different feature or variable.

    Returns:
    np.ndarray: A 1D numpy array containing the MaxLFQ normalization for each column in the input data_subset.
    """
    # Get the number of columns in the data subset
    num_columns = data_subset.shape[1]

    # Initialize matrices for the least squares computation
    AtA = np.zeros((num_columns, num_columns))  # product of transpose and matrix
    Atb = np.zeros((num_columns, 1))  # product of transpose and vector

    """
    For each pair of columns in the data subset, compute the median of the differences between the two columns. This median difference is used to update the AtA and Atb matrices, which are used in the least squares computation.  
    The AtA matrix represents the product of the transpose of matrix A and matrix A itself, while the Atb vector represents the product of the transpose of matrix A and vector b.  
    A matrix equation is then set up with 2*AtA on the left and a column of ones on the right, and a row of ones with a zero at the end is stacked below this matrix. This forms the stacked_matrix.  
    The mean of the data subset is computed, multiplied by the number of columns, and stacked on top of 2*Atb to form the stacked_vector.  
    The lstsq function is then used to solve the least squares problem defined by the stacked_matrix and stacked_vector. The solution to this problem gives the MaxLFQ normalization for each column in the input data subset.
    """

    # Iterate over each pair of columns in the data subset
    for column_index_1, column_1 in enumerate(data_subset.T[:-1]):
        for column_index_2, column_2 in enumerate(data_subset.T[column_index_1 + 1:], start=column_index_1 + 1):
            # Compute the median of the differences between the two columns
            median_difference = np.nanmedian(-column_1 + column_2)

            # If the median difference is not NaN
            if not np.isnan(median_difference):
                # Update the AtA and Atb matrices for the least squares computation
                AtA[column_index_1, column_index_2] = AtA[column_index_2, column_index_1] = -1
                # incrementing the diagonal elements of the AtA matrix. This is done for each pair of columns in the data subset where the median difference between the two columns is not NaN.
                AtA[column_index_1, column_index_1] += 1
                AtA[column_index_2, column_index_2] += 1
                Atb[column_index_1] -= median_difference
                Atb[column_index_2] += median_difference

    # Create a matrix with 2*AtA on the left and a column of ones on the right
    left_matrix = np.hstack([2 * AtA, np.ones((num_columns, 1))])

    # Create a row of ones with a zero at the end
    right_row = np.hstack([np.ones((1, num_columns)), np.zeros((1, 1))])

    # Stack the left_matrix on top of the right_row
    stacked_matrix = np.vstack([left_matrix, right_row])

    # Compute the mean of the data_subset, multiply it by the number of columns, and stack it on top of 2*Atb
    mean_data_subset = np.nanmean(data_subset)
    stacked_vector = np.vstack([2 * Atb, mean_data_subset * num_columns])

    # Solve the least squares problem
    least_squares_solution = linalg.lstsq(stacked_matrix, stacked_vector, lapack_driver='gelsy')

    # Return the solution up to the number of columns
    return least_squares_solution[0][:num_columns]


def maxLFQ(data_matrix: np.ndarray):
    """
    This function computes the MaxLFQ normalization for a given data matrix.

    Parameters:
    data_matrix (np.ndarray): A 2D numpy array where each row represents a different data point and each column represents a different feature or variable. The array should not contain any NaN values.

    Returns:
    dict: A dictionary containing the MaxLFQ normalization for each column in the input data_matrix and the group labels as a string.
    """
    # Check if all elements in the data matrix are NaN
    if np.all(np.isnan(data_matrix)):
        return {'estimate': np.nan, 'annotation': 'NA'}

    # Check if the data matrix has only one row
    if data_matrix.shape[0] == 1:
        return {'estimate': data_matrix[0, :].astype(float), 'annotation': ''}

    # Get the number of columns in the data matrix
    num_columns = data_matrix.shape[1]

    # Initialize the current group count and a list to store the group labels for each column
    current_group_count = 0
    group_labels = [np.nan] * num_columns

    # Iterate over each column in the data matrix
    for column_index, group_label in enumerate(group_labels):
        # If the group label for the current column has not been assigned
        if np.isnan(group_label):
            # Increment the current group count
            current_group_count += 1
            # Call the spread_group_label function to assign the current group count as the group label for the current column and any connected columns
            spread_group_label(column_index, data_matrix, group_labels, current_group_count)

    # Initialize a list to store the MaxLFQ normalization for each column
    maxLFQ_normalization = [np.nan] * num_columns

    # Iterate over each group
    for group_index in range(current_group_count):
        # Get the indices of the columns that belong to the current group
        group_indices = [index for index, label in enumerate(group_labels) if label == group_index + 1]

        # If the current group has only one column
        if len(group_indices) == 1:
            # Compute the median of the column and assign it as the MaxLFQ normalization for the column
            maxLFQ_normalization[group_indices[0]] = np.nanmedian(data_matrix[:, group_indices])
        else:
            # Call the compute_maxLFQ_normalization function to compute the MaxLFQ normalization for the columns in the current group
            for index in group_indices:
                maxLFQ_normalization[index] = compute_maxLFQ_normalization(data_matrix[:, group_indices])

    # If all elements in the MaxLFQ normalization are NaN
    if np.all(np.isnan(maxLFQ_normalization)):
        return {'estimate': maxLFQ_normalization, 'annotation': 'NA'}
    else:
        # Get the indices of the columns that have been quantified (i.e., their MaxLFQ normalization is not NaN)
        quantified_columns = [index for index, value in enumerate(maxLFQ_normalization) if not np.isnan(value)]

        # If all quantified columns belong to the same group
        if all(group_labels[i] == group_labels[quantified_columns[0]] for i in quantified_columns):
            return {'estimate': maxLFQ_normalization, 'annotation': ''}
        else:
            # Assign NaN as the group label for the columns that have not been quantified
            group_labels = [np.nan if np.isnan(value) else value for value in maxLFQ_normalization]

            # Return the MaxLFQ normalization for each column and the group labels as a string
            return {'estimate': maxLFQ_normalization, 'annotation': ';'.join(map(str, group_labels))}


def process(input_filename: str, output_filename: str, id_column: str, quant_columns: List[str], data_in_log_space=False):
    """
    This function reads the input file into a pandas DataFrame, computes the MaxLFQ normalization for each group of rows that share the same id_column value, and writes the updated DataFrame to the output file.
    """
    # Read the input file into a pandas DataFrame
    df = pd.read_csv(input_filename, sep="\t")

    # Remove rows where the id_column is NaN or empty
    df = df[df[id_column].notna()]
    df = df[df[id_column] != ""]

    # Iterate over each group of rows that share the same id_column value
    for group_id, group_df in df.groupby(id_column):
        # Select the quant_columns for the current group
        group_specific_df = group_df[quant_columns]

        # If the data is not in log space, apply a log2 transformation
        group_specific_df = np.log2(group_specific_df) if not data_in_log_space else group_specific_df

        # Compute the MaxLFQ normalization for the current group
        maxLFQ_result = maxLFQ(group_specific_df)

        # Update the quant_columns of the current group in the original DataFrame with the computed MaxLFQ normalization
        df.loc[df[id_column] == group_id, quant_columns] = maxLFQ_result['estimate']

    # Write the updated DataFrame to the output file
    df.to_csv(output_filename, sep="\t", index=False)
