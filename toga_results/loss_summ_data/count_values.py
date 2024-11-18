import sys
from collections import defaultdict

file_name = sys.argv[1]

possible_values = ['I', 'L', 'UL', 'PI', 'M', 'PG', 'PM']

try:
    column_counts = defaultdict(lambda: defaultdict(int))

    # Open and read the file
    with open(file_name, 'r') as f:
        header = f.readline().strip().split('\t')
        num_columns = len(header)

        # Initialize counts for each value in each column
        for col_idx in range(num_columns):
            for value in possible_values:
                column_counts[col_idx][value] = 0

        # Count occurrences of each value in each column
        for line in f:
            columns = line.strip().split('\t')
            for col_idx, value in enumerate(columns):
                if value in possible_values:
                    column_counts[col_idx][value] += 1

    # Print count and percentage of each value in each column
    print("Count and percentage of each value in each column:")
    for col_idx, col_name in enumerate(header):
        total = sum(column_counts[col_idx].values())  # Total count for the column
        print(f"\n{col_name}:")
        for value in possible_values:
            count = column_counts[col_idx][value]
            percentage = (count / total * 100) if total > 0 else 0
            print(f"  {value}: {count} ({percentage:.2f}%)")

