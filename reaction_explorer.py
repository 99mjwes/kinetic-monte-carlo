import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


row_of_interest = 146

threshold = 0.00001

count_data = pd.read_csv('reaction_counts.csv', sep=',', header=0)
# Extract reaction names and counts
headers = list(count_data.head())
names = headers[:-1:2]  # Exclude the last column (ReactionRate)

counts = count_data.iloc[row_of_interest, :-1:2].values.astype(float)  # Extract counts from the specified row
total_count = np.sum(counts)
sorted_count = np.sort(counts)[::-1]  # Sort counts in descending order
sorted_names = [name for _, name in sorted(zip(counts, names), reverse=True)]

i = 0
while i < len(sorted_count) and sorted_count[i] > threshold * total_count:
    print(f"{sorted_names[i]}: {sorted_count[i]}")
    i += 1