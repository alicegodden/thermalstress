# Title- average.py
# Author - Dr. Alice Godden
# put rows in numerical order
# Read the content of merged_output.txt
with open('merged_output.txt', 'r') as merged_file:
    lines = merged_file.readlines()

# Filter lines starting with a number
number_lines = [line.strip() for line in lines if line.strip().startswith(tuple(map(str, range(10))))]

# Sort the lines numerically
sorted_lines = sorted(number_lines, key=lambda x: int(x.split()[0]))

# Write the sorted lines to order_output.txt
with open('order_output.txt', 'w') as output_file:
    output_file.write('\n'.join(sorted_lines))


# calculate averages

from collections import defaultdict

# Read the content of order_output.txt
with open('order_output.txt', 'r') as order_file:
    lines = order_file.readlines()

# Group lines by their starting number
grouped_lines = defaultdict(list)
for line in lines:
    num, value = line.strip().split()
    grouped_lines[num].append(float(value))

# Calculate the average for each group
averages = {num: sum(values) / len(values) for num, values in grouped_lines.items()}

# Write the averages to a new file
with open('averages_output.txt', 'w') as output_file:
    for num, avg in averages.items():
        output_file.write(f"{num} {avg:.2f}\n")
