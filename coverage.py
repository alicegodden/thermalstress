# Title- coverage.py
# Author - Dr. Alice Godden
# plot the average coverage of each chromosome in a plot
# Import packages
import matplotlib.pyplot as plt


# Read data from averages_output.txt
with open('averages_output.txt', 'r') as input_file:
    lines = input_file.readlines()

# Extract x and y values
x = []
y = []
for line in lines:
    num, value = line.strip().split()
    x.append(int(num))
    y.append(float(value))

# Plot the data
plt.figure(figsize=(10, 6))
plt.plot(x, y, marker='o', linestyle='none', color='red', alpha=0.6)
plt.xlabel('Chromosome', weight='bold')
plt.ylabel('Coverage', weight='bold')
plt.title('Average Coverage per Chromosome', weight='bold')
plt.yticks(weight='bold')
plt.xticks(weight='bold')
plt.grid(linewidth=0.25)
# Set x-axis tick marks and labels
plt.xticks(range(min(x), max(x) + 1))
# set y axis range
plt.ylim(50, 60)  # Adjust the range as needed

# Save the plot as PNG
plt.savefig('coverage_plot_HS_gDNA_all_samples.png', dpi=600)
plt.savefig('coverage_plot_HS_gDNA_all_samples.pdf')

#show the plot
plt.show()
