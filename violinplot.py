# Import packages
import matplotlib.pyplot as plt
import seaborn as sns  # Make sure seaborn is installed (you can install it with: pip install seaborn)

# Read data from averages_output.txt
with open('order_output.txt', 'r') as input_file:
    lines = input_file.readlines()

# Extract x and y values
x = []
y = []
for line in lines:
    num, value = line.strip().split()
    x.append(int(num))
    y.append(float(value))

# Create a violin plot
plt.figure(figsize=(10, 6))
sns.violinplot(x=x, y=y, color='red')  # Specify x and y variables
plt.xlabel('Chromosome', weight='bold')
plt.ylabel('Coverage', weight='bold')
plt.title('Average Coverage per Chromosome', weight='bold')
plt.yticks(weight='bold')
plt.xticks(weight='bold')
plt.grid(linewidth=0.25)
# Set x-axis tick marks and labels
plt.xticks(range(0, 25))
# Set y-axis range
plt.ylim(0, 100)  # Adjust the range as needed

# Save the plot as PNG
plt.savefig('violinplot_coverage_HS_gDNA_all_samples.png', dpi=600)
plt.savefig('violinplot_coverage_HS_gDNA_all_samples.pdf')

# Show the plot
plt.show()
