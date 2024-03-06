import matplotlib.pyplot as plt

# Load data from file
with open('tstv_Te', 'r') as file:
    lines = file.readlines()

# Extract chromosome and data values
chromosomes = []
control_values = []
temperature_values = []

for line in lines[1:]:  # Skip the header line
    parts = line.strip().split('\t')
    chromosomes.append(int(parts[0]))
    control_values.append(float(parts[1]))
    temperature_values.append(float(parts[2]))

# Create the line plot
plt.figure(figsize=(10, 6))
plt.plot(chromosomes, control_values, marker='o', label='Control', color='mediumblue', alpha=0.6)
plt.plot(chromosomes, temperature_values, marker='D', label='Temperature', color='indianred', alpha=0.6)
plt.xlabel('Chromosome', fontsize=14, weight='bold')
plt.ylabel('Ts/Tv ratio', fontsize=14, weight='bold')
plt.title('Ts/Tv ratios Testes:Control vs Temperature', fontsize=16, weight='bold')
plt.xticks(range(1, 26))
plt.yticks(weight='bold')
plt.xticks(weight='bold')
plt.grid(True, linestyle='-', linewidth=0.5, alpha=0.5)

plt.legend()
plt.grid(True)


# save the file
output_file = 'male_tstvratio.png'  # Replace with your desired filename and extension
plt.savefig(output_file, dpi=600)  # Set dpi to 300 for high resolution
plt.show()


