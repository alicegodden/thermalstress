#Title: Water testing data and script
#Author: Dr. Alice M. Godden


import matplotlib.pyplot as plt

# Data for control group
control_data = {
    "Oxygen (ppm, mg/L)": [8, 8, 8],
    "KH (hardness, °dH)": [10, 11, 11],
    "pH meter": [7.4, 7.4, 7.4],
    "Temperature (°C)": [28, 28, 27.9],
    "Carbon dioxide (mg/L)": [12, 13, 13],
    "Conductivity": [398, 422, 370],
    "Light (lux)": [112, 122, 113],
    "Nitrites (mg/L)": [0, 0, 0],
}

# Data for temperature group
temperature_data = {
    "Oxygen (ppm, mg/L)": [6, 6, 5],
    "KH (hardness, °dH)": [11, 13, 13],
    "pH meter": [7.3, 7.5, 7.4],
    "Temperature (°C)": [34, 33.7, 33.4],
    "Carbon dioxide (mg/L)": [16.5, 15.5, 15.5],
    "Conductivity": [388, 391, 398],
    "Light (lux)": [118, 124, 114],
    "Nitrites (mg/L)": [0, 0, 0.05],
}

# Time points
time_points = ["Start", "1 week", "2 weeks"]

# Plotting the data
fig, axes = plt.subplots(nrows=2, ncols=4, figsize=(20, 10))
fig.tight_layout(pad=3.0)

# Iterate through all parameters in control_data (considering first 8 parameters)
for i, (parameter, values) in enumerate(list(control_data.items())[:8]):
    row = i // 4
    col = i % 4
    ax = axes[row, col]

    # Check if the parameter exists in temperature_data and if control_data has no None values
    if parameter in temperature_data and None not in values:
        ax.plot(time_points, values, label="Control", marker='o')
        ax.plot(time_points, temperature_data[parameter], label="Temperature", marker='D')
        ax.set_title(parameter)
        ax.set_xlabel("Time")
        ax.set_ylabel(parameter)
        ax.legend()
    else:
        # Handle cases where data might be missing or inconsistent
        if None in values:
            print(f"Control data for '{parameter}' contains None values and will not be plotted.")
        else:
            print(f"Temperature data for '{parameter}' is missing and will not be plotted.")

        # Ensure the plot area is cleared to avoid blank charts
        ax.axis('off')

# Hide the remaining subplot in the last row and column if there are less than 8 parameters
if len(control_data) < 8:
    for j in range(len(control_data), 8):
        axes[j // 4, j % 4].axis('off')

plt.savefig("comparison_line_graphs.png")
plt.show()
