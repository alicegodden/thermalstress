# Title- class_chart.py
# Author - Dr. Alice Godden
# Import packages

import re
import pandas as pd
import matplotlib.pyplot as plt

def extract_meinfo(line):
  """Extracts the MEINFO value from a line of VCF text."""

  pattern = r'MEINFO=(?P<meinfo>\w+)'
  match = re.search(pattern, line)
  if match:
    return match.group('meinfo')
  else:
    return None

def count_meinfo_words(meinfo_values, words):
  """Counts the number of times each word in a list of words appears in a list of MEINFO values."""

  meinfo_word_counts = {}
  for word in words:
    meinfo_word_counts[word] = 0

  for meinfo_value in meinfo_values:
    if meinfo_value is not None:
      for word in words:
        if word in meinfo_value:
          meinfo_word_counts[word] += 1

  return meinfo_word_counts

def plot_meinfo_word_counts(meinfo_word_counts):
  """Plots the counts of each MEINFO word with customizations."""

  meinfo_word_counts = pd.Series(meinfo_word_counts)
  meinfo_word_counts.plot(kind='bar')

  # Make text bold
  plt.xlabel('TE Class', fontweight='bold')  # Make x-axis label bold
  plt.ylabel('Count', fontweight='bold')  # Make y-axis label bold
  plt.title('Central', fontweight='bold')  # Make title bold

  # Rotate x-axis labels by 45 degrees
  plt.xticks(rotation=45, ha='right', fontweight = 'bold')  # Rotate x-axis labels by 45 degrees
  plt.yticks(fontweight = 'bold')
  plt.tight_layout()  # Adjust spacing for better layout




if __name__ == '__main__':
  # Read the VCF file.
  with open('C_retroseq_intersect_win150.vcf_filtered_CHROM.vcf', 'r') as f:
    vcf_lines = f.readlines()

  # Extract the MEINFO values from the VCF lines.
  meinfo_values = []
  for line in vcf_lines:
    meinfo = extract_meinfo(line)
    if meinfo is not None:
      meinfo_values.append(meinfo)

  # Count the number of times each word appears in the MEINFO values.
  words = ['LINE', 'LTR', 'DNA', 'SATELLITE', 'SINE', 'RC']
  meinfo_word_counts = count_meinfo_words(meinfo_values, words)

  # Plot the counts of each word with customizations
  plot_meinfo_word_counts(meinfo_word_counts)

  print(meinfo_word_counts)

# save the file
output_file = 'C_retroseq_intersect_win150_chrom_counts.png'  # Replace with your desired filename and extension
plt.savefig(output_file)
plt.show()
