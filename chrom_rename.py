# Title: Chromosome rename
# sub-title: For renaming NC chrom names to numeric IDs as in Ensembl
# Author: Dr. Alice M. Godden

def convert_vcf(input_file, output_file):
  """
  Converts NC accession codes to numeric IDs in a VCF file.

  Args:
      input_file (str): Path to the input VCF file.
      output_file (str): Path to write the converted VCF file.
  """
  conversion_map = {
      "NC_007112.7": "1",
      "NC_007113.7": "2",
      "NC_007114.7": "3",
      "NC_007115.7": "4",
      "NC_007116.7": "5",
      "NC_007117.7": "6",
      "NC_007118.7": "7",
      "NC_007119.7": "8",
      "NC_007120.7": "9",
      "NC_007121.7": "10",
      "NC_007122.7": "11",
      "NC_007123.7": "12",
      "NC_007124.7": "13",
      "NC_007125.7": "14",
      "NC_007126.7": "15",
      "NC_007127.7": "16",
      "NC_007128.7": "17",
      "NC_007129.7": "18",
      "NC_007130.7": "19",
      "NC_007131.7": "20",
      "NC_007132.7": "21",
      "NC_007133.7": "22",
      "NC_007134.7": "23",
      "NC_007135.7": "24",
      "NC_007136.7": "25"
  }

  with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    for line in infile:
      if line.startswith('#'):  # Copy header lines directly
        outfile.write(line)
      else:
        fields = line.strip().split('\t')  # Split line by tab
        try:
          fields[0] = conversion_map[fields[0]]  # Replace ID if found
        except KeyError:
          pass  # Skip lines with non-matching IDs silently
        outfile.write('\t'.join(fields) + '\n')

if __name__ == '__main__':
  convert_vcf('C_retroseq_intersect_win150.vcf_filtered.vcf', 'C_retroseq_intersect_win150.vcf_filtered_CHROM.vcf')  # Replace with your filenames
