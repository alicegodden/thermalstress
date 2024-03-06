# Title- filter.py
# Author - Dr. Alice Godden
# Import packages
def filter_vcf(input_vcf, output_vcf):
    with open(input_vcf, 'r') as infile, open(output_vcf, 'w') as outfile:
        header_lines = []
        data_lines = []

        # Read header lines and data lines
        for line in infile:
            if line.startswith('#'):
                header_lines.append(line)
            else:
                data_lines.append(line)

        # Write header lines to the output file
        for header_line in header_lines:
            outfile.write(header_line)

        # Write filtered data rows to the output file
        for data_line in data_lines:
            if 'HIGH' in data_line.lower():
                outfile.write(data_line)

    print(f"Filtered data has been written to {output_vcf}")


# Replace these with your actual input and output filenames
input_vcf = 'deletionsonly_snpeff_male_allchrs.eff.vcf'
output_vcf = 'highmod_male_del_allchrs.eff.vcf'

filter_vcf(input_vcf, output_vcf)
