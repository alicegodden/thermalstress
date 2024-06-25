# Title- filter_retroseq.py
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
            fields = data_line.strip().split('\t')
            format_field = fields[8]
            format_values = fields[9].split(':')

            gq_index = format_field.split(':').index('GQ')
            fl_index = format_field.split(':').index('FL')

            gq_value = float(format_values[gq_index])
            fl_value = int(format_values[fl_index])

            if gq_value >= 1000 and fl_value >= 8:
                outfile.write(data_line)

    print(f"Filtered data has been written to {output_vcf}")

# Replace these with your actual input and output filenames
input_vcf = 'FILTERmerged_with_header_USEME_win150_int_filteredHEADER.vcf'
output_vcf = 'HEADERFILTERPY_merged_with_header_USEME_win150_int_filtered.vcf'

filter_vcf(input_vcf, output_vcf)
