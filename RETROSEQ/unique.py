# Title: Removing non-unique rows from experimental retroseq vcf files
# Subtitle: Based on first 5 columns only
# Author: Dr. Alice M. Godden


import sys

def read_vcf(filename):
    with open(filename, 'r') as file:
        headers = []
        data = []
        for line in file:
            if line.startswith("#"):
                headers.append(line)
            else:
                columns = line.strip().split('\t')
                # Extract the first five columns: CHROM, POS, ID, REF, ALT
                key = '\t'.join(columns[:8])
                # Store the key and the full original line
                data.append((key, line.strip()))
        return headers, data

def check_uniqueness(data):
    unique_rows = {}
    duplicate_rows = []

    for key, full_line in data:
        if key in unique_rows:
            duplicate_rows.append(full_line)
        else:
            unique_rows[key] = full_line

    if not duplicate_rows:
        print("All rows in the experimental data are unique.")
    else:
	print(f"Warning: There are {len(duplicate_rows)} duplicate rows in the experimental data.")

    return list(unique_rows.items())  # Return key-full_line pairs
    
def compare_vcf(exp_data, ctrl_data):
    exp_dict = {key: full_line for key, full_line in exp_data}
    ctrl_set = {key for key, _ in ctrl_data}

    unique_to_exp = [full_line for key, full_line in exp_dict.items() if key not in ctrl_set]
    return unique_to_exp

def write_vcf(headers, unique_data, output_file):
    with open(output_file, 'w') as file:
        for header in headers:
            file.write(header)
        for line in unique_data:
            file.write(line + "\n")

def main():
    exp_file = "samplesFemaleTemperature_nomdups_call.out.vcf"
    ctrl_set = {key for key, _ in ctrl_data}

    unique_to_exp = [full_line for key, full_line in exp_dict.items() if key not in ctrl_set]
    return unique_to_exp

def write_vcf(headers, unique_data, output_file):
    with open(output_file, 'w') as file:
        for header in headers:
            file.write(header)
        for line in unique_data:
            file.write(line + "\n")

def main():
    exp_file = "samplesFemaleTemperature_nomdups_call.out.vcf"
    ctrl_file = "samplesFemaleControl_nomdups_call.out.vcf"
    output_file = "exp_unique_8_female_fullrow.vcf"

    print("Reading experimental VCF file...")
    exp_headers, exp_data = read_vcf(exp_file)
    print(f"Experimental file contains {len(exp_data)} data rows.")

    print("Checking for duplicate rows in the experimental data...")
    exp_data_unique = check_uniqueness(exp_data)

    print("Reading control VCF file...")
    _, ctrl_data = read_vcf(ctrl_file)
    print(f"Control file contains {len(ctrl_data)} data rows.")

    print("Comparing the files to find unique rows in the experimental file...")
    unique_to_exp = compare_vcf(exp_data_unique, ctrl_data)
    
    if unique_to_exp:
        print(f"Found {len(unique_to_exp)} unique rows in the experimental file.")
        # Print unique rows
        print("\nUnique rows in the experimental file:")
        for line in unique_to_exp:
            print(line)
    else:
	print("No unique rows found in the experimental file.")

    print(f"Writing unique rows to {output_file}...")
    write_vcf(exp_headers, unique_to_exp, output_file)
    print("Finished writing to file.")

if __name__ == "__main__":
    main()




