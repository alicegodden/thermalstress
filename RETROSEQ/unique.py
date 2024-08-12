# Title: Removing non-unique rows from experimental retroseq vcf files
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
                data.append(line.strip())
        return headers, data

def check_uniqueness(data):
    unique_rows = set(data)
    if len(unique_rows) == len(data):
        print("All rows in the experimental data are unique.")
    else:
        duplicates = len(data) - len(unique_rows)
        print(f"Warning: There are {duplicates} duplicate rows in the experimental data.")
    return unique_rows

def compare_vcf(exp_data, ctrl_data):
    exp_set = set(exp_data)
    ctrl_set = set(ctrl_data)

    unique_to_exp = exp_set - ctrl_set
    return unique_to_exp

def write_vcf(headers, unique_data, output_file):
    with open(output_file, 'w') as file:
        for header in headers:
            file.write(header)
        for line in unique_data:
            file.write(line + "\n")

def main():
    exp_file = "samplesMaleTemperature_nomdups_call.out.vcf"
    ctrl_file = "samplesMaleControl_nomdups_call.out.vcf"
    output_file = "exp_unique_male.vcf"

    print("Reading experimental vcf file...")
    exp_headers, exp_data = read_vcf(exp_file)
    print(f"Experimental file contains {len(exp_data)} data rows.")

    print("Checking for duplicate rows in the experimental data...")
    exp_data_unique = check_uniqueness(exp_data)

    print("Reading control VCF file...")
    _, ctrl_data = read_vcf(ctrl_file)
    print(f"Control file contains {len(ctrl_data)} data rows.")

    print("Comparing the files to find unique rows in the experimental vcf file...")
    unique_to_exp = compare_vcf(exp_data_unique, ctrl_data)

    if unique_to_exp:
        print(f"Found {len(unique_to_exp)} unique rows in the experimental vcf file.")
    else:
        print("No unique rows were found in the experimental vcf file.")

    print(f"Writing unique rows to {output_file}...")
    write_vcf(exp_headers, unique_to_exp, output_file)
    print("Finished writing to file.")

if __name__ == "__main__":
    main()
