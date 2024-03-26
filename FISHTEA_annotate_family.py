import csv

# Function to read TE_rmsk.gtf file and store start and end positions in a dictionary
def read_gtf_file(gtf_file):
    te_positions = {}
    with open(gtf_file, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            chromosome = parts[0]
            attributes = parts[8].split(';')
            te_id_parts = {}
            for attr in attributes:
                try:
                    key, value = attr.strip().split(' ')
                    te_id_parts[key] = value.strip('"')
                except ValueError:
                    continue
            te_id = te_id_parts.get('gene_id')
            family_id = te_id_parts.get('family_id')
            class_id = te_id_parts.get('class_id')
            if te_id and family_id and class_id:
                start = int(parts[3])
                end = int(parts[4])
                te_key = (te_id, family_id, class_id)
                if te_key in te_positions:
                    te_positions[te_key].append((chromosome, start, end))
                else:
                    te_positions[te_key] = [(chromosome, start, end)]
    return te_positions

# Function to match TE IDs from csv file with TE_rmsk.gtf file and create output
def match_and_write(csv_file, gtf_file, output_file):
    te_positions = read_gtf_file(gtf_file)
    with open(csv_file, 'r') as csvfile, open(output_file, 'w') as output:
        csv_reader = csv.reader(csvfile)
        csv_writer = csv.writer(output)
        next(csv_reader)  # Skip the header line
        for row in csv_reader:
            te_id_parts = row[0].split(':')
            if len(te_id_parts) != 3:
                print(f"Invalid TE name format in CSV file: {row[0]}")
                continue
            gene_id = te_id_parts[0]
            family_id = te_id_parts[1]
            class_id = te_id_parts[2]
            te_key = (gene_id, family_id, class_id)
            if te_key in te_positions:
                for chromosome, start, end in te_positions[te_key]:
                    csv_writer.writerow([chromosome, gene_id, family_id, class_id, start, end])


# Example usage
if __name__ == "__main__":
    csv_file = "Testes_Telescope_convert.csv"
    gtf_file = "GRCz11_Ensembl_rmsk_TE.gtf"
    output_file = "sigTEs_positions_telescope_family_TESTES.csv"
    match_and_write(csv_file, gtf_file, output_file)

