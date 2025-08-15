# Title: To truncate headers in de novo assemblies for input to RepeatMasker

# Define the directory where your FASTA files are located
INPUT_FOLDER="/gpfs/home/nfv16zpu/scratch/HS_RM/PKG-UOE.AG.ENQ-5178.A.01_24_LITE_R0893_S0001-S0024/rawFastQ/cat_raw_fq/assemblies_megahit/fasta"

# Loop through all files ending in .fasta in the input directory
for file in "$INPUT_FOLDER"/*.fasta; do

  # Check if the file exists to handle cases with no matching files
  if [ -f "$file" ]; then

    # Get the base filename without the path and extension
    base_name=$(basename "$file" .fasta)

    # Define the output file name and path
    output_file="$INPUT_FOLDER/${base_name}_short.fasta"

    echo "Processing $file -> $output_file"

    # Use awk to truncate headers
    awk '/^>/{
      gsub(">", "");
      print ">" substr($0, 1, 50);
      next
    }
    {
      print
    }' "$file" > "$output_file"

  fi
done

echo "All specified FASTA files have been processed."


