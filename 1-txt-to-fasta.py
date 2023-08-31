import os

# Declare input folder and output folder
input_dir = "./my_sequences"
output_dir= "./my_fasta"

# Check if the folder already exists
if not os.path.exists(output_dir):
    os.mkdir(output_dir)
    print(f"Folder '{output_dir}' created successfully.")
else:
    print(f"Folder '{output_dir}' already exists.")
    
# Loop through each TXT file in the input directory
forward = 0
reverse = 0
total = 0
NNNNN = 0
for filename in os.listdir(input_dir):
    if filename.endswith(".txt"):
        total += 1
        # Construct the full path to the input file
        input_path = os.path.join(input_dir, filename)

        # Read the DNA sequence from the input file
        with open(input_path) as f:
            sequence = f.read().strip()

        # Check if the content of the file is "NNNNN"
        if sequence == "NNNNN":
            NNNNN += 1
            continue  # Skip this file

        # Construct the description for the FASTA file (using the input file name)
        description = filename.replace(".txt", "")
        if "_F" in description:
            forward += 1
        elif "_R" in description:
            reverse += 1

        # Construct the output path for the FASTA file
        output_path = os.path.join(output_dir, filename.replace(".txt", ".fasta"))

        # Write the DNA sequence to the output FASTA file
        with open(output_path, "w") as f:
            f.write(">{}\n{}\n".format(description, sequence))
        
        print(f"Created file: {description}")
print(f"{forward} forward sequences")
print(f"{reverse} reverse sequences")
print(f"{total} total sequences")
print(f"{NNNNN} samples did not return any valid sequence")