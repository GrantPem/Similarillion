from collections import defaultdict
from Bio import SeqIO

def filter_fasta(input_file, output_file):
    # Dictionary to store sequences by sample name
    sequences = defaultdict(list)

    # Read input FASTA file and store sequences
    for record in SeqIO.parse(input_file, "fasta"):
        sample_name = record.id.rsplit('_', 1)[0]  # Extract sample name
        sequences[sample_name].append(record)

    # Filter sequences
    filtered_sequences = []
    for sample_name, records in sequences.items():
        # Find highest read count
        max_read_count = max(int(record.id.rsplit('_', 1)[1]) for record in records)
        
        # Collect unique sequences for each sample and read count
        unique_sequences = {}
        for record in records:
            read_count = int(record.id.rsplit('_', 1)[1])
            if read_count == max_read_count:
                sequence = str(record.seq)
                if sequence not in unique_sequences:
                    unique_sequences[sequence] = record
                elif len(record.seq) > len(unique_sequences[sequence].seq):
                    unique_sequences[sequence] = record

        filtered_sequences.extend(unique_sequences.values())

    # Remove sequences with read count of 2 or less
    filtered_sequences = [record for record in filtered_sequences if int(record.id.rsplit('_', 1)[1]) > 2]

    # Write filtered sequences to output file
    with open(output_file, "w") as out_handle:
        SeqIO.write(filtered_sequences, out_handle, "fasta")

# Usage example:
input_file = "sucA-Stone_Clipped.fasta"
output_file = "sucA-Stone_Clipped.fasta"
filter_fasta(input_file, output_file)
