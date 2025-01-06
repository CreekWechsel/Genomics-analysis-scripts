import sys
import re
import argparse

def load_rename_table(rename_table_file):
    """
    Load the rename table into a dictionary.
    
    Args:
        rename_table_file (str): Path to the rename table file.
    
    Returns:
        dict: A dictionary mapping old IDs to new IDs.
    """
    rename_dict = {}
    try:
        with open(rename_table_file, 'r') as file:
            for line in file:
                columns = line.strip().split("\t")
                if len(columns) >= 2:
                    old_id = columns[0]
                    new_id = columns[-1]
                    rename_dict[old_id] = new_id
    except FileNotFoundError:
        print(f"Error: File {rename_table_file} not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error loading rename table: {e}")
        sys.exit(1)
    
    return rename_dict

def rename_fasta(fasta_file, rename_dict):
    """
    Rename IDs in a FASTA file using the provided rename dictionary.
    
    Args:
        fasta_file (str): Path to the FASTA file to process.
        rename_dict (dict): Dictionary containing old-to-new ID mappings.
    """
    try:
        with open(fasta_file, 'r') as infile, open(f"{fasta_file}.renamed", 'w') as outfile:
            for line in infile:
                if line.startswith(">"):  # Only process header lines
                    for old_id, new_id in rename_dict.items():
                        pattern = r'\b{}\b'.format(re.escape(old_id))
                        line = re.sub(pattern, new_id, line)
                outfile.write(line)
    except FileNotFoundError:
        print(f"Error: File {fasta_file} not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error processing FASTA file: {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="Rename IDs in a FASTA file based on a rename table.")
    parser.add_argument('fasta_file', help="Input FASTA file")
    parser.add_argument('rename_table_file', help="Table mapping old IDs to new IDs")
    
    args = parser.parse_args()
    
    # Load rename table
    rename_dict = load_rename_table(args.rename_table_file)
    
    # Rename IDs in FASTA file
    rename_fasta(args.fasta_file, rename_dict)
    
    print(f"Renaming completed. Output saved to {args.fasta_file}.renamed")

if __name__ == "__main__":
    main()
