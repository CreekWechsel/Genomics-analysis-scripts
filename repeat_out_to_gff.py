import argparse
import sys
import os
import re

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Convert repeat raw formats (.dat, .out, .annot) to GFF3 format.",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="""
Example:
  python repeat_to_gff.py --in ./trf.dat --out trf.gff
  python repeat_to_gff.py --in ./RepeatMasker.out --out RepeatMasker.gff
  python repeat_to_gff.py --in ./Proteinmask.annot --out Proteinmask.gff
        """
    )
    parser.add_argument('--in', dest='input_file', required=True,
                        help='Set input file path (.dat, .out, or .annot)')
    parser.add_argument('--out', dest='output_file', required=True,
                        help='Set output GFF3 file name')
    parser.add_argument('--pre', dest='prefix', default='',
                        help='Set a prefix before repeat element ID')
    parser.add_argument('--verbose', dest='verbose', action='store_true',
                        help='Output running progress information to screen')
    parser.add_argument('--h', action='help',
                        help='Show this help message and exit')
    return parser.parse_args()

def create_marker():
    """Generate a unique marker starting from 1."""
    return 1

def dat_to_gff3(file_path, prefix, out_handle, verbose=False):
    """
    Convert TRF .dat file to GFF3 format.
    """
    if verbose:
        print(f"Processing .dat file: {file_path}")

    chr_name = None
    marker = create_marker()

    try:
        with open(file_path, 'r') as infile:
            for line_num, line in enumerate(infile, 1):
                line = line.strip()
                if line.startswith("Sequence:"):
                    match = re.match(r'^Sequence:\s+(\S+)', line)
                    if match:
                        chr_name = match.group(1)
                        if verbose:
                            print(f"Detected chromosome: {chr_name} at line {line_num}")
                    continue
                if not line:
                    continue
                fields = line.split()
                if len(fields) != 15:
                    if verbose:
                        print(f"Skipping line {line_num}: Expected 15 fields, got {len(fields)}")
                    continue
                start, end = fields[0], fields[1]
                period_size = fields[2]
                copy_number = fields[3]
                consensus_size = fields[4]
                percent_matches = fields[5]
                percent_indels = fields[6]
                score = fields[7]
                consensus = fields[13]
                repeat_seq = fields[14]

                if not chr_name:
                    if verbose:
                        print(f"Skipping line {line_num}: Chromosome name not detected yet.")
                    continue

                element_id = f"{prefix}_TR{marker}" if prefix else f"TR{marker}"
                strand = "+"

                attributes = (
                    f"ID={element_id};ConsensusSize={consensus_size};"
                    f"CopyNumber={copy_number};PercentMatches={percent_matches};"
                    f"PercentIndels={percent_indels};Consensus={consensus};RepeatSeq={repeat_seq}"
                )

                gff3_line = f"{chr_name}\tTRF\tTandemRepeat\t{start}\t{end}\t{score}\t{strand}\t.\t{attributes}\n"
                out_handle.write(gff3_line)
                if verbose and marker % 1000 == 0:
                    print(f"Processed {marker} entries")
                marker += 1
    except FileNotFoundError:
        sys.exit(f"Error: File not found - {file_path}")
    except Exception as e:
        sys.exit(f"Error processing .dat file: {e}")

def out_to_gff3(file_path, prefix, out_handle, verbose=False):
    """
    Convert RepeatMasker .out file to GFF3 format.
    """
    if verbose:
        print(f"Processing .out file: {file_path}")

    marker = create_marker()

    try:
        with open(file_path, 'r') as infile:
            for line_num, line in enumerate(infile, 1):
                line = line.strip()
                if not line or line.startswith("SW") or line.startswith("score"):
                    continue
                fields = line.split()
                if len(fields) < 14:
                    if verbose:
                        print(f"Skipping line {line_num}: Expected at least 14 fields, got {len(fields)}")
                    continue
                # Validate that the first field is numeric
                if not fields[0].isdigit():
                    if verbose:
                        print(f"Skipping line {line_num}: First field is not numeric.")
                    continue
                # Skip simple repeats
                class_family = fields[10]
                if re.search(r'Low|Simple|Satellite', class_family, re.IGNORECASE):
                    if verbose:
                        print(f"Skipping line {line_num}: Class/Family '{class_family}' is excluded.")
                    continue

                score = fields[0]
                perc_div = fields[1]
                perc_del = fields[2]
                perc_ins = fields[3]
                query = fields[4]
                q_begin = fields[5]
                q_end = fields[6]
                strand = fields[8]
                repeat = fields[9]
                class_family = fields[10]
                r_begin = fields[11]
                r_end = fields[12]

                start = q_begin
                end = q_end
                chr_name = query

                strand_symbol = "+" if strand == "+" else "-"

                target_positions = []
                for pos in [fields[11], fields[12], fields[13] if len(fields) > 13 else None]:
                    if pos and not re.search(r'[()]', pos):
                        try:
                            target_positions.append(int(pos))
                        except ValueError:
                            continue
                if len(target_positions) < 2:
                    if verbose:
                        print(f"Skipping line {line_num}: Not enough valid target positions.")
                    continue
                target_start, target_end = sorted(target_positions)[:2]

                element_id = f"{prefix}_TE{marker}" if prefix else f"TE{marker}"
                attributes = (
                    f"ID={element_id};Target={repeat} {target_start} {target_end};"
                    f"Class={class_family};PercDiv={perc_div};PercDel={perc_del};PercIns={perc_ins}"
                )

                gff3_line = f"{chr_name}\tRepeatMasker\tTransposon\t{start}\t{end}\t{score}\t{strand_symbol}\t.\t{attributes}\n"
                out_handle.write(gff3_line)
                if verbose and marker % 1000 == 0:
                    print(f"Processed {marker} entries")
                marker += 1
    except FileNotFoundError:
        sys.exit(f"Error: File not found - {file_path}")
    except Exception as e:
        sys.exit(f"Error processing .out file: {e}")

def annot_to_gff3(file_path, prefix, out_handle, verbose=False):
    """
    Convert RepeatProteinMask .annot file to GFF3 format.
    """
    if verbose:
        print(f"Processing .annot file: {file_path}")

    marker = create_marker()

    try:
        with open(file_path, 'r') as infile:
            for line_num, line in enumerate(infile, 1):
                line = line.strip()
                if not line or line.startswith("pValue"):
                    continue
                fields = line.split()
                if len(fields) < 11:
                    if verbose:
                        print(f"Skipping line {line_num}: Expected at least 11 fields, got {len(fields)}")
                    continue

                pvalue = fields[0]
                query = fields[3]
                start = fields[4]
                end = fields[5]
                strand = fields[6]
                target = fields[7]
                class_family = fields[8]
                try:
                    target_start = int(fields[9])
                    target_end = int(fields[10])
                except ValueError:
                    if verbose:
                        print(f"Skipping line {line_num}: Invalid target positions.")
                    continue
                target_start, target_end = sorted([target_start, target_end])

                element_id = f"{prefix}_TP{marker}" if prefix else f"TP{marker}"
                attributes = (
                    f"ID={element_id};Target={target} {target_start} {target_end};"
                    f"Class={class_family};pValue={pvalue}"
                )

                strand_symbol = strand

                gff3_line = f"{query}\tRepeatProteinMask\tTEprotein\t{start}\t{end}\t{pvalue}\t{strand_symbol}\t.\t{attributes}\n"
                out_handle.write(gff3_line)
                if verbose and marker % 1000 == 0:
                    print(f"Processed {marker} entries")
                marker += 1
    except FileNotFoundError:
        sys.exit(f"Error: File not found - {file_path}")
    except Exception as e:
        sys.exit(f"Error processing .annot file: {e}")

def main():
    args = parse_arguments()

    input_file = args.input_file
    output_file = args.output_file
    prefix = args.prefix.rstrip('_')
    verbose = args.verbose

    if not os.path.isfile(input_file):
        sys.exit(f"Error: Input file does not exist - {input_file}")

    try:
        with open(output_file, 'w') as out_handle:
            out_handle.write("##gff-version 3\n")
            
            if input_file.endswith('.dat'):
                dat_to_gff3(input_file, prefix, out_handle, verbose)
            elif input_file.endswith('.out'):
                out_to_gff3(input_file, prefix, out_handle, verbose)
            elif input_file.endswith('.annot'):
                annot_to_gff3(input_file, prefix, out_handle, verbose)
            else:
                sys.exit("Error: Unsupported file extension. Supported extensions are .dat, .out, .annot")
    except IOError as e:
        sys.exit(f"Error writing to output file: {e}")

    if verbose:
        print(f"Conversion completed. Output written to {output_file}")

if __name__ == "__main__":
    main()
