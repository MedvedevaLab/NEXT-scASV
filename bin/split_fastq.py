#!/usr/bin/env python3
import argparse
import gzip
import os
import sys
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description='Split FASTQ files based on barcodes')
    # Rename: sample_id -> sample (keep --sample_id as backwards-compatible alias)
    parser.add_argument('--sample', dest='sample', required=False, help='Sample (meta.json key / split_table sample column)')
    parser.add_argument('--group', required=True, help='Group to extract')
    parser.add_argument('--split_table', required=True, help='TSV file with barcode, sample, group columns')
    parser.add_argument('--input_dir', required=True, help='Input directory containing FASTQ files')
    parser.add_argument('--file_path', required=True, help='Path to FASTQ files (without .R1/R2.fastq.gz)')
    parser.add_argument('--outdir', required=True, help='Output directory')
    return parser.parse_args()

def read_split_table(split_table_path, sample, group):
    """Read split table and filter for specific sample and group"""
    try:
        df = pd.read_csv(split_table_path, sep='\t')
        filtered_df = df[(df['sample'] == sample) & (df['group'] == group)]
        
        if filtered_df.empty:
            print(f"No barcodes found for sample {sample}, group {group}")
            return []  # Return empty list instead of exiting
            
        return filtered_df['barcode'].tolist()
    except Exception as e:
        print(f"Error reading split table: {e}")
        sys.exit(1)

def process_fastq_files(input_dir, file_path, barcodes, sample, group, outdir):
    """Process FASTQ files and split them by barcode"""
    # Full paths to input FASTQ files
    r1_path = os.path.join(input_dir, f"{file_path}.R1.fastq.gz")
    r2_path = os.path.join(input_dir, f"{file_path}.R2.fastq.gz")
    
    # Check if files exist
    if not os.path.exists(r1_path) or not os.path.exists(r2_path):
        print(f"FASTQ files not found: {r1_path} or {r2_path}")
        sys.exit(1)
    
    # Set up output files
    output_prefix = f"{sample}-{group}"
    r1_out_path = os.path.join(outdir, f"{output_prefix}.R1.fastq")
    r2_out_path = os.path.join(outdir, f"{output_prefix}.R2.fastq")
    
    barcode_count = 0
    total_reads = 0
    
    # Create a set of barcodes for faster lookup
    barcode_set = set(barcodes)
    
    # Open output files
    with open(r1_out_path, 'w') as r1_out, open(r2_out_path, 'w') as r2_out:
        try:
            # Open input files
            with gzip.open(r1_path, 'rt') as r1_in, gzip.open(r2_path, 'rt') as r2_in:
                # Process files line by line
                while True:
                    # Read 4 lines at a time (FASTQ record)
                    r1_header = r1_in.readline().strip()
                    if not r1_header:
                        break  # End of file
                    
                    r1_seq = r1_in.readline().strip()
                    r1_plus = r1_in.readline().strip()
                    r1_qual = r1_in.readline().strip()
                    
                    r2_header = r2_in.readline().strip()
                    r2_seq = r2_in.readline().strip()
                    r2_plus = r2_in.readline().strip()
                    r2_qual = r2_in.readline().strip()
                    
                    total_reads += 1
                    
                    # Extract barcode and UMI
                    barcode = r1_seq[:16]  # First 16 bp is the barcode
                    umi = r1_seq[16:25]    # Next 9 bp is the UMI
                    
                    if barcode in barcode_set:
                        # Modify the headers to include barcode and UMI information
                        read_id = r1_header.split(' ')[0]
                        new_header = f"{read_id}_{barcode}_{umi}"
                        if ' ' in r1_header:
                            new_header += ' ' + r1_header.split(' ', 1)[1]
                        
                        # Write the read to output files with modified headers
                        r1_out.write(f"{new_header}\n{r1_seq}\n{r1_plus}\n{r1_qual}\n")
                        r2_out.write(f"{new_header}\n{r2_seq}\n{r2_plus}\n{r2_qual}\n")
                        barcode_count += 1
        
        except Exception as e:
            print(f"Error processing FASTQ files: {e}")
            sys.exit(1)
    
    print(f"Processed {total_reads} reads, extracted {barcode_count} reads with matching barcodes")
    
    # If no matching reads were found, remove empty output files
    if barcode_count == 0:
        os.remove(r1_out_path)
        os.remove(r2_out_path)
        print(f"No matching reads found for sample {sample}, group {group}")
        return False
    
    return True

def main():
    args = parse_args()
    
    print(f"Starting split for sample {args.sample}, group {args.group}")
    print(f"Input dir: {args.input_dir}")
    print(f"File path: {args.file_path}")
    print(f"Output dir: {args.outdir}")
    
    # Make sure output directory exists
    os.makedirs(args.outdir, exist_ok=True)
    
    # Get barcodes for this sample and group
    barcodes = read_split_table(args.split_table, args.sample, args.group)
    print(f"Found {len(barcodes)} barcodes for sample {args.sample}, group {args.group}")
    
    if not barcodes:
        print(f"No barcodes found. Skipping processing for {args.sample}, {args.group}")
        # Write an empty log file to indicate this sample/group was processed but had no data
        with open(os.path.join(args.outdir, f"{args.sample}_{args.group}_split.log"), 'w') as f:
            f.write(f"No barcodes found for sample {args.sample}, group {args.group}\n")
        sys.exit(0)  # Exit successfully as this is an expected condition
    
    # Process FASTQ files
    success = process_fastq_files(
        args.input_dir,
        args.file_path,
        barcodes,
        args.sample,
        args.group,
        args.outdir
    )
    
    if not success:
        # Write a log file with the failure info
        with open(os.path.join(args.outdir, f"{args.sample}_{args.group}_split.log"), 'w') as f:
            f.write(f"Failed to process sample {args.sample}, group {args.group}: No matching reads\n")
        sys.exit(0)  # Exit successfully to avoid pipeline failure
    else:
        # Write a success log file
        with open(os.path.join(args.outdir, f"{args.sample}_{args.group}_split.log"), 'w') as f:
            f.write(f"Successfully processed sample {args.sample}, group {args.group}\n")

if __name__ == "__main__":
    main()
