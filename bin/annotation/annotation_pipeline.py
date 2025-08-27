import argparse
import pandas as pd
import pickle
import os
from glob import glob
from annotation import annotate, plot_annotation_fraction

def main():
    parser = argparse.ArgumentParser(description="Annotate SNVs tables in a directory.")
    parser.add_argument('--pvalues_dir', required=True, help='Path to directory containing p-values TSV files.')
    parser.add_argument('--annotation_meta_dir', required=True, help='Path to annotation files.')
    args = parser.parse_args()

    # Load metadata once for all annotations
    print("Loading annotation metadata...")
    with open(args.annotation_meta_dir + '/qtl_info.pickle', 'rb') as f:
        qtl_info = pickle.load(f)
    phenotypes_info = pd.read_csv(args.annotation_meta_dir + '/phenotypes_info.csv')
    adastra_info = pd.read_csv(args.annotation_meta_dir + '/adastra_info.csv')
    print("Metadata loaded successfully")

    # Find all TSV files in the pvalues directory
    tsv_pattern = os.path.join(args.pvalues_dir, "*.tsv")
    tsv_files = glob(tsv_pattern)
    
    if not tsv_files:
        print(f"No TSV files found in {args.pvalues_dir}")
        return
    
    print(f"Found {len(tsv_files)} TSV files to process")
    
    # Process each TSV file
    for tsv_file in tsv_files:
        print(f"Processing {tsv_file}")
        basename = os.path.basename(tsv_file).replace('.tsv', '')
        
        # Annotate the file
        annotation = annotate(tsv_file, qtl_info, phenotypes_info, adastra_info)
        
        # Generate output paths
        out_annotated_path = f"{basename}_annotated.tsv"
        out_plot_path = f"{basename}_annotation_plot.png"
        
        annotation = annotation.fillna('-')
        annotation.to_csv(out_annotated_path, index=False, sep='\t', na_rep='-')
        plot_annotation_fraction(annotation, out_plot_path)
        print(f"Created {out_annotated_path} and {out_plot_path}")


if __name__ == '__main__':
    main() 