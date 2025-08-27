import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from glob import glob
from tqdm import tqdm


def annotate_with_gtex(row, qtl_info: dict, cis_trans: str, gene_or_pheno: int) -> str:
    if type(row['#chr']) == int:
        pos_id = 'chr' + str(row['#chr']) + '_' + str(row['end'])
    else:
        pos_id = row['#chr'] + '_' + str(row['end'])
    # pos_id = row['posID']
    if not pos_id in qtl_info[cis_trans]:
        return '-'
    return ';'.join(sorted(qtl_info[cis_trans][pos_id][gene_or_pheno]))


def annotate(snps_path, qtl_info, phenotypes_info, adastra_info):
    snps_positions = pd.read_table(snps_path)
    snps_positions['posID'] = snps_positions['#chr'].astype(str) + '_' + snps_positions['end'].astype(str)
    if snps_positions[~snps_positions.id.isna()].shape[0] == 0:
        print(f'no ids in {snps_path}, check the SNP annotation')
        return pd.DataFrame()
    # annotation = snps_positions.join(phenotypes_info, how='left', on='id')
    annotation = snps_positions.merge(phenotypes_info, how='left', left_on='id', right_on='RSID')
    annotation = annotation.merge(adastra_info, how='left', left_on='id', right_on='ID')
    
    annotation['eQTL_cis'] = annotation.apply(lambda x: annotate_with_gtex(x, qtl_info, 'cis', 0), axis=1)
    annotation['eQTL_cis_gene'] = annotation.apply(lambda x: annotate_with_gtex(x, qtl_info, 'cis', 1), axis=1)
    annotation['eQTL_trans'] = annotation.apply(lambda x: annotate_with_gtex(x, qtl_info, 'trans', 0), axis=1)
    annotation['eQTL_trans_gene'] = annotation.apply(lambda x: annotate_with_gtex(x, qtl_info, 'trans', 1), axis=1)
    
    return annotation.drop('posID', axis=1)


def annotate_proj(proj_path, out_path, qtl_info, phenotypes_info, adastra_info):
    #out_path = f'annotated'
    os.makedirs(out_path, exist_ok=True)
    snps_paths = glob(f'{proj_path}/*.tsv')
    for file in tqdm(snps_paths):
        print(file)
        annotation = annotate(file, qtl_info, phenotypes_info, adastra_info)
        filepath = file.split('/')[-1]
        if annotation.shape[0] > 0:
            annotation.to_csv(f'{out_path}/{filepath}', index=False, sep='\t', na_rep='-')
    
    
    
# Function to filter DataFrame based on threshold and return row count
def filter_by_threshold(df, threshold, colname):
    df_pval = df[df['fdr_comb_pval'] <= threshold].shape[0]
    df_pval_with_annot = df[(df[colname] != '-') & (df['fdr_comb_pval'] <= threshold)].shape[0]
    if df_pval == 0:
        return 0
    else:
        return df_pval_with_annot/df_pval

def plot_annotation_fraction(df, out_path):
    # Create lists to store thresholds and row counts
    thresholds = np.linspace(0, 5, 100)
    eqtl_cis_counts = [filter_by_threshold(df, 10**-threshold, 'eQTL_cis') for threshold in thresholds]
    ADASTRA_CL_counts = [filter_by_threshold(df, 10**-threshold, 'ADASTRA_CL') for threshold in thresholds]

    plt.plot(thresholds, ADASTRA_CL_counts, label = f'Adastra', c = '#6C8EBF', linewidth= 3)
    plt.plot(thresholds, eqtl_cis_counts, label = f'GTEx', c = '#D6B656', linewidth=3)
        
    plt.xlabel('-log10(FDR)')
    plt.ylabel('Fraction of annotated')
    plt.legend()
    plt.grid(True)
    plt.savefig(out_path)
    plt.close()
    # plt.show()