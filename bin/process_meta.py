import pandas as pd

df_for_pat= pd.read_csv('/mnt/flashgpu/lab2/Allele-specifity/data/input/seurat_metadata_B001-B0018.csv', )

df_for_pat.index = df_for_pat['Unnamed: 0'].str.split('_L').str.join('_L00')
df_for_pat = df_for_pat[df_for_pat.index.str.contains('RU')]

df_meta= pd.read_csv('/mnt/flashgpu/lab2/Allele-specifity/data/aida_freeze1_meta.csv')
df_meta = df_meta.set_index('Unnamed: 0')
df_meta = df_meta[df_meta.Country == 'RU']

con = pd.merge(df_meta, df_for_pat[['sample']], left_index= True, right_index= True, how= 'left')

con_pass = con[con['sample'].notna()]

print('before:', len(df_for_pat))
print('now freez1:', len(df_meta))
print('concat by people:', len(con_pass))

print(con[con['sample'].isna()]['orig.ident'].value_counts())

