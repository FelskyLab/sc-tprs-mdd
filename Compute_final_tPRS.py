from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
import os
import pandas as pd 

#Go to wd to save results
os.chdir('/external/rprshnas01/netdata_kcni/dflab/team/eoa/tPRS analysis/results/Computed_tPRSs/with_GTex') 

#Enter the path to gene list file, define output file names, and gene weights column
gene_list_file = '../../../Gene_lists/Processed_lists/Sig_MLithwick_avg_gene_list-M.xlsx' 
gene_expr_fpath = "/external/rprshnas01/external_data/uk_biobank/workspace/eadeniyi/Results/Predicted_with_gtex"
wts_col = 'logFC' #name of column containing gene weights in the gene list input file
#Define output file names
Genes_predicted_fn = 'Sig_MLithwick_Gtex_Genes_predicted_M.tsv'
tprs_output_fn = 'Sig_MLithwick_Gtex_final_tPRS_M.txt' 



#Read file from path
if gene_list_file.endswith(".csv"):
    gene_list = pd.read_csv(gene_list_file)
elif gene_list_file.endswith(".xlsx"):
    gene_list = pd.read_excel(gene_list_file, sheet_name=0)
else:
    raise ValueError("Unsupported file format. Please use a .csv or .xlsx file.")

#read in all predicted gene expression files chr 1:22
all_chr = []
for i in range(1,23):
    i = str(i)
    file_name = gene_expr_fpath + "/UKB_predicted_2024_chr" + i + ".txt"
    df = pd.read_csv(file_name, delimiter="\t", index_col='FID')
    all_chr.append(df)

#Function to transpose and process all gene expression files: removes IIDs (column 0), scales genes expr by column, transposes the df, 
#renames cols, extracts gene ids with version number, then drops rownames
def transpose_df(df):
    df = df.iloc[:, 1:]
    scaled_df = scaler.fit_transform(df)
    scaled_df = pd.DataFrame(scaled_df, columns=df.columns, index =df.index)
    transposed_df = scaled_df.transpose()
    transposed_df.columns = ["sample_" + str(i) for i in transposed_df.columns]
    transposed_df['ensembl_gene_id'] = transposed_df.index.str.replace(r'\.\d*', '', regex=True)
    transposed_df.reset_index(drop=True, inplace=True)
    return(transposed_df)

#Apply func to each predicted gene expression each expression df then bind them by rows. Merge the full expr df with the df containing avg dlpfc weights then drop non-matches
t_all_chr = [] #for transposed output
for df in all_chr:
    t_all_chr.append(transpose_df(df))
    
t_all_chr = pd.concat(t_all_chr, ignore_index=True) 
t_all_chr = t_all_chr.merge(gene_list, on='ensembl_gene_id', how="left")
t_all_chr.dropna(subset=[wts_col], inplace=True) #drop unmatched rows (i.e NAs in weights column)

#Write genes expressed info to file
genes_expressed = t_all_chr.loc[:, 'ensembl_gene_id':]
genes_expressed.to_csv(Genes_predicted_fn, sep="\t", index=False)

#Extract names of samples then weight the genes in the column by their effect. Write the result to file
sample_names = t_all_chr.filter(like='sample_').columns
with open(tprs_output_fn, 'w') as file:
    for col_name in sample_names:
        t_all_chr[col_name] = t_all_chr[col_name] * t_all_chr[wts_col]
        tPRS_val = col_name + "\t" + str(t_all_chr[col_name].sum())
        file.write(tPRS_val + '\n')
