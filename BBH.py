"""
@author: Antonietta Salerno
@date: 4/1/2022
@title: Comparative Genomics Assignment - Best Bi-directional Hits

Note: To run the following lines, BioPython installation is needed.
On the terminal run the following command: conda install biopython
"""
# Import packages
import pandas as pd
import numpy as np
import seaborn as sns 
import matplotlib.pyplot as plt
from Bio.Blast.Applications import NcbiblastpCommandline as cmd
import warnings
warnings.filterwarnings('ignore')

# 1. Run the BLAST query in the forward and reverse direction taking as input the 
    # FASTA files of the proteomes and returning as output a tab-separated file.
blastp_fwd = cmd(query="latimeria.fasta", db="python.fasta", outfmt=6, 
                 out="lat_pyt.tab", max_target_seqs=1)
blastp_bwd = cmd(query="python.fasta", db="latimeria.fasta", outfmt=6, 
                 out="pyt_lat.tab", max_target_seqs=1)
stderr, stdout = blastp_fwd()
stderr, stdout = blastp_bwd()

# 2. Read the files produced as pandas dataframes
pyt_lat = pd.read_csv("pyt_lat.tab", sep="\t", header=None)
lat_pyt = pd.read_csv("lat_pyt.tab", sep="\t", header=None)

# 3. Remove the columns not needed
pyt_lat = pyt_lat.drop(columns=[5,6,7,8,11])
lat_pyt = lat_pyt.drop(columns=[5,6,7,8,11])

# 4. Add headers to forward and reverse results dataframes
headers = ["query", "subject", "identity", "coverage", "qlength", 
           "bitscore", "E-value"]
pyt_lat.columns = headers
lat_pyt.columns = headers

# 5. Calculate the normalised bitscore and add the column to the datasets
pyt_lat['norm_bitscore'] = pyt_lat.bitscore.div(pyt_lat.qlength)
pyt_lat.loc[~np.isfinite(pyt_lat['norm_bitscore']), 'norm_bitscore'] = 0   
    # Solves division by zero problem
lat_pyt['norm_bitscore'] = lat_pyt.bitscore.div(lat_pyt.qlength)
lat_pyt.loc[~np.isfinite(lat_pyt['norm_bitscore']), 'norm_bitscore'] = 0   
    # Solves division by zero problem

# 6. Merge forward and reverse results by inner join, while keeping norm_bitscore 
    # and E-value for both directions
bbh = pd.merge(pyt_lat, lat_pyt[['query', 'subject',"norm_bitscore", "E-value"]],
                left_on='subject', right_on='query',
                how='inner')

# 7. Group duplicate RBH rows, taking the maximum value in each column
bbh = bbh.groupby(['query_x', 'subject_x']).max()

# 8. Sort the values by the normalised bitscore of one of the datasets and save 
    # into a .csv file = 90 hits
bbh.sort_values(by='norm_bitscore_x', ascending=False)
bbh.to_csv("BBH_python_latimeria.csv")

# 9. Explore the distribution of normalised bitscores
f, axes = plt.subplots(1, 2, figsize=(14, 7), sharex=True)
sns.despine(left=True)
sns.distplot(bbh.norm_bitscore_x, color="b", ax=axes[0], 
             axlabel="forward normalised bitscores")
sns.distplot(bbh.norm_bitscore_y, color="g", ax=axes[1], 
             axlabel="reverse normalised bitscores")
plt.savefig("density_norm_bitscore.png")

# 10. Filter by significance (E-value) based on the exploratory plots = 71 hits
bbh = bbh.drop(bbh[bbh["E-value_x"] > 0.001].index)
bbh = bbh.drop(bbh[bbh["E-value_y"] > 0.001].index)

# 11. Filter by quality of the alignment (normalised bitscore) 
    # based on the exploratory plots = 2 hits
bbh = bbh.drop(bbh[bbh["norm_bitscore_x"] < 300].index)
bbh = bbh.drop(bbh[bbh["norm_bitscore_y"] < 40].index)
bbh.to_csv("BBH_python_latimeria_orthologs.csv")
print(bbh.norm_bitscore_x)
