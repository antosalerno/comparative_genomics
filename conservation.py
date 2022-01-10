"""
@author: Antonietta Salerno
@date: 8/1/2022
@title: Comparative Genomics Assignment - Functional Conservation in MSA
Note: To run the following lines, BioPython installation is needed.
On the terminal run the following command: conda install biopython
"""

# Import packages
from Bio import AlignIO
import pandas as pd
import math
import scipy.stats as st
import matplotlib.pyplot as plt
from itertools import groupby
from operator import itemgetter

# Retrieve the ClustalW alignment file and parse it
filename = "species_tree.aln"
format = "clustal"
alignment = AlignIO.read(filename, format)

# Retrieve the MSA of a single column
def get_column(position):
    return [sequence[position] for sequence in alignment]

# Calculate the Shannon entropy for a single column (measure of conservation)
def get_entropy(column):
    entropy = 0
    if column.count('-')>20:
        entropy=+9
    entropy+= st.entropy(pd.Series(column).value_counts())
    return entropy

# Create a dictionary containing the Shannon entropy for every position in the MSA
sites = {}
for position in range(0, alignment.get_alignment_length()):
    col = get_column(position)
    sites[position] = get_entropy(col)

# Order the dictionary (most conserved sequences first: the smaller the entropy 
    # the higher the conservation) and save to csv
sites_ordered = sorted(sites.items(), key=lambda x: x[1])
sites_ordered = pd.DataFrame(sites_ordered, columns = ("position", "entropy"))
sites_ordered.to_csv("sites_ordered_entropy.csv")

# Check whether the site is actually conserved
most_conserved = sites_ordered["position"][0]
print("The most conserved position is ", most_conserved, 
      "and has the following amino acid composition: ", 
      ' '.join(get_column(int(most_conserved))))

# Getting regions of highest conservation
minimum = sites_ordered["entropy"] == sites_ordered["entropy"].min()
get_regions = list(sites_ordered[minimum]["position"])
ranges = []
for k, g in groupby(enumerate(get_regions), lambda ix:ix[0]-ix[1]):
    group = list(map(itemgetter(1), g))
    ranges.append((group[0], group[-1]))

# Plot sequence conservation regions by site for the whole MSA 
fig = plt.figure(figsize=(10,5))
ax = fig.add_subplot(111)
x, y = zip(*sites.items()) # unpack a list of pairs into two tuple
plt.plot(x,y)
for i in range(0, len(ranges)):
    plt.axvspan(ranges[i][0], ranges[i][1], color='lightblue', alpha=0.5)

plt.title("Per-site sequence conservation of MSA")
plt.xlabel("Sequence position")
plt.ylabel("Shannon entropy")
plt.show()
fig.savefig("sites_conservation.png", dpi = 300)
