# dna <- fasta2DNAbin(file="http://adegenet.r-forge.r-project.org/files/usflu.fasta")
# annot <- read.csv("http://adegenet.r-forge.r-project.org/files/usflu.annot.csv", header=TRUE, row.names=1)


# library(tidyverse)
library(adegenet)
library(ggtree)
library(ape)
# library(phangorn)


## Set up and load data
setwd('/Users/osipova/Documents/LabDocs/TA_human_PopGen/Phylogenetic_tutorial/')
dna <- fasta2DNAbin(file = 'label.subset.viruses.fa')
annotation <- read.csv('label.subset.annotation.tsv', sep='\t')




## Compute distance matrix
D <- dist.dna(dna, model = "TN93")
temp <- as.data.frame(as.matrix(D))
table.paint(temp, cleg=0, clabel.row=.5, clabel.col=.5)




## Build NJ tree
nj_tree <- nj(D)




## Different tree representations 
# circular
ggtree(nj_tree, branch.length="none", layout="circular") +
  geom_tiplab(size=3, color='black')+ 
  ggplot2::xlim(0, 30) 

# rectangular
ggtree(nj_tree, branch.length="none") +
  geom_tiplab(size=3, color='black') + 
  ggplot2::xlim(0, 30) 




## Final NJ tree
ggtree(nj_tree, branch.length="none") + 
  ggtitle("NJ tree") +
  geom_tiplab(size=5, color=annotation$color) + 
  ggplot2::xlim(0, 30) 



## UPGMA
upgma_tree <- as.phylo(hclust(D, method="average"))

ggtree(upgma_tree, branch.length="none") + 
      ggtitle("UPGMA tree") +
      geom_tiplab(size=5, color=annotation$color) + 
      ggplot2::xlim(0, 30) 

