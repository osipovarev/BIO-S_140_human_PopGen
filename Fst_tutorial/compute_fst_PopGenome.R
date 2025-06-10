devtools::install_github("pievos101/PopGenome")
library(PopGenome)

## Locate your work dir
setwd('/Users/osipova/Documents/LabDocs/TA_human_PopGen/Fst_tutorial/')

## Load the data
snp <- readData("popgenome-vcf", format="VCF")

## Define populations
pops <- get.individuals(snp)[[1]]
samples <- pops[seq(1, length(pops), 2)]

before <- samples[1:6]
after <- samples[7:12]
snp <- set.populations(snp, list(before, after))

## Transform object into object divided by sliding window
WINDOW = 10000
popgen_snp <- sliding.window.transform(snp, 
                                       width=WINDOW, jump=WINDOW, 
                                       type=2,
                                       whole.data=TRUE)
# Measurements per window
popgen_snp <- F_ST.stats(popgen_snp)
popgen_fst <- popgen_snp@nucleotide.F_ST[,1]

## Plot
# before_div <- win_snp@nuc.diversity.within[,1] # diversity among population before
# after_div <- win_snp@nuc.diversity.within[,2] # diversity among population after
positions <- seq(1, length(popgen_fst) * WINDOW, WINDOW)
plot(positions, popgen_fst, pch = 20, ylim = c(0, 1), frame.plot = FALSE)