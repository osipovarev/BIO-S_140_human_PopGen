####################################
### Admixture
setwd('/Users/osipova/Documents/LabDocs/TA_human_PopGen/Admixture_tutorial_Jonathan/')


# PLINK <- '/Users/osipova/Documents/LabDocs/TA_human_PopGen/software/plink/plink'
ADMIXTURE <- '/Users/osipova/Documents/LabDocs/TA_human_PopGen/software/admixture_macosx-1.3.0/admixture '

## run Admixture for k = 1,2,3,4,5
for (i in 1:2){
  system(paste0(ADMIXTURE, '--cv hapmap_r23_subset.bed ', i) )
}

## Load Admixture output
runs<-list()
for (i in 1:5){
  runs[[i]]<-read.table(paste0("hapmap_r23_subset.", i, ".Q"))
}

## Plot runs for all k
par(mfrow=c(4,1))
for (i in 1:4){
  barplot(t(as.matrix(runs[[i]])), col=rainbow(i), ylab="Ancestry", border=NA)
}




####################################
### Fst calculations

# # input the SNP data and the sample names
# snp_matrix   <- read.table("snp_matrix.txt")
# sample_names <- read.table("sample_names.txt")[[1]]
# 
# loci       <- snp_matrix[,1:2]
# colnames(loci) <- c("scaffold", "position")
# 
# # Turn the matrix on its side (rows = individuals, columns = loci)
# snp_matrix <- snp_matrix[,3:ncol(snp_matrix)]
# snp_matrix <- t(snp_matrix)
# 
# # add sample names
# row.names(snp_matrix) <- sample_names
# 
# snp <- new("genlight",
#            snp_matrix,
#            chromosome=loci$scaffold,
#            position=loci$position,
#            pop=as.factor(c(rep("B",6), rep("b",6))))


############### Fst with PopGenome #############

## Locate your work dir
setwd('/Users/osipova/Documents/LabDocs/TA_human_PopGen/Fst_tutorial/')

## add vcftools to yout PATH
Sys.setenv(PATH=paste(Sys.getenv("PATH"), "/Users/osipova/local/Homebrew/bin/", sep=":"))

## install and load packages
devtools::install_github("pievos101/PopGenome")
install.packages("stringr")
library(PopGenome)
library(stringr)

# Load the data
snp <- readData("popgenome-vcf", format="VCF")

# define populations
pops <- get.individuals(snp)[[1]]
samples <- pops[seq(1, length(pops), 2)]

before <- samples[1:6]
after <- samples[7:12]
snp <- set.populations(snp, list(before, after))

## Fst with PopGenome
# Transform object into object divided by sliding window
WINDOW = 10000
popgen_snp <- sliding.window.transform(snp, 
                                    width=WINDOW, jump=WINDOW, 
                                    type=2,
                                    whole.data=TRUE)
# Measurements per window
popgen_snp <- F_ST.stats(popgen_snp)

# A simple plot
popgen_fst <- popgen_snp@nucleotide.F_ST[,1]

# before_div <- win_snp@nuc.diversity.within[,1] # diversity among population before
# after_div <- win_snp@nuc.diversity.within[,2] # diversity among population after
positions <- seq(1, length(popgen_fst) * WINDOW, WINDOW)
plot(positions, popgen_fst, pch = 20, ylim = c(0, 1), frame.plot = FALSE)


## Fst with vcftools
VCF <- 'popgenome-vcf/corrected.head.chr5.indInd.vcf'
OUT_PREFIX <- 'before_vs_after_Fst_10kb'
WINDOW = 10000

system(str_c('vcftools',
              ' --vcf ', VCF, 
              ' --weir-fst-pop before.lst ',
              ' --weir-fst-pop after.lst ', 
              ' --fst-window-size ', WINDOW, 
              ' --out ', OUT_PREFIX))

## Load vcftools output
vcftools_fst <- read.table(paste0(OUT_PREFIX, '.windowed.weir.fst'), header = TRUE)
head(vcftools_fst)

vcftools_fst$COLOR <- with(vcftools_fst, ifelse(WEIGHTED_FST > 0.6, 'red', 'black'))
vcftools_fst$SIZE <- with(vcftools_fst, ifelse(WEIGHTED_FST > 0.6, 19, 20))

par(mfrow=c(1,1))
plot(vcftools_fst$BIN_START,
     vcftools_fst$WEIGHTED_FST, 
     pch = vcftools_fst$SIZE, 
     col = vcftools_fst$COLOR,
     ylim = c(0, 1), 
     xlab = 'position',
     ylab = 'Fst',
     frame.plot = FALSE,
     main = 'before vs after the Black Death')


