### Calculation of Fst the Black Death paper before vs after populations

## Add vcftools to your PATH
Sys.setenv(PATH=paste(Sys.getenv("PATH"), "/Users/osipova/local/Homebrew/bin/", sep=":"))

## Install and load packages
#install.packages("stringr")
library(stringr)

## Locate your work dir
setwd('/Users/osipova/Documents/LabDocs/TA_human_PopGen/Fst_tutorial/')




## Load VCF input
VCF <- 'popgenome-vcf/sampled.chr5.blackDeath.vcf'
WINDOW = 10000
OUT_PREFIX <- 'before_vs_after_Fst_10kb'

## compute window based Fst
system(str_c('vcftools',
             ' --vcf ', VCF,
             ' --weir-fst-pop before.lst ',
             ' --weir-fst-pop after.lst ',
             ' --fst-window-size ', WINDOW,
             ' --out ', OUT_PREFIX))

## Load vcftools output
fst <- read.table(paste0(OUT_PREFIX, '.windowed.weir.fst'), header = TRUE)
head(fst)


## Assign colors based on significance
fst$COLOR <- with(fst, ifelse(WEIGHTED_FST > 0.6, 'red', 'black'))
fst$SIZE <- with(fst, ifelse(WEIGHTED_FST > 0.6, 19, 20))

## Plot
par(mfrow=c(1,1))
plot(fst$BIN_START,
     fst$WEIGHTED_FST,
     pch = fst$SIZE,
     col = fst$COLOR,
     ylim = c(0, 1),
     xlab = 'position',
     ylab = 'Fst',
     frame.plot = FALSE,
     main = 'before vs after the Black Death')


