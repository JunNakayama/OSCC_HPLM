# ATAC-seq analysis
# Trim adapters, Alignment and Peakcall were were performed as described in Jason P. Smith, M. Ryan Corces, Jin Xu, Vincent P. Reuter, Howard Y. Chang, and Nathan C. Sheffield. PEPATAC: an optimized pipeline for ATAC-seq data analysis with serial alignments. NAR Genomics and Bioinformatics (2021). DOI: 10.1093/nargab/lqab101

# Differential analysis
library('DiffBind')

samples <- read.csv("config.csv")
db <- dba(sampleSheet=samples)
db <- dba.count(db)
norm <- dba.normalize(db, bRetrieve=TRUE)
db <- dba.contrast(db, minMembers=2)
db <- dba.analyze(db)
tmp <- dba.report(db, th=1, bCounts=TRUE) 
write.csv(tmp, file="DiffBind_report.csv")

