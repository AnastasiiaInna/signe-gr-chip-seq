---
title: 'Allele-specific ChiPseq analysis'
author: 'Anastasiia Hryhorzhevska'
output:
  pdf_document: default
---


```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

1) Creation of a union peak list/set with DiffBind: 
  a) Peak filtering
  b) Sample filtering: Some samples have very low peaks.
  c) union list
https://bioconductor.org/packages/release/bioc/html/DiffBind.html

Currently, peaks are called for each sample individually. Next, you need to determine which peaks are “valid”, e.g. a peak that was called in at least 10% of the samples is valid. Note: VALID peaks should first be assessed separately for DEX and VEHICLE conditions and for the analysis the union of these lists should  be used.

2) Differential analysis. Usually, in ChIPseq you have an input or IgG control  and you would use special tools for the differential analysis. This design is different and differential analysis can be carried out very similar to the corresponding RNAseq data analysis in DEseq or edgeR for treatment condition (veh vs. dex). You can also run differential analysis in DiffBind (which calls DeSeq)

3) Allele-specific binding: generate allele-specific count matrix and perform RasQUAL (R-package)

```{r, warning=FALSE, message=FALSE}
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install(version = "3.12")
}
  
if (!require(DiffBind))
  BiocManager::install("DiffBind")
library(DiffBind)  

library(GenomicFeatures)

if (!require(ChIPseeker))
  BiocManager::install("ChIPseeker")
library(ChIPseeker) 

if (!require(TxDb.Hsapiens.UCSC.hg19.knownGene))
  BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

if(!require(dplyr))
  install.packages("dplyr")
library(dplyr)

```

```{r, include=FALSE}
# Obtaining the sites significantly differentially bound (DB) between the dex samples and veh samples can be done in a five-step script:

# > data = dba(sampleSheet="samplesheeet.csv")
# > data = dba.count(data)
# > data = dba.contrast(data, categories=DBA_TREATMENT)
# > data = dba.analyze(data)
# > data.DB = dba.report(data)

```
1. Read in sampleSheet

```{r, warning=FALSE}
wd <- "/Users/anastasiia_hry/github/mpip/signe_chip-seq/"
setwd(wd)

result.dir.fn <- paste0(wd, "03_result/")

sample.sheet.fn <- paste0(wd, "01_sample_sheets/chipseq_diffbind_sample_sheet_local.csv")
```

```{r}
samples <- read.csv2(sample.sheet.fn, head = T)[, -11]
```
```{r}
head(samples, 5)
```

2. Create and save a DBA object

```{r, warning=FALSE}
dba.db <- dba(sampleSheet = samples)
# dba.save(dba.db, file = "Signe_ChIPseq", dir = result.dir.fn)
```
```{r, warning=FALSE}
# load(paste0(result.dir.fn, "dba_Signe_ChIPseq"))
dba.db
```
This shows how many peaks are in each peakset, as well as total number of unique peaks after merging overlapping ones (44'967) and the default binding matrix of 95 samples by the 24'124 sites that overlap in at least two of the samples.

3. Peaksets overlap rate

To determine which peaks are “valid”, we perform overlap peaksets analysis for each group (dex and veh) separately.

3.1. Dex group

```{r}
dba.overlap(dba.db, dba.db$masks$dex, mode=DBA_OLAP_RATE)
```
47 consensus peaksets were identified for dex group, 37'757 of which is the total number of unique sites, 21'482 - the number of unique sites appearing in at least two peaksets, 16'679 -  the number of sites overlapping in at least three peaksets, etc.

```{r}
dba.overlap(dba.db, dba.db$masks$dex & dba.db$masks$high, mode=DBA_OLAP_RATE)
```
24 consensus peaksets were identified for dex group, 24'355 of which is the total number of unique sites, 11'496 - the number of unique sites appearing in at least two peaksets, 7'477 -  the number of sites overlapping in at least three peaksets, etc.

3.2. Veh group

```{r, warning = FALSE}
dba.overlap(dba.db, dba.db$masks$veh, mode=DBA_OLAP_RATE)
```
48 consensus peaksets were identified for dex group, 21'486 of which is the total number of unique sites, 8'772 - the number of unique sites appearing in at least two peaksets, 5'760 -  the number of sites overlapping in at least three peaksets, etc.

```{r}
dba.overlap(dba.db, dba.db$masks$veh & dba.db$masks$high, mode=DBA_OLAP_RATE)
```

4. Create separate consensus peaksets for each of the two treatments (veh and dex), then take the union of these two peaksets as the overall consensus

Requiring that consensus peaks overlap in at least 20% of the samples in each group results in 3039 sites for the veh group and 11'128 sites for the dex group:

```{r, warning = FALSE}
dba.consensus <- dba.peakset(dba.db, consensus = c(DBA_TREATMENT), minOverlap = 0.2)
dba.consensus
```

```{r}
dba.plotVenn(dba.consensus, dba.consensus$masks$Consensus)
```

5. Merge peakset list

When forming the global binding matrix consensus peaksets, DiffBind first identifies all unique peaks amongst the relevant peaksets. As part of this process, it merges overlapping peaks,
replacing them with a single peak representing the narrowest region that covers all peaks that overlap by at least one base. 

Add a consensus peakset (derived from overlapping peaks in peaksets already present) :

```{r}
dba.consensus.merge <- dba.peakset(dba.consensus, peaks = dba.consensus$masks$Consensus, minOverlap = 1)
dba.consensus.merge
```

Save union peakset :

```{r}
n             <- length(dba.consensus.merge$peaks)
union.peakset <- dba.consensus.merge$peaks[[n]][,1:3]

write.table(union.peakset, file = paste0(result.dir.fn, "union_peakset.tsv"), sep = "\t", row.names = F, col.names = F, quote = F)
```

6. Generate a binding affinity matrix based on read count scores :

```{r, message=FALSE, warning=FALSE}
dba.count.db <- dba.count(dba.db, peaks = union.peakset, bRemoveDuplicates = F, score = DBA_SCORE_READS)  
save(dba.count.db, file = paste0(result.dir.fn, "dba_count.Rda"))

dba.count.norm.db <- dba.count(dba.db, peaks = union.peakset, bRemoveDuplicates = F)  
```

7. Establishing contrast

Next we need to let DiffBind know how we want to group our samples. In our case we will group based on treatment. We also have to set the minMembers parameter to 47 since we have 47 samples in dex and 48 in veh group.

```{r warning=FALSE}
# load(paste0(result.dir.fn, "dba_count.Rda"))
dba.contrast <- dba.contrast(dba.count.db, categories = c(DBA_CONDITION, DBA_TREATMENT), minMembers = 47)
dba.contrast 
```

8. Differential analysis :

8.1. Normalization

A windowing approach that identifies deferentially enriched regions (from csaw package)
```{r, message=FALSE, warning=FALSE}
# dba.contrast.norm <- dba.normalize(dba.contrast, normalize = DBA_NORM_LIB, library=DBA_LIBSIZE_PEAKREADS, method = DBA_ALL_METHODS)
dba.contrast.norm <- dba.normalize(dba.contrast, normalize = DBA_NORM_LIB, background = T, method = DBA_ALL_METHODS)
dba.dif.anal.norm <- dba.analyze(dba.contrast.norm, method = DBA_ALL_METHODS)
```
```{r}
dba.dif.anal.norm
```
Result :

FDR = 5% : 

DeSeq : 1809 of the 4781 are identified as being significantly differentially bound

EdgeR : 3285 of the 4781

```{r}
dba.plotMA(dba.dif.anal.norm, method = DBA_DESEQ2, contrast = 2, bNormalized=F)
```

```{r}
dba.plotMA(dba.dif.anal.norm, method = DBA_DESEQ2, contrast = 2)
```
```{r}
dba.plotMA(dba.dif.anal.norm, method = DBA_DESEQ2, contrast = 2, th = 0.1)
```

8.2. In-built normalization

```{r, warning=FALSE}
dba.dif.anal <- dba.analyze(dba.contrast, method = DBA_ALL_METHODS)
dba.dif.anal
```
Result :

FDR = 5% : 

DeSeq : 1766 of the 4781 are identified as being significantly differentially bound

EdgeR : 3270 of the 4781

```{r}
dba.plotMA(dba.dif.anal, method = DBA_DESEQ2, contrast = 2, bNormalized=T)
```

```{r}
dba.plotMA(dba.dif.anal, method = DBA_DESEQ2, contrast = 2, bNormalized = T, th = 0.1)
```

```{r}
dba.dif.anal <- dba.dif.anal.norm
```

```{r, echo=FALSE, warning=FALSE, message=FALSE, include=FALSE}
dba.plotHeatmap(dba.dif.anal, contrast = 2, attributes = c(DBA_CONDITION, DBA_TREATMENT))
```

```{r}
dba.plotVenn(dba.dif.anal, contrast = 2, method = DBA_ALL_METHODS, main = "Binding Site Overlaps, FDR = 5%")
```

```{r}
# dba.dif.anal.th.01 <- dba.dif.anal
dba.dif.anal$config$th <- 0.10
dba.plotVenn(dba.dif.anal, contrast = 2, method = DBA_ALL_METHODS, main = "Binding Site Overlaps, FDR = 10%")
```



```{r, message=FALSE, warning=FALSE}
# dba.dif.anal.overlap <- dba.overlap(dba.dif.anal, mode = DBA_OLAP_PEAKS, 
#                           contrast = 2, method=dba.dif.anal$config$AnalysisMethod, th = dba.dif.anal$config$th, 
#                           bUsePval = dba.dif.anal$config$bUsePval, 
#                           bCorOnly=TRUE, CorMethod="pearson", 
#                           DataType = dba.dif.anal$config$DataType)

dba.dif.anal.deseq <- dba.dif.anal$DESeq2
# dba.dif.anal.deseq@assays@data@listData$counts

dba.dif.anal.deseq <- dba.analyze(dba.contrast, method = DBA_DESEQ2)
dba.dif.anal.deseq$config$th <- 0.10
```

9. Retrieve the differentially bound sites

```{r}
dex.report <- dba.report(dba.dif.anal.deseq, contrast = 2, method = DBA_DESEQ2, th = dba.dif.anal.deseq$config$th, 
                                                bCounts = T, bCalled = T) # file = paste0(result.dir.fn, "dex_deseq_report.csv"))
write.csv2(dex.report, file = paste0(result.dir.fn, "dex_deseq_report.csv"), row.names = F)

# dex.report <- read.csv2(paste0(result.dir.fn, "dex_deseq_report.csv"))

head(dex.report)

```
- "Conc" shows the mean read concentration over all the samples (the default calculation uses log2 normalized read counts)
- "Conc_dex" - the mean concentration over the samples in dex group
- "Conc_veh" - the mean concentration over the samples in veh group. 
- "Fold" - the log fold changes (LFCs) between the two groups. A positive value indicates increased binding affinity in the DEX group, and a negative value indicates increased binding affinity in the VEH group. 
- normalized count data for individual samples

The number of differentially bound sites that have enriched in DEX samples :

```{r}
sum(dex.report$Fold > 0)
```

The number of differentially bound sites that have enriched in VEH samples :
```{r}
sum(dex.report$Fold < 0)
```
10. Peaks annotation

```{r}
#  This is the refeernce anno used: /ngs/references/human/GRCh38/Homo_sapiens.GRCh38.97.gtf
ref.annot.fn <- paste(wd, "Homo_sapiens.GRCh38.97.gtf", sep = "/")

txDb <- makeTxDbFromGFF(ref.annot.fn, format = "gtf",
                           organism = "Homo sapiens",
                           chrominfo = NULL,
                           dbxrefTag = "gene_id")
# annotate all peaks
dex.anno.report <- ChIPseeker::annotatePeak(dex.report, tssRegion = c(-1000, 1000), TxDb=txDb, verbose = FALSE) #, level = "transcript", overlap = "all" 
dex.anno.report
```

11. Export data

11.1. Export annotated data 
```{r}
annotation.df   <- as.data.frame(dex.anno.report@anno)
write.csv2(annotation.df, file = paste0(result.dir.fn, "dex_deseq_report_annotated.csv"), row.names = F)

anno.rasqual.df <- annotation.df[,c(113, 12:106)] # take only geneID and sampleIDs
write.csv2(anno.rasqual.df, file = paste0(result.dir.fn, "dex_deseq_anno_count_mtrx_for_rasqual.csv"), row.names = F)
```

11.2 Export only dex and veh enriched results

```{r warning=FALSE, include=FALSE}
count.mtrx.df <- as.data.frame(dex.report)
```

Write files for each set of significant regions identified by DeSeq, separating them based on the gain or loss of enrichment. 

DEX:
```{r, warning=FALSE}
dex.enrich.df <- count.mtrx.df  %>% filter(FDR < 0.1 & Fold > 0) 
write.csv2(dex.enrich.df, file = paste0(result.dir.fn, "dex_deseq_report_only_dex_echriched.csv"), row.names = F)

veh.enrich.df <- count.mtrx.df  %>% filter(FDR < 0.1 & Fold < 0) 
write.csv2(veh.enrich.df, file = paste0(result.dir.fn, "dex_deseq_report_only_veh_echriched.csv"), row.names = F)
```

11.3. Export count matrices for RASQUAL

DEX:
```{r}
count.mtrx.dex.rasqual <- read.csv2(paste0(result.dir.fn, "count_mtrx_dex_for_rasqual.csv"))
dex.samples <- c("PEAKSET_ID", c(read.table(paste0(result.dir.fn, "chipseq_dex_samples.txt")))[[1]])
count.mtrx.dex.rasqual.out <- count.mtrx.dex.rasqual[dex.samples]

write.table(count.mtrx.dex.rasqual[2:nrow(count.mtrx.dex.rasqual.out), ], paste0(result.dir.fn, "count_mtrx_dex_for_rasqual.txt"), sep="\t", row.names = F, col.names = F, quote = F)

write.table(count.mtrx.dex.rasqual[2:nrow(count.mtrx.dex.rasqual.out), ], paste0("/Users/anastasiia_hry/github/rasqual/data/signe_chipseq/", "count_mtrx_dex_for_rasqual.txt"), sep="\t", row.names = F, col.names = F, quote = F)
```

VEH:
```{r}
count.mtrx.veh.rasqual <- read.csv2(paste0(result.dir.fn, "count_mtrx_veh_for_rasqual_2.csv"))
 veh.samples <- c("PEAKSET_ID", c(read.table(paste0(result.dir.fn, "chipseq_veh_samples.txt")))[[1]])
 count.mtrx.veh.rasqual.out <- count.mtrx.dex.rasqual[veh.samples]
 
 write.table(count.mtrx.veh.rasqual[2:nrow(count.mtrx.veh.rasqual.out), ], paste0(result.dir.fn, "count_mtrx_veh_for_rasqual.txt"), sep="\t", row.names = F, col.names = F, quote = F)
 
 write.table(count.mtrx.veh.rasqual[2:nrow(count.mtrx.veh.rasqual.out), ], paste0("/Users/anastasiia_hry/github/rasqual/data/signe_chipseq/", "count_mtrx_veh_for_rasqual.txt"), sep="\t", row.names = F, col.names = F, quote = F)
```

12. Plotting

12.1 PCA

```{r}
dba.plotPCA(dba.dif.anal.deseq, method = DBA_DESEQ2, attributes = DBA_CONDITION, label = DBA_TREATMENT, contrast = 2)
```
12.2. Volcano 
```{r}
dba.plotVolcano(dba.dif.anal.deseq, method = DBA_DESEQ2, contrast = 2)
```
12.3. Boxplot
```{r}
dba.plotBox(dba.dif.anal.deseq, method = DBA_DESEQ2, contrast = 2)
```

