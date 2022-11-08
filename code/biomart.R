library(biomaRt)
library(tidyr)
library(dplyr)
library("DESeq2")
library("IsoformSwitchAnalyzeR")
packageVersion('IsoformSwitchAnalyzeR')



setwd("/run/user/1002/gvfs/smb-share:server=wfs-medizin.top.gwdg.de,share=ukps-all$/AG-Fornasiero/LAB/Public/Nisha/Transcriptomics/Jena_Nano_Illu/Salmon/Illumina/Transcript_cdna_salmon/")

files <- file.path("Salmon", list.files("Salmon"), "quant.sf")
names(files) <- list.files("Salmon")

salmonQuant <- importIsoformExpression(
    sampleVector = files
)


#biomaRt::listEnsemblArchives() # from here I retrieve the correct name for 103 version

biomaRt::listMarts(host='feb2021.archive.ensembl.org')

Mart102Human <- biomaRt::useMart(host='feb2021.archive.ensembl.org', 
                     biomart='ENSEMBL_MART_ENSEMBL', 
                     dataset='hsapiens_gene_ensembl')

Mart102Mouse <- biomaRt::useMart(host='feb2021.archive.ensembl.org', 
                     biomart='ENSEMBL_MART_ENSEMBL', 
                     dataset='mmusculus_gene_ensembl')
                     
MouseSymbols <- rownames(salmonQuant$counts) %>% unique()    


gtf = getBM(attributes = c("chromosome_name", "ensembl_transcript_id_version", "ensembl_gene_id",
                     "transcript_start", "transcript_end", 
                     "cdna_coding_start", "cdna_coding_end", 
                     "cds_start", "cds_end", 
                     "exon_chrom_start", "exon_chrom_end", 
                     "5_utr_start", "5_utr_end", 
                     "3_utr_start", "3_utr_end" ), 
      mart = Mart102Mouse)
      
write.csv(gtf,"/home/nhk/Desktop/mart_extract.csv")      












                 
