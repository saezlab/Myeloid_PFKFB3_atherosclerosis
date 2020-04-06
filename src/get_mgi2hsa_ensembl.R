### TITLE : Build ortholog dictionary for mgi (mmu) to human gene ensembl id
### AUTHOR : Perales-Paton, Javier - javier.perales@bioquant.uni-heidelberg.de

library(biomaRt);

#--- Create the BioMart
human = useMart(host="ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
mouse = useMart(host="ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")

#--- Attributes of interest
attributes = c("ensembl_gene_id", "mmusculus_homolog_ensembl_gene", "mmusculus_homolog_perc_id_r1")
orth.mouse1 = getBM(attributes, filters="with_mmusculus_homolog",values=TRUE,
                    mart = human, uniqueRows=TRUE)

orth.mouse_keys <- getBM(c("hgnc_symbol","ensembl_gene_id"),
                         filter="ensembl_gene_id",values=orth.mouse1$ensembl_gene_id,
                         mart = human, uniqueRows=TRUE)

orth.mouse_keys2 <- getBM(c("ensembl_gene_id","mgi_symbol"),
                          filter="ensembl_gene_id",
                          values=orth.mouse1$mmusculus_homolog_ensembl_gene,
                          mart = mouse, uniqueRows=TRUE)

# Table of orthologs
orth.mouse <- merge(orth.mouse1, orth.mouse_keys2, 
                    by.x="mmusculus_homolog_ensembl_gene",
                    by.y="ensembl_gene_id")
orth.mouse$mgi_symbol[orth.mouse$mgi_symbol==""] <- NA
orth.mouse <- na.omit(orth.mouse)

#-  Create dictionary (list of vectors)
mmu2hsa.dic <- split(orth.mouse$ensembl_gene_id, orth.mouse$mgi_symbol)
mmu2hsa.dic <- lapply(mmu2hsa.dic, function(z) unique(z))

# Create another one for gene symbols
orth.mouse <- merge(orth.mouse1,orth.mouse_keys,by.x="ensembl_gene_id",by.y="ensembl_gene_id")
orth.mouse <- merge(orth.mouse, orth.mouse_keys2, 
                    by.x="mmusculus_homolog_ensembl_gene",
                    by.y="ensembl_gene_id")
orth.mouse$mgi_symbol[orth.mouse$mgi_symbol==""] <- NA
orth.mouse <- na.omit(orth.mouse)

#-  Create dictionary (list of vectors)
mmu2hsa.dic2 <- split(orth.mouse$hgnc_symbol, orth.mouse$mgi_symbol)
mmu2hsa.dic2 <- lapply(mmu2hsa.dic2, function(z) unique(z))

hsa2mmu.dic <- split(orth.mouse$mgi_symbol, orth.mouse$hgnc_symbol)
hsa2mmu.dic <- lapply(hsa2mmu.dic, function(z) unique(z))

#--- Save it
if(!dir.exists("./data/Gene_annotation")) dir.create("./data/Gene_annotation");

saveRDS(mmu2hsa.dic, file = "./data/Gene_annotation/mgi2ensembl_mmu2hsa.rds")
saveRDS(mmu2hsa.dic2, file = "./data/Gene_annotation/mgi2hgnc_mmu2hsa.rds")
saveRDS(hsa2mmu.dic, file = "./data/Gene_annotation/hgnc2mgi_hsa2mmu.rds")

