Ligand activity analysis on in-vitro PHD2 experiment
================
Javier Perales-Patón - <javier.perales@bioquant.uni-heidelberg.de> -
ORCID: 0000-0003-0780-6683

Here we will perform a NicheNet ligand activity analysis with the
in-vitro signatures from PHD2-KO Macrophages medium conditioned with
fibroblasts (mouse). This consists on how ligands affect gene expression
in co-cultured through ligand-receptor interactions. For the analysis,
it is required to have a clear signature of genes that are responsive of
this interactions and related to the phenotype. Thus, we focused on gene
sets enriched in the phenotype of interest (pro-fibrotic fibroblasts)
and up-regulated ligands from the senders (Macrophage PHD2cKO).

The pipeline of a basic NicheNet analysis consist mainly of the
following steps:

  - 1.  Define a “sender/niche” cell population and a “receiver/target”
        cell population present in the expression data and determine
        which genes are expressed in both populations

  - 2.  Define a gene set of interest: these are the genes in the
        “receiver/target” cell population that are potentially
        affected by ligands expressed by interacting cells (e.g. genes
        differentially expressed upon cell-cell interaction)

  - 3.  Define a set of potential ligands: these are ligands that are
        expressed by the “sender/niche” cell population and bind a
        (putative) receptor expressed by the “receiver/target”
        population

  - 4.  Perform NicheNet ligand activity analysis: rank the potential
        ligands based on the presence of their target genes in the gene
        set of interest (compared to the background set of genes)

  - 5.  Infer top-predicted target genes of ligands that are top-ranked
        in the ligand activity
analysis

## Step 0: Load required packages, NicheNet’s ligand-target prior model and processed expression data of interacting cells

## set env

Define random seed for reproducible analysis.

``` r
# Seed number
set.seed(1234)
# Output directory
OUTDIR <- "./04_ligand_output/"
if(!dir.exists(OUTDIR)) dir.create(OUTDIR);

# Figures
FIGDIR <- paste0(OUTDIR, "/figures/")
knitr::opts_chunk$set(fig.path=FIGDIR)
knitr::opts_chunk$set(dev=c('png','tiff'))
# Data
DATADIR <- paste0(OUTDIR, "/data/")
if(!dir.exists(DATADIR)) dir.create(DATADIR);
```

### Load libraries

Essential packages for the analysis.

``` r
library(nichenetr)
library(tidyverse)
library(limma)

# pre-ranked GSEA
library(fgsea)

# GraphVis (final chunks)
source("../src/graphics.R")
library(DiagrammeR)
```

Ligand-target model:

This model denotes the prior potential that a particular ligand might
regulate the expression of a specific target
gene.

``` r
# ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
ligand_target_matrix = readRDS(("../data/nichenetr/ligand_target_matrix.rds"))
#ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns

# Because expr data is from mouse, we will convert it
colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols()
rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols()

ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]
```

Expression data of interacting cells: Macrophages.

``` r
v <- readRDS("./01_DGE_output/data/v.rds")
expression = v$E
sample_info = v$target # contains meta-information about the samples
```

## Step 1: Define expressed genes in sender and receiver cell populations

Our research question is to prioritize ligands releases by Macrophages
upon PHD2 knock-out perturbation, which triggers hypoxia response, that
is stimulating collagen production by fibroblasts.

First we have to define which genes are expressed in sender cells
(macrophages) and receiver cells (fibroblasts).

``` r
CPM_cutoff <- 10
N_samples <- c("Mac"=sum(grepl("^MC_", sample_info$group)), 
           "Fib"=sum(grepl("^Fib_", sample_info$group)))

expressed_genes_sender = rownames(expression)[rowSums(2^expression[,grep("^MC_", sample_info$group)] > CPM_cutoff) > N_samples["Mac"]*0.5]
# Here we use the DEG analysis to define important ligands
eBay <- readRDS("./01_DGE_output/data/eBay.rds")

topTab_sender <- topTable(eBay, coef="MC_PHD2", number=Inf) # FDR adjust
# Differentially expressed genes from step 01
DEup_sender <- topTab_sender$genes[topTab_sender$P.Value < 0.05 & topTab_sender$logFC > 0]

expressed_genes_receiver = rownames(expression)[rowSums(2^expression[,grep("^Fib_", sample_info$group)] > CPM_cutoff) > N_samples["Fib"]*0.5]

# Check the number of expressed genes: should be a 'reasonable' number of 
# total expressed genes in a cell type, e.g. between 5000-10000 (and not 500 or 20000)
length(expressed_genes_sender)
```

    ## [1] 8857

``` r
length(expressed_genes_receiver)
```

    ## [1] 9271

## Step 2: Define the gene set of interest and a background of genes

As gene set of interest, we consider the genes of which the expression
is possibly affected due to communication with other cells. The
definition of this gene set depends on your research question and is a
crucial step in the use of NicheNet. Because we here want to investigate
how Macrophages regulate the expression of Fibroblasts pro-fibrotic
signature, we will use the leading edge genes from the Matrisome gene
sets enriched in this phenoty as target genes for the analysis. And all
genes expressed in fibroblasts as background of
genes.

``` r
matrisomeDB <- read.table("../data/Matrisome/matrisome_mm_masterlist.tsv", sep="\t", header=TRUE, stringsAsFactors = FALSE)
matrisome <- split(matrisomeDB$Gene.Symbol, matrisomeDB$Category)
matrisome <- matrisome[names(matrisome)!="n/a"]

set.seed(123)
matrisome.res <- fgsea(pathways = matrisome, stats = eBay$t[,"Fib_PHD2"], nperm=10000)
# Define the downstream analysis based on leading edge genes from matrisome
geneset_oi <- unlist(matrisome.res[matrisome.res$padj<0.05, "leadingEdge"]) 
names(geneset_oi) <- NULL

head(geneset_oi)
```

    ## [1] "Col3a1"  "Col5a3"  "Col15a1" "Col28a1" "Col18a1" "Col27a1"

``` r
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
head(background_expressed_genes)
```

    ## [1] "0610009O20Rik" "0610010F05Rik" "0610010K14Rik" "0610012G03Rik"
    ## [5] "0610030E20Rik" "0610037L13Rik"

For the records, we show which genes we are trying to predict via ligand
analysis. This accounts for a total of n=98 genes involved in
    matrisome.

``` r
length(geneset_oi)
```

    ## [1] 98

``` r
print(geneset_oi)
```

    ##  [1] "Col3a1"  "Col5a3"  "Col15a1" "Col28a1" "Col18a1" "Col27a1" "Col5a1" 
    ##  [8] "Col1a1"  "Col5a2"  "Tinagl1" "Fbln5"   "Wisp2"   "Cyr61"   "Smoc2"  
    ## [15] "Thbs1"   "Nov"     "Igfbp4"  "Bmper"   "Mfge8"   "Mgp"     "Npnt"   
    ## [22] "Ltbp1"   "Emilin2" "Ctgf"    "Lama5"   "Smoc1"   "Thbs2"   "Aebp1"  
    ## [29] "Slit2"   "Igfbp6"  "Pcolce"  "Plxnc1"  "Sema7a"  "Sdc3"    "Clec4a1"
    ## [36] "Plxdc1"  "Clec12a" "Clec5a"  "Sema4d"  "Clec4d"  "Clec7a"  "Sema3c" 
    ## [43] "Clec4e"  "C1qa"    "Clec4n"  "Ovgp1"   "C1qc"    "Clec4a2" "C1qb"   
    ## [50] "Clec4a3" "Plxna2"  "Anxa2"   "Sema3a"  "Plxnd1"  "Clec2h"  "Sema4a" 
    ## [57] "Fcna"    "C1qtnf6" "Anxa6"   "Clec10a" "Clec9a"  "Sema4c"  "Igf1"   
    ## [64] "Angptl4" "Cxcl14"  "Ccl3"    "Igf2"    "Cx3cl1"  "Pdgfb"   "Tgfb1"  
    ## [71] "Ccl2"    "Tnf"     "Ccl4"    "Ccl12"   "Ccl11"   "Il16"    "Ccl5"   
    ## [78] "Ebi3"    "Cxcl10"  "Il1rn"   "S100a1"  "Pf4"     "Ccl9"    "Ccl7"   
    ## [85] "Cxcl2"   "Ccl6"    "Cxcl12"  "Tnfsf14" "Ccl17"   "Fstl1"   "Fgf13"  
    ## [92] "Fgf11"   "Tnfsf13" "S100a8"  "S100a13" "S100a16" "Clcf1"   "Osm"

``` r
# Save them
cat(geneset_oi, sep="\n", file=paste0(DATADIR, "/geneset_oi.txt"))
```

## Step 3: Define a set of potential ligands

As potentially active ligands, we will use ligands that are

1.  Expressed by Macrophages and up-regulated upon PHD2 knock-out
    perturbation.
2.  can bind a (putative) receptor expressed by malignant cells.
    Putative ligand-receptor links were gathered from NicheNet’s
    ligand-receptor data
sources.

<!-- end list -->

``` r
# lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
lr_network = readRDS(("../data/nichenetr/lr_network.rds"))
# Transform to mouse
lr_network = lr_network %>% mutate(from = convert_human_to_mouse_symbols(from), 
                   to = convert_human_to_mouse_symbols(to)) %>% 
drop_na()

ligands = lr_network %>% pull(from) %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_sender)

# NOTE: this is based in our experimental design
expressed_ligands <- intersect(expressed_ligands, DEup_sender)

receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = intersect(receptors,expressed_genes_receiver)

lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & 
                         to %in% expressed_receptors) 
head(lr_network_expressed)
```

    ## # A tibble: 6 x 4
    ##   from    to        source         database
    ##   <chr>   <chr>     <chr>          <chr>   
    ## 1 Pdgfa   Pdgfra    kegg_cytokines kegg    
    ## 2 Pdgfa   Pdgfrb    kegg_cytokines kegg    
    ## 3 Plekho2 Pdgfrb    kegg_cytokines kegg    
    ## 4 Hgf     Met       kegg_cytokines kegg    
    ## 5 Tnfsf12 Tnfrsf12a kegg_cytokines kegg    
    ## 6 Itgb1   Vcam1     kegg_cams      kegg

For the records, we show how many ligands (n=25) are considered for the
analysis. Since this is a short list, we show them all. We also show how
many receptors are considered as well.

``` r
## LIGANDS
print(expressed_ligands)
```

    ##  [1] "Il15"    "Pdgfa"   "Plekho2" "Hgf"     "Kitl"    "Tnfsf12" "Itgb1"  
    ##  [8] "Itga9"   "Alcam"   "Spp1"    "Calm1"   "Anxa1"   "Col18a1" "Gpi1"   
    ## [15] "Mmp13"   "Pf4"     "Plau"    "Rtn4"    "Tfpi"    "Glg1"    "Crlf2"  
    ## [22] "Mif"     "Flrt2"   "Pcdh7"   "Nptn"

``` r
cat(paste0("Expressed ligands accounts for n=",length(expressed_ligands)," genes","\n"),
    file=stdout())
```

    ## Expressed ligands accounts for n=25 genes

``` r
# Save them for the records
cat(expressed_ligands, sep="\n", file=paste0(DATADIR, "/expressed_ligands.txt"))

## RECEPTORS
cat(paste0("Expressed receptors accounts for n=",length(expressed_receptors)," genes","\n"),
    file=stdout())
```

    ## Expressed receptors accounts for n=238 genes

``` r
# Save them for the records
cat(expressed_receptors, sep="\n", file=paste0(DATADIR, "/expressed_receptors.txt"))
```

This ligand-receptor network contains the expressed ligand-receptor
interactions. As potentially active ligands for the NicheNet analysis,
we will consider the ligands from this network.

``` r
potential_ligands = lr_network_expressed %>% pull(from) %>% unique()
head(potential_ligands)
```

    ## [1] "Pdgfa"   "Plekho2" "Hgf"     "Tnfsf12" "Itgb1"   "Itga9"

## Step 4: Perform NicheNet’s ligand activity analysis on the gene set of interest

Now perform the ligand activity analysis: in this analysis, we will
calculate the ligand activity of each ligand, or in other words, we will
assess how well each Macrophage-derived ligands can predict the
Fibroblasts signature as compared to the background of expressed genes
(i.e. predict whether a gene belongs to the Fibroblast pro-fibrotic
program or not).

``` r
ligand_activities = predict_ligand_activities(geneset = geneset_oi, 
                          background_expressed_genes = background_expressed_genes, 
                          ligand_target_matrix = ligand_target_matrix, 
                          potential_ligands = potential_ligands)
```

Now, we want to rank the ligands based on their ligand activity. In
Nichenet validation, they showed that the pearson correlation
coefficient (PCC) between a ligand’s target predictions and the observed
transcriptional response was the most informative measure to define
ligand activity. Therefore, we will rank the ligands based on their
pearson correlation coefficient. This allows us to prioritize
Macrophage-released ligands by ranking that statistic.

``` r
ligand_activities %>% arrange(-pearson) 
```

    ## # A tibble: 22 x 4
    ##    test_ligand auroc   aupr pearson
    ##    <chr>       <dbl>  <dbl>   <dbl>
    ##  1 Spp1        0.623 0.0663  0.134 
    ##  2 Anxa1       0.626 0.0369  0.0862
    ##  3 Tnfsf12     0.638 0.0421  0.0787
    ##  4 Il15        0.631 0.0307  0.0604
    ##  5 Tfpi        0.626 0.0212  0.0525
    ##  6 Pf4         0.600 0.0218  0.0511
    ##  7 Mmp13       0.636 0.0183  0.0506
    ##  8 Pdgfa       0.607 0.0190  0.0436
    ##  9 Calm1       0.611 0.0186  0.0415
    ## 10 Col18a1     0.612 0.0175  0.0409
    ## # … with 12 more rows

Previous table reports that few top ligands are good predictors of the
gene set of interest. Then the PCC drops drastically. Becasue of this,
we will choose just the top
10.

``` r
best_upstream_ligands = ligand_activities %>% top_n(10, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
# We show all of them
print(best_upstream_ligands)
```

    ##  [1] "Spp1"    "Anxa1"   "Tnfsf12" "Il15"    "Tfpi"    "Pf4"     "Mmp13"  
    ##  [8] "Pdgfa"   "Calm1"   "Col18a1"

We see here that the performance metrics indicate that the 20 top-ranked
ligands can predict the Fibroblast signature reasonably, this implies
that ranking of the ligands might be accurate. However, it is possible
that for some gene sets, the target gene prediction performance of the
top-ranked ligands would not be much better than random prediction. In
that case, prioritization of ligands will be less trustworthy.

Additional note: we looked at the top 20 ligands here and will continue
the analysis by inferring target genes of these 20 ligands. However, the
choice of looking only at the 20 top-ranked ligands for further
biological interpretation is based on biological intuition and is quite
arbitrary. Therefore, users can decide to continue the analysis with a
different number of ligands. We recommend to check the selected cutoff
by looking at the distribution of the ligand activity values. Here, we
show the ligand activity histogram (the score for the 20th ligand is
indicated via the dashed line).

``` r
# show histogram of ligand activity scores
p_hist_lig_activity = ggplot(ligand_activities, aes(x=pearson)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  # geom_density(alpha=.1, fill="orange") +
  geom_vline(aes(xintercept=min(ligand_activities %>% 
                top_n(10, pearson) %>% pull(pearson))), 
         color="red", linetype="dashed", size=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()
p_hist_lig_activity
```

![](./04_ligand_output//figures/unnamed-chunk-13-1.png)<!-- -->

## Step 5: Infer target genes of top-ranked ligands and visualize in a heatmap

Now we will look into the regulatory potential scores between ligands
and target genes of interest. In this case, we will look at links
between top-ranked regulating ligands and genes. In the ligand-target
heatmaps, we show here regulatory potential scores for interactions
between the 20 top-ranked ligands and following target genes: genes that
belong to the gene set of interest and to the 250 most strongly
predicted targets of at least one of the 10 top-ranked ligands (the top
250 targets according to the general prior model, so not the top 250
targets for this dataset). Consequently, genes of your gene set that are
not a top target gene of one of the prioritized ligands, will not be
shown on the heatmap.

``` r
active_ligand_target_links_df = best_upstream_ligands %>% 
    lapply(get_weighted_ligand_target_links,
           geneset = geneset_oi, 
           ligand_target_matrix = ligand_target_matrix, 
           n = 250) %>% 
bind_rows()

nrow(active_ligand_target_links_df)
```

    ## [1] 51

``` r
head(active_ligand_target_links_df)
```

    ## # A tibble: 6 x 3
    ##   ligand target  weight
    ##   <chr>  <chr>    <dbl>
    ## 1 Spp1   Ccl11  0.00568
    ## 2 Spp1   Ccl3   0.00495
    ## 3 Spp1   Ccl4   0.00465
    ## 4 Spp1   Col3a1 0.00496
    ## 5 Spp1   Ctgf   0.00130
    ## 6 Spp1   Cxcl10 0.00417

For visualization purposes, we adapted the ligand-target regulatory
potential matrix as follows. Regulatory potential scores were set as 0
if their score was below a predefined threshold, which was here the 0.25
quantile of scores of interactions between the 10 top-ranked ligands and
each of their respective top targets (see the ligand-target network
defined in the data frame).

``` r
# Some got weight of NA
active_ligand_target_links_df <- na.omit(active_ligand_target_links_df)

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df,
                                      ligand_target_matrix = ligand_target_matrix,
                                      cutoff = 0.25)

nrow(active_ligand_target_links_df)
```

    ## [1] 51

``` r
head(active_ligand_target_links_df)
```

    ## # A tibble: 6 x 3
    ##   ligand target  weight
    ##   <chr>  <chr>    <dbl>
    ## 1 Spp1   Ccl11  0.00568
    ## 2 Spp1   Ccl3   0.00495
    ## 3 Spp1   Ccl4   0.00465
    ## 4 Spp1   Col3a1 0.00496
    ## 5 Spp1   Ctgf   0.00130
    ## 6 Spp1   Cxcl10 0.00417

``` r
# This is a circunvent somehow it is crashing the dimensionality of the matrix
if(any(!active_ligand_target_links_df$target %in% rownames(active_ligand_target_links))) {
  common <- intersect(active_ligand_target_links_df$target,
                        rownames(active_ligand_target_links))
  which_df <- active_ligand_target_links_df$target %in% common
  active_ligand_target_links_df <- active_ligand_target_links_df[which_df, ]
 active_ligand_target_links <- active_ligand_target_links[common,] 
 rm(common)
}

if(any(!active_ligand_target_links_df$ligand %in% colnames(active_ligand_target_links))) {
    common <- intersect(active_ligand_target_links_df$ligand,
                            colnames(active_ligand_target_links))
      which_df <- active_ligand_target_links_df$ligand %in% common
      active_ligand_target_links_df <- active_ligand_target_links_df[which_df,]
     active_ligand_target_links <- active_ligand_target_links[, common] 
}
```

The putatively active ligand-target links will now be visualized in a
heatmap. The order of the ligands accord to the ranking according to the
ligand activity
prediction.

``` r
order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets = active_ligand_target_links_df$target %>% unique()
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network = vis_ligand_target %>% 
    make_heatmap_ggplot("Prioritized Macrophage-ligands",
                "Responsive target genes from co-cultured Fibroblast", 
                color = "purple",legend_position = "top", x_axis_position = "top",
                legend_title = "Regulatory potential\n") + 
scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.005,0.01)) + # theme(axis.text.x = element_text(face = "italic")) +
theme(legend.key.width = unit(1.2, "cm"),
      legend.text = element_text(family=fontTXT, size=13),
      legend.title = element_text(family=fontTXT, size=14))

p_ligand_target_network
```

![](./04_ligand_output//figures/unnamed-chunk-16-1.png)<!-- -->

## Follow-up analysis 1: Ligand-receptor network inference for top-ranked ligands

One type of follow-up analysis is looking at which receptors of the
receiver cell population (here: fibroblasts) can potentially bind to the
prioritized ligands from the sender cell population (here: Macrophages).

So, we will now infer the predicted ligand-receptor interactions of the
top-ranked ligands and visualize these in a heatmap.

``` r
# get the ligand-receptor network of the top-ranked ligands
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

# get the weights of the ligand-receptor interactions as used in the NicheNet model
# weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
weighted_networks = readRDS(("../data/nichenetr/weighted_networks.rds"))
weighted_networks$lr_sig = weighted_networks$lr_sig %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
weighted_networks$gr = weighted_networks$gr %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
lr_network_top_df = weighted_networks$lr_sig %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

# convert to a matrix
lr_network_top_df = lr_network_top_df %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

# perform hierarchical clustering to order the ligands and receptors
dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
```

Show a heatmap of the ligand-receptor
interactions

``` r
vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% 
    make_heatmap_ggplot("Prioritized PHD2-KO Macrophage ligands",
            "Receptors expressed by Fibroblasts", 
            color = "mediumvioletred", x_axis_position = "top",
            legend_title = "Prior interaction potential")
p_ligand_receptor_network
```

![](./04_ligand_output//figures/unnamed-chunk-18-1.png)<!-- -->

## Follow-up analysis 2: Visualize expression of top-predicted ligands and their target genes in a combined heatmap

NicheNet only considers expressed ligands of sender cells, but does not
take into account their expression for ranking the ligands. This is why
we chose up-regulated ligands in the condition of under study (PHD2cKO).
Then the ranking is based on the potential that a ligand might regulate
the gene set of interest, given prior knowledge. Because it is also
useful to further look into expression of ligands and their target
genes. To conclude the analysis, we generate a combined figure showing
ligand activity, ligand expression, target gene expression and
ligand-target regulatory potential.

#### Load additional packages required for the visualization:

``` r
library(RColorBrewer)
library(cowplot)
library(ggpubr)
```

#### Prepare the ligand activity matrix

``` r
ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)

vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
```

``` r
p_ligand_pearson = vis_ligand_pearson %>% 
  make_heatmap_ggplot("Prioritized Macrophage-ligands","Ligand activity",
  color = "darkorange",legend_position = "top", x_axis_position = "top",
  legend_title = "Pearson Correlation Coef. \n(prediction ability)  ") +
#  scale_fill_gradient2(breaks = c(0, 0.05,0.10,0.15)) + # theme(axis.text.x = element_text(face = "italic")) +
  theme(legend.key.width = unit(1.2, "cm"),
    legend.text = element_text(family = fontTXT, size=13),
    legend.title = element_text(family = fontTXT, size=14))

p_ligand_pearson
```

![](./04_ligand_output//figures/unnamed-chunk-21-1.png)<!-- -->

#### Prepare expression of ligands (sender: macrophages)

``` r
expression_df_Mac = expression[order_ligands, grep("^MC_", v$targets$group)]
colnames(expression_df_Mac) <- gsub("^.*_PHD2_KO","PHD2cKO",colnames(expression_df_Mac))
colnames(expression_df_Mac) <- gsub("^.*_PHD2_WT","WT",colnames(expression_df_Mac))

vis_ligand_Mac_expression = expression_df_Mac
```

``` r
library(RColorBrewer)
color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
vis_ligand_Mac_expression = expression_df_Mac %>% t() %>% scale_quantile() %>% t()

p_ligand_Mac_scaled_expression = vis_ligand_Mac_expression  %>%
make_threecolor_heatmap_ggplot("Ligands","Sender (Macrophage)",
    low_color = color[1],mid_color = color[50], mid = 0.5,
    high_color = color[100], legend_position = "top", x_axis_position = "top" ,
    legend_title = "Scaled expression\n") +
theme(legend.key.width = unit(1.2, "cm"),
    legend.text = element_text(family = fontTXT, size=13),
    legend.title = element_text(family = fontTXT, size=14))
p_ligand_Mac_scaled_expression
```

![](./04_ligand_output//figures/unnamed-chunk-23-1.png)<!-- -->

#### Prepare expression of target genes (receiver: fibroblasts)

``` r
expression_df_target = expression[order_targets, grep("^Fib_", v$target$group)]
colnames(expression_df_target) <- gsub("3T3_PHD2_KO", "Stimulated", colnames(expression_df_target))
colnames(expression_df_target) <- gsub("3T3_PHD2_WT", "Control", colnames(expression_df_target))

vis_target_fib_expression_scaled = expression_df_target %>% t() %>% scale_quantile() 
colnames(vis_target_fib_expression_scaled) <- make.names(colnames(vis_target_fib_expression_scaled)) 
```

``` r
p_target_fib_scaled_expression = vis_target_fib_expression_scaled  %>%
make_threecolor_heatmap_ggplot("Receiver (Fibroblast)","Target",
    low_color = color[1],mid_color = color[50], mid = 0.5,
    high_color = color[100], legend_position = "top", x_axis_position = "top" ,
    legend_title = "Scaled Gene expression\n") +
theme(legend.key.width = unit(1.0, "cm"),
    legend.text = element_text(family = fontTXT, size=10),
    legend.title = element_text(family = fontTXT, size=14))
p_target_fib_scaled_expression
```

![](./04_ligand_output//figures/unnamed-chunk-25-1.png)<!-- -->

#### Combine the different heatmaps in one overview figure

``` r
figures_without_legend = plot_grid(
  p_ligand_pearson + 
      theme(legend.position = "none", 
        text = element_text(family=fontTXT, color="black"),
        axis.ticks = element_blank(),
        axis.title.x = element_text(family=fontTXT, size=14),
        axis.title.y = element_text(family=fontTXT, size=14),
        axis.text = element_text(family=fontTXT, color="black", size=14),
        plot.margin = unit(c(0, 0, 0, 0), "cm")), 
  p_ligand_Mac_scaled_expression + 
      theme(legend.position = "none", 
        legend.text = element_text(family = fontTXT, color="black", size=12),
        text = element_text(family=fontTXT, color="black", size=14),
        axis.ticks = element_blank(),
        axis.title.x = element_text(family=fontTXT, size=14),
        axis.text = element_text(family=fontTXT, color="black", size=14),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) + ylab(""),
  p_ligand_target_network + 
      theme(legend.position = "none",
        text = element_text(family = fontTXT, color="black"),
        axis.title.x = element_text(family=fontTXT, size=14),
        axis.ticks = element_blank(),
        axis.text = element_text(family=fontTXT, color="black", size=14),
        plot.margin = unit(c(0,0,0,0), "cm")) + ylab(""), 
  NULL,
  NULL,
  p_target_fib_scaled_expression + 
      theme(legend.position = "none",
        text = element_text(family = fontTXT, color="black"),
        axis.ticks = element_blank(),
        axis.title.y = element_text(family=fontTXT, size=14),
        axis.text = element_text(family=fontTXT, color="black", size=14),
        plot.margin = unit(c(0,0,0,0), "cm")) + xlab(""), 
  align = "hv",
  nrow = 2,  rel_widths = c(ncol(vis_ligand_pearson)+ 2.0, ncol(vis_ligand_Mac_expression) -0, ncol(vis_ligand_target)) + 2,
  rel_heights = c(nrow(vis_ligand_pearson), nrow(vis_target_fib_expression_scaled) + 1.50)) 

legends = plot_grid(
  as_ggplot(get_legend(p_ligand_pearson)),
  as_ggplot(get_legend(p_ligand_target_network)),
  as_ggplot(get_legend(p_ligand_Mac_scaled_expression)),
#  as_ggplot(get_legend(p_target_fib_scaled_expression)), # We remove it because it is the same as above
  nrow = 1, rel_widths = c(1.1,1,1),
  align = "h")

plot_grid(figures_without_legend, 
          legends, 
          rel_heights = c(10,2), nrow = 2, align = "hv")
```

![](./04_ligand_output//figures/combined_heatmaps_ligand-1.png)<!-- -->

## Assess performance of target prediction

For the top 20 ligands, we will now build a multi-ligand model that uses
all top-ranked ligands to predict whether a gene belongs to the
collagens program of not. This classification model will be trained via
cross-validation and returns a probability for every
gene.

``` r
# change rounds and folds here, to two rounds to reduce time: normally: do multiple rounds
k = 3 # 3-fold
n = 2 # 2 rounds

gs_gene_predictions_top20_list = seq(n) %>% lapply(assess_rf_class_probabilities, 
                             folds = k, 
                             geneset = geneset_oi, 
                             background_expressed_genes = background_expressed_genes, 
                             ligands_oi = best_upstream_ligands, 
                             ligand_target_matrix = ligand_target_matrix)
```

``` r
# get performance: auroc-aupr-pearson
target_prediction_performances_cv = gs_gene_predictions_top20_list %>% 
    lapply(classification_evaluation_continuous_pred_wrapper) %>% 
    bind_rows() %>% mutate(round=seq(1:nrow(.)))
```

What is the AUROC, AUPR and PCC of this model (averaged over
cross-validation rounds)?

``` r
target_prediction_performances_cv$auroc %>% mean()
```

    ## [1] 0.7149865

``` r
target_prediction_performances_cv$aupr %>% mean()
```

    ## [1] 0.07634632

``` r
target_prediction_performances_cv$pearson %>% mean()
```

    ## [1] 0.1607228

Evaluate now whether genes belonging to the gene set are more likely to
be top-predicted. We will look at the top 5% of predicted targets
here.

``` r
# get performance: how many collagen genes and non-collagen-genes among top 5% predicted targets
target_prediction_performances_discrete_cv = gs_gene_predictions_top20_list %>% 
    lapply(calculate_fraction_top_predicted, quantile_cutoff = 0.95) %>% 
    bind_rows() %>% ungroup() %>% mutate(round=rep(1:length(gs_gene_predictions_top20_list), each = 2))
```

What is the fraction of signature that belongs to the top 5% predicted
targets?

``` r
target_prediction_performances_discrete_cv %>% filter(true_target) %>% .$fraction_positive_predicted %>% mean()
```

    ## [1] 0.2258065

What is the fraction of non-collagen genes that belongs to the top 5%
predicted
targets?

``` r
target_prediction_performances_discrete_cv %>% filter(!true_target) %>% .$fraction_positive_predicted %>% mean()
```

    ## [1] 0.0493258

We see that the signature is enriched in the top-predicted target genes.
To test this, we will now apply a Fisher’s exact test for every
cross-validation round and report the average
p-value.

``` r
target_prediction_performances_discrete_fisher = gs_gene_predictions_top20_list %>% 
    lapply(calculate_fraction_top_predicted_fisher, quantile_cutoff = 0.95) 
target_prediction_performances_discrete_fisher %>% unlist() %>% mean()
```

    ## [1] 1.891916e-08

Finally, we will look at which genes from the signature are
well-predicted in every cross-validation round.

``` r
# get top predicted genes
top_predicted_genes = seq(length(gs_gene_predictions_top20_list)) %>% 
    lapply(get_top_predicted_genes, gs_gene_predictions_top20_list) %>% 
    reduce(full_join, by = c("gene","true_target"))
top_predicted_genes %>% filter(true_target)
```

    ## # A tibble: 27 x 4
    ##    gene   true_target predicted_top_target_roun… predicted_top_target_roun…
    ##    <chr>  <lgl>       <lgl>                      <lgl>                     
    ##  1 Ccl11  TRUE        TRUE                       TRUE                      
    ##  2 Ccl3   TRUE        TRUE                       TRUE                      
    ##  3 Col3a1 TRUE        TRUE                       TRUE                      
    ##  4 Cxcl2  TRUE        TRUE                       TRUE                      
    ##  5 Ccl2   TRUE        TRUE                       TRUE                      
    ##  6 Ccl4   TRUE        TRUE                       TRUE                      
    ##  7 Tnf    TRUE        TRUE                       TRUE                      
    ##  8 Ccl12  TRUE        TRUE                       NA                        
    ##  9 Il1rn  TRUE        TRUE                       TRUE                      
    ## 10 C1qb   TRUE        TRUE                       TRUE                      
    ## # … with 17 more rows

## SessionInfo

``` r
sessionInfo()
```

    ## R version 3.6.1 (2019-07-05)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 18.04.3 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
    ## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] e1071_1.7-2        ggpubr_0.2.2       magrittr_1.5      
    ##  [4] cowplot_1.0.0      RColorBrewer_1.1-2 DiagrammeR_1.0.5  
    ##  [7] extrafont_0.17     fgsea_1.10.1       Rcpp_1.0.2        
    ## [10] limma_3.40.6       forcats_0.4.0      stringr_1.4.0     
    ## [13] dplyr_0.8.3        purrr_0.3.2        readr_1.3.1       
    ## [16] tidyr_1.0.0        tibble_2.1.3       ggplot2_3.2.1     
    ## [19] tidyverse_1.2.1    nichenetr_0.1.0    rmarkdown_1.15    
    ## [22] nvimcom_0.9-82    
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] colorspace_1.4-1     ggsignif_0.6.0       ellipsis_0.3.0      
    ##  [4] class_7.3-15         htmlTable_1.13.1     base64enc_0.1-3     
    ##  [7] rstudioapi_0.10      fansi_0.4.0          prodlim_2019.11.13  
    ## [10] lubridate_1.7.4      xml2_1.2.2           codetools_0.2-16    
    ## [13] splines_3.6.1        knitr_1.24           zeallot_0.1.0       
    ## [16] Formula_1.2-3        jsonlite_1.6         pROC_1.16.1         
    ## [19] caret_6.0-85         broom_0.5.2          Rttf2pt1_1.3.8      
    ## [22] cluster_2.1.0        compiler_3.6.1       httr_1.4.1          
    ## [25] backports_1.1.4      assertthat_0.2.1     Matrix_1.2-17       
    ## [28] lazyeval_0.2.2       cli_1.1.0            acepack_1.4.1       
    ## [31] visNetwork_2.0.9     htmltools_0.3.6      tools_3.6.1         
    ## [34] igraph_1.2.4.1       gtable_0.3.0         glue_1.3.1          
    ## [37] reshape2_1.4.3       fastmatch_1.1-0      cellranger_1.1.0    
    ## [40] vctrs_0.2.0          gdata_2.18.0         nlme_3.1-141        
    ## [43] extrafontdb_1.0      iterators_1.0.12     timeDate_3043.102   
    ## [46] gower_0.2.1          xfun_0.9             rvest_0.3.4         
    ## [49] lifecycle_0.1.0      gtools_3.8.1         MASS_7.3-51.4       
    ## [52] scales_1.0.0         ipred_0.9-9          hms_0.5.1           
    ## [55] parallel_3.6.1       yaml_2.2.0           gridExtra_2.3       
    ## [58] rpart_4.1-15         latticeExtra_0.6-28  stringi_1.4.3       
    ## [61] foreach_1.4.7        randomForest_4.6-14  checkmate_1.9.4     
    ## [64] caTools_1.17.1.2     BiocParallel_1.18.1  lava_1.6.6          
    ## [67] rlang_0.4.0          pkgconfig_2.0.3      bitops_1.0-6        
    ## [70] evaluate_0.14        lattice_0.20-38      ROCR_1.0-7          
    ## [73] labeling_0.3         recipes_0.1.9        htmlwidgets_1.3     
    ## [76] tidyselect_0.2.5     plyr_1.8.4           R6_2.4.0            
    ## [79] gplots_3.0.1.1       generics_0.0.2       Hmisc_4.2-0         
    ## [82] pillar_1.4.2         haven_2.1.1          foreign_0.8-72      
    ## [85] withr_2.1.2          survival_2.44-1.1    nnet_7.3-12         
    ## [88] modelr_0.1.5         crayon_1.3.4         utf8_1.1.4          
    ## [91] fdrtool_1.2.15       KernSmooth_2.23-16   grid_3.6.1          
    ## [94] readxl_1.3.1         data.table_1.12.8    ModelMetrics_1.2.2.1
    ## [97] digest_0.6.21        stats4_3.6.1         munsell_0.5.0

``` r
{                                                                                                                                                                                                           
sink(file=paste0(OUTDIR,"/sessionInfo.txt"))
print(sessionInfo())
sink()
}
```
