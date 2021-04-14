### TITLE : Seurat-attached R functions
### AUTHOR : Perales-Paton, Javier - javier.perales@bioquant.uni-heidelberg.de
### DESCRIPTION : These functions are pretty similar to those used in the tutorial
###       but just a wrapper to facilitate readibility of the scripts.

# Get a SeuratObject from 10x data
# Also calculates the percentage of mitocondrial genes
getSeuratObject <- function(path, project_name, min.cells=3, min.features=200, mt.pattern="^MT-") {
  mat10x <-  Read10X(data.dir = path)
  Seurat.obj <- CreateSeuratObject(counts = mat10x, project = project_name, min.cells = min.cells, min.features = min.features)
  
  # Add mitocondrial genes
  Seurat.obj[["percent.mt"]] <- PercentageFeatureSet(Seurat.obj, pattern = mt.pattern)
  
  # Return object
  return(Seurat.obj)
}

# Save the typical violin plot for the lime QC from Seurat.
QCVln <- function(SeuratObject, outdir="./figs/QC") {
  png(paste0(outdir,"/QCvln_",Project(SeuratObject),".png"), width = 1600*3, height = 800*3, res=280)
  print(VlnPlot(SeuratObject, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
  dev.off()
}

# Save the Scatter plot to decide thresholds to remove doublets and death cells (mitochondrial genes)
# before subsetting
QCscatter <- function(SeuratObject, outdir="./figs/QC") {
  png(paste0(outdir,"/QCscatter_",Project(SeuratObject),".png"), width = 1600*3, height = 800*3, res=280)
  print(CombinePlots(plots = list(FeatureScatter(SeuratObject, feature1 = "nCount_RNA", feature2 = "percent.mt"),
                                  FeatureScatter(SeuratObject, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")))
  )
  dev.off()
}

# Scale data for all genes
Seurat_scaledata <- function(SeuratObject) {
  all_genes <- rownames(SeuratObject)
  SeuratObject <- ScaleData(SeuratObject, features = all_genes)
  
  return(SeuratObject)
}

# Show clustertree given a set of set of resolutions.
saveClusterTree <- function(SeuratObject, outdir="./figs/Cluster", prefix="RNA_snn_res.") {
  png(paste0(outdir,"/clustree_",Project(SeuratObject),".png"), width = 2000*3, height = 1000*3, res=280)
  print(clustree(SeuratObject, prefix= prefix))
  dev.off()
}

VarFeatPlot <- function(SeuratObject, Ngenes, outdir="./figs/Cluster") {
  png(filename = paste0(outdir,"/VariableGenes_",Project(SeuratObject),".png"), width = 800*3, height = 800*3, res=280)
  print(LabelPoints(VariableFeaturePlot(SeuratObject), points=VariableFeatures(S)[1:Ngenes],repel = TRUE))
  dev.off()
}

# Plot ModuleScore given a set of genes
PlotModuleScore <- function(SeuratObject, genes, reduction="tsne", tag=NULL) {
  cat("[INFO] : Intersect space with set of genes:\n", file = stdout())
  cat(paste0("\t #",sum(genes %in% rownames(SeuratObject))," out of ",length(genes),"\n",file=stdout()))
  isec <- intersect(rownames(SeuratObject), genes)
  if(length(isec)==0) {
    cat("[ERROR] : No genes are overlapping.\n",file = stdout())
  } else {
  
    SeuratObject <- AddModuleScore(SeuratObject, features=list(c(isec)), name = tag)
    # if(!is.null(tag)) {
    #   colnames(SeuratObject@meta.data)[which(colnames(tmp@meta.data)=="Cluster1")] <- tag
    # } else {
    #   tag <- "Cluster1"
    # }
    print(FeaturePlot(SeuratObject, features=grep(tag,colnames(SeuratObject@meta.data),value=TRUE), reduction=reduction))
  }
  
}

# Enhanced DoHeatmap function
DoHeatmap2 <- function(SeuratObject, GSC, assay="RNA", res=0.5, show_hr=TRUE) {
  library(ComplexHeatmap)
  
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  
  gg_color_hue2 <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 35, c = 100)[1:n]
  }
  
  
  genes <- unlist(geneIds(GSC))
  genes.cols <- unlist(sapply(names(GSC), function(abbn) rep(abbn, length(geneIds(GSC)[[abbn]]))))
  mat <- SeuratObject@assays[[assay]]@scale.data
  
  if(is.null(res)) {
    cl <- as.character(SeuratObject@meta.data[,"seurat_clusters"])
  } else {
    cl <- as.character(SeuratObject@meta.data[,paste0(assay,"_snn_res.",res)]) 
  }
  
  # Reorder
  ord <- order(as.numeric(cl), decreasing = FALSE)
  mat <- mat[,ord]
  cl <- cl[ord]
  
  cl.cols <- setNames(gg_color_hue(length(unique(cl))),unique(as.character(sort(as.numeric(cl)))))
  
  common_genes <- intersect(genes,rownames(mat))
  diff_genes <- setdiff(genes, rownames(mat))
  
  mat2 <- rbind(mat[common_genes,],
                matrix(NA, nrow=length(diff_genes), ncol=ncol(mat),
                       dimnames=list(diff_genes, colnames(mat)))
  )
  mat2 <- mat2[genes,]
  
  hc <- HeatmapAnnotation(df=data.frame("cluster"=cl),col = list("cluster"=cl.cols), show_annotation_name = FALSE,
                          show_legend = FALSE,
                          annotation_legend_param= list(legend_height = unit(8, "cm"),
                                                        grid_width = unit(5, "mm"),
                                                        title_gp=gpar(fontsize=16),
                                                        # at=c(-2.5,-2,-1,0,1,2,2.5),
                                                        # labels = c("","-2","-1","0","1","2",""),
                                                        labels_gp = gpar(fontsize = 14)))
  
  # hr <- rowAnnotation(df=data.frame(markers=genes.cols), show_annotation_name = FALSE)
  hr_classes <- gsub("\\-?[0-9]+$","",genes.cols)
  hr_classes <- gsub("^Modulated_","",hr_classes)
  hr_classes <- gsub("^Healthy_","",hr_classes)
  hr_classes <- gsub("^(CD8|Mixed|CXCR6)_","",hr_classes)
  hr_classes <- gsub("^(M2_?|ResLike|Inflammatory_?|TREM_high)-","",hr_classes)
  
  
  hr <- rowAnnotation(df=data.frame("type"=hr_classes), show_annotation_name = FALSE,
                      col=list("type"=setNames(gg_color_hue2(length(unique(hr_classes))),
                                               unique(hr_classes))),
                      annotation_legend_param= list(legend_height = unit(4, "cm"),
                                                    grid_width = unit(5, "mm"),
                                                    title_gp=gpar(fontsize=16),
                                                    # at=c(-2.5,-2,-1,0,1,2,2.5),
                                                    # labels = c("","-2","-1","0","1","2",""),
                                                    labels_gp = gpar(fontsize = 14)))
  
  f1 <-  circlize::colorRamp2(c(-2,0,+2), c("purple", "black", "yellow"))
  hp <- Heatmap(mat2, cluster_rows = FALSE, cluster_columns = FALSE,col = f1,
                name="Expression",
                top_annotation = hc, bottom_annotation = hc,
                split=factor(genes.cols, levels=unique(genes.cols)),row_title_rot = 0,row_gap = unit(1.5, "mm"), 
                column_split = factor(cl, levels=unique(cl)), column_title_rot=00, column_gap = unit(0.7, "mm"),
                # left_annotation = rowAnnotation(foo=anno_block(gpar(fill=table(genes.cols)[unique(genes.cols)]),
                #                                                labels=unique(genes.cols),
                #                                                labels_gp=gpar(col="white",fontsize=10))),
                heatmap_legend_param= list(legend_height = unit(4, "cm"),
                                           title_gp=gpar(fontsize=16),
                                           # at=c(-2.5,-2,-1,0,1,2,2.5),
                                           # labels = c("","-2","-1","0","1","2",""),
                                           labels_gp = gpar(fontsize = 15)),
                show_column_names = FALSE, row_names_side = "left",
                row_names_gp = gpar(fontsize=7))
  if(show_hr) {
    hh <- hp + hr 
  } else {
    hh <- hp
  }
  return(hh)
}


# Enhanced DoHeatmap function
DoHeatmap3 <- function(SeuratObject, GSC, assay="RNA", res=0.5, show_hr=TRUE, row_names_size=7, 
		       dat="scale", column_title_size=14, fontfamily="Arial", legend_nrow=1) {
  library(ComplexHeatmap)
  
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  
  gg_color_hue2 <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 35, c = 100)[1:n]
  }
  
  gg_color_hue3 <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 95, c = 100)[1:n]
  }
  
  
  genes <- unlist(geneIds(GSC))
  genes.cols <- unlist(sapply(names(GSC), function(abbn) rep(abbn, length(geneIds(GSC)[[abbn]]))))
  if(dat=="scale") {
	mat <- SeuratObject@assays[[assay]]@scale.data
  	f1 <-  circlize::colorRamp2(c(-2,0,+2), c("purple", "black", "yellow"))
  } else if(dat=="data") {
	mat <- as.matrix(SeuratObject@assays[[assay]]@data)
  	f1 <-  circlize::colorRamp2(c(-2,0,+2), c("purple", "black", "yellow"))
  } else {
	  stop("ERROR: data type is not supported")
  }
  
  if(is.null(res)) {
     cl <- as.character(SeuratObject@meta.data[,"seurat_clusters"])
  } else if(res=="Idents") {
     cl <- Idents(SeuratObject)
  } else {
    cl <- as.character(SeuratObject@meta.data[,paste0(assay,"_snn_res.",res)]) 
  }
  
  origin <- SeuratObject@meta.data$orig.ident
  
  # Reorder
  ord <- order(as.numeric(cl), decreasing = FALSE)
  mat <- mat[,ord]
  
  cl <- cl[ord]
  if(is.null(res)) {
  	cl.cols <- setNames(gg_color_hue(length(unique(cl))),unique(as.character(sort(as.numeric(cl)))))
  } else if(res=="Idents") {
  	cl.cols <- setNames(gg_color_hue(length(unique(cl))),levels((cl)))
  } else {
  	cl.cols <- setNames(gg_color_hue(length(unique(cl))),unique(as.character(sort(as.numeric(cl)))))
  }
  
  origin <- origin[ord]
  if(length(unique(origin)) > 2) {
	  origin.cols <- setNames(gg_color_hue3(length(unique(origin))),
                          unique(as.character(origin)))
  } else if (length(unique(origin))==2) {
	 origin.cols <- setNames(c("grey","red"), unique(as.character(origin)))
  }
		
  common_genes <- intersect(genes,rownames(mat))
  diff_genes <- setdiff(genes, rownames(mat))
  
  mat2 <- rbind(mat[common_genes,],
                matrix(NA, nrow=length(diff_genes), ncol=ncol(mat),
                       dimnames=list(diff_genes, colnames(mat)))
  )
  mat2 <- mat2[genes,]
  
  hc <- HeatmapAnnotation(df=data.frame("cluster"=as.character(cl), "origin"=origin, row.names=colnames(mat2)),
                          col = list("cluster"=cl.cols, "origin"=origin.cols),
                          show_annotation_name = TRUE,annotation_name_side = "left",
                          show_legend = c(TRUE,TRUE),
                          annotation_legend_param= list(
							cluster=list(legend_height = unit(8, "cm"),
								     grid_width = unit(5, "mm"),
								     title_gp=gpar(fontsize=16, 
										   fontfamily=fontfamily),
								     nrow=legend_nrow,
								     # at=c(-2.5,-2,-1,0,1,2,2.5),
								     # labels = c("","-2","-1","0","1","2",""),
								     labels_gp = gpar(fontsize = 14,
										      fontfamily = fontfamily)),
							"origin"=list(legend_height = unit(8, "cm"),
								     grid_width = unit(5, "mm"),
								     title_gp=gpar(fontsize=16,
										   fontfamily=fontfamily),
								     nrow=1,
								     # at=c(-2.5,-2,-1,0,1,2,2.5),
								     # labels = c("","-2","-1","0","1","2",""),
								     labels_gp = gpar(fontsize = 14,
										      fontfamily = fontfamily))

			  ))
  
  # hr <- rowAnnotation(df=data.frame(markers=genes.cols), show_annotation_name = FALSE)
  hr_classes <- gsub("\\-?[0-9]+$","",genes.cols)
  hr_classes <- gsub("^Modulated_","",hr_classes)
  hr_classes <- gsub("^Healthy_","",hr_classes)
  hr_classes <- gsub("^(CD8|Mixed|CXCR6)_","",hr_classes)
  hr_classes <- gsub("^(M2_?|ResLike|Inflammatory_?|TREM_high)-","",hr_classes)
  
  
  hr <- rowAnnotation(df=data.frame("type"=hr_classes), show_annotation_name = FALSE,
                      col=list("type"=setNames(gg_color_hue2(length(unique(hr_classes))),
                                               unique(hr_classes))),
                      annotation_legend_param= list(legend_height = unit(8, "cm"),
                                                    grid_width = unit(5, "mm"),
                                                    title_gp=gpar(fontsize=16,
								  fontfamily=fontfamily),
                                                    # at=c(-2.5,-2,-1,0,1,2,2.5),
                                                    # labels = c("","-2","-1","0","1","2",""),
                                                    labels_gp = gpar(fontsize = 14,
								     fontfamily=fontfamily)))
  # Choose main legend scale
  if(dat=="scale") {
  	f1 <-  circlize::colorRamp2(c(-2,0,+2), c("purple", "black", "yellow"))
  } else if(dat=="data") {
  	f1 <-  circlize::colorRamp2(c(min(as.vector(mat2)),
					max(as.vector(mat2))),
				      c("purple","yellow"))
  } else {
	  stop("ERROR: data type is not supported")
  }
  

  hp <- Heatmap(mat2, cluster_rows = FALSE, cluster_columns = FALSE,col = f1, na_col="black",
                name="Expression",
                top_annotation = hc,# bottom_annotation = hc,
                split=factor(genes.cols, levels=unique(genes.cols)),row_title_rot = 0,row_gap = unit(1.5, "mm"), 
                column_split = factor(cl, levels=unique(cl)), column_title_rot=00, column_title_gp = gpar(fontsize=column_title_size,fontfamily=fontfamily), column_gap = unit(0.7, "mm"),
                # left_annotation = rowAnnotation(foo=anno_block(gpar(fill=table(genes.cols)[unique(genes.cols)]),
                #                                                labels=unique(genes.cols),
                #                                                labels_gp=gpar(col="white",fontsize=10))),
                heatmap_legend_param= list(legend_height = unit(4, "cm"),
                                           title_gp=gpar(fontsize=16, fontfamily=fontfamily),
#                                            at=c(-2.5,-2,-1,0,1,2,2.5),
#                                            labels = c("","-2","-1","0","1","2",""),
# 					   grid_border="white",
                                           labels_gp = gpar(fontsize = 14, fontfamily=fontfamily)),
                show_column_names = FALSE, row_names_side = "left",
                row_names_gp = gpar(fontsize=row_names_size, fontfamily=fontfamily))
  if(show_hr) {
    hh <- hp + hr 
  } else {
    hh <- hp
  }
  return(hh)
}

## VlnPlot split by stim for integrated seurat object. Better aesthetics
VlnPlot.stim <- function(S, meta.feature, ylabTXT="", fontTXT="Arial") {
	dat <- S@meta.data[, c("stim",meta.feature)]
	colnames(dat) <- c("stim", "feat")

	# sample size
	sample_size = dat %>% group_by(stim) %>% summarize(num=n())
	# Plot
	dat <- dat %>%
	  left_join(sample_size) %>%
	  mutate(myaxis = paste0(stim, "\n", "n=", num))
	dat$myaxis <- factor(dat$myaxis, levels=sapply(levels(dat$stim), function(z) grep(z, unique(dat$myaxis), value=TRUE)))
	
	  ggplot(dat, aes(x=myaxis, y=feat, fill=stim)) +
	    geom_violin(width=1.0, lwd=1.2) +
	    geom_boxplot(width=0.1, 
			 outlier.size=0, 
			 color="grey60", fill="grey", alpha=0.4) +
	    #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.3) +
	    geom_jitter(shape=16, position=position_jitter(0.4)) +
	    scale_fill_manual(values=c("white","white")) +
	   # scale_fill_viridis(discrete = TRUE) +
	    xlab("") + ylab(ylabTXT) +
	    theme_cowplot() +
	    theme(
	      axis.text.x = element_text(family=fontTXT, size=28),
	      axis.text.y = element_text(family=fontTXT, size=24),
	      axis.title = element_text(family=fontTXT, size=22),
	      legend.position="none",
	    )
}


#' Save markers as individual CSV files 
#' 
#' @param x a list from loop-wilcox or genesorteR::sortGenes()
#' @param OUTDIR the output directory path
#' @examples
#' saveMarkers.CSV(x, "/output")
saveMarkers.CSV <- function(x, OUTDIR) {
	## Folder to store data
	MARKERS_OUTDIR <- paste0(OUTDIR,"/Markers")
	if(! dir.exists(MARKERS_OUTDIR)) dir.create(MARKERS_OUTDIR, recursive = TRUE)

	if(is.list(x)) {
		if(all(c("condGeneProb", "postClustProb", "specScore") %in% names(x))) {
			# Genesorter
			write.table(as.matrix(x$specScore), paste0(MARKERS_OUTDIR,"/specScore.tsv"),
				    sep="\t",row.names = TRUE,col.names = NA,quote=FALSE)

			write.table(as.matrix(x$condGeneProb), paste0(MARKERS_OUTDIR,"/condGeneProb.tsv"),
				    sep="\t",row.names = TRUE,col.names = NA,quote=FALSE)
		} else {
			# List of DEGs via wilcox
		    for(idx in names(x)) {
			    firstgofirst <- c("cluster", "gene")
			    remaining <- setdiff(colnames(x[[idx]]), firstgofirst)
			    col_names <- c(firstgofirst, remaining)
			write.table(x[[idx]][,col_names],
				    file = paste0(MARKERS_OUTDIR,"/cluster",idx,".tsv"),
		    		    sep="\t",col.names = TRUE, row.names = FALSE, quote=FALSE)
		    }
		}
	} else {
		stop("ERROR: Only list is allowed\n")
	}
}


#' Save markers as in Excel file
#'
#' @param up list with each FindMarkers() each
#' @param sg list object from genesorteR
#' @param SeuratObject with IDents() and $seurat_clusters
#' @param OUTDIR output directory path
#' @examples
#' saveMarkers.Excel(up, sg, SeuratObject, OUTDIR)

saveMarkers.Excel <- function(up=NULL, sg=NULL, SeuratObject, OUTDIR) {

	require(openxlsx)
	if(is.null(up) & is.null(sg)) stop("ERROR : either up or sg must be defined");

	## Folder to store data
	MARKERS_OUTDIR <- paste0(OUTDIR,"/Markers")
	if(! dir.exists(MARKERS_OUTDIR)) dir.create(MARKERS_OUTDIR, recursive = TRUE)

	## Excel file
	xlsx_file <- paste0(OUTDIR, "/", Project(SeuratObject), "_markers.xlsx")
	wb <- createWorkbook()

	NclusterNcells <- table(SeuratObject$seurat_clusters)

	IdentClust <- levels(SeuratObject$seurat_clusters)
	IdentSet <- levels(Idents(SeuratObject))

	# 1st sheet
	sheet1 <- data.frame(Cluster=paste0("cluster", IdentClust),
			     Ncells=as.numeric(NclusterNcells),
			     Perc=as.numeric(NclusterNcells)/sum(NclusterNcells)*100,
			     Assignment=rep("?", length(NclusterNcells)),
			     Markers=rep("?", length(NclusterNcells)))
	if(!all(IdentClust==IdentSet)) sheet1$Assignment <- IdentSet;

	addWorksheet(wb, sheetName = "Assignment")
	writeData(wb, sheet ="Assignment", x=sheet1, rowNames=FALSE)

	# GeneSorteR
	if(!is.null(sg)) {
		   specScore <- as.matrix(sg$specScore)
		   colnames(specScore) <- paste0("cluster", colnames(specScore))
		   addWorksheet(wb, sheetName = "genesorteR_specScore")
		   writeData(wb, sheet ="genesorteR_specScore", x=specScore, rowNames=TRUE)

		   condGeneProb <- as.matrix(sg$condGeneProb)
		   colnames(condGeneProb) <- paste0("cluster", colnames(condGeneProb))
		   addWorksheet(wb, sheetName = "genesorteR_condGeneProb")
		   writeData(wb, sheet ="genesorteR_condGeneProb", x=condGeneProb, rowNames=TRUE)
	}
	# Wilcox
	if(!is.null(up)) {
		for(idx in names(up)) {
			addWorksheet(wb, sheetName = paste0("wilcox_cluster",idx))

			firstgofirst <- c("cluster", "gene")
			remaining <- setdiff(colnames(up[[idx]]), firstgofirst)
			col_names <- c(firstgofirst, remaining)
			writeData(wb, sheet = paste0("wilcox_cluster",idx), 
				  x=up[[idx]][,c(col_names)], 
				  rowNames=FALSE)
		}
	}

	saveWorkbook(wb, file =  xlsx_file, overwrite = TRUE)

}

# source: split-violin at https://stackoverflow.com/questions/35717353/split-violin-plot-with-ggplot2
# source: quantiles in split-vilin at https://stackoverflow.com/questions/47651868/split-violin-plot-with-ggplot2-with-quantiles/47652563#47652563
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin,
  draw_group = function(self, data, ..., draw_quantiles = NULL) {
    # Original function by Jan Gleixner (@jan-glx)
    # Adjustments by Wouter van der Bijl (@Axeman)
    data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
    grp <- data[1, "group"]
    newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
    newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
    newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
    if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
      stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
      quantiles <- create_quantile_segment_frame(data, draw_quantiles, split = TRUE, grp = grp)
      aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
      aesthetics$alpha <- rep(1, nrow(quantiles))
      both <- cbind(quantiles, aesthetics)
      quantile_grob <- GeomPath$draw_panel(both, ...)
      ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
    }
    else {
      ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
    }
  }
)

create_quantile_segment_frame <- function(data, draw_quantiles, split = FALSE, grp = NULL) {
  dens <- cumsum(data$density) / sum(data$density)
  ecdf <- stats::approxfun(dens, data$y)
  ys <- ecdf(draw_quantiles)
  violin.xminvs <- (stats::approxfun(data$y, data$xminv))(ys)
  violin.xmaxvs <- (stats::approxfun(data$y, data$xmaxv))(ys)
  violin.xs <- (stats::approxfun(data$y, data$x))(ys)
  if (grp %% 2 == 0) {
    data.frame(
      x = ggplot2:::interleave(violin.xs, violin.xmaxvs),
      y = rep(ys, each = 2), group = rep(ys, each = 2)
    )
  } else {
    data.frame(
      x = ggplot2:::interleave(violin.xminvs, violin.xs),
      y = rep(ys, each = 2), group = rep(ys, each = 2)
    )
  }
}

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, position = position, 
        show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}


plotVln <- function(SeuratObject, gene=NULL, meta, stats, vlnsplit=TRUE, fontTXT="Arial", nCell.y=-03, pt.alpha=1, seed=12345) {
	set.seed(seed)
	# Theme
	theme_min <- theme_cowplot() +
		theme(
		      plot.title = element_text(size=20, hjust=0.5),
		      plot.margin = unit(c(0.3,0.3,0.3,0.3), "cm"),
		      axis.text.x = element_text(size=18, angle = 45, vjust=1, hjust=1, family=fontTXT),
		      axis.text.y = element_text(size=18, family=fontTXT),
		      axis.title = element_text(size=20, family=fontTXT),
		      legend.text = element_text(size=18, family=fontTXT),
		      legend.title = element_blank(),
		      legend.position="bottom",
		      )
	# Color palette
	cond_col <- c("grey30", "red")
	cond_bg <- c("grey60", "red")

	# Tidy data
	if(!is.null(gene)) {
		if(!gene %in% rownames(SeuratObject)) {
			stop(paste0("[ERROR]: Gene does NOT exist in dataset:", gene,"\n"))
	}

		dat <- data.frame(expr=SeuratObject@assays$RNA@data[gene,],
				  cell=Idents(SeuratObject),
				  stim=SeuratObject$stim)
		gtitle <- gene
	} else if(!is.null(meta)) {
		if(!meta %in% colnames(SeuratObject@meta.data)) {
			stop(paste0("[ERROR]: Meta does NOT exist in dataset:", meta,"\n"))
		}
		dat <- data.frame(expr=SeuratObject@meta.data[,meta],
				  cell=Idents(SeuratObject),
				  stim=SeuratObject$stim)
		gtitle <- meta
	}
	# If any stats to add,
	if(!is.null(stats)) {
		stats2 <- apply(stats,1,function(z) {
					# Make string with stats
					txt <- paste0(names(z),"=",z)
					# Rm empty stats if no test
					txt <- gsub("\n?^=", "", txt)
					# Rm NaN/NA stats if no test
					txt <- gsub(".*= +(NA|NaN)","",txt)
					return(txt)
				  })
		stats3 <- apply(stats2,2, function(z) paste(z,collapse="\n"))

		dat$tag <- stats3[dat$cell] 
		# For annotating significance
		dat$maxscore <- max(dat$expr)
		# Only one cell per key2, that is, per group by gene/meta. 
		# So plot is not overwritten per each observation
		dat$tag2 <- NA 
		dat$tag2[!duplicated(dat$cell)] <- dat$tag[!duplicated(dat$cell)]


	}

	# Fill missing groups
	tmp <- table((unique(dat[,c("cell","stim")])$cell))
 	stopifnot(all(tmp==length(unique(dat$stim))))

	give.n <- function(x, y=nCell.y){
		return(c(y = y,label = length(x)) # the actual tag
		) 
	}


	if(vlnsplit) {
	gg <- ggplot(dat,
       		aes(x=cell, y=expr, fill=stim, colour=stim)) +
		   geom_split_violin(width=0.5, scale="width",colour="black", alpha=0.3, draw_quantiles = 0.5) +
 		   geom_jitter(size=1, shape=16, alpha=pt.alpha, position=position_jitterdodge(0.6)) +
		   stat_summary(mapping=aes(colour=stim),
				fun.data = give.n, geom = "text", fun.y = median,
				size=4,
				position = position_dodge(width = 0.8)) +
		   scale_fill_manual(values=cond_bg) +
		   scale_colour_manual(values=cond_col) +
		   xlab("") + ggtitle(gtitle) + 
		   ylab("Expression Level") +
		   theme_min
	} else {
		gg <- ggplot(dat,
			     aes(x=cell, y=expr, fill=stim, colour=stim)) +
			   geom_violin(scale="width", colour="black", alpha=0.3, draw_quantiles = 0.5) +
 		   	geom_jitter(size=1, shape=16, alpha=pt.alpha, position=position_jitterdodge(0.6)) +
			   stat_summary(mapping=aes(colour=stim),
					fun.data = give.n, geom = "text", fun.y = median,
					size=4,
					position = position_dodge(width = 0.8)) +
			   scale_fill_manual(values=cond_bg) +
			   scale_colour_manual(values=cond_col) +
			   xlab("") + ggtitle(gtitle) + 
			   ylab("Expression Level") +
			   theme_min
	}
	if(!is.null(stats)) {
		gg <- gg + geom_text(aes(x=cell, y=maxscore*1.07, label=tag2), 
				     position = position_dodge(width=0),
				     hjust = 0.5, size=4, colour="black") +
			   coord_cartesian(clip="off")
	}
return(gg)

}
