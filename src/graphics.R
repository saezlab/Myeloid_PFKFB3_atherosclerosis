### TITLE : Prepare graphics for the project
### AUTHOR : Javier Perales-Paton, javier.perales@bioquant.uni-heidelberg.de


### FONT TEXT
if(!require(extrafont)) {
	install.packages("extrafont")
}

# We are going to use Calibri as first option for font text in graphics,
# we also provide an open-source alternative to Calibri, which is Carlito.

fontTXT <- NULL
if(any(grepl("Calibri", extrafont::fonts()))) {
	fontTXT <- "Calibri"
} else if(any(grepl("Carlito", extrafont::fonts()))) {
	fontTXT <- "Carlito"
} else {
	warning("The required font to reproduce the plots is not installed. Please root Readme.md for instructions")
}

if(is.null(fontTXT)){
	warning("Setting 'Arial' as font text by default.")
}

### PAR GRAPHICS

## VlnPlot
VlnPlot.stim <- function(S, meta.feature, ylabTXT="") {
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
