
#' Perform a two-sample Wilcoxon test of vectors of data and 
#' retrieve a tidy data.frame of statistics (W, pvalue, effect size).
#' 
#' @param xy a union of the two numeric vector samples as it would be split by gr
#' @param gr a factor vector to stratify samples from the two groups from xy.
wilcox.test_stats <- function(xy, gr, min.size=5, seed=101010) {

	stopifnot(length(xy)==length(gr))
	stopifnot(is.numeric(xy))
	stopifnot(is.factor(gr))
	stopifnot(length(levels(gr))==2)

	if(all(table(gr)>=min.size)) {
		   #NOTE: Non-parametric test. We keep same seed in both, which use stats::wilcox.test for stats.
		   set.seed(seed)
		   wout <- stats::wilcox.test(xy ~ gr)

		   #NOTE: We change Group1 so positive r effect sizes are related to test group
		   gr <- relevel(gr, ref=levels(gr)[2])
		   set.seed(seed)
		   weff <- rcompanion::wilcoxonR(x=xy, g=gr)


		   res <- setNames(c(wout$statistic, wout$p.value, weff),
				   c("W", "pvalue", "r"))

	} else {
		cat("[WARN] : Not enough sample size for wilcox test\n", file=stdout())
	res <- c("W"=NA, "pvalue"=NA, "r"=NA)
	}

	return(res)

}

#' Tag significance of a p-value or adjusted p-value
#'
#' @param p a numeric vector of p-values
tagSignif <- function(p) {
	stopifnot(is.numeric(p))
	if(any(na.omit(p)>1 | na.omit(p)<0)) stop("ERROR: expected p-values in a range 0,1");
	tags <- sapply(p, function(pval) {
			       if(!is.na(pval)) {
				       if(pval< 0.001) {
					       txt <- "***"
				       } else if (pval < 0.01) {
					       txt <- "**"
				       } else if (pval < 0.05) {
					       txt <- "*"
				       } else {
					       txt <- "ns"
				       }
			       } else {
					txt <- "nt"
			       }
			       return(txt)
	})

	return(tags)
}
