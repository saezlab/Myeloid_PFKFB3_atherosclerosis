#This is a function designed to compute progeny scores from log2foldchanges.
#Copyright (C) 2017  Aurelien Dugourd

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

#'\code{progenyScores}
#'
#'This function compute progeny pathway scores as a series of weighted sum.
#'For each pathway a score is computed as the weighted sum of the sample/contrast statistic using the progeny coeficients as weights.
#'The gene identifiers between the measurments and the progeny coeficient matrix should be coeherent.
#'
#'@param df a n*m data frame, where n is the number of omic features (genes). m isn't really important, as long as at least one column corespond to a sample or contrast statistic. One of the column should correspond to the gene symboles.
#'@param cm a progeny coeficient matrix. One of the column should be the gene symboles.
#'@param dfIndex an integer corresponding to the column number of the gene identifiers of df.
#'@param FCIndex an integer corresponding to the column number that contains the statistic to be consdered to compute the progeny scores.
#'@param cmIndex an integer corresponding to the column number of the gene identifiers of the weight matrix.
#'@return a named vector where each element is a scores and the names are the corresponding pathways.
progenyScores <- function(df, cm, dfIndex = 1, FCIndex = 3, cmIndex = 1) {
  names(df)[dfIndex] <- "X1"
  names(cm)[cmIndex] <- "X1"
  df <- df[complete.cases(df[,FCIndex]),]

  merged <- merge(df[,c(dfIndex,FCIndex)],cm)

  for (pathway in names(cm[,-cmIndex]))
  {
    merged[,pathway] <- merged[,2]*merged[,pathway]
  }

  progeny_scores <- colSums(merged[,c(3:length(merged[1,]))])
  names(progeny_scores) <- names(merged[,c(3:length(merged[1,]))])

  return(progeny_scores)
}

######################################################
######################################################

library(readr)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(ggrepel)
library(parallel)
library(cowplot)

#'\code{runProgeny}
#'
#'This function is designed to compute progeny pathway scores and assess there significance using a gene sampling based permutation strategy, for a series of experimental samples/contrasts.
#'
#'@param df A data.frame of n*m+1 dimension, where n is the number of omic features to be considered and m is the number of samples/contrasts.
#'The first column should be the identifiers of the omic features. These identifiers must be coherent with the identifers of the weight matrix.
#'@param weight_matrix A progeny coeficient matrix. the first column should be the identifiers of the omic features, and should be coherent with the identifiers provided in df.
#'@param k The number of permutations to be preformed to generate the null-distribution used to estimate significance of progeny scores. Default value is 10000.
#'@param nCores The number of cores that should be used for the parallelised generation of null distribution. If nothing is specified, the bumber of cores will be the minimum value between 1000 and the number of available cores minus 1.
#'@return This function returns a list of two elements. The first element is a dataframe of p*m+1 dimensions, where p is the number of progeny pathways, and m is the number of samples/contrasts.
#'Each cell represent the significance of a progeny pathway score for one sample/contrast. The signifcance ranges between -1 and 1. The significance is equal to x*2-1, x being the quantile of the progeny pathway score with respect to the null distribution.
#'Thus, this significance can be interpreted as the equivalent of 1-p.value (two sided test over an empirical distribution) with the sign indicating the direction of the regulation.
#'The sceond element is the null distribution list (a null distribution is generated for each sample/contrast).
runProgeny <- function(df,weight_matrix,k = 10000, nCores = 1000, z_score = T)
{
  scores <- list(0)
  for (i in 2:length(df[1,]))
  {
    scores[[i-1]] <- progenyScores(df[,c(1,i)],weight_matrix,dfIndex = 1,FCIndex = 2,cmIndex = 1)
  }
  scores <- as.data.frame(data.frame(scores))

  null_dist_list <- list(0)

  for (j in 2:length(df[1,]))
  {
    print(paste(j-1," ",sep=""))
    FC_sample <- df[,c(1,j)]
    names(FC_sample) <- c("ID","values")

    null_dist <- list(0)

    null_dist <- mclapply(1:k, function(i)
    {
      if (!(i %% 200))
      {
        print(paste(j-1,i,sep=" "))
      }
      FC_sample$values <- sample(FC_sample$values,size = length(FC_sample[,2]),replace = F)
      null_dist[[i]] <- progenyScores(FC_sample,weight_matrix,dfIndex = 1,FCIndex = 2,cmIndex = 1)
    }, mc.cores = min(c(nCores,detectCores()-1)))

    null_dist <- data.frame(null_dist)
    null_dist <- as.data.frame(null_dist)
    #print(head(null_dist))
    names(null_dist) <- c(1:k)
    null_dist_list[[j-1]] <- null_dist
  }

  score_probas <- scores

  if (z_score)
  {
    for (i in 1:length(scores[1,]))
    {
      null_dist <- null_dist_list[[i]]
      #View(as.data.frame(null_dist))
      for (j in 1:length(scores[,1]))
      {
        #percentile <- ecdf(null_dist[j,])
        #print(score_probas[j,i])
        #print(mean(null_dist[j, ]))
        #print(sd(null_dist[j, ]))
        score_probas[j,i] <- (score_probas[j,i] - mean(as.numeric(null_dist[j, ]))) / sd(as.numeric(null_dist[j, ]))
      }
    }

    score_probas$pathway <- row.names(score_probas)
    score_probas <- score_probas[,c(length(score_probas[1,]),1:(length(score_probas[1,])-1))]

    return(list(score_probas,null_dist_list))
  }
  else
  {
    for (i in 1:length(scores[1,]))
    {
      null_dist <- null_dist_list[[i]]
      for (j in 1:length(scores[,1]))
      {
        percentile <- ecdf(null_dist[j,])
        score_probas[j,i] <- percentile(scores[j,i])
      }
    }

    score_probas <- score_probas*2-1

    score_probas$pathway <- row.names(score_probas)
    score_probas <- score_probas[,c(length(score_probas[1,]),1:(length(score_probas[1,])-1))]

    return(list(score_probas,null_dist_list))
  }
}

#'\code{progenyScatter}
#'
#'This function generate a series of scatter plot with marginal distribution (in the form of an arrangeGrob object), for each progeny pathway and sample/contrast.
#'Each scatter plot has progeny weights as x axis and the gene level stat used to compute progeny score as y axis.
#'The marginal distribution of the gene level stats is displayed on the right of the plot to give a visual support of the significnce of each gene contributing to the progeny pathway score.
#'The colors green and red represent respectivelly positive and negative contribution of genes to the progeny pathway.
#'For each gene contribution, 4 cases are possible, as the combinaisons of the sign of the gene level stat and the sign of the gene level weigth.
#'Positive weight will lead to a positive(green)/negative(red) gene contribution if the gene level stat is positive/negative.
#'Negative weigth will lead to a negative(red)/positive(green) gene contribution if the gene level stat is positive/negative.
#'
#'@param df a n*m data frame, where n is the number of omic features (genes). m isn't really important, as long as at least one column corespond to a sample or contrast statistic. One of the column should correspond to the gene symboles.
#'@param cm a progeny coeficient matrix. One of the column should be the gene symboles.
#'@param dfID an integer corresponding to the column number of the gene identifiers of df.
#'@param weightID an integer corresponding to the column number of the gene identifiers of the weight matrix.
#'@param statname The neame of the stat used, to be displayed on the plot
#'@return The function returns a list of list of arrangeGrob object.The first level list elements correspond to samples/contrasts. The second level correspond to pathways.
#'The plots can be saved in a pdf format using the saveProgenyPlots function.
progenyScatter <- function(df,weight_matrix,dfID = 1, weightID = 1, statName = "gene stats", fontfamily="Arial")
{
  plot_list_contrasts <- list(0)
  for (i in 2:length(df[1,]))
  {
    plot_list_pathways <- list(0)
    for (j in 2:length(weight_matrix[1,]))
    {
      sub_df <- df[,c(dfID,i)]



      pathway_weights <- weight_matrix[,c(weightID,j)]
      names(sub_df) <- c("ID","stat")

      minstat <- min(sub_df$stat)
      maxstat <- max(sub_df$stat)
      histo <- ggplot(sub_df, aes(x = stat, fill = "blue")) + geom_density() + 
	      coord_flip() + scale_fill_manual( values = c("#00c5ff")) + xlim(minstat, maxstat) + 
	      theme_minimal() + 
	      theme(legend.position = "none", 
		    text=element_text(family=fontfamily, size=16),
		    axis.text.x = element_blank(), 
		    axis.ticks.x = element_blank(), 
		    axis.title.y = element_blank(), 
		    axis.text.y = element_blank(), 
		    axis.ticks.y = element_blank(), 
		    plot.background = element_rect(colour = "white"),
		    panel.grid.major = element_blank(), 
		    panel.grid.minor = element_blank())

      names(pathway_weights) <- c("ID","weight")
      pathway_weights <- pathway_weights[pathway_weights$weight != 0,]

      percentile <- ecdf(sub_df$stat)

      sub_df <- merge(sub_df,pathway_weights,by = "ID")

      sub_df$color <- "3"
      sub_df[(sub_df$weight > 0 & sub_df$stat > 0),"color"] <- "1"
      sub_df[(sub_df$weight > 0 & sub_df$stat < 0),"color"] <- "2"
      sub_df[(sub_df$weight < 0 & sub_df$stat > 0),"color"] <- "2"
      sub_df[(sub_df$weight < 0 & sub_df$stat < 0),"color"] <- "1"

      sub_df[(percentile(sub_df$stat) < .95 & percentile(sub_df$stat) > .05),1] <- NA

      print(paste("weights of ",names(weight_matrix)[j], sep = ""))

      title <- paste("weights of ",names(weight_matrix)[j], sep = "")

      scatterplot <- ggplot(sub_df, aes(x = weight, y = stat, color = color)) + geom_point() +
        # scale_colour_manual(values = c("#15ff00","#ff0000","#c9c9c9")) + #green and red
        scale_colour_manual(values = c("red","royalblue3","grey")) +
        geom_label_repel(aes(label = ID), family=fontfamily, size=5) +
        ylim(minstat, maxstat) +  
	theme_minimal() + cowplot::theme_cowplot() +  
	theme(legend.position = "none", plot.background = element_rect(colour = "white", fill="white"), 
	      text=element_text(family=fontfamily, size=16, color="black"),
	      axis.text=element_text(family=fontfamily, size=16, color="black")) + 
	geom_vline(xintercept = 0, linetype = 'dotted') + geom_hline(yintercept = 0, linetype = 'dotted') + labs(x = title, y = statName)

      lay <- t(as.matrix(c(1,1,1,1,2)))
      gg <- arrangeGrob(scatterplot, histo, nrow = 1, ncol = 2, layout_matrix = lay)

      #grid.arrange(gg)
      plot_list_pathways[[j-1]] <- gg
    }
    names(plot_list_pathways) <- names(weight_matrix[,-weightID])
    plot_list_contrasts[[i-1]] <- plot_list_pathways
  }
  return(plot_list_contrasts)
}

#######################################################################################################################################################################################################################################################
#######################################################################################################################################################################################################################################################
#######################################################################################################################################################################################################################################################

#'\code{saveProgenyPlots}
#'
#'This function is designed to save the plots (in pdf format) of a nested (2 level) list of arrangeGrob objects, such as the one returned by the progenyScatter function.
#'
#'@param plots a list of list of arrangeGrob object (such as the one returned by the progenyScatter function.).The first level list elements correspond to samples/contrasts. The second level correspond to pathways.
#'The plots can be saved in a pdf format using the saveProgenyPlots function.
#'@param contrast_names a vector of same length as the first level of the plot list corresponding to the names of each sample/contrast
#'@param dirpath the path to the directory where the plotsshould be saved
saveProgenyPlots <- function(plots, contrast_names, dirpath)
{
  i <- 1
  for (condition in plots)
  {
    dirname <- paste(dirpath,contrast_names[i], sep = "")
    dir.create(dirname, recursive = T, showWarnings = F)
    j <- 1
    for (pathway in condition)
    {
      filename <- paste(dirname,names(condition)[j],sep = "/")
      filename <- paste(filename,".pdf",sep = "")
      print(filename)
      ggsave(filename, pathway,device = "pdf", dpi = 300)
      j <- j+1
    }
    i <- i+1
  }
}

#######################################################################################################################################################################################################################################################
#######################################################################################################################################################################################################################################################
#######################################################################################################################################################################################################################################################

library(UniProt.ws)

progenyPathways <- function(tableTop_list,model_matrix, scores,omnipath, path, mapping = FALSE, all_nodes = FALSE)
{
  dir.create(path = path, recursive = T, showWarnings = F)
  if(mapping)
  {
    Up <- UniProt.ws(9606)
    columns(Up)
    kt <- "UNIPROTKB"
    keys <- unique(c(omnipath$source,omnipath$target))
    columns <- c("GENES")

    res <- select(Up, keys, columns, kt)

    res$GENES <- gsub("[ ].*","",res$GENES)
    names(res) <- c("source","source_gene")
    omnipath <- merge(omnipath,res, by = "source")
    names(res) <- c("target","target_gene")
    omnipath <- merge(omnipath,res, by = "target")
  } else
  {
    names(omnipath) <- gsub("source","source_gene",names(omnipath))
    names(omnipath) <- gsub("target","target_gene",names(omnipath))
  }

  #print(head(omnipath))
  names(model_matrix)[1] <- "ID"
  model_matrix$ID <- toupper(model_matrix$ID)

  scores_df <- scores[[1]]
  conditions <- names(scores_df)
  condition_list <- list(0)
  for (i in 2:length(conditions))
  {
    tableTop <- tableTop_list[[i-1]]
    names(tableTop)[1] <- "ID"
    pathway_list <- list(0)

    #print(length(row.names(scores_df)))
    for (j in 2:(length(row.names(scores_df))+1))
    {
      print(paste((j-1),row.names(scores_df)[j-1], sep = " "))
      network_and_attributes <- list(0)
      model_matrix_pathway <- model_matrix[model_matrix[,j] != 0,c(1,j)]
      #print(head(model_matrix_pathway))
      omnipath_pathway <- omnipath[omnipath$source_gene %in% model_matrix_pathway$ID & omnipath$target_gene %in% model_matrix_pathway$ID,]
      if (all_nodes == TRUE)
      {
        isolated_nodes <- as.data.frame(matrix(nrow = length(model_matrix_pathway[,1]), ncol = 6))
        isolated_nodes[,1] <- model_matrix_pathway[,1]
        names(isolated_nodes) <- names(omnipath_pathway)
        isolated_nodes <- isolated_nodes[!(isolated_nodes[,1] %in% omnipath_pathway[,1]) | !(isolated_nodes[,1] %in% omnipath_pathway[,2]),]
        omnipath_pathway <- as.data.frame(rbind(omnipath_pathway,isolated_nodes))
      }
      model_matrix_pathway <- model_matrix_pathway[model_matrix_pathway$ID %in% omnipath_pathway$source_gene | model_matrix_pathway$ID %in% omnipath_pathway$target_gene,]
      #print(head(model_matrix_pathway))
      tableTop$ID <- toupper(tableTop$ID)

      nodes <- merge(model_matrix_pathway,tableTop, by = "ID", all.x = TRUE)
      #print(head(nodes))
      nodes$adjusted_t <- nodes$t*nodes[,2]

      omnipath_pathway$sign <- omnipath_pathway$is_stimulation - omnipath_pathway$is_inhibition
      condname <- conditions[i]
      pathname <- row.names(scores_df)[j-1]
      network_name <- paste(path,paste(paste(condname,pathname, sep="_"),"network.csv",sep = "_"),sep = "")
      nodes_name <- paste(path,paste(paste(condname,pathname, sep="_"),"nodes.csv",sep = "_"),sep = "")
      write_csv(omnipath_pathway,network_name)
      write_csv(nodes,nodes_name)
      network_and_attributes[[1]] <- omnipath_pathway
      network_and_attributes[[2]] <- nodes
      pathway_list[[j-1]] <- network_and_attributes
    }
    names(pathway_list) <- names(model_matrix[,c(2:length(model_matrix[1,]))])
    condition_list[[i-1]] <- pathway_list
  }
  names(condition_list) <- names(scores_df)[c(2:length(names(scores_df)))]
  return(condition_list)
}

####################################

#'\code{plotProgenyPathway}
#'
#'In progress
#'
#'@param progenyPathway temp
#'@param filename temp
plotProgenyPathway <- function(progenyPathway, filename)
{
  if ( length(progenyPathway[[1]][,1] > 0))
  {
    names(progenyPathway[[2]])[2] <- "coeficients"

    network <- graph_from_data_frame(d=progenyPathway[[1]], vertices=progenyPathway[[2]], directed=T)

    myLines <- c(1,2,1)

    myEdgeColors <- c("red","grey","green")

    myColors <- unique(createLinearColors(V(network)$adjusted_t), 100)

    V(network)$sizes <- symmetricalScaling(abs(V(network)$coeficients), 15) + 5

    print(paste("saving netowork to :",filename,sep = ""))
    pdf(filename)
    plot.igraph(network, layout = layout_nicely, vertex.size = V(network)$sizes, vertex.color = myColors[symmetricalScaling(V(network)$adjusted_t, 100)], rescale = TRUE,
                edge.arrow.size=0.3,
                edge.curved=0.4,
                edge.lty = myLines[E(network)$sign+2],
                edge.color = myEdgeColors[E(network)$sign+2],
                vertex.label.color = "black",
                vertex.label.cex = 0.65)
    dev.off()
  }
  else
  {
    print("Network is empty... Skip it.")
  }
}

########################################

#'\code{plotProgenyPathwayResults}
#'
#'In progress
#'
#'@param progenyPathwaysResult temp
#'@param outpath temp
#'@return temp
plotProgenyPathwayResults <- function(progenyPathwaysResult, outpath)
{
  if (length(progenyPathwaysResult) > 1)
  {
    i <- 1
    for (condition in progenyPathwaysResult)
    {
      outDirName <- paste(outpath,paste("/",names(progenyPathwaysResult)[i], sep = ""), sep = "")
      dir.create(outDirName, recursive = T, showWarnings = F)
      j <- 1
      for (pathway in condition)
      {
        pathwayName <- names(pathway)[j]
        print(pathwayName)
        filename <- paste(outDirName,paste("/",paste(pathwayName,".pdf", sep = ""), sep = ""), sep = "")
        plotProgenyPathway(pathway,filename)

        j <- j+1
      }

      i <- i+1
    }
  }
  else
  {
    outDirName <- outpath
    dir.create(outDirName, recursive = T, showWarnings = F)
    condition <- progenyPathwaysResult[[1]]
    j <- 1
    for (pathway in condition)
    {
      pathwayName <- names(condition)[j]
      print(pathwayName)
      filename <- paste(outDirName,paste("/",paste(pathwayName,".pdf", sep = ""), sep = ""), sep = "")
      plotProgenyPathway(pathway,filename)

      j <- j+1
    }
  }
}


####################################

#'\code{progenyToTarget}
#'
#'This function is designed to generate a list of potential perturbation targets based on the significance of progeny scores. This list can give insights into what could be the starting point of perturbed signaling pathways.
#'The function simply isolate the perturbed pathway based on a threshold of pathway significance and genes associated with the pathway.
#'
#'@param progeny_scores a dataframe of dimension n*m+1 where n is the number of progeny pathways and m is the number of sample/contrast, where the first column is the pathway identifier.
#'@param progeny_mapping a mapping table that connect the progeny pathway names to the expectd targets.
#'@param threshold this threshold indicates which pathway should be considered as effectivelly deregulated. It is important to be aware that the significance scores are not corrected for multiple testing, hence the threshold should be decided with this in mind.
#'@return this function returns a list of dataframe (one for each sample/contrast). Each dataframe contain the suspected perturbation targets and the sign of the deregulation based on the progeny score.
progenyToTarget <- function(progeny_scores, progeny_mapping, threshold = 0.01)
{
  targets_list <- list(0)
  names(progeny_mapping) <- c("progeny","symbol","uniprot")
  for (i in 2:length(progeny_scores[1,]))
  {
    progeny_condition <- progeny_scores[,c(1,i)]
    progeny_condition <- progeny_condition[abs(progeny_condition[,2]) > 1-threshold,]
    progeny_condition[,2] <- ifelse(progeny_condition[,2] > 0, 1, -1)
    names(progeny_condition) <- c("progeny","score")
    progeny_condition <- merge(progeny_mapping, progeny_condition, by = "progeny")
    targets_list[[i-1]] <- progeny_condition
  }
  names(targets_list) <- names(progeny_scores[2:length(progeny_scores[1,])])
  return(targets_list)
}

######################################

#'\code{writeProgenyTargets}
#'
#'This function simply write as csv files the dataframes in the target list returned by the progenyToTarget function.
#'
#'@param targets_list a list of progeny target dataframe, such as the one returned by the progenyToTarget function.
#'@param outpath The path to the directory where the progeny target dataframes should be written.
writeProgenyTargets <- function(targets_list, outpath)
{
  i <- 1
  contrast_names <- names(targets_list)
  for(targets in targets_list)
  {
    dir.create(outpath, recursive = T, showWarnings = F)
    filename <- paste(contrast_names[i],"_targets.csv",sep = "")
    print(filename)
    filename <- paste("/",filename, sep = "")
    write_csv(targets,paste(outpath,filename, sep = ""))
    i <- i+1
  }
}

#######################################

#'\code{runProgenyFast}
#'
#'This function is designed to compute progeny pathway scores and assess there significance using a gene sampling based permutation strategy, for a series of experimental samples/contrasts.
#'
#'@param df A data.frame of n*m+1 dimension, where n is the number of omic features to be considered and m is the number of samples/contrasts.
#'The first column should be the identifiers of the omic features. These identifiers must be coherent with the identifers of the weight matrix.
#'@param weight_matrix A progeny coeficient matrix. the first column should be the identifiers of the omic features, and should be coherent with the identifiers provided in df.
#'@param k The number of permutations to be preformed to generate the null-distribution used to estimate significance of progeny scores. Default value is 10000.
#'@param z_score if true, the z-scores will be returned for the pathway activity estimations. Else, the function returns a normalised z-score value between -1 and 1.
#'@param get_nulldist if true, the null score distribution used for normalisation will be returned along with the actual normalised score data frame.
#'@return This function returns a list of two elements. The first element is a dataframe of p*m+1 dimensions, where p is the number of progeny pathways, and m is the number of samples/contrasts.
#'Each cell represent the significance of a progeny pathway score for one sample/contrast. The signifcance ranges between -1 and 1. The significance is equal to x*2-1, x being the quantile of the progeny pathway score with respect to the null distribution.
#'Thus, this significance can be interpreted as the equivalent of 1-p.value (two sided test over an empirical distribution) with the sign indicating the direction of the regulation.
#'The sceond element is the null distribution list (a null distribution is generated for each sample/contrast).
runProgenyFast <- function(df,weight_matrix,k = 10000, z_scores = T, get_nulldist = F)
{
  resList <- list()
  if(get_nulldist)
  {
    nullDist_list <- list()
  }

  for(i in 2:length(df[1,]))
  {
    current_df <- df[,c(1,i)]
    current_df <- current_df[complete.cases(current_df),]
    t_values <- current_df[,2]

    current_weights <- weight_matrix

    names(current_df)[1] <- "ID"
    names(current_weights)[1] <- "ID"

    common_ids <- merge(current_df, current_weights, by = "ID")
    common_ids <- common_ids$ID
    common_ids <- as.character(common_ids)

    row.names(current_df) <- current_df$ID
    current_df <- as.data.frame(current_df[common_ids,-1])

    row.names(current_weights) <- current_weights$ID
    current_weights <- as.data.frame(current_weights[common_ids,-1])

    current_mat <- as.matrix(current_df)
    current_weights <- t(current_weights)

    scores <- as.data.frame(current_weights %*% current_mat)

    null_dist_t <- replicate(k, sample(t_values,length(current_mat[,1]), replace = F))

    null_dist_scores <- current_weights %*% null_dist_t

    if(get_nulldist)
    {
      nullDist_list[[i-1]] <- null_dist_scores
    }

    if(z_scores)
    {
      scores$mean <- apply(null_dist_scores,1,mean)
      scores$sd <- apply(null_dist_scores,1,sd)
      resListCurrent <- (scores[,1]-scores[,2])/scores[,3]
      names(resListCurrent) <- names(weight_matrix[,-1])
      resList[[i-1]] <- resListCurrent
    }
    else
    {
      for(j in 1:length(weight_matrix[,-1]))
      {
        ecdf_function <- ecdf(null_dist_scores[j,])
        scores[j,1] <- ecdf_function(scores[j,1])
      }
      score_probas <- scores*2-1

      resListCurrent <- score_probas[,1]
      names(resListCurrent) <- names(weight_matrix[,-1])
      resList[[i-1]] <- resListCurrent
    }
  }
  names(resList) <- names(df[,-1])
  resDf <- as.data.frame(resList)
  if(get_nulldist)
  {
    names(nullDist_list) <- names(df[,-1])
    return(list(resDf, nullDist_list))
  }
  else
  {
    return(resDf)
  }

}

######################

#'\code{heatProgeny}
#'
#'This function is designed to generate two separate heatmap of the top 10 positive and top 10 negative responsive genes of a progeny pathway, displaying the measured expression of those genes.
#'
#'@param t_table A data.frame of n*m+1 dimension, where n is the number of omic features to be considered and m is the number of samples/contrasts.
#'The first column should be the identifiers of the omic features. These identifiers must be coherent with the identifers of the weight matrix.
#'@param weight_matrix A progeny coeficient matrix. the first column should be the identifiers of the omic features, and should be coherent with the identifiers provided in df.
#'@param pathway The name of the progeny pathway to generate heatmap for. The name should be same of those of the weight matrix
#'@param out_dir the output directory to save the heatmaps as pdf
#'@param height manual option for determining the output file height in inches.
#'@param width manual option for determining the output file width in inches.
heatProgeny <- function(t_table, model_mat, pathway, out_dir, height = 5, width = 5)
{
  row.names(t_table) <- t_table[,1]
  model_mat <- model_mat[order(model_mat[,pathway], decreasing = T),]
  model_mat <- model_mat[model_mat[,1] %in% row.names(t_table),]

  responsive_gene_up <- model_mat[model_mat[,pathway] > 0,1]
  if(length(responsive_gene_up) > 0)
  {
    responsive_gene_up <- responsive_gene_up[1:10]
  }

  model_mat <- model_mat[order(model_mat[,pathway], decreasing = F),]
  responsive_gene_down <- model_mat[model_mat[,pathway] < 0,1]
  if(length(responsive_gene_down) > 0)
  {
    responsive_gene_down <- responsive_gene_down[1:10]
  }

  t_table_up <- t_table[t_table[,1] %in% responsive_gene_up,-1]
  t_table_down <- t_table[t_table[,1] %in% responsive_gene_down,-1]
  t <- as.vector(t(t_table_up))

  if(min(t) < 0 & max(t) > 0)
  {
    row.names(t_table_up) <- toupper(row.names(t_table_up))
    spread <- max(t) - min(t)
    lower_part <- abs(round((min(t)/spread) * 100))
    upper_part <- round((1 - abs((min(t)/spread))) * 100)
    palette1 <- createLinearColors(t[t < 0],withZero = F , maximum = lower_part)
    palette2 <- createLinearColors(t[t > 0],withZero = F , maximum = upper_part)
    palette <- c(palette1,palette2)
    file_name_up <- paste(out_dir, paste("progeny_", paste(pathway, "_up.pdf", sep= ""), sep = ""), sep = "/")
    pheatmap(t_table_up, display_numbers = T, cluster_cols = F, color = palette, filename = file_name_up, height = height, width = width, cluster_rows = F)
  }

  t <- as.vector(t(t_table_down))

  if(min(t) < 0 & max(t) > 0)
  {
    row.names(t_table_down) <- toupper(row.names(t_table_down))
    spread <- max(t) - min(t)
    lower_part <- abs(round((min(t)/spread) * 100))
    upper_part <- round((1 - abs((min(t)/spread))) * 100)
    palette1 <- createLinearColors(t[t < 0],withZero = F , maximum = lower_part)
    palette2 <- createLinearColors(t[t > 0],withZero = F , maximum = upper_part)
    palette <- c(palette1,palette2)
    file_name_down <- paste(out_dir, paste("progeny_", paste(pathway, "_down.pdf", sep= ""), sep = ""), sep = "/")
    pheatmap(t_table_down,display_numbers = T, cluster_cols = F, color = palette, filename = file_name_down, height = height, width = width, cluster_rows = F)
  }
}
