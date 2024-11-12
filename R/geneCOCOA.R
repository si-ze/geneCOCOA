#' @rawNamespace import(data.table, except = c(last, first, between, transpose))
#' @rawNamespace import(psych, except=c(alpha, rescale, "%+%"))
#' @rawNamespace import(scales, except=c(discard, alpha))
#' @import dplyr purrr forcats
#' @import ggplot2 ggrepel ggpp stringr
#' @import msigdbr gemma.R
NULL




#
#' Main function of Gene-COCOA: performs statistical analysis.
#'
#' @param geneset_collection A list of lists. Each list holds gene symbols of a specific gene set. Can be generated with \code{\link{get_msigdb_genesets}}
#' @param GOI A string specifying the GOI gene symbol.
#' @param GOI_expr The GOI expression per sample. Should be generated with get_expr_info(expr, GOI)$GOI_expr.
#' @param expr_df Complete expression data. Should be generated with get_expr_info(expr, GOI)$expr_df.
#' @param samplesize How many predictors should be in each regression model? Should approximate (number of samples)/ 10.
#' @param nsims Bootstrapping: how many simulations should be run? At least 1k recommended.
#'
#' @return A list of data frames holding the results of the Gene-COCOA analysis.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # let CAD_disease: data frame holding expression data in coronary artery disease with rows=genes, columns=samples
#'
#' expr_info <- get_expr_info(expr=CAD_disease, GOI="PDGFD")
#' res <- get_stats(geneset_collection=get_msigigdb_genesets("HALLMARK"), GOI="PDGFD", GOI_expr=expr_info$GOI_expr, expr_df=expr_info$expr_df)
#'}
get_stats <- function(geneset_collection,
                      GOI,
                      GOI_expr,
                      expr_df,
                      samplesize=ncol((expr_df)/10),
                      nsims=1000,
                      verbose=TRUE) {

  message("STARTED calculating stats ", format(Sys.time(), "%H:%M:%S"))


  p_value_list <- list()
  low_power.p_value_list <- list()

  all_RMSEs.real <- list()
  all_RMSEs.random <- list()
  geom_mean_cor.real <- list()
  geom_mean_cor.random <- list()

  # DELETE_intercepts
  # intercepts <- list()
  # DELETE_mean_cor
  # mean_cors <- list()

  geom_mean_expr <- list()

  if (verbose) {
    cat(paste("time", "gene set", "set size", sep="\t"), "\n")
  }

  for (i in 1:length(geneset_collection)) {


    my_genes_in_set <- as.vector(unlist(geneset_collection[i]))
    my_genes_in_set <- gsub("-", "_", my_genes_in_set)
    tpm.genes_in_set <- expr_df[my_genes_in_set,]
    t.tpm.genes_in_set <- as.data.frame(t(tpm.genes_in_set))
    rownames(t.tpm.genes_in_set) <- colnames(tpm.genes_in_set)
    colnames(t.tpm.genes_in_set) <- rownames(tpm.genes_in_set)
    names(t.tpm.genes_in_set) <- rownames(tpm.genes_in_set)
    t.tpm.genes_in_set <- t.tpm.genes_in_set %>% purrr::discard(~all(is.na(.)))
    t.tpm.genes_in_set$GOI <- NULL


    my_set_name <- names(geneset_collection[i]) # get name of gene set
    my_set_size <- ncol(t.tpm.genes_in_set)



    # DELETE_intercepts
    # intercepts[[my_set_name]] <- get_intercept(GOI, GOI_expr, t.tpm.genes_in_set)
    # DELETE_mean_cor
    # mean_cors[[my_set_name]] <- get_mean_cor(GOI, GOI_expr, t.tpm.genes_in_set)
    geom_mean_cor.real[[my_set_name]] <- get_geom_mean_cor(GOI, GOI_expr, t.tpm.genes_in_set)

    tpm.my_rest_genes <- expr_df[!(rownames(expr_df) %in% my_genes_in_set),]
    t.tpm.my_rest_genes <- as.data.frame(t(tpm.my_rest_genes))
    rownames(t.tpm.my_rest_genes) <- colnames(tpm.my_rest_genes)
    colnames(t.tpm.my_rest_genes) <- rownames(tpm.my_rest_genes)
    names(t.tpm.my_rest_genes) <- rownames(tpm.my_rest_genes)
    t.tpm.my_rest_genes <- t.tpm.my_rest_genes %>% purrr::discard(~all(is.na(.)))
    t.tpm.my_rest_genes$GOI <- NULL

    if (verbose) {
      cat(paste(format(Sys.time(), "%H:%M:%S"), my_set_name, my_set_size, sep="\t"), "\n")
    }

    if (my_set_size>=samplesize) { # so that we don't sample from empty gene sets
    if (enough_power(samplesize = samplesize, my_set_size = my_set_size, nsims=nsims)==TRUE) {
      real_RMSEs <- compute_x_RMSEs(GOI,
                                    GOI_expr,
                                    my_t.tpms=t.tpm.genes_in_set,
                                    samplesize=samplesize, nsims=nsims)$RMSEs
      my_rand <- compute_x_RMSEs(GOI,
                                 GOI_expr,
                                 my_t.tpms=t.tpm.my_rest_genes,
                                 samplesize=samplesize, nsims=nsims, geom_mean_cors = TRUE)
      random_RMSEs <- my_rand$RMSEs
      geom_mean_cors <- my_rand$geom_mean_cors

      # t-test testing that the random model RMSE is greater than the real model RMSE
      my_p <- t.test(unlist(real_RMSEs),
                     unlist(random_RMSEs),
                     paired=FALSE, alternative="less")
      p_value_list[[my_set_name]] <- my_p$p.value

      all_RMSEs.real[[my_set_name]] <- real_RMSEs
      all_RMSEs.random[[my_set_name]] <- random_RMSEs

    }
    else {
      real_RMSEs <- compute_all_RMSEs(GOI,
                                    GOI_expr,
                                    my_t.tpms=t.tpm.genes_in_set,
                                    samplesize=samplesize)$RMSEs

      my_rand <- compute_x_RMSEs(GOI,
                                 GOI_expr,
                                 my_t.tpms=t.tpm.my_rest_genes,
                                 samplesize=samplesize,
                                 nsims=enough_power(samplesize = samplesize, my_set_size = my_set_size, nsims=nsims),
                                 geom_mean_cors = TRUE)
      random_RMSEs <- my_rand$RMSEs
      geom_mean_cors <- my_rand$geom_mean_cors


      low_power.p <- t.test(unlist(real_RMSEs),
                            unlist(random_RMSEs),
                            paired=FALSE, alternative="less")
      low_power.p_value_list[[paste(my_set_name,
                                    "_(", enough_power(samplesize = samplesize, my_set_size = my_set_size, nsims=nsims), ")",
                                    sep="")]] <- low_power.p$p.value

      all_RMSEs.real[[my_set_name]] <- real_RMSEs
      all_RMSEs.random[[my_set_name]] <- random_RMSEs
    }

    geom_mean_cors <- as.vector(unlist(geom_mean_cors))
    geom_mean_cors <- geom_mean_cors[!is.na(geom_mean_cors) & !is.infinite(geom_mean_cors)]
    geom_mean_cor.random[[my_set_name]] <- mean(geom_mean_cors)

    geom_mean_expr[[my_set_name]] <- get_geom_mean_expr(t.tpm.genes_in_set)

  } # close IF (not empty gene set)
  } # close FOR (iteration across gene sets)


  p_value_df <- data.frame("geneset"=names(p_value_list),
                           "p"=unlist(p_value_list),
                           "neglog10"=get_neg_log10(p_value_list))
  p_value_df$p.adj <- p.adjust(p_value_df$p, method = "BH")
  p_value_df$neglog10.adj <- get_neg_log10(p_value_df$p.adj)

  low_power.p_value_df <- data.frame("geneset"=names(low_power.p_value_list),
                                     "p"=unlist(low_power.p_value_list),
                                     "neglog10"=get_neg_log10(low_power.p_value_list))
  low_power.p_value_df$p.adj <- p.adjust(low_power.p_value_df$p, method = "BH")
  low_power.p_value_df$neglog10.adj <- get_neg_log10(low_power.p_value_df$p.adj)

  # transform "intercepts" from named list to dataframe
  # DELETE_intercepts
  # intercepts <- stack(intercepts) %>% rename("intercept"="values", "geneset"="ind")
  geom_mean_expr.df <- stack(geom_mean_expr) %>% rename("geom_mean_expr"="values", "geneset"="ind")
  # DELETE_mean_cor
  # mean_cors <- stack(mean_cors) %>% rename("mean_cor"="values", "geneset"="ind")

  geom_mean_cor.real.df <- stack(geom_mean_cor.real) %>% rename("GOI_cor.geneset"="values", "geneset"="ind")
  geom_mean_cor.random.df <- stack(geom_mean_cor.random) %>% rename("GOI_cor.background"="values", "geneset"="ind")


  geom_mean_cor.df <- merge(geom_mean_cor.real.df,  geom_mean_cor.random.df, by="geneset")
  geom_mean_cor.df$logFC <- log2(geom_mean_cor.df$GOI_cor.geneset/mean(geom_mean_cor.df$GOI_cor.background))


  message("ENDED calculating stats ", format(Sys.time(), "%H:%M:%S"))


  return(list("p_value_df"=p_value_df,
              "low_power.p_value_df"= low_power.p_value_df,
              "all_RMSEs.real"=all_RMSEs.real,
              "all_RMSEs.random"=all_RMSEs.random,
              # DELETE_intercepts
              # "intercepts"=intercepts,
              # DELETE_mean_cor
              # "mean_cor"=mean_cors,
              "geom_mean_cor.df"=geom_mean_cor.df,
              "geom_mean_expr.df"=geom_mean_expr.df
  ))

}


#' Clean a given expression data frame, i.e. remove all genes with too many missing values and all genes with 0 variance
#'
#' @param my_df data frame holding expression data with column=sample, row=gene
#'
#' @return cleaned data frame
clean_df <- function(my_df) {
  print(paste(nrow(my_df[rowSums(is.na(my_df)) > 0.25,]),
              " rows with more than 25% NaNs. Removing ...", sep=""))
  my_df <- my_df[rowSums(is.na(my_df)) <= 0.25,]

  print(paste(nrow(my_df[apply(my_df,1, var) == 0,]),
              " rows with variance==0. Removing ...", sep=""))
  my_df <- my_df[apply(my_df,1, var) != 0,]

  return(my_df)
}



# Note: this function either takes
#
#' Create Gene-COCOA input from expression data.
#'
#' @param expr This can either be a string holding a path to a table with expression data (tsv / csv) or a data frame already loaded into R.Columns should correspond to sample names, rows to gene names.
#' @param GOI String specifying the gene symbol of the GOI
#' @param clean Bool specifying whether the expression data should be cleaned (see \code{\link{clean_df}}) or not.
#'
#' @return Returns a list of two data frames with $GOIexpr = expression of gene of interest, $expr_df = expression info on all other genes.
#' @export
#'
#' @examples
#' \dontrun{
#' # let CAD_disease: data frame holding expression data in coronary artery disease with rows=genes, columns=samples
#'
#' expr_info <- get_expr_info(expr=CAD_disease, GOI="PDGFD")
#' }
get_expr_info <- function(expr, GOI, clean=TRUE) {
  if (is.data.frame(expr)) {
    get_expr_info.df(expr, GOI, clean)
  }
  else if (is.character(expr)) {
    get_expr_info.file(expr, GOI, clean)
  } else {
    stop("pass either a data frame containing the expression data or a string representing a path to the expression data")
  }
}



#' Takes a data frame holding expression data with rows = genes, columns = samples and transforms it to Gene-COCOA input format.
#'
#' @param expr_df
#' @param GOI
#' @param clean
#'
#' @keywords Internal
#'
get_expr_info.df <- function(expr_df, GOI, clean=TRUE) {
  rownames(expr_df) <- gsub("-", "_", rownames(expr_df))

  # get GOI expression as a column
  GOI_expr <- data.frame(t(expr_df[GOI,]))
  # remove GOI expression from expression matrix as to avoid autocorrelation
  expr_df <- expr_df[!rownames(expr_df)== GOI,]
  # clean expression matrix, if needed
  if (clean) {
    expr_df <- clean_df(expr_df)
  }

  return(list("expr_df"=expr_df,
              "GOI_expr"=GOI_expr))
}


#' Takes a path to a csv / tsv of expression values (rows=genes, columns=samples) and transforms it to Gene-COCOA input format.
#'
#' @param expr_file
#' @param GOI
#' @param clean
#'
#' @keywords Internal
#'
get_expr_info.file <- function(expr_file, GOI, clean=TRUE) {
  dt <- data.table::fread(expr_file)
  mat <- as.matrix(dt, rownames=1)
  expr_df <- as.data.frame(mat)
  rownames(expr_df) <- gsub("-", "_", rownames(expr_df))

  # get GOI expression as a column
  GOI_expr <- data.frame(t(expr_df[GOI,]))
  # remove GOI expression from expression matrix as to avoid autocorrelation
  expr_df <- expr_df[!rownames(expr_df)== GOI,]
  # clean expression matrix, if needed
  if (clean) {
    expr_df <- clean_df(expr_df)
  }

  return(list("expr_df"=expr_df,
              "GOI_expr"=GOI_expr))
}




#' Tests if a gene set is large enough to provide for enough distinguishable gene-combinations in a given number of (boostrap runs)
#'
#' @param samplesize The number of genes in one sample to be drawn.
#' @param my_set_size The size of the pool to draw from (size of this MSigDB gene set)
#' @param nsims Number of bootstrap runs
#'
#' @keywords Internal
#'
enough_power <- function(samplesize, my_set_size, nsims) {

  # if there there are more unique combinations possible than there are boot strap runs,
  # gene set has enough statistical power
  if (choose(n=my_set_size, k=samplesize)>nsims) {
    return(TRUE)
  } else { # low-power gene set
    return(choose(n=my_set_size, k=samplesize))
    # return number of bootstraps with distinct gene combinations
  }


}


#' Returns -log10 of given p-value list
#'
#' @param p_value_list
#'
#' @return
#' @keywords Internal
#'
get_neg_log10 <- function(p_value_list) {

  p_value_vec <- as.vector(unlist(p_value_list))

  # if the gene set is too small to rescale
  if (length(p_value_list)<2) {
    return(p_value_list)
  } else {
    # if there is a p=0, rescale the p value list to [min!=0,max]
    # else, keep as is
    pvalues_for_neglog <- if (min(!is.na(p_value_vec))==0)
      scales::rescale(p_value_vec,
                      to=c(sort(p_value_vec, FALSE)[2], max(!is.na(p_value_vec))))
    else p_value_vec

  }

  return((-1)*log10(pvalues_for_neglog))

}


#' Returns root mean square error of a given regression model.
#'
#' @param my_model
#'
#' @keywords Internal
#'
rmse <- function(my_model) {
  return(sqrt(mean(residuals(my_model)^2)))
}



#' Computes the intercept of each regression model
#'
#' @param GOI
#' @param GOI_expr
#' @param my_t.tpms
#'
#' @keywords Internal
get_intercept <- function(GOI, GOI_expr, my_t.tpms) {
  my_sample <- cbind.data.frame(GOI=GOI_expr,
                                my_t.tpms)
  my_model <- stats::lm(paste(GOI, "~ ."),data=my_sample)
  return(unlist(unname(coefficients(my_model)[1])))

}


#' Computes the geometric mean of the expression of a given gene set.
#'
#' @param GOI
#' @param GOI_expr
#' @param t.tpm.genes_in_set
#'
#' @keywords Internal
#'
get_geom_mean_expr <- function(t.tpm.genes_in_set) {
  # for a df, with genes as columns the psych::geometric.mean function calculates
  # the geometric mean per column, i.e. the mean gene expression across the samples
  geom_means_per_gene_list <- psych::geometric.mean(t.tpm.genes_in_set)

  # if less than 75% of the genes have a 0 expression, replace 0 with 0.01
  if ((length(geom_means_per_gene_list[geom_means_per_gene_list==0])/length(geom_means_per_gene_list))<0.75) {
    geom_means_per_gene_list[geom_means_per_gene_list==0] <- 0.01
  }
  # otherwise, the geometric mean for this dataset will return 0
  # if at least one value in the data set is 0, geometric mean returns 0


  # to obtain the base expression of the gene set,
  # take again the geom. mean, this time across all gene means
  base_mean_expr <- psych::geometric.mean(unlist(geom_means_per_gene_list))

  return(base_mean_expr)

}


#' Computes the arithmetic mean of the coexpression between each gene in a set and the GOI.
#'
#' @param GOI
#' @param GOI_expr
#' @param my_t.tpms
#'
#' @keywords Internal
#'
get_mean_cor <- function(GOI, GOI_expr, my_t.tpms) {
  my_sample <- cbind.data.frame(GOI=GOI_expr,
                                my_t.tpms)
  # get correlation matrix of given genes (including GOI)
  my_cor <- cor(my_sample)
  # Fisher-transform to z-values
  my_z <- psych::fisherz(my_cor)
  # get only lower triangle and cast to vector
  my_vector <- as.vector(my_z[lower.tri(my_z)])
  # keep only numericals in vector
  my_vector <- my_vector[!is.na(my_vector) & !is.infinite(my_vector)]
  # compute mean and transform back z->r
  my_mean_cor <- psych::fisherz2r(mean(my_vector))
  return(my_mean_cor)
}

#' Computes the geometric mean coexpression between a given gene set and the GOI.
#'
#' @param GOI
#' @param GOI_expr
#' @param my_t.tpms
#'
#' @keyowrds Internal
#'
get_geom_mean_cor <- function(GOI, GOI_expr, my_t.tpms) {
  my_sample <- cbind.data.frame(GOI=GOI_expr,
                                my_t.tpms)
  # get correlation matrix of given genes (including GOI)
  my_cor <- cor(my_sample)
  my_gm <- psych::geometric.mean(my_cor[lower.tri(my_cor)])

  return(my_gm)

}


#' Gets n random gene samples from a given gene set and computes regression models.
#'
#' @param GOI
#' @param GOI_expr
#' @param my_t.tpms
#' @param samplesize
#' @param nsims
#' @param geom_mean_cors
#'
#' @return Returns a list with $RMSEs=root mean square errors of all regression models, $geom_mean_cors=geometric mean correlation of each sample with GOI
#' @keywords Internal
#'
compute_x_RMSEs <- function(GOI,
                            GOI_expr,
                            my_t.tpms,
                            samplesize=10,
                            nsims=1000,
                            geom_mean_cors=FALSE) {

  my_RMSEs <- list()
  my_geom_mean_cors <- list()

  for (i in 1:nsims) {
    my_sample <- my_t.tpms[,sample(ncol(my_t.tpms), size=samplesize)]
    my_gene_names <- paste(names(my_sample), collapse=";")

    if (geom_mean_cors) { my_geom_mean_cors[[my_gene_names]] <- get_geom_mean_cor(GOI, GOI_expr, my_t.tpms = my_sample) }

    my_sample <- cbind.data.frame(GOI=GOI_expr,
                                  my_sample)
    my_model <- stats::lm(paste(GOI, "~ ."),data=my_sample)
    my_RMSEs[[my_gene_names]] <- rmse(my_model)

    #my_intercepts[[my_gene_names]] <- unlist(unname(coefficients(my_model)[1]))
  }

  return(list("RMSEs"=my_RMSEs, "geom_mean_cors"=my_geom_mean_cors))
}

#' Gets all combinations of <samplesize> genes from a given gene set and computes regression models.
#'
#' @param GOI
#' @param GOI_expr
#' @param my_t.tpms
#' @param samplesize
#' @param geom_mean_cors
#'
#' @return Returns a list with $RMSEs=root mean square errors of all regression models, $geom_mean_cors=geometric mean correlation of each sample with GOI
#' @keywords Internal
#'
compute_all_RMSEs <- function(GOI,
                              GOI_expr,
                              my_t.tpms,
                              samplesize=2,
                              geom_mean_cors=FALSE) {

  my_RMSEs <- list()
  my_geom_mean_cors <- list()

  # get all possible combinations of column indices (genes)
  allcomb <- combn(x=ncol(my_t.tpms), m=samplesize)
  # iterate over all
  for (i in 1:ncol(allcomb)) {
    my_sample <- my_t.tpms[,allcomb[,i]]
    my_gene_names <- paste(names(my_sample), collapse=";")

    if (geom_mean_cors) { my_geom_mean_cors[[my_gene_names]] <- get_geom_mean_cor(GOI, GOI_expr, my_t.tpms = my_sample) }

    my_sample <- cbind.data.frame(GOI=GOI_expr,
                                  my_sample)
    my_model <- stats::lm(paste(GOI, "~ ."),data=my_sample)
    my_RMSEs[[my_gene_names]] <- rmse(my_model)

  }
  return(list("RMSEs"=my_RMSEs, "geom_mean_cors"=my_geom_mean_cors))
}


#' Recommended visualisation of Gene-COCOA results.
#' Produces a ggplot with x=(GOI coexpression with gene set)/(GOI coexpression with random set), y=-log10(p.adj).
#' Dot size and alpha reflect the (geometric) mean expression of all genes in a respective gene set.
#'
#' @param mystats Output of \code{\link{get_stats}}.
#' @param sig_label_cutoff At which -log(p.adj) value should points be labelled? Default=2(=-log(p.adj=0.01).
#' @param sig_col_cutoff At which p.adj value should points be coloured (instead of grey)? Defaul p.adj=0.01.
#' @param sig_colour Which colour should significant dots have?
#' @param remove_outliers Should outliers (in terms of x-axis and y-axis) be removed. Default is FALSE; you should stick with this unless you get very cluttered plots.
#' @param filepath  Path + name of png to save If not provided, only the ggplot will be returned by this function.
#'
#' @return Returns a ggplot which can be stored in a variable and modified. If an absolute file name (path + file name) is provided, the ggplot with additionally be saved to a png.
#' @export
#'
#' @examples
#' \dontrun{
#' expr_info <- get_expr_info(expr=CAD_disease, GOI="PDGFD")
#' res <- get_stats(geneset_collection=get_msigigdb_genesets("HALLMARK"), GOI="PDGFD", GOI_expr=expr_info$GOI_expr, expr_df=expr_info$expr_df)
#' # show GeneCOCOA ggplot in "Plots" window
#' plot_volcano(res)
#'
#' # assign GeneCOCOA ggplot to variable
#' p <- plot_volcano(res)
#'
#' # assign GeneCOCOA ggplot to variable AND save GeneCOCOA ggplot to file
#' p <- plot_volcano(res, filepath="/path/to/my/GeneCOCOA.volcano.png")
#'}
plot_volcano <- function (mystats, sig_label_cutoff = 2, sig_col_cutoff =0.01, sig_colour = "dodgerblue4",
                          remove_outliers = FALSE, filepath = "")
{
  new_df <- merge(mystats$p_value_df, mystats$geom_mean_cor.df,
                  by = "geneset")
  new_df <- merge(new_df, mystats$geom_mean_expr.df, by = "geneset")
  new_df$label <- NA
  new_df$label[new_df$neglog10.adj >= sig_label_cutoff] <- sapply(new_df$geneset[new_df$neglog10.adj >= sig_label_cutoff],
                                                                  function(x) stringr::str_wrap(stringr::str_replace_all(x,"_", " "), width = 10))
  new_df$col <- "grey"
  new_df$col[new_df$p.adj < sig_col_cutoff] <- sig_colour
  new_df$alpha <- 0.5
  new_df$size <- scales::rescale(new_df$geom_mean_expr, to = c(5,20))

  if (remove_outliers) {
    new_df = new_df[!new_df$logFC %in% boxplot_stats(new_df$logFC)$out,]
    new_df = new_df[!new_df$neglog10.adj %in% boxplot_stats(new_df$neglog10.adj)$out,]
  }

  nudge_y_centre = diff(range(new_df$neglog10.adj))/2
  nudge_y = diff(range(new_df$neglog10.adj))/15
  nudge_x = diff(range(new_df$logFC))/15

  my_plot <- ggplot(new_df, aes(x = logFC, y = neglog10.adj,
                                label = label)) +
  geom_point(alpha = new_df$alpha, col = new_df$col,size = new_df$size) +
  geom_text_repel(aes(x = logFC,y = neglog10.adj, label = label),
                  size = 3, col = "black",
                  segment.color = sig_colour, segment.alpha = 0.4,
                  position = position_nudge_center(direction = "split"),
                  hjust = "outward", vjust = "outward") +
    # theme_minimal(base_size = 14) +
    theme_bw(base_size = 16) +
      theme(axis.text = element_text(colour = 'black'),
            axis.line = element_line(),
            panel.border = element_blank(),
            axis.ticks = element_line(colour = 'black'),
            strip.background = element_rect(colour = NA)) +
    scale_y_continuous(bquote(-log[10](P[adj]))) + # DELETE scale_y_continuous("-log10(p.adj)") +
    scale_x_continuous("impact",
                       limits = c(1.25 * (-max(abs(new_df$logFC))), 1.25 * max(abs(new_df$logFC))))

  if (nchar(filepath) > 0) {
    ggsave(filename = filepath, plot = my_plot, dpi = 300, height=10, width=10, units="in")
  }
  return(my_plot)
}


#' Comparative visualisation of the GeneCOCOA results for a GOI in two different conditions.
#' Produces a ggplot divering bar plot with x=-log10(p.adj) and y= top significant terms in chosen condition
#'
#' @param control_res Output of \code{\link{get_stats}} ("control" condition, left side of the plot).
#' @param treatment_res Output of \code{\link{get_stats}} ("treatment" condition, right side of the plot).
#' @param include_low_power Should gene sets which (for combinatorial reasons) could only be bootrstrapped <1000 times be included?
#' @param control_label How should the "control" group be labelled?
#' @param treatment_label How should the "treatment" group be labelled?
#' @param topN How many terms should be displayed? The terms are ordered by significance and only the top n terms plotted comparatively.
#' @param control_col Which colour should the "control" group bars be plotted in?
#' @param treatment_col Which colour should the "treatment" group bars be plotted in?
#' @param sort_by Should the bars be sorted by significance in the "treatment" or "control" group? (Default sort_by="treatment")
#' @param filepath  Path + name of png to save If not provided, only the ggplot will be returned by this function.
#'
#' @return Returns a ggplot which can be stored in a variable and modified. If an absolute file name (path + file name) is provided, the ggplot with additionally be saved to a png.
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # load complete Familial Hypercholesteroleamia data set
#' FH <- as.data.frame(get_dataset_expression("GSE6054"))
#' FH_disease <- FH %>% select(contains("FH"))
#' disease_input <- get_expr_info(expr=FH_disease, GOI="LDLR")
#' disease_res <- get_stats(geneset_list=hallmark_sets,
#'                          GOI="LDLR",
#'                          GOI_expr=disease_input$GOI_expr,
#'                          expr_df=disease_input$expr_df,
#'                          samplesize=2, nsims=1000)
#'
#'
#' FH_control <- FH %>% select(contains("Control"))
#' control_input <- get_expr_info(expr=FH_control, GOI="LDLR")
#' control_res <- get_stats(geneset_list=hallmark_sets,
#'                          GOI="LDLR",
#'                          GOI_expr=control_input$GOI_expr,
#'                          expr_df=control_input$expr_df,
#'                          samplesize=2, nsims=1000)
#'
#'
#' p <- plot_control_vs_treatment(control_res, treatment_res,
#'                                control_label="control", treatment_label="FH",
#'                                topN=5, sort_by="treament",
#'                                filepath="/path/to/my/LDLR.FH.diverging_bars.png")
#' p
#'}
plot_control_vs_treatment <- function(
  control_res, treatment_res, include_low_power=TRUE,
  control_label="control", treatment_label="treatment",
  control_col="#00000000", treatment_col="brown4",
  sort_by="treatment", topN=10,
  filepath="") {

  if (include_low_power) {
    control <- rbind(control_res$p_value_df, control_res$low_power.p_value_df)
    control$geneset <- gsub("_\\(\\d+)$", "", control$geneset)

    treatment <- rbind(treatment_res$p_value_df, treatment_res$low_power.p_value_df)
    treatment$geneset <- gsub("_\\(\\d+)$", "", treatment$geneset)

  }  else {
    control <- control_res$p_value_df
    treatment <- treatment_res$p_value_df
  }


  control$condition <- rep(control_label, nrow(control))
  control$neglog10.adj <- (-1)*control$neglog10.adj
  control <- control %>% arrange(p.adj)
  control$rank <- rownames(control)

  treatment$condition <- rep(treatment_label, nrow(treatment))
  treatment <- treatment %>% arrange(p.adj)
  treatment$rank <- rownames(treatment)



  if (sort_by=="control") {
    control$geneset <- factor(control$geneset, levels=control$geneset)
    treatment <- treatment[match(control$geneset, treatment$geneset),]
    treatment$geneset <- factor(treatment$geneset, levels=control$geneset)
  } else {
    treatment$geneset <- factor(treatment$geneset, levels=treatment$geneset)
    control <- control[match(treatment$geneset, control$geneset),]
    control$geneset <- factor(control$geneset, levels=treatment$geneset)
  }

  dat <- rbind(control[c(1:topN),], treatment[c(1:topN),])
  dat$geneset <- sapply(dat$geneset,  function(x) stringr::str_wrap(stringr::str_replace_all(x, "_", " "),width = 30))

  cols=c(control_col, treatment_col)
  cols <- setNames(cols, c(control_label, treatment_label))

  my_plot <- ggplot(dat, aes(y=fct_inorder(geneset), x=neglog10.adj, fill=condition)) +
    geom_bar(stat="identity", position="identity", colour="black", width=0.75) +
    theme_bw() +
    theme(axis.text = element_text(colour = 'black', size=16),
          axis.title.x=element_text(size=16),
          legend.title=element_text(size=16),
          legend.text=element_text(size=16),
          legend.position = "top",
          axis.line = element_line(),
          panel.border = element_blank(),
          axis.ticks = element_line(colour = 'black'),
          strip.background = element_rect(colour = NA)) +
    scale_fill_manual(values=cols) +
    xlab(bquote(-log[10](P[adj]))) +
    scale_y_discrete("", limits=rev)

  if (nchar(filepath) > 0) {
    ggsave(filename = filepath, plot = my_plot, height=5, width=10, units="in")
  }
  return(my_plot)
}



#' Get combination of genes which model the GOI best (based on root mean square error of regression model).
#' Both, samples from real gene sets and from random combinations are stored and compared.
#'
#' @param mystats Output of \code{\link{get_stats}}.
#' @param GOI A string specifying the GOI gene symbol.
#' @param expr The complete expression data frame.
#' @param output If set to "df", a data frame will be returned which provides an overview of all predictor combinations and their origins (gene set name or "RANDOM"), together with their RMSEs. If set to "cor", the coexpression matrix of predictors and GOI is returned. Otherwise, simply the best set is returned together with its origin and RMSE.
#'
#' @return If output=="df": data frame of rankings of gene combinations. If output=="cor": coexpression matrix between the best gene combinations and the GOI. Else: df row holding best gene combination and its RMSE.
#' @export
#'
#' @examples
#' \dontrun{
#' #' # let CAD_disease: data frame holding expression data in coronary artery disease with rows=genes, columns=samples
#'
#' expr_info <- get_expr_info(expr=CAD_disease, GOI="PDGFD")
#' res <- get_stats(geneset_collection=get_msigigdb_genesets("HALLMARK"), GOI="PDGFD", GOI_expr=expr_info$GOI_expr, expr_df=expr_info$expr_df)
#' gene_comb_ranking <- get_best_predictors(mystats=res, GOI="PDGFD", expr=CAD_disease, output="df")
#' }
get_best_predictors <- function(mystats, GOI, expr, output = FALSE) {
  # flatten the list of lists into dataframe with
  # rownames=combination of genes, column=RMSEs
  all_RMSEs.real.df <- as.data.frame(unlist(mystats$all_RMSEs.real))
  # the combination of genes is given in following format GENESETNAME.gene1;gene2;gene3;gene4
  # cast to column named "geneset.comb"
  all_RMSEs.real.df$geneset.comb <- rownames(all_RMSEs.real.df)
  names(all_RMSEs.real.df) <- c("RMSE", "geneset.comb")
  rownames(all_RMSEs.real.df) <- NULL
  # add new column holding only the gene set name (pattern up to first dot in geneset.comb)
  all_RMSEs.real.df <- all_RMSEs.real.df %>% mutate("geneset_name"=sub("[.].*", "", geneset.comb) )
  # add column holding only the combination of genes (everything after first dot)
  all_RMSEs.real.df <- all_RMSEs.real.df %>% mutate("comb"=sub("^.*?\\.", "", geneset.comb) )
  all_RMSEs.real.df$geneset.comb <- NULL
  rownames(all_RMSEs.real.df) <- NULL


  all_RMSEs.random.df <- as.data.frame(unlist(mystats$all_RMSEs.random))
  # the combination of genes is given in following format GENESETNAME.gene1;gene2;gene3;gene4
  # cast to column named "geneset.comb"
  all_RMSEs.random.df$geneset.comb <- rownames(all_RMSEs.random.df)
  names(all_RMSEs.random.df) <- c("RMSE", "geneset.comb")
  rownames(all_RMSEs.random.df) <- NULL
  # add column holding only the combination of genes (everything after first dot)
  all_RMSEs.random.df <- all_RMSEs.random.df %>% mutate("comb"=sub("^.*?\\.", "", geneset.comb) )
  # set gene set name to random
  all_RMSEs.random.df$geneset_name <- rep("RANDOM", nrow(all_RMSEs.random.df))
  all_RMSEs.random.df$geneset.comb <- NULL

  all_RMSEs <- rbind(all_RMSEs.real.df, all_RMSEs.random.df)

  best_row <- all_RMSEs[which.min(all_RMSEs$RMSE),]
  rownames(best_row) <- NULL

  best_set <- unlist(stringr::str_split(best_row, ";"))

  expr_best_set <- expr[rownames(expr) %in% c(GOI, best_set),]
  cor_best_set <-  cor(t(expr_best_set))


  if (output=="df") {
    return(all_RMSEs)
  } else if (output=="cor") {
    return(cor_best_set)
  } else {
    return(as.matrix(best_row))
  }

}



#' Plots p-value ranking (x-axis=adjusted -log10(p), y= gene set names, ranked).
#' This function works but \code{\link{plot_volcano}} is more informative.
#' @param my_stats Output of \code{\link{get_stats}}
#' @param filepath Path + name of png to save If not provided, only the ggplot will be returned by this function.
#'
#' @return Returns a ggplot which can be stored in a variable and modified. If an absolute file name (path + file name) is provided, the ggplot with additionally be saved to a png.
#' @export
#'
#'
plot_basic_p_ranking <- function (my_stats, filepath)
{
  ggplot(my_stats$p_value_df %>% arrange(neglog10.adj, rev(geneset)) %>%
           mutate(geneset = factor(geneset, levels = geneset)),
         aes(x = geneset, y = neglog10.adj, group = 1, label = formatC(p.adj,
                                                                       digits = 6))) + scale_x_discrete(name = "gene set") +
    scale_y_continuous(name = "-log10(p.adj)") + geom_point() +
    theme_minimal() +
    geom_line() + geom_text(size = 6,
                            hjust = "inward", vjust = "center") + coord_flip() +
    theme(
      axis.text = element_text(size = 18), axis.text.y = element_text(angle = 0),
      axis.title = element_text(size = 21), legend.title = element_text(size = 18),
      legend.text = element_text(size = 18))
  ggsave(filename = filepath,  width = 10, height = 25,
         units = "in", dpi = 300)
}



#' Returns a ggplot with x=-log10(p.adj), y=mean coexpression of the gene set with the GOI.
#' This function works but \code{\link{plot_volcano}} is more informative.
#'
#' @param mystats Output of \code{\link{get_stats}}.
#' @param sig_label_cutoff
#' @param sig_colour
#' @param remove_outliers
#' @param filepath  Path + name of png to save If not provided, only the ggplot will be returned by this function.
#'
#' @return Returns a ggplot which can be stored in a variable and modified. If an absolute file name (path + file name) is provided, the ggplot with additionally be saved to a png.
#' @export
#'
#'
plot_p_vs_geom_mean_cor <- function(mystats, sig_label_cutoff=0.05, sig_colour="dodgerblue4", remove_outliers=FALSE, filepath) {

  new_df <- merge(mystats$p_value_df, mystats$geom_mean_cor.df,
                  by="geneset")
  new_df <- merge(new_df, mystats$geom_mean_expr.df,
                  by = "geneset")
  new_df$label <- NA
  new_df$label[new_df$p.adj<sig_label_cutoff] <- sapply(new_df$geneset[new_df$p.adj<sig_label_cutoff],  function(x) stringr::str_wrap(stringr::str_replace_all(x, "_", " "),width = 10))
  new_df$col <- "grey"
  new_df$col[new_df$p.adj<sig_label_cutoff] <- sig_colour
  new_df$alpha <- 0.5
  new_df$size <- scales::rescale(new_df$geom_mean_expr, to = c(5,20))


  if (remove_outliers) {
    new_df = new_df[!new_df$GOI_cor.geneset %in% boxplot_stats(new_df$GOI_cor.geneset)$out,]
  }

  nudge_y_centre=diff(range(new_df$GOI_cor.geneset))/2
  nudge_y=diff(range(new_df$GOI_cor.geneset))/5

  my_plot <- ggplot(new_df, aes(x=GOI_cor.geneset, y=neglog10.adj,
                                label=label)) +
    geom_point(alpha=new_df$alpha, col=new_df$col, size=10) +
    geom_text_repel(aes(x=neglog10.adj, y=GOI_cor.geneset, label=label),
                    col=sig_colour,
                    min.segment.length = unit(0, 'lines'),
                    segment.color=sig_colour,
                    segment.alpha=0.5,
                    nudge_y = ifelse(new_df$GOI_cor.geneset >nudge_y_centre, nudge_y, -nudge_y)) +
    theme_minimal(base_size=14) +
    scale_y_continuous("Mean coexpression with GOI") +
    scale_x_continuous("Significane of coexpression with GOI [-log10(p.adj)]")

  # CHECK HOW TO ADD:
    # geom_text_repel(aes(x = logFC,y = neglog10.adj, label = label),
    #                 size = 3, col = "black",
    #                 segment.color = sig_colour, segment.alpha = 0.4,
    #                 position = position_nudge_center(direction = "split"),
    #                 hjust = "outward", vjust = "outward") +
    # theme_bw(base_size = 16) +
    # theme(axis.text = element_text(colour = 'black'),
    #       axis.line = element_line(),
    #       panel.border = element_blank(),
    #       axis.ticks = element_line(colour = 'black'),
    #       strip.background = element_rect(colour = NA)) +
    # scale_y_continuous(bquote(-log[10](P[adj]))) + # DELETE scale_y_continuous("-log10(p.adj)") +
    # scale_x_continuous("impact",
    #                    limits = c(1.25 * (-max(abs(new_df$logFC))), 1.25 * max(abs(new_df$logFC))))



  if (nchar(filepath)>0) {
    ggsave(filename=filepath, plot=my_plot, dpi=300)
  }

  return(my_plot)

}



#' Return a collection of curated gene sets from the MSigDB.
#'
#' @param geneset Name of gene set. Can be "HALLMARK", "POSITIONAL", "KEGG", "PID", "REACTOME", "WP" (for WikiPathways), "GOBP" (for GO: Biological Process), "GOCC" (for GO: Cellular Component), "GOMF" (for GO: Molecular Function), "ONCO" (for oncogenic signature). If you are unsure which collection of gene sets to choose for comparison, "HALLMARK" is recommended.
#' @param remove_prefix Bool. Should the collection prefix (e.g. "HALLMARK_") be removed from the gene set names. (FALSE might result in cluttered plots.)
#'
#' @return A list of named lists. Each sublist is named by a gene set, and its contents are the symbols of all genes in this set as strings.
#' @export
#'
#' @examples
#' \dontrun{
#' # assign gene set collection to variable and pass to get_stats
#' all_hallmark_gs <- get_msigdb_genesets("HALLMARK")
#' res <- get_stats(geneset_collection=all_hallmark_gs, GOI="PDGFD",  ...)
#'
#' # call directly in get_stats()
#' res <- get_stats(geneset_collection=get_msigdb_genesets("HALLMARK"), ...)
#' }
get_msigdb_genesets <- function(genesets="", remove_prefix=TRUE) {
  if (genesets=="HALLMARK") {
    gs.table<- msigdbr(species="Homo sapiens", category="H")
  } else if (genesets=="POSITIONAL") {
    gs.table<- msigdbr(species="Homo sapiens", category="C1")
  } else if (genesets=="BIOCARTA") {
    gs.table <- msigdbr(species="Homo sapiens", category="C2", subcategory = "BIOCARTA")
  } else if(genesets=="KEGG") {
    gs.table <- msigdbr(species="Homo sapiens", category="C2", subcategory = "KEGG")
  } else if(genesets=="PID") {
    gs.table <- msigdbr(species="Homo sapiens", category="C2", subcategory = "PID")
  } else if(genesets=="REACTOME") {
    gs.table <- msigdbr(species="Homo sapiens", category="C2", subcategory = "REACTOME")
  } else if(genesets=="WP") { # WIKI PATHWAYS
    gs.table <- msigdbr(species="Homo sapiens", category="C2", subcategory = "WIKIPATHWAYS")
  } else if(genesets=="GOBP") { # GO BIOLOGICAL PROCESS
    gs.table <- msigdbr(species="Homo sapiens", category="C5", subcategory = "BP")
  } else if(genesets=="GOCC") { # GO CELLULAR COMPONENT
    gs.table <- msigdbr(species="Homo sapiens", category="C5", subcategory = "CC")
  } else if(genesets=="GOMF") { # GO MOLECULAR FUNCTION
    gs.table <- msigdbr(species="Homo sapiens", category="C5", subcategory = "MF")
  } else if(genesets=="ONCO") { # ONCOGENIC SIGNATURE
    gs.table <- msigdbr(species="Homo sapiens", category="C8")
  } else {
    message("Unrecognised DB. Currently supported geneset arguments:")
    message("HALLMARK\t50 MSigDB hallmark gene sets.\n\t\tIf you are unsure which gene sets to compare, choose the hallmark gene sets.")
    message("POSITIONAL\tgene sets corresponding to human chromosome cytogenetic bands")
    message("KEGG\t\tgene sets derived from the KEGG pathway database")
    message("PID\t\tgene sets derived from the PID pathway database")
    message("REACTOME\tgene sets derived from the REACTOME pathway database")
    message("WP\t\tgene sets derived from the WikiPathways database")
    message("GOBP\t\tgene sets derived from the GO Biological Process ontology")
    message("GOCC\t\tgene sets derived from the GO Cellular Component ontology")
    message("GOMF\t\tgene sets derived from the GO Molecular Funciton ontology")
    message("ONCO\t\tgene sets that represent signatures of cellular pathways which are often disregulated in cancer.")
    return(NA)
  }

  gs.symbol_lists <- gs.table %>%
    split(f = as.factor(.$gs_name)) %>% # split into list of named tibbles (per individual gene set)
    lapply(FUN = function(x) subset(x, select=c("gene_symbol"))) # retain only the gene symbol column of each gene set

  if (remove_prefix) {
    names(gs.symbol_lists) <- sub(paste0('^', genesets, '_'), '', names(gs.symbol_lists))

  }

  return(gs.symbol_lists)

}





#' Returns a ggplot with x=-log10(p.adj), y=mean coexpression of the gene set with the GOI.
#' This function works but \code{\link{plot_volcano}} is more informative.
#'
#' @param treatment_df Expression matrix of treatment condition (gene sets associated more strongly with GOI in this condition will appear on the right-hand side of the plot).
#' @param control_df Expression matrix of treatment condition (gene sets associated more strongly with GOI in this condition will appear on the right-hand side of the plot).
#' @param GOI A string specifying the GOI gene symbol.
#' @param geneset_collection A list of lists. Each list holds gene symbols of a specific gene set. Can be generated with \code{\link{get_msigdb_genesets}}
#' @param output_dir  If the Laplace distribution fitting results should not only be returned as stats but also plotted, provide the output path for the plots. Default: working directory.
#' @param plot_ECDF  Should the ECDF of the Laplace fitting be plotted? Default: true.
#' @param plot_histo  Should the histogram of the Laplace fitting be plotted? Default: true.
#' @param histo_breaks How many breaks should the histogram have? Default: 100.
#' @return Returns a named list with "location" (the location estimate / mean of the fitted Laplace distribution), "scale" (the scale / standard deviation of the fitted Laplace distribution) and "ks.test" (the result of the Kolmogorov Smirnov test).
#' @export
estimate_laplace_parameters <- function(treatment_df=data.frame(), control_df=data.frame(),
                                        GOI, geneset_collection,
                                        output_dir=getwd(),
                                        plot_ECDF=TRUE,
                                        plot_histo=TRUE,
                                        histo_breaks=100) {
  # randomise the data frames
  disease_rownames <- rownames(treatment_df)
  treatment_df <- as.data.frame(lapply(treatment_df, sample))
  rownames(treatment_df) <- disease_rownames
  control_rownames <- rownames(control_df)
  control_df <- as.data.frame(lapply(control_df, sample))
  rownames(control_df) <- control_rownames


  expr_info.disease <- get_expr_info(expr=treatment_df, GOI=GOI)
  res_disease <- get_stats(geneset_collection=geneset_collection, GOI=GOI, GOI_expr=expr_info.disease$GOI_expr, expr_df=expr_info.disease$expr_df, samplesize=2, nsims=1000)

  expr_info.control <- get_expr_info(expr=control_df, GOI=GOI)
  res_control <- get_stats(geneset_collection=geneset_collection, GOI=GOI, GOI_expr=expr_info.control$GOI_expr, expr_df=expr_info.control$expr_df, samplesize=2, nsims=1000)


  disease <- rbind(as.data.frame(res_disease$p_value_df), as.data.frame(res_disease$low_power.p_value_df))
  control <- rbind(as.data.frame(res_control$p_value_df), as.data.frame(res_control$low_power.p_value_df))

  # remove low power suffix (number of samplings in brackets)
  disease$geneset <- gsub("_\\(\\d+)$", "", disease$geneset)
  control$geneset <- gsub("_\\(\\d+)$", "", control$geneset)


  diff_df <- merge(control[,c("geneset", "p", "p.adj")], disease[,c("geneset", "p", "p.adj")], by="geneset") %>% na.omit()
  names(diff_df) <- c("geneset", "p_control", "p_adj_control", "p_disease", "p_adj_disease")
  diff_df$p_ratio <- (diff_df$p_disease / diff_df$p_control) # if small, then "more significant" in disease. If large, "more significant" in control.
  diff_df$log10_p_ratio <- log10(diff_df$p_ratio) # if negative, "more significant" in disease. If positive, "more significant" in control.
  diff_df$neglog10_p_ratio <- (-1)*diff_df$log10_p_ratio # if positive, "more significant" in disease. If negative, "more significant" in control.



  laplace_fit <- fitdist(diff_df$neglog10_p_ratio, "laplace",
                         start = list(location = mean(diff_df$neglog10_p_ratio),
                                      scale = sd(diff_df$neglog10_p_ratio)))
  location_est <- laplace_fit$estimate["location"]
  scale_est <- laplace_fit$estimate["scale"]


  ks_test <- ks.test(rlaplace(n=500, location=laplace_fit$estimate["location"], scale=laplace_fit$estimate["scale"]),
                     rlaplace(n=500, location=0, scale=scale_est))
  ks_test_p <- as.character(as.numeric(ks_test$p.value))




  x_range <- range(c(rlaplace(500, location=0, scale=1),
                     rlaplace(500, location=laplace_fit$estimate["location"],
                              scale=laplace_fit$estimate["scale"])))

  y_range <- range(c(ecdf(rlaplace(n=500, location=0, scale=1))(x_range),
                     ecdf(rlaplace(n=500, location=laplace_fit$estimate["location"],
                                   scale=laplace_fit$estimate["scale"]))(x_range)))
  # Position legend based on data ranges
  legend_x <- min(x_range) + (0.05*diff(x_range)) # Slight offset from the left
  legend_y <- max(y_range) - (0.05*diff(y_range)) # Slight offset from the top

  if (plot_ECDF) {
    # Empirical cumulative distribution function of general and fitted Laplace distribution
    png(paste0(output_dir, "/Laplace_vs_custom.ECDF.png"), width=6, height=5, units="in", res=300)
    plot(ecdf(rlaplace(n=500, location=0, scale=1)), do.points=FALSE, lwd=2, xlim=c(-11, 11), col="grey60",
         main="")
    plot(ecdf(rlaplace(n=500, location=laplace_fit$estimate["location"], scale=laplace_fit$estimate["scale"])),
         add=TRUE, do.points=FALSE, lwd=2, col="darkred")
    # legend(-10.5, 0.9, legend=c("Laplace(m=0, sd=1)", "Laplace(data)"),
    #        col=c("grey60", "darkred"), lty=1, cex=0.8)
    legend(legend_x, legend_y, legend=c("Laplace(m=0, sd=1)", "Laplace(data)"),
           col=c("grey60", "darkred"), lty=1, cex=0.8, xjust=0, yjust=1)
    dev.off()
  }


  if (plot_histo) {
    png(paste0(output_dir, "/Laplace_vs_custom.histo.png", width=5, height=3, units="in", res=300))
    # calculate limits for legend and p-value annotation placements
    xlim <- range(diff_df$neglog10_p_ratio, na.rm = TRUE)
    ylim <- c(0, max(hist(diff_df$neglog10_p_ratio, breaks = histo_breaks, plot = FALSE)$density, na.rm = TRUE))
    # plot
    hist(diff_df$neglog10_p_ratio, prob = TRUE, main = "",
         xlab = bquote(-log[10](P[disease] / P[control])), col = "lightgrey", breaks = histo_breaks)
    curve(dlaplace(x, location = location_est, scale = scale_est), col = "darkorange2", lwd = 2, add = TRUE)
    curve(dlaplace(x, location = 0, scale = 1), col = "darkred", lwd = 2, add = TRUE)
    legend(
      x = xlim[1] + 0.05 * diff(xlim), # Slight offset from the left edge
      y = ylim[2] - 0.05 * diff(ylim), # Slight offset from the top edge
      legend = c("Laplace(data)", "Laplace(m=0, sd=1)"),
      col = c("darkorange2", "darkred"), lty = 1, cex = 0.8
    )
    text(
      x = xlim[2] - 0.05 * diff(xlim), # Slight offset from the right edge
      y = ylim[2] - 0.05 * diff(ylim), # Slight offset from the top edge
      labels = paste("P-value:", ks_test_p), # P-value text
      col = "black", cex = 0.8, adj = c(1, 1) # adj = c(1, 1) for upper right alignment
    )
    dev.off()
  }


  return ( list(
    "location"=location_est,
    "scale"=scale_est,
    "ks.test"=ks_test
  ))

} # END FUNCTION estimate_laplace_parameters



#' Differential analysis of two given GeneCOCOA results for conditions treatment and control. Returns a dataframe with "neglog10_p_ratio" indicating direction & strength of change in association between a gene set and the GOI (>0: more strongly associated in treatment, <0: more strongly associated in control), and "differential_p" the corresponding significance of the change in association.
#'
#' @param control_res Output of \code{\link{get_stats}} ("control" condition, left side of the plot).
#' @param treatment_res Output of \code{\link{get_stats}} ("treatment" condition, right side of the plot).
#' @param GOI A string specifying the GOI gene symbol.
#' @param geneset_collection A list of lists. Each list holds gene symbols of a specific gene set. Can be generated with \code{\link{get_msigdb_genesets}}
#' @param include_low_power Should gene sets which (for combinatorial reasons) could only be bootrstrapped <1000 times be included?
#' @param laplace_parameters Which Laplace distribution should be fitted to the data in order to derive p-values for the differential association values? Options: "Default": use precomputed Laplace-values (location=0.0076, scale=1.6614), "Laplace": location=0, scale=1, "Estimate": estimate parameters from given empirical data - cave: this will take longer than the other parameters!
#' @param treatment_df Only needed for laplace_parameters="Estimate". Expression matrix of treatment condition (gene sets associated more strongly with GOI in this condition will appear on the right-hand side of the plot).
#' @param control_df Only needed for laplace_parameters="Estimate". Expression matrix of treatment condition (gene sets associated more strongly with GOI in this condition will appear on the right-hand side of the plot).
#' @return Returns differential GeneCOCOA results: a data frame with with columns "geneset", "p_control", "p_adj_control", "p_disease". "p_adj_disease", "p_ratio", "log10_p_ratio", "neglog10_p_ratio", "differential_p"
#' @export
#'
#'
#'
#'
get_differential_GeneCOCOA_results <- function(treatment_res, control_res,
                                               GOI,
                                               geneset_collection,
                                               include_low_power=TRUE,
                                               laplace_parameters="Default", # options: "Default", "Laplace", "Estimate"
                                               treatment_df=data.frame(), # only needed for laplace_parameters="Estimate"
                                               control_df=data.frame() # only needed for laplace_parameters="Estimate"
)
{
  if (include_low_power) {
    disease <- rbind(as.data.frame(treatment_res$p_value_df), as.data.frame(treatment_res$low_power.p_value_df))
    control <- rbind(as.data.frame(control_res$p_value_df), as.data.frame(control_res$low_power.p_value_df))

    # remove low power suffix (number of samplings in brackets)
    disease$geneset <- gsub("_\\(\\d+)$", "", disease$geneset)
    control$geneset <- gsub("_\\(\\d+)$", "", control$geneset)

  } else {
    disease <- as.data.frame(treatment_res$p_value_df)
    control <- as.data.frame(control_res$p_value_df)
  }



  control <- control[match(disease$geneset, control$geneset),]
  disease$geneset <- factor(disease$geneset, levels=disease$geneset)
  control$geneset <- factor(control$geneset, levels=disease$geneset)

  diff_df <- merge(control[,c("geneset", "p", "p.adj")], disease[,c("geneset", "p", "p.adj")], by="geneset") %>% na.omit()
  names(diff_df) <- c("geneset", "p_control", "p_adj_control", "p_disease", "p_adj_disease")
  diff_df$p_ratio <- (diff_df$p_disease / diff_df$p_control) # if small (0-1), then "more significant" in disease. If large (several 1000), "more significant" in control.
  diff_df$log10_p_ratio <- log10(diff_df$p_ratio) # if negative, "more significant" in disease. If positive, "more significant" in control.
  diff_df$neglog10_p_ratio <- (-1)*diff_df$log10_p_ratio # if positive, "more significant" in disease. If negative, "more significant" in control.




  if (laplace_paramters=="Default") {
    location_est <- 0.007595472
    scale_est <-   1.661423
  } else if (laplace_parameters=="Laplace") {
    location_est <- 0
    scale_est <-   1
  } else if (laplace_parameters=="Estimate") {
    # run function


    # randomise data frames
    # run gene COCOA for


  } else {
    message("Please choose between values 'Default', 'Laplace', or 'Estimate'.")
    message("Data frame will be returned without p-values.")
  }


  diff_df$differential_p <- 2*(plaplace(abs(diff_df$neglog10_p_ratio), location = location_est, scale = scale_est, lower.tail = FALSE))

  return(diff_df)
}














#' Comparative visualisation of the GeneCOCOA results for a GOI in two different conditions.
#' Produces a ggplot divering bar plot with x=-log10(p.adj) and y= top significant terms in chosen condition
#'
#' @param diff_df Output of \code{\link{get_differential_GeneCOCOA_results}}
#' @param sig_label_cutoff At which differential p-value should points be labelled? Default: 0.05
#' @param output_dir  Output path for the plots. Default: working directory.
#' @param GOI A string specifying the GOI gene symbol.
#' @param condition A string specifying the treatment condition examined in contrast to control.
#'
#' @return Returns a list with two ggplots, a MA plot wand a volcano plot. Additionally saves them to files.
#' @export
#'

plot_differential_GeneCOCOA_results <- function(diff_df, sig_label_cutoff=0.05, output_dir=getwd(),
                                                GOI="GOI", condition="" # only needed for laplace_parameters="Estimate"
)
{
  

  plot_prefix <- png(paste0(output_dir, "/", condition, ".", GOI))
  
  plot_df <- diff_df
  plot_df$mean_p_adj <- rowMeans(plot_df[,c("p_adj_control", "p_adj_disease")])
  plot_df$min_p_adj <- apply(plot_df[,c("p_adj_control", "p_adj_disease")], 1, FUN = min)
  plot_df$label <- NA
  plot_df$label[plot_df$differential_p < sig_label_cutoff] <- plot_df$geneset[plot_df$differential_p < sig_label_cutoff] %>% tolower()
  plot_df$label <- gsub("dna_", "DNA_", plot_df$label)
  plot_df$label <- gsub("rna_", "RNA_", plot_df$label)
  plot_df$label <- gsub("ase_iii", "ase_III", plot_df$label)
  plot_df$label <- gsub("ase_ii", "ase_II", plot_df$label)
  plot_df$label <- gsub("ase_i", "ase_I", plot_df$label)
  plot_df$label <- ifelse(is.na(plot_df$label), NA,
                          paste0(toupper(substr(plot_df$label, 1, 1)), substr(plot_df$label, 2, nchar(plot_df$label))))
  
  plot_df$label <-
    sapply(plot_df$label, function(x) stringr::str_wrap(stringr::str_replace_all(x,"_", " "), width = 10))
  
  plot_df <- plot_df %>% mutate(col = ifelse(differential_p >= sig_label_cutoff, "n.s.", ifelse(neglog10_p_ratio>0, "Disease", "Control")))
  plot_df$col <- factor(plot_df$col, levels=c("n.s.", "Disease", "Control"))
  
  
  MA_plot <- ggplot(plot_df, aes(y=(-1)*log10(min_p_adj), x=neglog10_p_ratio)) +
    geom_point(alpha=0.5, aes(col = col, size=(-1)*log10(differential_p))) +
    # scale_x_continuous(limits=c(0,4), oob=scales::squish) +
    theme_bw() +
    geom_text_repel(aes(label=label)) +
    scale_colour_manual("Association tendency",
                        labels=c("n.s.", "Disease", "Control"), values=c("grey", "firebrick", "dodgerblue4")) +
    labs(size=bquote(-log[10](P[differential]))) +
    scale_x_continuous(bquote("Association ratio [" ~ -log[10](P[disease]/P[control]) ~ "]")) +
    scale_y_continuous(bquote("Conditional significance [" ~ -log[10] ~ "(min(" ~ P[adj] ~ "))]")) +
    theme(axis.title = element_text(size=14),
          axis.text = element_text(size=12))
  ggsave(paste0(plot_prefix, ".MA_plot.png"), plot=volcano_plot, width=80, height=50, units = "mm")
  ggsave(paste0(plot_prefix, ".MA_plot.svg"), plot=volcano_plot, width=80, height=50, units = "mm", device=svglite::svglite)
  
  
  
  
  
  plot_df$neglog10_diff_p <- (-1)*(log10(plot_df$differential_p))
  nudge_y_centre = diff(range(plot_df$neglog10_diff_p))/2
  nudge_y = diff(range(plot_df$neglog10_diff_p))/15
  nudge_x = diff(range(plot_df$logFC))/15
  
  # if (remove_outliers) {
  #   plot_df = plot_df[!plot_df$logFC %in% boxplot_stats(plot_df$neglog10_p_ratio)$out,]
  #   plot_df = plot_df[!plot_df$neglog10.adj %in% boxplot_stats(plot_df$neglog10.adj)$out,]
  # }
  
  
  volcano_plot <- ggplot(plot_df, aes(x = neglog10_p_ratio, y =(-1)*log10(differential_p),
                                      label = label)) +
    geom_point(alpha=0.5, size=2, aes(col = col)) +
    geom_text_repel(aes(x = neglog10_p_ratio,y = (-1)*log10(differential_p), label = label),
                    size = 5/.pt, col = "black",
                    segment.color = "black", segment.alpha = 0.4,
                    position = position_nudge_center(direction = "split"),
                    hjust = "outward", vjust = "outward",
                    lineheight=0.7) +
    scale_colour_manual("Association tendency",
                        labels=c("n.s.", "Disease", "Control"), values=c("grey", "firebrick", "dodgerblue4")) +
    theme_bw() +
    theme(text=element_text(size=7),
          axis.text = element_text(colour = 'black'),
          axis.line = element_line(),
          panel.border = element_blank(),
          axis.ticks = element_line(colour = 'black'),
          strip.background = element_rect(colour = NA)) +
    scale_y_continuous(bquote(-log[10](DifferentialP))) +
    scale_x_continuous(bquote(-log[10](P[disease]/P[control])))
  ggsave(paste0(plot_prefix, ".volcano_plot.png"), plot=volcano_plot, width=80, height=50, units = "mm")
  ggsave(paste0(plot_prefix, ".volcano_plot.svg"), plot=volcano_plot, width=80, height=50, units = "mm", device=svglite::svglite)
  
  return(list("MA_plot"=MA_plot, "volcano_plot"=volcano_plot))


}
















# # to do: add option to read from excel
# # to do: rewrite get_neg_log to include cases of several 0s
# # to do: add option to read gene set lists from
# # to do: check if get_best_predictors() works without the GOI, expr params




