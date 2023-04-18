#' @rawNamespace import(data.table, except = c(last, first, between, transpose))
#' @rawNamespace import(psych, except=c(alpha, rescale, "%+%"))
#' @rawNamespace import(scales, except=c(discard, alpha))
#' @import dplyr purrr
#' @import ggplot2 ggrepel ggpp stringr
#'
NULL


#
#' Main function of Gene-COCOA: performs statistical analysis.
#'
#' @param geneset_list A list of lists. Each list holds gene symbols of a specific gene set. Can be generated with \code{\link{get_msigdb_genesets}}
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
#' res <- get_stats(geneset_list=get_msigigdb_genesets("HALLMARK"), GOI="PDGFD", GOI_expr=expr_info$GOI_expr, expr_df=expr_info$expr_df)
#'}
get_stats <- function(geneset_list,
                      GOI,
                      GOI_expr,
                      expr_df,
                      samplesize=ncol((expr_df)/10),
                      nsims=1000) {


  p_value_list <- list()
  low_power.p_value_list <- list()

  all_RMSEs.real <- list()
  all_RMSEs.random <- list()
  geom_mean_cor.real <- list()
  geom_mean_cor.random <- list()

  intercepts <- list()
  mean_cors <- list()

  geom_mean_expr <- list()

  for (i in 1:length(geneset_list)) {


    my_genes_in_set <- as.vector(unlist(geneset_list[i]))
    my_genes_in_set <- gsub("-", "_", my_genes_in_set)
    tpm.genes_in_set <- expr_df[my_genes_in_set,]
    t.tpm.genes_in_set <- as.data.frame(t(tpm.genes_in_set))
    rownames(t.tpm.genes_in_set) <- colnames(tpm.genes_in_set)
    colnames(t.tpm.genes_in_set) <- rownames(tpm.genes_in_set)
    names(t.tpm.genes_in_set) <- rownames(tpm.genes_in_set)
    t.tpm.genes_in_set <- t.tpm.genes_in_set %>% purrr::discard(~all(is.na(.)))
    t.tpm.genes_in_set$GOI <- NULL


    my_set_name <- names(geneset_list[i]) # get name of gene set
    my_set_size <- ncol(t.tpm.genes_in_set)

    intercepts[[my_set_name]] <- get_intercept(GOI, GOI_expr, t.tpm.genes_in_set)
    mean_cors[[my_set_name]] <- get_mean_cor(GOI, GOI_expr, t.tpm.genes_in_set)
    geom_mean_cor.real[[my_set_name]] <- get_geom_mean_cor(GOI, GOI_expr, t.tpm.genes_in_set)

    tpm.my_rest_genes <- expr_df[!(rownames(expr_df) %in% my_genes_in_set),]
    t.tpm.my_rest_genes <- as.data.frame(t(tpm.my_rest_genes))
    rownames(t.tpm.my_rest_genes) <- colnames(tpm.my_rest_genes)
    colnames(t.tpm.my_rest_genes) <- rownames(tpm.my_rest_genes)
    names(t.tpm.my_rest_genes) <- rownames(tpm.my_rest_genes)
    t.tpm.my_rest_genes <- t.tpm.my_rest_genes %>% purrr::discard(~all(is.na(.)))
    t.tpm.my_rest_genes$GOI <- NULL




    # a gene set size of 17 is large enough to provide for more than 1000 unique samples
    # see binomial coefficient: while choose(16,13)=560, choose(17,13)=2380
    if (my_set_size>16) {
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
      my_sample <- cbind.data.frame(GOI=GOI_expr,
                                    t.tpm.genes_in_set)
      my_model <- stats::lm(paste(GOI, "~ ."),data=my_sample)
      real_RMSE <- rmse(my_model)

      my_rand <- compute_x_RMSEs(GOI,
                                 GOI_expr,
                                 my_t.tpms=t.tpm.my_rest_genes,
                                 samplesize=samplesize, nsims=nsims, geom_mean_cors = TRUE)
      random_RMSEs <- my_rand$RMSEs
      geom_mean_cors <- my_rand$geom_mean_cors

      low_power.p <- t.test(unlist(real_RMSEs),
                            unlist(random_RMSEs),
                            paired=FALSE, alternative="less")
      low_power.p_value_list[[paste(my_set_name,
                                    " (", my_set_size, ")",
                                    sep="")]] <- low_power.p$p.value


      all_RMSEs.real[[my_set_name]] <- real_RMSE
      all_RMSEs.random[[my_set_name]] <- random_RMSEs
    }


    geom_mean_cors <- as.vector(unlist(geom_mean_cors))
    geom_mean_cors <- geom_mean_cors[!is.na(geom_mean_cors) & !is.infinite(geom_mean_cors)]
    geom_mean_cor.random[[my_set_name]] <- mean(geom_mean_cors)


    geom_mean_expr[[my_set_name]] <- get_geom_mean_expr(t.tpm.genes_in_set)
  }


  p_value_df <- data.frame("geneset"=names(p_value_list),
                           "p"=unlist(p_value_list),
                           "neglog10"=get_neg_log10(unlist(p_value_list)))
  p_value_df$p.adj <- p.adjust(p_value_df$p, method = "BH")
  p_value_df$neglog10.adj <- get_neg_log10(unlist(p_value_df$p.adj))

  low_power.p_value_df <- data.frame("geneset"=names(low_power.p_value_list),
                                     "p"=unlist(low_power.p_value_list),
                                     "neglog10"=get_neg_log10(unlist(low_power.p_value_list)))

  # transform "intercepts" from named list to dataframe
  intercepts <- stack(intercepts) %>% rename("intercept"="values", "geneset"="ind")
  geom_mean_expr.df <- stack(geom_mean_expr) %>% rename("geom_mean_expr"="values", "geneset"="ind")
  mean_cors <- stack(mean_cors) %>% rename("mean_cor"="values", "geneset"="ind")

  geom_mean_cor.real.df <- stack(geom_mean_cor.real) %>% rename("GOI_cor.geneset"="values", "geneset"="ind")
  geom_mean_cor.random.df <- stack(geom_mean_cor.random) %>% rename("GOI_cor.background"="values", "geneset"="ind")


  geom_mean_cor.df <- merge(geom_mean_cor.real.df,  geom_mean_cor.random.df, by="geneset")
  geom_mean_cor.df$logFC <- log2(geom_mean_cor.df$GOI_cor.geneset/geom_mean_cor.df$GOI_cor.background)



  return(list("p_value_df"=p_value_df,
              "low_power.p_value_df"= low_power.p_value_df,
              "all_RMSEs.real"=all_RMSEs.real,
              "all_RMSEs.random"=all_RMSEs.random,
              "intercepts"=intercepts,
              "mean_cor"=mean_cors,
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


#' Returns -log10 of given p-value list
#'
#' @param p_value_list
#'
#' @return
#' @keywords Internal
#'
get_neg_log10 <- function(p_value_list) {
  # if the gene set is too small to rescale
  if (length(p_value_list)<2) {
    return(p_value_list)
  } else {
    # if there is a p=0, rescale the p value list to [min!=0,max]
    # else, keep as is
    pvalues_for_neglog <- if (min(p_value_list)==0) scales::rescale(
      p_value_list, to=c(sort(p_value_list, FALSE)[2], max(p_value_list))) else p_value_list

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


#' Recommended visualisation of Gene-COCOA results.
#' Produces a ggplot with x=(GOI coexpression with gene set)/(GOI coexpression with random set), y=-log10(p.adj).
#' Dot size and alpha reflect the (geometric) mean expression of all genes in a respective gene set.
#'
#' @param mystats Output of \code{\link{get_stats}}.
#' @param sig_label_cutoff At which p.adj-level should points be labelled? Default=0.01.
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
#' res <- get_stats(geneset_list=get_msigigdb_genesets("HALLMARK"), GOI="PDGFD", GOI_expr=expr_info$GOI_expr, expr_df=expr_info$expr_df)
#' # show GeneCOCOA ggplot in "Plots" window
#' plot_volcano(res)
#'
#' # assign GeneCOCOA ggplot to variable
#' p <- plot_volcano(res)
#'
#' # assign GeneCOCOA ggplot to variable AND save GeneCOCOA ggplot to file
#' p <- plot_volcano(res, filepath="/path/to/my/GeneCOCOA.volcano.png")
#'}
plot_volcano <- function(mystats, sig_label_cutoff=0.01, sig_colour="dodgerblue4", remove_outliers=FALSE, filepath="") {

  new_df <- merge(mystats$p_value_df, mystats$geom_mean_cor.df, by="geneset")
  new_df <- merge(new_df, mystats$geom_mean_expr.df, by="geneset")
  new_df$label <- NA
  new_df$label[new_df$p.adj<sig_label_cutoff] <- sapply(new_df$geneset[new_df$p.adj<sig_label_cutoff],  function(x) stringr::str_wrap(stringr::str_replace_all(x, "_", " "),width = 10))
  new_df$col <- "grey"
  new_df$col[new_df$p.adj<sig_label_cutoff] <- sig_colour
  new_df$alpha <- 0.7
  new_df$alpha[new_df$p.adj<sig_label_cutoff] <- 1
  new_df$alpha <- scales::rescale(new_df$geom_mean_expr, to=c(0,0.75))

  if (remove_outliers) {
    new_df = new_df[!new_df$logFC %in% boxplot_stats(new_df$logFC)$out,]
    new_df = new_df[!new_df$neglog10.adj %in% boxplot_stats(new_df$neglog10.adj)$out,]
  }

  nudge_y_centre=diff(range(new_df$neglog10.adj))/2
  nudge_y=diff(range(new_df$neglog10.adj))/15
  nudge_x=diff(range(new_df$logFC))/15

  my_plot <- ggplot(new_df, aes(x=logFC, y=neglog10.adj, label=label)) +
    geom_point(alpha=new_df$alpha, col=new_df$col, size=new_df$geom_mean_expr) +
    geom_text_repel(aes( x=logFC, y=neglog10.adj, label=label),
                    size=3,
                    col="black",
                    #min.segment.length = unit(0, 'lines'),
                    segment.color=sig_colour,
                    segment.alpha=0.4,
                    #nudge_y = ifelse(new_df$neglog10.adj >nudge_y_centre, nudge_y, -nudge_y),
                    position=position_nudge_center(#x = nudge_x,y = nudge_y,
                      direction = "split"),
                    hjust = "outward", vjust = "outward"

    ) +
    theme_minimal(base_size=14) +
    scale_y_continuous("-log10(p.adj)") +
    scale_x_continuous("log2(diffcoexp)", limits=c(1.25*(-max(abs(new_df$logFC))), 1.25*max(abs(new_df$logFC))))

  if (nchar(filepath)>0) {
    ggsave(filename=filepath, plot=my_plot, dpi=300)
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
#' res <- get_stats(geneset_list=get_msigigdb_genesets("HALLMARK"), GOI="PDGFD", GOI_expr=expr_info$GOI_expr, expr_df=expr_info$expr_df)
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



#
#' Function to plot multiple (long) p-value rankings. Obsolete.
#'
#'
#' @keywords Internal
#'
plot_multiple_p_plots <- function(subdf) ggplot(subdf, aes(x=reorder(geneset, p), y=p, group=1,
                                                           label=formatC(p, digits = 6)))+
  scale_x_discrete(limits=rev, name="") +
  scale_y_continuous(name ="p-value", limits=c(0,1.35)) +
  geom_line() + geom_point() +
  geom_text(size=6, hjust=-0.125, vjust=0.5) +
  coord_flip() +
  theme(
    plot_title = element_text(size=21, hjust = 0.5),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 21),
    legend.title =element_text(size = 18),
    legend.text =  element_text(size = 18)

  )


#' Plots p-value ranking (x-axis=adjusted -log10(p), y= gene set names, ranked).
#' This function works but \code{\link{plot_volcano}} is more informative.
#' @param my_stats Output of \code{\link{get_stats}}
#' @param filepath Path + name of png to save If not provided, only the ggplot will be returned by this function.
#'
#' @return Returns a ggplot which can be stored in a variable and modified. If an absolute file name (path + file name) is provided, the ggplot with additionally be saved to a png.
#' @export
#'
#'
plot_basic_p_ranking <- function(my_stats, filepath)  {


  ggplot(my_stats$p_value_df %>%
           arrange(neglog10.adj, rev(geneset)) %>%
           mutate(geneset = factor(geneset, levels = geneset)),
         aes(x=geneset,
             y=neglog10.adj, group=1, label=formatC(p.adj, digits = 6))) +
    scale_x_discrete(name="gene set") +
    scale_y_continuous(name ="-log10(p.adj)") +
    geom_point()  + geom_line() +
    ggthemes::theme_stata() +
    geom_text(size=6,  hjust="inward", vjust="center") +
    coord_flip() +
    theme(
      plot_title = element_text(size=21, hjust = 0.5),
      axis.text = element_text(size = 18),
      axis.text.y = element_text(angle = 0), # to have the labels (gene set names) horizontally
      axis.title = element_text(size = 21),
      legend.title =element_text(size = 18),
      legend.text =  element_text(size = 18)

    )
  ggsave(filename=filepath,
         path=getwd(),
         width=10, height=25, units="in", dpi=300)
}





#' Returns a ggplot with x=-log10(p.adj), y=mean intercept of gene set.
#' The intercept can be interpreted as "GOI expression if all genes in the gene set were not expressed".
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
plot_p_vs_intercept <- function(mystats, sig_label_cutoff=0.05, sig_colour="dodgerblue4", remove_outliers=FALSE, filepath) {

  new_df <- merge(mystats$p_value_df, mystats$intercepts, by="geneset")
  new_df$label <- NA
  new_df$label[new_df$p.adj<sig_label_cutoff] <- sapply(new_df$geneset[new_df$p.adj<sig_label_cutoff],  function(x) stringr::str_wrap(stringr::str_replace_all(x, "_", " "),width = 10))
  new_df$col <- "grey"
  new_df$col[new_df$p.adj<sig_label_cutoff] <- sig_colour
  new_df$alpha <- 0.7
  new_df$alpha[new_df$p.adj<sig_label_cutoff] <- 1


  if (remove_outliers) {
    new_df = new_df[!new_df$intercept %in% boxplot_stats(new_df$intercept)$out,]
  }

  nudge_y_centre=diff(range(new_df$intercept))/2
  nudge_y=diff(range(new_df$intercept))/5


  my_plot <- ggplot(new_df, aes(x=neglog10.adj, y=intercept, label=label)) +
    geom_point(alpha=new_df$alpha, col=new_df$col) +
    geom_text_repel(aes(x=neglog10.adj, y=intercept, label=label),
                    col=sig_colour,
                    min.segment.length = unit(0, 'lines'),
                    segment.color=sig_colour,
                    segment.alpha=0.5,
                    nudge_y = ifelse(new_df$intercept >nudge_y_centre, nudge_y, -nudge_y)) +
    theme_minimal(base_size=14) +
    scale_y_continuous("Intercept of regression model predicting GOI expression") +
    scale_x_continuous("Significane of coexpression with GOI [-log10(p.adj)]")

  if (nchar(filepath)>0) {
    ggsave(filename=filepath, plot=my_plot, dpi=300)
  }

  return(my_plot)

  # how to interpret the intercept:
  # the GOIs tpm equals the intercept value when none of the genes in the set are expressed

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

  new_df <- merge(mystats$p_value_df, mystats$geom_mean_cor.df, by="geneset")
  new_df$label <- NA
  new_df$label[new_df$p.adj<sig_label_cutoff] <- sapply(new_df$geneset[new_df$p.adj<sig_label_cutoff],  function(x) stringr::str_wrap(stringr::str_replace_all(x, "_", " "),width = 10))
  new_df$col <- "grey"
  new_df$col[new_df$p.adj<sig_label_cutoff] <- sig_colour
  new_df$alpha <- 0.7
  new_df$alpha[new_df$p.adj<sig_label_cutoff] <- 1

  if (remove_outliers) {
    new_df = new_df[!new_df$GOI_cor.geneset %in% boxplot_stats(new_df$GOI_cor.geneset)$out,]
  }

  nudge_y_centre=diff(range(new_df$GOI_cor.geneset))/2
  nudge_y=diff(range(new_df$GOI_cor.geneset))/5

  my_plot <- ggplot(new_df, aes(x=neglog10.adj, y=GOI_cor.geneset, label=label)) +
    geom_point(alpha=new_df$alpha, col=new_df$col) +
    geom_text_repel(aes(x=neglog10.adj, y=GOI_cor.geneset, label=label),
                    col=sig_colour,
                    min.segment.length = unit(0, 'lines'),
                    segment.color=sig_colour,
                    segment.alpha=0.5,
                    nudge_y = ifelse(new_df$GOI_cor.geneset >nudge_y_centre, nudge_y, -nudge_y)) +
    theme_minimal(base_size=14) +
    scale_y_continuous("Mean coexpression with GOI") +
    scale_x_continuous("Significane of coexpression with GOI [-log10(p.adj)]")

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
#' res <- get_stats(geneset_list=all_hallmark_gs, GOI="PDGFD",  ...)
#'
#' # call directly in get_stats()
#' res <- get_stats(geneset_list=get_msigdb_genesets("HALLMARK"), ...)
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






# to do: add option to read from excel
# to do: rewrite get_neg_log to include cases of several 0s
# to do: add option to read gene set lists from
# to do: check if get_best_predictors() works without the GOI, expr params




