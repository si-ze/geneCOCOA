
# GeneCOCOA

GeneCOCOA (**co**mparative **co**-expression **a**naylsis focussed on a **gene** of interest) has been developed as a method to integrate prior knowledge (e.g. Gene Ontology, MSigDB genesets, ...) with co-expression analysis while focussing on one gene-of-interest (GOI). 

After providing expression data from a specific experiment and defining a GOI, GeneCOCOA statistically ranks different gene sets/terms by strength of association with the GOI. GeneCOCOA thereby allows for functional annotation of experiment-specific expression data in a gene-centric manner. 



## Installation 


``` 
library(devtools)
devtools::install_github("si-ze/geneCOCOA")
``` 

## Usage
### Exploring curated data sets

In an example use case, we will explore the co-expression signature of the gene coding for the low-density lipoprotein receptor (*LDLR*) in familial hypercholesterolemia (FH). In this data set, we inspect co-expression patterns of *LDLR* with the [50 MSigDB Hallmark gene sets](https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp). 

``` 
# load MSigDB Hallmark gene sets 
hallmark_sets <- get_msigdb_genesets("HALLMARK")
``` 

GeneCOCOA offers the option to load public data sets using [gemma.R](https://github.com/PavlidisLab/gemma.R). We load the FH dataset with the GEO accession GSE6054. 


``` 
# load complete FH data set
FH <- as.data.frame(get_dataset_expression("GSE6054"))

# preprocessing:
# remove rows without gene symbols
FH <- FH[nchar(FH$GeneSymbol)>0,]
# remove columns not needed (we only want to keep the columns holding expression info)
FH$Probe = FH$GeneName = FH$NCBIid = NULL
# keep only official gene symbol / remove aliases 
# (e.g. "LDLR|LDLCQ2|Low-Density Lipoprotein Receptor" becomes "LDLR")
FH$GeneSymbol <- sub("[|].*", "",FH$GeneSymbol)
# if there are multiple rows per gene symbol, collapse
FH <-  group_by(FH, GeneSymbol) %>%
  summarize_at(c(1:(ncol(FH)-1)), sum) %>% as.data.frame()
# name rows by gene symbol
rownames(FH) = FH$GeneSymbol


``` 

Now, let's inspect the co-expression patterns of *LDLR* in only the patient samples (excluding healthy controls)


``` 
# select only the columns containing FH patient expression data
FH_disease <- FH %>% select(contains("FH"))

# prepare expression data as GeneCOCOA input
disease_input <- get_expr_info(expr=FH_disease, GOI="LDLR")

# run GeneCOCOA
disease_res <- get_stats(geneset_collection=hallmark_sets, 
                         GOI="LDLR", 
                         GOI_expr=disease_input$GOI_expr,
                         expr_df=disease_input$expr_df, 
                         samplesize=2, nsims=1000)

# plot GeneCOCOA results
plot_volcano(mystats=disease_res)
```

The resulting plot shows the -log10(*P*adj) plotted against the "direction" of co-expression (stronger or weaker than mean co-expression of all genes with *LDLR*). The size of the points in each plot reflects the relative mean expression level of each gene set.
The plot is returned as a ggplot2 object by the function. If you would additionally like to save it directly, you can pass a location via the optional `filepath` parameter, e.g. `filepath="./LDLR.FH.volcano.png"`.

<img align="center" width="45%" height="45%" src="https://github.com/si-ze/geneCOCOA/assets/129768077/616492c8-bc4f-41ae-b382-e8cc7a1a3bea">



### Differential GeneCOCOA (same GOI, same gene set collection, different conditions)
GeneCOCOA also offers the option to differentially test the association of one GOI with a collection of gene sets in two conditions.  

Therefore, we need to pass two results (produced by the `get_stats()` function). Please make sure that these were run with the **same GOI** and the **same geneset_collection**. 

For our example, we first need to generate GeneCOCOA results for both the FH disease and the FH control samples. 


``` 
library(devtools)
devtools::install_github("si-ze/geneCOCOA", lib="/mnt/f/Rlib")
library(geneCOCOA)
library(gemma.R)
library(dplyr)


# get data on familial hypercholesterolemia
FH <- as.data.frame(gemma.R::get_dataset_expression("GSE6054"))
FH <- FH[nchar(FH$GeneSymbol)>0,]
FH$Probe = FH$GeneName = FH$NCBIid = NULL
FH$GeneSymbol <- sub("[|].*", "",FH$GeneSymbol)
# summarise expression values on gene level
FH <-  group_by(FH, GeneSymbol) %>%
  summarize_at(c(1:(ncol(FH)-1)), sum) %>% as.data.frame()
# assign gene symbols as row names
rownames(FH) = FH$GeneSymbol

# get disease expression matrix (with gene symbols as row names)
FH_disease <- FH %>% select(contains("FH"))
# get control expression matrix (with gene symbols as row names)
FH_control <- FH %>% select(contains("Control"))


# define list of gene sets to be differentially analysed regarding their association with the GOI
hallmark_sets <- get_msigdb_genesets("HALLMARK")

# prepare input expression info for GeneCOCOA
expr_info.LDLR.FH_disease <- get_expr_info(expr=FH_disease, GOI="LDLR")
# run GeneCOCOA on control data set to get Condition P-values (indicating significance of association between GOI and gene set in disease)
res.LDLR.FH_disease <- get_stats(geneset_collection=hallmark_sets, GOI="LDLR", GOI_expr=expr_info.LDLR.FH_disease$GOI_expr, expr_df=expr_info.LDLR.FH_disease$expr_df, samplesize=2, nsims=1000)

# prepare input expression info for GeneCOCOA
expr_info.LDLR.FH_control <- get_expr_info(expr=FH_control, GOI="LDLR")
# run GeneCOCOA on control data set to get Condition P-values (indicating significance of association between GOI and gene set in control)
res.LDLR.FH_control <- get_stats(geneset_collection=hallmark_sets, GOI="LDLR", GOI_expr=expr_info.LDLR.FH_control$GOI_expr, expr_df=expr_info.LDLR.FH_control$expr_df, samplesize=2, nsims=1000)
```

Next, we pass the two results to the differential function. 

```
# feed the two Condition results into differential function 
differential_results <- get_differential_results(res.LDLR.FH_disease, res.LDLR.FH_control, 
                                      GOI="LDLR", geneset_collection="HALLMARK", 
                                      laplace_parameters = "Default")
```

This function returns a table: 
```
# > head(differential_results)
#               geneset    p_control p_adj_control    p_disease p_adj_disease
# 1        ADIPOGENESIS 4.912080e-03  2.824446e-02 8.827597e-01  1.000000e+00
# 2 ALLOGRAFT_REJECTION 9.613390e-01  9.999999e-01 2.021361e-01  4.893822e-01
# 3   ANDROGEN_RESPONSE 7.491303e-06  6.891999e-05 1.956240e-08  1.799740e-07
# 4        ANGIOGENESIS 9.966363e-01  9.966363e-01 5.295860e-01  9.666187e-01
# 5     APICAL_JUNCTION 4.700123e-02  1.801714e-01 1.495597e-08  1.719937e-07
# 6      APICAL_SURFACE 7.594417e-01  9.999999e-01 4.276494e-01  8.941761e-01
#        p_ratio log10_p_ratio neglog10_p_ratio differential_p
# 1 1.797120e+02     2.2545771       -2.2545771     0.25860745
# 2 2.102652e-01    -0.6772327        0.6772327     0.66827842
# 3 2.611348e-03    -2.5831353        2.5831353     0.21220528
# 4 5.313734e-01    -0.2746002        0.2746002     0.85154028
# 5 3.182039e-07    -6.4972946        6.4972946     0.02011863
# 6 5.631102e-01    -0.2494066        0.2494066     0.86455135
```

To visualise this, we use the following function, which outputs the results to images in the current working directory (but we can specify the output location via the  `output_dir` parameter - check the help page of `plot_differential_results` for usage. 


<img align="center" width="60%" height="60%" src="https://github.com/user-attachments/assets/3cfacb33-9323-4697-a2ef-890dc2b5276f">






