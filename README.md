
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

<img align="center" width="50%" height="50%" src="https://github.com/si-ze/geneCOCOA/assets/129768077/616492c8-bc4f-41ae-b382-e8cc7a1a3bea">



### Side-by-side comparison of two GeneCOCOA results (same GOI, different conditions)
GeneCOCOA also offers the option to visually compare the GeneCOCOA results for one GOI in two different conditions (diverging bar plot). 

Therefore, we need to pass two results (produced by the `get_stats()` function). Please make sure that these were run with the **same GOI** and the **same geneset_collection**. (The function can only visualise two GeneCOCOA results, not check the validity of the comparison.)

For our example, we first need to generate GeneCOCOA results for the FH control samples as well.


``` 
# select only the columns containing control expression data
FH_control <- FH %>% select(contains("Control"))

# prepare expression data as GeneCOCOA input
control_input <- get_expr_info(expr=FH_control, GOI="LDLR")

# run GeneCOCOA
control_res <- get_stats(geneset_collection=hallmark_sets, 
                         GOI="LDLR", 
                         GOI_expr=control_input$GOI_expr,
                         expr_df=control_input$expr_df, 
                         samplesize=2, nsims=1000)
```

Next, we pass the two results to `plot_control_vs_treatment()`. 


```
plot_control_vs_treatment(
  control_res=control_res,
  treatment_res=control_res,
  treatment_label="FH",
  control_label="control",
  topN <- 10,
  sort_by="treament"
)
```
The "treatment" and "control" labels and the respective colours can be customised. The function returns the top n terms (customisable via the `topN` parameter), either with respect to the control or the treatment condition (`sort_by` parameter). 
The resulting plot is returned as a ggplot2 object by the function. If you would additionally like to save it directly, you can pass a location via the optional `filepath` parameter, e.g. `filepath="./LDLR.FH.diverging_bars.png"`.


<img align="center" width="60%" height="60%" src="https://github.com/si-ze/geneCOCOA/assets/129768077/48d9caaf-b6bd-4921-8712-1b0672f34739">

