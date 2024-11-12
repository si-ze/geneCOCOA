# library(devtools)
# devtools::install_github("si-ze/geneCOCOA")
library(geneCOCOA)
library(gemma.R)
library(dplyr)


hallmark_sets <- get_msigdb_genesets("HALLMARK")
FH <- as.data.frame(gemma.R::get_dataset_expression("GSE6054"))


FH <- FH[nchar(FH$GeneSymbol)>0,]
FH$Probe = FH$GeneName = FH$NCBIid = NULL
FH$GeneSymbol <- sub("[|].*", "",FH$GeneSymbol)
FH <-  group_by(FH, GeneSymbol) %>%
  summarize_at(c(1:(ncol(FH)-1)), sum) %>% as.data.frame()
rownames(FH) = FH$GeneSymbol


FH_disease <- FH %>% select(contains("FH"))
FH_control <- FH %>% select(contains("Control"))

expr_info.LDLR.FH_disease <- get_expr_info(expr=FH_disease, GOI="LDLR")
res.LDLR.FH_disease <- get_stats(geneset_collection=hallmark_sets, GOI="LDLR", GOI_expr=expr_info.LDLR.FH_disease$GOI_expr, expr_df=expr_info.LDLR.FH_disease$expr_df, samplesize=2, nsims=1000)

expr_info.LDLR.FH_control <- get_expr_info(expr=FH_control, GOI="LDLR")
res.LDLR.FH_control <- get_stats(geneset_collection=hallmark_sets, GOI="LDLR", GOI_expr=expr_info.LDLR.FH_control$GOI_expr, expr_df=expr_info.LDLR.FH_control$expr_df, samplesize=2, nsims=1000)



test_diff <- get_differential_GeneCOCOA_results(res.LDLR.FH_disease, res.LDLR.FH_control, GOI="LDLR", geneset_collection="HALLMARK")
