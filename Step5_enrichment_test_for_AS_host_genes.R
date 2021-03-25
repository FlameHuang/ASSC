set.seed(1000)
### identify a funtion to perform KEGG pathway enrichment test in mice *****
run_KEGG_mice = function(candidate_gene, background_gene = NULL, gene_format, cutoff,
                    showCategory = 20, font.size = 10,title){
  diff_gene_ID<-clusterProfiler::bitr(candidate_gene, fromType = gene_format, toType ="ENTREZID", OrgDb="org.Mm.eg.db")

  if(is.null(background_gene)){
    ekegg<-enrichKEGG(gene=diff_gene_ID$ENTREZID, organism = 'mmu',
                      pAdjustMethod = "BH",
                      pvalueCutoff =cutoff,
                      qvalueCutoff = cutoff)
  }else{
    background_gene = clusterProfiler::bitr(background_gene, fromType = gene_format, toType ="ENTREZID", OrgDb="org.Mm.eg.db")
  ekegg<-enrichKEGG(gene=diff_gene_ID$ENTREZID, organism = 'mmu',
                    universe  = background_gene$ENTREZID,
                    pAdjustMethod = "BH",
                    pvalueCutoff =cutoff,
                    qvalueCutoff = cutoff)
  }
  #ekegg <- setReadable(ekegg, OrgDb = org.Mm.eg.db, keytype = "ENTREZID")
  ekegg.table = as.data.frame(ekegg)
  if(nrow(ekegg.table)>0){
    ekegg.table = ekegg.table[order(ekegg.table$p.adjust),]
    print(ekegg.table$Description)
    print(dotplot(ekegg, showCategory = showCategory, font.size = font.size, x='Count', title= title))
  }
  return(ekegg)
}


### identify a funtion to perform GO category enrichment test in mice ***
run_GO_mice = function(candidate_gene, background_gene=NULL, gene_format, ontology, cutoff,
                  showCategory,font.size,title){
  diff_gene_ID<-clusterProfiler::bitr(candidate_gene, fromType = gene_format, toType ="ENTREZID", OrgDb="org.Mm.eg.db")
  if(is.null(background_gene)){
    ego <-  simplify(enrichGO(gene = diff_gene_ID$ENTREZID,  OrgDb = org.Mm.eg.db,
                    keyType = 'ENTREZID', ont = ontology, readable = T,
                    pAdjustMethod = "BH", qvalueCutoff  = cutoff, pvalueCutoff  = cutoff))
  } else{
    background_gene = clusterProfiler::bitr(background_gene, fromType = gene_format, toType ="ENTREZID", OrgDb="org.Mm.eg.db")
    ego <-  simplify(enrichGO(gene = diff_gene_ID$ENTREZID,  OrgDb = org.Mm.eg.db,
                              universe = background_gene$ENTREZID,
                     keyType = 'ENTREZID', ont = ontology, readable = T,
                     pAdjustMethod = "BH", qvalueCutoff  = cutoff, pvalueCutoff  = cutoff))
  }

  if(nrow(ego@result)>0){
    print(dotplot(ego, showCategory = showCategory, font.size = font.size, x='Count',title=title))
  }
  return(ego)
}
