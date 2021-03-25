#### Alternative splicing (AS) site find ---------------------------
##### ********** we modify the DSU pipeline to optimize the identification of AFE ************
### detect the alternative splicing junctions (shared same start or same end) ------
ASSite <- function(junc, NT = 4){
  # NT: the number of threads to use
  # required package ------------------------
  usePackage("data.table")
  usePackage("parallel")
  usePackage("plyr")

  # Identify alternative splicing events ----
  juncs <- junc$junctions
  junc.names <- data.table(chr = as.character(mapply(function(x) x[1], strsplit(juncs, split = ":"))),
                           start = as.integer(mapply(function(x) x[2], strsplit(juncs, split = "[:-]"))),
                           end = as.integer(mapply(function(x) x[3], strsplit(juncs, split = "[:-]"))),
                           strand = as.character(mapply(function(x) x[3], strsplit(juncs, split = ":"))))
  #rownames(junc.names) <- juncs
  junc.names$names <- juncs
  junc.names <- junc.names[order(chr, start, end),][chr%in%c(1:22,'X','Y')]
  #### add the exon-intron boundary site (different for + and - strand)
  #### for - strand (end_base)
  junc.names[, first_boundary:=paste0(chr,':',end,':-')]
  ### using boundary:= failed !
  ##### for + strand (start base)
  junc.names[strand=='+'] = transform(junc.names[strand=='+'], first_boundary=paste0(chr,':',start,':+'))
  junc.names <- as.data.frame(junc.names)

  junc.as <- do.call(c, mclapply(unique(junc.names$chr), function(chr) {
    junc.chr <- junc.names[junc.names$chr == chr, ]
    #### notice: when applu the dlply to data.table, will return eorror results !!!!
    same.start <- plyr::dlply(junc.chr, c("start"), function(x) x)
    same.start <- lapply(same.start,function(x){ x$share = 'same_start'; return(x)})
    same.end <- plyr::dlply(junc.chr, c("end"), function(x) {
      if(nrow(x)>1) {
        return(x)
      } else {
        return(NULL)
      }
    })
    ### remove the NULL values
    same.end = same.end[!sapply(same.end, is.null)]
    same.end <- lapply(same.end,function(x){ x$share = 'same_end'; return(x)})
    return(c(same.start, same.end))
  }, mc.cores = NT))
  ## junction just has two alternative starts/ends:
  junc.as <- junc.as[sapply(junc.as, nrow) > 1]
  names(junc.as) <- sapply(junc.as, function(x) paste(x$names, collapse="@"))
  return(junc.as)
}
##### create indexed gtf (extract the first exon/intron) ------
GTFIndex <- function(gtfFile, NT = 1){
  if(!file.exists(gtfFile)) {
    stop("The GTF file does not exist! And you should provide the GTF file you used in STAR alignment.")
  }

  # required package ------------------------
  usePackage("data.table")
  usePackage("parallel")
  useBiocPackage('GenomicFeatures')
  #### contrust the database from GTF -----------
  txdb <- GenomicFeatures::makeTxDbFromGFF(file = gtfFile, format = "gtf")
  #### the intron regions to map gene id -------------
  intron <- intronicParts(txdb, linked.to.single.gene.only = TRUE)
  intron_tab <- as.data.table(intron)
  intron_tab[, names:=as.character(intron)]
  # AF/LE -----------------------------------
  txintron <- intronsByTranscript(txdb, use.names=TRUE)
  txintron <- txintron[S4Vectors::elementNROWS(txintron) > 0]
  s_s <- as.character(unlist(unique(strand(txintron))))
  s_l <- S4Vectors::elementNROWS(txintron)
  txintron_gr <- unlist(txintron)
  my.rank <- function(s, l) ifelse(s == "+", return(1:l), return(-(l:1)))
  GenomicRanges::mcols(txintron_gr)$intron_rank <- do.call(c, lapply(seq_along(s_s), function(i) my.rank(s_s[i], s_l[i])))
  GenomicRanges::mcols(txintron_gr)$tx_id <- names(txintron_gr)
  txintron_tab_rank <- as.data.table(txintron_gr)

  first_intron_tab <- txintron_tab_rank[abs(intron_rank) == 1, ]
  first_intron_tab[, junc:=paste0(seqnames, ":", start, "-", end, ":", strand)]
  setkey(first_intron_tab, tx_id)

  last_intron_tab <- txintron_tab_rank[ , .SD[abs(intron_rank) == max(abs(intron_rank)), ], by = "tx_id"]
  last_intron_tab[, junc:=paste0(seqnames, ":", start, "-", end, ":", strand)]
  setkey(last_intron_tab, tx_id)

  # ASS -------------------------------------
  exons<- as.data.table(unique(exons(txdb, columns = c("gene_id", "EXONNAME"))))
  ### find first exon
  exon_tab <- as.data.table(unique(exonicParts(txdb)))
  #exon_tab <- exon_tab[seqnames%in%c(1:22,'X','Y')][,exon_rank:=as.numeric(as.character(exon_rank))]
  exon_tab[,exon_rank:=as.numeric(mapply(function(x) x[1], strsplit(as.character(exon_rank),' ')))]
  #exon_tab[,'exon_rank']=as.numeric(I(exon_tab[,'exon_rank']))
  first_exon_tab <- exon_tab[abs(exon_rank) == 1, ][,gene_id:=as.character(gene_id)][order(gene_id)]
  first_exon_tab = transform(first_exon_tab, seqnames = as.character(seqnames), start=as.numeric(start), end = as.numeric(end))
  ### obtain the first exon boundary site
  ### end+1 for + strand
  first_exon_tab[,exon_boundary:=paste0(seqnames,':',end+1,':+')]
  ### start-1 for - strand
  first_exon_tab[strand=='-'] = transform(first_exon_tab[strand=='-'], exon_boundary=paste0(seqnames,':',start-1,':-'))
  ### we retain the gene with multiple transcripts
  retained_gene = first_exon_tab[duplicated(gene_id)]$gene_id
  first_exon_tab = first_exon_tab[gene_id%in%retained_gene]
  setkey(first_exon_tab, gene_id)
  ### classify the distinct exons and
  ### tandem exons in + strand (shared same start, with different ends)
  #first_exon_tab_pos = first_exon_tab[strand=='+',.(startl = length(start), startu = length(unique(start)), endl=length(end), endu = length(unique(end))),by=gene_id]
  #first_exon_tab_pos[,exon_type:=ifelse(startl == startu & endu < endl, 'Tandem_exon', 'Distinct_exon')]
  ### tandem exons in - strand (shared same end, with different starts)
  #first_exon_tab_neg = first_exon_tab[strand=='-',.(startl = length(start), startu = length(unique(start)), endl=length(end), endu = length(unique(end))),by=gene_id]
 # first_exon_tab_pos[,exon_type:=ifelse(startu < startl & endu == endl, 'Tandem_exon', 'Distinct_exon')]

  return(list(intron_tab = intron_tab, first_intron_tab = first_intron_tab, last_intron_tab = last_intron_tab,
              exons = exon_tab, first_exon_tab = first_exon_tab))
}
