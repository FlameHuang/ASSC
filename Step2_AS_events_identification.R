### using GTF as reference index
GTFIndex <- function(gtfFile, NT = 1){
  if(!file.exists(gtfFile)) {
    stop("The GTF file does not exist! And you should provide the GTF file you used in STAR alignment.")
  }

  # required package ------------------------
  usePackage("data.table")
  usePackage("parallel")
  useBiocPackage('GenomicFeatures')

  # AF/LE -----------------------------------
  txdb <- GenomicFeatures::makeTxDbFromGFF(file = gtfFile, format = "gtf")
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
  exons <- as.data.table(unique(exons(txdb, columns = c("gene_id", "EXONNAME"))))

  return(list(first_intron_tab = first_intron_tab, last_intron_tab = last_intron_tab, exons = exons))
}


#### Classifying all AS type based on filtered AS sites ----------------
#### optimization of identifying alternative the first exon (AFE) ---------
#tandem first exon (with different length)
#A5SS & AFE in + strand
#A3SS & AFE in - strand
####The canonical AFE:
#### the alternative exon are located in different transcripts from the sam gene
#### map the first exon's end with junction (+);  the first exon's start with junction (-)
identify_AStype <- function(cell.psi, sj_list, indexedGTF, NT = 1){
  ### cell_psi: the result table from betabinomial test
  ### sj_list: the list of AS junctions
  ### indexedGTF: index files' list from GTF
  # required package ------------------------
  usePackage("data.table")
  usePackage("parallel")

  junc_as2 <- as.list(cell.psi[event == 2, id])
  junc_as <- do.call(rbind, lapply(junc_as2, function(x) {
    unlist(strsplit(stringr::str_replace_all(x, "-(?=[[:digit:]])", "@"), "[:@]"))
  }))
  junc_as <- as.data.table(junc_as)
  colnames(junc_as) <- paste0(c("chr", "start", "end", "strand"), rep(1:2, c(4,4)))
  junc_as[, names:=cell.psi[event == 2, as.character(id)]]

  same_start_junc2 <- junc_as[start1 == start2, ]
  same_end_junc2 <- junc_as[end1 == end2, ]
  same_start_junc2 <- same_start_junc2[, .(names, chr1, start1, end1, end2, strand1)]
  colnames(same_start_junc2) <- c("names", "chr", "start", "end1", "end2", "strand")
  same_end_junc2 <- same_end_junc2[, .(names, chr1, start1, start2, end1, strand1)]
  colnames(same_end_junc2) <- c("names", "chr", "start1", "start2", "end", "strand")

  match_nage <- same_start_junc2[strand=="-", same_end_junc2[strand=="-",], by = c("names", "chr", "start", "end1", "end2", "strand")]
  match_posi <- same_start_junc2[strand=="+", same_end_junc2[strand=="+",], by = c("names", "chr", "start", "end1", "end2", "strand")]
  start_end_match <- rbind(match_posi, match_nage)
  start_end_match <- start_end_match[as.logical(start_end_match[,2]==start_end_match[,8]), ]

  # SE --------------------------------------
  SE <- start_end_match[start==start2 & end==end2 & end1 < start1,]
  SE_junc <- unique(c(SE[[1]], SE[[7]]))

  # MXE -------------------------------------
  MXE <- start_end_match[start < end1 & end1 < start2 & start2 < end2 & end2 < start1 & start1 < end, ]
  MXE_junc <- unique(c(MXE[[1]], MXE[[7]]))

  # AF/LE (special we need use sj_list)-----------------------------------
  first_intron_tab <- indexedGTF$first_intron_tab
  last_intron_tab <- indexedGTF$last_intron_tab
  first_exon_tab <- indexedGTF$first_exon_tab
  afle_index <- mcmapply(sj_list, FUN = function(x){
    if(sum(x$names %in% first_intron_tab[, junc]) > 0 ){
      res = "Pseudo_AFE"
      if(sum(x$first_boundary %in% first_exon_tab$exon_boundary)>1){
        if(x$strand == '+' & x$share == 'same_end'){
          res = 'Real_AFE'
        }
        if(x$strand == '-' & x$share == 'same_start') {
          res = 'Real_AFE'
        }
      }
    } else if (sum(x$names %in% last_intron_tab[, junc]) > 0 ) {
      res = "ALE"
    } else {
      res = "Non"
    }
    return(res)
  }, mc.cores = NT)
  #AFE_junc <- cell.psi[, as.character(id)][grep('AFE',afle_index)]
  R_AFE_junc <- cell.psi[id%in%names(afle_index[afle_index == "Real_AFE"]), as.character(id)]
  P_AFE_junc <- cell.psi[id%in%names(afle_index[afle_index == "Pseudo_AFE"]), as.character(id)]
  ALE_junc <- cell.psi[id%in%names(afle_index[afle_index == "ALE"]), as.character(id)]


  # A3/5SS ----------------------------------
  same_start_junc2_gr <- GRanges(seqnames = same_start_junc2[, chr],
                                 ranges = IRanges(start = as.integer(same_start_junc2[, end1]) + 1, end = as.integer(same_start_junc2[, end2])),
                                 strand = same_start_junc2[, strand],
                                 names = same_start_junc2[, names])
  same_end_junc2_gr <- GRanges(seqnames = same_end_junc2[, chr],
                               ranges = IRanges(start = as.integer(same_end_junc2[, start2]), end = as.integer(same_end_junc2[, start1]) - 1),
                               strand = same_end_junc2[, strand],
                               names = same_end_junc2[, names])
  tmp1 <- as.data.table(same_start_junc2_gr)
  tmp2 <- indexedGTF$exons
  tmp3 <- as.data.table(same_end_junc2_gr)

  same_start_ass <- merge(tmp1, tmp2, by = c("seqnames", "start", "strand"))[end.x < end.y, ]
  same_end_ass <- merge(tmp3, tmp2, by = c("seqnames", "end", "strand"))[start.x > start.y, ]
  A3SS_junc <- unique(c(same_start_ass[strand == "+", names], same_end_ass[strand == "-", names]))
  A5SS_junc <- unique(c(same_start_ass[strand == "-", names], same_end_ass[strand == "+", names]))

  AS_list <- list(SE = SE_junc, A3SS = A3SS_junc, A5SS = A5SS_junc,
                  Real_AFE = R_AFE_junc, Pseudo_AFE = P_AFE_junc, ALE = ALE_junc, MXE = MXE_junc)

  # Add AS type ------------------------
  as.type <-  function(x) {
    if(sum(mapply(function(y) is.element(x, y), AS_list)) > 0){
      AS = paste0(names(AS_list)[mapply(function(y) is.element(x, y), AS_list)], collapse = "|")
    } else {
      AS = "Other"
    }
    return(AS)
  }
  cell.psi[!is.na(event), AS:=mcmapply(as.type, cell.psi[!is.na(event), id], mc.cores = NT)]
  ### assign #A5SS & AFE in + and A3SS & AFE in - to tandem AFE

  cell.psi[grep('+',id)][AS=='A5SS|Real_AFE']$AS = 'Tandem_AFE'
  cell.psi[grep('+',id,invert = T)][AS=='A3SS|Real_AFE']$AS = 'Tandem_AFE'
  return(cell.psi)
}
