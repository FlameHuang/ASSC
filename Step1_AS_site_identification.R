#### AS site find ---------------------------
ASSite <- function(junc, NT = 1){

  # required package ------------------------
  usePackage("data.table")
  usePackage("parallel")
  usePackage("plyr")

  # Identify alternative splicing events ----
  juncs <- junc$junctions
  junc.names <- data.frame(chr = as.character(mapply(function(x) x[1], strsplit(juncs, split = ":"))),
                           start = as.integer(mapply(function(x) x[2], strsplit(juncs, split = "[:-]"))),
                           end = as.integer(mapply(function(x) x[3], strsplit(juncs, split = "[:-]"))),
                           strand = as.character(mapply(function(x) x[3], strsplit(juncs, split = ":"))))
  rownames(junc.names) <- juncs
  junc.names$names <- juncs
  junc.names <- junc.names[order(junc.names$chr, junc.names$start, junc.names$end),]

  junc.as <- do.call(c, mclapply(unique(junc.names$chr), function(chr) {
    junc.chr <- junc.names[junc.names$chr == chr, ]
    same.start <- plyr::dlply(junc.chr, c("start"), function(x) x)
    same.end <- plyr::dlply(junc.chr, c("end"), function(x) {
      if(nrow(x)>1) {
        return(x)
      } else {
        return(NULL)
      }
    })
    return(c(same.start, same.end[!sapply(same.end, is.null)]))
  }, mc.cores = NT))

  ## junction just has two alternative starts/ends:
  junc.as <- junc.as[sapply(junc.as, nrow) > 1]
  names(junc.as) <- sapply(junc.as, function(x) paste(x$names, collapse="@"))
  return(junc.as)
}
