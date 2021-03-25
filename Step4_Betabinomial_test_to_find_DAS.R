### modify the funtion to detect differential splicing events -------
BetabinomialDSS <- function(as.site, junc.tab, sampleInfo, cell1, cell2, NT = 1){
  if(!is.list(as.site)) {
    stop("Input of as.site must be the output list of ASSite!")
  }

  # required package ------------------------
  usePackage("data.table")
  usePackage("parallel")
  usePackage("VGAM")

  print(paste("Starting ", cell1, " vs ", cell2, " at ", Sys.time(), sep = ""))
  cell.psi <- parallel::mclapply(as.site, function(sjs){

    sjs.names <- sjs$names[order(as.numeric(sjs$end) - as.numeric(sjs$start))]
    #col.select
    tab0 <- as.matrix(t(junc.tab[eval(sjs.names), sampleInfo$ID[sampleInfo$Group %in% c(cell1, cell2)], with = F]))
    colnames(tab0) <- sjs.names

    if (sum(tab0) > 100 & sum((colSums(tab0)/sum(tab0)) > 0.05) == 2) {
      tab <- tab0[, (colSums(tab0)/sum(tab0)) > 0.05]
      event <- sum((colSums(tab)/sum(tab)) > 0.05)
      sjs.use <- sjs[sjs$names %in% colnames(tab), ]
      dist <- max(c(diff(as.numeric(sjs.use$start)), diff(as.numeric(sjs.use$end))))
      fit1 <- tryCatch(VGAM::vglm(tab ~ as.factor(sampleInfo[sampleInfo$Group %in% c(cell1, cell2), "Group"]), betabinomial, trace = FALSE, maxit = 60, subset = rowSums(tab) > 0), error = function(e) NULL)
      fit0 <- tryCatch(VGAM::vglm(tab ~ 1, betabinomial, trace = FALSE, maxit = 60, subset = rowSums(tab) > 0), error = function(e) NULL)
      beta.model <- fitted(fit1)
      cell1_PSI <- tryCatch(paste(round(beta.model[1], 3), 1 - round(beta.model[1], 3), sep= ";"), error = function(e) NA)
      cell2_PSI <- tryCatch(paste(round(beta.model[nrow(beta.model)], 3), 1 - round(beta.model[nrow(beta.model)], 3), sep = ";"), error=function(e) NA)
    } else if (sum(colSums(tab0) > 10 & colSums(tab0)/sum(tab0) > 0.05) > 2) {
      tab <- tab0[, (colSums(tab0)/sum(tab0)) > 0.05]
      event <- sum((colSums(tab)/sum(tab)) > 0.05)
      sjs.use <- sjs[sjs$names %in% colnames(tab), ]
      dist.all <- c(diff(as.numeric(sjs.use$start)), diff(as.numeric(sjs.use$end)))
      dist <- paste0(abs(dist.all[abs(dist.all) > 0]), collapse=";")
      fit1 <- tryCatch(VGAM::vglm(tab ~ as.factor(sampleInfo[sampleInfo$Group %in% c(cell1, cell2), "Group"]), dirmultinomial, trace = FALSE, maxit = 60, subset = rowSums(tab) > 0), error=function(e) NULL)
      fit0 <- tryCatch(VGAM::vglm(tab ~ 1, dirmultinomial, trace = FALSE, maxit = 60, subset = rowSums(tab) > 0), error=function(e) NULL)
      dir.model <- fitted(fit1)
      cell1_PSI <- tryCatch(paste0(round(dir.model[1, ], 3), collapse = ";"), error = function(e) NA)
      cell2_PSI <- tryCatch(paste0(round(dir.model[nrow(dir.model), ], 3), collapse = ";"), error = function(e) NA)
    } else {
      fit1 <- NULL
      fit0 <- NULL
    }
    if(!is.null(fit1) & !is.null(fit0)){
      pvalue <- tryCatch(1 - pchisq(2*(fit1@criterion$loglikelihood - fit0@criterion$loglikelihood), abs(fit0@df.residual - fit1@df.residual)), error = function(e) NA)
      return(
        data.frame(id = paste(colnames(tab), collapse = "@"),
                   event = event,
                   pvalue = pvalue,
                   cell1_PSI = cell1_PSI,
                   cell2_PSI = cell2_PSI,
                   cell1_sum = paste0(colSums(tab[sampleInfo[sampleInfo$Group %in% c(cell1, cell2), "Group"] == cell1, ]), collapse = ";"),
                   cell2_sum = paste0(colSums(tab[sampleInfo[sampleInfo$Group %in% c(cell1, cell2), "Group"] == cell2, ]), collapse = ";"),
                   dist = dist,
                   cell = paste0(cell1, "&", cell2)))
    } else {
      return(data.frame(
        id = paste(colnames(tab0), collapse="@"),
        event = NA,
        pvalue = NA,
        cell1_PSI = NA,
        cell2_PSI = NA,
        cell1_sum = NA,
        cell2_sum = NA,
        dist = NA,
        cell = paste0(cell1, "&", cell2)))
    }
  }, mc.cores = NT)

  cell.psi = do.call(rbind, cell.psi)
  # print(cell.psi)
  cell_psi <- as.data.table(cell.psi, keep.rownames = "SJ.cluster")
  cell_psi[, id:=as.character(id)]
  setkey(cell_psi, id)
  ### we extract the IC-PSI of the shortest junction
  cell_psi[,cell1_inclusion:=as.numeric(mapply(function(x) x[1], strsplit(as.character(cell1_PSI),';')))]
  cell_psi[,cell2_inclusion:=as.numeric(mapply(function(x) x[1], strsplit(as.character(cell2_PSI),';')))]
  ###  the shortest junction represent the inclusion isoform
  cell_psi[,inclusion_junc:=mapply(function(x) x[1], strsplit(as.character(id),'@'))]
  cell_psi[,FDR:=p.adjust(pvalue,method = 'BH')]
  cell_psi[,delta:=cell1_inclusion-cell2_inclusion]
  print(paste("Finished ", cell1, " vs ", cell2, " at ", Sys.time(), sep = ""))
  return(cell_psi)
}
