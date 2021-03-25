#### AS site read ---------------------------

# personal function
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dependencies = TRUE)
  require(p, character.only = TRUE)
}

# function for replace NA with 0 in data.table
DT.replace.NA  <- function(x) {
  for (j in seq_len(ncol(x))) set(x,which(is.na(x[[j]])),j,0)
}



readSJ <- function(file, input.rowsum.cutoff = NULL){
  if(!file.exists(file)) {
    stop("Input file does not exist!")
  }

  # required package ------------------------
  usePackage("data.table")

  # read files ------------------------------
  input <- fread(file)

  # replace missing value with 0 ------------
  DT.replace.NA(input)

  # filter by the rowsum
  input.rs <- rowSums(input[, -1])
  input.rowsum.cutoff <- ifelse(is.null(input.rowsum.cutoff),
                                max(unique(quantile(input.rs, probs = seq(0, 1, .05)))[2], 2*(ncol(input)-1)),
                                input.rowsum.cutoff)
  junc <- input[input.rs > input.rowsum.cutoff, ]
  setkey(junc, junctions)
  print(paste("Before filtering ", nrow(input), " After filtering ", nrow(junc), " at ", print(Sys.time()), collapse = ""))
  return(junc)
}
