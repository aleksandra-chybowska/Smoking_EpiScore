library(BSgenome.Hsapiens.UCSC.hg38) 
require(Biostrings)
require(parallel)

Find_CpG <- function(Genome, Cores){
  if (class(Genome) != "BSgenome") stop("Genome must be a BSgenome!")
  
  CpG <- mclapply(seqlevels(Genome), function(x) start(matchPattern("CG", Genome[[x]])), mc.cores = Cores)
  return(
    suppressWarnings(
      do.call(c, mclapply(1:length(seqlevels(Genome)), function(x) GRanges(names(Genome)[x], 
                                                                           IRanges(CpG[[x]], width = 2)
      ), mc.cores=Cores))
    )
  )
}

## Example:
hg38.CpG <- Find_CpG(Genome = BSgenome.Hsapiens.UCSC.hg38, Cores = 4)