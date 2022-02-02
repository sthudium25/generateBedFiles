packages = c('tidyverse', 'biomaRt', 'rtracklayer', 'tools', 'readxl', 'biomaRt')
suppressPackageStartupMessages(invisible(lapply(packages, library, character.only=T)))

# -------------------------------------------------------
## INPUTS
# type organism here: mouse or human; default is mouse
organism = 'human'

# type selection of genes: 'all' or 'select'; default is all
selection = 'all'

# If selecting for a list of genes, paste the file path here; allowable types are .txt, .xlxs, or .csv: 
# Make sure there are no column headers on your gene list
my_genes = 'path/to/my/genes'

# Select which genomic region you would like: 'genic', 'promoter', or 'both'; default is genic
region = 'genic'

# If retrieving promoter coordinates, set the upstream range (in basepairs) here: default is 2000
upstream = 2000

# ---------------------------------------------------------

# This function trims the 5' end of genes if there is a cds within a specified bp distance 
# DON'T ALTER WITH THIS FUNCTION
trimFivePrimeSide <- function(gr, remove)
{
  stopifnot(is(gr, "GenomicRanges"), is(remove, "GenomicRanges"))
  
  mapping <- as(findOverlaps(gr, remove), "IntegerList")
  
  new_start_on_plus_strand <- pmin(max(extractList(end(remove),
                                                   mapping)),
                                   end(gr)) + 1L
  new_end_on_minus_strand <- pmax(min(extractList(start(remove),
                                                  mapping)),
                                  start(gr)) - 1L
  
  on_plus_strand <- as.logical(strand(gr) == "+")
  on_minus_strand <- as.logical(strand(gr) == "-")
  has_hit <- lengths(mapping) != 0L
  
  new_start <- ifelse(on_plus_strand & has_hit,
                      new_start_on_plus_strand,
                      start(gr))
  new_end <- ifelse(on_minus_strand & has_hit,
                    new_end_on_minus_strand,
                    end(gr))
  ranges(gr) <- IRanges(new_start, new_end, names=names(gr))
  
  gr
}

loadTxDbFromOrganism <- function(organism) {
  if (organism == 'mouse') {
    library(TxDb.Mmusculus.UCSC.mm10.knownGene)
    txdb = TxDb.Mmusculus.UCSC.mm10.knownGene
  }
  if (organism == 'human') {
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
  }
  txdb
}

generateMartObject <- function(organism) {
  if (organism == 'mouse') {
    ensembl = useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
  } else if (organism == 'human') {
    ensembl = useMart(biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl')
  }
  ensembl
}

readInGeneList <- function(path) {
  ext = tools::file_ext(path)
  if (! (ext %in% c('txt', 'xlsx', 'csv'))) {
    print('Invalid genelist file type. Must be one of: .txt, .xlsx, or .csv')
    stop()
  } else if (ext == 'txt') {
    out = read_tsv(path, col_names = F)
  } else if (ext == 'xlsx') {
    out = read_xlsx(path, col_names = F)
  } else if (ext == 'csv') {
    out = read_csv(path, col_names = F)
  }
  out
}

txdb = loadTxDbFromOrganism(organism = organism)

mart = generateMartObject(organism)


createBEDs <- function(genes = c('all', 'select'), 
                       genomic_loc) {
  if (genes == 'all') {
    
  }
}

select_genes = readInGeneList(my_genes)

gr = genes(txdb)
tmp = as_tibble(gr)
query = getBM(attributes = c("entrezgene_id", 'hgnc_symbol'),
              filters = )
