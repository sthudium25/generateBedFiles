packages = c('tidyverse', 'biomaRt', 'rtracklayer', 'tools', 'readxl', 'glue')
suppressPackageStartupMessages(invisible(lapply(packages, library, character.only=T)))

# -------------------------------------------------------
## INPUTS
# type organism here: mouse or human; default is mouse
organism = 'human'

# type selection of genes: 'all' or 'select'; default is all
selection = 'all'

# If selecting for a list of genes, paste the file path here; allowable types are .txt, .xlxs, or .csv: 
# Make sure there are no column headers on your gene list
my_genes = 'Orf8_DEG_lists/Orf8overGFP_up'

# Select which genomic region you would like: 'genic', 'promoter', or 'both'; default is genic
region = 'both'

# If retrieving promoter coordinates, set the upstream range (in basepairs) here: default is 2000
upstream = 2000

# ---------------------------------------------------------
# DON'T ALTER THESE FUNCTIONS

# This function trims the 5' end of genes if there is a cds within a specified bp distance 
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

# This function loads the Transcript Database corresponding to your selected organism
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

# This function links to ensembl so that we can convert between different gene identifiers
generateMartObject <- function(organism) {
  if (organism == 'mouse') {
    ensembl = useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
  } else if (organism == 'human') {
    ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  }
  ensembl
}

# This funciton reads in your provided list of genes
readInGeneList <- function(path) {
  ext = tools::file_ext(path)
  if (! (ext %in% c('txt', 'xlsx', 'csv', ''))) {
    print('Invalid genelist file type. Must be one of: .txt, .xlsx, or .csv')
    stop()
  } else if (ext == 'txt' || ext == '') {
    out = read_tsv(path, col_names = F)
  } else if (ext == 'xlsx') {
    out = read_xlsx(path, col_names = F)
  } else if (ext == 'csv') {
    out = read_csv(path, col_names = F)
  }
  out
}

# This is a helper function to retrieve ensembl ids and gene names given a list of entrez ids
entrezToEnsemblPlusGeneSymbol <- function(organism, gr_df) {
  if (organism == 'human') {
    query = getBM(attributes = c("entrezgene_id", 'hgnc_symbol', 'ensembl_gene_id'),
                  filters = 'entrezgene_id',
                  values = gr_df$gene_id,
                  mart = mart)
  } else if (organism == 'mouse') {
    query = getBM(attributes = c("entrezgene_id", 'mgi_symbol', 'ensembl_gene_id'),
                  filters = 'entrezgene_id',
                  values = gr_df$gene_id,
                  mart = mart)
  }
  query
}

trimGenesWrapper <- function(txdb, genomic_loc, bed_df, promoter_width) {
  if (genomic_loc == 'promoter' || genomic_loc == 'both') {
    forBed_gr = makeGRangesFromDataFrame(bed_df, keep.extra.columns = T)
    geneFlanks = flank(forBed_gr, 
                       width = upstream)
    geneFlanksDF = as.data.frame(geneFlanks) 
    cdsCoords = cds(txdb)
    gnflankTrimmed = as_tibble(trimFivePrimeSide(geneFlanks, cdsCoords))
    return(gnflankTrimmed)
  } else {
    print('Promoter not selected')
  }
}

txdb = loadTxDbFromOrganism(organism = organism)

mart = generateMartObject(organism)

createBEDs <- function(txdb, 
                       genes = c('all', 'select'), 
                       genomic_loc = c('genic', 'promoter', 'both'),
                       organism = 'mouse',
                       upstream = 2000,
                       genelistPath) {
  curr_genome = txdb$user_genome[1]
  gr = suppressMessages(genes(txdb))
  tmp = as_tibble(gr)
  query = entrezToEnsemblPlusGeneSymbol(organism, gr_df = tmp)
  query$entrezgene_id = as.character(query$entrezgene_id)
  
  tmp = tmp %>% inner_join(query, by = c('gene_id' = 'entrezgene_id'))
  forBed = tmp %>%
    dplyr::select(-width, -gene_id) %>%
    arrange(seqnames)
  colnames(forBed)[5] = 'gene_symbol'
  
  gnflankTrimmed = trimGenesWrapper(txdb, genomic_loc = genomic_loc, 
                                    bed_df = forBed, promoter_width = upstream)
 
  if (genes == 'all' && (genomic_loc == 'genic' || genomic_loc == 'both')) {
    forBed %>%
      dplyr::select(-gene_symbol) %>%
      dplyr::rename(name = ensembl_gene_id, chrom = seqnames) %>%
      write_tsv(file = glue('outputs/{curr_genome}_allGenic.bed'))
  }
  
  if (genes == 'all' && (genomic_loc == 'promoter' || genomic_loc == 'both')) {
    gnflankTrimmed %>%
      dplyr::select(-width, -gene_symbol) %>%
      dplyr::rename(name = ensembl_gene_id, chrom = seqnames) %>%
      write_tsv(file=glue('outputs/{curr_genome}_allPromoter{upstream}.bed'))
  }
  
  if (genes == 'select') {
    select_genes = readInGeneList(genelistPath)
    select_genes = select_genes[,1][[1]]
    # Figure out which type of IDs are being used in the selected genes list
    selection_type = (grepl('ENSG', select_genes[1], fixed=T) || 
                        grepl('ENSMUS', select_genes[1], fixed=T))
    
    if (genomic_loc == 'genic' || genomic_loc == 'both') {
      if (selection_type) { 
        # this means we have ensembl ids in the user input gene list
        select_genic = forBed %>%
          filter(ensembl_gene_id %in% select_genes) %>%
          dplyr::select(seqnames, start, end, strand, gene_symbol) %>%
          dplyr::rename(name = gene_symbol, chrom = seqnames)
      } else {
        select_genic = forBed %>%
          filter(gene_symbol %in% select_genes) %>%
          dplyr::select(seqnames, start, end, strand, gene_symbol) %>%
          dplyr::rename(name = gene_symbol, chrom = seqnames)
      }
      select_genic %>%
        write_tsv(file = glue('outputs/{curr_genome}_{basename(genelistPath)}_Genic.bed'))
    }
    if (genomic_loc == 'promoter' || genomic_loc == 'both') {
      if (selection_type) { 
        # this means we have ensembl ids in the user input gene list
        select_promoter = gnflankTrimmed %>%
          filter(ensembl_gene_id %in% select_genes) %>%
          dplyr::select(seqnames, start, end, strand, gene_symbol) %>%
          dplyr::rename(name = gene_symbol, chrom = seqnames)
      } else {
        select_promoter = gnflankTrimmed %>%
          filter(gene_symbol %in% select_genes) %>%
          dplyr::select(seqnames, start, end, strand, gene_symbol) %>%
          dplyr::rename(name = gene_symbol, chrom = seqnames)
      }
      select_promoter %>%
        write_tsv(file = glue('outputs/{curr_genome}_{basename(genelistPath)}_Promoter{upstream}.bed'))
    }
  }
}

createBEDs(txdb = txdb, genes = selection, genomic_loc = region, 
           organism = organism, upstream = upstream, genelistPath = my_genes)






