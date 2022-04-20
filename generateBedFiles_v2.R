packages = c('tidyverse', 'biomaRt', 'rtracklayer', 'tools', 'readxl', 'glue')
suppressPackageStartupMessages(invisible(lapply(packages, library, character.only=T)))

# -------------------------------------------------------
## INPUTS
# Go to the directory where you would like output folders/files to be saved
setwd('')

# type organism here: mouse or human; default is mouse
organism = 'mouse'

# type selection of genes: 'all' or 'select'; default is all
selection = 'select'

# If selecting for a list of genes, paste the file path here; allowable types are .txt, .xlxs, or .csv: 
# Make sure there are no column headers on your gene list
# if your list is within the directory above, you can use the relative path 
my_genes = ''

# Select which genomic region you would like: 'genic', 'promoter', or 'both'; default is genic
region = 'both'

# If retrieving promoter coordinates, set the upstream range (in basepairs) here: default is 2000
upstream = 2000

# ---------------------------------------------------------
# DON'T ALTER THESE FUNCTIONS

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

makeOutputFolders <- function(genomic_loc) {
  folders = dir()
  if (genomic_loc == 'both') {
    if (!('BEDoutputs_promoter' %in% folders)) {
      dir.create('BEDoutputs_promoter')
    }
    if (!('BEDoutputs_genic' %in% folders)) {
      dir.create('BEDoutputs_genic')
    }
  } else {
    if (!(glue('BEDoutputs_{genomic_loc}') %in% folders)) {
      dir.create(glue('BEDoutputs_{genomic_loc}'))
    }
  }
  
}

addGeneNamesToGRangesObj <- function(organism, gr_obj) {
  tmp = as_tibble(gr_obj)
  query = entrezToEnsemblPlusGeneSymbol(organism, gr_df = tmp)
  query$entrezgene_id = as.character(query$entrezgene_id)
  
  tmp = tmp %>% inner_join(query, by = c('gene_id' = 'entrezgene_id'))
  forBed = tmp %>%
    dplyr::select(-width, -gene_id) %>%
    arrange(seqnames)
  colnames(forBed)[5] = 'gene_symbol'

  forBed
}

getPromoters <- function(txdb, upstream, organism) {
  gr = suppressWarnings(promoters(genes(txdb), upstream = upstream))
  tmp = addGeneNamesToGRangesObj(organism, gr)
  tmp
}

getCodingRegions <- function(txdb, organism) {
  gr <- GenomicRanges::reduce(cdsBy(txdb, 'gene'))
  gr.df <- as_tibble(gr) %>% 
    rename('gene_id' = 'group_name') %>%
    dplyr::select(-group) %>%
    relocate(gene_id, .after = strand)
  tmp <- addGeneNamesToGRangesObj(organism, gr.df)
  tmp
}

createBEDs <- function(txdb, 
                       genes = c('all', 'select'), 
                       genomic_loc = c('genic', 'promoter', 'both'),
                       organism = 'mouse',
                       upstream = 2000,
                       genelistPath) {
  makeOutputFolders(genomic_loc = genomic_loc)
  curr_genome = txdb$user_genome[1]

  forGenicBed <- getCodingRegions(txdb, organism)
  
  forPromoterBed <- getPromoters(txdb, upstream, organism)
  
  if (genes == 'all' && (genomic_loc == 'genic' || genomic_loc == 'both')) {
    forGenicBed %>%
      dplyr::select(-gene_symbol) %>%
      unique(.) %>%
      dplyr::rename(name = ensembl_gene_id, `#chrom` = seqnames) %>%
      write_tsv(file = glue('BEDoutputs_genic/{curr_genome}_allGenic.bed'))
  }
  
  if (genes == 'all' && (genomic_loc == 'promoter' || genomic_loc == 'both')) {
    forPromoterBed %>%
      dplyr::select(-gene_symbol) %>%
      unique(.) %>%
      dplyr::rename(name = ensembl_gene_id, `#chrom` = seqnames) %>%
      write_tsv(file=glue('BEDoutputs_promoter/{curr_genome}_allPromoter{upstream}.bed'))
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
        select_genic = forGenicBed %>%
          filter(ensembl_gene_id %in% select_genes) %>%
          dplyr::select(seqnames, start, end, strand, gene_symbol) %>%
          unique(.) %>%
          dplyr::rename(name = gene_symbol, `#chrom` = seqnames)
      } else {
        select_genic = forGenicBed %>%
          filter(gene_symbol %in% select_genes) %>%
          dplyr::select(seqnames, start, end, strand, gene_symbol) %>%
          unique(.) %>%
          dplyr::rename(name = gene_symbol, `#chrom` = seqnames)
      }
      select_genic %>%
        write_tsv(file = glue('BEDoutputs_genic/{basename(genelistPath)}_Genic.bed'))
    }
    if (genomic_loc == 'promoter' || genomic_loc == 'both') {
      if (selection_type) { 
        # this means we have ensembl ids in the user input gene list
        select_promoter = forPromoterBed %>%
          filter(ensembl_gene_id %in% select_genes) %>%
          dplyr::select(seqnames, start, end, strand, gene_symbol) %>%
          unique(.) %>%
          dplyr::rename(name = gene_symbol, `#chrom` = seqnames)
      } else {
        select_promoter = forPromoterBed %>%
          filter(gene_symbol %in% select_genes) %>%
          dplyr::select(seqnames, start, end, strand, gene_symbol) %>%
          unique(.) %>%
          dplyr::rename(name = gene_symbol, `#chrom` = seqnames)
      }
      select_promoter %>%
        write_tsv(file = glue('BEDoutputs_promoter/{basename(genelistPath)}_Promoter{upstream}.bed'))
    }
  }
}

# ---------------------------------------------------------

txdb <- loadTxDbFromOrganism(organism = organism)

mart <- generateMartObject(organism)

createBEDs(txdb = txdb, genes = selection, genomic_loc = region, 
           organism = organism, upstream = upstream, genelistPath = my_genes)





