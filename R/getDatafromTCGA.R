
#' @title Download multi-omic data used by PhenoDriver from TCGA
#' @description
#'   Download one or more of the following types of data: gene expression data (STAR-Count)
#'   and clinical data from TCGA, MC3 mutation data from
#'   https://gdc.cancer.gov/about-data/publications/mc3-2017
#' @param cancerType select cancer type for which analysis should be run.
#' \itemize{
#' \item{ TCGA-ACC }
#' \item{ TCGA-BLCA }
#' \item{ TCGA-BRCA }
#' \item{ TCGA-CESC }
#' \item{ TCGA-CHOL }
#' \item{ TCGA-COAD }
#' \item{ TCGA-DLBC }
#' \item{ TCGA-ESCA }
#' \item{ TCGA-GBM }
#' \item{ TCGA-HNSC }
#' \item{ TCGA-KICH }
#' \item{ TCGA-KIRC }
#' \item{ TCGA-KIRP }
#' \item{ TCGA-LAML }
#' \item{ TCGA-LGG }
#' \item{ TCGA-LIHC }
#' \item{ TCGA-LUAD }
#' \item{ TCGA-LUSC }
#' \item{ TCGA-MESO }
#' \item{ TCGA-OV }
#' \item{ TCGA-PAAD }
#' \item{ TCGA-PCPG }
#' \item{ TCGA-PRAD }
#' \item{ TCGA-READ }
#' \item{ TCGA-SARC }
#' \item{ TCGA-SKCM }
#' \item{ TCGA-STAD }
#' \item{ TCGA-TGCT }
#' \item{ TCGA-THCA }
#' \item{ TCGA-THYM }
#' \item{ TCGA-UCEC }
#' \item{ TCGA-UCS }
#' \item{ TCGA-UVM }
#' }
#' @param dataType indicate data type which need to be downloaded. It should be contained
#' one or more of the following item: 'mutation', 'expression', 'clinical', 'all'.
#' @param directory Directory/Folder where the data was downloaded. Default: GDCdownload
#' @param geneType Only indicated gene type will be reserved.
#' @param cleanDownload Should downloaded data be deleted?
#' @importFrom TCGAbiolinks GDCquery
#' @importFrom TCGAbiolinks colDataPrepare
#' @importFrom TCGAbiolinks TCGAquery_SampleTypes
#' @importFrom TCGAbiolinks GDCdownload
#' @importFrom TCGAbiolinks getGDCprojects
#' @importFrom plyr alply
#' @importFrom plyr adply
#' @importFrom purrr map_df
#' @importFrom utils read.delim
#' @importFrom downloader download
#' @importFrom R.utils gunzip
#' @importFrom readr cols
#' @importFrom readr read_tsv
#' @importFrom xfun is_windows
#' @importFrom magrittr %>%
#' @importFrom dplyr bind_cols
#'
#' @return A list with the results
#' @export getDatafromTCGA
#'
#' @examples
#' \dontrun{
#' cancerData <- getDatafromTCGA(cancerType = "BRCA",
#'                               dataType = 'all',
#'                               directory = "./data/",
#'                               cleanDownload = T)
#' }
#' @author Yan Li
getDatafromTCGA <- function(cancerType,
                            dataType,
                            directory = './GDCdownload/',
                            geneType = 'protein_coding',
                            cleanDownload = FALSE){
  
  dataType <- unique(dataType)
  if(sum(dataType %in% c('mutation', 'expression', 'clinical', 'all')) != length(dataType)) {
    message('Sorry! The legal "dataType" should consists with only one or more element in
            "mutation", "expression", "clinical" and "all"')
    return(NULL)
  }
  expression <- NULL
  mutation <- NULL
  clinical <- NULL
  
  if(sum(c('expression', 'clinical', 'all') %in% dataType) >= 1) {
    CancerProject <- paste0("TCGA-",cancerType)
    DataDirectory <- paste0(directory,"GDC_",gsub("-","_",CancerProject))
    
    query <- GDCquery(project = CancerProject,
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts")
    samplesDown <- query$results[[1]]$cases
    dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown,
                                      typesample = "TP")
    dataSmNT <- TCGAquery_SampleTypes(barcode = samplesDown,
                                      typesample = "NT")
    print(paste0(CancerProject," with ", length(dataSmTP),
                 " TP, ", length(dataSmNT), " NT"))
    
    if(sum(c('expression', 'all') %in% dataType) >= 1) {
      GDCdownload(query, directory = DataDirectory)
      files <- file.path(
        query$results[[1]]$project,
        gsub(" ","_",query$results[[1]]$data_category),
        gsub(" ","_",query$results[[1]]$data_type),
        gsub(" ","_",query$results[[1]]$file_id),
        gsub(" ","_",query$results[[1]]$file_name)
      )
      files <- file.path(DataDirectory, files)
      message('----------------')
      message(paste0("Start marge patients expression data of ", CancerProject))
      expMerge <- alply(files, 1, function(f) {
        a <- read.delim(file = f, skip = 6, header = F)
        if (!('all' %in% geneType)) a <- a[a[,3] %in% geneType,]
        return(a)
      }, .progress = "time")
      expMerge <- bind_cols(expMerge[[1]][,1],
                            expMerge[[1]][,2],
                            expMerge[[1]][,3],
                            expMerge %>% map_df(4))
      cases <- ifelse(grepl("TCGA|TARGET",query$results[[1]]$project %>% unlist()),
                      query$results[[1]]$cases,query$results[[1]]$sample.submitter_id)
      colnames(expMerge) <- c('Ensembl_id', 'Gene_Symbol', 'Gene_type', cases)
      expN <- expMerge[, c(1, 2, 3, which(colnames(expMerge) %in% dataSmNT))]
      expT <- expMerge[, c(1, 2, 3, which(colnames(expMerge) %in% dataSmTP))]
      expression <- list(Normal = expN, Tumor = expT)
    }
    if(sum(c('clinical', 'all') %in% dataType) >= 1) {
      clinical <- colDataPrepare(cases)
    }
  }
  
  if(sum(c('mutation', 'all') %in% dataType) >= 1) {
    mutation <- getMC3Data(directory)
    mutation <- mutation[mutation$project_id == CancerProject, ]
  }
  
  if (cleanDownload) unlink(DataDirectory, recursive = T)
  
  return(list(expressionData = expression,
              mutationData = mutation,
              clinicalData = clinical))
}


getMC3Data <- function(directory){
  DataDirectory <- paste0(directory, '/mc3.v0.2.8.PUBLIC.maf.gz')
  fpath <- "https://api.gdc.cancer.gov/data/1c8cfe5f-e52d-41ba-94da-f15ea1337efc"
  if(is_windows()) mode <- "wb" else  mode <- "w"
  message(rep("-",100))
  options(timeout = 10000) # set 10000 second to download the file, default is 60 seconds
  message("o Starting to download Public MAF from GDC")
  message("o More information at: https://gdc.cancer.gov/about-data/publications/mc3-2017")
  message("o Please, cite: Cell Systems. Volume 6 Issue 3: p271-281.e7, 28 March 2018 10.1016/j.cels.2018.03.002")
  if(!file.exists(gsub("\\.gz", "", DataDirectory))){
    download(fpath, DataDirectory, mode = mode)
    message("o Uncompressing file")
    gunzip(DataDirectory, remove = FALSE)
  }
  message("o Reading MAF")
  maf <- read_tsv(gsub("\\.gz", "", DataDirectory),progress = TRUE, col_types = cols())
  message("o Adding project_id information")
  project <- grep("TCGA",sort(getGDCprojects()$project_id),value = TRUE)
  df <- adply(
    project,
    .margins = 1,
    .fun = function(proj) {
      samples <- TCGAbiolinks:::getSubmitterID(proj)
      return(data.frame(proj,samples))
    }
  )
  maf$project_id <- df$proj[match(substr(maf$Tumor_Sample_Barcode,1,12),df$samples)] %>% as.character
  maf$project_id[is.na(maf$project_id)] <- 'TCGA-BRCA'
  message(rep("-",100))
  return(maf)
}
