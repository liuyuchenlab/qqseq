####qqseq
#自定义函数#
#' Title
#'
#' @param transcript_id esembl transcript_id
#' @param species species
#' @param max_retries he maximum number of retries allowed. Default is 20.
#'
#' @return sequences of transcript
#' @export
#'
#' @examples
#' result <- qqseq("ENSMUST00000127786", "mouse")
qqseq <- function(transcript_id, species, max_retries = 20) {
  # Load packages
  suppressMessages({
    if (!requireNamespace("biomaRt", quietly = TRUE)) stop("Package 'biomaRt' needed but not installed.")
    if (!requireNamespace("BSgenome", quietly = TRUE)) stop("Package 'BSgenome' needed but not installed.")
    if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) stop("Package 'BSgenome.Hsapiens.UCSC.hg38' needed but not installed.")
    if (!requireNamespace("BSgenome.Mmusculus.UCSC.mm39", quietly = TRUE)) stop("Package 'BSgenome.Mmusculus.UCSC.mm39' needed but not installed.")
    if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' needed but not installed.")
    if (!requireNamespace("Biostrings", quietly = TRUE)) stop("Package 'Biostrings' needed but not installed.")
  })

  # Set Ensembl dataset
  if (species == "human") {
    dataset <- "hsapiens_gene_ensembl"
    genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  } else if (species == "mouse") {
    dataset <- "mmusculus_gene_ensembl"
    genome <- BSgenome.Mmusculus.UCSC.mm39::BSgenome.Mmusculus.UCSC.mm39
  } else {
    stop("Unsupported species. Please use 'human' or 'mouse'.")
  }

  ensembl <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = dataset)

  # Function to get data
  get_biomart_data <- function(attributes, filters, values, mart, max_retries) {
    data <- NULL
    for (i in seq_len(max_retries)) {
      data <- tryCatch({
        biomaRt::getBM(attributes = attributes, filters = filters, values = values, mart = mart)
      }, error = function(e) {
        cat("Error in getBM (attempt ", i, "): ", e$message, "\n")
        NULL
      })
      if (!is.null(data) && nrow(data) > 0) break
    }
    if (is.null(data) || nrow(data) == 0) {
      stop("Data not found or error in retrieving data.")
    }
    return(data)
  }

  # Get exon and transcript info
  exons <- get_biomart_data(
    attributes = c("chromosome_name", "exon_chrom_start", "exon_chrom_end", "strand"),
    filters = "ensembl_transcript_id", values = transcript_id, mart = ensembl, max_retries = max_retries
  )

  transcript_info <- get_biomart_data(
    attributes = c("chromosome_name", "transcript_start", "transcript_end", "strand", "external_gene_name"),
    filters = "ensembl_transcript_id", values = transcript_id, mart = ensembl, max_retries = max_retries
  )

  # Extract gene name
  gene_name <- unique(transcript_info$external_gene_name)

  # Calculate exon and intron positions
  if (transcript_info$strand == 1) {  # Positive strand
    exons <- dplyr::arrange(exons, exon_chrom_start)  # Sort by start position in ascending order

    if (nrow(exons) < 2) {
      # Only one exon
      introns <- NULL
      cat("Only one exon, position:\n")
      print(exons)
    } else {
      introns <- data.frame(
        chromosome_name = rep(transcript_info$chromosome_name, nrow(exons) - 1),
        start = exons$exon_chrom_end[-nrow(exons)] + 1,  # From the end position of the first exon to the start position of the second exon
        end = exons$exon_chrom_start[-1] - 1  # From the start position of the second exon to the end position of the first exon
      )
      # Add strand column
      introns$strand <- exons$strand[-nrow(exons)]  # Assign the strand of the previous exon to the intron
      cat("Exon positions:\n")
      print(exons)
      cat("Intron positions:\n")
      print(introns)
    }
  } else {  # Negative strand
    exons <- dplyr::arrange(exons, desc(exon_chrom_start))  # Sort by start position in descending order

    if (nrow(exons) < 2) {
      # Only one exon
      introns <- NULL
      cat("Only one exon, position:\n")
      print(exons)
    } else {
      introns <- data.frame(
        chromosome_name = rep(transcript_info$chromosome_name, nrow(exons) - 1),
        start = exons$exon_chrom_end[-1] + 1,  # Until the second last exon
        end = exons$exon_chrom_start[-nrow(exons)] - 1  # From the second exon
      ) %>%
        dplyr::arrange(desc(start))  # Sort introns by start position in descending order
      # Add strand column
      introns$strand <- exons$strand[-nrow(exons)]  # Assign the strand of the previous exon to the intron
      cat("Exon positions:\n")
      print(exons)
      cat("Intron positions:\n")
      print(introns)
    }
  }

  # Function to extract sequences
  extract_sequences <- function(starts, ends, chromosome, strand, genome) {
    sequences <- mapply(function(start, end) {
      if (is.na(start) || is.na(end) || start > end) return(NULL)
      seq <- BSgenome::getSeq(genome, names = paste0("chr", chromosome), start = start, end = end)
      if (strand == -1) seq <- Biostrings::reverseComplement(seq)
      as.character(seq)
    }, starts, ends, SIMPLIFY = FALSE)
    return(sequences)
  }

  # Extract exon sequences
  exon_sequences <- extract_sequences(exons$exon_chrom_start, exons$exon_chrom_end, exons$chromosome_name[1], transcript_info$strand, genome)

  # Extract intron sequences (if any)
  intron_sequences <- if (!is.null(introns)) {
    extract_sequences(introns$start, introns$end, introns$chromosome_name[1], transcript_info$strand, genome)
  } else {
    NULL
  }

  # Assign IDs to exons and introns
  exons_with_id <- paste0("Exon_", seq_along(exon_sequences))
  introns_with_id <- if (!is.null(intron_sequences)) paste0("Intron_", seq_along(intron_sequences)) else NULL

  # Combine exons and introns and sort
  sequences <- setNames(exon_sequences, exons_with_id)
  if (!is.null(intron_sequences)) {
    sequences <- c(sequences, setNames(intron_sequences, introns_with_id))
  }

  sorted_sequences <- sequences[order(as.integer(gsub("Exon_|Intron_", "", names(sequences))))]

  # Convert sequences to data frame
  sequences_df <- data.frame(
    ID = names(sorted_sequences),
    Sequence = unlist(sorted_sequences),
    stringsAsFactors = FALSE
  )

  # Create file name with date, gene name, and transcript ID
  date_str <- format(Sys.Date(), "%Y-%m-%d")
  file_name <- paste0(transcript_id, "_", gene_name, "_", species, ".csv")

  # Save results to CSV file with the new name
  cat("Results saved to CSV file.\n")
  utils::write.csv(sequences_df, file = file_name, row.names = FALSE)

  return(sequences_df)
}
