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
  # 加载包
  suppressMessages({
    library(biomaRt)
    library(BSgenome)
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(BSgenome.Mmusculus.UCSC.mm39)
    library(dplyr)
    library(Biostrings)
  })

  # 设置Ensembl数据集
  if (species == "human") {
    dataset <- "hsapiens_gene_ensembl"
    genome <- BSgenome.Hsapiens.UCSC.hg38
  } else if (species == "mouse") {
    dataset <- "mmusculus_gene_ensembl"
    genome <- BSgenome.Mmusculus.UCSC.mm39
  } else {
    stop("Unsupported species. Please use 'human' or 'mouse'.")
  }

  ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = dataset)

  # 获取数据的通用函数
  get_biomart_data <- function(attributes, filters, values, mart, max_retries) {
    data <- NULL
    for (i in seq_len(max_retries)) {
      data <- tryCatch({
        getBM(attributes = attributes, filters = filters, values = values, mart = mart)
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

  # 获取外显子和转录本信息
  exons <- get_biomart_data(
    attributes = c("chromosome_name", "exon_chrom_start", "exon_chrom_end", "strand"),
    filters = "ensembl_transcript_id", values = transcript_id, mart = ensembl, max_retries = max_retries
  )

  transcript_info <- get_biomart_data(
    attributes = c("chromosome_name", "transcript_start", "transcript_end", "strand"),
    filters = "ensembl_transcript_id", values = transcript_id, mart = ensembl, max_retries = max_retries
  )

  # 计算外显子和内含子位置
  if (transcript_info$strand == 1) {  # 正链
    exons <- exons %>%
      arrange(exon_chrom_start)  # 按起始位置升序排列

    if (nrow(exons) < 2) {
      # 只有一个外显子的情况
      introns <- NULL
      cat("只有一个外显子，位置为：\n")
      print(exons)
    } else {
      introns <- data.frame(
        chromosome_name = rep(transcript_info$chromosome_name, nrow(exons) - 1),
        start = exons$exon_chrom_end[-nrow(exons)] + 1,  # 从第一个外显子的结束位置到第二个外显子的开始位置
        end = exons$exon_chrom_start[-1] - 1  # 从第二个外显子的开始位置到第一个外显子的结束位置
      )
      # 添加strand列
      introns$strand <- exons$strand[-nrow(exons)]  # 将前一个exon的strand赋给intron
      cat("外显子位置：\n")
      print(exons)
      cat("内含子位置：\n")
      print(introns)
    }
  } else {  # 负链
    exons <- exons %>%
      arrange(desc(exon_chrom_start))  # 按起始位置降序排列

    if (nrow(exons) < 2) {
      # 只有一个外显子的情况
      introns <- NULL
      cat("只有一个外显子，位置为：\n")
      print(exons)
    } else {
      introns <- data.frame(
        chromosome_name = rep(transcript_info$chromosome_name, nrow(exons) - 1),
        start = exons$exon_chrom_end[-1] + 1,  # 直到倒数第二个外显子
        end = exons$exon_chrom_start[-nrow(exons)] - 1  # 从第二个外显子开始

      ) %>%
        arrange(desc(start))  # 按内含子起始位置降序排列
      # 添加strand列
      introns$strand <- exons$strand[-nrow(exons)]  # 将前一个exon的strand赋给intron
      cat("........外显子位置：\n")
      print(exons)
      cat("........内含子位置：\n")
      print(introns)
    }
  }

  # 提取序列的通用函数
  extract_sequences <- function(starts, ends, chromosome, strand, genome) {
    sequences <- mapply(function(start, end) {
      if (is.na(start) || is.na(end) || start > end) return(NULL)
      seq <- getSeq(genome, names = paste0("chr", chromosome), start = start, end = end)
      if (strand == -1) seq <- reverseComplement(seq)
      as.character(seq)
    }, starts, ends, SIMPLIFY = FALSE)
    return(sequences)
  }

  # 提取外显子序列
  exon_sequences <- extract_sequences(exons$exon_chrom_start, exons$exon_chrom_end, exons$chromosome_name[1], transcript_info$strand, genome)

  # 提取内含子序列（如果有内含子）
  intron_sequences <- if (!is.null(introns)) {
    extract_sequences(introns$start, introns$end, introns$chromosome_name[1], transcript_info$strand, genome)
  } else {
    NULL
  }

  # 给外显子和内含子添加序号
  exons_with_id <- paste0("Exon_", seq_along(exon_sequences))
  introns_with_id <- if (!is.null(intron_sequences)) paste0("Intron_", seq_along(intron_sequences)) else NULL

  # 将外显子和内含子合并在一起并排序
  sequences <- setNames(exon_sequences, exons_with_id)
  if (!is.null(intron_sequences)) {
    sequences <- c(sequences, setNames(intron_sequences, introns_with_id))
  }

  sorted_sequences <- sequences[order(as.integer(gsub("Exon_|Intron_", "", names(sequences))))]

  # 将序列转换为数据框
  sequences_df <- data.frame(
    ID = names(sorted_sequences),
    Sequence = unlist(sorted_sequences),
    stringsAsFactors = FALSE
  )
  # 将结果保存到CSV文件，文件名即为转录本ID
  cat("...已将结果保存到CSV文件...\n")
  write.csv(sequences_df, file = paste0(transcript_id, ".csv"), row.names = FALSE)
  return(sequences_df)
}
