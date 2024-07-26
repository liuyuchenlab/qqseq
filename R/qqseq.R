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
  # 确认包已安装，否则停止执行并提示安装包
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop("Package 'biomaRt' is required but not installed.")
  }
  if (!requireNamespace("BSgenome", quietly = TRUE)) {
    stop("Package 'BSgenome' is required but not installed.")
  }
  if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
    stop("Package 'BSgenome.Hsapiens.UCSC.hg38' is required but not installed.")
  }
  if (!requireNamespace("BSgenome.Mmusculus.UCSC.mm39", quietly = TRUE)) {
    stop("Package 'BSgenome.Mmusculus.UCSC.mm39' is required but not installed.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required but not installed.")
  }
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Package 'Biostrings' is required but not installed.")
  }
  # 加载包
  suppressMessages(library(biomaRt))
  suppressMessages(library(BSgenome))
  suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
  suppressMessages(library(BSgenome.Mmusculus.UCSC.mm39))
  suppressMessages(library(dplyr))
  suppressMessages(library(Biostrings))
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

  # 获取转录本外显子信息
  exons <- NULL
  for (i in seq_len(max_retries)) {
    exons <- tryCatch({
      getBM(attributes = c("chromosome_name", "exon_chrom_start", "exon_chrom_end", "strand"),
            filters = "ensembl_transcript_id", values = transcript_id, mart = ensembl)
    }, error = function(e) {
      cat("Error in getBM for exons (attempt ", i, "): ", e$message, "\n")
      NULL
    })
    if (!is.null(exons) && nrow(exons) > 0) break
  }

  if (is.null(exons) || nrow(exons) == 0) {
    stop("Transcript not found or error in retrieving exon data.")
  }

  # 获取转录本信息
  transcript_info <- NULL
  for (i in seq_len(max_retries)) {
    transcript_info <- tryCatch({
      getBM(attributes = c("chromosome_name", "transcript_start", "transcript_end", "strand"),
            filters = "ensembl_transcript_id", values = transcript_id, mart = ensembl)
    }, error = function(e) {
      cat("Error in getBM for transcript_info (attempt ", i, "): ", e$message, "\n")
      NULL
    })
    if (!is.null(transcript_info) && nrow(transcript_info) > 0) break
  }

  if (is.null(transcript_info) || nrow(transcript_info) == 0) {
    stop("Transcript information not found or error in retrieving transcript data.")
  }

  # 计算外显子和内含子位置
  if (transcript_info$strand == 1) {  # 正链
    exons <- exons %>%
      arrange(exon_chrom_start)  # 按起始位置升序排列
    introns <- data.frame(
      chromosome_name = rep(transcript_info$chromosome_name, nrow(exons) - 1),
      start = exons$exon_chrom_end[-nrow(exons)] + 1,
      end = exons$exon_chrom_start[-1] - 1
    )
  } else {  # 负链
    exons <- exons %>%
      arrange(desc(exon_chrom_start))  # 按起始位置降序排列
    introns <- data.frame(
      chromosome_name = rep(transcript_info$chromosome_name, nrow(exons) - 1),
      start = exons$exon_chrom_start[-1] + 1,
      end = exons$exon_chrom_end[-nrow(exons)] - 1
    ) %>%
      arrange(desc(start))  # 按内含子起始位置降序排列
  }

  # 移除无效的内含子区域
  introns <- introns %>%
    filter(start <= end)

  # 提取外显子序列
  exon_sequences <- mapply(function(start, end) {
    seq <- getSeq(genome, names = paste0("chr", exons$chromosome_name[1]), start = start, end = end)
    if (transcript_info$strand == -1) {
      seq <- reverseComplement(seq)
    }
    as.character(seq)  # 转换为字符
  }, exons$exon_chrom_start, exons$exon_chrom_end, SIMPLIFY = FALSE)

  # 提取内含子序列
  intron_sequences <- mapply(function(start, end) {
    if (start > end) return(NULL)
    seq <- getSeq(genome, names = paste0("chr", introns$chromosome_name[1]), start = start, end = end)
    if (transcript_info$strand == -1) {
      seq <- reverseComplement(seq)
    }
    as.character(seq)  # 转换为字符
  }, introns$start, introns$end, SIMPLIFY = FALSE)

  # 给外显子和内含子添加序号
  exons_with_id <- paste0("Exon_", seq_along(exon_sequences))
  introns_with_id <- paste0("Intron_", seq_along(intron_sequences))

  # 将外显子和内含子合并在一起
  sequences <- c(
    setNames(exon_sequences, exons_with_id),
    setNames(intron_sequences, introns_with_id)
  )

  # 根据序列位置进行排序（将外显子按位置排序，将内含子按位置排序）
  sorted_sequences <- sequences[order(as.integer(gsub("Exon_|Intron_", "", names(sequences))))]

  # 将序列转换为数据框
  sequences_df <- data.frame(
    ID = names(sorted_sequences),
    Sequence = unlist(sorted_sequences),
    stringsAsFactors = FALSE
  )

  # 返回数据框
  return(sequences_df)
}
