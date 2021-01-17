#' @title DITHER
#'
#' @description
#' \code{DITHER} Calculate ITH based on the entropy of somatic mutation and CNA profiles in tumors.
#'
#' @details
#' This function calculates the ITH score in a bulk tumor, which is the average of the entropy of its somatic mutation profiles and the entropy of its CNA profiles. 
#'
#' @param input_data_mut A dataframe or matrix which at least includes two columns ("Sample" and "VAF").
#' @param input_data_cna A dataframe or matrix which at least includes two columns ("Sample" and "Segment_Mean").
#' @export
#' @return A dataframe with 2 columns:
#' \item{Sample}{Tumor sample ID.}
#' \item{DITHER_score}{The somatic mutation and CNA entropy-based ITH score in tumors.}
#' @author Lin Li <cpu_lilin@@126.com>, Xiaosheng Wang <xiaosheng.wang@@cpu.edu.cn>
#' @examples
#' path1 = system.file("extdata", "example_mDITHER.txt", package = "DITHER", mustWork = TRUE)
#' path2 = system.file("extdata", "example_cDITHER.txt", package = "DITHER", mustWork = TRUE)
#' input_data_mut = read.table(path1, stringsAsFactors = FALSE, sep = "\t", header = TRUE, quote = "")
#' input_data_cna = read.table(path2, stringsAsFactors = FALSE, sep = "\t", header = TRUE, quote = "")
#' DITHER(input_data_mut, input_data_cna)


DITHER <- function (input_data_mut, input_data_cna) {

	# Check arguments ------------------------------------------------------------
	if (missing(input_data_mut) || !class(input_data_mut) %in% c("matrix", "data.frame") || sum(c("Sample", "VAF") %in% colnames(input_data_mut)) != 2)
		stop("'input_data_mut' is missing or incorrect")
	if (missing(input_data_mut) || !class(input_data_mut) %in% c("matrix", "data.frame") || sum(c("Sample", "Segment_Mean") %in% colnames(input_data_cna)) != 2)
		stop("'input_data_cna' is missing or incorrect")
	
	#
	sample_mut <- unique(input_data_mut[ ,"Sample"])
	sample_cna <- unique(input_data_cna[ ,"Sample"])
	sample_all <- intersect(sample_mut, sample_cna)
	
	if (length(sample_all) == 0)
		stop("DITHER score cannot be calculated!\n\t\b\b---> Please use the function mDITHER() or cDITHER()")
	
	# mDITHER ------------------------------------------------------------
	sample_mut <- unique(input_data_mut[ ,"Sample"])
	min0 <- 0
	max0 <- 1
	step0 <- (max0 - min0) / 20

	mDITHER <- c()
	for (i in 1:length(sample_mut)) {
		sample_mut0 <- input_data_mut[input_data_mut$Sample == sample_mut[i],]
		a_all <- length(sample_mut0[ ,"VAF"])
		a <- list()
		entropy <- 0
		for (j in 1:20) {
			if (j < 20) {a[j] <- sum(sample_mut0[ ,"VAF"] >= round(min0 + step0 * (j - 1), 2) & sample_mut0[ ,"VAF"] < round(min0 + step0 * j, 2), na.rm = TRUE)}
			if (j == 20) {a[j] <- sum(sample_mut0[ ,"VAF"] >= round(min0 + step0 * (j - 1), 2) & sample_mut0[ ,"VAF"] <= round(min0 + step0 * j, 2), na.rm = TRUE)}
			if (a[j] == 0) {a[j] <- a_all}
			entropy <- entropy + (-a[[j]] / a_all) * log2(a[[j]] / a_all)
		}
		res <- data.frame(Sample = sample_mut[i], Entropy = entropy)
		mDITHER <- rbind(mDITHER, res)
	}
	colnames(mDITHER) <- c("Sample", "mDITHER_score")

	# cDITHER ------------------------------------------------------------.
	sample_cna <- unique(input_data_cna[ ,"Sample"])
	min0 <- -3
	max0 <- 3
	step0 <- (max0 - min0) / 20

	cDITHER <- c()
	for (i in 1:length(sample_cna)) {
		sample_cna0 <- input_data_cna[input_data_cna[ ,"Sample"] == sample_cna[i],]	
		sample_cna0 <- sample_cna0[sample_cna0[ ,"Segment_Mean"] <= 3 & sample_cna0[ ,"Segment_Mean"] >= -3,]
		a_all <- length(sample_cna0[ ,"Segment_Mean"])
		a <- list()
		entropy <- 0
		for (j in 1:20) {
			if (j < 20) {a[j] <- sum(sample_cna0[ ,"Segment_Mean"] >= round(min0 + step0 * (j - 1), 2) & sample_cna0[ ,"Segment_Mean"] < round(min0 + step0 * j, 2), na.rm = TRUE)}
			if (j == 20) {a[j] <- sum(sample_cna0[ ,"Segment_Mean"] >= round(min0 + step0 * (j - 1), 2) & sample_cna0[ ,"Segment_Mean"] <= round(min0 + step0 * j, 2), na.rm = TRUE)}
			if (a[j] == 0) {a[j] <- a_all}
			entropy <- entropy + (-a[[j]] / a_all) * log2(a[[j]] / a_all)
		}
		res <- data.frame(Sample = sample_cna[i], Entropy = entropy)
		cDITHER <- rbind(cDITHER, res)
	}
	colnames(cDITHER) <- c("Sample", "cDITHER_score")
	
	# DITHER ------------------------------------------------------------
	DITHER <- c()
	for (i in 1:length(sample_all)) {
		entropy <- mDITHER[mDITHER[ ,"Sample"] == sample_all[i], "mDITHER_score"] + cDITHER[cDITHER[ ,"Sample"] == sample_all[i], "cDITHER_score"]
		result <- data.frame(Sample = sample_all[i], Entropy = entropy / 2)
		DITHER <- rbind(DITHER, result)
	}
	colnames(DITHER) <- c("Sample", "DITHER_score")
	return(DITHER)
}


