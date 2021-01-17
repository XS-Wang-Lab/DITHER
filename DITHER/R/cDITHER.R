#' @title cDITHER
#'
#' @description
#' \code{cDITHER} Calculate ITH based on the entropy of CNA profiles in tumors.
#'
#' @details
#' This function calculates the ITH score in a bulk tumor based on the entropy of its CNA profiles.
#'
#' @param input_data_cna A dataframe or matrix which at least includes two columns ("Sample" and "Segment_Mean").
#' @export
#' @return A dataframe with 2 columns:
#' \item{Sample}{Tumor sample ID.}
#' \item{cDITHER_score}{The CNA entropy-based ITH score in tumors.}
#' @author Lin Li <cpu_lilin@@126.com>, Xiaosheng Wang <xiaosheng.wang@@cpu.edu.cn>
#' @examples
#' path = system.file("extdata", "example_cDITHER.txt", package = "DITHER", mustWork = TRUE)
#' input_data_cna = read.table(path, stringsAsFactors = FALSE, sep = "\t", header = TRUE, quote = "")
#' cDITHER(input_data_cna)


cDITHER <- function (input_data_cna) {

	# Check arguments ------------------------------------------------------------
	if (missing(input_data_cna) || !class(input_data_cna) %in% c("matrix", "data.frame") || sum(c("Sample", "Segment_Mean") %in% colnames(input_data_cna)) != 2)
		stop("'input_data_cna' is missing or incorrect")
	
	
	#
	sample_cna <- unique(input_data_cna[ ,"Sample"])
	min0 <- -3
	max0 <- 3
	step0 <- (max0 - min0) / 20

	# cDITHER ------------------------------------------------------------
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
	return(cDITHER)
}

