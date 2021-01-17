#' @title mDITHER
#'
#' @description
#' \code{mDITHER} Calculate ITH based on the entropy of somatic mutation profiles in tumors.
#'
#' @details
#' This function calculates the ITH score in a bulk tumor based on the entropy of its somatic mutation profiles.
#'
#' @param input_data_mut A dataframe or matrix which at least includes two columns ("Sample" and "VAF").
#' @export
#' @return A dataframe with 2 columns:
#' \item{Sample}{Tumor sample ID.}
#' \item{mDITHER_score}{The somatic mutation entropy-based ITH score in tumors.}
#' @author Lin Li <cpu_lilin@@126.com>, Xiaosheng Wang <xiaosheng.wang@@cpu.edu.cn>
#' @examples
#' path = system.file("extdata", "example_mDITHER.txt", package = "DITHER", mustWork = TRUE)
#' input_data_mut = read.table(path, stringsAsFactors = FALSE, sep = "\t", header = TRUE, quote = "")
#' mDITHER(input_data_mut)


mDITHER <- function (input_data_mut) {

	# Check arguments ------------------------------------------------------------
	if (missing(input_data_mut) || !class(input_data_mut) %in% c("matrix", "data.frame") || sum(c("Sample", "VAF") %in% colnames(input_data_mut)) != 2)
		stop("'input_data_mut' is missing or incorrect")
	
	
	#
	sample_mut <- unique(input_data_mut[ ,"Sample"])
	min0 <- 0
	max0 <- 1
	step0 <- (max0 - min0) / 20

	# mDITHER ------------------------------------------------------------
	mDITHER <- c()
	for (i in 1:length(sample_mut)) {
		sample_mut0 <- input_data_mut[input_data_mut[ ,"Sample"] == sample_mut[i],]
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
	return(mDITHER)
}

