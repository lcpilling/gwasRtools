#' Internal function to get -log10 p-value from test statistic
#' @param z Test statistic 
#' @param is_p Logistic. Default=FALSE. Is this a p-value? 
#' @noRd
P_neglog10 = function(z, 
                      is_p=FALSE, 
                      two_sided=TRUE) {
	if (!is.numeric(z))  stop("z needs to be numeric")
	if (is_p & two_sided)  z = z/2
	if (is_p)  z = abs(qnorm(z))
	neglog10_p = abs( ( pnorm(-abs(z), log.p=TRUE) + log(2) ) / log(10) )
	neglog10_p
}

