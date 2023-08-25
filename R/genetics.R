#' Estimate Lambda GC from p-values
#'
#' @description Estimate inflation of test statistics. Lambda GC compares the median test statistic against the expected median test statistic under the null hypothesis of no association. For well-powered quantitative traits with a known polygenic inheritance, we expect inflation of lambda GC, but for traits with no expected association, we expect lambda GC to be around 1.
#'
#' Original function by Jing Hua Zhao (https://jinghuazhao.github.io/)
#'
#' @return Returns a vector
#'
#' @author Jing Hua Zhao
#'
#' @name lambda_gc
#'
#' @param p a numeric vector of p-values
#'
#' @examples
#' lambda_gc(p)
#'
#' @export
#'

lambda_gc = function(p) {
	if (!is.numeric(p))  stop("p needs to be numeric")
	
	p = p[!is.na(p)]
	n = length(p)
	
	if (any(p > 1) | any(p < 0))  stop("These don't look like p-values")
	
	observed = qchisq(p,1,lower.tail=FALSE)
	expected = qchisq(1:n/n,1,lower.tail=FALSE)
	
	lambda = median(observed)/median(expected)
	return(lambda)
}


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

