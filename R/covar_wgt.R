covar_wgt <- function(x, wx, y, wy) {
	if (!(length(x) == length(wx) & length(x) == length(y) & length(y) == length(wy))) {
		stop("Vectors and their weights must be the same length.");
	}
	out <- .C("covariance_weighted", x=as.double(x), wx=as.double(wx), y=as.double(y), wy=as.double(wy), as.integer(length(x)), covar=as.double(0.0))
	return(out$covar)
}
