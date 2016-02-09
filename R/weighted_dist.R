weighted_dist <- function(m,w,p = 2) {
	# p = 2 = euclidean, else = minkowski, p = 1 = manhattan (abs)
	if (p < 0) {stop("p must be positive")}
	if (prod(dim(m) == dim(w))) {
		nrow = length(m[,1])
		ncol = length(m[1,])
		dist = matrix(rep(0, times = nrow*nrow), nrow = nrow, ncol=nrow)
		out <- .C("distance_wgt", y=as.double(m), w=as.double(w), nrow=as.integer(length(m[,1])), ncol=as.integer(length(m[1,])), exponent = as.double(p), out = as.double(dist) )
		
		# Format output
		output = matrix(out$dist, nrow = nrow, ncol = nrow)
		colnames(output) <- rownames(m)
		rownames(output) <- rownames(m) 
		return(as.dist(output))
	} else {
		stop("Weights matrix must be same dimension as data matrix.");
	}
}
