weighted_pca <- function(m,w,Nf=1,e=0.001,max_steps=10000) {
	if (prod(dim(m) == dim(w))) {
		nrow = length(m[,1])
		ncol = length(m[1,])
		A = matrix(rep(0, times = nrow*Nf), nrow = nrow, ncol=Nf)
		B = matrix(rep(0, times = ncol*Nf), nrow = Nf, ncol=ncol)
		out <- .C("PCA_WGT", y=as.double(m), w=as.double(w), Ng=as.integer(length(m[,1])), Nc=as.integer(length(m[1,])), e=as.double(e), Nf = as.integer(Nf), max_steps=as.integer(max_steps),A = as.double(A), B = as.double(B) )
		
		# Format output
		out$y = matrix(out$y, nrow = nrow, ncol = ncol)
		out$w = matrix(out$w, nrow = nrow, ncol = ncol)
		out$A = matrix(out$A, nrow = nrow)
		out$B = matrix(out$B, ncol = ncol)
		colnames(out$y) = colnames(m)
		rownames(out$y) = rownames(m)
		colnames(out$B) = colnames(m)
		rownames(out$A) = rownames(m)
		colnames(out$A) = paste("PC",1:length(out$A[1,]),sep="")
		rownames(out$B) = paste("PC",1:length(out$B[,1]),sep="")
		out$used_steps = out$max_steps
		out$max_steps = max_steps

		return(out)
	} else {
		stop("Weights matrix must be same dimension as data matrix.");
	}
}

weighted_pca_v2 <- function (m,w,method="covariance", Nf=1) {
	# calculate weighted covariance or correlation matrix
	c_matrix = matrix(0,nrow=length(m[,1]), ncol=length(m[,1]));
	.C(as.double(m),as.double(w),as.integer(length(m[,1])), as.integer(length(m[1,])), as.double(c_matrix));
	# do eigen decomposition
	decomp <- eigen(c_matrix, symmetric=TRUE)
	# Calculate Loadings
	loadings = t(decomp$vectors[,1:Nf]) %*% m;
	return(PCs = decomp$vectors[,1:NF], loadings = loadings);

}
