hidden__pca_fs <- function(expr_mat, pcs = c(1,2)) {
	pca <- prcomp(log(expr_mat+1)/log(2));
	if (length(pcs) > 1) {
		score <- Matrix::rowSums(abs(pca$x[,pcs]))
	} else {
		score <- abs(pca$x[,pcs])
	}
	names(score) = rownames(expr_mat);
	return(sort(-score))
}

irlbaPcaFS <- function(expr_mat, pcs=c(2,3)) {
	#require("irlba")
	#require("Matrix")
	norm <- expr_mat
	nz_genes <- which(Matrix::rowSums(norm) != 0)
	norm[nz_genes,] <- log(norm[nz_genes,] + 1)/log(2)
	# Create sparse Matrix

	gene_names <- rownames(norm);
	norm <- Matrix::Matrix(norm, sparse=TRUE)
	rownames(norm) <- gene_names
	#indices <- which(norm > 0, arr.ind=TRUE)
	#vals <- norm[indices]
	nc <- ncol(norm)

	#norm <- Matrix::sparseMatrix(i = indices[,1], j=indices[,2], x=vals)
	# Calculate the variance across genes without converting to a dense
	# matrix:
	expression_means <- Matrix::rowMeans(norm)
	expression_vars <- Matrix::rowMeans((norm - expression_means)^2)*(nc/(nc-1)) # sample variance.
	# Filter out genes that are constant across all cells:
	genes_to_keep <- expression_vars > 0
	norm <- norm[genes_to_keep,]
	expression_means <- expression_means[genes_to_keep]
	expression_vars <- expression_vars[genes_to_keep]
	# Here✬s how to take the top PCA loading genes, but using
	# sparseMatrix operations the whole time, using irlba.
	irlba_pca_res <- irlba::irlba(Matrix::t(norm), nu=0, center=expression_means, scale=sqrt(expression_vars), right_only=TRUE)$v
	row.names(irlba_pca_res) <- row.names(norm)
	if (length(pcs) > 1) {
		score <- Matrix::rowSums(abs(irlba_pca_res[, pcs]))
	} else {
		score <- abs(irlba_pca_res[, pcs])
	}
	names(score) = gene_names;
	return(sort(-score))
}


hidden__Grun_fs <- function(expr_mat, spikes) {
#	if (!is.factor(batches)) { batches <- factor(batches) }
#	batch_spikes <- sapply(levels(batches), function(b) {rowMeans(expr_mat[spikes,batches==b])})
	spike_T <- rowMeans(expr_mat[spikes,])
	beta2 <- apply(expr_mat, 2, function(c) {lm(c[spikes]~spike_T)})
	# Model 3
	beta3 <- expr_mat[spikes,]/spike_T
	fit_gamma <- function(x) {
		b = var(x)/mean(x)
		a = mean(x)/b
		return(c(a,1/b))
	}
	gammas <- apply(beta3, 1, fit_gamma)
	a_fit = lm(gammas[1,] ~ spike_T)
	b_fit = lm(gammas[2,] ~ spike_T)
	get_nb_params <- function(mu) {
		a = a_fit$coefficients[1]+a_fit$coefficients[2]*mu
		b = b_fit$coefficients[1]+b_fit$coefficients[2]*mu
		nb_mu = mu*a/b
		nb_size = a
		return(c(nb_mu, nb_size))
	}
	gammas <- log(gammas)/log(2)
	a_fit_l = lm(gammas[1,] ~ spike_T)
	b_fit_l = lm(gammas[2,] ~ spike_T)
	ka = a_fit_l$coefficients[1]
	kb = b_fit_l$coefficients[1]
	fa = a_fit_l$coefficients[2]
	fb = b_fit_l$coefficients[2]
	get_nb_params_l <- function(mu) {
		r = 2^(ka+fa*(kb-ka)/(1+fa-fb))*mu^((fa)/(1+fa-fb))
		return(r)
	}
	# deconvolution
	min_fun <- function(bio, ns, u_g, r_g) {
		u_bio=bio[1]
		r_bio=bio[2]
		diff <- function(n) {
			ms = 1:(2*n)
			deconvolve = sum(sapply(ms, function(m) {dnbinom(n,mu=m,size=get_nb_params_l(m))*dnbinom(m,mu=u_bio, size=r_bio)}));
			abs(dnbinom(n, mu=u_g, size=r_g) - deconvolve)
		}
		sum(sapply(ns, diff));
	}
	get_bio_disp <- function(g) {
	#	g = expr_mat[1,]
		u_g = mean(g);
		v_g = var(g);
		if (v_g <= u_g) {v_g = u_g+10^-10}
		r_g = u_g^2/(v_g-u_g)
		stuff = optim(c(u_g, r_g), min_fun, ns=round(g), u_g=u_g, r_g=r_g, method="BFGS")
		return(stuff$par[2])
	}
	bio_disp <- apply(expr_mat, 1, get_bio_disp)
	names(bio_disp) = rownames(expr_mat);
	return(bio_disp)
}

# How does GiniClust do it?
hidden__ginifs_simple <- function(expr_mat) {
	#require("reldist")
	ginis <- apply(expr_mat, 1, reldist::gini)
	d <- rowMeans(expr_mat>0)
	reg <- lm(ginis~d) # almost perfect linear relation in UMI data
	score <- reg$res
	return(sort(-score));
}

giniFS <- function(expr_mat, suppress.plot=TRUE) {
	# GiniClust
	expr_mat <- expr_mat[Matrix::rowSums(expr_mat) > 0,]
	#require("reldist")
	ginis <- apply(expr_mat, 1, reldist::gini)
	max_expr <- apply(expr_mat, 1, max)
	max_expr <- log(max_expr+1)/log(2)
	fit = loess(ginis~max_expr)
	outliers = abs(fit$residuals)
	outliers = outliers > quantile(outliers, probs=0.75)
	fit2 = loess(ginis[!outliers]~max_expr[!outliers])

	
	norm_ginis = rep(NA, times=length(ginis));
	norm_ginis[!outliers] = fit2$residuals;
	to_impute = which(is.na(norm_ginis))
	impute_loess <- function(i) {
		d = abs(max_expr-max_expr[i])
		closest = which(d[!outliers] == min(d[!outliers]))
		imputed = fit2$fitted[closest[1]]
		return(imputed)
	}
	fit_ginis = sapply(to_impute, impute_loess)
	norm_ginis[outliers] = ginis[outliers]-fit_ginis
	p = pnorm(norm_ginis, mean=mean(norm_ginis), sd=sd(norm_ginis), lower.tail=FALSE)
	names(p) = rownames(expr_mat)

	if (!suppress.plot) {
		my_cols <- colorRampPalette(c("grey75","black"))(20);
		plot(max_expr, ginis, pch=16, col=rev(my_cols)[cut(log(p), breaks=20)], xlab="Max Expression", ylab="Gini");
		legend("bottomleft", bty="n", c("Color indicates estimated p-value"), col="white", pch=16, cex=0.75)
		tmp <- order(max_expr);
		yes <-  rep(NA, times=length(ginis)); yes[!outliers]= fit2$fitted; yes[outliers] = fit_ginis
		lines(max_expr[tmp], yes[tmp], col="red", lwd=3)
	}

	return(sort(p));
}
corFS <- function(expr_mat, dir=c("both", "pos", "neg"), fdr=NULL) {
	# High memory
	expr_mat <- as.matrix(expr_mat);
	if (!is.null(fdr)) {
		cor_mat = Hmisc::rcorr(t(expr_mat), type="spearman")
		p_mat <- cor_mat$P
		cor_mat <- cor_mat$r
	} else {
		cor_mat <- cor(t(expr_mat), method="spearman")
	}
	diag(cor_mat) <- 0
	if (dir[1] == "both") {
		score = apply(cor_mat, 1, function(x) {sum(abs(min(x)), abs(max(x)))})
	} else if (dir[1] == "pos") {
		score = apply(cor_mat, 1, max)
	} else if (dir[1] == "neg") {
		score = abs(apply(cor_mat, 1, min))
	} else {
		stop("Unrecognized direction")
	}
	if (!is.null(fdr)) {
		sig <- matrix(p.adjust(p_mat, method="fdr"), ncol=ncol(p_mat))
		score <- score[apply(sig, 1, min, na.rm=T) < fdr]
	}
	rm(cor_mat)
	names(score) = rownames(expr_mat)
	return(sort(-score))
}

Consensus_FS <- function(counts, norm=NA, is.spike=rep(FALSE, times=nrow(counts)), pcs=c(2,3), include_cors=TRUE) {
	# Check input
	#if (!is.matrix(counts)) {
	#	counts <- as.matrix(counts);
	#}


	if (sum(dim(counts) != dim(norm)) > 0) {
		if (is.null(dim(norm))) {
			if (length(norm) < ncol(counts)) {
				# apply CPM
				sf <- Matrix::colSums(counts[!is.spike,])
			} else {
				sf <- norm
			}
			norm <- t(t(counts)/sf*median(sf));
		} else {
			stop("Error: counts and norm matrices must be of the same dimension");
		}
	}
	if (ncol(counts) < 1000) {
		row_vars <- (rowMeans(counts^2)-rowMeans(counts)^2);
	} else {
		row_vars <- Matrix::rowSums(counts);
	}
	invariant <- row_vars == 0;
	counts <- counts[!invariant,]
	norm <- norm[!invariant,]
	# apply FS methods
	# DANB
	fit <- NBumiFitModel(counts)
	DANB <- NBumiFeatureSelectionCombinedDrop(fit)
	DANB_var <- NBumiFeatureSelectionHighVar(fit)
	# HVG
	if (sum(is.spike == TRUE) > 10) {
		spikes <- rownames(norm)[is.spike]
		HVG <- BrenneckeGetVariableGenes(norm, spikes=spikes, fdr=2, suppress.plot=TRUE)
	} else {
		warning("Warning: insufficient spike-ins using all genes for HVG");
		HVG <- BrenneckeGetVariableGenes(norm, fdr=2, suppress.plot=TRUE)
	}
	# M3Drop
	m3drop <- M3DropFeatureSelection(norm, mt_method="fdr", mt_threshold=2, suppress.plot=TRUE)
	# Gini
	gini <- giniFS(norm)
	# pca
	pca_fs <- irlbaPcaFS(norm, pcs=pcs) 
	# cor
	if (include_cors) {
		cor_fs <- corFS(norm)
	} else {
		cor_fs <- rep(-1, times=nrow(norm));
		names(cor_fs) <- rownames(norm);
	}

	# sort by mean_rank of each gene.
	ranks <- 1:nrow(norm);
	ref_order <- rownames(counts); 
	out_table <- data.frame(
		DANB_drop = ranks[match(ref_order, DANB$Gene)],
		DANB_var = ranks[match(ref_order, names(DANB_var))],
		M3Drop = ranks[match(ref_order, m3drop$Gene)],
		HVG = ranks[match(ref_order, HVG$Gene)],
		PCA = ranks[match(ref_order, names(pca_fs))],
		Cor = ranks[match(ref_order, names(cor_fs))],
		Gini = ranks[match(ref_order, names(gini))]
		)
	if (!include_cors) {
		out_table$Cor <- rep(-1, times=nrow(out_table));
	}
	rownames(out_table) <- ref_order;
	consensus_score <- rowMeans(out_table, na.rm=TRUE);
	out_table <- out_table[order(consensus_score),]
	out_table$Cons <- 1:nrow(out_table) 

	return(out_table)
}
