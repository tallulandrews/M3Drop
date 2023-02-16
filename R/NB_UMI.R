## Cell-type DANB ##
# fit variances with contribution from global mean~var relationship for small clusters. 
#Copyright (c) 2015, 2016 Genome Research Ltd .
#Author : Tallulah Andrews <tallulandrews@gmail.com>
#This file is part of M3Drop.

#M3Drop is free software : you can redistribute it and/or modify it under
#the terms of the GNU General Public License as published by the Free Software
#Foundation; either version 2 of the License, or (at your option) any later
#version.

#This program is distributed in the hope that it will be useful, but WITHOUT
#ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
#FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

#You should have received a copy of the GNU General Public License along with
#this program . If not , see <http://www.gnu.org/licenses/>.

#### Fitting #####
hidden_calc_vals <- function(counts) {
        if (sum(counts < 0) >0) {stop("Expression matrix contains negative values! Please provide raw UMI counts!")}
	if ( sum(counts >= 1) != sum(counts > 0) ) {stop("Error: Expression matrix is not integers! Please provide raw UMI counts.")}
#        if (sum(!is.integer(counts)) >0) {stop("Expression matrix is not integers! Please provide a matrix (not data.frame) raw UMI counts!")}
	if (is.null(rownames(counts))) { rownames(counts) <- as.character(1:nrow(counts)) }

        tjs <- Matrix::rowSums(counts, na.rm=T) # Total molecules/gene
	no_detect <- sum(tjs <= 0)
	if (no_detect > 0) {
		stop(paste("Error: contains",no_detect, "undetected genes."))
	}
        tis <- Matrix::colSums(counts, na.rm=T) # Total molecules/cell
	if (sum(tis <= 0) > 0) {stop("Error: all cells must have at least one detected molecule.")}
        djs <- ncol(counts)-Matrix::rowSums(counts > 0, na.rm=T) # Observed Dropouts per gene
        dis <- nrow(counts)-Matrix::colSums(counts > 0, na.rm=T) # Observed Dropouts per cell
        nc <- length(counts[1,]) # Number of cells
        ng <- length(counts[,1]) # Number of genes
        total <- sum(tis, na.rm=T) # Total molecules sampled
        return(list(tis=tis, tjs=tjs, dis=dis, djs=djs, total=total,nc=nc,ng=ng));
}

NBumiConvertToInteger <- function(mat) {
        mat <- ceiling(as.matrix(mat))
        storage.mode(mat) <- "integer"
        mat <- mat[Matrix::rowSums(mat) > 0,]
        return(mat)
}

NBumiFitModel <- function(counts) {
	vals <- hidden_calc_vals(counts);
#	mus <- (vals$tjs) %*% t(vals$tis/vals$total)
	
#	max_size <- max(mus)^2;
	min_size <- 10^-10;
	my_rowvar <- sapply(1:nrow(counts), function(i){
				mu_is <- vals$tjs[i]*vals$tis/vals$total
				var(as.vector(unlist(counts[i,]))-mu_is)
			})
#	my_rowvar <- vector(length=nrow(counts));
#	for (i in 1:nrow(counts)) {
#		mu_is <- vals$tjs[i]*vals$tis/vals$total
#		my_rowvar[i] <- var(as.vector(unlist(counts[i,]))-mu_is)
#	}
		
	size <- vals$tjs^2*(sum(vals$tis^2)/vals$total^2)/((vals$nc-1)*my_rowvar-vals$tjs) # for this to work with sparse matrices might need to implement in C
	max_size <- 10*max(size);
	size[size < 0] <- max_size;
	size[size < min_size] <- min_size;
#	size[size > max_size] <- max_size;

	return(list(var_obs=my_rowvar, sizes=size, vals=vals))
}

NBumiFitBasicModel <- function(counts) {
	vals <- hidden_calc_vals(counts);
	mus <- vals$tjs/vals$nc
	gm <- Matrix::rowMeans(counts)
	v <- Matrix::rowSums((counts-gm)^2)/(ncol(counts)-1)
	#v <- matrixStats::rowVars(counts)
	errs <- v < mus;
	v[errs] <- mus[errs]+10^-10;
	size <- mus^2/(v-mus)
	max_size <- max(mus)^2;
	size[errs] <- max_size;
	my_rowvar <- vector(length=nrow(counts));
	for (i in 1:nrow(counts)) {
		my_rowvar[i] <- var(counts[i,]-mus[i])
	}
	#mus_mat <- matrix(rep(mus, times=vals$nc), ncol=vals$nc, byrow=FALSE)
	return(list(var_obs=my_rowvar, sizes=size, vals=vals))
}

NBumiCheckFit <- function(counts, fit, suppress.plot=FALSE) {
	vals <- fit$vals;
	size_mat <- matrix(rep(fit$sizes, times=vals$nc), ncol=vals$nc, byrow=FALSE)
#	vs <- vector(length=nrow(counts));
	row_ps <- vector(length=nrow(counts))
	col_ps <- rep(0, times=ncol(counts))
	for (i in 1:nrow(counts)) {
		mu_is <- vals$tjs[i]*vals$tis/vals$total
		p_is <- (1+mu_is/fit$size[i])^(-fit$size[i]);
		row_ps[i] <- sum(p_is);
		col_ps <- col_ps+p_is;
#		vs[i] <- mean(mu_is+mu_is*mu_is/fit$size[i])
	}
#	thing <- fit$var_obs
#	plot(thing, vs, log="xy", xlab="Observed", ylab="Expected", main="Gene-specific Variance")
#	abline(a=0, b=1, col="red")

#	exp_ps <- (1+fit$mus/size_mat)^(-size_mat)
	if (!suppress.plot) {
		plot(vals$djs, row_ps, xlab="Observed", ylab="Fit", main="Gene-specific Dropouts")
		abline(a=0, b=1, col="red")

		plot(vals$dis, col_ps, xlab="Observed", ylab="Expected", main="Cell-specific Dropouts")
		abline(a=0, b=1, col="red")
	}

	invisible(list(gene_error = sum((vals$djs-row_ps)^2), cell_error = sum((vals$dis-col_ps)^2), rowPs=row_ps, colPs=col_ps));
}

NBumiCheckFitFS <- function(counts, fit, suppress.plot=FALSE) {
	vals <- fit$vals;
	size_coeffs <- NBumiFitDispVsMean(fit, suppress.plot=TRUE) 
	smoothed_size <- exp( size_coeffs[1] + size_coeffs[2]*log(vals$tjs/vals$nc) )
#	size_mat <- matrix(rep(smoothed_size, times=vals$nc), ncol=vals$nc, byrow=FALSE)

# Matrix compatibility
#	vs <- vector(length=nrow(counts));
	row_ps <- vector(length=nrow(counts))
	col_ps <- rep(0, times=ncol(counts))
	for (i in 1:nrow(counts)) {
		mu_is <- vals$tjs[i]*vals$tis/vals$total
		p_is <- (1+mu_is/smoothed_size[i])^(-smoothed_size[i]);
		row_ps[i] <- sum(p_is);
		col_ps <- col_ps+p_is;
#		vs[i] <- mean(mu_is+mu_is*mu_is/fit$size[i])
	}
#	thing <- fit$var_obs;
#	plot(thing, vs, log="xy", xlab="Observed", ylab="Expected", main="Gene-specific Variance")
#	abline(a=0, b=1, col="red")

#	exp_ps <- (1+fit$mus/size_mat)^(-size_mat)
	if (!suppress.plot) {
		plot(vals$djs, row_ps, xlab="Observed", ylab="Fit", main="Gene-specific Dropouts")
		abline(a=0, b=1, col="red")

		plot(vals$dis, col_ps, xlab="Observed", ylab="Expected", main="Cell-specific Dropouts")
		abline(a=0, b=1, col="red")
	}

	invisible(list(gene_error = sum((vals$djs-row_ps)^2), cell_error = sum((vals$dis-col_ps)^2), rowPs=row_ps, colPs=col_ps));
}


NBumiCompareModels <- function(counts, size_factor=(Matrix::colSums(counts)/median(Matrix::colSums(counts)))) {
	if (max(counts) < max(size_factor)) {
		stop("Error: size factors are too large");
	}
	norm <- NBumiConvertToInteger(t(t(counts)/size_factor));
	counts <- counts[rownames(counts) %in% rownames(norm),]; # necessary for plotting
	fit_adjust <- NBumiFitModel(counts);
	fit_basic <- NBumiFitBasicModel(norm);
	check_adjust <- NBumiCheckFitFS(counts, fit_adjust, suppress.plot=TRUE)
	check_basic <- NBumiCheckFitFS(norm, fit_basic, suppress.plot=TRUE)
	nc = fit_adjust$vals$nc

	plot( fit_adjust$vals$tjs/nc, fit_adjust$vals$djs/nc, col="black", pch=16, xlab="Expression", log="x", ylab= "Dropout Rate", cex=0.75)
#	points( fit_adjust$vals$tjs/nc, fit_basic$vals$djs/nc, col="black", pch=16, cex=0.75)
	points( fit_adjust$vals$tjs/nc, check_adjust$rowPs/nc, col="goldenrod1", pch=16, cex=0.5 )
	points( fit_adjust$vals$tjs/nc, check_basic$rowPs/nc, col="purple", pch=16, cex=0.5 )

	err_adj <- sum(abs(check_adjust$rowPs/nc-fit_adjust$vals$djs/nc))
#	err_bas <- sum(abs(Matrix::rowSums(check_basic$exp_ps)/nc-fit_basic$vals$djs/nc))
	err_bas <- sum(abs(check_basic$rowPs/nc-fit_adjust$vals$djs/nc))
	legend("bottomleft", paste(c("Depth-Adjusted\nError:", "Basic\nError:"), round(c(err_adj, err_bas)), c("\n","\n")), col=c("goldenrod1","purple"), pch=16, bty="n", cex=0.75)
	out <- c(err_adj, err_bas); names(out) <- c("Depth-Adjusted", "Basic");
	return(list(errors=out, basic_fit=fit_basic, adjusted_fit=fit_adjust))
}

obsolete__fit_size_to_drop <- function(obs, mu_vec, max_size, min_size=10^-10, convergence=0.001) {
	step_size <- 1;
        size_fit <- 1;
        last <- 0;
        last2 <- 0;
        for (iter in 1:1000) {
        	if (size_fit < min_size) {size_fit <- min_size;}
                diff <- sum( (1+mu_vec/size_fit)^(-size_fit) ) - obs;
                if (abs(diff) < convergence) {return(size_fit)}
                if (diff > 0) {
                	size_fit <- size_fit+step_size;
                        if (last < 0) {
                                step_size <- step_size/2;
                        } else {
                                if (last2 > 0){step_size <- step_size*2;}
                        }
                } else {
                        size_fit <- size_fit-step_size;
                        if (last > 0) {
                                step_size <- step_size/2;
                        } else {
                                if (last2 < 0){step_size <- step_size*2;}
                        }
                }
                last2 <- last; last <- diff;
                if (size_fit > max_size) {return(max_size);}
       }
       warning("Fitting size did not converge.");
       return(size_fit);
}

bg__fit_size_to_var <- function(obs, mu_vec, max_size, min_size=10^-10, convergence=0.001) {
	step_size <- 1;
        size_fit <- 1;
        last <- 0;
        last2 <- 0;
        for (iter in 1:1000) {
        	if (size_fit < min_size) {size_fit <- min_size;}
		expect <- mu_vec+mu_vec*mu_vec/size_fit;
                diff <- sum( expect )/(length(mu_vec)-1) - obs;
                if (abs(diff) < convergence) {return(size_fit)}
                if (diff > 0) {
                	size_fit <- size_fit+step_size;
                        if (last < 0) {
                                step_size <- step_size/2;
                        } else {
                                if (last2 > 0){step_size <- step_size*2;}
                        }
                } else {
                        size_fit <- size_fit-step_size;
                        if (last > 0) {
                                step_size <- step_size/2;
                        } else {
                                if (last2 < 0){step_size <- step_size*2;}
                        }
                }
                last2 <- last; last <- diff;
                if (size_fit > max_size) {return(max_size);}
       }
       warning("Fitting size did not converge.");
       return(size_fit);
}



#### Dispersion vs Mean & Feature Selection ####
NBumiFitDispVsMean <- function(fit, suppress.plot=TRUE) {
	vals <- fit$vals;
	size_g <- fit$sizes
	forfit <- fit$size < max(size_g) & vals$tjs > 0 & size_g > 0
	higher <- log(vals$tjs/vals$nc)/log(2) > 4; # As per-Grun et al.
	if (sum(higher == TRUE) > 2000) {
		forfit = higher & forfit;
	}
	reg <- lm( log(size_g[forfit])~ log((vals$tjs/vals$nc)[forfit]) )
	if (!suppress.plot) {
		plot(log( (vals$tjs/vals$nc)[forfit] ), log(size_g[forfit]), xlab="Log Mean Expression", ylab="Log Size" )
		abline(reg, col="red")
	}
	return(reg$coefficients)
}

hidden_shift_size <- function(mu_all, size_all, mu_group, coeffs) {
	b <- log(size_all)-coeffs[2]*log(mu_all)
	size_group <- exp(coeffs[2]*log(mu_group)+b)
	return(size_group)
}

obsolete__nbFeatureSelectionHighVarDist2Med <- function(fit, window_size=1000) {
	vals <- fit$vals;
	mean_order <- order(vals$tjs);
	obs_mean <- vals$tjs[mean_order]/vals$nc
	fit_disp <- fit$sizes[mean_order]
	names(fit_disp) <- names(vals$tjs)
	keep <- fit_disp < max(fit_disp)
	fit_disp <- log(fit_disp[keep])
	obs_mean <- obs_mean[keep]
	if (window_size %% 2 == 0) {
		flank <- window_size/2;
	} else {
		flank <- (window_size+1)/2;
	}
	dist_from_med <- function(x) {
		low <- max(1, x-flank);
		high <- min(length(fit_disp), x+flank);
		return(fit_disp[x]-median(fit_disp[seq(from=low, to=high, by=1)]))
	}
	score <- sapply(1:length(fit_disp), dist_from_med)
	return(sort(score))
}

NBumiFeatureSelectionHighVar <- function(fit) {
	# Global mean-variance
	vals <- fit$vals;
	coeffs <- NBumiFitDispVsMean(fit, suppress.plot=TRUE);
	exp_size <- exp( coeffs[1] + coeffs[2]*log(vals$tjs/vals$nc) )
	res <- log(fit$size) - log(exp_size);
	return(sort(res))
}

obsolete__nbFeatureSelectionDropouts <- function(fit) {
	# Gene-specific variance, mean
	vals <- fit$vals;
	size_mat <- matrix(rep(fit$sizes, times=vals$nc), ncol=vals$nc, byrow=F)

	droprate_exp <- vector(length=vals$ng)
	droprate_exp_err <- vector(length=vals$ng)
	for (i in 1:vals$ng) {
		mu_is <- vals$tjs[i]*vals$tis/vals$total
		p_is <- (1+mu_is/vals$size[i])^(-vals$size[i]);
		p_var_is <- p_is*(1-p_is);
		droprate_exp[i] <- sum(p_is)/vals$nc;
		droprate_exp_err[i] <- sqrt(sum(p_var_is)/(vals$nc^2));
	}
#	Exp_p <- (1+fit$mus/size_mat)^(-size_mat)
#	Exp_p_var <- Exp_p*(1-Exp_p)

#	droprate_exp <- Matrix::rowSums(Exp_p)/vals$nc
	droprate_exp[droprate_exp < 1/vals$nc] <- 1/vals$nc # Is this necessary?
#	droprate_exp_err <- sqrt(Matrix::rowSums(Exp_p_var)/(vals$nc^2))
	droprate_obs <- vals$djs/vals$nc
	droprate_obs_err <- sqrt(droprate_obs*(1-droprate_obs)/vals$nc)

	diff <- droprate_obs-droprate_exp
	combined_err <- droprate_exp_err
	Zed <- diff/combined_err
        pvalue <- pnorm(Zed, lower.tail = FALSE)
	names(pvalue) <- names(vals$tjs)
	return(sort(pvalue))
}

NBumiFeatureSelectionCombinedDrop <- function(fit, ntop=NULL, method="fdr", qval.thresh=2, suppress.plot=TRUE) {
	# Global mean-variance, gene-specific mean
	vals <- fit$vals;

	coeffs <- NBumiFitDispVsMean(fit, suppress.plot=TRUE);
	exp_size <- exp( coeffs[1] + coeffs[2]*log(vals$tjs/vals$nc) )

	#size_mat <- matrix(rep(exp_size, times=vals$nc), ncol=vals$nc, byrow=F)
	#Exp_p <- (1+fit$mus/size_mat)^(-size_mat)
	#Exp_p_var <- Exp_p*(1-Exp_p)

	droprate_exp <- vector(length=vals$ng)
	droprate_exp_err <- vector(length=vals$ng)
	for (i in 1:vals$ng) {
		mu_is <- vals$tjs[i]*vals$tis/vals$total
		p_is <- (1+mu_is/exp_size[i])^(-exp_size[i]);
		p_var_is <- p_is*(1-p_is);
		droprate_exp[i] <- sum(p_is)/vals$nc;
		droprate_exp_err[i] <- sqrt(sum(p_var_is)/(vals$nc^2));
	}

	#droprate_exp <- Matrix::rowSums(Exp_p)/vals$nc
	droprate_exp[droprate_exp < 1/vals$nc] <- 1/vals$nc # Is this necessary?
	#droprate_exp_err <- sqrt(Matrix::rowSums(Exp_p_var)/(vals$nc^2))
	droprate_obs <- vals$djs/vals$nc
	droprate_obs_err <- sqrt(droprate_obs*(1-droprate_obs)/vals$nc)

	diff <- droprate_obs-droprate_exp
	combined_err <- sqrt(droprate_exp_err^2+droprate_obs_err^2)
	#alt_err <- sd(diff)/sqrt(vals$nc);
	Zed <- diff/combined_err
        pvalue <- pnorm(Zed, lower.tail = FALSE)
	names(pvalue) <- names(vals$tjs)
	reorder <- order(pvalue, droprate_exp-droprate_obs) # deal with ties (e.g. lots of zero p-values)

	out <- pvalue[reorder]
	diff <- diff[reorder]
	qval <- p.adjust(out, method=method)
	if (is.null(ntop)) {
		out <- out[qval < qval.thresh]
		diff <- diff[qval < qval.thresh]
		qval <- qval[qval < qval.thresh]
	} else {
		out <- out[1:ntop]
		diff <- diff[1:ntop]
		qval <- qval[1:ntop]
	}
	outTABLE <- data.frame(Gene=names(out), effect_size=diff, p.value=out, q.value=qval)

	if (!suppress.plot) {
		xes <- log10(vals$tjs/vals$nc);
        	dens.col <- densCols(xes, droprate_obs, colramp=colorRampPalette(c("grey75","black")))
		plot(xes, droprate_obs, col=dens.col, pch=16, xlab="", ylab="")
		title(ylab="Dropout Rate", line=2)
		title(xlab="log10(expression)", line=2)
		toplot = names(out)
		toplot = names(vals$tjs) %in% toplot
		points(xes[toplot], droprate_obs[toplot], col="darkorange", pch=16)
		points(xes, droprate_exp, col="dodgerblue", pch=16, cex=1)
	}
	return(outTABLE)
}

PoissonUMIFeatureSelectionDropouts <- function(fit) {
	# Poisson distribution
	vals <- fit$vals
	#Exp_p <- exp(-fit$mus)
	#Exp_p_var <- Exp_p*(1-Exp_p)
	droprate_exp <- vector(length=vals$ng)
	droprate_exp_err <- vector(length=vals$ng)
	for (i in 1:vals$ng) {
		mu_is <- vals$tjs[i]*vals$tis/vals$total
		p_is <- exp(-mu_is);
		p_var_is <- p_is*(1-p_is);
		droprate_exp[i] <- sum(p_is)/vals$nc;
		droprate_exp_err[i] <- sqrt(sum(p_var_is)/(vals$nc^2));
	}

	#droprate_exp <- Matrix::rowSums(Exp_p)/vals$nc
	droprate_exp[droprate_exp < 1/vals$nc] <- 1/vals$nc # Is this necessary?
	#droprate_exp_err <- sqrt(Matrix::rowSums(Exp_p_var)/(vals$nc^2))
	droprate_obs <- vals$djs/vals$nc
	droprate_obs_err <- sqrt(droprate_obs*(1-droprate_obs)/vals$nc)

	diff <- droprate_obs-droprate_exp
	combined_err <- droprate_exp_err
	Zed <- diff/combined_err
        pvalue <- pnorm(Zed, lower.tail = FALSE)
	names(pvalue) <- names(vals$tjs)
	return(sort(pvalue))
}

#PoissonUMIFeatureSelectionCV <- function(counts, fit) {
#	# Poisson distribution
#	vals <- fit$vals
#	g_means <- vals$tjs/vals$nc
#	CVexp <- 1/sqrt(g_means) # expected
#	CVobs <-  sqrt(apply(counts-fit$mus, 1, var))/g_means;
#
#	CVexp_log <- log(CVexp)
#	CVobs_log <- log(CVobs)
#	mean_log <- log(g_means)
#
#	res <- CVobs_log-CVexp_log;
#	res_sd <- sd(res);
#
#	# Variance of sample variance from http://mathworld.wolfram.com/SampleVarianceDistribution.html
#	# Kenney and Keeping 1951, p. 164; Rose and Smith 2002, p. 264
#	u4 <- g_means*(1+3*g_means) # central moment
#	u2 <- g_means
#	N <- vals$nc
#	V <- (N-1)^2/N^3*u4-(N-1)*(N-3)/N^3*u2^2
#
#	Zed <- (CVobs-CVexp)/sqrt(V)
#        pvalue <- pnorm(Zed, lower.tail = FALSE)
#	names(pvalue) <- names(vals$tjs)
#	return(sort(pvalue))
#}
#### Differential Expression #####

# HERE HERE HERE - editing for Matrix sparse-matrix compatibility #

unfinished__nbGroupDE <- function(counts, fit, groups) {
	# Global mean-variance, gene-specific variance & mean
	vals <- fit$vals;
	size_g <- fit$sizes
	group_specific_factor <- aggregate(t(counts), by=list(groups), sum)
	rownames(group_specific_factor) <- group_specific_factor[,1]
	group_specific_factor <- group_specific_factor[,-1]
	group_specific_factor <- t(group_specific_factor)
	group_specific_tjs <- group_specific_factor; 
	group_total <- Matrix::colSums(group_specific_tjs)
	group_specific_factor <- t(t(group_specific_tjs)/group_total) / (vals$tjs/vals$total) # Relative expression level in group vs across whole dataset.
	group_specific_factor[vals$tjs == 0,] <- rep(1, length(group_specific_factor[1,]))

	coeffs <- NBumiFitDispVsMean(fit, suppress.plot=TRUE);

	pvals <- vector(length=vals$ng)
	for (j in 1:vals$ng) {
		prob_null <- 0;
		prob_test <- 0;
		for (i in 1:vals$nc) {
			group <- which(colnames(group_specific_factor) == groups[i])
			mu_ji <- vals$tjs[j]*vals$tis[i]/vals$total
			p1 <- dnbinom(counts[j,i], mu=mu_ji, size=size_g[j], log=TRUE)
			group_mu <- mu_ji*group_specific_factor[j,group]
			shifted_size <- hidden_shift_size(mu_ji,size_g[j],group_mu,coeffs);
			p2 <- dnbinom(counts[j,i], mu=group_mu, size=shifted_size, log=TRUE)
			# Collect ps across i
			prob_null <- prob_null+p1
			prob_test <- prob_test+p2
		}
		# Perform log-likehood test for this gene
		D <- -2*(prob_null-prob_test)
		df <- length(unique(groups))-1
		pval <- pchisq(D, df=df, lower.tail=FALSE)
		pvals[j] <- pval;
	}
	# Format output.
	output <- cbind(group_specific_factor, pvals, p.adjust(pvals, method="fdr"));
	rownames(output) <- rownames(counts);
	colnames(output) <- c(colnames(group_specific_factor), "p.value","q.value")
	return(output);
}

hidden_nbImputeZeros <- function(counts, fit) {
	vals <- fit$vals
	for (i in 1:nrow(counts)) {
		mu_is <- vals$tjs[i]*vals$tis/vals$total
		low <- counts[i,] < mu_is
		counts[i,low] <- counts[i,low]+mu_is[low]
	}
	return(counts);
}

NBumiImputeNorm <- function(counts, fit, total_counts_per_cell=median(fit$vals$tis)) {
	# find p-value for current fit NB Umi
	# adjust parameters of the NB
	# determine value of equivalent p-value under new NB
	coeffs <- NBumiFitDispVsMean(fit, suppress.plot=TRUE);
	vals <- fit$vals;
	norm <- counts;
	normed_ti <- total_counts_per_cell;
	normed_mus <- vals$tjs/vals$total;
	# TODO modify so adjusts size to new mean expression level
	#normed_size <- fit$size;
	for (i in 1:nrow(counts)) {
		mu_is <- vals$tjs[i]*vals$tis/vals$total;
		p_orig <- pnbinom(counts[i,], mu=mu_is, size=fit$sizes[i]);
		new_size <- hidden_shift_size(mean(mu_is), fit$sizes[i], normed_mus[i]*normed_ti, coeffs)
		normed <- qnbinom(p_orig, mu=normed_mus[i]*normed_ti, size=new_size);
		norm[i,] <- normed;
	}
	return(norm);
}


NBumiConvertData <- function(input, is.log=FALSE, is.counts=FALSE, pseudocount=1) {
	type <- class(input)[1]
	lognorm <- NULL
	counts <- NULL
	if (type == "SCESet") {
		# Old scater
		lognorm <- scater::exprs(input)
		counts <- counts(input)

	} else if (type == "SingleCellExperiment") {
		# New scater
		c <- which(names(assays(input)) == "counts")
		ln <- which(names(assays(input)) == "logcounts")
		if (length(c) > 0) {
			counts <- assays(input)[[c]]
		} else if (length(ln) > 0) {
			lognorm <- assays(input)[[ln]]
		} else {
			stop("Error: Recognized SingleCellExperiment object but cannot find either counts or lognorm expression.")
		}
	} else if (type == "CellDataSet" | type == "ExpressionSet") {
		# monocle
		if (is.log) {
			lognorm <- Biobase::exprs(input)
		} else {
			counts <- Biobase::exprs(input)
		}
	} else if (type == "seurat") {
		# Seurat
		counts <- input@raw.data
	} else if (type == "matrix" | 
		   type == "data.frame" | 
		   type == "dgCMatrix" | 
		   type == "data.table" |
		   type == "DataTable" |
		   type == "DataFrame" |
		   type == "AnnotatedDataFrame") {
		if (type != "dgCMatrix") {
			input <- as.matrix(input)
		}

		if (is.log) {
			lognorm <- input;
		} else if (is.counts) {
			counts <- input
		} else {
			norm <- input;
		}
	} else {
		stop(paste("Error: Unrecognized input format :", type))
	}

	# Prefer raw counts to lognorm
	remove_undetected_genes <- function(mat) {
		no_detect <- Matrix::rowSums(mat > 0, na.rm=T) == 0;
		print(paste("Removing ",sum(no_detect), "undetected genes."))
		return(mat[!no_detect,])
	}
	# CPM transform raw counts
	if (!is.null(dim(counts))) {
		counts <- ceiling(counts)
		counts <- remove_undetected_genes(counts);
		return( counts )
	}
	# If normalized rescale by number of detected genes
	if (!is.null(dim(lognorm))) {
		norm <- 2^lognorm-pseudocount
	}
	if (!is.null(dim(norm))) {
		sf <- apply(norm, 2, function(x){min(x[x>0])}) #colSums(norm)
		sf <- 1/sf
		#detected <- colSums(norm > 0);
		#detected <- detected/median(detected);
		#libsize <- detected*median(sf)
		counts <- t( t(norm)/sf) #*libsize) )
		counts <- ceiling(counts)
		counts <- remove_undetected_genes(counts);
		return(counts)
	} 
}

