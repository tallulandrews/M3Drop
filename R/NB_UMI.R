#### Fitting #####
hidden_calc_vals <- function(counts) {
        if (sum(counts < 0) >0) {stop("Expression matrix contains negative values! Please provide raw UMI counts!")}
        if (sum(!is.integer(counts)) >0) {stop("Expression matrix is not integers! Please provide raw UMI counts!")}

        p = apply(counts,1,function(x){y = x[!is.na(x)]; sum(y==0)/length(y)});
        s = rowMeans(counts, na.rm=T);

        tis = colSums(counts, na.rm=T) # Total molecules/cell
        tjs = rowSums(counts, na.rm=T) # Total molecules/gene
        djs = rowSums(counts == 0, na.rm=T) # Observed Dropouts per gene
        dis = colSums(counts == 0, na.rm=T) # Observed Dropouts per cell
        nc = length(counts[1,]) # Number of cells
        ng = length(counts[,1]) # Number of genes
        total = sum(tis, na.rm=T) # Total molecules sampled
        return(list(tis=tis, tjs=tjs, djs=djs, dis=dis, total=total,nc=nc,ng=ng));
}

NBumiFitModel <- function(counts) {
	vals <- hidden_calc_vals(counts);
	mus <- (vals$tjs) %*% t(vals$tis/vals$total)
	
	convergence <- 0.001
	max_size <- max(mus)^2;
	min_size <- 10^-10;
	size <- sapply(1:vals$ng,function(j) {
		 bg__fit_size_to_var(var(counts[j,]-mus[j,]), mus[j,], max_size=max_size, min_size=10^-10, convergence=0.001)
		})
	return(list(mus=mus, sizes=size, vals=vals))

}

NBumiCheckFit <- function(counts, fit) {
	vals <- fit$vals;
	size_mat <- matrix(rep(fit$sizes, times=vals$nc), ncol=vals$nc, byrow=F)
	vs <- fit$mus+fit$mus*fit$mus/size_mat
	thing <- apply(counts-fit$mus, 1, var)
	plot(thing, rowMeans(vs), log="xy", xlab="Observed", ylab="Expected", main="Gene-specific Variance")
	abline(a=0, b=1, col="red")

	exp_ps <- (1+fit$mus/size_mat)^(-size_mat)
	plot(vals$djs, rowSums(exp_ps), xlab="Observed", ylab="Fit", main="Gene-specific Dropouts")
	abline(a=0, b=1, col="red")

	plot(vals$dis, colSums(exp_ps), xlab="Observed", ylab="Expected", main="Cell-specific Dropouts")
	abline(a=0, b=1, col="red")

#	plot(vals$tjs/vals$nc, rowMeans(fit$mus), xlab="Observed", ylab="Fit", main="Gene Mean Expression")
#	abline(a=0, b=1, col="red")

#	plot(vals$tis, colSums(fit$mus), xlab="Observed", ylab="Fit", main="Cell Count Total")
#	abline(a=0, b=1, col="red")
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

NBumiFeatureSelectionHighVar <- function(fit, window_size=1000) {
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
	score <- -1*sapply(1:length(fit_disp), dist_from_med)
	return(rev(sort(score)))
}

NBumiFeatureSelectionDropouts <- function(fit) {
	vals <- fit$vals;
	size_mat <- matrix(rep(fit$sizes, times=vals$nc), ncol=vals$nc, byrow=F)
	Exp_p <- (1+fit$mus/size_mat)^(-size_mat)
	Exp_p_var <- Exp_p*(1-Exp_p)

	droprate_exp <- rowSums(Exp_p)/vals$nc
	droprate_exp[droprate_exp < 1/vals$nc] <- 1/vals$nc # Is this necessary?
	droprate_exp_err <- sqrt(rowSums(Exp_p_var)/(vals$nc^2))
	droprate_obs <- vals$djs/vals$nc
	droprate_obs_err <- sqrt(droprate_obs*(1-droprate_obs)/vals$nc)

	diff <- droprate_obs-droprate_exp
	combined_err <- droprate_exp_err
	Zed <- diff/combined_err
        pvalue <- pnorm(Zed, lower.tail = FALSE)
	names(pvalue) <- names(vals$tjs)
	return(sort(pvalue))
}

PoissonUMIFeatureSelectionDropouts <- function(vals, fit) {
#	size_mat <-  matrix(rep(fit$sizes, times=vals$nc), ncol=vals$nc, byrow=F)
	Exp_p <- exp(-fit$mus)
	Exp_p_var <- Exp_p*(1-Exp_p)

	droprate_exp <- rowSums(Exp_p)/vals$nc
	droprate_exp[droprate_exp < 1/vals$nc] <- 1/vals$nc # Is this necessary?
	droprate_exp_err <- sqrt(rowSums(Exp_p_var)/(vals$nc^2))
	droprate_obs <- vals$djs/vals$nc
	droprate_obs_err <- sqrt(droprate_obs*(1-droprate_obs)/vals$nc)

	diff <- droprate_obs-droprate_exp
	combined_err <- droprate_exp_err
	Zed <- diff/combined_err
        pvalue <- pnorm(Zed, lower.tail = FALSE)
	names(pvalue) <- names(vals$tjs)
	return(sort(pvalue))
}
#### Differential Expression #####

NBumiGroupDE <- function(counts, fit, groups, mean2disp_coeffs) {
	vals <- fit$vals;
	size_g <- fit$sizes
	group_specific_factor <- aggregate(t(counts), by=list(groups), sum)
	rownames(group_specific_factor) <- group_specific_factor[,1]
	group_specific_factor <- group_specific_factor[,-1]
	group_specific_factor <- t(group_specific_factor)
	group_specific_tjs <- group_specific_factor; 
	group_total <- colSums(group_specific_tjs)
	group_specific_factor <- t(t(group_specific_tjs)/group_total) / (vals$tjs/vals$total) # Relative expression level in group vs across whole dataset.
	group_specific_factor[vals$tjs == 0,] <- rep(1, length(group_specific_factor[1,]))

	coeffs <- NBumiFitDispVsMean(fit, suppress.plot=TRUE);

	pvals <- vector(length=vals$ng)
	for (j in 1:vals$ng) {
		prob_null <- 0;
		prob_test <- 0;
		for (i in 1:vals$nc) {
			group <- which(colnames(group_specific_factor) == groups[i])
			p1 <- dnbinom(counts[j,i], mu=fit$mus[j,i], size=size_g[j], log=TRUE)
			group_mu <- fit$mus[j,i]*group_specific_factor[j,group]
			shifted_size <- hidden_shift_size(fit$mus[j,i],size_g[j],group_mu,coeffs);
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
