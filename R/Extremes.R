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

# DE Genes functions

hidden__test_DE_K_equiv_raw <- function (expr_mat, fit=NA) {
	gene_info <- bg__calc_variables(expr_mat);
	if (is.na(fit)[1]) {
		fit <- bg__fit_MM(gene_info$p, gene_info$s);
	}
	p_obs <- gene_info$p;
	N <- length(expr_mat[1,]);
	p_err <- gene_info$p_stderr;
	S_mean <- gene_info$s
	S_err <- gene_info$s_stderr
	K_err <- fit$Kerr;
	K_equiv <- p_obs*S_mean/(1-p_obs);
	K_equiv_err <- p_obs/(1-p_obs)*S_err
		
	Z <- (K_equiv - fit$K)/sqrt(K_equiv_err^2+K_err^2); # high = shifted right, low = shifted left
	pval <- pnorm(Z, lower.tail=FALSE)
	effect_size <- K_equiv/fit$K;
	return(list(pval = pval, fold_change = effect_size))
}

bg__test_DE_K_equiv <- function(gene_info, fit=NA) {
	p_obs <- gene_info$p;
	always_detected <- p_obs==0
	p_obs[p_obs==0] <- min(p_obs[p_obs > 0])/2 # Here so that K_equiv for p_obs==0 will be 0 but K_equiv_err will not throw errors for such genes.
	p_err <- gene_info$p_stderr;
	S_mean <- gene_info$s
	S_err <- gene_info$s_stderr
	K_err <- fit$Kerr;
	K_obs <- fit$K
	K_equiv <- p_obs*S_mean/(1-p_obs);
	K_equiv_err <- abs(K_equiv)*sqrt((S_err/S_mean)^2 + (p_err/p_obs)^2)


	K_equiv_log <- log(K_equiv)
	thing <- K_equiv-K_equiv_err; thing[thing <= 0 ] <- 10^-100
	K_equiv_err_log <- abs(log(thing)-K_equiv_log) # Testing to the right so this is what we are interested in
#	K_equiv_err_log = K_equiv_err/K_equiv # This does not hold when K_equiv_err =~ K_equiv which is particularly problematic for lowly expressed genes
	K_equiv_err_log[K_equiv-K_equiv_err <= 0 ] <- 10^10
	K_obs_log <- log(fit$K)
	K_err_log <- sd(K_equiv_log-K_obs_log)/sqrt(length(K_equiv_log)) 
		
	Z <- (K_equiv_log - K_obs_log)/sqrt(K_equiv_err_log^2+K_err_log^2); # high = shifted right, low = shifted left
	pval <- pnorm(Z, lower.tail=FALSE)
	pval[always_detected] <- 1;
	effect_size <- K_equiv/fit$K;
	return(list(pval = pval, fold_change = effect_size))
}
# Use the fact that errors of proportions are well define by converting S to proportion detected equivalents?
hidden__test_DE_P_equiv <- function (expr_mat,  fit=NA) {
	gene_info <- bg__calc_variables(expr_mat);
	if (is.na(fit)[1]) {
		fit <- bg__fit_MM(gene_info$p, gene_info$s);
	}
	p_obs <- gene_info$p;
	N <- length(expr_mat[1,]);
	p_err <- gene_info$p_stderr;
	S_mean <- gene_info$s
	S_err <- gene_info$s_stderr
	K_err <- fit$Kerr;
	p_equiv <- fit$predictions;
	propagated_err_p_equiv <- p_equiv*sqrt(((S_err+K_err)/(S_mean+fit$K))^2+(K_err/fit$K)^2)
	fitted_err_p_equiv <- fit$fitted_err
	Z <- (p_equiv - p_obs)/fitted_err_p_equiv; # low = shifted right, high = shifted left
	pval <- pnorm(Z, lower.tail=TRUE)
	effect_size <- p_obs/p_equiv;
	return(list(pval = pval, fold_change = effect_size))
}

# Use the fact that S as a function of P is more stable to noise for the main part of the curve
hidden__test_DE_S_equiv <- function (expr_mat, fit=NA, method="propagate") {
	gene_info <- bg__calc_variables(expr_mat);
	if (is.na(fit[1])) {
		fit <- bg__fit_MM(gene_info$p, gene_info$s);
	}
	p_obs <- gene_info$p;
	N <- length(expr_mat[1,]);
	p_err <- gene_info$p_stderr;
	S_mean <- gene_info$s
	S_err <- gene_info$s_stderr
	K_err <- fit$Kerr;
	S_equiv <- hidden__invert_MM(fit$K,p_obs);

	## Monte Carlo method to estimate error around S_equiv ##
	MC_err <- function (p_base) {
		p_rand <- rnorm(10000, mean = p_base, sd = p_err);
		p_rand <- p_rand[p_rand > 0 & p_rand < 1]
		K_rand <- rnorm(length(p_rand),fit$K,sd = K_err);
		K_rand[K_rand < 1] <- 1;
		S_equiv_rand <- hidden__invert_MM(K_rand, p_rand)
		sd(S_equiv_rand)
	}
	if (method == "MC") {
		S_equiv_err <- unlist(lapply(p_obs,MC_err))
	} else {
		S_equiv_err <- S_equiv*sqrt(2*(p_err/p_obs)^2+(K_err/fit$K)^2);
	}

	Z <- (S_equiv - S_mean)/sqrt(S_err^2+S_equiv_err^2); # low = shifted right, high = shifted left
	pval <- pnorm(Z, lower.tail=TRUE)*2
	effect_size <- (S_mean-S_equiv)/S_equiv;
	return(list(pval = pval, effect = effect_size))
}

bg__get_extreme_residuals <- function (expr_mat, fit=NA, v_threshold=c(0.05,0.95), percent = NA, fdr_threshold = 0.1, direction="right", suppress.plot = FALSE) {
	gene_info = bg__calc_variables(expr_mat);
	if (is.na(fit)[1]) {
		fit <- bg__fit_MM(gene_info$p, gene_info$s);
	}
	res <- bg__horizontal_residuals_MM_log10(fit$K, gene_info$p, gene_info$s)
	res <- res[gene_info$p < max(v_threshold) & gene_info$p > min(v_threshold)]

	if (is.na(percent)) {
		mu <- mean(res); sigma <- sd(res);
		# deal with potential bi-modality
		if (sum(res > mu-sigma & res < mu+sigma) < 0.5) { # should be 0.68 theoretically
			mu <- mean(res[res > quantile(res,0.33)]);
			sigma <- sd(res[res > quantile(res,0.33)]);
		}
	
		if (direction=="right") {
			pval <-pnorm((res-mu)/sigma, lower.tail=FALSE)
		} else {
			pval <- pnorm((res-mu)/sigma, lower.tail=TRUE)
		}
		qval <- p.adjust(pval, method="fdr");
		sig <- qval < fdr_threshold;

		# Plot fitted normal curve
		if (!suppress.plot) {
			hist(res, col="grey75", xlab="horizontal residuals", main="", prob=TRUE)
			curve(dnorm(x,mean=mu, sd=sigma), add=TRUE);
			if (direction=="right" & sum(sig) > 0) {
				abline(v=min(res[sig]), col="red");
			} else {
				abline(v=max(res[sig]), col="red");
			}
		}
		return(names(pval)[sig]);
	} else {
		if (direction=="right") {
			cut_off <- quantile(res,prob=1-percent);
			return(names(res)[res > cut_off]);
		} else {
			cut_off <- quantile(res,prob=percent);
			return(names(res)[res < cut_off]);
		}
	}
}
##### Assembled Analysis Chunks ####
M3DropDifferentialExpression <- function(expr_mat, mt_method="bon", mt_threshold=0.05, suppress.plot=FALSE) {
	BasePlot <- bg__dropout_plot_base(expr_mat, xlim = NA, suppress.plot=suppress.plot);
	MM <- bg__fit_MM(BasePlot$gene_info$p, BasePlot$gene_info$s);
	if (!suppress.plot) {
		sizeloc <- bg__add_model_to_plot(MM, BasePlot, lty=1, lwd=2.5, col="black",legend_loc = "topright");
	}
	DEoutput <- bg__test_DE_K_equiv(BasePlot$gene_info, fit=MM);
	sig <- which(p.adjust(DEoutput$pval, method=mt_method) < mt_threshold);
	DEgenes <- rownames(expr_mat)[sig];
	DEgenes <- DEgenes[!is.na(DEgenes)];
	if (!suppress.plot) {
		bg__highlight_genes(BasePlot, expr_mat, DEgenes);
	}
	
	TABLE <- data.frame(Gene = DEgenes, p.value = DEoutput$pval[sig], q.value= p.adjust(DEoutput$pval, method=mt_method)[sig])
	return(TABLE)
}

M3DropGetExtremes <- function(expr_mat, fdr_threshold = 0.1, percent = NA, v_threshold=c(0.05,0.95), suppress.plot=FALSE) {
	BasePlot <- bg__dropout_plot_base(expr_mat, xlim = NA, suppress.plot=suppress.plot);
	MM <- bg__fit_MM(BasePlot$p, BasePlot$s);
	if (!suppress.plot) {
		sizeloc <- bg__add_model_to_plot(MM, BasePlot, lty=1, lwd=2.5, col="black",legend_loc = "topright");
	}

	if (is.na(percent)) {
		shifted_right <- bg__get_extreme_residuals(expr_mat, fit=MM, v_threshold=v_threshold, fdr_threshold = fdr_threshold, direction="right", suppress.plot=TRUE)
		shifted_left  <- bg__get_extreme_residuals(expr_mat, fit=MM, v_threshold=v_threshold, fdr_threshold = fdr_threshold, direction="left",  suppress.plot=TRUE)
	} else {
		shifted_right <- bg__get_extreme_residuals(expr_mat, fit=MM, v_threshold=v_threshold, percent = percent, direction="right", suppress.plot=TRUE)
		shifted_left  <- bg__get_extreme_residuals(expr_mat, fit=MM, v_threshold=v_threshold, percent = percent, direction="left",  suppress.plot=TRUE)

	}
	if (!suppress.plot) {
		bg__highlight_genes(BasePlot, shifted_right, colour="orange");
		bg__highlight_genes(BasePlot, shifted_left, colour="purple");
	}
	return(list(left=shifted_left,right=shifted_right));
}

M3DropTestShift <- function(expr_mat, genes_to_test, name="", background=rownames(expr_mat), suppress.plot=FALSE) {
	BasePlot <- bg__dropout_plot_base(expr_mat, xlim = NA, suppress.plot=suppress.plot);
	MM <- bg__fit_MM(BasePlot$p, BasePlot$s);
	if (!suppress.plot) {
		sizeloc <- bg__add_model_to_plot(MM, BasePlot, lty=1, lwd=2.5, col="black",legend_loc = "topright");
		bg__highlight_genes(BasePlot, genes_to_test, colour="purple");
		title(main=name);
	}

	res <- bg__horizontal_residuals_MM_log10(MM$K, BasePlot$p, BasePlot$s)
	res[is.infinite(res)] <- NA;
	mu <- median(res[rownames(expr_mat) %in% as.character(background)], na.rm=TRUE); #sigma = sd(res, na.rm=TRUE);
	s_mu <- median(res[rownames(expr_mat) %in% as.character(genes_to_test)], na.rm=TRUE);
	pval<-suppressWarnings(wilcox.test(res[rownames(expr_mat) %in% as.character(genes_to_test)], res[rownames(expr_mat) %in% as.character(background)], na.rm=TRUE)$p.value)
	return(data.frame(sample=s_mu, background=mu, p.value=pval));
}
