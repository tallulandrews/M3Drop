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

hidden_get_K <- function(expr_mat) {
	stuff <- bg__calc_variables(expr_mat);
	fit <- bg__fit_MM(stuff$p,stuff$s);
	return(fit$K);
}

bg__get_mean2disp <- function(expr_mat) {
	cv2 <- matrixStats::rowVars(expr_mat)/((rowMeans(expr_mat, na.rm=T))^2)
	xes <- log(rowMeans(expr_mat, na.rm=T))/log(10)
	reg <- lm(log(cv2[xes > 0])~xes[xes>0])
	mean2disp_fun <- function(mu){
		cv2 <- exp(reg$coeff[1]+reg$coeff[2]*(log(mu)/log(10)))
		variance <- cv2*(mu^2)
		if (variance <= mu) {variance<-1.01*mu}
		disp <- mu^2/(variance-mu)
		return(1/disp)
	}
	return(mean2disp_fun)
}

bg__fitdispersion <- function(expr_mat) {
	V <- matrixStats::rowVars(expr_mat)
	mu <- rowMeans(expr_mat)
	xes <- log(mu)/log(10)
	V[V <= mu] <- mu[V <= mu]+10^-10;
	nb_size <- mu^2/(V-mu);

	reg <- lm(log(nb_size[xes>0])~xes[xes>0])
	return(reg$coefficients[2]);
}

hidden__cv2coeffs <- function(expr_mat) {
	cv2 <- matrixStats::rowVars(expr_mat)/((rowMeans(expr_mat, na.rm=T))^2)
	xes <- log(rowMeans(expr_mat, na.rm=T))/log(10)
	reg <- lm(log(cv2[xes > 0])~xes[xes>0])
	return(c(reg$coeff[1], reg$coeff[2]))
}

hidden_calc_p <- function(obs, mu, K, disp) {
	if (mu == 0 & obs != 0) {stop("Error:non-zero obs has zero mean")}
	if (obs == 0) {
		p <- 1-mu/(mu+K)
	} else {
		p <- dnbinom(obs, size=1/disp, mu=mu)
		if (p < 10^-200) {
			p = 10^-200
		}
	}
	return(p);
}
unfinished__m3dTraditionalDE <- function(expr_mat, groups, batches=rep(1, times=length(expr_mat[1,])), fdr=0.05) {
	# Batch-specific mean-variance
	# Check Input
	if (!is.factor(batches)) {
		batches <- factor(batches)	
	}
	if (!is.factor(groups)) {
		groups <- factor(groups)	
	}
	if (length(batches) != length(groups) | 
	    length(batches) != length(expr_mat[1,])) {
		stop("Error: length of groups and batches must match number of cells (columns of expr_mat)");
	}
	if (!is.matrix(expr_mat)) {
		expr_mat <- as.matrix(expr_mat);
	}
	# Fit Batches
	batch_levels <- levels(batches)
	Ks <- sapply(batch_levels, function(b){hidden_get_K(expr_mat[,batches == b])})
	DispFits <- sapply(batch_levels, function(b){bg__get_mean2disp(expr_mat[,batches == b])})

	Ms <- rowMeans(expr_mat, na.rm=T)
	Mis <- by(t(expr_mat), groups, colMeans)

	#### Move to C? ####
	AllOut <- sapply(1:length(expr_mat[,1]), function(g) {
#	for (g in 1:length(expr_mat[,1])) {
		probs <- sapply(1:length(expr_mat[1,]), function(i) {
				obs<-expr_mat[g,i]
				M <- Ms[g]
				b <- as.numeric(batches[i])
				group <- as.numeric(groups[i])
				Mi <- Mis[[group]][g]
				p1 <- hidden_calc_p(round(obs),M,Ks[b], DispFits[[b]](M))
				p2 <- hidden_calc_p(round(obs),Mi,Ks[b], DispFits[[b]](Mi))
				return(cbind(p1,p2))
			})
		D <- -2*(sum(log(probs[1,]))-sum(log(probs[2,])))
		df <- length(levels(groups))-1
		pval <- pchisq(D, df=df, lower.tail=FALSE)
		output <- c(as.vector(by(expr_mat[g,], groups, mean)),as.vector(by(expr_mat[g,], batches, mean)), pval)
	})
#		if (g == 1) {
#			AllOut <- output
#		} else {
#			AllOut = rbind(AllOut, output)
#		}
#	}
	#### --------- ####
	AllOut <- t(AllOut)
	rownames(AllOut) <- rownames(expr_mat)
	AllOut <- cbind(AllOut, p.adjust(AllOut[,length(AllOut[1,])], method="fdr"));		
	colnames(AllOut) <- c(levels(groups), levels(batches), "p.value", "q.value")
	AllOut <- AllOut[AllOut[,length(AllOut[1,])] < fdr,]
	return(AllOut);
}


unfinished__m3dTraditionalDEShiftDisp <- function(expr_mat, groups, batches=rep(1, times=length(expr_mat[1,])), fdr=0.05) {
	# Batch specific mean-variance, gene-specific variance.
	# Check Input
	if (!is.factor(batches)) {
		batches <- factor(batches)	
	}
	if (!is.factor(groups)) {
		groups <- factor(groups)	
	}
	if (length(batches) != length(groups) | 
	    length(batches) != length(expr_mat[1,])) {
		stop("Error: length of groups and batches must match number of cells (columns of expr_mat)");
	}
	if (!is.matrix(expr_mat)) {
		expr_mat <- as.matrix(expr_mat);
	}
	# Fit Batches
	batch_levels <- levels(batches)
	Ks <- sapply(batch_levels, function(b){hidden_get_K(expr_mat[,batches == b])})
	DispFits <- sapply(batch_levels, function(b){bg__fitdispersion(expr_mat[,batches == b])}) # Fit mean-variance for each batch

	Ms <- rowMeans(expr_mat, na.rm=T)
	Mis <- by(t(expr_mat), groups, colMeans)
	V <- matrixStats::rowVars(expr_mat)
	V[V<=Ms] <- Ms[V<=Ms] + 10^-10;
	nb_size <- Ms^2/(V-Ms); # gene-specific dataset-wide dispersion

	#### Move to C? ####
	AllOut <- sapply(1:length(expr_mat[,1]), function(g) {
#	for (g in 1:length(expr_mat[,1])) {
		probs <- sapply(1:length(expr_mat[1,]), function(i) {
				obs<-expr_mat[g,i]
				M <- Ms[g]
				b <- as.numeric(batches[i])
				group <- as.numeric(groups[i])
				Mi <- Mis[[group]][g]
				# Shift Dispersion
				slope <- DispFits[b]
				disp1 <- nb_size[g]
				if (disp1 <= 0) {disp1 = 10^-10}
				tmp_intercept <- log(disp1)-slope*log(M)
				disp2 <- exp(slope*log(Mi)+tmp_intercept)
				if (disp2 <= 0) {disp2 = 10^-10}
				# Calculate probability
				p1 <- hidden_calc_p(round(obs),M,Ks[b], 1/disp1)
				p2 <- hidden_calc_p(round(obs),Mi,Ks[b], 1/disp2)
				return(cbind(p1,p2))
			})
		D <- -2*(sum(log(probs[1,]))-sum(log(probs[2,])))
		df <- length(levels(groups))-1
		pval <- pchisq(D, df=df, lower.tail=FALSE)
		output <- c(as.vector(by(expr_mat[g,], groups, mean)),as.vector(by(expr_mat[g,], batches, mean)), pval)
	})
#		if (g == 1) {
#			AllOut <- output
#		} else {
#			AllOut = rbind(AllOut, output)
#		}
#	}
	#### --------- ####
	AllOut <- t(AllOut)
	rownames(AllOut) <- rownames(expr_mat)
	AllOut <- cbind(AllOut, p.adjust(AllOut[,length(AllOut[1,])], method="fdr"));		
	colnames(AllOut) <- c(levels(groups), levels(batches), "p.value", "q.value")
	AllOut <- AllOut[AllOut[,length(AllOut[1,])] < fdr,]
	return(AllOut);
}
