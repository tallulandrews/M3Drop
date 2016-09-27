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

hidden_get_mean2disp <- function(expr_mat) {
	cv2 <- hidden_rowVars(expr_mat)/((rowMeans(expr_mat, na.rm=T))^2)
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

hidden_calc_p <- function(obs, mu, K, mean2disp) {
	if (mu == 0 & obs != 0) {stop("Error:non-zero obs has zero mean")}
	if (obs == 0) {
		p <- 1-mu/(mu+K)
	} else {
		disp <- mean2disp(mu)
		p <- dnbinom(obs, size=1/disp, mu=mu)
		if (p < 10^-200) {
			p = 10^-200
		}
	}
	return(p);
}
M3DropTraditionalDE <- function(expr_mat, groups, batches=rep(1, times=length(expr_mat[1,]))) {
	# Check Input
	if (!is.factor(batches)) {
		batches <- factor(batches)	
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
	Ks <- vector(length=length(batch_levels))
	DispFun <- list(length=length(batch_levels))

	for (b in 1:length(batch_levels)) {
		Ks[b] <- hidden_get_K(expr_mat[,batches == batch_levels[b]]);
		DispFun[[b]] <- hidden_get_mean2disp(expr_mat[,batches == batch_levels[b]]);
	}

	Ms <- rowMeans(expr_mat, na.rm=T)
	Mis <- by(t(expr_mat), factor(groups), colMeans)
	for (g in 1:length(expr_mat[,1])) {
		probs <- sapply(1:length(expr_mat[1,]), function(i) {
				obs<-expr_mat[g,i]
				M <- Ms[g]
				b <- as.numeric(batches[i])
				group <- as.numeric(factor(groups)[i])
				Mi <- Mis[[group]][g]
				p1 <- hidden_calc_p(round(obs),M,Ks[b], DispFun[[b]])
				p2 <- hidden_calc_p(round(obs),Mi,Ks[b], DispFun[[b]])
				return(cbind(p1,p2))
			})
		D <- -2*(sum(log(probs[1,]))-sum(log(probs[2,])))
		df <- length(unique(groups))-1
		pval <- pchisq(D, df=df, lower.tail=FALSE)
		output <- c(as.vector(by(expr_mat[g,], factor(groups), mean)),as.vector(by(expr_mat[g,], factor(batches), mean)), pval)
		if (g == 1) {
			AllOut <- output
		} else {
			AllOut = rbind(AllOut, output)
		}
	}
	rownames(AllOut) <- rownames(expr_mat)
	AllOut <- cbind(AllOut, p.adjust(AllOut[,length(AllOut[1,])], method="fdr"));		
	colnames(AllOut) <- c(levels(factor(groups)), levels(batches), "p.value", "q.value")
	return(AllOut);
}
