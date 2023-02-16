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

# Normalization Functions
hidden__UQ <- function(x){quantile(x[x>0],0.75)};

bg__filter_cells <- function(expr_mat,labels=NA, suppress.plot=FALSE, min_detected_genes=NA) {
	num_detected <-  Matrix::colSums(expr_mat > 0, na.rm=TRUE);
	if (!is.na(min_detected_genes)) {
		low_quality <- num_detected < min_detected_genes;
	} else {
		num_zero <- Matrix::colSums(expr_mat == 0, na.rm=TRUE);
		cell_zero <- num_zero;
		mu <- mean(cell_zero);
		sigma <- sd(cell_zero);
		# Deal with bi-modal
		if (sum(cell_zero > mu-sigma & cell_zero < mu+sigma) < 0.5) { # should be 0.68 theoretically
			mu <- mean(cell_zero[cell_zero < median(cell_zero)]);
			sigma <- sd(cell_zero[cell_zero < median(cell_zero)]);
		}
		low_quality <- p.adjust(pnorm((cell_zero-mu)/sigma, lower.tail=FALSE), method="fdr") < 0.05;
		if (!suppress.plot) {
			hist(cell_zero, col="grey75", xlab="Number of zeros (per cell)", main="", prob=TRUE)
			curve(dnorm(x,mean=mu, sd=sigma), add=TRUE)
			if (sum(low_quality) > 0) {
				abline(v=min(cell_zero[low_quality]), col="red")
			}
		}
	}
	if (sum(low_quality) > 0) {
		if (length(labels)==length(expr_mat[1,])) {labels = labels[!low_quality]}
		expr_mat <- expr_mat[,!low_quality];
	}
	return(list(data = expr_mat, labels = labels));
}

hidden__normalize <- function(data) {
	# Combine UQ and detection rate adjusted normalization 
	# Stephanie Hick, Mingziang Teng, Rafael A Irizarry "On the widespread and critical impact of systematic single-cell RNA-Seq data" http://dx.doi.org/10.1101/025528 
	cell_zero <- Matrix::colSums(data == 0)/length(data[,1]);
	uq <- unlist(apply(data,2,hidden__UQ));
	normfactor <- (uq/median(uq)) * (median(cell_zero)/cell_zero); 
	data <- t(t(data)/normfactor);
	return(data);
}

M3DropCleanData <- function(expr_mat, labels = NA, is.counts=TRUE, suppress.plot=FALSE, pseudo_genes=NA, min_detected_genes=NA) {
	expr_mat[is.na(expr_mat)] = 0;
	if (length(pseudo_genes) > 1) {
		is_pseudo <- rownames(expr_mat) %in% as.character(pseudo_genes);
	        expr_mat <- expr_mat[!is_pseudo,];
	}

	data_list <- bg__filter_cells(expr_mat, labels, suppress.plot = suppress.plot, min_detected_genes=min_detected_genes);
	
        detected <- Matrix::rowSums(data_list$data > 0) > 3;
        expr_mat <- data_list$data[detected,];
        detected <- Matrix::rowSums(data_list$data> 0) > 3;
        expr_mat <- data_list$data[detected,];
	labels   <- data_list$labels

	spikes <- grep("ercc",rownames(expr_mat), ignore.case=TRUE)
	if (is.counts) {
		if (length(spikes) > 1) {
	                totreads <- Matrix::colSums(expr_mat[-c(spikes),])
		} else {
			totreads <- Matrix::colSums(expr_mat);
		}
                cpm <- t(t(expr_mat)/totreads)*1000000;
                lowExpr <- rowMeans(cpm) < 10^-5;
                cpm<-cpm[!lowExpr,];
                return(list(data=cpm, labels=labels));
        }

	lowExpr <- rowMeans(expr_mat) < 10^-5;
        data<-expr_mat[!lowExpr,];
        return(list(data=expr_mat, labels=labels));
}

#### Pearson Residuals ####
# February 16, 2023

NBumiPearsonResiduals <- function(counts, fits=NULL) {
	if (is.null(fits)) {
		fits <- NBumiFitModel(counts)
	}
	mus <- t(t(fits$vals$tjs/fits$vals$total)) %*% fits$vals$tis
	pearson <- (counts-mus)/sqrt(mus + mus^2/fits$sizes)
	return(pearson)
}

NBumiPearsonResidualsApprox <- function(counts, fits=NULL) {
	if (is.null(fits)) {
		vals <- hidden_calc_vals(counts);
	} else {
		vals <- fits$vals
	}
	mus <- t(t(vals$tjs/vals$total)) %*% vals$tis
	pearson <- (counts-mus)/sqrt(mus)
	return(pearson)
}
#############################
