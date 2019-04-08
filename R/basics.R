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

M3DropConvertData <- function(input, is.log=FALSE, is.counts=FALSE, pseudocount=1) {
	remove_undetected_genes <- function(mat) {
		no_detect <- Matrix::rowSums(mat > 0, na.rm=T) == 0;
		print(paste("Removing ",sum(no_detect), "undetected genes."))
		return(mat[!no_detect,])
	}
	type <- class(input)[1]
	lognorm <- NULL
	counts <- NULL
	if (type == "SCESet") {
		# Old scater
		lognorm <- scater::exprs(input)
		counts <- counts(input)

	} else if (type == "SingleCellExperiment") {
		# New scater
		c <- which(names(input@assays) == "counts")
		ln <- which(names(input@assays) == "logcounts")
		norm <- which(names(input@assays) == "normcounts")
		if (length(norm) > 0) {
			return(remove_undetected_genes(input));
		} else if (length(ln) > 0) {
			lognorm <- input@assays[[ln]]
		} else if (length(c) > 0) {
			counts <- input@assays[[c]]
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
			return(remove_undetected_genes(input));
		}
	} else {
		stop(paste("Error: Unrecognized input format :", type))
	}

	# Prefer log-norm to raw counts
	if (!is.null(dim(lognorm))) {
		norm <- 2^lognorm-pseudocount
		norm <- remove_undetected_genes(norm);
		return(norm)
	} 

	# CPM transform raw counts
	if (!is.null(dim(counts))) {
		sf <- Matrix::colSums(counts)
		norm <- t( t(counts)/sf * median(sf) )
		norm <- remove_undetected_genes(norm)
		return( norm )
	}
}

bg__calc_variables <- function(expr_mat) {
    if (class(expr_mat) != "matrix" & class(expr_mat) != "dgCMatrix" & class(expr_mat) != "Matrix") {
	warning("Warning: not a recognized matrix class, coercing to 'matrix'.")
	expr_mat <- as.matrix(expr_mat)
    }
    if (sum(is.na(expr_mat)) > 0) {
	stop("Error: Expression matrix contain NA values.");
    }
    # Calc variables
#    sum_neg <- sum(expr_mat < 0)
     sum_zero <- prod(dim(expr_mat)) - sum(expr_mat > 0)
#    sum_pos <- sum(expr_mat >= 0)
#    if (is.log) {
#	expr_mat <- is.log^expr_mat-1;
#    }

    lowest <- min(expr_mat)
    if (lowest < 0) {stop("Error: Expression matrix cannot contains negative values! Has the matrix been log-transformed?")}
    
    # Deal with strangely normalized data
    if (lowest > 0) {
        warning("Warning: No zero values (dropouts) detected will use minimum expression value instead.")
        # If no zeros in expression matrix convert minimum value into zero
        #expr_mat <- round(expr_mat, digits=2) # Round to accomodate errors
        min_val <- lowest+0.05 # Instead of rounding to accomodate errors
        expr_mat[expr_mat == min_val] <- 0;
    }
    if (sum_zero < 0.1*prod(dim(expr_mat))) {
        # Less than 10% zeros
        warning("Warning: Expression matrix contains few zero values (dropouts) this may lead to poor performance.")
    }
    
    p <- 1 - Matrix::rowSums(expr_mat > 0)/ncol(expr_mat)
    if (sum(p == 1) > 0) {
	warning(paste("Warning: Removing", sum(p==1),"undetected genes."))
	expr_mat <- expr_mat[p < 1,]
        p <- 1-Matrix::rowSums(expr_mat > 0)/ncol(expr_mat)
    }

    p_stderr <- sqrt(p*(1-p)/ncol(expr_mat))
    s <- rowMeans(expr_mat)
    s_stderr <- sqrt( (rowMeans(expr_mat^2) - s^2)/ncol(expr_mat) ) #sparse matrix friendly
#    s_stderr <- rowSds(expr_mat)/sqrt(ncol(expr_mat))
    names(s_stderr) <- rownames(expr_mat)
    return(list(s = s, p = p, s_stderr = s_stderr, p_stderr = p_stderr))
}

hidden__invert_MM <- function (K, p) {K*(1-p)/(p)}
bg__horizontal_residuals_MM_log10 <- function (K, p, s) {log(s)/log(10) - log(hidden__invert_MM(K,p))/log(10)}

hidden_getAUC <- function(gene, labels) {
	labels <- unlist(labels);
        ranked <- rank(gene);
	#MAT <- as.matrix(ranked)
        #x <- split(seq(ncol(MAT)), labels)
        #ms <- sapply(x, function(a) rowMeans(MAT[,a]))

        ms <- aggregate(ranked~unlist(labels),FUN=mean); #Get average score for each cluster
        posgroup <- as.character(unlist(ms[which(ms[,2]==max(ms[,2])),1])); #Get cluster with highest average score
        if (length(posgroup) > 1) {return (c(-1,-1,-1))} # Return negatives if there is a tie for cluster with highest average score (by definition this is not cluster specific)

        # Create 1/0 vector of truths for predictions, cluster with highest average score vs everything else
        truth <- labels == posgroup

        #Make predictions & get auc using RCOR package.
        pred <- ROCR::prediction(ranked,as.numeric(truth))
        val <- unlist(ROCR::performance(pred,"auc")@y.values)
        pval <- wilcox.test(gene[truth],gene[!truth])$p.value
        if (!exists("pval")) {pval <- 1}

        return(c(val,posgroup,pval))
}

hidden_fast_AUC_m3drop <- function(expression_vec, labels) {
	R = rank(expression_vec);
	labels <- unlist(labels);
	ms <- aggregate(R~unlist(labels),FUN=mean); #Get average score for each cluster
	posgroup <- as.character(unlist(ms[which(ms[,2]==max(ms[,2])),1])); #Get cluster with highest average score
	if (length(posgroup) > 1) {return (c(-1,-1,-1))} # Return negatives if there is a tie for cluster with highest average score (by definition this is not cluster specific)
	truth <- labels == posgroup

	out <- wilcox.test(expression_vec[truth], expression_vec[!truth], alternative="two.sided")
	
        N1 = sum(truth)
        N2 = sum(!truth);
        #U1 = sum(R[truth])-N1*(N1+1)/2
        U2 = sum(R[!truth])-N2*(N2+1)/2
        if (N1 == 0) {return(c(0,0,0))}
        if (N2 == 0) {return(c(1,1,1))}
        AUC = 1-U2/(N1*N2);  # smaller valued ranks (i.e. lower expression values for !true).
        # assumes large sample size
        #  https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/PASS/Confidence_Intervals_for_the_Area_Under_an_ROC_Curve.pdf
        # originally (Hanley and McNeil 1982)
	
        return(c(AUC, posgroup, out$p.value));
}

M3DropGetMarkers <- function(expr_mat, labels) {
	if (length(labels) != length(expr_mat[1,])) {
		stop("Length of labels does not match number of cells.")
	}
        aucs <- apply(expr_mat,1,hidden_fast_AUC_m3drop,labels=labels)
        auc_df <- data.frame(matrix(unlist(aucs), ncol=3, byrow=TRUE))
        rownames(auc_df) <- rownames(expr_mat)
        colnames(auc_df) <- c("AUC","Group", "pval")
	auc_df$Group <- as.character(auc_df$Group)
	if(sum(auc_df$Group == "-1") > 0) {
		auc_df$Group[auc_df$Group == "-1"] <- "Ambiguous";
	}
        auc_df[,1] <- as.numeric(as.character(auc_df[,1]))
        auc_df[,3] <- as.numeric(as.character(auc_df[,3]))
        auc_df <- auc_df[auc_df[,1] > 0,]
	auc_df <- auc_df[order(-auc_df$AUC),]
        return(auc_df);
}
