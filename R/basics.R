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

bg__calc_variables <- function(expr_mat) {
    # Calc variables
    sum_neg <- sum(expr_mat < 0)
    sum_zero <- sum(expr_mat == 0)
    sum_pos <- sum(expr_mat >= 0)
    if (sum_neg > 0) {stop("Expression matrix contains negative values! M3Drop requires an expression matrix that is not log-transformed.")}
    
    # Deal with strangely normalized data
    if (sum_zero == 0) {
        warning("Warning: No zero values (dropouts) detected will use minimum expression value instead.")
        # If no zeros in expression matrix convert minimum value into zero
        expr_mat <- round(expr_mat, digits=2) # Round to accomodate errors
        min_val <- min(expr_mat)
        expr_mat[expr_mat == min_val] <- 0;
    }
    if (sum_zero < 0.1*sum_pos) {
        # Less than 10% zeros
        warning("Warning: Expression matrix contains few zero values (dropouts) this may lead to poor performance.")
    }
    
    p <- apply(expr_mat,1,function(x){sum(x==0)/length(x)})
    p_stderr <- sqrt(p*(1-p)/ncol(expr_mat))
    s <- rowMeans(expr_mat)
    s_stderr <- unlist(apply(expr_mat,1,sd))/sqrt(ncol(expr_mat))
    return(list(s = s, p = p, s_stderr = s_stderr, p_stderr = p_stderr))
}

hidden__invert_MM <- function (K, p) {K*(1-p)/(p)}
bg__horizontal_residuals_MM_log10 <- function (K, p, s) {log(s)/log(10) - log(hidden__invert_MM(K,p))/log(10)}

hidden_getAUC <- function(gene, labels) {
        ranked <- rank(gene);
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

M3DropGetMarkers <- function(expr_mat, labels) {
	if (length(labels) != length(expr_mat[1,])) {
		stop("Length of labels does not match number of cells.")
	}
        aucs <- apply(expr_mat,1,hidden_getAUC,labels=labels)
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
