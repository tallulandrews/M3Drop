bg__calc_variables <- function(expr_mat) {
        # Calc variables
	if (sum(expr_mat < 0) >0) {stop("Expression matrix contains negative values! M3Drop requires an expression matrix that is not log-transformed.")}
	p = apply(expr_mat,1,function(x){y = x[!is.na(x)]; sum(y==0)/length(y)});
	s = rowMeans(expr_mat, na.rm=T);
	s_stderr = unlist(apply(expr_mat,1,sd))/sqrt(length(expr_mat[1,]));
	tmp = expr_mat; tmp[tmp == 0] = NA
	s_stderr_nozero = unlist(apply(tmp,1,sd, na.rm=T))/sqrt(rowSums(expr_mat>0));
	p_stderr = sqrt(p*(1-p)/length(expr_mat[1,]));
	names(s) = rownames(expr_mat);
	names(p) = rownames(expr_mat);
	return(list(s = s, p = p, s_stderr = s_stderr, s_stderr_nozero = s_stderr_nozero, p_stderr = p_stderr))
}

hidden__invert_MM <- function (K, p) {K*(1-p)/(p)}
bg__horizontal_residuals_MM_log10 <- function (K, p, s) {log(s)/log(10) - log(hidden__invert_MM(K,p))/log(10)}

hidden_getAUC <- function(gene, labels) {
        require("ROCR")
        ranked=rank(gene);
        ms = aggregate(ranked~unlist(labels),FUN=mean); #Get average score for each cluster
        posgroup = as.character(unlist(ms[which(ms[,2]==max(ms[,2])),1])); #Get cluster with highest average score
        if (length(posgroup) > 1) {return (c(-1,-1,-1))} # Return negatives if there is a tie for cluster with highest average score (by definition this is not cluster specific)

        # Create 1/0 vector of truths for predictions, cluster with highest average score vs everything else
        truth = labels == posgroup

        #Make predictions & get auc using RCOR package.
        pred=prediction(ranked,as.numeric(truth))
        val = unlist(performance(pred,"auc")@y.values)
        pval = wilcox.test(gene[truth],gene[!truth])$p.value
        if (!exists("pval")) {pval=NA}

        return(c(val,posgroup,pval))
}

M3Drop_getmarkers <- function(expr_mat, labels) {
	if (length(labels) != length(expr_mat[1,])) {
		stop("Length of labels does not match number of cells.")
	}
        aucs = apply(expr_mat,1,hidden_getAUC,labels=labels)
        auc_df <- data.frame(matrix(unlist(aucs), ncol=3, byrow=T))
        rownames(auc_df) = rownames(expr_mat)
        colnames(auc_df) = c("AUC","Group", "pval")
        auc_df[,1] = as.numeric(as.character(auc_df[,1]))
        auc_df[,3] = as.numeric(as.character(auc_df[,3]))
        auc_df = auc_df[auc_df[,1] > 0,]
        return(auc_df);

}

