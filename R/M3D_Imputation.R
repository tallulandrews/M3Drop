hidden_m3dropImputeZeros <- function(expr_mat) {
	BasePlot <- bg__dropout_plot_base(expr_mat, xlim = NA, suppress.plot=TRUE)
	MM <- bg__fit_MM(BasePlot$gene_info$p, BasePlot$gene_info$s);

	no_zero_p <- BasePlot$gene_info$p
	no_zero_p[no_zero_p==0] <- 1/(length(expr_mat[1,])*2);
	Expected_S <- MM$K*(1/no_zero_p -1);
#	Expected_S <- MM$K*(1/BasePlot$gene_info$p -1);
#	Expected_S[BasePlot$gene_info$p==0] = BasePlot$gene_info$s[BasePlot$gene_info$p==0];
#	new_mat <- (expr_mat+Expected_S)
	n_cells <- ncol(expr_mat)
	expect_mat <- matrix(rep(Expected_S, times=n_cells), ncol=n_cells, byrow=F)
	imp_mat <- (expect_mat+expr_mat)/2
	new_mat <- expr_mat
	new_mat[expr_mat < expect_mat] <- imp_mat[expr_mat < expect_mat]
	
	return(new_mat)
}
