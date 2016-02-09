count_zeros <-function(x) {sum(x == 0)}
num_zeros_wgt <-function(x,w) {sum(w[x==0])}

rowZero_wgt <- function(m,w) {
	w[is.na(m)] = 0;
	m[is.na(m)] = 0;
#	i = 1:length(m[,1])
#	unlist(lapply(i, function(j) {num_zeros_wgt(m[j,],w[j,])}))
	rowSums( (m==0)*w )
}

colZero_wgt <- function(m,w) {
	w[is.na(m)] = 0;
	m[is.na(m)] = 0;
#	i = 1:length(m[1,])
#	unlist(lapply(i, function(j) {num_zeros_wgt(m[,j],w[,j])}))
	colSums( (m==0)*w )
}

Calc_dropout_weights<- function(m) {
	nrow = dim(m)[1]
	ncol = dim(m)[2]
	drops_c = colSums(m==0);
	weight_c = min(drops_c)/drops_c; # down weight excess zeros
	weight_m = matrix(rep(weight_c,times=nrow),ncol=ncol,byrow=T);
	weight_m[m > 0] = 1;
	return(weight_m);
}

rowMeans_wgt <- function(m,w) {
	w[is.na(m)] = 0;
	m[is.na(m)] = 0;
#	unlist(apply(m*w,1,sum))/rowSums(w);
	return(rowMeans(m*w));
}

#covar_wgt <- function(x,wx,y,wy) {
#	w = wx*wy;
#	ux = sum(wx*x)/sum(wx)
#	uy = sum(wy*y)/sum(wy)
#	return((1/(sum(w)-1))*sum(w*(x-ux)*(y-uy)));
#}

var_wgt <- function(x,w) {
	#mu = sum(x*w)/sum(w);
	#s_sq = sum(w*(x-mu)*(x-mu))/(sum(w)-1)
	return(covar_wgt(x,w,x,w));
}

rowVar_wgt <- function(m,w) {
	w[is.na(m)] = 0;
	m[is.na(m)] = 0;
	# treat as frequency weights according to Wiki: https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Weighted_sample_variance
	i = 1:length(m[,1])
	unlist(lapply(i, function(j) {var_wgt(m[j,],w[j,])}))
}

colMeans_wgt <- function(m,w) {
	w[is.na(m)] = 0;
	m[is.na(m)] = 0;
	#unlist(apply(m*w,2,sum))/colSums(w);
	return(colMeans(m*w));
}

colVar_wgt <- function(m,w) {
	w[is.na(m)] = 0;
	m[is.na(m)] = 0;
	i = 1:length(m[1,])
	unlist(lapply(i, function(j) {var_wgt(m[,j],w[,j])}))
}

cor_pearson_matrix_wgt <- function(m,w) { 
	w[is.na(m)] = 0;
	m[is.na(m)] = 0;
	rs = matrix(nrow=length(m[,1]), ncol=length(m[,1]))
	ps = matrix(nrow=length(m[,1]), ncol=length(m[,1]))
	for (row_i in 1:length(m[,1])) {
		for (row_j in row_i:length(m[,1])) {
			weights = w[row_i,]*w[row_j,];
			r = covar_wgt(m[row_i,],w[row_i,],m[row_j,],w[row_j,])/(sqrt(var_wgt(m[row_i,],w[row_i,]))*sqrt(var_wgt(m[row_j,],w[row_j,])))
			t = r*sqrt((sum(weights)-2)/(1-r*r));

			p = pt(abs(t), sum(weights)-2, lower.tail=F)*2; # apparently df does not have to be an integer.

			rs[row_i,row_j] = r
			ps[row_i,row_j] = p
			rs[row_j,row_i] = r
			ps[row_j,row_i] = p
		}
	}
	return(list(coeff = rs, p.val = ps));
}

cor_spearman_matrix_wgt <- function(m,w) { # This doesn't work b/c weights completely change the ranking structure.
	w[is.na(m)] = 0;
	m[is.na(m)] = 0;
	rhos = matrix(nrow=length(m[,1]), ncol=length(m[,1]))
	ps = matrix(nrow=length(m[,1]), ncol=length(m[,1]))
	for (row_i in 1:length(m[,1])) {
		for (row_j in row_i:length(m[,1])) {
			ranks_i = rank(m[row_i,])
			ranks_j = rank(m[row_j,])
			weights = w[row_i,]*w[row_j,];
			rho = 1-(6*sum(weights*(ranks_i-ranks_j)^2))/(sum(weights)*(sum(weights)^2-1));
			z = sqrt((sum(weights)-3)/1.06)*0.5*log((1+rho)/(1-rho))
			p = pnorm(abs(z), lower.tail=F);
			rhos[row_i,row_j] = rho
			ps[row_i,row_j] = p
			rhos[row_j,row_i] = rho
			ps[row_j,row_i] = p
		}
	}
	return(list(coeff = rhos, p.val = ps));
}

kruskal_test_wgt <- function (x,w,g) {
	ranks = rank(x);
	r_bar = sum(ranks*w)/sum(w);
	denom = sum(w*(ranks-r_bar)^2);
	numer = 0;
	for (i in unique(g)) {
		ni = sum(w[g==i]);
		ri = sum(ranks[g==i]*w[g==i])/sum(w[g==i])
		numer = numer+ni*(ri-r_bar)^2
	}
	K = (sum(w)-1)*numer/denom
	pval = pchisq(K,length(unique(g))-1, lower.tail=F)
	return(list(estimate=K,p.value=pval))
}

demo <- function() {
	data_g1 = read.table("/nfs/users/nfs_t/ta6/Data/E-MTAB-2805_G1_singlecells_counts.txt", header=T)
	data_g2M = read.table("/nfs/users/nfs_t/ta6/Data/E-MTAB-2805_G2M_singlecells_counts.txt", header=T)
	data_S = read.table("/nfs/users/nfs_t/ta6/Data/E-MTAB-2805_S_singlecells_counts.txt", header=T)
	gene_data = data_g1[,1:4]

	merged = merge(data_g1,data_g2M,by=c("EnsemblGeneID","EnsemblTranscriptID","AssociatedGeneName","GeneLength"))
	merged = merge(merged,data_S,by=c("EnsemblGeneID","EnsemblTranscriptID","AssociatedGeneName","GeneLength"))
	groups = c(rep("g1", times=length(data_g1[1,])-4),rep("g2M", times=length(data_g2M[1,])-4), rep("s", times=length(data_S[1,])-4))
	clean <- function(x){
	        rownames(x) = x[,1]
	        x = x[,5:length(x[1,])];
	        filter <- apply(x, 1, function(y) {length(y[y>5])>=2})
	        x = x[filter,]
	        x = x[grep("^E",rownames(x)),]
	}
	merged = clean(merged)

	p = apply(merged, 2, function(x) {sum(x == 0)/length(x)})
	reads = colSums(merged);
	CPM = t(t(merged)/(reads/1000000));
	CPMnorm = t(t(CPM)*p/median(p));

	weights = Calc_dropout_weights(CPMnorm)
	plot(rowMeans(merged),sqrt(apply(merged,1,var))/rowMeans(merged), log="xy", xlab="Mean", ylab="Var/Mean")
	plot(rowMeans(CPMnorm),sqrt(apply(CPMnorm,1,var))/rowMeans(CPMnorm),log="xy", xlab="Mean", ylab="Var/Mean")
	plot(rowMeans_wgt(CPMnorm,weights),sqrt(rowVar_wgt(CPMnorm,weights))/rowMeans_wgt(CPMnorm,weights),log="xy",xlab="Mean", ylab="Var/Mean")

	p1 = vector(length=length(merged[,1]))
	p2 = vector(length=length(merged[,1]))
	for (i in 1:length(merged[,1])) {
		p1[i] = kruskal.test(unlist(CPMnorm[i,]),as.factor(groups))$p.value
		p2[i] = kruskal_test_wgt(CPMnorm[i,],weights[i,],groups)
	}
	thing = p.adjust(p2, method="fdr") < 0.0000000000005 & p.adjust(p1,method="fdr") < 0.0000000000005
	heatmap(log(as.matrix(CPMnorm[thing,]+1)), ColSideColors=colours[as.factor(groups)], Colv=NA)
}




