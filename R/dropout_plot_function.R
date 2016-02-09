# Modularize this stuff more sensibly
#  Plotting Functions
bg__dropout_plot_base <- function (expr_mat, xlim = NA, suppress.plot=FALSE) {
	require("RColorBrewer")
	
	gene_info = bg__calc_variables(expr_mat);

        xes = log(gene_info$s)/log(10);
        put_in_order = order(xes);
        fancy <- densCols(xes, gene_info$p, colramp=colorRampPalette(c("black","white")))
        dens <- col2rgb(fancy)[1,]+1L
#        colours <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F",
#                                    "#FCFF00", "#FF9400", "#FF3100"))(256) #rainbow
        colours <-  colorRampPalette(c("#000099", "#00FEFF", "#FCFF00"))(256) #blue->yellow
        dens.col = colours[dens]

        par(fg="black")
	if (!suppress.plot) {
		if (!(sum(is.na(xlim)))) {
	        	plot(xes,gene_info$p, main="", ylab="Dropout Proportion", xlab="log(expression)", col = dens.col,pch=16, xlim=xlim, ylim=c(0,1))
		} else {
	        	plot(xes,gene_info$p, main="", ylab="Dropout Proportion", xlab="log(expression)", col = dens.col,pch=16)
		}
	}
	invisible(list(P=gene_info$p, S=gene_info$s, xes=xes, data=expr_mat, order=put_in_order));
}

bg__add_model_to_plot <- function(fitted_model, base_plot, lty=1, lwd=1, col="black",legend_loc = "topright") {
	lines(base_plot$xes[base_plot$order],fitted_model$predictions[base_plot$order],lty=lty,lwd=lwd,col=col);
	par(fg=col)
	if (length(legend_loc) == 2) {
        	this_loc = legend(legend_loc[1], legend_loc[2], fitted_model$model, box.lty=lty, box.lwd=lwd, xjust=1)
	} else {
		this_loc = legend(legend_loc[1], fitted_model$model, box.lty=lty, box.lwd=lwd, xjust=1)
	}
	par(fg="black")
	invisible(this_loc)
}

bg__highlight_genes <- function (base_plot, genes, colour="purple", pch=16) {
	if(!is.numeric(genes) && !is.logical(genes)) {
		genes = match(as.character(genes), rownames(base_plot$data));
		nomatch = sum(is.na(genes));
		if (nomatch > 0) {warning(paste(nomatch, " genes could not be matched to data, they will not be highlighted."));}
		genes = genes[!is.na(genes)];
	}
	points(base_plot$xes[genes],base_plot$P[genes],col=colour, pch=pch)
}

bg__expression_heatmap <- function (genes, data, cell_labels=NA, gene_labels=NA, key_genes=NA, key_cells=NA) { 
	require("RColorBrewer")
	require("gplots")
	if(!is.numeric(genes)) {
		genes = match(genes, rownames(data));
		nomatch = sum(is.na(genes));
		if (nomatch > 0) {warning(paste(nomatch, " genes could not be matched to data, they will not be included in the heatmap."));}
		genes = genes[!is.na(genes)];
	}
	if (length(genes) < 1) {warning("No genes for heatmap.");return();}
	# Plot heatmap of expression
	heatcolours <- rev(brewer.pal(11,"RdBu"))
	col_breaks = c(-100,seq(-2,2,length=10),100)
	heat_data = as.matrix(data[genes,])
	heat_data = log(heat_data+1)/log(2);
	ColColors = rep("white", times=length(heat_data[1,]))
	RowColors = rep("white", times=length(heat_data[,1]))
	# remove row & column labels
	rownames(heat_data) = rep("", length(heat_data[,1]));
	if (!is.na(key_genes)) {
		rownames(heat_data)[rownames(data[genes,]) %in% key_genes] = rownames(data[genes,])[rownames(data[genes,]) %in% key_genes]; 
	}
	colnames(heat_data) = rep("", length(heat_data[1,]));
	if (!is.na(key_cells)) {
		colnames(heat_data)[colnames(data[genes,]) %in% key_cells] = colnames(data[genes,])[colnames(data[genes,]) %in% key_cells]; 
	}
	if (!is.na(cell_labels[1])) {
		colours = as.factor(cell_labels)
		palette = brewer.pal(max(3,length(unique(cell_labels))), "Set3");
		ColColors = palette[colours];	
		mylegend<- list(names = unique(cell_labels), fill = unique(ColColors));
	} 
	if (!is.na(gene_labels[1])) {
		# lowest factor level = grey (so 0-1 is striking)
		if (!is.numeric(gene_labels)) {
			colours = as.factor(gene_labels)
		} else {
			colours = gene_labels
		}
		palette = c("grey75",brewer.pal(max(3,length(unique(gene_labels))), "Set1"));
		RowColors = palette[colours];
	}
	# Custom Shit
	lwid=c(1,0.2,4)
	lhei=c(1,0.2,4)
	lmat=rbind(c(6,0,5),c(0,0,2),c(4,1,3))


	heatmap_output = heatmap.2(heat_data, ColSideColors = ColColors, RowSideColors = RowColors, col=heatcolours, breaks=col_breaks, scale="row",symbreaks=T, trace="none", dendrogram="column", key=FALSE, Rowv=TRUE, Colv=TRUE,lwid=lwid, lhei=lhei,lmat=lmat, hclustfun=function(x){hclust(x,method="ward.D2")})
	# Custom key
	par(fig = c(0, 1/(5.2),4/(5.2), 1), mar=c(4,1,1,1), new=TRUE)
	scale01 <- function(x, low = min(x), high = max(x)) {
        	x <- (x - low)/(high - low)
        	x
    	}
	par(mar=c(5,1,1,1))
	par(cex=0.75)
	par(mgp=c(2,1,0))
	key_breaks = seq(-2,2,length=10)
	key_col = heatcolours[2:(length(heatcolours)-1)]
	z = seq(min(key_breaks),max(key_breaks), by=min(diff(key_breaks)/4))
	image(z=matrix(z,ncol=1),col=key_col,breaks=key_breaks,xaxt="n",yaxt="n")
	par(usr = c(0, 1, 0, 1))
	lv <- pretty(key_breaks)
        xv <- scale01(as.numeric(lv), min(key_breaks),max(key_breaks))
        xargs <- list(at = xv, labels = lv)
	xargs$side <- 1
	do.call(axis, xargs)
	mtext(side = 1, "Expression Z-Score", line = par("mgp")[1], padj = 0.5, 
                cex = par("cex") * par("cex.lab"))

	# Legend
	par(fig = c(0/5.2, 1/(5.2),0/(5.2), 4/5.2), mar=c(0,0,0,0), new=TRUE)
	par(mar=c(0,0,0,0))
	if (!is.na(cell_labels[1])) {
		legend("left", mylegend$names, pt.bg = mylegend$fill,bg="white",col="black", pch=22, pt.cex=2.5, cex=1.25, bty="n",y.intersp = 2);
	}
	invisible(heatmap_output);
}

# Model-fitting/manipulation Functions
bg__calc_variables <- function(expr_mat) {
        # Calc variables
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

bg__invert_MM <- function (K, p) {K*(1-p)/(p)}
bg__horizontal_residuals_MM_log10 <- function (K, p, s) {log(s)/log(10) - log(bg__invert_MM(K,p))/log(10)}
#bg__num.zero <- function(x){sum(x==0)}
bg__fit_MM <- function (p,s) {
#        fit = nls(p ~ 1-(s/((krt+s))),data.frame(s=s),start=list(krt=3))
#	K_glm = glm(p ~ offset(-1*log(s)), family="binomial")
#	Kerr = summary(K_glm)$coeff[1,2];
#	Kcoeff = summary(K_glm)$coeff[1,1];
#	Kerr = exp(Kcoeff+Kerr)-exp(Kcoeff)
#        predicted = fitted(fit)
#        krt=summary(fit)$parameters[1,1]
#	return(list(K=krt,Kerr=Kerr,predictions=predicted, model=c("MMenton",paste("Krt =",round(krt,digits=3))),SSr=round(sum((residuals(fit))^2)),SAr=round(sum(abs(residuals(fit))))))
	if (length(p) != length(s)) {
		stop(print("Error: p and s not same length. Cannot fit Michaelis-Menten."))
	}
	require("bbmle")
	LL <- function(krt,sigma) {
		R = log(p)/log(10)-log((1-(s/((krt+s)))))/log(10) #log normal error
		R = suppressWarnings(dnorm(R,0,sigma,log=TRUE))
		-sum(R)
	}
	fit = mle2(LL,start=list(krt=3, sigma=0.25))
	thing = summary(fit)
	krt = fit@coef[1]
	res_err = attributes(summary(fit))$coef[2,1]
	Kerr = fit@coef[2]
	predicted = 1-(s/(krt+s))
	residuals = p-predicted
	return(list(K=krt,Kerr=Kerr,fitted_err = res_err,predictions=predicted, model=c("MMenton",paste("Krt =",round(krt,digits=3))),SSr=round(sum((residuals)^2)),SAr=round(sum(abs(residuals)))))

}
bg__fit_logistic <- function(p,s) {
	if (length(p) != length(s)) {
		stop(print("Error: p and s not same length. Cannot fit Logistic Regression."))
	}
        logistic = glm(p~log(s),family="binomial")
        predlog = fitted(logistic)
	return(list(predictions=predlog, B0 = logistic$coeff[1], B1=logistic$coeff[2] ,model=c( "Logistic", paste("Intercept =",round(logistic$coeff[1],digits=3)),paste("Coeff =",round(logistic$coeff[2],digits=3))),SSr=round(sum((fitted(logistic)-p)^2)),SAr=round(sum(abs(fitted(logistic)-p)))));
#	require("bbmle")
#	LL <- function(B0,B1,sigma) {
#		R = p-(1/(1+exp(-B0+B1*log(s)/log(10))))
#		R = suppressWarnings(dnorm(R,0,sigma,log=TRUE))
#		-sum(R)
#	}
#	fit = mle2(LL,start=list(B0=2, B1=-1, sigma=0.25))
#	thing = summary(fit)
#	B0 = attributes(summary(fit))$coef[1,1]
#	B1 = attributes(summary(fit))$coef[2,1]
#	res_err = attributes(summary(fit))$coef[3,1]
#	B0err = attributes(summary(fit))$coef[1,2]
#	B1err = attributes(summary(fit))$coef[2,2]
#	predicted = (1/(1+exp(-B0+B1*log(s)/log(10))))
#	residuals = p-predicted
#	return(list(B0=B0,B0err=B0err,B1=B1,B1err=B1err,fitted_err = res_err,predictions=predicted, model=c("Logistic",paste("Intercept =",round(B0,digits=3)),paste("B1 =",round(B1,digits=3))),SSr=round(sum((residuals)^2)),SAr=round(sum(abs(residuals)))))
}

bg__fit_ZIFA <- function(p,s) {
	if (length(p) != length(s)) {
		stop(print("Error: p and s not same length. Cannot fit double exponential."))
	}
#	doubleXfit = nls(p ~ exp(-lambda*s*s),data.frame(s=s),start=list(lambda=0.01), control=list(maxiter=100), algorithm="port", lower=list(lambda=0));
#	preddoubleX = fitted(doubleXfit);
#	lambda=summary(doubleXfit)$parameters[1,1];
#	return(list(predictions=preddoubleX, lambda=lambda, model=c("p ~ e^(-lambda*S^2)",paste("lambda =",signif(lambda,digits=2))),SSr = round(sum((residuals(doubleXfit))^2)),SAr = round(sum(abs(residuals(doubleXfit))))));
	require("bbmle")
	LL <- function(lambda,sigma) {
		R = p-exp(-lambda*s*s)
		R = suppressWarnings(dnorm(R,0,sigma,log=TRUE))
		-sum(R)
	}
	fit = mle2(LL,start=list(lambda=0.01, sigma=0.25))
	thing = summary(fit)
	lambda = attributes(summary(fit))$coef[1,1]
	res_err = attributes(summary(fit))$coef[2,1]
	Lerr = attributes(summary(fit))$coef[1,2]
	predicted = exp(-lambda*s*s)
	residuals = p-predicted
	return(list(lambda=lambda,Lerr=Lerr,fitted_err = res_err,predictions=predicted, model=c("p ~ e^(-lambda*S^2)",paste("lambda =",round(lambda,digits=2))),SSr=round(sum((residuals)^2)),SAr=round(sum(abs(residuals)))))
}

# Normalization Functions
#bg__UQ <- function(x){quantile(x[x>0],0.75)};
#bg__filter_genes <- function(data) {
#        # get rid of genes with 0 expression
##        filter <- apply(data, 1, function(x) length(x[x>5])>=2);
#	filter = rowSums(data > 5) >=2;
#        data = data[filter,];
#	return(data);
#}

bg__filter_cells <- function(expr_mat,labels=NA, suppress.plot=FALSE, threshold=NA) {
	num_detected =  colSums(expr_mat > 0);
	if (!is.na(threshold)) {
		low_quality = num_detected < threshold;
	} else {
		num_zero = colSums(expr_mat == 0);
		cell_zero = num_zero/length(expr_mat[,1]);
		mu = mean(cell_zero);
		sigma = sd(cell_zero);
		# Deal with bi-modal
		if (sum(cell_zero > mu-sigma & cell_zero < mu+sigma) < 0.5) { # should be 0.68 theoretically
			mu = mean(cell_zero[cell_zero < median(cell_zero)]);
			sigma = sd(cell_zero[cell_zero < median(cell_zero)]);
		}
		low_quality = p.adjust(pnorm((cell_zero-mu)/sigma, lower.tail=F), method="fdr") < 0.05;
		if (!suppress.plot) {
			hist(cell_zero, col="grey75", xlab="Number of zeros (per cell)", main="", prob=TRUE)
			curve(dnorm(x,mean=mu, sd=sigma), add=TRUE)
			if (sum(low_quality) > 0) {
				abline(v=min(cell_zero[low_quality]), col="red")
			}
		}
	}
	if (sum(low_quality) > 0) {
		expr_mat = expr_mat[,!low_quality];
		cell_zero = cell_zero[!low_quality];
		if (!is.na(labels)) {labels = labels[!low_quality]}
	}
	return(list(expr_mat = expr_mat, labels = labels));
}

#bg__normalize <- function(data) {
#	# Combine UQ and detection rate adjusted normalization 
#	# Stephanie Hick, Mingziang Teng, Rafael A Irizarry "On the widespread and critical impact of systematic single-cell RNA-Seq data" http://dx.doi.org/10.1101/025528 
#	cell_zero = colSums(data == 0)/length(data[,1]);
#	uq = unlist(apply(data,2,bg__UQ));
#	normfactor = (uq/median(uq)) * (median(cell_zero)/cell_zero); 
#	data = t(t(data)/normfactor);
#	return(data);
#}

# DE Genes functions

bg__test_DE_K_equiv <- function (expr_mat, fit=NA) {
	gene_info = bg__calc_variables(expr_mat);
	if (is.na(fit)) {
		fit = bg__fit_MM(gene_info$p, gene_info$s);
	}
	p_obs = gene_info$p;
	N = length(expr_mat[1,]);
	p_err = gene_info$p_stderr;
	S_mean = gene_info$s
	S_err = gene_info$s_stderr
	K_err = fit$Kerr;
	K_equiv = p_obs*S_mean/(1-p_obs);
	K_equiv_err = p_obs/(1-p_obs)*S_err
		
	Z = (K_equiv - fit$K)/sqrt(K_equiv_err^2+K_err^2); # high = shifted right, low = shifted left
	pval = pnorm(Z, lower.tail=F)
	effect_size = K_equiv/fit$K;
	return(list(pval = pval, fold_change = effect_size))
}
# Use the fact that errors of proportions are well define by converting S to proportion detected equivalents?
bg__test_DE_P_equiv <- function (expr_mat,  fit=NA) {
	gene_info = bg__calc_variables(expr_mat);
	if (is.na(fit)) {
		fit = bg__fit_MM(gene_info$p, gene_info$s);
	}
	p_obs = gene_info$p;
	N = length(expr_mat[1,]);
	p_err = gene_info$p_stderr;
	S_mean = gene_info$s
	S_err = gene_info$s_stderr
	K_err = fit$Kerr;
	p_equiv = fit$predictions;
	propagated_err_p_equiv = p_equiv*sqrt(((S_err+K_err)/(S_mean+fit$K))^2+(K_err/fit$K)^2)
	fitted_err_p_equiv = fit$fitted_err
#	Z = (p_equiv - p_obs)/sqrt(p_err^2+propagated_err_p_equiv^2); # low = shifted right, high = shifted left
	Z = (p_equiv - p_obs)/fitted_err_p_equiv; # low = shifted right, high = shifted left
	pval = pnorm(Z, lower.tail=T)
	effect_size = p_obs/p_equiv;
	return(list(pval = pval, fold_change = effect_size))
}

# Use the fact that S as a function of P is more stable to noise for the main part of the curve
bg__test_DE_S_equiv <- function (expr_mat, fit=NA, method="propagate") {
	gene_info = bg__calc_variables(expr_mat);
	if (is.na(fit[1])) {
		fit = bg__fit_MM(gene_info$p, gene_info$s);
	}
	p_obs = gene_info$p;
	N = length(expr_mat[1,]);
	p_err = gene_info$p_stderr;
	S_mean = gene_info$s
	S_err = gene_info$s_stderr
	K_err = fit$Kerr;
	S_equiv = bg__invert_MM(fit$K,p_obs);

	## Monte Carlo method to estimate error around S_equiv ##
	MC_err <- function (p_base) {
		p_rand = rnorm(10000, mean = p_base, sd = p_err);
		p_rand = p_rand[p_rand > 0 & p_rand < 1]
		K_rand = rnorm(length(p_rand),fit$K,sd = K_err);
		K_rand[K_rand < 1] = 1;
		S_equiv_rand = bg__invert_MM(K_rand, p_rand)
		sd(S_equiv_rand)
	}
	if (method == "MC") {
		S_equiv_err = unlist(lapply(p_obs,MC_err))
	} else {
		S_equiv_err = S_equiv*sqrt(2*(p_err/p_obs)^2+(K_err/fit$K)^2);
	}

	Z = (S_equiv - S_mean)/sqrt(S_err^2+S_equiv_err^2); # low = shifted right, high = shifted left
	pval = pnorm(Z, lower.tail=T)*2
	effect_size = (S_mean-S_equiv)/S_equiv;
	return(list(pval = pval, effect = effect_size))
}

bg__get_extreme_residuals <- function (expr_mat, fit=NA, v_threshold=c(0.05,0.95), perc_most_extreme = NA, fdr_threshold = 0.1, direction="right", suppress.plot = FALSE) {
	gene_info = bg__calc_variables(expr_mat);
	if (is.na(fit)) {
		fit = bg__fit_MM(gene_info$p, gene_info$s);
	}
	res = bg__horizontal_residuals_MM_log10(fit$K, gene_info$p, gene_info$s)
	res = res[gene_info$p < max(v_threshold) & gene_info$p > min(v_threshold)]

	if (is.na(perc_most_extreme)) {
		mu = mean(res); sigma = sd(res);
		# deal with potential bi-modality
		if (sum(res > mu-sigma & res < mu+sigma) < 0.5) { # should be 0.68 theoretically
			mu = mean(res[res > quantile(res,0.33)]);
			sigma = sd(res[res > quantile(res,0.33)]);
		}
	
		if (direction=="right") {
			pval =pnorm((res-mu)/sigma, lower.tail=F)
		} else {
			pval = pnorm((res-mu)/sigma, lower.tail=T)
		}
		qval = p.adjust(pval, method="fdr");
		sig = qval < fdr_threshold;

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
			cut_off = quantile(res,prob=1-perc_most_extreme);
			return(names(res)[res > cut_off]);
		} else {
			cut_off = quantile(res,prob=perc_most_extreme);
			return(names(res)[res < cut_off]);
		}
	}
}
##### Assembled Analysis Chunks ####
M3D_Clean_Data <- function(expr_mat, labels = NA, is.counts=TRUE, suppress.plot=FALSE, pseudo_genes=NA, min_detected_genes=NA) {
	if (length(pseudo_genes) > 1) {
		is_pseudo = rownames(expr_mat) %in% pseudo;
	        expr_mat = expr_mat[!is_pseudo,];
	}

	expr_mat = bg__filter_cells(expr_mat, labels, suppress.plot = suppress.plot, threshold=min_detected_genes);

        detected = rowSums(expr_mat > 0) > 3;
        expr_mat = expr_mat[detected,];

	spikes = grep("ercc",rownames(expr_mat), ignore.case=TRUE)
	if (is.counts) {
                totreads = colSums(expr_mat[-c(spikes),])
                cpm = t(t(expr_mat)/totreads)*1000000;
                lowExpr = rowMeans(cpm) < 10^-5;
                cpm=cpm[!lowExpr,];
                return(list(data=cpm, labels=labels));
        }

	lowExpr = rowMeans(expr_mat) < 10^-5;
        data=expr_mat[!lowExpr,];
        return(list(data=expr_mat, labels=labels));
}

M3D_Dropout_Models <- function(expr_mat, xlim=NA) {
	BasePlot = bg__dropout_plot_base(expr_mat, xlim = xlim);
	MM = bg__fit_MM(BasePlot$P, BasePlot$S);
	SCDE = bg__fit_logistic(BasePlot$P, BasePlot$S);
	ZIFA = bg__fit_ZIFA(BasePlot$P, BasePlot$S);
	sizeloc = bg__add_model_to_plot(MM, BasePlot, lty=1, lwd=2.5, col="black",legend_loc = "topright");
	sizeloc = bg__add_model_to_plot(SCDE, BasePlot, lty=2, lwd=2.5, col="grey35",legend_loc = c(sizeloc$rect$left+sizeloc$rect$w,sizeloc$rect$top-sizeloc$rect$h-0.05));
	sizeloc = bg__add_model_to_plot(ZIFA, BasePlot, lty=3, lwd=2.5, col="grey65",legend_loc = c(sizeloc$rect$left+sizeloc$rect$w,sizeloc$rect$top-sizeloc$rect$h-0.05));
	invisible(list(MMfit = MM, LogiFit = SCDE, ExpoFit = ZIFA));
}

M3D_Differential_Expression <- function(expr_mat, mt_method="bon", mt_threshold=0.05) {
	BasePlot = bg__dropout_plot_base(expr_mat, xlim = NA);
	MM = bg__fit_MM(BasePlot$P, BasePlot$S);
	sizeloc = bg__add_model_to_plot(MM, BasePlot, lty=1, lwd=2.5, col="black",legend_loc = "topright");
	DEoutput = bg__test_DE_K_equiv(expr_mat, fit=MM);

	sig = which(p.adjust(DEoutput$pval, method=mt_method) < mt_threshold);
	DEgenes = rownames(expr_mat)[sig];
	DEgenes = DEgenes[!is.na(DEgenes)];
	bg__highlight_genes(BasePlot, DEgenes);
	
	TABLE = data.frame(Gene = DEgenes, p.value = DEoutput$pval[sig], q.value= p.adjust(DEoutput$pval, method=mt_method)[sig])
	return(TABLE)
}

M3D_Expression_Heatmap <- function(Genes, expr_mat, cell_labels=NA, interesting_genes=NA, marker_genes=NA, outlier_cells=NA) {
	# Converted known DE genes into heatmap labels 
	gene_labels = rep(1, times = length(Genes));
	if (is.na(interesting_genes)) {
		gene_labels=NA
	}
 	if (is.list(interesting_genes)) {
                for (i in 1:length(interesting_genes)) {
                        gene_labels[Genes %in% interesting_genes[[i]]] = i+1;
                }
        } else {
                gene_labels[Genes %in% interesting_genes] = 2;
        }
	if (is.numeric(marker_genes) | is.logical(marker_genes)) {
		marker_genes = rownames(expr_mat)[marker_genes];
	}
	if (is.numeric(outlier_cells) | is.logical(outlier_cells)) {
		outlier_cells = rownames(expr_mat)[outlier_cells];
	}
	heatmap_output = bg__expression_heatmap(Genes, expr_mat, cell_labels=cell_labels, gene_labels=as.numeric(gene_labels), key_genes=as.character(marker_genes), key_cells=outlier_cells);
	invisible(heatmap_output);
}

M3D_Get_Extremes <- function(expr_mat, fdr_threshold = 0.1, percent = NA, v_threshold=c(0.05,0.95)) {
	BasePlot = bg__dropout_plot_base(expr_mat, xlim = NA);
	MM = bg__fit_MM(BasePlot$P, BasePlot$S);
	sizeloc = bg__add_model_to_plot(MM, BasePlot, lty=1, lwd=2.5, col="black",legend_loc = "topright");
	if (is.na(percent)) {
		shifted_right = bg__get_extreme_residuals(expr_mat, fit=MM, v_threshold=v_threshold, fdr_threshold = fdr_threshold, direction="right", suppress.plot=TRUE)
		shifted_left  = bg__get_extreme_residuals(expr_mat, fit=MM, v_threshold=v_threshold, fdr_threshold = fdr_threshold, direction="left",  suppress.plot=TRUE)
	} else {
		shifted_right = bg__get_extreme_residuals(expr_mat, fit=MM, v_threshold=v_threshold, perc_most_extreme = percent, direction="right", suppress.plot=TRUE)
		shifted_left  = bg__get_extreme_residuals(expr_mat, fit=MM, v_threshold=v_threshold, perc_most_extreme = percent, direction="left",  suppress.plot=TRUE)

	}
	bg__highlight_genes(BasePlot, shifted_right, colour="orange");
	bg__highlight_genes(BasePlot, shifted_left, colour="purple");
	return(list(left=shifted_left,right=shifted_right));
}

###### Extra Stuff for Comparisons ######
#from : http://www.nature.com/nmeth/journal/v10/n11/full/nmeth.2645.html#supplementary-information
Brennecke_getVariableGenes <- function(expr_mat, spikes=NA, suppress.plot=FALSE, fdr=0.1, minBiolDisp=0.5) {
        require(statmod)

        rowVars <- function(x) { unlist(apply(x,1,var))}

        colGenes = "black"
        colSp = "grey35"


        fullCountTable <- expr_mat;

        if (is.character(spikes)) {
                sp = rownames(fullCountTable) %in% spikes;
                countsSp <- fullCountTable[sp,];
                countsGenes <- fullCountTable[!sp,];
        } else if (is.numeric(spikes)) {
                countsSp <- fullCountTable[spikes,];
                countsGenes <- fullCountTable[-spikes,];
        } else {
                countsSp = fullCountTable;
                countsGenes = fullCountTable;
        }

        meansSp = rowMeans(countsSp)
        varsSp = rowVars(countsSp)
        cv2Sp = varsSp/meansSp^2
        meansGenes = rowMeans(countsGenes)
        varsGenes = rowVars(countsGenes)
        cv2Genes = varsGenes/meansGenes^2
        # Fit Model
        minMeanForFit <- unname( quantile( meansSp[ which( cv2Sp > 0.3 ) ], 0.80))
        useForFit <- meansSp >= minMeanForFit
        if (sum(useForFit) < 50) {
                warning("Too few spike-ins exceed minMeanForFit, recomputing using all genes.")
                meansAll = c(meansGenes, meansSp)
                cv2All = c(cv2Genes,cv2Sp)
                minMeanForFit <- unname( quantile( meansAll[ which( cv2All > 0.3 ) ], 0.80))
                useForFit <- meansSp >= minMeanForFit
        }
        if (sum(useForFit) < 50) {warning(paste("Only", sum(useForFit), "spike-ins to be used in fitting, may result in poor fit."))}
        fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/meansSp[useForFit] ), cv2Sp[useForFit] )
        a0 <- unname( fit$coefficients["a0"] )
        a1 <- unname( fit$coefficients["a1tilde"])

        # Test
        psia1theta <- a1
        minBiolDisp <- minBiolDisp^2
        m = ncol(countsSp);
        cv2th <- a0 + minBiolDisp + a0 * minBiolDisp
        testDenom <- (meansGenes*psia1theta + meansGenes^2*cv2th)/(1+cv2th/m)
        p <- 1-pchisq(varsGenes * (m-1)/testDenom,m-1)
        padj <- p.adjust(p,"BH")
        sig <- padj < fdr
        sig[is.na(sig)] <- FALSE
        if (!suppress.plot) {
                plot( meansGenes,cv2Genes, xaxt="n", yaxt="n", log="xy",
                        xlab = "average normalized read count",
                        ylab = "squared coefficient of variation (CV^2)", col="white")
                axis( 1, 10^(-2:5), c( "0.01", "0.1", "1", "10", "100", "1000",
                        expression(10^4), expression(10^5) ) )
                axis( 2, 10^(-2:3), c( "0.01", "0.1", "1", "10", "100","1000" ), las=2 )
                abline( h=10^(-2:1), v=10^(-1:5), col="#D0D0D0", lwd=2 )
                # Plot the genes, use a different color if they are highly variable
                points( meansGenes, cv2Genes, pch=20, cex=.2,
                        col = ifelse( padj < .1, "#C0007090", colGenes ) )
                # Add the technical noise fit
                xg <- 10^seq( -2, 6, length.out=1000 )
                lines( xg, (a1)/xg + a0, col="#FF000080", lwd=3 )
                # Add a curve showing the expectation for the chosen biological CV^2 thershold
                lines( xg, psia1theta/xg + a0 + minBiolDisp, lty="dashed", col="#C0007090", lwd=3)
        }
        return(names(meansGenes)[sig])
}

