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

### Functions ###
hidden_add_dropouts <- function(x,mu,K){
	p_drop <- 1 - mu/(mu+K);
	expect_pos <- mu*length(x)/sum(x>0);
	p_drop <- K/expect_pos
#	p_new_drop = (p_drop*length(x)-sum(x==0))/sum(x > 0)
	toss <- runif(n=length(x));
	x[toss < p_drop] <- 0;
	return(x);
}

hidden_amplification1 <- function (x, rounds=12, efficiency=0.97) {
	tot = sum(x);
	if (tot == 0) {return(x)} # Nothing to be done if all zero already
	# Amplify
	for (i in 1:rounds) {
		x <- sapply(x, function(y) {y+rbinom(1,size=y, prob=efficiency)})
	}
	amped = sum(x)
	# Downsample back to starting amounts
	x <- sapply(x, function(z) {rbinom(1,size=z, prob=1/((1+efficiency)^rounds))})
	return(x);
}

bg__default_mean2disp <- function(mu, coeffs=c(3.967816,-1.855054)){
		cv2 <- exp(coeffs[1]+coeffs[2]*(log(mu)/log(10)))
		variance <- cv2*(mu^2)
		disp <- mu^2/(variance-mu)
		return(1/disp)
}

bg__MakeSimData <- function(dispersion_fun=bg__default_mean2disp, n_cells=300, dispersion_factor=1, base_means=10^rnorm(25000,1,1), K=10.3) {
	# Make Simulated Matrix
        n_genes <- length(base_means);
        expr_mat <- sapply(1:n_genes, function(x){
                    base <- rnbinom(n_cells, 
			size=1/(dispersion_factor*dispersion_fun(base_means[x])), 
			mu=base_means[x])
                    if (!is.null(K)) {
                         base <- hidden_add_dropouts(base,base_means[x],K)
                    }
                    return(base)
                    })
        expr_mat <- t(expr_mat);
        rownames(expr_mat) <- 1:length(expr_mat[,1])
        colnames(expr_mat) <- 1:length(expr_mat[1,])
	return(expr_mat);
}

bg__MakeSimDE <- function(dispersion_fun=bg__default_mean2disp, fold_change=10, frac_change=0.1, n_cells=300, sub_pop=0.5, dispersion_factor=1, base_means=10^rnorm(25000,1,1), K=10.3){
        n_genes <- length(base_means);
        TP <- sample(1:n_genes, frac_change*n_genes)
        sub_pop <- round(sub_pop*n_cells)
        Pop_lab <- c(rep(1, times=n_cells-sub_pop),rep(2,times=sub_pop))
	# Base Population
	base <- bg__MakeSimData(dispersion_fun=dispersion_fun, 
			n_cells=sum(Pop_lab==1), 
			dispersion_factor=dispersion_factor, 
			base_means=base_means, K=K)
	changed_means <- base_means;
	changed_means[TP] <- base_means[TP]*fold_change;
	# Changed Subpopulation
	sub_pop <- bg__MakeSimData(dispersion_fun=dispersion_fun, 
			n_cells=sum(Pop_lab==2), 
			dispersion_factor=dispersion_factor, 
			base_means=changed_means, K=K)
	return(list(data=cbind(base,sub_pop), cell_labels=Pop_lab, TP=TP));
}


bg__MakeSimDVar <- function(dispersion_fun=bg__default_mean2disp, fold_change=10, frac_change=0.1, n_cells=300, sub_pop=0.5, dispersion_factor=1, base_means=10^rnorm(25000, 1,1), K=10.3) {
	# Make Simulated Matrix
        n_genes <- length(base_means);
        TP <- sample(1:n_genes, frac_change*n_genes)
        sub_pop <- round(sub_pop*n_cells)
        Pop_lab <- c(rep(1, times=n_cells-sub_pop),rep(2,times=sub_pop))
	# Whole Population
	base <- bg__MakeSimData(dispersion_fun=dispersion_fun, 
			n_cells=n_cells, 
			dispersion_factor=dispersion_factor, 
			base_means=base_means, K=K)
	# Changed Vals
	subpop <- bg__MakeSimData(dispersion_fun=dispersion_fun, 
			n_cells=sum(Pop_lab==2), 
			dispersion_factor=fold_change*dispersion_factor, 
			base_means=base_means[TP], K=K)
	base[TP,Pop_lab==2] <- subpop
	return(list(data=base, cell_labels=Pop_lab, TP=TP));
}

bg__MakeSimHVar <- function(dispersion_fun=bg__default_mean2disp, fold_change=10, frac_change=0.1, n_cells=300, dispersion_factor=1, base_means=10^rnorm(25000, 1,1), K=10.3) {
        n_genes <- length(base_means);
        TP <- sample(1:n_genes, frac_change*n_genes)
        Pop_lab <- rep(1, times=n_cells)
	# Whole Population
	base <- bg__MakeSimData(dispersion_fun=dispersion_fun, 
			n_cells=n_cells, 
			dispersion_factor=dispersion_factor, 
			base_means=base_means, K=K)
	# Changed Vals
	subpop <- bg__MakeSimData(dispersion_fun=dispersion_fun, 
			n_cells=n_cells, 
			dispersion_factor=fold_change*dispersion_factor, 
			base_means=base_means[TP], K=K)
	base[TP,] <- subpop
	return(list(data=base, cell_labels=Pop_lab, TP=TP));
}

bg__get_stats <- function(sig, TP, ngenes) {
	TPs=sum(sig %in% TP)
	FPs=length(sig)-TPs
	FNs=length(TP)-TPs
	TNs=ngenes-length(TP)-FPs;
	if ((TPs+FPs) == 0) {
		FDR <- 0;
	} else {
		FDR <- FPs/(TPs+FPs);
	}
	if ((FNs+TPs) == 0) {
		FNR <- 0;
	} else {
		FNR <- FNs/(FNs+TPs)
	}
	return(c(FDR,FNR));
}

##### Update #####


obsolete__calc_DE_stats <- function(expr_mat, TP, Observed_Means, mt_threshold=0.05, suppress.plot=TRUE) {
	#require("M3Drop")
        # Test DE
        mt_method="fdr"
	M3Drop_col="black"
	HVG_col="forestgreen"

        expr_mat = expr_mat[Matrix::rowSums(expr_mat) > 0,]
        expr_mat = expr_mat[Matrix::rowSums(expr_mat == 0) > 0,]
        TP = TP[TP %in% rownames(expr_mat)]

        DE = M3DropFeatureSelection(expr_mat, mt_method = mt_method, mt_threshold=2, suppress.plot=TRUE)
        HVG = BrenneckeGetVariableGenes(expr_mat, fdr=2, suppress.plot=TRUE)

        output = rbind(c(0,0,0),c(0,0,0));
	colnames(output) = c("AUC","FDR","FNR")
	rownames(output) = c("M3Drop","HVG")
        tmp = list(DE, HVG)

        # Plot Stuff - setup #
        act_means = rowMeans(expr_mat)
        thing = log(act_means)/log(10)
        bins = 10^seq(from = floor(min(thing)), to = ceiling(max(thing)), by = 1)

	for (bf in 1:9) { # Shift bins to make curve smoother
                binned = as.numeric(cut(rowMeans(expr_mat), bins*bf))
                names(binned) = rownames(expr_mat)
                plot_output = bins[2:length(bins)]*bf

	        for (i in 1:2) { # Do same thing for HVG & M3Drop
	                Diff = tmp[[i]]
	                sig = Diff[Diff$q.value < mt_threshold,]


	                stats = bg__get_stats(sig[,1], TP, ngenes=length(expr_mat[,1]));

	                #require("ROCR")
	                Diff$Truth = rep(0, times=length(Diff[,1]))
	                Diff[Diff[,1] %in% TP,]$Truth = 1;
	                pred <- ROCR::prediction(1-Diff$p.value, Diff$Truth)
	                val <- unlist(ROCR::performance(pred,"auc")@y.values)
	                output[i,]=c(val,stats);

	                # Plot Stuff - data #
	                get_bin_stats <- function(b){
	                        this_bin = names(binned)[binned==b];
	                        this_sig = this_bin[this_bin %in% sig[,1]];
	                        this_TP = this_bin[this_bin %in% TP];
	                        out = bg__get_stats(this_sig, this_TP, ngenes = length(expr_mat[,1]));
	                        c(out, length(this_sig), length(this_TP))
	                }
	                bin_stats = sapply(1:(length(bins)-1), get_bin_stats)
	                #colnames(bin_stats) = bins[2:length(bins)]
	                rownames(bin_stats) = c("FDR","FNR","Ncalled","nTrue")
	                plot_output = rbind(plot_output, bin_stats)
	        }
                if (bf==1) {
                        final_plot_output <- plot_output;
                } else {
                        final_plot_output = cbind(final_plot_output, plot_output);
                }
        }

	final_plot_output[1,] = final_plot_output[1,]/2
	if (!suppress.plot) {
#                par(mar=c(3.5,3.5,4,1))
                background = density(log(Observed_Means)/log(10))
                drawing=list()
                drawing$ylim<-c(0,max(c(background$y)));
#                drawing$xlim<-c(min(background$x),max(background$x));
                drawing$xlim<-c(min(background$x),max(c(background$x,ceiling(log(final_plot_output[1,])/log(10)))));
                plot(1, col="white", xlim=drawing$xlim, ylim=drawing$ylim, xaxt="n", yaxt="n", xlab="", ylab=""); polygon(background, col="grey75", border="grey75")


                convert_yvals = function(y) { # convert 0-1 to appropriate scale for plot
                        y*drawing$ylim[2]

                }
                final_plot_output <- final_plot_output[,order(final_plot_output[1,])]


                xes = log(final_plot_output[1,])/log(10)

                ticks = seq(from=1, to =length(xes), by=9)
		ticks = final_plot_output[1,ticks]*2
#                axis(1, at=xes[ticks], labels = final_plot_output[1,ticks])
                axis(1, at=log(ticks)/log(10), labels = ticks)

                axis(2, at=convert_yvals(seq(from=0, to=1, by=0.1)), labels = seq(from=0, to=1, by=0.1))
                title(xlab="Gene Expression (CPM)", ylab="FDR/FNR", line=2)
                abline(h=convert_yvals(mt_threshold), col="red", lty=2, xpd=F)
                # FDR
                lines(xes,convert_yvals(final_plot_output[2,]), col=M3Drop_col, lty=1, lwd=3)
                lines(xes,convert_yvals(final_plot_output[6,]), col=HVG_col, lwd=3)

                # FNR
                lines(xes,convert_yvals(final_plot_output[3,]), col=M3Drop_col, lwd=3, lty=2)
                lines(xes,convert_yvals(final_plot_output[7,]), col=HVG_col, lwd=3, lty=2)

                par(xpd=T)
                legend("top", inset=c(0,-0.22), c(expression(bold("Method :")), expression(bold("Measure :")), "M3Drop", "FDR","HVG", "FNR"), lty=c(1,1,1,1,1,2), col=c("white","white",M3Drop_col,"black",HVG_col,"black"), lwd=2.5, ncol=3, bty="n")
		par(xpd=F);
        }
	return(list(summary=output, per_expr=final_plot_output))
}

hidden__calc_DE_stats_simplified <- function(expr_mat, TP, Observed_Means, mt_threshold=0.05, suppress.plot=TRUE) {
	# Just ge AUC, FDR, FNR
	#require("M3Drop")
        # Test DE
        mt_method="fdr"

        expr_mat = expr_mat[Matrix::rowSums(expr_mat) > 0,]
        expr_mat = expr_mat[Matrix::rowSums(expr_mat == 0) > 0,]
        TP = TP[TP %in% rownames(expr_mat)]

        DE = M3DropFeatureSelection(expr_mat, mt_method = mt_method, mt_threshold=2, suppress.plot=TRUE)
        HVG = BrenneckeGetVariableGenes(expr_mat, fdr=2, suppress.plot=T)

        output = rbind(c(0,0,0),c(0,0,0));
	colnames(output) = c("AUC","FDR","FNR")
	rownames(output) = c("M3Drop","HVG")
        tmp = list(DE, HVG)

        # Plot Stuff - setup #
        act_means = rowMeans(expr_mat)
        thing = log(act_means)/log(10)
        bins = 10^seq(from = floor(min(thing)), to = ceiling(max(thing)), by = 1)

        for (i in 1:2) { # Do same thing for HVG & M3Drop
                Diff = tmp[[i]]
                sig = Diff[Diff$q.value < mt_threshold,]


                stats = bg__get_stats(sig[,1], TP, ngenes=length(expr_mat[,1]));

		# Need to somehow downsample so distribution of gene expression =~ observed distribution of gene expression
		# or do this when define means? - latter makes more sense b/c do it once for all sim
                #require("ROCR")
                Diff$Truth = rep(0, times=length(Diff[,1]))
                Diff[Diff[,1] %in% TP,]$Truth = 1;
                pred <- ROCR::prediction(1-Diff$p.value, Diff$Truth)
                val <- unlist(ROCR::performance(pred,"auc")@y.values)
                output[i,]=c(val,stats);

        }
	return(list(summary=output, per_expr=NULL))
}

hidden__calc_DE_stats_simplified_singular <- function(expr_mat, TP, DE) {
	qvals <- DE$q.value
	if (!is.numeric(qvals)) {
		qvals <- as.numeric(qvals)
	}
	pvals <- DE$p.value
	if (!is.numeric(pvals)) {
		pvals <- as.numeric(pvals)
	}

        sig = DE[qvals < 0.05,]

        stats = bg__get_stats(sig[,1], TP, ngenes=length(expr_mat[,1]));

        #require("ROCR")
        DE$Truth = rep(0, times=length(DE[,1]))
        DE[DE[,1] %in% TP,]$Truth = 1;
        pred <- ROCR::prediction(1-pvals, DE$Truth)
        val <- unlist(ROCR::performance(pred,"auc")@y.values)
        output=c(val,stats);

	return(list(summary=output, per_expr=NULL))
}

bg__var_vs_drop <- function(pop_size, fixed_mean, K=10.3, dispersion_from_mean=bg__default_mean2disp, suppress.plot=TRUE) {
        fc = seq(from=1, to=100, by=1)
        labels = c(rep(1, times=pop_size),rep(2,times=pop_size))
        lowmean_fun <- function(fc) {2*fixed_mean/(1+fc)}
        test = sapply(fc, function(f) {
                low_mean = lowmean_fun(f)
                high_mean = low_mean*f
                base = rnbinom(pop_size, size=1/dispersion_from_mean(low_mean), mu=low_mean);
                subpop = rnbinom(pop_size, size=1/dispersion_from_mean(high_mean), mu=high_mean)
                base = hidden_add_dropouts(base,low_mean,K)
                subpop = hidden_add_dropouts(subpop,high_mean*f,K)
                return(c(base,subpop))
        })
        var_btw_fun <- function(x) {
                a = aov(x~labels)
                sse_btw =unlist(summary(a))[3]
                var_btw = sse_btw/(2*pop_size-1)
                return(var_btw)
        }
        var_within_fun <- function(x) {
                a = aov(x~labels)
                sse_within =unlist(summary(a))[4]
                var_within = sse_within/(2*pop_size-1)
                return(var_within)
        }
        Vbtw = apply(test,2,var_btw_fun)
        Vwithin = apply(test,2,var_within_fun)
        vars = apply(test, 2, var)
        drops = apply(test, 2, function(x){sum(x==0)})/length(test[,1])
        if (!suppress.plot) {
                par(mar=c(3.5,3.5,4,1))
                # Var plot
                plot(fc, vars, pch=16, xlab="",ylab="", ylim=c(0,max(vars)), type="l")
                title(xlab="Fold Change", ylab="Variance", line=2)
                title(main=paste("mu = ",fixed_mean,", n = ",2*pop_size,sep=""))
                lines(fc, Vbtw, pch=16, xlab="",ylab="", col="blue")
                lines(fc, Vwithin, pch=16, xlab="",ylab="", col="red")
                legend("topleft", c("Total", "Btw","Within"), lty=1, col=c("black","blue","red"), bty="n")
                # Drop plot
                plot(fc, drops, pch=16, xlab="",ylab="")
                title(xlab="Fold Change", ylab="Dropout Rate", line=2)
                title(main=paste("mu = ",fixed_mean,", n = ",2*pop_size,sep=""))
        }
        var_r = cor(vars,fc)
        drop_r = cor(drops,fc)
        return(list(var_r=var_r, drop_r=drop_r, vars=vars, drops=drops, fc=fc, Vbtw=Vbtw, Vwithin=Vwithin));
}

hidden__calc_DE_stats_fpr<- function(expr_mat, TP, Observed_Means, mt_threshold=0.05, suppress.plot=TRUE) {
	#require("M3Drop")
        # Test DE
        mt_method="fdr"

        expr_mat = expr_mat[Matrix::rowSums(expr_mat) > 0,]
        expr_mat = expr_mat[Matrix::rowSums(expr_mat == 0) > 0,]
        TP = TP[TP %in% rownames(expr_mat)]

        DE = M3DropFeatureSelection(expr_mat, mt_method = mt_method, mt_threshold=2, suppress.plot=TRUE)
        HVG = BrenneckeGetVariableGenes(expr_mat, fdr=2, suppress.plot=T)

        output = c(0,0);
        tmp = list(DE, HVG)

        # Plot Stuff - setup #
        act_means = rowMeans(expr_mat)
        thing = log(act_means)/log(10)
        bins = 10^seq(from = floor(min(thing)), to = ceiling(max(thing)), by = 1)

        for (i in 1:2) { # Do same thing for HVG & M3Drop
                Diff = tmp[[i]]
                sig = Diff[Diff$q.value < mt_threshold,]
		val <- 1-sum(sig[,1] %in% TP)/length(sig[,1])
                output[i]=val;
        }
	return(list(summary=output, per_expr=NULL))
}

bg__fit_gamma <- function(x) {
	s = var(x)/mean(x)
	a = mean(x)/s
	return(list(shape=a, scale=s))
}
bg__shift_size <- function(mu_all, size_all, mu_group, coeffs) {
        b <- log(size_all)-coeffs[2]*log(mu_all)
        size_group <- exp(coeffs[2]*log(mu_group)+b)
        return(size_group)
}

NBumiSimulationTrifecta <- function(original_data, n_genes=25000, n_cells=250, sub_pop_prop=0.5) {
	# 31 Aug 2017
	n_cells1 = round(n_cells*2*sub_pop_prop)
	n_cells2 = round(n_cells*2*(1-sub_pop_prop))
#	require("M3Drop")
	fit <- NBumiFitModel(original_data);
	Tis = fit$vals$tis
	Tis_gamma = bg__fit_gamma(Tis)
	Mjs = fit$vals$tjs/fit$vals$nc
	Mjs_norm = c(mean(log(Mjs)/log(10)), sd(log(Mjs)/log(10)))
	mean2disp_coeffs = NBumiFitDispVsMean(fit, suppress.plot=TRUE)
	min_mean = 10^-5;

	g_means = rnorm(n_genes, mean=Mjs_norm[1], sd=Mjs_norm[2]);
	g_means = 10^g_means
	g_means[g_means < min_mean] = min_mean;
	g_means[g_means > max(Mjs)] = max(Mjs)
	g_means1 = g_means*(n_cells1);
	g_means2 = g_means*(n_cells2);

	c_depths1 = round(rgamma(n_cells1, shape=Tis_gamma$shape, scale=Tis_gamma$scale));
	c_depths2 = round(rgamma(n_cells2, shape=Tis_gamma$shape, scale=Tis_gamma$scale));
	l2fc <- rnorm(n_genes, sd=2)

	mus1 <- ((g_means1) %*% t(c_depths1)/sum(c_depths1)) # Fix 11 May 2017
	mus2 <- ((g_means2) %*% t(c_depths2)/sum(c_depths2)) # Fix 11 May 2017
	disp_size <- exp(log(rowMeans(cbind(mus1,mus2)))*mean2disp_coeffs[2]+mean2disp_coeffs[1])

	base <- sapply(1:n_genes, function(i) {sapply(mus1[i,], function(m) {rnbinom(1, mu=m, size=disp_size[i])})})

	shifted_size = bg__shift_size(rowMeans(mus2), disp_size, rowMeans(mus2)*2^l2fc, mean2disp_coeffs)
	de <- sapply(1:n_genes, function(i) {sapply(mus2[i,], function(m) {rnbinom(1, mu=m*2^l2fc[i], size=shifted_size[i])})})

	dv <- sapply(1:n_genes, function(i) {sapply(mus2[i,], function(m) {rnbinom(1, mu=m, size=disp_size[i]*2^l2fc[i])})})
	hv <- sapply(1:n_genes, function(i) {sapply(mus2[i,], function(m) {rnbinom(1, mu=m, size=disp_size[i]*2^l2fc[i])})})

	return(list(truth=l2fc, groups=c(rep(1, times=n_cells1), rep(2, times=n_cells2)), de = cbind(t(base),t(de)), dv = cbind(t(base),t(dv)), hv = cbind(t(dv),t(hv))))

}
M3DropSimulationTrifecta <- function(original_data, n_genes=25000, n_cells=250, sub_pop_prop=0.5) {
	n_cells1 = round(n_cells*2*sub_pop_prop)
	n_cells2 = round(n_cells*2*(1-sub_pop_prop))
#	require("M3Drop")
	tis = Matrix::colSums(original_data)
	norm = t(t(original_data)/tis*median(tis))
	Mjs = rowMeans(norm);
	vals <- bg__calc_variables(norm)
	fit = bg__fit_MM(vals$p, vals$s)
	K = fit$K
	Mjs_norm = c(mean(log(Mjs)/log(10)), sd(log(Mjs)/log(10)))
	mean2disp <- bg__get_mean2disp(norm);
	min_mean = 10^-5;
	
#	g_means = rgamma(n_genes, shape=Mjs_gamma$shape, scale=Mjs_gamma$scale)
	g_means = rnorm(n_genes, mean=Mjs_norm[1], sd=Mjs_norm[2])
	g_means = 10^g_means
	g_means[g_means < min_mean] = min_mean
	g_means[g_means > max(Mjs)] = max(Mjs)
	l2fc <- rnorm(n_genes, sd=2)
	
	# Make datasets
	base <- sapply(g_means, function(mu) {hidden_add_dropouts(rnbinom(n_cells1, mu=mu, size=1/mean2disp(mu)), mu, K)})
	de <- sapply(g_means*2^l2fc, function(mu) {hidden_add_dropouts(rnbinom(n_cells2, mu=mu, size=1/mean2disp(mu)), mu, K)})
	dv <- sapply(1:n_genes, function(i) {
			mu = g_means[i];
			hidden_add_dropouts(rnbinom(n_cells2, mu=mu, size=1/(mean2disp(mu))*2^l2fc[i]), mu, K)
			})
	hv <-  sapply(1:n_genes, function(i) {
                        mu = g_means[i];
                        hidden_add_dropouts(rnbinom(n_cells2, mu=mu, size=1/(mean2disp(mu))*2^l2fc[i]), mu, K)
                        })
	return(list(truth=l2fc, groups=c(rep(1, times=n_cells1), rep(2, times=n_cells2)), de = cbind(t(base),t(de)), dv = cbind(t(base),t(dv)), hv = cbind(t(dv),t(hv))))
}

#Fit_Simulation_Params <- function(original_data) {
#
#}
Make_Sim <- function(model_params, n_genes=25000, pop_params=data.frame(cells=c(125,125), size=c(1,1), hetero=c(0,0), distinct=c(0,0)), dV_params=list(mu=0, sd=1), dE_params=list(mu=0, sd=1)) {
# model_params from above function include:
#	slope & intercept of log-log relationship btw mean & variance
#	mean & sd of distribution of mean expression per gene
#	shape & scale of gamma distribution of total counts per cell (in 100,000s reads)
#	K of MM dropouts
# parameters for each cell population:
#	number of cells
#	relative size of the cells (affects total counts per cell)
#	heterogeneity of the population (multiply FC in variance)a  - setting to 0 = no change in variance.
#	distinctness of the population (multiply FC in mean expression); - setting to 0 = no change in mean expression.
# the "reference" population is never seen, (unless hetero=0, distinct=0, size=1)
# DE & DV are added to every population therefore actual distribution of FC is 
# sum of distribution of FC between two populations, therefore sd=1 is actually sd=2

# Does not support batch effects -> would require one population to have same dV/dE as another but with an additiona effect on top.

}
