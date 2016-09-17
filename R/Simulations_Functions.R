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

hidden_rowVars <- function(x) {
	unlist(apply(x, 1, var, na.rm = T))
}
hidden_colVars <- function(x) {
	unlist(apply(x, 2, var, na.rm = T))
}

M3DropMean2Disp <- function(mu, coeffs=c(3.967816,-1.855054)){
		cv2 <- exp(coeffs[1]+coeffs[2]*(log(mu)/log(10)))
		variance <- cv2*(mu^2)
		disp <- mu^2/(variance-mu)
		return(1/disp)
}

M3DropMakeSimData <- function(dispersion_fun, n_cells=300, dispersion_factor=1, base_means=10^runif(15000, min=-4, max=4), K=NULL) {
	# Make Simulated Matrix
        n_genes <- length(base_means);
        expr_mat <- sapply(1:n_genes, function(x){
                    base <- rnbinom(n_cells, 
			size<-dispersion_factor/dispersion_fun(base_means[x]), 
			mu<-base_means[x])
                    if (!is.null(K)) {
                         base <- add_dropouts(base,base_means[x],K)
                    }
                    return(base)
                    })
        expr_mat <- t(expr_mat);
        rownames(expr_mat) <- 1:length(expr_mat[,1])
        colnames(expr_mat) <- 1:length(expr_mat[1,])
	return(expr_mat);
}

M3DropMakeSimDE <- function(dispersion_func=M3DropMean2Disp, fold_change=10, frac_change=0.1, n_cells=300, sub_pop=0.5, dispersion_factor=1, base_means=10^rnorm(25000,1,1), K=10.3){
        n_genes <- length(base_means);
        TP <- sample(1:n_genes, frac_change*n_genes)
        sub_pop <- round(sub_pop*n_cells)
        Pop_lab <- c(rep(1, times=n_cells-sub_pop),rep(2,times=sub_pop))
	# Base Population
	base <- M3DropMakeSimData(dispersion_fun=disperion_func, 
			n_cells=sum(Pop_lab==1), 
			dispersion_factor=dispersion_factor, 
			base_means=base_means, K=K)
	changed_means <- base_means;
	changed_means[TP] <- base_means[TP]*fold_change;
	# Changed Subpopulation
	sub_pop <- M3DropMakeSimData(dispersion_fun=disperion_func, 
			n_cells=sum(Pop_lab==2), 
			dispersion_factor=dispersion_factor, 
			base_means=changed_means, K=K)
	return(list(data=cbind(base,sub_pop), cell_labels=Pop_lab, TP=TP));
}


M3DropMakeSimDVar <- function(dispersion_fun=M3DropMean2Disp, fold_change=10, frac_change=0.1, n_cells=300, sub_pop=0.5, dispersion_factor=1, base_means=10^rnorm(25000, 1,1), K=10.3) {
	# Make Simulated Matrix
        n_genes <- length(base_means);
        TP <- sample(1:n_genes, frac_change*n_genes)
        sub_pop <- round(sub_pop*n_cells)
        Pop_lab <- c(rep(1, times=n_cells-sub_pop),rep(2,times=sub_pop))
	# Whole Population
	base <- M3DropMakeSimData(dispersion_fun=disperion_func, 
			n_cells=n_cells, 
			dispersion_factor=dispersion_factor, 
			base_means=base_means, K=K)
	# Changed Vals
	subpop <- M3DropMakeSimData(dispersion_fun=disperion_func, 
			n_cells=sum(Pop_lab==2), 
			dispersion_factor=fold_change*dispersion_factor, 
			base_means=base_means[TP], K=K)
	base[TP,Pop_lab==2] <- subpop
	return(list(data=base, cell_labels=Pop_lab, TP=TP));
}

M3DropMakeSimHVar <- function(dispersion_fun=M3DropMean2Disp, fold_change=10, frac_change=0.1, n_cells=300, dispersion_factor=1, base_means=10^rnorm(25000, 1,1), K=10.3) {
        n_genes <- length(base_means);
        TP <- sample(1:n_genes, frac_change*n_genes)
        Pop_lab <- rep(1, times=n_cells)
	# Whole Population
	base <- M3DropMakeSimData(dispersion_fun=disperion_func, 
			n_cells=n_cells, 
			dispersion_factor=dispersion_factor, 
			base_means=base_means, K=K)
	# Changed Vals
	subpop <- M3DropMakeSimData(dispersion_fun=disperion_func, 
			n_cells=n_cells, 
			dispersion_factor=fold_change*dispersion_factor, 
			base_means=base_means[TP], K=K)
	base[TP,] <- subpop
	return(list(data=base, cell_labels=Pop_lab, TP=TP));
}

bg_get_stats <- function(sig, TP, ngenes) {
	TPs<-sum(sig %in% TP)
	FPs<-length(sig)-TPs
	FNs<-length(TP)-TPs
	TNs<-ngenes-length(TP)-FPs;
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

M3DropStatsForSim <- function(sim_mat, TP, Observed_Means, mt_threshold=0.05, suppress.plot=TRUE) {
	require("M3Drop")
        mt_method<-"fdr"

	# Remove Undetectable Genes
#       sim_mat = sim_mat[rowSums(sim_mat) > 0,]
        sim_mat <- sim_mat[rowSums(sim_mat == 0) > 0,]
        TP <- TP[TP %in% rownames(sim_mat)]

        DE <- M3Drop_Differential_Expression(sim_mat, mt_method = mt_method, mt_threshold=2, suppress.plot=TRUE)
        HVG <- BrenneckeGetVariableGenes(sim_mat, fdr=2, suppress.plot=TRUE)

        output <- rbind(c(0,0,0),c(0,0,0));
	colnames(output) <- c("AUC","FDR","FNR")
	rownames(output) <- c("M3Drop","HVG")
        tmp <- list(DE, HVG)

        # Plot Stuff - setup #
        act_means <- rowMeans(sim_mat)
        thing <- log(act_means)/log(10)
        bins <- 10^seq(from = floor(min(thing)), to = ceiling(max(thing)), by = 1)

	for (bf in 1:9) { # Shift bins to make curve smoother
                binned <- as.numeric(cut(rowMeans(sim_mat), bins*bf))
                names(binned) <- rownames(sim_mat)
                plot_output <- bins[2:length(bins)]*bf

	        for (i in 1:2) { # Do same thing for HVG & M3Drop
	                Diff <- tmp[[i]]
	                sig <- Diff[Diff$q.value < mt_threshold,]


	                stats <- get_stats(sig[,1], TP, ngenes=length(sim_mat[,1]));

	                Diff$Truth <- rep(0, times=length(Diff[,1]))
	                Diff[Diff[,1] %in% TP,]$Truth <- 1;
	                pred <- ROCR::prediction(1-Diff$p.value, Diff$Truth)
	                val <- unlist(ROCR::performance(pred,"auc")@y.values)
	                output[i,]<-c(val,stats);

	                # Plot Stuff - data #
	                get_bin_stats <- function(b){
	                        this_bin <- names(binned)[binned==b];
	                        this_sig <- this_bin[this_bin %in% sig[,1]];
	                        this_TP <- this_bin[this_bin %in% TP];
	                        out <- get_stats(this_sig, this_TP, ngenes = length(sim_mat[,1]));
	                        c(out, length(this_sig), length(this_TP))
	                }
	                bin_stats <- sapply(1:(length(bins)-1), get_bin_stats)
	                rownames(bin_stats) <- c("FDR","FNR","Ncalled","nTrue")
	                plot_output <- rbind(plot_output, bin_stats)
	        }
                if (bf==1) {
                        final_plot_output <- plot_output;
                } else {
                        final_plot_output <- cbind(final_plot_output, plot_output);
                }

        }

	final_plot_output[1,] <- final_plot_output[1,]/2
	if (!suppress.plot) {
                background <- density(log(Observed_Means)/log(10))
                drawing<-list()
                drawing$ylim<-c(0,max(c(background$y)));
                drawing$xlim<-c(min(background$x),max(c(background$x,ceiling(log(final_plot_output[1,])/log(10)))));
                plot(1, col="white", xlim=drawing$xlim, ylim=drawing$ylim, xaxt="n", yaxt="n", xlab="", ylab=""); polygon(background, col="grey75", border="grey75")


                convert_yvals <- function(y) { # convert 0-1 to appropriate scale for plot
                        y*drawing$ylim[2]

                }
                final_plot_output <- final_plot_output[,order(final_plot_output[1,])]


                xes <- log(final_plot_output[1,])/log(10)

                ticks <- seq(from=1, to =length(xes), by=9)
		ticks <- final_plot_output[1,ticks]*2
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

bg_var_vs_drop <- function(pop_size, fixed_mean, suppress.plot=TRUE) {
	# Relationship between Fold Change and Var/Dropouts for fixed mean
        fc <- seq(from=1, to=100, by=1)
        labels <- c(rep(1, times=pop_size),rep(2,times=pop_size))
        lowmean_fun <- function(fc) {2*fixed_mean/(1+fc)}
        test <- sapply(fc, function(f) {
                low_mean <- lowmean_fun(f)
                high_mean <- low_mean*f
                base <- rnbinom(pop_size, size=1/dispersion_from_mean(low_mean), mu=low_mean);
                subpop <- rnbinom(pop_size, size=1/dispersion_from_mean(high_mean), mu=high_mean)
                base <- add_dropouts(base,low_mean,K)
                subpop <- add_dropouts(subpop,high_mean*f,K)
                return(c(base,subpop))
        })
        var_btw_fun <- function(x) {
                a <- aov(x~labels)
                sse_btw <-unlist(summary(a))[3]
                var_btw <- sse_btw/(2*pop_size-1)
                return(var_btw)
        }
        var_within_fun <- function(x) {
                a <- aov(x~labels)
                sse_within <-unlist(summary(a))[4]
                var_within <- sse_within/(2*pop_size-1)
                return(var_within)
        }
        Vbtw <- apply(test,2,var_btw_fun)
        Vwithin <- apply(test,2,var_within_fun)
        vars <- apply(test, 2, var)
        drops <- apply(test, 2, function(x){sum(x==0)})/length(test[,1])
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
        var_r <- cor(vars,fc)
        drop_r <- cor(drops,fc)
        return(list(var_r=var_r, drop_r=drop_r, vars=vars, drops=drops, fc=fc, Vbtw=Vbtw, Vwithin=Vwithin));
}

