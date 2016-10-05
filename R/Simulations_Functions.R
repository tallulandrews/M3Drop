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

hidden_rowVars <- function(x) {
	unlist(apply(x, 1, var, na.rm = T))
}
hidden_colVars <- function(x) {
	unlist(apply(x, 2, var, na.rm = T))
}

bg__mean2disp <- function(mu, coeffs=c(3.967816,-1.855054)){
		cv2 <- exp(coeffs[1]+coeffs[2]*(log(mu)/log(10)))
		variance <- cv2*(mu^2)
		disp <- mu^2/(variance-mu)
		return(1/disp)
}

M3DropMakeSimData <- function(dispersion_fun=bg__mean2disp, n_cells=300, dispersion_factor=1, base_means=10^rnorm(25000,1,1), K=10.3) {
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

M3DropMakeSimDE <- function(dispersion_fun=bg__mean2disp, fold_change=10, frac_change=0.1, n_cells=300, sub_pop=0.5, dispersion_factor=1, base_means=10^rnorm(25000,1,1), K=10.3){
        n_genes <- length(base_means);
        TP <- sample(1:n_genes, frac_change*n_genes)
        sub_pop <- round(sub_pop*n_cells)
        Pop_lab <- c(rep(1, times=n_cells-sub_pop),rep(2,times=sub_pop))
	# Base Population
	base <- M3DropMakeSimData(dispersion_fun=dispersion_fun, 
			n_cells=sum(Pop_lab==1), 
			dispersion_factor=dispersion_factor, 
			base_means=base_means, K=K)
	changed_means <- base_means;
	changed_means[TP] <- base_means[TP]*fold_change;
	# Changed Subpopulation
	sub_pop <- M3DropMakeSimData(dispersion_fun=dispersion_fun, 
			n_cells=sum(Pop_lab==2), 
			dispersion_factor=dispersion_factor, 
			base_means=changed_means, K=K)
	return(list(data=cbind(base,sub_pop), cell_labels=Pop_lab, TP=TP));
}


M3DropMakeSimDVar <- function(dispersion_fun=bg__mean2disp, fold_change=10, frac_change=0.1, n_cells=300, sub_pop=0.5, dispersion_factor=1, base_means=10^rnorm(25000, 1,1), K=10.3) {
	# Make Simulated Matrix
        n_genes <- length(base_means);
        TP <- sample(1:n_genes, frac_change*n_genes)
        sub_pop <- round(sub_pop*n_cells)
        Pop_lab <- c(rep(1, times=n_cells-sub_pop),rep(2,times=sub_pop))
	# Whole Population
	base <- M3DropMakeSimData(dispersion_fun=dispersion_fun, 
			n_cells=n_cells, 
			dispersion_factor=dispersion_factor, 
			base_means=base_means, K=K)
	# Changed Vals
	subpop <- M3DropMakeSimData(dispersion_fun=dispersion_fun, 
			n_cells=sum(Pop_lab==2), 
			dispersion_factor=fold_change*dispersion_factor, 
			base_means=base_means[TP], K=K)
	base[TP,Pop_lab==2] <- subpop
	return(list(data=base, cell_labels=Pop_lab, TP=TP));
}

M3DropMakeSimHVar <- function(dispersion_fun=bg__mean2disp, fold_change=10, frac_change=0.1, n_cells=300, dispersion_factor=1, base_means=10^rnorm(25000, 1,1), K=10.3) {
        n_genes <- length(base_means);
        TP <- sample(1:n_genes, frac_change*n_genes)
        Pop_lab <- rep(1, times=n_cells)
	# Whole Population
	base <- M3DropMakeSimData(dispersion_fun=dispersion_fun, 
			n_cells=n_cells, 
			dispersion_factor=dispersion_factor, 
			base_means=base_means, K=K)
	# Changed Vals
	subpop <- M3DropMakeSimData(dispersion_fun=dispersion_fun, 
			n_cells=n_cells, 
			dispersion_factor=fold_change*dispersion_factor, 
			base_means=base_means[TP], K=K)
	base[TP,] <- subpop
	return(list(data=base, cell_labels=Pop_lab, TP=TP));
}

bg__get_stats <- function(sig, TP, ngenes) {
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

bg__var_vs_drop <- function(pop_size, fixed_mean, K=10.3, suppress.plot=TRUE) {
	# Relationship between Fold Change and Var/Dropouts for fixed mean
        fc <- seq(from=1, to=100, by=1)
        labels <- c(rep(1, times=pop_size),rep(2,times=pop_size))
        lowmean_fun <- function(fc) {2*fixed_mean/(1+fc)}
        test <- sapply(fc, function(f) {
                low_mean <- lowmean_fun(f)
                high_mean <- low_mean*f
                base <- rnbinom(pop_size, size=1/bg__mean2disp(low_mean), mu=low_mean);
                subpop <- rnbinom(pop_size, size=1/bg__mean2disp(high_mean), mu=high_mean)
                base <- hidden_add_dropouts(base,low_mean,K)
                subpop <- hidden_add_dropouts(subpop,high_mean*f,K)
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

