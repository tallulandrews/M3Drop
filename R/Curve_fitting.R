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

bg__fit_MM <- function (p,s) {
	if (length(p) != length(s)) {
		stop(print("Error: p and s not same length. Cannot fit Michaelis-Menten."))
	}
	#require("bbmle")
	LL <- function(krt,sigma) {
		R <- p-(1-(s/((krt+s))))
		R <- suppressWarnings(dnorm(R,0,sigma,log=TRUE))
		-sum(R)
	}
	fit <- mle2(LL,start=list(krt=3, sigma=0.25))
	thing <- summary(fit)
	krt <- fit@coef[1]
	res_err <- fit@coef[2];
	Kerr <- max(fit@coef[2],attributes(summary(fit))$coef[1,2]);
	predicted <- 1-(s/(krt+s))
	residuals <- p-predicted
	return(list(K=krt,Kerr=Kerr,fitted_err = res_err,predictions=predicted, model=c("MMenten",paste("K =",round(krt,digits=2))),SSr=round(sum((residuals)^2)),SAr=round(sum(abs(residuals)))))
}

hidden__fit_MM_lognormal<-function(p,s){
	# This consistently underestimates K
	if (length(p) != length(s)) {
		stop(print("Error: p and s not same length. Cannot fit Michaelis-Menten."))
	}
	#require("bbmle")
	p_c <- p[p<1 & p>0]
	s_c <- s[p<1 & p>0]
	LL <- function(krt,sigma) {
		if (krt >0) {
			obs_Ks <- p_c/(1-p_c)*s_c
			R <- log(obs_Ks)-log(krt);
			Qs <- quantile(R,prob=c(0.25,0.5,0.75))
			IQR <- Qs[2]-Qs[1];

			thing <- densCols(p_c, log(s_c)/log(10), colramp = colorRampPalette(c("black","white")));
			dens <- Matrix::colSums(col2rgb(thing))
			thresh <- quantile(dens, prob=0.05);

			R <- R[dens > thresh]
			R <- suppressWarnings(dnorm(R,0,sigma,log=TRUE))
			-sum(R)
		} else {
			10^100
		}
	}
	fit <- mle2(LL,start=list(krt=6, sigma=0.25))
	thing <- summary(fit)
	krt <- fit@coef[1]
	res_err <- fit@coef[2];
	Kerr <- max(fit@coef[2],attributes(summary(fit))$coef[1,2]);
	predicted <- 1-(s/(krt+s))
	residuals <- p-predicted
	return(list(K=krt,Kerr=Kerr,fitted_err = res_err,predictions=predicted, model=c("MMenten",paste("K =",round(krt,digits=3))),SSr=round(sum((residuals)^2)),SAr=round(sum(abs(residuals)))))
}

hidden__fit_MM_logistic <- function(p,s){
	s_nozero <- s[s>0];
	p_nozero <- p[s>0];
	fit <- glm(p_nozero ~ offset(-1*log(s_nozero)), family="binomial")
	Kerr <- summary(fit)$coeff[1,2];
	Kcoeff <- summary(fit)$coeff[1,1];
	krt <- exp(Kcoeff)
	Kerr <- Kerr*krt
	predicted <- rep(0, times=length(s));
	predicted[s>0] <- fitted(fit)
	residuals <- p-predicted
	return(list(K=krt,Kerr=Kerr,predictions=predicted, model=c("MMenten",paste("K =",round(krt,digits=3))),SSr=round(sum((residuals)^2)),SAr=round(sum(abs(residuals)))))
}

bg__fit_logistic <- function(p,s) {
	if (length(p) != length(s)) {
		stop(print("Error: p and s not same length. Cannot fit Logistic Regression."))
	}
	s_nozero <- s[s>0];
	p_nozero <- p[s>0];
        logistic <- suppressWarnings(glm(p_nozero~log(s_nozero),family="binomial")) #warns that is not binary data
        predlog <- fitted(logistic)
	fullpredictions <- rep(0, times=length(s));
	fullpredictions[s>0] <- predlog
	res <- fullpredictions-p;
	return(list(predictions=fullpredictions, B0 = logistic$coeff[1], B1=logistic$coeff[2] ,model=c( "Logistic", paste("Intercept =",round(logistic$coeff[1],digits=3)),paste("Coeff =",round(logistic$coeff[2],digits=3))),SSr=round(sum(res^2)),SAr=round(sum(abs(res)))));
}

bg__fit_ZIFA <- function(p,s) {
	if (length(p) != length(s)) {
		stop(print("Error: p and s not same length. Cannot fit double exponential."))
	}
	p_nozero <- p; p_nozero[p == 0] = min(p[p>0])/10
	fit <- lm(log(p_nozero)~-1+s^2)
	lambda <- -fit$coeff[1]
	res_err <- attributes(summary(fit))$coef[2,1]
	Lerr <- summary(fit)$coeff[1,2]
	predicted <- exp(-lambda*s*s)
	residuals <- p-predicted
	return(list(lambda=lambda,Lerr=Lerr,fitted_err = res_err,predictions=predicted, model=c("p ~ e^(-lambda*S^2)",paste("lambda =",signif(lambda,digits=2))),SSr=round(sum((residuals)^2)),SAr=round(sum(abs(residuals)))))
}

M3DropDropoutModels <- function(expr_mat, xlim=NA, suppress.plot=FALSE) {
	BasePlot <- bg__dropout_plot_base(expr_mat, xlim = xlim, suppress.plot=suppress.plot);
	MM <- bg__fit_MM(BasePlot$gene_info$p, BasePlot$gene_info$s);
	SCDE <- bg__fit_logistic(BasePlot$gene_info$p, BasePlot$gene_info$s);
	ZIFA <- bg__fit_ZIFA(BasePlot$gene_info$p, BasePlot$gene_info$s);
	if (!suppress.plot) {
	  	sizeloc <- bg__add_model_to_plot(MM, BasePlot, lty=1, lwd=2.5, col="black",legend_loc = "topright");
		sizeloc <- bg__add_model_to_plot(SCDE, BasePlot, lty=2, lwd=2.5, col="magenta3",legend_loc = c(sizeloc$rect$left+sizeloc$rect$w,sizeloc$rect$top-sizeloc$rect$h-0.05));
		sizeloc <- bg__add_model_to_plot(ZIFA, BasePlot, lty=3, lwd=2.5, col="red",legend_loc = c(sizeloc$rect$left+sizeloc$rect$w,sizeloc$rect$top-sizeloc$rect$h-0.05));
	}
	invisible(list(MMFit = MM, LogiFit = SCDE, ExpoFit = ZIFA));
}

