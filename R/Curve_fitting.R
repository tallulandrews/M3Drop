bg__fit_MM <- function (p,s) {
#        fit = nls(p ~ 1-(s/((krt+s))),data.frame(s=s),start=list(krt=3))
#	K_glm = glm(p ~ offset(-1*log(s)), family="binomial")
#	Kerr = summary(K_glm)$coeff[1,2];
#	Kcoeff = summary(K_glm)$coeff[1,1];
#	Kerr = exp(Kcoeff+Kerr)-exp(Kcoeff)
#        predicted = fitted(fit)
#        krt=summary(fit)$parameters[1,1]
#	return(list(K=krt,Kerr=Kerr,predictions=predicted, model=c("MMenten",paste("Krt =",round(krt,digits=3))),SSr=round(sum((residuals(fit))^2)),SAr=round(sum(abs(residuals(fit))))))
	if (length(p) != length(s)) {
		stop(print("Error: p and s not same length. Cannot fit Michaelis-Menten."))
	}
	#require("bbmle")
	LL <- function(krt,sigma) {
		R = p-(1-(s/((krt+s))))
		R = suppressWarnings(dnorm(R,0,sigma,log=TRUE))
		-sum(R)
	}
	LL_s <- function(krt,sigma) {
		p_nozero = p[p<1 & s > 0 & p>0];
		s_nozero = s[p<1 & s > 0 & p>0];
		R = log(s_nozero)/log(10)-log(krt*p_nozero/(1-p_nozero))/log(10); #horizontal residuals are log-normally distributed
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
	return(list(K=krt,Kerr=Kerr,fitted_err = res_err,predictions=predicted, model=c("MMenten",paste("Krt =",round(krt,digits=3))),SSr=round(sum((residuals)^2)),SAr=round(sum(abs(residuals)))))

}
bg__fit_logistic <- function(p,s) {
	if (length(p) != length(s)) {
		stop(print("Error: p and s not same length. Cannot fit Logistic Regression."))
	}
	s_nozero = s[s>0];
	p_nozero = p[s>0];
        logistic = suppressWarnings(glm(p_nozero~log(s_nozero),family="binomial")) #warns that is not binary data
        predlog = fitted(logistic)
	fullpredictions = rep(0, times=length(s));
	fullpredictions[s>0] = predlog
	res = fullpredictions-p;
	return(list(predictions=fullpredictions, B0 = logistic$coeff[1], B1=logistic$coeff[2] ,model=c( "Logistic", paste("Intercept =",round(logistic$coeff[1],digits=3)),paste("Coeff =",round(logistic$coeff[2],digits=3))),SSr=round(sum(res^2)),SAr=round(sum(abs(res)))));
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
	#require("bbmle")
#	LL <- function(lambda) {
#		R = p-exp(-lambda*s*s)
#		R = suppressWarnings(dnorm(R,0,0.1,log=TRUE))
#		-sum(R)
#	}
#	fit = mle2(LL,start=list(lambda=180/(mean(s)^2)))
	p_nozero = p; p_nozero[p == 0] = min(p[p>0])/10
	fit = lm(log(p_nozero)~-1+s^2)
	lambda = -fit$coeff[1]
	res_err = attributes(summary(fit))$coef[2,1]
	Lerr = summary(fit)$coeff[1,2]
	predicted = exp(-lambda*s*s)
	residuals = p-predicted
	return(list(lambda=lambda,Lerr=Lerr,fitted_err = res_err,predictions=predicted, model=c("p ~ e^(-lambda*S^2)",paste("lambda =",signif(lambda,digits=2))),SSr=round(sum((residuals)^2)),SAr=round(sum(abs(residuals)))))
}

M3Drop_Dropout_Models <- function(expr_mat, xlim=NA, suppress.plot=FALSE) {
	BasePlot = bg__dropout_plot_base(expr_mat, xlim = xlim, suppress.plot=suppress.plot);
	MM = bg__fit_MM(BasePlot$p, BasePlot$s);
	SCDE = bg__fit_logistic(BasePlot$p, BasePlot$s);
	ZIFA = bg__fit_ZIFA(BasePlot$p, BasePlot$s);
	if (!suppress.plot) {
	  	sizeloc = bg__add_model_to_plot(MM, BasePlot, lty=1, lwd=2.5, col="black",legend_loc = "topright");
		sizeloc = bg__add_model_to_plot(SCDE, BasePlot, lty=2, lwd=2.5, col="magenta3",legend_loc = c(sizeloc$rect$left+sizeloc$rect$w,sizeloc$rect$top-sizeloc$rect$h-0.05));
		sizeloc = bg__add_model_to_plot(ZIFA, BasePlot, lty=3, lwd=2.5, col="red",legend_loc = c(sizeloc$rect$left+sizeloc$rect$w,sizeloc$rect$top-sizeloc$rect$h-0.05));
	}
	invisible(list(MMFit = MM, LogiFit = SCDE, ExpoFit = ZIFA));
}

