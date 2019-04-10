NBumiHVG <- function(counts, fit, fdr_thresh=0.05, suppress.plot=FALSE, method=c("DANB", "basic")) {
	# v of v : https://math.stackexchange.com/questions/72975/variance-of-sample-variance
	# Moments of NB: http://mathworld.wolfram.com/NegativeBinomialDistribution.html

	n <- ncol(counts);
	
	# DANB model
	if (method[1] == "DANB") {
		#fit <- NBumiFitModel(counts);
		mu_obs <- fit$vals$tjs/n
		v_obs <- mu_obs + mu_obs^2/fit$size
	} else {
		mu_obs <- Matrix::rowMeans(counts)
		v_obs <- Matrix::rowSums((counts-mu_obs)^2)/(n-1)
	}
	
	# var = mu + disp * mu^2

	tmp <- mu_obs^2
	disp <- glm((v_obs-mu_obs)~tmp+0)$coef[1]
	v_fitted <- mu_obs+disp*mu_obs^2

	p <- mu_obs/v_fitted
	r <- mu_obs*p/(1-p)	


	mu4 <- r*(1-p)*(6-6*p+p^2+3*r-3*p*r)/(p^4)
	sigma2 <- r*(1-p)/(p^2)

	#v_of_v <- mu4/n - (sigma2^2*(n-3)/(n*(n-1)) #https://math.stackexchange.com/questions/72975/variance-of-sample-variance

	v_of_v <- mu4*(n-1)^2/n^3 - (sigma2^2*(n-3)*(n-1))/(n^3) #http://mathworld.wolfram.com/SampleVarianceDistribution.html

	z <- (v_obs - sigma2)/sqrt(v_of_v)
	p <- pnorm(z, lower.tail=FALSE)
	q <- p.adjust(p, method="fdr")
        eff <- v_obs-sigma2
	tab <- data.frame(rownames(counts), eff, p, q)
	colnames(tab) <- c("Gene", "effect.size", "p.value", "q.value")
	tab <- tab[!is.na(tab$p.value),]
	tab <- tab[order(-tab$q.value, tab$effect.size, decreasing=TRUE),]
	if (!suppress.plot) {
		plot(mu_obs, v_obs, cex=0.75, pch=16, xlab="mean", ylab="variance",log="xy")
		points(mu_obs[q < fdr_thresh], v_obs[q < fdr_thresh], col="red", pch=16)
		# Lines
		reorder <- order(mu_obs)
		lines(mu_obs[reorder], sigma2[reorder], col="grey80", lwd=2, lty=1)
		lines(mu_obs[reorder], sigma2[reorder]+sqrt(v_of_v[reorder])*qnorm(fdr_thresh, lower.tail=FALSE), col="grey80", lwd=2, lty=2)
	}
	return(tab[tab$q.value < fdr_thresh,])
}
