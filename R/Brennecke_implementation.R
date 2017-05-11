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

#This contains code by Brennecke et al. published in http://www.nature.com/nmeth/journal/v10/n11/full/nmeth.2645.html#supplementary-information
BrenneckeGetVariableGenes <- function(expr_mat, spikes=NA, suppress.plot=FALSE, fdr=0.1, minBiolDisp=0.5, fitMeanQuantile=0.8) {
        #require(statmod)

        rowVars <- function(x) { unlist(apply(x,1,var, na.rm=TRUE))}

        colGenes <- "black"
        colSp <- "blue"


        fullCountTable <- expr_mat;

        if (is.character(spikes)) {
                sp <- rownames(fullCountTable) %in% spikes;
                countsSp <- fullCountTable[sp,];
                countsGenes <- fullCountTable[!sp,];
        } else if (is.numeric(spikes)) {
                countsSp <- fullCountTable[spikes,];
                countsGenes <- fullCountTable[-spikes,];
        } else {
                countsSp <- fullCountTable;
                countsGenes <- fullCountTable;
        }

        meansSp <- rowMeans(countsSp, na.rm=TRUE)
        varsSp <- rowVars(countsSp)
        cv2Sp <- varsSp/meansSp^2
        meansGenes <- rowMeans(countsGenes, na.rm=TRUE)
        varsGenes <- rowVars(countsGenes)
        cv2Genes <- varsGenes/meansGenes^2
        # Fit Model
        minMeanForFit <- unname( quantile( meansSp[ which( cv2Sp > 0.3 ) ], fitMeanQuantile))
        useForFit <- meansSp >= minMeanForFit
        if (sum(useForFit, na.rm=TRUE) < 20) {
                warning("Too few spike-ins exceed minMeanForFit, recomputing using all genes.")
                meansAll <- c(meansGenes, meansSp)
                cv2All <- c(cv2Genes,cv2Sp)
                minMeanForFit <- unname( quantile( meansAll[ which( cv2All > 0.3 ) ], 0.80))
                useForFit <- meansSp >= minMeanForFit
        }
        if (sum(useForFit, na.rm=TRUE) < 30) {warning(paste("Only", sum(useForFit), "spike-ins to be used in fitting, may result in poor fit."))}
        fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/meansSp[useForFit] ), cv2Sp[useForFit] )
        a0 <- unname( fit$coefficients["a0"] )
        a1 <- unname( fit$coefficients["a1tilde"])
	res <- cv2Genes - (a0 + a1/meansGenes);

        # Test
        psia1theta <- a1
        minBiolDisp <- minBiolDisp^2
        m <- ncol(countsSp);
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
		# Plot/highlight the spike-ins if they are different from the genes
		if (length(meansSp) < length(meansGenes)) {
			points(meansSp, cv2Sp, pch=20, cex=.5, col=colSp)
		}
                # Add the technical noise fit
                xg <- 10^seq( -2, 6, length.out=1000 )
                lines( xg, (a1)/xg + a0, col="#FF000080", lwd=3 )
                # Add a curve showing the expectation for the chosen biological CV^2 thershold
                lines( xg, psia1theta/xg + a0 + minBiolDisp, lty="dashed", col="#C0007090", lwd=3)
        }
	TABLE <- data.frame(Gene = names(meansGenes)[sig], effect.size=res[sig], p.value = p[sig], q.value= padj[sig])
	TABLE <- TABLE[order(-TABLE[,2]),];
        return(TABLE)
}
