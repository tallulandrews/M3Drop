# Normalization Functions
hidden__UQ <- function(x){quantile(x[x>0],0.75)};

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

hidden__normalize <- function(data) {
	# Combine UQ and detection rate adjusted normalization 
	# Stephanie Hick, Mingziang Teng, Rafael A Irizarry "On the widespread and critical impact of systematic single-cell RNA-Seq data" http://dx.doi.org/10.1101/025528 
	cell_zero = colSums(data == 0)/length(data[,1]);
	uq = unlist(apply(data,2,hidden__UQ));
	normfactor = (uq/median(uq)) * (median(cell_zero)/cell_zero); 
	data = t(t(data)/normfactor);
	return(data);
}

M3D_Clean_Data <- function(expr_mat, labels = NA, is.counts=TRUE, suppress.plot=FALSE, pseudo_genes=NA, min_detected_genes=NA) {
	if (length(pseudo_genes) > 1) {
		is_pseudo = rownames(expr_mat) %in% as.character(pseudo_genes);
	        expr_mat = expr_mat[!is_pseudo,];
	}

	data_list = bg__filter_cells(expr_mat, labels, suppress.plot = suppress.plot, threshold=min_detected_genes);
	
        detected = rowSums(data_list$expr_mat > 0) > 3;
        expr_mat = data_list$expr_mat[detected,];
	labels   = data_list$labels

	spikes = grep("ercc",rownames(expr_mat), ignore.case=TRUE)
	if (is.counts) {
		if (length(spikes) > 1) {
	                totreads = colSums(expr_mat[-c(spikes),])
		} else {
			totreads = colSums(expr_mat);
		}
                cpm = t(t(expr_mat)/totreads)*1000000;
                lowExpr = rowMeans(cpm) < 10^-5;
                cpm=cpm[!lowExpr,];
                return(list(data=cpm, labels=labels));
        }

	lowExpr = rowMeans(expr_mat) < 10^-5;
        data=expr_mat[!lowExpr,];
        return(list(data=expr_mat, labels=labels));
}
