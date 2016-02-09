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
