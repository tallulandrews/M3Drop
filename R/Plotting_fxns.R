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

# Modularize this stuff more sensibly
#  Plotting Functions
bg__dropout_plot_base <- function (expr_mat, xlim = NA, suppress.plot=FALSE) {
	#require("RColorBrewer")
	
	gene_info <- bg__calc_variables(expr_mat);

        xes <- log(gene_info$s)/log(10);
        put_in_order <- order(xes);
        fancy <- densCols(xes, gene_info$p, colramp=colorRampPalette(c("black","white")))
        dens <- col2rgb(fancy)[1,]+1L
#        colours <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F",
#                                    "#FCFF00", "#FF9400", "#FF3100"))(256) #rainbow
        colours <-  colorRampPalette(c("#000099", "#00FEFF", "#FCFF00"))(256) #blue->yellow
        dens.col <- colours[dens]

	if (!suppress.plot) {
        	par(fg="black")
		if (!(sum(is.na(xlim)))) {
	        	plot(xes,gene_info$p, main="", ylab="", xlab="", col = dens.col,pch=16, xlim=xlim, ylim=c(0,1))
		} else {
	        	plot(xes,gene_info$p, main="", ylab="", xlab="", col = dens.col,pch=16, ylim=c(0,1))
		}
		title(ylab="Dropout Rate", line=2)
		title(xlab="log10(expression)", line=2)
	}
	invisible(list(gene_info = gene_info, xes=xes, order=put_in_order));
}

bg__add_model_to_plot <- function(fitted_model, base_plot, lty=1, lwd=1, col="black",legend_loc = "topright") {
	lines(base_plot$xes[base_plot$order],fitted_model$predictions[base_plot$order],lty=lty,lwd=lwd,col=col);
	par(fg=col)
	if (length(legend_loc) == 2) {
        	this_loc <- legend(legend_loc[1], legend_loc[2], fitted_model$model, box.lty=lty, box.lwd=lwd, xjust=1)
	} else {
		this_loc <- legend(legend_loc[1], fitted_model$model, box.lty=lty, box.lwd=lwd, xjust=1)
	}
	par(fg="black")
	invisible(this_loc)
}

bg__highlight_genes <- function (base_plot, data, genes, colour="purple", pch=16) {
	if(!is.numeric(genes) && !is.logical(genes)) {
		genes <- match(as.character(genes), rownames(data));
		nomatch <- sum(is.na(genes));
		if (nomatch > 0) {warning(paste(nomatch, " genes could not be matched to data, they will not be highlighted."));}
		if (nomatch == length(genes)) {invisible(cbind(c(NA,NA),c(NA,NA)))}
		genes <- genes[!is.na(genes)];
	}
	points(base_plot$xes[genes],base_plot$gene_info$p[genes],col=colour, pch=pch)
	invisible(cbind(base_plot$gene_info$s[genes],base_plot$gene_info$p[genes]));
}

bg__expression_heatmap <- function (genes, expr_mat, cell_labels=NA, gene_labels=NA, key_genes=genes, key_cells=NA) { 
	#require("RColorBrewer")
	#require("gplots")
	if(!is.numeric(genes)) {
		new_genes <- match(genes, rownames(expr_mat));
		nomatch <- sum(is.na(new_genes));
		if (nomatch > 0) {warning(paste("Warning: ",nomatch, " genes could not be matched to data, they will not be included in the heatmap."));}
		genes <- new_genes[!is.na(new_genes)];
	}
	if (length(genes) < 1) {stop("Error: No genes for heatmap.");return();}
	# Plot heatmap of expression
	heatcolours <- rev(brewer.pal(11,"RdBu"))
	col_breaks <- c(-100,seq(-2,2,length=10),100)
	heat_data <- as.matrix(expr_mat[genes,])
	heat_data <- log(heat_data+1)/log(2);
	ColColors <- rep("white", times=length(heat_data[1,]))
	RowColors <- rep("white", times=length(heat_data[,1]))
	# remove row & column labels
	rownames(heat_data) <- rep("", length(heat_data[,1]));
	if (!is.na(key_genes[1])) {
		rownames(heat_data)[rownames(expr_mat[genes,]) %in% key_genes] <- rownames(expr_mat[genes,])[rownames(expr_mat[genes,]) %in% key_genes]; 
	}
	if(length(unique(colnames(heat_data))) < length(heat_data[1,])) {
		colnames(heat_data) <- 1:length(colnames(heat_data));
	}
	if (!is.na(key_cells[1])) {
		colnames(heat_data)[colnames(expr_mat[genes,]) %in% key_cells] <- colnames(expr_mat[genes,])[colnames(expr_mat[genes,]) %in% key_cells]; 
	}
	if (!is.na(cell_labels[1])) {
		cell_labels <- as.character(cell_labels);
		colours <- as.factor(cell_labels)
		palette <- brewer.pal(max(3,length(unique(cell_labels))), "Set3");
		ColColors <- palette[colours];	
		mylegend<- list(names = unique(cell_labels), fill = unique(ColColors));
	} 
	if (!is.na(gene_labels[1])) {
		# lowest factor level = grey (so 0-1 is striking)
		if (!is.numeric(gene_labels)) {
			colours <- as.factor(gene_labels)
		} else {
			colours <- gene_labels
		}
		palette <- c("grey75",brewer.pal(max(3,length(unique(gene_labels))), "Set1"));
		RowColors <- palette[colours];
	}
	# Custom Shit
	lwid<-c(1,0.2,4)
	lhei<-c(1,0.2,4)
	lmat<-rbind(c(6,0,5),c(0,0,2),c(4,1,3))


	if (dim(heat_data)[1] < 10000) {
		heatmap_output <- suppressWarnings(heatmap.2(heat_data, ColSideColors = ColColors, RowSideColors = RowColors, col=heatcolours, breaks=col_breaks, scale="row",symbreaks=TRUE, trace="none", dendrogram="column", key=FALSE, Rowv=TRUE, Colv=TRUE,lwid=lwid, lhei=lhei,lmat=lmat, hclustfun=function(x){hclust(x,method="ward.D2")}))
	} else {
		heatmap_output = suppressWarnings(heatmap.2(heat_data, ColSideColors = ColColors, RowSideColors = RowColors, col=heatcolours, breaks=col_breaks, scale="row",symbreaks=TRUE, trace="none", dendrogram="column", key=FALSE, Rowv=FALSE, Colv=TRUE,lwid=lwid, lhei=lhei,lmat=lmat, hclustfun=function(x){hclust(x,method="ward.D2")}))
	}
	# Custom key
	par(fig = c(0, 1/(5.2),4/(5.2), 1), mar=c(4,1,1,1), new=TRUE)
	scale01 <- function(x, low = min(x), high = max(x)) {
        	x <- (x - low)/(high - low)
        	x
    	}
	par(mar=c(5,1,1,1))
	par(cex=0.75)
	par(mgp=c(2,1,0))
	key_breaks <- seq(-2,2,length=10)
	key_col <- heatcolours[2:(length(heatcolours)-1)]
	z <- seq(min(key_breaks),max(key_breaks), by=min(diff(key_breaks)/4))
	image(z=matrix(z,ncol=1),col=key_col,breaks=key_breaks,xaxt="n",yaxt="n")
	par(usr = c(0, 1, 0, 1))
	lv <- pretty(key_breaks)
        xv <- scale01(as.numeric(lv), min(key_breaks),max(key_breaks))
        xargs <- list(at = xv, labels = lv)
	xargs$side <- 1
	do.call(axis, xargs)
	mtext(side = 1, "Expression Z-Score", line = par("mgp")[1], padj = 0.5, 
                cex <- par("cex") * par("cex.lab"))

	# Legend
	par(fig = c(0/5.2, 1/(5.2),0/(5.2), 4/5.2), mar=c(0,0,0,0), new=TRUE)
	par(mar=c(0,0,0,0))
	if (!is.na(cell_labels[1])) {
		legend("left", mylegend$names, pt.bg = mylegend$fill,bg="white",col="black", pch=22, pt.cex=2.5, cex=1.25, bty="n",y.intersp = 2);
	}
	invisible(heatmap_output);
}

M3DropExpressionHeatmap <- function(genes, expr_mat, cell_labels=NA, interesting_genes=NA, key_genes=genes, key_cells=NA) {
	# Converted known DE genes into heatmap labels 
	gene_labels <- rep(1, times = length(genes));
	if (is.na(interesting_genes[1])) {
		gene_labels<-NA
	}
 	if (is.list(interesting_genes)) {
                for (i in 1:length(interesting_genes)) {
                        gene_labels[genes %in% interesting_genes[[i]]] <- i+1;
                }
        } else {
                gene_labels[genes %in% interesting_genes] <- 2;
        }
	if (is.numeric(key_genes) | is.logical(key_genes)) {
		key_genes <- rownames(expr_mat)[key_genes];
	}
	if (is.numeric(key_cells) | is.logical(key_cells)) {
		key_cells <- rownames(expr_mat)[key_cells];
	}
	if (is.factor(genes)) {
		genes <- as.character(genes);
	}
	if (!is.vector(genes)) {
		is.gene <- grepl("gene",colnames(genes), ignore.case=TRUE)
		if (sum(is.gene) == 1) {
			genes <- unlist(genes[,is.gene]);
		} else {
			stop("Error: please provide a vector of gene names not a table.")
		}
	}
	heatmap_output <- bg__expression_heatmap(genes, expr_mat, cell_labels=cell_labels, gene_labels=as.numeric(gene_labels), key_genes=as.character(key_genes), key_cells=key_cells);
	invisible(heatmap_output);
}

M3DropGetHeatmapCellClusters <- function (heatout, k) {
        dendro<-heatout$colDendrogram
        curr_k <- 1;
        dendro_list <- list(dendro)
        dendro_heights <- attr(dendro, "height")
        while( curr_k < k ){
                to_split <- which(dendro_heights == max(dendro_heights))
                to_split_dendro <- dendro_list[[to_split]]
                to_split_height <-  dendro_heights[to_split]

                children <- as.list(to_split_dendro)
                for (i in 1:length(children)) {
                        dendro_heights <- c(dendro_heights,attr(children[[i]],"height"))
                        dendro_list[[length(dendro_list)+1]] <- children[[i]]
                }
                # Remove to split
                dendro_list[to_split] <- NULL
                dendro_heights <- dendro_heights[-to_split]
                curr_k <- curr_k-1+length(children)
        }
        # Make group vector
        names_orig_order <- labels(dendro)[order(heatout$colInd)]
        groups <- rep(0, times=length(names_orig_order))
        for (i in 1:length(dendro_list)) {
                groups[names_orig_order %in% labels(dendro_list[[i]])] <- i
        }
	names(groups) <- names_orig_order;
        return(groups);
}
