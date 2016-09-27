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

# By: Tallulah Andrews, Senior Ph.D Student, 
# MRC Functional Genomics Unit, University of Oxford.
#
# Using three steps takes three sets and creates a labelled three-way Venn Diagram with areas approximately proportional to the counts.
# steps:
# table = ConvertSetstoTable(set1, set2, set3, names)
# circles = plot.venn.diagram(table)
# PlaceLabels(circles, table, 1/0) -> last argument indicates whether to place set names within their circle

#Based on Original Code from: http://tolstoy.newcastle.edu.au/R/help/03a/1115.html
# By: David J. States, M.D., PhD. Professor of Human Genetics, Director of Bioinformatics, University of Mishigan School of Medicine, Ann Arbor, MI, USA
# Date Downloaded: 10 Apr 2014
# Areas of circles and areas of pairwise overlaps proportional to counts, hopes for the best for the 3-way overlap

venn.overlap <-
function(r, a, b, target = 0)
{
#
# calculate the overlap area for circles of radius a and b
# with centers separated by r
# target is included for the root finding code
#
        pi <- acos(-1)
        if(r >= a + b) {
                return( - target)
        }
        if(r < a - b) { #b completely overlapped by a
                return(pi * b * b - target)
        }
        if(r < b - a) { #a complely overlapped by b
                return(pi * a * a - target)
        }
        s <- (a + b + r)/2
        triangle.area <- sqrt(s * (s - a) * (s - b) * (s - r))
        h <- (2 * triangle.area)/r
        aa <- 2 * atan(sqrt(((s - r) * (s - a))/(s * (s - b))))
        ab <- 2 * atan(sqrt(((s - r) * (s - b))/(s * (s - a))))
        sector.area <- aa * (a * a) + ab * (b * b)
        overlap <- sector.area - 2 * triangle.area
        return(overlap - target)
} 

plot.venn.diagram <-
function(d)
{
#
# Draw Venn diagrams with proportional overlaps
# d$table = 3 way table of overlaps
# d$labels = array of character string to use as labels
#
pi <- acos(-1)
csz <- 0.1
# Normalize the data
n <- length(dim(d$table))
c1 <- vector(length = n)
c1[1] <- sum(d$table[2, , ])
c1[2] <- sum(d$table[, 2, ])
c1[3] <- sum(d$table[, , 2])
n1 <- c1
#
c2 <- matrix(nrow = n, ncol = n, 0)
c2[1, 2] <- sum(d$table[2, 2, ])
c2[2, 1] <- c2[1, 2]
c2[1, 3] <- sum(d$table[2, , 2])
c2[3, 1] <- c2[1, 3]
c2[2, 3] <- sum(d$table[, 2, 2])
c2[3, 2] <- c2[2, 3]
n2 <- c2
#
c3 <- d$table[2, 2, 2]
n3 <- c3
c2 <- c2/sum(c1)
c3 <- c3/sum(c1)
c1 <- c1/sum(c1)
n <- length(c1)
# Radii are set so the area is proporitional to number of counts
pi <- acos(-1)
r <- sqrt(c1/pi)
r12 <- uniroot(venn.overlap, interval = c(max(r[1] - r[2], r[2] - r[1],
0) + 0.01, r[1] + r[2] - 0.01), a = r[1], b = r[
2], target = c2[1, 2])$root
r13 <- uniroot(venn.overlap, interval = c(max(r[1] - r[3], r[3] - r[1],
0) + 0.01, r[1] + r[3] - 0.01), a = r[1], b = r[
3], target = c2[1, 3])$root
r23 <- uniroot(venn.overlap, interval = c(max(r[2] - r[3], r[3] - r[2],
0) + 0.01, r[2] + r[3] - 0.01), a = r[2], b = r[
3], target = c2[2, 3])$root
s <- (r12 + r13 + r23)/2
x <- vector()
y <- vector()
x[1] <- 0
y[1] <- 0
x[2] <- r12
y[2] <- 0
angle <- 2 * atan(sqrt(((s - r12) * (s - r13))/(s * (s - r13))))
x[3] <- r13 * cos(angle)
y[3] <- r13 * sin(angle)
xc <- cos(seq(from = 0, to = 2 * pi, by = 0.01)) #Resolution of Circles
yc <- sin(seq(from = 0, to = 2 * pi, by = 0.01)) #Resolution of Circles
cmx <- sum(x * c1)
cmy <- sum(y * c1)
x <- x - cmx
y <- y - cmy
rp<-sqrt(x*x + y*y)
frame()
par(usr = c(-1, 1, -1, 1), pty = "s")
box()
for(i in 1:3) { #Draw Circles
lines(xc * r[i] + x[i], yc * r[i] + y[i])
}
return(list(r = r, x = x, y = y))
} 

# Calculate intersection to robustly place labels in appropriate locations
circle_intersection <- function (x1,y1,r1,x2,y2,r2) {
	if (x1 == x2) {
		y <- (r2*r2-r1*r1+x1*x1-x2*x2+y1*y1-y2*y2)/(2*y1-2*y2);
		a <- 1; b <- 2*x1; c <- (y-y1)*(y-y1)-r1*r1+x1*x1;
		xi1 <- (-b+sqrt(b*b-4*a*c))/(2*a)
                xi2 <- (-b-sqrt(b*b-4*a*c))/(2*a)
 		return(rbind(c(xi1,y),c(xi2,y)));
	} else {
		A <- (r2*r2-r1*r1+x1*x1-x2*x2+y1*y1-y2*y2)/(2*x1-2*x2)
		B <- (2*y2-2*y1)/(2*x1-2*x2)

		a <- (B*B)+1; b <- -2*B*x1+2*A*B-2*y1; c <- A*A-2*x1*A+x1*x1+y1*y1-r1*r1;
		yi1 <- (-b+sqrt(b*b-4*a*c))/(2*a)
		yi2 <- (-b-sqrt(b*b-4*a*c))/(2*a)
		xi1 <- A+B*yi1
		xi2 <- A+B*yi2
		return(rbind(c(xi1,yi1),c(xi2,yi2)));
	}
}


circle_line_intersection<-function (xc,yc,r,x1,y1,x2,y2) {
	if (x2 == x1) {
		xi <- x1;
		a <- 1; b <- -2*yc; c <- (xi-xc)*(xi-xc)+yc*yc-r*r;
	        yi1 <- (-b+sqrt(b*b-4*a*c))/(2*a)
	        yi2 <- (-b-sqrt(b*b-4*a*c))/(2*a)
		return(rbind(c(xi, yi1), c(xi, yi2)))
	} else {
		M <- (y2-y1)/(x2-x1);
		B <- y1-M*x1;
		a <- M*M+1; b <- 2*M*(B-yc)-2*xc; c <- xc*xc+(B-yc)*(B-yc)-r*r;
	        xi1 <- (-b+sqrt(b*b-4*a*c))/(2*a)
	        xi2 <- (-b-sqrt(b*b-4*a*c))/(2*a)
		yi1 <- M*xi1+B
		yi2 <- M*xi2+B
		return(rbind(c(xi1,yi1),c(xi2,yi2)));
	}
}


#By: Tallulah Andrews
#Date: 4th April 2014
# Places labels in each section of the venn diagram based on the interection points of the three circles

PlaceLabels <- function(circles, d, labelsincircles) {
  circles$left <- circles$x-circles$r
  circles$right <- circles$x+circles$r
  circles$top <- circles$y+circles$r
  circles$bottom <- circles$y-circles$r


  # Better labelling - Does not cope with one circle completely within another circle
  topcircle <- which(circles$y == max(circles$y))
  rightcircle <- which(circles$x == max(circles$x))
  leftcircle <- which(circles$x == min(circles$x))

  # Names for each Circle
  # Outside most distant edge of the circle
  if (!labelsincircles) {
	  text(circles$x[topcircle],circles$top[topcircle], d$labels[topcircle],pos=3)
	  text(circles$right[rightcircle],circles$y[rightcircle], d$labels[rightcircle],pos=4)
	  text(circles$left[leftcircle],circles$y[leftcircle], d$labels[leftcircle],pos=2)
  }

  # Intersection Points
  intersectionLR <- circle_intersection(circles$x[leftcircle],circles$y[leftcircle],circles$r[leftcircle],circles$x[rightcircle],circles$y[rightcircle],circles$r[rightcircle])
  intersectionLT <- circle_intersection(circles$x[leftcircle],circles$y[leftcircle],circles$r[leftcircle],circles$x[topcircle],circles$y[topcircle],circles$r[topcircle])
  intersectionTR <- circle_intersection(circles$x[topcircle],circles$y[topcircle],circles$r[topcircle],circles$x[rightcircle],circles$y[rightcircle],circles$r[rightcircle])


  onewayvalues <- c(d$table[2, 1, 1], d$table[1, 2, 1], d$table[1, 1, 2]) #1, 2, 3
  twowayvalues <- c(d$table[2,2,1],d$table[2,1,2],d$table[1,2,2]) #1&2, 1&3, 2&3 ->ceiling(a*b/2)
  threewayvalue <- d$table[2,2,2];

  # Names in circles above their one-way value
  if (labelsincircles) {onewayvalues <- paste(d$labels, onewayvalues, sep="\n");}

  # 3-way intersection
  topLRintersection <- intersectionLR[which(intersectionLR[,2] == max(intersectionLR[,2])),];
  leftTRintersection <- intersectionTR[which(intersectionTR[,1] == min(intersectionTR[,1])),];
  rightTLintersection <- intersectionLT[which(intersectionLT[,1] == max(intersectionLT[,1])),];
  boundarypoints <- rbind(topLRintersection, leftTRintersection, rightTLintersection) #R of LT intersection
  midpoint <- colMeans(boundarypoints); text(midpoint[1],midpoint[2], threewayvalue)


  #2-way intersection
  # LR - take advantage of fact that L & R circles has centres at the same height/
  textx <- mean(intersectionLR[,1]); # mid-point betweed the intersection points
  lowestT <- min(leftTRintersection[2],rightTLintersection[2], topLRintersection[2]);
  if (circles$x[topcircle] > rightTLintersection[1] & circles$x[topcircle] < leftTRintersection[1]) {lowestT <- min(circles$bottom[topcircle], lowestT);}
  texty <- mean(c(min(intersectionLR[,2]), lowestT))
  text(textx, texty, twowayvalues[ceiling(leftcircle*rightcircle/2)])

  # LT
  intersectpoint <- circle_line_intersection(circles$x[rightcircle], circles$y[rightcircle], circles$r[rightcircle], rightTLintersection[1], rightTLintersection[2], intersectionLT[which(intersectionLT[,1] != max(intersectionLT[,1])),1], intersectionLT[which(intersectionLT[,1] != max(intersectionLT[,1])),2])

# if intersect points are outside of the area wedge use the inner circle-pairwise intersection instead
  midpointv <- rbind(intersectionLT[which(intersectionLT[,1] == min(intersectionLT[,1])),], intersectpoint[which(intersectpoint[,1] == min(intersectpoint[,1])),]);
  if (midpointv[2,1] > rightTLintersection[1]) {midpointv[2,] <- rightTLintersection;}
  midpoint <- colMeans(midpointv); 
  text(midpoint[1], midpoint[2], twowayvalues[ceiling(leftcircle*topcircle/2)])

  # RT
  intersectpoint <- circle_line_intersection(circles$x[leftcircle], circles$y[leftcircle], circles$r[leftcircle], leftTRintersection[1], leftTRintersection[2], intersectionTR[which(intersectionTR[,1] !=  min(intersectionTR[,1])),1], intersectionTR[which(intersectionTR[,1] !=  min(intersectionTR[,1])),2])
  midpointv <- rbind(intersectionTR[which(intersectionTR[,1] !=  min(intersectionTR[,1])),], intersectpoint[which(intersectpoint[,1] == max(intersectpoint[,1])),]);
  if (midpointv[2,1] < leftTRintersection[1]) {midpointv[2,] <- leftTRintersection;}
  midpoint <- colMeans(midpointv);
  text(midpoint[1], midpoint[2], twowayvalues[ceiling(rightcircle*topcircle/2)])
  #1-way intersection
  #L 
#if distance from  leftTRintersection to the centre is smaller than the radius. else put in centre
  if (dist(rbind(leftTRintersection,c(circles$x[leftcircle], circles$y[leftcircle])), method="euclidean") < circles$r[leftcircle]) {
	  intersectpoint <- circle_line_intersection(circles$x[leftcircle], circles$y[leftcircle], circles$r[leftcircle], intersectionTR[which(intersectionTR[,1] != leftTRintersection[1]),1], intersectionTR[which(intersectionTR[,1] != leftTRintersection[1]),2], leftTRintersection[1], leftTRintersection[2])
	  midpoint <- colMeans(rbind(leftTRintersection, intersectpoint[which(intersectpoint[,1] == min(intersectpoint[,1])),]))
	  text(midpoint[1], midpoint[2], onewayvalues[leftcircle])
  } else {
	  text(circles$x[leftcircle], circles$y[leftcircle],onewayvalues[leftcircle]);
  }
  #R 
#if distance from  rightTLintersection to the centre is smaller than the radius. else put in centre
  if (dist(rbind(rightTLintersection,c(circles$x[rightcircle], circles$y[rightcircle])), method="euclidean") < circles$r[rightcircle]) {
	  intersectpoint <- circle_line_intersection(circles$x[rightcircle], circles$y[rightcircle], circles$r[rightcircle], intersectionLT[which(intersectionLT[,1] != rightTLintersection[1]),1], intersectionLT[which(intersectionLT[,1] != rightTLintersection[1]),2] , rightTLintersection[1], rightTLintersection[2])
	  midpoint <- colMeans(rbind(rightTLintersection, intersectpoint[which(intersectpoint[,1] == max(intersectpoint[,1])),]))
	  text(midpoint[1], midpoint[2], onewayvalues[rightcircle])
  } else {
	text(circles$x[rightcircle], circles$y[rightcircle],onewayvalues[rightcircle]);
  }
  #T
#if distance from  topLRintersection to the centre is smaller than the radius. else put in centre
  if (dist(rbind(topLRintersection,c(circles$x[topcircle], circles$y[topcircle])), method="euclidean") < circles$r[topcircle]) {
  	intersectpoint <- circle_line_intersection(circles$x[topcircle], circles$y[topcircle], circles$r[topcircle], intersectionLR[which(intersectionLR[,2] != topLRintersection[2]),1], intersectionLR[which(intersectionLR[,2] != topLRintersection[2]),2], topLRintersection[1], topLRintersection[2])
  	midpoint <- colMeans(rbind(topLRintersection, intersectpoint[which(intersectpoint[,2] == max(intersectpoint[,2])),]))
  	text(midpoint[1], midpoint[2], onewayvalues[topcircle])
  } else {
	text(circles$x[topcircle], circles$y[topcircle],onewayvalues[topcircle]);
  }
}



# Based on code from: http://research.stowers-institute.org/efg/R/Math/VennDiagram.htm

ConvertSetstoTable <- function(set1, set2, set3, names)
{stopifnot( length(names) == 3)
  # Form universe as union of all three sets
  universe <- sort( unique( c(set1, set2, set3) ) )
  InOut <- matrix(0, nrow=length(universe), ncol=3) #each element in/out of each set
  colnames(InOut) <- names

  for (i in 1:length(universe)) {
    InOut[i,1] <- universe[i] %in% set1
    InOut[i,2] <- universe[i] %in% set2
    InOut[i,3] <- universe[i] %in% set3
  }
  d <-list()
  d$table<-table(InOut[,1],InOut[,2], InOut[,3])
  d$labels <- names;
  return(d)
}

# Top-level function
M3DropThreeSetVenn <- function(set1,set2,set3,names){
  table <- ConvertSetstoTable(as.character(set1), as.character(set2), as.character(set3), names=names)
  nicer_table <- table
  if (max(table$table) > 5*mean(table$table)) {
    bigones <- table$table > 5*mean(table$table);
    nicer_table$table[bigones] <- table$table[bigones]/5
  }
  if (sum(table$table == 0) > 1) {
    smallest_notzero <- min(table$table[table$table > 0])
    nicer_table$table[table$table == 0] <- smallest_notzero;
  }
  circles <- plot.venn.diagram(nicer_table);
  PlaceLabels(circles, table, labelsincircles=TRUE)
}
