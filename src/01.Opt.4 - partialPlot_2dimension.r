partialPlot_2dimension <-
    function (x, pred.data, x1.var, x2.var, which.class, w, plot=TRUE,
              n.pt = min(length(unique(pred.data[, x1.var])), 51),
              xlab=deparse(substitute(x1.var)), ylab=deparse(substitute(x2.var)), zlab="log(fraction of votes",
              main=paste("Partial Dependence on", deparse(substitute(x1.var)),"and", deparse(substitute(x2.var))),
              ...)
{
    classRF <- x$type != "regression"
    if (is.null(x$forest))
        stop("The randomForest object must contain the forest.\n")
    x1.var <- substitute(x1.var)
    x2.var <- substitute(x2.var)
    x1name <- if (is.character(x1.var)) x1.var else {
        if (is.name(x1.var)) deparse(x1.var) else {
            eval(x1.var)
        }
    }
    x2name <- if (is.character(x2.var)) x2.var else {
        if (is.name(x2.var)) deparse(x2.var) else {
            eval(x2.var)
        }
    }
    x1v <- pred.data[, x1name]
    x2v <- pred.data[, x2name]
    n <- nrow(pred.data)
    if (missing(w)) w <- rep(1, n)
    if (classRF) {
        if (missing(which.class)) {
            focus <- 1
        }
        else {
            focus <- charmatch(which.class, colnames(x$votes))
            if (is.na(focus))
                stop(which.class, "is not one of the class labels.")
        }
    	}
	
			if (is.ordered(x1v)) x1v <- as.numeric(x1v)
			if (is.ordered(x2v)) x2v <- as.numeric(x2v)

			x1.pt <- seq(min(x1v), max(x1v), length = n.pt)
			x2.pt <- seq(min(x2v), max(x2v), length = n.pt)
			y.pt <- matrix(NA, length(x1.pt), length(x2.pt))

			for (i in seq(along = x1.pt)) {
				for (j in seq(along = x2.pt)) {
					x.data <- pred.data
					x.data[, x1name] <- rep(x1.pt[i], n)
					x.data[, x2name] <- rep(x2.pt[j], n)
					if (classRF) {
							pr <- predict(x, x.data, type = "prob")
							y.pt[i,j] <- weighted.mean(log(ifelse(pr[, focus] == 0, .Machine$double.eps, pr[, focus])) - rowMeans(log(ifelse(pr == 0, .Machine$double.eps, pr))), w, na.rm=TRUE)
					} else {
							y.pt[i,j] <- weighted.mean(predict(x, x.data), w, na.rm=TRUE)
					}
				}
			}
			
			#if (plot) plot(x.pt, y.pt, type = "l", xlab=xlab, ylab=ylab,
			#									 main = main, ...)

			if (plot) persp(x=x1.pt, y=x2.pt, z=y.pt, phi = 45, theta = 45, xlab=xlab, ylab=ylab, zlab=zlab)

			return(list(x1.pt, x2.pt, y.pt))
}


#partial_output <- partialPlot_2dimension(x=popfit_final, pred.data=x_data, x1.var="dis", x2.var="YEARPOP", zlab="log(Pop. Density)", n.pt=20)
#
###	Using output re-call persp:
#persp(x=partial_output[[1]], y=partial_output[[2]], z=partial_output[[3]], phi = 30, theta = 30, xlab="dis", ylab="YEARPOP", zlab="log(Pop. Density)", main="Bi-variate Partial Plot", sub="Other Exp. Variables Held Constant at Median")
#
###	Or use rgl:
#require(rgl)
#ylim <- 10*range(partial_output[[3]])
#ylen <- ylim[2] - ylim[1] + 1
#colorlut <- terrain.colors(ylen) # height color lookup table
#col <- colorlut[ 10*partial_output[[3]]-ylim[1]+1 ] # assign colors to heights for each point
#rgl.surface(x=1:length(partial_output[[1]]), z=1:length(partial_output[[2]]), y=10*partial_output[[3]], color=col, back="lines")
#axes3d()
#title3d(xlab="dis", zlab="YEARPOP", ylab="log(Pop. Density)")
#
