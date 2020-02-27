MiHC.stat <-
function(Zs, hs=hs, Ws=Ws) {

	d <- length(Zs)
	pvs <- 1-pchisq(Zs^2, df=1)
	i0 <- which(pvs < 0.00000001)
	i1 <- which(pvs > 0.99999999)
	pvs[i0] <- 0.00000001
	pvs[i1] <- 0.99999999
	Is <- order(order(pvs))
	HCs <- (Is/d-pvs)/sqrt(pvs*(1-pvs)/d)
	HCs.ord <- HCs[order(HCs, decreasing=TRUE)]
	if (!is.null(Ws)) {
		wHCs <- Ws*(Is/d-pvs)/sqrt(pvs*(1-pvs)/d)
		wHCs.ord <- wHCs[order(wHCs, decreasing=TRUE)]
	}
	if (is.null(Ws)) {
		maxHC <- sapply(as.list(hs), function(x) sum(HCs.ord[1:x], na.rm=TRUE))
	}
	if (!is.null(Ws)) {
		maxHC.u <- sapply(as.list(hs), function(x) sum(HCs.ord[1:x], na.rm=TRUE))
		maxHC.w <- sapply(as.list(hs), function(x) sum(wHCs.ord[1:x], na.rm=TRUE))
		maxHC <- c(maxHC.u, maxHC.w)	
	}
	return(maxHC)
}
