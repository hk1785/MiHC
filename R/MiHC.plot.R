MiHC.plot <-
function(MiHC.out, leg.loc="bottomright", pdf.filename=NULL) {

	if (length(MiHC.out$graph.els) == 5) {
		otu.ids <- MiHC.out$graph.els$otu.ids
		exp.HC <- MiHC.out$graph.els$exp.HC
		obs.HC <- MiHC.out$graph.els$obs.HC
		exp.wHC <- MiHC.out$graph.els$exp.wHC
		obs.wHC <- MiHC.out$graph.els$obs.wHC		
		HCs <- abs(exp.HC-obs.HC)
		wHCs <- abs(exp.wHC-obs.wHC)
		d <- length(otu.ids)
		if (!is.null(pdf.filename)) {
			pdf(file=pdf.filename, width=4.5, height=9.5)
			layout(matrix(c(1,2), nrow=2, byrow=TRUE))
			par(mai=c(1,0.8,0.4,0.2))
			par(oma=c(0,0,1,0))

			uHC.pvs <- round(c(MiHC.out$simes.pv, MiHC.out$ind.pvs[,1], MiHC.out$ada.pvs[1]),3)
			ind.0 <- which(uHC.pvs == 0)
			ind.1 <- which(uHC.pvs == 1)
			uHC.pvs[ind.0] <- "<.001"
			uHC.pvs[ind.1] <- ">.999"
			if (paste(rownames(MiHC.out$ind.pvs), collapse="") == c("13579")) {
				xlab.v <- paste("\n\nObserved quantiles\nuHC(1):", uHC.pvs[2], ", uHC(3):", uHC.pvs[3], ", uHC(5):", uHC.pvs[4], ", uHC(7):", uHC.pvs[5], ", uHC(9):", uHC.pvs[6], "\nSimes:", uHC.pvs[1], ", uHC(A):", uHC.pvs[7], sep="")
			} else {
				xlab.v <- paste("uHC(A):", uHC.pvs[length(uHC.pvs)], sep="")
			}
			plot(obs.HC, exp.HC, pch=20, panel.first=grid(lty=2, lwd=1), xlab=xlab.v, ylab="Expected quantiles", main="uHC", cex.lab=0.7, cex.axis=0.7, cex.main=0.9, col="blue2")
	
			abline(0, 1, col="red2")
			abline(v=obs.HC[which(rank(HCs)==d)], col="gray1")
			abline(v=obs.HC[which(rank(HCs)==d-1)], col="gray10")
			abline(v=obs.HC[which(rank(HCs)==d-2)], col="gray20")
			abline(v=obs.HC[which(rank(HCs)==d-3)], col="gray30")
			abline(v=obs.HC[which(rank(HCs)==d-4)], col="gray40")
			abline(v=obs.HC[which(rank(HCs)==d-5)], col="gray50")
			abline(v=obs.HC[which(rank(HCs)==d-6)], col="gray60")
			abline(v=obs.HC[which(rank(HCs)==d-7)], col="gray70")
			abline(v=obs.HC[which(rank(HCs)==d-8)], col="gray80")
			abline(v=obs.HC[which(rank(HCs)==d-9)], col="gray10")
			legend(leg.loc, adj=-0.05, legend=c(otu.ids[which(rank(HCs)==d)], otu.ids[which(rank(HCs)==d-1)], otu.ids[which(rank(HCs)==d-2)],
			otu.ids[which(rank(HCs)==d-3)], otu.ids[which(rank(HCs)==d-4)], otu.ids[which(rank(HCs)==d-5)], 
			otu.ids[which(rank(HCs)==d-6)], otu.ids[which(rank(HCs)==d-7)], otu.ids[which(rank(HCs)==d-8)], otu.ids[which(rank(HCs)==d-9)]),
			col=c("gray1", "gray10", "gray20", "gray30", "gray40", "gray50", "gray60", "gray70", "gray80", "gray90"), title="10 most infl OTUs", pch=15, bty='n', cex=0.5)
	
			wHC.pvs <- round(c(MiHC.out$simes.pv, MiHC.out$ind.pvs[,2], MiHC.out$ada.pvs[2]),3)
			ind.0 <- which(wHC.pvs == 0)
			ind.1 <- which(wHC.pvs == 1)
			wHC.pvs[ind.0] <- "<.001"
			wHC.pvs[ind.1] <- ">.999"
			if (paste(rownames(MiHC.out$ind.pvs), collapse="") == c("13579")) {
				xlab.v <- paste("\n\nObserved quantiles\nwHC(1):", wHC.pvs[2], ", wHC(3):", wHC.pvs[3], ", wHC(5):", wHC.pvs[4], ", wHC(7):", wHC.pvs[5], ", wHC(9):", wHC.pvs[6], "\nSimes:", wHC.pvs[1], ", wHC(A):", wHC.pvs[7], sep="")
			} else {
				xlab.v <- paste("wHC(A):", wHC.pvs[length(wHC.pvs)], sep="")
			}
			plot(obs.wHC, exp.wHC, pch=20, panel.first=grid(lty=2, lwd=1), xlab=xlab.v, ylab="Expected quantiles", main="wHC", cex.lab=0.7, cex.axis=0.7, cex.main=0.9, col="blue2")

			abline(0, 1, col="red2")
			abline(v=obs.wHC[which(rank(wHCs)==d)], col="gray1")
			abline(v=obs.wHC[which(rank(wHCs)==d-1)], col="gray10")
			abline(v=obs.wHC[which(rank(wHCs)==d-2)], col="gray20")
			abline(v=obs.wHC[which(rank(wHCs)==d-3)], col="gray30")
			abline(v=obs.wHC[which(rank(wHCs)==d-4)], col="gray40")
			abline(v=obs.wHC[which(rank(wHCs)==d-5)], col="gray50")
			abline(v=obs.wHC[which(rank(wHCs)==d-6)], col="gray60")
			abline(v=obs.wHC[which(rank(wHCs)==d-7)], col="gray70")
			abline(v=obs.wHC[which(rank(wHCs)==d-8)], col="gray80")
			abline(v=obs.wHC[which(rank(wHCs)==d-9)], col="gray10")
			legend(leg.loc, adj=-0.05, legend=c(otu.ids[which(rank(wHCs)==d)], otu.ids[which(rank(wHCs)==d-1)], otu.ids[which(rank(wHCs)==d-2)],
			otu.ids[which(rank(wHCs)==d-3)], otu.ids[which(rank(wHCs)==d-4)], otu.ids[which(rank(wHCs)==d-5)], 
			otu.ids[which(rank(wHCs)==d-6)], otu.ids[which(rank(wHCs)==d-7)], otu.ids[which(rank(wHCs)==d-8)], otu.ids[which(rank(wHCs)==d-9)]),
			col=c("gray1", "gray10", "gray20", "gray30", "gray40", "gray50", "gray60", "gray70", "gray80", "gray90"), title="10 most infl OTUs", pch=15, bty='n', cex=0.5)
	
			MiHC.pv <- round(MiHC.out$ada.pvs[3],3)
			if (MiHC.pv == 0) MiHC.pv <- "<.001"
			if (MiHC.pv == 1) MiHC.pv <- ">.999"
			mtext(paste("*MiHC: ", MiHC.pv, sep=""), side=3, adj=0.91, cex=0.9, line=-1, font=2, outer=TRUE)
			graphics.off()
		} else {
			layout(matrix(c(1,2), nrow=2, byrow=TRUE))
			par(oma=c(0,0,1,0))

			uHC.pvs <- round(c(MiHC.out$simes.pv, MiHC.out$ind.pvs[,1], MiHC.out$ada.pvs[1]),3)
			ind.0 <- which(uHC.pvs == 0)
			ind.1 <- which(uHC.pvs == 1)
			uHC.pvs[ind.0] <- "<.001"
			uHC.pvs[ind.1] <- ">.999"
			if (paste(rownames(MiHC.out$ind.pvs), collapse="") == c("13579")) {
				xlab.v <- paste("\n\nObserved quantiles\nuHC(1):", uHC.pvs[2], ", uHC(3):", uHC.pvs[3], ", uHC(5):", uHC.pvs[4], ", uHC(7):", uHC.pvs[5], ", uHC(9):", uHC.pvs[6], "\nSimes:", uHC.pvs[1], ", uHC(A):", uHC.pvs[7], sep="")
			} else {
				xlab.v <- paste("uHC(A):", uHC.pvs[length(uHC.pvs)], sep="")
			}
			plot(obs.HC, exp.HC, pch=20, panel.first=grid(lty=2, lwd=1), xlab=xlab.v, ylab="Expected quantiles", main="uHC", cex.lab=0.7, cex.axis=0.7, cex.main=0.9, col="blue2")
			
			abline(0, 1, col="red2")
			abline(v=obs.HC[which(rank(HCs)==d)], col="gray1")
			abline(v=obs.HC[which(rank(HCs)==d-1)], col="gray10")
			abline(v=obs.HC[which(rank(HCs)==d-2)], col="gray20")
			abline(v=obs.HC[which(rank(HCs)==d-3)], col="gray30")
			abline(v=obs.HC[which(rank(HCs)==d-4)], col="gray40")
			abline(v=obs.HC[which(rank(HCs)==d-5)], col="gray50")
			abline(v=obs.HC[which(rank(HCs)==d-6)], col="gray60")
			abline(v=obs.HC[which(rank(HCs)==d-7)], col="gray70")
			abline(v=obs.HC[which(rank(HCs)==d-8)], col="gray80")
			abline(v=obs.HC[which(rank(HCs)==d-9)], col="gray10")
			legend(leg.loc, adj=-0.05, legend=c(otu.ids[which(rank(HCs)==d)], otu.ids[which(rank(HCs)==d-1)], otu.ids[which(rank(HCs)==d-2)],
			otu.ids[which(rank(HCs)==d-3)], otu.ids[which(rank(HCs)==d-4)], otu.ids[which(rank(HCs)==d-5)], 
			otu.ids[which(rank(HCs)==d-6)], otu.ids[which(rank(HCs)==d-7)], otu.ids[which(rank(HCs)==d-8)], otu.ids[which(rank(HCs)==d-9)]),
			col=c("gray1", "gray10", "gray20", "gray30", "gray40", "gray50", "gray60", "gray70", "gray80", "gray90"), title="10 most infl OTUs", pch=15, bty='n', cex=0.5)
	
			wHC.pvs <- round(c(MiHC.out$simes.pv, MiHC.out$ind.pvs[,2], MiHC.out$ada.pvs[2]),3)
			ind.0 <- which(wHC.pvs == 0)
			ind.1 <- which(wHC.pvs == 1)
			wHC.pvs[ind.0] <- "<.001"
			wHC.pvs[ind.1] <- ">.999"
			if (paste(rownames(MiHC.out$ind.pvs), collapse="") == c("13579")) {
				xlab.v <- paste("\n\nObserved quantiles\nwHC(1):", wHC.pvs[2], ", wHC(3):", wHC.pvs[3], ", wHC(5):", wHC.pvs[4], ", wHC(7):", wHC.pvs[5], ", wHC(9):", wHC.pvs[6], "\nSimes:", wHC.pvs[1], ", wHC(A):", wHC.pvs[7], sep="")
			} else {
				xlab.v <- paste("wHC(A):", wHC.pvs[length(wHC.pvs)], sep="")
			}
			plot(obs.wHC, exp.wHC, pch=20, panel.first=grid(lty=2, lwd=1), xlab=xlab.v, ylab="Expected quantiles", main="wHC", cex.lab=0.7, cex.axis=0.7, cex.main=0.9, col="blue2")

			abline(0, 1, col="red2")
			abline(v=obs.wHC[which(rank(wHCs)==d)], col="gray1")
			abline(v=obs.wHC[which(rank(wHCs)==d-1)], col="gray10")
			abline(v=obs.wHC[which(rank(wHCs)==d-2)], col="gray20")
			abline(v=obs.wHC[which(rank(wHCs)==d-3)], col="gray30")
			abline(v=obs.wHC[which(rank(wHCs)==d-4)], col="gray40")
			abline(v=obs.wHC[which(rank(wHCs)==d-5)], col="gray50")
			abline(v=obs.wHC[which(rank(wHCs)==d-6)], col="gray60")
			abline(v=obs.wHC[which(rank(wHCs)==d-7)], col="gray70")
			abline(v=obs.wHC[which(rank(wHCs)==d-8)], col="gray80")
			abline(v=obs.wHC[which(rank(wHCs)==d-9)], col="gray10")
			legend(leg.loc, adj=-0.05, legend=c(otu.ids[which(rank(wHCs)==d)], otu.ids[which(rank(wHCs)==d-1)], otu.ids[which(rank(wHCs)==d-2)],
			otu.ids[which(rank(wHCs)==d-3)], otu.ids[which(rank(wHCs)==d-4)], otu.ids[which(rank(wHCs)==d-5)], 
			otu.ids[which(rank(wHCs)==d-6)], otu.ids[which(rank(wHCs)==d-7)], otu.ids[which(rank(wHCs)==d-8)], otu.ids[which(rank(wHCs)==d-9)]),
			col=c("gray1", "gray10", "gray20", "gray30", "gray40", "gray50", "gray60", "gray70", "gray80", "gray90"), title="10 most infl OTUs", pch=15, bty='n', cex=0.5)
	
			MiHC.pv <- round(MiHC.out$ada.pvs[3],3)
			if (MiHC.pv == 0) MiHC.pv <- "<.001"
			if (MiHC.pv == 1) MiHC.pv <- ">.999"
			mtext(paste("*MiHC: ", MiHC.pv, sep=""), side=3, adj=0.91, cex=0.9, line=-1, font=2, outer=TRUE)
		}
	}
	if (length(MiHC.out$graph.els) == 3) {
		otu.ids <- MiHC.out$graph.els$otu.ids
		exp.HC <- MiHC.out$graph.els$exp.HC
		obs.HC <- MiHC.out$graph.els$obs.HC	
		HCs <- abs(exp.HC-obs.HC)
		d <- length(otu.ids)
		if (!is.null(pdf.filename)) {
			pdf(file=pdf.filename, width=4.5, height=5.5)
			layout(matrix(c(1), nrow=1, byrow=TRUE))
			par(oma=c(0,0,1,0))
			uHC.pvs <- round(c(MiHC.out$simes.pv, MiHC.out$ind.pvs[,1], MiHC.out$ada.pvs[1]),3)
			ind.0 <- which(uHC.pvs == 0)
			ind.1 <- which(uHC.pvs == 1)
			uHC.pvs[ind.0] <- "<.001"
			uHC.pvs[ind.1] <- ">.999"
			if (paste(rownames(MiHC.out$ind.pvs), collapse="") == c("13579")) {
				xlab.v <- paste("\n\nObserved quantiles\nuHC(1):", uHC.pvs[2], ", uHC(3):", uHC.pvs[3], ", uHC(5):", uHC.pvs[4], ", uHC(7):", uHC.pvs[5], ", uHC(9):", uHC.pvs[6], "\nSimes:", uHC.pvs[1], ", uHC(A):", uHC.pvs[7], sep="")
			} else {
				xlab.v <- paste("uHC(A):", uHC.pvs[length(uHC.pvs)], sep="")
			}
			plot(obs.HC, exp.HC, pch=20, panel.first=grid(lty=2, lwd=1), xlab=xlab.v, ylab="Expected quantiles", main="uHC", cex.lab=0.7, cex.axis=0.7, cex.main=0.9, col="blue2")

			abline(0, 1, col="red2")
			abline(v=obs.HC[which(rank(HCs)==d)], col="gray1")
			abline(v=obs.HC[which(rank(HCs)==d-1)], col="gray10")
			abline(v=obs.HC[which(rank(HCs)==d-2)], col="gray20")
			abline(v=obs.HC[which(rank(HCs)==d-3)], col="gray30")
			abline(v=obs.HC[which(rank(HCs)==d-4)], col="gray40")
			abline(v=obs.HC[which(rank(HCs)==d-5)], col="gray50")
			abline(v=obs.HC[which(rank(HCs)==d-6)], col="gray60")
			abline(v=obs.HC[which(rank(HCs)==d-7)], col="gray70")
			abline(v=obs.HC[which(rank(HCs)==d-8)], col="gray80")
			abline(v=obs.HC[which(rank(HCs)==d-9)], col="gray10")
			legend(leg.loc, adj=-0.05, legend=c(otu.ids[which(rank(HCs)==d)], otu.ids[which(rank(HCs)==d-1)], otu.ids[which(rank(HCs)==d-2)],
			otu.ids[which(rank(HCs)==d-3)], otu.ids[which(rank(HCs)==d-4)], otu.ids[which(rank(HCs)==d-5)], 
			otu.ids[which(rank(HCs)==d-6)], otu.ids[which(rank(HCs)==d-7)], otu.ids[which(rank(HCs)==d-8)], otu.ids[which(rank(HCs)==d-9)]),
			col=c("gray1", "gray10", "gray20", "gray30", "gray40", "gray50", "gray60", "gray70", "gray80", "gray90"), title="10 most infl OTUs", pch=15, bty='n', cex=0.5)
			
			MiHC.pv <- round(MiHC.out$ada.pvs[1],3)
			if (MiHC.pv == 0) MiHC.pv <- "<.001"
			if (MiHC.pv == 1) MiHC.pv <- ">.999"
			mtext(paste("*uHC(A): ", MiHC.pv, sep=""), side=3, adj=0.91, cex=0.9, line=-1, font=2, outer=TRUE)
			graphics.off()
		} else {
			layout(matrix(c(1), nrow=1, byrow=TRUE))
			par(oma=c(0,0,1,0))
			uHC.pvs <- round(c(MiHC.out$simes.pv, MiHC.out$ind.pvs[,1], MiHC.out$ada.pvs[1]),3)
			ind.0 <- which(uHC.pvs == 0)
			ind.1 <- which(uHC.pvs == 1)
			uHC.pvs[ind.0] <- "<.001"
			uHC.pvs[ind.1] <- ">.999"
			if (paste(rownames(MiHC.out$ind.pvs), collapse="") == c("13579")) {
				xlab.v <- paste("\n\nObserved quantiles\nuHC(1):", uHC.pvs[2], ", uHC(3):", uHC.pvs[3], ", uHC(5):", uHC.pvs[4], ", uHC(7):", uHC.pvs[5], ", uHC(9):", uHC.pvs[6], "\nSimes:", uHC.pvs[1], ", uHC(A):", uHC.pvs[7], sep="")
			} else {
				xlab.v <- paste("uHC(A):", uHC.pvs[length(uHC.pvs)], sep="")
			}
			plot(obs.HC, exp.HC, pch=20, panel.first=grid(lty=2, lwd=1), xlab=xlab.v, ylab="Expected quantiles", main="uHC", cex.lab=0.7, cex.axis=0.7, cex.main=0.9, col="blue2")
			
			abline(0, 1, col="red2")
			abline(v=obs.HC[which(rank(HCs)==d)], col="gray1")
			abline(v=obs.HC[which(rank(HCs)==d-1)], col="gray10")
			abline(v=obs.HC[which(rank(HCs)==d-2)], col="gray20")
			abline(v=obs.HC[which(rank(HCs)==d-3)], col="gray30")
			abline(v=obs.HC[which(rank(HCs)==d-4)], col="gray40")
			abline(v=obs.HC[which(rank(HCs)==d-5)], col="gray50")
			abline(v=obs.HC[which(rank(HCs)==d-6)], col="gray60")
			abline(v=obs.HC[which(rank(HCs)==d-7)], col="gray70")
			abline(v=obs.HC[which(rank(HCs)==d-8)], col="gray80")
			abline(v=obs.HC[which(rank(HCs)==d-9)], col="gray10")
			legend(leg.loc, adj=-0.05, legend=c(otu.ids[which(rank(HCs)==d)], otu.ids[which(rank(HCs)==d-1)], otu.ids[which(rank(HCs)==d-2)],
			otu.ids[which(rank(HCs)==d-3)], otu.ids[which(rank(HCs)==d-4)], otu.ids[which(rank(HCs)==d-5)], 
			otu.ids[which(rank(HCs)==d-6)], otu.ids[which(rank(HCs)==d-7)], otu.ids[which(rank(HCs)==d-8)], otu.ids[which(rank(HCs)==d-9)]),
			col=c("gray1", "gray10", "gray20", "gray30", "gray40", "gray50", "gray60", "gray70", "gray80", "gray90"), title="10 most infl OTUs", pch=15, bty='n', cex=0.5)
			
			MiHC.pv <- round(MiHC.out$ada.pvs[1],3)
			if (MiHC.pv == 0) MiHC.pv <- "<.001"
			if (MiHC.pv == 1) MiHC.pv <- ">.999"
			mtext(paste("*uHC(A): ", MiHC.pv, sep=""), side=3, adj=0.91, cex=0.9, line=-1, font=2, outer=TRUE)
		}
	}
}
