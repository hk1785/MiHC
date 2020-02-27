MiHC <-
function(y, covs=NULL, otu.tab, tree, model, hs=c(1,3,5,7,9), W=TRUE, comp=FALSE, CLR=FALSE, opt.ncl=30, n.perm=5000) {

	if (!comp) {
		otu.tab <- t(apply(otu.tab,1,function(x)x/sum(x)))
	}
	com.otu.tab <- otu.tab
	if (CLR) {
		otu.tab <- t(apply(com.otu.tab, 1, clr))
	}
	n <- length(y)
	p <- ncol(otu.tab)
	
    if (is.null(covs)) {
        fit <- glm(y ~ 1, family = model)
		f.y <- fitted.values(fit)
        r <- y - f.y
    } else {
		fit <- glm(y ~ ., family = model, as.data.frame(covs))
		f.y <- fitted.values(fit)
        r <- y - f.y
    }

	r.s <- list()
	for (j in 1:n.perm) {
		r.s[[j]] <- r[shuffle(n)]
	}
	
	n.Zs <- as.numeric(r %*% otu.tab)
	n.Z0s <- lapply(r.s, function(x) as.numeric(x %*% otu.tab))
	U0s <- lapply(apply(sapply(r.s, function(x) return(x %*% otu.tab)),1,list),unlist)
	se <- sapply(U0s, sd)
	Zs <- n.Zs/se
	Z0s <- lapply(n.Z0s, function(x) x/se)
	
	if (W) {
		CDis <- cophenetic(tree)	
		for (j in 1:nrow(CDis)) {
			ind.t <- which(CDis[j,] == 0)
			ind.d <- which(ind.t == j)
			ind.r <- ind.t[-ind.d]
			CDis[j,ind.r] <- min(CDis[j,-ind.t])/2
		}	
		
		if (opt.ncl == "gmx") {
			asw <- numeric(p-1)
			for (j in 2:(p-1)) {
				asw[j] <- pam(CDis, j)$silinfo$avg.width
			}
			k.best <- which.max(asw)
			clust.pam <- pam(CDis, k.best, cluster.only=TRUE)
		}
		if (opt.ncl == "fmx") {
			asw <- numeric(p-1)
			j <- 2
			asw[j] <- pam(CDis, j)$silinfo$avg.width
			while (asw[j-1] < asw[j] & j <= p-1) {
				j <- j + 1
				asw[j] <- pam(CDis, j)$silinfo$avg.width
			}		
			k.best <- j-1
			clust.pam <- pam(CDis, k.best, cluster.only=TRUE)
		}
		if (opt.ncl != "gmx" & opt.ncl != "fmx") {
			asw <- numeric(p-1)
			for (j in 2:opt.ncl) {
				asw[j] <- pam(CDis, j)$silinfo$avg.width
			}
			k.best <- which.max(asw)
			clust.pam <- pam(CDis, k.best, cluster.only=TRUE)
		}
		Ws <- rep(NA, p)
		for (j in 1:p) {
			ind.1 <- match(colnames(CDis)[j], names(clust.pam))
			ind.2 <- which(clust.pam == as.numeric(clust.pam[ind.1]))
			ind.3 <- which(ind.2 == ind.1)
			ind <- ind.2[-ind.3]
			if (length(ind) == 0) {
				Ws[j] <- 1
			} else {
				inv.Ds <- 1/(CDis[j,ind])
				abs.Zs <- abs(Zs[ind])
				Ws[j] <- sum(inv.Ds*abs.Zs)/sum(inv.Ds) + 1
			}
		}
		Ws <- Ws/sum(Ws)

		ahc.o <- as.list(MiHC.stat(Zs=Zs, hs=hs, Ws=Ws))
		ahc.p <- lapply(apply(sapply(Z0s, function(x) MiHC.stat(Zs=x, hs=hs, Ws=Ws)), 1, list), unlist)
		pvs <- mapply(function(x, y) (length(which(x < y))+0.01)/(n.perm+0.01), ahc.o, ahc.p)
	
		ind.pvs <- 1-pchisq(Zs^2, df=1)
		simes.pv <- round(min(length(ind.pvs)*ind.pvs/rank(ind.pvs)),4)
		simes.pv0s.ori <- unlist(lapply(Z0s, function(x) {xx <- 1-pchisq(x^2, df=1); min(length(xx)*xx/rank(xx))}))
		simes.pv0s <- rep(NA, n.perm)
		for (j in 1:n.perm) {	
			simes.pv0s[j] <- (length(which(simes.pv0s.ori[-j] < simes.pv0s.ori[j]))+0.01)/(n.perm+0.01)
		}
		simes.pv <- (length(which(simes.pv0s < simes.pv))+0.01)/(n.perm+0.01)
		
		l.hs <- length(hs)
		Tu <- min(pvs[1:l.hs], simes.pv)
		T0u <- rep(NA, n.perm)
		Tw <- min(pvs[(l.hs+1):(l.hs*2)], simes.pv)
		T0w <- rep(NA, n.perm)
		for (j in 1:n.perm) {
			T0u.o <- lapply(ahc.p[1:l.hs], function(x) x[j])
			T0u.p <- lapply(ahc.p[1:l.hs], function(x) x[-j])
			T0u[j] <- min(mapply(function(x, y) (length(which(x < y))+0.01)/(n.perm+0.01), T0u.o, T0u.p))
		}
		for (j in 1:n.perm) {
			T0w.o <- lapply(ahc.p[(l.hs+1):(l.hs*2)], function(x) x[j])
			T0w.p <- lapply(ahc.p[(l.hs+1):(l.hs*2)], function(x) x[-j])
			T0w[j] <- min(mapply(function(x, y) (length(which(x < y))+0.01)/(n.perm+0.01), T0w.o, T0w.p))
		}
		T0u.minp <- apply(cbind(T0u, simes.pv0s),1,min)
		T0w.minp <- apply(cbind(T0w, simes.pv0s),1,min)
		pv.opt.u <- (length(which(T0u.minp < Tu))+0.01)/(n.perm+0.01)
		pv.opt.w <- (length(which(T0w.minp < Tw))+0.01)/(n.perm+0.01)
	
		T <- min(pvs, simes.pv)
		T0 <- apply(cbind(T0u, T0w, simes.pv0s),1,min)
		pv.opt <- (length(which(T0 < T))+0.01)/(n.perm+0.01)
		
		names(simes.pv) <- "Simes"
		item.pvs <- cbind(pvs[1:l.hs], pvs[(l.hs+1):(l.hs*2)])
		colnames(item.pvs) <- c("uHC(h)", "wHC(h)")
		rownames(item.pvs) <- as.character(hs)
		
		ada.pvs <- c(pv.opt.u, pv.opt.w, pv.opt)
		names(ada.pvs) <- c("uHC(A)", "wHC(A)", "MiHC")
		
		otu.ids <- colnames(otu.tab)
		otu.ids <- sub("New.ReferenceOTU", "N", otu.ids)
		otu.ids <- sub("New.CleanUp.ReferenceOTU", "NC", otu.ids)
		pvs <- ind.pvs
		i0 <- which(pvs < 0.00000001)
		i1 <- which(pvs > 0.99999999)
		pvs[i0] <- 0.00000001
		pvs[i1] <- 0.99999999
		Is <- order(order(pvs))
		exp.HC <- (Is/p)/sqrt(pvs*(1-pvs)/p)
		obs.HC <- pvs/sqrt(pvs*(1-pvs)/p)
		exp.wHC <- Ws*(Is/p)/sqrt(pvs*(1-pvs)/p)
		obs.wHC <- Ws*pvs/sqrt(pvs*(1-pvs)/p)
		return(list(graph.els = list(otu.ids = otu.ids, exp.HC = exp.HC, obs.HC = obs.HC, exp.wHC = exp.wHC, obs.wHC = obs.wHC), opt.nclust = k.best, simes.pv = simes.pv, ind.pvs = item.pvs, ada.pvs = ada.pvs))
	}
	
	if (!W) {
		ahc.o <- as.list(MiHC.stat(Zs=Zs, hs=hs, Ws=NULL))
		ahc.p <- lapply(apply(sapply(Z0s, function(x) MiHC.stat(Zs=x, hs=hs, Ws=NULL)), 1, list), unlist)
		pvs <- mapply(function(x, y) (length(which(x < y))+0.01)/(n.perm+0.01), ahc.o, ahc.p)
	
		ind.pvs <- 1-pchisq(Zs^2, df=1)
		simes.pv <- min(length(ind.pvs)*ind.pvs/rank(ind.pvs))
		simes.pv0s.ori <- unlist(lapply(Z0s, function(x) {xx <- 1-pchisq(x^2, df=1); min(length(xx)*xx/rank(xx))}))
		simes.pv0s <- rep(NA, n.perm)
		for (j in 1:n.perm) {	
			simes.pv0s[j] <- (length(which(simes.pv0s.ori[-j] < simes.pv0s.ori[j]))+0.01)/(n.perm+0.01)
		}
		simes.pv <- (length(which(simes.pv0s < simes.pv))+0.01)/(n.perm+0.01)
		
		l.hs <- length(hs)
		Tu <- min(pvs, simes.pv)
		T0u <- rep(NA, n.perm)
		for (j in 1:n.perm) {
			T0u.o <- lapply(ahc.p, function(x) x[j])
			T0u.p <- lapply(ahc.p, function(x) x[-j])
			T0u[j] <- min(mapply(function(x, y) (length(which(x < y))+0.01)/(n.perm+0.01), T0u.o, T0u.p))
		}
		T0u.minp <- apply(cbind(T0u, simes.pv0s),1,min)
		pv.opt.u <- (length(which(T0u.minp < Tu))+0.01)/(n.perm+0.01)

		names(simes.pv) <- "Simes"
		item.pvs <- as.matrix(pvs)
		colnames(item.pvs) <- "uHC(h)"
		rownames(item.pvs) <- as.character(hs)
		
		ada.pvs <- pv.opt.u
		names(ada.pvs) <- "uHC(A)"
		
		otu.ids <- colnames(otu.tab)
		otu.ids <- sub("New.ReferenceOTU", "N", otu.ids)
		otu.ids <- sub("New.CleanUp.ReferenceOTU", "NC", otu.ids)
		pvs <- ind.pvs
		i0 <- which(pvs < 0.00000001)
		i1 <- which(pvs > 0.99999999)
		pvs[i0] <- 0.00000001
		pvs[i1] <- 0.99999999
		Is <- order(order(pvs))
		exp.HC <- (Is/p)/sqrt(pvs*(1-pvs)/p)
		obs.HC <- pvs/sqrt(pvs*(1-pvs)/p)
		return(list(graph.els = list(otu.ids = otu.ids, exp.HC = exp.HC, obs.HC = obs.HC), simes.pv = simes.pv, ind.pvs = item.pvs, ada.pvs = ada.pvs))
	}
}
