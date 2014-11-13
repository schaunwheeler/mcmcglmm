# File-Name:       mcmcglmm.R           
# Date:            2012-09-21                     
# Author:          Schaun Wheeler
# Email:           schaun.wheeler@gmail.com                                      
# Purpose:         To expand the capabilities of the MCMCglmm package,
#                  particularly as it relates to issues of cross-validation.
# Data Used:       none
# Packages Used:   MCMCglmm
# Output File:     
# Data Output:     
# Machine:         Schaun Wheeler's Dell Precision T7500 and MacBook Pro

# Copyright (c) 2012, under the Simplified BSD License.  
# For more information on FreeBSD see: 
# http://www.opensource.org/licenses/bsd-license.php
# All rights reserved.

# THIS PACKAGE IS NO LONGER ACTIVELY MAINTAINED. IT ONLY EVER WORKED FOR
# THE "GAUSSIAN" FAMILY, AND EVEN THAT WAS NOT TESTED MUCH.

# Load libraries
library(MCMCglmm)

# Load functions
SplitData <- function(data, percent = .8, ignore = NULL){

	facs <- sapply(data,is.factor)
	data[,facs] <- lapply(data[,facs],as.character)
	
	chars <- sapply(data,is.character)
	ignore <- colnames(data) %in% ignore
	
	look <- chars & !ignore

	num <- round(nrow(data) * percent, 0)
	
	rows <- sample(1:nrow(data), num)
	rowind <- 1:nrow(data) %in% rows

	big <- data[rowind,]
	small <- data[!rowind,]
	
	bigind <- 1:nrow(big)
	smallind <- 1:nrow(small)
	
	bigvals <- rep(NA,ncol(big))
	smallvals <- rep(NA,ncol(small))
	
	bigvals[look] <- lapply(big[,look],function(x)sort(unique(x)))
	smallvals[look] <- lapply(small[,look],function(x)sort(unique(x)))
	
	matches <- lapply(1:length(smallvals), 
										function(x)smallvals[[x]] %in% bigvals[[x]])
	
	misses <- which(sapply(matches, function(x)1-mean(x)) > 0)
	
	missvals <- lapply(1:length(smallvals), 
										function(x)smallvals[[x]][!(smallvals[[x]] %in% bigvals[[x]])])

	if(length(misses) > 0){
		for(i in misses){
			for(j in 1:length(missvals[[i]])){
			pulls <- smallind[small[,i] == missvals[[i]][j]]
			take <- ifelse(length(pulls) == 1, pulls, 
										 try(sample(pulls, 1), silent = T))
			if(is.numeric(take)){
			big <- rbind(big,small[take,])
			small <- small[-take,]
			}
			}
		}
		}
	list("large" = big,
			 "small" = small)	
}

mcmcglmm <- function(fixed, random=NULL, rcov=~units, family="gaussian", 
										 mev=NULL, data, start=NULL, prior="InvW-pe", tune=NULL, 
										 pedigree=NULL, nodes="ALL", scale=TRUE, nitt=13000, 
										 thin=10, burnin=3000, pr=TRUE, pl=FALSE, verbose=TRUE, 
										 DIC=TRUE, singular.ok=FALSE, saveX=TRUE, saveZ=TRUE, 
										 saveXL=TRUE, slice=FALSE, ginverse=NULL,...){
	
	chars <- sapply(data,is.character)
	data[,chars] <- lapply(data[,chars],as.factor)
	datalevels <- sapply(data,levels)

	if(is.null(random)){
		pr = FALSE
	}
	
	split.direct.sum <- function (x){
		if (is.na(x)) {
			return(NULL)
		}	else {
			openB <- gregexpr("\\(", x)[[1]]
			closeB <- gregexpr("\\)", x)[[1]]
			true_openB <- openB
			true_closeB <- closeB
			for (i in 1:length(openB)) {
				dist <- outer(openB, closeB, function(x, y) {
					y - x
				})
				dist[which(dist < 0)] <- Inf
				dist <- which(dist == min(dist), arr.ind = T)[1, 
																											]
				true_openB[i] <- openB[dist[1]]
				true_closeB[i] <- closeB[dist[2]]
				openB <- openB[-dist[1]]
				closeB <- closeB[-dist[2]]
			}
			plus <- gregexpr("\\+", x)[[1]]
			internals <- matrix(mapply(function(x, y) {
				(plus > x & plus < y)
			}, x = true_openB, y = true_closeB), length(plus), length(true_openB))
			rterms <- strsplit(x, "")[[1]]
			rterms[plus[which(rowSums(internals) != 0)]] <- "leaveMCMCleave"
			rterms <- paste(rterms, collapse = "")
			rterms <- strsplit(rterms, " *\\+ *")[[1]]
			rterms <- gsub("leaveMCMCleave", "+", rterms)
			return(rterms)
		}
	}
	
  if(is.character(prior)){
	
	fixvar <- split.direct.sum(as.character(fixed)[3])
	ranvar <- split.direct.sum(as.character(random)[2])

	f <- sapply(fixvar,function(x){
		if(grepl(":",x)){
			frt <- gsub("(\\w+):(\\w+)", "\\1", x)
			bck <- gsub("(\\w+):(\\w+)", "\\2", x)
			levfrt <- length(datalevels[[frt]])
			levbck <- length(datalevels[[bck]])
			y <- levfrt * levbck
		}else{
			y <- length(datalevels[[x]])
		}
		if(y == 0){
			y <- 1
		}
		y
	})
    
    if(is.null(fixvar)){
      f <- 0
    }
	
  r <- sapply(ranvar,function(x){
    if(grepl("(us|idh|cor)\\(1 [+] ", x)){
      y <- 2
    }
    if(grepl("(us|idh|cor)\\((?!1 [+] )", x, perl = T)){
      y <- length(unlist(datalevels[gsub("(us|idh|cor)\\((\\w+)\\):\\w+", "\\2", x)]))
      if(length(y) == 0){
        y <- 1
        }
      }
    if(!grepl("(us|idh|cor)", x)){
      y <- 1
      }
    y
	})
	
	names(r)<-NULL
  
	B = list(mu = rep(0,sum(f)+1), V = diag(sum(f)+1)*(1e10))
	
	if(grepl("\\bInvW", prior)){
	  R = list(V = 1, nu = 1)
	  G = lapply(r, function(x){l1 <- list(V = diag(x), nu = x)
	  													if(grepl("-pe\\b", prior)){
	  														l1 <- c(l1,
	  																		list(alpha.mu=rep(0,x), 
	  																		alpha.V=diag(x)*1000))
	  														l1}})
	}
	if(grepl("\\bInvG", prior)){
	  R = list(V = 1, nu = 0.002)  
	  G = lapply(r, function(x){l1 <- list(V = diag(x), nu = x-1+0.002)
	  													if(grepl("-pe\\b", prior)){
	  														l1 <- c(l1,
	  																		list(alpha.mu=rep(0,x), 
	  																				 alpha.V=diag(x)*1000))
	  														l1}})
	}
	
  if(grepl("ordinal|categorical|multinomial", family)){
     R = list(fix = 1, V = 1, nu = 1) 
  }
  
	if(!is.null(random)){
	names(G) <- paste("G",1:length(ranvar),sep="")
	}
  prior = list("B" = B, "R" = R, "G" = G)
  
}

	m1 <- MCMCglmm(fixed=fixed, random=random, rcov=rcov, family=family, 
	               mev=mev, data=data, start=start, prior=prior, tune=tune, 
	               pedigree=pedigree, nodes=nodes, scale=scale, nitt=nitt, 
	               thin=thin, burnin=burnin, pr=pr, pl=pl, verbose=verbose, 
	               DIC=DIC, singular.ok=singular.ok, saveX=saveX, saveZ=saveZ, 
	               saveXL=saveXL, slice=slice, ginverse=ginverse)
	
	m1$datalevels <- datalevels
	
	m1
}
										 
QuickSummary <- function(out, prob = .95, m = 2, rnd = 4){
	out.mean <- colMeans(out$Sol)
	out.hpd <- HPDinterval(out$Sol, prob = prob)
	out.types <- ifelse(colMeans(out$Sol) > 0,
											colMeans(out$Sol < 0),
											colMeans(out$Sol > 0))
	raw.typem <- ifelse(colMeans(out$Sol) > 0,
											colMeans(t(apply(out$Sol, MARGIN = 1, 
																			 function(x){
																			 	x[x < 0] <- 0
																			 	out.mean/m - x})) > 0),
											colMeans(t(apply(out$Sol, MARGIN = 1, 
																			 function(x){
																			 	x[x > 0] <- 0
																			 	out.mean/m - x})) < 0))
	out.typem <- raw.typem - out.types

	round(cbind("mean" = out.mean,
							out.hpd,
							"type_S" = out.types,
							"type_M" = out.typem), digits = rnd)
}

PredictNew <- function (object, newdata = NULL, marginal = NULL, type = "terms", 
												interval = "none", level = 0.95, index = NULL, ...) {
	
	# Do checks to make sure all the inputs are in order.
	
	if (type %in% c("response", "terms") == FALSE) {
		stop("type must be response or terms")
	}
	if (interval %in% c("none", "confidence", "prediction", "all") == 
		FALSE) {
		stop("interval must be none, confidence or prediction")
	}
	
	# Load function to identify model components
	
	split.direct.sum <- function (x){
		if (is.na(x)) {
			return(NULL)
		}	else {
			openB <- gregexpr("\\(", x)[[1]]
			closeB <- gregexpr("\\)", x)[[1]]
			true_openB <- openB
			true_closeB <- closeB
			for (i in 1:length(openB)) {
				dist <- outer(openB, closeB, function(x, y) {
					y - x
				})
				dist[which(dist < 0)] <- Inf
				dist <- which(dist == min(dist), arr.ind = T)[1, 
																											]
				true_openB[i] <- openB[dist[1]]
				true_closeB[i] <- closeB[dist[2]]
				openB <- openB[-dist[1]]
				closeB <- closeB[-dist[2]]
			}
			plus <- gregexpr("\\+", x)[[1]]
			internals <- matrix(mapply(function(x, y) {
				(plus > x & plus < y)
			}, x = true_openB, y = true_closeB), length(plus), length(true_openB))
			rterms <- strsplit(x, "")[[1]]
			rterms[plus[which(rowSums(internals) != 0)]] <- "leaveMCMCleave"
			rterms <- paste(rterms, collapse = "")
			rterms <- strsplit(rterms, " *\\+ *")[[1]]
			rterms <- gsub("leaveMCMCleave", "+", rterms)
			return(rterms)
		}
	}
	
	# Identify which random variables should be marginalized
	
	rcomponents <- split.direct.sum(as.character(object$Random$formula)[2])
	mcomponents <- split.direct.sum(as.character(marginal)[2])
	if (any(mcomponents %in% rcomponents == FALSE)) {
		stop("marginal formula does not correspond to model formula")
	}
	marginal <- rep(rep(as.numeric(rcomponents %in% mcomponents), 
											object$Random$nrt), object$Random$nfl)
	
	if (any(marginal == 0) & dim(object$Sol)[2] == dim(object$X)[2]) {
		stop("posterior distribution of random effects not saved: pass pr=TRUE to MCMCglmm")
	}
	if (any(marginal == 0) & is.null(object$Z)) {
		stop("random effect design matrix not saved: pass saveZ=TRUE to MCMCglmm")
	}
	if (is.null(object$X) & is.null(newdata)) {
		stop("either specify saveX=TRUE when model fitting or pass a data frame")
	}
	if (is.null(object$Random$nfl) == FALSE) {
		st <- c(1, cumsum(rep(object$Random$nrl, object$Random$nfl)) + 
			1)
		st <- st[-length(st)]
		end <- cumsum(rep(object$Random$nrl, object$Random$nfl))
		comp <- rep(1:length(object$Random$nfl), object$Random$nfl)
		keep <- unlist(mapply(st[which(marginal == 0)], end[which(marginal == 
			0)], FUN = ":"))
	} else {
		keep <- NULL
	}
	
	# Keep solutions only for non-marginalized variables
	
	object$Sol <- object$Sol[, c(1:object$Fixed$nfl, object$Fixed$nfl + 
		keep), drop = FALSE]
	W <- cBind(object$X, object$Z)
	W <- W[, c(1:object$Fixed$nfl, object$Fixed$nfl + keep), 
				 drop = FALSE]
	
	# Insert new data if new data exists
	
	if(!is.null(newdata)){
		chars <- sapply(newdata,is.character)
		newdata[,chars] <- lapply(newdata[,chars, drop = FALSE],as.factor)
		
		vars.o <- paste(as.character(object$Fixed[[1]])[-c(1:2)],
										as.character(object$Random[[1]]), collapse = " ")
		vars.o <- gsub("~|(us|idh|cor)\\(|[+]|\\):|\\b1\\b"," ", vars.o)
		vars.o <- unlist(strsplit(vars.o, split = "\\s+"))
		vars.o <- unique(vars.o[vars.o != ""])
		
		if(any(vars.o %in% colnames(newdata)) == F){
			stop("'newdata' is missing variables needed for the model")
		}
		
		facs <- sapply(newdata,is.factor)
		facs.o <- vars.o[vars.o %in% names(facs)[facs]]
		
		for(i in 1:length(facs.o)){
			newdata[,facs.o[i]] <- factor(newdata[,facs.o[i]], 
																		levels = sort(unique(c(levels(newdata[,facs.o[i]]), 
																							 object$datalevels[[facs.o[i]]]))),
																		labels = object$datalevels[[facs.o[i]]])
		}
		
		fixef <- sparse.model.matrix(
      as.formula(paste("~",as.character(object$Fixed[[1]])[-c(1:2)])), 
      newdata)
		
		rterms <- split.direct.sum(as.character(object$Random[[1]])[2])
		
		ranef <- lapply(rterms,function(x, df = newdata){
			covms <- grepl("\\w{2,3}\\([[:print:]]+\\):",x)
			ints <- grepl("\\w{2,3}\\((1 [+] )?([[:print:]]+)\\):([[:print:]]+)", x)
			if(covms == T & ints == T){
				full <- sparse.model.matrix(as.formula(
					gsub("\\w{2,3}\\(1 [+] ([[:print:]]+)\\):([[:print:]]+)", 
							 "~ 0 + \\1 : \\2", x)),df)
				binary <- full!=0
				matching <- vector("logical",length(colnames(df)))
				for(j in 1:length(colnames(df))){
					matching[j] <- grepl(paste(":",colnames(df)[j],sep=""), x)
				}
				matchvar <- colnames(df)[matching]
				firstvar <- gsub("\\w{2,3}\\(1 [+] ([[:print:]]+)\\):([[:print:]]+)", 
												 "\\1", x)
				colnames(binary) <- paste(matchvar, "(Intercept)", matchvar, 
																	sort(object$datalevels[[matchvar]]),sep=".")
				colnames(full) <- paste(matchvar, firstvar, matchvar, 
																sort(object$datalevels[[matchvar]]),sep=".")
				
				out <- cBind(binary,full)
			}
			if(covms == T & ints == F){
				full <- sparse.model.matrix(as.formula(
					gsub("\\w{2,3}\\(([[:print:]]+)\\):([[:print:]]+)", 
							 "~ 0 + \\1 : \\2", x)),df)
				matching <- vector("logical",length(colnames(df)))
				for(j in 1:length(colnames(df))){
					matching[j] <- grepl(paste(":",colnames(df)[j],sep=""), x)
				}
				matchvar <- colnames(df)[matching]
				firstvar <- gsub("\\w{2,3}\\(1 [+] ([[:print:]]+)\\):([[:print:]]+)", 
												 "\\1", rterms[i])
				colnames(full) <- paste(matchvar, firstvar, matchvar, 
																sort(unique(as.character(object$datalevels[[matchvar]]))),sep=".")
				
				out <- full
			}
			if(covms == F & ints == F){
				matchvar <- colnames(df)[colnames(df) %in% x]
				full <- sparse.model.matrix(as.formula(paste("~ 0 +", x, sep= "")),df)
				colnames(full) <- paste(x,sort(unique(as.character(object$datalevels[[matchvar]]))),sep=".")
				
				out <- full
			}
			out
		})
		
		ranef <- do.call("cBind",ranef)
		
		Wn <- cBind(fixef,ranef)
		
		object$X <- fixef[,match(colnames(object$X),colnames(fixef))]
		
		object$Z <- ranef[,match(colnames(object$Z),colnames(ranef))]
		
		object$error.term <- rep(object$error.term[1],nrow(Wn))
		
		W <- Wn[,match(colnames(W),colnames(Wn))]
		
		W <- W[, c(1:object$Fixed$nfl, object$Fixed$nfl + keep), 
					 drop = FALSE]
	}
	
	# Determine variance for posterior predictive simulation
	
	if ((type == "response" & any(object$family != "gaussian" & 
		object$family != "cengaussian")) | interval == "prediction") {
		vpred <- matrix(0, dim(object$X)[1], sum(object$Random$nfl[which(marginal == 
			1)]^2, na.rm = T) + sum(object$Residual$nfl^2))
		cnt <- 0
		if (any(marginal == 1)) {
			st <- st[which(marginal == 1)]
			end <- end[which(marginal == 1)]
			comp <- comp[which(marginal == 1)]
			for (i in 1:length(st)) {
				for (j in 1:length(st)) {
					if (comp[i] == comp[j]) {
						cnt <- cnt + 1
						vpred[, cnt] <- diag(object$Z[, st[i]:end[i]] %*% 
							t(object$Z[, st[j]:end[j]]))
					}
				}
			}
		}
		comp <- rep(1:length(object$Residual$nfl), object$Residual$nfl)
		for (i in 1:length(comp)) {
			for (j in 1:length(comp)) {
				if (comp[i] == comp[j]) {
#					cnt <- cnt + 1             # old code doesn't identify residual column
					cnt <- sum(object$Random$nfl[which(marginal == 1)]^2, na.rm = T) + 1
					vpred[, cnt][which(object$error.term == i & 
						object$error.term == j)] <- 1
				}
			}
		}

# Original code: fails whenever marginal = TRUE.
# 		keep <- which(rep(rep(as.numeric(rcomponents %in% mcomponents), 
# 													object$Random$nrt), object$Random$nfl^2) == 1)
# 		keep <- c(keep, which(rep(rep(rep(1, length(object$Residual$nrt)), 
# 																	object$Residual$nrt), object$Residual$nfl^2) == 1))

# Replacement code
		keep <- rep(rep(as.numeric(rcomponents %in% mcomponents), 
										object$Random$nrt), object$Random$nfl^2) == 1
		keep <- c(keep, rep(rep(rep(1, length(object$Residual$nrt)), 
														object$Residual$nrt), object$Residual$nfl^2) == 1)
		keep <- which(keep)

		if(ncol(vpred) == 1){
			keep <- 1
			}
		
		vpred <- vpred[, keep, drop = FALSE]
		
# Original code: non-conformable matrices when marginal = TRUE
		postvar <- t(apply(object$VCV[, keep, drop = FALSE], 
											 1, function(x) {(vpred %*% x)}))
		
#Replacement code		
# 		postvar <- t(apply(t(object$VCV[, keep, drop = FALSE]), 
# 											 2, function(x) {
# 											 	(vpred[,ifelse(ncol(vpred)>1,keep,1), 
# 											 				 drop = FALSE] %*% x)}))
	}
	
	# Posterior predicitons
	
	post.pred <- t(apply(object$Sol, 1, function(x){(W %*% x)@x}))

	if (interval == "prediction") {
		post.pred <- matrix(rnorm(prod(dim(post.pred)), post.pred, 
															sqrt(postvar)), dim(post.pred)[1], dim(post.pred)[2])
	}

	# Only relevant if predictions are called in the response scale
	
	if (type == "response") {
		if (any(object$family %in% c("poisson", "cenpoisson", 
																 "multinomial", "categorical", "gaussian", "cengaussian", 
																 "ordinal") == FALSE)) {
			stop("sorry - prediction on data scale not implemented for this family")
		}
		if (any(object$family %in% c("poisson", "cenpoisson"))) {
			keep <- which(object$family %in% c("poisson", "cenpoisson"))
			if (interval == "prediction") {
				post.pred[, keep, drop = FALSE] <- exp(post.pred[, keep, drop = FALSE])
			}
			else {
				post.pred[, keep, drop = FALSE] <- exp(post.pred[, keep, drop = FALSE] + 
					0.5 * postvar[, keep, drop = FALSE])
			}
		}
		if (any(object$family %in% c("multinomial", "categorical"))) {
			c2 <- (16 * sqrt(3)/(15 * pi))^2
			keep <- which(object$family %in% c("multinomial", 
																				 "categorical"))
			if (interval == "prediction") {
				post.pred[, keep, drop = FALSE] <- plogis(post.pred[, keep, drop = FALSE])
			}
			else {
				post.pred[, keep] <- plogis(post.pred[, keep, drop = FALSE]/sqrt(1 + 
					c2 * postvar[, keep, drop = FALSE]))
			}
		}
		if (any(object$family %in% c("ordinal"))) {
			keep <- which(object$family %in% c("ordinal"))
			CP <- cbind(-Inf, 0, object$CP, Inf)
			q <- matrix(0, dim(post.pred)[1], length(keep))
			if (interval == "prediction") {
				for (i in 2:(dim(CP)[2] - 1)) {
					q <- q + (pnorm(CP[, i + 1] - post.pred[, keep, drop = FALSE]) - 
						pnorm(CP[, i] - post.pred[, keep, drop = FALSE])) * (i - 
						1)
				}
			}
			else {
				for (i in 2:(dim(CP)[2] - 1)) {
					q <- q + (pnorm(CP[, i + 1] - post.pred[, keep], 
													0, sqrt(postvar[, keep] + 1)) - pnorm(CP[, 
																																	 i] - post.pred[, keep], 0, sqrt(postvar[, 
																																	 																				keep] + 1))) * (i - 1)
				}
			}
			post.pred[, keep] <- q
			rm(q)
		}
	}
	
	# Predictions
	
	pred <- matrix(colMeans(post.pred), dim(post.pred)[2], 1)
	if (!grepl("none|all", interval)) {
		pred <- cbind(pred, coda::HPDinterval(mcmc(post.pred), 
																					prob = level))
		colnames(pred) <- c("fit", "lwr", "upr")
	}
	if(interval == "all"){
		pred <- t(post.pred)
	}
	
	if(!is.null(index)) {
		prednames <- index
	} else {
		prednames <- 1:dim(pred)[1]
	}
		rownames(pred) <- prednames
	return(pred)
}
