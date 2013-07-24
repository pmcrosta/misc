## by peter m crosta: pmcrosta at gmail
### Functions for oaxac decompositions in R given sample means and regression coefficients as inputs.

setwd('~/CSVs')
# read in means file
meansraw <- read.csv('means.csv', header=T, stringsAsFactors=F)

# read in coeffs file
coeffsraw <- read.csv('coeffs.csv', header=T, stringsAsFactors=F)

meansraw <- meansraw[meansraw$X!="",]
coeffsraw <- coeffsraw[coeffsraw$X!="",]

coefnames <- coeffsraw$X
sumnames <- meansraw$X

coefcol <- colnames(coeffsraw)
sumcol <- colnames(meansraw)

oaxaca <- function(mods) UseMethod("oaxaca")

oaxaca.default <- function(mods) {
    ## Input is a charactervector. Likely taken from a list
    mod1 <- mods[1]
    mod2 <- mods[2]
    ### get model coefficients
    m1 <- as.numeric(coeffsraw[, which(mod1==coefcol)])
    m2 <- as.numeric(coeffsraw[, which(mod2==coefcol)])

    names(m1)<-coefnames
    names(m2)<-coefnames

    m1 <- na.omit(m1)
    m2 <- na.omit(m2)

    ## get means
    s1 <- as.numeric(meansraw[, which(mod1==sumcol)])
    s2 <- as.numeric(meansraw[, which(mod2==sumcol)])

    names(s1)<-sumnames
    names(s2)<-sumnames

    s1 <- na.omit(s1)
    s2 <- na.omit(s2)

    if (substr(mod1,1,1)=="r") depvar <- "readdep" else if (substr(mod1,1,1)=="m") depvar <- "mathdep"

    k <- NROW(m1)
    k1 <- NROW(m2)
    if (k != k1) stop("models are not the same size")


    W <- list("W = 1 (Oaxaca, 1973)" = diag(1, k),
              "W = 0 (Blinder, 1973)" = diag(0, k),
              "W = 0.5 (Reimers 1983)" = diag(0.5, k))

    mnames <- names(m1)[-which(names(m1)=="Constant")]
    X1 <- c(s1[mnames], 1)
    X2 <- c(s2[mnames], 1)

    Q <- sapply(W, function(w) crossprod(X1 - X2, w %*% m1 + (diag(1, k) - w) %*% m2))
    U <- sapply(W, function(w) (crossprod(X1, diag(1, k) - w) + crossprod(X2, w)) %*% (m1 - m2))

    attcoeff <- c("_Iattdec_2", "_Iattdec_3", "_Iattdec_4", "_Iattdec_5", "_Iattdec_6", "_Iattdec_7", "_Iattdec_8", "_Iattdec_9", "_Iattdec_10")
    
    #Qa <- sapply(W, function(w) crossprod(X1['lnattend'] - X2['lnattend'], w[1] %*% m1['lnattend'] + (diag(1, 1) - w[1]) %*% m2['lnattend']))
    #Ua <- sapply(W, function(w) (crossprod(X1['lnattend'], diag(1, 1) - w[1]) + crossprod(X2['lnattend'], w[1])) %*% (m1['lnattend'] - m2['lnattend']))

    Qa <- Ua <- 0
    for (ii in attcoeff) {
        Qa <- Qa + sapply(W, function(w) crossprod(X1[ii] - X2[ii], w[1] %*% m1[ii] + (diag(1, 1) - w[1]) %*% m2[ii]))
        Ua <- Ua + sapply(W, function(w) (crossprod(X1[ii], diag(1, 1) - w[1]) + crossprod(X2[ii], w[1])) %*% (m1[ii] - m2[ii]))
    }

    KK<-length(X1)
    Qdetail <- sapply(W, function(w) sapply(1:KK, function(jj) crossprod(X1[jj] - X2[jj], w[1] %*% m1[jj] + (diag(1, 1) - w[1]) %*% m2[jj])))
    Udetail <- sapply(W, function(w) sapply(1:KK, function(jj) (crossprod(X1[jj], diag(1, 1) - w[1]) + crossprod(X2[jj], w[1])) %*% (m1[jj] - m2[jj])))
    rownames(Qdetail) <- rownames(Udetail) <- c(mnames, "Constant")
    
    # Var-covar matrices
    setwd('~/outregs')
    VB1 <- as.matrix(read.table(paste(mod1, 'mat.txt', sep='')))
    VB2 <- as.matrix(read.table(paste(mod2, 'mat.txt', sep='')))

    VX1 <- as.matrix(read.table(paste(mod1, 'cor.txt', sep='')))
    VX2 <- as.matrix(read.table(paste(mod2, 'cor.txt', sep='')))

    dep1 <- s1[depvar]
    dep2 <- s2[depvar]
    Ddep <- dep1-dep2
    
    ans <- list(R = Ddep,
              VR = crossprod(X1, VB1) %*% X1 + crossprod(m1, VX1) %*% m1 + sum(diag(VX1 %*% VB1)) +
              crossprod(X2, VB2) %*% X2 + crossprod(m2, VX2) %*% m2 + sum(diag(VX2 %*% VB2)))

    
    VQ <- sapply(W, function(w)
              sum(diag((VX1 + VX2) %*% (w %*% tcrossprod(VB1, w) + (diag(1, k) - w) %*% tcrossprod(VB2, diag(1, k) - w)))) +
              crossprod(X1 - X2, w %*% tcrossprod(VB1, w) + (diag(1, k) - w) %*% tcrossprod(VB2, diag(1, k) - w)) %*% (X1 - X2) +
              crossprod(w %*% m1 + (diag(1, k) - w) %*% m2, VX1 + VX2) %*% (w %*% m1 + (diag(1, k) - w) %*% m2))
    VU <- sapply(W, function(w)
              sum(diag((crossprod(diag(1, k) - w, VX1) %*% (diag(1, k) - w) + crossprod(w, VX2) %*% w) %*% (VB1 + VB2))) +
              crossprod(crossprod(diag(1, k) - w, X1) + crossprod(w, X2), VB1 + VB2) %*% (crossprod(diag(1, k) - w, X1) + crossprod(w, X2))+
              crossprod(m1 - m2, crossprod(diag(1, k) - w, VX1) %*% (diag(1, k) - w) + crossprod(w, VX2) %*% w) %*% (m1 - m2))

    
    ans$Q <- Q
    ans$U <- U
    ans$VQ <- VQ
    ans$VU <- VU
    ans$W <- W
    ans$Qa <- Qa
    ans$Ua <- Ua
    ans$call <- match.call()
    ans$mod1 <- mod1
    ans$mod2 <- mod2
    ans$Qdetail <- Qdetail
    ans$Udetail <- Udetail
    class(ans) <- "oaxaca"
    ans
    
}


print.oaxaca <- function(x) {
  se <- sqrt(x$VQ)
  zval <- x$Q / se
  decomp <- cbind(Explained = x$Q, StdErr = se, "z-value" = zval, "Pr(>|z|)" = 2*pnorm(-abs(zval)))
  cat("\nBlinder-Oaxaca decomposition\n\nCall:\n")
  print(x$call)
  cat(x$mod1, ' ', x$mod2, "\n")
  decomp <- matrix(nrow = 1, ncol = 4)
  decomp[1, c(1, 2)] <- c(x$R, sqrt(x$VR))
  decomp[1, 3] <- c(decomp[, 1] / decomp[, 2])
  decomp[1, 4] <- c(2 * pnorm(-abs(decomp[, 3])))
  colnames(decomp) <- c("Difference", "StdErr", "z-value", "Pr(>|z|)")
  rownames(decomp) <- "Mean"
  cat("\n")
  print(zapsmall(decomp))
  cat("\nLinear decomposition:\n")
  for (i in 1:3) {
    cat(paste("\nWeight: ", names(x$W[i]), "\n"))
    decomp <- cbind(Difference = c(x$Q[i], x$U[i]), StdErr = c(sqrt(x$VQ[i]), sqrt(x$VU[i])))
    decomp <- cbind(decomp, "z-value" = decomp[, 1] / decomp[, 2])
    decomp <- cbind(decomp, "Pr(>|z|)" = 2 * pnorm(-abs(decomp[, 3])))
    decomp <- cbind(decomp, c(x$Qa[i], x$Ua[i]))
    colnames(decomp) <- c("Difference", "StdErr", "z-value", "Pr(>|z|)", "attendance")
    rownames(decomp) <- c("Explained", "Unexplained")
    print(zapsmall(decomp))
  }
  cat("\n")
}

write.oaxaca <- function(x) {
    if (class(x) != "oaxaca") stop("Input needs to be of class oaxaca")
  #cat(x$mod1, ' ', x$mod2, "\n")
  decomp <- matrix(nrow = 3, ncol = 5)
  decomp[1, c(1, 2)] <- c(x$R, sqrt(x$VR))
  decomp[1, 3] <- c(decomp[1, 1] / decomp[1, 2])
  decomp[1, 4] <- c(2 * pnorm(-abs(decomp[1, 3])))
  colnames(decomp) <- c("Difference", "StdErr", "z-value", "Pr(>|z|)", paste(x$mod1, x$mod2, sep=''))
  if (substr(x$mod1,1,1)=="r") depvar <- "Reading" else if (substr(x$mod1,1,1)=="m") depvar <- "Math"
  rownames(decomp) <- c(depvar, "Explained", "Unexplained")
    i<-3
    decomp[2,1] <- x$Q[i]
    decomp[3,1] <- x$U[i]
    decomp[2,2] <- sqrt(x$VQ[i])
    decomp[3,2] <- sqrt(x$VU[i])
    decomp[2,3] <- decomp[2, 1] / decomp[2, 2]
    decomp[3,3] <- decomp[3, 1] / decomp[3, 2]
    decomp[2,4] <- 2 * pnorm(-abs(decomp[2, 3]))
    decomp[3,4] <- 2 * pnorm(-abs(decomp[3, 3]))
    decomp[2,5] <- x$Qa[i]
    decomp[3,5] <- x$Ua[i]
    zapsmall(decomp)    
}

detail.oaxaca <- function(x, depth=6) {
  if (class(x) != "oaxaca") stop("Input needs to be of class oaxaca")
  ##only deal with Cotton/Reimers
  Qdet <- sort(x$Qdetail[,3])
  Udet <- sort(x$Udetail[,3])
  Qh <- head(Qdet, depth)
  Qt <- tail(Qdet, depth)
  Uh <- head(Udet, depth)
  Ut <- tail(Udet, depth)
  Qsum <- cbind(Qh, names(Qt), Qt)
  Usum <- cbind(Uh, names(Ut), Ut)
  ret.list <- list(Qsum, Usum)
}

### Below is code specific to the analysis
### which decompositions do we want to do?
### coefcol[-1] lists the options:
### math, reading
### pooled, fixed effects
### racedum1-5, econdum1-4
### remember: race: amerindian, asian, black, hispanic, white
### remember: econ: not poor, free lunch, reduced lunch, other dis 

## For pooled and FE models: white-black, white-hispanic, hispanic-black
##              notpoor-free, notpoor-reduced, notpoor-other

## This becomes: 5-3, 5-4, 4-3; 1-2, 1-3, 1-4
grouprace <- list(
 c('mpa4racedum5C', 'mpa4racedum3C'),
 c('rpa4racedum5C', 'rpa4racedum3C'),
 c('mfa4racedum5C', 'mfa4racedum3C'),
 c('rfa4racedum5C', 'rfa4racedum3C'),

 c('mpa4racedum5C', 'mpa4racedum4C'),
 c('rpa4racedum5C', 'rpa4racedum4C'),
 c('mfa4racedum5C', 'mfa4racedum4C'),
 c('rfa4racedum5C', 'rfa4racedum4C'),

 c('mpa4racedum4C', 'mpa4racedum3C'),
 c('rpa4racedum4C', 'rpa4racedum3C'),
 c('mfa4racedum4C', 'mfa4racedum3C'),
 c('rfa4racedum4C', 'rfa4racedum3C'))

groupecon <- list(
 c('mpa4econdum1C', 'mpa4econdum2C'),
 c('rpa4econdum1C', 'rpa4econdum2C'),
 c('mfa4econdum1C', 'mfa4econdum2C'),
 c('rfa4econdum1C', 'rfa4econdum2C'),

 c('mpa4econdum1C', 'mpa4econdum3C'),
 c('rpa4econdum1C', 'rpa4econdum3C'),
 c('mfa4econdum1C', 'mfa4econdum3C'),
 c('rfa4econdum1C', 'rfa4econdum3C'),

 c('mpa4econdum1C', 'mpa4econdum4C'),
 c('rpa4econdum1C', 'rpa4econdum4C'),
 c('mfa4econdum1C', 'mfa4econdum4C'),
 c('rfa4econdum1C', 'rfa4econdum4C'))

oax.race<-lapply(grouprace, oaxaca)
oax.econ<-lapply(groupecon, oaxaca)
 
## oaxaca print out function
lapply(oax.race, write.oaxaca)
lapply(oax.econ, write.oaxaca)

## print out detailed decomposition
race.det <- lapply(oax.race, detail.oaxaca)
names(race.det) <- sapply(grouprace, function(x) paste(x, collapse=''))
race.det

econ.det <- lapply(oax.econ, detail.oaxaca)
names(econ.det) <- sapply(groupecon, function(x) paste(x, collapse=''))
econ.det
