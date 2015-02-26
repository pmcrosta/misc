#######################
# INPUT DATASETS: ~/projects/virginia/original_data/Aug2010_Transfer/stata
###############################################

## working on chagall with rstudio
library(foreign)
library(lsa)
library(irlba)
library(Cairo)
library(igraph)
library(descr)
library(R.utils)
library(snow)
library(maptools)
library(plyr)
library(cluster)
library(ggdendro)
library(Hmisc)
library(sda)
library(scatterplot3d)

lunique <- function(x) length(unique(x))
pdensity <- function(x) plot(density(x, na.rm=T))
cum.na <- function(x) {
  ## cumulative sum function that counts NA as 0 in sum
  y <- x
  x[which(is.na(x))] <- 0 
  x <- cumsum(x)
  x[which(is.na(y))] <- NA
  return(x)
} 

initialload <- function() {
  setwd("~/projects/virginia/original_data/Aug2010_Transfer/stata/")
  
  ## read in stata data files so they are in local R environment
  ## years go from 04 to 09
  ## types are classes, fa, graduate_awards, nsc_transfers, recommendations, student_profiles
  
  yearvec <- paste(0, 4:9, sep="")
  filetype <- c("classes", "fa", "graduate_awards", "nsc_transfers", "recommendations", "student_profiles")  
  
  if (!exists("allcourse")) {
    for (ii in yearvec) {
      for (jj in filetype) {
        assign(paste("y", ii, "_", jj, sep=''), read.dta(paste("y", ii, "_", jj, ".dta", sep='')))
      }
    }
  }
  
  
  ## append course data frames for all represented years
  appcom <- paste("allcourse <- rbind(", paste("y", yearvec, "_classes ", collapse=',', sep=''), ")", sep='')
  eval(parse(text=appcom))
  rm(list=unlist(paste("y", yearvec, "_classes", sep='')))
  
  ## append other data frames.
  appcom <- paste("allfa <- rbind(", paste("y", yearvec, "_fa ", collapse=',', sep=''), ")", sep='')
  eval(parse(text=appcom))
  rm(list=unlist(paste("y", yearvec, "_fa", sep='')))
  
  appcom <- paste("allgrad <- rbind(", paste("y", yearvec, "_graduate_awards ", collapse=',', sep=''), ")", sep='')
  eval(parse(text=appcom))
  rm(list=unlist(paste("y", yearvec, "_graduate_awards", sep='')))
  
  appcom <- paste("allnsc <- rbind(", paste("y", yearvec, "_nsc_transfers ", collapse=',', sep=''), ")", sep='')
  eval(parse(text=appcom))
  rm(list=unlist(paste("y", yearvec, "_nsc_transfers", sep='')))
  
  appcom <- paste("allrecom <- rbind(", paste("y", yearvec, "_recommendations ", collapse=',', sep=''), ")", sep='')
  eval(parse(text=appcom))
  rm(list=unlist(paste("y", yearvec, "_recommendations", sep='')))
  
  appcom <- paste("allstudprof <- rbind(", paste("y", yearvec, "_student_profiles ", collapse=',', sep=''), ")", sep='')
  eval(parse(text=appcom))
  rm(list=unlist(paste("y", yearvec, "_student_profiles", sep='')))
  
  rm(appcom, ii, jj)
  gc()
  save.image("Aug2010_Transfer.Rdata")
}

setwd("~/projects/virginia/original_data/Aug2010_Transfer/Rdata/")

if (!exists("allcourse")) load("Aug2010_Transfer_Mod.Rdata")

if (!exists("allcourse")) {
  load("Aug2010_Transfer.Rdata")
  
  ##Looks like there are 161,343 student observations
  ## in course file, identify cohort of student
  allcourse$cohort <- substr(allcourse$id, 1, 2)
  
  ## create factor of term variable, in order
  termvec <- c("SU00", "FA00", paste(rep(c("SP", "SU", "FA"),9), 0, sapply(1:9, rep, 3), sep=""), "SP10", "SU10", "FA10")
  allcourse$term.fac <- factor(allcourse$term, levels=termvec)
  
  ## create cleaner grade variable.
  allcourse$grade2 <- substr(allcourse$grade, 1,1)
  allcourse$grade <- NULL
  
  # A, B, C, D, F get numbers for GPA calculation
  allcourse$grade.num <- NA
  allcourse$grade.num[allcourse$grade2=="A"] <- 4
  allcourse$grade.num[allcourse$grade2=="B"] <- 3
  allcourse$grade.num[allcourse$grade2=="C"] <- 2
  allcourse$grade.num[allcourse$grade2=="D"] <- 1
  allcourse$grade.num[allcourse$grade2=="F"] <- 0
  
  ## pass/fail boolean
  allcourse$pass <- NA
  allcourse$pass[allcourse$grade2 %in% c("A", "B", "C", "D", "S", "P")] <- TRUE
  allcourse$pass[allcourse$grade2 %in% c("F", "I", "N", "U", "W")] <- FALSE
  
  ## remove duplicate courses. 
  allcourse.sort <- allcourse[order(allcourse$id, allcourse$course, -as.numeric(allcourse$term.fac)),]
  duplist3 <- duplicated(allcourse.sort[,c(1,6)]) # id, course
  allcourse <- allcourse.sort
  faclevels <- levels(allcourse$term.fac)
  rm(allcourse.sort)
  gc()
  save.image("Aug2010_Transfer_Mod.Rdata")
}

######################### SVD ################################
## this section uses SVD to build and test a model of determining when a student
## is on a course of study. The basic idea is that the student will become more and more
## similar to a given completer who has completed a certain program. Eventually over time,
## the probability will be high enough to claim an area of concentration, and thus a program of study.
## The analysis is built on an SVD of the course-id matrix from completers. We will determine
## appropriate thresholds for the dissimilarity metric and the metric itself
## by training the model on portion of completers and testing with another portion.

## individual courses in a binary incidence matrix.
## make factors out of these so tables all have the same dimensions
allgrad$award <- factor(allgrad$award)
allgrad$award_curriculum <- factor(allgrad$award_curriculum)
allgrad$cip <- factor(allgrad$cip)
allcourse$hegis_name <- factor(allcourse$hegis_name)
allcourse$hegis_code <- factor(allcourse$hegis_code)
allcourse$program_major <- factor(allcourse$program_major)

## remove CSC from graduate data set
allgrad <- subset(allgrad, award != "CSC")

## to get a cleaner course X transcript matrix, we are only going to use those who
## completed a single program of study
dups <- allgrad$id[duplicated(allgrad$id)]
allgrad <- subset(allgrad, !(id %in% dups))

## really need to combine award with award_curriculum.
allgrad$progofstudy <- gsub(pattern="[0-9]",replacement="",paste0(allgrad$award, allgrad$award_curriculum))

## First step: determine which completers to use, drawn from those 
completers <- unique(allgrad$id)

## Create the cross-tabulation matrix for completers only. wmat2 is a binary course incidence matrix. 
load("wmat2.Rdata")
if (!exists("wmat2")) wmat2 <- xtabs(pass ~ course + id, data=allcourse[!duplist3 & allcourse$id %in% completers,])

## another binary incidence matrix of college that awarded degree
collegemat <- table(factor(allgrad$college),allgrad$id)
wmat2.bk <- wmat2
wmat2 <- rbind(wmat2.bk,collegemat)

## decision to use 5-fold cross validation.
## divide sample into five groups
set.seed(7402313)

# Generating random indices 
kk <- 5 # number of folds
group <- sample(rep(seq_len(kk), length.out=length(completers)))

# lapply over them:
indicies <- lapply(seq_len(kk), function(a) list(
  test_matrix_indices = which(group==a),
  train_matrix_indices = which(group!=a)
))
#str(indicies)

## Weighting options:
#  1) Use wmat2 as it stands
#  2) Use global weight idf (gw_idf)
#  3) Use global weight idf adjusted for global freq (gw_gfidf)
#  4) Use gw_entropy global weightas.integer(trainpct*length(completers))

GW_IDF <- function (m) 
{
  df = rowSums(lw_bintf(m), na.rm = TRUE)+1
  return((log2(ncol(m)/df) + 1))
}
GW_GFIDF <- function (m) 
{
  gf = rowSums(m, na.rm = TRUE)
  df = rowSums(lw_bintf(m), na.rm = TRUE)+1
  return(gf/df)
}

wtfunc <- function(inp, wt) {
  if (wt==1) {
    return(inp)
  } else if (wt==2) {
    return(inp*gwidf)
  } else if (wt==3) {
    return(inp*gwgfidf)
  } else if (wt==4) {
    return(inp*gwentropy)
  }
}

gwidf <- GW_IDF(wmat2)
#gwgfidf <- GW_GFIDF(wmat2)
#gwentropy <- gw_entropy(wmat2)

### This is pretty much all the preprocessing stuff. The rest is modeling and eval
gc()
save.image("ready_for_svd.Rdata")

setwd("~/projects/virginia/original_data/Aug2010_Transfer/Rdata/")
load("ready_for_svd.Rdata")

## &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& ##
# choose weight (1 or 2), build weighted matrix

wtscheme <- 2 #1 or 2
wmat3 <- wtfunc(wmat2, wtscheme)

## compute SVD
load(paste("svd-wt", wtscheme, ".Rdata", sep=''))
if (!exists("sings")) {
  sings <- vector("list", 5)
  for (ii in 1:5) {
    sings[[ii]] <- svd(wmat3[,indicies[[ii]]$train_matrix_indices]) ## 64 seconds
  }
  save(sings, file=paste("svd-wt", wtscheme, "-inst.Rdata", sep=''))
}

dims <- c(100, 250, 500, 750, 1000, 1500, 2000, 2500, dim(wmat3)[1])

## initialize results list of matrices
reslist <- vector("list", length(dims))
names(reslist) <- dims

timestamp()
for (ndim in dims) { ## Looping over dimensions
  reslist[[as.character(ndim)]] <- vector("list", length=5)
  reslist[[as.character(ndim)]][[1]] <- matrix(data=NA, nrow=5, ncol=9, dimnames=list(1:5, c("recall", "wtrecall", "precision", "wtprecision", "f1", "wtf1", "sens1", "sens2", "sens3")))
  reslist[[as.character(ndim)]][[2]] <- vector("list", length=5)
  reslist[[as.character(ndim)]][[3]] <- vector("list", length=5)
  reslist[[as.character(ndim)]][[4]] <- vector("list", length=5)
  reslist[[as.character(ndim)]][[5]] <- vector("list", length=5)
  for (ii in 1:5) { ## Looping over validation folds ####
    ## pull out orthogonal SVD matrices and vector of singular values
    Uk <- sings[[ii]]$u[,1:ndim]
    Vk <- sings[[ii]]$v[,1:ndim]
    Sk <- sings[[ii]]$d[1:ndim]
    SkVk <- diag(Sk) %*% t(Vk)
    rownames(Uk) <- rownames(wmat3)
    colnames(SkVk) <- colnames(wmat3[,indicies[[ii]]$train_matrix_indices])
    
    ## store student ids of current train and test set
    train <- allgrad$id[indicies[[ii]]$train_matrix_indices]
    test <- allgrad$id[indicies[[ii]]$test_matrix_indices]
    
    ## tabulate train and test award curricula
    trainawards <- table(as.character(allgrad$progofstudy[indicies[[ii]]$train_matrix_indices]))
    testawards <- table(as.character(allgrad$progofstudy[indicies[[ii]]$test_matrix_indices]))
    awardcur <- as.character(allgrad$progofstudy[indicies[[ii]]$test_matrix_indices])
    
    ## compute centroids for awards and variance
    award.centroid <- centroids(t(SkVk), allgrad$progofstudy[indicies[[ii]]$train_matrix_indices], lambda.var=0, var.groups=TRUE, verbose=F)
    award.centroid$means <- award.centroid$means[,1:length(trainawards)]
    
    # sort(colSums(award.centroid$variances)) ## variance of centroid  
    
    ## create default query, all entries zero; transcript length matrix (num courses); top courses by program major; distance to award_curriculum medioids 
    q.default <- matrix(data=0, nrow=nrow(wmat3), ncol=1, dimnames=list(rownames(wmat3), "id"))
    translenmat <- matrix(data=0, nrow=length(test), ncol=1, dimnames=list(test, "ncourse"))
    toplist <- vector(mode="list", length=length(test))
    names(toplist) <- test
    
    ## create temporary course matrix for test set
    tmpcour <- subset(allcourse, subset=!duplist3 & allcourse$id %in% test & allcourse$pass, select=c("id", "course", "program_major"))
    
    loopdex <-na.omit(unique(tmpcour$id)) #same as test generally when using all 32 terms  
    ## Looping over test set members to compute distances ####
    for (jj in 1:length(loopdex)) {
      vicki <- loopdex[jj]
      vicki.q <- q.default
      partialcourse <- as.character(na.omit(unique(tmpcour$course[tmpcour$id==vicki])))
      
      ## ignore observation if student has no passed courses in the term
      if (length(partialcourse)==0) next		
      
      ## ignore observation if no new information added in this term
      translenmat[vicki,1] <- length(partialcourse)
      
      ## generate base query
      parttran <-	table(partialcourse)
      vicki.q[partialcourse,1] <- parttran
      
      ## apply weighting scheme to query vector
      vicki.q <-  wtfunc(vicki.q, wtscheme)
      
      ## project partial transcript into course space
      #	and/or use distance metric to determine proximity to completers		
      vicki.q <- t(Uk) %*% vicki.q
      simmatsvd <- sapply(1:dim(award.centroid$means)[2], function(x) cosine(as.vector(vicki.q), award.centroid$means[,x]))
      names(simmatsvd) <- colnames(award.centroid$means)
      
      # skip if simmatsvd is all Nan
      if (all(is.nan(simmatsvd))) next
      
      toplist[[vicki]] <- sort(simmatsvd, decreasing=T)      
    } ## End loop over test set members ####
    
    ## pull out top 1, 2, 3 matches
    lasttermmed <- sapply(toplist, function(x) unlist(x)[1])
    lasttermmed2 <- lapply(toplist, function(x) unlist(x)[c(1,2)])
    lasttermmed3 <- lapply(toplist, function(x) unlist(x)[c(1,2,3)])
    
    ## identify predicted class; establish true positives for each student
    firststring <- sapply(strsplit(names(lasttermmed),split='\\.'), function(x) x[[2]])
    TP <- awardcur==firststring
    
    ## proportions on these will be aggregate sensitivity@n numbers
    predmatch <- sapply(1:length(test), function(x) substring(names(lasttermmed)[x], first=nchar(test[1])+2) == awardcur[x])
    predmatch2 <- sapply(1:length(test), function(x) any(names(lasttermmed2[[x]]) %in% awardcur[x]))
    predmatch3 <- sapply(1:length(test), function(x) any(names(lasttermmed3[[x]]) %in% awardcur[x]))
    
    ## create the matrix used to compute precision and recall
    TPmat <- sapply(names(testawards), function(x) sum(TP[awardcur==x]))
    evalmat <- cbind(TPmat, testawards)
    evalmat <- merge(evalmat, data.frame(table(firststring), stringsAsFactors=F), all.x=T, by.x="row.names", by.y="firststring")
    evalmat$Freq[is.na(evalmat$Freq)] <-0
    evalmat$precision <- evalmat$TPmat/evalmat$Freq
    evalmat$precision[is.nan(evalmat$precision)] <- 1
    evalmat$recall <- evalmat$TPmat/evalmat$testawards
    
    ## compute average over classes precision, recall, weighted, f1, sensitivity
    recall <- mean(evalmat$recall, na.rm=T)
    wtrecall <- wtd.mean(evalmat$recall, weights=evalmat$testawards)
    precision <- mean(evalmat$precision, na.rm=T)
    wtprecision <- wtd.mean(evalmat$precision, weights=evalmat$testawards)
    f1 <- (2*recall*precision)/(recall+precision)
    wtf1 <- (2*wtrecall*wtprecision)/(wtrecall+wtprecision)
    
    sens1 <- prop.table(table(predmatch))["TRUE"]
    sens2 <- prop.table(table(predmatch2))["TRUE"]
    sens3 <- prop.table(table(predmatch3))["TRUE"]
    
    reslist[[as.character(ndim)]][[1]][ii,] <- c(recall, wtrecall, precision, wtprecision, f1, wtf1, sens1, sens2, sens3)
    reslist[[as.character(ndim)]][[2]][[ii]] <- award.centroid
    reslist[[as.character(ndim)]][[3]][[ii]] <- translenmat
    reslist[[as.character(ndim)]][[4]][[ii]] <- toplist
    reslist[[as.character(ndim)]][[5]][[ii]] <- evalmat
    ##recall is TPR or sensitivity. precision is positive predicted value 
  } ## End loop over validation folds ####
  timestamp()  
} ## End loop over dimensions ####
timestamp()
#close(pb)

## Save and clean up ####
setwd("~/working/virginia/code/")
save(reslist, sings, file=paste("course-wt", wtscheme, "-svd-inst", ".Rdata", sep=''))
rm(SkVk, Uk, Vk, evalmat, q.default, tmpcour, translenmat, vicki.q, Sk, TP, TPmat, award.centroid,
   awardcur, centtest, completers, conttab, coursedistlist, dims, f1, firststring,
   ii, jj, kk, lasttermmed, lasttermmed2, lasttermmed3, loopdex, ndim, partialcourse,
   parttran, pb, precision, predmatch, predmatch2, predmatch3, recall,
   sens1, sens2, sens3, simmatsvd, test, testawards, toplist, train, trainawards,
   vicki, wtf1, wtprecision, wtrecall)

## to get variances???
#sort(colSums(reslist[["3549"]][[2]][[1]]$variances))


########## SUMMARIZE RESULTS ##############
## bring in result data sets
setwd("~/working/virginia/code/")

wtscheme <- 2 #1 or 2
load(paste("course-wt", wtscheme, "-svd", ".Rdata", sep=''))

## output summary statistics. Average over the folds for each SVD
outmat <- do.call("rbind", lapply(lapply(reslist, function(x) x[[1]]), function(y) colMeans(y)))
outmat

##### Partial Transcripts ######
wtshceme <- 2
ndim <- 1000
wmat3 <- wtfunc(wmat2, wtscheme)

## initialize results list of matrices
partreslist <- vector("list", length=5)
names(partreslist) <- 1:5

for (ii in 1:5) 
{ ## Looping over validation folds ####
  timestamp()
  Uk <- sings[[ii]]$u[,1:ndim]
  Vk <- sings[[ii]]$v[,1:ndim]
  Sk <- sings[[ii]]$d[1:ndim]
  SkVk <- diag(Sk) %*% t(Vk)
  rownames(Uk) <- rownames(wmat3)
  colnames(SkVk) <- colnames(wmat3[,indicies[[ii]]$train_matrix_indices])
  
  award.centroid <- reslist[[as.character(ndim)]][[2]][[ii]]
  tranlength <- reslist[[as.character(ndim)]][[3]][[ii]] 
  
  ## store student ids of current train and test set
  train <- allgrad$id[indicies[[ii]]$train_matrix_indices]
  test <- allgrad$id[indicies[[ii]]$test_matrix_indices]
  
  ## tabulate train and test award curricula
  trainawards <- table(as.character(allgrad$progofstudy[indicies[[ii]]$train_matrix_indices]))
  testawards <- table(as.character(allgrad$progofstudy[indicies[[ii]]$test_matrix_indices]))
  awardcur <- as.character(allgrad$progofstudy[indicies[[ii]]$test_matrix_indices])
  
  ## create default query, all entries zero; transcript length matrix (num courses); top courses by program major; distance to award_curriculum medioids 
  q.default <- matrix(data=0, nrow=nrow(wmat3), ncol=1, dimnames=list(rownames(wmat3), "id"))
  translenmat <- matrix(data=0, nrow=length(test), ncol=length(faclevels), dimnames=list(test, paste(1:length(faclevels))))
  toplist <- vector(mode="list", length=length(test))
  names(toplist) <- test
  
  for (kk in 1:length(faclevels)) {
    tmpcour <- subset(allcourse, subset=!duplist3 & allcourse$term.fac %in% faclevels[1:kk] & allcourse$id %in% test & allcourse$pass & !is.na(allcourse$pass), select=c("id", "course", "program_major", "term.fac"))
    loopdex <-na.omit(unique(tmpcour$id))
  
    ## Looping over test set members to compute distances ####
    for (jj in 1:length(loopdex)) {
      vicki <- loopdex[jj]
      vicki.q <- q.default
      partialcourse <- as.character(na.omit(unique(tmpcour$course[tmpcour$id==vicki])))
      
      ## ignore observation if student has no passed courses in the term
      if (length(partialcourse)==0) next  	
      
      ## generate base query
      parttran <-	table(partialcourse)
      vicki.q[partialcourse,1] <- parttran
      
      ## ignore observation if no new information added in this term
      translenmat[vicki,kk] <- length(partialcourse)
      if (kk>1) {
        if (translenmat[vicki,kk]==translenmat[vicki,kk-1]) next
      }
      
      ## apply weighting scheme to query vector
      vicki.q <-  wtfunc(vicki.q, wtscheme)
      
      ## project partial transcript into course space
      #	and/or use distance metric to determine proximity to completers		
      vicki.q <- t(Uk) %*% vicki.q
      simmatsvd <- sapply(1:dim(award.centroid$means)[2], function(x) cosine(as.vector(vicki.q), award.centroid$means[,x]))
      names(simmatsvd) <- colnames(award.centroid$means)
      
      # skip if simmatsvd is all Nan
      if (all(is.nan(simmatsvd))) next
      
      toplist[[vicki]][[kk]] <- sort(simmatsvd, decreasing=T)      
    } ## End loop over test set members ####
  } ## End loop over fac levels ####
  partreslist[[ii]][[1]] <- toplist ## list of test, each having 1 to 32 lists of top programs
  partreslist[[ii]][[2]] <- translenmat ## matrix of test x 32
  partreslist[[ii]][[3]] <- matrix(data=NA, nrow=max(translenmat, na.rm=T), ncol=9, dimnames=list(1:max(translenmat, na.rm=T), c("recall", "wtrecall", "precision", "wtprecision", "f1", "wtf1", "sens1", "sens2", "sens3")))
  
  tmpcour$one <- 1
  tranlengthterm <- with(tmpcour, tapply(one, INDEX=list(id, term.fac), FUN=sum, na.rm=T))
  tranlengthtermcum <- t(apply(tranlengthterm, MARGIN=1, FUN=cum.na)) 
  
  ## COMPUTATIONS for evaluation metrics. do for each number of courses 
  for (mm in 1:max(tranlengthtermcum, na.rm=T)) {
    ## a lot of code to pull out the most recent max highest term only
    who <- which(tranlengthtermcum>=1 & tranlengthtermcum<=mm , arr.ind=T)
    who.reorder <- who[order(who[,"row"], -who[,"col"]),]
    who.dup <- duplicated(who.reorder[,"row"])
    who <- who.reorder[!who.dup,]
    if (dim(who)[1]==0) next
    
    notnull <- unlist(sapply(1:nrow(who), function(x) !is.null(toplist[[who[x,"row"]]][[who[x,"col"]]][1])))
    who <- who[notnull,]
    awardcur.sub <- awardcur[who[,"row"]]
    testawards.sub <- table(awardcur.sub)
    
    ## pull out top 1, 2, 3 matches
    lasttermmed <- sapply(1:nrow(who), function(x) toplist[[who[x,"row"]]][[who[x,"col"]]][1])
    lasttermmed2 <- lapply(1:nrow(who), function(x) toplist[[who[x,"row"]]][[who[x,"col"]]][c(1,2)])
    lasttermmed3 <- lapply(1:nrow(who), function(x) toplist[[who[x,"row"]]][[who[x,"col"]]][c(1,2,3)])
    
    ## identify predicted class; establish true positives for each student
    firststring <- names(lasttermmed)
    TP <- awardcur.sub==firststring
    
    ## proportions on these will be aggregate recall@n numbers
    predmatch <- sapply(1:nrow(who), function(x) firststring[x] == awardcur.sub[x])
    predmatch2 <- sapply(1:nrow(who), function(x) any(names(lasttermmed2[[x]]) %in% awardcur.sub[x]))
    predmatch3 <- sapply(1:nrow(who), function(x) any(names(lasttermmed3[[x]]) %in% awardcur.sub[x]))
    
    ## create the matrix used to compute precision and recall
    TPmat <- sapply(names(testawards.sub), function(x) sum(TP[awardcur.sub==x]))
    evalmat <- cbind(TPmat, testawards.sub)
    evalmat <- merge(evalmat, data.frame(table(firststring), stringsAsFactors=F), all.x=T, by.x="row.names", by.y="firststring")
    evalmat$Freq[is.na(evalmat$Freq)] <-0
    #<- evalmat$TPmat[is.na(evalmat$TPmat)] <- evalmat$testawards.sub[is.na(evalmat$testawards.sub)] <- 0
    evalmat$precision <- evalmat$TPmat/evalmat$Freq
    evalmat$precision[is.nan(evalmat$precision)] <- 1
    evalmat$recall <- evalmat$TPmat/evalmat$testawards
    
    ## compute average over classes precision, recall, weighted, f1, sensitivity
    recall <- mean(evalmat$recall, na.rm=T)
    wtrecall <- wtd.mean(evalmat$recall, weights=evalmat$testawards.sub)
    precision <- mean(evalmat$precision, na.rm=T)
    wtprecision <- wtd.mean(evalmat$precision, weights=evalmat$testawards.sub)
    f1 <- (2*recall*precision)/(recall+precision)
    wtf1 <- (2*wtrecall*wtprecision)/(wtrecall+wtprecision)
    
    sens1 <- prop.table(table(predmatch))["TRUE"]
    sens2 <- prop.table(table(predmatch2))["TRUE"]
    sens3 <- prop.table(table(predmatch3))["TRUE"]
    
    partreslist[[ii]][[3]][mm,] <- c(recall, wtrecall, precision, wtprecision, f1, wtf1, sens1, sens2, sens3)    
  }
  timestamp()
} ## End loop over folds ####

## need to save toplist for test set in fold, translenmat
save(partreslist, file="course-wt2-svd100-partcourse.Rdata") 

## summarize partial results over folds
partmean <- (partreslist[[1]][[3]][1:46,] + partreslist[[2]][[3]][1:46,] + partreslist[[3]][[3]][1:46,] + partreslist[[4]][[3]][1:46,] + partreslist[[5]][[3]][1:46,])/5
partmean     

## average compelte transcript length
median(reslist[[as.character(ndim)]][[3]][[1]])

## Section to focus on how classifier gets confused ####
setwd("~/projects/virginia/original_data/Aug2010_Transfer/Rdata/")
load("ready_for_svd.Rdata")
wtscheme <- 2 #1 or 2
wmat3 <- wtfunc(wmat2, wtscheme)
load(paste("svd-wt", wtscheme, ".Rdata", sep=''))
setwd("~/working/virginia/code/")
load(paste("course-wt", wtscheme, "-svd", ".Rdata", sep=''))


ii <- 1 ## first validation fold
ndim <- 1000
Uk <- sings[[ii]]$u[,1:ndim]
Vk <- sings[[ii]]$v[,1:ndim]
Sk <- sings[[ii]]$d[1:ndim]
SkVk <- diag(Sk) %*% t(Vk)
rownames(Uk) <- rownames(wmat3)
colnames(SkVk) <- colnames(wmat3[,indicies[[ii]]$train_matrix_indices])
award.centroid <- reslist[[as.character(ndim)]][[2]][[ii]]
tranlength <- reslist[[as.character(ndim)]][[3]][[ii]] 
## store student ids of current train and test set
train <- allgrad$id[indicies[[ii]]$train_matrix_indices]
test <- allgrad$id[indicies[[ii]]$test_matrix_indices]
## tabulate train and test award curricula
trainawards <- table(as.character(allgrad$progofstudy[indicies[[ii]]$train_matrix_indices]))
testawards <- table(as.character(allgrad$progofstudy[indicies[[ii]]$test_matrix_indices]))
awardcur <- as.character(allgrad$progofstudy[indicies[[ii]]$test_matrix_indices])
awardcurtrain <- as.character(allgrad$progofstudy[indicies[[ii]]$train_matrix_indices])

## here is the evalmat for the first fold (ii).
reslist[[as.character(ndim)]][[5]][[ii]]

## Let's poke into a few programs to see exactly how they are misclassified
toplist <- reslist[[as.character(ndim)]][[4]][[ii]]
lasttermmed <- sapply(toplist, function(x) unlist(x)[1])
firststring <- sapply(strsplit(names(lasttermmed),split='\\.'), function(x) x[[2]])
 
## recall checks
cbind(sort(table(firststring[awardcur=="881/Science"]), decreasing=T))
cbind(sort(table(firststring[awardcur=="648/Liberal Arts"]), decreasing=T))
#cbind(sort(table(firststring[awardcur=="882/Social Sciences"]), decreasing=T))
cbind(sort(table(firststring[awardcur=="831/Engineering"]), decreasing=T))
cbind(sort(table(firststring[awardcur=="909/Automotive"]), decreasing=T))
cbind(sort(table(firststring[awardcur=="290/Information Technology"]), decreasing=T))

## precision checks
cbind(sort(table(awardcur[firststring=="118/Dental Hygiene"]), decreasing=T))
cbind(sort(table(awardcur[firststring=="831/Engineering"]), decreasing=T))
cbind(sort(table(awardcur[firststring=="882/Social Sciences"]), decreasing=T))
cbind(sort(table(awardcur[firststring=="290/Information Technology"]), decreasing=T))








cbind(table(allgrad$progofstudy[allgrad$id %in% unmatch1]), table(allgrad$progofstudy[allgrad$id %in% unmatch2]), table(allgrad$progofstudy[allgrad$id %in% unmatch3]), table(allgrad$progofstudy[allgrad$id %in% test]))

incorpredFP <- sapply(unique(lt_mrg$progofstudy), function(x) sort(table(lt_mrg$award_pred[lt_mrg$progofstudy==x]), decreasing=T))
xx <- lapply(lapply(incorpredFP, function(x) paste(names(x), " (", x, ")", sep="")), function(y) paste(y, collapse=", "))
xx[unlist(lapply(xx, function(x) length(x) > 0))]

## prettified list of award_curr and ways they were misclassified.
paste(names(unlist(xx[unlist(lapply(xx, function(x) length(x) > 0))])),paste(unlist(xx[unlist(lapply(xx, function(x) length(x) > 0))])), sep=": ")


## look at top two predictions
ag_onetest <- allgrad[allgrad$id %in% intersect(onetest, unmatch2), c("id", "progofstudy")]
lt_onetest <- data.frame(dist=lasttermmed[which(test %in% intersect(onetest, unmatch2))])
lt_onetest$id <- substring(rownames(lt_onetest), first=1, last=12)
lt_onetest$award_pred <- substring(rownames(lt_onetest), first=nchar(test[1])+2)
lt_mrg <- merge(ag_onetest, lt_onetest, by="id")
lt_mrg$progofstudy <- as.character(lt_mrg$progofstudy)

incorpredFP <- sapply(unique(lt_mrg$progofstudy), function(x) sort(table(lt_mrg$award_pred[lt_mrg$progofstudy==x]), decreasing=T))
xx <- lapply(lapply(incorpredFP, function(x) paste(names(x), " (", x, ")", sep="")), function(y) paste(y, collapse=", "))
xx[unlist(lapply(xx, function(x) length(x) > 0))]

## prettified list of award_curr and ways they were misclassified.
paste(names(unlist(xx[unlist(lapply(xx, function(x) length(x) > 0))])),paste(unlist(xx[unlist(lapply(xx, function(x) length(x) > 0))])), sep=": ")


progvec <- allgrad$progofstudy
names(progvec) <- allgrad$award_curriculum

newcentroidnames <- progvec[dimnames(award.centroid$means)[[2]]]

dimnames(award.centroid$means) <- list(NULL, newcentroidnames)

## hierarchical cluster of award medioids; dendrograms; will prob split up for presentation
centroid.dist <- dist(t(award.centroid$means))
centroid.hcl <- hclust(centroid.dist, method="ward")
groups <- cutree(centroid.hcl, k=174)
plot(centroid.hcl, labels=substring(colnames(award.centroid$means), 5), cex=.75, xlab="", ylab="",ann=F)
ggdendrogram(centroid.hcl, rotate=T, size=1.62)


  
##plot first and second dimension.
plot(award.centroid$means[1,], award.centroid$means[2,], type="n", xlab="Dimension 1", ylab="Dimension 2", main="")
xyc <- pointLabel(award.centroid$means[1,], award.centroid$means[2,], labels=substring(colnames(award.centroid$means), 5), cex=.75, offset=0)
ex1 <- names(head(sort(award.centroid$means[2,]), 10))
ex2 <- names(tail(sort(award.centroid$means[2,]), 12))
# highlight some 2nd dimension clusters
lines(predict(ellipsoidhull(as.matrix(data.frame(xyc)[which(colnames(award.centroid$means) %in% ex1),]))), col="red")
lines(predict(ellipsoidhull(as.matrix(data.frame(xyc)[which(colnames(award.centroid$means) %in% ex2),]))), col="red")

## Plot of d1 and d2
d1 <- 2
d2 <- 3
plot(award.centroid$means[d1,], award.centroid$means[d2,], type="n", xlab=paste("Dimension",d1), ylab=paste("Dimension",d2), main="")
xyc <- pointLabel(award.centroid$means[d1,], award.centroid$means[d2,], labels=substring(colnames(award.centroid$means), 5), cex=.75)
ex1 <- names(head(sort(award.centroid$means[d1,]), 9))
ex2 <- names(tail(sort(award.centroid$means[d1,]), 7))
# highlight some 10th and 20th dimension clusters
lines(predict(ellipsoidhull(as.matrix(data.frame(xyc)[which(colnames(award.centroid$means) %in% ex1),]))), col="red")


 


  
 
  
  ## other idea is to standardize toplist?
  
  lasttermscale1 <- sapply(test, function(chuck) { 
    scale(toplist[[chuck]][[length(toplist[[chuck]])]])[1]
  })  
  
  # density distributions showing strength of relationship between unmatched and matched lastterm's max scaled score
  # thepoint is that correct matchers have a stronger relationship (compared to the other program medioids) than
  # studends who were incorrectly matched. This was subsetted to only those with one program and only
  # considers the top match (as opposed to the top two or three)
  plot(density(lasttermscale1[which(test %in% intersect(onetest, unmatch1))]), xlim=c(0,15), ylim=c(0,.5), main="", lty=2, axes=F,  xlab="", ylab="")
  par(new=T)
  plot(density(lasttermscale1[which(test %in% intersect(onetest, test[predmatch]))]), xlim=c(0,15), ylim=c(0,.5), main="", col="red", xlab="", ylab="Density")
  legend("topright", legend=c("Matched", "Unmatched"), lty=c(1,2), col=c("red", "black"))
  
  # histograms of match and unmatch distributions
  hist(lasttermscale1[which(test %in% intersect(onetest, test[predmatch]))], nclass=90)
  hist(lasttermscale1[which(test %in% intersect(onetest, unmatch1))], nclass=90)
  
  # difference in mean of match and unmatch distributions
  mean(lasttermscale1[which(test %in% intersect(onetest, unmatch1))])
  mean(lasttermscale1[which(test %in% intersect(onetest, test[predmatch]))])
  
  # difference in quantiles of match and unmatch distributions
  quantile(lasttermscale1[which(test %in% intersect(onetest, unmatch1))], probs=seq(.1,.9,.1))
  quantile(lasttermscale1[which(test %in% intersect(onetest, test[predmatch]))], probs=seq(.1,.9,.1))
  
  
  ####################### timing. ##########################
  ####################### timing. ##########################
  ####################### timing. ##########################
  ####################### timing. ##########################
  ####################### timing. ##########################
  ## The test set sample size is very small and that informative for looking at timing (a lot of 
  ## programs have only one student). I think we need to bring in the full sample of completers,
  ## and predict their program in each term. This should also help us see program differences more clearly.
  load("course-wt2-svd1000COMPLETERS.Rdata")
  if (dim(translenmat)[1] != length(completers)) {
    wtscheme <- 2 
    wmat3 <- wtfunc(wmat2, wtscheme)
    q.default <- matrix(data=0, nrow=nrow(wmat3), ncol=1, dimnames=list(rownames(wmat3), "id"))
    translenmat <- matrix(data=0, nrow=length(completers), ncol=length(faclevels), dimnames=list(completers, paste(1:length(faclevels))))
    coursedistlist <- vector(mode="list", length=length(completers))
    names(coursedistlist) <- completers
    toplist <- vector(mode="list", length=length(completers))
    names(toplist) <- completers
    pb <- txtProgressBar(min = 0, max = length(completers)*length(faclevels), style = 3)
    timestamp()
    
    ## loop through the enrollment terms.
    for (ii in 1:length(faclevels)) {
      tmpcour <- subset(allcourse, subset=!duplist3 & allcourse$term.fac %in% faclevels[1:ii] & allcourse$id %in% completers & allcourse$pass, select=c("id", "course", "hegis_name"))
      loopdex <-na.omit(unique(tmpcour$id))
      for (jj in 1:length(loopdex)) {
        vicki <- loopdex[jj]
        vicki.q <- q.default
        partialcourse <- as.character(na.omit(unique(tmpcour$course[tmpcour$id==vicki])))
        if (length(partialcourse)==0) next  	
        translenmat[vicki,ii] <- length(partialcourse)
        if (ii>1) {	
          if (translenmat[vicki,ii]==translenmat[vicki,ii-1]) next
        }
        parttran <-	table(partialcourse)
        vicki.q[partialcourse,1] <- parttran
        vicki.q <-  wtfunc(vicki.q, wtscheme)
        vicki.q <- t(Uk) %*% vicki.q
        simmatsvd <- sapply(1:dim(award.medioid)[2], function(x) cosine(as.vector(vicki.q), award.medioid[,x]))
        names(simmatsvd) <- colnames(award.medioid)
        if (all(is.nan(simmatsvd))) next
        toplist[[vicki]][[ii]] <- sort(simmatsvd, decreasing=T)      
        coursedistlist[[vicki]][[ii]] <- sort(table(tmpcour$hegis_name[tmpcour$id==vicki]), decreasing=T)
      }
      setTxtProgressBar(pb, ii*length(completers))
    }
    timestamp()
    close(pb)
    save(translenmat, coursedistlist, toplist, ndim, award.medioid, file=paste("course-wt", wtscheme, "-svd", ndim, "COMPLETERS.Rdata", sep=''))
  }
  
  ## what is predicted field in last term
  lasttermmed <- sapply(completers, function(chuck) { 
    unlist(toplist[[chuck]][[length(toplist[[chuck]])]][1])
  })  
  
  testfields <- sapply(completers, function(x) allgrad[allgrad$id %in% x, "progofstudy"])
  predmatch <- sapply(1:length(completers), function(x) substring(names(lasttermmed)[x], first=nchar(test[1])+2) %in% as.character(testfields[[x]]))
  
  ## Use predicted program as "Actual" program in determining timing to concentration.
  predfield <- substring(names(lasttermmed), first=nchar(test[1])+2)
  names(predfield) <- completers
  
  ## used scaled cosine scores as a magnitude
  scalelist <- lapply(toplist, function(y) lapply(y, function(x) if(!is.null(x)) scale(x)))
  
  scaledpred <- lapply(completers, function(chuck) {
    prog <- predfield[chuck]
    if (prog=="") return(NA)
    unlist(lapply(scalelist[[chuck]], function(x) x[prog,]))})
  names(scaledpred) <- predfield
  
  ## this loop says for each predicted field, what are the standardized cosine dist in each term for each student
  ## then plot the trend for each student within a field to get an idea of program trajectory
  setwd('~/projects/virginia/virginia_PC/graphics/')
  
  for (ii in 1:lunique(predfield)) {
    png(filename=paste("scalecos", ii, substring(unique(predfield)[ii], first=5), ".png", sep=''))
    plot(1, xlim=c(0,16), ylim=c(0,16), type="n", main=unique(predfield)[ii], xlab="Term", ylab="Scaled cosine distance")
    l_ply(scaledpred[which(predfield==unique(predfield)[ii])], function(x) lines(x))
    dev.off()
  }
  
  ## predmax is the maximum cosine distance in each term, regardless of actual or predicted field
  predmax <- lapply(completers, function(chuck) 
    unlist(lapply(toplist[[chuck]], function(x) if (!is.null(x)) max(x))))
  names(predmax) <- completers
  
  for (ii in 1:lunique(predfield)) {
    png(filename=paste("cos", ii, substring(unique(predfield)[ii], first=5), ".png", sep=''))
    plot(1, xlim=c(0,16), ylim=c(0,1), type="n", main=unique(predfield)[ii], xlab="Term", ylab="Cosine distance")
    l_ply(predmax[which(predfield==unique(predfield)[ii])], function(x) lines(x))
    dev.off()
  }
  
  ## multiply transcript length by max cos dist each term to get some weighted confidence
  predconf2 <- lapply(completers, function(chuck)
    predmax[[chuck]]*translenmat[chuck,which(unlist(lapply(toplist[[chuck]], function(x) !is.null(x))))])
  
  ## use magnitude of biggest jump in weighted cosine distance as indicator of concentrating
  predconfdiffmag <- unlist(lapply(predconf2, function(x) max(diff(x))))
  names(predconfdiffmag) <- completers
  ## location of biggest jump
  predconfdiff2 <- as.numeric(unlist(lapply(predconf2, function(x) names(which.max(diff(x))))))
  names(predconfdiff2) <- completers
  
  ## magnitude of biggest jump
  predconfdiffmagunw <- unlist(lapply(predmax, function(x) max(diff(x))))
  names(predconfdiffmagunw) <- completers
  ## location of biggest jump
  predconfdiff2unw <- as.numeric(unlist(lapply(predmax, function(x) which.max(diff(x)))))+1
  names(predconfdiff2unw) <- completers
  
  pdensity(predconfdiffmag)
  pdensity(predconfdiffmagunw)
  
  chuck.terms <- sapply(completers, function(chuck) which(unlist(lapply(toplist[[chuck]], function(x) !is.null(x)))))
  ## term of concentration relative to start term
  concterm <- sapply(completers, function(chuck) which(which(unlist(lapply(toplist[[chuck]], function(x) !is.null(x)))) %in% predconfdiff2[chuck]))
  cbind(sort(tapply(concterm, predfield, mean)))
  conctermunw <- predconfdiff2unw
  names(conctermunw) <- completers
  cbind(sort(tapply(conctermunw, predfield, mean)))
  
  
  ## what are the most common courses that push students into concentration land? use weighted and unw concterm
  for (ii in lunique(predfield)) {
    chucks <- names(predfield[predfield==unique(predfield)[ii]])
    c1 <- vector()
    c2 <- vector()
    for (chuck in chucks) {
      c1 <- c(c1, names(table(allcourse$course[allcourse$id==chuck & allcourse$pass & allcourse$term.fac %in% faclevels[chuck.terms[[chuck]][1:conctermunw[chuck]]]])))
      c2 <- c(c2, names(table(allcourse$course[allcourse$id==chuck & allcourse$pass & allcourse$term.fac %in% faclevels[chuck.terms[[chuck]][1:conctermunw[chuck]-1]]])))
    }
    coursecomp <- merge(data.frame(table(c2)), data.frame(table(c1)), by.x="c2", by.y="c1", all=T)
    coursecomp$Freq.x[is.na(coursecomp$Freq.x)] <- 0
    coursecomp$diff <- coursecomp$Freq.y-coursecomp$Freq.x
    coursecomp <- coursecomp[order(coursecomp$diff, decreasing=T),]
    rbind(head(coursecomp), tail(coursecomp))
  }
  
  ## Predicted field in each term
  predmaxfield <- lapply(completers, function(chuck) 
    unlist(lapply(toplist[[chuck]], function(x) names(x)[1])))
  names(predmaxfield) <- completers
  
  ## what is the predicted field in the term of the max jump of the weighted cosine dist?
  prefieldatjump <- sapply(completers, function(chuck) predmaxfield[[chuck]][concterm[[chuck]]])
  prop.table(table(prefieldatjump==predfield))
  
  ## what is the predicted field in the term of the max jump of the unweighted cosine dist?
  prefieldatjumpunw <- sapply(completers, function(chuck) predmaxfield[[chuck]][conctermunw[[chuck]]])
  prop.table(table(prefieldatjumpunw==predfield))
  
  nterms <- sapply(completers, function(chuck) length(which(unlist(lapply(toplist[[chuck]], function(x) !is.null(x))))))
  
  boxplot(nterms~conctermunw, main="Last enrolled term vs. concentration term", xlab="Concentration term", ylab="Last term")
  boxplot(nterms~concterm, main="Last enrolled term vs. concentration term", xlab="Concentration term", ylab="Last term")
  
  ncourse <- unlist(lapply(completers, function(chuck) sum(coursedistlist[[chuck]][[predconfdiff2[chuck]]])))
  names(ncourse) <- completers
  boxplot(ncourse~concterm, main="Number of courses vs. concentration term", xlab="Concentration term", ylab="Number of courses in concentration term")
  cbind(sort(tapply(ncourse, predfield, mean)))
  
  ncourseunw <- unlist(lapply(completers, function(chuck) sum(coursedistlist[[chuck]][[chuck.terms[[chuck]][conctermunw[chuck]]]])))
  names(ncourseunw) <- completers
  boxplot(ncourseunw~conctermunw, main="Number of courses vs. concentration term", xlab="Concentration term", ylab="Number of courses in concentration term")
  cbind(sort(tapply(ncourseunw, predfield, mean)))
  
  
  ## predicted field each term
  whenpred <- sapply(completers, function(chuck) match(predfield[[chuck]], predmaxfield[[chuck]]))
  predfact <- factor(predfield, levels=names(sort(tapply(whenpred, predfield, median))))
  boxplot(whenpred~predfact)
  
  ncoursewhen <- unlist(lapply(completers, function(chuck) {
    chuck.terms <- which(unlist(lapply(toplist[[chuck]], function(x) !is.null(x))))
    sum(coursedistlist[[chuck]][[chuck.terms[whenpred[chuck]]]])}))
  names(ncoursewhen) <- completers
  
  boxplot(nterms~whenpred, main="Last enrolled term vs. concentration term", xlab="Concentration term", ylab="Last term")
  boxplot(ncoursewhen~whenpred, main="Number of courses vs. concentration term", xlab="Concentration term", ylab="Number of courses in concentration term")
  cbind(sort(tapply(whenpred, predfield, mean)))
  
  
  iis <- c(1, 2, 8, 10)
  for (ii in iis) {
    chucks <- names(predfield[predfield==unique(predfield)[ii]])
    c1 <- vector()
    c2 <- vector()
    for (chuck in chucks) {
      c1 <- c(c1, names(table(allcourse$course[allcourse$id==chuck & allcourse$pass & allcourse$term.fac %in% faclevels[chuck.terms[[chuck]][1:whenpred[chuck]]]])))
      c2 <- c(c2, names(table(allcourse$course[allcourse$id==chuck & allcourse$pass & allcourse$term.fac %in% faclevels[chuck.terms[[chuck]][1:whenpred[chuck]-1]]])))
    }
    coursecomp <- merge(data.frame(table(c2)), data.frame(table(c1)), by.x="c2", by.y="c1", all=T)
    coursecomp$Freq.x[is.na(coursecomp$Freq.x)] <- 0
    coursecomp$diff <- coursecomp$Freq.y-coursecomp$Freq.x
    coursecomp <- coursecomp[order(coursecomp$diff, decreasing=T),]
    print(rbind(head(coursecomp), tail(coursecomp)))
  }
  
  
  ################ MATT requests top 3 predicted programs for all students ##########
  # choose weight, build weighted matrix
  wtscheme <- 2 
  wmat3 <- wtfunc(wmat2, wtscheme)
  sing5 <- irlba(wmat3[,colnames(wmat3) %in% completers], nu=1000, nv=1000)
  wmat4 <- wtfunc(wmat1, wtscheme)
  
  ndim <- 1000
  Uk <- sing5$u[,1:ndim]
  Vk <- sing5$v[,1:ndim]
  Sk <- sing5$d[1:ndim]
  SkVk <- diag(Sk) %*% t(Vk)
  rownames(Uk) <- rownames(wmat3)
  colnames(SkVk) <- colnames(wmat3[,colnames(wmat3) %in% completers])
  awardcur <- allgrad$progofstudy[allgrad$id %in% completers]
  
  award.medioid <- sapply(unique(as.character(awardcur)), function(x)
  { 
    # this gets us the medioids in ndim-space for each award_curriculum in train
    ugh <- try(rowMeans(SkVk[,c(allgrad$id[allgrad$progofstudy==x & allgrad$id %in% completers])]), silent=TRUE)
    if (class(ugh)=="try-error") ugh <- SkVk[,c(allgrad$id[allgrad$progofstudy==x & allgrad$id %in% completers])]
    return(ugh)
  })
  
  
  ## only use last enrollment term - this is a most recent transcript info analysis
  tmpcour <- subset(allcourse, subset=!duplist3 & allcourse$pass, select=c("id", "course", "program_major"))
  loopdex <-na.omit(unique(tmpcour$id))
  
  ## create default query, all entries zero; transcript length matrix (num courses); top courses by program major; distance to award_curriculum medioids 
  q.default <- matrix(data=0, nrow=nrow(wmat4), ncol=1, dimnames=list(rownames(wmat4), "id"))
  translenmat <- matrix(data=0, nrow=length(loopdex), ncol=1, dimnames=list(loopdex, 1))
  coursedistlist <- vector(mode="list", length=length(loopdex))
  names(coursedistlist) <- loopdex
  toplist <- vector(mode="list", length=length(loopdex))
  names(toplist) <- loopdex
  
  ## jj is index of ids for test members). 
  for (jj in 1:length(loopdex)) {
    vicki <- loopdex[jj]
    vicki.q <- q.default
    partialcourse <- as.character(na.omit(unique(tmpcour$course[tmpcour$id==vicki])))
    
    ## ignore observation if student has no passed courses in the term
    if (length(partialcourse)==0) next  	
    
    ## ignore observation if no new information added in this term
    translenmat[vicki,1] <- length(partialcourse)
    #if (ii>1) {	
    #    if (translenmat[vicki,ii]==translenmat[vicki,ii-1]) next
    #  }
    
    ## generate base query
    parttran <-	table(partialcourse)
    vicki.q[partialcourse,1] <- parttran
    
    ## apply weighting scheme to query vector
    vicki.q <-  wtfunc(vicki.q, wtscheme)
    
    ## project partial transcript into course space
    #	and/or use distance metric to determine proximity to completers		
    vicki.q <- t(Uk) %*% vicki.q
    simmatsvd <- sapply(1:dim(award.medioid)[2], function(x) cosine(as.vector(vicki.q), award.medioid[,x]))
    names(simmatsvd) <- colnames(award.medioid)
    
    # skip if simmatsvd is all Nan
    if (all(is.nan(simmatsvd))) next
    
    toplist[[vicki]][[1]] <- sort(simmatsvd, decreasing=T)      
    coursedistlist[[vicki]][[1]] <- sort(table(tmpcour$program_major[tmpcour$id==vicki]), decreasing=T)
  }
  
  save(toplist, translenmat, file="PredictionsForEntireSample.Rdata")
  
  ## top 3 programs
  aa <- lapply(toplist, function(x) x[[1]][1:3])
  
  prognames <- lapply(aa, names)
  progfr <- do.call("rbind", prognames)
  colnames(progfr) <- c("Top1", "Top2", "Top3")
  
  progdist <- lapply(aa, function(x) as.numeric(x)) 
  distfr <- do.call("rbind", progdist)
  colnames(distfr) <- c("Top1", "Top2", "Top3")
  
  write.csv(progfr, file="VA_Top3programs.csv")
  write.csv(distfr, file="VA_Top3programs_dist.csv")
  write.csv(translenmat, file="VA_Top3programs_tranlen.csv")
  
  
