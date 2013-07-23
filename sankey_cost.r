## Code to produce JSON data set for budget Sankey diagram
## Simple one with just Fund -> Purpose
## read in data set and create summary table.
costs <- read.csv("./sankey/budg-11.csv")
yum <- tapply(costs$Closing.Balance,INDEX=list(costs$FUND.DESC,costs$PURP.1.DESC),FUN=sum, na.rm=T)
yum[is.na(yum)] <-0

library("RJSONIO")

## Create nodes list
aa <- c(rownames(yum), colnames(yum))
names(aa) <- rep("name", length(aa))
cc <- fromJSON("./sankey/energy.json")
nn <- NULL
for(ii in 1:length(aa)) nn <- c(nn,list(aa[ii]))
ee <- cc
ee[["nodes"]] <- nn

## Create links list
library(bipartite)
ff <- web2edges(yum1, return=TRUE)
rownames(ff) <- NULL
ff <- data.frame(ff)
ff$row <- ff$row-1
ff$col <- ff$col-1
colnames(ff) <- c("source", "target", "value")

ll <- NULL
for(ii in 1:nrow(ff)) ll <- c(ll,list(unlist(ff[ii,])))
ee[["links"]] <- ll

## save to file.
cat(toJSON(ee),file="./sankey/budg.json")

########################
## Let's kick it up a notch with FUND -> PURPOSE -> OBJECT.DESC 
# Need to build on codes above
yum2 <- tapply(costs$Closing.Balance,INDEX=list(costs$PURP.1.DESC,costs$OBJECT.DESC),FUN=sum, na.rm=T)

aa <- c(rownames(yum), colnames(yum), colnames(yum2))
names(aa) <- rep("name", length(aa))
cc <- fromJSON("./sankey/energy.json")
nn <- NULL
for(ii in 1:length(aa)) nn <- c(nn,list(aa[ii]))
gg <- cc
gg[["nodes"]] <- nn

nodevec <- as.character(unlist(gg$nodes))
hh <- web2edges(yum2, return=TRUE)
rownames(hh) <- NULL
hh <- data.frame(hh)
hh$row <- hh$row-1+length(rownames(yum))
hh$col <- hh$col-1+length(rownames(yum))
colnames(hh) <- c("source", "target", "value")

ww <- rbind(ff,hh)

ll <- NULL
for(ii in 1:nrow(ww)) ll <- c(ll,list(unlist(ww[ii,])))
gg[["links"]] <- ll

cat(toJSON(gg),file="./sankey/cpcc2.json")

########################
## Unwieldly.
## Maybe just for academic departments?
## FUND -> PURPOSE -> OBJECT.DESC -> UNIT.DESC
yum3 <- tapply(costs$Closing.Balance,INDEX=list(costs$OBJECT.DESC,costs$UNIT.DESC),FUN=sum, na.rm=T)

facsub <- as.character(unique(costs$UNIT.DESC[costs$OBJECT.DESC=="FT Faculty" | costs$OBJECT.DESC=="PT Faculty"]))
facsub <- as.character(unique(costs$UNIT.DESC[costs$UNIT.1==3]))
acadunits <- read.csv("/projects/CBD/crosta/data/CostsPerCredit0708-0910.csv")
## only include curriculum budget codes as last stop.
facsub <- as.character(unique(costs$UNIT.DESC[costs$UNIT %in% unique(acadunits$UnitId)]))
yum3 <- yum3[,facsub]

aa <- c(rownames(yum), colnames(yum), colnames(yum2), colnames(yum3))
names(aa) <- rep("name", length(aa))
cc <- fromJSON("/home/crosta/public_html/sankey/energy.json")
nn <- NULL
for(ii in 1:length(aa)) nn <- c(nn,list(aa[ii]))
qq <- cc
qq[["nodes"]] <- nn

nodevec <- as.character(unlist(gg$nodes))
zz <- web2edges(yum3, return=TRUE)
rownames(zz) <- NULL
zz <- data.frame(zz)
zz$row <- zz$row-1+length(rownames(yum))+length(rownames(yum2))
zz$col <- zz$col-1+length(rownames(yum))+length(rownames(yum2))
colnames(zz) <- c("source", "target", "value")

vv <- rbind(ff,hh,zz)

ll <- NULL
for(ii in 1:nrow(vv)) ll <- c(ll,list(unlist(vv[ii,])))
qq[["links"]] <- ll

cat(toJSON(qq),file="./sankey/cpcc3.json")

########### ISOLATE spending on FT and PT faculty  ###############
faculty <- subset(costs, subset=OBJECT.DESC=="FT Faculty" | OBJECT.DESC=="PT Faculty")
#faculty <- subset(costs, subset=(OBJECT.DESC=="FT Faculty" | OBJECT.DESC=="PT Faculty") & UNIT %in% unique(acadunits$UnitId))
facsum <- with(faculty, tapply(Closing.Balance, INDEX=list(as.character(UNIT.DESC), as.character(OBJECT.DESC)), FUN=sum))
facsum <- data.frame(facsum)
facsum$total <- rowSums(facsum[,c(1,2)], na.rm=T)
facsum$pctFT <- 100*facsum$FT.Faculty/facsum$total
facsum$pctPT <- 100*facsum$PT.Faculty/facsum$total
hist(facsum$pctFT)
plot(density(facsum$pctFT, na.rm=T))
#############################################################


