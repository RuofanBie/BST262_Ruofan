### R code from vignette source 'BayesMendel.Rnw'

###################################################
### code chunk number 1: options
###################################################
options(width=65)


###################################################
### code chunk number 2: package
###################################################
library(BayesMendel)


###################################################
### code chunk number 3: package
###################################################

# Change future risk to be calculated in intervals of 2 y instead of the default of 5 y.
# Leave all other parameters as set.
myparams <- brcaparams(age.by=2)

# Run BRCAPRO with family history information for example family
out = brcapro(family=brca.fam)
slotNames(out)
out@probs
out@family



###################################################
### code chunk number 4: package
###################################################
# Turn off age imputation
out <- brcapro(family=brca.fam, imputeAges=FALSE, imputeRelatives=FALSE)

# Calculate risks with imputed ages
out = brcapro(family=brca.fam, imputeAges=TRUE, imputeRelatives=TRUE)
# When age imputation is done, the original
#family (with NA inputs re-coded to
#unaffected, age = 1) is returned by brcapro
out@family

# Can also impute ages, but not relatives.
out = brcapro(family=brca.fam, imputeAges=TRUE, imputeRelatives=FALSE)



###################################################
### code chunk number 5: package
###################################################
myparams <- brcaparams(penetrance = BRCApenet.Italian.2008)
out <- brcapro(family=brca.fam, params=myparams)


###################################################
### code chunk number 6: package
###################################################
out <- brcapro(family=brca.fam, race="Hispanic")


###################################################
### code chunk number 7: package
###################################################
# Add the testing results for BRCA1 and BRCA2
BRCA1 <- BRCA2 <- TestOrder <- rep(0,nrow(brca.fam))
germline.testing <- data.frame(BRCA1,BRCA2,TestOrder)
germline.testing[2,] <- c(2,0,1)
out <- brcapro(family=brca.fam, germline.testing=germline.testing)


###################################################
### code chunk number 8: package
###################################################
# Add the testing results for breast cancer markers
marker.testing <- data.frame(matrix(rep(0,nrow(brca.fam)*5),ncol=5))
colnames(marker.testing) <- c("ER","CK14","CK5.6","PR","HER2")
brca.fam[1,"AffectedBreast"] <- 1 
marker.testing[1,"ER"] <- 2 
out <- brcapro(family=brca.fam, germline.testing=germline.testing, marker.testing=marker.testing)


###################################################
### code chunk number 9: package
###################################################
# Add the information for oophorectomy
Oophorectomy <- c(1,rep(0,(nrow(brca.fam)-1)))
AgeOophorectomy <- c(30,rep(1,(nrow(brca.fam)-1)))
oophorectomy <- data.frame(Oophorectomy,AgeOophorectomy)
out <- brcapro(family=brca.fam, germline.testing=germline.testing, marker.testing=marker.testing, oophorectomy=oophorectomy)


###################################################
### code chunk number 10: package
###################################################
# Add the information for mastectomy
Mastectomy <- c(1,rep(0,(nrow(brca.fam)-1)))
AgeMastectomy <- c(57,rep(1,(nrow(brca.fam)-1)))
mastectomy <- data.frame(Mastectomy,AgeMastectomy)
out <- brcapro(family=brca.fam, mastectomy=mastectomy)


###################################################
### code chunk number 11: package
###################################################

# Change future risk to be calculated up to age 95 instead of the default 85.
# Leave all other parameters as set.
myparams <- MMRparams(age.to=95)

# Run MMRpro with family history information for example family
out = MMRpro(family=MMR.fam, params=myparams)



###################################################
### code chunk number 12: package
###################################################

## The counselee's father tested negative for MLH1 and MSH2.
## No testing for MSH6 was done.
MLH1 <- MSH2 <- MSH6 <- TestOrder <- rep(0, nrow(MMR.fam))
germline.testing = data.frame(MLH1, MSH2, MSH6, TestOrder)
germline.testing[3,] <- c(2,2,0,1)  

out <- MMRpro(family=MMR.fam, germline.testing = germline.testing)



###################################################
### code chunk number 13: package
###################################################

## Now let's say the counselee's sister has a colorectal tumor

MMR.fam[7, "AffectedColon"] <- 1

## The counselee's sister's tumor was found to be MSI high.
## Add in this MSI result.

MSI <- location <- rep(0, nrow(MMR.fam))
marker.testing <- data.frame(MSI, location)
marker.testing[7, "MSI"] <- 1

out <- MMRpro(family = MMR.fam, marker.testing = marker.testing)



###################################################
### code chunk number 14: package
###################################################

# Change the output for future risk to be calculated
# in age intervals of 1 year up to
# age 65 instead of the default 5 years.
# Leave all other parameters as set.
myparams <- pancparams(age.by=1, age.to=65)

# Run PancPRO with family history information for example family
pancpro(family=panc.fam, params=myparams)



###################################################
### code chunk number 15: package
###################################################

# Change likelihood ratio for single melanomas
# among noncarriers from default 0.702 to 0.80
myparams <- melaparams(spm.lr.noncarrier=0.80)

# Run PancPRO with family history information for example family
melapro(family=mela.fam, params=myparams)



###################################################
### code chunk number 16: package
###################################################

# The counselee's sister was tested for
# germline mutations in P16, and one was found.
# Proband was also tested, but no mutation was found.
P16 <- TestOrder <- rep(0, nrow(mela.fam))
germline.testing = data.frame(P16, TestOrder)
germline.testing[4,] <- c(1,1)
germline.testing[1,] <- c(2,2)

out <- melapro(family=mela.fam, germline.testing = germline.testing)



###################################################
### code chunk number 17: package
###################################################
pdf("brcafamplot.pdf")
brca.fam$Death <- rbinom(nrow(brca.fam), 1, 0.2)
myfamily <- new("BayesMendel", family=brca.fam, counselee.id=1)
plot.BayesMendel(myfamily, cex=0.2)
dev.off()


###################################################
### code chunk number 18: package
###################################################

pdf("mmrfamplot.pdf")
MMR.fam$Death <- rbinom(nrow(MMR.fam), 1, 0.2)
mmrpro.out <- MMRpro(family=MMR.fam, counselee.id=1)
plot(mmrpro.out, cex=0.2)
dev.off()


