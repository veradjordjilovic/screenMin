#Data analysis for mediation project

library(statgenGWAS)
library(hdi)
data(dropsMarkers)  #41723 markers 
data(dropsPheno)    #We will use the outcome "yield" from this
unique(dropsPheno$Experiment)  #code names for the different trials (locations):
#W=watered, R=rainfed


###### obtaining p-values with statgenGWAS #######

set.seed(123)


## Add genotypes as row names of dropsMarkers and drop Ind column.
rownames(dropsMarkers) <- dropsMarkers[["Ind"]]
dropsMarkers <- dropsMarkers[colnames(dropsMarkers) != "Ind"]

data(dropsMap)

## Add genotypes as row names of dropsMap.
rownames(dropsMap) <- dropsMap[["SNP.names"]]
## Rename Chomosome and Position columns.
colnames(dropsMap)[match(c("Chromosome", "Position"), colnames(dropsMap))] <- c("chr", "pos")

## Create a gData object containing map and marker information.
gDataDrops <- createGData(geno = dropsMarkers, map = dropsMap)


## Rename Variety_ID to genotype.
colnames(dropsPheno)[colnames(dropsPheno) == "Variety_ID"] <- "genotype"
## Select relevant columns and convert data to a list.
dropsPhenoList <- split(x = dropsPheno[c("genotype", "grain.yield",
                                         "grain.number", "seed.size",
                                         "anthesis", "silking", "plant.height",
                                         "tassel.height", "ear.height")], 
                        f = dropsPheno[["Experiment"]])
## Add phenotypic data to gDataDrops.
gDataDrops <- createGData(gData = gDataDrops, pheno = dropsPhenoList)



## Remove duplicate SNPs from gDataDrops. 
gDataDropsDedup <- codeMarkers(gDataDrops, impute = FALSE, verbose = TRUE)   #36624 SNPs left



## Run single trait GWAS for traits 'grain.yield'  for trial Bol12R. (computation time: <10seconds)
GWASDrops <- runSingleTraitGwas(gData = gDataDropsDedup,
                                trials = c("Mur13W","Mur13R","Kar13W"), #select trial you want to analyze
                                traits = c("grain.yield"))


## Plot a QQ-plot of GWAS Drops.
plot(GWASDrops, plotType = "qq", trait = "grain.yield",trial="Mur13W")  #strong signal
plot(GWASDrops, plotType = "qq", trait = "grain.yield",trial="Mur13R")  #ok signal
plot(GWASDrops, plotType = "qq", trait = "grain.yield",trial="Kar13W")   #ok signal
if(FALSE){
  plot(GWASDrops, plotType = "qq", trait = "grain.yield",trial="Bol12R")
  plot(GWASDrops, plotType = "qq", trait = "grain.yield",trial="Kar12R")
  plot(GWASDrops, plotType = "qq", trait = "grain.yield",trial="Bol12W")
  plot(GWASDrops, plotType = "qq", trait = "grain.yield",trial="Kar12W")
  plot(GWASDrops, plotType = "qq", trait = "grain.yield",trial="Kar13R")   #weak
  plot(GWASDrops, plotType = "qq", trait = "grain.yield",trial="Gai13W")   #weak
}

Result = GWASDrops$GWAResult
summary(Result)

ResultMur13W = Result[[1]]
pvalsMur13W = ResultMur13W$pValue
length(pvalsMur13W)
hist(pvalsMur13W)
min(pvalsMur13W)

ResultMur13R = Result[[2]]
pvalsMur13R = ResultMur13R$pValue
min(pvalsMur13R)


ResultKar13W = Result[[3]]
pvalsKar13W = ResultKar13W$pValue
min(pvalsKar13W)

# ScreenMin analysis
p1 <-  pvalsKar13W
# p1 <- pvalsMur13R
# p1 <- pvalsMur13W
# p2 <-  pvalsMur13R
p2 <- pvalsMur13W

alpha <- 0.05
sum(p1 <= alpha/length(p1)) 
sum(p2 <= alpha/length(p2))
which(p1 <= alpha/length(p1))
which(p2<= alpha/length(p2))
minp <- pmin(p1, p2)
maxp <- pmax(p1, p2)

######################################
# Default ScreenMin Bonferroni method
######################################
S <- (which(pmin(p1,p2) < alpha/length(p1)))
length(S)
p.adjust(maxp[S], method="bonferroni")
SM <- (which(pmax(p1, p2) < alpha/ length(S)))
# discoveries
(SM <- intersect(S, SM))
length(SM)

#####################################
# Adaptive ScreenMin threshold
#####################################
sortmin <- sort(minp)
temp <- sortmin *(1:length(p1))
ind.temp <- which(temp>alpha)
if (length(ind.temp)>0) {k.prim <- ind.temp[1]
if (sortmin[k.prim]*(k.prim-1) <= alpha) {k <- k.prim}else{
  k <- k.prim-1}} else {
    k <- m}
adapt.threshold <- alpha/k
#discoveries
(adapt.S <- which(pmax(p1,p2) <= adapt.threshold))
length(adapt.S)
maxp[adapt.S]


#### PFER analysis
threshold <- FpControl(p1, p2, 1)
ind_pfer <- which(pmax(p1, p2) <= threshold)
pval_adjusted <- maxp[ind_pfer] * sum(minp <= threshold)

### HDMT analysis (not reported in the manuscript)
require(HDMT)
nullprop <- null_estimation(cbind(p1, p2), lambda = 0.9)
out1
library(help="HDMT")
fwercut0 <- fwer_est(nullprop$alpha10,nullprop$alpha01,nullprop$alpha00,nullprop$alpha1,
                     nullprop$alpha2,cbind(p1, p2),alpha=0.05,exact=0)
(hdmt<- which(pmax(p1,p2) <= fwercut0))
length(hdmt)

fwercut1 <- fwer_est(nullprop$alpha10,nullprop$alpha01,nullprop$alpha00,nullprop$alpha1,
                     nullprop$alpha2,cbind(p1, p2),alpha=0.05,exact=1)
(hdmt1<- which(pmax(p1,p2) <= fwercut1))
length(hdmt1)




###### A look at the two markers that Vera found:  #####

Markers = gDataDropsDedup$marker
dim(Markers)
#Vera found significant result for the 12652-th and 17925-th SNP, among the 36624 SNPs:
colnames(Markers)[12652] #"PUT-163a-148986271-678"
colnames(Markers)[17925] #"PZE-104137686"

SNPinfo =  read_csv("7b-InfoSNP_50K_41722.csv")
View(SNPinfo[1:20,])

##positions on B73 reference genome V2:
SNPinfo[SNPinfo$SNP.names=="PUT-163a-148986271-678",]  #chr3, position 155960553, allele1 A, allele2 C
SNPinfo[SNPinfo$SNP.names=="PZE-104137686",]     #chr4, position 224625164, allele1 T, allele2 C


Phen = gDataDropsDedup$pheno
which(unique(dropsPheno$Experiment)=="Kar13W")  #20
which(unique(dropsPheno$Experiment)=="Mur13W")  #24


dataKar13W = dropsPheno[dropsPheno$Experiment=="Kar13W",]
dataMur13W = dropsPheno[dropsPheno$Experiment=="Mur13W",]
#View(dataKar13W)
#sort phenotype data alphabetically by genotype name:
dataKar13W = dataKar13W[order(dataKar13W$genotype),  ]
dataMur13W = dataMur13W[order(dataMur13W$genotype),  ]

#the found markers are indeed strongy correlated with yield:
cor(dataKar13W$grain.yield, Markers[,12652])   #correlation of 0.3311227
cor(dataMur13W$grain.yield, Markers[,12652])   # corr of 0.3067632
cor(dataKar13W$grain.yield, Markers[,17925])  #corr of -0.3148384
cor(dataMur13W$grain.yield, Markers[,17925])  #corr of -0.3644031




save(pvalsMur13W,pvalsMur13R,pvalsKar13W, file = "pvalsVera.RData")

#save.image(file='workspmediation.RData')
#load("~/onderzoek/Fred van Eeuwijk/INVITA analyses/workspmediation.RData")


