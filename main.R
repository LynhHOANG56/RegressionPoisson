#################### Regression - d'occurrence de déclarations en ALD ############################
#################### 
x0 <- read.csv("data_lagged.csv", header=T) ; dim(x0); 
head(x0)


#################### Regression Lineaire
##### Model 1 - non spatial, 2- Autoregresse spatial, 3- Durbin spatial
##### Les meilleurs modeles
regressor1 = lm(formula = occurrence ~ 
                  TauxBac + TauxChomage + Salaire + TauxOuvrier + Club, data = x0)
summary(regressor1)            

regressor2 = lm(formula = occurrence ~ 
                  TauxBac + TauxChomage + TauxOuvrier + Club + yLag, data = x0)
summary(regressor2)  

regressor3 = lm(formula = occurrence ~ 
                  TauxBac + TauxOuvrier + Club + yLag
                  + x0$SalaireLag + x0$TauxOuvrierLag + x0$ClubLag, data = x0)
summary(regressor3)

#################### Regression Poisson
##### Model A - non spatial, B- Autoregresse spatial, C- Durbin spatial
exposure <- x0$occurrence/x0$Pop2017
regressorA = glm(formula = occurrence ~ 
                   TauxBac + TauxChomage + Salaire + TauxOuvrier + TauxClub
                 + offset(log(exposure)),
                 data = x0, family = poisson(link = log))
summary(regressorA)

regressorB = glm(formula = occurrence ~ 
                   TauxBac + TauxChomage + Salaire + TauxOuvrier + TauxClub + yLag
                 + offset(log(exposure)),
                 data = x0, family = poisson(link = log))
summary(regressorB)

regressorC = glm(formula = occurrence ~ 
                   TauxBac + TauxChomage + Salaire + TauxOuvrier + TauxClub + yLag
                 + x0$TauxBacLag + x0$TauxChomageLag + x0$SalaireLag + x0$TauxOuvrierLag + x0$TauxClubLag
                 + offset(log(exposure)),
                 data = x0, family = poisson(link = log))
summary(regressorC)

######################## residual
###### 1,2,3 
z1 <- residuals(regressor1); 
z2 <- residuals(regressor2); 
z3 <- residuals(regressor3); 

###### A,B,C 
zA <- residuals(regressorA); 
zB <- residuals(regressorB); 
zC <- residuals(regressorC); 

residus <- cbind(x0$Code, x0$indD, z1,z2,z3,zA,zB,zC)
colnames(residus) <- c("Dep", "indD", "Model1","Model2","Model3", "ModelA","ModelB","ModelC")
resStand <- cbind(x0$Code, x0$indD,scale(z1),scale(z2),scale(z3), scale(zA), scale(zB), scale(zC))
colnames(resStand) <- c("Dep", "indD", "Model1","Model2","Model3", "ModelA","ModelB","ModelC")

boxplot(residus[,-c(1,2)], main="Les résidus des modèles 1-2-3 et A-B-C")
boxplot(resStand[,-c(1,2)], main="Les résidus standardisés des modèles 1-2-3 et A-B-C")

###### Plot carte de Residus/Residus standardisés
###### La France 
library(cartography)
library(raster)
res <- residus
jj = 8 # 3,4, 5,6,7, 8
z <- res[,jj]
formes <- getData(name="GADM", country="FRA", level=2)
concordance <- res[res[,2], jj]; concordance
formes$z <- concordance
mycols <- carto.pal(pal1 = "blue.pal", n1 = 9, pal2 = "orange.pal", n2 = 9)
ech <-  c(-3.45, -3.4, -2.9, -2.36, -2, -1.5, -1, 0, 1,1.5,2.3,2.5,3.7,4.5,4.77,4.9)
spplot(formes, "z", col.regions=mycols,
       at = ech,#couleurs(30), 
       main=list(label="Résidus standardisés, Durbin-spatial modèle",cex=1.4))

###### L'Ile-de-France
France <- getData(name="GADM", country="FRA", level=2)
FrancePartie<- subset(France, NAME_1=="Île-de-France")
#plot(FrancePartie, main="Carte de la France")
nI <- FrancePartie$NAME_2; nI
iI <- c()
for (j in 1:length(nI)){
  iI[j] <- which(x0$Nom==FrancePartie$NAME_2[j])
}
idI <- x0$indD[iI]; idI; iI
concordanceI <- res[iI, jj]; concordanceI #
FrancePartie$z <- concordanceI
mycols <- carto.pal(pal1 = "blue.pal", n1 = 9, pal2 = "orange.pal", n2 = 9)
spplot(FrancePartie, "z", col.regions=mycols,
       at =ech ,#couleurs(30), 
       main=list(label="Ile-de-France",cex=1.8))


######################################################
## Test Moran I
install.packages("ape")
library(ape)
occ.distsA <- as.matrix(dist(cbind(x0$TauxBac, x0$TauxChomage, x0$Salaire, x0$TauxOuvrier, x0$TauxClub)))
z <- zA
occ.distsB <- as.matrix(dist(cbind(x0$TauxBac, x0$TauxChomage, x0$Salaire, x0$TauxOuvrier, x0$TauxClub, x0$yLag)))
z <- zB
occ.distsC <- as.matrix(dist(cbind(x0$TauxBac, x0$TauxChomage, x0$Salaire, x0$TauxOuvrier,
                                   x0$ylab, x0$TauxClub,
                                   x0$TauxBacLag, x0$TauxChomageLag, x0$SalaireLag,x0$TauxOuvrierLag,
                                   x0$TauxClubLag)))
z <- zC

MoranTest <- function(z, occ.dists){
  occ.dists; dim(occ.dists)
  occ.dists.inv <- 1/occ.dists
  diag(occ.dists.inv) <- 0
  occ.dists.inv[1:5, 1:5]
  return(Moran.I(z, occ.dists.inv))
}
MoranTest(zA, occ.distsA)
MoranTest(zB, occ.distsB)
MoranTest(zC, occ.distsC)
