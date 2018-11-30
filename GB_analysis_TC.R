
################################################

# Tom's additions start here
#setwd("C:/Users/cristto/Desktop/Documents/Manuscripts and Data/Leno")
getwd()
setwd("C:/Users/davisao2/Desktop/open source GIS/ebird r code")
library(ape)# includes moran's i function moran.i
require(cooccur)# for cooccurrence matrices
library(fpc)#includes Clujaccard
require(pgirmess)
library(readr)
library(reshape2)
library(vegan)#


GBeBirdDEC2016 <- read_csv("GBeBirdDEC2016TerrestrialHabitat.csv")

str(GBeBirdDEC2016)


#create data frame with all the observers that have conducted surveys and the number of surveys at each location####
LocalityObs<-dcast(GBeBirdDEC2016,
                   `LOCALITY ID`+ 
                     LATITUDE+
                     LONGITUDE+
                     Habitat~
                     `OBSERVER ID`,
                   value.var="SAMPLING EVENT IDENTIFIER",  fun.aggregate = function(x) length(unique(x))
)

LocalityObsMinutes<-dcast(GBeBirdDEC2016,
                          `LOCALITY ID`+ 
                            LATITUDE+
                            LONGITUDE+
                            Habitat~
                            `OBSERVER ID`,
                          value.var="DURATION MINUTES",  fun.aggregate = sum
)


#calculate Observer diversity indices for each locality, Richness, Shannon and Simpson####

LocalityObs$ObserverRichness<-rowSums(LocalityObs[,5:284]>0)
LocalityObs$TotalSurveys<-rowSums(LocalityObs[,5:284])
LocalityObs$ObserverShannonDiversity<-apply(LocalityObs[,5:284], 1, diversity, index = "shannon")
LocalityObs$ObserverSimpsonDiversity<-apply(LocalityObs[,5:284], 1, diversity, index = "simpson")

#Create Data Frame of all species seen at each locality and the number of surveys they turn up in.####
LocalitySps<-dcast(GBeBirdDEC2016,
                   `LOCALITY ID`+ 
                     LATITUDE+
                     LONGITUDE+
                     Habitat~
                     `COMMON NAME`,
                   value.var="SAMPLING EVENT IDENTIFIER",  fun.aggregate = function(x) length(unique(x))
)

#Calculate Species Diversity indices for each locality, Species Richness, Shannon and Simpson Diversity####
LocalitySps$SpeciesRichness<-rowSums(LocalitySps[,5:282]>0)
LocalitySps$SpeciesShannonDiversity<-apply(LocalitySps[,5:282], 1, diversity, index = "shannon")
LocalitySps$SpeciesSimpsonDiversity<-apply(LocalitySps[,5:282], 1, diversity, index = "simpson")

#Create Data Frame of observer effort at each location####
LocalityObsMinutes[is.na(LocalityObsMinutes)]<-0
LocalityObsMinutes$TotalMinutes<-rowSums(LocalityObsMinutes[,5:284], na.rm=TRUE)
LocalityObsMinutes$ObserverMinutesShannonDiversity<-apply(LocalityObsMinutes[,5:284], 1, diversity, index = "shannon")
LocalityObsMinutes$ObserverMinutesSimpsonDiversity<-apply(LocalityObsMinutes[,5:284], 1, diversity, index = "simpson")

#Merge diversity data####
ObsDiversity<-as.data.frame(c(LocalityObs[,c(1:4, 285:287)]))
SpsDiversity<-as.data.frame(c(LocalitySps[,c(1:4, 283:285)]))
ObsMinutesDiversity<-as.data.frame(c(LocalityObsMinutes[,c(1:4, 285:287)]))
AllDiversity<-merge(ObsDiversity, SpsDiversity)
AllDiversity<-merge(AllDiversity, ObsMinutesDiversity)
#write.csv(AllDiversity, "AllDiversityGBDEC2016.csv")

##Habitat diversity measures####
WaterDiversity<-AllDiversity[AllDiversity$Habitat=="Water",]
PineDiversity<-AllDiversity[AllDiversity$Habitat=="Pine",]
WetlandDiversity<-AllDiversity[AllDiversity$Habitat=="Wetland",]
SandDiversity<-AllDiversity[AllDiversity$Habitat=="Sand",]
UrbanDiversity<-AllDiversity[AllDiversity$Habitat=="Urban",]
GrassDiversity<-AllDiversity[AllDiversity$Habitat=="Grass",]
HWTCDiversity<-AllDiversity[AllDiversity$Habitat=="HWTC",]

### Here I am creating a species by sample matrix for species accumulation curves and ordinations
## Add the species x sample matrix to the AllDiversity data
SpecSamp<-merge(AllDiversity,LocalitySps)

str(SpecSamp)
SpecSampMat<-as.matrix(SpecSamp[,15:292])

## Total Species Accumulation Curve by surveys, minutes, or localities

cum_surveys<-mean(SpecSamp$TotalSurvey)*seq(1,nrow(SpecSamp))
cum_minutes<-mean(SpecSamp$TotalMinutes)*seq(1,nrow(SpecSamp))
rc<-specaccum(SpecSampMat,method="exact")
cum_ind<-rc$individuals
cum_rich<-rc$richness
cum_samp<-rc$sites
rich_sd<-rc$sd
rich_low<-cum_rich-1.96*rich_sd
rich_hi<-cum_rich+1.96*rich_sd
plot(cum_rich~cum_samp,type="o",xlab="Number of Localities",ylab="Species Richness",
     xlim=c(0,520),ylim=c(0,370),cex.axis=1.2,cex.lab=1.2)
segments(cum_samp,rich_low,cum_samp,rich_hi,col="grey",lwd=2)
points(cum_samp,cum_rich,pch=21)
lines(cum_samp,cum_rich,col="black")
rich_est<-specpool(SpecSampMat)
points(500,rich_est$boot,pch=16,cex=1.2)
segments(500,rich_est$boot-(1.96*rich_est$boot.se),500,rich_est$boot+(1.96*rich_est$boot.se))
points(510,rich_est$jack1,pch=16,cex=1.2)
segments(510,rich_est$jack1-(1.96*rich_est$jack1.se),510,rich_est$jack1+(1.96*rich_est$jack1.se))
points(520,rich_est$chao,pch=16,cex=1.2)
segments(520,rich_est$chao-(1.96*rich_est$chao.se),520,rich_est$chao+(1.96*rich_est$chao.se))
#abline(a=rich_est$chao,b=0,lty=2)
#abline(a=rich_est$jack1,b=0,lty=3)
#abline(a=rich_est$boot,b=0,lty=4)
text(473,350,"Chao")
text(460,330,"Jackknife")
text(457,300,"Bootstrap")

### Species Rarefaction Curves by Habitat
WaterSpecSamp<-SpecSamp[SpecSamp$Habitat=="Water",]
WaterSpecSampMat<-as.matrix(WaterSpecSamp[,15:292])
PineSpecSamp<-SpecSamp[SpecSamp$Habitat=="Pine",]
PineSpecSampMat<-as.matrix(PineSpecSamp[,15:292])
WetlandSpecSamp<-SpecSamp[SpecSamp$Habitat=="Wetland",]
WetlandSpecSampMat<-as.matrix(WetlandSpecSamp[,15:292])
SandSpecSamp<-SpecSamp[SpecSamp$Habitat=="Sand",]
SandSpecSampMat<-as.matrix(SandSpecSamp[,15:292])
UrbanSpecSamp<-SpecSamp[SpecSamp$Habitat=="Urban",]
UrbanSpecSampMat<-as.matrix(UrbanSpecSamp[,15:292])
GrassSpecSamp<-SpecSamp[SpecSamp$Habitat=="Grass",]
GrassSpecSampMat<-as.matrix(GrassSpecSamp[,15:292])
HWTCSpecSamp<-SpecSamp[SpecSamp$Habitat=="HWTC",]
HWTCSpecSampMat<-as.matrix(HWTCSpecSamp[,15:292])

# water
wat_rc<-specaccum(WaterSpecSampMat,method="exact")
cum_ind<-wat_rc$individuals
cum_rich<-wat_rc$richness
cum_samp<-wat_rc$sites
rich_sd<-wat_rc$sd
rich_low<-cum_rich-1.96*rich_sd
rich_hi<-cum_rich+1.96*rich_sd
plot(cum_rich~cum_samp,type="l",xlab="Number of Localities",ylab="Bird Species Richness",
     xlim=c(0,120),ylim=c(0,250),col="blue",lwd=2,pch=16,cex.axis=1.2,cex.lab=1.2)
text(40,110,"Water",cex=1.2)
# pine
pine_rc<-specaccum(PineSpecSampMat,method="exact")
cum_ind<-pine_rc$individuals
cum_rich<-pine_rc$richness
cum_samp<-pine_rc$sites
rich_sd<-pine_rc$sd
rich_low<-cum_rich-1.96*rich_sd
rich_hi<-cum_rich+1.96*rich_sd
lines(cum_rich~cum_samp,col="darkgreen",lwd=2)
#points(cum_rich~cum_samp,col="darkgreen",pch=16)
text(85,150,"Pine",cex=1.2)
# wetlands
wet_rc<-specaccum(WetlandSpecSampMat,method="exact")
cum_ind<-wet_rc$individuals
cum_rich<-wet_rc$richness
cum_samp<-wet_rc$sites
rich_sd<-wet_rc$sd
rich_low<-cum_rich-1.96*rich_sd
rich_hi<-cum_rich+1.96*rich_sd
lines(cum_rich~cum_samp,col="purple",lwd=2)
#points(cum_rich~cum_samp,col="purple",pch=16)
text(53,185,"Wetlands",cex=1.2)
# sand
sand_rc<-specaccum(SandSpecSampMat,method="exact")
cum_ind<-sand_rc$individuals
cum_rich<-sand_rc$richness
cum_samp<-sand_rc$sites
rich_sd<-sand_rc$sd
rich_low<-cum_rich-1.96*rich_sd
rich_hi<-cum_rich+1.96*rich_sd
lines(cum_rich~cum_samp,col="sandybrown",lwd=2)
#points(cum_rich~cum_samp,col="sandybrown",pch=16)
text(100,195,"Sand",cex=1.2)
# urban
urban_rc<-specaccum(UrbanSpecSampMat,method="exact")
cum_ind<-urban_rc$individuals
cum_rich<-urban_rc$richness
cum_samp<-urban_rc$sites
rich_sd<-urban_rc$sd
rich_low<-cum_rich-1.96*rich_sd
rich_hi<-cum_rich+1.96*rich_sd
lines(cum_rich~cum_samp,col="black",lwd=2)
#points(cum_rich~cum_samp,col="black",pch=16)
text(100,165,"Urban",cex=1.2)
# grass
grass_rc<-specaccum(GrassSpecSampMat,method="exact")
cum_ind<-grass_rc$individuals
cum_rich<-grass_rc$richness
cum_samp<-grass_rc$sites
rich_sd<-grass_rc$sd
rich_low<-cum_rich-1.96*rich_sd
rich_hi<-cum_rich+1.96*rich_sd
lines(cum_rich~cum_samp,col="limegreen",lwd=2)
#points(cum_rich~cum_samp,col="limegreen",pch=16)
text(100,230,"Grass",cex=1.2)
# HWTC
hwtc_rc<-specaccum(HWTCSpecSampMat,method="exact")
cum_ind<-hwtc_rc$individuals
cum_rich<-hwtc_rc$richness
cum_samp<-hwtc_rc$sites
rich_sd<-hwtc_rc$sd
rich_low<-cum_rich-1.96*rich_sd
rich_hi<-cum_rich+1.96*rich_sd
lines(cum_rich~cum_samp,col="magenta",lwd=2)
#points(cum_rich~cum_samp,col="magenta",pch=16)
text(23,125,"HWTC",cex=1.2)




### General linear model of species richness by habitat using a Poisson error distribution
m1<-glm(SpeciesRichness ~ Habitat + log(TotalSurveys),family="poisson",data=SpecSamp)
summary(m1)

# note this uses a log scale for the y axis
boxplot(SpeciesRichness ~ Habitat, data=SpecSamp,log="y")

plot(SpecSamp$TotalSurveys,SpecSamp$SpeciesRichness,log="xy",xlab="Number of Transect Surveys",
     ylab="Oberved Number of Bird Species")
curve(exp(2.394+0.479*log(x)),add=TRUE)

# Do some ordinations by habitat type

# transform abundance date using quarter power
# use Bray-Curtis dissimilarities
SpecSampMatTrans<-SpecSampMat^0.25

# NMDS does not converge
# ord.m1<-metaMDS(SpecSampMatTrans,dist="bray",k=3,autotransform=FALSE)

# Use Distance-based redundancy analysis for ordination with habitat as a constraining variable
ord.m2<-dbrda(SpecSampMatTrans~SpecSamp$Habitat+log(SpecSamp$TotalSurveys),dist="bray")
summary(ord.m2)
anova(ord.m2,by="term") # uses permuation test to get p-values for differences in species composition by habitat and sampling efforts

# get habitat labels and plot ordination with habitat centroids and ellipses
hab_type<-SpecSamp$Habitat
site_scores<-scores(ord.m2,display="wa")
row.names(site_scores)<-hab_type
hab_names<-c("Grass","HWTC","Pine","Sand","Urban","Water","Wetland")
hab_col<-c("limegreen","magenta","darkgreen","sandybrown","black","blue","purple")
hab_centroids<-scores(ord.m2,display="cn")
plot(site_scores[,1],site_scores[,2],xlim=c(-1.5,1.5),ylim=c(-1.5,1.5),xlab="Axis 1",ylab="Axis 2",pch=3,col="grey",
     cex.axis=1.2,cex=0.7,cex.lab=1.2)
points(hab_centroids[,1],hab_centroids[,2],pch=15,col=hab_col,cex=1.3)
cen_label_x<-hab_centroids[,1]
cen_label_y<-hab_centroids[,2]+0.15
cen_label_x[7]<-hab_centroids[7,1]+0.1
cen_label_y[7]<-hab_centroids[7,2]-0.1
ordiellipse(ord.m2,hab_type,col=hab_col,lwd=1,kind="se",conf=0.99)
text(cen_label_x,cen_label_y,hab_names,col=hab_col,cex=1.2)
points(0,0,pch=19,col="black",cex=1.2)
arrows(0,0,-0.22,0.96,lwd=2)
text(-0.25,1.05,"Number of Surveys",cex=1.2)


#### Tom's edits stop here - I'll pick up spatial analysis next

#Moran's I####
#Distance weights####
#All Data###
All.dists <- as.matrix(dist(cbind(AllDiversity$LONGITUDE, AllDiversity$LATITUDE)))
All.dists.inv <- 1/All.dists
diag(All.dists.inv) <- 0
All.dists.inv[is.infinite(All.dists.inv)] <- 0

All.dists.inv[1:5, 1:5]
#Water
Water.dists<-as.matrix(dist(cbind(AllDiversity$LONGITUDE[AllDiversity$Habitat=="Water"], AllDiversity$LATITUDE[AllDiversity$Habitat=="Water"])))
Water.dists.inv <- 1/Water.dists
diag(Water.dists.inv) <- 0
Water.dists.inv[is.infinite(Water.dists.inv)] <- 0
Water.dists.inv[1:5, 1:5]
#Pine
Pine.dists<-as.matrix(dist(cbind(AllDiversity$LONGITUDE[AllDiversity$Habitat=="Pine"], AllDiversity$LATITUDE[AllDiversity$Habitat=="Pine"])))
Pine.dists.inv <- 1/Pine.dists
diag(Pine.dists.inv) <- 0
Pine.dists.inv[is.infinite(Pine.dists.inv)] <- 0
Pine.dists.inv[1:5, 1:5]
#Wetland
Wetland.dists<-as.matrix(dist(cbind(AllDiversity$LONGITUDE[AllDiversity$Habitat=="Wetland"], AllDiversity$LATITUDE[AllDiversity$Habitat=="Wetland"])))
Wetland.dists.inv <- 1/Wetland.dists
diag(Wetland.dists.inv) <- 0
Wetland.dists.inv[is.infinite(Wetland.dists.inv)] <- 0
Wetland.dists.inv[1:5, 1:5]
#Sand
Sand.dists<-as.matrix(dist(cbind(AllDiversity$LONGITUDE[AllDiversity$Habitat=="Sand"], AllDiversity$LATITUDE[AllDiversity$Habitat=="Sand"])))
Sand.dists.inv <- 1/Sand.dists
diag(Sand.dists.inv) <- 0
Sand.dists.inv[is.infinite(Sand.dists.inv)] <- 0
Sand.dists.inv[1:5, 1:5]
#Urban
Urban.dists<-as.matrix(dist(cbind(AllDiversity$LONGITUDE[AllDiversity$Habitat=="Urban"], AllDiversity$LATITUDE[AllDiversity$Habitat=="Urban"])))
Urban.dists.inv <- 1/Urban.dists
diag(Urban.dists.inv) <- 0
Urban.dists.inv[is.infinite(Urban.dists.inv)] <- 0
Urban.dists.inv[1:5, 1:5]
#Grass
Grass.dists<-as.matrix(dist(cbind(AllDiversity$LONGITUDE[AllDiversity$Habitat=="Grass"], AllDiversity$LATITUDE[AllDiversity$Habitat=="Grass"])))
Grass.dists.inv <- 1/Grass.dists
diag(Grass.dists.inv) <- 0
Grass.dists.inv[is.infinite(Grass.dists.inv)] <- 0
Grass.dists.inv[1:5, 1:5]
#HWTC
HWTC.dists<-as.matrix(dist(cbind(AllDiversity$LONGITUDE[AllDiversity$Habitat=="HWTC"], AllDiversity$LATITUDE[AllDiversity$Habitat=="HWTC"])))
HWTC.dists.inv <- 1/HWTC.dists
diag(HWTC.dists.inv) <- 0
HWTC.dists.inv[is.infinite(HWTC.dists.inv)] <- 0
HWTC.dists.inv[1:5, 1:5]



#Moran's I of Geo Distances and Observer and Species Diversity
#All Data####
set.seed(1981);Moran.I(AllDiversity$ObserverRichness, All.dists.inv)
set.seed(1981);Moran.I(AllDiversity$ObserverShannonDiversity, All.dists.inv)
set.seed(1981);Moran.I(AllDiversity$ObserverSimpsonDiversity, All.dists.inv)
set.seed(1981);Moran.I(AllDiversity$SpeciesRichness, All.dists.inv)
set.seed(1981);Moran.I(AllDiversity$SpeciesShannonDiversity, All.dists.inv)
set.seed(1981);Moran.I(AllDiversity$SpeciesSimpsonDiversity, All.dists.inv)
set.seed(1981);Moran.I(AllDiversity$TotalMinutes, All.dists.inv)
set.seed(1981);Moran.I(AllDiversity$TotalSurveys, All.dists.inv)
set.seed(1981);Moran.I(AllDiversity$ObserverMinutesShannonDiversity, All.dists.inv)
set.seed(1981);Moran.I(AllDiversity$ObserverMinutesSimpsonDiversity, All.dists.inv)

summary(AllDiversity$ObserverRichness)

#By Habitat####
#Water

set.seed(1981);Moran.I(WaterDiversity$ObserverRichness, Water.dists.inv)
set.seed(1981);Moran.I(WaterDiversity$ObserverShannonDiversity, Water.dists.inv)
set.seed(1981);Moran.I(WaterDiversity$ObserverSimpsonDiversity, Water.dists.inv)
set.seed(1981);Moran.I(WaterDiversity$SpeciesRichness, Water.dists.inv)
set.seed(1981);Moran.I(WaterDiversity$SpeciesShannonDiversity, Water.dists.inv)
set.seed(1981);Moran.I(WaterDiversity$SpeciesSimpsonDiversity, Water.dists.inv)
set.seed(1981);Moran.I(WaterDiversity$TotalMinutes, Water.dists.inv)
set.seed(1981);Moran.I(WaterDiversity$TotalSurveys, Water.dists.inv)
set.seed(1981);Moran.I(WaterDiversity$ObserverMinutesShannonDiversity, Water.dists.inv)
set.seed(1981);Moran.I(WaterDiversity$ObserverMinutesSimpsonDiversity, Water.dists.inv)
summary(WaterDiversity$ObserverRichness)

#Pine
set.seed(1981);Moran.I(PineDiversity$ObserverRichness, Pine.dists.inv)
set.seed(1981);Moran.I(PineDiversity$ObserverShannonDiversity, Pine.dists.inv)
set.seed(1981);Moran.I(PineDiversity$ObserverSimpsonDiversity, Pine.dists.inv)
set.seed(1981);Moran.I(PineDiversity$SpeciesRichness, Pine.dists.inv)
set.seed(1981);Moran.I(PineDiversity$SpeciesShannonDiversity, Pine.dists.inv)
set.seed(1981);Moran.I(PineDiversity$SpeciesSimpsonDiversity, Pine.dists.inv)
set.seed(1981);Moran.I(PineDiversity$TotalMinutes, Pine.dists.inv)
set.seed(1981);Moran.I(PineDiversity$TotalSurveys, Pine.dists.inv)
set.seed(1981);Moran.I(PineDiversity$ObserverMinutesShannonDiversity, Pine.dists.inv)
set.seed(1981);Moran.I(PineDiversity$ObserverMinutesSimpsonDiversity, Pine.dists.inv)

summary(PineDiversity$ObserverRichness)

#Wetland
set.seed(1981);Moran.I(WetlandDiversity$ObserverRichness, Wetland.dists.inv)
set.seed(1981);Moran.I(WetlandDiversity$ObserverShannonDiversity, Wetland.dists.inv)
set.seed(1981);Moran.I(WetlandDiversity$ObserverSimpsonDiversity, Wetland.dists.inv)
set.seed(1981);Moran.I(WetlandDiversity$SpeciesRichness, Wetland.dists.inv)
set.seed(1981);Moran.I(WetlandDiversity$SpeciesShannonDiversity, Wetland.dists.inv)
set.seed(1981);Moran.I(WetlandDiversity$SpeciesSimpsonDiversity, Wetland.dists.inv)
set.seed(1981);Moran.I(WetlandDiversity$TotalMinutes, Wetland.dists.inv)
set.seed(1981);Moran.I(WetlandDiversity$TotalSurveys, Wetland.dists.inv)
set.seed(1981);Moran.I(WetlandDiversity$ObserverMinutesShannonDiversity, Wetland.dists.inv)
set.seed(1981);Moran.I(WetlandDiversity$ObserverMinutesSimpsonDiversity, Wetland.dists.inv)
summary(WetlandDiversity$ObserverRichness)
#Sand
set.seed(1981);Moran.I(SandDiversity$ObserverRichness, Sand.dists.inv)
set.seed(1981);Moran.I(SandDiversity$ObserverShannonDiversity, Sand.dists.inv)
set.seed(1981);Moran.I(SandDiversity$ObserverSimpsonDiversity, Sand.dists.inv)
set.seed(1981);Moran.I(SandDiversity$SpeciesRichness, Sand.dists.inv)
set.seed(1981);Moran.I(SandDiversity$SpeciesShannonDiversity, Sand.dists.inv)
set.seed(1981);Moran.I(SandDiversity$SpeciesSimpsonDiversity, Sand.dists.inv)
set.seed(1981);Moran.I(SandDiversity$TotalMinutes, Sand.dists.inv)
set.seed(1981);Moran.I(SandDiversity$TotalSurveys, Sand.dists.inv)
set.seed(1981);Moran.I(SandDiversity$ObserverMinutesShannonDiversity, Sand.dists.inv)
set.seed(1981);Moran.I(SandDiversity$ObserverMinutesSimpsonDiversity, Sand.dists.inv)
summary(SandDiversity$ObserverRichness)
#Urban
set.seed(1981);Moran.I(UrbanDiversity$ObserverRichness, Urban.dists.inv)
set.seed(1981);Moran.I(UrbanDiversity$ObserverShannonDiversity, Urban.dists.inv)
set.seed(1981);Moran.I(UrbanDiversity$ObserverSimpsonDiversity, Urban.dists.inv)
set.seed(1981);Moran.I(UrbanDiversity$SpeciesRichness, Urban.dists.inv)
set.seed(1981);Moran.I(UrbanDiversity$SpeciesShannonDiversity, Urban.dists.inv)
set.seed(1981);Moran.I(UrbanDiversity$SpeciesSimpsonDiversity, Urban.dists.inv)
set.seed(1981);Moran.I(UrbanDiversity$TotalMinutes, Urban.dists.inv)
set.seed(1981);Moran.I(UrbanDiversity$TotalSurveys, Urban.dists.inv)
set.seed(1981);Moran.I(UrbanDiversity$ObserverMinutesShannonDiversity, Urban.dists.inv)
set.seed(1981);Moran.I(UrbanDiversity$ObserverMinutesSimpsonDiversity, Urban.dists.inv)
summary(UrbanDiversity$ObserverRichness)
#Grass
set.seed(1981);Moran.I(GrassDiversity$ObserverRichness, Grass.dists.inv)
set.seed(1981);Moran.I(GrassDiversity$ObserverShannonDiversity, Grass.dists.inv)
set.seed(1981);Moran.I(GrassDiversity$ObserverSimpsonDiversity, Grass.dists.inv)
set.seed(1981);Moran.I(GrassDiversity$SpeciesRichness, Grass.dists.inv)
set.seed(1981);Moran.I(GrassDiversity$SpeciesShannonDiversity, Grass.dists.inv)
set.seed(1981);Moran.I(GrassDiversity$SpeciesSimpsonDiversity, Grass.dists.inv)
set.seed(1981);Moran.I(GrassDiversity$TotalMinutes, Grass.dists.inv)
set.seed(1981);Moran.I(GrassDiversity$TotalSurveys, Grass.dists.inv)
set.seed(1981);Moran.I(GrassDiversity$ObserverMinutesShannonDiversity, Grass.dists.inv)
set.seed(1981);Moran.I(GrassDiversity$ObserverMinutesSimpsonDiversity, Grass.dists.inv)
summary(GrassDiversity$ObserverRichness)
#HWTC
set.seed(1981);Moran.I(HWTCDiversity$ObserverRichness, HWTC.dists.inv)
set.seed(1981);Moran.I(HWTCDiversity$ObserverShannonDiversity, HWTC.dists.inv)
set.seed(1981);Moran.I(HWTCDiversity$ObserverSimpsonDiversity, HWTC.dists.inv)
set.seed(1981);Moran.I(HWTCDiversity$SpeciesRichness, HWTC.dists.inv)
set.seed(1981);Moran.I(HWTCDiversity$SpeciesShannonDiversity, HWTC.dists.inv)
set.seed(1981);Moran.I(HWTCDiversity$SpeciesSimpsonDiversity, HWTC.dists.inv)
set.seed(1981);Moran.I(HWTCDiversity$TotalMinutes, HWTC.dists.inv)
set.seed(1981);Moran.I(HWTCDiversity$TotalSurveys, HWTC.dists.inv)
set.seed(1981);Moran.I(HWTCDiversity$ObserverMinutesShannonDiversity, HWTC.dists.inv)
set.seed(1981);Moran.I(HWTCDiversity$ObserverMinutesSimpsonDiversity, HWTC.dists.inv)
summary(HWTCDiversity$ObserverRichness)


###########
#Moran's I of Observer Community Distances and Species Diversity
#All Data####
set.seed(1981);Moran.I(AllDiversity$ObserverRichness, AllJaccardObsDistance)
set.seed(1981);Moran.I(AllDiversity$ObserverShannonDiversity, AllJaccardObsDistance)
set.seed(1981);Moran.I(AllDiversity$ObserverSimpsonDiversity, AllJaccardObsDistance)
set.seed(1981);Moran.I(AllDiversity$SpeciesRichness, AllJaccardObsDistance)
set.seed(1981);Moran.I(AllDiversity$SpeciesShannonDiversity, AllJaccardObsDistance)
set.seed(1981);Moran.I(AllDiversity$SpeciesSimpsonDiversity, AllJaccardObsDistance)
set.seed(1981);Moran.I(AllDiversity$TotalMinutes, AllJaccardObsDistance)
set.seed(1981);Moran.I(AllDiversity$TotalSurveys, AllJaccardObsDistance)

set.seed(1981);Moran.I(AllDiversity$ObserverMinutesShannonDiversity, AllJaccardObsDistance)
set.seed(1981);Moran.I(AllDiversity$ObserverMinutesSimpsonDiversity, AllJaccardObsDistance)

summary(AllDiversity$ObserverRichness)

#By Habitat####
#Water

set.seed(1981);Moran.I(WaterDiversity$ObserverRichness, WaterJaccardObsDistance)
set.seed(1981);Moran.I(WaterDiversity$ObserverShannonDiversity, WaterJaccardObsDistance)
set.seed(1981);Moran.I(WaterDiversity$ObserverSimpsonDiversity, WaterJaccardObsDistance)
set.seed(1981);Moran.I(WaterDiversity$SpeciesRichness, WaterJaccardObsDistance)
set.seed(1981);Moran.I(WaterDiversity$SpeciesShannonDiversity, WaterJaccardObsDistance)
set.seed(1981);Moran.I(WaterDiversity$SpeciesSimpsonDiversity, WaterJaccardObsDistance)
set.seed(1981);Moran.I(WaterDiversity$TotalMinutes, WaterJaccardObsDistance)
set.seed(1981);Moran.I(WaterDiversity$TotalSurveys, WaterJaccardObsDistance)
set.seed(1981);Moran.I(WaterDiversity$ObserverMinutesShannonDiversity, WaterJaccardObsDistance)
set.seed(1981);Moran.I(WaterDiversity$ObserverMinutesSimpsonDiversity, WaterJaccardObsDistance)
summary(WaterDiversity$ObserverRichness)

#Pine
set.seed(1981);Moran.I(PineDiversity$ObserverRichness, PineJaccardObsDistance)
set.seed(1981);Moran.I(PineDiversity$ObserverShannonDiversity, PineJaccardObsDistance)
set.seed(1981);Moran.I(PineDiversity$ObserverSimpsonDiversity, PineJaccardObsDistance)
set.seed(1981);Moran.I(PineDiversity$SpeciesRichness, PineJaccardObsDistance)
set.seed(1981);Moran.I(PineDiversity$SpeciesShannonDiversity, PineJaccardObsDistance)
set.seed(1981);Moran.I(PineDiversity$SpeciesSimpsonDiversity, PineJaccardObsDistance)
set.seed(1981);Moran.I(PineDiversity$TotalMinutes, PineJaccardObsDistance)
set.seed(1981);Moran.I(PineDiversity$TotalSurveys, PineJaccardObsDistance)
set.seed(1981);Moran.I(PineDiversity$ObserverMinutesShannonDiversity, PineJaccardObsDistance)
set.seed(1981);Moran.I(PineDiversity$ObserverMinutesSimpsonDiversity, PineJaccardObsDistance)


#Wetland
set.seed(1981);Moran.I(WetlandDiversity$ObserverRichness, WetlandJaccardObsDistance)
set.seed(1981);Moran.I(WetlandDiversity$ObserverShannonDiversity, WetlandJaccardObsDistance)
set.seed(1981);Moran.I(WetlandDiversity$ObserverSimpsonDiversity, WetlandJaccardObsDistance)
set.seed(1981);Moran.I(WetlandDiversity$SpeciesRichness, WetlandJaccardObsDistance)
set.seed(1981);Moran.I(WetlandDiversity$SpeciesShannonDiversity, WetlandJaccardObsDistance)
set.seed(1981);Moran.I(WetlandDiversity$SpeciesSimpsonDiversity, WetlandJaccardObsDistance)
set.seed(1981);Moran.I(WetlandDiversity$TotalMinutes, WetlandJaccardObsDistance)
set.seed(1981);Moran.I(WetlandDiversity$TotalSurveys, WetlandJaccardObsDistance)
set.seed(1981);Moran.I(WetlandDiversity$ObserverMinutesShannonDiversity, WetlandJaccardObsDistance)
set.seed(1981);Moran.I(WetlandDiversity$ObserverMinutesSimpsonDiversity, WetlandJaccardObsDistance)
summary(WetlandDiversity$ObserverRichness)
#Sand
set.seed(1981);Moran.I(SandDiversity$ObserverRichness, SandJaccardObsDistance)
set.seed(1981);Moran.I(SandDiversity$ObserverShannonDiversity, SandJaccardObsDistance)
set.seed(1981);Moran.I(SandDiversity$ObserverSimpsonDiversity, SandJaccardObsDistance)
set.seed(1981);Moran.I(SandDiversity$SpeciesRichness, SandJaccardObsDistance)
set.seed(1981);Moran.I(SandDiversity$SpeciesShannonDiversity, SandJaccardObsDistance)
set.seed(1981);Moran.I(SandDiversity$SpeciesSimpsonDiversity, SandJaccardObsDistance)
set.seed(1981);Moran.I(SandDiversity$TotalMinutes, SandJaccardObsDistance)
set.seed(1981);Moran.I(SandDiversity$TotalSurveys, SandJaccardObsDistance)
set.seed(1981);Moran.I(SandDiversity$ObserverMinutesShannonDiversity, SandJaccardObsDistance)
set.seed(1981);Moran.I(SandDiversity$ObserverMinutesSimpsonDiversity, SandJaccardObsDistance)
summary(SandDiversity$ObserverRichness)
#Urban
set.seed(1981);Moran.I(UrbanDiversity$ObserverRichness, UrbanJaccardObsDistance)
set.seed(1981);Moran.I(UrbanDiversity$ObserverShannonDiversity, UrbanJaccardObsDistance)
set.seed(1981);Moran.I(UrbanDiversity$ObserverSimpsonDiversity, UrbanJaccardObsDistance)
set.seed(1981);Moran.I(UrbanDiversity$SpeciesRichness, UrbanJaccardObsDistance)
set.seed(1981);Moran.I(UrbanDiversity$SpeciesShannonDiversity, UrbanJaccardObsDistance)
set.seed(1981);Moran.I(UrbanDiversity$SpeciesSimpsonDiversity, UrbanJaccardObsDistance)
set.seed(1981);Moran.I(UrbanDiversity$TotalMinutes, UrbanJaccardObsDistance)
set.seed(1981);Moran.I(UrbanDiversity$TotalSurveys, UrbanJaccardObsDistance)
set.seed(1981);Moran.I(UrbanDiversity$ObserverMinutesShannonDiversity, UrbanJaccardObsDistance)
set.seed(1981);Moran.I(UrbanDiversity$ObserverMinutesSimpsonDiversity, UrbanJaccardObsDistance)
summary(UrbanDiversity$ObserverRichness)
#Grass
set.seed(1981);Moran.I(GrassDiversity$ObserverRichness, GrassJaccardObsDistance)
set.seed(1981);Moran.I(GrassDiversity$ObserverShannonDiversity, GrassJaccardObsDistance)
set.seed(1981);Moran.I(GrassDiversity$ObserverSimpsonDiversity, GrassJaccardObsDistance)
set.seed(1981);Moran.I(GrassDiversity$SpeciesRichness, GrassJaccardObsDistance)
set.seed(1981);Moran.I(GrassDiversity$SpeciesShannonDiversity, GrassJaccardObsDistance)
set.seed(1981);Moran.I(GrassDiversity$SpeciesSimpsonDiversity, GrassJaccardObsDistance)
set.seed(1981);Moran.I(GrassDiversity$TotalMinutes, GrassJaccardObsDistance)
set.seed(1981);Moran.I(GrassDiversity$TotalSurveys, GrassJaccardObsDistance)
set.seed(1981);Moran.I(GrassDiversity$ObserverMinutesShannonDiversity, GrassJaccardObsDistance)
set.seed(1981);Moran.I(GrassDiversity$ObserverMinutesSimpsonDiversity, GrassJaccardObsDistance)
summary(GrassDiversity$ObserverRichness)
#HWTC
set.seed(1981);Moran.I(HWTCDiversity$ObserverRichness, HWTCJaccardObsDistance)
set.seed(1981);Moran.I(HWTCDiversity$ObserverShannonDiversity, HWTCJaccardObsDistance)
set.seed(1981);Moran.I(HWTCDiversity$ObserverSimpsonDiversity, HWTCJaccardObsDistance)
set.seed(1981);Moran.I(HWTCDiversity$SpeciesRichness, HWTCJaccardObsDistance)
set.seed(1981);Moran.I(HWTCDiversity$SpeciesShannonDiversity, HWTCJaccardObsDistance)
set.seed(1981);Moran.I(HWTCDiversity$SpeciesSimpsonDiversity, HWTCJaccardObsDistance)
set.seed(1981);Moran.I(HWTCDiversity$TotalMinutes, HWTCJaccardObsDistance)
set.seed(1981);Moran.I(HWTCDiversity$TotalSurveys, HWTCJaccardObsDistance)
set.seed(1981);Moran.I(HWTCDiversity$ObserverMinutesShannonDiversity, HWTCJaccardObsDistance)
set.seed(1981);Moran.I(HWTCDiversity$ObserverMinutesSimpsonDiversity, HWTCJaccardObsDistance)
summary(HWTCDiversity$ObserverRichness)


###########
#Moran's I of Sps Community Distances and Species Diversity
#All Data####
set.seed(1981);Moran.I(AllDiversity$ObserverRichness, AllJaccardSpsDistance)
set.seed(1981);Moran.I(AllDiversity$ObserverShannonDiversity, AllJaccardSpsDistance)
set.seed(1981);Moran.I(AllDiversity$ObserverSimpsonDiversity, AllJaccardSpsDistance)
set.seed(1981);Moran.I(AllDiversity$SpeciesRichness, AllJaccardSpsDistance)
set.seed(1981);Moran.I(AllDiversity$SpeciesShannonDiversity, AllJaccardSpsDistance)
set.seed(1981);Moran.I(AllDiversity$SpeciesSimpsonDiversity, AllJaccardSpsDistance)
set.seed(1981);Moran.I(AllDiversity$TotalMinutes, AllJaccardSpsDistance)
set.seed(1981);Moran.I(AllDiversity$TotalSurveys, AllJaccardSpsDistance)

set.seed(1981);Moran.I(AllDiversity$ObserverMinutesShannonDiversity, AllJaccardSpsDistance)
set.seed(1981);Moran.I(AllDiversity$ObserverMinutesSimpsonDiversity, AllJaccardSpsDistance)

summary(AllDiversity$ObserverRichness)

#By Habitat####
#Water

set.seed(1981);Moran.I(WaterDiversity$ObserverRichness, WaterJaccardSpsDistance)
set.seed(1981);Moran.I(WaterDiversity$ObserverShannonDiversity, WaterJaccardSpsDistance)
set.seed(1981);Moran.I(WaterDiversity$ObserverSimpsonDiversity, WaterJaccardSpsDistance)
set.seed(1981);Moran.I(WaterDiversity$SpeciesRichness, WaterJaccardSpsDistance)
set.seed(1981);Moran.I(WaterDiversity$SpeciesShannonDiversity, WaterJaccardSpsDistance)
set.seed(1981);Moran.I(WaterDiversity$SpeciesSimpsonDiversity, WaterJaccardSpsDistance)
set.seed(1981);Moran.I(WaterDiversity$TotalMinutes, WaterJaccardSpsDistance)
set.seed(1981);Moran.I(WaterDiversity$TotalSurveys, WaterJaccardSpsDistance)
set.seed(1981);Moran.I(WaterDiversity$ObserverMinutesShannonDiversity, WaterJaccardSpsDistance)
set.seed(1981);Moran.I(WaterDiversity$ObserverMinutesSimpsonDiversity, WaterJaccardSpsDistance)
summary(WaterDiversity$ObserverRichness)

#Pine
set.seed(1981);Moran.I(PineDiversity$ObserverRichness, PineJaccardSpsDistance)
set.seed(1981);Moran.I(PineDiversity$ObserverShannonDiversity, PineJaccardSpsDistance)
set.seed(1981);Moran.I(PineDiversity$ObserverSimpsonDiversity, PineJaccardSpsDistance)
set.seed(1981);Moran.I(PineDiversity$SpeciesRichness, PineJaccardSpsDistance)
set.seed(1981);Moran.I(PineDiversity$SpeciesShannonDiversity, PineJaccardSpsDistance)
set.seed(1981);Moran.I(PineDiversity$SpeciesSimpsonDiversity, PineJaccardSpsDistance)
set.seed(1981);Moran.I(PineDiversity$TotalMinutes, PineJaccardSpsDistance)
set.seed(1981);Moran.I(PineDiversity$TotalSurveys, PineJaccardSpsDistance)
set.seed(1981);Moran.I(PineDiversity$ObserverMinutesShannonDiversity, PineJaccardSpsDistance)
set.seed(1981);Moran.I(PineDiversity$ObserverMinutesSimpsonDiversity, PineJaccardSpsDistance)


#Wetland
set.seed(1981);Moran.I(WetlandDiversity$ObserverRichness, WetlandJaccardSpsDistance)
set.seed(1981);Moran.I(WetlandDiversity$ObserverShannonDiversity, WetlandJaccardSpsDistance)
set.seed(1981);Moran.I(WetlandDiversity$ObserverSimpsonDiversity, WetlandJaccardSpsDistance)
set.seed(1981);Moran.I(WetlandDiversity$SpeciesRichness, WetlandJaccardSpsDistance)
set.seed(1981);Moran.I(WetlandDiversity$SpeciesShannonDiversity, WetlandJaccardSpsDistance)
set.seed(1981);Moran.I(WetlandDiversity$SpeciesSimpsonDiversity, WetlandJaccardSpsDistance)
set.seed(1981);Moran.I(WetlandDiversity$TotalMinutes, WetlandJaccardSpsDistance)
set.seed(1981);Moran.I(WetlandDiversity$TotalSurveys, WetlandJaccardSpsDistance)
set.seed(1981);Moran.I(WetlandDiversity$ObserverMinutesShannonDiversity, WetlandJaccardSpsDistance)
set.seed(1981);Moran.I(WetlandDiversity$ObserverMinutesSimpsonDiversity, WetlandJaccardSpsDistance)
summary(WetlandDiversity$ObserverRichness)
#Sand
set.seed(1981);Moran.I(SandDiversity$ObserverRichness, SandJaccardSpsDistance)
set.seed(1981);Moran.I(SandDiversity$ObserverShannonDiversity, SandJaccardSpsDistance)
set.seed(1981);Moran.I(SandDiversity$ObserverSimpsonDiversity, SandJaccardSpsDistance)
set.seed(1981);Moran.I(SandDiversity$SpeciesRichness, SandJaccardSpsDistance)
set.seed(1981);Moran.I(SandDiversity$SpeciesShannonDiversity, SandJaccardSpsDistance)
set.seed(1981);Moran.I(SandDiversity$SpeciesSimpsonDiversity, SandJaccardSpsDistance)
set.seed(1981);Moran.I(SandDiversity$TotalMinutes, SandJaccardSpsDistance)
set.seed(1981);Moran.I(SandDiversity$TotalSurveys, SandJaccardSpsDistance)
set.seed(1981);Moran.I(SandDiversity$ObserverMinutesShannonDiversity, SandJaccardSpsDistance)
set.seed(1981);Moran.I(SandDiversity$ObserverMinutesSimpsonDiversity, SandJaccardSpsDistance)
summary(SandDiversity$ObserverRichness)
#Urban
set.seed(1981);Moran.I(UrbanDiversity$ObserverRichness, UrbanJaccardSpsDistance)
set.seed(1981);Moran.I(UrbanDiversity$ObserverShannonDiversity, UrbanJaccardSpsDistance)
set.seed(1981);Moran.I(UrbanDiversity$ObserverSimpsonDiversity, UrbanJaccardSpsDistance)
set.seed(1981);Moran.I(UrbanDiversity$SpeciesRichness, UrbanJaccardSpsDistance)
set.seed(1981);Moran.I(UrbanDiversity$SpeciesShannonDiversity, UrbanJaccardSpsDistance)
set.seed(1981);Moran.I(UrbanDiversity$SpeciesSimpsonDiversity, UrbanJaccardSpsDistance)
set.seed(1981);Moran.I(UrbanDiversity$TotalMinutes, UrbanJaccardSpsDistance)
set.seed(1981);Moran.I(UrbanDiversity$TotalSurveys, UrbanJaccardSpsDistance)
set.seed(1981);Moran.I(UrbanDiversity$ObserverMinutesShannonDiversity, UrbanJaccardSpsDistance)
set.seed(1981);Moran.I(UrbanDiversity$ObserverMinutesSimpsonDiversity, UrbanJaccardSpsDistance)
summary(UrbanDiversity$ObserverRichness)
#Grass
set.seed(1981);Moran.I(GrassDiversity$ObserverRichness, GrassJaccardSpsDistance)
set.seed(1981);Moran.I(GrassDiversity$ObserverShannonDiversity, GrassJaccardSpsDistance)
set.seed(1981);Moran.I(GrassDiversity$ObserverSimpsonDiversity, GrassJaccardSpsDistance)
set.seed(1981);Moran.I(GrassDiversity$SpeciesRichness, GrassJaccardSpsDistance)
set.seed(1981);Moran.I(GrassDiversity$SpeciesShannonDiversity, GrassJaccardSpsDistance)
set.seed(1981);Moran.I(GrassDiversity$SpeciesSimpsonDiversity, GrassJaccardSpsDistance)
set.seed(1981);Moran.I(GrassDiversity$TotalMinutes, GrassJaccardSpsDistance)
set.seed(1981);Moran.I(GrassDiversity$TotalSurveys, GrassJaccardSpsDistance)
set.seed(1981);Moran.I(GrassDiversity$ObserverMinutesShannonDiversity, GrassJaccardSpsDistance)
set.seed(1981);Moran.I(GrassDiversity$ObserverMinutesSimpsonDiversity, GrassJaccardSpsDistance)
summary(GrassDiversity$ObserverRichness)
#HWTC
set.seed(1981);Moran.I(HWTCDiversity$ObserverRichness, HWTCJaccardSpsDistance)
set.seed(1981);Moran.I(HWTCDiversity$ObserverShannonDiversity, HWTCJaccardSpsDistance)
set.seed(1981);Moran.I(HWTCDiversity$ObserverSimpsonDiversity, HWTCJaccardSpsDistance)
set.seed(1981);Moran.I(HWTCDiversity$SpeciesRichness, HWTCJaccardSpsDistance)
set.seed(1981);Moran.I(HWTCDiversity$SpeciesShannonDiversity, HWTCJaccardSpsDistance)
set.seed(1981);Moran.I(HWTCDiversity$SpeciesSimpsonDiversity, HWTCJaccardSpsDistance)
set.seed(1981);Moran.I(HWTCDiversity$TotalMinutes, HWTCJaccardSpsDistance)
set.seed(1981);Moran.I(HWTCDiversity$TotalSurveys, HWTCJaccardSpsDistance)
set.seed(1981);Moran.I(HWTCDiversity$ObserverMinutesShannonDiversity, HWTCJaccardSpsDistance)
set.seed(1981);Moran.I(HWTCDiversity$ObserverMinutesSimpsonDiversity, HWTCJaccardSpsDistance)
summary(HWTCDiversity$ObserverRichness)

###########
#Mantel Tests of Observer and Species Community Indices####
#Jaccard Observer Distance Matrices####
#Create a function to complete the jaccard function on all pairs of columns in selected data#####
pairedcolumns <- function(x,fun) #This function works on a data frame x using whichever other function you select
{
  n <- ncol(x)##find out how many columns are in the data frame
  
  foo <- matrix(0,n,n)
  for ( i in 1:n)
  {
    for (j in 1:n)
    {
      foo[i,j] <- fun(x[,i],x[,j])
    }
  }
  colnames(foo)<-rownames(foo)<-colnames(x)
  return(foo)
}

#Create columnwise data on How many surveys each observers conducted at each Locality####
AllObserverSurveys<-dcast(GBeBirdDEC2016,
                          `OBSERVER ID`~`LOCALITY ID`, value.var ="SAMPLING EVENT IDENTIFIER" )
WaterObserverSurveys<-dcast(GBeBirdDEC2016[GBeBirdDEC2016$Habitat=="Water",],
                            `OBSERVER ID`~`LOCALITY ID`, value.var ="SAMPLING EVENT IDENTIFIER" )
PineObserverSurveys<-dcast(GBeBirdDEC2016[GBeBirdDEC2016$Habitat=="Pine",],
                           `OBSERVER ID`~`LOCALITY ID`, value.var ="SAMPLING EVENT IDENTIFIER" )
WetlandObserverSurveys<-dcast(GBeBirdDEC2016[GBeBirdDEC2016$Habitat=="Wetland",],
                              `OBSERVER ID`~`LOCALITY ID`, value.var ="SAMPLING EVENT IDENTIFIER" )
SandObserverSurveys<-dcast(GBeBirdDEC2016[GBeBirdDEC2016$Habitat=="Sand",],
                           `OBSERVER ID`~`LOCALITY ID`, value.var ="SAMPLING EVENT IDENTIFIER" )
UrbanObserverSurveys<-dcast(GBeBirdDEC2016[GBeBirdDEC2016$Habitat=="Urban",],
                            `OBSERVER ID`~`LOCALITY ID`, value.var ="SAMPLING EVENT IDENTIFIER" )
GrassObserverSurveys<-dcast(GBeBirdDEC2016[GBeBirdDEC2016$Habitat=="Grass",],
                            `OBSERVER ID`~`LOCALITY ID`, value.var ="SAMPLING EVENT IDENTIFIER" )
HWTCObserverSurveys<-dcast(GBeBirdDEC2016[GBeBirdDEC2016$Habitat=="HWTC",],
                           `OBSERVER ID`~`LOCALITY ID`, value.var ="SAMPLING EVENT IDENTIFIER" )
#Create Columnwise data on whether each observer visited each locality####

AllObserverVisits<-AllObserverSurveys
AllObserverVisits[,2:ncol(AllObserverVisits)]<-AllObserverVisits[,2:ncol(AllObserverVisits)]>0
#Water
WaterObserverVisits<-WaterObserverSurveys
WaterObserverVisits[,2:ncol(WaterObserverVisits)]<-WaterObserverVisits[,2:ncol(WaterObserverVisits)]>0
#Pine
PineObserverVisits<-PineObserverSurveys
PineObserverVisits[,2:ncol(PineObserverVisits)]<-PineObserverVisits[,2:ncol(PineObserverVisits)]>0
#Wetland
WetlandObserverVisits<-WetlandObserverSurveys
WetlandObserverVisits[,2:ncol(WetlandObserverVisits)]<-WetlandObserverVisits[,2:ncol(WetlandObserverVisits)]>0
#Sand
SandObserverVisits<-SandObserverSurveys
SandObserverVisits[,2:ncol(SandObserverVisits)]<-SandObserverVisits[,2:ncol(SandObserverVisits)]>0
#Urban
UrbanObserverVisits<-UrbanObserverSurveys
UrbanObserverVisits[,2:ncol(UrbanObserverVisits)]<-UrbanObserverVisits[,2:ncol(UrbanObserverVisits)]>0
#Grass
GrassObserverVisits<-GrassObserverSurveys
GrassObserverVisits[,2:ncol(GrassObserverVisits)]<-GrassObserverVisits[,2:ncol(GrassObserverVisits)]>0
#HWTC
HWTCObserverVisits<-HWTCObserverSurveys
HWTCObserverVisits[,2:ncol(HWTCObserverVisits)]<-HWTCObserverVisits[,2:ncol(HWTCObserverVisits)]>0

#Create Jaccard Observer Similarity Matrix for All data and each habitat####
AllJaccardObsSimilarity<-pairedcolumns(AllObserverVisits[,2:ncol(AllObserverVisits)], clujaccard)
#Water
WaterJaccardObsSimilarity<-pairedcolumns(WaterObserverVisits[,2:ncol(WaterObserverVisits)], clujaccard)
#Pine
PineJaccardObsSimilarity<-pairedcolumns(PineObserverVisits[,2:ncol(PineObserverVisits)], clujaccard)
#Wetland
WetlandJaccardObsSimilarity<-pairedcolumns(WetlandObserverVisits[,2:ncol(WetlandObserverVisits)], clujaccard)
#Sand
SandJaccardObsSimilarity<-pairedcolumns(SandObserverVisits[,2:ncol(SandObserverVisits)], clujaccard)
#Urban
UrbanJaccardObsSimilarity<-pairedcolumns(UrbanObserverVisits[,2:ncol(UrbanObserverVisits)], clujaccard)
#Grass
GrassJaccardObsSimilarity<-pairedcolumns(GrassObserverVisits[,2:ncol(GrassObserverVisits)], clujaccard)
#HWTC
HWTCJaccardObsSimilarity<-pairedcolumns(HWTCObserverVisits[,2:ncol(HWTCObserverVisits)], clujaccard)  

#Create Jaccard Observer Difference Matrix for All data and each habitat####
AllJaccardObsDistance<-1- AllJaccardObsSimilarity
#Water
WaterJaccardObsDistance<-1-WaterJaccardObsSimilarity
#Pine
PineJaccardObsDistance<-1-PineJaccardObsSimilarity
#Wetland
WetlandJaccardObsDistance<-1-WetlandJaccardObsSimilarity
#Sand
SandJaccardObsDistance<-1-SandJaccardObsSimilarity
#Urban
UrbanJaccardObsDistance<-1-UrbanJaccardObsSimilarity
#Grass
GrassJaccardObsDistance<-1-GrassJaccardObsSimilarity
#HWTC
HWTCJaccardObsDistance<-1-HWTCJaccardObsSimilarity  
##Create Coordinates in Observer Space####
WaterObscoords<-as.data.frame(cmdscale(WaterJaccardObsDistance)) #THis generates 2 dimensional coordinates for the observer communities.
plot(WaterObscoords)
PineObscoords<-as.data.frame(cmdscale(PineJaccardObsDistance)) #THis generates 2 dimensional coordinates for the observer communities.
plot(PineObscoords)
WetlandObscoords<-as.data.frame(cmdscale(WetlandJaccardObsDistance)) #THis generates 2 dimensional coordinates for the observer communities.
plot(WetlandObscoords)
SandObscoords<-as.data.frame(cmdscale(SandJaccardObsDistance)) #THis generates 2 dimensional coordinates for the observer communities.
plot(SandObscoords)
UrbanObscoords<-as.data.frame(cmdscale(UrbanJaccardObsDistance)) #THis generates 2 dimensional coordinates for the observer communities.
plot(UrbanObscoords)
GrassObscoords<-as.data.frame(cmdscale(GrassJaccardObsDistance)) #THis generates 2 dimensional coordinates for the observer communities.
plot(GrassObscoords)
HWTCObscoords<-as.data.frame(cmdscale(HWTCJaccardObsDistance)) #THis generates 2 dimensional coordinates for the observer communities.
plot(HWTCObscoords)
AllObscoords<-as.data.frame(cmdscale(AllJaccardObsDistance)) #THis generates 2 dimensional coordinates for the observer communities.
plot(AllObscoords, main= "Observer community space in eBird Grand Bahama 1988-2016")
AllObscoords$LOCALITY.ID<-rownames(AllObscoords)
colnames(AllObscoords)<-c("Obs.x.coord", "Obs.y.coord", "LOCALITY.ID")



#Jaccard Species Distance Matrices####

#Create columnwise data on How many surveys each Species occurred in at each Locality####
AllSpeciesSurveys<-dcast(GBeBirdDEC2016,
                         `COMMON NAME`~`LOCALITY ID`, value.var ="SAMPLING EVENT IDENTIFIER" )
WaterSpeciesSurveys<-dcast(GBeBirdDEC2016[GBeBirdDEC2016$Habitat=="Water",],
                           `COMMON NAME`~`LOCALITY ID`, value.var ="SAMPLING EVENT IDENTIFIER" )
PineSpeciesSurveys<-dcast(GBeBirdDEC2016[GBeBirdDEC2016$Habitat=="Pine",],
                          `COMMON NAME`~`LOCALITY ID`, value.var ="SAMPLING EVENT IDENTIFIER" )
WetlandSpeciesSurveys<-dcast(GBeBirdDEC2016[GBeBirdDEC2016$Habitat=="Wetland",],
                             `COMMON NAME`~`LOCALITY ID`, value.var ="SAMPLING EVENT IDENTIFIER" )
SandSpeciesSurveys<-dcast(GBeBirdDEC2016[GBeBirdDEC2016$Habitat=="Sand",],
                          `COMMON NAME`~`LOCALITY ID`, value.var ="SAMPLING EVENT IDENTIFIER" )
UrbanSpeciesSurveys<-dcast(GBeBirdDEC2016[GBeBirdDEC2016$Habitat=="Urban",],
                           `COMMON NAME`~`LOCALITY ID`, value.var ="SAMPLING EVENT IDENTIFIER" )
GrassSpeciesSurveys<-dcast(GBeBirdDEC2016[GBeBirdDEC2016$Habitat=="Grass",],
                           `COMMON NAME`~`LOCALITY ID`, value.var ="SAMPLING EVENT IDENTIFIER" )
HWTCSpeciesSurveys<-dcast(GBeBirdDEC2016[GBeBirdDEC2016$Habitat=="HWTC",],
                          `COMMON NAME`~`LOCALITY ID`, value.var ="SAMPLING EVENT IDENTIFIER" )
#Create Columnwise data on whether each Species visited each locality####

AllSpeciesVisits<-AllSpeciesSurveys
AllSpeciesVisits[,2:ncol(AllSpeciesVisits)]<-AllSpeciesVisits[,2:ncol(AllSpeciesVisits)]>0
#Water
WaterSpeciesVisits<-WaterSpeciesSurveys
WaterSpeciesVisits[,2:ncol(WaterSpeciesVisits)]<-WaterSpeciesVisits[,2:ncol(WaterSpeciesVisits)]>0
#Pine
PineSpeciesVisits<-PineSpeciesSurveys
PineSpeciesVisits[,2:ncol(PineSpeciesVisits)]<-PineSpeciesVisits[,2:ncol(PineSpeciesVisits)]>0
#Wetland
WetlandSpeciesVisits<-WetlandSpeciesSurveys
WetlandSpeciesVisits[,2:ncol(WetlandSpeciesVisits)]<-WetlandSpeciesVisits[,2:ncol(WetlandSpeciesVisits)]>0
#Sand
SandSpeciesVisits<-SandSpeciesSurveys
SandSpeciesVisits[,2:ncol(SandSpeciesVisits)]<-SandSpeciesVisits[,2:ncol(SandSpeciesVisits)]>0
#Urban
UrbanSpeciesVisits<-UrbanSpeciesSurveys
UrbanSpeciesVisits[,2:ncol(UrbanSpeciesVisits)]<-UrbanSpeciesVisits[,2:ncol(UrbanSpeciesVisits)]>0
#Grass
GrassSpeciesVisits<-GrassSpeciesSurveys
GrassSpeciesVisits[,2:ncol(GrassSpeciesVisits)]<-GrassSpeciesVisits[,2:ncol(GrassSpeciesVisits)]>0
#HWTC
HWTCSpeciesVisits<-HWTCSpeciesSurveys
HWTCSpeciesVisits[,2:ncol(HWTCSpeciesVisits)]<-HWTCSpeciesVisits[,2:ncol(HWTCSpeciesVisits)]>0

#Create Jaccard Species Similarity Matrix for All data and each habitat####
AllJaccardSpsSimilarity<-pairedcolumns(AllSpeciesVisits[,2:ncol(AllSpeciesVisits)], clujaccard)
#Water
WaterJaccardSpsSimilarity<-pairedcolumns(WaterSpeciesVisits[,2:ncol(WaterSpeciesVisits)], clujaccard)
#Pine
PineJaccardSpsSimilarity<-pairedcolumns(PineSpeciesVisits[,2:ncol(PineSpeciesVisits)], clujaccard)
#Wetland
WetlandJaccardSpsSimilarity<-pairedcolumns(WetlandSpeciesVisits[,2:ncol(WetlandSpeciesVisits)], clujaccard)
#Sand
SandJaccardSpsSimilarity<-pairedcolumns(SandSpeciesVisits[,2:ncol(SandSpeciesVisits)], clujaccard)
#Urban
UrbanJaccardSpsSimilarity<-pairedcolumns(UrbanSpeciesVisits[,2:ncol(UrbanSpeciesVisits)], clujaccard)
#Grass
GrassJaccardSpsSimilarity<-pairedcolumns(GrassSpeciesVisits[,2:ncol(GrassSpeciesVisits)], clujaccard)
#HWTC
HWTCJaccardSpsSimilarity<-pairedcolumns(HWTCSpeciesVisits[,2:ncol(HWTCSpeciesVisits)], clujaccard)  


#Create Jaccard Species Difference Matrix for All data and each habitat####
AllJaccardSpsDistance<-1- AllJaccardSpsSimilarity
#Water
WaterJaccardSpsDistance<-1-WaterJaccardSpsSimilarity
#Pine
PineJaccardSpsDistance<-1-PineJaccardSpsSimilarity
#Wetland
WetlandJaccardSpsDistance<-1-WetlandJaccardSpsSimilarity
#Sand
SandJaccardSpsDistance<-1-SandJaccardSpsSimilarity
#Urban
UrbanJaccardSpsDistance<-1-UrbanJaccardSpsSimilarity
#Grass
GrassJaccardSpsDistance<-1-GrassJaccardSpsSimilarity
#HWTC
HWTCJaccardSpsDistance<-1-HWTCJaccardSpsSimilarity  

##Create Coordinates in Species Space####
WaterSpscoords<-as.data.frame(cmdscale(WaterJaccardSpsDistance)) #THis generates 2 dimensional coordinates for the Spserver communities.
plot(WaterSpscoords)
PineSpscoords<-as.data.frame(cmdscale(PineJaccardSpsDistance)) #THis generates 2 dimensional coordinates for the Spserver communities.
plot(PineSpscoords)
WetlandSpscoords<-as.data.frame(cmdscale(WetlandJaccardSpsDistance)) #THis generates 2 dimensional coordinates for the Spserver communities.
plot(WetlandSpscoords)
SandSpscoords<-as.data.frame(cmdscale(SandJaccardSpsDistance)) #THis generates 2 dimensional coordinates for the Spserver communities.
plot(SandSpscoords)
UrbanSpscoords<-as.data.frame(cmdscale(UrbanJaccardSpsDistance)) #THis generates 2 dimensional coordinates for the Spserver communities.
plot(UrbanSpscoords)
GrassSpscoords<-as.data.frame(cmdscale(GrassJaccardSpsDistance)) #THis generates 2 dimensional coordinates for the Spserver communities.
plot(main = "Grass habitat species community space in eBird data for Grand Bahama Island 1988-2016",
     GrassSpscoords)



HWTCSpscoords<-as.data.frame(cmdscale(HWTCJaccardSpsDistance)) #THis generates 2 dimensional coordinates for the Spserver communities.
plot(HWTCSpscoords)

AllSpscoords<-as.data.frame(cmdscale(AllJaccardSpsDistance))#THis generates 2 dimensional coordinates for the Species data communities.
AllSpscoords$LOCALITY.ID<-rownames(AllSpscoords)

colnames(AllSpscoords)<-c("Sps.x.coord", "Sps.y.coord", "LOCALITY.ID")
Diversityin3spaces<-merge(merge(AllDiversity, AllSpscoords), AllObscoords)

#Create matrix of locality id and coordinates in Observer and Species Space
LocalitynamesObs<-cbind(as.matrix(rownames(AllObscoords)), AllObscoords)
LocalitynamesSps<-cbind(as.matrix(rownames(AllSpscoords)), AllSpscoords)
colnames(LocalitynamesObs)<-c("LOCALITY.ID", "obs.x.coord", "obs.y.coord")
  colnames(LocalitynamesSps)<-c("LOCALITY.ID", "sps.x.coord", "sps.y.coord")
AllDiversity
DiversityCoordinates<-merge(merge(AllDiversity, LocalitynamesObs), LocalitynamesSps)

#Mantel Test of Geographic Distance and Observer Similarity####
#Spearman
set.seed(1981);mantel(All.dists.inv, AllJaccardObsDistance, method = "spearman")
set.seed(1981);mantel(Water.dists.inv, WaterJaccardObsDistance, method = "spearman")
set.seed(1981);mantel(Pine.dists.inv, PineJaccardObsDistance, method = "spearman")
set.seed(1981);mantel(Wetland.dists.inv, WetlandJaccardObsDistance, method = "spearman")
set.seed(1981);mantel(Sand.dists.inv, SandJaccardObsDistance, method = "spearman")
set.seed(1981);mantel(Urban.dists.inv, UrbanJaccardObsDistance, method = "spearman")
set.seed(1981);mantel(Grass.dists.inv, GrassJaccardObsDistance, method = "spearman")
set.seed(1981);mantel(HWTC.dists.inv, HWTCJaccardObsDistance, method = "spearman")
#Pearson
set.seed(1981);mantel(All.dists.inv, AllJaccardObsDistance)
set.seed(1981);mantel(Water.dists.inv, WaterJaccardObsDistance)
set.seed(1981);mantel(Pine.dists.inv, PineJaccardObsDistance)
set.seed(1981);mantel(Wetland.dists.inv, WetlandJaccardObsDistance)
set.seed(1981);mantel(Sand.dists.inv, SandJaccardObsDistance)
set.seed(1981);mantel(Urban.dists.inv, UrbanJaccardObsDistance)
set.seed(1981);mantel(Grass.dists.inv, GrassJaccardObsDistance)
set.seed(1981);mantel(HWTC.dists.inv, HWTCJaccardObsDistance)

#Mantel Test of Geographic Distance and Species Similarity####
#pearson
set.seed(1981);mantel(All.dists.inv, AllJaccardSpsDistance)
set.seed(1981);mantel(Water.dists.inv, WaterJaccardSpsDistance)
set.seed(1981);mantel(Pine.dists.inv, PineJaccardSpsDistance)
set.seed(1981);mantel(Wetland.dists.inv, WetlandJaccardSpsDistance)
set.seed(1981);mantel(Sand.dists.inv, SandJaccardSpsDistance)
set.seed(1981);mantel(Urban.dists.inv, UrbanJaccardSpsDistance)
set.seed(1981);mantel(Grass.dists.inv, GrassJaccardSpsDistance)
set.seed(1981);mantel(HWTC.dists.inv, HWTCJaccardSpsDistance)
#Spearman
set.seed(1981);mantel(All.dists.inv, AllJaccardSpsDistance, method = "spearman")
set.seed(1981);mantel(Water.dists.inv, WaterJaccardSpsDistance, method = "spearman")
set.seed(1981);mantel(Pine.dists.inv, PineJaccardSpsDistance, method = "spearman")
set.seed(1981);mantel(Wetland.dists.inv, WetlandJaccardSpsDistance, method = "spearman")
set.seed(1981);mantel(Sand.dists.inv, SandJaccardSpsDistance, method = "spearman")
set.seed(1981);mantel(Urban.dists.inv, UrbanJaccardSpsDistance, method = "spearman")
set.seed(1981);mantel(Grass.dists.inv, GrassJaccardSpsDistance, method = "spearman")
set.seed(1981);mantel(HWTC.dists.inv, HWTCJaccardSpsDistance, method = "spearman")

#Mantel Test of Observer and Species Similarity####
#Pearson
set.seed(1981);mantel(AllJaccardObsDistance, AllJaccardSpsDistance)
set.seed(1981);mantel(WaterJaccardObsDistance, WaterJaccardSpsDistance)
set.seed(1981);mantel(PineJaccardObsDistance, PineJaccardSpsDistance)
set.seed(1981);mantel(WetlandJaccardObsDistance, WetlandJaccardSpsDistance)
set.seed(1981);mantel(SandJaccardObsDistance, SandJaccardSpsDistance)
set.seed(1981);mantel(UrbanJaccardObsDistance, UrbanJaccardSpsDistance)
set.seed(1981);mantel(GrassJaccardObsDistance, GrassJaccardSpsDistance)
set.seed(1981);mantel(HWTCJaccardObsDistance, HWTCJaccardSpsDistance)

#Spearman
set.seed(1981);mantel(AllJaccardObsDistance, AllJaccardSpsDistance, method = "spearman")
set.seed(1981);mantel(WaterJaccardObsDistance, WaterJaccardSpsDistance, method = "spearman")
set.seed(1981);mantel(PineJaccardObsDistance, PineJaccardSpsDistance, method = "spearman")
set.seed(1981);mantel(WetlandJaccardObsDistance, WetlandJaccardSpsDistance, method = "spearman")
set.seed(1981);mantel(SandJaccardObsDistance, SandJaccardSpsDistance, method = "spearman")
set.seed(1981);mantel(UrbanJaccardObsDistance, UrbanJaccardSpsDistance, method = "spearman")
set.seed(1981);mantel(GrassJaccardObsDistance, GrassJaccardSpsDistance, method = "spearman")
set.seed(1981);mantel(HWTCJaccardObsDistance, HWTCJaccardSpsDistance, method = "spearman")

#Create data frames with species in columns and the number of surveys they were found in per locality####
#All HabitatTypes
AllHabitatSpsSurveysAtLocality<-t(AllSpeciesSurveys) #transpose the locality surveys by species
colnames(AllHabitatSpsSurveysAtLocality)<-AllHabitatSpsSurveysAtLocality[1,]# use the row of species names for column names
AllHabitatSpsSurveysAtLocality<-AllHabitatSpsSurveysAtLocality[-1,] #remove the row of species names

#Water HabitatTypes
WaterHabitatSpsSurveysAtLocality<-as.data.frame(t(WaterSpeciesSurveys)) #transpose the locality surveys by species
colnames(WaterHabitatSpsSurveysAtLocality)<-WaterHabitatSpsSurveysAtLocality[1,]# use the row of species names for column names
WaterHabitatSpsSurveysAtLocality<-WaterHabitatSpsSurveysAtLocality[-1,] #remove the row of species names

#Pine HabitatTypes
PineHabitatSpsSurveysAtLocality<-t(PineSpeciesSurveys) #transpose the locality surveys by species
colnames(PineHabitatSpsSurveysAtLocality)<-PineHabitatSpsSurveysAtLocality[1,]# use the row of species names for column names
PineHabitatSpsSurveysAtLocality<-PineHabitatSpsSurveysAtLocality[-1,] #remove the row of species names

#Wetland HabitatTypes
WetlandHabitatSpsSurveysAtLocality<-t(WetlandSpeciesSurveys) #transpose the locality surveys by species
colnames(WetlandHabitatSpsSurveysAtLocality)<-WetlandHabitatSpsSurveysAtLocality[1,]# use the row of species names for column names
WetlandHabitatSpsSurveysAtLocality<-WetlandHabitatSpsSurveysAtLocality[-1,] #remove the row of species names

#Sand HabitatTypes
SandHabitatSpsSurveysAtLocality<-t(SandSpeciesSurveys) #transpose the locality surveys by species
colnames(SandHabitatSpsSurveysAtLocality)<-SandHabitatSpsSurveysAtLocality[1,]# use the row of species names for column names
SandHabitatSpsSurveysAtLocality<-SandHabitatSpsSurveysAtLocality[-1,] #remove the row of species names

#Urban HabitatTypes
UrbanHabitatSpsSurveysAtLocality<-t(UrbanSpeciesSurveys) #transpose the locality surveys by species
colnames(UrbanHabitatSpsSurveysAtLocality)<-UrbanHabitatSpsSurveysAtLocality[1,]# use the row of species names for column names
UrbanHabitatSpsSurveysAtLocality<-UrbanHabitatSpsSurveysAtLocality[-1,] #remove the row of species names

#Grass HabitatTypes
GrassHabitatSpsSurveysAtLocality<-t(GrassSpeciesSurveys) #transpose the locality surveys by species
colnames(GrassHabitatSpsSurveysAtLocality)<-GrassHabitatSpsSurveysAtLocality[1,]# use the row of species names for column names
GrassHabitatSpsSurveysAtLocality<-GrassHabitatSpsSurveysAtLocality[-1,] #remove the row of species names

#HWTC HabitatTypes
HWTCHabitatSpsSurveysAtLocality<-t(HWTCSpeciesSurveys) #transpose the locality surveys by species
colnames(HWTCHabitatSpsSurveysAtLocality)<-HWTCHabitatSpsSurveysAtLocality[1,]# use the row of species names for column names
HWTCHabitatSpsSurveysAtLocality<-HWTCHabitatSpsSurveysAtLocality[-1,] #remove the row of species names





#Create data with Species in Columns and whether or not they were found at a locality or by an observer####
SpsAtAllLoc<-t(AllSpeciesVisits)
SpsAtWaterLoc<-t(WaterSpeciesVisits)
SpsAtPineLoc<-t(PineSpeciesVisits)
SpsAtWetlandLoc<-t(WetlandSpeciesVisits)
SpsAtSandLoc<-t(SandSpeciesVisits)
SpsAtUrbanLoc<-t(UrbanSpeciesVisits)
SpsAtGrassLoc<-t(GrassSpeciesVisits)
SpsAtHWTCLoc<-t(HWTCSpeciesVisits)

#Create data with species in columns and whether or not they were seen by an observer####
AllSpsObsSurveys<-dcast(GBeBirdDEC2016,
                         `OBSERVER ID`~`COMMON NAME`, value.var ="SAMPLING EVENT IDENTIFIER" )
WaterSpsObsSurveys<-dcast(GBeBirdDEC2016[GBeBirdDEC2016$Habitat=="Water",],
                           `OBSERVER ID`~`COMMON NAME`, value.var ="SAMPLING EVENT IDENTIFIER" )
PineSpsObsSurveys<-dcast(GBeBirdDEC2016[GBeBirdDEC2016$Habitat=="Pine",],
                          `OBSERVER ID`~`COMMON NAME`, value.var ="SAMPLING EVENT IDENTIFIER" )
WetlandSpsObsSurveys<-dcast(GBeBirdDEC2016[GBeBirdDEC2016$Habitat=="Wetland",],
                             `OBSERVER ID`~`COMMON NAME`, value.var ="SAMPLING EVENT IDENTIFIER" )
SandSpsObsSurveys<-dcast(GBeBirdDEC2016[GBeBirdDEC2016$Habitat=="Sand",],
                          `OBSERVER ID`~`COMMON NAME`, value.var ="SAMPLING EVENT IDENTIFIER" )
UrbanSpsObsSurveys<-dcast(GBeBirdDEC2016[GBeBirdDEC2016$Habitat=="Urban",],
                           `OBSERVER ID`~`COMMON NAME`, value.var ="SAMPLING EVENT IDENTIFIER" )
GrassSpsObsSurveys<-dcast(GBeBirdDEC2016[GBeBirdDEC2016$Habitat=="Grass",],
                           `OBSERVER ID`~`COMMON NAME`, value.var ="SAMPLING EVENT IDENTIFIER" )
HWTCSpsObsSurveys<-dcast(GBeBirdDEC2016[GBeBirdDEC2016$Habitat=="HWTC",],
                          `OBSERVER ID`~`COMMON NAME`, value.var ="SAMPLING EVENT IDENTIFIER" )
rownames(AllSpsObsSurveys)<-AllSpsObsSurveys[,1]; AllSpsObsSurveys<-AllSpsObsSurveys[,-1]
rownames(WaterSpsObsSurveys)<-WaterSpsObsSurveys[,1]; WaterSpsObsSurveys<-WaterSpsObsSurveys[,-1]
rownames(PineSpsObsSurveys)<-PineSpsObsSurveys[,1]; PineSpsObsSurveys<-PineSpsObsSurveys[,-1]
rownames(WetlandSpsObsSurveys)<-WetlandSpsObsSurveys[,1]; WetlandSpsObsSurveys<-WetlandSpsObsSurveys[,-1]
rownames(SandSpsObsSurveys)<-SandSpsObsSurveys[,1]; SandSpsObsSurveys<-SandSpsObsSurveys[,-1]
rownames(UrbanSpsObsSurveys)<-UrbanSpsObsSurveys[,1]; UrbanSpsObsSurveys<-UrbanSpsObsSurveys[,-1]
rownames(GrassSpsObsSurveys)<-GrassSpsObsSurveys[,1]; GrassSpsObsSurveys<-GrassSpsObsSurveys[,-1]
rownames(HWTCSpsObsSurveys)<-HWTCSpsObsSurveys[,1]; HWTCSpsObsSurveys<-HWTCSpsObsSurveys[,-1]

AllSpsObsSeen<-AllSpsObsSurveys>0
WaterSpsObsSeen<-WaterSpsObsSurveys>0
PineSpsObsSeen<-PineSpsObsSurveys>0
WetlandSpsObsSeen<-WetlandSpsObsSurveys>0
SandSpsObsSeen<-SandSpsObsSurveys>0
UrbanSpsObsSeen<-UrbanSpsObsSurveys>0
GrassSpsObsSeen<-GrassSpsObsSurveys>0
HWTCSpsObsSeen<-HWTCSpsObsSurveys>0

##Species names need to be made without spaces for this to work with species names####
HWTCObscooccur<-cooccur(HWTCSpsObsSeen)
#Plot localities and their habitat type####
#Create a custom color scale
library(RColorBrewer)
myColors <- brewer.pal(7,"Set1")
names(myColors) <- levels(DiversityCoordinates$Habitat)
colScale <- scale_colour_manual(name = "Habitat",values = myColors)

plotGeoHabitat <- ggplot(DiversityCoordinates,
            aes(x=LONGITUDE,y=LATITUDE,colour = Habitat)) + 
  geom_point()+
  labs(title = "eBird survey geographic locations on Grand Bahama 1988-2016")
plotGeoHabitat

plotGeoObsRich <- ggplot(DiversityCoordinates,
                  aes(x=LONGITUDE,y=LATITUDE,colour = ObserverRichness)) + 
  geom_point()+
  scale_colour_gradientn(colours=c("light blue", "black"))+   
  labs(title = "Observer Richness at eBird survey geographic locations on Grand Bahama 1988-2016")
plotGeoObsRich

plotGeoSpsRich <- ggplot(DiversityCoordinates,
                         aes(x=LONGITUDE,y=LATITUDE,colour = SpeciesRichness)) + 
  geom_point()+
  scale_colour_gradientn(colours=c("light blue", "red"))+   
  labs(title = "Species Richness at eBird survey geographic locations on Grand Bahama 1988-2016",
       x="Longitude", y = "Latitude")
plotGeoSpsRich



plotObsSpaceHabitat <- ggplot(DiversityCoordinates,
                  aes(x=obs.x.coord,y=obs.y.coord,colour = Habitat)) + 
  geom_point(aes(shape=factor(Habitat)))+
     labs(title = "eBird surveys by habitat in observer community space on Grand Bahama 1988-2016",
        x="x coordinate in jaccard dissimilarity space", y = "y coordinate in jaccard dissimilarity space" )
plotObsSpaceHabitat

plotSpsSpaceHabitat <- ggplot(DiversityCoordinates,
                       aes(x=jitter(sps.x.coord),
                           y=sps.y.coord,colour = Habitat)) + 
  geom_point(aes(shape=factor(Habitat)))+
  labs(title = "eBird surveys by habitat in observer community space on Grand Bahama 1988-2016",
       x="x coordinate in jaccard dissimilarity space", y = "y coordinate in jaccard dissimilarity space" )

plotSpsSpaceHabitat

summary(DiversityCoordinates$obs.x.coord)
plot(DiversityCoordinates$LONGITUDE, DiversityCoordinates$LATITUDE, col=DiversityCoordinates$Habitat, xlab="Longitude", ylab="Latitude", main="eBird Localities")
legend(-79,26.8, legend=c("Water", "Pine", "Wetland","Sand", "Urban", "Grass", "HWTC"),
       col=DiversityCoordinates$Habitat, lty=3, cex=0.8)
summary(DiversityCoordinates$Habitat)
#Plot corellogram of observer richness for distance in different spaces
plot(correlog(as.matrix(
  cbind(Diversityin3spaces[,3], Diversityin3spaces[,2])), 
  Diversityin3spaces$ObserverRichness), 
  main="Moran's I correlogram for observer richness over geographic distance")

plot(correlog(as.matrix(
  cbind(Diversityin3spaces[,3], Diversityin3spaces[,2])), 
  Diversityin3spaces$SpeciesRichness), 
  main="Moran's I correlogram for species richness over geographic distance")

plot(correlog(as.matrix(
  cbind(Diversityin3spaces[,3], Diversityin3spaces[,2])), 
  Diversityin3spaces$SpeciesShannonDiversity), 
  main="Moran's I correlogram for Shannon's Species Diversity over geographic distance")

plot(correlog(as.matrix(
  cbind(Diversityin3spaces[,3], Diversityin3spaces[,2])), 
  Diversityin3spaces$TotalMinutes), 
  main="Moran's I correlogram for observer effort in minutes over geographic distance")

plot(correlog(as.matrix(
  cbind(Diversityin3spaces[,3], Diversityin3spaces[,2])), 
  Diversityin3spaces$TotalSurveys), 
  main="Moran's I correlogram for observer effort in number of surveys over geographic distance")

plot(correlog(as.matrix(
  cbind(Diversityin3spaces[,17], Diversityin3spaces[,16])), 
  Diversityin3spaces$ObserverRichness), 
  main="Moran's I correlogram Observer Richness in Observer community space")

plot(correlog(as.matrix(
  cbind(Diversityin3spaces[,17], Diversityin3spaces[,16])), 
  Diversityin3spaces$SpeciesRichness), 
  main="Moran's I correlogram for species richness in Observer community space")

plot(correlog(as.matrix(
  cbind(Diversityin3spaces[,17], Diversityin3spaces[,16])), 
  Diversityin3spaces$SpeciesShannonDiversity), 
  main="Moran's I correlogram for Shannon's Species Diversity in Observer community space")

plot(correlog(as.matrix(
  cbind(Diversityin3spaces[,17], Diversityin3spaces[,16])), 
  Diversityin3spaces$TotalMinutes), 
  main="Moran's I correlogram for observer effort in minutes in Observer community space")

plot(correlog(as.matrix(
  cbind(Diversityin3spaces[,17], Diversityin3spaces[,16])), 
  Diversityin3spaces$TotalSurveys), 
  main="Moran's I correlogram for observer effort in number of surveys in Observer community space")
PineDiversityin3spaces<-Diversityin3spaces[Diversityin3spaces$Habitat=="Pine",]
