#libraries
library(raster)
library(rgdal)
#read dataset in a shapefile 
shape <- readOGR(dsn = ".", layer = "edaf_puntos_sii")
covar <- readRDS('covariates_1km_10PCA_6OGC_scaled.rds')

lim <- getData('GADM', country='MEX', level=1)

lim <- spTransform(lim, CRS(projection(covar)))

plot(lim)

covar <- crop(covar, drawExtent())

#covar <- as(covar, 'SpatialPixelsDataFrame')


proj4string(shape)<-crs('+proj=lcc +lat_1=17.5 +lat_2=29.5 +lat_0=12 +lon_0=-102 +x_0=2500000 +y_0=0 +ellps=WGS84 +units=m +no_defs')
#proj4string(shape)<-crs(projection(covar))

shape <- spTransform(shape, CRS(projection(covar)))

shape <- shape[covar,]

#Mexico limit 
#country limit from the global administrative areas project\
#https://gadm.org/
 

#visualize the points
#figure 1
x <- shape
#saturate for improving visualization of variability below 30% 
x$CO[x$CO>30] <- 30
#x <- x[is.na(x@data$CO)==FALSE,]
bubble(x, "CO",
        panel=function(...) {
          sp.polygons(lim, fill='white')
          sp:::panel.bubble(...)
        }) 
#and histograms
hist(shape$CO)
#log transformed
hist(log1p(shape$CO))


#generate a year only colum
year <- strsplit(as.character(shape$FECHA), '/')
y <- numeric()
for (i in 1:length(year)){
y[i] <- year[[i]][3]} 
shape$year <- as.numeric(y)
#select only those years after 1998 
shape@data <- shape@data[shape@data$year > 1998,]
shape@data <- na.omit(shape@data)
#bulk density funciton
estimateBD <- function(SOC, method="Saini1996"){
OM <- SOC * 1.724
if(method=="Saini1996"){BD <- 1.62 - 0.06 * OM}
if(method=="Drew1973"){BD <- 1 / (0.6268 + 0.0361 * OM)}
if(method=="Jeffrey1979"){BD <- 1.482 - 0.6786 * (log(OM))}
if(method=="Grigal1989"){BD <- 0.669 + 0.941 * exp(1)^(-0.06 * OM)}
if(method=="Adams1973"){BD <- 100 / (OM /0.244 + (100 - OM)/2.65)}
if(method=="Honeyset_Ratkowsky1989"){BD <- 1/(0.564 + 0.0556 * OM)}
return(BD)
}
shape@data$BLD <- estimateBD(shape@data$CO, method="Saini1996")
#coarse fragments data
shape$CRF <- as.numeric(shape$PEDREG)
shape$CRF[shape$CRF==1] <-  0
shape$CRF[shape$CRF==2] <-  20
shape$CRF[shape$CRF==3] <-  40
shape$CRF[shape$CRF==4] <-  60
shape$CRF[shape$CRF==5] <-  80

#algorthms for quantitative pedology
library(aqp)
sp4=shape@data
sp4$IDPROF <- paste0("IDPROF_", sp4$COORD_Y, "_", sp4$COORD_X)
#generate a soil profile collection object
depths(sp4) <- IDPROF  ~ LIM_SUP + LIM_INF
site(sp4) <- ~ COORD_X + COORD_Y
coordinates(sp4) <- ~ COORD_X + COORD_Y
proj4string(sp4)<-crs(projection(shape))
library(GSIF)
try(SOC <- mpspline(sp4, 'CO', d = t(c(0,100))))
try(BLD <- mpspline(sp4, 'BLD', d = t(c(0,100))))
try(CRF <- mpspline(sp4, 'CRF', d = t(c(0,100))))
dat <- data.frame(id = sp4@site$IDPROF,
X = sp4@sp@coords[,1],
Y = sp4@sp@coords[,2],
SOC = SOC$var.std[,1],
BLD = BLD$var.std[,1],
CRF = CRF$var.std[,1])
head(dat)


#agg <- slab(sp4, fm= ~ CO + BLD + CRF)
library(lattice)

xyplot(top ~ p.q50 | variable, data=agg, ylab='Profundidad (cm)',
             xlab='Valor medio de la variable dentro de los 25th y 75th percentiles',
             lower=agg$p.q25, upper=agg$p.q75, ylim=c(105,-2),
             panel=panel.depth_function,
             alpha=0.25, sync.colors=TRUE,
             par.settings=list(superpose.line=list(col=c('darkgray'), lwd=2)),
             prepanel=prepanel.depth_function,
             cf=agg$contributing_fraction, cf.col='black', cf.interval=5, 
             layout=c(3, 1), strip=strip.custom(bg=grey(0.8)),
             scales=list(x=list(tick.number=4, cex=1.5,alternating=3, relation='free'), y=list(cex=1.5))
             )

#Remove zero values

#Remove zero values

dat$SOC[dat$SOC==0] <- NA
dat <- na.omit(dat)

#Calculate SOC stocks
OCSKGM <- OCSKGM(ORCDRC = dat$SOC*10, BLD = dat$BLD*1000,
CRFVOL = dat$CRF, HSIZE = 100)
dat$OCSKGM <- OCSKGM
dat$meaERROR <- attr(OCSKGM,"measurementError")
dat <- dat[dat$OCSKGM>0,]
summary(dat)

#lis is a covariates files list names b


library(landmap)
library(rgdal)
library(geoR)
library(plotKML)
library(raster)
library(glmnet)
library(xgboost)
library(kernlab)
library(deepnet)
library(mlr)
library(mapview)

#covar <- covar[1:25]

#dat <- read.csv('soc_12_2020.xls - soc.csv')

coordinates(dat) <- ~ X + Y

proj4string(dat)<-crs(projection(covar))


m <- train.spLearner(dati["OCSKGM"], covariates=covar, lambda = 1, spc=FALSE,oblique.coords = FALSE, buffer.dist = FALSE, parallel=TRUE)

#saveRDS(m, file='trainSL_SOC_prof_effec.rds')

#m <- readRDS('trainSL_SOC_prof_effec.rds')

mp <- predict(m)

#saveRDS(m, file='trainSL_SOC_prof_effec_predicted.rds')

#model <- readRDS('trainSL_SOC_prof_effec.rds')

summary(m@spModel$learner.model$super.model$learner.model)

#prediction <- readRDS('trainSL_SOC_prof_effec_predicted.rds')

spplot(mp$pred[,c("response","q.lwr","q.upr")])

dim(dat)[1]

mean(dat$soc.t.ha)
sd(dat$soc.t.ha)

cellStats(expm1(stack(mp$pred[,c("response")])), mean)
cellStats(expm1(stack(mp$pred[,c("response")])), sd)


