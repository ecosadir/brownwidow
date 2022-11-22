# Overlap Vignette
# Brown and Carnaval 2019

library(devtools)
library(humboldt)
library(raster)
library(sp)
library(rgdal)
library(dplyr)
library(adehabitatHR)
library(spatstat)
library(tcltk)


#### LOAD RASTERS ####
# Set up is: Long , Lat , X1 , X2 , X... for both species
# Load rasters
rasternames <- list.files("/path", 
                          pattern="*.asc", 
                          full.names=TRUE)
rasterapply <- lapply(rasternames, raster)

for (i in 1:length(rasternames)) {
  assign(paste("Bio", i, sep = ""), 
         raster(rasternames[i]))
  }

bb<-paste("Bio",1:19, sep="") # create character string for each raster object name
foo <- lapply(bb, get) # creates a list of the objects from character strings
envs.all <- stack(foo)
proj <- CRS(as.character("+proj=longlat +datum=WGS84")) #add projection
proj4string(envs.all)<- proj


#### OCCURRENCES ####
# occurence data
all <- read.csv("/path")

colnames(all) <- c("sp", "x", "y")

# Filter by species
occs.geo <- all %>% 
  filter(sp == "Latrodectus geometricus")


occs.hesp <- all %>% 
  filter(sp == "Latrodectus hesperus")


occs.var <- all %>% 
  filter(sp == "Latrodectus variolus")


occs.mac <- all %>% 
  filter(sp == "Latrodectus mactans")





#### ENVIRONMENTAL SPACE ####
#convert one of the rasters to a point dataframe to sample.  Use any raster input.
env.points<-rasterToPoints(Bio3, fun=NULL, spatial=FALSE)

# Subset only the x and y data
env.sampling.res<- env.points[,1:2]

# Extract values to points from rasters
RAST_VAL<-data.frame(extract(envs.all, env.sampling.res))

# merge sampled data to input
Env1<-cbind(env.sampling.res,RAST_VAL)
Env2<-Env1

# remove NAs and make sure all variables are imported as numbers
Env1<-humboldt.scrub.env(Env1)
Env2<-humboldt.scrub.env(Env2)

#save the file as '.csv' for future analyses 
write.csv(Env1, file = "/path")
write.csv(Env2, file = "/path")



#### GEOMETRIUS HESPERUS ####


# Select top variables
reduc.vars.geohesp <- humboldt.top.env(
  Env1,
  Env2,
  occs.geo,
  occs.hesp,
  rarefy.dist = 0,
  rarefy.units = "km",
  env.reso=0.0416669,
  learning.rt1 = 0.01,
  learning.rt2 = 0.01,
  e.var=c(3:21),
  pa.ratio = 5,
  steps1 = 50,
  steps2 = 50,
  method = "contrib",
  contrib.greater = 10) # keep variables that contribute > 5% 


write.csv(reduc.vars.geohesp, file = "variables_geohesp.csv")

## Adjust the number of variables input for e.vars after reduction to only important variables
num.var.e.geohesp<-ncol(reduc.vars.geohesp$env1)


# # # NICHE OVERLAP # # # 


## run it first with full environmental for backgroud tests and equivalence statistic 
## (total equivalence or divergence in current distributions)
full.geohesp <- humboldt.doitall(inname="doitall_full_extent_geohesp", 
                                 env1=reduc.vars.geohesp$env1, 
                                 env2=reduc.vars.geohesp$env2, 
                                 sp1=occs.geo, 
                                 sp2=occs.hesp, 
                                 rarefy.dist=0, 
                                 rarefy.units="km", 
                                 env.reso=0.416669, 
                                 reduce.env=0, 
                                 reductype="PCA", 
                                 non.analogous.environments="YES", 
                                 correct.env=T, 
                                 env.trim=T,  
                                 env.trim.type="MCP", 
                                 trim.buffer.sp1=500, 
                                 trim.buffer.sp2=500, 
                                 pcx=1, 
                                 pcy=2, 
                                 col.env= e.var, 
                                 e.var=c(3:num.var.e.geohesp), 
                                 R=100, 
                                 kern.smooth=1, 
                                 e.reps=500, 
                                 b.reps=500, 
                                 nae="YES",
                                 thresh.espace.z=0.0001, 
                                 p.overlap=T, 
                                 p.boxplot=T, 
                                 p.scatter=T, 
                                 run.silent=F, 
                                 ncores="all")



# # # NICHE DIVERGENCE # # # 

## run it a second time with a trimmed, shared-espace. Here the equivalence statistic tests for niche evolution or niche 
## divergence. For comparing results, change only the following model parameters: reduce.env, 
## non.analogous.environmental, env.trim
shared_ae.geohesp <- humboldt.doitall(inname="shared_espace_ae_geohesp", 
                                      env1=reduc.vars.geohesp$env1, 
                                      env2=reduc.vars.geohesp$env2, 
                                      sp1=occs.geo, 
                                      sp2=occs.hesp, 
                                      rarefy.dist=0, 
                                      rarefy.units="km", 
                                      env.reso=0.416669, 
                                      reduce.env=0, 
                                      reductype="PCA", 
                                      non.analogous.environments="NO", 
                                      nae.window=2,
                                      correct.env=T, 
                                      env.trim=T, 
                                      env.trim.type="MCP", 
                                      trim.buffer.sp1=500, 
                                      trim.buffer.sp2=500, 
                                      pcx=1,
                                      pcy=2, 
                                      col.env=e.var, 
                                      e.var=c(3:num.var.e.geohesp), 
                                      R=100, 
                                      kern.smooth=1, 
                                      e.reps=500, 
                                      b.reps=500, 
                                      nae="NO",
                                      thresh.espace.z=0.0001, 
                                      p.overlap=T, 
                                      p.boxplot=T, 
                                      p.scatter=T,
                                      run.silent=F, 
                                      ncores="all")





# # # # EXTRA PRACTICE JUNK # # # # 
# HUMBOLDT.G2E
## Converts geographic space to environmental space
g2e.geohesp <- humboldt.g2e(
  Env1,
  Env2,
  occs.geo,
  occs.hesp,
  reduce.env = 0,
  reductype = "PCA",
  non.analogous.environments = "NO",
  nae.window = 5,
  env.trim = T,
  e.var = c(3:8),
  col.env = e.var,
  env.trim.type = "MCP",
  trim.buffer.sp1 = 750,
  trim.buffer.sp2 = 750,
  pcx = 1,
  pcy = 2,
  rarefy.dist = 0,
  rarefy.units = "km",
  env.reso = 0.0416669,
  kern.smooth = 1,
  R = 100,
  run.silent = F)



# CREATE GRID E-SPACE
## store espace scores for sp1 and environments 1,2 and both environments combined output from humboldt.g2e
scores.env1 <- g2e.geohesp$scores.env1[1:2]
scores.env2 <- g2e.geohesp$scores.env2[1:2]
scores.env12 <- rbind(g2e.geohesp$scores.env1[1:2],g2e.geohesp$scores.env2[1:2])
scores.sp1 <- g2e.geohesp$scores.sp1[1:2]
scores.sp2 <- g2e.geohesp$scores.sp2[1:2]



Geospace <- humboldt.grid.espace(scores.env12,scores.env1,scores.sp1,kern.smooth=1,R=100)

Hespspace <- humboldt.grid.espace(scores.env12,scores.env1,scores.sp1,kern.smooth=1,R=100)



# BACKGROUND TEST
## Runs a modified niche background test
### If the observed values of the niche similarity measures obtained from the two original populations are significantly 
### higher than expected from this null distribution, then the null hypothesis that similarity between species is solely 
### due to differences in habitat availability is rejected. Or in other words, a significant value suggests that the two 
### species are more divergent than expected solely on the habitat availability (pending a significant equivalency test). 
### A non-significant test suggest that most occupied environments are similar among environments.


### This test measures the power of the equivalency test to detect differences based on the available e-space. 
### If both the equivalency test and background test are non-significant, this means the two species occupied 
### environmental spaces are not significantly different and resulting niche 'equivalence' is likely a result of 
### the limited environmental space present in habitat(s).

bg.test.geohesp <- humboldt.background.stat(g2e.geohesp, rep = 100, sim.dir = 1, env.reso=0.0416669,
                                             kern.smooth = 1, R = 100, correct.env = F,
                                             thresh.espace.z = 1e-04, force.equal.sample = F,
                                             run.silent.bak = F, ncores = "All")


# BACKGROUND STAT
##Runs a modified niche background statistic(see Warren et al 2008) based on two species' occurrence density grids. 

### A significant value suggests that the two species are more divergent than expected solely on the habitat 
### availability (pending a significant equivalence statistic). 
### A non-significant test suggest that most occupied environments are similar among environments. 



### This test measures the power of the equivalence statistic to detect differences based on the available e-space. 
### If both the equivalence statistic and background statistic are non-significant, this means the two species occupied 
### environmental spaces are not significantly different and resulting niche 'equivalence' is likely a result of the 
### limited environmental space present in habitat(s). Basically in these situations, there is limited power for the 
### equivalence statistic to actually detect and significant differences among taxa, even if they existed. 
### However, conversely it also doesn't provide any evidence that they in fact on not equivalent- simply there is 
### little power to detect it in input environmental data.


### If both background statistic are non-significant and equivalence statistic are non-significant, try to increase 
### the spatial extent of input climate data. If the equivalence statistic is significant and this is not, this means 
### the two species occupied environmental space is significantly different despite the existence of largely similar habitats. 
### This is strong evidence of niche divergence.

bg.stat.geohesp <- humboldt.background.stat(
  g2e.geohesp,
  rep = 100,
  sim.dir = 1,
  env.reso=0.0416669,
  kern.smooth = 1,
  R = 100,
  correct.env = T,
  thresh.espace.z = 0.001,
  force.equal.sample = F,
  run.silent.bak = F,
  ncores = "All")


# EQUIVALENCY STAT

equiv.stat.geohesp <- humboldt.equivalence.stat(Geospace, 
                                                Hespspace,
                                                rep = 100,
                                                correct.env = T,
                                                kern.smooth = 1,
                                                nae = "YES",
                                                thresh.espace.z = 0.001,
                                                run.silent.equ = F,
                                                ncores = "All")













#### GEOMETRIUS MACTANS ####

# Select top variables
reduc.vars.geomac <- humboldt.top.env(
  Env1,
  Env2,
  occs.geo,
  occs.mac,
  rarefy.dist = 0,
  rarefy.units = "km",
  env.reso=0.0416669,
  learning.rt1 = 0.01,
  learning.rt2 = 0.01,
  e.var=c(3:21),
  pa.ratio = 5,
  steps1 = 50,
  steps2 = 50,
  method = "contrib",
  contrib.greater = 10) # keep variables that contribute > 5% 


write.csv(reduc.vars.geomac, file = "variables_geomac.csv")



## Adjust the number of variables input for e.vars after reduction to only important variables
num.var.e.geomac<-ncol(reduc.vars.geomac$env1)


# # # NICHE OVERLAP # # # 


## run it first with full environmental for backgroud tests and equivalence statistic 
## (total equivalence or divergence in current distributions)
full.geomac <- humboldt.doitall(inname="doitall_full_extent_geomac", 
                                 env1=reduc.vars.geomac$env1, 
                                 env2=reduc.vars.geomac$env2, 
                                 sp1=occs.geo, 
                                 sp2=occs.mac, 
                                 rarefy.dist=0, 
                                 rarefy.units="km", 
                                 env.reso=0.416669, 
                                 reduce.env=0, 
                                 reductype="PCA", 
                                 non.analogous.environments="YES", 
                                 correct.env=T, 
                                 env.trim=T,  
                                 env.trim.type="MCP", 
                                 trim.buffer.sp1=500, 
                                 trim.buffer.sp2=500, 
                                 pcx=1, 
                                 pcy=2, 
                                 col.env= e.var, 
                                 e.var=c(3:num.var.e.geomac), 
                                 R=100, 
                                 kern.smooth=1, 
                                 e.reps=500, 
                                 b.reps=500, 
                                 nae="YES",
                                 thresh.espace.z=0.0001, 
                                 p.overlap=T, 
                                 p.boxplot=T, 
                                 p.scatter=T, 
                                 run.silent=F, 
                                 ncores="all")



# # # NICHE DIVERGENCE # # # 

## run it a second time with a trimmed, shared-espace. Here the equivalence statistic tests for niche evolution or niche 
## divergence. For comparing results, change only the following model parameters: reduce.env, 
## non.analogous.environmental, env.trim
shared_ae.geomac <- humboldt.doitall(inname="shared_espace_ae_geomac", 
                                      env1=reduc.vars.geomac$env1, 
                                      env2=reduc.vars.geomac$env2, 
                                      sp1=occs.geo, 
                                      sp2=occs.mac, 
                                      rarefy.dist=0, 
                                      rarefy.units="km", 
                                      env.reso=0.416669, 
                                      reduce.env=0, 
                                      reductype="PCA", 
                                      non.analogous.environments="NO", 
                                      nae.window=2,
                                      correct.env=T, 
                                      env.trim=T, 
                                      env.trim.type="MCP", 
                                      trim.buffer.sp1=500, 
                                      trim.buffer.sp2=500, 
                                      pcx=1,
                                      pcy=2, 
                                      col.env=e.var, 
                                      e.var=c(3:num.var.e.geomac), 
                                      R=100, 
                                      kern.smooth=1, 
                                      e.reps=500, 
                                      b.reps=500, 
                                      nae="NO",
                                      thresh.espace.z=0.0001, 
                                      p.overlap=T, 
                                      p.boxplot=T, 
                                      p.scatter=T,
                                      run.silent=F, 
                                      ncores="all")















#### GEOMETRIUS VARIOLUS ####

# Select top variables
reduc.vars.geovar <- humboldt.top.env(
  Env1,
  Env2,
  occs.geo,
  occs.var,
  rarefy.dist = 0,
  rarefy.units = "km",
  env.reso=0.0416669,
  learning.rt1 = 0.01,
  learning.rt2 = 0.01,
  e.var=c(3:21),
  pa.ratio = 5,
  steps1 = 50,
  steps2 = 50,
  method = "contrib",
  contrib.greater = 10) # keep variables that contribute > 5% 


write.csv(reduc.vars.geovar, file = "variables_geovar.csv")

## Adjust the number of variables input for e.vars after reduction to only important variables
num.var.e.geovar<-ncol(reduc.vars.geovar$env1)


# # # NICHE OVERLAP # # # 


## run it first with full environmental for backgroud tests and equivalence statistic 
## (total equivalence or divergence in current distributions)
full.geovar <- humboldt.doitall(inname="doitall_full_extent_geovar", 
                                env1=reduc.vars.geovar$env1, 
                                env2=reduc.vars.geovar$env2, 
                                sp1=occs.geo, 
                                sp2=occs.var, 
                                rarefy.dist=0, 
                                rarefy.units="km", 
                                env.reso=0.416669, 
                                reduce.env=0, 
                                reductype="PCA", 
                                non.analogous.environments="YES", 
                                correct.env=T, 
                                env.trim=T,  
                                env.trim.type="MCP", 
                                trim.buffer.sp1=500, 
                                trim.buffer.sp2=500, 
                                pcx=1, 
                                pcy=2, 
                                col.env= e.var, 
                                e.var=c(3:num.var.e.geovar), 
                                R=100, 
                                kern.smooth=1, 
                                e.reps=500, 
                                b.reps=500, 
                                nae="YES",
                                thresh.espace.z=0.0001, 
                                p.overlap=T, 
                                p.boxplot=T, 
                                p.scatter=T, 
                                run.silent=F, 
                                ncores="all")



# # # NICHE DIVERGENCE # # # 

## run it a second time with a trimmed, shared-espace. Here the equivalence statistic tests for niche evolution or niche 
## divergence. For comparing results, change only the following model parameters: reduce.env, 
## non.analogous.environmental, env.trim
shared_ae.geovar <- humboldt.doitall(inname="shared_espace_ae_geovar", 
                                     env1=reduc.vars.geovar$env1, 
                                     env2=reduc.vars.geovar$env2, 
                                     sp1=occs.geo, 
                                     sp2=occs.var, 
                                     rarefy.dist=0, 
                                     rarefy.units="km", 
                                     env.reso=0.416669, 
                                     reduce.env=0, 
                                     reductype="PCA", 
                                     non.analogous.environments="NO", 
                                     nae.window=2,
                                     correct.env=T, 
                                     env.trim=T, 
                                     env.trim.type="MCP", 
                                     trim.buffer.sp1=500, 
                                     trim.buffer.sp2=500, 
                                     pcx=1,
                                     pcy=2, 
                                     col.env=e.var, 
                                     e.var=c(3:num.var.e.geovar), 
                                     R=100, 
                                     kern.smooth=1, 
                                     e.reps=500, 
                                     b.reps=500, 
                                     nae="NO",
                                     thresh.espace.z=0.0001, 
                                     p.overlap=T, 
                                     p.boxplot=T, 
                                     p.scatter=T,
                                     run.silent=F, 
                                     ncores="all")















#### HESPERUS VARIOLUS ####

# Select top variables
reduc.vars.hespvar <- humboldt.top.env(
  Env1,
  Env2,
  occs.hesp,
  occs.var,
  rarefy.dist = 0,
  rarefy.units = "km",
  env.reso=0.0416669,
  learning.rt1 = 0.01,
  learning.rt2 = 0.01,
  e.var=c(3:21),
  pa.ratio = 5,
  steps1 = 50,
  steps2 = 50,
  method = "contrib",
  contrib.greater = 10) # keep variables that contribute > 5% 

write.csv(reduc.vars.hespvar, file = "variables_hespvar.csv")


## Adjust the number of variables input for e.vars after reduction to only important variables
num.var.e.hespvar <- ncol(reduc.vars.hespvar$env1)


# # # NICHE OVERLAP # # # 


## run it first with full environmental for backgroud tests and equivalence statistic 
## (total equivalence or divergence in current distributions)
full.hespvar <- humboldt.doitall(inname="doitall_full_extent_hespvar", 
                                 env1=reduc.vars.hespvar$env1, 
                                 env2=reduc.vars.hespvar$env2, 
                                 sp1=occs.hesp, 
                                 sp2=occs.var, 
                                 rarefy.dist=0, 
                                 rarefy.units="km", 
                                 env.reso=0.416669, 
                                 reduce.env=0, 
                                 reductype="PCA", 
                                 non.analogous.environments="YES", 
                                 correct.env=T, 
                                 env.trim=T,  
                                 env.trim.type="MCP", 
                                 trim.buffer.sp1=500, 
                                 trim.buffer.sp2=500, 
                                 pcx=1, 
                                 pcy=2, 
                                 col.env= e.var, 
                                 e.var=c(3:num.var.e.hespvar), 
                                 R=100, 
                                 kern.smooth=1, 
                                 e.reps=500, 
                                 b.reps=500, 
                                 nae="YES",
                                 thresh.espace.z=0.0001, 
                                 p.overlap=T, 
                                 p.boxplot=T, 
                                 p.scatter=T, 
                                 run.silent=F, 
                                 ncores="all")



# # # NICHE DIVERGENCE # # # 

## run it a second time with a trimmed, shared-espace. Here the equivalence statistic tests for niche evolution or niche 
## divergence. For comparing results, change only the following model parameters: reduce.env, 
## non.analogous.environmental, env.trim
shared_ae.hespvar <- humboldt.doitall(inname="shared_espace_ae_hespvar", 
                                      env1=reduc.vars.hespvar$env1, 
                                      env2=reduc.vars.hespvar$env2, 
                                      sp1=occs.hesp, 
                                      sp2=occs.var, 
                                      rarefy.dist=0, 
                                      rarefy.units="km", 
                                      env.reso=0.416669, 
                                      reduce.env=0, 
                                      reductype="PCA", 
                                      non.analogous.environments="NO", 
                                      nae.window=2,
                                      correct.env=T, 
                                      env.trim=T, 
                                      env.trim.type="MCP", 
                                      trim.buffer.sp1=500, 
                                      trim.buffer.sp2=500, 
                                      pcx=1,
                                      pcy=2, 
                                      col.env=e.var, 
                                      e.var=c(3:num.var.e.hespvar), 
                                      R=100, 
                                      kern.smooth=1, 
                                      e.reps=500, 
                                      b.reps=500, 
                                      nae="NO",
                                      thresh.espace.z=0.0001, 
                                      p.overlap=T, 
                                      p.boxplot=T, 
                                      p.scatter=T,
                                      run.silent=F, 
                                      ncores="all")











#### HESPERUS MACTANS ####

## Select top variables
reduc.vars.hespmac <- humboldt.top.env(
  Env1,
  Env2,
  occs.hesp,
  occs.mac,
  rarefy.dist = 0,
  rarefy.units = "km",
  env.reso=0.0416669,
  learning.rt1 = 0.01,
  learning.rt2 = 0.01,
  e.var=c(3:21),
  pa.ratio = 5,
  steps1 = 50,
  steps2 = 50,
  method = "contrib",
  contrib.greater = 10) # keep variables that contribute > 5% 

write.csv(reduc.vars.hespmac, file = "variables_hespmac.csv")

## Adjust the number of variables input for e.vars after reduction to only important variables
num.var.e<-ncol(reduc.vars.hespmac$env1)


# # # NICHE OVERLAP # # # 


## run it first with full environmental for backgroud tests and equivalence statistic 
## (total equivalence or divergence in current distributions)
full.hespmac <- humboldt.doitall(inname="doitall_full_extent_hespmac", 
                       env1=reduc.vars.hespmac$env1, 
                       env2=reduc.vars.hespmac$env2, 
                       sp1=occs.hesp, 
                       sp2=occs.mac, 
                       rarefy.dist=0, 
                       rarefy.units="km", 
                       env.reso=0.416669, 
                       reduce.env=0, 
                       reductype="PCA", 
                       non.analogous.environments="YES", 
                       correct.env=T, 
                       env.trim=T,  
                       env.trim.type="MCP", 
                       trim.buffer.sp1=500, 
                       trim.buffer.sp2=500, 
                       pcx=1, 
                       pcy=2, 
                       col.env= e.var, 
                       e.var=c(3:num.var.e), 
                       R=100, 
                       kern.smooth=1, 
                       e.reps=500, 
                       b.reps=500, 
                       nae="YES",
                       thresh.espace.z=0.0001, 
                       p.overlap=T, 
                       p.boxplot=T, 
                       p.scatter=T, 
                       run.silent=F, 
                       ncores="all")



# # # NICHE DIVERGENCE # # # 

## run it a second time with a trimmed, shared-espace. Here the equivalence statistic tests for niche evolution or niche 
## divergence. For comparing results, change only the following model parameters: reduce.env, 
## non.analogous.environmental, env.trim
shared_ae.hespmac <- humboldt.doitall(inname="shared_espace_ae_hespmac", 
                            env1=reduc.vars.hespmac$env1, 
                            env2=reduc.vars.hespmac$env2, 
                            sp1=occs.hesp, 
                            sp2=occs.mac, 
                            rarefy.dist=0, 
                            rarefy.units="km", 
                            env.reso=0.416669, 
                            reduce.env=0, 
                            reductype="PCA", 
                            non.analogous.environments="NO", 
                            nae.window=2,
                            correct.env=T, 
                            env.trim=T, 
                            env.trim.type="MCP", 
                            trim.buffer.sp1=500, 
                            trim.buffer.sp2=500, 
                            pcx=1,
                            pcy=2, 
                            col.env=e.var, 
                            e.var=c(3:num.var.e), 
                            R=100, 
                            kern.smooth=1, 
                            e.reps=500, 
                            b.reps=500, 
                            nae="NO",
                            thresh.espace.z=0.0001, 
                            p.overlap=T, 
                            p.boxplot=T, 
                            p.scatter=T,
                            run.silent=F, 
                            ncores="all")


#### MACTANS VARIOLUS ####


# Select top variables
reduc.vars.macvar <- humboldt.top.env(
  Env1,
  Env2,
  occs.mac,
  occs.var,
  rarefy.dist = 0,
  rarefy.units = "km",
  env.reso=0.0416669,
  learning.rt1 = 0.01,
  learning.rt2 = 0.01,
  e.var=c(3:21),
  pa.ratio = 5,
  steps1 = 50,
  steps2 = 50,
  method = "contrib",
  contrib.greater = 10) # keep variables that contribute > 5% 

write.csv(reduc.vars.macvar, file = "variables_macvar.csv")


## Adjust the number of variables input for e.vars after reduction to only important variables
num.var.e.macvar <- ncol(reduc.vars.macvar$env1)


# # # NICHE OVERLAP # # # 


## run it first with full environmental for backgroud tests and equivalence statistic 
## (total equivalence or divergence in current distributions)
full.macvar <- humboldt.doitall(inname="doitall_full_extent_macvar", 
                                 env1=reduc.vars.macvar$env1, 
                                 env2=reduc.vars.macvar$env2, 
                                 sp1=occs.mac, 
                                 sp2=occs.var, 
                                 rarefy.dist=0, 
                                 rarefy.units="km", 
                                 env.reso=0.416669, 
                                 reduce.env=0, 
                                 reductype="PCA", 
                                 non.analogous.environments="YES", 
                                 correct.env=T, 
                                 env.trim=T,  
                                 env.trim.type="MCP", 
                                 trim.buffer.sp1=500, 
                                 trim.buffer.sp2=500, 
                                 pcx=1, 
                                 pcy=2, 
                                 col.env= e.var, 
                                 e.var=c(3:num.var.e.macvar), 
                                 R=100, 
                                 kern.smooth=1, 
                                 e.reps=500, 
                                 b.reps=500, 
                                 nae="YES",
                                 thresh.espace.z=0.0001, 
                                 p.overlap=T, 
                                 p.boxplot=T, 
                                 p.scatter=T, 
                                 run.silent=F, 
                                 ncores="all")



# # # NICHE DIVERGENCE # # # 

## run it a second time with a trimmed, shared-espace. Here the equivalence statistic tests for niche evolution or niche 
## divergence. For comparing results, change only the following model parameters: reduce.env, 
## non.analogous.environmental, env.trim
shared_ae.macvar <- humboldt.doitall(inname="shared_espace_ae_macvar", 
                                      env1=reduc.vars.macvar$env1, 
                                      env2=reduc.vars.macvar$env2, 
                                      sp1=occs.mac, 
                                      sp2=occs.var, 
                                      rarefy.dist=0, 
                                      rarefy.units="km", 
                                      env.reso=0.416669, 
                                      reduce.env=0, 
                                      reductype="PCA", 
                                      non.analogous.environments="NO", 
                                      nae.window=2,
                                      correct.env=T, 
                                      env.trim=T, 
                                      env.trim.type="MCP", 
                                      trim.buffer.sp1=500, 
                                      trim.buffer.sp2=500, 
                                      pcx=1,
                                      pcy=2, 
                                      col.env=e.var, 
                                      e.var=c(3:num.var.e.macvar), 
                                      R=100, 
                                      kern.smooth=1, 
                                      e.reps=500, 
                                      b.reps=500, 
                                      nae="NO",
                                      thresh.espace.z=0.0001, 
                                      p.overlap=T, 
                                      p.boxplot=T, 
                                      p.scatter=T,
                                      run.silent=F, 
                                      ncores="all")

