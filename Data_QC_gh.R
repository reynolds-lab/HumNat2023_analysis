#### Geographic Analysis QC ####
#Austin W Reynolds
#Mark N Grote
#2021.03.24

### Table of Contents ###
#1. Dependencies
#2. User defined variables
#3. Migration Distances
#4. Ancestry
#5. Create Longform Hierarchical Dataset
#6. Censored Migration Data
#7. Derived variables for Modeling
#8. Create Dataset for Migration Model
#9. Create Dataset for Locality Model
#10. Write datasets to file

#this analysis workbook requires two input files:
#"latlondata" which is a csv file with three columns containing place names, latitude and longitude.
#see example file "placedata.csv" for format
#"demographicdata" which is another csv containing demographic information for probands and their families
#see example file "pedigreedata.csv" for format



### 1. Read in dependencies ###
library(geosphere)
library(maptools)
library(data.table)
library(tidyverse)


### 2. User defined variables ###
latlonpath <- "~/Dropbox/Cederberg_project/Datasets/combined_placedata.csv" #path to Latitude and Longitude data for birthplaces of probands and parents/grandparents of probands
demographicpath <- "~/Dropbox/Cederberg_project/Datasets/combined_pedigreedata_newancestry.csv" #path to raw demographic data for probands and their families
sibs <- c() #vector of Subject_IDs for one member of each pair of siblings in the demographic data. These were removed due to redundancy in parent and grandparent data.
maxdistpath <- "~/Dropbox/Cederberg_project/Scripts/maxdist2line.r" #path to function to calculate max distance
outputpath <- "~/Dropbox/Cederberg_project/Datasets/" #path to the directory where final datasets will be written

#source the long function that calculates max distance across a polygon
source(maxdistpath)

#input boundary coordinates for each polygon here
#these are used when an individual's birthplace is a larger region (province, state, country) instead of a point location (town, city, village)
# for our analysis, these were regions that were named as birthplaces by a sizeable number of individuals in the datasets.
regional.polygons <- list(
  Bosmanland = matrix(c(22.389911,-30.483715, 
                        19.950946,-30.937101, 
                        19.028094,-30.398468, 
                        18.632586,-29.828259, 
                        18.742450,-29.168475, 
                        19.885028,-29.072499, 
                        20.643084,-28.889897, 
                        21.346209,-28.793662, 
                        21.730731,-28.976433, 
                        22.389911,-30.483715),ncol=2,byrow = TRUE),
  CDB = matrix(c(18.679861,-31.793486,
                 18.602957,-32.170894,
                 18.773245,-32.500428,
                 18.872122,-32.542114,
                 18.970999,-32.639306,
                 19.031423,-32.782588,
                 19.042410,-32.907194,
                 19.031423,-33.031624,
                 19.047903,-33.192663,
                 19.069876,-33.312100,
                 19.229177,-33.435957,
                 19.399465,-33.440541,
                 19.712576,-33.321280,
                 19.641165,-32.778374,
                 19.509329,-32.491567,
                 19.053396,-31.750534,
                 18.679861,-31.793486),ncol=2,byrow = TRUE),
  NMB = matrix(c(20.045721,-22.017536, 
                 20.018674,-28.437526, 
                 19.864866,-28.447187, 
                 19.820920,-28.505129, 
                 19.733030,-28.495475, 
                 19.678098,-28.543740, 
                 19.612180,-28.514784, 
                 19.557248,-28.563041, 
                 19.491331,-28.630563, 
                 19.458372,-28.698043, 
                 19.326536,-28.755847, 
                 19.249631,-28.765478, 
                 19.271604,-28.880981, 
                 19.238645,-28.919454, 
                 19.150754,-28.967524, 
                 19.062864,-28.957912, 
                 18.953000,-28.871361, 
                 18.755247,-28.852118, 
                 18.557493,-28.871361, 
                 18.414670,-28.880981, 
                 18.304807,-28.861740, 
                 18.194944,-28.861740, 
                 18.063108,-28.852118, 
                 17.920286,-28.746215, 
                 17.733518,-28.746215, 
                 17.557737,-28.669128, 
                 17.491819,-28.688405, 
                 17.447873,-28.601630, 
                 17.403928,-28.466504, 
                 17.436887,-28.360213, 
                 17.370969,-28.234458, 
                 17.294065,-28.186051, 
                 17.261106,-28.215098, 
                 17.228147,-28.176367, 
                 17.206174,-28.098863, 
                 17.151243,-28.040699, 
                 17.085325,-28.011605, 
                 16.909543,-28.069785, 
                 16.843625,-28.108554, 
                 16.799680,-28.234458, 
                 16.733762,-28.292517, 
                 16.755735,-28.350544, 
                 16.744748,-28.456846, 
                 16.667844,-28.456846, 
                 16.579954,-28.485819, 
                 16.437131,-28.572689, 
                 16.415159,-28.601630, 
                 15.668088,-27.875729, 
                 15.338498,-27.447583, 
                 15.140745,-26.763009, 
                 14.942991,-26.113752, 
                 14.855100,-25.361634, 
                 14.811155,-24.884180, 
                 14.591428,-24.564846, 
                 14.459592,-24.024124, 
                 14.415647,-22.286795, 
                 19.952756,-22.083335),ncol=2,byrow = TRUE),
  Richtersveld = matrix(c(17.008642,-29.555820,
                          17.843603,-29.508025,
                          17.761205,-29.268712,
                          17.849096,-29.019231,
                          17.634862,-28.773960,
                          17.568944,-28.692075,
                          17.426122,-28.730617,
                          17.426122,-28.576363,
                          17.332738,-28.484667,
                          17.409643,-28.388059,
                          17.327245,-28.233304,
                          17.244848,-28.242982,
                          17.189916,-28.218784,
                          17.145971,-28.097708,
                          17.085546,-28.054087,
                          17.047094,-28.049239,
                          16.964696,-28.068629,
                          16.882299,-28.088016,
                          16.871313,-28.184896,
                          16.843847,-28.165527,
                          16.816381,-28.218784,
                          16.810888,-28.281689,
                          16.777929,-28.272014,
                          16.722997,-28.281689,
                          16.750463,-28.478977,
                          16.701025,-28.483805,
                          16.453832,-28.618909,
                          17.008642,-29.555820),ncol=2,byrow = TRUE)
)
#input region names here
#used in for loop for assigning max/min migration distances to individuals in these regions
regional.names<-c("Bosmanland","CDB","NMB","Richtersveld")




### 3. Migration ###

#Define function to assign distances
distassn <- function(df,dm){
  #function to assign calculated distances from a matrix (dm)
  #to place names in a dataframe (df)
  for (i in 1:nrow(df)){
    df$distance[i] <- dm[df$V2[i],df$V3[i]]
  }
  return(df)
}

#read in latlon dataset
latlondata<-read.csv(latlonpath,stringsAsFactors = F)

#create matrix of distances between all included locations
x<-latlondata[,c(2,3)] #remove placenames for convienence
distancematrix<-distm(x,x)/1000 # calculate distances in km between all locations
#assign colnames and rownames of distancematrix to the placenames from latlondata
colnames(distancematrix)<-latlondata$birthplace
rownames(distancematrix)<-latlondata$birthplace

#read in demographic dataset
demographicdata<-read.csv(demographicpath,stringsAsFactors = F,colClasses = c("character"))

#remove siblings to avoid double counting birthplaces 
demographicdata <- demographicdata[which(!demographicdata$Subject_ID %in% sibs),]

#Downstream data wrangling cannot accept duplicate subject IDs
#True duplicates will need to be detected and removed

#create dataframes for parent/offspring distances
#proband-mother distance
sm<-as.data.frame(cbind(demographicdata$Subject_ID,demographicdata$birthplace,demographicdata$m_birthplace,demographicdata$birthyear,demographicdata$m_birthyear,demographicdata$eg,demographicdata$m_eg),stringsAsFactors = F)
#sm$V2[which(!sm$V2%in%rownames(distancematrix))] #these are checks to see if any names are mispelled. If everything is good, the results should all be NA
sm$V2[which(!sm$V2%in%rownames(distancematrix))]<-"Unknown" #assign unknown birthplaces (NA) to string "Unknown" for convenience
#sm$V3[which(!sm$V3%in%rownames(distancematrix))]
sm$V3[which(!sm$V3%in%rownames(distancematrix))]<-"Unknown"

#proband-father distance
sf<-as.data.frame(cbind(demographicdata$Subject_ID,demographicdata$birthplace,demographicdata$f_birthplace,demographicdata$birthyear,demographicdata$f_birthyear,demographicdata$eg,demographicdata$f_eg),stringsAsFactors = F)
#sf$V2[which(!sf$V2%in%rownames(distancematrix))]
sf$V2[which(!sf$V2%in%rownames(distancematrix))]<-"Unknown"
#sf$V3[which(!sf$V3%in%rownames(distancematrix))]
sf$V3[which(!sf$V3%in%rownames(distancematrix))]<-"Unknown"

#mother-maternal grandmother distance
mmgm<-as.data.frame(cbind(demographicdata$Subject_ID,demographicdata$m_birthplace,demographicdata$mgm_birthplace,demographicdata$m_birthyear,demographicdata$mgm_birthyear,demographicdata$m_eg,demographicdata$mgm_eg),stringsAsFactors = F)
#mmgm$V2[which(!mmgm$V2%in%rownames(distancematrix))]
mmgm$V2[which(!mmgm$V2%in%rownames(distancematrix))]<-"Unknown"
#mmgm$V3[which(!mmgm$V3%in%rownames(distancematrix))]
mmgm$V3[which(!mmgm$V3%in%rownames(distancematrix))]<-"Unknown"

#mother-maternal grandfather distance
mmgf<-as.data.frame(cbind(demographicdata$Subject_ID,demographicdata$m_birthplace,demographicdata$mgf_birthplace,demographicdata$m_birthyear,demographicdata$mgf_birthyear,demographicdata$m_eg,demographicdata$mgf_eg),stringsAsFactors = F)
#mmgf$V2[which(!mmgf$V2%in%rownames(distancematrix))]
mmgf$V2[which(!mmgf$V2%in%rownames(distancematrix))]<-"Unknown"
#mmgf$V3[which(!mmgf$V3%in%rownames(distancematrix))]
mmgf$V3[which(!mmgf$V3%in%rownames(distancematrix))]<-"Unknown"

#father-paternal grandmother distance
fpgm<-as.data.frame(cbind(demographicdata$Subject_ID,demographicdata$f_birthplace,demographicdata$pgm_birthplace,demographicdata$f_birthyear,demographicdata$pgm_birthyear,demographicdata$f_eg,demographicdata$pgm_eg),stringsAsFactors = F)
#fpgm$V2[which(!fpgm$V2%in%rownames(distancematrix))]
fpgm$V2[which(!fpgm$V2%in%rownames(distancematrix))]<-"Unknown"
#fpgm$V3[which(!fpgm$V3%in%rownames(distancematrix))]
fpgm$V3[which(!fpgm$V3%in%rownames(distancematrix))]<-"Unknown"

#father-paternal grandfather distance
fpgf<-as.data.frame(cbind(demographicdata$Subject_ID,demographicdata$f_birthplace,demographicdata$pgf_birthplace,demographicdata$f_birthyear,demographicdata$pgf_birthyear,demographicdata$f_eg,demographicdata$pgf_eg),stringsAsFactors = F)
#fpgf$V2[which(!fpgf$V2%in%rownames(distancematrix))]
fpgf$V2[which(!fpgf$V2%in%rownames(distancematrix))]<-"Unknown"
#fpgf$V3[which(!fpgf$V3%in%rownames(distancematrix))]
fpgf$V3[which(!fpgf$V3%in%rownames(distancematrix))]<-"Unknown"

#assign distances to parent-offspring pairs
sm<-distassn(sm,distancematrix)
sf<-distassn(sf,distancematrix)
mmgm<-distassn(mmgm,distancematrix)
mmgf<-distassn(mmgf,distancematrix)
fpgm<-distassn(fpgm,distancematrix)
fpgf<-distassn(fpgf,distancematrix)

#rename columns
colnames(sm)<-c("subjectID","subject_birthplace","m_birthplace","subject_birthyear","m_birthyear","subject_eg","m_eg","sm_distance")
colnames(sf)<-c("subjectID","subject_birthplace","f_birthplace","subject_birthyear","f_birthyear","subject_eg","f_eg","sf_distance")
colnames(mmgm)<-c("subjectID","m_birthplace","mgm_birthplace","m_birthyear","mgm_birthyear","m_eg","mgm_eg","mmgm_distance")
colnames(mmgf)<-c("subjectID","m_birthplace","mgf_birthplace","m_birthyear","mgf_birthyear","m_eg","mgf_eg","mmgf_distance")
colnames(fpgm)<-c("subjectID","f_birthplace","pgm_birthplace","f_birthyear","pgm_birthyear","f_eg","pgm_eg","fpgm_distance")
colnames(fpgf)<-c("subjectID","f_birthplace","pgf_birthplace","f_birthyear","pgf_birthyear","f_eg","pgf_eg","fpgf_distance")

#combine parent/offspring distance dataframes
d1<-data.table(sm, key=c("subjectID"))
d2<-data.table(sf, key=c("subjectID"))
d3<-data.table(mmgm, key=c("subjectID"))
d4<-data.table(mmgf, key=c("subjectID"))
d5<-data.table(fpgm, key=c("subjectID"))
d6<-data.table(fpgf, key=c("subjectID"))
mergeddf<-Reduce(merge,list(d1,d2,d3,d4,d5,d6))

#there are some duplicate columns from this merge denoted with ".x",".y",etc...
#We can fix by renaming a few and removing the rest
names(mergeddf)[names(mergeddf) == "subject_birthplace.x"]<-"subject_birthplace"
names(mergeddf)[names(mergeddf) == "subject_birthyear.x"]<-"subject_birthyear"
names(mergeddf)[names(mergeddf) == "subject_eg.x"]<-"subject_eg"
tidymergeddf <- mergeddf %>% select(-contains(".")) #remove duplicate columns


### 4. Ancestry ###
ancestry <- demographicdata[ , c("Subject_ID", "KHS", "BTU","EUR")] #isolate columns of ancestry proportions from ADMIXTURE (Alexander et al. 2009 - Genome Research)
colnames(ancestry)[1] <- "subjectID" #change naming convention for merging
ancestry$subjectID <- as.character(ancestry$subjectID) #Cast subjectID from numeric to character
tidymergeddf<-merge(tidymergeddf, ancestry, by = "subjectID") #merge ancestry with tidymergeddf 

#Cast all ancestry columns as numeric
tidymergeddf$KHS<-as.numeric(tidymergeddf$KHS)
tidymergeddf$BTU<-as.numeric(tidymergeddf$BTU)
tidymergeddf$EUR<-as.numeric(tidymergeddf$EUR)

#define overall max value of any ancestry (for boundary correction)
overall.max.anc <- max(c(tidymergeddf$KHS,tidymergeddf$BTU,tidymergeddf$EUR),na.rm = T)
#define boundary correction
bound.correction <- (1 - overall.max.anc)/2
#bound.correction != 0
#calculate the estimated ratios of African ancestry to Eurasian ancestry for each proband
#value for the proband is associated with all members of the probands family, as they are not measured
tidymergeddf$logratio_KHS_EUR <- log((tidymergeddf$KHS + bound.correction) / (tidymergeddf$EUR + bound.correction))
tidymergeddf$logratio_BTU_EUR <- log((tidymergeddf$BTU + bound.correction) / (tidymergeddf$EUR + bound.correction))


### 5. Create Longform Hierarchical Dataset ###
#Adding in region id to tidymergeddf in the same way we add ancestry
region <- demographicdata[ , c("Subject_ID", "Region")] #isolate relevant columns
colnames(region)[1] <- "subjectID" #change naming convention for merging
region$subjectID <- as.character(region$subjectID) #make column character, was double
tidymergeddf<-merge(tidymergeddf, region, by = "subjectID") #merge 

#make col for "self" distances
tidymergeddf$self_distance<-rep(NA,nrow(tidymergeddf))

#conveniently reorder all columns 
tidymergeddf<-as.data.frame(tidymergeddf[,c("subjectID",
                                            #m = mother, f = father, mgm = maternal grandmother, mgf = maternal grandfather, ...
                                            "subject_birthyear","m_birthyear","f_birthyear","mgm_birthyear","pgm_birthyear","mgf_birthyear","pgf_birthyear",
                                            "subject_birthplace","m_birthplace","f_birthplace","mgm_birthplace","pgm_birthplace","mgf_birthplace","pgf_birthplace",
                                            #sm = proband to mother, sf = proband to father, mmgm = mother to maternal grandmother, ...
                                            "self_distance","sm_distance","sf_distance","mmgm_distance","fpgm_distance","mmgf_distance","fpgf_distance",
                                            #eg = self-identified ethnic group
                                            "subject_eg","m_eg","f_eg","mgm_eg","pgm_eg","mgf_eg","pgf_eg",
                                            #KHS = KhoeSan, BTU = Bantu, EUR = Eurasian (European + East and South Asian)
                                            "KHS","BTU","EUR","logratio_KHS_EUR","logratio_BTU_EUR","Region")])

#convert to long form hierarchical data structure
#making 7 repeats of everything because of the pedigree structure of our data (proband, mother, father, maternal grandmother, ...)
longformids<-rep(tidymergeddf$subjectID,each=7)
longformrelations<-rep(c("self","m","f","mgm","pgm","mgf","pgf"),nrow(tidymergeddf))
longformsubjectbirthyears<-rep(tidymergeddf$subject_birthyear,each=7)
longformKHS<-rep(tidymergeddf$KHS,each=7)
longformBTU<-rep(tidymergeddf$BTU,each=7)
longformEUR<-rep(tidymergeddf$EUR,each=7)
longformlogratio_KHS_EUR<-rep(tidymergeddf$logratio_KHS_EUR,each=7)
longformlogratio_BTU_EUR<-rep(tidymergeddf$logratio_BTU_EUR,each=7)
longformRegion<-rep(tidymergeddf$Region,each=7)
longformsex<-rep(c(NA,"Female","Male","Female","Female","Male","Male"),nrow(tidymergeddf))
transposedbirthyears<-melt(t(tidymergeddf[,c("subject_birthyear","m_birthyear","f_birthyear","mgm_birthyear","pgm_birthyear","mgf_birthyear","pgf_birthyear")]))
transposedbirthplaces<-melt(t(tidymergeddf[,c("subject_birthplace","m_birthplace","f_birthplace","mgm_birthplace","pgm_birthplace","mgf_birthplace","pgf_birthplace")]))
transposeddistances<-melt(t(tidymergeddf[,c("self_distance","sm_distance","sf_distance","mmgm_distance","fpgm_distance","mmgf_distance","fpgf_distance")]))
transposedegs<-melt(t(tidymergeddf[,c("subject_eg","m_eg","f_eg","mgm_eg","pgm_eg","mgf_eg","pgf_eg")]))
longformdf<-as.data.frame(cbind(longformids,longformrelations,as.numeric(as.character(transposedbirthyears$value)),longformsubjectbirthyears,as.character(transposedbirthplaces$value),transposeddistances$value,as.character(transposedegs$value),longformsex,longformKHS,longformBTU,longformEUR,longformlogratio_KHS_EUR,longformlogratio_BTU_EUR,longformRegion),stringsAsFactors = F)
colnames(longformdf)<-c("subjectID","subjectRelation","birthyear","subject_birthyear","birthplace","distance","eg","sex","KHS","BTU","EUR","logratio_KHS_EUR","logratio_BTU_EUR","Region")

#convert distance and year columns to numeric
longformdf$distance<-as.numeric(longformdf$distance)
longformdf$birthyear<-as.numeric(longformdf$birthyear)
longformdf$subject_birthyear<-as.numeric(longformdf$subject_birthyear)
longformdf$KHS<-as.numeric(longformdf$KHS)
longformdf$BTU<-as.numeric(longformdf$BTU)
longformdf$EUR<-as.numeric(longformdf$EUR)
longformdf$logratio_KHS_EUR<-as.numeric(longformdf$logratio_KHS_EUR)
longformdf$logratio_BTU_EUR<-as.numeric(longformdf$logratio_BTU_EUR)


### 6. Censored Migration Data ###
#add new columns to longformdf for max/min distance and whether they needed to be calculated in "censored"
longformdf$censored<-0
longformdf$maxdist<-NA
longformdf$mindist<-NA

#loop through all possible polygons
#polygons and their respective names must be aligned in regional.polygons/names
for (i in 1:length(regional.polygons)){
  #print(names(regional.polygons[i]))
  rctpolygon <- as.matrix(regional.polygons[[i]]) #assign rctpolygon to the ith polygon in the list
  #print(regional.names[i])
  region.name <- regional.names[i] # assign region.name to the ith name in the list
  
  #assign the max/min/censored values for the inds born in region.name
  for (i in 1:nrow(longformdf)){#for every individual in the dataset
    if(longformdf[i,"birthplace"]==region.name&longformdf[i,"subjectRelation"]%in%c("f","m")){#if a mother or father was born in region.name
      #print(longformdf[i,])
      x<-longformdf[which(longformdf$subjectID==longformdf[i,"subjectID"]&longformdf$subjectRelation=="self"),] #assign the self individual to x
      longformdf[i,"censored"]<-1 #add the binary censored marker to column
      #print(x)
      if(is.na(latlondata[which(latlondata$birthplace==x$birthplace),"Longitude"])==FALSE){#if the latlon of birthplace of self isn't missing
        #print(latlondata[which(latlondata$birthplace==x$birthplace),])
        y<-latlondata[which(latlondata$birthplace==x$birthplace),c(2,3)] #assign birthplace lat/lon of self to y
        #longformdf[i,"mindist"]<-dist2Line(y,rctpolygon)[1]/1000
        ## CASE C for f,m ##
        if(point.in.polygon(y[,1],y[,2],rctpolygon[,1],rctpolygon[,2])>0){ #if the birthplace y is in the polygon
          longformdf[i,"mindist"]<-0 #make mindistance of self 0
          longformdf[i,"maxdist"]<-maxdist2Line(y,rctpolygon)[1]/1000 #and max distance of self the max dist between birthplace of self and polygon border
        } else{ 
          ## CASE A for f,m ##
          longformdf[i,"mindist"]<-dist2Line(y,rctpolygon)[1]/1000 #make mindistance of self the min dist between birthplace of self and polygon border
          longformdf[i,"maxdist"]<-maxdist2Line(y,rctpolygon)[1]/1000 #and max distance of self the max dist between birthplace of self and polygon border
        }
      } else{ #if the latlon of birthplace of self is missing
        ## CASE B/E for f,m ##
        longformdf[i,"mindist"]<-0
        longformdf[i,"maxdist"]<-NA
      }
    } else if(longformdf[i,"birthplace"]==region.name&longformdf[i,"subjectRelation"]%in%c("mgm","mgf")){#if a mgm or mgf was born in region.name
      x<-longformdf[which(longformdf$subjectID==longformdf[i,"subjectID"]&longformdf$subjectRelation=="m"),] #assign the m individual to x
      longformdf[i,"censored"]<-1
      #print(x)
      if(is.na(latlondata[which(latlondata$birthplace==x$birthplace),"Longitude"])==FALSE){#if the latlon of birthplace of m isn't missing
        #print(latlondata[which(latlondata$birthplace==x$birthplace),])
        y<-latlondata[which(latlondata$birthplace==x$birthplace),c(2,3)] #assign birthplace lat/lon of m to y
        #longformdf[i,"mindist"]<-dist2Line(y,rctpolygon)[1]/1000
        ## CASE C for mgf,mgm ##
        if(point.in.polygon(y[,1],y[,2],rctpolygon[,1],rctpolygon[,2])>0){ #if the birthplace y is in the polygon
          longformdf[i,"mindist"]<-0 #make mindistance of m 0
          longformdf[i,"maxdist"]<-maxdist2Line(y,rctpolygon)[1]/1000 #and max distance of m the max dist between birthplace of m and polygon border
        } else{ 
          ## CASE A for mgf,mgm ##
          longformdf[i,"mindist"]<-dist2Line(y,rctpolygon)[1]/1000 #make mindistance of m the min dist between birthplace of m and polygon border
          longformdf[i,"maxdist"]<-maxdist2Line(y,rctpolygon)[1]/1000 #and max distance of m the max dist between birthplace of m and polygon border
        }
      } else{ #if the latlon of birthplace of m is missing
        ## CASE B/E for mgf,mgm ##
        longformdf[i,"mindist"]<-0
        longformdf[i,"maxdist"]<-NA
      }
    } else if(longformdf[i,"birthplace"]==region.name&longformdf[i,"subjectRelation"]%in%c("pgm","pgf")){#if a pgm or pgf was born in region.name
      x<-longformdf[which(longformdf$subjectID==longformdf[i,"subjectID"]&longformdf$subjectRelation=="f"),] #assign the f individual to x
      longformdf[i,"censored"]<-1
      #print(x)
      if(is.na(latlondata[which(latlondata$birthplace==x$birthplace),"Longitude"])==FALSE){#if the latlon of birthplace of f isn't missing
        #print(latlondata[which(latlondata$birthplace==x$birthplace),])
        y<-latlondata[which(latlondata$birthplace==x$birthplace),c(2,3)] #assign birthplace lat/lon of f to y
        #longformdf[i,"mindist"]<-dist2Line(y,rctpolygon)[1]/1000
        ## CASE C for pgf,pgm ##
        if(point.in.polygon(y[,1],y[,2],rctpolygon[,1],rctpolygon[,2])>0){ #if the birthplace y is in the polygon
          longformdf[i,"mindist"]<-0 #make mindistance of f 0
          longformdf[i,"maxdist"]<-maxdist2Line(y,rctpolygon)[1]/1000 #and max distance of f the max dist between birthplace of f and polygon border
        } else{ 
          ## CASE A for pgf,pgm ##
          longformdf[i,"mindist"]<-dist2Line(y,rctpolygon)[1]/1000 #make mindistance of f the min dist between birthplace of f and polygon border
          longformdf[i,"maxdist"]<-maxdist2Line(y,rctpolygon)[1]/1000 #and max distance of f the max dist between birthplace of f and polygon border
        }
      } else{ #if the latlon of birthplace of f is missing
        ## CASE B/E for pgf,pgm ##
        longformdf[i,"mindist"]<-0
        longformdf[i,"maxdist"]<-NA
      }
    }
  }
  #Now fix the few instances where the m or f was marked as being born in your region and the gf or gm had an actual location
  for (i in 1:nrow(longformdf)){
    if(is.na(longformdf[i,"maxdist"])==TRUE &  #if the row is missing maxdist
       is.na(longformdf[i,"distance"])==TRUE & #and if the row is missing distance
       longformdf[i,"subjectRelation"]%in%c("pgm","pgf") & #and if they are pgm or pgf
       is.na(latlondata[which(latlondata$birthplace==longformdf[i,"birthplace"]),"Longitude"])==FALSE & #and if their birthplace has a non-NA value of longitude
       longformdf$birthplace[which(longformdf$subjectID==longformdf[i,"subjectID"]&longformdf$subjectRelation=="f")]==region.name){ #and if their son "f" was born in region.name
      y<-latlondata[which(latlondata$birthplace==longformdf[i,"birthplace"]),c(2,3)] #assign a dataframe containing the grandparent birthplace latlong coordinates
      longformdf[i,"censored"] <- 1 #add 1 to censored column
      #print(longformdf[i,])
      if(point.in.polygon(y[,1],y[,2],rctpolygon[,1],rctpolygon[,2])>0){ #if the birthplace y is in the polygon
        longformdf[i,"mindist"]<-0 #make mindistance of pgp 0
        longformdf[i,"maxdist"]<-maxdist2Line(y,rctpolygon)[1]/1000 #and max distance of pgp the max dist between birthplace of pgp and polygon border
      } else{ 
        longformdf[i,"mindist"]<-dist2Line(y,rctpolygon)[1]/1000 #make mindistance of pgp the min dist between birthplace of pgp and polygon border
        longformdf[i,"maxdist"]<-maxdist2Line(y,rctpolygon)[1]/1000 #and max distance of pgp the max dist between birthplace of pgp and polygon border
      }
      #print(longformdf[i,])
    } else if(longformdf[i,"subjectRelation"]%in%c("mgm","mgf") & #or if they are missing maxdist value and mgm or mgf
              is.na(latlondata[which(latlondata$birthplace==longformdf[i,"birthplace"]),"Longitude"])==FALSE & #and if their birthplace has a non-NA value of longitude
              longformdf$birthplace[which(longformdf$subjectID==longformdf[i,"subjectID"]&longformdf$subjectRelation=="m")]==region.name){ #and if their daughter "m" was born in region.name
      y<-latlondata[which(latlondata$birthplace==longformdf[i,"birthplace"]),c(2,3)] #assign a dataframe containing the grandparent birthplace latlong coordinates
      longformdf[i,"censored"] <- 1 #add 1 to censored column
      #print(longformdf[i,])
      if(point.in.polygon(y[,1],y[,2],rctpolygon[,1],rctpolygon[,2])>0){ #if the birthplace y is in the polygon
        longformdf[i,"mindist"]<-0 #make mindistance of m 0
        longformdf[i,"maxdist"]<-maxdist2Line(y,rctpolygon)[1]/1000 #and max distance of mgp the max dist between birthplace of mgp and polygon border
      } else{ 
        longformdf[i,"mindist"]<-dist2Line(y,rctpolygon)[1]/1000 #make mindistance of mgp the min dist between birthplace of mgp and polygon border
        longformdf[i,"maxdist"]<-maxdist2Line(y,rctpolygon)[1]/1000 #and max distance of mgp the max dist between birthplace of mgp and polygon border
      }
      #print(longformdf[i,])
    }
  }
}

#finally loop through every row that isn't self and still has no min value to set them all to zero
for (i in 1:nrow(longformdf)){#for every individual in the dataset
  if(is.na(longformdf[i,"mindist"])==TRUE &  #if the row is missing maxdist
     is.na(longformdf[i,"distance"])==TRUE & #and if the row is missing distance
     longformdf[i,"subjectRelation"]!="self" #and if the ind is not the proband
  ){
    ## CASE B/E for individuals where birthplace is not known for them or their parent/offspring so no distance was calculated##
    longformdf[i,"mindist"]<-0
    longformdf[i,"maxdist"]<-NA
  }
}

#make new column 'migrated' which is a binary if they migrated or not (1 if migrated, 0 if no)
#and another column 'case' to denote the pattern 
#Case D/2 - complete cases, migrated a known distance
longformdf$migrated<-ifelse(longformdf$distance>0,1,0)
longformdf$case<-ifelse(longformdf$distance>0,2,NA)
#Case D/1 - complete cases, did not migrate
longformdf$case[which(longformdf$distance==0)]<-1
#Case B,E/5 - no information, completely missing
longformdf$case[which(is.na(longformdf$distance) & is.na(longformdf$maxdist))]<-5
#Case C/3 - left-censored, because the unknown migration distance is below known maximum (maxdist)
longformdf$case[which(is.na(longformdf$distance) & longformdf$mindist==0 & longformdf$maxdist>0)]<-3
longformdf$migrated[which(is.na(longformdf$distance) & longformdf$mindist==0 & longformdf$maxdist>0)]<-1
#Case A/4 - interval-censored, because the unknown migration distance is between an upper and a lower bound
longformdf$case[which(is.na(longformdf$distance) & longformdf$mindist>0 & longformdf$maxdist>0)]<-4
longformdf$migrated[which(is.na(longformdf$distance) & longformdf$mindist>0 & longformdf$maxdist>0)]<-1



### 7. Derived variables for Modeling ###
#Revert "Unknowns" back to NA so they are not elvauated in the model
is.na(longformdf$birthplace)<-longformdf$birthplace=="Unknown"

#make longform df into a data.table to combine proband and family sex information
longformdf <- data.table(longformdf, key = c("subjectID","subjectRelation"))

#create a data.table containing the proband sex along with subjectID and subject relation to combine with longformdf
assignmentdf <- demographicdata[c("Subject_ID","Sex")]
assignmentdf$subjectRelation <- "self" #add subject relation column
colnames(assignmentdf) <- c("subjectID","sex","subjectRelation")
assignmentdf$sex <- factor(assignmentdf$sex, levels = c("M","F"), labels = c("Male","Female")) #change sex coding to that used in longformdf
assignmentdf <- data.table(assignmentdf, key = c("subjectID","subjectRelation"))

#combine proband sex and family sex information
longformdf<-merge(longformdf,assignmentdf, all = TRUE)
longformdf$sex <- ifelse(!is.na(longformdf$sex.x), longformdf$sex.x, as.character(longformdf$sex.y))

#remove intermediate sex columns
longformdf <- longformdf %>% select(-contains(".")) #remove duplicate columns

#change longformdf back into a dataframe
longformdf <- as.data.frame(longformdf)

#make the male binary column
longformdf$male <- ifelse(longformdf$sex=="Male",1,0)

#Make a column of "1"s as a placeholder for use in constructing model matrices
longformdf$ones<-1 

#prepare a dataset to identify locality codes
#create a dataframe containing the birthplace of the proband and proband's mother and father
mf<-demographicdata[,c("Subject_ID","birthplace","m_birthplace","f_birthplace")]
mf$relation <- "self"

#create a dataframe containing the birthplace of proband's mother with birthplaces of maternal grandmother and maternal grandfather
mgmmgf<-demographicdata[,c("Subject_ID","m_birthplace","mgm_birthplace","mgf_birthplace")]
mgmmgf$relation <- "m"
colnames(mgmmgf)<-colnames(mf)

#create a dataframe containing the birthplace of proband's father with birthplaces of paternal grandmother and paternal grandfather
pgmpgf<-demographicdata[,c("Subject_ID","f_birthplace","pgm_birthplace","pgf_birthplace")]
pgmpgf$relation <- "f"
colnames(pgmpgf)<-colnames(mf)

#combine the birthplace dataframes and add a new column locality
parents<-rbind(mf,mgmmgf,pgmpgf)
parents$locality <- NA
colnames(parents) <- c("subjectID","subject_birthplace","m_bp","f_bp","subjectRelation","locality")

#identify and assign locality codes
#Identify and assign "same-same" locality pattern, where the child shares a birthplace with both the child's mother and father
parents[which(parents$subject_birthplace==parents$m_bp & parents$subject_birthplace==parents$f_bp),"locality"]<- "ss"
#Identify and assign "matrilocal", where the child shares a birthplace with the child's mother but not the child's father
parents[which(parents$subject_birthplace==parents$m_bp & parents$subject_birthplace!=parents$f_bp),"locality"]<- "matrilocal"
#Identify and assign "patrilocal", where the child shares a birthplace with the child's father but not the child's mother
parents[which(parents$subject_birthplace!=parents$m_bp & parents$subject_birthplace==parents$f_bp),"locality"]<- "patrilocal"
#Identify and assign "neolocal", where the child does not share a birthplace with either the child's mother or father
parents[which(parents$subject_birthplace!=parents$m_bp & parents$subject_birthplace!=parents$f_bp),"locality"]<- "neolocal"
#Identify and assign "smuf", where the child shares a birthplace with the child's mother but the father's birthplace is unknown
parents[which(parents$subject_birthplace==parents$m_bp & is.na(parents$f_bp)==TRUE),"locality"]<- "smuf"
#Identify and assign "dmuf", where the child does not share a birthplace with the child's mother but the father's birthplace is unknown
parents[which(parents$subject_birthplace!=parents$m_bp & is.na(parents$f_bp)==TRUE),"locality"]<- "dmuf"
#Identify and assign "dmuf", where the child shares a birthplace with the child's father but the mother's birthplace is unknown
parents[which(parents$subject_birthplace==parents$f_bp & is.na(parents$m_bp)==TRUE),"locality"]<- "umsf"
#Identify and assign "dmuf", where the child does not share a birthplace with the child's father but the mother's birthplace is unknown
parents[which(parents$subject_birthplace!=parents$f_bp & is.na(parents$m_bp)==TRUE),"locality"]<- "umdf"
#Identify and assign "umuf", where the birthplace of both of the child's parents is unknown OR the child's birthplace is unknown
parents[which(is.na(parents$f_bp) & is.na(parents$m_bp)==TRUE),"locality"]<- "umuf"
parents[which(is.na(parents$subject_birthplace) ==TRUE),"locality"]<- "umuf"
parents <- parents[,c("subjectID","subjectRelation","locality","m_bp","f_bp")]

#now add locality column to longformdf
longformdf <- full_join(longformdf,parents)
#define grandparents as umuf
longformdf[which(longformdf$subjectRelation %in% c("mgm","mgf","pgm","pgf")),"locality"]<- "umuf"

#add city size information to dataset
longformdf$size_categorical <- sapply(longformdf$birthplace, function(x) latlondata$designation[match(x, latlondata$birthplace)])
longformdf$size_categorical[which(longformdf$size_categorical %in% c("R","N"))] <- "U"
longformdf$size_categorical <- factor(longformdf$size_categorical, levels = c("C","F","T","U"), labels = c("City","Farm","Town","Unknown"))



### NOW THE DATASET SPLITS INTO TWO ###

### 8. Create Dataset for Migration Model ###
#Remove probands from longformdf because migration status as defined in this study is only known for proband's parents and grandparents
df.migration<-subset(longformdf, subset = subjectRelation != "self")

# Remove records completely missing for migration info. 
df.migration <- subset(df.migration, subset = case != 5)
df.migration <- droplevels(df.migration) 

# Remove records missing ancestry data
# We only check one ratio here because if they are missing information in one ancestry category, they will be missing all information
complete.vec <- complete.cases(df.migration[,"logratio_KHS_EUR"])
df.migration<-df.migration[complete.vec,]
df.migration <- droplevels(df.migration) 

# Make subjectID.numeric and birthplace.numeric, to ensure consecutive integers for nested references in random effect vectors 
df.migration$subjectID.numeric <- as.integer(as.factor(df.migration$subjectID))
df.migration$birthplace.numeric <- as.integer(as.factor(df.migration$birthplace))
# Check. should be true
# max(df.migration$subjectID.numeric) == length(unique(df.migration$subjectID.numeric))  
# max(df.migration$birthplace.numeric) == length(unique(df.migration$birthplace.numeric))  

#adding a new column, generation, for parent or grandparent status
df.migration$generation<-ifelse(df.migration$subjectRelation=="m"|df.migration$subjectRelation=="f",1,0)

# Center subject birthyear on the median of all subject birthyears
df.migration$subject_birthyear_centered <- df.migration$subject_birthyear - median(df.migration$subject_birthyear)
df.migration$subject_birthyear_centered_centennial <- df.migration$subject_birthyear_centered/100



### 9. Create Dataset for Locality Model ###
# Remove records completely missing for migration info. 
df.locality <- subset(longformdf, subset = locality %in% c("ss","neolocal","matrilocal","patrilocal"))
df.locality <- droplevels(df.locality)

# Remove records missing ancestry data
# We only check one ratio here because if they are missing information in one ancestry category, they will be missing all information
complete.vec <- complete.cases(df.locality[,"logratio_KHS_EUR"])
df.locality<-df.locality[complete.vec,]
df.locality <- droplevels(df.locality)

# Make subjectID.numeric and vectors for mother and father's birthplace, to ensure consecutive integers for nested references in random effect vectors 
df.locality$subjectID.numeric <- as.integer(as.factor(df.locality$subjectID))
df.locality$m_bp.numeric <- as.integer(as.factor(df.locality$m_bp))
df.locality$f_bp.numeric <- as.integer(as.factor(df.locality$f_bp))
# Check. should be true
# max(df.locality$subjectID.numeric) == length(unique(df.locality$subjectID.numeric))
# max(df.locality$m_bp.numeric) == length(unique(df.locality$m_bp.numeric))
# max(df.locality$f_bp.numeric) == length(unique(df.locality$f_bp.numeric))

# Redeclare locality as a factor in our desired order
# The LAST level in the vector of levels is the reference level in the Stan model.
# Our preference is to have "ss" as the reference level. 
df.locality$locality <- factor(df.locality$locality, levels = c("neolocal","matrilocal","patrilocal","ss"))
# Now, the first of three effects for each condition (parent, region, etc.) is neolocal : ss and so on...
#check
# table(df.locality$locality,as.numeric(df.locality$locality))

#make column for numeric locality
df.locality$locality.numeric <- as.integer(as.factor(df.locality$locality))

# adding a new column, generation, for parent or grandparent status
# The dataset df.locality contains only fathers, mothers, and probands (subjectRelation == "self")
# This is because we can only determine locality by comparing a child's birthplace with the birthplaces of BOTH parents
# Because we do not have birthplace information for great-grandparents, grandparental records are removed
# The determination of locality will be associated with the child's record
# The locality code on a proband record refers to the locality pattern of the proband's parents
# The locality code on a mother/father record refers to the locality pattern of the maternal/paternal grandparents
df.locality$generation<-ifelse(df.locality$subjectRelation=="m"|df.locality$subjectRelation=="f",0,1)

# Center subject birthyear on the median of all subject birthyears
df.locality$subject_birthyear_centered <- df.locality$subject_birthyear - median(df.locality$subject_birthyear)
df.locality$subject_birthyear_centered_centennial <- df.locality$subject_birthyear_centered/100


### 10. Write datasets to file ###
# Write QC'd longformdf before it was split into dataframes for the models
write.table(longformdf,
            file = paste(outputpath,"longformdf.csv",sep = ""),
            sep = ",",
            col.names = T,
            row.names = F,
            quote = F
)

# Write the dataframe for the migration model
write.table(df.migration,
            file = paste(outputpath,"df.migration.csv",sep = ""),
            sep = ",",
            col.names = T,
            row.names = F,
            quote = F
)

# Write the dataframe for the locality model
write.table(df.locality,
            file = paste(outputpath,"df.locality.csv",sep = ""),
            sep = ",",
            col.names = T,
            row.names = F,
            quote = F
)
