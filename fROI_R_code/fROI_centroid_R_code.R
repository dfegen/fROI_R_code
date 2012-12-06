#############################################################################
#############################################################################
#read in the data:

library("RSQLite")

####to read in a .sqlite file
con<-dbConnect(drv="SQLite",dbname="FROIAtlas.sqlite")
tables<-dbListTables(con)
#exclude sqlite_sequence
tables<-tables[tables != "sqlite_sequence"]
lDataFrames<-vector("list",length=length(tables))

for (i in seq(along=tables)) {
  lDataFrames[[i]] <- dbGetQuery(conn=con, statement=paste("SELECT * FROM '", tables[[i]], "'", sep=""))
}


number_of_rows<-sapply(lDataFrames,NROW)
foci_table_rows<-number_of_rows[1]
study_table_rows<-number_of_rows[2]

foci_table<-lDataFrames[1]
study_table<-lDataFrames[2]

#all columns read into dframe as FACTORS, so need to convert to numeric
foci.dframe<-data.frame(matrix(unlist(foci_table),nrow=foci_table_rows,byrow=F))
study.dframe<-data.frame(matrix(unlist(study_table),nrow=study_table_rows,byrow=F))


names(foci.dframe)<-c('PMID','X','Y','Z','Hemisphere','Key')
foci.dframe$PMID<-as.numeric(as.character(foci.dframe$PMID))
foci.dframe$X<-as.numeric(as.character(foci.dframe$X))
foci.dframe$Y<-as.numeric(as.character(foci.dframe$Y))
foci.dframe$Z<-as.numeric(as.character(foci.dframe$Z))
foci.dframe$Hemisphere<-as.factor(foci.dframe$Hemisphere)
foci.dframe$Key<-as.factor(foci.dframe$Key)

names(study.dframe)<-c('Author','Year','SampleSize','CoordSpace','PMID','Journal','FROI')






#############################################################################
#############################################################################
#analyze the data:

library("car")
library("scatterplot3d")
library("MASS")
library("plyr")

#convert Tal to MNI (and un-Brett Matthew Brett transforms) - but study.df has Tal/MNI
# --> for now will just use on Tal & MNI w/o converting





#segment into left and right - need to make sure L/R naming is consistent in dframe
foci.right.dframe<-subset(foci.dframe,Hemisphere=='Right' | Hemisphere=='R')
foci.left.dframe<-subset(foci.dframe,Hemisphere=='Left' | Hemisphere=='L')


#calculate the mean
mean.foci.right<-as.vector(apply(foci.right.dframe[,2:4],MARGIN=2,FUN=mean))
mean.foci.left<-as.vector(apply(foci.left.dframe[,2:4],MARGIN=2,FUN=mean))

#calculate the median
median.foci.right<-as.vector(apply(foci.right.dframe[,2:4],MARGIN=2,FUN=median))
median.foci.left<-as.vector(apply(foci.left.dframe[,2:4],MARGIN=2,FUN=median))

#calculate the huber mean - in MASS package use 'huber' or 'hubers'???
huber.foci.right<-apply(foci.right.dframe[,2:4],MARGIN=2,FUN=huber)
huber.foci.right<-c(huber.foci.right[[1]][[1]],huber.foci.right[[2]][[1]],huber.foci.right[[3]][[1]])

hubers.foci.right<-apply(foci.right.dframe[,2:4],MARGIN=2,FUN=hubers)
hubers.foci.right<-c(hubers.foci.right[[1]][[1]],hubers.foci.right[[2]][[1]],hubers.foci.right[[3]][[1]])

huber.foci.left<-apply(foci.left.dframe[,2:4],MARGIN=2,FUN=huber)
huber.foci.left<-c(huber.foci.left[[1]][[1]],huber.foci.left[[2]][[1]],huber.foci.left[[3]][[1]])

hubers.foci.left<-apply(foci.left.dframe[,2:4],MARGIN=2,FUN=hubers)
hubers.foci.left<-c(hubers.foci.left[[1]][[1]],hubers.foci.left[[2]][[1]],hubers.foci.left[[3]][[1]])


###calculate Euclidean distance from centroid measure:
#formula for Euclidean distance:
## result<-sqrt(sum((vector1-vector2)^2))

foci.right.dframe.Edist<-ddply(foci.right.dframe,.(X,Y,Z),transform,Eucldist_huber=sqrt(sum((c(X-huber.foci.right[1],Y-huber.foci.right[2],Z-huber.foci.right[3]))^2)))
foci.left.dframe.Edist<-ddply(foci.left.dframe,.(X,Y,Z),transform,Eucldist_huber=sqrt(sum((c(X-huber.foci.left[1],Y-huber.foci.left[2],Z-huber.foci.left[3]))^2)))
#note: ddply() automically sorts the output from smallest to largest according to X -> no way to stop this with ddply()









###put individual coordinates and centroid measures in one dataframe:
#right
results.foci.right.dframe<-foci.right.dframe.Edist[,2:4]
length<-nrow(results.foci.right.dframe)

total.results.right<-rbind(results.foci.right.dframe,mean.foci.right,median.foci.right,huber.foci.right)
total.results.right[,"Coord"]<-as.factor(c(rep('individual',length),'mean','median','huber'))

length<-nrow(total.results.right)
total.results.right[,"Hemi"]<-as.factor(c(rep('right',length)))

#left
results.foci.left.dframe<-foci.left.dframe.Edist[,2:4]
length<-nrow(results.foci.left.dframe)

total.results.left<-rbind(results.foci.left.dframe,mean.foci.left,median.foci.left,huber.foci.left)
total.results.left[,"Coord"]<-as.factor(c(rep('individual',length),'mean','median','huber'))

length<-nrow(total.results.left)
total.results.left[,"Hemi"]<-as.factor(c(rep('left',length)))

#combine LEFT and RIGHT
total.results<-rbind(total.results.right,total.results.left)









###plot
levels(total.results$Coord) <- sub("individual", "black", levels(total.results$Coord))
levels(total.results$Coord) <- sub("mean", "green", levels(total.results$Coord))
levels(total.results$Coord) <- sub("median", "red", levels(total.results$Coord))
levels(total.results$Coord) <- sub("huber", "blue", levels(total.results$Coord))



scatterplot3d(total.results$X,total.results$Y,total.results$Z,color=total.results$Coord,pch=19,main="fROI: MT",xlab="X",ylab="Y",zlab="Z")
legend("bottomright",bty="o",bg="gray",title="Centroid:",cex=0.5,c("mean","median","huber"),fill=c("green","red","blue"))
#we could denote outliers as some color...(or same color and instead label them)
#be nice to make interactive -> when move cursor over point it tells you what study it was from
#--> for 2D scatterplots, functions identify() or locator() can be used to do this...
#would be nice for each ROI to give back some metric of how close all the points are


#note there seems to be different near overlap of colors in the 1st vs. 2 & 3 plots... -> something wrong with color scheme of these two below plots
library(rgl)
plot3d(total.results$X,total.results$Y,total.results$Z,xlab="X",ylab="Y",zlab="Z",col=total.results$Coord)

scatter3d(total.results$X,total.results$Y,total.results$Z,xlab="X",ylab="Y",zlab="Z",point.col=total.results$Coord,surface=FALSE)
#scatter3d() seems more to plot 3D regression
#this program plots spheres, flag can be used to plot points instead
#---> for each point size of sphere could be # of subjects -> could have button to change what the size of sphere represents

#python seems to have better 3D-scatterplots -> more "ggplot2" like

#what about plotting 3D distributions - can't find much on the internet:
#   https://stat.ethz.ch/pipermail/r-help/2012-March/305278.html
#   functions persp() and wireframe()

########
#notes:

#robust centroid
#find huber mean of coordinates:
#  for each coordinate get Euclidean distance to the center (sum of squared differences to the center)
#  then weight it somehow as inverse of the distance
#  then put a Gaussian around the center of each point -> ones with low weight multiply by weight
#  then sum across all points -> will give us the distribution, should be centered around the centroid
#  ? is how going to threshold? could pick arbitrary 0.5 -> or could take all the values and normalize them so they sum 1 and then try diff thresholds until make sense
#  --> could threshold at the median, anything below the median...
