# DAPC code

# Load necessary packages
library(ade4)
library("poppr")

# Choose the directory where your data file is

setwd("~/Dropbox/UPRRP_Maestria/Investigacion_Iguana/Iguana_Source_Manuscript/Iguana_source_manuscript_data/Source_manuscript_microsatellites/Manuscript_Analysis_MSat/DAPC_Test_20181205/DAPC_Iguana_test")
###### L6 subset #########
# Load datafile
Loci_L6 <- read.genalex("loci6_sampling_sites_20180822.txt", sep = "\t") 


# Transform for accurate reading
Loci_L6 <- genclone2genind(Loci_L6)

Loci_L6 

######optim.a to determin the n.pca to choose for L6########

dapc2 <- dapc(Loci_L6, n.da =100, n.pca = 50)
temp <- optim.a.score(dapc2) 
#based on dapc2 and temp n.pca should for DAPC should be between 10-20

dapc3 <- dapc(Loci_L6, n.da =100, n.pca = 15)
dapc3


#### Sequential K-means algorithm to identify minimum number of clusters
set.seed(3363)
GroupInds_L6 <- find.clusters(Loci_L6, max.n.clust = 40, n.pca = 100) # choose a n.pca value that allows you to encapsulate over 90% of the cumulative variance.  
                                                                      

# Choose number of clusters you are interested in to see
6
GroupInds_L6

#assignment table sampling sites vs cluster

x<-read.genalex("loci6_sampling_sites_20180822.txt" , sep = "\t")
x
table(pop(x), GroupInds_L6$grp)

table.value(table(pop(x), GroupInds_L6$grp), col.lab=paste("Cluster", 1:6), row.lab=paste("ori", 1:10))


# Next step, new pca cut off

dapcGrouping_L6 <- dapc(Loci_L6, GroupInds_L6$grp, n.pca = 15, n.da = 9) # n.pca Based on alpha value, here between 10 and 20. n.da based on discriminant functions saved in dapc3

dapcGrouping_L6



# This code allows you to see which individual is assigned to which population, as its not possible to see that in the standard plot


assignplot(dapcGrouping_L6,subset=1:169) #  number of individuals


# Adjust colors for clusters in plot
myCol <- funky(6)



#Scatter plot of L6

L6_Scatterplot<- scatter(dapcGrouping_L6, posi.da="bottomleft", bg="white", 
        pch=15:21, cstar=0, col=myCol, solid = .6, cex = 2, clab=0, 
        scree.pca=FALSE, posi.pca="bottomright",
        leg=TRUE, txt.leg=paste("Cluster",1:6))


mtext("L6 DAPC clusters", side = 3, line = -1, outer = FALSE ) #title in the margins of plot
title("L6 DAPC clusters")



#structure like figure for L6


L6_Barplot<- compoplot(dapcGrouping_L6 , posi= list(x=-5,y=-.01),
          txt.leg=paste("Cluster", 1:6), lab="",
          ncol=2, xlab="individuals", col=funky(6))



# which are the most ???admixed??? individuals?



temp <- which(apply(dapcGrouping_L6$posterior,1, function(e) all(e<0.9)))
temp


L6_AdmixturePlot<-compoplot(dapcGrouping_L6, subset=temp, posi=list(x=-3,y=-.05),
          txt.leg=paste("Cluster", 1:6),
          ncol=2, col=funky(6))







################# W6 subset analysis #############


# Choose the directory where your data file is

setwd("~/Dropbox/UPRRP_Maestria/Investigacion_Iguana/Iguana_Source_Manuscript/Iguana_source_manuscript_data/Source_manuscript_microsatellites/Manuscript_Analysis_MSat/DAPC_Test_20181205/DAPC_Iguana_test")

###### W6 subset #########
# Load datafile
Loci_W6 <- read.genalex("W6_sampling_sites_20180822.txt", sep = "\t") 


# Transform for accurate reading
Loci_W6 <- genclone2genind(Loci_W6)

Loci_W6 

######optim.a to determin the n.pca to choose for W6########
dapc3 <- dapc(Loci_W6, n.da =100, n.pca = 50)
temp2 <- optim.a.score(dapc3) 

 
#based on dapc2 and temp n.pca should for DAPC should be between 1-10

dapc4 <- dapc(Loci_L6, n.da =100, n.pca = 5)
dapc4
#results  5 discriminant functions 

#### Sequential K-means algorithm to identify minimum number of clusters
set.seed(3365)
GroupInds_W6 <- find.clusters(Loci_W6, max.n.clust = 40, n.pca = 100) # choose a n.pca value, I normally take the # of PC's at between 80-95% of Cumulative variance 

# Choose number of clusters you are interested in to see based on BIC valuees
3
GroupInds_W6

#assignment table sampling sites vs cluster
y<-read.genalex("W6_sampling_sites_20180822.txt", sep = "\t") 
table(pop(y), GroupInds_W6$grp)

table.value(table(pop(y), GroupInds_W6$grp), col.lab=paste("Cluster", 1:3), row.lab=paste("ori", 1:10))


# Next step, new pca cut off based on alpha value using the optim.a function

dapcGrouping_W6 <- dapc(Loci_W6, GroupInds_W6$grp, n.pca = 5, n.da = 5) 

dapcGrouping_W6



# This code allows you to see which individual is assigned to which population, as its not possible to see that in the standard plot


assignplot(dapcGrouping_W6, subset=1:41) # number of individuals in data set


# Adjust colors for clusters in plot
#myCol <- c("orange","blue","pink","brown","lightblue","purple") 
myCol <- funky(3)



#Scatter plot of W6

W6_Scatterplot<- scatter(dapcGrouping_W6, posi.da="bottomright", bg="white", 
                         pch=15:21, cstar=0, col=myCol, solid = .6, cex = 2, clab=0, scree.pca=FALSE, posi.pca="bottomright",
                         leg=TRUE, txt.leg=paste("Cluster",1:3))

#mtext("W6 DAPC clusters", side = 3, line = 0, outer = FALSE)
#title("W6 DAPC clusters ")



#structure like figure for W6


W6_Barplot<- compoplot(dapcGrouping_W6 ,posi= list(x=4,y=-.01),
                       txt.leg=paste("Cluster", 1:3), lab="",
                       ncol=2, xlab="individuals", col=funky(3))


# which are the most ???admixed??? individuals?



temp <- which(apply(dapcGrouping_W6$posterior,1, function(e) all(e<0.9)))
temp


W6_AdmixturePlot<-compoplot(dapcGrouping_W6, subset=temp, posi= list(x=4,y=-.01),
                            txt.leg=paste("Cluster", 1:3),
                            ncol=2, col=funky(3))

