rm(list = ls())
#### Initial data preparation ####
library(readxl)
library(tidyverse)
library(ggplot2)
library(vegan)
library(reshape2)
library(RColorBrewer)
library(clValid)
library(stats)
library(plotly)
library(RColorBrewer)
library(RSelenium)#for saving plotly as svg
library(rcompanion)
library(FSA)
library(ggpubr)
library(factoextra) 
library(ecodist)
library(party) 
library(caret)
library(Rmisc)
library(nlme)
library(multcomp)
library(emmeans)
library(MuMIn)
library(effects)

fish_data <- read_xlsx("Manuscript_data.xlsx","Fish_data")

#setting all numrical columns to be numeric
fish_data$TL_mm <- as.numeric(fish_data$TL_mm)
fish_data$SL_mm <- as.numeric(fish_data$SL_mm)
fish_data$gape_height_mm <- as.numeric(fish_data$gape_height_mm)
fish_data$gape_width_mm <- as.numeric(fish_data$gape_width_mm)
fish_data$stomach_weight_g <- as.numeric(fish_data$stomach_weight_g)
fish_data$caudal_fin_area    <- as.numeric(fish_data$caudal_fin_area)          
fish_data$caudal_fin_height_mm    <- as.numeric(fish_data$caudal_fin_height_mm)          
fish_data$caudal_fin_AR    <- as.numeric(fish_data$caudal_fin_AR)          
fish_data$left_pectoral_fin_area    <- as.numeric(fish_data$left_pectoral_fin_area)          
fish_data$left_pectoral_fin_length_mm    <- as.numeric(fish_data$left_pectoral_fin_length_mm)          
fish_data$left_pectoral_fin_AR    <- as.numeric(fish_data$left_pectoral_fin_AR)          
fish_data$right_pectoral_fin_area    <- as.numeric(fish_data$right_pectoral_fin_area)          
fish_data$right_pectoral_fin_length_mm    <- as.numeric(fish_data$right_pectoral_fin_length_mm)          
fish_data$right_pectoral_fin_AR    <- as.numeric(fish_data$right_pectoral_fin_AR) 

#There are a small number of missing values for TL. We have sufficient data to estimate the relationship between SL and TL for each species,
#and therefore predict the missing TL values

#visualising the relationship:
ggplot(fish_data,aes(x=TL_mm,y=SL_mm))+geom_point()+facet_wrap(~Species)+geom_smooth(method="lm")

#fitting linear models for each species
fitted_models = fish_data %>% 
  group_by(Species) %>% 
  do(model = lm(as.numeric(TL_mm)~as.numeric(SL_mm), data = .,na.action = na.exclude))

#replacing the missing TL values with predictions from SL
for(i in 1:nrow(fitted_models)){
  sp_NA_index <- which(fish_data$Species==fitted_models$Species[i]&is.na(fish_data$TL_mm)==TRUE)
  SL_pred_vals <- fish_data$SL_mm[sp_NA_index]
  fish_data$TL_mm[sp_NA_index] <- predict(fitted_models$model[[i]],data.frame(SL_mm=SL_pred_vals))
}

#### Creating size classes ####
#Here, we split the entire size range of the community into four equal size bins
seq(min(fish_data$TL_mm,na.rm = T),max(fish_data$TL_mm,na.rm = T),length.out=4)
#visualising the bins
ggplot(fish_data,aes(x=TL_mm))+geom_histogram(bins=100)+geom_vline(xintercept = 93,linetype="dashed")+
  geom_vline(xintercept = 268.5,linetype="dashed")+geom_vline(xintercept = 444,linetype="dashed")+
  geom_vline(xintercept = 619.5,linetype="dashed")+geom_vline(xintercept = 795,linetype="dashed")+
  facet_wrap(~Species)

#designating the size classes
fish_data$size_class <- NA
fish_data[which(fish_data$TL_mm <=268.5),]$size_class <- "Size_class_1"
fish_data[which(fish_data$TL_mm >268.5 & fish_data$TL_mm <=444),]$size_class <- "Size_class_2"
fish_data[which(fish_data$TL_mm >444 & fish_data$TL_mm <=619.5),]$size_class <- "Size_class_3"
fish_data[which(fish_data$TL_mm >619.5 & fish_data$TL_mm <=795),]$size_class <- "Size_class_4"

#adding the species ID
fish_data$sp_size_class <- paste(fish_data$Species,fish_data$size_class,sep="_")

#resetting the size class for individuals with missing length data to NA
fish_data[is.na(fish_data$SL_mm)==TRUE,]$sp_size_class <- NA

#### Calculating %IRI for each species-size class combination ####
prey_data <- read_xlsx("Manuscript_data.xlsx","Prey_data")

# Removing rows with no value in the prey number column
prey_data$Number <- as.numeric(prey_data$Number)
prey_data$Total_weight_g <- as.numeric(prey_data$Total_weight_g)
prey_data <- prey_data%>%
  drop_na(Number,Total_weight_g)

#Calculating $IRI by Species and size class
prey_data <- subset(prey_data,prey_data$Number>=1) #removing rows with values less than 1

#filtering out rarely-consumed prey and those difficult to quantify
table(prey_data$Prey_group)
prey_data <- subset(prey_data,prey_data$Prey_group!="unid."& prey_data$Prey_group!="Parasitic Worm" & 
                        prey_data$Prey_group!="Jelly" & prey_data$Prey_group!="Octopus" & prey_data$Prey_group!="Squid" & 
                        prey_data$Prey_group!="Unid. Crustacean")

#separating themisto and other amphipods into different groups
prey_data$Prey_group[which(prey_data$Best_guess_group=="Themisto")] <- "Themisto"
prey_data$Prey_group[which(prey_data$Prey_group=="Amphipod")] <- "Other Amphipods"

#removing net feeding events from the data
prey_data <- prey_data[-which(prey_data$Comments=="Net feed"),]

# grouping brittlestars, worms, bivalves, cumacean, decapod and gastropod into a misc. benthos group
prey_data$Prey_group[which(prey_data$Prey_group=="Bivalve")] <- "Misc. benthos"
prey_data$Prey_group[which(prey_data$Prey_group=="Brittle star")] <- "Misc. benthos"
prey_data$Prey_group[which(prey_data$Prey_group=="Cumacean")] <- "Misc. benthos"
prey_data$Prey_group[which(prey_data$Prey_group=="Decapod")] <- "Misc. benthos"
prey_data$Prey_group[which(prey_data$Prey_group=="Gastropod")] <- "Misc. benthos"
prey_data$Prey_group[which(prey_data$Prey_group=="Worm")] <- "Misc. benthos"

#rounding up prey numbers, as in some cases the number was estimated from the weight of a subsample and therefore includes decimals
prey_data$Number <- ceiling(prey_data$Number) 

# transferring the species-size class information from the fish dataset by matching unique IDs
prey_data$sp_size_class <- NA
prey_data$size_class <- NA

#adding size classes to data
for(i in 1:nrow(prey_data)){
  index <- which(fish_data$Unique_ID==prey_data$Unique_ID[i])
  prey_data$sp_size_class[i] <- fish_data$sp_size_class[index]
  prey_data$size_class[i] <- fish_data$size_class[index]
}

#There are four fish which have neither a TL nor an SL value. These are excluded
prey_data <- prey_data[is.na(prey_data$sp_size_class)==F,] 
fish_data <- fish_data[is.na(fish_data$TL_mm)==F,]

#There are two species-size cass combinations with very few (<5) stomach records. These are removed. 
prey_data <- subset(prey_data,prey_data$sp_size_class!="SSI_Size_class_1"&prey_data$sp_size_class!="PGE_Size_class_3")

# #adding a unique Fish ID number based on the original unique ID (just in case the current format screws things up)
# prey_data <- prey_data%>%
#   group_by(Unique_ID)%>%
#   dplyr::mutate(Fish_ID = cur_group_id())%>%
#   ungroup()

#There are various cases with duplicates of the same prey group in the same stomach, due to differences in digestive state. Here we merge them
prey_data <- prey_data%>%
  group_by(Species,sp_size_class,size_class,Unique_ID,Prey_group)%>%
  dplyr::summarise(prey_count = sum(Number),prey_weight = sum(Total_weight_g),Unique_ID = Unique_ID)

prey_data <- prey_data%>%
  group_by(Unique_ID)%>%
  dplyr::distinct() 

#adding a column identifying how many unique stomachs there are for each fish species and size class
prey_data <- prey_data%>%
  group_by(Species,sp_size_class,size_class)%>%
  dplyr::mutate(n_fish=n_distinct(Unique_ID))

#checking how many individuals there are per species and size class
prey_data%>%
  group_by(Species,sp_size_class,size_class)%>%
  dplyr::summarise(individuals = mean(n_fish))%>%
  print(n=30)

#now summarising the data to identify how often each prey occurs, and the total weight and counts of each prey item
prey_data_summarised <- prey_data%>%
  group_by(Species,sp_size_class,size_class,Prey_group)%>%
  dplyr::summarise(stom_count = (n()),sum_weight = sum(prey_weight),sum_prey_count = sum(prey_count))

#Estimating the frequency of occurance of each prey group
prey_data_summarised$FO <- NA

for(i in 1:nrow(prey_data_summarised)){
  index <- which(prey_data$Species==prey_data_summarised$Species[i]&prey_data$sp_size_class==prey_data_summarised$sp_size_class[i]&prey_data$Prey_group==prey_data_summarised$Prey_group[i])
  prey_data_summarised$FO[i] <- prey_data_summarised$stom_count[i]/prey_data$n_fish[index][1]
}

#now calculating %IRI
prey_data_summarised <- prey_data_summarised%>%
  group_by(Species,sp_size_class)%>%
  dplyr::mutate(W = sum_weight/sum(sum_weight),N = sum_prey_count/sum(sum_prey_count))%>%
  dplyr::mutate(IRI = (N+W)*FO)%>%
  dplyr::mutate(IRI_perc = IRI/sum(IRI))


#converting the %IRI values from long to wide format
wide_data_IRI <- prey_data_summarised[,c(1,2,3,4,12)] %>%
  spread(Prey_group,IRI_perc)

wide_data_IRI[is.na(wide_data_IRI)] <- 0

# #renaming the species size class for better plotting:
# split <- str_split(wide_data_IRI$sp_size_class,"_")
# for(i in 1:nrow(wide_data_IRI)){
#   wide_data_IRI$sp_size_class_short[i] <- paste(split[[i]][1],split[[i]][4],sep=" ")
# }

# Calculating a similarity matrix (by calculating the inverse of the dissimilarity matrix)
dissim_IRI_bray <- vegdist(wide_data_IRI[,-c(1,2,3)],method = "bray",binary=F,diag=F,na.rm=T)
dissim <- as.matrix(dissim_IRI_bray)
dissim <- 1-dissim
diag(dissim) <- 0
dissim[upper.tri(dissim)] <- NA
#plotting the similarity matrix between species-size classes
melt(dissim)%>%
  ggplot(aes(x=Var1,y=Var2,fill=value))+geom_tile()+
  scale_fill_gradientn(colours = c('darkgrey',brewer.pal(11,"Spectral")[6], brewer.pal(11,"Spectral")[11]),
                       values = c(0, 0.0001, 1),na.value="white",limits=c(0,1))+
  labs(x="Species & size class",y="Species & size class",fill="Similarity")+scale_x_continuous(breaks = seq(1:20),labels=paste(wide_data_IRI$sp_size_class))+
  scale_y_continuous(breaks = seq(1:20),labels=paste(wide_data_IRI$sp_size_class))+
  theme(axis.title = element_text(size=18),axis.text.x = element_text(size=13,colour="black",angle=45,hjust=1),axis.text.y = element_text(size=13,colour="black"),
        legend.position = c(0.1,0.8),legend.text = element_text(size=13),legend.title = element_text(size=15),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),panel.border = element_rect(colour = "black",fill=NA))+
  geom_text(aes(label = substr(as.character(value), 1, 5)))

#now converting back to long format to plot the %IRI values
long_data_IRI <- wide_data_IRI%>%
  gather(Prey_group,IRI,Fish:Themisto,factor_key = T)

IRI_plot <- melt(prey_data_summarised[,c(1,2,3,4,12)]) #%IRI
IRI_plot$value <- IRI_plot$value+0.00000001 #needed to prevent some values coming up to fewer decimal places
options(scipen = 10) #needed to stop some values coming up in scientific notation

IRI_plot%>%
  ggplot(aes(x=sp_size_class,y=Prey_group,fill=value))+geom_tile()+
  scale_fill_gradientn(colours = c(brewer.pal(11,"Spectral")[6], brewer.pal(11,"Spectral")[11]),
                       values = c(0, 1),na.value="darkgrey")+
  labs(x="Species",y="Species")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),panel.border = element_rect(colour = "black",fill=NA))+
  geom_text(aes(label = substr(as.numeric(value), 1, 7)))+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

#### Running cluster analysis on diets ####
#substituting NA for 0:
wide_data_IRI[is.na(wide_data_IRI)] <- 0

#setting up a dataframe for the results
clust_results <- as.data.frame(wide_data_IRI$sp_size_class)
names(clust_results) <- "size_class"
clust_results$size_class <- as.factor(clust_results$size_class)

#run cluster dentification
clust_IRI <- hclust(dissim_IRI_bray,method = "complete",members = clust_results$size_class)

#plot dendogram
plot(clust_IRI)

#checking cluster validity at different levels of dissimilarity
h_seq <- seq(0.1,0.9,0.01)
frame <- as.data.frame(h_seq)
frame$index <- NA
for(i in 1:length(h_seq)){
  groups <- cutree(clust_IRI,h=frame$h_seq[i])
  frame$index[i] <- dunn(dissim_IRI_bray,groups)
}
plot(frame$h_seq, frame$index)
#need to choose a suitable dissimilarity. The highest levels are naturally at the lowest dissimilarity
#I need to find the right balance between the validity of the cluster and the ecological interpretability of the clusters
#plotting the dendogram with lines at different dissimilarity levels:
plot(clust_IRI)
abline(a=0.3,b=0,col="orange")
abline(a=0.35,b=0,col="blue")
abline(a=0.4,b=0,col="red")
abline(a=0.45,b=0,col="green")
abline(a=0.5,b=0,col="black")
#values below approximately 0.45 split the two PGE size classes into separate clusters. 
#This is because they show differences in the proportion of mysids versus other shrimps in their diets.
#For our purposes, this distinction is not necessary
#Values below 0.48 split the group of largely fish feeders into two (one with 3 individuals (larger NOR and largest SSI, the other with smaller NOR and larger NOS). 
#This is because the latter have a slightly lower proportion of fish and a larger proportion of “other amphipods” in their diets.
#Otherwise they have quite similar diets to the other fish feeders. 
#We have selected a cutoff of 0.5 as this retains the most general, and ecologically interpretable, groupings.
groups <- cutree(clust_IRI,h=0.5)
clust_results$cluster <- groups

#Some size classes are not present in the diet data despite being in the main data. 
#This is dealt with by subsetting my data to exclude all individuals from size classes not present in the diet data
cluster <- as.data.frame(clust_results)
fish_data_clusters <- fish_data[which(fish_data$sp_size_class %in% clust_results$size_class==TRUE),]

for(i in 1:nrow(fish_data_clusters)){
  index <- which(fish_data_clusters$sp_size_class[i]==clust_results$size_class)
  fish_data_clusters$cluster[i] <- clust_results$cluster[index]
}

# running simper analysis on cluster results
wide_data_IRI_frame <- as.data.frame(as.matrix(wide_data_IRI[,-c(1,2,3)]))
clust_simper <- (vegan::simper(wide_data_IRI_frame,cluster$cluster,permutations = 99,ordered=T))
summary(clust_simper)

#### plotting bipartite food webs with size classes ####

prey_data_summarised$cluster <- NA
prey_data_summarised$Prey_group[which(prey_data_summarised$Prey_group=="Isopod")] <- "Isopods"

#linking diet cluster to data, to use as colours:
for(i in 1:nrow(prey_data_summarised)){
  index <- which(prey_data_summarised$sp_size_class[i]==clust_results$size_class)
  prey_data_summarised$cluster[i] <- clust_results$cluster[index]
}

#reordering the dataset to group clusters together
prey_data_summarised <- (prey_data_summarised[order(prey_data_summarised$cluster, prey_data_summarised$sp_size_class ),])
selected_size_class <- prey_data_summarised

#automating  source and target links
#identifying unique nodes for network
all_links <- as.data.frame(cbind(selected_size_class$sp_size_class,selected_size_class$Prey_group,
                                 selected_size_class$IRI_perc,selected_size_class$cluster))
names(all_links) <- c("target","source","value","group_colour")

all_links$target_num <- all_links$target
all_links$source_num <- all_links$source
all_links$link_colour <- NA

#converting the link information to numeric to reference the node names
all_links_num <- all_links

all_links_num$link_colour[which(all_links_num$source=="Krill")] <- "#1B9E777D"
all_links_num$link_colour[which(all_links_num$source=="Mysid")] <- "#E6AB027D"
all_links_num$link_colour[which(all_links_num$source=="Isopods")] <- "#D95F027D"
all_links_num$link_colour[which(all_links_num$source=="Misc. benthos")] <- "#D95F027D"
all_links_num$link_colour[which(all_links_num$source=="Themisto")] <- "#7570B37D"
all_links_num$link_colour[which(all_links_num$source=="Other Amphipods")] <- "#7570B37D"
all_links_num$link_colour[which(all_links_num$source=="Fish")] <- "#E7298A7D"
all_links_num$link_colour[which(all_links_num$source=="Notocrangon")] <- "#E6AB027D"

all_links_num$group_colour[which(all_links_num$group_colour==1)] <- "#1B9E77"
all_links_num$group_colour[which(all_links_num$group_colour==2)] <- "#D95F02"
all_links_num$group_colour[which(all_links_num$group_colour==3)] <- "#E7298A"
all_links_num$group_colour[which(all_links_num$group_colour==4)] <- "#7570B3"
all_links_num$group_colour[which(all_links_num$group_colour==5)] <- "#E6AB02"

node_names <- c(unique(all_links$target),unique(all_links$source))

for(i in 1:nrow(all_links_num)){
  index_target <- which(node_names==all_links_num$target[i])
  index_source <- which(node_names==all_links_num$source[i])
  
  all_links_num$target_num[i] <- index_target-1
  all_links_num$source_num[i] <- index_source-1
  
}

target_colours <- unique(all_links$target)
source_colours <- unique(all_links$source)

#setting the colours for the fish
for(i in 1:length(target_colours)){
  index_target <- which(all_links_num$target==target_colours[i])
  target_colours[i] <- all_links_num$group_colour[index_target][1]
}

#setting the colours for the prey
source_colours[which(source_colours=="Krill")] <- "#1B9E77"
source_colours[which(source_colours=="Mysid")] <- "#E6AB02"
source_colours[which(source_colours=="Themisto")] <- "#7570B3"
source_colours[which(source_colours=="Fish")] <- "#E7298A"
source_colours[which(source_colours=="Isopods")] <- "#D95F02"
source_colours[which(source_colours=="Other Amphipods")] <- "#7570B3"
source_colours[which(source_colours=="Misc. benthos")] <- "#D95F02"
source_colours[which(source_colours=="Notocrangon")] <- "#E6AB02"

node_colours <- c(target_colours,source_colours)

#removing links that are < 1%
all_links_num <- all_links_num[-which(all_links_num$value<0.01),]

fig <- plot_ly(type = "sankey",arrangement="freeform",
               domain=list(
                 x=c(0,1),
                 y=c(0,1)
               ),
               orientation = "h",
               
               node = list(label = node_names,
                           color = node_colours,
                           pad = 15,thickness = 30,line = list(color = "black",width = 0.5)),
               link = list(target = all_links_num$target_num, #switching these so that the targets are on the right.
                           source = all_links_num$source_num,
                           value =  all_links_num$value,
                           color=all_links_num$link_colour)
)

fig%>%
  layout(margin = list(r = 200,l = 200,t=10,b=50))%>%
  config(fig, toImageButtonOptions = list(format= 'svg', # one of png, svg, jpeg, webp
                                          
                                          filename= 'custom_image',
                                          
                                          height= 500,
                                          
                                          width= 700,
                                          
                                          scale= 1 ))  




#### Adding length-standardised variables ####

fish_data_clusters$species_cluster <- paste(fish_data_clusters$Species,fish_data_clusters$cluster,sep="_")

for(i in 1:length(unique(fish_data_clusters$species_cluster))){
  subset_data <- fish_data_clusters[which(fish_data_clusters$species_cluster==unique(fish_data_clusters$species_cluster)[i]),]
  
  GH_model <- lm(log(gape_height_mm)~log(TL_mm),data=subset_data)
  GW_model <- lm(log(gape_width_mm)~log(TL_mm),data=subset_data)
  caudal_area_model <- lm(log(as.numeric(caudal_fin_area))~log(TL_mm),data=subset_data)
  caudal_height_model <- lm(log(as.numeric(caudal_fin_height_mm))~log(TL_mm),data=subset_data)
  mean_pec_length_model <- lm(log(as.numeric(mean_pec_fin_length_mm))~log(TL_mm),data=subset_data)
  mean_pec_area_model <- lm(log(as.numeric(mean_pec_fin_area))~log(TL_mm),data=subset_data)
  TL_mean <- mean(na.omit(subset_data$TL_mm))
  
  subset_data$gape_width_mm_stand <- subset_data$gape_width_mm*(TL_mean/subset_data$TL_mm)^coef(GW_model)[2]
  subset_data$gape_height_mm_stand <- subset_data$gape_height_mm*(TL_mean/subset_data$TL_mm)^coef(GH_model)[2]
  subset_data$caudal_area_stand <- as.numeric(subset_data$caudal_fin_area)*(TL_mean/subset_data$TL_mm)^coef(caudal_area_model)[2]
  subset_data$caudal_fin_height_mm <- as.numeric(subset_data$caudal_fin_height_mm)*(TL_mean/subset_data$TL_mm)^coef(caudal_height_model)[2]
  subset_data$mean_pec_fin_length_mm <- as.numeric(subset_data$mean_pec_fin_length_mm)*(TL_mean/subset_data$TL_mm)^coef(mean_pec_length_model)[2]
  subset_data$mean_pec_area_stand <- as.numeric(subset_data$mean_pec_fin_area)*(TL_mean/subset_data$TL_mm)^coef(mean_pec_area_model)[2]
  
  subset_data$gape_area_stand <- pi*((subset_data$gape_height_mm_stand/2)*(subset_data$gape_width_mm_stand/2))
  subset_data$pec_fin_AR_stand <- (subset_data$mean_pec_fin_length_mm^2)/subset_data$mean_pec_area_stand
  subset_data$caudal_fin_AR_stand <- (subset_data$caudal_fin_height_mm^2)/subset_data$caudal_area_stand
  
  #adding back to the normal data
  if(i==1){
    fish_data_clusters_stand <- subset_data
  } else if(i>1){
    fish_data_clusters_stand <- rbind(fish_data_clusters_stand,subset_data)
  }
}

#### Plotting general trait-size relationships and doing pairwise comparisons ####
fish_data_clusters_stand$diet_guild <- NA
fish_data_clusters_stand$gape_area_mm <- as.numeric(fish_data_clusters_stand$gape_area_mm)

axis_text_size <- 22
axis_title_size <- 25

#gape
gape_tl_plot <- ggplot(fish_data_clusters_stand,aes(x=log10(TL_mm),y=log10(fish_data_clusters_stand$gape_area_mm),colour=Species))+geom_point(size=4,alpha=0.7)+
  geom_smooth(method="lm",aes(fill=Species),se=FALSE,lwd=1.5)+
  scale_color_manual(values=c("#990F26","#99600F","#54990F","#0F8299","#3D0F99","#333333","#FF9233","#CC7A88","#999999"))+
  #scale_fill_manual(values=c("#990F26","#99600F","#54990F","#0F8299","#3D0F99","#333333","#FF9233","#CC7A88","#999999"))+
  theme_bw()+labs(y="",x="")+
  theme(axis.title.y = element_text(vjust=2),axis.title.x = element_blank(),
        axis.text = element_text(size=axis_text_size,colour="black"),
        axis.title = element_text(size=axis_title_size),legend.position = c(0.7,0.1),
        legend.text = element_text(size=20),legend.title = element_blank())+
  guides(colour = guide_legend(nrow = 3))+
  scale_y_continuous(limits=c(1,4.1),breaks=seq(1,4,1),labels = c("1.0","2.0","3.0","4.0"))+scale_x_continuous(limits=c(1.95,3.05),breaks=seq(2,3,0.5))

#comparing the species pairwise
kruskal_mass <- kruskal.test(log10(fish_data_clusters_stand$gape_area_mm)~as.factor(fish_data_clusters_stand$Species))
kruskal_mass #there is a significant difference

#running dunn test with bonferroni correction (most conservative)
test_dunn <- dunnTest(log10(fish_data_clusters_stand$gape_area_mm)~fish_data_clusters_stand$Species,method = "bonferroni") 
test_dunn = test_dunn$res
test_dunnlist <- cldList(P.adj ~ Comparison, data = test_dunn, threshold = 0.05)
test_dunnlist

gape_boxplot <- ggplot(fish_data_clusters_stand,aes(x=fct_reorder(Species,-log10(gape_area_mm),median),y=log10(gape_area_mm),group=Species))+
  geom_boxplot(width=0.5,lwd=1,fill="#56a6d1")+theme_bw()+labs(y=bquote("Gape area (log10 mm^2)"),x="")+
  theme(axis.title.y = element_text(vjust=2),axis.title.x = element_blank(),
        axis.text = element_text(size=axis_text_size,colour="black"),axis.title = element_text(size=axis_title_size),
        plot.title = element_text(size=20),legend.position = "none")+
  scale_y_continuous(limits=c(1,4.1),breaks=seq(1,4,1),labels = c("1.0","2.0","3.0","4.0"))+
  stat_summary(geom = 'text', (aes(fontface=3)), label = c("b","c","d","a","c","b","a","a","c"), fun= max, hjust=0.5,vjust = -0.5, size = 10)

#caudal fin
caudal_AR_plot <- ggplot(fish_data_clusters_stand,aes(x=log10(TL_mm),y=(caudal_fin_AR),colour=Species))+geom_point(size=4,alpha=0.7)+
  scale_color_manual(values=c("#990F26","#99600F","#54990F","#0F8299","#3D0F99","#333333","#FF9233","#CC7A88","#999999"))+theme_bw()+
  geom_smooth(method="lm",formula = y~poly(x,2),se=F,lwd=1.5)+labs(y=bquote(""),x="")+
  theme(axis.title.y = element_text(vjust=2),axis.title.x = element_blank(),
        axis.text = element_text(size=axis_text_size,colour="black"),
        axis.title = element_text(size=axis_title_size),legend.title = element_blank(),
        plot.title = element_text(size=20),legend.position="none")+
  scale_x_continuous(limits=c(1.95,3.05),breaks=seq(2,3,0.5))+scale_y_continuous(limits=c(0.75,3.25),breaks=seq(0,3,0.5))

kruskal_mass <- kruskal.test((fish_data_clusters_stand$caudal_fin_AR)~as.factor(fish_data_clusters_stand$Species))
kruskal_mass #there is a significant difference

#running dunn test with bonferroni correction (most conservative)
test_dunn <- dunnTest(fish_data_clusters_stand$caudal_fin_AR~fish_data_clusters_stand$Species,method = "bonferroni") 
test_dunn = test_dunn$res
test_dunnlist <- cldList(P.adj ~ Comparison, data = test_dunn, threshold = 0.05)
test_dunnlist

caudal_boxplot <- ggplot(fish_data_clusters_stand,aes(x=fct_reorder(Species,-caudal_fin_AR,median),y=caudal_fin_AR,group=Species))+
  geom_boxplot(width=0.5,lwd=1,fill="#56a6d1")+theme_bw()+  labs(y=bquote("Caudal fin AR"),x="")+
  theme(axis.title.y = element_text(vjust=2),axis.title.x = element_blank(),
        axis.text = element_text(size=axis_text_size,colour="black"),
        axis.title = element_text(size=axis_title_size),legend.position = "none")+
  scale_y_continuous(limits=c(0.75,3.25),breaks=seq(0,3,0.5))+
  stat_summary(geom = 'text', (aes(fontface=3)), label = c("a","c","cd","b","bc","d","bc","a","c"), fun= max, hjust=0.5,vjust = -0.5, size = 10)

#pectoral fin
kruskal_mass <- kruskal.test((fish_data_clusters_stand$mean_pec_fin_AR)~as.factor(fish_data_clusters_stand$Species))
kruskal_mass #there is a significant difference
#running dunn test with bonferroni correction (most conservative)
fish_data_clusters_stand$mean_pec_fin_AR <- as.numeric(fish_data_clusters_stand$mean_pec_fin_AR)

test_dunn <- dunnTest(as.numeric(fish_data_clusters_stand$mean_pec_fin_AR)~fish_data_clusters_stand$Species,method = "bonferroni") 
test_dunn = test_dunn$res
test_dunnlist <- cldList(P.adj ~ Comparison, data = test_dunn, threshold = 0.05)
test_dunnlist

pectoral_plot <- ggplot(fish_data_clusters_stand,aes(x=log10(TL_mm),y=(as.numeric(mean_pec_fin_AR)),colour=Species))+geom_point(size=4,alpha=0.7)+
  scale_color_manual(values=c("#990F26","#99600F","#54990F","#0F8299","#3D0F99","#333333","#FF9233","#CC7A88","#999999"))+theme_bw()+
  geom_smooth(aes(group=Species),method="lm",formula = y~poly(x,2),se=F,lwd=1.5)+
  labs(y=bquote(""),x="Total length (log10 mm)")+
  theme(axis.text = element_text(size=axis_text_size,colour="black"),axis.title = element_text(size=axis_title_size),
        legend.position = "none",legend.title = element_text(size=18),
        legend.text = element_text(size=18),strip.text = element_text(size=20))+guides(colour = guide_legend(nrow = 2))+
  scale_x_continuous(limits=c(1.95,3.05),breaks=seq(2,3,0.5))+scale_y_continuous(limits=c(0.45,2),breaks=seq(0.5,2,0.5))

pectoral_boxplot <- ggplot(fish_data_clusters_stand,aes(x=fct_reorder(Species,-pec_fin_AR_stand,median),y=pec_fin_AR_stand,group=Species))+
  geom_boxplot(width=0.5,lwd=1,fill="#56a6d1")+theme_bw()+  labs(y=bquote("Pectoral fin AR"),x="Species")+
  theme(axis.title.y = element_text(vjust=2),axis.text = element_text(size=22,colour="black"),axis.title = element_text(size=25),
        plot.title = element_text(size=20),legend.position = "none")+
  scale_y_continuous(limits=c(0.45,2),breaks=seq(0.5,2,0.5))+
  stat_summary(geom = 'text', (aes(fontface=3)), label = c("b","cd","cd","a","d","e","a","bc","e"), fun= max, hjust=0.5,vjust = -0.5, size = 10)


ggarrange(gape_boxplot,gape_tl_plot,caudal_boxplot, caudal_AR_plot,pectoral_boxplot,pectoral_plot,
          labels = c("a", "b", "c","d","e","f"),font.label = list(size=30),
          ncol = 2, nrow = 3)

#### Violin plotting and analysis of trait distributions across dietary guilds ####
#adding feeding guild names
fish_data_clusters_stand$diet_guild <- NA
fish_data_clusters_stand$rel_gape_area_log10 <- log10(fish_data_clusters_stand$gape_area_mm/fish_data_clusters_stand$TL_mm)

table(fish_data_clusters_stand$sp_size_class,fish_data_clusters_stand$cluster)

fish_data_clusters_stand$diet_guild <- NA
fish_data_clusters_stand$diet_guild[which(fish_data_clusters_stand$cluster=="1")] <- "Primarily Krill"
fish_data_clusters_stand$diet_guild[which(fish_data_clusters_stand$cluster=="2")] <- "Primarily Benthos"
fish_data_clusters_stand$diet_guild[which(fish_data_clusters_stand$cluster=="3")] <- "Primarily Fish"
fish_data_clusters_stand$diet_guild[which(fish_data_clusters_stand$cluster=="4")] <- "Themisto and Krill"
fish_data_clusters_stand$diet_guild[which(fish_data_clusters_stand$cluster=="5")] <- "Benthic shrimps"

fish_data_clusters_stand$diet_guild <- factor(fish_data_clusters_stand$diet_guild, levels=c("Primarily Krill", "Primarily Benthos", "Primarily Fish","Themisto and Krill","Benthic shrimps"))

axis_text_size <- 22
axis_title_size <- 25

#gape area
kruskal_gape <- kruskal.test(fish_data_clusters_stand$gape_area_stand~as.factor(fish_data_clusters_stand$cluster))
kruskal_gape #there is a significant difference
#running dunn test with bonferroni correction (most conservative)
test_dunn <- dunnTest(log10(fish_data_clusters_stand$gape_area_mm)~as.factor(fish_data_clusters_stand$diet_guild),method = "bonferroni") 
test_dunn = test_dunn$res
test_dunnlist <- cldList(P.adj ~ Comparison, data = test_dunn, threshold = 0.05)
test_dunnlist

fish_data_clusters_stand%>%
  filter(!is.na(log10(gape_area_mm)))%>%
  group_by(diet_guild)%>%
  dplyr::summarise(total=n())

gape_plot <- ggplot(fish_data_clusters_stand,aes(x=fct_reorder(as.factor(diet_guild),-log10(gape_area_mm)),y=log10(gape_area_mm),fill=as.factor(cluster)))+
  geom_boxplot(width=0.5,aes(group = as.factor(cluster)),lwd=1)+
  scale_fill_manual(values=c("#1B9E77","#D95F02","#E7298A","#7570B3","#E6AB02"))+
  theme_bw()+labs(y="Gape area (log10 mm^2)",x="")+
  theme(axis.title.y = element_text(vjust=2),axis.text = element_text(size=axis_text_size,colour="black"),axis.title.y.left = element_text(size=axis_title_size),
        axis.title.x = element_text(size=axis_title_size,vjust = -0.1),legend.position = "none")+
  scale_x_discrete(labels=c("Fish", "Krill","Benthic \n shrimps","Themisto\n & krill","Benthos"))+
  scale_y_continuous(limits=c(1.5,4.5),breaks=seq(1.5,4.5,0.5))+
  stat_summary(geom = 'text', (aes(fontface=3)), label = c("a","b","b","c","d"), fun= max, hjust=0.5,vjust = -0.5, size = 10)+
  stat_summary(geom = 'text', (aes(fontface=3)), label = c("\n\n\n(160)","\n\n\n(315)","\n\n\n(85)","\n\n\n(160)","\n\n\n(103)"), fun= max, hjust=0.5,vjust = -0.5, size = 8)

#standardised caudal AR
kruskal_caudal <- kruskal.test(fish_data_clusters_stand$caudal_fin_AR_stand~as.factor(fish_data_clusters_stand$cluster))
kruskal_caudal #there is a significant difference
#running dunn test with bonferroni correction (most conservative)
test_dunn <- dunnTest(((fish_data_clusters_stand$caudal_fin_AR))~as.factor(fish_data_clusters_stand$diet_guild),method = "bonferroni") 
test_dunn = test_dunn$res
test_dunnlist <- cldList(P.adj ~ Comparison, data = test_dunn, threshold = 0.05)
test_dunnlist

fish_data_clusters_stand%>%
  filter(!is.na(caudal_fin_AR))%>%
  dplyr::group_by(diet_guild)%>%
  dplyr::summarise(total=n())

caudal_plot <- ggplot(fish_data_clusters_stand,aes(x=fct_reorder(as.factor(diet_guild),-(caudal_fin_AR)),y=(caudal_fin_AR) ,fill=as.factor(cluster)))+
  geom_boxplot(width=0.5,aes(group = as.factor(cluster)),lwd=1)+
  scale_fill_manual(values=c("#1B9E77","#D95F02","#E7298A","#7570B3","#E6AB02"))+
  theme_bw()+labs(y=bquote("Caudal fin AR"),x="")+
  theme(axis.title.y = element_text(vjust=2),axis.text = element_text(size=axis_text_size,colour="black"),axis.title.y.left = element_text(size=axis_title_size),
        axis.title.x = element_text(size=axis_title_size,vjust = -0.1),legend.position = "none")+
  scale_x_discrete(labels=c("Fish", "Krill","Themisto\n & krill","Benthos","Benthic \n shrimps"))+
  scale_y_continuous(limits=c(0.9,3.3),breaks=seq(1,3,0.5))+
  stat_summary(geom = 'text', (aes(fontface=3)), label = c("a","a","b","b","c"), fun= max, hjust=0.5,vjust = -0.5, size = 12)+
  stat_summary(geom = 'text', (aes(fontface=3)), label = c("\n\n\n(146)","\n\n\n(280)","\n\n\n(137)","\n\n\n(88)","\n\n\n(74)"), fun= max, hjust=0.5,vjust = -0.5, size = 8)

#standardised pectoral AR
kruskal_pec <- kruskal.test(fish_data_clusters_stand$pec_fin_AR_stand~as.factor(fish_data_clusters_stand$cluster))
kruskal_pec #there is a significant difference
#running dunn test with bonferroni correction (most conservative)
test_dunn <- dunnTest(fish_data_clusters_stand$mean_pec_fin_AR~as.factor(fish_data_clusters_stand$diet_guild),method = "bonferroni") 
test_dunn = test_dunn$res
test_dunnlist <- cldList(P.adj ~ Comparison, data = test_dunn, threshold = 0.05)
test_dunnlist

fish_data_clusters_stand%>%
  filter(!is.na(mean_pec_fin_AR))%>%
  group_by(diet_guild)%>%
  dplyr::summarise(total=n())

pectoral_plot <- ggplot(fish_data_clusters_stand,aes(x=fct_reorder(as.factor(diet_guild),-(mean_pec_fin_AR)),y=mean_pec_fin_AR,fill=as.factor(cluster)))+
  geom_boxplot(width=0.5,aes(group = as.factor(cluster)),lwd=1)+
  scale_fill_manual(values=c("#E7298A","#D95F02","#1B9E77","#7570B3","#E6AB02"))+
  theme_bw()+labs(y=bquote("Pectoral fin AR"),x="Feeding guild")+
  theme(axis.title.y = element_text(vjust=2),axis.text = element_text(size=axis_text_size,colour="black"),axis.title.y.left = element_text(size=axis_title_size),
        axis.title.x = element_text(size=axis_title_size,vjust = -0.1),legend.position = "none")+
  scale_x_discrete(labels=c("Fish", "Krill", "Benthos","Themisto\n & krill","Benthic \n shrimps"))+
  scale_y_continuous(limits=c(0.5,2),breaks=seq(0.5,2.5,0.5))+
  stat_summary(geom = 'text', (aes(fontface=3)), label = c("a","a","b","c","d"), fun= max, hjust=0.5,vjust = -0.5, size = 12)+
  stat_summary(geom = 'text', (aes(fontface=3)), label = c("\n\n\n(147)","\n\n\n(296)","\n\n\n(102)","\n\n\n(143)","\n\n\n(83)"), fun= max, hjust=0.5,vjust = -0.5, size = 8)

ggarrange(gape_plot, caudal_plot, pectoral_plot,
          labels = c("a", "b", "c"),font.label = list(size=30),
          ncol = 1, nrow = 3)

#### PCA of diet groups - by size class: ####
#omitting the species_cluster combos with very few individuals
fish_data_clusters_stand$rel_gape_area_log10 <- (fish_data_clusters_stand$gape_area_mm)/(fish_data_clusters_stand$TL_mm)
fish_data_clusters_stand$TL_log10 <- log10(fish_data_clusters_stand$TL_mm)
fish_data_clusters_stand$gape_area_log10 <- log10(fish_data_clusters_stand$gape_area_mm)
fish_data_clusters_stand$rel_caudal_AR <- log10(fish_data_clusters_stand$caudal_fin_AR/fish_data_clusters_stand$TL_mm)
fish_data_clusters_stand$rel_pectoral_AR <- log10(fish_data_clusters_stand$mean_pec_fin_AR/fish_data_clusters_stand$TL_mm)

fish_data_2 <- fish_data_clusters_stand%>%
  dplyr::select(c(size_class,cluster,gape_area_stand,caudal_fin_AR_stand,
                  pec_fin_AR_stand))

fish_data_2$size_class[which(fish_data_2$size_class=="Size_class_1")] <- 1
fish_data_2$size_class[which(fish_data_2$size_class=="Size_class_2")] <- 2
fish_data_2$size_class[which(fish_data_2$size_class=="Size_class_3")] <- 3
fish_data_2$size_class[which(fish_data_2$size_class=="Size_class_4")] <- 4

fish_data_2$size_class <- as.numeric(fish_data_2$size_class)

fish_data_2 <- na.omit(fish_data_2) #removing rows with missing values, to allow me to plot properly

fish_data_2$cluster <- as.factor(fish_data_2$cluster)

levels(fish_data_2$cluster) <- c("Krill", "Benthos", "Fish","Themisto sp. & krill","Benthic shrimps")

names(fish_data_2) <- c("Size_class","Cluster","Gape area","Caudal AR",
                        "Pectoral AR")


species_PCA <- prcomp(fish_data_2[,c(3,4,5)],
                      center = TRUE,
                      scale = TRUE)

# ggplot(fish_data_clusters_stand,aes(x=log10(TL_mm),y=log10(mean/TL_mm),colour=Species))+geom_point()+geom_smooth()+facet_wrap(~Species)
PCA_plot <- fviz_pca_biplot(species_PCA,axes=c(1,2),pointsize=4,geom = "point",labelsize=8,alpha.ind=0.4,pch=16,mean.point=F,addEllipses = TRUE,ellipse.level=0.8,ellipse.alpha = 0, 
                            col.var="black",habillage = as.factor(fish_data_2$Cluster),title = "",legend.title="Feeding guild",
                            palette = c("#1B9E77","#D95F02","#E7298A","#7570B3","#E6AB02"))+
  theme(axis.text = element_text(size=22,colour="black"),axis.title = element_text(size=25),
        legend.title = element_text(size=22),legend.position=c(0.5,-0.1),legend.text = element_text(size=20),
        plot.margin = unit(c(0, 0, 0.6, 0), 
                           "inches"))+guides(colour=guide_legend(ncol=5))

species_PCA
get_eigenvalue(species_PCA)

#### Random Forest on traits versus diet clusters by species####
fish_data_clusters_stand$size_class[which(fish_data_clusters_stand$size_class=="Size_class_1")] <- 1
fish_data_clusters_stand$size_class[which(fish_data_clusters_stand$size_class=="Size_class_2")] <- 2
fish_data_clusters_stand$size_class[which(fish_data_clusters_stand$size_class=="Size_class_3")] <- 3
fish_data_clusters_stand$size_class[which(fish_data_clusters_stand$size_class=="Size_class_4")] <- 4
fish_data_clusters_stand$size_class <- as.numeric(fish_data_clusters_stand$size_class)

fish_data_2 <- fish_data_clusters_stand%>%
  dplyr::select(c(size_class,cluster,gape_area_stand,caudal_fin_AR_stand,
                  pec_fin_AR_stand))

fish_data_2 <- na.omit(fish_data_2)

names(fish_data_2) <- c("Size_class","Cluster","Gape_area","Caudal_AR","Pectoral_AR")
fish_data_2[,c(3,4,5)] <- scale(fish_data_2[,c(3,4,5)])

#setting up random forest model from party package
my_cforest_control <- cforest_control(teststat = "quad",
                                      testtype = "Univ", 
                                      mincriterion = 0, 
                                      mtry=sqrt(4),
                                      ntree = 10,
                                      replace = FALSE)

TSS_values <- as.data.frame(matrix(nrow=100,ncol=2))
names(TSS_values) <- c("TSS","SD")

gini_values <- as.data.frame(matrix(nrow=100,ncol=3))
sense_spec_table <- as.data.frame(matrix(nrow=100,ncol=5))
colnames(sense_spec_table) <- c("Krill", "Benthos","Fish","Themisto & krill","Benthic shrimps")

start <- Sys.time()
for(i in 1:100){
  sample <- sample(c(TRUE, FALSE), nrow(fish_data_2), replace=TRUE, prob=c(0.7,0.3))
  train <- fish_data_2[sample,c(2,3,4,5)]
  test <- fish_data_2[!sample,c(2,3,4,5)]
  
  my_cforest <- cforest(as.factor(Cluster) ~ ., data = train, controls = my_cforest_control)
  predictions <- predict(my_cforest, newdata=test,OOB=TRUE, type = "response")
  conf_matrix <- confusionMatrix(predictions,as.factor(test$Cluster))
  class_1_TSS <- (conf_matrix$byClass[1,1]+conf_matrix$byClass[1,2]-1)
  sense_spec_table[i,1] <- class_1_TSS
  class_2_TSS <- (conf_matrix$byClass[2,1]+conf_matrix$byClass[2,2]-1)
  sense_spec_table[i,2] <- class_2_TSS
  class_3_TSS <- (conf_matrix$byClass[3,1]+conf_matrix$byClass[3,2]-1)
  sense_spec_table[i,3] <- class_3_TSS
  class_4_TSS <- (conf_matrix$byClass[4,1]+conf_matrix$byClass[4,2]-1)
  sense_spec_table[i,4] <- class_4_TSS
  class_5_TSS <- (conf_matrix$byClass[5,1]+conf_matrix$byClass[5,2]-1)
  sense_spec_table[i,5] <- class_5_TSS
  # class_6_TSS <- (conf_matrix$byClass[6,1]+conf_matrix$byClass[5,2]-1)
  # sense_spec_table[i,6] <- class_6_TSS
  TSS_values[i,1] <- mean(c(class_1_TSS,class_2_TSS,class_3_TSS,class_4_TSS,class_5_TSS))
  TSS_values[i,2] <- sd(c(class_1_TSS,class_2_TSS,class_3_TSS,class_4_TSS,class_5_TSS))
  varimp_cforest <-  varimp(my_cforest,conditional = TRUE)
  gini <- as.data.frame(varimp_cforest)
  colnames(gini_values) <- rownames(gini)
  gini_values[i,] <- gini$varimp_cforest
  
  print(i)
}

stop <- Sys.time()
stop-start

mean(TSS_values[,1],na.rm = T)
mean(TSS_values[,2],na.rm=T)

gini_values <- gini_values[,-1]

gini_plotting <- as.data.frame(matrix(nrow=3,ncol=4))
colnames(gini_plotting) <- c("Trait","Mean_Importance", "lower_CI","upper_CI")
gini_plotting[,1] <- colnames(gini_values)
gini_values <- na.exclude(gini_values)

gini_relative <- gini_values
for(i in 1:nrow(gini_values)){ #figuring out the proportion of importance of each trait for each run
  for(a in 1:ncol(gini_values)){
    gini_relative[i,a] <- gini_values[i,a]/sum(gini_values[i,])
  }
}

for (i in 1:ncol(gini_values)){
  gini_plotting[i,2] <- mean(gini_relative[,i])
  gini_plotting[i,3] <- CI(gini_relative[,i],ci=0.95)[3]
  gini_plotting[i,4] <- CI(gini_relative[,i],ci=0.95)[1]
  
}

gini_plotting <- gini_plotting[order(gini_plotting$Mean_Importance),]
gini_plotting$Trait[1] <- "Pectoral AR"
gini_plotting$Trait[2] <- "Caudal AR"
gini_plotting$Trait[3] <- "Gape area"

gini_plotting$Trait = factor(gini_plotting$Trait, levels=gini_plotting$Trait)

my_palette = (brewer.pal(10, "Spectral")[seq(1,10,1)])

RF_plot <- ggplot(data=gini_plotting,aes(y=Mean_Importance,x=as.factor(Trait)))+geom_point(size=3)+
  theme_bw()+labs()+
  theme(axis.text = element_text(size=22,colour="black"),axis.title.y = element_text(size=25),axis.title.x = element_blank())+
  scale_y_continuous(name = "Relative importance",breaks=seq(0,1,0.2),limits=c(0,1.01),expand=c(0,0))+
  geom_errorbar(aes(y=as.numeric(Trait),ymin=lower_CI,ymax=upper_CI),width=0.2,lwd=0.9)+
  annotate("text", x=0.8, y=1.3, label= "TSS = 0.77 ± 0.12",cex=8) 

ggarrange(PCA_plot, RF_plot,
          labels = c("a", "b"),font.label = list(size=30),
          ncol = 1, nrow = 2,heights = c(1, 0.4))


#### Linear models of prey mass vs predator mass * diet guild ####
prey_data <- read_xlsx("Manuscript_data.xlsx","Prey_data")

#removing rows with no prey weight:
prey_data <- prey_data[-which(prey_data$Total_weight_g=="NA"),]

#subsetting to only keep rows flagged to use for PPMR:
prey_data <- prey_data[-which(prey_data$Use_for_PPMR=="No"),]

#adding a unique fish ID:
prey_data$Unique_ID <- paste(prey_data$Date_sampled_lab,prey_data$Event,prey_data$Species,prey_data$Temp_ID, sep="_")

#removing some prey groups and aggregating others
prey_data <- subset(prey_data,prey_data$Prey_group!="Misc")

prey_data$Prey_group[which(prey_data$Prey_group=="Bivalve")] <- "Misc. benthos"
prey_data$Prey_group[which(prey_data$Prey_group=="Worm")] <- "Misc. benthos"
prey_data$Prey_group[which(prey_data$Prey_group=="Gastropod")] <- "Misc. benthos"
prey_data$Prey_group[which(prey_data$Prey_group=="Cumacean")] <- "Misc. benthos"
prey_data$Prey_group[which(prey_data$Prey_group=="Brittle star")] <- "Misc. benthos"
prey_data$Prey_group[which(prey_data$Prey_group=="Decapod")] <- "Other crustacean"
prey_data$Prey_group[which(prey_data$Prey_group=="Unid. Crustacean")] <- "Other crustacean"

#calculating average prey mass for each group
#removing rows with no prey counts
prey_data <- prey_data[-which(prey_data$Number=="NA"),]
prey_data$Number <- as.numeric(prey_data$Number)

prey_data$total_prey_mass_log10 <- log10(as.numeric(prey_data$Total_weight_g))

#Keeping only the individuals which are also present in the main data (i.e. have length and weight estimates)
prey_data <- prey_data[which(prey_data$Unique_ID%in%fish_data_clusters_stand$Unique_ID==TRUE),]

#averaging prey size by prey group
prey_data$Species <- as.factor(prey_data$Species)
prey_data$Prey_group <- as.factor(prey_data$Prey_group)

prey_data2 <- prey_data%>%
  dplyr::group_by(Unique_ID,Species,Prey_group)%>%
  drop_na(Number)%>%
  dplyr::summarise(mean_prey_mass_log10 = weighted.mean(as.numeric(total_prey_mass_log10),Number))

#matching fish sizes and clusters to prey sizes
prey_data2$predator_mass_log10 <- NA
prey_data2$cluster <- NA

fish_data_clusters_stand$fish_weight_g <- as.numeric(fish_data_clusters_stand$fish_weight_g)

for(i in 1:nrow(prey_data2)){
  index <- which(fish_data_clusters_stand$Unique_ID==prey_data2$Unique_ID[i])
  prey_data2$predator_mass_log10[i] <- log10(fish_data_clusters_stand$fish_weight_g)[index]
  prey_data2$cluster[i] <- fish_data_clusters_stand$cluster[index]
  
}

#identifying the diet group
prey_data2$diet_group <- NA

for(i in 1:nrow(prey_data2)){
  if(prey_data2$cluster[i]==1){
    prey_data2$diet_group[i] <- "Krill"
  } else if (prey_data2$cluster[i]==2){
    prey_data2$diet_group[i] <- "Benthos"
  } else if (prey_data2$cluster[i]==3){
    prey_data2$diet_group[i] <- "Fish"
  } else if (prey_data2$cluster[i]==4){
    prey_data2$diet_group[i] <- "Themisto sp. and krill"
  } else if (prey_data2$cluster[i]==5){
    prey_data2$diet_group[i] <- "Benthic shrimps"
  }
}

prey_data2$diet_group <- as.factor(prey_data2$diet_group)

#dropping any NA predator masses
prey_data2 <- prey_data2%>%
  drop_na(predator_mass_log10)

#fitting ful model with different random effect structures
lmec = lmeControl(msMaxIter = 500, msMaxEval = 500)
mod_gls <- gls(mean_prey_mass_log10~predator_mass_log10*diet_group,data=prey_data2,method = "REML",na.action = na.omit) 
mod_lme <- lme(mean_prey_mass_log10~predator_mass_log10*diet_group,random = list(~1|Prey_group),data=prey_data2,method = "REML",control=lmec,na.action = na.omit) 
mod_lme2 <- lme(mean_prey_mass_log10~predator_mass_log10*diet_group,random = list(~predator_mass_log10|Prey_group),data=prey_data2,method = "REML",control=lmec,na.action = na.omit) 

anova(mod_gls,mod_lme,mod_lme2)

#best model includes prey group as a random effect with random intercept but fixed slope
windows(record=T)
plot(mod_lme)
qqnorm(residuals(mod_lme))
qqline(residuals(mod_lme))

summary(mod_lme)
anova(mod_lme)

#checking the distribution of residuals
plot(factor(prey_data2$diet_group),resid(mod_lme,type="pearson"))
E <- resid(mod_lme)
coplot(E ~ predator_mass_log10 | factor(diet_group), data = prey_data2)

plot(factor(prey_data2$Prey_group),resid(mod_lme,type="pearson"))
E <- resid(mod_lme)
coplot(E ~ predator_mass_log10 | factor(Prey_group), data = prey_data2)

#adding variance structures
lme_ident_diet <- lme(mean_prey_mass_log10~predator_mass_log10*diet_group,random = list(~1|Prey_group),
                      weights = varIdent(form= ~ 1 | diet_group),data=prey_data2,method = "REML",control=lmec,na.action = na.omit) 

lme_ident_prey <- lme(mean_prey_mass_log10~predator_mass_log10*diet_group,random = list(~1|Prey_group),
                      weights = varIdent(form= ~ 1 | Prey_group),data=prey_data2,method = "REML",control=lmec,na.action = na.omit) 

lme_fix <- lme(mean_prey_mass_log10~predator_mass_log10*diet_group,random = list(~1|Prey_group),
               weights = varFixed(~predator_mass_log10),data=prey_data2,method = "REML",control=lmec,na.action = na.omit)

lme_exp <- lme(mean_prey_mass_log10~predator_mass_log10*diet_group,random = list(~1|Prey_group),
               weights = varExp(form=~predator_mass_log10),data=prey_data2,method = "REML",control=lmec,na.action = na.omit)

lme_const <- lme(mean_prey_mass_log10~predator_mass_log10*diet_group,random = list(~1|Prey_group),
                 weights = varConstPower(form=~predator_mass_log10),data=prey_data2,method = "REML",control=lmec,na.action = na.omit)

lme_comb_fix_mass_ident_diet <- lme(mean_prey_mass_log10~predator_mass_log10*diet_group,random = list(~1|Prey_group),
                                    weights = varComb(varFixed(~predator_mass_log10),varIdent(form= ~ 1 | diet_group)),
                                    data=prey_data2,method = "REML",control=lmec,na.action = na.omit)

lme_comb_fix_mass_ident_prey <- lme(mean_prey_mass_log10~predator_mass_log10*diet_group,random = list(~1|Prey_group),
                                    weights = varComb(varFixed(~predator_mass_log10),varIdent(form= ~ 1 | Prey_group)),
                                    data=prey_data2,method = "REML",control=lmec,na.action = na.omit)

lme_comb_exp_mass_ident_diet <-lme(mean_prey_mass_log10~predator_mass_log10*diet_group,random = list(~1|Prey_group),
                                   weights = varComb(varExp(form=~predator_mass_log10),varIdent(form= ~ 1 | diet_group)),
                                   data=prey_data2,method = "REML",control=lmec,na.action = na.omit)


lme_comb_exp_mass_ident_prey <- lme(mean_prey_mass_log10~predator_mass_log10*diet_group,random = list(~1|Prey_group),
                                    weights = varComb(varExp(form=~predator_mass_log10),varIdent(form= ~ 1 | Prey_group)),
                                    data=prey_data2,method = "REML",control=lmec,na.action = na.omit)


lme_comb_const_mass_ident_diet <- lme(mean_prey_mass_log10~predator_mass_log10*diet_group,random = list(~1|Prey_group),
                                      weights = varComb(varConstPower(form=~predator_mass_log10),varIdent(form= ~ 1 | diet_group)),
                                      data=prey_data2,method = "REML",control=lmec,na.action = na.omit)

lme_comb_const_mass_ident_prey <- lme(mean_prey_mass_log10~predator_mass_log10*diet_group,random = list(~1|Prey_group),
                                      weights = varComb(varConstPower(form=~predator_mass_log10),varIdent(form= ~ 1 | Prey_group)),
                                      data=prey_data2,method = "REML",control=lmec,na.action = na.omit)

lme_comb_ident_diet_prey <- lme(mean_prey_mass_log10~predator_mass_log10*diet_group,random = list(~1|Prey_group),
                                weights = varComb(varIdent(form=~1|diet_group),varIdent(form= ~ 1 | Prey_group)),
                                data=prey_data2,method = "REML",control=lmec,na.action = na.omit)


anova_frame <- anova(mod_lme,lme_ident_diet,lme_ident_prey,lme_fix,lme_exp,lme_const,
                     lme_comb_fix_mass_ident_diet,lme_comb_fix_mass_ident_prey,lme_comb_exp_mass_ident_diet,
                     lme_comb_exp_mass_ident_prey,lme_comb_const_mass_ident_prey,lme_comb_const_mass_ident_diet,
                     lme_comb_ident_diet_prey)

anova_frame[order(anova_frame$AIC,decreasing = TRUE),]

#best variance structure using AIC is fixed variance for body size and identity variance structure for prey group
r.squaredGLMM(lme_comb_fix_mass_ident_prey)

summary(lme_comb_fix_mass_ident_prey)
anova(lme_comb_fix_mass_ident_prey)

prey_data2$diet_group <- factor(prey_data2$diet_group,levels=c("Krill","Benthos", "Fish","Themisto sp. and krill","Benthic shrimps"))

#refining the fixed effects
lme_fix_inter <- lme(mean_prey_mass_log10~predator_mass_log10*diet_group,random = list(~1|Prey_group),
                     weights = varComb(varFixed(~predator_mass_log10),varIdent(form= ~ 1 | Prey_group)),data=prey_data2,method = "ML",control=lmec,na.action = na.omit)

lme_fix_add <- lme(mean_prey_mass_log10~predator_mass_log10+diet_group,random = list(~1|Prey_group),
                   weights = varComb(varFixed(~predator_mass_log10),varIdent(form= ~ 1 | Prey_group)),data=prey_data2,method = "ML",control=lmec,na.action = na.omit)

lme_fix_mass <- lme(mean_prey_mass_log10~predator_mass_log10,random = list(~1|Prey_group),
                    weights = varComb(varFixed(~predator_mass_log10),varIdent(form= ~ 1 | Prey_group)),data=prey_data2,method = "ML",control=lmec,na.action = na.omit)

lme_fix_diet <- lme(mean_prey_mass_log10~diet_group,random = list(~1|Prey_group),
                    weights = varComb(varFixed(~predator_mass_log10),varIdent(form= ~ 1 | Prey_group)),data=prey_data2,method = "ML",control=lmec,na.action = na.omit)

lme_fix_null <- lme(mean_prey_mass_log10~1,random = list(~1|Prey_group),
                    weights = varComb(varFixed(~predator_mass_log10),varIdent(form= ~ 1 | Prey_group)),data=prey_data2,method = "ML",control=lmec,na.action = na.omit)

anova_frame <- anova(lme_fix_inter,lme_fix_add,lme_fix_mass,lme_fix_diet,lme_fix_null)
anova_frame[order(anova_frame$AIC,decreasing = TRUE),]

#best model by AIC removes the interaction between predator mass and dietary group
anova(lme_fix_add)
summary(lme_fix_add)

#plotting a partial residuals plot, with individual lines for each group (need to be estimated manually)
partial_model <- data.frame(effect("predator_mass_log10", lme_fix_add, xlevels=list(predator_mass_log10=seq(min(prey_data2$predator_mass_log10),max(prey_data2$predator_mass_log10),length=100))))
prey_data2$resids <- resid(lme_fix_add) + summary(lme_fix_add)$tTable[1] + summary(lme_fix_add)$tTable[2]*prey_data2$predator_mass_log10

intercepts <- summary(lme_fix_add)[[4]]$fixed
intercepts <- intercepts[-2]

abline <- data.frame(matrix(nrow=5,ncol=1))
names(abline)=c("Group")
abline$Group <- levels(prey_data2$diet_group)
abline$xmin <- NA
abline$xmax <- NA
abline$intercept[1] <- intercepts[1]
abline$intercept[2:5] <- intercepts[1]+intercepts[2:5]
abline$slope <- 0.5619859 

for(i in 1:nrow(abline)){
  index <- which(prey_data2$diet_group==abline$Group[i])
  abline$xmin[i] <- min(prey_data2$predator_mass_log10[index])
  abline$xmax[i] <- max(prey_data2$predator_mass_log10[index])
}

abline$ymin <- abline$intercept+(abline$xmin*abline$slope)
abline$ymax <- abline$ymin+((abline$xmax-abline$xmin)*abline$slope)

ggplot(partial_model,aes(x=predator_mass_log10,y=fit))+
  geom_point(prey_data2,mapping=aes(y=resids,x=predator_mass_log10,colour=diet_group),size=4,alpha=0.7)+
  geom_smooth(method="lm",size=1.2,alpha=0.8,colour="black",linetype="dashed")+
  geom_ribbon(aes(ymin=lower, ymax=upper),alpha=0.2,show.legend = F)+theme_classic()+
  geom_segment(data=abline,aes(x=xmin,xend=xmax,y=ymin,yend=ymax),inherit.aes = FALSE,lwd=1.5,colour=c("#1B9E77","#D95F02","#E7298A","#7570B3","#E6AB02"))+
  scale_color_manual(name = "Feeding guild",labels=c("Krill","Benthos", "Fish","Themisto sp. & krill","Benthic shrimps"),values = c("#1B9E77","#D95F02","#E7298A","#7570B3","#E6AB02"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size=26,colour="black"),axis.title.x = element_text(size=30),
        axis.title.y = element_text(size=30,vjust=1.5),legend.position = "bottom",legend.title = element_text(size=27),legend.text = element_text(size=24),legend.spacing.y = unit(0.2, 'cm'))+
  labs(x="Predator mass (log10 g)",y="Prey mass (log10 g)")+scale_x_continuous(limits=c(1,4))+guides(colour = guide_legend(nrow = 3,override.aes=list(fill=NA)))


