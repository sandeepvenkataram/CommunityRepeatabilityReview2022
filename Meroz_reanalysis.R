################
##
## Sandeep Venkataram
## Reanalysis of data from Meroz et al 2021
##
################


library(reshape2)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(reshape2)
library(rstatix)
library(RColorBrewer)

merozFreqs <- read.table(file="Meroz_freqs.txt",sep="\t",header=TRUE)
speciesNames <- names(merozFreqs)[3:18]
merozFreqs <- merozFreqs[merozFreqs$coexistence,]
generations <- unique(merozFreqs$Generation)
cultures <- unique(merozFreqs$ident)
cultureTypes <- unique(merozFreqs$sample)
plyr::count(merozFreqs$sample[merozFreqs$transfer==0])

#test whether initial frequency (t0) has any correlation with species frequency at the next timepoint. If not, we can conduct randomization for null repeatability expectation without controlling for t0 frequency.
merozFreqsConsecTimepoints = data.frame(ident = character(), sample = character(), sample_kind = character(), species = character(), t0Freq = numeric(), t1Freq = numeric())
for(ident in unique(merozFreqs$ident)){
  if(sum(merozFreqs$ident == ident & merozFreqs$transfer == 0)>0 & sum(merozFreqs$ident == ident & merozFreqs$transfer == 2)>0){
    zeroVals = merozFreqs[merozFreqs$ident == ident & merozFreqs$transfer==0,3:18]
    oneVals = merozFreqs[merozFreqs$ident == ident & merozFreqs$transfer==2,3:18]
    goodPositions = zeroVals > 0 | oneVals > 0
    merozFreqsConsecTimepoints<-rbind(merozFreqsConsecTimepoints, data.frame(ident = ident, sample = merozFreqs$sample[merozFreqs$ident == ident][1], sample_kind = merozFreqs$sample_kind[merozFreqs$ident == ident][1], species = speciesNames[goodPositions], t0Freq = zeroVals[goodPositions], t1Freq = oneVals[goodPositions]))
  }
}

#cor.test(merozFreqsConsecTimepoints$t0Freq, merozFreqsConsecTimepoints$t1Freq)
#ggplot(merozFreqsConsecTimepoints)+geom_point(aes(x=t0Freq,y=t1Freq))+theme_classic()+geom_abline(aes(slope=1,intercept=0))+xlim(c(0,1))+ylim(c(0,1))


#functions to estimate parallelism from Meroz data and expectations under neutrality


calculate_s <- function(x,y){
  return(1 - sqrt(sum((x-y)^2)/2)) #normalize by sqrt2 so that values go from 0 to 1
}

numCommunities = 3
pair_expectations = sapply(c(1:100000),function(x){
  p1 <- runif(numCommunities)
  p2 <- 1-p1
  return(mean(sapply(c(1:numCommunities),function(y){
    return(mean(sapply(c(1:numCommunities),function(z){
      return(calculate_s(c(p1[y],p2[y]),c(p1[z],p2[z])))
    })))
  })))
})

trio_expectations = sapply(c(1:100000),function(x){
  p1 <- runif(numCommunities)
  p2 <- sapply(p1,function(y){return(runif(1,min=0,max=1-y))})
  p3 <- 1 - p1 - p2
  return(mean(sapply(c(1:numCommunities),function(y){
    return(mean(sapply(c(1:numCommunities),function(z){
      return(calculate_s(c(p1[y],p2[y],p3[y]),c(p1[z],p2[z],p3[z])))
    })))
  })))
})




outputDF<-data.frame(species_richness = character(), cultureType = character(), generation = numeric(), parallelism = numeric(), culture_pair = character())

#calculate parallelism for each pair of replicates for each community composition at each timepoint

for(generation in generations){
  for(cultureType in cultureTypes){
    myData <-merozFreqs[merozFreqs$sample == cultureType & merozFreqs$Generation == generation,]
    
    if(NROW(myData)<=1){
      next
    }
    for(pop1ID in c(1:(NROW(myData)))){
      pop1Dat = myData[pop1ID,3:18]
      for(pop2ID in c(1:NROW(myData))){
        pop2Dat = myData[pop2ID,3:18]
        parallelism <- calculate_s(pop1Dat,pop2Dat)
        outputDF<-rbind(outputDF,data.frame(species_richness = myData$sample_kind[1],cultureType = cultureType, generation = generation, parallelism = parallelism, culture_pair = paste(myData$ident[pop1ID],myData$ident[pop2ID],sep=",")))
      }
    }
  }
}

#calculate means and SEM of parallelism for each community composition at each timepoint, also grand mean for each species richness type

meanCountDF<-data.frame(species_richness = character(), cultureType = character(), generation = numeric(), meanVal = numeric(), SEM = numeric())
grandMeanDF <- data.frame(species_richness = character(), generation = numeric(), meanVal = numeric(), SEM = numeric())
for(generation in generations){
  for(richness in unique(outputDF$species_richness)){
    myValsAll <- c()
    for(cultureType in unique(outputDF$cultureType)){
      myVals <- outputDF$parallelism[outputDF$species_richness == richness & outputDF$generation == generation & outputDF$cultureType == cultureType]
      
      if(NROW(myVals)>0){
        myValsAll <- c(myValsAll,mean(myVals))
        meanCountDF <- rbind(meanCountDF, data.frame(species_richness = richness, cultureType = cultureType, generation = generation, meanVal = mean(myVals), SEM = sd(myVals)/sqrt(NROW(myVals))))  
      }
    }
    grandMeanDF<-rbind(grandMeanDF,data.frame(species_richness=richness, generation = generation, meanVal = mean(myValsAll), SEM = sd(myValsAll)/sqrt(NROW(myValsAll))))
  }
}


#paired t test between min s and final s for communities with min <s> less than 0.9

pairedTestDF<-data.frame(community = character(), species_richness = character(), R = character(), ttestP = character(), ttestT = character(), ttestDF = character(), tinit = numeric(), tfinal = numeric())
for(community in unique(meanCountDF$cultureType)){
  myMeanVals = meanCountDF[meanCountDF$cultureType==community,]
  minVal = min(myMeanVals$meanVal[myMeanVals$generation < max(myMeanVals$generation)])
  finalVal = myMeanVals$meanVal[max(which(!is.na(myMeanVals$meanVal)))]
  if(minVal < 0.9 & minVal != finalVal){ #if the minimum value is sufficiently low and is not the final value
    myGen = myMeanVals$generation[myMeanVals$meanVal == minVal]
    myInitComparisonVals = outputDF[outputDF$cultureType==community & outputDF$generation==myGen,]
    myFinalComparisonVals = outputDF[outputDF$cultureType==community & outputDF$generation==max(myMeanVals$generation),]
    mergedData <- merge(myInitComparisonVals,myFinalComparisonVals,by="culture_pair")
    if(NROW(mergedData) > 1 & sum(!is.na(mergedData$parallelism.x) & !is.na(mergedData$parallelism.y))>1){
      ttestRes <-  t.test(mergedData$parallelism.x,mergedData$parallelism.y,paired=TRUE)
      pairedTestDF<-rbind(pairedTestDF, data.frame(community = community, species_richness = myMeanVals$species_richness[1], R = cor(mergedData$parallelism.x,mergedData$parallelism.y), ttestP =ttestRes$p.value, ttestT = ttestRes$statistic, ttestDF = ttestRes$parameter,tinit = myGen, tfinal = max(myMeanVals$generation) ))
    }
  }
}

#these are communities to highlight as significant
mySignifCommunities <- pairedTestDF[which(p.adjust(pairedTestDF$ttestP,method="BH")<0.05),] #communities with significant increases 



#make plots!
myColors <- RColorBrewer::brewer.pal(4,"Paired")

p1<-ggplot()+
  geom_line(data = meanCountDF[meanCountDF$species_richness=="Pair" & ! meanCountDF$cultureType %in% mySignifCommunities$community,], aes(x=generation,y=meanVal,group=cultureType, col = cultureType %in% mySignifCommunities$community),size=1)+
  geom_line(data = meanCountDF[meanCountDF$species_richness=="Pair" & meanCountDF$cultureType %in% mySignifCommunities$community,], aes(x=generation,y=meanVal,group=cultureType, col = cultureType %in% mySignifCommunities$community),size=1)+
  geom_line(data = grandMeanDF[grandMeanDF$species_richness == "Pair",], aes(x=generation,y=meanVal),col="red",size=2)+
  theme_classic()+
  theme(legend.position="none",plot.title=element_text(size=18,color="black"),axis.text=element_text(size=14,color="black"),axis.title=element_text(size=16),panel.grid=element_blank())+
  xlab("Generation")+
  ylab("Repeatability, <s>")+
  geom_tile(aes(x=max(outputDF$generation/2),y=mean(pair_expectations),width=max(outputDF$generation)+40,height=sd(pair_expectations)*2/sqrt(NROW(pair_expectations))),col="grey45",alpha=.5)+
  ggtitle("Two-species communities")+
  scale_color_manual(values=c(myColors[1],myColors[2]))



p2<-ggplot()+
  geom_line(data = meanCountDF[meanCountDF$species_richness=="Trio" & ! meanCountDF$cultureType %in% mySignifCommunities$community,], aes(x=generation,y=meanVal,group=cultureType, col = cultureType %in% mySignifCommunities$community), size=1)+
  geom_line(data = meanCountDF[meanCountDF$species_richness=="Trio" &  meanCountDF$cultureType %in% mySignifCommunities$community,], aes(x=generation,y=meanVal,group=cultureType, col = cultureType %in% mySignifCommunities$community), size=1)+
  geom_line(data = grandMeanDF[grandMeanDF$species_richness == "Trio",], aes(x=generation,y=meanVal),col="red",size=2)+
  theme_classic()+
  theme(legend.position="none",plot.title=element_text(size=18,color="black"),axis.text=element_text(size=14,color="black"),axis.title=element_text(size=16),panel.grid=element_blank())+
  xlab("Generation")+
  ylab("Repeatability, <s>")+
  geom_tile(aes(x=max(outputDF$generation/2),y=mean(trio_expectations),width=max(outputDF$generation)+40,height=sd(trio_expectations)*2/sqrt(NROW(trio_expectations))),col="grey45",alpha=.5)+
  ggtitle("Three-species communities")+
  scale_color_manual(values=c(myColors[1],myColors[2]))



p3<-ggarrange(p1,p2,nrow=1,labels=c("A","B"))
ggsave(p3,file="../Figures/Meroz_timecourse.png",width=200,height=100,units="mm")
ggsave(p3,file="../Figures/Meroz_timecourse.svg",width=200,height=100,units="mm")


## repeating measures ANOVA tests reported in manuscript text and methods


anova_test(data=outputDF[outputDF$species_richness=="Pair",],dv=parallelism,wid = culture_pair, within=generation,between=cultureType)
anova_test(data=outputDF[outputDF$species_richness=="Pair" & outputDF$generation>70,],dv=parallelism,wid = culture_pair, within=generation,between=cultureType)


anova_test(data=outputDF[outputDF$species_richness=="Trio",],dv=parallelism,wid = culture_pair, within=generation,between=cultureType)
anova_test(data=outputDF[outputDF$species_richness=="Trio" & outputDF$generation>70,],dv=parallelism,wid = culture_pair, within=generation,between=cultureType)

#is parallelism different between two and three species communities?
anova_test(data=outputDF,dv=parallelism,wid = culture_pair, within=generation,between=species_richness)
