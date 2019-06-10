##Practicing TPC

rm(list=ls())

##Install packages
# load packages
library(nls.multstart)
library(broom)
library(purrr)
library(tidyverse)
library(nlstools)
library(nls2)
library(ggforce)
library(gridExtra)



#load data

photo.data <- read.csv("GalapagosRates.csv")
photo.data$X <- NULL
#View(photo.data)
#glimpse(photo.data)

# remove the NAs from the data
#photo.data<-photo.data[-which(is.na(photo.data$umol.cm2.hr)),]

# remove the three organisms that are wrong (had too much messy respiration files)
#remove<-c('Egala_Bart_1','Egala_Ibbet_3','Egala_Botel_2', 'Egala_Corm_10')

#bad.ID<-which(photo.data$Organism.ID %in% remove)
#photo.data<-photo.data[-bad.ID,]



mydata <- photo.data
## clean some of the bad data
#bad<-which(mydata$umol.cm2.hr< -40)
bad1<-c(which(mydata$Species=='Ptube' & mydata$umol.cm2.hr>70), which(mydata$Species=='Tcoc' & mydata$umol.cm2.hr>75),
        which(mydata$Species=='Lsemi' & mydata$Temp.Cat<20 & mydata$umol.cm2.hr>25),
        which(mydata$Species=='Secu' & mydata$umol.cm2.hr>300 & mydata$Light_Dark=='Dark'), which(mydata$umol.cm2.hr< -40))
mydata<-mydata[-bad1,]

Sp.info<-read.csv('SpeciesList.csv')

#join with the mydata
mydata<-left_join(mydata,Sp.info)

# read in the species info 
# remove respiration rates that are less then 0 and make them 0

mydata$umol.cm2.hr[mydata$Light_Dark=='Dark']<-ifelse(mydata$umol.cm2.hr[mydata$Light_Dark=='Dark']<0,0,mydata$umol.cm2.hr[mydata$Light_Dark=='Dark'])

#mydata$log.rate <- log(mydata$umol.cm2.hr+1)  #logging and adding 0.1 because a log of zero does not exist

# make GP by adding light to dark (NP to R)

#light will be assigned NP for net photosynthesis 
mydata$rate.type <-ifelse(mydata$Light_Dark=='Light', "NP", "R")
mydata$rate.type<-as.factor(mydata$rate.type)
View(mydata)


# add a fragment ID
mydata$Fragment.ID<-factor(paste0(mydata$Organism.ID, "_", mydata$rate.type))

#rename dataframe to work from Photo.R2 to preserve Photo.R
#make column for GP and group by fragment ID and temp to keep R and NP together
Photo.R2 <- mydata %>%
  filter(Functional.group=='Producer') %>%  # only do this for the producers since they are the only ones that photosynthesize
  group_by(Organism.ID, Temp.Cat, Species) %>% 
  summarize(rates = sum(umol.cm2.hr), Temp.C=mean(Temp.C)) %>%
  mutate(rate.type="GP", Light_Dark="Light") %>%  
  rename(umol.cm2.hr=rates) %>%
  mutate(Fragment.ID=paste0(Organism.ID, "_", rate.type)) %>%
  left_join(.,Sp.info)

#rename light darks
#Photo.R2$light_dark <- ifelse(Photo.R2$light_dark=="L", "light", "dark")

#make the negative gphoto values zeros
Photo.R2$umol.cm2.hr[Photo.R2$umol.cm2.hr<0]<-0

# merge together with mydata

mydata<-bind_rows(mydata,Photo.R2)
  
# convert temp to K
mydata$K<-mydata$Temp.C + 273.15

# define the Sharpe-Schoolfield equation
schoolfield_high <- function(lnc, E, Eh, Th, temp, Tc) {
  Tc <- 273.15 + Tc
  k <- 8.62e-5
  boltzmann.term <- lnc + log(exp(E/k*(1/Tc - 1/temp))) #units are eV/K, electrovolts/Kelvin
  inactivation.term <- log(1/(1 + exp(Eh/k*(1/Th - 1/temp))))
  return(boltzmann.term + inactivation.term)
}


# something weird is happening with the locations
mydata$Organism.ID<-droplevels(mydata$Organism.ID)

#filter out the NP
mydata<-mydata %>%
  filter(rate.type !='NP')%>%
  droplevels()

#log rate
mydata$log.rate<-log(mydata$umol.cm2.hr+1)
#mydata$Location<-as.character(mydata$Location)
#mydata$Organism.ID<-as.character(mydata$Organism.ID)
mydata$Phylum<-as.character(mydata$Phylum)
mydata$Habitat<-as.character(mydata$Habitat)
mydata$CommonName<-as.character(mydata$CommonName)
mydata$Species.Fullname<-as.character(mydata$Species.Fullname)
mydata$Species<-as.character(mydata$Species)

fits <- mydata %>%
  group_by(Fragment.ID, rate.type, Species) %>%
  nest() %>%
  mutate(fit = purrr::map(data, ~ nls_multstart(log.rate ~ schoolfield_high(lnc, E, Eh, Th, temp = K, Tc = 26),
                                                data = .x,
                                                iter = 1000,
                                                start_lower = c(lnc = -10, E = 0.1, Eh = 0.2, Th = 285),
                                                start_upper = c(lnc = 10, E = 2, Eh = 5, Th = 330),
                                                supp_errors = 'Y',
                                                na.action = na.omit,
                                                lower = c(lnc = -10, E = 0, Eh = 0, Th = 0))))



#broom, models over and over again and purr for lists
#make dplyr code and export the slope, intercept and R2 values
#get r2, extract predit and pull out and join lists, shows o vs predcited and r2
PredictedAcli_1 <- predict(fits$fit[[1]])
ObservedAcli_1 <- mydata$umol.cm2.hr[mydata$Fragment.ID == "Acali_1_R"]

po <- lm(PredictedAcli_1 ~ ObservedAcli_1)
summary(po)
plot(ObservedAcli_1,PredictedAcli_1)
abline(po)
legend("topleft", bty="n", legend=paste("r2 =", format(summary(po)$adj.r.squared, digits=4)))



# look at a single fit
summary(fits$fit[[1]])

# look at output object
select(fits, Fragment.ID, data, fit)  

# get summary info
info <- fits %>%
  unnest(fit %>% map(glance))

# get params
params <- fits %>%
  unnest(fit %>% map(tidy))

# get confidence intervals
CI <- fits %>% 
  unnest(fit %>% map(~ confint2(.x) %>%
                       data.frame() #%>%
                     #     rename(., conf.low = X2.5.., conf.high = X97.5..)
  )) %>%
  group_by(., Fragment.ID) %>%
  mutate(., term = c('lnc', 'E', 'Eh', 'Th')) %>%
  ungroup()

colnames(CI)[3:4]<-c("conf.low", "conf.high") # rename columns

# merge parameters and CI estimates
params <- merge(params, CI, by = intersect(names(params), names(CI)))

# get predictions
preds <- fits %>%
  unnest(fit %>% map(augment))

#select(info, fragment.ID, logLik, AIC, BIC, deviance, df.residual)

# new data frame of predictions, do this to set a sequence to make a smooth curve with your prediction points
new_preds <- mydata %>%  do(., data.frame(K = seq(min(.$K), max(.$K), length.out = 150), stringsAsFactors = FALSE)) #setting a specific sequence so you can have a smooth curve

# max and min for each curve
max_min <- mydata %>% group_by(Fragment.ID) %>%
  dplyr::summarise(., min_K = min(K), max_K = max(K)) %>%
  ungroup()

# create new predictions
preds2 <- fits %>%
  unnest(fit %>% map(augment, newdata = new_preds)) %>%
  merge(., max_min, by = "Fragment.ID") %>%
  group_by(., Fragment.ID) %>%
  filter(., K > unique(min_K) & K < unique(max_K)) %>%
  dplyr::rename(., ln.rate = .fitted) %>%
  ungroup()

#want to do ggplot where we look at fragments individually
ggplot() +
  geom_point(aes(K - 273.15, log.rate, col = Species), size = 2, mydata) +
  geom_line(aes(K - 273.15, ln.rate, col = Species, group = Fragment.ID), alpha = 0.5, preds2) +
  facet_wrap_paginate(~ Fragment.ID, labeller = labeller(.multi_line = FALSE), nrow=5, ncol =5,  page = 9, scales = 'free') + # put the plots on multiple pages
  #facet_wrap(~ Fragment.ID, labeller = labeller(.multi_line = FALSE)) +
  #scale_colour_manual(values = c('green4', 'blue', 'black')) +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  ylab(expression("Log " *mu*"mol (" *cm^-2 *hr^-1*")")) +
  xlab('Temperature (C)') +
  #geom_hline(yintercept=0, color = "red") +
  theme(legend.position = c(0.91, 0.85))+
  labs(color = "Rate Type")  +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "none") 

# plot all P, R and C values in TPCs
ggplot() +
  geom_point(aes(K - 273.15, log.rate, col =  rate.type, group = Species), size = 2, mydata) +
  geom_line(aes(K - 273.15, ln.rate, col = rate.type, group = Fragment.ID), alpha = 0.5, preds2) +
  facet_wrap(~ Species, labeller = labeller(.multi_line = FALSE), scales = 'free') +
  #scale_colour_manual(values = c('green4', 'blue', 'black')) +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  ylab(expression("Log " *mu*"mol (" *cm^-2 *hr^-1*")")) +
  xlab('Temperature (?C)') +
  #geom_hline(yintercept=0, color = "red") +
  #theme(legend.position = c(0.91, 0.85))+
 # theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave(filename = "Output/TPCcurvesGalapagos.png", device = "png", width = 10, height = 8)

# function for calculating Topt

get_topt <- function(E, Th, Eh){
  return((Eh*Th)/(Eh + (8.62e-05 *Th*log((Eh/E) - 1))))
}


# calc topts for all 
Topt_data <- params %>%
  select(Fragment.ID, term, estimate,rate.type,Species) %>%
  spread(term, estimate) %>%
  mutate(Topt = get_topt(E, Th, Eh)) %>%
  group_by(., Species, Fragment.ID, rate.type)

#get temerature back in celcius not K
Topt_data$Topt <- Topt_data$Topt - 273.15 


data.summary<-Topt_data %>%
  group_by(Species, rate.type) %>% #tells to group by these two factors
  #dplyr::summarise_if(is.numeric, list(mean = mean, se=sd(.)/sqrt(n())), na.rm=TRUE) #calculates mean and s.e.
  dplyr::summarise_if(is.numeric, mean, na.rm=TRUE) #calculates mean and s.e.

data.summary


ggplot(data.summary, aes(x=reorder(Species, -mean), y=mean, group = rate.type)) +
  theme_bw()+
  theme(legend.title=element_text(colour="black", size=12), axis.text.x=element_text(face="bold", color="black", size=12), axis.text.y=element_text(face="bold", color="black", size=13), axis.title.x = element_text(color="black", size=18, face="bold"), axis.title.y = element_text(color="black", size=18, face="bold"),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_point(position="dodge", size=2) +
  scale_y_continuous(expand = c(0,0.2)) + scale_x_discrete(expand = c(0.1,0.1)) + #adjust space around graph to left and right of discrtet values (GP and C)
  theme(legend.text=element_text(size=rel(1))) + #makes legend elements larger
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), position=position_dodge(width=0.9), width=0.1) +
  labs(x="Species", y="Thermal Optimum (°C)") +
  facet_wrap(~rate.type, scales = 'free')+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(filename = "Output/Topt_graph.png", device = "png", width = 10, height = 5)


# Join the parameter data with the species info
Topt_data<-left_join(Topt_data, Sp.info)

Topt_data.clean<-Topt_data[complete.cases(Topt_data),] # remove all the NAs

# pull out only the Respiration data
Topt_data.clean.R<- Topt_data.clean %>%
  filter(rate.type=='R', Eh<20, Functional.group!='Producer')

pca<-prcomp(sqrt(sqrt(Topt_data.clean.R[,c(4:8)])), center = TRUE, scale. = TRUE)

# bring together the PCA with the species info
pca.all<-cbind(pca$x, Topt_data.clean.R[, c(3,9:12)])

# get the loadings
df_out_r <- as.data.frame(pca$rotation)
df_out_r$feature <- row.names(df_out_r)

#pca.all<-pca.all[-18,]
# PCA by habitat
habplot<-ggplot(pca.all, color =Habitat)+
  geom_point(aes(x = PC1, y = PC2, color = Habitat))+
  stat_ellipse(aes(x = PC1, y = PC2, color = Habitat))+
  geom_segment(data = df_out_r, aes(x = 0, y = 0,xend=PC1*3, yend=PC2*3))+
  geom_text(data = df_out_r, aes(x=PC1*3, y=PC2*3,label = feature))+
  theme_bw()+
  xlim(-8,8)+
  ylim(-4,4)+
  ggtitle('Habitat')

# PCA by functional group
Funplot<- ggplot(pca.all, color =Functional.group)+
  geom_point(aes(x = PC1, y = PC2, color = Functional.group))+
  stat_ellipse(aes(x = PC1, y = PC2, color = Functional.group))+
  geom_segment(data = df_out_r, aes(x = 0, y = 0,xend=PC1*3, yend=PC2*3))+
  geom_text(data = df_out_r, aes(x=PC1*3, y=PC2*3,label = feature))+
  theme_bw()+
  xlim(-8,8)+
  ylim(-4,4)+
   ggtitle('Functional Group')

 # PCA by phyla
Phylplot<-ggplot(pca.all, color =Phylum)+
   geom_point(aes(x = PC1, y = PC2, color = Phylum))+
   stat_ellipse(aes(x = PC1, y = PC2, color = Phylum))+
   geom_segment(data = df_out_r, aes(x = 0, y = 0,xend=PC1*3, yend=PC2*3))+
   geom_text(data = df_out_r, aes(x=PC1*3, y=PC2*3,label = feature))+
   theme_bw()+
  xlim(-8,8)+
  ylim(-4,4)+
   ggtitle('Phylum')

pcaplots<-grid.arrange(habplot, Funplot, Phylplot)
ggsave(filename = 'Output/pcaplots.png', , plot = pcaplots, device = 'png', width = 5, height = 10, dpi = 250) 
