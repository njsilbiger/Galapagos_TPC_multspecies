rm(list=ls())

##Install packages
# load packages
#library(nls.multstart)
library(broom)
library(purrr)
library(tidyverse)
#library(nlstools)
#library(nls2)
library(ggforce)
library(gridExtra)
library(brms)
library(tidybayes)
library(modelr)
library(bayesplot)


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

# something weird is happening with the locations
mydata$Organism.ID<-droplevels(mydata$Organism.ID)

#filter out the NP
mydata<-mydata %>%
  filter(rate.type !='NP')%>%
  droplevels()

#log rate
mydata$log.rate<-log(mydata$umol.cm2.hr+1)


Acali<-mydata[mydata$Species=='Lsemi',]

options("scipen"=100,digits=12) # stan doesnt like scientific notation. This fixes that

# fit the model
# create a generated quantities block for Topt
stanvars <- stanvar(scode = "real Topt;
   Topt =  (b_Eh[1] * b_Th[1]) / (b_Eh[1] + (0.0000862 * b_Th[1] * log((b_Eh[1] / b_E[1]) - 1)));", 
                    block = "genquant")
fit1<-brm(
  bf(log.rate ~ lnc + log(exp(E/.0000862*(1/299.15 - 1/K))) +log(1/(1 + exp(Eh/.0000862*(1/Th - 1/K)))), 
            lnc ~ 1 + (1|Organism.ID), Eh~1, Th~1 , E~1 ,  nl = TRUE ),
            data = Acali,
            family = student(),
            prior = c(
              prior(normal(0,10), nlpar = "E", lb = 0),  # set the priors
              prior(normal(0,10), nlpar = "Eh", lb = 0),
              prior(normal(320, 10), nlpar = "Th", lb = 0),
              prior(normal(0, 10), nlpar = "lnc", lb = 0)
            ), control = list(adapt_delta = 0.99, max_treedepth = 20), # force stan to take smaller steps to reduce divergent errors
  cores = 4, chains = 4, seed = 126, iter = 3000, warmup = 2000,stanvars = stanvars) 
 # use multiple cores for parallel computing
# pull out all the draws and calculate Topt with appropriate error
# fit1.data<-as.data.frame(fit1) 
# #Topt function
# get_topt <- function(E, Th, Eh){
#   return((Eh*Th)/(Eh + (8.62e-05 *Th*log((Eh/E) - 1))))
# }
# 
#  dataList =  standata(fit1)
# # try with stan model with Topt in generated qualtity
# fit <- stan(file = 'GalapagosTPCStan.stan', data = dataList, iter=2000, warmup=1000, cores = 3, control = list(adapt_delta = 0.99, max_treedepth = 20))
# 
# ppc_intervals(
#   y = Acali$log.rate,
#   yrep = posterior_predict(fit),
#   x = Acali$K,
#   prob = 0.5
# ) +
#   labs(
#     x = "Weight (1000 lbs)",
#     y = "MPG",
#     title = "50% posterior predictive intervals \nvs observed miles per gallon",
#     subtitle = "by vehicle weight"
#   ) +
#   panel_bg(fill = "gray95", color = NA) +
#   grid_lines(color = "white")
# 
# ### plot the within species variance by species... which ones species have more variability in tpc? 
# 
# # split the underscores for clean names and gather the data so there is a columnen for Organism ID for easier plotting
# 
# 
# fit1.data$Topt<-get_topt(E = fit1.data$b_E_Intercept, Eh = fit1.data$b_Eh_Intercept, Th = fit1.data$b_Th_Intercept)
# 

# plot the population level effects
pt1<-marginal_effects(fit1)

 #population level 
pt2<-ggplot(as.data.frame(pt1$K))+ # pull pur the fitted data
  geom_line(aes(K-273.15, estimate__))+ # convert temp to celcius
  geom_ribbon(aes(K-273.15, ymin = lower__, ymax = upper__), alpha=0.3)+
  geom_point(data = Acali, aes(x = K-273.15, y = log.rate)) +# add the raw data
  theme_bw()+
  xlab(expression(paste('Temperature (',~degree,'C)')))+
  ylab(expression(paste('Log respiration rate (', mu, "mol cm"^-2, 'hr'^-1,")")))+
  theme(text = element_text(size=18))

# plot the individual effects
conditions <- data.frame(Organism.ID = unique(Acali$Organism.ID))
rownames(conditions) <- unique(Acali$Organism.ID)
me_loss <- marginal_effects(
  fit1, conditions = conditions, 
  re_formula = NULL, method = "predict"
)
#plot(me_loss, ncol = 5, points = TRUE)

# individual level
pt3<-ggplot(as.data.frame(me_loss$K))+ # pull pur the fitted data
  geom_line(aes(K-273.15, estimate__, group = cond__, color = cond__, lwd = 1))+ # convert temp to celcius
  geom_point(data = Acali, aes(x = K-273.15, y = log.rate)) +# add the raw data
  theme_bw()+
  scale_color_brewer(palette = "Set2")+
  xlab(expression(paste('Temperature (',~degree,'C)')))+
  ylab(expression(paste('Log respiration rate (', mu, "mol cm"^-2, 'hr'^-1,")")))+
  theme(text = element_text(size=18),legend.position = "none")

# plot the individual effects
# conditions <- data.frame(Species = unique(Acali$Species))
# rownames(conditions) <- unique(Acali$Species)
# me_loss <- marginal_effects(
#   fit1, conditions = conditions, 
#   re_formula = NULL, method = "predict"
# )
# plot(me_loss, ncol = 2, points = TRUE)

# plot all the fits on the same plot with error
Acali %>%
  group_by(Organism.ID) %>%
  data_grid(K = seq(min(K),max(K)+10, by = 0.1)) %>% # make the predictions go to 10 degrees plus 
  #data_grid(K = seq_range(K+10, n = 150)) %>%
  add_fitted_draws(fit1) %>%
  ggplot(aes(x = K-273.15, y = log.rate, color = ordered(Organism.ID))) +
  stat_lineribbon(aes(y = .value), .width = c(0.95)) + # 95%CI
  geom_point(data = Acali) +
  scale_fill_brewer(palette = "Greys") +
  scale_color_brewer(palette = "Set2")+
  theme_bw()+ xlab(expression(paste('Temperature (',~degree,'C)')))+
  ylab(expression(paste('Log respiration rate (', mu, "mol cm"^-2, 'hr'^-1,")")))+
  theme(text = element_text(size=18),legend.position = "none")


# calculate 95%Ci
params<-fit1 %>%
  gather_draws( b_E_Intercept, b_Eh_Intercept, b_lnc_Intercept, b_Th_Intercept, sd_Organism.ID__lnc_Intercept, sigma, Topt) %>%
  median_qi()

# find Tmax
Tmax<-Acali %>%
  group_by(Organism.ID, .draw) %>%
  data_grid(K = seq(min(K),max(K)+10, by = 0.01)) %>%
  add_predicted_draws(fit1) 

Tmax<-
  Tmax[which(Tmax$K>params[params$.variable=='Topt',2] && Tmax$value<1),]

# 
# fit1.data %>%
#   median_qi(Topt)

plot(fit1)

#plot(density(fit1.data$Topt))

# PP check
pp_check(fit1)

