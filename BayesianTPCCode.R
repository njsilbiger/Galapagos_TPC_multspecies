### TPC code for Galapagos inverts collect by JB and MB
### Data analysis by Nyssa Silbiger
### Last edited on 8/2/2019
#########################################################



rm(list=ls())

##Load packages ##########

library(broom)
library(purrr)
library(tidyverse)
library(ggforce)
library(gridExtra)
library(brms)
library(tidybayes)
library(modelr)
library(bayesplot)
library(cowplot)


#load data ########################

photo.data <- read.csv("Data/GalapagosRates.csv")
photo.data$X <- NULL
#View(photo.data)
#glimpse(photo.data)

# species list
Sp.info<-read.csv('Data/SpeciesList.csv')

# list of files to remove 
BadData<-read.csv('Data/quality_control_inverts.csv')

## remove the bad data after QC #################
badrows<-which(photo.data$ID %in% BadData$ID)

photo.data<-photo.data[-badrows,]

# remove the NAs from the data
#photo.data<-photo.data[-which(is.na(photo.data$umol.cm2.hr)),]


# Clean Data for Analysis #################
mydata <- photo.data
## clean some of the bad data
bad<-which(mydata$umol.cm2.hr< 0)
mydata<-mydata[-bad,]

## might need to add this back
 bad1<-c( which(mydata$Species=='Bgran' & mydata$umol.cm2.hr>70), 
         which(mydata$Species=='Secu' & mydata$umol.cm2.hr>700), which(mydata$Species=='Ltube' & mydata$umol.cm2.hr>31),
         which(mydata$Species=='Tcoc' & mydata$umol.cm2.hr<5),which(mydata$Species=='Egala2' & mydata$umol.cm2.hr<6 & mydata$Temp.C>30),
         which(mydata$Species=='Cfusc' & mydata$umol.cm2.hr>400)
         )
 mydata<-mydata[-bad1,]

ggplot(mydata) +
  geom_point(aes(x = Temp.C, y = umol.cm2.hr, group =Species ))+
  facet_wrap(~Species, scales = "free_y")


#join with the mydata
mydata<-left_join(mydata,Sp.info)

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

#### Run Bayesian analysis ##############

# Make a function that runs the analysis 

BayeTPC<-function(mydata = mydata, SpeciesName,  PlotDiagnostics = TRUE, PlotResults = TRUE){
Species.selected<-mydata[mydata$Species==SpeciesName,]
spec.name<-unique(Species.selected$Species.Fullname)
y<-Species.selected$log.rate

options("scipen"=100,digits=12) # stan doesnt like scientific notation. This fixes that

# fit the model
# create a generated quantities block for Topt and for posterior predictive checks
stanvars <- stanvar(scode = "real Topt; vector[N] y_new;
vector[N] nlp_lnc = X_lnc * b_lnc;
  vector[N] nlp_Eh = X_Eh * b_Eh;
  vector[N] nlp_Th = X_Th * b_Th;
  vector[N] nlp_E = X_E * b_E;
Topt =  (b_Eh[1] * b_Th[1]) / (b_Eh[1] + (0.0000862 * b_Th[1] * log((b_Eh[1] / b_E[1]) - 1)));
for (n in 1:N) {  
nlp_lnc[n] += r_1_lnc_1[J_1[n]] * Z_1_lnc_1[n];
 y_new[n] = student_t_rng(nu, nlp_lnc[n] + log(exp(nlp_E[n] / 0.0000862 * (1 / 299.15 - 1 / C_1[n]))) + log(1 / (1 + exp(nlp_Eh[n] / 0.0000862 * (1 / nlp_Th[n] - 1 / C_1[n])))), sigma);};",
                    block = "genquant")


fit1<-brm(
  bf(log.rate ~ lnc + log(exp(E/.0000862*(1/299.15 - 1/K))) +log(1/(1 + exp(Eh/.0000862*(1/Th - 1/K)))), 
            lnc ~ 1 + (1|Organism.ID), Eh~1, Th~1 , E~1 ,  nl = TRUE ),
            data = Species.selected,
            family = student(),
            prior = c(
              prior(normal(0,10), nlpar = "E", lb = 0),  # set the priors
              prior(normal(0,10), nlpar = "Eh", lb = 0),
              prior(normal(320, 10), nlpar = "Th", lb = 0),
              prior(normal(0, 10), nlpar = "lnc", lb = 0)
            ), control = list(adapt_delta = 0.99, max_treedepth = 20), # force stan to take smaller steps to reduce divergent errors
  cores = 4, chains = 4, seed = 126, iter = 3000, warmup = 2000,stanvars = stanvars) 


if (PlotDiagnostics==TRUE){
## assess the fits
posterior <-  posterior_samples(fit1)
#str(posterior)
post <- as.data.frame(fit1)
post %>% select( nu, sigma,b_lnc_Intercept, b_Eh_Intercept,b_Th_Intercept,b_E_Intercept) %>% cor() # correlation matrix
(coef_mean <- post %>% select(nu, sigma,b_lnc_Intercept, b_Eh_Intercept,b_Th_Intercept,b_E_Intercept) %>% summarise_all(mean) %>% as.numeric) # parameter means

# convergence diagnostics -------------------------------------------------
# plotting one at a time
#plot(fit1,'b_lnc_Intercept')
# plot the traceplots
pdf(paste0('Output/diagnostics/traceplots_',spec.name,'.pdf'))
plot(fit1, ask = FALSE, newpage = FALSE) # all of them
dev.off()

## posterior predictive checks
y_rep <- as.matrix(fit1, pars = "y_new")
dim(y_rep)
check1<-ppc_dens_overlay(y, y_rep[1:200, ]) # comparing density of y with densities of y over 200 posterior draws.
check2<-ppc_stat(y = y, yrep = y_rep, stat = "mean") # compare estimates of summary statistics
check3<-ppc_scatter_avg(y = y, yrep = y_rep) # observed vs predicted with 1:1 line

ppcheckfig<-plot_grid(check1, check2, check3, labels = 'AUTO', nrow = 1 )
# add title
title <- ggdraw() + draw_label(spec.name, fontface='italic')
ppcheckf<-plot_grid(title, ppcheckfig,nrow = 2 , rel_heights=c(0.1, 1))
  
ggsave(filename = paste0('Output/diagnostics/ppchecks_',spec.name,'.pdf'),plot = ppcheckf, width = 9, height = 6 )

}

### plot the results
# plot the population level effects
#pt1<-marginal_effects(fit1)

 #population level 
pt2<-ggplot(as.data.frame(pt1$K))+ # pull pur the fitted data
  geom_line(aes(K-273.15, estimate__))+ # convert temp to celcius
  geom_ribbon(aes(K-273.15, ymin = lower__, ymax = upper__), alpha=0.3)+
  geom_point(data = Species.selected, aes(x = K-273.15, y = log.rate)) +# add the raw data
  theme_bw()+
  xlab(expression(paste('Temperature (',~degree,'C)')))+
  ylab(expression(paste('Log respiration rate (', mu, "mol cm"^-2, 'hr'^-1,")")))+
  theme(text = element_text(size=18), title = element_text(face="italic"))
  

# plot the individual effects
# conditions <- data.frame(Organism.ID = unique(Species.selected$Organism.ID))
# rownames(conditions) <- unique(Species.selected$Organism.ID)
# me_loss <- marginal_effects(
#   fit1, conditions = conditions, 
#   re_formula = NULL, method = "predict"
# )
#plot(me_loss, ncol = 5, points = TRUE)

# individual level
# pt3<-ggplot(as.data.frame(me_loss$K))+ # pull pur the fitted data
#   geom_line(aes(K-273.15, estimate__, group = cond__, color = cond__, lwd = 1))+ # convert temp to celcius
#   geom_point(data = Species.selected, aes(x = K-273.15, y = log.rate)) +# add the raw data
#   theme_bw()+
#   scale_color_brewer(palette = "Set2")+
#   xlab(expression(paste('Temperature (',~degree,'C)')))+
#   ylab(expression(paste('Log respiration rate (', mu, "mol cm"^-2, 'hr'^-1,")")))+
#   theme(text = element_text(size=18),legend.position = "none", title = element_text(face="italic"))+
#   ggtitle(spec.name)

# plot the individual effects
# conditions <- data.frame(Species = unique(Acali$Species))
# rownames(conditions) <- unique(Acali$Species)
# me_loss <- marginal_effects(
#   fit1, conditions = conditions, 
#   re_formula = NULL, method = "predict"
# )
# plot(me_loss, ncol = 2, points = TRUE)

# plot all the fits on the same plot with error
pt4<-Species.selected %>%
  group_by(Organism.ID) %>%
  data_grid(K = seq(min(K),max(K)+10, by = 0.1)) %>% # make the predictions go to 10 degrees plus 
  #data_grid(K = seq_range(K+10, n = 150)) %>%
  add_fitted_draws(fit1) %>%
  ggplot(aes(x = K-273.15, y = log.rate, color = ordered(Organism.ID))) +
  stat_lineribbon(aes(y = .value), .width = c(0.95)) + # 95%CI
  geom_point(data = Species.selected) +
  scale_fill_brewer(palette = "Greys") +
  scale_color_brewer(palette = "Set2")+
  theme_bw()+ xlab(expression(paste('Temperature (',~degree,'C)')))+
  ylab(expression(paste('Log respiration rate (', mu, "mol cm"^-2, 'hr'^-1,")")))+
  theme(text = element_text(size=18),legend.position = "none", title = element_text(face="italic"))
  #ggtitle(spec.name)

TPC_plots1<-plot_grid(pt2, pt4, labels = 'AUTO')
title <- ggdraw() + draw_label(spec.name, fontface='italic')
TPC_plots<-plot_grid(title, TPC_plots1, nrow = 2 , rel_heights=c(0.1, 1))
ggsave(filename = paste0('Output/TPC_plots/TPC_',spec.name,'.pdf'), plot = TPC_plots, width = 8, height = 5)

# calculate 95%Ci
params<-fit1 %>%
  gather_draws( b_E_Intercept, b_Eh_Intercept, b_lnc_Intercept, b_Th_Intercept, sd_Organism.ID__lnc_Intercept, sigma, Topt) %>%
  median_qi()

return(params)
}

#Run Bayesian TPC for all species, export plots, and save the 95% CI of each parameter
Param.output<-mydata %>%
  filter(Species != 'Acali') %>% # do everything but Acali for now
  group_by(Species) %>%
  do(BayeTPC(mydata = ., SpeciesName = unique(.$Species)))


# find Tmax
## UNCOMMENT THIS WHEN IT STOPS CRASHING
# Tmax<-Acali %>%
#   group_by(Organism.ID, .draw) %>%
#   data_grid(K = seq(min(K),max(K)+10, by = 0.01)) %>%
#   add_predicted_draws(fit1) 

# Tmax<-
#   Tmax[which(Tmax$K>params[params$.variable=='Topt',2] && Tmax$value<1),]

# 
# fit1.data %>%
#   median_qi(Topt)

