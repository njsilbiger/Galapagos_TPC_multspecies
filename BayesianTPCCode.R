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

## remove the bad data after QC ################# (i.e. citter kept touching O2 probe giving unreliable measurements)
badrows<-which(photo.data$ID %in% BadData$ID)

photo.data<-photo.data[-badrows,]

# remove the NAs from the data
#photo.data<-photo.data[-which(is.na(photo.data$umol.cm2.hr)),]


# Clean Data for Analysis #################
mydata <- photo.data
## clean some of the bad data
bad<-which(mydata$umol.cm2.hr< 0)
mydata<-mydata[-bad,]

## removing clear outliers in the data. Probably because the critter touched the probed
 bad1<-c( which(mydata$Species=='Bgran' & mydata$umol.cm2.hr>70), 
         which(mydata$Species=='Secu' & mydata$umol.cm2.hr>700), which(mydata$Species=='Ltube' & mydata$umol.cm2.hr>31),
         which(mydata$Species=='Tcoc' & mydata$umol.cm2.hr<5),which(mydata$Species=='Egala2' & mydata$umol.cm2.hr<10 & mydata$Temp.C>30 & mydata$Temp.C<40),
         which(mydata$Species=='Cfusc' & mydata$umol.cm2.hr>400)
         )
 mydata<-mydata[-bad1,]

ggplot(mydata) +
  geom_point(aes(x = Temp.C, y = umol.cm2.hr, group =Organism.ID ))+
  geom_line(aes(x = Temp.C, y = umol.cm2.hr, group =Organism.ID ))+
  facet_wrap(~Species, scales = "free_y")


# Egala was done twice, because the first time was methods testing.  Remove the first time and keep the second
mydata<-mydata %>%
  filter(Species != "Egala")

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
  dplyr::filter(Functional.group=='Producer') %>%  # only do this for the producers since they are the only ones that photosynthesize
  dplyr::group_by(Organism.ID, Temp.Cat, Species) %>% 
  dplyr::summarize(rates = sum(umol.cm2.hr), Temp.C=mean(Temp.C)) %>%
  dplyr::mutate(rate.type="GP", Light_Dark="Light") %>%  
  dplyr::rename(umol.cm2.hr=rates) %>%
  dplyr::mutate(Fragment.ID=paste0(Organism.ID, "_", rate.type)) %>%
  dplyr::left_join(.,Sp.info)

# put in the full species names for light
algaenames<-mydata %>%
  dplyr::filter(Functional.group=='Producer' & Light_Dark=='Dark') %>%
  dplyr::select(Species.Fullname, Species)%>%
  unique()

for(i in 1:length(algaenames)){
mydata[which(mydata$Species==algaenames[i,2]),'Species.Fullname']<-algaenames[i,1]
}

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

BayeTPC<-function(mydata = mydata, SpeciesName, rate, PlotDiagnostics = TRUE, PlotResults = TRUE){
#Species.selected<-mydata[mydata$Species==SpeciesName & mydata$rate.type==rate,]
Species.selected<-mydata
spec.name<-unique(Species.selected$Species.Fullname)
y<-Species.selected$log.rate
LD<-unique(Species.selected$rate.type) # light or dark for the algae

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
  cores = 4, chains = 4, seed = 126, iter = 3000, warmup = 2000,stanvars = stanvars, silent = TRUE) 



if (PlotDiagnostics==TRUE){
## assess the fits
posterior <-  posterior_samples(fit1)
#str(posterior)
#post <- as.data.frame(fit1)
#post %>% select( nu, sigma,b_lnc_Intercept, b_Eh_Intercept,b_Th_Intercept,b_E_Intercept) %>% cor() # correlation matrix
#(coef_mean <- post %>% select(nu, sigma,b_lnc_Intercept, b_Eh_Intercept,b_Th_Intercept,b_E_Intercept) %>% summarise_all(mean) %>% as.numeric) # parameter means

# convergence diagnostics -------------------------------------------------
# plotting one at a time
#plot(fit1,'b_lnc_Intercept')
# plot the traceplots
pdf(paste0('Output/diagnostics/traceplots_',spec.name,LD,'.pdf'))
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
  
ggsave(filename = paste0('Output/diagnostics/ppchecks_',spec.name,LD,'.pdf'),plot = ppcheckf, width = 9, height = 6 )

}

if (PlotResults == TRUE){
### plot the results
# plot the population level effects
pt1<-marginal_effects(fit1)

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
ggsave(filename = paste0('Output/TPC_plots/TPC_',spec.name,LD,'.pdf'), plot = TPC_plots, width = 8, height = 5)
}

# only print params if it exists
#if (exists(x = 'fit1')) {
# calculate 95%Ci
params<-fit1 %>%
  gather_draws( b_E_Intercept, b_Eh_Intercept, b_lnc_Intercept, b_Th_Intercept, sd_Organism.ID__lnc_Intercept, sigma, Topt) %>%
  median_qi()

return(params)
#} else {
#  params<-as.data.frame(x = matrix(NA,nrow = 1, ncol = 7))
#  colnames(params)<-c('.variable',' \.value','.lower', '.upper', '.width', '.point','.interval')
#  return(params)
#}

# print when each species is done
print(paste(spec.name, 'is done'))
rm(fit1) # clear the fit

}

p<-as.data.frame(x = matrix(NA,nrow = 1, ncol = 7))
colnames(p)<-c('.variable','.value','.lower', '.upper', '.width', '.point','.interval')

#Same functions, but skips over species with an error
BayeTPC_noerror <- possibly(BayeTPC, otherwise = p)

#Run Bayesian TPC for all species, export plots, and save the 95% CI of each parameter
Param.output<-mydata %>%
#  filter(Species == 'Padi') %>% # do everything but Acali for now
  group_by(Species, rate.type) %>%
  do(BayeTPC_noerror(mydata = ., SpeciesName = unique(.$Species)))

write.csv(Param.output, file = 'Data/paramoutput.csv', row.names = NULL) # write out the results
  # nest() %>%
  # mutate(params =  data %>%
  #          map(BayeTPC_noerror(mydata = data, SpeciesName = unique(.$Species), 
  #                              rate = unique(.$rate.type))))


# remove the species and rates that did not converge
not.converge<-which(is.na(Param.output$.value))
Param.output<-Param.output[-not.converge,]

# Put Topt aand Th in celcius
Param.output[Param.output$.variable=='Topt', c(4:6)] <-Param.output[Param.output$.variable=='Topt', c(4:6)]-273.15
Param.output[Param.output$.variable=='b_Th_Intercept', c(4:6)] <-Param.output[Param.output$.variable=='b_Th_Intercept', c(4:6)]-273.15

# add a column for plotting names
prettynames<-data.frame(.variable = unique(Param.output$.variable))
prettynames$varnames<-c("E (eV)","Eh","b(Tc)","Th (°C)", "Within species variance","Observation error", "Thermal Optimum (°C)")

# add species full name
sp.fullnames<- mydata %>%
  select(Species.Fullname, Species)%>%
  unique() %>%
  filter(Species.Fullname != '<NA>') # remove the NA

Sp.info<-left_join(Sp.info, sp.fullnames)

Param.output<-left_join(Param.output,prettynames) %>% # make the names easier for plotting
  left_join(.,Sp.info) %>% # join with the species info
  filter(Species != 'Egala') # remove the first Egala since this was done twice

# Plot the Topt in decending order for R colored by functional group
paramplots<-Param.output %>%
  filter(rate.type == 'R')%>%
  group_by(.variable) %>%
  nest() %>%
  mutate(plot = map2(data, .variable, ~ggplot(., aes(x=reorder(Species, -.value), y=.value, color = Functional.group)) +
  theme_bw()+
  theme(legend.title=element_text(colour="black", size=12), axis.text.x=element_text(face="bold", color="black", size=12), axis.text.y=element_text(face="bold", color="black", size=13), axis.title.x = element_text(color="black", size=18, face="bold"), axis.title.y = element_text(color="black", size=18, face="bold"),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_point(position="dodge", size=2) +
  scale_y_continuous(expand = c(0,0.2)) + scale_x_discrete(expand = c(0.1,0.1)) + #adjust space around graph to left and right of discrtet values (GP and C)
  theme(legend.text=element_text(size=rel(1))) + #makes legend elements larger
  geom_errorbar(aes(ymax=.upper, ymin=.lower), position=position_dodge(width=0.9), width=0.1) +
  labs(x="", y=unique(.$varnames)) +
  theme(legend.position="none",axis.text.x = element_text(angle = 90, hjust = 1))))
#save the results
#map2(paste0("Output/resultplots/",paramplots$.variable, "_fungroup.pdf"), paramplots$plot, ggsave)

#extract legend
legend <- get_legend(
    paramplots$plot[[1]] + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom", legend.title = element_blank())
)

# put everything together
p1<-plot_grid(paramplots$plot[[1]],paramplots$plot[[2]], paramplots$plot[[3]],
          paramplots$plot[[5]], paramplots$plot[[6]], paramplots$plot[[7]])
results_by_FunGroup<-plot_grid(p1,legend, ncol = 1, rel_heights = c(1, .1))

ggsave(filename = 'Output/resultplots/results_by_FunGroup.pdf', results_by_FunGroup, width = 10, height = 8)

## Same plot, but color by intertidal/subtidal

# Plot the Topt in decending order for R colored by functional group
paramplots2<-Param.output %>%
  filter(rate.type == 'R')%>%
  group_by(.variable) %>%
  nest() %>%
  mutate(plot = map2(data, .variable, ~ggplot(., aes(x=reorder(Species, -.value), y=.value, color = Habitat)) +
                       theme_bw()+
                       theme(legend.title=element_text(colour="black", size=12), axis.text.x=element_text(face="bold", color="black", size=12), axis.text.y=element_text(face="bold", color="black", size=13), axis.title.x = element_text(color="black", size=18, face="bold"), axis.title.y = element_text(color="black", size=18, face="bold"),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
                       geom_point(position="dodge", size=2) +
                       scale_y_continuous(expand = c(0,0.2)) + scale_x_discrete(expand = c(0.1,0.1)) + #adjust space around graph to left and right of discrtet values (GP and C)
                       theme(legend.text=element_text(size=rel(1))) + #makes legend elements larger
                       geom_errorbar(aes(ymax=.upper, ymin=.lower), position=position_dodge(width=0.9), width=0.1) +
                       labs(x="", y=unique(.$varnames)) +
                       theme(legend.position="none",axis.text.x = element_text(angle = 90, hjust = 1))))

p2<-plot_grid(paramplots2$plot[[1]],paramplots2$plot[[2]], paramplots2$plot[[3]],
              paramplots2$plot[[5]], paramplots2$plot[[6]], paramplots2$plot[[7]])
#extract legend
legend <- get_legend(
  paramplots2$plot[[1]] + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom", legend.title = element_blank())
)
results_by_habitat<-plot_grid(p2,legend, ncol = 1, rel_heights = c(1, .1))

ggsave(filename = 'Output/resultplots/results_by_habitat.pdf', results_by_habitat, width = 10, height = 8)

# make a PCA
# spread data

Param.spread<-Param.output %>% 
  select(-c(.lower,.upper,.width,.point,.interval,varnames))%>%
  filter(rate.type=="R") %>%
  spread(.variable, .value)

pca<-prcomp(Param.spread[,c(8:10,14)], center = TRUE, scale. = TRUE)

# bring together the PCA with the species info
pca.all<-cbind(pca$x, Param.spread[, c(1:7)])

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
  xlim(-6,6)+
  ylim(-3,3)+
  ggtitle('Habitat')


fungplot<-ggplot(pca.all, color =Functional.group)+
  geom_point(aes(x = PC1, y = PC2, color = Functional.group))+
 # stat_ellipse(aes(x = PC1, y = PC2, color = Functional.group))+
  geom_segment(data = df_out_r, aes(x = 0, y = 0,xend=PC1*3, yend=PC2*3))+
  geom_text(data = df_out_r, aes(x=PC1*3, y=PC2*3,label = feature))+
  theme_bw()+
  xlim(-6,6)+
  ylim(-3,3)+
  ggtitle('Functional group')


phylagplot<-ggplot(pca.all, color =Phylum)+
  geom_point(aes(x = PC1, y = PC2, color = Phylum))+
  # stat_ellipse(aes(x = PC1, y = PC2, color = Functional.group))+
  geom_segment(data = df_out_r, aes(x = 0, y = 0,xend=PC1*3, yend=PC2*3))+
  geom_text(data = df_out_r, aes(x=PC1*3, y=PC2*3,label = feature))+
  theme_bw()+
  xlim(-6,6)+
  ylim(-3,3)+
  ggtitle('Phylum')

pcaplot<-plot_grid(habplot, fungplot, phylagplot, nrow =1)
ggsave(filename = 'Output/resultplots/pcaplots.pdf', plot = pcaplot, width = 12, height = 3)


### bring in functional data
fundata<-read.csv('Data/Inverts_List_Physio_12_Jan.csv')
fundata$X<-NULL

View(left_join(Param.output, fundata))
# make a dendrogram https://jcoliver.github.io/learn-r/008-ggplot-dendrograms-and-heatmaps.html
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

