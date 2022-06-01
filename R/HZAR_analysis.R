library(quantreg)
library(MCMCpack)
library(hzar)
library(geodist)
library(geosphere)
library(tidyverse)

## The following code follows Generate_HZAR_input_file.R script and runs HZAR (tests a number of cline models) 
## from the admixture index results. Much of the code to run HZAR is adapted from the Manakin tutorial 
## (Supp. Dataset 1) from Derryberry et al. (2014) MER 14:652-663. This also builds Fig. S8.

##    FILES REQUIRED:
##          input_hzar_avg object from the Generate_HZAR_input_file.R script [TO RUN HZAR]
##          pure_allele_dict object from Admixture_index_analysis.R script [TO BUILD FIGURE S8]
##          ../../4_Data_visualization/data_files_input_into_scripts/admixture_data.txt [TO BUILD FIGURE S8]
##          ../3_Admixture_index/admix_index_input_files/seq_data.txt [TO BUILD FIGURE S8]

##    STRUCTURE OF CODE:
##              (1) Make empty lists and input observed data
##              (2) Build helper function
##              (3) Compile and prepare for fitting
##              (4) Build initial chains for models
##              (5) Randomize initial value for each fit
##              (6) Run a chain of 3 runs for every fit request
##              (7) Check for model convergence
##              (8) Aggregate data for analysis
##              (9) Import data files and tidy them
##              (10) Calculate distances for individual samples (for plotting)
##              (11) Build plot (Fig. S8)



# (1) Make lists and input observed data -----------------------------------------------

# For average frequency data (freq$avg)
if(length(apropos("^freq$",ignore.case=FALSE)) == 0 ||
   !is.list(freq) ) freq <- list()

freq$avg <- list()
freq$avg$obs <- list()
freq$avg$models <- list()
freq$avg$fitRs <- list()
freq$avg$runs <- list()
freq$avg$analysis <- list()

# Create observed data points for average:
freq$avg$obs <-
  hzar.doMolecularData1DPops(input_hzar_avg$distance,
                             input_hzar_avg$avg_pop_sys,
                             input_hzar_avg$no_samples)



# (2) Build helper function ---------------------------------------------------

# We want to run multiple models on the same dataset (in this case, the avg freqs)
# Scaling: none (fixed min 0 and max 1); fixed (min and max fitted to obs min and max); free (min and max free params)
# Tails: none (no tails); right / left (just one tail); mirror; both (two tails independent)
freq.loadavgmodel <- function(scaling, tails, id=paste(scaling, tails, sep="."))
  freq$avg$models[[id]] <<- hzar.makeCline1DFreq(freq$avg$obs, scaling, tails)

# freq.loadavgmodel("free", "none", "modelI") # model 1
# freq.loadavgmodel("free" , "left", "modelI") # model 2
# freq.loadavgmodel("free" , "right", "modelI") # model 3
# freq.loadavgmodel("free" , "mirror", "modelI") # model 4
# freq.loadavgmodel("free" , "both", "modelI") # model 5
freq.loadavgmodel("fixed", "none", "modelI") # model 6; this is the model with the best support

# Modify models to look at distance we observed data within
# Data gathered from 0 to 641km; set to -30 and 660; this will adjust the cline center param
freq$avg$models <- sapply(freq$avg$models,
                          hzar.model.addBoxReq,
                          -30 , 660,
                          simplify=FALSE)

# Check that everything looks ok
print(freq$avg$models)



# (3) Compile and prepare for fitting -----------------------------------------

freq$avg$fitRs$init <- sapply(freq$avg$models,
                              hzar.first.fitRequest.old.ML,
                              obsData=freq$avg$obs,
                              verbose=FALSE,
                              simplify=FALSE)

# To speed up computation:
if(require(doMC)){
  registerDoMC()
} else {
  registerDoSEQ();
}

chainLength=1e5

# Each model will run off separate seeds:
# sample(1:1000, 6, replace=FALSE)
mainSeed = list(A=c(596,528,124,978,544,99),
                B=c(528,124,978,544,99,596),
                C=c(124,978,544,99,596,528),
                D=c(177,712,660,642,395,488),
                E=c(2,769,147,177,841,295),
                F=c(165, 272, 622, 738, 213, 656))


freq$avg$fitRs$init$modelI$mcmcParam$chainLength <- chainLength #1e5
freq$avg$fitRs$init$modelI$mcmcParam$burnin <- chainLength %/% 10 #1e4
freq$avg$fitRs$init$modelI$mcmcParam$seed[[1]] <- mainSeed$F # CHANGE THIS EACH TIME

# Check that the above worked
print(freq$avg$fitRs$init)



# (4) Build initial chains for models -----------------------------------------

freq$avg$runs$init <- list()

freq$avg$runs$init$modelI <-
  hzar.doFit(freq$avg$fitRs$init$modelI)
# Plot the trace:
# plot(hzar.mcmc.bindLL(freq$avg$runs$init$modelI))

# Compile a new set of fit requests based on initial chains
freq$avg$fitRs$chains <- lapply(freq$avg$runs$init,
                                hzar.next.fitRequest)

# Replicate each fit request 3 times, keeping original seeds
# while switching to a new seed channel
freq$avg$fitRs$chains <- hzar.multiFitRequest(freq$avg$fitRs$chains,
                                              each=3,
                                              baseSeed=NULL)


# (5) Randomize initial value for each fit ---------------------------------------

####################### SET CENTER VALUES ########################
## runif(15,-30,660) center for models 1-6
# model 1:
# freq$avg$fitRs$chains[[1]]$modelParam$init["center"]= 432.006881
# freq$avg$fitRs$chains[[2]]$modelParam$init["center"]= 100.356347
# freq$avg$fitRs$chains[[3]]$modelParam$init["center"]= 587.144966

# model 2:
# freq$avg$fitRs$chains[[1]]$modelParam$init["center"]= 458.839561
# freq$avg$fitRs$chains[[2]]$modelParam$init["center"]= 588.279165
# freq$avg$fitRs$chains[[3]]$modelParam$init["center"]= 28.742727

# model 3:
# freq$avg$fitRs$chains[[1]]$modelParam$init["center"]= 158.790457
# freq$avg$fitRs$chains[[2]]$modelParam$init["center"]= 649.814984
# freq$avg$fitRs$chains[[3]]$modelParam$init["center"]= 203.995312

# model 4:
# freq$avg$fitRs$chains[[1]]$modelParam$init["center"]= 189.459069
# freq$avg$fitRs$chains[[2]]$modelParam$init["center"]= 320.329683
# freq$avg$fitRs$chains[[3]]$modelParam$init["center"]= 249.764026

# model 5:
# freq$avg$fitRs$chains[[1]]$modelParam$init["center"]= 452.036011
# freq$avg$fitRs$chains[[2]]$modelParam$init["center"]= 178.492941
# freq$avg$fitRs$chains[[3]]$modelParam$init["center"]= 6.501785

# model 6: 
freq$avg$fitRs$chains[[1]]$modelParam$init["center"]= 518.74816
freq$avg$fitRs$chains[[2]]$modelParam$init["center"]= 108.20439
freq$avg$fitRs$chains[[3]]$modelParam$init["center"]= 616.63816

##################### SET WIDTH VALUES #######################
## runif(15,0,690) width for models 1-6
# model 1: 
# freq$avg$fitRs$chains[[1]]$modelParam$init["width"]= 283.09967
# freq$avg$fitRs$chains[[2]]$modelParam$init["width"]= 436.41024
# freq$avg$fitRs$chains[[3]]$modelParam$init["width"]= 84.82553

# model 2:
# freq$avg$fitRs$chains[[1]]$modelParam$init["width"]= 8.662236
# freq$avg$fitRs$chains[[2]]$modelParam$init["width"]= 516.097896
# freq$avg$fitRs$chains[[3]]$modelParam$init["width"]= 171.830732

# model 3:
# freq$avg$fitRs$chains[[1]]$modelParam$init["width"]= 93.505144
# freq$avg$fitRs$chains[[2]]$modelParam$init["width"]= 497.208883
# freq$avg$fitRs$chains[[3]]$modelParam$init["width"]= 485.961920

# model 4:
# freq$avg$fitRs$chains[[1]]$modelParam$init["width"]= 434.640996
# freq$avg$fitRs$chains[[2]]$modelParam$init["width"]= 428.257560
# freq$avg$fitRs$chains[[3]]$modelParam$init["width"]= 639.039660

# model 5:
# freq$avg$fitRs$chains[[1]]$modelParam$init["width"]= 182.189841
# freq$avg$fitRs$chains[[2]]$modelParam$init["width"]= 218.605990
# freq$avg$fitRs$chains[[3]]$modelParam$init["width"]= 20.161982

# model 6:
freq$avg$fitRs$chains[[1]]$modelParam$init["width"]= 256.992895
freq$avg$fitRs$chains[[2]]$modelParam$init["width"]= 303.897667
freq$avg$fitRs$chains[[3]]$modelParam$init["width"]= 320.772177

####################### SET pMIN VALUES ########################
## runif(15,0,1) pMin
# freq$avg$fitRs$chains[[1]]$modelParam$init["pMin"]= 0.63840100
# freq$avg$fitRs$chains[[2]]$modelParam$init["pMin"]= 0.56778522
# freq$avg$fitRs$chains[[3]]$modelParam$init["pMin"]= 0.15172399

# model 2: free // left
# freq$avg$fitRs$chains[[1]]$modelParam$init["pMin"]= 0.91416327
# freq$avg$fitRs$chains[[2]]$modelParam$init["pMin"]= 0.66370784
# freq$avg$fitRs$chains[[3]]$modelParam$init["pMin"]= 0.39050743

# model 3: free // right
# freq$avg$fitRs$chains[[1]]$modelParam$init["pMin"]= 0.40052113
# freq$avg$fitRs$chains[[2]]$modelParam$init["pMin"]= 0.25350542
# freq$avg$fitRs$chains[[3]]$modelParam$init["pMin"]= 0.37046092

# model 4: free // mirror
# freq$avg$fitRs$chains[[1]]$modelParam$init["pMin"]= 0.92339812
# freq$avg$fitRs$chains[[2]]$modelParam$init["pMin"]= 0.02262226
# freq$avg$fitRs$chains[[3]]$modelParam$init["pMin"]= 0.32464770

# model 5: free // both
# freq$avg$fitRs$chains[[1]]$modelParam$init["pMin"]= 0.23246689
# freq$avg$fitRs$chains[[2]]$modelParam$init["pMin"]= 0.69235535
# freq$avg$fitRs$chains[[3]]$modelParam$init["pMin"]= 0.24519266

####################### SET pMAX VALUES ########################
## runif(15,0,1) pMax for models 1-5
# freq$avg$fitRs$chains[[1]]$modelParam$init["pMax"]= 0.8207061
# freq$avg$fitRs$chains[[2]]$modelParam$init["pMax"]= 0.7225998
# freq$avg$fitRs$chains[[3]]$modelParam$init["pMax"]= 0.9503700

# model 2:
# freq$avg$fitRs$chains[[1]]$modelParam$init["pMax"]= 0.3551041
# freq$avg$fitRs$chains[[2]]$modelParam$init["pMax"]= 0.9698223
# freq$avg$fitRs$chains[[3]]$modelParam$init["pMax"]= 0.7361401

# model 3:
# freq$avg$fitRs$chains[[1]]$modelParam$init["pMax"]= 0.5489139
# freq$avg$fitRs$chains[[2]]$modelParam$init["pMax"]= 0.1495655
# freq$avg$fitRs$chains[[3]]$modelParam$init["pMax"]= 0.8675297

# model 4:
# freq$avg$fitRs$chains[[1]]$modelParam$init["pMax"]= 0.8374373
# freq$avg$fitRs$chains[[2]]$modelParam$init["pMax"]= 0.8774964
# freq$avg$fitRs$chains[[3]]$modelParam$init["pMax"]= 0.2848403

# model 5:
# freq$avg$fitRs$chains[[1]]$modelParam$init["pMax"]= 0.5586895
# freq$avg$fitRs$chains[[2]]$modelParam$init["pMax"]= 0.4446208
# freq$avg$fitRs$chains[[3]]$modelParam$init["pMax"]= 0.4682641

####################### SET deltaL VALUES ########################
## runif(3,0,690) deltaL: models 2, 5
# freq$avg$fitRs$chains[[1]]$modelParam$init["deltaL"]= 290.2459
# freq$avg$fitRs$chains[[2]]$modelParam$init["deltaL"]= 242.7785
# freq$avg$fitRs$chains[[3]]$modelParam$init["deltaL"]= 351.2703

# freq$avg$fitRs$chains[[1]]$modelParam$init["deltaL"]= 152.6738
# freq$avg$fitRs$chains[[2]]$modelParam$init["deltaL"]= 343.9682
# freq$avg$fitRs$chains[[3]]$modelParam$init["deltaL"]= 671.5616

####################### SET tauL VALUES ########################
## runif(3,0,1) tauL for models 2, 5
# freq$avg$fitRs$chains[[1]]$modelParam$init["tauL"]= 0.3205238
# freq$avg$fitRs$chains[[2]]$modelParam$init["tauL"]= 0.9736836
# freq$avg$fitRs$chains[[3]]$modelParam$init["tauL"]= 0.6674259

# freq$avg$fitRs$chains[[1]]$modelParam$init["tauL"]= 0.1867361
# freq$avg$fitRs$chains[[2]]$modelParam$init["tauL"]= 0.2796495
# freq$avg$fitRs$chains[[3]]$modelParam$init["tauL"]= 0.3601411

####################### SET deltaR VALUES ########################
## runif(3,0,690) deltaR for models 3, 5
# freq$avg$fitRs$chains[[1]]$modelParam$init["deltaR"]=472.6632
# freq$avg$fitRs$chains[[2]]$modelParam$init["deltaR"]=390.6439
# freq$avg$fitRs$chains[[3]]$modelParam$init["deltaR"]=545.4608

# freq$avg$fitRs$chains[[1]]$modelParam$init["deltaR"]=55.74919
# freq$avg$fitRs$chains[[2]]$modelParam$init["deltaR"]=315.64708
# freq$avg$fitRs$chains[[3]]$modelParam$init["deltaR"]=330.91392

####################### SET tauR VALUES ########################
## runif(3,0,1) tauR for models 3, 5
# freq$avg$fitRs$chains[[1]]$modelParam$init["tauR"]=0.9146836
# freq$avg$fitRs$chains[[2]]$modelParam$init["tauR"]=0.2064311
# freq$avg$fitRs$chains[[3]]$modelParam$init["tauR"]= 0.2435238

# freq$avg$fitRs$chains[[1]]$modelParam$init["tauR"]=0.6868013
# freq$avg$fitRs$chains[[2]]$modelParam$init["tauR"]=0.2202784
# freq$avg$fitRs$chains[[3]]$modelParam$init["tauR"]= 0.8867558

####################### SET tauM VALUES ########################
## runif(3,0,1) tauM for model 4
# freq$avg$fitRs$chains[[1]]$modelParam$init["tauM"]=0.6743654
# freq$avg$fitRs$chains[[2]]$modelParam$init["tauM"]=0.2145613
# freq$avg$fitRs$chains[[3]]$modelParam$init["tauM"]= 0.3934755

####################### SET deltaM VALUES ######################
# ## runif(3,0,690) deltaM for model 4
# freq$avg$fitRs$chains[[1]]$modelParam$init["deltaM"]=279.9647
# freq$avg$fitRs$chains[[2]]$modelParam$init["deltaM"]=481.7391
# freq$avg$fitRs$chains[[3]]$modelParam$init["deltaM"]=269.2668



# (6) Run a chain of 3 runs for every fit request -----------------------------

freq$avg$runs$chains <- hzar.doChain.multi(freq$avg$fitRs$chains,
                                           doPar=TRUE,
                                           inOrder = FALSE,
                                           count=3)


# (7) Check for model convergence ---------------------------------------------

summary(do.call(mcmc.list,
                lapply(freq$avg$runs$chains[1:3], # model I: chains 1-3; model II: chains 4-6, etc.
                       function(x)hzar.mcmc.bindLL(x[[3]]))))


# (8) Aggregate data for analysis -------------------------------

freq$avg$analysis$initDGs <- list()
freq$avg$analysis$initDGs$modelI <- hzar.dataGroup.add(freq$avg$runs$init$modelI)
freq$avg$analysis$oDG <- hzar.make.obsDataGroup(freq$avg$analysis$initDGs)
freq$avg$analysis$oDG <- hzar.copyModelLabels(freq$avg$analysis$initDGs,
                                              freq$avg$analysis$oDG)
freq$avg$analysis$oDG <- hzar.make.obsDataGroup(lapply(freq$avg$runs$chains,
                                                       hzar.dataGroup.add),
                                                freq$avg$analysis$oDG)
# print(summary(freq$avg$analysis$oDG$data.groups))

# Model selection
print(freq$avg$analysis$AICcTable <- hzar.AICc.hzar.obsDataGroup(freq$avg$analysis$oDG))

## Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(freq$avg$analysis$oDG$data.groups$modelI,
                         names(freq$avg$analysis$oDG$data.groups$modelI$data.param)))

## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(freq$avg$analysis$oDG$data.groups$modelI))


# (9) Import data and tidy ----------------------------------------------

# admixture_data.txt requires lat and long values; please contact author if needed
dat <- read_tsv("admixture_data.txt", col_names = TRUE)
seq <- read_tsv("../3_Admixture_index/admix_index_input_files/seq_data.txt", col_names = TRUE, col_types = cols(.default = "c"))

seq %>% 
  gather(locus, value, -population, -sample_ID, -read_data) %>% 
  filter(value != 'N') %>% 
  inner_join(pure_allele_dict, by = 'locus') %>% 
  group_by(sample_ID, population) %>% 
  summarize(frac_gent_ind = sum(value==pure_gent)/n(),
            frac_sys_ind = sum(value==pure_sys)/n()) %>% 
  filter(population != 'alterna', population !='elapsoides') %>% 
  left_join(dat %>% 
              group_by(SW_onedeglong) %>% 
              summarize(avg_long = mean(long),
                        avg_lat = mean(lat)) %>% 
              rename(population = SW_onedeglong), by = 'population') -> freq_data_inds

dat %>% 
  filter(SW_onedeglong!="elapsoides",
         SW_onedeglong!="alterna") %>% 
  dplyr::select(sample_ID, lat, long) %>% 
  left_join(freq_data_inds, by = "sample_ID") -> dist_data



# (10) Calculate distances for individual samples -------------------------

coords_inds <-
  dist_data %>% 
  dplyr::select(lat, long)

distance_inds <- geodist(coords_inds, measure="geodesic")
distance_inds <- as.data.frame(distance_inds[,1]) %>% 
  summarize(distance = `distance_inds[, 1]`/1000)

plot_data <- cbind(distance_inds, as.data.frame(freq_data_inds))

plot_data %>% 
  mutate(recip_dist = (641.9518-distance)) %>% 
  mutate(color = case_when(population == "ks_1" ~ "#fcae60",
                           population == "ks_2" ~ "#89d062",
                           population == "ks_3" ~ "#35b779",
                           population == "ks_4" ~ "#22908c",
                           population == "ks_5" ~ "#443a83",
                           population == "pure_gent" ~ "#b67431",
                           population == "pure_sys" ~ "#a8a2ca")) -> plot_data


# (11) Build and export plot (Fig. S8) ------------------------------------

# Build base plot with HZAR curve and averages per group as black points
hzar.plot.fzCline(freq$avg$analysis$oDG$data.groups$modelI, pch = 19, xlim = c(0,642),
                  xlab = "Distance (km)", ylab = "Admixture index")

# Add each individual's data point based on distance
points(plot_data$recip_dist, plot_data$frac_sys_ind, pch = 19, col=plot_data$color)

# Export the plot 11x6
ggsave("FigS8_HZAR.pdf", width=11, height=6)
