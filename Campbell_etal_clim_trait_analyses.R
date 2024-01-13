#############################################
# R-code for trait analyses in:             #
#                                           #
# Campbell C, Granath G, Rydin H.           #
# 'Trait variation across species           # 
#  distribution boundaries in Sphagnum'.    #
#                                           #
# In review American Journal of Botany.     #
#                                           #
#############################################

# Tested on R version 4.3.2, macOS 

#Libraries####
library(sf)
library(tidyverse)
library(lme4)
library(lmerTest)
library(vegan)
library(terra)
library(gstat)
library(broom)
library(cowplot)
library(egg)
library(openxlsx)
library(fractaldim)

theme_set(theme_cowplot())

# Data ####
#__Environmental ####
stk <- list.files( path="./data/env/", #
                   pattern=".grd",recursive = TRUE,
                   full.names = TRUE)
stk <- rast( stk)
# Correlation coefs between climate variables
cors <- layerCor(stk,  "pearson" ,"pearson", na.rm=TRUE)

names( stk) <- c( "ADD","AWB","SDD","SWB")
env <- read.csv( 
  list.files(  "./data/sample_data/",
               pattern ="clim_traits_patch_site_data.csv",
               full.names = TRUE,
               recursive = TRUE),stringsAsFactors = FALSE) %>%
  as_tibble() %>%
  
  st_as_sf(coords = c( "lon", "lat"),
           crs = st_crs(4326)) %>%
  st_transform(st_crs(stk)) 
env <- bind_cols( env, extract(stk, vect(env)) %>%
                    tibble() %>%
                    dplyr::select(-ID)) %>%
  st_transform(3006)


#__maxent model ####
# average output from maxent model

sp_pts <- read_delim(list.files("./data/maxent_model",
                                pattern = ".csv",
                                full.names = TRUE)) %>%
  st_as_sf(crs = st_crs(stk),
           coords = c("x", "y")) %>%
  st_transform(3006)
lfs <- list.files("./data/maxent_model/",
                  pattern = "asc",
                  full.names = TRUE)[c(1,3)]

max_sp <- rast( lfs)
names(max_sp) <- paste0("p_",
                        substr(gsub("_avg","",
                                    gsub("Sphagnum_","",names(max_sp))),1,1))
# calculate bray curtis dissimilarity
# verified by vegan::vegdist()
disim <-  (2 * min( max_sp$p_c,
                    max_sp$p_l,
                    na.rm = TRUE))/sum(max_sp,na.rm = TRUE)

names(disim) <- "disim" 
max_sp <- c(max_sp, disim)
crs(max_sp) <- crs(stk)
max_sp <- project(max_sp,
                  crs("epsg:3006"))

env <- bind_cols( env, extract(max_sp, vect(env)) %>%
                    tibble() %>%
                    dplyr::select(-ID)) %>%
  st_transform(3006)

#___response curves####
# read files and change the response names
maxent.marg <- list.files(path="./data/maxent_model",
                pattern = "\\.dat$", full.names = TRUE) 
read_file <- function(x) {
  out <- read_csv(x)
  species = if(grepl("cusp",x)) {"S.cuspidatum"} else {"S.lindbergii"}
  out$variable <- if(out$variable[1]=="ann.dd 10") {"ADD"} else 
    if(out$variable[1]=="ann.wb 10") {"AWB"} else 
      if(out$variable[1]=="seas.dd 10") {"SDD"} else {
        "SWB"
      }
  
  cbind(species=species, out)
}
maxent.marg <- lapply(maxent.marg, read_file) 

clim.sites = unique(st_drop_geometry(env[,c("species","ADD", "AWB", "SDD", "SWB")]))

# merge site data with the maxent output
maxent.marg <- lapply(maxent.marg, function (x) {
  x$samp <- NA
  dat <- clim.sites[clim.sites$species == x$species[1],
                    which(colnames(clim.sites) == x$variable[1])]
for (i in 1:NROW(dat)) {
    closestto <- which.min(abs(x$x - unlist(dat[i,])))
    x$samp[closestto] <- x$species[1]
    #print(paste(x$samp[closestto], x$variable[1], x$species[1], sep=" "))
  }
  return(x)
}
)
maxent.marg <- maxent.marg %>%
  bind_rows()


#__Canopy variables ####
ras <- list.files( path="./data/rasters/",
                   pattern=".asc",recursive = TRUE,full.names=TRUE)#raster files

ras.xyz <- grep( "li|cu",
                 ras,
                 value = TRUE)#lindbergii cuspidatum only
ras.xyz <- grep( "xyz_08",ras.xyz, value = TRUE)
can_data <- lapply(ras.xyz, function(xyz){
  
  xyz <- rast(xyz)
  nam <- gsub( "_xyz_08_","",names( xyz)) #get name of canopy
  names( xyz) <- "z"
  
  pt.xyz <- terra::as.points( xyz) #convert to points
  
  v_xyz <- gstat::variogram( z~1,data = st_as_sf( pt.xyz)) #get variogram
  
  #calculate spherical model
  vg_xyz <- fit.variogram( v_xyz,
                           vgm( psill = max( v_xyz$gamma),
                                model = "Sph",
                                range = v_xyz$dist[which( v_xyz$gamma==max( v_xyz$gamma))],
                                nugget = min( v_xyz$gamma[1:3])))  
  
  #from spherical model calculate surface roughness
  #fractal dimension
  zz <- tibble ( sample_id = nam,
                 # Calculate Ra 
                 Ra = 2*sqrt( vg_xyz$psill[2]),
                 # distance between 
                 Sr = ( 2*vg_xyz$range[2])/( 2*sqrt( vg_xyz$psill[2])), #distance between features,
                 #fractal dimension
                 Fd = fractaldim::fd.estim.variogram( bind_cols( st_coordinates( st_as_sf( pt.xyz)),
                                                                 z=pt.xyz$z))$fd)
  
  return(zz)}) %>%
  bind_rows()
can_data <- can_data %>%
  filter(Sr<3000)

#__Shoot data #### 
df.cl_stem <- read.csv( list.files(  "./data/sample_data/",
                                     pattern = "clim_traits_shoot_data.csv",
                                     full.names = TRUE,
                                     recursive = TRUE),stringsAsFactors = FALSE) %>%
  tibble() %>%
  left_join( env %>%
               st_drop_geometry() %>%
               dplyr::select(sample_id, n)) %>%
  dplyr::select( #-cap.tube.number,
    #-tube.wet.weight,
    -water.mas) %>%
  tibble( ) %>%
  
  #add env data
  left_join( env) %>%
  dplyr::select( -grep("Spreading|Pendent|^P.S",colnames(.)))  %>%#remove branch measurements
  mutate(CS = ( CPg)/( STm/SL),
         nF=NF/SL,
         nS=NSlf/SL,
         STg = STm/SL) %>%
  dplyr::select( sample_id,
                 stem_id,
                 p_c,
                 p_l,
                 AWB,
                 ADD,
                 SWB,
                 SDD,
                 species,
                 HWT,
                 disim,
                 site,
                 wc,
                 n,
                 CPg,
                 STg,
                 nF,
                 nS,
                 Pbr,
                 Sbr,
                 CS) 

#___add two canopy traits####
# STBD and CPBD
df.cl <- left_join(df.cl_stem,
                   df.cl_stem %>%
                     group_by(sample_id,species,n) %>%
                     summarise(STg = mean(STg, na.rm = TRUE),
                               CPg = mean(CPg, na.rm = TRUE))%>%
                     ungroup() %>%
                     mutate(STBD = STg* n,
                            CPBD = CPg * n) %>%
                     left_join(can_data) %>%
                     dplyr::select(-STg, -CPg))%>%
  mutate(species = gsub("[.]",". ",species))

# Water Content  ####

mod <- lmer( data=df.cl,( 100*wc)~species-1+(species)/HWT+  #model slope comparison
               ( 1|site)+
               ( 1|sample_id))
mod.test <- lmer( data=df.cl,( 100*wc)~species*scale(HWT)+  #model slope comparison
               ( 1|site)+
               ( 1|sample_id))

# F test
anova(mod.test,type=2)
summary(mod.test)

# CIs
wc.ci <- confint( mod,level=0.95,method="boot")

model <- bind_cols(`Fixed effects` = rownames(summary(mod)$coefficients),
                   coeffs = as_tibble(summary(mod)$coefficients),
                  `2.5%` = unname(wc.ci[4:7,1]),
                  `97.5%` = unname(wc.ci[4:7,2]))

# Save data frame for plotting
df.cl.plot <- df.cl

# Height above the Water level####
# HWT and range differences
hwt.means <- df.cl %>% 
  group_by(sample_id, species, site) %>%
  summarise(meanHWT = mean(HWT, na.rm = TRUE))

mod.test <- lmer( data=hwt.means, meanHWT ~ species +
                    ( 1|site) )
anova(mod.test,type=2)

# range
hwt.means.sites <- df.cl %>% 
  group_by(site, species) %>%
  summarise(rangeHWT = max(HWT, na.rm = TRUE) - min(HWT, na.rm = TRUE)) 
t.test(rangeHWT ~ species, data=hwt.means.sites,
       alternative = "two.sided", var.equal = FALSE)

# Arrange trait data set####
df.cl <- df.cl %>%
  group_by(sample_id,species,site) %>%
    summarise(across(where(is.numeric),function(q){
      mean(q,na.omit=TRUE)
    })) %>%
  ungroup()

st_mods_df <- split( df.cl %>%
                   dplyr::select(species, site,AWB, ADD, SWB, SDD, #disim,
                                 CPg,STg,CS,nF,nS,Sbr, Pbr,
                                 n,Ra,Sr,Fd,STBD,CPBD) %>%
                     distinct() %>%
                   pivot_longer(cols = c(CPg,STg,CS,nF,nS,Sbr, Pbr,
                                         n,Ra,Sr,Fd,STBD,CPBD),
                                names_to = "trait",
                                values_to = "tvalue"),df.cl %>%
                   dplyr::select(species, site,AWB, ADD, SWB, SDD, #disim,
                                 CPg,STg,CS,nF,nS,Sbr, Pbr,
                                 n,Ra,Sr,Fd,STBD,CPBD) %>%
                     distinct() %>%
                   pivot_longer(cols = c(CPg,STg,CS,nF,nS,Sbr, Pbr,
                                         n,Ra,Sr,Fd,STBD,CPBD),
                                names_to = "trait",
                                values_to = "tvalue") %>%
                   .$trait)

# Trait models####
# For correct statistical Type II test we need to use a centralised
# continuous predictor, but for the figure we need the model to
# be fitted with non-centralised predictor. Hence we fit the same model
# twice.
# First for anova testing
smodsLmer.anova <-  lapply( st_mods_df,function( q){
  return( list(ADD = lmer( data=q,tvalue ~ species * scale(ADD,scale=F)+( 1|site)),
               SDD = lmer( data=q,tvalue ~ species * scale(SDD, scale=F)+( 1|site)),
               AWB = lmer( data=q,tvalue ~ species * scale(AWB, scale=F)+( 1|site)),
               SWB = lmer( data=q,tvalue ~ species * scale(SWB, scale=F)+( 1|site))
#               disim = lmer( data=q,tvalue~species * disim + (  1|site)))
  ))
})

# Now for coefs
smodsLmer <-  lapply( st_mods_df,function( q){
  return( list(ADD = lmer( data=q,tvalue ~ species * ADD+( 1|site)),
               SDD = lmer( data=q,tvalue ~ species * SDD+( 1|site)),
               AWB = lmer( data=q,tvalue ~ species * AWB+( 1|site)),
               SWB = lmer( data=q,tvalue ~ species * SWB+( 1|site))
               #               disim = lmer( data=q,tvalue~species * disim + (  1|site)))
  ))
})
# Extract coefs for plotting
straitModCoef <- lapply(names(smodsLmer), function(q){
lapply(names(smodsLmer[[q]]), function(w){
    return(tibble(species = c("S. cuspidatum", "S. lindbergii"),
          climate = w,
          trait = q,
          intercet = c(smodsLmer[[q]][[w]]@beta[1],sum(smodsLmer[[q]][[w]]@beta[1:2])),
          slope    = c(smodsLmer[[q]][[w]]@beta[3],sum(smodsLmer[[q]][[w]]@beta[3:4]))))
  }) %>%
  bind_rows()
}) %>%
  bind_rows()

# Anovas for traits
sanvSum <- lapply(names(smodsLmer.anova), function(q){
  mods <- lapply(names(smodsLmer.anova[[q]]), function(w){
    df.adj = if(VarCorr(smodsLmer.anova[[q]][[w]])[1]$site[1] == 0) {"Kenward-Roger"
      } else {"Satterthwaite"}
  anv <- anova(smodsLmer.anova[[q]][[w]], type=2, ddf = df.adj)
  anv <- tidy(anv) %>%
        mutate(trait = q,
               climate = w) %>%
    mutate(statistic = round(statistic, 3)) %>%
    arrange(term)
  
    return(anv)
  }) %>%
    bind_rows()
  
  }) %>%
    bind_rows() %>%
  mutate(p.adj = p.adjust(p.value, "BH"))

# Change the names in the vector to remove the "scale"-stuff
sanvSum$term[grepl("ADD", sanvSum$term) & grepl(":", sanvSum$term)] <- "species:ADD"
sanvSum$term[grepl("SDD", sanvSum$term) & grepl(":", sanvSum$term)] <- "species:SDD"
sanvSum$term[grepl("AWB", sanvSum$term) & grepl(":", sanvSum$term)] <- "species:AWB"
sanvSum$term[grepl("SWB", sanvSum$term) & grepl(":", sanvSum$term)] <- "species:SWB"
#sanvSum$term[grepl("disim", sanvSum$term) & grepl(":", sanvSum$term)] <- "species:disim"
sanvSum$term[grepl("ADD", sanvSum$term) & !(grepl(":", sanvSum$term))] <- "ADD"
sanvSum$term[grepl("SDD", sanvSum$term) & !(grepl(":", sanvSum$term))] <- "SDD"
sanvSum$term[grepl("AWB", sanvSum$term) & !(grepl(":", sanvSum$term))] <- "AWB"
sanvSum$term[grepl("SWB", sanvSum$term) & !(grepl(":", sanvSum$term))] <- "SWB"
#sanvSum$term[grepl("disim", sanvSum$term) & !(grepl(":", sanvSum$term))] <- "disim"

slabel <- sanvSum  %>%
  mutate(label = gsub("species","S",gsub("ADD|SDD|AWB|SWB", #|disim",
                                         "C",gsub(":","x",term)))) %>%
  #mutate(climate = w)%>%
  dplyr::select(climate,trait,p.adj, label) %>%
  pivot_wider(names_from = label, values_from = p.adj) %>%
  
  mutate(label = paste0("C = ",ifelse(C<0.001,"***, ",
                                      ifelse(C<0.01,"**, ",
                        #                     ifelse(C<0.05,"*, ", "NS, "))),
                                              ifelse(C<0.05,"*, ", paste(round(C,2),", ",sep="")))),
                        "S = ", ifelse(S<0.001,"***, ",
                                       ifelse(S<0.01,"**, ",
                                             # ifelse(S<0.05,"*, ", "NS, "))),
                                              ifelse(S<0.05,"*, ", paste(round(S,2),", ",sep="")))),
                        
                        "S x C = ", ifelse(SxC<0.001,"*** ",
                                           ifelse(SxC<0.01,"** ",
                                                  #ifelse(SxC<0.05,"* ", "NS ")))))
                                                  ifelse(SxC<0.05,"*", round(SxC,2))))))


#___save model output for appendix####
 # write.xlsx(
 #   list(statistic = sanvSum %>%
 #          mutate(statistic = round(statistic, 1)) %>%
 #     dplyr::select(-p.value,-p.adj,
 #                   -sumsq,-meansq, -DenDF) %>%
 #     pivot_wider(names_from = trait,
 #                 values_from = statistic) %>%
 #       relocate(climate) %>%
 #       relocate(CPg,STg,CS,nF,nS,Sbr,Pbr,
 #                n,Ra,Sr,Fd,STBD,CPBD, .after = c(NumDF)),
 #     DenDF = sanvSum %>%
 #       mutate(DenDF = round(DenDF, 1)) %>%
 #       dplyr::select(-p.value,-statistic,
 #                     -sumsq,-meansq, -p.adj) %>%
 #       pivot_wider(names_from = trait,
 #                   values_from = DenDF) %>%
 #       relocate(climate) %>%
 #       relocate(CPg,STg,CS,nF,nS,Sbr,Pbr,
 #                n,Ra,Sr,Fd,STBD,CPBD, .after = NumDF),
 #     p.value = sanvSum %>%
 #       mutate_at(.vars = vars(p.adj),
 #     .funs = ~ case_when(round(.x, digits = 2) == 0 ~ round(.x, digits = 5),
 #                         TRUE ~ round(.x, digits = 2))) %>%
 #      dplyr::select(-p.value,-statistic,
 #                     -sumsq,-meansq, -DenDF) %>%
 #       pivot_wider(names_from = trait,
 #                   values_from = p.adj) %>%
 #       relocate(climate) %>%
 #       relocate(CPg,STg,CS,nF,nS,Sbr,Pbr,
 #                n,Ra,Sr,Fd,STBD,CPBD, .after = NumDF))
 #   ,
 #     "./results/Appendix_S2_v2.xlsx")

# PCA ordination ####
df.clSum <- (df.cl %>%
dplyr::select(site,species,n,CPg,STg,nF,nS,Pbr,Sbr,CS,Ra,Sr,Fd,STBD,CPBD) %>%
  pivot_longer(-c(site, species)) %>%
  pivot_wider(names_from = name, values_from = value,
              values_fn = function(q){mean(q, na.rm = TRUE)})) %>%
  na.omit(.)
df.clSum = data.frame(df.clSum)

# we use capital N for capitulum density in ms so need to change here
colnames(df.clSum)[3] <- "N"

# data for PCA
pca.dat = df.clSum %>%
  select(-c(site, species)) %>%
  na.omit(.)

# climate variables data
clim.sites.pca = data.frame(unique(st_drop_geometry(
  env[,c("site", "species","ADD", "AWB", "SDD", "SWB")])))

# fix a space in the species names
clim.sites.pca$species <- 
  ifelse(clim.sites.pca$species == "S.cuspidatum", "S. cuspidatum",
         "S. lindbergii")

# order the two data frames in the same order
clim.sites.pca$order <- paste(clim.sites.pca[,1],clim.sites.pca[,2], sep="")
df.clSum$order <- paste(df.clSum[,1], df.clSum[,2], sep="")

clim.sites.pca = clim.sites.pca[order(match(
  clim.sites.pca$order,df.clSum$order)
),]

# run pca
pca <- rda(pca.dat, scale = TRUE) 
sumpca = summary(pca)
varexp1 = round((sumpca$cont$importance[2, "PC1"])*100,digits=0)
varexp2 = round((sumpca$cont$importance[2, "PC2"])*100,digits=0)

# save scores
pca.scores <- plot(pca, scaling = 3)
pca_data <- data.frame(sample=df.clSum$site, 
                       site.type = ifelse(df.clSum$site %in% 
                                            c("S1","S2","S9", "S10"), "allo",
                                          "symp"),
                       species = df.clSum$species,
                       x=pca.scores$sites[,1],
                       y=pca.scores$sites[,2])
# save loadings
loadingscores <- data.frame(x=pca.scores$species[,1],
                            y=pca.scores$species[,2],
                            labs = colnames(pca.dat))

# fit env variables onto PCA
ef <- envfit (pca ~ ADD+AWB+SDD+SWB, 
              data = clim.sites.pca, perm = 1999)
ef
arrow_heads <- ef$vectors$arrows
env.var.scores <- data.frame(scores(ef, display="vectors", arrow.mul=2.5),
                             variables=c("ADD","AWB","DDS","WBS"))

