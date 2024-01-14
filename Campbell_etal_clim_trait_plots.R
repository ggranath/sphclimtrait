#############################################
# R-code for plotting figures in:           #
#                                           #
# Campbell C, Granath G, Rydin H.           #
# 'Trait variation across species           # 
#  distribution boundaries in Sphagnum'.    #
#                                           #
# In review.                                #
#                                           #
#############################################

# Get data and output from analyses
# source("./Campbell_etal_clim_trait_analyses.R")

# Libraries#####
library(ggplot2)
library(colorRamps)
library(tidyterra)
library(ggspatial)

# Fig 2. MaxEnt figures and sampling sites####
env$sympsites <- ifelse(env$site %in% c("S1", "S2", "S9", "S10"), "allop", "symp")

dist_plots <- lapply( names( max_sp), 
                      function(q){
                        species <- ifelse( q == "p_c", "Sphagnum cuspidatum",
                                           ifelse( q == "p_l", "Sphagnum lindbergii","Dissimilarity"))
                        
                        p1 <- ggplot( )+
                          geom_spatraster(data=max_sp[[q]],
                                          show.legend = TRUE)+
                          scale_fill_gradientn( colors=matlab.like2( n=100), 
                                                na.value = "white", limits=c(0,1))+
                          scale_x_continuous(expand=c(0,0)) + 
                          scale_y_continuous(expand=c(0,0), 
                                             limits=c(6110949, 7800555)) +
                          theme(legend.position = "bottom",axis.line=element_blank(),axis.text.x=element_blank(),
                                axis.text.y=element_blank(),axis.ticks=element_blank(),
                                axis.title.x=element_blank(),
                                axis.title.y=element_blank(),
                                panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                                panel.grid.minor=element_blank(),
                                plot.background=element_blank(),
                                plot.margin = unit(c(0,0,0,0),"cm")) +
                          guides(fill = "none")+
                          coord_sf(ylim = c(57, 71), 
                                   clip = "off")
                        #      if(species == "Sphagnum cuspidatum") {
                        #        guides(fill = guide_colorbar(barheight = 1, barwidth = 10))
                        #      } else {guides(fill = "none")
                        #          }
                        
                        
                        
                        if(species != "Dissimilarity"){
                          p1 <- p1+
                            geom_sf(data = sp_pts %>%
                                      filter(sp == species),
                                    size = 0.1)
                          
                          
                        } else {
                          p1 <- p1 +
                            geom_sf(data = env, aes(color=sympsites),
                                    size = 3, shape=1)+
                            scale_color_manual(values = c('allop' = "deeppink",  'symp'= "black"),
                                               labels=c('allopatric', 'sympatric'))+
                            guides(color = "none") +
                            annotation_scale(width_hint =0.5, 
                                             location ="tl", 
                                             pad_y=unit(0.3, "cm"),
                                             pad_x=unit(0, "cm"),
                                             text_cex=1.3 ) 
                        }
                        
                        return(p1)         })

maxent_plot <- plot_grid( dist_plots[[1]]+
                            guides(fill = guide_legend(title = "Suitability probability (a,b) \nBC similarity (c)"))+
                            theme(legend.title = element_text(size = 10),
                                  legend.justification='left',
                                  legend.box.margin = margin(0,0,0,-1.5,"cm")),
                          dist_plots[[2]] +
                            guides(fill = "none"), 
                          dist_plots[[3]]+
                            guides(color = guide_legend(title = "Sampled \nsites",
                                                        override.aes = list(fill = NA)),
                                   linetype = guide_legend(override.aes = list(fill = NA))) +
                            theme(legend.title = element_text(size = 10), 
                                  legend.justification='left',
                                  legend.box.margin = margin(0,0,0,-1.55,"cm"),
                                  legend.key = element_rect(fill = "white"),
                                  legend.background=element_blank()
                            ),
                          ncol = 3,
                          labels = c( "( a)","( b)","( c)"),align = "h", axis = 'b') 
maxent_plot
ggsave2("./plots/fig2.png",
        plot = maxent_plot,
        height= 70*1.4,
        width = 160*1.4,
        dpi = 400, bg="white",
        units = "mm")

# Fig 3. MaxEnt response curves####

mylimits <- function(x) range(scales::breaks_extended()(x))

mybreaks <- function(n = 3) {
  function(x) {
    breaks <- mylimits(x)
    seq(breaks[1], breaks[2], length.out = n)  
  }
}

fig3 <- ggplot(data= maxent.marg, aes(x=x, y=y, colour = species)) +
  geom_line() +
  scale_colour_manual( name="species",
                       breaks=c( "S.cuspidatum",
                                 "S.lindbergii"),
                       values = c( "green",
                                   "orange"),
                       labels = c("S. cuspidatum", "S. lindbergii")) +
  #scale_x_continuous(breaks = mybreaks(n = 3), limits = mylimits) +
  ylab("Probability of suitability") + 
  xlab(NULL)+
  facet_wrap(~ variable, scale="free_x",
             labeller = as_labeller(c(
               `ADD` = 'Annual~Degree~Days~"("*degree*"C"~y^{-1}*")"', 
               `SDD` = 'Degree~Day~Seasonality~"("*degree*"C"~y^{-1}*")"',
               `AWB` = 'Annual~Water~Balance~"(mm "*y^{-1}*")"',
               `SWB` = 'Water~Balance~Seasonality~"(mm "*y^{-1}*")"'),
               label_parsed),
             strip.position = "bottom") +
  theme(
    strip.placement = "outside",   
    strip.background = element_blank(),
    legend.position = "bottom",
    legend.justification = "center",
    legend.text = element_text( face="italic"),
  ) +
  geom_point(data= maxent.marg[!(is.na(maxent.marg$samp)),], 
             aes(x=x, y=y), shape=19, color="black")
fig3
ggsave2("./plots/fig3.png",
        plot = fig5,
        width = 170,
        dpi = 140, bg="white",
        units = "mm")

# Fig. 4 water content ####

# End points in WC graph to draw lines confined by data range
Xs = do.call(data.frame, aggregate(HWT~species, df.cl.plot, FUN = function (x) c(min(x), max(x))))
betas = c(model$Estimate[grep( "tum$",model$`Fixed effects`)],model$Estimate[grep( "gii$",model$`Fixed effects`)],
          model$Estimate[grep( "tum:HWT",model$`Fixed effects`)], model$Estimate[grep( "gii:HWT",model$`Fixed effects`)])
ys = c(cbind(1, t(Xs[1,2:3])) %*% betas[c(1,3)],
       cbind(1, t(Xs[2,2:3])) %*% betas[c(2,4)])

# extract mean and standard error inside the ggplot function
calculate_mean_and_se <- function(data) {
  n <- length(data)
  y <- mean(data)
  ymin <- y - sd(data) / sqrt(n)
  ymax <- y + sd(data) / sqrt(n)
  data.frame(y=y, ymin=ymin, ymax=ymax)
}

p.wc <- ggplot( data=df.cl.plot, aes( x=( HWT),y=wc*100, colour = species))+
  stat_summary(geom = "linerange", fun.data = calculate_mean_and_se)+
  stat_summary(geom = "point", fun.data = calculate_mean_and_se, shape=1, size=2)+
  #stat_summary( fun.data = "mean_cl_boot", shape = 1)+
  scale_colour_manual( name="species",
                       breaks=c( "S. cuspidatum",
                                 "S. lindbergii"),
                       values = c( "green",
                                   "orange"))+
  geom_segment(aes(x = Xs[1,2], xend = Xs[1,3], y = ys[1], yend = ys[2]),
               colour = "green") +
  geom_segment(aes(x = Xs[2,2], xend = Xs[2,3], y = ys[3], yend = ys[4]),
               colour = "orange") +
#  geom_abline( slope=model$Estimate[grep( "tum:HWT",model$`Fixed effects`)],
#               intercept=model$Estimate[grep( "tum$",model$`Fixed effects`)],
#               colour="green")+
#  
#  geom_abline( slope=model$Estimate[grep( "gii:HWT",model$`Fixed effects`)],
#               intercept=model$Estimate[grep( "gii$",model$`Fixed effects`)],               
#               colour="orange")+
  xlab( "Height above Water Table (cm)")+
  ylab( "Water Content (%)")+
  scale_x_continuous(n.breaks=10, limits= c(-2, 16)) +
  scale_y_continuous(n.breaks=10, limits= c(0, 4000)) +
  guides(colour = guide_legend(title  = "Species"))+
  theme( legend.position = "bottom",
         legend.text = element_text( face="italic", size = 10),
         axis.title = element_text(size=12),
         axis.text = element_text(size=10),
         legend.title = element_text(size = 12))+
  NULL

p.wc
ggsave2("./plots/height_above_water_table_ver2.png",
        plot = p.wc,
        height = 110,
        width = 140,
        dpi = 300,bg="white",
        units = "mm")

# Fig. 5 shoot traits plot ####
# this is for placing the stats text in the graph
spos <- tibble(expand.grid(climate = unique(straitModCoef$climate),
                          trait = unique(straitModCoef$trait),
                          stringsAsFactors = FALSE)) %>%
  left_join(tibble(climate = c("ADD","SDD","AWB","SWB", "disim"),
                   x = c(1750,260,125,42.5, .65))) %>%
  left_join(tibble(trait =c("CPg","STg","CS","nF","nS","Sbr", "Pbr"),
                   y = c(0.04,0.018,6,21,65,1.4,1.3)))

# Filter data for trait plot
trait.plot.dat <- df.cl.plot %>%
  dplyr::select(species, AWB, ADD, SWB, SDD,# disim,
                CPg,STg,CS,nF,nS,Sbr, Pbr) %>%
  #filter(climate != "disim") %>%
  distinct() %>%
  pivot_longer(cols = c(ADD, AWB, SWB, SDD),
               names_to  = "climate",
               values_to = "cvalue") %>%
  pivot_longer(cols = c(CPg,STg,CS,nF,nS,Sbr, Pbr),
               names_to = "trait",
               values_to = "tvalue")

# End points in trait plot to draw regression lines confined by data range
rr = df.cl.plot %>% 
  pivot_longer(cols = c(CPg,STg,CS,nF,nS,Sbr, Pbr),
               names_to = "trait",
               values_to = "trait.dat") %>%
  pivot_longer(cols = c(AWB,ADD,SWB,SDD),
               names_to = "climate",
               values_to = "climate.dat") %>%
  group_by(species, climate, trait)  %>%
  summarise(xmin = min(climate.dat, na.rm = TRUE), xmax= max(climate.dat, na.rm = TRUE),
            .groups = "keep") %>%
  arrange(factor(trait, levels = c("CPg","CS","nF","nS","Pbr","Sbr","STg")), 
          factor(climate, levels = c("ADD", "SDD", "AWB", "SWB")),
          factor(species, levels = c("S. cuspidatum", "S. lindbergii")))
  # get the coefficients
betas = straitModCoef %>%
  filter(trait %in% c("CPg","STg","CS","nF","nS","Sbr", "Pbr"))
  # take end points x coefs to get outet values of the line
ranges = data.frame(ymin = NA, ymax = NA)
for ( i in 1:NROW(rr)) {
  ranges[i,] =t( 
    cbind(1, t(rr[i,4:5])) %*% unlist(c(betas[i,4:5])))
}
line.lims = cbind(betas[,1:3], rr[,4:5], ranges)
  
  # plot the graph
p <- ggplot(data = trait.plot.dat, aes(x=cvalue, y = tvalue
                         )) +
  stat_summary(geom = "linerange", fun.data = calculate_mean_and_se,
               aes(colour = species))+
  stat_summary(geom = "point", fun.data = calculate_mean_and_se, 
               aes(colour = species), shape=1, size=2)+
  #stat_summary(fun.data = "mean_cl_boot",aes(colour = species))+
  #geom_smooth( method="lm",aes( linetype=species),colour="black",se=FALSE)+
  #geom_abline(data = straitModCoef %>%
  #              filter(trait %in% c("CPg","STg","CS","nF","nS","Sbr", "Pbr")),
  #            aes( lty = species,
  #            slope = slope,
  #            intercept = intercet))+
  geom_segment(data=line.lims, aes(x = xmin, xend = xmax, y = ymin, 
                                   yend = ymax, lty = species)) +
  scale_colour_manual( name="species",
                       breaks=c( "S. cuspidatum",
                                 "S. lindbergii",
                                 "label"),
                       values = c( "green",
                                   "orange"))+
  facet_grid( trait~climate,scales="free",
              labeller = labeller( climate=c( "ADD"="Annual Degree Days",
                                              "SDD"="Degree Day Seasonality",
                                              "AWB"="Annual Water Balance",
                                              "SWB"="Water Balance Seasonality"),
                                            #  "disim" = "Disimilarity"),
                                   trait=c( "CS"="Capitulum:Stem mass",
                                            "Pbr"="Pendent branch",
                                            "Sbr"="Spreading branch",
                                            "nF"="Fascicle density",
                                            "nS"="Stem leaf density",
                                            "CPg"="Capitulum mass",
                                            "STg"="Stem mass")))+
  scale_y_continuous(n.breaks = 6) +
  theme_bw( )+
  theme( legend.position = "bottom",
         strip.text = element_text( face="bold"),
         strip.background = element_blank( ),
         #  legend.title = element_text( 12),
         legend.text = element_text( face = "italic"),
         axis.text = element_text( ),
         axis.title = element_blank( ),
         panel.grid = element_blank( ))
p
p.trait_stm <-  p+
  geom_text( data = spos %>%
               filter(trait %in% c("CPg","STg","CS","nF","nS","Sbr", "Pbr")) %>%
               left_join(slabel %>%
                           dplyr::select(-c(S,C,SxC))) %>%
               mutate(x=x*c(0.77,0.98,0.68,0.96)),
             aes(x = x, y= y, label = label), size = 2)

ggsave2("./plots/stem_traits vs climate_ver_jan11_v2.png",
        plot = p.trait_stm,
        width = 180,
        height = 270,
        dpi = 300,
        units = "mm")

# Fig. 6 canopy traits plot####
# Filter data for canopy plot
spos <- tibble(expand.grid(climate = unique(straitModCoef$climate),
                           trait = unique(straitModCoef$trait),
                           stringsAsFactors = FALSE)) %>%
  left_join(tibble(climate = c("ADD","SDD","AWB","SWB", "disim"),
                   x = c(1750,260,125,42.5, .65))) %>%
  left_join(tibble(trait =c("CPBD", "Fd","n","Ra","Sr","STBD"),
                   y = c(0.04, 2.3, 1.7, 0.6, 1750, 0.024)))


canopy.plot.dat <- df.cl.plot %>%
  dplyr::select(species, AWB, ADD, SWB,  SDD, #disim,
                n,STBD,CPBD,Ra,Sr,Fd)  %>%
  distinct(.)%>%
  pivot_longer(cols = c(ADD, AWB, SWB, SDD), #, disim
               names_to  = "climate",
               values_to = "cvalue") %>%
  pivot_longer(cols = c( n,STBD,CPBD,Ra,Sr,Fd),
               names_to = "trait",
               values_to = "tvalue")

# End points in trait plot to draw regression lines confined by data range
rr.can = df.cl.plot %>% 
  pivot_longer(cols = c(n,STBD,CPBD,Ra,Sr,Fd),
               names_to = "trait",
               values_to = "trait.dat") %>%
  pivot_longer(cols = c(AWB,ADD,SWB,SDD), #, disim
               names_to = "climate",
               values_to = "climate.dat") %>%
  group_by(species, climate, trait)  %>%
  summarise(xmin = min(climate.dat, na.rm = TRUE), xmax= max(climate.dat, na.rm = TRUE),
            .groups = "keep") %>%
  arrange(factor(trait, levels = c("CPBD", "Fd","n","Ra","Sr","STBD")), 
          factor(climate, levels = c("ADD", "SDD", "AWB", "SWB")), # , "disim"
          factor(species, levels = c("S. cuspidatum", "S. lindbergii")))
  # get coefs
betas.can = straitModCoef %>%
  filter(trait %in% c("CPBD", "Fd","n","Ra","Sr","STBD"))
  
# take end points x coefs to get line
ranges.can = data.frame(ymin = NA, ymax = NA)
for ( i in 1:NROW(rr.can)) {
  ranges.can[i,] =t( 
    cbind(1, t(rr.can[i,4:5])) %*% unlist(c(betas.can[i,4:5])))
}
line.lims.can = cbind(betas.can[,1:3], rr.can[,4:5], ranges.can)

# plot canopy trait graph
p.can <- ggplot(data = canopy.plot.dat , aes(x=cvalue, y = tvalue)) +
  stat_summary(geom = "linerange", fun.data = calculate_mean_and_se,
               aes(colour = species))+
  stat_summary(geom = "point", fun.data = calculate_mean_and_se, 
               aes(colour = species), shape=1, size=2) +
  
  #stat_summary(fun.data = "mean_cl_boot",aes(colour = species))+
      #geom_smooth( method="lm",aes( linetype=species),colour="black",se=FALSE)+
  #geom_abline(data = ctraitModCoef,
  #            aes( lty = species,
  #                 slope = slope,
  #                 intercept = intercet))+
  geom_segment(data=line.lims.can, aes(x = xmin, xend = xmax, y = ymin, 
                                   yend = ymax, lty = species)) +
  scale_colour_manual( name="species",
                       breaks=c( "S. cuspidatum",
                                 "S. lindbergii",
                                 "label"),
                       values = c( "green",
                                   "orange"))+
  facet_grid( trait~climate,scales="free",
              labeller = labeller( climate=c( "ADD"="Ann Degree Days",
                                              "SDD"="Degree Day Season",
                                              "AWB"="Ann Water Balance",
                                              "SWB"="Water Balance Season"),
                                              #"disim" = "Dissimilarity"),
                                   trait=c( "n"="Capitulum density",
                                            "CPBD"="Capitulum bulk density",
                                            "STBD"="Stem bulk density",
                                            "Ra"="Surface roughness",
                                            "Sr"="Roughness scale",
                                            "Fd"="Fractal dimension")))+
  scale_y_continuous(expand = expansion(mult = c(0,0.05)), n.breaks=6) +
  theme_bw( )+
  theme( legend.position = "bottom",
         strip.text = element_text( face="bold"),
         strip.background = element_blank( ),
         #  legend.title = element_text( 12),
         legend.text = element_text( face = "italic"),
         axis.text = element_text( ),
         axis.title = element_blank( ),
         panel.grid = element_blank( ))
p.can

# this is to ad coordinates for added text
# this can be removed now
cpos <- tibble(expand.grid(climate = unique(straitModCoef$climate),
                           trait = unique(straitModCoef$trait),
                           stringsAsFactors = FALSE)) %>%
  left_join(tibble(climate = c("ADD","SDD","AWB","SWB"), #, "disim"),
                   cvalue = c(1750,260,125,42.5)))  %>% #, .7))) %>%
  left_join(tibble(trait =c("n","CPBD","STBD","Ra","Sr","Fd"),
                   tvalue = c( 1.8,
                               0.045,
                               0.03,
                               0.65,
                               2000,
                               2.32)))
# # add stats text to plot
# p.trait_can <- p.can +
#   #geom_text( data = cpos %>%
#    #            left_join(clabel %>%
#     #                       dplyr::select(-c(S,C,SxC))),
#      #        aes(x = cvalue, y= tvalue, 
#       #           label = label), 
#        #      size = 3)
# geom_text( data = spos %>%
#              left_join(slabel %>%
#                          dplyr::select(-c(S,C,SxC))),
#            aes(x = -Inf, y = Inf, label = label), 
#            hjust = -0.025, vjust = 1.05, size=2.5)

p.trait_can <- p.can +
  geom_text( data = spos %>%
               filter(trait %in% c("CPBD", "Fd","n","Ra","Sr","STBD")) %>%
               left_join(slabel %>%
                           dplyr::select(-c(S,C,SxC))) %>%
               mutate(x=x*c(0.77,0.98,0.68,0.96)),
             aes(x = x, y= y, label = label), size = 2)

ggsave2("./plots/canopy_traits vs climate_t2_jan13.png",
        plot = p.trait_can,
        width = 190,
        height = 280,
        dpi = 300,
        units = "mm")


#Fig 7. PCA####
# to adjust variable labels
adj.var.place = data.frame(y = c(-0.1,     0,  0.05, -0.1, -0.1,  0.06,
                                 0.1,  0.1,     0, -0.03, -0.03,
                                 -0.08, 0), 
                           x =c(0.05,    0.15,  0.15,    0,    0,  0.08,
                                0,-0.05,-0.15,  0.1,  -0.15,
                                0.2, 0.2))
# make plot fig7
fig7 = ggplot(pca_data, aes(x=x, y=y, color=species, shape=site.type)) +
  geom_text(aes(label=sample), size=3) +
  geom_point(pca_data[pca_data$site.type == "allo",], 
             mapping=aes(x=x, y=y, color=species, 
                         shape=site.type), 
             size=6,  inherit.aes = F) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  xlim(c(-1.5, 2)) +
  ylim(c(-1.5, 1.5)) +
  xlab(paste("PC1 (",varexp1,"%)", sep="")) + 
  ylab(paste("PC2 (",varexp2,"%)", sep="")) + 
  geom_text(data=loadingscores, aes(y=y, x=x, label=labs), 
            inherit.aes = F, parse=FALSE, 
            nudge_y=adj.var.place$y, nudge_x=adj.var.place$x) + 
  geom_text(data=env.var.scores, aes(y=PC2, x=PC1, label=variables), 
            inherit.aes = F, parse=FALSE, 
            nudge_y=c(0.07,-0.05,0,0), nudge_x=c(0,0,0.2,-0.2)) + 
  scale_shape_manual( name=NULL, breaks= c("allo"),
                      labels = c("allopatric sites"),
                      values = c(1,2),
                      guide = guide_legend(override.aes = 
                                             list(size = 5))) +
  scale_colour_manual( name="species",
                       breaks=c( "S. cuspidatum",
                                 "S. lindbergii"),
                       values = c( "green",
                                   "orange"),
                       guide = guide_legend(
                         label.theme = element_text(angle = 0, face = "italic"),
                         override.aes = list(
                           label = c("S", "S"),size=4, shape=""))) +
  theme(legend.position = "bottom",
        legend.box="vertical", 
        legend.justification = "center",
        legend.margin=margin()) +
  geom_segment(data = loadingscores, inherit.aes = FALSE,
               aes(x = 0, xend = x, 
                   y = 0, yend = y),
               arrow = arrow(length = unit(0.025, "npc"), type = "open"), 
               lwd = 0.5, color="red") +
  geom_segment(data = env.var.scores, inherit.aes = FALSE,
               aes(x = 0, xend = PC1, 
                   y = 0, yend = PC2),
               arrow = arrow(length = unit(0.025, "npc"), type = "open"), 
               lwd = 0.5, color= "black") 
fig7
ggsave2("./plots/fig7_pca.png",
        plot = fig7,
        height = 150,
        width = 130,
        dpi = 300,bg="white",
        units = "mm")


