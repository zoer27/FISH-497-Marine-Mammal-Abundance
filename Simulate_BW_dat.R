#Simulation code for Blue whale abundance estimation
library(tidyverse)
library(dsims)
library(sf)
library(rnaturalearth)


#get map for plotting
US_map<-ne_countries(country = c("canada", "united states of america", "mexico"))
plot(US_map)
head(US_map)
#convert to utm zone 10
CC_coast_sf<-US_map %>% 
  st_crop(xmin = -129, xmax = -115, ymin = 30, ymax = 50) %>% 
  st_transform(crs = 32610)

ggplot(CC_coast_sf) + geom_sf()

#get survey area
eez_sf<-read_sf("eez/eez.shp")
head(eez_sf)
ggplot(eez_sf) + geom_sf() + coord_sf(ylim = c(30, 50), xlim = c(-129, -115))

CC_eez<-eez_sf %>% st_crop(xmin = -129, xmax = -115, ymin = 30, ymax = 50)
ggplot() + geom_sf(data = CC_coast_sf, fill = NA) + geom_sf(data =CC_eez)

#convert to utm zone 10
CC_sf<-st_transform(CC_eez, crs = 32610)
CC_sf

region<-make.region(region.name = "California current", 
                    shape = CC_sf)
plot(region)

#define relative density surface
density <- make.density(region = region,
                        x.space = 10000,
                        constant = 1)

plot(density, region, scale = 0.001)

density<- add.hotspot(object = density,
                      centre = c(800000, 4000000),
                      sigma = 300000, 
                      amplitude = 4)

plot(density, region, scale = 0.001)

#define the zigzag design
zigzag.design <- make.design(region = region, 
                             design = "eszigzag",
                             spacing = 100000,
                             edge.protocol = "minus",
                             design.angle = 0,
                             bounding.shape = "convex.hull",
                             truncation = 6000)
z.survey <- generate.transects(zigzag.design)
p2<-plot(region, z.survey)
p2 + geom_sf(data = CC_coast_sf)
#save image
ggsave("images/simulated_transects.png")

tracklines_sf<-z.survey@samplers
ggplot() + geom_sf()
st_write(tracklines_sf,"data/tracklines.shp", append = FALSE)
st_write(region@region, "data/region.shp", append = FALSE)
#define the popualtion
#covariates <- list(size = list(distribution = "ztruncpois", mean = 0.5))
pop.desc <- make.population.description(region = region,
                                        density = density,
                                        #covariates = covariates,
                                        N = 800,
                                        fixed.N = TRUE)
# Create the detectability description
# Define the covariate parameters on the log scale
#cov.param <- list(size = log(1.01))

detect <- make.detectability(key.function = "hn",
                             scale.param = 1500,
                             #cov.param = cov.param,
                             truncation = 6000)
plot(detect, pop.desc)

#set up some analyses (just to set up simualation but not important for this use)
analyses <- make.ds.analysis(dfmodel = list(~1, ~1),
                             key = c("hn", "hr"),
                             truncation = 6000,
                             er.var = "R2",
                             criteria = "AIC")

#set up simulation
sim.zigzag <- make.simulation(reps = 999,
                              design = zigzag.design,
                              population.description = pop.desc,
                              detectability = detect,
                              ds.analysis = analyses)
# Generate a single instance of a survey: a population, set of transects 
# and the resulting distance data
eg.zigzag.survey <- run.survey(sim.zigzag)

# Plot it to view a summary
plot(eg.zigzag.survey, region)

#save simulated data
sim_dat<-eg.zigzag.survey@dist.data

head(sim_dat)



final_sim_dat<-sim_dat %>% select(individual, distance, x, y, Effort, Region.Label, Sample.Label, Area) %>% 
  filter(!is.na(individual))

write_csv(final_sim_dat, "data/simulated_bw.csv")
