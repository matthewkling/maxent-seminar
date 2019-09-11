
library(tidyverse)
library(raster)
library(data.table)
library(dismo)

select <- dplyr::select


# load full consortium of cali herbaria dataset and climate data
cch <- read_csv("data/cch.csv")
clim <- raster("data/maxtemp.tif")

# some spatial wrangling
coordinates(cch) <- c("longitude", "latitude")
projection(cch) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
cch <- spTransform(cch, crs(clim))
cch$clim <- extract(clim, cch)

# focal species
spp <- "Sequoia sempervirens"

# prep three datasets each representing a climate distribution: 
# region-wide, cch-wide, focal species presences
region <- clim %>% 
   rasterToPoints() %>% 
   as.data.frame() %>% 
   rename(clim = maxtemp) %>% 
   mutate(dataset = "region") %>%
   sample_n(10000)
herbaria <-  as.data.frame(cch) %>%
   sample_n(10000) %>%
   rename(x=longitude, y=latitude) %>%
   select(x, y, clim) %>%
   mutate(dataset = "herbaria")
pres <- as.data.frame(cch) %>%
   rename(x=longitude, y=latitude) %>%
   filter(species == spp) %>%
   select(x, y, clim) %>%
   mutate(dataset = "presence")
d <- bind_rows(region, herbaria, pres)

# probability densities for each dataset
mn <- min(d$clim, na.rm=T)
mx <- max(d$clim, na.rm=T)
dregion <- density(d$clim[d$dataset=="region"], from=mn, to=mx)
dherb <- density(d$clim[d$dataset=="herbaria" & !is.na(d$clim)], from=mn, to=mx)
dpres <- density(d$clim[d$dataset=="presence"], from=mn, to=mx)
densities <- data.frame(clim=dregion$x, 
                        region=dregion$y, 
                        herbaria=dherb$y, 
                        presence=dpres$y) %>%
   mutate(pres_region = presence / region,
          pres_herb = presence / herbaria,
          herb_region = herbaria / region)

dd <- densities %>%
   gather(stat, value, -clim) %>%
   mutate(type = ifelse(grepl("_", stat), "ratio", "pdf"))


# map of presences on climate background
mapdata <- clim %>% rasterToPoints() %>% as.data.frame() %>% rename(clim=maxtemp)
map <- ggplot() +
   geom_raster(data=mapdata, aes(x, y, fill=clim)) +
   geom_point(data=pres, aes(x, y), color="black") +
   scale_fill_gradientn(colours=c("cyan", "dodgerblue", "darkorchid1", "red", "yellow")) +
   theme_void() +
   theme(legend.position=c(.7, .75),
         legend.title=element_text(size=20)) +
   labs(fill="summer max\ntemp, °C")
ggsave("figures/map.png", map, width=8, height=10, units="in")


# fit three maxent models
md <- filter(d, dataset %in% c("presence", "herbaria")) %>% 
   mutate(pres = as.integer(dataset=="presence")) %>%
   filter(is.finite(clim))
fit <- maxent(select(md, clim), md$pres)
md$pred <- predict(fit, md)
md$pred <- md$pred

md2 <- filter(d, dataset %in% c("presence", "region")) %>% 
   mutate(pres = as.integer(dataset=="presence")) %>%
   filter(is.finite(clim))
fit2 <- maxent(select(md2, clim), md2$pres)
md2$pred <- predict(fit2, md2)
md2$pred <- md2$pred

md3 <- filter(d, dataset %in% c("herbaria", "region")) %>% 
   mutate(pres = as.integer(dataset=="herbaria")) %>%
   filter(is.finite(clim))
fit3 <- maxent(select(md3, clim), md3$pres)
md3$pred <- predict(fit3, md3)
md3$pred <- md3$pred

# combined dataset of model predictions
pmd <- bind_rows(md %>% filter(dataset=="herbaria"),
                 md2 %>% filter(dataset=="region"),
                 md3 %>% filter(dataset=="region") %>% 
                    mutate(dataset="botanists"))

# merge with densities dataset
dd <- pmd %>%
   select(clim, dataset, pred) %>%
   rename(stat = dataset, value = pred) %>%
   mutate(type = "model",
          stat = case_when(stat=="botanists" ~ "herb_region",
                           stat=="herbaria" ~ "pres_herb",
                           stat=="region" ~ "pres_region")) %>%
   bind_rows(dd) %>%
   mutate(type = factor(type, levels=c("pdf", "ratio", "model")),
          stat = factor(stat, levels=c("presence", "region", "herbaria",
                                       "pres_region", "pres_herb", "herb_region"),
                        labels=c("presence", "land", "museums",
                                 "pres_land", "pres_museum", "museum_land")))

# plot density and maxent curves
p <- ggplot(dd %>% filter(clim > 10), 
            aes(clim, value, color=stat)) +
   facet_grid(type ~ ., scales="free") +
   geom_line(size=1) +
   theme_minimal() +
   theme(text=element_text(size=20)) +
   labs(x="summer max temp, °C",
        color=NULL, y=NULL)
ggsave("figures/curves.png", p, width=10, height=7, units="in")


# plot maxent range predictions
mapdata$pred_museum <- predict(fit, mapdata)
mapdata$pred_land <- predict(fit2, mapdata)
map <- mapdata %>% gather(stat, value, pred_museum, pred_land) %>%
   ggplot(aes(x, y, fill=value)) +
   facet_grid(.~stat) +
   geom_raster() +
   scale_fill_viridis_c() +
   theme_void() +
   theme(legend.position=c(.9, .75),
         text=element_text(size=20)) +
   labs(fill="maxent\nprediction")
ggsave("figures/maxent_maps.png", map, width=10, height=7, units="in")
