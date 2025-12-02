sapply(c('dplyr', 'ggplot2', 'lubridate', 'forcats', 
         'magrittr', 'cmdstanr', 'rethinking', 'cowplot', 
         'bipartite', 'igraph', 'iNEXT', 'hypervolume', 'readxl'), 
       library, character.only = T)

extrafont::loadfonts(device = 'win')

source('functions_mod_diagnostics.r')


path_alien_native <- '/Users/andres/Documents/github_repos/hawaii_projects/fecal_songbirds/native_plants_alien.xlsx'

alien_native1 <- read_xlsx(path_alien_native, sheet = 1, col_names = T)
alien_native2 <- read_xlsx(path_alien_native, sheet = 2, col_names = T)
alien_native2$native <- ifelse(alien_native2$native == 'A', F, T)

alien_native <- 
  full_join(alien_native1, 
            alien_native2, 
            by = c('sp_plant', 'native')) |> unique()

df <- readRDS('FecalSongbirds_CORRECTED.rds')

codes_birds <- readRDS('birds_code.rds')

codes_plants <- readRDS('plants_code.rds')
codes_plants2 <- read_xlsx('plants_codes2.xlsx', sheet = 1, col_names = T)

plant_traits <- readRDS('plant_traits.rds')

bird_traits <- readRDS('bird_traits.rds')

df <- df[, c("Site", "Sample_Date", "Month", "Bird_sp", 
             "PlantSpecies", "SeedCount")]

colnames(df) <- c('site', 'date', 'month', 'bird', 'plant', 'seeds')

df <- df[df$plant != 'none', ]

df$code <- df %$% paste0(bird, '_', plant)

df <- df[!grepl('[0-9]', df$code), ]

df <- df[!is.na(df$seeds), ]

df$plant <- toupper(df$plant)

codes_plants$code <- toupper(codes_plants$code)

missing_plants <- unique(df$plant)[!(unique(df$plant) %in% unique(codes_plants2$Acronym))]


plants_codes_network <- codes_plants2[codes_plants2$Acronym %in% unique(df$plant), ]

codes_plants$species <- 
  gsub('^(.*)(\\s)([a-z]*)(\\s)(.*)$', '\\1\\2\\3', codes_plants$species)

colnames(plants_codes_network) <- c('species', 'code', 'origin')

plants_codes_full <- full_join(codes_plants, plants_codes_network, by = c('code', 'species'))

plants_codes_full  <- plants_codes_full[, 1:4]

plants_codes_full$genus <- gsub('^(.*)(\\s)(.*)$', '\\1', plants_codes_full$species)

codes_temp <- 
  cbind(
  plants_codes_full[which(is.na(plants_codes_full$family)), c(1:2, 5)], 
  tribble(~family,       ~order, 
          'Urticaceae',  'Rosales', 
          'Araliaceae',  'Apiales', 
          'Moraceae',     'Rosales', 
          'Moraceae',     'Rosales', 
          'Moraceae',     'Rosales',
          'Verbenaceae',  'Lamiales', 
          'Euphorbiaceae', 'Malpighiales', 
          'Passifloaraceae', 'Malpighiales', 
          'Dipentodontaceae', 'Huerteales', 
          'Araliaceae', 'Apiales', 
          'Boraginaceae', 'Boraginales', 
          'Euphorbiaceae', 'Malpighiales', 
          'Cucurbitaceae', 'Cucurbitales')
) |> as_tibble()

plants_codes_full <- 
  full_join(plants_codes_full[-which(is.na(plants_codes_full$family)), ], 
            codes_temp, by = c('code', 'species', 'genus', 'family', 'order')) 
plants_codes_full <- 
  rbind(plants_codes_full, 
      tibble(code = missing_plants, 
             species = NA, 
             family = 'Arecaceae', 
             order = 'Arecales', 
             genus = NA))

names(plant_traits)[!(names(plant_traits) %in% plants_codes_full$species)] <- 'Sambucus nigra'

plant_traits$`Sambucus nigra`$sp <- 'Sambucus nigra'

codes_birds$Acronym <- toupper(codes_birds$Acronym)


interactions <- unique(df$code)

rarefaction_dat <- split(df, df$site)

rarefaction_dat <- 
  lapply(rarefaction_dat, FUN = 
           function(x) {
             
             sapply(interactions, FUN = 
                      function(i) {
                        length(which(x$code == i))
                      })
             
           })

rarefaction_plot <- iNEXT(rarefaction_dat, q = 0, datatype = 'abundance')


ggiNEXT(rarefaction_plot) +
  theme_minimal() +
  labs(x = 'Interaction events', y = 'Interaction richness')

df <- split(df, list(df$site, df$month))

df <- df[unlist(lapply(df, function(x) nrow(x) > 0), use.names = F)]

par(mfrow = c(1, 2))
plot(density(unlist(lapply(df, function(x) length(unique(x$bird))), use.names = T)), 
     main = '', xlab = 'Bird species per network')
plot(density(unlist(lapply(df, function(x) length(unique(x$plant))), use.names = T)), 
     main = '', xlab = 'plant species per network')
par(mfrow = c(1, 1))

N_networks <- 
  do.call('rbind', 
          lapply(df, FUN = 
                   function(x) {
                     tibble(N_birds = length(unique(x$bird)), 
                            N_plants = length(unique(x$plant)), 
                            site = x$site[1], 
                            month = x$month[1], 
                            N = length(unique(x$date))) # N sampling
                   }))

final_nets_indx <- 
  lapply(df, FUN = 
         function(x) {
           i <- length(unique(x$bird))
           j <- length(unique(x$plant))
           
           i >= 2 & j >=2
         }) |> unlist()

df <- df[final_nets_indx]

matrix_networks <- 
  lapply(df, FUN = 
           function(x) {
             x$code <- toupper(x$code)
             d_temp <- unique(x[, c('bird', 'plant')])
             d_temp <- as.matrix(d_temp)
             d_temp <- expand.grid(d_temp[,1], d_temp[, 2])
             d_temp <- as_tibble(unique(d_temp))
             colnames(d_temp) <- c('bird', 'plant')
             d_temp$code <- d_temp %$% paste(bird, plant, sep = '_')
             d_temp <- full_join(x, d_temp, by = c('bird', 'plant', 'code'))
             indx <- which(is.na(d_temp$seeds))
             d_temp$seeds[indx] <- 0
             d_temp$site[indx] <- d_temp$site[1]
             d_temp$month[indx] <- d_temp$month[1]
             d_temp$year <- year(d_temp$date)
             
             d_temp <-
               d_temp |>
               unique() |>
               group_by(year, code) |>
               mutate(seeds = sum(seeds)) |>
               unique() |>
               ungroup() |>
               group_by(code) |>
               mutate(seeds = mean(seeds)) |>
               ungroup() |>
               dplyr::select(site, month, code, bird,
                             plant, seeds) |>
               unique() |>
               mutate(freq_rel = seeds/sum(seeds))
             d_temp <- d_temp[order(d_temp$seeds, decreasing = T), ]

             birds <- unique(d_temp$bird)
             plants <- unique(d_temp$plant)

             sapply(plants, FUN =
                      function(i) {
                        sapply(birds, FUN =
                                 function(j) {
                                   d_temp[d_temp$bird == j &
                                            d_temp$plant == i,]$freq_rel
                                 })
                      })
             
           })


matrix_networks_plots <- 
  lapply(matrix_networks, FUN = 
           function(x) {
             
             w <- t(x)
             w <- w[which(w > 0)]
             #net <- graph_from_biadjacency_matrix(x)
             net <- graph_from_incidence_matrix(x) # to be run in mac
             V(net)$size <- igraph::degree(net) * 5
             V(net)$label <- ''
             V(net)$color <- c(rep('seagreen', nrow(x)), rep('tan1', ncol(x)))
             E(net)$weight <- w
             net
           })


par(mfrow = c(8, 4), mar = c(1, 1, 1, 1))
for (i in seq_along(matrix_networks)) {
  plot(matrix_networks_plots[[i]], 
       layout = layout_in_circle, 
       edge.width = E(matrix_networks_plots[[i]])$weight * 10)
  text(x = 0, y = 1.1, lab = names(matrix_networks_plots)[i])
}
par(mfrow = c(1, 1))

metrics_networks <- 
  do.call('rbind', 
          lapply(seq_along(matrix_networks), FUN = 
                   function(x) {
                     t <- matrix_networks[[x]]
                     
                     modu <- computeModules(t)
                     
                     tibble(size = nrow(t) + ncol(t), 
                            N_birds = nrow(t), 
                            N_plants = ncol(t), 
                            asimetry = (N_birds - N_plants) / (N_plants + N_birds),
                            nestedness = networklevel(t, index = 'weighted NODF'),
                            modularity = modu@likelihood,
                            H2 = networklevel(t, index = 'H2'), 
                            site = names(matrix_networks)[x])
                   }))

par(mfrow = c(2, 2), mar = c(4, 4, 1, 1))
for (i in 1:4) {
  plot(metrics_networks[[i]], metrics_networks$nestedness, 
       xlab = colnames(metrics_networks)[i], ylab = 'nestedness')
}

for (i in 1:4) {
  plot(metrics_networks[[i]], metrics_networks$H2, 
       xlab = colnames(metrics_networks)[i], ylab = 'H2')
}

for (i in 1:4) {
  plot(metrics_networks[[i]], metrics_networks$modularity, 
       xlab = colnames(metrics_networks)[i], ylab = 'Modularity')
}
par(mfrow = c(1, 1))

metrics_networks$month <- gsub('(.*)([0-9][0-9])$', '\\2', 
                               metrics_networks$site)

metrics_networks$site2 <- gsub('([A-Z]*)(.*)$', '\\1', 
                               metrics_networks$site)

metrics_networks$month <- as.numeric(metrics_networks$month)


names(bird_traits) <- toupper(names(bird_traits))

bird_traits <- 
  lapply(bird_traits, FUN = 
           function(x) {
             x$sp <- toupper(x$sp)
             x
           })

plant_traits <- 
  lapply(plant_traits, FUN = 
         function(x) {
           sp <- unique(x$sp)
           indx <- which(plants_codes_full$species == sp)
           cod <- unique(plants_codes_full$code[indx])
           
           if (length(cod) == 1) {
             x$code <- cod
           } else {
             x$code <- cod[2]
           }
           x
         })

names(plant_traits) <- 
  unlist(lapply(plant_traits, FUN = function(x) unique(x$code)))

(indx <- do.call('rbind', df)$plant %in% sort(names(plant_traits)))

plants_imputation <- 
  do.call('rbind', df)$plant[!indx] |> unique() # those are the plants to
                                              # conduct imputation
df_relative_abu <- 
  lapply(df, FUN = 
         function(x) {
           
           # Here we calculated the relative abundance of plant and bird
           # species in the network. Then, we use it to calculate the 
           # number of individuals per species in a simulated community 
           # of 2e3 individuals
           
           x$year <- year(x$date)
           
           bird <- 
             x |> 
             dplyr::select(site, year, month, bird) |> 
             group_by(bird, year) |> 
             mutate(bird_abu = length(month)) |> 
             ungroup() |> unique() |> 
             mutate(bird_abu = (bird_abu/sum(bird_abu))*100, 
                    ind_perc = round((2e3 * bird_abu)/100), 
                    t = sum(ind_perc))
           
           plant <- 
             x |> 
             dplyr::select(site, year, month, plant) |> 
             group_by(plant, year) |> 
             mutate(plant_abu = length(month)) |> 
             ungroup() |> unique() |> 
             mutate(plant_abu = (plant_abu/sum(plant_abu))*100, 
                    ind_perc = round((2e3 * plant_abu)/100), 
                    t = sum(ind_perc)) 
           
           list(bird = bird, 
                plant = plant)
         })

df_relative_abu <- 
  lapply(df_relative_abu, FUN = 
         function(x) {
           j <- lapply(x, FUN = 
                         function(i) {
                           
                           # this is necessary to guarantee that the summation 
                           # of individuals in plants and birds is equal to 2e3.
                           # We need this to bind the two matrices of traits
                           
                           tot <- unique(i$t)
                           
                           if (tot < 2e3) {
                             i <- i[order(i$ind_perc, decreasing = T), ]
                             diff <- 2e3 - i$t[1]
                             i$ind_perc[1] <- i$ind_perc[1] + diff
                           }
                           
                           if (tot > 2e3) {
                             i <- i[order(i$ind_perc, decreasing = T), ]
                             diff <- i$t[1] - 2e3
                             i$ind_perc[1] <- i$ind_perc[1] - diff
                           }
                           
                           i$t2 <- sum(i$ind_perc)
                           i 
                         })
           j
         })


fun_imputation <- 
  function(sp_code = plants_imputation[1], sim_ind = 1e3) {
    d <- plants_codes_full
    d2 <- d[d$code != sp_code, ]
    tax <- d[which(d$code == sp_code), ]
    
    # Here we detect the species matching genus, family or order 
    # of the species to conduct imputation, and extract indexes to 
    # index the list of simulated traits 
    df_gen <- na.omit(d2[d2$genus == tax$genus, ])
    indx_gen_trait <- names(plant_traits) %in% df_gen$code
    df_fam <- na.omit(d2[d2$family == tax$family, ])
    indx_fam_trait <- names(plant_traits) %in% df_fam$code
    df_ord <- na.omit(d2[d2$order == tax$order, ])
    indx_ord_trait <- names(plant_traits) %in% df_ord$code
    
    traits_imputation <- # function to conduct trait imputation 
      function(df) { # list of species and posterior predictions of the traits
        imputed_traits <- 
          lapply(seq_along(df[[1]]), FUN = 
                 function(j) {
                   trait <- do.call('cbind', # binding columns of same traits
                                    lapply(df, function(i) i[, j]))
                   
                   trait <- apply(trait, 1, mean) # imputation trait i
                   matrix(trait[1:sim_ind], # simulated individuals
                          ncol = 1)
                 })
        
        # binding imputed traits on a single matrix
        imputed_traits <- do.call('cbind', imputed_traits) 
        colnames(imputed_traits) <- colnames(df[[1]])
        imputed_traits
      } # end of function `traits_imputation`
    
    # The following statements control the imputation depending on the 
    # available species at genus, family or oderd level.
    if (nrow(df_gen) > 0 & sum(indx_gen_trait) > 0) {
      message(paste0('Imputing across genus: ', unique(df_gen$genus)))
      message('...')
      
      d_trait <- lapply(plant_traits[indx_gen_trait], # selecting species of same genus
                        function(x) x[, 2:7]) 
      
      imputation_out <- traits_imputation(d_trait)
      imputation_out
      
    } else if (nrow(df_fam) > 0 & sum(indx_fam_trait) > 0) {
      message(paste0('Imputing across family: ', unique(df_fam$family)))
      message('...')
      
      d_trait <- lapply(plant_traits[indx_fam_trait], # selecting species of same family
                        function(x) x[, 2:7]) 
      
      imputation_out <- traits_imputation(d_trait)
      imputation_out
      
    } else if (nrow(df_ord) > 0 & sum(indx_ord_trait) > 0) {
      message(paste0('Imputing across order: ', unique(df_ord$order)))
      
      d_trait <- lapply(plant_traits[indx_ord_trait], # selecting species of same order
                        function(x) x[, 2:7]) 
      
      imputation_out <- traits_imputation(d_trait)
      imputation_out
      
    } else {
      message(paste0('There is not genus, family or orders to impute the species ', 
                    sp_code, '. Imputing across all species.'))
      message('...')
      
      d_trait <- lapply(plant_traits, function(x) x[, 2:7]) # using all species
      
      imputation_out <- traits_imputation(d_trait)
      imputation_out
      
    }
  }

extract_functional_space <- # function to estimate the hypervolume
  function(month_site = 'EKA.01') {
    
    # first, we separete data of birds and plants
    plants <- df_relative_abu[[month_site]]$plant
    birds <- df_relative_abu[[month_site]]$bird
    
    # This is to detect if we have to impute one or more plant 
    no_imputation <- sum(plants$plant %in% plants_imputation) == 0
    
    # the estimations will consider the structure of the community
    HV_birds <- 
      lapply(seq_along(birds$bird), FUN = 
               function(i) {
                 # extracting the traits of N individuals of
                 # the species i depending on its relative abundance 
                 
                 perc <- birds$ind_perc[i]
                 
                 sp <- birds$bird[i]
                 
                 indx <- which(names(bird_traits) == sp)
                 
                 traits <- bird_traits[[indx]]
                 
                 as.matrix(traits[1:perc, 2:4])
                 
               })
    
    HV_birds <- do.call('rbind', HV_birds)
    
    indx_imputation <- plants$plant %in% plants_imputation
    
    message(paste0('We have ', sum(indx_imputation), ' sp/spp to impute'))
    
    if (no_imputation) {
      HV_plants <-
        lapply(seq_along(plants$plant), FUN =
                 function(i) {
                   
                   # number of individuals to simulate
                   perc <- plants$ind_perc[i]
                   
                   # code of the plant
                   sp <- plants$plant[i]
                   
                   # select the spp in the posterior predictive data 
                   indx <- which(names(plant_traits) == sp)
                   
                   traits <- plant_traits[[indx]]
                   
                   # bind the simulated individuals in a 
                   as.matrix(traits[1:perc, 2:7])
                   
                 })
      
      HV_plants <- do.call('rbind', HV_plants)
    } else {
      plants_imp <- plants[indx_imputation, ] # select the species to impute
      
      plants_imp <- lapply(1:nrow(plants_imp), FUN = 
                             function(i) {
                               
                               # number of individuals to simulate 
                               N_inds <- plants_imp$ind_perc[i] 
                               # code of the plant
                               code_plant <- plants_imp$plant[i]
                               
                               # imputation
                               fun_imputation(code_plant, sim_ind = N_inds)
                             })
      
      plants_imp <- do.call('rbind', plants_imp)
      
      plants_no_imp <- plants[!indx_imputation, ] # select spp that do not need imputation
      
      plants_no_imp <-
        lapply(seq_along(plants_no_imp$plant), FUN =
                 function(i) {
                   
                   # number of individuals to simulate
                   perc <- plants_no_imp$ind_perc[i]
                   
                   # code of the plant
                   sp <- plants_no_imp$plant[i]
                   
                   # select the spp in the posterior predictive data 
                   indx <- which(names(plant_traits) == sp)
                   
                   traits <- plant_traits[[indx]]
                   
                   # bind the simulated individuals in a 
                   as.matrix(traits[1:perc, 2:7])
                   
                 })
      
      plants_no_imp <- do.call('rbind', plants_no_imp)
      
      # binding plants with and without imputation
      HV_plants <- rbind(plants_no_imp, plants_imp) 
      
    }
    
    # We bind the birds and plants matrices to generate the trait matrix
    # for the network.
    HV_network <- cbind(HV_birds, HV_plants)

    # We reduce the dimensionality of the functional space of the network to
    # improve hypervolume calculation (only gape, fruit diameter and nutritional
    # traits are considered)
    HV_network <- HV_network[, -c(1:2, 5:6)]

    # Here we index plants matrix to estimate hypervolumes of the
    # morphological and nutritional spaces
    HV_plants_morfo <- HV_plants[, 1:3]
    HV_plants_nut <- HV_plants[, 4:6]

    # Hypervolume estimations
    message('Calculating network hypervolume')
    set.seed(5)
    HV_networks <- hypervolume(HV_network, method = 'svm', svm.gamma = 0.5)

    message('Calculating plants hypervolume (total)')
    set.seed(5)
    HV_plants <- hypervolume(HV_plants, method = 'svm', svm.gamma = 0.5)

    message('Calculating plants hypervolume (morphological)')
    set.seed(5)
    HV_plants_morfo <- hypervolume(HV_plants_morfo, method = 'svm', svm.gamma = 0.5)

    message('Calculating plants hypervolume (nutritional)')
    set.seed(5)
    HV_plants_nut <- hypervolume(HV_plants_nut, method = 'svm', svm.gamma = 0.5)

    message('Calculating birds hypervolume')
    set.seed(5)
    HV_birds <- hypervolume(HV_birds, method = 'svm', svm.gamma = 0.5)

    message('Hypervolumes calculation done')
    # extract volumes
    tibble(site = month_site,
           HV_network = HV_networks@Volume,
           HV_plant = HV_plants@Volume,
           HV_plants_morfo = HV_plants_morfo@Volume,
           HV_plants_nut = HV_plants_nut@Volume,
           HV_bird = HV_birds@Volume)
  }


t1 <- Sys.time()
hypervolumes <-
  lapply(seq_along(df_relative_abu), FUN =
             function(x) {
               i <- names(df_relative_abu)[x]
               message(paste('Running site', x,
                             'from', length(names(df_relative_abu)),
                             'sites'))
               extract_functional_space(i)
             })
Sys.time() - t1

hypervolumes <- do.call('rbind', hypervolumes)

hypervolumes$site <- names(df_relative_abu)

rainfall <- readRDS('AVG_climate_data.rds')[[1]]
temperature <- readRDS('AVG_climate_data.rds')[[2]]

sites <- 
  paste(gsub('^([A-Z]*)(\\.)([0-9]*)$', '\\1', hypervolumes$site), 
      as.numeric(gsub('^([A-Z]*)(\\.)([0-9]*)$', '\\3', hypervolumes$site)), 
      sep = '_')

hypervolumes$site <- sites

metrics_networks$site <- metrics_networks %$% paste(site2, month, sep = '_')

colnames(rainfall)[1] <- 'z_rainfall'

rainfall$z_rainfall <- as.vector(scale(rainfall$z_rainfall))

rainfall$site <- 
  rainfall %$% paste(sites, month, sep = '_')

colnames(temperature)[1] <- 'z_temperature'
temperature$z_temperature <- as.vector(scale(temperature$z_temperature))

temperature$site <- 
  temperature %$% paste(sites, month, sep = '_')

dat <- full_join(metrics_networks, hypervolumes, by = 'site')

dat <- left_join(dat, 
                 temperature[, c("site", "z_temperature")], 
                 by = 'site')

left_join(dat, 
          rainfall[, c("site", "z_rainfall")], 
          by = 'site') |> 
  apply(2, function(x) sum(is.na(x)))


dat <- 
  left_join(dat, 
            rainfall[, c("site", "z_rainfall")], 
            by = 'site')

#saveRDS(dat, 'data_for_models.rds')

sessionInfo()
