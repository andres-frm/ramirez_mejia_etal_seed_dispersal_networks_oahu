pks <- 
  c('dplyr', 'magrittr', 'lubridate', 'lubridate',
    'readxl', 'tibble', 'tidyr', 'ggplot2', 
    'ggfortify', 'terra', 'raster', 'sf', 'patchwork', 
    'cmdstanr', 'rethinking', 'animation')

sapply(pks, FUN = function(x) library(x, character.only = T))

extrafont::loadfonts(device = 'win')


path_bird_fecal <- 
  '/Users/andres/Documents/github_repos/hawaii_projects/data/FecalSongbirds_v1.xls'

path_bird_morpho <- 
  '/Users/andres/Documents/github_repos/hawaii_projects/bird_morphology/BirdMorphol_JMG-Master_2019May17.xlsx'

bird_morpho <- read_xlsx(path_bird_morpho, 
                         sheet = 1, col_names = T)

date_birds <- bird_morpho$date

bird_morpho <- 
  lapply(seq_along(bird_morpho), FUN = 
           function(i) {
             
             x <- bird_morpho[[i]]
             
             n <- sum(grep('^([0-9])', x))
             
             if (n > 0) {
               df <- tibble(v = as.numeric(x))
               colnames(df) <- colnames(bird_morpho)[i]
               df
             } else {
               bird_morpho[, i]
             }
             
           })

bird_morpho <- as_tibble(do.call('cbind', bird_morpho))

bird_morpho$date <- date_birds

bird_morpho <- 
  bird_morpho[, c("site", "date", 'net', 'species', 
                  "band_number", 'age', 'sex', 'tail', 
                  'tarsus', 'wing', "total_culmen", 
                  "nares_to_tip", "width", "depth", "gape", 
                  "mass")]

bird_morpho <- split(bird_morpho, bird_morpho$species)

bird_morpho <- bird_morpho[unlist(lapply(bird_morpho, FUN = function(x) nrow(x) > 0))]

bird_morpho <- do.call('rbind', bird_morpho)

bird_seed <- read_xls(path_bird_fecal, sheet = 1)

bird_seed <- bird_seed[!is.na(bird_seed$Month), ]

bird_seed$Month <- 
  sapply(bird_seed$Month, 
         function(x) {
           if (x < 10) paste('0', x, sep = '')
           else as.character(x)
         })


bird_seed$Day <- 
  sapply(bird_seed$Day, 
         function(x) {
           if (nchar(x) == 1) paste0('0', x)
           else x
         })

bird_seed$Sample_Date <- 
  bird_seed %$% paste(Year, Month, Day, sep = '-')


bird_seed$Sample_Date <- as.Date(bird_seed$Sample_Date)

codes_birds <- read_xls(path_bird_fecal, sheet = 2)

codes_plants <- read_xls(path_bird_fecal, sheet = 3)

coordinates_sites <- read_xls(path_bird_fecal, sheet = 4)

bird_morpho$date <- as.Date(bird_morpho$date)

dates_morpho <- unique(bird_morpho[, c('site', 'date')])

date_seeds <- unique(bird_seed[, c("Site", "Sample_Date")])

colnames(date_seeds) <- c('site', 'date')

date_seeds$dat_seed <- 1

dates_1 <- unique(bird_morpho[, c('site', 'date')])

dates_1$dat_trait <- 1

dates_1$site[grep('^PAL$', dates_1$site)] <- 'PAH'

dates_1 <- dates_1[dates_1$site != 'HOUSE', ]

date_join <- left_join(date_seeds, dates_1, by = c('site', 'date'))

# using band number to asses compatibility of dates

bird_morpho$band_number <- 
  sapply(bird_morpho$band_number, FUN = 
           function(x) {
             if (is.na(x)) {
               NA
             } else {
               x <- as.character(x)
               gsub('^(....)(.*)', '\\1-\\2', x)
             }
           })

code_traits <- unique(bird_morpho[, c('site', "date", "band_number")])
code_seeds <- unique(bird_seed[, c("Site", "Sample_Date", "Band_Num")])
colnames(code_seeds) <- colnames(code_traits)

compare_bands <- 
  function(date1, site1, date2, site2) {
    seed <- code_seeds[code_seeds$date == date1 , ]
    trait <- code_traits[code_traits$date == date2, ]
    
    list(bands =
           list(trait = sort(trait$band_number), site_trait = unique(trait$site), 
                seeds = sort(seed$band_number)), site_seed = unique(seed$site), 
         prop_same = mean(seed$band_number %in% trait$band_number),
         prop_total = mean(seed$band_number %in% code_traits$band_number))
    
  }


date_join$correction <- date_join$date

t <- date_join
t$day <- day(t$date)
t$filter <- 
  sapply(seq_along(t$day), FUN = 
           function(x) {
             i <- t$day[x] <= 12
             j <- is.na(t$dat_trait[x])
             (i+j) == 2
           })

t$year <- year(t$date)

t1 <- t[(t$filter), ] %$% aggregate(date ~ site + year, FUN = length)

t1 <- t1[order(t1$site), ]

t1$perc <- (t1$date/sum(t1$date))*100

t <- t[!(t$filter), ]

t$cod <- t %$% paste(site, date, sep = '_')

bird_seed$cod <- bird_seed %$% paste(Site, Sample_Date, sep = '_')

bird_seed1 <- split(bird_seed, bird_seed$cod)

indx <- names(bird_seed1) %in% t$cod

bird_seed1 <- do.call('rbind', bird_seed1[indx])

# data for constructing interaction networks
#saveRDS(bird_seed1, 'FecalSongbirds_CORRECTED.rds')



path_environmental <- 
  '/Users/andres/Documents/github_repos/hawaii_projects'

coords <- 
  tribble(~lon,            ~lat, 
          -158.08451944444, 21.443305555556,
          -158.19316666667, 21.536819444444, 
          -157.87335833333, 21.376041666667,
          -158.1447694444, 21.506827777778,
          -158.1799,       21.536472222222, 
          -157.81091388889, 21.338386111111,
          -158.04074722222, 21.632116666667)

coords$sites <- coordinates_sites$Acronym[-8]

folderRain <- 
  paste0(path_environmental, 
         '/rainfall_temperature/HCDP_data/rainfall/data_map/', 
         '2015')

mapsRain <- dir(folderRain)

path_rain <- paste0(folderRain, '/', mapsRain[1])

rainfall <- rast(path_rain)

folderTEMP <- 
  paste0(path_environmental, 
         '/rainfall_temperature/HCDP_data/temperature/data_map/', 
         '2015')

mapsTEMP <- dir(folderTEMP)

path_TEMP <- paste0(folderTEMP, '/', mapsTEMP[1])

temperature <- rast(path_TEMP)


par(mfrow = c(1, 2), mar = c(1, 1, 4, 1))
plot(rainfall, xlim = c(-158.3, -157.6), ylim = c(21.2, 21.75), 
     main = 'Rainfall (January 2014)')
coords %$% 
  points(lon, lat, col = 'red', pch = 16, cex = 0.5)
plot(temperature, xlim = c(-158.3, -157.6), ylim = c(21.2, 21.75), 
     main = 'Temperature (January 2014)')
coords %$% 
  points(lon, lat, col = 'red', pch = 16, cex = 0.5)
par(mfrow = c(1, 1))

rainfall_data <- 
  lapply(2014:2017, FUN = 
           function(x) {
             
             folder <- 
               paste0(path_environmental, 
                      '/rainfall_temperature/HCDP_data/rainfall/data_map/', 
                      x)
             
             maps <- dir(folder)
             
             path <- paste0(folder, '/', maps)
             
             dat <- lapply(path, FUN = 
                             function(i) {
                               rainfall <- rast(i)
                               d <- extract(rainfall, coords[, -3])
                               month <- gsub('^(.*)(..)(.tif)$', '\\2', i)
                               date <- as.Date(paste0(x, '-', month, '-', '01'))
                               d$date <- date
                               d$month <- as.numeric(month)
                               d$year <- x
                               colnames(d)[2] <- 'rainfall_mm'
                               d$sites <- coordinates_sites$Acronym[-8]
                               d[, -1]
                             })
             
             dat <- do.call('rbind', dat)
             dat
           })

rainfall_data <- as_tibble(do.call('rbind', rainfall_data))

rainfall_data$z_rainfall <- as.vector(scale(rainfall_data$rainfall_mm))

rainfall_data <- split(rainfall_data, rainfall_data$sites)

rainfall_data <- 
  do.call('rbind', 
          lapply(rainfall_data, FUN = 
                   function(i) {
                     i <- i[i$year != 2014, ]
                     i$month_id <- 1:nrow(i)
                     i
                   }))

temp_data <- 
  lapply(2014:2017, FUN = 
           function(x) {
             
             folder <- 
               paste0(path_environmental, 
                      '/rainfall_temperature/HCDP_data/temperature/data_map/', 
                      x)
             
             maps <- dir(folder)
             
             path <- paste0(folder, '/', maps)
             
             dat <- lapply(path, FUN = 
                             function(i) {
                               temp <- rast(i)
                               d <- extract(temp, coords[, -3])
                               month <- gsub('^(.*)(..)(.tif)$', '\\2', i)
                               date <- as.Date(paste0(x, '-', month, '-', '01'))
                               d$date <- date
                               d$month <- as.numeric(month)
                               d$year <- x
                               colnames(d)[2] <- 'temperature'
                               d$sites <- coordinates_sites$Acronym[-8]
                               d[, -1]
                             })
             
             dat <- do.call('rbind', dat)
             dat
           })

temp_data <- as_tibble(do.call('rbind', temp_data))

temp_data$z_temperature <- as.vector(scale(temp_data$temperature))

temp_data <- split(temp_data, temp_data$sites)

temp_data <- 
  do.call('rbind', 
          lapply(temp_data, FUN = 
                   function(i) {
                     i <- i[i$year != 2014, ]
                     i$month_id <- 1:nrow(i)
                     i
                   }))

# Data for climatic models

# saveRDS(list(temperature = temp_data, 
#              rainfall = rainfall_data), 'rainfall_temperature.rds')

temp_data <- readRDS('rainfall_temperature.rds')[[1]]

rainfall_data <- readRDS('rainfall_temperature.rds')[[2]]

dat_rainfall <- 
  list(
    N = nrow(rainfall_data),
    M = 12,
    N_year = 3,
    L = max(as.numeric(as.factor(rainfall_data$sites))), 
    t = rainfall_data$month_id,
    year_id = as.numeric(as.factor(rainfall_data$year)),
    month = rainfall_data$month,
    level = as.numeric(as.factor(rainfall_data$sites)),
    rainfall = rainfall_data$rainfall_mm,
    period = 12
  )


cat(file = 'climate_model.stan', 
    '
    data {
      int N;                // Total observations (252)
      int M;                // N months (12)
      int L;                // levels of factor (sites) (7)
      int N_year;           // N year
      array[N] int year_id;  // years
      array[N] int t;       // number of time points (1, ..., 36)
      array[N] int month;   // month indices (1, ..., 12)
      array[N] int level;   // factor (grouping variable) (1, ..., 7)
      vector[N] rainfall;   // rainfall data
      int period;           // period of circular time (12)
    }
    
    parameters {
      real alpha;
      real<lower = 0> sigma;                // noise for the likelihood
      vector<lower = 0>[L] length_scale;    // length scale (smooth term)
      vector<lower = 0>[L] sigma_f;         // noise for GP
      array[L] vector[M] eta;    // latent variables for each month
      vector[N_year] z_year;
      real mu_year;
      real<lower = 0> sigma_year;
    }
    
    transformed parameters {
      array[L] vector[M] f;
      for (l in 1:L) {
        matrix[M, M] K;
        matrix[M, M] L_K;
        
        for (i in 1:(M-1)) {
          for (j in (i+1):M) {
            real distance = abs(i - j);
            real periodic_distance = fmin(distance, period - distance);
            // periodic kernel
            K[i, j] = sigma_f[l]^2 * exp(-2 * square(sin(pi()*periodic_distance / period))/square(length_scale[l]));
            K[j, i] = K[i, j];
          }
          K[i, i] = sigma_f[l]^2 + 1e-9;
        }
        K[M, M] = sigma_f[l]^2 + 1e-9;
        
        // cholesky decomposition
        L_K = cholesky_decompose(K);
        
        // transforme the latent variable eta to the GP
        f[l] = L_K * eta[l];
      }
    
      vector[N_year] year;
      year = mu_year + z_year * sigma_year;
    
    }
    
    model {
      sigma ~ cauchy(0, 1);
      length_scale ~ inv_gamma(5, 5);
      sigma_f ~ cauchy(0, 1);
      for (l in 1:L) {
        eta[l] ~ normal(0, 1);
      }
    
      alpha ~ normal(5, 1);
    
      z_year ~ normal(0, 1);
      mu_year ~ normal(0, 0.25);
      sigma_year ~ exponential(1);
    
      for (n in 1:N) {
        rainfall[n] ~ lognormal(alpha + f[level[n]][month[n]] + 
                                year[year_id[n]], sigma);
      }
    }
    
    generated quantities {
      vector[N] y_pred;
      
      for (i in 1:N) {
        y_pred[i] = lognormal_rng(alpha + f[level[i]][month[i]] + 
                                  year[year_id[i]], sigma);
      }
    }
    ')


file <- paste0(getwd(), '/climate_model.stan')

fit_climate <- cmdstan_model(file, compile = T)

mod_temperature <-
  fit_climate$sample(
    data = dat_temperature,
    chains = 3,
    parallel_chains = 3,
    iter_sampling = 1e4,
    iter_warmup = 1e3,
    thin = 10,
    seed = 123
  )

summary_rain <- mod_rainfall$summary()

summary_rain %$% plot(rhat ~ ess_bulk)
summary_rain %$% points(rhat ~ ess_tail, col = 'red')




ppcheck_rainfall <- mod_rainfall$draws('y_pred', format = 'matrix')

plot(density(dat_rainfall$rainfall), main = '', 
     ylim = c(0, 0.008), xlim = c(-40, 900), xlab = 'Cumulative rainfall')
for (i in 1:100) lines(density(ppcheck_rainfall[i, ]), lwd = 0.1)
lines(density(dat_rainfall$rainfall), col = 'red', lwd = 2)



post_rainfall <- mod_rainfall$draws(c('alpha', 'f', 'year', 'sigma'), format = 'df')

post_rainfall <- 
  list(year = post_rainfall[, grep('^year', colnames(post_rainfall))],
       alpha = post_rainfall$alpha, 
       f = post_rainfall[, grep('^f', colnames(post_rainfall))], 
       sigma = post_rainfall$sigma)

dat_rainfall <- as_tibble(do.call('cbind', dat_rainfall[6:9]))

site_month <- unique(dat_rainfall[, c('level', "month")])

est_rainfall2 <- 
  lapply(1:12, FUN = 
           function(x) {
             
             cc <- paste0('^(.*)(,', x, '\\])$')
             
             mu <- mean(dat_rainfall$rainfall)
             sigma <- sd(dat_rainfall$rainfall)
             
             indx <- grep(cc, colnames(post_rainfall$f))
             
             post <- apply(post_rainfall$f[, indx, drop = T], 1, mean)
             a <- post_rainfall$alpha
             year <- apply(post_rainfall$year, 1, mean)
             sigma <- post_rainfall$sigma
             
             mu_est <- a + post + year
             
             mu_est <- exp(mu_est)
             
             tibble(mu = mean(mu_est) , 
                    li_mu = quantile(mu_est, 0.025), 
                    ls_mu = quantile(mu_est, 0.975),
                    month = x)
             
           })

est_rainfall2 <- do.call('rbind', est_rainfall2)

for (i in 2:3) {
  est_rainfall2[[i]] <- 
    (est_rainfall2[[i]] - mean(est_rainfall2$mu)) / 
    sd(est_rainfall2$mu)
}

est_rainfall2$mu <- as.vector(scale(est_rainfall2$mu))

est_rainfall <- 
  lapply(1:nrow(site_month), FUN = 
           function(x) {
             site <- site_month$level[x]
             month <- site_month$month[x]
             
             f_GP <- 
               paste0('f[', site, ',', month, ']')
             
             post <- post_rainfall$f[, f_GP, drop = T]
             a <- post_rainfall$alpha
             year <- apply(post_rainfall$year, 1, mean)
             sigma <- post_rainfall$sigma
             
             mu_est <- a + post + year
             mu_est <- mu_est
             
             set.seed(555)
             post_pred <- rlnorm(length(post), mu_est, sigma)
             
             tibble(mu = mean(exp(mu_est)), 
                    li_mu = quantile(exp(mu_est), 0.025), 
                    ls_mu = quantile(exp(mu_est), 0.975), 
                    li_pred = quantile(post_pred, 0.025), 
                    ls_pred = quantile(post_pred, 0.975), 
                    sites = site, 
                    month = month)
             
           })

est_rainfall <- do.call('rbind', est_rainfall)

est_rainfall$sites <- as.factor(est_rainfall$sites)
est_rainfall$sites <- factor(est_rainfall$sites, 
                             labels = levels(as.factor(rainfall_data$sites)))

dat_rainfall$level <- as.factor(dat_rainfall$level)
dat_rainfall$sites <- factor(dat_rainfall$level, 
                             labels = levels(as.factor(rainfall_data$sites)))



ggplot() +
  # geom_ribbon(data = est_rainfall,
  #             aes(month, ymin = li_pred, ymax = ls_pred, fill = sites),
  #             alpha = 0.2) +
  geom_ribbon(data = est_rainfall, 
              aes(month, ymin = li_mu, ymax = ls_mu, fill = sites), 
              alpha = 0.2) +
  geom_line(data = est_rainfall, 
            aes(month, mu, color = sites)) +
  geom_point(data = dat_rainfall, 
             aes(month, rainfall, color = sites), size = 0.8) +
  facet_wrap(~sites, scales = 'free_y') +
  scale_x_continuous(breaks = 1:12) +
  labs(y = 'Rainfall (mm)', x = 'Years') +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        text = element_text(size = 14), 
        legend.position = 'none')


dat_temperature <- 
  list(
    N = nrow(temp_data),
    M = 12,
    N_year = 3,
    L = max(as.numeric(as.factor(temp_data$sites))), 
    t = temp_data$month,
    year_id = as.numeric(as.factor(temp_data$year)),
    month = temp_data$month,
    level = as.numeric(as.factor(temp_data$sites)),
    temperature = temp_data$temperature,
    period = 12
  )


cat(file = 'climate_model.stan', 
    '
    data {
      int N;                // Total observations (252)
      int M;                // N months (12)
      int L;                // levels of factor (sites) (7)
      int N_year;           // N year
      array[N] int year_id;  // years
      array[N] int t;       // number of time points (1, ..., 36)
      array[N] int month;   // month indices (1, ..., 12)
      array[N] int level;   // factor (grouping variable) (1, ..., 7)
      vector[N] temperature;   // rainfall data
      int period;           // period of circular time (12)
    }
    
    parameters {
      real alpha;
      real<lower = 0> sigma;                // noise for the likelihood
      vector<lower = 0>[L] length_scale;    // length scale (smooth term)
      vector<lower = 0>[L] sigma_f;         // noise for GP
      array[L] vector[M] eta;    // latent variables for each month
      vector[N_year] z_year;
      real mu_year;
      real<lower = 0> sigma_year;
    }
    
    transformed parameters {
      array[L] vector[M] f;
      for (l in 1:L) {
        matrix[M, M] K;
        matrix[M, M] L_K;
        
        for (i in 1:(M-1)) {
          for (j in (i+1):M) {
            real distance = abs(i - j);
            real periodic_distance = fmin(distance, period - distance);
            // periodic kernel
            K[i, j] = sigma_f[l]^2 * exp(-2 * square(sin(pi()*periodic_distance / period))/square(length_scale[l]));
            K[j, i] = K[i, j];
          }
          K[i, i] = sigma_f[l]^2 + 1e-9;
        }
        K[M, M] = sigma_f[l]^2 + 1e-9;
        
        // cholesky decomposition
        L_K = cholesky_decompose(K);
        
        // transforme the latent variable eta to the GP
        f[l] = L_K * eta[l];
      }
    
      vector[N_year] year;
      year = mu_year + z_year * sigma_year;
    
    }
    
    model {
      sigma ~ cauchy(0, 1);
      length_scale ~ inv_gamma(5, 5);
      sigma_f ~ cauchy(0, 1);
      for (l in 1:L) {
        eta[l] ~ normal(0, 1);
      }
    
      alpha ~ normal(20, 5);
    
      z_year ~ normal(0, 1);
      mu_year ~ normal(0, 0.25);
      sigma_year ~ exponential(1);
    
      for (n in 1:N) {
        temperature[n] ~ student_t(7, alpha + f[level[n]][month[n]] + 
                                      year[year_id[n]], sigma);
      }
    }
    
    generated quantities {
      vector[N] y_pred;
      
      for (i in 1:N) {
        y_pred[i] = student_t_rng(7, alpha + f[level[i]][month[i]] + 
                                  year[year_id[i]], sigma);
      }
    }
    ')

file <- paste0(getwd(), '/climate_model.stan')

fit_climate <- cmdstan_model(file, compile = T)

mod_temperature <-
  fit_climate$sample(
    data = dat_temperature,
    chains = 3,
    parallel_chains = 3,
    iter_sampling = 1e4,
    iter_warmup = 1e3,
    thin = 10,
    seed = 123
  )

summary_temp <- mod_temperature$summary()

summary_temp %$% plot(rhat ~ ess_bulk)
summary_temp %$% points(rhat ~ ess_tail, col = 'red')




ppcheck_temp <- mod_temperature$draws('y_pred', format = 'matrix')

plot(density(dat_temperature$temperature), main = '', ylim = c(0, 0.2), 
     xlab = 'Average monthly temperature')
for (i in 1:400) lines(density(ppcheck_temp[i, ]), lwd = 0.1)
lines(density(dat_temperature$temperature), col = 'red', lwd = 2)



post_temp <- mod_temperature$draws(c('alpha', 'f', 'year', 'sigma'), format = 'df')

post_temp <- 
  list(f = post_temp[, grep('^f', colnames(post_temp))],
       year = post_temp[, grep('^year', colnames(post_temp))],
       alpha = post_temp$alpha,
       sigma = post_temp$sigma)

dat_temperature <- as_tibble(do.call('cbind', dat_temperature[1:7]))

est_temp2 <- 
  lapply(1:12, FUN = 
           function(x) {
             
             cc <- paste0('^(.*)(,', x, '\\])$')
             
             indx <- grep(cc, colnames(post_rainfall$f))
             
             post <- apply(post_temp$f[, indx, drop = T], 1, mean)
             a <- post_temp$alpha
             year <- apply(post_temp$year, 1, mean)
             sigma <- post_temp$sigma
             
             mu_est <- a + post + year
             
             tibble(mu = mean(mu_est) , 
                    li_mu = quantile(mu_est, 0.025), 
                    ls_mu = quantile(mu_est, 0.975),
                    month = x)
             
           })

est_temp2 <- do.call('rbind', est_temp2)

for (i in 2:3) {
  est_temp2[[i]] <- 
    (est_temp2[[i]] - mean(est_temp2$mu)) / 
    sd(est_temp2$mu)
}

est_temp2$mu <- as.vector(scale(est_temp2$mu))

est_rainfall2$class <- 'Rainfall'
est_temp2$class <- 'Temperature'

# saveRDS(rbind(est_rainfall2, 
#               est_temp2), 'climate_for_plot.rds')

est_temp <- 
  lapply(1:nrow(site_month), FUN = 
           function(x) {
             site <- site_month$level[x]
             month <- site_month$month[x]
             
             f_GP <- 
               paste0('f[', site, ',', month, ']')
             
             post <- post_temp$f[, f_GP, drop = T]
             sigma <- post_temp$sigma
             a <- post_temp$alpha
             year <- apply(post_temp$year, 1, mean)
             
             mu_est <- a + post + year
             
             set.seed(555)
             post_pred <- rstudent(7, length(post), mu_est, sigma)
             
             tibble(mu = mean(mu_est), 
                    li_mu = quantile(mu_est, 0.025), 
                    ls_mu = quantile(mu_est, 0.975), 
                    li_pred = quantile(post_pred, 0.025), 
                    ls_pred = quantile(post_pred, 0.975), 
                    sites = site, 
                    month = month)
             
           })

est_temp <- do.call('rbind', est_temp)

est_temp$sites <- as.factor(est_temp$sites)
est_temp$sites <- factor(est_temp$sites, 
                         labels = levels(as.factor(temp_data$sites)))

# saveRDS(list(rainfall = est_rainfall, 
#              temperature = est_temp), 'AVG_climate_data.rds')


ggplot() +
  # geom_ribbon(data = est_temp, 
  #             aes(month, ymin = li_pred, ymax = ls_pred, fill = sites), 
  #             alpha = 0.2) +
  geom_ribbon(data = est_temp, 
              aes(month, ymin = li_mu, ymax = ls_mu, fill = sites), 
              alpha = 0.2) +
  geom_line(data = est_temp, 
            aes(month, mu, color = sites)) +
  geom_point(data = temp_data, 
             aes(month, temperature, color = sites), size = 0.8) +
  facet_wrap(~sites, scales = 'free') +
  labs(y = 'Temperature (°C)', x = 'Years') +
  scale_x_continuous(breaks = 1:12) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        text = element_text(size = 14), 
        legend.position = 'none')


sessionInfo()
