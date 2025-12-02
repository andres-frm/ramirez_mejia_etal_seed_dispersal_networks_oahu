
sapply(c('dplyr', 'ggplot2', 'lubridate', 'forcats', 
         'magrittr', 'cmdstanr', 'rethinking', 'cowplot', 
         'readxl', 'patchwork', 'tidyr'), 
       library, character.only = T)

extrafont::loadfonts(device = 'win')

source('functions_mod_diagnostics.r')


d <- read_xlsx('SeedTraps.xlsx', col_names = T, sheet = 1)

d <- d[, c("TrapID", "Site", "Date", "Month", "MonthJef", 
           "Year", "Species", '#Seeds', "#seed/#seed per fruit", 
           "#fruitsTotal", "TrapTipped", "Dispersed? (Y/N)")]

d <- d[d$TrapTipped != 'Y', ]

d$`#seed/#seed per fruit`[which(d$`#seed/#seed per fruit` == 'NA')] <- 0

d$Date <- d %$% paste(day(Date), MonthJef, Year, sep = '-')
d <- d[-grep('NA', d$Date), ]

d$Date <- dmy(d$Date)

sampling_effort <- 
  d[, c("TrapID", "Site", "Date")] |> 
  unique() |> 
  group_by(Date, Site) |> 
  transmute(SampEff = length(TrapID)) |> 
  unique()

sampling_effort$month <- month(sampling_effort$Date)

sampling_effort <- 
  sampling_effort |> 
  group_by(Site, month) |> 
  transmute(SampEff = mean(SampEff)) |> 
  unique()

d$Species[which(d$`Dispersed? (Y/N)` == 'N' |
                  is.na(d$`Dispersed? (Y/N)`))] <- 'NA'

d <- 
  d[, c("TrapID", "Site", "Date", "Species", 
      "#seed/#seed per fruit", "#fruitsTotal")] |> 
  unique()

d$month <- month(d$Date)

d$date2 <- d %$% paste(year(Date), month, sep = '-')

# fruit production data 

fruit_abundance <- 
  d |> 
  group_by(Site, date2) |> 
  transmute(abun = sum(`#fruitsTotal`), 
            month = month) |> 
  ungroup() |> 
  group_by(Site, month) |> 
  transmute(abun = mean(abun)) |> 
  unique()

# richness of fruit production fruits
S_fruiting <- 
  d |> 
  filter(Species != 'NA') |> 
  select(Site, date2, Species, month) |> 
  unique() |> 
  group_by(Site, date2) |> 
  transmute(S = length(Species), 
            month = month) |> 
  unique() |> 
  group_by(Site, month) |> 
  transmute(S = mean(S)) |> 
  unique()

phenology <- full_join(S_fruiting, fruit_abundance, by = c('Site', 'month'))

phenology <- full_join(phenology, sampling_effort, by = c('Site', 'month'))

phenology$SampEff_std <- phenology$SampEff / mean(phenology$SampEff)

phenology$Site <- as.factor(phenology$Site)

dist_mat <- readRDS('mat_distance_dites.rds')

dat <- lapply(phenology, function(x) x)

dat$Site <- as.numeric(dat$Site)

dat$abun <- log(dat$abun)

dat$dist_sites <- dist_mat

dat$S <- round(dat$S)

dat$N <- length(dat$Site)
dat$N_sites <- max(dat$Site)
dat$N_months <- 12



cat(file = 'S_phenology.stan', 
    '
    functions {
  vector GP_periodic(int period,   // periodicity 
                     real gamma,   // smoothing term of the GP
                     real sigma,   // noise parameter of the GP
                     vector eta){  // latent variable per month
                       
                       int M = period;
                       matrix[M, M] K;
                       matrix[M, M] L_K;
                       
                       for (i in 1:(M-1)) {
                         for (j in (i+1):M) {
                           
                           real distance = abs(i - j);
                           real periodic_distance = fmin(distance, period - distance);
                           K[i, j] = sigma^2 * exp(-2*square(sin(pi()*periodic_distance/period))/gamma^2);
                           K[j, i] = K[i, j]; // filling the lower triangle 
                         }
                         K[i, i] = sigma^2 + 1e-9; // small values at the diagonal for stability
                       }
                       
                       K[M, M] = sigma^2 + 1e-9; // small values at the begining and end of the matrix
                       return(cholesky_decompose(K) * eta);
                     }
                     
  matrix GP_quadratic(matrix x,
                      real eta,
                      real rho, 
                      real delta) {
                        
                        int N = dims(x)[1];
                        matrix[N, N] K;
                        matrix[N, N] L_K;
                        
                        for (i in 1:(N-1)) {
                          K[i, i] = eta + delta; // small values to the diagonal
                          for (j in (i+1):N) {
                            K[i, j] = eta * exp(-rho * square(x[i, j]));
                            K[j, i] = K[i, j]; // filling low triangle
                          }
                        }
                      
                      K[N, N] = eta + delta;
                      L_K = cholesky_decompose(K);
                      return L_K;
                      }
}

data {
  int N;
  int N_sites;
  int N_months;
  array[N] int Site;
  array[N] int month;
  array[N] int S;
  vector[N] abun;
  vector[N] SampEff_std;
  matrix[N_sites, N_sites] dist_sites;
}

parameters {
  
  real alpha;
  //real<lower = 0> sigma;
  
  // periodic GP
  vector<lower = 0>[N_sites] gamma_f;
  vector<lower = 0>[N_sites] sigma_f;
  matrix[N_months, N_sites] eta_f;
  
  // Quadratic GP
  vector[N_sites] z_sites;
  real<lower = 0> eta;
  real<lower = 0> rho;
  
}

transformed parameters {
  // periodic GP
  matrix[N_months, N_sites] f;
  
  for (i in 1:N_sites) {
    f[, i] = GP_periodic(
      12, 
      gamma_f[i],
      sigma_f[i],
      eta_f[,i]
    );
  }
  
  // Quadratic GP
  
  vector[N_sites] theta;
  matrix[N_sites, N_sites] L_K_theta;
  L_K_theta = GP_quadratic(dist_sites, 
                           eta, 
                           rho, 
                           0.001);
  theta = L_K_theta * z_sites;
}

model {
  alpha ~ normal(1.5, 0.25);
  //sigma ~ exponential(1);
  
  to_vector(eta_f) ~ normal(1.5, 0.25);
  gamma_f ~ inv_gamma(5, 5);
  sigma_f ~ cauchy(0, 1);
  
  z_sites ~ normal(0, 1);
  eta ~ exponential(4); 
  rho ~ exponential(1);
  
  // likelihood
  
  for (i in 1:N) {
   S[i] ~ poisson(exp(alpha + f[month[i], Site[i]] + 
                           theta[Site[i]] + log(SampEff_std[i]))); 
  }
}

generated quantities {
  vector[N] mu;
  array[N] int ppcheck;
  
  for (i in 1:N) {
    mu[i] = exp(alpha + f[month[i], Site[i]] + theta[Site[i]]) .*SampEff_std[i];
  }
  ppcheck = poisson_rng(mu);
  
}
    ')

file <- paste0(getwd(), '/S_phenology.stan')

fit_S <- cmdstan_model(file, compile = T)

mod_S <- 
  fit_S$sample(
    data = dat, 
    iter_sampling = 2e3, 
    iter_warmup = 500, 
    chains = 3, 
    parallel_chains = 3, 
    thin = 3, 
    seed = 5
  )



summary_S <- mod_S$summary()

mod_diagnostics(mod_S, summary_S)



ppcheck_S <- mod_S$draws('ppcheck', format = 'matrix')
plot(density(dat$S), ylim = c(0, 0.20), xlim = c(-2, 20), 
     xlab = 'S (fruiting plants)', main = '')
for (i in 1:200) lines(density(ppcheck_S[i, ]), lwd = 0.1) 
lines(density(dat$S), col = 'red', lwd = 2)



mu_S <- mod_S$draws('mu', format = 'df')

mu_S <- mu_S[, grep('mu', colnames(mu_S))]

mu_S <- lapply(mu_S, FUN = 
                 function(x) {
                   tibble(mu = mean(x), 
                          li = quantile(x, 0.025), 
                          ls = quantile(x, 0.975))
                 })

mu_S <- do.call('rbind', mu_S)

colnames(mu_S) <- paste0('S_', colnames(mu_S))

phenology <- cbind(phenology, mu_S)

phenology$date <- dmy(paste0('01-', phenology$month, '-2025'))

inter_data <- readRDS('data_for_models.rds')
inter_data <- unique(inter_data[, c("site2", "month")])
inter_data$cod <- inter_data %$% paste0(site2, '_', month)

phenology$cod <- phenology %$% paste0(Site, '_', month)

impute_dates <- !(inter_data$cod %in% phenology$cod)

inter_data$site2 <- as.factor(inter_data$site2)
inter_data$site <- as.numeric(inter_data$site2)

post_S <- mod_S$draws(c('alpha', 'f', 'theta'), format = 'df')
post_S <- lapply(c('alpha', 'f', 'theta'), FUN = 
                   function(x) {
                     post_S[, grep(x, colnames(post_S))]
                   })

names(post_S) <- c('alpha', 'f', 'theta')

pred_S <- 
  lapply(1:7, FUN = 
           function(i) {
             
             indx <- which(i == dat$Site)[1]
             S <- phenology$Site[indx]
             
             site <- 
               lapply(1:12, FUN = 
                        function(j) {
                          
                          GP <- paste0('f[', j, ',', i, ']')
                          t <- paste0('theta[', i, ']')
                          
                          cod <- paste0(S, '_', j)
                          indx_SEff <- which(cod == phenology$cod)
                          indx_months <- which(i == phenology$month)
                          
                          est <- 
                            with(post_S, 
                                 {
                                   alpha[, 1, drop = T] +
                                     f[, GP, drop = T] +
                                     theta[, t, drop = T]
                                 })
                          
                          if (length(indx_SEff) > 0) {
                            est <- exp(est) * dat$SampEff_std[indx_SEff]
                          } else {
                            est <- exp(est) * median(dat$SampEff_std[indx_months])
                          }
                          
                          est <- 
                            tibble(Site = i, 
                                   date = dmy(paste0('01-', j, '-2025')), 
                                   S_mu = sample(est, 1e3))
                          est
                        })
             
             site <- do.call('rbind', site)
             
             site$Site <- S
             site
           })

pred_S <- do.call('rbind', pred_S)

pred_S$z_s <- (pred_S$S_mu - mean(dat$S)) / sd(dat$S)

est_S <- 
  lapply(1:7, FUN = 
         function(i) {
           
           indx <- which(i == dat$Site)[1]
           S <- phenology$Site[indx]
           
           site <- 
             lapply(1:12, FUN = 
                      function(j) {
                        
                        GP <- paste0('f[', j, ',', i, ']')
                        t <- paste0('theta[', i, ']')
                        
                        cod <- paste0(S, '_', j)
                        indx_SEff <- which(cod == phenology$cod)
                        indx_months <- which(i == phenology$month)
                        
                        est <- 
                          with(post_S, 
                               {
                                 alpha[, 1, drop = T] +
                                   f[, GP, drop = T] +
                                   theta[, t, drop = T]
                               })
                        
                        if (length(indx_SEff) > 0) {
                          est <- exp(est) * dat$SampEff_std[indx_SEff]
                        } else {
                          est <- exp(est) * median(dat$SampEff_std[indx_months])
                        }
                        
                        est <- 
                          tibble(Site = i, 
                                 date = dmy(paste0('01-', j, '-2025')), 
                                 S_mu = mean(est), 
                                 S_li = quantile(est, 0.025), 
                                 S_ls = quantile(est, 0.975))
                        est
                      })
           
           site <- do.call('rbind', site)
           
           site$Site <- S
           site
         })

est_S <- do.call('rbind', est_S)



ggplot() +
  geom_line(data = phenology, aes(date, S_mu, color = 'No imputation')) +
  geom_point(data = phenology, aes(date, S_mu, color = 'No imputation'), 
             size = 1.5) +
  geom_ribbon(data = phenology, 
              aes(date, S_mu, ymin = S_li, ymax = S_ls, 
                  fill = 'No imputation'), 
              alpha = 0.3) +
  geom_line(data = est_S, aes(date, S_mu, color = 'Imputation'), lwd = 0.1) +
  geom_point(data = est_S, aes(date, S_mu, color = 'Imputation'), size = 0.2) +
  geom_ribbon(data = est_S, 
              aes(date, S_mu, ymin = S_li, ymax = S_ls, fill = 'Imputation'), 
              alpha = 0.3) +
  scale_x_date(date_labels = '%b', date_breaks = '2 month') + 
  facet_wrap(~Site) +
  labs(y = 'S (fruiting plants)') +
  theme_bw()

cat(file = 'abun_phenology.stan', 
    '
    functions {
  vector GP_periodic(int period,   // periodicity 
                     real gamma,   // smoothing term of the GP
                     real sigma,   // noise parameter of the GP
                     vector eta){  // latent variable per month
                       
                       int M = period;
                       matrix[M, M] K;
                       matrix[M, M] L_K;
                       
                       for (i in 1:(M-1)) {
                         for (j in (i+1):M) {
                           
                           real distance = abs(i - j);
                           real periodic_distance = fmin(distance, period - distance);
                           K[i, j] = sigma^2 * exp(-2*square(sin(pi()*periodic_distance/period))/gamma^2);
                           K[j, i] = K[i, j]; // filling the lower triangle 
                         }
                         K[i, i] = sigma^2 + 1e-9; // small values at the diagonal for stability
                       }
                       
                       K[M, M] = sigma^2 + 1e-9; // small values at the begining and end of the matrix
                       return(cholesky_decompose(K) * eta);
                     }
                     
  matrix GP_quadratic(matrix x,
                      real eta,
                      real rho, 
                      real delta) {
                        
                        int N = dims(x)[1];
                        matrix[N, N] K;
                        matrix[N, N] L_K;
                        
                        for (i in 1:(N-1)) {
                          K[i, i] = eta + delta; // small values to the diagonal
                          for (j in (i+1):N) {
                            K[i, j] = eta * exp(-rho * square(x[i, j]));
                            K[j, i] = K[i, j]; // filling low triangle
                          }
                        }
                      
                      K[N, N] = eta + delta;
                      L_K = cholesky_decompose(K);
                      return L_K;
                      }
}

data {
  int N;
  int N_sites;
  int N_months;
  array[N] int Site;
  array[N] int month;
  array[N] int S;
  vector[N] abun;
  vector[N] SampEff_std;
  matrix[N_sites, N_sites] dist_sites;
}

parameters {
  
  real alpha;
  real<lower = 0> sigma;
  
  // periodic GP
  vector<lower = 0>[N_sites] gamma_f;
  vector<lower = 0>[N_sites] sigma_f;
  matrix[N_months, N_sites] eta_f;
  
  // Quadratic GP
  vector[N_sites] z_sites;
  real<lower = 0> eta;
  real<lower = 0> rho;
  
}

transformed parameters {
  // periodic GP
  matrix[N_months, N_sites] f;
  
  for (i in 1:N_sites) {
    f[, i] = GP_periodic(
      12, 
      gamma_f[i],
      sigma_f[i],
      eta_f[,i]
    );
  }
  
  // Quadratic GP
  
  vector[N_sites] theta;
  matrix[N_sites, N_sites] L_K_theta;
  L_K_theta = GP_quadratic(dist_sites, 
                           eta, 
                           rho, 
                           0.001);
  theta = L_K_theta * z_sites;
}

model {
  alpha ~ normal(0, 1);
  sigma ~ exponential(1);
  
  to_vector(eta_f) ~ normal(0, 0.5);
  gamma_f ~ inv_gamma(5, 5);
  sigma_f ~ cauchy(0, 1);
  
  z_sites ~ normal(0, 1);
  eta ~ exponential(4); 
  rho ~ exponential(1);
  
  // likelihood
  
  for (i in 1:N) {
   abun[i] ~ student_t(6, alpha + f[month[i], Site[i]] + 
                          theta[Site[i]] + SampEff_std[i], 
                       sigma); 
  }
}

generated quantities {
  vector[N] mu;
  array[N] real ppcheck;
  
  for (i in 1:N) {
    mu[i] = alpha + f[month[i], Site[i]] + 
            theta[Site[i]] + SampEff_std[i];
  }
  ppcheck = student_t_rng(7, mu, sigma);
  
}

    ')

dat$abun <- as.vector(scale(dat$abun))

file <- paste0(getwd(), '/abun_phenology.stan')

fit_ABU <- cmdstan_model(file, compile = T)

mod_ABU <- 
  fit_ABU$sample(
    data = dat, 
    iter_sampling = 2e3, 
    iter_warmup = 500, 
    chains = 3, 
    parallel_chains = 3, 
    thin = 3, 
    seed = 5
  )



summary_ABU <- mod_ABU$summary()

mod_diagnostics(mod_ABU, summary_ABU)



ppcheck_ABU <- mod_ABU$draws('ppcheck', format = 'matrix')
plot(density(dat$abun), ylim = c(0, 0.6), xlim = c(-5, 5), 
     xlab = 'Fruit abundance', main = '')
for (i in 1:200) lines(density(ppcheck_ABU[i, ]), lwd = 0.1) 
lines(density(dat$abun), col = 'red', lwd = 2)



mu_ABU <- mod_ABU$draws('mu', format = 'df')

mu_ABU <- mu_ABU[, grep('mu', colnames(mu_ABU))]

mu_ABU <- lapply(mu_ABU, FUN = 
                 function(x) {
                   tibble(mu = mean(x), 
                          li = quantile(x, 0.025), 
                          ls = quantile(x, 0.975))
                 })

mu_ABU <- do.call('rbind', mu_ABU)

colnames(mu_ABU) <- paste0('ABU_', colnames(mu_ABU))

phenology <- cbind(phenology, mu_ABU)

post_ABU <- mod_ABU$draws(c('alpha', 'f', 'theta'), format = 'df')
post_ABU <- lapply(c('alpha', 'f', 'theta'), FUN = 
                   function(x) {
                     post_ABU[, grep(x, colnames(post_ABU))]
                   })

names(post_ABU) <- c('alpha', 'f', 'theta')

pred_ABU <- 
  lapply(1:7, FUN = 
           function(i) {
             
             indx <- which(i == dat$Site)[1]
             S <- phenology$Site[indx]
             
             site <- 
               lapply(1:12, FUN = 
                        function(j) {
                          
                          GP <- paste0('f[', j, ',', i, ']')
                          t <- paste0('theta[', i, ']')
                          
                          cod <- paste0(S, '_', j)
                          indx_SEff <- which(cod == phenology$cod)
                          indx_months <- which(i == phenology$month)
                          
                          est <- 
                            with(post_ABU, 
                                 {
                                   alpha[, 1, drop = T] +
                                     f[, GP, drop = T] +
                                     theta[, t, drop = T]
                                 })
                          
                          if (length(indx_SEff) > 0) {
                            est <- est + dat$SampEff_std[indx_SEff]
                          } else {
                            est <- est + median(dat$SampEff_std[indx_months])
                          }
                          
                          est <- 
                            tibble(Site = i, 
                                   date = dmy(paste0('01-', j, '-2025')), 
                                   ABU_mu = sample(est, 1e3))
                          est
                        })
             
             site <- do.call('rbind', site)
             
             site$Site <- S
             site
           })

pred_ABU <- do.call('rbind', pred_ABU)

est_ABU <- 
  lapply(1:7, FUN = 
           function(i) {
             
             indx <- which(i == dat$Site)[1]
             S <- phenology$Site[indx]
             
             site <- 
               lapply(1:12, FUN = 
                        function(j) {
                          
                          GP <- paste0('f[', j, ',', i, ']')
                          t <- paste0('theta[', i, ']')
                          
                          cod <- paste0(S, '_', j)
                          indx_SEff <- which(cod == phenology$cod)
                          indx_months <- which(i == phenology$month)
                          
                          est <- 
                            with(post_ABU, 
                                 {
                                   alpha[, 1, drop = T] +
                                     f[, GP, drop = T] +
                                     theta[, t, drop = T]
                                 })
                          
                          if (length(indx_SEff) > 0) {
                            est <- est + dat$SampEff_std[indx_SEff]
                          } else {
                            est <- est + median(dat$SampEff_std[indx_months])
                          }
                          
                          est <- 
                            tibble(Site = i, 
                                   date = dmy(paste0('01-', j, '-2025')), 
                                   ABU_mu = mean(est), 
                                   ABU_li = quantile(est, 0.025), 
                                   ABU_ls = quantile(est, 0.975))
                          est
                        })
             
             site <- do.call('rbind', site)
             
             site$Site <- S
             site
           })

est_ABU <- do.call('rbind', est_ABU)


ggplot() +
  geom_line(data = phenology, aes(date, ABU_mu, color = 'No imputation')) +
  geom_point(data = phenology, aes(date, ABU_mu, color = 'No imputation'), 
             size = 1.5) +
  geom_ribbon(data = phenology, 
              aes(date, ABU_mu, ymin = ABU_li, ymax = ABU_ls, 
                  fill = 'No imputation'), 
              alpha = 0.3) +
  geom_line(data = est_ABU, aes(date, ABU_mu, color = 'Imputation'), lwd = 0.1) +
  geom_point(data = est_ABU, aes(date, ABU_mu, color = 'Imputation'), size = 0.2) +
  geom_ribbon(data = est_ABU, 
              aes(date, ABU_mu, ymin = ABU_li, ymax = ABU_ls, fill = 'Imputation'), 
              alpha = 0.3) +
  scale_x_date(date_labels = '%b', date_breaks = '2 month') + 
  facet_wrap(~Site) +
  labs(y = 'Fruit abundance') +
  theme_bw()


est_ABU$month <- month(est_ABU$date)
est_S$month <- month(est_S$date)

# data for global models

# saveRDS(list(abun_phenology = est_ABU, 
#              S_phenology = est_S), 'phenology_data.rds')
# 
# 
# saveRDS(list(abun_phenology = pred_ABU, 
#              S_phenology = pred_S), 'phenology_data_PREDICTED.rds')

sessionInfo()
