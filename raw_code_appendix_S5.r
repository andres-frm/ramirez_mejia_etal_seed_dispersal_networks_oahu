sapply(c('dplyr', 'ggplot2', 'lubridate', 'forcats', 
         'magrittr', 'cmdstanr', 'rethinking', 'cowplot', 
         'bipartite', 'igraph', 'iNEXT', 'hypervolume', 
         'readxl', 'patchwork', 'ggfortify', 'tidyr', 'gtools', 
         'ggsankey', 'ggridges'), 
       library, character.only = T)

extrafont::loadfonts(device = 'win')

source('functions_mod_diagnostics.r')



phenology <- readRDS('phenology_data.rds')

alien_native <- readRDS('alien_native_prop.rds')

colnames(alien_native)[1:2] <- c('prop_alien', 'site2')

dist_mat <- readRDS('mat_distance_dites.rds')

d <- readRDS('data_for_models.rds')

d <- full_join(d, alien_native, by = c('site2', 'month'))

d <- d[, -8]

colnames(d)[9] <- 'site'

phenology <- lapply(phenology, 
                    function(x) {
                      colnames(x) <- tolower(colnames(x))
                      x
                    })

d <- 
  left_join(d, 
          phenology$abun_phenology[, c('site', 'month', 'abu_mu')], 
          by = c('site', 'month'))

d <- 
  left_join(d, 
            phenology$S_phenology[, c('site', 'month', 's_mu')], 
            by = c('site', 'month'))

colnames(d)[grep('mu', colnames(d))] <- c('fruit_abun', 'fruit_S')

dat <- lapply(d, function(x) x)

dat$site <- as.numeric(as.factor(dat$site))

dat$HV_network <- as.vector(scale(sqrt(dat$HV_network)))
dat$size <- as.vector(scale(log(dat$size)))
dat$fruit_S <- as.vector(scale(dat$fruit_S))
dat$HV_bird <- as.vector(scale(dat$HV_bird))
dat$HV_plants_nut <- as.vector(scale(dat$HV_plants_nut)) 
dat$HV_plants_morfo <- as.vector(scale(sqrt(dat$HV_plants_morfo)))
dat$HV_plant <- as.vector(scale(dat$HV_plant))
dat$nestedness <- as.vector(scale(dat$nestedness))

dat$H2[which(dat$H2 >= 0.99)] <- 0.9

names(dat)[1] <- 'network_size'

dat$N <- length(dat$network_size)
dat$N_sites <- max(dat$site)
dat$N_month <- max(dat$month)
dat$dist_sites <- dist_mat

codes_sites <- 
  unique(tibble(cod = dat$site, 
                label = d$site))



cat(file = 'Nestedness_global.stan',
    "

    functions {
      vector GP_periodic(int period,        // periodicity
                         real gamma,        // smoothing term of the GP
                         real sigma,        // noise paramether
                         vector eta) {      // latent variable for each month

                         int M = period;
                         matrix[M, M] K;
                         matrix[M, M] L_K;

                         for (i in 1:(M - 1)) {
                           for (j in (i+1):M) {
                               real distance = abs(i - j);
                               real periodic_distance = fmin(distance, period - distance);
                               K[i, j] = sigma^2 * exp(-2 * square(sin(pi()*periodic_distance/period))/gamma^2);
                               K[j, i] = K[i, j];   // filling the lower triangle
                            }
                             K[i, i] = sigma^2 + 1e-9;   // small values to guarante stability
                          }
                           K[M, M] = sigma^2 + 1e-9;    // small values to guarante stability
                           return cholesky_decompose(K) * eta;

                        }

      matrix GP_quadratic(matrix x,
                          real eta,
                          real rho,
                          real delta) {

                          int N = dims(x)[1];
                          matrix[N, N] K;
                          matrix[N, N] L_K;

                          for (i in 1:(N-1)) {
                            K[i, i] = eta + delta;
                            for (j in (i+1):N) {
                              K[i, j] = square(eta) * exp(-rho * square(x[i, j]));
                              K[j, i] = K[i, j];
                            }
                          }

                          K[N, N] = eta + delta;
                          L_K = cholesky_decompose(K);
                          return L_K;
                          }
    }

    data{
      int N;
      int N_sites;
      int N_month;
      array[N] int network_size;
      vector[N] modularity;
      vector[N] nestedness;
      vector[N] H2;
      array[N] int month;
      array[N] int site;
      vector[N] HV_network;
      vector[N] HV_plant;
      vector[N] HV_plants_morfo;
      vector[N] HV_plants_nut;
      vector[N] HV_bird;
      vector[N] z_temperature;
      vector[N] z_rainfall;
      vector[N] fruit_S;
      vector[N] fruit_abun;
      matrix[N_sites, N_sites] dist_sites;
    }

    parameters {


      //////////////////////////// Fruit abundance //////////////////////
      ////////////////////////////////////////////////////////////////

      // GP periodic
      real<lower = 0> gamma_fruitS;
      real<lower = 0> sigma_f_fruitS;
      vector[N_month] eta_fruitS;
      // GP quadratic
      vector[N_sites] z_sites_fruitS;
      real<lower = 0> eta_site_fruitS;
      real<lower = 0> rho_site_fruitS;

      // Pars linear model
      real beta_temp_fruitS;
      real beta_rain_fruitS;
      // noise for likelihood function
      real<lower = 0> sigma_fruitS;


      //////////////////////////// Fruit abundance //////////////////////
      ////////////////////////////////////////////////////////////////

      // GP periodic
      real<lower = 0> gamma_fruitAB;
      real<lower = 0> sigma_f_fruitAB;
      vector[N_month] eta_fruitAB;
      // GP quadratic
      vector[N_sites] z_sites_fruitAB;
      real<lower = 0> eta_site_fruitAB;
      real<lower = 0> rho_site_fruitAB;

      // Pars linear model

      real beta_temp_fruitAB;
      real beta_rain_fruitAB;
      // noise for likelihood function
      real<lower = 0> sigma_fruitAB;


      //////////////////////////// Fru. Nut. HV //////////////////////
      ////////////////////////////////////////////////////////////////

      // GP periodic
      real<lower = 0> gamma_nut_HV;
      real<lower = 0> sigma_f_nut_HV;
      vector[N_month] eta_nut_HV;
      // GP quadratic
      vector[N_sites] z_sites_nut_HV;
      real<lower = 0> eta_site_nut_hv;
      real<lower = 0> rho_site_nut_hv;

      // Pars linear model

      real beta_temp_nut_HV;
      real beta_rain_nut_HV;
      real beta_fruitAB_nut_HV;
      real beta_fruitS_nut_HV;
      real beta_BIRD_nut_HV;
      // noise for likelihood function
      real<lower = 0> sigma_nut_HV;

      //////////////////////////// Fru. morpho. HV //////////////////////
      ////////////////////////////////////////////////////////////////
      // GP periodic
      real<lower = 0> gamma_mor_HV;
      real<lower = 0> sigma_f_mor_HV;
      vector[N_month] eta_mor_HV;
      // GP quadratic
      vector[N_sites] z_sites_mor_HV;
      real<lower = 0> eta_site_mor_hv;
      real<lower = 0> rho_site_mor_hv;

      // Pars linear model

      real beta_temp_mor_HV;
      real beta_rain_mor_HV;
      real beta_fruitAB_mor_HV;
      real beta_fruitS_mor_HV;
      real beta_BIRD_mor_HV;
      // noise for likelihood function
      real<lower = 0> sigma_mor_HV;
      real nu;

      // //////////////////////////// Plants HV //////////////////////
      // ////////////////////////////////////////////////////////////////
      // GP periodic
      real<lower = 0> gamma_PLAN_HV;
      real<lower = 0> sigma_f_PLAN_HV;
      vector[N_month] eta_PLAN_HV;
      // GP quadratic
      vector[N_sites] z_sites_PLAN_HV;
      real<lower = 0> eta_site_PLAN_hv;
      real<lower = 0> rho_site_PLAN_hv;

      // Pars linear model
      real beta_temp_PLAN_HV;
      real beta_rain_PLAN_HV;
      real beta_fruitAB_PLAN_HV;
      real beta_fruitS_PLAN_HV;
      real beta_nut_PLAN_HV;
      real beta_mor_PLAN_HV;
      real beta_BIRD_PLAN_HV;
      // noise for likelihood function
      real<lower = 0> sigma_PLAN_HV;

      //////////////////////////// Birsds HV //////////////////////
      // ////////////////////////////////////////////////////////////////
      // GP periodic
      real<lower = 0> gamma_BIRD_HV;
      real<lower = 0> sigma_f_BIRD_HV;
      vector[N_month] eta_BIRD_HV;
      // GP quadratic
      vector[N_sites] z_sites_BIRD_HV;
      real<lower = 0> eta_site_BIRD_hv;
      real<lower = 0> rho_site_BIRD_hv;

      // Pars linear model
      real beta_temp_BIRD_HV;
      real beta_rain_BIRD_HV;
      real beta_fruitAB_BIRD_HV;
      real beta_fruitS_BIRD_HV;
      // noise for likelihood function
      real<lower = 0> sigma_BIRD_HV;

      //////////////////////////// Network metrics //////////////////////
      ////////////////////////////////////////////////////////////////
      // //////////////////////////// Nestedness //////////////////////
      // ////////////////////////////////////////////////////////////////
      // GP periodic
      real<lower = 0> gamma_NSS;
      real<lower = 0> sigma_f_NSS;
      vector[N_month] eta_NSS;
      // GP quadratic
      vector[N_sites] z_sites_NSS;
      real<lower = 0> eta_site_NSS;
      real<lower = 0> rho_site_NSS;

      // Pars linear model
      real beta_temp_NSS;
      real beta_rain_NSS;
      real beta_fruitAB_NSS;
      real beta_fruitS_NSS;
      real beta_nut_NSS;
      real beta_mor_NSS;
      real beta_PLAN_NSS;
      real beta_BIRD_NSS;
      // vector[N_sites] beta_NW_NSS;
      // noise for likelihood function
      real<lower = 0> sigma_NSS;

    }

    transformed parameters {

    //////////////////////////// Fruit abundance //////////////////////
    ////////////////////////////////////////////////////////////////
      // GP periodic
      vector[N_month] f_fruitAB;

      f_fruitAB = GP_periodic(
              12,
              gamma_fruitAB,
              sigma_f_fruitAB,
              eta_fruitAB
        );

      // GP quadratic
      vector[N_sites] alpha_fruitAB;
      matrix[N_sites, N_sites] L_K_fruitAB;
      L_K_fruitAB = GP_quadratic(dist_sites,
                                eta_site_fruitAB,
                                rho_site_fruitAB, 0.001);
      alpha_fruitAB = L_K_fruitAB * z_sites_fruitAB;


    //////////////////////////// Fruit richness //////////////////////
    ////////////////////////////////////////////////////////////////
      // GP periodic
      vector[N_month] f_fruitS;

      f_fruitS = GP_periodic(
              12,
              gamma_fruitS,
              sigma_f_fruitS,
              eta_fruitS
        );

      // GP quadratic
      vector[N_sites] alpha_fruitS;
      matrix[N_sites, N_sites] L_K_fruitS;
      L_K_fruitS = GP_quadratic(dist_sites,
                                eta_site_fruitS,
                                rho_site_fruitS, 0.001);
      alpha_fruitS = L_K_fruitS * z_sites_fruitS;

    //////////////////////////// Fru. Nut. HV //////////////////////
    ////////////////////////////////////////////////////////////////
      // GP periodic
      vector[N_month] f_nut_HV;

      f_nut_HV = GP_periodic(
              12,
              gamma_nut_HV,
              sigma_f_nut_HV,
              eta_nut_HV
        );

      // GP quadratic
      vector[N_sites] alpha_nut_HV;
      matrix[N_sites, N_sites] L_K_nut_HV;
      L_K_nut_HV = GP_quadratic(dist_sites,
                                eta_site_nut_hv,
                                rho_site_nut_hv, 0.001);
      alpha_nut_HV = L_K_nut_HV * z_sites_nut_HV;

    //   //////////////////////////// Fru. morpho. HV //////////////////////
    //   ////////////////////////////////////////////////////////////////
      // GP periodic
      vector[N_month] f_mor_HV;

      f_mor_HV = GP_periodic(
              12,
              gamma_mor_HV,
              sigma_f_mor_HV,
              eta_mor_HV
        );

      // GP quadratic
      vector[N_sites] alpha_mor_HV;
      matrix[N_sites, N_sites] L_K_mor_HV;
      L_K_mor_HV = GP_quadratic(dist_sites,
                                eta_site_mor_hv,
                                rho_site_mor_hv, 0.001);
      alpha_mor_HV = L_K_mor_HV * z_sites_mor_HV;
    //
    //   //////////////////////////// Plants HV //////////////////////
    //   ////////////////////////////////////////////////////////////////
      // GP periodic
      vector[N_month] f_PLAN_HV;

      f_PLAN_HV = GP_periodic(
              12,
              gamma_PLAN_HV,
              sigma_f_PLAN_HV,
              eta_PLAN_HV
        );

      // GP quadratic
      vector[N_sites] alpha_PLAN_HV;
      matrix[N_sites, N_sites] L_K_PLAN_HV;
      L_K_PLAN_HV = GP_quadratic(dist_sites,
                                eta_site_PLAN_hv,
                                rho_site_PLAN_hv, 0.001);
      alpha_PLAN_HV = L_K_PLAN_HV * z_sites_PLAN_HV;

      ////////////////////////// Birsds HV //////////////////////
      //////////////////////////////////////////////////////////////
      // GP periodic
      vector[N_month] f_BIRD_HV;

      f_BIRD_HV = GP_periodic(
              12,
              gamma_BIRD_HV,
              sigma_f_BIRD_HV,
              eta_BIRD_HV
        );

      // GP quadratic
      vector[N_sites] alpha_BIRD_HV;
      matrix[N_sites, N_sites] L_K_BIRD_HV;
      L_K_BIRD_HV = GP_quadratic(dist_sites,
                                eta_site_BIRD_hv,
                                rho_site_BIRD_hv, 0.001);
      alpha_BIRD_HV = L_K_BIRD_HV * z_sites_BIRD_HV;

      // //////////////////////////// Nestedness //////////////////////
      // ////////////////////////////////////////////////////////////////
      // GP periodic
      vector[N_month] f_NSS;

      f_NSS = GP_periodic(
              12,
              gamma_NSS,
              sigma_f_NSS,
              eta_NSS
        );

      // GP quadratic
      vector[N_sites] alpha_NSS;
      matrix[N_sites, N_sites] L_K_NSS;
      L_K_NSS = GP_quadratic(dist_sites,
                              eta_site_NSS,
                              rho_site_NSS, 0.001);
      alpha_NSS = L_K_NSS * z_sites_NSS;

    }

    model {

       //////////////////////////// Fruit abundance //////////////////////
       ////////////////////////////////////////////////////////////////
       // priors for periodic GP
       eta_fruitAB ~ normal(0, 0.5);
       gamma_fruitAB ~ inv_gamma(5, 5);
       sigma_f_fruitAB ~ cauchy(0, 1);
       // GP quadratic
       z_sites_fruitAB ~ normal(0, 0.25);
       eta_site_fruitAB ~ exponential(4);
       rho_site_fruitAB ~ exponential(1);
       // Pars linear model
       beta_temp_fruitAB ~ normal(0, 0.5);
       beta_rain_fruitAB ~ normal(0, 0.5);
       // noise parametert of the likelihood
       sigma_fruitAB ~ exponential(1);

       //////////////////////////// Fruit richness //////////////////////
       ////////////////////////////////////////////////////////////////
       // priors for periodic GP
       eta_fruitS ~ normal(0, 0.5);
       gamma_fruitS ~ inv_gamma(5, 5);
       sigma_f_fruitS ~ cauchy(0, 1);
       // GP quadratic
       z_sites_fruitS ~ normal(0, 0.25);
       eta_site_fruitS ~ exponential(4);
       rho_site_fruitS ~ exponential(1);
       // Pars linear model
       beta_temp_fruitS ~ normal(0, 0.5);
       beta_rain_fruitS ~ normal(0, 0.5);
       // noise parametert of the likelihood
       sigma_fruitS ~ exponential(1);

       //////////////////////////// Fru. Nut. HV //////////////////////
       ////////////////////////////////////////////////////////////////
       // priors for periodic GP
       eta_nut_HV ~ normal(0, 0.5);
       gamma_nut_HV ~ inv_gamma(5, 5);
       sigma_f_nut_HV ~ cauchy(0, 1);
       // GP quadratic
       z_sites_nut_HV ~ normal(0, 0.25);
       eta_site_nut_hv ~ exponential(4);
       rho_site_nut_hv ~ exponential(1);
       // Pars linear model
       beta_temp_nut_HV ~ normal(0, 0.5);
       beta_rain_nut_HV ~ normal(0, 0.5);
       beta_fruitAB_nut_HV ~ normal(0, 0.5);
       beta_fruitS_nut_HV ~ normal(0, 0.5);
       beta_BIRD_nut_HV ~ normal(0, 0.5);
       // noise parametert of the likelihood
       sigma_nut_HV ~ exponential(1);
       //
       //////////////////////////// Fru. morpho. HV //////////////////////
       // ////////////////////////////////////////////////////////////////
       // priors for periodic GP
       eta_mor_HV ~ normal(0, 0.5);
       gamma_mor_HV ~ inv_gamma(5, 5);
       sigma_f_mor_HV ~ cauchy(0, 1);
       // GP quadratic
       z_sites_mor_HV ~ normal(0, 0.25);
       eta_site_mor_hv ~ exponential(4);
       rho_site_mor_hv ~ exponential(1);
       // Pars linear model
       beta_temp_mor_HV ~ normal(0, 0.5);
       beta_rain_mor_HV ~ normal(0, 0.5);
       beta_fruitAB_mor_HV ~ normal(0, 0.5);
       beta_fruitS_mor_HV ~ normal(0, 0.5);
       beta_BIRD_mor_HV ~ normal(0, 0.5);
       // noise parametert of the likelihood
       sigma_mor_HV ~ exponential(1);

       // //////////////////////////// Plants HV //////////////////////
       // ////////////////////////////////////////////////////////////////
       // priors for periodic GP
       eta_PLAN_HV ~ normal(0, 0.5);
       gamma_PLAN_HV ~ inv_gamma(5, 5);
       sigma_f_PLAN_HV ~ cauchy(0, 1);
       // GP quadratic
       z_sites_PLAN_HV ~ normal(0, 0.5);
       eta_site_PLAN_hv ~ exponential(4);
       rho_site_PLAN_hv ~ exponential(1);
       // Pars linear model
       beta_temp_PLAN_HV ~ normal(0, 0.5);
       beta_rain_PLAN_HV ~ normal(0, 0.5);
       beta_fruitAB_PLAN_HV ~ normal(0, 0.5);
       beta_fruitS_PLAN_HV ~ normal(0, 0.5);
       beta_nut_PLAN_HV ~ normal(0, 0.5);
       beta_mor_PLAN_HV ~ normal(0, 0.5);
       beta_BIRD_PLAN_HV ~ normal(0, 0.5);
       // noise parametert of the likelihood
       sigma_PLAN_HV ~ exponential(1);

       ////////////////////////// Birds HV //////////////////////
       //////////////////////////////////////////////////////////////
       // priors for periodic GP
       eta_BIRD_HV ~ normal(0, 0.5);
       gamma_BIRD_HV ~ inv_gamma(5, 5);
       sigma_f_BIRD_HV ~ cauchy(0, 1);
       // GP quadratic
       z_sites_BIRD_HV ~ normal(0, 0.5);
       eta_site_BIRD_hv ~ exponential(4);
       rho_site_BIRD_hv ~ exponential(1);
       // Pars linear model
       beta_temp_BIRD_HV ~ normal(0, 0.5);
       beta_rain_BIRD_HV ~ normal(0, 0.5);
       beta_fruitAB_BIRD_HV ~ normal(0, 0.5);
       beta_fruitS_BIRD_HV ~ normal(0, 0.5);
       // noise parametert of the likelihood
       sigma_BIRD_HV ~ exponential(1);

      // //////////////////////////// Nestedness //////////////////////
      // ////////////////////////////////////////////////////////////////
      // priors for periodic GP
      eta_NSS ~ normal(0, 0.5);
      gamma_NSS ~ inv_gamma(5, 5);
      sigma_f_NSS ~ cauchy(0, 1);
      // GP quadratic
      z_sites_NSS ~ normal(0, 0.25);
      eta_site_NSS ~ exponential(4);
      rho_site_NSS ~ exponential(1);
      // Pars linear model
      beta_temp_NSS ~ normal(0, 0.5);
      beta_rain_NSS ~ normal(0, 0.5);
      beta_fruitAB_NSS ~ normal(0, 0.5);
      beta_fruitS_NSS ~ normal(0, 0.5);
      beta_nut_NSS ~ normal(0, 0.5);
      beta_mor_NSS ~ normal(0, 0.5);
      beta_PLAN_NSS ~ normal(0, 0.5);
      beta_BIRD_NSS ~ normal(0, 0.5);
      // beta_NW_NSS ~ normal(0, 0.5);
      // noise parametert of the likelihood
      sigma_NSS ~ exponential(1);


      //////////////////////////// fruit abundance //////////////////////
       ////////////////////////////////////////////////////////////////
       for (i in 1:N) {
         fruit_abun[i] ~ student_t(7, alpha_fruitAB[site[i]] + f_fruitAB[month[i]] +
                                    beta_temp_fruitAB * z_temperature[i] +
                                    beta_rain_fruitAB * z_rainfall[i],
                                    sigma_fruitAB);
       }

      //////////////////////////// fruit richness //////////////////////
       ////////////////////////////////////////////////////////////////
       for (i in 1:N) {
         fruit_S[i] ~ student_t(7, alpha_fruitS[site[i]] + f_fruitS[month[i]] +
                                    beta_temp_fruitS * z_temperature[i] +
                                    beta_rain_fruitS * z_rainfall[i],
                                    sigma_fruitS);
       }

       //////////////////////////// Fru. Nut. HV //////////////////////
       ////////////////////////////////////////////////////////////////
       for (i in 1:N) {
         HV_plants_nut[i] ~ student_t(7, alpha_nut_HV[site[i]] + f_nut_HV[month[i]] +
                                   beta_temp_nut_HV * z_temperature[i] +
                                   beta_rain_nut_HV * z_rainfall[i] +
                                   beta_BIRD_nut_HV * HV_bird[i] +
                                   beta_fruitAB_nut_HV * fruit_abun[i] +
                                   beta_fruitS_nut_HV * fruit_S[i],
                                   sigma_nut_HV);
       }
       //
       // //////////////////////////// Fru. morpho. HV //////////////////////
       // ////////////////////////////////////////////////////////////////
       for (i in 1:N) {
         HV_plants_morfo[i] ~ student_t(15, alpha_mor_HV[site[i]] + f_mor_HV[month[i]] +
                                     beta_temp_mor_HV * z_temperature[i] +
                                     beta_rain_mor_HV * z_rainfall[i] +
                                     beta_fruitAB_mor_HV * fruit_abun[i] +
                                     beta_fruitS_mor_HV * fruit_S[i] +
                                     beta_BIRD_mor_HV * HV_bird[i],
                                     sigma_mor_HV);
       }

       // //////////////////////////// Plants HV //////////////////////
       // ////////////////////////////////////////////////////////////////
       for (i in 1:N) {
         HV_plant[i] ~ student_t(7, alpha_PLAN_HV[site[i]] + f_PLAN_HV[month[i]] +
                                    beta_temp_PLAN_HV * z_temperature[i] +
                                    beta_rain_PLAN_HV * z_rainfall[i] +
                                    beta_fruitAB_PLAN_HV * fruit_abun[i] +
                                    beta_fruitS_PLAN_HV * fruit_S[i] +
                                    beta_nut_PLAN_HV * HV_plants_nut[i] +
                                    beta_mor_PLAN_HV * HV_plants_morfo[i] +
                                    beta_BIRD_PLAN_HV * HV_bird[i],
                                    sigma_PLAN_HV);
       }

       ////////////////////////// Birds HV //////////////////////
       //////////////////////////////////////////////////////////////
       for (i in 1:N) {
         HV_bird[i] ~ student_t(7, alpha_BIRD_HV[site[i]] + f_BIRD_HV[month[i]] +
                                    beta_temp_BIRD_HV * z_temperature[i] +
                                    beta_rain_BIRD_HV * z_rainfall[i] +
                                    beta_fruitAB_BIRD_HV * fruit_abun[i] +
                                    beta_fruitS_BIRD_HV * fruit_S[i],
                                    sigma_BIRD_HV);
       }

      // //////////////////////////// Nestedness //////////////////////
      // ////////////////////////////////////////////////////////////////
      for (i in 1:N) {
        nestedness[i] ~ student_t(7, alpha_NSS[site[i]] + f_NSS[month[i]] +
                                       beta_temp_NSS * z_temperature[i] +
                                       beta_rain_NSS * z_rainfall[i] +
                                       beta_fruitAB_NSS * fruit_abun[i] +
                                       beta_fruitS_NSS * fruit_S[i] +
                                       beta_nut_NSS * HV_plants_nut[i] +
                                       beta_mor_NSS * HV_plants_morfo[i] +
                                       beta_PLAN_NSS * HV_plant[i] +
                                       beta_BIRD_NSS * HV_bird[i],
                                       sigma_NSS);
      }

    }

    generated quantities {

       //////////////////////////// Fruit abundance //////////////////////
       ////////////////////////////////////////////////////////////////
       array[N] real ppcheck_fruitAB;
       vector[N] mu_fruitAB;

       for (i in 1:N) {
         mu_fruitAB[i] = alpha_fruitAB[site[i]] + f_fruitAB[month[i]] +
                         beta_temp_fruitAB * z_temperature[i] +
                         beta_rain_fruitAB * z_rainfall[i];
       }

       ppcheck_fruitAB = student_t_rng(7, mu_fruitAB, sigma_fruitAB);


       //////////////////////////// Fruit richness //////////////////////
       ////////////////////////////////////////////////////////////////
       array[N] real ppcheck_fruitS;
       vector[N] mu_fruitS;

       for (i in 1:N) {
         mu_fruitS[i] = alpha_fruitS[site[i]] + f_fruitS[month[i]] +
                        beta_temp_fruitS * z_temperature[i] +
                        beta_rain_fruitS * z_rainfall[i];
       }

       ppcheck_fruitS = student_t_rng(7, mu_fruitS, sigma_fruitS);


       //////////////////////////// Fru. Nut. HV //////////////////////
       ////////////////////////////////////////////////////////////////
       array[N] real ppcheck_nut_HV;
       vector[N] mu_nut_HV;

       for (i in 1:N) {
         mu_nut_HV[i] = alpha_nut_HV[site[i]] + f_nut_HV[month[i]] +
                        beta_temp_nut_HV * z_temperature[i] +
                        beta_rain_nut_HV * z_rainfall[i] +
                        beta_fruitAB_nut_HV * fruit_abun[i] +
                        beta_fruitS_nut_HV * fruit_S[i] +
                        beta_BIRD_nut_HV * HV_bird[i];
       }

       ppcheck_nut_HV = student_t_rng(7, mu_nut_HV, sigma_nut_HV);
       //
       // //////////////////////////// Fru. morpho. HV //////////////////////
       // ////////////////////////////////////////////////////////////////
       array[N] real ppcheck_mor_HV;
       vector[N] mu_mor_HV;

       for (i in 1:N) {
         mu_mor_HV[i] = alpha_mor_HV[site[i]] + f_mor_HV[month[i]] +
                        beta_temp_mor_HV * z_temperature[i] +
                        beta_rain_mor_HV * z_rainfall[i] +
                        beta_fruitAB_mor_HV * fruit_abun[i] +
                        beta_fruitS_mor_HV * fruit_S[i] +
                        beta_BIRD_mor_HV * HV_bird[i];
       }

       ppcheck_mor_HV = student_t_rng(15, mu_mor_HV, sigma_mor_HV);

       // //////////////////////////// Plants HV //////////////////////
       // ////////////////////////////////////////////////////////////////
       array[N] real ppcheck_PLAN_HV;
       vector[N] mu_PLAN_HV;

       for (i in 1:N) {
         mu_PLAN_HV[i] = alpha_PLAN_HV[site[i]] + f_PLAN_HV[month[i]] +
                         beta_temp_PLAN_HV * z_temperature[i] +
                         beta_rain_PLAN_HV * z_rainfall[i] +
                         beta_fruitAB_PLAN_HV * fruit_abun[i] +
                         beta_fruitS_PLAN_HV * fruit_S[i] +
                         beta_nut_PLAN_HV * HV_plants_nut[i] +
                         beta_mor_PLAN_HV * HV_plants_morfo[i] +
                         beta_BIRD_PLAN_HV * HV_bird[i];
       }

       ppcheck_PLAN_HV = student_t_rng(7, mu_PLAN_HV, sigma_PLAN_HV);

       ////////////////////////// Birds HV //////////////////////
       //////////////////////////////////////////////////////////////
       array[N] real ppcheck_BIRD_HV;
       vector[N] mu_BIRD_HV;

       for (i in 1:N) {
         mu_BIRD_HV[i] = alpha_BIRD_HV[site[i]] + f_BIRD_HV[month[i]] +
                         beta_temp_BIRD_HV * z_temperature[i] +
                         beta_rain_BIRD_HV * z_rainfall[i] +
                         beta_fruitAB_BIRD_HV * fruit_abun[i] +
                         beta_fruitS_BIRD_HV * fruit_S[i];
       }

       ppcheck_BIRD_HV = student_t_rng(7, mu_BIRD_HV, sigma_BIRD_HV);

      ////////////////////////// Nestedness //////////////////////
      //////////////////////////////////////////////////////////////
      array[N] real ppcheck_NSS;
      vector[N] mu_NSS;

      for (i in 1:N) {
        mu_NSS[i] = alpha_NSS[site[i]] + f_NSS[month[i]] +
                    beta_temp_NSS * z_temperature[i] +
                    beta_rain_NSS * z_rainfall[i] +
                    beta_fruitAB_NSS * fruit_abun[i] +
                    beta_fruitS_NSS * fruit_S[i] +
                    beta_nut_NSS * HV_plants_nut[i] +
                    beta_mor_NSS * HV_plants_morfo[i] +
                    beta_PLAN_NSS * HV_plant[i] +
                    beta_BIRD_NSS * HV_bird[i];
      }

      ppcheck_NSS = student_t_rng(7, mu_NSS, sigma_NSS);

    }

    "
    )


file <- paste0(getwd(), '/Nestedness_global.stan')
fit_nestednessG <- cmdstan_model(file, compile = T)

mod_nestednessG <-
  fit_nestednessG$sample(
    data = dat,
    iter_sampling = 2e3,
    iter_warmup = 500,
    thin = 5,
    chains = 3,
    parallel_chains = 3,
    seed = 5
  )

summary_nestednessG <- mod_nestednessG$summary()

par(mfrow = c(3, 3), mar = c(4, 4, 2, 1))
mod_diagnostics(mod_nestednessG, 
                summary_nestednessG[grep('fruitAB', summary_nestednessG$variable), ])
title('Fruit abundance')
mod_diagnostics(mod_nestednessG, 
                summary_nestednessG[grep('fruitS', summary_nestednessG$variable), ])
title('Fruit richness')
mod_diagnostics(mod_nestednessG, 
                summary_nestednessG[grep('nut_HV', summary_nestednessG$variable), ])
title('Frut. Nut. HV')
mod_diagnostics(mod_nestednessG, 
                summary_nestednessG[grep('mor_HV', summary_nestednessG$variable), ])
title('Frut. Morpho. HV')
mod_diagnostics(mod_nestednessG, 
                summary_nestednessG[grep('PLAN_HV', summary_nestednessG$variable), ])
title('Plants HV')
mod_diagnostics(mod_nestednessG, 
                summary_nestednessG[grep('BIRD_HV', summary_nestednessG$variable), ])
title('Birds HV')
mod_diagnostics(mod_nestednessG,
                summary_nestednessG[grep('NSS', summary_nestednessG$variable), ])
title('nestednessG')
par(mfrow = c(1, 1))



par(mfrow = c(3, 3), mar = c(4, 4, 1, 1))
ppcheck_nestednessG <- mod_nestednessG$draws('ppcheck_fruitAB', format = 'matrix')
plot(density(dat$fruit_abun), main = '', xlab = 'Fruit abundance',
     ylim = c(0, 1))
for (i in 1:100) lines(density(ppcheck_nestednessG[i, ]), lwd = 0.1)
lines(density(dat$fruit_abun), col = 'red', lwd = 2)

ppcheck_nestednessG <- mod_nestednessG$draws('ppcheck_fruitS', format = 'matrix')
plot(density(dat$fruit_S), main = '', xlab = 'Fruit richness',
     ylim = c(0, 0.7))
for (i in 1:100) lines(density(ppcheck_nestednessG[i, ]), lwd = 0.1)
lines(density(dat$fruit_S), col = 'red', lwd = 2)

ppcheck_nestednessG <- mod_nestednessG$draws('ppcheck_nut_HV', format = 'matrix')
plot(density(dat$HV_plants_nut), main = '', xlab = 'Frut. Nut. HV',
     ylim = c(0, 0.5))
for (i in 1:100) lines(density(ppcheck_nestednessG[i, ]), lwd = 0.1)
lines(density(dat$HV_plants_nut), col = 'red', lwd = 2)

ppcheck_nestednessG <- mod_nestednessG$draws('ppcheck_mor_HV', format = 'matrix')
plot(density(dat$HV_plants_morfo), main = '', xlab = 'Frut. Morpho. HV',
     ylim = c(0, 0.5))
for (i in 1:100) lines(density(ppcheck_nestednessG[i, ]), lwd = 0.1)
lines(density(dat$HV_plants_morfo), col = 'red', lwd = 2)

ppcheck_nestednessG <- mod_nestednessG$draws('ppcheck_PLAN_HV', format = 'matrix')
plot(density(dat$HV_plant), main = '', xlab = 'Plants HV',
     ylim = c(0, 0.7))
for (i in 1:100) lines(density(ppcheck_nestednessG[i, ]), lwd = 0.1)
lines(density(dat$HV_plant), col = 'red', lwd = 2)

ppcheck_nestednessG <- mod_nestednessG$draws('ppcheck_BIRD_HV', format = 'matrix')
plot(density(dat$HV_bird), main = '', xlab = 'Birds HV',
     ylim = c(0, 0.5))
for (i in 1:100) lines(density(ppcheck_nestednessG[i, ]), lwd = 0.1)
lines(density(dat$HV_bird), col = 'red', lwd = 2)

ppcheck_nestednessG <- mod_nestednessG$draws('ppcheck_NSS', format = 'matrix')
plot(density(dat$nestedness), main = '', xlab = 'Nestedness',
     ylim = c(0, 0.5))
for (i in 1:200) lines(density(ppcheck_nestednessG[i, ]), lwd = 0.1)
lines(density(dat$nestedness), col = 'red', lwd = 2)
par(mfrow = c(1, 1))


pars_nestednessG <- summary_nestednessG[grep('^beta', summary_nestednessG$variable), ]$variable

betas_nestednessG <- mod_nestednessG$draws(pars_nestednessG, format = 'df')

betas_nestednessG <- betas_nestednessG[, -c((ncol(betas_nestednessG)-2):ncol(betas_nestednessG))]

betas_nestednessG <- 
  lapply(1:ncol(betas_nestednessG), 
         function(i) {
           x <- betas_nestednessG[[i]]
           var <- colnames(betas_nestednessG)[[i]]
           
           df <- 
             tibble(var = var, 
                    mu = mean(x), 
                    li = quantile(x, 0.025), 
                    ls = quantile(x, 0.975),
                    p_pos = mean(x > 0), 
                    p_neg = mean(x < 0))
           
           df$effect <- df$p_pos > 0.7 | df$p_neg > 0.7
           
           df
           
         })


betas_nestednessG <- do.call('rbind', betas_nestednessG)

betas_nestednessG_LABS <- 
  betas_nestednessG[grep('NSS', betas_nestednessG$var), ]


betas_nestednessG_z <- betas_nestednessG

betas_nestednessG_z$var <- 
  gsub('^(beta_)([A-Za-z]*)(_)(.*)$', '\\2 -> \\4', betas_nestednessG_z$var)

colnames(betas_nestednessG_z) <- 
  c('Directional effect', 'Mean', '0.025%', '0.975%', 'P(effect > 0)', 
    'P(effect < 0)', 'effect')

knitr::kable(betas_nestednessG_z, 
             digits = 2, 
             caption = 'Average values, 95% CIs, and probabilities of directional effects for the relationships assessed in the generative model. An effect is marked as TRUE if P(effect > 0) > 0.7 or P(effect < 0) > 0.7 or.')


betas_nestednessG <- betas_nestednessG[betas_nestednessG$effect == T, ]

betas_nestednessGEFFECT <- betas_nestednessG[grep('NSS', betas_nestednessG$var), ]

pars_nestednessG <- 
  summary_nestednessG[c(grep('^beta', summary_nestednessG$variable), 
                       grep('^alpha', summary_nestednessG$variable), 
                       grep('^f', summary_nestednessG$variable), 
                       grep('^sigma', summary_nestednessG$variable)), ]$variable

post_nestednessG <- mod_nestednessG$draws(pars_nestednessG, format = 'df')

post_nestednessG <- 
  lapply(c('^beta', '^alpha', '^f', '^sigma'), FUN = 
           function(x) {
             post_nestednessG[, grep(x, colnames(post_nestednessG))]
           })

names(post_nestednessG) <- c('beta', 'alpha', 'f', 'sigma')

estimated_nestednessG_effect <-
  lapply(1, FUN =
           function(x) {
             
             
             f_ <- 'f_NSS['
             s <- grep('(.*)(NSS)(.*)', colnames(post_nestednessG$alpha))
             est <- lapply(1:12, FUN =
                             function(j) {
                               
                               f_ <- paste0(f_, j, ']')
                               
                               effect <-
                                 apply(post_nestednessG$alpha[, s], 1, mean) +
                                 post_nestednessG$f[[f_]] +
                                 post_nestednessG$beta$beta_temp_NSS * mean(dat$z_temperature) +
                                 post_nestednessG$beta$beta_rain_NSS * mean(dat$z_rainfall) +
                                 post_nestednessG$beta$beta_fruitAB_NSS * mean(dat$fruit_abun) +
                                 post_nestednessG$beta$beta_fruitS_NSS * mean(dat$fruit_S) +
                                 post_nestednessG$beta$beta_nut_NSS * mean(dat$HV_plants_nut) +
                                 post_nestednessG$beta$beta_mor_NSS * mean(dat$HV_plants_morfo) +
                                 post_nestednessG$beta$beta_PLAN_NSS * mean(dat$HV_plant) +
                                 post_nestednessG$beta$beta_BIRD_NSS * mean(dat$HV_bird)
                               set.seed(123)
                               effect <- sample(effect, 1e3)

                               tibble(est = effect,
                                      month = j, 
                                      class = 'mNSSu')
                               
                             })
             
             est <- do.call('rbind', est)
             est
             
           })

estimated_nestednessG_effect <- 
  do.call('rbind', estimated_nestednessG_effect)


month_abbr <- na.omit(ymd(paste0("2023-", estimated_nestednessG_effect$month, "-01"), 
                          label = TRUE, abbr = TRUE))

estimated_nestednessG_effect$month2 <- month_abbr



estimated_nestednessG_effect |> 
  group_by(month2) |> 
  transmute(mu = mean(est), 
            li = quantile(est, 0.025), 
            ls = quantile(est, 0.975)) |> 
  unique() |> 
  ggplot() +
  geom_ribbon(aes(month2, ymin = li, ymax = ls), alpha = 0.5) +
  geom_line(aes(month2, mu)) +
  scale_x_date(date_labels = '%b', date_breaks = '2 month')


cat(file = 'Modularity_global.stan',
    "
    functions {
      vector GP_periodic(int period,        // periodicity
                         real gamma,        // smoothing term of the GP
                         real sigma,        // noise paramether
                         vector eta) {      // latent variable for each month

                         int M = period;
                         matrix[M, M] K;
                         matrix[M, M] L_K;

                         for (i in 1:(M - 1)) {
                           for (j in (i+1):M) {
                               real distance = abs(i - j);
                               real periodic_distance = fmin(distance, period - distance);
                               K[i, j] = sigma^2 * exp(-2 * square(sin(pi()*periodic_distance/period))/gamma^2);
                               K[j, i] = K[i, j];   // filling the lower triangle
                            }
                             K[i, i] = sigma^2 + 1e-9;   // small values to guarante stability
                          }
                           K[M, M] = sigma^2 + 1e-9;    // small values to guarante stability
                           return cholesky_decompose(K) * eta;

                        }

      matrix GP_quadratic(matrix x,
                          real eta,
                          real rho,
                          real delta) {

                          int N = dims(x)[1];
                          matrix[N, N] K;
                          matrix[N, N] L_K;

                          for (i in 1:(N-1)) {
                            K[i, i] = eta + delta;
                            for (j in (i+1):N) {
                              K[i, j] = square(eta) * exp(-rho * square(x[i, j]));
                              K[j, i] = K[i, j];
                            }
                          }

                          K[N, N] = eta + delta;
                          L_K = cholesky_decompose(K);
                          return L_K;
                          }
    }

    data{
      int N;
      int N_sites;
      int N_month;
      array[N] int network_size;
      vector[N] modularity;
      vector[N] nestedness;
      vector[N] H2;
      array[N] int month;
      array[N] int site;
      vector[N] HV_network;
      vector[N] HV_plant;
      vector[N] HV_plants_morfo;
      vector[N] HV_plants_nut;
      vector[N] HV_bird;
      vector[N] z_temperature;
      vector[N] z_rainfall;
      vector[N] fruit_S;
      vector[N] fruit_abun;
      matrix[N_sites, N_sites] dist_sites;
    }

    parameters {

      //////////////////////////// Fruit abundance //////////////////////
      ////////////////////////////////////////////////////////////////

      // GP periodic
      real<lower = 0> gamma_fruitS;
      real<lower = 0> sigma_f_fruitS;
      vector[N_month] eta_fruitS;
      // GP quadratic
      vector[N_sites] z_sites_fruitS;
      real<lower = 0> eta_site_fruitS;
      real<lower = 0> rho_site_fruitS;

      // Pars linear model
      real beta_temp_fruitS;
      real beta_rain_fruitS;
      // noise for likelihood function
      real<lower = 0> sigma_fruitS;


      //////////////////////////// Fruit abundance //////////////////////
      ////////////////////////////////////////////////////////////////

      // GP periodic
      real<lower = 0> gamma_fruitAB;
      real<lower = 0> sigma_f_fruitAB;
      vector[N_month] eta_fruitAB;
      // GP quadratic
      vector[N_sites] z_sites_fruitAB;
      real<lower = 0> eta_site_fruitAB;
      real<lower = 0> rho_site_fruitAB;

      // Pars linear model

      real beta_temp_fruitAB;
      real beta_rain_fruitAB;
      // noise for likelihood function
      real<lower = 0> sigma_fruitAB;


      //////////////////////////// Fru. Nut. HV //////////////////////
      ////////////////////////////////////////////////////////////////

      // GP periodic
      real<lower = 0> gamma_nut_HV;
      real<lower = 0> sigma_f_nut_HV;
      vector[N_month] eta_nut_HV;
      // GP quadratic
      vector[N_sites] z_sites_nut_HV;
      real<lower = 0> eta_site_nut_hv;
      real<lower = 0> rho_site_nut_hv;

      // Pars linear model

      real beta_temp_nut_HV;
      real beta_rain_nut_HV;
      real beta_fruitAB_nut_HV;
      real beta_fruitS_nut_HV;
      real beta_BIRD_nut_HV;
      // noise for likelihood function
      real<lower = 0> sigma_nut_HV;

      //////////////////////////// Fru. morpho. HV //////////////////////
      ////////////////////////////////////////////////////////////////
      // GP periodic
      real<lower = 0> gamma_mor_HV;
      real<lower = 0> sigma_f_mor_HV;
      vector[N_month] eta_mor_HV;
      // GP quadratic
      vector[N_sites] z_sites_mor_HV;
      real<lower = 0> eta_site_mor_hv;
      real<lower = 0> rho_site_mor_hv;

      // Pars linear model

      real beta_temp_mor_HV;
      real beta_rain_mor_HV;
      real beta_fruitAB_mor_HV;
      real beta_fruitS_mor_HV;
      real beta_BIRD_mor_HV;
      // noise for likelihood function
      real<lower = 0> sigma_mor_HV;
      real nu;

      // //////////////////////////// Plants HV //////////////////////
      // ////////////////////////////////////////////////////////////////
      // GP periodic
      real<lower = 0> gamma_PLAN_HV;
      real<lower = 0> sigma_f_PLAN_HV;
      vector[N_month] eta_PLAN_HV;
      // GP quadratic
      vector[N_sites] z_sites_PLAN_HV;
      real<lower = 0> eta_site_PLAN_hv;
      real<lower = 0> rho_site_PLAN_hv;

      // Pars linear model
      real beta_temp_PLAN_HV;
      real beta_rain_PLAN_HV;
      real beta_fruitAB_PLAN_HV;
      real beta_fruitS_PLAN_HV;
      real beta_nut_PLAN_HV;
      real beta_mor_PLAN_HV;
      real beta_BIRD_PLAN_HV;
      // noise for likelihood function
      real<lower = 0> sigma_PLAN_HV;

      //////////////////////////// Birsds HV //////////////////////
      // ////////////////////////////////////////////////////////////////
      // GP periodic
      real<lower = 0> gamma_BIRD_HV;
      real<lower = 0> sigma_f_BIRD_HV;
      vector[N_month] eta_BIRD_HV;
      // GP quadratic
      vector[N_sites] z_sites_BIRD_HV;
      real<lower = 0> eta_site_BIRD_hv;
      real<lower = 0> rho_site_BIRD_hv;

      // Pars linear model
      real beta_temp_BIRD_HV;
      real beta_rain_BIRD_HV;
      real beta_fruitAB_BIRD_HV;
      real beta_fruitS_BIRD_HV;
      // noise for likelihood function
      real<lower = 0> sigma_BIRD_HV;

      //////////////////////////// Network metrics //////////////////////
      ////////////////////////////////////////////////////////////////
      // //////////////////////////// Modularity //////////////////////
      // ////////////////////////////////////////////////////////////////
      // GP periodic
      real<lower = 0> gamma_MODU;
      real<lower = 0> sigma_f_MODU;
      vector[N_month] eta_MODU;
      // GP quadratic
      vector[N_sites] z_sites_MODU;
      real<lower = 0> eta_site_MODU;
      real<lower = 0> rho_site_MODU;

      // Pars linear model
      real beta_temp_MODU;
      real beta_rain_MODU;
      real beta_fruitAB_MODU;
      real beta_fruitS_MODU;
      real beta_nut_MODU;
      real beta_mor_MODU;
      real beta_PLAN_MODU;
      real beta_BIRD_MODU;
      // noise for likelihood function
      real<lower = 0> sigma_MODU;

    }

    transformed parameters {

        //////////////////////////// Fruit abundance //////////////////////
    ////////////////////////////////////////////////////////////////
      // GP periodic
      vector[N_month] f_fruitAB;

      f_fruitAB = GP_periodic(
              12,
              gamma_fruitAB,
              sigma_f_fruitAB,
              eta_fruitAB
        );

      // GP quadratic
      vector[N_sites] alpha_fruitAB;
      matrix[N_sites, N_sites] L_K_fruitAB;
      L_K_fruitAB = GP_quadratic(dist_sites,
                                eta_site_fruitAB,
                                rho_site_fruitAB, 0.001);
      alpha_fruitAB = L_K_fruitAB * z_sites_fruitAB;


    //////////////////////////// Fruit richness //////////////////////
    ////////////////////////////////////////////////////////////////
      // GP periodic
      vector[N_month] f_fruitS;

      f_fruitS = GP_periodic(
              12,
              gamma_fruitS,
              sigma_f_fruitS,
              eta_fruitS
        );

      // GP quadratic
      vector[N_sites] alpha_fruitS;
      matrix[N_sites, N_sites] L_K_fruitS;
      L_K_fruitS = GP_quadratic(dist_sites,
                                eta_site_fruitS,
                                rho_site_fruitS, 0.001);
      alpha_fruitS = L_K_fruitS * z_sites_fruitS;

    //////////////////////////// Fru. Nut. HV //////////////////////
    ////////////////////////////////////////////////////////////////
      // GP periodic
      vector[N_month] f_nut_HV;

      f_nut_HV = GP_periodic(
              12,
              gamma_nut_HV,
              sigma_f_nut_HV,
              eta_nut_HV
        );

      // GP quadratic
      vector[N_sites] alpha_nut_HV;
      matrix[N_sites, N_sites] L_K_nut_HV;
      L_K_nut_HV = GP_quadratic(dist_sites,
                                eta_site_nut_hv,
                                rho_site_nut_hv, 0.001);
      alpha_nut_HV = L_K_nut_HV * z_sites_nut_HV;

    //   //////////////////////////// Fru. morpho. HV //////////////////////
    //   ////////////////////////////////////////////////////////////////
      // GP periodic
      vector[N_month] f_mor_HV;

      f_mor_HV = GP_periodic(
              12,
              gamma_mor_HV,
              sigma_f_mor_HV,
              eta_mor_HV
        );

      // GP quadratic
      vector[N_sites] alpha_mor_HV;
      matrix[N_sites, N_sites] L_K_mor_HV;
      L_K_mor_HV = GP_quadratic(dist_sites,
                                eta_site_mor_hv,
                                rho_site_mor_hv, 0.001);
      alpha_mor_HV = L_K_mor_HV * z_sites_mor_HV;
    //
    //   //////////////////////////// Plants HV //////////////////////
    //   ////////////////////////////////////////////////////////////////
      // GP periodic
      vector[N_month] f_PLAN_HV;

      f_PLAN_HV = GP_periodic(
              12,
              gamma_PLAN_HV,
              sigma_f_PLAN_HV,
              eta_PLAN_HV
        );

      // GP quadratic
      vector[N_sites] alpha_PLAN_HV;
      matrix[N_sites, N_sites] L_K_PLAN_HV;
      L_K_PLAN_HV = GP_quadratic(dist_sites,
                                eta_site_PLAN_hv,
                                rho_site_PLAN_hv, 0.001);
      alpha_PLAN_HV = L_K_PLAN_HV * z_sites_PLAN_HV;

      ////////////////////////// Birsds HV //////////////////////
      //////////////////////////////////////////////////////////////
      // GP periodic
      vector[N_month] f_BIRD_HV;

      f_BIRD_HV = GP_periodic(
              12,
              gamma_BIRD_HV,
              sigma_f_BIRD_HV,
              eta_BIRD_HV
        );

      // GP quadratic
      vector[N_sites] alpha_BIRD_HV;
      matrix[N_sites, N_sites] L_K_BIRD_HV;
      L_K_BIRD_HV = GP_quadratic(dist_sites,
                                eta_site_BIRD_hv,
                                rho_site_BIRD_hv, 0.001);
      alpha_BIRD_HV = L_K_BIRD_HV * z_sites_BIRD_HV;

      // //////////////////////////// Modularity //////////////////////
      // ////////////////////////////////////////////////////////////////
      // GP periodic
      vector[N_month] f_MODU;

      f_MODU = GP_periodic(
              12,
              gamma_MODU,
              sigma_f_MODU,
              eta_MODU);

      // GP quadratic
      vector[N_sites] alpha_MODU;
      matrix[N_sites, N_sites] L_K_MODU;
      L_K_MODU = GP_quadratic(dist_sites,
                              eta_site_MODU,
                              rho_site_MODU, 0.001);
      alpha_MODU = L_K_MODU * z_sites_MODU;

    }

    model {

              //////////////////////////// Fruit abundance //////////////////////
       ////////////////////////////////////////////////////////////////
       // priors for periodic GP
       eta_fruitAB ~ normal(0, 0.5);
       gamma_fruitAB ~ inv_gamma(5, 5);
       sigma_f_fruitAB ~ cauchy(0, 1);
       // GP quadratic
       z_sites_fruitAB ~ normal(0, 0.25);
       eta_site_fruitAB ~ exponential(4);
       rho_site_fruitAB ~ exponential(1);
       // Pars linear model
       beta_temp_fruitAB ~ normal(0, 0.5);
       beta_rain_fruitAB ~ normal(0, 0.5);
       // noise parametert of the likelihood
       sigma_fruitAB ~ exponential(1);

       //////////////////////////// Fruit richness //////////////////////
       ////////////////////////////////////////////////////////////////
       // priors for periodic GP
       eta_fruitS ~ normal(0, 0.5);
       gamma_fruitS ~ inv_gamma(5, 5);
       sigma_f_fruitS ~ cauchy(0, 1);
       // GP quadratic
       z_sites_fruitS ~ normal(0, 0.25);
       eta_site_fruitS ~ exponential(4);
       rho_site_fruitS ~ exponential(1);
       // Pars linear model
       beta_temp_fruitS ~ normal(0, 0.5);
       beta_rain_fruitS ~ normal(0, 0.5);
       // noise parametert of the likelihood
       sigma_fruitS ~ exponential(1);

       //////////////////////////// Fru. Nut. HV //////////////////////
       ////////////////////////////////////////////////////////////////
       // priors for periodic GP
       eta_nut_HV ~ normal(0, 0.5);
       gamma_nut_HV ~ inv_gamma(5, 5);
       sigma_f_nut_HV ~ cauchy(0, 1);
       // GP quadratic
       z_sites_nut_HV ~ normal(0, 0.25);
       eta_site_nut_hv ~ exponential(4);
       rho_site_nut_hv ~ exponential(1);
       // Pars linear model
       beta_temp_nut_HV ~ normal(0, 0.5);
       beta_rain_nut_HV ~ normal(0, 0.5);
       beta_fruitAB_nut_HV ~ normal(0, 0.5);
       beta_fruitS_nut_HV ~ normal(0, 0.5);
       beta_BIRD_nut_HV ~ normal(0, 0.5);
       // noise parametert of the likelihood
       sigma_nut_HV ~ exponential(1);
       //
       //////////////////////////// Fru. morpho. HV //////////////////////
       // ////////////////////////////////////////////////////////////////
       // priors for periodic GP
       eta_mor_HV ~ normal(0, 0.5);
       gamma_mor_HV ~ inv_gamma(5, 5);
       sigma_f_mor_HV ~ cauchy(0, 1);
       // GP quadratic
       z_sites_mor_HV ~ normal(0, 0.25);
       eta_site_mor_hv ~ exponential(4);
       rho_site_mor_hv ~ exponential(1);
       // Pars linear model
       beta_temp_mor_HV ~ normal(0, 0.5);
       beta_rain_mor_HV ~ normal(0, 0.5);
       beta_fruitAB_mor_HV ~ normal(0, 0.5);
       beta_fruitS_mor_HV ~ normal(0, 0.5);
       beta_BIRD_mor_HV ~ normal(0, 0.5);
       // noise parametert of the likelihood
       sigma_mor_HV ~ exponential(1);

       // //////////////////////////// Plants HV //////////////////////
       // ////////////////////////////////////////////////////////////////
       // priors for periodic GP
       eta_PLAN_HV ~ normal(0, 0.5);
       gamma_PLAN_HV ~ inv_gamma(5, 5);
       sigma_f_PLAN_HV ~ cauchy(0, 1);
       // GP quadratic
       z_sites_PLAN_HV ~ normal(0, 0.5);
       eta_site_PLAN_hv ~ exponential(4);
       rho_site_PLAN_hv ~ exponential(1);
       // Pars linear model
       beta_temp_PLAN_HV ~ normal(0, 0.5);
       beta_rain_PLAN_HV ~ normal(0, 0.5);
       beta_fruitAB_PLAN_HV ~ normal(0, 0.5);
       beta_fruitS_PLAN_HV ~ normal(0, 0.5);
       beta_nut_PLAN_HV ~ normal(0, 0.5);
       beta_mor_PLAN_HV ~ normal(0, 0.5);
       beta_BIRD_PLAN_HV ~ normal(0, 0.5);
       // noise parametert of the likelihood
       sigma_PLAN_HV ~ exponential(1);

       ////////////////////////// Birds HV //////////////////////
       //////////////////////////////////////////////////////////////
       // priors for periodic GP
       eta_BIRD_HV ~ normal(0, 0.5);
       gamma_BIRD_HV ~ inv_gamma(5, 5);
       sigma_f_BIRD_HV ~ cauchy(0, 1);
       // GP quadratic
       z_sites_BIRD_HV ~ normal(0, 0.5);
       eta_site_BIRD_hv ~ exponential(4);
       rho_site_BIRD_hv ~ exponential(1);
       // Pars linear model
       beta_temp_BIRD_HV ~ normal(0, 0.5);
       beta_rain_BIRD_HV ~ normal(0, 0.5);
       beta_fruitAB_BIRD_HV ~ normal(0, 0.5);
       beta_fruitS_BIRD_HV ~ normal(0, 0.5);
       // noise parametert of the likelihood
       sigma_BIRD_HV ~ exponential(1);

      // //////////////////////////// Modularity //////////////////////
      // ////////////////////////////////////////////////////////////////
      vector[N] mu_MODU;
      vector[N] p1_MODU;
      vector[N] p2_MODU;
      // priors for periodic GP
      eta_MODU ~ normal(0, 1);
      gamma_MODU ~ inv_gamma(5, 5);
      sigma_f_MODU ~ cauchy(0, 1);
      // GP quadratic
      z_sites_MODU ~ normal(0, 0.25);
      eta_site_MODU ~ exponential(4);
      rho_site_MODU ~ exponential(1);
      // Pars linear model
      beta_fruitAB_MODU ~ normal(0, 0.055);
      beta_fruitS_MODU ~ normal(0, 0.055);
      beta_nut_MODU ~ normal(0, 0.055);
      beta_mor_MODU ~ normal(0, 0.055);
      beta_PLAN_MODU ~ normal(0, 0.055);
      beta_BIRD_MODU ~ normal(0, 0.055);
      // beta_NW_MODU ~ normal(0, 0.055);
      // noise parametert of the likelihood
      sigma_MODU ~ exponential(2);

      //////////////////////////// fruit abundance //////////////////////
       ////////////////////////////////////////////////////////////////
       for (i in 1:N) {
         fruit_abun[i] ~ student_t(7, alpha_fruitAB[site[i]] + f_fruitAB[month[i]] +
                                    beta_temp_fruitAB * z_temperature[i] +
                                    beta_rain_fruitAB * z_rainfall[i],
                                    sigma_fruitAB);
       }

      //////////////////////////// fruit richness //////////////////////
       ////////////////////////////////////////////////////////////////
       for (i in 1:N) {
         fruit_S[i] ~ student_t(7, alpha_fruitS[site[i]] + f_fruitS[month[i]] +
                                    beta_temp_fruitS * z_temperature[i] +
                                    beta_rain_fruitS * z_rainfall[i],
                                    sigma_fruitS);
       }

       //////////////////////////// Fru. Nut. HV //////////////////////
       ////////////////////////////////////////////////////////////////
       for (i in 1:N) {
         HV_plants_nut[i] ~ student_t(7, alpha_nut_HV[site[i]] + f_nut_HV[month[i]] +
                                   beta_temp_nut_HV * z_temperature[i] +
                                   beta_rain_nut_HV * z_rainfall[i] +
                                   beta_BIRD_nut_HV * HV_bird[i] +
                                   beta_fruitAB_nut_HV * fruit_abun[i] +
                                   beta_fruitS_nut_HV * fruit_S[i],
                                   sigma_nut_HV);
       }
       //
       // //////////////////////////// Fru. morpho. HV //////////////////////
       // ////////////////////////////////////////////////////////////////
       for (i in 1:N) {
         HV_plants_morfo[i] ~ student_t(15, alpha_mor_HV[site[i]] + f_mor_HV[month[i]] +
                                     beta_temp_mor_HV * z_temperature[i] +
                                     beta_rain_mor_HV * z_rainfall[i] +
                                     beta_fruitAB_mor_HV * fruit_abun[i] +
                                     beta_fruitS_mor_HV * fruit_S[i] +
                                     beta_BIRD_mor_HV * HV_bird[i],
                                     sigma_mor_HV);
       }

       // //////////////////////////// Plants HV //////////////////////
       // ////////////////////////////////////////////////////////////////
       for (i in 1:N) {
         HV_plant[i] ~ student_t(7, alpha_PLAN_HV[site[i]] + f_PLAN_HV[month[i]] +
                                    beta_temp_PLAN_HV * z_temperature[i] +
                                    beta_rain_PLAN_HV * z_rainfall[i] +
                                    beta_fruitAB_PLAN_HV * fruit_abun[i] +
                                    beta_fruitS_PLAN_HV * fruit_S[i] +
                                    beta_nut_PLAN_HV * HV_plants_nut[i] +
                                    beta_mor_PLAN_HV * HV_plants_morfo[i] +
                                    beta_BIRD_PLAN_HV * HV_bird[i],
                                    sigma_PLAN_HV);
       }

       ////////////////////////// Birds HV //////////////////////
       //////////////////////////////////////////////////////////////
       for (i in 1:N) {
         HV_bird[i] ~ student_t(7, alpha_BIRD_HV[site[i]] + f_BIRD_HV[month[i]] +
                                    beta_temp_BIRD_HV * z_temperature[i] +
                                    beta_rain_BIRD_HV * z_rainfall[i] +
                                    beta_fruitAB_BIRD_HV * fruit_abun[i] +
                                    beta_fruitS_BIRD_HV * fruit_S[i],
                                    sigma_BIRD_HV);
       }
      // //////////////////////////// Modularity //////////////////////
      // ////////////////////////////////////////////////////////////////
      for (i in 1:N) {
        mu_MODU[i] = inv_logit(alpha_MODU[site[i]] + f_MODU[month[i]] +
                     beta_temp_MODU * z_temperature[i] +
                     beta_rain_MODU * z_rainfall[i] +
                     beta_fruitAB_MODU * fruit_abun[i] +
                     beta_fruitS_MODU * fruit_S[i] +
                     beta_nut_MODU * HV_plants_nut[i] +
                     beta_mor_MODU * HV_plants_morfo[i] +
                     beta_PLAN_MODU * HV_plant[i] +
                     beta_BIRD_MODU * HV_bird[i]);
      }

      p1_MODU = mu_MODU * sigma_MODU;
      p2_MODU = (1 - mu_MODU) * sigma_MODU;

      modularity ~ beta(p1_MODU, p2_MODU);

    }

    generated quantities {

       //////////////////////////// Fruit abundance //////////////////////
       ////////////////////////////////////////////////////////////////
       array[N] real ppcheck_fruitAB;
       vector[N] mu_fruitAB;

       for (i in 1:N) {
         mu_fruitAB[i] = alpha_fruitAB[site[i]] + f_fruitAB[month[i]] +
                         beta_temp_fruitAB * z_temperature[i] +
                         beta_rain_fruitAB * z_rainfall[i];
       }

       ppcheck_fruitAB = student_t_rng(7, mu_fruitAB, sigma_fruitAB);


       //////////////////////////// Fruit richness //////////////////////
       ////////////////////////////////////////////////////////////////
       array[N] real ppcheck_fruitS;
       vector[N] mu_fruitS;

       for (i in 1:N) {
         mu_fruitS[i] = alpha_fruitS[site[i]] + f_fruitS[month[i]] +
                        beta_temp_fruitS * z_temperature[i] +
                        beta_rain_fruitS * z_rainfall[i];
       }

       ppcheck_fruitS = student_t_rng(7, mu_fruitS, sigma_fruitS);


       //////////////////////////// Fru. Nut. HV //////////////////////
       ////////////////////////////////////////////////////////////////
       array[N] real ppcheck_nut_HV;
       vector[N] mu_nut_HV;

       for (i in 1:N) {
         mu_nut_HV[i] = alpha_nut_HV[site[i]] + f_nut_HV[month[i]] +
                        beta_temp_nut_HV * z_temperature[i] +
                        beta_rain_nut_HV * z_rainfall[i] +
                        beta_fruitAB_nut_HV * fruit_abun[i] +
                        beta_fruitS_nut_HV * fruit_S[i] +
                        beta_BIRD_nut_HV * HV_bird[i];
       }

       ppcheck_nut_HV = student_t_rng(7, mu_nut_HV, sigma_nut_HV);
       //
       // //////////////////////////// Fru. morpho. HV //////////////////////
       // ////////////////////////////////////////////////////////////////
       array[N] real ppcheck_mor_HV;
       vector[N] mu_mor_HV;

       for (i in 1:N) {
         mu_mor_HV[i] = alpha_mor_HV[site[i]] + f_mor_HV[month[i]] +
                        beta_temp_mor_HV * z_temperature[i] +
                        beta_rain_mor_HV * z_rainfall[i] +
                        beta_fruitAB_mor_HV * fruit_abun[i] +
                        beta_fruitS_mor_HV * fruit_S[i] +
                        beta_BIRD_mor_HV * HV_bird[i];
       }

       ppcheck_mor_HV = student_t_rng(15, mu_mor_HV, sigma_mor_HV);

       // //////////////////////////// Plants HV //////////////////////
       // ////////////////////////////////////////////////////////////////
       array[N] real ppcheck_PLAN_HV;
       vector[N] mu_PLAN_HV;

       for (i in 1:N) {
         mu_PLAN_HV[i] = alpha_PLAN_HV[site[i]] + f_PLAN_HV[month[i]] +
                         beta_temp_PLAN_HV * z_temperature[i] +
                         beta_rain_PLAN_HV * z_rainfall[i] +
                         beta_fruitAB_PLAN_HV * fruit_abun[i] +
                         beta_fruitS_PLAN_HV * fruit_S[i] +
                         beta_nut_PLAN_HV * HV_plants_nut[i] +
                         beta_mor_PLAN_HV * HV_plants_morfo[i] +
                         beta_BIRD_PLAN_HV * HV_bird[i];
       }

       ppcheck_PLAN_HV = student_t_rng(7, mu_PLAN_HV, sigma_PLAN_HV);

       ////////////////////////// Birds HV //////////////////////
       //////////////////////////////////////////////////////////////
       array[N] real ppcheck_BIRD_HV;
       vector[N] mu_BIRD_HV;

       for (i in 1:N) {
         mu_BIRD_HV[i] = alpha_BIRD_HV[site[i]] + f_BIRD_HV[month[i]] +
                         beta_temp_BIRD_HV * z_temperature[i] +
                         beta_rain_BIRD_HV * z_rainfall[i] +
                         beta_fruitAB_BIRD_HV * fruit_abun[i] +
                         beta_fruitS_BIRD_HV * fruit_S[i];
       }

       ppcheck_BIRD_HV = student_t_rng(7, mu_BIRD_HV, sigma_BIRD_HV);


      // ////////////////////////// Modularity //////////////////////
      // //////////////////////////////////////////////////////////////
      array[N] real ppcheck_MODU;
      vector[N] mu_MODU;
      vector[N] p1_MODU;
      vector[N] p2_MODU;
      real variance_fixed;
      real variance_GP_temp;
      real variance_GP_space;
      real variance_total;
      real prop_GP_temp;
      real prop_GP_space;
      real prop_fixed;

      for (i in 1:N) {
        mu_MODU[i] = alpha_MODU[site[i]] + f_MODU[month[i]] +
                    beta_temp_MODU * z_temperature[i] +
                    beta_rain_MODU * z_rainfall[i] +
                    beta_fruitAB_MODU * fruit_abun[i] +
                    beta_fruitS_MODU * fruit_S[i] +
                    beta_nut_MODU * HV_plants_nut[i] +
                    beta_mor_MODU * HV_plants_morfo[i] +
                    beta_PLAN_MODU * HV_plant[i] +
                    beta_BIRD_MODU * HV_bird[i];
      }

      // variance decomposition
      variance_fixed = variance(mu_MODU);
      variance_GP_space = eta_site_MODU^2;
      variance_GP_temp = sigma_f_MODU^2;
      variance_total = variance_fixed + variance_GP_temp +
                       variance_GP_space + ((pi()^2)/3);
      prop_GP_temp = variance_GP_temp / variance_total;
      prop_GP_space = variance_GP_space / variance_total;
      prop_fixed = variance_fixed / variance_total;

      mu_MODU = inv_logit(mu_MODU);
      p1_MODU = mu_MODU * sigma_MODU;
      p2_MODU = (1 - mu_MODU) * sigma_MODU;

      ppcheck_MODU = beta_rng(p1_MODU, p2_MODU);

    }
    "
    )


file <- paste0(getwd(), '/Modularity_global.stan')
fit_MODU_G <- cmdstan_model(file, compile = T)

mod_MODU_G <-
  fit_MODU_G$sample(
    data = dat,
    iter_sampling = 2e3,
    iter_warmup = 500,
    thin = 3,
    chains = 3,
    parallel_chains = 3,
    seed = 5
  )

summary_MODU_G <- mod_MODU_G$summary()

na.omit(summary_MODU_G[summary_MODU_G$ess_bulk < 500,]) |> print(n = 1e3)

par(mfrow = c(3, 3), mar = c(4, 4, 2, 1))
mod_diagnostics(mod_MODU_G, 
                summary_MODU_G[grep('fruitAB', summary_MODU_G$variable), ])
title('Fruit abundance')
mod_diagnostics(mod_MODU_G, 
                summary_MODU_G[grep('fruitS', summary_MODU_G$variable), ])
title('Fruit richness')
mod_diagnostics(mod_MODU_G, 
                summary_MODU_G[grep('nut_HV', summary_MODU_G$variable), ])
title('Frut. Nut. HV')
mod_diagnostics(mod_MODU_G, 
                summary_MODU_G[grep('mor_HV', summary_MODU_G$variable), ])
title('Frut. Morpho. HV')
mod_diagnostics(mod_MODU_G, 
                summary_MODU_G[grep('PLAN_HV', summary_MODU_G$variable), ])
title('Plants HV')
mod_diagnostics(mod_MODU_G, 
                summary_MODU_G[grep('BIRD_HV', summary_MODU_G$variable), ])
title('Birds HV')

mod_diagnostics(mod_MODU_G,
                summary_MODU_G[grep('MODU', summary_MODU_G$variable), ])
title('Modularity')
par(mfrow = c(1, 1))


par(mfrow = c(3, 3), mar = c(4, 4, 1, 1))
ppcheck_MODU_G <- mod_MODU_G$draws('ppcheck_fruitAB', format = 'matrix')
plot(density(dat$fruit_abun), main = '', xlab = 'Fruit abundance',
     ylim = c(0, 1))
for (i in 1:100) lines(density(ppcheck_MODU_G[i, ]), lwd = 0.1)
lines(density(dat$fruit_abun), col = 'red', lwd = 2)

ppcheck_MODU_G <- mod_MODU_G$draws('ppcheck_fruitS', format = 'matrix')
plot(density(dat$fruit_S), main = '', xlab = 'Fruit richness',
     ylim = c(0, 1))
for (i in 1:100) lines(density(ppcheck_MODU_G[i, ]), lwd = 0.1)
lines(density(dat$fruit_S), col = 'red', lwd = 2)

ppcheck_MODU_G <- mod_MODU_G$draws('ppcheck_nut_HV', format = 'matrix')
plot(density(dat$HV_plants_nut), main = '', xlab = 'Frut. Nut. HV',
     ylim = c(0, 1))
for (i in 1:100) lines(density(ppcheck_MODU_G[i, ]), lwd = 0.1)
lines(density(dat$HV_plants_nut), col = 'red', lwd = 2)

ppcheck_MODU_G <- mod_MODU_G$draws('ppcheck_mor_HV', format = 'matrix')
plot(density(dat$HV_plants_morfo), main = '', xlab = 'Frut. Morpho. HV',
     ylim = c(0, 1))
for (i in 1:100) lines(density(ppcheck_MODU_G[i, ]), lwd = 0.1)
lines(density(dat$HV_plants_morfo), col = 'red', lwd = 2)

ppcheck_MODU_G <- mod_MODU_G$draws('ppcheck_PLAN_HV', format = 'matrix')
plot(density(dat$HV_plant), main = '', xlab = 'Plants HV',
     ylim = c(0, 1))
for (i in 1:100) lines(density(ppcheck_MODU_G[i, ]), lwd = 0.1)
lines(density(dat$HV_plant), col = 'red', lwd = 2)

ppcheck_MODU_G <- mod_MODU_G$draws('ppcheck_BIRD_HV', format = 'matrix')
plot(density(dat$HV_bird), main = '', xlab = 'Birds HV',
     ylim = c(0, 1))
for (i in 1:100) lines(density(ppcheck_MODU_G[i, ]), lwd = 0.1)
lines(density(dat$HV_bird), col = 'red', lwd = 2)

ppcheck_MODU_G <- mod_MODU_G$draws('ppcheck_MODU', format = 'matrix')
plot(density(dat$modularity), main = '', xlab = 'Modularity',
     ylim = c(0, 5))
for (i in 1:100) lines(density(ppcheck_MODU_G[i, ]), lwd = 0.1)
lines(density(dat$modularity), col = 'red', lwd = 2)
par(mfrow = c(1, 1))


pars_MODU_G <- summary_MODU_G[grep('^beta', summary_MODU_G$variable), ]$variable

betas_MODU_G <- mod_MODU_G$draws(pars_MODU_G, format = 'df')

betas_MODU_G <- betas_MODU_G[, -c((ncol(betas_MODU_G)-2):ncol(betas_MODU_G))]

betas_MODU_G <- 
  lapply(1:ncol(betas_MODU_G), 
         function(i) {
           x <- betas_MODU_G[[i]]
           var <- colnames(betas_MODU_G)[[i]]
           
           df <- 
             tibble(var = var, 
                    mu = mean(x), 
                    li = quantile(x, 0.025), 
                    ls = quantile(x, 0.975),
                    p_pos = mean(x > 0), 
                    p_neg = mean(x < 0))
           
           df$effect <- df$p_pos > 0.7 | df$p_neg > 0.7
           
           df
           
         })


betas_MODU_G <- do.call('rbind', betas_MODU_G)

betas_MODU_G_LABS <- betas_MODU_G[c(grep('PLAN_MODU', betas_MODU_G$var), 
                                grep('mor_MODU', betas_MODU_G$var)), ]

betas_MODU_GEFFECT <- betas_MODU_G[grep('MODU', betas_MODU_G$var), ]


betas_MODU_G_z <- betas_MODU_G

betas_MODU_G_z$var <- 
  gsub('^(beta_)([A-Za-z]*)(_)(.*)$', '\\2 -> \\4', betas_MODU_G$var)

colnames(betas_MODU_G_z) <- 
  c('Directional effect', 'Mean', '0.025%', '0.975%', 'P(effect > 0)', 
    'P(effect < 0)', 'effect')

knitr::kable(betas_MODU_G_z, 
             digits = 2, 
             caption = 'Average values, 95% CIs, and probabilities of directional effects for the relationships assessed in the generative model. An effect is marked as TRUE if P(effect > 0) > 0.7 or P(effect < 0) > 0.7 or.')


betas_MODU_G <- betas_MODU_G[betas_MODU_G$effect == T, ]

pars_MODU_G <- 
  summary_MODU_G[c(grep('^beta', summary_MODU_G$variable), 
                 grep('^alpha', summary_MODU_G$variable), 
                 grep('^f', summary_MODU_G$variable), 
                 grep('^sigma', summary_MODU_G$variable)), ]$variable

post_MODU_G <- mod_MODU_G$draws(pars_MODU_G, format = 'df')

post_MODU_G <- 
  lapply(c('^beta', '^alpha', '^f', '^sigma'), FUN = 
           function(x) {
             post_MODU_G[, grep(x, colnames(post_MODU_G))]
           })

names(post_MODU_G) <- c('beta', 'alpha', 'f', 'sigma')


post_MODU_G$f$`f_fruitAB[1]`
post_MODU_G$alpha$`alpha_fruitAB[1]`
post_MODU_G$beta$beta_temp_fruitS
betas_MODU_GEFFECT
estimated_MODU_G_effect <- 
  lapply(1, FUN = 
           function(x) {
             
             
             f_ <- 'f_MODU['
             
             mu_MODU <- mean(dat$modularity)
             sd_MODU <- sd(dat$modularity)
             
             s <- grep('(.*)(MODU)(.*)$', colnames(post_MODU_G$alpha))
             
             est <- lapply(1:12, FUN =
                             function(j) {
                               
                               f_ <- paste0(f_, j, ']')
                               
                               m <-
                                 apply(post_MODU_G$alpha[, s], 1, mean) +
                                 post_MODU_G$f[[f_]] +
                                 post_MODU_G$beta$beta_temp_MODU * mean(dat$z_temperature) +
                                 post_MODU_G$beta$beta_rain_MODU * mean(dat$z_rainfall) +
                                 post_MODU_G$beta$beta_fruitAB_MODU * mean(dat$fruit_abun) +
                                 post_MODU_G$beta$beta_fruitS_MODU * mean(dat$fruit_S) +
                                 post_MODU_G$beta$beta_nut_MODU * mean(dat$HV_plants_nut) +
                                 post_MODU_G$beta$beta_mor_MODU * mean(dat$HV_plants_morfo) +
                                 post_MODU_G$beta$beta_PLAN_MODU * mean(dat$HV_plant) +
                                 post_MODU_G$beta$beta_BIRD_MODU * mean(dat$HV_bird)
                               
                               m <- inv_logit(m)
                               set.seed(123)
                               m <- sample(m, 1e3)
                               
                               tibble(est = (m - mu_MODU) / sd_MODU,
                                      month = j, 
                                      class = 'MODU')
                               
                             })
             est <- do.call('rbind', est)
             
             est
             
           })

estimated_MODU_G_effect <- 
  do.call('rbind', estimated_MODU_G_effect)

month_abbr <- na.omit(ymd(paste0("2023-", estimated_MODU_G_effect$month, "-01"), 
                          label = TRUE, abbr = TRUE))

estimated_MODU_G_effect$month1 <- month_abbr

estimated_MODU_G_effect |> 
  group_by(month1) |> 
  transmute(mu = mean(est), 
            li = quantile(est, 0.025), 
            ls = quantile(est, 0.975)) |> 
  unique() |> 
  ggplot() +
  geom_ribbon(aes(month1, ymin = li, ymax = ls), alpha = 0.5) +
  geom_line(aes(month1, mu)) +
  scale_x_date(date_labels = '%b', date_breaks = '2 month')



cat(file = 'H2_global.stan',
    "

    functions {
      vector GP_periodic(int period,        // periodicity
                         real gamma,        // smoothing term of the GP
                         real sigma,        // noise paramether
                         vector eta) {      // latent variable for each month

                         int M = period;
                         matrix[M, M] K;
                         matrix[M, M] L_K;

                         for (i in 1:(M - 1)) {
                           for (j in (i+1):M) {
                               real distance = abs(i - j);
                               real periodic_distance = fmin(distance, period - distance);
                               K[i, j] = sigma^2 * exp(-2 * square(sin(pi()*periodic_distance/period))/gamma^2);
                               K[j, i] = K[i, j];   // filling the lower triangle
                            }
                             K[i, i] = sigma^2 + 1e-9;   // small values to guarante stability
                          }
                           K[M, M] = sigma^2 + 1e-9;    // small values to guarante stability
                           return cholesky_decompose(K) * eta;

                        }

      matrix GP_quadratic(matrix x,
                          real eta,
                          real rho,
                          real delta) {

                          int N = dims(x)[1];
                          matrix[N, N] K;
                          matrix[N, N] L_K;

                          for (i in 1:(N-1)) {
                            K[i, i] = eta + delta;
                            for (j in (i+1):N) {
                              K[i, j] = square(eta) * exp(-rho * square(x[i, j]));
                              K[j, i] = K[i, j];
                            }
                          }

                          K[N, N] = eta + delta;
                          L_K = cholesky_decompose(K);
                          return L_K;
                          }
    }

    data{
      int N;
      int N_sites;
      int N_month;
      array[N] int network_size;
      vector[N] modularity;
      vector[N] nestedness;
      vector[N] H2;
      array[N] int month;
      array[N] int site;
      vector[N] HV_network;
      vector[N] HV_plant;
      vector[N] HV_plants_morfo;
      vector[N] HV_plants_nut;
      vector[N] HV_bird;
      vector[N] z_temperature;
      vector[N] z_rainfall;
      vector[N] fruit_S;
      vector[N] fruit_abun;
      matrix[N_sites, N_sites] dist_sites;
    }

    parameters {

      //////////////////////////// Fruit abundance //////////////////////
      ////////////////////////////////////////////////////////////////

      // GP periodic
      real<lower = 0> gamma_fruitS;
      real<lower = 0> sigma_f_fruitS;
      vector[N_month] eta_fruitS;
      // GP quadratic
      vector[N_sites] z_sites_fruitS;
      real<lower = 0> eta_site_fruitS;
      real<lower = 0> rho_site_fruitS;

      // Pars linear model
      real beta_temp_fruitS;
      real beta_rain_fruitS;
      // noise for likelihood function
      real<lower = 0> sigma_fruitS;


      //////////////////////////// Fruit abundance //////////////////////
      ////////////////////////////////////////////////////////////////

      // GP periodic
      real<lower = 0> gamma_fruitAB;
      real<lower = 0> sigma_f_fruitAB;
      vector[N_month] eta_fruitAB;
      // GP quadratic
      vector[N_sites] z_sites_fruitAB;
      real<lower = 0> eta_site_fruitAB;
      real<lower = 0> rho_site_fruitAB;

      // Pars linear model

      real beta_temp_fruitAB;
      real beta_rain_fruitAB;
      // noise for likelihood function
      real<lower = 0> sigma_fruitAB;


      //////////////////////////// Fru. Nut. HV //////////////////////
      ////////////////////////////////////////////////////////////////

      // GP periodic
      real<lower = 0> gamma_nut_HV;
      real<lower = 0> sigma_f_nut_HV;
      vector[N_month] eta_nut_HV;
      // GP quadratic
      vector[N_sites] z_sites_nut_HV;
      real<lower = 0> eta_site_nut_hv;
      real<lower = 0> rho_site_nut_hv;

      // Pars linear model

      real beta_temp_nut_HV;
      real beta_rain_nut_HV;
      real beta_fruitAB_nut_HV;
      real beta_fruitS_nut_HV;
      real beta_BIRD_nut_HV;
      // noise for likelihood function
      real<lower = 0> sigma_nut_HV;

      //////////////////////////// Fru. morpho. HV //////////////////////
      ////////////////////////////////////////////////////////////////
      // GP periodic
      real<lower = 0> gamma_mor_HV;
      real<lower = 0> sigma_f_mor_HV;
      vector[N_month] eta_mor_HV;
      // GP quadratic
      vector[N_sites] z_sites_mor_HV;
      real<lower = 0> eta_site_mor_hv;
      real<lower = 0> rho_site_mor_hv;

      // Pars linear model

      real beta_temp_mor_HV;
      real beta_rain_mor_HV;
      real beta_fruitAB_mor_HV;
      real beta_fruitS_mor_HV;
      real beta_BIRD_mor_HV;
      // noise for likelihood function
      real<lower = 0> sigma_mor_HV;
      real nu;

      // //////////////////////////// Plants HV //////////////////////
      // ////////////////////////////////////////////////////////////////
      // GP periodic
      real<lower = 0> gamma_PLAN_HV;
      real<lower = 0> sigma_f_PLAN_HV;
      vector[N_month] eta_PLAN_HV;
      // GP quadratic
      vector[N_sites] z_sites_PLAN_HV;
      real<lower = 0> eta_site_PLAN_hv;
      real<lower = 0> rho_site_PLAN_hv;

      // Pars linear model
      real beta_temp_PLAN_HV;
      real beta_rain_PLAN_HV;
      real beta_fruitAB_PLAN_HV;
      real beta_fruitS_PLAN_HV;
      real beta_nut_PLAN_HV;
      real beta_mor_PLAN_HV;
      real beta_BIRD_PLAN_HV;
      // noise for likelihood function
      real<lower = 0> sigma_PLAN_HV;

      //////////////////////////// Birsds HV //////////////////////
      // ////////////////////////////////////////////////////////////////
      // GP periodic
      real<lower = 0> gamma_BIRD_HV;
      real<lower = 0> sigma_f_BIRD_HV;
      vector[N_month] eta_BIRD_HV;
      // GP quadratic
      vector[N_sites] z_sites_BIRD_HV;
      real<lower = 0> eta_site_BIRD_hv;
      real<lower = 0> rho_site_BIRD_hv;

      // Pars linear model
      real beta_temp_BIRD_HV;
      real beta_rain_BIRD_HV;
      real beta_fruitAB_BIRD_HV;
      real beta_fruitS_BIRD_HV;
      // noise for likelihood function
      real<lower = 0> sigma_BIRD_HV;

      //////////////////////////// Network metrics //////////////////////
      ////////////////////////////////////////////////////////////////
      // //////////////////////////// H2 //////////////////////
      // ////////////////////////////////////////////////////////////////
      // GP periodic
      real<lower = 0> gamma_H2;
      real<lower = 0> sigma_f_H2;
      vector[N_month] eta_H2;
      // GP quadratic
      vector[N_sites] z_sites_H2;
      real<lower = 0> eta_site_H2;
      real<lower = 0> rho_site_H2;

      // Pars linear model
      real beta_temp_H2;
      real beta_rain_H2;
      real beta_fruitAB_H2;
      real beta_fruitS_H2;
      real beta_nut_H2;
      real beta_mor_H2;
      real beta_PLAN_H2;
      real beta_BIRD_H2;
      // noise for likelihood function
      real<lower = 0> sigma_H2;

    }

    transformed parameters {

        //////////////////////////// Fruit abundance //////////////////////
    ////////////////////////////////////////////////////////////////
      // GP periodic
      vector[N_month] f_fruitAB;

      f_fruitAB = GP_periodic(
              12,
              gamma_fruitAB,
              sigma_f_fruitAB,
              eta_fruitAB
        );

      // GP quadratic
      vector[N_sites] alpha_fruitAB;
      matrix[N_sites, N_sites] L_K_fruitAB;
      L_K_fruitAB = GP_quadratic(dist_sites,
                                eta_site_fruitAB,
                                rho_site_fruitAB, 0.001);
      alpha_fruitAB = L_K_fruitAB * z_sites_fruitAB;


    //////////////////////////// Fruit richness //////////////////////
    ////////////////////////////////////////////////////////////////
      // GP periodic
      vector[N_month] f_fruitS;

      f_fruitS = GP_periodic(
              12,
              gamma_fruitS,
              sigma_f_fruitS,
              eta_fruitS
        );

      // GP quadratic
      vector[N_sites] alpha_fruitS;
      matrix[N_sites, N_sites] L_K_fruitS;
      L_K_fruitS = GP_quadratic(dist_sites,
                                eta_site_fruitS,
                                rho_site_fruitS, 0.001);
      alpha_fruitS = L_K_fruitS * z_sites_fruitS;

    //////////////////////////// Fru. Nut. HV //////////////////////
    ////////////////////////////////////////////////////////////////
      // GP periodic
      vector[N_month] f_nut_HV;

      f_nut_HV = GP_periodic(
              12,
              gamma_nut_HV,
              sigma_f_nut_HV,
              eta_nut_HV
        );

      // GP quadratic
      vector[N_sites] alpha_nut_HV;
      matrix[N_sites, N_sites] L_K_nut_HV;
      L_K_nut_HV = GP_quadratic(dist_sites,
                                eta_site_nut_hv,
                                rho_site_nut_hv, 0.001);
      alpha_nut_HV = L_K_nut_HV * z_sites_nut_HV;

    //   //////////////////////////// Fru. morpho. HV //////////////////////
    //   ////////////////////////////////////////////////////////////////
      // GP periodic
      vector[N_month] f_mor_HV;

      f_mor_HV = GP_periodic(
              12,
              gamma_mor_HV,
              sigma_f_mor_HV,
              eta_mor_HV
        );

      // GP quadratic
      vector[N_sites] alpha_mor_HV;
      matrix[N_sites, N_sites] L_K_mor_HV;
      L_K_mor_HV = GP_quadratic(dist_sites,
                                eta_site_mor_hv,
                                rho_site_mor_hv, 0.001);
      alpha_mor_HV = L_K_mor_HV * z_sites_mor_HV;
    //
    //   //////////////////////////// Plants HV //////////////////////
    //   ////////////////////////////////////////////////////////////////
      // GP periodic
      vector[N_month] f_PLAN_HV;

      f_PLAN_HV = GP_periodic(
              12,
              gamma_PLAN_HV,
              sigma_f_PLAN_HV,
              eta_PLAN_HV
        );

      // GP quadratic
      vector[N_sites] alpha_PLAN_HV;
      matrix[N_sites, N_sites] L_K_PLAN_HV;
      L_K_PLAN_HV = GP_quadratic(dist_sites,
                                eta_site_PLAN_hv,
                                rho_site_PLAN_hv, 0.001);
      alpha_PLAN_HV = L_K_PLAN_HV * z_sites_PLAN_HV;

      ////////////////////////// Birsds HV //////////////////////
      //////////////////////////////////////////////////////////////
      // GP periodic
      vector[N_month] f_BIRD_HV;

      f_BIRD_HV = GP_periodic(
              12,
              gamma_BIRD_HV,
              sigma_f_BIRD_HV,
              eta_BIRD_HV
        );

      // GP quadratic
      vector[N_sites] alpha_BIRD_HV;
      matrix[N_sites, N_sites] L_K_BIRD_HV;
      L_K_BIRD_HV = GP_quadratic(dist_sites,
                                eta_site_BIRD_hv,
                                rho_site_BIRD_hv, 0.001);
      alpha_BIRD_HV = L_K_BIRD_HV * z_sites_BIRD_HV;

      // //////////////////////////// H2 //////////////////////
      // ////////////////////////////////////////////////////////////////
      // GP periodic
      vector[N_month] f_H2;

      f_H2 = GP_periodic(
              12,
              gamma_H2,
              sigma_f_H2,
              eta_H2);

      // GP quadratic
      vector[N_sites] alpha_H2;
      matrix[N_sites, N_sites] L_K_H2;
      L_K_H2 = GP_quadratic(dist_sites,
                              eta_site_H2,
                              rho_site_H2, 0.001);
      alpha_H2 = L_K_H2 * z_sites_H2;

    }

    model {

              //////////////////////////// Fruit abundance //////////////////////
       ////////////////////////////////////////////////////////////////
       // priors for periodic GP
       eta_fruitAB ~ normal(0, 0.5);
       gamma_fruitAB ~ inv_gamma(5, 5);
       sigma_f_fruitAB ~ cauchy(0, 1);
       // GP quadratic
       z_sites_fruitAB ~ normal(0, 0.25);
       eta_site_fruitAB ~ exponential(4);
       rho_site_fruitAB ~ exponential(1);
       // Pars linear model
       beta_temp_fruitAB ~ normal(0, 0.5);
       beta_rain_fruitAB ~ normal(0, 0.5);
       // noise parametert of the likelihood
       sigma_fruitAB ~ exponential(1);

       //////////////////////////// Fruit richness //////////////////////
       ////////////////////////////////////////////////////////////////
       // priors for periodic GP
       eta_fruitS ~ normal(0, 0.5);
       gamma_fruitS ~ inv_gamma(5, 5);
       sigma_f_fruitS ~ cauchy(0, 1);
       // GP quadratic
       z_sites_fruitS ~ normal(0, 0.25);
       eta_site_fruitS ~ exponential(4);
       rho_site_fruitS ~ exponential(1);
       // Pars linear model
       beta_temp_fruitS ~ normal(0, 0.5);
       beta_rain_fruitS ~ normal(0, 0.5);
       // noise parametert of the likelihood
       sigma_fruitS ~ exponential(1);

       //////////////////////////// Fru. Nut. HV //////////////////////
       ////////////////////////////////////////////////////////////////
       // priors for periodic GP
       eta_nut_HV ~ normal(0, 0.5);
       gamma_nut_HV ~ inv_gamma(5, 5);
       sigma_f_nut_HV ~ cauchy(0, 1);
       // GP quadratic
       z_sites_nut_HV ~ normal(0, 0.25);
       eta_site_nut_hv ~ exponential(4);
       rho_site_nut_hv ~ exponential(1);
       // Pars linear model
       beta_temp_nut_HV ~ normal(0, 0.5);
       beta_rain_nut_HV ~ normal(0, 0.5);
       beta_fruitAB_nut_HV ~ normal(0, 0.5);
       beta_fruitS_nut_HV ~ normal(0, 0.5);
       beta_BIRD_nut_HV ~ normal(0, 0.5);
       // noise parametert of the likelihood
       sigma_nut_HV ~ exponential(1);
       //
       //////////////////////////// Fru. morpho. HV //////////////////////
       // ////////////////////////////////////////////////////////////////
       // priors for periodic GP
       eta_mor_HV ~ normal(0, 0.5);
       gamma_mor_HV ~ inv_gamma(5, 5);
       sigma_f_mor_HV ~ cauchy(0, 1);
       // GP quadratic
       z_sites_mor_HV ~ normal(0, 0.25);
       eta_site_mor_hv ~ exponential(4);
       rho_site_mor_hv ~ exponential(1);
       // Pars linear model
       beta_temp_mor_HV ~ normal(0, 0.5);
       beta_rain_mor_HV ~ normal(0, 0.5);
       beta_fruitAB_mor_HV ~ normal(0, 0.5);
       beta_fruitS_mor_HV ~ normal(0, 0.5);
       beta_BIRD_mor_HV ~ normal(0, 0.5);
       // noise parametert of the likelihood
       sigma_mor_HV ~ exponential(1);

       // //////////////////////////// Plants HV //////////////////////
       // ////////////////////////////////////////////////////////////////
       // priors for periodic GP
       eta_PLAN_HV ~ normal(0, 1);
       gamma_PLAN_HV ~ inv_gamma(5, 5);
       sigma_f_PLAN_HV ~ cauchy(0, 1);
       // GP quadratic
       z_sites_PLAN_HV ~ normal(0, 0.5);
       eta_site_PLAN_hv ~ exponential(4);
       rho_site_PLAN_hv ~ exponential(1);
       // Pars linear model
       beta_temp_PLAN_HV ~ normal(0, 0.5);
       beta_rain_PLAN_HV ~ normal(0, 0.5);
       beta_fruitAB_PLAN_HV ~ normal(0, 0.5);
       beta_fruitS_PLAN_HV ~ normal(0, 0.5);
       beta_nut_PLAN_HV ~ normal(0, 0.5);
       beta_mor_PLAN_HV ~ normal(0, 0.5);
       beta_BIRD_PLAN_HV ~ normal(0, 0.5);
       // noise parametert of the likelihood
       sigma_PLAN_HV ~ exponential(1);

       ////////////////////////// Birds HV //////////////////////
       //////////////////////////////////////////////////////////////
       // priors for periodic GP
       eta_BIRD_HV ~ normal(0, 0.5);
       gamma_BIRD_HV ~ inv_gamma(5, 5);
       sigma_f_BIRD_HV ~ cauchy(0, 1);
       // GP quadratic
       z_sites_BIRD_HV ~ normal(0, 0.5);
       eta_site_BIRD_hv ~ exponential(4);
       rho_site_BIRD_hv ~ exponential(1);
       // Pars linear model
       beta_temp_BIRD_HV ~ normal(0, 0.5);
       beta_rain_BIRD_HV ~ normal(0, 0.5);
       beta_fruitAB_BIRD_HV ~ normal(0, 0.5);
       beta_fruitS_BIRD_HV ~ normal(0, 0.5);
       // noise parametert of the likelihood
       sigma_BIRD_HV ~ exponential(1);

      // //////////////////////////// H2 //////////////////////
      // ////////////////////////////////////////////////////////////////
      vector[N] mu_H2;
      vector[N] p1_H2;
      vector[N] p2_H2;
      // priors for periodic GP
      eta_H2 ~ normal(0, 1);
      gamma_H2 ~ inv_gamma(5, 5);
      sigma_f_H2 ~ cauchy(0, 1);
      // GP quadratic
      z_sites_H2 ~ normal(0, 0.25);
      eta_site_H2 ~ exponential(4);
      rho_site_H2 ~ exponential(1);
      // Pars linear model
      beta_fruitAB_H2 ~ normal(0, 0.055);
      beta_fruitS_H2 ~ normal(0, 0.055);
      beta_nut_H2 ~ normal(0, 0.055);
      beta_mor_H2 ~ normal(0, 0.055);
      beta_PLAN_H2 ~ normal(0, 0.055);
      beta_BIRD_H2 ~ normal(0, 0.055);
      // beta_NW_H2 ~ normal(0, 0.055);
      // noise parametert of the likelihood
      sigma_H2 ~ exponential(2);

      //////////////////////////// fruit abundance //////////////////////
       ////////////////////////////////////////////////////////////////
       for (i in 1:N) {
         fruit_abun[i] ~ student_t(7, alpha_fruitAB[site[i]] + f_fruitAB[month[i]] +
                                    beta_temp_fruitAB * z_temperature[i] +
                                    beta_rain_fruitAB * z_rainfall[i],
                                    sigma_fruitAB);
       }

      //////////////////////////// fruit richness //////////////////////
       ////////////////////////////////////////////////////////////////
       for (i in 1:N) {
         fruit_S[i] ~ student_t(7, alpha_fruitS[site[i]] + f_fruitS[month[i]] +
                                    beta_temp_fruitS * z_temperature[i] +
                                    beta_rain_fruitS * z_rainfall[i],
                                    sigma_fruitS);
       }

       //////////////////////////// Fru. Nut. HV //////////////////////
       ////////////////////////////////////////////////////////////////
       for (i in 1:N) {
         HV_plants_nut[i] ~ student_t(7, alpha_nut_HV[site[i]] + f_nut_HV[month[i]] +
                                   beta_temp_nut_HV * z_temperature[i] +
                                   beta_rain_nut_HV * z_rainfall[i] +
                                   beta_BIRD_nut_HV * HV_bird[i] +
                                   beta_fruitAB_nut_HV * fruit_abun[i] +
                                   beta_fruitS_nut_HV * fruit_S[i],
                                   sigma_nut_HV);
       }
       //
       // //////////////////////////// Fru. morpho. HV //////////////////////
       // ////////////////////////////////////////////////////////////////
       for (i in 1:N) {
         HV_plants_morfo[i] ~ student_t(15, alpha_mor_HV[site[i]] + f_mor_HV[month[i]] +
                                     beta_temp_mor_HV * z_temperature[i] +
                                     beta_rain_mor_HV * z_rainfall[i] +
                                     beta_fruitAB_mor_HV * fruit_abun[i] +
                                     beta_fruitS_mor_HV * fruit_S[i] +
                                     beta_BIRD_mor_HV * HV_bird[i],
                                     sigma_mor_HV);
       }

       // //////////////////////////// Plants HV //////////////////////
       // ////////////////////////////////////////////////////////////////
       for (i in 1:N) {
         HV_plant[i] ~ student_t(7, alpha_PLAN_HV[site[i]] + f_PLAN_HV[month[i]] +
                                    beta_temp_PLAN_HV * z_temperature[i] +
                                    beta_rain_PLAN_HV * z_rainfall[i] +
                                    beta_fruitAB_PLAN_HV * fruit_abun[i] +
                                    beta_fruitS_PLAN_HV * fruit_S[i] +
                                    beta_nut_PLAN_HV * HV_plants_nut[i] +
                                    beta_mor_PLAN_HV * HV_plants_morfo[i] +
                                    beta_BIRD_PLAN_HV * HV_bird[i],
                                    sigma_PLAN_HV);
       }

       ////////////////////////// Birds HV //////////////////////
       //////////////////////////////////////////////////////////////
       for (i in 1:N) {
         HV_bird[i] ~ student_t(7, alpha_BIRD_HV[site[i]] + f_BIRD_HV[month[i]] +
                                    beta_temp_BIRD_HV * z_temperature[i] +
                                    beta_rain_BIRD_HV * z_rainfall[i] +
                                    beta_fruitAB_BIRD_HV * fruit_abun[i] +
                                    beta_fruitS_BIRD_HV * fruit_S[i],
                                    sigma_BIRD_HV);
       }
      // //////////////////////////// H2 //////////////////////
      // ////////////////////////////////////////////////////////////////
      for (i in 1:N) {
        mu_H2[i] = inv_logit(alpha_H2[site[i]] + f_H2[month[i]] +
                     beta_temp_H2 * z_temperature[i] +
                     beta_rain_H2 * z_rainfall[i] +
                     beta_fruitAB_H2 * fruit_abun[i] +
                     beta_fruitS_H2 * fruit_S[i] +
                     beta_nut_H2 * HV_plants_nut[i] +
                     beta_mor_H2 * HV_plants_morfo[i] +
                     beta_PLAN_H2 * HV_plant[i] +
                     beta_BIRD_H2 * HV_bird[i]);
      }

      p1_H2 = mu_H2 * sigma_H2;
      p2_H2 = (1 - mu_H2) * sigma_H2;

      H2 ~ beta(p1_H2, p2_H2);

    }

    generated quantities {

       //////////////////////////// Fruit abundance //////////////////////
       ////////////////////////////////////////////////////////////////
       array[N] real ppcheck_fruitAB;
       vector[N] mu_fruitAB;

       for (i in 1:N) {
         mu_fruitAB[i] = alpha_fruitAB[site[i]] + f_fruitAB[month[i]] +
                         beta_temp_fruitAB * z_temperature[i] +
                         beta_rain_fruitAB * z_rainfall[i];
       }

       ppcheck_fruitAB = student_t_rng(7, mu_fruitAB, sigma_fruitAB);


       //////////////////////////// Fruit richness //////////////////////
       ////////////////////////////////////////////////////////////////
       array[N] real ppcheck_fruitS;
       vector[N] mu_fruitS;

       for (i in 1:N) {
         mu_fruitS[i] = alpha_fruitS[site[i]] + f_fruitS[month[i]] +
                        beta_temp_fruitS * z_temperature[i] +
                        beta_rain_fruitS * z_rainfall[i];
       }

       ppcheck_fruitS = student_t_rng(7, mu_fruitS, sigma_fruitS);


       //////////////////////////// Fru. Nut. HV //////////////////////
       ////////////////////////////////////////////////////////////////
       array[N] real ppcheck_nut_HV;
       vector[N] mu_nut_HV;

       for (i in 1:N) {
         mu_nut_HV[i] = alpha_nut_HV[site[i]] + f_nut_HV[month[i]] +
                        beta_temp_nut_HV * z_temperature[i] +
                        beta_rain_nut_HV * z_rainfall[i] +
                        beta_fruitAB_nut_HV * fruit_abun[i] +
                        beta_fruitS_nut_HV * fruit_S[i] +
                        beta_BIRD_nut_HV * HV_bird[i];
       }

       ppcheck_nut_HV = student_t_rng(7, mu_nut_HV, sigma_nut_HV);
       //
       // //////////////////////////// Fru. morpho. HV //////////////////////
       // ////////////////////////////////////////////////////////////////
       array[N] real ppcheck_mor_HV;
       vector[N] mu_mor_HV;

       for (i in 1:N) {
         mu_mor_HV[i] = alpha_mor_HV[site[i]] + f_mor_HV[month[i]] +
                        beta_temp_mor_HV * z_temperature[i] +
                        beta_rain_mor_HV * z_rainfall[i] +
                        beta_fruitAB_mor_HV * fruit_abun[i] +
                        beta_fruitS_mor_HV * fruit_S[i] +
                        beta_BIRD_mor_HV * HV_bird[i];
       }

       ppcheck_mor_HV = student_t_rng(15, mu_mor_HV, sigma_mor_HV);

       // //////////////////////////// Plants HV //////////////////////
       // ////////////////////////////////////////////////////////////////
       array[N] real ppcheck_PLAN_HV;
       vector[N] mu_PLAN_HV;

       for (i in 1:N) {
         mu_PLAN_HV[i] = alpha_PLAN_HV[site[i]] + f_PLAN_HV[month[i]] +
                         beta_temp_PLAN_HV * z_temperature[i] +
                         beta_rain_PLAN_HV * z_rainfall[i] +
                         beta_fruitAB_PLAN_HV * fruit_abun[i] +
                         beta_fruitS_PLAN_HV * fruit_S[i] +
                         beta_nut_PLAN_HV * HV_plants_nut[i] +
                         beta_mor_PLAN_HV * HV_plants_morfo[i] +
                         beta_BIRD_PLAN_HV * HV_bird[i];
       }

       ppcheck_PLAN_HV = student_t_rng(7, mu_PLAN_HV, sigma_PLAN_HV);

       ////////////////////////// Birds HV //////////////////////
       //////////////////////////////////////////////////////////////
       array[N] real ppcheck_BIRD_HV;
       vector[N] mu_BIRD_HV;

       for (i in 1:N) {
         mu_BIRD_HV[i] = alpha_BIRD_HV[site[i]] + f_BIRD_HV[month[i]] +
                         beta_temp_BIRD_HV * z_temperature[i] +
                         beta_rain_BIRD_HV * z_rainfall[i] +
                         beta_fruitAB_BIRD_HV * fruit_abun[i] +
                         beta_fruitS_BIRD_HV * fruit_S[i];
       }

       ppcheck_BIRD_HV = student_t_rng(7, mu_BIRD_HV, sigma_BIRD_HV);


      // ////////////////////////// H2 //////////////////////
      // //////////////////////////////////////////////////////////////
      array[N] real ppcheck_H2;
      vector[N] mu_H2;
      vector[N] p1_H2;
      vector[N] p2_H2;
      real variance_fixed;
      real variance_GP_temp;
      real variance_GP_space;
      real variance_total;
      real prop_GP_temp;
      real prop_GP_space;
      real prop_fixed;

      for (i in 1:N) {
        mu_H2[i] = alpha_H2[site[i]] + f_H2[month[i]] +
                    beta_temp_H2 * z_temperature[i] +
                    beta_rain_H2 * z_rainfall[i] +
                    beta_fruitAB_H2 * fruit_abun[i] +
                    beta_fruitS_H2 * fruit_S[i] +
                    beta_nut_H2 * HV_plants_nut[i] +
                    beta_mor_H2 * HV_plants_morfo[i] +
                    beta_PLAN_H2 * HV_plant[i] +
                    beta_BIRD_H2 * HV_bird[i];
      }

      // variance decomposition

      variance_fixed = variance(mu_H2);
      variance_GP_space = eta_site_H2^2;
      variance_GP_temp = sigma_f_H2^2;
      variance_total = variance_fixed + variance_GP_temp +
                       variance_GP_space + ((pi()^2)/3);
      prop_GP_temp = variance_GP_temp / variance_total;
      prop_GP_space = variance_GP_space / variance_total;
      prop_fixed = variance_fixed / variance_total;

      mu_H2 = inv_logit(mu_H2);
      p1_H2 = mu_H2 * sigma_H2;
      p2_H2 = (1 - mu_H2) * sigma_H2;

      ppcheck_H2 = beta_rng(p1_H2, p2_H2);

    }
    "
    )


file <- paste0(getwd(), '/H2_global.stan')
fit_H2_G <- cmdstan_model(file, compile = T)

mod_H2_G <-
  fit_H2_G$sample(
    data = dat,
    iter_sampling = 2e3,
    iter_warmup = 500,
    thin = 3,
    chains = 3,
    parallel_chains = 3,
    seed = 5
  )


summary_H2_G <- mod_H2_G$summary()

par(mfrow = c(3, 3), mar = c(4, 4, 2, 1))

mod_diagnostics(mod_H2_G, 
                summary_H2_G[grep('fruitAB', summary_H2_G$variable), ])
title('Fruit abundance')
mod_diagnostics(mod_H2_G, 
                summary_H2_G[grep('fruitS', summary_H2_G$variable), ])
title('Fruit richness')
mod_diagnostics(mod_H2_G, 
                summary_H2_G[grep('nut_HV', summary_H2_G$variable), ])
title('Frut. Nut. HV')
mod_diagnostics(mod_H2_G, 
                summary_H2_G[grep('mor_HV', summary_H2_G$variable), ])
title('Frut. Morpho. HV')
mod_diagnostics(mod_H2_G, 
                summary_H2_G[grep('PLAN_HV', summary_H2_G$variable), ])
title('Plants HV')
mod_diagnostics(mod_H2_G, 
                summary_H2_G[grep('BIRD_HV', summary_H2_G$variable), ])
title('Birds HV')
mod_diagnostics(mod_H2_G,
                summary_H2_G[grep('H2', summary_H2_G$variable), ])
title('H2')
par(mfrow = c(1, 1))


par(mfrow = c(3, 3), mar = c(4, 4, 1, 1))
ppcheck_H2_G <- mod_H2_G$draws('ppcheck_fruitAB', format = 'matrix')
plot(density(dat$fruit_abun), main = '', xlab = 'Fruit abundance',
     ylim = c(0, 1))
for (i in 1:100) lines(density(ppcheck_H2_G[i, ]), lwd = 0.1)
lines(density(dat$fruit_abun), col = 'red', lwd = 2)

ppcheck_H2_G <- mod_H2_G$draws('ppcheck_fruitS', format = 'matrix')
plot(density(dat$fruit_S), main = '', xlab = 'Fruit richness',
     ylim = c(0, 1))
for (i in 1:100) lines(density(ppcheck_H2_G[i, ]), lwd = 0.1)
lines(density(dat$fruit_S), col = 'red', lwd = 2)
ppcheck_H2_G <- mod_H2_G$draws('ppcheck_nut_HV', format = 'matrix')
plot(density(dat$HV_plants_nut), main = '', xlab = 'Frut. Nut. HV',
     ylim = c(0, 1))
for (i in 1:100) lines(density(ppcheck_H2_G[i, ]), lwd = 0.1)
lines(density(dat$HV_plants_nut), col = 'red', lwd = 2)

ppcheck_H2_G <- mod_H2_G$draws('ppcheck_mor_HV', format = 'matrix')
plot(density(dat$HV_plants_morfo), main = '', xlab = 'Frut. Morpho. HV',
     ylim = c(0, 1))
for (i in 1:100) lines(density(ppcheck_H2_G[i, ]), lwd = 0.1)
lines(density(dat$HV_plants_morfo), col = 'red', lwd = 2)

ppcheck_H2_G <- mod_H2_G$draws('ppcheck_PLAN_HV', format = 'matrix')
plot(density(dat$HV_plant), main = '', xlab = 'Plants HV',
     ylim = c(0, 1))
for (i in 1:100) lines(density(ppcheck_H2_G[i, ]), lwd = 0.1)
lines(density(dat$HV_plant), col = 'red', lwd = 2)

ppcheck_H2_G <- mod_H2_G$draws('ppcheck_BIRD_HV', format = 'matrix')
plot(density(dat$HV_bird), main = '', xlab = 'Birds HV',
     ylim = c(0, 1))
for (i in 1:100) lines(density(ppcheck_H2_G[i, ]), lwd = 0.1)
lines(density(dat$HV_bird), col = 'red', lwd = 2)

ppcheck_H2_G <- mod_H2_G$draws('ppcheck_H2', format = 'matrix')
plot(density(dat$H2), main = '', xlab = 'H2',
     ylim = c(0, 2))
for (i in 1:100) lines(density(ppcheck_H2_G[i, ]), lwd = 0.1)
lines(density(dat$H2), col = 'red', lwd = 2)
par(mfrow = c(1, 1))


pars_H2_G <- summary_H2_G[grep('^beta', summary_H2_G$variable), ]$variable

betas_H2_G <- mod_H2_G$draws(pars_H2_G, format = 'df')

betas_H2_G <- betas_H2_G[, -c((ncol(betas_H2_G)-2):ncol(betas_H2_G))]

betas_H2_G <- 
  lapply(1:ncol(betas_H2_G), 
         function(i) {
           x <- betas_H2_G[[i]]
           var <- colnames(betas_H2_G)[[i]]
           
           df <- 
             tibble(var = var, 
                    mu = mean(x), 
                    li = quantile(x, 0.025), 
                    ls = quantile(x, 0.975),
                    p_pos = mean(x > 0), 
                    p_neg = mean(x < 0))
           
           df$effect <- df$p_pos > 0.7 | df$p_neg > 0.7
           
           df
           
         })


betas_H2_G <- do.call('rbind', betas_H2_G)

betas_H2_G_LABS <- betas_H2_G[c(grep('PLAN_H2', betas_H2_G$var), 
                            grep('mor_H2', betas_H2_G$var), 
                            grep('nut_H2', betas_H2_G$var)), ]

betas_H2_GEFFECT <- betas_H2_G[grep('H2', betas_H2_G$var), ]


betas_H2_G_z <- betas_H2_G

betas_H2_G_z$var <-
  gsub('^(beta_)([A-Za-z]*)(_)(.*)$', '\\2 -> \\4', betas_H2_G_z$var)

colnames(betas_H2_G_z) <-
  c('Directional effect', 'Mean', '0.025%', '0.975%', 'P(effect > 0)',
    'P(effect < 0)', 'effect')

knitr::kable(betas_H2_G_z, 
             digits = 2, 
             caption = 'Average values, 95% CIs, and probabilities of directional effects for the relationships assessed in the generative model. An effect is marked as TRUE if P(effect > 0) > 0.7 or P(effect < 0) > 0.7 or.')


pars_H2_G <- 
  summary_H2_G[c(grep('^beta', summary_H2_G$variable), 
               grep('^alpha', summary_H2_G$variable), 
               grep('^f', summary_H2_G$variable), 
               grep('^sigma', summary_H2_G$variable)), ]$variable

post_H2_G <- mod_H2_G$draws(pars_H2_G, format = 'df')

post_H2_G <- 
  lapply(c('^beta', '^alpha', '^f', '^sigma'), FUN = 
           function(x) {
             post_H2_G[, grep(x, colnames(post_H2_G))]
           })

names(post_H2_G) <- c('beta', 'alpha', 'f', 'sigma')



estimated_H2_G_effect <- 
  lapply(1, FUN = 
           function(x) {
             
             
             f_ <- 'f_H2['
             
             s <- grep('alpha_H2', colnames(post_H2_G$alpha))
             
             mu_H2_G <- mean(dat$H2)
             sd_H2_G <- sd(dat$H2)
             
             est <- lapply(1:12, FUN =
                             function(j) {
                               
                               f_ <- paste0(f_, j, ']')
                               
                               m <-
                                 apply(post_H2_G$alpha[, s], 1, mean) +
                                 post_H2_G$f[[f_]] +
                                 post_H2_G$beta$beta_temp_H2 * mean(dat$z_temperature) +
                                 post_H2_G$beta$beta_rain_H2 * mean(dat$z_rainfall) +
                                 post_H2_G$beta$beta_fruitAB_H2 * mean(dat$fruit_abun) +
                                 post_H2_G$beta$beta_fruitS_H2 * mean(dat$fruit_S) +
                                 post_H2_G$beta$beta_nut_H2 * mean(dat$HV_plants_nut) +
                                 post_H2_G$beta$beta_mor_H2 * mean(dat$HV_plants_morfo) +
                                 post_H2_G$beta$beta_PLAN_H2 * mean(dat$HV_plant) +
                                 post_H2_G$beta$beta_BIRD_H2 * mean(dat$HV_bird)
                               
                               m <- inv_logit(m)
                               
                               set.seed(123)
                               m <- sample(m, 1e3)
                               
                               tibble(est = (m - mu_H2_G) / sd_H2_G, 
                                      month = j, 
                                      class = 'H2')
                               
                             })
             est <- do.call('rbind', est)
             
             est
             
           })

estimated_H2_G_effect <- 
  do.call('rbind', estimated_H2_G_effect)


month_abbr <- na.omit(ymd(paste0("2023-", estimated_H2_G_effect$month, "-01"), 
                          label = TRUE, abbr = TRUE))

estimated_H2_G_effect$month1 <- month_abbr


estimated_H2_G_effect |> 
  group_by(month1) |> 
  transmute(mu = mean(est), 
            li = quantile(est, 0.025), 
            ls = quantile(est, 0.975)) |> 
  unique() |> 
  ggplot() +
  geom_ribbon(aes(month1, ymin = li, ymax = ls), alpha = 0.5) +
  geom_line(aes(month1, mu)) +
  scale_x_date(date_labels = '%b', date_breaks = '2 month')



global_climate <- readRDS('climate_for_plot.rds')

global_climate$class <- rep(c('Rainfall', 'Temperature'), each = 12)

global_climate$month1 <- ymd(paste0('2025-',
                                    global_climate$month, 
                                    '-01'))

plot_global_climate <- 
  global_climate |> 
  mutate(type = 'Climate') |> 
  ggplot(aes(month1, mu, ymin = li_mu, ymax = ls_mu)) +
  geom_line(aes(color = class), linewidth = 0.7, alpha = 0.5) +
  # geom_point(data = 
  #              rbind(global_climate[global_climate$class == 'Temperature', ] |> 
  #                      filter(mu == min(mu) | mu == max(mu)),
  #                    global_climate[global_climate$class == 'Rainfall', ] |> 
  #                      filter(mu == min(mu) | mu == max(mu))), 
  #           mapping = aes(month1, mu, color = class)) +
  geom_ribbon(aes(fill = class), alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = 3, linewidth = 0.25) +
  geom_vline(xintercept = c(ymd(c('2025-04-01', '2025-11-01'))), 
             linetype = 2, linewidth = 0.3) +
  annotate(geom = 'text', label = 'Rainy season', 
           x = ymd('2025-07-15'), y = -2) +
  scale_x_date(date_labels = '%b', date_breaks = '2 month') +
  labs(y = 'z-scores', x = NULL) +
  theme_classic() + 
  facet_wrap(~ type) +
  scale_color_manual(values = c("steelblue", "tan1")) +
  scale_fill_manual(values = c("steelblue", "tan1")) +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.5, 0.95),
        legend.direction = 'horizontal',
        strip.background = element_blank(),
        axis.line = element_line(linewidth = 0.25),
        axis.ticks = element_line(linewidth = 0.25),
        #text = element_text(family = 'Times New Roman'),
        legend.background = element_blank(),
        #axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text = element_text(size = 8.5),
        #title = element_blank(),
        legend.margin = margin(c(0, 0, 0, 0)))


plot_global_climate

global_phenology <- 
  lapply(1:2, FUN = 
         function(x) {
           
           
           
           f_ <- c('f_fruitAB[', 'f_fruitS[')
           
           alpha <- c('alpha_fruitAB', 'alpha_fruitS')
           
           s <- grep(alpha[x], colnames(post_H2_G$alpha))
           
           est <- lapply(1:12, FUN =
                           function(j) {
                             
                             class = c('Fruit production', 
                                       'Fruit richness')
                             
                             f_ <- paste0(f_[x], j, ']')
                             
                             m <-
                               apply(post_H2_G$alpha[, s], 1, mean) +
                               post_H2_G$f[[f_]] +
                               post_H2_G$beta$beta_temp_fruitAB * mean(dat$z_temperature) +
                               post_H2_G$beta$beta_rain_fruitAB * mean(dat$z_rainfall) 
                               # post_H2_G$beta$beta_fruitAB_H2 * mean(dat$fruit_abun) +
                               # post_H2_G$beta$beta_fruitS_H2 * mean(dat$fruit_S) +
                               # post_H2_G$beta$beta_nut_H2 * mean(dat$HV_plants_nut) +
                               # post_H2_G$beta$beta_mor_H2 * mean(dat$HV_plants_morfo) +
                               # post_H2_G$beta$beta_PLAN_H2 * mean(dat$HV_plant) +
                               # post_H2_G$beta$beta_BIRD_H2 * mean(dat$HV_bird)
                             set.seed(123)
                             m <- sample(m, 1e3)
                             
                             tibble(est = m, 
                                    month = j, 
                                    class = class[x])
                             
                           })
           est <- do.call('rbind', est)
           
           est
           
         })

global_phenology <- 
  do.call('rbind', global_phenology)


month_abbr <- na.omit(ymd(paste0("2023-", global_phenology$month, "-01"), 
                          label = TRUE, abbr = TRUE))

global_phenology$month1 <- month_abbr

global_phenology |> 
  group_by(class, month1) |> 
  transmute(mu = mean(est), 
            li = quantile(est, 0.1), 
            ls = quantile(est, 0.9), 
            type = 'Fruiting phenology', 
            month = month) |> 
  unique() |> 
  filter(class == 'Fruit richness')

dat_mod_pheno <- readRDS('data_phenology_model.rds')

plot_global_phenology <- 
  global_phenology |> 
  group_by(class, month1) |> 
  transmute(mu = mean(est), 
            li = quantile(est, 0.1), 
            ls = quantile(est, 0.9), 
            type = 'Fruiting phenology') |> 
  unique() |> 
  ggplot(aes(month1, mu, ymin = li, ymax = ls)) +
  geom_line(aes(color = class), linewidth = 0.7, alpha = 0.5) +
  geom_ribbon(aes(fill = class), alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = 3, linewidth = 0.25) +
  geom_vline(xintercept = c(ymd(c('2023-04-01', '2023-11-01'))), 
             linetype = 2, linewidth = 0.3) +
  scale_x_date(date_labels = '%b', date_breaks = '2 month') +
  labs(y = 'z-scores', x = NULL) +
  theme_classic() + 
  facet_wrap(~ type) +
  scale_color_manual(values = c("tan1", "steelblue")) +
  scale_fill_manual(values = c("tan1", "steelblue")) +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.5, 0.95),
        legend.direction = 'horizontal',
        strip.background = element_blank(),
        axis.line = element_line(linewidth = 0.25),
        axis.ticks = element_line(linewidth = 0.25),
        #text = element_text(family = 'Times New Roman'),
        #axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text = element_text(size = 8.5),
        #title = element_blank(),
        legend.margin = margin(c(0, 0, 0, 0)))


plot_global_phenology

## species richness

# maximum species number
max_phenologyS_std <- 
  global_phenology |> 
    filter(class == 'Fruit richness' &
             month >= 7 & month <= 8) # July and August

max_phenologyS <- 
  max_phenologyS_std$est + (mean(d$fruit_S) * sd(d$fruit_S))

median(max_phenologyS) # average maximum 
quantile(max_phenologyS, 0.025) # low boundary CI
quantile(max_phenologyS, 0.975) # upper boundary CI

min_phenologyS_std <- 
  global_phenology |> 
  filter(month >= 2 & month <=3) |> 
    filter(class == 'Fruit richness')

min_phenologyS <- 
  min_phenologyS_std$est + mean(d$fruit_S) * sd(d$fruit_S)

median(min_phenologyS) # average maximum 
quantile(min_phenologyS, 0.025) # low boundary CI
quantile(min_phenologyS, 0.975) # upper boundary CI

set.seed(1)
mean(max_phenologyS > # P(max > high)
       sample(min_phenologyS, length(max_phenologyS)))



# maximum fruit production
max_phenologyABU_std <- 
  global_phenology |> 
  filter(class == 'Fruit production' &
           month >= 11 & month <= 12) # November and December

max_phenologyABU <- 
  exp(max_phenologyABU_std$est + 
  (mean(dat_mod_pheno$abun) * sd(dat_mod_pheno$abun)))

median(max_phenologyABU) # average maximum
quantile(max_phenologyABU, 0.025) # lower boundary CI
quantile(max_phenologyABU, 0.975) # upper boundary CI

min_phenologyABU_std <- 
  global_phenology |> 
  filter(month >= 5 & month <=6) |> 
  filter(class == 'Fruit production')

min_phenologyABU <- 
  exp(min_phenologyABU_std$est + 
  (mean(dat_mod_pheno$abun) * sd(dat_mod_pheno$abun)))

median(min_phenologyABU) # Average minimum
quantile(min_phenologyABU, 0.025) # lower boundary CI
quantile(min_phenologyABU, 0.975) # upper boundary CI

tot_fru_increase <- 100 - (min_phenologyABU * 100) / max_phenologyABU

median(tot_fru_increase) # Average contrast
quantile(tot_fru_increase, 0.025) # lower boundary CI
quantile(tot_fru_increase, 0.975) # upper boundary CI
mean(max_phenologyABU > min_phenologyABU) # P(higher > lower)


global_hypervolumes <- 
  lapply(1:4, FUN = 
           function(x) {
             
             class <- c("Fruits' nutr. hypervolume", 
                        "Fruits' morph. hypervolume", 
                        "Fruits' total hypervolume", 
                        "Birds' hypervolume")
             
             f_ <- c('f_nut_HV[', 'f_mor_HV[', 'f_PLAN_HV[', 'f_BIRD_HV[')
             
             alpha <- c('alpha_nut_HV', 'alpha_mor_HV', 
                        'alpha_PLAN_HV', 'alpha_BIRD_HV')
             
             s <- grep(alpha[x], colnames(post_H2_G$alpha))
             
             if (f_[x] == 'f_nut_HV[') {
               
               est <- lapply(1:12, FUN =
                               function(j) {
                                 
                                 f_ <- paste0(f_[x], j, ']')
                                 
                                 m <-
                                   apply(post_H2_G$alpha[, s], 1, mean) +
                                   post_H2_G$f[[f_]] +
                                   post_H2_G$beta$beta_temp_nut_HV * mean(dat$z_temperature) +
                                   post_H2_G$beta$beta_rain_nut_HV * mean(dat$z_rainfall) +
                                   post_H2_G$beta$beta_fruitAB_nut_HV * mean(dat$fruit_abun) +
                                   post_H2_G$beta$beta_fruitS_nut_HV * mean(dat$fruit_S) +
                                   #post_H2_G$beta$beta_nut_H2 * mean(dat$HV_plants_nut) +
                                   # post_H2_G$beta$beta_mor_H2 * mean(dat$HV_plants_morfo) +
                                   # post_H2_G$beta$beta_PLAN_H2 * mean(dat$HV_plant) +
                                   post_H2_G$beta$beta_BIRD_nut_HV * mean(dat$HV_bird)
                                 set.seed(123)
                                 m <- sample(m, 1e3)
                                 
                                 tibble(est = m, 
                                        month = j, 
                                        class = class[x])
                                 
                               })
               est <- do.call('rbind', est)
               
               est
             } else if (f_[x] == 'f_mor_HV[') {
               
               est <- lapply(1:12, FUN =
                               function(j) {
                                 
                                 f_ <- paste0(f_[x], j, ']')
                                 
                                 m <-
                                   apply(post_H2_G$alpha[, s], 1, mean) +
                                   post_H2_G$f[[f_]] +
                                   post_H2_G$beta$beta_temp_mor_HV * mean(dat$z_temperature) +
                                   post_H2_G$beta$beta_rain_mor_HV * mean(dat$z_rainfall) +
                                   post_H2_G$beta$beta_fruitAB_mor_HV * mean(dat$fruit_abun) +
                                   post_H2_G$beta$beta_fruitS_mor_HV * mean(dat$fruit_S) +
                                   #post_H2_G$beta$beta_nut_H2 * mean(dat$HV_plants_nut) +
                                   # post_H2_G$beta$beta_mor_H2 * mean(dat$HV_plants_morfo) +
                                   # post_H2_G$beta$beta_PLAN_H2 * mean(dat$HV_plant) +
                                   post_H2_G$beta$beta_BIRD_mor_HV * mean(dat$HV_bird)
                                 set.seed(123)
                                 m <- sample(m, 1e3)
                                 
                                 tibble(est = m, 
                                        month = j, 
                                        class = class[x])
                                 
                               })
               est <- do.call('rbind', est)
               
               est
             } else if (f_[x] == 'f_PLAN_HV[') {
               
               est <- lapply(1:12, FUN =
                               function(j) {
                                 
                                 f_ <- paste0(f_[x], j, ']')
                                 
                                 m <-
                                   apply(post_H2_G$alpha[, s], 1, mean) +
                                   post_H2_G$f[[f_]] +
                                   post_H2_G$beta$beta_temp_PLAN_HV * mean(dat$z_temperature) +
                                   post_H2_G$beta$beta_rain_PLAN_HV * mean(dat$z_rainfall) +
                                   post_H2_G$beta$beta_fruitAB_PLAN_HV * mean(dat$fruit_abun) +
                                   post_H2_G$beta$beta_fruitS_PLAN_HV * mean(dat$fruit_S) +
                                   post_H2_G$beta$beta_nut_PLAN_HV * mean(dat$HV_plants_nut) +
                                   post_H2_G$beta$beta_mor_PLAN_HV * mean(dat$HV_plants_morfo) +
                                   #post_H2_G$beta$beta_PLAN_HV * mean(dat$HV_plant) +
                                   post_H2_G$beta$beta_BIRD_PLAN_HV * mean(dat$HV_bird)
                                 set.seed(123)
                                 m <- sample(m, 1e3)
                                 
                                 tibble(est = m, 
                                        month = j, 
                                        class = class[x])
                                 
                               })
               est <- do.call('rbind', est)
               
               est
             } else if (f_[x] == 'f_BIRD_HV[') {
               
               est <- lapply(1:12, FUN =
                               function(j) {
                                 
                                 f_ <- paste0(f_[x], j, ']')
                                 
                                 m <-
                                   apply(post_H2_G$alpha[, s], 1, mean) +
                                   post_H2_G$f[[f_]] +
                                   post_H2_G$beta$beta_temp_BIRD_HV * mean(dat$z_temperature) +
                                   post_H2_G$beta$beta_rain_BIRD_HV * mean(dat$z_rainfall) +
                                   post_H2_G$beta$beta_fruitAB_BIRD_HV * mean(dat$fruit_abun) +
                                   post_H2_G$beta$beta_fruitS_BIRD_HV * mean(dat$fruit_S) 
                                   # post_H2_G$beta$beta_nut_PLAN_HV * mean(dat$HV_plants_nut) +
                                   # post_H2_G$beta$beta_mor_PLAN_HV * mean(dat$HV_plants_morfo) +
                                   #post_H2_G$beta$beta_PLAN_HV * mean(dat$HV_plant) +
                                   #post_H2_G$beta$beta_BIRD_PLAN_HV * mean(dat$HV_bird)
                                 set.seed(123)
                                 m <- sample(m, 1e3)
                                 
                                 tibble(est = m, 
                                        month = j, 
                                        class = class[x])
                                 
                               })
               est <- do.call('rbind', est)
               
               est
             }
             
           })

global_hypervolumes <- 
  do.call('rbind', global_hypervolumes)


month_abbr <- na.omit(ymd(paste0("2023-", global_hypervolumes$month, "-01"), 
                          label = TRUE, abbr = TRUE))

global_hypervolumes$month1 <- month_abbr

plot_global_hypervolume <- 
  lapply(unique(global_hypervolumes$class), FUN = 
           function(x) {
             p <- 
             global_hypervolumes |> 
               filter(class == x) |> 
               group_by(class, month1) |> 
               transmute(mu = mean(est), 
                         li = quantile(est, 0.1), 
                         ls = quantile(est, 0.9)) |> 
               unique() |> 
               ggplot(aes(month1, mu, ymin = li, ymax = ls)) +
               geom_line(aes(), linewidth = 0.7, alpha = 0.5, color = 'lightblue3') +
               geom_ribbon(aes(), alpha = 0.3, fill = 'lightblue3') +
               geom_hline(yintercept = 0, linetype = 3, linewidth = 0.25) +
               geom_vline(xintercept = c(ymd(c('2023-04-01', '2023-11-01'))), 
                          linetype = 2, linewidth = 0.3) +
               scale_x_date(date_labels = '%b', date_breaks = '2 month') +
               labs(y = 'z-scores', x = NULL) +
               facet_wrap(~class) +
               #lims(y = c(-0.8, 1)) +
               theme_classic() + 
               #scale_color_manual(values = c("tomato", "steelblue", "goldenrod")) +
               #scale_fill_manual(values = c("tomato", "steelblue", "goldenrod")) +
               theme(panel.grid = element_blank(),
                     legend.title = element_blank(),
                     legend.position = 'top',
                     strip.background = element_blank(),
                     axis.line = element_line(linewidth = 0.25),
                     axis.ticks = element_line(linewidth = 0.25),
                     #text = element_text(family = 'Times New Roman'),
                     axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
                     axis.text = element_text(size = 8.5),
                     #title = element_blank(),
                     legend.margin = margin(c(0, 0, 0, 0)))
             p
           })



plot_global_hypervolume <- 
  plot_grid(plot_global_hypervolume[[1]] +
            lims(y = c(-0.75, 0.95)), 
          plot_global_hypervolume[[2]] +
            theme(axis.title.y = element_blank()), 
          plot_global_hypervolume[[3]] +
            lims(y = c(-0.15, 0.19)) +
            theme(axis.title.y = element_blank()), 
          plot_global_hypervolume[[4]] +
            lims(y = c(-0.5, 0.48)) +
            theme(axis.title.y = element_blank()), 
          nrow = 1, 
          rel_widths = c(1, rep(0.95, 3)), 
          labels = c('i', 'ii', 'iii', 'iv'), 
          #label_fontfamily = 'Times New Roman', 
          label_fontface = 'plain', 
          label_size = 12, 
          label_y = 0.825, 
          label_x = c(0.3, 0.25, 0.23, 0.25))

global_hypervolumes <- split(global_hypervolumes, global_hypervolumes$class)


plot_global_hypervolume

max_nut_HV <- 
  global_hypervolumes$`Fruits' nutr. hypervolume` |> 
  filter(month1 == '2023-12-01') |> 
  mutate(est = mean(d$HV_plants_nut) + est * sd(d$HV_plants_nut))

min_nut_HV <- 
  global_hypervolumes$`Fruits' nutr. hypervolume` |> 
  filter(month1 == '2023-04-01') |> 
  mutate(est = mean(d$HV_plants_nut) + est * sd(d$HV_plants_nut))

cont_nut_MAX_MIN <- 100 - (min_nut_HV$est * 100) / max_nut_HV$est

plot(density(cont_nut_MAX_MIN))

mean(cont_nut_MAX_MIN) # mean contrast 
quantile(cont_nut_MAX_MIN, 0.025) # lower boundary CI
quantile(cont_nut_MAX_MIN, 0.975) # lower boundary CI


max_mor_HV <- 
  global_hypervolumes$`Fruits' morph. hypervolume` |> 
  filter(month1 >= '2023-01-01' & month1 <= '2023-03-01' |
           month1 == '2023-12-01') |> 
  mutate(est = mean(sqrt(d$HV_plants_morfo)) + 
           est * sd(sqrt(d$HV_plants_morfo)))

min_mor_HV <- 
  global_hypervolumes$`Fruits' morph. hypervolume` |> 
  filter(month1 >= '2023-04-01' & month1 <= '2023-11-01') |> 
  mutate(est = mean(sqrt(d$HV_plants_morfo)) + 
           est * sd(sqrt(d$HV_plants_morfo)))

cont_mor_MAX_MIN <- 100 - (min_mor_HV$est * 100) / max_mor_HV$est

mean(cont_mor_MAX_MIN) # mean contrast 
quantile(cont_mor_MAX_MIN, 0.025) # lower boundary CI
quantile(cont_mor_MAX_MIN, 0.975) # upper boundary CI
mean(cont_mor_MAX_MIN > 0) # P(higher > lower)



max_PLAN_HV <- 
  global_hypervolumes$`Fruits' total hypervolume` |> 
  filter(month1 >= '2023-02-01' & month1 <= '2023-03-01') |> 
  mutate(est = mean(d$HV_plants_morfo) + 
           est * sd(d$HV_plants_morfo))

min_PLAN_HV <- 
  global_hypervolumes$`Fruits' total hypervolume` |> 
  filter(month1 >= '2023-06-01' & month1 <= '2023-07-01') |> 
  mutate(est = mean(d$HV_plants_morfo) + 
           est * sd(d$HV_plants_morfo))

cont_PLAN_MAX_MIN <- 100 - (min_PLAN_HV$est * 100) / max_PLAN_HV$est

mean(cont_PLAN_MAX_MIN) # mean contrast
quantile(cont_PLAN_MAX_MIN, 0.025) # lower boundary CI
quantile(cont_PLAN_MAX_MIN, 0.975) # upper boundary CI
mean(cont_PLAN_MAX_MIN > 0) # P(higher > lower)


max_BIRD_HV <- 
  global_hypervolumes$`Birds' hypervolume` |> 
  filter(month1 == '2023-04-01') |> 
  mutate(est = mean(d$HV_plants_morfo) + 
           est * sd(d$HV_plants_morfo))

min_BIRD_HV <- 
  global_hypervolumes$`Birds' hypervolume` |> 
  filter(month1 >= '2023-08-01' & month1 <= '2023-09-01') |> 
  mutate(est = mean(d$HV_plants_morfo) + 
           est * sd(d$HV_plants_morfo))

cont_BIRD_MAX_MIN <- 100 - (min_BIRD_HV$est * 100) / max_BIRD_HV$est

mean(cont_BIRD_MAX_MIN) # mean contrast 
quantile(cont_BIRD_MAX_MIN, 0.025) # lower boundary CI
quantile(cont_BIRD_MAX_MIN, 0.975) # upper boundary CI
mean(cont_BIRD_MAX_MIN > 0) # P(higher > lower)

colnames(estimated_nestednessG_effect)[4] <- 'month1'

global_time_NM <- 
  rbind(estimated_MODU_G_effect,
        estimated_nestednessG_effect, 
        estimated_H2_G_effect)

global_time_NM$class <- as.factor(global_time_NM$class)
levels(global_time_NM$class)
global_time_NM$class <- factor(global_time_NM$class, 
                               labels = c("H2'", 'Nestedness', 'Modularity'))


#==== plot network metrics ==== ####

plot_global_NM <- 
  lapply(levels(global_time_NM$class), FUN = 
         function(x) {
           
           p <- 
           global_time_NM |> 
             filter(class == x) |> 
             group_by(class, month1) |> 
             transmute(mu = mean(est), 
                       li = quantile(est, 0.1), 
                       ls = quantile(est, 0.9)) |> 
             unique() |> 
             ggplot(aes(month1, mu, ymin = li, ymax = ls)) +
             geom_line(aes(), linewidth = 0.7, alpha = 0.5, color = 'lightblue3') +
             geom_ribbon(aes(), alpha = 0.3, fill = 'lightblue3') +
             geom_hline(yintercept = 0, linetype = 3, linewidth = 0.25) +
             geom_vline(xintercept = c(ymd(c('2023-04-01', '2023-11-01'))), 
                        linetype = 2, linewidth = 0.3) +
             scale_x_date(date_labels = '%b', date_breaks = '2 month') +
             labs(y = 'z-scores', x = NULL) +
             facet_wrap(~class) +
             theme_classic() + 
             #scale_color_manual(values = c("tomato", "steelblue", "goldenrod")) +
             #scale_fill_manual(values = c("tomato", "steelblue", "goldenrod")) +
             theme(panel.grid = element_blank(),
                   legend.title = element_blank(),
                   legend.position = 'top',
                   strip.background = element_blank(),
                   axis.line = element_line(linewidth = 0.25),
                   axis.ticks = element_line(linewidth = 0.25),
                   #text = element_text(family = 'Times New Roman'),
                   #axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
                   axis.text = element_text(size = 8.5),
                   #title = element_blank(),
                   legend.margin = margin(c(0, 0, 0, 0)))
             p
         })

plot_global_network_metrics <- 
  plot_grid(plot_global_NM[[1]], 
            plot_global_NM[[2]] +
              theme(axis.title.y = element_blank()),
            plot_global_NM[[3]] + 
              theme(axis.title.y = element_blank()) +
              lims(y = c(-0.65, 0.8)), 
            nrow = 1, 
            rel_widths = c(1, 0.95, 0.95), 
            labels = c('i', 'ii', 'ii'), 
            #label_fontfamily = 'Times New Roman', 
            label_fontface = 'plain', 
            label_size = 12, 
            label_y = 0.825, 
            label_x = c(0.21, 0.185, 0.185))

global_time_NM <- split(global_time_NM, global_time_NM$class)


plot_global_network_metrics

max_NSS <- 
  global_time_NM$Nestedness |>
  filter(month == 12) |> # December
  mutate(est = mean(d$nestedness) + 
           est * sd(d$nestedness))

min_NSS <- 
  global_time_NM$Nestedness |>
  filter(month >= 8 & month <= 9) |> # August and September
  mutate(est = mean(d$nestedness) + 
           est * sd(d$nestedness))

max_NSS <- 
  global_time_NM$Nestedness |>
  filter(month == 12) |> 
  mutate(est = mean(d$nestedness) + 
           est * sd(d$nestedness))

min_NSS <- 
  global_time_NM$Nestedness |>
  filter(month >= 8 & month <= 9) |> 
  mutate(est = mean(d$nestedness) + 
           est * sd(d$nestedness))

cont_NSS_MIN_MAX <- 
  100 - ((min_NSS$est * 100) / max_NSS$est)

mean(cont_NSS_MIN_MAX) # average difference
quantile(cont_NSS_MIN_MAX, 0.025) # low boundary of CI
quantile(cont_NSS_MIN_MAX, 0.975) # upper boundary of CIU
mean(cont_NSS_MIN_MAX > 0) # P(higher month > lower month)


max_MODU <- 
  global_time_NM$Modularity |>
  filter(month >= 10 & month <= 11) |> 
  mutate(est = mean(d$modularity) + 
           est * sd(d$modularity))

min_MODU <- 
  global_time_NM$Modularity |>
  filter(month >= 4 & month <= 5) |> 
  mutate(est = mean(d$modularity) + 
           est * sd(d$modularity))

unique(min_MODU$month)

cont_MODU_MIN_MAX <- max_MODU$est - min_MODU$est
#100 - ((min_MODU$est * 100) / max_MODU$est)


mean(cont_MODU_MIN_MAX) * 100 # average difference
quantile(cont_MODU_MIN_MAX, 0.025) * 100 # low boundary of CI
quantile(cont_MODU_MIN_MAX, 0.975) * 100 # upper boundary of CIU
mean(cont_MODU_MIN_MAX > 0) # P(higher month > lower month)



max_h2 <- 
  global_time_NM$`H2'` |>
  filter(month == 4 | month == 9) |> 
  mutate(est = mean(d$H2) + 
           est * sd(d$H2))

min_h2 <- 
  global_time_NM$`H2'` |>
  filter(month == 12 | month == 7) |> 
  mutate(est = mean(d$H2) + 
           est * sd(d$H2))

cont_h2_MIN_MAX <- max_h2$est - min_h2$est
  #100 - ((min_h2$est * 100) / max_h2$est)

mean(cont_h2_MIN_MAX) * 100 # average difference
quantile(cont_h2_MIN_MAX, 0.025) * 100 # low boundary of CI
quantile(cont_h2_MIN_MAX, 0.975) * 100 # upper boundary of CIU
mean(cont_h2_MIN_MAX > 0) # P(higher month > lower month)


plot_grid(plot_grid(NULL, 
                    plot_global_climate +
                      theme(legend.background = element_blank()), 
                    plot_global_phenology +
                      theme(axis.title.y = element_blank(), 
                            legend.background = element_blank()),
                    NULL, nrow = 1, 
                    rel_widths = c(0.2, 1, 0.975, 0.2), 
                    labels = c('', 'i', 'ii', ''), 
                    #label_fontfamily = 'Times New Roman', 
                    label_fontface = 'plain', 
                    label_size = 12, 
                    label_y = 0.85, 
                    label_x = c(0, 0.15, 0.13, 0)), 
          plot_global_hypervolume, 
          plot_global_network_metrics, 
          ncol = 1, 
          labels = paste0('(', letters, ')'), 
          #label_fontfamily = 'Times New Roman', 
          label_fontface = 'plain', 
          label_size = 12, 
          label_x = c(0, 0, 0))

betas_MODU_G$p <- 
  sapply(1:nrow(betas_MODU_G), FUN = 
         function(x) {
           pos <- betas_MODU_G$p_pos[x]
           neg <- betas_MODU_G$p_neg[x]
           
           if (pos > neg) pos
           else neg
           
         })

betas_MODU_G$factor <- gsub('^(beta_)([a-zA-Z]*)(_)(.*)$', '\\2', 
                            betas_MODU_G$var)

betas_MODU_G$response <- gsub('^(beta_)([a-zA-Z]*)(_)(.*)$', '\\4', 
                              betas_MODU_G$var)

betas_H2_G$p <- 
  sapply(1:nrow(betas_H2_G), FUN = 
           function(x) {
             pos <- betas_H2_G$p_pos[x]
             neg <- betas_H2_G$p_neg[x]
             
             if (pos > neg) pos
             else neg
             
           })

betas_H2_G$factor <- gsub('^(beta_)([a-zA-Z]*)(_)(.*)$', '\\2', 
                          betas_H2_G$var)

betas_H2_G$response <- gsub('^(beta_)([a-zA-Z]*)(_)(.*)$', '\\4', 
                            betas_H2_G$var)

betas_nestednessG$p <- 
  sapply(1:nrow(betas_nestednessG), FUN = 
           function(x) {
             pos <- betas_nestednessG$p_pos[x]
             neg <- betas_nestednessG$p_neg[x]
             
             if (pos > neg) pos
             else neg
             
           })

betas_nestednessG$factor <- gsub('^(beta_)([a-zA-Z]*)(_)(.*)$', '\\2', 
                                 betas_nestednessG$var)

betas_nestednessG$response <- gsub('^(beta_)([a-zA-Z]*)(_)(.*)$', '\\4', 
                            betas_nestednessG$var)


effects_net <- 
  rbind(betas_H2_G[betas_H2_G$p >= 0.7,], 
      betas_MODU_G[grep('MODU$', betas_MODU_G$var),],
      betas_nestednessG[grep('NSS$', betas_nestednessG$var),])

effects_net <- effects_net[, c("factor", "response", 'p')]

colnames(effects_net) <- c('from', 'to', 'weight')

effects_net$from[grep('^nut$', effects_net$from)] <- 'nut_HV'
effects_net$from[grep('^mor$', effects_net$from)] <- 'mor_HV'
effects_net$from[grep('^BIRD$', effects_net$from)] <- 'BIRD_HV'
effects_net$from[grep('^PLAN$', effects_net$from)] <- 'PLAN_HV'

vertex_df <- 
  tibble(var = unique(c(effects_net$from, effects_net$to)), 
       type = c('Climate', 'Climate', 'Phenology', 'Phenology', 
                'Plants', 'Plants', 'Birds', 'Plants', 'Network', 
                'Network', 'Network'), 
       color = c('lightblue3', 'lightblue3', 'seagreen', 'seagreen', 
                 'tomato3', 'tomato3', 'tan1', 'tomato3', 'cyan4', 
                 'cyan4', 'cyan4'))

#effects_net <- effects_net[effects_net$weight >= 0.8, ]

net_plot <- 
  graph_from_data_frame(effects_net, directed = T, vertices = vertex_df)


sp <- lapply(1:2, FUN = 
               function(x) {
                 shortest_paths(net_plot, from = x, to = 11, 
                                weights = E(net_plot)$weight)
               })

(shortest_path <- lapply(1:2, FUN = 
                           function(x) {
                             sp[[x]]$vpath[[1]]
                           }))

all_paths_nestednes <- lapply(1:2, FUN = 
                                function(x) {
                                  all_simple_paths(net_plot, from = x, to = 11)
                                  })

path_lengths_nestednes <- 
  lapply(all_paths_nestednes, FUN = 
           function(i) {
             sapply(i, function(path) {
               sum(E(net_plot, path = path)$weight)
             })
           })
all_paths_nestednes

nestedness_paths <- 
  lapply(all_paths_nestednes, FUN = 
         function(z) {
           d <- 
             lapply(z, FUN = 
                      function(x) {
                        
                        x <- as.numeric(x)
                        
                        edges <- c('Temperature', 'Rainfall', 'Fr. Production', 'Fr. S.', 
                                   'Nutr. ifs', 'Morph. ifs', 'Birds. ifs', 
                                   'Fr. Tot. ifs', "H2'", 'Modularity', 'Nestedness')
                        
                        x <- edges[x]
                        
                        v <- matrix(NA, nrow = 1, ncol = 5)
                        
                        for (i in seq_along(x)) {
                          v[1, i] <- x[i]
                        }
                        
                        v
                      })
           d <- as_tibble(do.call('rbind', d))
           colnames(d) <- paste0('Stage ', 1:5)
           d
         })

nestedness_paths <- do.call('rbind', nestedness_paths)

#nestedness_paths <- 
  nestedness_paths |> 
  make_long(`Stage 1`, `Stage 2`, `Stage 3`, `Stage 4`, `Stage 5`) |> 
  ggplot(aes(
  x = x,
  next_x = next_x,
  node = node,
  next_node = next_node,
  fill = node,
  label = node
)) +
  geom_sankey(flow.alpha = 0.6, node.color = "gray30") +
  geom_sankey_label(size = 3, color = "black", fill = "white") +
  theme_sankey(base_size = 14) +
  #ggtitle("Nestedness") +
  #scale_fill_brewer(palette = "Set2") +
  theme(axis.title = element_blank(), 
        legend.position = 'none')


sp <- lapply(1:2, FUN = 
               function(x) {
                 shortest_paths(net_plot, from = x, to = 10, 
                                weights = E(net_plot)$weight)
               })

(shortest_path_MODU <- lapply(1:2, FUN = 
                           function(x) {
                             sp[[x]]$vpath[[1]]
                           }))

all_paths_MODU <- lapply(1:2, FUN = 
                                function(x) {
                                  all_simple_paths(net_plot, from = x, to = 10)
                                })

all_paths_MODU


MODU_paths <- 
  lapply(all_paths_MODU, FUN = 
           function(z) {
             d <- 
               lapply(z, FUN = 
                        function(x) {
                          
                          x <- as.numeric(x)
                          
                          edges <- c('Temperature', 'Rainfall', 'Fr. Production', 'Fr. S.', 
                                     'Nutr. ifs', 'Morph. ifs', 'Birds. ifs', 
                                     'Fr. Tot. ifs', "H2'", 'Modularity', 'Nestedness')
                          
                          x <- edges[x]
                          
                          v <- matrix(NA, nrow = 1, ncol = 5)
                          
                          for (i in seq_along(x)) {
                            v[1, i] <- x[i]
                          }
                          
                          v
                        })
             d <- as_tibble(do.call('rbind', d))
             colnames(d) <- paste0('Stage ', 1:5)
             d
           })

MODU_paths <- do.call('rbind', MODU_paths)



#MODU_paths <- 
  MODU_paths |> 
  make_long(`Stage 1`, `Stage 2`, `Stage 3`, `Stage 4`, `Stage 5`) |> 
    ggplot(aes(
  x = x,
  next_x = next_x,
  node = node,
  next_node = next_node,
  fill = node,
  label = node
)) +
  geom_sankey(flow.alpha = 0.6, node.color = "gray30") +
  geom_sankey_label(size = 3, color = "black", fill = "white") +
  theme_sankey(base_size = 14) +
  #ggtitle("Modularity") +
  #scale_fill_brewer(palette = "Set2") +
  theme(axis.title = element_blank(), 
        legend.position = 'none')




sp <- lapply(1:2, FUN = 
               function(x) {
                 shortest_paths(net_plot, from = x, to = 9, 
                                weights = E(net_plot)$weight)
               })

(shortest_path_H2 <- lapply(1:2, FUN = 
                                function(x) {
                                  sp[[x]]$vpath[[1]]
                                }))
# 

all_paths_H2 <- lapply(1:2, FUN = 
                           function(x) {
                             all_simple_paths(net_plot, from = x, to = 9)
                           })

all_paths_H2


H2_paths <- 
  lapply(all_paths_H2, FUN = 
           function(z) {
             d <- 
               lapply(z, FUN = 
                        function(x) {
                          
                          x <- as.numeric(x)
                          
                          edges <- c('Temperature', 'Rainfall', 'Fr. Production', 'Fr. S.', 
                                     'Nutr. ifs', 'Morph. ifs', 'Birds. ifs', 
                                     'Fr. Tot. ifs', "H2'", 'Modularity', 'Nestedness')
                          
                          x <- edges[x]
                          
                          v <- matrix(NA, nrow = 1, ncol = 5)
                          
                          for (i in seq_along(x)) {
                            v[1, i] <- x[i]
                          }
                          
                          v
                        })
             d <- as_tibble(do.call('rbind', d))
             colnames(d) <- paste0('Stage ', 1:5)
             d
           })

H2_paths <- do.call('rbind', H2_paths)

total_paths <- rbind(nestedness_paths, MODU_paths, H2_paths)

total_paths <- 
  total_paths |> 
  make_long(`Stage 1`, `Stage 2`, `Stage 3`, `Stage 4`, `Stage 5`)

# total_paths$alpha <- 0
# total_paths$alpha[which(sum(is.na(total_paths$node)) > 0)] <- 0.2

H2_paths <- 
  H2_paths |> 
  make_long(`Stage 1`, `Stage 2`)

#total_paths_plot <- 
  ggplot(H2_paths, aes(
    x = x,
    next_x = next_x,
    node = node,
    next_node = next_node,
    fill = node,
    label = node
  )) +
  geom_sankey(flow.alpha = 0.6, node.color = "gray30") +
  geom_sankey_label(size = 3, color = "black", fill = "white",) +
  theme_sankey(base_size = 18) +
  #ggtitle("H2") +
  #scale_fill_brewer(palette = "Set2") +
  theme(axis.title = element_blank(), 
        legend.position = 'none')



 ggplot(total_paths, aes(
    x = x,
    next_x = next_x,
    node = node,
    next_node = next_node,
    fill = node,
    label = node
  )) +
  geom_sankey(flow.alpha = 0.6, node.color = "gray30") +
  geom_sankey_label(size = 3, color = "black", fill = "white",) +
  theme_sankey(base_size = 18) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(11, "Spectral")) +
  #ggtitle("H2") +
  #scale_fill_brewer(palette = "Set2") +
  theme(axis.title = element_blank(), 
        legend.position = 'none')

causal_effects <- 
  function(pars = post_H2_G, # object with posterior distributions of the BN model
           alpha_1, f_1, # site and GP parameters for equation 1
           alpha_2, f_2, # site and GP parameters for equation 2
           alpha_3, f_3, # site and GP parameters for equation 3
           alpha_4, f_4, # site and GP parameters for equation 4
           beta_1, beta_2, beta_3, beta_4, # effects to assess at equation i
           var, # exposure or predictor variable to assess
           y_obs, # outcome variable in its natural scale
           level_1 = T, # direct effect (x -> y)
           level_2 = F, # indirect effect (x -> z -> y)
           level_3 = F, # indirect effect (x -> z -> k -> y)
           level_4 = F  # indirect effect (x -> z -> k -> l -> y)
           ) {
    
    if (level_1) {
      
      # minimum and maximum values of predictor variable 
      x_1 <- quantile(var, c(0, 1))
      
      # indexing site and GP parameters EQ 1
      a <- grep(alpha_1, colnames(pars$alpha))
      f_ <- grep(f_1, colnames(pars$f))
      
      effect <- # simulating the intervention
        sapply(x_1, FUN = 
               function(x) {
                 
                 with(pars, {
                   # EQ 1
                   apply(alpha[, a], 1, mean) +
                     apply(f[, f_], 1, mean) +
                     beta[[beta_1]] * x
                 })
               })
      
      colnames(effect) <- c('min', 'max')
      
      # going from z-scores to the scale of the variable
      if ((grepl('MODU', beta_1) +
           grepl('H2', beta_1)) == 0) {
        effect <- apply(effect, 2, function(z) mean(y_obs) + z * sd(y_obs))
      } else {
        effect <- apply(effect, 2, function(z) inv_logit(z))
      }
      
      effect <- as_tibble(effect)
      #contrast
      effect$contrast <- effect %$% (max - min)
      # p(effect)
      effect$p_effect <- mean(effect$contrast > 0)
      # % change
      effect$relative_change <- ((effect$contrast * 100) / mean(y_obs))
      effect$effect <- paste0(beta_1)
    } 
    
    if (level_2) {
      # minimum and maximum values of predictor variable 
      x_1 <- quantile(var, c(0, 1))
      
      # indexing site and GP parameters EQ 1
      a_1 <- grep(alpha_1, colnames(pars$alpha))
      f_1 <- grep(f_1, colnames(pars$f))
      # indexing site and GP parameters EQ 2
      a_2 <- grep(alpha_2, colnames(pars$alpha))
      f_2 <- grep(f_2, colnames(pars$f))
      
      effect <- # simulating the intervention
        sapply(x_1, FUN = 
                 function(x) {
                   with(pars, {
                     # EQ 1
                     effect1 <- 
                       apply(alpha[, a_1], 1, mean) +
                       apply(f[, f_1], 1, mean) +
                       beta[[beta_1]] * x
                     
                     # EQ 2
                     effect2 <- 
                       apply(alpha[, a_2], 1, mean) +
                       apply(f[, f_2], 1, mean) +
                       beta[[beta_2]] * effect1
                     
                     effect2
                   })
                 })
      
      colnames(effect) <- c('min', 'max')
      
      # going from z-scores to the scale of the variable
      if ((grepl('MODU', beta_2) +
           grepl('H2', beta_2)) == 0) {
        effect <- apply(effect, 2, function(z) mean(y_obs) + z * sd(y_obs))
      } else {
        effect <- apply(effect, 2, function(z) inv_logit(z))
      }
      effect <- as_tibble(effect)
      #contrast
      effect$contrast <- effect %$% (max - min)
      # p(effect)
      effect$p_effect <- mean(effect$contrast > 0)
      # % change
      effect$relative_change <- ((effect$contrast * 100) / mean(y_obs))
      effect$effect <- paste0(beta_1, '->', beta_2)
    }
    
    if (level_3) {
      # minimum and maximum values of predictor variable 
      x_1 <- quantile(var, c(0, 1))
      
      # indexing site and GP parameters EQ 1
      a_1 <- grep(alpha_1, colnames(pars$alpha))
      f_1 <- grep(f_1, colnames(pars$f))
      # indexing site and GP parameters EQ 2
      a_2 <- grep(alpha_2, colnames(pars$alpha))
      f_2 <- grep(f_2, colnames(pars$f))
      # indexing site and GP parameters EQ 3
      a_3 <- grep(alpha_3, colnames(pars$alpha))
      f_3 <- grep(f_3, colnames(pars$f))
      
      effect <- # simulating the intervention
        sapply(x_1, FUN = 
                 function(x) {
                   with(pars, {
                     # EQ 1
                     effect1 <- 
                       apply(alpha[, a_1], 1, mean) +
                       apply(f[, f_1], 1, mean) +
                       beta[[beta_1]] * x
                     
                     # EQ 2
                     effect2 <- 
                       apply(alpha[, a_2], 1, mean) +
                       apply(f[, f_2], 1, mean) +
                       beta[[beta_2]] * effect1
                     
                     # EQ 3
                     effect3 <- 
                       apply(alpha[, a_3], 1, mean) +
                       apply(f[, f_3], 1, mean) +
                       beta[[beta_3]] * effect2
                     
                     effect3
                   })
                 })
      
      colnames(effect) <- c('min', 'max')
      
      # going from z-scores to the scale of the variable
      if ((grepl('MODU', beta_3) +
           grepl('H2', beta_3)) == 0) {
        effect <- apply(effect, 2, function(z) mean(y_obs) + z * sd(y_obs))
      } else {
        effect <- apply(effect, 2, function(z) inv_logit(z))
      }
      effect <- as_tibble(effect)
      #contrast
      effect$contrast <- effect %$% (max - min)
      # p(effect)
      effect$p_effect <- mean(effect$contrast > 0)
      # % change
      effect$relative_change <- ((effect$contrast * 100) / mean(y_obs))
      effect$effect <- paste0(beta_1, '->', beta_2, '->', beta_3)
      
    }
    
    if (level_4) {
      # minimum and maximum values of predictor variable 
      x_1 <- quantile(var, c(0, 1))
      
      # indexing site and GP parameters EQ 1
      a_1 <- grep(alpha_1, colnames(pars$alpha))
      f_1 <- grep(f_1, colnames(pars$f))
      # indexing site and GP parameters EQ 2
      a_2 <- grep(alpha_2, colnames(pars$alpha))
      f_2 <- grep(f_2, colnames(pars$f))
      # indexing site and GP parameters EQ 3
      a_3 <- grep(alpha_3, colnames(pars$alpha))
      f_3 <- grep(f_3, colnames(pars$f))
      # indexing site and GP parameters EQ 4
      a_4 <- grep(alpha_4, colnames(pars$alpha))
      f_4 <- grep(f_4, colnames(pars$f))
      
      effect <- # simulating the intervention
        sapply(x_1, FUN = 
                 function(x) {
                   with(pars, {
                     # EQ 1
                     effect1 <- 
                       apply(alpha[, a_1], 1, mean) +
                       apply(f[, f_1], 1, mean) +
                       beta[[beta_1]] * x
                     
                     # EQ 2
                     effect2 <- 
                       apply(alpha[, a_2], 1, mean) +
                       apply(f[, f_2], 1, mean) +
                       beta[[beta_2]] * effect1
                     
                     # EQ 3
                     effect3 <- 
                       apply(alpha[, a_3], 1, mean) +
                       apply(f[, f_3], 1, mean) +
                       beta[[beta_3]] * effect2
                     
                     # EQ 4
                     effect4 <- 
                       apply(alpha[, a_4], 1, mean) +
                       apply(f[, f_4], 1, mean) +
                       beta[[beta_4]] * effect2
                     
                     effect4
                   })
                 })
      
      colnames(effect) <- c('min', 'max')
      
      # going from z-scores to the scale of the variable
      if ((grepl('MODU', beta_4) +
           grepl('H2', beta_4)) == 0) {
        effect <- apply(effect, 2, function(z) mean(y_obs) + z * sd(y_obs))
      } else {
        effect <- apply(effect, 2, function(z) inv_logit(z))
      }
      effect <- as_tibble(effect)
      #contrast
      effect$contrast <- effect %$% (max - min)
      # p(effect)
      effect$p_effect <- mean(effect$contrast > 0)
      # % change
      effect$relative_change <- ((effect$contrast * 100) / mean(y_obs))
      effect$effect <- paste0(beta_1, '->', beta_2, '->', beta_3, '->', beta_4)
    }
    
    effect
  }

knitr::kable(betas_H2_G_z[grep('BIRD_HV$', betas_H2_G_z[[1]]), ], 
             digits = 2, 
             caption = "Factors affecting birds' functional space. Average values, 95% CIs, and probabilities of directional effects for the relationships assessed in the generative model. An effect is marked as TRUE if P(effect > 0) > 0.7 or P(effect < 0) > 0.7 or.")

# fruit production - temperature 

((betas_H2_G[grep('BIRD_HV$', betas_H2_G$var), ]$mu[3] - 
  abs(betas_H2_G[grep('BIRD_HV$', betas_H2_G$var), ]$mu[1])) /
  abs(betas_H2_G[grep('BIRD_HV$', betas_H2_G$var), ]$mu[1])) * 100

# fruit production - rainfall 
((abs(betas_H2_G[grep('BIRD_HV$', betas_H2_G$var), ]$mu[3]) - 
    abs(betas_H2_G[grep('BIRD_HV$', betas_H2_G$var), ]$mu[2])) /
    abs(betas_H2_G[grep('BIRD_HV$', betas_H2_G$var), ]$mu[2])) * 100


effect_fruitAB_bird <- 
  causal_effects(
  pars = post_H2_G, 
  alpha_1 = 'alpha_BIRD_HV',
  f_1 = 'f_BIRD_HV', 
  beta_1 = 'beta_fruitAB_BIRD_HV', 
  var = dat$fruit_abun,
  y_obs = d$HV_bird,
  level_1 = T
)


quantile(effect_fruitAB_bird$relative_change, c(0.025, 0.5, 0.975))

unique(effect_fruitAB_bird$p_effect)

effect_rain_bird <- 
  causal_effects(
    pars = post_H2_G, 
    alpha_1 = 'alpha_BIRD_HV',
    f_1 = 'f_BIRD_HV', 
    beta_1 = 'beta_rain_BIRD_HV', 
    var = dat$z_rainfall,
    y_obs = d$HV_bird,
    level_1 = T
  )


quantile(effect_rain_bird$relative_change, c(0.025, 0.5, 0.975))

unique(effect_rain_bird$p_effect)

effect_temp_bird <- 
  causal_effects(
    pars = post_H2_G, 
    alpha_1 = 'alpha_BIRD_HV',
    f_1 = 'f_BIRD_HV', 
    beta_1 = 'beta_temp_BIRD_HV', 
    var = dat$z_rainfall,
    y_obs = d$HV_bird,
    level_1 = T
  )


quantile(effect_temp_bird$relative_change, c(0.025, 0.5, 0.975))

unique(effect_temp_bird$p_effect)

effect_rain_fruitAB_bird <- 
  causal_effects(
    pars = post_H2_G, 
    alpha_1 = 'alpha_fruitAB',
    f_1 = 'f_fruitAB', 
    beta_1 = 'beta_rain_fruitAB',
    alpha_2 = 'alpha_BIRD_HV',
    f_2 = 'f_BIRD_HV', 
    beta_2 = 'beta_fruitAB_BIRD_HV', 
    var = dat$z_rainfall,
    y_obs = d$HV_bird,
    level_1 = F, 
    level_2 = T
  )


quantile(effect_rain_fruitAB_bird$relative_change, c(0.025, 0.5, 0.975))

unique(effect_rain_fruitAB_bird$p_effect)

df_causal_BIRD <- 
  rbind(effect_fruitAB_bird,
        effect_temp_bird,
        effect_rain_fruitAB_bird) |> 
  group_by(effect) |> 
  transmute(mu = median(relative_change), 
            li = quantile(relative_change, 0.025), 
            ls = quantile(relative_change, 0.975)) |> 
  unique() |> ungroup() |> 
  mutate(x = 
           c('Fruit production', 
             'Temperature', 
             '*Rainfall      \n (fruit production)'), 
         y = "Birds' hypervolume") 

causal_plot_BIRD_HV <- 
  rbind(effect_fruitAB_bird,
      effect_temp_bird,
      effect_rain_fruitAB_bird) |> 
  group_by(effect) |> 
  transmute(mu = median(relative_change), 
            li = quantile(relative_change, 0.025), 
            ls = quantile(relative_change, 0.975)) |> 
  unique() |> ungroup() |> 
  mutate(x = 
           c('Fruit production', 
             'Temperature', 
             '*Rainfall      \n (fruit production)'), 
         y = "Birds' Hypervolume") |> 
  ggplot() +
  geom_hline(yintercept = 0, linetype = 3, linewidth = 0.25) +
  geom_errorbar(aes(x, ymin = li, ymax = ls), 
                width = 0, linewidth = 1.25, alpha = 0.5, 
                color = '#1e81b0') +
  geom_point(aes(x, mu), color = '#e28743', size = 2) +
  labs(y = 'Relative causal effect (%)', x = '') +
  facet_wrap(~y) +
  theme_classic() + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.line = element_line(linewidth = 0.25),
        axis.ticks = element_line(linewidth = 0.25),
        axis.title.x = element_blank(),
        #text = element_text(family = 'Times New Roman'),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text = element_text(size = 8.5))


causal_plot_BIRD_HV

knitr::kable(betas_H2_G_z[grep('nut_HV$', betas_H2_G_z[[1]]), ], 
             digits = 2, 
             caption = "Factors affecting fruits' nutritional functional space. Average values, 95% CIs, and probabilities of directional effects for the relationships assessed in the generative model. An effect is marked as TRUE if P(effect > 0) > 0.7 or P(effect < 0) > 0.7 or.")

# temperature - fruitAB
((abs(betas_H2_G[grep('nut_HV$', betas_H2_G$var), ]$mu[1]) - 
    abs(betas_H2_G[grep('nut_HV$', betas_H2_G$var), ]$mu[3])) /
  abs(betas_H2_G[grep('nut_HV$', betas_H2_G$var), ]$mu[3])) * 100



effect_temp_nut <- 
  causal_effects(
    pars = post_H2_G, 
    alpha_1 = 'alpha_nut_HV',
    f_1 = 'f_nut_HV', 
    beta_1 = 'beta_temp_nut_HV', 
    var = dat$z_temperature,
    y_obs = d$HV_plants_nut,
    level_1 = T
  )


quantile(effect_temp_nut$relative_change, 
         c(0.025, 0.5, 0.975))

unique(effect_temp_nut$p_effect)


effect_fruitAB_nut <- 
  causal_effects(
    pars = post_H2_G, 
    alpha_1 = 'alpha_nut_HV',
    f_1 = 'f_nut_HV', 
    beta_1 = 'beta_fruitAB_nut_HV', 
    var = dat$fruit_abun,
    y_obs = d$HV_plants_nut,
    level_1 = T
  )


quantile(effect_fruitAB_nut$relative_change, 
         c(0.025, 0.5, 0.975))

unique(effect_fruitAB_nut$p_effect)

df_causal_nut <- 
  rbind(effect_temp_nut,
        effect_fruitAB_nut) |> 
  group_by(effect) |> 
  transmute(mu = median(relative_change), 
            li = quantile(relative_change, 0.025), 
            ls = quantile(relative_change, 0.975)) |> 
  unique() |> ungroup() |> 
  mutate(x = 
           c('Temperature', 
             'Fruit production'), 
         y = "Fruits' nutr. hypervolume")

causal_plot_nut_HV <- 
  rbind(effect_temp_nut,
        effect_fruitAB_nut) |> 
  group_by(effect) |> 
  transmute(mu = median(relative_change), 
            li = quantile(relative_change, 0.025), 
            ls = quantile(relative_change, 0.975)) |> 
  unique() |> ungroup() |> 
  mutate(x = 
           c('Temperature', 
             'Fruit production'), 
         y = "Fruits' nutr. Hypervolume") |> 
  ggplot() +
  geom_hline(yintercept = 0, linetype = 3, linewidth = 0.25) +
  geom_errorbar(aes(x, ymin = li, ymax = ls), 
                width = 0, linewidth = 1.25, alpha = 0.5, 
                color = '#1e81b0') +
  geom_point(aes(x, mu), color = '#e28743', size = 2) +
  labs(y = 'Relative causal effect (%)', x = '') +
  facet_wrap(~y) +
  theme_classic() + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.line = element_line(linewidth = 0.25),
        axis.ticks = element_line(linewidth = 0.25),
        #text = element_text(family = 'Times New Roman'),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text = element_text(size = 8.5))


causal_plot_nut_HV

knitr::kable(betas_H2_G_z[grep('mor_HV$', betas_H2_G_z[[1]]), ], 
             digits = 2, 
             caption = "Factors affecting fruits' nutritional functional space. Average values, 95% CIs, and probabilities of directional effects for the relationships assessed in the generative model. An effect is marked as TRUE if P(effect > 0) > 0.7 or P(effect < 0) > 0.7 or.")


effect_fruitAB_mor <- 
  causal_effects(
    pars = post_H2_G, 
    alpha_1 = 'alpha_mor_HV',
    f_1 = 'f_mor_HV', 
    beta_1 = 'beta_fruitAB_mor_HV', 
    var = dat$fruit_abun,
    y_obs = d$HV_plants_morfo,
    level_1 = T
  )


quantile(effect_fruitAB_mor$relative_change, 
         c(0.025, 0.5, 0.975))

unique(effect_fruitAB_mor$p_effect)


effect_rain_fruitAB_mor <- 
  causal_effects(
    pars = post_H2_G, 
    alpha_1 = 'alpha_fruitAB',
    f_1 = 'f_fruitAB', 
    beta_1 = 'beta_rain_fruitAB', 
    alpha_2 = 'alpha_mor_HV',
    f_2 = 'f_mor_HV', 
    beta_2 = 'beta_fruitAB_mor_HV', 
    var = dat$z_rainfall,
    y_obs = d$HV_plants_morfo,
    level_1 = F, 
    level_2 = T 
  )


quantile(effect_rain_fruitAB_mor$relative_change, 
         c(0.025, 0.5, 0.975))

unique(effect_rain_fruitAB_mor$p_effect)

df_causal_mor <- 
  rbind(effect_fruitAB_mor,
        effect_rain_fruitAB_mor) |> 
  group_by(effect) |> 
  transmute(mu = median(relative_change), 
            li = quantile(relative_change, 0.025), 
            ls = quantile(relative_change, 0.975)) |> 
  unique() |> ungroup() |> 
  mutate(x = 
           c('Fruit production', 
             '*Rainfall      \n (fruit production)'), 
         y = "Fruits' morph. hypervolume") 

causal_plot_mor_HV <- 
  rbind(effect_fruitAB_mor,
        effect_rain_fruitAB_mor) |> 
  group_by(effect) |> 
  transmute(mu = median(relative_change), 
            li = quantile(relative_change, 0.025), 
            ls = quantile(relative_change, 0.975)) |> 
  unique() |> ungroup() |> 
  mutate(x = 
           c('Fruit production', 
             '*Rainfall      \n (fruit production)'), 
         y = "Fruits' morph. Hypervolume") |> 
  ggplot() +
  geom_hline(yintercept = 0, linetype = 3, linewidth = 0.25) +
  geom_errorbar(aes(x, ymin = li, ymax = ls), 
                width = 0, linewidth = 1.25, alpha = 0.5, 
                color = '#1e81b0') +
  geom_point(aes(x, mu), color = '#e28743', size = 2) +
  labs(y = 'Relative causal effect (%)', x = '') +
  facet_wrap(~y) +
  theme_classic() + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.line = element_line(linewidth = 0.25),
        axis.ticks = element_line(linewidth = 0.25),
        #text = element_text(family = 'Times New Roman'),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text = element_text(size = 8.5))


causal_plot_mor_HV

knitr::kable(betas_H2_G_z[grep('PLAN_HV$', betas_H2_G_z[[1]]), ], 
             digits = 2, 
             caption = "Factors affecting fruits' nutritional functional space. Average values, 95% CIs, and probabilities of directional effects for the relationships assessed in the generative model. An effect is marked as TRUE if P(effect > 0) > 0.7 or P(effect < 0) > 0.7 or.")

zzz <- 
  sapply(c(3, 4, 7), FUN = 
         function(x) {
           ((abs(betas_H2_G_z[grep('PLAN_HV$', betas_H2_G_z[[1]]), ]$Mean[2]) -
              abs(betas_H2_G_z[grep('PLAN_HV$', betas_H2_G_z[[1]]), ]$Mean[x])) /
             abs(betas_H2_G_z[grep('PLAN_HV$', betas_H2_G_z[[1]]), ]$Mean[x])) * 100
         })

mean(zzz); sqrt(sd(zzz)/length(zzz))

tibble(var = betas_H2_G_z[grep('PLAN_HV$', betas_H2_G_z[[1]]), ][[1]][c(3, 4, 7)], 
       `Importance of rainfall` = zzz)


zzz <- 
  sapply(c(2, 4, 7), FUN = 
         function(x) {
           ((abs(betas_H2_G_z[grep('PLAN_HV$', betas_H2_G_z[[1]]), ]$Mean[3]) -
              abs(betas_H2_G_z[grep('PLAN_HV$', betas_H2_G_z[[1]]), ]$Mean[x])) /
             abs(betas_H2_G_z[grep('PLAN_HV$', betas_H2_G_z[[1]]), ]$Mean[x])) * 100
         })

mean(zzz); sqrt(sd(zzz)/length(zzz))

tibble(var = betas_H2_G_z[grep('PLAN_HV$', betas_H2_G_z[[1]]), ][[1]][c(2, 4, 7)], 
       `Importance of fruit production` = zzz)



effect_rain_PLAN <- 
  causal_effects(
    pars = post_H2_G, 
    alpha_1 = 'alpha_PLAN_HV',
    f_1 = 'f_PLAN_HV', 
    beta_1 = 'beta_rain_PLAN_HV', 
    var = dat$z_rainfall,
    y_obs = d$HV_plant,
    level_1 = T
  )


quantile(effect_rain_PLAN$relative_change, 
         c(0.025, 0.5, 0.975))

unique(effect_rain_PLAN$p_effect)


effect_fruitAB_PLAN <- 
  causal_effects(
    pars = post_H2_G, 
    alpha_1 = 'alpha_PLAN_HV',
    f_1 = 'f_PLAN_HV', 
    beta_1 = 'beta_fruitAB_PLAN_HV', 
    var = dat$fruit_abun,
    y_obs = d$HV_plant,
    level_1 = T
  )


quantile(effect_fruitAB_PLAN$relative_change, 
         c(0.025, 0.5, 0.975))

unique(effect_fruitAB_PLAN$p_effect)


effect_rain_fruitAB_PLAN <- 
  causal_effects(
    pars = post_H2_G, 
    alpha_1 = 'alpha_fruitAB',
    f_1 = 'f_fruitAB', 
    beta_1 = 'beta_rain_fruitAB', 
    alpha_2 = 'alpha_PLAN_HV',
    f_2 = 'f_PLAN_HV', 
    beta_2 = 'beta_fruitAB_PLAN_HV', 
    var = dat$z_rainfall,
    y_obs = d$HV_plants_morfo,
    level_1 = F, 
    level_2 = T 
  )


quantile(effect_rain_fruitAB_PLAN$relative_change, 
         c(0.025, 0.5, 0.975))

unique(effect_rain_fruitAB_PLAN$p_effect)

df_causal_PLAN <- 
  rbind(effect_rain_PLAN, 
        effect_fruitAB_PLAN,
        effect_rain_fruitAB_PLAN) |> 
  group_by(effect) |> 
  transmute(mu = median(relative_change), 
            li = quantile(relative_change, 0.025), 
            ls = quantile(relative_change, 0.975)) |> 
  unique() |> ungroup() |> 
  mutate(x = 
           c('Rainfall', 
             'Fruit production', 
             '*Rainfall      \n (fruit production)'), 
         y = "Fruits' total hypervolume")

causal_plot_PLAN_HV <- 
  rbind(effect_rain_PLAN, 
        effect_fruitAB_PLAN,
        effect_rain_fruitAB_PLAN) |> 
  group_by(effect) |> 
  transmute(mu = median(relative_change), 
            li = quantile(relative_change, 0.025), 
            ls = quantile(relative_change, 0.975)) |> 
  unique() |> ungroup() |> 
  mutate(x = 
           c('Rainfall', 
             'Fruit production', 
             '*Rainfall      \n (fruit production)'), 
         y = "Fruits' total Hypervolume") |> 
  ggplot() +
  geom_hline(yintercept = 0, linetype = 3, linewidth = 0.25) +
  geom_errorbar(aes(x, ymin = li, ymax = ls), 
                width = 0, linewidth = 1.25, alpha = 0.5, 
                color = '#1e81b0') +
  geom_point(aes(x, mu), color = '#e28743', size = 2) +
  labs(y = 'Relative causal effect (%)', x = '') +
  facet_wrap(~y) +
  theme_classic() + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.line = element_line(linewidth = 0.25),
        axis.ticks = element_line(linewidth = 0.25),
        # text = element_text(family = 'Times New Roman'),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text = element_text(size = 8.5))

causal_plot_PLAN_HV

knitr::kable(betas_nestednessG_z[grep('NSS$', betas_nestednessG_z[[1]]), ], 
             digits = 2, 
             caption = "Factors affecting fruits' nutritional functional space. Average values, 95% CIs, and probabilities of directional effects for the relationships assessed in the generative model. An effect is marked as TRUE if P(effect > 0) > 0.7 or P(effect < 0) > 0.7 or.")

zzz <- 
  sapply(c(1:3, 6,8), FUN = 
         function(x) {
           ((abs(betas_nestednessG_z[grep('NSS$', betas_nestednessG_z[[1]]), ]$Mean[7]) -
              abs(betas_nestednessG_z[grep('NSS$', betas_nestednessG_z[[1]]), ]$Mean[x])) /
             abs(betas_nestednessG_z[grep('NSS$', betas_nestednessG_z[[1]]), ]$Mean[x])) * 100
         })

mean(zzz); sqrt(sd(zzz)/length(zzz))

tibble(var = betas_nestednessG_z[grep('NSS$', betas_nestednessG_z[[1]]), ][[1]][c(1:3, 6,8)], 
       `Importance of fruit's functional space` = zzz)



effect_temp_NSS <- 
  causal_effects(
    pars = post_nestednessG, 
    alpha_1 = 'alpha_NSS',
    f_1 = 'f_NSS', 
    beta_1 = 'beta_temp_NSS', 
    var = dat$z_temperature,
    y_obs = d$nestedness,
    level_1 = T
  )


quantile(effect_temp_NSS$relative_change, 
         c(0.025, 0.5, 0.975))

unique(effect_temp_NSS$p_effect)


effect_rain_NSS <- 
  causal_effects(
    pars = post_nestednessG, 
    alpha_1 = 'alpha_NSS',
    f_1 = 'f_NSS', 
    beta_1 = 'beta_rain_NSS', 
    var = dat$z_rainfall,
    y_obs = d$nestedness,
    level_1 = T
  )


quantile(effect_rain_NSS$relative_change, 
         c(0.025, 0.5, 0.975))

unique(effect_rain_NSS$p_effect)


effect_fruitAB_NSS <- 
  causal_effects(
    pars = post_nestednessG, 
    alpha_1 = 'alpha_NSS',
    f_1 = 'f_NSS', 
    beta_1 = 'beta_fruitAB_NSS', 
    var = dat$fruit_abun,
    y_obs = d$nestedness,
    level_1 = T
  )



quantile(effect_fruitAB_NSS$relative_change, 
         c(0.025, 0.5, 0.975))

unique(effect_fruitAB_NSS$p_effect)


effect_mor_NSS <- 
  causal_effects(
    pars = post_nestednessG, 
    alpha_1 = 'alpha_NSS',
    f_1 = 'f_NSS', 
    beta_1 = 'beta_mor_NSS', 
    var = dat$HV_plants_morfo,
    y_obs = d$nestedness,
    level_1 = T
  )



quantile(effect_mor_NSS$relative_change, 
         c(0.025, 0.5, 0.975))

unique(effect_mor_NSS$p_effect)


effect_PLAN_NSS <- 
  causal_effects(
    pars = post_nestednessG, 
    alpha_1 = 'alpha_NSS',
    f_1 = 'f_NSS', 
    beta_1 = 'beta_PLAN_NSS', 
    var = dat$HV_plant,
    y_obs = d$nestedness,
    level_1 = T
  )



quantile(effect_PLAN_NSS$relative_change, 
         c(0.025, 0.5, 0.975))

unique(effect_PLAN_NSS$p_effect)


effect_BIRD_NSS <- 
  causal_effects(
    pars = post_nestednessG, 
    alpha_1 = 'alpha_NSS',
    f_1 = 'f_NSS', 
    beta_1 = 'beta_BIRD_NSS', 
    var = dat$HV_bird,
    y_obs = d$nestedness,
    level_1 = T
  )


quantile(effect_BIRD_NSS$relative_change, 
         c(0.025, 0.5, 0.975))

unique(effect_BIRD_NSS$p_effect)

direct_nss <- rbind(effect_temp_NSS, 
                    effect_rain_NSS, 
                    effect_fruitAB_NSS, 
                    effect_mor_NSS, 
                    effect_PLAN_NSS, 
                    effect_BIRD_NSS)

direct_nss$effect <- as.factor(direct_nss$effect)

causal_plot_NSS_DIRECT <- 
  direct_nss |> 
    ggplot(aes(x = relative_change, y = effect, fill = after_stat(x))) +
    geom_density_ridges(scale = 1, rel_min_height = 0.005,
                        alpha = 0.5,
                        color = 'lightblue3',
                        fill = 'lightblue', 
                        linewidth = 0.25) +
    # geom_density_ridges_gradient(scale = 1.15, rel_min_height = 0.01, 
    #                              linewidth = 0.1) +
    # scale_fill_viridis_c(name = "Direction", option = "C") +
    #scale_fill_viridis_c(name = "", option = "C") +
    geom_vline(xintercept = 0, linetype = 3, linewidth = 0.35, color = 'tan1') +
    scale_y_discrete(labels = c(expression('Bifs'%->%"N"),
                                expression('FP'%->%"N"),
                                expression('Mifs'%->%"N"),
                                expression('FTifs'%->%"N"),
                                expression('R'%->%"N"),
                                expression('T'%->%"N"))) +
    labs(x = 'Relative causal effect (%)', x = '') +
    theme_classic() + 
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          axis.line = element_line(linewidth = 0.25),
          axis.ticks = element_line(linewidth = 0.25),
          legend.position = 'none',
          #legend.position = c(0.925, 0.8),
          legend.background = element_blank(),
          legend.key.size = unit(0.25, 'cm'),
          legend.text = element_text(size = 5),
          legend.title = element_text(size = 6),
          axis.title.y = element_blank(),
          axis.text = element_text(size = 8.5))
  


causal_plot_NSS_DIRECT

p_direct_nss <- 
  sapply(levels(direct_nss$effect), FUN = 
         function(i) {
           sapply(levels(direct_nss$effect), FUN = 
                    function(j) {
                      i <- direct_nss[direct_nss$effect == i, ]$relative_change
                      j <- direct_nss[direct_nss$effect == j, ]$relative_change
                      mean(i < j)
                    })
         })

colnames(p_direct_nss) <- c('Bifs', 'Fp', 'Mifs', 'FTifs', 'R', 'T')
rownames(p_direct_nss) <- c('Bifs', 'Fp', 'Mifs', 'FTifs', 'R', 'T')


knitr::kable(lower.tri.remove(p_direct_nss), 
             digits = 3, 
             caption = 'Probability that a causal direct effect on the row is lower than the effect on the columns to explaing network nestedness')

# Bottom-up effects (temperature)

# temp    fruitS  PLAN_HV NSS    
# temp    nut_HV  PLAN_HV NSS    
# temp    BIRD_HV PLAN_HV NSS    
# temp    BIRD_HV NSS    
# temp NSS 

bottom_up_tem_NSS <- 
  lapply(1, FUN = 
           function(x) {
             
             f_nut <- grep('f_nut_HV', colnames(post_nestednessG$f))
             alpha_nut <- grep('alpha_nut_HV', colnames(post_nestednessG$alpha))
             
             f_Frs <- grep('f_fruitS', colnames(post_nestednessG$f))
             alpha_Frs <- grep('alpha_fruitS', colnames(post_nestednessG$alpha))
             
             f_FT <- grep('f_PLAN_HV', colnames(post_nestednessG$f))
             alpha_FT <- grep('alpha_PLAN_HV', colnames(post_nestednessG$alpha))
             
             f_BIRD <- grep('f_BIRD_HV', colnames(post_nestednessG$f))
             alpha_BIRD <- grep('alpha_BIRD_HV', colnames(post_nestednessG$alpha))
             
             f_NSS <- grep('f_NSS', colnames(post_nestednessG$f))
             alpha_NSS <- grep('alpha_NSS', colnames(post_nestednessG$alpha))
             
             contrast <- 
               lapply(c(min, max), FUN = 
                        function(the_fun) {
                          est <- 
                            with(post_nestednessG, 
                                 {
                                   # temp -> nut
                                   mu_nut <- 
                                     apply(alpha[, alpha_nut], 1, mean) +
                                     apply(f[, f_nut], 1, mean) +
                                     beta$beta_temp_nut_HV * the_fun(dat$z_temperature) 
                                   
                                   # temp -> Fr. S
                                   mu_FrS <-
                                     apply(alpha[, alpha_Frs], 1, mean) +
                                     apply(f[, f_Frs], 1, mean) +
                                     beta$beta_temp_fruitS * the_fun(dat$z_temperature)
                                   
                                   # temp -> Birds ifs
                                   mu_BIRD <-
                                     apply(alpha[, alpha_BIRD], 1, mean) +
                                     apply(f[, f_BIRD], 1, mean) +
                                     beta$beta_temp_fruitS * the_fun(dat$z_temperature)
                                   
                                   # FT ifs
                                   
                                   mu_FT <-
                                     apply(alpha[, alpha_FT], 1, mean) +
                                     apply(f[, f_FT], 1, mean) +
                                     beta$beta_temp_PLAN_HV * the_fun(dat$z_temperature) +
                                     beta$beta_BIRD_PLAN_HV * mu_BIRD +
                                     beta$beta_nut_PLAN_HV * mu_nut +
                                     beta$beta_fruitS_PLAN_HV * mu_FrS
                                   
                                   # temp -> nut -> modu
                                   # temp -> Fr. S -> FT ifs -> modu
                                   # temp -> Birds ifs -> FT ifs -> modu
                                   
                                   mu_NSS <-
                                     apply(alpha[, alpha_NSS], 1, mean) +
                                     apply(f[, f_NSS], 1, mean) +
                                     beta$beta_temp_NSS * the_fun(dat$z_temperature) +
                                     beta$beta_BIRD_NSS * mu_BIRD +
                                     beta$beta_nut_NSS * mu_nut +
                                     beta$beta_fruitS_NSS * mu_FrS +
                                     beta$beta_PLAN_NSS * mu_FT
                                   
                                   mean(d$nestedness) + mu_NSS * sd(d$nestedness)
                                   
                                 }
                            )
                          est
                        })
             
             
             contrast <- as.data.frame(do.call('cbind', contrast))
             colnames(contrast) <- c('min', 'max')
             contrast$contrast <- contrast$max - contrast$min
             contrast$p_effect <- mean(contrast$contrast > 0)
             contrast$relative_change <- ((contrast$contrast * 100) / mean(d$nestedness))
             contrast$effect <- 'total_temp'
             as_tibble(contrast)
             
           })[[1]]


quantile(bottom_up_tem_NSS$relative_change, c(0.025, 0.5, 0.975))

mean(bottom_up_tem_NSS$relative_change < 0)

# Bottom-up effects (rainfall)

# rain    fruitAB nut_HV  PLAN_HV NSS    
# rain    fruitAB mor_HV  PLAN_HV NSS    
# rain    fruitAB mor_HV  NSS    
# rain    fruitAB BIRD_HV PLAN_HV NSS    
# rain    fruitAB BIRD_HV NSS    
# rain    fruitAB PLAN_HV NSS    
# rain    fruitAB NSS    
# rain    fruitS  PLAN_HV NSS    
# rain    BIRD_HV PLAN_HV NSS    
# rain    BIRD_HV NSS    
# rain    PLAN_HV NSS    
# rain NSS 

bottom_up_rain_NSS <- 
  lapply(1, FUN = 
           function(x) {
             
             f_mor <- grep('f_mor_HV', colnames(post_nestednessG$f))
             alpha_mor <- grep('alpha_mor_HV', colnames(post_nestednessG$alpha))
             
             f_nut <- grep('f_nut_HV', colnames(post_nestednessG$f))
             alpha_nut <- grep('alpha_nut_HV', colnames(post_nestednessG$alpha))
             
             f_Frs <- grep('f_fruitS', colnames(post_nestednessG$f))
             alpha_Frs <- grep('alpha_fruitS', colnames(post_nestednessG$alpha))
             
             f_Fp <- grep('f_fruitAB', colnames(post_nestednessG$f))
             alpha_Fp <- grep('alpha_fruitAB', colnames(post_nestednessG$alpha))
             
             f_FT <- grep('f_PLAN_HV', colnames(post_nestednessG$f))
             alpha_FT <- grep('alpha_PLAN_HV', colnames(post_nestednessG$alpha))
             
             f_BIRD <- grep('f_BIRD_HV', colnames(post_nestednessG$f))
             alpha_BIRD <- grep('alpha_BIRD_HV', colnames(post_nestednessG$alpha))
             
             f_NSS <- grep('f_NSS', colnames(post_nestednessG$f))
             alpha_NSS <- grep('alpha_NSS', colnames(post_nestednessG$alpha))
             
             contrast <- 
               lapply(c(min, max), FUN = 
                        function(the_fun) {
                          est <- 
                            with(post_nestednessG, 
                                 {
                                   # rain -> Fp
                                   mu_Fp <-
                                     apply(alpha[, alpha_Fp], 1, mean) +
                                     apply(f[, f_Fp], 1, mean) +
                                     beta$beta_rain_fruitAB * the_fun(dat$z_rainfall)
                                   
                                   # rain -> nut
                                   # Fp -> nut
                                   mu_nut <- 
                                     apply(alpha[, alpha_nut], 1, mean) +
                                     apply(f[, f_nut], 1, mean) +
                                     beta$beta_rain_nut_HV * the_fun(dat$z_rainfall) +
                                     beta$beta_fruitAB_nut_HV * mu_Fp
                                   
                                   # rain -> morph
                                   # Fp -> morph
                                   mu_morp <- 
                                     apply(alpha[, alpha_mor], 1, mean) +
                                     apply(f[, f_mor], 1, mean) +
                                     beta$beta_rain_mor_HV * the_fun(dat$z_rainfall) +
                                     beta$beta_fruitAB_mor_HV * mu_Fp
                                   
                                   # rain -> Birds ifs
                                   mu_BIRD <-
                                     apply(alpha[, alpha_BIRD], 1, mean) +
                                     apply(f[, f_BIRD], 1, mean) +
                                     beta$beta_rain_fruitS * the_fun(dat$z_rainfall)
                                   
                                   # rain -> Fr. S
                                   mu_FrS <-
                                     apply(alpha[, alpha_Frs], 1, mean) +
                                     apply(f[, f_Frs], 1, mean) +
                                     beta$beta_rain_fruitS * the_fun(dat$z_rainfall)
                                   
                                   # FT ifs
                                   mu_FT <-
                                     apply(alpha[, alpha_FT], 1, mean) +
                                     apply(f[, f_FT], 1, mean) +
                                     beta$beta_rain_PLAN_HV * the_fun(dat$z_rainfall) +
                                     beta$beta_BIRD_PLAN_HV * mu_BIRD +
                                     beta$beta_nut_PLAN_HV * mu_nut +
                                     beta$beta_fruitS_PLAN_HV * mu_FrS +
                                     beta$beta_fruitAB_PLAN_HV * mu_Fp +
                                     beta$beta_mor_PLAN_HV * mu_morp
                                   
                                   # temp -> nut -> modu
                                   # temp -> Fr. S -> FT ifs -> modu
                                   # temp -> Birds ifs -> FT ifs -> modu
                                   
                                   mu_NSS <-
                                     apply(alpha[, alpha_NSS], 1, mean) +
                                     apply(f[, f_NSS], 1, mean) +
                                     beta$beta_rain_NSS * the_fun(dat$z_temperature) +
                                     beta$beta_BIRD_NSS * mu_BIRD +
                                     beta$beta_mor_NSS * mu_morp +
                                     beta$beta_nut_NSS * mu_nut +
                                     beta$beta_fruitS_NSS * mu_FrS +
                                     beta$beta_fruitAB_NSS * mu_Fp +
                                     beta$beta_PLAN_NSS * mu_FT
                                   
                                   mean(d$nestedness) + mu_NSS * sd(d$nestedness)
                                   
                                 }
                            )
                          est
                        })
             
             
             contrast <- as.data.frame(do.call('cbind', contrast))
             colnames(contrast) <- c('min', 'max')
             contrast$contrast <- contrast$max - contrast$min
             contrast$p_effect <- mean(contrast$contrast < 0)
             contrast$relative_change <- ((contrast$contrast * 100) / mean(d$nestedness))
             contrast$effect <- 'total_Rain'
             as_tibble(contrast)
             
           })[[1]]

quantile(bottom_up_rain_NSS$relative_change, c(0.025, 0.5, 0.975))

mean(bottom_up_rain_NSS$relative_change < 0)

causal_plot_NSS_total <- 
  rbind(bottom_up_rain_NSS,
      bottom_up_tem_NSS) |> 
  ggplot(aes(x = relative_change, y = effect, fill = after_stat(x))) +
  geom_density_ridges(scale = 0.9, rel_min_height = 0.005,
                      alpha = 0.5,
                      color = 'lightblue3',
                      fill = 'lightblue', 
                      linewidth = 0.25) +
  # geom_density_ridges_gradient(scale = 1.15, rel_min_height = 0.01, 
  #                              linewidth = 0.1) +
  # scale_fill_viridis_c(name = "Direction", option = "C") +
  #scale_fill_viridis_c(name = "", option = "C") +
  geom_vline(xintercept = 0, linetype = 3, linewidth = 0.35, color = 'tan1') +
  scale_y_discrete(labels = c(expression('R (bottom-up)'%->%"N"),
                              expression('T (bottom-up)'%->%"N"))) +
  labs(x = 'Relative causal effect (%)', x = '') +
  theme_classic() + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.line = element_line(linewidth = 0.25),
        axis.ticks = element_line(linewidth = 0.25),
        legend.position = 'none',
        #legend.position = c(0.925, 0.8),
        legend.background = element_blank(),
        legend.key.size = unit(0.25, 'cm'),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 6),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 8.5))

causal_plot_NSS_total

mean(effect_rain_NSS$relative_change > bottom_up_rain_NSS$relative_change)

quantile(effect_rain_NSS$relative_change - bottom_up_rain_NSS$relative_change, c(0.025, 0.5, 0.975))

mean(effect_temp_NSS$relative_change > bottom_up_tem_NSS$relative_change)

quantile(effect_temp_NSS$relative_change - bottom_up_tem_NSS$relative_change, 
         c(0.025, 0.5, 0.975))

causal_plot_NSS_total_CONT <- 
  rbind(tibble(val = effect_rain_NSS$relative_change - 
               bottom_up_rain_NSS$relative_change, 
             contrast = 'Direct vs bottom-up effect (R)'), 
      tibble(val = effect_temp_NSS$relative_change - 
               bottom_up_tem_NSS$relative_change, 
             contrast = 'Direct vs bottom-up effect (T)')) |> 
  ggplot(aes(x = val), 
         color = 'ligth') +
  geom_density(alpha = 0.5,
               color = 'lightblue3',
               fill = 'lightblue', 
               linewidth = 0.25) +
  facet_wrap(~ contrast, ncol = 2, scales = 'free') +
  geom_vline(xintercept = 0, linetype = 3, linewidth = 0.35, color = 'tan1') +
  labs(x = 'Relative causal effect (%)', y = 'Density') +
  theme_classic() + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.line = element_line(linewidth = 0.25),
        axis.ticks = element_line(linewidth = 0.25),
        legend.position = 'none',
        #legend.position = c(0.925, 0.8),
        legend.background = element_blank(),
        legend.key.size = unit(0.25, 'cm'),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 6),
        axis.text = element_text(size = 8.5))

causal_plot_NSS_total_CONT

GP_nestedness <- mod_nestednessG$draws(c('sigma_NSS', 
                                         'sigma_f_NSS', 
                                         'gamma_NSS', 
                                         'eta_site_NSS'), 
                                       format = 'df')

total_sigma_NSS <- 
  GP_nestedness$sigma_NSS^2 + 
  GP_nestedness$sigma_f_NSS^2 +
  GP_nestedness$eta_site_NSS^2
prop_GP_temp_NSS <- GP_nestedness$sigma_f_NSS^2 / total_sigma_NSS
prop_GP_site_NSS <- GP_nestedness$eta_site_NSS^2 / total_sigma_NSS
prop_PRED_NSS <- GP_nestedness$sigma_NSS^2 / total_sigma_NSS

GP_NSS <- 
  tibble(component = c('Tiemporal processes', 
                     'Spatial processes', 
                     'Fized effects'), 
       mu = c(mean(prop_GP_temp_NSS), 
              mean(prop_GP_site_NSS), 
              mean(prop_PRED_NSS)), 
       sd = c(sd(prop_GP_temp_NSS), 
              sd(prop_GP_site_NSS), 
              sd(prop_PRED_NSS)))

knitr::kable(GP_NSS, 
             digits = 3, 
             caption = 'Variance decomposition of the components explaning intra-annual pattern of nestedness in seed dispersal networks. Average and SD values of proportions are provided')


mean(GP_nestedness$gamma_NSS); sd(GP_nestedness$gamma_NSS)

knitr::kable(betas_MODU_G_z[grep('MODU$', betas_MODU_G_z[[1]]), ], 
             digits = 2, 
             caption = "Factors affecting fruits' nutritional functional space. Average values, 95% CIs, and probabilities of directional effects for the relationships assessed in the generative model. An effect is marked as TRUE if P(effect > 0) > 0.7 or P(effect < 0) > 0.7 or.")

zzz <- 
  sapply(c(5, 7), FUN = 
         function(x) {
           ((abs(betas_MODU_G_z[grep('MODU$', betas_MODU_G_z[[1]]), ]$Mean[2]) -
              abs(betas_MODU_G_z[grep('MODU$', betas_MODU_G_z[[1]]), ]$Mean[x])) /
             abs(betas_MODU_G_z[grep('MODU$', betas_MODU_G_z[[1]]), ]$Mean[x])) * 100
         })

mean(zzz); sqrt(sd(zzz)/length(zzz))

tibble(var = betas_MODU_G_z[grep('MODU$', betas_MODU_G_z[[1]]), ][[1]][c(5, 7)], 
       `Importance of rainfall` = zzz)



effect_rain_MODU <- 
  causal_effects(
    pars = post_MODU_G, 
    alpha_1 = 'alpha_MODU',
    f_1 = 'f_MODU', 
    beta_1 = 'beta_rain_MODU', 
    var = dat$z_rainfall,
    y_obs = d$modularity,
    level_1 = T
  )


quantile(effect_rain_MODU$relative_change, 
         c(0.025, 0.5, 0.975))

unique(effect_rain_MODU$p_effect)


effect_nut_MODU <- 
  causal_effects(
    pars = post_MODU_G, 
    alpha_1 = 'alpha_MODU',
    f_1 = 'f_MODU', 
    beta_1 = 'beta_nut_MODU', 
    var = dat$HV_plants_nut,
    y_obs = d$modularity,
    level_1 = T
  )


quantile(effect_nut_MODU$relative_change, 
         c(0.025, 0.5, 0.975))

unique(effect_nut_MODU$p_effect)


effect_PLAN_MODU <- 
  causal_effects(
    pars = post_MODU_G, 
    alpha_1 = 'alpha_MODU',
    f_1 = 'f_MODU', 
    beta_1 = 'beta_PLAN_MODU', 
    var = dat$HV_plant,
    y_obs = d$modularity,
    level_1 = T
  )


quantile(effect_PLAN_MODU$relative_change, 
         c(0.025, 0.5, 0.975))

unique(effect_PLAN_MODU$p_effect)

causal_plot_MODU_DIRECT <- 
  rbind(effect_rain_MODU,
        effect_nut_MODU, 
        effect_PLAN_MODU) |>
  ggplot(aes(x = relative_change, y = effect, fill = after_stat(x))) +
  geom_density_ridges(scale = 1, rel_min_height = 0.005,
                      alpha = 0.5,
                      color = 'lightblue3',
                      fill = 'lightblue', 
                      linewidth = 0.25) +
  # geom_density_ridges_gradient(scale = 1.15, rel_min_height = 0.01, 
  #                              linewidth = 0.1) +
  # scale_fill_viridis_c(name = "Direction", option = "C") +
  #scale_fill_viridis_c(name = "", option = "C") +
  geom_vline(xintercept = 0, linetype = 3, linewidth = 0.35, color = 'tan1') +
  scale_y_discrete(labels = c(expression('Nifs'%->%"M"),
                              expression('FTifs'%->%"M"),
                              expression('R'%->%"M"))) +
  lims(x = c(-100, 300)) +
  labs(x = 'Relative causal effect (%)', x = '') +
  theme_classic() + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.line = element_line(linewidth = 0.25),
        axis.ticks = element_line(linewidth = 0.25),
        legend.position = 'none',
        #legend.position = c(0.925, 0.8),
        legend.background = element_blank(),
        legend.key.size = unit(0.25, 'cm'),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 6),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 8.5))

causal_plot_MODU_DIRECT

causal_plot_MODU_DIRECT

# Bottom-up effect (temperature)

# temp -> nut -> modu
# temp -> Fr. S -> FT ifs -> modu
# temp -> Birds ifs -> FT ifs -> modu

bottom_up_tem_MODU <- 
  lapply(1, FUN = 
         function(x) {
           
           f_nut <- grep('f_nut_HV', colnames(post_MODU_G$f))
           alpha_nut <- grep('alpha_nut_HV', colnames(post_MODU_G$alpha))
           
           f_Frs <- grep('f_fruitS', colnames(post_MODU_G$f))
           alpha_Frs <- grep('alpha_fruitS', colnames(post_MODU_G$alpha))
           
           f_FT <- grep('f_PLAN_HV', colnames(post_MODU_G$f))
           alpha_FT <- grep('alpha_PLAN_HV', colnames(post_MODU_G$alpha))
           
           f_BIRD <- grep('f_BIRD_HV', colnames(post_MODU_G$f))
           alpha_BIRD <- grep('alpha_BIRD_HV', colnames(post_MODU_G$alpha))
           
           f_MODU <- grep('f_MODU', colnames(post_MODU_G$f))
           alpha_MODU <- grep('alpha_MODU', colnames(post_MODU_G$alpha))
           
           contrast <- 
             lapply(c(min, max), FUN = 
                      function(the_fun) {
                        est <- 
                          with(post_MODU_G, 
                               {
                                 # temp -> nut
                                 mu_nut <- 
                                   apply(alpha[, alpha_nut], 1, mean) +
                                   apply(f[, f_nut], 1, mean) +
                                   beta$beta_temp_nut_HV * the_fun(dat$z_temperature) 
                                 
                                 # temp -> Fr. S
                                 mu_FrS <-
                                   apply(alpha[, alpha_Frs], 1, mean) +
                                   apply(f[, f_Frs], 1, mean) +
                                   beta$beta_temp_fruitS * the_fun(dat$z_temperature)
                                 
                                 # temp -> Birds ifs
                                 mu_BIRD <-
                                   apply(alpha[, alpha_BIRD], 1, mean) +
                                   apply(f[, f_BIRD], 1, mean) +
                                   beta$beta_temp_fruitS * the_fun(dat$z_temperature)
                                 
                                 # FT ifs
                                 
                                 mu_FT <-
                                   apply(alpha[, alpha_FT], 1, mean) +
                                   apply(f[, f_FT], 1, mean) +
                                   beta$beta_temp_PLAN_HV * the_fun(dat$z_temperature) +
                                   beta$beta_BIRD_PLAN_HV * mu_BIRD +
                                   beta$beta_nut_PLAN_HV * mu_nut +
                                   beta$beta_fruitS_PLAN_HV * mu_FrS
                                 
                                 # temp -> nut -> modu
                                 # temp -> Fr. S -> FT ifs -> modu
                                 # temp -> Birds ifs -> FT ifs -> modu
                                 
                                 mu_MODU <-
                                   apply(alpha[, alpha_MODU], 1, mean) +
                                   apply(f[, f_MODU], 1, mean) +
                                   beta$beta_temp_MODU * the_fun(dat$z_temperature) +
                                   beta$beta_BIRD_MODU * mu_BIRD +
                                   beta$beta_nut_MODU * mu_nut +
                                   beta$beta_fruitS_MODU * mu_FrS +
                                   beta$beta_PLAN_MODU * mu_FT
                                 
                                 inv_logit(mu_MODU)
                                 
                               }
                          )
                        est
                      })
           
           
           contrast <- as.data.frame(do.call('cbind', contrast))
           colnames(contrast) <- c('min', 'max')
           contrast$contrast <- contrast$max - contrast$min
           contrast$p_effect <- mean(contrast$contrast < 0)
           contrast$relative_change <- ((contrast$contrast * 100) / mean(dat$modularity))
           contrast$effect <- 'total_temp'
           as_tibble(contrast)
           
         })[[1]]

quantile(bottom_up_tem_MODU$relative_change, c(0.025, 0.5, 0.975))

mean(bottom_up_tem_MODU$contrast < 0)

# Bottom-up effect (rainfall)

# rain -> fruitAB -> nut_HV -> PLAN_HV MODU   
# rain ->  fruitAB -> nut_HV  -> MODU   
# rain  ->  fruitAB -> mor_HV ->  PLAN_HV -> MODU   
# rain -> fruitAB -> BIRD_HV -> PLAN_HV -> MODU   
# rain -> fruitAB -> PLAN_HV -> MODU   
# rain -> fruitS -> PLAN_HV MODU   
# rain -> BIRD_HV -> PLAN_HV -> MODU   
# rain -> PLAN_HV -> MODU   
# rain -> MODU

bottom_up_rain_MODU <- 
  lapply(1, FUN = 
           function(x) {
             
             f_mor <- grep('f_mor_HV', colnames(post_MODU_G$f))
             alpha_mor <- grep('alpha_mor_HV', colnames(post_MODU_G$alpha))
             
             f_nut <- grep('f_nut_HV', colnames(post_MODU_G$f))
             alpha_nut <- grep('alpha_nut_HV', colnames(post_MODU_G$alpha))
             
             f_Frs <- grep('f_fruitS', colnames(post_MODU_G$f))
             alpha_Frs <- grep('alpha_fruitS', colnames(post_MODU_G$alpha))
             
             f_Fp <- grep('f_fruitAB', colnames(post_MODU_G$f))
             alpha_Fp <- grep('alpha_fruitAB', colnames(post_MODU_G$alpha))
             
             f_FT <- grep('f_PLAN_HV', colnames(post_MODU_G$f))
             alpha_FT <- grep('alpha_PLAN_HV', colnames(post_MODU_G$alpha))
             
             f_BIRD <- grep('f_BIRD_HV', colnames(post_MODU_G$f))
             alpha_BIRD <- grep('alpha_BIRD_HV', colnames(post_MODU_G$alpha))
             
             f_MODU <- grep('f_MODU', colnames(post_MODU_G$f))
             alpha_MODU <- grep('alpha_MODU', colnames(post_MODU_G$alpha))
             
             contrast <- 
               lapply(c(min, max), FUN = 
                        function(the_fun) {
                          est <- 
                            with(post_MODU_G, 
                                 {
                                   # rain -> Fp
                                   mu_Fp <-
                                     apply(alpha[, alpha_Fp], 1, mean) +
                                     apply(f[, f_Fp], 1, mean) +
                                     beta$beta_rain_fruitAB * the_fun(dat$z_rainfall)
                                   
                                   # rain -> nut
                                   # Fp -> nut
                                   mu_nut <- 
                                     apply(alpha[, alpha_nut], 1, mean) +
                                     apply(f[, f_nut], 1, mean) +
                                     beta$beta_rain_nut_HV * the_fun(dat$z_rainfall) +
                                     beta$beta_fruitAB_nut_HV * mu_Fp
                                   
                                   # rain -> morph
                                   # Fp -> morph
                                   mu_morp <- 
                                     apply(alpha[, alpha_mor], 1, mean) +
                                     apply(f[, f_mor], 1, mean) +
                                     beta$beta_rain_mor_HV * the_fun(dat$z_rainfall) +
                                     beta$beta_fruitAB_mor_HV * mu_Fp
                                   
                                   # rain -> Birds ifs
                                   mu_BIRD <-
                                     apply(alpha[, alpha_BIRD], 1, mean) +
                                     apply(f[, f_BIRD], 1, mean) +
                                     beta$beta_rain_fruitS * the_fun(dat$z_rainfall)
                                   
                                   # rain -> Fr. S
                                   mu_FrS <-
                                     apply(alpha[, alpha_Frs], 1, mean) +
                                     apply(f[, f_Frs], 1, mean) +
                                     beta$beta_rain_fruitS * the_fun(dat$z_rainfall)
                                   
                                   # FT ifs
                                   mu_FT <-
                                     apply(alpha[, alpha_FT], 1, mean) +
                                     apply(f[, f_FT], 1, mean) +
                                     beta$beta_rain_PLAN_HV * the_fun(dat$z_rainfall) +
                                     beta$beta_BIRD_PLAN_HV * mu_BIRD +
                                     beta$beta_nut_PLAN_HV * mu_nut +
                                     beta$beta_fruitS_PLAN_HV * mu_FrS +
                                     beta$beta_fruitAB_PLAN_HV * mu_Fp +
                                     beta$beta_mor_PLAN_HV * mu_morp
                                   
                                   # temp -> nut -> modu
                                   # temp -> Fr. S -> FT ifs -> modu
                                   # temp -> Birds ifs -> FT ifs -> modu
                                   
                                   mu_MODU <-
                                     apply(alpha[, alpha_MODU], 1, mean) +
                                     apply(f[, f_MODU], 1, mean) +
                                     beta$beta_rain_MODU * the_fun(dat$z_temperature) +
                                     beta$beta_BIRD_MODU * mu_BIRD +
                                     beta$beta_mor_MODU * mu_morp
                                     beta$beta_nut_MODU * mu_nut +
                                     beta$beta_fruitS_MODU * mu_FrS +
                                     beta$beta_fruitAB_MODU * mu_Fp +
                                     beta$beta_PLAN_MODU * mu_FT
                                   
                                   inv_logit(mu_MODU)
                                   
                                 }
                            )
                          est
                        })
             
             
             contrast <- as.data.frame(do.call('cbind', contrast))
             colnames(contrast) <- c('min', 'max')
             contrast$contrast <- contrast$max - contrast$min
             contrast$p_effect <- mean(contrast$contrast > 0)
             contrast$relative_change <- ((contrast$contrast * 100) / mean(dat$modularity))
             contrast$effect <- 'total_Rain'
             as_tibble(contrast)
             
           })[[1]]


quantile(bottom_up_rain_MODU$relative_change, c(0.025, 0.5, 0.975))

mean(bottom_up_rain_MODU$contrast > 0)


causal_plot_MODU_total <- 
  rbind(bottom_up_rain_MODU, 
      bottom_up_tem_MODU) |> 
  ggplot(aes(x = relative_change, y = effect, fill = after_stat(x))) +
  geom_density_ridges(scale = 0.9, rel_min_height = 0.005,
                      alpha = 0.5,
                      color = 'lightblue3',
                      fill = 'lightblue', 
                      linewidth = 0.25) +
  # geom_density_ridges_gradient(scale = 1.15, rel_min_height = 0.01, 
  #                              linewidth = 0.1) +
  # scale_fill_viridis_c(name = "Direction", option = "C") +
  #scale_fill_viridis_c(name = "", option = "C") +
  geom_vline(xintercept = 0, linetype = 3, linewidth = 0.35, color = 'tan1') +
  scale_y_discrete(labels = c(expression('R (bottom-up)'%->%"M"),
                              expression('T (bottom-up)'%->%"M"))) +
  labs(x = 'Relative causal effect (%)', x = '') +
  theme_classic() + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.line = element_line(linewidth = 0.25),
        axis.ticks = element_line(linewidth = 0.25),
        legend.position = 'none',
        #legend.position = c(0.925, 0.8),
        legend.background = element_blank(),
        legend.key.size = unit(0.25, 'cm'),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 6),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 8.5))

causal_plot_MODU_total

mean(effect_rain_MODU$relative_change > bottom_up_rain_MODU$relative_change)

quantile(effect_rain_MODU$relative_change - bottom_up_rain_MODU$relative_change, 
         c(0.025, 0.5, 0.975))

causal_plot_MODU_total_CONT <- 
  tibble(val = effect_rain_MODU$relative_change - 
         bottom_up_rain_MODU$relative_change, 
       contrast = 'Direct vs bottom-up effect (R)') |> 
  ggplot(aes(x = val), 
         color = 'ligth') +
  geom_density(alpha = 0.5,
               color = 'lightblue3',
               fill = 'lightblue', 
               linewidth = 0.25) +
  facet_wrap(~ contrast) +
  geom_vline(xintercept = 0, linetype = 3, linewidth = 0.35, color = 'tan1') +
  labs(x = 'Relative causal effect (%)', y = 'Density') +
  theme_classic() + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.line = element_line(linewidth = 0.25),
        axis.ticks = element_line(linewidth = 0.25),
        legend.position = 'none',
        #legend.position = c(0.925, 0.8),
        legend.background = element_blank(),
        legend.key.size = unit(0.25, 'cm'),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 6),
        axis.text = element_text(size = 8.5))

causal_plot_MODU_total_CONT

GP_modu <- mod_MODU_G$draws(c('gamma_MODU', 'prop_GP_temp', 
                              'prop_fixed', 'prop_GP_space'), format = 'df')

GP_MODU <- 
  do.call('rbind', 
        lapply(2:4, FUN = 
                 function(x) {
                   tibble(component = c('Temporal processes',
                                        'Fixed effects', 
                                        'Spatial processes')[x-1], 
                          mu = median(GP_modu[[x]]), 
                          sd = sd(GP_modu[[x]]))
                 }))

knitr::kable(GP_MODU, 
             digits = 3, 
             caption = 'Variance decomposition of the components explaning intra-annual pattern of modularity in seed dispersal networks. Average and SD of proportions are provided')


mean(GP_modu$gamma_MODU); sd(GP_modu$gamma_MODU)

knitr::kable(betas_H2_G_z[grep('H2$', betas_H2_G_z[[1]]), ], 
             digits = 2, 
             caption = "Factors affecting fruits' nutritional functional space. Average values, 95% CIs, and probabilities of directional effects for the relationships assessed in the generative model. An effect is marked as TRUE if P(effect > 0) > 0.7 or P(effect < 0) > 0.7 or.")

((abs(betas_H2_G[grep('H2$', betas_H2_G$var), ]$mu[2]) -
    abs(betas_H2_G[grep('H2$', betas_H2_G$var), ]$mu[1])) /
    abs(betas_H2_G[grep('H2$', betas_H2_G$var), ]$mu[1])) * 100



effect_rain_H2 <- 
  causal_effects(
    pars = post_H2_G, 
    alpha_1 = 'alpha_H2',
    f_1 = 'f_H2', 
    beta_1 = 'beta_rain_H2', 
    var = dat$z_rainfall,
    y_obs = d$H2,
    level_1 = T
  )


quantile(effect_rain_H2$relative_change, 
         c(0.025, 0.5, 0.975))

unique(effect_rain_H2$p_effect)


effect_temp_H2 <- 
  causal_effects(
    pars = post_H2_G, 
    alpha_1 = 'alpha_H2',
    f_1 = 'f_H2', 
    beta_1 = 'beta_temp_H2', 
    var = dat$z_temperature,
    y_obs = d$H2,
    level_1 = T
  )


quantile(effect_temp_H2$relative_change, 
         c(0.025, 0.5, 0.975))

unique(effect_temp_H2$p_effect)

causal_plot_H2 <- 
  rbind(effect_rain_H2,
        effect_temp_H2) |>
  ggplot(aes(x = relative_change, y = effect, fill = after_stat(x))) +
  geom_density_ridges(scale = 1, rel_min_height = 0.001,
                               alpha = 0.5,
                      color = 'lightblue3',
                      fill = 'lightblue', 
                      linewidth = 0.25) +
  # geom_density_ridges_gradient(scale = 1.15, rel_min_height = 0.01, 
  #                              linewidth = 0.1) +
  # scale_fill_viridis_c(name = "Direction", option = "C") +
  #scale_fill_viridis_c(name = "", option = "C") +
  geom_vline(xintercept = 0, linetype = 3, linewidth = 0.35, color = 'tan1') +
  scale_y_discrete(labels = c(expression('R'%->%"H2'"), 
                              expression('T'%->%"H2'"))) +
  #lims(y = c(-100, 224)) +
  labs(x = 'Relative causal effect (%)', x = '') +
  theme_classic() + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.line = element_line(linewidth = 0.25),
        axis.ticks = element_line(linewidth = 0.25),
        legend.position = 'none',
        #legend.position = c(0.925, 0.8),
        legend.background = element_blank(),
        legend.key.size = unit(0.25, 'cm'),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 6),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 8.5))

causal_plot_H2

GP_modu_H2 <- mod_H2_G$draws(c('gamma_H2', 'prop_GP_temp', 
                              'prop_fixed', 'prop_GP_space'), format = 'df')

GP_H2 <- do.call('rbind', 
        lapply(2:4, FUN = 
                 function(x) {
                   tibble(component = c('Temporal processes',
                                        'Fixed effects', 
                                        'Spatial processes')[x-1], 
                          mu = median(GP_modu_H2[[x]]), 
                          sd = sd(GP_modu_H2[[x]]))
                 }))


knitr::kable(GP_H2, 
             digits = 3, 
             caption = 'Variance decomposition of the components explaning intra-annual pattern of specialization in seed dispersal networks. Average and SD of proportions are provided')


mean(GP_modu_H2$gamma_H2); sd(GP_modu_H2$gamma_H2)

sessionInfo()
