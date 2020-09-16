cat('\f')
rm(list = ls())

setwd('/Users/robertladwig/Documents/DSI/b-odem')

library(tidyverse)
library(lubridate)
library(ggdark)
library(deSolve)
library(reshape2)
library(raster)
library(zoo)

lake.list <- c('Allequash', 'BigMuskellunge', 'Crystal', 'Fish', 'Mendota',
               'Monona', 'Sparkling', 'Trout', 'Wingra')

# for (lake.id in lake.list){
lake.id = lake.list[5]
  setwd(paste0(lake.id))
  years <- 1979:2019
  

  input <- readr::read_csv(
    paste0('input.txt'), 
    col_names=c('datetime', 'thermocline_depth', 'temperature_epi',
                'temperature_hypo', 'temperature_total', 'volume_total',
                'volume_epi', 'volume_hypo', 'area_thermocline',
                'area_surface', 'upper_meta', 'lower_meta',
                'year', 'day_of_year', 'z.max', 'wind','airtemp'),
    col_types=cols(datetime=col_datetime(), year=col_integer(), day_of_year=col_integer(), .default=col_double()))
  inyr <- filter(input, year %in% years)


  # obs <- read.table(
  #   paste0('observed.txt'),
  #   header=FALSE,
  #   sep=' ',
  #   as.is=TRUE) %>%
  #   t() %>%
  #   as_tibble(.name_repair='minimal') %>%
  #   setNames(., nm=c('dateint', 'DO_tot', 'DO_epi', 'DO_hypo')) %>%
  #   mutate(date = zoo::as.Date(dateint, origin='1979-04-01')) %>% # just guessing at origin and therefore at dates
  #   select(date, everything())
  obs <- read.table(
    paste0('observed.txt'),
    header=FALSE,
    sep=' ',
    as.is=TRUE) %>%
    t() %>%
    as_tibble(.name_repair='minimal') %>%
    setNames(., nm=c('dateint', 'DO_tot', 'DO_epi', 'DO_hypo')) %>%
    mutate(date = zoo::as.Date(dateint, origin='1979-04-01')) 
  obsyr <- filter(obs, lubridate::year(date) %in% years)

  
  idx1 = which(!is.na(inyr$thermocline_depth))[1]
  idx2 = rev(which(!is.na(inyr$thermocline_depth)))[1]
  inyr = inyr #[idx1:idx2,]
  inyr$strat <- ifelse(is.na(inyr$thermocline_depth),0,1)
  strat.pos <- c()
  for (ii in 1:length(inyr$strat)){
    if (inyr$strat[ii] == 1 && inyr$strat[ii-1] == 0){
      strat.pos <- append(strat.pos, ii)
    }
  }
  
  idy = match(obsyr$date,zoo::as.Date(inyr$datetime))
  idx = idy[!is.na(obsyr$DO_tot)]
  idxx = idy[!is.na(obsyr$DO_epi)]
  idz = match(idx, idy)
  idzz =match(idxx, idy)
  DO_obs_epi = rep(NA, length(inyr$datetime))
  DO_obs_epi[idxx] = obsyr$DO_epi[idzz]
  DO_obs_hyp = rep(NA, length(inyr$datetime))
  DO_obs_hyp[idxx] = obsyr$DO_hypo[idzz]
  DO_obs_tot = rep(NA, length(inyr$datetime))
  DO_obs_tot[idx] = obsyr$DO_tot[idz]
  
  dates = inyr$datetime
  
  library(lazyeval)
  inyr <- inyr[,-c(1)] %>% mutate_each( funs_( interp( ~replace(., is.na(.),0) ) ) )
  
  library(LakeMetabolizer)
  
  simdata <- tibble(
    DO_obs_epi = DO_obs_epi * 1000,
    DO_obs_hyp = DO_obs_hyp * 1000,
    DO_obs_tot = DO_obs_tot * 1000,
    day = seq(1, nrow(inyr))
  )
  strat.pos <- c()
  for (ii in 1:length(inyr$strat)){
    if (inyr$strat[ii] == 1 && inyr$strat[ii-1] == 0){
      strat.pos <- append(strat.pos, ii)
    }
  }
  shiftdays <- rep(0, nrow(inyr))
  shiftdays[strat.pos] <- 1
  
  dummyinput <- data.frame(
    theta0 = 1.08^(inyr$temperature_total - 20),
    theta1 = 1.08^(inyr$temperature_epi - 20),
    theta2 = 1.08^(inyr$temperature_hypo - 20),
    k600t = k600.2.kGAS.base(k.cole.base(inyr$wind),temperature = inyr$temperature_total, gas = "O2"),
    o2satt = o2.at.sat.base(temp = inyr$temperature_total, altitude = 450) * 1000,
    k600 = k600.2.kGAS.base(k.cole.base(inyr$wind),temperature = inyr$temperature_epi, gas = "O2"),
    o2sat = o2.at.sat.base(temp = inyr$temperature_epi, altitude = 450) * 1000,
    volume_epi = inyr$volume_epi,
    volume_tot = inyr$volume_total,
    area_epi = inyr$area_surface,
    volume_hyp = inyr$volume_hypo,
    area_hyp = inyr$area_thermocline,
    tddepth = inyr$thermocline_depth,
    wtr_epi = inyr$temperature_epi,
    wtr_hyp = inyr$temperature_hypo,
    wtr_tot = inyr$temperature_total,
    khalf = rep(500, nrow(inyr)), 
    stratified = inyr$strat,
    airtemp = inyr$airtemp,
    shiftdays = shiftdays
  )
  
  inyr$volume_epi[which(inyr$volume_epi==0)] = NA
  f1 <- function(dat) {
    N <- length(dat)
    na.pos <- which(is.na(dat))
    if (length(na.pos) %in% c(0, N)) {
      return(dat)
    }
    non.na.pos <- which(!is.na(dat))
    intervals  <- findInterval(na.pos, non.na.pos,
                               all.inside = TRUE)
    left.pos   <- non.na.pos[pmax(1, intervals)]
    right.pos  <- non.na.pos[pmin(N, intervals+1)]
    left.dist  <- na.pos - left.pos
    right.dist <- right.pos - na.pos
    
    dat[na.pos] <- ifelse(left.dist <= right.dist,
                          dat[left.pos], dat[right.pos])
    return(dat)
  }
  inyr$volume_epi = f1(inyr$volume_epi)

  inyr$volume_hypo[which(inyr$volume_hypo==0)] = NA
  inyr$volume_hypo = f1(inyr$volume_hypo)
  
  dummyinput$k600t[which(inyr$airtemp <= 0 & inyr$temperature_total <= 4)] = 1e-5
  dummyinput$volume_epi = ifelse(dummyinput$stratified == 0, dummyinput$volume_tot, dummyinput$volume_epi)
  dummyinput$volume_hyp = ifelse(dummyinput$stratified == 0, dummyinput$volume_tot, dummyinput$volume_hyp)
  # dummyinput$volchange_epi = c(diff(inyr$volume_epi),0)/c(inyr$volume_epi)#/inyr$volume_epi
  # dummyinput$volchange_hypo = c(diff(inyr$volume_hyp),0)/c(inyr$volume_hyp)#/inyr$volume_hypo
  dummyinput$max_depth = rep(25, nrow(inyr))
  dummyinput$volchange = 0
  for (i in 2:(nrow(inyr))){
    dummyinput$volchange[i] = inyr$volume_epi[i] - inyr$volume_epi[i-1]
  }

  
  ic_cond <- data.frame(    DO_epi_init = 15 * 1000 * dummyinput$volume_tot[1], 
                            DO_hyp_init = 15 * 1000 * dummyinput$volume_tot[1],
                            DO_tot_init = 15 * 1000 * dummyinput$volume_tot[1])
  
  
  times = seq(1, nrow(dummyinput))
  theta0 <- approxfun(x = times, y = dummyinput$theta0, method = "linear", rule = 2) 
  theta1 <- approxfun(x = times, y = dummyinput$theta1, method = "linear", rule = 2) 
  theta2 <- approxfun(x = times, y = dummyinput$theta2, method = "linear", rule = 2)
  k600t <- approxfun(x = times, y = dummyinput$k600t, method = "linear", rule = 2) 
  o2satt <- approxfun(x = times, y = dummyinput$o2satt, method = "linear", rule = 2) 
  k600 <- approxfun(x = times, y = dummyinput$k600, method = "linear", rule = 2) 
  o2sat <- approxfun(x = times, y = dummyinput$o2sat, method = "linear", rule = 2) 
  # volume_epi <- approxfun(x = times, y = dummyinput$volume_epi, method = "linear", rule = 2) 
  # volume_tot <- approxfun(x = times, y = dummyinput$volume_tot, method = "linear", rule = 2) 
  # volume_hyp <- approxfun(x = times, y = dummyinput$volume_hyp, method = "linear", rule = 2) 
  volume_epi <- approxfun(x = times, y = inyr$volume_epi, method = "linear", rule = 2) 
  volume_tot <- approxfun(x = times, y = inyr$volume_total, method = "linear", rule = 2) 
  volume_hyp <- approxfun(x = times, y = inyr$volume_hypo, method = "linear", rule = 2) 
  area_epi <- approxfun(x = times, y = dummyinput$area_epi, method = "linear", rule = 2) 
  area_hyp <- approxfun(x = times, y = dummyinput$area_hyp, method = "linear", rule = 2) 
  tddepth <- approxfun(x = times, y = dummyinput$tddepth, method = "linear", rule = 2) 
  wtr_epi <- approxfun(x = times, y = dummyinput$wtr_epi, method = "linear", rule = 2) 
  wtr_hyp <- approxfun(x = times, y = dummyinput$wtr_hyp, method = "linear", rule = 2) 
  wtr_tot <- approxfun(x = times, y = dummyinput$wtr_tot, method = "linear", rule = 2) 
  khalf <- approxfun(x = times, y = dummyinput$khalf, method = "linear", rule = 2) 
  stratified <- approxfun(x = times, y = dummyinput$stratified, method = "linear", rule = 2) 
  airtemp <- approxfun(x = times, y = dummyinput$airtemp, method = "linear", rule = 2) 
  shiftdays <- approxfun(x = times, y = dummyinput$shiftdays, method = "linear", rule = 2) 
  # volchange_epi <- approxfun(x = times, y = dummyinput$volchange_epi, method = "linear", rule = 2) 
  # volchange_hypo <- approxfun(x = times, y = dummyinput$volchange_hypo, method = "linear", rule = 2) 
  max_depth <- approxfun(x = times, y = dummyinput$max_depth, method = "linear", rule = 2) 
  volchange <- approxfun(x = times, y = dummyinput$volchange, method = "linear", rule = 2) 
  
  run_odem <- function(t, y, parms){
    
    NEP = params[1]
    SED = params[2]
    
    fNEP = 0
    fSED = 0
    fATM =  0
    fENTR1 = 0
    fENTR2 = 0
    
    if (floor(stratified(t)) == 0) {
      
      fNEP = NEP * theta0(t)
      fSED = -  SED *  (y[1]/(khalf(t) + y[1])) * theta0(t) / (volume_tot(t)/area_epi(t))
      fATM =  k600t(t)  *  (o2sat(t) - max(y[1]/volume_tot(t),0.01))  / (volume_tot(t)/area_epi(t)) 
      fENTR1 = 0
      fENTR2 = 0
      
      dDOdt_tot = ( fNEP - 
        fSED +
        fATM) 
      
      if(abs(dDOdt_tot)>abs(y[1])){
        if(dDOdt_tot < 0){
          flux_tot = - y[1];
        }else{
          flux_tot = dDOdt_tot  
        }
      }else{
        flux_tot = dDOdt_tot
      }
      
      flux_epi = flux_tot 
      flux_hypo = flux_tot
      

      
    } else if (stratified(t) == 1){
    
      
      if (volchange(t) >= 0){
        x_do1 = y[3] 
      } else {
        x_do1 = y[2] 
      }
      
      # if (delvol_hypo >= 0){
      #   x_do2 = y[2] / volume_epi(t)
      # } else {
      #   x_do2 = y[3] / volume_hyp(t)
      # }
      
      if (y[2] <= 0){
        atm_do = o2sat(t) 
      } else {
        atm_do = (o2sat(t) - y[2]/volume_epi(t))
      }
      
      fNEP = NEP *theta1(t)
      fSED = - SED*  (y[3]/(khalf(t) + y[3])) * theta2(t)  / max(volume_hyp(t)/area_hyp(t),1)#(max_depth(t) - tddepth(t))
      fATM =  k600(t) *  (o2sat(t) - max(y[2]/volume_epi(t),0.01)) / tddepth(t)
      fENTR1 = volchange(t) * x_do1 / volume_epi(t)
      fENTR2 =  - volchange(t) * x_do1 / volume_hyp(t)
      
      dDOdt_epi= fNEP-
        fATM +
        fENTR1
    
    if(abs(dDOdt_epi)>abs(y[2])){
      if(dDOdt_epi < 0){
        flux_epi = - y[2];
      }else{
        flux_epi = dDOdt_epi 
      }
    }else{
      flux_epi = dDOdt_epi 
    }
      dDOdt_hyp =  fSED +
        fENTR2
    
    if(abs(dDOdt_hyp)>abs(y[3])){
      if(dDOdt_hyp < 0){
        flux_hypo = - y[3];
      }else{
        flux_hypo = dDOdt_hyp 
      }
    }else{
      flux_hypo = dDOdt_hyp 
    }
    
    # flux_tot = (flux_epi / volume_epi(t) + flux_hypo / volume_hyp(t))/2 * volume_tot(t)
    flux_tot = (flux_epi + flux_hypo )/2 

  
    
    }
    write.table(matrix(c(fNEP, fSED, fATM, fENTR1, fENTR2,
                         flux_tot, flux_epi, flux_hypo, t), nrow=1), '/Users/robertladwig/Documents/DSI/odem/output.txt', append = TRUE,
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    return(list(c(flux_tot, flux_epi, flux_hypo)))
  }

  params <- c(0.1, 0.1)#c(300, 2500)
  
  file.remove('/Users/robertladwig/Documents/DSI/odem/output.txt')
  
  out <- deSolve::ode(times = seq(1, nrow(dummyinput), by = 1), y = as.numeric(ic_cond), func = run_odem, parms = params,
             method = "rk4")
  
  vol_tot = dummyinput$volume_tot*1000
  vol_epi = dummyinput$volume_epi*1000
  vol_hypo = dummyinput$volume_hyp*1000
  # vol_tot = inyr$volume_total*1000
  # vol_epi = inyr$volume_epi*1000
  # vol_hypo = inyr$volume_hypo*1000
  # 
  t = 100
  plot(out[1:t,1], out[1:t,2]/(vol_tot[1:t]))
  plot(out[1:t,1], out[1:t,3]/(vol_epi[1:t]))
  plot(out[1:t,1], out[1:t,4]/(vol_hypo[1:t]))
  
  plot(out[1:t,1], out[1:t,2])
  plot(out[1:t,1], out[1:t,3])
  plot(out[1:t,1], out[1:t,4])
  
  output = data.frame('time' = out[,1],
                      'total_conc' = out[,2]/vol_tot,
                      'epi_conc' = out[,3]/vol_epi,
                      'hypo_conc' = out[,4]/vol_hypo,
                      'total_mass' = out[,2],
                      'epi_mass' = out[,3],
                      'hypo_mass' = out[,4],
                      'total_vol' = vol_tot,
                      'epi_vol' = vol_epi,
                      'hypo_vol' = vol_hypo)
  ggplot(output)+
  geom_point(aes(time, total_conc, col = 'total')) +
  geom_point(aes(time, epi_conc, col = 'epi')) +
  geom_point(aes(time, hypo_conc, col = 'hypo')) +
    ylim(0,500) +
    xlim(0,50) +
    ylab('Conc in mg/m3')
  ggplot(output)+
    geom_point(aes(time, total_mass, col = 'total')) +
    geom_point(aes(time, epi_mass, col = 'epi')) +
    geom_point(aes(time, hypo_mass, col = 'hypo')) +
    # xlim(0,100) +
    ylab('Mass in mg')
  ggplot(output)+
    geom_point(aes(time, total_vol, col = 'total')) +
    geom_point(aes(time, epi_vol, col = 'epi')) +
    geom_point(aes(time, hypo_vol, col = 'hypo')) +
    ylab('Volume in m3')
  
  
  output <- read.table('/Users/robertladwig/Documents/DSI/odem/output.txt')
  output <- data.frame('fNEP' = output[,1], 
                       'fSED' = output[,2], 'fATM' = output[,3], 'fENTR1' = output[,4],
                       'fENTR2' = output[,5], 'flux_tot'= output[,6], 'flux_epi'= output[,7], 
                       'flux_hypo'= output[,8], 'timestep' = output[,9])
  
  ggplot(output) +
    geom_line(aes(timestep, fATM)) +
    xlim(0,5000) 
  
  ggplot(output) +
    geom_line(aes(timestep, flux_epi)) +
    xlim(0,5000) 
  
  ggplot(output) +
    geom_line(aes(timestep, fENTR1)) 
  
  melt_output = reshape2::melt(output, id.vars = 'timestep')
  
  ggplot(melt_output, aes(timestep, value)) +
    geom_point() +
    facet_wrap(~ variable, scales = 'free')
  
# }
