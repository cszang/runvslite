if (!require(VSLiteR)) {
  devtools::install_github("suztolwinskiward/VSLiteR")
}
if (!require(DEoptim)) {
  install.packages("DEoptim")
}
if (!require(dplyr)) {
  install.packages("dplyr")
}
if (!require(magrittr)) {
  install.packages("magrittr")
}

matlab_month_format <- function(.climate, .var) {
  ..climate <- .climate %>% dplyr::select(year, month, !!.var)
  t(as.matrix(tidyr::spread(..climate, month, !!.var) %>% dplyr::select(-year)))
}

make_vsinput_historic <- function(.rwl, .climate, restrict = NULL) {
  vars <- names(.climate)[!names(.climate) %in% c("year", "month")]
  nvars <- length(vars)
  .chron <- .rwl %>%
    detrend(method = "Spline", nyrs = 32) %>% 
    chron() %>% 
    na.omit()
  .chron <- data.frame(
    year = as.numeric(rownames(.chron)),
    rwi = .chron$std
  )
  tree_years <- .chron$year
  climate_years <- unique(.climate$year)
  common_years <- inner_join(data.frame(year = tree_years),
                               data.frame(year = climate_years)) %>% .$year
  if (!is.null(restrict)) {
    restrict_years <- min(restrict):max(restrict)
    if (all(restrict_years %in% common_years)) {
      common_years <- restrict_years  
    } else {
      stop("Please specify `restrict` as integer vector of length 2 giving range of years for calibration.")
    }
  }
  .chron <- .chron %>% filter(year %in% common_years)
  .climate <- .climate %>% filter(year %in% common_years)
  out <- list()
  out$trw <- .chron$rwi
  out$syear <- min(common_years)
  out$eyear <- max(common_years)
  for (i in 1:nvars) {
    .var <- vars[i]
    clim_var <- matlab_month_format(.climate, .var)
    out[[.var]] <- clim_var
  }
  out
}

make_vsinput_transient <- function(.climate) {
  vars <- names(.climate)[!names(.climate) %in% c("year", "month", "rcp")]
  nvars <- length(vars)
  out <- list()
  out$syear <- min(.climate$year)
  out$eyear <- max(.climate$year)
  for (i in 1:nvars) {
    .var <- vars[i]
    clim_var <- matlab_month_format(.climate, .var)
    out[[.var]] <- clim_var
  }
  out
}

vs_params <- function(trw, .temp, .prec, .syear, .eyear, .phi, iter = 200) {
  diffcor_optim <- function(x) {
    vs <- VSLite(syear = .syear,
                 eyear = .eyear,
                 phi = .phi,
                 T = .temp,
                 P = .prec,
                 T1 = x[1],
                 T2 = x[2],
                 M1 = x[3],
                 M2 = x[4])
    model <- t(vs$trw)
    corr <- 1 - cor(model, trw)
    if (is.na(corr)) corr <- Inf
    corr
  }
  lower <- c(0, 9, 0.01, 0.1)
  upper <- c(8.5, 20, 0.03, 0.5)
  opti <- DEoptim(diffcor_optim, lower, upper,
                  control = DEoptim.control(itermax = iter))
  out <- opti$optim$bestmem
  names(out) <- c("T1", "T2", "M1", "M2")
  out
}

vs_run_forward <- function(.vs_params, .temp, .prec, .syear, .eyear, .phi, 
                           return_pretty = TRUE) {
  out <- VSLite(syear = .syear, eyear = .eyear, phi = .phi, T = .temp, P = .prec,
                T1 = .vs_params["T1"], T2 = .vs_params["T2"], 
                M1 = .vs_params["M1"], M2 = .vs_params["M2"])
  if (return_pretty) {
    data.frame(year = .syear:.eyear,
               trw = t(out$trw))
  }  else out
}
