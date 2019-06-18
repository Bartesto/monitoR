# precursor code for WIN data clean up for annual reporting

library(tidyverse)
library(lubridate)
library(readxl)
library(janitor)

# vars - to be deleted
path <- "../nutrient_monitoring/swan_nutrient_WIN"
reportingYear <- 2019


cleanR <- function(path, reportingYear = 2019, export = TRUE){
  #setup info
  #these are the fields to import - need verification
  fields <- c("collection_method_name", "sample_type_code", "collect_date",
              "area_ref", "alkalinity_tot_ca_co3_mg_l",
              "c_sol_org_doc_doc_as_npoc_mg_l", "chlorophyll_a_by_vol_mg_l",
              "n_sum_sol_org_don_mg_l", "n_sum_sol_ox_n_ox_n_ton_mg_l",
              "n_tot_tn_p_tn_mg_l", "nh3_n_nh4_n_sol_mg_l", "o_do_in_situ_mg_l",
              "p_tot_tp_p_tp_mg_l", "po4_p_sol_react_srp_frp_mg_l",
              "salinity_mg_l", "secchi_depth_m", "si_o2_si_sol_react_mg_l",
              "tss_mg_l", "temperature_in_situ_deg_c", "p_h_no_units")

  #function to replace all '<' values LOR/2
  conv_less <- function(x){
    ifelse(str_detect(x, "^<"), as.numeric(sapply(str_split(x, "<"),
                                                  function(x) x[2]))/2,
           as.numeric(x))
  }

  report_start <- paste0(reportingYear - 6, "-06-01")
  report_fin <- paste0(reportingYear - 1, "-06-01")

  files <- list.files(path = path, pattern = ".xlsx",
                          full.names = TRUE)

  #do one to add to (handling missing rows)
  # temp_read <- readxl::read_excel(files[1])
  #
  # skip_rows <- NULL
  # search_string <- "Site Ref"
  # max_rows_to_search <- 10
  #
  # while (length(skip_rows) == 0) {
  #   skip_rows <- which(stringr::str_detect(temp_read[1:max_rows_to_search,][[1]],search_string)) - 0
  # }
  #
  # d1 <- readxl::read_excel(files[1], skip = skip_rows, .name_repair = ~ janitor::make_clean_names)

  tidy_all <- d1 %>%
    select(fields) %>%
    mutate(chlorophyll_a_by_vol_mg_l = conv_less(chlorophyll_a_by_vol_mg_l),
           n_sum_sol_ox_n_ox_n_ton_mg_l = conv_less(n_sum_sol_ox_n_ox_n_ton_mg_l),
           nh3_n_nh4_n_sol_mg_l = conv_less(nh3_n_nh4_n_sol_mg_l),
           po4_p_sol_react_srp_frp_mg_l = conv_less(po4_p_sol_react_srp_frp_mg_l))

  #control for only 1 file or append all data to first
  if(length(files) > 1){
    for(i in seq_along(files)[-1]){
      d <- readxl::read_excel(files[i], .name_repair = ~ janitor::make_clean_names)
      t_d <- d %>%
        select(fields) %>%
        mutate(chlorophyll_a_by_vol_mg_l = conv_less(chlorophyll_a_by_vol_mg_l),
               n_sum_sol_ox_n_ox_n_ton_mg_l = conv_less(n_sum_sol_ox_n_ox_n_ton_mg_l),
               nh3_n_nh4_n_sol_mg_l = conv_less(nh3_n_nh4_n_sol_mg_l),
               po4_p_sol_react_srp_frp_mg_l = conv_less(po4_p_sol_react_srp_frp_mg_l))
      tidy_all <- dplyr::bind_rows(tidy_all, t_d)
    }
  } else {
    tidy_all <- tidy_all
  }

  #add reporting zones plot order for months
  tidy_all <- tidy_all %>%
    mutate(emz = case_when(area_ref == "BLA" | area_ref == "ARM" | area_ref == "HEA" | area_ref == "NAR" ~ "lower",
                           area_ref == "NIL" | area_ref == "STJ" | area_ref == "MAY" | area_ref == "RON" ~ "middle",
                           area_ref == "KIN" | area_ref == "SUC" | area_ref == "WMP" | area_ref == "MSB" ~ "upper",
                           TRUE ~ "nrz")) %>%
    filter(collect_date >= report_start & collect_date < report_fin ) %>%
    mutate(mth = month(collect_date),
           pord = case_when(mth == 1 ~ 8,
                            mth == 2 ~ 9,
                            mth == 3 ~ 10,
                            mth == 4 ~ 11,
                            mth == 5 ~ 12,
                            mth == 6 ~ 1,
                            mth == 7 ~ 2,
                            mth == 8 ~ 3,
                            mth == 9 ~ 4,
                            mth == 10 ~ 5,
                            mth == 11 ~ 6,
                            mth == 12 ~ 7))

  #setup datasets for different variables
  #1 surface and bottom grab samples
  seven_vars <- c("n_tot_tn_p_tn_mg_l", "nh3_n_nh4_n_sol_mg_l",
                  "n_sum_sol_ox_n_ox_n_ton_mg_l", "n_sum_sol_org_don_mg_l",
                  "p_tot_tp_p_tp_mg_l", "po4_p_sol_react_srp_frp_mg_l",
                  "si_o2_si_sol_react_mg_l")

  seven <- tidy_all %>%
    filter(emz != "nrz" & collection_method_name == "Grab sample" & sample_type_code != "ST_STAND") %>%
    gather(xvar, value, seven_vars) %>%
    group_by(emz, pord, sample_type_code, xvar) %>%
    summarise(lower = quantile(value, probs = 0.10, na.rm = TRUE),
              upper = quantile(value, probs = 0.90, na.rm = TRUE),
              median = median(value, na.rm = TRUE))

  #2 surface and bottom insitu samples
  four_vars <- c("o_do_in_situ_mg_l", "salinity_mg_l",
                 "temperature_in_situ_deg_c", "p_h_no_units")

  four <- tidy_all %>%
    filter(emz != "nrz" & collection_method_name == "Insitu" & sample_type_code != "ST_STAND") %>%
    gather(xvar, value, four_vars) %>%
    group_by(emz, pord, sample_type_code, xvar) %>%
    summarise(lower = quantile(value, probs = 0.10, na.rm = TRUE),
              upper = quantile(value, probs = 0.90, na.rm = TRUE),
              median = median(value, na.rm = TRUE))

  #3 surface only grab samples
  three_vars <- c("tss_mg_l", "alkalinity_tot_ca_co3_mg_l",
             "c_sol_org_doc_doc_as_npoc_mg_l")

  three <- tidy_all %>%
    filter(emz != "nrz" & collection_method_name == "Grab sample" & sample_type_code == "ST_SURF") %>%
    gather(xvar, value, three_vars) %>%
    group_by(emz, pord, sample_type_code, xvar) %>%
    summarise(lower = quantile(value, probs = 0.10, na.rm = TRUE),
              upper = quantile(value, probs = 0.90, na.rm = TRUE),
              median = median(value, na.rm = TRUE))

  #4 chlorophyl a - surface and integrated
  chla <- tidy_all %>%
    filter(emz != "nrz" & collection_method_name != "Insitu" & sample_type_code != "ST_BOTT") %>%
    gather(xvar, value, chlorophyll_a_by_vol_mg_l) %>%
    group_by(emz, pord, sample_type_code, xvar) %>%
    summarise(lower = quantile(value, probs = 0.10, na.rm = TRUE),
              upper = quantile(value, probs = 0.90, na.rm = TRUE),
              median = median(value, na.rm = TRUE))

  #5 secchi is what it is
  secchi_1 <- tidy_all %>%
    filter(emz != "nrz") %>%
    select(emz, pord, collect_date, sample_type_code,chlorophyll_a_by_vol_mg_l) %>%
    gather(xvar, value, chlorophyll_a_by_vol_mg_l) %>%
    group_by(emz, pord, sample_type_code, xvar) %>%
    summarise(lower = quantile(value, probs = 0.10, na.rm = TRUE),
              upper = quantile(value, probs = 0.90, na.rm = TRUE),
              median = median(value, na.rm = TRUE))


}

tidy_all %>%
  select(emz, pord, area_ref, collect_date, sample_type_code,secchi_depth_m) %>%
  gather(xvar, value, secchi_depth_m) %>%
  drop_na(value)
