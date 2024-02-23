# Library -----------------------------------------------------------------
library( tidyverse)
library( ARTofR)
library( hrbrthemes)
library( viridis)
library( sjlabelled)
library( lubridate)
# Cohort construction -----------------------------------------------------
# comtemporary cohort
comtemporary_comparison_cohort <- comtemporary_func() %>% mutate( infection_status = case_when( BT_index %in% c( "BT_infection", "infection") ~ "Infected", TRUE ~ "Uninfected"),
                                                                  Infection_settings = case_when( is.na( Infection_settings) ~ "Uninfected", TRUE ~ Infection_settings))

#historical cohort construction
historical_comparison_cohort <- historical_2017_2019_func() %>% mutate( infection_status = case_when( BT_index %in% c( "BT_infection", "infection") ~ "Infected", TRUE ~ "Uninfected"),
                                                                   Infection_settings = case_when( is.na( Infection_settings) ~ "Uninfected", TRUE ~ Infection_settings))
#no breakthrough cohort
no_breakthrough_comtemporary_cohort <- no_breakthrough_comtemporary_func() %>% mutate( infection_status = case_when( BT_index %in% c( "BT_infection", "infection") ~ "Infected", TRUE ~ "Uninfected"),
                                                                  Infection_settings = case_when( is.na( Infection_settings) ~ "Uninfected", TRUE ~ Infection_settings))
#breakthough cohort
breakthrough_comtemporary_cohort <- breakthrough_comtemporary_func() %>% mutate( infection_status = case_when( BT_index %in% c( "BT_infection", "infection") ~ "Infected", TRUE ~ "Uninfected"),
                                                                                       Infection_settings = case_when( is.na( Infection_settings) ~ "Uninfected", TRUE ~ Infection_settings))
#community setting cohort
community_comtemporary_cohort <- community_comtemporary_func() %>% mutate( infection_status = case_when( BT_index %in% c( "BT_infection", "infection") ~ "Infected", TRUE ~ "Uninfected"),
                                                                           Infection_settings = case_when( is.na( Infection_settings) ~ "Uninfected", TRUE ~ Infection_settings))

# hospital setting cohort
hospital_comtemporary_cohort <- hospital_comtemporary_func() %>% mutate( infection_status = case_when( BT_index %in% c( "BT_infection", "infection") ~ "Infected", TRUE ~ "Uninfected"),
                                                                           Infection_settings = case_when( is.na( Infection_settings) ~ "Uninfected", TRUE ~ Infection_settings))

# Baseline covariates and control for confounding (weighting) -----------------------------------------------------
# link other covariates dataframe
covariate_generate_func <- function( cohort){
  
  start.time <- Sys.time()
  # link to GP medication within 1 years before index date
  GP_medication_long <- 
    cohort %>% 
    select( eid, index_date) %>% 
    left_join( combined_gp_medication, by = "eid") %>%
    filter( issue_date < index_date , issue_date >= index_date - 365) %>%  
    left_join( all_medication_dmd_code, by = "dmd_code") %>% 
    left_join( medication_dictionary_dpa, by = c("top_level" = "Top_level_medication_classes")) 
  
  GP_medication_wide <- 
    GP_medication_long %>% 
    group_by( eid, Manual_category) %>% 
    filter( row_number() == 1) %>% 
    ungroup() %>% 
    mutate( medication = "1") %>% 
    pivot_wider( id_cols = eid,names_from = Manual_category, values_from = medication) %>% 
    janitor::clean_names()%>% 
    ungroup()
  
  #link to emis diagnosis (all history)
  GP_emis_diagnosis_long <- 
    cohort %>% 
    select( eid, index_date) %>% 
    left_join( emis_gp_clinical_extraction, by = "eid") %>% 
    filter( event_dt < index_date ) %>%  
    left_join( all_clinical_emis_code, by = c( "code" = "SnomedCTConceptId")) %>% 
    select( eid, second_level)
  
  #link to TPP diagnosis (all history)
  GP_tpp_diagnosis_long <- 
    cohort %>% 
    select( eid, index_date) %>% 
    left_join( tpp_gp_clinical_extraction, by = "eid") %>% 
    filter( event_dt < index_date ) %>%  
    left_join( all_clinical_tpp_code, by = c( "code" = "readcode_new")) %>% 
    select( eid, second_level)
  
  combined_diagnosis_wide <- 
    bind_rows( GP_emis_diagnosis_long, GP_tpp_diagnosis_long) %>% 
    group_by( eid, second_level) %>% 
    filter( row_number() == 1) %>% 
    ungroup() %>% 
    mutate( diagnosis = "1") %>% 
    pivot_wider( id_cols = eid, names_from = second_level, values_from = diagnosis) %>% 
    janitor::clean_names()%>% 
    ungroup()
  
  #link to hes phecode within 3 years before index date
  previous_clinical_sequelae <- 
    cohort %>% 
    select( eid, index_date) %>% 
    left_join( hes_phecode, by = "eid") %>% 
    filter( admidate < index_date, admidate >= index_date - 365 * 3) %>% 
    group_by( eid, phe_code) %>% 
    summarise( n = n()) %>% 
    ungroup()
  
  phecode_covariates_list <- 
    previous_clinical_sequelae %>% 
    group_by( phe_code ) %>% 
    summarise( persons_count = n_distinct( eid)) %>% 
    mutate( perc = persons_count/ length(unique( cohort$eid)) * 100) %>% 
    filter( perc > 0.5)
    
  previous_clinical_sequelae_wide <- 
    previous_clinical_sequelae %>% 
    mutate( outcome_index = "1") %>% 
    pivot_wider( id_cols = eid, names_from = phe_code, values_from = outcome_index) %>%   # attention: count frequency
    select( eid, all_of(phecode_covariates_list$phe_code)) %>% 
    ungroup()
  
  #link to hes number within 1 years before index date
  hes_admission_numb <- 
    cohort %>% 
    select( eid, index_date) %>% 
    left_join( hesin, by = "eid") %>% 
    filter( admidate < index_date, admidate >= index_date - 365) %>%
    group_by( eid) %>% 
    summarise( hes_count  = n_distinct( admidate)) %>% 
    ungroup()
  
  
  #link to vaccination
  GP_vaccination <- 
    cohort %>% 
    select( eid, index_date) %>% 
    left_join( vaccine_diagnosis_source, by = "eid") %>% 
    filter( event_dt < index_date) %>%
    group_by( eid) %>% 
    summarise( vaccination_count  = n_distinct( event_dt))
  
  all_output_df <- 
    cohort %>% 
    left_join( baseline_df, by = c("eid")) %>% #baseline charateristics
    left_join( select( baseline_df_addtional, eid, Diastolic_BP, Systolic_BP, smoking_status, alcohol_status, PA_status), by = "eid") %>% 
    left_join( GP_medication_wide, by = "eid")  %>% 
    left_join( combined_diagnosis_wide, by = "eid") %>% 
    left_join( previous_clinical_sequelae_wide, by = "eid") %>% 
    left_join( GP_vaccination, by = "eid") %>% 
    left_join( hes_admission_numb, by = "eid") %>%
    mutate( age_year = year( index_date) - birth_year) %>% 
    mutate( across( setdiff( c( names(GP_medication_wide), 
                                names(combined_diagnosis_wide),
                                phecode_covariates_list$phe_code), 
                             "eid"), ~ case_when( is.na(.) ~ "0", TRUE ~ .))) %>% 
    mutate( hes_count = case_when( is.na( hes_count) ~ 0, TRUE ~ as.numeric(hes_count)),
            vaccination_binary = case_when( is.na( vaccination_count) ~ "Not or partially vaccinated",  TRUE ~ "Fully vaccinated"),
            IMD = case_when( is.na( IMD) ~ mean( IMD, na.rm = TRUE), TRUE ~ IMD),
            bmi = case_when( is.na( bmi) ~ mean( bmi, na.rm = TRUE), TRUE ~ bmi),
            ethnicity = replace_na( ethnicity, "NA"),
            ethnicity = case_when( ethnicity %in% c( "White") ~ "White",
                                   TRUE ~ "Other ethnic groups")) %>% 
    mutate( calendar_scale_year = as.character( year( index_date)),
            calendar_scale_week = as.character( isoweek( index_date)),
            calendar_scale_week_num = isoweek( index_date)) 
  
  print( Sys.time() - start.time)
  
  return( list( GP_medication_wide = GP_medication_wide, 
                combined_diagnosis_wide = combined_diagnosis_wide,
                phecode_covariates_list = phecode_covariates_list,
                previous_clinical_sequelae_wide = previous_clinical_sequelae_wide,
                hes_admission_numb = hes_admission_numb,
                GP_vaccination = GP_vaccination,
                all_output_df = all_output_df))

  
  
}
comtemporary_comparison_cohort_covariates_df_list <- covariate_generate_func( cohort = comtemporary_comparison_cohort)
historical_comparison_cohort_covariates_df_list <- covariate_generate_func( cohort = historical_comparison_cohort)
no_breakthrough_comparison_cohort_covariates_df_list <- covariate_generate_func( cohort = no_breakthrough_comtemporary_cohort)
breakthrough_comparison_cohort_covariates_df_list <- covariate_generate_func( cohort = breakthrough_comtemporary_cohort)
community_comtemporary_cohort_df_list <- covariate_generate_func( cohort = community_comtemporary_cohort)
hospital_comtemporary_cohort_df_list <- covariate_generate_func( cohort = hospital_comtemporary_cohort)
test_negative_comparison_cohort_covariates_df_list <- covariate_generate_func( cohort = test_negative_comparison_cohort)
# baseline charateristics for different cohorts
minimal_covariate_list <- 
  c("age_year", "sex", "ethnicity", "IMD", "bmi",
    "Diastolic_BP", "Systolic_BP", "smoking_status", "alcohol_status", "PA_status",
    "vaccination_binary",
    "lipid_lowering_drugs",
    "ra_sinh",
    "antihypertensives",
    "anticoagulants",
    "antithrombotics",
    "proton_pump_inhibitors",
    "antidiabetes",
    "antidepressants",
    "systemic_glucocorticoids",
    "immunosuppressants_excl_corticosteroids",
    "antineoplastic_agents",
    "ch_cancer",
    "ch_cancer_meta",
    "ch_diabetes",
    "ch_diabetes_comp",
    "ch_heart_cong",
    "ch_mi",
    "ch_cerebrovasc",
    "ch_peripheralvasc",
    "ch_liver_mild",
    "ch_liver_mod",
    "ch_pulmonary",
    "ch_renal",
    "ch_pepticulcer",
    "ch_rheumatol",
    "ch_dementia",
    "ch_hemiplegia",
    "ch_aids",
    "hes_count"
    #"calendar_scale_week", "calendar_scale_year"
    )

maximal_covariate_list <- c( minimal_covariate_list, finalised_Covariate_phecode_list)

additional_covariate_list_for_maximal <- setdiff( Covariate_phecode_list, maximal_covariate_list)

breakthrough_maximal_covariate_list <- setdiff(  maximal_covariate_list, "vaccination_binary")

# historical only
historical_minimal_covariate_list <- 
  c("age_year", "sex", "ethnicity", "IMD", "bmi",
    "Diastolic_BP", "Systolic_BP", "smoking_status", "alcohol_status", "PA_status",
    #"vaccination_binary",
    "lipid_lowering_drugs",
    "ra_sinh",
    "antihypertensives",
    "anticoagulants",
    "antithrombotics",
    "proton_pump_inhibitors",
    "antidiabetes",
    "antidepressants",
    "systemic_glucocorticoids",
    "immunosuppressants_excl_corticosteroids",
    "antineoplastic_agents",
    "ch_cancer",
    "ch_cancer_meta",
    "ch_diabetes",
    "ch_diabetes_comp",
    "ch_heart_cong",
    "ch_mi",
    "ch_cerebrovasc",
    "ch_peripheralvasc",
    "ch_liver_mild",
    "ch_liver_mod",
    "ch_pulmonary",
    "ch_renal",
    "ch_pepticulcer",
    "ch_rheumatol",
    "ch_dementia",
    "ch_hemiplegia",
    "ch_aids",
    "hes_count",
    "hes_count_after",
    "GP_count_after" 
    #"calendar_scale_week", "calendar_scale_year"
  )

historical_finalised_Covariate_phecode_list <- 
  historical_comparison_cohort_covariates_df_list$phecode_covariates_list %>% 
  filter( perc >= 1) %>% 
  pull( phe_code)


historical_maximal_covariate_list <- c( historical_minimal_covariate_list, historical_finalised_Covariate_phecode_list)

# weighting
library( survey)
library( tableone)

Exact_Weighting_process_func <- function( input_data , 
                                          proportion , 
                                          model_covariate , 
                                          exposure_value,
                                          control_value) {
  
  set.seed(1)
  
  sample <- sample_frac( input_data, 1)
  
  print(paste("Number of missing value in all covariates are",sum(summarise( sample, across( all_of( model_covariate), ~ sum( is.na(.)))))))

  sample <- 
    sample%>% 
    mutate( exposure = case_when( .data[["infection_status"]] == exposure_value ~ 1, .data[["infection_status"]] == control_value ~ 0)) %>% 
    filter( !is.na(exposure))
  
  print( table(sample$exposure))
  
  PS_model <- glm( reformulate( termlabels = model_covariate , response = "exposure"), family  = binomial(link = "logit"), data = sample)
  
  # Weighted data with probability
  PS_probability <- mutate( sample, 
                            P_treatment = PS_model$fitted.values, 
                            P_control = 1 - PS_model$fitted.values,
                            PS_weights = case_when( exposure == 1 ~ 1/P_treatment, TRUE ~ 1/P_control))
  
  PS_probability_trim <-  
    PS_probability %>% 
    mutate( P_treatment = case_when( P_treatment >= quantile( P_treatment, 0.99) ~ quantile( P_treatment, 0.99),
                                     P_treatment <= quantile( P_treatment, 0.01) ~ quantile( P_treatment, 0.01)))
    #filter( P_treatment <= quantile( P_treatment, 0.99), P_treatment >= quantile( P_treatment, 0.01)) 

  
  ## Weighted data format
  #PS_IPTW <- survey::svydesign( ids = ~ 1, data = PS_probability, weights = ~ PS_weights)
  #PS_IPTW_trim <- survey::svydesign( ids = ~ 1, data = PS_probability_trim, weights = ~ PS_weights)
  
  return( list( PS_probability = PS_probability, 
                PS_probability_trim = PS_probability_trim
                #PS_IPTW = PS_IPTW,
                #PS_IPTW_trim = PS_IPTW_trim
                ))
  
}


Infection_comtemporary_HES_weighted <- Exact_Weighting_process_func( input_data = output_df, 
                                                                     proportion = 1, 
                                                                     model_covariate = historical_minimal_covariate_list, 
                                                                     exposure_value = "Infected",
                                                                     control_value = "Uninfected")

Infection_comtemporary_GP_weighted <- Exact_Weighting_process_func( input_data = linked_comtemporary_GP_cohort, 
                                                                     proportion = 1, 
                                                                     model_covariate = minimal_covariate_list, 
                                                                     exposure_value = "Infected",
                                                                     control_value = "Uninfected")

# Wrap up all primary analyses pipeline --------------------------------------------------------------------
HES_outcome_name_list <- unique( tidy_target_list$Outcome_name) %>% na.omit() %>%  set_names()
GP_outcome_name_list <- unique( psy_code_all_code_longer$Outcome_level) %>% na.omit() %>%  set_names()

# define function to construct cohorts
Combine_HES_GP_weighted_cohort_construction <- function( df,history, look_window, follow_up_window, covariates_df, covariate_list){
  
  # history outcome HES diagnosis
  all_outcome_HES <- 
    hesin_diag %>% 
    filter( diag_icd10 %in% tidy_target_list$raw_codes) %>% 
    left_join( select( hesin, eid, ins_index, admidate), by = c("eid", "ins_index")) %>% 
    select( eid, diag_icd10, admidate) %>% 
    rename( codes = diag_icd10)
  
  # history outcome GP medication
  all_outcome_GP <- 
    combined_covid19_emis_gp_scripts %>% 
    filter( dmd_code  %in% psy_code_all_code_longer$Code_value) %>% 
    select( eid, dmd_code, issue_date) %>% 
    rename( codes = dmd_code, admidate = issue_date)
  
  all_outcome <- bind_rows( all_outcome_HES, all_outcome_GP)
  
  history_outcome_wide <- 
    df %>% 
    select( eid, index_date) %>% 
    left_join( all_outcome, by = c( "eid")) %>% 
    filter( !is.na( codes), admidate < index_date, admidate >= index_date - look_window) %>% 
    group_by( eid) %>% 
    filter( row_number() == 1) %>% 
    ungroup() %>% 
    mutate( outcome_history_index = "1") %>% 
    janitor::clean_names()%>% 
    select( eid, outcome_history_index) %>% 
    ungroup()
  
  # incident outcome HES diagnosis
  incident_outcome_wide <- 
    df %>% 
    select( eid, index_date) %>% 
    left_join( all_outcome, by = c( "eid")) %>% 
    filter( admidate > index_date) %>% 
    arrange( admidate) %>% 
    group_by( eid) %>% 
    filter( row_number() == 1) %>% 
    ungroup() %>% 
    select( eid, admidate)
  
  incident_outcome_follow_up <- 
    df %>% 
    left_join( incident_outcome_wide, by = "eid") %>% 
    left_join( select( death, eid, date_of_death), by = "eid") %>% 
    mutate( ending_date = as.Date( "2021-09-30")) %>% 
    mutate( follow_up_end_date = case_when( is.na( date_of_death) ~ pmin( index_date + follow_up_window, ending_date), 
                                            TRUE ~ pmin( date_of_death, index_date + follow_up_window, ending_date)),
            incident_outcome = case_when( admidate >= index_date & admidate <= follow_up_end_date ~ 1, TRUE ~ 0),
            follow_up_days = as.numeric( pmin( follow_up_end_date, admidate, na.rm = TRUE) - index_date)) 
  
  # link covariates and historical outcome
  output_df <- 
    df %>% 
    left_join( history_outcome_wide, by = "eid")  %>% 
    left_join( select( incident_outcome_follow_up, eid, incident_outcome, follow_up_days), by = "eid") %>% 
    left_join( select(covariates_df$all_output_df, eid, all_of(covariate_list)), by = "eid") %>%   
    mutate( outcome_history_index = case_when( is.na(outcome_history_index) ~ "0", TRUE ~ outcome_history_index)) %>% 
    filter( outcome_history_index %in% history) 
  
  Infection_comtemporary_weighted <- Exact_Weighting_process_func( input_data = output_df, 
                                                                   proportion = 1, 
                                                                   model_covariate = covariate_list, 
                                                                   exposure_value = "Infected",
                                                                   control_value = "Uninfected")

  
  return( Infection_comtemporary_weighted)
  
  
  
}
# attention: HES outcome needs addtional censors
HES_weighted_cohort_construction <- function( df, target_outcome, history, look_window, follow_up_window, covariates_df, covariate_list){
  
  # history outcome HES diagnosis
  specific_outcome_code_list <- 
    filter( tidy_target_list, Outcome_name %in% target_outcome) %>% 
    pull( raw_codes)
  
  all_outcome_HES <- 
    hesin_diag %>% 
    filter( diag_icd10 %in% specific_outcome_code_list) 
  
  history_outcome_HES_wide <- 
    df %>% 
    select( eid, index_date) %>% 
    left_join( select( hesin, eid, ins_index, admidate), by = "eid")%>% 
    left_join( all_outcome_HES, by = c( "eid", "ins_index")) %>% 
    filter( !is.na( diag_icd10), admidate < index_date, admidate >= index_date - look_window) %>% 
    group_by( eid) %>% 
    filter( row_number() == 1) %>% 
    ungroup() %>% 
    mutate( outcome_history_index = "1") %>% 
    janitor::clean_names()%>% 
    select( eid, outcome_history_index) %>% 
    ungroup()
  
  # incident outcome HES diagnosis
  incident_outcome_HES_wide <- 
    df %>% 
    select( eid, index_date) %>% 
    left_join( select( hesin, eid, ins_index, admidate), by = "eid")%>% 
    left_join( all_outcome_HES, by = c( "eid", "ins_index")) %>% 
    filter( !is.na( diag_icd10), admidate > index_date) %>% 
    arrange( admidate) %>% 
    group_by( eid) %>% 
    filter( row_number() == 1) %>% 
    ungroup() %>% 
    select( eid, admidate)

  
  incident_outcome_HES_follow_up <- 
    df %>% 
    left_join( incident_outcome_HES_wide, by = "eid") %>% 
    left_join( select( death, eid, date_of_death), by = "eid") %>% 
    left_join( Infection_update, by = "eid") %>% 
    mutate( follow_up_end_date = case_when( is.na( date_of_death) ~ pmin( index_date + follow_up_window, as.Date( ending_date), infection_update_end_date, na.rm = TRUE), 
                                            TRUE ~ pmin( date_of_death, index_date + follow_up_window, as.Date( ending_date), infection_update_end_date, na.rm = TRUE)),
            incident_outcome = case_when( admidate >= index_date & admidate <= follow_up_end_date ~ 1, TRUE ~ 0),
            follow_up_days = as.numeric( pmin( follow_up_end_date, admidate, na.rm = TRUE) - index_date)) 
  
  
  # link covariates and historical outcome
  output_df <- 
    df %>% 
    left_join( history_outcome_HES_wide, by = "eid")  %>% 
    left_join( select( incident_outcome_HES_follow_up, eid, incident_outcome, follow_up_days), by = "eid") %>% 
    left_join( select(covariates_df$all_output_df, eid, all_of(covariate_list)), by = "eid") %>%   
    mutate( outcome_history_index = case_when( is.na(outcome_history_index) ~ "0", TRUE ~ outcome_history_index)) %>% 
    filter( outcome_history_index %in% history) 
  
  Infection_comtemporary_HES_weighted <- Exact_Weighting_process_func( input_data = output_df, 
                                                                       proportion = 1, 
                                                                       model_covariate = covariate_list, 
                                                                       exposure_value = "Infected",
                                                                       control_value = "Uninfected")

  
  return( Infection_comtemporary_HES_weighted)
  
  
}
GP_weighted_cohort_construction <- function( df, Outcomes , history, look_window, follow_up_window, covariates_df, covariate_list){
  
  # history outcome GP medication
  specific_outcome_code_list <- 
    filter( psy_code_all_code_longer, Outcome_level == Outcomes) %>% 
    pull( Code_value)
  
  all_outcome_GP <- 
    combined_covid19_emis_gp_scripts %>% 
    filter( dmd_code  %in% specific_outcome_code_list) 
  
  history_outcome_GP_wide <- 
    df %>% 
    select( eid, index_date) %>% 
    left_join( all_outcome_GP, by = c( "eid")) %>% 
    filter(  issue_date < index_date, issue_date >= index_date - look_window) %>% 
    group_by( eid) %>% 
    filter( row_number() == 1) %>% 
    ungroup() %>% 
    mutate( outcome_history_index = "1") %>% 
    janitor::clean_names()%>% 
    select( eid, outcome_history_index) %>% 
    ungroup()
  
  # incident outcome HES diagnosis
  incident_outcome_GP_wide <- 
    df %>% 
    select( eid, index_date) %>% 
    left_join( all_outcome_GP, by = c( "eid")) %>% 
    filter(  issue_date >= index_date) %>% 
    arrange( issue_date) %>% 
    group_by( eid) %>% 
    filter( row_number() == 1) %>% 
    ungroup() %>% 
    select( eid, issue_date)
  
  incident_outcome_GP_follow_up <- 
    df %>% 
    left_join( incident_outcome_GP_wide, by = "eid") %>% 
    left_join( select( death, eid, date_of_death), by = "eid") %>% 
    mutate( ending_date = as.Date( "2021-09-30")) %>% 
    mutate( follow_up_end_date = case_when( is.na( date_of_death) ~ pmin( index_date + follow_up_window, as.Date( ending_date)), 
                                            TRUE ~ pmin( date_of_death, index_date + follow_up_window, as.Date( ending_date))),
            incident_outcome = case_when( issue_date >= index_date & issue_date <= follow_up_end_date ~ 1, TRUE ~ 0),
            follow_up_days = as.numeric( pmin( follow_up_end_date, issue_date, na.rm = TRUE) - index_date))
  
  
  # link covariates and historical outcome
  output_df <- 
    df %>% 
    left_join( history_outcome_GP_wide, by = "eid")  %>% 
    left_join( select( incident_outcome_GP_follow_up, eid, incident_outcome, follow_up_days), by = "eid") %>% 
    left_join( select(covariates_df$all_output_df, eid, all_of(covariate_list)), by = "eid") %>%   
    mutate( outcome_history_index = case_when( is.na(outcome_history_index) ~ "0", TRUE ~ outcome_history_index)) %>% 
    filter( outcome_history_index %in% history) 
  
  Infection_comtemporary_GP_weighted <- Exact_Weighting_process_func( input_data = output_df, 
                                                                       proportion = 1, 
                                                                       model_covariate = covariate_list, 
                                                                       exposure_value = "Infected",
                                                                       control_value = "Uninfected")
  
  return( Infection_comtemporary_GP_weighted)
  
  
}

# define function to preduce HR
Survival_weighted_Model_func <- function( df = HES_weighted_cohort_df_list$Psychiatric_diagnoses){
  
  temp_df <- df$PS_probability_trim 
  
  formula_input <- as.formula( paste( paste( "survival::Surv( follow_up_days, incident_outcome )", "~"), paste( c("exposure"), collapse="+")))
  
  output_HR <- survival::coxph( formula = formula_input, data = temp_df, weights = PS_weights) %>% temp_func()
  
  output_summary_statistics <- 
    temp_df %>% 
    group_by( exposure) %>% 
    summarise( sum_follow = sum( as.numeric(follow_up_days * PS_weights)),
               median_follow = median( follow_up_days),
               #total_number = n(),
               cases = sum( incident_outcome * PS_weights),
               #proportion = cases/total_number * 1000,
               rate = cases/ sum_follow * 1000  * 365) %>% 
    rowwise() %>% 
    mutate( epiR::epi.conf( as.matrix( cbind(.data[["cases"]], .data[["sum_follow"]])), ctype = "inc.rate", method = "exact", N = 1000, design = 1, conf.level = 0.95)* 1000 * 365) %>% 
    select( - est) %>% 
    pivot_wider( names_from = exposure, values_from = c( median_follow, sum_follow, cases, rate, lower, upper)) %>% 
    mutate( ARD_lower = rate_0 * (output_HR$lower_95 - 1),
            ARD = rate_0 * (output_HR$exp_coef - 1),
            ARD_higer = rate_0 * (output_HR$upper_95 - 1)) %>% 
    mutate( ARD_CI = paste(format(round(ARD, 2), nsmall = 2), " (",
                           format(round(ARD_lower, 2), nsmall = 2), " to ",
                           format(round(ARD_higer, 2), nsmall = 2), ")",sep = ""))
  
  output <- bind_cols( output_HR, output_summary_statistics)
  
  return( output)
  
  
}
Survival_weighted_setting_Model_func <- function( df){
  
  temp_df <- df$PS_probability_trim
  
  temp_df <- temp_df %>% mutate( Infection_settings = factor( Infection_settings, levels = c( "Uninfected", "community_infection", "hospital_infection")))
  
  formula_input <- as.formula( paste( paste( "survival::Surv( follow_up_days, incident_outcome )", "~"), paste( c("Infection_settings"), collapse="+")))
  
  output_HR <- survival::coxph( formula = formula_input, data = temp_df, weights = PS_weights) %>% temp_func() %>% 
    select( vars, exp_coef, lower_95, upper_95) %>% 
    pivot_wider( names_from = vars, values_from = c( exp_coef, lower_95, upper_95))
  
  output_summary_statistics <- 
    temp_df %>% 
    group_by( Infection_settings) %>% 
    summarise( sum_follow = sum( as.numeric(follow_up_days * PS_weights) ),
               median_follow = median( follow_up_days),
               #total_number = n(),
               cases = sum( incident_outcome * PS_weights),
               #proportion = cases/total_number * 1000,
               rate = cases/ sum_follow * 1000  * 365) %>% 
    rowwise() %>% 
    mutate( epiR::epi.conf( as.matrix( cbind(.data[["cases"]], .data[["sum_follow"]])), ctype = "inc.rate", method = "exact", N = 1000, design = 1, conf.level = 0.95)* 1000 * 365) %>% 
    select( - est) %>%  
    pivot_wider( names_from = Infection_settings, values_from = c( median_follow, sum_follow, cases, rate, lower, upper)) %>% 
    mutate( community_infection_ARD_lower = rate_Uninfected * (output_HR$lower_95_Infection_settingscommunity_infection - 1),
            community_infection_ARD = rate_Uninfected * (output_HR$exp_coef_Infection_settingscommunity_infection - 1),
            community_infection_ARD_higer = rate_Uninfected * (output_HR$upper_95_Infection_settingscommunity_infection - 1)) %>% 
    mutate( hospital_infection_ARD_lower = rate_Uninfected * (output_HR$lower_95_Infection_settingshospital_infection - 1),
            hospital_infection_ARD = rate_Uninfected * (output_HR$exp_coef_Infection_settingshospital_infection - 1),
            hospital_infection_ARD_higer = rate_Uninfected * (output_HR$upper_95_Infection_settingshospital_infection - 1))
  
  
  output <- bind_cols( output_HR, output_summary_statistics)
  
  return( output)
  
  
}
Survival_weighted_breakthrough_Model_func <- function( df){
  
  temp_df <- df$PS_probability_trim
  
  temp_df <- temp_df %>% mutate( infection_type = factor( BT_index, levels = c( "comtemporary_control", "infection", "BT_infection")))
  
  formula_input <- as.formula( paste( paste( "survival::Surv( follow_up_days, incident_outcome )", "~"), paste( c("infection_type"), collapse="+")))
  
  output_HR <- survival::coxph( formula = formula_input, data = temp_df, weights = PS_weights) %>% temp_func() %>% 
    select( vars, exp_coef, lower_95, upper_95) %>% 
    pivot_wider( names_from = vars, values_from = c( exp_coef, lower_95, upper_95))
  
  output_summary_statistics <- 
    temp_df %>% 
    group_by( infection_type) %>% 
    summarise( sum_follow = sum( as.numeric(follow_up_days * PS_weights) ),
               median_follow = median( follow_up_days),
               #total_number = n(),
               cases = sum( incident_outcome * PS_weights),
               #proportion = cases/total_number * 1000,
               rate = cases/ sum_follow * 1000  * 365) %>% 
    rowwise() %>% 
    mutate( epiR::epi.conf( as.matrix( cbind(.data[["cases"]], .data[["sum_follow"]])), ctype = "inc.rate", method = "exact", N = 1000, design = 1, conf.level = 0.95)* 1000 * 365) %>% 
    select( - est) %>% 
    pivot_wider( names_from = infection_type, values_from = c( median_follow, sum_follow, cases, rate, lower, upper)) %>% 
    mutate( normal_infection_ARD_lower = rate_comtemporary_control * (output_HR$lower_95_infection_typeinfection - 1),
            normal_infection_ARD = rate_comtemporary_control * (output_HR$exp_coef_infection_typeinfection - 1),
            normal_infection_ARD_higer = rate_comtemporary_control * (output_HR$upper_95_infection_typeinfection - 1)) %>% 
    mutate( BT_infection_ARD_lower = rate_comtemporary_control * (output_HR$lower_95_infection_typeBT_infection - 1),
            BT_infection_ARD = rate_comtemporary_control * (output_HR$exp_coef_infection_typeBT_infection - 1),
            BT_infection_ARD_higer = rate_comtemporary_control * (output_HR$upper_95_infection_typeBT_infection - 1))
  
  
  output <- bind_cols( output_HR, output_summary_statistics)
  
  return( output)
  
}

#produce weighted incident cohort
Combine_HES_GP_weighted_cohort_df_list <- Combine_HES_GP_weighted_cohort_construction( df = comtemporary_comparison_cohort, 
                                                                                       history = "0",
                                                                                       look_window = 365*2, 
                                                                                       follow_up_window = 365,
                                                                                       covariates_df = comtemporary_comparison_cohort_covariates_df_list,
                                                                                       covariate_list = maximal_covariate_list)

HES_weighted_cohort_df_list <- map( HES_outcome_name_list[1] , 
                                    HES_weighted_cohort_construction, 
                                    df = comtemporary_comparison_cohort, 
                                    look_window = 365*2, 
                                    follow_up_window = 365,
                                    history = "0",
                                    covariates_df = comtemporary_comparison_cohort_covariates_df_list,
                                    covariate_list = maximal_covariate_list)

GP_weighted_cohort_df_list <- map( GP_outcome_name_list[2], 
                                   GP_weighted_cohort_construction, 
                                   df = comtemporary_comparison_cohort, 
                                   look_window = 365*2, 
                                   follow_up_window = 365,
                                   history = "0",
                                   covariates_df = comtemporary_comparison_cohort_covariates_df_list,
                                   covariate_list = maximal_covariate_list)


#produce weighted prevalent cohort
Prevalent_Combine_HES_GP_weighted_cohort_df_list <- Combine_HES_GP_weighted_cohort_construction( df = comtemporary_comparison_cohort, 
                                                                                                 history = c( "0", "1"),
                                                                                                 look_window = 365*2, 
                                                                                                 follow_up_window = 365,
                                                                                                 covariates_df = comtemporary_comparison_cohort_covariates_df_list,
                                                                                                 covariate_list = maximal_covariate_list)

Prevalent_HES_weighted_cohort_df_list <- map( HES_outcome_name_list[1], 
                                              HES_weighted_cohort_construction, 
                                              df = comtemporary_comparison_cohort, 
                                              look_window = 365*2, 
                                              follow_up_window = 365,
                                              history = c( "0", "1"),
                                              covariates_df = comtemporary_comparison_cohort_covariates_df_list,
                                              covariate_list = maximal_covariate_list)


Prevalent_GP_weighted_cohort_df_list <- map( GP_outcome_name_list[2], 
                                             GP_weighted_cohort_construction, 
                                             df = comtemporary_comparison_cohort, 
                                             look_window = 365*2, 
                                             follow_up_window = 365,
                                             history = c( "0", "1"),
                                             covariates_df = comtemporary_comparison_cohort_covariates_df_list,
                                             covariate_list = maximal_covariate_list)

# produce HR 
Combine_HES_GP_infection_HR <- Survival_weighted_Model_func( Combine_HES_GP_weighted_cohort_df_list) %>% mutate( study_outcome = "any related outcome", .before = vars)

HES_infection_HR <- 
  map_df( HES_weighted_cohort_df_list, 
          Survival_weighted_Model_func, 
          .id = "study_outcome")

GP_infection_HR <- 
  map_df( GP_weighted_cohort_df_list, 
          Survival_weighted_Model_func, 
          .id = "study_outcome")

Prevalent_Combine_HES_GP_infection_HR <- Survival_weighted_Model_func( Prevalent_Combine_HES_GP_weighted_cohort_df_list) %>% mutate( study_outcome = "any related outcome", .before = vars)

Prevalent_HES_infection_HR <- 
  map_df( Prevalent_HES_weighted_cohort_df_list, 
          Survival_weighted_Model_func, 
          .id = "study_outcome")

Prevalent_GP_infection_HR <- 
  map_df( Prevalent_GP_weighted_cohort_df_list, 
          Survival_weighted_Model_func, 
          .id = "study_outcome")

# Wrap up all primary analyses pipeline (historical control) --------------------------------------------------------------------
Hist_combine_HES_GP_weighted_cohort_construction <- function( df,history, look_window, follow_up_window, covariates_df, covariate_list){
  
  # history outcome HES diagnosis
  all_outcome_HES <- 
    hesin_diag %>% 
    filter( diag_icd10 %in% tidy_target_list$raw_codes) %>% 
    left_join( select( hesin, eid, ins_index, admidate), by = c("eid", "ins_index")) %>% 
    select( eid, diag_icd10, admidate) %>% 
    rename( codes = diag_icd10)
  
  # history outcome GP medication
  all_outcome_GP <- 
    combined_covid19_emis_gp_scripts %>% 
    filter( dmd_code  %in% psy_code_all_code_longer$Code_value) %>% 
    select( eid, dmd_code, issue_date) %>% 
    rename( codes = dmd_code, admidate = issue_date)
  
  all_outcome <- bind_rows( all_outcome_HES, all_outcome_GP)
  
  history_outcome_wide <- 
    df %>% 
    select( eid, index_date) %>% 
    left_join( all_outcome, by = c( "eid")) %>% 
    filter( !is.na( codes), admidate < index_date, admidate >= index_date - look_window) %>% 
    group_by( eid) %>% 
    filter( row_number() == 1) %>% 
    ungroup() %>% 
    mutate( outcome_history_index = "1") %>% 
    janitor::clean_names()%>% 
    select( eid, outcome_history_index) %>% 
    ungroup()
  
  # incident outcome HES diagnosis
  incident_outcome_wide <- 
    df %>% 
    select( eid, index_date) %>% 
    left_join( all_outcome, by = c( "eid")) %>% 
    filter( admidate > index_date) %>% 
    arrange( admidate) %>% 
    group_by( eid) %>% 
    filter( row_number() == 1) %>% 
    ungroup() %>% 
    select( eid, admidate)
  
  incident_outcome_follow_up <- 
    df %>% 
    left_join( incident_outcome_wide, by = "eid") %>% 
    mutate( ending_date = case_when( infection_status == "Infected" ~ as.Date( "2021-09-30"),
                                     TRUE ~ as.Date( "2018-09-30"))) %>% 
    left_join( select( death, eid, date_of_death), by = "eid") %>% 
    mutate( follow_up_end_date = case_when( is.na( date_of_death) ~ pmin( index_date + follow_up_window, ending_date), 
                                            TRUE ~ pmin( date_of_death, index_date + follow_up_window, ending_date)),
            incident_outcome = case_when( admidate >= index_date & admidate <= follow_up_end_date ~ 1, TRUE ~ 0),
            follow_up_days = as.numeric( pmin( follow_up_end_date, admidate, na.rm = TRUE) - index_date)) 
  
  # link covariates and historical outcome
  output_df <- 
    df %>% 
    select( -hes_count_after, -GP_count_after) %>% 
    left_join( history_outcome_wide, by = "eid")  %>% 
    left_join( select( incident_outcome_follow_up, eid, incident_outcome, follow_up_days), by = "eid") %>% 
    left_join( select(covariates_df$all_output_df, eid, all_of(covariate_list)), by = "eid") %>%   
    mutate( outcome_history_index = case_when( is.na(outcome_history_index) ~ "0", TRUE ~ outcome_history_index)) %>% 
    filter( outcome_history_index %in% history) 
  
  Infection_comtemporary_weighted <- Exact_Weighting_process_func( input_data = output_df, 
                                                                   proportion = 1, 
                                                                   model_covariate = covariate_list, 
                                                                   exposure_value = "Infected",
                                                                   control_value = "Uninfected")
  
  
  return( Infection_comtemporary_weighted)
  
  
  
}
Hist_HES_weighted_cohort_construction <- function( df, Outcomes = "Psychiatric_diagnoses", history, look_window, follow_up_window, covariates_df, covariate_list){
  
  # history outcome HES diagnosis
  specific_outcome_code_list <- 
    filter( tidy_target_list, Outcome_name == Outcomes) %>% 
    pull( raw_codes)
  
  all_outcome_HES <- 
    hesin_diag %>% 
    filter( diag_icd10 %in% specific_outcome_code_list) 
  
  history_outcome_HES_wide <- 
    df %>% 
    select( eid, index_date) %>% 
    left_join( select( hesin, eid, ins_index, admidate), by = "eid")%>% 
    left_join( all_outcome_HES, by = c( "eid", "ins_index")) %>% 
    filter( !is.na( diag_icd10), admidate < index_date, admidate >= index_date - look_window) %>% 
    group_by( eid) %>% 
    filter( row_number() == 1) %>% 
    ungroup() %>% 
    mutate( outcome_history_index = "1") %>% 
    janitor::clean_names()%>% 
    select( eid, outcome_history_index) %>% 
    ungroup()
  
  # incident outcome HES diagnosis
  incident_outcome_HES_wide <- 
    df %>% 
    select( eid, index_date) %>% 
    left_join( select( hesin, eid, ins_index, admidate), by = "eid")%>% 
    left_join( all_outcome_HES, by = c( "eid", "ins_index")) %>% 
    filter( !is.na( diag_icd10), admidate >= index_date) %>% 
    arrange( admidate) %>% 
    group_by( eid) %>% 
    filter( row_number() == 1) %>% 
    ungroup() %>% 
    select( eid, admidate)
  
  incident_outcome_HES_follow_up <- 
    df %>% 
    left_join( incident_outcome_HES_wide, by = "eid") %>% 
    left_join( select( death, eid, date_of_death), by = "eid") %>% 
    mutate( follow_up_end_date = case_when( is.na( date_of_death) ~ pmin( index_date + follow_up_window, as.Date( ending_date)), 
                                            TRUE ~ pmin( date_of_death, index_date + follow_up_window, as.Date( ending_date))),
            incident_outcome = case_when( admidate >= index_date & admidate <= follow_up_end_date ~ 1, TRUE ~ 0),
            follow_up_days = as.numeric( pmin( follow_up_end_date, admidate, na.rm = TRUE) - index_date)) 
  
  
  # link covariates and historical outcome
  output_df <- 
    df %>% 
    select( -hes_count_after, -GP_count_after) %>% 
    left_join( history_outcome_HES_wide, by = "eid")  %>% 
    left_join( select( incident_outcome_HES_follow_up, eid, incident_outcome, follow_up_days), by = "eid") %>% 
    left_join( select(covariates_df$all_output_df, eid, all_of(covariate_list)), by = "eid") %>%   
    mutate( outcome_history_index = case_when( is.na(outcome_history_index) ~ "0", TRUE ~ outcome_history_index)) %>% 
    filter( outcome_history_index %in% history) 
  
  
  Infection_comtemporary_HES_weighted <- Exact_Weighting_process_func( input_data = output_df, 
                                                                       proportion = 1, 
                                                                       model_covariate = covariate_list, 
                                                                       exposure_value = "Infected",
                                                                       control_value = "Uninfected")
  
  
  return( Infection_comtemporary_HES_weighted)
  
  
}
Hist_GP_weighted_cohort_construction <- function( df, Outcomes , history, look_window, follow_up_window, covariates_df, covariate_list){
  
  # history outcome GP medication
  specific_outcome_code_list <- 
    filter( psy_code_all_code_longer, Outcome_level == Outcomes) %>% 
    pull( Code_value)
  
  all_outcome_GP <- 
    combined_covid19_emis_gp_scripts %>% 
    filter( dmd_code  %in% specific_outcome_code_list) 
  
  history_outcome_GP_wide <- 
    df %>% 
    select( eid, index_date) %>% 
    left_join( all_outcome_GP, by = c( "eid")) %>% 
    filter(  issue_date < index_date, issue_date >= index_date - look_window) %>% 
    group_by( eid) %>% 
    filter( row_number() == 1) %>% 
    ungroup() %>% 
    mutate( outcome_history_index = "1") %>% 
    janitor::clean_names()%>% 
    select( eid, outcome_history_index) %>% 
    ungroup()
  
  # incident outcome HES diagnosis
  incident_outcome_GP_wide <- 
    df %>% 
    select( eid, index_date) %>% 
    left_join( all_outcome_GP, by = c( "eid")) %>% 
    filter(  issue_date >= index_date) %>% 
    arrange( issue_date) %>% 
    group_by( eid) %>% 
    filter( row_number() == 1) %>% 
    ungroup() %>% 
    select( eid, issue_date)
  
  incident_outcome_GP_follow_up <- 
    df %>% 
    left_join( incident_outcome_GP_wide, by = "eid") %>% 
    mutate( ending_date = case_when( infection_status == "Infected" ~ as.Date( "2021-09-30"),
                                     TRUE ~ as.Date( "2018-09-30"))) %>% 
    left_join( select( death, eid, date_of_death), by = "eid") %>% 
    mutate( follow_up_end_date = case_when( is.na( date_of_death) ~ pmin( index_date + follow_up_window, as.Date( ending_date)), 
                                            TRUE ~ pmin( date_of_death, index_date + follow_up_window, as.Date( ending_date))),
            incident_outcome = case_when( issue_date >= index_date & issue_date <= follow_up_end_date ~ 1, TRUE ~ 0),
            follow_up_days = as.numeric( pmin( follow_up_end_date, issue_date, na.rm = TRUE) - index_date))
  
  
  # link covariates and historical outcome
  output_df <- 
    df %>% 
    select( -hes_count_after, -GP_count_after) %>% 
    left_join( history_outcome_GP_wide, by = "eid")  %>% 
    left_join( select( incident_outcome_GP_follow_up, eid, incident_outcome, follow_up_days), by = "eid") %>% 
    left_join( select(covariates_df$all_output_df, eid, all_of(covariate_list)), by = "eid") %>%   
    mutate( outcome_history_index = case_when( is.na(outcome_history_index) ~ "0", TRUE ~ outcome_history_index)) %>% 
    filter( outcome_history_index %in% history) 
  
  Infection_comtemporary_GP_weighted <- Exact_Weighting_process_func( input_data = output_df, 
                                                                      proportion = 1, 
                                                                      model_covariate = covariate_list, 
                                                                      exposure_value = "Infected",
                                                                      control_value = "Uninfected")
  
  return( Infection_comtemporary_GP_weighted)
  
  
}

table( historical_comparison_cohort$infection_status)

#produce weighted incident cohort
Hist_Combine_HES_GP_weighted_cohort_df_list <- Hist_combine_HES_GP_weighted_cohort_construction( df = historical_comparison_cohort, 
                                                                                                 history = "0",
                                                                                                 look_window = 365*2, 
                                                                                                 follow_up_window = 365,
                                                                                                 covariates_df = historical_comparison_cohort_covariates_df_list,
                                                                                                 covariate_list = historical_minimal_covariate_list)

table( historical_comparison_cohort$GP_count_after)

Hist_HES_weighted_cohort_df_list <- map( HES_outcome_name_list[1], 
                                         Hist_HES_weighted_cohort_construction, 
                                         df = historical_comparison_cohort, 
                                         look_window = 365*2, 
                                         follow_up_window = 365,
                                         history = c("0"),
                                         covariates_df = historical_comparison_cohort_covariates_df_list,
                                         covariate_list = historical_minimal_covariate_list)


Hist_GP_weighted_cohort_df_list <- map( GP_outcome_name_list[2], 
                                        Hist_GP_weighted_cohort_construction, 
                                        df = historical_comparison_cohort, 
                                        look_window = 365*2, 
                                        follow_up_window = 365,
                                        history = c( "0"),
                                        covariates_df = historical_comparison_cohort_covariates_df_list,
                                        covariate_list = historical_minimal_covariate_list)


#produce weighted prevalent cohort
Hist_Prevalent_Combine_HES_GP_weighted_cohort_df_list <- Hist_combine_HES_GP_weighted_cohort_construction( df = historical_comparison_cohort, 
                                                                                                           history = c("0", "1"),
                                                                                                           look_window = 365*2, 
                                                                                                           follow_up_window = 365,
                                                                                                           covariates_df = historical_comparison_cohort_covariates_df_list,
                                                                                                           covariate_list = historical_minimal_covariate_list)

Hist_Prevalent_HES_weighted_cohort_df_list <- map( HES_outcome_name_list[1], 
                                                   Hist_HES_weighted_cohort_construction, 
                                                   df = historical_comparison_cohort, 
                                                   look_window = 365*2, 
                                                   follow_up_window = 365,
                                                   history = c("0","1"),
                                                   covariates_df = historical_comparison_cohort_covariates_df_list,
                                                   covariate_list = historical_minimal_covariate_list)


Hist_Prevalent_GP_weighted_cohort_df_list <- map( GP_outcome_name_list[2], 
                                                  Hist_GP_weighted_cohort_construction, 
                                                  df = historical_comparison_cohort, 
                                                  look_window = 365*2, 
                                                  follow_up_window = 365,
                                                  history = c( "0", "1"),
                                                  covariates_df = historical_comparison_cohort_covariates_df_list,
                                                  covariate_list = historical_minimal_covariate_list)


# produce HR 
Hist_Combine_HES_GP_infection_HR <- Survival_weighted_Model_func( Hist_Combine_HES_GP_weighted_cohort_df_list) %>% mutate( study_outcome = "any related outcome", .before = vars)

Hist_HES_infection_HR <- 
  map_df( Hist_HES_weighted_cohort_df_list, 
          Survival_weighted_Model_func, 
          .id = "study_outcome")

Hist_GP_infection_HR <- 
  map_df( Hist_GP_weighted_cohort_df_list, 
          Survival_weighted_Model_func, 
          .id = "study_outcome")

Hist_Prevalent_Combine_HES_GP_infection_HR <- Survival_weighted_Model_func( Hist_Prevalent_Combine_HES_GP_weighted_cohort_df_list) %>% mutate( study_outcome = "any related outcome", .before = vars)

Hist_Prevalent_HES_infection_HR <- 
  map_df( Hist_Prevalent_HES_weighted_cohort_df_list, 
          Survival_weighted_Model_func, 
          .id = "study_outcome")

Hist_Prevalent_GP_infection_HR <- 
  map_df( Hist_Prevalent_GP_weighted_cohort_df_list, 
          Survival_weighted_Model_func, 
          .id = "study_outcome")

# Wrap up all primary analyses pipeline for test negetive control ---------
TN_Combine_HES_GP_weighted_cohort_df_list <- Combine_HES_GP_weighted_cohort_construction( df = test_negative_comparison_cohort , 
                                                                                          history = "0",
                                                                                          look_window = 365*2, 
                                                                                          follow_up_window = 365,
                                                                                          covariates_df = test_negative_comparison_cohort_covariates_df_list ,
                                                                                          covariate_list = maximal_covariate_list)

TN_HES_weighted_cohort_df_list <- map(HES_outcome_name_list[1] , 
                                      HES_weighted_cohort_construction, 
                                      df = test_negative_comparison_cohort, 
                                      look_window = 365*2, 
                                      follow_up_window = 365,
                                      history = "0",
                                      covariates_df = test_negative_comparison_cohort_covariates_df_list,
                                      covariate_list = maximal_covariate_list)



TN_GP_weighted_cohort_df_list <- map( GP_outcome_name_list[2], 
                                      GP_weighted_cohort_construction, 
                                      df = test_negative_comparison_cohort, 
                                      look_window = 365*2, 
                                      follow_up_window = 365,
                                      history = "0",
                                      covariates_df = test_negative_comparison_cohort_covariates_df_list,
                                      covariate_list = maximal_covariate_list)



# produce HR 
TN_Combine_HES_GP_infection_HR <- Survival_weighted_Model_func( TN_Combine_HES_GP_weighted_cohort_df_list) %>% mutate( study_outcome = "any related outcome", .before = vars)

TN_HES_infection_HR <- 
  map_df( TN_HES_weighted_cohort_df_list, 
          Survival_weighted_Model_func, 
          .id = "study_outcome")

TN_GP_infection_HR <- 
  map_df( TN_GP_weighted_cohort_df_list, 
          Survival_weighted_Model_func, 
          .id = "study_outcome")

# Wrap up non-breakthrough and breakthrough pipeline --------------------------------------------------------------------
#produce weighted incident cohort
bt_func <- function(index_label){
  
  no_breakthrough_HES_GP_weighted_cohort_df_list <- Combine_HES_GP_weighted_cohort_construction( df = no_breakthrough_comtemporary_cohort, 
                                                                                                 history = index_label,
                                                                                                 look_window = 365*2, 
                                                                                                 follow_up_window = 365,
                                                                                                 covariates_df = no_breakthrough_comparison_cohort_covariates_df_list ,
                                                                                                 covariate_list = breakthrough_maximal_covariate_list)
  
  no_breakthrough_HES_weighted_cohort_df_list <- map( HES_outcome_name_list[1] , 
                                                      HES_weighted_cohort_construction, 
                                                      df = no_breakthrough_comtemporary_cohort, 
                                                      look_window = 365*2, 
                                                      follow_up_window = 365,
                                                      history = index_label,
                                                      covariates_df = no_breakthrough_comparison_cohort_covariates_df_list,
                                                      covariate_list = breakthrough_maximal_covariate_list)
  
  no_breakthrough_GP_weighted_cohort_df_list <- map( GP_outcome_name_list[2], 
                                                     GP_weighted_cohort_construction, 
                                                     df = no_breakthrough_comtemporary_cohort, 
                                                     look_window = 365*2, 
                                                     follow_up_window = 365,
                                                     history = index_label,
                                                     covariates_df = no_breakthrough_comparison_cohort_covariates_df_list,
                                                     covariate_list = breakthrough_maximal_covariate_list)
  
  
  
  breakthrough_HES_GP_weighted_cohort_df_list <- Combine_HES_GP_weighted_cohort_construction( df = breakthrough_comtemporary_cohort, 
                                                                                              history = index_label,
                                                                                              look_window = 365*2, 
                                                                                              follow_up_window = 365,
                                                                                              covariates_df = breakthrough_comparison_cohort_covariates_df_list ,
                                                                                              covariate_list = breakthrough_maximal_covariate_list)
  
  
  breakthrough_HES_weighted_cohort_df_list <- map( HES_outcome_name_list[1], 
                                                   HES_weighted_cohort_construction, 
                                                   df = breakthrough_comtemporary_cohort, 
                                                   look_window = 365*2, 
                                                   follow_up_window = 365,
                                                   history = index_label,
                                                   covariates_df = breakthrough_comparison_cohort_covariates_df_list,
                                                   covariate_list = breakthrough_maximal_covariate_list)
  
  
  breakthrough_GP_weighted_cohort_df_list <- map( GP_outcome_name_list[2], 
                                                  GP_weighted_cohort_construction, 
                                                  df = breakthrough_comtemporary_cohort, 
                                                  look_window = 365*2, 
                                                  follow_up_window = 365,
                                                  history = index_label,
                                                  covariates_df = breakthrough_comparison_cohort_covariates_df_list,
                                                  covariate_list = breakthrough_maximal_covariate_list)
  
  return( output = list( no_breakthrough_HES_GP_weighted_cohort_df_list = no_breakthrough_HES_GP_weighted_cohort_df_list,
                         no_breakthrough_HES_weighted_cohort_df_list = no_breakthrough_HES_weighted_cohort_df_list,
                         no_breakthrough_GP_weighted_cohort_df_list = no_breakthrough_GP_weighted_cohort_df_list,
                         breakthrough_HES_GP_weighted_cohort_df_list = breakthrough_HES_GP_weighted_cohort_df_list,
                         breakthrough_HES_weighted_cohort_df_list = breakthrough_HES_weighted_cohort_df_list,
                         breakthrough_GP_weighted_cohort_df_list = breakthrough_GP_weighted_cohort_df_list))
  
  
}

incident_bt <- bt_func( index_label = c("0"))
prevalent_bt <- bt_func( index_label = c("0", "1"))

# produce HR 
no_breakthrough_HES_infection_HR <- 
  map_df( incident_bt$no_breakthrough_HES_weighted_cohort_df_list, 
          Survival_weighted_Model_func, 
          .id = "study_outcome")

no_breakthrough_GP_infection_HR <- 
  map_df( incident_bt$no_breakthrough_GP_weighted_cohort_df_list, 
          Survival_weighted_Model_func, 
          .id = "study_outcome")

no_breakthrough_HES_GP_infection_HR <- Survival_weighted_Model_func( incident_bt$no_breakthrough_HES_GP_weighted_cohort_df_list) %>% mutate( study_outcome = "any related outcome", .before = vars)


breakthrough_HES_infection_HR <- 
  map_df( breakthrough_HES_weighted_cohort_df_list, 
          Survival_weighted_Model_func, 
          .id = "study_outcome")


breakthrough_GP_infection_HR <- 
  map_df( incident_bt$breakthrough_GP_weighted_cohort_df_list, 
          Survival_weighted_Model_func, 
          .id = "study_outcome")

breakthrough_HES_GP_infection_HR <- Survival_weighted_Model_func( incident_bt$breakthrough_HES_GP_weighted_cohort_df_list) %>% mutate( study_outcome = "any related outcome", .before = vars)

prevalent_no_breakthrough_HES_GP_infection_HR <- Survival_weighted_Model_func( prevalent_bt$no_breakthrough_HES_GP_weighted_cohort_df_list) %>% mutate( study_outcome = "any related outcome", .before = vars)

prevalent_no_breakthrough_HES_infection_HR <- 
  map_df( prevalent_bt$no_breakthrough_HES_weighted_cohort_df_list, 
          Survival_weighted_Model_func, 
          .id = "study_outcome")

prevalent_no_breakthrough_GP_infection_HR <- 
  map_df( prevalent_bt$no_breakthrough_GP_weighted_cohort_df_list, 
          Survival_weighted_Model_func, 
          .id = "study_outcome")

prevalent_breakthrough_HES_GP_infection_HR <- Survival_weighted_Model_func( prevalent_bt$breakthrough_HES_GP_weighted_cohort_df_list) %>% mutate( study_outcome = "any related outcome", .before = vars)

prevalent_breakthrough_HES_infection_HR <- 
  map_df( prevalent_bt$breakthrough_HES_weighted_cohort_df_list, 
          Survival_weighted_Model_func, 
          .id = "study_outcome")

prevalent_breakthrough_GP_infection_HR <- 
  map_df( prevalent_bt$breakthrough_GP_weighted_cohort_df_list, 
          Survival_weighted_Model_func, 
          .id = "study_outcome")

# Wrap up community and hospital pipeline --------------------------------------------------------------------
#produce weighted incident cohort
community_func <- function(index_label = "0"){
  
  community_HES_GP_weighted_cohort_df_list <- Combine_HES_GP_weighted_cohort_construction( df = community_comtemporary_cohort, 
                                                                                                 history = index_label,
                                                                                                 look_window = 365*2, 
                                                                                                 follow_up_window = 365,
                                                                                                 covariates_df = community_comtemporary_cohort_df_list ,
                                                                                                 covariate_list = maximal_covariate_list)
  
  community_HES_weighted_cohort_df_list <- map( HES_outcome_name_list[1] , 
                                                HES_weighted_cohort_construction, 
                                                df = community_comtemporary_cohort, 
                                                look_window = 365*2, 
                                                follow_up_window = 365,
                                                history = index_label,
                                                covariates_df = community_comtemporary_cohort_df_list,
                                                covariate_list = maximal_covariate_list)
  
  community_GP_weighted_cohort_df_list <- map( GP_outcome_name_list[2], 
                                                     GP_weighted_cohort_construction, 
                                                     df = community_comtemporary_cohort, 
                                                     look_window = 365*2, 
                                                     follow_up_window = 365,
                                                     history = index_label,
                                                     covariates_df = community_comtemporary_cohort_df_list,
                                                     covariate_list = maximal_covariate_list)
  
  
  
  hospital_HES_GP_weighted_cohort_df_list <- Combine_HES_GP_weighted_cohort_construction( df = hospital_comtemporary_cohort, 
                                                                                              history = index_label,
                                                                                              look_window = 365*2, 
                                                                                              follow_up_window = 365,
                                                                                              covariates_df = hospital_comtemporary_cohort_df_list ,
                                                                                              covariate_list = maximal_covariate_list)
  
  
  hospital_HES_weighted_cohort_df_list <- map( HES_outcome_name_list[1] , 
                                                   HES_weighted_cohort_construction, 
                                                   df = hospital_comtemporary_cohort, 
                                                   look_window = 365*2, 
                                                   follow_up_window = 365,
                                                   history = index_label,
                                                   covariates_df = hospital_comtemporary_cohort_df_list,
                                                   covariate_list = maximal_covariate_list)
  
  
  hospital_GP_weighted_cohort_df_list <- map( GP_outcome_name_list[2], 
                                                  GP_weighted_cohort_construction, 
                                                  df = hospital_comtemporary_cohort, 
                                                  look_window = 365*2, 
                                                  follow_up_window = 365,
                                                  history = index_label,
                                                  covariates_df = hospital_comtemporary_cohort_df_list,
                                                  covariate_list = breakthrough_maximal_covariate_list)
  
  return( output = list( community_HES_GP_weighted_cohort_df_list = community_HES_GP_weighted_cohort_df_list,
                         community_HES_weighted_cohort_df_list = community_HES_weighted_cohort_df_list,
                         community_GP_weighted_cohort_df_list = community_GP_weighted_cohort_df_list,
                         hospital_HES_GP_weighted_cohort_df_list = hospital_HES_GP_weighted_cohort_df_list,
                         hospital_HES_weighted_cohort_df_list = hospital_HES_weighted_cohort_df_list,
                         hospital_GP_weighted_cohort_df_list = hospital_GP_weighted_cohort_df_list))
  
  
}
incident_community <- community_func( index_label = c("0"))
prevalent_community <- community_func( index_label = c("0", "1"))

# produce HR 
community_HES_GP_infection_HR <- Survival_weighted_Model_func( incident_community$community_HES_GP_weighted_cohort_df_list) %>% mutate( study_outcome = "any related outcome", .before = vars)

community_HES_infection_HR <- 
  map_df( incident_community$community_HES_weighted_cohort_df_list, 
          Survival_weighted_Model_func, 
          .id = "study_outcome")

community_GP_infection_HR <- 
  map_df( incident_community$community_GP_weighted_cohort_df_list, 
          Survival_weighted_Model_func, 
          .id = "study_outcome")

hospital_HES_GP_infection_HR <- Survival_weighted_Model_func( incident_community$hospital_HES_GP_weighted_cohort_df_list) %>% mutate( study_outcome = "any related outcome", .before = vars)

hospital_HES_infection_HR <- 
  map_df( incident_community$hospital_HES_weighted_cohort_df_list, 
          Survival_weighted_Model_func, 
          .id = "study_outcome")

hospital_GP_infection_HR <- 
  map_df( incident_community$hospital_GP_weighted_cohort_df_list, 
          Survival_weighted_Model_func, 
          .id = "study_outcome")

write.csv( community_HES_GP_infection_HR, "D:/Document_sync/OneDrive - Nexus365/DPhil_training/Collaboration/Yunhe Wang/Covid_Mental_health/submission/Revision for NHB/Results/community_HES_GP_infection_HR.csv", row.names = TRUE)
write.csv( community_HES_infection_HR, "D:/Document_sync/OneDrive - Nexus365/DPhil_training/Collaboration/Yunhe Wang/Covid_Mental_health/submission/Revision for NHB/Results/community_HES_infection_HR.csv", row.names = TRUE)
write.csv( community_GP_infection_HR, "D:/Document_sync/OneDrive - Nexus365/DPhil_training/Collaboration/Yunhe Wang/Covid_Mental_health/submission/Revision for NHB/Results/community_GP_infection_HR.csv", row.names = TRUE)

write.csv( hospital_HES_GP_infection_HR, "D:/Document_sync/OneDrive - Nexus365/DPhil_training/Collaboration/Yunhe Wang/Covid_Mental_health/submission/Revision for NHB/Results/hospital_HES_GP_infection_HR.csv", row.names = TRUE)
write.csv( hospital_HES_infection_HR, "D:/Document_sync/OneDrive - Nexus365/DPhil_training/Collaboration/Yunhe Wang/Covid_Mental_health/submission/Revision for NHB/Results/hospital_HES_infection_HR.csv", row.names = TRUE)
write.csv( hospital_GP_infection_HR, "D:/Document_sync/OneDrive - Nexus365/DPhil_training/Collaboration/Yunhe Wang/Covid_Mental_health/submission/Revision for NHB/Results/hospital_GP_infection_HR.csv", row.names = TRUE)

prevalent_community_HES_GP_infection_HR <- Survival_weighted_Model_func( prevalent_community$community_HES_GP_weighted_cohort_df_list) %>% mutate( study_outcome = "any related outcome", .before = vars)

prevalent_community_HES_infection_HR <- 
  map_df( prevalent_community$community_HES_weighted_cohort_df_list, 
          Survival_weighted_Model_func, 
          .id = "study_outcome")

prevalent_community_GP_infection_HR <- 
  map_df( prevalent_community$community_GP_weighted_cohort_df_list, 
          Survival_weighted_Model_func, 
          .id = "study_outcome")

prevalent_hospital_HES_GP_infection_HR <- Survival_weighted_Model_func( prevalent_community$hospital_HES_GP_weighted_cohort_df_list) %>% mutate( study_outcome = "any related outcome", .before = vars)

prevalent_hospital_HES_infection_HR <- 
  map_df( prevalent_community$hospital_HES_weighted_cohort_df_list, 
          Survival_weighted_Model_func, 
          .id = "study_outcome")

prevalent_hospital_GP_infection_HR <- 
  map_df( prevalent_community$hospital_GP_weighted_cohort_df_list, 
          Survival_weighted_Model_func, 
          .id = "study_outcome")

write.csv( prevalent_community_HES_GP_infection_HR, "D:/Document_sync/OneDrive - Nexus365/DPhil_training/Collaboration/Yunhe Wang/Covid_Mental_health/submission/Revision for NHB/Results/prevalent_community_HES_GP_infection_HR.csv", row.names = TRUE)
write.csv( prevalent_community_HES_infection_HR, "D:/Document_sync/OneDrive - Nexus365/DPhil_training/Collaboration/Yunhe Wang/Covid_Mental_health/submission/Revision for NHB/Results/prevalent_community_HES_infection_HR.csv", row.names = TRUE)
write.csv( prevalent_community_GP_infection_HR, "D:/Document_sync/OneDrive - Nexus365/DPhil_training/Collaboration/Yunhe Wang/Covid_Mental_health/submission/Revision for NHB/Results/prevalent_community_GP_infection_HR.csv", row.names = TRUE)

write.csv( prevalent_hospital_HES_GP_infection_HR, "D:/Document_sync/OneDrive - Nexus365/DPhil_training/Collaboration/Yunhe Wang/Covid_Mental_health/submission/Revision for NHB/Results/prevalent_hospital_HES_GP_infection_HR.csv", row.names = TRUE)
write.csv( prevalent_hospital_HES_infection_HR, "D:/Document_sync/OneDrive - Nexus365/DPhil_training/Collaboration/Yunhe Wang/Covid_Mental_health/submission/Revision for NHB/Results/prevalent_hospital_HES_infection_HR.csv", row.names = TRUE)
write.csv( prevalent_hospital_GP_infection_HR, "D:/Document_sync/OneDrive - Nexus365/DPhil_training/Collaboration/Yunhe Wang/Covid_Mental_health/submission/Revision for NHB/Results/prevalent_hospital_GP_infection_HR.csv", row.names = TRUE)


