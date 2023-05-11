setwd('...')
source('codes/000 - functions and paths.R')
library(mediation)
REGION.COMPONENT = -3
NUM_BOOT_MEDIATION = 10000
path.result <- create_when_absent('results/300-mediation')

# Reads the data ----------------------------------------------------------
# Network efficiency
data.efficiency <- read_csv(path.df.cross.sectional) %>% 
  dplyr::select(id, region, name_efficiency = measure, efficiency = value) %>% 
  mutate(name_efficiency = ifelse(name_efficiency == 'NE', 'global efficiency', ifelse(name_efficiency == 'NLE', 'generalised local efficiency', 'local efficiency'))) %>% 
  rbind(
    read_csv('results/400-network PCA/component score.csv') %>% 
      mutate(
        id,
        region = REGION.COMPONENT,
        name_efficiency = 'generalised local efficiency', # it's only component
        efficiency = component,
        .keep = 'none'
      )
  ) %>% 
  # Add covariates
  left_join(
    data.efficiency <- read_csv(path.df.cross.sectional) %>% 
      dplyr::select(id, diagnosis = diagnosis.scan, age.scan, sex, education, APOE4) %>% 
      distinct(),
    by = 'id'
  ) 


# Cognition and blood markers
data.response <- read_csv(path.df.cross.sectional) %>% 
  dplyr::select(id, age.blood, diagnosis.blood, MMSE, MOCA = 'MOCA-B', Aβ42:GFAP) %>% 
  pivot_longer(MMSE:MOCA, names_to = 'name_cognition', values_to = 'cognition') %>% 
  pivot_longer(Aβ42:GFAP, names_to = 'name_marker', values_to = 'marker') %>% 
  distinct()

data.all <- data.efficiency %>% left_join(data.response) %>% 
  mutate(APOE4 = case_when(
    APOE4 %in% c(24, 34, 44) ~ 'carrier',
    APOE4 %in% c(22, 23, 33) ~ 'non-carrier',
    TRUE ~ NA_character_
  )) %>% mutate(APOE4 = factor(APOE4, levels = c('non-carrier', 'carrier'))) %>% 
  mutate(sex = factor(sex, levels = c('f', 'm'))) 

data.nest <- data.all %>% 
  filter(region<0) %>% 
  mutate(name_efficiency = if_else(region == REGION.COMPONENT, str_c('component of ', name_efficiency), name_efficiency)) %>% 
  group_by(name_efficiency, name_cognition, name_marker) %>% nest()



# Run mediation analyses on multiple regions ------------------------------
path.mediation.result <- str_c(path.result, '/RDS_mediation_MTC') %>% create_when_absent()

X.label = 'marker'
M.label = 'efficiency'
Y.label = 'cognition'
future::plan(future::multisession, workers = 8)   
result.mtc <- furrr::future_pmap_dfr(
  list(data.nest$name_cognition, data.nest$name_marker, data.nest$name_efficiency, data.nest$data), 
  function(COGNITION, MARKER, TOPOLOGY, df) {
    # print(df)
    df_scaled <- 
      df %>%
      mutate(age = (age.scan + age.blood) / 2) %>% 
      mutate(across(c(age, marker, efficiency, education, cognition), scale)) %>% 
      mutate(across(everything(), as.vector))
    
    # Fit the model for mediator and X
    model.mediator <- lm(efficiency ~ marker + age + sex + education + APOE4, data = df_scaled)
    coef.mediator <- broom::tidy(model.mediator)
    a <- coef.mediator %>% filter(term == X.label) %>% pull('estimate')
    a.se <- coef.mediator %>% filter(term == X.label) %>% pull('std.error')
    a.pval <- coef.mediator %>% filter(term == X.label) %>% pull('p.value')
    
    # Fit regression model for the total effect (Y, X and M)
    model.total <- lm(cognition ~ marker + efficiency + age + sex + education + APOE4, data = df_scaled)
    coef.total <- broom::tidy(model.total)
    b <- coef.total %>% filter(term == M.label) %>% pull('estimate')
    b.se <-  coef.total %>% filter(term == M.label) %>% pull('std.error')
    b.pval <-  coef.total %>% filter(term == M.label) %>% pull('p.value')
    
    # Conduct mediation analysis
    mediate.result <- mediation::mediate(model.mediator, model.total,  treat = X.label, mediator = M.label, sims = NUM_BOOT_MEDIATION, boot = T)
    mediate.sum <- summary(mediate.result)
    saveRDS(mediate.result, str_c(path.mediation.result, '/', COGNITION, '_', str_replace(MARKER, '/', '_'), '_', TOPOLOGY, '.rds'))
    
    # ACME stands for average causal mediation effects. This is the indirect effect of the X on the Y that goes through the M.
    acme <- mediate.sum$d0
    acme.pval <- mediate.sum$d0.p
    acme.low <- mediate.sum$d0.ci[[1]]
    acme.up <- mediate.sum$d0.ci[[2]]
    
    # ADE stands for average direct effects. It describes the direct effect of the IV on the DV.
    ade <- mediate.sum$z0
    ade.pval <- mediate.sum$z0.p
    ade.low <- mediate.sum$z0.ci[[1]]
    ade.up <- mediate.sum$z0.ci[[2]]
    
    # Total Effect stands for the total effect (direct + indirect) of the IV onto the DV.
    total <- mediate.sum$tau.coef
    total.pval <- mediate.sum$tau.p
    total.low <- mediate.sum$tau.ci[[1]]
    total.up <- mediate.sum$tau.ci[[2]]
    
    # Prop. Mediated describes the proportion of the effect of the IV on the DV that goes through the mediator.
    
    # Save the result of mediation
    return(tibble(
      cognition = COGNITION, marker = MARKER, efficiency = TOPOLOGY,
      a, a.se, a.pval,
      b, b.se, b.pval,
      acme, acme.pval, acme.low, acme.up,
      ade, ade.pval, ade.low, ade.up,
      total, total.pval, total.low, total.up
    ))
  }
)

pval.star <- function(p.value) {
  case_when(
    p.value < 0.0001 ~ '****',
    p.value < 0.001 ~ '***',
    p.value < 0.01 ~ '**',
    p.value < 0.05 ~ '*',
    p.value < 0.1 ~ '.',
    TRUE ~ 'NS'
  )
}


result <- result.mtc %>% mutate(
  report = str_c(
    'a=', round(a, 3), pval.star(a.pval), '(', round(a.low, 3), ',', round(a.up), ')\n', 
    'b=', round(b, 3), pval.star(b.pval), '(', round(b.low, 3), ',', round(b.up), ')\n',
    #'b=', round(b, 3), pval.star(b.pval), '(', round(b.se, 3), ')\n', 
    #'c\'=', round(ade, 3), pval.star(ade.pval), '(', round(ade.se, 3), ')\n', 
    'c\'=', round(ade, 3), pval.star(ade.pval), '(', round(ade.low, 3), ',', round(ade.up), ')\n',
    'ab=', round(acme, 3), pval.star(acme.pval), '(', round(acme.low, 3), ',', round(acme.up), ')'
  )
)
write_csv(result, str_c(path.result, '/mediation_MTC.csv'))
write_csv(result %>% filter(acme.pval<=0.05), str_c(path.result, '/mediation_significant_acme_MTC.csv'))

