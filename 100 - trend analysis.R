setwd('...') # Your working directory
source('codes/000 - functions and paths.R')
path.result <- create_when_absent('results/100-trend analyses')
path.fit.anova <- str_c(path.result, '/fit_ANOVA') %>% create_when_absent()

data.marker <- read_csv(path.df.blood) %>% 
  select(id, age=age.blood, diagnosis, sex, education, APOE4, Aβ42:GFAP) %>% 
  mutate(APOE4 = case_when(
    APOE4 %in% c(24, 34, 44) ~ 'carrier',
    APOE4 %in% c(22, 23, 33) ~ 'non-carrier',
    TRUE ~ NA_character_
  )) %>% 
  pivot_longer(Aβ42:GFAP, names_to = 'name_response', values_to = 'response') %>% 
  mutate(diagnosis = factor(diagnosis, levels = c('NC', 'SCD', 'MCI', 'AD'))) %>% 
  mutate(region = -2, .after = id)

data.image <- read_csv(path.df.image) %>% 
  select(id, region, age=age.scan, diagnosis, sex, education, APOE4, name_response = measure, response = value) %>% 
  mutate(name_response = case_when(
    name_response == 'NE' ~ 'global efficiency',
    name_response == 'NLE' ~ 'generalised local efficiency',
    name_response == 'NLEgretna' ~ 'local efficiency',
    TRUE ~ NA_character_
  )) %>% 
  mutate(APOE4 = case_when(
    APOE4 %in% c(24, 34, 44) ~ 'carrier',
    APOE4 %in% c(22, 23, 33) ~ 'non-carrier',
    TRUE ~ NA_character_
  )) %>% 
  mutate(diagnosis = factor(diagnosis, levels = c('NC', 'SCD', 'MCI', 'AD')))

data.all <- rbind(data.marker, data.image) 

data.all.nest <- data.all %>% 
  group_by(region, name_response) %>% nest() 


# Run trend ANCOVA for each measures ----------------------------------------------

future::plan(future::multisession, workers = 8)   
ancov.coef <- furrr::future_pmap_dfr(
  list(data.all.nest$region, data.all.nest$name_response, data.all.nest$data),
  function(REGION, RESPONSE, tbl) {
    path.rds <- str_c(path.fit.anova, '/', str_replace(RESPONSE, '/', '-'), '_', REGION, '.rds')
    if (file.exists(path.rds)) {
      model <- readRDS(path.rds)
    } else {
      tbl <- tbl %>% na.omit() %>% 
        mutate(sex = factor(sex, levels = c("f", "m"), labels = c("Female", "Male"))) %>% 
        mutate(APOE4 = factor(APOE4, levels = c("non-carrier", "carrier")) ) 
      
      # Create a trend variable
      tbl$trend <- as.numeric(tbl$diagnosis) - 1
      
      # Create polynomial contrasts
      contrasts(tbl$diagnosis) <- contr.poly(4)
      model <- lm(response ~ diagnosis + age + sex + education + APOE4, data = tbl)
      saveRDS(model, path.rds)      
    }
    model.coef <- broom::tidy(model)
    sum.model <- summary(model)
    return(
      tibble(
        region = REGION,
        var = RESPONSE,
        r.adj = sum.model$adj.r.squared,
        p.value.L = model.coef %>% filter(term == 'diagnosis.L') %>% pull(p.value),
        t.L = model.coef %>% filter(term == 'diagnosis.L') %>% pull(statistic),
        p.value.Q = model.coef %>% filter(term == 'diagnosis.Q') %>% pull(p.value),
        t.Q = model.coef %>% filter(term == 'diagnosis.Q') %>% pull(statistic)
      )
    )
  }
)

# Adjust p-value for regional statistics -----------------------------
ancov.coef.global <- 
  ancov.coef %>% filter(region < 0) %>% 
  mutate(p.adj.L = p.value.L, p.adj.Q = p.value.Q) %>% 
  mutate(`full name` = var)

ancov.coef.regional <- 
  ancov.coef %>% filter(region > 0) %>% 
  group_by(var) %>% 
  mutate(p.adj.L = p.adjust(p.value.L, method = 'fdr')) %>% 
  mutate(p.adj.Q = p.adjust(p.value.Q, method = 'fdr')) %>% 
  left_join(
    bna.description %>% select(region, 'full name')
  ) %>% ungroup()

result.anconv <- 
  rbind(ancov.coef.global, ancov.coef.regional) %>% 
  mutate(report.L = str_c(
    't=', round(t.L, 3), ', ',
    'p FDR=', format(p.adj.L, scienticif = T, digits = 3)
  )) %>% 
  mutate(report.Q = str_c(
    't=', round(t.Q, 3), ', ',
    'p FDR=', format(p.adj.Q, scienticif = T, digits = 3)
  ))
write_csv(result.anconv, str_c(path.result, '/result_trend_report_all.csv'))
write_csv(result.anconv %>% dplyr::filter(p.adj.L < 0.05 | p.adj.Q < 0.05), str_c(path.result, '/result_trend_report_signif.csv'))


# Pairwise comparison  -----------------------------------
path.pairwise <- create_when_absent(str_c(path.result, '/pairwise'))
path.pairwise.comparison <-  create_when_absent(str_c(path.pairwise, '/comparison'))
trend.significant <- read_csv(str_c(path.result, '/result_trend_report_signif.csv'))
pwalk(
  list(trend.significant$region, trend.significant$var, trend.significant$report.L, trend.significant$report.Q, trend.significant$r.adj, trend.significant$`full name`), 
  function(REGION, VAR, REPORT.L, REPORT.Q, R.adj, FULLNAME) {
    model <- str_c(path.fit.anova, '/', str_replace(VAR, '/', '-'), '_', REGION, '.rds') %>% readRDS()
    tbl <- model$model
    # Pairwise comparisons
    df.residuals <- 
      tibble(
        residuals = resid(lm(response ~ age + sex + education + APOE4, data = tbl)),
        diagnosis = tbl$diagnosis
      ) %>% na.omit()
    
    pairwise.results <- 
      df.residuals %>% 
      rstatix::pairwise_t_test(residuals ~ diagnosis, pool.sd = FALSE, p.adjust.method = "bonferroni")
    write_csv(pairwise.results, str_c(path.pairwise.comparison, '/', str_replace(VAR, '/', '-'), '_', REGION, '.csv'))
  }
)