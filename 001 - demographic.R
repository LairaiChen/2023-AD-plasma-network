setwd('...') # Your working directory
source('codes/000 - functions and paths.R')
path.result <- str_c('results/100-demographic/') %>% create_when_absent()

# Related functions ---------
pvalue <- function(x, ...) {
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x) # data with group info
  g <- factor(rep(1:length(x), times=sapply(x, length))) # e.g. 11222
  if (is.numeric(y)) {
    # For numeric variables, perform a standard 2-sample t-test or ANOVA
    if (length(unique(g)) == 2) {
      p <- t.test(y ~ g)$p.value
    } else {
      p <- broom::tidy(aov(y ~ g))$p.value[[1]]
    }
  } else {
    # For categorical variables, perform a chi-squared test of independence
    p <- chisq.test(table(y, g))$p.value
  }
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
}

statistic.and.pval <- function(x, ...) {
  y <- unlist(x) # data with group info
  g <- factor(rep(1:length(x), times=sapply(x, length))) # e.g. 11222
  if (is.numeric(y)) {
    # For numeric variables, perform a standard 2-sample t-test or ANOVA
    if (length(unique(g)) == 2) {
      temp <- t.test(y ~ t)
      p <- temp$p.value
      statistic <- temp$statistic %>% unname()
    } else {
      temp <- broom::tidy(aov(y ~ g))
      p <- temp$p.value[[1]]
      statistic <- temp$statistic[[1]]
    }
  } else {
    # For categorical variables, perform a chi-squared test of independence
    temp <- chisq.test(table(y, g))
    p <- temp$p.value
    statistic <- temp$statistic %>% unname()
  }
  # Format the p-value, using an HTML entity for the less-than sign.
  str_c(round(statistic, 3), '(', case_when(
    p < .0001 ~ '****',
    p < .001 ~ '***',
    p < .01 ~ '**',
    p < .05 ~ '*',
    p < .1 ~ '&',
    TRUE ~ 'N.S'
  ), ')')
}
# Demographic of blood set ------------------------------------------------------------

df <- read_csv(path.df.blood) %>%
  rename(age = age.blood) %>% 
  mutate(diagnosis = factor(diagnosis, levels = c('NC', 'SCD', 'MCI', 'AD'))) %>% 
  mutate(APOE4 = case_when(
    APOE4 %in% c(24, 34, 44) ~ 'carrier',
    APOE4 %in% c(22, 23, 33) ~ 'non-carrier',
    TRUE ~ NA_character_
  ))

table1::table1(
  ~sex+age+sex+education+MMSE+`MOCA-B`+APOE4+Aβ42+Aβ40+pTau181+NfL+GFAP+`pTau/Aβ42`+`Aβ42/Aβ40`|diagnosis, 
  data = df,
  extra.col = list(`Statistics`=statistic.and.pval)
  ) %>% 
  as.data.frame(demo.table) %>% 
  write_csv(str_c(path.result, '/demo_blood.csv'))

  
# Demographic of image set -------
df <- read_csv(path.df.image) %>% 
  filter(region < 0) %>%
  filter(value > 0) %>% 
  filter(instance %in% c(1,2)) %>% 
  distinct() %>% 
  pivot_wider(names_from = 'measure', values_from = 'value') %>% 
  mutate(diagnosis = factor(diagnosis, levels = c('NC', 'SCD', 'MCI', 'AD'))) %>% 
  select(age=age.scan,sex, instance, APOE4, `Global efficiency` = NE, `Local efficiency` = NLEgretna, `generalized Local efficiency` = NLE, diagnosis) %>% 
  mutate(APOE4 = case_when(
    APOE4 %in% c(24, 34, 44) ~ 'carrier',
    APOE4 %in% c(22, 23, 33) ~ 'non-carrier',
    TRUE ~ NA_character_
  )) %>% 
  mutate(APOE4 = factor(APOE4, levels = c('non-carrier', 'carrier')))


table1::table1(
  ~age+sex+APOE4+`Global efficiency`+`Local efficiency`+`generalized Local efficiency`|instance + diagnosis, 
  data = df %>% mutate(instance = instance - 1) %>% mutate(instance = ifelse(instance == 0, 'baseline', str_c('follow-up ', instance))) %>% na.omit(),
  extra.col = list(`Statistics`=statistic.and.pval)
) %>% 
  as.data.frame(demo.table) %>% 
  write_csv(str_c(path.result, '/demo_image.csv'))

# Demographic of cross-sectional set ---------------------------------------------

df <- read_csv(path.df.cross.sectional) %>% 
  filter(region < 0) %>% 
  distinct() %>% 
  pivot_wider(names_from = 'measure', values_from = 'value') %>% 
  select(age=age.scan, sex, education, education:GFAP, `Global efficiency` = NE, `Local efficiency` = NLEgretna, `generalized Local efficiency` = NLE, diagnosis = diagnosis.blood) %>% 
  mutate(diagnosis = factor(diagnosis, levels = c('NC', 'SCD', 'MCI', 'AD'))) %>% 
  mutate(APOE4 = case_when(
    APOE4 %in% c(24, 34, 44) ~ 'carrier',
    APOE4 %in% c(22, 23, 33) ~ 'non-carrier',
    TRUE ~ NA_character_
  ))


table1::table1(
  ~sex+age+sex+education+MMSE+`MOCA-B`+APOE4+Aβ42+Aβ40+pTau181+NfL+GFAP+`pTau/Aβ42`+`Aβ42/Aβ40`+`Global efficiency`+`Local efficiency`+`generalized Local efficiency`|diagnosis, 
  data = df,
  extra.col = list(`Statistics`=statistic.and.pval)
) %>% 
  as.data.frame(demo.table) %>% 
  write_csv(str_c(path.result, '/demo_cross_sectional.csv'))
