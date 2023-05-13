setwd('...')
source('codes/000 - functions and paths.R')
path.result <- create_when_absent('results/401-network and (cognition and marker)')
REGION.COMPONENT = -3

# Reads the data while accounting for the effects of covariates -----------

# Network efficiency
data.efficiency <- read_csv(path.df.cross.sectional) %>% 
  select(id, region, name_efficiency = measure, efficiency = value) %>% 
  mutate(name_efficiency = ifelse(name_efficiency == 'NE', 'global efficiency', ifelse(name_efficiency == 'NLE', 'generalised local efficiency', 'local efficiency'))) %>% 
  rbind(
    read_csv('results/400-network PCA/component score.csv') %>% 
      mutate(
        id,
        region = REGION.COMPONENT,
        name_efficiency = 'generalised local efficiency', 
        efficiency = component,
        .keep = 'none'
      )
  ) %>% 
  # Add covariates
  left_join(
    data.efficiency <- read_csv(path.df.cross.sectional) %>% 
      select(id, diagnosis.scan, age.scan, sex, education, APOE4) %>% 
      distinct(),
    by = 'id'
  ) %>% 
  mutate(APOE4 = case_when(
    APOE4 %in% c(24, 34, 44) ~ 'carrier',
    APOE4 %in% c(22, 23, 33) ~ 'non-carrier',
    TRUE ~ NA_character_
  )) %>% mutate(APOE4 = factor(APOE4, levels = c('non-carrier', 'carrier'))) %>% 
  mutate(sex = factor(sex, levels = c('f', 'm'))) %>% 
  group_by(region, name_efficiency) %>%
  do(data.frame(., efficiency.adj = rmv(., efficiency ~ age.scan + sex + education + APOE4))) %>%
  ungroup()


# Cognition and blood markers
data.response <- read_csv(path.df.cross.sectional) %>% 
  select(id, age.blood, diagnosis.blood, diagnosis.scan, sex, education, APOE4, MMSE, MOCA = 'MOCA-B', Aβ42:GFAP) %>% 
  mutate(APOE4 = case_when(
    APOE4 %in% c(24, 34, 44) ~ 'carrier',
    APOE4 %in% c(22, 23, 33) ~ 'non-carrier',
    TRUE ~ NA_character_
  )) %>% 
  pivot_longer(MMSE:GFAP, names_to = 'name_response', values_to = 'response') %>%
  distinct() %>% 
  mutate(sex = factor(sex, levels = c('f', 'm'))) %>% 
  mutate(APOE4 = factor(APOE4, levels = c('non-carrier', 'carrier'))) %>% 
  group_by(name_response) %>% 
  do(data.frame(., response.adj = rmv(., response ~ age.blood + sex + education + APOE4))) %>%
  ungroup()


# Correlations between fitted efficiency and markers ----------------------

# Global and component
data.whole.brain <- 
  data.efficiency %>% filter(region < 0) %>% 
  mutate(name_efficiency = if_else(region == REGION.COMPONENT, str_c('component of ', name_efficiency), name_efficiency)) %>% 
  select(id, name_efficiency, efficiency = efficiency.adj) %>% 
  left_join(
    data.response %>% select(id, name_response, response = response.adj, diagnosis = diagnosis.blood) %>% 
      mutate(diagnosis = factor(diagnosis, levels = c('NC', 'SCD', 'MCI', 'AD'))) %>% 
      mutate(name_response = factor(name_response, levels = c("Aβ40","Aβ42","Aβ42/Aβ40","GFAP","NfL","pTau181","pTau/Aβ42","MMSE","MOCA")))
  )

size.point = 1.5
size.line = 0.5
size.cor = 2
path.plot <- create_when_absent(str_c(path.result, '/plot_cor'))
walk(unique(data.whole.brain$name_efficiency), function(EFFICIENCY) {
  tbl <- data.whole.brain %>% filter(name_efficiency == EFFICIENCY) %>% na.omit()
  
  p <- ggplot(data = tbl, aes(x=efficiency, y=response)) +
    facet_wrap(~name_response, scales = 'free', ncol = 3) +
    geom_point(aes(fill = diagnosis, color = diagnosis), 
               size = size.point) +
    ggsci::scale_color_jco() +
    ggnewscale::new_scale_color() +
    geom_smooth(aes(group = 1),fill="#F2EAE6", method=lm, se=TRUE, size = size.line)+
    ggpubr::stat_cor(aes(group = 1), method = "pearson", size = size.cor) +
    ggsci::scale_color_lancet() +
    theme_classic() 
  
  ggsave(str_c(path.plot, '/', EFFICIENCY, '.pdf'), p , width = 11, height = 10)

})