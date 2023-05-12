setwd('...')
source('codes/000 - functions and paths.R')
path.result <- create_when_absent('results/300-marker and cognition')

data.marker <- read_csv(path.df.blood) %>% 
  select(id, age=age.blood, diagnosis, sex, education, APOE4, Aβ42:GFAP, MMSE, MOCA = 'MOCA-B') %>% 
  mutate(APOE4 = case_when(
    APOE4 %in% c(24, 34, 44) ~ 'carrier',
    APOE4 %in% c(22, 23, 33) ~ 'non-carrier',
    TRUE ~ NA_character_
  )) %>% 
  mutate(sex = factor(sex, levels = c("f", "m"), labels = c("Female", "Male"))) %>% 
  mutate(APOE4 = factor(APOE4, levels = c("non-carrier", "carrier")) ) %>% 
  pivot_longer(Aβ42:GFAP, names_to = 'name_marker', values_to = 'marker') %>% 
  pivot_longer(MMSE:MOCA, names_to = 'name_cognition', values_to = 'cognition') %>% 
  mutate(diagnosis = factor(diagnosis, levels = c('NC', 'SCD', 'MCI', 'AD'))) %>% 
  na.omit() %>% 
  # Adjust cognition performances
  group_by(name_marker, name_cognition) %>% 
  do(data.frame(., cognition.adj = rmv(., cognition ~ age + sex + education + APOE4))) %>%
  ungroup()


# Correlation and visualization -----------------------------------------------------------
size.point = #Custome
size.line = #Custome
size.cor = #Custome
data.plot <- data.marker %>% select(cognition = cognition.adj, marker, diagnosis, name_marker, name_cognition)
p <- ggplot(data.plot, aes(x=marker, y=cognition)) +
  geom_point(aes(fill = diagnosis, color = diagnosis), 
             size = size.point) +
  ggsci::scale_color_jco() +
  ggnewscale::new_scale_color() +
  geom_smooth(aes(group = 1),fill="#F2EAE6", method=lm, se=TRUE, size = size.line)+
  ggpubr::stat_cor(aes(group = 1), method = "pearson", size = size.cor, label.x = 0) +
  #ggsci::scale_color_lancet() +
  theme_classic() +
  facet_wrap(name_marker ~ name_cognition, scales = 'free', ncol = 4)

ggsave(str_c(path.result, '/cor.pdf'), p , width = 11, height = 10)
