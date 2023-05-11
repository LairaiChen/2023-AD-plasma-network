setwd('...') # Your working directory
source('codes/000 - functions and paths.R')
path.result <- create_when_absent('results/200-network PCA')

# Retrieves progression-related regions -----------------------------------
region.progression.related <- read_csv('results/100-trend analyses/result_trend_report_all.csv') %>% 
  filter(region > 0 & p.adj.L <= 0.05) %>% 
  pull(region)

# Computes components for progression-related regions ---------------------
efficiency.data <- read_csv(path.df.cross.sectional) %>% 
  filter(measure == 'NLE') %>% 
  filter(region %in% region.progression.related) %>% 
  select(region, value, id) %>% 
  distinct() %>% 
  pivot_wider(names_from = 'region', values_from = 'value') %>% 
  na.omit()

# Gets regional measures in blood-marker related ROIs
#	matrix; Data matrix (each row is an observation, each column is a variable)
pc.parallel <- psych::fa.parallel(efficiency.data %>% select(-id), fa='pc', n.iter=10000, show.legend=FALSE)
df.parallel <- 
  tibble(
    x = 1:length(pc.parallel$fa.values),
    `Principle values` = pc.parallel$pc.values,
    `Simulated data` = pc.parallel$pc.sim,  
    `Resampled data` = pc.parallel$pc.simr) %>% 
  pivot_longer(cols = `Principle values`:`Resampled data`,
               values_to = 'pc',
               names_to = 'type')
plot.parallel <- ggplot(df.parallel, aes(x, pc)) + 
  geom_point(aes(colour = type, shape = type)) + 
  geom_line(aes(colour = type)) + 
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'red') +
  geom_vline(xintercept = as.numeric(pc.parallel$ncomp), linetype = 'dashed', color = 'blue') + 
  labs(
    title = str_c('Parallel analysis suggests ', pc.parallel$ncomp, ' components'),
    x = 'Principle component numbers',
    y = 'Eigen values'
  )  + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + # Mediate the title
  ggsci::scale_color_lancet()
plot.parallel

ggsave(str_c(path.result, '/scree.pdf'), plot.parallel, width=6, height=5)

# Do pca and saves the scores and weights
NCOMP = pc.parallel$ncomp
pca.result <- psych::principal(efficiency.data %>% select(-id), nfactors=NCOMP, rotate='geominT', scores=TRUE)
# Visualize the weights
pca.weights <- 
  as.matrix(pca.result$weights) %>% 
  reshape::melt.matrix() %>% 
  set_names(c("x", "y", "value")) %>% 
  mutate(y = factor(y, levels = c(str_c('RC', 1:NCOMP)))) %>% 
  rename(region = x) %>% 
  left_join(
    bna.description %>% 
      mutate(acronym = stringr::str_extract(acronym, pattern = '^[a-zA-Z]+_[LR]')) %>%
      mutate(acronym = str_replace(acronym, pattern = "^([a-zA-Z]+)_(L|R)$", replacement = "\\1(\\2)")) %>% 
      mutate(acronym = str_c(region, ', ', acronym)) %>% 
      select(region, acronym, `full name`)
  ) 
acronym.order <- pca.weights %>% select(region, acronym) %>% distinct() %>% arrange(region) %>% pull(acronym)
pca.weights <- pca.weights %>% mutate(acronym = factor(acronym, levels = acronym.order))

plot.pca.loading <- ggplot(pca.weights, aes(x = y, y = acronym, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#3636B3", high = "#AC3636", mid = 'white') +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20)) +
  labs(
    x = 'Principle component',
    y = 'Loading weights'
  ) +
  theme_classic() +
  ggsci::scale_color_lancet()
ggsave(str_c(path.result, '/loading.pdf'), plot.pca.loading, width=6, height=13)

data.pca <- 
  cbind(id=efficiency.data$id, pca.result$scores) %>% 
  as_tibble() %>% 
  pivot_longer(-id, names_to = 'name_component', values_to = 'component') %>% 
  mutate(name_component = factor(name_component, levels = c(str_c('PC', 1:NCOMP))))

write_csv(data.pca, str_c(path.result, '/component score.csv'))
