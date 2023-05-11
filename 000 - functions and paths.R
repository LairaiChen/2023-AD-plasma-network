setwd('...') # Your working directory
library(tidyverse)
library(lubridate)
library(tidyverse)
library(lubridate)
path.df.blood <- '...' # Plasma markers
path.df.image <- '...' # Network markers
path.df.cross.sectional <- '...' # Participants with both plasma markers and network markers
bna.description <- read_csv('data/bna_246_description_BNV.csv')
# Files ----------------------
create_when_absent <- function(path) {
  if(!dir.exists(path)) dir.create(path)
  return(path)
}

get_file_name <- function(x) {
  strsplit(x, '\\.')[[1]][1]
}

write_connetome_lairai <- function(network, path.save, .bnv = T) {
  if(.bnv) path.save <- str_replace(path.save, '\\.txt', '\\.edge')
  write.table(network, path.save, row.names = F, col.names = F, quote = F)
}


# Statistics -------------------------------------------------------

rmv <- function(data, formula) {
  # Convert the formula to a formula object
  formula <- as.formula(formula)
  
  # Fit a linear model with the given formula using the input data
  model <- lm(formula, data = data)
  
  # Extract the residuals
  residuals <- residuals(model)
  
  # Adjust the residuals to have the same mean as the original response variable
  y_mean <- mean(model$model[[1]])
  residuals_adj <- residuals + y_mean
  
  # Return the adjusted residuals as a vector
  return(residuals_adj)
}

