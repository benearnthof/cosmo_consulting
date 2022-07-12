library(rcosmo)
library(INLA) 
require(FITSio)
library(data.table)
library(Matrix)
library(ggplot2)
list.of.packages <- c("gapminder", "gt", "tidyverse")

# install required packages, if necessary, and load them --
{
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  lapply(list.of.packages, require, character.only = TRUE)
}