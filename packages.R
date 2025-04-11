#load packages

# R version 4.4.1
require(renv)#1.1.4

renv::restore()
renv::status()

require(dplyr)#1.1.4
require(tidyr)#1.3.0
require(taxize)#0.9.99
require(vroom)#1.5.5
require(ape)#5.5
require(PhyloMeasures)#2.1
require(ggplot2)#3.4.4
require(corrplot)#0.91
require(tm)#0.7.8
require(wordcloud)#2.6
require(cowplot)#1.1.1
require(rstan)#2.35.0.900
require(brms)#2.21.6
require(bayesplot)#1.8.1
require(xtable)#1.8.4
require(ggtree)#3.0.4
require(patchwork)#1.1.2
require(caper)#1.0.1
require(lubridate)#1.8.0
require(parameters)#0.21.6
require(knitr)#1.4.9
require(kableExtra)#1.3.4
require(tidybayes)#3.0.6
require(egg); packageVersion("egg")#0.4.5
require(forcats); packageVersion("forcats")#0.5.1
require(magick); packageVersion("magick")#2.7.3

# renv::install("phangorn@2.8.0")
# renv::install("treeio@1.16.2")
# renv::install("BiocVersion@3.13.1")
