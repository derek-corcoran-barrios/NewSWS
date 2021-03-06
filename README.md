Documentation for the SQS analysis for the manuscript Increase in
ecological complexity of marine ecosystems during the Paleozoic-Mesozoic
transition
================

## How to use this guide

this guide will guide you on taking the `SQS.R` file in this same
repository and you will be able to get the same results shown in this
document

### Inputs

The only input needed is the csv in this same repository called
`EcoData2019.csv`.

### Packages needed

The packages needed for running this script are *tidyverse*,
*PBSmodelling*, *foreach* and *doParallel*

### The script and how to change it

The script `SQS.R` is commented, this script creates two functions:

  - **RichnessWhile:** This function takes a database, and it samples it
    until it reaches a certain value of *Good’s U*. This function takes
    two parameters:
      - *Data:* A data frame with the data
      - *quorum:* The desired *Good’s U* as a threshold

This function returns a list with a Graph of the trajectory of the
*Good’s U* sample by sample, the Ensemble of genus, a data frame with
the *GoodsU*, genus Richness, genus richness for temperate genus, genus
richness for tropical genus, ecocode richness, ecocode richness for
temperate genus and ecocode richness for tropical genus.

  - **ResampleWhile:** This function takes a database, and it samples it
    until it reaches a certain value of *Good’s U*. through several
    simulations, this functions loops through the *RichnessWhile*
    functions several times and it uses the geometric mean or median to
    calculate the reported Genera richness and ecocode richness. This
    function takes four parameters:
      - *DF:* A data frame with the data
      - *Quorum:* The desired *Good’s U* as a threshold
      - *nsim:* The number of simulations to run
      - *method:* a character, it can be *“Geometric”* or *“Median”*,
        depending on that the reported genera richness and ecocode
        richness will be based on the median or geometric mean of all
        simulations
      - *ncores:* It defaults to NULL, but if not it will use the
        provided number of cores to run the model.

This function returns a data frame with the information of every
simulation, plus the recorded genus Richness, genus richness for
temperate genus, genus richness for tropical genus, ecocode richness,
ecocode richness for temperate genus and ecocode richness for tropical
genus calculated as the geometric mean or median of all simulations.

### The code used for the problem

The code was run in an *m5.12xlarge* AWS linux instance with 48 cores.

Besides using SQS we used bootstrapping to get a confidence interval on
the estimates, if you go to line 148 of the code you can change the
number of bootsraps which was set to 500.

In line 154 we used `set.seed` to make the results reproducible. Every
use of *ResampleWhile* was set to 500 simulations a quorum of 0.4, using
40 cores.

All datasets were saved in `RDS` and `csv` format.

### Visualization

The `Graph.r` script generates the main graph for this result which is
shown in the following figure:

``` r
Results <- read_rds("Results.rds")

Results2 <- read_rds("Results2.rds")

Results3 <- Results %>% mutate(Type = case_when(Meassurement %in% c("Ecos", "EcosTemp", "EcosTrop") ~ "Ecological entities",Meassurement %in% c("Richness", "RichnessTemp", "RichnessTrop") ~ "Richness")) %>% mutate(Meassurement = case_when(Meassurement %in% c("Ecos", "Richness") ~ "Total", Meassurement %in% c("EcosTrop", "RichnessTrop") ~ "Tropical", Meassurement %in% c("EcosTemp", "RichnessTemp") ~ "Temperate")) %>% mutate(Type = relevel(factor(Type), "Richness")) %>% mutate(Meassurement = relevel(factor(Meassurement), "Total"))


Results3 <- bind_rows(Results3)
saveRDS(Results3, "Results3.rds")


ggplot(Results3, aes(y = Mean, x = Time, group = Meassurement)) + geom_ribbon(aes(ymax = Upper, ymin = Lower, fill = Meassurement), alpha = 0.3) + geom_path(aes(color = Meassurement)) + geom_point(aes(color = Meassurement)) + theme_bw() + facet_grid(Type ~ ., scales = "free_y") + theme(legend.position = "bottom") + scale_x_reverse() + xlab("Time [Ma]") + ylab("") + geom_vline(xintercept = 250, lty = 2, color = "red") + geom_vline(xintercept = 200, lty = 2, color = "red")
```

![](README_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->
