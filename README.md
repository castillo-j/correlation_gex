Correlation analysis
================

## Correlation analysis using large gene expression datasets

Analysis of correlations between mutational indicators (e.g., tumor
mutational burden) and gene expression across large-scale datasets. This
repository contains code, data processing workflows, and visualization
tools for exploring these relationships.

# Structure the R project as follow:

```
├── README.md
├── data              <- Raw data, avoid editing directly.
├── data_interim      <- Intermidiate data that has be transformed.
├── data_processed    <- Final data sets for making plots and sharing.
├── scrips
|   ├── 1_process.r     <- Script to format raw data for correlation analysis (varies by input)
|   ├── 2_correlation.r <- Correlation script
|   ├── 3_plots.r       <- Script for plots and table outputs
├── results             <- Generated plots, tables, ppts
```

## Including Code

You can include R code in the document as follows:

``` r
summary(cars)
```

    ##      speed           dist       
    ##  Min.   : 4.0   Min.   :  2.00  
    ##  1st Qu.:12.0   1st Qu.: 26.00  
    ##  Median :15.0   Median : 36.00  
    ##  Mean   :15.4   Mean   : 42.98  
    ##  3rd Qu.:19.0   3rd Qu.: 56.00  
    ##  Max.   :25.0   Max.   :120.00

## Including Plots

You can also embed plots, for example:

![](README_files/figure-gfm/pressure-1.png)<!-- -->

Note that the `echo = FALSE` parameter was added to the code chunk to
prevent printing of the R code that generated the plot.
