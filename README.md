jlong
================

The `R` package **jlong** provides a main function `jlong` to connect
two linear mixed-effects models via a joint model method.

## Installation

You can install the development version at
[GitHub](https://github.com/ytalau/jlong) using the following code,

``` r
if (!require(devtools)) install.packages("devtools")
devtools::install_github("https://github.com/ytalau/jlong")
```

``` r
library(jlong)
```

## Simulated Data

### Submodel/Biomarker Data

There are some baseline variables contained in the submodel data:

- `biomk`: biomarker value
- `id`: unique subject identifier
- `time`: the visit time
- `gender`: 0 corresponds to female and 1 corresponds to male.
- `baseline_age`: age at time = 0

``` r
head(example_dat$biomk)
```

    ##       biomk id time gender baseline_age
    ## 1  8.375319  1    0      0     74.94612
    ## 2  9.325451  1    1      0     74.94612
    ## 3  7.637588  1    2      0     74.94612
    ## 4  3.784104  1    3      0     74.94612
    ## 5 -2.337580  1    4      0     74.94612
    ## 6 -6.781847  1    5      0     74.94612

There are additional baseline variables contained in the primary model
data:

- `score`: clinical outcome score
- `id`: unique subject identifier which needs to match the `id` in the
  submodel
- `time`: the visit time which could be different from the visit time in
  submodel
- `gender`: 0 corresponds to female and 1 corresponds to male.
- `baseline_age`: age at time = 0
- `baseline_edu` number of years of education completed at time = 0
- `gene`: a binary variable indicating the existence of certain genes

``` r
head(example_dat$primary)
```

    ##      score id time gender baseline_age baseline_edu gene
    ## 1 13.66661  1    0      0     74.94612           18    1
    ## 2 17.74032  1    1      0     74.94612           18    1
    ## 3 16.34908  1    2      0     74.94612           18    1
    ## 4 18.06208  1    3      0     74.94612           18    1
    ## 5 20.19826  1    4      0     74.94612           18    1
    ## 6 19.96913  1    5      0     74.94612           18    1

## True Model: Primary model has random slope and random intercept

The `jlong` function is used to jointly fit two linear mixed models.
There are five argument needed for the `jlong` to function.

- `data.primod`: a data frame / data.table for the primary model
- `data.submod`: a data frame / data.table for the biomarker model
- `formula.primod`: an object of class “formula” (or one that can be
  coerced to that class): a symbolic description of the primary model to
  be fitted; the syntax of the formula is similar to the functions in
  `lme4` package
- `formula.submod`: an object of class “formula” (or one that can be
  coerced to that class): a symbolic description of the sub model to be
  fitted similar to the functions in `lme4` package
- `var.shared`: a string specifying the variable name shared in both
  models

``` r
res2 <- jlong(data.primod = example_dat$primary,
                  data.submod = example_dat$biomk,
                  formula.primod = score ~ gender + baseline_age + time +
                      baseline_edu + traj + traj:time + gene + (time | id),
                  formula.submod = biomk ~ gender + baseline_age + time +
                      (time | id),
                  var.shared = "traj")

summary(res2)
```

    ## 
    ## Call:
    ## jlong(formula.primod = score ~ gender + baseline_age + time + 
    ##     baseline_edu + traj + traj:time + gene + (time | id), formula.submod = biomk ~ 
    ##     gender + baseline_age + time + (time | id), data.primod = example_dat$primary, 
    ##     data.submod = example_dat$biomk, var.shared = "traj")
    ## 
    ## Submodel Coefficients:
    ##                Estimate    Std.Err   Wald   Pr(>z)   
    ## (Intercept)  17.0291237  5.7322582 8.8254 0.002971 **
    ## gender1      -1.1167113  1.1780046 0.8986 0.343145   
    ## baseline_age  0.0018813  0.0798789 0.0006 0.981210   
    ## time          0.0934905  0.0925551 1.0203 0.312444   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Primary Model Coefficients:
    ##               Estimate   Std.Err     Wald    Pr(>z)    
    ## (Intercept)  26.153467  2.926162  79.8844 < 2.2e-16 ***
    ## gender1      -0.282431  0.527825   0.2863 0.5925909    
    ## baseline_age -0.102168  0.035020   8.5111 0.0035298 ** 
    ## baseline_edu  0.377730  0.099078  14.5349 0.0001376 ***
    ## gene1        -0.678778  0.514809   1.7385 0.1873345    
    ## time         -0.132861  0.065816   4.0751 0.0435200 *  
    ## traj          3.136236  0.285822 120.3996 < 2.2e-16 ***
    ## time:traj    -0.146780  0.079891   3.3755 0.0661722 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## 
    ## Submodel Random Effects:
    ##             (Intercept)       time
    ## (Intercept)  34.0681412 -0.1166917
    ## time         -0.1166917  0.6980659
    ## 
    ## Primary Model Random Effect:
    ##             (Intercept)      time
    ## (Intercept)   2.7331717 0.1536701
    ## time          0.1536701 0.2879756

### Simpler Model: Drop the random slope in the primary model

Here we drop the random slope in the primary model, which means the
`time` variable in the primary data set is now a time-dependent
variable.

``` r
res2 <- jlong(data.primod = example_dat$primary,
                  data.submod = example_dat$biomk,
                  formula.primod = score ~ gender + baseline_age + time +
                      baseline_edu + traj + gene + (1 | id),
                  formula.submod = biomk ~ gender + baseline_age + time +
                      (time | id),
                  var.shared = "traj")

summary(res2)
```

    ## 
    ## Call:
    ## jlong(formula.primod = score ~ gender + baseline_age + time + 
    ##     baseline_edu + traj + gene + (1 | id), formula.submod = biomk ~ 
    ##     gender + baseline_age + time + (time | id), data.primod = example_dat$primary, 
    ##     data.submod = example_dat$biomk, var.shared = "traj")
    ## 
    ## Submodel Coefficients:
    ##                Estimate    Std.Err   Wald   Pr(>z)   
    ## (Intercept)  17.0291237  5.7322582 8.8254 0.002971 **
    ## gender1      -1.1167113  1.1780046 0.8986 0.343145   
    ## baseline_age  0.0018813  0.0798789 0.0006 0.981210   
    ## time          0.0934905  0.0925551 1.0203 0.312444   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Primary Model Coefficients:
    ##               Estimate   Std.Err    Wald    Pr(>z)    
    ## (Intercept)  28.366496  2.873588 97.4456 < 2.2e-16 ***
    ## gender1      -0.893070  0.550239  2.6343  0.104577    
    ## baseline_age -0.100342  0.032417  9.5815  0.001965 ** 
    ## baseline_edu  0.257904  0.101493  6.4572  0.011050 *  
    ## gene1        -0.597922  0.564361  1.1225  0.289387    
    ## traj          2.782110  0.334494 69.1786 < 2.2e-16 ***
    ## time         -0.144654  0.119213  1.4723  0.224976    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## 
    ## Submodel Random Effects:
    ##             (Intercept)       time
    ## (Intercept)  34.0681412 -0.1166917
    ## time         -0.1166917  0.6980659
    ## 
    ## Primary Model Random Effect:
    ##             (Intercept)
    ## (Intercept)     5.35149
