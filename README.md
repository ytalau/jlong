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
    ## 1  8.114381  1    0      1     66.04639
    ## 2 11.696321  1    1      1     66.04639
    ## 3 12.048433  1    2      1     66.04639
    ## 4 15.178319  1    3      1     66.04639
    ## 5 19.413096  1    4      1     66.04639
    ## 6 17.537775  1    5      1     66.04639

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
    ## 1 33.41925  1    0      1     66.04639           20    0
    ## 2 34.63046  1    1      1     66.04639           20    0
    ## 3 35.67740  1    2      1     66.04639           20    0
    ## 4 29.83353  1    3      1     66.04639           20    0
    ## 5 30.35575  1    4      1     66.04639           20    0
    ## 6 30.35968  1    5      1     66.04639           20    0

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
                  formula.primod = score ~ gender + baseline_age + 
                      baseline_edu + gene + time  + gene:time  + traj +  
                  traj:time  + (time | id),
                  formula.submod = biomk ~ gender + baseline_age + time +
                      (time | id),
                  var.shared = "traj")

summary(res2)
```

    ## 
    ## Call:
    ## jlong(formula.primod = score ~ gender + baseline_age + baseline_edu + 
    ##     gene + time + gene:time + traj + traj:time + (time | id), 
    ##     formula.submod = biomk ~ gender + baseline_age + time + (time | 
    ##         id), data.primod = example_dat$primary, data.submod = example_dat$biomk, 
    ##     var.shared = "traj")
    ## 
    ## Submodel Coefficients:
    ##               Estimate   Std.Err   Wald   Pr(>z)   
    ## (Intercept)  13.563242  4.346967 9.7354 0.001808 **
    ## gender1      -0.945353  1.138385 0.6896 0.406294   
    ## baseline_age  0.034940  0.063278 0.3049 0.580833   
    ## time          0.180725  0.076431 5.5911 0.018052 * 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Primary Model Coefficients:
    ##               Estimate   Std.Err     Wald    Pr(>z)    
    ## (Intercept)  30.833937  2.908782 112.3661 < 2.2e-16 ***
    ## gender1       0.601532  0.448534   1.7986 0.1798861    
    ## baseline_age -0.093815  0.036089   6.7577 0.0093345 ** 
    ## baseline_edu  0.171857  0.068273   6.3363 0.0118290 *  
    ## gene1         0.140187  0.845063   0.0275 0.8682443    
    ## time         -0.377745  0.105108  12.9159 0.0003258 ***
    ## gene1:time   -0.179873  0.202056   0.7925 0.3733504    
    ## traj          3.804733  0.742007  26.2925 2.934e-07 ***
    ## time:traj     0.022440  0.160807   0.0195 0.8890160    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## 
    ## Submodel Random Effects:
    ##             (Intercept)      time
    ## (Intercept)   28.856943 0.2311660
    ## time           0.231166 0.4216604
    ## 
    ## Primary Model Random Effect:
    ##             (Intercept)       time
    ## (Intercept)   1.7656466 -0.1609224
    ## time         -0.1609224  0.3365222

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
    ##               Estimate   Std.Err   Wald   Pr(>z)   
    ## (Intercept)  13.563242  4.346967 9.7354 0.001808 **
    ## gender1      -0.945353  1.138385 0.6896 0.406294   
    ## baseline_age  0.034940  0.063278 0.3049 0.580833   
    ## time          0.180725  0.076431 5.5911 0.018052 * 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Primary Model Coefficients:
    ##               Estimate   Std.Err     Wald    Pr(>z)    
    ## (Intercept)  32.099681  3.151288 103.7589 < 2.2e-16 ***
    ## gender1       0.240999  0.495376   0.2367  0.626615    
    ## baseline_age -0.112178  0.035840   9.7967  0.001748 ** 
    ## baseline_edu  0.198433  0.084138   5.5621  0.018353 *  
    ## gene1        -0.329853  0.747013   0.1950  0.658805    
    ## traj          3.831815  0.697400  30.1888  3.92e-08 ***
    ## time         -0.466354  0.189960   6.0270  0.014088 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## 
    ## Submodel Random Effects:
    ##             (Intercept)      time
    ## (Intercept)   28.856943 0.2311660
    ## time           0.231166 0.4216604
    ## 
    ## Primary Model Random Effect:
    ##             (Intercept)
    ## (Intercept)    3.091165
