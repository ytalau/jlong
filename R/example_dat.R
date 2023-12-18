#' Longitudinal study of childhood growth data
#'
#' The \code{example_dat} data frame has 3537 rows and 6 columns.
#'
#' @details Data contains body mass index (BMI) z-scores for 486 subjects
#' from their outpatient visits.
#' Subjects were categorized as obese or non-obese depending on their average
#' BMI z-scores. This dataset can be used to predict obesity in adolescence
#' (age 12--17) and obesity in adulthood (age 20-29)
#' using the childhood (age 2--7) BMI z-scores and preadolescence (age 7--12)
#' BMI z-scores, respectively. Note that the true long term BMI z-scores
#' cannot be measured due to finite measurement and thus they are subject to
#' measurement error.
#'
#'
#' @format This data frame contains the following columns:
#'
#' \describe{
#' \item{id}{identifer of child}
#' \item{time}{an indicator of age intervals (1 = age 2--7 for childhood
#'  BMI z-scores and age 12--17 for obesity in adolescence,
#'  2 = age 7--12 for preadolescence
#'  BMI z-scores and age 20--29 for obesity in adulthood)}
#' \item{gender}{an indicator of gender (1 = male, 0 = female)}
#' \item{age}{child's age}
#' \item{bmi_z}{BMI z-score}
#' \item{obes}{an indicator of obesity (1 = obese, 0 = none)}
#' }
#'
#' @keywords datasets
#' @examples
#' data(bmidata)
#' \dontrun{
#' library(data.table)
#' library(geepack)
#' dat_cal <- bmidat[order(id, age)]
#' dat_cal <- dat_cal[, `:=`(
#'     W_star = bmi_z - shift(bmi_z, n = 1, fill = 0, type = "lag"),
#'     age_diff = age - shift(age, n = 1, fill = 0, type = "lag")),
#'     by = .(id, time)
#'     ][, .(age_diff = tail(age_diff, .N - 1),
#'           W_star = tail(W_star, .N - 1)), by = .(id, time)]
#'
#' dat_cal <- dat_cal[, `:=`(W_star2 := W_star * shift(W_star, n = 1, fill = NA,
#'                                                     type = "lag"),
#'                           t1 = shift(age_diff, n = 1, fill = NA,
#'                                      type = "lag"), t2 = age_diff,
#'                           t3 = shift(age_diff, n = 1, fill = NA,
#'                                      type = "lag") + age_diff),
#'                    by = .(id, time)]
#'
#' dat_cal <- dat_cal[!is.na(W_star2), ]
#' mome <- function(theta, W = dat_cal$W_star,
#'                  t = dat_cal$age_diff,
#'                  W2 = dat_cal[!is.na(W_star2), W_star2],
#'                  t2 = unique(dat_cal[!is.na(W_star2), .(t1, t2, t3)])) {
#'     sigma <- theta[1]
#'     rho <- theta[2]
#'     mo1 <- sum(W^2 - 2 * sigma^2 * (1 - rho^t))
#'     mo2 <- sum(W2 - sigma^2 * (rho^t2[, 2] - 1 - rho^t2[, 3] + rho^t2[, 1]))
#'     c(mo1, mo2)
#' }
#'
#' res <- nleqslv::nleqslv(c(0.5, 0.5), fn = mome)
#'
#' bmidat <- bmidat[, me_sd := sqrt(res$x[1] *
#'                                      sum(res$x[2]^abs(outer(age, age, FUN = "-"))))/.N,
#'                  by = .(id, time)]
#' me_sd2 <- round(mean(bmidat$me_sd), 3)^2
#'
#' ## fitting models
#'
#' bmidat <- unique(bmidat[, W := mean(bmi_z), by = .(id, time)
#'                         ][, .(id, time, obes, W, me_sd)])
#'
#' init <- geepack::geeglm(obes ~ time + time:W -1,
#'                         family = binomial, data = bmidat,
#'                         corstr = "ar1", id = id)
#'
#' cres <- eiv::eivgmm(obes ~ time + time:W - 1,
#'                     data = bmidat, me.var = c("time1:W", "time2:W"),
#'                     mcov = list(matrix(c(me_sd2, 0, 0, 0), 2, 2),
#'                                 matrix(c(0, 0, 0, me_sd2), 2, 2)),
#'                     time.var = "time", id.var = "id",
#'                     start = init$coefficients,
#'                     modify_inv = 1,
#'                     finsam_cor = 0)
#' summary(cres)
#' }
"example_dat"