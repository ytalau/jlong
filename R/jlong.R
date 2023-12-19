#' Fit Linear Mixed-Effects Model with Covariate from another Linear Mixed-Effects Model
#'
#' @description This function connects two linear mixed-effects models via a joint model.
#'
#' @param formula.primod A two-sided formula for the primary response variable
#' @param formula.submod A two-sided formula for the sub model (biomarker data model)
#' @param data.primod Primary response data
#' @param data.submod Sub model data
#' @param var.shared The name of the variable that both models share but remember
#'              to put this name in the formula for the response variable
#' @param ...
#'
#' @return An object of class \code{jlong}
#' @export jlong
#'
jlong <- function(formula.primod, formula.submod, data.primod, data.submod,
                  var.shared, ...) {
    call <- match.call()
    dat.sub <- data.process2(formula = formula.submod, data = data.submod)

    ## as placeholder
    data.primod[, var.shared] <- 1
    dat.pri <- data.process2(formula = formula.primod, data = data.primod,
                             var.shared = var.shared)

    mod <- estMod2(data.sub = dat.sub, data.pri = dat.pri)
    m <- length(c(mod$beta_1r, mod$beta_v))
    convergence_message <-
        switch(unique(mod$convergence),
               "Convergence of function values has been achieved",
               paste("The relative distance between consecutive solutions",
                     "is smaller than specified xtol value without convergence",
                     "of function values"),
               "No better point found without convergence of function values",
               "Iteration limit maxit exceeded without convergence")
    coeffs.submod <- c(mod$beta_1r, mod$beta_v)
    names(coeffs.submod) <- c(colnames(dat.sub$tind[[1]]),
                              colnames(dat.sub$tdep[[1]]))
    coeffs.primod <- c(mod$beta_u, mod$beta_xt, mod$beta_z)
    names(coeffs.primod) <- c(colnames(dat.pri$tind[[1]]),
                              colnames(dat.pri$tind_shared[[1]]),
                              colnames(dat.pri$tdep[[1]]))
    Sigma_gm <- mod$Sigma_gm
    Sigma_b <- mod$Sigma_b
    dimnames(Sigma_gm) <- list(c("(Intercept)", dat.sub$random.vars),
                                   c("(Intercept)", dat.sub$random.vars))
    dimnames(Sigma_b) <- list(c("(Intercept)", dat.pri$random.vars),
                                  c("(Intercept)", dat.pri$random.vars))
    fit <- list(call = call,
                inf.submod = mod$ses[1:m],
                inf.primod = mod$se[(m + 1) : length(mod$ses)],
                formula.primod = formula.primod,
                formula.submod = formula.submod,
                coefficients.submod = coeffs.submod,
                coefficients.primod = coeffs.primod,
                sigma.submod = sqrt(mod$hsigma_e2),
                sigma.primod = sqrt(mod$hsigma_2),
                rcov.submod = Sigma_gm,
                rcov.primod = Sigma_b,
                convergence_code = mod$convergence,
                convergence_message = convergence_message)
    class(fit) <- "jlong"
    fit
}


print.jlong <- function(obj, ...) {
    cat("\nCall:\n")
    dput(obj$call)
    ## coefficients
    cat("\nSubmodel Coefficients:\n")
    print(obj$coefficients.submod)
    cat("\nPrimary Model Coefficients:\n")
    print(obj$coefficients.primod)
    cat("\n")
    cat(obj$convergence_message)
    invisible(obj)

}

#' @export
summary.jlong <- function(obj, ...) {
    coef.matrix.submod <- data.frame(Estimate = unname(obj$coefficients.submod),
                              Std.Err =  unname(obj$inf.submod))
    coef.matrix.primod <- data.frame(Estimate = unname(obj$coefficients.primod),
                                     Std.Err =  unname(obj$inf.primod))
    coef.matrix.submod$wald <- (coef.matrix.submod$Estimate /
                                    coef.matrix.submod$Std.Err)^2
    coef.matrix.primod$wald <- (coef.matrix.primod$Estimate /
                                    coef.matrix.primod$Std.Err)^2
    coef.matrix.submod$chi.squared <- 1 - pchisq(coef.matrix.submod$wald,
                                                 df = 1)
    coef.matrix.primod$chi.squared <- 1 - pchisq(coef.matrix.primod$wald,
                                                 df = 1)
    colnames(coef.matrix.submod) <- c("Estimate", "Std.Err", "Wald", "Pr(>z)")
    colnames(coef.matrix.primod) <- c("Estimate", "Std.Err", "Wald", "Pr(>z)")
    rownames(coef.matrix.submod) <- names(obj$coefficients.submod)
    rownames(coef.matrix.primod) <- names(obj$coefficients.primod)
    out <- list(formula = obj$formula,
                call = obj$call,
                rcov.submod = obj$rcov.submod,
                rcov.primod = obj$rcov.primod,
                formula.primod = obj$formula.primod,
                formula.submod = obj$formula.submod,
                convergence_code = obj$convergence_code,
                convergence_message = obj$convergence_message,
                coef.matrix.submod = coef.matrix.submod,
                coef.matrix.primod = coef.matrix.primod)
    class(out) <- "summary.jlong"
    out
}

#' @export
print.summary.jlong <- function(obj, ...) {
    cat("\nCall:\n")
    dput(obj$call)
    ## coefficients
    cat("\nSubmodel Coefficients:\n")
    printCoefmat(obj$coef.matrix.submod)
    cat("\nPrimary Model Coefficients:\n")
    printCoefmat(obj$coef.matrix.primod)
    cat("\n")
    cat("\nSubmodel Random Effects:\n")
    print(obj$rcov.submod)
    cat("\nPrimary Model Random Effect:\n")
    print(obj$rcov.primod)
    cat("\n")
    invisible(obj)
}
