
jointLong2 <- function(formula.submod, data.submod,
                      formula.primod, data.primod, var.shared, ...) {

    ## get submod data process done first
    dat.sub <- data.process2(formula = formula.submod, data = data.submod)

    ## as placeholder
    data.primod[, var.shared] <- 1
    dat.pri <- data.process2(formula = formula.primod, data = data.primod,
                            var.shared = var.shared)

    mod <- estMod2(data.sub = dat.sub, data.pri = dat.pri)
    m <- length(c(mod$beta_1r, mod$beta_v))
    ## TODO temporary
    df <- data.frame(vars = c(colnames(dat.sub$tind[[1]]),
                              colnames(dat.sub$tdep[[1]]),
                              "sigma_e",
                              paste0("Sigma_gm",
                                     1:length(mod$Sigma_gm)),
                              colnames(dat.pri$tind[[1]]),
                              colnames(dat.pri$tind_shared[[1]]),
                              colnames(dat.pri$tdep[[1]]),
                              "sigma",
                              paste0("Sigma_b",
                                     1:length(mod$Sigma_b))),
                     vals = c(mod$beta_1r, mod$beta_v,
                              sqrt(mod$hsigma_e2),
                              matrix(mod$Sigma_gm),
                              mod$beta_u, mod$beta_xt, mod$beta_z,
                              sqrt(mod$hsigma_2),
                              matrix(mod$Sigma_b)),
                     ses = c(mod$ses[1:m],
                             rep(NA, length(mod$Sigma_gm) + 1),
                             mod$se[(m + 1) : length(mod$ses)],
                             rep(NA, length(mod$Sigma_b) + 1)),
                     conv = mod$convergence)

    df

}


S_beta_1r <- function(matR, hX, Sigma_gamma_inv, beta_1r) {
    t(matR) %*% Sigma_gamma_inv %*% (hX - matR %*% beta_1r)
}

dS_beta_1r <- function(matR, Sigma_gamma_inv) {
    t(matR) %*% Sigma_gamma_inv %*% (-matR)
}

S_beta_u <- function(matU, Sigma_b_inv, meanvec) {
    t(matU) %*% Sigma_b_inv %*% meanvec
}

dS_beta_u <- function(matU, Sigma_b_inv) {
    t(matU) %*% Sigma_b_inv %*% (-matU)
}

dS_beta_uxt <- function(matUX, Sigma_b_inv, beta_xt, merr.sub.mat,
                        mmX) {
    out <- t(matUX) %*% Sigma_b_inv %*% (-matUX)
    idx <- ncol(matUX) - (length(beta_xt)-1):0
    coef_mats <- kronecker(merr.sub.mat, Sigma_b_inv)

    mat_idx <- which(mmX != 0, arr.ind = T)
    nE <- nrow(mat_idx)
    varfun <- function(i, j, k, l) Sigma_b_inv[i, k] * merr.sub.mat[l, j]
    d_corr_term_vec <- mapply(varfun,
                              i = rep(mat_idx[, 1], each = nE),
                              j = rep(mat_idx[, 2], times = nE),
                              k = mat_idx[, 1], l = mat_idx[, 2])
    d_corr_term <- matrix(d_corr_term_vec, nrow = nE, ncol = nE, byrow = T)
    out[idx, idx] <- out[idx, idx] + d_corr_term
    out

}

S_beta_xt <- function(mathX, Sigma_b_inv, beta_xt,
                      meanvec, merr.sub.mat, mmX) {
    # TODO:more than 2 random effects and Xhat is only
    # involved in some of them
    mmX[mmX != 0] <- beta_xt
    cor_term <- Sigma_b_inv %*% mmX %*% merr.sub.mat
    t(mathX) %*% Sigma_b_inv %*% meanvec + cor_term[mmX != 0]

}




S_beta_z <- function(matZ, heta, A_blck, Sigma_b_inv, Y,
                     invM, projM, meanvec, beta_z, sigma_2) {
    t(matZ) %*% invM %*% (Y - A_blck %*% heta - matZ %*% beta_z) / sigma_2 +
    t(matZ) %*% t(projM) %*% Sigma_b_inv %*% meanvec
}

## placeholder
dS_beta_z <- function(matZ, Sigma_b_inv, invM, projM,  sigma_2) {
    - t(matZ) %*% invM %*% invM %*% matZ  / sigma_2 +
        t(matZ) %*% t(projM) %*% Sigma_b_inv %*% projM %*% matZ
}


S_beta_v <- function(matV, hX, D_blck, Sigma_gamma_inv, Sigma_b_inv,
                     sigma_e2, W, invM, projM, meanvec, meanvec_sub,
                     beta_v, beta_xt
                     ) {
    t(matV) %*% invM %*% (W - D_blck %*% hX - matV %*% beta_v) / sigma_e2 +
        t(matV) %*% t(projM) %*% Sigma_gamma_inv %*% meanvec_sub +
        (t(matV) %*% t(projM))[1,2] %*% t(beta_xt) %*% Sigma_b_inv %*% meanvec

}

## placeholder
dS_beta_v <- function(matV, Sigma_gamma_inv, invM, projM, sigma_e2,
                      beta_xt, Sigma_b_inv) {

    - t(matV) %*% invM %*% invM %*% matV / sigma_e2 +
        t(matV) %*% t(projM) %*% Sigma_gamma_inv %*% projM %*% matV -
        (t(matV) %*% t(projM))[1,2] %*% t(beta_xt) %*% Sigma_b_inv %*%
         beta_xt %*% (t(matV) %*% t(projM))[1,2]
}

S_Sigma <- function(mmX, beta_xt, meanvec, merr, hsigma_2, blck_mat) {
    if (is.null(beta_xt)) {
        b.merr <- lapply(1:n, function(i) 0)
    } else {
        mmX[mmX != 0] <- beta_xt
        b.merr <-  mmX %*% merr %*% t(mmX)
    }
    meanvec %*% t(meanvec) - b.merr - c(hsigma_2) *
        solve(t(blck_mat) %*% blck_mat)
}


jointSolve <- function(par, matR, matV = NULL, matZ = NULL,
                       matX, matU, pms.sub, pms.pri, Y, W,
                       n, merr.sub.mat, A_blck, D_blck,
                       Sigma_gamma_inv, Sigma_b_inv,
                       hsigma_e2, hsigma_2,
                       idx_beta_v, idx_beta_1r, mmX,
                       idx_beta_z, idx_beta_u, idx_beta_xt) {

    beta_1r <- as.matrix(par[1:idx_beta_1r])
    beta_u <- as.matrix(par[(idx_beta_v + 1):idx_beta_u])
    beta_xt <- as.matrix(par[(idx_beta_u + 1):idx_beta_xt])

    ## the structure of Xhat should be given at the first place
    if (is.null(matV[[1]])) {
        hX <- lapply(1:n, function(i) pms.sub$projM[[i]] %*% W[[i]])
    } else {
        beta_v <- as.matrix(par[(idx_beta_1r + 1):idx_beta_v])
        hX <- lapply(1:n, function(i) pms.sub$projM[[i]] %*%
                         (W[[i]] - matV[[i]] %*% beta_v))
    }
    # TODO: how can we remove the hard index
    mathX <- lapply(1:n, function(i) {
         matX[[i]] * hX[[i]][2, ]
    })

    ## need to setup the matrix for measurement error
    if (is.null(matZ[[1]])) {
        heta <- lapply(1:n, function(i) pms.pri$projM[[i]] %*% Y[[i]])
    } else {
        beta_z <- as.matrix(par[(idx_beta_xt + 1):idx_beta_z])
        heta <- lapply(1:n, function(i) pms.pri$projM[[i]] %*%
                         (Y[[i]] - matZ[[i]] %*% beta_z))
    }

    meanvec <- lapply(1:n, function(i) heta[[i]] - matU[[i]] %*% beta_u -
                         mathX[[i]] %*% beta_xt)
    meanvec_sub <- lapply(1:n, function(i) hX[[i]] - matR[[i]] %*% beta_1r)

    S_1r <- lapply(1:n, function(i)
        S_beta_1r(matR = matR[[i]], beta_1r = beta_1r, hX = hX[[i]],
                  Sigma_gamma_inv = Sigma_gamma_inv))

    S_u <- lapply(1:n, function(i)
        S_beta_u(matU = matU[[i]], meanvec = meanvec[[i]],
                  Sigma_b_inv = Sigma_b_inv))

    S_xt <- lapply(1:n, function(i)
        S_beta_xt(mathX = mathX[[i]], Sigma_b_inv = Sigma_b_inv,
                  beta_xt = beta_xt, meanvec = meanvec[[i]],
                  mmX = mmX, merr.sub.mat = merr.sub.mat[[i]]))

    if (is.null(matZ[[1]])) {
        S_z <- NULL
    } else {
        S_z <- lapply(1:n, function(i)
            S_beta_z(matZ = matZ[[i]], heta = heta[[i]],
                     A_blck = A_blck[[i]], Sigma_b_inv = Sigma_b_inv,
                     Y = Y[[i]], sigma_2 = hsigma_2,
                     invM = pms.pri$invM[[i]],
                     projM = pms.pri$projM[[i]],
                     meanvec = meanvec[[i]], beta_z = beta_z))
        S_z <- Reduce(`+`, S_z)
    }

    if (is.null(matV[[1]])) {
        S_v <- NULL
    } else {
        S_v <- lapply(1:n, function(i)
            S_beta_v(matV = matV[[i]], beta_xt = beta_xt,
                     hX = hX[[i]], D_blck = D_blck[[i]],
                     Sigma_gamma_inv = Sigma_gamma_inv,
                     Sigma_b_inv = Sigma_b_inv,
                     sigma_e2 = hsigma_e2,
                     W = W[[i]], invM = pms.sub$invM[[i]],
                     projM = pms.sub$projM[[i]],
                     meanvec = meanvec[[i]],
                     meanvec_sub = meanvec_sub[[i]], beta_v = beta_v))
    }
    S_v <- Reduce(`+`, S_v)
    S_1r <- Reduce(`+`, S_1r)
    S_u <- Reduce(`+`, S_u)
    S_xt <- Reduce(`+`, S_xt)
    c(S_1r, S_v , S_u , S_xt , S_z)

}



cMat <- function(beta_1r, beta_u, beta_xt,
                 beta_v = NULL, beta_z = NULL,
                 matR, matV = NULL, matZ = NULL,
                 matX, matU, pms.sub, pms.pri, Y, W,
                 n, merr.sub.mat, A_blck, D_blck,
                 hX, mathX, heta, mmX, hsigma_e2,
                 hsigma_2,
                 Sigma_gamma_inv, Sigma_b_inv) {

    meanvec <- lapply(1:n, function(i) heta[[i]] - matU[[i]] %*% beta_u -
                          mathX[[i]] %*% beta_xt)
    meanvec_sub <- lapply(1:n, function(i) hX[[i]] - matR[[i]] %*% beta_1r)

    if (is.null(beta_v)) {
        S_v <- NULL
        dS_v <- NULL
    } else {
        S_v <- lapply(1:n, function(i)
            S_beta_v(matV = matV[[i]], beta_xt = beta_xt,
                     hX = hX[[i]], D_blck = D_blck[[i]],
                     Sigma_gamma_inv = Sigma_gamma_inv,
                     Sigma_b_inv = Sigma_b_inv,
                     W = W[[i]], invM = pms.sub$invM[[i]],
                     projM = pms.sub$projM[[i]],
                     meanvec = meanvec[[i]], sigma_e2 = hsigma_e2,
                     meanvec_sub = meanvec_sub[[i]], beta_v = beta_v))
        dS_v <- lapply(1:n, function(i)
            dS_beta_v(matV = matV[[i]], sigma_e2 = hsigma_e2,
                      Sigma_gamma_inv = Sigma_gamma_inv,
                      Sigma_b_inv = Sigma_b_inv,
                      beta_xt = beta_xt,
                      invM = pms.sub$invM[[i]],
                      projM = pms.sub$projM[[i]]))
    }

    if (is.null(beta_z)) {
        S_z <- NULL
        dS_z <- NULL
    } else {
        S_z <- lapply(1:n, function(i)
            S_beta_z(matZ = matZ[[i]], heta = heta[[i]],
                     A_blck = A_blck[[i]], Sigma_b_inv = Sigma_b_inv,
                     Y = Y[[i]], sigma_2 = hsigma_2,
                     invM = pms.pri$invM[[i]],
                     projM = pms.pri$projM[[i]],
                     meanvec = meanvec[[i]], beta_z = beta_z))
        dS_z <- lapply(1:n, function(i)
            dS_beta_z(matZ = matZ[[i]], Sigma_b_inv = Sigma_b_inv,
                      sigma_2 = hsigma_2,
                      invM = pms.pri$invM[[i]], projM = pms.pri$projM[[i]]))
    }

    cvec <- lapply(1:n, function(i) {
        matrix(c(S_beta_1r(matR = matR[[i]], beta_1r = beta_1r, hX = hX[[i]],
                           Sigma_gamma_inv = Sigma_gamma_inv),
                 S_v[[i]],
                 S_beta_u(matU = matU[[i]], meanvec = meanvec[[i]],
                          Sigma_b_inv = Sigma_b_inv),
                 S_beta_xt(mathX = mathX[[i]], Sigma_b_inv = Sigma_b_inv,
                           mmX = mmX,
                           beta_xt = beta_xt, meanvec = meanvec[[i]],
                           merr.sub.mat = merr.sub.mat[[i]]),
                 S_z[[i]]))
    })
    cmat <- Reduce(`+`,
                   lapply(1:n, function(i) cvec[[i]] %*% t(cvec[[i]]))) / n

    dS_1r <- lapply(1:n, function(i)
        dS_beta_1r(matR = matR[[i]],
                   Sigma_gamma_inv = Sigma_gamma_inv))


    dS_uxt <- lapply(1:n, function(i)
        dS_beta_uxt(matUX = cbind(matU[[i]], mathX[[i]]),
                    mmX = mmX,
                    Sigma_b_inv = Sigma_b_inv, beta_xt = beta_xt,
                    merr.sub.mat = merr.sub.mat[[i]]))


    tmp <- list(Reduce(`+`, dS_1r),
                Reduce(`+`, dS_v),
                Reduce(`+`, dS_uxt),
                Reduce(`+`, dS_z))
    tmp <- tmp[!sapply(tmp,is.null)]
    dmat <- as.matrix(Matrix::bdiag(tmp) / n)
    list(cmat, solve(dmat) %*% cmat %*% t(solve(dmat)))
}

## second step: need to calculate pm

PMcals <- function(blck_mat, clsz, n) {
    projM <- lapply(1:n,
                    function(i) solve(t(blck_mat[[i]])%*% blck_mat[[i]]) %*%
                        t(blck_mat[[i]]))
    invM <- lapply(1:n,
                   function(i) diag(clsz[[i]]) - blck_mat[[i]] %*% projM[[i]])
    diagM <- as.matrix(Matrix::bdiag(invM))
    list(projM = projM, invM = invM, diagM = diagM)
}


estSigma2 <- function(data, pms, merr = NULL) {
    ## estimate beta here
    if(!is.null(data$tdep[[1]])) {
        merr_tdep <- lapply(1:data$n,
                            function(i) matrix(0, ncol = ncol(data$tdep[[i]]),
                                               nrow = ncol(data$tdep[[i]])))
        cidx_tdep <- 1
        if (!(is.null(data$cidx_tdep) | identical(integer(0), data$cidx_tdep))) {
            cidx_tdep <- data$cidx_tdep
            merr_tdep <- lapply(1:data$n, function(i) {
                merr_tdep[[i]][cidx_tdep, cidx_tdep] <- merr[[i]]
                merr_tdep[[i]]
            })
        }
    }
    if (is.null(data$tdep[[1]])) {
        beta_tdep <- NULL
    } else {
        beta_tdep_res <- nlm(estTimeDep,
                             p = matrix(rep(0, ncol(data$tdep[[1]]))),
                             rm = data$response, pm = pms$invM,
                             dm = data$tdep, dat_merr_tdep = data$merr_tdep,
                             merr_tdep = merr_tdep, n = data$n)
        beta_tdep <- matrix(beta_tdep_res$estimate)
    }
    p <- ncol(data$blck_mat[[1]])
    if(!is.null(beta_tdep)) {
        hsigma2_est <- lapply(1:data$n,
                              function(i) t(data$response[[i]] -
                                                data$tdep[[i]]%*% beta_tdep) %*%
                                  pms$invM[[i]] %*%
                                  (data$response[[i]] -
                                       data$tdep[[i]]%*% beta_tdep))
    } else {
        hsigma2_est <- lapply(1:data$n, function(i) t(data$response[[i]]) %*%
                                  pms$invM[[i]] %*% (data$response[[i]]))
    }
    c(Reduce(`+`, hsigma2_est) /(Reduce(`+`, data$clsz) - data$n * p))
}

## next one to correct for!
estSigma_ranef <- function(beta, beta_xt, lvl2est, dm,
                            hsigma_2, blck_mat, n, merr, mmX = NULL) {
    ## forget about the design matrix at this point?
    if (is.null(beta_xt)) {
        b.merr <- lapply(1:n, function(i) 0)
    } else {
        mmX[mmX != 0] <- beta_xt
        b.merr <- lapply(1:n, function(i) mmX %*% merr[[i]] %*%
                             t(mmX))
    }

    Sigma_est <- lapply(1:n, function(i) (lvl2est[[i]] - dm[[i]] %*% beta) %*%
                            t(lvl2est[[i]] - dm[[i]] %*% beta) -
                            c(hsigma_2) *
                            solve(t(blck_mat[[i]]) %*% blck_mat[[i]]) -
                            b.merr[[i]])
    #out <- Reduce(`+`, Sigma_est) / n
    out <- as.matrix(Matrix::nearPD(Reduce(`+`, Sigma_est) / n)$mat)
    out
}

estTimeDep <- function(par, dm, rm, pm, merr_tdep, n, dat_merr_tdep) {
    sum(sapply(1:n, function(i)
        t(rm[[i]] - dm[[i]] %*% par) %*% pm[[i]] %*% (rm[[i]] - dm[[i]] %*% par)))
}



estMod2 <- function(data.sub, data.pri, ...) {
    # if(is.null(merr)) merr <- lapply(1:data$n, function(i) 0)
    pms.sub <- PMcals(blck_mat = data.sub$blck_mat,
                      clsz = data.sub$clsz, n = data.sub$n)
    pms.pri <- PMcals(blck_mat = data.pri$blck_mat,
                      clsz = data.pri$clsz, n = data.pri$n)

    ## needs to keep this

    hsigma_e2 <- estSigma2(data = data.sub, pms = pms.sub)


    nX <- nrow(data.pri$tind[[1]])
    ncX <- ncol(data.pri$tind[[1]])
    merr <- lapply(1:data.sub$n, function(i)
        solve(t(data.sub$blck_mat[[i]]) %*%
                  data.sub$blck_mat[[i]])[2, 2] * c(hsigma_e2))

    mmX <- as.matrix(Reduce(`+`, data.pri$tind_shared))
    merr.sub.mat <- lapply(1:data.sub$n, function(i) {
        tmp <- matrix(0, nrow = ncX, ncol = ncX)
        mm <- data.pri$tind_shared[[i]]
        mm <- matrix(mm[mmX != 0])
        mm %*% t(mm) * merr[[i]]
    })
    hsigma_2 <- estSigma2(data = data.pri, pms = pms.pri, merr = merr)

    idx_beta_1r <- ncol(data.sub$tind[[1]])
    idx_beta_v <- idx_beta_1r + ifelse(is.null(ncol(data.sub$tdep[[1]])), 0,
                         ncol(data.sub$tdep[[1]]))

    idx_beta_u <- idx_beta_v + ncol(data.pri$tind[[1]])
    idx_beta_xt <- idx_beta_u + ncol(data.pri$tind_shared[[1]])
    idx_beta_z <- idx_beta_xt + ifelse(is.null(ncol(data.pri$tdep[[1]])), 0,
                                       ncol(data.pri$tdep[[1]]))


    res <- nleqslv::nleqslv(x = rep(1, idx_beta_z),
                            fn = jointSolve,
                            idx_beta_v =  idx_beta_v,
                            idx_beta_1r = idx_beta_1r,
                            idx_beta_z = idx_beta_z,
                            idx_beta_u = idx_beta_u,
                            idx_beta_xt = idx_beta_xt,
                            Y = data.pri$response,
                            W = data.sub$response,
                            matR = data.sub$tind, matV = data.sub$tdep,
                            matU = data.pri$tind, matX = data.pri$tind_shared,
                            matZ = data.pri$tdep,
                            pms.sub = pms.sub, pms.pri = pms.pri,
                            n = data.sub$n,
                            mmX = mmX, hsigma_e2 = hsigma_e2,
                            hsigma_2 = hsigma_2,
                            A_blck = data.pri$blck_mat,
                            D_blck = data.sub$blck_mat,
                            merr.sub.mat = merr.sub.mat,
                            Sigma_gamma_inv = diag(nrow(data.sub$tind[[1]])),
                            Sigma_b_inv = diag(nX)
    )
    par <- res$x
    ## estimate random effects
    beta_1r <- as.matrix(par[1:idx_beta_1r])
    beta_u <- as.matrix(par[(idx_beta_v + 1):idx_beta_u])
    beta_xt <- as.matrix(par[(idx_beta_u + 1):idx_beta_xt])

    ## the structure of Xhat should be given at the first place
    if (is.null(data.sub$tdep[[1]])) {
        hX <- lapply(1:data.sub$n, function(i) pms.sub$projM[[i]] %*%
                         data.sub$response[[i]])
    } else {
        beta_v <- as.matrix(par[(idx_beta_1r + 1):idx_beta_v])

        hX <- lapply(1:data.sub$n, function(i) pms.sub$projM[[i]] %*%
                         (data.sub$response[[i]] -
                              data.sub$tdep[[i]] %*% beta_v))
    }

    mathX <- lapply(1:data.sub$n, function(i) {
        data.pri$tind_shared[[i]] * hX[[i]][2, ]
    })

    if (is.null(data.pri$tdep[[1]])) {
        heta <- lapply(1:data.pri$n, function(i) pms.pri$projM[[i]] %*%
                           data.pri$response[[i]])
    } else {
        beta_z <- as.matrix(par[(idx_beta_xt + 1):idx_beta_z])
        heta <- lapply(1:data.pri$n, function(i) pms.pri$projM[[i]] %*%
                           (data.pri$response[[i]] -
                                data.pri$tdep[[i]] %*% beta_z))
    }


    hSigma_gamma <- estSigma_ranef(beta = beta_1r, beta_xt = NULL,
                                    lvl2est = hX, dm = data.sub$tind,
                                    hsigma_2 = hsigma_e2,
                                    blck_mat = data.sub$blck_mat,
                                    n = data.sub$n, mmX = NULL,
                                    merr = NULL)

    hSigma_b <- estSigma_ranef(beta = rbind(beta_u, beta_xt),
                    beta_xt = beta_xt,
                    lvl2est = heta,
                    dm = lapply(1:data.pri$n, function(i)
                        cbind(data.pri$tind[[i]], mathX[[i]])),
                    hsigma_2 = hsigma_2, mmX = mmX,
                    blck_mat = data.pri$blck_mat, n = data.pri$n,
                    merr = merr.sub.mat)

    res <- nleqslv::nleqslv(x = res$x,
                            fn = jointSolve,
                            idx_beta_v =  idx_beta_v,
                            idx_beta_1r = idx_beta_1r,
                            idx_beta_z = idx_beta_z,
                            idx_beta_u = idx_beta_u,
                            idx_beta_xt = idx_beta_xt,
                            Y = data.pri$response,
                            W = data.sub$response,
                            mmX = mmX, hsigma_e2 = hsigma_e2,
                            hsigma_2 = hsigma_2,
                            matR = data.sub$tind, matV = data.sub$tdep,
                            matU = data.pri$tind, matX = data.pri$tind_shared,
                            matZ = data.pri$tdep,
                            pms.sub = pms.sub, pms.pri = pms.pri,
                            n = data.sub$n,
                            A_blck = data.pri$blck_mat,
                            D_blck = data.sub$blck_mat,
                            merr.sub.mat = merr.sub.mat,
                            Sigma_gamma_inv = solve(hSigma_gamma),
                            Sigma_b_inv = solve(hSigma_b)
    )
    par <- res$x
    beta_1r <- as.matrix(par[1:idx_beta_1r])
    beta_u <- as.matrix(par[(idx_beta_v + 1):idx_beta_u])
    beta_xt <- as.matrix(par[(idx_beta_u + 1):idx_beta_xt])
    #
    ## the structure of Xhat should be given at the first place
    if (is.null(data.sub$tdep[[1]])) {
        hX <- lapply(1:data.sub$n, function(i) pms.sub$projM[[i]] %*%
                         data.sub$response[[i]])
        beta_v <- NULL
    } else {
        beta_v <- as.matrix(par[(idx_beta_1r + 1):idx_beta_v])
        hX <- lapply(1:data.sub$n, function(i) pms.sub$projM[[i]] %*%
                         (data.sub$response[[i]] -
                              data.sub$tdep[[i]] %*% beta_v))
    }

    mathX <- lapply(1:data.sub$n, function(i) {
        data.pri$tind_shared[[i]] * hX[[i]][2, ]
    })

    if (is.null(data.pri$tdep[[1]])) {
        heta <- lapply(1:data.pri$n, function(i) pms.pri$projM[[i]] %*%
                           data.pri$response[[i]])
        beta_z <- NULL
    } else {
        beta_z <- as.matrix(par[(idx_beta_xt + 1):idx_beta_z])
        heta <- lapply(1:data.pri$n, function(i) pms.pri$projM[[i]] %*%
                           (data.pri$response[[i]] -
                                data.pri$tdep[[i]] %*% beta_z))
    }


    infs <- cMat(beta_1r = beta_1r, beta_u = beta_u, beta_xt = beta_xt,
                 beta_v = beta_v, beta_z = beta_z,
                 matR = data.sub$tind, matV = data.sub$tdep,
                 matU = data.pri$tind, matX = data.pri$tind_shared,
                 matZ = data.pri$tdep, Y = data.pri$response,
                 W = data.sub$response,
                 pms.sub = pms.sub, pms.pri = pms.pri,
                 n = data.sub$n, mmX = mmX,
                 hsigma_e2 = hsigma_e2,
                 hsigma_2 = hsigma_2,
                 A_blck = data.pri$blck_mat,
                 D_blck = data.sub$blck_mat,
                 merr.sub.mat = merr.sub.mat,
                 Sigma_gamma_inv = solve(hSigma_gamma),
                 Sigma_b_inv = solve(hSigma_b),
                 heta = heta,
                 hX = hX, mathX = mathX)
    ses <- sqrt(diag(infs[[2]]) / data.sub$n)
    hSigma_gamma <- estSigma_ranef(beta = beta_1r, beta_xt = NULL,
                                    lvl2est = hX, dm = data.sub$tind,
                                    hsigma_2 = hsigma_e2,
                                    blck_mat = data.sub$blck_mat,
                                    n = data.sub$n, mmX = NULL,
                                    merr = NULL)

    hSigma_b <- estSigma_ranef(beta = rbind(beta_u, beta_xt),
                    beta_xt = beta_xt,
                    lvl2est = heta,
                    dm = lapply(1:data.pri$n, function(i)
                        cbind(data.pri$tind[[i]], mathX[[i]])),
                    hsigma_2 = hsigma_2, mmX = mmX,
                    blck_mat = data.pri$blck_mat, n = data.pri$n,
                    merr = merr.sub.mat)

    list(beta_1r = beta_1r, beta_u = beta_u, beta_xt = beta_xt,
         beta_v = beta_v, beta_z = beta_z,
         hsigma_2 = hsigma_2, hsigma_e2 = hsigma_e2,
         convergence = res$termcd,
         hX = hX, heta = heta, ses = ses,
         Sigma_gm = hSigma_gamma,
         Sigma_b = hSigma_b)
}

