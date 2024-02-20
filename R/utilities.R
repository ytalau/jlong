#' @importFrom stats as.formula model.matrix model.response model.extract
#' model.frame terms pchisq printCoefmat rnorm glm coef
#' @export
extract.formula <- function(f.obj) {

    stopifnot(inherits(f.obj, "formula"))

    rhs <- if(length(f.obj) > 2) f.obj[[3L]] else f.obj[[2L]]
    lhs <- if(length(f.obj) > 2) f.obj[[2L]] else NULL

    fixed <- if(length(rhs) > 2) rhs[[2L]] else NULL
    random <- if(length(rhs) > 2) rhs[[3L]][[2L]] else rhs[[2L]][[2L]]

    random_formula <- random[[2]]
    id.var <- random[[3]]

    list(lhs = lhs, fixed = fixed, random = random_formula, id.var = id.var)
}


## for fixed effect divided into time-depdnent or time-independent
#' @export
split.fixed <- function(model.data) {
    ## for each data frame
    tind <- NULL
    tdep <- NULL
    idx <- apply(model.data, 2, sd) == 0
    if (sum(idx) == 0) {
        tdep <- model.data
    } else if (sum(idx) == ncol(model.data)) {
        ## redefine....matR
        tind <- unique(model.data)
    } else {
        tdep  <- model.data[, !idx, drop = F]
        tind <- unique(model.data[, idx, drop = F])
    }
    list(tind = tind, tdep = tdep)
}

#' @export
getClsz <- function(formula, data) {

    fm <- as.formula(formula)
    parts <- extract.formula(fm)

    sdat <- split(data, f = as.formula(paste("~", as.character(parts$id.var))))

    sapply(1:length(sdat), function(i) nrow(sdat[[i]]))
}

#' @importFrom stats as.formula model.matrix model.response model.extract
#' model.frame terms pchisq printCoefmat rnorm glm coef
#' @importFrom data.table is.data.table
#' @export
data.process2 <- function(formula, data, var.shared = NULL) {

    fm <- as.formula(formula)
    parts <- extract.formula(fm)
    ## this comes from the package data table...
    if(is.data.table(data)) {
        sdat <- split(data, by =as.character(parts$id.var))
    } else {
        sdat <- base::split(data,
                            f = as.formula(paste("~", as.character(parts$id.var))))
    }

    n <- length(sdat)
    random.vars <- as.character(all.vars(parts$random))

    fixed.fm <- as.formula(paste0(parts$lhs, "~",
                                  paste(parts$fixed, collapse = "+")))

    mf <- lapply(1:n, function(i) model.frame(fixed.fm, data = sdat[[i]]))

    response <- lapply(1:n, function(i) as.matrix(model.response(mf[[i]])))

    mdat <- lapply(1:n, function(i) model.matrix(fixed.fm, data = sdat[[i]]))
    var.names <- colnames(mdat[[1]])
    mdat.slps <- NULL
    blck_mat <- lapply(1:n, function(i) cbind(rep(1, nrow(mdat[[i]])),
                                              mdat[[i]][, random.vars]))
    fmatrix <- attr(terms(fixed.fm), "factors")
    assign_id <- attr(mdat[[1]], "assign")

    ## work on this part so that we can extract an additional design matrix
    if (!identical(random.vars, character(0))) {
        if(length(random.vars) > 1) {
            ## this would guarantee it returns a list; under development
            idx <- apply(attr(terms(fixed.fm),
                              "factors")[random.vars, ] == 1, 1,
                         which, simplify = FALSE)

        } else {
            ## grab the row of random effect variables
            idx <- fmatrix[random.vars, ]
            ## grab the covariates "interact" with the random effect
            name_idx <- names(idx)[which(idx == 1)]
            selectidx <- c(attr(terms(fixed.fm), "intercept") == 1, !idx)
            selectidx_names <- names(selectidx)[selectidx == F]
            tmpidx <- apply(as.matrix(fmatrix[, setdiff(name_idx, random.vars)]) == 1,
                            1, any)
            colidx <- which(tmpidx)
            colidx2 <- which(colnames(mdat[[1]]) %in% random.vars)
            slp.idx <- setdiff(colidx, colidx2)
            slp.vars <- setdiff(names(colidx), random.vars)
            if(identical(slp.idx, integer(0))) {
                mdat.slps <- lapply(1:length(mdat), function(i)
                    matrix(1, dimnames = list(NULL, c(random.vars))))
            } else {
                mdat.names <- var.names[which(idx == 1) + 1L]
                mdat.slps <- lapply(1:length(mdat), function(i)
                    matrix(cbind(1, unique(mdat[[i]][, slp.idx])),
                           nrow = 1, dimnames = list(NULL, mdat.names)))
            }

            mdat <- lapply(1:n, function(i) mdat[[i]][, selectidx])
        }
    }

    design <- lapply(1:n, function(i) split.fixed(mdat[[i]]))

    tdep <- lapply(1:n, function(i) design[[i]][["tdep"]])
    idx_shared <- which(fmatrix[rownames(fmatrix) == var.shared, ] == 1)
    add.var.shared <- var.names[assign_id %in% idx_shared]

    sf_names <- setdiff(colnames(design[[1]][["tind"]]), add.var.shared)
    int_names <- intersect(colnames(design[[1]][["tind"]]), add.var.shared)

    tind <- lapply(1:n, function(i)
        design[[i]][["tind"]][, sf_names, drop = F])
    tind_shared <- lapply(1:n, function(i)
        design[[i]][["tind"]][, int_names, drop = F])

    if (!is.null(mdat.slps)) {
        if (is.null(var.shared)) {
            tind_shared <- NULL
            tind <- lapply(1:n, function(i) {
                tmp <- as.matrix(Matrix::bdiag(tind[[i]], mdat.slps[[i]]))
                colnames(tmp) <- c(colnames(tind[[i]]), colnames(mdat.slps[[i]]))
                tmp
            })
        } else {
            if (ncol(mdat.slps[[1]]) == 1) {
                sf_idx <- TRUE
                int_idx <- FALSE
            } else {
                sf_idx <- c(TRUE, !slp.vars %in% var.shared)
                int_idx <- c(FALSE, slp.vars %in% var.shared)
            }

            mdat.slps_shared <- lapply(1:n, function(i)
                mdat.slps[[i]][, which(int_idx), drop = F])
            mdat.slps <- lapply(1:n, function(i)
                mdat.slps[[i]][, which(sf_idx), drop = F])
            tind <- lapply(1:n, function(i) {
                tmp <- as.matrix(Matrix::bdiag(tind[[i]], mdat.slps[[i]]))
                colnames(tmp) <- c(colnames(tind[[i]]), colnames(mdat.slps[[i]]))
                tmp
            })
            tind_shared <- lapply(1:n, function(i) {
                tmp <- as.matrix(Matrix::bdiag(tind_shared[[i]],
                                               mdat.slps_shared[[i]]))
                colnames(tmp) <- c(colnames(tind_shared[[i]]),
                                   colnames(mdat.slps_shared[[i]]))
                tmp
            })
        }

}


    if(is.null(var.shared)) {
        cidx_tind <- 1
        cidx_tdep <- NULL
        merr_tdep <- NULL
    } else {
        cidx_tind <- which(grepl(var.shared, colnames(tind[[1]])))
        if(is.null(tdep[[1]])){
            cidx_tdep <- NULL
            merr_tdep <- NULL
        } else {
            cidx_tdep <- which(grepl(var.shared, colnames(tdep[[1]])))
            merr_tdep <- lapply(1:n, function(i)
                matrix(0,nrow = nrow(tdep[[i]]), ncol = ncol(tdep[[i]])))
        }

    }
    if(!(is.null(cidx_tdep) |identical(integer(0), cidx_tdep))) {
        ## TODO: redesign for multiple interaction terms.
        merr_tdep_names <- setdiff(
            names(which(fmatrix[, colnames(tdep[[1]])[cidx_tdep]] == 1)),
            var.shared)
        merr_tdep <- lapply(1:n, function(i)
            mdat[[i]][, merr_tdep_names, drop = F])

    }
    clsz <- lapply(1:n, function(i) nrow(blck_mat[[i]]))

    list(response = response, tind = tind, tdep = tdep,
         tind_shared = tind_shared, random.vars = random.vars,
         blck_mat = blck_mat,
         n = n, clsz = clsz, cidx_tdep = cidx_tdep,
         cidx_tind = cidx_tind,
         merr_tdep = merr_tdep)
}
