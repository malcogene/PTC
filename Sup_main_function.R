
setClass(
  "binary",
  slot=c(
    trueclass  = "vector",
    predclass = "vector",
    predprob = "numeric"
  ),
  validity=function(obj)
  {
    if(length(obj@trueclass) != length(obj@predclass) | length(obj@trueclass) != length(obj@predprob) )
      return("ERROR: length NOT equal!")
    if( length(obj@predprob[obj@predprob < 0 | obj@predprob >1 ]) >=1 )
      return("Slot `predprob` must be a closed inteval [0,1]")
    return(TRUE)
  }
)

setGeneric("brier", function(obj) standardGeneric("brier"))
setMethod("brier", "binary", function(obj) mean((obj@trueclass - obj@predprob)^2))

setGeneric("APRFMscore", function(obj) standardGeneric("APRFMscore"))
setMethod("APRFMscore", "binary", function(obj) {
  tp <- sum( (obj@predclass==1) & (obj@trueclass == 1) )
  tn <- sum( (obj@predclass==0) & (obj@trueclass == 0) )
  fp <- sum( (obj@predclass==1) & (obj@trueclass == 0) )
  fn <- sum( (obj@predclass==0) & (obj@trueclass == 1) )
  ACC <-  mean(obj@predclass == obj@trueclass)
  precision <- sum(obj@predclass & obj@trueclass) / sum(obj@predclass)
  recall <- sum(obj@predclass & obj@trueclass) / sum(obj@trueclass)
  F1 <- 2 * precision * recall / (precision + recall)
  MCC <- (tp * tn - fp * fn) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  
  cat(sprintf("accuracy : %f" , ACC))
  cat("\n")
  cat(sprintf("precision : %f" , precision))
  cat("\n")
  cat(sprintf("recall : %f", recall))
  cat("\n")
  cat(sprintf("F-score (F1 score) is the harmonic mean of precision and recall: %f", F1))
  cat("\n")
  cat(sprintf("MCC : %f" , MCC))
  cat("\n")
  cat("\n")
  return(list(ACC=ACC, precision=precision, recall=recall, F1=F1, MCC=MCC ))
})

setGeneric("F.beta.score", function(obj) standardGeneric("F.beta.score"))
setMethod("F.beta.score", "binary", function(obj) {
  precision <- sum(obj@predclass & obj@trueclass) / sum(obj@predclass)
  recall <- sum(obj@predclass & obj@trueclass) / sum(obj@trueclass)
  f.beta<- function(beta) { (beta^2 +1) * ( ( precision * recall) / ( ( (beta^2)*precision)  + recall)    ) }
  f.beta
})


setGeneric("AUC", function(obj) standardGeneric("AUC"))
setMethod("AUC", "binary", function(obj) { a <- pROC::roc(obj@trueclass, obj@predprob,  direction = "<")
a <- as.numeric(gsub(".*: ", '', a$auc))
a } )



browseGEO = function(GEOid) {
  url <- paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", GEOid)
  browseURL(url)
}


browseGEOiter = function(x) {
  for(i in 1:length(x)) {
    url <- paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", x[i])
    browseURL(url) }
}


.epi.tests<- function (dat, conf.level = 0.95) 
{
    elements <- list()
    elements <- within(elements, {
        N. <- 1 - ((1 - conf.level)/2)
        z <- qnorm(N., mean = 0, sd = 1)
        .funincrisk <- function(cdat, conf.level) {
            N. <- 1 - ((1 - conf.level)/2)
            a <- cdat[, 1]
            n <- cdat[, 2]
            b <- n - a
            p <- a/n
            a. <- ifelse(a == 0, a + 1, a)
            b. <- ifelse(b == 0, b + 1, b)
            low <- a./(a. + (b. + 1) * (1/qf(1 - N., 2 * a., 
                2 * b. + 2)))
            up <- (a. + 1)/(a. + 1 + b./(1/qf(1 - N., 2 * b., 
                2 * a. + 2)))
            low <- ifelse(a == 0, 0, low)
            up <- ifelse(a == n, 1, up)
            rval <- data.frame(est = p, lower = low, upper = up)
            rval
        }
        a <- dat[1]
        b <- dat[3]
        c <- dat[2]
        d <- dat[4]
        M1 <- a + c
        M0 <- b + d
        N1 <- a + b
        N0 <- c + d
        total <- a + b + c + d
        tdat <- as.matrix(cbind(M1, total))
        trval <- .funincrisk(tdat, conf.level)
        tp <- trval$est
        tp.low <- trval$lower
        tp.up <- trval$upper
        tprev <- data.frame(est = tp, lower = tp.low, upper = tp.up)
        tdat <- as.matrix(cbind(N1, total))
        trval <- .funincrisk(tdat, conf.level)
        ap <- trval$est
        ap.low <- trval$lower
        ap.up <- trval$upper
        aprev <- data.frame(est = ap, lower = ap.low, upper = ap.up)
        tdat <- as.matrix(cbind(a, M1))
        trval <- .funincrisk(tdat, conf.level)
        se <- trval$est
        se.low <- trval$lower
        se.up <- trval$upper
        sensitivity <- data.frame(est = se, lower = se.low, upper = se.up)
        tdat <- as.matrix(cbind(d, M0))
        trval <- .funincrisk(tdat, conf.level)
        sp <- trval$est
        sp.low <- trval$lower
        sp.up <- trval$upper
        specificity <- data.frame(est = sp, lower = sp.low, upper = sp.up)
        tdat <- as.matrix(cbind(a, N1))
        trval <- .funincrisk(tdat, conf.level)
        ppv <- trval$est
        ppv.low <- trval$lower
        ppv.up <- trval$upper
        pv.positive <- data.frame(est = ppv, lower = ppv.low, 
            upper = ppv.up)
        tdat <- as.matrix(cbind(d, N0))
        trval <- .funincrisk(tdat, conf.level)
        npv <- trval$est
        npv.low <- trval$lower
        npv.up <- trval$upper
        pv.negative <- data.frame(est = npv, lower = npv.low, 
            upper = npv.up)
        lrpos <- (a/M1)/(1 - (d/M0))
        lrpos.low <- exp(log(lrpos) - z * sqrt((1 - se)/(M1 * 
            se) + (sp)/(M0 * (1 - sp))))
        lrpos.up <- exp(log(lrpos) + z * sqrt((1 - se)/(M1 * 
            se) + (sp)/(M0 * (1 - sp))))
        lr.positive <- data.frame(est = lrpos, lower = lrpos.low, 
            upper = lrpos.up)
        lrneg <- (1 - (a/M1))/(d/M0)
        lrneg.low <- exp(log(lrneg) - z * sqrt((se)/(M1 * (1 - 
            se)) + (1 - sp)/(M0 * (sp))))
        lrneg.up <- exp(log(lrneg) + z * sqrt((se)/(M1 * (1 - 
            se)) + (1 - sp)/(M0 * (sp))))
        lr.negative <- data.frame(est = lrneg, lower = lrneg.low, 
            upper = lrneg.up)
        tdat <- as.matrix(cbind((a + d), total))
        trval <- .funincrisk(tdat, conf.level)
        da <- trval$est
        da.low <- trval$lower
        da.up <- trval$upper
        diag.acc <- data.frame(est = da, lower = da.low, upper = da.up)
        dOR.p <- (a * d)/(b * c)
        lndOR <- log(dOR.p)
        lndOR.var <- 1/a + 1/b + 1/c + 1/d
        lndOR.se <- sqrt(1/a + 1/b + 1/c + 1/d)
        lndOR.l <- lndOR - (z * lndOR.se)
        lndOR.u <- lndOR + (z * lndOR.se)
        dOR.se <- exp(lndOR.se)
        dOR.low <- exp(lndOR.l)
        dOR.up <- exp(lndOR.u)
        diag.or <- data.frame(est = dOR.p, lower = dOR.low, upper = dOR.up)
        ndx <- 1/(se - (1 - sp))
        ndx.1 <- 1/(se.low - (1 - sp.low))
        ndx.2 <- 1/(se.up - (1 - sp.up))
        ndx.low <- min(ndx.1, ndx.2)
        ndx.up <- max(ndx.1, ndx.2)
        nnd <- data.frame(est = ndx, lower = ndx.low, upper = ndx.up)
        c.p <- se - (1 - sp)
        c.1 <- se.low - (1 - sp.low)
        c.2 <- se.up - (1 - sp.up)
        c.low <- min(c.1, c.2)
        c.up <- max(c.1, c.2)
        youden <- data.frame(est = c.p, lower = c.low, upper = c.up)
    })
    rval <- list(aprev = elements$aprev, tprev = elements$tprev, 
        se = elements$sensitivity, sp = elements$specificity, 
        diag.acc = elements$diag.acc, diag.or = elements$diag.or, 
        nnd = elements$nnd, youden = elements$youden, ppv = elements$pv.positive, 
        npv = elements$pv.negative, plr = elements$lr.positive, 
        nlr = elements$lr.negative)
    r1 <- with(elements, c(a, b, N1))
    r2 <- with(elements, c(c, d, N0))
    r3 <- with(elements, c(M1, M0, M0 + M1))
    tab <- as.data.frame(rbind(r1, r2, r3))
    colnames(tab) <- c("   Disease +", "   Disease -", "     Total")
    rownames(tab) <- c("Test +", "Test -", "Total")
    tab <- format.data.frame(tab, digits = 3, justify = "right")
    out <- list(conf.level = conf.level, elements = elements, 
        rval = rval, tab = tab)
    class(out) <- "epi.tests"
    return(out)
}


summary.epi.tests <- function (obj) 
{
    names(obj$rval) <- c("Apparent Prevalence", "True Prevalence", "Sensitivity (Recall / TPR)", "Specificity", "Accuracy", "Diagnostic Odds Ratio (OR)", "Number Needed to Diagnose (NND)", "Youden's index", "Precision (PPV)", "Negative Predictive Value (NPV)", "Likelihood Ratio Positive (LR+)", "Likelihood Ratio Negative (LR-)") 
    out <- do.call(rbind, obj$rval)
    return(out)
}




wrap.PDS =
  function (data,
            allgenes,
            syms,
            pathwaynames,
            normals = NULL,
            ranks = NULL,
            attempts = 100,
            maximize_stability = TRUE,
            logfile = "",
            samplings = NULL,
            min_exp = 4,
            min_std = 0.4)
  {
    cat(
      file = logfile,
      append = FALSE,
      "robust_score_bydist. min_exp=",
      min_exp,
      ", min_std=",
      min_std,
      "\n"
    )
    
    if (!is.null(seed <- getOption("seed"))) set.seed(seed)  # .w.
    
    data[data < min_exp] = min_exp
    n <- ncol(data)
    if (is.null(normals)) {
      normals <- rep(TRUE, n)
      start <- "by pca"
    }
    else {
      start <- "by ranks"
    }
    if (is.null(ranks))
      ranks <- !normals
    ranks <- rank(ranks)
    if ((length(normals) != n) || (length(ranks) != n)) {
      stop("invalid dimentions")
    }
    l <- length(syms)
    nn <- sum(normals)
    m <- floor(0.8 * (n - nn)) + nn
    if (is.null(samplings)) {
      samplings <- matrix(0, attempts, m)
      w <- which(!normals)
      for (a in 1:attempts) {
        samplings[a,] <- sort(c(w[sample(n - nn, m - nn)],
                                which(normals)))
      }
    }
    s <- NULL
    ind <- NULL
    for (i in 1:l) {
      pathway <- syms[[i]]
      pathwayindata <-
        pathifier:::.getpathway(pathway, allgenes, data)
      k1 = sum(pathwayindata$isin)
      if (k1 < 3) {
        si <- NULL
        cat(file = logfile,
            append = TRUE,
            "skipping pathway ",
            i,
            " k1=",
            k1,
            "\n")
      }
      else {
        x <- pathwayindata$x
        pathway <- pathway[pathwayindata$isin]
        xm <- colMeans(x[normals,])
        xs <- apply(x[normals,], 2, sd)
        xs[xs < min_std] = min_std
        if (0 %in% xs) {
          si <- NULL
          cat(file = logfile,
              append = TRUE,
              "skipping pathway ",
              i,
              " (0 in xs)\n")
        }
        else {
          z <- (x - matrix(rep(xm, each = n), nrow = n)) / (matrix(rep(xs,
                                                                       each = n), nrow = n))
          t <- prcomp(z)
          k2 = max(sum(t$sdev > 1.1), 4)
          k2 = min(k2, k1, 0.75 * dim(x)[1], sum(t$sdev >
                                                   0.25))
          if (k2 < 3) {
            si <- NULL
            cat(file = logfile,
                append = TRUE,
                "skipping pathway ",
                i,
                " k2=",
                k2,
                "\n")
          }
          else {
            pca <- t$x[, 1:k2]
            res <-
              pathifier:::.score_all_pathways_helper(pca,
                                                     ranks,
                                                     samplings,
                                                     i,
                                                     attempts,
                                                     maximize_stability,
                                                     logfile,
                                                     start = start)
            if (is.null(res)) {
              si <- NULL
              cat(file = logfile,
                  append = TRUE,
                  "skipping pathway ",
                  i,
                  "\n")
            }
            else {
              ind <- c(ind, i)
              si <- list(
                res$score,
                pathway,
                res$sig,
                res$origsig,
                res$k,
                res$thecurve$s,
                res$thecurve$tag,
                res$z,
                res$isin,
                xm,
                xs,
                t$center,
                t$rotation,
                k2,
                pca
              )
            }
          }
        }
      }
      s <- rbind(s, si)
    }
    cat(
      file = logfile,
      append = TRUE,
      length(ind),
      "pathways processed with start=",
      start,
      "\n"
    )
    rownames(s) <- pathwaynames[ind]
    list(
      scores = s[, 1],
      genesinpathway = s[, 2],
      newmeanstd = s[,
                     3],
      origmeanstd = s[, 4],
      pathwaysize = s[, 5],
      curves = s[,
                 6],
      curves_order = s[, 7],
      z = s[, 8],
      compin = s[, 9],
      xm = s[, 10],
      xs = s[, 11],
      center = s[, 12],
      rot = s[,
              13],
      pctaken = s[, 14],
      pca = s[, 15],
      samplings = samplings,
      sucess = ind,
      logfile = logfile
    )
  }


GA.solution <-
  function(x,
           y,
           weights,
           foldid,
           family,
           popSize,
           maxiter,
           seed,
           ...) {
    e <- new.env()
    e$GA.cvFitList = NULL
    e$GA.sol =NULL
    e$Crude3Dmatrix = NULL


    GA_f = function(x, y, weights, foldid, family, alpha, ...) {
      set.seed(seed)
      fit = glmnet::cv.glmnet(
        x = x,
        y = y,
        weights = weights,
        foldid = foldid,
        family = family,
        alpha = alpha,
        ...

      )

      print(fit)
      print(fit$lambda.min)
      print(fit$cvm)
      print(min(fit$cvm))
      e$Crude3Dmatrix <-
        rbind(e$Crude3Dmatrix, c(alpha, fit$lambda.min, min(fit$cvm)))
      e$GA.cvFitList <- c(e$GA.cvFitList, list(fit))
      return(min(fit$cvm))
    }
    monitor <- function(obj) {
      alpha_x = obj@population
      fitness = obj@fitness
      plot(alpha_x, fitness, main = paste("Generation =", obj@iter))
      rug(alpha_x, col = 2)
      a = cbind(alpha_x, fitness)
      print(a)
    }
    fitness <-
      function(alpha, ...)
        - GA_f(x, y, weights, foldid, family, alpha, ...)
    GA <<-
      GA::ga(
        type = "real-valued",
        fitness = fitness,
        min = 0,
        max = 1,
        popSize = popSize,
        maxiter = maxiter,
        seed = seed,
        monitor = monitor,
        ...
      )
    plot(GA)
    summary(GA)
    e$GA.sol <- GA@solution
    GA.results <<-list(e$GA.cvFitList, e$GA.sol, e$Crude3Dmatrix)
    return(GA.results)
  }


cv.merged <-
  function(ematMerged,
           sampleMetadata,
           weights,
           alphas,
           nFolds = 10,
           foldid = NA,
           nRepeats = 3,
           yName = 'class',
           clinVarColnames = NA,
           GlobalOp = NA,
           seed = 1234,
           popSize=30,
           maxiter=5,
           type.min="lambda.min",
           ...) {
    args = list(...)
    if (!is.null(args[['family']]) & args[['family']] == 'cox') {
      y = as.matrix(sampleMetadata[colnames(ematMerged), yName])
      colnames(y) = c('time', 'status')
    } else {
      y = sampleMetadata[colnames(ematMerged), yName]
    }

    if (is.na(clinVarColnames[1])) {
      x = scale(t(ematMerged), center = TRUE, scale = FALSE)
    } else {
      clinVarTmp = data.frame(lapply(sampleMetadata[colnames(ematMerged), clinVarColnames], factor))
      clinDummy = model.matrix( ~ 0 + ., data = clinVarTmp)
      x = cbind(scale(t(ematMerged), center = TRUE, scale = FALSE), clinDummy)
    }


    # if (is.na(foldid[1]) & is.na(GlobalOp) ) {
    if (is.na(GlobalOp) ) {
      if (is.na(foldid[1])) {

        cvFitList = list()
        for (ii in 1:nRepeats) {
          foldid = sample(rep(seq(nFolds), length = ncol(ematMerged)))
          cvFitList[[ii]] = foreach(alpha = alphas) %do% {
            cv.glmnet(
              x,
              y,
              weights = weights[colnames(ematMerged)],
              foldid = foldid,
              alpha = alpha,
              ...
            )
          }
        }
        return(cvFitList) } else {
          cvFitList = foreach(alpha = alphas) %do% {
            cv.glmnet(
              x,
              y,
              weights = weights[colnames(ematMerged)],
              foldid = foldid[colnames(ematMerged)],
              alpha = alpha,
              ...
            )
          }
          return(cvFitList) }

    } else if (GlobalOp == "epsgo") {

      if (is.na(foldid[1])) foldid = sample(rep(seq(nFolds), length = ncol(ematMerged)))

      y.classes <- y
      bounds <- t(data.frame(alpha=c(0, 1)))
      colnames(bounds)<-c("lower","upper")

      epsgo.cvFitList <- c060::epsgo(Q.func=".w.tune.glmnet.interval",
                                     parms.coding="none",
                                     seed = seed,
                                     show="none",
                                     bounds = bounds,
                                     x = x, y = y.classes,
                                     foldid=foldid,
                                     weights =weights,
                                     type.min=type.min,
                                     ...
      )

      return(epsgo.cvFitList)

    } else if (GlobalOp == "GA") {

      if (is.na(foldid[1])) foldid = sample(rep(seq(nFolds), length = ncol(ematMerged)))

      GA.results <- GA.solution (
        x = x,
        y = y,
        weights = weights,
        foldid = foldid,
        seed =seed,
        popSize=popSize,
        maxiter=maxiter,
        ...
      )
      return(GA.results)

    } else { }
  }



.w.tune.glmnet.interval <-
  function (parms, x, y, weights, offset = NULL, lambda = NULL,
            type.measure = c("mse", "deviance", "class", "auc", "mae"),
            seed = 12345, nfolds = 10, foldid = foldid, grouped = TRUE,
            type.min = c("lambda.min", "lambda.1se"), family, verbose = FALSE,
            ...)
  {
    alpha <- parms[1]
    names(alpha) <- NULL
    if (verbose)
      print(paste("alpha=", alpha))
    set.seed(seed)
    cv <- cv.glmnet(x = x, y = y, family = family, alpha = alpha,
                    offset = NULL, lambda = NULL, type.measure = type.measure, weights=weights, nfolds = nfolds, foldid = foldid, grouped = grouped, keep=T, intercept=T, standardize = F ) 
    opt.lambda <- ifelse(type.min == "lambda.min", cv$lambda.min,
                         cv$lambda.1se)
    q.val <- cv$cvm[which(cv$lambda == opt.lambda)]
    fit <- glmnet(x = x, y = y, family = family, alpha = alpha,
                  lambda = opt.lambda)
    ret <- list(q.val = q.val, model = list(alpha = alpha, lambda = opt.lambda,
                                            nfolds = nfolds, cvreg = cv, fit = fit))
    return(ret)
  }






PDS.predictValidationData <-
  function(ematList, studyMetadata, sampleMetadata, discoverySampleNames, classesTrain,
           alpha, lambda, weights, batchColname='study', covariateName=NA, className='class',
           familyName='binomial', predType='response', intercept=TRUE) {

    discoveryStudyNames = studyMetadata[studyMetadata[,'discovery'], 'study']
    validationStudyNames = studyMetadata[studyMetadata[,'validation'], 'study']

    predsPDSList = foreach(validationStudyName=validationStudyNames) %do% {
      idxValidation = sampleMetadata[,'study']==validationStudyName & (sampleMetadata[,className] %in% classesTrain)
      validationSampleNames = sampleMetadata[idxValidation, 'sample']

      ematListNow = ematList[c(discoveryStudyNames, validationStudyName)]
      ematMergedDiscVal = mergeStudyData(ematListNow, sampleMetadata, batchColname=batchColname, covariateName=covariateName)
      ematMergedDisc = ematMergedDiscVal[,discoverySampleNames]

      exp <- ematMergedDiscVal
      exp.merged = exp [which(rownames(exp) %in% gene_ids),] ; dim(exp)
      nor <- rep(F, dim(ematMergedDiscVal)[2])
      nor[which(colnames(ematMergedDisc) %in% sampleMetadata$sample[which(sampleMetadata$class == "C" )])] <- T
      normal.merged <-nor

      PDS.result <-wrap.PDS(exp.merged, rownames(exp.merged), g.CPDB.KEGG$entrez_gene_ids, g.CPDB.KEGG$pathway, normal.merged, attempts = 2, min_exp = -3, min_std = 0.4)

      # PDS.result <- foreach(i=1:length(rownames(g.CPDB.KEGG)), .combine = list, .multicombine = T)  %dopar%  (wrap.PDS(exp.merged, rownames(exp.merged), g.CPDB.KEGG$entrez_gene_ids[i], g.CPDB.KEGG$pathway[i], normal.merged, attempts = 2, min_exp = -3, min_std = 0.4)))

      PDSmatrix <- do.call(rbind.data.frame, PDS.result$scores)
      print(dim(PDSmatrix))
      d <-as.matrix(PDSmatrix)
      colnames(d)<-colnames(ematMergedDiscVal)

      fitResult = glmnet(t(d[,discoverySampleNames]), sampleMetadata[discoverySampleNames, className], alpha=alpha, lambda=lambda, weights=weights[discoverySampleNames], family=familyName, standardize=FALSE, intercept=intercept)
      preds = predict(fitResult, newx=t(d[,validationSampleNames]), s=lambda, type=predType)
    }
    names(predsPDSList) = validationStudyNames
    return(predsPDSList)}





theme_complete_bw =
  function (base_size = 11,
            base_family = "")
  {
    theme_grey(base_size = base_size, base_family = base_family) %+replace%
      theme(
        panel.background = element_rect(fill = "white", colour = NA),
        panel.border = element_rect(fill = NA,colour = "grey20"),
        panel.grid.major = element_line(colour = "grey92"),
        panel.grid.minor = element_line(colour = "grey92",size = 0.25),
        strip.background = element_rect(fill = "grey100", colour = "grey20"),
        legend.key = element_rect(fill = "white", colour = NA),
        complete = TRUE
      )
  }




wrapper.MetaDE.filter =
  function (x, DelPerc) {
    Mean.rank <- sapply(x, function(z) rank(apply(z[[1]], 1,
                                                  mean, na.rm = T)))
    mean.r.mv <- rowMeans(Mean.rank, na.rm = T)
    mean.r.mv <- mean.r.mv[order(mean.r.mv, decreasing = T)]
    Gene_mv <- names(mean.r.mv)[which(mean.r.mv > quantile(mean.r.mv,
                                                           DelPerc[1]))]
    Gene_mv = Gene_mv[-which(Gene_mv=="")]  # wrapper
    SD.rank <- sapply(x, function(z) rank(apply(z[[1]][Gene_mv,], 1, mean, na.rm = T)))
    mean.r.sd <- rowMeans(SD.rank, na.rm = T)
    mean.r.sd <- mean.r.sd[order(mean.r.sd, decreasing = T)]
    final.genes <- names(mean.r.sd)[which(mean.r.sd > quantile(mean.r.sd,
                                                               DelPerc[2]))]
    K <- length(x)
    for (i in 1:K) {
      x[[i]][[1]] <- x[[i]][[1]][final.genes, ]
    }
    return(x)
  }



wrapper.posterior.mean=
  function (data, varname, nsamp, permute = 0)
  {
    # data=Data; varname="MSIstatus"; nsamp=2; permute = 0

    N <- length(GEDM(data))
    dataZ <- function(data, varname, nsamp) {
      N <- length(GEDM(data))
      dum <- GEDM(data)
      cl <- clinical(data)
      for (i in 1:N) {
        group <- grep(varname, names(clinical(data)[[i]]))
        lev <- levels(as.factor(as.numeric(as.factor(clinical(data)[[i]][,
                                                                         group]))))
        dum[[i]] <- dum[[i]][, c(sample(which(as.numeric(as.factor(clinical(data)[[i]][,
                                                                                       group])) == lev[1]), nsamp), sample(which(as.numeric(as.factor(clinical(data)[[i]][,
                                                                                                                                                                          group])) == lev[2]), nsamp))]
        vars = list("var1", "var2")
        tmp = data.frame(cl = c(rep(0, nsamp), rep(1, nsamp)))
        #  names = colnames(dum[[i]]))
        # names(tmp) <- c(varname, "names")
        # names(vars) = names(tmp)
        cl[[i]] = tmp
      }
      newdata <- new("MetaArray", GEDM = dum, clinical = cl,
                     datanames = datanames(data))
      return(newdata)
    }
    sampl.data <- dataZ(data, varname, nsamp)
    sampl.data <- as.list(sampl.data)
    arg <- list()
    for (i in c(1:N)) {
      cl <- sampl.data[[i]]$clinical
      row.names(cl) <- cl$names
      ex <- as.matrix(sampl.data[[i]]$GEDM)
      arg[[i]] <- new("ExpressionSet", exprs = ex, phenoData = new("AnnotatedDataFrame",
                                                                   data = cl))
    }
    merged <- mergeExprs2(arg, datanames(data))
    z.stat <- metaArray::Zscore(merged, pheno = rep(1, N), permute = permute)
    class(z.stat) <- c(class(z.stat), "posterior.mean")
    return(z.stat)
  }




.w.pds.plotExpressionHeatmapMerged =
  function(fitResult, lambda, ematMerged, sampleMetadata, annoNames, annoLevels, annoColors, clusterSamplesTogether=FALSE, geneIdOrder=NA, className='class', classLevels=NA, metaAnalysisName='metaAnalysis', width=8, height=8,color = colorRampPalette(rev(brewer.pal(n=7, name='RdBu')))(100), geneSymbol = F, ...) {
    coefResult = coef(fitResult, s=lambda)
    coefDf = makeCoefDf(coefResult)
    geneIds = coefDf[coefDf[,'geneId']!='(Intercept)', 'geneId']
    if(geneSymbol) { geneSymbols = annotate::getSYMBOL(geneIds, 'org.Hs.eg')
    geneTexts = sprintf('%s (%s)', geneSymbols, geneIds) } else {
    geneTexts = geneIds }
    names(geneTexts) = geneIds
    emat = ematMerged[geneIds,]
    
    # order the samples
    if (clusterSamplesTogether) {
      d = dist(t(emat))
      co = cba::order.optimal(d, hclust(d)$merge)
      emat = emat[,co$order]
    } else {
      if (is.na(classLevels[1])) {
        classLevels = unique(sampleMetadata[colnames(ematMerged), className])}
      ematSmallList = foreach(classLevel=classLevels) %do% {
        x = emat[,colnames(emat) %in% sampleMetadata[sampleMetadata[,className]==classLevel, 'sample']]
        d = dist(t(x))
        co = cba::order.optimal(d, hclust(d)$merge)
        x = x[,co$order]}
      emat = do.call(cbind, ematSmallList)}
    
    # order the genes
    if (is.na(geneIdOrder[1])) {
      d = dist(emat)
      co = cba::order.optimal(d, hclust(d)$merge)
      emat = emat[co$order,]
      rownames(emat) = geneTexts[co$order]
    } else {
      emat = emat[geneIdOrder,]
      rownames(emat) = geneTexts[geneIdOrder]}
    
    # scale the matrix
    emat = t(scale(t(emat)))
    emat[emat>3] = 3
    emat[emat<(-3)] = -3
    
    annotation = sampleMetadata[colnames(ematMerged), annoNames, drop=FALSE]
    for (annoName in annoNames) {
      if (!is.na(annoLevels[[annoName]][1])) {
        annotation[,annoName] = factor(annotation[,annoName], levels=annoLevels[[annoName]])}}
    
    pdf(file=sprintf('%s_lambda_%.3g_heatmap.pdf', metaAnalysisName, lambda), width=width, height=height)
    pheatmap::pheatmap(emat, color,
                       breaks=seq(from=-3, to=3, length.out=101), cluster_rows=FALSE, cluster_cols=FALSE, treeheight_row=0,
                       treeheight_col=0, show_colnames=FALSE, border_color=NA, annotation=annotation, annotation_colors=annoColors, ...)
    dev.off()}




.w.writeConfusionCrossValidation =
  function(cvFit, lambda, ematMerged, sampleMetadata, className='class',
           classLevels=NA, metaAnalysisName='metaAnalysis') {

    if (is.na(classLevels[1])) {
      classLevels = names(cvFit$glmnet.fit$beta)}

    if (class(cvFit$glmnet.fit)[1] == "lognet" ) {
      cvProbs = cvFit$fit.preval[,which.min(abs(cvFit$lambda - lambda)), drop=F]
      cvProbs <- cbind(1-cvProbs, cvProbs) } else {
        cvProbs <- cvFit$fit.preval[, , which.min(abs(cvFit$lambda - lambda))]
      }
    rownames(cvProbs) <- colnames(ematMerged)
    colnames(cvProbs) <- cvFit$glmnet.fit$classnames
    cvProb <<- cvProbs
    preds = colnames(cvProbs)[apply(cvProbs, MARGIN=1, function(x) which.max(x))]
    predsFactor = factor(preds, levels=classLevels)
    trueClasses =  factor(sampleMetadata[colnames(ematMerged), className], levels=classLevels)
    cv.preds <<- preds
    cv.predsFactor <<- predsFactor
    cv.trueClasses <<- trueClasses
    confus.cv <<- table(trueClasses, predsFactor)
    write.csv(confus.cv, file=sprintf('%s_cv_lambda_%.3g_confusion.csv', metaAnalysisName, lambda), quote=FALSE)
  }



.w.plotClassProbsCrossValidation =
  function(cvFit, lambda, sampleMetadata, discoveryStudyNames, discoverySampleNames,
           className, classesTrain, metaAnalysisName, size=2, width=8, height=12, ggplotArgs=NA, ddf=NULL, ...) {
    if (class(cvFit$glmnet.fit)[1] == "lognet" ) {
      cvProbs = cvFit$fit.preval[,which.min(abs(cvFit$lambda - lambda)), drop=F]
      cvProbs <- cbind(1-cvProbs, cvProbs) } else {
        cvProbs = cvFit$fit.preval[, , which.min(abs(cvFit$lambda - lambda))]
      }

    pList = list()
    ddf=NULL; ddf=list()
    for (discoveryStudyName in discoveryStudyNames) {
      discoverySampleNamesNow = discoverySampleNames[sampleMetadata[discoverySampleNames, 'study']==discoveryStudyName]
      df = data.frame(cvProbs[sampleMetadata[discoverySampleNames, 'study']==discoveryStudyName,])
      colnames(df) = cvFit$glmnet.fit$classnames
      df[,'study'] = discoveryStudyName
      df[,'sample'] = discoverySampleNamesNow
      df[,'trueClass'] = factor(sampleMetadata[discoverySampleNamesNow, className], levels=classesTrain)
      df[,'trueClassProb'] = apply(df, MARGIN=1, function(x) as.numeric(x[x['trueClass']]))

      df = df[order(df[,'trueClass'], -df[,'trueClassProb']),]
      df = do.call(rbind, lapply(classesTrain, function(x) df[df[,'trueClass']==x,]))
      ddf[[discoveryStudyName]] <<- df
      idxTmp = c()
      for (classTrain in classesTrain) {
        if (any(df[,'trueClass']==classTrain)) {
          idxTmp = c(idxTmp, 1:(sum(df[,'trueClass']==classTrain)))}}
      df[,'idx'] = idxTmp

      dfMolten = melt(df, measure.vars=classesTrain, variable.name='probClass', value.name='prob')
      p = ggplot(dfMolten) + facet_grid(study ~ trueClass, scales='free_x', space='free_x', ...) +
        geom_point(aes_string(x='idx', y='prob', color='probClass', shape='probClass' ), size=size) + geom_hline(yintercept=.5, linetype="dashed", color = "red", size=.5) +
        labs( x="", y='Probability') + theme(legend.title=element_blank()) + scale_x_continuous( breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))
      if (!is.na(ggplotArgs[1])) {
        for (ggplotArg in ggplotArgs) {
          p = p + ggplotArg }}
      pList[[discoveryStudyName]] = p}

    g = do.call(arrangeGrob, c(pList, list(nrow=length(discoveryStudyNames))))
    ggsave(sprintf('%s_cv_class_probs.pdf', metaAnalysisName), plot=g, width=width, height=height)}




.w.writeConfusionValidation =
  function(predsList, lambda, sampleMetadata, className='class',
           classLevels=NA, metaAnalysisName='metaAnalysis', familyName) {
    if (is.na(classLevels[1])) {
      classLevels = colnames(predsList[[1]])}

    if (familyName == "binomial" ) {
      predsProb = do.call(rbind, lapply(predsList, function(x) x))
      predsProb = cbind(1-predsProb, predsProb)
      colnames(predsProb) = c("C", "E")

    } else {
      predsProb = do.call(rbind, lapply(predsList, function(x) x[,,1]))

    }
    valProb <<- predsProb
    predsClass = colnames(predsProb)[apply(predsProb, MARGIN=1, function(x) which.max(x))]
    predsFactor <<- factor(predsClass, levels=classLevels)
    trueClasses = factor(sampleMetadata[rownames(predsProb), className], levels=classLevels)
    val.trueClasses <<- trueClasses
    confus <<- table(trueClasses, predsFactor)
    write.csv(confus, file=sprintf('%s_val_lambda_%.3g_confusion.csv', metaAnalysisName, lambda), quote=FALSE)}


.w.writeConfusionValidationEach =
  function(predsList, lambda, sampleMetadata, className='class',
           classLevels=NA, metaAnalysisName='metaAnalysis', valProb.ind = NULL, familyName) {

    if (is.na(classLevels[1])) {
      classLevels = colnames(predsList[[1]])}
    each.val=list()
    valProb.ind=list()
    confus.each=list()

    for (validationStudyName in names(predsList)) {
      if (familyName == "binomial" ) {
        predsProb = predsList[[validationStudyName]]
        predsProb = cbind(1-predsProb, predsProb)
        colnames(predsProb) = c("C", "E")
      } else {
        predsProb = predsList[[validationStudyName]][,,1]
      }
      predsClass = colnames(predsProb)[apply(predsProb, MARGIN=1, function(x) which.max(x))]
      predsFactor = factor(predsClass, levels=classLevels)
      trueClasses = factor(sampleMetadata[rownames(predsProb), className], levels=classLevels)
      each.val[[validationStudyName]] <-data.frame(predsFactor=predsFactor, trueClasses=trueClasses, predsProb= predsProb[,2])

      confus.each[[validationStudyName]] <- table(trueClasses, predsFactor)
      valProb.ind[[validationStudyName]] <- data.frame(predsProb=predsProb[,2], trueClasses=trueClasses)
      write.csv(confus.each[[validationStudyName]], file=sprintf('%s_val_%s_lambda_%.3g_confusion.csv', metaAnalysisName, validationStudyName, lambda), quote=FALSE)

    }
    each.val <<- each.val
    valProb.ind <<-valProb.ind
    confus.each <<- confus.each 
    print("ok")
  }



.w.plotClassProbsValidation <-
  function(predsList, sampleMetadata, className, classesTrain,
           size=2, width=8, height=9, ggplotArgs=NA, familyName, ...) {
    pList = list()
    for (validationStudyName in names(predsList)) {
      if (familyName == "binomial" ) {
        d = predsList[[validationStudyName]]
        d = cbind(1-d, d)
        colnames(d) = c("C", "E")
        df = data.frame(d)
      } else {
        df = data.frame(predsList[[validationStudyName]][,,1])
      }
      df[,'study'] = sampleMetadata[rownames(df), 'study']
      df[,'sample'] = rownames(df)
      df[,'trueClass'] = factor(sampleMetadata[rownames(df), className], levels=classesTrain)
      df[,'trueClassProb'] = apply(df, MARGIN=1, function(x) as.numeric(x[x['trueClass']]))

      df = df[order(df[,'trueClass'], -df[,'trueClassProb']),]
      df = do.call(rbind, lapply(classesTrain, function(x) df[df[,'trueClass']==x,]))

      idxTmp = c()
      for (classTrain in classesTrain) {
        if (any(df[,'trueClass']==classTrain)) {
          idxTmp = c(idxTmp, 1:(sum(df[,'trueClass']==classTrain)))}}
      df[,'idx'] = idxTmp
      rownames(df) = NULL

      dfMolten = melt(df, measure.vars=classesTrain, variable.name='probClass', value.name='prob')
      p = ggplot(dfMolten) + facet_grid(study ~ trueClass, scales='free_x', space='free_x') +
        geom_point(aes_string(x='idx', y='prob', color='probClass', shape='probClass'), size=size) +  geom_hline(yintercept=.5, linetype="dashed", color = "red", size=.5) +
        labs(x='Sample', y='Probability') + theme(legend.title=element_blank())
      if (!is.na(ggplotArgs[1])) {
        for (ggplotArg in ggplotArgs) {
          p = p + ggplotArg}}
      pList[[validationStudyName]] = p}

    g = do.call(arrangeGrob, c(pList, list(nrow=length(predsList))))
    ggsave(sprintf('%s_val_class_probs.pdf', metaAnalysisName), plot=g, width=width, height=height)}


# .w.  validationStudyNames <- val.list
.w.predictValidationData <- 
  function(ematList, studyMetadata, sampleMetadata, discoverySampleNames, classesTrain,
           alpha, lambda, weights, batchColname='study', covariateName=NA, className='class',
           familyName='binomial', predType='response', intercept=TRUE, val.list=NA) {
    # .w. val.list
    if(!is.na(val.list)) {discoveryStudyNames = studyMetadata[studyMetadata[,'discovery'], 'study'] 
    validationStudyNames <- val.list } else {   
      discoveryStudyNames = studyMetadata[studyMetadata[,'discovery'], 'study']
      validationStudyNames = studyMetadata[studyMetadata[,'validation'], 'study'] }
    
    predsList = foreach(validationStudyName=validationStudyNames) %do% {
      if(!is.na(val.list)) validationStudyName <- unlist(validationStudyName)
      idxValidation = sampleMetadata[,'study'] %in% validationStudyName & (sampleMetadata[,className] %in% classesTrain)
      validationSampleNames = sampleMetadata[idxValidation, 'sample']
      
      ematListNow = ematList[c(discoveryStudyNames, validationStudyName)]
      ematMergedDiscVal = mergeStudyData(ematListNow, sampleMetadata, batchColname=batchColname, covariateName=covariateName)
      ematMergedDisc = ematMergedDiscVal[,discoverySampleNames]
      
      fitResult = glmnet(t(ematMergedDisc), sampleMetadata[discoverySampleNames, className], alpha=alpha, lambda=lambda,
                         weights=weights[discoverySampleNames], family=familyName, standardize=FALSE, intercept=intercept)
      preds = predict(fitResult, newx=t(ematMergedDiscVal[,validationSampleNames]), s=lambda, type=predType)}
    
    names(predsList) = validationStudyNames
    return(predsList)}




.w.PDS.predictValidationData <- function(ematList, studyMetadata, sampleMetadata, discoverySampleNames, classesTrain,
                                         alpha, lambda, weights, batchColname='study', covariateName=NA, className='class',
                                         familyName='binomial', predType='response', intercept=TRUE, ...) {
  
  discoveryStudyNames = studyMetadata[studyMetadata[,'discovery'], 'study']
  validationStudyNames = studyMetadata[studyMetadata[,'validation'], 'study']
  
  predsPDSList = foreach(validationStudyName=validationStudyNames) %do% {
    idxValidation = sampleMetadata[,'study']==validationStudyName & (sampleMetadata[,className] %in% classesTrain)
    validationSampleNames = sampleMetadata[idxValidation, 'sample']
    
    ematListNow = ematList[c(discoveryStudyNames, validationStudyName)]
    ematMergedDiscVal = .w.mergeStudyData(ematListNow, sampleMetadata, batchColname=batchColname, covariateName=covariateName, ...)
    ematMergedDisc = ematMergedDiscVal[,discoverySampleNames]
    
    exp <- ematMergedDiscVal
    exp.merged = exp [which(rownames(exp) %in% gene_ids),] ; dim(exp)
    nor <- rep(F, dim(ematMergedDiscVal)[2])
    nor[which(colnames(ematMergedDisc) %in% sampleMetadata$sample[which(sampleMetadata$class == "C" )])] <- T
    normal.merged <-nor
    
    PDS.result <-wrap.PDS(exp.merged, rownames(exp.merged), pw.in$entrez_gene_ids, pw.in$pathway, normal.merged, attempts = 2, min_exp = -3, min_std = 0.4)
    
    PDSmatrix <- do.call(rbind.data.frame, PDS.result$scores)
    print(dim(PDSmatrix))
    d <-as.matrix(PDSmatrix)
    colnames(d)<-colnames(ematMergedDiscVal)
    
    fitResult = glmnet(t(d[,discoverySampleNames]), sampleMetadata[discoverySampleNames, className], alpha=alpha, lambda=lambda, weights=weights[discoverySampleNames], family=familyName, standardize=FALSE, intercept=intercept)
    preds = predict(fitResult, newx=t(d[,validationSampleNames]), s=lambda, type=predType)
  }
  names(predsPDSList) = validationStudyNames
  return(predsPDSList)}





.w.mergeStudyData <- function(ematList, sampleMetadata, batchColname='study', covariateName=NA,
                              batchCorrection=TRUE, parPrior=TRUE,  merge.mathod=c("combat", "reference.combat", "tdm", "c.bind"),  ref.batch = 1) {
  merge.mathod <-match.arg(merge.mathod)  # .w.
  
  # merge data and perform cross-study normalization
  geneIds = Reduce(intersect, lapply(ematList, function(x) rownames(x)))
  ematList2 = foreach(studyName=names(ematList)) %do% {ematNow = ematList[[studyName]][geneIds,]}
  if (batchCorrection) {
    # if both one-color and two-color data is present, ComBat can fail catastrophically, if data is not scaled beforehand
    ematListScaled = lapply(ematList2, function(emat) (emat - mean(emat)) / sd(emat))
    ematMerged = do.call(cbind, ematListScaled)
    if (is.na(covariateName)) {
      covariateInfo = model.matrix(~rep_len(1, ncol(ematMerged)))
    } else {
      covariateInfo = model.matrix(~sampleMetadata[colnames(ematMerged), covariateName])}
    
    if (length(unique(sampleMetadata[colnames(ematMerged), batchColname]))>1) {
      if(merge.mathod == "combat") {
        ematMergedNorm = sva::ComBat(ematMerged, batch=sampleMetadata[colnames(ematMerged), batchColname],
                                     mod=covariateInfo, par.prior=parPrior) }
      if(merge.mathod == "reference.combat") {
        ematMergedNorm = .w.ComBat(ematMerged, batch=sampleMetadata[colnames(ematMerged), batchColname],
                                   mod=covariateInfo, par.prior=parPrior, ref.batch =  ref.batch) }
      if(merge.mathod == "tdm") {
        ematMergedNorm = .w.tdm(ematList2[[1]], ematList2[[2]])
      }
      if(merge.mathod == "c.bind") {
        ematMergedNorm = cbind(ematList2[[1]], ematList2[[2]])
      }
      
       
    } else {
      ematMergedNorm = ematMerged}
    
    
    batch <<- c(batch, list(batch=unique(sampleMetadata[colnames(ematMerged), batchColname])))
    batch.expr <<-c(batch.expr, list(ematMergedNorm=ematMergedNorm)) 
    
    
    return(ematMergedNorm)
  } else {
    return(do.call(cbind, ematList2))}}





.w.dual.predictValidationData <- function(ematList, studyMetadata, sampleMetadata, discoverySampleNames, classesTrain,
                                          alpha, lambda, weights, batchColname='study', covariateName=NA, className='class',
                                          familyName='binomial', predType='response', intercept=TRUE, val.list=NA, type.m =c("m", "PDS"), pw.in = NULL, gene_ids=NULL, min_std = 0.4, attempts = 2, min_exp = -3, ...) {
  type.m <- match.arg(type.m)
  
  # .w. val.list
  if(!is.na(val.list)) {discoveryStudyNames = studyMetadata[studyMetadata[,'discovery'], 'study']     
  validationStudyNames <- val.list } else {   
    discoveryStudyNames = studyMetadata[studyMetadata[,'discovery'], 'study']
    validationStudyNames = studyMetadata[studyMetadata[,'validation'], 'study'] }
  
  if(type.m == "m") {
    predsList = foreach(validationStudyName=validationStudyNames) %do% {
      if(!is.na(val.list)) validationStudyName <- unlist(validationStudyName)
      idxValidation = sampleMetadata[,'study'] %in% validationStudyName & (sampleMetadata[,className] %in% classesTrain)
      validationSampleNames = sampleMetadata[idxValidation, 'sample']
      
      ematListNow = ematList[c(discoveryStudyNames, validationStudyName)]
      ematMergedDiscVal = .w.mergeStudyData(ematListNow, sampleMetadata, batchColname=batchColname, covariateName=covariateName,...)
      ematMergedDisc = ematMergedDiscVal[,discoverySampleNames]
      
      fitResult = glmnet(t(ematMergedDisc), sampleMetadata[discoverySampleNames, className], alpha=alpha, lambda=lambda,
                         weights=weights[discoverySampleNames], family=familyName, standardize=FALSE, intercept=intercept)
      preds = predict(fitResult, newx=t(ematMergedDiscVal[,validationSampleNames]), s=lambda, type=predType)}
    
    names(predsList) = validationStudyNames
    return(predsList)
    
  } else {      # .w.
    predsPDSList = foreach(validationStudyName=validationStudyNames) %do% {
      if(!is.na(val.list)) validationStudyName <- unlist(validationStudyName)
      idxValidation = sampleMetadata[,'study'] %in% validationStudyName & (sampleMetadata[,className] %in% classesTrain)
      validationSampleNames = sampleMetadata[idxValidation, 'sample'] 
      
      ematListNow = ematList[c(discoveryStudyNames, validationStudyName)]
      ematMergedDiscVal = .w.mergeStudyData(ematListNow, sampleMetadata, batchColname=batchColname, covariateName=covariateName, ...)
      ematMergedDisc = ematMergedDiscVal[,discoverySampleNames]
      exp <- ematMergedDiscVal
      exp.merged = exp [which(rownames(exp) %in% gene_ids),] ; dim(exp)
      nor <- rep(F, dim(ematMergedDiscVal)[2])
      nor[which(colnames(ematMergedDisc) %in% sampleMetadata$sample[which(sampleMetadata$class == "C" )])] <- T
      normal.merged <-nor
      PDS.result <-wrap.PDS(exp.merged, rownames(exp.merged), pw.in$entrez_gene_ids, pw.in$pathway, normal.merged, attempts = attempts, min_exp = min_exp, min_std = min_std)
      PDSmatrix <- do.call(rbind.data.frame, PDS.result$scores)
      print(dim(PDSmatrix))
      d <-as.matrix(PDSmatrix)
      colnames(d)<-colnames(ematMergedDiscVal)
      fitResult = glmnet(t(d[,discoverySampleNames]), sampleMetadata[discoverySampleNames, className], alpha=alpha, lambda=lambda, weights=weights[discoverySampleNames], family=familyName, standardize=FALSE, intercept=intercept)
      preds = predict(fitResult, newx=t(d[,validationSampleNames]), s=lambda, type=predType)
    }
    names(predsPDSList) = validationStudyNames
    predsList <- predsPDSList
    return(predsList)
  }
}



makeMatchSampleMapping = function(metadata, subStudyNames, matchColname) {
  metadataNow = metadata[metadata[,'study'] %in% subStudyNames, c('study', 'sample', matchColname)]
  metadataNow = metadataNow[order(metadataNow[,'study'], decreasing=is.unsorted(subStudyNames)), c('sample', matchColname)]
  headFunc = function(x) x[1]
  mappingDf = metadataNow %>% group_by_(.dots=list(matchColname)) %>% summarise_each(funs(headFunc)) %>%
    data.frame(check.names=FALSE)
  mapping = mappingDf[,'sample']
  names(mapping) = mappingDf[,matchColname]
  return(mapping)}


mergeMatchStudyData = function(ematAtomicList, studyMetadataAtomic, sampleMetadataAtomic, matchColname,
                               mergeFunc=function(x) mean(x, na.rm=TRUE)) {
  ematList = list()
  sampleMetadataList = list()

  for (matchStudyName in unique(studyMetadataAtomic[,'matchStudy'])) {
    if (sum(studyMetadataAtomic[,'matchStudy']==matchStudyName)==1) {
      ematList[[matchStudyName]] = ematAtomicList[[matchStudyName]]
      sampleMetadataList[[matchStudyName]] = sampleMetadataAtomic[sampleMetadataAtomic[,'study']==matchStudyName,]

    } else if (sum(studyMetadataAtomic[,'matchStudy']==matchStudyName)>1) {
      atomicStudyNames = studyMetadataAtomic[studyMetadataAtomic[,'matchStudy']==matchStudyName, 'study']
      edfListNow = list()
      for (atomicStudyName in atomicStudyNames) {
        edf = data.frame(rownames(ematAtomicList[[atomicStudyName]]), ematAtomicList[[atomicStudyName]])
        rownames(edf) = NULL
        colnames(edf) = c('geneId', sampleMetadataAtomic[colnames(edf)[2:ncol(edf)], matchColname])
        edfListNow[[atomicStudyName]] = edf}

      edfBound = suppressWarnings(rbind_all(edfListNow))
      edfMerged = edfBound %>% group_by_(.dots=list('geneId')) %>% summarise_each(funs(mergeFunc)) %>%
        data.frame(check.names=FALSE)
      rownames(edfMerged) = edfMerged[,'geneId']
      edfMerged = edfMerged[,-1]

      mapping = makeMatchSampleMapping(sampleMetadataAtomic, atomicStudyNames, matchColname)
      colnames(edfMerged) = mapping[colnames(edfMerged)]
      ematList[[matchStudyName]] = as.matrix(edfMerged)

      idx = (sampleMetadataAtomic[,'study'] %in% atomicStudyNames) &
        (sampleMetadataAtomic[,'sample'] %in% colnames(edfMerged))
      sampleMetadataList[[matchStudyName]] = sampleMetadataAtomic[idx,]
      sampleMetadataList[[matchStudyName]][,'study'] = matchStudyName}}

  headFunc = function(x) x[1]
  studyMetadata = studyMetadataAtomic %>% group_by_(.dots=list('matchStudy')) %>% summarise_each(funs(headFunc)) %>%
    data.frame(check.names=FALSE)
  studyMetadata = studyMetadata[,colnames(studyMetadata)!='study']
  colnames(studyMetadata)[colnames(studyMetadata)=='matchStudy'] = 'study'
  rownames(studyMetadata) = studyMetadata[,'study']
  studyMetadata[,'matchStudy'] = studyMetadata[,'study']

  sampleMetadata = suppressWarnings(rbind_all(sampleMetadataList)) %>% data.frame(check.names=FALSE)
  rownames(sampleMetadata) = sampleMetadata[,'sample']

  result = list(ematList, studyMetadata, sampleMetadata)
  names(result) = c('ematList', 'studyMetadata', 'sampleMetadata')
  return(result)}


mergeStudyData = function(ematList, sampleMetadata, batchColname='study', covariateName=NA,
                          batchCorrection=TRUE, parPrior=TRUE) {
  # merge data and perform cross-study normalization
  geneIds = Reduce(intersect, lapply(ematList, function(x) rownames(x)))
  ematList2 = foreach(studyName=names(ematList)) %do% {ematNow = ematList[[studyName]][geneIds,]}
  if (batchCorrection) {
    # if both one-color and two-color data is present, ComBat can fail catastrophically, if data is not scaled beforehand
    ematListScaled = lapply(ematList2, function(emat) (emat - mean(emat)) / sd(emat))
    ematMerged = do.call(cbind, ematListScaled)
    if (is.na(covariateName)) {
      covariateInfo = model.matrix(~rep_len(1, ncol(ematMerged)))
    } else {
      covariateInfo = model.matrix(~sampleMetadata[colnames(ematMerged), covariateName])}

    if (length(unique(sampleMetadata[colnames(ematMerged), batchColname]))>1) {
      ematMergedNorm = ComBat(ematMerged, batch=sampleMetadata[colnames(ematMerged), batchColname],
                              mod=covariateInfo, par.prior=parPrior)
    } else {
      ematMergedNorm = ematMerged}

    return(ematMergedNorm)
  } else {
    return(do.call(cbind, ematList2))}}


makeGlmnetArgs = function(metadata, foldidColname='study') {
  # construct foldid and weights for glmnet
  foldid = as.numeric(factor(metadata[,foldidColname], labels=1:length(unique(metadata[,foldidColname]))))
  names(foldid) = rownames(metadata)
  weights = length(unique(foldid)) /
    do.call("c", sapply(sapply(unique(foldid), function(x) sum(foldid==x)), function(n) rep_len(n, n), simplify=FALSE))
  names(weights) = rownames(metadata)
  return(list(foldid=foldid, weights=weights))}


crossValidateMerged = function(ematMerged, sampleMetadata, weights, alphas, nFolds=10, foldid=NA, nRepeats=3,
                               yName='class', clinVarColnames=NA, ...) {
  args = list(...)
  if (!is.null(args[['family']]) & args[['family']]=='cox') {
    y = as.matrix(sampleMetadata[colnames(ematMerged), yName])
    colnames(y) = c('time', 'status')
  } else {
    y = sampleMetadata[colnames(ematMerged), yName]}

  if (is.na(clinVarColnames[1])) {
    x = scale(t(ematMerged), center=TRUE, scale=FALSE)
  } else {
    clinVarTmp = data.frame(lapply(sampleMetadata[colnames(ematMerged), clinVarColnames], factor))
    clinDummy = model.matrix(~ 0 + ., data=clinVarTmp)
    x = cbind(scale(t(ematMerged), center=TRUE, scale=FALSE), clinDummy)}

  if (is.na(foldid[1])) {
    cvFitList = list()
    for (ii in 1:nRepeats) {
      foldid = sample(rep(seq(nFolds), length=ncol(ematMerged)))
      cvFitList[[ii]] = foreach(alpha=alphas) %do% {
        cv.glmnet(x, y, weights=weights[colnames(ematMerged)], foldid=foldid, alpha=alpha, standardize=FALSE, ...)}}
  } else {
    cvFitList = foreach(alpha=alphas) %do% {
      cv.glmnet(x, y, weights=weights[colnames(ematMerged)], foldid=foldid[colnames(ematMerged)], alpha=alpha,
                standardize=FALSE, ...)}}
  return(cvFitList)}


predictValidationData = function(ematList, studyMetadata, sampleMetadata, discoverySampleNames, classesTrain,
                                 alpha, lambda, weights, batchColname='study', covariateName=NA, className='class',
                                 familyName='binomial', predType='response', intercept=TRUE) {

  discoveryStudyNames = studyMetadata[studyMetadata[,'discovery'], 'study']
  validationStudyNames = studyMetadata[studyMetadata[,'validation'], 'study']

  predsList = foreach(validationStudyName=validationStudyNames) %do% {
    idxValidation = sampleMetadata[,'study']==validationStudyName & (sampleMetadata[,className] %in% classesTrain)
    validationSampleNames = sampleMetadata[idxValidation, 'sample']

    ematListNow = ematList[c(discoveryStudyNames, validationStudyName)]
    ematMergedDiscVal = mergeStudyData(ematListNow, sampleMetadata, batchColname=batchColname, covariateName=covariateName)
    ematMergedDisc = ematMergedDiscVal[,discoverySampleNames]

    fitResult = glmnet(t(ematMergedDisc), sampleMetadata[discoverySampleNames, className], alpha=alpha, lambda=lambda,
                       weights=weights[discoverySampleNames], family=familyName, standardize=FALSE, intercept=intercept)
    preds = predict(fitResult, newx=t(ematMergedDiscVal[,validationSampleNames]), s=lambda, type=predType)}

  names(predsList) = validationStudyNames
  return(predsList)}


calcPredictionAuc = function(predsList, sampleMetadata, className='class') {
  nLambda = ncol(predsList[[1]])
  auc = matrix(nrow=nLambda, ncol=length(predsList))
  colnames(auc) = names(predsList)
  for (validationStudyName in names(predsList)) {
    pred = prediction(predsList[[validationStudyName]],
                      matrix(rep(sampleMetadata[rownames(predsList[[validationStudyName]]), className], nLambda), ncol=nLambda),
                      label.ordering=levels(sampleMetadata[,className]))
    auc[,validationStudyName] = sapply(performance(pred, 'auc')@y.values, function(x) x[[1]])}
  return(auc)}


plotCvError = function(cvFit, metaAnalysisName='metaAnalysis', size=0.4, width=5, height=3, ggplotArgs=NA) {
  df = data.frame(log(cvFit$lambda), cvFit$cvm, cvFit$cvlo, cvFit$cvup)
  colnames(df) = c('logLambda', 'cvm', 'cvlo', 'cvup')
  p = ggplot(data=df) + geom_pointrange(aes_string(x='logLambda', y='cvm', ymax='cvup', ymin='cvlo'), size=size) +
    geom_vline(xintercept=log(cvFit$lambda.min), color='blue', linetype='dashed') + labs(x='log(lambda)', y=cvFit$name)
  if (!is.na(ggplotArgs[1])) {
    for (ggplotArg in ggplotArgs) {
      p = p + ggplotArg}}
  ggsave(sprintf('%s_cv_lambda_error.pdf', metaAnalysisName), plot=p, width=width, height=height)}


makeCoefDf = function(coefResult, decreasing=TRUE, classLevels=NA) {
  if (is.list(coefResult)) {
    coefResultNonzero = foreach(coefSparse=coefResult) %do% {
      x = data.frame(rownames(coefSparse)[(coefSparse@i)+1], coefSparse[(coefSparse@i)+1], stringsAsFactors=FALSE)
      colnames(x) = c('geneId', 'coefficient')
      return(x)}
    names(coefResultNonzero) = names(coefResult)
    if (!is.na(classLevels[1])) {
      coefResult = coefResult[classLevels]
      coefResultNonzero = coefResultNonzero[classLevels]}

    for (ii in 1:length(coefResult)) {
      colnames(coefResultNonzero[[ii]])[2] = names(coefResult)[ii]}
    coefDf = Reduce(function(x, y) merge(x, y, by='geneId', all=TRUE), coefResultNonzero)
    idx = do.call(order, c(coefDf[,2:ncol(coefDf)], list(decreasing=decreasing)))
    coefDf = coefDf[idx,]
    coefDf[is.na(coefDf)] = 0

  } else {
    coefDf = data.frame(names(coefResult[(coefResult@i)+1,]), coefResult[(coefResult@i)+1,], stringsAsFactors=FALSE)
    colnames(coefDf) = c('geneId', 'coefficient')
    coefDf = coefDf[order(coefDf[,'coefficient'], decreasing=decreasing),]}
  rownames(coefDf) = NULL
  return(coefDf)}



plotCoefficients = function(fitResult, lambda, classLevels=NA, decreasing=FALSE, geneIdOrder=NA,
                            metaAnalysisName='metaAnalysis', width=4, height=10, ggplotArgs=NA) {
  coefResult = coef(fitResult, s=lambda)
  coefDf = makeCoefDf(coefResult, decreasing=decreasing, classLevels=classLevels)
  coefDf = coefDf[coefDf[,'geneId']!='(Intercept)',]

  if (!is.na(geneIdOrder[1])) {
    rownames(coefDf) = coefDf[,'geneId']
    coefDf = coefDf[geneIdOrder,]
    rownames(coefDf) = NULL}

  if (ncol(coefDf)==2) {
    geneSymbols = getSYMBOL(coefDf[,'geneId'], 'org.Hs.eg')
    coefDf[,'geneId'] = factor(coefDf[,'geneId'], levels=rev(coefDf[,'geneId']),
                               labels=sprintf('%s (%s)', rev(geneSymbols), rev(coefDf[,'geneId'])))
    p = ggplot(data=coefDf) + geom_bar(aes_string(x='geneId', y='coefficient'), stat='identity')

  } else {
    if (is.na(classLevels[1])) {
      classLevels = colnames(coefDf)[2:ncol(coefDf)]}
    coefDfMolten = melt(coefDf, id.vars='geneId', variable.name='class', value.name='coefficient')
    coefDfMolten[,'class'] = factor(coefDfMolten[,'class'], levels=classLevels)

    geneIds = coefDf[,'geneId']
    geneSymbols = getSYMBOL(geneIds, 'org.Hs.eg')
    coefDfMolten[,'geneId'] = factor(coefDfMolten[,'geneId'], levels=rev(geneIds),
                                     labels=sprintf('%s (%s)', rev(geneSymbols), rev(geneIds)))
    p = ggplot(data=coefDfMolten) + facet_wrap(as.formula('~ class'), ncol=ncol(coefDf)-1) +
      geom_bar(aes_string(x='geneId', y='coefficient', fill='class'), stat='identity') + guides(fill=FALSE)}

  p = p + coord_flip() + labs(x='', y='Coefficient')
  if (!is.na(ggplotArgs[1])) {
    for (ggplotArg in ggplotArgs) {
      p = p + ggplotArg}}
  ggsave(filename=sprintf('%s_lambda_%.3g_gene_coef.pdf', metaAnalysisName, lambda), plot=p, width=width, height=height)}



writeCoefficients = function(fitResult, lambda, metaAnalysisName='metaAnalysis', decreasing=TRUE) {
  coefResult = coef(fitResult, s=lambda)
  coefDf = makeCoefDf(coefResult, decreasing=decreasing)
  coefDf1 = coefDf[c(which(coefDf[,'geneId']=='(Intercept)'), which(coefDf[,'geneId']!='(Intercept)')),]
  write.csv(coefDf1, file=sprintf('%s_lambda_%.3g_coef.csv', metaAnalysisName, lambda), quote=FALSE, row.names=FALSE)}


writeConfusionCrossValidation = function(cvFit, lambda, ematMerged, sampleMetadata, className='class',
                                         classLevels=NA, metaAnalysisName='metaAnalysis') {
  if (is.na(classLevels[1])) {
    classLevels = names(cvFit$glmnet.fit$beta)}

  cvProbs = cvFit$fit.preval[,,which.min(abs(cvFit$lambda - lambda))]
  rownames(cvProbs) = colnames(ematMerged)
  colnames(cvProbs) = names(cvFit$glmnet.fit$beta)
  preds = colnames(cvProbs)[apply(cvProbs, MARGIN=1, function(x) which.max(x))]
  predsFactor = factor(preds, levels=classLevels)
  trueClasses = factor(sampleMetadata[colnames(ematMerged), className], levels=classLevels)
  confus = table(trueClasses, predsFactor)
  write.csv(confus, file=sprintf('%s_cv_lambda_%.3g_confusion.csv', metaAnalysisName, lambda), quote=FALSE)}


writeConfusionValidation = function(predsList, lambda, sampleMetadata, className='class',
                                    classLevels=NA, metaAnalysisName='metaAnalysis') {
  if (is.na(classLevels[1])) {
    classLevels = colnames(predsList[[1]])}

  predsProb = do.call(rbind, lapply(predsList, function(x) x[,,1]))
  predsClass = colnames(predsProb)[apply(predsProb, MARGIN=1, function(x) which.max(x))]
  predsFactor = factor(predsClass, levels=classLevels)
  trueClasses = factor(sampleMetadata[rownames(predsProb), className], levels=classLevels)
  confus = table(trueClasses, predsFactor)
  write.csv(confus, file=sprintf('%s_val_lambda_%.3g_confusion.csv', metaAnalysisName, lambda), quote=FALSE)}


writeConfusionValidationEach = function(predsList, lambda, sampleMetadata, className='class',
                                        classLevels=NA, metaAnalysisName='metaAnalysis') {
  if (is.na(classLevels[1])) {
    classLevels = colnames(predsList[[1]])}

  for (validationStudyName in names(predsList)) {
    predsProb = predsList[[validationStudyName]][,,1]
    predsClass = colnames(predsProb)[apply(predsProb, MARGIN=1, function(x) which.max(x))]
    predsFactor = factor(predsClass, levels=classLevels)
    trueClasses = factor(sampleMetadata[rownames(predsProb), className], levels=classLevels)
    confus = table(trueClasses, predsFactor)
    write.csv(confus, file=sprintf('%s_val_%s_lambda_%.3g_confusion.csv', metaAnalysisName, validationStudyName,
                                   lambda), quote=FALSE)}}





.w.plotCoefficients = function(fitResult, lambda, classLevels=NA, decreasing=FALSE, geneIdOrder=NA,
                            metaAnalysisName='metaAnalysis', width=4, height=10, ggplotArgs=NA, PDS=F) {
  coefResult = coef(fitResult, s=lambda)
  coefDf = makeCoefDf(coefResult, decreasing=decreasing, classLevels=classLevels)
  coefDf = coefDf[coefDf[,'geneId']!='(Intercept)',]
  
  if (!is.na(geneIdOrder[1])) {
    rownames(coefDf) = coefDf[,'geneId']
    coefDf = coefDf[geneIdOrder,]
    rownames(coefDf) = NULL}
  
  if (ncol(coefDf)==2) {
    geneSymbols = getSYMBOL(coefDf[,'geneId'], 'org.Hs.eg')
    if (PDS)  { coefDf[,'geneId'] = factor(coefDf[,'geneId'], levels=rev(coefDf[,'geneId']),
                                         labels=sprintf('%s',  rev(coefDf[,'geneId']))) } else {
    
    coefDf[,'geneId'] = factor(coefDf[,'geneId'], levels=rev(coefDf[,'geneId']),
                               labels=sprintf('%s (%s)', rev(geneSymbols), rev(coefDf[,'geneId']))) }
    
    p = ggplot(data=coefDf) + geom_bar(aes_string(x='geneId', y='coefficient'), stat='identity')
    
  } else {
    if (is.na(classLevels[1])) {
      classLevels = colnames(coefDf)[2:ncol(coefDf)]}
    coefDfMolten = melt(coefDf, id.vars='geneId', variable.name='class', value.name='coefficient')
    coefDfMolten[,'class'] = factor(coefDfMolten[,'class'], levels=classLevels)
    
    geneIds = coefDf[,'geneId']
    geneSymbols = getSYMBOL(geneIds, 'org.Hs.eg')
    if (PDS)  { coefDf[,'geneId'] = factor(coefDf[,'geneId'], levels=rev(coefDf[,'geneId']),
                                           labels=sprintf('%s',  rev(coefDf[,'geneId']))) } else {
    coefDfMolten[,'geneId'] = factor(coefDfMolten[,'geneId'], levels=rev(geneIds),
                                     labels=sprintf('%s (%s)', rev(geneSymbols), rev(geneIds))) }

    p = ggplot(data=coefDfMolten) + facet_wrap(as.formula('~ class'), ncol=ncol(coefDf)-1) +
      geom_bar(aes_string(x='geneId', y='coefficient', fill='class'), stat='identity') + guides(fill=FALSE)}
  
  p = p + coord_flip() + labs(x='', y='Coefficient')
  if (!is.na(ggplotArgs[1])) {
    for (ggplotArg in ggplotArgs) {
      p = p + ggplotArg}}
  ggsave(filename=sprintf('%s_lambda_%.3g_gene_coef.pdf', metaAnalysisName, lambda), plot=p, width=width, height=height)}





