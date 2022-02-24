## Function for pre-processing data
processSurvivalTable <- function(careerData, type, PhDorEMBL, censorYear = 2021) {

  if (!type %in% c("AcPI", "AcOt", "IndR", "NonRes", "NonSci")) stop("Position type not supported") #types of career role
  if (!PhDorEMBL %in% c("PhD","EMBL")) stop("Start time must be either PhD or EMBL") # whether time after PhD or time after EMBL

  #define variable names based on job types and PhD/Postdoc
  posVar <- ifelse(type == "AcPI", "was_GL", paste0("CV_was", type)) #position
  beginVar <- ifelse(PhDorEMBL == "PhD", "phd_year_if_known", "to_year") #pre or post-doc
  timeVar <- ifelse(type == "AcPI", paste0(PhDorEMBL,"toPI"),
                    paste0("CV_",PhDorEMBL, "to", type))

  #check if all the column are present
  stopifnot(c(posVar, beginVar, timeVar) %in% colnames(careerData))

  #calculate values for Cox models
  survT <- select(careerData, all_of(c("unique_ID",posVar, beginVar, timeVar))) %>%
    mutate(ifPos = ifelse(.data[[posVar]] %in% "y", TRUE, FALSE)) %>%
    mutate(timeToPos = ifelse(ifPos, .data[[timeVar]], censorYear - .data[[beginVar]])) %>%
    select(unique_ID, ifPos, timeToPos) %>%
    filter(!is.na(ifPos), !is.na(timeToPos))

  #sainity check
  ## Number of negative time
  negNum <- nrow(survT[survT$timeToPos <0, ])
  totalNum <- nrow(survT)
  eventNum <- nrow(survT[survT$ifPos,])

  #print summarise information
  print(paste0("Number of subjects with negative time: ", negNum))
  print("Negative values are changed to zero")
  survT <- mutate(survT, timeToPos = ifelse(timeToPos <0, 0, timeToPos))
  print(paste0("Total number of subjects used in the analysis: ", totalNum))
  print(paste0("Number of events (being found in a position): ", eventNum))

  return(survT)
}


#Function for cox regression
com <- function(response, time, endpoint, scale =FALSE) {

  if (scale) {
    #calculate z-score
    response <- (response - mean(response, na.rm = TRUE))/sd(response, na.rm=TRUE)
  }
  surv <- coxph(Surv(time, endpoint) ~ response)


  tibble(p = summary(surv)[[7]][,5],
         HR = summary(surv)[[7]][,2],
         lower = summary(surv)[[8]][,3],
         higher = summary(surv)[[8]][,4])
}

#Function to format floats
formatNum <- function(i, limit = 0.01, digits =1, format="e") {
  r <- sapply(i, function(n) {
    if (n < limit) {
      formatC(n, digits = digits, format = format)
    } else {
      format(n, digits = digits)
    }
  })
  return(r)
}


# Function for Kaplan-Meier plot, without faceting
km_noFacet <- function(testTab, strata, titlePlot = "", titleLegend = "Cohort", pval = NULL,
                       stat = "median", maxTime =NULL, showP = TRUE, showTable = FALSE,
                       ylab = "Fraction", xlab = "Time (years)", colList = NULL,
                       table_ratio = c(0.7,0.3)) {

  # function for km plot
  survS <- tibble(time = testTab$timeToPos,
                  endpoint = testTab$ifPos)

  response <- testTab[[strata]]

  # set maximal time limit
  if (!is.null(maxTime)) {
    survS <- mutate(survS, endpoint = ifelse(time > maxTime, FALSE, endpoint),
                    time = ifelse(time > maxTime, maxTime, time))
  } else maxTime <- max(survS$time)

  # change group to factor
  if (!is.factor(response)) {
    survS$group <- factor(response)
  } else {
    survS$group <- response
  }

  #calculate log-rank p value
  if (nlevels(survS$group) > 2) {
    sdf <- survdiff(Surv(survS$time,survS$endpoint) ~ survS$group)
    p <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
  } else {
    p <- com(survS$group, survS$time, survS$endpoint)$p
  }

  if(p< 1e-16) {
    pAnno <- bquote(italic("P")~"< 1e-16")
  } else {
    pval <- formatNum(p, digits = 1)
    pAnno <- bquote(italic("P")~"="~.(pval))
  }

  if (!showP) pAnno <- ""

  if (is.null(colList)) {
    colorPal <- brewer.pal(length(levels(survS$group)), "Set2")
  } else {
    colorPal <- colList[1:length(levels(survS$group))]
  }

  p <- ggsurvplot(survfit(Surv(time, endpoint) ~ group, data = survS),
                  data = survS, pval = FALSE,  conf.int = FALSE,
                  fun = "event",
                  legend = ifelse(showTable, "none","right"),
                  ylab = ylab, xlab = xlab, title = titlePlot,
                  pval.coord = c(0,0.1), risk.table = showTable, xlim = c(0, maxTime), break.time.x=5,
                  legend.labs = levels(survS$group), palette = colorPal, legend.title = titleLegend,
                  risk.table.title = "Number at risk", fontsize=4,
                  ggtheme = theme_classic() + theme(axis.text = element_text(size=15),
                                                    axis.title = element_text(size=15),
                                                    plot.title = element_text(size=15, face ="bold"),
                                                    legend.text = element_text(size=12),
                                                    legend.title = element_text(size=12)))

  if (!showTable) {
    p <- p$plot + annotate("text",label=pAnno, x = 0.1, y=0.1, hjust =0, size =5) +
      geom_hline(yintercept = 1, lty = "dashed")
    return(p)
  } else {
    #construct a gtable
    pp <- p$plot + annotate("text",label=pAnno, x = 0.1, y=0.8, hjust =0, size=5) +
      geom_hline(yintercept = 1, lty = "dashed") +
      theme(axis.title.y = element_text(vjust =-10))
    pt <- p$table + ylab("") + xlab("") + theme(plot.title = element_text(hjust=0, size =10))
    p <- plot_grid(pp,pt, rel_heights = table_ratio, nrow =2, align = "v")
    return(p)
  }
}

# Function for Kaplan-Meier plot
km <- function(testTab, strata, titlePlot = "", titleLegend = "Cohort", pval = NULL, maxTime =NULL,
               ylab = "Fraction", xlab = "Time (years)", colList = NULL, facetBy = "type_pre_postdoc") {

  # function for km plot
  survS <- tibble(time = testTab$timeToPos,
                  endpoint = testTab$ifPos,
                  facetCol = testTab[[facetBy]])

  response <- testTab[[strata]]

  # set maximal time limit
  if (!is.null(maxTime)) {
    survS <- mutate(survS, endpoint = ifelse(time > maxTime, FALSE, endpoint),
                    time = ifelse(time > maxTime, maxTime, time))
  } else maxTime <- max(survS$time)

  # change group to factor
  if (!is.factor(response)) {
    survS$group <- factor(response)
  } else {
    survS$group <- response
  }


  if (is.null(colList)) {
    colorPal <- brewer.pal(length(levels(survS$group)), "Set2")
  } else {
    colorPal <- colList[1:length(levels(survS$group))]
  }

  p <- ggsurvplot(survfit(Surv(time, endpoint) ~ group, data = survS),
                  data = survS, pval = FALSE,  conf.int = FALSE,
                  fun = "event", facet.by = "facetCol",
                  ylab = ylab, xlab = xlab, title = titlePlot,
                  xlim = c(0, maxTime), break.time.x=5,
                  legend.labs = levels(survS$group), short.panel.labs = TRUE,
                  palette = colorPal, legend.title = titleLegend, fontsize=4) +
    geom_hline(yintercept = 1, lty = "dashed") +
    theme_classic() +
    theme(axis.text = element_text(size=15), axis.title = element_text(size=15),
          plot.title = element_text(size=15, face ="bold"), legend.text = element_text(size=12),
          legend.title = element_text(size=12), legend.position = c(0.85,0.75),
          strip.background  = element_rect(fill = "grey80", color=NA), strip.text = element_text(size=15, face = "bold"))

  return(p)
}


#function for pair-wise cox regression test
pairwiseTest <- function(testTab, strata) {
  combo <- combn(levels(testTab[[strata]]),2)
  resTab <- lapply(seq(ncol(combo)), function(i) {
    testVar <- combo[2,i]
    refVar <- combo[1,i]
    eachTab <- filter(testTab, .data[[strata]] %in% c(testVar, refVar)) %>%
      mutate(group = factor(.data[[strata]], levels = c(refVar, testVar)))
    res <- com(eachTab$group, eachTab$timeToPos, eachTab$ifPos)

    tibble(`Test` = testVar,
           `Reference` = refVar,
           `Hazard ratio` = formatC(res$HR, digits=2, format = "f"),
           `95% CI` = paste0(formatC(res$lower, digits=2, format = "f"), " - ", formatC(res$higher, digits=2, format = "f")),
           `P value` = formatC(res$p, digits=2))
  }) %>% bind_rows()
  return(resTab)
}

#Function for KM plot using cohort as strata
plotCohort <- function(careerData, type, PhDorEMBL = "EMBL", titlePlot = "", showTable = FALSE) {

  #get survival table
  survT <- processSurvivalTable(careerData, type,PhDorEMBL)

  testTab <- select(survT,unique_ID, ifPos, timeToPos) %>%
    left_join(select(careerData, unique_ID, type_pre_postdoc, cohort), by = "unique_ID") %>%
    filter(!is.na(cohort)) %>%
    mutate(cohort = factor(cohort),
           type_pre_postdoc = ifelse(type_pre_postdoc == "predoc","PhD alumni",
                                     ifelse(type_pre_postdoc == "postdoc","Postdoc alumni",NA))) %>%
    mutate(type_pre_postdoc = factor(type_pre_postdoc, levels = c("PhD alumni","Postdoc alumni")))

  print("Total subject number stratified by pre/post-doc")
  print(table(testTab$type_pre_postdoc))
  print("Event number stratified by pre/post-doc")
  print(table(filter(testTab, ifPos)$type_pre_postdoc))


  plotOut <- km(testTab, "cohort",  maxTime = 25, titlePlot = type,
                xlab = sprintf("Time after %s (years)", PhDorEMBL),
                ylab = paste0("Probability of being found as ", type)) +
    theme(plot.title = element_text(size=18, hjust=0.5))



  tabOut <-  lapply(levels(testTab$type_pre_postdoc), function(n) {
    eachTab <- filter(testTab, type_pre_postdoc == n)
    pairwiseTest(eachTab, "cohort") %>%
      mutate(Alumni = str_remove(n," alumni"))
  }) %>% bind_rows()

  return(list(plot = plotOut, table = tabOut))

}


# Function for calculating Harrel's C
calcHarrelsC <- function(survT, careerDataSub, preOrPost) {

  combine <- ifelse(preOrPost == "both", TRUE, FALSE)

  df_survan <- select(survT,unique_ID, timeToPos, ifPos) %>%
    left_join(careerDataSub, by = "unique_ID") %>%
    filter(timeToPos >0)

  SurvObject <- Surv(df_survan$timeToPos, df_survan$ifPos)
  options(na.action='na.pass')

  FeatureList <- list(
    publications = dplyr::select(df_survan, starts_with("pubs", ignore.case = TRUE)),
    group_publications = dplyr::select(df_survan, starts_with("group_pubs", ignore.case = TRUE)),
    nationality = model.matrix(~ 0 + nationality, df_survan),
    gender =  model.matrix(~ 0 + gender, df_survan)[,1, drop=FALSE],
    group_leaderseniority = model.matrix(~ 0 + groupleader_seniority, df_survan)[,1, drop=FALSE],
    time_cohort =  model.matrix(~ 0 + cohort + phd_year_if_known + from_year + to_year, df_survan))

  if (combine) {
    FeatureList[["type_pre_postdoc"]] <- model.matrix(~ 0 + type_pre_postdoc, df_survan)[,1, drop=FALSE]
  }

  FeatureList$all <- do.call(cbind, FeatureList)
  sapply(FeatureList, ncol)

  # remove nas
  isna <- sapply(FeatureList, function(f) apply(f,1, function(r) any(is.na(r))))
  isna <- apply(isna,1, any)
  #sum(!isna)
  FeatureList <- lapply(FeatureList, function(mat)  mat[!isna,,drop=F])
  SurvObject <- SurvObject[!isna]
  #unique(sapply(FeatureList, nrow))
  #nrow(SurvObject)

  #stratified CV to include same proportion of uncensored cases
  set.seed(12345)
  uncensored <- SurvObject[,"status"] == 1
  cv_ix <- dismo::kfold(1:nrow(SurvObject), k=10, by=uncensored)

  # fit a cox model and predict on left-out fold
  fit <- lapply(unique(cv_ix), function(i){

    # fit cox model with ridge penalty
    c_all = lapply(FeatureList, function(x) {
      if(ncol(x)<2) {
        coxph(SurvObject[cv_ix!=i,] ~ as.matrix(x[cv_ix!=i,]))
      }
      else {
        #coxph(SurvObject[cv_ix!=i,] ~ as.matrix(x[cv_ix!=i,]))
        cv.glmnet(as.matrix(x[cv_ix!=i,]),SurvObject[cv_ix!=i,], family="cox", alpha=0)
      }
    })

    coefs <- lapply(c_all, function(l) if(class(l) == "cv.glmnet") coef(l, l$lambda.min) else coef(l))
    p_all = mapply(function(x,y) as.matrix(x[cv_ix==i,]) %*% y, FeatureList, coefs)
    p_all <- do.call(cbind, p_all)
    colnames(p_all) <- names(FeatureList)

    # calculate CI on left-out fold
    CI=apply(-p_all,2, Hmisc::rcorr.cens, SurvObject[cv_ix==i])[1,]
    list(CI=CI,c=c, c_all=c_all)
  })

  # get cross-validated CI
  concordanceCV <- sapply(fit, function(l) l$CI)
  colnames(concordanceCV) <- unique(cv_ix)
  df_hc <- melt(concordanceCV, varnames = c("predictor", "cv_idx"))
  df_hc <- df_hc %>% dplyr::rename(CI = value)

  # calculate summary statistics across folds
  summaryHc <- aggregate(df_hc$CI,
                         by = list(predictors = df_hc$predictor),
                         FUN = function(x) c(mean = mean(x), sd = sd(x),
                                             n = length(x)))
  summaryHc <- do.call(data.frame, summaryHc)
  summaryHc$se <- summaryHc$x.sd / sqrt(summaryHc$x.n)

  return(summaryHc)

}


# Function for calculating Harrel's C
calcHarrelsC.bootstrap <- function(survT, careerDataSub, preOrPost, nSample=100, testRatio = 0.3) {

  combine <- ifelse(preOrPost == "both", TRUE, FALSE)

  df_survan <- select(survT,unique_ID, timeToPos, ifPos) %>%
    left_join(careerDataSub, by = "unique_ID") %>%
    filter(timeToPos >0)

  SurvObject <- Surv(df_survan$timeToPos, df_survan$ifPos)
  options(na.action='na.pass')

  FeatureList <- list(
    publications = dplyr::select(df_survan, starts_with("pubs", ignore.case = TRUE)),
    group_publications = dplyr::select(df_survan, starts_with("group_pubs", ignore.case = TRUE)),
    nationality = model.matrix(~ 0 + nationality, df_survan),
    gender =  model.matrix(~ 0 + gender, df_survan)[,1, drop=FALSE],
    group_leaderseniority = model.matrix(~ 0 + groupleader_seniority, df_survan)[,1, drop=FALSE],
    time_cohort =  model.matrix(~ 0 + cohort + phd_year_if_known + from_year + to_year, df_survan))

  if (combine) {
    FeatureList[["type_pre_postdoc"]] <- model.matrix(~ 0 + type_pre_postdoc, df_survan)[,1, drop=FALSE]
  }

  FeatureList$all <- do.call(cbind, FeatureList)
  sapply(FeatureList, ncol)

  # remove nas
  isna <- sapply(FeatureList, function(f) apply(f,1, function(r) any(is.na(r))))
  isna <- apply(isna,1, any)
  #sum(!isna)
  FeatureList <- lapply(FeatureList, function(mat)  mat[!isna,,drop=F])
  SurvObject <- SurvObject[!isna]
  #unique(sapply(FeatureList, nrow))
  #nrow(SurvObject)

  #stratified CV to include same proportion of uncensored cases
  set.seed(12345)
  sampleCensored <- seq(nrow(SurvObject))[SurvObject[,"status"] != 1]
  sampleUncensored <- seq(nrow(SurvObject))[SurvObject[,"status"] == 1]
  bootsSample <- lapply(seq(nSample), function(i) {
    c(sample(sampleCensored, length(sampleCensored)*testRatio),
      sample(sampleUncensored, length(sampleUncensored)*testRatio))
  }) #ensure the balance of censored and uncensored samples

  # fit a cox model and predict on left-out fold
  fit <- lapply(bootsSample, function(i){

    # fit cox model with ridge penalty
    c_all = lapply(FeatureList, function(x) {
      if(ncol(x)<2) {
        coxph(SurvObject[-i,] ~ as.matrix(x[-i,]))
      }
      else {
        #coxph(SurvObject[cv_ix!=i,] ~ as.matrix(x[cv_ix!=i,]))
        cv.glmnet(as.matrix(x[-i,]),SurvObject[-i,], family="cox", alpha=0, nfolds =3)
      }
    })

    coefs <- lapply(c_all, function(l) if(class(l) == "cv.glmnet") coef(l, l$lambda.min) else coef(l))
    p_all = mapply(function(x,y) as.matrix(x[i,]) %*% y, FeatureList, coefs)
    p_all <- do.call(cbind, p_all)
    colnames(p_all) <- names(FeatureList)

    # calculate CI on left-out fold
    CI=apply(-p_all,2, Hmisc::rcorr.cens, SurvObject[i])[1,]
    list(CI=CI,c=c, c_all=c_all)
  })

  # get cross-validated CI
  concordanceCV <- sapply(fit, function(l) l$CI)
  colnames(concordanceCV) <- seq(nSample)
  df_hc <- melt(concordanceCV, varnames = c("predictor", "cv_idx"))
  df_hc <- dplyr::rename(df_hc,CI = value)

  # calculate summary statistics across folds
  summaryHc <- aggregate(df_hc$CI,
                         by = list(predictors = df_hc$predictor),
                         FUN = function(x) c(mean = mean(x), sd = sd(x),
                                             n = length(x)))
  summaryHc <- do.call(data.frame, summaryHc)
  summaryHc$se <- summaryHc$x.sd / sqrt(summaryHc$x.n)

  return(summaryHc)
}

# Function for calculating Harrel's C
calcHarrelsC.noCV <- function(survT, careerDataSub, preOrPost) {

  combine <- ifelse(preOrPost == "both", TRUE, FALSE)

  df_survan <- select(survT,unique_ID, timeToPos, ifPos) %>%
    left_join(careerDataSub, by = "unique_ID") %>%
    filter(timeToPos >0)

  SurvObject <- Surv(df_survan$timeToPos, df_survan$ifPos)
  options(na.action='na.pass')

  FeatureList <- list(
    publications = dplyr::select(df_survan, starts_with("pubs", ignore.case = TRUE)),
    group_publications = dplyr::select(df_survan, starts_with("group_pubs", ignore.case = TRUE)),
    nationality = model.matrix(~ 0 + nationality, df_survan),
    gender =  model.matrix(~ 0 + gender, df_survan)[,1, drop=FALSE],
    group_leaderseniority = model.matrix(~ 0 + groupleader_seniority, df_survan)[,1, drop=FALSE],
    time_cohort =  model.matrix(~ 0 + cohort + phd_year_if_known + from_year + to_year, df_survan))

  if (combine) {
    FeatureList[["type_pre_postdoc"]] <- model.matrix(~ 0 + type_pre_postdoc, df_survan)[,1, drop=FALSE]
  }

  FeatureList$all <- do.call(cbind, FeatureList)
  sapply(FeatureList, ncol)

  # remove nas
  isna <- sapply(FeatureList, function(f) apply(f,1, function(r) any(is.na(r))))
  isna <- apply(isna,1, any)
  #sum(!isna)
  FeatureList <- lapply(FeatureList, function(mat)  mat[!isna,,drop=F])
  SurvObject <- SurvObject[!isna]
  #unique(sapply(FeatureList, nrow))
  #nrow(SurvObject)


  # fit cox model with ridge penalty
  c_all = lapply(FeatureList, function(x) {
    if(ncol(x)<2) {
      coxph(SurvObject ~ as.matrix(x))
    }
    else {
      #coxph(SurvObject[cv_ix!=i,] ~ as.matrix(x[cv_ix!=i,]))
      cv.glmnet(as.matrix(x),SurvObject, family="cox", alpha=0, nfolds =10)
    }
  })

  coefs <- lapply(c_all, function(l) if(class(l) == "cv.glmnet") coef(l, l$lambda.min) else coef(l))
  p_all = mapply(function(x,y) as.matrix(x) %*% y, FeatureList, coefs)
  p_all <- do.call(cbind, p_all)
  colnames(p_all) <- names(FeatureList)

  # calculate CI on left-out fold
  CI=apply(-p_all,2, Hmisc::rcorr.cens, SurvObject)

  summaryHc <- t(CI) %>%
    as_tibble(rownames = "predictors") %>%
    dplyr::rename(x.mean = "C Index", x.sd = S.D., x.n = n) %>%
    mutate(se = x.sd / sqrt(x.n)) %>%
    select(predictors, x.mean, x.sd, x.n, se)

  return(summaryHc)
}


#function for plotting HarrelsC index
plotHarrelsC <- function(careerData, type, preOrPost = "both", method = "CV") {

  PhDorEMBL <- ifelse(preOrPost == "both", "PhD","EMBL")

  #get survival table
  if (preOrPost != "both") {
    careerDataSub <- filter(careerData, type_pre_postdoc == preOrPost)
  } else careerDataSub <- careerData

  survT <- processSurvivalTable(careerDataSub, type, PhDorEMBL)

  #calculate
  if (method == "bootstrapping") {
    summaryHc <- calcHarrelsC.bootstrap(survT, careerDataSub, preOrPost)
  } else if (method == "CV") {
    summaryHc <- calcHarrelsC(survT, careerDataSub, preOrPost)
  } else if (method == "noCV") {
    summaryHc <- calcHarrelsC.noCV(survT, careerDataSub, preOrPost)
  }

  limits <- aes(ymax = summaryHc$x.mean + 1.96*summaryHc$se,
                ymin = summaryHc$x.mean - 1.96*summaryHc$se)
  # barplot
  nameMap <- c(time_cohort = "time & cohort",
               type_pre_postdoc = "predoc/postdoc",
               group_leaderseniority = "group leader seniority",
               nationality = "nationality",
               gender = "gender",
               publications = "publications",
               group_publications = "group publications",
               all = "all")


  plotTab <- summaryHc %>% mutate(predictors = nameMap[as.character(predictors)]) %>%
    mutate(predictors =factor(predictors, nameMap))

  plotTitle <- ifelse(preOrPost == "both", "All subjects",
                      ifelse(preOrPost == "predoc","PhD alumni", "Postdoc alumni"))

  p <- ggplot(plotTab, aes(x=predictors, y=x.mean, fill=x.mean, group=predictors))+
    geom_hline(yintercept = 0.5, color = "grey50") +
    geom_bar(stat = "identity") +
    coord_cartesian(ylim = c(0.30,0.90)) +theme_bw()+
    geom_errorbar(limits, position = "dodge", width = 0.25)+
    scale_fill_gradient(low = "grey80",high="grey10") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=18),
          axis.text.y = element_text(size=18),
          axis.title.y = element_text(size=18),
          strip.background = element_rect(fill="white"),
          strip.text = element_blank(),
          legend.position= "none",
          plot.title = element_text(size=20)) +
    ylab("Harrell's C-index") + xlab("") +
    ggtitle(sprintf("%s (%s)", type, plotTitle))

  return(list(plot = p, table = summaryHc %>% mutate(role = type, preOrPostDoc = preOrPost)))
}

#Function for KM plot using publication as strata
plotPublication <- function(careerData, type, PhDorEMBL = "EMBL", showTable = TRUE) {

  #get survival table
  survT <- processSurvivalTable(careerData, type,PhDorEMBL)

  testTab <- survT %>%
    left_join(select(careerData, unique_ID, type_pre_postdoc, cohort, pubs_FIRST_ra_only_TOTAL), by = "unique_ID") %>%
    filter(!is.na(timeToPos), !is.na(pubs_FIRST_ra_only_TOTAL), !is.na(type_pre_postdoc)) %>%
    mutate(cohort = factor(cohort),
           type_pre_postdoc = ifelse(type_pre_postdoc == "predoc","PhD alumni", "Postdoc alumni")) %>%
    mutate(numPub = ifelse(pubs_FIRST_ra_only_TOTAL > 1, "2+", pubs_FIRST_ra_only_TOTAL)) %>%
    mutate(type_pre_postdoc = factor(type_pre_postdoc, levels = c("PhD alumni","Postdoc alumni")),
           numPub = factor(numPub, levels = c("0", "1", "2+")))

  print("Total subject number stratified by pre/post-doc")
  print(table(testTab$type_pre_postdoc))
  print("Event number stratified by pre/post-doc")
  print(table(filter(testTab, ifPos)$type_pre_postdoc))


  plotOut <- km(testTab, "numPub",  maxTime = 25, titlePlot = type,
                titleLegend = "Number of\nfirst-author\npublications",
                xlab = sprintf("Time after %s (years)", PhDorEMBL),
                ylab = paste0("Probability of being found as ", type)) +
    theme(legend.position = c(0.86,0.65),
          plot.title = element_text(size=18, hjust=0.5))


  tabOut <-  lapply(levels(testTab$type_pre_postdoc), function(n) {
    eachTab <- filter(testTab, type_pre_postdoc == n)
    pairwiseTest(eachTab, "numPub") %>%
      mutate(Alumni = str_remove(n," alumni"))
  }) %>% bind_rows()

  return(list(plot = plotOut, table = tabOut))

}
