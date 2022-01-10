#' productPipeline function for when a new bsp object has been supplied
#'
#' Example: post burn-in
#'
#' function to advance a simulated breeding product pipeline forward by one generation. See Gaynor et al. 2017 for the general idea.
#'
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @param newbsp the new list of breeding scheme parameters
#' @param bsp the original parameter set used
#' @param SP the AlphaSimR SimParam object
#' @return A records object that has new records created by advancing by a generation
#'
#' @details The breeding program product pipeline will have been set by initializeFunc. This function moves the breeding program along by one generation and saves all the resulting phenotypes to the records object.
#'
#' @examples
#' bsp <- specifyPipeline()
#' bsp <- specifyPopulation(bsp)
#' initList <- initializeFunc(bsp)
#' SP <- initList$SP
#' bsp <- initList$bsp
#' records <- initList$records
#' records <- productPipeline(records, bsp, SP)
#' records <- popImprov1(records, bsp, SP)
#'
#' @export
productPipelinePostBurnIn <- function(records, bsp, SP){
  # Calculate the selection criterion. selCritPipeAdv has to be given in bsp
  candidates <- records$F1@id
  selCrit <- bsp$selCritPipeAdv(records, candidates, bsp, SP)

  # Make summary for the incoming F1s
  year <- max(records$stageOutputs$year)+1 # Add a year relative to last year
  nF1 <- bsp$nCrosses * bsp$nProgeny
  nGenoRec <- nInd(records$F1)
  # Analyze the most-recent F1s
  newF1Idx <- nGenoRec - nF1 + 1:nF1
  id <- records$F1[newF1Idx]@id
  records$stageOutputs <- records$stageOutputs %>%
      dplyr::bind_rows(stageOutputs(id=id, f1=records$F1, selCrit=selCrit, stage=0, year=year, bsp=bsp))
  # Will be added to the phenotype records
  toAdd <- list()
  for(stage in bsp$stageNames){
    # Make a summary for this stage
    id <- last(records[[stage]]) %>% .[!. %in% bsp$checks@id]
    records$stageOutputs<-records$stageOutputs %>%
      dplyr::bind_rows(stageOutputs(id=id, f1=records$F1, selCrit=selCrit,
                                            stage=which(bsp$stageNames==stage),
                                            year=year, bsp=bsp))
    if(which(bsp$stageNames==stage) == 1){
      # Stage 1 different: no phenotypes but full Pop-class
      # Use phenotypes to select the F1 going into Stage 1?
      if(bsp$phenoF1toStage1){ # Use phenotypes to choose what goes to Stage 1
        phenoF1 <- setPheno(records$F1[newF1Idx], varE=bsp$errVarPreStage1, onlyPheno=T, simParam=SP)
        indToAdv <- records$F1@id[nGenoRec - nF1 + (phenoF1 %>% order(decreasing=T))[1:bsp$nEntries[stage]] %>% sort]
      } else {
        # Do the F1 have genotypic values that could be used?
        if(selCrit[newF1Idx] %>% is.na %>% all){ # Choose at random
          indToAdv <- records$F1@id[nGenoRec - nF1 + sort(sample(nF1, bsp$nEntries[stage]))]
        } else { # Use selCrit
          indToAdv <- records$F1@id[nGenoRec - nF1 + (selCrit[newF1Idx] %>%
                                                        order(decreasing=T))[1:bsp$nEntries[stage]] %>% sort]
        }
      }
    } else { # Beyond stage 1
      # Don't allow checks to be advanced: use 1:bsp$nEntries[stage-1]
      id <- last(records[[bsp$stageNames[which(bsp$stageNames==stage) -1]]])$id %>% .[!. %in% bsp$checks@id]
      selCritToAdv <- selCrit[id]
      indToAdv<-selCritToAdv %>% sort(.,decreasing = T) %>% .[1:bsp$nEntries[stage]] %>% names
    }

    entries <- records$F1[indToAdv]
    varE <- bsp$gxyVar + (bsp$gxlVar + bsp$gxyxlVar + bsp$errVars[stage] / bsp$nReps[stage]) / bsp$nLocs[stage]
    # reps=1 because varE is computed above
    entries <- setPheno(entries, varE=varE, reps=1, simParam=SP)
    phenoRec <- phenoRecFromPop(entries, bsp, stage)
    # If provided, add checks to the population
    if(!is.null(bsp$checks) & bsp$nChks[stage] > 0){
      varE <- bsp$gxyVar + (bsp$gxlVar + bsp$gxyxlVar + bsp$errVars[stage] / bsp$chkReps[stage]) / bsp$nLocs[stage]
      chkPheno <- setPheno(bsp$checks[1:bsp$nChks[stage]], varE=varE, reps=1, simParam=SP)
      chkRec <- phenoRecFromPop(chkPheno, bsp, stage, checks=T)
      phenoRec <- dplyr::bind_rows(phenoRec, chkRec)
    }
    toAdd <- c(toAdd, list(phenoRec))
  }#END 1:nStages
  for(stage in bsp$stageNames){
    records[[stage]] <- c(records[[stage]], toAdd[which(bsp$stageNames==stage)])
  }

  # Remove old records if needed
  if (length(records[[2]]) > bsp$nCyclesToKeepRecords) records <- removeOldestCyc(records, bsp)

  return(records)
}


#' productSelCritBLUP function
#'
#' function to select parents among individuals with phenotypes, assuming individual effects are IID
#'
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @param candidates Character vector of ids of the candidates to be parents
#' @param bsp The breeding scheme parameter list
#' @param SP The AlphaSimR SimParam object (not used, here for uniformity)
#' @return An IID BLUP of the trait of the candidates
#' @details Accesses all individuals in \code{records} to pick the highest among candidates. If candidates do not have records, a random sample is returned
#'
#' @export
productSelCritBLUP <- function(records, candidates, bsp, SP){
  phenoDF <- dataframePhenoRec(records, bsp)
  # Candidates don't have phenotypes so return random vector
  if (!any(candidates %in% phenoDF$id)){
    crit <- runif(length(candidates))
  } else {
    crit <- iidPhenoEval(phenoDF)
    crit <- crit[candidates]
  }
  names(crit) <- candidates
  return(crit)
}

#' dataframePhenoRec replaces framePhenoRec function post-burn in
#'
#' function to make a data.frame to be used as a source of data to analyze the phenotypic \code{records}
#'
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @param bsp The breeding scheme parameter list
#' @return A data.frame of phenotypic records with four columns: 1. The id of individuals; 2. The trial type of the phenotype record; 3. The year the observation was recorded; 4. The phenotypic value
#' @details \code{records} is a list of lists of populations and is primarily useful for maintaining the phenotypic observations across years and stages. For analysis, you need just the phenotypes in a matrix with relevant independent values
#'
#' @examples
#' phenoDF <- framePhenoRec(records, bsp)
#'
#' @export
dataframePhenoRec <- function(records, bsp){
  stages2frame<-names(records)[!names(records) %in% c("F1","stageOutputs")]
  if(!is.null(bsp$burnInBSP)){ stages_during_burn_in<-bsp$burnInBSP$stageNames } else {
    stages_during_burn_in<-bsp$stageNames }
  current_stages<-bsp$stageNames
  new_stages<-current_stages[!current_stages %in% stages_during_burn_in]

  allPheno<-tibble(Stage=names(records[stages2frame]),
                   Records=records[stages2frame]) %>%
    dplyr::mutate(Records=purrr::map2(Stage,Records,function(Stage,Records){
      if(Stage %in% stages_during_burn_in){
        phenodf<-tibble(year=1:length(Records),df=Records) %>% tidyr::unnest(df) }
      if(length(new_stages)>0){
        if(Stage %in% new_stages){
          phenodf<-tibble(year=(bsp$maxYearBurnInStage+1):(bsp$maxYearBurnInStage+length(Records)),
                          df=Records) %>% tidyr::unnest(df) } }
      return(phenodf)})) %>%
    dplyr::select(-Stage) %>%
    tidyr::unnest(Records)
  return(allPheno)
}
