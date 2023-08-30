#' @title Results preprocessing
#' @description
#' Enables to combine simulation results stored in different files.
#' @param fileNamesPattern a [regex](regular expression) matching names of the
#' files that should be combined
#' @param path optionally a path to the directory in which these files are
#' stored (if it is different than the current working directory)
#' @param saveToFile optionally the name of the (.RData) file to save the
#' combined results to
#' @param mergeResultsWithConditions a logical value indicating whether
#' variables describing specific values of the data-generating model parameters
#' in each simulation condition should be appended to objects storing simulation
#' results (if not, only the simulation condition number will be stored in these
#' objects, allowing such a merge to be performed in the future)
#' @param keepFileNames a logical value indicating whether information about
#' the source file should be stored with the combined results
#' @param trimIterationsNumber a logical value indicating whether to trim the
#' number of iterations within each simulation condition to the smallest number
#' of iterations among all conditions in the results
#' @returns a list of combined objects
#' @details
#' Please note that the function assigns new (ordinal) numbers to the simulation
#' conditions and iterations to make them unique among all combined files.
#' To emphasize this fact the column storing ordinal number of condition
#' (condition id) is named *condition* in the combined results (instead of
#' *cond*) and the column storing ordinal number of interation (iteration id)
#' is named *iter* in the combined results (instead of *i*).
#' @importFrom dplyr .data
#' @export
combine_results <- function(fileNamesPattern, path = "./",
                            saveToFile = NA_character_,
                            mergeResultsWithConditions = FALSE,
                            keepFileNames = FALSE,
                            trimIterationsNumber = FALSE) {
  stopifnot(is.character(fileNamesPattern), length(fileNamesPattern) == 1,
            !is.na(fileNamesPattern),
            is.character(path), length(path) == 1, !is.na(path),
            dir.exists(path),
            is.character(saveToFile), length(saveToFile) == 1,
            is.logical(mergeResultsWithConditions),
            length(mergeResultsWithConditions) == 1,
            !is.na(mergeResultsWithConditions),
            is.logical(keepFileNames), length(keepFileNames) == 1,
            !is.na(keepFileNames),
            is.logical(trimIterationsNumber), length(trimIterationsNumber) == 1,
            !is.na(trimIterationsNumber))
  files <- list.files(path = path, pattern = fileNamesPattern, full.names = TRUE)
  files <- files[grep("\\.RData", files)]
  if (length(files) == 0L) {
    stop("There are no files matching the provided expression.")
  }
  message("There are ", length(files), " file(s) matching:\n- ",
          paste(sub(paste0(sub("/$", "", path), "/"), "", files),
                collapse = ",\n- "), ".\n")

  message("Reading files...")
  allObjects <- vector(mode = "character", length = 0L)
  for (f in files) {
    objects <- load(f)
    missingObjects <- setdiff(c("conditions", "modelSummaries", "measurePars",
                                "structPars"), objects)
    if (length(missingObjects) > 0L) {
      warning("In file '", f, "' there are missing some objects storing simulation results: `",
              paste0(missingObjects, collapse = "`, `"), "`",
              call. = FALSE, immediate. = TRUE)
    }
    if ("conditions" %in% missingObjects && mergeResultsWithConditions) {
      stop("With argument `mergeResultsWithConditions=TRUE` there must be `conditions` object in each file.")
    }
    unexpectedObjects <- setdiff(objects, c("conditions", "modelSummaries",
                                            "measurePars", "structPars"))
    if (length(unexpectedObjects) > 0L) {
      warning("In file '", f, "' there are some additional objects that were not expected: `",
              paste0(unexpectedObjects, collapse = "`, `"), "`",
              call. = FALSE, immediate. = TRUE)
    }
    allObjects <- union(allObjects, objects)
    if ("conditions" %in% objects) {
      if (!("cond" %in% names(get("conditions")))) {
        conditions <- cbind(get("conditions"), cond = 1L:nrow(get("conditions")))
      }
    }
    for (o in objects) {
      assign(paste0(o, "Temp"),
             dplyr::bind_rows(get0(paste0(o, "Temp"),
                                   ifnotfound = data.frame()),
                              dplyr::mutate(get(o),
                                            file = f)))
    }
    rm(list = objects)
  }
  for (o in 1:length(allObjects)) {
    assign(allObjects[o],
           get(paste0(allObjects[o], "Temp")))
  }
  rm(list = paste0(allObjects, "Temp"))

  message("Ensuring consistency of iteration numbers and condition numbers across files...")
  if ("conditions" %in% allObjects) {
    conditions <- enumerate_conditions(get("conditions"))
    if (mergeResultsWithConditions) {
      condVarsToMerge <- names(conditions)
    } else {
      condVarsToMerge <- c("cond", "file", "condition")
    }
    for (o in setdiff(allObjects, "conditions")) {
      assign(o, dplyr::left_join(get(o), conditions[, condVarsToMerge],
                                 by = c("cond", "file"),
                                 suffix = c(".res", ".gen")))
      condVarsToMerge2 <- condVarsToMerge
      duplicatedNames <- setdiff(condVarsToMerge, names(get(o)))
      condVarsToMerge2[condVarsToMerge2 %in% duplicatedNames] <-
        paste0(condVarsToMerge2[condVarsToMerge2 %in% duplicatedNames], ".gen")
      assign(o, get(o)[, c("condition",
                           setdiff(condVarsToMerge2, c("cond", "condition")),
                           setdiff(names(get(o)), condVarsToMerge2))])
    }
    conditions <- conditions[, !(names(conditions) %in% c("cond", "file"))]
    conditions <- dplyr::distinct(conditions)
  }
  iterations <- enumerate_iterations(lapply(setdiff(allObjects, "conditions"),
                                            function(x) get(x, pos = parent.frame(n = 2))))
  for (o in setdiff(allObjects, "conditions")) {
    assign(o, dplyr::left_join(get(o), iterations,
                               by = c("condition", "i", "file")))
    assign(o, get(o)[, c("condition", "file", "iter",
                         setdiff(names(get(o)), c("condition", "file",
                                                  "iter", "i")))])
  }
  if (trimIterationsNumber) {
    maxIter <- min(dplyr::summarise(dplyr::group_by(iterations, .data$condition),
                                    iter = max(.data$iter),
                                    .groups = "drop")$iter)
    for (o in setdiff(allObjects, "conditions")) {
      assign(o, get(o)[get(o)$iter <= maxIter, ])
    }
  }
  if (!keepFileNames) {
    for (o in allObjects) {
      assign(o, get(o)[, names(get(o)) != "file"])
    }
  }
  if (!is.na(saveToFile) && saveToFile != "") {
    message("Writing combined results to the .RData file...")
    save(list = allObjects,
         file = paste0(sub("\\.([Rr][Dd]ata|rdt)$", "", saveToFile),
                       ".RData"))
  }
  message("Done.")
  results <- lapply(allObjects,
                    function(x) get(x, pos = parent.frame(n = 2)))
  names(results) <- allObjects
  invisible(results)
}
#' @title Results preprocessing
#' @description
#' Provides unique numbering of conditions among multiple files that store
#' simulation results, even if they were ran using different sets of simulation
#' conditions.
#' @param conditions a data frame with combined simulation conditions
#' descriptions
#' @returns provided data frame with added column *condition*
#' @export
enumerate_conditions <- function(conditions) {
  uniqueConditions <- dplyr::distinct(conditions[, setdiff(names(conditions),
                                                           c("cond", "file"))])
  uniqueConditions$condition <- 1L:nrow(uniqueConditions)
  conditions <- dplyr::left_join(conditions, uniqueConditions,
                                 by = setdiff(names(uniqueConditions),
                                              "condition"))
  return(conditions)
}
#' @title Results preprocessing
#' @description
#' Provides unique numbering of iterations among multiple files that store
#' simulation results.
#' @param resultsList a list of data frames with combined simulation results
#' @returns a data frame with columns *condition*, *i*, *file* and *iter*
#' @importFrom dplyr .data
#' @export
enumerate_iterations <- function(resultsList) {
  iterations <- lapply(resultsList,
                       function(x) return(x[, c("condition", "i", "file")]))
  iterations <- dplyr::distinct(dplyr::bind_rows(iterations))
  iterations <- dplyr::arrange(iterations, .data$condition, .data$file, .data$i)
  iterations <- dplyr::mutate(dplyr::group_by(iterations, .data$condition),
                              iter = 1L:dplyr::n())
  iterations <- dplyr::ungroup(iterations)
  return(iterations)
}
