# Find example in depression study that had most assymetric confidence intervals, while still having estimates for all 4 DBs
# Exclude psychotherapy and ECT
exclude <- c("Psychotherapy", "Electroconvulsive therapy")

results <- read.csv("c:/data/Depression/DepressionResults.csv")
results <- results[results$analysisId == 3, ] # On-treatment, matching
results <- results[!(results$targetName %in% exclude) & !(results$comparatorName %in% exclude), ]
results <- results[results$outcomeType == "hoi", ]
results$ratio <- 1 - (log(results$rr) - log(results$ci95lb)) / (log(results$ci95ub) - log(results$rr))
results$absRatio <- abs(results$ratio)
meanRatio <- aggregate(absRatio ~ targetId + comparatorId + outcomeId, results, mean)
results$count <- 1
dbCount <- aggregate(count ~ targetId + comparatorId + outcomeId, results[!is.na(results$seLogRr), ], sum)
meanRatio <- merge(meanRatio, dbCount)
meanRatio <- meanRatio[meanRatio$count == 4, ]
meanRatio <- meanRatio[order(-meanRatio$absRatio), ]

index <- 1
example <- results[results$targetId == meanRatio$targetId[index] & results$comparatorId == meanRatio$comparatorId[index] & results$outcomeId == meanRatio$outcomeId[index], ]
example

# temp <- example$targetId
# example$targetId <- example$comparatorId
# example$comparatorId <- temp


# Fetch patient-level data from S3 storage ----------------------------------
ArchiveR::listFiles("mschuemi/Epi421")
localFolder <- "C:/home/Research/SmallCountMetaAnalysis/Example"

dbs <- data.frame(s3Folders = c("mschuemi/Epi421/PopEstDepression_Ccae/cmOutput/", 
                                "mschuemi/Epi421/PopEstDepression_Mdcd/cmOutput/",
                                "mschuemi/Epi421/PopEstDepression_Mdcr/cmOutput/",
                                "mschuemi/Epi421/PopEstDepression_Optum/cmOutput/"),
                  database = c("CCAE",
                               "MDCD",
                               "MDCR",
                               "Optum"),
                  stringsAsFactors = FALSE)

for (database in dbs$database) {
  writeLines(paste("Getting file for", database))
  s3Folder <- dbs$s3Folders[dbs$database == database]
  ArchiveR::restoreFile(paste0(s3Folder, "outcomeModelReference.rds"), localFolder = localFolder)
  omr <- readRDS(file.path(localFolder, "outcomeModelReference.rds")) 
  unlink(file.path(localFolder, "outcomeModelReference.rds"))
  strataFile <- omr$strataFile[omr$targetId == example$targetId[1] &
                                   omr$comparatorId == example$comparatorId[1] &
                                   omr$outcomeId == example$outcomeId[1] &
                                   omr$analysisId == 3]
  strataFile <- basename(strataFile)
  localFile <- file.path(localFolder, sprintf("StratPop_%s.rds", database))
  ArchiveR::restoreFile(paste0(s3Folder, strataFile), localFile = localFile)
  
  # Cleanup
  x <- readRDS(localFile)
  x <- data.frame(rowId = x$rowId,
                  stratumId = x$stratumId,
                  x = x$treatment,
                  y = as.integer(x$outcomeCount != 0),
                  time = x$survivalTime)
  saveRDS(x, localFile)
}

