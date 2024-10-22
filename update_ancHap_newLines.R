library(tidyverse)
library(parallel)
library(MASS)
library(foreach)
library(doParallel)
source("/mnt/updateAncHap_newLines/strParser.R")
#query data from CSW
library(bigrquery)
library(entVaultR)
library(tidyverse)
library(aws.s3)

#### pre-settings ####
# genetic map, bin map

#1. load a bin index map; usually can find from Danny's haploProb data
generateNewPredMap <- function(){
  #myNewCsvRaw = get_object("tzuo/Project/AncestorHaplotypes/newAncestorHaploCaller/supportingData/predictionMap_v6.csv", bucket="breeding-scratch-space")
  predictionMap <- read.csv("/mnt/data/predictionMap_v6.csv",stringsAsFactors = FALSE)
  row.names(predictionMap) <- predictionMap$markerName
  newPredictionMap <- as.character()
  for(chr in 1:10){
    tmp <- predictionMap %>% filter(genMapChr == chr)
    tmpDf <- tmp[1:(nrow(tmp)-1),]
    newPredictionMap <- rbind(newPredictionMap,tmpDf)
  }
  row.names(newPredictionMap) <- newPredictionMap$markerName
  newPredictionMap$index <- 1:nrow(newPredictionMap)
  return(newPredictionMap)
}

newPredictionMap <- generateNewPredMap()

chrStartEnd <- newPredictionMap %>% group_by(genMapChr) %>% summarise(st = min(index),ed = max(index))
totalBins = 17116

#2. Genetic map
genMap <- MONv6 <- read.csv("/mnt/data/Maize__MONv6_0__Map_2.txt",stringsAsFactors = F)

genMap$predictMap_index <- NA
for(i in 1:nrow(genMap)){
  chr <- genMap$genMapChr[i]
  pos <- round(genMap$genMapPos[i],1)
  genMap$predictMap_index[i] <- subset(newPredictionMap,genMapChr == chr & genMapPos == pos)$index
}



# get lines with LSH data
# s3://bay-breeding-np-genome-analytics-635590664227/applications/longSharedHaplotypes/Maize/data/Maize__MONv6_0__LshFpGenos.HetSoft.index

Sys.setenv("AWS_ACCESS_KEY_ID" = "",
           "AWS_SECRET_ACCESS_KEY" = "",
           "AWS_DEFAULT_REGION" = "us-east-1",
           "AWS_SESSION_TOKEN" = "")


lshLines <- s3read_using(FUN=read.delim, object="applications/longSharedHaplotypes/Maize/data/Maize__MONv6_0__LshFpGenos.HetSoft.index",
                         bucket = "bay-breeding-np-genome-analytics-635590664227",header=F)

system("aws s3 cp /mnt/updateAncHap_newLines/toS3 s3://bay-breeding-np-genome-analytics-635590664227/data/pedigreeBasedHaplotypes/ancestralHaplotypes/ancestralCall_20241018/ --recursive")
#saveRDS(lshLines,"/mnt/updateAncHap_newLines/lshLines_20241017.rds")

# get lines with parental haplotype data
bq_auth(use_oob = TRUE) 

parHapLines = as.data.frame(pullParMetaData(minMk = 15000))
#saveRDS(parHapLines,"/mnt/updateAncHap_newLines/parHapLines_20241017.rds")

#parHapLines = parHapLines[!duplicated(parHapLines$pedigree),] # same pedigree can have different germplasm_group_ids, but same parental calls.

ancHap = readRDS("/mnt/ancestralHaploCaller_20230303/results/strCondensed_ancestralHaploData_pbh_20230303.rds")


# lines to be added
newLines = parHapLines$pedigree[parHapLines$pedigree %in% lshLines$V1]
newLines = newLines[!newLines %in% ancHap$lineName] # 10/16/2024. 6692 unique lines are new, 4330/6692 has been updated in July. 
parHapLines = parHapLines %>% filter(pedigree %in% newLines)
dupLines = parHapLines %>% filter(pedigree %in% parHapLines$pedigree[duplicated(parHapLines$pedigree)]) %>% arrange(pedigree)

# lines have added
ancHap_new = readRDS("/mnt/updateAncHap_newLines/ancHap_newLines.rds")
tmpLines = colnames(ancHap_new)

parHapLines = parHapLines %>% filter(!pedigree %in% tmpLines)

# fetch credentials
# option1: use the line below if you authenticate with user credentials
bq_auth(use_oob = TRUE) # check view and manage your data in Google BigQuery

# get ancestral haplotype for all lines
numCores <- detectCores()
numCores
registerDoParallel(numCores)


rowsWithIssue = which(str_detect(parHapLines$origin_used,"lineage, dtype")) # 195 has strange origin
newPrtData <- foreach(i = 1:nrow(parHapLines),.combine=rbind) %dopar% {
  if(i %in% rowsWithIssue){return(i)} 
  origin = parHapLines$origin_used[i]
  if(str_detect(origin,"\\\\")){origin = str_replace_all(origin,"\\\\","/")} #JIDB1845-DFK-B*3\\JIDB1845-XIW
  gid = parHapLines$germplasm_group_id[i]
  qLine = parHapLines$pedigree[i]
  totalBins = 17116
  
  outFile = paste0("/mnt/updateAncHap_newLines/tmpData2/",qLine,"__",gid,".csv")
  if(!file.exists(outFile)){
    
    
    # step0: pull parental haplotypes
    parHap = as.data.frame(pullParHap(origin,gid))
    
    if(nrow(parHap) != totalBins){return(i)}
    # initialize ancestral haplotypes and ancestralHapProb
    parHap$AncestralCall = parHap$AncestralCallProb = NA
    parHap$processed_parental_hap_cp = parHap$processed_parental_hap
    parHap$processed_parental_hap_prob_cp = parHap$processed_parental_hap_prob
    
    # set the less confidence region as the line itself and prob as 1
    lowProbBins = which(parHap$processed_parental_hap_prob < 0.9 & parHap$processed_parental_hap != "Heterozygous")
    if(length(lowProbBins) > 0){
      parHap$processed_parental_hap[lowProbBins] = qLine
      parHap$processed_parental_hap_prob[lowProbBins] = 1
    } # replace low confidence bin as the line itself
    
    # step1: pull/decode ancHap for all parents from CSW; and pull lsh 
    allPrts = unique(parHap$processed_parental_hap)
    if(!all(allPrts[!allPrts %in% c("Heterozygous","Other",qLine)] %in% ancHap$lineName)){return(i)}
    prtAncHap = as.data.frame(pullPrtAncHap(allPrts))
    ancHapTable = prtAncHap %>% dplyr::select(Line,Chromosome,Position,AncestralCall) %>% 
      spread(key = "Line",value = "AncestralCall") %>% arrange(Chromosome,Position)
    
    # step1 assign parent's ancHap to this line
    for(prt in allPrts){
      if(prt == "Heterozygous"){next} # will be worked next
      
      idx = which(parHap$processed_parental_hap == prt)
      
      if(prt == qLine | prt == "Other"){ # line itself
        parHap$AncestralCall[idx] = qLine
        parHap$AncestralCallProb[idx] = 1
      }else{ # other parental haplotypes
        tmpAnc = prtAncHap %>% filter(Line == prt) %>% arrange(Chromosome, Position)
        if(nrow(tmpAnc) == totalBins){
          parHap$AncestralCall[idx] = tmpAnc$AncestralCall[idx]
          parHap$AncestralCallProb[idx] = round(tmpAnc$AncestralCallProb[idx] * parHap$processed_parental_hap_prob[idx],4)
        }else{
          print(paste(prt,"Having issue with the number of rows"))
        }
      }
      
    }
    
    allPrts = allPrts[!allPrts %in% c("Heterozygous",qLine)]
    
    
    # step2: deal with het regions
    pedInfo = parsePedToDf(origin)
    # create a layers_new data frame which is used to call the functions writen before
    layers_new = data.frame(pathLevel = 2, directParents = paste(pedInfo$Parent,collapse = ";"),Parents = paste(allPrts,collapse = ":") )
    row.names(layers_new) = qLine
    hetIdx = which(parHap$processed_parental_hap == "Heterozygous")
    
    
    parsedHaploData <- parseHetRegions(qLine,layers_new,parHap$processed_parental_hap,ancHapTable,c(qLine,allPrts))
    prtHaplos <- parsedHaploData$prtHaplo
    astHaplos <- parsedHaploData$astHaplo
    parHap$processed_parental_hap <- prtHaplos
    
    
    # step3: using LSH to fill in with others
    lshHaplos = lshToFillOthers(parHap,qLine,allPrts)

    if(length(lshHaplos) > 0){
      prtHaplos[as.numeric(names(lshHaplos))] <- lshHaplos
      parHap$processed_parental_hap[as.numeric(names(lshHaplos))] <- lshHaplos
      parHap$processed_parental_hap_prob[as.numeric(names(lshHaplos))] <- 0.9
      
    }
    
    # step3: convert to ast hap and prob
    
    ancHapProb = prtAncHap %>% dplyr::select(Line,Chromosome,Position,AncestralCallProb) %>% 
      spread(key = "Line",value = "AncestralCallProb") %>% arrange(Chromosome,Position)
    
    astHaploData <- convertPrtToAst(qLine,prtHaplos,astHaplos,ancHapTable,allPrts,ancHapProb)
    parHap$AncestralCall <- astHaploData$astHaplo
    parHap$AncestralCallProb <- parHap$processed_parental_hap_prob * astHaploData$astProb
    
    # change column names
    parHap$processed_parental_hap = parHap$processed_parental_hap_cp
    parHap$processed_parental_hap_prob = parHap$processed_parental_hap_prob_cp
    parHap = parHap[,1:11]
    parHap = parHap[,c(1:9,11,10)]
    colnames(parHap) = c("Line","Chromosome","Position","PreprocessedParentalHap","PreprocessedParentalHapProb","ProcessedParentalHap",
                      "ProcessedParentalHapProb","midasGermplasmID","originUsed","AncestralCall","AncestralCallProb")
    parHap$Version = "maize_MONv6_20241018"
    
    # step5: save the data
    
    write.csv(tmp,outFile)
    return(i)
  }
}



#### adjust the column names ####
tmpFiles = list.files("/mnt/updateAncHap_newLines/tmpData2/")
allLines = as.character(sapply(tmpFiles,function(x){unlist(strsplit(x,"\\."))[1]}))
for(line in allLines){
  fileName = paste0(line,".rds")
  if(file.exists(paste0("/mnt/updateAncHap_newLines/tmpData2/",fileName))){
    tmp = readRDS(paste0("/mnt/updateAncHap_newLines/tmpData2/",fileName))
    tmp$processed_parental_hap = tmp$processed_parental_hap_cp
    tmp$processed_parental_hap_prob = tmp$processed_parental_hap_prob_cp
    tmp = tmp[,1:11]
    tmp = tmp[,c(1:9,11,10)]
    colnames(tmp) = c("Line","Chromosome","Position","PreprocessedParentalHap","PreprocessedParentalHapProb","ProcessedParentalHap",
                      "ProcessedParentalHapProb","midasGermplasmID","originUsed","AncestralCall","AncestralCallProb")
    tmp$Version = "maize_MONv6_20241018"
    write.csv(tmp,paste0("/mnt/updateAncHap_newLines/toS3/",line,".csv"))
  }
}


# collect all ancHap together
tmpFiles = list.files("/mnt/updateAncHap_newLines/tmpData2/")
allLines = as.character(sapply(tmpFiles,function(x){unlist(strsplit(x,"\\."))[1]}))
ancHap = as.data.frame(matrix(NA,nrow = 17116,ncol=length(allLines)))
colnames(ancHap) = allLines
qLines = as.character()
for(line in allLines){
  file = paste0(line,".rds")
  if(file %in% tmpFiles){
    tmp = readRDS(paste0("/mnt/updateAncHap_newLines/tmpData2/",line,".rds"))
    ancHap[,line] = tmp$AncestralCall
  }else{
    qLines = c(qLines,line)
  }

}

eLines = allLines[!allLines %in% qLines]
ancHap = ancHap[,eLines]
saveRDS(ancHap,"/mnt/updateAncHap_newLines/ancHap_newLines.rds")






#### analyze the data ####
pct_other = group_id
pct_other$numOfSelf_parHap = NA
pct_other$numOfSelf_ancHap = NA

for(i in 1:nrow(pct_other)){
  origin = pct_other$origin_used[i]
  gid = pct_other$germplasm_pct_other[i]
  qLine = pct_other$pedigree[i]
  outFile = paste0("/mnt/updateAncHap_newLines/tmpData/",qLine,".rds")
  if(file.exists(outFile)){
    tmp = readRDS(outFile)
    pct_other$numOfSelf_parHap[i] = length(which(tmp$processed_parental_hap == "Other" | tmp$processed_parental_hap == qLine))
    pct_other$numOfSelf_ancHap[i] = length(which(tmp$AncestralCall == "Other" | tmp$AncestralCall == qLine))
  }
}
saveRDS(pct_other,"/mnt/updateAncHap_newLines/line_info.RDS")

pct_other = readRDS("/mnt/updateAncHap_newLines/line_info.RDS")





#### function ####

list2df <- function (ll){
  ll <- ll[!unlist(lapply(ll, is.null))]
  n <- length(ll)
  out <- ll[[1]]
  if (n > 1) {
    for (i in 2:n) {
      out <- rbind(out, ll[[i]])
    }
  }
  out
}


parsePedToDf = function(ped){
  peds2 <- pedparser(ped)
  peds.list <- lapply(peds2, parentsFromOrigin1rec)
  
  peds.contribs <- lapply(peds.list, parentalContrib1)
  
  out = data.frame(1, names(peds.contribs[[1]]), 
                   peds.contribs[[1]], stringsAsFactors = F)
  names(out) <- c("Line", "Parent", "Contribution")
  return(out)
}

parsePedToDf_old = function(ped){
  peds2 <- pedparser(ped)
  peds.list <- lapply(peds2, parentsFromOrigin1rec)
  
  peds.contribs <- lapply(peds.list, parentalContrib1)
  out.list <- lapply(names(peds.contribs), function(x) {y <- data.frame(x, names(peds.contribs[[x]]), 
                                                                        peds.contribs[[x]], stringsAsFactors = F); 
  names(y) <- c("Line", "Parent", "Contribution"); 
  y})
  
  out <- list2df(out.list)
  return(out)
}


# SELECT pedigree,origin_used,germplasm_group_id FROM `bcs-breeding-datasets.breeding_genomics.parental_haplotypes_metadata_maize_monv6_0_v2` 
# where crop='Maize' and marker_count >= 15000 and passed_qc

pullParMetaData = function(minMk = 10000,project = "bcs-breeding-datasets"){
  # pull parental haplotypes meta data
  sql<- paste0("SELECT ped.pedigree,origin_used,germplasm_group_id,mk.mcs_name, mk.platform_name", 
               " FROM breeding_genomics.parental_haplotypes_metadata_maize_monv6_0_v2 as ped,",
               "unnest(marker_call_sets) as mk",
               " where crop = 'Maize' and ped.marker_count >= ",minMk, " and passed_qc")
  tb <- bq_project_query(project, sql)
  df <- bq_table_download(tb,bigint="integer64") 
  return(df)
}

pullParHap = function(origin,gid){
  # pull parental haplotypes for a single line based on its origin_used and germplasm_group_id from CSW
  project = "bcs-breeding-datasets"
  sql<- paste0("select pedigree,chromosome,position,preprocessed_parental_hap, preprocessed_parental_hap_prob,
               processed_parental_hap,processed_parental_hap_prob, midas_germplasm_id, origin_used", 
               " from breeding_genomics.parental_haplotypes_maize_monv6_0_v2",
               " where origin_used ='",origin,"'", "and germplasm_group_id='",gid,"'")
  tb <- bq_project_query(project, sql)
  df <- bq_table_download(tb,bigint="integer64") 
  df = df %>% arrange(chromosome,position)
  return(df)
}

pullPrtAncHap = function(lines){
# pull ancestral haplotypes for all parents from CSW
  project = "bcs-breeding-datasets"
  sql<- paste0("select Line,Chromosome,Position,AncestralCall,AncestralCallProb", 
               " from breeding_genomics.ancestral_haplotypes_maize_v2",
               " where Line in (%s)")
  sql <- sprintf(sql,toString(sprintf("'%s'",lines)))
  tb <- bq_project_query(project, sql)
  df <- bq_table_download(tb,bigint="integer64") 
  df = df %>% arrange(Line,Chromosome,Position) # sort by Line, Chr and Pos
  return(df)
}

pullLshHap = function(qLine,prtLines){
  # pull ancestral haplotypes for all parents from CSW
  project = "bcs-breeding-datasets"
  sql<- paste0("select line1,line2,start_index,end_index", 
               " from breeding_genomics.long_shared_haplotypes_maize",
               " where line1 in (%s) and line2 in (%s)")
  sql <- sprintf(sql,toString(sprintf("'%s'",qLine)),toString(sprintf("'%s'",prtLines)))
  tb <- bq_project_query(project, sql)
  df <- bq_table_download(tb,bigint="integer64") 
  return(as.data.frame(df))
}

assignHetOtherToHet <- function(line,prtHaplo,hetRegionData,allAstLines = allPBHLines){
  # assign regions that are ("other" in PBH & "Heterozygous") in LSH to "line_het"
  if(line %in% unique(hetRegionData$V1) & line %in% allAstLines){
    selfRegions <- which(prtHaplo == line)
    subHetRegionData <- hetRegionData %>% filter(V1 == line)
    subHetRegions <- unlist(apply(subHetRegionData,1,function(x){return(x[11]:x[12])}))
    selfHetRegions <- selfRegions[selfRegions %in% subHetRegions]
    if(length(selfHetRegions) > 0){prtHaplo[selfHetRegions] <- paste0(line,"_het")}
  }
  return(prtHaplo)
}

parseHetRegions <- function(line,layers_new,prtHaplo,prtAstData,allAstLines){
  # return parsed het regions at both parental and ancestral levels, because:
  # regions could be het at parental level but not necessary at the ancestral level; 
  
  
  # step 1. get direct parents and all parents used in Parental calls
  pathLevel <- layers_new[line,"pathLevel"]
  prts <- sort(unlist(strsplit(layers_new[line,"directParents"],";")))
  numOfPrt <- length(prts)
  allPrts <- sort(unlist(strsplit(layers_new[line,"Parents"],";")))
  
  if(numOfPrt == 0){return(list(prtHaplo=prtHaplo,astHaplo=prtHaplo))} # none direct parent is FP'ed
  
  # step 2. assign regions that are ("other" in PBH & "Heterozygous") in LSH to "line_het"
  #prtHaplo <- assignHetOtherToHet(line,prtHaplo,hetRegionData,allPBHLines)
  
  # step 3. skip lines with pathLevel ==  1; no parents
  if(pathLevel == 1){return(list(prtHaplo=prtHaplo,astHaplo=prtHaplo))} 
  
  # step 5. parse "Heterozygous" regions
  astHaplo <- prtHaplo
  if(line %in% allAstLines){
    prts <- prts[prts %in% allAstLines]
    numOfAstPrt <- length(prts)
    hetRegions <- which(prtHaplo == "Heterozygous")
    numOfHetRegion <- length(hetRegions)
    if(numOfHetRegion == 0){return(list(prtHaplo=prtHaplo,astHaplo=prtHaplo))} # if no "Heterozygous", skip;
    hetHaplos <- rep(paste0(line,"_het"),numOfHetRegion)
    hetHaploPrt <- hetHaplos
    
    #stp 5_0: skip lines that are questionable (with direct parents not in all parents, or without directParents or all directParents are not FP'ed)
    if(numOfPrt == 0 | numOfAstPrt == 0){prtHaplo[hetRegions] <- hetHaplos;return(list(prtHaplo=prtHaplo,astHaplo=prtHaplo))}
    
    #step 5_1: only one parent
    if(numOfPrt == 1){
      if(numOfAstPrt == 1){
        hetHaploData <- parseHetRegionOneParent(hetHaplos,prtAstData[hetRegions,prts],prts)
        hetHaplos <- hetHaploData$hetHaplos
        hetHaploPrt <- hetHaploData$hetHaploPrt
      }
    }
    
    #step 5_2: two parents
    if(numOfPrt == 2){
      # step 5_2_1: both parents are not FP'ed
      #if(numOfAstPrt == 0){} # do nothing
      
      # step 5_2_2: one parent is not FP'ed
      if(numOfAstPrt == 1){hetHaploData <- parseHetRegionOneParent(hetHaplos,prtAstData[hetRegions,prts],prts)}
      
      # step 5_2_3: both parents are present
      if(numOfAstPrt == 2){hetHaploData <- parseHetRegionTwoParents(hetHaplos,prtAstData[hetRegions,prts],prts)}
      
      hetHaplos <- hetHaploData$hetHaplos
      hetHaploPrt <- hetHaploData$hetHaploPrt
    }
    
    #step 5_3: more than 2 parents
    if(numOfPrt >= 3){ # at least three parents
      
      # step 5_3_1: no more than one missing
      if(numOfPrt - numOfAstPrt <= 1) {
        tmpAstTable <- prtAstData[hetRegions,prts]
        numOfUniqHaplos <- apply(tmpAstTable,1,function(x){y <- unique(unlist(strsplit(x,"&|\\|"))); return(length(y))})
        combHaplos <- apply(tmpAstTable,1,function(x){y <- sort(unique(unlist(strsplit(x,"&|\\|")))); return(paste(y,collapse = "&"))})
        
        # no missing prt
        if(numOfPrt == numOfAstPrt){
          pos <- which(numOfUniqHaplos == 1) # all prts have the same haplo;
          if(length(pos) > 0 ){hetHaplos[pos] <- combHaplos[pos];hetHaploPrt[pos] <- paste(prts,collapse = "&")}
          pos <- which(numOfUniqHaplos == 2) # two unique haplos
          if(length(pos) > 0 ){hetHaplos[pos] <- combHaplos[pos];hetHaploPrt[pos] <- paste(prts,collapse = "&")}
          #pos <- which(numOfUniqHaplos > 2) #more than 2 unique haplos; Do nothing
        }
        
        # missing 1 prt
        if(numOfPrt - numOfAstPrt == 1){ 
          pos <- which(numOfUniqHaplos == 1) # 1 unique haplo
          if(length(pos) > 0 ){hetHaplos[pos] <- paste(combHaplos[pos],"&?",sep="");hetHaploPrt[pos] <- paste(paste(prts,collapse = "|"),"&?",sep="") }
          #pos <- which(numOfUniqHaplos >= 2) # more than 1 unique haplo; Do nothing
        }
        
      }
      
      # step 5_3_2: more than one missing; Do nothing
    }
    
    # assing new values to heterozygous regions
    prtHaplo[hetRegions] <- hetHaploPrt
    astHaplo[hetRegions] <- hetHaplos
  }
  return(list(prtHaplo=prtHaplo,astHaplo=astHaplo))
}

parseHetRegionOneParent <- function(hetHaplos,prtHaplos,prts){
  # hetHaplos is the regions where the query line is "Heterozygous";
  # prtHaplo is the only parent's haplotype at the same region
  # parse hetHaplo based on the only parent's haplotype (prtHaplo)
  haploPatterns <- unique(prtHaplos)
  hetHaploPrt <- hetHaplos
  for(pat in haploPatterns){
    pos <- which(prtHaplos == pat)
    if(str_detect(pat,"\\|") ){ #5_2_2_1: parental haplo is het and with a complciated pattern;
      #layers_new$type1_1[i] <- layers_new$type1_1[i] + length(pos)
    }else if(str_detect(pat,"_het")){ #t_2_2_2: parental haplo is line_prt
      #layers_new$type1_2[i] <- layers_new$type1_2[i] + length(pos)
    }else if(str_detect(pat,"&")){ #5_2_2_3: parental haplo is A&B
      hetHaplos[pos] <- paste0("(",paste(sort(unique(unlist(strsplit(pat,"&")))),collapse="|"),")&?")
      hetHaploPrt[pos] <- paste0(prts,"&?")
      #layers_new$type1_3[i] <- layers_new$type1_3[i] + length(pos)
    }else { #5_2_2_4: parental haplo is A
      hetHaplos[pos] <- paste0(pat,"&?")
      hetHaploPrt[pos] <- paste0(prts,"&?")
      #layers_new$type1_4[i] <- layers_new$type1_4[i] + length(pos)
    }
  }
  return(list(hetHaplos=hetHaplos,hetHaploPrt = hetHaploPrt))
}

parseHetRegionTwoParents <- function(hetHaplos,tmpAstTable,prts){
  # hetHaplos is the regions where the query line is "Heterozygous";
  # tmpAst is the data frame with two parents' haplotype data
  # parse hetHaplo based on the parent's haplotype (prtHaplo)
  
  tmpAstTable$pattern <- apply(tmpAstTable,1,function(x){return(paste(x,collapse =";"))})
  haploPatterns <- unique(tmpAstTable$pattern)
  hetHaploPrt <- hetHaplos
  
  for(pat in haploPatterns){
    pos <- which(tmpAstTable$pattern == pat)
    astHaplos <- sort(unlist(strsplit(pat,";")))
    
    if(str_detect(pat,"\\|")){ # 2_2_1: at least one parental haplo is het and with a complciated pattern;
      
    }else if(str_detect(pat,"_het")){ #2_2_2: at least one parental haplo is line_prt
      
    }else if(str_count(pat,"&") == 2){ # both parents are het
      if(astHaplos[1] == astHaplos[2]){ # parents are the same haplos
        hetHaplos[pos] <- astHaplos[1]
        hetHaploPrt[pos] <- paste(prts,collapse = "&")
      }else{
        
      }
    }else if (str_count(pat,"&") == 1){ # one is homo and the other is het
      astHaplos1 <- sort(unique(unlist(strsplit(astHaplos[1],"&"))))
      astHaplos2 <- sort(unique(unlist(strsplit(astHaplos[2],"&"))))
      hetHaploPrt[pos] <- paste(prts,collapse = "&")
      if(all(astHaplos1 %in% astHaplos2) | all(astHaplos2 %in% astHaplos1)) { # (A & B, A)
        hetHaplos[pos] <- ifelse(length(astHaplos1) > length(astHaplos2), astHaplos[1], astHaplos[2])
      }else{  # 3-2: (A & B, C)
        
        hetHaplos[pos] <- ifelse(length(astHaplos1) > length(astHaplos2),
                                 paste0("(",paste(astHaplos1,collapse="|"),")&",astHaplos[2]),
                                 paste0("(",paste(astHaplos2,collapse="|"),")&",astHaplos[1]))
      }
      
    }else if(str_count(pat,"&") == 0){ # both are homo
      hetHaploPrt[pos] <- paste(prts,collapse = "&")
      if(astHaplos[1] == astHaplos[2]){ #1.Same homo; This happens sometimes in the transition.
        hetHaplos[pos] <- astHaplos[1] 
      }else{ #2. Regions are homo in both parents, but with diffrent haplo
        hetHaplos[pos] <- paste(astHaplos[1],astHaplos[2],sep="&")
      }
    }else{
      
    }
  }
  return(list(hetHaplos=hetHaplos,hetHaploPrt = hetHaploPrt))
}

convertPrtToAst <- function(line,prtHaplos,astHaplos,astData,allPrts,astHaploProb){
  # This function takes a line's parental calls and its parents's ancestral calls together,
  # to trace ancestral calls and prob for a give line.
  # eg. for a given bin, line1's parental call is A, and then A's ancestral haplotype is O, then line1's astCall is O;
  
  astProb <- rep(1,length(prtHaplos))
  for(haplo in unique(prtHaplos)){
    posIndex <- which(prtHaplos == haplo)
    
    # 1. handle heterogyzous cases (A&B, A&B&C,A&?, A|B&C,); no changes on the astHaplo since it is already parsed.
    if(str_detect(haplo,"&")){
      hetPrts <- unlist(strsplit(haplo,"&|\\|"))
      hetPrts <- hetPrts[hetPrts %in% colnames(astHaploProb)]
      if(length(hetPrts) > 1){
        astProb[posIndex] <- apply(astHaploProb[posIndex,hetPrts],1,mean)
      }else if(length(hetPrts) == 1){
        astProb[posIndex] <- astHaploProb[posIndex,hetPrts]
      }else{
        # no changes
      }
    }
    
    # 2. other extreme cases
    if(!haplo %in% c("Other",allPrts)){next;} # include het,& and any other haplos
    
    # 3. typical cases
    if(length(posIndex) > 0){
      if(haplo == "Other"){
        astHaplos[posIndex] <- line
      }else{
        posIndexHaplo <- as.character(astData[posIndex,haplo])
        posIndexIndex <- str_detect(posIndexHaplo,"&|_het|Heterozygous")
        if(length(posIndexIndex) > 0){
          astHaplos[posIndex[!posIndexIndex]] <- as.character(astData[posIndex[!posIndexIndex],haplo])
          astProb[posIndex[!posIndexIndex]] <- astHaploProb[posIndex[!posIndexIndex],haplo]
          astHaplos[posIndex[posIndexIndex]] <- line
          
        }else{
          astHaplos[posIndex] <- posIndexHaplo
          astProb[posIndex] <- astHaploProb[posIndex,haplo]
        }
        
      }
    }
    #print(haplo)
    #print(table(astHaplos))
  }
  return(list(astHaplos=astHaplos, astProb = astProb))
  
}

lshToFillOthers = function(lineProbData,line,allPrts,totalBins = 17116){
  lshHaplos = as.character()
  lineProbData$haplo <- lineProbData$AncestralCall
  lineProbData$index <- 1:totalBins
  lineProbData$genMapChr = lineProbData$chromosome
  lineProbData$genMapPos = lineProbData$position
  lineProbDataHaplo <- lineProbData %>% filter(haplo == line)
  if(nrow(lineProbDataHaplo) > 0){
    otherHaploTable <- convertOtherToHaploTable(lineProbDataHaplo)
    if(nrow(otherHaploTable) > 0){
      otherHaploTable <- otherHaploTable %>% filter((end_index - start_index) >= 10)
      if(nrow(otherHaploTable) > 0){
        lsh <- pullLshHap(line,allPrts)
        lsh <- lsh %>%
          mutate(start_index_map = genMap$predictMap_index[lsh$start_index],
                 end_index_map = genMap$predictMap_index[lsh$end_index])
        
        lshHaplos <- findLineShareOtherSimpliefied(line,otherHaploTable,lsh,allPrts,lineProbData)
        lshHaplos <- lshHaplos[lshHaplos != line] # remove itself
      }
    }
  }
  return(lshHaplos)
}

convertOtherToHaploTable <- function(lineProbDataHaplo){
  row.names(lineProbDataHaplo) <- lineProbDataHaplo$index
  haploTable <- data.frame()
  
  for (chr in unique(lineProbDataHaplo$genMapChr)){
    markerIndex <- subset(lineProbDataHaplo,genMapChr == chr)$index
    if(length(markerIndex) <= 10){next;}
    startIndex <- markerIndex[1]
    currentIndex <- markerIndex[1]
    for (i in 2:length(markerIndex)){
      if(markerIndex[i] - currentIndex != 1){
        haploTable <- rbind(haploTable,data.frame(
          chromosome = chr,
          start_index = startIndex,
          end_index = markerIndex[i - 1],
          start_pos = lineProbDataHaplo[as.character(startIndex),'genMapPos'],
          end_pos = lineProbDataHaplo[as.character(markerIndex[i-1]),'genMapPos'],
          stringsAsFactors = F
          
        ))
        startIndex <- markerIndex[i]
      }
      
      if(i == length(markerIndex)){ # last data
        haploTable <- rbind(haploTable,data.frame(
          chromosome = chr,
          start_index = startIndex,
          end_index = markerIndex[i],
          start_pos = lineProbDataHaplo[as.character(startIndex),'genMapPos'],
          end_pos = lineProbDataHaplo[as.character(markerIndex[i-1]),'genMapPos'],
          stringsAsFactors = F
        ))
      }
      currentIndex <- markerIndex[i]
    }
    
  }
  return(haploTable)
}

findLineShareOtherSimpliefied <- function(line,otherHaploTable,lsh,allPrts,lineProbData){
  # Fill the "other" regions based on LSH
  # prt line must have prob. >= 0.1 
  
  allIndex <- as.numeric()
  haplos <- as.character()
  
  for(i in 1:nrow(otherHaploTable)){
    st <- otherHaploTable$start_index[i]
    ed <- otherHaploTable$end_index[i]
    
    subLsh <- lsh %>%  filter(start_index_map <= st, end_index_map >= ed)
    
    if(nrow(subLsh) > 0){
      sharedLines <- subLsh$line2
      haplos <- c(haplos,rep(sharedLines[1],ed-st+1)) # pick the first matched line
      
    }else{
      haplos <- c(haplos,rep(line,,ed-st+1))
    }
    
    allIndex = c(allIndex,st:ed)
  }
  names(haplos) <- allIndex
  return(haplos)
}
