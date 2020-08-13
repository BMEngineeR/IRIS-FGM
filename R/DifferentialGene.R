#' @include generics.R
#' @include Classes.R
NULL

#' @param object BRIC object
#'
#' @param SimpleResult marker gene only output log fold change (LFC), p-value, and adjusted p-value.
#' @param FDR a number to specify the threshold of FDR, default by 0.05
#' @name FindMarkers
#' @importFrom DEsingle DEsingle DEtype
#'
.findMarkers <- function(object,SimpleResult = T, FDR = 0.05){
  # two group number as factor.
  message("select condition to compare")
  message(paste0(c(1:ncol(object@MetaInfo))," : ",c(colnames(object@MetaInfo)),"\n"))
  ident.index <-  readline(prompt="select index of cell condition: ")
  ident.index <- as.numeric(ident.index)
  tmp.ident <- object@MetaInfo[,ident.index]
  names(tmp.ident) <- rownames(object@MetaInfo)
  label.used <- colnames(object@MetaInfo)[ident.index]
  # create index table
  tmp.group.table <- data.frame(index = 1:length(unique(tmp.ident)), condition = as.character(sort(unique(tmp.ident))),stringsAsFactors = F)
  tmp.group.table <- rbind(tmp.group.table, c(nrow(tmp.group.table)+1,"rest of all"))
  # select groups to compare
  message("select index (left) of first group to compare : ")
  message(paste0(tmp.group.table$index[1:nrow(tmp.group.table)-1], " : ", tmp.group.table$condition[1:nrow(tmp.group.table)-1],"\n"))
  group.1.idx <- readline("input first group index : ")
  message("select index (left) of second group to compare : ")
  message(paste0(tmp.group.table$index[tmp.group.table$index != group.1.idx], " : ", tmp.group.table$condition[tmp.group.table$index != group.1.idx],"\n"))
  group.2.idx <- readline(prompt="select index of group 2: ")
  group.1 <- tmp.group.table[tmp.group.table$index == group.1.idx, 2]
  group.2 <- tmp.group.table[tmp.group.table$index == group.2.idx, 2]
  tmp.expression.table <- object@LTMG@LTMG_discrete
  if(group.2 == "rest of all") {
    new.condition <- as.character(tmp.ident)
    # set group in new condition
    new.condition[tmp.ident == group.1] <- group.1
    new.condition[tmp.ident != group.1] <- "all"
    results <- DEsingle(counts = tmp.expression.table, group = as.factor(new.condition))
  } else {
    new.condition <- vector(mode = "logical", length = length(tmp.ident))
    new.condition[tmp.ident == group.1] <- "TRUE"
    new.condition[tmp.ident == group.2] <- "TRUE"
    tmp.new.condition <- as.factor(as.character(tmp.ident)[new.condition == "TRUE"])
    tmp.expression.table <- tmp.expression.table[,new.condition == "TRUE"]
    results <- DEsingle(counts = tmp.expression.table, group = as.factor(tmp.new.condition))
  }
  results.classified <- DEtype(results = results, threshold = FDR)
  if (SimpleResult == TRUE) {
    gene.name <- rownames(results.classified)
    results.classified <- cbind(log(results.classified$foldChange),results.classified$pvalue,results.classified$pvalue.adj.FDR )
    colnames(results.classified) <- c("LFC","pval","pvalue.adj.FDR")
    rownames(results.classified) <- gene.name
    results.classified <- as.data.frame(results.classified)
    results.classified <- results.classified[results.classified$pvalue.adj.FDR< FDR,]
  }
  if(grepl("MC", label.used, ignore.case = T)){
    object@BiCluster@MarkerGene <- results.classified
  } else {
    object@LTMG@MarkerGene <- results.classified
  }

  return(object)
}

#' @rdname FindMarkers
#' @export
setMethod("FindMarkers", "BRIC", .findMarkers)

