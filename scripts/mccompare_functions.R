# From AdriTara TEammo_functions:
teReadLib <- function(libPath, libIdentifier = NULL, species, strain){
  if(file.exists(libPath)){
    cat("Loading TE library:", libPath, "\n")
    teLib <- Biostrings::readDNAStringSet(libPath)
    #species <- libPath %>% strsplit("/", fixed = TRUE) %>% lapply("[[", 4) %>% unlist()
    #strain <- libPath %>% strsplit("/", fixed = TRUE) %>% lapply("[[", 5) %>% unlist()
    teLib@metadata$species <- species
    teLib@metadata$strain <- strain
    teLib@ranges@NAMES <- gsub(pattern = "#DNA", replacement = "#TIR", teLib@ranges@NAMES)
    teLib@ranges@NAMES <- gsub(pattern = "#Unknown", replacement = "#Unclassified", teLib@ranges@NAMES)
    teLib@ranges@NAMES <- gsub(pattern = "?", replacement = "", teLib@ranges@NAMES, fixed = TRUE)
    teLib@metadata$libID <- libIdentifier
    if(length(teLib) > 0){#prevent error when there are no sequences in the fasta file
      repType <- teLib@ranges@NAMES %>%
        strsplit("#") %>% lapply("[[", 2) %>% unlist() %>% #separate header by the # and select the second part (eliminate rnd-X...)
        strsplit(" ") %>% lapply("[[", 1) %>% unlist() #separate the second part of the header by the space and take the repetitive family type
      teOrder <- repType %>% strsplit("/", fixed = TRUE) %>% lapply("[[", 1) %>% unlist()
      teType <- ifelse(!is.na(match(teOrder, teClassification$Order))
      , yes = teClassification$AppName[match(teOrder, teClassification$Order)]
      , no = repType)
      teLib@ranges@metadata$repetitiveType <- teType
      cat(libPath, " successfully loaded\n")
    } else {
      cat("There are no sequences in the selected library: ", libPath, " -------\n")
      teLib <- NULL
    }
  } else {
    cat(libPath, " does not exist\n")
    teLib <- NULL
  }
  teLib #return the library
}

tePlotLib <- function(teLibList, libType = "Please identify the library plot!!", libIdentifier = "Please identify the library"){
  if(class(teLibList) == "list"){
    # cat("This is executed")
    inputDataComb <- lapply(teLibList, function(teLib){
      if(length(teLib) >= 1){
        inputData <- data.frame(
          element = teLib@ranges@metadata$repetitiveType
          , width = teLib@ranges@width
          , species = teLib@metadata$species, strain = teLib@metadata$strain
          , libID = teLib@metadata$libID
        )
      }else{
        cat("There are empty libraries", teLib, "\n")
        inputData <- data.frame(element = character(), width = numeric(), species = character(), strain = character(), libID = character())
      }
      inputData  
    }) %>% do.call(rbind, .)
    # cat("Data for plotting is generated")
    p <- inputDataComb %>%
      ggplot(aes(x = element)) + theme_classic() +
      geom_bar(aes(fill = libID), position = "dodge") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            ,plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
      labs(x = "", y = "Count", title = paste0(unique(inputDataComb$species), "--", unique(inputDataComb$strain))) +
      scale_fill_discrete(name = "Library type")
    ggplotly(p) %>% layout(bargap = 0.3, legend = list(orientation = 'h', x = 0.1, y = 1))
  } else {
    teLib <- teLibs
    inputData <- data.frame(
      element = teLib@ranges@metadata$repetitiveType
      , width = teLib@ranges@width
      , species = teLib@metadata$species, strain = teLib@metadata$strain
    )
    p2 <- inputData %>%
      ggplot(aes(x = element)) + theme_classic() +
      geom_bar(aes(fill = element), position = "dodge", show.legend = FALSE) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "bottom"
            ,plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
      labs(x = "", y = "Count", title = paste0(unique(inputData$species), "--", unique(inputData$strain)), subtitle = libType)
    p2
  }
}

teCurationStats <- function(teLibList){
  libNames <- c("teRepeatLib" = "Raw lib", "teCleanLib" = "Clean lib"
                , "teAutoCuratedLib" = "Curated seqs", "teCompleteModelsLib" = "Complete models"
                , "teIncompleteModelsLib" = "Incomplete models"
                , "teAutoCuratedTmp" = "Remaining to MI"
                , "te80Lib" = "808080 lib", "teNon80TmpLib" = "Remaining to MI"
                , "teNon80Lib" = "Non 808080"
                , "te80RecovLib" = "808080 recovered"
                , "teStandbyLib" = "Stand by lib"
                , "te70FiltLib" = "Filtered 707070"
                , "teIncompleteRecovLib" = "Incomplete 808080 recovered"
                , "teIncompleteNon80Tmp" = "Remaining to MI"
                , "teFinalNR" = "Final lib"
                )
  
  teStats <- sapply(teLibList, function(lib){
    list(
      nSeqs = names(lib) %>% length()
      , nComplete = names(lib)[!grepl(names(lib), pattern = "_inc|_unconfirmed")] %>% length()
      , nIncomplete = names(lib)[grepl(names(lib), pattern = "_inc|_unconfirmed")] %>% length()
    ) %>% do.call(cbind, .)
  }, USE.NAMES = TRUE, simplify = FALSE) %>% do.call(rbind.data.frame, .)
  rownames(teStats) <- libNames[rownames(teStats)]
  teStats
}


# Mine

# Load the blast results comparing curated libraries
LoadBlastComparison <- function(
  blast_out,
  blast_query_lib,
  blast_subject_lib,
  species,
  strain,
  comparison) {

  blast <- data.table::fread(blast_out, sep = "\t")

  colnames(blast) <- c(
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"
  )

  # Idenitfy TEs with no hits in either direction
  query_only <- data.frame(setdiff(names(blast_query_lib), blast$qseqid)) # Seqs in dean with no hit at all
  colnames(query_only) <- c("qseqid")

  subject_only <- data.frame(setdiff(names(blast_subject_lib), blast$sseqid)) # Seqs in marta with no hits from dean
  colnames(subject_only) <- c("sseqid")

  blast <- rbind(blast, query_only, fill = TRUE)
  blast <- rbind(blast, subject_only, fill = TRUE)

  blast$species <- species
  blast$strain <- strain
  blast$compare <- comparison

  # Split the query and ubject names 
  blast <- blast %>%
  separate_wider_regex(qseqid,
    c(q_id = "[^#]+", "#",
      q_class = "[^/]+", "/?",
      q_family = "[^/]+", "/?",
      q_subfamily = ".*"),
      cols_remove = FALSE
    ) %>%
  separate_wider_regex(sseqid,
    c(s_id = "[^#]+", "#",
      s_class = "[^/]+", "/?",
      s_family = "[^/]+", "/?",
      s_subfamily = ".*"),
      cols_remove = FALSE
    )

  blast
}