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


loadPlotLibs <- function(species, strain) {
    teLibs <- list()
    fileNames <- list()

    # Define the files to load
    fileNames$teOutDir <- paste("./data/MCH_output", species, strain, sep = "/")
    fileNames$teManualOutDir <- paste0(fileNames$teOutDir,"/1_MI_MCH")
    fileNames$te80OutDir <- gsub(fileNames$teNon80Lib, pattern = "complete_non808080.fa", replacement = "2_MI_MCH")
    fileNames$te80TrimDir <- paste0(fileNames$te80OutDir, "/trimming")
    fileNames$teIncompleteOutDir <- paste0(fileNames$teManualOutDir, "/3_MI_MCH")
    fileNames$teIncompleteTrimDir <- paste0(fileNames$teIncompleteOutDir, "/trimming")
    fileNames$teFinalMchDir <- paste0(fileNames$teOutDir, "/MCH_final")
    fileNames$teFinalRenamed <- paste0(fileNames$teOutDir, "/", species, "_", strain, "_curated-TE-library.fa")

    fileNames$teCleanLib <- paste0(fileNames$teOutDir, "/", strain, "-clean_families.fa")
    fileNames$teAutoCuratedLib <- paste0(fileNames$teOutDir, "/curated_sequences_NR.fa")
    fileNames$teIncompleteModelsLib <- paste0(fileNames$teManualOutDir, "/incomplete_models.fa")
    fileNames$te80Lib <- paste0(fileNames$teManualOutDir, "/complete_808080.fa")
    fileNames$teCompleteModelsLib <- paste0(fileNames$teManualOutDir, "/complete_models.fa")
    fileNames$teNon80Lib <- gsub(fileNames$teCompleteModelsLib, pattern = "complete_models.fa", replacement = "complete_non808080.fa")
    fileNames$te80RecovLib <- gsub(fileNames$teCompleteModelsLib, pattern = "complete_models.fa", replacement = "complete_808080_recovered.fa")
    fileNames$teStandbyLib <- paste0(fileNames$te80OutDir, "/standby_sequences.fa")
    fileNames$te70Lib <- gsub(pattern = "standby_sequences.fa", replacement = "standby_707070.fa", x = fileNames$teStandbyLib)
    fileNames$teIncompleteNon80 <- gsub(fileNames$teIncompleteModelsLib, pattern = "incomplete_models.fa", replacement = "incomplete_non808080.fa")
    fileNames$teIncompleteRecovLib <- gsub(x = fileNames$teIncompleteNon80, pattern = "incomplete_non808080.fa", replacement = "incomplete_808080_recovered.fa")
    fileNames$te70FiltLib <- gsub(pattern = "standby_707070.fa", replacement = "standby_707070_filtered.fa",fileNames$te70Lib)
    fileNames$teFinalNR <- paste0(fileNames$teFinalMchDir, "/final_curated_NR.fa")

    fileNames$teNon80TmpLib <- paste0(fileNames$te80OutDir, "/tmp_non808080.fa")
    fileNames$te80BlastLog <- gsub(fileNames$teCompleteModelsLib, pattern = "complete_models.fa", replacement = "blast808080.log")
    fileNames$te80AssignLog <- gsub(fileNames$teCompleteModelsLib, pattern = "complete_models.fa", replacement = "assign_families.log")
    fileNames$te80MchLog <- gsub(fileNames$teCompleteModelsLib, pattern = "complete_models.fa", replacement = "mchelper_manual.log")
    fileNames$te70InputDb <- gsub(pattern = "standby_sequences.fa", replacement = "allDatabases.clustered_merged.fa", x = fileNames$teStandbyLib)
    fileNames$te70FilterLog <- gsub(fileNames$te70Lib, pattern = "standby_707070.fa", replacement = "707070_filter.log")
    fileNames$teIncompleteNewfam <- gsub(x = fileNames$teIncompleteNon80, pattern = "incomplete_non808080.fa", replacement = "3_MI_MCH/incomplete_newfam.fa")
    fileNames$teIncompleteNon80Tmp <- gsub(pattern = "incomplete_newfam.fa", replacement = "tmp_incomplete_non808080.fa", fileNames$teIncompleteNewfam)
    
    

    # Then load them

    teLibs$teCleanLib <- teReadLib(fileNames$teCleanLib, libIdentifier = "Clean lib", species = species, strain = strain)
    teLibs$teAutoCuratedLib <- teReadLib(fileNames$teAutoCuratedLib, libIdentifier = "Curated sequences", species = species, strain = strain)
    teLibs$teCompleteModelsLib <- teReadLib(fileNames$teCompleteModelsLib, libIdentifier = "Complete models", species = species, strain = strain)
    teLibs$teIncompleteModelsLib <- teReadLib(fileNames$teIncompleteModelsLib, libIdentifier = "Incomplete models", species = species, strain = strain)
    teLibs$te80Lib <- teReadLib(fileNames$te80Lib, libIdentifier = "808080 lib", species = species, strain = strain)
    teLibs$teNon80Lib <- teReadLib(fileNames$teNon80Lib, libIdentifier = "Non 80-80-80 lib", species = species, strain = strain)
    teLibs$te80RecovLib <- teReadLib(fileNames$te80RecovLib, libIdentifier = "Complete 808080 recovered", species = species, strain = strain)
    teLibs$teStandbyLib <- teReadLib(fileNames$teStandbyLib, libIdentifier = "Stand by lib", species = species, strain = strain)
    teLibs$te70Lib <- teReadLib(fileNames$te70Lib, libIdentifier = "707070 lib", species = species, strain = strain)
    teLibs$teIncompleteRecovLib <- teReadLib(fileNames$teIncompleteRecovLib, libIdentifier = "Incomplete 808080 recovered", species = species, strain = strain)
    teLibs$te70FiltLib <- teReadLib(fileNames$te70FiltLib, libIdentifier = "Standby 707070 filtered", species = species, strain = strain)
    teLibs$teFinalNR <- teReadLib(fileNames$teFinalNR, libIdentifier = "Final curated TE library", species = species, strain = strain)

  # Plot em
    tePlotLibTotal(
      list(
        teLibs$teCleanLib,
        teLibs$teAutoCuratedLib,
        teLibs$teCompleteModelsLib,
        teLibs$teIncompleteModelsLib,
        #teLibs$te80Lib,
        #teLibs$teNon80Lib,
        teLibs$te80RecovLib,
        #teLibs$teStandbyLib,
        #teLibs$te70Lib,
        teLibs$teIncompleteRecovLib,
        teLibs$te70FiltLib,
        teLibs$teFinalNR
      )
    )
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
    ggplotly(p) %>% layout(bargap = 0.3)
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

tePlotLibTotal <- function(teLibList, libType = "Please identify the library plot!!", libIdentifier = "Please identify the library"){
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
      ggplot(aes(x = fct_infreq(libID))) + theme_classic() +
      geom_bar(aes(fill = libID), position = "dodge") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
      labs(x = "", y = "Count", title = paste0(unique(inputDataComb$species), "--", unique(inputDataComb$strain))) +
      scale_fill_discrete(name = "Library type")
    ggplotly(p) %>% layout(bargap = 0.3)
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
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "bottom",
        plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
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
  comparison,
  reduce_hits = TRUE) {

  blast <- data.table::fread(blast_out, sep = "\t")

  colnames(blast) <- c(
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"
  )
  # Reduce to best hit per query
  if (reduce_hits) {
    blast <- blast %>%
      group_by(qseqid) %>%
      slice_max(order_by = bitscore, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      as.data.frame()
  }

  # Idenitfy TEs with no hits in either direction
  query_only <- data.frame(setdiff(names(blast_query_lib), blast$qseqid)) # Seqs in query with no hit at all
  colnames(query_only) <- c("qseqid")

  subject_only <- data.frame(setdiff(names(blast_subject_lib), blast$sseqid)) # Seqs in marta with no hits from dean
  colnames(subject_only) <- c("sseqid")

  blast <- bind_rows(blast, query_only) %>%
    bind_rows(subject_only)

  blast$species <- species
  blast$strain <- strain
  blast$compare <- comparison

  # Split the query and ubject names 
  blast <- blast %>%
  separate_wider_regex(qseqid,
    c(q_id = "[^#]+", "#",
      q_class = "[^/]+", "/?",
      q_family = "[^/]*", "/?",
      q_subfamily = ".*"),
      cols_remove = FALSE
    ) %>%
  separate_wider_regex(sseqid,
    c(s_id = "[^#]+", "#",
      s_class = "[^/]+", "/?",
      s_family = "[^/]*", "/?",
      s_subfamily = ".*"),
      cols_remove = FALSE
    )
  
  blast <- blast %>%
  mutate(
    seq_match = case_when(
      between(pident, 95, 100) ~ "Perfect, 95-100",
      between(pident, 80, 94.999) ~ "Present, 80",
      between(pident, 70, 79.999) ~ "Present, 70",
      (pident <= 69.999 | is.na(pident)) & is.na(sseqid) ~ "Missing from subject lib",
      (pident <= 69.999 | is.na(pident)) & is.na(qseqid) ~ "Missing from query lib"
    ),
    seq_match_score = case_when(
      between(pident, 95, 100) ~ 4,
      between(pident, 80, 94.999) ~ 3,
      between(pident, 70, 79.999) ~ 2,
      (pident <= 69.999 | is.na(pident)) & is.na(sseqid) ~ 1,
      (pident <= 69.999 | is.na(pident)) & is.na(qseqid) ~ 0
    ),
    class_match = case_when(
        is.na(sseqid) ~ "Missing from subject lib",
        is.na(qseqid) ~ "Missing from query lib",
        q_class == s_class & q_family == s_family & q_subfamily == s_subfamily ~ "Class, Family, Subfamily",
        q_class == s_class & q_family == s_family & q_subfamily != s_subfamily ~ "Class, Family",
        q_class == s_class & q_family != s_family & q_subfamily != s_subfamily ~ "Class",
        q_class != s_class & q_family != s_family & q_subfamily != s_subfamily ~ "None"
      ),
    class_match_score = case_when(
        is.na(sseqid) ~ 0,
        is.na(qseqid) ~ 0,
        q_class == s_class & q_family == s_family & q_subfamily == s_subfamily ~ 4,
        q_class == s_class & q_family == s_family & q_subfamily != s_subfamily ~ 3,
        q_class == s_class & q_family != s_family & q_subfamily != s_subfamily ~ 2,
        q_class != s_class & q_family != s_family & q_subfamily != s_subfamily ~ 1
      )
  )

  blast <- blast %>%
    mutate(
      across(where(is.character), na_if, ""),
      across(where(is.character), replace_na, "None")
    )

  blast <- blast %>%
  mutate(
    q_classification = paste(q_class, q_family, q_subfamily, sep = "/"),
    s_classification = paste(s_class, s_family, s_subfamily, sep = "/")
  )

  blast
}

PlotBlastTileMatches <- function(blast_input) {
  df1 <- blast_input %>%
    #filter(
      #!seq_match == "Missing from query lib",
      #!seq_match == "Missing from subject lib"
    #) %>%
    group_by(
      q_classification,
      s_classification,
      seq_match, class_match,
      class_match_score, seq_match_score) %>%
    summarise(
      n = n()
      ) %>%
    ungroup() %>%
    group_by(
      q_classification) %>%
    mutate(
      percentage = n / sum(n) * 100
    ) %>%
    ungroup() %>%
    mutate(
      class_match = factor(class_match, levels = c("Missing from query lib", "Missing from subject lib", "None", "Class", "Class, Family", "Class, Family, Subfamily")),
      seq_match = factor(seq_match, levels = c("Missing from query lib", "Missing from subject lib", "Present, 70", "Present, 80", "Perfect, 95-100"))
    )

    p1 <- df1 %>%
      ggplot(
        aes(
          x = reorder(q_classification, class_match_score),
          y = reorder(paste0(seq_match_score, s_classification), class_match_score),
          fill = percentage,
          text = paste(
            "Query:", q_classification,
            "\nSubject:", s_classification,
            "\nSequence match:", seq_match,
            "\nClassification match:", class_match,
            "\n% of n that classification agreed:", percentage
            )
          )
        ) +
      geom_tile() +
      theme_bw() +
      theme(
        strip.text.y = element_text(angle = 0),
        strip.text.x = element_text(angle = 90),
        axis.text.x = element_text(angle = 90)
      ) +
      scale_fill_viridis_c(option = "plasma") +
      labs(fill = "% of n assigned")

  ggplotly(p1, tooltip = "text")
}

PlotBlastViolMatches <- function(blast_input) {
  blast_input %>%
  mutate(pident = ifelse(is.na(pident), 0, pident)) %>%
  ggplot(aes(q_classification, pident)) +
  geom_violin() +
  geom_point(aes(fill = factor(class_match_score)),
    size = 3, stroke = 0.5, color = "black", shape = 21) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  facet_wrap(~ s_classification, scales = "free_x") +
  scale_fill_viridis_d(option = "plasma")
}

PlotBlastBarMatches <- function(blast_input) {
  df <- blast_input %>%
    group_by(
        q_classification,
        seq_match
        ) %>%
      summarize(
        n = n()
        ) %>%
      mutate(
        seq_match = factor(seq_match, levels = c("Missing from query lib", "Missing from subject lib", "Present, 70", "Present, 80", "Perfect, 95-100"))
      )
  bp1 <- df %>%
    ggplot(aes(fill=seq_match, y=n, x=reorder(q_classification, n))) + 
    geom_bar(position="stack", stat="identity") +
    theme_bw() +
    scale_fill_manual(values = palette_seq_match) +
    coord_flip()

  bp2 <- df %>%
    ggplot(aes(fill=seq_match, y=n, x=reorder(q_classification, n))) + 
    geom_bar(position="fill", stat="identity") +
    theme_bw() +
    scale_fill_manual(values = palette_seq_match) +
    coord_flip()
  
  # Now plot for the missing from query broken down by classification:
  df2 <- blast_input %>%
    filter(seq_match == "Missing from query lib") %>%
    group_by(
      s_classification,
      seq_match
    ) %>%
    summarize(
      n = n()
    )

  bp3 <- df2 %>%
    ggplot(aes(fill=seq_match, y=n, x=reorder(s_classification, n))) + 
    geom_bar(position="stack", stat="identity", show.legend = FALSE) +
    theme_bw() +
    scale_fill_manual(values = palette_seq_match) +
    coord_flip()
  
  df3 <- blast_input %>%
    filter(seq_match == "Missing from subject lib") %>%
    group_by(
      q_classification,
      seq_match
    ) %>%
    summarize(
      n = n()
    )

  bp4 <- df3 %>%
    ggplot(aes(fill=seq_match, y=n, x=reorder(q_classification, n))) + 
    geom_bar(position="stack", stat="identity", show.legend = FALSE) +
    theme_bw() +
    scale_fill_manual(values = palette_seq_match) +
    coord_flip()
  
  bp1 + bp2 + bp3 + bp4 +
  plot_layout(guides = 'collect', ncol = 1) +
  plot_annotation(tag_levels = 'A')
}


summarize_libraries <- function(lib_list) {
  summary_df <- data.frame()
  
  for (lib_name in names(lib_list)) {
    lib <- lib_list[[lib_name]]
    
    length_val <- length(lib)
    
    # Optionally, parse species and person from name
    parts <- strsplit(lib_name, "_")[[1]]
    species <- parts[1]
    what <- parts[2]
    
    summary_df <- rbind(
      summary_df,
      data.frame(
        name = lib_name,
        species = species,
        what = what,
        length = length_val,
        stringsAsFactors = FALSE
      )
    )
  }
  
  return(summary_df)
}



checkSeqFate <- function(
  species,
  strain,
  MCH_output_dir,
  query_ids
) {
  teLibs <- list()
    fileNames <- list()

    # Define the files to load
    fileNames$teOutDir <- paste(MCH_output_dir, species, strain, sep = "/")
    fileNames$teManualOutDir <- paste0(fileNames$teOutDir,"/1_MI_MCH")
    fileNames$te80OutDir <- paste0(fileNames$teOutDir,"/1_MI_MCH/2_MI_MCH")
    fileNames$te80TrimDir <- paste0(fileNames$te80OutDir, "/trimming")
    fileNames$teIncompleteOutDir <- paste0(fileNames$teManualOutDir, "/3_MI_MCH")
    fileNames$teIncompleteTrimDir <- paste0(fileNames$teIncompleteOutDir, "/trimming")
    fileNames$teFinalMchDir <- paste0(fileNames$teOutDir, "/MCH_final")
    fileNames$teFinalRenamed <- paste0(fileNames$teOutDir, "/", species, "_", strain, "_curated-TE-library.fa")

    fileNames$teCleanLib <- paste0(fileNames$teOutDir, "/", strain, "-clean_families.fa")
    fileNames$teAutoCuratedLib <- paste0(fileNames$teOutDir, "/curated_sequences_NR.fa")
    fileNames$teAutoCuratedLibRed <- paste0(fileNames$teOutDir, "/curated_sequences_R.fa")
    fileNames$teIncompleteModelsLib <- paste0(fileNames$teManualOutDir, "/incomplete_models.fa")
    fileNames$te80Lib <- paste0(fileNames$teManualOutDir, "/complete_808080.fa")
    fileNames$teCompleteModelsLib <- paste0(fileNames$teManualOutDir, "/complete_models.fa")
    fileNames$teNon80Lib <- gsub(fileNames$teCompleteModelsLib, pattern = "complete_models.fa", replacement = "complete_non808080.fa")
    fileNames$te80RecovLib <- gsub(fileNames$teCompleteModelsLib, pattern = "complete_models.fa", replacement = "complete_808080_recovered.fa")
    fileNames$teStandbyLib <- paste0(fileNames$te80OutDir, "/standby_sequences.fa")
    fileNames$te70Lib <- gsub(pattern = "standby_sequences.fa", replacement = "standby_707070.fa", x = fileNames$teStandbyLib)
    fileNames$teIncompleteNon80 <- gsub(fileNames$teIncompleteModelsLib, pattern = "incomplete_models.fa", replacement = "incomplete_non808080.fa")
    fileNames$teIncompleteRecovLib <- gsub(x = fileNames$teIncompleteNon80, pattern = "incomplete_non808080.fa", replacement = "incomplete_808080_recovered.fa")
    fileNames$te70FiltLib <- gsub(pattern = "standby_707070.fa", replacement = "standby_707070_filtered.fa",fileNames$te70Lib)
    fileNames$teFinalNR <- paste0(fileNames$teFinalMchDir, "/final_curated_NR.fa")
    fileNames$teFinalR <- paste0(fileNames$teFinalMchDir, "/final_curated_R.fa")

    fileNames$teNon80TmpLib <- paste0(fileNames$te80OutDir, "/tmp_non808080.fa")
    fileNames$te80BlastLog <- gsub(fileNames$teCompleteModelsLib, pattern = "complete_models.fa", replacement = "blast808080.log")
    fileNames$te80AssignLog <- gsub(fileNames$teCompleteModelsLib, pattern = "complete_models.fa", replacement = "assign_families.log")
    fileNames$te80MchLog <- gsub(fileNames$teCompleteModelsLib, pattern = "complete_models.fa", replacement = "mchelper_manual.log")
    fileNames$te70InputDb <- gsub(pattern = "standby_sequences.fa", replacement = "allDatabases.clustered_merged.fa", x = fileNames$teStandbyLib)
    fileNames$te70FilterLog <- gsub(fileNames$te70Lib, pattern = "standby_707070.fa", replacement = "707070_filter.log")
    fileNames$teIncompleteNewfam <- gsub(x = fileNames$teIncompleteNon80, pattern = "incomplete_non808080.fa", replacement = "3_MI_MCH/incomplete_newfam.fa")
    fileNames$teIncompleteNon80Tmp <- gsub(pattern = "incomplete_newfam.fa", replacement = "tmp_incomplete_non808080.fa", fileNames$teIncompleteNewfam)

    fileNames$teFirstMIDiscard <- paste0(fileNames$teManualOutDir, "/discarded_sequences.fa")
    fileNames$teCompleteDiscard <- paste0(fileNames$teIncompleteOutDir, "/discarded_non808080.fa")
    fileNames$teIncompleteDiscard <- paste0(fileNames$te80OutDir, "/discarded_non808080.fa")


    # Then load them

    teLibs$teCleanLib <- teReadLib(fileNames$teCleanLib, libIdentifier = "Clean lib", species = species, strain = strain)
    teLibs$teAutoCuratedLib <- teReadLib(fileNames$teAutoCuratedLib, libIdentifier = "Curated sequences", species = species, strain = strain)
    teLibs$teAutoCuratedLibRed <- teReadLib(fileNames$teAutoCuratedLibRed, libIdentifier = "Curated sequences", species = species, strain = strain)
    teLibs$teCompleteModelsLib <- teReadLib(fileNames$teCompleteModelsLib, libIdentifier = "Complete models", species = species, strain = strain)
    teLibs$teIncompleteModelsLib <- teReadLib(fileNames$teIncompleteModelsLib, libIdentifier = "Incomplete models", species = species, strain = strain)
    teLibs$te80Lib <- teReadLib(fileNames$te80Lib, libIdentifier = "808080 lib", species = species, strain = strain)
    teLibs$teNon80Lib <- teReadLib(fileNames$teNon80Lib, libIdentifier = "Non 80-80-80 lib", species = species, strain = strain)
    teLibs$te80RecovLib <- teReadLib(fileNames$te80RecovLib, libIdentifier = "Complete 808080 recovered", species = species, strain = strain)
    teLibs$teStandbyLib <- teReadLib(fileNames$teStandbyLib, libIdentifier = "Stand by lib", species = species, strain = strain)
    teLibs$te70Lib <- teReadLib(fileNames$te70Lib, libIdentifier = "707070 lib", species = species, strain = strain)
    teLibs$teIncompleteRecovLib <- teReadLib(fileNames$teIncompleteRecovLib, libIdentifier = "Incomplete 808080 recovered", species = species, strain = strain)
    teLibs$te70FiltLib <- teReadLib(fileNames$te70FiltLib, libIdentifier = "Standby 707070 filtered", species = species, strain = strain)
    teLibs$teFinalNR <- teReadLib(fileNames$teFinalNR, libIdentifier = "Final curated TE library", species = species, strain = strain)
    teLibs$teFinalR <- teReadLib(fileNames$teFinalR, libIdentifier = "Final curated TE library Red", species = species, strain = strain)
    teLibs$teFirstMIDiscard <- teReadLib(fileNames$teFirstMIDiscard, libIdentifier = "First MI Discard", species = species, strain = strain)
    teLibs$teCompleteDiscard <- teReadLib(fileNames$teCompleteDiscard, libIdentifier = "808080 Discard", species = species, strain = strain)
    teLibs$teIncompleteDiscard <- teReadLib(fileNames$teIncompleteDiscard, libIdentifier = "Incomplete Discard", species = species, strain = strain)

  # Check presence of each query in each library
  query_ids_clean <- sub("#.*", "", query_ids)
  results <- lapply(teLibs, function(lib) {
    seq_ids <- names(lib)

    seq_ids_clean <- sub("#.*", "", seq_ids)

    query_ids_clean %in% seq_ids_clean
  })

  # Convert results into a data frame
  presence_df <- as.data.frame(results)
  presence_df <- cbind(seqID = query_ids_clean, presence_df)
  presence_df
}


checkSeqFate2 <- function(
  species,
  strain,
  MCH_output_dir,
  query_ids  # list of character vectors OR data frame with alt names
) {
  teLibs <- list()
  fileNames <- list()

  # ---------- define all filenames ----------
  # Define the files to load
    fileNames$teOutDir <- paste(MCH_output_dir, species, strain, sep = "/")
    fileNames$teManualOutDir <- paste0(fileNames$teOutDir,"/1_MI_MCH")
    fileNames$te80OutDir <- paste0(fileNames$teOutDir,"/1_MI_MCH/2_MI_MCH")
    fileNames$te80TrimDir <- paste0(fileNames$te80OutDir, "/trimming")
    fileNames$teIncompleteOutDir <- paste0(fileNames$teManualOutDir, "/3_MI_MCH")
    fileNames$teIncompleteTrimDir <- paste0(fileNames$teIncompleteOutDir, "/trimming")
    fileNames$teFinalMchDir <- paste0(fileNames$teOutDir, "/MCH_final")
    fileNames$teFinalRenamed <- paste0(fileNames$teOutDir, "/", species, "_", strain, "_curated-TE-library.fa")

    fileNames$teCleanLib <- paste0(fileNames$teOutDir, "/", strain, "-clean_families.fa")
    fileNames$teAutoCuratedLib <- paste0(fileNames$teOutDir, "/curated_sequences_NR.fa")
    fileNames$teAutoCuratedLibRed <- paste0(fileNames$teOutDir, "/curated_sequences_R.fa")
    fileNames$teIncompleteModelsLib <- paste0(fileNames$teManualOutDir, "/incomplete_models.fa")
    fileNames$te80Lib <- paste0(fileNames$teManualOutDir, "/complete_808080.fa")
    fileNames$teCompleteModelsLib <- paste0(fileNames$teManualOutDir, "/complete_models.fa")
    fileNames$teNon80Lib <- gsub(fileNames$teCompleteModelsLib, pattern = "complete_models.fa", replacement = "complete_non808080.fa")
    fileNames$te80RecovLib <- gsub(fileNames$teCompleteModelsLib, pattern = "complete_models.fa", replacement = "complete_808080_recovered.fa")
    fileNames$teStandbyLib <- paste0(fileNames$te80OutDir, "/standby_sequences.fa")
    fileNames$te70Lib <- gsub(pattern = "standby_sequences.fa", replacement = "standby_707070.fa", x = fileNames$teStandbyLib)
    fileNames$teIncompleteNon80 <- gsub(fileNames$teIncompleteModelsLib, pattern = "incomplete_models.fa", replacement = "incomplete_non808080.fa")
    fileNames$teIncompleteRecovLib <- gsub(x = fileNames$teIncompleteNon80, pattern = "incomplete_non808080.fa", replacement = "incomplete_808080_recovered.fa")
    fileNames$te70FiltLib <- gsub(pattern = "standby_707070.fa", replacement = "standby_707070_filtered.fa",fileNames$te70Lib)
    fileNames$teFinalNR <- paste0(fileNames$teFinalMchDir, "/final_curated_NR.fa")
    fileNames$teFinalR <- paste0(fileNames$teFinalMchDir, "/final_curated_R.fa")

    fileNames$teNon80TmpLib <- paste0(fileNames$te80OutDir, "/tmp_non808080.fa")
    fileNames$te80BlastLog <- gsub(fileNames$teCompleteModelsLib, pattern = "complete_models.fa", replacement = "blast808080.log")
    fileNames$te80AssignLog <- gsub(fileNames$teCompleteModelsLib, pattern = "complete_models.fa", replacement = "assign_families.log")
    fileNames$te80MchLog <- gsub(fileNames$teCompleteModelsLib, pattern = "complete_models.fa", replacement = "mchelper_manual.log")
    fileNames$te70InputDb <- gsub(pattern = "standby_sequences.fa", replacement = "allDatabases.clustered_merged.fa", x = fileNames$teStandbyLib)
    fileNames$te70FilterLog <- gsub(fileNames$te70Lib, pattern = "standby_707070.fa", replacement = "707070_filter.log")
    fileNames$teIncompleteNewfam <- gsub(x = fileNames$teIncompleteNon80, pattern = "incomplete_non808080.fa", replacement = "3_MI_MCH/incomplete_newfam.fa")
    fileNames$teIncompleteNon80Tmp <- gsub(pattern = "incomplete_newfam.fa", replacement = "tmp_incomplete_non808080.fa", fileNames$teIncompleteNewfam)

    fileNames$teFirstMIDiscard <- paste0(fileNames$teManualOutDir, "/discarded_sequences.fa")
    fileNames$teCompleteDiscard <- paste0(fileNames$teIncompleteOutDir, "/discarded_non808080.fa")
    fileNames$teIncompleteDiscard <- paste0(fileNames$te80OutDir, "/discarded_non808080.fa")


    # Then load them

    teLibs$teCleanLib <- teReadLib(fileNames$teCleanLib, libIdentifier = "Clean lib", species = species, strain = strain)
    teLibs$teAutoCuratedLib <- teReadLib(fileNames$teAutoCuratedLib, libIdentifier = "Curated sequences", species = species, strain = strain)
    teLibs$teAutoCuratedLibRed <- teReadLib(fileNames$teAutoCuratedLibRed, libIdentifier = "Curated sequences", species = species, strain = strain)
    teLibs$teCompleteModelsLib <- teReadLib(fileNames$teCompleteModelsLib, libIdentifier = "Complete models", species = species, strain = strain)
    teLibs$teIncompleteModelsLib <- teReadLib(fileNames$teIncompleteModelsLib, libIdentifier = "Incomplete models", species = species, strain = strain)
    teLibs$te80Lib <- teReadLib(fileNames$te80Lib, libIdentifier = "808080 lib", species = species, strain = strain)
    teLibs$teNon80Lib <- teReadLib(fileNames$teNon80Lib, libIdentifier = "Non 80-80-80 lib", species = species, strain = strain)
    teLibs$te80RecovLib <- teReadLib(fileNames$te80RecovLib, libIdentifier = "Complete 808080 recovered", species = species, strain = strain)
    teLibs$teStandbyLib <- teReadLib(fileNames$teStandbyLib, libIdentifier = "Stand by lib", species = species, strain = strain)
    teLibs$te70Lib <- teReadLib(fileNames$te70Lib, libIdentifier = "707070 lib", species = species, strain = strain)
    teLibs$teIncompleteRecovLib <- teReadLib(fileNames$teIncompleteRecovLib, libIdentifier = "Incomplete 808080 recovered", species = species, strain = strain)
    teLibs$te70FiltLib <- teReadLib(fileNames$te70FiltLib, libIdentifier = "Standby 707070 filtered", species = species, strain = strain)
    teLibs$teFinalNR <- teReadLib(fileNames$teFinalNR, libIdentifier = "Final curated TE library", species = species, strain = strain)
    teLibs$teFinalR <- teReadLib(fileNames$teFinalR, libIdentifier = "Final curated TE library Red", species = species, strain = strain)
    teLibs$teFirstMIDiscard <- teReadLib(fileNames$teFirstMIDiscard, libIdentifier = "First MI Discard", species = species, strain = strain)
    teLibs$teCompleteDiscard <- teReadLib(fileNames$teCompleteDiscard, libIdentifier = "808080 Discard", species = species, strain = strain)
    teLibs$teIncompleteDiscard <- teReadLib(fileNames$teIncompleteDiscard, libIdentifier = "Incomplete Discard", species = species, strain = strain)

  # ---------- normalise query IDs ----------
  clean_names <- function(x) sub("#.*", "", x)

  # If data frame, turn to list of character vectors
  if (is.data.frame(query_ids)) {
    query_list <- split(query_ids, seq_len(nrow(query_ids)))
    query_list <- lapply(query_list, function(row) clean_names(unlist(row)))
  } else if (is.list(query_ids)) {
    query_list <- lapply(query_ids, clean_names)
  } else {
    query_list <- as.list(clean_names(query_ids))
  }

  # ---------- check presence ----------
  results <- lapply(teLibs, function(lib) {
    lib_ids_clean <- clean_names(names(lib))
    sapply(query_list, function(possible_names) {
      any(possible_names %in% lib_ids_clean)
    })
  })

  # ---------- make data frame ----------
  presence_df <- as.data.frame(results)
  presence_df <- cbind(
  seqID = sapply(query_list, function(x) paste(x, collapse = "|")), 
  presence_df
)
  presence_df
}

# Get TE names for those missing from query and subject
GetMissingTEs <- function(
  genome_name = "Genome",
  names_relate_file_Query,
  blast_Genome_Lib1_vs_Lib2,
  blast_Genome_Lib2_vs_Lib1,
  blast_Genome_Lib2_vs_MCH,
  lib1_name = "Lib1",
  lib2_name = "Lib2"
){
names_relate <- data.table::fread(names_relate_file_Query, sep = "\t") %>%
  select(seqID, newIDs)

Lib2_MCH_names <- blast_Genome_Lib2_vs_MCH %>%
      filter(qseqid != "None") %>%
      group_by(sseqid) %>%
      slice_max(order_by = bitscore, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
  select(
    qseqid, sseqid
  ) %>%
  as.data.frame()

  # Missing from subject
  missing_from_Lib1 <- blast_Genome_Lib2_vs_Lib1 %>%
    filter(
      seq_match == "Missing from subject lib") %>%
      select(qseqid) %>%
    as.data.frame()
  
  missing_from_Lib1_recip <- blast_Genome_Lib1_vs_Lib2 %>%
    filter(
      seq_match == "Missing from query lib") %>%
      select(qseqid = sseqid) %>%
    as.data.frame()
  
  missing_from_Lib1 <- bind_rows(missing_from_Lib1, missing_from_Lib1_recip) %>%
    distinct()

  missing_from_Lib1_final <- missing_from_Lib1 %>%
    left_join(Lib2_MCH_names, by = c("qseqid" = "qseqid")) %>%
    filter(!is.na(sseqid)) %>%
    select(sseqid)

  # Missing from query
  missing_from_Lib2 <- blast_Genome_Lib2_vs_Lib1 %>%
    filter(
      seq_match == "Missing from query lib") %>%
      select(sseqid) %>%
    as.data.frame()
  
  missing_from_Lib2_recip <- blast_Genome_Lib1_vs_Lib2 %>%
    filter(
      seq_match == "Missing from subject lib") %>%
      select(sseqid = qseqid) %>%
    as.data.frame()
  
  missing_from_Lib2 <- bind_rows(missing_from_Lib2, missing_from_Lib2_recip) %>%
    distinct()

  # Now join in the original MCH names - seqID
  missing_from_Lib2_final <- missing_from_Lib2 %>%
    left_join(names_relate, by = c("sseqid" = "newIDs")) %>%
    filter(!is.na(seqID)) %>%
    select(seqID)

  write.csv(
    missing_from_Lib1_final,
    paste0("results/missing_from_Lib1_", lib1_name, "_", genome_name, ".csv"),
    row.names = FALSE, quote = FALSE
  )

  write.csv(
    missing_from_Lib2_final,
    paste0("results/missing_from_Lib2_", lib2_name, "_", genome_name, ".csv"),
    row.names = FALSE, quote = FALSE
)
}