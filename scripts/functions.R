################################################################################
# GENERAL FUNCTIONS
################################################################################

# Auto installer
install.packages.auto <- function(x) { 
  x <- as.character(substitute(x)) 
  if(isTRUE(x %in% .packages(all.available = TRUE))) { 
    eval(parse(text = sprintf("require(\"%s\")", x)))
  } else { 
    # Update installed packages - this may mean a full upgrade of R, which in turn
    # may not be warrented. 
    #update.install.packages.auto(ask = FALSE) 
    eval(parse(text = sprintf("install.packages(\"%s\", dependencies = TRUE, repos = \"https://cloud.r-project.org/\")", x)))
  }
  if(isTRUE(x %in% .packages(all.available = TRUE))) { 
    eval(parse(text = sprintf("require(\"%s\")", x)))
  } else {
    if (!requireNamespace("BiocManager"))
      install.packages("BiocManager")
    BiocManager::install() # this would entail updating installed packages, which in turned may not be warrented
    
    # Code for older versions of R (<3.5.0)
    # source("http://bioconductor.org/biocLite.R")
    # Update installed packages - this may mean a full upgrade of R, which in turn
    # may not be warrented.
    # biocLite(character(), ask = FALSE) 
    eval(parse(text = sprintf("BiocManager::install(\"%s\")", x)))
    eval(parse(text = sprintf("require(\"%s\")", x)))
  }
}


################################################################################
# (GENERAL) LINEAR/LOGISTIC MODELING
################################################################################
# Function to grep data from glm()/lm()

### CONTINUOUS TRAITS
GLM.CON <- function(fit, DATASET, x_name, y, verbose=c(TRUE,FALSE)){
  cat("Analyzing in dataset '", DATASET ,"' the association of '", x_name ,"' with '", y ,"' .\n")
  if (nrow(summary(fit)$coefficients) == 1) {
    output = c(DATASET, x_name, y, rep(NA,8))
    cat("Model not fitted; probably singular.\n")
  }else {
    cat("Collecting data.\n\n")
    effectsize = summary(fit)$coefficients[2,1]
    SE = summary(fit)$coefficients[2,2]
    OReffect = exp(summary(fit)$coefficients[2,1])
    CI_low = exp(effectsize - 1.96 * SE)
    CI_up = exp(effectsize + 1.96 * SE)
    tvalue = summary(fit)$coefficients[2,3]
    pvalue = summary(fit)$coefficients[2,4]
    R = summary(fit)$r.squared
    R.adj = summary(fit)$adj.r.squared
    sample_size = nrow(model.frame(fit))
    N = study.samplesize
    Perc_Miss = 100 - ((sample_size * 100)/N)
    
    output = c(DATASET, x_name, y, effectsize, SE, OReffect, CI_low, CI_up, tvalue, pvalue, R, R.adj, N, sample_size, Perc_Miss)
    
    if (verbose == TRUE) {
      cat("We have collected the following and summarize it in an object:\n")
      cat("Dataset...................:", DATASET, "\n")
      cat("Score/Exposure/biomarker..:", x_name, "\n")
      cat("Trait/outcome.............:", y, "\n")
      cat("Effect size...............:", round(effectsize, 6), "\n")
      cat("Standard error............:", round(SE, 6), "\n")
      cat("Odds ratio (effect size)..:", round(OReffect, 3), "\n")
      cat("Lower 95% CI..............:", round(CI_low, 3), "\n")
      cat("Upper 95% CI..............:", round(CI_up, 3), "\n")
      cat("T-value...................:", round(tvalue, 6), "\n")
      cat("P-value...................:", signif(pvalue, 8), "\n")
      cat("R^2.......................:", round(R, 6), "\n")
      cat("Adjusted r^2..............:", round(R.adj, 6), "\n")
      cat("Sample size of AE DB......:", N, "\n")
      cat("Sample size of model......:", sample_size, "\n")
      cat("Missing data %............:", round(Perc_Miss, 6), "\n")
    } else {
      cat("Collecting data in summary object.\n")
    }
  }
  return(output)
  print(output)
}

### BINARY TRAITS
GLM.BIN <- function(fit, DATASET, x_name, y, verbose=c(TRUE,FALSE)){
  cat("Analyzing in dataset '", DATASET ,"' the association of '", x_name ,"' with '", y ,"' ...\n")
  if (nrow(summary(fit)$coefficients) == 1) {
    output = c(DATASET, x_name, y, rep(NA,9))
    cat("Model not fitted; probably singular.\n")
  }else {
    cat("Collecting data...\n")
    effectsize = summary(fit)$coefficients[2,1]
    SE = summary(fit)$coefficients[2,2]
    OReffect = exp(summary(fit)$coefficients[2,1])
    CI_low = exp(effectsize - 1.96 * SE)
    CI_up = exp(effectsize + 1.96 * SE)
    zvalue = summary(fit)$coefficients[2,3]
    pvalue = summary(fit)$coefficients[2,4]
    dev <- fit$deviance
    nullDev <- fit$null.deviance
    modelN <- length(fit$fitted.values)
    R.l <- 1 - dev / nullDev
    R.cs <- 1 - exp(-(nullDev - dev) / modelN)
    R.n <- R.cs / (1 - (exp(-nullDev/modelN)))
    sample_size = nrow(model.frame(fit))
    N = study.samplesize
    Perc_Miss = 100 - ((sample_size * 100)/N)
    
    output = c(DATASET, x_name, y, effectsize, SE, OReffect, CI_low, CI_up, zvalue, pvalue, R.l, R.cs, R.n, N, sample_size, Perc_Miss)
    if (verbose == TRUE) {
      cat("We have collected the following and summarize it in an object:\n")
      cat("Dataset...................:", DATASET, "\n")
      cat("Score/Exposure/biomarker..:", x_name, "\n")
      cat("Trait/outcome.............:", y, "\n")
      cat("Effect size...............:", round(effectsize, 6), "\n")
      cat("Standard error............:", round(SE, 6), "\n")
      cat("Odds ratio (effect size)..:", round(OReffect, 3), "\n")
      cat("Lower 95% CI..............:", round(CI_low, 3), "\n")
      cat("Upper 95% CI..............:", round(CI_up, 3), "\n")
      cat("Z-value...................:", round(zvalue, 6), "\n")
      cat("P-value...................:", signif(pvalue, 8), "\n")
      cat("Hosmer and Lemeshow r^2...:", round(R.l, 6), "\n")
      cat("Cox and Snell r^2.........:", round(R.cs, 6), "\n")
      cat("Nagelkerke's pseudo r^2...:", round(R.n, 6), "\n")
      cat("Sample size of AE DB......:", N, "\n")
      cat("Sample size of model......:", sample_size, "\n")
      cat("Missing data %............:", round(Perc_Miss, 6), "\n")
    } else {
      cat("Collecting data in summary object.\n")
    }
  }
  return(output)
  print(output)
}


################################################################################
# REGIONAL ASSOCIATION PLOTTING
################################################################################


# RACER singleRegionalAssocPlot
singlePlotRACER2 <- function (assoc_data, chr, build = "hg19", set = "protein_coding", 
                              plotby, gene_plot = NULL, snp_plot = NULL, start_plot = NULL, 
                              end_plot = NULL, label_lead = FALSE, 
                              grey_colors = FALSE, 
                              cred_set = FALSE, 
                              gene_track_h = 3, gene_name_s = 2.5) {
  if (missing(assoc_data)) {
    stop("Please provide a data set to plot.")
  }
  else if (missing(chr)) {
    stop("Please specify which chromosome you wish to plot.")
  }
  else if (missing(plotby)) {
    stop("Please specify the method by which you wish to plot.")
  }
  else if (plotby == "gene") {
    if (is.null(gene_plot)) {
      stop("Please specify a gene to plot by.")
    }
  }
  else if (plotby == "snp") {
    if (is.null(snp_plot)) {
      stop("Please specify a snp to plot by.")
    }
  }
  else if (plotby == "coord") {
    if (is.null(start_plot) | is.null(end_plot)) {
      stop("Please specify start coordinate for plot.")
    }
  }
  else {
    message("All inputs are go.")
  }
  reqs = c("CHR", "POS", "LOG10P")
  cols = colnames(assoc_data)
  if (sum(reqs %in% cols) == 3) {
  }
  else {
    stop("Association Data Set is missing a required column, please format your data set using formatRACER.R.")
  }
  reqs_2 = c("LD", "LD_BIN")
  if (sum(reqs_2 %in% cols) == 2) {
  }
  else {
    message("Association Data Set is missing LD data, the resulting plot won't have LD information, but you can add it using the ldRACER.R function.")
  }
  `%>%` <- magrittr::`%>%`
  if (build == "hg38") {
    utils::data(hg38)
    chr_in = chr
    colnames(hg38) = c("GENE_ID", "CHR", "TRX_START", "TRX_END", 
                       "LENGTH", "GENE_NAME", "TYPE")
    gene_sub = hg38[hg38$CHR == chr_in, ]
  }
  else if (build == "hg19") {
    utils::data(hg19)
    chr_in = chr
    colnames(hg19) = c("GENE_ID", "CHR", "TRX_START", "TRX_END", 
                       "LENGTH", "GENE_NAME", "TYPE")
    gene_sub = hg19[hg19$CHR == chr_in, ]
  }
  if (set == "protein_coding") {
    gene_sub = gene_sub[gene_sub$TYPE == "protein_coding", 
    ]
  }
  else {
    gene_sub = gene_sub
  }
  if (sum(is.null(plotby)) == 1) {
    stop("Please specify a method by which to plot.")
  }
  if (sum(is.null(plotby)) == 0) {
    message("Plotting by...")
    if ((plotby == "coord") == TRUE) {
      message("coord")
      start = start_plot
      end = end_plot
    }
    else if ((plotby == "gene") == TRUE) {
      message(paste("gene:", gene_plot))
      if (sum(is.null(gene_plot)) == 0) {
        p = subset(gene_sub, gene_sub$GENE_NAME == gene_plot)
        start = min(p$TRX_START) - 5e+05
        end = max(p$TRX_END) + 5e+05
      }
      else {
        message("No gene specified.")
      }
    }
    else if ((plotby == "snp") == TRUE) {
      message(paste("snp", snp_plot))
      q = assoc_data[assoc_data$RS_ID == snp_plot, ]
      w = q$POS
      w = as.numeric(as.character(w))
      start = w - 5e+05
      end = w + 5e+05
    }
  }
  gene_sub = subset(gene_sub, gene_sub$TRX_START > (start - 
                                                      50000))
  gene_sub = subset(gene_sub, gene_sub$TRX_END < (end + 50000))
  gene_sub = gene_sub[, c(3, 4, 6)]
  gene_sub = reshape2::melt(gene_sub, id.vars = "GENE_NAME")
  gene_sub$y_value = as.numeric(as.factor(gene_sub$GENE_NAME))
  plot_lab = subset(gene_sub, gene_sub$variable == "TRX_END")
  message("Reading in association data")
  in.dt <- as.data.frame(assoc_data)
  in.dt$POS = as.numeric(as.character(in.dt$POS))
  in.dt$LOG10P = as.numeric(as.character(in.dt$LOG10P))
  in.dt$CHR = as.numeric(as.character(in.dt$CHR))
  in.dt = dplyr::filter(in.dt, .data$CHR == chr_in)
  in.dt = dplyr::filter(in.dt, .data$POS > start) %>% dplyr::filter(.data$POS < 
                                                                      end)
  if (label_lead == TRUE) {
    message("Determining lead SNP")
    lsnp_row = which(in.dt$LABEL == "LEAD")
    label_data = in.dt[lsnp_row, ]
    if (dim(label_data)[1] == 0) {
      lsnp_row = in.dt[in.dt$LOG10P == max(in.dt$LOG10P), ]
      label_data = lsnp_row[1, ]
    }
  }
  
  if (cred_set == TRUE) {
    message("Collecting posterior probabilities")
    ppsnp_row = which(in.dt$Posterior_Prob >= 0)
    pp_data = in.dt[ppsnp_row, ]
    if (dim(pp_data)[1] == 0) {
      ppsnp_row = in.dt[in.dt$LOG10P == max(in.dt$LOG10P), ]
      pp_data = ppsnp_row[1, ]
    }
  }
  
  
  message("Generating Plot")
  if ("LD" %in% colnames(in.dt) && "LD_BIN" %in% colnames(in.dt)) {
    c = ggplot2::ggplot(gene_sub, ggplot2::aes_string(x = "value", 
                                                      y = "y_value")) + ggplot2::geom_line(ggplot2::aes_string(group = "GENE_NAME"), 
                                                                                           size = 2) + ggplot2::theme_minimal() + 
      ggplot2::geom_text(data = plot_lab, 
                         ggplot2::aes_string(x = "value", y = "y_value", label = "GENE_NAME"), 
                         hjust = -0.1, vjust = 0.3,
                         size = gene_name_s) + 
      ggplot2::theme(axis.title.y = ggplot2::element_text(color = "transparent", size = 28),
                     axis.text.y = ggplot2::element_blank(), 
                     axis.ticks.y = ggplot2::element_blank(),
                     panel.border = element_blank(), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line.x = element_line(colour = "#595A5C")) + 
      ggplot2::xlab(paste0("chromosome ", 
                           chr_in, " position")) + ggplot2::coord_cartesian(xlim = c(start, 
                                                                                     end), ylim = c(0, (max(gene_sub$y_value) + 0.25)))
    b = ggplot2::ggplot(in.dt, ggplot2::aes_string(x = "POS", 
                                                   y = "LOG10P", color = "LD_BIN")) + ggplot2::geom_point() + 
      ggplot2::scale_colour_manual(values = c(`1.0-0.8` = "#DC0000FF", # "#DB003F", #"red", 
                                              `0.8-0.6` = "#F39B7FFF", # "#F59D10", #"darkorange1", 
                                              `0.6-0.4` = "#00A087FF", # "#49A01D", #"green1", 
                                              `0.4-0.2` = "#4DBBD5FF", # "#1290D9", #"skyblue1", 
                                              `0.2-0.0` = "#3C5488FF", # "#4C81BF", #"navyblue", 
                                              `NA` = "#A2A3A4" # "grey"
      ), 
      drop = FALSE) + 
      labs(color = bquote(LD~r^2)) + 
      ggplot2::theme(panel.border = element_blank(), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line.x = element_blank(),
                     axis.line.y = element_line(colour = "#595A5C"),
                     axis.title.x=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank()) +  
      # ggplot2::xlab("Chromosome Position") + 
      # ggplot2::ylab("-log10(p-value)") + 
      ggplot2::ylab(bquote(-log[10]~(p-value))) + 
      ggplot2::coord_cartesian(xlim = c(start, end), ylim = c(min(in.dt$LOG10P), 
                                                              max(in.dt$LOG10P)))
  }
  else {
    c = ggplot2::ggplot(gene_sub, ggplot2::aes_string(x = "value", 
                                                      y = "y_value")) + ggplot2::geom_line(ggplot2::aes_string(group = "GENE_NAME"), 
                                                                                           size = 2) + ggplot2::theme_minimal() + 
      ggplot2::geom_text(data = plot_lab, 
                         ggplot2::aes_string(x = "value", y = "y_value", label = "GENE_NAME"), 
                         hjust = -0.1, vjust = 0.3,
                         size = gene_name_s) + 
      ggplot2::theme(axis.title.y = ggplot2::element_blank(), #ggplot2::element_text(color = "white", size = 28),  
                     axis.text.y = ggplot2::element_blank(), 
                     axis.ticks.y = ggplot2::element_blank(),
                     panel.border = element_blank(), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line.x = element_line(colour = "#595A5C")) + 
      ggplot2::xlab(paste0("chromosome ", 
                           chr_in, " position")) + ggplot2::coord_cartesian(xlim = c(start, 
                                                                                     end), ylim = c(0, (max(gene_sub$y_value) + 0.25)))
    b = ggplot2::ggplot(in.dt, ggplot2::aes_string(x = "POS", 
                                                   y = "LOG10P")) + ggplot2::geom_point() + ggplot2::theme(panel.border = element_blank(), 
                                                                                                           panel.grid.major = element_blank(),
                                                                                                           panel.grid.minor = element_blank(), 
                                                                                                           axis.line.x = element_blank(),
                                                                                                           axis.line.y = element_line(colour = "#595A5C"),
                                                                                                           axis.title.x=element_blank(),
                                                                                                           axis.text.x=element_blank(),
                                                                                                           axis.ticks.x=element_blank()) + 
      # ggplot2::xlab("Chromosome Position") + 
      # ggplot2::ylab("-log10(p-value)") + 
      ggplot2::ylab(bquote(-log[10]~(p-value))) + 
      ggplot2::coord_cartesian(xlim = c(start, end), ylim = c(min(in.dt$LOG10P), 
                                                              max(in.dt$LOG10P)))
  }
  if (label_lead == TRUE) {
    b = b + geom_point(data = label_data, 
                       aes_string(x = "POS", 
                                  y = "LOG10P"), 
                       color = "#8491B4FF", # "#9A3480", #"purple")
                       fill = "#8491B4FF", # "#9A3480",
                       size = 4, shape = 23) 
    
    b = b + geom_text(data = label_data, 
                      aes_string(label = "RS_ID"), 
                      color = "black", 
                      size = 4, hjust = 1.25)
  }
  if (grey_colors == TRUE) {
    b = b + geom_point(color = "#A2A3A4", fill = "#A2A3A4")
  }
  
  if (cred_set == TRUE) {
    
    b = b + geom_point(data = pp_data, 
                       aes_string(x = "POS", 
                                  y = "LOG10P")) +
      geom_point(aes(colour = Posterior_Prob)) +
      scale_colour_gradient(
        low = "#F39B7FFF", # "#132B43",
        high = "#DC0000FF", # "#56B1F7",
        space = "Lab",
        na.value = "#A2A3A4", # "grey50",
        guide = "colourbar",
        aesthetics = "colour"
      ) + scale_fill_gradient(
        low = "#F39B7FFF", # "#132B43",
        high = "#DC0000FF", # "#56B1F7",
        space = "Lab",
        na.value = "#A2A3A4", # "grey50",
        guide = "colourbar",
        aesthetics = "colour"
      )
    
  }
  
  ggpubr::ggarrange(b, c, heights = c(gene_track_h, 1), nrow = 2, ncol = 1, 
                    common.legend = TRUE, legend = "right")
}

# RACER mirrorAssociationPlot
mirrorPlotRACER2 <- function (assoc_data1, assoc_data2, chr, build = "hg19", set = "protein_coding", 
          name1 = "Association Dataset #1", name2 = "Association Dataset #2", 
          plotby, gene_plot = NULL, snp_plot = NULL, start_plot = NULL, 
          end_plot = NULL, label_lead = FALSE) {
  reqs = c("CHR", "POS", "LOG10P")
  cols_1 = colnames(assoc_data1)
  cols_2 = colnames(assoc_data2)
  if (sum(reqs %in% cols_1) == 3) {
  }
  else {
    stop("Association Data Set #1 is missing a required column.")
  }
  if (sum(reqs %in% cols_2) == 3) {
  }
  else {
    stop("Association Data Set #2 is missing a required column.")
  }
  if (build == "hg38") {
    utils::data(hg38)
    chr_in = chr
    colnames(hg38) = c("GENE_ID", "CHR", "TRX_START", "TRX_END", 
                       "LENGTH", "GENE_NAME", "TYPE")
    gene_sub = hg38[hg38$CHR == chr_in, ]
  }
  else if (build == "hg19") {
    utils::data(hg19)
    chr_in = chr
    colnames(hg19) = c("GENE_ID", "CHR", "TRX_START", "TRX_END", 
                       "LENGTH", "GENE_NAME", "TYPE")
    gene_sub = hg19[hg19$CHR == chr_in, ]
  }
  if (set == "protein_coding") {
    gene_sub = gene_sub[gene_sub$TYPE == "protein_coding", 
    ]
  }
  else {
    gene_sub = gene_sub
  }
  `%>%` <- magrittr::`%>%`
  if ((sum(is.null(plotby)) == 0) == TRUE) {
    message("Plotting by...")
    if ((plotby == "coord") == TRUE) {
      message("coord")
      start = start_plot
      end = end_plot
    }
    else if ((plotby == "gene") == TRUE) {
      message(paste("gene:", gene_plot))
      if (sum(is.null(gene_plot)) == 0) {
        print(head(gene_sub))
        p = subset(gene_sub, gene_sub$GENE_NAME == gene_plot)
        start = min(p$TRX_START) - 5e+05
        end = max(p$TRX_END) + 5e+05
      }
      else {
        stop("No gene specified.")
      }
    }
    else if ((plotby == "snp") == TRUE) {
      message(paste("snp", snp_plot))
      q = assoc_data1[assoc_data1$RS_ID == snp_plot, ]
      w = q$POS
      w = as.numeric(as.character(w))
      start = w - 5e+05
      end = w + 5e+05
    }
  }
  else {
    stop("Please specify a parameter to plotby.")
  }
  gene_sub = subset(gene_sub, gene_sub$TRX_START > (start - 
                                                      5000))
  gene_sub = subset(gene_sub, gene_sub$TRX_END < (end + 5000))
  gene_sub = gene_sub[!duplicated(gene_sub$GENE_ID), ]
  gene_sub = gene_sub[, c(3, 4, 6)]
  gene_sub = reshape2::melt(gene_sub, id.vars = "GENE_NAME")
  gene_sub$y_value = as.numeric(as.factor(gene_sub$GENE_NAME))
  plot_lab = subset(gene_sub, gene_sub$variable == "TRX_END")
  message("Reading in association data")
  in.dt <- as.data.frame(assoc_data1)
  in.dt$POS = as.numeric(as.character(in.dt$POS))
  in.dt$LOG10P = as.numeric(as.character(in.dt$LOG10P))
  in.dt$CHR = as.numeric(as.character(in.dt$CHR))
  in.dt = dplyr::filter(in.dt, .data$CHR == chr_in)
  in.dt = dplyr::filter(in.dt, .data$POS > start) %>% dplyr::filter(.data$POS < 
                                                                      end)
  if (label_lead == TRUE) {
    lsnp_row_1 = which(in.dt$LABEL == "LEAD")
    label_data_1 = in.dt[lsnp_row_1, ]
    if (dim(label_data_1)[1] == 0) {
      lsnp_row_1 = in.dt[in.dt$LOG10P == max(in.dt$LOG10P), 
      ]
      label_data_1 = lsnp_row_1[1, ]
    }
  }
  in.dt.2 <- as.data.frame(assoc_data2)
  in.dt.2$POS = as.numeric(as.character(in.dt.2$POS))
  in.dt.2$LOG10P = as.numeric(as.character(in.dt.2$LOG10P))
  in.dt.2$CHR = as.numeric(as.character(in.dt.2$CHR))
  in.dt.2 = dplyr::filter(in.dt.2, .data$CHR == chr_in)
  in.dt.2 = dplyr::filter(in.dt.2, .data$POS > start) %>% dplyr::filter(.data$POS < 
                                                                          end)
  if (label_lead == TRUE) {
    lsnp_row_2 = which(in.dt.2$LABEL == "LEAD")
    label_data_2 = in.dt.2[lsnp_row_2, ]
    if (dim(label_data_2)[1] == 0) {
      lsnp_row_2 = in.dt.2[in.dt.2$LOG10P == max(in.dt.2$LOG10P), 
      ]
      label_data_2 = lsnp_row_2[1, ]
    }
  }
  len1 = nchar(trunc(max(in.dt$LOG10P)))
  len2 = nchar(trunc(max(in.dt.2$LOG10P)))
  scaleFUN0 <- function(x) sprintf("%.0f", x)
  scaleFUN1 <- function(x) sprintf("%.1f", x)
  scaleFUN2 <- function(x) sprintf("%.2f", x)
  scaleFUN3 <- function(x) sprintf("%.3f", x)
  scaleFUN4 <- function(x) sprintf("%.4f", x)
  message("Generating plot.")
  if ("LD" %in% cols_1 && "LD_BIN" %in% cols_1) {
    a = ggplot2::ggplot(data = in.dt, ggplot2::aes_string(x = "POS", 
                                                          y = "LOG10P", color = "LD_BIN")) + ggplot2::geom_point() + 
      ggplot2::scale_colour_manual(values = c(`1.0-0.8` = "#DC0000FF", # "#DB003F", #"red", 
                                              `0.8-0.6` = "#F39B7FFF", # "#F59D10", #"darkorange1", 
                                              `0.6-0.4` = "#00A087FF", # "#49A01D", #"green1", 
                                              `0.4-0.2` = "#4DBBD5FF", # "#1290D9", #"skyblue1", 
                                              `0.2-0.0` = "#3C5488FF", # "#4C81BF", #"navyblue",
                                              `NA` = "#A2A3A4" # "grey"
                                              ), drop = FALSE) + 
      labs(color = bquote(LD~r^2)) + 
      ggplot2::theme_minimal() + 
      ggplot2::xlab(paste0("Chromosome ", chr_in, " Position")) + 
      ggplot2::ylab("-log10(p-value)") + ggplot2::scale_y_reverse() + 
      ggplot2::theme(axis.title.x = ggplot2::element_blank(), 
                     axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
                     panel.border = element_blank(), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank()) + 
      ggplot2::theme(legend.position = "none") + 
      ggplot2::xlim(start, end) + ggplot2::ggtitle(paste0(name1)) + 
      theme(plot.title = element_text(size = 10, vjust = -1)) + theme(plot.margin = margin(5.5, 5.5, -3, 5.5))
  }
  else {
    message("No LD information for dataset #1.")
    a = ggplot2::ggplot(in.dt, ggplot2::aes_string(x = "POS", 
                                                   y = "LOG10P")) + ggplot2::geom_point() + ggplot2::theme_minimal() + 
      ggplot2::xlab(paste0("Chromosome ", chr_in, " Position")) + 
      ggplot2::ylab("-log10(p-value)") + ggplot2::scale_y_reverse() + 
      ggplot2::theme(axis.title.x = ggplot2::element_blank(), 
                     axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
                     panel.border = element_blank(), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank()) + 
      ggplot2::theme(legend.position = "none") + 
      ggplot2::xlim(start, end) + ggplot2::ggtitle(paste0(name1)) + 
      theme(plot.title = element_text(size = 10, vjust = -1)) + 
      theme(plot.margin = margin(5.5, 5.5, -3, 5.5))
  }
  if ("LD" %in% cols_2 && "LD_BIN" %in% cols_2) {
    b = ggplot2::ggplot(data = in.dt.2, ggplot2::aes_string(x = "POS", 
                                                            y = "LOG10P", color = "LD_BIN")) + ggplot2::geom_point() + 
      ggplot2::scale_colour_manual(values = c(`1.0-0.8` = "#DC0000FF", # "#DB003F", #"red", 
                                              `0.8-0.6` = "#F39B7FFF", # "#F59D10", #"darkorange1", 
                                              `0.6-0.4` = "#00A087FF", # "#49A01D", #"green1", 
                                              `0.4-0.2` = "#4DBBD5FF", # "#1290D9", #"skyblue1", 
                                              `0.2-0.0` = "#3C5488FF", # "#4C81BF", #"navyblue",
                                              `NA` = "#A2A3A4" # "grey"
                                              ), drop = FALSE) + 
      labs(color = bquote(LD~r^2)) + 
      ggplot2::theme_minimal() + 
      ggplot2::xlab(paste0("Chromosome ", chr_in, " Position (Mbp)")) + 
      ggplot2::ylab("-log10(p-value)") + ggplot2::theme(legend.position = "bottom") + 
      ggplot2::xlim(start, end) + ggplot2::ylim(min(in.dt.2$LOG10P), 
                                                max(in.dt.2$LOG10P)) + 
      ggplot2::theme(axis.title.x = ggplot2::element_blank(), 
                     axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
                     panel.border = element_blank(), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank()) + 
      ggplot2::ggtitle(paste0(name2)) + theme(plot.title = element_text(size = 10, 
                                                                        vjust = -1))
  }
  else {
    b = ggplot2::ggplot(in.dt.2, ggplot2::aes_string(x = "POS", 
                                                     y = "LOG10P")) + ggplot2::geom_point() + ggplot2::theme_minimal() + 
      ggplot2::xlab(paste0("Chromosome ", chr_in, " Position (Mbp)")) + 
      ggplot2::ylab("-log10(p-value)") + ggplot2::theme(legend.position = "bottom") + 
      ggplot2::xlim(start, end) + ggplot2::ylim(min(in.dt.2$LOG10P), 
                                                max(in.dt.2$LOG10P)) + 
      ggplot2::theme(axis.title.x = ggplot2::element_blank(), 
                     axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
                     panel.border = element_blank(), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank()) + 
      ggplot2::ggtitle(paste0(name2)) + theme(plot.title = element_text(size = 10, 
                                                                        vjust = -1))
  }
  c = ggplot2::ggplot(gene_sub, ggplot2::aes_string(x = "value", 
                                                    y = "y_value")) + 
    ggplot2::geom_line(ggplot2::aes_string(group = "GENE_NAME"), size = 2) + 
    ggplot2::theme_minimal() + 
    ggplot2::geom_text(data = plot_lab, 
                       ggplot2::aes_string(x = "value", y = "y_value", label = "GENE_NAME"), 
                                             hjust = -0.1, vjust = 0.3, size = 2.5) + ggplot2::xlim(start, end) + 
    ggplot2::theme(axis.title.y = ggplot2::element_text(color = "white", size = 28), 
                   axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
                   panel.border = element_blank(), 
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   axis.line.x = element_line(colour = "#595A5C")) + 
    ggplot2::xlab(paste0("Chromosome ", chr_in, " Position")) + 
    ggplot2::ylim(0, (max(gene_sub$y_value) + 1))
  if (len1 == len2) {
    a = a + scale_y_reverse(labels = scaleFUN0) + ggplot2::theme(axis.title.y = ggplot2::element_blank(), #ggplot2::element_text(color = "white", size = 28),  
                                                                 axis.text.y = ggplot2::element_blank(), 
                                                                 axis.ticks.y = ggplot2::element_blank(),
                                                                 panel.border = element_blank(), 
                                                                 panel.grid.major = element_blank(),
                                                                 panel.grid.minor = element_blank(),
                                                                 axis.line.x = element_line(colour = "#595A5C"))
    b = b + scale_y_continuous(labels = scaleFUN0) + ggplot2::theme(axis.title.y = ggplot2::element_blank(), #ggplot2::element_text(color = "white", size = 28),  
                                                                    axis.text.y = ggplot2::element_blank(), 
                                                                    axis.ticks.y = ggplot2::element_blank(),
                                                                    panel.border = element_blank(), 
                                                                    panel.grid.major = element_blank(),
                                                                    panel.grid.minor = element_blank(),
                                                                    axis.line.x = element_line(colour = "#595A5C"))
  }
  else if (len1 > len2) {
    a = a + scale_y_reverse(labels = scaleFUN1)
    diff = len1 - len2
    if (diff == 1) {
      b = b + scale_y_continuous(labels = scaleFUN2)
    }
    else if (diff == 2) {
      b = b + scale_y_continuous(labels = scaleFUN3)
    }
    else if (diff == 3) {
      b = b + scale_y_continuous(labels = scaleFUN4)
    }
  }
  else if (len2 > len1) {
    b = b + scale_y_continuous(labels = scaleFUN1)
    diff = len2 - len1
    if (diff == 1) {
      a = a + scale_y_reverse(labels = scaleFUN2)
    }
    else if (diff == 2) {
      a = a + scale_y_reverse(labels = scaleFUN3)
    }
    else if (diff == 3) {
      a = a + scale_y_reverse(labels = scaleFUN4)
    }
  }
  if (label_lead == TRUE) {
    a = a + geom_point(data = label_data_1, aes_string(x = "POS", 
                                                       y = "LOG10P"), color = "purple")
    a = a + geom_text(data = label_data_1, aes_string(label = "RS_ID"), 
                      color = "black", size = 3, hjust = 1.25)
    b = b + geom_point(data = label_data_2, aes_string(x = "POS", 
                                                       y = "LOG10P"), color = "purple")
    b = b + geom_text(data = label_data_2, aes_string(label = "RS_ID"), 
                      color = "black", size = 3, hjust = 1.25)
  }
  ggpubr::ggarrange(a, b, c, heights = c(2, 2, 1), nrow = 3, 
                    ncol = 1, common.legend = TRUE, legend = "right")
}

# RACER scatterAssociationPlot
scatterPlotRACER2 <- function (assoc_data1, assoc_data2, chr, name1 = "Association Dataset #1", 
                                  name2 = "Association Dataset #2", region_start, region_end, 
                                  ld_df = NULL, label = FALSE) 
{
  reqs = c("CHR", "POS", "LOG10P", "RS_ID")
  cols_1 = colnames(assoc_data1)
  cols_2 = colnames(assoc_data2)
  if (sum(reqs %in% cols_1) == 4) {
  }
  else {
    stop("Association Data Set #1 is missing a required column.")
  }
  if (sum(reqs %in% cols_2) == 4) {
  }
  else {
    stop("Association Data Set #2 is missing a required column.")
  }
  `%>%` <- magrittr::`%>%`
  message("Reading in association data")
  in.dt <- as.data.frame(assoc_data1)
  in.dt$POS = as.numeric(as.character(in.dt$POS))
  in.dt$LOG10P = as.numeric(as.character(in.dt$LOG10P))
  in.dt$CHR = as.numeric(as.character(in.dt$CHR))
  in.dt = dplyr::filter(in.dt, .data$CHR == chr)
  in.dt = dplyr::filter(in.dt, .data$POS > region_start) %>% 
    dplyr::filter(.data$POS < region_end)
  in.dt.2 <- as.data.frame(assoc_data2)
  in.dt.2$POS = as.numeric(as.character(in.dt.2$POS))
  in.dt.2$LOG10P = as.numeric(as.character(in.dt.2$LOG10P))
  in.dt.2$CHR = as.numeric(as.character(in.dt.2$CHR))
  in.dt.2 = dplyr::filter(in.dt.2, .data$CHR == chr)
  in.dt.2 = dplyr::filter(in.dt.2, .data$POS > region_start) %>% 
    dplyr::filter(.data$POS < region_end)
  if (ld_df > 0) {
    if (ld_df == 1) {
      in.dt.final = dplyr::select(in.dt, "RS_ID", "LOG10P", 
                                  "LD", "LD_BIN")
      colnames(in.dt.final) = c("RS_ID", "LOG10P1", "LD", 
                                "LD_BIN")
      in.dt.2.final = dplyr::select(in.dt.2, "RS_ID", "LOG10P")
      colnames(in.dt.2.final) = c("RS_ID", "LOG10P2")
    }
    else if (ld_df == 2) {
      in.dt.2.final = dplyr::select(in.dt.2, "RS_ID", "LOG10P", 
                                    "LD", "LD_BIN")
      colnames(in.dt.2.final) = c("RS_ID", "LOG10P2", "LD", 
                                  "LD_BIN")
      in.dt.final = dplyr::select(in.dt, "RS_ID", "LOG10P")
      colnames(in.dt.final) = c("RS_ID", "LOG10P1")
    }
  }
  else {
    in.dt.final = dplyr::select(in.dt, "RS_ID", "LOG10P")
    colnames(in.dt.final) = c("RS_ID", "LOG10P1")
    in.dt.2.final = dplyr::select(in.dt.2, "RS_ID", "LOG10P")
    colnames(in.dt.2.final) = c("RS_ID", "LOG10P2")
  }
  df_plot = merge(in.dt.final, in.dt.2.final, by = "RS_ID")
  lab.in = df_plot[which.max(df_plot$LOG10P1 + df_plot$LOG10P2), 
  ]
  message(paste0("Generating plot for ", name1, " vs. ", name2, "."))
  if (ld_df > 0) {
    ggplot2::ggplot(df_plot, aes_string(x = "LOG10P1", y = "LOG10P2", 
                                        color = "LD_BIN")) + ggplot2::geom_point(size = 4) + 
      ggplot2::xlab(bquote("-"~log[10]~"("~italic(p)~"-value) for "~.(name1))) + 
      ggplot2::ylab(bquote("-"~log[10]~"("~italic(p)~"-value) for "~.(name2))) + 
      ggplot2::theme_bw() + ggplot2::scale_colour_manual(values = c(`1.0-0.8` = "#E55738", 
                                                                    `0.8-0.6` = "#F59D10", `0.6-0.4` = "#49A01D", 
                                                                    `0.4-0.2` = "#1396D8", `0.2-0.0` = "#4C81BF", 
                                                                    `NA` = "#A2A3A4"), 
                                                         drop = FALSE) + ggplot2::geom_point(data = lab.in, 
                                                                                             color = "#9A3480", size = 4) + geom_text(data = lab.in, aes_string(label = "RS_ID"), 
                                                                                                                                      color = "black", size = 5, hjust = 1.25) +
      theme(
        axis.text.x = element_text(size = 14, face = "plain"),
        axis.text.y = element_text(size = 14, face = "plain"),  
        axis.title.x = element_text(size = 16, face = "plain"),
        axis.title.y = element_text(size = 16, face = "plain"))
  }
  else {
    ggplot2::ggplot(df_plot, aes_string(x = "-LOG10P1", y = "-LOG10P2")) + 
      ggplot2::geom_point(size = 4) + 
      ggplot2::xlab(bquote("-"~log[10]~"("~italic(p)~"-value) for "~.(name1))) + 
      ggplot2::ylab(bquote("-"~log[10]~"("~italic(p)~"-value) for "~.(name2))) + 
      ggplot2::theme_bw() + 
      ggplot2::geom_point(data = lab.in, 
                          color = "#9A3480", size = 4) + geom_text(data = lab.in, aes_string(label = "RS_ID"), 
                                                                   color = "black", size = 5, hjust = 1.25) +
      theme(
        axis.text.x = element_text(size = 14, face = "plain"),
        axis.text.y = element_text(size = 14, face = "plain"),  
        axis.title.x = element_text(size = 16, face = "plain"),
        axis.title.y = element_text(size = 16, face = "plain"))
  }
}

#' Create a regional association mirror plot with ggplot2.
mirrorPlotRACER_NEW <- function (assoc_data1, assoc_data2, chr, build = "hg19", set = "protein_coding", 
          name1 = "Association Dataset #1", name2 = "Association Dataset #2", 
          plotby, gene_plot = NULL, snp_plot = NULL, start_plot = NULL, 
          end_plot = NULL, label_lead = FALSE) 
  {
    reqs = c("CHR", "POS", "LOG10P")
    cols_1 = colnames(assoc_data1)
    cols_2 = colnames(assoc_data2)
    if (sum(reqs %in% cols_1) == 3) {
    }
    else {
      stop("Association Data Set #1 is missing a required column.")
    }
    if (sum(reqs %in% cols_2) == 3) {
    }
    else {
      stop("Association Data Set #2 is missing a required column.")
    }
    if (build == "hg38") {
      utils::data(hg38)
      chr_in = chr
      colnames(hg38) = c("GENE_ID", "CHR", "TRX_START", "TRX_END", 
                         "LENGTH", "GENE_NAME", "TYPE")
      gene_sub = hg38[hg38$CHR == chr_in, ]
    }
    else if (build == "hg19") {
      utils::data(hg19)
      chr_in = chr
      colnames(hg19) = c("GENE_ID", "CHR", "TRX_START", "TRX_END", 
                         "LENGTH", "GENE_NAME", "TYPE")
      gene_sub = hg19[hg19$CHR == chr_in, ]
    }
    if (set == "protein_coding") {
      gene_sub = gene_sub[gene_sub$TYPE == "protein_coding", 
      ]
    }
    else {
      gene_sub = gene_sub
    }
    `%>%` <- magrittr::`%>%`
    if ((sum(is.null(plotby)) == 0) == TRUE) {
      message("Plotting by...")
      if ((plotby == "coord") == TRUE) {
        message("coord")
        start = start_plot
        end = end_plot
      }
      else if ((plotby == "gene") == TRUE) {
        message(paste("gene:", gene_plot))
        if (sum(is.null(gene_plot)) == 0) {
          p = subset(gene_sub, gene_sub$GENE_NAME == gene_plot)
          start = min(p$TRX_START) - 5e+05
          end = max(p$TRX_END) + 5e+05
        }
        else {
          stop("No gene specified.")
        }
      }
      else if ((plotby == "snp") == TRUE) {
        message(paste("snp", snp_plot))
        q = assoc_data1[assoc_data1$RS_ID == snp_plot, ]
        w = q$POS
        w = as.numeric(as.character(w))
        start = w - 5e+05
        end = w + 5e+05
      }
    }
    else {
      stop("Please specify a parameter to plotby.")
    }
    gene_sub = subset(gene_sub, gene_sub$TRX_START > (start - 
                                                        5000))
    gene_sub = subset(gene_sub, gene_sub$TRX_END < (end + 5000))
    gene_sub = gene_sub[!duplicated(gene_sub$GENE_ID), ]
    gene_sub = gene_sub[, c(3, 4, 6)]
    gene_sub = reshape2::melt(gene_sub, id.vars = "GENE_NAME")
    gene_sub$y_value = as.numeric(as.factor(gene_sub$GENE_NAME))
    plot_lab = subset(gene_sub, gene_sub$variable == "TRX_END")
    
    message("Reading in association data.")
    message("- dataset #1")
    in.dt <- as.data.frame(assoc_data1)
    in.dt$POS = as.numeric(in.dt$POS)
    in.dt$LOG10P = as.numeric(in.dt$LOG10P)
    in.dt$CHR = as.numeric(as.character(in.dt$CHR))
    in.dt = dplyr::filter(in.dt, .data$CHR == chr_in)
    in.dt = dplyr::filter(in.dt, .data$POS > start) %>% dplyr::filter(.data$POS < 
                                                                        end)
    # bug fixing
    # print(head(in.dt))
    if (label_lead == TRUE) {
      lsnp_row_1 = which(in.dt$LABEL == "LEAD")
      label_data_1 = in.dt[lsnp_row_1, ]
      if (dim(label_data_1)[1] == 0) {
        lsnp_row_1 = in.dt[in.dt$LOG10P == max(in.dt$LOG10P), 
        ]
        label_data_1 = lsnp_row_1[1, ]
      }
    }
    message("- dataset #2")
    in.dt.2 <- as.data.frame(assoc_data2)
    in.dt.2$POS = as.numeric(in.dt.2$POS)
    in.dt.2$LOG10P = as.numeric(in.dt.2$LOG10P)
    in.dt.2$CHR = as.numeric(as.character(in.dt.2$CHR))
    in.dt.2 = dplyr::filter(in.dt.2, .data$CHR == chr_in)
    in.dt.2 = dplyr::filter(in.dt.2, .data$POS > start) %>% dplyr::filter(.data$POS < 
                                                                            end)
    # bug fixing
    # print(head(in.dt.2))
    if (label_lead == TRUE) {
      lsnp_row_2 = which(in.dt.2$LABEL == "LEAD")
      label_data_2 = in.dt.2[lsnp_row_2, ]
      if (dim(label_data_2)[1] == 0) {
        lsnp_row_2 = in.dt.2[in.dt.2$LOG10P == max(in.dt.2$LOG10P), 
        ]
        label_data_2 = lsnp_row_2[1, ]
      }
    }
    
    y.max <- trunc(max(in.dt.2$LOG10P, na.rm = TRUE))
    y.max.2 <- trunc(max(in.dt$LOG10P, na.rm = TRUE))
    
    if (y.max > y.max.2) {
      y.max.lim = y.max + 1
    } else {
      y.max.lim = y.max.2 + 1
    }
    rm(y.max, y.max.2)
    
    len1 = nchar(trunc(max(in.dt$LOG10P, na.rm = TRUE)))
    # print(len1)
    len2 = nchar(trunc(max(in.dt.2$LOG10P, na.rm = TRUE)))
    # print(len2)
    scaleFUN0 <- function(x) sprintf("%.0f", x)
    scaleFUN1 <- function(x) sprintf("%.1f", x)
    scaleFUN2 <- function(x) sprintf("%.2f", x)
    scaleFUN3 <- function(x) sprintf("%.3f", x)
    scaleFUN4 <- function(x) sprintf("%.4f", x)
    message("Generating plot.")
    if ("LD" %in% cols_1 && "LD_BIN" %in% cols_1) {
      a = ggplot2::ggplot(data = in.dt, ggplot2::aes_string(x = "POS", 
                                                            y = "LOG10P", color = "LD_BIN")) + ggplot2::geom_point() + 
        ggplot2::scale_colour_manual(values = c(`1.0-0.8` = "#E55738", 
                                                `0.8-0.6` = "#F59D10", `0.6-0.4` = "#49A01D", 
                                                `0.4-0.2` = "#1396D8", `0.2-0.0` = "#4C81BF", 
                                                `NA` = "#A2A3A4"), drop = FALSE) + ggplot2::theme_bw() + 
        ggplot2::xlab(paste0("Chromosome ", chr_in, " Position")) + 
        ggplot2::ylab("-log10(p-value)") + ggplot2::scale_y_reverse() + 
        ggplot2::theme(axis.title.x = ggplot2::element_blank(), 
                       axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank()) + 
        ggplot2::theme(legend.position = "none") + ggplot2::xlim(start, 
                                                                 end) + ggplot2::ggtitle(paste0(name1)) + theme(plot.title = element_text(size = 10, 
                                                                                                                                          vjust = -1)) + theme(plot.margin = margin(5.5, 5.5, 
                                                                                                                                                                                    -3, 5.5))
      }
    else {
      message("No LD information for dataset #1.")
      a = ggplot2::ggplot(in.dt, ggplot2::aes_string(x = "POS", 
                                                     y = "LOG10P")) + ggplot2::geom_point() + ggplot2::theme_bw() + 
        ggplot2::xlab(paste0("Chromosome ", chr_in, " Position")) + 
        ggplot2::ylab("-log10(p-value)") + ggplot2::scale_y_reverse() + 
        ggplot2::theme(axis.title.x = ggplot2::element_blank(), 
                       axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank()) + 
        ggplot2::theme(legend.position = "none") + ggplot2::xlim(start, 
                                                                 end) + ggplot2::ggtitle(paste0(name1)) + theme(plot.title = element_text(size = 10, 
                                                                                                                                          vjust = -1)) + theme(plot.margin = margin(5.5, 5.5, 
                                                                                                                                                                                    -3, 5.5))
    }
    if ("LD" %in% cols_2 && "LD_BIN" %in% cols_2) {
      b = ggplot2::ggplot(data = in.dt.2, ggplot2::aes_string(x = "POS", 
                                                              y = "LOG10P", color = "LD_BIN")) + ggplot2::geom_point() + 
        ggplot2::scale_colour_manual(values = c(`1.0-0.8` = "#E55738", 
                                                `0.8-0.6` = "#F59D10", `0.6-0.4` = "#49A01D", 
                                                `0.4-0.2` = "#1396D8", `0.2-0.0` = "#4C81BF", 
                                                `NA` = "#A2A3A4"), drop = FALSE) + ggplot2::theme_bw() + 
        ggplot2::xlab(paste0("Chromosome ", chr_in, " Position (Mbp)")) + 
        ggplot2::ylab("-log10(p-value)") + ggplot2::theme(legend.position = "bottom") + 
        ggplot2::xlim(start, end) + ggplot2::ylim(min(in.dt.2$LOG10P), 
                                                  max(in.dt.2$LOG10P)) + ggplot2::theme(axis.title.x = ggplot2::element_blank(), 
                                                                                        axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank()) + 
        ggplot2::ggtitle(paste0(name2)) + theme(plot.title = element_text(size = 10, 
                                                                          vjust = -1))
    }
    else {
      b = ggplot2::ggplot(in.dt.2, ggplot2::aes_string(x = "POS", 
                                                       y = "LOG10P")) + ggplot2::geom_point() + ggplot2::theme_bw() + 
        ggplot2::xlab(paste0("Chromosome ", chr_in, " Position (Mbp)")) + 
        ggplot2::ylab("-log10(p-value)") + ggplot2::theme(legend.position = "bottom") + 
        ggplot2::xlim(start, end) + ggplot2::ylim(min(in.dt.2$LOG10P), 
                                                  max(in.dt.2$LOG10P)) + ggplot2::theme(axis.title.x = ggplot2::element_blank(), 
                                                                                        axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank()) + 
        ggplot2::ggtitle(paste0(name2)) + theme(plot.title = element_text(size = 10, 
                                                                          vjust = -1))
    }
    c = ggplot2::ggplot(gene_sub, ggplot2::aes_string(x = "value", 
                                                      y = "y_value")) + ggplot2::geom_line(ggplot2::aes_string(group = "GENE_NAME"), 
                                                                                           size = 2) + ggplot2::theme_bw() + ggplot2::geom_text(data = plot_lab, 
                                                                                                                                                ggplot2::aes_string(x = "value", y = "y_value", label = "GENE_NAME"), 
                                                                                                                                                hjust = -0.1, vjust = 0.3, size = 2.5) + ggplot2::xlim(start, 
                                                                                                                                                                                                       end) + ggplot2::theme(axis.title.y = ggplot2::element_text(color = "white", 
                                                                                                                                                                                                                                                                  size = 28), axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank()) + 
      ggplot2::xlab(paste0("Chromosome ", chr_in, " Position")) + 
      ggplot2::ylim(0, (max(gene_sub$y_value) + 1))
    
    if (len1 == len2) {
      a = a + scale_y_reverse(labels = scaleFUN0)
      b = b + scale_y_continuous(labels = scaleFUN0)
    }
    else if (len1 > len2) {
      a = a + scale_y_reverse(labels = scaleFUN1)
      diff = len1 - len2
      if (diff == 1) {
        b = b + scale_y_continuous(labels = scaleFUN2)
      }
      else if (diff == 2) {
        b = b + scale_y_continuous(labels = scaleFUN3)
      }
      else if (diff == 3) {
        b = b + scale_y_continuous(labels = scaleFUN4)
      }
    }
    else if (len2 > len1) {
      b = b + scale_y_continuous(labels = scaleFUN1)
      diff = len2 - len1
      if (diff == 1) {
        a = a + scale_y_reverse(labels = scaleFUN2)
      }
      else if (diff == 2) {
        a = a + scale_y_reverse(labels = scaleFUN3)
      }
      else if (diff == 3) {
        a = a + scale_y_reverse(labels = scaleFUN4)
      }
    }
    if (label_lead == TRUE) {
      a = a + geom_point(data = label_data_1, aes_string(x = "POS", 
                                                         y = "LOG10P"), color = "purple")
      a = a + geom_text(data = label_data_1, aes_string(label = "RS_ID"), 
                        color = "black", size = 3, hjust = 1.25)
      b = b + geom_point(data = label_data_2, aes_string(x = "POS", 
                                                         y = "LOG10P"), color = "purple")
      b = b + geom_text(data = label_data_2, aes_string(label = "RS_ID"), 
                        color = "black", size = 3, hjust = 1.25)
    }
    ggpubr::ggarrange(a, b, c, heights = c(2, 2, 1), nrow = 3, 
                      ncol = 1, common.legend = TRUE, legend = "right")
}


################################################################################
# MANHATTAN PLOTTING
################################################################################
# based on the qqman package by Stephen Turner


manhattan_uu <- function (x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", col = c("gray10", "gray60"), chrlabs = NULL, suggestiveline = -log10(1e-05), 
                          genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE, 
                          annotatePval = NULL, annotateTop = TRUE, ...) {
  CHR = BP = P = index = NULL
  if (!(chr %in% names(x))) 
    stop(paste("Column", chr, "not found!"))
  if (!(bp %in% names(x))) 
    stop(paste("Column", bp, "not found!"))
  if (!(p %in% names(x))) 
    stop(paste("Column", p, "not found!"))
  if (!(snp %in% names(x))) 
    warning(paste("No SNP column found. OK unless you're trying to highlight."))
  if (!is.numeric(x[[chr]])) 
    stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
  if (!is.numeric(x[[bp]])) 
    stop(paste(bp, "column should be numeric."))
  if (!is.numeric(x[[p]])) 
    stop(paste(p, "column should be numeric."))
  if (!is.null(x[[snp]])) 
    d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]], 
                   pos = NA, index = NA, SNP = x[[snp]], stringsAsFactors = FALSE)
  else d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]], 
                      pos = NA, index = NA)
  d <- d[order(d$CHR, d$BP), ]
  if (logp) {
    d$logp <- -log10(d$P)
  }
  else {
    d$logp <- d$P
  }
  d$index = rep.int(seq_along(unique(d$CHR)), times = tapply(d$SNP, 
                                                             d$CHR, length))
  nchr = length(unique(d$CHR))
  if (nchr == 1) {
    d$pos = d$BP
    xlabel = paste("Chromosome", unique(d$CHR), "position")
  }
  else {
    lastbase = 0
    ticks = NULL
    for (i in unique(d$index)) {
      if (i == 1) {
        d[d$index == i, ]$pos = d[d$index == i, ]$BP
      }
      else {
        lastbase = lastbase + max(d[d$index == (i - 1), 
                                    "BP"])
        d[d$index == i, "BP"] = d[d$index == i, "BP"] - 
          min(d[d$index == i, "BP"]) + 1
        d[d$index == i, "pos"] = d[d$index == i, "BP"] + 
          lastbase
      }
    }
    ticks <- tapply(d$pos, d$index, quantile, probs = 0.5)
    xlabel = "Chromosome"
    labs <- unique(d$CHR)
  }
  xmax = ceiling(max(d$pos) * 1.03)
  xmin = floor(max(d$pos) * -0.03)
  def_args <- list(xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", 
                   las = 1, pch = 20, xlim = c(xmin, xmax), ylim = c(0, 
                                                                     ceiling(max(d$logp))), xlab = xlabel, ylab = expression(-log[10](italic(p))))
  dotargs <- list(...)
  do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% 
                                            names(dotargs)]))
  if (!is.null(chrlabs)) {
    if (is.character(chrlabs)) {
      if (length(chrlabs) == length(labs)) {
        labs <- chrlabs
      }
      else {
        warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
      }
    }
    else {
      warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
    }
  }
  if (nchr == 1) {
    axis(1, ...)
  }
  else {
    axis(1, at = ticks, labels = labs, ...)
  }
  col = rep_len(col, max(d$index))
  if (nchr == 1) {
    with(d, points(pos, logp, pch = 20, col = col[1], ...))
  }
  else {
    icol = 1
    for (i in unique(d$index)) {
      points(d[d$index == i, "pos"], d[d$index == i, "logp"], 
             col = col[icol], pch = 20, ...)
      icol = icol + 1
    }
  }
  if (suggestiveline) 
    # abline(h = suggestiveline, col = "blue") # original
    abline(h = suggestiveline, col = "#595A5C", lty = "dashed") # original
  if (genomewideline) 
    # abline(h = genomewideline, col = "red") # original
    abline(h = genomewideline, col = "#E55738", lty = "dashed") # original
  if (!is.null(highlight)) {
    if (any(!(highlight %in% d$SNP))) 
      warning("You're trying to highlight SNPs that don't exist in your results.")
    d.highlight = d[which(d$SNP %in% highlight), ]
    # with(d.highlight, points(pos, logp, col = "green3", pch = 20, ...)) # original
    with(d.highlight, points(pos, logp, col = "#9FC228", pch = 20, ...))
  }
  if (!is.null(annotatePval)) {
    if (logp) {
      topHits = subset(d, P <= annotatePval)
    }
    else topHits = subset(d, P >= annotatePval)
    par(xpd = TRUE)
    if (annotateTop == FALSE) {
      if (logp) {
        with(subset(d, P <= annotatePval), textxy(pos, 
                                                  -log10(P), offset = 0.625, labs = topHits$SNP, 
                                                  cex = 0.45), ...)
      }
      else with(subset(d, P >= annotatePval), textxy(pos, 
                                                     P, offset = 0.625, labs = topHits$SNP, cex = 0.45), 
                ...)
    }
    else {
      topHits <- topHits[order(topHits$P), ]
      topSNPs <- NULL
      for (i in unique(topHits$CHR)) {
        chrSNPs <- topHits[topHits$CHR == i, ]
        topSNPs <- rbind(topSNPs, chrSNPs[1, ])
      }
      if (logp) {
        textxy(topSNPs$pos, -log10(topSNPs$P), offset = 0.625, 
               labs = topSNPs$SNP, cex = 0.5, ...)
      }
      else textxy(topSNPs$pos, topSNPs$P, offset = 0.625, 
                  labs = topSNPs$SNP, cex = 0.5, ...)
    }
  }
  par(xpd = FALSE)
}

