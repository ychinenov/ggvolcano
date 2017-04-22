ggvolcano <- function(gene.id = NULL,
                      p.val = NULL,
                      fold.c = NULL,
                      labels = "none",
                      de.l = list(),
                      top.list = NULL,
                      cutoff.p = 0.01,
                      cutoff.fc = 2,
                      xlim = vector(),
                      ylim = vector(),
                      logfold.c = T,
                      FDR = T,
                      binhex = F) {
  
  
  
  #_____________________________________________________________________
  # This function producess pretty volcano plot using ggplot framework
  # the following arguments are currently available
  #gene.id=NULL, #vector of gene names
  #p.val=NULL,   #vector of pvalues
  #fold.c=NULL,  # fector of fold changes
  #labels=NULL,   # vector of gene names to be show on the plot. If NULL upper 15 sorted by p.val will be labeled
  #                 labels="none" - no labels, labels="auto" top 20 genes sorted by p value (or adjusted pvalue)
  #de.l=substitute(list()), #a list of gene lists for which dots will be colored; if NULL all genes will be colored according to cutoffs
  #top.list=NULL   # when lists of gene to color is given labels top N genes, colors by group
  #cutoff.p=0.01,
  #cutoff.fc=2,
  #xlim=vector()  #by default will use min(log2(fold.c),max(lof2(fold.c))), otherwise expects a vector c(min,max)
  #xlim=vector()  #by default
  #logfold.c=T,   #by default assumes that fold change is log2 transformed, otherwise fold.c will be log2 transformed
  #FDR=T,         #by default assumes that p values are FDR, otherwise p.addjust() will be applied with method="fdr"
  #binhex=F){     # if binhex=T geom_hex will be used to bin nearby data points to avoid overplotting
  
  
  #______Checking input data
  
  if (is.null(gene.id)) {
    warning("***** No gene IDs were provided.*****")
    no.id = TRUE
  } else{
    no.id <-FALSE
    gene.id <- as.character(gene.id)
  } #we can leave without geneIDs but numbers are ugly
  
  if (length(p.val) != length(fold.c)) {
    warning("*****The number of provided p values and fold changes are not equal.Full Stop!****")
    stop()
  }
  
  if (length(gene.id) != length(p.val)) {
    warning("*****The number of p values/fold changes and gene IDs are not equal.No gene ID will be used !****")
    no.id = TRUE
  } else{
    no.id <- FALSE
  }
  
  if(is.null(labels)){labels="none"}
  
  # check and instal required packages
  
  if (!"ggplot2" %in% installed.packages()) {
    install <- readline(prompt = "Package 'ggplot2' is required but not found. Install? (y/n).")
    if (install == "y") {
       install.packages(c("ggplot2"), dep = TRUE)
       require("ggrepel")
    } else {
       stop("Too bad!")
    }
  } else {
    require("ggplot2")
  }
  
  if (!"ggrepel" %in% installed.packages()) {
    install <- readline(prompt = "Package 'ggrepel' is required but not found. Install? (y/n).")
    if (install == "y") {
      install.packages(c("ggrepel"), dep = TRUE)
      require("ggrepel")
    } else {
      stop("Too bad!")
    }
  } else {
    require("ggrepel")
  }
  
  #My prefered theme for ggplots
  t2 <- theme(
    plot.margin = unit(c(2, 5, 2, 2), "lines"),
    axis.title.x = element_text( face = "bold", color = "black", size = 14), 
    axis.title.y = element_text( face = "bold", color = "black", size = 14),
    axis.text = element_text(face = "bold", color = "black", size = 14),
    plot.title = element_text(face = "bold", color = "black", size = 14),
    legend.title = element_text(face = "bold", color = "black", size = 14),
    legend.text = element_text(face = "bold", colour = "black", size = 14),
    legend.key = element_rect(colour = 'white', fill = 'white', size = 0.5),
    legend.direction = 'horizontal',
    legend.position = "top",
    legend.spacing.x = unit(c(2), "mm"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA),
    panel.background = element_blank(),
    axis.line = element_line(size = 0.4)
  )
  
  #Transforming vectors of fold changes and p values is needed
  
  if (FDR) {
    p.val <- -log10(p.val)
  } else
  {
    p.val <- -log10(p.adjust(p.val, method = "fdr"))
  }
  
  if (logfold.c == FALSE) {
    fold.c <- log2(fold.c)
  }
  
  #Expressions for X and Y axis labels
  ylegend = expression(-log[10] ~ p.val)
  xlegend <- expression(log[2] ~ fold ~ change)
  
  #Building basic framework: combine gene.id, fold.c and p.val into a single dataframe
  df <- cbind.data.frame(gene.id, fold.c, p.val, stringsAsFactors = F)
  df <- df[order(df$p.val, decreasing = T), ]
  
  #Create a table of labels depending on user's input
  if (labels == "none") {
    labels.t <-data.frame(
                 gene.id = character(),
                 fold.c = double(),
                 p.val = double(),
                 list_name = character())
  } else if (labels == "auto") {
    labels.t <- df[1:20, ]
    labels.t$list_name <- "User"
  }   else {
    labels <- unique(labels)
    labels.t <- df[df$gene.id %in% labels,]
       if(nrow(labels.t)==0){warning("No matching genes found! Labels will be set to labels=\"auto\"")
           labels.t <- df[1:20, ]
           labels.t$list_name <- "User"} else {
    labels.t$list_name <- "User"}
  }
  
  # determining max Y and quanttile 75% for plot scaling
  maxY <- max(p.val)
  quant4 <- quantile(p.val)[[4]]
 
  #seting up xlim/ylim if not provided
  if(length(xlim) !=2){xlim <- c(min(fold.c),max(fold.c))}
  if(length(ylim) !=2){ylim <- c(-0.4,max(p.val)) } else if(ylim[1]<0){warning("-log10(p.val) cannot be negative. Reset to 0:max(p.val)"); ylim <-c(-0.4,max(p.val))} 
    
    
   #Plotting basic volcano plot here when no subsets to visualize were provided
  if (length(eval(de.l)) == 0 | no.id == T) {
    p <- ggplot() + xlim(xlim) + ylim(ylim) + t2 +
    {
      if (binhex == F)
      {
        geom_point(mapping = aes( x = fold.c, y = p.val, col = (abs(fold.c) >= log2(cutoff.fc) & p.val > -log10(cutoff.p))), alpha = 0.5, size = 6 )
      } else
      {
        geom_hex(mapping = aes(x = fold.c, y = p.val), col = "gray85", fill = "gray80",alpha = 0.5, size = 6, bins = (length(gene.id) / 250))
      }
    } +
      geom_hline(yintercept = -log10(cutoff.p), col = "purple", linetype = "dashed", size = 1 ) +
      geom_vline(xintercept = c(-log2(cutoff.fc), log2(cutoff.fc)), col = "blue",linetype = "dashed", size = 1) +
      geom_text(aes(label = paste0("p=", cutoff.p), x = max(fold.c) - 0.15, y = -log10(cutoff.p) + 0.15),  size = 5) +
      geom_text(aes(label = paste0(c( paste("-", cutoff.fc, "fold"),paste0("+", cutoff.fc, "fold"))), x = c(-log2(cutoff.fc) - 1.2, log2(cutoff.fc) + 1.2), y = (-0.4), size = 5)) +
        {
        if (binhex == T)
        {
          geom_point(data = df, mapping = aes(x = fold.c, y = p.val), col = ifelse((abs(fold.c) >= log2(cutoff.fc) &  p.val > -log10(cutoff.p)), "red", NA), alpha = 0.5, size = 6)
        } else
        {
          scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "gray"))
        }
      } +
      xlab(label = expression(log[2] ~ fold ~ change)) +
      ylab(label = expression(-log[10] ~ p ~ value)) +
      geom_point(data = labels.t, mapping = aes(x = fold.c, y = p.val), size = 7, shape = 1, stroke = 1.2) +
      geom_label_repel(mapping = aes(x = fold.c, y = p.val, label = gene.id), 
                       data = labels.t, 
                       force =3,
                       nudge_x = ifelse(sign(labels.t$fold.c) == -1, -5, 5), 
                       nudge_y = ifelse(labels.t$p.val > quantile(p.val)[[4]], 2, -2),
                       point.padding = unit(0.7, 'lines')) +
      theme(legend.position = "none")
    
  } else {
    #extract subset names from user input
    nm <- deparse(substitute(de.l))
    nm <- gsub("[\\(\\)]", "", regmatches(nm, gregexpr("\\(.*?\\)", nm)))
    nm <- unlist(strsplit(nm, ", "))
    de.l <- eval(de.l) # evaluate the expression
    names(de.l) <- nm # and put names from a function call into actual list
    
    
    # creating a long form dataframe from input list
    a <- list()# This list contains dataframes. Each dataframe contains genes from input input list along with pvalues/fold changes
    for (i in 1:length(de.l)) {
      #this really needs cleaning.... eventually
      a[[i]] <- df[df$gene.id %in% de.l[[i]],]
      a[[i]]$list_name <- names(de.l[i])
    }
    a <-  do.call(rbind.data.frame, a)
    a$list_name <- as.factor(a$list_name)
    
    # If user wants to label top (by pvalue) genes of each list these genes need to be added to the dataframe of labels. If lebels="none"
    # a new table will be created that contains top genes from each list
    if (!is.null(top.list)) {
      temp <- a[ave(a$p.val, a$list_name, FUN = seq_along) <= top.list, ] #may require presorting
      if (labels == "none") { 
        labels.t <- temp
      } else{
        labels.t <- rbind.data.frame(temp, labels.t)
      }
    }
    p <- ggplot() + xlim(xlim) + ylim(ylim) + t2 +
    {
      if (binhex == F)
      {
        geom_point(mapping = aes(x = fold.c, y = p.val), alpha = 0.5, size = 6, col = "gray80")
      } else
      {
        geom_hex(mapping = aes(x = fold.c, y = p.val), col = "gray85", fill = "gray80", alpha = 0.5, size = 6, bins = (length(gene.id)/250))
      }
    } + #Is there a "smart" way to calculate bin number?
      
      geom_hline(yintercept = -log10(cutoff.p), col = "purple", linetype = "dashed", size = 1) +
      geom_vline(xintercept = c(-log2(cutoff.fc), log2(cutoff.fc)), col = "blue", linetype = "dashed", size = 1 ) +
      geom_text(aes(label = paste0("p=", cutoff.p), x = max(fold.c) - 0.15, y = -log10(cutoff.p) + 0.15),  size = 5) +
      geom_text(aes(label = paste0(c(paste("-", cutoff.fc, "fold"), paste0("+", cutoff.fc, "fold"))), 
                    x = c(-log2(cutoff.fc) - 1.2, log2(cutoff.fc) + 1.2) ,y = (-0.4),size = 10)) +
      xlab(label = expression(log[2] ~ fold ~ change)) +
      ylab(label = expression(-log[10] ~ p ~ value)) +
      geom_point(data = a,mapping = aes(x = fold.c, y = p.val, col = list_name),alpha = 0.5, size = 6) +
      {
        if (binhex == F) # this conditional needs to be collapsed to a single geom. I don't think I need two different appearances.
        {
          geom_point(data = labels.t, mapping = aes(x = fold.c, y = p.val, col = list_name), size = 7, shape = 1, stroke = 1.2)
        } else
        {
          geom_point(data = labels.t, mapping = aes(x = fold.c, y = p.val, col = list_name), size = 7, shape = 1, stroke = 1.2)
        }
      } +
      geom_label_repel(mapping = aes(x = fold.c, y = p.val, label = gene.id, col = factor(list_name)),
                      data = labels.t,
                      nudge_x = ifelse(sign(labels.t$fold.c) == -1, -5, 5), 
                      nudge_y = ifelse(labels.t$p.val > quantile(p.val)[[4]], 2, -2), 
                      force = 1,
                      point.padding = unit(0.7, 'lines'), segment.size = 0.2, show.legend = FALSE) +
      
      {
        if ("ggsci" %in% installed.packages())  #conditional color choice
        {
          scale_color_npg()
        } else if ("RColorBrewer" %in% installed.packages()) {
          scale_color_brewer(palette = "Dark2")
        } else
        {
          scale_color_hue(l = 50)
        }
      } +
      guides(size = F) + #removing legend label for hexagonal bins
      guides(col = guide_legend(title = "Gene set:", keywidth = 2))
  }
  
  return(p)
}

#p <- ggvolcano(gene.id = toptag[[2]]$genes,
#                      p.val = toptag[[2]]$`FDR_-1*WT 1*Fcr2bKO`,
#                      fold.c = toptag[[2]]$`logFC_-1*WT 1*Fcr2bKO`,
#                      labels = c("Klf2","Klf7","Klf5","Klf9"),
#                      de.l = list(EMT_hallmark,TNF_halmark,INFLAMMATION_halmark),
#                      top.list = 5,
#                      cutoff.p = 0.01,
#                      cutoff.fc = 2,
#                      xlim = c(),
#                      ylim = c(),
#                      logfold.c = T,
#                      FDR = T,
#                      binhex = F) 

#rm(gene.id,p.val,fold.c,labels,de.l,top.list,cutoff.p,cutoff.fc,xlim,ylim,logfold.c,FDR,binhex)

