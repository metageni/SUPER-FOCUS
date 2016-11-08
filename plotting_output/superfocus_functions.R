#!/usr/bin/Rscript
# superfocus_functions.R
# Figures for functions from SUPER-FOCUS output.
# GitHub: https://github.com/metageni/SUPER-FOCUS/
# Ref: Silva GGZ, Green K., B. E. Dutilh, and R. A. Edwards: SUPER-FOCUS: A tool for agile functional analysis of shotgun metagenomic data.
#      (Bioinformatics. 2015 Oct 9. pii: btv584.)
#
# Author: Daniel A Cuevas (dcuevas08.at.gmail.com)
# Created on 26 Oct 2016
# Updated on 08 Nov 2016

# Import necessary packages
# These may need to be installed first
if ("getopt" %in% rownames(installed.packages()) == F) {
    install.packages("getopt")
}
if ("ggplot2" %in% rownames(installed.packages()) == F) {
    install.packages("ggplot2")
}
if ("reshape2" %in% rownames(installed.packages()) == F) {
    install.packages("reshape2")
}
if ("plyr" %in% rownames(installed.packages()) == F) {
    install.packages("plyr")
}
if ("gridExtra" %in% rownames(installed.packages()) == F) {
    install.packages("gridExtra")
}
suppressMessages(require("getopt"))
suppressMessages(require("ggplot2"))
suppressMessages(require("reshape2"))
suppressMessages(require("plyr"))
suppressMessages(require("gridExtra"))

# Suppress warning messages
options(warn=-1)

#################################################################
# UTILITY FUNCTIONS
#################################################################
# My color map using tableau colors
# http://tableaufriction.blogspot.com/2012/11/finally-you-can-use-tableau-data-colors.html
getColorMap <- function() {
    colors <- c("yellow", "blue", "orange", "green", "red", "purple", "brown",
                "light grey", "light red", "green-orange", "olive",
                "light blue", "light green", "grey", "pink", "orange-yellow",
                "dark blue", "aqua", "dark green", "black")
    colorMap <- table(colors)
    colorMap["yellow"]       <- "#ffd94a"
    colorMap["blue"]         <- "#1f77b4"
    colorMap["orange"]       <- "#ff7f0e"
    colorMap["green"]        <- "#2ca02c"
    colorMap["red"]          <- "#d62728"
    colorMap["light grey"]   <- "#c7c7c7"
    colorMap["purple"]       <- "#9467bd"
    colorMap["brown"]        <- "#8c564b"
    colorMap["light red"]    <- "#ff9896"
    colorMap["green-orange"] <- "#86b4a9"
    colorMap["olive"]        <- "#bcbd22"
    colorMap["light blue"]   <- "#aec7e8"
    colorMap["light green"]  <- "#98df8a"
    colorMap["grey"]  <- "#c7c7c7"
    colorMap["pink"]  <- "#e377c2"
    colorMap["orange-yellow"]  <- "#ffc156"
    colorMap["dark blue"]  <- "#000bca"
    colorMap["aqua"]  <- "#17becf"
    colorMap["dark green"]  <- "#006400"
    colorMap["dark red"]  <- "#8b0000"
    colorMap["black"]  <- "#000000"
    return(colorMap)
}

#### Plotting variables ###
# Load color map
cm <- getColorMap()
cmplots <- c(cm[["blue"]], cm[["red"]], cm[["green"]], cm[["purple"]], cm[["orange"]], cm[["black"]], cm[["grey"]])
# Set theme
my.theme <-
    theme(axis.text=element_text(colour="black", size=12),
          axis.title=element_text(face="bold", size=15),
          axis.title.x=element_text(margin=margin(t=10, b=5)),
          axis.ticks=element_blank(),
          axis.line.x=element_line(colour="black"),
          axis.line.y=element_line(colour="black"),
          legend.key=element_rect(fill=NA),
          legend.text=element_text(size=12),
          plot.title=element_text(size=16, face="bold"),
          panel.background=element_blank(),
          panel.margin=unit(3, "mm"),
          panel.grid.major.x=element_line(colour="#e7e7e7"),
          panel.grid.major.y=element_blank(),
          panel.grid.minor=element_blank())


# summarySE
# Gives count, mean, standard deviation, standard error of the mean,
# and confidence interval (default 95%).
#   data:          a data frame.
#   measurevar:    the name of a column that contains the
#                  variable to be summariezed
#   groupvars:     a vector containing names of columns
#                  that contain grouping variables
#   na.rm:         a boolean that indicates whether to ignore NA's
#   conf.interval: the percent range of the
#                  confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {

    # New version of length which can handle NA's:
    # if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
		if (na.rm) {
            sum(!is.na(x))
        }
		else {
            length(x)
        }
	}

	# This does the summary. For each group's data frame, return a vector with
	# N, mean, and sd
    datac <- ddply(data,
                   groupvars,
                   .drop=.drop,
                   .fun = function(xx, col) {
                       c(N = length2(xx[[col]], na.rm=na.rm),
                         mean = mean(xx[[col]], na.rm=na.rm),
                         sd = sd(xx[[col]], na.rm=na.rm))
                   },
                   measurevar)

    # Rename the "mean" column
	datac <- rename(datac, c("mean" = measurevar))

    # Calculate standard error of the mean
	datac$se <- datac$sd / sqrt(datac$N)

	# Confidence interval multiplier for standard error
	# Calculate t-statistic for confidence interval:
	# e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
	ciMult <- qt(conf.interval/2 + .5, datac$N-1)
	datac$ci <- datac$se * ciMult
	return(datac)
}


#################################################################
# ARGUMENT PARSING
#################################################################
spec <- matrix(c(
    "sf_dir",     "d", 1, "character",    "SUPER-FOCUS file directory (required)",
    "help",       "h", 0, "logical",      "This help message"
    ), ncol=5, byrow=T)

opt <- getopt(spec)

# Check if help flag was given
if (!is.null(opt$help)) {
    cat(paste(getopt(spec, usage=T), "\n"))
    q(status=1)
}
# Check for super-focus directory
if (is.null(opt$sf_dir)) {
    cat("\nSUPER-FOCUS directory not specified. Use the '-d' option.\n\n")
    cat(paste(getopt(spec, usage=T), "\n"))
    q(status=1)
}


#################################################################
# DATA PROCESSING
#################################################################
# Find the data file
# Uses the suffix '_____results__all_levels_and_function.xls'
fp <- list.files(opt$sf_dir, pattern="_____results__all_levels_and_function.xls")[[1]]
fp <- paste(opt$sf_dir, fp, sep="/")

# Load data, skip first two lines
data <- read.delim(fp, skip=2)

# Change column name
names(data)[names(data) == "Subsystem.Level.3"] <- "SS3"
names(data)[names(data) == "SEED.Function"] <- "f"

# Determine number of samples
nSamples <- (ncol(data) - 4) / 2

# Expand data for plotting
# Grab only function and values
pl.data <- melt(data[, c(3, 4, seq(ncol(data) - nSamples + 1, ncol(data)))],
                id.vars=c("SS3", "f"),
                variable.name="sample",
                value.name="rel_abundance")

# Rename sample names
pl.data$newsample <- apply(subset(pl.data, select="sample"), 1, function(x) {
    strsplit(as.character(x[1]), "..Relative")[[1]][1]
})

# Replace underscores "_" with space " " in function names
pl.data$f <- as.factor(gsub("_", " ", pl.data$f))
pl.data$f <- as.factor(gsub("^(.{85})(.*)$", "\\1\n\\2", pl.data$f))

# Create new function names with subsystems
pl.data$fullf <- apply(subset(pl.data, select=c("SS3", "f")), 1, function(x) {
    paste(x[1], x[2], sep=":\n")
})

pl.data$fullf <- as.factor(pl.data$fullf)

#################################################################
# PLOT ALL DATA
#################################################################
# Find maximum relative abundance of each function
stats <- aggregate(list(rel_abundance=pl.data$rel_abundance), by=list(fullf=pl.data$fullf), FUN=max)

# Sort data
all.sorder <- order(-stats$rel_abundance)
pl.data$fullf <- factor(pl.data$fullf, levels=stats$fullf[rev(all.sorder)])

# Top 50 functions
top <- stats$fullf[all.sorder][seq(1, 30)]

# Plot top 30 functions
pl <- ggplot(pl.data[pl.data$fullf %in% top,], aes(x=fullf, y=rel_abundance, fill=newsample)) +
    geom_bar(stat="identity", width=0.5, position=position_dodge(width=0.5)) +
    my.theme +
    scale_y_continuous(expand=c(0, 0)) +
    coord_flip() +
    scale_fill_manual("sample", values=cmplots) +
    ggtitle("Top 30 most occurring functions") + xlab("") + ylab("Relative Abundance (%)")

ggsave(paste(opt$sf_dir, "/all_top_functions.png", sep=""),
       plot=pl,
       width=40,
       height=40,
       units="cm",
       dpi=300)

#################################################################
# PLOT INDIVIDUAL
#################################################################
for (s in unique(pl.data$newsample)) {
    # Subset of specific sample name
    subdata <- subset(pl.data, newsample == s)

    #adata <- aggregate(subdata$rel_abundance, by=list(subdata$fullf), FUN=mean)

    # Get top 30 functions based on relative abundance
    top30 <- order(-subdata$rel_abundance)[seq(1, 30)]
    subdata <- subdata[top30,]
    tmp <- droplevels(subdata$fullf)
    tmp <- droplevels(subdata$SS3)

    # Order from highest to lowest relative abundance
    subdata.order <- rev(order(-subdata$rel_abundance))
    subdata$fullf <- factor(subdata$fullf, levels=subdata$fullf[subdata.order])
    subdata$SS3 <- factor(subdata$SS3, levels=subdata$SS3[subdata.order])

    pl <- ggplot(subdata, aes(x=fullf, y=rel_abundance)) +
    geom_bar(stat="identity", width=0.5, position=position_dodge(width=0.5), fill="#444444") +
    my.theme +
    scale_y_continuous(expand=c(0, 0)) +
    coord_flip() +
    scale_fill_manual("sample", values=cmplots) +
    ggtitle(paste(s, "\nTop 30 most occurring functions", sep="")) + xlab("") + ylab("Relative Abundance (%)")
    ggsave(paste(opt$sf_dir, "/", s, "_top_functions.png", sep=""),
       plot=pl,
       width=40,
       height=40,
       units="cm",
       dpi=300)

    # No longer plotting subsystems bar chart
#    # Plot subsystems of top 50 functions
#    pl <- ggplot(subdata, aes(x=SS3, y=rel_abundance)) +
#    geom_bar(stat="identity", position="dodge", fill="#444444") +
#    my.theme +
#    scale_y_continuous(expand=c(0, 0)) +
#    coord_flip() +
#    scale_fill_manual("sample", values=cmplots) +
#    ggtitle(paste(s, "\nSubsystems of top 30 most occurring functions", sep="")) + xlab("") + ylab("relative abundance (%)")
#    ggsave(paste(opt$sf_dir, "/", s, "_top_functions_ss3.png", sep=""),
#       plot=pl,
#       width=40,
#       height=40,
#       units="cm",
#       dpi=300)

    # Plot table of subsystems, functional roles, and relative abundances
    tabledata <- subset(subdata, select=c("SS3", "f", "rel_abundance"))
    tabledata$f <- factor(tabledata$f, levels=tabledata$f[rev(subdata.order)])
    tabledata$SS3 <- factor(tabledata$SS3, levels=tabledata$SS3[rev(subdata.order)])
    names(tabledata) <- c("Subsystem", "Functional Role", "Rel. Abundance (%)")
    tt <- ttheme_default(core=list(fg_params=list(hjust=1, x=1)),
                         colhead=list(fg_params=list(hjust=1, x=1)))
    tbl <- tableGrob(tabledata, rows=NULL, theme=tt)
    plt <- arrangeGrob(tbl, ncol=1, as.table=T, top=paste(s, "\nSubsystems of top 30 most occurring functions", sep=""))
    ggsave(paste(opt$sf_dir, "/", s, "_top_functions_table.png", sep=""),
           plot=plt,
           width=55,
           height=25,
           units="cm",
           dpi=300)

}

