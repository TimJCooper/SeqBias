#!/usr/bin/env Rscript
#Package Version: 3.1

######################################################################################################
# Author(s): T.J.Cooper
# Updated: 20/12/2017
# Plots sequence bias profiles on a per-strand and combined basis
######################################################################################################

suppressMessages(require(ggplot2))
suppressMessages(require(reshape2))
suppressMessages(require(gridExtra))
suppressMessages(require(grid))
args <- commandArgs(trailing = TRUE)
strand <- list("Watson", "Crick", "Combined")
outfile <- file.path(getwd(), "Plots", args[4])
pdf(sprintf("%s.pdf", outfile), width = 7, height = 9)
plot_list <- list()
for (i in 1:3) {
  DF <- read.table(args[i], sep = "\t", header = TRUE)
  cDF <- melt(DF[, c("A", "T", "G", "C", "Pos")], id = "Pos", variable.name = "Base", value.name = "Frequency")
  plot_list[[i]] <-
    ggplot(data = cDF, aes(x = Pos, y = Frequency, group = Base, colour = Base)) +
    geom_line() +
    scale_color_manual(values = c("#FFD000", "#246FAB", "#CA0017", "#78AB46")) +
    ggtitle(strand[[i]]) +
    theme(
      text = element_text(size = 7),
      axis.line.x = element_line(color = "black", size = 0.2),
      axis.line.y = element_line(color = "black", size = 0.2)
    ) +
    xlab("±bp") +
    ylab("Base Frequency")
  if (max(cDF[, "Pos"]) < 50) {
    plot_list[[i]] <- plot_list[[i]] + geom_point(size = 1.0)
  }
}
grid.arrange(grobs = plot_list, nrow = 3, top = textGrob(args[4], gp = gpar(fontsize = 8, font = 8)))
