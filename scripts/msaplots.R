

#' Plot a multiple sequence alignment
#' 
#' Scans through all the SNP and gap files in the lists of files 
#' and for each sequence plots the regions that are covered as 
#' rectangles and indicate SNPs and ambiguous sites.
#' The input files should be produced by msasummary.py
#' 
#' @param snpFiles    List of filenames containing SNP and ambiguous site positions
#'                    for each sequence, produced by msasummary.py
#' @param gapFiles    List of filenames containing gap locations for each sequence, 
#'                    produced by msasummary.py. Files should be in the same order as
#'                    \code{snpFiles}
#' @param seqlen      The length of the reference sequence / alignment.
#' @param names       Names of sequences, in the same order as the filenames. If NULL
#'                    names will be taken from the files.
#' @param barwidth    Height of the rectangles. Values > 1 means sequences will overlap.
#' @param plotNames   Add sequence names at the right of the plot
#' @param plotGrid    Plot vertical grid lines
#' @param plotXAxis   Plot an X-axis
#' @param plotStats   Add coverage and SNP statistics to the sequence names
#' @param seqCol      Colour of the rectangles. If this is a vector of the same length 
#'                    as the number of sequences a different colour will be used for each
#'                    sequence.
#' @param seqBorder   Border colour of the rectangles. By default not plotted. Can also be
#'                    a vector.
#' @param bgCol       Colour of a thinner background rectangle for each sequence.
#' @param snpCol      Colour of SNP positions
#' @param ambCol      Colour of ambiguous sites
#' @param gridCol     Colour of the grid lines
#' @param add         If TRUE do not start a new plot
#' @param cex.axis 
#' @param ...         Extra options passed to plot()  
#' 
#' @return            A list of the number of SNPs at each site in the alignment. Includes
#'                    ambiguous sites.
#' 
plotAlignment <- function(snpFiles, gapFiles, seqlen, names=NULL, barwidth=0.7, 
                          plotNames=TRUE, plotGrid=TRUE, plotXAxis=TRUE, plotStats=TRUE, 
                          seqCol=mPal(oxCols$blue1), seqBorder=NA, bgCol=mPal(oxCols$gray3), 
                          snpCol=mPal(oxCols$red2), ambCol=mPal(oxCols$orange2), gridCol=mPal(oxCols$gray3),
                          add=FALSE, cex.axis=0.6, ...) {
  
  allSNPs  <- list()
  nrSeqs   <- length(snpFiles)
  barwidth <- barwidth/2
  
  if (length(seqCol) < nrSeqs) {
      seqCol <- rep(seqCol, nrSeqs)
  }
  
  if (length(seqBorder) < nrSeqs) {
      seqBorder <- rep(seqBorder, nrSeqs)
  }
  
  if (add == FALSE) {
      plot(1, type='n', ylim=c(0,nrSeqs+1), xlim=c(1,seqlen), xaxs='i', yaxs='i', 
           yaxt='n', xlab="", ylab="", axes=FALSE, ...)
  }
  
  if (plotXAxis) {
      axis(1, at=c(0, axTicks(1), seqlen), cex.axis=cex.axis, xpd=TRUE)
  }
  
  
  labels=c()
  for (i in 1:nrSeqs) {
    
      #############
      # Read data #
      #############
    
      print(snpFiles[i])
      
      # Read SNPs
      snps     <- read.csv(snpFiles[i], skip=1, colClasses = c("numeric", "character", "character"))
      snpIdxs  <- which(snps[, 3] %in% c("A","C","G","T"))
      ambIdxs  <- setdiff(1:nrow(snps), snpIdxs)
      nrSnps   <- length(snpIdxs)
      nrAmbig  <- length(ambIdxs)
      
      for (j in 1:nrow(snps)) {
          if (!is.na(snps[j,1]) && length(snps[j,2]) > 0) {
              snp <- paste(snps[j,c(2,1,3)], collapse=" ")
              if (snp %in% names(allSNPs)) {
                  allSNPs[[snp]] <- allSNPs[[snp]] + 1
              } else {
                  allSNPs[[snp]] <- 1
              }
          }
      }

      print(gapFiles[i])
    
      # Read gaps and complement
      # (probably needs a little more robust testing)
      gaps <- read.csv(gapFiles[i], skip=1)
      if (nrow(gaps) == 0) {
          coverage <- 1
          contigs  <- data.frame(start=1, end=seqlen)
      } else {
          coverage <- 1 - sum(gaps$length)/seqlen
          
          contigs <- data.frame(start=(gaps$end[1:(nrow(gaps)-1)]+1), end=(gaps$start[2:nrow(gaps)]-1))
          if (gaps$start[1] > 1)
            contigs <- rbind(c(1, gaps$start[1]-1), contigs)
          if (gaps$end[nrow(gaps)] != seqlen)
            contigs <- rbind(contigs, c(gaps$end[nrow(gaps)]+1, seqlen))
      }
    
    
      ########
      # Plot #
      ########
      
      # Background
      rect(1, i-barwidth/2, seqlen, i+barwidth/2, col=bgCol, border=NA)
      
      # Rectangles
      for (j in 1:nrow(contigs)) {
          rect(contigs$start[j], i-barwidth, contigs$end[j], i+barwidth, col=seqCol[i], border=seqBorder[i], lwd=0.5, xpd=TRUE)
      }
      
      # SNPs and ambiguous sites
      if (nrSnps > 0)  segments(snps$position[snpIdxs], i-barwidth, snps$position[snpIdxs], i+barwidth, col=snpCol, lwd=1.5)
      if (nrAmbig > 0) segments(snps$position[ambIdxs], i-barwidth, snps$position[ambIdxs], i+barwidth, col=ambCol, lwd=1.5)
      
      # Labels
      if (is.null(names)) {
          seqid <- gsub("\\.","-",colnames(snps[3]))
      } else {
          seqid <- names[i]
      }
      
      if (plotStats) {
          seqid <- paste0(seqid, " (", nrSnps, "/", nrAmbig,", ", round(coverage*100,2), "%)")
      }
      labels <- c(labels, seqid)
      
  }
  
  if (plotGrid) {
      abline(v=seq(1000,seqlen,by=1000), lty=3, lwd=0.5, col=gridCol)
  }
  
  # Labels
  if (plotNames) {
      axis(4, at=1:nrSeqs, labels=labels, las=1, lwd.ticks=NA, lwd=NA, cex.axis=cex.axis)
  }
  
  
  return( allSNPs )
  
}




plotSNPHist <- function(snpTable, cutoff=0, ylim=NULL, plotGrid=TRUE, 
                        col=mPal(oxCols$gray6), labelCol=mPal(oxCols$red2), gridCol=mPal(oxCols$gray3), 
                        cex.axis=0.8) {
  
    x <- y <- c()
    if (length(snpTable) > 0) {
        for (snp in names(snpTable)) {
            x <- c(x, as.numeric(strsplit(snp, split=" ")[[1]][2]))
            y <- c(y, snpTable[[snp]])
        }
        
        if (is.null(ylim)) {
            ylim   <- range(pretty(c(0, 2+round(max(unlist(snpTable))))))
        }
    } else {
        ylim <- c(0,1)
    }
    
    plot(1, type='n', ylim=ylim, xlim=c(1,seqlen), xaxs='i', yaxs='i', 
         yaxt='n', xlab="", ylab="", bty='n', axes=FALSE)

    yticks <- unique(floor(axTicks(2)))
    
    if (plotGrid) {
        abline(h=yticks[-1], lty=3, lwd=0.5, col=gridCol)
        #grid(nx=0, ny=NULL, lwd=0.5, col=gridCol)
    }
    
    # SNPs
    #print(col)
    if (length(x) > 0) {
        segments(x, rep(0, length(y)), x, y, col=col, lwd=1.5)
    }
    
    axis(1, at=c(0, axTicks(1), seqlen), cex.axis=cex.axis, xpd=TRUE)
    axis(2, las=1, cex.axis=cex.axis, at=yticks)
    
    
    # Label SNPs > cutoff
    if (max(y) > cutoff) {
        text(x[y > cutoff], y[y > cutoff]+(0.15*ylim[2]), names(snpTable)[y > cutoff], offset = c(0,3), srt=45, cex=cex.axis, col=labelCol, xpd=TRUE)
    }
  
}
