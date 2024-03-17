.cumsum0 <- function(x, left = TRUE, right = FALSE, n = NULL) {
    xx <- c(0, cumsum(as.numeric(x)))
    if (!left) 
        xx <- xx[-1]
    if (!right) 
        xx <- head(xx, -1)
    names(xx) <- n
    xx
}

get_gains <- function(gains, pos, sample_ids, h){
    gains_new <- gains
    for (i in 1:length(pos)){
        if (i==1){
            start <- 0
        }
        else {
            start <- pos[i-1]
        }

        if (i<44){
            end <- pos[i]
            chrs_1 <- h[(h>=start)&(h<=end)]
        }
        else{
            chrs_1 <- h[h>=start]
        }

        v_sum <- gains[names(chrs_1)] #sum(gains[names(chrs_1)])
        v_len <- 1 #length(gains[names(chrs_1)])
        v_prop <- v_sum/(v_len*length(sample_ids))

        gains_new[names(chrs_1)] <- v_prop
    }
    return(gains_new)
}

get_gains_regions <- function(gains_regions, sample_ids){
    gains_regions_new <- gains_regions/length(sample_ids)
    return(gains_regions_new)
}

composite_cnv_dist_gr <- function(high_sample_ids, ctrl_ids, anno, minfi.data, savename,
                         ylim, cols, chr, chrX, chrY, centromere, main){
    
    ########### Gains
    id_n = 1
    # message(id_n)
    colors = c("green", "skyblue")
    pdf(savename, width=20)
    for (ids in list(high_sample_ids)){
        
        for (sample in ids){
            message(paste("processing sample", sample))
            x <- CNV.fit(minfi.data[sample], minfi.data[ctrls], anno)
            x <- CNV.bin(x)
            x <- CNV.detail(x)
            x <- CNV.segment(x)

            bin.ratio <- x@bin$ratio - x@bin$shift #  Adjust for shift, different per sample
            detail.ratio <- x@detail$ratio - x@bin$shift

            bin.ratio[bin.ratio < ylim[1]] <- ylim[1] # Clip values 
            bin.ratio[bin.ratio > ylim[2]] <- ylim[2]
            
            if (sample==ids[1]){
                if (sample==high_sample_ids[1]){
                    message("Starting")
                    if (is.null(main)) 
                        main <- x@name
                    if (chr[1] == "all") {
                        chr <- x@anno@genome$chr
                    } else {
                        chr <- intersect(chr, x@anno@genome$chr)
                    }

                    plot(NA, xlim = c(0, sum(as.numeric(x@anno@genome[chr, "size"])) - 
                    0), ylim = ylim, xaxs = "i", xaxt = "n", yaxt = "n", xlab = NA, 
                    ylab = NA, main = main)

                    chr.cumsum0 <- .cumsum0(x@anno@genome[chr, "size"], n = chr)

                    if (!chrX & is.element("chrX", names(chr.cumsum0))) 
                            chr.cumsum0["chrX"] <- NA
                    if (!chrY & is.element("chrY", names(chr.cumsum0))) 
                            chr.cumsum0["chrY"] <- NA

                    abline(v = .cumsum0(x@anno@genome[chr, "size"], right = TRUE), 
                            col = "grey")

                    if (centromere) {
                        abline(v = .cumsum0(x@anno@genome[chr, "size"]) + x@anno@genome[chr, 
                            "pq"], col = "grey", lty = 2) 
                    }

                    axis(1, at = .cumsum0(x@anno@genome[chr, "size"]) + x@anno@genome[chr, 
                        "size"]/2, labels = x@anno@genome[chr, "chr"], las = 2)
                    if (all(ylim == c(-1.25, 1.25))) {
                        axis(2, at = round(seq(-1.2, 1.2, 0.2), 1), las = 2)
                    } else {
                        axis(2, las = 2)
                }
                }
                gains <- bin.ratio
                gains[gains>=0.4] = 1
                gains[gains<0.4] = 0
                
                gains_regions <- detail.ratio
                gains_regions[gains_regions>=0.4] = 1
                gains_regions[gains_regions<0.4] = 0
            }
            else {
                # gains
                gains_tmp <- bin.ratio
                gains_tmp[gains_tmp>=0.4] = 1
                gains_tmp[gains_tmp<0.4] = 0

                gains_regions_tmp <- detail.ratio
                gains_regions_tmp[gains_regions_tmp>=0.4] = 1
                gains_regions_tmp[gains_regions_tmp<0.4] = 0
                
                gains <- gains + gains_tmp
                gains_regions <- gains_regions + gains_regions_tmp
            
            }
        }
        
        # chr
        h = chr.cumsum0[as.vector(seqnames(x@anno@bins))] + values(x@anno@bins)$midpoint
        names(h) <- names(gains)
        pq_pos = .cumsum0(x@anno@genome[chr, "size"]) + x@anno@genome[chr, 
                            "pq"]
        chr_pos <- .cumsum0(x@anno@genome[x@anno@genome$chr, "size"])
        pos = c(rbind(chr_pos, pq_pos))

        gains_new <- get_gains(gains, pos, ids, h)
        
        
        lines(chr.cumsum0[as.vector(seqnames(x@anno@bins))] + values(x@anno@bins)$midpoint, 
                    gains_new, type = "h", pch = 16, cex = 0.4, col=colors[id_n])
        
        # gene regions
        
        gains_regions_new <- get_gains_regions(gains_regions, ids)
        # detail.ratio.above <- (gains_regions_new > 0 & gains_regions_new < 0.85) | 
        #                         gains_regions_new < -0.85
        
        lines(start(x@anno@detail) + (end(x@anno@detail) - start(x@anno@detail)) /2
              + chr.cumsum0[as.vector(seqnames(x@anno@detail))], 
                gains_regions_new, type = "p", pch = 21, col="black", bg=colors[id_n])
                #bg=alpha(colors[id_n], 0.5))
        
        if (id_n==1){ # Label only once
            labels <- values(x@anno@detail)$name
            labels <- gsub("\\/", "-", labels)
            print(labels)
            
            text(start(x@anno@detail) + (end(x@anno@detail) - start(x@anno@detail)) /2
                 + chr.cumsum0[as.vector(seqnames(x@anno@detail))], gains_regions_new, # ifelse(detail.ratio.above, detail.ratio, NA), 
                 labels = labels, adj = c(0, 
                    0.5), srt = 90, col = "darkblue")

            # text(start(x@anno@detail) + (end(x@anno@detail) - start(x@anno@detail)) /2
            #      + chr.cumsum0[as.vector(seqnames(x@anno@detail))], gains_regions_new, # ifelse(detail.ratio.above, NA, detail.ratio), 
            #     labels = labels, adj = c(1, 0.5), srt = 90, col = "darkblue")
        }
        
        id_n <- 2
        
        }

    ########### Losses
    id_n = 1
    colors = c("red", "blue")
    for (ids in list(high_sample_ids)){
        for (sample in ids){
            message(paste("processing sample", sample))
            x <- CNV.fit(minfi.data[sample], minfi.data[ctrls], anno)
            x <- CNV.bin(x)
            x <- CNV.detail(x)
            x <- CNV.segment(x)

            bin.ratio <- x@bin$ratio - x@bin$shift #  Adjust for shift, different per sample
            detail.ratio <- x@detail$ratio - x@bin$shift

            bin.ratio[bin.ratio < ylim[1]] <- ylim[1] # Clip values 
            bin.ratio[bin.ratio > ylim[2]] <- ylim[2]
            
            if (sample==ids[1]){
                
                gains <- bin.ratio
                gains[gains>-0.4] = 0
                gains[gains<=-0.4] = 1
                
                gains_regions <- detail.ratio
                gains_regions[gains_regions>-0.4] = 0
                gains_regions[gains_regions<=-0.4] = 1
            }
            else {
                # gains
                gains_tmp <- bin.ratio
                gains_tmp[gains_tmp>-0.4] = 0
                gains_tmp[gains_tmp<=-0.4] = 1
                
                gains_regions_tmp <- detail.ratio
                gains_regions_tmp[gains_regions_tmp>-0.4] = 0
                gains_regions_tmp[gains_regions_tmp<=-0.4] = 1
    
                gains <- gains + gains_tmp
                gains_regions <- gains_regions + gains_regions_tmp
            
            }
        }
        
        # chr
        h = chr.cumsum0[as.vector(seqnames(x@anno@bins))] + values(x@anno@bins)$midpoint
        names(h) <- names(gains)
        pq_pos = .cumsum0(x@anno@genome[chr, "size"]) + x@anno@genome[chr, 
                            "pq"]
        chr_pos <- .cumsum0(x@anno@genome[x@anno@genome$chr, "size"])
        pos = c(rbind(chr_pos, pq_pos))

        gains_new <- get_gains(gains, pos, ids, h)
        
        lines(chr.cumsum0[as.vector(seqnames(x@anno@bins))] + values(x@anno@bins)$midpoint, 
                    -gains_new, type = "h", pch = 16, cex = 0.4, col=colors[id_n])
        
        # gene regions
        gains_regions_new <- get_gains_regions(gains_regions, ids)
        # detail.ratio.above <- (gains_regions_new > 0 & gains_regions_new < 0.85) | 
        #                         gains_regions_new < -0.85
        
        lines(start(x@anno@detail) + (end(x@anno@detail) - start(x@anno@detail)) /2
              + chr.cumsum0[as.vector(seqnames(x@anno@detail))], 
                -gains_regions_new, type = "p", pch = 21, col="black", bg=colors[id_n])
                # bg=alpha(colors[id_n], 0.5))
        
        id_n <- 2
        }

    
        legend("topright", inset=c(-0.2,0), bg="transparent", legend=c("Low neuronal (gains)",
                                                          "High neuronal (gains)",
                                                          "Low neuronal (losses)",
                                                          "High neuronal (losses)"), 
           fill=c("skyblue", "green", "blue", "red"), title="Group")
    
    dev.off()
    }


composite_cnv_dist_gr_1 <- function(high_sample_ids, low_sample_ids, cell_culture, h_cortex, neun, ctrl_ids, 
                                anno, minfi.data, savename, ylim, cols, chr, chrX, chrY, centromere, main){
    ########### Gains
    id_n = 1
    colors = c("green", "skyblue", "gray", "pink", "red")
    pdf(savename, width=20)
    for (ids in list(high_sample_ids, low_sample_ids, cell_culture, h_cortex, neun)){
        for (sample in ids){
            message(paste("processing sample", sample))
            x <- CNV.fit(minfi.data[sample], minfi.data[ctrls], anno)
            x <- CNV.bin(x)
            x <- CNV.detail(x)
            x <- CNV.segment(x)

            bin.ratio <- x@bin$ratio - x@bin$shift #  Adjust for shift, different per sample
            detail.ratio <- x@detail$ratio - x@bin$shift

            bin.ratio[bin.ratio < ylim[1]] <- ylim[1] # Clip values 
            bin.ratio[bin.ratio > ylim[2]] <- ylim[2]
            
            if (sample==ids[1]){
                if (sample==high_sample_ids[1]){
                    message("Starting")
                    if (is.null(main)) 
                        main <- x@name
                    if (chr[1] == "all") {
                        chr <- x@anno@genome$chr
                    } else {
                        chr <- intersect(chr, x@anno@genome$chr)
                    }

                    plot(NA, xlim = c(0, sum(as.numeric(x@anno@genome[chr, "size"])) - 
                    0), ylim = ylim, xaxs = "i", xaxt = "n", yaxt = "n", xlab = NA, 
                    ylab = NA, main = main)

                    chr.cumsum0 <- .cumsum0(x@anno@genome[chr, "size"], n = chr)

                    if (!chrX & is.element("chrX", names(chr.cumsum0))) 
                            chr.cumsum0["chrX"] <- NA
                    if (!chrY & is.element("chrY", names(chr.cumsum0))) 
                            chr.cumsum0["chrY"] <- NA

                    abline(v = .cumsum0(x@anno@genome[chr, "size"], right = TRUE), 
                            col = "grey")

                    if (centromere) {
                        abline(v = .cumsum0(x@anno@genome[chr, "size"]) + x@anno@genome[chr, 
                            "pq"], col = "grey", lty = 2) 
                    }

                    axis(1, at = .cumsum0(x@anno@genome[chr, "size"]) + x@anno@genome[chr, 
                        "size"]/2, labels = x@anno@genome[chr, "chr"], las = 2)
                    if (all(ylim == c(-1.25, 1.25))) {
                        axis(2, at = round(seq(-1.2, 1.2, 0.2), 1), las = 2)
                    } else {
                        axis(2, las = 2)
                }
                }
                gains <- bin.ratio
                gains[gains>=0.4] = 1
                gains[gains<0.4] = 0
                
                gains_regions <- detail.ratio
                gains_regions[gains_regions>=0.4] = 1
                gains_regions[gains_regions<0.4] = 0
            }
            else {
                # gains
                gains_tmp <- bin.ratio
                gains_tmp[gains_tmp>=0.4] = 1
                gains_tmp[gains_tmp<0.4] = 0

                gains_regions_tmp <- detail.ratio
                gains_regions_tmp[gains_regions_tmp>=0.4] = 1
                gains_regions_tmp[gains_regions_tmp<0.4] = 0
                
                gains <- gains + gains_tmp
                gains_regions <- gains_regions + gains_regions_tmp
            
            }
        }
        
        # chr
        h = chr.cumsum0[as.vector(seqnames(x@anno@bins))] + values(x@anno@bins)$midpoint
        names(h) <- names(gains)
        pq_pos = .cumsum0(x@anno@genome[chr, "size"]) + x@anno@genome[chr, 
                            "pq"]
        chr_pos <- .cumsum0(x@anno@genome[x@anno@genome$chr, "size"])
        pos = c(rbind(chr_pos, pq_pos))

        gains_new <- get_gains(gains, pos, ids, h)
        
        
        lines(chr.cumsum0[as.vector(seqnames(x@anno@bins))] + values(x@anno@bins)$midpoint, 
                    gains_new, type = "h", pch = 16, cex = 0.4, col=colors[id_n])
        
        # gene regions
        
        gains_regions_new <- get_gains_regions(gains_regions, ids)
        # detail.ratio.above <- (gains_regions_new > 0 & gains_regions_new < 0.85) | 
        #                         gains_regions_new < -0.85
        
        lines(start(x@anno@detail) + (end(x@anno@detail) - start(x@anno@detail)) /2
              + chr.cumsum0[as.vector(seqnames(x@anno@detail))], 
                gains_regions_new, type = "p", pch = 21, col="black", bg=colors[id_n])
                #bg=alpha(colors[id_n], 0.5))
        
        if (id_n==1){ # Label only once
            labels <- values(x@anno@detail)$name
            labels <- gsub("\\/", "-", labels)
            print(labels)
            
            text(start(x@anno@detail) + (end(x@anno@detail) - start(x@anno@detail)) /2
                 + chr.cumsum0[as.vector(seqnames(x@anno@detail))], gains_regions_new, # ifelse(detail.ratio.above, detail.ratio, NA), 
                 labels = labels, adj = c(0, 
                    0.5), srt = 90, col = "darkblue")

            # text(start(x@anno@detail) + (end(x@anno@detail) - start(x@anno@detail)) /2
            #      + chr.cumsum0[as.vector(seqnames(x@anno@detail))], gains_regions_new, # ifelse(detail.ratio.above, NA, detail.ratio), 
            #     labels = labels, adj = c(1, 0.5), srt = 90, col = "darkblue")
        }
        
        id_n <- id_n + 1
        
        }

    ########### Losses
    id_n = 1
    for (ids in list(high_sample_ids, low_sample_ids, cell_culture, h_cortex, neun)){
        for (sample in ids){
            message(paste("processing sample", sample))
            x <- CNV.fit(minfi.data[sample], minfi.data[ctrls], anno)
            x <- CNV.bin(x)
            x <- CNV.detail(x)
            x <- CNV.segment(x)

            bin.ratio <- x@bin$ratio - x@bin$shift #  Adjust for shift, different per sample
            detail.ratio <- x@detail$ratio - x@bin$shift

            bin.ratio[bin.ratio < ylim[1]] <- ylim[1] # Clip values 
            bin.ratio[bin.ratio > ylim[2]] <- ylim[2]
            
            if (sample==ids[1]){
                
                gains <- bin.ratio
                gains[gains>-0.4] = 0
                gains[gains<=-0.4] = 1
                
                gains_regions <- detail.ratio
                gains_regions[gains_regions>-0.4] = 0
                gains_regions[gains_regions<=-0.4] = 1
            }
            else {
                # gains
                gains_tmp <- bin.ratio
                gains_tmp[gains_tmp>-0.4] = 0
                gains_tmp[gains_tmp<=-0.4] = 1
                
                gains_regions_tmp <- detail.ratio
                gains_regions_tmp[gains_regions_tmp>-0.4] = 0
                gains_regions_tmp[gains_regions_tmp<=-0.4] = 1
    
                gains <- gains + gains_tmp
                gains_regions <- gains_regions + gains_regions_tmp
            
            }
        }
        
        # chr
        h = chr.cumsum0[as.vector(seqnames(x@anno@bins))] + values(x@anno@bins)$midpoint
        names(h) <- names(gains)
        pq_pos = .cumsum0(x@anno@genome[chr, "size"]) + x@anno@genome[chr, 
                            "pq"]
        chr_pos <- .cumsum0(x@anno@genome[x@anno@genome$chr, "size"])
        pos = c(rbind(chr_pos, pq_pos))

        gains_new <- get_gains(gains, pos, ids, h)
        
        lines(chr.cumsum0[as.vector(seqnames(x@anno@bins))] + values(x@anno@bins)$midpoint, 
                    -gains_new, type = "h", pch = 16, cex = 0.4, col=colors[id_n])
        
        # gene regions
        gains_regions_new <- get_gains_regions(gains_regions, ids)
        # detail.ratio.above <- (gains_regions_new > 0 & gains_regions_new < 0.85) | 
        #                         gains_regions_new < -0.85
        
        lines(start(x@anno@detail) + (end(x@anno@detail) - start(x@anno@detail)) /2
              + chr.cumsum0[as.vector(seqnames(x@anno@detail))], 
                -gains_regions_new, type = "p", pch = 21, col="black", bg=colors[id_n])
                # bg=alpha(colors[id_n], 0.5))
        
        id_n <- id_n + 1
        }

    
        legend("topright", inset=c(0), bg="transparent", legend=c("Low neuronal",
                                                          "High neuronal",
                                                          "cell_culture",
                                                          "healthy cortex",
                                                          "NeuN+"),
           fill=colors, title="Group")
    
    dev.off()
    }
