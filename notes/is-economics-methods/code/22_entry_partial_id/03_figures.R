# ============================================================
# Chapter 18 - figures: (1) 2-firm multiplicity phase diagram,
# (2) Ciliberto-Tamer criterion Q(delta) and the identified set.
# ============================================================
suppressPackageStartupMessages({library(data.table)})
CT<-readRDS("ch22_ct.rds")

## ---------- Figure 1: multiplicity phase diagram (2-firm game) ----------
# firm i enters iff base + eps_i - delta * y_{-i} >= 0. Symmetric base, delta.
base<-0.5; delta<-1.1
lo <- -base; hi <- delta - base                 # eps thresholds: want-in-alone / want-in-together
svglite::svglite("../../assets/fig/fig_22_multiplicity.svg", width=6.4, height=6.0)
par(mar=c(4.4,4.6,2.4,1), family="sans")
plot(NA, xlim=c(-2.2,2.2), ylim=c(-2.2,2.2), xlab=expression(epsilon[1]),
     ylab=expression(epsilon[2]), axes=FALSE)
axis(1); axis(2)
# regions
rect(-2.2,-2.2,lo,lo, col="#edf2f7", border=NA)                 # (0,0)
rect(hi,hi,2.2,2.2,   col="#c6dbef", border=NA)                 # (1,1)
rect(hi,-2.2,2.2,lo,  col="#bcd4e6", border=NA)                 # (1,0) unique (firm1 in regardless)
rect(-2.2,hi,lo,2.2,  col="#bcd4e6", border=NA)                 # (0,1) unique
rect(lo,lo,hi,hi,     col="#f6c9c0", border=NA)                 # multiplicity: {(1,0),(0,1)}
# also the L-shaped unique (1,0)/(0,1) strips adjacent to the square are (1,0)/(0,1) too;
# keep the picture simple: label the key regions
abline(v=c(lo,hi), h=c(lo,hi), col="#718096", lty=3)
text(-1.5,-1.5,"(0,0)\nneither",cex=0.95,col="#2d3748")
text(1.5,1.5,"(1,1)\nboth",cex=0.95,col="#2d3748")
text(1.5,-1.5,"(1,0)",cex=0.95,col="#1a365d"); text(-1.5,1.5,"(0,1)",cex=0.95,col="#1a365d")
text(0.05,0.05,"multiple eq\n{(1,0),(0,1)}",cex=0.92,col="#c53030",font=2)
title(main="Entry game with 2 firms: a region of multiple equilibria",
      cex.main=1.0,font.main=1,col.main="#1a365d")
box(col="#a0aec0")
dev.off(); cat("wrote fig_22_multiplicity.svg\n")

## ---------- Figure 2: CT criterion Q(delta) and identified set ----------
svglite::svglite("../../assets/fig/fig_22_identifiedset.svg", width=8.2, height=4.4)
par(mar=c(4.4,4.6,2.4,1), family="sans")
plot(CT$grid, CT$Q*1e3, type="o", pch=16, lwd=2, col="#1a365d",
     xlab=expression(paste("competitive effect  ", delta)),
     ylab="criterion  Q(δ) × 1000", axes=FALSE)
axis(1); axis(2)
inset<-CT$inset
rect(min(inset)-0.05, -1, max(inset)+0.05, 0.4, col="#c6f6d533", border=NA)  # identified set band
abline(v=1.1, col="#c53030", lwd=2, lty=2)                                    # truth
text(1.1, max(CT$Q*1e3)*0.7, " truth δ=1.1", col="#c53030", pos=4, cex=0.9)
text(mean(inset), 0.9, paste0("identified set [", min(inset), ", ", max(inset), "]"),
     col="#2f855a", cex=0.92, font=2)
title(main="Partial identification: the criterion floors over an interval, not a point",
      cex.main=0.98, font.main=1, col.main="#1a365d")
dev.off(); cat("wrote fig_22_identifiedset.svg\n")
