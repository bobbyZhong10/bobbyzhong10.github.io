R <- readRDS("nova_results.rds")
dir.create("../../assets/fig", recursive=TRUE, showWarnings=FALSE)

# Axis labels follow the text's notation R^2_{D~U|X} and R^2_{Y~U|D,X}. The tilde and bar are plain
# strings: plotmath's %~% (U+223C) has no glyph in the svg() font and printed as an empty box.
z <- matrix(R$grid$adjusted, nrow=81, ncol=81)
contour_plot <- function() {
  par(mar=c(5.1,5.1,4.1,2.1))
  filled.contour(seq(0,0.30,length.out=81), seq(0,0.30,length.out=81), z,
    color.palette=function(n) hcl.colors(n,"Blues 3",rev=TRUE),
    xlab=expression(partial~R[D*" ~ "*U*" | "*X]^2),
    ylab=expression(partial~R[Y*" ~ "*U*" | "*D*", "*X]^2),
    key.title=title(main="Adjusted\neffect",cex.main=0.8))
}
svg("../../assets/fig/fig_23_robustness_contour.svg", width=7.2, height=5.6)
contour_plot()
dev.off()
# The SVG is larger than the 300 KB threshold of inline_figures.py, so the rendered page embeds
# this PNG sibling instead; draw it from the same function so the two never drift apart.
png("../../assets/fig/fig_23_robustness_contour.png", width=7.2, height=5.6, units="in", res=216,
    type="cairo")
contour_plot()
dev.off()

svg("../../assets/fig/fig_23_identification_region.svg", width=8.2, height=4.6)
plot(R$region$bias_cap, R$region$lower, type="l", lwd=2, col="#c53030",
     ylim=range(0,R$region$lower,R$region$upper,R$truth), xlab="Allowed absolute omitted-variable bias",
     ylab="Identified set for treatment effect")
lines(R$region$bias_cap,R$region$upper,lwd=2,col="#2b6cb0")
polygon(c(R$region$bias_cap,rev(R$region$bias_cap)),
        c(R$region$lower,rev(R$region$upper)),col="#bee3f855",border=NA)
abline(h=0,lty=2,col="#4a5568"); abline(h=R$truth,lty=3,lwd=2,col="#4a7c59")
legend("topleft",c("Lower bound","Upper bound","True effect"),
       col=c("#c53030","#2b6cb0","#4a7c59"),lty=c(1,1,3),lwd=2,bty="n")
dev.off()
cat(sprintf("breakdown bias for including zero = %.3f | for including the true effect = %.3f | lower bound at the largest cap = %.3f\n",
            abs(R$naive[1]), R$naive[1]-R$truth, min(R$region$lower)))
