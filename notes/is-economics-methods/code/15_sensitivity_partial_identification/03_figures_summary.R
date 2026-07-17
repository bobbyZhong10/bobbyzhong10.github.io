R <- readRDS("nova_results.rds")
dir.create("../../assets/fig", recursive=TRUE, showWarnings=FALSE)

z <- matrix(R$grid$adjusted, nrow=81, ncol=81)
svg("../../assets/fig/fig_15_robustness_contour.svg", width=7.2, height=5.6)
filled.contour(seq(0,0.30,length.out=81), seq(0,0.30,length.out=81), z,
  color.palette=function(n) hcl.colors(n,"Blues 3",rev=TRUE),
  xlab=expression(partial~R^2[D %~% U %.% X]),
  ylab=expression(partial~R^2[Y %~% U %.% (D+X)]),
  key.title=title(main="Adjusted\neffect",cex.main=0.8))
dev.off()

svg("../../assets/fig/fig_15_identification_region.svg", width=8.2, height=4.6)
plot(R$region$bias_cap, R$region$lower, type="l", lwd=2, col="#c53030",
     ylim=range(R$region$lower,R$region$upper,R$truth), xlab="Allowed absolute omitted-variable bias",
     ylab="Identified region for treatment effect")
lines(R$region$bias_cap,R$region$upper,lwd=2,col="#2b6cb0")
polygon(c(R$region$bias_cap,rev(R$region$bias_cap)),
        c(R$region$lower,rev(R$region$upper)),col="#bee3f855",border=NA)
abline(h=0,lty=2,col="#4a5568"); abline(h=R$truth,lty=3,lwd=2,col="#4a7c59")
legend("topright",c("Lower bound","Upper bound","True effect"),
       col=c("#c53030","#2b6cb0","#4a7c59"),lty=c(1,1,3),lwd=2,bty="n")
dev.off()
cat(sprintf("breakdown bias for including zero = %.3f\n", abs(R$naive[1])))
