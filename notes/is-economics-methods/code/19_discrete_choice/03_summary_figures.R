# ============================================================
# Chapter 15 - naive-logit elasticity, summary table, and figures.
# ============================================================
suppressPackageStartupMessages({library(data.table)})
set.seed(15)
D <- readRDS("ch19_data.rds"); E <- readRDS("ch19_estimates.rds")
grid<-D$grid; truth<-D$truth
rbar<-E$rbar; pbar<-E$pbar

## naive plain-logit own elasticity & diversion (uses biased alpha)
b<-D$est$logit
shares_plain<-function(pv,par){V<-c(0,par[1:4]+par[5]*rbar+par[6]*pv);e<-exp(V-max(V));e/sum(e)}
div_plain<-function(par){
  s0<-shares_plain(pbar,par);p1<-pbar;p1[1]<-p1[1]+1;s1<-shares_plain(p1,par);ds<-s1-s0
  list(own=(ds[2]/s0[2])/(1/pbar[1]), div=setNames(-ds[-2]/ds[2],c("outside","B","C","D")), s0=s0)
}
naive<-div_plain(b)
cat("NAIVE plain logit: own_elas A =",round(naive$own,3),
    " s_A =",round(naive$s0[2],3),"\n")
cat("  diversion:", round(naive$div,3),"\n")
cat("CF   plain logit: own_elas A =",round(E$div$plain$own_elas,3),"\n")
cat("TRUTH own_elas A =",round(E$div$truth$own_elas,3)," mixed =",round(E$div$mixed$own_elas,3),"\n")

## ---------- Figure 1: diversion from A (IIA failure and its fix) ----------
tr<-E$div$truth$div; pl<-E$div$plain$div; nl<-E$div$nested$div; mx<-E$div$mixed$div
M<-rbind(TRUTH=tr, `plain logit`=pl, nested=nl, mixed=mx)
cols<-c("#1a365d","#c53030","#2b6cb0","#6b8e23")     # navy, red, blue, olive
svglite::svglite("../../assets/fig/fig_19_substitution.svg", width=8.2, height=4.6)
par(mar=c(4,4.3,2.4,1), family="sans")
bp<-barplot(M, beside=TRUE, col=cols, border=NA, ylim=c(0,0.42),
    names.arg=c("outside","B (same nest)","C","D"),
    ylab="diversion ratio from A", cex.names=1.0, cex.axis=0.9)
title(main="Diversion from A after a $1 price increase", cex.main=1.1, font.main=1, col.main="#1a365d")
legend("topright", rownames(M), fill=cols, border=NA, bty="n", cex=0.95)
abline(h=0)
dev.off()
cat("wrote fig_19_substitution.svg\n")

## ---------- Figure 2: recovered heterogeneity in price sensitivity ----------
am<-E$est_m["alpha_bar"]; sm<-E$est_m["sigma_a"]
svglite::svglite("../../assets/fig/fig_19_heterogeneity.svg", width=8.2, height=4.4)
par(mar=c(4.2,4.3,2.4,1), family="sans")
xs<-seq(0,0.65,length=400)
plot(xs, dnorm(xs,truth$alpha_bar,truth$sigma_a), type="l", lwd=2.4, col="#1a365d",
     xlab=expression(paste("price sensitivity  ", alpha[i], "  (utils / $)")), ylab="density",
     main="", axes=FALSE)
axis(1,cex.axis=0.9); axis(2,cex.axis=0.9)
lines(xs, dnorm(xs,am,sm), lwd=2.4, col="#c53030", lty=2)
abline(v=-D$est$logit["price_coef"], col="#718096", lwd=2, lty=3)   # naive point
abline(v=-D$est$cf["price_coef"],   col="#2b6cb0", lwd=2, lty=3)    # cf point
legend("topright", bty="n", cex=0.9,
  legend=c("truth  N(0.30, 0.12)","mixed logit estimate","plain logit point (naive)","control-function point"),
  col=c("#1a365d","#c53030","#718096","#2b6cb0"), lwd=c(2.4,2.4,2,2), lty=c(1,2,3,3))
title(main="Mixed logit recovers the price-sensitivity heterogeneity plain logit averages away",
      cex.main=0.98, font.main=1, col.main="#1a365d")
dev.off()
cat("wrote fig_19_heterogeneity.svg\n")

## ---------- Figure 3: nesting-structure schematic ----------
svglite::svglite("../../assets/fig/fig_19_nest.svg", width=7.6, height=3.6)
par(mar=c(0.5,0.5,0.5,0.5)); plot.new(); plot.window(xlim=c(0,10),ylim=c(0,6))
box(col="white")
# root
rect(4.0,5,6.0,5.8,border="#1a365d",lwd=2); text(5,5.4,"one order",col="#1a365d",cex=1)
# outside + two nests
draw<-function(x,y,lab,col){rect(x-1.0,y-0.45,x+1.0,y+0.45,border=col,lwd=2);text(x,y,lab,col=col,cex=0.92)}
draw(1.4,3.4,"outside\n(no order)","#718096")
draw(5,3.4,"nest 1\ncomfort","#c53030")
draw(8.6,3.4,"nest 2\nasian","#2b6cb0")
segments(5,5,1.4,3.85,col="#718096"); segments(5,5,5,3.85,col="#c53030"); segments(5,5,8.6,3.85,col="#2b6cb0")
# leaves
draw(4.0,1.4,"A","#c53030"); draw(6.0,1.4,"B","#c53030")
draw(7.7,1.4,"C","#2b6cb0"); draw(9.5,1.4,"D","#2b6cb0")
segments(5,2.95,4.0,1.85,col="#c53030");segments(5,2.95,6.0,1.85,col="#c53030")
segments(8.6,2.95,7.7,1.85,col="#2b6cb0");segments(8.6,2.95,9.5,1.85,col="#2b6cb0")
text(5,0.3,expression(paste("within-nest correlation = 1 - ",lambda^2,";  estimated ",lambda,"=0.875, correlation=0.235")),cex=0.92,col="#2d3748")
dev.off()
cat("wrote fig_19_nest.svg\n")

## save headline numbers
saveRDS(list(naive=naive), "ch19_naive.rds")
cat("\n=== HEADLINE NUMBERS ===\n")
cat("naive alpha", round(D$est$alpha_naive,3), "-> cf", round(D$est$alpha_cf,3),
    "-> mixed", round(E$est_m["alpha_bar"],3), "(truth 0.30)\n")
cat("naive own-elas A", round(naive$own,3), "-> mixed", round(E$div$mixed$own_elas,3),
    "(truth", round(E$div$truth$own_elas,3),")\n")
