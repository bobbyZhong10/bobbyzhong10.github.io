D <- readRDS("meridian_marketplace.rds"); R <- readRDS("meridian_results.rds")
dir.create("../../assets/fig",recursive=TRUE,showWarnings=FALSE)

svg("../../assets/fig/fig_26_market_interference.svg",width=8.2,height=4.6)
par(mar=c(6,4,1,1))
bp<-barplot(R$table$estimate,names.arg=c("Individual","Market-period\ncluster","Switchback"),
            col=c("#c53030","#2b6cb0","#4a7c59"),border=NA,
            ylab="Estimated treatment contrast",ylim=c(0,0.13))
arrows(bp,R$table$estimate-1.96*R$table$se,bp,R$table$estimate+1.96*R$table$se,
       angle=90,code=3,length=0.05)
abline(h=D$truth["direct"],lty=3,lwd=2,col="#805ad5")
abline(h=D$truth["overall"],lty=2,lwd=2,col="#1a365d")
legend("topright",c("Direct effect","Global rollout effect"),
       col=c("#805ad5","#1a365d"),lty=c(3,2),lwd=2,bty="n")
dev.off()

svg("../../assets/fig/fig_26_switchback.svg",width=8.2,height=4.4)
sb <- unique(D$switchback[D$switchback$market<=6,c("market","period","d")])
sb <- sb[order(sb$market,sb$period),]
S6 <- matrix(sb$d,nrow=20)
S7 <- rbind(S6,S6[20,])  # repeat period 20 so the last block is drawn two periods wide
matplot(1:21,S7+rep((1:6-3.5)*0.03,each=21),type="s",lty=1,lwd=2,
        col=hcl.colors(6,"Dark 3"),xlab="Period",ylab="Treatment status",
        xlim=c(1,24),ylim=c(-0.15,1.15),yaxt="n")
axis(2,at=c(0,1),labels=c("Control","Treated"))
legend("right",paste("Market",1:6),col=hcl.colors(6,"Dark 3"),lty=1,lwd=2,bty="n")
dev.off()
cat(sprintf("individual %.3f | cluster %.3f | switchback %.3f\n",
            R$table$estimate[1],R$table$estimate[2],R$table$estimate[3]))
