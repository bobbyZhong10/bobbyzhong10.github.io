set.seed(2604)
N <- 6000
ai_access <- rbinom(N, 1, 0.5)
workflow_training <- rbinom(N, 1, 0.5)
baseline_skill <- rnorm(N)

# A 2x2 factorial design identifies whether workflow redesign complements AI.
log_productivity <- 2 + 0.05*ai_access + 0.12*workflow_training +
  0.18*ai_access*workflow_training + 0.25*baseline_skill + rnorm(N, sd=0.55)
dat <- data.frame(ai_access, workflow_training, baseline_skill, log_productivity)
fit <- lm(log_productivity ~ ai_access*workflow_training + baseline_skill, dat)
co <- coef(summary(fit))
effects <- c(
  AI_without_training=coef(fit)["ai_access"],
  AI_with_training=coef(fit)["ai_access"]+coef(fit)["ai_access:workflow_training"],
  complementarity=coef(fit)["ai_access:workflow_training"]
)
saveRDS(list(data=dat,fit=fit,effects=effects),"res26_complements.rds")

dir.create("../../assets/fig",recursive=TRUE,showWarnings=FALSE)
means <- with(dat,tapply(log_productivity,list(workflow_training,ai_access),mean))
svg("../../assets/fig/fig_26_complements.svg",width=8.2,height=4.6)
matplot(c(0,1),t(means),type="b",pch=16,lty=1,lwd=2,
        col=c("#2b6cb0","#c53030"),xaxt="n",xlab="AI access",
        ylab="Mean log productivity")
axis(1,at=c(0,1),labels=c("No access","Access"))
legend("topleft",c("No workflow training","Workflow training"),
       col=c("#2b6cb0","#c53030"),pch=16,lty=1,lwd=2,bty="n")
dev.off()
cat(sprintf("AI effect without training %.3f | with training %.3f | interaction %.3f (SE %.3f)\n",
            effects[1],effects[2],effects[3],co["ai_access:workflow_training","Std. Error"]))
