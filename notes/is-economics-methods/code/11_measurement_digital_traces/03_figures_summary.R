dat <- readRDS("lyra_labels.rds")
out <- readRDS("lyra_results.rds")
dir.create("../../assets/fig", recursive = TRUE, showWarnings = FALSE)

svg("../../assets/fig/fig_11_label_error.svg", width = 8.2, height = 4.6)
set.seed(11); ii <- sample.int(nrow(dat), 1200)
cols <- ifelse(dat$treatment[ii] == 1, "#c5303070", "#2b6cb070")
plot(dat$true_satisfaction[ii], dat$llm_satisfaction[ii], pch = 16, cex = 0.55,
     col = cols, xlab = "Gold-standard satisfaction", ylab = "LLM-generated label")
abline(0, 1, lty = 2, lwd = 2, col = "#4a5568")
legend("topleft", c("Control", "Treatment"), pch = 16,
       col = c("#2b6cb0", "#c53030"), bty = "n")
dev.off()

svg("../../assets/fig/fig_11_logging_funnel.svg", width = 8.2, height = 4.6)
par(mar = c(6, 4, 1, 1))
bp <- barplot(out$treatment, names.arg = c("Oracle", "LLM label\nnaive", "Gold only", "PPI\nrectified"),
              col = c("#4a7c59", "#c53030", "#805ad5", "#2b6cb0"), border = NA,
              ylab = "Treatment coefficient", ylim = c(0, max(out$treatment + 2*out$se)))
arrows(bp, out$treatment - 1.96*out$se, bp, out$treatment + 1.96*out$se,
       angle = 90, code = 3, length = 0.05)
abline(h = 0.40, lty = 2, lwd = 2, col = "#1a365d")
dev.off()
cat(sprintf("naive %.3f | gold %.3f | PPI %.3f\n",
            out$treatment[2], out$treatment[3], out$treatment[4]))
