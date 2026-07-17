R <- readRDS("aurora_results.rds")
tab <- R$table

dir.create("../../assets/fig", recursive = TRUE, showWarnings = FALSE)
svg("../../assets/fig/fig_10_estimand_contrast.svg", width = 8.2, height = 4.6)
par(mar = c(6, 4, 1, 1), family = "sans")
cols <- c("#1a365d", "#2b6cb0", "#c53030")
bp <- barplot(tab$estimate, names.arg = c("Difference\nin means", "Baseline\nadjusted", "Post-treatment\ncontrol"),
              col = cols, border = NA, ylab = "Coefficient on early access",
              ylim = range(c(0, tab$estimate + 2 * tab$se, R$truth)))
arrows(bp, tab$estimate - 1.96 * tab$se, bp, tab$estimate + 1.96 * tab$se,
       angle = 90, code = 3, length = 0.05)
abline(h = R$truth["total"], col = "#4a7c59", lty = 2, lwd = 2)
abline(h = R$truth["direct"], col = "#805ad5", lty = 3, lwd = 2)
legend("topright", c("True total effect", "Structural direct effect"),
       col = c("#4a7c59", "#805ad5"), lty = c(2, 3), lwd = 2, bty = "n")
dev.off()

svg("../../assets/fig/fig_10_assignment_dag.svg", width = 8.2, height = 4.2)
plot.new(); plot.window(xlim = c(0, 10), ylim = c(0, 6))
node <- function(x, y, label, fill = "#edf2f7") {
  rect(x - 0.85, y - 0.42, x + 0.85, y + 0.42, col = fill, border = "#1a365d", lwd = 1.5)
  text(x, y, label, cex = 0.9)
}
arr <- function(x1, y1, x2, y2, col = "#4a5568") arrows(x1, y1, x2, y2, length = 0.08, col = col, lwd = 1.5)
node(1.5, 4.6, "Randomized\nearly access", "#bee3f8")
node(5.0, 4.6, "Post-treatment\nengagement", "#fed7d7")
node(8.5, 4.6, "Sales", "#c6f6d5")
node(5.0, 1.5, "Latent service\ncapability", "#faf089")
arr(2.4, 4.6, 4.1, 4.6); arr(5.9, 4.6, 7.6, 4.6)
# Route the direct-effect arrow below the mediator box instead of through it.
segments(2.25, 4.15, 3.0, 3.25, col = "#4a5568", lwd = 1.5)
segments(3.0, 3.25, 7.0, 3.25, col = "#4a5568", lwd = 1.5)
arr(7.0, 3.25, 7.75, 4.15)
arr(5.0, 1.95, 5.0, 4.1); arr(5.7, 1.8, 8.0, 4.15)
text(5, 0.55, "Engagement is a mediator and a collider descendant, not a baseline control", cex = 0.9, col = "#c53030")
dev.off()

cat(sprintf("difference %.3f | baseline adjusted %.3f | post-treatment control %.3f\n",
            tab$estimate[1], tab$estimate[2], tab$estimate[3]))
