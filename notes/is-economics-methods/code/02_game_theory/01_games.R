# 第 2 章 运行案例：Aster 与 Boreal 的若干博弈
# 全部数字均可手算；本脚本把手算过程原样誊写并自检。
# 运行时请 cd 进本目录。

## 6.1 定价博弈：占优与纳什均衡 --------------------------------------------
# 收益矩阵 A[Aster 行动, Boreal 行动]，值为 Aster 收益；对称博弈 B = t(A)
A <- matrix(c(10, 2, 14, 6), nrow = 2, byrow = TRUE,
            dimnames = list(c("High", "Low"), c("High", "Low")))
High_dominated <- all(A["Low", ] > A["High", ])   # Low 严格支配 High
cat("High 被 Low 严格支配（Aster）:", High_dominated, " -> NE=(Low,Low)=(6,6)\n")

## 6.2 兼容性博弈：纯与混合策略均衡 ----------------------------------------
# (S1,S1)=(3,1) (S2,S2)=(1,3) mismatch=(0,0)；A 偏好 S1，B 偏好 S2
# 纯策略 NE：逐格无单方面有利偏离
# 混合：A 以 p 选 S1，B 以 q 选 S1
p <- 3/4        # 由 B 无差异 p*1 = (1-p)*3 解出
q <- 1/4        # 由 A 无差异 3q = 1-q 解出
coord_prob <- p*q + (1-p)*(1-q)
EU_A <- 3*q     # A 期望收益（均衡下两纯策略相等）
cat(sprintf("混合 NE: p=%.2f q=%.2f 协调概率=%.3f A 期望=%.2f\n",
            p, q, coord_prob, EU_A))
stopifnot(abs(coord_prob - 3/8) < 1e-12)

## 6.3 序贯博弈：逆向归纳与 SPE --------------------------------------------
# A 先动，B 观察后匹配（匹配得 1 或 3，均 > 0）；A 预见 -> 选 S1 得 3
aster_choice <- if (3 > 1) "S1" else "S2"
cat("SPE: Aster ->", aster_choice, "; Boreal 匹配 -> 收益 (3, 1)\n")

## 6.4 不完全信息：Bayesian Nash equilibrium -------------------------------
# Boreal: Normal 型(prob beta) 占优 Low；Committed 型(1-beta) 占优 High
# Aster 期望收益：面对 High 选 High=10/Low=7；面对 Low 选 High=3/Low=6
beta <- 2/5
EU_high <- beta*3 + (1 - beta)*10     # = 10 - 7 beta
EU_low  <- beta*6 + (1 - beta)*7      # = 7 - beta
aster_bne <- ifelse(EU_high >= EU_low, "High", "Low")
cat(sprintf("BNE beta=%.1f: EU_high=%.1f EU_low=%.1f -> Aster=%s (门槛 beta*=0.5)\n",
            beta, EU_high, EU_low, aster_bne))

## 6.5 重复博弈：合作门槛 --------------------------------------------------
R <- 10; T <- 14; P <- 6
delta_star <- (T - R) / (T - P)
coop  <- function(d) R / (1 - d)
cheat <- function(d) T + d * P / (1 - d)
cat(sprintf("grim 门槛 delta*=%.2f；在门槛处 coop=%.1f cheat=%.1f\n",
            delta_star, coop(delta_star), cheat(delta_star)))
stopifnot(abs(coop(delta_star) - cheat(delta_star)) < 1e-9)

## 4.1 Comparative statics checks ------------------------------------------
# Higher continuation value makes cooperation strictly easier.
deltas <- c(0.4, 0.5, 0.7, 0.9)
slack <- coop(deltas) - cheat(deltas)
cat("cooperation slack by delta:", paste(round(slack, 2), collapse=", "), "\n")
stopifnot(all(diff(slack) > 0), abs(slack[deltas == delta_star]) < 1e-9)
