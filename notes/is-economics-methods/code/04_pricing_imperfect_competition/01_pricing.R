# 第 4 章 运行案例：Cirrus 市场的定价与竞争
# 全部数字均可手算；本脚本把手算过程原样誊写并自检。
# 运行时请 cd 进本目录。

a <- 12; b <- 1; cost <- 2       # 需求 p = 12 - q，边际成本 2

## 6.1 单一价格垄断与 Lerner index ---------------------------------------
qm <- (a - cost)/(2*b); pm <- a - b*qm; profit_m <- (pm - cost)*qm
qc <- (a - cost)/b                          # 竞争产量（p = c）
CS_m <- 0.5*(a - pm)*qm; DWL <- 0.5*(qc - qm)*(pm - cost)
lerner <- (pm - cost)/pm; inv_elas <- qm/pm  # |eps| = 1*p/q, 倒数 = q/p
cat(sprintf("垄断: qm=%g pm=%g profit=%g CS=%g DWL=%g; Lerner=%.4f 1/|e|=%.4f\n",
            qm, pm, profit_m, CS_m, DWL, lerner, inv_elas))
stopifnot(abs(lerner - inv_elas) < 1e-9)

## 6.2 三级价格歧视（逆弹性规则） ----------------------------------------
q1 <- (12 - cost)/2; p1 <- 12 - q1
q2 <- (8  - cost)/2; p2 <- 8  - q2
cat(sprintf("三级PD: 细分1 p=%g |e|=%.4f markup=%.4f；细分2 p=%g |e|=%.4f markup=%.4f\n",
            p1, p1/q1, (p1-cost)/p1, p2, p2/q2, (p2-cost)/p2))

## 6.3 二级价格歧视（分层菜单 = screening） -----------------------------
thH <- 10; thL <- 6; mu <- 0.5
# 一口价基准（两类型模型内部）：单一方案至多 25
single_both <- optimize(function(q) thL*q - q^2/2, c(0, 12), maximum = TRUE)$objective  # 服务两类
single_high <- mu * (thH*thH - thH^2/2)      # 只服务高类型 q=10,T=50,概率 1/2
uniform_bench <- max(single_both, single_high)
cat(sprintf("二级PD 一口价基准: 服务两类=%.0f 只服务高类型=%.0f -> 上限=%.0f\n",
            single_both, single_high, uniform_bench))
qH <- thH                                    # 顶端不扭曲
qL <- thL - (mu/(1 - mu))*(thH - thL)        # 低端扭曲
rentH <- (thH - thL)*qL
TL <- thL*qL - qL^2/2                          # 低类型 IR 绑定
TH <- (thH*qH - qH^2/2) - rentH                # 高类型 IC 绑定
profit2 <- mu*TH + (1 - mu)*TL
profit_fb <- mu*thH^2/2 + (1 - mu)*thL^2/2     # 完全信息（一级）
uH_own <- thH*qH - qH^2/2 - TH
uH_mim <- thH*qL - qL^2/2 - TL
uL_own <- thL*qL - qL^2/2 - TL
cat(sprintf("二级PD: qL=%g qH=%g rentH=%g TL=%g TH=%g profit=%g (完全信息=%g)\n",
            qL, qH, rentH, TL, TH, profit2, profit_fb))
stopifnot(abs(uH_own - uH_mim) < 1e-9, abs(uL_own) < 1e-9, profit2 > uniform_bench)

## 6.4 寡头竞争 -----------------------------------------------------------
q_c <- (a - cost)/(3*b); Q <- 2*q_c
p_cournot <- a - b*Q; profit_each <- (p_cournot - cost)*q_c
p_bertrand <- cost                             # 无差别价格战
t <- 1; p_hotelling <- cost + t; profit_hot <- (p_hotelling - cost)*0.5
cat(sprintf("Cournot p=%.4f 每家利润=%.4f 总=%.4f；Bertrand p=%g；Hotelling p=%g 利润=%g\n",
            p_cournot, profit_each, 2*profit_each, p_bertrand, p_hotelling, profit_hot))
stopifnot(pm > p_cournot, p_cournot > p_bertrand)

## 4.8 Comparative statics checks ------------------------------------------
# Linear monopoly pass-through is 1/2; symmetric Cournot is n/(n+1).
n_grid <- 1:8
pass_through <- n_grid/(n_grid + 1)
cat("linear-demand pass-through, monopoly to Cournot n=8:",
    paste(round(pass_through, 3), collapse=", "), "\n")
stopifnot(abs(pass_through[1] - 0.5) < 1e-12, all(diff(pass_through) > 0))
