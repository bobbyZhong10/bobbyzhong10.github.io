# 第 3 章 运行案例：Vireo 广告位拍卖与若干机制
# 全部数字均可手算；本脚本把手算过程原样誊写并自检。
# 运行时请 cd 进本目录。

## 6.1 效率配置与 VCG 支付 -------------------------------------------------
v   <- c(10, 6, 4)           # 每次点击估值（已排序）
ctr <- c(100, 40, 0)         # 广告位点击率（第 3 位 = 落空）
# 效率配置：广告主 i -> 广告位 i
# VCG 支付 = 他人在其缺席时的最大价值 - 他人在选定配置下的实得价值
pay1 <- (v[2]*ctr[1] + v[3]*ctr[2]) - (v[2]*ctr[2] + v[3]*ctr[3])
pay2 <- (v[1]*ctr[1] + v[3]*ctr[2]) - (v[1]*ctr[1] + v[3]*ctr[3])
VCG_rev <- pay1 + pay2
util1 <- v[1]*ctr[1] - pay1
util2 <- v[2]*ctr[2] - pay2
cat(sprintf("VCG: pay1=%d (%.1f/click) pay2=%d (%.1f/click) revenue=%d\n",
            pay1, pay1/ctr[1], pay2, pay2/ctr[2], VCG_rev))
cat(sprintf("     util1=%d util2=%d (均非负 -> 真话占优)\n", util1, util2))
stopifnot(VCG_rev == 680, util1 >= 0, util2 >= 0)

## 6.2 GSP 对照 -----------------------------------------------------------
# 天真如实报价：顶部位付第二名报价，次位付第三名报价
gsp_naive <- v[2]*ctr[1] + v[3]*ctr[2]
cat(sprintf("GSP 天真真话收入=%d；局部无嫉妒均衡复现 VCG=%d\n", gsp_naive, VCG_rev))

## 6.3 收入等价与保留价 ---------------------------------------------------
# 两广告主估值独立均匀[0,1]
n <- 2
spa <- (n - 1)/(n + 1)                 # 第二价格：次高值期望
fpa <- (n - 1)/n * n/(n + 1)           # 第一价格：压价 b(v)=(n-1)/n v，最高报价期望
cat(sprintf("收入等价: 第二价格=%.4f 第一价格=%.4f\n", spa, fpa))
stopifnot(abs(spa - fpa) < 1e-12)
# 最优保留价 r* = 1/2（解 r = (1-F)/f = 1-r）；n=2 期望收入
rev_reserve <- 5/12
cat(sprintf("最优保留价 r*=0.5，期望收入=%.4f（对比无保留价 %.4f）\n", rev_reserve, spa))

## 3.1 Akerlof 逆向选择自检 -----------------------------------------------
# 质量 theta~U[0,1]，卖家估值 theta，买家估值 1.5 theta
# 价格 p 下卖者 theta<=p，均质 p/2，买家 WTP = 1.5*(p/2) = 0.75p < p
p_grid <- seq(0.1, 1, by = 0.1)
wtp <- 1.5 * (p_grid/2)
cat("Akerlof: WTP < 价格 恒成立 ->", all(wtp < p_grid), " 市场崩解至 theta=0\n")

## 3.4 道德风险自检 -------------------------------------------------------
# 高努力成本 5，高业绩概率 0.8（低努力 0.2），有限责任 w>=0
# IC: 0.6 (wH - wL) >= 5 -> wH-wL >= 25/3；wL=0 -> wH=25/3；租 = 0.8 wH - 5
wH <- 25/3; wL <- 0
rent <- 0.8*wH - 5
cat(sprintf("道德风险: wH-wL>=%.3f, wH=%.3f, 代理人租=%.3f\n", 25/3, wH, rent))
stopifnot(abs(0.6*(wH - wL) - 5) < 1e-9)

## 4.6 Comparative statics checks ------------------------------------------
# First-price equilibrium b(v)=(n-1)v/n: shading falls as bidder count rises.
n_grid <- 2:8; v <- 1
shading <- v - ((n_grid - 1)/n_grid)*v
cat("first-price shading by bidder count:", paste(round(shading, 3), collapse=", "), "\n")
stopifnot(all(diff(shading) < 0))
