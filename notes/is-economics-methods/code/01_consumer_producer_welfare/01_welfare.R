# 第 1 章 运行案例：Lumen 的 Cobb-Douglas 消费者与福利度量
# 全部数字均可手算；本脚本把手算过程原样誊写并自检。
# 运行时请 cd 进本目录。

a  <- 2/5      # good 1 的支出份额
w  <- 100      # 预算
p2 <- 1        # outside good 价格
p1_0 <- 2      # premium content 初始价格
p1_1 <- 4      # 涨价后价格

## 6.1 UMP：Walrasian demand 与 indirect utility -----------------------------
x1 <- function(p1) a * w / p1
x2 <- function(p2) (1 - a) * w / p2
u  <- function(p1) x1(p1)^a * x2(p2)^(1 - a)   # indirect utility v(p,w)

u0 <- u(p1_0); u1 <- u(p1_1)
cat("x1_0 =", x1(p1_0), " x2 =", x2(p2), " u0 =", round(u0, 4), "\n")
cat("x1_1 =", x1(p1_1), " u1 =", round(u1, 4), "\n")
cat("支出 p1*x1 涨价前 =", p1_0 * x1(p1_0), " 涨价后 =", p1_1 * x1(p1_1), "\n")

## 6.2 EMP：支出函数与 Hicksian demand（Shephard 自检） ----------------------
e <- function(p1, uu) uu * (p1 / a)^a * (p2 / (1 - a))^(1 - a)
cat("e(p0,u0) 应为 100 =", round(e(p1_0, u0), 6), "\n")
cat("e(p1,u1) 应为 100 =", round(e(p1_1, u1), 6), "\n")
# Shephard: h1 = de/dp1，在 (p0,u0) 应等于 Walrasian x1=20
h1_shephard <- function(p1, uu) a * uu * (1 / a)^a * (p2 / (1 - a))^(1 - a) * p1^(a - 1)
cat("h1(p0,u0) =", round(h1_shephard(p1_0, u0), 4), " (应=20)\n")

## 6.3 Slutsky 分解（自身价格，初始点） --------------------------------------
total  <- -a * w / p1_0^2               # dx1/dp1
income <- -x1(p1_0) * (a / p1_0)        # -x1 * dx1/dw
subst  <- total - income                # 替代效应
cat("Slutsky: total =", total, " subst =", subst, " income =", income, "\n")

## 6.4 三个福利数字 ----------------------------------------------------------
dCS <- 40 * log(2)                      # 积分 Walrasian 需求
EV  <- w - e(p1_0, u1)                  # 旧价格计价
CV  <- e(p1_1, u0) - w                  # 新价格计价
cat("EV =", round(EV, 4), " dCS =", round(dCS, 4), " CV =", round(CV, 4), "\n")
stopifnot(EV < dCS, dCS < CV)           # 次序 EV < dCS < CV

## 6.5 生产者对偶小例 -------------------------------------------------------
# y = z1^.5 z2^.5, r=(4,1), y=10 -> c=2y sqrt(r1 r2); Shephard z_i
r1 <- 4; r2 <- 1; y <- 10
cost <- 2 * y * sqrt(r1 * r2)
z1 <- y * sqrt(r2 / r1); z2 <- y * sqrt(r1 / r2)
cat("cost =", cost, " z1 =", z1, " z2 =", z2, " 回代 =", r1 * z1 + r2 * z2, "\n")

## 6.5 不确定性小例 ---------------------------------------------------------
# u(w)=sqrt(w), 50-50 彩票 {64,144}
uw <- sqrt
Ew  <- 0.5 * 64 + 0.5 * 144
Eu  <- 0.5 * uw(64) + 0.5 * uw(144)
CE  <- Eu^2
rp  <- Ew - CE
A_Ew <- 1 / (2 * Ew)
varw <- 0.5 * (64 - Ew)^2 + 0.5 * (144 - Ew)^2
rp_approx <- 0.5 * A_Ew * varw
cat("E[w] =", Ew, " CE =", CE, " risk premium =", rp,
    " 近似 =", round(rp_approx, 4), "\n")

## 4.4 Comparative statics checks ------------------------------------------
# CRRA u(w)=w^(1-rho)/(1-rho): local risk premium rises with rho and variance.
rho_grid <- c(0.5, 1, 2)
risk_var <- 100
rp_local <- 0.5 * (rho_grid / Ew) * risk_var
cat("local risk premia by CRRA rho:", paste(round(rp_local, 3), collapse=", "), "\n")
stopifnot(all(diff(rp_local) > 0))
