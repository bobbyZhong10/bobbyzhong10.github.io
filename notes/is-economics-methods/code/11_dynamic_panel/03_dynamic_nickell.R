# ============================================================================
# 03_dynamic_nickell.R  --  the dynamic panel: two wrong numbers bracket the
# truth. Pooled OLS of y on its lag is biased UP, within/FE is biased DOWN
# (Nickell), and Arellano-Bond / Blundell-Bond GMM recover rho. Also: the
# Nickell bias shrinks like 1/T, and difference GMM goes weak as rho -> 1.
# ============================================================================

suppressMessages(library(plm))
N <- readRDS("nimbus.rds"); truth <- N$truth; gen_dyn <- N$gen_dyn

d <- pdata.frame(N$dyn_main, index = c("id", "t"))

## pooled OLS and within, both regressing y on lag(y)
pool <- plm(y ~ lag(y, 1), data = d, model = "pooling")
fe   <- plm(y ~ lag(y, 1), data = d, model = "within")

## Arellano-Bond difference GMM: lagged LEVELS instrument differenced equation
ab <- pgmm(y ~ lag(y, 1) | lag(y, 2:99), data = d,
           effect = "individual", model = "twosteps", transformation = "d")
## Blundell-Bond system GMM: add level equation with lagged-difference instruments
bb <- pgmm(y ~ lag(y, 1) | lag(y, 2:99), data = d,
           effect = "individual", model = "twosteps", transformation = "ld")

rho_ab <- coef(ab)["lag(y, 1)"]; rho_bb <- coef(bb)["lag(y, 1)"]
cat("=== dynamic panel: persistence rho (truth 0.70) ===\n")
cat(sprintf("pooled OLS     rho = %.4f   (biased UP: alpha in y and its lag)\n",
            coef(pool)["lag(y, 1)"]))
cat(sprintf("within / FE    rho = %.4f   (biased DOWN: Nickell)\n",
            coef(fe)["lag(y, 1)"]))
cat(sprintf("Arellano-Bond  rho = %.4f   (difference GMM)\n", rho_ab))
cat(sprintf("Blundell-Bond  rho = %.4f   (system GMM)\n", rho_bb))

## AB specification tests: Sargan overid + AR(1) present, AR(2) absent
sg  <- sargan(ab)
ar1 <- mtest(ab, order = 1); ar2 <- mtest(ab, order = 2)
cat(sprintf("\nAB tests: Sargan p = %.3f ; AR(1) p = %.4f (expect small); AR(2) p = %.3f (expect large)\n",
            sg$p.value, ar1$p.value, ar2$p.value))

## ---- Nickell bias vs 1/T: FE bias shrinks as T grows -----------------------
## FE within (fast) shows the O(1/T) attenuation; AB (capped, collapsed
## instruments) stays consistent at every T.
cat("\n=== FE (Nickell) bias of rho vs sample length T (truth 0.70) ===\n")
set.seed(303)
Ts <- c(4, 8, 16, 32)
nick <- sapply(Ts, function(TT){
  est <- replicate(40, {
    dd <- pdata.frame(gen_dyn(N = 300, T = TT, rho = 0.70), index = c("id","t"))
    fe <- coef(plm(y ~ lag(y,1), dd, model = "within"))["lag(y, 1)"]
    ab <- coef(pgmm(y ~ lag(y,1) | lag(y, 2:4), dd, effect="individual",
                    model="onestep", transformation="d", collapse=TRUE))["lag(y, 1)"]
    c(fe = fe, ab = ab)
  })
  c(fe = mean(est["fe.lag(y, 1)",]), ab = mean(est["ab.lag(y, 1)",]))
})
colnames(nick) <- Ts
nick <- rbind(nick, nickell_approx = 0.70 - (1+0.70)/(Ts-1))
print(round(nick, 4))

## ---- rho near 1: difference GMM weakens, system GMM rescues ----------------
dp <- pdata.frame(N$dyn_persist, index = c("id","t"))
ab9 <- pgmm(y ~ lag(y,1) | lag(y,2:99), dp, effect="individual",
            model="twosteps", transformation="d")
bb9 <- pgmm(y ~ lag(y,1) | lag(y,2:99), dp, effect="individual",
            model="twosteps", transformation="ld")
cat(sprintf("\n=== rho near 1 (truth 0.90): difference vs system GMM ===\n"))
cat(sprintf("difference GMM rho = %.4f (SE %.4f)  <- weak instruments\n",
            coef(ab9)["lag(y, 1)"], sqrt(diag(vcov(ab9)))["lag(y, 1)"]))
cat(sprintf("system GMM     rho = %.4f (SE %.4f)  <- rescued\n",
            coef(bb9)["lag(y, 1)"], sqrt(diag(vcov(bb9)))["lag(y, 1)"]))

saveRDS(list(pool=coef(pool)["lag(y, 1)"], fe=coef(fe)["lag(y, 1)"],
             ab=rho_ab, bb=rho_bb, sargan_p=sg$p.value,
             ar1_p=ar1$p.value, ar2_p=ar2$p.value, nick=nick,
             ab9=coef(ab9)["lag(y, 1)"], bb9=coef(bb9)["lag(y, 1)"],
             se_ab9=sqrt(diag(vcov(ab9)))["lag(y, 1)"],
             se_bb9=sqrt(diag(vcov(bb9)))["lag(y, 1)"]), "res_dynamic.rds")
