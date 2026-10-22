---
title: "Instrumental Variables & System Estimation"
subtitle: "Instruments, 2SLS and Systems of Equations"
seriesline: "Foundations of Information Systems Economics · Chapter 9"
series: "Foundations of Information Systems Economics: Empirical Methods for Information Systems Economics"
sidebartitle: "Foundations of Information Systems Economics"
sidebarsub: "Theory · Methods · Applications"
chapterlabel: "Chapter 9 · Instrumental Variables & System Estimation"
---

## Introduction

Bramble is a fictional two-sided platform: households in a community use it to order food from home bakeries, small farms, and independent kitchens, and the platform operates in 1,500 local markets. It faces the old problem of every two-sided platform: more sellers draw buyers, and more buyers draw sellers. The growth team cares about only one half of it: if the platform recruits 100 more sellers, how many more active buyers will they bring? That number decides how much the platform should spend on onboarding subsidies for sellers. Section 1 shows that regressing the number of buyers on the number of sellers gives an answer almost twice the truth.

The problem is simultaneity. The buyer equation and the seller equation hold at the same time, so every pair of seller and buyer counts observed in a market is the intersection of the two equations, and any local shock that pushes both sides moves the two together. This chapter shows how to pull a single equation out of those intersections, and whether several equations can be estimated more precisely together. The first thread is single-equation IV (instrumental variables) and 2SLS: the onboarding fee, which appears only in the seller equation, moves the number of sellers without directly affecting households' orders, so it serves as an instrument in the buyer equation, moving the market equilibrium along that equation. The second thread is system estimation: when the errors of several equations are correlated, SUR (seemingly unrelated regressions) and 3SLS (three-stage least squares) borrow information from the other equations to gain efficiency, and the price is that one misspecified equation drags the others down with it. The third thread runs through the whole chapter: when the instruments are weak, conventional 2SLS inference breaks down, while the Anderson-Rubin set stays honest.

The chapter proceeds as follows. Section 2 writes down Bramble's two structural equations, the reduced form, and the estimands of the five data blocks. Section 3, the most detailed in the chapter, takes up in turn identification in simultaneous equations, the consistency and asymptotic normality of 2SLS, tests of endogeneity and overidentification, weak instruments, SUR and feasible generalized least squares, system 2SLS and 3SLS, and finally the trade-off between robustness and efficiency. Section 4 covers computation, diagnostics, and software, Section 5 examines three anchoring papers, Section 6 runs a full walkthrough on Bramble, and Section 7 discusses failure modes. The chapter uses the sandwich and two-step estimation of Chapter 5, the GMM and J test of Chapter 7, and the robust inference of Chapter 8; the random effects and Hausman-Taylor estimators of Chapter 10 return to this chapter's feasible generalized least squares and 2SLS.

## 1 A Cross-Side Effect That Nearly Doubles

::: {.case}
Bramble's block S1 is one month of data on 1,500 local markets, with buyers counted in thousands and sellers in hundreds. Three policies vary across markets: an ad budget that headquarters reallocates at random (mean USD 40 thousand, standard deviation USD 10 thousand), a first-order coupon (mean USD 6, standard deviation USD 1.5), and the sellers' monthly onboarding fee (mean USD 30, standard deviation USD 7). City programs pair deeper coupons with fee reductions, and the correlation between the two is $-0.7$.

Because this is a simulation, we know the answer: holding advertising, the coupon, and market size fixed, 100 more sellers bring 1,200 more buyers, so the coefficient on the number of sellers in the buyer equation is 1.2. In the other direction, 1,000 more buyers bring 25 more sellers, and each USD 1 cut in the onboarding fee brings 10 more sellers.
:::

The analyst starts with the most direct approach: an OLS regression of the number of buyers on the number of sellers, controlling for the ad budget, the coupon, and market size. The coefficient is 2.27 (robust standard error 0.015), with a 95% interval of $[2.239, 2.297]$, which reads as "100 more sellers bring more than two thousand buyers," almost twice the truth. Going by this number, the platform would set its onboarding subsidies too high. This is not sampling error: on one million simulated markets, the same regression gives 2.264, 88.6% above the truth, a number that can be taken as its probability limit.

The reason is that the same local shock pushes both sides. When a market catches a wave of local demand, say because a new school opens nearby, buyers increase and sellers follow; in the DGP, the correlation between the errors of the two equations is 0.8. The number of sellers is therefore positively correlated with the error of the buyer equation, and OLS credits the sellers with the shock-driven comovement as well.

Change the estimator and the bias disappears. The onboarding fee appears only in the seller equation: it changes sellers' cost of joining but does not directly affect whether households order. With the fee as the instrument for the number of sellers, 2SLS gives 1.14 (robust standard error 0.102), and the interval $[0.941, 1.339]$ covers 1.2. New questions follow. The standard error is about seven times that of OLS: can the estimate be made more precise? Estimating the two equations together by 3SLS gives 1.212 (standard error 0.086), a standard error about 15% smaller. Where does this precision come from, and what does it cost? In another set of markets, a single pilot moved the onboarding fee by just a few dozen cents. Can the 2SLS interval still be trusted there? These are the three questions this chapter answers.

## 2 The Economic Model and the Estimand

### 2.1 The Buyer Equation, the Seller Equation, and the Reduced Form

Write the number of active buyers in market $i$ (thousands) as $\mathrm{buyers}_i$ and the number of active sellers (hundreds) as $\mathrm{sellers}_i$. The two structural equations are

$$\begin{aligned}
\mathrm{buyers}_i &= 14.4 + 1.2\,\mathrm{sellers}_i + 0.10\,\mathrm{ad}_i \\
&\quad + 0.40\,\mathrm{coupon}_i + 2.0\,\mathrm{lpop}_i + \varepsilon_{iB}, \\
\mathrm{sellers}_i &= 9.0 + 0.25\,\mathrm{buyers}_i - 0.10\,\mathrm{fee}_i \\
&\quad + 0.6\,\mathrm{lpop}_i + \varepsilon_{iS},
\end{aligned}$$

Here $\mathrm{ad}_i$ is the ad budget (thousands of USD), $\mathrm{coupon}_i$ and $\mathrm{fee}_i$ are the coupon and the onboarding fee (USD), and $\mathrm{lpop}_i$ is standardized log population. The first is the buyer equation, which describes how households respond to the number of sellers: more sellers mean more dishes to choose from and shorter deliveries. The second is the seller equation, which describes how sellers respond to the number of buyers. Every USD 1,000 of advertising brings 100 buyers, that is, USD 10 per buyer; each USD 1 of extra coupon depth brings 400 more buyers. The local shocks $(\varepsilon_{iB}, \varepsilon_{iS})$ have standard deviations of 2.5 and 0.8 and a correlation of 0.8, and they are unrelated to the three policies and to market size.

The coefficients of the structural equations are causal quantities. The 1.2 answers this question: if the platform raised the number of sellers in a market exogenously by 100 (by filling in with its own first-party sellers, say), holding advertising, the coupon, and population fixed, how many more buyers would there be? It is not the total change in the market equilibrium. The two equations amplify each other: more sellers bring buyers, and more buyers in turn bring sellers, with a feedback multiplier of $1/(1 - 1.2 \times 0.25) = 1.43$. A USD 1 cut in the onboarding fee adds 10 sellers directly through the seller equation; in equilibrium there are 14.3 more sellers and 171 more buyers.

Write the system in general form. Let $\mathbf{y}_i$ collect the $M$ endogenous variables ($M = 2$ in Bramble) and $\mathbf{z}_i$ all $r$ exogenous variables (here the constant, advertising, the coupon, the onboarding fee, and population, so $r = 5$). Wooldridge writes $G$ for the number of equations; in this book $G$ denotes a cohort and $J$ the $J$ statistic, so we write the number of equations as $M$.

$$\mathbf{y}_i'\Gamma + \mathbf{z}_i'\Delta + \boldsymbol\varepsilon_i' = 0,$$

In $\Gamma$ ($M \times M$), column $m$ holds the coefficients of the endogenous variables in equation $m$, with the coefficient of the left-hand variable normalized to $-1$; $\Delta$ ($r \times M$) holds the coefficients of the exogenous variables. Bramble's two matrices are as follows. In both, the columns are the buyer equation and the seller equation, in that order; the rows of $\Gamma$ are the numbers of buyers and sellers, and the rows of $\Delta$ are the constant, advertising, the coupon, the onboarding fee, and population:

$$\Gamma = \begin{pmatrix} -1 & 0.25 \\ 1.2 & -1 \end{pmatrix}, \qquad \Delta = \begin{pmatrix} 14.4 & 9.0 \\ 0.10 & 0 \\ 0.40 & 0 \\ 0 & -0.10 \\ 2.0 & 0.6 \end{pmatrix}.$$

Here $\Gamma$ is the matrix of structural coefficients on the endogenous variables and has nothing to do with the $\Gamma_0$ of two-step estimation in Chapters 5 and 6. Likewise, $\Omega = \mathbb{E}[\boldsymbol\varepsilon_i\boldsymbol\varepsilon_i']$ is the variance matrix of the equation errors, not the HAC long-run variance $\hat\Omega^{HAC}$ of Chapter 8. When $\Gamma$ is invertible, solving the system gives the reduced form (Wooldridge 2010 equation (9.14)):

$$\begin{aligned}
\mathbf{y}_i' &= \mathbf{z}_i'\Pi + \mathbf{v}_i', \\
\Pi &= -\Delta\Gamma^{-1}, \qquad \mathbf{v}_i' = -\boldsymbol\varepsilon_i'\Gamma^{-1}.
\end{aligned}$$

$\Pi$ ($r \times M$), like the $\pi$ of Chapter 7, is a reduced-form parameter: it describes how far each exogenous variable pushes each endogenous variable in equilibrium, and equation-by-equation OLS estimates it consistently. The reduced-form error $\mathbf{v}_i$ ($M \times 1$) is a linear combination of the structural errors, with variance $\Omega_{\mathbf{v}} = (\Gamma^{-1})'\Omega\,\Gamma^{-1}$. The reduced-form error of the number of sellers contains $\varepsilon_{iB}$; in the DGP the two have a correlation of 0.935, and this is where the OLS bias comes from. In the $\Pi$ computed from the two matrices above, the coefficients of the onboarding fee on the numbers of sellers and buyers are $-0.143$ and $-0.171$ (hundreds and thousands per dollar), and their ratio is exactly 1.2: the equilibrium movement driven by the fee happens entirely along the buyer equation.

The estimands of this chapter are the coefficients of the structural equations: the 1.2 in block S1's buyer equation and the 0.25 in its seller equation. They cannot be read directly off the reduced form; they rest on exclusion restrictions: the onboarding fee is not in the buyer equation, and advertising and the coupon are not in the seller equation. The buyer equation excludes only one exogenous variable and has one endogenous variable on the right-hand side, so it is exactly identified; the seller equation excludes two, one more than it needs, so it is overidentified. Here the coefficients are the same in every market. If the cross-side effect differed across markets, 2SLS would estimate a weighted average only under certain conditions; Section 3.2 spells out those conditions and the weights.

### 2.2 The Other Four Data Blocks

Block S1m uses the same markets and policies as S1, but the true seller equation has an extra term: the coefficient on the number of buyers varies with the coupon, changing by $-0.0045$ for every USD 1 by which the coupon is deeper than USD 6, so the seller equation contains an interaction between the coupon and the number of buyers. Households drawn in by coupons order less, so each USD 1 of extra coupon depth makes the cross-side response smaller by 1.8%, and a coupon one standard deviation deeper makes it smaller by 2.7%. The analyst still estimates the linear seller equation of S1. The omitted term carries the coupon into the error of the seller equation, so the coupon is no longer a valid instrument there; the buyer equation, by contrast, is entirely correct.

Block S2 holds two other outcomes for the same markets: the average basket value (USD) and the number of orders per active buyer per month. The platform varies the free-delivery minimum (mean USD 35, standard deviation USD 5); a higher minimum raises the basket value and lowers the number of orders:

$$\begin{aligned}
\mathrm{basket}_i &= 18 + 0.50\,\mathrm{shipmin}_i + 25\,\mathrm{premium}_i + \varepsilon_{i1}, \\
\mathrm{orders}_i &= 3.6 - 0.04\,\mathrm{shipmin}_i + 0.5\,\mathrm{loyalty}_i + \varepsilon_{i2},
\end{aligned}$$

Here $\mathrm{premium}_i$ is the share of sellers with a premium badge, and $\mathrm{loyalty}_i$ indicates whether the market has a loyalty program (about 41% of markets do). The two equations are not simultaneous, and every right-hand variable is exogenous; the equations are linked only through their errors. A local income shock raises both the basket value and the number of orders; the error standard deviations are USD 4 and 0.45 orders, and the correlation is 0.6. This is the SUR setting. The estimands are the six coefficients and one cross-equation combination: monthly revenue per buyer is the basket value times the number of orders, and the first-order effect of raising the minimum by USD 1 is $\mathbb{E}[\mathrm{orders}] \times 0.50 + \mathbb{E}[\mathrm{basket}] \times (-0.04) = -0.518$ dollars, where $\mathbb{E}[\mathrm{basket}] = 43$ dollars and $\mathbb{E}[\mathrm{orders}] = 2.40$.

Block S3 is the weak-instrument version: the structural equations and errors are the same as in S1, the coupon and the onboarding fee are independent, and the fee varies only in a pilot, around USD 30 with a standard deviation of USD 0.349, that is, 34.9 cents. Its push on the number of sellers is drowned out by the local shocks. The concentration parameter $\mu^2$ of Chapter 7 measures this: the sample size times the variance of the number of sellers explained by the onboarding fee (after partialling out the other exogenous variables of the equation), divided by the variance of the reduced-form error of the number of sellers. In S3, $\mu^2 = 1$; with a single instrument, the first-stage F statistic has a mean of about $1 + \mu^2$.

Block S4 is one quarter of data on another set of markets, and its buyer equation is overidentified. Three seller-side programs are randomized separately: the onboarding fee, payout days (uniform on 1 to 7 days), and a referral bonus (uniform on USD 20 to 80). In the seller equation, each extra day of payout delay costs 30 sellers, and each extra USD 1 of referral bonus brings 2 sellers. In S4 the referral bonus also has a side effect: referred sellers bring their existing customers onto the platform, so every extra USD 10 of bonus adds 500 buyers directly. The bonus therefore violates the exclusion restriction in the buyer equation; the control block S4v has no such side effect. In both blocks the errors vary with the referral bonus, and the largest conditional variance is about 32 times the smallest.

The five data blocks ask the same thing: under what conditions a single structural equation can be pulled out of equilibrium data, and when estimating several equations together is worth it. Section 3 begins by turning identification into theorems.

## 3 Identification and Large-Sample Theory

This section follows the format of Chapters 5 to 7. It first uses Bramble's numbers or figures to show what is to be proved, then states numbered assumptions and theorems, followed by proof sketches of two to six steps; each step names the assumptions it uses, and regularity conditions are only cited to their source. After each theorem, an "In words" paragraph answers four questions: under which assumptions it estimates which quantity, what it estimates when the specification is wrong, who sets the weights, and what Bramble's numbers say. Assumption labels follow the format of Wooldridge (2010), with SEM.1 and WI.1 coined in this chapter; the theorem numbers are this chapter's own. Section 3.1 covers identification in simultaneous equations, Sections 3.2 to 3.5 cover single-equation 2SLS, its tests, and weak instruments, and Sections 3.6 to 3.8 cover system estimation.

### 3.1 Simultaneity and Identification

Start with a figure. Partial out advertising, the coupon, and market size from both the number of buyers and the number of sellers, and plot one point for every third market (figure below). The cloud lies along a steep line, and the OLS slope is 2.27. The true buyer equation has slope 1.2 and is much flatter; in this plane, the seller equation has slope $1/0.25 = 4$ (not drawn). The cloud shows how the local shocks scatter the intersections of the two equations: shocks on the buyer side move the intersection along the seller equation, shocks on the seller side move it along the buyer equation, and OLS fits a mixture of the two movements, not either equation. Yet when the markets are grouped into five bins by the onboarding fee (with the other variables likewise partialled out), the five group means fall almost exactly on the buyer equation: the fee moves only the seller equation, so the intersection can only slide along the buyer equation.

![Simultaneity in block S1. The horizontal axis is the number of sellers (hundreds) and the vertical axis the number of buyers (thousands), both with the ad budget, the coupon, and market size partialled out; the gray points show every third market, 500 in all. The blue dashed line is the OLS fit through the cloud (slope 2.27), the dark gray solid line is the true buyer equation (slope 1.2), and the orange solid line is 2SLS with the onboarding fee as the instrument (slope 1.14); the orange circles are the group means by quintile of the partialled-out fee.](assets/fig/fig_09_simultaneity.svg)

This is the question that Working (1927) raised: whether the "statistical demand curve" traced out by a scatter of prices and quantities is really a demand curve. As Imbens (2014) §3.5 restates it, when the shocks to the two equations are uncorrelated, the slope of the regression of equilibrium quantity on price is a weighted average of the supply and demand slopes, with weights from the residual variances of the two equations; when the shocks are correlated an extra term appears, so the statistical demand curve is in general not the demand function. Bramble's two shocks are correlated, so Bramble falls in the second case. For the textbook derivation of how excluded shifters identify an equation, see Wooldridge (2010) §9.2.1, pp. 242-244.

::: {.assumption}
**SEM.1 (Wooldridge 2010 equation (9.2), p. 242; for the invertibility of $\Gamma$, see §9.2.2, p. 246)** $\mathbb{E}[\mathbf{z}_i\boldsymbol\varepsilon_i'] = 0$, $\mathbb{E}[\mathbf{z}_i\mathbf{z}_i']$ is nonsingular, and $\Gamma$ is invertible. In each equation the coefficient of the left-hand variable is normalized to $-1$, and the remaining prior information consists of exclusion restrictions: certain exogenous or endogenous variables do not appear in certain equations.
:::

The exogeneity in SEM.1 is stronger than what a single equation requires: the onboarding fee must be uncorrelated not only with the error of the seller equation but also with that of the buyer equation. The latter is the other half of what the instrument needs, besides "the fee does not directly affect households' orders." Exclusion restrictions are economic assumptions; the data cannot choose them for the researcher. The identification question is this: the reduced form $\Pi$ can always be estimated, but can a given structural equation be recovered from it? Take equation $m$ and let $\Pi_{(m)}$ be the submatrix of $\Pi$ formed by the exogenous variables excluded from equation $m$ (rows) and the endogenous variables on the right-hand side of that equation (columns).

::: {.theorem}
**Theorem 1 (identification in simultaneous equations; Wooldridge 2010 equations (9.9) and (9.10), p. 244; Theorem 9.1, p. 245; Theorems 9.2 and 9.3, pp. 248-249)** Suppose SEM.1 holds.

\(i) With exclusion restrictions only, equation $m$ is identified if and only if $\Pi_{(m)}$ has full column rank; this is called the rank condition.

\(ii) A necessary condition for it is the order condition: the number of exogenous variables excluded from the equation is at least the number of endogenous variables on its right-hand side (Theorem 9.1). When the two are equal the equation is exactly identified, and the excess of the former over the latter is the number of overidentifying restrictions.

\(iii) More generally, if the prior information on equation 1 is a set of linear restrictions (that is, $R_1$ times column 1 of $[\Gamma', \Delta']'$ equals zero), then it is identified if and only if $\mathrm{rank}\big(R_1[\Gamma', \Delta']'\big) = M - 1$ (Theorem 9.2); a necessary condition is that the number of restrictions be at least $M - 1$ (Theorem 9.3).
:::

The proof sketch has three steps. Step one (SEM.1): $\mathbb{E}[\mathbf{z}_i\mathbf{z}_i']$ is nonsingular, so equation-by-equation OLS estimates $\Pi$ and $\Omega_{\mathbf{v}}$ consistently, and they are therefore identified (after equation (9.14)). Step two (SEM.1 and the exclusion restrictions): the right-hand side of equation $m$ consists of its own exogenous variables and a few endogenous variables, whose projection coefficients on $\mathbf{z}_i$ are the corresponding columns of $\Pi$. Single-equation identification requires $\mathbb{E}[\mathbf{z}_i\mathbf{x}_i']$ to have full column rank, where $\mathbf{x}_i$ collects all right-hand variables of the equation (Section 3.2 calls this 2SLS.2(b)). Each of the equation's own exogenous variables takes up one unit vector, so the rest of the requirement falls entirely on the excluded rows: $\Pi_{(m)}$ must have full column rank (equations (9.9) and (9.10)), which is (i). The matrix $\Pi_{(m)}$ has as many rows as there are excluded variables and as many columns as there are right-hand endogenous variables, and with too few rows it cannot have full column rank, which is (ii). Step three (general restrictions): postmultiplying the whole system by any nonsingular $M \times M$ matrix $C$ leaves the reduced form unchanged (equation (9.15)), so identifying the equation requires the restrictions to rule out every nontrivial $C$: $R_1[\Gamma', \Delta']'$ times column 1 of $C$ must be zero, and this must allow that column to be only a multiple of the first unit vector, that is, the null space of this matrix must be one-dimensional. The normalization then fixes the multiple (equations (9.18) and (9.19)), which is (iii).

In words: under SEM.1 and the exclusion restrictions, the structural coefficients are functions of the reduced form. The buyer equation is exactly identified, and its coefficient is simply the ratio of the onboarding fee's coefficients in the two reduced forms, $-0.171$ to $-0.143$, which is exactly 1.2. If the exclusion restriction is wrong, say because sellers in high-fee cities also charge higher prices and households therefore order less, the ratio still exists, but it equals 1.2 plus "the fee's direct effect on buyers divided by its equilibrium effect on sellers," and what is identified is a number with no economic meaning. An exactly identified equation has no spare restrictions to test, so the theorem itself raises no alarm. Who sets the weights: under exact identification there are no weights to choose; under overidentification, the same coefficient can be computed once from each excluded variable, and how to combine these is left to the 2SLS of Section 3.2 and the 3SLS of Section 3.7. Bramble's seller equation excludes advertising and the coupon; its $\Pi_{(m)} = (0.143, 0.571)'$ has rank 1, with one restriction to spare. The buyer equation's $\Pi_{(m)} = -0.143$ also has rank 1. Here (iii) reduces to (i): the buyer equation's only prior restriction is that the fee does not appear, $R_1$ picks the fee's row out of $[\Gamma', \Delta']'$, and $R_1[\Gamma', \Delta']' = (0, -0.10)$, with rank $1 = M - 1$. If the fee did not affect sellers at all (set its coefficient in the seller equation to zero), the buyer equation's order condition would still hold, but $\Pi_{(m)}$ would become zero, with rank 0, and the equation would not be identified: the order condition only counts, while the rank condition checks whether the excluded variables actually move the endogenous variables.

Beyond exclusion restrictions there are two other kinds of identifying information (Wooldridge 2010 §9.4). One is cross-equation restrictions: if a coefficient is equal in two equations, the coefficient identified in one equation can be substituted into the other, and an equation that was not identified may become identified as a result. In estimation, 3SLS writes the restriction directly into the system (equations (9.31) to (9.34), pp. 256-257). The other is covariance restrictions: if the errors of two equations are assumed uncorrelated, that is, if the off-diagonal elements of $\Omega$ are zero, the residuals of an identified equation can serve as instruments for the other equation (equations (9.35) to (9.37), pp. 257-258). In a two-sided market like Bramble, where local shocks act on both sides at once, this diagonal assumption is hard to defend; Wooldridge also points out that in truly simultaneous microeconomic systems it is hard to find reasons for the errors to be uncorrelated. In linear simultaneous equations, the identifying information discussed in this chapter comes from three sources: exclusion restrictions, cross-equation restrictions, and covariance restrictions. Each is an economic assumption.

### 3.2 Consistency of 2SLS

The seller equation has two excluded exogenous variables, advertising and the coupon, each of which could identify the coefficient on the number of buyers on its own. 2SLS first regresses the number of buyers on all exogenous variables; this is the first stage, and its coefficients on advertising and the coupon are 0.119 and 0.657. The fitted values then replace the number of buyers in a second-step regression. These fitted values contain only the variation driven by the exogenous variables and are unrelated to the local shocks. The result is 0.219 (robust standard error 0.018), against a truth of 0.25.

Write the structural equation as $y_i = \mathbf{x}_i'\beta + \varepsilon_i$, where $\mathbf{x}_i$ is the $k \times 1$ vector of regressors (including the constant and the exogenous variables of the equation) and $\mathbf{z}_i$ is the $r \times 1$ vector of instruments (including the exogenous variables in $\mathbf{x}_i$). The single-equation first stage is written $\mathbf{x}_i' = \mathbf{z}_i'\Pi_{\mathrm{fs}} + \mathbf{v}_i'$, where $\Pi_{\mathrm{fs}}$ ($r \times k$) holds the linear projection coefficients and the first-stage error $\mathbf{v}_i$ is $k \times 1$; the exogenous part of $\mathbf{x}_i$ is also in $\mathbf{z}_i$, so the corresponding elements of $\mathbf{v}_i$ are zero. $\Pi_{\mathrm{fs}}$ and the reduced form $\Pi$ are not the same matrix: the former projects the regressors of one equation on all instruments, the latter all endogenous variables on all exogenous variables. Stack the observations so that $\mathbf{y}$, $\mathbf{X}$, and $\mathbf{Z}$ each have $n$ rows, and let $\mathbf{P}_{\mathbf{Z}} = \mathbf{Z}(\mathbf{Z}'\mathbf{Z})^{-1}\mathbf{Z}'$; the first-stage fitted values are $\hat{\mathbf{X}} = \mathbf{P}_{\mathbf{Z}}\mathbf{X}$, so that

$$\hat\beta_{2SLS} = (\hat{\mathbf{X}}'\mathbf{X})^{-1}\hat{\mathbf{X}}'\mathbf{y} = (\hat{\mathbf{X}}'\hat{\mathbf{X}})^{-1}\hat{\mathbf{X}}'\mathbf{y}.$$

The second equality holds because $\mathbf{P}_{\mathbf{Z}}$ is symmetric and idempotent. When $r = k$ this is the IV estimator $\hat\beta_{IV} = (\mathbf{Z}'\mathbf{X})^{-1}\mathbf{Z}'\mathbf{y}$.

::: {.assumption}
**2SLS.1 (Wooldridge 2010 §5.2.1, p. 99)** $\mathbb{E}[\mathbf{z}_i\varepsilon_i] = 0$.

**2SLS.2 (p. 99)** (a) $\mathrm{rank}\,\mathbb{E}[\mathbf{z}_i\mathbf{z}_i'] = r$; (b) $\mathrm{rank}\,\mathbb{E}[\mathbf{z}_i\mathbf{x}_i'] = k$.
:::

2SLS.1 is exogeneity plus the exclusion restriction: the excluded instruments are uncorrelated with the structural error, and it is only through $\mathbf{x}_i$ that they affect $y_i$; otherwise their direct effect remains in the error. 2SLS.2(b) is the single-equation rank condition, and Theorem 1 (i) is what it looks like in simultaneous equations; with a single endogenous regressor, it amounts to at least one excluded instrument having a nonzero first-stage coefficient, which can be tested directly (Wooldridge 2010 equation (5.20), p. 98). A necessary condition for it is $r \geq k$.

::: {.theorem}
**Theorem 2 (consistency of 2SLS; Wooldridge 2010 Theorem 5.1, p. 100; identification equation (5.21))** Under 2SLS.1 and 2SLS.2 (random sampling),

$$\begin{aligned}
\beta &= \big(\mathbb{E}[\Pi_{\mathrm{fs}}'\mathbf{z}_i\mathbf{x}_i']\big)^{-1}\mathbb{E}[\Pi_{\mathrm{fs}}'\mathbf{z}_iy_i], \\
\Pi_{\mathrm{fs}} &= \big(\mathbb{E}[\mathbf{z}_i\mathbf{z}_i']\big)^{-1}\mathbb{E}[\mathbf{z}_i\mathbf{x}_i'],
\end{aligned}$$

and $\hat\beta_{2SLS} \xrightarrow{p} \beta$.
:::

The proof sketch has three steps. Step one (2SLS.2(a)): $\Pi_{\mathrm{fs}}$ is the coefficient of the linear projection of $\mathbf{x}_i$ on $\mathbf{z}_i$, and the sample first stage estimates it consistently. Step two (2SLS.1): premultiply both sides of the structural equation by $\Pi_{\mathrm{fs}}'\mathbf{z}_i$ and take expectations; the error term vanishes because $\mathbb{E}[\mathbf{z}_i\varepsilon_i] = 0$, which gives equation (5.21). Step three (2SLS.2(b)): $\mathbb{E}[\Pi_{\mathrm{fs}}'\mathbf{z}_i\mathbf{x}_i'] = \mathbb{E}[\mathbf{x}_i\mathbf{z}_i']\big(\mathbb{E}[\mathbf{z}_i\mathbf{z}_i']\big)^{-1}\mathbb{E}[\mathbf{z}_i\mathbf{x}_i']$ is invertible if and only if $\mathbb{E}[\mathbf{z}_i\mathbf{x}_i']$ has rank $k$; the law of large numbers and Slutsky's theorem take the sample counterparts to their population values.

In words: under 2SLS.1 and 2SLS.2, 2SLS estimates the coefficient $\beta$ of the structural equation, provided it is the same in every market. With a single instrument it is the ratio of the reduced-form coefficient to the first-stage coefficient; the $\mathrm{Cov}(z, y)/\mathrm{Cov}(z, x)$ of Section 2.1 of Chapter 7 is this same ratio. If 2SLS.1 fails, for example because the instrument has a direct effect, the ratio converges to $\beta$ plus "the covariance of the instrument with the error divided by the covariance of the instrument with the regressor" (Wooldridge 2010 equation (5.36), p. 108): the denominator is the strength of the first stage, so the weaker the instrument, the larger the bias from the same small violation of the exclusion restriction; it can even exceed the bias of OLS. Who sets the weights: with several instruments, 2SLS combines them into one through the first-stage fitted values, using the first-stage coefficients as weights, and these are set by the data rather than by the researcher. In the seller equation the first-stage coefficients of advertising and the coupon are 0.119 and 0.657 (per thousand dollars and per dollar), and they decide how much each instrument contributes to the fitted values. When the coefficient varies across markets, the weighted-average reading needs two additional conditions: the instrument pushes the regressor in the same direction in every market, and the market-level effects are unrelated to the instrument. When both hold, 2SLS estimates a weighted average of the market-level effects with nonnegative weights that depend on how far the instrument moves each market's regressor, and markets the instrument cannot move do not enter the average; when they fail, weights can be negative, and 2SLS need not lie within the range of the market-level effects. With several instruments, the same-direction condition must hold jointly over the values of all the instruments; it is not enough for each instrument to push in a consistent direction on its own (Section 4.1 of Chapter 18). With a common coefficient none of these conditions is needed. What it estimates under heterogeneous treatment effects is the subject of Chapter 18 on LATE (local average treatment effect), where monotonicity and independence are these same two conditions. In Bramble the effect is the same in every market: the buyer equation's 2SLS aims at 1.2, with a sample estimate of 1.140, and the seller equation's aims at 0.25, with a sample estimate of 0.219 and an interval $[0.184, 0.254]$ that covers the truth.

### 3.3 Asymptotic Normality, Relative Efficiency, and the Robust Variance

In 1,000 Monte Carlo replications of block S1, the buyer-equation 2SLS estimates average 1.195 with a standard deviation of 0.096, and the nominal 95% interval built from the robust standard error has a coverage rate of 0.942. In a single sample, the standard error from the homoskedastic formula, 0.1010, is almost the same as the robust one, 0.1017, because the errors of S1 are homoskedastic. The theorem below shows where these standard errors come from.

::: {.assumption}
**2SLS.3 (Wooldridge 2010 §5.2.2, p. 101)** $\mathbb{E}[\varepsilon_i^2\mathbf{z}_i\mathbf{z}_i'] = \sigma^2\,\mathbb{E}[\mathbf{z}_i\mathbf{z}_i']$, where $\sigma^2 = \mathbb{E}[\varepsilon_i^2]$.
:::

2SLS.3 is the homoskedasticity assumption of 2SLS, and $\mathbb{E}[\varepsilon_i^2 \mid \mathbf{z}_i] = \sigma^2$ is a sufficient condition for it (equation (5.23)). When the regressors are endogenous, an assumption on the conditional variance given $\mathbf{x}_i$ makes no sense, so the conditioning is on the instruments.

::: {.theorem}
**Theorem 3 (asymptotic normality and efficiency of 2SLS; Wooldridge 2010 Theorem 5.2 and equations (5.24) to (5.27), pp. 101-102; Theorem 5.3, p. 103; robust variance equation (5.34), p. 106)**

\(i) Under 2SLS.1, 2SLS.2, and the usual moment conditions, $\sqrt n(\hat\beta_{2SLS} - \beta)$ is asymptotically normal with mean zero and variance

$$\begin{aligned}
&\big(\mathbb{E}[\Pi_{\mathrm{fs}}'\mathbf{z}_i\mathbf{z}_i'\Pi_{\mathrm{fs}}]\big)^{-1} \\
&\quad \times \mathbb{E}[\varepsilon_i^2\,\Pi_{\mathrm{fs}}'\mathbf{z}_i\mathbf{z}_i'\Pi_{\mathrm{fs}}] \\
&\quad \times \big(\mathbb{E}[\Pi_{\mathrm{fs}}'\mathbf{z}_i\mathbf{z}_i'\Pi_{\mathrm{fs}}]\big)^{-1},
\end{aligned}$$

estimated by $(\hat{\mathbf{X}}'\hat{\mathbf{X}})^{-1}\big(\sum_i\hat\varepsilon_i^2\hat{\mathbf{x}}_i\hat{\mathbf{x}}_i'\big)(\hat{\mathbf{X}}'\hat{\mathbf{X}})^{-1}$ (equation (5.34)), where $\hat{\mathbf{x}}_i'$ is the row of $\hat{\mathbf{X}}$ for observation $i$ and $\hat\varepsilon_i = y_i - \mathbf{x}_i'\hat\beta_{2SLS}$.

\(ii) If 2SLS.3 also holds, the variance collapses to $\sigma^2\big(\mathbb{E}[\Pi_{\mathrm{fs}}'\mathbf{z}_i\mathbf{z}_i'\Pi_{\mathrm{fs}}]\big)^{-1}$ (equation (5.24)), estimated by $\hat\sigma^2(\hat{\mathbf{X}}'\hat{\mathbf{X}})^{-1}$ with $\hat\sigma^2 = (n - k)^{-1}\sum_i\hat\varepsilon_i^2$ (equations (5.26) and (5.27)); moreover, among all IV estimators that use linear combinations of $\mathbf{z}_i$ as instruments, 2SLS has the smallest asymptotic variance (Theorem 5.3).
:::

The proof sketch has four steps. Step one (2SLS.1 and 2SLS.2): substitute $y_i = \mathbf{x}_i'\beta + \varepsilon_i$; then $\sqrt n(\hat\beta_{2SLS} - \beta)$ equals $(\hat{\mathbf{X}}'\mathbf{X}/n)^{-1}$ times the transpose of the sample first-stage coefficients times $n^{-1/2}\sum_i\mathbf{z}_i\varepsilon_i$, and the first two factors converge in probability to their population counterparts. Step two (2SLS.1 and the moment conditions): the central limit theorem gives $n^{-1/2}\sum_i\mathbf{z}_i\varepsilon_i \xrightarrow{d} \mathcal{N}(0, \mathbb{E}[\varepsilon_i^2\mathbf{z}_i\mathbf{z}_i'])$; multiplying by the convergent matrices gives (i), and the sample counterparts give equation (5.34). Step three (2SLS.3): replacing the $\mathbb{E}[\varepsilon_i^2\mathbf{z}_i\mathbf{z}_i']$ in the filling by $\sigma^2\mathbb{E}[\mathbf{z}_i\mathbf{z}_i']$ collapses the sandwich, which is the variance of (ii); the residuals must be computed with $\mathbf{x}_i$, not with the fitted values (equation (5.25)). Step four (2SLS.1 through 2SLS.3): if another set of linear combinations of $\mathbf{z}_i$ is used as instruments, the inverse of the asymptotic variance of 2SLS minus the inverse of the asymptotic variance of that estimator equals $\sigma^{-2}$ times the second moment of "the residual from the linear projection of $\Pi_{\mathrm{fs}}'\mathbf{z}_i$ on that set of instruments"; a second-moment matrix is positive semidefinite, so the variance of 2SLS is no larger (p. 103).

In words: under 2SLS.1 and 2SLS.2, 2SLS still estimates $\beta$, and its variance is the sandwich of Section 3.4 of Chapter 5, whose bread is the explanatory power of the instruments for the regressors and whose filling is the variance of the product of instruments and error. 2SLS.3 decides only whether the sandwich can collapse; it has nothing to do with what is estimated. In the language of Chapter 7, 2SLS is GMM with weight $(\mathbf{Z}'\mathbf{Z}/n)^{-1}$, and under homoskedasticity this weight is the efficient weight (Theorem 3 of Chapter 7), which is what Theorem 5.3 says; under heteroskedasticity with overidentification, two-step GMM is more efficient. When the specification is wrong, for example when an instrument violates the exclusion restriction, the variance formula can be computed just the same; it merely describes the spread around a wrong center. Who sets the weights: the bread $\mathbb{E}[\Pi_{\mathrm{fs}}'\mathbf{z}_i\mathbf{z}_i'\Pi_{\mathrm{fs}}]$ contains only the part of the variation in the regressors that the instruments explain; weak instruments make the bread thin and the variance large, and this is the price 2SLS pays for consistency. In Bramble's buyer equation, the onboarding fee explains only a small part of the variation in the number of sellers, and the 2SLS standard error of 0.102 is about seven times that of OLS; the robust and homoskedastic standard errors, 0.1017 and 0.1010, are almost the same and agree with the Monte Carlo standard deviation of 0.096.

### 3.4 Tests of Endogeneity and of Overidentification

The standard error of 2SLS is several times that of OLS; if the number of sellers is in fact exogenous, OLS is better. Can this be tested? Conversely, the seller equation has two instruments and needs only one; can the extra one be used to check the instruments? The answer to both questions is a test in regression form.

::: {.theorem}
**Theorem 4 (tests of endogeneity and overidentification; Wooldridge 2010 §6.2 equations (6.13) to (6.15), pp. 126-127; §6.3.1 equation (6.22), p. 131; §6.3.2 equation (6.31), p. 135; §8.5.2 equation (8.58), p. 228)** Suppose 2SLS.1 and 2SLS.2 hold.

\(i) Let $\hat{\mathbf{v}}_i$ be the first-stage residuals of the endogenous regressors on $\mathbf{z}_i$ (only the components for the endogenous regressors), and regress $y_i$ on $(\mathbf{x}_i, \hat{\mathbf{v}}_i)$ by OLS:

$$y_i = \mathbf{x}_i'\beta + \hat{\mathbf{v}}_i'\rho + e_i.$$

The coefficient on $\mathbf{x}_i$ is numerically equal to $\hat\beta_{2SLS}$; under the null hypothesis $\rho = 0$ (exogenous regressors), for $\hat\rho$ the robust $t$ or Wald statistic has the usual asymptotic distribution, with no correction needed for the first stage.

\(ii) If $r > k$ and 2SLS.3 holds, regress the 2SLS residuals on all instruments; $n$ times the uncentered $R^2$ of this regression (the Sargan statistic) is, under the null hypothesis $\mathbb{E}[\mathbf{z}_i\varepsilon_i] = 0$, distributed as $\chi^2_{r-k}$.

\(iii) If $r > k$, 2SLS.3 is not needed: at the two-step GMM estimate with weight $\hat S = n^{-1}\sum_i\hat\varepsilon_i^2\mathbf{z}_i\mathbf{z}_i'$, Hansen's $J = n\,\bar g_n(\hat\beta)'\hat S^{-1}\bar g_n(\hat\beta)$ is distributed as $\chi^2_{r-k}$, where $\bar g_n(\beta) = n^{-1}\sum_i\mathbf{z}_i(y_i - \mathbf{x}_i'\beta)$.
:::

The proof sketch has four steps. Step one (2SLS.1 and linear projection): project the structural error linearly on the first-stage error, $\varepsilon_i = \mathbf{v}_i'\rho + e_i$; $e_i$ is uncorrelated with both $\mathbf{v}_i$ and $\mathbf{z}_i$, and hence with $\mathbf{x}_i$ (equations (6.13) and (6.14)); the regressors are exogenous if and only if $\rho = 0$ (equation (6.22)). Step two: using $\hat{\mathbf{v}}_i$ in place of $\mathbf{v}_i$ brings in a generated regressor, the case of Section 3.6 of Chapter 5; by the Frisch-Waugh-Lovell theorem, the coefficient on $\mathbf{x}_i$ is exactly 2SLS (Wooldridge 2010 Problem 5.1). When $\rho = 0$, the derivative of the second-step score with respect to the first-step parameters has expectation zero, that is, the $\Gamma_0 = 0$ of Theorem 6 of Chapter 5, so the first step's estimation error does not enter the asymptotic distribution of $\hat\rho$, which is (i). Step three (Theorem 4 of Chapter 7): (iii) is the special case of the J test of Chapter 7 for the moment function $\mathbf{z}_i(y_i - \mathbf{x}_i'\beta)$; estimating $\beta$ uses up $k$ directions, leaving $r - k$. Step four (2SLS.3): under homoskedasticity $S_0 = \sigma^2\mathbb{E}[\mathbf{z}_i\mathbf{z}_i']$, estimated by $\hat\sigma^2\mathbf{Z}'\mathbf{Z}/n$; efficient GMM is then 2SLS, and $J$ becomes $\hat{\boldsymbol\varepsilon}'\mathbf{P}_{\mathbf{Z}}\hat{\boldsymbol\varepsilon}/\hat\sigma^2$, that is, $n$ times the uncentered $R^2$, which is (ii); without 2SLS.3, this denominator no longer estimates $S_0$ consistently, and the $\chi^2$ limit fails with it.

In words: (i) estimates $\rho$, the projection coefficient of the structural error on the first-stage error; it tests whether OLS and 2SLS have the same probability limit, which is the principle of Hausman (1978) for comparing two estimators. This regression form is asymptotically equivalent to the Durbin-Wu-Hausman test, which compares the two estimators directly (Durbin 1954; Wu 1973; Wooldridge 2010 pp. 130-131). It presupposes valid instruments: if the instrument itself is wrong, the comparison is between OLS and an estimator that is just as inconsistent. In Bramble's buyer equation, $\hat\rho = 1.26$ (robust standard error 0.035) and $t = 35.53$: the number of sellers is clearly endogenous. Tests (ii) and (iii) check whether the $r$ moment conditions can all be satisfied by the same $\beta$, not whether the instruments are valid: if several instruments go wrong in the same direction, the test cannot detect it (Wooldridge 2010 pp. 135-136). When effects vary across markets, different instruments can aim at different weighted averages; each instrument is then right in its own sense, and $J$ may still reject. Who sets the weights: $J$ weights by $\hat S^{-1}$, so the more precisely estimated a moment, the more its departure from zero counts; Sargan uses $(\hat\sigma^2\mathbf{Z}'\mathbf{Z}/n)^{-1}$, which is the efficient weight only under 2SLS.3. In the seller equation of S1, $J = 1.690$ ($p = 0.194$) and Sargan is 1.904 ($p = 0.168$), with one degree of freedom; neither rejects. In S4v the errors vary with the referral bonus, and across 2,000 samples with valid instruments the nominal 5% Sargan test rejects 7.20% of the time and $J$ 5.35%.

### 3.5 Weak Instruments

Block S3's pilot moved the onboarding fee by only 34.9 cents. In this sample the first-stage F statistic is 0.019, in the lower tail of its sampling distribution (with $\mu^2 = 1$ the mean of F is about 2); 2SLS gives $-4.80$ (robust standard error 53.50), and the 95% interval $[-109.66, 100.07]$ covers 1.2 but is useless. Repeated sampling tells the story: across 2,000 replications with $\mu^2 = 1$, the median of 2SLS is 1.544, between the truth 1.2 and the OLS probability limit 2.408, and the nominal 95% Wald interval covers only 82.45%. The OLS probability limit of S3 is even higher than the 2.264 of S1: the fee barely moves, so almost all the remaining variation in the number of sellers comes from the local shocks.

::: {.intuition}
The trouble with weak instruments lies in the denominator of a ratio. With a single instrument, 2SLS is "the fee's effect on buyers" divided by "the fee's effect on sellers." The pilot moved the fee by only a few dozen cents, so the denominator is almost all noise, and when it happens to come close to zero, the ratio lands far from the truth. The noise in the numerator and the denominator comes from the same local shocks and is correlated, so the center of the ratio is pulled toward OLS. The Anderson-Rubin test sidesteps the division: given a candidate effect, it asks only whether "the buyer-equation errors computed with this effect are unrelated to the onboarding fee," a question that involves only addition and subtraction and can be answered whatever the strength of the instruments.
:::

::: {.assumption}
**WI.1 (the local-to-zero setup of Staiger and Stock 1997; following the assumptions of their NBER Technical Working Paper 151, hereafter TWP 151, p. 5)** There is a single endogenous regressor; the first-stage coefficients of the excluded instruments go to zero at the rate $n^{-1/2}$, with $\mu^2$ remaining bounded; the instruments are exogenous (2SLS.1), and the errors are homoskedastic.
:::

This setup does not claim that the data are really generated this way. It is an approximation: only when the first-stage signal is of the same order as the sampling noise does the asymptotic distribution keep the features that weak instruments show in finite samples (Wooldridge 2010 p. 111).

::: {.theorem}
**Theorem 5 (2SLS and the Anderson-Rubin test under weak instruments; Staiger and Stock 1997; theorem and equation numbers follow Theorem 1, Theorem 3(a), and equation (2.24) of TWP 151, p. 8 and pp. 15-16)** Under WI.1:

\(i) $\hat\beta_{2SLS}$ is inconsistent and converges in distribution to a nonnormal random variable whose center is pulled toward the probability limit of OLS; the conventional $t$ statistic is not standard normal; the first-stage F statistic converges to a noncentral $\chi^2$ divided by the number of excluded instruments, with noncentrality parameter $\mu^2$.

\(ii) For the null hypothesis $H_0: \beta = \beta^{\circ}$, subtract from $y_i$ the product of $\beta^{\circ}$ and the endogenous regressor, regress the result on the excluded instruments (controlling for the exogenous variables of the equation), and denote by $\mathrm{AR}(\beta^{\circ})$ the F statistic for the joint hypothesis that these coefficients are zero; under the null it converges to a $\chi^2$ divided by the number of excluded instruments, whatever $\mu^2$. The set of all $\beta^{\circ}$ not rejected by the AR test covers the truth with asymptotic probability equal to the nominal level.
:::

Here $\beta^{\circ}$ is the value under the null hypothesis; in this book $\beta_0$ denotes an intercept.

The proof sketch has three steps. Step one (WI.1): in the sample covariance between the excluded instruments and the regressor, multiplied by $\sqrt n$, the signal part is bounded and the noise part is asymptotically normal, so the two are of the same order and the first stage remains random in the limit; the limit of the F statistic is therefore a noncentral $\chi^2$. Step two (WI.1 and correlated errors): $\hat\beta_{2SLS} - \beta$ is the ratio of "the covariance of the instrument with the structural error" to "the covariance of the instrument with the regressor," and the two noises are linked through the correlation between the structural error and the first-stage error; the weaker the signal, the more the denominator looks like pure noise and the closer the ratio comes to the bias of OLS (TWP 151 Theorem 1). Step three (the null hypothesis only): when $\beta^{\circ}$ equals the truth, $y_i$ minus $\beta^{\circ}$ times the regressor is the structural error itself, whose sample covariance with the instruments does not involve the first stage, and the central limit theorem gives the limit directly (Theorem 3(a)). The argument is the same as for Theorem 6 of Chapter 7. With the moment function $\mathbf{z}_i(y_i - \mathbf{x}_i'\beta^{\circ})$, the $\mathcal{S}_n(\beta^{\circ})$ of Chapter 7 is the heteroskedasticity-robust AR; replace its $\hat S_n(\beta^{\circ})$ by the homoskedastic estimate (the residual variance times $n^{-1}\sum_i\mathbf{z}_i\mathbf{z}_i'$), partial out the exogenous variables of the equation, and divide by the number of excluded instruments, and the result is the classical AR in F form; the two differ only by degrees of freedom and a finite-sample correction of the denominator (for the general formulation, see Andrews, Stock and Sun 2019, pp. 22-23 of the 2018 draft).

In words: under WI.1, 2SLS no longer estimates any definite quantity; its sampling distribution is governed by the noise of the first stage and centered between the truth and OLS, and which way it leans is set by $\mu^2$: the smaller $\mu^2$, the closer to OLS. The AR set does not estimate a point; it collects every $\beta^{\circ}$ compatible with "the excluded instruments are uncorrelated with the structural error." When the specification is wrong and an instrument violates the exclusion restriction, AR under overidentification is also testing the overidentifying restrictions, and the set shrinks and may even be empty (Section 5.1 has an example); under exact identification the set is instead centered on the wrong value. Who sets the weights: AR needs no first-stage weights and uses only the reduced form. But it tests the moments of all excluded instruments jointly, so under overidentification with strong instruments it is wider than the Wald interval. In Bramble's S3 sample, $\mathrm{AR}(1.2) = 0.417$ ($p = 0.518$) and $\mathrm{AR}(0) = 0.074$ ($p = 0.785$), neither of which rejects; the AR set is the two rays $(-\infty, 2.191] \cup [2.718, \infty)$, which cover the truth. Across 2,000 replications with $\mu^2 = 1$, the coverage rate of the AR set is 0.953 and that of the Wald interval 0.8245.

An unbounded set is an honest answer, not a defect of the method. The shape of the set comes from a quadratic inequality: requiring $\mathrm{AR}(\beta^{\circ})$ not to exceed the critical value is equivalent to requiring a quadratic in $\beta^{\circ}$ not to exceed zero, and the coefficient of the quadratic term has the same sign as "the first-stage F minus the same critical value." When the first-stage F is above the critical value, the solution set is a bounded interval (possibly empty under overidentification); when it is below, the solution set is two rays or the whole real line. As the absolute value of $\beta^{\circ}$ goes to infinity, $\mathrm{AR}(\beta^{\circ})$ tends to the first-stage F statistic, that is, to a test of "the onboarding fee does not move the number of sellers." The first-stage F of the S3 sample is 0.019, far below the 5% critical value of 3.848, so the data cannot rule out that the fee has no effect on sellers at all; if that were so, every effect large enough in absolute value would be compatible with the data, so the set has to be unbounded. The only excluded stretch is the short one between 2.191 and 2.718, which contains the OLS estimate 2.418. Reporting a finite Wald interval here amounts to claiming knowledge that the data do not contain. With $\mu^2 = 1$, 69.9% of the AR sets are two rays, 13.5% are the whole real line, and only 16.6% are bounded.

How weak is weak? The most common answer is the rule of thumb that "a first-stage F statistic (written $F_{\mathrm{fs}}$ below) under 10 means weak," which Section 3.6 of Chapter 7 mentioned. Stock and Yogo (2005) interpret it as a test (following §4.3 and Tables 1 and 2 of NBER Technical Working Paper 284, pp. 29-30 and pp. 58-59). With the bias of 2SLS relative to OLS as the criterion and at least three excluded instruments, the 5% critical values for testing the weak-instrument null hypothesis "the maximal relative bias exceeds 10%" all lie between 9.08 and 11.52, so requiring F greater than 10 amounts roughly to this test. With the actual rejection rate of the Wald test as the criterion and a single instrument, guaranteeing that a nominal 5% test rejects at most 10% of the time requires F above 16.38; F greater than 10 guarantees only a rate of at most 15%, with a corresponding critical value of 8.96. All of these critical values assume homoskedasticity. When the errors are heteroskedastic, serially correlated, or clustered by group, use the effective F statistic $F_{\mathrm{eff}}$ of Montiel Olea and Pflueger (2013); for its simplified test with a bias tolerance of 10%, the critical value is at most 23.1, and with a single instrument it is 23.11 (following their April 2013 working-paper version, p. 3 and Table 1).

With a single instrument, the tF of Lee, McCrary, Moreira and Porter (2022) takes a different approach: instead of first screening samples by F, it lets the critical value of the $t$ statistic, $c_{\mathrm{tF}}(F_{\mathrm{fs}})$, vary smoothly with the first-stage F (following the March 2022 revision of NBER Working Paper 29124, §II.A to §II.B and Table 3, pp. 14-17). At the 5% level, the critical value returns to 1.96 only when $F_{\mathrm{fs}}$ reaches about 104.7; at $F_{\mathrm{fs}} = 10$ the standard error must be multiplied by 1.751; and when $F_{\mathrm{fs}}$ is at most 3.84 (the 5% critical value of $\chi^2_1$) the interval is the whole real line; if AR uses the same critical value, this is exactly the condition under which the AR set is unbounded, and the F-form value 3.848 used above differs only by a degrees-of-freedom correction. The first-stage F of the buyer equation of S1 is 174.7, so the conventional interval needs no adjustment; the F of the S3 sample is only 0.019, and tF, like AR, gives an unbounded set.

### 3.6 SUR and Feasible Generalized Least Squares

The two equations of block S2 have no endogenous regressors, so equation-by-equation OLS is consistent. But their errors are correlated: a local income shock in a market raises both the basket value and the number of orders, and in the sample the residuals of the two equations have correlation 0.605. The residuals of the other equation carry information about this equation's error; can that information be borrowed?

For unit $i$, stack the $M$ equations as $\mathbf{y}_i = \mathbf{X}_i\beta + \boldsymbol\varepsilon_i$, where $\mathbf{y}_i$ is $M \times 1$ and $\mathbf{X}_i$ is $M \times k$. For SUR, $\mathbf{X}_i$ is the block-diagonal matrix formed from the regressors $\mathbf{x}_{im}'$ of the individual equations, and $k = \sum_m k_m$. OLS on the whole system is called system OLS (SOLS), and for SUR it is just equation-by-equation OLS. GLS (generalized least squares) weights by $\Omega^{-1}$; since $\Omega$ is unknown, FGLS (feasible GLS) first estimates $\hat\Omega = n^{-1}\sum_i\hat{\boldsymbol\varepsilon}_i\hat{\boldsymbol\varepsilon}_i'$ from the residuals of equation-by-equation OLS (equation (7.40)) and then weights by it (equation (7.45); for SUR this approach goes back to Zellner 1962):

$$\hat\beta_{FGLS} = \Big(\sum_{i=1}^n\mathbf{X}_i'\hat\Omega^{-1}\mathbf{X}_i\Big)^{-1}\sum_{i=1}^n\mathbf{X}_i'\hat\Omega^{-1}\mathbf{y}_i.$$

::: {.assumption}
**SOLS.1 (Wooldridge 2010 §7.3.2, p. 167)** $\mathbb{E}[\mathbf{X}_i'\boldsymbol\varepsilon_i] = 0$.

**SOLS.2 (p. 168)** $\mathbb{E}[\mathbf{X}_i'\mathbf{X}_i]$ is nonsingular.

**SGLS.1 (Wooldridge 2010 §7.4, p. 173)** $\mathbb{E}[\mathbf{X}_i \otimes \boldsymbol\varepsilon_i] = 0$, that is, every element of $\mathbf{X}_i$ is uncorrelated with every element of $\boldsymbol\varepsilon_i$.

**SGLS.2 (p. 174)** $\Omega$ is positive definite, and $\mathbb{E}[\mathbf{X}_i'\Omega^{-1}\mathbf{X}_i]$ is nonsingular.

**SGLS.3 (Wooldridge 2010 §7.5.2, p. 180)** $\mathbb{E}[\mathbf{X}_i'\Omega^{-1}\boldsymbol\varepsilon_i\boldsymbol\varepsilon_i'\Omega^{-1}\mathbf{X}_i] = \mathbb{E}[\mathbf{X}_i'\Omega^{-1}\mathbf{X}_i]$.
:::

SOLS.1 requires only that each equation's regressors be uncorrelated with that equation's own error; SGLS.1 is stronger, requiring them to be uncorrelated with the error of every equation; this is the price FGLS pays for borrowing the residuals of other equations. In panel data, where the equations become the periods of the same unit, SGLS.1 is strict exogeneity. SGLS.3 is system homoskedasticity, and $\mathbb{E}[\boldsymbol\varepsilon_i\boldsymbol\varepsilon_i' \mid \mathbf{X}_i] = \Omega$ is a sufficient condition for it (equation (7.53)).

::: {.theorem}
**Theorem 6 (SOLS and FGLS; Wooldridge 2010 Theorems 7.1 and 7.2, pp. 169-171; Theorem 7.3, p. 179; Theorem 7.4, p. 181; §7.5.3 equation (7.55), p. 182)**

\(i) Under SOLS.1 and SOLS.2, SOLS is consistent and asymptotically normal, with a sandwich variance.

\(ii) Under SGLS.1 and SGLS.2, with $\hat\Omega$ taken from equation (7.40) (consistent under SOLS.2), FGLS is asymptotically equivalent to GLS with known $\Omega$: $\sqrt n$ times the difference between the two goes to zero in probability; FGLS is asymptotically normal, and its fully robust variance estimator is equation (7.52).

\(iii) If $\hat\Omega$ converges to some positive definite matrix $\Omega_*$, not necessarily equal to $\Omega$, then as long as SGLS.1 holds and $\mathbb{E}[\mathbf{X}_i'\Omega_*^{-1}\mathbf{X}_i]$ is nonsingular, FGLS is still consistent and asymptotically normal and is asymptotically equivalent to GLS weighted by $\Omega_*$, and its variance must use equation (7.52).

\(iv) If SGLS.3 also holds, $\mathrm{Avar}\big(\sqrt n(\hat\beta_{FGLS} - \beta)\big) = \big(\mathbb{E}[\mathbf{X}_i'\Omega^{-1}\mathbf{X}_i]\big)^{-1}$, estimated by $\big(\sum_i\mathbf{X}_i'\hat\Omega^{-1}\mathbf{X}_i\big)^{-1}$ (equation (7.54)); under SGLS.1, SOLS.2, SGLS.2, and equation (7.53), FGLS is more efficient than SOLS (pp. 181-182); more generally, it is asymptotically efficient among GMM estimators that use only the moment conditions of SGLS.1 (p. 222).
:::

The proof sketch has three steps. Step one (SOLS.1 and SOLS.2): $\beta = \big(\mathbb{E}[\mathbf{X}_i'\mathbf{X}_i]\big)^{-1}\mathbb{E}[\mathbf{X}_i'\mathbf{y}_i]$, and the law of large numbers and the central limit theorem give (i). Step two (SGLS.1): for any positive definite $\Omega_*$, $\mathrm{vec}\,\mathbb{E}[\mathbf{X}_i'\Omega_*^{-1}\boldsymbol\varepsilon_i] = \mathbb{E}[\boldsymbol\varepsilon_i' \otimes \mathbf{X}_i']\,\mathrm{vec}\,\Omega_*^{-1} = 0$, so the weighted moment conditions hold whatever the weight, and GLS is consistent; this is the source of (iii). When $\hat\Omega^{-1}$ is replaced by its limit, the difference is a bounded quantity times a quantity that goes to zero (equations (7.46) to (7.48)), which gives the asymptotic equivalence in (ii) and (iii). Step three (SGLS.3): the filling equals the bread and the sandwich collapses, which is the variance of (iv). The efficiency is a special case of Theorem 8 in Section 3.7: take the regressors of all equations of the system together and use them as instruments in every equation, and the resulting 3SLS is FGLS (p. 222). Hence, under SGLS.1 and $\mathbb{E}[\boldsymbol\varepsilon_i\boldsymbol\varepsilon_i' \mid \mathbf{X}_i] = \Omega$, FGLS is asymptotically efficient among GMM estimators that use only the moment conditions of SGLS.1, and SOLS is one of them.

In words: under SGLS.1 and SGLS.2, FGLS estimates the same $\beta$ that SOLS estimates, and the two differ only in precision. Even if the structure of $\hat\Omega$ is misspecified, FGLS remains consistent as long as the regressors are exogenous with respect to the error of every equation; only the conventional variance is wrong, and the robust equation (7.52) fixes it. The random effects estimator of Chapter 10 is a case in point: the equations become the $T$ periods of one unit, $\Omega$ is the $T \times T$ variance matrix of the composite error, strict exogeneity with respect to the composite error (individual effect included) is the SGLS.1 of that setting, and any positive definite $\hat\Omega$ gives a consistent estimate. When the structure of $\Omega$ is correctly specified, FGLS is asymptotically equivalent to GLS, and with equation (7.53) added it is asymptotically efficient among the estimators that use only the moment conditions of SGLS.1. When the specification is wrong and only SOLS.1 holds, so that the regressors of one equation are correlated with the error of another, FGLS borrows that other equation's error and may be inconsistent, while equation-by-equation OLS is unaffected. This is the chapter's first encounter with the trade-off between robustness and efficiency, which Section 3.8 meets again in system IV. Who sets the weights: $\hat\Omega^{-1}$ is set by the data, and the stronger the correlation of the errors, the more useful the residuals of the other equation. In Bramble's S2, the coefficients on the free-delivery minimum, which appears in both equations, gain only a little, with standard-error ratios of 0.944 and 0.983; the premium-seller share and the loyalty program, each of which appears in only one equation, both have their standard errors cut to 0.796 times those of OLS, 20.4% smaller.

::: {.theorem}
**Theorem 7 (when FGLS equals equation-by-equation OLS; Wooldridge 2010 Theorem 7.5, p. 185; Theorem 7.6, p. 186)** (i) If $\hat\Omega$ is a diagonal matrix, FGLS and equation-by-equation OLS are identical; if $\Omega$ is diagonal, the two are asymptotically equivalent. (ii) If every equation has the same regressors, $\mathbf{x}_{i1} = \cdots = \mathbf{x}_{iM}$, then whatever $\hat\Omega$ is, FGLS and equation-by-equation OLS are identical.
:::

The proof sketch has two steps. Step one ($\hat\Omega$ diagonal): $\hat\Omega^{-1}$ is diagonal, so the blocks of the equations are not connected, the weighting merely multiplies each equation by a constant, and the normal equations separate into equation-by-equation OLS; when $\Omega$ is diagonal, GLS is OLS, and FGLS is asymptotically equivalent to GLS (p. 186). Step two (identical regressors): write the common regressors as $\mathbf{x}_i$, so that $\mathbf{X}_i = I_M \otimes \mathbf{x}_i'$ and $\sum_i\mathbf{X}_i'\hat\Omega^{-1}\mathbf{X}_i = \hat\Omega^{-1} \otimes \sum_i\mathbf{x}_i\mathbf{x}_i'$; after inversion, $\hat\Omega^{-1}$ and $\hat\Omega$ cancel, and the coefficients of equation $m$ reduce to $\big(\sum_i\mathbf{x}_i\mathbf{x}_i'\big)^{-1}\sum_i\mathbf{x}_iy_{im}$ (Wooldridge 2010 Problem 7.5).

In words: with identical regressors, SUR contributes nothing to the point estimates, which estimate the same $\beta$ as equation-by-equation OLS; but the cross-equation covariance still matters, since any quantity that uses the coefficients of both equations at once needs it (Wooldridge 2010 p. 186). When the specification is wrong, for example when coefficients that are truly zero are put into every equation, the point estimates remain consistent and only the efficiency that comes from "different regressors" is lost. Who sets the weights: $\hat\Omega$ cancels out of the point estimates but stays in the variance, where it determines the covariance between the coefficients of the two equations. When S2 puts all four variables into both equations, the largest difference between SUR and OLS is $4.3 \times 10^{-12}$. The platform cares about revenue: raising the minimum by USD 1 changes monthly revenue per buyer by an estimated $-0.430$ dollars (truth $-0.518$), with standard error 0.169; the estimation errors of the two minimum coefficients are positively correlated (correlation 0.605), and if the two equations are treated as unrelated and the cross-equation covariance is dropped, the standard error is 0.139, an understatement of 18.1%. Both standard errors treat the sample means as fixed weights and do not account for the sampling error of the means themselves. The interval $[-0.762, -0.098]$ covers the truth.

### 3.7 System 2SLS, 3SLS, and System GMM

Back to block S1. Both equations have endogenous regressors, so SUR cannot be used; equation-by-equation 2SLS is consistent but makes no use of the correlation between the two equations' errors, and in the sample the 2SLS residuals of the two equations have a correlation of 0.858. 3SLS (Zellner and Theil 1962) combines 2SLS with FGLS: first run equation-by-equation 2SLS, estimate $\hat\Omega$ from the residuals, then weight by $\hat\Omega$ in system GMM. The buyer equation's coefficient changes from 1.140 with 2SLS to 1.212, and its standard error falls from 0.101 to 0.086.

Write the system for unit $i$ as $\mathbf{y}_i = \mathbf{X}_i\beta + \boldsymbol\varepsilon_i$, and let $\mathbf{Z}_i$ be its matrix of instruments, with $M$ rows and as many columns as the system has moment conditions ($Mr$ with common instruments); when every equation uses all exogenous variables as instruments, $\mathbf{Z}_i = I_M \otimes \mathbf{z}_i'$ (equation (9.27)). Stack the units' $\mathbf{y}_i$, $\mathbf{X}_i$, and $\mathbf{Z}_i$ into $\mathbf{y}$, $\mathbf{X}$, and $\mathbf{Z}$; the GMM estimator for the moment conditions $\mathbb{E}[\mathbf{Z}_i'(\mathbf{y}_i - \mathbf{X}_i\beta)] = 0$ is then (equation (8.28))

$$\hat\beta = \big(\mathbf{X}'\mathbf{Z}\hat W\mathbf{Z}'\mathbf{X}\big)^{-1}\mathbf{X}'\mathbf{Z}\hat W\mathbf{Z}'\mathbf{y}.$$

The weight $\hat W = (\mathbf{Z}'\mathbf{Z}/n)^{-1}$ gives system 2SLS (equations (8.32) and (8.33)), and $\hat W = \big(n^{-1}\sum_i\mathbf{Z}_i'\hat\Omega\mathbf{Z}_i\big)^{-1}$ gives 3SLS (equations (8.39) and (8.40)).

::: {.assumption}
**SIV.1 (Wooldridge 2010 §8.2, p. 210)** $\mathbb{E}[\mathbf{Z}_i'\boldsymbol\varepsilon_i] = 0$.

**SIV.2 (p. 210)** $\mathrm{rank}\,\mathbb{E}[\mathbf{Z}_i'\mathbf{X}_i] = k$.

**SIV.3 (Wooldridge 2010 §8.3.1, p. 215)** $\hat W \xrightarrow{p} W$, with $W$ symmetric and positive definite.

**SIV.5 (Wooldridge 2010 §8.3.4, p. 220)** $\mathbb{E}[\mathbf{Z}_i'\boldsymbol\varepsilon_i\boldsymbol\varepsilon_i'\mathbf{Z}_i] = \mathbb{E}[\mathbf{Z}_i'\Omega\mathbf{Z}_i]$.
:::

Wooldridge's SIV.4 says that the weight is the efficient weight; here it is written directly into the theorem. SIV.5 is the system version of 2SLS.3, and $\mathbb{E}[\boldsymbol\varepsilon_i\boldsymbol\varepsilon_i' \mid \mathbf{Z}_i] = \Omega$ is a sufficient condition for it (equation (8.41)); it restricts the conditional variance, not $\Omega$ itself, so the variances and covariances of the equations can be arbitrary.

Two more estimators go by the name 3SLS (Wooldridge 2010 §8.4, pp. 222-226). The GIV (generalized instrumental variables) estimator first applies a GLS-like transformation by $\hat\Omega^{-1/2}$ to the system, transforms the instruments in the same way, and then runs system 2SLS (equations (8.46) and (8.47)). Traditional 3SLS does not transform the first stage: it first takes the fitted values $\hat{\mathbf{X}}_i$ and then estimates the transformed equations with $\hat\Omega^{-1/2}\hat{\mathbf{X}}_i$ as instruments (equation (8.50)). The consistency of both requires $\mathbb{E}[\mathbf{Z}_i'\Omega^{-1}\boldsymbol\varepsilon_i] = 0$ (equations (8.48) and (8.51)), which is stronger than SIV.1: unless $\Omega$ is diagonal, it generally requires the instruments of each equation to be uncorrelated with the error of every equation. This is the system-IV counterpart of the distinction between SGLS.1 and SOLS.1 in Section 3.6. When every equation uses the same set of instruments (that is, any variable that is exogenous in one equation is exogenous in all of them), this condition follows directly from SIV.1, and the three versions of 3SLS coincide (Theorem 8 (iv)).

::: {.theorem}
**Theorem 8 (system IV; Wooldridge 2010 Theorems 8.1 and 8.2, pp. 215-216; equation (8.33), p. 217; Theorem 8.3, p. 218; Theorem 8.4 and equation (8.45), p. 221; §8.4.2, p. 225; equation (8.58), p. 228)**

\(i) Under SIV.1 through SIV.3, the GMM estimator is consistent and asymptotically normal, with the sandwich variance of Theorem 2 of Chapter 7 (equation (8.29)).

\(ii) With $W = S_0^{-1}$, where $S_0 = \mathbb{E}[\mathbf{Z}_i'\boldsymbol\varepsilon_i\boldsymbol\varepsilon_i'\mathbf{Z}_i]$, GMM has the smallest variance in this family (Theorem 8.3) and is called the minimum chi-square estimator, with $S_0$ estimated from the residuals of an initial estimate (equation (8.36)); when there are more moment conditions than $k$, the $J$ statistic at the efficient weight follows a $\chi^2$ distribution under SIV.1, with degrees of freedom equal to the number of moment conditions minus $k$ (equation (8.58)).

\(iii) If SIV.5 also holds, the 3SLS weight consistently estimates $S_0^{-1}$, 3SLS is efficient GMM, and its variance is estimated by $\big(\mathbf{X}'\mathbf{Z}\{\sum_i\mathbf{Z}_i'\hat\Omega\mathbf{Z}_i\}^{-1}\mathbf{Z}'\mathbf{X}\big)^{-1}$ (equation (8.45)); without SIV.5, 3SLS is in general less efficient than the minimum chi-square estimator, and equation (8.45) no longer holds.

\(iv) When $\mathbf{Z}_i = I_M \otimes \mathbf{z}_i'$, system 2SLS is equation-by-equation 2SLS, and 3SLS is identical to the GIV estimator and to traditional 3SLS.
:::

The proof sketch has four steps. Step one (SIV.1 through SIV.3): the moment function $\mathbf{Z}_i'(\mathbf{y}_i - \mathbf{X}_i\beta)$ is linear in $\beta$, so Theorems 1 and 2 of Chapter 7 apply as they stand, and SIV.2 guarantees that the Jacobian $-\mathbb{E}[\mathbf{Z}_i'\mathbf{X}_i]$ has full column rank, which is (i). Step two (SIV.1 through SIV.3 and the efficient weight): the efficient weight $W = S_0^{-1}$ is Theorem 3 of Chapter 7, and $J$ is Theorem 4 of Chapter 7, which is (ii). Step three (SIV.5): $S_0 = \mathbb{E}[\mathbf{Z}_i'\Omega\mathbf{Z}_i]$, and the 3SLS weight is the inverse of a consistent estimate of it, so 3SLS is efficient GMM and the variance collapses to equation (8.45), which is (iii). Step four (algebra): when $\mathbf{Z}_i = I_M \otimes \mathbf{z}_i'$, $\mathbf{Z}'\mathbf{Z}$ is block diagonal, and the normal equations of system 2SLS separate equation by equation (p. 217). The equality of 3SLS and GIV is an algebraic result (cited in Wooldridge 2010 p. 225). GIV and traditional 3SLS differ only in the first stage, which the former in effect estimates by FGLS and the latter by OLS; since the first stages of all equations have the same regressors, the two give the same result (Theorem 7 (ii)), so all three coincide (p. 225), which is (iv).

In words: under SIV.1 through SIV.3, system 2SLS, 3SLS, and the minimum chi-square estimator all estimate the same $\beta$; they differ only in precision, and precision comes from the weight. The 3SLS weight is set by $\hat\Omega$ and is the efficient weight under SIV.5; with heteroskedastic errors it is still consistent, but its conventional standard errors are no longer correct, so either use the robust sandwich or switch to the minimum chi-square estimator. When the specification is wrong, that is, when the moment conditions of some equation fail, the system estimator converges to the pseudo-true value that brings the weighted moments closest to zero, and the failing moments drag in the other equations through the off-diagonal blocks of the weight, which Section 3.8 discusses in detail. The system of S1 has 10 moment conditions and 9 parameters, so the system $J$ has only one degree of freedom, which is the seller equation's overidentifying restriction. The errors of S1 are homoskedastic, so 3SLS is efficient here: the buyer equation's standard error falls from 0.101 to 0.086, a ratio of 0.853. The two homoskedastic standard errors have different denominators (2SLS divides by $n - k$, the $\hat\Omega$ of 3SLS divides by $n$); with $n$ used for both, the ratio is 0.855. The ratio of the population asymptotic standard errors is 0.857, and the ratio of the two sampling standard deviations across 1,000 Monte Carlo replications is 0.837, so all four ratios put the gain at about 15%. The robust 3SLS standard error, 0.082, is close to the conventional one. In panel data, where the $T$ periods of one unit are treated as $T$ equations, system 2SLS is pooled 2SLS when each period uses the same number of instruments and $\mathbf{Z}_i$ is stacked as in Wooldridge (2010) equation (8.22) (p. 217). The Hausman-Taylor estimator of Chapter 10 applies this to the equation after quasi-demeaning (defined in Chapter 10): its consistency and asymptotic normality come from Theorems 2 and 3 and from this theorem, its order and rank conditions are checked through 2SLS.2(b) and SIV.2 (Theorem 1 (i) and (ii) gives their form in simultaneous equations), and the strength of its first stage is judged as in Section 3.5.

### 3.8 The Trade-Off Between Robustness and Efficiency

The gain of 3SLS has a price. In block S1m the buyer equation is entirely correct, and only the seller equation omits the interaction between the coupon and the number of buyers. Across 1,000 Monte Carlo replications (figure below), the buyer-equation 2SLS estimates average 1.191, and the nominal 95% interval has a coverage rate of 0.957; the 3SLS estimates average 0.793, with coverage of only 0.001.

![Monte Carlo sampling distributions of the coefficient on the number of sellers in the buyer equation in blocks S1 and S1m, 1,000 samples each: the blue line is 2SLS on the buyer equation alone, the orange line is 3SLS on the two equations together, and the vertical dashed line is the truth, 1.2. In the left panel both equations are correctly specified and the 3SLS distribution is narrower; in the right panel the seller equation omits the interaction between the coupon and the number of buyers, 2SLS stays centered on the truth, and 3SLS shifts left as a whole.](assets/fig/fig_09_contagion.svg)

::: {.theorem}
**Theorem 9 (when 2SLS and 3SLS coincide, and the consequences of one misspecified equation; Wooldridge 2010 §9.3.2, pp. 254-255; §9.3.1, p. 252)** Suppose $\mathbf{Z}_i = I_M \otimes \mathbf{z}_i'$.

\(i) When every equation is exactly identified, 3SLS and equation-by-equation 2SLS are identical; they are also identical when $\hat\Omega$ is diagonal.

\(ii) When only one equation is overidentified and all the others are exactly identified, the 3SLS of the overidentified equation equals its 2SLS; the 3SLS of the exactly identified equations generally differs from 2SLS.

\(iii) If some equation is misspecified so that its moment conditions fail, 3SLS and system GMM are in general inconsistent for all parameters; for any equation that is correctly specified and has exogenous instruments, 2SLS remains consistent.
:::

The proof sketch has four steps. Step one (exact identification): the sample moments of every equation can be solved exactly to zero, so the weight plays no role; when $\hat\Omega$ is diagonal, the 3SLS weight is block diagonal and differs from that of system 2SLS only by a constant multiple in each block, which is (i). Step two (the other equations exactly identified, SIV.2): put the overidentified equation first, split the sample moments by equation into two blocks, $\bar g_{(1)}$ and $\bar g_{(2)}$, and partition the 3SLS weight $\hat W$ accordingly. The sample moments of the other equations are linear in their own parameters with an invertible coefficient matrix, so whatever the parameters of the first equation, the other parameters can always be chosen to set $\bar g_{(2)}$ to any value. Minimizing first over $\bar g_{(2)}$ gives $\bar g_{(2)} = -\hat W_{(22)}^{-1}\hat W_{(21)}\bar g_{(1)}$, and the criterion reduces to $\bar g_{(1)}'\big[(\hat W^{-1})_{(11)}\big]^{-1}\bar g_{(1)}$. For 3SLS, $\hat W^{-1} = \hat\Omega \otimes n^{-1}\sum_i\mathbf{z}_i\mathbf{z}_i'$, whose upper-left block is $[\hat\Omega]_{11}\,n^{-1}\sum_i\mathbf{z}_i\mathbf{z}_i'$, so the remaining criterion differs from the 2SLS criterion of the first equation only by a constant multiple, and the first equation's 3SLS is its 2SLS. Step three (overidentification): for the other equations, 3SLS sets $\bar g_{(2)}$ equal to $-\hat W_{(22)}^{-1}\hat W_{(21)}\bar g_{(1)}$ rather than to zero. The first equation is overidentified, so $\bar g_{(1)}$ is generally not zero at the estimate; as long as the off-diagonal elements of $\hat\Omega^{-1}$ are not zero, this term is nonzero. Equation-by-equation 2SLS, by contrast, sets $\bar g_{(2)}$ to zero, so the two differ, which is (ii) (an algebraic result cited in Wooldridge 2010 pp. 254-255). Intuitively, the only extra information the system provides is that equation's overidentifying restriction, which 2SLS has already used up in estimating that equation, so it can help only the other equations. Step four (misspecification): the 3SLS pseudo-true value sets a weighted combination of all moment conditions to zero, and the off-diagonal blocks of the weight contain the off-diagonal elements of $\Omega^{-1}$, so the failing moment of the seller equation enters the first-order conditions of the buyer equation, and the buyer equation's pseudo-true value departs from $\beta$. The weight of equation-by-equation 2SLS is block diagonal, so the buyer equation uses only its own moments and is not drawn in, which is (iii).

In words: the gain and the risk of 3SLS travel the same road. The gain comes only from the overidentifying restrictions of other equations, passed on to the exactly identified equations through the correlation of the errors. That is why, in S1, the 3SLS and 2SLS estimates of the seller equation differ by only $2.4 \times 10^{-10}$, and the entire gain falls on the exactly identified buyer equation. Which coefficient it reaches depends on how strongly the instruments involved in that restriction are correlated with the coefficient's source of identification: the seller equation's overidentifying restriction comes from advertising and the coupon, the city programs give the coupon and the onboarding fee a correlation of $-0.7$, and the fee is the source of identification for the coefficient on the number of sellers in the buyer equation. If the two were independent, the ratio of the population standard errors for the sellers coefficient would be exactly 1: 3SLS would not improve it at all, though it would still improve the precision of the two coefficients on advertising and the coupon in the buyer equation. When the specification is wrong, the same road carries the error across. In S1m, what fails is the seller equation's coupon moment. The buyer equation's coupon coefficient is identified by the coupon's own moment and is contaminated under both city-program designs: on one million markets, its 3SLS is 0.738 (correlation $-0.7$) and 0.732 (correlation zero), while its 2SLS is near the truth, 0.40, in both. The sellers coefficient is identified by the onboarding fee and is drawn in only when the coupon and the fee are correlated: with correlation $-0.7$, its 2SLS is about 1.20 and its 3SLS about 0.80; with correlation zero, 2SLS and 3SLS are 1.202 and 1.201, and it is no longer contaminated. So the contamination is not transferred from one coefficient to another: it always lands on the coupon coefficient and spreads to the sellers coefficient through the correlation between the coupon and the fee; in a misspecified system, other coefficients may shift along with them. Who sets the weights: $\Omega^{-1}$ is set by the data, and the researcher decides which equations go into the system. The seller equation's own $J$ is 67.89 in the S1m sample ($p = 1.7 \times 10^{-16}$) and 1.690 in S1; it is the first line of defense for spotting the problem.

This is the trade-off described in Wooldridge (2010) §9.3.1: when all equations are right, system methods are asymptotically more efficient; when only one equation matters, 2SLS is consistent as long as that equation is correctly specified and its instruments are exogenous. In practice both should be reported: when every equation is correctly specified, 3SLS and 2SLS have the same probability limit, so the difference between them is itself a diagnostic. In S1 the difference for the buyer equation is 0.072, of the same order as the 3SLS standard error; in S1m it is $-0.399$, far larger than any standard error. The efficiency of system estimation, then, comes at the cost of robustness.

## 4 Estimation and Practice

### 4.1 Computation: The Traps of Doing the Two Stages by Hand

The name 2SLS comes from two regression steps, but running the two steps by hand sets two traps (Wooldridge 2010 §5.1.2, p. 97). The first is the standard error: the residuals used by the second-step OLS are $y_i - \hat{\mathbf{x}}_i'\hat\beta$, while the correct residuals are $y_i - \mathbf{x}_i'\hat\beta$, and the $\hat\sigma^2$ computed from the wrong residuals is inconsistent. The second is leaving the equation's exogenous variables out of the first step: regressing the endogenous variable only on the excluded instruments and then putting the fitted values into the second step gives an estimator that is in general inconsistent (p. 97), because the first-step residuals are correlated with the exogenous variables of the equation. Software 2SLS commands handle both points correctly; when coding by hand, use $(\hat{\mathbf{X}}'\hat{\mathbf{X}})^{-1}\hat{\mathbf{X}}'\mathbf{y}$ directly and compute the residuals with the original $\mathbf{X}$. Nor does 3SLS need three steps: the single formula of equation (8.40) suffices, and the name is a historical holdover (p. 220). Its conventional variance, equation (8.45), relies on SIV.5; when the errors are heteroskedastic or must be clustered by unit, use the GMM sandwich or switch to the minimum chi-square estimator.

### 4.2 Diagnostics and Reporting

An IV estimate should come with at least the following:

- First stage: the coefficients and standard errors of the excluded instruments, and $F_{\mathrm{fs}}$; with heteroskedastic, serially correlated, or clustered errors, report $F_{\mathrm{eff}}$ and compare it with the critical values of Montiel Olea and Pflueger (2013), not with 10.
- With a single endogenous regressor, Andrews, Stock and Sun (2019; p. 3 of the 2018 draft) recommend judging instrument strength by the effective F; with a single instrument, they recommend reporting the AR interval however large F is. If $t$ inference is still reported, at least use the tF critical values of Lee et al. (2022) (Section 3.5).
- Under overidentification, report the robust $J$ with its degrees of freedom; the homoskedastic Sargan test has the correct rejection rate only when the errors are homoskedastic.
- Endogeneity: $\hat\rho$ from the control function regression, with its robust $t$.
- System estimation: equation-by-equation 2SLS side by side with 3SLS, and the $J$ of each equation reported separately; when the two differ sharply, suspect first that some equation is misspecified.

### 4.3 The Estimators at a Glance

The table below lines up the estimators of this chapter. LIML (limited information maximum likelihood) estimates the structural equation together with the reduced form of the endogenous regressors under a normal likelihood, using the same set of instruments.

| Estimator | Moment conditions and weight | Consistency requires | Efficient when | S1 buyer equation, coefficient on sellers |
|---|---|---|---|---|
| OLS | $\mathbb{E}[\mathbf{x}_i\varepsilon_i] = 0$ | Exogenous regressors | Exogenous and homoskedastic | 2.268 (0.015) |
| 2SLS | $\mathbb{E}[\mathbf{z}_i\varepsilon_i] = 0$, weight $(\mathbf{Z}'\mathbf{Z}/n)^{-1}$ | 2SLS.1, 2SLS.2 | 2SLS.3 | 1.140 (0.102) |
| Two-step GMM | As above, weight $\hat S^{-1}$ | As above | Always efficient for these moments | Exactly identified, equals 2SLS |
| LIML | As above | As above | Asymptotically equivalent to 2SLS; small median bias under weak instruments | Exactly identified, equals 2SLS |
| 3SLS | Moments of all equations, weight $\big(n^{-1}\sum_i\mathbf{Z}_i'\hat\Omega\mathbf{Z}_i\big)^{-1}$ | Every equation correctly specified | SIV.5 | 1.212 (0.086) |
| Minimum chi-square | As above, weight $\hat S^{-1}$ | As above | Always efficient for these moments | Asymptotically equivalent to 3SLS under SIV.5 |

To choose among them, first ask which number you want. When the structural coefficient is the same in every market and the instruments are strong and valid, 2SLS is enough; if you want more precision and are willing to answer for the specification of every equation, use 3SLS or the minimum chi-square estimator; when the instruments are weak, no point estimate is reliable, and inference should switch to AR or tF. LIML's advantage under weak instruments rests on a finding of Staiger and Stock: it is approximately median unbiased (the introduction of TWP 151). Under exact identification it is the same as 2SLS, and in the S3 sample both are $-4.795$; the Fuller correction ($b = 1$) gives 2.280 (standard error 0.460), not far from the OLS estimate of 2.418, yet the 95% interval built from this standard error does not contain the truth, 1.2. When the first-stage F is only 0.019, no point estimate helps; only the AR set remains honest.

### 4.4 Software

In R, single-equation 2SLS uses `ivreg::ivreg(y ~ x | z)`; `summary(fit, vcov. = function(x) sandwich::vcovHC(x, type = "HC0"), diagnostics = TRUE)` gives the robust weak-instrument F and the Wu-Hausman test, while its Sargan test is always the homoskedastic version; writing only `vcov. = sandwich::vcovHC` uses the default HC3, which differs from the HC0 of this chapter. `fixest::feols(y ~ w | x ~ z)` goes with `fitstat(fit, ~ ivf + ivwald + wh + sargan + cd)`: `ivf` is the homoskedastic first-stage F, and `ivwald` uses the variance specified at fitting, so with one endogenous regressor it can serve as the robust first-stage F. The packages used in this chapter offer neither the effective F nor tF, so the chapter's code computes them itself; the CRAN package `ivDiag` provides both (`eff_F()` and `tF()`), and tF can also be looked up in the tables the authors have published. The AR set, CLR (conditional likelihood ratio), LIML, and Fuller use `ivmodel` (`AR.test`, `CLR`, `LIML`, `Fuller`). SUR, system 2SLS, and 3SLS use `systemfit(eqs, method = "SUR", inst = ~ z1 + z2, data = d)`, where `method` can be changed to `"2SLS"` or `"3SLS"`, and `methodResidCov = "noDfCor"` makes $\hat\Omega$ divide by $n$, consistent with equation (7.40) and equation (8.38); `systemfit` has no cluster variance, which `sandwich::vcovCL(fit, cluster = rep(id, M), type = "HC0", cadjust = FALSE, residFit = FALSE)` supplies, with the cluster vector stacked by equation. `residFit = FALSE` cannot be omitted: without it, the score uses residuals computed from the fitted regressors, which is the trap of Section 4.1, and the IV standard errors are wrong; for the 3SLS of S1, the standard error of the sellers coefficient computed this way is much larger than the correct robust value of 0.082. For system GMM there is `momentfit`, with `sysMomentModel` and `gmmFit`. The Stata counterparts are `ivregress` (followed by `estat firststage`, `estat overid`, `estat endogenous`), `sureg`, and `reg3`.

## 5 Anchoring Papers

Each of the three anchoring papers corresponds to one thread of the chapter: building instruments from natural experiments, simultaneous equations on a two-sided platform, and honest inference when the first stage is not strong enough. Each paper is laid out in five parts: paper, method, data, results, and limitations.

### 5.1 Angrist and Krueger (1991)

::: {.case}
Paper: "Does Compulsory School Attendance Affect Schooling and Earnings?" Quarterly Journal of Economics. The question is the return to education: by what percentage does an extra year of schooling raise earnings? People of higher ability get more schooling and also earn more, so OLS is commonly believed to overstate the return.

Method: quarter of birth as the instrument. School districts set the school-entry age by date of birth, so children born early in the year are older when they start school and have completed fewer grades when they reach the age at which they may legally drop out; quarter of birth is assumed to be unrelated to ability and to affect earnings only through years of schooling. The authors then expand the instruments with interactions of quarter of birth with year of birth and with state of birth.

Data: men from the 1970 Census (born 1920 to 1929) and from the 5% Public Use Microdata Sample of the 1980 Census (born 1930 to 1939 and 1940 to 1949), grouped by decade of birth; the two 1980 Census samples each have between three hundred thousand and five hundred thousand men (Wooldridge 2010 p. 111). The reanalysis by Staiger and Stock below uses only the men in the 1980 Census born 1930 to 1949 with earnings in 1979 (TWP 151 p. 29).

Results: the authors conclude that the IV estimates of the return are quite close to OLS, so the bias due to ability is small. Staiger and Stock reestimate the 1980 Census part (TWP 151 Table 11): among the 329,509 men born 1930 to 1939, OLS is 0.0632; with only the three quarter-of-birth dummies as instruments, the first-stage F is 30.53, 2SLS is 0.0990 (standard error 0.0207), and the AR interval is $[0.052, 0.153]$; adding age controls and quarter-of-birth by year-of-birth interactions, which brings the instruments to 28, the first-stage F falls to 1.613 and 2SLS becomes 0.0600 (standard error 0.0290), almost on top of OLS, while the AR interval widens to $[-0.441, 0.490]$.

Limitations: with many weak instruments, 2SLS is pulled toward OLS, as Theorem 5 (i) predicts, so in these specifications "IV is close to OLS" may be merely an artifact of weak instruments; Bound, Jaeger and Baker (1995) pointed out the problems of 2SLS in this setting (Wooldridge 2010 p. 111). Staiger and Stock conclude that, by the statistics they consider reliable, the return to education is higher than in the original paper, but the intervals are wider than the original 2SLS intervals; for the cohort born 1940 to 1949, the overidentifying restrictions are rejected and the AR sets of the first three specifications are empty, so they rely only on the 1930 to 1939 cohort.
:::

This study made the "natural experiment" the standard narrative of IV: the instrument comes from an accidental detail of institutions rather than from another economic variable. It also became the starting point of the weak-instruments literature: even a large sample of more than three hundred thousand people cannot make up for a weak first stage.

### 5.2 Chu and Manchanda (2016)

::: {.case}
Paper: "Quantifying Cross and Direct Network Effects in Online Consumer-to-Consumer Platforms," Marketing Science. The question is the same as Bramble's: on a C2C (consumer-to-consumer) platform such as Taobao, how many buyers do more sellers attract, and how many sellers do more buyers attract?

Method: starting from a utility model, the authors write one growth equation for buyers and one for sellers, and each side's installed base appears in the other side's equation. The two equations are estimated jointly, with errors specified as correlated bivariate normal to absorb common shocks that hit both sides at once. Identification rests on three things: each side's own excluded variables (product variety enters only the buyer equation, buyer quality only the seller equation); external instruments (national indices of consumer expectations, confidence, and satisfaction for the buyer base, and an entrepreneur confidence index and components of the purchasing managers' index for the seller base); and the control function of Petrin and Train (2010).

Data: daily data on the installed bases and new users of both sides of Taobao, from May 11, 2003, to December 31, 2012; the installed base is the cumulative number of registered users (pp. 876-877).

Results: the cross-side network effects on both sides are significantly positive and asymmetric, with the effect of the seller base on buyer growth far larger than the effect of the buyer base on seller growth; direct network effects are small, with only a small positive effect on the buyer side; the cross-side effect of sellers on buyers first rises and then falls over time before stabilizing. The authors report that the OLS and 2SLS estimates differ only slightly (Table 5, p. 883).

Limitations: the authors acknowledge that within the same population it is hard to find instruments that affect only one side (p. 882); the paper reports no first-stage F or overidentification test, so the validity of the excluded variables rests on economic argument alone. The joint estimation uses maximum likelihood rather than 3SLS, but its dependence on the specification of both equations is the same as described in Section 3.8: one wrong equation drags down the other.
:::

The paper's value is that it puts Bramble's problem on a real platform: each side needs a variable that moves only its own side before the cross-side effects can be separated out of equilibrium data. It also shows how hard such variables are to find in real data.

### 5.3 Lee, McCrary, Moreira and Porter (2022)

::: {.case}
Paper: "Valid t-Ratio Inference for IV," American Economic Review. This section follows the March 2022 revision of its NBER Working Paper 29124.

Method: with one instrument and one endogenous regressor, the conventional practice is to check first whether the first-stage F exceeds 10 and then multiply the standard error by $\pm 1.96$. Using the asymptotic representation of Stock and Yogo (2005), the authors construct a critical value $c_{\mathrm{tF}}(F_{\mathrm{fs}})$ that varies smoothly with F and guarantees that the probability of $t^2$ exceeding it is at most the nominal level; restated as an adjustment factor for the standard error, it multiplies the conventional standard error by a number that decreases as F increases (Section 3.5). As long as the first stage and 2SLS use the same robust variance, it applies equally to heteroskedasticity and clustering.

Data: of the 757 papers in the American Economic Review from 2013 to 2019, 123 use IV, and 61 of these use only exactly identified single-instrument specifications, 1,311 specifications in all.

Results: the weighted quartiles of the first-stage F are 14.23, 45.84, and 225, so most specifications exceed 10. After the tF adjustment, the confidence intervals of a quarter of the specifications are at least 49% wider; among specifications with F above 10 and $|t|$ above 1.96, which by convention would be treated as significant, 24% are no longer significant at the 5% level. The authors also prove that when both intervals are bounded, the tF interval has a shorter expected length than the AR interval.

Limitations: it applies only to a single instrument; it controls the unconditional rejection rate, not the rejection rate conditional on the observed F; with several instruments one still has to return to methods of the AR kind.
:::

Read together, the three papers line up as follows: quarter of birth is the model of IV as a natural experiment, and the reanalysis by Staiger and Stock shows that the same data tell an entirely different story when the instruments become weak; Chu and Manchanda take up Bramble's own problem on a two-sided platform, with identification coming from each side's excluded variables; Lee et al. show that a first-stage F above 10 does not guarantee that conventional $t$ inference is credible.

## 6 A Full Walkthrough on the Bramble Data

The code below runs under R 4.6.1. It lists only the statements that generate the data, estimate, and summarize the Monte Carlo; several statements are merged onto one line, and statements used only for checks are omitted. Every number in the text comes from the output of an actual run, and the result lines in the code that start with `#` copy what the run printed. Unless stated otherwise, 2SLS standard errors are HC0 robust standard errors (equation (5.34)), and 3SLS standard errors are those of equation (8.45). For a few small functions only the purpose is described: `first_stage_F` and `eff_F` compute the first-stage F and the effective F, `cf_test` runs the control function regression of Theorem 4 (i), `hansen_j` and `sargan` compute the two overidentification statistics, `ar_test`, `ar_ci`, and `ar_covers` compute the AR statistic, the closed-form AR set, and whether it contains a given value, and `rank_condition` computes the rank of a matrix from its singular values.

### 6.1 The DGP

Blocks S1 and S1m share one generating function: the two structural equations are solved jointly market by market, and S1m adds only the interaction coefficient `kappa`. Seeds are fixed by block.

```r
P1 <- list(n = 1500, ad0 = 40, ad_sd = 10, cp0 = 6, cp_sd = 1.5, fee0 = 30, fee_sd = 7, r = -0.7,
           aB = 14.4, bBS = 1.2, gad = 0.10, gcp = 0.40, lB = 2.0, sB = 2.5,
           aS = 9.0, bSB = 0.25, dfee = -0.10, lS = 0.6, sS = 0.8, rho = 0.8)
gen_S1 <- function(par, n = par$n, kappa = 0) {
  log_pop <- rnorm(n)                                                 # standardized market size
  ad_budget <- par$ad0 + par$ad_sd * rnorm(n)                         # USD thousands, reallocated at random
  u <- rnorm(n); w <- rnorm(n)
  coupon <- par$cp0 + par$cp_sd * u                                   # first-order coupon, USD
  fee <- par$fee0 + par$fee_sd * (par$r * u + sqrt(1 - par$r^2) * w)  # onboarding fee, USD per month
  Om <- matrix(c(par$sB^2, par$rho * par$sB * par$sS, par$rho * par$sB * par$sS, par$sS^2), 2)
  e <- matrix(rnorm(2 * n), n) %*% chol(Om)                           # local shocks hit both sides
  cB <- par$aB + par$gad * ad_budget + par$gcp * coupon + par$lB * log_pop + e[, 1]
  slope <- par$bSB + kappa * (coupon - par$cp0)                       # S1m: kappa = -0.0045
  cS <- par$aS + par$dfee * fee + par$lS * log_pop + e[, 2]
  buyers <- (cB + par$bBS * cS) / (1 - par$bBS * slope)               # both equations solved market by market
  data.frame(buyers = buyers, n_sellers = cS + slope * buyers, ad_budget = ad_budget, coupon = coupon,
             fee = fee, log_pop = log_pop)
}
set.seed(9101); S1 <- gen_S1(P1)
set.seed(9102); S1m <- gen_S1(P1, kappa = -0.0045)
Gam <- matrix(c(-1, 1.2, 0.25, -1), 2)                                # columns: buyer and seller equation
Del <- cbind(c(14.4, 0.10, 0.40, 0, 2.0), c(9.0, 0, 0, -0.10, 0.6))   # rows: const, ad_budget, coupon, fee, log_pop
Pi <- -Del %*% solve(Gam)                                             # reduced form, eq. (9.14)
dimnames(Pi) <- list(c("const", "ad_budget", "coupon", "fee", "log_pop"), c("buyers", "n_sellers"))
round(Pi[2:4, ], 4)
c(buyer = rank_condition(Pi["fee", "n_sellers", drop = FALSE])$rank,
  seller = rank_condition(Pi[c("ad_budget", "coupon"), "buyers", drop = FALSE])$rank)
#            buyers n_sellers
# ad_budget  0.1429    0.0357
# coupon     0.5714    0.1429
# fee       -0.1714   -0.1429
#  buyer seller
#      1      1
```

The reduced form agrees with Section 2.1: for each USD 1 increase in the onboarding fee, sellers drop by 0.143 and buyers drop by 0.171 (in hundreds and thousands); the identifying submatrices of both equations have rank 1, with the buyer equation using the fee's row and the seller equation the two rows of advertising and the coupon.

### 6.2 Block S1: OLS, 2SLS, and 3SLS

`tsls` is the 2SLS of equation (5.22), with a choice of the robust variance of equation (5.34) or the homoskedastic one of equation (5.27); `sys_gls` solves the weighted normal equations of the system, and `three_sls` estimates $\hat\Omega$ from the residuals of equation-by-equation 2SLS and then calls it.

```r
tsls <- function(y, X, Z, vcov = c("robust", "homo")) {               # eqs. (5.22), (5.34) and (5.27)
  vcov <- match.arg(vcov)
  Xh <- qr.fitted(qr(Z), X); qx <- qr(Xh)                             # first stage: X_hat = P_Z X
  b <- qr.coef(qx, y); A <- chol2inv(qr.R(qx))                        # b = (X_hat'X_hat)^-1 X_hat'y
  e <- drop(y - X %*% b)                                              # residuals with X, not X_hat
  V <- if (vcov == "robust") A %*% crossprod(Xh * e) %*% A else sum(e^2) / (length(y) - ncol(X)) * A
  list(b = b, V = V)
}
sys_gls <- function(ylist, L, R, Oi) {                                # A[m, l] = [Omega^-1]_ml L_m'R_l, b = A^-1 c, V = A^-1
  M <- seq_along(ylist)
  A <- do.call(rbind, lapply(M, function(m) do.call(cbind, lapply(M, function(l) Oi[m, l] * crossprod(L[[m]], R[[l]])))))
  cc <- unlist(lapply(M, function(m) Reduce(`+`, lapply(M, function(l) Oi[m, l] * crossprod(L[[m]], ylist[[l]])))))
  V <- solve((A + t(A)) / 2); list(b = drop(V %*% cc), V = V)
}
three_sls <- function(ylist, Xlist, Z) {                              # common instruments: eqs. (8.39), (8.40), (8.45)
  E <- sapply(seq_along(ylist), function(m) ylist[[m]] - Xlist[[m]] %*% tsls(ylist[[m]], Xlist[[m]], Z, "homo")$b)
  Oi <- solve(crossprod(E) / nrow(Z))                                 # Omega_hat^-1 from equation-by-equation 2SLS
  c(sys_gls(ylist, lapply(Xlist, function(X) qr.fitted(qr(Z), X)), Xlist, Oi), list(Omega = solve(Oi)))
}
Z  <- cbind(const = 1, ad_budget = S1$ad_budget, coupon = S1$coupon, fee = S1$fee, log_pop = S1$log_pop)
XB <- cbind(const = 1, n_sellers = S1$n_sellers, ad_budget = S1$ad_budget, coupon = S1$coupon, log_pop = S1$log_pop)
XS <- cbind(const = 1, buyers = S1$buyers, fee = S1$fee, log_pop = S1$log_pop)
ci <- function(b, s) sprintf("[%.3f, %.3f]", b - qnorm(0.975) * s, b + qnorm(0.975) * s)
o  <- lm(S1$buyers ~ XB - 1); so <- sqrt(sandwich::vcovHC(o, type = "HC0")[2, 2])
tB <- tsls(S1$buyers, XB, Z); tBh <- tsls(S1$buyers, XB, Z, "homo"); tS <- tsls(S1$n_sellers, XS, Z)
th <- three_sls(list(S1$buyers, S1$n_sellers), list(XB, XS), Z)
sB <- sqrt(tB$V[2, 2]); s3 <- sqrt(th$V[2, 2])
cat(sprintf("OLS  %.4f (HC0 %.4f) %s\n", coef(o)[2], so, ci(coef(o)[2], so)))
cat(sprintf("2SLS %.4f (HC0 %.4f, homoskedastic %.4f) %s\n", tB$b[2], sB, sqrt(tBh$V[2, 2]), ci(tB$b[2], sB)))
cat(sprintf("3SLS %.4f (eq. 8.45 %.4f) %s; SE ratio to homoskedastic 2SLS %.3f\n", th$b[2], s3, ci(th$b[2], s3),
            s3 / sqrt(tBh$V[2, 2])))
cat(sprintf("seller equation: 2SLS %.4f (HC0 %.4f); max |3SLS - 2SLS| %.1e\n", tS$b[2], sqrt(tS$V[2, 2]),
            max(abs(th$b[6:9] - tS$b))))
eqs <- list(buyer = buyers ~ n_sellers + ad_budget + coupon + log_pop, seller = n_sellers ~ buyers + fee + log_pop)
ivB <- ivreg::ivreg(buyers ~ n_sellers + ad_budget + coupon + log_pop | ad_budget + coupon + fee + log_pop, data = S1)
sf3 <- systemfit::systemfit(eqs, method = "3SLS", inst = ~ ad_budget + coupon + fee + log_pop, data = S1,
                            methodResidCov = "noDfCor")
cat(sprintf("hand vs package: 2SLS %.1e, 3SLS %.1e\n", max(abs(tB$b - coef(ivB))), max(abs(th$b - coef(sf3)))))
# OLS  2.2678 (HC0 0.0148) [2.239, 2.297]
# 2SLS 1.1401 (HC0 0.1017, homoskedastic 0.1010) [0.941, 1.339]
# 3SLS 1.2123 (eq. 8.45 0.0862) [1.043, 1.381]; SE ratio to homoskedastic 2SLS 0.853
# seller equation: 2SLS 0.2186 (HC0 0.0178); max |3SLS - 2SLS| 2.4e-10
# hand vs package: 2SLS 9.2e-14, 3SLS 2.0e-09
```

The hand-written 2SLS differs from `ivreg` by $9.2 \times 10^{-14}$ and the 3SLS from `systemfit` by $2.0 \times 10^{-9}$, both mere floating-point error. The 3SLS and 2SLS of the seller equation differ by only $2.4 \times 10^{-10}$, which is Theorem 9 (ii): the overidentified seller equation gains nothing from the system. The diagnostics follow.

```r
WB <- Z[, c("const", "ad_budget", "coupon", "log_pop")]; WS <- Z[, c("const", "fee", "log_pop")]
fB <- first_stage_F(S1$n_sellers, Z[, "fee", drop = FALSE], WB)
fS <- first_stage_F(S1$buyers, Z[, c("ad_budget", "coupon")], WS)
cat(sprintf("first stage: buyer equation F_fs %.1f, F_eff %.1f, pi_fee %.4f; seller equation F_fs %.2f, F_eff %.2f\n",
            fB$F, eff_F(S1$n_sellers, Z[, "fee", drop = FALSE], WB), fB$pi, fS$F,
            eff_F(S1$buyers, Z[, c("ad_budget", "coupon")], WS)))
cf <- cf_test(S1$buyers, XB, Z, endog = "n_sellers")                  # y on (X, v_hat): eq. (6.15)
J <- hansen_j(S1$n_sellers, XS, Z); sg <- sargan(S1$n_sellers, XS, Z)
cat(sprintf("control function: rho %.4f (HC0 %.4f), t %.2f; seller equation: J %.3f (p %.3f), Sargan %.3f (p %.3f)\n",
            cf$rho, cf$se, cf$t, J$J, J$p, sg$stat, sg$p))
# first stage: buyer equation F_fs 174.7, F_eff 166.6, pi_fee -0.1378; seller equation F_fs 57.39, F_eff 56.76
# control function: rho 1.2595 (HC0 0.0354), t 35.53; seller equation: J 1.690 (p 0.194), Sargan 1.904 (p 0.168)
```

The instruments of both equations are strong: the first-stage F of the buyer equation is 174.7, with a sample fee coefficient of $-0.138$ close to the population value $-0.143$; that of the seller equation is 57.39, above 19.93, the strictest Stock-Yogo critical value for two instruments, and its effective F is 56.76. The control function's $t$ value of 35.53 shows that the number of sellers is clearly endogenous, and the single overidentifying restriction of the seller equation is not rejected.

### 6.3 Block S1m: One Misspecified Equation

`fit_S1` computes, for one sample, both the 2SLS of the buyer equation and the 3SLS of the two equations; the same function is applied first to the S1m sample, then to two Monte Carlo runs of 1,000 replications each, and finally to one million markets to compare the two city-program designs.

```r
fit_S1 <- function(d) {                                               # buyer equation: 2SLS (HC0) and 3SLS (eq. 8.45)
  Z  <- cbind(1, d$ad_budget, d$coupon, d$fee, d$log_pop)
  XB <- cbind(1, d$n_sellers, d$ad_budget, d$coupon, d$log_pop); XS <- cbind(1, d$buyers, d$fee, d$log_pop)
  t2 <- tsls(d$buyers, XB, Z); t3 <- three_sls(list(d$buyers, d$n_sellers), list(XB, XS), Z)
  c(tsls = t2$b[[2]], se_tsls = sqrt(t2$V[2, 2]), three = t3$b[[2]], se_three = sqrt(t3$V[2, 2]),
    coupon_tsls = t2$b[[4]], coupon_three = t3$b[[4]])
}
f <- fit_S1(S1m)
Jm <- hansen_j(S1m$n_sellers, cbind(1, S1m$buyers, S1m$fee, S1m$log_pop),
               cbind(1, S1m$ad_budget, S1m$coupon, S1m$fee, S1m$log_pop))
cat(sprintf("S1m: 2SLS %.4f %s; 3SLS %.4f %s; seller equation J %.2f (p %.1e)\n", f[["tsls"]],
            ci(f[["tsls"]], f[["se_tsls"]]), f[["three"]], ci(f[["three"]], f[["se_three"]]), Jm$J, Jm$p))
mc <- function(kappa, seed, R = 1000) { set.seed(seed); t(replicate(R, fit_S1(gen_S1(P1, kappa = kappa)))) }
cover <- function(est, se) mean(abs(est - 1.2) <= qnorm(0.975) * se)
for (k in list(list("S1", 0, 9121), list("S1m", -0.0045, 9122))) {
  s <- mc(k[[2]], k[[3]])
  cat(sprintf("MC %-3s 2SLS mean %.4f sd %.4f cover %.3f | 3SLS mean %.4f sd %.4f cover %.3f\n", k[[1]],
              mean(s[, "tsls"]), sd(s[, "tsls"]), cover(s[, "tsls"], s[, "se_tsls"]),
              mean(s[, "three"]), sd(s[, "three"]), cover(s[, "three"], s[, "se_three"])))
}
for (r in c(0, -0.7)) {                                               # one million markets per design
  set.seed(if (r == 0) 9113 else 9114)
  g <- fit_S1(gen_S1(modifyList(P1, list(r = r)), n = 1e6, kappa = -0.0045))
  cat(sprintf("corr(coupon, fee) = %4.1f: n_sellers 2SLS %.4f, 3SLS %.4f | coupon 2SLS %.4f, 3SLS %.4f\n",
              r, g[["tsls"]], g[["three"]], g[["coupon_tsls"]], g[["coupon_three"]]))
}
# S1m: 2SLS 1.2120 [1.033, 1.390]; 3SLS 0.8131 [0.666, 0.960]; seller equation J 67.89 (p 1.7e-16)
# MC S1  2SLS mean 1.1954 sd 0.0961 cover 0.942 | 3SLS mean 1.1972 sd 0.0804 cover 0.932
# MC S1m 2SLS mean 1.1909 sd 0.0887 cover 0.957 | 3SLS mean 0.7925 sd 0.1236 cover 0.001
# corr(coupon, fee) =  0.0: n_sellers 2SLS 1.2017, 3SLS 1.2013 | coupon 2SLS 0.4026, 3SLS 0.7320
# corr(coupon, fee) = -0.7: n_sellers 2SLS 1.1940, 3SLS 0.7967 | coupon 2SLS 0.4045, 3SLS 0.7384
```

The S1 results for the two estimators match Section 3.7: both estimators average about 1.2, the 3SLS sampling standard deviation of 0.080 is smaller than the 2SLS one of 0.096, and the coverage rates of 0.932 and 0.942 are both close to the nominal level. In the S1m sample, the 3SLS interval $[0.666, 0.960]$ excludes the truth while the 2SLS interval $[1.033, 1.390]$ covers it; the Monte Carlo shows that this is no accident: the 3SLS average falls 0.407 below the truth. The two city-program designs on one million markets are discussed in Section 3.8.

### 6.4 Block S2: SUR

```r
P2 <- list(n = 1500, ship0 = 35, ship_sd = 5, prem0 = 0.30, prem_ship = 0.008, prem_sd = 0.06, loy_ship = 0.08,
           loy_cut = 0.25, b1 = c(18, 0.50, 25), b2 = c(3.6, -0.04, 0.5), s1 = 4, s2 = 0.45, rho = 0.6)
gen_S2 <- function(par, n = par$n) {
  free_ship_min <- par$ship0 + par$ship_sd * rnorm(n)                 # free-delivery minimum, USD
  premium_share <- par$prem0 + par$prem_ship * (free_ship_min - par$ship0) + par$prem_sd * rnorm(n)
  loyalty <- as.integer(par$loy_ship * (free_ship_min - par$ship0) + rnorm(n) > par$loy_cut)
  Om <- matrix(c(par$s1^2, par$rho * par$s1 * par$s2, par$rho * par$s1 * par$s2, par$s2^2), 2)
  e <- matrix(rnorm(2 * n), n) %*% chol(Om)                           # local income shocks raise both outcomes
  data.frame(free_ship_min, premium_share, loyalty,
             basket_usd = par$b1[1] + par$b1[2] * free_ship_min + par$b1[3] * premium_share + e[, 1],
             orders_pb = par$b2[1] + par$b2[2] * free_ship_min + par$b2[3] * loyalty + e[, 2])
}
set.seed(9103); S2 <- gen_S2(P2)
sur_fgls <- function(ylist, Xlist) {                                  # eqs. (7.40), (7.45) and (7.54)
  E <- sapply(seq_along(ylist), function(m) ylist[[m]] - qr.fitted(qr(Xlist[[m]]), ylist[[m]]))
  Oi <- solve(crossprod(E) / length(ylist[[1]]))                      # Omega_hat from equation-by-equation OLS
  c(sys_gls(ylist, Xlist, Xlist, Oi), list(Omega = solve(Oi)))
}
y  <- list(S2$basket_usd, S2$orders_pb)
Xs <- cbind(1, S2$free_ship_min, S2$premium_share, S2$loyalty)       # the same regressors in both equations
Xd <- list(Xs[, 1:3], Xs[, c(1, 2, 4)])                               # each equation its own regressors
same <- sur_fgls(y, list(Xs, Xs)); dif <- sur_fgls(y, Xd)
ols <- c(qr.coef(qr(Xs), y[[1]]), qr.coef(qr(Xs), y[[2]]))
se_ols <- sqrt(c(dif$Omega[1, 1] * diag(solve(crossprod(Xd[[1]]))), dif$Omega[2, 2] * diag(solve(crossprod(Xd[[2]])))))
cat(sprintf("same regressors: max |SUR - OLS| %.1e\n", max(abs(same$b - ols))))
cat("different regressors, SE of SUR / SE of OLS:", sprintf("%.4f", sqrt(diag(dif$V)) / se_ols), "\n")
w <- c(0, mean(S2$orders_pb), 0, 0, 0, mean(S2$basket_usd), 0, 0)   # sample means, treated as fixed
Vb <- same$V; Vb[1:4, 5:8] <- 0; Vb[5:8, 1:4] <- 0                   # as if the two equations were unrelated
est <- sum(w * same$b); s_full <- sqrt(drop(t(w) %*% same$V %*% w)); s_blk <- sqrt(drop(t(w) %*% Vb %*% w))
cat(sprintf("revenue per buyer, minimum USD 1 higher: %.4f (SE %.4f; without the cross-equation covariance %.4f) %s\n",
            est, s_full, s_blk, ci(est, s_full)))
# same regressors: max |SUR - OLS| 4.3e-12
# different regressors, SE of SUR / SE of OLS: 0.9997 0.9439 0.7964 0.9930 0.9830 0.7964 
# revenue per buyer, minimum USD 1 higher: -0.4296 (SE 0.1694; without the cross-equation covariance 0.1388) [-0.762, -0.098]
```

The six standard-error ratios are, in order, for the constant, the minimum, and the premium-seller share of the basket equation, and for the constant, the minimum, and the loyalty program of the orders equation; the OLS standard errors use the diagonal elements of the same $\hat\Omega$, so the comparison reflects only the difference in weighting. The constants hardly improve (0.9997 and 0.993), and the two variables that appear in only one equation improve the most. The revenue-effect line confirms Theorem 7: with identical regressors the point estimate is the same as OLS, yet the cross-equation covariance changes the standard error from 0.139 to 0.169; both treat the sample means as fixed weights.

### 6.5 Block S3: Weak Instruments

`fee_sd_for` uses the reduced form to compute the standard deviation of the onboarding fee that makes the concentration parameter equal a given value; `weak_stats` returns, for one sample, the 2SLS estimate, whether its robust Wald interval covers the truth, the first-stage F, and the shape of the AR set.

```r
Om <- matrix(c(P1$sB^2, P1$rho * P1$sB * P1$sS, P1$rho * P1$sB * P1$sS, P1$sS^2), 2)
s2v <- (t(solve(Gam)) %*% Om %*% solve(Gam))[2, 2]                   # reduced-form error variance of n_sellers
fee_sd_for <- function(mu2, n = 1500) sqrt(mu2 * s2v / (n * Pi["fee", "n_sellers"]^2))
P3 <- modifyList(P1, list(r = 0, fee_sd = fee_sd_for(1)))            # concentration parameter 1
set.seed(9104); S3 <- gen_S1(P3)
W3 <- cbind(1, S3$ad_budget, S3$coupon, S3$log_pop); zf <- cbind(S3$fee)
t3 <- tsls(S3$buyers, cbind(1, S3$n_sellers, W3[, -1]), cbind(W3, zf))
f3 <- first_stage_F(S3$n_sellers, zf, W3); ar <- ar_ci(S3$buyers, S3$n_sellers, zf, W3)
a12 <- ar_test(S3$buyers, S3$n_sellers, zf, W3, 1.2); a0 <- ar_test(S3$buyers, S3$n_sellers, zf, W3, 0)
cat(sprintf("fee sd USD %.3f; F_fs %.4f, F_eff %.4f; 2SLS %.3f (HC0 %.2f) %s\n", P3$fee_sd, f3$F, eff_F(S3$n_sellers, zf, W3),
            t3$b[2], sqrt(t3$V[2, 2]), ci(t3$b[2], sqrt(t3$V[2, 2]))))
cat(sprintf("AR set: %s, (-Inf, %.3f] and [%.3f, Inf); AR(1.2) %.3f (p %.3f), AR(0) %.3f (p %.3f)\n",
            ar$type, ar$lo, ar$hi, a12$stat, a12$p, a0$stat, a0$p))
weak_stats <- function(d, b = 1.2) {                                  # one sample: 2SLS, its HC0 Wald interval, F_fs, AR set
  W <- cbind(1, d$ad_budget, d$coupon, d$log_pop); rs <- function(v) v - qr.fitted(qr(W), v)
  y <- rs(d$buyers); x <- rs(d$n_sellers); z <- rs(d$fee); df2 <- length(y) - ncol(W) - 1
  est <- sum(z * y) / sum(z * x); e <- y - est * x
  se <- sqrt(sum(z^2 * e^2)) / abs(sum(z * x))                        # HC0, eq. (5.34)
  p <- sum(z * x) / sum(z^2); Ffs <- p^2 * sum(z^2) / (sum((x - p * z)^2) / df2)
  a <- ar_ci(d$buyers, d$n_sellers, cbind(d$fee), W)
  c(est = est, F = Ffs, wald = abs(est - b) <= qnorm(0.975) * se, ar = ar_covers(a, b),
    type = match(a$type, c("bounded", "two_rays", "real_line")))
}
grid <- c(0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128)
for (j in seq_along(grid)) {
  set.seed(if (grid[j] == 1) 9124 else 94000 + j)
  s <- t(replicate(if (grid[j] == 1) 2000 else 1000, weak_stats(gen_S1(modifyList(P3, list(fee_sd = fee_sd_for(grid[j])))))))
  cat(sprintf("mu2 %6.2f: mean F %6.2f  Wald %.4f  AR %.4f  median 2SLS %.3f  AR bounded %.3f\n", grid[j],
              mean(s[, "F"]), mean(s[, "wald"]), mean(s[, "ar"]), median(s[, "est"]), mean(s[, "type"] == 1)))
}
# fee sd USD 0.349; F_fs 0.0195, F_eff 0.0182; 2SLS -4.795 (HC0 53.50) [-109.662, 100.071]
# AR set: two_rays, (-Inf, 2.191] and [2.718, Inf); AR(1.2) 0.417 (p 0.518), AR(0) 0.074 (p 0.785)
# mu2   0.25: mean F   1.28  Wald 0.7250  AR 0.9390  median 2SLS 2.096  AR bounded 0.080
# mu2   0.50: mean F   1.51  Wald 0.7660  AR 0.9510  median 2SLS 1.860  AR bounded 0.116
# mu2   1.00: mean F   1.94  Wald 0.8245  AR 0.9530  median 2SLS 1.544  AR bounded 0.166
# mu2   2.00: mean F   3.10  Wald 0.8740  AR 0.9370  median 2SLS 1.376  AR bounded 0.291
# mu2   4.00: mean F   5.12  Wald 0.8780  AR 0.9490  median 2SLS 1.232  AR bounded 0.524
# mu2   8.00: mean F   8.98  Wald 0.9170  AR 0.9430  median 2SLS 1.210  AR bounded 0.795
# mu2  16.00: mean F  17.30  Wald 0.9180  AR 0.9480  median 2SLS 1.220  AR bounded 0.980
# mu2  32.00: mean F  32.87  Wald 0.9500  AR 0.9680  median 2SLS 1.201  AR bounded 1.000
# mu2  64.00: mean F  65.45  Wald 0.9340  AR 0.9300  median 2SLS 1.203  AR bounded 1.000
# mu2 128.00: mean F 129.22  Wald 0.9430  AR 0.9380  median 2SLS 1.202  AR bounded 1.000
```

![The Monte Carlo of block S3. The horizontal axis is the average first-stage F in each design (log scale, with a dash-dotted line at 10), and the concentration parameter of the ten designs runs from 0.25 to 128. The left panel shows the coverage rates of the nominal 95% Wald interval (HC0) and of the AR set, with a horizontal solid line at 0.95 and two dotted lines above and below it marking the Monte Carlo error band for 1,000 replications; the right panel shows the median of the 2SLS estimates, with the interquartile range shaded, and the two dashed lines are the truth, 1.2, and the probability limit of OLS.](assets/fig/fig_09_weak_iv.svg)

The average F is about $1 + \mu^2$, as Section 2.2 said. The coverage of the Wald interval rises as the first stage strengthens: at $\mu^2 = 0.25$ it is only 0.725, in the two designs with average F of 8.98 and 17.30 it is 0.917 and 0.918, and only at $\mu^2 = 32$ does it return to 0.950. The coverage of the AR set lies between 0.930 and 0.968 in all ten designs, whatever the strength; the Monte Carlo standard error of a coverage rate is about 0.0069 with 1,000 replications and about 0.0049 with 2,000. The median of 2SLS falls all the way from 2.096 to near the truth, hugging the OLS probability limit at the weak end. The share of bounded AR sets rises with strength, to about half (0.524) at $\mu^2 = 4$, and all sets are bounded once $\mu^2$ is at least 32: whether the set is bounded is itself a report on the strength of the instruments.

### 6.6 Block S4: Overidentification Tests

`gen_S4` generates the three randomized seller-side programs, with the error standard deviation varying exponentially with the referral bonus; `overid` computes 2SLS, $J$, and Sargan for one sample and can be told to use only some of the instruments.

```r
P4 <- modifyList(P1, list(aS = 9.2, pay_lo = 1, pay_hi = 7, ref_lo = 20, ref_hi = 80, dpay = -0.30, dref = 0.02,
                          eta = 0.5, theta = 0.05))
gen_S4 <- function(par, n = par$n) {
  log_pop <- rnorm(n); ad_budget <- par$ad0 + par$ad_sd * rnorm(n); coupon <- par$cp0 + par$cp_sd * rnorm(n)
  fee <- par$fee0 + par$fee_sd * rnorm(n)                             # three seller-side programs, randomized
  payout_days <- runif(n, par$pay_lo, par$pay_hi); referral_bonus <- runif(n, par$ref_lo, par$ref_hi)
  e <- matrix(rnorm(2 * n), n) %*% chol(matrix(c(par$sB^2, rep(par$rho * par$sB * par$sS, 2), par$sS^2), 2))
  zr <- (referral_bonus - 50) / (60 / sqrt(12)); t2 <- 2 * par$eta * sqrt(3)
  h <- exp(par$eta * zr) / sqrt(sinh(t2) / t2)                        # error scale rises with the referral bonus
  cB <- par$aB + par$gad * ad_budget + par$gcp * coupon + par$lB * log_pop + par$theta * (referral_bonus - 50) + e[, 1] * h
  cS <- par$aS + par$dfee * fee + par$dpay * payout_days + par$dref * referral_bonus + par$lS * log_pop + e[, 2] * h
  buyers <- (cB + par$bBS * cS) / (1 - par$bBS * par$bSB)
  data.frame(buyers, n_sellers = cS + par$bSB * buyers, ad_budget, coupon, log_pop, fee, payout_days, referral_bonus)
}
overid <- function(d, inst = c("fee", "payout_days", "referral_bonus")) {
  W <- cbind(1, d$ad_budget, d$coupon, d$log_pop); X <- cbind(1, d$n_sellers, W[, -1]); Z <- cbind(W, as.matrix(d[inst]))
  t2 <- tsls(d$buyers, X, Z); J <- hansen_j(d$buyers, X, Z); sg <- sargan(d$buyers, X, Z)
  c(est = t2$b[[2]], se = sqrt(t2$V[2, 2]), J = J$J, J_p = J$p, sargan = sg$stat, sargan_p = sg$p)
}
set.seed(9105); S4 <- gen_S4(P4); set.seed(9106); S4v <- gen_S4(modifyList(P4, list(theta = 0)))
for (nm in c("S4v", "S4")) {
  a <- overid(get(nm)); a2 <- overid(get(nm), c("fee", "payout_days"))
  cat(sprintf("%-3s 2SLS %.4f %s, J %.3f (p %.3g), Sargan %.3f (p %.3g) | fee and payout only: 2SLS %.4f %s, J %.3f (p %.3f)\n",
              nm, a[["est"]], ci(a[["est"]], a[["se"]]), a[["J"]], a[["J_p"]], a[["sargan"]], a[["sargan_p"]],
              a2[["est"]], ci(a2[["est"]], a2[["se"]]), a2[["J"]], a2[["J_p"]]))
}
c05 <- qchisq(0.95, 2); c01 <- qchisq(0.99, 2)
set.seed(9125); sv <- t(replicate(2000, overid(gen_S4(modifyList(P4, list(theta = 0))))))
set.seed(9126); si <- t(replicate(1000, overid(gen_S4(P4))))
cat(sprintf("valid (2000 samples): reject at 5%%: J %.4f, Sargan %.4f; at 1%%: J %.4f, Sargan %.4f | invalid (1000): J %.3f, mean 2SLS %.4f\n",
            mean(sv[, "J"] > c05), mean(sv[, "sargan"] > c05), mean(sv[, "J"] > c01), mean(sv[, "sargan"] > c01),
            mean(si[, "J"] > c05), mean(si[, "est"])))
# S4v 2SLS 1.1998 [1.118, 1.282], J 0.445 (p 0.8), Sargan 0.514 (p 0.774) | fee and payout only: 2SLS 1.2018 [1.115, 1.288], J 0.443 (p 0.506)
# S4  2SLS 1.5117 [1.435, 1.589], J 140.379 (p 3.29e-31), Sargan 186.218 (p 3.66e-41) | fee and payout only: 2SLS 1.1691 [1.050, 1.289], J 0.534 (p 0.465)
# valid (2000 samples): reject at 5%: J 0.0535, Sargan 0.0720; at 1%: J 0.0120, Sargan 0.0180 | invalid (1000): J 1.000, mean 2SLS 1.5185
```

In S4v all three instruments are valid: 2SLS gives 1.200, on target, and $J = 0.445$ ($p = 0.800$). In S4 the referral bonus violates the exclusion restriction: 2SLS drifts to 1.512, a bias of 0.312, and the interval $[1.435, 1.589]$ excludes the truth; $J = 140.38$ (the 5% critical value with two degrees of freedom is 5.99), an unambiguous rejection. Drop the referral bonus and use only the onboarding fee and payout days, and 2SLS returns to 1.169, with $J = 0.534$ ($p = 0.465$) on the one remaining degree of freedom. In repeated sampling, $J$ rejects in all 1,000 invalid samples, and the 2SLS average of 1.518 is 0.318 above the truth; in the 2,000 valid samples, the nominal 5% $J$ test rejects 5.35% of the time and the homoskedastic Sargan test 7.20%, and at the nominal 1% level the two reject 1.2% and 1.8%. Because the errors vary with the referral bonus, Sargan's denominator misestimates $S_0$, and the test rejects too often.

### 6.7 The Grand Reconciliation

| Block | Estimator | Estimate | 95% interval or set | True value | Covered |
|---|---|---|---|---|---|
| S1 | OLS, buyer equation | 2.268 | $[2.239, 2.297]$ | 1.2 | No |
| S1 | 2SLS, buyer equation | 1.140 | $[0.941, 1.339]$ | 1.2 | Yes |
| S1 | 3SLS, buyer equation | 1.212 | $[1.043, 1.381]$ | 1.2 | Yes |
| S1 | 2SLS, seller equation | 0.219 | $[0.184, 0.254]$ | 0.25 | Yes |
| S1m | 2SLS, buyer equation | 1.212 | $[1.033, 1.390]$ | 1.2 | Yes |
| S1m | 3SLS, buyer equation | 0.813 | $[0.666, 0.960]$ | 1.2 | No |
| S2 | SUR, revenue effect (dollars) | $-0.430$ | $[-0.762, -0.098]$ | $-0.518$ | Yes |
| S3 | 2SLS, robust Wald interval | $-4.80$ | $[-109.66, 100.07]$ | 1.2 | Yes |
| S3 | AR set |  | $(-\infty, 2.191] \cup [2.718, \infty)$ | 1.2 | Yes |
| S4v | 2SLS, three instruments | 1.200 | $[1.118, 1.282]$ | 1.2 | Yes |
| S4 | 2SLS, three instruments | 1.512 | $[1.435, 1.589]$ | 1.2 | No |
| S4 | 2SLS, referral bonus dropped | 1.169 | $[1.050, 1.289]$ | 1.2 | Yes |

The three intervals in the table that miss the truth all do so by design: OLS has simultaneity bias, the 3SLS of S1m is infected by the seller equation, and S4 uses the referral bonus, which violates the exclusion restriction. The Wald interval of S3 does cover the truth, but it is too wide to carry any information; conclusions about coverage have to be read from the Monte Carlo: with $\mu^2 = 1$, the Wald interval covers 82.45% and the AR set 95.3%.

The reconciliation of this section comes down to one sentence: every estimator aims at the quantity defined by its moment conditions, and when things go wrong, either the moment conditions themselves fail or the sampling distribution under weak instruments has been taken for a normal one.

## 7 Failure Modes and Robustness

A failure of the exclusion restriction is the most fundamental threat to IV, and also the hardest to detect. Under exact identification it is completely untestable: the $\Pi_{(m)}$ of the buyer equation of S1 has a single row, and whether the onboarding fee has a direct effect is a question the data cannot answer. Under overidentification, $J$ raises the alarm only when the instruments disagree with one another; if several instruments go wrong in the same direction, it stays silent. There are three practical responses. The first is to treat the exclusion restriction as an economic assumption that needs defending, spelling out the channel, and the only channel, through which the instrument affects the outcome. The second is to look at the reduced form: when the first stage is not weak and the instrument has no effect on the outcome, the structural effect cannot be large either. The third, under overidentification, is to drop instruments one at a time and reestimate: when the referral bonus is dropped in S4 of Section 6.6, 2SLS returns to near the truth and $J$ no longer rejects.

Weak instruments are the most common failure, and the danger is that they look perfectly normal. When the first-stage F is below the rule-of-thumb threshold, 2SLS is biased toward OLS and the Wald interval systematically misses the truth; in the Monte Carlo of S3, coverage is far below the nominal level. The response is to report the first stage up front, using the effective F under heteroskedasticity or clustering, and with a single instrument to report the AR set directly and accept that it may be unbounded. Do not pick specifications by whether F exceeds 10 and then report conventional intervals: Lee et al. point out that inference done this way is still distorted, and Andrews, Stock and Sun (2019) discuss the distortions that come from selectively dropping specifications. Young (2022) replicates 1,309 regressions estimated by 2SLS in 30 papers in the journals of the American Economic Association and finds that when errors are not independent and identically distributed and a few observations have high leverage, F-based weak-instrument pretests carry almost no information about the actual rejection rates and biases of tests, and that IV tests also have low power. With many weak instruments things are worse: the 28 instruments of Section 5.1 pushed 2SLS onto OLS. The many-moment bias of Section 3.5 of Chapter 7 and instrument proliferation in the dynamic panels of Chapter 11 are the same problem.

Heteroskedastic errors throw homoskedastic tests off. In S4v the errors vary with one of the instruments, the homoskedastic Sargan test rejects too often, and the robust $J$ stays close to the nominal level. The rule is simple: use $J$ for overidentification tests, the robust $t$ for endogeneity tests, equation (5.34) for the standard errors of 2SLS, and the sandwich or the minimum chi-square estimator for system estimation.

The risk of system estimation is contagion. The efficiency of 3SLS comes from the overidentifying restrictions of the other equations, and the same channel also brings in their specification errors, which is why the 3SLS of S1m lands far from the truth. There are two lines of defense. The first is to report the $J$ of each equation separately: the $J$ of the seller equation of S1m is far above the critical value. The second is to set 3SLS side by side with equation-by-equation 2SLS: the two have the same probability limit under correct specification, and in S1m their difference far exceeds any standard error. When only one equation matters, 2SLS is the more robust choice; if one equation in the system cannot be trusted, it is better to replace it with an unrestricted reduced form or to leave it out of the system.

Heterogeneous effects change what 2SLS estimates. Bramble's cross-side effect is the same in every market; on a real platform, new sellers may attract buyers in big cities and in small towns to different degrees. Then, under the two conditions of Section 3.2 (the instrument pushes every market in the same direction, and the effects are unrelated to the instrument), 2SLS estimates a weighted average whose weights are set by how far the instrument moves the number of sellers in each market; different instruments aim at different averages, and $J$ may reject for that reason even though each instrument is valid in its own sense. Chapter 18 writes this weighted average as a LATE; with several instruments, the same-direction condition must hold jointly over the values of all the instruments (Section 4.1 of Chapter 18).

String these failure modes together, and the robustness of this chapter's methods comes down, in the end, to three things: whether the instruments affect the outcome only through the endogenous regressors, whether the first stage is strong enough, and whether every equation in the system is correctly specified. The first rests on economic argument and overidentification tests, the second on first-stage diagnostics and robust inference of the AR kind, and the third on the equation-by-equation $J$ and the comparison of 2SLS with 3SLS. Efficiency always comes last.

## 8 Further Reading

::: {.readings}
Required reading, in suggested reading order:

- Wooldridge (2010, *Econometric Analysis of Cross Section and Panel Data*, 2nd ed.), §5.1 to §5.2, §6.2 to §6.3, §7.3 to §7.7, §8.2 to §8.5, and §9.1 to §9.4. A unified treatment of 2SLS, tests of endogeneity and overidentification, SUR and FGLS, system IV, and simultaneous equations; the source of Theorems 1 to 4 and 6 to 9 of this chapter.
- Staiger and Stock (1997, Econometrica) and Stock and Yogo (2005). The local-to-zero asymptotic theory of weak instruments and the source of the critical values for the "first-stage F"; this chapter follows the NBER technical working paper versions of both.
- Andrews, Stock and Sun (2019, Annual Review of Economics). A survey of the weak-instruments literature, focusing on the heteroskedastic and clustered cases, with practical recommendations.
- Lee, McCrary, Moreira and Porter (2022, American Economic Review). The tF critical values for a single instrument and their implications for published research.

Further reading:

- Working (1927, Quarterly Journal of Economics) and Haavelmo (1943, Econometrica). The identification problem of the demand curve and the probabilistic formulation of simultaneous equations, the two starting points of the simultaneous-equations tradition in econometrics; Imbens (2014, Statistical Science) §3.5 gives an accessible restatement of Working's problem.
- Anderson and Rubin (1949, Annals of Mathematical Statistics). The source of the AR test and confidence set.
- Sargan (1958, Econometrica) and Hansen (1982, Econometrica). The two sources of overidentification tests.
- Durbin (1954), Wu (1973, Econometrica), and Hausman (1978, Econometrica). Tests of endogeneity and the testing principle of comparing two estimators.
- Zellner (1962, Journal of the American Statistical Association) and Zellner and Theil (1962, Econometrica). The original papers on SUR and 3SLS.
- Montiel Olea and Pflueger (2013, Journal of Business & Economic Statistics). The effective F under heteroskedasticity and its critical values.
- Angrist and Krueger (1991, Quarterly Journal of Economics) and Bound, Jaeger and Baker (1995, Journal of the American Statistical Association). The quarter-of-birth instrument and the weak-instrument problem it exposed.
- Young (2022, European Economic Review). The quality of inference in published IV studies.
- Chu and Manchanda (2016, Marketing Science) and Rysman (2004, Review of Economic Studies). Simultaneous-equations estimation of cross-side network effects in two-sided markets.
- Petrin and Train (2010, Journal of Marketing Research) and Rossi (2014, Marketing Science). The control function in the marketing literature, and a critique of the quality of IV in marketing applications.
- Hansen (2014) §15.1 to §15.7. Another textbook treatment of IV, including Bekker asymptotics and identification failure.
:::

::: {.apa-refs}
- Anderson, T. W., & Rubin, H. (1949). Estimation of the parameters of a single equation in a complete system of stochastic equations. *The Annals of Mathematical Statistics, 20*(1), 46-63. https://doi.org/10.1214/aoms/1177730090
- Andrews, I., Stock, J. H., & Sun, L. (2019). Weak instruments in instrumental variables regression: Theory and practice. *Annual Review of Economics, 11*(1), 727-753. https://doi.org/10.1146/annurev-economics-080218-025643
- Angrist, J. D., & Krueger, A. B. (1991). Does compulsory school attendance affect schooling and earnings? *The Quarterly Journal of Economics, 106*(4), 979-1014. https://doi.org/10.2307/2937954
- Bound, J., Jaeger, D. A., & Baker, R. M. (1995). Problems with instrumental variables estimation when the correlation between the instruments and the endogenous explanatory variable is weak. *Journal of the American Statistical Association, 90*(430), 443-450. https://doi.org/10.1080/01621459.1995.10476536
- Chu, J., & Manchanda, P. (2016). Quantifying cross and direct network effects in online consumer-to-consumer platforms. *Marketing Science, 35*(6), 870-893. https://doi.org/10.1287/mksc.2016.0976
- Durbin, J. (1954). Errors in variables. *Review of the International Statistical Institute, 22*(1/3), 23-32. https://doi.org/10.2307/1401917
- Haavelmo, T. (1943). The statistical implications of a system of simultaneous equations. *Econometrica, 11*(1), 1-12. https://doi.org/10.2307/1905714
- Hansen, B. E. (2014). *Econometrics* [Unpublished manuscript]. Department of Economics, University of Wisconsin.
- Hansen, L. P. (1982). Large sample properties of generalized method of moments estimators. *Econometrica, 50*(4), 1029-1054. https://doi.org/10.2307/1912775
- Hausman, J. A. (1978). Specification tests in econometrics. *Econometrica, 46*(6), 1251-1271. https://doi.org/10.2307/1913827
- Imbens, G. W. (2014). Instrumental variables: An econometrician's perspective. *Statistical Science, 29*(3), 323-358. https://doi.org/10.1214/14-sts480
- Lee, D. S., McCrary, J., Moreira, M. J., & Porter, J. (2022). Valid t-ratio inference for IV. *American Economic Review, 112*(10), 3260-3290. https://doi.org/10.1257/aer.20211063
- Montiel Olea, J. L., & Pflueger, C. (2013). A robust test for weak instruments. *Journal of Business & Economic Statistics, 31*(3), 358-369. https://doi.org/10.1080/00401706.2013.806694
- Petrin, A., & Train, K. (2010). A control function approach to endogeneity in consumer choice models. *Journal of Marketing Research, 47*(1), 3-13. https://doi.org/10.1509/jmkr.47.1.3
- Rossi, P. E. (2014). Invited paper: Even the rich can make themselves poor: A critical examination of IV methods in marketing applications. *Marketing Science, 33*(5), 655-672. https://doi.org/10.1287/mksc.2014.0860
- Rysman, M. (2004). Competition between networks: A study of the market for Yellow Pages. *The Review of Economic Studies, 71*(2), 483-512. https://doi.org/10.1111/0034-6527.00512
- Sargan, J. D. (1958). The estimation of economic relationships using instrumental variables. *Econometrica, 26*(3), 393-415. https://doi.org/10.2307/1907619
- Staiger, D., & Stock, J. H. (1997). Instrumental variables regression with weak instruments. *Econometrica, 65*(3), 557-586. https://doi.org/10.2307/2171753
- Stock, J. H., & Yogo, M. (2005). Testing for weak instruments in linear IV regression. In D. W. K. Andrews & J. H. Stock (Eds.), *Identification and inference for econometric models: Essays in honor of Thomas Rothenberg* (pp. 80-108). Cambridge University Press. https://doi.org/10.1017/cbo9780511614491.006
- Wooldridge, J. M. (2010). *Econometric analysis of cross section and panel data* (2nd ed.). MIT Press.
- Working, E. J. (1927). What do statistical "demand curves" show? *The Quarterly Journal of Economics, 41*(2), 212-235. https://doi.org/10.2307/1883501
- Wu, D.-M. (1973). Alternative tests of independence between stochastic regressors and disturbances. *Econometrica, 41*(4), 733-750. https://doi.org/10.2307/1914093
- Young, A. (2022). Consistency without inference: Instrumental variables in practical application. *European Economic Review, 147*, Article 104112. https://doi.org/10.1016/j.euroecorev.2022.104112
- Zellner, A. (1962). An efficient method of estimating seemingly unrelated regressions and tests for aggregation bias. *Journal of the American Statistical Association, 57*(298), 348-368. https://doi.org/10.1080/01621459.1962.10480664
- Zellner, A., & Theil, H. (1962). Three-stage least squares: Simultaneous estimation of simultaneous equations. *Econometrica, 30*(1), 54-78. https://doi.org/10.2307/1911287
:::
