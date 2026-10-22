---
title: "Quantile Regression"
subtitle: "Beyond the Conditional Mean"
seriesline: "Foundations of Information Systems Economics · Chapter 13"
series: "Foundations of Information Systems Economics: Empirical Methods for Information Systems Economics"
sidebartitle: "Foundations of Information Systems Economics"
sidebarsub: "Theory · Methods · Applications"
chapterlabel: "Chapter 13 · Quantile Regression"
---

## Introduction

Marigold is a fictional e-commerce platform that hosts online shops for independent merchants. The platform has built an automatic pricing tool, AutoPrice, that raises prices when demand is strong and cuts them when demand is weak. To evaluate it, the platform picks half of its merchants at random and switches the tool on for them, and one month later it compares the two groups' monthly revenue. The two groups' average revenue is almost identical, and the first report writes itself: AutoPrice does nothing. But line up each group's merchants from lowest to highest revenue and compare them at the same position, and the picture changes completely: at the bottom of the line, the merchants with the tool earn less; at the top, they earn more. A tool whose average effect is zero may be quietly moving revenue from small merchants to large ones.

This chapter covers quantile regression (QR) and a set of causal tools built around it. A mean regression answers "when $X$ changes by one unit, by how much does the conditional mean of $Y$ move?"; QR answers "when $X$ changes by one unit, by how much does the conditional distribution of $Y$ move at its $\tau$-th position?" ($\tau \in (0,1)$ is the quantile level; Section 2.1 gives the details). On the surface this merely swaps the mean squared error of Chapter 5 for another loss, but underneath lie three difficulties that never arise with means. First, quantiles are not additive: the difference between two distributions at the same position is not a quantile of the individual differences, and a shift in a distribution is not anyone's effect. Second, conditional quantiles do not average back to the unconditional quantile: an effect obtained after controlling for covariates cannot be averaged back to the population by weighting with the distribution of the covariates, as a mean effect can. Third, in a panel the merchant fixed effects cannot be removed by subtracting each merchant's mean, as they can in a regression for the mean. The chapter is organized around these three points.

Section 2 states the estimand of each of Marigold's seven data blocks. Section 3 is the most detailed in the chapter: it first gives the properties of the quantile function and the generalized inverse, then treats QR as a nonsmooth M-estimator in the sense of Chapter 5, covering identification, consistency, and the density-weighted sandwich, and explains why Wooldridge (2010) Theorem 12.3, used in Section 3.4 of Chapter 5, fails here; it then takes up, in turn, the comparison of distributions in randomized experiments and under selection on observables, whether individual ranks are preserved, unconditional quantiles, endogenous treatment and instrumental variables, and merchant fixed effects in panels, and ends with the two-group, two-period design that matches by rank. Section 4 covers computation and inference, Section 5 reads three anchoring papers, Section 6 runs every tool on the Marigold data, and Section 7 discusses failure modes. The bootstrap of Chapter 8, the instrumental variables of Chapter 9, and the within transformation of Chapter 10 are all used here; Chapter 20, on more general difference-in-differences, returns to the quantile tools of this chapter.

## 1 A Tool with an Average Effect of Zero

::: {.case}
In all, 5,002 of Marigold's merchants took part in the AutoPrice experiment, split by complete randomization into two groups of exactly 2,501 each. The outcome is monthly revenue in the month after the tool is switched on, in thousands of USD. The control group averages 11.49 and the AutoPrice group 11.21; the control group's median is 10.04, that is, monthly revenue of around ten thousand dollars.

Because this is a simulation, we know the answer: the average effect of AutoPrice is exactly zero, but it stretches the revenue distribution out toward both ends. The 10th percentile of the AutoPrice group's revenue distribution falls 1.652 thousand dollars below the control group's, and the 90th percentile lies 2.772 thousand dollars above it.
:::

The analyst first does the most routine thing: an OLS regression of monthly revenue on an AutoPrice indicator. The coefficient is -0.284, with HC2 standard error 0.221 and a 95% interval [-0.717, 0.150] that covers zero; after adding the previous year's revenue and its interaction with treatment, the coefficient becomes 0.071, with interval [-0.239, 0.381]. The conclusion looks clear: AutoPrice has no effect on revenue.

Ask the question another way: line up each group's merchants from lowest to highest and compare them at the same position. At the 10th percentile, the AutoPrice group is lower by 1.747 thousand dollars, with interval [-1.946, -1.548]; at the 90th percentile it is higher by 2.750 thousand dollars, with interval [1.768, 3.733]. Both intervals are far from zero, and both cover the true values. The bottom loses, the top gains, and the mean lets the two ends cancel (figure below).

![Block Q1: the effect of AutoPrice on monthly revenue (thousands of USD). Dots and vertical bars are the differences in quantiles on the Q1 data at each level and their pointwise 95% intervals, at levels from 5% to 95% in steps of 5 percentage points; the solid line is the true difference in quantiles, which is the same for Q1 and Q1r; the dashed and dash-dotted lines are the true quantiles of the individual effects in Q1 and Q1r (read the horizontal axis as the level $u$); the horizontal dotted line and the light band are the OLS mean effect and its 95% interval. The vertical axis is truncated at -10 and 10 thousand dollars.](assets/fig/fig_13_qte.svg)

At this point it is tempting to say something stronger: "AutoPrice makes small merchants earn more than seventeen hundred dollars less a month, and large merchants more than twenty-seven hundred dollars more." That sentence is about what particular merchants gain and lose, whereas the estimates above only compared two distributions. To measure the gap between the two, Marigold has a parallel world, Q1r: the merchants, the covariates, the treatment assignment, and the control group's revenue are all the same as in Q1, and the population distribution of the AutoPrice group's revenue is also the same as in Q1; only the revenue with the tool is paired with the revenue without it in a different way. The two worlds have exactly the same differences in quantiles and, if we look only at revenue and treatment status $(Y_i, D_i)$, exactly the same distribution of the experimental data, yet their individual effects are entirely different: the merchants whose revenue falls make up 72.6% in Q1 and 57.6% in Q1r; the 10th percentile of the individual effects is -1.648 thousand dollars in Q1 and -7.979 in Q1r; the 90th percentiles are 2.772 and 8.884. Using only $(Y_i, D_i)$, no estimator can tell these two worlds apart.

This is where the chapter's question comes from: what exactly does the difference between the two groups' quantiles estimate, and under what conditions can it be read as "who gains and who loses"? And how do covariates, endogenous selection, panel structure, and before-and-after policy comparisons change the answer?

## 2 The Economic Model and the Estimand

### 2.1 Quantiles, Quantile Treatment Effects, and Individual Effects

Let $F_Y(y) = \Pr(Y \leq y)$ be the distribution function of monthly revenue and $f_Y$ its density. The $\tau$-th quantile is defined as

$$F_Y^{-1}(\tau) = \inf\{y : F_Y(y) \geq \tau\}, \qquad \tau \in (0,1),$$

that is, the smallest value such that "at least a share $\tau$ of merchants have revenue no greater than it." When the distribution function is strictly increasing, this is the ordinary inverse; when the distribution function has jumps or flat stretches, it is still defined and is called the generalized inverse. In this chapter, $\tau$ without a subscript denotes only a quantile level; the individual effect is still written $\tau_i$, always with a unit subscript; average effects are written ATE, ATT (the average treatment effect on the treated), and LATE, never a bare $\tau$. The $\tau$ in the regression equations of Chapter 8 is a treatment-effect parameter and means something different. The conditional quantile given covariates $X_i$ is written $\mathrm{Quant}_\tau(Y_i \mid X_i)$; it is the generalized inverse of the conditional distribution function $F_{Y \mid X}$, and the conditional density is written $f_{Y \mid X}$.

Let $Y_i(1)$ and $Y_i(0)$ be the monthly revenue of merchant $i$ with and without AutoPrice, its two potential outcomes, and let $D_i$ be the randomly assigned adoption indicator; we observe $Y_i = D_iY_i(1) + (1 - D_i)Y_i(0)$. The individual effect is $\tau_i = Y_i(1) - Y_i(0)$. The quantile treatment effect (QTE) is the difference between the distributions of the two potential outcomes at the same position:

$$\mathrm{QTE}(\tau) = F_{Y(1)}^{-1}(\tau) - F_{Y(0)}^{-1}(\tau).$$

It generally differs from the quantile of the distribution of individual effects, $F_{\tau_i}^{-1}(u)$. Here the individual effect $\tau_i$ and the quantile level appear together, so the level is written $u \in (0,1)$ for the moment; it is the same quantity as $\tau$ everywhere else in the chapter, and $\tau_i$ always carries a unit subscript. Chapter 20 uses the definition of the generalized inverse when it writes down Wald-CIC, and there too the quantile level is written $u$, which is this chapter's $\tau$. The condition that links the QTE to the distribution of individual effects is called rank invariance: every merchant sits at the same position in the distributions of both potential outcomes,

$$Y_i(1) = F_{Y(1)}^{-1}\big(F_{Y(0)}(Y_i(0))\big).$$

The DGP of Q1 satisfies it. Let $X_i$ be the standardized log revenue of the previous year and $\varepsilon_i$ an independent standard normal shock:

$$\log Y_i(0) = \log 10 + 0.4\,X_i + 0.3\,\varepsilon_i, \qquad Y_i(1) = c\,Y_i(0)^{1.4},$$

where the constant $c$ makes $\mathbb{E}[Y_i(1)] = \mathbb{E}[Y_i(0)]$, so the ATE is exactly zero. $Y_i(1)$ is a strictly increasing function of $Y_i(0)$, so ranks are invariant; the exponent exceeds 1, so merchants with high revenue gain and merchants with low revenue lose. Q1r keeps the distribution of $Y_i(1)$ but makes it only partly correlated with $Y_i(0)$: the logs of the two are paired through a normal copula (a function that joins two marginal distributions into a joint distribution) with correlation coefficient 0.5, and rank invariance no longer holds. The estimands of Q1 and Q1r are $\mathrm{QTE}(\tau)$ at five levels and the ATE. The chapter's five levels are the 10th, 25th, 50th, 75th, and 90th percentiles (that is, $\tau = 0.1$, $\tau = 0.9$, and so on); the figures use 19 levels from 5% to 95% in steps of 5 percentage points.

### 2.2 The Other Five Data Blocks

Block Q2 is about a shift in the distribution of a covariate. Of 6,002 merchants, 1,201 use Fulfilled by Marigold (FBM), under which the platform ships orders on the merchant's behalf; this is written $W_i = 1$, and in this chapter $W_i$ indicates only whether FBM is used and has nothing to do with the time-invariant regressor $W_i$ of Chapter 10. Log monthly revenue (in USD) is $\log Y_i = \log 8000 + 0.5\,W_i + 0.5\,X_i + 0.35\,\varepsilon_i$, where $X_i$ is the standardized log number of listings, independent of $W_i$. Given $(W_i, X_i)$, FBM shifts every conditional quantile of log revenue up by 0.5, a pure location shift. What the platform cares about, however, is something else: if the share of FBM rises a little from one fifth, how far does the $\tau$-th quantile of the platform-wide revenue distribution move? This quantity is called the unconditional quantile partial effect (UQPE); it varies with $\tau$, and at $\tau = 0.9$ it is 0.665, not 0.5.

Block Q3 is about endogenous treatment. AdBoost is paid promotion that merchants adopt voluntarily, and the merchants who adopt had better prospects to begin with. The platform sent an ad credit voucher at random to 3,001 of 6,002 merchants ($Z_i = 1$). Merchants come in three types: 30% adopt whether or not they have a voucher (always-takers), 40% adopt only with a voucher (compliers), and 30% do not adopt even with a voucher (never-takers). The potential outcomes are

$$\log Y_i(d) = \log 8000 + 0.10\,d + \sigma_d\,\Phi^{-1}(U_i) + 0.3\,X_i,$$

where $\sigma_0 = 0.5$, $\sigma_1 = 0.6$, $\Phi$ is the standard normal distribution function, $X_i$ is the standardized log follower count, and $U_i \sim \mathrm{U}(0,1)$ is the merchant's rank in the distribution of potential revenue given its follower count; it is positively correlated with the propensity to adopt, so adoption $D_i$ is endogenous. There are two different estimands here. One is the structural quantile effect for all merchants, $\beta_D(\tau) = 0.10 + 0.1\,\Phi^{-1}(\tau)$, which increases with $\tau$; the other is the compliers' local quantile treatment effect (LQTE), which also increases with $\tau$ but is closer to 0.10 at both ends. The compliers' ranks are concentrated in the middle and symmetric about the median, so the two are equal only at the median, where both are 0.10.

Block Q4 is about panels. Each of $N = 200$ merchants is observed for $\mathcal{T} = 5$ or 20 months, where $N$ denotes the number of merchants and $\mathcal{T}$ the number of panel periods. $X_{it} \in (0,1)$ is the share of listings in promotions, and larger merchants promote more. Log monthly revenue is $Y_{it} = \log 8000 + \alpha_i + 0.4\,X_{it} + 0.3\,(1 + X_{it})\,\varepsilon_{it}$, so that

$$\begin{aligned}
\mathrm{Quant}_\tau(Y_{it} \mid X_{it}, \alpha_i) &= \alpha_i + \log 8000 + 0.3\,\Phi^{-1}(\tau) + \beta(\tau)\,X_{it}, \\
\beta(\tau) &= 0.4 + 0.3\,\Phi^{-1}(\tau).
\end{aligned}$$

The merchant effect $\alpha_i$ only shifts the location; promotions act more strongly on the upper part of the distribution, and $\beta(0.9) = 0.7845$.

Blocks Q5 and Q5f are about before-and-after policy comparisons. One-Click Checkout launches in one region ($S_i = 1$) for the year-end peak season ($T_i = 1$) but not in another region ($S_i = 0$). Here $T_i \in \{0,1\}$ is the period indicator; it is unrelated to the number of periods observed for unit $i$ in an unbalanced panel, which Chapter 10 writes as $T_i$. The data are repeated cross-sections, with 2,500 merchants in each region-by-quarter cell. Untreated monthly revenue (in thousands of USD) is

$$\begin{aligned}
Y_i(0) &= h(U_i, T_i), \\
h(U_i, t) &= 10\exp\big(0.05\,t + 0.5\,(1 + 0.4\,t)\,U_i\big),
\end{aligned}$$

$U_i$ follows $\mathcal{N}(0, 1)$ in the control region and $\mathcal{N}(0.6, 0.8^2)$ in the launch region, the same in both seasons: merchants in the launch region are larger. Here $U_i$ is not normalized to a uniform distribution, and its distribution differs between the two regions; the CIC of Section 3.9 requires only that it be one-dimensional. Peak-season growth, $h(U_i, 1)/h(U_i, 0) = \exp(0.05 + 0.2\,U_i)$, increases with rank, so large merchants grow more and the time effect is nonlinear in revenue. One-Click Checkout raises revenue by 12%. In Q5, $D_i = S_iT_i$, and the estimand is the QTT (quantile treatment effect on the treated) of the launch region in the peak season,

$$\mathrm{QTT}(\tau) = F_{Y(1) \mid D = 1}^{-1}(\tau) - F_{Y(0) \mid D = 1}^{-1}(\tau),$$

which increases with $\tau$ and is 1.920 thousand dollars at the median. Q5f is the voluntary-adoption version, with 10,000 merchants per cell: in the launch region, 10% adopted during the beta test, 60% start adopting in the peak season, and 30% never adopt; in the control region, 15% have been using an older version of One-Click Checkout all along, a share that is the same in both seasons. The merchants who start adopting in the peak season are called switchers. The estimand is the switchers' $\mathrm{LQTE}_{sw}(\tau)$, 2.047 thousand dollars at the median; the switchers' average effect (LATE) is 2.350.

The seven estimands have one thing in common: each is a difference between distributions at some position, not anyone's effect. What differs is "which distribution": conditional or unconditional; all merchants, compliers, switchers, or the treated group. Section 3 explains, one by one, under what assumptions each can be identified from the data.

## 3 Identification and Large-Sample Theory

This section follows the plan of Chapter 5: Marigold's numbers or a small example first show what is to be proved, then come numbered assumptions and theorems, followed by proof sketches of two to six steps; each step names the assumptions it uses, and regularity conditions are only cited to their source. After each theorem, an "In words" paragraph answers four questions: under which assumptions it estimates which quantity, what it estimates when the specification is wrong, who sets the weights, and what Marigold's numbers say. Assumptions and theorems are numbered within this chapter; the parentheses give the numbers in Wooldridge (2010), Hansen (2014), or the original papers, and when the numbers of a working-paper version are cited, the version is named. Sections 3.1 to 3.4 are the general theory of QR, and Sections 3.5 to 3.9 cover five causal designs.

### 3.1 The Quantile Function, the Generalized Inverse, and the Check Function

Start with a distribution that takes only three values: $Y \in \{1, 2, 3\}$, with $F_Y$ equal to 0.2, 0.5, and 1 at the three points. By definition, $F_Y^{-1}(0.5) = 2$ and $F_Y^{-1}(0.6) = 3$: as soon as the level moves a little past 0.5, the quantile jumps to the next value. Conversely, $F_Y(F_Y^{-1}(0.3)) = F_Y(2) = 0.5$, which is larger than 0.3: at the jump, the distribution function overshoots its level. A sample quantile is the generalized inverse of the empirical distribution function, and the empirical distribution function is a staircase like this everywhere, so the inequalities below are the norm in samples.

Quantiles have a second identity: they solve a minimization problem. Define the check function

$$\rho_\tau(u) = \big(\tau - \mathbf{1}\{u < 0\}\big)\,u,$$

which for $u > 0$ has slope $\tau$ and for $u < 0$ has slope $\tau - 1$, so it looks like a tilted check mark. Take 500 random draws from an exponential distribution and find the value that minimizes $\sum_i \rho_{0.3}(Y_i - \theta)$: this $\theta$ lies between the 150th and 151st order statistics. The reason is simple: move $\theta$ up a little, and the loss of each observation below it rises by $1 - \tau$ times the move, while the loss of each observation above it falls by $\tau$ times the move; when exactly three tenths of the observations lie below, the increases and decreases cancel.

::: {.theorem}
**Theorem 1 (the generalized inverse and the check function; Wacker 2023 Lemma 1 (g)(3), (h)(3), (i), and (l), pp. 2-3; Athey and Imbens 2006 working-paper version, after equation (3.8), p. 9; Wooldridge 2010 equation (12.104), p. 450; Hansen 2014 equation (9.14), p. 225)** Let $F_Y$ be a distribution function and $\tau \in (0,1)$.

\(i) For all $y$, $F_Y^{-1}(\tau) \leq y$ if and only if $\tau \leq F_Y(y)$; hence $F_Y(F_Y^{-1}(\tau)) \geq \tau$ and $F_Y^{-1}(F_Y(y)) \leq y$; when $F_Y$ is continuous at $F_Y^{-1}(\tau)$, the first inequality holds with equality.

\(ii) If $\mathbb{E}\lvert Y \rvert < \infty$, the set of minimizers of $\mathbb{E}[\rho_\tau(Y - \theta)]$ is $\{\theta : \Pr(Y < \theta) \leq \tau \leq F_Y(\theta)\}$, and $F_Y^{-1}(\tau)$ is its smallest element; when there is no interval on which $F_Y$ is identically $\tau$, it is the unique minimizer.

\(iii) If $g$ is strictly increasing and continuous, then $F_{g(Y)}^{-1}(\tau) = g\big(F_Y^{-1}(\tau)\big)$; the same holds for conditional distributions, $\mathrm{Quant}_\tau(g(Y_i) \mid X_i) = g\big(\mathrm{Quant}_\tau(Y_i \mid X_i)\big)$.
:::

The proof sketch has three steps. Step one ($F_Y$ nondecreasing and right-continuous): the set $\{y : F_Y(y) \geq \tau\}$ is a closed right half-line whose left endpoint is $F_Y^{-1}(\tau)$, which is the equivalence in (i); taking $y = F_Y^{-1}(\tau)$ gives the first inequality, and taking $\tau = F_Y(y)$ gives the second; when $F_Y$ is continuous at the left endpoint, its values to the left of that point are all below $\tau$, and taking limits gives the reverse inequality. Step two ($\mathbb{E}\lvert Y \rvert < \infty$): $\mathbb{E}[\rho_\tau(Y - \theta)]$ is a convex function of $\theta$, with right derivative $F_Y(\theta) - \tau$ and left derivative $\Pr(Y < \theta) - \tau$; a convex function attains its minimum where the left derivative is nonpositive and the right derivative nonnegative, and step one then picks out the smallest such point. Step three ($g$ strictly increasing and continuous): $\{g(Y) \leq g(y)\} = \{Y \leq y\}$, so $F_{g(Y)}(g(y)) = F_Y(y)$, and among the values with $F_{g(Y)}(z) \geq \tau$, the smallest $z$ is the value of $g$ at $F_Y^{-1}(\tau)$; with a conditional distribution in place of the distribution, the argument goes through word for word. Wooldridge (2010) p. 453 states this property for the median: the median passes through monotone functions and the expectation does not; step three is the case of a general $\tau$.

In words: a quantile has two equivalent definitions, as the generalized inverse of the distribution function or as the minimizer of the expected check function. The latter puts quantiles into the framework of Chapter 5: the estimand is the minimizer of some population loss. There is no model here that has to be correctly specified; the population distribution itself determines the answer. When the distribution has a flat stretch at level $\tau$, the minimizer is an interval, and the generalized inverse takes its left endpoint by convention, in samples as well. The weights are set by the two slopes of the check function: the larger $\tau$, the cheaper the observations below and the costlier those above, so the minimizer is pushed toward the upper end of the distribution until exactly a share $\tau$ of the probability mass lies to its left. Property (iii) is a real advantage of quantiles over means. At the 90th percentile of Q1, the true quantile of the AutoPrice group is 21.751 thousand dollars and that of the control group is 18.980, and the difference of their logs, 0.136, is the QTE of log monthly revenue at this level; means have no such property, and the average effect on log revenue cannot be computed from the average effect on revenue.

### 3.2 Conditional Quantile Regression: Identification and Consistency

On block Q2, a QR of log monthly revenue on the FBM indicator and log listings gives an FBM coefficient of 0.473 at $\tau = 0.1$, 0.492 at the median, and 0.523 at $\tau = 0.9$, all near the true value 0.5; the 95% intervals at all 19 levels cover 0.5. QR generalizes the univariate problem of Theorem 1 (ii) to regression:

$$\hat\beta(\tau) = \arg\min_{b \in \mathcal{B}} \frac{1}{n}\sum_{i=1}^n \rho_\tau(Y_i - X_i'b),$$

that is, the sample problem of Wooldridge (2010) equation (12.105) (Koenker and Bassett 1978). It is an M-estimator in the sense of Chapter 5, whose loss for a single observation is $q(w, b) = \rho_\tau(Y - X'b)$.

::: {.assumption}
**QR.1 (linear conditional quantile; Wooldridge 2010 equation (12.103), p. 450)** For the given $\tau$, $\mathrm{Quant}_\tau(Y_i \mid X_i) = X_i'\beta(\tau)$, and the first element of $X_i$ is 1.

**QR.2 (identification; p. 451)** $\beta(\tau)$ is the unique minimizer of $\mathbb{E}[\rho_\tau(Y_i - X_i'b) - \rho_\tau(Y_i)]$ on $\mathcal{B}$. A sufficient condition: given $X_i$, the conditional density of $Y_i$ is continuous and positive at $X_i'\beta(\tau)$, and $\mathbb{E}[X_iX_i']$ is positive definite.

**QR.3 (regularity conditions; the conditions of Wooldridge 2010 Theorem 12.2, p. 404)** $\{(Y_i, X_i)\}$ is a random sample, $\mathcal{B}$ is a compact set containing $\beta(\tau)$, and $\mathbb{E}\lVert X_i \rVert < \infty$.
:::

QR.1 is the only economic assumption: the conditional quantile is linear in the parameters. It holds separately for each $\tau$ and does not require the slopes of different quantiles to be the same; the error can be heteroskedastic and arbitrarily distributed. Subtracting $\rho_\tau(Y_i)$ in QR.2 does not change the minimizer, but it means the expectation requires only a first moment of $X_i$; $Y_i$ need not have a mean. The sufficient condition means that there are data near the true quantile (a positive density) and no perfect collinearity among the regressors. QR.3 corresponds to M.2 of Chapter 5: the check function is continuous everywhere, and Wooldridge (2010) p. 449 notes that the consistency of LAD comes from Theorem 12.2, and the same argument covers a general $\tau$.

::: {.theorem}
**Theorem 2 (consistency of QR; Wooldridge 2010 §12.10.1, p. 451, via Theorem 12.2, p. 404; Koenker and Bassett 1978)** Under QR.1 through QR.3, $\hat\beta(\tau) \xrightarrow{p} \beta(\tau)$.
:::

The proof sketch has three steps. Step one (QR.1, Theorem 1 (ii)): given $X_i$, $\mathbb{E}[\rho_\tau(Y_i - \theta) \mid X_i]$ is minimized at $\theta = \mathrm{Quant}_\tau(Y_i \mid X_i) = X_i'\beta(\tau)$, so $\beta(\tau)$ minimizes the conditional expected loss at every $X_i$, and after taking expectations it is a minimizer of the population problem (equation (12.104)). Step two (QR.2): the minimizer is unique, which is M.1 of Chapter 5. Step three (QR.3): $\lvert \rho_\tau(a) - \rho_\tau(a') \rvert \leq \lvert a - a' \rvert$, so $\rho_\tau(Y - X'b) - \rho_\tau(Y)$ is continuous in $b$ and dominated by $\lVert X \rVert \sup_{b \in \mathcal{B}}\lVert b \rVert$, which has finite expectation; all the conditions of Theorem 2 of Chapter 5 are met.

In words: under QR.1 through QR.3, QR estimates the parameters of the conditional quantile: when $X$ changes by one unit, how far the $\tau$-th position of the conditional distribution moves. It says "the 90th percentile of the revenue distribution moves up," not "the merchant who was at the 90th percentile earns more"; the latter needs the rank condition of Section 3.5. When the conditional quantile is not linear, QR converges to the pseudo-true value of Theorem 4 in Section 3.4. The weights are set by the check function: observations above the fitted line are penalized at $\tau$ and those below at $1 - \tau$, and the sample first-order condition puts roughly $n\tau$ observations below the fitted hyperplane; a distant observation counts only through its direction, not its distance, so QR is insensitive to extreme values of $Y$. This robustness concerns only the conditional quantile: Wooldridge (2010) p. 451 warns that when the parameters of the conditional mean are wanted, LAD is not robust, and it estimates the same thing as OLS only under extra assumptions such as a symmetric conditional distribution. On Q2 the conditional QR coefficients of FBM are near 0.5 at all five levels, and the reading is: of two merchants with the same number of listings, the one that uses FBM has its whole conditional distribution shifted to the right by 0.5 log points.

### 3.3 Asymptotic Normality: The Density-Weighted Sandwich

Q2 again. With the same data and the same regressor, the standard error of the FBM coefficient is 0.0232 at $\tau = 0.1$, 0.0144 at the median, and 0.0214 at $\tau = 0.9$: the two ends of the distribution are harder to estimate than the middle. The reason is the density of the data near the quantile: near the median the merchants are dense, and a slight shift of the fitted line changes the sign of many residuals; in the tails they are sparse, and the fitted line must move a long way before it meets the next observation.

Theorem 3 of Section 3.4 of Chapter 5 cannot be used here. That theorem requires the loss of a single observation to be twice differentiable and the expected Hessian to be positive definite (conditions M.3(b) and (d)), whereas the check function is not differentiable at zero and has second derivative zero everywhere else, so the Hessian is the zero matrix. Wooldridge (2010) points out explicitly, on p. 407 and again on p. 454, that LAD and QR in general fall outside Theorem 12.3. The remedy is to change the order of differentiation: first take expectations over the data, smoothing out the nonsmooth score, and then differentiate.

The QR score is $s_i(b) = X_i\big(\mathbf{1}\{Y_i < X_i'b\} - \tau\big)$; it is undefined at $Y_i = X_i'b$, but at the true value this happens with probability zero. Taking expectations given $X_i$,

$$\mathbb{E}[s_i(b) \mid X_i] = X_i\big(F_{Y \mid X}(X_i'b \mid X_i) - \tau\big)$$

(equation (12.111)), which is a smooth function of $b$ whose Jacobian at $\beta(\tau)$ is

$$A_0(\tau) = \mathbb{E}\big[f_{Y \mid X}(X_i'\beta(\tau) \mid X_i)\,X_iX_i'\big].$$

Here $A_0(\tau)$ is a density-weighted matrix, not the expected Hessian of $q$: the second derivative of the check function is zero wherever it exists, and Wooldridge (2010) Theorem 12.3, used in Chapter 5, does not apply. The filling of the sandwich is written $B_0(\tau) = \mathbb{E}[s_i(\beta(\tau))\,s_i(\beta(\tau))']$.

::: {.assumption}
**QR.4 (Wooldridge 2010 §12.10.2, pp. 454-455; Angrist, Chernozhukov and Fernández-Val 2006 working-paper version, conditions A.4 and A.5, p. 12)** In addition to QR.1 through QR.3: (a) given $X_i$, the conditional density of $Y_i$, $f_{Y \mid X}$, exists, is bounded, and is uniformly continuous in $y$, uniformly in $x$ over the support; (b) $A_0(\tau)$ is positive definite; (c) for some $\delta > 0$, $\mathbb{E}\lVert X_i \rVert^{2+\delta} < \infty$.
:::

QR.4(a) validates the smoothing obtained by taking expectations first (step two of the proof) and makes the residuals that are exactly zero near the solution negligible (step one); (b) says that $A_0(\tau)$ is invertible, replacing the positive definite Hessian of M.3(d) in Chapter 5; (c) is what the central limit theorem and the empirical-process arguments require of the regressors.

::: {.theorem}
**Theorem 3 (asymptotic normality of QR; Wooldridge 2010 equations (12.111) to (12.115), pp. 454-455; Hansen 2014 Theorem 9.6.1, p. 226; under misspecification, see Angrist, Chernozhukov and Fernández-Val 2006 working-paper version, Theorem 5 and Corollary 1, pp. 13-14)** Under QR.1 through QR.4,

$$\begin{aligned}
\sqrt{n}\big(\hat\beta(\tau) - \beta(\tau)\big) &= -A_0(\tau)^{-1}\,\frac{1}{\sqrt n}\sum_{i=1}^n s_i(\beta(\tau)) + o_p(1) \\
&\xrightarrow{d} \mathcal{N}\big(0,\ A_0(\tau)^{-1}B_0(\tau)A_0(\tau)^{-1}\big),
\end{aligned}$$

and $B_0(\tau) = \tau(1 - \tau)\,\mathbb{E}[X_iX_i']$ (equation (12.115)). If $f_{Y \mid X}(X_i'\beta(\tau) \mid X_i)$ does not vary with $X_i$, the asymptotic variance simplifies to $\tau(1 - \tau)\big(\mathbb{E}[X_iX_i']\big)^{-1}$ divided by the square of this density (equation (12.119)).
:::

The proof sketch has four steps. Step one (QR.4(a)(c)): at the solution of the sample problem, with probability 1 exactly $k$ residuals are zero ($k$ is the dimension of $\beta$; QR.4(a) makes the probability that more than $k$ points lie on a common hyperplane zero), and they make the sample score deviate from zero by at most $k\max_i\lVert X_i \rVert$, which tends to zero after division by $\sqrt n$, so the first-order condition holds approximately (Wooldridge 2010 p. 454). Step two (QR.1, QR.4(a)): first take conditional expectations to get equation (12.111), which is smooth in $b$ and zero at $\beta(\tau)$ (QR.1 makes $F_{Y \mid X}(X_i'\beta(\tau) \mid X_i) = \tau$), with Jacobian $A_0(\tau)$ (equation (12.112)). Step three (QR.4(b)(c)): the difference between the sample score and its expectation is asymptotically the same at $\hat\beta(\tau)$ and at $\beta(\tau)$; this is a stochastic equicontinuity argument, for which Wooldridge draws on Huber (1967) and Newey and McFadden (1994) §7; combined with step two, and using the invertibility of $A_0(\tau)$, it gives the linear representation (equation (12.113)). Step four (QR.1, QR.4(c)): the $s_i(\beta(\tau))$ are independent and identically distributed, and QR.1 gives $\mathbb{E}[\mathbf{1}\{Y_i < X_i'\beta(\tau)\} \mid X_i] = \tau$, so the mean is zero and the variance is $\tau(1 - \tau)\,\mathbb{E}[X_iX_i']$ (equation (12.115)); the central limit theorem completes the proof.

In words: under QR.1 through QR.4, $\hat\beta(\tau)$ fluctuates normally around the parameters of the conditional quantile at the $\sqrt n$ rate; the shape is still the sandwich of Chapter 5, except that the bread has become the density-weighted $A_0(\tau)$. Under misspecification, the step $B_0(\tau) = \tau(1 - \tau)\mathbb{E}[X_iX_i']$ no longer holds, because it uses the fact that QR.1 makes the conditional mean of the indicator exactly $\tau$. Theorem 5 in the working-paper version of Angrist, Chernozhukov and Fernández-Val (2006), NBER WP 10428, proves that under misspecification the estimator converges to the pseudo-true value $\beta_*(\tau)$ and is still asymptotically normal, but the filling must be replaced by $\mathbb{E}[(\tau - \mathbf{1}\{Y_i < X_i'\beta_*(\tau)\})^2X_iX_i']$; the common standard errors computed with $\tau(1-\tau)$ are therefore no longer robust. The weights are set by the conditional density in $A_0(\tau)$: where the data near the true quantile are dense, the curvature is large and the estimate precise; where they are sparse, the curvature is small and the estimate poor. The reciprocal of the density is therefore called the sparsity. In Q2 the error is independent of the regressors, so the simplified formula holds, and the standard errors vary only with $\tau(1 - \tau)$ and with the density at the quantile: 0.0232 at $\tau = 0.1$ and 0.0144 at the median, so the ends are hard to estimate and the middle is easy.

### 3.4 Misspecification: A Weighted Best Linear Approximation

On block Q4, treat the panel as a cross-section and run a pooled QR of log revenue on the promotion share, ignoring the merchant effects. In 1,000 Monte Carlo replications, the average pooled QR estimate at $\tau = 0.9$ is 1.823 for $\mathcal{T} = 5$ and 1.810 for 20, against a true value of 0.7845. Adding periods does not help: the probability limit of pooled QR ($N \to \infty$) does not change with the number of periods; it consistently estimates something else. Given only the promotion share, the conditional quantile of revenue also absorbs the differences in merchant size, and it need not be linear. When the linear model is misspecified, what does QR estimate?

::: {.theorem}
**Theorem 4 (QR under misspecification; Angrist, Chernozhukov and Fernández-Val 2006 working-paper version, Theorem 1, p. 4; cited in Wooldridge 2010 p. 453)** Suppose the conditional density given $X$ exists, and that the conditional quantile $\mathrm{Quant}_\tau(Y \mid X)$ and the population QR coefficient $\beta_*(\tau) = \arg\min_b \mathbb{E}[\rho_\tau(Y - X'b)]$ are both unique. Then

$$\beta_*(\tau) = \arg\min_b\ \mathbb{E}\Big[\bar f_\tau(X, b)\,\big(X'b - \mathrm{Quant}_\tau(Y \mid X)\big)^2\Big],$$

where

$$\begin{aligned}
\bar f_\tau(X, b) &= \int_0^1 2(1 - v)\,f_{Y \mid X}(m_v \mid X)\,dv, \\
m_v &= v\,X'b + (1 - v)\,\mathrm{Quant}_\tau(Y \mid X).
\end{aligned}$$

$m_v$ runs along the segment between the fitted value and the true quantile as $v$ varies, and $\bar f_\tau$ is a weighted average of the conditional density over this segment. When QR.1 holds, $\beta_*(\tau) = \beta(\tau)$.
:::

The proof sketch has two steps. Step one (the conditional density exists): write $\bar\rho_\tau(\theta \mid X) = \mathbb{E}[\rho_\tau(Y - \theta) \mid X]$; its first derivative in $\theta$ is $F_{Y \mid X}(\theta \mid X) - \tau$, which is zero at the true quantile (Theorem 1 (ii)), and its second derivative is the conditional density; a Taylor expansion with integral remainder gives $\bar\rho_\tau(X'b \mid X) - \bar\rho_\tau(\mathrm{Quant}_\tau(Y \mid X) \mid X) = \tfrac12\,\bar f_\tau(X, b)\,\big(X'b - \mathrm{Quant}_\tau(Y \mid X)\big)^2$. Step two (the two uniqueness conditions): subtracting the term that does not involve $b$, $\bar\rho_\tau(\mathrm{Quant}_\tau(Y \mid X) \mid X)$, leaves the minimizer unchanged, and taking expectations over $X$ gives the result. The original paper's weight is half of $\bar f_\tau$, which gives the same minimizer.

In words: if the conditional quantile is linear, QR estimates $\beta(\tau)$ itself; if it is not, QR estimates $\beta_*(\tau)$, a weighted best linear approximation to the conditional quantile function, for the same reason that NLS approximates the conditional mean in Theorem 5 of Chapter 5. The weights come in two layers: the distribution of $X$ decides how large a share each region gets, and $\bar f_\tau$ is what Angrist, Chernozhukov and Fernández-Val call the importance weight; the denser the conditional distribution near the true quantile, the more the approximation error there matters. In their census example the importance weights are fairly flat, and QR nearly coincides with the minimum-distance fit weighted only by the distribution of $X$ (Chapter 7). The pooled QR on Q4 is a case in point: it faithfully approximates the conditional quantile of revenue given only the promotion share, but that conditional quantile mixes in the differences in merchant size, so 1.823 is far from 0.7845, and no number of periods can correct it. Conversely, when the only regressor is a treatment indicator taking the values 0 and 1, the model is saturated and cannot be misspecified: on Q1, a QR on $D_i$ alone gives a coefficient that is exactly the difference between the two groups' sample quantiles.

### 3.5 Randomized Experiments, Selection on Observables, and Rank Invariance

Back to Section 1: Q1 and Q1r have exactly the same QTE, yet their distributions of individual effects differ. Making sense of this takes three assumptions.

::: {.assumption}
**QTE.1 (random assignment)** $(Y_i(0), Y_i(1)) \perp D_i$ and $0 < \Pr(D_i = 1) < 1$.

**QTE.2 (selection on observables; Firpo 2007)** $(Y_i(0), Y_i(1)) \perp D_i \mid X_i$, and the propensity score $e(x) = \Pr(D_i = 1 \mid X_i = x)$ satisfies $0 < e(x) < 1$.

**QTE.3 (rank invariance; Heckman, Smith and Clements 1997)** There exists $U_i \sim \mathrm{U}(0,1)$ such that $Y_i(d) = F_{Y(d)}^{-1}(U_i)$, $d = 0, 1$.
:::

QTE.1 is the design of Q1. QTE.2 replaces random assignment with "as good as randomly assigned, given the covariates"; the first half is called unconfoundedness and the second half overlap, and Chapter 15 discusses both in detail. QTE.3 concerns how the two potential outcomes are paired, and since each merchant reveals only one of them in the data, it cannot be tested from $(Y_i, D_i)$ alone; it allows effects to differ across merchants but requires that treatment not reorder the merchants.

::: {.theorem}
**Theorem 5 (QTE and individual effects; Firpo 2007; Heckman, Smith and Clements 1997; Bitler, Gelbach and Hoynes 2006 working-paper version, §4.1, pp. 16-17)**

\(i) Under QTE.1, $F_{Y \mid D = d} = F_{Y(d)}$, so $\mathrm{QTE}(\tau) = F_{Y \mid D = 1}^{-1}(\tau) - F_{Y \mid D = 0}^{-1}(\tau)$. Under QTE.2, $F_{Y(1)}(y) = \mathbb{E}\big[D_i\mathbf{1}\{Y_i \leq y\}/e(X_i)\big]$ and $F_{Y(0)}(y) = \mathbb{E}\big[(1 - D_i)\mathbf{1}\{Y_i \leq y\}/(1 - e(X_i))\big]$, and the QTE is identified by the difference of their generalized inverses.

\(ii) With QTE.3 added, $\tau_i = \mathrm{QTE}(U_i)$: the effect of the merchant with rank $\tau$ is exactly $\mathrm{QTE}(\tau)$. If $\mathrm{QTE}(\cdot)$ is continuous and strictly increasing, then $F_{\tau_i}^{-1}(u) = \mathrm{QTE}(u)$.

\(iii) Without QTE.3, the two marginal distributions do not determine the distribution of $\tau_i$. But if at some $\tau$ we have $\mathrm{QTE}(\tau) < 0$, then $\Pr(\tau_i < 0) > 0$; if at some $\tau$ we have $\mathrm{QTE}(\tau) > 0$, then $\Pr(\tau_i > 0) > 0$.
:::

The proof sketch has three steps. Step one (QTE.1 or QTE.2): under random assignment, $\Pr(Y_i \leq y \mid D_i = d) = \Pr(Y_i(d) \leq y \mid D_i = d) = F_{Y(d)}(y)$, and then take generalized inverses; under QTE.2, by iterated expectations, $\mathbb{E}[D_i\mathbf{1}\{Y_i \leq y\}/e(X_i)] = \mathbb{E}\big[e(X_i)\Pr(Y_i(1) \leq y \mid X_i)/e(X_i)\big] = F_{Y(1)}(y)$, where $e(X_i) > 0$ makes the division meaningful, and the untreated side is analogous. Step two (QTE.3): $\tau_i = F_{Y(1)}^{-1}(U_i) - F_{Y(0)}^{-1}(U_i) = \mathrm{QTE}(U_i)$; when $\mathrm{QTE}(\cdot)$ is continuous and strictly increasing, Theorem 1 (iii) says that the quantile of $\mathrm{QTE}(U_i)$ at level $u$ is $\mathrm{QTE}(u)$. Step three (by contradiction): if $\tau_i \geq 0$ almost surely, then $Y_i(1) \geq Y_i(0)$ and $F_{Y(1)} \leq F_{Y(0)}$, so that for all $\tau$, $\mathrm{QTE}(\tau) \geq 0$; the other direction is analogous. For the first half of (iii), Q1 and Q1r give a counterexample: the $(Y_i, D_i)$ of the two worlds have the same distribution, while the distributions of $\tau_i$ differ.

In words: in a randomized experiment, the difference between the two groups' quantiles estimates the difference between the two marginal distributions at the same position, and this step uses only random assignment; under selection on observables, it is the difference between two distributions reweighted by the inverse of the propensity score, and if the propensity score model is misspecified, the weights are wrong and what is estimated is the difference between two improperly weighted distributions. Who sets the weights: in a randomized experiment each quantile is determined by its own group's sample alone; under selection on observables, by $1/e(X_i)$. Reading the QTE as "who gains and who loses" also requires rank invariance, an assumption about counterfactual pairing that the data cannot test: Q1 and Q1r have the same QTE, yet the merchants whose revenue falls make up 72.6% and 57.6% respectively. Even when rank invariance holds, the QTE need not be the quantile of the individual effects: the individual effect in Q1 is a U-shaped function of revenue without the tool, with a minimum of -1.663 thousand dollars where revenue without the tool is 5.821 thousand dollars, so at the 25th percentile the QTE is -1.606 while the 25th percentile of the individual effects is -1.562; only at the 75th and 90th percentiles do the two coincide. Nor are we entirely in the dark without rank invariance: the significantly negative QTE at the 10th percentile of Q1 is strong evidence that some merchants are indeed hurt, though it cannot say how many.

Covariates add a little information. The conditional distributions of Q1 and Q1r given the previous year's revenue are not the same, so neither are the distributions of $(Y_i, D_i, X_i)$: in a QR of log monthly revenue on $D_i$, $X_i$, and $D_iX_i$, the coefficient on $D_i$ at the 90th percentile is 0.041 in Q1 and 0.308 in Q1r. Some information about the pairing is hidden in the covariates, but by themselves they say nothing about the distribution of individual effects unless further assumptions are added.

### 3.6 Unconditional Quantiles and Counterfactual Distributions

On block Q2, a recentered influence function (RIF) regression gives 0.723 at the 90th percentile, against a true UQPE of 0.665; conditional QR gives 0.523, which aims at 0.5. Neither number is wrong; they answer two different questions. Switch one merchant from not using FBM to using it, and its log revenue moves up by 0.5; but how far the 90th percentile of the platform-wide revenue distribution moves depends on how many merchants of each group lie near that quantile, and also on how dense the merchants are there.

The RIF is the influence function of the quantile plus the quantile itself:

$$\mathrm{RIF}(Y;\tau) = F_Y^{-1}(\tau) + \frac{\tau - \mathbf{1}\{Y \leq F_Y^{-1}(\tau)\}}{f_Y(F_Y^{-1}(\tau))},$$

and its expectation is exactly $F_Y^{-1}(\tau)$. When the share of FBM is $p$, the platform-wide distribution is $F_Y(y; p) = p\,F_{Y \mid W = 1}(y) + (1 - p)\,F_{Y \mid W = 0}(y)$, and the UQPE is defined as the derivative of the quantile when the share rises slightly, $\mathrm{UQPE}(\tau) = \partial F_Y^{-1}(\tau; p)/\partial p$.

::: {.assumption}
**RIF.1 (a change in the share does not change the conditional distributions; the counterfactual of Firpo, Fortin and Lemieux 2009 working-paper version, Corollary 3, p. 13)** When the share of FBM changes from $p$ to $p'$, given $W_i$ the joint distribution of $(Y_i, X_i)$ does not change, so $F_Y(y; p') = p'F_{Y \mid W = 1}(y) + (1 - p')F_{Y \mid W = 0}(y)$.

**RIF.2 (density)** $f_Y$ is continuous and positive in a neighborhood of $F_Y^{-1}(\tau)$.
:::

RIF.1 is the core of counterfactuals of this kind: the merchants who newly join FBM are no different from the existing FBM merchants, the policy does not change the revenue distribution of the FBM merchants themselves, and there is no general equilibrium feedback. It is an assumption about the policy, not about the data.

::: {.theorem}
**Theorem 6 (the unconditional quantile partial effect; Firpo, Fortin and Lemieux 2009 working-paper version, Corollary 3 and equation (15), pp. 13-14)** Under RIF.1 and RIF.2,

$$\begin{aligned}
\mathrm{UQPE}(\tau) &= \mathbb{E}[\mathrm{RIF}(Y_i;\tau) \mid W_i = 1] - \mathbb{E}[\mathrm{RIF}(Y_i;\tau) \mid W_i = 0] \\
&= \frac{F_{Y \mid W = 0}(F_Y^{-1}(\tau)) - F_{Y \mid W = 1}(F_Y^{-1}(\tau))}{f_Y(F_Y^{-1}(\tau))}.
\end{aligned}$$

If $W_i$ is independent of $X_i$, then in the population OLS of $\mathrm{RIF}(Y_i;\tau)$ on $(1, W_i, X_i)$, the coefficient on $W_i$ is $\mathrm{UQPE}(\tau)$.
:::

The proof sketch has three steps. Step one (RIF.1): $F_Y(y; p)$ is linear in $p$, with derivative $F_{Y \mid W = 1}(y) - F_{Y \mid W = 0}(y)$. Step two (RIF.2, Theorem 1 (i)): the distribution is continuous near the quantile, so $F_Y(F_Y^{-1}(\tau; p); p) = \tau$ holds with equality; differentiating with respect to $p$ by the implicit function theorem gives $f_Y \cdot \partial F_Y^{-1}/\partial p + F_{Y \mid W = 1} - F_{Y \mid W = 0} = 0$, and solving gives the second equality. Step three (the definition of the RIF): $\mathbb{E}[\mathrm{RIF}(Y_i;\tau) \mid W_i = w]$ equals $F_Y^{-1}(\tau)$ plus $\tau - F_{Y \mid W = w}(F_Y^{-1}(\tau))$ divided by $f_Y(F_Y^{-1}(\tau))$, and differencing the two groups gives the first equality; when $W_i$ is independent of $X_i$, by the Frisch-Waugh-Lovell theorem the OLS coefficient on $W_i$ equals $\mathrm{Cov}(W_i, \mathrm{RIF})/\mathrm{Var}(W_i)$, which is the difference between the two groups' RIF means.

In words: under RIF.1 and RIF.2, RIF regression estimates the effect of a small change in the share on the unconditional quantile; it asks "what happens to the revenue distribution of the whole platform," whereas conditional QR asks "how do two merchants with the same number of listings compare." The two agree for the mean (the RIF of the mean is $Y$ itself) and part ways at quantiles. The OLS of the RIF on $(1, W_i, X_i)$ is a linear probability approximation to $\Pr(Y_i > F_Y^{-1}(\tau) \mid W_i, X_i)$, and when $W_i$ and $X_i$ are correlated it is in general only an approximation; Firpo, Fortin and Lemieux also report a logit version for comparison. The weights are set by the unconditional density in the denominator: the sparser the merchants near the quantile, the farther a given number of merchants switching groups pushes the quantile. In Q2, $W_i$ and $X_i$ are independent; RIF regression gives 0.723 at the 90th percentile, with interval [0.619, 0.826] covering 0.665; conditional QR gives 0.523, with interval [0.481, 0.565], which covers its own target 0.5 but not 0.665. It is not a wrong estimate; it answers a different question. At the median the UQPE is 0.496, almost the same as the conditional effect; away from the median, the two move apart in opposite directions: at the low end the UQPE is below 0.5 (0.346 at $\tau = 0.1$), and at the high end it is above 0.5. The reason lies in the numerator and the denominator of the second equality of Theorem 6. The numerator is the difference between the two groups' distribution functions at the quantile, that is, the probability that a merchant crosses this quantile when it switches groups; it is largest near the median and shrinks toward both ends, but shrinks more at the low end: the low end of the distribution consists almost entirely of merchants who do not use FBM, so the probability of crossing a low quantile is small; at the high end the share of FBM merchants is high, and the probability of crossing is much larger. The denominator is the density there, which differs little between the two ends. So the effect is small at the low end and large at the high end (figure below).

![Block Q2: two kinds of FBM effect on log monthly revenue. Circles are the FBM coefficients in the RIF regression and triangles the FBM coefficients in conditional QR; vertical bars are pointwise 95% intervals (from a bootstrap of the whole procedure for the RIF, from nid standard errors for conditional QR), at levels from 5% to 95%; the solid line is the true UQPE, and the dashed line is the true conditional effect, 0.5 at every level.](assets/fig/fig_13_uqpe.svg)

Another route is to estimate the whole conditional distribution first and then integrate it over a counterfactual distribution of the covariates. Distribution regression (Foresi and Peracchi 1995; Chernozhukov, Fernández-Val and Melly 2013) runs a separate binary regression at each threshold $y$,

$$F_{Y \mid X}(y \mid x) = \Phi\big(x'\gamma(y)\big),$$

which is the probit link, the one used on Q2 (logit replaces $\Phi$ with the logistic distribution function $\Lambda$), and here $x$ contains the FBM indicator, the number of listings, and their interaction. Setting each merchant's FBM indicator to 1 and then to 0 and averaging over the sample gives the two groups' distributions $\bar F_1$ and $\bar F_0$; the counterfactual distribution after the share rises by 0.1 is $(p + 0.1)\bar F_1 + (0.9 - p)\bar F_0$, and inverting it by linear interpolation between thresholds gives the counterfactual quantiles. It suits changes that are finite and not too small, and the theory of inference is in Theorem 5.2 and Corollary 5.3 of the arXiv version of Chernozhukov, Fernández-Val and Melly (2013). On Q2, with the share rising from twenty to thirty percent, the true discrete effect per unit of share at the 90th percentile is 0.632 and distribution regression gives 0.736; at the median the truth is 0.505 and the estimate 0.510. There are only point estimates here, and from this one estimate alone there is no telling how much of the gap between 0.736 and the truth is sampling error.

### 3.7 Endogenous Treatment and Instrumental Variables

On block Q3, a QR of log monthly revenue on AdBoost adoption and follower count gives an adoption coefficient of 0.466 at the median, with interval [0.434, 0.498], while the truth is only 0.10: the adopting merchants ranked higher to begin with. The vouchers were sent at random, and merchants with a voucher have an adoption rate 0.399 higher; this is an instrumental variable in the sense of Section 3.2 of Chapter 9. The question is what an instrumental variable identifies at a quantile.

A key assumption of instrumental variable quantile regression (IV-QR) is rank similarity.

::: {.assumption}
**IVQR.1 (structural quantile model; the linear case of condition A1 in the arXiv version of the Chernozhukov and Hansen 2013 survey, p. 5)** Given $X_i$, for each $d$, $Y_i(d) = d\,\beta_D(U_{id}) + X_i'\beta_X(U_{id})$, $U_{id} \mid X_i \sim \mathrm{U}(0,1)$, and the right-hand side is strictly increasing in $U_{id}$.

**IVQR.2 (independence of the instrument; A2)** Given $X_i$, $U_{id}$ is independent of $Z_i$.

**IVQR.3 (selection; A3)** $D_i$ is some unknown function of $(Z_i, X_i, V_i)$, where $V_i$ is an unobserved vector.

**IVQR.4 (rank similarity; A4)** Given $(X_i, Z_i, V_i)$, $U_{i0}$ and $U_{i1}$ have the same distribution.
:::

IVQR.4 is weaker than rank invariance ($U_{i0} = U_{i1}$): it allows ranks to slip randomly between the two states, but the slippage must not be systematically related to selection. When deciding whether to adopt, a merchant may know roughly where it ranks, but it cannot foresee exactly how its rank will change after adoption. IVQR.3 places almost no restriction on selection; the correlation between $V_i$ and $U_{id}$ is the source of endogeneity.

::: {.theorem}
**Theorem 7 (the IV-QR moment condition; Chernozhukov and Hansen 2005; arXiv version of the Chernozhukov and Hansen 2013 survey, Theorem 1 (ii) and equation (2.3), p. 6)** Under IVQR.1 through IVQR.4, for each $\tau \in (0,1)$,

$$\Pr\big(Y_i \leq D_i\,\beta_D(\tau) + X_i'\beta_X(\tau) \mid X_i, Z_i\big) = \tau.$$

Hence at $b = \beta_D(\tau)$, the conditional quantile of $Y_i - D_ib$ given $(X_i, Z_i)$ at level $\tau$ is $X_i'\beta_X(\tau)$, which does not depend on $Z_i$.
:::

The proof sketch has three steps. Step one (IVQR.1, IVQR.3, IVQR.4): write $U_i = U_{iD_i}$; the observed outcome is $Y_i = D_i\,\beta_D(U_i) + X_i'\beta_X(U_i)$; given $(X_i, Z_i, V_i)$, $D_i$ is already determined, and rank similarity makes $U_i$ distributed as $U_{i0}$. Step two (IVQR.2, IVQR.1): integrating over $V_i$, $\Pr(U_i \leq \tau \mid X_i, Z_i) = \Pr(U_{i0} \leq \tau \mid X_i, Z_i) = \Pr(U_{i0} \leq \tau \mid X_i) = \tau$. Step three (the strict monotonicity of IVQR.1): $\{Y_i \leq D_i\beta_D(\tau) + X_i'\beta_X(\tau)\} = \{U_i \leq \tau\}$, which has probability $\tau$.

In words: under IVQR.1 through IVQR.4, IV-QR estimates the structural quantile effect for all merchants, $\beta_D(\tau)$: how much revenue changes when a merchant of rank $\tau$ is switched from not adopting to adopting. This is not the compliers' effect; rank similarity makes the merchants moved by the instrument representative of all merchants. Theorem 7 gives only a moment condition; to solve it uniquely for $\beta_D(\tau)$, the instrument must actually move adoption, and a full-rank condition must guarantee that the solution is unique (the identification conditions of Chernozhukov and Hansen 2005). Otherwise the moment condition only marks out an identified set, and with a weak instrument the acceptance set can be very wide, which is why the confidence sets below are reported by inverting over a grid. Estimation uses the inverse QR of Chernozhukov and Hansen (2006, 2008): at each point $b$ of a grid, run a QR of $Y_i - D_ib$ on $(1, X_i, Z_i)$ and look for the point at which the coefficient on $Z_i$, as a function of $b$, is zero; the grid points at which the coefficient on $Z_i$ has a $t$ statistic no larger than the critical value form the confidence set, which does not rely on a strong first stage and follows the same idea as the AR set of Section 3.5 of Chapter 9. When rank similarity fails, this moment condition no longer corresponds to any structural parameter, and IV-QR estimates a quantity that is hard to interpret. The weights are set by the variation that the instrument moves; there is no explicit weighting. On Q3, IV-QR at the median is 0.130, with confidence set [0.048, 0.218] covering 0.10; at the 25th and 75th percentiles the accepted grid points each have a gap in the middle: the coefficient on $Z_i$ is a continuous piecewise linear function of $b$, and the nid standard error also varies with $b$, so the $t$ statistic is continuous but not monotone and, near the critical value, can cross it more than once. In that case the convex hull of the acceptance set is reported, [-0.070, 0.132] and [0.128, 0.288] respectively, covering 0.033 and 0.167; at the 90th percentile IV-QR is 0.202, and [0.072, 0.350] covers 0.228.

The compliers' effect takes another route, one that does not need rank similarity. The subscript $C$ denotes compliers, and $D_i(z)$ is the potential treatment when $Z_i = z$.

::: {.assumption}
**LQ.1 (independence and exclusion of the instrument; NBER working-paper version of Abadie, Angrist and Imbens 2002, Assumption 1)** $Z_i$ is randomly assigned, independent of $(Y_i(0), Y_i(1), D_i(0), D_i(1))$, and affects the outcome only through $D_i$.

**LQ.2 (first stage and monotonicity; same source)** $\Pr(D_i(1) = 1) > \Pr(D_i(0) = 1)$ and $D_i(1) \geq D_i(0)$: no merchant would adopt without a voucher yet not adopt with one.
:::

::: {.theorem}
**Theorem 8 (quantile treatment effects for compliers; IZA working-paper version of Frölich and Melly 2013, Theorem 1 and equation (3), p. 12; NBER working-paper version of Abadie, Angrist and Imbens 2002, Lemma 2 and Theorem 2)** Under LQ.1 and LQ.2, the share of compliers is $\Pr(C) = \Pr(D_i = 1 \mid Z_i = 1) - \Pr(D_i = 1 \mid Z_i = 0)$, and the distribution of the compliers' $Y_i(1)$ is identified by a Wald-type ratio:

$$\begin{aligned}
F_{Y(1) \mid C}(y) = \frac{1}{\Pr(C)}\Big[&\Pr(Y_i \leq y,\ D_i = 1 \mid Z_i = 1) \\
&- \Pr(Y_i \leq y,\ D_i = 1 \mid Z_i = 0)\Big],
\end{aligned}$$

$F_{Y(0) \mid C}$ replaces $D_i = 1$ by $D_i = 0$ and swaps the order of the two terms in the numerator, so $\mathrm{LQTE}(\tau) = F_{Y(1) \mid C}^{-1}(\tau) - F_{Y(0) \mid C}^{-1}(\tau)$. Equivalently, for any function $g$, $\mathbb{E}[\kappa_i\,g(Y_i, D_i, X_i)] = \Pr(C)\,\mathbb{E}[g(Y_i, D_i, X_i) \mid C]$, where

$$\kappa_i = 1 - \frac{D_i(1 - Z_i)}{1 - \Pr(Z_i = 1 \mid X_i)} - \frac{(1 - D_i)Z_i}{\Pr(Z_i = 1 \mid X_i)}.$$
:::

The proof sketch has three steps. Step one (LQ.2): monotonicity divides the merchants into three types, compliers, always-takers, and never-takers; a merchant that adopts with $Z_i = 0$ can only be an always-taker, and one that does not adopt with $Z_i = 1$ can only be a never-taker. Step two (LQ.1): random assignment makes the share of each type, and its outcome distribution, the same in the two instrument groups; $\Pr(Y_i \leq y, D_i = 1 \mid Z_i = 1)$ is the sum of an always-taker part and a complier part, while $\Pr(Y_i \leq y, D_i = 1 \mid Z_i = 0)$ has only the always-taker part; subtracting leaves only the compliers, and dividing by the compliers' share gives $F_{Y(1) \mid C}$, and the non-adopting side is analogous. $\kappa_i$ expresses the same subtraction as weights: observations whose $Z_i$ and $D_i$ disagree receive negative weights, which cancel out the noncompliers. Step three (Theorem 1): take generalized inverses of the two distribution functions; the sample Wald-type distribution function need not be monotone, but the generalized inverse is defined all the same.

In words: under LQ.1 and LQ.2, the LQTE estimates the difference between the compliers' two marginal distributions; it needs no rank similarity, and the price is that it concerns only the compliers, the forty percent of merchants moved by the vouchers; to read it as individual effects still requires the rank condition of Section 3.5. If the instrument is invalid, for example if the vouchers also affect revenue through some other channel, the subtraction leaves in a component from the noncompliers. The weights are set by the instrument: who the compliers are depends on which voucher is used, and a different instrument means a different set of compliers. On Q3 the LQTE at the median is 0.099, with interval [0.006, 0.193] covering 0.10; the Wald estimate of the LATE is 0.0996, with interval [0.0177, 0.1815]. IV-QR and the LQTE aim at the same number only at the median, a consequence of the symmetric design of Q3: at $\tau = 0.1$ the structural effect is -0.028 and the compliers' effect 0.004; at $\tau = 0.9$ they are 0.228 and 0.196, and the LQTE estimate of 0.271 (interval [0.151, 0.391]) aims at the latter. Away from the median, IV-QR cannot be taken as an estimate of the LQTE, or vice versa. Chapter 18 starts from LATE and treats compliers systematically.

### 3.8 Panel QR and Merchant Fixed Effects

On block Q4, pooled QR is far from the truth (Section 3.4). Chapter 10's remedy is within demeaning: expectations are linear, so subtracting each merchant's mean removes $\alpha_i$, and putting in a dummy variable for each merchant gives exactly the same result as the within transformation (Theorem 3 (iii) of Chapter 10). Quantiles have neither of these conveniences.

::: {.assumption}
**PQR.1 (conditional quantile with merchant fixed effects; Wooldridge 2010 equation (12.126), p. 460; the location shift of Canay 2011)** $\mathrm{Quant}_\tau(Y_{it} \mid X_{i1}, \ldots, X_{i\mathcal{T}}, \alpha_i) = \alpha_i + X_{it}'\beta(\tau)$, where $X_{it}$ contains a constant and the correlation between $\alpha_i$ and $X_{it}$ is unrestricted; Canay's estimator further requires $\alpha_i$ to be the same at every $\tau$, shifting only the location, whereas dummy-variable QR allows it to vary with $\tau$.

**PQR.2 (sampling; accepted manuscript of Kato, Galvao and Montes-Rojas 2012, conditions (A1) and (A2), p. 6)** Different merchants are independent of one another; for a given merchant, $(Y_{it}, X_{it})$ is independent and identically distributed across periods; $X_{it}$ has moments of sufficiently high order.
:::

::: {.theorem}
**Theorem 9 (panel QR; Wooldridge 2010 §12.10.3, p. 460; Koenker 2004; accepted manuscript of Kato, Galvao and Montes-Rojas 2012, Theorems 3.1 and 3.2, pp. 7-8; Canay 2011 Lemma 3.2 and Theorem 4.1, pp. 373-374)**

\(i) With the number of periods $\mathcal{T}$ fixed, estimating each $\alpha_i$ as a parameter together with $\beta(\tau)$ (one dummy variable per merchant) is in general inconsistent.

\(ii) Under PQR.1, PQR.2, and identification conditions, dummy-variable QR is consistent when $N/\mathcal{T}^{a} \to 0$, where the exponent $a$ is determined by the moment conditions on $X_{it}$; asymptotic normality centered at zero requires $N^2(\log N)^3/\mathcal{T} \to 0$.

\(iii) Under PQR.1 and Canay's conditions, $\mathcal{T} \geq 2$ is enough to identify $\beta(\tau)$; Canay's two-step estimator (first obtain $\hat\alpha_i$ from the within estimator of the mean model, then run a pooled QR on $Y_{it} - \hat\alpha_i$) is consistent and asymptotically normal when $N/\mathcal{T}^{a'} \to 0$ (for some $a' > 1$), but Besstremyannaya and Golovan (2019) point out that this rate condition is not strong enough.
:::

The proof sketch has three steps. Step one (PQR.1): quantiles are not additive, so the conditional quantile of $Y_{it} - \bar Y_i$ is not the conditional quantile of $Y_{it}$ minus a quantity involving only $\alpha_i$; Wooldridge (2010) p. 460 points out that no simple transformation removes $\alpha_i$. Step two (PQR.2): when $\alpha_i$ is treated as a parameter, each $\hat\alpha_i$ uses only $\mathcal{T}$ observations, and its error does not shrink with $N$; the first-order condition of the check function is nonlinear in $\alpha_i$, so this error does not cancel in the average but passes into $\hat\beta(\tau)$, creating a bias that vanishes only with $\mathcal{T}$; this is the incidental parameters problem discussed in Section 3.8 of Chapter 12, and the expansion of Kato, Galvao and Montes-Rojas gives the order of this bias, and from it the rate condition of (ii). Step three (the location shift of PQR.1): $\alpha_i$ is the same at every quantile, so it is also the merchant effect in the conditional mean; the within estimator of Chapter 10 estimates the mean slope consistently, and $\hat\alpha_i$ is the within-merchant average of $Y_{it} - X_{it}'\hat\beta$; its error vanishes only as $\mathcal{T} \to \infty$, so (iii) likewise requires $\mathcal{T}$ to grow relative to $N$.

In words: under PQR.1, both estimators aim at $\beta(\tau)$, the slope of the conditional quantile given the merchant effect; if the merchant effect also changes the shape of the distribution (for example, if large merchants' revenue is more dispersed), Canay's estimator no longer aims at it, while dummy-variable QR still does, though it depends more on a long panel. No asymptotic theorem covers Q4: $N = 200$ while $\mathcal{T}$ is only 5 or 20, far from $N^2(\log N)^3/\mathcal{T} \to 0$, and Canay's rate condition is itself disputed, so Q4 is only a finite-sample Monte Carlo illustration, not a verification of any theorem. Over 1,000 replications at $\tau = 0.9$, the average dummy-variable QR estimate is 0.6091 at $\mathcal{T} = 5$ and 0.7573 at $\mathcal{T} = 20$; Canay's two-step estimates average 0.6903 and 0.7581; the truth is 0.7845, while the pooled QR average lies far above the truth (Section 3.4). The bias is concentrated in the tails: at the median, the dummy-variable QR average at $\mathcal{T} = 5$ is 0.4008, Canay's is 0.4022, and the truth is 0.4. The weights are set by each merchant's own number of periods: in dummy-variable QR each $\alpha_i$ is located only by its own $\mathcal{T}$ months; with short panels the tail quantiles rest on very few observations, which is where the bias comes from. A single sample says nothing about bias: in this particular draw, dummy-variable QR at $\mathcal{T} = 5$ gives 0.778 and Canay at $\mathcal{T} = 20$ gives 0.851, so here the former is the closer to the truth, the reverse of the Monte Carlo averages; the standard errors of the two Monte Carlo means are only 0.0039 and 0.0017, so the gap between the averages is not simulation noise.

### 3.9 Matching by Rank: Quantile Effects for the Treated and for Switchers

On block Q5, the most direct approach is quantile DiD, a DiD run level by level: at each level, subtract the difference between the control region's quantiles in the two seasons from the difference between the launch region's quantiles in the two seasons. At the median it gives 3.575, with interval [2.862, 4.289], while the true QTT is only 1.920. Athey and Imbens (2006) give the model under which quantile DiD is valid: the untreated outcome is additively separable in group and time, and the unobservable is independent of both group and time (working-paper version, Theorem 3.2); in Q5, peak-season growth increases with size and merchants in the launch region are larger, so quantile DiD credits One-Click Checkout with the launch region's inherently faster peak-season growth.

Changes-in-changes (CIC) takes a different approach. It matches by rank, but ranks are matched only between the two periods of the control region; across regions, it is the value of the outcome that is matched, not the rank. Let $Y_{st}$ be the outcome in the cell of region $s$ and period $t$ (the first subscript is the region, the second the period), and let $F_{Y_{st}}$ be its distribution function.

::: {.assumption}
**CIC.1 (the model; Athey and Imbens 2006 working-paper version, Assumption 3.1, p. 8)** The untreated outcome is $Y_i(0) = h(U_i, T_i)$, where $U_i$ is a one-dimensional unobservable.

**CIC.2 (strict monotonicity; Assumption 3.2)** For $t = 0, 1$, $h(\cdot, t)$ is strictly increasing.

**CIC.3 (time invariance; Assumption 3.3)** $U_i \perp T_i \mid S_i$: within each region, the distribution of $U_i$ is the same in both periods.

**CIC.4 (support; Assumption 3.4)** The support of $U_i$ in the launch region is contained in its support in the control region.
:::

Taken together, CIC.1 through CIC.3 say that merchants with the same $U_i$ change over time in the same way, whichever region they are in; the distributions of $U_i$ in the two regions can differ arbitrarily, so two such merchants generally hold different ranks in their own regions. The model replaces parallel trends with "people with the same $U_i$ change in the same way over time." CIC.3 rules out compositional change within a region, such as a wave of new merchants pouring into the control region in the peak season.

::: {.theorem}
**Theorem 10 (CIC identification and the QTT; Athey and Imbens 2006 working-paper version, Theorem 3.1, equations (3.9) and (3.15), and the quantile effects that follow them, pp. 10-11)** Under CIC.1 through CIC.4, write $Q_{d=0}(y) = F_{Y_{01}}^{-1}\big(F_{Y_{00}}(y)\big)$; then

$$\begin{aligned}
F_{Y_{11}(0)}(y) &= F_{Y_{10}}\Big(F_{Y_{00}}^{-1}\big(F_{Y_{01}}(y)\big)\Big), \\
\mathrm{QTT}(\tau) &= F_{Y_{11}}^{-1}(\tau) - Q_{d=0}\big(F_{Y_{10}}^{-1}(\tau)\big).
\end{aligned}$$
:::

The proof sketch has three steps. Step one (CIC.1 through CIC.3): in the control region $F_{Y_{0t}}(y) = F_{U \mid S = 0}\big(h^{-1}(y, t)\big)$ (working-paper version, equation (3.10)), so $Q_{d=0}$ maps the outcome at rank $\tau$ in the base period to the outcome at the same rank in the later period, $Q_{d=0}(y) = h\big(h^{-1}(y, 0), 1\big)$. Step two (CIC.2, CIC.4): a launch-region merchant with base-period outcome $y$ and a control-region merchant with base-period outcome $y$ have the same $U_i$, and the support condition guarantees that such a control merchant exists; without the launch, its later-period outcome would be $h(U_i, 1) = Q_{d=0}(y)$, so $Y_{11}(0)$ has the same distribution as $Q_{d=0}(Y_{10})$. Step three (Theorem 1 (iii)): $Q_{d=0}$ is strictly increasing, so the quantile of $Q_{d=0}(Y_{10})$ at level $\tau$ is $Q_{d=0}(F_{Y_{10}}^{-1}(\tau))$, and subtracting it from $F_{Y_{11}}^{-1}(\tau)$ gives the QTT.

In words: under CIC.1 through CIC.4, CIC estimates the QTT of the launch region. It still compares two distributions, the launch region's actual later-period distribution and the counterfactual distribution carried over by rank; to read it as individual effects, we must add the rank condition of Section 3.5. CIC does not depend on the measurement scale of the outcome: take logs of revenue, and $Q_{d=0}$ is transformed along with it, leaving the conclusions about the QTT unchanged; DiD does, since parallel trends in revenue and parallel trends in log revenue are two different assumptions, and DiDs on the two scales can be far apart (the Athey and Imbens application in Section 5). When CIC.3 fails, for example when a wave of new merchants pours into the control region in the peak season, $Q_{d=0}$ treats the compositional change as a time effect. As for the weights, there are no explicit ones here: each treated merchant is matched to the control merchants with the same base-period outcome, hence the same $U_i$, and takes the later-period value at their rank in the control group; the control group's distributions in the two periods determine the shape of the mapping. On Q5, CIC at the median gives 1.744, with interval [0.610, 2.879] covering 1.920; quantile DiD's 3.575 aims in the population at 3.909, exceeding the QTT by 1.989, and the excess is the part due to the launch region's merchants being larger and growing faster in the peak season. At the five levels, CIC's intervals all cover the QTT; quantile DiD covers only at the 10th percentile ([0.543, 1.530] against 0.937), where its population bias is smallest, at 0.704, and its interval is wide enough.

![Block Q5: the QTT of One-Click Checkout on monthly revenue (thousands of USD) in the launch region. Circles are CIC estimates and triangles quantile DiD estimates; vertical bars are pointwise 95% bootstrap intervals, resampling within each region-by-quarter cell, at levels from 5% to 95%; the solid line is the true QTT.](assets/fig/fig_13_cic.svg)

At the 19 plotting levels, the pointwise intervals of CIC cover the truth at 18; the one miss is at the 20th percentile, with estimate 0.586, interval [0.027, 1.146], and truth 1.198; one miss among 19 levels is in line with the nominal level of pointwise 95% intervals. Quantile DiD covers at only 5.

In Q5f adoption is voluntary: in the launch region some merchants adopted long ago and some never adopt. de Chaisemartin and D'Haultfœuille (2018, hereafter dCDH) extend CIC to such fuzzy designs, in which the treatment rate changes only partly within a group, with switchers as the target. Carry out the mapping separately for each treatment status $d$,

$$Q_d(y) = F_{Y_{01} \mid D = d}^{-1}\big(F_{Y_{00} \mid D = d}(y)\big),$$

that is, among the control-region merchants with treatment status $d$, map the outcome at rank $\tau$ in the base period to the outcome at the same rank in the later period; for $d = 0$ this is the $Q_{d=0}$ above.

The idea behind the identification is this: the later-period adopters in the launch region consist of the early adopters and the switchers; first use the $d = 1$ version of $Q_d$ to carry the early adopters' base-period distribution over to the later period, then subtract it, by shares, from the distribution of the later-period adopters, and what remains is the distribution of the switchers' $Y_i(1)$; the non-adopting side is analogous.

::: {.assumption}
**FZ.1 (treatment participation; dCDH 2018 Assumption 3)** $D_i = \mathbf{1}\{V_i \geq v_{S_iT_i}\}$, with the threshold varying by region and period, and $V_i \perp T_i \mid S_i$.

**FZ.2 (first stage; Assumption 1)** The adoption rate in the launch region rises, and rises more than in the control region.

**FZ.3 (stable adoption rate in the control group; Assumption 2)** $0 < \mathbb{E}(D_{01}) = \mathbb{E}(D_{00}) < 1$.

**FZ.4 (monotonicity and time invariance; Assumption 7)** $Y_i(d) = h_d(U_{id}, T_i)$, with $h_d$ strictly increasing in $U_{id}$; given the region and the base-period treatment status, $U_{id}$ is independent of $T_i$.

**FZ.5 (data conditions; Assumption 8)** The outcome has the same support, a closed interval, in each of the eight treatment-by-region-by-period cells, and the distribution function of each cell is continuous and strictly increasing on the support.
:::

::: {.theorem}
**Theorem 11 (the switchers' LQTE; dCDH 2018 Theorem 3)** Under FZ.1 through FZ.5, let

$$\begin{aligned}
&\big[\Pr(D_{11} = d) - \Pr(D_{10} = d)\big]\,F_{\mathrm{CIC},d}(y) \\
&\quad = \Pr(D_{11} = d)\,F_{Y_{11} \mid D = d}(y) \\
&\qquad - \Pr(D_{10} = d)\,F_{Q_d(Y_{10}) \mid D = d}(y),
\end{aligned}$$

where $F_{Q_d(Y_{10}) \mid D = d}$ is the distribution of the launch-region merchants with base-period treatment status $d$ after they are mapped by $Q_d$. Then $\mathrm{LQTE}_{sw}(\tau) = F_{\mathrm{CIC},1}^{-1}(\tau) - F_{\mathrm{CIC},0}^{-1}(\tau)$, and the switchers' LATE equals

$$\frac{\mathbb{E}(Y_{11}) - \mathbb{E}[Q_{D_i}(Y_i) \mid S_i = 1, T_i = 0]}{\mathbb{E}(D_{11}) - \mathbb{E}(D_{10})}.$$
:::

The proof sketch has three steps. Step one (FZ.1, FZ.3): the control region's adoption rate does not change, and with $V_i \perp T_i \mid S_i$, the cell with treatment status $d$ contains the same merchants in both periods; FZ.4 makes $Q_d$ map the outcome at rank $\tau$ among these merchants in the base period to the outcome at the same rank in the later period. Step two (FZ.4, FZ.5): a launch-region merchant with base-period treatment status $d$ and outcome $y$ has the same $U_{id}$ as a control-region merchant with the same status and the same outcome, and if it keeps its original status its later-period outcome is $Q_d(y)$, so $F_{Q_d(Y_{10}) \mid D = d}$ is the distribution of these merchants' later-period $Y_i(d)$. Step three (FZ.2): the launch-region merchants who adopt in the later period consist of switchers and early adopters; subtracting, by shares, the reconstructed distribution of the early adopters from the distribution of the later-period adopters and dividing by the switchers' share gives the distribution of the switchers' $Y_i(1)$, by the same subtraction as in Theorem 8; the non-adopting side (with a negative denominator) gives the distribution of the switchers' $Y_i(0)$; taking generalized inverses and subtracting gives the LQTE, and taking means and subtracting gives the LATE.

In words: under FZ.1 through FZ.5, fuzzy CIC estimates the quantile treatment effects of the launch region's switchers; like the compliers of Theorem 8, this is the effect for the merchants moved by the policy. It requires the control region's adoption rate not to change: in Q5f the rates in the two periods are 0.1515 and 0.1475, while in the population both are 0.15, so the sample difference is sampling error; the launch region rises from 0.1027 to 0.705, a first stage of 0.602. When FZ.3 fails, the control cells no longer contain the same merchants in both periods, compositional change seeps into the mapping, and dCDH show that only bounds can then be given. The weights are set by the shares: the switchers' distribution is a difference of two distributions weighted by adoption shares, and the smaller the share, the more the subtraction amplifies noise. At the median of Q5f, the switchers' LQTE estimate is 2.158, with interval [1.225, 3.091] covering 2.047; in the tails there are only point estimates, 1.372 against a truth of 1.045 at the 10th percentile and 6.100 against 4.012 at the 90th. A distribution function obtained by subtraction is least stable in the tails, and without intervals there is no way to judge how much of these gaps is sampling error. The same mapping gives a switchers' LATE of 3.016, with interval [2.078, 3.953] covering 2.350; Chapter 20 calls it Wald-CIC and uses the same $Q_d$ to discuss fuzzy DiD.

## 4 Estimation and Practice

### 4.1 Computation and the Convention for Sample Quantiles

The sample problem of QR is a linear program: split each residual into its positive and negative parts; the objective is their weighted sum, and the constraints are linear. The Barrodale-Roberts simplex method, the default in `quantreg` (`method = "br"`), is fast with a few thousand observations; with many observations or a sparse design matrix, use the Frisch-Newton interior point method (`"fn"`, whose sparse version is `rq.fit.sfn`), which is how the QR with one dummy variable per merchant in Section 6 is solved. The solution of a linear program lies at a vertex, with exactly $k$ residuals equal to zero. When the only regressor is a binary variable and $n\tau$ in some group happens to be an integer, the minimizer is an interval and the software reports only one of its endpoints; the two experimental groups of Marigold each have an odd number of merchants to avoid this nonuniqueness. Newton's method is of no use: the objective is not differentiable, and the Hessian is zero everywhere.

The convention for sample quantiles must also be stated. Throughout this chapter, a sample quantile is the generalized inverse of the empirical distribution function, that is, the smallest order statistic with $j/n \geq \tau$. R's `quantile(type = 1)` usually gives the same value, but when $n\tau$ is slightly larger than an integer because of floating-point error, it goes one order statistic too far, so Section 6 computes quantiles from the definition; only the RIF step uses `type = 7`, to stay consistent with `rifreg`, and distribution regression interpolates linearly between thresholds.

### 4.2 Inference

The options of `summary.rq(se = )` correspond to different ways of estimating the variance of Theorem 3:

| `se =` | How the bread $A_0(\tau)$ is estimated | When it is valid |
|---|---|---|
| `"iid"` | One common density of the residuals at zero | The error density at the quantile does not vary with $X_i$ (equations (12.119) to (12.121)) |
| `"nid"` | A local sparsity for each observation, with the bandwidth of Hall and Sheather (1988) | QR.1 holds; heteroskedasticity allowed |
| `"ker"` | Kernel weighting of the residuals, of the form of equation (12.118) | Same as above |
| `"rank"` | No bread is estimated; intervals come from inverting a rank test (the default when $n < 1001$) | Errors independent and identically distributed, or a relaxed version of this |
| `"boot"` | No formula; the pairs bootstrap (`bsmethod = "xy"`), which can resample by unit | Does not rely on the $\tau(1-\tau)$ filling |

The filling of the first three is $\tau(1-\tau)X'X/n$, that is, $B_0(\tau)$ under correct specification; when the specification may be wrong, the warning after Theorem 3 applies, and the bootstrap is safer. The bandwidth $h_n$ must satisfy $h_n \to 0$ and $\sqrt n\,h_n \to \infty$; the example Wooldridge (2010) p. 457 gives is a constant multiple of $n^{-1/3}$, and the Hall-Sheather bandwidth is of this order. The QR objective is not smooth, but the estimator itself is asymptotically normal, and the pairs bootstrap of Section 3.3 of Chapter 8 still applies to it (Wooldridge 2010 p. 457); with panel and clustered data, resample by unit, or use the cluster-robust estimators of equations (12.124) and (12.125).

The other estimators each have their own form of inference. If RIF regression treats the RIF as data, the OLS sandwich misses the estimation error of the two steps, the quantile and the density: at the 90th percentile of Q2, the bootstrap standard error of the whole procedure is 0.0527, while treating the RIF as data gives 0.0435. The IV-QR confidence set relies on inverting over a grid, and the grid must be fine enough and wide enough; when the acceptance set touches the edge of the grid, widen the range; Section 6 uses [-0.30, 0.60] with a step of 0.002. CIC and fuzzy CIC resample within each region-by-quarter cell, consistent with how the repeated cross-sections are sampled. The intervals in this chapter are all pointwise: the probability that all 19 levels are covered at once is below 95%, and drawing conclusions about the whole curve requires uniform confidence bands of the kind in Angrist, Chernozhukov and Fernández-Val (2006).

### 4.3 The Estimators at a Glance

| Estimand | Key assumption | Estimator | Marigold |
|---|---|---|---|
| Parameters of the conditional quantile $\beta(\tau)$ | QR.1 | `rq` | Q2, 90th percentile: 0.523 |
| QTE | Random assignment or unconfoundedness | Difference of the two groups' quantiles (possibly weighted) | Q1, 90th percentile: 2.750 |
| Distribution of individual effects | Add rank invariance | Same as above | Cannot be tested from the data |
| UQPE | A change in the share does not change the conditional distributions | RIF regression | Q2, 90th percentile: 0.723 |
| Structural effect $\beta_D(\tau)$ | Rank similarity | Inverse QR | Q3, median: 0.130 |
| LQTE | Monotonicity, exclusion | Wald-type distribution functions or $\kappa_i$ weighting | Q3, median: 0.099 |
| Panel $\beta(\tau)$ | Canay needs a location shift; both need $\mathcal{T}$ large relative to $N$ | Canay two-step, dummy-variable QR | Q4: Monte Carlo only |
| QTT | CIC.1 through CIC.4 | CIC | Q5, median: 1.744 |
| $\mathrm{LQTE}_{sw}$ | FZ.1 through FZ.5 | Fuzzy CIC | Q5f, median: 2.158 |

To choose a row, first answer "which position of which distribution do I want to know?", and only then ask whether the data offer a randomized experiment, an instrument, a panel, or a before-and-after policy comparison.

### 4.4 Software

`quantreg` 6.1 provides `rq`, `summary.rq(se = )`, `rq.fit.br`, the sparse `rq.fit.sfn`, and `boot.rq`; `conquer`, which fits QR on large samples with a smoothed check function, can serve as a complement. `rifreg` 1.1.0 implements RIF regression, and it gives standard errors that account for the two-step error only with `bootstrap = TRUE`. In `qte` 2.0.0, `unc_qte` computes the QTE under random assignment and under unconfoundedness, and `cic` computes CIC (write `gt_type = "qtt"` and supply `probs`). When checking results against these packages, first align the quantile conventions, and do not expect agreement to the last digit: `unc_qte` finds quantiles by numerical optimization, close to the generalized inverse but not identical digit for digit; the new `cic()` uses another quantile convention and differs visibly from this chapter's CIC on Q5; the deprecated `CiC()`, when the cells are of equal size, takes the wrong order statistic at a few levels because of floating-point error. In Section 6, IV-QR, distribution regression, Canay's two-step estimator, and fuzzy CIC are all coded as functions directly from their definitions.

## 5 Anchoring Papers

Each of the three anchoring papers corresponds to one thread of this chapter: the comparison of distributions and individual effects in randomized experiments, conditional versus unconditional quantiles, and the two-group, two-period design that matches by rank. In IS and marketing journals we have not yet found a verifiable causal application of QR, QTE, or CIC, so all three come from economics journals, and each is laid out in five parts: paper, method, data, results, and limitations.

### 5.1 Bitler, Gelbach and Hoynes (2006)

::: {.case}
Paper: "What Mean Impacts Miss: Distributional Effects of Welfare Reform Experiments," American Economic Review. The numbers cited follow the published version (identical to the final pre-publication version the authors posted, dated August 23, 2005); the subsection and page numbers cited in Theorem 5 follow its NBER Working Paper 10121 (November 2003).

Method: a randomized experiment of Connecticut's Jobs First welfare reform; the paper compares the treatment and control groups' income distributions by the difference at the same quantile, that is, the QTE of Theorem 5 (i) of this chapter; the two groups' earnings before random assignment are unbalanced, and the authors adjust by inverse propensity score weighting, with intervals from a bootstrap that uses individuals as blocks.

Data: the experiment's public-use sample has 4,803 cases in all, 2,396 assigned to Jobs First and 2,407 remaining in the old AFDC system; income is recorded quarterly, and the unit is the person-quarter.

Results: over 16 quarters, the average effect of Jobs First on total quarterly income (earnings plus transfers) is 136 dollars, against a control-group baseline of 2,609 dollars, about 5%. Now take quarterly earnings before the time limit takes effect (quarters 1 to 7): 48% of the treatment group's person-quarters and 55% of the control group's are zero, so the QTE is identically zero at nearly half of the quantiles; the treatment group is higher from the 49th to the 82nd percentile, the two are equal from the 83rd to the 87th, and the control group is higher from the 88th to the 97th (the authors report only up to the 97th percentile). The point estimates of the earnings QTE range roughly from minus 300 to 500 dollars, while the average effect on earnings over the same period is only 82 dollars (standard error 52 dollars).

Limitations: the authors stress that the QTE is not necessarily the effect for identifiable individuals: if Jobs First reverses ranks, the difference between the medians of the two distributions does not give the effect for the person originally at the median; the working-paper version further notes that ranks are probably not preserved at the low end of the distribution (some respond to the incentives, some do not). So the QTE can be read only as a difference between two distributions. They list conclusions that do not depend on ranks, for example that a negative QTE somewhere implies that at least someone is hurt, which is Theorem 5 (iii) of this chapter. With the person-quarter as the unit, the correlation across quarters for the same person is handled only in the bootstrap.
:::

This paper rereads an experiment with a small mean effect as "some gain and some lose," with the same structure as Q1, and its caution matches this chapter's: it makes clear how the distribution moves, but does not say who moves.

### 5.2 Firpo, Fortin and Lemieux (2009)

::: {.case}
Paper: "Unconditional Quantile Regressions," Econometrica. The numbers cited follow the published version (Table I on p. 964, the sample description on pp. 962-963); the proposition numbers follow its NBER Technical Working Paper 339 (July 2007), since the published version renumbered the propositions.

Method: RIF regression, that is, Theorem 6 of this chapter: regressing the recentered influence function on the covariates by OLS gives the effect of a small change in the distribution of the covariates on the unconditional quantiles.

Data: male wages from the 1983 to 1985 CPS, 266,956 observations; the outcome is the log hourly wage, and the focus is union coverage.

Results: OLS gives a union coverage coefficient of 0.179 (standard error 0.002). The unconditional quantile effects at the 10th, 50th, and 90th percentiles are 0.195, 0.337, and -0.135 (the published text still uses the working-paper numbers; the values here are taken from the table); conditional QR at the same three positions gives 0.288, 0.195, and 0.088, declining monotonically across quantiles. The authors' explanation is that unions raise the location of the conditional distribution and compress its dispersion: at the low end the two forces point the same way, at the high end in opposite ways.

Limitations: RIF regression gives a first-order approximation for a small change in the distribution of the covariates, and when union coverage changes a lot, the approximation error depends on how far apart the new and old distributions are; it assumes that the change does not alter the wage distribution given the covariates, so there are no general equilibrium effects. The OLS of the RIF is a linear approximation to a conditional probability, and the authors compare it with a logit version; the two are very close.
:::

At the 90th percentile, the conditional effect 0.088 and the unconditional effect -0.135 have opposite signs: holding the other covariates fixed, union coverage moves the 90th percentile of the conditional wage distribution up; yet a higher union coverage rate pushes down the upper end of the whole wage distribution. This is the same divergence between the conditional and unconditional effects as in Q2, only starker.

### 5.3 Athey and Imbens (2006)

::: {.case}
Paper: "Identification and Inference in Nonlinear Difference-in-Differences Models," Econometrica. The numbering of assumptions, theorems, and equations follows its NBER Technical Working Paper 280 (September 2002); the published version moved this application to the online supplement (Athey and Imbens 2006, Supplementary Material §3), and the numbers cited follow the supplement.

Method: CIC, that is, Theorem 10 of this chapter; with a discrete outcome, its discrete version and bounds are used, and compared with two DiDs, one in levels and one in logs.

Data: the increase in workers' compensation benefits in Kentucky studied by Meyer, Viscusi and Durbin (1995): the higher benefit cap benefited only high-earning workers (the treatment group), and low-earning workers serve as the control group; the outcome is the number of weeks receiving benefits, counted in whole weeks and heavily right-skewed. Meyer, Viscusi and Durbin themselves found that the effect is significant when the outcome is log weeks and not significant when it is weeks.

Results: for the treatment group, the average effects on log weeks of the log DiD and of discrete CIC (with a conditional independence assumption added) are 0.191 and 0.183, both with standard error 0.07 in the table, and the authors regard the two as comparable; at the 75th percentile, the discrete CIC effect is 2 weeks (standard error 0.77). For reference, the treatment group's mean of log weeks before the increase is 1.38, and its 75th percentile is 8 weeks. The authors do not report the effect of the levels DiD on log weeks: a counterfactual built by shifting levels would produce negative weeks, whose logs cannot be taken. At the 25th percentile, the levels DiD gives -0.77 weeks (standard error 0.58), whereas the lower bound from discrete CIC is 0; the authors use this contrast to show that the choice of model can change the conclusion.

Limitations: weeks are discrete, so CIC needs a conditional independence assumption for point identification and otherwise gives only bounds; CIC requires the distribution of the unobservables within each group to be the same in both periods, and the authors argue that the characteristics of injured workers would not change in the short run because of the policy change.
:::

This application shows CIC's selling point most clearly: the DiD conclusion depends on whether the outcome is measured in weeks or in logs, while CIC is invariant to monotone transformations (Theorem 1 (iii) and Theorem 10), so the researcher need not decide in advance on which scale parallel trends holds.

Taken together: Bitler, Gelbach and Hoynes show that in a randomized experiment a difference in quantiles is a comparison of distributions, not anyone's effect; Firpo, Fortin and Lemieux show that the effect controlling for covariates and the effect on the whole distribution can have opposite signs; Athey and Imbens show that in before-and-after policy comparisons, matching by rank depends less on the scale than shifting by means. The discipline common to all three is the main thread of this chapter: first say clearly which position of which distribution you want to look at.

## 6 A Full Walkthrough on Marigold Data

The code below runs under R 4.6.1. It lists only the statements that generate the data, estimate, compute variances, and summarize the Monte Carlo; several statements are merged onto one line, and statements used only for checks are omitted. Every number in the text comes from the output of an actual run, and the result lines in the code that start with `#` summarize the numbers printed by the run. The code uses `quantreg`, `sandwich`, and `SparseM`. The helper functions not shown are: `gen_inverse(support, cdf, tau)`, the generalized inverse of an arbitrary distribution function (not necessarily monotone); `fe_dummy_qr`, QR with one dummy variable per merchant, solved by the sparse Frisch-Newton interior point method; `dist_reg` and `dr_cdf`, the threshold-by-threshold probit and the distribution function averaged over the sample; `cic_lqte`, which computes the switchers' LQTE and LATE from the $F_{\mathrm{CIC},d}$ of Theorem 11, with a bootstrap within cells; `.muffle_nonunique`, which suppresses only the message of `quantreg` that the solution may be nonunique.

### 6.1 The DGP

Q1 and Q1r share the merchants, the covariates, the treatment assignment, and $Y_i(0)$; Q1r uses new normal draws to re-pair $Y_i(1)$ with $Y_i(0)$, leaving the marginal distributions unchanged. The panel of Q4 is also written out in full; the other blocks are generated from the formulas of Section 2.2.

```r
taus <- c(0.1, 0.25, 0.5, 0.75, 0.9)
taus_fine <- round(seq(0.05, 0.95, by = 0.05), 2)               # the 19 levels of the figures
P1 <- list(n = 5002, n_treat = 2501, m = log(10), b = 0.4, s = 0.3, k = 1.4, rho_r = 0.5)
P1$v <- P1$b^2 + P1$s^2                                        # variance of log Y(0)
P1$logc <- (1 - P1$k) * P1$m + (1 - P1$k^2) * P1$v / 2          # makes E[Y(1)] = E[Y(0)]: ATE = 0
gen_Q1 <- function(P) {                                         # AutoPrice, monthly revenue in USD thousands
  n <- P$n
  x <- rnorm(n); e <- rnorm(n)
  D <- sample(rep(0:1, c(n - P$n_treat, P$n_treat)))           # complete randomisation, 2,501 per arm
  ly0 <- P$m + P$b * x + P$s * e
  y0 <- exp(ly0); y1 <- exp(P$logc + P$k * ly0)                # Y(1) = c Y(0)^k: ranks are kept
  data.frame(merchant = seq_len(n), x = x, D = D, y = ifelse(D == 1, y1, y0), y0_true = y0, y1_true = y1)
}
gen_Q1r <- function(P, Q1) {                                    # same merchants and marginals, Y(1) re-paired
  xi <- rnorm(nrow(Q1))
  z0 <- (log(Q1$y0_true) - P$m) / sqrt(P$v)
  y1 <- exp(P$logc + P$k * (P$m + sqrt(P$v) * (P$rho_r * z0 + sqrt(1 - P$rho_r^2) * xi)))
  data.frame(merchant = Q1$merchant, x = Q1$x, D = Q1$D, y = ifelse(Q1$D == 1, y1, Q1$y0_true),
             y0_true = Q1$y0_true, y1_true = y1)
}
gen_Q4 <- function(P, Tn) {                                     # merchant panel, Tn months, log revenue
  N <- P$N
  id <- rep(seq_len(N), each = Tn)
  a <- rnorm(N, 0, P$sa)[id]                                    # merchant effect alpha_i
  x <- plogis(P$cx * a / P$sa + rnorm(N * Tn))                  # promotion share, larger for larger merchants
  e <- rnorm(N * Tn)
  data.frame(id = id, t = rep(seq_len(Tn), N), x = x, y = P$mu + a + P$beta * x + P$sig * (1 + P$gam * x) * e,
             alpha_true = a)
}
P2  <- list(n = 6002, n_fbm = 1201, p = 0.2, mu = log(8000), g = 0.5, bx = 0.5, s0 = 0.35, s1 = 0.35)
P3  <- list(n = 6002, n_z1 = 3001, mu0 = log(8000), dmu = 0.10, s0 = 0.5, s1 = 0.6, bx = 0.3, rho = 0.5,
            pi0 = 0.3, pi1 = 0.7)
P4  <- list(N = 200, mu = log(8000), sa = 0.6, cx = 0.5, beta = 0.4, sig = 0.3, gam = 1)
P5  <- list(nc = 2500, A = 10, g = 0.05, s = 0.5, kap = 0.4, mu1 = 0.6, sd1 = 0.8, delta = 0.12)
P5f <- modifyList(P5, list(nc = 10000, rho = 0.4, p_at1 = 0.10, p_sw1 = 0.60, p_at0 = 0.15))
# gen_Q2, gen_Q3, gen_Q5 and gen_Q5f draw the other blocks from these parameters as described in section 2.2
set.seed(1301); Q1 <- gen_Q1(P1); set.seed(1302); Q1r <- gen_Q1r(P1, Q1)
set.seed(1303); Q2 <- gen_Q2(P2); set.seed(1304); Q3 <- gen_Q3(P3)
set.seed(1305); Q4_T5 <- gen_Q4(P4, 5); set.seed(1306); Q4_T20 <- gen_Q4(P4, 20)
set.seed(1307); Q5 <- gen_Q5(P5); set.seed(1308); Q5f <- gen_Q5f(P5f)
```

In Q3, `rho` is the correlation between rank and the propensity to adopt, and `pi0` and `pi1` are the adoption rates without and with a voucher; in Q5, `delta` is the revenue increase brought by One-Click Checkout; the three shares of Q5f are, respectively, the shares of early adopters and of peak-season switchers in the launch region, and the adoption rate in the control region.

### 6.2 Blocks Q1 and Q1r: Means, Differences in Quantiles, and Conditional QR

`emp_quantile` takes the generalized inverse by definition; `qte_marginal` bootstraps by resampling with replacement within each of the two groups.

```r
library(quantreg)
emp_quantile <- function(y, tau) {          # F_n^{-1}(tau) = inf{y : F_n(y) >= tau}, the generalized inverse
  s <- sort(as.numeric(y)); n <- length(s)
  j <- ceiling(n * tau)
  j <- j - ((j - 1) / n >= tau); j <- j + (j / n < tau)        # the smallest j with j / n >= tau, exactly
  s[pmin(pmax(j, 1), n)]
}
qte_marginal <- function(y, d, taus, B = 999, level = 0.95) {   # difference of the two arms' quantiles
  y1 <- y[d == 1]; y0 <- y[d == 0]
  est <- emp_quantile(y1, taus) - emp_quantile(y0, taus)
  bs <- vapply(seq_len(B), function(b) emp_quantile(y1[sample.int(length(y1), replace = TRUE)], taus) -
                 emp_quantile(y0[sample.int(length(y0), replace = TRUE)], taus), numeric(length(taus)))
  se <- apply(matrix(bs, nrow = length(taus)), 1, sd)           # bootstrap that resamples each arm
  z <- qnorm(1 - (1 - level) / 2)
  data.frame(tau = taus, est = est, lo = est - z * se, hi = est + z * se, se = se)
}
ols_hc2 <- function(fit) {                                      # OLS coefficient on D, HC2 standard error
  est <- coef(fit)[["D"]]; se <- sqrt(sandwich::vcovHC(fit, type = "HC2")["D", "D"])
  round(c(est = est, se = se, lo = est - qnorm(0.975) * se, hi = est + qnorm(0.975) * se), 3)
}
qte_true <- function(t) { q0 <- exp(P1$m + sqrt(P1$v) * qnorm(t)); exp(P1$logc) * q0^P1$k - q0 }
cqte_true <- function(t) P1$logc + (P1$k - 1) * (P1$m + P1$s * qnorm(t))   # in logs, at X = 0
qlog_true <- function(t) P1$logc + (P1$k - 1) * (P1$m + sqrt(P1$v) * qnorm(t))   # in logs, marginal
cqr_D <- function(d) sapply(taus, function(t)                   # QR of log Y on D, X, D X: coefficient on D
  summary(rq(log(y) ~ D * x, tau = t, data = d), se = "nid")$coefficients["D", 1])
ols_hc2(lm(y ~ D, data = Q1)); ols_hc2(lm(y ~ D * x, data = Q1))
set.seed(13101); q1 <- qte_marginal(Q1$y, Q1$D, taus_fine)
set.seed(13102); q1r <- qte_marginal(Q1r$y, Q1r$D, taus_fine)
round(cbind(q1[q1$tau %in% taus, ], truth = qte_true(taus)), 3)
cover <- function(q) q$lo <= qte_true(taus_fine) & qte_true(taus_fine) <= q$hi
c(Q1 = sum(cover(q1)), Q1r = sum(cover(q1r))); taus_fine[!cover(q1)]
round(q1r[q1r$tau %in% taus, c("est", "lo", "hi")], 3)
round(rbind(Q1 = cqr_D(Q1), Q1r = cqr_D(Q1r), truth_Q1 = cqte_true(taus), marginal = qlog_true(taus)), 3)
# OLS on D: -0.284 (SE 0.221) [-0.717, 0.150]; with X and D x X: 0.071 (SE 0.158) [-0.239, 0.381]
# QTE (Q1): -1.747 -1.829 -1.240 -0.451 2.750; SE 0.102 0.150 0.199 0.342 0.501; truth -1.652 -1.606 -1.131 0.210 2.772
#   lower -1.946 -2.124 -1.630 -1.123 1.768; upper -1.548 -1.534 -0.850 0.220 3.733
# intervals covering the true QTE on the 19-level grid: Q1 14, Q1r 19; Q1 misses at tau 0.35 0.55 0.60 0.65 0.70
# QTE (Q1r): -1.782 -1.648 -1.128 -0.129 2.820; lower -1.996 -1.955 -1.555 -0.912 1.515; upper -1.567 -1.340 -0.702 0.654 4.125
# QR of log Y on D, X, D x X, coefficient on D: Q1 -0.284 -0.198 -0.116 -0.041 0.041; Q1r -0.532 -0.325 -0.114 0.096 0.308
# truth for Q1 at X = 0: -0.274 -0.201 -0.120 -0.039 0.034; marginal QTE in logs -0.376 -0.255 -0.120 0.015 0.136
```

At the five levels, the intervals of Q1 all cover the truth. At the 19 plotting levels, Q1 covers at 14, missing at the 35th, 55th, 60th, 65th, and 70th percentiles, while Q1r covers at all 19. This is the result of one draw: the intervals are pointwise, and the sample quantiles at neighboring levels use nearby order statistics, so their errors are highly correlated and misses, once they occur, tend to come in runs; over this stretch, the two groups' sample quantiles each strayed from the truth in opposite directions. Q1r shares its control group with Q1, so the two sets of differences in quantiles differ only through the sampling of the treatment group.

The conditional QR coefficient and the marginal difference in quantiles are not the same thing, even in a randomized experiment: in Q1, on the log scale, the true conditional effect at $\tau = 0.1$ given $X_i = 0$ is -0.274, while the true marginal difference in log quantiles is -0.376. In mean regression, covariates only improve precision; in QR they change the estimand.

### 6.3 Block Q2: RIF Regression, Conditional QR, and Distribution Regression

```r
rif_quantile <- function(y, tau, bw = NULL) {   # RIF(Y; tau) = q + (tau - 1{Y <= q}) / f_Y(q), as rifreg 1.1.0 does
  q <- quantile(y, tau, type = 7, names = FALSE)
  dens <- if (is.null(bw)) density(y, from = q, to = q) else density(y, bw = bw, from = q, to = q)
  f <- dens$y[1]                                                # Gaussian kernel, Silverman bandwidth by default
  structure(q + (tau - as.numeric(y <= q)) / f, quantile = q, density = f, bw = dens$bw)
}
M <- cbind(1, fbm = Q2$fbm, x = Q2$x)
rif_coef <- function(ly, Mq) qr.coef(Mq, vapply(taus, function(t) rif_quantile(ly, t), numeric(length(ly))))[2, ]
rif_est <- rif_coef(Q2$logy, qr(M))                             # OLS of RIF on (1, W, X): coefficient on W
set.seed(13201)
boot <- vapply(1:999, function(b) { i <- sample.int(nrow(Q2), replace = TRUE); rif_coef(Q2$logy[i], qr(M[i, ])) },
               numeric(length(taus)))                           # bootstrap of the whole procedure
rif_se <- apply(boot, 1, sd)
rif_se_known <- sapply(taus, function(t)                        # HC0 standard error treating the RIF as data
  sqrt(sandwich::vcovHC(lm(rif_quantile(Q2$logy, t) ~ fbm + x, data = Q2), type = "HC0")["fbm", "fbm"]))
cqr <- t(sapply(taus, function(t) summary(rq(logy ~ fbm + x, tau = t, data = Q2), se = "nid")$coefficients["fbm", 1:2]))
uqpe_true <- function(t) {                       # (F_0(q) - F_1(q)) / f_Y(q) at the unconditional quantile q
  sd0 <- sqrt(P2$bx^2 + P2$s0^2); sd1 <- sqrt(P2$bx^2 + P2$s1^2)
  Fm <- function(q) P2$p * pnorm(q, P2$mu + P2$g, sd1) + (1 - P2$p) * pnorm(q, P2$mu, sd0)
  q <- uniroot(function(q) Fm(q) - t, P2$mu + c(-10, 10), tol = 1e-13)$root
  fm <- P2$p * dnorm(q, P2$mu + P2$g, sd1) + (1 - P2$p) * dnorm(q, P2$mu, sd0)
  (pnorm(q, P2$mu, sd0) - pnorm(q, P2$mu + P2$g, sd1)) / fm
}
round(rbind(rif = rif_est, se = rif_se, se_rif_known = rif_se_known, uqpe = sapply(taus, uqpe_true),
            cqr = cqr[, 1], cqr_se = cqr[, 2]), 4)
round(rif_est[5] + c(-1, 1) * qnorm(0.975) * rif_se[5], 3); round(cqr[5, 1] + c(-1, 1) * qnorm(0.975) * cqr[5, 2], 3)
th <- quantile(Q2$logy, seq(0.005, 0.995, by = 0.005), type = 1, names = FALSE)   # 199 thresholds
gam <- dist_reg(Q2$logy, cbind(fbm = Q2$fbm, x = Q2$x, fbm_x = Q2$fbm * Q2$x), th, link = "probit")
F1 <- dr_cdf(gam, cbind(fbm = 1, x = Q2$x, fbm_x = Q2$x), "probit")   # every merchant with FBM
F0 <- dr_cdf(gam, cbind(fbm = 0, x = Q2$x, fbm_x = 0), "probit")      # every merchant without
p <- mean(Q2$fbm)
qinv <- function(Fv, t) approx(cummax(Fv), th, xout = t, ties = "ordered")$y
dr <- (sapply(taus, qinv, Fv = (p + 0.1) * F1 + (1 - p - 0.1) * F0) - sapply(taus, qinv, Fv = p * F1 + (1 - p) * F0)) / 0.1
round(dr, 3)
# RIF coefficient on FBM: 0.3666 0.4067 0.4852 0.5555 0.7225; bootstrap SE 0.0251 0.0230 0.0213 0.0299 0.0527
# SE treating the RIF as data: 0.0228 0.0194 0.0181 0.0243 0.0435; true UQPE 0.3457 0.4108 0.4961 0.5883 0.6646
# conditional QR coefficient on FBM: 0.4734 0.4871 0.4919 0.5121 0.5228; nid SE 0.0232 0.0164 0.0144 0.0165 0.0214
# tau = 0.9: RIF [0.619, 0.826]; conditional QR [0.481, 0.565]
# distribution regression, (q_cf - q) / 0.1 for the FBM share raised by 0.1: 0.386 0.471 0.510 0.622 0.736
```

At the five levels, the intervals of RIF regression all cover the UQPE; the intervals of conditional QR all cover the conditional effect 0.5 and, except at the median, none covers the UQPE. The standard error that treats the RIF as data is smaller than the bootstrap one at every level; what it misses is the estimation error of the two steps, the quantile and the density. The `cummax` in distribution regression is merely a safeguard: distribution functions estimated threshold by threshold are not guaranteed to be monotone in finite samples, and taking the cumulative maximum turns them into a valid distribution function.

### 6.4 Block Q3: IV-QR, the LQTE, and the Wald Estimate

```r
ivqr_grid <- function(y, d, z, X, tau, grid, level = 0.95, se = "nid") {   # inverse QR (Chernozhukov-Hansen)
  M <- cbind("(Intercept)" = 1, X, z = z); jz <- ncol(M); crit <- qnorm(1 - (1 - level) / 2)
  res <- vapply(grid, function(b) {                             # QR of Y - D b on (1, X, Z) at each candidate b
    cf <- .muffle_nonunique(summary(.muffle_nonunique(rq(y - d * b ~ M - 1, tau = tau)), se = se))$coefficients
    cf[jz, 1:2]
  }, numeric(2))
  tstat <- res[1, ] / res[2, ]; acc <- which(abs(tstat) <= crit)  # accepted b: |t| of the Z coefficient <= 1.96
  list(est = grid[which.min(abs(tstat))], ci = range(grid[acc]), connected = all(diff(acc) == 1),
       n_accepted = length(acc))
}
lqte_iv <- function(y, d, z, taus) {             # complier cdfs as Wald ratios, then generalized inverses
  n1 <- sum(z == 1); n0 <- sum(z == 0); fs <- mean(d[z == 1]) - mean(d[z == 0])
  supp <- sort(unique(y))
  cnt <- function(sel) if (any(sel)) findInterval(supp, sort(y[sel])) else numeric(length(supp))
  F1 <- (cnt(z == 1 & d == 1) / n1 - cnt(z == 0 & d == 1) / n0) / fs      # Y(1) among compliers
  F0 <- (cnt(z == 0 & d == 0) / n0 - cnt(z == 1 & d == 0) / n1) / fs      # Y(0) among compliers
  data.frame(tau = taus, lqte = gen_inverse(supp, F1, taus) - gen_inverse(supp, F0, taus), first_stage = fs)
}
wald <- function(ly, d, z) (mean(ly[z == 1]) - mean(ly[z == 0])) / (mean(d[z == 1]) - mean(d[z == 0]))
grid <- round(seq(-0.30, 0.60, by = 0.002), 3)
naive <- t(sapply(taus, function(t) summary(.muffle_nonunique(rq(logy ~ D + x, tau = t, data = Q3)),
                                            se = "nid")$coefficients["D", 1:2]))
iv <- lapply(taus, function(t) ivqr_grid(Q3$logy, Q3$D, Q3$Z, cbind(x = Q3$x), t, grid))
lq <- lqte_iv(Q3$logy, Q3$D, Q3$Z, taus)
i1 <- which(Q3$Z == 1); i0 <- which(Q3$Z == 0)
set.seed(13301)
boot <- vapply(1:999, function(b) {                             # bootstrap within each Z arm
  i <- c(i1[sample.int(length(i1), replace = TRUE)], i0[sample.int(length(i0), replace = TRUE)])
  c(lqte_iv(Q3$logy[i], Q3$D[i], Q3$Z[i], taus)$lqte, wald(Q3$logy[i], Q3$D[i], Q3$Z[i]))
}, numeric(length(taus) + 1))
se <- apply(boot, 1, sd)
round(c(first_stage = lq$first_stage[1], naive = naive[3, 1], se = naive[3, 2],
        lo = naive[3, 1] - qnorm(0.975) * naive[3, 2], hi = naive[3, 1] + qnorm(0.975) * naive[3, 2]), 3)
round(t(sapply(iv, function(r) c(r$est, r$ci, r$connected, r$n_accepted))), 3)
round(rbind(lqte = lq$lqte, se = se[1:5], lo = lq$lqte - qnorm(0.975) * se[1:5], hi = lq$lqte + qnorm(0.975) * se[1:5]), 3)
w <- wald(Q3$logy, Q3$D, Q3$Z); round(c(wald = w, se = se[6], lo = w - qnorm(0.975) * se[6], hi = w + qnorm(0.975) * se[6]), 4)
# first stage 0.399; naive QR at the median 0.466 (SE 0.016) [0.434, 0.498]
# IV-QR (estimate, lower, upper, accepted set connected, accepted grid points):
#   tau 0.1  0.002 -0.182 0.124 1 154;  tau 0.25 0.018 -0.070 0.132 0 101;  tau 0.5 0.130 0.048 0.218 1 86
#   tau 0.75 0.202  0.128 0.288 0  80;  tau 0.9  0.202  0.072 0.350 1 140
# LQTE 0.067 0.017 0.099 0.163 0.271; SE 0.067 0.045 0.048 0.057 0.061
#   lower -0.064 -0.071 0.006 0.050 0.151; upper 0.198 0.106 0.193 0.275 0.391
# Wald 0.0996 (SE 0.0418) [0.0177, 0.1815]
```

The five IV-QR confidence sets all cover $\beta_D(\tau)$ at their respective levels, and the five LQTE intervals all cover the compliers' true values; at the 25th and 75th percentiles the accepted grid points are not connected, and each convex hull contains one grid point more than the accepted set.

### 6.5 Block Q4: A Monte Carlo for Panel QR

```r
canay_qr <- function(y, X, id, tau) {             # Canay (2011): alpha_i from the within fit, then pooled QR
  X <- as.matrix(X)
  bmean <- qr.coef(qr(X - apply(X, 2, function(col) ave(col, id))), y - ave(y, id))   # within estimator
  ah <- ave(drop(y - X %*% bmean), id)                          # alpha_hat_i: merchant mean of Y - X'b
  fit <- .muffle_nonunique(rq.fit(cbind(1, X), y - ah, tau = tau, method = "br"))
  list(beta = setNames(fit$coefficients, c("(Intercept)", colnames(X))), beta_mean = bmean, alpha_hat = ah)
}
est_panel <- function(d) {                        # slope at every tau: pooled QR, one dummy per merchant, Canay
  X <- cbind(x = d$x)
  rbind(pooled = sapply(taus, function(t) .muffle_nonunique(rq.fit(cbind(1, d$x), d$y, tau = t, method = "br"))$coefficients[2]),
        dummies = sapply(taus, function(t) fe_dummy_qr(d$y, X, d$id, t)$beta[["x"]]),   # sparse Frisch-Newton
        canay = sapply(taus, function(t) canay_qr(d$y, X, d$id, t)$beta[["x"]]))
}
beta_true <- P4$beta + P4$sig * P4$gam * qnorm(taus)              # 0.4 + 0.3 z_tau
round(rbind(T5 = est_panel(Q4_T5)[, 5], T20 = est_panel(Q4_T20)[, 5]), 3)   # single samples, tau = 0.9
reps <- parallel::mclapply(1:1000, function(r) {
  set.seed(134000 + r)
  c(est_panel(gen_Q4(P4, 5)), est_panel(gen_Q4(P4, 20)))        # estimator fastest, then tau, then panel length
}, mc.cores = 8, mc.set.seed = FALSE)
A <- array(do.call(rbind, reps), dim = c(1000, 3, 5, 2),
           dimnames = list(NULL, c("pooled", "dummies", "canay"), taus, c("T5", "T20")))
mc_mean <- apply(A, c(2, 3, 4), mean)
round(mc_mean[, , "T5"], 4); round(mc_mean[, , "T20"], 4); round(beta_true, 4)
sdA <- apply(A, c(2, 3, 4), sd)[, "0.9", ]; round(sdA, 4); round(sdA / sqrt(1000), 4)   # sd and Monte Carlo SE of the mean
# single samples at tau = 0.9: panel length 5  pooled 1.891 dummies 0.778 Canay 0.875
#                              panel length 20 pooled 1.942 dummies 0.792 Canay 0.851
# Monte Carlo means, panel length 5:  pooled 1.3306 1.4495 1.5743 1.7076 1.8230; dummies 0.1964 0.2422 0.4008 0.5628 0.6091
#                                     Canay 0.1140 0.2499 0.4022 0.5553 0.6903
# Monte Carlo means, panel length 20: pooled 1.3339 1.4474 1.5701 1.6968 1.8102; dummies 0.0416 0.2083 0.3999 0.5919 0.7573
#                                     Canay 0.0425 0.2109 0.3991 0.5873 0.7581
# truth beta(tau) 0.0155 0.1977 0.4000 0.6023 0.7845
# tau = 0.9, sd (length 5, 20): pooled 0.2095 0.1437; dummies 0.1242 0.0641; Canay 0.1154 0.0552
#   Monte Carlo SE of the mean: pooled 0.0066 0.0045; dummies 0.0039 0.0020; Canay 0.0036 0.0017
```

At $\mathcal{T} = 5$ both fixed effects estimators shrink toward the middle at both ends: they are biased upward at the 10th percentile, biased downward at the 90th, and nearly unbiased at the median; at $\mathcal{T} = 20$ both averages are already close to the truth. This agrees with the mechanism of Section 3.8, but it is only a finite-sample illustration for this one design and these two panel lengths.

### 6.6 Blocks Q5 and Q5f: CIC, Quantile DiD, and Switchers

`qq_map` computes $Q_d$ in integer arithmetic: the empirical distribution function takes the value $k_0/n_0$, and the generalized inverse takes, among the values satisfying $j\,n_0 \geq n_1k_0$, the smallest $j$, so that when cell sizes are equal or proportional there is no misalignment from floating-point error.

```r
qq_map <- function(yeval, y_d00, y_d01) {         # Q_d(y) = F^{-1}_{Y_d01}(F_{Y_d00}(y)) in exact integer arithmetic
  s00 <- sort(y_d00); s01 <- sort(y_d01); n0 <- length(s00); n1 <- length(s01)
  k0 <- findInterval(yeval, s00)                                # n0 F_{Y_d00}(y)
  target <- as.double(n1) * k0
  j <- ceiling(target / n0); j <- j + (j * n0 < target) - ((j - 1) * n0 >= target)   # smallest j with j n0 >= n1 k0
  s01[pmax(j, 1)]
}
cic_qtt <- function(y00, y01, y10, taus, y11) {   # QTT(tau) = F^{-1}_{Y11}(tau) - Q_{d=0}(F^{-1}_{Y10}(tau))
  cf <- qq_map(emp_quantile(y10, taus), y00, y01); q11 <- emp_quantile(y11, taus)
  data.frame(tau = taus, qtt = q11 - cf, q11 = q11, cf = cf)
}
qdid <- function(y00, y01, y10, y11, taus) emp_quantile(y11, taus) - emp_quantile(y10, taus) -
  (emp_quantile(y01, taus) - emp_quantile(y00, taus))           # quantile DiD, level by level
cells <- lapply(list(c00 = c(0, 0), c01 = c(0, 1), c10 = c(1, 0), c11 = c(1, 1)),
                function(st) Q5$y[Q5$S == st[1] & Q5$T == st[2]])
both <- function(cl) c(cic_qtt(cl$c00, cl$c01, cl$c10, taus_fine, cl$c11)$qtt, qdid(cl$c00, cl$c01, cl$c10, cl$c11, taus_fine))
est <- both(cells)
set.seed(13501)
boot <- vapply(1:999, function(b) both(lapply(cells, function(v) v[sample.int(length(v), replace = TRUE)])),
               numeric(2 * length(taus_fine)))                  # bootstrap within each region-quarter cell
se <- apply(boot, 1, sd); lo <- est - qnorm(0.975) * se; hi <- est + qnorm(0.975) * se
qtt_true <- function(t) P5$delta * P5$A * exp(P5$g + P5$s * (1 + P5$kap) * (P5$mu1 + P5$sd1 * qnorm(t)))
k <- which(taus_fine %in% taus); kq <- length(taus_fine) + k; tr <- qtt_true(taus_fine)
round(rbind(cic = est[k], lo = lo[k], hi = hi[k], qdid = est[kq], qdid_lo = lo[kq], qdid_hi = hi[kq], truth = tr[k]), 3)
c(cic = sum(lo[1:19] <= tr & tr <= hi[1:19]), qdid = sum(lo[20:38] <= tr & tr <= hi[20:38]))
round(c(est[4], lo[4], hi[4], tr[4]), 3)                        # tau = 0.2
set.seed(13502); lq <- cic_lqte(Q5f, "y", "D", "S", "T", 0.5, B = 999)   # switchers' LQTE at the median
round(c(est = lq$est, lo = lq$ci[1], hi = lq$ci[2], first_stage = lq$first_stage,
        wald = lq$wald, wald_lo = lq$ci_wald[1], wald_hi = lq$ci_wald[2]), 3)
round(sapply(taus, function(t) cic_lqte(Q5f, "y", "D", "S", "T", t, B = 0)$est), 3)   # point estimates only
round(with(Q5f, tapply(D, list(S = S, T = T), mean)), 4)
# CIC  0.477 0.802 1.744 1.852 3.422; lower -0.132 0.134 0.610 0.185 -0.313; upper 1.086 1.471 2.879 3.519 7.157
# QDID 1.037 1.902 3.575 4.898 6.937; lower 0.543 1.352 2.862 3.663 5.004; upper 1.530 2.452 4.289 6.134 8.871
# true QTT 0.937 1.316 1.920 2.801 3.935
# intervals covering the true QTT on the 19-level grid: CIC 18, QDID 5; CIC at tau = 0.2: 0.586 [0.027, 1.146], truth 1.198
# switchers' LQTE at the median 2.158 [1.225, 3.091]; first stage 0.602; LATE (Wald-CIC) 3.016 [2.078, 3.953]
# switchers' LQTE, point estimates: 1.372 2.065 2.158 3.030 6.100
# adoption rates by (S, T): (0, 0) 0.1515  (0, 1) 0.1475  (1, 0) 0.1027  (1, 1) 0.7050
```

The CIC interval for the launch region's treated group at the 90th percentile, [-0.313, 7.157], is very wide: merchants at the upper end are scarce, and there $Q_{d=0}$ has to be pinned down by the few high-revenue merchants at the upper end of the control region. The quantile DiD intervals are no wider than CIC's; its problem is bias, not noise.

### 6.7 The Grand Reconciliation

| Block | Estimator | Estimate | 95% interval | Target and true value | Covered |
|---|---|---|---|---|---|
| Q1 | OLS mean effect | -0.284 | [-0.717, 0.150] | ATE 0 | Yes |
| Q1 | Difference in quantiles, 10th percentile | -1.747 | [-1.946, -1.548] | QTE -1.652 | Yes |
| Q1 | Difference in quantiles, 90th percentile | 2.750 | [1.768, 3.733] | QTE 2.772 | Yes |
| Q1r | Difference in quantiles, 90th percentile | 2.820 | [1.515, 4.125] | QTE 2.772 | Yes |
| Q2 | RIF regression, 90th percentile | 0.723 | [0.619, 0.826] | UQPE 0.665 | Yes |
| Q2 | Conditional QR, 90th percentile | 0.523 | [0.481, 0.565] | UQPE 0.665; conditional effect 0.5 | No; yes |
| Q3 | Naive QR, median | 0.466 | [0.434, 0.498] | $\beta_D$ 0.10 | No |
| Q3 | IV-QR, median | 0.130 | [0.048, 0.218] | $\beta_D$ and LQTE both 0.10 | Yes |
| Q3 | IV-QR, 25th percentile (convex hull) | 0.018 | [-0.070, 0.132] | $\beta_D$ 0.033 | Yes |
| Q3 | LQTE, median | 0.099 | [0.006, 0.193] | LQTE 0.10 | Yes |
| Q3 | Wald estimate | 0.0996 | [0.0177, 0.1815] | LATE 0.10 | Yes |
| Q4 | Dummy-variable QR, $\mathcal{T} = 5$, Monte Carlo mean | 0.6091 | Not applicable | $\beta(0.9)$ 0.7845 | Not applicable |
| Q4 | Canay, $\mathcal{T} = 20$, Monte Carlo mean | 0.7581 | Not applicable | $\beta(0.9)$ 0.7845 | Not applicable |
| Q5 | CIC, median | 1.744 | [0.610, 2.879] | QTT 1.920 | Yes |
| Q5 | Quantile DiD, median | 3.575 | [2.862, 4.289] | QTT 1.920 | No |
| Q5 | CIC, 20th percentile | 0.586 | [0.027, 1.146] | QTT 1.198 | No |
| Q5f | Switchers' LQTE, median | 2.158 | [1.225, 3.091] | 2.047 | Yes |
| Q5f | Switchers' LATE | 3.016 | [2.078, 3.953] | 2.350 | Yes |

The "No" entries in the table are of two kinds. The misses of conditional QR against the UQPE, of naive QR against $\beta_D$, and of quantile DiD against the QTT are by design: each estimator aims at something else, and luck has nothing to do with it. The miss of CIC at the 20th percentile is an ordinary single miss of a pointwise interval across 19 levels, and at the five levels CIC covers every time. The two Q4 rows belong to neither kind and read "Not applicable": what must be compared there is the average location of the estimators, not one sample, so there are no intervals.

The reconciliation of this section comes down to one sentence: every estimator aims at the position of its own distribution, and where things go wrong, either that position was read as some individual's effect, or the position of one distribution was taken for the position of another.

## 7 Failure Modes and Robustness

The most common misreading is to take a shift in a distribution for individual effects. "AutoPrice makes small merchants earn less and large merchants earn more" requires rank invariance, which concerns a counterfactual pairing that can never be seen: Q1 and Q1r have identical QTEs, yet the shares of merchants whose revenue falls are far apart (Section 1). There are three practical responses: report only "how far the $\tau$-th position of the distribution moves"; when conclusions at the individual level are needed, report conclusions that do not depend on the pairing, such as the sign result of Theorem 5 (iii), or bounds on the distribution of individual effects (the range given by all the pairings that the two marginal distributions allow; Heckman, Smith and Clements 1997 discuss this kind of problem starting from the Fréchet-Hoeffding bounds of the joint distribution); with covariates, check whether the QTEs by subgroup tell the same story, bearing in mind that this is only circumstantial evidence.

The second kind of error is mixing up conditional and unconditional quantities. A QR coefficient is a shift in the conditional distribution; "FBM makes the upper end of the platform's revenue distribution rise more" calls for RIF or distribution regression, not a direct reading of the QR coefficient. The same holds in randomized experiments: in mean regression, adding covariates only improves precision, but in QR, adding covariates changes the estimand; in Q1, on the log scale, the conditional and marginal effects at the 10th percentile differ (Section 6.2). If you want the marginal QTE, compare the two groups' marginal distributions, or use the weighting of Firpo (2007); do not substitute the conditional QR coefficient for it.

The third kind concerns inference. The `nid`, `ker`, and `iid` standard errors all assume that the conditional quantile is correctly specified, so use the bootstrap when the specification is in doubt; at tail quantiles the density is low and the standard errors are large, so pay special attention to the width of the intervals at the 90th and 95th percentiles; the choice of bandwidth affects the estimate of the sparsity, and conclusions should not flip with the bandwidth. Pointwise intervals cannot be used to say "the whole curve is significant"; the 5 misses among the 19 plotting levels of Q1 do not mean that the estimator is biased. Panel and clustered data must be resampled by unit.

The fourth kind concerns endogenous treatment. IV-QR and the LQTE aim at different people: the former relies on rank similarity to extend its conclusion to everyone, while the latter concerns only the compliers, and in Q3 the two coincide only at the median. Rank similarity cannot be tested directly, and it fails if the adoption decision depends on the change in rank after adoption (as when merchants try the product first and then decide). With a weak instrument, use the confidence set obtained by inverting over a grid, not a Wald-type point estimate plus or minus standard errors; when the acceptance set has gaps, report the convex hull and say so.

The fifth kind concerns panels. With short panels, fixed effects QR has sizable bias in the tails; the number of periods must be large enough relative to the number of units for the asymptotic theory to apply, and the Monte Carlo of Q4 is only an illustration. Canay's two-step estimator relies on the merchant effects only shifting the location; if the merchant effects also change the dispersion, it misses the target, while dummy-variable QR pays the price of needing a longer panel; the location-scale approach of Machado and Santos Silva (2019), in which location and scale move together, is another option, mentioned here but not developed.

The sixth kind concerns before-and-after policy comparisons. CIC requires the distribution of the unobservables within each group to be the same in both periods, and compositional change within a group will be taken for a time effect; when the support condition fails, the highest or lowest stretch of the launch region has no control merchants to match, and the QTT cannot be point identified over that stretch; with a discrete outcome only bounds are available. The fuzzy version further requires the control region's adoption rate not to change, and the switchers' distribution, obtained by subtraction, is least stable in the tails. With panel data, Callaway and Li (2019) give the QTT under another set of assumptions.

Last come computation and conventions. QR curves estimated separately at different levels can cross, producing the absurd result that the 90th percentile lies below the 75th somewhere; the rearrangement of Chernozhukov, Fernández-Val and Galichon (2010) sorts them back into monotone order. The convention for sample quantiles and the interpolation used by software packages can make the numbers differ in the last few decimal places, so when checking results, align the conventions first.

String these failure modes together, and the robustness of this chapter comes down to one question: which position of which distribution do you want to look at? The conditional distribution or the platform-wide one; all merchants, compliers, switchers, or the treated group; a shift in the distribution, or some individual's gain or loss. Answer this question clearly, and the estimator and the assumptions follow; leave it unclear, and even the most refined standard error merely attaches a correct confidence interval to a wrong question.

## 8 Further Reading

::: {.readings}
Required reading, in suggested reading order:

- Wooldridge (2010, *Econometric Analysis of Cross Section and Panel Data*, 2nd ed.) §12.10. The textbook treatment of quantile estimation: consistency via Theorem 12.2, the asymptotic distribution in equations (12.111) to (12.115), and the panel section in equations (12.122) to (12.126); the source of Theorems 2 and 3 and of Section 3.8 of this chapter.
- Koenker and Hallock (2001, Journal of Economic Perspectives). An introductory survey of QR and of how to read its coefficients.
- Angrist, Chernozhukov and Fernández-Val (2006, Econometrica). What QR estimates under misspecification, and inference on the whole QR process.
- Firpo, Fortin and Lemieux (2009, Econometrica). The RIF and the unconditional quantile partial effect.
- Athey and Imbens (2006, Econometrica). Identification and inference for changes-in-changes, and the conditions under which quantile DiD is valid.

Further reading:

- Koenker and Bassett (1978, Econometrica). The original paper on QR and the check function.
- Heckman, Smith and Clements (1997, Review of Economic Studies) and Bitler, Gelbach and Hoynes (2006, American Economic Review). The difference between differences in quantiles and the distribution of individual effects, and what can still be said without rank invariance.
- Firpo (2007, Econometrica). Efficient semiparametric estimation of the QTE under unconfoundedness.
- Chernozhukov and Hansen (2005, Econometrica; 2006, 2008, Journal of Econometrics) and their survey Chernozhukov and Hansen (2013, Annual Review of Economics). IV-QR, rank similarity, and inference robust to weak instruments.
- Abadie, Angrist and Imbens (2002, Econometrica) and Frölich and Melly (2013, Journal of Business & Economic Statistics). Quantile treatment effects for compliers, $\kappa_i$ weighting, and the unconditional version.
- Foresi and Peracchi (1995, Journal of the American Statistical Association) and Chernozhukov, Fernández-Val and Melly (2013, Econometrica). Distribution regression and inference on counterfactual distributions.
- Chernozhukov, Fernández-Val and Galichon (2010, Econometrica). Rearrangement for quantile curves that cross.
- Hall and Sheather (1988, Journal of the Royal Statistical Society: Series B (Methodological)). The bandwidth for estimating the sparsity.
- Koenker (2004, Journal of Multivariate Analysis), Canay (2011, Econometrics Journal), Kato, Galvao and Montes-Rojas (2012, Journal of Econometrics), and Besstremyannaya and Golovan (2019, Econometrics Journal). Fixed effects in panel QR, their asymptotic theory, and the critique of Canay's rate condition.
- de Chaisemartin and D'Haultfœuille (2018, Review of Economic Studies) and Callaway and Li (2019, Quantitative Economics). The switchers' effects in fuzzy designs, and the QTT under panel DiD.
- Wacker (2023). A complete, self-contained compilation of the properties of generalized inverses; the source of Theorem 1 (i) of this chapter.
- Hansen (2014) §9.5 and §9.6. Another textbook account of the asymptotic distributions of LAD and QR, with proof sketches.
:::

::: {.apa-refs}
- Abadie, A., Angrist, J., & Imbens, G. (2002). Instrumental variables estimates of the effect of subsidized training on the quantiles of trainee earnings. *Econometrica, 70*(1), 91-117. https://doi.org/10.1111/1468-0262.00270
- Angrist, J., Chernozhukov, V., & Fernández-Val, I. (2006). Quantile regression under misspecification, with an application to the U.S. wage structure. *Econometrica, 74*(2), 539-563. https://doi.org/10.1111/j.1468-0262.2006.00671.x
- Athey, S., & Imbens, G. W. (2006). Identification and inference in nonlinear difference-in-differences models. *Econometrica, 74*(2), 431-497. https://doi.org/10.1111/j.1468-0262.2006.00668.x
- Besstremyannaya, G., & Golovan, S. (2019). Reconsideration of a simple approach to quantile regression for panel data. *The Econometrics Journal, 22*(3), 292-308. https://doi.org/10.1093/ectj/utz012
- Bitler, M. P., Gelbach, J. B., & Hoynes, H. W. (2006). What mean impacts miss: Distributional effects of welfare reform experiments. *American Economic Review, 96*(4), 988-1012. https://doi.org/10.1257/aer.96.4.988
- Callaway, B., & Li, T. (2019). Quantile treatment effects in difference in differences models with panel data. *Quantitative Economics, 10*(4), 1579-1618. https://doi.org/10.3982/QE935
- Canay, I. A. (2011). A simple approach to quantile regression for panel data. *The Econometrics Journal, 14*(3), 368-386. https://doi.org/10.1111/j.1368-423X.2011.00349.x
- Chernozhukov, V., Fernández-Val, I., & Galichon, A. (2010). Quantile and probability curves without crossing. *Econometrica, 78*(3), 1093-1125. https://doi.org/10.3982/ECTA7880
- Chernozhukov, V., Fernández-Val, I., & Melly, B. (2013). Inference on counterfactual distributions. *Econometrica, 81*(6), 2205-2268. https://doi.org/10.3982/ECTA10582
- Chernozhukov, V., & Hansen, C. (2005). An IV model of quantile treatment effects. *Econometrica, 73*(1), 245-261. https://doi.org/10.1111/j.1468-0262.2005.00570.x
- Chernozhukov, V., & Hansen, C. (2006). Instrumental quantile regression inference for structural and treatment effect models. *Journal of Econometrics, 132*(2), 491-525. https://doi.org/10.1016/j.jeconom.2005.02.009
- Chernozhukov, V., & Hansen, C. (2008). Instrumental variable quantile regression: A robust inference approach. *Journal of Econometrics, 142*(1), 379-398. https://doi.org/10.1016/j.jeconom.2007.06.005
- Chernozhukov, V., & Hansen, C. (2013). Quantile models with endogeneity. *Annual Review of Economics, 5*(1), 57-81. https://doi.org/10.1146/annurev-economics-080511-110952
- de Chaisemartin, C., & D'Haultfœuille, X. (2018). Fuzzy differences-in-differences. *The Review of Economic Studies, 85*(2), 999-1028. https://doi.org/10.1093/restud/rdx049
- Firpo, S. (2007). Efficient semiparametric estimation of quantile treatment effects. *Econometrica, 75*(1), 259-276. https://doi.org/10.1111/j.1468-0262.2007.00738.x
- Firpo, S., Fortin, N. M., & Lemieux, T. (2009). Unconditional quantile regressions. *Econometrica, 77*(3), 953-973. https://doi.org/10.3982/ECTA6822
- Foresi, S., & Peracchi, F. (1995). The conditional distribution of excess returns: An empirical analysis. *Journal of the American Statistical Association, 90*(430), 451-466. https://doi.org/10.1080/01621459.1995.10476537
- Frölich, M., & Melly, B. (2013). Unconditional quantile treatment effects under endogeneity. *Journal of Business & Economic Statistics, 31*(3), 346-357. https://doi.org/10.1080/07350015.2013.803869
- Hall, P., & Sheather, S. J. (1988). On the distribution of a studentized quantile. *Journal of the Royal Statistical Society: Series B (Methodological), 50*(3), 381-391. https://doi.org/10.1111/j.2517-6161.1988.tb01735.x
- Hansen, B. E. (2014). *Econometrics* [Unpublished manuscript]. Department of Economics, University of Wisconsin.
- Heckman, J. J., Smith, J., & Clements, N. (1997). Making the most out of programme evaluations and social experiments: Accounting for heterogeneity in programme impacts. *The Review of Economic Studies, 64*(4), 487-535. https://doi.org/10.2307/2971729
- Huber, P. J. (1967). The behavior of maximum likelihood estimates under nonstandard conditions. In L. M. Le Cam & J. Neyman (Eds.), *Proceedings of the Fifth Berkeley Symposium on Mathematical Statistics and Probability* (Vol. 1, pp. 221-233). University of California Press.
- Kato, K., Galvao, A. F., Jr., & Montes-Rojas, G. V. (2012). Asymptotics for panel quantile regression models with individual effects. *Journal of Econometrics, 170*(1), 76-91. https://doi.org/10.1016/j.jeconom.2012.02.007
- Koenker, R. (2004). Quantile regression for longitudinal data. *Journal of Multivariate Analysis, 91*(1), 74-89. https://doi.org/10.1016/j.jmva.2004.05.006
- Koenker, R., & Bassett, G., Jr. (1978). Regression quantiles. *Econometrica, 46*(1), 33-50. https://doi.org/10.2307/1913643
- Koenker, R., & Hallock, K. F. (2001). Quantile regression. *Journal of Economic Perspectives, 15*(4), 143-156. https://doi.org/10.1257/jep.15.4.143
- Machado, J. A. F., & Santos Silva, J. M. C. (2019). Quantiles via moments. *Journal of Econometrics, 213*(1), 145-173. https://doi.org/10.1016/j.jeconom.2019.04.009
- Meyer, B. D., Viscusi, W. K., & Durbin, D. L. (1995). Workers' compensation and injury duration: Evidence from a natural experiment. *American Economic Review, 85*(3), 322-340.
- Newey, W. K., & McFadden, D. (1994). Large sample estimation and hypothesis testing. In R. F. Engle & D. L. McFadden (Eds.), *Handbook of econometrics* (Vol. 4, pp. 2111-2245). Elsevier. https://doi.org/10.1016/S1573-4412(05)80005-4
- Wacker, P. (2023). *Please, not another note about generalized inverses* [Preprint]. arXiv. https://arxiv.org/abs/2306.06989
- Wooldridge, J. M. (2010). *Econometric analysis of cross section and panel data* (2nd ed.). MIT Press.
:::
