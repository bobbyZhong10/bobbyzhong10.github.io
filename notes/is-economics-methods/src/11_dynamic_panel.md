---
title: "Dynamic Panels"
subtitle: "From Fixed Effects to Arellano-Bond"
seriesline: "Foundations of Information Systems Economics · Chapter 11"
series: "Foundations of Information Systems Economics: Empirical Methods for Information Systems Economics"
sidebartitle: "Foundations of Information Systems Economics"
sidebarsub: "Theory · Methods · Applications"
chapterlabel: "Chapter 11 · Dynamic Panels"
---

## Introduction

Nimbus is the fictional online marketplace of Chapter 10. The question there was static: when the same merchant spends USD 100 more on promotion, how much higher is its GMV that month? This chapter's question adds a time dimension: how long does a merchant's good month last? The platform needs the answer to decide whether a support subsidy should be paid once or several months in a row, and to judge how much of the GMV from a big sales event will stick. Regress this month's GMV directly on last month's, and persistence comes out at 0.946; add merchant fixed effects, and the estimate crashes to 0.425. Both models are perfectly standard, yet the two answers differ by more than a factor of two, while the true value known from the simulation, 0.70, sits squarely between them. The first model mistakes merchant talent for persistence, and the second, in scrubbing that talent away, manufactures a new correlation of its own.

This chapter covers the problems that arise when a lagged dependent variable, or any other regressor that responds to past shocks, enters a panel model. The FE, FD, and strict exogeneity of Chapter 10 are the starting point: a lagged dependent variable violates strict exogeneity by its very nature, so in short panels FE has a bias that does not vanish as the number of merchants grows, known as Nickell bias. The main remedy differences the merchant effect away first and then uses earlier levels as instruments; these are the estimators of Anderson-Hsiao and Arellano-Bond. When persistence approaches 1 these instruments become weak, and Blundell and Bond add a levels equation, at the cost of an assumption about how merchants enter the sample. Dynamic panel GMM is not a push-button remedy either: the instruments can be too weak, or so numerous that they pull the estimate back toward the biased estimator it was meant to avoid, while leaving the overidentification test toothless. Serial correlation tests and overidentification tests are therefore part of the analysis. Panel data really do hand the researcher a second dimension, but time never gives away identification for free.

The chapter proceeds as follows. Section 1 looks at the two wrong numbers that straddle the truth, and Section 2 writes down the estimand of the dynamic model and Nimbus's data blocks. Section 3, the most detailed in the chapter, takes up in turn sequential exogeneity and the lagged dependent variable, Nickell bias, Anderson-Hsiao and Arellano-Bond, inference and specification tests for two-step GMM, system GMM (a two-equation system, differences and levels, in the sense of Section 3.7 of Chapter 9) and the initial conditions, instrument proliferation, predetermined regressors, and bias-corrected FE. Section 4 covers the choice of estimator, a reporting checklist, and software; Section 5 reads three anchoring papers; Section 6 runs the whole toolkit on Nimbus; and Section 7 discusses failure modes. The chapter uses the sandwich of Section 3.4 of Chapter 5, the GMM theory of Chapter 7 (two-step GMM, the Hansen J test, the Windmeijer correction, and many-moment bias), the cluster-robust variance of Section 3.4 of Chapter 8, and the FE, FD, and strict exogeneity tests of Chapter 10.

## 1 Two Wrong Numbers That Straddle the Truth

::: {.case}
Block D1 is Nimbus's inherited dynamic panel: 500 merchants, 8 months each, with GMV measured by a standardized log index. In the simulation, the GMV of merchant $i$ in month $t$ satisfies $Y_{it} = 0.70\,Y_{i,t-1} + \alpha_i + \varepsilon_{it}$: seventy percent of this month's excess over the long-run level carries over into next month; $\alpha_i$ is the merchant's operating talent, and high-talent merchants also have high long-run levels. Management wants to know the persistence $\rho$: the closer it is to 1, the more durable the effect of a shock, and the more worthwhile it becomes to take a long view of the platform's subsidies and campaigns.
:::

The analyst first runs the most direct regression, a pooled OLS of this month's GMV on last month's. The estimated persistence is 0.946, with a merchant-clustered 95% interval of [0.938, 0.955]. This implies that GMV is as sticky as a random walk, with last month's good showing carrying over almost entirely into this month. She then realizes that top merchants post high GMV month after month, and that merchant-specific talent may have inflated persistence. So she adds merchant fixed effects and estimates again, and the estimate is nearly cut in half: 0.425, with an interval of [0.388, 0.461]. One number says GMV is extremely sticky, the other that shocks fade fairly fast. Because this is a simulation, we know the true value is 0.70: both are wrong, and the truth is caught between them.

Neither deviation is sampling error. As the number of merchants goes to infinity, pooled OLS converges to 0.955 and FE to 0.416; each sample estimate is close to its own limit and far from the truth. The biases point in opposite directions, and each has its own source. Pooled OLS runs high because talent lifts both this month's GMV and its lag, and OLS misreads the correlation across merchants created by talent as persistence. FE runs low for a far subtler reason: the demeaning step that removes talent creates a negative correlation between the transformed lag and the transformed error. This is Nickell bias, and Section 3.2 will compute its exact size.

::: {.intuition}
An honest estimate of persistence must do two things at once: remove merchant talent, and keep the step that removes it from damaging the dynamic structure. Pooled OLS fails at the first, FE at the second. It is no luck that the truth, 0.70, is straddled by 0.946 and 0.425; it is the inevitable consequence of two biases pointing in opposite directions. The bracket is useful in itself: an estimator that claims to be consistent but lands outside the two has most likely gone wrong somewhere. Later in the chapter we will see that Arellano-Bond difference GMM uses the model's own past values as instruments, sidesteps both biases at once, and pulls persistence back to around 0.70.
:::

This chapter answers three questions: why FE is biased, by how much, and why the bias does not vanish with the number of merchants (Sections 3.1 and 3.2); whether differencing plus instruments estimates the same $\rho$, and at what cost (Sections 3.3 to 3.5); and when this approach quietly breaks down (Sections 3.5 to 3.8).

## 2 The Economic Model and the Estimand

### 2.1 The Dynamic Model and State Dependence

The basic dynamic panel model puts a lagged dependent variable into the unobserved effects model of Chapter 10:

$$Y_{it} = \rho\, Y_{i,t-1} + X_{it}'\beta + \alpha_i + \varepsilon_{it}, \qquad t = 1, \ldots, T.$$

The unit $i = 1, \ldots, N$ is a merchant, and the notation follows Chapter 10: $\alpha_i$ is the merchant effect, $\varepsilon_{it}$ the idiosyncratic error, $X_{it}$ the other time-varying regressors, and $\Delta$ the first difference. The timing convention is that of Wooldridge (2010) §11.6.2 and Nickell (1981): each merchant contributes $Y_{i0}, Y_{i1}, \ldots, Y_{iT}$, $T + 1$ months in all; the first month supplies only a lag, and the number of regression periods is $T$. Block D1 has 8 months per merchant, so $T = 7$. The convention looks trivial, but Section 3.2 will show that if $T$ is taken to be the number of months, the same bias formula gives a different number.

The special case without $X_{it}$ is the AR(1) panel $Y_{it} = \rho\,Y_{i,t-1} + \alpha_i + \varepsilon_{it}$, and blocks D1, D1p, D2, and D4 all take this form. When $|\rho| < 1$, a merchant's GMV fluctuates around its own long-run level $\alpha_i/(1-\rho)$: talent decides where the merchant settles, and $\rho$ decides how fast it comes back after a deviation. The estimand of this chapter is $\rho$, and, when $X_{it}$ is present, $\beta$ as well; the latter is the same-month effect of promotion, and the long-run effect is $\beta/(1-\rho)$.

Reading these two parameters as causal quantities rests on what Wooldridge (2010) equation (11.61) calls dynamic completeness: given the merchant's talent and last month's GMV, earlier history no longer helps predict this month. $\rho$ answers a question of state dependence: holding talent fixed, does selling well last month in itself make this month's sales better? This is what pooled OLS confounds. The correlation of GMV across months has two sources, the time-invariant heterogeneity of talent and genuine state dependence; pooled OLS adds the two together, and FE tries to strip off the first but damages the second in the process. So management's question, "how long does a good month last?", is a question about state dependence, not about the raw autocorrelation of GMV: a one-off platform subsidy raises this month's GMV, and $k$ months later a share $\rho^k$ of the increase remains; a high-talent merchant sells well every month, and that has nothing to do with the subsidy.

::: {.intuition}
Picture two merchants, one with high talent and a high long-run GMV, the other with low talent and a low long-run level. Even if $\rho = 0$, so that last month's fortunes have no effect at all on this month, the high-talent merchant is high last month and high again this month; pool all merchant-months, and this month's GMV is still highly correlated with last month's. This correlation comes from "who is who," not from "what happened last month." State dependence asks about the latter: if a merchant sold unexpectedly well last month, above its own long-run level, will it also be above its long-run level this month? To answer that, each merchant's own long-run level has to be removed, and how it is removed is where this chapter's trouble lies.
:::

### 2.2 Nimbus's Data Blocks

| Block | Merchants and months | Design | Target and true value | Used in |
|---|---|---|---|---|
| D1 | 500, 8 | AR(1) plus merchant effect, stationary start (inherited) | $\rho = 0.70$ | Section 1, Sections 3.1 to 3.4 and 3.8 |
| D1p | 500, 8 | As D1, with higher persistence (inherited) | $\rho = 0.90$ | Sections 3.2 to 3.5 |
| D1T | 300; 4, 8, 16, 32 | The design of D1, 40 samples per length (inherited) | $\rho = 0.70$ | Section 3.2 |
| D2 | 70, 16 | 70 key accounts, monthly GMV in thousands of USD | $\rho = 0.70$ | Sections 3.3 and 3.6 |
| Twin of D2 | 70, 16 | As D2, with MA(1) errors, 200 samples | Instruments invalid | Section 3.6 |
| D3 | 4,000, 8 | Lagged GMV plus promotion that responds to last month's shock | $\rho = 0.6$, $\beta = 0.5$ | Sections 3.1 and 3.7 |
| D4 | 2,000, 8 | A cohort of merchants onboarded in the same month, starting at 80% of their long-run level | $\rho = 0.6$ | Sections 3.4 and 3.5 |
| D5 | 500; 8, 16, 32 | The design of D1, 1,000 samples per length | $\rho = 0.70$ | Section 3.8 |

The inherited blocks keep standardized units; in the new blocks $Y_{it}$ is monthly GMV in thousands of USD and $X_{it}$ is monthly promotion spending in hundreds of USD. Each new block changes one thing only, putting one assumption in the spotlight. D2 consists of Nimbus's key accounts: long-run monthly GMV averages USD 40,000, with a standard deviation of USD 5,000 across merchants and a standard deviation of USD 2,000 for the monthly shocks; few merchants and many months is the situation in which the number of instruments easily exceeds the number of merchants. Its twin uses the same design except that a shock carries over into the following month, $\varepsilon_{it} = e_{it} + 0.6\,e_{i,t-1}$, where $e_{it}$ is the innovation, so the level two periods back is no longer a valid instrument. A D3 merchant whose GMV comes in USD 1,000 below expectations in some month spends USD 60 more on promotion the next month; each extra USD 100 of promotion brings USD 500 more GMV that month and USD 1,250 more in the long run. D4 consists of 2,000 merchants onboarded in the same month, who reach only 80% of their long-run level in their first month (long-run monthly GMV averages USD 15,000, with a standard deviation of USD 3,000); the correlation between the first month's deviation from the long-run level and talent is $-0.433$, so the more talented the merchant, the farther it starts from its own target. D4 also has a twin sample: the same random draws, except that every merchant starts at its long-run level.

## 3 Identification and Large-Sample Theory

This section first uses Nimbus's numbers to show what is to be proved, then states numbered assumptions and theorems, followed by proof sketches of two to six steps, each step naming the assumptions it uses. Finally, an "In words" paragraph answers four questions: under which assumptions it estimates which quantity, what it estimates when the specification is wrong, who sets the weights, and what Nimbus's numbers say. Wooldridge (2010) §11.6 has numbered equations but no numbered theorems. For Nickell bias, it gives only the order $O(1/T)$ and the sign (p. 323), not the exact formula; the Windmeijer correction and bias correction are not in the book, and the initial conditions of system GMM receive a single sentence (p. 374). The theorem boxes therefore cite equation numbers, pages, and the original papers, and the theorem numbering is this chapter's own. Except in Section 3.8, the large-sample theory is always $N \to \infty$ with $T$ fixed.

### 3.1 Why Strict Exogeneity Must Fail: Sequential Exogeneity and the Lagged Dependent Variable

FE.1 of Section 3.1 of Chapter 10 requires the current error to be unrelated to the regressors of every month, past, present, and future. A lagged dependent variable cannot meet this: $Y_{it}$ is the regressor of month $t + 1$, and by definition it contains $\varepsilon_{it}$. In block D1 the variance of $\varepsilon_{it}$ is 1, so $\mathbb{E}[Y_{it}\,\varepsilon_{it}] = 1$: next month's regressor fully "knows" this month's shock. Block P3 of Chapter 10 showed the same kind of failure: merchants whose sales were poor last month promote more this month, and the lead test at the end of Section 3.4 of Chapter 10 rejected strict exogeneity with a $p$-value of practically zero. A lagged dependent variable is the extreme form of this feedback: it is not just a variable that may respond to past shocks; it is built directly out of them. We need an exogeneity condition that is weaker and better suited to the dynamic model.

::: {.assumption}
**SX.1 (sequential exogeneity; Wooldridge 2010 equation (11.51), p. 368)** $\mathbb{E}[\varepsilon_{it} \mid X_{it}, X_{i,t-1}, \ldots, X_{i1}, \alpha_i] = 0$, $t = 1, \ldots, T$.

**DC.1 (dynamic completeness; Wooldridge 2010 equations (11.60) and (11.61), p. 371)** In the AR(1) model $Y_{it} = \rho\,Y_{i,t-1} + \alpha_i + \varepsilon_{it}$, $\mathbb{E}[\varepsilon_{it} \mid Y_{i,t-1}, Y_{i,t-2}, \ldots, Y_{i0}, \alpha_i] = 0$, $t = 1, \ldots, T$.
:::

SX.1 requires only that the current error be unrelated to the regressors of the current and earlier periods, and allows future regressors to respond to today's shock. Equivalently (equation (11.52)), given the current regressors and the merchant effect, past regressors no longer affect this month's conditional mean. The conditioning set stops at period $t$, which is far more natural than the "all periods" of strict exogeneity: how a merchant spends next month may well depend on how this month's sales went. DC.1 is the special case of SX.1 with $X_{it} = Y_{i,t-1}$, and it says that the first lag already captures all of the dynamics. It has an often overlooked implication: the error has no serial correlation, because $\varepsilon_{i,t-1} = Y_{i,t-1} - \rho\,Y_{i,t-2} - \alpha_i$ lies in the conditioning set. If GMV really needs two lags, or if shocks themselves persist (one stockout affecting two months), DC.1 fails, and the $m_2$ test of Section 3.4 watches for this.

::: {.theorem}
**Theorem 1 (moment conditions under sequential exogeneity; Wooldridge 2010 §11.6.1 and §11.6.2, equations (11.56), (11.57), and (11.59), pp. 370-372)** (i) Under SX.1, $\mathbb{E}[X_{ij}\,\varepsilon_{it}] = 0$ for $j \le t$; the first-differenced equation $\Delta Y_{it} = \Delta X_{it}'\beta + \Delta\varepsilon_{it}$ satisfies

$$\begin{aligned}
&\mathbb{E}[X_{ij}\,\Delta\varepsilon_{it}] = 0, \\
&j = 1, \ldots, t - 1, \quad t = 2, \ldots, T,
\end{aligned}$$

so that the number of available instruments grows with $t$. (ii) Under DC.1, the AR(1) model is the special case of (i) with $X_{it} = Y_{i,t-1}$, and the moment conditions are $\mathbb{E}[Y_{i,t-s}\,\Delta\varepsilon_{it}] = 0$, $s \ge 2$; moreover, $\mathbb{E}[\varepsilon_{it}\,\varepsilon_{i,t-1}] = 0$. (iii) But FE.1 of Chapter 10 fails: $\mathbb{E}[Y_{it}\,\varepsilon_{it}] = \mathbb{E}[\varepsilon_{it}^2] > 0$.
:::

The proof sketch has four steps. Step one (SX.1): for $j \le t$, $X_{ij}$ is in the conditioning set, and iterated expectations give $\mathbb{E}[X_{ij}\varepsilon_{it}] = \mathbb{E}\big[X_{ij}\,\mathbb{E}(\varepsilon_{it} \mid X_{it}, \ldots, X_{i1}, \alpha_i)\big] = 0$. Step two (step one): $\Delta\varepsilon_{it} = \varepsilon_{it} - \varepsilon_{i,t-1}$, and for $j \le t - 1$ each of the two terms is orthogonal to $X_{ij}$, which is equation (11.57); $j = t$ is excluded, because SX.1 allows $X_{it}$ to be correlated with $\varepsilon_{i,t-1}$. Step three (DC.1): take $X_{it} = Y_{i,t-1}$; the conditioning set of DC.1 is exactly that of SX.1, $X_{ij} = Y_{i,j-1}$ with $j \le t-1$, and in terms of the lag order $s = t - j + 1$ this is $Y_{i,t-s}$ with $s \ge 2$; $\varepsilon_{i,t-1}$ is a function of the conditioning set and is therefore uncorrelated with $\varepsilon_{it}$. Step four (DC.1): the model itself gives $\mathbb{E}[Y_{it}\varepsilon_{it}] = \rho\,\mathbb{E}[Y_{i,t-1}\varepsilon_{it}] + \mathbb{E}[\alpha_i\varepsilon_{it}] + \mathbb{E}[\varepsilon_{it}^2]$; the first two terms are zero by DC.1, and since $Y_{it}$ is the regressor of period $t+1$, FE.1 fails.

In words: under SX.1, the differenced equation with lagged regressors as instruments estimates the $\beta$ of the conditional mean under sequential exogeneity (equation (11.52)); under DC.1 it estimates the state dependence $\rho$, not the raw autocorrelation of GMV. If the dynamics are misspecified, say the error is in fact serially correlated, then at least the moment condition of (ii) with $s = 2$ fails (under MA(1) errors those with $s \ge 3$ still hold; under autoregressive errors all of them fail), $Y_{i,t-2}$ is correlated with $\Delta\varepsilon_{it}$ through a common shock, and the estimator converges to a quantity contaminated by the serial correlation; the twin block of Section 3.6 is this case. The question of the weights has only half an answer here: the moment conditions say which variables may serve as instruments, with more instruments available in later months; how the moment conditions are weighted against one another is decided by the weight matrix of Section 3.3. Back to Nimbus: in block D1, $\mathbb{E}[Y_{it}\,\varepsilon_{it}] = 1$, so FE, which relies on FE.1, has no reason to be consistent here, and it gives 0.425; in block D3 the coefficient of promotion's response to last month's shock is $-0.6$, promotion likewise satisfies only SX.1, and FE estimates the same-month effect of promotion at 0.346, against a true value of 0.5.

### 3.2 Nickell Bias: How Demeaning Mixes the Error into the Lag

FE gives 0.425 on D1. Repeat the same design 40 times on panels of 4, 8, 16, and 32 months (block D1T, 300 merchants each time), and the average FE estimates are 0.082, 0.417, 0.570, and 0.640: the shorter the panel, the worse the bias, and at 32 months it still falls some way short of the truth. The exact size of this bias was worked out by Nickell (1981).

::: {.intuition}
The source of the bias is a hidden consequence of the demeaning step. Within demeaning of the dynamic model turns the lag into $Y_{i,t-1} - \bar Y_{i,-1}$, where $\bar Y_{i,-1} = T^{-1}\sum_{t=1}^T Y_{i,t-1}$ is the mean of this merchant's lags over the periods, and the error into $\varepsilon_{it} - \bar\varepsilon_i$. For FE to be unbiased, the two must be uncorrelated; but they are correlated, and the culprit is the two means. $\bar Y_{i,-1}$ contains GMV from month 1 to month $T-1$, and each month's GMV contains that month's shock, so $\bar Y_{i,-1}$ picks up a component of $\bar\varepsilon_i$, and the two means are positively correlated. Through the $-\bar Y_{i,-1}$ term, the demeaned lag is negatively correlated with the demeaned error, and a regressor negatively correlated with the error drags the coefficient down: this is why the bias is downward. The magnitude of the contamination is $1/T$: the mean is an average of $T$ numbers, a single shock makes up only a $1/T$ share of it, and the more months there are, the more it is diluted; more merchants do not help, because the contamination happens inside each merchant's own mean.
:::

::: {.assumption}
**NK.1 (the setting of Nickell 1981)** (a) DC.1 holds; (b) $\varepsilon_{it}$ is independent and identically distributed across $i$ and $t$, with mean zero and variance $\sigma_\varepsilon^2$, and independent of $\alpha_i$; (c) $|\rho| < 1$, and the initial value is drawn from the stationary distribution: $Y_{i0} - \alpha_i/(1-\rho)$ is independent of $\alpha_i$ and of later shocks, with variance $\sigma_\varepsilon^2/(1-\rho^2)$.
:::

::: {.theorem}
**Theorem 2 (Nickell bias; Nickell 1981 equation (18); for the general $O(1/T)$ result, see Wooldridge 2010 equation (10.76), p. 323)** Under NK.1, with $\zeta_T = (1 - \rho^T)/(1 - \rho)$, the FE (within) estimator satisfies

$$\begin{aligned}
&\mathrm{plim}_{N \to \infty}(\hat\rho_{FE} - \rho) \\
&\quad = -\frac{1+\rho}{T-1}\Big(1 - \frac{\zeta_T}{T}\Big) \\
&\quad \times \Big[1 - \frac{2\rho}{(1-\rho)(T-1)}\Big(1 - \frac{\zeta_T}{T}\Big)\Big]^{-1}.
\end{aligned}$$

The bias is negative and does not vanish with $N$; as $T \to \infty$, $T \cdot \mathrm{plim}_{N \to \infty}(\hat\rho_{FE} - \rho) \to -(1+\rho)$.
:::

The proof sketch has four steps. Step one (DC.1): FE is the pooled OLS of the demeaned $Y_{it}$ on the demeaned $Y_{i,t-1}$ (Theorem 3 of Chapter 10), and the numerator of $\hat\rho_{FE} - \rho$ is $N^{-1}\sum_i\big(\sum_t Y_{i,t-1}\varepsilon_{it} - T\,\bar Y_{i,-1}\bar\varepsilon_i\big)$; each $Y_{i,t-1}$ depends only on earlier shocks, so the first term has expectation zero. Step two (NK.1(b) and (c)): the second term is not zero. For $j \le t - 1$, $\mathrm{Cov}(Y_{i,t-1}, \varepsilon_{ij}) = \rho^{t-1-j}\sigma_\varepsilon^2$, and otherwise it is zero; summing over $t$ and $j$ gives $\mathbb{E}[\bar Y_{i,-1}\bar\varepsilon_i] = \sigma_\varepsilon^2\big(1 - \zeta_T/T\big)/\big(T(1-\rho)\big) > 0$, so the numerator converges to $-\sigma_\varepsilon^2(1 - \zeta_T/T)/(1-\rho)$. Step three (NK.1(c) and the law of large numbers): the denominator is the sum of squares of the demeaned lag, which converges to a positive number that depends only on $\rho$, $T$, and $\sigma_\varepsilon^2$; in the ratio $\sigma_\varepsilon^2$ cancels, leaving the expression in the theorem. Step four (expansion): as $T \to \infty$, $\zeta_T \to 1/(1-\rho)$, both brackets are $1 + O(1/T)$, and the bias equals $-(1+\rho)/T + O(T^{-2})$.

The first factor of the exact expression is often singled out as the leading term of Nickell bias:

$$\mathrm{plim}_{N \to \infty}\,(\hat\rho_{FE} - \rho) \approx -\,\frac{1 + \rho}{T - 1}.$$

Strictly speaking it is only one factor, and it agrees with the exact bias only to order $1/T$; expanded in $1/T$, the first-order term is $-(1+\rho)/T$. The $T$ in the formula is the number of regression periods. In block D1, $\rho = 0.70$ and $T = 7$; the exact limit is 0.416, and the first factor alone gives 0.417, almost the same; the first-order term, $\rho - (1+\rho)/T$, gives 0.457, clearly too high. The near coincidence of the first factor is an accident of these parameter values: in D1p, with $\rho = 0.90$, the first factor gives 0.583 against an exact limit of 0.558, and the two have already come apart. If $T$ is taken to be the 8 months and plugged into the first factor, the resulting $-(1+\rho)/7$ is exactly the first-order term, and it too lands at 0.457. In short panels the three are further apart: at 4 months ($T = 3$) the exact limit is 0.080, the first factor gives $-0.150$, and the first-order term gives 0.133. So use the exact expression in short panels, and count $T$ carefully.

In words: under NK.1, FE estimates not $\rho$ but $\rho$ plus a negative bias determined by $\rho$ and $T$ alone; more merchants only make the estimate settle more precisely in the wrong place, and a platform panel with many merchants and few months is the worst case. What happens when the specification is wrong: the exact expression relies on a stationary start and iid errors, and with a nonstationary start or heteroskedastic errors the size of the bias changes, but the mechanism, "demeaning mixes the current shock into the lag," does not, and FE is still generally inconsistent; beyond contemporaneous exogeneity (equation (10.74), which holds automatically under DC.1), the $O(1/T)$ result of Wooldridge (2010) equation (10.76) requires only that the process be stable and weakly dependent. The weights are set by within-unit variation: each merchant's contribution is weighted by its own sum of squared demeaned lags, and the bias comes from the mean over $T$ months inside each unit, so it shrinks at rate $1/T$, not $1/N$. Back to Nimbus: D1's sample estimate is 0.425, and its merchant-clustered interval [0.388, 0.461] covers the limit 0.416 but not 0.70; in D1p, with $\rho = 0.90$, the limit is 0.558 and the sample estimate 0.535, so the larger $\rho$, the larger the bias; D1T's Monte Carlo means differ from the exact limits 0.080, 0.416, 0.572, and 0.641 by at most 0.0024 (figure below).

![FE estimates of the persistence $\rho$ (true value 0.70, horizontal dashed line) against the number of months each merchant is observed; the horizontal axis is the number of months, $T + 1$, and the number of regression periods $T$ is one month fewer. The three curves are the exact expression for the probability limit of FE and two approximations: the solid line is the exact Nickell formula, the long-dashed line keeps only the first factor, $\rho - (1+\rho)/(T-1)$, and the dotted line is $\rho - (1+\rho)/T$, that is, the bias of order $1/T$ that the Hahn-Kuersteiner correction removes; the vertical axis is cut off a little below zero. Dots are the Monte Carlo means of FE at 4, 8, 16, and 32 months (40 samples per length, 300 merchants each), and triangles are the means of one-step Arellano-Bond on the same samples (each lag order merged into one column of instruments, lags 2 to 4); squares and diamonds are the Monte Carlo means of the two bias corrections of Section 3.8 at 8, 16, and 32 months (1,000 samples per length, 500 merchants each).](assets/fig/fig_11_nickell.svg)

::: {.warning}
The most dangerous thing about Nickell bias is that it is silent. The FE regression runs just fine, coefficients and standard errors all present, with no error message of any kind. An unwitting analyst would take the 0.425 and report "GMV persistence is moderate," when the truth is 0.70, a fairly high persistence, and the two conclusions give opposite advice on whether subsidies should be a long-term commitment. More insidiously, the more merchants there are, the more precise that wrong 0.425 looks: the standard errors grow ever smaller, dressing up a biased estimate as ever more credible.
:::

Compare pooled OLS. It leaves talent untouched, so the lag enters the regression together with the $\alpha_i$ in the error, and its probability limit can be written down directly:

$$\begin{aligned}
&\mathrm{plim}_{N \to \infty}\,\hat\rho_{POLS} \\
&\quad = \rho + (1-\rho)\,\frac{\mathrm{Var}(\alpha_i)}{\mathrm{Var}(\alpha_i) + \sigma_\varepsilon^2(1-\rho)/(1+\rho)}.
\end{aligned}$$

The fraction is the share of the variance of the long-run level $\alpha_i/(1-\rho)$ in the total variance of GMV (numerator and denominator have both been multiplied by $(1-\rho)^2$): the larger the differences in talent, the closer pooled OLS is to 1. In D1, $\mathrm{Var}(\alpha_i) = \sigma_\varepsilon^2 = 1$, and the variance of the long-run level is far larger than that of the fluctuations around it; the limit is 0.955, with the sample estimate 0.946 right beside it, and the limit for D1p is 0.995.

### 3.3 Differencing Plus Lagged Levels: Anderson-Hsiao and Arellano-Bond

The first step of the remedy is to switch to a different way of removing $\alpha_i$, using first differencing in place of demeaning:

$$\Delta Y_{it} = \rho\, \Delta Y_{i,t-1} + \Delta\varepsilon_{it}, \qquad t = 2, \ldots, T.$$

Differencing likewise removes $\alpha_i$, but it does not solve the endogeneity: $\Delta Y_{i,t-1}$ contains $\varepsilon_{i,t-1}$, and $\Delta\varepsilon_{it}$ also contains $\varepsilon_{i,t-1}$; Hansen (2014) §19.3 computes their covariance as $-\sigma_\varepsilon^2$. OLS run directly on the differenced equation is the FD estimator of Chapter 10. On the sample of block D2 from Section 3.6, it gives $-0.165$, against a limit of $-(1-\rho)/2 = -0.150$: even the sign is wrong.

The key move is to find an instrument for $\Delta Y_{i,t-1}$, and the answer hides in the model's own history. $Y_{i,t-2}$ is clearly correlated with $\Delta Y_{i,t-1} = Y_{i,t-1} - Y_{i,t-2}$; it is determined by information at $t-2$ and earlier, and under DC.1 it is uncorrelated with $\Delta\varepsilon_{it}$. This gives the moment condition

$$\mathbb{E}[Y_{i,t-2}\, \Delta\varepsilon_{it}] = 0.$$

This is the idea of Anderson and Hsiao (1981, 1982): pooled IV on the differenced equation, with $Y_{i,t-2}$ as the instrument, or $\Delta Y_{i,t-2}$, which uses one period fewer. Arellano and Bond (1991) push the idea to its limit: the later the month, the more valid lagged levels there are, and all of them are used. The differenced equation for period $t$ can use $Y_{i0}, \ldots, Y_{i,t-2}$, and the moment conditions stack into

$$\mathbb{E}[Y_{ij}\, \Delta\varepsilon_{it}] = 0, \qquad j \le t - 2.$$

::: {.intuition}
Why is the level of two months ago a clean instrument? The error of the differenced equation, $\Delta\varepsilon_{it}$, contains only this month's and last month's shocks. GMV two months ago was already settled: it knows everything that came before, but not what was to happen last month and this month; as long as shocks do not persist, it is unrelated to $\Delta\varepsilon_{it}$. It is, however, related to last month's change: a merchant that sold above its own long-run level two months ago most likely fell back last month, because deviations shrink when $\rho < 1$. So past levels predict the next change, and the instrument's strength comes from this pull back toward the long-run level; the closer $\rho$ is to 1, the weaker the pull, and the weaker the instrument.
:::

Write $\Delta\mathbf{Y}_i = (\Delta Y_{i2}, \ldots, \Delta Y_{iT})'$ and $\Delta\mathbf{Y}_{i,-1} = (\Delta Y_{i1}, \ldots, \Delta Y_{i,T-1})'$; the instrument matrix $\mathbf{Z}_i$ is a block-diagonal $(T-1) \times r$ matrix whose row for period $t$ holds $(Y_{i0}, \ldots, Y_{i,t-2})$, with zeros elsewhere, for $r = T(T-1)/2$ columns in all (Wooldridge 2010 equation (11.59)). In D1, $T = 7$ and $r = 21$.

::: {.assumption}
**AB.1 (rank condition; Wooldridge 2010 §11.6.2, p. 372 and p. 374)** $\mathbb{E}[\mathbf{Z}_i'\Delta\mathbf{Y}_{i,-1}] \neq 0$.
:::

AB.1 says that past levels help predict last month's change. The closer $\rho$ is to 1, the more GMV resembles a random walk and the harder it is for past levels to predict the next change; Wooldridge (2010) p. 374 points out that at $\rho = 1$ this correlation vanishes and both estimators fail. This is the subject of Section 3.5.

::: {.theorem}
**Theorem 3 (Anderson-Hsiao and Arellano-Bond; Wooldridge 2010 §11.6.2, p. 372; Hansen 2014 §19.3, equations (19.3) and (19.4); Anderson and Hsiao 1982; Arellano and Bond 1991)** Under DC.1, AB.1, and regularity conditions (finite fourth moments, independence across merchants): (i) the Anderson-Hsiao estimator, that is, pooled IV on the differenced equation with $Y_{i,t-2}$ or $\Delta Y_{i,t-2}$ as the instrument, is consistent and asymptotically normal under the corresponding rank condition (the instrument used is correlated with $\Delta Y_{i,t-1}$), with its variance clustered by merchant; (ii) for any positive definite $r \times r$ weight matrix $W$, the Arellano-Bond estimator

$$\begin{aligned}
\hat\rho_{AB} &= \arg\min_{\rho}\ \bar g_N(\rho)'\,W\,\bar g_N(\rho), \\
\bar g_N(\rho) &= \frac{1}{N}\sum_{i=1}^N \mathbf{Z}_i'\big(\Delta\mathbf{Y}_i - \rho\,\Delta\mathbf{Y}_{i,-1}\big)
\end{aligned}$$

is consistent and asymptotically normal; (iii) the one-step estimator $\hat\rho^{(1)}_{AB}$ takes $W^{(1)} = \big(\sum_i \mathbf{Z}_i'\mathbf{H}\mathbf{Z}_i\big)^{-1}$, where $\mathbf{H}$ is the matrix with 2 on the main diagonal, $-1$ on the diagonals next to it, and zeros elsewhere, of size $(T-1) \times (T-1)$; if $\varepsilon_{it}$ is also conditionally homoskedastic, it is the efficient GMM estimator for this set of instruments.
:::

The proof sketch has four steps. Step one (DC.1, Theorem 1 (ii)): every element of $\mathbf{Z}_i'\Delta\boldsymbol\varepsilon_i$ has mean zero, so the sample moments have expectation zero at the true value. Step two (AB.1): the moment conditions are linear in $\rho$, and their expectation at any candidate value equals the difference between the true and candidate values times $\mathbb{E}[\mathbf{Z}_i'\Delta\mathbf{Y}_{i,-1}]$, which is not zero, so the true value is the unique solution. Step three (regularity conditions, $T$ fixed): the sample moments are an average over $N$ independent merchants, so the law of large numbers and the central limit theorem give consistency and asymptotic normality, and the variance is the GMM sandwich of Chapter 7; Anderson-Hsiao is the special case with a single moment shared across periods. Step four (conditional homoskedasticity): the variance of $\Delta\varepsilon_{it}$ is $2\sigma_\varepsilon^2$, the covariance between adjacent periods is $-\sigma_\varepsilon^2$, and it is zero two periods apart, so $\mathrm{Var}(\Delta\boldsymbol\varepsilon_i) = \sigma_\varepsilon^2\mathbf{H}$, and $W^{(1)}$ is proportional to the inverse of the variance of the moments (the 3SLS form of Wooldridge 2010 p. 372).

In words: under DC.1, Anderson-Hsiao and Arellano-Bond estimate the same $\rho$, with no assumption about $\alpha_i$ and none about the initial conditions; the cost is that they use only the variation in differences, and the instruments, "past levels predicting current changes," have limited predictive power. When the error is serially correlated, some of the moment conditions fail and the estimators are generally inconsistent, converging to a point that depends on the form of the serial correlation. The weights are set by $W$: Anderson-Hsiao uses a single column of instruments, pooled across periods; one-step Arellano-Bond uses $\mathbf{H}$ to account for the correlation between adjacent differenced errors, and is efficient only when the errors are homoskedastic; the weights of the two-step estimator come from the data (Section 3.4). Back to Nimbus: on D1, Anderson-Hsiao gives 0.680 (standard error 0.062), with an interval [0.559, 0.802] that covers 0.70; one-step Arellano-Bond, with 21 instruments, gives 0.707 (0.048), and the standard error of the former is 1.283 times that of the latter, which is the payoff from using every lag. With $\Delta Y_{i,t-2}$ as the instrument, Anderson-Hsiao uses one month fewer and a weaker instrument, and gives 0.477 (0.145). On D1p the two versions of Anderson-Hsiao give 1.000 (0.277) and 0.312 (0.261), and the interval of the latter, [$-0.201$, 0.824], does not cover 0.90. This one miss in a single sample says nothing about bias. When the estimator is repeated 1,000 times under the design of D1p, the share of nominal 95% intervals that cover 0.90 is 0.933, a little under 95%, and every interval that misses lies below 0.90; the median estimate is 0.895, within Monte Carlo error of 0.90. The real problem is that the distribution is extremely dispersed: the interquartile range is 0.813, and the middle 95% of the estimates stretch all the way from 0.091 to 3.913. An exactly identified IV estimator has no finite mean, so only the median and quantiles are reported here. The instrument, then, is weak at $\rho = 0.90$: its standard error of 0.261 is much larger than the 0.145 at $\rho = 0.70$.

### 3.4 Inference for Two-Step GMM: The Windmeijer Correction, $m_2$, and Hansen J

The two-step estimator gives 0.7045 on D1. By the conventional formula its standard error is 0.0465; with the Windmeijer correction it is 0.0548, so the former is only 1/1.177 of the latter. On D1p the absolute gap is larger: the conventional standard error of two-step difference GMM is 0.1185, against 0.137 corrected. The specification tests are all there too: the test for first-order serial correlation in the differenced residuals gives $m_1 = -12.21$, $p < 10^{-4}$; the second-order test gives $m_2 = 0.964$, $p = 0.335$; Hansen J is 22.97 with 20 degrees of freedom, $p = 0.290$.

::: {.theorem}
**Theorem 4 (variance and specification tests of two-step GMM; Arellano and Bond 1991; Windmeijer 2005; Roodman 2009a §2.4, §2.5, and §3.5, equations (18), (27), and (28); for the general theory of two-step GMM and the J statistic, see Chapter 7)** Under DC.1, AB.1, and regularity conditions, let $\Delta\hat{\boldsymbol\varepsilon}_i$ be the differenced residuals of the one-step estimator, $g_i = \mathbf{Z}_i'\Delta\hat{\boldsymbol\varepsilon}_i$, and $\hat S = N^{-1}\sum_i g_i g_i'$.

\(i) The two-step estimator $\hat\rho^{(2)}_{AB}$ takes $W = \hat S^{-1}$ and is asymptotically efficient for this set of instruments, even when the errors are heteroskedastic.

\(ii) The conventional two-step variance treats $\hat S$ as known; Windmeijer's corrected variance adds the term that comes from the dependence of $\hat S$ on the one-step estimate. The two are asymptotically equivalent, and in finite samples the conventional variance is too small (the simulations of Windmeijer 2005).

\(iii) $m_l$ is the autocovariance of order $l$ of the differenced residuals, $N^{-1/2}\sum_i\sum_t \Delta\hat\varepsilon_{it}\Delta\hat\varepsilon_{i,t-l}$, divided by its standard error. When the differenced errors have no autocorrelation of order $l$, $m_l \xrightarrow{d} \mathcal{N}(0,1)$; under DC.1 this holds for $m_2$, while $\Delta\varepsilon_{it}$ and $\Delta\varepsilon_{i,t-1}$ share $\varepsilon_{i,t-1}$, so $m_1$ should be significantly negative.

\(iv) $J = N\,\bar g_N(\hat\rho^{(2)}_{AB})'\,\hat S^{-1}\,\bar g_N(\hat\rho^{(2)}_{AB})$; with $r$ fixed and all moment conditions true, $J \xrightarrow{d} \chi^2_{r-k}$, where $k$ is the number of parameters.
:::

The proof sketch has four steps. Step one (Theorem 3, Chapter 7): the one-step estimator is consistent, $\hat S$ consistently estimates the variance of the moments, and weighting by its inverse is efficient GMM. Step two (step one): the two-step estimator depends on the data not only through $\Delta\mathbf{Y}_i$ but also, through $\hat S$, on the one-step estimate; a first-order expansion of the two-step estimator around the one-step estimate produces an extra term proportional to the error of the one-step estimate, which the conventional formula drops and Windmeijer puts back, with its variance and covariances (Roodman 2009a equation (18)); this term is negligible as $N \to \infty$, so the two variances are asymptotically equivalent. Step three (DC.1): $\Delta\varepsilon_{it}\Delta\varepsilon_{i,t-2} = (\varepsilon_{it} - \varepsilon_{i,t-1})(\varepsilon_{i,t-2} - \varepsilon_{i,t-3})$ involves no common shock, so it has mean zero when the errors are uncorrelated; the central limit theorem applied to sums over merchants needs only a large $N$; replacing errors by residuals changes the variance, and the standard error of Arellano and Bond (1991) accounts for this (Roodman 2009a equations (27) and (28)). For $l = 1$ the expectation is $-\mathbb{E}[\varepsilon_{i,t-1}^2] < 0$. Step four (DC.1, $r$ fixed): the overidentification theorem of Chapter 7.

In words: under DC.1 the two-step and one-step estimators estimate the same $\rho$; the two-step estimator only changes the weights, and the Windmeijer correction changes only the standard error, not the point estimate. $m_2$ tests one implication of DC.1 (no serial correlation in the errors), and $J$ tests all $r - k$ overidentifying restrictions; passing both tests is only a necessary condition. When the specification is wrong, a rejection by $m_2$ says that $Y_{i,t-2}$ is no longer clean and the instruments must start from more distant lags; a rejection by $J$ says that some moment condition fails, but its power depends on how $r$ compares with $N$ (Section 3.6). The weights are set by $\hat S$: the larger the estimated variance of a moment, the lighter its weight in the two-step estimator and in $J$, and Roodman (2009b) warns that the moment that fails most badly may end up with the lightest weight. Back to Nimbus: D1's two-step estimate is 0.7045, with a Windmeijer interval of [0.597, 0.812]; $m_2$ has a $p$-value of 0.335 and Hansen J a $p$-value of 0.29, consistent with a correctly specified model. D1p's conventional standard error of 0.1185 is likewise too small, and whenever this chapter reports a two-step estimate it pairs it with the Windmeijer standard error; for the new data blocks the main estimate is one-step with robust standard errors, and the two-step estimate serves as a reference. Even under a correct specification a $p$-value can be small: all of D4's Arellano-Bond moment conditions hold, yet the $p$-value of its two-step Hansen J is 0.0519, right at the edge of 5%.

### 3.5 Persistence Near 1: Weak Instruments and Blundell-Bond System GMM

In D1p, $\rho = 0.90$. Two-step difference GMM gives 0.809, with Windmeijer standard error 0.137; two-step system GMM gives 0.911, with Windmeijer standard error 0.025, so the standard error of difference GMM is 5.46 times that of system GMM. The conventional standard error of system GMM is only 0.017, likewise too small. The reason lies in AB.1: write the model in terms of deviations from the long-run level, $\Delta Y_{i,t-1} = (\rho - 1)\big(Y_{i,t-2} - \alpha_i/(1-\rho)\big) + \varepsilon_{i,t-1}$; past levels can predict only the first term, and when $\rho$ is close to 1 the factor $\rho - 1$ squeezes that term to very little; what remains is the unpredictable shock. Blundell and Bond (1998), building on Arellano and Bover (1995), add a new set of moment conditions: turn things around and use past changes as instruments for the levels equation.

::: {.assumption}
**BB.1 (mean-stationary initial conditions; Blundell and Bond 1998; statement following Roodman 2009b, working-paper version, equation (10), and Wooldridge 2010 p. 374)** $\mathbb{E}\big[\big(Y_{i0} - \alpha_i/(1-\rho)\big)\,\alpha_i\big] = 0$.
:::

BB.1 says that a merchant's distance from its own long-run level when it enters the sample is unrelated to its talent. It does not require merchants to be at their long-run level already, only that the direction and size of the deviation be uncorrelated with $\alpha_i$; Wooldridge (2010) p. 374 puts it as the initial value being drawn from the steady-state distribution, which is a sufficient condition. In platform data it does not hold automatically: newly onboarded merchants are all ramping up, and high-talent merchants have higher targets and start farther from them; block D4 is built this way.

::: {.intuition}
Picture two merchants onboarded in the same month, one of ordinary talent and one of outstanding talent. Both reach only eighty percent of their own long-run level in the first month; the outstanding one has the higher target and falls further short of it, so it also ramps up faster over the following months. In the first months of the sample, GMV growth is then positively correlated with talent, and talent is part of the levels-equation error $\alpha_i + \varepsilon_{it}$: past changes are no longer clean instruments for the levels equation. If the two started equally far from their own targets, growth would be unrelated to talent. What BB.1 needs is not that merchants have already reached their long-run level, but that their distance from it is unrelated to talent.
:::

::: {.theorem}
**Theorem 5 (the level moment conditions of system GMM; Arellano and Bover 1995; Blundell and Bond 1998; Wooldridge 2010 p. 374; derivation following equations (9) and (10) of the working-paper version of Roodman 2009b)** Under DC.1 and $|\rho| < 1$, for $t = 2, \ldots, T$,

$$\begin{aligned}
&\mathbb{E}\big[\Delta Y_{i,t-1}(\alpha_i + \varepsilon_{it})\big] \\
&\quad = -(1-\rho)\,\rho^{t-2}\,\mathbb{E}\Big[\Big(Y_{i0} - \frac{\alpha_i}{1-\rho}\Big)\alpha_i\Big].
\end{aligned}$$

\(i) Under BB.1 the moment conditions of the levels equation, $\mathbb{E}[\Delta Y_{i,t-1}(\alpha_i + \varepsilon_{it})] = 0$, all hold, and system GMM, which stacks them with the differenced moment conditions of Theorem 3, is consistent and asymptotically normal. (ii) When BB.1 fails and $\rho \neq 0$, all of the level moment conditions fail and system GMM is generally inconsistent, while difference GMM is unaffected.
:::

The proof sketch has three steps. Step one (the model): $\alpha_i = (1-\rho)\cdot\alpha_i/(1-\rho)$, so the deviation from the long-run level satisfies $Y_{it} - \alpha_i/(1-\rho) = \rho\big(Y_{i,t-1} - \alpha_i/(1-\rho)\big) + \varepsilon_{it}$, and iterating gives $Y_{i,t-2} - \alpha_i/(1-\rho) = \rho^{t-2}\big(Y_{i0} - \alpha_i/(1-\rho)\big) + \sum_{j=1}^{t-2}\rho^{t-2-j}\varepsilon_{ij}$. Step two (DC.1): $\alpha_i$ is in the conditioning set, so $\mathbb{E}[\varepsilon_{ij}\alpha_i] = 0$; $\Delta Y_{i,t-1}$ is a function of the conditioning set of period $t$, so $\mathbb{E}[\Delta Y_{i,t-1}\varepsilon_{it}] = 0$. Step three (steps one and two): multiply $\Delta Y_{i,t-1} = (\rho - 1)\big(Y_{i,t-2} - \alpha_i/(1-\rho)\big) + \varepsilon_{i,t-1}$ by $\alpha_i + \varepsilon_{it}$ and take expectations; only the term in the initial deviation remains, giving the expression in the theorem; BB.1 makes it zero, which is (i), and otherwise (ii) follows.

In words: system GMM estimates the same $\rho$, but it bets on one more assumption, about the initial conditions. When the assumption holds, the levels equation uses past changes to predict current levels, and this set of instruments remains informative when $\rho$ is close to 1; when it is wrong, system GMM is biased instead, because the violation decays only as $\rho^{t-2}$, and the more persistent the series, the slower the decay. As Roodman (2009b) puts it, system GMM is most likely to fail exactly where it is most useful. The weights are set by the relative information in the two sets of moment conditions: the weaker the instruments of the differenced equation, the more the levels equation weighs in the estimate, and the larger the bias once the level moments are invalid. Back to Nimbus: D1p has a stationary start, and the interval [0.862, 0.960] around system GMM's 0.911 covers 0.90. D4's merchants start at 80% of their long-run level, and one-step system GMM gives 0.816, with an interval [0.798, 0.834] far from 0.6; in a large sample of twenty thousand merchants the one-step estimate is 0.817 and the two-step estimate 0.806 (the two-step estimate on this block is 0.805), so the bias does not vanish with the number of merchants, and the one-step and two-step estimators each have their own pseudo-true value. One-step difference GMM gives 0.593, with an interval [0.535, 0.651] that covers 0.6; it is a consistent estimator whose limit is 0.6 itself, and in the large sample of twenty thousand merchants the one-step and two-step estimates are both 0.607, a gap that is only sampling error. Here $r = 34$, far fewer than the 2,000 merchants, so Hansen J has power, and it rejects the moment conditions of system GMM with $p < 10^{-4}$. With the same random draws but a start at the long-run level, system GMM gives 0.618, with an interval [0.575, 0.660] that covers 0.6 and a Hansen J $p$-value of 0.231. One more thing to watch out for: Bun and Windmeijer (2010) prove (as stated in their 2009 working-paper version) that in a covariance-stationary AR(1) panel where the merchant effect and the idiosyncratic error have equal variances, the expected concentration parameter in each period's cross-section is the same for the levels equation and the differenced equation, so the levels equation suffers from weak instruments as well. D1p happens to be the case in which both variances are 1: system GMM's greater precision does not mean its instruments are strong.

### 3.6 Instrument Proliferation: The More Instruments, the More GMM Resembles the Estimator It Meant to Avoid

Block D2 has 70 key accounts and 16 months. One-step system GMM with month effects, using every available lag in every month, has 134 instruments, more than there are merchants: the estimate is 0.820, and its interval [0.760, 0.879] does not cover 0.70; the $J$ of this one-step fit is 70.00 with 118 degrees of freedom, $p = 0.9999$, which looks impeccable. Collapse the instruments (each lag order merged into one column across months) down to 30, and the estimate becomes 0.629, with an interval [0.530, 0.729] that covers the truth, and $J$ has a $p$-value of 0.121. The figure below adds instruments step by step, from lag 2 alone to every lag: the uncollapsed estimate rises steadily as instruments are added, moving toward the 0.929 of OLS in levels; FE is 0.549.

::: {.intuition}
Overfitting has an extreme version: in two-stage least squares with as many instruments as observations, the first-stage $R^2$ is 1, the fitted values are the endogenous regressor itself, and the second stage degenerates into OLS (Roodman 2009b). With fewer instruments, the first stage still fits some of the part of the endogenous regressor that should have been purged, namely the part correlated with the error; the more instruments, the more of it gets fitted, and the estimator slides bit by bit toward the estimator that uses no instruments. The failure of the test is the other side of the same coin: once the instruments are numerous enough to "explain" away every merchant's sample moments, $J$ can no longer measure how far the moment conditions are from zero.
:::

![Block D2 (70 key accounts, 16 months): estimates of $\rho$ with 95% intervals from one-step system GMM with month effects (left) and the one-step $J$ statistic's $p$-value (right), against the number of instruments $r$, with the lag depth of the instruments running from lag 2 alone to every lag. Dots joined by lines are uncollapsed instruments (one column per lag per month), and triangles joined by lines are collapsed instruments (one column per lag). The three horizontal lines in the left panel are the true value 0.70 (dashed), OLS in levels with month dummies, 0.929, and FE, 0.549 (dotted); the vertical dash-dotted line in both panels is $r = N = 70$. The horizontal dotted line in the right panel marks the 5% significance level, and in the gray region $r \ge N$, $J$ is identically $N$ and the $p$-value depends only on the degrees of freedom.](assets/fig/fig_11_proliferation.svg)

::: {.theorem}
**Theorem 6 (two consequences of instrument proliferation; the statement and proof are this chapter's own; for the phenomenon, see the working-paper version of Roodman 2009b, §III, and the simulations of Bowsher 2002)** Let $\mathbf{g}$ be the matrix with $g_i'$ in row $i$, of size $N \times r$, let $\boldsymbol\iota_N$ be the vector of all 1s, and let $J = \boldsymbol\iota_N'\mathbf{g}(\mathbf{g}'\mathbf{g})^{+}\mathbf{g}'\boldsymbol\iota_N$, where $(\cdot)^{+}$ is the Moore-Penrose generalized inverse. This $J$ is evaluated at the one-step estimate, with both the moments and the weight built from the one-step residuals; it is the $J$ that a one-step fit reports. The $J$ of Theorem 4 is evaluated at the two-step estimate, which minimizes the same quadratic form under the same weight, so it never exceeds this $J$. The definitions of the two $J$ statistics, the ordering between them, and (i) below rely only on the algebra of the moment vector. They therefore hold unchanged when $g_i$ stacks the differenced and level moments of system GMM (with month-effect instruments), as it does for block D2 in this section; the $\chi^2_{r-k}$ limit of Theorem 4 (iv) also applies there, again provided that $r$ is fixed and all the moment conditions hold.

\(i) If $r \ge N$ and $\mathbf{g}$ has rank $N$, then $J = N$, whatever the data; the $p$-value $\Pr(\chi^2_{r-k} > N)$ is determined by $r$, $k$, and $N$ alone. In this case the $J$ of Theorem 4 is at most $N$, and its $p$-value is at least this probability.

\(ii) Let $\boldsymbol\Delta$ be the $(T-1) \times T$ first-difference matrix, and $\mathbf{H} = \boldsymbol\Delta\boldsymbol\Delta'$. If one-step difference GMM has so many instruments that they span all the differenced observations (the stacked instrument matrix is an invertible square matrix; this is a limiting case, equivalent to giving each differenced observation its own instrument column, which the block-diagonal Arellano-Bond instruments cannot reach), the estimator is the GLS of the differenced equation with error variance $\mathbf{H}$, and this GLS is exactly FE.
:::

The proof sketch has three steps. Step one (algebra): $\mathbf{g}(\mathbf{g}'\mathbf{g})^{+}\mathbf{g}'$ is the orthogonal projection onto the column space of $\mathbf{g}$, which lies in $\mathbb{R}^N$; when the rank is $N$ the projection is the identity, $J = \boldsymbol\iota_N'\boldsymbol\iota_N = N$, and (i) follows; rank $N$ requires $r \ge N$. Step two (algebra): when the stacked instrument matrix $\mathbf{Z}$ is invertible, the one-step weight gives $\mathbf{Z}(\mathbf{Z}'\mathbf{H}\mathbf{Z})^{-1}\mathbf{Z}' = \mathbf{H}^{-1}$ (here $\mathbf{H}$ is block diagonal by merchant), and the estimator becomes $\big(\Delta\mathbf{Y}_{-1}'\mathbf{H}^{-1}\Delta\mathbf{Y}_{-1}\big)^{-1}\Delta\mathbf{Y}_{-1}'\mathbf{H}^{-1}\Delta\mathbf{Y}$, which is GLS. Step three (algebra): $\boldsymbol\Delta\boldsymbol\iota_T = 0$ and $\boldsymbol\Delta$ has rank $T-1$, so $\boldsymbol\Delta'(\boldsymbol\Delta\boldsymbol\Delta')^{-1}\boldsymbol\Delta = I_T - T^{-1}\boldsymbol\iota_T\boldsymbol\iota_T'$, which is exactly within demeaning; the GLS is the FE of Chapter 10, and (ii) follows. When the instruments fall short of this limit, the first stage gets ever closer to the endogenous regressor itself as instruments are added, and the estimator slides step by step toward the limit: this is what the many-moment bias of Chapter 7 looks like in a dynamic panel.

In words: whether $J$ carries any information depends on how $r$ compares with $N$. When $r \ge N$, it is an algebraic fact that the $J$ of a one-step fit is identically $N$ and that the two-step $J$ is at most $N$; the $p$-value is at least $\Pr(\chi^2_{r-k} > N)$, and with few parameters the test cannot reject, so "Hansen does not reject" is not evidence. When $r$ is below $N$ but close to it, the power of the test is also badly diminished. The estimate itself is dragged away too: the more instruments, the more GMM resembles the estimator that uses no instruments. Difference GMM is pulled toward FE, and so downward with FE's Nickell bias; system GMM contains both the differenced equation and the levels equation, and its limit is a combination of the two equations' uninstrumented estimators; on block D2 most of the weight falls on OLS in levels, so it is biased upward. When the specification is wrong, cutting instruments cannot rescue the estimate: the twin block's errors are MA(1), so $Y_{i,t-2}$ and the instruments of the levels equation are themselves invalid; the mean estimate is 0.889 with all instruments and 1.016 collapsed, and neither is on target for $\rho$. The twin block is useful for studying the power of the test, not for estimation. The weights are set by the number of instruments the researcher chooses: lag depth and collapsing decide how closely the first stage can fit the endogenous regressor. Back to Nimbus: in the 500 Monte Carlo replications of D2, the average estimate is 0.809 with all instruments and 0.697 collapsed, and OLS in levels averages 0.927; in 99.0% of the samples the full-instrument estimate is closer to OLS than the collapsed one, and the Hansen $p$-value with all instruments exceeds 0.99 in every sample.

The power of the test has to be read from the twin block, alongside the size; the test here is the two-step $J$ of Theorem 4. The test with all instruments never rejects, and neither does the test with the middle set (lags 2 to 4, $r = 68$, slightly fewer than the 70 merchants), whether the specification is correct or wrong: when $r$ is close to $N$, the test fails as well. The collapsed test rejects in 5.6% of the correctly specified samples, in line with the nominal 5%, and in all 200 samples of the twin block. If the $J$ of a one-step fit were used instead, the collapsed test would also reject in 26.6% of the correctly specified samples: the one-step weight of system GMM is not the optimal weight, and the one-step estimate does not minimize the quadratic form of Theorem 6, so this $J$ is too large and does not follow $\chi^2_{r-k}$.

The $p$-value curve in the figure is based on the $J$ of the one-step fits, and it illustrates Theorem 6 (i). As the uncollapsed instruments go from lag 2 alone to lags 2 to 4 ($r = 68$), the estimate rises from 0.731 to 0.748, and $J$ has a $p$-value of 0.139, still in the normal range; at lags 2 to 6, $r = 89$ already exceeds the 70 merchants, $J$ equals exactly 70 with 73 degrees of freedom, and the $p$-value is 0.578; at lags 2 to 9, $r = 113$, with 97 degrees of freedom and a $p$-value of 0.982; with every lag, the degrees of freedom are 118 and the $p$-value is 0.9999. $J$ does not move, the $p$-value climbs with the degrees of freedom all the way to 1, and the data say nothing new over these steps. The collapsed line is much steadier: across the six lag depths the estimates all lie between 0.611 and 0.636, and the $p$-value swings between 0.044 and 0.430; $J$ is not pinned at $N$ and is still responding to the data. This line, however, is also based on the $J$ of the one-step fits, which is too large for the collapsed set, so the 0.044 at lags 2 to 9 is not evidence that a moment condition fails. Whether the test can be trusted must be judged from the size and power of the two-step $J$ in the previous paragraph. So the value of cutting instruments lies not only in the point estimate but also in giving the test something to say again.

### 3.7 Predetermined Regressors and Feedback

Block D3 puts lagged GMV and promotion into the same equation, $Y_{it} = 0.6\,Y_{i,t-1} + 0.5\,X_{it} + \alpha_i + \varepsilon_{it}$, and merchants step up promotion after a month with a bad shock. FE estimates the promotion effect at 0.346, with an interval [0.321, 0.371] that does not cover 0.5, and persistence at 0.306. Promotion responds to last month's shock, so we first state the exogeneity that regressors of this kind satisfy.

::: {.assumption}
**SX.2 (sequential exogeneity with a lagged dependent variable; Wooldridge 2010 around equations (11.62) and (11.63), p. 372; Arellano and Bond 1991)** In the model $Y_{it} = \rho\,Y_{i,t-1} + X_{it}'\beta + \alpha_i + \varepsilon_{it}$, $\mathbb{E}[\varepsilon_{it} \mid Y_{i,t-1}, \ldots, Y_{i0}, X_{it}, \ldots, X_{i1}, \alpha_i] = 0$, $t = 1, \ldots, T$.
:::

An $X_{it}$ that satisfies SX.2 is called a predetermined regressor, the predetermined variable of Chapter 10: this month's promotion may be set using all the information up to last month, including last month's shock, but it cannot anticipate this month's shock.

::: {.theorem}
**Theorem 7 (predetermined regressors; Wooldridge 2010 §11.6.1 equations (11.53) to (11.57), §11.6.2 equations (11.62) and (11.63), pp. 369-372)** Under SX.2: (i) for the differenced equation $\Delta Y_{it} = \rho\,\Delta Y_{i,t-1} + \Delta X_{it}'\beta + \Delta\varepsilon_{it}$, the valid instruments in period $t$ are $(Y_{i0}, \ldots, Y_{i,t-2})$ and $(X_{i1}, \ldots, X_{i,t-1})$, and under the corresponding rank condition difference GMM with them as instruments is consistent and asymptotically normal; (ii) if $X_{it}$ responds to $\varepsilon_{i,t-1}$, $X_{it}$ itself is not a valid instrument in period $t$, and treating $X$ as strictly exogenous, so that $X_{it}$ also enters the instruments, makes the estimator inconsistent; (iii) FE is inconsistent for both $\rho$ and $\beta$, and the bias in $\beta$ comes both from feedback and from the Nickell bias of the lagged dependent variable, transmitted through the correlation between $X_{it}$ and $Y_{i,t-1}$.
:::

The proof sketch has three steps. Step one (SX.2): for $j \le t - 1$, $X_{ij}$ and $Y_{i,j-1}$ both lie in the conditioning sets of $\varepsilon_{it}$ and $\varepsilon_{i,t-1}$, so, as in Theorem 1, they are orthogonal to $\Delta\varepsilon_{it}$, and the argument of Theorem 3 then gives (i). Step two (feedback): $\mathbb{E}[X_{it}\Delta\varepsilon_{it}] = -\mathbb{E}[X_{it}\varepsilon_{i,t-1}] \neq 0$, so this moment condition fails, and GMM that uses it converges to the pseudo-true value that Chapter 7 describes under misspecification, which is (ii). Step three (the Frisch-Waugh-Lovell theorem): the $\beta$ of FE equals the coefficient from regressing the demeaned $Y_{it}$ on the residual left once the $Y_{i,t-1}$ part is partialled out of the demeaned $X_{it}$; feedback makes the demeaned $X_{it}$ correlated with the demeaned error (the mechanism of Theorem 4 (iii) of Chapter 10); the mechanism of Theorem 2 makes the demeaned lag correlated with the error; and since $X_{it}$ and $Y_{i,t-1}$ both move with $\alpha_i$ and past shocks, both correlations are transmitted into $\beta$, which is (iii).

In words: under SX.2, difference GMM with lagged $Y$ and lagged $X$ as instruments estimates the same-month effect of promotion, $\beta$, and the persistence $\rho$; the long-run effect is $\beta/(1-\rho) = 1.25$, so every USD 100 of promotion brings USD 1,250 of GMV in the long run. Misspecification can go two ways: treating $X$ as strictly exogenous uses invalid moments; treating $X$ as endogenous and using only lags two or more periods back keeps every moment valid but throws information away. The weights are set by the researcher's classification of $X$: the classification decides which moments enter the estimation. Back to Nimbus: with promotion treated as predetermined, the one-step estimate of $\beta$ is 0.513, interval [0.462, 0.564], and that of $\rho$ is 0.610; the Hansen J $p$-value is 0.668, and $m_2$ has a $p$-value of 0.605; the two-step estimate is 0.508, with a Windmeijer standard error of 0.0264. Treated as endogenous, $\beta$ is 0.472, interval [0.361, 0.583], with a standard error of 0.057, about twice that under the predetermined treatment. Treated as strictly exogenous, $\beta$ is 0.835, interval [0.780, 0.890], and Hansen J rejects with $p < 10^{-4}$; here $r = 54$, far fewer than the 4,000 merchants, so the test has power. The bias of FE can be taken apart: in a large sample of a hundred thousand merchants with the same random draws, the $\beta$ of FE is 0.367 with feedback and still 0.466 with feedback switched off, the latter due entirely to the lagged dependent variable; the two large-sample estimates of $\rho$ are 0.313 and 0.355. So "adding a lagged dependent variable to control for dynamics" and then running FE biases even a regressor that had no problem to begin with.

### 3.8 Bias-Corrected FE: Hahn-Kuersteiner and the Split-Panel Jackknife

There is one more route: since the bias of FE has a formula, subtract it. On D1 the analytical correction of Hahn-Kuersteiner pushes 0.425 up to 0.628, but its interval [0.586, 0.670] still does not cover 0.70. The split-panel jackknife (written half-panel jackknife in the code output and the figure legends) needs no formula: cut the 7 regression periods into a first half and a second half; 7 cannot be split evenly, so take both cuts, the first 3 and the last 4, and the first 4 and the last 3; the FE estimates on the half panels are 0.0908 and 0.1970, and 0.1973 and 0.0392. Following Dhaene and Jochmans (2015), the two halves of each cut are weighted by their shares of the 7 periods (3/7 and 4/7), and the two cuts are then averaged, giving 0.1405, so $2 \times 0.4246 - 0.1405 = 0.709$. The half panels have only three or four periods each and a much larger Nickell bias, and that is what lets them estimate the bias of the full panel.

::: {.assumption}
**HK.1 (Conditions 1 to 3 of Hahn and Kuersteiner 2002, numbered as in the 2001 working-paper version, univariate case)** NK.1 holds and the errors have moments of high enough order; $N$ and $T$ go to infinity together, with $N/T$ converging to a positive finite constant.
:::

::: {.theorem}
**Theorem 8 (bias correction; Hahn and Kuersteiner 2002, numbered as in the 2001 working-paper version, Theorem 1, Theorem 2, and equation (6); Dhaene and Jochmans 2015, numbered as in the accepted manuscript, equations (2.2), (2.3), and (2.4) and Theorem 2.1)** Under HK.1:

\(i) $\sqrt{NT}(\hat\rho_{FE} - \rho) \xrightarrow{d} \mathcal{N}\big(-\sqrt{\lim N/T}\,(1+\rho),\ 1-\rho^2\big)$;

\(ii) $\hat\rho_{HK} = \hat\rho_{FE} + (1 + \hat\rho_{FE})/T = \frac{T+1}{T}\hat\rho_{FE} + \frac{1}{T}$ satisfies $\sqrt{NT}(\hat\rho_{HK} - \rho) \xrightarrow{d} \mathcal{N}(0, 1-\rho^2)$;

\(iii) Cut the $T$ regression periods in time into a first and a second half, each using all merchants, and write the FE estimates on the two halves as $\hat\rho_{FE,1}$ and $\hat\rho_{FE,2}$. When $T$ is even, $\hat\rho_{SPJ} = 2\hat\rho_{FE} - (\hat\rho_{FE,1} + \hat\rho_{FE,2})/2$; when $T$ is odd, one half is longer than the other, the $(\hat\rho_{FE,1} + \hat\rho_{FE,2})/2$ in the formula is replaced by the average weighted by each half's share of the $T$ periods, and the result is then averaged over the two cuts, short-then-long and long-then-short (equations (2.2) and (2.4)). $\hat\rho_{SPJ}$ removes the bias of order $1/T$, and its asymptotic variance is the same as that of FE.
:::

The proof sketch has three steps. Step one (Theorem 2, HK.1): the bias of FE is $-(1+\rho)/T + O(T^{-2})$; multiplied by $\sqrt{NT}$ it gives $-\sqrt{N/T}(1+\rho)$, which does not go to zero, so the asymptotic distribution is centered away from zero, which is (i); the higher-order terms go to zero once multiplied by $\sqrt{NT}$. Step two (step one): HK subtracts the estimated bias $-(1+\hat\rho_{FE})/T$; the error of this plug-in is $(\hat\rho_{FE} - \rho)/T$, which goes to zero when multiplied by $\sqrt{NT}$, and this gives (ii). Step three (step one): a half with a given number of periods has a bias of about $-(1+\rho)$ divided by that number; weighted by their shares of the $T$ periods, the two halves together come to exactly $-2(1+\rho)/T$, and this holds for odd $T$ as well (with unequal halves, an equally weighted average would leave an extra $O(T^{-3})$ term). The bias of $2\hat\rho_{FE}$ is also $-2(1+\rho)/T$; subtracting the two removes the term of order $1/T$, which is (iii). This step uses only the fact that the bias is inversely proportional to the number of periods, not the constant of proportionality, so the jackknife does not depend on the AR(1) bias formula.

In words: both corrections still estimate $\rho$, but they remove the bias only under double asymptotics in which $N$ and $T$ grow at the same rate; with $T$ fixed they remove only the part of order $1/T$, and the remaining $O(T^{-2})$ can still be sizable in short panels. When the specification is wrong: the HK expression is that of an AR(1) with no other regressors and iid errors, and with regressors or heteroskedasticity the form of the bias changes. With regressors, the LSDVC of Kiviet (1995) can be used (corrected LSDV, where LSDV stands for least squares dummy variables, the regression on unit dummies, which is the FE of Theorem 3 (iii) of Chapter 10). Kiviet derives an approximation to the FE bias for balanced panels with strictly exogenous regressors; the unknown parameters in it are first replaced by a consistent initial estimate, Arellano-Bond for example, and the estimated bias is then subtracted. The jackknife requires the first and second halves to share the same bias structure, and in a panel with a nonstationary start, like D4, the two halves are not comparable. The weights are the same as those of FE: only within-unit variation is used. Back to Nimbus: at $T = 7$ the probability limit of HK is 0.6185 and that of the jackknife 0.6747. Block D5 repeats the design 1,000 times at each of 8, 16, and 32 months: the average HK estimates are 0.619, 0.677, and 0.694, and the jackknife's 0.675, 0.709, and 0.706; for nominal 95% intervals, the coverage of HK is only 2.2%, 46.9%, and 84.7%, that of the jackknife (with FE standard errors) 62.5%, 74.7%, and 79.0%, and that of FE zero at all three lengths. At 8 months the root mean squared error of the jackknife, 0.037, is clearly smaller than that of HK, 0.084, but its average estimate still falls 0.025 below 0.70, so the bias has not been removed completely; at 32 months the two have similar root mean squared errors. In short panels, bias correction can push the point estimate toward the truth, but it cannot deliver credible intervals.

## 4 Estimation and Practice

### 4.1 The Estimators at a Glance

| Estimator | Assumptions needed | What it estimates | D1 (truth 0.70) | D1p (truth 0.90) |
|---|---|---|---|---|
| Pooled OLS | Merchant effect unrelated to the lag (impossible) | Persistence plus the share of talent | 0.946 | 0.995 |
| FE | Strict exogeneity (impossible) | $\rho$ plus Nickell bias | 0.425 | 0.535 |
| FE + Hahn-Kuersteiner | NK.1; $N$ and $T$ of the same order | Removes the bias of order $1/T$ | 0.628 | 0.755 |
| Split-panel jackknife | Stationarity; $N$ and $T$ of the same order | Removes the bias of order $1/T$ | 0.709 | 0.820 |
| Anderson-Hsiao ($Y_{i,t-2}$) | DC.1 | $\rho$ | 0.680 | 1.000 |
| One-step difference GMM | DC.1; efficiency also needs homoskedasticity | $\rho$ | 0.707 | 0.789 |
| Two-step difference GMM | DC.1 | $\rho$ | 0.7045 | 0.809 |
| Two-step system GMM | DC.1 plus BB.1 | $\rho$ | 0.701 | 0.911 |

Which to choose depends first on two questions: is the error still serially correlated once the dynamics are fully specified, and is merchants' distance from their own long-run level when they enter the sample unrelated to their talent? The first decides which lags can be used, the second whether system GMM can be used. One-step or two-step is a third and much less important choice: on D1 one-step gives 0.707 and two-step 0.7045, hardly different. The two-step estimator is asymptotically more efficient when the errors are heteroskedastic, and the simulations of Windmeijer (2005) reported by Roodman (2009a) show that the two-step estimator with corrected standard errors does slightly better than the cluster-robust one-step estimator. With many instruments and few merchants, the two-step weight matrix is poorly estimated, and once $r$ exceeds $N$ it becomes singular; the one-step estimator is then more stable, which is why block D2 reports only one-step estimates. When $T$ is long (several dozen periods or more), the bias of FE, about $-(1+\rho)/T$, is already small, and bias-corrected FE is often simpler and more reliable than GMM; platform panels with short $T$ and large $N$ are where GMM is most at home.

### 4.2 A Reporting Checklist for a Dynamic Panel Analysis

1. Report pooled OLS and FE first, as upper and lower bounds on persistence; Roodman (2009a) relays Bond's advice that a credible estimate should lie between the two or near them.
2. Report the number of instruments $r$ and the number of merchants $N$; by default collapse the instruments or limit the lag depth, and show how the estimates and tests change when instruments are cut (Roodman 2009b). When $r$ is close to or above $N$, the $J$ statistic's $p$-value carries no information.
3. For one-step estimates report robust standard errors; for two-step estimates report only Windmeijer standard errors.
4. Report $m_1$ and $m_2$; $m_1$ should be significantly negative, a non-rejection by $m_2$ is only a necessary condition, and when $m_2$ rejects, start the instruments from more distant lags.
5. Report the Hansen J statistic of the two-step fit, its degrees of freedom, and its $p$-value, even when the main estimate is one-step (Roodman 2009a §2.5); Roodman (2009b) warns that a $p$-value as high as 0.25 should already raise concern, and one close to 1 is often a sign of too many instruments. With system GMM, run a separate Difference-in-Hansen test on the instruments of the levels equation: it is the difference between the two Hansen J statistics with and without this set of instruments, and if the other moments hold, then under the null that this set is also valid, it is approximately $\chi^2$ with degrees of freedom equal to the number of instruments in the set.
6. Report difference GMM and system GMM side by side; when they differ by more than sampling error, go back and examine the initial conditions.
7. Include time effects to absorb shocks common to all merchants; otherwise the "independence across merchants" on which $m_l$ and $J$ rely does not hold.
8. For each regressor, state whether it is strictly exogenous, predetermined, or endogenous, and which instruments were chosen on that basis.

The many-moment bias of Chapter 7 matters especially here: the number of instruments grows quadratically with $T$, and once $T$ is even moderately large, $r$ is of the same order as $N$.

### 4.3 Software

In R, `plm` 2.6-7 provides `pgmm`, which covers the GMM estimators of this chapter: `pgmm(y ~ lag(y, 1) | lag(y, 2:99), data, effect = "individual", model = "onestep", transformation = "d")` is one-step difference GMM, `model = "twosteps"` gives two-step, `transformation = "ld"` gives system GMM, `collapse = TRUE` collapses the instruments, and `lag(y, 2:4)` limits the lag depth. Note that the default of `effect` is `"twoways"`, that is, with time effects; `vcovHC(fit)` gives the Windmeijer variance for a two-step estimate and the robust variance for a one-step estimate; `sargan(fit)` (the name is inherited from Sargan) gives the Hansen J of Theorem 4 (iv) for a fit with `model = "twosteps"` and the $J$ of Theorem 6 for a fit with `model = "onestep"`. The latter is too large for one-step system GMM (Section 3.6), so when the main estimate is one-step, take the Hansen J from a fit with `model = "twosteps"`. Finally, `mtest(fit, order = 2, vcov = vcovHC(fit))` is $m_2$. `sargan(fit, "onestep")` returns a quadratic form that is not scaled by the error variance and cannot be read as a Sargan statistic, so this chapter does not use it. FE and pooled OLS use `plm` with the `within` and `pooling` models, and cluster variances use `vcovHC(fit, method = "arellano", type = "HC0", cluster = "group")`. The Stata counterparts are `xtabond2` (Roodman 2009a) and the official `xtabond` and `xtdpdsys`. Alternative estimators include the maximum likelihood of Moral-Benito, Allison and Williams (2019), whose finite-sample performance may be better when $N$ is small and the series is persistent, and the LSDVC of Kiviet (1995).

## 5 Anchoring Papers

Each of the three anchoring papers corresponds to one thread of the chapter: where difference GMM and its specification tests come from, how instrument proliferation makes conclusions look sound, and how a management application handles predetermined regressors and serial correlation. Each paper is laid out in five parts: paper, method, data, results, and limitations.

### 5.1 Arellano and Bond (1991)

::: {.case}
Paper: "Some Tests of Specification for Panel Data: Monte Carlo Evidence and an Application to Employment Equations," The Review of Economic Studies 58(2):277-297. The difference GMM and serial correlation tests of Sections 3.3 and 3.4 of this chapter both come from it.

Method: first-difference the individual effect out of the dynamic equation, use every available lagged level as an instrument for the differenced equation, with the instrument set growing over time (Theorem 3), and combine these moment conditions by one-step and two-step GMM; the paper also proposes the accompanying tests for first-order and second-order serial correlation in the differenced residuals ($m_1$ and $m_2$), which judge the specification together with the overidentification test. The treatment in Wooldridge (2010) §11.6.2 takes this paper as its blueprint.

Data: an employment panel of UK firms. As Roodman (2009a) describes it, the sample is 140 UK firms surveyed once a year from 1976 to 1984, an unbalanced panel; the equation is a dynamic regression of a firm's employment on two of its own lags, on wages, capital, and industry output (with their lags), and on year dummies; adjusting employment is costly, hence the lags.

Results: Roodman (2009a) demonstrates this procedure on the same data: the coefficient on lagged employment is 1.045 by pooled OLS and 0.733 by the within estimator, and by Bond's account a credible estimate should lie between the two or near them; one-step difference GMM (robust standard errors) gives 0.686 (0.145), with 41 instruments and 611 observations, and the second-order serial correlation test has a $p$-value of 0.606. He also relays a finding from the Monte Carlo of the original paper: on simulated panels of 7 periods and 100 units, the serial correlation test is better than the overidentification test at detecting lagged instruments spoiled by autocorrelation, but when the autocorrelation falls to 0.2 it rejects only half the time.

Limitations: when persistence approaches 1, lagged levels are weak instruments, which is what Blundell and Bond later remedied; the instruments grow quadratically with $T$, and many-instrument bias and test failure are the problems most in need of caution in practice; the serial correlation tests have limited power against weak autocorrelation.
:::

What makes this paper exemplary is that it gives not just an estimator but a whole suite of diagnostics that make the estimator testable and falsifiable. The two $p$-values on D1, 0.335 for $m_2$ and 0.29 for Hansen J, come from following this workflow.

### 5.2 Roodman (2009b)

::: {.case}
Paper: "A Note on the Theme of Too Many Instruments," Oxford Bulletin of Economics and Statistics 71(1):135-158. The numbers below follow the May 2008 revision of its Center for Global Development Working Paper 125.

Method: the paper sets out four consequences of instrument proliferation: the instruments overfit the endogenous regressor and pull the estimate toward the estimator that uses no instruments; the optimal weight matrix is poorly estimated; two-step standard errors are too small; and the Hansen test is weakened. It explains what the Blundell-Bond initial condition actually requires (BB.1 in this chapter), and reviews and simulates two ways that researchers have used to cut instruments: limiting the lag depth, and collapsing, which replaces the period-by-period moment conditions $\mathbb{E}[Y_{i,t-l}\Delta\varepsilon_{it}] = 0$ (one for each $l$ and $t$) with one moment condition for each lag order $l$; the author designed collapsing independently when writing xtabond2.

Data: a Monte Carlo (100 units, $\rho$ of 0.8, 500 simulated panels per setting, with an adjustable correlation between the initial deviation and the individual effect), and replications of two early growth studies: Forbes (2000), which uses difference GMM to study income inequality and growth in 45 countries over five-year periods, and Levine, Loayza and Beck (2000), which uses system GMM to study financial development and growth.

Results: in the Monte Carlo, when the initial condition of system GMM is violated, the Hansen test with all instruments never rejects on 20-period panels, with an average $p$-value of 1.000, while the leanest instrument set, collapsed and using only the first lag, almost always rejects, with an average $p$-value of 0.017. The Forbes replication uses 80 instruments and 138 observations, with a Hansen $p$-value of 1.00; after collapsing to 30 instruments, the $t$-statistic of the inequality coefficient falls from 2.12 to 1.09 and loses significance. The replication of Levine et al. uses 75 instruments, 77 countries, and 353 observations, and the Difference-in-Hansen test of the system GMM instruments has a $p$-value of 0.75; after collapsing, the $p$-value of the overall Hansen test falls to 0.03, and that of the test on the lagged-growth instruments to 0.001. The paper also relays the simulations of Bowsher (2002): with 100 units, 15 periods, and 91 instruments, the Hansen test never rejects at the 5% or 10% significance level.

Limitations: there is no agreed safe number of instruments, and the author points out that both the simulations and the replications suggest that pushing $r$ below $N$ does not guarantee a valid test; cutting instruments increases variance; the replications depend on the data and specifications available, and their conclusion is that "the original results may come from undetected endogeneity"; they do not prove the opposite causal relation.
:::

This paper turned "the Hansen test does not reject" from a pass into a number that needs defending. The $J$ of block D2's one-step fit, identically 70, is the extreme case of the weakened test it describes: once the instruments outnumber the units, the $J$ of a one-step fit equals the number of units as a matter of algebra, the two-step $J$ cannot exceed it, and with few parameters the test cannot reject.

### 5.3 Chung, Kim and Park (2019)

::: {.case}
Paper: "How Do Sales Efforts Pay Off? Dynamic Panel Data Analysis in the Nerlove-Arrow Framework," Management Science 65(11):5197-5218. The description below follows the paper's abstract and its June 2017 Harvard Business School working-paper version (Working Paper 17-108).

Method: within the Nerlove-Arrow goodwill framework, the authors estimate the short-run and long-run returns to pharmaceutical sales representatives' visits to physicians: prescriptions depend on last period's prescriptions (a lagged dependent variable), and the number of visits is a predetermined regressor (visit plans take past prescriptions into account). The authors point out that when the error is serially correlated the usual dynamic panel moment conditions fail, while the usual serial correlation tests have weak power, so researchers may use invalid moment conditions without realizing it; they give moment conditions that are more robust for both the lagged dependent variable and the predetermined regressor, and illustrate them with theoretical examples and Monte Carlo simulations.

Data: visits by pharmaceutical sales representatives to different types of physicians and physician-level sales outcomes (monthly prescriptions in the working-paper version), forming a physician-level dynamic panel.

Results: the long-run persistence of the effect of visits is stronger for specialists, and the current-period marginal effect is larger for general practitioners; this distinction has direct implications for planning visit routes.

Limitations: the more robust moment conditions usually give up the lags closest to the current period, and the instruments become weaker as a result, which is the trade-off of Section 3.7; the conclusions may also depend on the form of goodwill decay under the Nerlove-Arrow specification.
:::

The paper speaks directly to two points in this chapter: the response of visits to past prescriptions is the feedback of D3, and the low power of the serial correlation tests is the empirical version of the point in Section 3.4 that "a non-rejection by $m_2$ is only a necessary condition."

Taken together: Arellano and Bond give a way to build instruments from inside the model, with diagnostics to match; Roodman shows that with many instruments these diagnostics fail without anyone noticing; and Chung, Kim and Park show how predetermined regressors and serial correlation appear together in real marketing and management problems, and how the moment conditions should be adjusted for them.

## 6 A Full Walkthrough on the Nimbus Data

The code below runs under R 4.6.1 and `plm` 2.6-7, with the random seed fixed block by block (D1 and D1p keep 808 and 7372, D1T keeps 303), and the result lines in the code that start with `#` copy what the run printed. Standard errors are clustered by merchant, without a small-sample factor (CR0); GMM standard errors are robust standard errors for one-step estimates and Windmeijer standard errors for two-step estimates. Section 6.1 first lists the helper functions shared by all blocks: `fe_ar1` (FE, returning the estimate, the standard error, and the number of regression periods `T`), `pooled_ar1`, `anderson_hsiao` (IV on the differenced equation, with $Y_{i,t-2}$ or $\Delta Y_{i,t-2}$ as the instrument), `gmm_summary` (which summarizes a `pgmm` fit: the estimate, the robust or Windmeijer standard error, the uncorrected two-step standard error, the number of instruments $r$, Hansen J, $m_1$, and $m_2$), `ci95`, and a few functions they all call; then comes the data-generating process, with the functions for Nickell bias and the two bias corrections right after the code that generates D1 and D1p.

### 6.1 The Data-Generating Process

The functions shared by the blocks of Section 6 follow. Merchant means average over the regression periods only, the cluster-robust variance sums by merchant (CR0), and functions whose names start with a dot are internal.

```r
library(plm)

.unit_index <- function(id) match(id, unique(id))                  # 1..G in order of appearance
.demean <- function(M, g) {                                         # subtract unit means (rows of M)
  M <- as.matrix(M)
  M - (rowsum(M, g) / as.vector(rowsum(rep(1, nrow(M)), g)))[g, , drop = FALSE]
}
.cl_vcov <- function(Z, res, g, bread = solve(crossprod(Z))) {
  S <- rowsum(Z * as.vector(res), g)
  V <- bread %*% crossprod(S) %*% bread
  dimnames(V) <- list(colnames(Z), colnames(Z))
  V
}
.panel_key <- function(id, t, pad = 0) {
  if (any(t != round(t))) stop("t must be integer-valued (months)")
  span <- max(t) - min(t) + pad + 1
  match(id, unique(id)) * span + (t - min(t))
}
.check_panel <- function(y, id, t) {
  if (length(y) != length(id) || length(y) != length(t)) stop("y, id and t must have the same length")
  if (anyNA(id) || anyNA(t)) stop("id and t must not contain NA")
  if (anyDuplicated(.panel_key(id, t))) stop("(id, t) pairs must be unique")
  invisible(TRUE)
}
ci95 <- function(est, se) c(est - qnorm(0.975) * se, est + qnorm(0.975) * se)
lag_within <- function(v, id, t, k = 1L) {
  key <- .panel_key(id, t, pad = abs(k))
  v[match(key - k, key)]
}
fe_ar1 <- function(y, id, t, X = NULL) {
  .check_panel(y, id, t)
  yl <- lag_within(y, id, t, 1L)
  R <- cbind(rho = yl)
  if (!is.null(X)) {
    X <- as.matrix(X); if (is.null(colnames(X))) colnames(X) <- paste0("x", seq_len(ncol(X)))
    R <- cbind(R, X)
  }
  ok <- !is.na(y) & stats::complete.cases(R)
  g <- .unit_index(id[ok])
  Yd <- .demean(y[ok], g); Rd <- .demean(R[ok, , drop = FALSE], g); colnames(Rd) <- colnames(R)
  A <- solve(crossprod(Rd)); b <- drop(A %*% crossprod(Rd, Yd)); names(b) <- colnames(R)
  e <- drop(Yd - Rd %*% b)
  V <- .cl_vcov(Rd, e, g, A)
  Tg <- tabulate(g)
  list(est = b, se = sqrt(diag(V)), V = V, n = length(e), G = max(g),
       T = if (length(unique(Tg)) == 1L) Tg[1] else NA_integer_)
}
pooled_ar1 <- function(y, id, t, time_fe = FALSE) {
  .check_panel(y, id, t)
  yl <- lag_within(y, id, t, 1L); ok <- !is.na(y) & !is.na(yl)
  Z <- if (time_fe) stats::model.matrix(~ yl + factor(t), data.frame(yl = yl[ok], t = t[ok])) else cbind(1, yl = yl[ok])
  colnames(Z)[2] <- "yl"
  A <- solve(crossprod(Z)); b <- drop(A %*% crossprod(Z, y[ok])); e <- drop(y[ok] - Z %*% b)
  V <- .cl_vcov(Z, e, .unit_index(id[ok]), A)
  list(est = unname(b[2]), se = sqrt(V[2, 2]), n = length(e))
}
anderson_hsiao <- function(y, id, t, inst = c("level", "diff")) {
  inst <- match.arg(inst)
  .check_panel(y, id, t)
  dy <- y - lag_within(y, id, t, 1L); dyl <- lag_within(dy, id, t, 1L)
  z <- if (inst == "level") lag_within(y, id, t, 2L) else lag_within(dy, id, t, 2L)
  ok <- !is.na(dy) & !is.na(dyl) & !is.na(z)
  zx <- sum(z[ok] * dyl[ok])
  est <- sum(z[ok] * dy[ok]) / zx
  e <- dy[ok] - est * dyl[ok]
  s <- rowsum(z[ok] * e, .unit_index(id[ok]))
  se <- sqrt(sum(s^2)) / abs(zx)
  list(est = est, se = se, ci = ci95(est, se), n = sum(ok), inst = inst)
}
n_instruments <- function(fit) {
  if (!inherits(fit, "pgmm")) stop("n_instruments: fit must be a plm::pgmm object")
  as.integer(ncol(fit$W[[1L]]))
}
gmm_summary <- function(fit, par = "lag(y, 1)") {
  V <- suppressWarnings(plm::vcovHC(fit))
  est <- unname(coef(fit)[par]); se <- unname(sqrt(V[par, par]))
  two <- identical(fit$args$model, "twosteps")
  J <- plm::sargan(fit, "twosteps")
  m1 <- suppressWarnings(plm::mtest(fit, order = 1L, vcov = V))
  m2 <- suppressWarnings(plm::mtest(fit, order = 2L, vcov = V))
  list(est = est, se = se, se_conv = if (two) unname(sqrt(stats::vcov(fit)[par, par])) else NA_real_,
       ci = ci95(est, se), r = n_instruments(fit), N = length(fit$W),
       J = unname(J$statistic), J_df = unname(J$parameter), J_p = unname(J$p.value),
       m1 = unname(m1$statistic), m1_p = unname(m1$p.value),
       m2 = unname(m2$statistic), m2_p = unname(m2$p.value))
}
.stack <- function(N, TT, cols) {
  out <- data.frame(id = rep(seq_len(N), each = TT), t = rep(seq_len(TT), times = N))
  for (nm in names(cols)) out[[nm]] <- if (is.matrix(cols[[nm]])) as.vector(t(cols[[nm]])) else rep(cols[[nm]], each = TT)
  out
}
```

The generating code of the inherited blocks D1 and D1p is kept as it was. The `T` in this code is the number of months per merchant, the $T + 1$ of the text; the start is drawn from the stationary distribution, and then 50 months are discarded as burn-in to wash out the influence of the initial condition.

```r
## ---- dynamic panel: AR(1) with merchant effects ----------------------------
# stationary start; burn-in to erase initial-condition transient
gen_dyn <- function(N, T, rho, burn = 50, sigma_a = 1, sigma_e = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  alpha <- rnorm(N, 0, sigma_a)
  Ttot  <- T + burn
  y <- matrix(0, N, Ttot)
  y[,1] <- alpha/(1-rho) + rnorm(N, 0, sigma_e/sqrt(1-rho^2))  # stationary draw
  for (s in 2:Ttot) y[,s] <- rho*y[,s-1] + alpha + rnorm(N, 0, sigma_e)
  y <- y[, (burn+1):Ttot]                              # keep last T
  data.frame(id = rep(1:N, each = T), t = rep(1:T, times = N),
             y = as.vector(t(y)))
}
dyn_main <- gen_dyn(N = 500, T = 8, rho = 0.70, seed = 808)    # main
dyn_persist <- gen_dyn(N = 500, T = 8, rho = 0.90, seed = 7372) # rho near 1
```

Nickell's exact bias, the Hahn-Kuersteiner correction, and the split-panel jackknife (cut into two halves) are written as a few functions; here `T` is the number of regression periods. `spj_parts` weights the two halves following equations (2.2) and (2.4) of Dhaene and Jochmans (2015), each half's weight being its share of the regression periods.

```r
nickell_bias <- function(rho, T) {             # exact plim of rho_hat_FE - rho; T = regression periods
  if (any(T < 2)) stop("nickell_bias: T (regression periods) must be at least 2")
  if (any(abs(rho) >= 1)) stop("nickell_bias: |rho| must be below 1")
  a <- (1 - rho^T) / (1 - rho)
  -(1 + rho) / (T - 1) * (1 - a / T) / (1 - 2 * rho / ((1 - rho) * (T - 1)) * (1 - a / T))
}
hk_correct <- function(rho_fe, T) ((T + 1) / T) * rho_fe + 1 / T     # Hahn-Kuersteiner
spj_parts <- function(y, id, t) {              # FE on all T periods and on two halves in time
  .check_panel(y, id, t)
  yl <- lag_within(y, id, t, 1L); ok <- !is.na(y) & !is.na(yl)
  periods <- sort(unique(t[ok])); TT <- length(periods)
  if (TT < 2) stop("spj_fe: need at least two regression periods")
  fe_on <- function(keep) {
    g <- .unit_index(id[keep]); Yd <- .demean(y[keep], g); Ld <- .demean(yl[keep], g)
    sum(Ld * Yd) / sum(Ld^2)
  }
  fe <- fe_on(ok)
  cuts <- unique(c(floor(TT / 2), ceiling(TT / 2)))                 # odd T: two cuts
  halves <- do.call(rbind, lapply(cuts, function(k) {
    s1 <- periods[seq_len(k)]
    c(S1 = fe_on(ok & t %in% s1), S2 = fe_on(ok & !(t %in% s1)))
  }))
  weights <- cbind(S1 = cuts / TT, S2 = (TT - cuts) / TT)             # shares of the T periods, DJ (2.2)
  half_mean <- mean(rowSums(weights * halves))                        # averaged over the cuts, DJ (2.4)
  list(fe = fe, halves = halves, weights = weights, half_mean = half_mean, spj = 2 * fe - half_mean, T = TT)
}
spj_fe <- function(y, id, t) spj_parts(y, id, t)$spj
```

The parameters and generating functions of the new blocks follow. `gen_D2` generates a stationary AR(1): each merchant's long-run level (in thousands of USD) is normal with mean 40 and standard deviation 5, $\alpha_i$ equals $1 - \rho$ times the long-run level, the start is drawn from the stationary distribution, and 50 months are discarded; `gen_D2_ma` only changes the error to $e_{it} + 0.6\,e_{i,t-1}$. In the code, `eps_x` is promotion's own random disturbance, `start_share` is the share of the long-run level reached in the first month, and `w` is the random deviation in the first month; they appear only in the code.

```r
P_D2 <- list(N = 70, T = 16, rho = 0.70, m_mean = 40, m_sd = 5, sigma_e = 2, burn = 50, R_mc = 500,
             lags = c(2, 3, 4, 6, 9, 99), mid_lag = 4)
P_D2_ma <- modifyList(P_D2, list(ma_coef = 0.6, R_mc = 200))   # MA(1) errors: Y[t-2] is an invalid instrument
P_D3 <- list(N = 4000, T = 8, rho = 0.6, beta = 0.5, a_mean = 3, a_sd = 1, x_mean = 6, x_alpha = 0.4,
             feedback = -0.6, eps_x_sd = 0.5, sigma_e = 1, burn = 50)
P_D4 <- list(N = 2000, T = 8, rho = 0.6, m_mean = 15, m_sd = 3, start_share = 0.8, sigma_e = 1)
gen_D2 <- function(P) {                        # key accounts: a stationary AR(1), Y in thousands of USD
  m <- stats::rnorm(P$N, P$m_mean, P$m_sd); alpha <- (1 - P$rho) * m
  Tt <- P$T + P$burn; Y <- matrix(0, P$N, Tt)
  Y[, 1] <- m + stats::rnorm(P$N, 0, P$sigma_e / sqrt(1 - P$rho^2))
  for (s in 2:Tt) Y[, s] <- P$rho * Y[, s - 1] + alpha + stats::rnorm(P$N, 0, P$sigma_e)
  .stack(P$N, P$T, list(y = Y[, (P$burn + 1):Tt], m_true = m, alpha_true = alpha))
}
gen_D2_ma <- function(P) {                     # the twin: MA(1) errors, Y[t-2] is no longer a valid instrument
  m <- stats::rnorm(P$N, P$m_mean, P$m_sd); alpha <- (1 - P$rho) * m
  Tt <- P$T + P$burn; Y <- matrix(0, P$N, Tt)
  e_prev <- stats::rnorm(P$N, 0, P$sigma_e)
  Y[, 1] <- m + stats::rnorm(P$N, 0, P$sigma_e * sqrt((1 + 2 * P$rho * P$ma_coef + P$ma_coef^2) / (1 - P$rho^2)))
  for (s in 2:Tt) {
    e <- stats::rnorm(P$N, 0, P$sigma_e)
    Y[, s] <- P$rho * Y[, s - 1] + alpha + e + P$ma_coef * e_prev
    e_prev <- e
  }
  .stack(P$N, P$T, list(y = Y[, (P$burn + 1):Tt], m_true = m, alpha_true = alpha))
}
gen_D3 <- function(P) {                        # GMV with a lag and a promotion budget that reacts to last month
  alpha <- stats::rnorm(P$N, P$a_mean, P$a_sd)
  Tt <- P$T + P$burn; Y <- X <- E <- matrix(0, P$N, Tt)
  E[, 1] <- stats::rnorm(P$N, 0, P$sigma_e)
  X[, 1] <- P$x_mean + P$x_alpha * alpha + stats::rnorm(P$N, 0, P$eps_x_sd)
  Y[, 1] <- (alpha + P$beta * (P$x_mean + P$x_alpha * alpha)) / (1 - P$rho) + E[, 1]
  for (s in 2:Tt) {
    E[, s] <- stats::rnorm(P$N, 0, P$sigma_e)
    X[, s] <- P$x_mean + P$x_alpha * alpha + P$feedback * E[, s - 1] + stats::rnorm(P$N, 0, P$eps_x_sd)
    Y[, s] <- P$rho * Y[, s - 1] + P$beta * X[, s] + alpha + E[, s]
  }
  k <- (P$burn + 1):Tt
  .stack(P$N, P$T, list(y = Y[, k], x = X[, k], eps_true = E[, k], alpha_true = alpha))
}
gen_D4 <- function(P) {                        # an onboarding cohort that starts at 80% of its long-run level
  m <- stats::rnorm(P$N, P$m_mean, P$m_sd); alpha <- (1 - P$rho) * m
  w <- stats::rnorm(P$N, 0, P$sigma_e / sqrt(1 - P$rho^2))
  Y <- Ys <- matrix(0, P$N, P$T)
  Y[, 1] <- P$start_share * m + w; Ys[, 1] <- m + w             # Ys: the same draws, mean-stationary start
  for (s in 2:P$T) {
    e <- stats::rnorm(P$N, 0, P$sigma_e)
    Y[, s] <- P$rho * Y[, s - 1] + alpha + e; Ys[, s] <- P$rho * Ys[, s - 1] + alpha + e
  }
  .stack(P$N, P$T, list(y = Y, y_stat_true = Ys, m_true = m, alpha_true = alpha))
}
set.seed(1102); dD2 <- gen_D2(P_D2)
set.seed(1103); dD3 <- gen_D3(P_D3)
set.seed(1104); dD4 <- gen_D4(P_D4)
```

### 6.2 Blocks D1 and D1p: The Hero Turns Traitor, and GMM Recovers

On the dynamic block, regress this month's GMV on last month's: first the two wrong estimators, then difference GMM and system GMM. The designs of D1 and D1p have no shocks common to all merchants, so the inherited specification has no month effects (`effect = "individual"`).

```r
d  <- pdata.frame(dyn_main, index = c("id", "t"))
pool <- plm(y ~ lag(y, 1), d, model = "pooling")
fe   <- plm(y ~ lag(y, 1), d, model = "within")
```

```r
ab <- pgmm(y ~ lag(y, 1) | lag(y, 2:99), d, effect = "individual",
           model = "twosteps", transformation = "d")    # difference GMM (Arellano-Bond)
bb <- pgmm(y ~ lag(y, 1) | lag(y, 2:99), d, effect = "individual",
           model = "twosteps", transformation = "ld")   # system GMM (Blundell-Bond)
round(c(pooled = coef(pool)[[2]], FE = coef(fe)[[1]], AB = coef(ab)[[1]], BB = coef(bb)[[1]]), 4)
# pooled     FE     AB     BB
# 0.9461 0.4246 0.7045 0.7005
```

The two wrong numbers of Section 1 are confirmed here: pooled OLS, at 0.946, is pushed up by talent; FE, at 0.425, is pressed down by Nickell bias; and two-step difference GMM and system GMM both land near 0.70. The same FE that was the hero recovering the truth in Section 6.2 of Chapter 10 becomes a villain here, and the gap is entirely due to the lagged dependent variable violating strict exogeneity. Below, the estimators of Sections 3.2 to 3.8 are all run on D1 and D1p.

```r
one_block <- function(dd) {
  pd <- pdata.frame(dd, index = c("id", "t"))
  gm <- function(model, tf) pgmm(y ~ lag(y, 1) | lag(y, 2:99), pd, effect = "individual", model = model,
                                 transformation = tf)
  f <- fe_ar1(dd$y, dd$id, dd$t)                                  # within estimator, CR0 by merchant
  list(pool = pooled_ar1(dd$y, dd$id, dd$t), fe = f, hk = hk_correct(f$est[["rho"]], f$T),   # f$T = 7
       spj = spj_fe(dd$y, dd$id, dd$t), ah = anderson_hsiao(dd$y, dd$id, dd$t, "level"),
       ah_diff = anderson_hsiao(dd$y, dd$id, dd$t, "diff"), ab1 = gmm_summary(gm("onestep", "d")),
       ab2 = gmm_summary(gm("twosteps", "d")), bb2 = gmm_summary(gm("twosteps", "ld")))
}
for (b in c("dyn_main", "dyn_persist")) {
  x <- one_block(get(b))
  cat(sprintf("%s: pooled %.4f (%.4f), FE %.4f (%.4f), FE + HK %.4f, half-panel jackknife %.4f\n", b, x$pool$est,
              x$pool$se, x$fe$est[["rho"]], x$fe$se[["rho"]], x$hk, x$spj))
  cat(sprintf("  Anderson-Hsiao: Y[t-2] %.4f (%.4f), dY[t-2] %.4f (%.4f); AB one-step %.4f (%.4f), r = %d\n",
              x$ah$est, x$ah$se, x$ah_diff$est, x$ah_diff$se, x$ab1$est, x$ab1$se, x$ab1$r))
  cat(sprintf("  two-step AB %.4f (Windmeijer %.4f, conventional %.4f), BB %.4f (Windmeijer %.4f, conventional %.4f)\n",
              x$ab2$est, x$ab2$se, x$ab2$se_conv, x$bb2$est, x$bb2$se, x$bb2$se_conv))
  cat(sprintf("  two-step AB tests: Hansen J %.2f on %d df (p %.3f), m1 %.2f, m2 %.3f (p %.3f)\n", x$ab2$J,
              as.integer(x$ab2$J_df), x$ab2$J_p, x$ab2$m1, x$ab2$m2, x$ab2$m2_p))
}
# dyn_main: pooled 0.9461 (0.0043), FE 0.4246 (0.0188), FE + HK 0.6281, half-panel jackknife 0.7087
#   Anderson-Hsiao: Y[t-2] 0.6803 (0.0621), dY[t-2] 0.4774 (0.1446); AB one-step 0.7067 (0.0484), r = 21
#   two-step AB 0.7045 (Windmeijer 0.0548, conventional 0.0465), BB 0.7005 (Windmeijer 0.0343, conventional 0.0265)
#   two-step AB tests: Hansen J 22.97 on 20 df (p 0.290), m1 -12.21, m2 0.964 (p 0.335)
# dyn_persist: pooled 0.9950 (0.0015), FE 0.5354 (0.0162), FE + HK 0.7548, half-panel jackknife 0.8202
#   Anderson-Hsiao: Y[t-2] 1.0002 (0.2771), dY[t-2] 0.3116 (0.2613); AB one-step 0.7892 (0.1210), r = 21
#   two-step AB 0.8094 (Windmeijer 0.1367, conventional 0.1185), BB 0.9113 (Windmeijer 0.0250, conventional 0.0169)
#   two-step AB tests: Hansen J 25.23 on 20 df (p 0.193), m1 -6.49, m2 1.416 (p 0.157)
```

The result lines match Section 3 one by one: FE and pooled OLS stop near their own probability limits, 0.416 and 0.955 (Theorem 2 and the expression of Section 3.2); Anderson-Hsiao and one-step difference GMM both aim at 0.70, the latter more precisely (Theorem 3); the Windmeijer standard errors of the two-step estimates are larger than the conventional ones (Theorem 4); and at $\rho = 0.90$ the standard error of difference GMM balloons while system GMM stays precise (Theorem 5).

![Blocks D1 ($\rho = 0.70$, left) and D1p ($\rho = 0.90$, right), 500 merchants and 8 months each: point estimates and 95% intervals of eight estimators, namely pooled OLS, FE, Hahn-Kuersteiner-corrected FE, the split-panel jackknife (point estimate only; written half-panel jackknife in the legend), Anderson-Hsiao (with $Y_{i,t-2}$ as the instrument), one-step Arellano-Bond, two-step Arellano-Bond, and two-step Blundell-Bond. Intervals are clustered by merchant, and the GMM intervals use one-step robust or two-step Windmeijer standard errors; hollow diamonds are the probability limits of pooled OLS and FE as the number of merchants goes to infinity (FE by the exact Nickell formula), and the vertical dashed line is the true value; the number at the right end of each row is that row's point estimate. The right panel has a wider horizontal axis, both to fit the Anderson-Hsiao interval, whose upper end reaches 1.543, and to leave room for the column of numbers on the right.](assets/fig/fig_11_estimators.svg)

### 6.3 Nickell Bias Fading with $T$, and the Bias Corrections

Fix $\rho = 0.70$, let the number of months per merchant grow from 4 to 32, and look at the average estimates of FE and of one-step difference GMM (collapsed, lags 2 to 4); the third row is the exact Nickell limit.

```r
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
nick <- rbind(nick, nickell_exact = 0.70 + nickell_bias(0.70, Ts - 1))   # T = months - 1 regression periods
print(round(nick, 4))
#                    4      8     16     32
# fe            0.0821 0.4169 0.5703 0.6400
# ab            0.7050 0.6864 0.6995 0.7014
# nickell_exact 0.0797 0.4162 0.5723 0.6413
```

| Months per merchant | 4 | 8 | 16 | 32 |
|---|---|---|---|---|
| Average FE estimate | 0.082 | 0.417 | 0.570 | 0.640 |
| Average one-step Arellano-Bond estimate | 0.705 | 0.686 | 0.700 | 0.701 |
| Exact Nickell limit of FE | 0.080 | 0.416 | 0.572 | 0.641 |

FE climbs steadily toward 0.70 as the number of months grows, and the exact limit traces this path closely; one-step difference GMM holds near 0.70 at every length, with no $1/T$ drift. These are only means over 40 samples of 300 merchants each per length; they are compatible with "GMM is consistent" but do not by themselves prove consistency, which comes from Theorem 3. The 0.417 in the 8-month cell of the table is the mean of 40 samples (300 merchants each), while the 0.425 of the main sample is a single draw of 500 merchants; under the same design a small gap between the two is normal.

Block D5 compares the two bias corrections under the same design, with 1,000 samples of 500 merchants each per length.

```r
d5 <- sapply(c(8, 16, 32), function(Tobs) {                 # months per merchant; T = Tobs - 1
  TT <- Tobs - 1
  set.seed(switch(as.character(Tobs), "8" = 1151, "16" = 1152, "32" = 1153))
  m <- t(vapply(seq_len(1000), function(r) {
    dd <- gen_dyn(N = 500, T = Tobs, rho = 0.70)
    f <- fe_ar1(dd$y, dd$id, dd$t)
    c(fe = f$est[["rho"]], se = f$se[["rho"]], spj = spj_fe(dd$y, dd$id, dd$t))
  }, numeric(3)))
  hk <- hk_correct(m[, "fe"], TT); se_hk <- (TT + 1) / TT * m[, "se"]   # the correction is linear in FE
  cover <- function(est, se) mean(abs(est - 0.70) <= qnorm(0.975) * se)
  c(T_obs = Tobs, fe = mean(m[, "fe"]), hk = mean(hk), spj = mean(m[, "spj"]),
    rmse_hk = sqrt(mean((hk - 0.70)^2)), rmse_spj = sqrt(mean((m[, "spj"] - 0.70)^2)),
    cover_fe = cover(m[, "fe"], m[, "se"]), cover_hk = cover(hk, se_hk), cover_spj = cover(m[, "spj"], m[, "se"]))
})
print(as.data.frame(round(t(d5), 4)), row.names = FALSE)
#  T_obs     fe     hk    spj rmse_hk rmse_spj cover_fe cover_hk cover_spj
#      8 0.4163 0.6186 0.6754  0.0839   0.0373        0    0.022     0.625
#     16 0.5723 0.6772 0.7095  0.0252   0.0173        0    0.469     0.747
#     32 0.6410 0.6940 0.7057  0.0090   0.0100        0    0.847     0.790
```

At all three lengths, the average estimates of FE, HK, and the jackknife sit right on their own exact limits (FE 0.4162, 0.5723, and 0.6413; HK 0.6185, 0.6771, and 0.6942; jackknife 0.6747, 0.7092, and 0.7059), just as Section 3.8 said: the corrections remove only the part of order $1/T$. Coverage tells a different story: both corrections push the point estimate toward the truth, but their intervals fall far short of the nominal 95%.

### 6.4 Block D2: Instrument Proliferation

```r
sys_fit <- function(pd, L, coll, model = "onestep") {   # system GMM with month effects, instruments lags 2:L
  f <- as.formula(sprintf("y ~ lag(y, 1) | lag(y, 2:%d)", L))
  suppressWarnings(pgmm(f, pd, effect = "twoways", model = model, transformation = "ld", collapse = coll))
}
pd2 <- pdata.frame(dD2, index = c("id", "t"))
for (coll in c(FALSE, TRUE)) {
  s <- suppressWarnings(gmm_summary(sys_fit(pd2, 99, coll)))   # every available lag, uncollapsed or collapsed
  cat(sprintf("r = %3d: %.4f [%.4f, %.4f]; Hansen J %.2f on %d df, p %.4f\n", s$r, s$est, s$ci[1], s$ci[2],
              s$J, as.integer(s$J_df), s$J_p))
}
cat(sprintf("OLS in levels (month dummies) %.4f; FE %.4f\n", pooled_ar1(dD2$y, dD2$id, dD2$t, time_fe = TRUE)$est,
            fe_ar1(dD2$y, dD2$id, dD2$t)$est[["rho"]]))
# r = 134: 0.8197 [0.7600, 0.8793]; Hansen J 70.00 on 118 df, p 0.9999
# r =  30: 0.6292 [0.5296, 0.7288]; Hansen J 20.30 on 14 df, p 0.1210
# OLS in levels (month dummies) 0.9292; FE 0.5494
```

The 134 instruments consist of 105 for the differenced equation, 14 for the levels equation, and 15 month dummies; $r = 134$ exceeds the 70 merchants, and $J$ equals exactly 70 (Theorem 6 (i)). The lag-depth path in the figure of Section 3.6 comes from the same sample: the lag depth runs to 2, 3, 4, 6, and 9 and then to every lag, with each depth fitted once uncollapsed and once collapsed.

```r
for (coll in c(FALSE, TRUE)) for (L in P_D2$lags) {   # the lag-depth path drawn in 3.6
  s <- suppressWarnings(gmm_summary(sys_fit(pd2, L, coll)))
  cat(sprintf("%-9s lags 2:%-2d r = %3d: %.4f; Hansen J %6.2f on %3d df, p %.4f\n", if (coll) "collapsed" else "full",
              L, s$r, s$est, s$J, as.integer(s$J_df), s$J_p))
}
# full      lags 2:2  r =  43: 0.7311; Hansen J  34.85 on  27 df, p 0.1427
# full      lags 2:3  r =  56: 0.7312; Hansen J  54.21 on  40 df, p 0.0662
# full      lags 2:4  r =  68: 0.7478; Hansen J  63.11 on  52 df, p 0.1391
# full      lags 2:6  r =  89: 0.7690; Hansen J  70.00 on  73 df, p 0.5778
# full      lags 2:9  r = 113: 0.8032; Hansen J  70.00 on  97 df, p 0.9824
# full      lags 2:99 r = 134: 0.8197; Hansen J  70.00 on 118 df, p 0.9999
# collapsed lags 2:2  r =  17: 0.6128; Hansen J   0.62 on   1 df, p 0.4300
# collapsed lags 2:3  r =  18: 0.6355; Hansen J   2.16 on   2 df, p 0.3398
# collapsed lags 2:4  r =  19: 0.6212; Hansen J   4.77 on   3 df, p 0.1897
# collapsed lags 2:6  r =  21: 0.6196; Hansen J  10.05 on   5 df, p 0.0739
# collapsed lags 2:9  r =  24: 0.6114; Hansen J  15.91 on   8 df, p 0.0437
# collapsed lags 2:99 r =  30: 0.6292; Hansen J  20.30 on  14 df, p 0.1210
```

A single sample shows only one draw; the Monte Carlo below repeats the design 500 times, and 200 times on the twin block, to look at the size and power of the tests.

```r
Jp <- function(f) unname(suppressWarnings(sargan(f))$p.value)
Jp3 <- function(pr, model) c(full = Jp(sys_fit(pr, 99, FALSE, model)), mid = Jp(sys_fit(pr, P_D2$mid_lag, FALSE, model)),
                             coll = Jp(sys_fit(pr, 99, TRUE, model)))
rej <- function(p) sprintf("full %.3f, lags 2:4 %.3f, collapsed %.3f", mean(p[, 1] < 0.05), mean(p[, 2] < 0.05),
                           mean(p[, 3] < 0.05))
set.seed(1112)
mc <- t(vapply(seq_len(P_D2$R_mc), function(r) {
  dr <- gen_D2(P_D2); pr <- pdata.frame(dr, index = c("id", "t"))
  c(full = coef(sys_fit(pr, 99, FALSE))[[1]], collapsed = coef(sys_fit(pr, 99, TRUE))[[1]],
    ols = pooled_ar1(dr$y, dr$id, dr$t, time_fe = TRUE)$est, fe = fe_ar1(dr$y, dr$id, dr$t)$est[["rho"]],
    p1 = Jp3(pr, "onestep"), p2 = Jp3(pr, "twosteps"))   # J of the one-step fits; Hansen J of the two-step fits
}, numeric(10)))
closer <- abs(mc[, "full"] - mc[, "ols"]) < abs(mc[, "collapsed"] - mc[, "ols"])
cat(sprintf("means: full %.4f, collapsed %.4f, OLS in levels %.4f, FE %.4f; full closer to OLS in %.1f%%\n",
            mean(mc[, "full"]), mean(mc[, "collapsed"]), mean(mc[, "ols"]), mean(mc[, "fe"]), 100 * mean(closer)))
cat(sprintf("size, one-step J: %s\nsize, two-step J: %s; full p > 0.99 in %.1f%%\n", rej(mc[, 5:7]), rej(mc[, 8:10]),
            100 * mean(mc[, "p2.full"] > 0.99)))
set.seed(1132)                                   # the MA(1) twin: every Hansen test should reject
mis <- t(vapply(seq_len(P_D2_ma$R_mc), function(r) {
  pr <- pdata.frame(gen_D2_ma(P_D2_ma), index = c("id", "t"))
  c(p1 = Jp3(pr, "onestep"), p2 = Jp3(pr, "twosteps"))
}, numeric(6)))
cat(sprintf("power, one-step J: %s\npower, two-step J: %s\n", rej(mis[, 1:3]), rej(mis[, 4:6])))
# means: full 0.8089, collapsed 0.6967, OLS in levels 0.9271, FE 0.5699; full closer to OLS in 99.0%
# size, one-step J: full 0.000, lags 2:4 0.074, collapsed 0.266
# size, two-step J: full 0.000, lags 2:4 0.000, collapsed 0.056; full p > 0.99 in 100.0%
# power, one-step J: full 0.000, lags 2:4 0.235, collapsed 1.000
# power, two-step J: full 0.000, lags 2:4 0.000, collapsed 1.000
```

`sargan(f)` builds $\sum_i g_i g_i'$ from the one-step residuals, uses its (generalized) inverse as the weight, and evaluates the moments at the fit's estimate: for a one-step fit this gives the $J$ of Theorem 6, and for a two-step fit the Hansen J of Theorem 4 (iv), which is also what `xtabond2` reports after a robust one-step estimate (Roodman 2009a §2.5). The `gmm_summary` printouts in Sections 6.2, 6.4, and 6.5 label the statistic Hansen J (or Hansen p) for one-step and two-step fits alike, and the figure in Section 3.6 also labels the $J$ of the one-step fits as Hansen J. Wherever the text reports a Hansen J alongside a one-step estimate, it is the $J$ of Theorem 6. These result lines supply the numbers in the two interpretive paragraphs after Theorem 6 in Section 3.6. As the number of merchants goes to infinity, the limit of OLS in levels is 0.928 by the expression of Section 3.2, and that of FE is 0.572 by Theorem 2; the Monte Carlo averages, 0.927 and 0.570, sit right beside them. The full-instrument estimates lie between the truth and OLS, and the collapsed average, 0.697, is almost exactly the truth.

### 6.5 Blocks D3 and D4: Predetermined Regressors and Initial Conditions

```r
fp <- function(p) if (p < 1e-4) "< 1e-4" else sprintf("%.4f", p)
f3 <- fe_ar1(dD3$y, dD3$id, dD3$t, X = cbind(x = dD3$x))          # within estimator of (rho, beta)
cat(sprintf("FE                 rho %.4f, beta %.4f [%.4f, %.4f]\n", f3$est[["rho"]], f3$est[["x"]],
            ci95(f3$est[["x"]], f3$se[["x"]])[1], ci95(f3$est[["x"]], f3$se[["x"]])[2]))
pd3 <- pdata.frame(dD3, index = c("id", "t"))
sets <- c(predetermined = "1:99", endogenous = "2:99", strictly_exogenous = "0:99")   # lags of X as instruments
for (k in names(sets)) {
  g <- pgmm(as.formula(sprintf("y ~ lag(y, 1) + x | lag(y, 2:99) + lag(x, %s)", sets[[k]])), pd3,
            effect = "individual", model = "onestep", transformation = "d")
  b <- gmm_summary(g, "x"); r <- gmm_summary(g, "lag(y, 1)")
  cat(sprintf("%-18s rho %.4f, beta %.4f [%.4f, %.4f]; r = %d, Hansen p %s, m2 p %s\n", k, r$est, b$est,
              b$ci[1], b$ci[2], b$r, fp(b$J_p), fp(b$m2_p)))
}
# FE                 rho 0.3055, beta 0.3462 [0.3210, 0.3714]
# predetermined      rho 0.6097, beta 0.5128 [0.4617, 0.5638]; r = 48, Hansen p 0.6677, m2 p 0.6046
# endogenous         rho 0.5765, beta 0.4720 [0.3607, 0.5834]; r = 42, Hansen p 0.7708, m2 p 0.7923
# strictly_exogenous rho 0.5183, beta 0.8348 [0.7798, 0.8897]; r = 54, Hansen p < 1e-4, m2 p < 1e-4
```

Under the strictly exogenous treatment even $m_2$ rejects: once the estimate is biased, serial correlation is left in the residuals. Of the three instrument sets, only the first is consistent while making full use of the information.

```r
pd4 <- pdata.frame(dD4, index = c("id", "t"))
for (v in c("y", "y_stat_true")) for (tf in c("d", "ld")) {     # y: 80% start; y_stat_true: mean-stationary twin
  f <- as.formula(sprintf("%s ~ lag(%s, 1) | lag(%s, 2:99)", v, v, v))
  s <- gmm_summary(pgmm(f, pd4, effect = "twoways", model = "onestep", transformation = tf), sprintf("lag(%s, 1)", v))
  cat(sprintf("%-11s %-2s %.4f [%.4f, %.4f]; r = %d, Hansen p %s, m2 p %.3f\n", v, if (tf == "d") "AB" else "BB",
              s$est, s$ci[1], s$ci[2], s$r, fp(s$J_p), s$m2_p))
}
# y           AB 0.5927 [0.5347, 0.6508]; r = 27, Hansen p 0.0515, m2 p 0.646
# y           BB 0.8160 [0.7983, 0.8336]; r = 34, Hansen p < 1e-4, m2 p 0.227
# y_stat_true AB 0.6009 [0.5552, 0.6467]; r = 27, Hansen p 0.0717, m2 p 0.615
# y_stat_true BB 0.6178 [0.5752, 0.6604]; r = 34, Hansen p 0.2306, m2 p 0.585
```

The two `y` rows and the two `y_stat_true` rows use the same random draws; the only difference is where merchants start in their first month. Difference GMM is unaffected by this and aims at 0.6 both times; system GMM drifts to 0.816 with the nonstationary start and comes back to 0.618 with the stationary one.

### 6.6 The Grand Reconciliation

| Block | Estimator | Estimate | 95% interval | Target and true value | Covered |
|---|---|---|---|---|---|
| D1 | Pooled OLS | 0.946 | [0.938, 0.955] | $\rho$ 0.70 (probability limit 0.955) | No |
| D1 | FE | 0.425 | [0.388, 0.461] | $\rho$ 0.70 (probability limit 0.416) | No |
| D1 | FE + Hahn-Kuersteiner | 0.628 | [0.586, 0.670] | $\rho$ 0.70 (probability limit 0.6185) | No |
| D1 | Anderson-Hsiao, $Y_{i,t-2}$ | 0.680 | [0.559, 0.802] | $\rho$ 0.70 | Yes |
| D1 | One-step difference GMM | 0.707 | [0.612, 0.802] | $\rho$ 0.70 | Yes |
| D1 | Two-step difference GMM | 0.7045 | [0.597, 0.812] | $\rho$ 0.70 | Yes |
| D1 | Two-step system GMM | 0.701 | [0.633, 0.768] | $\rho$ 0.70 | Yes |
| D1p | FE | 0.535 | [0.504, 0.567] | $\rho$ 0.90 (probability limit 0.558) | No |
| D1p | Anderson-Hsiao, $\Delta Y_{i,t-2}$ | 0.312 | [$-0.201$, 0.824] | $\rho$ 0.90 | No |
| D1p | Two-step difference GMM | 0.809 | [0.542, 1.077] | $\rho$ 0.90 | Yes |
| D1p | Two-step system GMM | 0.911 | [0.862, 0.960] | $\rho$ 0.90 | Yes |
| D2 | One-step system GMM, 134 instruments | 0.820 | [0.760, 0.879] | $\rho$ 0.70 | No |
| D2 | One-step system GMM, collapsed to 30 instruments | 0.629 | [0.530, 0.729] | $\rho$ 0.70 | Yes |
| D3 | FE, $\beta$ | 0.346 | [0.321, 0.371] | $\beta$ 0.5 (large-sample estimate 0.367) | No |
| D3 | Difference GMM, promotion as predetermined | 0.513 | [0.462, 0.564] | $\beta$ 0.5 | Yes |
| D3 | Difference GMM, promotion as strictly exogenous | 0.835 | [0.780, 0.890] | $\beta$ 0.5 | No |
| D4 | One-step difference GMM | 0.593 | [0.535, 0.651] | $\rho$ 0.6 | Yes |
| D4 | One-step system GMM | 0.816 | [0.798, 0.834] | $\rho$ 0.6 (large-sample one-step estimate 0.817) | No |
| D4 twin | One-step system GMM | 0.618 | [0.575, 0.660] | $\rho$ 0.6 | Yes |

Almost every row of the table that fails to cover the truth does so by design: pooled OLS, FE, and HK aim at their own probability limits; the full instrument set of D2 is pulled toward OLS in levels; the FE of D3 and the GMM that treats promotion as strictly exogenous use moments that do not hold; and the system GMM of D4 bets on the wrong initial conditions. These biases come from the identifying assumptions, not from luck. The only exception is Anderson-Hsiao with $\Delta Y_{i,t-2}$ as the instrument on D1p, which should aim at 0.90 but this time misses it. A nominal 95% interval is bound to miss now and then, and across 1,000 repetitions of the same design its coverage is 0.933, so this one miss is no surprise; the real problem is that the instrument is weak and the estimate far too imprecise, with a standard error of 0.261.

The reconciliation of this section can be summed up as follows: the FE that recovered the true promotion effect in Chapter 10 is pressed down to 0.425 by Nickell bias in the dynamic panel, and together with pooled OLS, inflated to 0.946, it straddles the true persistence; GMM with differencing plus lagged levels recovers $\rho$ at 0.70, at a cost in precision that is highest when $\rho$ is close to 1; system GMM restores the precision, provided the initial conditions are right; and with many instruments the estimate is pulled back toward the biased side while the Hansen test still waves it through.

## 7 Failure Modes and Robustness

The problems in the simulation are manufactured, and in real panels they are better hidden. Below we lay out the most common failure modes and the responses to them.

Nickell bias lurks in any short panel with a lagged dependent variable, and no number of merchants can help. In any setup with a lagged outcome on the right-hand side and a $T$ that is not large, the FE estimate of persistence cannot be trusted directly. There is a hidden variant. Even if what you care about is not persistence but the effect of some covariate, as long as a lagged dependent variable is included in the model to control for dynamics, the bias will seep along the correlation into that covariate's coefficient: in D3, even with feedback switched off, the promotion effect estimated by FE in the large sample of a hundred thousand merchants is only 0.466. "Adding a lagged dependent variable to control for dynamics a bit" is not a harmless move.

When the dynamics are misspecified, the lagged instruments fail. If the error is serially correlated, say because shocks persist, $Y_{i,t-2}$ is no longer clean. $m_2$ is the gatekeeper, but its power is limited: in Arellano and Bond's own simulations, the test rejects only half the time when the autocorrelation is weak, and Chung, Kim and Park stress this point as well. The response is to make the dynamics more flexible (add one more lag) and to check whether the estimate is stable when the instruments start from more distant lags.

When persistence approaches 1, the instruments of difference GMM become weak, the estimate is imprecise, and the weak-instrument bias also points in the direction of FE. System GMM can restore precision, provided that merchants' distance from their long-run level when they enter the sample is unrelated to talent. This is the condition platform data are most likely to violate: newly onboarded merchants, newly launched products, and newly opened cities are all ramping up, start systematically below their long-run levels, and the higher the talent, the farther they start from the target, as in D4. In such samples system GMM can be worse than difference GMM. The basic robustness exercise is to report both; when they differ markedly, run a Difference-in-Hansen test on the instruments of the levels equation, or drop the earliest months and estimate again.

Instrument proliferation is the most common pitfall. The number of instruments grows quadratically with $T$, and too many of them bias GMM toward the estimator that uses no instruments (difference GMM toward FE, system GMM toward OLS in levels) while robbing the Hansen test of power; when the instruments outnumber the units, the $J$ of a one-step fit is identically the number of units and the two-step $J$ cannot exceed it. A model with hundreds of instruments and a $p$-value close to 1 is not grounds for reassurance but an alarm of overfitting. The response is to collapse, limit the lag depth, report the number of instruments, and show how sensitive the estimates and tests are to the number of instruments.

The classification of regressors decides whether the instruments are valid. Treating a regressor that responds to past shocks as strictly exogenous uses invalid moments; treating it as endogenous throws information away. The key question is: could the value of this variable be influenced by past outcomes? On a platform, promotion budgets, prices, inventory, and advertising almost always are.

Bias-corrected FE is a useful complement when $T$ is moderate, but in short panels it removes only part of the bias and its intervals cannot be trusted; the corrections rely on the AR(1) bias structure or on the stationarity of the panel, and call for caution when there are other regressors or the start is nonstationary.

Finally, inference. Always use Windmeijer standard errors for two-step estimates; $m_l$ and $J$ both assume independence across merchants, and time effects can absorb shocks common to all merchants; when merchants are few, these asymptotic tests are themselves unreliable.

Put these failure modes together, and the robustness of dynamic panels comes down to three judgments: whether the dynamics are fully specified, that is, whether the error is still serially correlated; what exogeneity the regressors and the initial conditions satisfy, which decides which moments are valid; and how many instruments are used, which decides how far the estimate is from the biased estimator that uses no instruments and whether the tests still have power. Once these three judgments are thought through clearly, the panel's gift of differencing out heterogeneity is truly in hand; if they are not, demeaning and differencing, the two knives that scrub out contamination for you, will also quietly wound the dynamic parameter.

## 8 Further Reading

::: {.readings}
Required reading, in suggested reading order:

- Wooldridge (2010, *Econometric Analysis of Cross Section and Panel Data*, 2nd ed.), §11.6. Sequential exogeneity, models with a lagged dependent variable, and instruments for the differenced equation; the source of the equation numbers in Theorems 1, 3, and 7 of this chapter.
- Arellano and Bond (1991, *Review of Economic Studies*). The founding work of difference GMM and the serial correlation tests, the starting point of empirical work with dynamic panels.
- Blundell and Bond (1998, *Journal of Econometrics*). System GMM and the initial conditions: why it is needed when persistence approaches 1, and when it cannot be used; read together with Arellano and Bover (1995).
- Roodman (2009b, *Oxford Bulletin of Economics and Statistics*) and Roodman (2009a, *Stata Journal*). The former covers instrument proliferation and collapsing; the latter is the most practical operating guide to difference GMM and system GMM, including the derivations of the Windmeijer correction and the serial correlation tests.

Further reading:

- Nickell (1981, *Econometrica*). The original derivation of the dynamic FE bias, the source of Theorem 2 of this chapter.
- Anderson and Hsiao (1981, *Journal of the American Statistical Association*; 1982, *Journal of Econometrics*). The original idea of using lagged values as instruments.
- Windmeijer (2005, *Journal of Econometrics*). The finite-sample correction of two-step GMM standard errors.
- Bun and Windmeijer (2010, *Econometrics Journal*). The levels equation of system GMM can also suffer from weak instruments.
- Bowsher (2002, *Economics Letters*). Simulation evidence that overidentification tests in dynamic panels lose power when there are many instruments.
- Hahn and Kuersteiner (2002, *Econometrica*), Dhaene and Jochmans (2015, *Review of Economic Studies*), and Kiviet (1995, *Journal of Econometrics*). Analytical and jackknife bias corrections for FE.
- Chung, Kim and Park (2019, *Management Science*). More robust moment conditions under predetermined regressors and serial correlation, with a marketing application.
- Moral-Benito, Allison and Williams (2019, *Applied Economics*). Maximum likelihood for dynamic panels, an alternative to Arellano-Bond.
- Hansen (2014) §19.3 (p. 329). A concise treatment of dynamic panel regression in less than a page.
:::

::: {.apa-refs}
- Anderson, T. W., & Hsiao, C. (1981). Estimation of dynamic models with error components. *Journal of the American Statistical Association, 76*(375), 598-606. https://doi.org/10.1080/01621459.1981.10477691
- Anderson, T. W., & Hsiao, C. (1982). Formulation and estimation of dynamic models using panel data. *Journal of Econometrics, 18*(1), 47-82. https://doi.org/10.1016/0304-4076(82)90095-1
- Arellano, M., & Bond, S. (1991). Some tests of specification for panel data: Monte Carlo evidence and an application to employment equations. *The Review of Economic Studies, 58*(2), 277-297. https://doi.org/10.2307/2297968
- Arellano, M., & Bover, O. (1995). Another look at the instrumental variable estimation of error-components models. *Journal of Econometrics, 68*(1), 29-51. https://doi.org/10.1016/0304-4076(94)01642-D
- Blundell, R., & Bond, S. (1998). Initial conditions and moment restrictions in dynamic panel data models. *Journal of Econometrics, 87*(1), 115-143. https://doi.org/10.1016/S0304-4076(98)00009-8
- Bowsher, C. G. (2002). On testing overidentifying restrictions in dynamic panel data models. *Economics Letters, 77*(2), 211-220. https://doi.org/10.1016/S0165-1765(02)00130-1
- Bun, M. J. G., & Windmeijer, F. (2010). The weak instrument problem of the system GMM estimator in dynamic panel data models. *The Econometrics Journal, 13*(1), 95-126. https://doi.org/10.1111/j.1368-423X.2009.00299.x
- Chung, D. J., Kim, B., & Park, B. G. (2019). How do sales efforts pay off? Dynamic panel data analysis in the Nerlove-Arrow framework. *Management Science, 65*(11), 5197-5218. https://doi.org/10.1287/mnsc.2018.3189
- Dhaene, G., & Jochmans, K. (2015). Split-panel jackknife estimation of fixed-effect models. *The Review of Economic Studies, 82*(3), 991-1030. https://doi.org/10.1093/restud/rdv007
- Forbes, K. J. (2000). A reassessment of the relationship between inequality and growth. *American Economic Review, 90*(4), 869-887. https://doi.org/10.1257/aer.90.4.869
- Hahn, J., & Kuersteiner, G. (2002). Asymptotically unbiased inference for a dynamic panel model with fixed effects when both n and T are large. *Econometrica, 70*(4), 1639-1657. https://doi.org/10.1111/1468-0262.00344
- Hansen, B. E. (2014). *Econometrics* [Unpublished manuscript]. Department of Economics, University of Wisconsin.
- Kiviet, J. F. (1995). On bias, inconsistency, and efficiency of various estimators in dynamic panel data models. *Journal of Econometrics, 68*(1), 53-78. https://doi.org/10.1016/0304-4076(94)01643-E
- Levine, R., Loayza, N., & Beck, T. (2000). Financial intermediation and growth: Causality and causes. *Journal of Monetary Economics, 46*(1), 31-77. https://doi.org/10.1016/S0304-3932(00)00017-9
- Moral-Benito, E., Allison, P., & Williams, R. (2019). Dynamic panel data modelling using maximum likelihood: An alternative to Arellano-Bond. *Applied Economics, 51*(20), 2221-2232. https://doi.org/10.1080/00036846.2018.1540854
- Nickell, S. (1981). Biases in dynamic models with fixed effects. *Econometrica, 49*(6), 1417-1426. https://doi.org/10.2307/1911408
- Roodman, D. (2009a). How to do xtabond2: An introduction to difference and system GMM in Stata. *The Stata Journal, 9*(1), 86-136. https://doi.org/10.1177/1536867X0900900106
- Roodman, D. (2009b). A note on the theme of too many instruments. *Oxford Bulletin of Economics and Statistics, 71*(1), 135-158. https://doi.org/10.1111/j.1468-0084.2008.00542.x
- Windmeijer, F. (2005). A finite sample correction for the variance of linear efficient two-step GMM estimators. *Journal of Econometrics, 126*(1), 25-51. https://doi.org/10.1016/j.jeconom.2004.02.005
- Wooldridge, J. M. (2010). *Econometric analysis of cross section and panel data* (2nd ed.). MIT Press.
:::
