---
title: "Panel Data & Dynamic Panels"
subtitle: "From Fixed Effects to Arellano-Bond"
seriesline: "Foundations of Information Systems Economics · Chapter 7"
series: "Foundations of Information Systems Economics: Empirical Methods for Information Systems Economics"
sidebartitle: "Foundations of Information Systems Economics"
sidebarsub: "Theory · Methods · Applications"
chapterlabel: "Chapter 7 · Panel Data & Dynamic Panels"
---

## Introduction

Nimbus wants to know how long a merchant's good month will last. Regress this month's sales directly on last month's, and persistence comes out at 0.946; add merchant fixed effects, and the estimate crashes to 0.425. Both models are perfectly standard, yet the two answers differ by a factor of two, while the true value known from the simulation, 0.70, sits squarely between them. This is not statistics being coy: the first model mistakes merchant talent for persistence, and the second, in scrubbing that talent away, manufactures a new correlation of its own.

The appeal of panel data lies precisely in its ability to follow the same unit over time and so difference out unobservable differences that do not vary over time. For static models, fixed effects are often a reliable workhorse; but once a lagged outcome appears on the right-hand side, the demeaned lag becomes correlated with the demeaned error, and Nickell bias shows up in short panels. The distinctions among strict exogeneity, sequential exogeneity, and the random effects assumption determine which transformation and which estimator are entitled to work.

This chapter starts from the split between 0.946 and 0.425, explains how Arellano-Bond builds instruments for the differenced equation out of earlier levels, and why Blundell-Bond adds a levels equation when persistence is high. Dynamic panel GMM is not a push-button remedy: the instruments can be too weak, or so numerous that they "explain" everything. Serial correlation tests and overidentification tests are therefore part of the analysis. Panel data really do hand the researcher a second dimension, but time never gives away identification for free.

## 1 Two Wrong Numbers That Straddle the Truth

::: {.case}
Nimbus is a fictional e-commerce platform. We tell this case with simulated data, and the payoff is that the data generating process is fully known, so every estimator can be reconciled against the truth, a teaching condition that real data can never provide. The platform observes 500 merchants across 8 months of monthly GMV (here a standardized log-GMV measure). Management wants to know something that matters operationally: how "sticky" a merchant's GMV is, that is, to what extent a good month carries over into the next. We capture this with a persistence parameter $\rho$; the closer $\rho$ is to 1, the more durable the effect of a shock, and the more worthwhile it becomes to take a long view of the platform's subsidies and campaigns.

The analyst runs the most direct regression, this month's GMV on last month's (pooled OLS). The estimated persistence is $\hat\rho = 0.946$, a frighteningly high number, implying that GMV is almost as sticky as a random walk, with last month's good showing carrying over almost entirely into this month. She grows wary at once: merchants differ innately, top merchants post high GMV month after month, and this merchant-specific operating talent might have inflated persistence. So she adds merchant fixed effects to difference out each merchant's fixed level and estimates again. This time $\hat\rho = 0.425$, nearly cut in half, implying that GMV is not so sticky after all and that shocks fade fairly fast.
:::

Two numbers, one saying GMV is extremely sticky (0.946), the other saying not so much (0.425), differing by more than a factor of two. Which do we trust? Because this is a simulation, we know the true value is $\rho = 0.70$. The answer is unexpected: both are wrong, and the truth happens to be caught between them.

This is no coincidence but two systematic biases pushing the estimate in opposite directions. Pooled OLS runs high because it never removes merchant talent: a high-talent merchant posts high GMV both this month and last, the talent lifts both the dependent variable and its lag at once, and OLS misreads this "talent-induced high correlation" as "persistence," pushing $\hat\rho$ up to 0.946. Fixed effects run low for a far subtler reason, which is the heart of this chapter: the demeaning step that differences out merchant talent creates a negative correlation between the transformed lag and the transformed error, systematically pressing persistence down, down to 0.425. This downward bias is called Nickell bias, and Section 3 will derive it cleanly.

::: {.intuition}
An honest estimate of persistence must both remove merchant talent and avoid introducing Nickell bias, and pooled OLS and fixed effects each manage one and botch the other. Pooled OLS leaves talent untouched and is contaminated by it; fixed effects remove talent but damage the dynamic structure in the removing. That the truth 0.70 is straddled by 0.946 and 0.425 is not luck but the inevitable consequence of these two biases pointing in opposite directions. Later in the chapter we will see that Arellano-Bond dynamic panel GMM sidesteps both pitfalls at once with a clever device, pulling $\hat\rho$ back to around 0.70. To understand why it is clever, we first have to understand why fixed effects are a hero in a static panel and why they turn traitor in a dynamic one.
:::

The lesson of Section 1 can be summed up as follows: in panel data, pooled OLS that ignores individual heterogeneity overstates persistence, while fixed effects that remove heterogeneity understate it through Nickell bias, and the truth is caught between them; to estimate the dynamic parameter correctly we can neither ignore heterogeneity nor use the fixed effects that introduce Nickell bias, which leads us to a GMM designed specifically for dynamic panels.

## 2 The Economic Model and the Estimand

Before rolling up our sleeves we set up the panel model and, in the simplest static case, install fixed effects as the hero, because all of the trouble in dynamic panels is the story of this hero turning traitor.

### 2.1 The Panel Model and Individual Heterogeneity

The basic panel model writes the outcome as the sum of observable covariates, an individual effect, and an idiosyncratic error:

$$y_{it} = x_{it}'\beta + \alpha_i + \varepsilon_{it}.$$

The subscript $i$ is the unit (a merchant in Nimbus) and $t$ is the time point (a month). $x_{it}$ is a time-varying covariate and $\beta$ is the coefficient we want. The key new character is $\alpha_i$, the individual effect, which captures everything time-invariant and not in $x$: a merchant's operating talent, the natural scale of its category, the capability of its founding team. $\varepsilon_{it}$ is an idiosyncratic shock that fluctuates over time. $\alpha_i$ is unobservable, and it is almost surely correlated with $x_{it}$, which is exactly the core difficulty of identification: a high-talent merchant has both high GMV (high $\alpha_i$) and a greater willingness to spend on promotions (high $x_{it}$), so stuffing $\alpha_i$ into the error term and running OLS directly contaminates $\beta$ with this correlation.

The way out that panels offer is to use the fact that the same merchant is observed repeatedly, and difference $\alpha_i$ away. The most common device is the within transformation, also called demeaning: for each unit, subtract its own time average across periods. Because $\alpha_i$ does not vary over time, it equals its own time mean, and one subtraction removes it:

$$y_{it} - \bar y_i = (x_{it} - \bar x_i)'\beta + (\varepsilon_{it} - \bar\varepsilon_i).$$

Running OLS on the demeaned data is the fixed effects (FE) estimator, also called the within estimator. $\alpha_i$ is wiped out completely, and however it correlates with $x$ it no longer does any harm. This is the precise cashing-in of the panel's gift: as long as the unobservable heterogeneity does not vary over time, fixed effects are immune to it, even if it correlates arbitrarily with the covariates.

### 2.2 Fixed Effects as Hero: A Static Example

Let us stage this gift on Nimbus. Suppose promotion intensity $x_{it}$ has a true effect $\beta = 0.5$ on GMV, and merchant talent $\alpha_i$ is positively correlated with promotion (high-talent merchants promote more). Here three estimators give three different numbers. Pooled OLS leaves $\alpha_i$ in the error and yields $\hat\beta = 0.951$, badly inflated by the positive correlation between talent and promotion. Fixed effects difference $\alpha_i$ away and yield $\hat\beta = 0.506$, hugging the true value of 0.5. In between sits a random effects (RE) estimator that gives 0.739; Section 3.2 explains what it is and why it too is biased. Fixed effects are the only winner here, and they deliver on the panel's promise.

::: {.intuition}
Why do fixed effects win? In a phrase: they only compare each merchant with itself. They do not ask "do merchants that promote more have higher GMV" (a comparison contaminated by talent, because merchants that promote more tend to be more talented) but rather "in the months when a given merchant promotes more, is its GMV higher than in its own months of less promotion." A single merchant's talent is the same constant across different months, and a within comparison cancels it automatically, leaving the net effect of promotion. This "compare with itself" idea is the entire source of the power of fixed effects, and also their entire cost: they use only within (over-time) variation and throw away all between (across-unit) variation, so they cannot estimate any time-invariant covariate and are not precise when sample variation is scarce.
:::

### 2.3 The Dynamic Estimand

Fixed effects triumph in the static case, but management's real question is dynamic: how sticky is GMV. This requires putting the outcome's own lag on the right-hand side:

$$y_{it} = \rho\, y_{i,t-1} + \alpha_i + \varepsilon_{it}.$$

The target estimand is the persistence $\rho$. On the surface, this just treats $y_{i,t-1}$ as one more covariate, so copying the fixed effects demeaning from Section 2.1 ought to work. Section 1 has already spoiled the ending: it does not, because the demeaning step is exactly where things go wrong here, and $\hat\rho$ is pressed down to 0.425. Why the fixed effects that were invincible in the static case fail the moment a lagged dependent variable appears is the heart of the next section's identification analysis.

The main points of Section 2 can be summed up as follows: the panel model splits the outcome into covariates, an individual effect $\alpha_i$, and an idiosyncratic error, and the unobservable heterogeneity $\alpha_i$ is the primary threat to identification; fixed effects difference $\alpha_i$ away by within demeaning and identify $\beta$ purely by "comparing with itself," a hero in the static panel; but once the target estimand becomes the dynamic persistence $\rho$ and a lagged dependent variable appears on the right-hand side, this hero turns traitor.

## 3 Identification: The Fixed-Effects Gift and the Nickell Curse

This is the most analytically detailed section of the chapter. It first clarifies the exogeneity assumption on which fixed effects rest (3.1) and the choice between fixed and random effects together with the Hausman test (3.2); then the distinction between strict and sequential exogeneity and why a lagged dependent variable breaks the former (3.3); then it derives Nickell bias and exposes the failure of fixed effects in the dynamic model (3.4); and finally it lays out the remedy: differencing plus lagged levels as instruments (3.5).

### 3.1 What Fixed Effects Need: Strict Exogeneity

For fixed effects to estimate $\beta$ consistently, it takes more than the mechanical step of "$\alpha_i$ being differenced away"; it also takes an assumption about the error. What remains after demeaning is $(\varepsilon_{it} - \bar\varepsilon_i)$, and for it to be uncorrelated with the demeaned covariate $(x_{it} - \bar x_i)$ requires a condition far stronger than cross-sectional exogeneity.

::: {.assumption}
**(strict exogeneity)** The error is uncorrelated with the covariates in every period, conditional on the individual effect:

$$\mathbb{E}[\varepsilon_{it} \mid x_{i1}, x_{i2}, \ldots, x_{iT}, \alpha_i] = 0, \quad t = 1, \ldots, T.$$
:::

Note that the conditioning set here is $x$ from every period, from $x_{i1}$ all the way to $x_{iT}$, not just the contemporaneous $x_{it}$. It requires the contemporaneous error $\varepsilon_{it}$ to be unrelated to covariates in all periods past, present, and future. Why so strong? Because within demeaning stirs all periods together: $\bar x_i$ is the average of $x$ across all periods, $\bar\varepsilon_i$ is the average of $\varepsilon$ across all periods, the demeaned $(x_{it} - \bar x_i)$ contains $x$ from every period, and the demeaned $(\varepsilon_{it} - \bar\varepsilon_i)$ contains $\varepsilon$ from every period. As soon as $x$ in any period correlates with $\varepsilon$ in any period, the two demeaned quantities correlate and fixed effects are biased. Strict exogeneity is precisely what seals off all these cross-period correlations at once.

It is worth spelling out what strict exogeneity rules out, because that leads straight to the trouble in dynamic panels. The most typical case it rules out is a covariate that responds to past shocks (feedback): if this month's promotion $x_{it}$ was partly raised because last month's GMV came in surprisingly low ($\varepsilon_{i,t-1}$ negative), then $x_{it}$ correlates with $\varepsilon_{i,t-1}$ and strict exogeneity breaks. Worse still, a lagged dependent variable $y_{i,t-1}$ is by its very nature such a covariate that "responds to past shocks," which is the cue for Section 3.3.

### 3.2 Fixed or Random Effects: The Hausman Test

Fixed effects allow $\alpha_i$ to correlate arbitrarily with $x$, at the cost of throwing away between variation, lower precision, and an inability to estimate time-invariant covariates. Is there a way not to pay this cost? Yes, provided you dare make a stronger assumption: that $\alpha_i$ is in fact uncorrelated with $x$.

::: {.assumption}
**(random effects assumption)** The individual effect is uncorrelated with the covariates in every period:

$$\mathbb{E}[\alpha_i \mid x_{i1}, \ldots, x_{iT}] = 0.$$
:::

If this assumption holds, $\alpha_i$ is just a harmless random disturbance that has slipped into the error term and is no longer a threat to identification. We can then use the random effects estimator, which does not do full within demeaning but a partial demeaning (quasi-demeaning), subtracting $\theta_{RE}\, \bar y_i$ rather than the whole $\bar y_i$, where $\theta_{RE} = 1 - \sqrt{\sigma_\varepsilon^2 / (\sigma_\varepsilon^2 + T\sigma_\alpha^2)}$ is determined by the variance components (here $\theta_{RE}$ is a local notation for the weight in the RE transformation, unrelated to the $\theta$ the book uses for structural parameters). This partial demeaning is a GLS that uses both within and between variation, and so is more efficient than fixed effects, which use only within variation. The cost is that if the assumption that "$\alpha_i$ is unrelated to $x$" is wrong, random effects have not scrubbed the contamination clean and the estimate is biased. In the static Nimbus example RE gives 0.739, wedged between the contaminated OLS (0.951) and the clean FE (0.506): its partial demeaning shaved off some of the talent contamination, but because $\alpha_i$ really is correlated with promotion, it did not shave it clean.

How do we choose between the robustness of fixed effects and the efficiency of random effects? The Hausman test provides a data-driven criterion.

::: {.theorem}
**(Hausman test)** When the random effects assumption holds, both fixed and random effects are consistent, but random effects are more efficient, and the two should be close; when the assumption fails, only fixed effects are consistent, and the two diverge significantly. The test statistic

$$H = (\hat\beta_{FE} - \hat\beta_{RE})'\, \big[\mathrm{Var}(\hat\beta_{FE}) - \mathrm{Var}(\hat\beta_{RE})\big]^{-1}\, (\hat\beta_{FE} - \hat\beta_{RE}) \xrightarrow{d} \chi^2_k,$$

with degrees of freedom $k$ equal to the number of time-varying covariates. A large $H$ rejects the random effects assumption and calls for fixed effects.
:::

Its logic is elegant: it works with the difference between the two estimators, $\hat\beta_{FE} - \hat\beta_{RE}$. Under the null ($\alpha_i$ uncorrelated with $x$) both are correct, the difference is just sampling noise, and $H$ is small; when the null is false, random effects are biased and fixed effects are not, the difference contains a systematic bias, and $H$ is large. That the variance is taken as the difference of the two, $\mathrm{Var}(\hat\beta_{FE}) - \mathrm{Var}(\hat\beta_{RE})$, rests on a Hausman lemma: under the null, the sampling error of the efficient estimator (RE) is uncorrelated with "the difference of the two," so the variance of the difference equals exactly the difference of the two variances. In the static Nimbus example $H$ reaches 3362, the $p$-value is effectively zero, random effects are decisively rejected, and the correct choice is fixed effects, consistent with the setup in which $\alpha_i$ really is correlated with promotion.

### 3.3 Strict or Sequential: The Original Sin of the Lagged Dependent Variable

Now put the lagged dependent variable back and watch how it breaks the strict exogeneity of Section 3.1. The dynamic model is $y_{it} = \rho\, y_{i,t-1} + \alpha_i + \varepsilon_{it}$. Ask a simple question: does the lagged dependent variable $y_{i,t-1}$ satisfy strict exogeneity?

It does not, and it cannot in principle. Strict exogeneity requires $y_{i,t-1}$ to be unrelated to the errors in every period, yet $y_{i,t-1}$ by definition equals $\rho\, y_{i,t-2} + \alpha_i + \varepsilon_{i,t-1}$, and it plainly contains $\varepsilon_{i,t-1}$. So $y_{i,t-1}$ is perfectly correlated with last period's error $\varepsilon_{i,t-1}$, and strict exogeneity breaks at the root. This is not a matter of poor specification but the original sin of a lagged dependent variable: a variable that is the outcome's own past value must carry the past error.

Here we need a weaker notion of exogeneity, better suited to the dynamic model.

::: {.assumption}
**(sequential exogeneity / predetermined)** The error is uncorrelated with contemporaneous and past covariates, but may correlate with future covariates:

$$\mathbb{E}[\varepsilon_{it} \mid x_{i1}, \ldots, x_{it}, \alpha_i] = 0.$$
:::

Sequential exogeneity only requires the contemporaneous error to be unrelated to $x$ in the current and earlier periods, not the future. It is weaker than strict exogeneity: it allows covariates to respond to past shocks (future $x$ may correlate with the current error). A lagged dependent variable satisfies sequential exogeneity ($y_{i,t-1}$ is determined by information at $t-1$ and earlier, and is unrelated to the current shock $\varepsilon_{it}$) but not strict exogeneity (it correlates with $\varepsilon_{i,t-1}$). This distinction is the watershed: the consistency of fixed effects requires strict exogeneity, while a dynamic model has only sequential exogeneity, so fixed effects must run into trouble in a dynamic model. Exactly how large the trouble is, the next section quantifies.

### 3.4 Nickell Bias: How the Hero Turns Traitor

The failure of fixed effects in the dynamic model has a precise name and a precise magnitude, due to Nickell (1981). First we explain where the bias comes from, then give the formula.

::: {.intuition}
The source of the bias is a hidden consequence of the demeaning step. Within demeaning of the dynamic model gives $(y_{it} - \bar y_i) = \rho\,(y_{i,t-1} - \bar y_i) + (\varepsilon_{it} - \bar\varepsilon_i)$. Look at the demeaned lag $(y_{i,t-1} - \bar y_i)$ and the demeaned error $(\varepsilon_{it} - \bar\varepsilon_i)$; they ought to be uncorrelated for fixed effects to be unbiased. But they are correlated, and the culprit is those two means. $\bar y_i$ is the average of $y$ across all periods for this merchant, and it contains $y_{it}$, which contains $\varepsilon_{it}$, so $\varepsilon_{it}$ has leaked into $\bar y_i$. Then the demeaned lag $(y_{i,t-1} - \bar y_i)$, through the $-\bar y_i$ term, indirectly contains a $-\varepsilon_{it}/T$ component, which is negatively correlated with the $\varepsilon_{it}$ component in the demeaned error. A regressor negatively correlated with the error drags the coefficient down, and this is the downward direction of Nickell bias. The key is that the magnitude of this contamination is $1/T$: the mean is an average of $T$ numbers, and a single $\varepsilon_{it}$ makes up only a $1/T$ share of $\bar y_i$, so the larger $T$ is, the more the contamination is diluted and the smaller the bias.
:::

Nickell (1981) worked out the limit of this bias (as $N \to \infty$ with $T$ fixed), whose leading term is

$$\mathrm{plim}_{N \to \infty}\,(\hat\rho_{FE} - \rho) \approx -\,\frac{1 + \rho}{T - 1}.$$

Three things are worth reading off this formula. First, the bias is negative, fixed effects systematically understate persistence, consistent with intuition. Second, the bias is $O(1/T)$, smaller the larger $T$ is, vanishing as $T \to \infty$, so Nickell bias is a short-panel problem and irrelevant in long panels. Third, the bias does not vanish with $N$: no number of merchants can rescue a short panel, because this is mean contamination inside each unit, and adding units without adding periods is of no help. This point is especially counterintuitive and especially important, for it means that in a typical platform panel with large $N$ and small $T$ (thousands of merchants, a mere handful of months), fixed effects can be wildly wrong about persistence, and the enormous sample size helps not at all. In Nimbus $T = 8$ and $\rho = 0.70$, and the leading term predicts a bias of about $-(1.7)/7 = -0.24$, so fixed effects should press $\rho$ down to around 0.46; the actual estimate is 0.425, right in direction and magnitude (the actual bias is slightly larger than the leading term, because the full Nickell formula has higher-order negative terms that are not negligible at short $T$).

::: {.warning}
The most dangerous thing about Nickell bias is that it is silent. The fixed effects regression runs just fine, coefficients and standard errors all present, with no error message of any kind. An unwitting analyst would take the 0.425 and report "GMV persistence is moderate," when the truth is 0.70, persistence is quite high, and the two conclusions give opposite advice about whether subsidies should be invested for the long term. More insidiously, adding sample size (more merchants) not only fails to help but makes that wrong 0.425 look ever more "precise," the standard errors ever smaller, dressing up a biased estimate as ever more credible. This is the echo, in dynamic panels, of that line at the end of Chapter 5: a beautiful, precise estimate may be nothing more than a systematic bias computed very stably.
:::

### 3.5 The Remedy: Differencing Plus Lagged Levels

The root of Nickell bias is that demeaning mixes $\varepsilon_{it}$ into the lag. The first step of the remedy is to switch to a different way of removing $\alpha_i$, using first differencing in place of demeaning:

$$\Delta y_{it} = \rho\, \Delta y_{i,t-1} + \Delta \varepsilon_{it},$$

where $\Delta y_{it} = y_{it} - y_{i,t-1}$. Differencing likewise removes the time-invariant $\alpha_i$. But it does not solve the endogeneity: the differenced lag $\Delta y_{i,t-1} = y_{i,t-1} - y_{i,t-2}$ contains $\varepsilon_{i,t-1}$, and the differenced error $\Delta\varepsilon_{it} = \varepsilon_{it} - \varepsilon_{i,t-1}$ also contains $\varepsilon_{i,t-1}$, so the two are still correlated, and running OLS on the differenced equation directly is still biased.

The key move is to find an instrument for this endogenous $\Delta y_{i,t-1}$. What variable is both correlated with $\Delta y_{i,t-1}$ and uncorrelated with $\Delta\varepsilon_{it}$? The answer hides in the model's own history: sufficiently distant lagged levels. Take $y_{i,t-2}$; it is clearly correlated with $\Delta y_{i,t-1} = y_{i,t-1} - y_{i,t-2}$ ($y_{i,t-2}$ is right there inside it). What about its correlation with $\Delta\varepsilon_{it} = \varepsilon_{it} - \varepsilon_{i,t-1}$? $y_{i,t-2}$ is determined by information at $t-2$ and earlier, containing $\varepsilon_{i,t-2}$ and earlier errors, while $\Delta\varepsilon_{it}$ contains only $\varepsilon_{it}$ and $\varepsilon_{i,t-1}$, and when the errors have no serial correlation the two are uncorrelated. So $y_{i,t-2}$ is a valid instrument, giving the moment condition

$$\mathbb{E}[y_{i,t-2}\, \Delta\varepsilon_{it}] = 0.$$

This is the idea of Anderson and Hsiao, and the starting point of dynamic panel GMM. Sequential exogeneity earns its keep here: precisely because $y_{i,t-1}$ is only predetermined rather than strictly exogenous, its earlier lags $y_{i,t-2}, y_{i,t-3}, \ldots$ are clean enough to serve as instruments. Section 4 will see Arellano-Bond push this idea to its limit, using every valid lagged level.

The main points of this section can be summed up as follows: the consistency of fixed effects depends on strict exogeneity, that the error be unrelated to covariates in every period; the random effects assumption is stronger (the individual effect uncorrelated with covariates) but more efficient, and the Hausman test chooses between the two; a lagged dependent variable by nature satisfies only the weaker sequential exogeneity and violates strict exogeneity, causing fixed effects, when demeaning, to mix the current error into the lag and produce Nickell bias, which is downward in direction, $O(1/T)$ in magnitude, and does not vanish with $N$; the remedy is to switch to differencing to remove $\alpha_i$, then use a sufficiently distant lagged level $y_{i,t-2}$ as an instrument, whose validity is guaranteed by sequential exogeneity.

## 4 Estimation: Dynamic Panel GMM

Identification has laid out the road of differencing plus lagged levels; this section explains how to make it into an efficient estimator. The main line is how Arellano-Bond wring out every valid instrument (4.1), how to test whether these instruments hold up (4.2), how Blundell-Bond come to the rescue when persistence approaches 1 (4.3), and the practical matters of too many instruments and clustering (4.4).

### 4.1 Arellano-Bond: Building Instruments From Inside the Model

Anderson-Hsiao use only a single lag $y_{i,t-2}$, wasting the earlier lags. The insight of Arellano and Bond (1991) is that as time advances, more and more valid instruments become available, and all of them should be used. Look at the instruments for the differenced equation in each period: at $t = 3$ (the first period that can be differenced and instrumented), the valid instrument for $\Delta y_{i2}$ is $y_{i1}$; at $t = 4$, the valid instruments for $\Delta y_{i3}$ are $y_{i1}, y_{i2}$; and on to $t = T$, where the instruments are $y_{i1}, \ldots, y_{i,T-2}$. The valid lagged levels grow linearly with the period, and all of them are stacked into one system of moment conditions:

$$\mathbb{E}[y_{is}\, \Delta\varepsilon_{it}] = 0, \quad s \leq t - 2.$$

That the instruments are "built from inside the model" is the most beautiful part of this method: no external instrument is needed, the model's own historical values are the instruments. Solve these moment conditions with the GMM machinery of Chapter 5, weighting by the efficient weight matrix in two-step feasible GMM, and you get the Arellano-Bond estimator, usually called difference GMM. In Nimbus it gives $\hat\rho = 0.705$, almost dead-on the true value of 0.70, cleanly pulling out the truth between pooled OLS's 0.946 and fixed effects' 0.425.

### 4.2 Do the Instruments Hold Up: Overidentification and Serial Correlation Tests

Difference GMM uses far more instruments than parameters and is heavily overidentified, which is both the source of its precision and the reason for two tests that must be run.

The first is the overidentification test, the Hansen $J$ of Chapter 5 (the older version under homoskedasticity is called Sargan). It tests whether all the lagged instruments are jointly valid: if they are all clean, the $\rho$ estimated from different subsets of instruments should agree, and $J$ is small; if some instruments are in fact endogenous (say the error is serially correlated, so $y_{i,t-2}$ is no longer clean), $J$ grows and rejects. In Nimbus the Sargan $p$-value is 0.29, no rejection, and the instruments pass their check-up.

The second is the Arellano-Bond serial correlation test, which zeroes in on one fatal hazard: if the error $\varepsilon_{it}$ is serially correlated, then $y_{i,t-2}$ correlates with $\Delta\varepsilon_{it}$ and the instrument fails. The cleverness of the test is to look at the serial correlation of the differenced residual $\Delta\hat\varepsilon_{it}$. Differencing itself manufactures first-order serial correlation ($\Delta\varepsilon_{it}$ and $\Delta\varepsilon_{i,t-1}$ share $\varepsilon_{i,t-1}$), so first-order correlation AR(1) appearing is normal and expected. The real worry is second-order correlation AR(2): if $\Delta\varepsilon_{it}$ and $\Delta\varepsilon_{i,t-2}$ are correlated, it means the level error $\varepsilon_{it}$ is itself serially correlated and the instrument $y_{i,t-2}$ is dirty. So the correct diagnosis is: AR(1) significant (as it should be), AR(2) insignificant (only then are the instruments clean). In Nimbus the AR(1) $p$-value is 0.000 (significant as expected) and the AR(2) $p$-value is 0.33 (insignificant), and the validity of the instruments is supported.

### 4.3 Blundell-Bond: The Remedy When Persistence Approaches 1

Difference GMM has a fatal weak spot, exposed when persistence $\rho$ approaches 1: the lagged levels become weak instruments. The reasoning is this: when $\rho$ approaches 1, the level series of $y$ approaches a random walk, and the past level of a random walk $y_{i,t-2}$ has almost no predictive power for its difference $\Delta y_{i,t-1}$ (the increments of a random walk are unpredictable). Once the instruments are weak, difference GMM repeats the weak-identification story of Chapter 5: biased estimates and exploding standard errors.

The remedy of Blundell and Bond (1998) is to add a new set of moment conditions and enlarge the system. On top of "using lagged levels as instruments for the differenced equation," they add "using lagged differences as instruments for the levels equation":

$$\mathbb{E}[\Delta y_{i,t-1}\,(\alpha_i + \varepsilon_{it})] = 0.$$

That is, let $\Delta y_{i,t-1}$ serve as an instrument for the undifferenced levels equation. This new set of moment conditions remains valid near a random walk (differences still predict levels), rescuing identification. The cost is that it requires an additional assumption: the initial condition satisfies mean stationarity, roughly that the extent to which each unit's starting point deviates from its long-run mean is unrelated to $\alpha_i$. The estimator that uses both sets of moment conditions together is called system GMM. In Nimbus we build a near-random-walk variant with $\rho = 0.90$: difference GMM gives $\hat\rho = 0.809$ with a standard error as high as 0.119 (biased and imprecise), while system GMM gives 0.911 with a standard error of only 0.017, greatly improving both precision and accuracy, which is exactly where system GMM earns its place over difference GMM.

### 4.4 Too Many Instruments, and Clustering

Dynamic panel GMM has a practical trap that recurs: too many instruments. Following the ladder of Section 4.1, the differenced equation runs from $t = 3$ to $t = T$, with $t - 2$ lagged instruments in period $t$, for a total of $(T-1)(T-2)/2$, growing quadratically with the period, into the hundreds or thousands once $T$ is even moderately large. Too many instruments do two kinds of harm. First, the many-instruments bias mentioned in Chapter 5: when the number of instruments is large relative to the sample size, the finite-sample bias of GMM springs back toward the fixed effects direction, partly cancelling the gains of GMM. Second, too many instruments make the overidentification test lose power; a model overfit by thousands of instruments will almost surely fail to be rejected by the Sargan test, and that beautiful $p$-value is a danger signal rather than reassurance. The remedy is to limit the number of instruments: use only a few finite lags (say $y_{i,t-2}$ through $y_{i,t-4}$), or "collapse" the instruments from different lags of the same variable into fewer columns. Reporting the instrument count and running robustness checks on it is basic hygiene for dynamic panel empirical work.

On clustering, panels have within-group correlation by nature, the errors of the same merchant across months are usually positively correlated, and standard errors must be clustered at the individual level, a direct application of Chapter 6. Default non-clustered standard errors would, as in Chapter 6, systematically understate uncertainty. In the static Nimbus example the default fixed effects standard error is 0.0175 and the merchant-clustered one is 0.0179, not much different (because the errors in this setup are nearly independent), but in real panels serial correlation is often sizeable, and clustering is standard rather than optional.

The main points of this section can be summed up as follows: Arellano-Bond difference GMM uses every valid lagged level from the model itself as an instrument, pulling persistence out of the gap between pooled OLS and fixed effects; the validity of its instruments is guarded by the Sargan overidentification test and the AR(2) serial correlation test; when persistence approaches 1 the lagged levels become weak instruments, and Blundell-Bond system GMM adds the "lagged difference as instrument for the levels equation" moment condition to the rescue, at the cost of an initial-condition stationarity assumption; in practice one must beware the bias springback and test failure caused by too many instruments, and cluster the standard errors at the individual level.

## 5 Anchoring Papers

A method only stands once it lands in real research. Three anchoring papers: one founds difference GMM and its accompanying tests, one uses system GMM to patch the shortfall near a random walk, and one applies dynamic panel GMM to the empirics of information technology and firm performance. Each is laid out along five elements: paper, method, data, results, limitations. Nickell (1981), the source of the bias derivation, was covered in detail in Section 3 and is not listed separately here.

### 5.1 Arellano and Bond (1991)

::: {.case}
Paper: "Some Tests of Specification for Panel Data: Monte Carlo Evidence and an Application to Employment Equations," The Review of Economic Studies 58(2):277-297. It founds difference GMM, the workhorse method for dynamic panels, and everything in Sections 4.1 through 4.2 of this chapter comes from it.

Method: Arellano and Bond first-difference the individual effect out of the dynamic panel equation, then use every sufficiently distant lagged level as an instrument for the differenced equation, with the instrument set growing linearly as time advances, and combine these moment conditions efficiently by two-step GMM. They also propose the accompanying diagnostics: the $m_1$ and $m_2$ tests for first- and second-order serial correlation of the differenced residuals (expecting AR(1) significant, AR(2) insignificant), and the Sargan overidentification test for instrument validity. This combination of "building instruments from inside the model plus built-in diagnostics" became the standard workflow of dynamic panel empirical work.

Data: an employment panel of UK firms, about 140 firms spanning 1976 to 1984, used to estimate a dynamic employment equation, that is, how a firm's current employment depends on last period's employment, wages, and capital.

Results: the paper uses Monte Carlo to show the bias improvement of difference GMM over fixed effects in short panels, and also gives concrete dynamic estimates on the employment equation, demonstrating how the $m_2$ and Sargan tests judge whether a specification is credible in practice. It has since become one of the most heavily cited methodological papers in economics and adjacent fields.

Limitations: the paper also plants the seed of the trouble later remedied by Blundell-Bond: when persistence approaches 1, the lagged level instruments become weak and difference GMM is biased and imprecise. The many-instruments bias and test failure brought on by the rapid proliferation of instruments with $T$ are also the aspects of this method most in need of caution in practice.
:::

The exemplary significance of this paper is that it gives not just an estimator but a whole suite of diagnostics that make the estimator testable and falsifiable. The AR(2) and Sargan tests of Section 4.2 are its legacy, and those two $p$-values in Nimbus (0.33 and 0.29) are produced by following its very workflow.

### 5.2 Blundell and Bond (1998)

::: {.case}
Paper: "Initial Conditions and Moment Restrictions in Dynamic Panel Data Models," Journal of Econometrics 87(1):115-143. It patches the fatal shortfall of difference GMM when persistence approaches 1, and the system GMM of Section 4.3 of this chapter comes from it.

Method: Blundell and Bond point out that when the autoregressive parameter approaches 1, or when the variance of the individual effect is large relative to the variance of the idiosyncratic error, the lagged level instruments of difference GMM become weak and estimation deteriorates badly. Beyond the moment conditions of the differenced equation, they append a set of moment conditions for the levels equation, using lagged differences as instruments for the levels equation, forming system GMM. This new set of moment conditions relies on an initial-condition assumption (mean stationarity), whose meaning and reasonableness they argue carefully.

Data: a Monte Carlo system demonstrates the enormous improvement of system GMM over difference GMM near a random walk, accompanied by an empirical application (such as a firm-level production function or an employment panel) illustrating its use.

Results: system GMM greatly reduces bias and improves precision in the high-persistence case, becoming the standard upgrade to difference GMM. The great volume of empirical work using dynamic panels today reports both difference GMM and system GMM by default, a norm established by this paper.

Limitations: the extra efficiency of system GMM comes from that initial-condition stationarity assumption, and once it fails, the newly added level moment conditions are themselves invalid and instead introduce bias. It is also not a plain linear IV; it needs the full GMM machinery, and it too suffers from too many instruments. Whether or not to use system GMM is at bottom whether one is willing to shoulder that stationarity assumption for the sake of efficiency.
:::

Blundell-Bond and Arellano-Bond make a pair: the former more robust with fewer assumptions, the latter more efficient with more assumptions, and in practice the two are usually reported together and cross-checked against each other. The $\rho = 0.90$ variant in Nimbus (difference GMM standard error 0.119 versus system GMM's 0.017) stages the division of labor between this pair very plainly.

### 5.3 Tambe and Hitt (2012)

::: {.case}
Paper: "The Productivity of Information Technology Investments: New Evidence from IT Labor Data," Information Systems Research 23(3, Part 1 of 2):599-617. It applies dynamic panel GMM to the classic IS question of the return on IT investment, and shares a methodological core with Nimbus: using the time dimension of a panel and constructing instruments from lagged values to deal with the endogeneity of the regressors.

Method: Tambe and Hitt estimate a firm-level production function, regressing output on IT inputs, non-IT capital, and labor, with the core difficulty being the endogeneity of IT investment (high-productivity firms both invest more in IT and produce more output, and both reverse causality and unobservable productivity shocks would contaminate OLS). They line up OLS, fixed effects, and several panel-based endogeneity-robust estimators for comparison, including Arellano-Bond-type dynamic panel GMM (using lagged levels as instruments to handle the endogenous regressor), placing the estimate of the IT return under different identification assumptions to cross-check one another.

Data: a newly constructed firm-level panel of IT labor inputs, matching IT personnel expenditure data with financial data, covering roughly 1987 to 2006, and, unlike earlier studies focused only on Fortune 1000 large firms, including many medium-sized and service-sector firms.

Results: the IT output elasticity measured after correcting for endogeneity falls only about ten percent (the IT labor elasticity in the dynamic panel GMM column is about 0.041), implying that the endogeneity bias in past studies of IT value was in fact not large; the IT return is significantly lower for medium-sized firms than for large firms and is realized more slowly in large firms; and the marginal output of IT spending from 2000 to 2006 exceeds that of any prior period, hinting that IT-driven process innovation continued after the internet boom.

Limitations: the main thread of the paper is the measurement of IT productivity, and dynamic panel GMM is one of several robustness tools it uses against endogeneity rather than the sole workhorse; production function estimation also relies on the exogeneity of inputs and on functional-form specification, and dynamic panel GMM too suffers from the too-many-instruments and weak-instruments problems this chapter keeps stressing.
:::

The significance of this paper for the chapter is that it returns dynamic panel GMM to a real IS research question: whether IT investment is worth it, the difficulty being that IT investment is not sprinkled down at random. It demonstrates that when the covariate is endogenous and there is panel structure, building instruments from the model's own lagged values is a ready-made path in IS empirical work, and it also reminds us of the prudent practice of using it alongside other identification strategies for mutual corroboration.

Taken together, the three make the meaning of anchoring clear: Arellano-Bond found difference GMM and its accompanying diagnostics, Blundell-Bond use system GMM to patch the shortfall near a random walk, and Tambe and Hitt apply the toolkit to measuring the return on IT investment. They share the chapter's core, that the time dimension of a panel can both difference out unobservable heterogeneity and build instruments from history, handling dynamics and endogeneity together.

## 6 A Full Walkthrough on the Nimbus Data

Now we run the earlier tools on Nimbus from start to finish. The code below uses R 4.5.3 with the random seed fixed at 737, and every number cited in the text comes from the actual output of these runs. Dynamic panel GMM uses `plm::pgmm`.

### 6.1 The Data Generating Process

We design two blocks of data and one seed. The static block installs fixed effects: $y_{it} = \beta x_{it} + \alpha_i + \varepsilon_{it}$ with $\beta = 0.5$, and merchant talent $\alpha_i$ positively correlated with promotion $x_{it}$ (high-talent merchants promote more), which biases pooled OLS and random effects and leaves fixed effects unbiased. The dynamic block stages Nickell and GMM: $y_{it} = \rho\, y_{i,t-1} + \alpha_i + \varepsilon_{it}$ with $\rho = 0.70$, $N = 500$ merchants and $T = 8$ months, plus a near-random-walk variant with $\rho = 0.90$.

```r
set.seed(737)
## Static block: alpha positively correlated with x
alpha <- rnorm(600); a <- alpha[rep(1:600, each = 6)]
x <- 0.7 * a + rnorm(3600)
y <- 0.5 * x + a + rnorm(3600)                    # beta = 0.5
## Dynamic block: AR(1) plus merchant effect, stationary start
rho <- 0.70; al <- rnorm(500)
Y <- matrix(0, 500, 8 + 50)                       # 50-period burn-in
Y[,1] <- al/(1-rho) + rnorm(500, 0, 1/sqrt(1-rho^2))
for (s in 2:58) Y[,s] <- rho * Y[,s-1] + al + rnorm(500)
```

The dynamic block uses a 50-period burn-in to wash out the transient of the initial condition and ensure entry into stationarity. The true values $\beta = 0.5$, $\rho = 0.70$, and the $\rho = 0.90$ of the near-random-walk variant are the targets for all the estimates that follow.

### 6.2 Fixed Effects as Hero: Static

On the static block we compare three estimators.

```r
library(plm)
d <- pdata.frame(static, index = c("id", "t"))
pool <- plm(y ~ x, d, model = "pooling")
re   <- plm(y ~ x, d, model = "random")
fe   <- plm(y ~ x, d, model = "within")
phtest(fe, re)                                     # Hausman
```

Pooled OLS gives $\hat\beta = 0.951$, badly inflated by the positive correlation between talent and promotion. Random effects partially demean and give $\hat\beta = 0.739$, shaving off some contamination but not cleanly (because $\alpha_i$ really is correlated with $x$). Fixed effects within-demean and give $\hat\beta = 0.506$, hugging the true value. The Hausman statistic reaches 3362 with a $p$-value of zero, random effects are decisively rejected, and fixed effects are confirmed as the choice. Fixed effects are the undisputed hero of the static panel.

### 6.3 The Hero Turns Traitor: Dynamics and Nickell

On the dynamic block, regress this month's GMV on last month's.

```r
d  <- pdata.frame(dyn_main, index = c("id", "t"))
pool <- plm(y ~ lag(y, 1), d, model = "pooling")
fe   <- plm(y ~ lag(y, 1), d, model = "within")
```

Pooled OLS gives $\hat\rho = 0.946$ (inflated by merchant talent) and fixed effects give $\hat\rho = 0.425$ (pressed down by Nickell), with the truth of 0.70 caught between them. The two wrong numbers of Section 1 are confirmed here. The same fixed effects that were a hero in Section 6.2 become a villain here, and the gap is entirely due to the lagged dependent variable violating strict exogeneity.

The $1/T$ character of Nickell bias can be demonstrated directly. Fix $\rho = 0.70$, let the panel length $T$ grow from 4 to 32, and watch the fixed effects estimate of $\rho$:

| $T$ | 4 | 8 | 16 | 32 |
|---|---|---|---|---|
| Fixed effects $\hat\rho$ | 0.082 | 0.417 | 0.570 | 0.640 |
| Arellano-Bond $\hat\rho$ | 0.705 | 0.686 | 0.700 | 0.701 |
| Nickell leading term | 0.133 | 0.457 | 0.587 | 0.645 |

The fixed effects estimate climbs steadily toward the true value of 0.70 as $T$ grows, as low as 0.082 at $T = 4$ and back to 0.640 at $T = 32$, and the leading-term prediction of $-(1+\rho)/(T-1)$ (third row) traces this climbing path quite accurately (the 0.417 in the $T = 8$ cell of the table comes from this independent resampled experiment, and together with the 0.425 from the main sample they are two draws from the same DGP, so a slight discrepancy is normal). Arellano-Bond, meanwhile, holds steady near 0.70 at every $T$, with no $1/T$ drift, and this is the intuitive picture of "adding periods can rescue fixed effects, but the truly consistent one is GMM." Note that fixed effects run a bit below even the leading term at short $T$, because the higher-order terms of the full Nickell formula are also pulling down in short panels.

### 6.4 Arellano-Bond and Blundell-Bond

Using lagged levels as instruments, run difference GMM and system GMM.

```r
ab <- pgmm(y ~ lag(y,1) | lag(y, 2:99), d, effect = "individual",
           model = "twosteps", transformation = "d")   # difference GMM
bb <- pgmm(y ~ lag(y,1) | lag(y, 2:99), d, effect = "individual",
           model = "twosteps", transformation = "ld")  # system GMM
```

Difference GMM gives $\hat\rho = 0.705$ and system GMM gives 0.701, both landing precisely on the true value of 0.70, cleanly pulling out the truth from the gap between pooled OLS and fixed effects. The validity of the instruments survives testing: Sargan overidentification $p = 0.29$ (no rejection), AR(1) $p = 0.000$ (significant as expected), AR(2) $p = 0.33$ (insignificant, instruments clean). Only a result with this full complement of diagnostics is a persistence estimate fit to write into a report.

The value of system GMM over difference GMM only shows up when persistence approaches 1. Use the near-random-walk variant with $\rho = 0.90$:

```r
ab9 <- pgmm(y ~ lag(y,1) | lag(y,2:99), dp, model="twosteps", transformation="d")
bb9 <- pgmm(y ~ lag(y,1) | lag(y,2:99), dp, model="twosteps", transformation="ld")
```

Difference GMM gives $\hat\rho = 0.809$ with standard error 0.119, both biased (far from 0.90) and imprecise (large standard error), exactly the symptom of lagged levels degenerating into weak instruments. System GMM gives 0.911 with standard error 0.017, greatly improving both accuracy and precision. Adding the "lagged difference as instrument for the levels equation" moment condition rescues identification near a random walk.

![Left: three estimates of the promotion effect in the static panel (truth 0.50). Pooled OLS (0.95) and random effects (0.74) are inflated by merchant talent, and only fixed effects (0.51) recover the truth. Right: four estimates of persistence in the dynamic panel (truth 0.70). Pooled OLS (0.95) runs high, fixed effects (0.42) are pressed down by Nickell bias, the two straddle the truth, and both Arellano-Bond and Blundell-Bond GMM recover 0.70.](assets/fig/fig_07_estimators.svg)

![Left: Nickell bias fading with panel length $T$. Fixed effects (blue) badly understate persistence at short $T$ and climb toward the truth of 0.70 as $T$ grows, with the $-(1+\rho)/(T-1)$ leading term (grey dashed) tracing the path; Arellano-Bond (navy) holds steady at 0.70 across all $T$ with no $1/T$ drift. Right: when persistence approaches 1 (truth 0.90), difference GMM (gold) is biased with an extremely wide confidence interval because the lagged levels degenerate into weak instruments, while system GMM (navy) rescues precision through the extra moment conditions of the levels equation.](assets/fig/fig_07_nickell.svg)

### 6.5 The Full Reconciliation

Putting the key numbers side by side:

| Quantity | Estimator | Result | Target |
|---|---|---|---|
| Promotion effect $\beta$ | pooled OLS / RE / FE | 0.951 / 0.739 / 0.506 | truth 0.50, only FE correct |
| Hausman | FE vs RE | $\chi^2 = 3362$, $p \approx 0$ | reject RE, use FE |
| Persistence $\rho$ | pooled OLS / FE | 0.946 / 0.425 | truth 0.70 straddled |
| Persistence $\rho$ | difference GMM / system GMM | 0.705 / 0.701 | truth 0.70, both recovered |
| AB diagnostics | Sargan / AR(1) / AR(2) | 0.29 / 0.000 / 0.33 | instruments valid |
| $\rho \to 1$ | difference / system GMM | 0.809 (SE 0.119) / 0.911 (SE 0.017) | truth 0.90 |

The reconciliation of this section can be summed up as follows: in the static panel, fixed effects are the only estimator to recover the true promotion effect, and Hausman decisively rejects random effects; the same fixed effects, in the dynamic panel, are pressed down by Nickell bias and, together with the inflated pooled OLS, straddle the true persistence; Nickell bias fades with $T$ at rate $1/T$ and the leading term predicts it well, while Arellano-Bond is consistent at every $T$; both difference and system GMM recover persistence to 0.70 and pass the Sargan and AR(2) diagnostics; when persistence approaches 1 the lagged level instruments of difference GMM become weak, and system GMM rescues precision through the level moment conditions.

## 7 Failure Modes and Robustness

The setups in the simulation are manufactured, and in real panels they hide deeper. This section lays out the most common failure modes and the responses to them.

Fixed or random effects is the first choice to answer head-on. Random effects are more efficient and can also estimate time-invariant covariates, but they wager on the strong assumption that "the individual effect is uncorrelated with the covariates," which in observational data usually does not hold (units that select into a state often differ in unobservable dimensions to begin with). The Hausman test gives a data-driven criterion, but it is not a cure-all: a failure to reject by Hausman does not mean random effects are correct, only perhaps that the test lacks power; and the difference of variances that Hausman uses can be non-positive-definite in finite samples, making the statistic negative and uninterpretable. The safe default is that, unless there is strong reason to believe the individual effect is exogenous and time-invariant covariates truly need estimating, use fixed effects and put robustness ahead of efficiency.

Strict exogeneity is the most easily overlooked weak spot of fixed effects. It is violated not only by a lagged dependent variable but by any covariate that "responds to past shocks": promotion adjusted to last period's performance, price adjusted to last period's demand, investment adjusted to last period's profit; these feedbacks all make the covariate correlate with past errors, breaking strict exogeneity and biasing fixed effects. The key to judging is to ask: could the value of this covariate be influenced by past outcomes. If so, treat it as a predetermined variable and use an instrument under sequential exogeneity (a lagged value) rather than putting it straight into a fixed effects regression.

Nickell bias lurks in any short panel with a lagged dependent variable, and large sample size is no help. Any setup with a lagged outcome on the right-hand side and a not-large $T$ cannot trust the fixed effects estimate of persistence directly and should use dynamic panel GMM. A hidden variant is that even if what you care about is not persistence itself but the effect of some covariate, as long as a lagged dependent variable is included in the model to control for dynamics, Nickell bias will seep along the correlation into that covariate's coefficient, so the seemingly harmless move of "adding a lagged dependent variable to control for dynamics a bit" biases the entire fixed effects regression.

Dynamic panel GMM has a few things of its own to watch. Too many instruments is the most common pit: the instrument count grows quadratically with $T$, and too many make GMM lean toward fixed effects and rob the Sargan test of power; a model with thousands of instruments will almost surely fail to be rejected by Sargan, and that is not reassurance but an alarm of overfitting. The response is to limit the lag depth or collapse the instruments, and report the instrument count and run robustness checks on it. When persistence approaches 1 difference GMM has weak instruments and must rely on system GMM, but the extra moment conditions of system GMM depend on initial-condition stationarity, and when that assumption fails system GMM is instead worse, so reporting both GMMs and seeing whether they agree is a basic robustness exercise. Serial correlation in the error makes the lagged instruments fail, and the AR(2) test is the gatekeeper; once it rejects, one must use more distant lags as instruments or reexamine the specification.

Stringing these failure modes together, the robustness of panel methods rests at bottom on two judgments: first, whether your unobservable heterogeneity really is time-invariant and your covariates really are strictly exogenous, which decides whether fixed effects can be used; second, whether your model is dynamic and whether $T$ is long enough, which decides whether Nickell bias matters and whether dynamic panel GMM should be brought in. Once these two judgments are thought through clearly, the panel's gift of differencing out heterogeneity is truly in hand; if they are not, demeaning is a double-edged sword that may scrub out contamination for you or quietly carve a $1/T$ scar into your dynamic parameter.

## 8 Further Reading

::: {.readings}
Required reading, in suggested order:

- Wooldridge (2010, *Econometric Analysis of Cross Section and Panel Data*, 2nd ed.), Chapters 10 to 11. The clearest textbook treatment of fixed effects, random effects, and strict versus sequential exogeneity, the source of the framework for the static panel part of this chapter.
- Nickell (1981, *Econometrica*). The original derivation of the dynamic fixed effects bias, the source of that $-(1+\rho)/(T-1)$ in Section 3.4 of this chapter, short and essential.
- Arellano and Bond (1991, *Review of Economic Studies*). The founding work of difference GMM and the AR(2) and Sargan diagnostics, the starting point of dynamic panel empirical work.
- Blundell and Bond (1998, *Journal of Econometrics*). System GMM, to understand why nothing else will do when persistence approaches 1, read together with Arellano-Bond.

Further reading:

- Anderson and Hsiao (1981, *JASA*). The original idea of using lagged levels as instruments, the direct forerunner of Arellano-Bond, read it to see where the instruments come from.
- Hausman (1978, *Econometrica*). The original source of specification testing, the origin of this chapter's Hausman test, of a piece with the QMLE tests of Chapter 5.
- Mundlak (1978, *Econometrica*). Writing the individual effect as a function of the within means of the covariates, building a bridge between fixed and random effects (correlated random effects), another angle for understanding what the Hausman test tests.
- Arellano and Bover (1995, *Journal of Econometrics*). The other half of the theoretical foundation for the level moment conditions of system GMM, read together with Blundell-Bond.
- Roodman (2009, *Stata Journal*). The most practical operating guide to dynamic panel GMM, giving the clearest account of too many instruments, collapsing instruments, and reporting both GMMs, an essential checklist before getting your hands dirty.
- Baltagi (*Econometric Analysis of Panel Data*). The authoritative monograph on panel econometrics, covering static to dynamic and linear to nonlinear in full, a desk reference.
:::

::: {.apa-refs}
- Anderson, T. W., & Hsiao, C. (1981). Estimation of dynamic models with error components. *Journal of the American Statistical Association, 76*(375), 598-606. https://doi.org/10.1080/01621459.1981.10477691
- Arellano, M., & Bond, S. (1991). Some tests of specification for panel data: Monte Carlo evidence and an application to employment equations. *The Review of Economic Studies, 58*(2), 277-297. https://doi.org/10.2307/2297968
- Arellano, M., & Bover, O. (1995). Another look at the instrumental variable estimation of error-components models. *Journal of Econometrics, 68*(1), 29-51. https://doi.org/10.1016/0304-4076(94)01642-D
- Baltagi, B. H. (2021). *Econometric analysis of panel data* (6th ed.). Springer. https://doi.org/10.1007/978-3-030-53953-5
- Blundell, R., & Bond, S. (1998). Initial conditions and moment restrictions in dynamic panel data models. *Journal of Econometrics, 87*(1), 115-143. https://doi.org/10.1016/S0304-4076(98)00009-8
- Hausman, J. A. (1978). Specification tests in econometrics. *Econometrica, 46*(6), 1251-1271. https://doi.org/10.2307/1913827
- Mundlak, Y. (1978). On the pooling of time series and cross section data. *Econometrica, 46*(1), 69-85. https://doi.org/10.2307/1913646
- Nickell, S. (1981). Biases in dynamic models with fixed effects. *Econometrica, 49*(6), 1417-1426. https://doi.org/10.2307/1911408
- Roodman, D. (2009). How to do xtabond2: An introduction to difference and system GMM in Stata. *The Stata Journal, 9*(1), 86-136. https://doi.org/10.1177/1536867X0900900106
- Tambe, P., & Hitt, L. M. (2012). The productivity of information technology investments: New evidence from IT labor data. *Information Systems Research, 23*(3, Part 1 of 2), 599-617. https://doi.org/10.1287/isre.1110.0398
- Wooldridge, J. M. (2010). *Econometric analysis of cross section and panel data* (2nd ed.). MIT Press.
:::
