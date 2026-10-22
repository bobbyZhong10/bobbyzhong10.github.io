---
title: "Linear Panels: FE, RE, FD & CRE"
subtitle: "Unobserved Effects and the Within Transformation"
seriesline: "Foundations of Information Systems Economics · Chapter 10"
series: "Foundations of Information Systems Economics: Empirical Methods for Information Systems Economics"
sidebartitle: "Foundations of Information Systems Economics"
sidebarsub: "Theory · Methods · Applications"
chapterlabel: "Chapter 10 · Linear Panels: FE, RE, FD & CRE"
---

## Introduction

Nimbus is a fictional online marketplace whose merchants decide each month how much to spend on promotions inside the marketplace. The merchant growth team wants one plain number: when the same merchant spends USD 100 more on promotion, how much higher, on average, is its GMV (gross merchandise value) that month? The data at hand are monthly records on the same set of merchants. The appeal of panel data is that they observe the same unit repeatedly, which lets us sweep out unobservable differences that do not vary over time. For static models, fixed effects (FE) are often a reliable workhorse.

This chapter covers the static side of linear panels: when the individual effect is correlated with the covariates, why pooled OLS and random effects (RE) are contaminated, how FE, first differences (FD), and correlated random effects (CRE) deal with that correlation, and under what conditions even FE goes wrong. Section 1 starts with a wrong number that is almost double the truth, and Section 2 writes down the model and the estimand of each data block. Section 3, the most detailed in the chapter, gives in turn the assumptions and theorems for pooled OLS, RE, FE, and FD and the tests of strict exogeneity, and then covers how RE relates to FE, the Mundlak regression and Chamberlain's approach, the classical and robust Hausman tests, the Hausman-Taylor estimator, random slopes, and unbalanced panels. Section 4 covers estimation and software, Section 5 reads three anchoring papers, Section 6 runs the whole toolkit on Nimbus, and Section 7 discusses failure modes. The chapter uses the sandwich variance of Section 3.4 of Chapter 5 and the cluster-robust variance of Section 3.4 of Chapter 8; RE, as feasible GLS, uses the GLS of Chapter 9; Chamberlain's approach uses the minimum distance (MD) of Chapter 7; and Hausman-Taylor uses the IV of Chapter 9. Dynamic panels, with a lagged dependent variable on the right-hand side, are left to Chapter 11; the relation between two-way FE (TWFE) and DiD is taken up again in Chapter 19.

## 1 A Promotion Effect Nearly Twice the Truth

::: {.case}
Block P1 is Nimbus's inherited static panel: 600 merchants, 6 months each, with the variables kept in standardized units. In the simulation the true effect of promotion on GMV is $\beta = 0.5$; merchants with more operating talent also promote more, and the regression coefficient of promotion on talent is 0.7. Talent is unobserved, and it raises both GMV and promotion.
:::

The analyst pools the 3,600 merchant-months, runs OLS, gets 0.951 with a merchant-clustered 95% interval of [0.908, 0.994], and reads it as the effect of promotion; yet the true value in the simulation is only 0.5, so the estimate is almost double. Switch to RE and the answer falls to 0.739; switch to FE and the answer is 0.506, with an interval [0.471, 0.541] that covers the truth. The Hausman statistic comparing FE and RE is 3362.02. Same data, three estimators, three answers. The differences come not from sampling error but from the different assumptions each estimator makes about the correlation between talent and promotion. Pooled OLS errs not in its arithmetic but in its comparison: it compares merchants who promote a lot with merchants who promote little, and merchants who promote a lot are better at business to begin with, so nearly half of the GMV gap is the work of talent. The value of a panel is that it lets us make a different comparison.

This chapter answers three questions: why pooled OLS is biased, and by how much (Section 3.1); why RE corrects only part of the bias (Section 3.2); and what exactly FE estimates, and when it too goes wrong. The answer to the third question is spread across Nimbus's other data blocks: when promotion responds to last month's shocks (feedback), FE gives 0.574 rather than 0.5; when slopes vary across merchants and the merchants whose promotion fluctuates more also have steeper slopes, it gives 0.683, which targets a weighted average; and when merchants leave the platform because of poor performance, it gives 0.444.

## 2 The Economic Model and the Estimand

### 2.1 The Panel Model and Individual Heterogeneity

The basic panel model writes the outcome as the sum of observable covariates, an individual effect, and an idiosyncratic error:

$$Y_{it} = X_{it}'\beta + \alpha_i + \varepsilon_{it}.$$

The subscript $i = 1, \ldots, N$ is the unit (a merchant in Nimbus) and $t = 1, \ldots, T$ is the time point (a month). $X_{it}$ is a $k \times 1$ vector of time-varying covariates and $\beta$ is the coefficient we want. The key new character is $\alpha_i$, the unobserved effect, which captures every individual characteristic that is time-invariant and not in $X$: a merchant's operating talent, the natural scale of its category, the capability of its founding team. $\varepsilon_{it}$ is the idiosyncratic error, which varies across merchants and months. $\alpha_i$ is unobservable, and it is almost surely correlated with $X_{it}$, which is the core difficulty of identification: a high-talent merchant has high GMV and is also more willing to spend on promotions. If $\alpha_i$ is left in the error term and OLS is run directly, the estimate of $\beta$ is contaminated by this correlation. Folding it into the error term gives $\nu_{it} = \alpha_i + \varepsilon_{it}$, the composite error. In a balanced panel every merchant has $T$ periods; when merchants leave partway through or join late, the panel is unbalanced; the number of periods in which unit $i$ is observed is written $T_i$, and $T_i = T$ in a balanced panel. All the large-sample theory in this chapter lets the number of merchants $N$ go to infinity with the number of months $T$ fixed, which matches the shape of data on a few thousand merchants observed for anywhere from a few months to a few years.

The names RE and FE originally came from whether $\alpha_i$ is treated as a random variable or as a parameter to be estimated; Wooldridge (2010) §10.2.1 argues that when many units are drawn at random from a population, $\alpha_i$ should almost always be treated as a random variable drawn from the population together with $Y_{it}$ and $X_{it}$, and that the only thing that matters is whether it is correlated with $X_{it}$. This chapter adopts that view: FE does not mean that $\alpha_i$ is nonrandom, but that it is allowed to be arbitrarily correlated with the regressors; RE assumes the two are unrelated; CRE sits between them and writes out the form of the correlation.

Panels offer a way out: because the same merchant is observed repeatedly, $\alpha_i$ can be removed. The most common device is the within transformation, also called demeaning: for each unit, subtract its own average over time, $\bar Y_i = T_i^{-1}\sum_t Y_{it}$ and $\bar X_i = T_i^{-1}\sum_t X_{it}$, where the sums run only over unit $i$'s $T_i$ observed periods. Because $\alpha_i$ does not vary over time, it equals its own time mean, and the subtraction removes it:

$$Y_{it} - \bar Y_i = (X_{it} - \bar X_i)'\beta + (\varepsilon_{it} - \bar\varepsilon_i).$$

Write $\ddot Y_{it} = Y_{it} - \bar Y_i$ and $\ddot X_{it} = X_{it} - \bar X_i$. OLS on the demeaned data is the FE estimator, also called the within estimator. $\alpha_i$ is wiped out completely, and however it correlates with $X$ it no longer causes any bias.

::: {.intuition}
Why does FE win? In a phrase: it compares each merchant only with itself. It does not ask "Do merchants that promote more have higher GMV?" (a comparison contaminated by talent, because merchants that promote more tend to be more talented) but rather "In the months when a given merchant promotes more, is its GMV higher than in its own months of less promotion?" A single merchant's talent is the same constant in different months, and a within comparison cancels it automatically. This idea is the whole source of FE's power, and also its whole cost: it uses only within (over-time) variation and throws away all between (across-unit) variation, so it cannot estimate any time-invariant covariate and is imprecise when within-unit variation is small.
:::

### 2.2 The Estimand and Nimbus's Seven Data Blocks

The estimand of this chapter is "the effect of promotion once merchant talent is held fixed":

$$\mathbb{E}[Y_{it} \mid \mathbf{X}_i, \alpha_i] = X_{it}'\beta + \alpha_i,$$

where $\mathbf{X}_i$ stacks unit $i$'s regressors from all periods. $\beta$ answers this question: when the same merchant spends one more unit on promotion, how much higher is its GMV that month, on average? The conditioning set contains $\mathbf{X}_i$ from all periods, which means that promotion in other months has no further effect on this month's GMV; Section 3.1 shows how strong this condition is. Reading $\beta$ as a causal effect rests on the condition that, given talent, the month-to-month swings in promotion are unrelated to the GMV shocks of every month; it does not require promotion to be randomly assigned across merchants, and that is the advantage of a panel over a cross-section. In the new data blocks $Y_{it}$ is monthly GMV in thousands of USD and $X_{it}$ is monthly promotion spending in hundreds of USD, so $\beta = 0.5$ means that each USD 100 of promotion brings USD 500 of GMV.

Slopes need not be the same for every merchant. When merchant $i$ has its own slope $b_i$, management usually wants the average slope $\beta = \mathbb{E}[b_i]$: if every merchant spent USD 100 more, how much more GMV would there be on average? Section 3.8 shows that what FE converges to is not necessarily this average but a different weighted average, $\beta_*$. Time-invariant variables have their own estimand too: in block P4, "Nimbus Certified merchant" status $W_i$ raises monthly GMV by $\gamma = 1.5$, that is, USD 1,500, and FE cannot estimate it.

| Block | Merchants and months | Design | Target and true value | Used in |
|---|---|---|---|---|
| P1 | 600, 6 | Talent positively correlated with promotion, standardized units | $\beta = 0.5$ | Sections 3.1 to 3.6 |
| P2 | 500, 6 | The same monthly shocks built into an iid error and into a random-walk error | $\beta = 0.5$ | Section 3.4 |
| P3 | 1,000, 6 | Merchants hit by a bad shock last month promote more this month | $\beta = 0.5$ | Sections 3.4 and 3.5 |
| P4 | 2,000, 6 | Certified status is time-invariant and correlated with talent | $\gamma = 1.5$ | Section 3.7 |
| P5, P5b | 800, 8 | Slopes vary with the within-unit fluctuation of promotion; in P5b only with its level | Average slope 0.5; FE probability limit in P5 0.66125 | Section 3.8 |
| P6 | 4,000, 8 | Exit depends on talent, or on the current month's shock | $\beta = 0.5$ | Sections 3.5 and 3.9 |
| P7 | 1,000, 5 | Errors heteroskedastic and serially correlated; RE.1 holds but RE.3 fails | The Hausman null is true | Section 3.6 |

Each new block changes one thing only: it makes one assumption fail and leaves the rest as they were. P2 changes the serial correlation of the error, and P3 changes how promotion is decided: merchants whose sales were poor last month promote more this month, which is what any merchant with an operations team would do. P4 adds a time-invariant regressor, P5 makes promotion work differently for different merchants, P6 lets merchants leave the platform, and P7 changes only the second moments of the error, not the first. That way, whenever an estimator goes wrong, we know which assumption is responsible.

## 3 Identification and Large-Sample Theory

This section first uses Nimbus's numbers to show what is to be proved, then states numbered assumptions and theorems, followed by proof sketches of two to six steps, each step naming the assumptions it uses. Finally, an "In words" paragraph answers four questions: under which assumptions it estimates which quantity, what it estimates when the specification is wrong, who sets the weights, and what Nimbus's numbers say. The two parts of Wooldridge (2010) on panels have numbered assumptions but no numbered theorems, so the theorem boxes cite assumption labels, sections, and equation numbers; the theorem numbering is this chapter's own. The large-sample theory is always $N \to \infty$ with $T$ fixed.

### 3.1 Strict Exogeneity and Pooled OLS

Pooled OLS on P1 gives 0.951, and its probability limit can be computed directly: talent $\alpha_i$ has variance 1, and promotion equals 0.7 times $\alpha_i$ plus independent noise with variance 1, so $\mathrm{Cov}(X_{it}, \alpha_i)/\mathrm{Var}(X_{it}) = 0.7/(0.7^2 + 1)$, and pooled OLS converges to 0.9698. The estimate differs from it only by sampling error.

::: {.assumption}
**POLS.1 through POLS.3 (Wooldridge 2010 §7.8.1, p. 192, applied to the composite error $\nu_{it}$)** (1) $\mathbb{E}[X_{it}\nu_{it}] = 0$, $t = 1, \ldots, T$; in this chapter's model it in effect requires $\mathbb{E}[X_{it}\varepsilon_{it}] = 0$ plus the $\mathbb{E}[X_{it}\alpha_i] = 0$ of equation (10.22). (2) $\sum_t \mathbb{E}[X_{it}X_{it}']$ has full rank. (3) Conditional homoskedasticity, with errors uncorrelated across periods.

**FE.1 (strict exogeneity; Wooldridge 2010 §10.5.1, p. 301, that is, equation (10.14))** $\mathbb{E}[\varepsilon_{it} \mid \mathbf{X}_i, \alpha_i] = 0$, $t = 1, \ldots, T$.
:::

POLS.1 requires only contemporaneous uncorrelatedness, at the price of having to assume that talent is unrelated to promotion, which almost never holds in observational data. FE.1 is the reverse: it allows $\alpha_i$ to be arbitrarily correlated with $X$, but the conditioning set is $\mathbf{X}_i$ from all periods, so the current error must be unrelated to the covariates of every period, past, present, and future. Why so strong? Because within demeaning mixes all periods together: $\bar X_i$ contains $X$ from every period and $\bar\varepsilon_i$ contains $\varepsilon$ from every period, so as soon as $X$ in any period correlates with $\varepsilon$ in any period, the two demeaned quantities correlate. Strict exogeneity seals off all these cross-period correlations at once. Conditioning on $\alpha_i$ matters a great deal: the version that does not condition on talent, $\mathbb{E}[Y_{it} \mid \mathbf{X}_i] = X_{it}'\beta$, almost surely fails, because talent affects the promotion choice in every month, so promotion in other months "predicts" this month's GMV through talent (Wooldridge 2010 §10.2.2 gives the same example with farm inputs).

The classic case that strict exogeneity rules out is feedback: if this month's promotion was partly raised because last month's GMV came in surprisingly low ($\varepsilon_{i,t-1} < 0$), then $X_{it}$ correlates with $\varepsilon_{i,t-1}$ and strict exogeneity no longer holds; block P3 is built this way. A lagged dependent variable $Y_{i,t-1}$ is by its very nature a covariate that "responds to past shocks," and this is the cue for Section 3.1 of Chapter 11.

::: {.theorem}
**Theorem 1 (pooled OLS; Wooldridge 2010 Theorem 7.7 and equation (10.22))** Under POLS.1 and POLS.2, pooled OLS is consistent and asymptotically normal; with POLS.3 added, the usual OLS variance applies. If POLS.1 fails because $\mathbb{E}[X_{it}\alpha_i] \neq 0$, but $\mathbb{E}[X_{it}\varepsilon_{it}] = 0$, then

$$\begin{aligned}
\mathrm{plim}\,\hat\beta_{POLS} &= \beta \\
&\quad + \Big(\sum_{t}\mathbb{E}[X_{it}X_{it}']\Big)^{-1} \\
&\quad \times \sum_{t}\mathbb{E}[X_{it}\alpha_i].
\end{aligned}$$
:::

The proof sketch has three steps. Step one (POLS.2): $\hat\beta_{POLS} - \beta$ equals the inverse of the sample second-moment matrix times the sample average of $X_{it}\nu_{it}$, and the law of large numbers makes each converge to its population moment. Step two (POLS.1): when $\mathbb{E}[X_{it}\nu_{it}] = 0$ the second factor converges to zero, and the central limit theorem, applied to sums over units, gives asymptotic normality, with a variance that is the panel version of the sandwich of Section 3.4 of Chapter 5. Step three: when POLS.1 fails, $\mathbb{E}[X_{it}\nu_{it}] = \mathbb{E}[X_{it}\alpha_i]$, and substituting into step one gives the probability limit.

In words: under POLS.1, pooled OLS estimates $\beta$; when talent and promotion are correlated, it estimates $\beta$ plus the coefficient from projecting talent on promotion, and on P1 this bias is 0.4698; the cluster interval of the estimate 0.951 covers the probability limit 0.9698 but not the truth. Who sets the weights: pooled OLS puts every merchant-month together with equal weight, treating differences across merchants and differences within a merchant alike, and in P1 the differences in promotion across merchants come mainly from talent. Even when POLS.1 holds, $\alpha_i$ makes the composite errors of the same merchant correlated across periods, so the usual standard error is wrong: on P1 it is 0.0172, against 0.0220 when clustered by merchant.

### 3.2 RE: Feasible GLS and Quasi-Demeaning

RE gives 0.739, which lies between pooled OLS and FE. It is not an arbitrary compromise between the two but a GLS: it leaves $\alpha_i$ in the error and weights the data by the special correlation structure of the composite error.

::: {.assumption}
**RE.1 (Wooldridge 2010 §10.4.1, p. 292)** (a) FE.1; (b) $\mathbb{E}[\alpha_i \mid \mathbf{X}_i] = \mathbb{E}[\alpha_i] = 0$.

**RE.2 (p. 293)** $\mathbb{E}[\mathbf{X}_i'\Omega^{-1}\mathbf{X}_i]$ has full rank, where $\Omega = \mathbb{E}[\boldsymbol\nu_i\boldsymbol\nu_i']$ is a $T \times T$ matrix.

**RE.3 (p. 294)** (a) $\mathbb{E}[\boldsymbol\varepsilon_i\boldsymbol\varepsilon_i' \mid \mathbf{X}_i, \alpha_i] = \sigma_\varepsilon^2 I_T$; (b) $\mathbb{E}[\alpha_i^2 \mid \mathbf{X}_i] = \sigma_\alpha^2$.
:::

RE.1(b) is RE's entire wager: the individual effect is unrelated to the covariates in all periods. When it holds, $\alpha_i$ is just a harmless disturbance that has slipped into the error term. RE.3 is an auxiliary assumption that gives the variance matrix $\Omega$ of the composite error the RE structure

$$\Omega = \sigma_\varepsilon^2 I_T + \sigma_\alpha^2\boldsymbol\iota_T\boldsymbol\iota_T',$$

where $\boldsymbol\iota_T$ is the vector with every entry equal to 1: the correlation between the composite errors of the same merchant in any two periods is $\sigma_\alpha^2/(\sigma_\alpha^2 + \sigma_\varepsilon^2)$, and it does not decay as the gap between the periods grows. GLS under this structure has a simple form, quasi-demeaning: subtract part of the unit mean from each period, in the proportion

$$\varphi = 1 - \sqrt{\frac{\sigma_\varepsilon^2}{\sigma_\varepsilon^2 + T\sigma_\alpha^2}}.$$

$\varphi$ is the quasi-demeaning parameter. Wooldridge writes it as $\lambda$; in this book $\lambda_t$ is the time effect, so we write $\varphi$ instead.

::: {.intuition}
Quasi-demeaning can be seen as a knob between pooled OLS and FE. At $\varphi = 0$ nothing is subtracted, and we are back to pooled OLS; at $\varphi = 1$ the whole unit mean is subtracted, which is FE. RE turns the knob to a position set by the data: the larger the variance of the individual effect relative to the idiosyncratic error, and the more periods each unit has, the purer the "talent" component of the unit mean, and the more of it should be subtracted. When RE.1(b) and RE.3 both hold, this position makes the estimate most precise; when RE.3 fails, RE is still consistent, just no longer the most precise (P7 in Section 3.6); when RE.1(b) fails, the part of the mean left unsubtracted brings talent back into the regression.
:::

::: {.theorem}
**Theorem 2 (RE; Wooldridge 2010 §10.4.1 and §10.7.2, equations (10.34) and (10.78) to (10.81))** (i) Under RE.1 and RE.2, with $\hat\Omega$ built from estimates of the variance components (when RE.3 fails, it is enough that these converge to some limit), the feasible GLS estimator $\hat\beta_{RE}$ is consistent and asymptotically normal whether or not RE.3 holds; (ii) if RE.3 also holds, it is efficient in the class of estimators consistent under $\mathbb{E}[\boldsymbol\nu_i \mid \mathbf{X}_i] = 0$, and the usual variance estimator applies; (iii) $\hat\beta_{RE}$ equals the pooled OLS of $Y_{it} - \hat\varphi\bar Y_i$ on $(1 - \hat\varphi,\ X_{it} - \hat\varphi\bar X_i)$, where $\hat\varphi$ plugs the estimated variance components into the formula for $\varphi$.
:::

The proof sketch has three steps. Step one (RE.1): RE.1 implies $\mathbb{E}[\boldsymbol\nu_i \mid \mathbf{X}_i] = 0$, so the regressors are strictly exogenous with respect to the composite error, and the GLS argument of Chapter 9 shows that weighting by any positive definite $\hat\Omega$ then leaves consistency intact; RE.2 guarantees invertibility. Step two (RE.3): when the variance really has the RE structure, feasible GLS is asymptotically equivalent to true GLS, and GLS is efficient. Step three (the algebra of RE.3): $\Omega^{-1/2}$ is proportional to $I_T - \varphi\,T^{-1}\boldsymbol\iota_T\boldsymbol\iota_T'$, and premultiplying by it subtracts $\varphi$ times the unit mean from each period, which gives (iii).

In words: under RE.1, RE estimates $\beta$; when RE.1(b) fails, it converges to a weighted average of $\beta$ and the probability limit of the between estimator (Section 3.5), and on P1 this limit is 0.7537, against an estimate of 0.739. The weights are set by $\varphi$: $\varphi = 0$ is pooled OLS and $\varphi = 1$ is FE; $\varphi$ tends to 1 as $T$ and $\sigma_\alpha^2/\sigma_\varepsilon^2$ grow, so when $T$ is large or the variance of the individual effect is large, RE automatically moves toward FE. On P1, $\hat\varphi$ is only 0.3785. $\hat\sigma_\alpha^2$ comes from the variance of the between-regression residuals minus $\hat\sigma_\varepsilon^2/T$; when talent is correlated with promotion, $\bar X_i$ in the between regression soaks up part of $\alpha_i$, the residual variance is too small, and $\hat\sigma_\alpha^2$ is only 0.2520, far below the true value of 1; with the true variance components, $\varphi$ would be 0.6220. So when RE.1(b) fails, even the variance components that RE uses for its weights are wrong.

The practical implication of Theorem 2 (i) is that when RE.3 fails, RE is still consistent; only the usual variance is wrong. A conditional variance of the composite error that varies with $\mathbf{X}_i$, or serial correlation in the idiosyncratic error, makes RE.3 fail, and in monthly data this is almost the norm. The remedy of Wooldridge (2010) §10.4.2 is the same as for FE: treat the quasi-demeaned regression as a pooled OLS and use a variance clustered by unit. On P1, RE's classical standard error is 0.0171, and 0.0182 after clustering.

### 3.3 FE: Using Only Within-Unit Variation

FE gives 0.506, and its cluster interval [0.471, 0.541] covers the truth. Running OLS with a dummy variable for each of the 600 merchants gives exactly the same promotion coefficient, 0.506302.

::: {.assumption}
**FE.2 (Wooldridge 2010 §10.5.1, p. 304)** $\sum_t \mathbb{E}[\ddot X_{it}\ddot X_{it}']$ has full rank.

**FE.3 (Wooldridge 2010 §10.5.2, p. 304)** $\mathbb{E}[\boldsymbol\varepsilon_i\boldsymbol\varepsilon_i' \mid \mathbf{X}_i, \alpha_i] = \sigma_\varepsilon^2 I_T$.
:::

FE.2 rules out time-invariant regressors: any variable that is constant over time for every unit is identically zero after demeaning. It does not rule out interactions of such variables with time: whether the effect of certified status changes from month to month can be estimated within FE, except that the level effect in any given period cannot be estimated, only the change relative to a base period (Wooldridge 2010 equation (10.43)). FE.3 is the same as RE.3(a) and affects only the variance.

::: {.theorem}
**Theorem 3 (FE; Wooldridge 2010 §10.5, equations (10.48), (10.53), (10.56), (10.57), and (10.59))** (i) Under FE.1 and FE.2, $\hat\beta_{FE}$ (the pooled OLS of $\ddot Y_{it}$ on $\ddot X_{it}$) is unbiased, consistent, and asymptotically normal, with no restriction at all on the correlation between $\alpha_i$ and $X$; (ii) with FE.3 added, $\mathrm{Avar}\big(\sqrt N(\hat\beta_{FE} - \beta)\big) = \sigma_\varepsilon^2\big(\mathbb{E}[\ddot{\mathbf{X}}_i'\ddot{\mathbf{X}}_i]\big)^{-1}$, a consistent estimator of $\sigma_\varepsilon^2$ is $\mathrm{SSR}/[N(T-1) - k]$, and FE is efficient in the class of estimators consistent under FE.1; (iii) in an OLS on $N$ unit dummies and $X_{it}$, the coefficient on $X_{it}$ and the residuals are the same as FE's; (iv) without FE.3, the cluster-robust variance $(\ddot{\mathbf{X}}'\ddot{\mathbf{X}})^{-1}\big(\sum_i \ddot{\mathbf{X}}_i'\hat{\ddot{\boldsymbol\varepsilon}}_i\hat{\ddot{\boldsymbol\varepsilon}}_i'\ddot{\mathbf{X}}_i\big)(\ddot{\mathbf{X}}'\ddot{\mathbf{X}})^{-1}$ applies.
:::

The proof sketch has four steps. Step one (FE.1): $\mathbb{E}[\ddot\varepsilon_{it} \mid \mathbf{X}_i] = \mathbb{E}[\varepsilon_{it} \mid \mathbf{X}_i] - \mathbb{E}[\bar\varepsilon_i \mid \mathbf{X}_i] = 0$, and $\ddot X_{it}$ is a function of $\mathbf{X}_i$, so the demeaned equation satisfies strict exogeneity and pooled OLS on it is unbiased and consistent; FE.2 guarantees invertibility. Step two (FE.3): the variance of the demeaned error is $\sigma_\varepsilon^2(1 - 1/T)$ and the correlation between any two periods is $-1/(T-1)$ (equation (10.52)), but $\ddot{\mathbf{X}}_i'\ddot{\boldsymbol\varepsilon}_i = \ddot{\mathbf{X}}_i'\boldsymbol\varepsilon_i$, so the negative correlation does not enter the asymptotic variance; it only changes the degrees of freedom to $N(T-1) - k$. Step three: in a partitioned regression (the Frisch-Waugh-Lovell theorem), taking residuals on the unit dummies is demeaning, so (iii) is an algebraic identity. Step four (FE.1 and FE.2): without FE.3 the asymptotic variance is the sandwich of Section 3.4 of Chapter 5, and its estimate summed by unit is the cluster-robust variance of Section 3.4 of Chapter 8, which requires $T$ to be small relative to $N$.

In words: under FE.1, FE estimates $\beta$ however talent and promotion are correlated, and on P1 it is 0.506. When strict exogeneity fails, it converges to something else, and Section 3.4 gives that limit. The weights are set by within-unit variation:

$$\hat\beta_{FE} = \Big(\sum_i\sum_t \ddot X_{it}\ddot X_{it}'\Big)^{-1}\sum_i\sum_t \ddot X_{it}\ddot Y_{it},$$

in which each merchant's contribution is weighted by its own $\sum_t \ddot X_{it}\ddot X_{it}'$. A merchant whose promotion never changes gets zero weight, and one whose promotion swings up and down weighs the most; when everyone has the same slope this is only a matter of efficiency, but when slopes differ across merchants it decides whose average FE estimates (Section 3.8). The dummy-variable regression gives the same $\beta$, but each of the 600 merchant intercepts $\hat\alpha_i$ rests on only 6 months of data and is not consistent as $N$ grows; their average, however, is a consistent estimator of $\mathbb{E}[\alpha_i]$, and it is the FE "intercept" that software reports. One more small but common trap: if you run OLS on the demeaned data yourself, the denominator of the error variance is $NT - k$ instead of $N(T-1) - k$, so the classical standard errors you report are systematically too small, and the smaller $T$, the larger the understatement.

### 3.4 FD, Choosing Between FE and FD, and Testing Strict Exogeneity

Block P2 lets the same set of monthly shocks enter in two ways: either independent across months (iid), or accumulating into a random walk, as when a bad review or a stockout permanently changes a merchant's situation. Across 1,000 replications, the standard deviation of the FE estimate is 0.0239 with iid errors, against 0.0310 for FD; with random-walk errors it is the other way round, 0.0195 for FD and 0.0318 for FE. These are Monte Carlo sampling standard deviations, not the standard errors of a single sample.

::: {.assumption}
**FD.1 (Wooldridge 2010 §10.6.1, p. 316)** The same as FE.1.

**FD.2 (p. 317)** $\sum_{t=2}^T \mathbb{E}[\Delta X_{it}\Delta X_{it}']$ has full rank.

**FD.3 (p. 318)** With $e_{it} = \Delta\varepsilon_{it}$, $\mathbb{E}[\mathbf{e}_i\mathbf{e}_i' \mid \mathbf{X}_i, \alpha_i] = \sigma_e^2 I_{T-1}$.
:::

FD.3 amounts to saying that $\varepsilon_{it} = \varepsilon_{i,t-1} + e_{it}$ is a random walk, with $e_{it}$ mutually uncorrelated innovations. It is the opposite extreme from FE.3: FE.3 says the error has no serial correlation, and FD.3 says its serial correlation is as strong as a unit root. With only two periods, FE and FD are numerically identical (Wooldridge 2010 §10.7.1), so the choice arises only when $T > 2$.

::: {.theorem}
**Theorem 4 (FD, and the choice between FE and FD; Wooldridge 2010 §10.6 and §10.7.1, equations (10.65), (10.67), (10.71), and (10.75) to (10.77))** (i) Under FD.1 and FD.2, $\hat\beta_{FD}$ (the pooled OLS of $\Delta Y_{it}$ on $\Delta X_{it}$) is consistent; with FD.3 added, it is efficient in the class of estimators consistent under FE.1, and the usual variance applies. (ii) Under FE.3 the first-order autocorrelation of the FD residuals is $-0.5$, and under FD.3 it is 0. (iii) If there is only contemporaneous exogeneity, $\mathbb{E}[X_{it}\varepsilon_{it}] = 0$, and $(X_{it}, \varepsilon_{it})$ is stationary and weakly dependent as a time series (the idiosyncratic error has no unit root), the inconsistency of FE is $O(1/T)$, while the inconsistency of FD does not depend on $T$ and is generally not zero.
:::

The proof sketch has three steps. Step one (FD.1 and FD.2): differencing removes $\alpha_i$, and FE.1 implies $\mathbb{E}[\Delta X_{it}\Delta\varepsilon_{it}] = 0$, that is, the differenced equation satisfies POLS.1. Step two (FD.3): the differenced error $e_{it}$ is homoskedastic and uncorrelated, so the differenced equation satisfies POLS.3, and by Theorem 1 the usual variance applies; efficiency comes from the fact that the differenced error already has the spherical structure that GLS requires. If instead $\varepsilon_{it}$ is uncorrelated, $\mathrm{Corr}(e_{it}, e_{i,t-1}) = -0.5$, which is (ii). Step three (contemporaneous exogeneity and weak dependence): the asymptotic bias of FE is $\big(T^{-1}\sum_t\mathbb{E}[\ddot X_{it}\ddot X_{it}']\big)^{-1}T^{-1}\sum_t\mathbb{E}[\ddot X_{it}\varepsilon_{it}]$, whose numerator equals $-\mathbb{E}[\bar X_i\bar\varepsilon_i]$; under weak dependence $\mathrm{Var}(\bar X_i)$ and $\mathrm{Var}(\bar\varepsilon_i)$ are both $O(1/T)$, so by the Cauchy-Schwarz inequality the numerator is $O(1/T)$ (equation (10.76)). The numerator for FD is $-\big(\mathbb{E}[X_{it}\varepsilon_{i,t-1}] + \mathbb{E}[X_{i,t-1}\varepsilon_{it}]\big)$ (Wooldridge 2010 p. 323), which under stationarity does not depend on $T$; under sequential exogeneity the second term is zero, leaving only $-\mathbb{E}[X_{it}\varepsilon_{i,t-1}]$, and this is the case in P3. If the idiosyncratic error is a random walk, this advantage of FE no longer holds (p. 324).

In words: under strict exogeneity FE and FD estimate the same $\beta$, and the serial correlation of the error decides which is more precise; P2's two sets of Monte Carlo standard deviations trade places neatly at the two poles (0.0239 against 0.0310, 0.0318 against 0.0195), and real data mostly lie somewhere in between. The autocorrelation of the FD residuals can help decide: in one sample of P2 it is $-0.537$ with iid errors and $-0.026$ with a random walk; on P1 it is $-0.4906$ (cluster standard error 0.0157), compatible with the $-0.5$ of FE.3. Under both kinds of error, the Monte Carlo average of the cluster standard error is close to the true standard deviation (0.0234 for FE and 0.0310 for FD with iid errors, 0.0315 for FE and 0.0200 for FD with a random walk), so choosing the wrong estimator costs precision, not the validity of inference. When strict exogeneity fails, the two aim at different things, and the gap between them is itself an alarm. In block P3, feedback makes $\mathrm{Cov}(X_{i,t+1}, \varepsilon_{it}) = -0.6$ and the within-merchant variance of promotion is 1.36, so the probability limit of FE is $0.5 + 0.6/(6 \times 1.36) = 0.5735$ and that of FD is 0.7206; the sample estimates 0.574 and 0.724 each have an interval that covers its own limit, and neither covers 0.5. As for the weights, FE compares each period with the unit mean, while FD compares only adjacent periods: the correlation created by feedback appears only between adjacent periods, so in FE it is diluted by the mean over $T$ periods, while in FD it enters every difference undiluted.

FE and FD both bet on FE.1, but fortunately it can be tested. The most direct way comes from Wooldridge (2010) p. 325: add next period's promotion, the lead, to the FE regression

$$Y_{it} = X_{it}'\beta + X_{i,t+1}'\delta + \alpha_i + \varepsilon_{it},$$

using only the months $t = 1, \ldots, T-1$. Under strict exogeneity $\delta = 0$: next period's promotion should not explain this period's GMV. If it is significant, there are two common reasons: next month's promotion responded to this month's shock (feedback, as in P3), or customers anticipated next month's promotion and put off buying (anticipation), in which case the model has left out the dynamic effects of promotion. In either case, an FE regression with only the current month's promotion cannot be trusted. In P3, feedback makes $\hat\delta = -0.4262$, with a $p$-value that is numerically zero; in P1, $\hat\delta = -0.0342$ and $p = 0.0952$: not rejected at the 5% level, but rejected at the 10% level. P1's design has no feedback, and when the null is true a $p$-value below 0.10 has a one-in-ten chance anyway, so this result can only be read as "no strong evidence"; it falls far short of proving strict exogeneity. The autocorrelation test on the FD residuals (equation (10.71) of Wooldridge 2010 §10.6.3) tests FD.3, not strict exogeneity.

### 3.5 How RE Relates to FE: Between, Mundlak, and Chamberlain

There is a third estimator: average each merchant's 6 months into one row and run OLS of $\bar Y_i$ on $\bar X_i$; this is the between estimator. On P1 it is 1.550; it uses only variation across merchants and suffers the heaviest contamination from talent. The RE of Section 3.2 is exactly a weighted average of it and FE:

$$\hat\beta_{RE} = \omega\,\hat\beta_{BE} + (1 - \omega)\,\hat\beta_{FE},$$

$$\omega = \frac{(1-\hat\varphi)^2\,T\sum_i(\bar X_i - \bar X)^2}{\sum_i\sum_t \ddot X_{it}^2 + (1-\hat\varphi)^2\,T\sum_i(\bar X_i - \bar X)^2}.$$

This is the form for a single regressor, with $\bar X$ the full-sample mean; with several regressors $\omega$ becomes a matrix (Wooldridge 2010 p. 332 attributes this result to Hausman and Taylor 1981). In P1 the within sum of squares of promotion, $\sum_i\sum_t \ddot X_{it}^2$, is 3099.1 and the between part, $T\sum_i(\bar X_i - \bar X)^2$, is 2302.53; multiplying the latter by $(1-\hat\varphi)^2$ gives $\omega = 0.2230$, and plugging this back in reproduces RE up to a difference on the order of $10^{-15}$. This formula puts the two faces of RE side by side: when RE.1 through RE.3 hold, the larger the weight on between, the more precise RE is relative to FE, because it also uses the variation across merchants; but the larger this weight, the greater the danger, because variation across merchants is where talent does its work.

Mundlak (1978) turned this relationship into a regression. Rather than assume that $\alpha_i$ is unrelated to $X$, write the correlation out:

::: {.assumption}
**CRE.1 (Mundlak 1978; Wooldridge 2010 §10.7.3, before equation (10.88))** $\alpha_i = \xi_0 + \bar X_i'\xi + a_i$, where $a_i$ has mean zero and is uncorrelated with $\mathbf{X}_i$.
:::

Substituting into the model gives

$$Y_{it} = \xi_0 + X_{it}'\beta + \bar X_i'\xi + a_i + \varepsilon_{it}.$$

::: {.theorem}
**Theorem 5 (the Mundlak regression; Wooldridge 2010 §10.7.3, equation (10.88); for the unbalanced case of (i), see Wooldridge 2019)** Run pooled OLS or RE on $(1, X_{it}, \bar X_i)$: (i) the coefficient on $X_{it}$ is numerically identical to $\hat\beta_{FE}$, in balanced and unbalanced panels alike, with $\bar X_i$ the mean over unit $i$'s $T_i$ observed periods; (ii) in a balanced panel, $\hat\xi = \hat\beta_{BE} - \hat\beta_{FE}$; (iii) with month dummies $\lambda_t$ added, the coefficient on $X_{it}$ in a balanced panel is identical to the TWFE estimator; in an unbalanced panel the identity holds only if the unit means of the month dummies are added as well.
:::

The proof sketch has three steps. Step one (reparameterization): $X_{it} = \ddot X_{it} + \bar X_i$, so the regression on $(1, X_{it}, \bar X_i)$ and the regression on $(1, \ddot X_{it}, \bar X_i)$ span the same space; the coefficient on $X_{it}$ is the coefficient on $\ddot X_{it}$, and the coefficient on $\bar X_i$ becomes $\beta + \xi$. Step two (orthogonality): for each unit $\sum_t \ddot X_{it} = 0$, so $\ddot X_{it}$ is orthogonal in the sample to any unit-level variable, and its coefficient is the coefficient from regressing $Y_{it}$ on $\ddot X_{it}$ alone, which is FE; this uses only means over the observed periods and does not depend on whether the panel is balanced, which gives (i). Likewise, the coefficient on $\bar X_i$ comes from the regression of $Y_{it}$ on $(1, \bar X_i)$, which in a balanced panel is the between regression, so $\hat\beta + \hat\xi = \hat\beta_{BE}$, which is (ii). Step three: repeat the first two steps with the month dummies treated as regressors too; what is needed is each regressor's own unit mean. In a balanced panel the unit means of the month dummies are constants absorbed by the intercept, while in an unbalanced panel they are not, which gives (iii).

In words: under FE.1 and CRE.1, the $\beta$ that the Mundlak regression estimates is FE's $\beta$, and RE amounts to forcing $\xi = 0$. This also explains RE's efficiency (Wooldridge 2010 p. 332): $\bar X_i$ is highly correlated with $X_{it}$, and when $\xi$ really is zero, taking $\bar X_i$ out of the equation reduces collinearity; when $\xi$ is not zero, taking it out creates an omitted variable. So "RE or FE" becomes a question about $\xi$, which Section 3.6 tests with a cluster-robust test. Even if CRE.1 is wrong ($\alpha_i$ not linear in $\bar X_i$), no harm is done: identity (i) is an algebraic fact, $\hat\beta$ is still FE, and $\hat\xi$ is merely a linear projection coefficient. The weights: the weight $\omega$ that RE puts on between is set by $\varphi$ and by the ratio of between to within variation; on P1, between gets a weight of 0.2230 and itself sits far off at 1.550, so RE stops at 0.739. On P1, $\hat\xi = 1.0435$, the difference between 1.550 and 0.506, with a merchant-clustered standard error of 0.0369; the coefficient of the two-way version, 0.504978, is identical to TWFE's; on the unbalanced panel of block P6 the one-way identity holds exactly as well, while the two-way regression that leaves out the month means is off by 0.00101. This two-way Mundlak regression is what the ETWFE (extended TWFE) of Chapter 19 points back to.

Mundlak lets $\alpha_i$ be correlated with the regressors only through $\bar X_i$. Chamberlain (1982) relaxes this, writing down the linear projection of $\alpha_i$ on the regressors from all periods, $\alpha_i = \xi_0 + \sum_{t=1}^T X_{it}'\xi_t + a_i$ (the projection always exists; it is not an assumption), of which Mundlak is the special case $\xi_t \equiv \xi/T$. Substituting into the model gives

$$Y_{it} = \xi_0 + X_{i1}'\xi_1 + \cdots + X_{it}'(\xi_t + \beta) + \cdots + X_{iT}'\xi_T + a_i + \varepsilon_{it}.$$

This is the reduced form for period $t$: in the regression of $Y_{it}$ on $(1, X_{i1}, \ldots, X_{iT})$, only the coefficient on the period-$t$ regressor carries an extra $\beta$. Stack the reduced-form coefficients of the $T$ periods into $\pi$; they are linear functions of $\theta = (\xi_0, \beta, \xi_1, \ldots, \xi_T)$. These restrictions also give a more systematic test of strict exogeneity.

::: {.theorem}
**Theorem 6 (Chamberlain's approach; Wooldridge 2010 §11.1.2, equations (11.5) to (11.10); for the derivation of minimum distance, see §14.6.2)** Under FE.1 and a rank condition, the optimal MD estimator of the reduced form $\pi = h(\theta)$ is consistent and asymptotically normal; the overidentification statistic $\chi^2$ has $(T-1)(1 + Tk) - k$ degrees of freedom, and it tests precisely FE.1, without any assumption about the error variance.
:::

The proof sketch has three steps. Step one (FE.1): the projection error $a_i$ is uncorrelated with $X$ in every period, and FE.1 makes $\varepsilon_{it}$ uncorrelated with those regressors as well, so OLS of each period's outcome on all of $\mathbf{X}_i$ estimates $\pi$ consistently, and the variance computed by unit is robust to any heteroskedasticity and serial correlation. Step two (rank condition): $\pi$ has $T(1 + Tk)$ elements and $\theta$ only $1 + k + Tk$, and minimum distance pulls $\hat\pi$ back onto the restricted set $\pi = h(\theta)$, with the inverse of the variance of $\hat\pi$ as the optimal weight. Step three: the number of restrictions equals the difference between the two, and the overidentification statistic follows a $\chi^2$ with the corresponding degrees of freedom. Steps two and three are derived in Section 3.8 of Chapter 7, and only the conclusions are used here.

In words: under FE.1, Chamberlain's MD estimator and FE estimate the same $\beta$, 0.5056 on P1 (standard error 0.0174); when strict exogeneity fails, $\pi$ no longer satisfies those restrictions, and the estimator aims at a quantity contaminated by feedback, 0.5558 on P3, whose interval [0.5354, 0.5762] does not cover 0.5 either. The weights are set by MD's optimal weight matrix, the inverse of the variance of the reduced-form estimation error, which keeps its inference valid under heteroskedasticity and serial correlation. The test numbers are clear: on P1 the overidentification statistic is 30.65 with 34 degrees of freedom, $p = 0.633$; on P3 it is 1492.75, again with 34 degrees of freedom. P1's six $\hat\xi_t$ lie between 0.158 and 0.192 and sum to 1.053, almost the same as Mundlak's $\hat\xi = 1.0435$, which looks compatible with Mundlak's restriction that forces the coefficients of all periods to be equal (this chapter does not test that restriction formally); the value of Chamberlain's approach is that it need not assume this in advance.

### 3.6 The Hausman Test: Classical and Robust

With the Mundlak regression in hand, the Hausman test comes in two versions, and the focus here is the size of the test (the rejection rate when the null is true). In block P7 the null is true: talent is mean independent of promotion, so RE.1 holds, but the variance of talent changes with the level of promotion, and the monthly errors are heteroskedastic and serially correlated, so RE.3 fails. Across 2,000 replications, the classical Hausman test at a nominal 5% rejects at a rate of 0.1690, and the robust test based on the Mundlak regression at a rate of 0.0510 (figure below). If the true size is 5%, the Monte Carlo standard error of a rejection rate over 2,000 replications is 0.0049.

![Monte Carlo for block P7 (2,000 samples, null true: RE.1 holds but RE.3 fails): actual rejection rates of the classical Hausman test and of the merchant-clustered robust test based on the Mundlak regression, at nominal levels from 1% to 20%. The black line is where the rejection rate equals the nominal level, and the gray band is the Monte Carlo 95% error band around it.](assets/fig/fig_10_hausman_size.svg)

::: {.intuition}
The Hausman test is like cross-examining two witnesses. FE is a cautious witness whose testimony is credible whether or not talent and promotion are correlated, only not very precise; RE is a precise but partisan witness who tells the truth only when talent and promotion are unrelated. If their testimonies differ by more than sampling noise, the second witness's partisanship is the problem. To judge how much counts as "more than," we need to know how much the difference between the two testimonies jitters, and the classical formula computes this with a shortcut: if RE is the most efficient witness, the variance of the difference is exactly the difference of their variances. Once the premise of the shortcut fails, the computed jitter is wrong, and the test's ruler is bent.
:::

::: {.theorem}
**Theorem 7 (the Hausman test; Wooldridge 2010 §10.7.3, equations (10.82) to (10.88), Lemma 14.1; Hausman 1978, Arellano 1993)** Compare only the coefficients on regressors that vary across both units and time, $k$ of them; regressors that vary only over time, such as month dummies, cannot enter the comparison, and $\bar X_i$ must not be perfectly collinear with the intercept and other time-constant regressors, or the variance matrix of the difference is singular. The $\mathrm{Avar}(\hat\beta)$ below is shorthand for $\mathrm{Avar}\big(\sqrt N(\hat\beta - \beta)\big)/N$, that is, the form in Theorem 3 (ii) divided by $N$. (i) Under RE.1 through RE.3 and FE.2,

$$\begin{aligned}
H &= (\hat\beta_{FE} - \hat\beta_{RE})' \\
&\quad \times \big[\widehat{\mathrm{Avar}}(\hat\beta_{FE}) - \widehat{\mathrm{Avar}}(\hat\beta_{RE})\big]^{-1} \\
&\quad \times (\hat\beta_{FE} - \hat\beta_{RE}) \xrightarrow{d} \chi^2_k;
\end{aligned}$$

\(ii) under RE.1 and FE.2, whether or not RE.3 holds, the cluster-robust Wald statistic for $\xi = 0$ in the Mundlak regression converges in distribution to $\chi^2_k$. Both test that $\bar X_i$ and $\alpha_i$ are uncorrelated, and both presuppose FE.1.
:::

The proof sketch has four steps. Step one (RE.1): both estimators are consistent, their difference converges to zero, and $\sqrt N$ times the difference is asymptotically normal. Step two (RE.1 through RE.3): RE is efficient in the class of estimators that are consistent under $\mathbb{E}[\boldsymbol\nu_i \mid \mathbf{X}_i] = 0$ (Theorem 2 (ii)), a class that includes FE, and $\hat\beta_{RE} + c(\hat\beta_{FE} - \hat\beta_{RE})$ is consistent for any $c$, with its variance minimized at $c = 0$, which requires $\mathrm{Cov}(\hat\beta_{RE}, \hat\beta_{FE} - \hat\beta_{RE}) = 0$. Step three (step two): hence $\mathrm{Avar}(\hat\beta_{FE} - \hat\beta_{RE}) = \mathrm{Avar}(\hat\beta_{FE}) - \mathrm{Avar}(\hat\beta_{RE})$, and the quadratic form gives (i). Step four (RE.1(b) and Theorem 5): RE.1(b) implies $\xi = 0$, the cluster-robust variance of the Mundlak regression does not depend on RE.3 (Section 3.4 of Chapter 8), and the Wald statistic gives (ii).

In words: the Hausman test is aimed not at "whether RE is right" but at the moment condition of equation (10.87), namely whether the unit mean is correlated with the individual effect; when strict exogeneity fails, FE and RE are both biased, and the test loses this meaning. When RE.3 fails, the covariance in step two is generally not zero, and the two variances in the classical formula are no longer the true asymptotic variances either, so $\mathrm{Avar}(\hat\beta_{FE} - \hat\beta_{RE}) \neq \mathrm{Avar}(\hat\beta_{FE}) - \mathrm{Avar}(\hat\beta_{RE})$: the limit of the classical statistic is not $\chi^2_k$, and whether it is too large or too small cannot be told in advance. On P7, across 2,000 replications, the standard deviation of the difference between FE and RE is 0.0105, while the classical formula gives only 0.0073 on average, a ratio of 1.45, so the statistic is inflated by roughly the square of that ratio: the classical $H$ averages 2.09 and the robust statistic averages 0.98, against a mean of 1 for $\chi^2_1$. The weights are set by the auxiliary assumption: the classical version quietly puts RE.3 into the null, while the robust version requires only RE.1. On P1 both clearly reject RE: 3362.02 for the classical test and 798.95 for the robust one, both with negligible $p$-values. When the null fails, both statistics diverge, and since they use different variance estimators, their values need not be close.

Once time-invariant controls are put into RE, equation (10.87) tests whether the part of $\bar X_i$ that the controls cannot explain is correlated with $\alpha_i$: with rich enough controls, RE can be close to FE even if $\bar X_i$ itself is correlated with talent (Wooldridge 2010 p. 331).

### 3.7 Hausman-Taylor: Time-Invariant Regressors

In block P4, "Nimbus Certified merchant" status $W_i$ raises monthly GMV by USD 1,500. FE's demeaning wipes it out together with $\alpha_i$; RE can estimate it, but more talented merchants apply for certification more often, so RE.1(b) fails, and RE gives 2.913. The idea of Hausman and Taylor (1981) is to find instruments inside the model: some time-varying regressors are unrelated to $\alpha_i$, and their unit means can serve as instruments for $W_i$. In P4, the share of orders shipped through Nimbus fulfillment, $X_{1it}$, depends mainly on how far the merchant is from a warehouse and is unrelated to talent, while a record of Nimbus fulfillment is in turn a prerequisite for applying for certification; promotion, $X_{2it}$, is correlated with talent. The group subscripts 1 and 2 come first and the unit and period last: $X_{1it}$ is the first group of regressors for unit $i$ in period $t$.

::: {.assumption}
**HT.1 (Wooldridge 2010 equation (11.30))** In the model $Y_{it} = X_{1it}'\beta_1 + X_{2it}'\beta_2 + W_i'\gamma + \alpha_i + \varepsilon_{it}$, $\mathbb{E}[\varepsilon_{it} \mid W_i, \mathbf{X}_i, \alpha_i] = 0$.

**HT.2 (equation (11.32))** $W_i = (W_{1i}', W_{2i}')'$, $\mathbb{E}[W_{1i}\alpha_i] = 0$, and $\mathbb{E}[X_{1it}\alpha_i] = 0$; $X_{2it}$ and $W_{2i}$ may be arbitrarily correlated with $\alpha_i$.

**HT.3 (Wooldridge 2010 p. 361, the order condition when $\bar X_{1i}$ is the instrument, plus a rank condition)** The number of $X_{1it}$ is at least the number of $W_{2i}$, and $\bar X_{1i}$ has explanatory power for $W_{2i}$.
:::

Two extreme cases help explain HT: when all regressors are unrelated to $\alpha_i$, HT reduces to RE; when there is no exogenous time-varying regressor at all, $\gamma$ cannot be identified, and we are back to FE. HT.2 is the core of identification, and it is an exclusion restriction that cannot be tested on its own. The general order condition after equation (11.34) uses the $X_{1it}$ of all periods as instruments and requires only that the number of $X_{1it}$ times $T$ be at least the number of $W_{2i}$; this chapter uses only the unit means, so the condition is correspondingly stricter. On P4 the rank condition of HT.3 shows up in the first stage: in the unit-level regression of $W_i$ on $\bar X_{1i}$ the coefficient is 0.2656 and the robust $F$ statistic is 1566.1.

::: {.intuition}
Why can a merchant's average fulfillment share serve as an instrument for certified status? FE has already used the month-to-month swings in the fulfillment share to estimate its own coefficient, but the long-run differences between merchants go unused: merchants close to a warehouse have a persistently higher fulfillment share and also find it easier to get certified, while distance to a warehouse has nothing to do with operating talent. This cross-merchant variation thus acts like a natural experiment, isolating the part of the difference between certified and uncertified merchants that is "unrelated to talent." The instrument's persuasiveness comes entirely from that phrase, "unrelated to talent," and the data themselves cannot vouch for it.
:::

::: {.theorem}
**Theorem 8 (Hausman-Taylor; Wooldridge 2010 §11.3, equations (11.29) to (11.34); for estimation, see the REIV (random effects IV, that is, IV on the quasi-demeaned equation) of §11.2)** Under HT.1 through HT.3, IV on the quasi-demeaned equation with $(\ddot X_{1it}, \ddot X_{2it}, W_{1i}, \bar X_{1i})$ as instruments gives consistent and asymptotically normal estimates of $(\beta_1, \beta_2, \gamma)$.
:::

The proof sketch has three steps. Step one (HT.1): within each unit $\sum_t \ddot X_{it} = 0$, so $\sum_t \ddot X_{it}\alpha_i = 0$, and the correlation between the demeaned regressors and the composite error reduces to $\sum_t \ddot X_{it}\varepsilon_{it}$, whose expectation is zero by strict exogeneity; they can serve as their own instruments and identify $\beta$. Step two (HT.2): $W_{1i}$ and $\bar X_{1i}$ are both uncorrelated with $\alpha_i + \varepsilon_{it}$ and are valid instruments; HT.3 makes them sufficient to identify $\gamma$. Step three: quasi-demeaning changes only the weights, not the moment conditions, and the IV theory of Chapter 9 gives consistency and asymptotic normality.

In words: under HT.1 through HT.3, Hausman-Taylor estimates the structural parameter $\gamma$, 1.449 on P4 (classical standard error 0.197, interval [1.063, 1.836]), or about USD 1,449; if the exclusion restriction is wrong, say because the fulfillment share also reflects talent, the estimate aims at a quantity contaminated by talent, and the data cannot reveal this directly. The weights come from the instruments: $\gamma$ is identified only by the part of the variation in $\bar X_{1i}$ that predicts certification, and the weaker the first stage, the less stable the estimate. On the same data RE gives 2.913 (interval [2.780, 3.045] with classical standard errors, [2.743, 3.083] when clustered by merchant) and pooled OLS gives 2.8257; both count talent as part of certification. P4 is exactly identified (one $X_{1it}$ and one $W_{2i}$), so HT's $\hat\beta_2 = 0.4972$ is numerically the same as FE's, while RE's 0.7781 is pushed up by talent.

### 3.8 Random Slopes: Whose Average Does FE Estimate?

In block P5 each merchant has its own promotion slope $b_i$, with an average of 0.5; the merchants whose promotion fluctuates most are also the ones for whom it works best (shops that run flash sales, for example). FE gives 0.683, and its interval [0.658, 0.709] does not cover 0.5. FE has not got it wrong; it is estimating a different average.

::: {.assumption}
**RS.1 (Wooldridge 2010 equations (11.82) and (11.83))** $Y_{it} = X_{it}'b_i + \alpha_i + \varepsilon_{it}$, $\mathbb{E}[\varepsilon_{it} \mid \mathbf{X}_i, \alpha_i, b_i] = 0$, and FE.2 holds.

**RS.2 (equation (11.86))** $\mathbb{E}[b_i \mid \ddot X_{it}] = \mathbb{E}[b_i] = \beta$, $t = 1, \ldots, T$.
:::

RS.2 is weaker than it looks. Split promotion into a merchant's long-run level and its month-to-month swings: RS.2 allows the slope to be arbitrarily correlated with the long-run level (large merchants may be better at using promotion) and rules out only a correlation between the slope and the swings (Wooldridge 2010 p. 383). P5 violates the latter: merchants with larger swings also have steeper slopes.

::: {.theorem}
**Theorem 9 (FE under random slopes; Wooldridge 2010 §11.7.3, equations (11.84) to (11.86); Wooldridge 2005a)** (i) Under RS.1,

$$\begin{aligned}
\mathrm{plim}\,\hat\beta_{FE} &= \beta_* \\
&= \Big(\sum_t\mathbb{E}[\ddot X_{it}\ddot X_{it}']\Big)^{-1} \\
&\quad \times \sum_t\mathbb{E}[\ddot X_{it}\ddot X_{it}'\,b_i],
\end{aligned}$$

and with a single regressor $\beta_* = \mathbb{E}\big[b_i\sum_t\ddot X_{it}^2\big]\big/\mathbb{E}\big[\sum_t\ddot X_{it}^2\big]$, the average slope weighted by within-unit variance; (ii) if RS.2 also holds, $\beta_* = \beta$.
:::

The proof sketch has three steps. Step one (RS.1): write $b_i = \beta + (b_i - \beta)$; after demeaning, $\ddot Y_{it} = \ddot X_{it}'b_i + \ddot\varepsilon_{it}$, with $\ddot\varepsilon_{it}$ uncorrelated with $\ddot X_{it}$. Step two: substituting into the FE formula, the numerator converges to $\sum_t\mathbb{E}[\ddot X_{it}\ddot X_{it}'b_i]$, which gives (i). Step three (RS.2): by iterated expectations, $\mathbb{E}[\ddot X_{it}\ddot X_{it}'b_i] = \mathbb{E}[\ddot X_{it}\ddot X_{it}'\,\mathbb{E}(b_i \mid \ddot X_{it})] = \mathbb{E}[\ddot X_{it}\ddot X_{it}']\beta$, which gives (ii).

In words: under RS.1, FE estimates $\beta_*$, which equals the average slope only when RS.2 holds. The weights are the within-unit variances: the population FE weights of P5's four types of merchant (the figure shows the sample weights; see Section 6.5) are 3.3%, 13.3%, 30.0%, and 53.3%, although each type makes up a quarter of merchants; the type with the largest fluctuations has the steepest slope, so $\beta_* = 0.66125$, and the interval of the sample estimate 0.683 covers it. In P5b the slope varies only with the level of promotion and is unrelated to within-unit fluctuations; FE gives 0.502, and its interval [0.470, 0.534] covers the average slope. If management asks "What happens if every merchant spends USD 100 more?", P5's FE overstates the answer by more than thirty percent; if the question is "What happens if the merchants who already promote often promote a little more?", FE's weights actually fit the question better. Estimators are not right or wrong in themselves; they either match the question or do not.

The test of §11.7.4 adds interactions of $X_{it}$ with unit-level variables to the FE regression: with the $\bar X_i$ interaction, the $p$-value is numerically zero on P5b but 0.835 on P5, where the test sees no problem; with the within-unit variance interaction instead, the $p$-value on P5 is $3.7 \times 10^{-7}$. The two tests answer different questions: the $\bar X_i$ interaction tests whether the slope varies with the level of promotion, and it rejects on P5b, but this kind of heterogeneity satisfies RS.2, and FE still aims at the average slope; what matters for FE is whether the slope varies with within-unit fluctuations, which only the within-unit variance interaction can detect, and P5 is that case. Wooldridge (2010) p. 386 also warns that the $\bar X_i$ interaction test cannot detect slope heterogeneity unrelated to $\bar X_i$ and is not in itself designed to judge whether FE is consistent. The mean group estimator, which runs a regression merchant by merchant and then averages, aims at the average slope; on P5 it is 0.545 with interval [0.5075, 0.5825], which just misses 0.5 (the average of the true slopes in the sample is 0.5021); a nominal 95% interval is supposed to miss about once in twenty anyway.

![Blocks P5 (left) and P5b (right): the horizontal axis is the within-unit standard deviation of a merchant's monthly promotion spending (four types, USD 50, 100, 150, and 200, a quarter each); gray dots are the true slopes $b_i$ of one merchant in every four, blue dots are the sample mean slopes of each type, with area proportional to the type's FE weight in the sample, and the weights are printed below; the dashed line is the average slope $\mathbb{E}[b_i]$, the solid line is the probability limit of FE, and the blue band is the 95% interval of FE. The vertical axis is in thousands of USD of GMV per hundred USD of promotion.](assets/fig/fig_10_slopes.svg)

### 3.9 Unbalanced Panels and Attrition

Block P6 has 4,000 merchants and data for the first 8 months, with two attrition rules applied to the same balanced panel. Under one, exit depends on talent: a merchant whose talent is one standard deviation below average ($\alpha_i = -1.5$) leaves with probability 18.4% each month, one a standard deviation above average with probability only 1.8%, and in the end 75.4% of the merchant-months remain. Under the other, exit depends on the shock: a merchant leaves the platform in the month its GMV falls below USD 9,500 and reports nothing for that month, which leaves 72.3%.

::: {.assumption}
**Assumption 19.8 (Wooldridge 2010 §19.9.1, p. 829)** Let $s_{it} = 1$ indicate that unit $i$ is observed in period $t$, and $T_i = \sum_{t=1}^T s_{it}$. (a) $\mathbb{E}[\varepsilon_{it} \mid \mathbf{X}_i, \mathbf{s}_i, \alpha_i] = 0$; (b) $\sum_t\mathbb{E}[s_{it}\ddot X_{it}\ddot X_{it}']$ has full rank; (c) $\mathbb{E}[\boldsymbol\varepsilon_i\boldsymbol\varepsilon_i' \mid \mathbf{X}_i, \mathbf{s}_i, \alpha_i] = \sigma_\varepsilon^2 I_T$.
:::

In an unbalanced panel the unit mean is taken over the observed periods only: $\bar X_i = T_i^{-1}\sum_t s_{it}X_{it}$. Condition (a) allows selection to be arbitrarily correlated with $\alpha_i$ and $\mathbf{X}_i$, and rules out only selection that is correlated with the idiosyncratic error.

::: {.theorem}
**Theorem 10 (unbalanced panels; Wooldridge 2010 §19.9.1, equations (19.96) to (19.99))** (i) Under Assumption 19.8(a)(b), FE on an unbalanced panel is consistent and asymptotically normal, and the cluster-robust variance applies; with (c) added, the usual variance applies as well; units with $T_i = 1$ contribute no information. (ii) RE further requires $\mathbb{E}[\alpha_i \mid \mathbf{X}_i, \mathbf{s}_i] = \mathbb{E}[\alpha_i]$, and the quasi-demeaning parameter varies by unit, $\varphi_i = 1 - \sqrt{\sigma_\varepsilon^2/(\sigma_\varepsilon^2 + T_i\sigma_\alpha^2)}$.
:::

The proof sketch has two steps. Step one (19.8(a)): $s_{it}\ddot X_{it}$ is a function of $(\mathbf{X}_i, \mathbf{s}_i)$, so iterated expectations give $\mathbb{E}[s_{it}\ddot X_{it}\varepsilon_{it}] = 0$, and 19.8(b) guarantees invertibility, which is (i). Step two: RE's composite error contains $\alpha_i$, and selection that is correlated with $\alpha_i$, even if it is unrelated to $\varepsilon$, breaks $\mathbb{E}[\boldsymbol\nu_i \mid \mathbf{X}_i, \mathbf{s}_i] = 0$, so an extra condition is needed; when $T_i$ differs across units, so does the GLS transformation, which is (ii).

In words: under Assumption 19.8, FE on an unbalanced panel still estimates $\beta$. On the P6 panel with exit driven by talent, FE is 0.4916, and its interval [0.4713, 0.5119] covers 0.5; with exit driven by the shock, what remains are the months with favorable shocks, selection is correlated with $\varepsilon_{it}$, and FE is 0.4437, with an interval [0.4241, 0.4632] that does not cover 0.5; in large samples it converges to about 0.4435, so the bias does not vanish with the sample size. The weights: in an unbalanced panel each merchant is weighted by the number of months in which it is observed and by its within-unit variation, and those who remain tend to be the more talented merchants (average talent $-0.934$ among those who leave, 0.702 among those who stay), which does not affect FE but does affect RE: on the same panel RE gives 0.8669, with $\varphi_i$ ranging from 0.095 to 0.398 as $T_i$ varies; RE is so far off not mainly because the panel is unbalanced but because talent in P6 is correlated with promotion to begin with, so RE.1(b) fails, and the robust Hausman statistic from the Mundlak regression is as high as 5149.8. The Mundlak identity holds exactly here too, provided that $\bar X_i$ is averaged over the observed periods only.

Two simple diagnostics come from Wooldridge (2010) §19.9. Under Assumption 19.8, FE is also consistent on any balanced subsample, so estimate FE once on the unbalanced panel and once on its balanced part: a large difference is an alarm. The second is to add to the FE regression the indicator $s_{i,t+1}$ of whether the unit is still observed next period: if it is significant, staying in the sample is related to the current shock.

## 4 Estimation and Practice

### 4.1 The Estimators at a Glance

| Estimator | Assumption on $\alpha_i$ | Assumption on $\varepsilon$ | What it estimates | Estimate on P1 |
|---|---|---|---|---|
| Pooled OLS | Uncorrelated with $X_{it}$ | Contemporaneous exogeneity | Slope plus the projection of talent | 0.951 |
| Between | Uncorrelated with $\bar X_i$ | Strict exogeneity | Comparison across units | 1.550 |
| RE | Mean independent of all of $\mathbf{X}_i$ | Strict exogeneity; RE.3 matters only for efficiency | Weighted average of between and FE | 0.739 |
| FE | Arbitrary | Strict exogeneity | Slope weighted by within-unit variance | 0.506 |
| FD | Arbitrary | Strict exogeneity | Slope of changes between adjacent periods | 0.518 |
| Mundlak (CRE) | Linear projection on $\bar X_i$ | Strict exogeneity | $\beta$ equal to FE, plus $\xi$ | 0.506 |
| Chamberlain MD | Linear projection on all of $\mathbf{X}_i$ | Strict exogeneity, testable | Same as FE | 0.5056 |
| Hausman-Taylor | Some regressors uncorrelated with it | Strict exogeneity | Also estimates $\gamma$ | Block P4 |

The choice depends first on two questions, "Is talent correlated with the regressors?" and "Are the regressors strictly exogenous?", not on which regression happens to run. In observational data the default answer to the first question is that they are correlated, so the default estimator is FE or FD; to estimate a time-invariant variable, use Hausman-Taylor or CRE and defend the exclusion restrictions they involve. Models with a lagged dependent variable on the right-hand side are not in this table; their problems are taken up in Chapter 11.

### 4.2 Variances, Tests, and Software

Here are the two kinds of standard error for four estimators, side by side on P1:

| Estimator | Classical standard error valid under | P1 classical standard error | Cluster standard error valid under | P1 cluster standard error |
|---|---|---|---|---|
| Pooled OLS | POLS.1 through POLS.3 ($\alpha_i$ makes them almost impossible to satisfy) | 0.0172 | POLS.1 and POLS.2 | 0.0220 |
| RE | RE.1 through RE.3 | 0.0171 | RE.1 and RE.2 | 0.0182 |
| FE | FE.1 through FE.3 | 0.0175 | FE.1 and FE.2 | 0.0179 |
| FD | FD.1 through FD.3 | 0.0174 | FD.1 and FD.2 | 0.0212 |

P1's idiosyncratic errors are iid, so FE.3 holds and FD.3 fails, and the classical standard error of FD is therefore about a fifth too small; the classical standard error of pooled OLS ignores the within-unit correlation created by $\alpha_i$. Unless stated otherwise, the standard errors in this chapter are clustered by merchant, without a small-sample factor (CR0); Section 3.4 of Chapter 8 explains why the default choice for panels is clustering rather than heteroskedasticity-robust standard errors, and Section 4.2 of Chapter 8 covers the problems with few clusters. For the Hausman test, prefer the robust version based on the Mundlak regression: it is robust to heteroskedasticity and serial correlation, can include time effects and time-invariant controls, and can test just a subset of the coefficients. For strict exogeneity, use the lead test or Chamberlain's overidentification test; for selective attrition, use the balanced-subsample comparison and the selection-indicator test of Section 3.9.

Put in order, these checks give roughly the following steps for a panel analysis:

1. Report FE and FD together, with standard errors clustered by unit, and see whether the two differ by more than sampling error.
2. Run the lead test and, when $T$ is not large, Chamberlain's overidentification test as well, to judge whether strict exogeneity is credible.
3. Estimate the Mundlak regression and report $\hat\xi$ with its robust Wald test; with time effects, use the two-way version.
4. When slopes may vary across units, add interactions of $X_{it}$ with the unit mean and with the within-unit variance (to judge whether FE aims at the average slope, look at the within-unit variance interaction), and say which units carry FE's weight.
5. When the panel is unbalanced, look at the distribution of $T_i$, compare FE on the balanced subsample, and test next period's selection indicator.

In R, the package `plm` 2.6-7 provides `plm()`, which covers the chapter's main estimators: `model` takes `"pooling"`, `"within"`, `"random"`, `"fd"`, or `"between"`, and `effect = "twoways"` gives TWFE; `phtest(fe, re)` is the classical Hausman test, and `phtest(formula, data, method = "aux", vcov = ...)` adds $\bar X_i$ to the RE regression and runs a cluster Wald test, which on P1 gives 798.95, the same as this chapter's hand-coded robust test; `pwartest` and `pwfdtest` test for serial correlation in the FE and FD errors; Hausman-Taylor is written `plm(y ~ x1 + x2 + w | x1 | x2, model = "random", random.method = "ht", inst.method = "baltagi")`; and `pvcm(model = "within")` is mean group. `fixest::feols(y ~ x | id + t, cluster = ~id)` is faster with many fixed effects, and `sandwich::vcovCL` and `clubSandwich::vcovCR` provide cluster variances, the latter with small-sample corrections such as CR2. The Stata counterparts are `xtreg, fe`, `xtreg, re`, `hausman`, and `xthtaylor`.

### 4.3 From FE to TWFE: The Bridge to DiD

The model with month effects $\lambda_t$, $Y_{it} = X_{it}'\beta + \alpha_i + \lambda_t + \varepsilon_{it}$, is estimated by TWFE, which removes the unit means and the month means at the same time. Theorem 5 (iii) shows that this is equivalent to the two-way Mundlak regression

$$Y_{it} = \xi_0 + X_{it}'\beta + \bar X_i'\xi + \lambda_t + a_i + \varepsilon_{it},$$

where $\bar X_i = T_i^{-1}\sum_t s_{it}X_{it}$ averages over the observed periods only; it is estimated by pooled OLS, and in an unbalanced panel the unit means of the month dummies must be added too. In a balanced panel, the month dummies can also be replaced by the cross-sectional means of the regressors in each month without changing the coefficient on $X_{it}$, and this is the form Chapter 19 uses. When $X_{it}$ is a treatment indicator, TWFE is the most common DiD regression. When treatment starts at different times for different units, the TWFE coefficient is a weighted average of many comparisons, some of them even with negative weights, and Imai and Kim (2021) prove that it generally does not equal multi-period DiD; Chapter 19 takes this up again, returning to this chapter's two-way Mundlak regression to build ETWFE, which lets the treatment effect vary freely by cohort and by period.

## 5 Anchoring Papers

Each of the three anchoring papers corresponds to one thread of the chapter: using a panel to tackle the classic IS problem that "unobservable differences between firms are correlated with their inputs," finding instruments inside the model to estimate time-invariant effects, and putting this toolkit to work in an IS study of a platform. Each paper is laid out in five parts: paper, method, data, results, and limitations.

### 5.1 Tambe and Hitt (2012)

::: {.case}
Paper: "The Productivity of Information Technology Investments: New Evidence from IT Labor Data," Information Systems Research 23(3, Part 1 of 2):599-617. The question is how much output IT investment actually brings firms.

Method: on a firm-level panel, the authors estimate the marginal product of IT inputs with production-function specifications. The core difficulty is this chapter's $\alpha_i$ problem, plus shocks that vary over time: high-productivity firms both invest more in IT and produce more output. The authors handle the endogeneity of IT spending with GMM-based estimators.

Data: newly collected IT labor data measure firms' IT inputs and form a firm-level panel whose sample includes both Fortune 500 firms and midsize firms.

Results: after correcting for endogeneity, the IT coefficient is only about 10% lower than the unadjusted estimate, implying that the endogeneity bias in past studies may not be large; the IT return is significantly lower for midsize firms than for Fortune 500 firms; in large firms the contribution of IT is realized more slowly, while for midsize firms the short-run and long-run contributions are similar; the marginal product of IT spending from 2000 to 2006 exceeds that of any earlier period; and IT productivity is higher in manufacturing.

Limitations: the closeness of the coefficients from the two sets of estimators, like a Hausman test that does not reject, shows only that the answers are similar under their respective identifying assumptions, and does not prove that IT inputs are exogenous; the validity of the GMM instruments has to be defended with the tests that Chapter 11 will cover; and the production function also depends on how inputs are measured and on functional form.
:::

What this paper brings to the chapter is that it puts "comparing estimators under different identifying assumptions" back into a real IS question: IT investment is not scattered at random, and the difficulty is that it is correlated with unseen differences between firms.

### 5.2 Hausman and Taylor (1981)

::: {.case}
Paper: "Panel Data and Unobservable Individual Effects," Econometrica. The question is whether the effect of a time-invariant variable can be estimated when the individual effect is correlated with some of the regressors. The description of the application below follows the account in Wooldridge (2010) §11.3 (p. 361).

Method: the regressors are split into a group unrelated to the individual effect and a group related to it, each with time-varying and time-invariant members; the demeaned time-varying regressors serve as their own instruments, and the unit means of the time-varying regressors unrelated to the individual effect serve as instruments for the correlated time-invariant regressors (Theorem 8 of this chapter); estimation is by REIV.

Data: an individual-level wage panel with two periods, in which years of schooling do not change between the periods. The time-varying regressors unrelated to the individual effect are work experience, an indicator of poor health, and an indicator of unemployment in the previous year; race and union status serve as time-invariant regressors unrelated to the individual effect; and the only time-invariant regressor correlated with the individual effect is years of schooling.

Results: the estimated return to schooling is almost twice the pooled OLS estimate; once some time-varying regressors are allowed to be correlated with the individual effect, the estimate is larger still.

Limitations: Wooldridge points out that the identifying assumptions are not especially convincing, and that assuming that work experience and union status are unrelated to unobserved ability, and then using this to identify the return to schooling, is quite a stretch; an overidentification test can check only the surplus exclusion restrictions.
:::

The paper turned a shortcoming of FE, that it cannot estimate time-invariant variables, into an instrumental-variables problem, at the cost of a set of exclusion restrictions that must be defended one by one.

### 5.3 Susarla, Oh and Tan (2012)

::: {.case}
Paper: "Social Networks and the Diffusion of User-Generated Content: Evidence from YouTube," Information Systems Research. The question is how social networks shape the diffusion of user-generated content. The data and results below follow its July 2008 working-paper version (Oh, Susarla and Tan 2008).

Method: a video-level panel of repeated observations, with the growth in view counts as the dependent variable. Users may self-select into different subscription networks according to taste, and there may be unobservable associations between users and videos that are correlated with the regressors. The authors first apply $\rho$-differencing to remove serial correlation, then use exclusion restrictions to handle self-selection into subscription networks, and finally estimate by Hausman-Taylor: exogenous time-varying variables instrument the endogenous time-varying variables, and the means of the exogenous time-varying variables instrument the endogenous time-invariant variables, the same structure as Theorem 8 of this chapter.

Data: 4,106 videos uploaded by 913 users in a music and blogging interest group on YouTube, observed every 5 days over 2 months, 11 observation points in all, with friends, subscribers, and subscriptions recorded as well.

Results: the friend network matters more in the later stages of a video's diffusion, while centrality in the subscription network dominates in the early stages; ties to friends outside the group have about three times the effect on diffusion of ties within the group.

Limitations: a sample of one group over two months limits extrapolation; Hausman-Taylor's identification rests on which variables are unrelated to the associations between videos and users, and these exclusion restrictions can only be defended, not fully tested; the specifications of the working paper and the published version may differ.
:::

Taken together: Tambe and Hitt judge the size of the bias from firm heterogeneity by how close the estimators are under different identifying assumptions; Hausman and Taylor give a way to estimate time-invariant effects within the FE framework; and Susarla, Oh, and Tan apply that method to platform data, inheriting all of its exclusion restrictions as well.

## 6 A Full Walkthrough on the Nimbus Data

The code below runs under R 4.6.1, with the random seed fixed block by block (P1 keeps 737); every number in the text comes from the output of an actual run, and the result lines in the code that start with `#` copy what the run printed. Unless stated otherwise, standard errors are clustered by merchant (CR0).

### 6.1 The Data-Generating Process

The static block P1 makes the case for FE: $Y_{it} = \beta X_{it} + \alpha_i + \varepsilon_{it}$ with $\beta = 0.5$, and merchant talent $\alpha_i$ positively correlated with promotion $X_{it}$, which biases pooled OLS and RE and leaves FE unbiased.

```r
set.seed(737)

## ---- static panel: FE vs RE ------------------------------------------------
gen_static <- function(N = 600, T = 6, beta = 0.5, ca = 0.7) {
  alpha <- rnorm(N, 0, 1)                              # merchant quality (unobs.)
  id <- rep(1:N, each = T); tt <- rep(1:T, times = N)
  a  <- alpha[id]
  x  <- ca * a + rnorm(N*T, 0, 1)                      # promotions, corr with alpha
  y  <- beta * x + a + rnorm(N*T, 0, 1)
  data.frame(id, t = tt, y, x)
}
static <- gen_static()
```

The true value $\beta = 0.5$ is the target for all the P1 estimates. The parameters of the new blocks and three of their generating functions follow. The three that are not listed work as follows: in P2, promotion is AR(1); in P4, certification is determined jointly by distance to the warehouse, talent, and a random term, while the fulfillment share depends only on distance to the warehouse; in P7, the standard deviation of talent grows with $\lvert\bar X_i - 6\rvert$, and the error is a heteroskedastic AR(1).

```r
P2 <- list(N = 500, T = 6, alpha_sd = 1.5, x_mean = 6, x_alpha = 0.4, rho_x = 0.5, beta = 0.5, mu0 = 10,
           sigma_e = 1, R_mc = 1000)
P3 <- list(N = 1000, T = 6, alpha_sd = 1.5, x_mean = 6, x_alpha = 0.4, feedback = -0.6, x_sd = 1, beta = 0.5,
           mu0 = 10, sigma_e = 1)
P4 <- list(N = 2000, T = 6, alpha_sd = 2, x1_mean = 5, x1_mu = 1.0, x1_sd = 0.6, x2_mean = 6, x2_alpha = 0.3,
           x2_sd = 1, w_mu = 1.2, w_alpha = 0.35, w_cut = 0.3, beta1 = 0.3, beta2 = 0.5, gamma = 1.5, mu0 = 8,
           sigma_e = 1.5)
P5 <- list(N = 800, T = 8, alpha_sd = 1.5, mu_mean = 10, mu_alpha = 0.4, mu_sd = 0.5, s_types = c(0.5, 1, 1.5, 2),
           s_probs = rep(0.25, 4), avg_slope = 0.5, kappa = 0.15, eta_sd = 0.1, design = "within_sd", mu0 = 8,
           sigma_e = 1)
P5b <- modifyList(P5, list(kappa = 0.3, design = "level"))
P6 <- list(N = 4000, T = 8, alpha_sd = 1.5, x_mean = 6, x_alpha = 0.5, beta = 0.5, mu0 = 9, sigma_e = 1.5,
           haz0 = -1.5, haz_alpha = -0.4, floor = 9.5)
P7 <- list(N = 1000, T = 5, x_mean = 6, mu_sd = 1, rho_x = 0.6, alpha_base = 0.3, alpha_slope = 1.0,
           eps_base = 0.5, eps_slope = 0.5, rho_e = 0.6, beta = 0.5, mu0 = 22, R_mc = 2000)
.panel_index <- function(N, T) list(id = rep(seq_len(N), each = T), t = rep(seq_len(T), times = N))
gen_P3 <- function(P) {
  N <- P$N; T <- P$T; ix <- .panel_index(N, T); id <- ix$id
  alpha <- rnorm(N, 0, P$alpha_sd)
  e <- matrix(rnorm(N * (T + 1), 0, P$sigma_e), N, T + 1)          # column 1 is month 0
  u <- matrix(rnorm(N * T, 0, P$x_sd), N, T)
  x <- P$x_mean + P$x_alpha * alpha + P$feedback * e[, 1:T, drop = FALSE] + u
  eps <- e[, 2:(T + 1), drop = FALSE]
  x <- as.vector(t(x)); eps <- as.vector(t(eps))
  data.frame(id, t = ix$t, x, y = P$mu0 + P$beta * x + alpha[id] + eps, alpha_true = alpha[id], eps_true = eps)
}

gen_P5 <- function(P) {
  N <- P$N; T <- P$T; ix <- .panel_index(N, T); id <- ix$id
  alpha <- rnorm(N, 0, P$alpha_sd)
  mu <- P$mu_mean + P$mu_alpha * alpha + rnorm(N, 0, P$mu_sd)
  s <- sample(P$s_types, N, replace = TRUE, prob = P$s_probs)
  eta <- rnorm(N, 0, P$eta_sd)
  Es2 <- sum(P$s_probs * P$s_types^2)
  b <- switch(P$design,
              within_sd = P$avg_slope + P$kappa * (s^2 - Es2) + eta,
              level = P$avg_slope + P$kappa * (mu - P$mu_mean) + eta,
              stop("gen_P5: design must be within_sd or level"))
  x <- mu[id] + s[id] * rnorm(N * T)
  eps <- rnorm(N * T, 0, P$sigma_e)
  data.frame(id, t = ix$t, x, y = P$mu0 + b[id] * x + alpha[id] + eps, b_true = b[id], s_true = s[id],
             mu_true = mu[id], alpha_true = alpha[id], eps_true = eps)
}

#' P6: a balanced panel and two attrition rules (absorbing: a merchant that leaves never returns; month 1 is
#' always observed).
gen_P6 <- function(P) {
  N <- P$N; T <- P$T; ix <- .panel_index(N, T); id <- ix$id
  alpha <- rnorm(N, 0, P$alpha_sd)
  x <- P$x_mean + P$x_alpha * alpha[id] + rnorm(N * T)
  eps <- rnorm(N * T, 0, P$sigma_e)
  y <- P$mu0 + P$beta * x + alpha[id] + eps
  leave <- matrix(runif(N * T) < pnorm(P$haz0 + P$haz_alpha * alpha), N, T); leave[, 1] <- FALSE
  low <- matrix(y < P$floor, N, T, byrow = TRUE); low[, 1] <- FALSE
  stay <- function(m) as.vector(apply(m, 1, function(v) as.integer(cumsum(v) == 0)))
  data.frame(id, t = ix$t, x, y, alpha_true = alpha[id], eps_true = eps, s_alpha = stay(leave), s_eps = stay(low))
}
# gen_P2, gen_P4 and gen_P7 draw blocks P2, P4 and P7 as described in the text
set.seed(1002); d2 <- gen_P2(P2); set.seed(1003); d3 <- gen_P3(P3); set.seed(1004); d4 <- gen_P4(P4)
set.seed(1005); d5 <- gen_P5(P5); set.seed(1015); d5b <- gen_P5(P5b); set.seed(1006); d6 <- gen_P6(P6)
set.seed(1007); d7 <- gen_P7(P7)
```

Below are the functions shared by all the blocks of Section 6: unit means average over the observed periods only, the cluster-robust variance sums by unit (CR0), and the Wald test takes the quadratic form. `fd_fit`, `fd_ar1`, `chamberlain_md`, `mg_fit`, `slope_test`, and `re_fit` are helper functions written the same way, each described where it is used.

```r
.as_X <- function(X, nm = "x") {
  if (is.null(dim(X))) X <- matrix(as.numeric(X), ncol = 1L, dimnames = list(NULL, nm))
  X <- as.matrix(X)
  if (is.null(colnames(X))) colnames(X) <- paste0(nm, seq_len(ncol(X)))
  storage.mode(X) <- "double"
  X
}
.unit_index <- function(id) match(id, unique(id))                  # 1..G in order of appearance
.unit_means <- function(M, g) {                                       # G x k matrix of unit means (row j = unit j)
  M <- as.matrix(M)
  rowsum(M, g) / as.vector(rowsum(rep(1, nrow(M)), g))
}
# CR0 cluster-robust variance: bread * sum_g (Z_g' u_g)(Z_g' u_g)' * bread
.cl_vcov <- function(Z, u, g, bread = solve(crossprod(Z))) {
  S <- rowsum(Z * as.vector(u), g)
  V <- bread %*% crossprod(S) %*% bread
  dimnames(V) <- list(colnames(Z), colnames(Z))
  V
}
.wald <- function(b, V) {
  stat <- drop(crossprod(b, solve(V, b)))
  list(stat = stat, df = length(b), p = pchisq(stat, length(b), lower.tail = FALSE))
}
ci95 <- function(est, se) c(est - qnorm(0.975) * se, est + qnorm(0.975) * se)
within_tf <- function(y, X, id) {
  X <- .as_X(X); g <- .unit_index(id)
  m <- .unit_means(cbind(y, X), g)[g, , drop = FALSE]
  Xd <- X - m[, -1, drop = FALSE]; colnames(Xd) <- colnames(X)
  list(y = as.numeric(y) - m[, 1], X = Xd)
}
fe_fit <- function(y, X, id) {
  X <- .as_X(X); g <- .unit_index(id); w <- within_tf(y, X, id); Xd <- w$X
  A <- solve(crossprod(Xd)); b <- drop(A %*% crossprod(Xd, w$y)); names(b) <- colnames(X)
  e <- drop(w$y - Xd %*% b)
  n <- length(e); G <- max(g); k <- ncol(X)
  s2 <- sum(e^2) / (n - G - k)
  V_classic <- s2 * A; dimnames(V_classic) <- list(colnames(X), colnames(X))
  list(b = b, V = .cl_vcov(Xd, e, g, A), V_classic = V_classic, sigma2_eps = s2, resid = e, n = n, G = G)
}
mundlak_fit <- function(y, X, id, time_fe = FALSE, t = NULL) {
  X <- .as_X(X); g <- .unit_index(id)
  Xbar <- .unit_means(X, g)[g, , drop = FALSE]; colnames(Xbar) <- paste0("xbar_", colnames(X))
  Z <- cbind("(Intercept)" = 1, X, Xbar)
  if (isTRUE(time_fe)) {
    if (is.null(t)) stop("mundlak_fit: time_fe = TRUE needs the period index t")
    lev <- sort(unique(t))[-1]
    D <- outer(t, lev, "==") * 1; colnames(D) <- paste0("t_", lev)
    Dbar <- .unit_means(D, g)[g, , drop = FALSE]; colnames(Dbar) <- paste0("tbar_", lev)
    Z <- cbind(Z, D, Dbar)
  }
  q <- qr(Z, tol = 1e-9)
  keep <- sort(q$pivot[seq_len(q$rank)])
  core <- c("(Intercept)", colnames(X), colnames(Xbar))
  if (!all(match(core, colnames(Z)) %in% keep)) stop("mundlak_fit: X or Xbar is collinear")
  dropped <- colnames(Z)[-keep]
  Z <- Z[, keep, drop = FALSE]
  A <- solve(crossprod(Z)); coef <- drop(A %*% crossprod(Z, y)); names(coef) <- colnames(Z)
  e <- as.numeric(y) - drop(Z %*% coef)
  list(b_x = coef[colnames(X)], b_xbar = coef[colnames(Xbar)], V = .cl_vcov(Z, e, g, A), coef = coef,
       n = length(e), G = max(g), dropped = dropped)
}
hausman_classic <- function(b_fe, V_fe, b_re, V_re) {
  d <- as.numeric(b_fe) - as.numeric(b_re)
  .wald(d, as.matrix(V_fe) - as.matrix(V_re))
}
hausman_robust <- function(fit) {
  nm <- names(fit$b_xbar)
  .wald(fit$b_xbar, fit$V[nm, nm, drop = FALSE])
}
lead_test <- function(y, X, id, t) {
  X <- .as_X(X); o <- order(id, t); y <- as.numeric(y)[o]; X <- X[o, , drop = FALSE]; id <- id[o]; t <- t[o]
  n <- length(y); has_lead <- c(id[-1] == id[-n] & t[-1] == t[-n] + 1, FALSE)
  Xl <- rbind(X[-1, , drop = FALSE], NA); colnames(Xl) <- paste0(colnames(X), "_lead")
  k <- which(has_lead)
  f <- fe_fit(y[k], cbind(X[k, , drop = FALSE], Xl[k, , drop = FALSE]), id[k])
  dn <- colnames(Xl)
  w <- .wald(f$b[dn], f$V[dn, dn, drop = FALSE])
  structure(w$p, stat = w$stat, df = w$df, delta = f$b[dn], se = sqrt(diag(f$V)[dn]), n = length(k))
}
```

### 6.2 FE as Hero: Static

On the static block we compare three estimators.

```r
library(plm)
d <- pdata.frame(static, index = c("id", "t"))
pool <- plm(y ~ x, d, model = "pooling")
re   <- plm(y ~ x, d, model = "random")
fe   <- plm(y ~ x, d, model = "within")
phtest(fe, re)                                     # Hausman
```

`phtest` reports $\chi^2 = 3362$, with a $p$-value below the smallest number the software can display. Next we add between, FD, and the dummy-variable regression, and check the weighted average of Section 3.5.

```r
cl <- function(m) vcovHC(m, method = "arellano", type = "HC0", cluster = "group")   # clustered by merchant (CR0)
be <- plm(y ~ x, d, model = "between"); fd <- plm(y ~ x, d, model = "fd")
lsdv <- lm(y ~ x + factor(id), data = static)                                      # one dummy per merchant
show <- function(nm, m, clu = TRUE) cat(sprintf("%-8s %.4f  SE %.4f%s\n", nm, coef(m)[["x"]], sqrt(vcov(m)["x", "x"]),
                                         if (clu) sprintf("  cluster %.4f", sqrt(cl(m)["x", "x"])) else ""))
show("pooled", pool); show("RE", re); show("FE", fe); show("between", be, FALSE); show("FD", fd); show("LSDV", lsdv, FALSE)
ec <- ercomp(re)                                         # Swamy-Arora variance components and varphi
xw <- static$x - ave(static$x, static$id); xb <- tapply(static$x, static$id, mean)
S_W <- sum(xw^2); S_B <- 6 * sum((xb - mean(xb))^2)       # within and between variation of X
w <- (1 - ec$theta)^2 * S_B / (S_W + (1 - ec$theta)^2 * S_B)
cat(sprintf("varphi %.4f, sigma2_alpha %.4f, sigma2_eps %.4f; weight on between %.4f; w BE + (1 - w) FE - RE = %.1e\n",
            ec$theta, ec$sigma2[["id"]], ec$sigma2[["idios"]], w, w * coef(be)[["x"]] + (1 - w) * coef(fe)[["x"]] - coef(re)[["x"]]))
# pooled   0.9511  SE 0.0172  cluster 0.0220
# RE       0.7390  SE 0.0171  cluster 0.0182
# FE       0.5063  SE 0.0175  cluster 0.0179
# between  1.5498  SE 0.0327
# FD       0.5183  SE 0.0174  cluster 0.0212
# LSDV     0.5063  SE 0.0175
# varphi 0.3785, sigma2_alpha 0.2520, sigma2_eps 0.9517; weight on between 0.2230; w BE + (1 - w) FE - RE = 2.1e-15
```

In the static Nimbus example the default FE standard error is 0.0175 and the merchant-clustered one is 0.0179, not much different (because the errors in this setup are nearly independent), but in real panels serial correlation is often sizeable, and clustering is the default rather than an option. The two standard errors of pooled OLS differ more (0.0172 against 0.0220), because its composite error contains $\alpha_i$.

```r
mk <- mundlak_fit(static$y, static$x, static$id)
mk_tw <- mundlak_fit(static$y, static$x, static$id, time_fe = TRUE, t = static$t)
tw <- plm(y ~ x, d, model = "within", effect = "twoways")              # two-way FE
hc <- hausman_classic(coef(fe), vcov(fe), coef(re)["x"], vcov(re)["x", "x"]); hr <- hausman_robust(mk)
lt <- lead_test(static$y, static$x, static$id, static$t)
cat(sprintf("Mundlak: beta %.6f (FE %.6f), xi %.4f (BE - FE %.4f), cluster SE %.4f and %.4f\n", mk$b_x, coef(fe),
            mk$b_xbar, coef(be)[["x"]] - coef(fe), sqrt(mk$V["x", "x"]), sqrt(mk$V["xbar_x", "xbar_x"])))
cat(sprintf("two-way Mundlak %.6f, two-way FE %.6f\n", mk_tw$b_x, coef(tw)))
cat(sprintf("Hausman: classical %.2f (p %.2g), robust %.2f (p %.2g)\n", hc$stat, hc$p, hr$stat, hr$p))
cat(sprintf("lead test: delta %.4f (SE %.4f), Wald %.2f, p %.4f\n", attr(lt, "delta"), attr(lt, "se"), attr(lt, "stat"), lt))
ar <- fd_ar1(fd_fit(static$y, static$x, static$id, static$t))         # AR(1) coefficient of the FD residuals, cluster SE
ch <- chamberlain_md(static$y, static$x, static$id, static$t)           # Chamberlain by minimum distance (chapter 7)
cat(sprintf("FD residuals: rho %.4f (SE %.4f); Chamberlain beta %.4f (SE %.4f), overid %.2f on %d df, p %.3f\n",
            ar$rho, ar$se, ch$beta, ch$se, ch$stat, ch$df, ch$p))
# Mundlak: beta 0.506302 (FE 0.506302), xi 1.0435 (BE - FE 1.0435), cluster SE 0.0179 and 0.0369
# two-way Mundlak 0.504978, two-way FE 0.504978
# Hausman: classical 3362.02 (p 0), robust 798.95 (p 9.1e-176)
# lead test: delta -0.0342 (SE 0.0205), Wald 2.78, p 0.0952
# FD residuals: rho -0.4906 (SE 0.0157); Chamberlain beta 0.5056 (SE 0.0174), overid 30.65 on 34 df, p 0.633
```

The result lines match the theorems of Section 3 one by one: the promotion coefficient of the Mundlak regression agrees with FE to the sixth decimal place, $\hat\xi$ equals the between estimate minus FE, and the two-way version equals TWFE (Theorem 5); the classical and robust Hausman tests both reject RE (Theorem 7); the autocorrelation of the FD residuals is close to $-0.5$, Chamberlain's overidentification test does not reject, and the $p$-value of the lead test falls between 5% and 10% (Sections 3.4 and 3.5).

![Point estimates and 95% intervals of six estimators on block P1 (the inherited static panel, 600 merchants, 6 months): pooled OLS, between, RE, FE, FD, and Chamberlain minimum distance. Intervals are clustered by merchant, except that the between estimator, with one row per merchant, uses classical standard errors; hollow diamonds are each estimator's probability limit as the number of merchants goes to infinity, and the dashed line is the true value 0.5.](assets/fig/fig_10_estimators.svg)

The six estimators line up clearly: between, which uses only variation across merchants, is the furthest off, pooled OLS comes next, and RE stops partway from between toward FE; the three estimators that rely only on strict exogeneity all land near the truth, and each interval covers its own probability limit.

### 6.3 Blocks P2 and P3: FE Versus FD, and Feedback

```r
for (k in c("iid", "rw")) {                                  # P2, one sample
  y <- d2[[paste0("y_", k)]]
  ff <- fe_fit(y, d2$x, d2$id); fdf <- fd_fit(y, d2$x, d2$id, d2$t)
  cat(sprintf("P2 %-3s: FE %.4f (cluster SE %.4f), FD %.4f (%.4f), FD residual rho %.3f\n", k, ff$b,
              sqrt(ff$V[1, 1]), fdf$b[["x"]], sqrt(fdf$V["x", "x"]), fd_ar1(fdf)$rho))
}
set.seed(1012)
mc <- t(vapply(seq_len(P2$R_mc), function(r) {               # P2, Monte Carlo
  d <- gen_P2(P2)
  out <- numeric(0)
  for (k in c("iid", "rw")) {
    y <- d[[paste0("y_", k)]]
    ff <- fe_fit(y, d$x, d$id); fdf <- fd_fit(y, d$x, d$id, d$t)
    out <- c(out, ff$b[1], sqrt(ff$V[1, 1]), fdf$b["x"], sqrt(fdf$V["x", "x"]))
  }
  out
}, numeric(8)))
colnames(mc) <- c("fe_iid", "se_fe_iid", "fd_iid", "se_fd_iid", "fe_rw", "se_fe_rw", "fd_rw", "se_fd_rw")
cat(sprintf("P2 Monte Carlo (%d samples): SD iid FE %.4f FD %.4f; SD random walk FE %.4f FD %.4f\n", nrow(mc),
            sd(mc[, "fe_iid"]), sd(mc[, "fd_iid"]), sd(mc[, "fe_rw"]), sd(mc[, "fd_rw"])))
f3 <- fe_fit(d3$y, d3$x, d3$id); fd3 <- fd_fit(d3$y, d3$x, d3$id, d3$t)   # P3, feedback
lt3 <- lead_test(d3$y, d3$x, d3$id, d3$t); ch3 <- chamberlain_md(d3$y, d3$x, d3$id, d3$t)
cat(sprintf("P3: FE %.4f [%.4f, %.4f], FD %.4f [%.4f, %.4f]\n", f3$b, ci95(f3$b, sqrt(f3$V[1, 1]))[1],
            ci95(f3$b, sqrt(f3$V[1, 1]))[2], fd3$b[["x"]], ci95(fd3$b[["x"]], sqrt(fd3$V["x", "x"]))[1],
            ci95(fd3$b[["x"]], sqrt(fd3$V["x", "x"]))[2]))
cat(sprintf("P3: lead test delta %.4f, p %.2g; Chamberlain overid %.1f on %d df\n", attr(lt3, "delta"), lt3,
            ch3$stat, ch3$df))
# P2 iid: FE 0.4883 (cluster SE 0.0220), FD 0.4906 (0.0319), FD residual rho -0.537
# P2 rw : FE 0.4860 (cluster SE 0.0301), FD 0.4796 (0.0211), FD residual rho -0.026
# P2 Monte Carlo (1000 samples): SD iid FE 0.0239 FD 0.0310; SD random walk FE 0.0318 FD 0.0195
# P3: FE 0.5743 [0.5511, 0.5975], FD 0.7235 [0.6977, 0.7494]
# P3: lead test delta -0.4262, p 2.5e-293; Chamberlain overid 1492.8 on 34 df
```

In the P2 Monte Carlo the averages of all four estimates are near 0.501, and the coverage rates of the cluster intervals lie between 0.946 and 0.953; they differ only in precision: with iid errors the standard deviation of FD is 1.30 times that of FE, and with a random walk that of FE is 1.63 times that of FD. Each of the two P3 estimates aims at its own probability limit from Section 3.4, 0.5735 and 0.7206 respectively.

### 6.4 Block P4: Hausman-Taylor

```r
pd4 <- pdata.frame(d4, index = c("id", "t"))
ht <- plm(y ~ x1 + x2 + w | x1 | x2, data = pd4, model = "random", random.method = "ht", inst.method = "baltagi")
re4 <- plm(y ~ x1 + x2 + w, data = pd4, model = "random")
po4 <- plm(y ~ x1 + x2 + w, data = pd4, model = "pooling")
fe4 <- plm(y ~ x1 + x2, data = pd4, model = "within")
row4 <- function(nm, m, v) {                                    # classical (non-clustered) SEs from vcov()
  e <- coef(m)[[v]]; s <- sqrt(vcov(m)[v, v])
  cat(sprintf("%-6s %-2s %.4f  SE %.4f  [%.3f, %.3f]\n", nm, v, e, s, e - qnorm(0.975) * s, e + qnorm(0.975) * s))
}
row4("HT", ht, "w"); row4("RE", re4, "w"); row4("HT", ht, "x2"); row4("RE", re4, "x2")
cat(sprintf("pooled OLS gamma %.4f; FE beta1 %.4f, beta2 %.4f; varphi HT %.3f, RE %.3f\n", coef(po4)[["w"]],
            coef(fe4)[["x1"]], coef(fe4)[["x2"]], ercomp(ht)$theta, ercomp(re4)$theta))
u <- d4[d4$t == 1, c("id", "w")]; u$x1bar <- tapply(d4$x1, d4$id, mean)[as.character(u$id)]
fs <- lm(w ~ x1bar, data = u)                                    # unit-level first stage of W_i on Xbar_1i
cat(sprintf("first stage: coefficient %.4f, robust F %.1f\n", coef(fs)[["x1bar"]],
            coef(fs)[["x1bar"]]^2 / sandwich::vcovHC(fs, type = "HC0")["x1bar", "x1bar"]))
# HT     w  1.4493  SE 0.1973  [1.063, 1.836]
# RE     w  2.9130  SE 0.0676  [2.780, 3.045]
# HT     x2 0.4972  SE 0.0150  [0.468, 0.527]
# RE     x2 0.7781  SE 0.0155  [0.748, 0.809]
# pooled OLS gamma 2.8257; FE beta1 0.3291, beta2 0.4972; varphi HT 0.710, RE 0.504
# first stage: coefficient 0.2656, robust F 1566.1
```

The standard errors here are the classical ones from `vcov()`, not cluster standard errors. Certified merchants make up 42.9% of the sample; the correlation of certification with $\bar X_{1i}$ is 0.554, and its correlation with talent is 0.330.

### 6.5 Blocks P5 and P5b: Random Slopes

```r
for (blk in c("P5", "P5b")) {
  dd <- if (blk == "P5") d5 else d5b
  f <- fe_fit(dd$y, dd$x, dd$id)                           # FE with cluster SE
  m <- mg_fit(dd$y, dd$x, dd$id)                           # mean group: average of merchant-by-merchant OLS slopes
  s2hat <- ave((dd$x - ave(dd$x, dd$id))^2, dd$id) * 8 / 7  # within variance of X_i, one value per merchant
  t_xbar <- slope_test(dd$y, dd$x, dd$id, ave(dd$x, dd$id))  # FE with X_it (Z_i - mean Z) added, cluster Wald
  t_s2 <- slope_test(dd$y, dd$x, dd$id, s2hat)
  cat(sprintf("%-3s: FE %.4f [%.4f, %.4f]; MG %.4f [%.4f, %.4f]; slope tests p (Xbar) %.3g (within var.) %.3g\n", blk,
              f$b, ci95(f$b, sqrt(f$V[1, 1]))[1], ci95(f$b, sqrt(f$V[1, 1]))[2], m$est, ci95(m$est, m$se)[1],
              ci95(m$est, m$se)[2], t_xbar$p, t_s2$p))
}
# P5 : FE 0.6834 [0.6582, 0.7086]; MG 0.5450 [0.5075, 0.5825]; slope tests p (Xbar) 0.835 (within var.) 3.71e-07
# P5b: FE 0.5018 [0.4698, 0.5339]; MG 0.5061 [0.4659, 0.5462]; slope tests p (Xbar) 9.51e-53 (within var.) 0.88
```

In the sample, the FE weights of P5's four types of merchant are 3.6%, 12.8%, 27.7%, and 55.9%. Weighted by each merchant's own within-unit sum of squares, the average of the true slopes is 0.6674, and FE's 0.683 is a little more than one standard error away from it; the equally weighted average is 0.5021.

### 6.6 Block P6: An Unbalanced Panel

```r
fe6 <- function(dd) {
  f <- fe_fit(dd$y, dd$x, dd$id); e <- f$b[[1]]; s <- sqrt(f$V[1, 1])
  sprintf("%.4f [%.4f, %.4f]", e, e - qnorm(0.975) * s, e + qnorm(0.975) * s)
}
sa <- d6[d6$s_alpha == 1, ]; se6 <- d6[d6$s_eps == 1, ]       # attrition on alpha_i; attrition on the shock
cat(sprintf("P6 FE: full %s; attrition on alpha %s (kept %.3f); on the shock %s (kept %.3f)\n", fe6(d6), fe6(sa),
            nrow(sa) / nrow(d6), fe6(se6), nrow(se6) / nrow(d6)))
pa <- pdata.frame(sa, index = c("id", "t"))
mk6 <- mundlak_fit(sa$y, sa$x, sa$id); mk6_tw <- mundlak_fit(sa$y, sa$x, sa$id, time_fe = TRUE, t = sa$t)
naive <- coef(lm(y ~ x + ave(x, id) + factor(t), data = sa))[["x"]]   # two-way without averaged month dummies
re6 <- plm(y ~ x, data = pa, model = "random")
cat(sprintf("unbalanced: Mundlak %.6f vs FE %.6f; two-way Mundlak %.6f vs TWFE %.6f; naive %.6f\n", mk6$b_x,
            coef(plm(y ~ x, data = pa, model = "within")), mk6_tw$b_x,
            coef(plm(y ~ x, data = pa, model = "within", effect = "twoways")), naive))
cat(sprintf("RE %.4f with varphi_i from %.3f to %.3f; robust Hausman %.1f\n", coef(re6)[["x"]],
            min(ercomp(re6)$theta), max(ercomp(re6)$theta), hausman_robust(mk6)$stat))
# P6 FE: full 0.4974 [0.4805, 0.5143]; attrition on alpha 0.4916 [0.4713, 0.5119] (kept 0.754); on the shock 0.4437 [0.4241, 0.4632] (kept 0.723)
# unbalanced: Mundlak 0.491571 vs FE 0.491571; two-way Mundlak 0.491178 vs TWFE 0.491178; naive 0.492192
# RE 0.8669 with varphi_i from 0.095 to 0.398; robust Hausman 5149.8
```

Under attrition on the shock, the FE bias is $-0.056$: the estimate falls 5.64 standard errors below the truth, which is not sampling error. After attrition on talent, merchants are observed for 6.03 months on average, and 9.7% of merchants are left with only 1 month; these merchants contribute nothing to FE. The two-way Mundlak regression without the month means differs from TWFE by 0.00101, a small amount, but it shows which set of regressors the identity depends on.

### 6.7 Block P7: The Size of the Hausman Test

```r
both_tests <- function(d) {
  fe <- fe_fit(d$y, d$x, d$id); re <- re_fit(d$y, d$x, d$id)       # re_fit: Swamy-Arora RE, classical vcov
  hc <- hausman_classic(fe$b, fe$V_classic, re$b["x"], re$V_classic["x", "x"])
  hr <- hausman_robust(mundlak_fit(d$y, d$x, d$id))
  c(fe = unname(fe$b), re = unname(re$b["x"]), H_classic = hc$stat, p_classic = hc$p, H_robust = hr$stat,
    p_robust = hr$p, vd_classic = unname(fe$V_classic[1, 1] - re$V_classic["x", "x"]), varphi = re$varphi)
}
s1 <- both_tests(d7)
cat(sprintf("P7 sample: H classical %.3f (p %.3f), robust %.3f (p %.3f)\n", s1[["H_classic"]], s1[["p_classic"]],
            s1[["H_robust"]], s1[["p_robust"]]))
set.seed(1017)
mc7 <- t(vapply(seq_len(P7$R_mc), function(r) both_tests(gen_P7(P7)), numeric(8)))
rej <- function(p, a) mean(p < a)
cat(sprintf("P7 Monte Carlo (%d samples), rejection rate at 1%%, 5%%, 10%%: classical %.4f %.4f %.4f; robust %.4f %.4f %.4f\n",
            nrow(mc7), rej(mc7[, "p_classic"], 0.01), rej(mc7[, "p_classic"], 0.05), rej(mc7[, "p_classic"], 0.10),
            rej(mc7[, "p_robust"], 0.01), rej(mc7[, "p_robust"], 0.05), rej(mc7[, "p_robust"], 0.10)))
cat(sprintf("SD of FE - RE across samples %.4f; classical formula, averaged %.4f\n", sd(mc7[, "fe"] - mc7[, "re"]),
            mean(sqrt(mc7[, "vd_classic"]))))
# P7 sample: H classical 1.572 (p 0.210), robust 0.597 (p 0.440)
# P7 Monte Carlo (2000 samples), rejection rate at 1%, 5%, 10%: classical 0.0795 0.1690 0.2530; robust 0.0075 0.0510 0.0940
# SD of FE - RE across samples 0.0105; classical formula, averaged 0.0073
```

In this sample neither test rejects; a single sample does not reveal the problem, and size can only be read from the Monte Carlo. At all three levels, 1%, 5%, and 10%, the rejection rate of the classical test is far above the nominal level, while that of the robust test stays within the Monte Carlo error band; the last result line shows the mechanism described in Section 3.6.

### 6.8 The Grand Reconciliation

| Block | Estimator | Estimate | 95% interval | Target and true value | Covered |
|---|---|---|---|---|---|
| P1 | Pooled OLS | 0.951 | [0.908, 0.994] | $\beta$ 0.5 (probability limit 0.9698) | No |
| P1 | RE | 0.739 | [0.703, 0.775] | $\beta$ 0.5 (probability limit 0.7537) | No |
| P1 | FE | 0.506 | [0.471, 0.541] | $\beta$ 0.5 | Yes |
| P1 | FD | 0.518 | [0.477, 0.560] | $\beta$ 0.5 | Yes |
| P1 | Chamberlain MD | 0.5056 | [0.4716, 0.5396] | $\beta$ 0.5 | Yes |
| P1 | Hausman | Classical 3362.02, robust 798.95 | | Reject RE | |
| P3 | FE | 0.574 | [0.551, 0.597] | $\beta$ 0.5 (probability limit 0.5735) | No |
| P3 | FD | 0.724 | [0.698, 0.749] | $\beta$ 0.5 (probability limit 0.7206) | No |
| P4 | Hausman-Taylor $\gamma$, classical standard error | 1.449 | [1.063, 1.836] | $\gamma$ 1.5 | Yes |
| P4 | RE $\gamma$, classical standard error | 2.913 | [2.780, 3.045] | $\gamma$ 1.5 | No |
| P5 | FE | 0.683 | [0.658, 0.709] | Average slope 0.5 ($\beta_*$ 0.66125) | No |
| P5 | Mean group | 0.545 | [0.5075, 0.5825] | Average slope 0.5 | No |
| P5b | FE | 0.502 | [0.470, 0.534] | Average slope 0.5 | Yes |
| P6 | FE, attrition on talent | 0.4916 | [0.4713, 0.5119] | $\beta$ 0.5 | Yes |
| P6 | FE, attrition on the shock | 0.4437 | [0.4241, 0.4632] | $\beta$ 0.5 | No |
| P7 | 5% size of the classical and robust Hausman tests | 0.1690 and 0.0510 | | Nominal 0.05 | |

Almost every row of the table that fails to cover the truth does so by design: the pooled OLS and RE of P1, the FE and FD of P3, the FE of P5, and the FE of P6 under attrition on the shock each cover their own probability limit or pseudo-true value; P4's RE counts talent as part of certification and is biased upward. These biases come from the identifying assumptions, not from luck. The only exception is P5's mean group, which should aim at the average slope and this time does not cover it, the kind of miss that a nominal 95% interval makes about once in twenty.

The reconciliation of this section can be summed up as follows: in the static panel FE is the estimator that recovers the true promotion effect, and Hausman decisively rejects RE; but once strict exogeneity fails, slopes are correlated with the within-unit fluctuations of promotion (P5; P5b, whose slopes vary only with the level of promotion, is unaffected), or attrition is related to the shock, FE aims at a different quantity (a problem with any one of the three is enough), and in a simulation this quantity can be computed in advance.

## 7 Failure Modes and Robustness

The problems in the simulation are manufactured, and in real panels they hide deeper. Below we lay out the most common failure modes and the responses to them.

The first choice to confront head-on is FE or RE. When RE.1 through RE.3 hold, RE is more efficient and can also estimate time-invariant covariates, but it wagers on the strong assumption that "the individual effect is uncorrelated with the covariates," which in observational data usually does not hold (units that select into a state often differ in unobservable dimensions to begin with). The Hausman test gives a data-driven criterion, but it is not a cure-all: a failure to reject does not mean RE is correct; the test may simply lack power, and Wooldridge (2010) p. 334 warns that accepting RE in that case may be a Type II error; the difference of variances that the classical version uses can be non-positive-definite in finite samples, making the statistic negative and uninterpretable, and when RE.3 fails its size is wrong as well: on P7 the nominal 5% test rejects 16.9% of the time. The safe practice is to use the robust version based on the Mundlak regression and to report $\hat\xi$ as a coefficient with economic meaning; use FE by default, putting robustness ahead of efficiency, and turn to Hausman-Taylor or CRE only when a time-invariant variable really must be estimated, defending the exclusion restrictions one by one.

Strict exogeneity is the most easily overlooked weak spot of FE. It is violated not only by a lagged dependent variable but by any covariate that "responds to past shocks": promotion adjusted to last period's performance, price adjusted to last period's demand, investment adjusted to last period's profit; these feedbacks all make the covariate correlate with past errors and bias FE. The key question is: could the value of this covariate be influenced by past outcomes? If so, treat it as a predetermined variable and use an instrument under sequential exogeneity (a lagged value) rather than putting it straight into an FE regression. There are three practical diagnostics: the lead test, Chamberlain's overidentification test, and whether FE and FD differ by more than sampling error; on P3 the gap between 0.574 and 0.724 is such an alarm. The dynamic case, with sequential exogeneity and a lagged dependent variable, is taken up in Chapter 11.

Slope heterogeneity is the second easily overlooked problem. As long as the slope is unrelated to the within-unit fluctuations of the regressor (RS.2), FE still aims at the average slope, as in P5b; once the two are correlated, FE does not "go wrong"; it just quietly switches estimand. What FE reports is the average slope weighted by within-unit variance, dominated by the merchants whose promotion fluctuates most. If the policy question concerns the average response of all merchants, check where the weights fall: add to the FE regression interactions of $X_{it}$ with unit-level variables (the unit mean, the within-unit variance, size), of which the within-unit variance interaction is the one that matters for FE; or report an estimator such as mean group that aims directly at the average slope. Imai and Kim (2019) reach the same conclusion for a binary treatment: the unit FE estimator converges to an average of unit effects weighted by the within-unit variance of the treatment.

Attrition harms FE only when it is correlated with the idiosyncratic error. If merchants leave one after another because of low talent, FE is unaffected; if they leave because of one bad month, the remaining sample is systematically tilted toward good shocks, and on P6 the bias is about a tenth of the true value and does not vanish however large the sample. The selection-indicator test of Section 3.9 can help here; when a refreshment sample is available, Hirano, Imbens, Ridder and Rubin (2001) give a model that can identify more general attrition mechanisms.

FE can remove only time-invariant heterogeneity. If a merchant's "talent" itself changes, with some shops growing and others declining, and promotion follows these trends, demeaning is not enough. The random trend model of Wooldridge (2010) gives each unit its own linear trend; it can be estimated by differencing first and then applying FE (§11.7.1), or by detrending unit by unit (§11.7.2), at the cost of needing at least three periods, even less within-unit variation, and less precise estimates.

Measurement error is another weak spot of FE and FD. Demeaning and differencing remove talent, but they also remove most of the persistent variation in the regressor, so measurement noise makes up a larger share of the within-unit variation that remains, and the attenuation caused by classical measurement error can therefore be more severe than in pooled OLS; the more persistent the regressor and the less persistent the noise, the bigger the problem (Wooldridge 2010 §11.5, equation (11.44)). Promotion spending is booked by the platform's own systems, so this problem is minor here; if the regressor comes from a survey or a sample estimate, it must be taken seriously.

Finally, inference. In a panel the errors of the same unit are almost always correlated, so standard errors should be clustered by unit by default, with a large enough number of units; if treatment in a TWFE regression starts at different times, the meaning of the coefficient has to be reexamined with the methods of Chapter 19.

Put these failure modes together, and the robustness of panel methods comes down to three judgments: whether your unobservable heterogeneity really is time-invariant, whether your covariates really are strictly exogenous, and whether the average you want is the one that FE automatically delivers. The first two decide whether FE can be used, and the third decides how its number should be read.

## 8 Further Reading

::: {.readings}
Required reading, in suggested reading order:

- Wooldridge (2010, *Econometric Analysis of Cross Section and Panel Data*, 2nd ed.), Chapters 10 to 11. The clearest textbook treatment of FE, RE, and strict versus sequential exogeneity, and the source of this chapter's assumption labels; unbalanced panels are in §19.9.
- Mundlak (1978, *Econometrica*). Writing the individual effect as a function of the within means of the covariates builds a bridge between FE and RE (CRE), and gives another angle on what the Hausman test is testing.
- Hausman (1978, *Econometrica*). The original source of specification testing and of this chapter's Hausman test, of a piece with the QMLE tests of Chapter 6.
- Wooldridge (2019, *Journal of Econometrics*). The Mundlak and Chamberlain approaches on unbalanced panels, with unit means averaged over the observed periods only.

Further reading:

- Chamberlain (1982, *Journal of Econometrics*). Panel models in a multivariate regression framework, strict exogeneity stated as a conditional restriction, and minimum distance estimation.
- Hausman and Taylor (1981, *Econometrica*). Finding instruments inside the model to estimate the coefficients on time-invariant regressors.
- Arellano (1993, *Journal of Econometrics*). The Hausman test written as a robust Wald test in an augmented regression.
- Wooldridge (2005a, *Review of Economics and Statistics*). Conditions under which FE-type estimators remain consistent for the average slope when slopes are heterogeneous.
- Imai and Kim (2019, *American Journal of Political Science*) and Imai and Kim (2021, *Political Analysis*). Causal interpretations of unit FE and TWFE: the weights, the assumption of no dynamics, and how TWFE differs from multi-period DiD.
- Hirano, Imbens, Ridder and Rubin (2001, *Econometrica*). Identification under attrition when a refreshment sample is available.
- Baltagi (2021, *Econometric Analysis of Panel Data*, 6th ed.). The authoritative monograph on panel econometrics, covering static and dynamic, linear and nonlinear models in full; a desk reference.
:::

::: {.apa-refs}
- Arellano, M. (1993). On the testing of correlated effects with panel data. *Journal of Econometrics, 59*(1-2), 87-97. https://doi.org/10.1016/0304-4076(93)90040-C
- Baltagi, B. H. (2021). *Econometric analysis of panel data* (6th ed.). Springer. https://doi.org/10.1007/978-3-030-53953-5
- Chamberlain, G. (1982). Multivariate regression models for panel data. *Journal of Econometrics, 18*(1), 5-46. https://doi.org/10.1016/0304-4076(82)90094-X
- Hausman, J. A. (1978). Specification tests in econometrics. *Econometrica, 46*(6), 1251-1271. https://doi.org/10.2307/1913827
- Hausman, J. A., & Taylor, W. E. (1981). Panel data and unobservable individual effects. *Econometrica, 49*(6), 1377-1398. https://doi.org/10.2307/1911406
- Hirano, K., Imbens, G. W., Ridder, G., & Rubin, D. B. (2001). Combining panel data sets with attrition and refreshment samples. *Econometrica, 69*(6), 1645-1659. https://doi.org/10.1111/1468-0262.00260
- Imai, K., & Kim, I. S. (2019). When should we use unit fixed effects regression models for causal inference with longitudinal data? *American Journal of Political Science, 63*(2), 467-490. https://doi.org/10.1111/ajps.12417
- Imai, K., & Kim, I. S. (2021). On the use of two-way fixed effects regression models for causal inference with panel data. *Political Analysis, 29*(3), 405-415. https://doi.org/10.1017/pan.2020.33
- Mundlak, Y. (1978). On the pooling of time series and cross section data. *Econometrica, 46*(1), 69-85. https://doi.org/10.2307/1913646
- Oh, J.-H., Susarla, A., & Tan, Y. (2008). *Examining the diffusion of user-generated content in online social networks* [Working paper]. University of Washington. http://hdl.handle.net/2451/28035
- Susarla, A., Oh, J.-H., & Tan, Y. (2012). Social networks and the diffusion of user-generated content: Evidence from YouTube. *Information Systems Research, 23*(1), 23-41. https://doi.org/10.1287/isre.1100.0339
- Tambe, P., & Hitt, L. M. (2012). The productivity of information technology investments: New evidence from IT labor data. *Information Systems Research, 23*(3, Part 1 of 2), 599-617. https://doi.org/10.1287/isre.1110.0398
- Wooldridge, J. M. (2005a). Fixed-effects and related estimators for correlated random-coefficient and treatment-effect panel data models. *The Review of Economics and Statistics, 87*(2), 385-390. https://doi.org/10.1162/0034653053970320
- Wooldridge, J. M. (2010). *Econometric analysis of cross section and panel data* (2nd ed.). MIT Press.
- Wooldridge, J. M. (2019). Correlated random effects models with unbalanced panels. *Journal of Econometrics, 211*(1), 137-150. https://doi.org/10.1016/j.jeconom.2018.12.010
:::
