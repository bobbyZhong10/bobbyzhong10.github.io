---
title: "Synthetic Control and Panel Factor Models"
subtitle: "Building the Counterfactual from Donors"
seriesline: "Foundations of Information Systems Economics · Chapter 22"
series: "Foundations of Information Systems Economics: Empirical Methods for Information Systems Economics"
sidebartitle: "Foundations of Information Systems Economics"
sidebarsub: "Theory · Methods · Applications"
chapterlabel: "Chapter 22 · Synthetic Control and Panel Factor Models"
---

## Introduction

The Aster marketplace is piloting a fee reform in a single city. Here there is no threshold, and no natural "identical twin" city to serve as a control. Synthetic control takes a different route: it lets several untreated cities, combined with weights, jointly play the role of the counterfactual city that does not exist. Whether the pre-treatment path can be pieced together closely enough decides whether the gap that opens after treatment is credible.

The regression discontinuity design (RDD) of Chapter 21 relies on continuity near a threshold and speaks only to a local effect at the cutoff. Synthetic control (hereafter SC) likewise draws its identifying power from concrete institutional detail: it relies on the pre-treatment fit and a donor pool, and answers only the counterfactual for the unit that was intervened on. A single-city pilot does not come with a synthetic control automatically. This chapter follows Aster's single-city fee reform and asks when a weighted set of cities can truly stand in for "Aster without the reform." It then develops the same weighting idea in the order of Section 3: its unification with DiD and regression; synthetic DiD, which weights units and periods at once; a treated city that lies outside the convex hull of the donor pool; completing the counterfactual with a panel factor model; staggered adoption across several cities; inference for SC; and choosing the pilot cities before the pilot begins.

## 1 Two Numbers That Do Not Add Up

::: {.case}
Aster is a fictional two-sided online marketplace. We use simulated data to follow one design line. The payoff is that the data-generating process is fully known, so every estimator can be checked against the truth.

This design line is the single-city fee reform. Aster rolls out a seller fee reform in one metro only, with the last month before the reform denoted $T_0$, while twenty other metros carry on as usual. Management wants to know the effect of this city-level policy on the city's log GMV. The most natural calculation is a two-period, two-group DiD: the change in the treated city from before to after the reform, minus the average change in the twenty donor cities over the same window. This number is 0.244.
:::

This number is not credible, and earlier chapters have already shown where the problem lies: the fee reform's 0.244 goes wrong because parallel trends, the assumption discussed in Chapter 19, does not hold. In this DGP the cities' growth trends differ from the start: the treated city's pre-reform monthly trend slope in log GMV is 0.131, while the donor cities average only 0.117. The simple average of twenty cities is not the path the treated city would have followed without the reform, and using it as the counterfactual counts the pre-existing trend difference as part of the policy effect. The true post-reform average effect is only 0.095, so the 0.244 is about 2.6 times the truth.

The same platform has a second number that does not add up, and it shows that the problem is not confined to single-city designs. Aster later launches a seller support program, picking 15 of its 80 metros and rolling the program out in three waves of 5, with a true effect of 0.06 everywhere. The static TWFE of Section 4.3 of Chapter 19 gives 0.099; switching to the Callaway-Sant'Anna estimator recommended in the same chapter (hereafter CS) makes the result higher still: 0.139. Changing the estimator cannot rescue it, because the problem lies in who gets selected: the platform picks the cities hit hardest by the demand cycle, these cities respond more strongly to the same cycle to begin with, and in the recovery they naturally grow faster than the cities that were never selected.

SC offers a way around the trap: instead of using the simple average of twenty cities as the counterfactual, pick a convex combination of donors that reproduces the treated city's log-GMV trajectory before the reform as closely as possible. The closer the pre-reform fit, the more the gap that opens after the reform between the treated city's actual path and this synthetic version looks like the net effect of the policy. The second number calls for a generalization of this idea: build a synthetic control for each selected city, or estimate directly how each city responds to the common shocks.

This idea is the thread running through the chapter. Why it holds, and where it fails, are the questions the next seven sections answer.

## 2 The Economic Model and the Estimand

Before estimating anything, we need to be clear about two things: which causal quantity SC is meant to estimate, and what notation writes it down. The framework is the potential outcomes of Chapter 15.

There is a single treated unit here (denoted $i = 1$, the reform city Aster-Metro), treated for $t > T_0$; the donor cities $j = 2, \dots, J+1$ are never treated. The treatment indicator is $D_{it} = \mathbf{1}\{t \geq G_i\}$, where $G_1 = T_0 + 1$ and the donors have $G_j = \infty$ (following the notation of Section 2 of Chapter 19); $T_0$ denotes the last untreated period, month 24 in this example, and $\mathcal{T}$ denotes the total number of periods. In potential outcomes,

$$Y_{it} = Y_{it}(0) + \tau_{it} D_{it}.$$

The SC estimand is the ATT path of the treated unit, an effect curve that varies over time:

$$\tau_{1t} = Y_{1t}(1) - Y_{1t}(0) = \underbrace{Y_{1t}}_{\text{observed}} - \underbrace{Y_{1t}(0)}_{\text{counterfactual}}, \qquad t > T_0.$$

$Y_{1t}$ is observed, and all of the difficulty falls on $Y_{1t}(0)$, the log GMV the reform city would have realized without the reform. The whole machinery of SC does one thing: it constructs this unobserved counterfactual path. Averaging the effect path over the post-period gives $ATT_{\mathrm{post}} = (\mathcal{T} - T_0)^{-1}\sum_{t > T_0}\tau_{1t}$, and this is what the true value of 0.095 in Section 1 refers to.

Aster's data come in four blocks, which this chapter calls S1 to S4, each with its own estimand. S1 is the single-city fee reform above, with 36 months and 20 donors; the estimand is $ATT_{\mathrm{post}}$. S2 moves the same reform to Aster-Capital, the metro with the highest intrinsic level: its intrinsic level of 6.20 exceeds that of every donor (the highest donor is at 5.97). The estimand is again $ATT_{\mathrm{post}}$, with a true value of 0.079. The log GMV in S1 and S2 is an index without currency units, and its level cannot be compared with that of the two blocks below. S3 is the support program of Section 1: 80 metros and 48 months, with the three waves launching in months 29, 33, and 37; the outcome is the log of monthly GMV in millions of USD. Because there are several treatment dates, we switch to the potential outcomes $Y_{it}(g)$ of Chapter 19, indexed by the first treated period; the estimands are the $ATT(g,t) = \mathbb{E}[Y_{t}(g) - Y_{t}(0) \mid G = g]$ of Chapter 19 and their average over all treated city-months. In the DGP every $ATT(g,t)$ equals 0.06. S4 changes the question: the platform wants to pick 5 of 40 metros for a 6-month commission-cut pilot, with only the history of the previous 24 months in hand; the outcome is again the log of GMV in millions of USD. The estimand is the average effect over all cities, $\tau^{\mathrm{pop}}_t = N^{-1}\sum_{i=1}^{N}\tau_{it}$ ($N = 40$), with a true value of 0.04. It is not the pilot cities' own effect but "what would happen if every city cut commissions," and this distinction decides whom the weights in Section 3.8 should target.

To summarize: the target of SC is the ATT path of a single treated unit, $\tau_{1t}$, a curve that opens over time. It is local, local to the unit that was intervened on; S3 and S4 change the target to an average over several units. The next section covers the assumptions under which SC identifies these quantities.

## 3 Identification

Identification asks whether, under a set of assumptions about the counterfactual, the desired causal quantity can be written as a function of the distribution of the observable data; the answer depends only on the design and the institutional background, not on which estimator is used or how large the sample is. This section follows the format of Chapters 19 and 20: first look at Aster's numbers or a figure, then state numbered assumptions, theorems, and proof sketches that name the assumptions used at each step, and finally answer whose effect is estimated, who determines the weights, and where the assumptions are most likely to fail. Sections 3.1 to 3.4 cover a single treated unit, Sections 3.5 and 3.6 turn to several treated units and staggered adoption, and Sections 3.7 and 3.8 cover inference and experimental design. The assumption labels A1 to A12 are shared across the whole chapter.

### 3.1 Weights and Identification in Synthetic Control: Matching Unobserved Factor Loadings with a Convex Combination

Start with SC's result on S1. Over the 24 months before the reform the synthetic city tracks the treated city closely: the pre-period RMSPE (the root mean squared prediction error before treatment) is only 0.017, comparable to the standard deviation of the noise itself, 0.02. After the reform the two lines separate, with an average gap of 0.086 against a truth of 0.095 (the figure is in Section 6.3). The DiD of Chapter 19 borrows the counterfactual from the change in the control group; SC borrows the entire path of a weighted set of donors. Whether the fit is close can be seen; the question is under what model "fits closely" implies "estimates accurately."

To turn "if it fits, it is credible" into an identification result, we need a model of how the untreated outcome is generated. What SC assumes is interactive fixed effects, also called a factor model:

::: {.assumption}
**A1 (factor model; Abadie, Diamond and Hainmueller 2010, hereafter ADH, equation (1) with the unit effect $\alpha_i$ added, which the paper folds into the loadings; unnumbered in the paper)** The untreated potential outcome is

$$Y_{it}(0) = \alpha_i + \lambda_t + \boldsymbol\theta_t'\mathbf{W}_i + \boldsymbol\mu_i'\mathbf{f}_t + u_{it}.$$

Here $u_{it}$ has mean zero.

**A2 (no anticipation and no interference; the observation rule of §2.2 of ADH 2010, unnumbered there)** Observed outcomes satisfy the $Y_{it} = Y_{it}(0) + \tau_{it} D_{it}$ of Section 2: the treated city's outcome for $t \leq T_0$ is $Y_{1t}(0)$, and every donor's outcome in every period is its own $Y_{jt}(0)$.

**A3 (perfect pre-period fit; equation (2) of the same paper)** There exist $w^*_j \geq 0$ with $\sum_j w^*_j = 1$ such that $\sum_j w^*_j Y_{jt} = Y_{1t}$ for all $t \leq T_0$ and $\sum_j w^*_j \mathbf{W}_j = \mathbf{W}_1$.

**A4 (regularity conditions; Appendix B of the same paper, unnumbered there)** $u_{it}$ is independent across units and over time and mean-independent of $(\mathbf{W}_i, \boldsymbol\mu_i)$; for some even $p$, $\mathbb{E}|u_{jt}|^p$ is finite; the factors are bounded, $|f_{tk}| \leq f_{\max}$; and the average of $\mathbf{f}_t\mathbf{f}_t'$ over any number of the final pre-periods has smallest eigenvalue no less than $\xi > 0$.
:::

Here $\alpha_i$ is the city's intrinsic level and $\lambda_t$ is a time effect common to all cities, both with the same names and meanings as in Chapter 19; $\mathbf{W}_i$ are observable city covariates (the $\mathbf{W}_i$ of Chapter 21, fixed before treatment), whose coefficients $\boldsymbol\theta_t$ may vary over time (here $\boldsymbol\theta_t$ follows the notation of ADH and is unrelated to the $\theta$ used throughout the book as a structural parameter and to the event-study coefficients $\theta_e$ of Chapter 19); $\mathbf{f}_t$ is a set of $r$ time-varying common factors, with $r$ the number of factors (unrelated to the number of moment conditions in GMM); $\boldsymbol\mu_i$ collects city $i$'s loadings on these factors; and $u_{it}$ is idiosyncratic noise. This model is far more general than the additive two-way fixed effects $\alpha_i + \lambda_t$ of Chapter 19. Additive fixed effects are a special case: when all cities have the same loadings $\boldsymbol\mu_i$ and the covariate coefficients $\boldsymbol\theta_t$ do not vary over time, $\boldsymbol\mu_i'\mathbf{f}_t$ is absorbed into $\lambda_t$, the covariate term is absorbed into $\alpha_i$, and the cities' trends are parallel to one another. In general the loadings can differ, meaning that different cities respond with different magnitudes to the same time-varying common shock, so their trends are naturally not parallel. Aster's DGP is built exactly this way, and this is the root of both the different trend slopes of the two cities in Section 1 and the failure of DiD.

The other three assumptions each handle one thing. The no-anticipation half of A2 requires that the treated city not react before the reform launches because it foresees the reform; otherwise the pre-period is contaminated, $T_0$ is misdated, and the fitted loadings are contaminated too. No interference requires that the donors not be affected by spillovers from the reform: if lower fees draw sellers in from neighboring cities and lead buyers to order across cities, the donors' outcomes are no longer a clean $Y(0)$. A3 is SC's own condition: it requires the treated city's pre-period outcomes and covariates to lie in the convex hull of the donors. A4 is a technical condition, and its most substantive part is that the factors must vary enough at the end of the pre-period: when the factors barely move, the pre-period path cannot reveal the loadings.

::: {.theorem}
**Theorem 1 (bias bound for SC; ADH 2010 Appendix B, pp. 503-504, unnumbered in the paper)** Suppose A1 to A4 hold, and fold $\alpha_i$ into $\boldsymbol\mu_i'\mathbf{f}_t$ as the loading on a factor that always equals 1 (raising $r$ by one). Take $\hat w_j = w^*_j$ and, for $t > T_0$, form the gap

$$\hat\tau_{1t} = Y_{1t} - \sum_{j} \hat w_j Y_{jt}$$

Then $\hat\tau_{1t} - \tau_{1t} = R_{1t} + R_{2t} + R_{3t}$, where $\mathbb{E}[R_{2t}] = \mathbb{E}[R_{3t}] = 0$ and

$$\mathbb{E}|R_{1t}| \leq C(p)^{1/p}\,\frac{f_{\max}^{2} r}{\xi}\,J^{1/p}\max\Big\{\frac{\bar m_p^{1/p}}{T_0^{1-1/p}},\ \frac{\bar\sigma}{T_0^{1/2}}\Big\},$$

where $C(p)$ depends only on $p$, and $\bar m_p$ and $\bar\sigma^2$ are the maxima across donors of the $p$-th absolute moment and of the variance of donor noise, each averaged over the pre-period. The bias is therefore no larger than the right-hand side, and it goes to zero as $T_0$ grows.
:::

The proof sketch has four steps (Appendix B of the paper). Step one (A2): after the reform we observe $Y_{1t}(1)$ for the treated city and $Y_{jt}(0)$ for the donors, so $\hat\tau_{1t} - \tau_{1t} = Y_{1t}(0) - \sum_j w^*_j Y_{jt}(0)$, and the error lies entirely in the counterfactual. Step two (A1): substituting the factor model,

$$\begin{aligned}Y_{1t}(0) - \sum_j w_j Y_{jt}(0) &= \boldsymbol\theta_t'\underbrace{\Big(\mathbf{W}_1 - \sum_j w_j \mathbf{W}_j\Big)}_{\text{match directly}} \\ &\quad + \underbrace{\Big(\alpha_1 - \sum_j w_j \alpha_j\Big) + \Big(\boldsymbol\mu_1 - \sum_j w_j \boldsymbol\mu_j\Big)'\mathbf{f}_t}_{\text{match indirectly}} \\ &\quad + \sum_j w_j (u_{1t} - u_{jt}).\end{aligned}$$

The time effect $\lambda_t$ cancels from the difference because $\sum_j w_j = 1$. Step three (A3 and A4): write the same expression for every pre-period. A3 makes the left-hand side and the covariate term zero, so the pre-period factor matrix times the difference in loadings equals minus the difference in pre-period noise; the eigenvalue condition of A4 makes this equation solvable, so the difference in loadings is a linear combination of pre-period noise. Substituting back into the expression for $t > T_0$ splits the error into three pieces: $R_{1t}$ comes from the donors' pre-period noise, $R_{2t}$ from the treated city's pre-period noise, and $R_{3t}$ is post-period noise. Step four (A4): the noise is independent over time, so $R_{2t}$ and $R_{3t}$ have mean zero; $R_{1t}$ does not, because $\mathbf{w}^*$ itself was chosen with the donors' pre-period noise in view. Apply the Cauchy-Schwarz, Hölder, and Rosenthal inequalities to it in turn: the bounded factors give $f_{\max}^{2} r/\xi$; across the $J$ donors, $p$-th moments give $J^{1/p}$; and averaging over $T_0$ periods gives the powers of $T_0$.

The intuition is that the pre-reform outcome path is itself a fingerprint of the factor loadings: a city's exposure to each common factor is imprinted in the trajectory of how its log GMV has grown and fluctuated with the seasons over the years. If two cities track each other month by month over a long enough pre-period, their exposures are nearly the same, even though we have never observed them directly.

Theorem 1 uses the exact-fit weights; in practice, the weights are determined by a constrained minimum-distance problem. Stack the covariates and pre-period outcomes into a predictor vector $\mathbf{X}_i = (\mathbf{W}_i, Y_{i1}, \dots, Y_{iT_0})$, where the boldface $\mathbf{X}_i$ refers specifically to SC's stacked predictors and is a different object from the scalar running variable $X_i$ of Chapter 21,

$$\begin{aligned}\hat{\mathbf{w}} = \arg\min_{\mathbf{w}}\ & \Big(\mathbf{X}_1 - \textstyle\sum_j w_j \mathbf{X}_j\Big)' \boldsymbol\Omega \Big(\mathbf{X}_1 - \textstyle\sum_j w_j \mathbf{X}_j\Big) \\ \text{s.t. } & \sum_{j} w_j = 1,\ w_j \geq 0,\end{aligned}$$

where $\boldsymbol\Omega$ is a positive-semidefinite predictor weighting matrix. The two constraints are SC's signature: the weights are nonnegative and sum to one, which together are called the simplex constraints. They guarantee that the synthetic city is a convex combination of donors, lying inside the donors' convex hull, with no extrapolation. This is a deliberate choice of SC relative to ordinary regression: the implicit weights of ordinary regression can be negative or larger than one, which amounts to allowing extrapolation into regions with no donor support, whereas SC uses these two constraints to confine itself to the support of the data, at the cost that the convex hull is sometimes too small to reproduce the treated city.

In words: under A1 to A4, the SC gap estimates the period-$t$ effect for this one treated city, $\tau_{1t}$; the weights are set by the pre-period fit, and the researcher chooses only the predictors and the weighting matrix. Two forces pull in opposite directions in the bias bound: the longer the pre-period, the more firmly the loadings are pinned down; the more donors, the easier it is to assemble weights that happen to fit the noise as well, and $J^{1/p}$ is the price of this overfitting, which is one reason not to enlarge the donor pool indiscriminately. S1 has $T_0 = 24$ and $J = 20$. A hand-written SC that uses only the 24 pre-period outcomes and sets $\boldsymbol\Omega$ to the identity uses 9 donors to push the pre-period RMSPE down to 0.013, with an estimate of 0.076. This is lower still than the 0.017 of tidysynth (which matches only 5 predictors, whereas the hand-written SC matches 24 outcomes period by period), and both are below the noise standard deviation of 0.02. That the hand-written SC fits better than the noise itself shows that the weights have fitted the noise too, and the part coming from donor noise causes bias, namely $R_{1t}$. A3 is most likely to fail when the treated unit falls outside the convex hull (S2 in Section 3.4). Proposition 1 of Ferman and Pinto (2021) shows that with an imperfect fit, even as $T_0 \to \infty$, the SC weights do not converge to the weights that recover the loadings; the bias equals the residual imbalance in levels plus the imbalance in loadings times the post-period factors, and in general it does not vanish as the pre-period lengthens. A2 is most likely to fail through spillovers: sellers moving into the treated city depress the donors' GMV, and the gap then counts the spillover as part of the effect.

### 3.2 SC, DiD, and Regression: Different Constraints on One Weighting Problem

On the same S1 data, DiD gives 0.244, the hand-written SC gives 0.076, and the truth is 0.095. The two numbers are so far apart that the two methods look entirely different; Doudchenko and Imbens (2016) point out that they are the same expression with different constraints. Write the treated city's counterfactual in period $t$ as a linear combination of the donors' outcomes in the same period plus an intercept (equation (4.1) of the paper):

$$\hat Y_{1t}(0) = \omega_0 + \sum_{j=2}^{J+1} w_j Y_{jt}.$$

§4.2 of Doudchenko and Imbens (2016) lists four constraints: NO-INTERCEPT ($\omega_0 = 0$), ADDING-UP ($\sum_j w_j = 1$), NON-NEGATIVITY ($w_j \geq 0$), and CONSTANT-WEIGHTS (all $w_j$ equal).

::: {.theorem}
**Proposition 1 (a unified representation of three counterfactuals; Doudchenko and Imbens 2016 equation (4.1) and §4.2, unnumbered in the paper; for the weights in (iii) see Abadie 2021 §4; (iii) and (iv) are derived here)**

(i) SC imposes the first three constraints, with the weights given by the minimum-distance problem of Section 3.1; when the predictors are only the pre-period outcomes and $\boldsymbol\Omega$ is the identity, SC is then least squares for $Y_{1t}$ on the pre-period, without an intercept and subject to the simplex constraints.

(ii) DiD imposes the last three constraints, so that $w_j = 1/J$, and pre-period least squares sets the intercept to $\hat\omega_0 = T_0^{-1}\sum_{t \leq T_0}\big(Y_{1t} - J^{-1}\sum_j Y_{jt}\big)$; the average post-period gap is then the 2x2 DiD of Chapter 19.

(iii) Regressing $Y_{jt}$ on $\tilde{\mathbf{X}}_j = (1, \mathbf{X}_j')'$ by OLS across the donors and predicting at $\tilde{\mathbf{X}}_1$ is, when $\tilde{\mathbf{X}}_0\tilde{\mathbf{X}}_0'$ is invertible, the same as taking $\omega_0 = 0$ and the weights $\mathbf{w}^{\mathrm{reg}} = \tilde{\mathbf{X}}_0'(\tilde{\mathbf{X}}_0\tilde{\mathbf{X}}_0')^{-1}\tilde{\mathbf{X}}_1$ (where $\tilde{\mathbf{X}}_0$ has the donors' $\tilde{\mathbf{X}}_j$ as its columns); these weights satisfy ADDING-UP and balance $\mathbf{X}$ exactly, but they can be negative or larger than one.

(iv) Under A1 and A2, the bias of (ii) equals $\big(\boldsymbol\mu_1 - J^{-1}\sum_j\boldsymbol\mu_j\big)'\big(\bar{\mathbf{f}}_{\mathrm{post}} - \bar{\mathbf{f}}_{\mathrm{pre}}\big)$ plus the analogous covariate term, where $\bar{\mathbf{f}}_{\mathrm{pre}}$ and $\bar{\mathbf{f}}_{\mathrm{post}}$ are the means of the factors over the pre-period and the post-period; differences in $\alpha_i$ are absorbed entirely by the intercept.
:::

The proof sketch has three steps. Step one (algebra): (i) and (ii) substitute their respective constraints into pre-period least squares; under uniform weights, the first-order condition for the intercept gives the average pre-period gap, and the average post-period gap minus it is the double difference (§5.1 of Doudchenko and Imbens 2016). Step two (algebra): the first row of $\tilde{\mathbf{X}}_0$ is $\boldsymbol\iota_J'$, the transpose of the $J$-dimensional vector of all 1s, so $\boldsymbol\iota_J'\mathbf{w}^{\mathrm{reg}}$ equals the first element of $\tilde{\mathbf{X}}_1$, which is 1; moreover $\tilde{\mathbf{X}}_0\mathbf{w}^{\mathrm{reg}} = \tilde{\mathbf{X}}_1$, so the balance is exact, but no constraint rules out negative weights. Step three (A1 and A2): substituting the factor model into (ii), $\lambda_t$ cancels because the weights sum to one, the differences in $\alpha_i$ are time-invariant and absorbed by the intercept, and what remains is the difference in loadings times the change in the factors from the pre-period to the post-period.

In words: all three methods estimate the same $\tau_{1t}$, and they differ only in how the weights are set. The two-way fixed effects model of DiD, $Y_{it} = \alpha_i + \lambda_t + \beta\, D_{it} + \varepsilon_{it}$, assumes additively separable fixed effects: it gives every donor equal weight and absorbs level differences with the intercept, betting that the donors' average loadings equal the treated city's. As soon as different cities respond differently to common shocks, that is, as soon as there is a term like $\boldsymbol\mu_i'\mathbf{f}_t$ with heterogeneous loadings, the bias in (iv) is nonzero. In S1 the treated city's pre-period monthly trend slope exceeds the donor average by 0.014, and DiD overestimates by 0.149, mostly because of this trend gap. SC drops CONSTANT-WEIGHTS and lets the data set the weights, betting that a convex combination can balance the loadings term by term, and this pulls the estimate back to 0.076. Regression is the most flexible, at the price of extrapolation through negative weights; the cross-donor regression of (iii) cannot even be solved in S1: the predictors include 24 pre-period outcomes, there are only 20 donors, and $\tilde{\mathbf{X}}_0\tilde{\mathbf{X}}_0'$ is not invertible. Doudchenko and Imbens relax all four constraints because each is hard to sustain in practice; once they are relaxed, the regression cannot be solved when there are more donors than pre-periods, which is exactly what elastic net regularization is for. In their regression, the treated city is regressed on the donors over the pre-period (equation (4.2) of the paper), and in S1 the 24 periods outnumber the 20 donors, so it can be run. In short: DiD tolerates only time-invariant confounding, while SC can tolerate time-varying confounding that enters through the factor loadings; SC also turns the choice of control group into a data-driven step, with weights that are computed and can be audited. Its NO-INTERCEPT constraint, in turn, is most costly when the treated unit's level differs from that of every donor (Section 3.4).

### 3.3 Synthetic DiD: Weighting Units and Periods at Once

S1 has a third number: synthetic DiD (Arkhangelsky, Athey, Hirshberg, Imbens and Wager 2021, hereafter AAHIW, with the estimator written SDID) gives 0.057, far below the 0.244 of DiD but also below SC, and farther from the truth than SC. To see how the three relate, first write SDID down. It estimates two sets of weights at once: SC-style unit weights $\hat\omega_i$ (making the weighted donors' pre-period path parallel to the treated unit's, up to a constant) and time weights $\hat\nu_t$ (giving more weight to the pre-periods that better predict the post-period), and then runs a weighted two-way fixed effects regression:

$$\big(\hat\tau^{\mathrm{sdid}}, \cdot\big) = \arg\min_{\tau,\, \alpha_i,\, \lambda_t}\sum_{i,t}\big(Y_{it} - \alpha_i - \lambda_t - \tau D_{it}\big)^2\, \hat\omega_i\, \hat\nu_t.$$

Here $\tau$ is the coefficient on $D_{it}$, the treatment effect parameter, not a quantile level; $\hat\omega_i$ is unrelated to the aggregation weights $\omega_e^{g,e'}$ of Chapter 19. For the treated units $\hat\omega_i$ is $1/N_1$ ($N_1$ is the number of treated units, and the total number of units is $N = N_1 + J$), and for each post-period $\hat\nu_t$ is $1/(\mathcal{T} - T_0)$. The donors' unit weights solve equation (2.1) of the paper, an SC problem with an intercept $\omega_0$ and a ridge penalty $\zeta_{\mathrm{sdid}}$, where the intercept plays the same role as the $\omega_0$ of Section 3.2; the pre-period time weights solve equation (2.3). Unlike pure SC, SDID keeps the unit fixed effects, so it does not require an exact match in levels, only parallel trends after weighting; unlike DiD, it does not give all donors and all periods equal weight.

::: {.theorem}
**Proposition 2 (SDID nests DiD and SC; stated here from AAHIW 2021 equations (1.1) to (1.3), (2.4), (2.5), and footnote 5, unnumbered in the paper; no identifying assumption is needed)** For any given unit weights $\omega$ (nonnegative on the donors and summing to one) and time weights $\nu$ (nonnegative on the pre-periods and summing to one), the solution of the regression above is the weighted double difference

$$\begin{aligned}\hat\tau(\omega, \nu) &= \Big(\bar Y_{1,\mathrm{post}} - \sum_{t \leq T_0}\nu_t Y_{1t}\Big) \\ &\quad - \sum_{j}\omega_j\Big(\bar Y_{j,\mathrm{post}} - \sum_{t \leq T_0}\nu_t Y_{jt}\Big),\end{aligned}$$

where $\bar Y_{i,\mathrm{post}}$ is unit $i$'s post-period mean (with several treated units, the first term becomes their average). Hence: (i) with $\omega_j = 1/J$ and $\nu_t = 1/T_0$, it is the unweighted two-way fixed effects regression, that is, DiD; with a single treated unit and a single adoption date, it also equals the TWFE coefficient and the 2x2 DiD of Chapter 19. (ii) Dropping $\alpha_i$ and the time weights and using SC weights without an intercept, the regression reduces to equation (1.3) of the paper, whose solution is the average post-period gap $\bar Y_{1,\mathrm{post}} - \sum_j \hat w_j \bar Y_{j,\mathrm{post}}$, that is, SC.
:::

The proof sketch has three steps. Step one (algebra, weighted least squares with the weights fixed): $\alpha_i$ and $\lambda_t$ can be partialled out one at a time, and $\tau$ is determined only by the weighted means of the treated block and the other three blocks, which gives the expression above, the "adjusted outcomes" of SDID in equation (2.4) together with (2.5) of the paper (the post-period mean minus the time-weighted pre-period mean). Step two (algebra, (i)): with uniform weights substituted in, the adjusted outcome becomes the DiD version in (2.5), the post-period mean minus the pre-period mean, which on a balanced panel coincides term by term with the 2x2 formula; with a single treated unit and a single adoption date, the TWFE coefficient is the same double difference. Step three (algebra, (ii)): without $\alpha_i$, the $\lambda_t$ of each post-period only needs to fit the weighted means of the two groups, the first-order condition for $\tau$ gives the average post-period gap, and the pre-period observations affect only $\lambda_t$; footnote 5 says that the SC unit weights are those of equation (2.1) with the intercept removed and only SC's own tiny regularization kept (penalty $10^{-6}\hat\sigma$, where $\hat\sigma$ is the standard deviation of the first differences of the donors' pre-period outcomes).

Numerical confirmation on Aster: in `synthdid`, `did_estimate` (with both sets of weights uniform) gives 0.244, and its absolute difference from the 2x2 DiD of Section 1 is $6.7 \times 10^{-16}$; the absolute difference between the TWFE coefficient from `fixest` and the 2x2 DiD is $5.6 \times 10^{-16}$. The three numbers are the same number: the DiD and TWFE of Chapter 19 are both special cases of SDID with uniform weights.

Whether SDID is accurate is a question for the asymptotic theory of AAHIW.

::: {.assumption}
**A5 (AAHIW 2021 Assumptions 1 to 4, arXiv v4)** Let $Y_{it} = \alpha_i + \lambda_t + \boldsymbol\mu_i'\mathbf{f}_t + \tau_{it}D_{it} + u_{it}$, with the systematic component and the treatment assignment treated as fixed. (1) The units' noise vectors $(u_{i1}, \ldots, u_{i\mathcal{T}})$ are i.i.d. normal, with covariance-matrix eigenvalues bounded and bounded away from zero; (2) $N_1(\mathcal{T} - T_0) \to \infty$, $J$ and $T_0$ both go to infinity, $T_0/J$ is bounded and bounded away from zero, and $J$ grows faster than $N_1(\mathcal{T} - T_0)\max\{N_1, \mathcal{T} - T_0\}\log^2 J$; (3) the matrix formed by the systematic component has only a few large singular values; (4) the oracle weights (the data-independent weights obtained by replacing the objective functions of equations (2.1) and (2.3) with their expectations) are sufficiently dispersed and fit well enough on the donors and the pre-periods that the bias of the double difference they produce, with respect to the systematic component, is $o\big((N_1(\mathcal{T} - T_0))^{-1/2}\big)$.
:::

::: {.theorem}
**Theorem 2 (asymptotic normality of SDID; AAHIW 2021 Theorem 1, arXiv v4)** Suppose A5 holds and $(N_1(\mathcal{T} - T_0))^{1/2}\log J = o(\zeta_{\mathrm{sdid}}^2)$. Let $\tau$ be the average of $\tau_{it}$ over the treated block. Then

$$\begin{aligned}\hat\tau^{\mathrm{sdid}} - \tau &= \frac{1}{N_1}\sum_{i:\, G_i < \infty}\Big(\bar u_{i,\mathrm{post}} - \mathbf{u}_{i,\mathrm{pre}}'\boldsymbol\psi\Big) \\ &\quad + o_p\big((N_1(\mathcal{T} - T_0))^{-1/2}\big),\end{aligned}$$

where $\boldsymbol\psi$ is the population coefficient from regressing the average post-period noise on the pre-period noise (zero when the noise has no serial correlation). Hence $\hat\tau^{\mathrm{sdid}}$ is asymptotically normal, with an asymptotic variance equal to the oracle variance one would obtain if the systematic component and the noise covariance were known in advance.
:::

The proof sketch has four steps (§4 and §9 of AAHIW 2021, arXiv v4). Step one (Proposition 2): for any weights, $\hat\tau(\omega, \nu) - \tau$ splits into the weighted double difference of the systematic component (the bias) and the weighted double difference of the noise, equation (4.4) of the paper. Step two (A5 (4)): substitute the oracle weights first; the bias is assumed to be of higher order. As long as either the unit weights or the time weights balance the systematic component, the bias is small; this is what the paper calls double robustness. Step three (A5 (1) and (2)): under the oracle weights the noise term is a linear combination of independent rows; the donor side is negligible because the weights are dispersed, and what remains is the treated units' post-period noise minus its linear prediction from their pre-period noise, whose normality gives the limiting distribution. Step four (A5 (2) and (3), and the rate of $\zeta_{\mathrm{sdid}}$): the difference between the estimated and the oracle weights has only a higher-order effect on $\hat\tau$.

In words: SDID estimates the average effect over the treated block (treated units times post-periods), with equal weight on each treated cell; the unit weights are set by the fit to pre-period trends, and the time weights by which pre-periods best predict the post-period. Compared with DiD it adds two layers of localization, and compared with SC it adds a unit fixed effect, so it is not troubled by level differences: on S2 in Section 3.4, the unit fixed effect absorbs a level difference of about 1.0 between the treated city and the weighted donors (the pre-period difference weighted by the time weights). Theorem 2 requires $N_1(\mathcal{T} - T_0)$ to go to infinity, whereas S1 has only 1 treated city, 12 post-period months, and 20 donors, far from this framework. With a single treated unit the jackknife variance is not defined either, so S1 can use only the placebo variance, computed by repeatedly drawing one donor at random to serve as the treated unit, 500 times: SDID gives 0.057 with a standard error of 0.035, and its 95% interval [-0.012, 0.125] covers the truth of 0.095, yet the error of the point estimate is -0.038, larger than that of SC. The time weights are one reason: they put all of their weight on a single period, month 24, so the base period of the double difference is that one month alone, and that month's noise enters the estimate in full. The loadings of the treated city in S1 already lie in the convex hull of the donors' loadings, and SC's pre-period fit has already reached the scale of the noise (Section 3.1), so SDID's two layers of insurance go unused. SDID is not better on every dataset; its value shows when A3 fails, which is the situation of the next subsection.

### 3.4 Imperfect Fit: Augmented SC and Penalized SC

Theorem 1 requires a perfect fit, which is exactly what S2 cannot provide. Aster-Capital's intrinsic level of 6.20 is above that of every donor (the highest is 5.97), so it lies outside the convex hull: when the treated city's characteristics are so extreme that no convex combination can reach them, the convex-weight constraint refuses to extrapolate, and SC can only leave a visible residual imbalance. SC fits along the edge of the hull with only 3 donors, and its pre-period RMSPE is 0.111, 8.7 times that of the same procedure on S1; its post-reform estimate is -0.119, with even the sign wrong, against a truth of 0.079. This convex-hull problem differs from interpolation bias: the latter arises from using very dissimilar donors to interpolate a possibly nonlinear response surface over a wide range, can occur even when the treated unit lies inside the convex hull, and is mitigated by restricting the donor pool to similar units. One way out of the convex-hull problem is to change the metric: in Section 5.1, ADH (2010) match per-capita cigarette sales rather than total sales, and on a per-capita basis the difference in size no longer pushes California out of the convex hull.

![The fee reform in Aster-Capital (S2). Top panel: the gray lines are the 20 donors, the orange line is Aster-Capital's actual log GMV, and the other three lines are the counterfactuals built by SC (blue dashed), augmented SC (green), and SDID (purple dotted); the treated city's intrinsic level is above that of every donor, yet its observed path lies among the gray lines: what falls outside the convex hull is the combination of level and loadings. Bottom panel: the three gaps and the true effect (orange). SC's gap fails to stay at zero even before the reform, first positive and then negative, and falls below zero after the reform; the other two fluctuate around zero before the reform and roughly track the true effect after it.](assets/fig/fig_22_hull.svg)

Augmented SC (Ben-Michael, Feller and Rothstein 2021, hereafter BMFR, with the estimator written ASC) targets exactly this situation. It adds a bias-correction term: an outcome model $\hat m_t(\cdot)$ (commonly a ridge-penalized regression) fits the residual imbalance, which is then subtracted from the estimate,

$$\hat Y_{1t}^{\text{aug}}(0) = \sum_j \hat w_j Y_{jt} + \underbrace{\Big(\hat m_t(\mathbf{X}_1) - \sum_j \hat w_j \hat m_t(\mathbf{X}_j)\Big)}_{\text{correction for residual pre-period imbalance}}.$$

The idea is akin to that of the AIPW of Section 4.5 of Chapter 15: on top of the weighting, an outcome model corrects the residual bias once more; as BMFR themselves put it, the approach is analogous to standard doubly robust estimation.

::: {.assumption}
**A6 (linear factor model and sub-Gaussian noise; BMFR 2021 equation (21), unnumbered in the paper, arXiv v3)** $Y_{it}(0) = \boldsymbol\mu_i'\mathbf{f}_t + u_{it}$, where the factors are bounded and orthonormal over the pre-period ($T_0^{-1}\sum_{t \leq T_0}\mathbf{f}_t\mathbf{f}_t' = \mathbf{I}$); the $u_{it}$ are independent, mean-zero, and sub-Gaussian (with tails no heavier than the normal's) with scale $\sigma$, and uncorrelated with treatment assignment.
:::

::: {.theorem}
**Lemma 1 (ridge ASC is a weighting estimator; BMFR 2021 Lemma 1)** Let $\hat m_t$ be the ridge regression, across the donors, of the period-$t$ outcome on the centered pre-period outcomes $\mathbf{X}_j$, with penalty $\zeta_{\mathrm{ridge}}$; let $\mathbf{X}_0$ have the donors' $\mathbf{X}_j$ as its columns, and let $\hat{\mathbf{w}}^{\mathrm{sc}}$ be the SC weights. Then $\hat Y_{1t}^{\text{aug}}(0) = \sum_j \hat w^{\mathrm{aug}}_j Y_{jt}$, with

$$\hat w^{\mathrm{aug}}_j = \hat w^{\mathrm{sc}}_j + \big(\mathbf{X}_1 - \mathbf{X}_0\hat{\mathbf{w}}^{\mathrm{sc}}\big)'\big(\mathbf{X}_0\mathbf{X}_0' + \zeta_{\mathrm{ridge}}\mathbf{I}\big)^{-1}\mathbf{X}_j,$$

which are the weights that, subject to $\sum_j w_j = 1$, minimize $(2\zeta_{\mathrm{ridge}})^{-1}\|\mathbf{X}_1 - \mathbf{X}_0\mathbf{w}\|^2 + \frac{1}{2}\|\mathbf{w} - \hat{\mathbf{w}}^{\mathrm{sc}}\|^2$. When SC fits perfectly, the two sets of weights coincide.

**Theorem 3 (error bound for ridge ASC; BMFR 2021 Theorem 1)** Under A6, for any $\Delta > 0$, with probability at least $1 - 6\exp(-\Delta^2/2) - \exp\big(-2(\log 2 + J\log 5)\Delta^2\big)$, $\big|Y_{1t}(0) - \sum_j \hat w^{\mathrm{aug}}_j Y_{jt}\big|$ is no larger than the sum of the four terms in equation (24) of the paper: (a) the residual pre-period imbalance, shrunk by ridge along the singular directions of $\mathbf{X}_0$; (b) the additional approximation error from extrapolation; (c) SC's own approximation error, the cost of balancing noisy outcomes rather than loadings; and (d) the post-period noise, $\Delta\sigma\big(1 + \|\hat{\mathbf{w}}^{\mathrm{aug}}\|_2\big)$. As $\zeta_{\mathrm{ridge}} \to \infty$, ASC reverts to SC.
:::

The $\Delta$ in Theorem 3 and in Theorem 7 of Section 3.6 is just an arbitrary positive number in a probability bound, and the $\Delta_{\mathcal{T}}$ of Section 3.7 is an error bound; neither is related to the difference operator of Chapter 19.

The proof sketch has three steps. Step one (the closed form of ridge): the correction term is the inner product of the pre-period imbalance $\mathbf{X}_1 - \mathbf{X}_0\hat{\mathbf{w}}^{\mathrm{sc}}$ with the ridge coefficients, which are linear in the donors' outcomes; combining the two gives $\hat w^{\mathrm{aug}}_j$, and writing out the first-order conditions of the constrained problem verifies that these weights solve the penalized problem. Step two (A6): the error of any weighting estimator equals the loading imbalance times $\mathbf{f}_t$ plus noise, and the loading imbalance can in turn be written as the imbalance in pre-period outcomes minus the imbalance in pre-period noise (equations (22) and (23) of the paper); the former gives (a), and the latter gives (c) through the SC part and (b) through the part passed on by the ridge correction. Step three (the sub-Gaussian tails of A6): tail bounds on these noise terms give (b), (c), (d), and the probability.

ASC fixes the bias by extrapolating; the penalized SC of Abadie and L'Hour (2021) fixes something else. SC on S2 uses only 3 donors: the treated city lies outside the convex hull, and the solution lies on a boundary face of the hull (Section 4.1). When the treated unit lies inside the convex hull, the weights that fit perfectly are often not unique: which set should we pick?

::: {.assumption}
**A7 (general position; the condition of Abadie and L'Hour 2021 Theorem 1, unnumbered in the paper)** In the matrix whose rows are the donors' $(\mathbf{X}_j', 1, \|\mathbf{X}_1 - \mathbf{X}_j\|^2)$, every submatrix formed by any subset of rows has full rank; this rules out two donors with identical predictors.
:::

::: {.theorem}
**Theorem 4 (uniqueness and sparsity of penalized SC; Abadie and L'Hour 2021 Theorem 1, p. 1820)** Penalized SC solves

$$\begin{aligned}\min_{\mathbf{w}}\ & \Big\|\mathbf{X}_1 - \sum_j w_j\mathbf{X}_j\Big\|^2 \\ & + \zeta_{\mathrm{pen}}\sum_j w_j\big\|\mathbf{X}_1 - \mathbf{X}_j\big\|^2, \\ \text{s.t. } & \sum_j w_j = 1,\ w_j \geq 0.\end{aligned}$$

Under A7, when $\zeta_{\mathrm{pen}} > 0$ the solution is unique and has at most $\dim(\mathbf{X}_1) + 1$ nonzero weights.
:::

§2 of Abadie and L'Hour (2021) also points out two limits: as $\zeta_{\mathrm{pen}} \to 0$ the estimator picks, among all SC solutions, the one with the smallest sum of pairwise distances, and as $\zeta_{\mathrm{pen}} \to \infty$ it reduces to nearest-neighbor matching. The proof sketch of Theorem 4 has four steps (see the paper's appendix). Step one (convexity): the objective is convex and the feasible set is the simplex, so a minimum exists and the solution set is convex. Step two (strict convexity of the distance): if two solutions had different fitted values, averaging them would strictly lower the objective, so all solutions have the same $\sum_j w_j\mathbf{X}_j$; $\zeta_{\mathrm{pen}} > 0$ then makes their penalty terms equal as well. Step three (A7): the KKT multipliers are the same for all solutions, so the set of donors that can receive positive weight is the same too; if this set had $\dim(\mathbf{X}_1) + 2$ or more donors, the linear system given by the first-order conditions would contradict the full-rank condition, so there are at most $\dim(\mathbf{X}_1) + 1$ nonzero weights. Step four (A7): on these rows the matrix has full row rank, and the difference between two solutions is mapped to zero on them (the fitted values, the sums, and the penalty terms are all equal), so it must be zero, which gives uniqueness.

In words: ASC and penalized SC both still estimate this treated city's $\tau_{1t}$; they differ in how they respond to the failure of A3. ASC allows negative weights and limited extrapolation, trading a little interpretability for smaller bias; Theorem 3 spells out the cost: extrapolation reduces (a) but increases (b) and (d). On S2, ridge ASC uses 9 negative weights (the smallest -0.069, the largest 0.489, summing to one), and the pre-period RMSPE falls from SC's 0.111 to 0.027. Penalized SC does not solve the extrapolation problem: the penalty only pushes the weights toward donors close to the treated unit, but even the nearest donor's level is below Aster-Capital's, so the penalty cannot reach it. Leave-one-out cross-validation as in §4.1 of Abadie and L'Hour (2021) selects $\zeta_{\mathrm{pen}} = 0.0018$, which still uses SC's 3 donors and gives an estimate of -0.119, almost the same as SC and no sparser; only when the penalty rises to about 0.032 does the number of donors fall to 2. Under 20 different orderings of the donors the solved weights are exactly the same, consistent with Theorem 4; S2 lies outside the convex hull, so the unpenalized solution is already unique, and A7 holds almost surely for continuous data. Finally, with an imperfect fit SC need not beat DiD, a corollary of Proposition 1 of Ferman and Pinto (2021): when SC cannot match the level, NO-INTERCEPT lets the level difference enter the estimate in full, whereas the intercept of DiD and the unit fixed effects of SDID absorb it. The table below reconciles the estimates on S2:

| Estimator | Weights | Unit fixed effects | Estimate | Error |
|---|---|---|---|---|
| SC | Nonnegative, sum to one, 3 donors | None | -0.1193 | -0.1986 |
| Penalized SC (leave-one-out cross-validation) | Same, 3 donors | None | -0.1190 | -0.1983 |
| Ridge ASC | Sum to one, 9 negative weights | None (outcome-model correction) | 0.044 | -0.035 |
| SDID | Unit plus time weights, 19 donors | Yes (level difference of about 1.0) | 0.031 | -0.048 |
| DiD | Uniform | Yes | 0.019 | -0.060 |

The truth is 0.079. No estimator fully fixes the convex-hull problem; of the two closest to the truth, one relies on extrapolation and the other on fixed effects.

### 3.5 Interactive Fixed Effects and Matrix Completion

The second number in Section 1 that did not add up comes from S3. The top panel of the figure below lists eight estimators: the intervals of the DiD family (TWFE at 0.099, and CS at about 0.139 with either control group; the not-yet-treated control group also consists mostly of never-selected cities) all exclude the truth of 0.06, while the intervals of the factor and SC families all include it. The reason is selection: the 15 selected cities have an average loading of 0.208 on the demand-cycle factor, while the 65 never-selected cities have only 0.113; among the selected cities, the loadings are almost unrelated to the wave (a correlation of -0.113, with average loadings of 0.212, 0.210, and 0.203 in the three waves). So the problem lies not in the order of launch but in the fact that the selected and the never-selected cities respond differently to the same cyclical shock, and the parallel trends of Chapter 19 fails head-on. With a single treated city SC sidesteps the factors; with many cities and many treatment dates, another route is to estimate them directly: the interactive fixed effects estimator (hereafter IFE) estimates a finite number of factors and loadings, and matrix completion (hereafter MC) requires only that the counterfactual matrix be approximately low-rank.

![The staggered support program in S3. Top panel: point estimates and 95% intervals of the eight estimators, with the DiD family in orange, the factor and SC families in blue, and the truth of 0.06 as a dashed line. Bottom panel: across 24 redrawn datasets, the number of times each estimator's interval covers the truth; TWFE and CS never do, IFE and multisynth always do, and the two matrix completion variants cover it only 6 and 14 times.](assets/fig/fig_22_staggered.svg)

::: {.assumption}
**A8 (functional form, strict exogeneity, and low rank; Liu, Wang and Xu 2024 Assumptions 1 to 3, arXiv v3)** (1) $Y_{it}(0) = \alpha_i + \lambda_t + L_{it} + u_{it}$ with $L_{it} = \boldsymbol\mu_i'\mathbf{f}_t$; (2) $u_{it}$ is independent of the treatment status of all units in all periods and of $(\alpha_j, \lambda_s, L_{js})$; (3) the matrix $\mathbf{L} = [L_{it}]$ has rank $r$ much smaller than $\min(N, \mathcal{T})$.

**A9 (Athey, Bayati, Doudchenko, Imbens and Khosravi 2021 Assumptions 1 and 2, arXiv v5)** (1) The noise is independent of $\mathbf{L}$, and its elements are mutually independent and sub-Gaussian; (2) given $\mathbf{L}$, the units' adoption dates are independent of one another and of the noise.
:::

Here $\mathbf{L}$ is a low-rank matrix, unrelated to the HAC lag truncation $L$ of Chapter 8 and to the upper bound $L$ of the event window in Chapter 19. A8 (2) and A9 (2) are the key: they allow treatment to be arbitrarily correlated with the loadings, so a platform that deliberately picks the cities with high cyclical exposure does not by itself violate them (although taking the top 15 by rank ties the cities' adoption dates to one another, so A9 (2) holds only approximately); what they forbid is treatment correlated with idiosyncratic shocks, for example a platform that selects a city because it sees an unexpected fall in that city's GMV this month, which is the kind of selection behind the Ashenfelter dip of Chapter 19. A8 (1) also rules out anticipation and spillovers; (3) is a feasibility condition: if the loadings are to be estimated from the data, there cannot be too many factors.

::: {.theorem}
**Theorem 5 (consistency of the IFE counterfactual estimator; Liu, Wang and Xu 2024 Proposition 3, with regularity conditions of the same kind as Bai 2009 Assumptions A to D)** Suppose A8 and the regularity conditions of Bai (2009) hold (the second-moment matrices of the factors and loadings converge to positive-definite matrices, the noise is weakly dependent and has eighth moments, and so on). Estimate $\alpha_i$, $\lambda_t$, and the $r$ factors and loadings using only the untreated cells, impute $\hat Y_{it}(0)$ for the treated cells, and average $Y_{it} - \hat Y_{it}(0)$ over the treated cells. Then as $N, \mathcal{T} \to \infty$, this average converges in probability to the average effect over the treated cells.

**Theorem 6 (error bound for matrix completion; Athey, Bayati, Doudchenko, Imbens and Khosravi 2021 Theorem 2)** Let $\mathcal{O}$ be the set of untreated cells. MC solves

$$\begin{aligned}\hat{\mathbf{L}} = \arg\min_{\mathbf{L}}\ & \frac{1}{|\mathcal{O}|}\sum_{(i,t) \in \mathcal{O}}\big(Y_{it} - L_{it}\big)^2 + \zeta_{\mathrm{nn}}\|\mathbf{L}\|_* \\ \text{s.t. } & \|\mathbf{L}\|_{\max} \leq L_{\max},\end{aligned}$$

where $\|\mathbf{L}\|_*$ is the nuclear norm (the sum of the singular values), $L_{\max}$ is an upper bound on the absolute values of the elements, and the paper folds the two-way fixed effects into $\mathbf{L}$ as well. Under A9, if $\mathrm{rank}(\mathbf{L}) = r$, $\mathcal{T}$ is at least some constant multiple of $\log(N + \mathcal{T})$, and $\zeta_{\mathrm{nn}}$ is of the order given in the paper, then with probability at least $1 - 2(N + \mathcal{T})^{-2}$, $\|\mathbf{L} - \hat{\mathbf{L}}\|_F/\sqrt{N\mathcal{T}}$ is bounded by a quantity that increases with $r$, $L_{\max}$, the noise scale, and $\log(N + \mathcal{T})$, and that goes to zero only when $Np_c^2$ and $\mathcal{T}p_c^2$ both grow faster than the logarithmic factors; $p_c$ is the smallest probability that a unit is still untreated in the last period.
:::

The proof sketch has four steps. Step one (A8 (1) and (2)): the observations in the untreated cells are $Y_{it}(0)$, and their noise is unrelated to who is treated, so fitting the model to the untreated cells alone brings in no selection bias; this is the logic of the BJS imputation of Section 4.5 of Chapter 19, with the model changed from $\alpha_i + \lambda_t$ to $\alpha_i + \lambda_t + L_{it}$. Step two (A8 (3) and Bai's regularity conditions): the factor space estimated from the untreated cities is consistent (Bai 2009 Proposition 1: the space spanned by the estimated factors converges to the true factor space); the loadings of the selected cities are estimated from their own pre-periods, more accurately the longer the pre-period, which is why Remark 2 of Xu (2017) requires the number of control units and $T_0$ to go to infinity together (that paper has no numbered theorems). Step three (A8 (1)): the imputed values are subtracted from the observations and averaged over the treated cells, and consistency carries over from the first two steps. Step four (A9): MC does not estimate $r$ but uses a nuclear-norm penalty instead, set above the spectral norm of the noise; it then uses the structure of the adoption dates to control how well the observed cells reflect the whole matrix, which yields an RMSE bound (Appendix A.2 of the paper).

In words: IFE and MC both estimate the average effect over the treated cells, with equal weight on each of the 240 treated city-months in S3; the weights do not appear explicitly but are implied by how the model imputes $Y(0)$. IFE stands or falls with the number of factors: the cross-validation of `fect` selects $r = 2$, matching the two factors of the DGP, with an estimate of 0.051 and an interval [0.028, 0.073] that contains the truth; across 24 redrawn datasets it covers the truth every time, selecting $r = 2$ in 23 of them. The MC theorem must be read with care. Theorem 6 controls the average error over the whole matrix, of which the 240 treated cells are only a small block; the nuclear-norm penalty shrinks the singular values toward zero, so the imputation is biased toward simpler structure, that is, toward the two-way fixed effects fit, and the bias on the treated block does not vanish along with the average error. On the S3 dataset shown in this chapter, MC refitted with a tightened convergence tolerance gives 0.067, with an interval [0.047, 0.087] that contains the truth, removing 0.82 of the TWFE bias. But this dataset was picked from 24 draws that change only the random seed, by the criterion that "every conclusion holds under both MC procedures and the default MC covers most reliably," so it is no accident that it is the draw most favorable to MC (it ranks first under both MC procedures). Across the 24 draws MC overestimates the truth every time, with an average estimate of 0.080, removing on average only 0.55 of the TWFE bias (the mean of the per-draw shares); the refitted interval covers the truth in only 14 of them, and the default settings of `fect` in only 6. The more serious violation in S3 concerns the independence of the noise: the idiosyncratic shocks contain a persistent AR(1) component with coefficient 0.8 and stationary standard deviation 0.04, twice the 0.02 of the independent part. Assumption C of Bai allows weak serial correlation, so A8 still holds, consistent with IFE covering the truth in all 24 draws; A9 (1) requires the elements of the noise to be independent, which fails here, so Theorem 6 gives no guarantee at all for S3. The performance of MC on S3 is only an empirical result, which is one more reason not to generalize from it.

### 3.6 SC under Staggered Adoption

IFE and MC first estimate a global model; another route keeps SC's transparent weights and builds a synthetic control for each selected city. The multisynth row in the figure of Section 3.5 follows this route: an estimate of 0.053, close to the truth, but an interval [-0.018, 0.135] wide enough to cover zero. The difficulty is that good fits for each city need not average into a good fit, and vice versa.

::: {.assumption}
**A10 (Ben-Michael, Feller and Rothstein 2022 Assumptions 1, 2b, and 3, arXiv v2)** (1) No anticipation: $Y_{it}(g) = Y_{it}(0)$ for $t < g$; (2) $Y_{it}(0) = \boldsymbol\mu_i'\mathbf{f}_t + u_{it}$, where the factors are bounded and $u_{it}$ has mean zero, is independent across units and over time, and is independent of the adoption date $G_i$; (3) $u_{it}$ is sub-Gaussian.
:::

On selection, A10 (2) is the same as A8 (2): adoption dates may be arbitrarily correlated with the loadings, so a platform picking cities by cyclical exposure is fully allowed; what is not allowed is picking them by idiosyncratic shocks.

::: {.theorem}
**Theorem 7 (error bound for staggered SC; Ben-Michael, Feller and Rothstein 2022 Theorem 2)** Suppose A10 holds, each selected city fits weights on the simplex using the same number of pre-adoption periods, the weights are independent of the post-period noise, and each city's pre-period factor matrix is nonsingular and normalized (the paper's regularity conditions, of the same kind as A4). Let $q^{\mathrm{sep}}$ be the root mean square of the cities' pre-period fit errors and $q^{\mathrm{pool}}$ the pre-period fit error of the selected cities' average path. Then for any $\Delta > 0$, the error of the estimated average effect at event time $e$ is, with probability at least $1 - 6\exp(-\Delta^2/2)$, no larger than $\|\bar{\mathbf{f}}_e\|_2\, q^{\mathrm{pool}} + S_e\, q^{\mathrm{sep}}$ plus two terms for approximation error and noise; here $\bar{\mathbf{f}}_e$ is the average of the (projected) factors across adoption dates and $S_e$ is their dispersion, with $S_e = 0$ when all cities adopt at the same time.
:::

The proof sketch has three steps (see the paper's appendix). Step one (A10 (1)): pre-adoption outcomes are $Y(0)$, and the estimated average effect at event time $e$ is the average of the cities' SC gaps, which also equals a single SC gap for the "average selected city" (the two readings of equation (5) of the paper). Step two (A10 (2)): substituting the factor model, the loading imbalance is expressed through the pre-period fit errors; when averaging over adoption dates, the common part is multiplied by the fit of the average path, $q^{\mathrm{pool}}$, and the part that disperses because adoption dates differ is multiplied by the city-by-city fit, $q^{\mathrm{sep}}$. Step three (A10 (3)): the approximation error from pre-period noise and the post-period noise are controlled with sub-Gaussian tail bounds.

Theorem 7 leads directly to partially pooled SC (`augsynth::multisynth`): it minimizes $\pi_{\mathrm{pool}}(\tilde q^{\mathrm{pool}})^2 + (1 - \pi_{\mathrm{pool}})(\tilde q^{\mathrm{sep}})^2$ plus a ridge penalty on the weights (equation (6) of the paper, where the tilde denotes fit errors normalized by their values under city-by-city SC; the paper writes the pooling parameter as $\nu$, and we write it as $\pi_{\mathrm{pool}}$ instead to avoid confusion with the time weights of SDID). Here $\pi_{\mathrm{pool}} = 0$ is city-by-city SC and $\pi_{\mathrm{pool}} = 1$ is a single SC for the average path. SDID itself handles only simultaneous adoption; the approach in §8 of the appendix to AAHIW 2021 (arXiv v4) is to split the sample by adoption date: each wave forms a simultaneous-adoption subsample with the never-selected cities alone, SDID is run once on each, and the results are averaged with weights equal to each wave's number of treated cells.

In words: both estimate the average post-adoption effect for the selected cities. In S3 the effect is 0.06 everywhere, so any weighted average of the $ATT(g,t)$ of Chapter 19 is 0.06, and the comparison here concerns bias alone. By default, multisynth averages only the first 12 months after each city's adoption, with equal weight per city (180 of the 240 treated city-months), while SDID by wave uses every post month. The weights of multisynth trade off, through $\pi_{\mathrm{pool}}$, a close fit for each city against a close fit for the average path; on S3 it selects $\pi_{\mathrm{pool}} = 0.343$. By wave, the three averages are 0.025, 0.061, and 0.073: the SC for an individual city is very noisy, and only the average comes close to the truth. Across the 24 draws its interval covers the truth every time, with an average estimate of 0.062. SDID by wave gives 0.054 with interval [0.031, 0.076], and 0.046, 0.049, and 0.073 for the three waves; although its bootstrap resamples all 80 cities, the treated side has only 15 selected cities, so the reliability of the interval is limited. Both methods also require cities with similar loadings in the donor pool, and in S3 these are provided by the 20 cities, chosen at random, that cannot be selected because of existing contracts. What is truly violated is the independence over time in A10 (2): the persistent shocks of S3 leave the bound of Theorem 7 without a guarantee here, and the fact that multisynth covers the truth in all 24 draws is an empirical result for this design; its intervals come from the bootstrap and do not rely on Theorem 7.

### 3.7 Inference: Permutation, Conformal Inference, and Prediction Intervals

With only one treated unit, where does a standard error come from? The in-space placebo test of Section 4.2 takes each donor in turn as the treated unit, and the permutation p-value for S1 is the smallest that 21 units can give, 0.048. It answers whether the treated city's gap is the most extreme within a family of fake effects; its resolution is capped by the number of units, and it requires the units to be exchangeable. The two routes of this subsection, conformal inference and prediction intervals, move the source of randomness from the unit dimension to the time dimension.

![Monthly effects in S1: the dots are the SC estimates, the vertical bars are the `scpi` prediction intervals (with 5% error allowed for each of the in-sample and out-of-sample parts, for coverage of at least 90%), and the orange line is the true effect path. All 12 monthly intervals cover that month's truth.](assets/fig/fig_22_scpi.svg)

::: {.assumption}
**A11 (Chernozhukov, Wüthrich and Zhu 2021 Assumptions 1 to 3, arXiv v10)** (1) There is a proxy $P_t$ with $\mathbb{E}[P_t] = \mathbb{E}[Y_{1t}(0)]$ such that $Y_{1t}(0) = P_t + u_t$, where $u_t$ has mean zero, is stationary, and has a distribution unaffected by the reform; (2) one of two alternatives: $u_t$ is i.i.d. (Assumption 2.1 of the paper), or $u_t$ is stationary and strong mixing (Assumption 2.2); (3) the estimator of $P_t$ is consistent under the null: there exist $\Delta_{\mathcal{T}}$ and $\Delta'_{\mathcal{T}}$ going to zero such that, with probability at least $1 - \Delta'_{\mathcal{T}}$, both the root mean squared error of the estimator and its pointwise errors in the post-period are no larger than $\Delta_{\mathcal{T}}$.
:::

::: {.theorem}
**Theorem 8 (conformal inference; Chernozhukov, Wüthrich and Zhu 2021, with (i) Theorem D.1 of Appendix D and (ii) Theorem 1)** Test $H_0$: the post-period effect path equals a given value (for example, zero throughout). Under $H_0$, treat all $\mathcal{T}$ periods as untreated and fit $P_t$, obtaining residuals $\hat u_t$; form the statistic $\mathrm{TS}(\hat{\mathbf{u}}) = (\mathcal{T} - T_0)^{-1/2}\sum_{t > T_0}|\hat u_t|$ and the p-value

$$\hat p_{\mathrm{cf}} = \frac{1}{|\Pi|}\sum_{\pi \in \Pi}\mathbf{1}\{\mathrm{TS}(\hat{\mathbf{u}}_\pi) \geq \mathrm{TS}(\hat{\mathbf{u}})\},$$

where $\Pi$ is either the set of all permutations, $\Pi_{\mathrm{all}}$, or the set of all cyclic shifts, $\Pi_{\to}$. (i) If the data are exchangeable under the null and the estimator of $P_t$ is invariant to permutations of the data, then $\Pr(\hat p_{\mathrm{cf}} \leq \alpha) \leq \alpha$: the test is exactly valid in finite samples, whatever the choice of $\Pi$. (ii) If A11 holds, with $\Pi_{\mathrm{all}}$ paired with Assumption 2.1 and $\Pi_{\to}$ with Assumption 2.2, the density of the statistic is bounded, and the number of post-periods is fixed, then $|\Pr(\hat p_{\mathrm{cf}} \leq \alpha) - \alpha| \leq C\big(\big((\mathcal{T} - T_0)/T_0\big)^{1/4}\log\mathcal{T} + \Delta_{\mathcal{T}} + \sqrt{\Delta_{\mathcal{T}}} + \Delta'_{\mathcal{T}}\big)$, where the constant $C$ does not depend on $\mathcal{T}$.

**Lemma 2 (prediction intervals; Cattaneo, Feng and Titiunik 2021 Lemma 1, arXiv v4)** Split $\hat\tau_{1t} - \tau_{1t}$ into an in-sample part (the error from estimating the weights on the pre-period) and an out-of-sample part (new post-period noise and misspecification error). If $[a_{\mathrm{in}}, b_{\mathrm{in}}]$ covers the former with conditional probability at least $1 - \alpha_{\mathrm{in}}$ and $[a_{\mathrm{out}}, b_{\mathrm{out}}]$ covers the latter with conditional probability at least $1 - \alpha_{\mathrm{out}}$, then $\mathrm{PI}_{1t} = [\hat\tau_{1t} - b_{\mathrm{in}} - b_{\mathrm{out}},\ \hat\tau_{1t} - a_{\mathrm{in}} - a_{\mathrm{out}}]$ covers, with conditional probability at least $1 - \alpha_{\mathrm{in}} - \alpha_{\mathrm{out}}$, the effect $\tau_{1t}$ (the paper adds a probability slack that goes to zero as $T_0$ grows).
:::

Here $\alpha$ without a subscript is the significance level, and with the subscripts in and out it denotes the error level of each part; none of these is related to the unit effect $\alpha_i$. The symbol $\pi$ denotes a single permutation and is unrelated to the $\pi_{\mathrm{pool}}$ of Section 3.6. The proof sketch of Theorem 8 has three steps. Step one ($H_0$): the null hypothesis specifies the effects in the treated periods, and once these are subtracted, in all $\mathcal{T}$ periods the treated city is at $Y(0)$, so all periods can be used for the fit; residuals from fitting on the pre-period and extrapolating to the post-period lack this symmetry, and this step is the key to the method. Step two (exchangeability and permutation invariance, or A11 (2)): SC fitted on all periods is invariant to permutations of the data, so under exchangeability the distribution of the residual vector is invariant to any permutation, and the rank of the observed statistic among the $|\Pi|$ permuted statistics is uniformly distributed, which proves (i); i.i.d. noise allows arbitrary reordering, whereas stationary weak dependence preserves only the block structure in time, so only cyclic shifts are allowed. Step three (A11 (3)): estimation error moves $\hat u_t$ away from $u_t$, and consistency keeps this deviation at the order of $\Delta_{\mathcal{T}}$, at the cost of the bound in (ii). The proof sketch of Lemma 2 has two steps. Step one (algebra): if $\tau_{1t}$ falls outside $\mathrm{PI}_{1t}$, at least one part of the error must fall outside its own interval. Step two (the two coverage conditions): the union bound gives this event a probability of at most $\alpha_{\mathrm{in}} + \alpha_{\mathrm{out}}$. The difficulty lies in the two intervals themselves: the in-sample part is approximated by simulation following Theorems 1 and 2 of the paper, and the out-of-sample part requires an assumption about the distribution of the post-period noise; this chapter uses a Gaussian bound.

In words: the conformal test asks whether this post-period effect path is compatible with zero; its randomness comes from time rather than from units, so it does not need many donors. What it needs is A11 (1): the reform only shifts the outcome and does not change the distribution of the shocks; if the reform also changes volatility, the p-value tests the joint hypothesis of a zero effect and an unchanged distribution. The shocks in S1 are i.i.d. by construction and satisfy Assumption 2.1, and 9,999 random permutations give $\hat p_{\mathrm{cf}} = 0.082$; assuming only Assumption 2.2 and using the 36 cyclic shifts, the p-value is 0.222 (0.111 when the statistic uses the sum of squares instead). At the 10% level the former rejects and the latter does not; the difference lies not in the data but in whether you are willing to assume independent shocks or only weakly dependent ones. Under cyclic shifts the p-value moves in steps of 1/36 and cannot hit a level of exactly 10%, so the test is slightly conservative; two further features cost power: the fit under the null must absorb an effect that ramps up month by month, and the statistics of adjacent shifts are close to one another. Prediction intervals answer a different question: where each month's effect lies. On S1, allowing 0.05 error for each of the in-sample and out-of-sample parts, all 12 monthly intervals cover that month's truth, with an average length of 0.141; intervals that use only the in-sample part cover just 9 months, so the out-of-sample uncertainty cannot be left out. When the shocks are heavy-tailed or serially correlated, the out-of-sample Gaussian bound is too narrow, and coverage may fall below the nominal level.

### 3.8 SC for Experimental Design

The previous seven subsections all repair things after the fact: treatment has already happened, and we look for a credible set of controls. S4 turns the problem around: the pilot has not started yet, and the platform can choose the pilot cities itself. Pick 5 at random, and with bad luck you get a set of cities that are atypical in level or trend. Abadie and Zhao (published online in 2026 in the Review of Economics and Statistics) propose picking with SC: choose a set of pilot cities whose weighted average reproduces the historical average of all 40 cities, then build a synthetic control from the remaining cities that also reproduces this average; both sets of weights lie on the simplex, and the two sides share no cities (the constraints of design problem (6) in the paper).

![The pilot design in S4. Top panel: the average log GMV of the 40 cities (gray), the synthetic average of the 5 pilot cities (orange), and the synthetic control from the remaining cities (blue dashed); the three lines coincide over the 24 months of history, and the orange line rises after the pilot starts. Bottom panel: the distributions of the estimates of the three approaches over 500 replicated pilot periods, with the truth of 0.04 as a dashed line; the RMSE is 0.0087 for the synthetic design, 0.077 for random selection followed by DiD, and 0.212 for random selection followed by a comparison of means.](assets/fig/fig_22_design.svg)

::: {.assumption}
**A12 (Abadie and Zhao Assumptions 1 to 3, arXiv v5)** (1) The potential outcomes in both the untreated and the treated state follow a factor model, and the factor coefficients in the treated state may differ; (2) over the fitting period ($T_E$ months of history), the smallest eigenvalue of the factor matrix divided by $T_E$ is no less than $\xi > 0$, and the noise is sub-Gaussian (with scale $\sigma$) and independent over time; (3) the weights on the pilot side, $\mathbf{w}^{\mathrm{tr}}$, and on the control side, $\mathbf{w}^{\mathrm{co}}$, both exactly reproduce the average predictors of all cities over the fitting period.
:::

::: {.theorem}
**Theorem 9 (bias bound for the synthetic design; Abadie and Zhao Theorem 1, arXiv v5 equation (15))** Suppose A12 holds, $\mathbf{w}^{\mathrm{tr}}$ and $\mathbf{w}^{\mathrm{co}}$ solve design problem (6) of the paper, and the estimator is $\hat\tau_t = \sum_i w^{\mathrm{tr}}_i Y_{it} - \sum_i w^{\mathrm{co}}_i Y_{it}$. Then for every pilot period $t$,

$$\big|\mathbb{E}[\hat\tau_t - \tau^{\mathrm{pop}}_t]\big| \leq \frac{f_{\max}\,(f^{\mathrm{I}}_{\max} + f_{\max})\,r}{\xi}\,\sqrt{2\log(2N)}\,\frac{\sigma}{\sqrt{T_E}},$$

where $f^{\mathrm{I}}_{\max}$ is the upper bound on the factor coefficients in the treated state, and $N$ is the total number of candidate cities (written $J$ in the paper). If A12 (3) holds only approximately (Assumption 4 of the paper: a fit error no larger than $d$), the right-hand side gains a term proportional to $d$ (equation (16) of the paper).
:::

The proof sketch has three steps (see the paper's appendix). Step one (A12 (3)): the weighted averages of both sides equal the overall average over the fitting period, so all that remains of $\hat\tau_t - \tau^{\mathrm{pop}}_t$ is each side's loading imbalance times the corresponding factor coefficients, plus noise. Step two (A12 (1) and (2)): as in step three of the proof of Theorem 1, solve for the loading imbalance as a linear combination of the fitting-period noise, with the eigenvalue condition guaranteeing solvability. Step three (A12 (2)): the weights are chosen with the fitting-period noise in view, and the bias comes from the correlation between the two; the maximum of sub-Gaussian noise over $N$ cities gives $\sqrt{2\log(2N)}$, and the longer the fitting period, the smaller this term.

In words: the synthetic design estimates the average effect over all 40 cities, with equal weight on each city, not the pilot cities' own effect; the weights are set from history before the experiment begins, when no outcomes exist yet, so the researcher has no opportunity to choose weights based on outcomes. The design for S4 picks Metro-01, 09, 14, 16, and 27 (weights 0.185, 0.203, 0.273, 0.082, and 0.256), with a historical fit error of 0.0034, which is not zero, so what applies is equation (16) with its $d$ term. Over 500 replicated pilot periods, the RMSE of the synthetic design is 0.0087, while picking 5 cities at random and running DiD gives 0.077, 8.8 times as large. Comparing means directly after random selection gives an RMSE as high as 0.212, mostly noise from differences in intrinsic levels across cities, which DiD removes, so the factor of 8.8 is the fair comparison. The bias of the synthetic design (0.0062) is conditional on this particular history and is not the expected bias bounded by Theorem 9, but its sources are of the same kind, the overfitting term and the approximate-fit term: with only 24 fitting months and 40 candidate cities, overfitting the noise is its main risk. Its point of failure is A12 (1): if the pilot period brings a shock not seen in the history, or the effect varies with some city characteristic that was not fitted, the synthetic average no longer represents the whole.

To summarize: identification in all eight subsections rests on assumptions about the counterfactual, with the factor structure, the fit, the noise, and the selection mechanism each accounting for a part; very few of them can be tested directly, and what the pre-period fit, placebos, and resampling can test is only the shadows they cast.

## 4 Estimation

Identification has made clear what SC estimates; this section covers how to estimate it: Section 4.1 covers SC's weight optimization, Section 4.2 its permutation inference, and Section 4.3 lists implementations of the estimators of Section 3.

### 4.1 SC: Weight Optimization and Fit

The estimation core of SC is the simplex-constrained minimum-distance problem of Section 3.1. In practice it proceeds in a few steps: first choose the predictors (covariates plus the values of several pre-period outcomes; a common practice is the pre-period mean plus lags at a few specific months), then optimize the weights $\hat{\mathbf{w}}$ over the pre-period so that the synthetic city is as close as possible to the treated city on these predictors, and then use these weights to construct the synthetic path $\sum_j \hat w_j Y_{jt}$ over all periods; the treated city's actual path minus this path is the period-by-period gap $\hat\tau_{1t}$. When the treated city's predictors lie outside the convex hull of the donors, the solution of this constrained quadratic problem lies on a face of the hull's boundary, the weights are sparse, and the number of nonzero donors usually does not exceed the number of predictors plus one. This sparsity does not come from $\ell_1$ shrinkage (on the simplex, $\sum_j |w_j| = \sum_j w_j = 1$ is always constant) but from the vertex geometry of the simplex: the optimum of the quadratic objective tends to land on a low-dimensional face with only a few nonzero coordinates. When the predictors lie inside the convex hull, the weights that fit exactly form a whole set, and the optimizer returns only one point of it, which need not be sparse. When the weights are sparse, the synthetic city is pieced together from a few nameable donors, which is more transparent than a twenty-city average. In R, `Synth` and `tidysynth` are the common implementations, and `MSCMT` fixes some numerical issues in the optimization of `Synth`.

The weighting matrix $\boldsymbol\Omega$ is usually left to the data as well: `Synth` and `tidysynth` run a nested optimization, in which the inner loop solves for the weights given $\boldsymbol\Omega$ and the outer loop picks the diagonal elements of $\boldsymbol\Omega$ to minimize the fit error of the pre-period outcomes. If every pre-period outcome is used as a predictor, the remaining covariates receive zero weight in this optimization (Kaul, Klößner, Pfeifer and Schieler 2022); the hand-written SC of Section 3.1 likewise uses only the 24 pre-period outcomes and simply sets $\boldsymbol\Omega$ to the identity.

### 4.2 SC Inference: Permutation and Placebo

For inference with a single treated unit, the classic approach relies not on analytic standard errors but on permutation and placebos; the conformal test and the prediction intervals of Section 3.7 are two other routes.

The core is the in-space placebo: take each donor city in turn and pretend it is the "treated unit," run SC on it, compute a gap curve, and obtain a whole family of "fake effects." If the treated city truly has an effect, its post-period gap should look abnormally large within this family. To quantify this "abnormality," use the post/pre MSPE ratio:

$$\mathcal{R}_i = \frac{\mathrm{MSPE}^{\text{post}}_i}{\mathrm{MSPE}^{\text{pre}}_i}, \qquad \mathrm{MSPE}^{\text{pre}}_i = \frac{1}{T_0}\sum_{t \leq T_0}\big(Y_{it} - \hat Y_{it}(0)\big)^2.$$

The denominator is the pre-period fit error and the numerator is the post-period gap. A unit that truly has an effect fits well in the pre-period (small denominator) and has a large post-period gap (large numerator), so its ratio is large. Rank the treated unit's $\mathcal{R}_1$ among all units' $\{\mathcal{R}_i\}$; the permutation p-value is the share of units whose ratio is no smaller than the treated unit's, $\hat p_{\mathrm{perm}} = \frac{1}{J+1}\sum_{i=1}^{J+1} \mathbf{1}\{\mathcal{R}_i \geq \mathcal{R}_1\}$. The ratio is already normalized by the pre-period fit: a placebo with a poor pre-period fit has a large denominator, so other things equal its ratio is smaller, and ADH use the ratio precisely so that no separate exclusion threshold is needed. Only when gap curves are compared directly do placebos with a poor pre-period fit have to be dropped first (the figure in Section 6.3). The permutation p-value of Section 6.3 is computed from the ratios without any exclusion.

Another diagnostic is the in-time placebo: set the treatment date to some period before the actual reform and re-estimate. A design that holds up should show no gap around this fake, not-yet-reformed date. It tests not the size of the effect but the timing: if a gap appears before the actual reform, the fit is capturing something else, not the policy.

### 4.3 Implementations of the Modern Estimators

The estimators of Sections 3.3 to 3.8 all have ready-made implementations, but the default settings are not always appropriate. The implementations and settings used in this chapter are as follows:

| Estimator | Implementation | How the weights or tuning are set | Inference |
|---|---|---|---|
| SC | `tidysynth`; hand-written `sc_weights` | $\boldsymbol\Omega$ by nested optimization, or $\boldsymbol\Omega = \mathbf{I}$ | In-space and in-time placebos |
| Penalized SC | Hand-written `penalized_sc` | $\zeta_{\mathrm{pen}}$ chosen by leave-one-out cross-validation | Not done |
| Ridge ASC | `augsynth`, `progfunc = "ridge"` | $\zeta_{\mathrm{ridge}}$ chosen by the package's cross-validation | Conformal |
| SDID | `synthdid` | Unit and time weights, $\zeta_{\mathrm{sdid}}$ from equation (2.2) of the paper | Placebo variance |
| IFE | `fect`, `method = "ife"` | $r$ cross-validated between 0 and 4 | Bootstrap |
| MC | `fect`, `method = "mc"` | $\zeta_{\mathrm{nn}}$ cross-validated, then refitted with a tightened tolerance | Bootstrap |
| multisynth | `augsynth::multisynth` | $\pi_{\mathrm{pool}}$ by the heuristic of §4.2 of Ben-Michael, Feller and Rothstein (2022) | Bootstrap |
| SDID by wave | `synthdid`, estimated wave by wave and then weighted | Each wave against the never-selected cities | City bootstrap |
| Conformal test | Hand-written `conformal_p` | SC fitted on all periods under the null | Random permutations or cyclic shifts |
| Prediction intervals | `scpi` | Simplex weights, Gaussian bound | In-sample plus out-of-sample |
| Synthetic design | `sc_weights` plus local search | Design problem (6) of the paper | Replicated pilot periods |

Four settings materially change the numbers. First, with a single treated unit the jackknife and bootstrap of `synthdid` return no result, so the SDID intervals for S1 and S2 can come only from the placebo variance, which requires the treated city and the donors to have identically distributed noise (Section 5.2). Second, the penalty of ridge ASC is set by cross-validation: on S1 it selects 263.6, the correction term is almost zero, and ASC estimates 0.076, the same as the hand-written SC, with the smallest of its 9 negative weights only $-8.7 \times 10^{-6}$; on S2 it selects 0.752, and only then does the correction term matter. This is exactly the last sentence of Theorem 3: the larger $\zeta_{\mathrm{ridge}}$, the closer ASC is to SC. Third, after cross-validation the MC of `fect` keeps a looser convergence tolerance for the final fit and does not fully converge on S3; fixing the penalty and refitting with the tolerance tightened to $10^{-7}$ moves the estimate (Section 6.5). Fourth, the inner solver of `tidysynth`, `ipop`, does not reach the exact optimum under the package's settings; tightening its tolerance `margin_ipop` moves the `tidysynth` numbers throughout the chapter, including the pre-period RMSPE, the gaps, the MSPE ratios, the donor weights, and the diagonal of $\boldsymbol\Omega$, while leaving the permutation ranking unchanged. Tightening the tolerance alone does not make the weights sparse either: the sparse exact solution requires re-solving at the stored $\boldsymbol\Omega$.

To summarize: SC estimation is a simplex-constrained minimum-distance weight optimization whose lifeline is the pre-period fit. For inference with a single treated unit, the in-space and in-time placebos locate the randomness in reassigning units and dates, while the conformal test and prediction intervals of Section 3.7 locate it in the time dimension. Augmented SC and synthetic DiD are two modern ways out when the fit is inadequate (Sections 3.3 and 3.4). Software defaults change the numbers and should be reported along with the results.

## 5 Anchoring Papers

The three anchoring papers correspond to three threads of this chapter: the founding application of SC, the re-estimation of the same data with synthetic DiD, and an IS study that uses a factor model under staggered adoption. Each is examined under five headings (paper, method, data, results, limitations), with a focus on how the authors defend their assumptions.

### 5.1 Abadie, Diamond and Hainmueller (2010)

::: {.case}
The paper: Abadie, A., Diamond, A., & Hainmueller, J. (2010). Synthetic control methods for comparative case studies: Estimating the effect of California's tobacco control program. Journal of the American Statistical Association, 105(490), 493-505. It is the founding application of the synthetic control method and the direct prototype of this chapter.

Method: in 1988 California passed Proposition 99, which raised the cigarette excise tax by 25 cents per pack and came with a tobacco control program. To estimate this policy's effect on per-capita cigarette sales, the difficulty is that there is only one treated unit, California, and no single state can serve as California's counterfactual on its own. SC's solution is to piece together a "synthetic California" from a convex combination of other states, making it reproduce California's per-capita cigarette sales and several predictors as closely as possible before 1988, so that after the policy the gap that opens between synthetic California and real California is the policy effect.

Data: a state-level panel of U.S. states from 1970 to 2000, with per-capita cigarette sales and a set of predictors. The donor pool consists of the remaining U.S. states; excluding the states that launched large tobacco control programs of their own or sharply raised cigarette taxes over the same period (such as Massachusetts, Arizona, Oregon, and Florida) leaves 38 states.

Results: synthetic California is a weighted combination of five states, Colorado, Connecticut, Montana, Nevada, and Utah (weights 0.164, 0.069, 0.199, 0.234, and 0.334), with zero weight on all other states. After the policy, real California's per-capita cigarette sales stay persistently below those of synthetic California: lower on average by almost 20 packs per person per year over 1989 to 2000, about 25%, and lower by about 26 packs per person per year by 2000. In the in-space placebo, which takes each donor state in turn as the treated unit, California's post/pre MSPE ratio is about 130, the largest among the 39 states, and the permutation p-value is 1/39 = 0.026.

Limitations: with only one treated unit, inference rests entirely on permutation logic, with no analytic standard error; the credibility of the conclusion hinges on the quality of the pre-period fit and the cleanliness of the donor pool. The authors discuss no interference: anti-smoking sentiment spilling over into other states, or tobacco companies shifting advertising budgets from other states to California, would depress the donors' sales and lead to an underestimate of the effect.
:::

Nearly all later SC studies follow this paper's template (convex-combination weights, pre-period fit, and in-space placebo inference), and the Aster case in Section 6 of this chapter is no exception.

### 5.2 Arkhangelsky, Athey, Hirshberg, Imbens and Wager (2021) Revisit Proposition 99

::: {.case}
The paper: Arkhangelsky, D., Athey, S., Hirshberg, D. A., Imbens, G. W., & Wager, S. (2021). Synthetic difference-in-differences. American Economic Review, 111(12), 4088-4118. The numbers in this section are taken from §2 and Table 1 of its arXiv version v4.

Method: on the same Proposition 99 data, the paper reports five estimators side by side: SDID, SC, DID, MC, and DIFP (the SC with an intercept of Doudchenko and Imbens and of Ferman and Pinto). All standard errors use the placebo method: draw one control state at random, treat it as the treated state, recompute the estimator, and use the dispersion of these estimates as the variance (Algorithm 4 of the paper).

Data: 39 states (California and 38 control states), 1970 to 2000, with 19 pre-treatment years and 12 post-treatment years.

Results: for the average effect on per-capita cigarette sales over the 12 post-treatment years, SDID gives -15.6 (8.4), SC -19.6 (9.9), DID -27.3 (17.7), MC -20.2 (11.5), and DIFP -11.1 (9.5). The authors argue that pre-treatment trends are clearly not parallel, so DID very likely overstates the effect; SC's -19.6 is comparable to ADH's "almost 20 packs on average"; SDID and DIFP are clearly smaller in absolute value than SC, which shows that time weights and unit fixed effects can materially change the conclusion, and when SC and SDID disagree the authors lean toward SDID. SDID is more flexible than DID yet has a smaller standard error, because the weighting improves the local fit; DID and SC both have highly influential states (such as New Hampshire), while SDID does not.

Limitations: with only one treated state the jackknife is not defined; the placebo variance relies on homoskedasticity across units, and if the treated state's noise differs from that of the control states, the variance cannot be learned from the control states. The placebo standard errors come from random draws and change with the seed when recomputed, so the values cited here are the paper's.
:::

This paper is the empirical counterpart of Proposition 2 in Section 3.3: changing only the weights and the fixed effects moves the estimate from -27.3 to -15.6. Aster's S1 offers a reminder in the other direction: with the truth known, SDID's 0.057 is farther from 0.095 than the hand-written SC's 0.076; the authors' preference for SDID comes from simulations and asymptotic theory, not from a guarantee for every dataset.

### 5.3 Song, Agarwal and Wen (2026, ISR)

::: {.case}
The paper: Song, F., Agarwal, A., & Wen, W. (2026). The impact of generative AI on collaborative open-source software development: Evidence from GitHub Copilot. Information Systems Research. Advance online publication. The numbers in this section are taken from the main text, Table 2, and Appendix Table G.6 of its arXiv version v4.

Method: GitHub projects start using Copilot in different months, which is staggered adoption. The treated group consists of projects in which developers use Copilot and whose development environment supports it; the control group consists of projects whose development environment does not support it and in which no one uses it. Estimation uses the generalized synthetic control (gsynth) of Xu (2017): a factor model with project and month fixed effects is estimated on the control projects only, and the loadings of each treated project are then estimated from its pre-treatment data to impute the counterfactual. The number of factors is chosen by cross-validation: each pre-treatment period of the treated projects is held out in turn, the loadings are estimated from the remaining periods and used to predict it, and the candidate value with the smallest prediction error is taken; the value selected is 0.

Data: a project-month panel from January 2021 to December 2022, with 7,637 projects (4,491 treated and 3,146 control) and 139,329 project-months.

Results: in Table 2, the ATT on the log number of merged pull requests is 0.059 (standard error 0.007), and the ATT on the log time to merge is 0.080 (0.021). Code contributions increase, and coordination also takes longer. The appendix runs TWFE on the full sample (9,244 projects) and obtains 0.029 (0.008) and 0.036 (0.018): the same direction, with smaller magnitudes.

Limitations: with 0 factors, gsynth reduces to imputing the counterfactual with two-way fixed effects alone, a version of the imputation of Section 4.5 of Chapter 19, and identification falls back on parallel trends. That cross-validation selects 0 factors shows only that factors estimated from the control projects do not improve the leave-one-out prediction of the treated projects' pre-treatment periods; it cannot prove that the two groups do not differ in their factor exposures, and the paper itself notes that factors able to absorb more confounding are not necessarily selected. Developers choose for themselves whether to use Copilot, and the authors run robustness checks using comparisons within the same development environment, propensity score matching combined with DiD, and instrumental variables.
:::

Taken together: ADH (2010) established SC's workflow and permutation inference; AAHIW show that the choice of weights and fixed effects can by itself move the estimate substantially; and the study by Song, Agarwal and Wen shows how factor models are used in IS research, and also that when the number of factors selected is 0 this machinery falls back to DiD, and the difference between its estimate and TWFE mixes sample differences with the weighting problem discussed in Chapter 19.

## 6 A Full Walkthrough on Aster Data

We now run the earlier tools from start to finish on Aster's four data blocks. The code below uses R 4.6.1, and every number cited in the text comes from the actual output of this chapter's code; only the main calls are shown here, with random seeds and data preparation omitted. SC uses `tidysynth`, and the other software is listed in the table of Section 4.3; `sc_weights`, `penalized_sc`, and `conformal_p` are functions written for this chapter, which compute the SC weights, the penalized SC weights, and the conformal p-value of Section 3.7, respectively.

### 6.1 The Data-Generating Process for S1

The design parameters are as follows: fix `set.seed(12021)`, with 20 donor metros plus 1 treated metro (Aster-Metro), 36 months, and the reform launching in month 25, the first month after the text's $T_0$, so the post-period is months 25 through 36, 12 months in all (the `T0` in the code below refers to this first treated month, 25). The data are generated by a factor model, two common factors plus city-specific loadings, with very small noise ($\sigma = 0.02$). The treated city's loadings are a convex combination of the loadings of four donors (numbers 3, 7, 12, and 18), with weights 0.35/0.30/0.20/0.15, which guarantees that its loadings lie inside the convex hull of the donors' loadings and that SC can reproduce its factor exposure.

```r
set.seed(12021)
n_donors <- 20; Tt <- 36; T0 <- 25; t_idx <- 1:Tt
f1 <- 2.0 * (1 - exp(-t_idx / 12)) + 0.05 * t_idx  # common factor 1: platform adoption trend (concave)
f2 <- sin(2 * pi * t_idx / 12) + 0.02 * t_idx       # common factor 2: seasonal term
tau_path <- ifelse(t_idx >= T0, 0.12 * (1 - exp(-(t_idx - T0 + 1) / 3)), 0)  # true effect: smooth ramp
```

Because the cities have different loadings and their trends are naturally not parallel, the parallel trends of DiD fails: the treated city's pre-period slope is 0.131 and the donor average is 0.117. The naive two-period, two-group DiD gives 0.244, about 2.6 times the true post-period average effect of 0.095, which confirms the number from Section 1 that did not add up.

### 6.2 SC Estimation for S1: Weights, Fit, and the Gap

Build the synthetic control with `tidysynth`, taking as predictors the pre-period mean plus the lags at months 6, 12, 18, and 24, with the weights optimized over months 1 through 24:

```r
library(tidysynth)
sc <- panel %>%
  synthetic_control(outcome = log_gmv, unit = metro, time = month,
                    i_unit = "Aster-Metro", i_time = 25, generate_placebos = TRUE) %>%
  generate_predictor(time_window = 1:24, mean_pre = mean(log_gmv)) %>%
  generate_predictor(time_window = 6,  gmv_m06 = log_gmv) %>%
  generate_predictor(time_window = 12, gmv_m12 = log_gmv) %>%
  generate_predictor(time_window = 18, gmv_m18 = log_gmv) %>%
  generate_predictor(time_window = 24, gmv_m24 = log_gmv) %>%
  generate_weights(optimization_window = 1:24) %>%
  generate_control()
```

The pre-period RMSPE is 0.017, comparable to the noise standard deviation of 0.02. The average post-period gap of 0.086 is SC's effect estimate, already quite close to the truth of 0.095 and far better than the naive DiD's 0.244; the gap in the last month is 0.103, against a truth of 0.118. The donors with the largest weights in the synthetic city are Donor-10 (0.240), Donor-02 (0.155), Donor-04 (0.108), Donor-08 (0.059), Donor-12 (0.039), and Donor-13 (0.037), yet these weights are not sparse: 19 donors have weights above 1%. This is not because the weighting problem has multiple solutions, but because the solver `ipop` does not reach the exact optimum (Section 4.3): the exact solution of the problem is unique, with only a few nonzero donors; the hand-written SC that uses only pre-period outcomes uses 9 donors, and the SC for S2 uses 3. The hand-written SC is sparse because it matches 24 noisy outcomes: 20 donors span a convex hull of at most 19 dimensions, the noise pushes the treated city's 24-dimensional outcome vector outside the hull, and the weights land on a boundary face of the hull (Section 4.1).

::: {.warning}
The weights SC estimates do not correspond to the four donors used to generate the data (numbers 3, 7, 12, and 18; see Section 6.1). The reason is that with only two factors plus a level, many different convex combinations can reproduce the treated city's factor exposure, so the weights that reproduce the factor exposure are not unique (whereas the solution of the weighting problem is). This is not a bug. What SC delivers is a good pre-period fit and a credible counterfactual, not an exact recovery of the weights, so do not read the weights as "the true control group."
:::

### 6.3 Placebo Inference for S1

Inference relies on placebos. The in-space placebo takes each donor in turn as the treated unit and computes the post/pre MSPE ratio. The treated city's ratio is as high as 32.777 (`tidysynth` counts month 25 as part of the pre-period; recomputing on the window of Section 4.2 while also switching to the hand-written SC gives 42.975, with the same ranking), first among all 21 units, with a permutation p-value of 0.048. Note that 0.048 is the strongest result 21 units can give ($1/21 = 0.048$). For inference with a single treated unit, all of the weight rests on this ranking.

The in-time placebo sets the reform date 8 months before the actual reform and re-estimates:

```r
sc_time <- panel %>%
  synthetic_control(outcome = log_gmv, unit = metro, time = month,
                    i_unit = "Aster-Metro", i_time = 25 - 8, generate_placebos = FALSE) %>%
  generate_predictor(time_window = 1:16, mean_pre = mean(log_gmv)) %>%
  generate_predictor(time_window = 6,  gmv_m06 = log_gmv) %>%
  generate_predictor(time_window = 12, gmv_m12 = log_gmv) %>%
  generate_weights(optimization_window = 1:16) %>%
  generate_control()
```

In the window after this fake, not-yet-reformed date and before the actual reform, the average gap is -0.019, but this number moves with the solver tolerance (Section 4.3), so it counts only as supporting evidence.

![Synthetic control for Aster's single-city fee reform. Top panel: the actual log-GMV trajectory of the treated city Aster-Metro (orange) and the synthetic city (blue dashed), with a vertical line marking the month the reform launched (the first month after $T_0$); before the reform the two lines nearly coincide (pre-period RMSPE 0.017), and after the reform a gap opens. Bottom panel: the gap plot. The treated city's gap (orange) stays near zero before the reform and opens to about 0.10 after it; the gray background lines are the in-space placebo gaps obtained by taking each donor in turn as the treated unit, with those whose pre-period RMSPE is more than 5 times the treated city's dropped, which leaves 3, and the treated city's average post-reform gap is larger than that of each of them.](assets/fig/fig_22_synth.svg)

The remaining S1 numbers in Sections 3.3 and 3.7 and the placebo intervals of Section 6.7 come from the calls below; `sc_estimate` is the SC in `synthdid`, the estimator of Proposition 2 (ii) in Section 3.3:

```r
library(synthdid); library(scpi)
# Ysd: donors in rows 1-20, Aster-Metro in row 21; Y0: months x donors
did  <- did_estimate(Ysd, 20, 24)                    # uniform weights: the 2x2 DiD and TWFE
sqrt(vcov(did, method = "placebo", replications = 500))
sdid <- synthdid_estimate(Ysd, 20, 24)
sqrt(vcov(sdid, method = "placebo", replications = 500))
sc   <- sc_estimate(Ysd, 20, 24)                     # no intercept, no time weights
sqrt(vcov(sc, method = "placebo", replications = 500))
w_h0 <- sc_weights(y1, Y0, diag(36))                 # SC fitted on all 36 months under H0
conformal_p(y1, drop(Y0 %*% w_h0), T0 = 24, scheme = "iid", n_perm = 9999)
conformal_p(y1, drop(Y0 %*% w_h0), T0 = 24, scheme = "moving_block")
sdat <- scdata(df, id.var = "metro", time.var = "month", outcome.var = "log_gmv",
               period.pre = 1:24, period.post = 25:36, unit.tr = "Aster-Metro", unit.co = donors)
scpi(sdat, w.constr = list(name = "simplex"), sims = 200, e.method = "gaussian")
```

### 6.4 S2: A Treated City Outside the Convex Hull

The setup of S2 is the same as that of S1 (two factors, 20 donors, 24 pre-period months); the donors are a new set drawn from the same distribution, keeping only the names, and the treated city becomes Aster-Capital: its two loadings sit in the middle of the donors' range (0.9 and 0), but its intrinsic level is 0.232 above that of the highest donor, and in the noise-free generating parameters its $\ell_1$ distance to the donors' convex hull is 0.298. The main calls are as follows:

```r
library(augsynth); library(synthdid)
# y1: Aster-Capital; Y0: months x donors; months 1-24 are the pre-period
w_sc <- sc_weights(y1[1:24], Y0[1:24, ], diag(24))           # SC on outcomes, Omega = I
pan  <- subset(S2$panel, select = -unit)                      # augsynth clashes with an integer column named unit
asc  <- augsynth(log_gmv ~ treat_it, metro, month, pan, progfunc = "ridge", scm = TRUE)
Ysd  <- rbind(Y[-1, ], Y[1, ])                                # synthdid: donors first, treated last
sdid <- synthdid_estimate(Ysd, 20, 24)
se   <- sqrt(vcov(sdid, method = "placebo", replications = 500))
pen  <- penalized_sc(y1[1:24], Y0[1:24, ], pen = 10^-2.75)    # chosen by leave-one-out CV over 21 grid values
```

The estimates and errors are in the table and figure of Section 3.4. SC puts its weight on Donor-12, Donor-17, and Donor-13 (0.433, 0.402, and 0.165), the three donors with the highest intrinsic levels. The penalty of ridge ASC is selected at 0.752, and the correction term matters; its conformal p-value is 0.204: extrapolation removes most of the bias but does not buy enough precision to reject a zero effect. The placebo standard error of SDID is 0.037, and its interval [-0.042, 0.104] covers the truth of 0.079. The penalty of penalized SC is chosen by leave-one-out cross-validation as in §4.1 of Abadie and L'Hour (2021), with each donor taking a turn as the treated unit; if one instead uses the rule of §4.2 of the paper, holding out a stretch at the end of the treated city's own pre-period for validation, the penalty is 0.178, only 2 donors remain, and the estimate is -0.160. Neither rule fixes the sign: the penalty governs uniqueness and sparsity, not extrapolation.

### 6.5 S3: Staggered Adoption

The design parameters of S3 are as follows: 80 metros and 48 months; 20 cities that cannot be selected because of existing contracts are chosen at random, the other cities are ranked by their cyclical loading plus noise with standard deviation 0.05, and the top 15 are rolled out by rank in three waves of 5. Within this group the loadings are similar, the differences are swamped by the ranking noise, and the waves are close to random. The idiosyncratic shocks are an independent part (standard deviation 0.02) plus an AR(1) with coefficient 0.8 and stationary standard deviation 0.04. The main calls are as follows:

```r
library(fixest); library(did); library(fect); library(augsynth); library(synthdid)
twfe <- feols(log_gmv ~ treat_it | unit + month, data = p, cluster = ~unit)
p$gname <- ifelse(is.finite(p$cohort), p$cohort, 0)          # did codes never-treated as 0
gt   <- att_gt("log_gmv", "month", "unit", gname = "gname", data = p,
               control_group = "nevertreated", est_method = "reg", biters = 999)
cs   <- aggte(gt, type = "simple")
cs_ny <- aggte(att_gt("log_gmv", "month", "unit", gname = "gname", data = p,
                      control_group = "notyettreated", est_method = "reg", biters = 999), type = "simple")
ife  <- fect(log_gmv ~ treat_it, data = p, index = c("unit", "month"), method = "ife",
             force = "two-way", CV = TRUE, r = c(0, 4), se = TRUE, nboots = 500)
mc0  <- fect(log_gmv ~ treat_it, data = p, index = c("unit", "month"), method = "mc",
             force = "two-way", CV = TRUE, se = TRUE, nboots = 500)        # fect default
mc   <- fect(log_gmv ~ treat_it, data = p, index = c("unit", "month"), method = "mc",
             force = "two-way", CV = FALSE, lambda = mc0$lambda.cv, tol = 1e-7,
             max.iteration = 50000, se = TRUE, nboots = 500)               # refit at the CV penalty
ms   <- multisynth(log_gmv ~ treat_it, unit, month, p)
summary(ms)                                                                 # bootstrap intervals
# SDID by wave: each wave against the never-treated metros, weighted by its treated metro-months
never <- which(is.infinite(cohort))
est   <- sapply(c(29, 33, 37), function(g)
  as.numeric(synthdid_estimate(Y[c(never, which(cohort == g)), ], length(never), g - 1)))
cells <- sapply(c(29, 33, 37), function(g) sum(cohort == g) * (48 - g + 1))
sum(cells * est) / sum(cells)                                              # metro bootstrap, B = 500
library(scpi)
sm   <- scdataMulti(p, id.var = "metro", time.var = "month", outcome.var = "log_gmv",
                    treatment.var = "treat_it", effect = "time")           # never-treated donors
scpi(sm, w.constr = list(name = "simplex"), sims = 200, e.method = "gaussian") # one interval per event month
```

The table below sets this S3 dataset side by side with the 24 draws (the latter include the former; the `fect` bootstrap uses 100 replications in the 24 draws and 500 for this dataset), with a truth of 0.06 in all cases:

| Estimator | This estimate | 95% interval | Mean estimate over 24 draws | Times the interval covers the truth |
|---|---|---|---|---|
| TWFE | 0.099 | [0.071, 0.127] | 0.102 | 0 |
| CS (never-selected as controls) | 0.139 | [0.105, 0.173] | 0.149 | 0 |
| CS (not-yet-selected as controls) | 0.139 | [0.105, 0.173] | Not redrawn | Not redrawn |
| IFE | 0.051 | [0.028, 0.073] | 0.063 | 24 |
| MC (`fect` default) | 0.074 | [0.053, 0.094] | 0.088 | 6 |
| MC (refitted with tightened tolerance) | 0.067 | [0.047, 0.087] | 0.080 | 14 |
| multisynth | 0.053 | [-0.018, 0.135] | 0.062 | 24 |
| SDID by wave | 0.054 | [0.031, 0.076] | Not redrawn | Not redrawn |
| `scpi` | 0.046 | 2 of 12 monthly prediction intervals cover 0.06 | Not redrawn | Not redrawn |

Reading this table requires keeping three things apart. First, the bias of the DiD family comes from the design: in all 24 draws TWFE and CS overestimate, and their intervals never cover 0.06, because the platform picks cities by cyclical exposure. Second, this draw is unusually favorable to MC: both MC estimates are the lowest of the 24 draws; the default procedure removes 0.65 of the TWFE bias this time, against a mean of the per-draw shares of only 0.34 across the 24 draws, and the corresponding numbers for the refit are 0.82 and 0.55. Looking at this one dataset alone, one would conclude that MC is as reliable as IFE. Third, `scpi` fits each city separately ($\pi_{\mathrm{pool}} = 0$), so the error of the average path, $q^{\mathrm{pool}}$ (Theorem 7), becomes a common bias, and the intervals are too narrow to allow for it. For IFE and multisynth, and the figure, see Sections 3.5 and 3.6.

### 6.6 S4: Pilot Design

The design parameters of S4 are as follows: 40 metros and 24 months of history; the common demand shock is an AR(1) with coefficient 0.7 and innovation standard deviation 0.05, and the idiosyncratic noise has standard deviation 0.02; the pilot lasts 6 months, and once piloted, every city has an effect of 0.04; evaluation uses 500 replicated pilot periods that redraw only the shocks of the pilot period. Design problem (6) of the paper searches over all pilot sets; here we solve two simplex problems for each candidate set and run a local search across sets:

```r
# X: 24 history months x 40 metros; xb: the 40-metro average history
fit_set <- function(S) {                                     # S: candidate pilot metros
  C <- setdiff(1:40, S)
  w <- sc_weights(xb, X[, S, drop = FALSE], diag(24))        # pilot side reproduces xb
  v <- sc_weights(xb, X[, C, drop = FALSE], diag(24))        # control side, from the other metros
  sum((xb - X[, S, drop = FALSE] %*% w)^2) + sum((xb - X[, C, drop = FALSE] %*% v)^2)
}
# local search over S (swap, add or drop one metro) from a greedy start and 100 random starts
```

Of the 101 starting points, 6 reach the best set found, which is Metro-01, 09, 14, 16, and 27 of Section 3.8, and the control side uses 25 cities. Performance over the 500 replicated pilot periods is as follows:

| Design and estimator | RMSE | Bias | Standard deviation |
|---|---|---|---|
| Synthetic design | 0.0087 | 0.0062 | 0.0061 |
| Random 5, DiD | 0.077 | -0.004 | 0.077 |
| Random 5, difference in means | 0.212 | 0.019 | 0.211 |

The random designs are nearly unbiased and lose on variance; the synthetic design's bias and standard deviation are of the same order, and both are small. In replicated pilot period 1, the synthetic design estimates 0.053. The figure is in Section 3.8.

### 6.7 The Grand Reconciliation

| Block | Estimator | Estimate | Interval or test | Truth |
|---|---|---|---|---|
| S1 | Naive DiD (that is, SDID with uniform weights, and TWFE) | 0.244 | [-0.959, 1.447] | 0.095 |
| S1 | SC, `tidysynth` | 0.086 | Permutation p-value 0.048 | 0.095 |
| S1 | SC, $\boldsymbol\Omega = \mathbf{I}$ | 0.076 | Conformal p-value 0.082 or 0.222 | 0.095 |
| S1 | SC, `synthdid` | 0.078 | [-0.326, 0.483] | 0.095 |
| S1 | SDID | 0.057 | [-0.012, 0.125] | 0.095 |
| S2 | SC | -0.119 | | 0.079 |
| S2 | Ridge ASC | 0.044 | | 0.079 |
| S2 | SDID | 0.031 | [-0.042, 0.104] | 0.079 |
| S3 | TWFE | 0.099 | [0.071, 0.127] | 0.06 |
| S3 | IFE | 0.051 | [0.028, 0.073] | 0.06 |
| S4 | Synthetic design (replicated pilot period 1) | 0.053 | RMSE 0.0087 | 0.04 |

The intervals for S1 and S2 are all built from placebo standard errors: the intervals of DiD and of the SC in `synthdid` are so wide as to be uninformative, possibly because these estimators fit the placebo units poorly; the weighting of SDID probably improves on this, consistent with the explanation of AAHIW in Section 5.2.

The reconciliation of this section can be summarized as follows. On S1, the naive DiD of 0.244 is inflated to 2.6 times the truth because trends are not parallel; synthetic control uses the tight pre-period fit of 0.017 (the solver `ipop` does not reach the exact optimum; see Section 4.3 for a note on the tolerance) to pull the estimate back to 0.086, close to the truth of 0.095; the in-space placebo ranking supports the conclusion that this gap comes from the reform, and the in-time placebo counts only as supporting evidence. On S2, SC outside the convex hull gets even the sign wrong, and extrapolation and fixed effects each recover part of the effect. On S3, selection enters through factor exposure, the DiD family fails in every draw, IFE is stable, and MC's good showing belongs to this one draw only. On S4, the pilot is picked in advance from history, and random selection followed by DiD has 8.8 times its RMSE.

## 7 Failure Modes and Robustness

In a simulation the identifying assumptions are built in, but in real research they can fail at any time. This section lays out the most common ways SC fails and the actionable responses.

The first kind of failure is a poor pre-period fit. If the synthetic city cannot track the treated city before the reform, the loading-matching argument collapses and the post-period gap cannot be interpreted. Check the pre-period RMSPE, and if the fit is poor, consider augmented SC, which patches the residual imbalance with an outcome model. The second kind is a convex-hull violation: the treated unit's characteristics are too extreme and lie outside the donors' convex hull, no convex combination can reproduce it, and under the simplex constraint SC can give only a suboptimal fit that clings to the edge of the hull. On S2 such a fit gets even the sign wrong. Report ASC, SDID, and SC side by side; the gaps among the three are themselves a diagnostic. Interpolation bias can arise inside the convex hull too (Section 3.4), and the response is to keep only similar donors or to use penalized SC. The third kind is contamination of the donor pool: a donor that is in fact treated as well, or touched by spillovers from the intervention on the treated unit, or hit by a large shock of its own that the treated unit does not experience, is no longer clean. In a platform setting this is a particular concern: lower fees may draw sellers in from neighboring cities and prompt buyers to order across cities (the interference of Section 3.1). The response is to exclude the obviously entangled donors from the pool, or to redefine the market boundary.

SC has two subtler problems. One is too few or too many donors: with too few you cannot assemble a good convex combination, and with too many you tend to overfit noise in the pre-period, piecing together a synthetic city that resembles the treated city only in the pre-period and falls apart in the post-period. The other is overfitting the pre-period path itself: a long and flexible fit may be matching idiosyncratic noise rather than the factor structure. The in-time placebo is exactly the diagnostic against this: if a gap already appears between the fake date and the actual reform, overfitting is the most likely cause. Finally, inference with a single treated unit has no analytic standard error. The conformal test moves the randomness to the time dimension, at the cost of taking a stand on the time structure of the shocks: whether the p-value on S1 is 0.082 or 0.222 depends on exactly this, and both should be reported. With too few donors, permutation inference runs out of power; this is an inborn limitation of SC, to be admitted honestly in a paper rather than covered up with a precise-looking asterisk.

With many treated units, failure takes other forms. A factor model stands or falls with the number of factors: too many overfit, too few leave non-parallel trends behind, and selecting 0 falls back to imputation with two-way fixed effects (Section 5.3); report the selected $r$ and re-estimate at neighboring values of $r$. MC does not select $r$ but shrinks the treated block toward the two-way fixed effects fit, and on S3 it overestimates in every draw. Persistent shocks do not prevent IFE from identifying the effect, but they strip the theorems for MC and staggered SC of their guarantees. Be explicit about the selection mechanism: if treated units are picked by factor exposure, a factor model can handle it and DiD cannot; if they are picked by idiosyncratic shocks, neither can rescue the estimate.

With staggered adoption, an analysis with multisynth should report $\pi_{\mathrm{pool}}$ and include a sensitivity analysis over its values. An ex ante synthetic design can fail too: with a short fitting period and many candidate cities, the chosen set may fit the history only by chance; the bound of Theorem 9 grows with $N$ and shrinks as $T_E$ lengthens. Abadie and Zhao leave a blank period after the fitting period and use its residuals for a permutation test (their Theorem 2); S4 in this chapter is fitted on all 24 months of history and does not take this step.

Taken together, these failure modes teach one lesson: the credibility of SC and its extensions lies not in the estimators themselves but in the design and the institutions, and their success or failure hinges on whether the pre-period fit is really good, whether the donors are really clean, and how treatment was assigned. The pre-period RMSPE, the in-space and in-time placebos, the side-by-side comparison of several estimators, and checks on redrawn data are all evidence offered on these questions, and none of them can substitute for a substantive judgment about the institutional background.

## 8 Further Reading

::: {.readings}
Required reading, in suggested order:

- Abadie, Diamond and Hainmueller (2010, Journal of the American Statistical Association). The founding application of synthetic control, California's Proposition 99, and the direct prototype of this chapter; focus on the in-space placebo and the weight optimization.
- Abadie (2021, Journal of Economic Literature). The authoritative survey of the SC method, which maps out the convex hull, loading matching, inference, and extensions; the first read for an overview.
- Arkhangelsky, Athey, Hirshberg, Imbens and Wager (2021, American Economic Review). Synthetic DiD, which blends SC and DiD and provides formal inference; the re-estimation of Proposition 99 in §2 and Theorem 1 are the sources of Sections 3.3 and 5.2 of this chapter.
- Doudchenko and Imbens (2016, NBER Working Paper). Writes SC, DiD, and regression as different constraints on one weighting problem; the source of Section 3.2 of this chapter.
- Arkhangelsky and Imbens (2024, Econometrics Journal). A survey of causal panel models that places SC, factor models, and DiD in a single framework.

Further reading:

- Abadie and Gardeazabal (2003, American Economic Review). The methodological origin of synthetic control: the effect of terrorism on the economy of the Basque Country.
- Ben-Michael, Feller and Rothstein (2021, Journal of the American Statistical Association). Augmented SC: bias correction when the pre-period fit is imperfect.
- Ferman and Pinto (2021, Quantitative Economics). Why SC is generally biased when the fit is imperfect, and when it does worse than DiD.
- Abadie and L'Hour (2021, Journal of the American Statistical Association). Uniqueness and sparsity of penalized SC.
- Kaul, Klößner, Pfeifer and Schieler (2022, Journal of Business & Economic Statistics). Why covariates play no role when all pre-period outcomes are used as predictors.
- Bai (2009, Econometrica), Xu (2017, Political Analysis), and Liu, Wang and Xu (2024, American Journal of Political Science). The asymptotic theory of interactive fixed effects, gsynth, and the counterfactual estimators in `fect`.
- Athey, Bayati, Doudchenko, Imbens and Khosravi (2021, Journal of the American Statistical Association). A unified representation of matrix completion and its error bounds.
- Ben-Michael, Feller and Rothstein (2022, Journal of the Royal Statistical Society Series B). Partially pooled SC for staggered adoption.
- Chernozhukov, Wüthrich and Zhu (2021, Journal of the American Statistical Association) and Cattaneo, Feng and Titiunik (2021, Journal of the American Statistical Association). Conformal inference and prediction intervals; Cattaneo, Feng, Palomba and Titiunik (2025, Journal of Statistical Software) documents `scpi`, whose `scdataMulti` handles staggered adoption (see Section 6.5).
- Abadie and Zhao (published online in 2026 in the Review of Economics and Statistics). Using SC to choose pilot units ex ante.
- Song, Agarwal and Wen (published online in 2026 in Information Systems Research) and Lambrecht, Tucker and Zhang (2024, Journal of Marketing Research). Two applications in IS and marketing: the former uses gsynth for staggered adoption, and the latter uses SDID in a field experiment with TV advertising rolled out region by region.
:::

::: {.apa-refs}
- Abadie, A. (2021). Using synthetic controls: Feasibility, data requirements, and methodological aspects. *Journal of Economic Literature, 59*(2), 391-425. https://doi.org/10.1257/jel.20191450
- Abadie, A., Diamond, A., & Hainmueller, J. (2010). Synthetic control methods for comparative case studies: Estimating the effect of California's tobacco control program. *Journal of the American Statistical Association, 105*(490), 493-505. https://doi.org/10.1198/jasa.2009.ap08746
- Abadie, A., & Gardeazabal, J. (2003). The economic costs of conflict: A case study of the Basque Country. *American Economic Review, 93*(1), 113-132. https://doi.org/10.1257/000282803321455188
- Abadie, A., & L'Hour, J. (2021). A penalized synthetic control estimator for disaggregated data. *Journal of the American Statistical Association, 116*(536), 1817-1834. https://doi.org/10.1080/01621459.2021.1971535
- Abadie, A., & Zhao, J. (2026). Synthetic controls for experimental design. *The Review of Economics and Statistics*. Advance online publication. https://doi.org/10.1162/rest.a.1837
- Arkhangelsky, D., Athey, S., Hirshberg, D. A., Imbens, G. W., & Wager, S. (2021). Synthetic difference-in-differences. *American Economic Review, 111*(12), 4088-4118. https://doi.org/10.1257/aer.20190159
- Arkhangelsky, D., & Imbens, G. (2024). Causal models for longitudinal and panel data: A survey. *The Econometrics Journal, 27*(3), C1-C61. https://doi.org/10.1093/ectj/utae014
- Athey, S., Bayati, M., Doudchenko, N., Imbens, G., & Khosravi, K. (2021). Matrix completion methods for causal panel data models. *Journal of the American Statistical Association, 116*(536), 1716-1730. https://doi.org/10.1080/01621459.2021.1891924
- Bai, J. (2009). Panel data models with interactive fixed effects. *Econometrica, 77*(4), 1229-1279. https://doi.org/10.3982/ECTA6135
- Ben-Michael, E., Feller, A., & Rothstein, J. (2021). The augmented synthetic control method. *Journal of the American Statistical Association, 116*(536), 1789-1803. https://doi.org/10.1080/01621459.2021.1929245
- Ben-Michael, E., Feller, A., & Rothstein, J. (2022). Synthetic controls with staggered adoption. *Journal of the Royal Statistical Society Series B: Statistical Methodology, 84*(2), 351-381. https://doi.org/10.1111/rssb.12448
- Cattaneo, M. D., Feng, Y., Palomba, F., & Titiunik, R. (2025). scpi: Uncertainty quantification for synthetic control methods. *Journal of Statistical Software, 113*(1), 1-38. https://doi.org/10.18637/jss.v113.i01
- Cattaneo, M. D., Feng, Y., & Titiunik, R. (2021). Prediction intervals for synthetic control methods. *Journal of the American Statistical Association, 116*(536), 1865-1880. https://doi.org/10.1080/01621459.2021.1979561
- Chernozhukov, V., Wüthrich, K., & Zhu, Y. (2021). An exact and robust conformal inference method for counterfactual and synthetic controls. *Journal of the American Statistical Association, 116*(536), 1849-1864. https://doi.org/10.1080/01621459.2021.1920957
- Doudchenko, N., & Imbens, G. W. (2016). *Balancing, regression, difference-in-differences and synthetic control methods: A synthesis* (NBER Working Paper No. 22791). National Bureau of Economic Research. https://doi.org/10.3386/w22791
- Ferman, B., & Pinto, C. (2021). Synthetic controls with imperfect pretreatment fit. *Quantitative Economics, 12*(4), 1197-1221. https://doi.org/10.3982/QE1596
- Kaul, A., Klößner, S., Pfeifer, G., & Schieler, M. (2022). Standard synthetic control methods: The case of using all preintervention outcomes together with covariates. *Journal of Business & Economic Statistics, 40*(3), 1362-1376. https://doi.org/10.1080/07350015.2021.1930012
- Lambrecht, A., Tucker, C., & Zhang, X. (2024). TV advertising and online sales: A case study of intertemporal substitution effects for an online travel platform. *Journal of Marketing Research, 61*(2), 248-270. https://doi.org/10.1177/00222437231180171
- Liu, L., Wang, Y., & Xu, Y. (2024). A practical guide to counterfactual estimators for causal inference with time-series cross-sectional data. *American Journal of Political Science, 68*(1), 160-176. https://doi.org/10.1111/ajps.12723
- Song, F., Agarwal, A., & Wen, W. (2026). The impact of generative AI on collaborative open-source software development: Evidence from GitHub Copilot. *Information Systems Research*. Advance online publication. https://doi.org/10.1287/isre.2024.1361
- Xu, Y. (2017). Generalized synthetic control method: Causal inference with interactive fixed effects models. *Political Analysis, 25*(1), 57-76. https://doi.org/10.1017/pan.2016.2
:::
