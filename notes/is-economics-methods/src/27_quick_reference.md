---
title: "Quick Reference & Tools"
subtitle: "Formula Cards, a Method Decision Tree, and a Cross-Software Command Map"
seriesline: "Foundations of Information Systems Economics · Chapter 27"
series: "Foundations of Information Systems Economics: Empirical Methods for Information Systems Economics"
sidebartitle: "Foundations of Information Systems Economics"
sidebarsub: "Theory · Methods · Applications"
chapterlabel: "Chapter 27 · Quick Reference & Tools"
---

## Introduction

The moment you actually need a reference sheet is usually not before an exam, but after the regression has already run. Should the standard errors cluster by user, by store, or by market? Can a staggered event study still be run as plain TWFE? When the instrument is weak, which interval is still trustworthy? The numbers on the screen tend to arrive faster than these questions, and they are all too willing to tempt you into writing up a conclusion on the spot. This chapter sits on your desk for one reason: to buy that moment one extra pause.

We organize everything by research question first, not by software menu. The formula cards help you check the estimand and the notation; the decision tree starts from your goal, your data structure, and where the identifying variation comes from, and points you at the relevant chapter; the R, Python, and Stata command tables answer "how do I run this"; and the pitfalls checklist presses the harder question, "why might that implementation still be wrong". See a marketplace experiment, for instance, and the path does not stop at a difference in means: it reminds you of interference, switchback designs, and the rollout estimand. See AI adoption, and it will not let you treat access, use, and workflow redesign as one and the same treatment.

The danger of a reference sheet is exactly its convenience. A formula has no room for the full institutional background, a command will not check the exclusion restriction for you, and a decision tree cannot decide on your behalf whether parallel trends is credible. So what this chapter offers is an entrance and a set of signposts for the way back, not a verdict on which method to use: once you have a candidate, return to the corresponding chapter and check the economic model, the identifying assumptions, the estimation details, and the failure modes.

Notation is uniform throughout: individuals are $i$, products $j$, markets or time $t$, firms $f$, cohorts $g$; potential outcomes are written $Y_i(d)$, treatment is $D_i$, the instrument is $Z_i$, covariates are $X_i$, and the propensity score is $e(x)$. If this chapter does its job, it will not help you pick a fashionable method faster, it will help you discover faster that the question in front of you has not yet been defined clearly.

## 1 Formula Cards

Each card gives one core method: the name, the key formula, and a one-line gloss. Formulas are in English and LaTeX throughout, and the cards are grouped by the book's five Parts. To understand the conditions and limits under which a formula holds, return to the chapter noted in parentheses.

### Part I Microeconomic Theory Foundations

Part I is the pure-theory part: the main text demonstrates by hand in base R, with no dedicated estimation package. The cards below give the primitives that later identification and structural estimation borrow again and again.

::: {.theorem data-label="Consumer duality (Roy / Shephard / the four duality identities) · Chapter 1"}
$$x_\ell(p, w) = -\frac{\partial v(p, w)/\partial p_\ell}{\partial v(p, w)/\partial w}, \qquad \frac{\partial e(p, \bar u)}{\partial p_\ell} = h_\ell(p, \bar u)$$

$$e(p, v(p, w)) = w, \qquad h_\ell(p, \bar u) = x_\ell(p, e(p, \bar u))$$

In words: the indirect utility $v$ and the expenditure function $e$ are inverses of each other, differentiating with respect to price reads Walrasian demand $x$ and Hicksian demand $h$ straight off the value function, and the four objects are pinned to one another pairwise.
:::

::: {.theorem data-label="Slutsky equation · Chapter 1"}
$$\underbrace{\frac{\partial x_\ell(p, w)}{\partial p_k}}_{\text{total}} = \underbrace{\frac{\partial h_\ell(p, \bar u)}{\partial p_k}}_{\text{substitution}}\; \underbrace{-\; x_k(p, w)\,\frac{\partial x_\ell(p, w)}{\partial w}}_{\text{income}}$$

In words: the observable price slope of demand splits into a compensated substitution effect and a purchasing-power income effect, and the substitution piece corresponds to a negative semidefinite Slutsky matrix.
:::

::: {.theorem data-label="Welfare measures CV / EV · Chapter 1"}
$$\mathrm{CV} = e(p^1, u^0) - e(p^1, u^1) = e(p^1, u^0) - w, \qquad \mathrm{EV} = e(p^0, u^0) - e(p^0, u^1) = w - e(p^0, u^1)$$

In words: the money-metric welfare of a price change is the difference between two evaluations of the expenditure function, CV is priced at the new prices and EV at the old, and the usual ordering is $\mathrm{EV} \leq \Delta\mathrm{CS} \leq \mathrm{CV}$.
:::

::: {.theorem data-label="Arrow-Pratt risk aversion · Chapter 1"}
$$A(w) = -\frac{u''(w)}{u'(w)}, \qquad \pi \approx \tfrac{1}{2}A(w)\,\mathrm{Var}(w)$$

In words: the curvature of the Bernoulli utility measures risk aversion, and for a small risk the risk premium is roughly half the coefficient of absolute risk aversion times the variance.
:::

::: {.theorem data-label="Nash equilibrium · Chapter 2"}
$$u_i(s_i^*, s_{-i}^*) \geq u_i(s_i, s_{-i}^*) \quad \text{for all } s_i,\ \text{all } i$$

In words: every player's strategy is a best response to the others', so no one has an incentive to deviate unilaterally; in a finite game, once mixed strategies are allowed, existence is guaranteed.
:::

::: {.theorem data-label="Repeated-game cooperation threshold (grim trigger) · Chapter 2"}
$$\frac{R}{1-\delta} \geq T + \frac{\delta P}{1-\delta} \quad\Longleftrightarrow\quad \delta \geq \frac{T - R}{T - P}$$

In words: the discounted value of cooperating forever is at least that of defecting once and then being punished forever, so the cooperation threshold is set by the ratio of the temptation to the severity of punishment; the folk theorem says that when $\delta$ is large enough almost any feasible, individually rational payoff can be supported.
:::

::: {.theorem data-label="VCG mechanism · Chapter 3"}
$$t_i = \max_{a}\sum_{j \neq i} \hat v_j(a) \;-\; \sum_{j \neq i} \hat v_j(a^*)$$

In words: each agent pays the externality it imposes on everyone else, so private self-interest and total social welfare differ only by a constant that is independent of the agent's own report, and truth-telling becomes a dominant strategy (the revelation principle guarantees that we need only study direct, truthful mechanisms of this kind).
:::

::: {.theorem data-label="Revenue equivalence and the Myerson reserve · Chapter 3"}
$$b(v) = \frac{n-1}{n}\, v, \qquad r^\ast:\ r = \frac{1 - F(r)}{f(r)}$$

In words: any auction with risk-neutral bidders, independent private values, that awards the object to the highest valuation and leaves the lowest type zero rent earns the same expected revenue; the optimal reserve is where the virtual value hits zero, which is $r^\ast = 1/2$ under the uniform distribution.
:::

::: {.theorem data-label="Lerner index / monopoly markup · Chapter 4"}
$$\underbrace{p(q) + q\,\frac{dp}{dq}}_{MR} = \underbrace{c'(q)}_{MC}, \qquad \frac{p - c'}{p} = \frac{1}{|\varepsilon|}$$

In words: the monopolist sets marginal revenue equal to marginal cost, the markup equals the reciprocal of the price elasticity of demand, and the less elastic demand is the greater the market power; this inverse-elasticity relation is the starting point for recovering markups from demand in Chapter 21.
:::

::: {.theorem data-label="Cournot and Bertrand · Chapter 4"}
$$\text{Cournot}:\ p > c \qquad\qquad \text{Bertrand}:\ p = c = mc$$

In words: the dimension of competition determines market power. Cournot competition in quantities leaves a markup, while Bertrand competition in prices drives price down to marginal cost when products are homogeneous and costs identical (the Bertrand paradox); capacity constraints, product differentiation, or repeated play can push price back up.
:::

### Part II Econometric Core

Part II is the machinery of estimation and inference: every downstream reduced-form and structural method draws on the sandwich variance, clustering, GMM, and nonparametric smoothing collected here.

::: {.theorem data-label="Extremum estimator (MLE / NLS / GMM unified framework) · Chapter 5"}
$$\hat\theta = \arg\max_{\theta \in \Theta} Q_n(\theta), \qquad \hat\theta_{MLE} = \arg\max_\theta \frac{1}{n} \sum_{i=1}^n \ln f(z_i; \theta)$$

$$\hat\theta_{GMM} = \arg\min_\theta\ \bar g_n(\theta)'\, W\, \bar g_n(\theta), \qquad \bar g_n(\theta) = \frac{1}{n}\sum_{i=1}^n g(w_i, \theta)$$

In words: MLE, NLS, and GMM are all special cases of extremizing one objective function, identification requires $Q_0(\theta)$ to have a unique peak at the truth, and how steep that peak is measures the strength of identification.
:::

::: {.theorem data-label="Sandwich variance and the information equality · Chapter 5"}
$$\sqrt{n}(\hat\theta - \theta_0) \xrightarrow{d} \mathcal{N}\big(0,\ A^{-1} B A^{-1}\big), \quad A = -\,\mathbb{E}[\nabla^2 q], \ B = \mathbb{E}[\nabla q\, \nabla q']$$

In words: the variance of any extremum estimator is a sandwich, the bread is curvature and the filling is gradient noise; when the specification is exactly right $A = B = I(\theta_0)$ and it collapses to the Cramér-Rao bound $I^{-1}$. Report the sandwich by default, not $I^{-1}$.
:::

::: {.theorem data-label="GMM efficient weighting and the Hansen J · Chapter 5"}
$$\mathrm{Var}_{\text{eff}} = (D'S^{-1}D)^{-1}, \qquad n\, \bar g_n(\hat\theta)'\, \hat S^{-1}\, \bar g_n(\hat\theta) \xrightarrow{d} \chi^2_{q-k}$$

In words: efficient GMM takes the weight $W = S^{-1}$ to reach minimum variance ($D$ is the Jacobian of the moment conditions, $S$ their covariance), and under overidentification the objective value gives the Hansen J test, whose degrees of freedom are the number of moments minus the number of parameters.
:::

::: {.theorem data-label="Delta method · Chapter 6"}
$$\sqrt{n}\big(g(\hat\theta) - g(\theta_0)\big) \xrightarrow{d} \mathcal{N}\big(0,\ \nabla g(\theta_0)'\, \Sigma\, \nabla g(\theta_0)\big)$$

In words: a first-order Taylor expansion of a nonlinear transform at the truth gives standard errors for derived quantities like ratios, elasticities, and WTP; when the denominator is near zero (the Fieller problem) the approximation breaks down, and you should switch to percentile-t or Fieller intervals.
:::

::: {.theorem data-label="Cluster-robust variance and the Moulton factor · Chapter 6"}
$$\hat V^{CR} = (X'X)^{-1} \left( \sum_{g=1}^{G} X_g'\, \hat u_g\, \hat u_g'\, X_g \right) (X'X)^{-1}, \qquad \sqrt{1 + (\bar n - 1)\,\rho}$$

In words: adding the within-cluster cross-covariances back into the filling of the sandwich corrects the understatement caused by the iid assumption; the asymptotics run on the number of clusters $G$ rather than the number of observations, and the Moulton factor is the multiple by which within-cluster correlation inflates the standard error.
:::

::: {.theorem data-label="HAC (Newey-West) · Chapter 6"}
$$\hat\Omega^{HAC} = \sum_{t} \hat u_t^2\, x_t x_t' + \sum_{\ell=1}^{L} w_\ell \sum_{t=\ell+1}^{T} \hat u_t \hat u_{t-\ell}\big(x_t x_{t-\ell}' + x_{t-\ell} x_t'\big),\quad w_\ell = 1 - \frac{\ell}{L+1}$$

In words: the time-series version of clustering, using Bartlett triangular weights to add several lagged autocovariances into the variance estimate, handling serial correlation and heteroskedasticity.
:::

::: {.theorem data-label="The testing trinity Wald / LR / LM · Chapter 6"}
$$W = \frac{(\hat\theta - \theta_0)^2}{\widehat{\mathrm{Var}}(\hat\theta)}, \qquad LR = 2(\ell_u - \ell_r), \qquad LM = s(\hat\theta_r)'\, I(\hat\theta_r)^{-1}\, s(\hat\theta_r)$$

In words: the same restriction tested from three angles, distance, height, and slope, asymptotically equivalent but diverging in finite samples; Wald is not invariant to reparameterization, so when in doubt prefer LR.
:::

::: {.theorem data-label="Within (FE) and Nickell bias · Chapter 7"}
$$y_{it} - \bar y_i = (x_{it} - \bar x_i)'\beta + (\varepsilon_{it} - \bar\varepsilon_i), \qquad \mathrm{plim}_{N \to \infty}\,(\hat\rho_{FE} - \rho) \approx -\,\frac{1 + \rho}{T - 1}$$

In words: within-demeaning differences out the individual fixed effect and identifies $\beta$ by comparing each unit with itself, which requires strict exogeneity; once the right-hand side contains a lagged dependent variable, demeaning mixes the error into the lagged term and produces the Nickell bias, which is downward, of order $O(1/T)$, and does not vanish as $N$ grows.
:::

::: {.theorem data-label="Dynamic panel GMM (Arellano-Bond / Blundell-Bond) · Chapter 7"}
$$\text{difference GMM}:\ \mathbb{E}[y_{is}\, \Delta\varepsilon_{it}] = 0,\ s \leq t-2 \qquad \text{system GMM}:\ \mathbb{E}[\Delta y_{i,t-1}\,(\alpha_i + \varepsilon_{it})] = 0$$

In words: first-difference to remove $\alpha_i$, then use deep lagged levels as instruments to estimate the persistence $\rho$; when $\rho$ is close to 1 lagged levels are weak instruments, and Blundell-Bond restores efficiency by adding lagged differences as instruments for the levels equation, at the cost of a mean-stationarity assumption.
:::

::: {.theorem data-label="Kernel regression (Nadaraya-Watson / local linear) · Chapter 8"}
$$\hat m(x_0) = \frac{\sum_{i} y_i\, K\!\big(\frac{x_0 - x_i}{h}\big)}{\sum_{i} K\!\big(\frac{x_0 - x_i}{h}\big)}, \qquad \min_{a, b}\ \sum_{i} K\!\Big(\frac{x_0 - x_i}{h}\Big)\,\big(y_i - a - b\,(x_i - x_0)\big)^2$$

In words: a distance-decaying weighted local average, or local weighted least squares, estimates the whole curve $m(x) = \mathbb{E}[y \mid x]$, with the bandwidth $h$ as the smoothing knob; local linear corrects the boundary bias of Nadaraya-Watson.
:::

::: {.theorem data-label="Bias-variance and the optimal bandwidth · Chapter 8"}
$$\text{Bias} \approx \frac{h^2}{2}\, m''(x_0)\!\int\! u^2 K(u)\, du, \quad \text{Var} \approx \frac{1}{nh}\, \frac{\sigma^2}{f(x_0)}\!\int\! K(u)^2 du, \quad h^\ast \propto n^{-1/5}$$

In words: bias grows with $h^2$ and variance falls with $1/(nh)$, the MSE is U-shaped, and the optimal one-dimensional rate $n^{-4/5}$ is slower than the parametric $n^{-1}$; at the optimal bandwidth bias and standard deviation are of the same order, so you must undersmooth for the confidence interval to cover correctly.
:::

::: {.theorem data-label="Robinson partially linear (double residual) · Chapter 8"}
$$y_i - \mathbb{E}[y_i \mid x_i] = \beta\,\big(z_i - \mathbb{E}[z_i \mid x_i]\big) + \varepsilon_i, \qquad \sqrt{n}\,(\hat\beta - \beta) \xrightarrow{d} \mathcal{N}(0, V)$$

In words: partial out the nuisance function $g(x)$ nonparametrically, then regress the two sets of residuals to estimate the $\beta$ you care about, so that the slow convergence of $g$ does not drag down the $\sqrt{n}$ rate for $\beta$; this is precisely the forerunner of DML in Chapter 16.
:::

### Part III Causal Identification Toolbox

Part III gives one family of estimands per module and is the workhorse of reduced-form causal research. The identifying assumptions are all about counterfactuals and mostly untestable, so the cards give only the conclusion in closed form.

::: {.theorem data-label="Selection-bias decomposition · Chapter 9"}
$$\underbrace{\mathbb{E}[Y_i(1) \mid D_i = 1] - \mathbb{E}[Y_i(0) \mid D_i = 0]}_{\text{observed difference}} = \underbrace{\mathbb{E}[Y_i(1) \mid D_i = 1] - \mathbb{E}[Y_i(0) \mid D_i = 1]}_{\text{ATT}} + \underbrace{\mathbb{E}[Y_i(0) \mid D_i = 1] - \mathbb{E}[Y_i(0) \mid D_i = 0]}_{\text{selection bias}}$$

In words: the naive difference in group means equals the ATT plus a selection-bias term, the latter being the difference the two groups would have shown in the untreated state anyway, which can be nonzero even under homogeneous treatment effects, and randomization exists precisely to zero it out.
:::

::: {.theorem data-label="Strong ignorability (unconfoundedness + overlap) · Chapter 9"}
$$\{Y_i(1), Y_i(0)\} \perp D_i \mid X_i, \qquad 0 < e(X_i) < 1, \quad e(X) \equiv \Pr(D = 1 \mid X)$$

In words: after controlling for $X$, treatment is independent of the potential outcomes (everything that affects both adoption and outcome is in $X$), and every type of unit has positive probability of both treatment states; the former is untestable, the latter can be checked.
:::

::: {.theorem data-label="IPW and AIPW (doubly robust) · Chapter 9"}
$$\hat\tau_{HT} = \frac{1}{N} \sum_i \frac{D_i Y_i}{\hat e(X_i)} - \frac{1}{N} \sum_i \frac{(1 - D_i) Y_i}{1 - \hat e(X_i)}$$

$$\hat\tau_{AIPW} = \frac{1}{N} \sum_{i} \left\{ \hat\mu_1(X_i) - \hat\mu_0(X_i) + \frac{D_i\,(Y_i - \hat\mu_1(X_i))}{\hat e(X_i)} - \frac{(1 - D_i)\,(Y_i - \hat\mu_0(X_i))}{1 - \hat e(X_i)} \right\}$$

In words: IPW reconstructs the two counterfactual means by weighting on the inverse of the probability of being observed; AIPW layers a further outcome model $\hat\mu_d$ on top, and is consistent as long as one of the propensity and outcome models is correctly specified (doubly robust).
:::

::: {.theorem data-label="Research-design checklist · Chapter 10"}
$$\mathcal D=\{\text{unit},\text{treatment},\text{version},\text{timing},\text{exposure},\text{outcome},\text{population},\text{assignment}\}$$

In words: before running a regression, fix the eight design objects, and in particular distinguish assignment, actual adoption, and algorithmic exposure; when adding a control, recheck whether the estimand has changed and whether the variable is realized after treatment.
:::

::: {.theorem data-label="Prediction-powered rectification · Chapter 11"}
$$\widehat\beta_{PPI}=\widehat\beta_{pred,all}+\big(\widehat\beta_{gold,labeled}-\widehat\beta_{pred,labeled}\big)$$

In words: a large volume of machine predictions supplies precision, and a random gold sample estimates the bias between prediction and truth in the target estimating equation; validity comes from the validation design, not from treating model accuracy as a guarantee of unbiasedness.
:::

::: {.theorem data-label="IV assumptions and Wald / LATE · Chapter 12"}
$$\beta_{IV} = \frac{\mathbb{E}[Y_i \mid Z_i = 1] - \mathbb{E}[Y_i \mid Z_i = 0]}{\mathbb{E}[D_i \mid Z_i = 1] - \mathbb{E}[D_i \mid Z_i = 0]} = \mathbb{E}\big[Y_i(1) - Y_i(0) \mid \text{complier}\big]$$

In words: under the four conditions of relevance, random assignment, exclusion, and monotonicity, the Wald ratio (reduced form divided by first stage) identifies the LATE for the complier subpopulation, and only relevance is directly testable.
:::

::: {.theorem data-label="MTE (marginal treatment effect) · Chapter 12"}
$$D_i = \mathbf{1}\{U_{D,i} \leq P(Z_i)\}, \qquad \Delta^{MTE}(u) = \mathbb{E}\big[Y_i(1) - Y_i(0) \mid U_{D,i} = u\big]$$

$$\text{ATE} = \int_0^1 \Delta^{MTE}(u)\, du, \qquad \text{ATT} = \int_0^1 \Delta^{MTE}(u)\, \omega_{ATT}(u)\, du$$

In words: the MTE is the effect for the marginal population whose "resistance" $U_D = u$ sits exactly at the threshold, ATE, ATT, LATE, and OLS are all integrals of it under different weights, and you need $P(Z)$ to have enough support to trace out the whole curve.
:::

::: {.theorem data-label="2×2 DiD · Chapter 13"}
$$\Delta_{DD} = \big(\mathbb{E}[Y_{i2} \mid D=1] - \mathbb{E}[Y_{i1} \mid D=1]\big) - \big(\mathbb{E}[Y_{i2} \mid D=0] - \mathbb{E}[Y_{i1} \mid D=0]\big) = ATT$$

In words: the before-after change in the treated group minus the contemporaneous change in the control group identifies the ATT under parallel trends, no anticipation, and SUTVA, with each of the three assumptions blocking one route of bias.
:::

::: {.theorem data-label="Goodman-Bacon decomposition · Chapter 13"}
$$\mathrm{plim}\; \hat\beta_{TWFE} = VWATT + VWCT - \Delta ATT$$

In words: in a staggered design the static TWFE is a variance-weighted average of all the 2×2 comparisons, the $\Delta ATT$ term comes from the forbidden comparison that uses already-treated units as controls, and when effects grow over time it systematically pulls the estimate down or even flips its sign.
:::

::: {.theorem data-label="Modern staggered estimators (CS / SA / BJS) · Chapter 13"}
$$ATT(g,t) = \mathbb{E}\big[Y_t - Y_{g-1} \mid G = g\big] - \mathbb{E}\big[Y_t - Y_{g-1} \mid C = 1\big]$$

In words: first break the effect down to the finest group-time $ATT(g,t)$, using only clean controls in each cell (never-treated or not-yet-treated), then aggregate under transparent weights; Callaway-Sant'Anna estimates directly, Sun-Abraham uses a saturated interaction regression, and Borusyak-Jaravel-Spiess uses imputation to fill in the counterfactual.
:::

::: {.theorem data-label="Sharp / fuzzy RDD · Chapter 14"}
$$\tau_{\mathrm{RD}} = \lim_{x\downarrow c}\mathbb{E}[Y_i \mid X_i = x] - \lim_{x\uparrow c}\mathbb{E}[Y_i \mid X_i = x] = \mathbb{E}[\tau_i \mid X_i = c]$$

$$\tau_{\mathrm{FRD}}(c) = \frac{\lim_{x\downarrow c}\mathbb{E}[Y_i \mid X_i = x] - \lim_{x\uparrow c}\mathbb{E}[Y_i \mid X_i = x]}{\lim_{x\downarrow c}\Pr(D_i = 1 \mid X_i = x) - \lim_{x\uparrow c}\Pr(D_i = 1 \mid X_i = x)}$$

In words: at the cutoff of the running variable the conditional expectations of the potential outcomes are continuous, so the jump in the outcome is the local effect; in the sharp case treatment flips deterministically, in the fuzzy case only the probability jumps, and the effect equals the outcome jump divided by the treatment jump, which is the LATE for compliers at the cutoff. Estimate with local linear plus bias-corrected robust intervals.
:::

::: {.theorem data-label="Synthetic control weights · Chapter 14"}
$$\hat{\mathbf{w}} = \arg\min_{\mathbf{w}}\ \Big(\mathbf{X}_1 - \textstyle\sum_j w_j \mathbf{X}_j\Big)' \boldsymbol\Omega \Big(\mathbf{X}_1 - \textstyle\sum_j w_j \mathbf{X}_j\Big) \quad\text{s.t. } \sum_{j} w_j = 1,\ w_j \geq 0$$

In words: for a single treated unit, find a set of non-negative weights summing to one in the donor pool so that the synthetic unit's pre-period path tracks the treated unit, and the post-period gap $Y_{1t} - \sum_j \hat w_j Y_{jt}$ is the effect; the simplex constraint ensures the convex combination does not extrapolate.
:::

::: {.theorem data-label="OVB sensitivity and the robustness value · Chapter 15"}
$$R^2_{D\sim U\mid X},\qquad R^2_{Y\sim U\mid D,X},\qquad \Theta(B)=[\widehat\tau-B,\widehat\tau+B]$$

In words: two partial $R^2$s express the strength of an unobserved confounder with the treatment and with the outcome; the robustness value gives the minimum equal-strength confounding needed to overturn the conclusion, and allowing a bias bound $B$ expands the point estimate into an identified region.
:::

::: {.theorem data-label="DML: Neyman orthogonality and cross-fitting · Chapter 16"}
$$\psi(W; \theta, \ell, m) = \big(Y - \ell(X) - \theta (D - m(X))\big)\big(D - m(X)\big), \qquad \hat\theta = \frac{\sum_i \tilde D_i\, \tilde Y_i}{\sum_i \tilde D_i^2}$$

In words: residualize both $Y$ and $D$ on $X$ ($\tilde Y, \tilde D$) and then regress, the orthogonal score is first-order insensitive to the nuisance, and paired with cross-fitting, as long as $\|\hat\ell - \ell_0\| \cdot \|\hat m - m_0\| = o_p(n^{-1/2})$ you get $\sqrt{n}$ normality, so the nuisance can be estimated by any ML method.
:::

::: {.theorem data-label="Causal forest / CATE · Chapter 16"}
$$\hat\tau(x) = \frac{\sum_i \alpha_i(x)\, (D_i - \bar D_x)(Y_i - \bar Y_x)}{\sum_i \alpha_i(x)\, (D_i - \bar D_x)^2}$$

In words: the random forest splits on treatment-effect heterogeneity, the forest weights $\alpha_i(x)$ are a data-adaptive similarity, and a local residual-on-residual regression over them yields $\tau(x)$; honesty (separating the splitting sample from the estimation sample) buys pointwise asymptotic normality and confidence intervals.
:::

::: {.theorem data-label="Power / MDE / sample size · Chapter 17"}
$$\text{power} \approx \Phi\!\left(\frac{\tau}{\text{SE}} - z_{1-\alpha/2}\right), \quad \text{MDE} = (z_{1-\alpha/2} + z_{1-\beta})\, \sigma\sqrt{\tfrac{2}{n}}, \quad n = \frac{2\sigma^2 (z_{1-\alpha/2} + z_{1-\beta})^2}{\tau^2}$$

In words: whether an experiment can detect an effect is jointly determined by sample size, variance, and the true effect, the MDE is the smallest effect reliably detectable at a given power, and sample size grows in inverse proportion to the square of the target effect; ex-post observed power is a monotone transform of the p-value and carries no information.
:::

::: {.theorem data-label="CUPED (variance reduction) · Chapter 17"}
$$Y_i^{cv} = Y_i - \theta (X_i - \bar X), \quad \theta = \frac{\text{Cov}(Y, X)}{\text{Var}(X)}, \quad \text{Var}(Y^{cv}) = \text{Var}(Y)\,(1 - \rho^2)$$

In words: residualizing the outcome on a pre-experiment covariate reduces variance by a factor $1 - \rho^2$ (equivalent to ANCOVA), and the stronger the correlation the more the effective sample size multiplies; the covariate must come from before the experiment, otherwise it is a bad control.
:::

::: {.theorem data-label="Direct and overall effects under interference · Chapter 18"}
$$\tau_D(g)=\mathbb E[Y_i(1,g)-Y_i(0,g)],\qquad \tau_O=\mathbb E[Y_i(1,g_1)-Y_i(0,g_0)]$$

In words: individual randomization typically identifies the direct contrast at a fixed saturation $g$; a platform-wide rollout changes both a unit's own treatment and its exposure environment at once, and identifying the corresponding overall effect requires cluster, two-stage saturation, or switchback variation.
:::

### Part IV Structural Methods

Part IV consumes the theory of Part I and the identification of Part III, structuring a model to study substitution, counterfactuals, and welfare. Here "identification" often rests on demand instruments, conduct assumptions, or equilibrium selection, things that cannot be tested from the first-order conditions themselves.

::: {.theorem data-label="Logit choice probability and IIA · Chapter 19"}
$$s_{ij} = \frac{e^{V_{ij}}}{1 + \sum_{k=1}^{J} e^{V_{ik}}}, \qquad \frac{s_{ij}}{s_{ik}} = e^{V_{ij} - V_{ik}}$$

In words: a Type-I extreme value shock gives closed-form shares, and the choice probability is the normalization of the utility exponentials; the price is IIA, the ratio of any two probabilities is independent of a third, which forces proportional substitution (red-bus/blue-bus). Coefficients are identified only up to a scale normalization, and what is comparable are ratios such as WTP.
:::

::: {.theorem data-label="Nested logit (GEV) · Chapter 19"}
$$s_{ij} = \frac{e^{V_{ij}/\lambda_g}\Big(\sum_{k \in g} e^{V_{ik}/\lambda_g}\Big)^{\lambda_g - 1}}{\sum_{h}\Big(\sum_{k \in h} e^{V_{ik}/\lambda_h}\Big)^{\lambda_h}}$$

In words: it allows the shocks of products in the same nest to be correlated, with the nest parameter $\lambda_g \in (0,1]$, $\lambda_g = 1$ collapsing back to logit and $\lambda_g \to 0$ making within-nest products nearly perfect substitutes; the nesting structure must be specified in advance, and getting the nests wrong ruins everything.
:::

::: {.theorem data-label="Mixed / random-coefficients logit · Chapter 19"}
$$s_{ij}(\theta) = \int \frac{e^{\,\beta\, x_{ij}}}{\sum_k e^{\,\beta\, x_{ik}}}\, dF(\beta \mid \theta)$$

In words: letting the coefficients $\beta_i$ vary across the population by a distribution $F$ and integrating over the taste distribution can approximate any RUM and break IIA; the random-coefficient distribution is identified by rich choice-set and characteristic variation, and the repeated choices of a panel can sharpen identification substantially through cross-situation correlation, though this is not in principle necessary.
:::

::: {.theorem data-label="Berry (1994) inversion · Chapter 20"}
$$\ln s_{jt} - \ln s_{0t} = \delta_{jt} = x_{jt}\beta - \alpha p_{jt} + \xi_{jt}$$

In words: invert the observed market shares into the mean utility $\delta_{jt}$, which under logit is just the difference in log shares, turning the problem into a linear IV problem with an endogenous price, and wherever shares appear an instrument is needed.
:::

::: {.theorem data-label="BLP contraction and GMM · Chapter 20"}
$$\delta_t^{(h+1)} = \delta_t^{(h)} + \ln S_t^{\text{obs}} - \ln s_t\big(\delta_t^{(h)}, \theta_2\big), \qquad \hat\theta_2 = \arg\min_{\theta_2}\ \big(Z' \xi(\theta_2)\big)' W \big(Z' \xi(\theta_2)\big)$$

In words: under random coefficients the shares have no analytic inversion, so for a fixed $\theta_2$ a contraction mapping is iterated to solve for $\delta_t$ (modulus below 1, converging from any starting value), and GMM is then run on the orthogonality of instruments with $\xi$; the inner tolerance must be set very tight (say $10^{-13}$) to keep it from seeping into the outer loop.
:::

::: {.theorem data-label="Multiproduct Bertrand FOC and markup (ℋ∘Δ)⁻¹ · Chapter 21"}
$$mc = \mathbf{p} - \big(\mathcal{H} \odot \Delta(\mathbf{p})\big)^{-1} s(\mathbf{p})$$

In words: given demand, the multiproduct Bertrand first-order condition backs marginal cost out exactly, where $\mathcal{H}$ is the ownership matrix, $\Delta$ the negative matrix of demand price derivatives, and $\odot$ the Hadamard product; cost is inverted rather than estimated, which is both convenient and dangerous, since a misspecified demand or conduct propagates through the FOC and is amplified systematically.
:::

::: {.theorem data-label="Diversion ratio and GUPPI / UPP · Chapter 21"}
$$D_{j \to k} = \frac{\partial s_k / \partial p_j}{-\,\partial s_j / \partial p_j}, \qquad GUPPI_j = \frac{(p_k - mc_k)\, D_{j \to k}}{p_j}$$

In words: the diversion ratio is how much of the volume $j$ loses to a price increase goes to $k$ (equal to the second-choice share), and is the core demand output for merger analysis; GUPPI combines it with the markup into a dimensionless first-pass merger screen, and UPP then subtracts an efficiency offset.
:::

::: {.theorem data-label="Merger simulation and log-sum welfare · Chapter 21"}
$$\Delta CS = \int \frac{1}{\alpha_i}\Big[\log\big(1 + \textstyle\sum_j e^{V_{ij}^{\text{post}}}\big) - \log\big(1 + \textstyle\sum_j e^{V_{ij}^{\text{pre}}}\big)\Big]\, dF_i$$

In words: a merger flips the corresponding cross-terms of the ownership matrix from 0 to 1, and re-solving the FOC with $mc$ and demand held fixed gives the post-merger prices; the change in consumer surplus is the difference in the log-sum expected maximum utility before and after, divided by the marginal utility of price $\alpha_i$.
:::

::: {.theorem data-label="Entry games and Bresnahan-Reiss · Chapter 22"}
$$\Pi_{fm} = X_m \beta + Z_{fm} \alpha - \delta\, N_{-f,m} + \xi_m + \varepsilon_{fm}, \qquad a_{fm} = \mathbf{1}\{\Pi_{fm} \geq 0\}$$

In words: a firm enters when its expected profit is non-negative, with $\delta > 0$ the intensity of competition; Bresnahan-Reiss sidesteps the multiplicity of the game's equilibria by running an ordered probit on the number of entrants only, and uses the ratio of adjacent per-capita entry thresholds $s_{N+1}/s_N$ to gauge how competition intensifies with the number of incumbents.
:::

::: {.theorem data-label="Partial identification (Ciliberto-Tamer) · Chapter 22"}
$$H_1(\theta, X) \leq P(y \mid X) \leq H_2(\theta, X), \qquad Q(\theta) = \sum_y \Big[ \big(P(y) - H_1\big)_-^2 + \big(P(y) - H_2\big)_+^2 \Big]$$

In words: when the entry game has multiple equilibria the model gives only an interval for the outcome probability ($H_1$ the lower bound from the unique-equilibrium probability, $H_2$ the upper bound from some equilibrium), the parameters where the criterion function is zero form the identified set, and no assumption need be imposed on the unknown equilibrium-selection mechanism.
:::

::: {.theorem data-label="Bellman and the Rust contraction · Chapter 23"}
$$V(x) = \mathbb{E}_\varepsilon \max_{i}\Big\{ u(x, i) + \varepsilon_i + \beta\, \mathbb{E}\big[V(x') \mid x, i\big] \Big\}$$

$$(T_\theta V)(x) = \gamma + \log\Big[ \exp\big(u(x,0) + \beta\, \mathbb{E}[V(x') \mid x, 0]\big) + \exp\big(u(x,1) + \beta\, \mathbb{E}[V(x') \mid x, 1]\big) \Big]$$

In words: the value function is its own fixed point, and under conditional independence plus Type-I EV the Bellman equation is written as the log-sum operator $T_\theta$ (the leading term $\gamma$ is the Euler-Mascheroni constant, arising from the expected maximum $\mathbb{E}_\varepsilon\max$ and not affecting the choice probability), which is a contraction with modulus $\beta < 1$ and iterates to a unique solution; NFXP maximizes the likelihood in the outer loop and solves the fixed point in the inner loop.
:::

::: {.theorem data-label="Hotz-Miller CCP inversion · Chapter 23"}
$$v(x, 1) - v(x, 0) = \log P(1 \mid x) - \log P(0 \mid x)$$

In words: under Type-I EV the value difference between two actions maps one-to-one to the ratio of choice probabilities, so the dynamic value difference can be recovered from observed conditional choice probabilities without solving the fixed point, which is the key to two-step estimation; the discount factor $\beta$ is generally not identified and must be fixed or pinned down by a separate exclusion restriction.
:::

### Part V Topics in Digital Economics

Part V brings the earlier three legs together onto the phenomena of the digital economy: two-sided platforms, consumer search, algorithms, and data. The methods are mostly combinations and applications of the earlier chapters.

::: {.theorem data-label="Two-sided pricing (Armstrong / Rochet-Tirole) · Chapter 24"}
$$p_i = f_i - \alpha_j\, n_j + \frac{n_i}{\partial n_i / \partial u_i}$$

In words: the platform's optimal price for one group of members equals its own connection cost, minus the cross-side external benefit $\alpha_j n_j$ this group brings to the other, plus the usual monopoly markup; the middle term can price one side below cost or even subsidize it, and what matters is the price structure (how the total price is split across the two sides), not the price level.
:::

::: {.theorem data-label="Network effects (Katz-Shapiro) and tipping · Chapter 24"}
$$P(q) = r(q) + v(q), \qquad P'(q^\ast) = r'(q^\ast) + v'(q^\ast) < 0$$

In words: the marginal user's willingness to pay equals a decreasing stand-alone value plus an increasing network value, the fulfilled-expectations demand can be non-monotone, and where it crosses the price line at several points there are multiple equilibria; a stable equilibrium requires the demand curve to cross the price line from above, otherwise you get tipping and a critical mass.
:::

::: {.theorem data-label="Switching costs / state dependence · Chapter 24"}
$$u_{ijt} = x_{jt}'\beta + \xi_{ij} + \delta\, \mathbf{1}\{y_{i,t-1} = j\} + \varepsilon_{ijt}$$

In words: $\delta$ is genuine structural state dependence (a switching cost, monetized as $\delta / |\alpha|$ with $\alpha$ the price coefficient), which must be separated from the spurious state dependence created by persistent preferences $\xi_{ij}$; failing to control for heterogeneity systematically overstates $\delta$, and the initial-conditions problem must also be handled.
:::

::: {.theorem data-label="Sequential search (Weitzman reservation value) · Chapter 25"}
$$c_j = \mathbb{E}\big[\max(u_j - z_j, 0)\big]$$

In words: each option is assigned a reservation value $z_j$ by its own search-cost equation, the optimal rule reduces to a set of static indices, inspect in decreasing order of $z_j$ and stop once the realized utility exceeds every unopened index; the higher the search cost the lower the $z$ and the shallower the search.
:::

::: {.theorem data-label="Limited consideration sets and identifying search costs · Chapter 25"}
$$\Pr(\text{buy } j \mid \text{choice set} = \text{inspected set}_i)$$

In words: a full-information logit mistakes a "cheap item never inspected" for one "seen and rejected", attenuating the price coefficient toward zero; replacing the choice set with the actually-inspected set is directionally right but overstates because stopping is endogenous, and the correct fix is to model the search process explicitly, using randomized ordering (position affects purchase only through "whether it was inspected") to pin down search costs and preferences separately.
:::

::: {.theorem data-label="Algorithmic pricing and Q-learning (collusion) · Chapter 26"}
$$Q(s_t, a_t) \leftarrow (1 - \alpha)\, Q(s_t, a_t) + \alpha\Big[ r_t + \delta\, \max_{a'} Q(s_{t+1}, a') \Big], \qquad \Delta = \frac{\bar{p} - p^N}{p^M - p^N}$$

In words: with the rival's previous-period price as the state and the current profit as the reward, Q-learning bootstraps its updates of the action value; the collusion index $\Delta$ normalizes the average price between the competitive price $p^N$ and the monopoly price $p^M$, and in simulation the algorithm can spontaneously learn supra-competitive pricing, but "simulations collude" is not the same as "reality colludes".
:::

::: {.theorem data-label="AI adoption and organizational complementarity · Chapter 26"}
$$Y=F(K,L,A,O),\qquad F_{AO}>0,\qquad \tau_{AO}=(\bar Y_{11}-\bar Y_{01})-(\bar Y_{10}-\bar Y_{00})$$

In words: the cross-marginal product of AI capital $A$ with organizational capital $O$ (workflow, skills, data readiness, and the like) is positive; the cleanest empirical design is a factorial experiment that randomizes AI access and the organizational intervention separately, and the interaction contrast $\tau_{AO}$ identifies the policy-relevant complementarity.
:::

::: {.theorem data-label="Generated regressors: attenuation and correction / PPI · Chapter 26"}
$$\mathrm{plim}\; \hat\beta_{\text{OLS}} = \beta \cdot \frac{\sigma_{x^\ast}^2}{\sigma_{x^\ast}^2 + \sigma_u^2} = \beta \cdot \lambda, \qquad \hat\beta_{\text{corrected}} = \frac{\hat\beta_{\text{OLS}}}{\hat\lambda}$$

In words: putting a machine-predicted variable into the right-hand side as if it were the truth, classical measurement error attenuates the coefficient toward zero by the factor $\lambda$ (reliability), and this does not vanish as the sample grows; a small randomly drawn set of gold labels makes $\lambda$ or the systematic bias estimable and thus correctable (PPI uses a rectifier to subtract the bias), and the standard errors must fold in the uncertainty of the generation step.
:::

## 2 A Method Decision Tree

The first step in choosing a method is not to ask "which model is the most advanced" but "what do I want to answer, what does the data look like, and where does the variation come from". The tree below starts from the research goal, follows the data structure to the family of methods you should use, and labels each exit with the corresponding chapter. It is only an entry-level guide: the real choice always turns on a substantive judgment about the treatment-assignment mechanism and the data-generating process, and on the same data a different credible assumption can send you to a different family of methods.

First tell apart whether the goal is description and prediction, a causal effect, or counterfactuals and welfare:

- **You only want to describe, predict, or characterize a relationship curve** (no causal claim)
    - Want the data to decide the functional shape and avoid specification error → kernel estimation / local linear / splines / partially linear (Chapter 8)
    - The goal is itself a prediction-policy problem where "predicting well improves the decision" → supervised ML and its inference (Chapter 16), and guard against attenuation when a machine-generated variable enters a regression (Chapter 26)
- **You want the causal effect of some treatment** (treatment effect). First look at how treatment is assigned:
    - First use Chapter 10 to fix the unit, treatment, version, timing, exposure, outcome, population, and assignment; when a variable comes from platform logs or machine labels, first run the measurement audit of Chapter 11
    - Treatment is randomly assigned (RCT / online A/B test) → difference in means and randomization inference in the potential-outcomes framework (Chapter 9); experimental power, MDE, CUPED, sequential testing (Chapter 17); when treatment spills over to other units or the goal is a market-wide rollout, see interference and marketplace experiments (Chapter 18)
    - Observational data with non-random treatment, then look for the identifying handle:
        - You believe the confounders are all observed (selection on observables) → matching / IPW / AIPW (Chapter 9); high-dimensional covariates or complex nuisance → DML (Chapter 16)
        - Treatment is endogenous but there is an instrument that is relevant, exogenous, and satisfies exclusion and monotonicity → IV / LATE; want the marginal population and extrapolation weights → MTE and control functions (Chapter 12)
        - A continuous running variable that triggers treatment only past a cutoff → sharp RDD; a threshold that only changes the probability of treatment → fuzzy RDD (Chapter 14)
        - Panel data where units are treated at different times, with pre- and post-periods and a control group → DiD / event study; staggered with time-heterogeneous effects → Callaway-Sant'Anna / Sun-Abraham / BJS (Chapter 13)
        - Only one (or very few) treated aggregates, plus a set of untreated donors and a long enough pre-period → synthetic control or synthetic DiD (Chapter 14)
    - You care about heterogeneous effects (who benefits most) or the optimal assignment policy (whom to treat) → causal forest / GATES / policy tree (Chapter 16)
    - The core identifying assumption holds only approximately and you want to quantify how large a departure would overturn the conclusion → sensitivity, negative controls, bounds, and the breakdown frontier (Chapter 15)
    - One party's treatment changes the outcomes of other units, or the goal is a market-wide rollout → cluster / two-stage / switchback marketplace experiments (Chapter 18)
- **You want counterfactuals, substitution patterns, market power, or welfare** (structural modeling). First look at what data you have:
    - Individual-level choice data (who chose what among several options) → discrete choice, moving from logit to nested logit to mixed logit as the substitution pattern grows richer (Chapter 19)
    - Only market-level shares and prices, with endogenous prices → BLP demand: Berry inversion plus contraction plus demand instruments (Chapter 20)
    - Demand already estimated, and you want to recover markups, marginal costs, simulate a merger, or test conduct → supply-side first-order conditions and merger simulation (Chapter 21)
    - You want to explain how many firms are in a market and the entry-exit decision → the Bresnahan-Reiss threshold model; with multiple equilibria, switch to partial identification (Ciliberto-Tamer) (Chapter 22)
    - The decision is dynamic, where today's choice changes tomorrow's state (replacing equipment, adopting technology, optimal stopping) → single-agent Rust NFXP or Hotz-Miller CCP (Chapter 23)
- **Landing on a specific phenomenon of the digital economy** (combining the tools above)
    - Two-sided platform pricing, network effects and tipping, switching costs and lock-in → Chapter 24
    - Consumer search frictions, limited consideration sets, advertising and information goods → Chapter 25
    - Algorithmic pricing, AI adoption, organizational complements, inference with machine-generated variables and the economics of data → Chapter 26

Whichever path you take, the underlying machinery of estimation and inference lives in Part II: the unified framework of extremum estimation and GMM (Chapter 5), the choice of standard errors, that is, the cluster / HAC / bootstrap / few-clusters first-aid kit (Chapter 6), and panel fixed effects together with the Nickell bias of dynamic panels (Chapter 7). The language of counterfactuals and welfare comes from Part I: consumer surplus, CV/EV, and market power (Chapters 1 and 4).

For a quick lookup by data shape, use the table below:

| Data / situation | First-choice method | Chapter |
|---|---|---|
| Randomized experiment, online A/B test | difference in means + randomization inference; power, CUPED, sequential | Chapters 9, 17 |
| Platform logs / machine labels | measurement audit, random validation, PPI | Chapter 11 |
| Observational data, confounders observed | matching / IPW / AIPW; DML in high dimensions | Chapters 9, 16 |
| Endogenous treatment + valid instrument | IV / LATE / MTE / control function | Chapter 12 |
| running variable + cutoff | sharp / fuzzy RDD | Chapter 14 |
| Panel + staggered treatment timing | DiD / event study; CS / SA / BJS when staggered | Chapter 13 |
| Single treated aggregate + donor pool | synthetic control / synthetic DiD | Chapter 14 |
| Assumptions hold only approximately / hidden bias | sensitivity, negative controls, partial ID | Chapter 15 |
| Want heterogeneous effects / policy | causal forest / policy learning | Chapter 16 |
| Spillover / marketplace rollout | cluster, two-stage, switchback | Chapter 18 |
| Individual choice data, want substitution elasticities | logit / nested / mixed logit | Chapter 19 |
| Market shares + endogenous prices | BLP demand | Chapter 20 |
| Demand in hand, want markups / mergers | supply-side FOC / merger simulation | Chapter 21 |
| Entry / market structure | Bresnahan-Reiss / partial identification | Chapter 22 |
| Dynamic single-agent decisions | Rust NFXP / Hotz-Miller CCP | Chapter 23 |

## 3 Cross-Software Command Map

Each row of the table below is one method, giving the main commands in R, Python, and Stata, with R as the primary. All commands are real packages and functions. A blank cell means there is no widely used standard package, and the corresponding chapter mostly implements the formula by hand (base R's `optim`, `glm`, and matrix algebra suffice), so we leave it blank rather than invent an API. Package and function names may change across versions, so check the respective documentation before use.

| Method (chapter) | R | Python | Stata |
|---|---|---|---|
| OLS / high-dimensional fixed effects | `fixest::feols`, `lfe::felm`, `plm::plm` | `linearmodels` PanelOLS, `pyfixest`, `statsmodels` OLS | `reg`, `reghdfe`, `areg`, `xtreg` |
| MLE / GLM (logit, probit, Poisson) (5) | `glm`, `stats4::mle`, `bbmle::mle2` | `statsmodels` Logit/GLM, `scipy.optimize` | `logit`, `probit`, `glm`, `ml` |
| GMM (5) | `gmm::gmm`, `momentfit` | `linearmodels` IVGMM, `statsmodels` GMM | `gmm` |
| Robust / cluster standard errors (6) | `sandwich::vcovHC`, `sandwich::vcovCL` | `statsmodels` cov_type='cluster', `linearmodels` | `vce(robust)`, `vce(cluster)` |
| HAC standard errors (6) | `sandwich::NeweyWest` | `statsmodels` cov_type='HAC' | `newey` |
| Bootstrap (6) | `boot::boot` | `scipy.stats.bootstrap`, `arch` | `bootstrap` |
| Delta method (6) | `car::deltaMethod`, `msm::deltamethod` |  | `nlcom` |
| Panel fixed effects / RE (7) | `plm::plm`, `fixest::feols`, `plm::phtest` | `linearmodels` PanelOLS/RandomEffects | `xtreg`, `xtreg, re` |
| Dynamic panel GMM (7) | `plm::pgmm` | `pydynpd` | `xtabond2`, `xtdpdgmm`, `xtabond` |
| Nonparametric / kernel regression (8) | `np::npreg`, `KernSmooth::locpoly`, `smooth.spline`, `mgcv::gam` | `statsmodels` KernelReg/lowess | `npregress`, `lpoly` |
| Matching (9) | `MatchIt::matchit`, `Matching::Match` |  | `teffects psmatch`, `psmatch2`, `kmatch` |
| IPW / AIPW (9) | `WeightIt::weightit`, `AIPW`, `drtmle` | `DoubleML` IRM, `econml` | `teffects ipw`, `teffects aipw`, `teffects ipwra` |
| Validation / PPI (11) | base estimating equations, `ppi` implementations | `ppi_py` |  |
| IV / 2SLS / LATE (12) | `fixest::feols`, `AER::ivreg`, `estimatr::iv_robust` | `linearmodels` IV2SLS | `ivregress 2sls`, `ivreg2` |
| MTE (12) | `ivmte`, `localIV` |  | `mtefe` |
| DiD / event study (13) | `fixest::feols`, `fixest::sunab`, `did::att_gt` | `differences`, `pyfixest` | `csdid`, `eventstudyinteract`, `reghdfe` |
| Goodman-Bacon decomposition (13) | `bacondecomp::bacon` |  | `bacondecomp` |
| Staggered imputation (13) | `didimputation::did_imputation`, `did2s::did2s` |  | `did_imputation`, `did2s` |
| RDD (14) | `rdrobust::rdrobust`, `rddensity::rddensity` | `rdrobust`, `rddensity` | `rdrobust`, `rddensity` |
| Synthetic control / SDID (14) | `Synth::synth`, `tidysynth`, `augsynth`, `synthdid` | `pysyncon`, `SparseSC` | `synth`, `synth_runner`, `sdid` |
| OVB sensitivity / honest DiD (15) | `sensemakr`, `HonestDiD` |  |  |
| DML (16) | `DoubleML`, `grf` | `DoubleML`, `econml` | `ddml`, `pdslasso` |
| Causal forest (16) | `grf::causal_forest` | `econml` CausalForestDML |  |
| Policy learning (16) | `policytree::policy_tree` | `econml` policy |  |
| Power / MDE (17) | `pwr::pwr.t.test`, `DeclareDesign` | `statsmodels.stats.power` | `power twomeans` |
| Sequential test boundaries (17) | `gsDesign` |  |  |
| Marketplace / switchback (18) | `DeclareDesign` or hand-coded by assignment |  |  |
| Logit / conditional logit (19) | `mlogit::mlogit`, `survival::clogit`, `apollo` | `xlogit`, `pylogit`, `biogeme` | `clogit`, `asclogit`, `mlogit` |
| Nested logit (19) | `mlogit::mlogit`, `apollo`, `gmnl::gmnl` | `biogeme`, `pylogit` | `nlogit` |
| Mixed / RC logit (19) | `mlogit::mlogit`, `mixl`, `gmnl::gmnl`, `apollo` | `xlogit`, `biogeme` | `mixlogit`, `cmxtmixlogit` |
| BLP demand (20) | `BLPestimatoR` | `pyblp` |  |
| Markups / merger simulation (21) | `antitrust` | `pyblp` | `mergersim` |
| Entry / structural partial identification (22) |  |  |  |
| Dynamic structural Rust / CCP (23) |  | `ruspy` |  |
| RE / dynamic logit (24) | `lme4::glmer` |  | `xtlogit, re` |
| Consumer search (25) | `survival::clogit` |  |  |
| Measurement error / generated regressors (26) | `simex` | `ppi_py` |  |

## 4 Common Pitfalls Checklist

Below we gather the most common errors from each chapter's "failure modes" section, grouped by method, for a line-by-line self-check before and after you get your hands dirty. Each item names only the problem and the direction of the fix; the full argument is in the corresponding chapter.

### Part II Econometric Core

- Specification error (Chapter 5): distinguish "mean right, distribution wrong" (QMLE is still consistent, fix the standard errors with the sandwich) from "mean wrong" (even consistency is gone, and robust standard errors are meaningless); report sandwich standard errors by default rather than Fisher information (the latter requires the information-matrix equality to hold); complete separation in logit makes the MLE diverge, and the software may not warn you.
- Weak identification (Chapter 5): when the Jacobian or the curvature is near zero the variance explodes, a large sample cannot save you, and the judgment must come from "where the identifying variation comes from" rather than from significance stars.
- The level of the standard errors (Chapter 6): heteroskedasticity-robust standard errors do nothing about within-cluster correlation, so cluster at the level of treatment assignment; the asymptotics run on the number of clusters $G$, and a small $G$ (empirically below 30 to 40, and worse when few clusters are treated) over-rejects, so upgrade to $t_{G-1}$, CR2, or the wild cluster bootstrap.
- Delta method (Chapter 6): when the denominator is near zero (the Fieller problem) the symmetric normal approximation fails, so switch to percentile-t or Fieller intervals; it does not apply where the function is non-differentiable or the gradient is zero.
- Bootstrap (Chapter 6): on correlated data a naive bootstrap amounts to assuming iid, so use the block or wild cluster bootstrap; boundary and extreme-value statistics are not asymptotically normal, so use subsampling.
- Generated regressors (Chapter 6): two-step estimation understates the standard errors, so use the Murphy-Topel variance or bootstrap the whole pipeline.
- Panel Nickell bias (Chapter 7): in a short panel, a lagged dependent variable together with fixed effects produces a downward bias of order $O(1/T)$ that does not vanish as $N$ grows, and the regression runs without complaint; switch to dynamic panel GMM.
- Too many instruments in dynamic panels (Chapter 7): the number of instruments grows with the square of $T$, dragging GMM back toward FE and destroying the power of the Sargan test (almost never rejecting is a warning sign, not a reassurance), so limit the lag depth or collapse the instruments, and report the instrument count.
- Nonparametric smoothing and inference (Chapter 8): too large a bandwidth smooths away real bends and too small chases noise; the boundary and sparse regions are the least trustworthy (kernel regression has boundary bias, so use local linear); at the optimal bandwidth bias and standard deviation are of the same order, so you must undersmooth for the confidence interval to cover; beyond three or four continuous regressors the curse of dimensionality bites, so move to semiparametrics.

### Part III Causal Identification Toolbox

- Selection on observables / matching (Chapter 9): the CIA is untestable, and balance and placebo checks are only indirect circumstantial support; when overlap is poor the IPW weights explode, and trimming changes the estimand; the bootstrap is invalid for fixed-number NN matching, so use the Abadie-Imbens variance; do not control for downstream variables that respond to treatment (bad controls).
- Research design (Chapter 10): first fix the unit, treatment, version, timing, exposure, outcome, population, and assignment; post-treatment controls can both change the estimand and open a collider path; the arrows in a DAG come from institutional assumptions, not facts that a correlation matrix verifies automatically.
- Measurement (Chapter 11): high predictive accuracy does not guarantee downstream inference; the gold sample must represent the target population and carry reliable labels; the bias of a machine variable differs on the left-hand side and the right-hand side; PPI fixes measurement, not treatment endogeneity.
- IV / LATE (Chapter 12): of the four assumptions only relevance is testable, and exclusion is the one most in need of argument and the most hidden; weak instruments are the number-one threat, so report the first-stage F (the effective F under heteroskedasticity or clustering) and base the main inference on Anderson-Rubin; feeding the fitted values of a nonlinear first stage straight into the second stage is a forbidden regression; the LATE is only the effect for compliers, so do not extrapolate it as the ATE.
- DiD (Chapter 13): under staggered timing the weights of static TWFE are contaminated by forbidden comparisons, biasing downward or even flipping sign when effects grow, which cannot be "patched a little", so run a Goodman-Bacon decomposition or switch to a modern estimator; flat pre-trends are supportive evidence of low power rather than proof; throwing in a lagged dependent variable to "absorb the trend" brings the Nickell bias; cluster at the level of treatment assignment, and few treated clusters are the most dangerous.
- RDD (Chapter 14): do not use high-order global polynomials (they give extreme weight to distant observations, are sensitive to the order, cover poorly, and suffer the Runge phenomenon), use local linear with a bandwidth; test manipulation and sorting with the McCrary density test plus covariate continuity, and use a donut for integer heaping; the effect is a local quantity at the cutoff; hanging several programs on the same threshold treats their combined force as one effect.
- Synthetic control (Chapter 14): if the pre-period fit is poor, do not use it (switch to augmented SC); the treated unit must lie inside the convex hull of the donors, otherwise extrapolation and interpolation bias; keep the donor pool clean and remove donors contaminated by spillover; with a single treated unit only permutation inference is available, and its resolution is capped by the number of units; the weights are not unique, so do not read them as "the true control group".
- Sensitivity / partial ID (Chapter 15): coefficient stability puts no constraint on unobserved confounding; the robustness value has no universal size threshold and must be compared against observed benchmarks; a negative-control failure cannot be buried in the appendix; the identified region and the confidence interval express assumption uncertainty and sampling uncertainty respectively.
- DML / causal forest (Chapter 16): DML is not a master key, all its guarantees rest on the product-rate condition, and when the sample is too small or the nuisance too complex the point estimate is biased and the interval undercovers; when unconfoundedness and overlap fail, ML cannot save you either; ML will fit heterogeneity that does not exist, so validate through calibration, GATES, and RATE; the regret guarantee of policy learning is only relative to a restricted policy class.
- Digital experiments (Chapter 17): peeking (continuous monitoring and stopping the moment it is significant) poisons everything downstream, so use always-valid or group-sequential boundaries; interference is subtle (it does not violate randomization and the p-value looks fine, yet it precisely estimates a quantity you did not want), and the fix is to change the design (cluster, switchback) rather than to test; use an SRM chi-square as the first guardrail; handle multiple comparisons with Bonferroni or FDR and pre-register; ex-post observed power carries no information.
- Interference / marketplace experiments (Chapter 18): an individual A/B test may correctly estimate the direct effect yet cannot answer the rollout effect; cluster leakage, carryover, and the wrong level of inference are the three big threats; the short-run experimental effect does not automatically equal the long-run equilibrium effect.

### Part IV Structural Methods

- Discrete choice (Chapter 19): logit coefficients are identified only up to a scale normalization and cannot be compared in magnitude across models or samples, what is comparable are ratios and elasticities; misspecified IIA does not trigger an automatic error, and under separation a finite MLE may not even exist; getting the nests wrong in nested logit or forgetting to instrument the within-nest share pushes $\lambda$ toward 1; the variance terms of mixed logit are identified by rich choice-set variation, which a panel can strengthen substantially but is not necessary; the second-step standard errors of a control function must account for the first step; ordinary SML retains simulation bias at a fixed number of draws $R$.
- BLP demand (Chapter 20): the validity of the price instrument is the foundation (a cost shifter must enter only cost, and BLP instruments degenerate into weak instruments in large markets); random coefficients collapsing toward zero quietly pulls the model back to IIA, so report the standard errors of $\sigma$ and the instrument strength; too loose a contraction tolerance seeps into the outer GMM (MPEC sidesteps this problem); zero shares leave the log share undefined, and discarding them casually introduces selection bias.
- Supply side / mergers (Chapter 21): conduct cannot be identified from the first-order conditions themselves, two conduct assumptions that both yield positive costs cannot be told apart from cost alone, so report cost plausibility plus a Rivers-Vuong test and do not default to Bertrand; a misspecified demand (logit IIA understating near-neighbor diversion) propagates through the FOC and is amplified systematically into markups and merger price increases; an efficiency defense passes through to price only if marginal cost falls; merger simulation is static and short-run, and its extrapolation is limited by the Lucas critique.
- Entry / partial identification (Chapter 22): multiple equilibria break point identification, and forcing in an equilibrium-selection mechanism makes the estimate drift with the assumption, so be honest about the dependence or switch to partial identification; the information structure (complete versus incomplete information) must be argued from institutional detail and is not a technical option to switch at will; the identified set can be so wide as to carry no information, yet the width itself is also information.
- Dynamic structural (Chapter 23): the discount factor $\beta$ is generally not identified, and treating fixing it as harmless is a common overconfidence, so report how the conclusion changes at several values of $\beta$ or identify it formally with an exclusion restriction; conditional independence rules out persistent unobserved heterogeneity, and a violation biases the estimate, so model types explicitly (finite mixtures); counterfactual validity is the Lucas critique; the CCP first step is noisy on rarely visited states, so smooth or gather more data.

### Part V Topics in Digital Economics

- Platforms / switching costs (Chapter 24): the bias in switching costs is determinedly upward, almost never understating as long as persistent heterogeneity is left uncontrolled, so when you see a large switching cost first ask "how much preference has been mistaken for lock-in"; the initial-conditions problem is especially severe in a short panel, so when the panel is not long default to Wooldridge CRE; peer and network effects have the reflection problem, which must be broken with network structure.
- Search / advertising (Chapter 25): a full-information model attenuates systematically toward zero in a search market and a large sample cannot save it, so when the search logs show shallow search it is time to change the model; getting the sequential versus fixed-sample protocol wrong makes the likelihood wrong, so estimate both and test with the prevalence of return visits; the exogeneity of the ranking cannot be taken for granted, and needs a ranking experiment or the quasi-random perturbation of an algorithm change; advertising effects are the easiest to be exaggerated by selection, and the observational-correlation ROI should be discounted heavily.
- AI / technology adoption / algorithms / data (Chapter 26): machine-prediction error is often non-classical, and the gold subsample must represent the target population; AI access, actual use, and workflow redesign are different treatment margins; an observed-use interaction usually suffers post-treatment selection, and organizational complementarity is best studied with factorial assignment; short-run productivity does not equal the long-run skill and J-curve payoff; the appearance of algorithmic collusion in simulation is not collusion in reality, and needs field evidence to pin down.

Strung together, these pitfalls point to a single lesson that runs through the whole book: the credibility of a method lies not in the sophistication of the regression specification or the algorithm, but in whether the identifying assumptions are laid on the table and in how the conclusion moves when those assumptions are loosened. A reference sheet can help you recall the formulas, the commands, and the common traps, but stating the assumptions clearly and reporting sensitivity honestly always sends you back to the main text and to your own judgment about the institutional background.
