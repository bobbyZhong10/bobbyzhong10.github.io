---
title: "Quick Reference & Tools"
subtitle: "Formula Cards, a Method Decision Tree, and a Cross-Software Command Map"
seriesline: "Foundations of Information Systems Economics · Chapter 35"
series: "Foundations of Information Systems Economics: Empirical Methods for Information Systems Economics"
sidebartitle: "Foundations of Information Systems Economics"
sidebarsub: "Theory · Methods · Applications"
chapterlabel: "Chapter 35 · Quick Reference & Tools"
---

## Introduction

The moment you actually need a reference sheet is usually not before an exam, but after the regression has already run. Should the standard errors be clustered by user, by store, or by market? Can an event study with staggered treatment still be run as plain TWFE? When the instrument is weak, which interval is still trustworthy? The results on the screen usually arrive faster than these questions, and they readily tempt you to write down a conclusion on the spot. This chapter sits on your desk for one reason: to make you pause once more at that moment.

Here you look up methods by research question first, not commands by software menu. The formula cards help you check the estimand and the notation; the decision tree starts from your goal, your data structure, and where the identifying variation comes from, and points you to the relevant chapter; the R, Python, and Stata command tables answer "how do I run this"; and the pitfalls checklist presses the harder question: "why might this implementation still be wrong?" Every card and every exit of the tree names its chapter and section, so that you can go straight back to the source. For a marketplace experiment, for instance, the path does not stop at a difference in means: it reminds you of interference, switchback designs, and the rollout estimand. For AI adoption, it will not let you treat access, use, and workflow redesign as one and the same treatment.

The danger of a reference sheet comes precisely from its convenience. A formula has no room for the full institutional background, a command will not check the exclusion restriction for you, and a decision tree cannot decide on your behalf whether parallel trends is credible. So what this chapter offers is an entry point and a set of signposts for the way back, not a verdict on which method to use: once you have a candidate, return to the corresponding chapter and check the economic model, the identifying assumptions, the estimation details, and the failure modes.

The prose of this chapter uses a uniform notation: individuals are $i$, products $j$, markets or time $t$, firms $f$, cohorts $g$; potential outcomes are written $Y_i(d)$, treatment is $D_i$, the instrument is $Z_i$, covariates are $X_i$, and the propensity score is $e(x)$. If this chapter does its job, it will not help you pick a fashionable method faster; it will help you discover faster that the question in front of you has not yet been clearly defined.

## 1 Formula Cards

Each card gives one core method: its name, the key formula, and a one-line gloss; the "Chapter N §x.y" at the end of the label is the source. Formulas are in English and LaTeX throughout, and the cards are grouped by the book's five Parts. A card only restates the theorems and results of the section it cites, and assumption labels (such as M.1 or A5) keep that chapter's own numbering; for the conditions and limits under which a formula holds, go back to the section named on the card. The card formulas keep the notation of the chapter they cite (the Robinson card's $z_i$ and $x_i$, for example) rather than being forced into the uniform notation above.

### Part I Microeconomic Foundations

Part I consists of theory chapters: the main text works through its hand calculations and simulation demonstrations in base R, with no dedicated estimation package. The cards below give the primitives that the later chapters on identification and structural estimation borrow again and again.

::: {.theorem data-label="Consumer duality (Roy / Shephard / the four duality identities) · Chapter 1 §3.3 to §3.4"}
$$x_\ell(p, w) = -\frac{\partial v(p, w)/\partial p_\ell}{\partial v(p, w)/\partial w}, \qquad \frac{\partial e(p, \bar u)}{\partial p_\ell} = h_\ell(p, \bar u)$$

$$e(p, v(p, w)) = w, \qquad h_\ell(p, \bar u) = x_\ell(p, e(p, \bar u))$$

In words: when $u$ is continuous and C.3 holds (Theorem 4), $v$ and $e$ are inverses of each other at fixed prices; add C.4(b), and the price slopes of the two value functions are Walrasian demand $x$ and Hicksian demand $h$ (Theorem 5, a local result).
:::

::: {.theorem data-label="Slutsky equation · Chapter 1 §3.5"}
$$\underbrace{\frac{\partial x_\ell(p, w)}{\partial p_k}}_{\text{total}} = \underbrace{\frac{\partial h_\ell(p, \bar u)}{\partial p_k}}_{\text{substitution}}\; \underbrace{-\; x_k(p, w)\,\frac{\partial x_\ell(p, w)}{\partial w}}_{\text{income}}$$

In words: with $u$ continuous and C.3 and C.4(b) in force (Theorem 6), the price slope of demand splits into a substitution term and an income term; if $h$ is continuously differentiable, the observable Slutsky matrix $S(p, w)$ equals the Hessian of $e$, so it is symmetric and negative semidefinite, with $S(p, w)\,p = 0$.
:::

::: {.theorem data-label="Welfare measures CV / EV · Chapter 1 §4.1"}
$$\mathrm{CV} = e(p^1, u^0) - e(p^1, u^1) = e(p^1, u^0) - w, \qquad \mathrm{EV} = e(p^0, u^0) - e(p^0, u^1) = w - e(p^0, u^1)$$

In words: CV, EV, and $\Delta\mathrm{CS}$ (the loss in consumer surplus) are all signed as losses, the opposite of MWG's sign; CV is priced at the new prices and EV at the old ones. With $u$ continuous, C.3 and C.4(b) holding, and only $p_1$ changing (Theorem 10), CV and EV are areas under Hicksian demand; for a normal good, $\mathrm{EV} \leq \Delta\mathrm{CS} \leq \mathrm{CV}$ for price increases and decreases alike, the ordering reverses for an inferior good, and with no income effect the three are equal.
:::

::: {.theorem data-label="Arrow-Pratt risk aversion · Chapter 1 §4.3"}
$$A(w) = -\frac{u''(w)}{u'(w)}, \qquad \pi(w) \approx \tfrac{1}{2}A(w)\operatorname{Var}(\widetilde\varepsilon)$$

In words: within expected utility (U.1, U.2, Theorem 15), requiring at every $x$ that $A(x; u_2) \geq A(x; u_1)$ is equivalent to requiring that, for any distribution of payoffs, the certainty equivalent under $u_2$ be no higher than under $u_1$ (Theorem 16); for a small zero-mean risk $\widetilde\varepsilon$, the premium is roughly half of $A(w)$ times the variance.
:::

::: {.theorem data-label="Nash equilibrium · Chapter 2 §3.2"}
$$u_i(s_i^*, s_{-i}^*) \geq u_i(s_i, s_{-i}^*) \quad \text{for all } s_i,\ \text{all } i$$

In words: every player's strategy is a best response to the others', so no one has an incentive to deviate unilaterally; in a finite game (G.1(a)), an equilibrium is guaranteed to exist once mixed strategies are allowed (Theorem 2), and a continuous game under G.2 has a pure-strategy equilibrium (Chapter 2 §3.3, Theorem 3).
:::

::: {.theorem data-label="Repeated-game cooperation threshold (grim trigger) · Chapter 2 §4.1"}
$$\frac{R}{1-\delta} \geq T + \frac{\delta P}{1-\delta} \quad\Longleftrightarrow\quad \delta \geq \frac{T - R}{T - P}$$

In words: when $T > R > P$, the punishment profile is a Nash equilibrium of the stage game, and play is repeated infinitely often (G.3), grim trigger is an SPE if and only if the discounted value of cooperating forever is at least that of defecting once and then being punished forever (Theorem 16), so the cooperation threshold is set by the ratio of the temptation to the severity of punishment; the folk theorem (G.1(a), G.3, Theorem 17) says that, in two-player games or in games with more players that satisfy a dimensionality condition, any feasible payoff profile strictly above the minmax is the average payoff of some SPE when $\delta$ is large enough.
:::

::: {.theorem data-label="VCG mechanism · Chapter 3 §3.6"}
$$t_i = \max_{a}\sum_{j \neq i} v_j'(a) \;-\; \sum_{j \neq i} v_j'(a^*)$$

In words: with quasilinear private values (Theorem 8), each agent pays the externality it imposes on everyone else, so private self-interest and total social welfare differ only by a constant that is independent of the agent's own report, and truthful reporting becomes a (weakly) dominant strategy, at the price of a budget that is in general not balanced (the revelation principle guarantees that we need only study direct mechanisms with truthful reporting, Chapter 3 §3.5, Theorem 6).
:::

::: {.theorem data-label="Revenue equivalence and the Myerson reserve · Chapter 3 §4.1 to §4.3"}
$$b(v) = \frac{n-1}{n}\, v, \qquad \Psi(r^\ast) = r^\ast - \frac{1 - F(r^\ast)}{f(r^\ast)} = 0$$

In words: with risk-neutral bidders and independent private values (I.1, I.2), if the Bayesian Nash equilibria of two auctions give each bidder the same probability of winning at every profile of values and give the lowest type the same expected utility, the two earn the same expected revenue (Theorem 12); when the seller values the object at zero, the virtual value $\Psi$ is increasing, and bidders are symmetric (I.1'), a second-price or first-price auction with reserve price $r^\ast$ maximizes revenue, and $r^\ast$ does not depend on the number of bidders (Theorem 13). With values uniform on $[0, 1]$, the left formula is the symmetric equilibrium bid in a first-price auction (Theorem 11), and $r^\ast = 1/2$.
:::

::: {.theorem data-label="Lerner index / monopoly markup · Chapter 4 §2"}
$$\underbrace{p(q) + q\,\frac{dp}{dq}}_{MR} = \underbrace{c'(q)}_{MC}, \qquad \frac{p - c'}{p} = \frac{1}{\lvert\mathcal{E}\rvert}$$

In words: under D.1 (Theorem 1), the monopolist sets marginal revenue equal to marginal cost, and at the monopoly price the Lerner index equals the reciprocal of the absolute price elasticity of demand (not a causal relation); with positive marginal cost, demand is elastic there; this inverse-elasticity relation is the starting point for recovering markups from demand in Chapter 29.
:::

::: {.theorem data-label="Cournot and Bertrand · Chapter 4 §4.1 to §4.2"}
$$\text{Cournot}:\ p'(Q^\ast)\,Q^\ast/n + p(Q^\ast) = c,\ \ p(Q^\ast) > c \qquad\qquad \text{Bertrand}:\ p_1 = p_2 = c$$

In words: the dimension of competition determines market power. With a common marginal cost $c$ (D.2), Cournot competition in quantities has a pure-strategy equilibrium under D.1 and O.1, and the symmetric equilibrium price is above $c$ and, as $n \to \infty$, tends to $c$ (Theorem 7); in Bertrand competition in prices with a homogeneous product and $c > 0$, the unique equilibrium is $p_1 = p_2 = c$ (Theorem 8, the Bertrand paradox). Capacity constraints, product differentiation, or repeated play can push price back up.
:::

### Part II Econometric Core

Part II is the machinery of estimation and inference: every downstream reduced-form and structural method draws on the M-estimation, sandwich variances, clustering, GMM, panel methods, and nonparametric smoothing collected here. Part II always writes estimation as minimization: an objective to be maximized (such as a log-likelihood) enters with a minus sign.

::: {.theorem data-label="The M-estimation framework · Chapter 5 §3.1"}
$$\hat\theta = \arg\min_{\theta \in \Theta} Q_n(\theta), \qquad Q_n(\theta) = \frac{1}{n}\sum_{i=1}^n q(w_i,\theta), \qquad q_{MLE}(w,\theta) = -\ln f(w;\theta)$$

$$\hat\theta_{GMM} = \arg\min_\theta\ \bar g_n(\theta)'\, W\, \bar g_n(\theta), \qquad \bar g_n(\theta) = \frac{1}{n}\sum_{i=1}^n g(w_i, \theta)$$

In words: OLS, NLS, Poisson PML, and MLE differ only in the per-observation loss $q$, while GMM is a quadratic form in sample moments, and all of them are extremum estimators; MLE and GMM are developed in Chapters 6 and 7. Identification requires the population objective $Q_0(\theta) = \mathbb{E}[q(w,\theta)]$ to have a unique minimum (M.1, Chapter 5 §3.2); under misspecification $\hat\theta$ converges to the pseudo-true value picked by $q$ itself, $\theta_*$, not to a structural parameter.
:::

::: {.theorem data-label="Asymptotic normality and the sandwich · Chapter 5 §3.4"}
$$\sqrt{n}(\hat\theta - \theta_0) \xrightarrow{d} \mathcal{N}\big(0,\ A_0^{-1} B_0 A_0^{-1}\big), \qquad A_0 = \mathbb{E}[\nabla^2_\theta q(w,\theta_0)], \quad B_0 = \mathbb{E}[s(w,\theta_0)\,s(w,\theta_0)']$$

In words: under M.1 to M.3 (Theorem 3), with $s = \nabla_\theta q'$ the score, the bread is the curvature and the filling is the noise in the score. The variance collapses only when $B_0 = \sigma^2 A_0$, and the model's default standard errors are right only then, so report the sandwich by default.
:::

::: {.theorem data-label="NLS · Chapter 5 §3.5"}
$$\hat\theta = \arg\min_\theta \tfrac{1}{2}\sum_{i=1}^n \big(y_i - m(x_i,\theta)\big)^2, \qquad A_0 = \mathbb{E}[\nabla_\theta m'\nabla_\theta m], \quad B_0 = \mathbb{E}[\varepsilon^2\,\nabla_\theta m'\nabla_\theta m]$$

In words: following the convention of Theorem 4, take $q = (y - m(x,\theta))^2/2$; only then do the $A_0$ and $B_0$ above hold. Under NLS.1 and NLS.2 the estimator is consistent and asymptotically normal, and only when homoskedasticity (NLS.3) is added does $B_0 = \sigma^2 A_0$ hold (Theorem 4); with a misspecified mean it converges to the value that minimizes $\mathbb{E}[(\mathbb{E}[y \mid x] - m(x,\theta))^2]$, namely $\theta_*$, and inference needs the sandwich with the sample Hessian as bread (Theorem 5).
:::

::: {.theorem data-label="Poisson PML · Chapter 5 §3.7"}
$$q(w,\theta) = \exp(x'\theta) - y\,x'\theta, \qquad A_0 = \mathbb{E}[\exp(x'\theta_0)\,xx'], \quad B_0 = \mathbb{E}[\mathrm{Var}(y \mid x)\,xx']$$

In words: as long as $\mathbb{E}[y \mid x] = \exp(x'\theta_0)$ (PML.1, PML.2), the estimator is consistent with no assumption on the variance at all, and $y$ need not be a count (Theorem 7). Use it for non-negative outcomes with many zeros; the coefficients of a log(1+y) regression change with the units of measurement (Chapter 5 §3.8).
:::

::: {.theorem data-label="Two-step M-estimation · Chapter 5 §3.6"}
$$\mathrm{Avar}\big(\sqrt n(\hat\theta - \theta_0)\big) = A_0^{-1}\,\mathrm{Var}\big(s_i(\theta_0,\gamma_*) + \Gamma_0\,\psi_i(\gamma_*)\big)\,A_0^{-1}, \qquad \Gamma_0 = \mathbb{E}[\nabla_\gamma s(w,\theta_0,\gamma_*)]$$

In words: under TS.1 to TS.3, the first step's estimation error $\psi_i$ enters the second step's filling through $\Gamma_0$ (Theorem 6); when $\Gamma_0 = 0$ you may treat $\hat\gamma$ as known; otherwise, use this correction or bootstrap the whole two-step procedure.
:::

::: {.theorem data-label="MLE and the information matrix equality · Chapter 6 §3.3"}
$$A_0 = B_0 \equiv I(\theta_0), \qquad \sqrt{n}(\hat\theta_{MLE} - \theta_0) \xrightarrow{d} \mathcal{N}\big(0,\ I(\theta_0)^{-1}\big)$$

In words: taking $q = -\ln f$, under ML.1 to ML.3 the conditional information matrix equality collapses the sandwich to $I(\theta_0)^{-1}$ (Theorem 3); this is the variance bound for a class of moment estimators, and the CMLE attains it (Chapter 6 §3.4). That efficiency requires the whole conditional distribution to be correctly specified.
:::

::: {.theorem data-label="QMLE · Chapter 6 §3.6"}
$$\theta_* = \arg\min_\theta\ \mathbb{E}\big[\mathrm{KL}\big(f_0(\cdot \mid x) \,\Vert\, f(\cdot \mid x;\theta)\big)\big], \qquad \sqrt n(\hat\theta - \theta_*) \xrightarrow{d} \mathcal{N}(0, A_*^{-1}B_*A_*^{-1})$$

In words: with a misspecified likelihood, the MLE estimates the pseudo-true value closest in the KL sense, and in general $A_* \neq B_*$; the only valid variance is the sandwich with the sample Hessian as bread, the LM test must take its robust form, and LR is no longer asymptotically $\chi^2$ (Chapter 6 Theorem 6). Poisson is the happy case: with a correctly specified conditional mean it remains consistent (Chapter 5 §3.7).
:::

::: {.theorem data-label="The testing trinity Wald / LR / LM · Chapter 6 §3.5"}
$$\mathrm{W} = \frac{(\hat\theta - \theta^{\circ})^2}{\widehat{\mathrm{Var}}(\hat\theta)}, \qquad \mathrm{LR} = 2(\ell_u - \ell_r), \qquad \mathrm{W},\ \mathrm{LR},\ \mathrm{LM} \xrightarrow{d} \chi^2_\nu$$

In words: the W above is the single-restriction form ($\nu = 1$, with null value $\theta^{\circ}$); under ML.1 to ML.3, TR.1, and $H_0$, the three statistics share the same asymptotic distribution and local power (Theorem 5 (i)); with a misspecified likelihood only the sandwich Wald and the robust LM remain valid, and LR fails (Theorem 5 (ii)). Wald is not invariant to reparameterization, so when the parameterization is in doubt use LR (Chapter 8 §4.1).
:::

::: {.theorem data-label="GMM and efficient weighting · Chapter 7 §3.3"}
$$\mathrm{Avar}\big(\sqrt n(\hat\theta - \theta_0)\big) = (D_0'WD_0)^{-1}D_0'WS_0WD_0(D_0'WD_0)^{-1}, \qquad W = S_0^{-1}:\ (D_0'S_0^{-1}D_0)^{-1}$$

In words: the sandwich under GMM.1 to GMM.3 (Theorem 2) minus $(D_0'S_0^{-1}D_0)^{-1}$ is positive semidefinite, with equality at $W = S_0^{-1}$; two-step GMM computes $\hat S$ from a first-step estimate and attains this bound (Theorem 3). $D_0$ is the Jacobian of the moment conditions and $S_0$ their covariance.
:::

::: {.theorem data-label="Hansen J test · Chapter 7 §3.4"}
$$J = n\,\bar g_n(\hat\theta^{(2)})'\,\hat S^{-1}\,\bar g_n(\hat\theta^{(2)}) \xrightarrow{d} \chi^2_{r-k}$$

In words: with more moments than parameters ($r > k$), it tests whether the moment conditions are mutually compatible (Theorem 4); it can only detect contradictions and has no power when all moments are wrong in the same direction. Use the Windmeijer correction for two-step standard errors (Chapter 7 §3.5), and switch to the S confidence set under weak moments (Chapter 7 §3.6).
:::

::: {.theorem data-label="Classical minimum distance · Chapter 7 §3.7"}
$$\hat\theta_{MD} = \arg\min_\theta\ \big(\hat\pi - h(\theta)\big)'\hat\Sigma_\pi^{-1}\big(\hat\pi - h(\theta)\big), \qquad \sqrt n(\hat\theta_{MD} - \theta_0) \xrightarrow{d} \mathcal{N}\big(0, (D_0'\Sigma_\pi^{-1}D_0)^{-1}\big)$$

In words: first estimate the reduced form $\hat\pi$, then find the structural parameters whose $h(\theta)$ comes closest to it; under MD.1 and MD.2 the optimal weight gives the smallest variance, the residual quadratic form is a $\chi^2_{r-k}$ overidentification test (Theorem 7), and the estimator is asymptotically equivalent to the corresponding two-step GMM (Chapter 7 §3.8). $\hat\Sigma_\pi$ must be the full variance of $\hat\pi$.
:::

::: {.theorem data-label="Delta method · Chapter 8 §3.2"}
$$\sqrt{n}\big(g(\hat\theta) - g(\theta_0)\big) \xrightarrow{d} \mathcal{N}\big(0,\ \nabla g(\theta_0)'\, \Sigma\, \nabla g(\theta_0)\big)$$

In words: when $g$ is continuously differentiable at the truth with a nonzero gradient, a first-order Taylor expansion gives standard errors for derived quantities such as ratios, elasticities, and WTP; when the denominator is near zero (the Fieller problem), the approximation breaks down and you should switch to a Fieller interval.
:::

::: {.theorem data-label="Cluster-robust variance and the Moulton factor · Chapter 8 §3.4"}
$$\hat V^{CR} = (X'X)^{-1} \left( \sum_{g=1}^{G} X_g'\, \hat u_g\, \hat u_g'\, X_g \right) (X'X)^{-1}, \qquad \sqrt{1 + (\bar n - 1)\,\rho}$$

In words: adding the within-cluster cross-covariances back into the filling of the sandwich corrects the understatement caused by the iid assumption; the Moulton factor is the multiple by which within-cluster correlation inflates the standard error. The asymptotics are in the number of clusters $G$, not the number of observations; with few clusters, use $t(G-1)$, CR2, or the wild cluster bootstrap (Chapter 8 §4.2).
:::

::: {.theorem data-label="HAC (Newey-West) · Chapter 8 §3.5"}
$$\hat\Omega^{HAC} = \sum_{t} \hat u_t^2\, x_t x_t' + \sum_{\ell=1}^{L} \omega_\ell \sum_{t=\ell+1}^{T} \hat u_t \hat u_{t-\ell}\big(x_t x_{t-\ell}' + x_{t-\ell} x_t'\big),\quad \omega_\ell = 1 - \frac{\ell}{L+1}$$

In words: the time-series counterpart of clustering, which uses Bartlett triangular weights to add several lagged autocovariances to the variance estimate and so handles serial correlation and heteroskedasticity.
:::

::: {.theorem data-label="2SLS and its robust variance · Chapter 9 §3.3"}
$$\hat\beta_{2SLS} = (\hat{\mathbf{X}}'\hat{\mathbf{X}})^{-1}\hat{\mathbf{X}}'\mathbf{y}, \qquad \widehat{\mathrm{Avar}} = (\hat{\mathbf{X}}'\hat{\mathbf{X}})^{-1}\Big(\sum_i\hat\varepsilon_i^2\hat{\mathbf{x}}_i\hat{\mathbf{x}}_i'\Big)(\hat{\mathbf{X}}'\hat{\mathbf{X}})^{-1}$$

In words: under 2SLS.1 and 2SLS.2 it is consistent and asymptotically normal, and the residual $\hat\varepsilon_i = y_i - \mathbf{x}_i'\hat\beta_{2SLS}$ uses the original regressors, not the fitted values (Theorem 3); under homoskedasticity (2SLS.3), it is the most efficient among estimators that use linear combinations of $\mathbf{z}_i$ as instruments. With weak instruments, 2SLS is pulled toward OLS, so switch to the Anderson-Rubin set (Chapter 9 §3.5).
:::

::: {.theorem data-label="SUR and FGLS · Chapter 9 §3.6"}
$$\hat\beta_{FGLS} = \Big(\sum_{i=1}^n\mathbf{X}_i'\hat\Omega^{-1}\mathbf{X}_i\Big)^{-1}\sum_{i=1}^n\mathbf{X}_i'\hat\Omega^{-1}\mathbf{y}_i$$

In words: under SGLS.1 and SGLS.2 FGLS is consistent and asymptotically equivalent to GLS with known $\Omega$, and its variance takes the fully robust form; under the homoskedasticity-type conditions of Theorem 6 (iv), it is more efficient than equation-by-equation OLS. When $\hat\Omega$ is diagonal or every equation has the same regressors, FGLS is equation-by-equation OLS (Theorem 7).
:::

::: {.theorem data-label="System 2SLS, 3SLS, and system GMM · Chapter 9 §3.7"}
$$\hat\beta = \arg\min_\beta\ \Big(\sum_i \mathbf{Z}_i'(\mathbf{y}_i - \mathbf{X}_i\beta)\Big)'\,W\,\Big(\sum_i \mathbf{Z}_i'(\mathbf{y}_i - \mathbf{X}_i\beta)\Big), \qquad W = S_0^{-1}$$

In words: under SIV.1 to SIV.3 system GMM is consistent, the variance is smallest at $W = S_0^{-1}$, and under SIV.5, 3SLS is exactly that efficient GMM (Theorem 8). The price is contagion: when one equation is misspecified, 3SLS is generally inconsistent for all parameters, while 2SLS on a correctly specified equation stays consistent (Chapter 9 §3.8).
:::

::: {.theorem data-label="FE, RE, and FD · Chapter 10 §3.3"}
$$\hat\beta_{FE}:\ \ddot Y_{it}\ \text{on}\ \ddot X_{it}, \qquad \hat\beta_{RE}:\ Y_{it} - \hat\varphi\bar Y_i\ \text{on}\ (1 - \hat\varphi,\ X_{it} - \hat\varphi\bar X_i), \qquad \hat\beta_{FD}:\ \Delta Y_{it}\ \text{on}\ \Delta X_{it}$$

In words: FE is consistent under FE.1 (strict exogeneity) and FE.2 without restricting the correlation between $\alpha_i$ and the regressors, and its standard errors are clustered by unit (Theorem 3); RE additionally requires $\mathbb{E}[\alpha_i \mid \mathbf{X}_i] = \mathbb{E}[\alpha_i]$ (RE.1, Chapter 10 §3.2). When only contemporaneous exogeneity holds and the series are stationary and weakly dependent (no unit root in the idiosyncratic error), the inconsistency of FE is $O(1/T)$, while that of FD does not vanish as $T$ grows (Theorem 4, Chapter 10 §3.4).
:::

::: {.theorem data-label="Mundlak (CRE) and the robust Hausman test · Chapter 10 §3.6"}
$$Y_{it}\ \text{on}\ (1, X_{it}, \bar X_i):\ \text{coefficient on } X_{it} = \hat\beta_{FE}, \qquad H_0:\ \xi = 0$$

In words: in pooled OLS or RE on $(1, X_{it}, \bar X_i)$, the coefficient on $X_{it}$ numerically equals the FE estimate (Theorem 5, Chapter 10 §3.5); the cluster-robust Wald test that the coefficient on $\bar X_i$ satisfies $\xi = 0$ is $\chi^2_k$ under RE.1 and FE.2 whether or not RE.3 holds; this is the robust Hausman test (Theorem 7 (ii)). The classical version needs all of RE.1 to RE.3 and FE.2.
:::

::: {.theorem data-label="Nickell bias · Chapter 11 §3.2"}
$$T \cdot \mathrm{plim}_{N \to \infty}(\hat\rho_{FE} - \rho) \to -(1+\rho) \qquad (T \to \infty)$$

In words: under NK.1, within-demeaning mixes the error into the lagged term, so the FE estimate of persistence has a downward bias of order $1/T$ that does not vanish as $N$ grows (Theorem 2); through correlation, the coefficients on the other covariates in the model are biased as well. For bias-corrected FE see Chapter 11 §3.8.
:::

::: {.theorem data-label="Arellano-Bond and Blundell-Bond · Chapter 11 §3.5"}
$$\text{difference GMM}:\ \mathbb{E}[Y_{i,t-s}\,\Delta\varepsilon_{it}] = 0,\ s \geq 2 \qquad \text{system GMM}:\ \mathbb{E}[\Delta Y_{i,t-1}(\alpha_i + \varepsilon_{it})] = 0$$

In words: first-difference to remove $\alpha_i$, then use deep lagged levels as instruments (DC.1, AB.1, Chapter 11 §3.3); when $\rho$ is close to 1, the instruments are weak, and system GMM adds moment conditions in levels, which rest on BB.1 (the initial deviation from the long-run level is uncorrelated with $\alpha_i$); when BB.1 fails and $\rho \neq 0$, system GMM is generally inconsistent, while difference GMM is unaffected (Theorem 5).
:::

::: {.theorem data-label="The Windmeijer correction and instrument proliferation · Chapter 11 §3.4"}
$$\hat V_c = \hat V_2 + \hat D_W\hat V_2 + \hat V_2\hat D_W' + \hat D_W\hat V_1\hat D_W', \qquad J = N\,\bar g_N(\hat\rho^{(2)}_{AB})'\,\hat S^{-1}\,\bar g_N(\hat\rho^{(2)}_{AB})$$

In words: the usual two-step variance treats $\hat S$ as known and is too small in finite samples, and the Windmeijer correction adds the term that arises because $\hat S$ depends on the one-step estimate (Theorem 4; the general formula is in Chapter 7 §3.5); $m_2$ should be insignificant and $m_1$ significantly negative. With $r \geq N$ instruments, a one-step fit gives $J = N$ regardless of the data, and the two-step $J$ is at most $N$ (Chapter 11 §3.6), so collapse the instruments and report their number.
:::

::: {.theorem data-label="APE and interaction terms · Chapter 12 §3.3"}
$$\widehat{\mathrm{APE}}_j = \frac{1}{n}\sum_{i=1}^n F'(X_i'\hat\beta)\,\hat\beta_j, \qquad \frac{\partial^2 P(x)}{\partial x_j\,\partial x_{j'}} = \beta_{jj'}F'(x'\beta) + (\beta_j + \beta_{jj'}x_{j'})(\beta_{j'} + \beta_{jj'}x_j)\,F''(x'\beta)$$

In words: probit and logit coefficients are meaningful only up to a scale normalization; what is comparable is the APE (consistent under BR.1 and BR.2, Theorem 3); the cross-partial effect is generally nonzero even when $\beta_{jj'} = 0$, its sign can change with $x$, and its significance cannot be read off the $\beta_{jj'}$ coefficient's $z$ statistic.
:::

::: {.theorem data-label="Control function and IV probit · Chapter 12 §3.5"}
$$\Pr(Y_i = 1 \mid X_i, Z_i, v_i) = \Phi\big((X_i'\beta + \beta_D D_i + \rho v_i)/\sigma_e\big), \qquad \mathrm{ASF}(x, d) = \Phi(x'\beta + \beta_D d)$$

In words: under CF.1 and CF.2, first run the first-stage OLS of a continuous $D_i$ and keep the residual $\hat v_i$, then add it to the probit: this consistently estimates the scaled coefficients, and the $\hat v_i$ term's $t$ statistic tests exogeneity; the APE comes from averaging the ASF over $\hat v_i$ and differentiating (Theorem 5).
:::

::: {.theorem data-label="FE logit and incidental parameters · Chapter 12 §3.8"}
$$\Pr(Y_{i2} = 1 \mid Y_{i1} + Y_{i2} = 1, X_{i1}, X_{i2}, \alpha_i) = \Lambda\big((X_{i2} - X_{i1})'\beta\big)$$

In words: conditioning on $\sum_t Y_{it}$ removes $\alpha_i$, and the conditional MLE is consistent for fixed $T$, but the APE is not identified (Theorem 8); the dummy-variable MLE that treats $\alpha_i$ as parameters is inconsistent for fixed $T$, and when $N$ and $T$ grow at the same rate it has a bias of order $1/T$ that can be corrected analytically (Theorem 9). When only the APE is needed, CRE probit works (Chapter 12 §3.7).
:::

::: {.theorem data-label="The check function and QR asymptotics · Chapter 13 §3.3"}
$$\rho_\tau(u) = \big(\tau - \mathbf{1}\{u < 0\}\big)\,u, \qquad \sqrt{n}\big(\hat\beta(\tau) - \beta(\tau)\big) \xrightarrow{d} \mathcal{N}\big(0,\ A_0(\tau)^{-1}B_0(\tau)A_0(\tau)^{-1}\big)$$

$$A_0(\tau) = \mathbb{E}\big[f_{Y \mid X}(X_i'\beta(\tau) \mid X_i)\,X_iX_i'\big], \qquad B_0(\tau) = \tau(1 - \tau)\,\mathbb{E}[X_iX_i']$$

In words: QR is M-estimation with the check function as loss and is asymptotically normal under QR.1 to QR.4; the bread $A_0(\tau)$ is the density-weighted Jacobian of the expected score, not $\mathbb{E}[\nabla^2 q]$ (the check function's second derivative is zero almost everywhere), and the density is low in the tails, so tail quantiles are estimated imprecisely (Theorem 3); under misspecification the QR coefficient is the density-weighted best linear approximation (Chapter 13 §3.4).
:::

::: {.theorem data-label="QTE and rank invariance · Chapter 13 §3.5"}
$$\mathrm{QTE}(\tau) = F_{Y(1)}^{-1}(\tau) - F_{Y(0)}^{-1}(\tau)$$

In words: under random assignment (QTE.1) it is identified by the difference of the two groups' marginal quantiles, and under selection on observables (QTE.2) by propensity-score weighting; it is a shift of the distribution, and reading it as the individual effect for the units ranked at $\tau$ also requires rank invariance (QTE.3), without which only sign results remain (Theorem 5).
:::

::: {.theorem data-label="RIF regression and unconditional quantile effects · Chapter 13 §3.6"}
$$\mathrm{RIF}(Y;\tau) = F_Y^{-1}(\tau) + \frac{\tau - \mathbf{1}\{Y \leq F_Y^{-1}(\tau)\}}{f_Y(F_Y^{-1}(\tau))}$$

In words: the expectation of the RIF is exactly the quantile; under RIF.1 and RIF.2 the UQPE is the difference in mean RIF between the two groups of the binary policy variable $W_i$ (in Chapter 13, whether a merchant uses FBM), and when $W_i$ is independent of $X_i$, the OLS coefficient of the RIF on $(1, W_i, X_i)$ is the UQPE (Theorem 6). The conditional QR coefficient answers a different question.
:::

::: {.theorem data-label="Kernel regression (Nadaraya-Watson / local linear) · Chapter 14 §3.1 to §3.3"}
$$\hat m(x_0) = \frac{\sum_{i} y_i\, K\!\big(\frac{x_0 - x_i}{h}\big)}{\sum_{i} K\!\big(\frac{x_0 - x_i}{h}\big)}, \qquad \min_{a, b}\ \sum_{i} K\!\Big(\frac{x_0 - x_i}{h}\Big)\,\big(y_i - a - b\,(x_i - x_0)\big)^2$$

In words: a distance-decaying weighted local average, or local weighted least squares, estimates the whole curve $m(x) = \mathbb{E}[y \mid x]$, with the bandwidth $h$ as the smoothing knob; local linear corrects the boundary bias of Nadaraya-Watson.
:::

::: {.theorem data-label="Bias-variance and the optimal bandwidth · Chapter 14 §3.2"}
$$\text{Bias} \approx \frac{h^2}{2}\, m''(x_0)\!\int\! u^2 K(u)\, du, \quad \text{Var} \approx \frac{1}{nh}\, \frac{\sigma^2}{f(x_0)}\!\int\! K(u)^2 du, \quad h^\ast \propto n^{-1/5}$$

In words: the bias grows with $h^2$ and the variance is of order $1/(nh)$, so the MSE is U-shaped, and the optimal one-dimensional rate $n^{-4/5}$ is slower than the parametric $n^{-1}$; at the optimal bandwidth, the bias and the standard deviation are of the same order, so you must undersmooth or correct the bias explicitly for the confidence interval to cover correctly (Chapter 14 §4.1 to §4.2).
:::

::: {.theorem data-label="Robinson partially linear (double residual) · Chapter 14 §4.3"}
$$y_i - \mathbb{E}[y_i \mid x_i] = \beta\,\big(z_i - \mathbb{E}[z_i \mid x_i]\big) + \varepsilon_i, \qquad \sqrt{n}\,(\hat\beta - \beta) \xrightarrow{d} \mathcal{N}(0, V)$$

In words: partial out the nuisance function $g(x)$ nonparametrically, then regress one set of residuals on the other to estimate the $\beta$ of interest, so that the slow convergence of $g$ does not drag down $\beta$'s $\sqrt{n}$ rate. Here $V = \mathbb{E}[\varepsilon^2 e_z^2]/(\mathbb{E}[e_z^2])^2$ and $e_z = z - \mathbb{E}[z \mid x]$. This double-residual regression is precisely the forerunner of DML in Chapter 24.
:::

### Part III Causal Identification Toolbox

Part III gives one family of estimands per module and is the workhorse of reduced-form causal research. The identifying assumptions all concern counterfactuals and are mostly untestable, so the cards give only the expressions that state the results.

::: {.theorem data-label="Selection-bias decomposition · Chapter 15 §3.1"}
$$\underbrace{\mathbb{E}[Y_i(1) \mid D_i = 1] - \mathbb{E}[Y_i(0) \mid D_i = 0]}_{\text{observed difference}} = \underbrace{\mathbb{E}[Y_i(1) \mid D_i = 1] - \mathbb{E}[Y_i(0) \mid D_i = 1]}_{\text{ATT}} + \underbrace{\mathbb{E}[Y_i(0) \mid D_i = 1] - \mathbb{E}[Y_i(0) \mid D_i = 0]}_{\text{selection bias}}$$

In words: the naive difference in group means equals the ATT plus a selection-bias term. The latter is the difference the two groups would show in the untreated state anyway; it can be nonzero even when treatment effects are homogeneous, and randomization exists precisely to set it to zero.
:::

::: {.theorem data-label="Strong ignorability (unconfoundedness + overlap) · Chapter 15 §3.4 to §3.5"}
$$\{Y_i(1), Y_i(0)\} \perp D_i \mid X_i, \qquad 0 < e(X_i) < 1, \quad e(X) \equiv \Pr(D = 1 \mid X)$$

In words: after controlling for $X$, treatment is independent of the potential outcomes (everything that affects both adoption and outcome is in $X$), and every type of unit has a positive probability of being in either treatment state; the former is untestable, the latter can be checked.
:::

::: {.theorem data-label="IPW and AIPW (doubly robust) · Chapter 15 §4.4 to §4.5"}
$$\hat\tau_{HT} = \frac{1}{N} \sum_i \frac{D_i Y_i}{\hat e(X_i)} - \frac{1}{N} \sum_i \frac{(1 - D_i) Y_i}{1 - \hat e(X_i)}$$

$$\hat\tau_{AIPW} = \frac{1}{N} \sum_{i} \left\{ \hat\mu_1(X_i) - \hat\mu_0(X_i) + \frac{D_i\,(Y_i - \hat\mu_1(X_i))}{\hat e(X_i)} - \frac{(1 - D_i)\,(Y_i - \hat\mu_0(X_i))}{1 - \hat e(X_i)} \right\}$$

In words: IPW reconstructs the two counterfactual means by weighting with the inverse of the probability of being observed; AIPW adds an outcome model $\hat\mu_d$ on top and is consistent as long as either the propensity-score model or the outcome model is correctly specified (doubly robust).
:::

::: {.theorem data-label="Research-design checklist · Chapter 16 §2.1, §3.3, §4.3"}
$$\mathcal D=\{\text{unit},\text{treatment},\text{version},\text{timing},\text{exposure},\text{outcome},\text{population},\text{assignment}\}$$

In words: before running a regression, fix the eight design objects, and in particular distinguish assignment, actual adoption, and algorithmic exposure; when adding a control, recheck whether the estimand has changed and whether the variable is realized after treatment.
:::

::: {.theorem data-label="Prediction-powered rectification · Chapter 17 §3.3"}
$$\widehat\theta_{PPI}=\widehat\theta_{pred,all}+\big(\widehat\theta_{gold,labeled}-\widehat\theta_{pred,labeled}\big)$$

In words: a large volume of machine predictions supplies precision when they are informative (with nearly uninformative predictions, basic PPI is even wider than gold-only, and PPI++ is needed), and a random gold sample estimates the bias between prediction and truth in the target estimating equation; the expression above is the simplified form for linear regression, and in general the rectifier is built inside the same estimating equation. Validity comes from the validation design, not from treating model accuracy as a guarantee of unbiasedness.
:::

::: {.theorem data-label="IV assumptions and Wald / LATE · Chapter 18 §3.4"}
$$\beta_{IV} = \frac{\mathbb{E}[Y_i \mid Z_i = 1] - \mathbb{E}[Y_i \mid Z_i = 0]}{\mathbb{E}[D_i \mid Z_i = 1] - \mathbb{E}[D_i \mid Z_i = 0]} = \mathbb{E}\big[Y_i(1) - Y_i(0) \mid \text{complier}\big]$$

In words: under the four conditions of relevance, random assignment, exclusion, and monotonicity, the Wald ratio (reduced form divided by first stage) identifies the LATE for the complier subpopulation, and only relevance is directly testable.
:::

::: {.theorem data-label="MTE (marginal treatment effect) · Chapter 18 §3.6"}
$$D_i = \mathbf{1}\{U_{D,i} \leq P(Z_i)\}, \qquad \Delta^{MTE}(u) = \mathbb{E}\big[Y_i(1) - Y_i(0) \mid U_{D,i} = u\big]$$

$$\text{ATE} = \int_0^1 \Delta^{MTE}(u)\, du, \qquad \text{ATT} = \int_0^1 \Delta^{MTE}(u)\, \omega_{ATT}(u)\, du$$

In words: the MTE is the effect for the marginal population whose "resistance" $U_D = u$ sits exactly at the threshold; ATE, ATT, LATE, and OLS are all integrals of it with different weights, and $P(Z)$ needs enough support to trace out the whole curve.
:::

::: {.theorem data-label="2×2 DiD · Chapter 19 §3.1"}
$$\Delta_{DD} = \big(\mathbb{E}[Y_{i2} \mid D=1] - \mathbb{E}[Y_{i1} \mid D=1]\big) - \big(\mathbb{E}[Y_{i2} \mid D=0] - \mathbb{E}[Y_{i1} \mid D=0]\big) = ATT$$

In words: the before-after change in the treated group minus the contemporaneous change in the control group identifies the ATT under parallel trends (A1), no anticipation (A2), and SUTVA (A3) (Theorem 1), and each of the three assumptions blocks one source of bias.
:::

::: {.theorem data-label="Goodman-Bacon decomposition · Chapter 19 §4.3"}
$$\mathrm{plim}\; \hat\beta_{TWFE} = VWATT + VWCT - \Delta ATT$$

In words: in a staggered design, the static TWFE estimate is exactly a weighted average of all the 2×2 comparisons, with non-negative weights summing to one (Theorem 4); the $\Delta ATT$ term comes from the forbidden comparison that uses already-treated units as controls, and when effects grow over time it systematically pulls the estimate down or even flips its sign.
:::

::: {.theorem data-label="Group-time ATT and CS / SA / BJS · Chapter 19 §3.6"}
$$ATT(g,t) = \mathbb{E}\big[Y_t - Y_{g-\kappa-1} \mid G = g\big] - \mathbb{E}\big[Y_t - Y_{g-\kappa-1} \mid C = 1\big]$$

In words: under A2', A3, A4, and the never-treated version of parallel trends (A5), the group-time effect is identified, and under A5' the control group becomes the not-yet-treated units (Theorem 2); $\kappa$ is the number of anticipation periods allowed. CS estimates cell by cell and then aggregates with transparent weights, Sun-Abraham uses a saturated interaction regression, and BJS fits the fixed effects on untreated observations only and then imputes the counterfactual, with an implicit parallel trends assumption that covers the whole pre-period (Chapter 19 §4.5).
:::

::: {.theorem data-label="HonestDiD (Rambachan-Roth) · Chapter 19 §6.9"}
$$\big|(b_{e+1} - b_e) - (b_e - b_{e-1})\big| \leq M, \qquad \big|b_{e+1} - b_e\big| \leq \bar M \max_{s \leq -2} \big|b_{s+1} - b_s\big|$$

In words: instead of assuming that parallel trends holds exactly, it restricts the counterfactual trend difference $b_e$ to a smoothness (SD) or relative-magnitudes (RM) set and gives robust intervals that are uniformly valid over that set; the $M$ or $\bar M$ at which the interval first includes zero is the breakdown value. Report both, and state the yardstick each uses (Chapter 19 §7).
:::

::: {.theorem data-label="Fuzzy DiD and Wald-DID · Chapter 20 §3.1"}
$$W_{DID} = \frac{DID_Y}{DID_D}, \qquad DID_R = \mathbb{E}(R_{11}) - \mathbb{E}(R_{10}) - \big[\mathbb{E}(R_{01}) - \mathbb{E}(R_{00})\big]$$

In words: under A1 to A4, $W_{DID}$ is a weighted combination of the switchers' effects in the treatment and control groups, and it equals $\mathrm{LATE}_{sw}$ only when A5 (equal switcher effects in the two groups) or A6 (a stable control-group treatment rate) is added (Theorem 1). A4 requires the effect for the incumbent treated to be stable over time and is often violated.
:::

::: {.theorem data-label="Wald-TC and Wald-CIC · Chapter 20 §3.3 to §3.4"}
$$W_{TC} = \frac{\mathbb{E}(Y_{11}) - \mathbb{E}\big(Y_{10} + \delta_{D_{10}}\big)}{\mathbb{E}(D_{11}) - \mathbb{E}(D_{10})}, \qquad W_{CIC} = \frac{\mathbb{E}(Y_{11}) - \mathbb{E}\big[Q_{D_i}(Y_i) \mid S_i = 1, T_i = 0\big]}{\mathbb{E}(D_{11}) - \mathbb{E}(D_{10})}$$

In words: the yardstick is taken within groups defined by baseline treatment status, so A4 is no longer needed. Under A1, A2, A3', and A6, $W_{TC} = \mathrm{LATE}_{sw}$ (Theorem 2); under A1, A2, A6, A7, and A8, $W_{CIC} = \mathrm{LATE}_{sw}$ and the LQTE is identified too (Theorem 3). When A6 fails, only the TC or CIC bounds remain (Theorem 4).
:::

::: {.theorem data-label="Dose DiD: ATT(d|d) and ACRT · Chapter 20 §3.5"}
$$\mathrm{ATT}(d \mid d) = \mathbb{E}[\Delta Y \mid D = d] - \mathbb{E}[\Delta Y \mid D = 0], \qquad \frac{\partial\, \mathbb{E}[\Delta Y \mid D = d]}{\partial d} = \mathrm{ACRT}(d \mid d) + \frac{\partial\, \mathrm{ATT}(d \mid d')}{\partial d'}\Big|_{d' = d}$$

In words: A9 (each dose group's $Y(0)$ moves in parallel with the zero-dose group's) identifies the height of the dose-response curve; the slope of the curve still contains a selection bias, and identifying $\mathrm{ACRT}(d)$ requires strong parallel trends (A10) (Theorem 5). In a two-period dose design, TWFE puts negative weights on treated units with below-average doses (Chapter 20 §4.3).
:::

::: {.theorem data-label="Stayers and the heterogeneous adoption design · Chapter 20 §3.6"}
$$\mathrm{WAS} = \frac{\mathbb{E}[\Delta Y] - \lim_{d \downarrow 0} \mathbb{E}[\Delta Y \mid D_2 \leq d]}{\mathbb{E}[D_2]}$$

In words: when every unit receives a positive dose, quasi-untreated units with doses near zero serve as the yardstick, and under A12 and bounded effects, the WAS is identified (Theorem 6 (ii)); with stayers, use stayers with the same baseline dose instead (A11, Theorem 6 (i)).
:::

::: {.theorem data-label="Switching treatments and DID_ℓ · Chapter 20 §3.7"}
$$\mathrm{AVSQ}_{i,\ell} = \mathbb{E}\big[Y_{i,G_i - 1 + \ell} - Y_{i,G_i - 1 + \ell}(D_{i1}, \ldots, D_{i1})\big]$$

In words: when treatment can switch on and off and take several values, compare each unit with its own path of "staying at its period 1 level"; under A13 (no anticipation) and A14 (parallel status-quo outcomes for units with the same initial level), $DID_\ell$ is conditionally unbiased for $\mathrm{AVSQ}_\ell$ and asymptotically normal as the number of units grows (Theorem 7). In a binary absorbing design, it is numerically identical to the CS event study with not-yet-treated controls (Chapter 19 §4.5).
:::

::: {.theorem data-label="Sharp RD · Chapter 21 §3.1, §4.2 to §4.3"}
$$\tau_{\mathrm{RD}} = \lim_{x\downarrow c}\mathbb{E}[Y_i \mid X_i = x] - \lim_{x\uparrow c}\mathbb{E}[Y_i \mid X_i = x] = \mathbb{E}[\tau_i \mid X_i = c]$$

In words: with positive density on both sides of the cutoff and continuous conditional expectations of the potential outcomes (A1), the jump in the outcome is the average effect at the cutoff (Theorem 1). Estimate with local linear regression and the MSE-optimal bandwidth $h \propto n^{-1/5}$, base inference on bias-corrected robust intervals (Chapter 21 §4.2), and do not use high-order global polynomials.
:::

::: {.theorem data-label="Fuzzy RD · Chapter 21 §3.5"}
$$\tau_{\mathrm{FRD}}(c) = \frac{\lim_{x\downarrow c}\mathbb{E}[Y_i \mid X_i = x] - \lim_{x\uparrow c}\mathbb{E}[Y_i \mid X_i = x]}{\lim_{x\downarrow c}\Pr(D_i = 1 \mid X_i = x) - \lim_{x\uparrow c}\Pr(D_i = 1 \mid X_i = x)}$$

In words: when the threshold only makes the probability of treatment jump, the effect is the outcome jump divided by the treatment jump; under A1, A6, and independent gains (A6′), it is the average effect of everyone at the cutoff, and under A6 and monotonicity at the cutoff (A7) it is the LATE of compliers (Theorem 4).
:::

::: {.theorem data-label="Regression kink · Chapter 21 §3.6"}
$$\tau_{\mathrm{RK}} = \frac{\displaystyle\lim_{x\downarrow c}\frac{d}{dx}\mathbb{E}[Y_i \mid X_i = x] - \lim_{x\uparrow c}\frac{d}{dx}\mathbb{E}[Y_i \mid X_i = x]}{\displaystyle\lim_{x\downarrow c} b'(x) - \lim_{x\uparrow c} b'(x)}$$

In words: when a policy only bends its slope at the threshold, the ratio of the two kinks identifies, under A8, the average of the marginal effects of the unit types at the cutoff weighted by their relative density (Theorem 5); the boundary bias of a slope is one order larger than that of a level, so kink designs default to local quadratic fits.
:::

::: {.theorem data-label="Local randomization · Chapter 21 §3.3, §7"}
$$p_{\mathrm{F}} = \Pr\big(|T(\mathbf{Z})| \geq |T(\mathbf{z}^{\mathrm{obs}})|\big)$$

In words: treat a window around the cutoff as a small experiment in which the number of treated units in the window is fixed and every combination of which units are treated is equally likely (fixed-margins randomization); under this premise, A3, and the sharp null, the Fisher p-value is exact for any sample size in the window; adding A4 (a locally constant effect) lets you invert the test into an exact confidence set (Proposition 1). The price is that $Y(0)$ must be flat within the window, and you should report a range of windows.
:::

::: {.theorem data-label="Discrete running variables and honest intervals · Chapter 21 §3.4"}
$$\Big[\hat\tau \pm \mathrm{cv}_{1-\alpha}(r_{\sup})\,\hat\sigma_{\mathrm{NN}}/\sqrt{N_h}\Big]$$

In words: when the score takes finitely many values, the variance clustered by score value does not converge to a constant, so clustering is no remedy (Theorem 3 (i)); over the whole function class of A5 (second derivative bounded by $M$), this interval's coverage is asymptotically at least $1 - \alpha$ (Theorem 3 (ii)). $M$ cannot be estimated from the data, so run a sensitivity analysis over it.
:::

::: {.theorem data-label="Multiple cutoffs and geographic RD · Chapter 21 §3.7"}
$$\tau_{\mathrm{pool}} = \sum_{c \in \mathrm{Supp}(c_i)} \varpi(c)\,\tau_{\mathrm{RD}}(c), \qquad \varpi(c) = \frac{f_{X \mid c_i = c}(c)\Pr(c_i = c)}{\sum_{c'} f_{X \mid c_i = c'}(c')\Pr(c_i = c')}$$

In words: normalizing and pooling several cutoffs estimates an average weighted by the density at each cutoff and its share (Theorem 6, A9); a geographic RD identifies the effect at every point on the boundary, so the estimand is a curve (Theorem 7, A10), and administrative boundaries often switch several treatments at once.
:::

::: {.theorem data-label="Extrapolation away from the cutoff · Chapter 21 §3.8"}
$$\tau_{\mathrm{RD}}(\bar x, \underline{c}) = m_{\underline{c}}(\bar x) - \Big[m_{\bar{c}}(\bar x) + \lim_{x\uparrow \underline{c}} m_{\underline{c}}(x) - m_{\bar{c}}(\underline{c})\Big]$$

In words: RD estimates a local quantity at the cutoff; with a low and a high cutoff, if the difference between the two groups' $Y(0)$ does not vary with the score (A11), the effect for the low-cutoff group can be extrapolated to the region between the two cutoffs (Theorem 8). This constant bias is untestable; it can only be falsified by checking parallelism below the cutoff.
:::

::: {.theorem data-label="Synthetic control weights and the bias bound · Chapter 22 §3.1"}
$$\hat{\mathbf{w}} = \arg\min_{\mathbf{w}}\ \Big(\mathbf{X}_1 - \textstyle\sum_j w_j \mathbf{X}_j\Big)' \boldsymbol\Omega \Big(\mathbf{X}_1 - \textstyle\sum_j w_j \mathbf{X}_j\Big) \quad\text{s.t. } \sum_{j} w_j = 1,\ w_j \geq 0$$

In words: for a single treated unit, find a set of non-negative weights in the donor pool that sum to one, so that the synthetic unit's pre-period path tracks the treated unit; the post-period gap $Y_{1t} - \sum_j \hat w_j Y_{jt}$ is the effect. Under the factor model and a perfect pre-period fit (A1 to A4), the bias has an upper bound, and the bound shrinks to zero as the pre-period lengthens but grows with the number of donors (Theorem 1); the simplex constraint keeps the synthetic unit inside the convex hull, with no extrapolation.
:::

::: {.theorem data-label="Synthetic DiD · Chapter 22 §3.3"}
$$\hat\tau(\omega, \nu) = \Big(\bar Y_{1,\mathrm{post}} - \sum_{t \leq T_0}\nu_t Y_{1t}\Big) - \sum_{j}\omega_j\Big(\bar Y_{j,\mathrm{post}} - \sum_{t \leq T_0}\nu_t Y_{jt}\Big)$$

In words: a double difference that weights units and periods at once; $\omega_j = 1/J$ and $\nu_t = 1/T_0$ give DiD, and dropping the unit fixed effect $\alpha_i$ and the time weights while using SC weights without an intercept gives SC (Proposition 2). Under A5, it is asymptotically normal (Theorem 2); with a single treated unit, the intervals of `synthdid` can only come from the placebo variance (Chapter 22 §4.3).
:::

::: {.theorem data-label="Augmented SC and penalized SC · Chapter 22 §3.4"}
$$\hat Y_{1t}^{\text{aug}}(0) = \sum_j \hat w_j Y_{jt} + \Big(\hat m_t(\mathbf{X}_1) - \sum_j \hat w_j \hat m_t(\mathbf{X}_j)\Big), \qquad \zeta_{\mathrm{pen}}\sum_j w_j\big\|\mathbf{X}_1 - \mathbf{X}_j\big\|^2$$

In words: when the pre-period fit is imperfect or the treated unit lies outside the convex hull, ASC uses a ridge outcome model to subtract the remaining imbalance; its weights can be negative, and its error bound is given in Theorem 3 (A6). Penalized SC adds the right-hand term to the objective to penalize donors far from the treated unit, which eases interpolation bias; under A7 the solution is unique, with at most $\dim(\mathbf{X}_1) + 1$ nonzero weights (Theorem 4).
:::

::: {.theorem data-label="IFE and matrix completion · Chapter 22 §3.5, §7"}
$$\hat{\mathbf{L}} = \arg\min_{\mathbf{L}}\ \frac{1}{|\mathcal{O}|}\sum_{(i,t) \in \mathcal{O}}\big(Y_{it} - L_{it}\big)^2 + \zeta_{\mathrm{nn}}\|\mathbf{L}\|_*$$

In words: with many treated units, estimate the factor structure on untreated cells only and fill in $\hat Y_{it}(0)$: IFE is consistent under A8 as $N$ and $\mathcal{T}$ grow (Theorem 5); MC replaces the choice of the number of factors with a nuclear-norm penalty, and its error bound requires that each unit's probability of still being untreated in the last period not be too small (Theorem 6, A9). For IFE, report the selected $r$ and re-estimate at neighboring values of $r$.
:::

::: {.theorem data-label="Synthetic control under staggered adoption · Chapter 22 §3.6, §7"}
$$\big\|\bar{\mathbf{f}}_e\big\|_2\, q^{\mathrm{pool}} + S_e\, q^{\mathrm{sep}}$$

In words: with staggered adoption, the error bound for the average effect at event time $e$ depends both on the pre-period fit of the treated units' average path, $q^{\mathrm{pool}}$, and on each unit's own fit, $q^{\mathrm{sep}}$ (Theorem 7, A10); multisynth trades the two off, so report its tuning parameter and run a sensitivity analysis on it.
:::

::: {.theorem data-label="Conformal inference and prediction intervals · Chapter 22 §3.7"}
$$\hat p_{\mathrm{cf}} = \frac{1}{|\Pi|}\sum_{\pi \in \Pi}\mathbf{1}\{\mathrm{TS}(\hat{\mathbf{u}}_\pi) \geq \mathrm{TS}(\hat{\mathbf{u}})\}$$

In words: inference with a single treated unit is not limited to permutations across donors. The conformal test puts the randomness in the time dimension and is exact in finite samples when the data are exchangeable under the null (Theorem 8 (i)), and approximately valid under A11 (Theorem 8 (ii)); the prediction intervals of `scpi` add up the in-sample and out-of-sample uncertainty (Lemma 2). Neither is limited by the permutation resolution set by the number of donors.
:::

::: {.theorem data-label="Synthetic control for experimental design · Chapter 22 §3.8"}
$$\hat\tau_t = \sum_i w^{\mathrm{tr}}_i Y_{it} - \sum_i w^{\mathrm{co}}_i Y_{it}$$

In words: before launch, use historical data to choose weights for a treatment group and a control group at once, so that both track the population; under A12, the bound on the estimator's bias grows slowly with the number of candidate units $N$ and shrinks as the fitting period $T_E$ lengthens (Theorem 9). With a short fitting period and many candidates, the chosen combination may track history only by chance.
:::

::: {.theorem data-label="OVB sensitivity and the robustness value · Chapter 23 §2.2 to §2.4"}
$$R^2_{D\sim U\mid X},\qquad R^2_{Y\sim U\mid D,X},\qquad \widehat\Theta(B)=[\widehat\tau-B,\widehat\tau+B]$$

In words: two partial $R^2$s express the strength of an unobserved confounder's association with the treatment and with the outcome; the robustness value gives the minimum equal-strength confounding needed to overturn the conclusion, and allowing a bias bound $B$ expands the point estimate into an identified set.
:::

::: {.theorem data-label="DML: Neyman orthogonality and cross-fitting · Chapter 24 §3.2"}
$$\psi(W; \theta, \ell, m) = \big(Y - \ell(X) - \theta (D - m(X))\big)\big(D - m(X)\big), \qquad \hat\theta = \frac{\sum_i \tilde D_i\, \tilde Y_i}{\sum_i \tilde D_i^2}$$

In words: residualize both $Y$ and $D$ on $X$ ($\tilde Y, \tilde D$) and then regress one residual on the other. The orthogonal score is first-order insensitive to the nuisance, and with cross-fitting, as long as $\|\hat m - m_0\|(\|\hat m - m_0\| + \|\hat\ell - \ell_0\|) = o_p(n^{-1/2})$, the estimator is $\sqrt{n}$ asymptotically normal, so the nuisance can be estimated by any ML method.
:::

::: {.theorem data-label="Causal forest / CATE · Chapter 24 §3.4"}
$$\hat\tau(x) = \frac{\sum_i \alpha_i(x)\, (D_i - \bar D_x)(Y_i - \bar Y_x)}{\sum_i \alpha_i(x)\, (D_i - \bar D_x)^2}$$

In words: the random forest splits on treatment-effect heterogeneity, the forest weights $\alpha_i(x)$ are a data-adaptive similarity, and a local residual-on-residual regression with these weights yields $\tau(x)$; honesty (separating the splitting sample from the estimation sample) delivers pointwise asymptotic normality and confidence intervals.
:::

::: {.theorem data-label="Power / MDE / sample size · Chapter 25 §3.1"}
$$\text{power} \approx \Phi\!\left(\frac{\tau}{\text{SE}} - z_{1-\alpha/2}\right), \quad \text{MDE} = (z_{1-\alpha/2} + z_{1-\beta})\, \sigma\sqrt{\tfrac{2}{n}}, \quad n = \frac{2\sigma^2 (z_{1-\alpha/2} + z_{1-\beta})^2}{\tau^2}$$

In words: $n$ is the sample size per arm (two equal arms) and $\sigma$ the outcome standard deviation; whether an experiment can detect an effect is jointly determined by sample size, variance, and the true effect; the MDE is the smallest effect reliably detectable at a given power, and the required sample size grows in inverse proportion to the square of the target effect. An observed power computed after the fact (power backed out from the observed effect) cannot vindicate a null result (Chapter 25 §7).
:::

::: {.theorem data-label="CUPED (variance reduction) · Chapter 25 §3.2"}
$$Y_i^{cv} = Y_i - \theta (X_i - \bar X), \quad \theta = \frac{\text{Cov}(Y, X)}{\text{Var}(X)}, \quad \text{Var}(Y^{cv}) = \text{Var}(Y)\,(1 - \rho^2)$$

In words: residualizing the outcome on a pre-experiment covariate reduces the variance by the factor $1 - \rho^2$ (equivalent to ANCOVA), and the stronger the correlation, the more the effective sample size multiplies; the covariate must be measured before the experiment, or it becomes a bad control.
:::

::: {.theorem data-label="Direct and indirect effects under interference · Chapter 26 §2.2, §3.2 to §3.5, §4.3"}
$$\tau_D(g)=\mathbb E[Y_i(1,g)-Y_i(0,g)],\qquad \tau_I(d)=\mathbb E[Y_i(d,g_1)-Y_i(d,g_0)]$$

In words: $\tau_D(g)$ holds exposure fixed and changes a unit's own treatment, and individual randomization identifies it at a given saturation; $\tau_I(d)$ holds own treatment fixed and changes the neighbors' exposure, and two-stage randomization identifies it by comparing units whose own treatment is $d$ across saturations (Chapter 26 §3.3, §4.3; untreated units give $\tau_I(0)$). The total effect changes own treatment and exposure together, and cluster randomization identifies it by comparing treated and control markets (Chapter 26 §3.2); the overall effect is a different object: the difference in population means between two assignment policies. When a whole market must share one policy, use a switchback design (Chapter 26 §3.5).
:::

### Part IV Structural Methods

Part IV draws on the theory of Part I and the identification of Part III and builds structural models to study substitution, counterfactuals, and welfare. Here "identification" often rests on demand instruments, conduct assumptions, or equilibrium selection, none of which can be tested from the first-order conditions themselves.

::: {.theorem data-label="Logit choice probability and IIA · Chapter 27 §3.1 to §3.3"}
$$s_{ij} = \frac{e^{V_{ij}}}{1 + \sum_{k=1}^{J} e^{V_{ik}}}, \qquad \frac{s_{ij}}{s_{ik}} = e^{V_{ij} - V_{ik}}$$

In words: a Type-I extreme value shock gives closed-form shares, and the choice probabilities are the normalized exponentials of utility; the price is IIA: the ratio of any two probabilities does not depend on a third option, which forces proportional substitution (red-bus/blue-bus). Coefficients are identified only up to a scale normalization; what is comparable are ratios such as WTP.
:::

::: {.theorem data-label="Nested logit (GEV) · Chapter 27 §3.5"}
$$s_{ij} = \frac{e^{V_{ij}/\lambda_g}\Big(\sum_{k \in g} e^{V_{ik}/\lambda_g}\Big)^{\lambda_g - 1}}{\sum_{h}\Big(\sum_{k \in h} e^{V_{ik}/\lambda_h}\Big)^{\lambda_h}}$$

In words: it allows the shocks of products in the same nest to be correlated, with the nest parameter $\lambda_g \in (0,1]$: $\lambda_g = 1$ reduces to logit, and $\lambda_g \to 0$ makes products within a nest nearly perfect substitutes; the nesting structure must be specified in advance, and getting the nests wrong biases the substitution pattern.
:::

::: {.theorem data-label="Mixed / random-coefficients logit · Chapter 27 §3.5"}
$$s_{ij}(\theta) = \int \frac{e^{\,\beta\, x_{ij}}}{\sum_k e^{\,\beta\, x_{ik}}}\, dF(\beta \mid \theta)$$

In words: letting the coefficients $\beta_i$ vary across the population according to a distribution $F$ and integrating over this taste distribution can approximate any RUM and break IIA; the random-coefficient distribution is identified by rich variation in choice sets and characteristics, and the repeated choices in a panel can substantially strengthen identification through correlation across choice situations, though panel data are not necessary in principle.
:::

::: {.theorem data-label="Berry (1994) inversion · Chapter 28 §3.1"}
$$\ln s_{jt} - \ln s_{0t} = \delta_{jt} = x_{jt}\beta - \alpha p_{jt} + \xi_{jt}$$

In words: invert the observed market shares into the mean utility $\delta_{jt}$, which under logit is just the difference in log shares; this turns the problem into a linear IV problem with an endogenous price, and an instrument is needed wherever shares appear.
:::

::: {.theorem data-label="BLP contraction and GMM · Chapter 28 §3.4"}
$$\delta_t^{(h+1)} = \delta_t^{(h)} + \ln S_t^{\text{obs}} - \ln s_t\big(\delta_t^{(h)}, \theta_2\big), \qquad \hat\theta_2 = \arg\min_{\theta_2}\ \big(Z' \xi(\theta_2)\big)' W \big(Z' \xi(\theta_2)\big)$$

In words: under random coefficients the shares cannot be inverted analytically, so for a fixed $\theta_2$ a contraction mapping is iterated to solve for $\delta_t$ (its modulus is below 1, so it converges from any starting value), and GMM is then run on the orthogonality between the instruments and $\xi$ (Chapter 28 §4.3); the inner tolerance must be set extremely tight so that its error does not seep into the outer loop.
:::

::: {.theorem data-label="Multiproduct Bertrand FOC and markup (ℋ∘Δ)⁻¹ · Chapter 29 §3.1, §3.3, §3.5"}
$$mc = \mathbf{p} - \big(\mathcal{H} \odot \Delta(\mathbf{p})\big)^{-1} s(\mathbf{p})$$

In words: given demand, the multiproduct Bertrand first-order condition backs out marginal cost exactly, where $\mathcal{H}$ is the ownership matrix, $\Delta$ the negative of the matrix of demand price derivatives, and $\odot$ the Hadamard product. Cost is inverted rather than estimated: some set of costs can be backed out under any conduct assumption, so conduct itself cannot be identified from the first-order conditions (Chapter 29 §3.3), while a misspecified demand passes systematically through $\Delta$ into costs, merger simulations, and welfare (Chapter 29 §3.5).
:::

::: {.theorem data-label="Diversion ratio and GUPPI / UPP · Chapter 29 §3.2"}
$$D_{j \to k} = \frac{\partial s_k / \partial p_j}{-\,\partial s_j / \partial p_j}, \qquad GUPPI_j = \frac{1}{p_j} \sum_{k \in \mathcal{J}_{f'}} (p_k - mc_k)\, D_{j \to k}$$

In words: the diversion ratio measures how much of the volume lost by $j$ after a price increase goes to $k$ (equal to the second-choice share) and is the core demand output for merger analysis; GUPPI combines it with the markups of the merger partner's products ($\mathcal{J}_{f'}$) into a dimensionless first-pass merger screen, and UPP then subtracts an efficiency offset (Chapter 29 §4.3).
:::

::: {.theorem data-label="Merger simulation and log-sum welfare · Chapter 29 §4.2, §4.4"}
$$\Delta CS = \int \frac{1}{\alpha_i}\Big[\log\big(1 + \textstyle\sum_j e^{V_{ij}^{\text{post}}}\big) - \log\big(1 + \textstyle\sum_j e^{V_{ij}^{\text{pre}}}\big)\Big]\, dF_i$$

In words: a merger flips the corresponding cross-terms of the ownership matrix from 0 to 1, and re-solving the FOC with $mc$ and demand held fixed gives the post-merger prices; $\Delta CS$ is the change in the log-sum from before to after the merger, divided by $\alpha_i$ (the marginal utility of income, that is, the price coefficient, which must not vary with income), and is signed as a gain, the opposite of Chapter 1's $\Delta\mathrm{CS}$.
:::

::: {.theorem data-label="Entry games and Bresnahan-Reiss · Chapter 30 §2, §4.1"}
$$\Pi_{fm} = X_m \beta + Z_{fm} \alpha - \delta\, N_{-f,m} + \xi_m + \varepsilon_{fm}, \qquad a_{fm} = \mathbf{1}\{\Pi_{fm} \geq 0\}$$

In words: a firm enters when its profit from entering is non-negative, with $\delta > 0$ the intensity of competition; Bresnahan-Reiss sidesteps the game's multiple equilibria by running an ordered probit on the number of entrants alone, and uses the ratio of adjacent per-firm entry thresholds $s_{N+1}/s_N$ to gauge how competition intensifies with the number of incumbents.
:::

::: {.theorem data-label="Partial identification (Ciliberto-Tamer) · Chapter 30 §3.4, §4.4"}
$$H_1(\theta, X) \leq P(y \mid X) \leq H_2(\theta, X), \qquad Q(\theta) = \sum_y \Big[ \big(P(y) - H_1\big)_-^2 + \big(P(y) - H_2\big)_+^2 \Big]$$

In words: when the entry game has multiple equilibria, the model gives only an interval for the outcome probability ($H_1$ the lower bound from the unique-equilibrium probability, $H_2$ the upper bound from some equilibrium), the parameters at which the criterion function is zero form the identified set, and no assumption needs to be imposed on the unknown equilibrium-selection mechanism.
:::

::: {.theorem data-label="Bellman and the Rust contraction · Chapter 31 §2, §3.1 to §3.2"}
$$V(x) = \mathbb{E}_\varepsilon \max_{i}\Big\{ u(x, i) + \varepsilon_i + \beta\, \mathbb{E}\big[V(x') \mid x, i\big] \Big\}$$

$$(T_\theta V)(x) = \log\Big[ \exp\big(u(x,0) + \beta\, \mathbb{E}[V(x') \mid x, 0]\big) + \exp\big(u(x,1) + \beta\, \mathbb{E}[V(x') \mid x, 1]\big) \Big]$$

In words: the value function is its own fixed point, and under conditional independence plus Type-I EV, the expected maximum has a closed log-sum form (up to Euler's constant), so the Bellman equation becomes the operator $T_\theta$, a contraction with modulus $\beta < 1$ that iterates to a unique solution; NFXP maximizes the likelihood in the outer loop and solves the fixed point in the inner loop.
:::

::: {.theorem data-label="Hotz-Miller CCP inversion · Chapter 31 §4.3"}
$$v(x, 1) - v(x, 0) = \log P(1 \mid x) - \log P(0 \mid x)$$

In words: under Type-I EV the value difference between two actions is in one-to-one correspondence with the ratio of choice probabilities, so the dynamic value difference can be recovered from observed conditional choice probabilities without solving the fixed point; this is the key to two-step estimation; the discount factor $\beta$ is generally not identified and must be fixed or pinned down by a separate exclusion restriction (Chapter 31 §3.3).
:::

### Part V Topics in the Digital Economy

Part V brings the three earlier legs to bear together on the phenomena of the digital economy: two-sided platforms, consumer search, algorithms, and data. The methods are mostly combinations and applications of the earlier chapters.

::: {.theorem data-label="Two-sided pricing (Armstrong / Rochet-Tirole) · Chapter 32 §2.1"}
$$p_i = f_i - \alpha_j\, n_j + \frac{n_i}{\partial n_i / \partial u_i}$$

In words: the platform's optimal price for one group of members equals the cost of connecting it, minus the cross-side external benefit $\alpha_j n_j$ that this group brings to the other, plus the usual monopoly markup; the middle term can push one side's price below cost or even turn it into a subsidy, and what matters is the price structure (how the total price is split across the two sides), not the price level.
:::

::: {.theorem data-label="Network effects (Katz-Shapiro) and tipping · Chapter 32 §2.2"}
$$P(q) = r(q) + v(q), \qquad P'(q^\ast) = r'(q^\ast) + v'(q^\ast) < 0$$

In words: the marginal user's willingness to pay equals a decreasing stand-alone value plus an increasing network value; fulfilled-expectations demand can therefore be non-monotone, and where it crosses the price line at several points there are multiple equilibria. A stable equilibrium requires the demand curve to cross the price line from above; otherwise tipping and a critical mass appear.
:::

::: {.theorem data-label="Switching costs / state dependence · Chapter 32 §2.6"}
$$u_{ijt} = x_{jt}'\beta + \alpha_{ij} + \delta\, \mathbf{1}\{y_{i,t-1} = j\} + \varepsilon_{ijt}$$

In words: $\delta$ is genuine structural state dependence (a switching cost, monetized as $\delta / |\beta|$, with $\beta$ the utility response to the fee), which must be separated from the spurious state dependence created by persistent preferences $\alpha_{ij}$; failing to control for heterogeneity systematically overstates $\delta$, and the initial-conditions problem must also be handled.
:::

::: {.theorem data-label="Sequential search (Weitzman reservation value) · Chapter 33 §2.2"}
$$c_j = \mathbb{E}\big[\max(u_j - z_j, 0)\big]$$

In words: each option is assigned a reservation value $z_j$ by its own search-cost equation, and the optimal rule reduces to a set of static indices: inspect in decreasing order of $z_j$, and stop once the realized utility exceeds every unopened index. The higher the search cost, the lower the $z$ and the shallower the search.
:::

::: {.theorem data-label="Limited consideration sets and identifying search costs · Chapter 33 §3.2 to §3.4, §4.2"}
$$\Pr(\text{buy } j \mid \text{choice set} = \mathcal{S}_i) = \frac{\exp(-\alpha p_j)}{\sum_{k \in \mathcal{S}_i} \exp(-\alpha p_k)}, \qquad j \in \mathcal{S}_i$$

In words: a full-information logit mistakes a "cheap item never inspected" for one "seen and rejected", which attenuates the price coefficient $\alpha$ toward zero. Replacing the choice set with the set actually inspected, $\mathcal{S}_i$ (the formula above), points in the right direction but overstates, because stopping is endogenous. The correct fix is to model the search process explicitly: randomized ordering (position affects purchase only through "whether it was inspected") makes the order of inspection exogenous, so that once the stopping rule is modeled, search costs and preferences are identified jointly.
:::

::: {.theorem data-label="Algorithmic pricing and Q-learning (collusion) · Chapter 34 §2.1"}
$$Q(s_t, a_t) \leftarrow (1 - \alpha)\, Q(s_t, a_t) + \alpha\Big[ r_t + \delta\, \max_{a'} Q(s_{t+1}, a') \Big], \qquad \Delta = \frac{\bar{p} - p^N}{p^M - p^N}$$

In words: with all firms' previous-period prices as the state and the current profit as the reward, Q-learning updates the action value by bootstrapping; the collusion index $\Delta$ normalizes the average price between the competitive price $p^N$ and the monopoly price $p^M$, and in simulations the algorithms can learn supra-competitive pricing on their own, but "colluding in simulation" is not the same as "colluding in reality".
:::

::: {.theorem data-label="AI adoption and organizational complementarity · Chapter 34 §2.3, §3.5, §4.6"}
$$Y=F(K,L,A,O),\qquad F_{AO}>0,\qquad \tau_{AO}=(\bar Y_{11}-\bar Y_{01})-(\bar Y_{10}-\bar Y_{00})$$

In words: the cross-marginal product of AI capital $A$ with organizational capital $O$ (workflow, skills, data readiness, and the like) is positive; the cleanest empirical design is a factorial experiment that randomizes AI access and the organizational intervention separately, and the factorial contrast $\tau_{AO}$ identifies the policy-relevant complementarity.
:::

::: {.theorem data-label="Generated regressors: attenuation and correction / PPI · Chapter 34 §3.1, §4.3 to §4.4"}
$$\mathrm{plim}\; \hat\beta_{\text{OLS}} = \beta \cdot \frac{\sigma_{x^\ast}^2}{\sigma_{x^\ast}^2 + \sigma_u^2} = \beta \cdot \lambda, \qquad \hat\beta_{\text{corrected}} = \frac{\hat\beta_{\text{OLS}}}{\hat\lambda}$$

In words: when a machine-predicted variable is put on the right-hand side as if it were the truth, classical measurement error attenuates the coefficient toward zero by the factor $\lambda$ (the reliability), and the attenuation does not vanish as the sample grows; a small, randomly drawn set of gold labels makes $\lambda$ or the systematic bias estimable and hence correctable (PPI uses a rectifier to subtract the bias), and the standard errors must incorporate the uncertainty of the generation step.
:::

## 2 A Method Decision Tree

The first step in choosing a method is not to ask "which model is the most advanced" but "what do I want to answer, what do the data look like, and where does the variation come from". The tree below starts from the research goal, follows the data structure to the family of methods you should use, and labels each exit with the corresponding chapter and section. It is only an introductory guide: the real choice always turns on a substantive judgment about the treatment-assignment mechanism and the data-generating process, and with the same data, a different credible assumption can lead to a different family of methods.

First decide whether the goal is description and prediction, a causal effect, or counterfactuals and welfare; whichever it is, estimation and inference first require the right machinery from Part II:

- **First settle the estimation machinery** (Part II, used by every design)
    - The mean is a nonlinear function of the parameters → NLS (Chapter 5 §3.5); non-negative outcomes with many zeros → Poisson PML rather than log(1+y) (Chapter 5 §3.7 to §3.8)
    - Censored, truncated, or duration data → MLE that writes the observation mechanism into the likelihood, with QMLE and the sandwich when the specification is in doubt (Chapter 6 §3.2, §3.6)
    - Several moment conditions restricting a few parameters → GMM and the J test (Chapter 7 §3.3 to §3.4); reduced-form estimates already in hand → minimum distance (Chapter 7 §3.7)
    - Simultaneous equations, or several equations with correlated errors → 2SLS, SUR, 3SLS, and system GMM (Chapter 9 §3.3, §3.6 to §3.7)
    - Panel data with time-invariant individual heterogeneity → FE, RE, FD, or Mundlak (CRE) with a robust Hausman test (Chapter 10 §3.2 to §3.6); a lagged dependent variable on the right-hand side → Arellano-Bond or Blundell-Bond (Chapter 11 §3.3, §3.5)
    - A binary, ordered, or multinomial outcome → probit and logit with APEs, a control function for an endogenous regressor, FE logit or CRE probit for panels (Chapter 12 §3.3, §3.5, §3.7 to §3.9)
    - You care about the whole distribution rather than the mean → conditional QR, QTE, RIF regression (Chapter 13 §3.3, §3.5, §3.6); IV-QR for an endogenous treatment (Chapter 13 §3.7), and, with a control group, CIC to compare the treated group's distributions before and after a policy (Chapter 13 §3.9)
    - Standard errors: clustering, HAC, the bootstrap, and fixes for few clusters (Chapter 8 §3.3 to §3.5, §4.2)
- **You only want to describe, predict, or characterize a relationship curve** (no causal claim)
    - You want the data to decide the functional form and to avoid specification error → kernel estimation / local linear / splines / partially linear (Chapter 14 §3.1 to §3.4, §4.3)
    - You want to make or evaluate decisions with machine predictions or labels → first check for selective labels (Chapter 17 §2.5), and do not use the same gold labels both to train the model and to estimate the correction (Chapter 17 §3.4); guard against attenuation when a machine-generated variable enters a regression (Chapter 34 §3.1)
- **You want the causal effect of some treatment** (treatment effect). First look at how treatment is assigned:
    - First use Chapter 16 §2.1 to fix the unit, treatment, version, timing, exposure, outcome, population, and assignment; when a variable comes from platform logs or machine labels, first run the measurement audit of Chapter 17 §3.2 to §3.3 (checking the machine variable's errors against a random gold sample)
    - Treatment is randomly assigned (RCT / online A/B test) → difference in means and randomization inference (Chapter 15 §4.1); power, MDE, CUPED, sequential testing (Chapter 25 §3.1 to §3.3)
        - Random assignment but incomplete take-up → the difference in means estimates the ITT; for the effect of take-up itself, use assignment as the instrument to estimate the LATE of compliers (Chapter 18 §3.4)
    - Observational data with non-random treatment; next, look for the source of identification:
        - You believe the confounders are all observed (selection on observables) → matching / IPW / AIPW (Chapter 15 §4.3 to §4.5); high-dimensional covariates or complex nuisance → DML (Chapter 24 §3.2)
        - Treatment is endogenous but there is an instrument that is relevant, exogenous, and satisfies exclusion and monotonicity → IV / LATE (Chapter 18 §3.4); you want the marginal population and extrapolation weights → MTE and control functions (Chapter 18 §3.6 to §3.7)
        - A running variable that triggers treatment only once it crosses a cutoff → sharp RD (Chapter 21 §3.1); a threshold that only changes the probability of treatment → fuzzy RD (Chapter 21 §3.5); a threshold that only changes the slope of a policy → regression kink (Chapter 21 §3.6)
            - The running variable takes only a few values → honest intervals (Chapter 21 §3.4) or local randomization between the nearest support points on either side of the cutoff (Chapter 21 §3.3)
            - Several cutoffs or a geographic boundary → the weights of normalized pooling and geographic RD (Chapter 21 §3.7); you want the effect away from the cutoff → multi-cutoff extrapolation (Chapter 21 §3.8)
        - Panel data with pre- and post-periods and a control group → 2×2 DiD and event studies (Chapter 19 §3.1, §4.2)
            - Staggered timing → CS, Sun-Abraham, BJS, or $DID_\ell$ (Chapter 19 §4.5): static TWFE is biased when effects vary over time (Chapter 19 §4.3), and when effects differ across cohorts, the event-study coefficients of dynamic TWFE contaminate one another, so even pre-period coefficients can be nonzero (Chapter 19 §4.4); parallel trends holds only approximately → HonestDiD (Chapter 19 §6.9)
            - Only some units in the treated group take up treatment (fuzzy, partial take-up) → Wald-DID; when the effect for the incumbent treated may change, Wald-TC or Wald-CIC (Chapter 20 §3.1, §3.3 to §3.4)
            - Dose or continuous treatment → with zero-dose units, estimate $\mathrm{ATT}(d \mid d)$, and the slope also requires strong parallel trends (Chapter 20 §3.5); with no zero-dose units → HAD (Chapter 20 §3.6)
            - Treatment can switch on and off or take several values (non-absorbing) → $DID_\ell$ (Chapter 20 §3.7)
        - Only one (or very few) treated units, plus a set of untreated donors and a long enough pre-period → synthetic control (Chapter 22 §3.1)
            - Poor pre-period fit, or the treated unit outside the convex hull → augmented SC; donors very different from the treated unit → penalized SC (Chapter 22 §3.4)
            - Inference → in-space and in-time placebos (Chapter 22 §4.2), or the conformal test and `scpi` prediction intervals (Chapter 22 §3.7)
            - Not yet launched and choosing pilot units → synthetic design (Chapter 22 §3.8)
        - Several treated units or staggered adoption, with the counterfactual built from donor weights or a factor structure → multisynth (Chapter 22 §3.6) or synthetic DiD (Chapter 22 §3.3); trends driven by a factor structure → IFE or matrix completion (Chapter 22 §3.5)
    - You care about heterogeneous effects (who benefits most) or the optimal assignment policy (whom to treat) → causal forest / GATES / policy tree (Chapter 24 §3.4 to §3.6)
    - The core identifying assumption holds only approximately and you want to quantify how large a departure would overturn the conclusion → sensitivity, negative controls, bounds, and the breakdown frontier (Chapter 23 §2.2 to §2.5, §3.1 to §3.2)
    - One unit's treatment changes the outcomes of other units, or the goal is a market-wide rollout → cluster / two-stage / switchback marketplace experiments (Chapter 26 §3.2 to §3.5)
- **You want counterfactuals, substitution patterns, market power, or welfare** (structural modeling). First look at what data you have:
    - Individual-level choice data (who chose what among several options) → discrete choice, moving from logit to nested logit to mixed logit as the substitution pattern grows richer (Chapter 27 §3.1, §3.5)
    - Only market-level shares and prices, with endogenous prices → BLP demand: the Berry inversion plus the contraction plus demand instruments (Chapter 28 §3.1, §3.4)
    - Demand already estimated, and you want to recover markups and marginal costs, simulate a merger, or test conduct → supply-side first-order conditions and merger simulation (Chapter 29 §3.1, §3.3 to §3.4, §4.2, §4.5)
    - You want to explain how many firms a market has and their entry and exit decisions → if only the total intensity of competition matters, model the number of entrants with the Bresnahan-Reiss threshold model or Berry's model of heterogeneous firms (Chapter 30 §4.1, §4.2); for firm-specific competitive effects, use Seim if you accept incomplete information (Chapter 30 §4.3), or Ciliberto-Tamer partial identification if you make no assumption on equilibrium selection (Chapter 30 §3.4, §4.4)
    - The decision is dynamic: today's choice changes tomorrow's state (replacing equipment, adopting technology, optimal stopping) → single-agent Rust NFXP or Hotz-Miller CCP (Chapter 31 §4.1, §4.3)
- **Specific phenomena of the digital economy** (combining the tools above)
    - Two-sided platform pricing, network effects and tipping, switching costs and lock-in → Chapter 32 §2.1, §2.2, §2.6
    - Consumer search frictions, limited consideration sets, advertising and experience goods → Chapter 33 §2.2, §2.4 to §2.5, §3.1, §3.3, §3.6
    - Algorithmic pricing, AI adoption, organizational complements, inference with machine-generated variables, and the economics of data → Chapter 34 §2.1 to §2.3, §3.1 to §3.5

The language of counterfactuals and welfare comes from Part I: consumer surplus, CV/EV, and market power (Chapter 1 §4.1, Chapter 4 §2).

For a quick lookup by data shape, use the table below:

| Data / situation | First-choice method | Chapter |
|---|---|---|
| Non-negative outcomes with many zeros; censoring or durations | PPML; MLE, QMLE | Chapters 5, 6 |
| Several moment conditions; simultaneous equations | GMM, MD; 2SLS, 3SLS | Chapters 7, 9 |
| Panel + individual effects; lagged dependent variable | FE, RE, FD, CRE; AB, BB | Chapters 10, 11 |
| Binary or multinomial outcomes; distributional effects | APE, control function; QR, QTE, RIF | Chapters 12, 13 |
| Randomized experiment, online A/B test | difference in means + randomization inference; power, CUPED, sequential testing | Chapters 15, 25 |
| Platform logs / machine labels | measurement audit, random validation sample, PPI | Chapter 17 |
| Observational data, confounders observed | matching / IPW / AIPW; DML in high dimensions | Chapters 15, 24 |
| Endogenous treatment + valid instrument | IV / LATE / MTE / control function | Chapter 18 |
| Running variable + cutoff | sharp, fuzzy, kink; honest or local randomization for a discrete score; multiple cutoffs and extrapolation | Chapter 21 |
| Panel + staggered treatment timing | DiD / event study; CS / SA / BJS / $DID_\ell$ when staggered; HonestDiD | Chapter 19 |
| Partial take-up, dose, or switching treatment | Wald-DID / TC / CIC; $\mathrm{ATT}(d \mid d)$, HAD; $DID_\ell$ | Chapter 20 |
| Few treated units + donor pool | SC, ASC, penalized SC, SDID; conformal, `scpi` | Chapter 22 |
| Many treated units + factor structure | multisynth, IFE, matrix completion | Chapter 22 |
| Assumptions hold only approximately / hidden bias | sensitivity, negative controls, partial ID | Chapter 23 |
| Want heterogeneous effects / policy | causal forest / policy learning | Chapter 24 |
| Spillover / marketplace rollout | cluster, two-stage, switchback | Chapter 26 |
| Individual choice data, want substitution elasticities | logit / nested / mixed logit | Chapter 27 |
| Market shares + endogenous prices | BLP demand | Chapter 28 |
| Demand in hand, want markups / mergers | supply-side FOC / merger simulation | Chapter 29 |
| Entry / market structure | Bresnahan-Reiss / partial identification | Chapter 30 |
| Dynamic single-agent decisions | Rust NFXP / Hotz-Miller CCP | Chapter 31 |

## 3 Cross-Software Command Map

Each row of the table below is one method, giving the main commands in R, Python, and Stata, with R as the primary language. Every R function was confirmed to exist, one by one, in R 4.6.1, and every Python and Stata command was checked against its official documentation or package page. A blank cell means we could not verify a widely used standard implementation; the corresponding chapter mostly implements the formula by hand (base R's `optim`, `glm`, and matrix algebra suffice), and we leave the cell blank rather than invent an API. Package and function names may change across versions, so check the respective documentation before use.

| Method (chapter) | R | Python | Stata |
|---|---|---|---|
| OLS and high-dimensional fixed effects (10) | `fixest::feols`, `plm::plm` | `linearmodels` PanelOLS, `pyfixest` | `reghdfe`, `xtreg` |
| NLS and the sandwich (5) | `nls`, `sandwich::sandwich` | `scipy.optimize.least_squares` | `nl` |
| Poisson PML (5) | `fixest::fepois`, `glm` | `pyfixest`, `statsmodels` Poisson | `ppmlhdfe`, `poisson` |
| Two-step M-estimation, bootstrap both steps together (5) | `boot::boot` | `scipy.stats.bootstrap` | `bootstrap` |
| MLE, Tobit, and durations (6) | `maxLik::maxLik`, `censReg::censReg`, `survival::survreg` | `statsmodels` GenericLikelihoodModel, `lifelines` | `ml`, `tobit`, `streg` |
| LR and Wald tests (6, 8) | `lmtest::lrtest`, `lmtest::waldtest` |  | `lrtest`, `test` |
| GMM and the J test (7) | `gmm::gmm`, `momentfit::gmmFit` | `linearmodels` IVGMM | `gmm` |
| Classical minimum distance (7) |  |  |  |
| Robust and cluster standard errors (8) | `sandwich::vcovHC`, `sandwich::vcovCL`, `clubSandwich::vcovCR` | `statsmodels` `cov_type="cluster"` | `vce(robust)`, `vce(cluster)` |
| HAC (8) | `sandwich::NeweyWest` | `statsmodels` `cov_type="HAC"` | `newey` |
| Bootstrap and wild cluster bootstrap (8) | `boot::boot`, `fwildclusterboot::boottest` | `scipy.stats.bootstrap`, `wildboottest` | `bootstrap`, `boottest` |
| Delta method (8) | `car::deltaMethod` |  | `nlcom` |
| 2SLS and weak instruments (9) | `ivreg::ivreg`, `fixest::feols`, `ivmodel::AR.test` | `linearmodels` IV2SLS | `ivregress 2sls`, `ivreg2`, `weakiv` |
| SUR, 3SLS (9) | `systemfit::systemfit` | `linearmodels` SUR, IV3SLS | `sureg`, `reg3` |
| FE, RE, FD, CRE, and Hausman (10) | `plm::plm`, `plm::phtest` | `linearmodels` PanelOLS, RandomEffects, FirstDifferenceOLS | `xtreg`, `hausman`, `xthtaylor` |
| Dynamic panel GMM (11) | `plm::pgmm`, `plm::mtest` | `pydynpd` | `xtabond2`, `xtabond`, `xtdpdsys` |
| Binary response and APEs (12) | `glm`, `marginaleffects::avg_slopes` | `statsmodels` Probit, get_margeff | `probit`, `logit`, `margins` |
| Control function and IV probit (12) |  |  | `ivprobit` |
| FE logit, bias correction, and CRE (12) | `survival::clogit`, `alpaca::feglm`, `alpaca::biasCorr`, `bife::bife` | `statsmodels` ConditionalLogit | `clogit`, `xtlogit, fe`, `xtprobit, re` |
| Ordered and multinomial response (12) | `ordinal::clm`, `MASS::polr`, `nnet::multinom`, `mlogit::hmftest` | `statsmodels` OrderedModel, MNLogit | `oprobit`, `mlogit`, `suest` |
| Quantile regression (13) | `quantreg::rq` | `statsmodels` QuantReg | `qreg`, `sqreg` |
| RIF, QTE, and CIC (13) | `rifreg::rifreg`, `qte::unc_qte`, `qte::cic` |  | `rifhdreg`, `cic` |
| Nonparametric and kernel regression (14) | `np::npreg`, `KernSmooth::locpoly`, `mgcv::gam` | `statsmodels` KernelReg | `npregress`, `lpoly` |
| Matching (15) | `MatchIt::matchit`, `Matching::Match` |  | `teffects psmatch`, `teffects nnmatch` |
| IPW / AIPW (15) | `WeightIt::weightit`, `AIPW::AIPW` | `doubleml` | `teffects ipw`, `teffects aipw` |
| Validation / PPI (17, 34) |  | `ppi_py` |  |
| IV / LATE and MTE (18) | `ivreg::ivreg`, `fixest::feols` | `linearmodels` IV2SLS | `ivregress 2sls`, `mtefe` |
| 2×2 DiD and event studies (19) | `fixest::feols`, `fixest::sunab` | `pyfixest` | `reghdfe`, `eventstudyinteract` |
| Staggered DiD: CS, BJS, two-stage, ETWFE, $DID_\ell$ (19) | `did::att_gt`, `did::aggte`, `didimputation::did_imputation`, `did2s::did2s`, `etwfe::etwfe`, `DIDmultiplegtDYN::did_multiplegt_dyn` | `csdid`, `differences`, `pyfixest` | `csdid`, `did_imputation`, `did2s`, `jwdid`, `did_multiplegt_dyn` |
| Goodman-Bacon decomposition (19) | `bacondecomp::bacon` |  | `bacondecomp` |
| HonestDiD (19) | `HonestDiD::createSensitivityResults`, `HonestDiD::createSensitivityResults_relativeMagnitudes` |  | `honestdid` |
| Fuzzy DiD (20) | `Rfuzzydid::fuzzydid` |  | `fuzzydid` |
| Doses, stayers, and HAD (20) | `contdid::cont_did`, `DIDmultiplegtSTAT::did_multiplegt_stat`, `DIDHAD::did_had` |  | `did_multiplegt_stat`, `did_had` |
| TWFE weight diagnostics (20) | `TwoWayFEWeights::twowayfeweights` |  | `twowayfeweights` |
| Sharp and fuzzy RD, kink (21) | `rdrobust::rdrobust`, `rdrobust::rdplot`, `rddensity::rddensity` | `rdrobust`, `rddensity` | `rdrobust`, `rddensity` |
| Local randomization, honest intervals, multiple cutoffs (21) | `rdlocrand::rdrandinf`, `rdlocrand::rdwinselect`, `RDHonest::RDHonest`, `rdmulti::rdmc`, `rdmulti::rdms` | `rdlocrand`, `rdmulti` | `rdrandinf`, `rdmc` |
| Synthetic control and SDID (22) | `tidysynth::synthetic_control`, `Synth::synth`, `synthdid::synthdid_estimate` | `pysyncon` | `synth`, `sdid` |
| ASC, multisynth, IFE, and MC (22) | `augsynth::augsynth`, `augsynth::multisynth`, `fect::fect`, `gsynth::gsynth` |  |  |
| SC prediction intervals and inference (22) | `scpi::scest`, `scpi::scpi` | `scpi_pkg` | `scpi` |
| Synthetic control for experimental design (22) |  |  |  |
| OVB sensitivity (23) | `sensemakr::sensemakr` | `PySensemakr` | `sensemakr` |
| DML (24) | `DoubleML::DoubleMLPLR`, `DoubleML::DoubleMLIRM` | `doubleml`, `econml` | `ddml` |
| Causal forest and policy learning (24) | `grf::causal_forest`, `policytree::policy_tree` | `econml` CausalForestDML, DRPolicyTree |  |
| Power / MDE (25) | `power.t.test` | `statsmodels` TTestIndPower | `power twomeans` |
| Sequential testing boundaries (25) | `gsDesign::gsDesign` |  | `gsdesign` |
| Marketplace / switchback (26) | `DeclareDesign::declare_model`, `DeclareDesign::diagnose_design` |  |  |
| Logit / conditional logit (27) | `mlogit::mlogit`, `survival::clogit` | `xlogit`, `biogeme` | `clogit`, `cmclogit` |
| Nested logit (27) | `mlogit::mlogit` | `biogeme` | `nlogit` |
| Mixed / RC logit (27) | `mlogit::mlogit` | `xlogit` | `cmmixlogit`, `mixlogit` |
| BLP demand (28) | `BLPestimatoR::estimateBLP` | `pyblp` |  |
| Markups and merger simulation (29) | `antitrust::logit`, `antitrust::upp.bertrand` | `pyblp` |  |
| Entry and structural partial identification (30) |  |  |  |
| Dynamic structural Rust / CCP (31) |  | `ruspy` |  |
| RE and CRE logit (32) | `lme4::glmer` |  | `xtlogit, re` |
| Consumer search (33) | `survival::clogit` |  |  |
| Measurement error and generated regressors (34) | `simex::simex` | `ppi_py` |  |

## 4 Common Pitfalls Checklist

Below we gather the most common errors from each chapter's "failure modes" section, grouped by method, for an item-by-item self-check before and after you run the analysis. Each item names only the problem and the direction of the fix; the full argument is in the corresponding chapter.

### Part II Econometric Core

- A misspecified mean (Chapter 5): the point estimate does not "blow up"; it targets the best approximation picked by the loss. Make the mean more flexible, run a RESET-type test, and put NLS and PPML side by side (a large gap between them is an alarm); base inference on the fully robust variance with the sample Hessian as bread. Report the units behind any log(1+y) coefficient, and do not read it as a percentage.
- Wrong standard errors (Chapters 5 and 6): the model's default standard errors are right only when $B_0 = \sigma^2 A_0$ or the whole likelihood is correctly specified, so report the sandwich by default; OPG often performs badly in finite samples; LR relies on the information matrix equality, so when the specification is in doubt use the robust Wald or LM test instead. A suspiciously narrow standard error often signals a wrong formula.
- Two-step estimation and the observation mechanism (Chapters 5 and 6): treating the first step as known generally makes the naive standard errors wrong, in a direction that depends on the case, so correct them with Chapter 5 Theorem 6 or bootstrap the two steps together; censoring must be independent of duration given the covariates (CENS.1), and sampling only the people still subscribed on a given date creates length-biased sampling.
- Moment conditions (Chapter 7): randomization does not guarantee the exclusion restrictions; J only detects contradictions among moments, so when it rejects, first ask whether the cause is misspecification or effect heterogeneity (under the latter, the estimate targets a precision-weighted average); with many moments and a small sample, report Windmeijer standard errors under linear moment conditions with CUE as a cross-check, use the S confidence set for weak moments, and look first to the design for a remedy.
- The level of the standard errors (Chapter 8): heteroskedasticity-robust standard errors do nothing about within-cluster correlation, so cluster at the level of treatment assignment; the asymptotics are in the number of clusters $G$, and with few clusters (as a rule of thumb, fewer than thirty or forty) tests over-reject, so upgrade to $t(G-1)$, CR2, or the wild cluster bootstrap.
- Delta method and bootstrap (Chapter 8): when the denominator is near zero, the symmetric interval breaks down, so switch to a Fieller interval (percentile-t is unreliable there too); on correlated data, a naive bootstrap amounts to assuming iid, so use the block or wild cluster bootstrap; correct for multiple testing.
- IV and system estimation (Chapter 9): with exact identification, the exclusion restriction is completely untestable; with weak instruments, 2SLS is pulled toward OLS, so report the effective F and the AR set, and do not pick specifications by whether F clears 10; use the robust J for overidentification; a misspecified equation contaminates 3SLS, so report each equation's $J$ alongside equation-by-equation 2SLS; with heterogeneous effects, 2SLS estimates a weighted average under the two conditions of Chapter 9 §3.2; with several instruments, the condition that they push in the same direction must hold over their joint values.
- Panels (Chapter 10): default to FE, and use the robust Mundlak regression instead of the classical Hausman test; when covariates respond to past shocks, strict exogeneity fails, so run a lead test and compare FE with FD; when slopes are correlated with the regressor's within-unit variation, FE weights units by their within variance and quietly changes the estimand; attrition correlated with idiosyncratic shocks, and measurement error amplified by demeaning, both hurt FE.
- Dynamic panels (Chapter 11): in a short panel, a lagged dependent variable plus fixed effects produces the Nickell bias, and the regression runs without any warning; the power of $m_2$ is limited; on platforms where new merchants or new cities are still ramping up, BB.1 fails and system GMM can be worse, so report it alongside difference GMM; the number of instruments grows with the square of $T$, and with more instruments than units, the $J$ of a one-step fit equals the number of units and the two-step $J$ cannot exceed it, so a $p$-value near 1 is an alarm: collapse, limit the lag depth, and report the instrument count.
- Binary and multinomial response (Chapter 12): complete separation makes the MLE diverge, and the software may not warn you; the scale of the coefficients is set by the latent variance, so comparing coefficients across groups compares scales, and the sign of an interaction coefficient is not the sign of the cross-partial effect, so report APEs; a control function requires a continuous $D_i$, and plugging fitted probabilities of a binary endogenous regressor into a probit is inconsistent; FE logit gives no APE; the classical Hausman-McFadden test is distorted in finite samples, so use the version with stacked scores or estimate a nested logit directly.
- Quantiles (Chapter 13): a QTE is a shift of the distribution, not any person's gain or loss, and the latter requires rank invariance; a conditional QR coefficient is not a shift of the platform-wide distribution, for which you need RIF or distribution regression; intervals at tail quantiles are wide, and pointwise intervals do not license "the whole curve is significant"; QR curves fitted at different levels can cross, so rearrange them back into monotone order.
- Nonparametric smoothing and inference (Chapter 14): too large a bandwidth smooths away real bends, and too small a bandwidth chases noise; the boundary and sparse regions are the least trustworthy (kernel regression has boundary bias, so use local linear); at the optimal bandwidth, the bias and the standard deviation are of the same order, so you must undersmooth or correct the bias explicitly for the confidence interval to cover; beyond three or four continuous regressors the curse of dimensionality bites, so move to semiparametric methods.

### Part III Causal Identification Toolbox

- Selection on observables / matching (Chapter 15): the CIA is untestable, and balance and placebo checks are only indirect evidence; when overlap is poor, the IPW weights explode, and trimming changes the estimand; the bootstrap is invalid for NN matching with a fixed number of matches, so use the Abadie-Imbens variance; do not control for downstream variables that respond to treatment (bad controls).
- Research design (Chapter 16): first fix the unit, treatment, version, timing, exposure, outcome, population, and assignment; post-treatment controls can both change the estimand and open a collider path; the arrows in a DAG come from institutional assumptions, not from facts that a correlation matrix verifies automatically.
- Measurement (Chapter 17): high predictive accuracy does not guarantee downstream inference; the gold sample must represent the target population and carry reliable labels; the bias of a machine variable differs depending on whether it is on the left-hand or the right-hand side; PPI fixes measurement, not treatment endogeneity.
- IV / LATE (Chapter 18): of the four assumptions only relevance is testable, and exclusion is the one most in need of argument and the hardest to see; weak instruments are the number-one threat, so report the first-stage F (the effective F under heteroskedasticity or clustering) and base the main inference on Anderson-Rubin; using the fitted values of a nonlinear first stage such as a probit directly as a regressor is consistent only if that first stage is correctly specified, and plugging in a nonlinear transform of them is the forbidden regression, which is in general inconsistent; the LATE is only the effect for compliers, so do not extrapolate it as the ATE.
- DiD (Chapter 19): under staggered timing, the weights of static TWFE are contaminated by forbidden comparisons, which bias the estimate downward or even flip its sign when effects grow, so run a Goodman-Bacon decomposition or switch to a modern estimator; flat pre-trends are low-power supporting evidence rather than proof, and proceeding only after a pre-test passes distorts the subsequent estimates, so use HonestDiD to report the departure needed to overturn the result; throwing in a lagged dependent variable to "absorb the trend" brings the Nickell bias; cluster at the level of treatment assignment; few treated clusters are the most dangerous case.
- Non-binary DiD (Chapter 20): a changing treatment rate in the control group is the main threat to fuzzy DiD, so first report the treatment rate in every cell and test A6, and when it fails, report only the TC or CIC bounds; when Wald-DID diverges systematically from Wald-TC and Wald-CIC, the most likely cause is a changing effect for the incumbent treated. In a dose design, do not read the slope as a causal response: that reading needs strong parallel trends, and pre-trends have no power against its failure; with no zero-dose units, run the quasi-untreated test first. A long-horizon $DID_\ell$ needs parallel trends over a longer window, and an insignificant placebo may only mean low power.
- RDD (Chapter 21): the number-one threat is manipulation and sorting, so use `rddensity` to look for bunching in the density and check that predetermined covariates pass smoothly through the threshold; when there is heaping at the cutoff, use a donut and report a range of radii; report how the estimate moves with the bandwidth, do not use global polynomials of order three or higher, and use local quadratic fits by default for kink designs; when the score takes only a few values, clustering by score value makes things worse, so report honest intervals with a sensitivity analysis over $M$, or use local randomization, and remember that $M$ cannot be estimated from the data; when the threshold triggers several treatments at once (compound treatment), the jump is their combined effect, and only institutional detail can rule this out; pooling several cutoffs requires stating the weights, and extrapolation rests on an untestable constant bias.
- Synthetic control (Chapter 22): when the pre-period fit is poor, do not force an interpretation of the gap; switch to augmented SC and report it alongside SDID and SC. A convex-hull violation is an extrapolation problem: when the treated unit lies outside the donors' convex hull, no convex combination can reproduce it. Interpolation bias is a different matter: even when the treated unit is inside the hull, donors far from it interpolate a possibly nonlinear response surface over a wide range and create bias, and the remedy is to keep only similar donors or to use penalized SC (Chapter 22 §3.4). "Non-unique weights" covers different cases: when the predictors lie outside the hull, the solution sits on a face of the hull's boundary and the weights are sparse, while inside the hull the weights that fit exactly form a whole set and the optimizer returns just one point of it (Chapter 22 §4.1); the dense weights of `tidysynth` come from the solver not reaching the exact optimum, not from non-uniqueness (Chapter 22 §4.3); either way, the weights need not recover the donors that generated the data and must not be read as "the true control group". Keep the donor pool clean and remove donors touched by spillover. With a single treated unit, the resolution of permutation inference is capped by the number of donors, but the conformal test and the `scpi` prediction intervals are also available (Chapter 22 §3.7); the conformal p-value depends on what you assume about the time structure of the shocks, so report several versions.
- Sensitivity / partial ID (Chapter 23): coefficient stability puts no constraint on unobserved confounding; the robustness value has no universal size threshold and must be compared against observed benchmarks; a negative-control failure must not be buried in the appendix; the identified set and the confidence interval express assumption uncertainty and sampling uncertainty respectively.
- DML / causal forest (Chapter 24): DML is not a master key: all its guarantees rest on the product-rate condition, and when the sample is too small or the nuisance too complex, the point estimate is biased and the interval undercovers; when unconfoundedness and overlap fail, ML cannot save you either; ML will fit heterogeneity that does not exist, so validate through calibration, GATES, and RATE; the regret guarantee of policy learning is only relative to a restricted policy class.
- Digital experiments (Chapter 25): peeking (monitoring continuously and stopping the moment the result is significant) poisons everything downstream, so use always-valid or group-sequential boundaries; interference is subtle (it does not violate randomization and the p-value looks fine, yet the experiment precisely estimates a quantity you did not want), and the fix is to change the design (cluster, switchback) rather than to test; use an SRM chi-square test as the first guardrail; handle multiple comparisons with Bonferroni or FDR and pre-register; ex-post observed power carries no information.
- Interference / marketplace experiments (Chapter 26): an individual A/B test may correctly estimate the direct effect yet cannot recover the rollout effect; cluster leakage, carryover, and the wrong level of inference are the three big threats; the short-run experimental effect does not automatically equal the long-run equilibrium effect.

### Part IV Structural Methods

- Discrete choice (Chapter 27): logit coefficients are identified only up to a scale normalization and cannot be compared in magnitude across models or samples; what is comparable are ratios and elasticities; wrongly imposing IIA does not raise an error automatically, and under separation a finite MLE may not even exist; getting the nests wrong in nested logit sends the substitution pattern in the wrong direction, and forgetting to instrument the endogenous within-nest share pushes $\lambda$ toward 0, faking a strong within-nest correlation that does not exist (Chapter 28 §4.2); the variance terms of mixed logit are identified by rich choice-set variation; a panel can strengthen this identification substantially but is not necessary; the second-step standard errors of a control function must account for the first step; standard SML retains simulation bias at a fixed number of draws $R$.
- BLP demand (Chapter 28): the validity of the price instrument is the foundation (a cost shifter must enter only through cost, and BLP instruments degenerate into weak instruments in large markets); random coefficients that collapse toward zero quietly pull the model back to IIA, so report the standard errors of $\sigma$ and the instrument strength; too loose a contraction tolerance lets error seep into the outer GMM (MPEC sidesteps this problem); zero shares leave the log share undefined, and discarding them casually introduces selection bias.
- Supply side / mergers (Chapter 29): conduct cannot be identified from the first-order conditions themselves: two conduct assumptions that both yield positive costs cannot be told apart from costs alone, so report cost plausibility together with a Rivers-Vuong test, and do not default to Bertrand; a misspecified demand (logit IIA understating near-neighbor diversion) propagates through the FOC and is systematically amplified in markups and merger price increases; an efficiency defense passes through to price only if marginal cost falls; merger simulation is static and short-run, and its extrapolation is limited by the Lucas critique.
- Entry / partial identification (Chapter 30): multiple equilibria break point identification, and forcing in an equilibrium-selection mechanism makes the estimate drift with the assumption, so be open about this dependence or switch to partial identification; the information structure (complete versus incomplete information) must be argued from institutional detail and is not a technical option to switch at will; the identified set can be so wide that it carries no information, yet its width is itself informative.
- Dynamic structural models (Chapter 31): the discount factor $\beta$ is generally not identified, and assuming that fixing it is harmless is a common form of overconfidence, so report how the conclusion changes at several values of $\beta$ or identify it formally with an exclusion restriction; conditional independence rules out persistent unobserved heterogeneity, and a violation biases the estimates, so model types explicitly (finite mixtures); the validity of counterfactuals is exactly what the Lucas critique questions; the CCP first step is noisy on rarely visited states, so smooth the estimates or gather more data.

### Part V Topics in the Digital Economy

- Platforms / switching costs (Chapter 32): the bias in estimated switching costs is unambiguously upward: as long as persistent heterogeneity is left uncontrolled, they are almost never understated, so when you see a large switching cost, first ask "how much preference has been mistaken for lock-in?"; the initial-conditions problem is especially severe in a short panel, so when the panel is not long, default to Wooldridge's CRE approach; peer and network effects face the reflection problem, which must be broken using network structure.
- Search / advertising (Chapter 33): in a search market, a full-information model's price coefficient is attenuated systematically toward zero, and a large sample cannot save it; when the search logs show shallow search, it is time to change the model; getting the sequential versus fixed-sample protocol wrong makes the likelihood wrong, so estimate both and test between them with the prevalence of return visits; the exogeneity of the ranking cannot be taken for granted and needs a ranking experiment or the quasi-random perturbation from an algorithm change; advertising effects are the most easily exaggerated by selection, and an ROI based on observational correlations should be discounted heavily.
- AI / technology adoption / algorithms / data (Chapter 34): machine-prediction error is often non-classical, and the gold subsample must represent the target population; AI access, actual use, and workflow redesign are different treatment margins; an interaction with observed use usually suffers from post-treatment selection, and organizational complementarity is best studied with factorial assignment; short-run productivity does not equal the long-run payoff in skills and along the J-curve; algorithmic collusion appearing in simulation is not collusion in reality, and establishing the latter requires field evidence.

Taken together, these pitfalls point to a single lesson that runs through the whole book: the credibility of a method lies not in the sophistication of the regression specification or the algorithm, but in whether the identifying assumptions are laid on the table and in how the conclusion moves when those assumptions are loosened. A reference sheet can help you recall the formulas, the commands, and the common traps, but stating the assumptions clearly and reporting sensitivity honestly always sends you back to the main text and to your own judgment about the institutional background.
