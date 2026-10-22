---
title: "Binary & Multinomial Response Models"
subtitle: "From Latent Utility to Marginal Effects"
seriesline: "Foundations of Information Systems Economics · Chapter 12"
series: "Foundations of Information Systems Economics: Empirical Methods for Information Systems Economics"
sidebartitle: "Foundations of Information Systems Economics"
sidebarsub: "Theory · Methods · Applications"
chapterlabel: "Chapter 12 · Binary & Multinomial Response Models"
---

## Introduction

Helix is a fictional B2B developer platform. Its data scientist wants to know whether client firms with higher engagement are more likely to turn on the paid Pro tier, so she runs OLS with "adopted or not" as the outcome. The software reports no error, and the coefficients even come with standard errors, but when she prints the fitted values firm by firm, one firm's "adoption probability" is 121% and another's is -48%. Stranger still, the slope of the same regression almost exactly equals the average effect of engagement on the adoption probability. Which number is wrong, which is right, and why?

This chapter covers models whose outcome takes only a few values: binary response (adopt or not), ordered response (satisfaction on a scale of 1 to 5 stars), and multinomial response (three plans: Free, Pro, and Enterprise). They share one latent-variable framework: a firm has an unobserved net utility and adopts only when it crosses a threshold. This gives the chapter its main thread: the scale of the coefficients $\beta$ is set by the variance of the latent error, and the coefficients themselves rarely carry economic meaning; what can be reported and compared across models is the average partial effect of a regressor on the probability. On the estimation side, all of these models are conditional MLEs of the kind studied in Chapter 6, and their consistency and asymptotic normality follow directly from Theorem 1 in Section 3.1 of Chapter 6 and Theorem 3 in Section 3.3 of Chapter 6; the more basic framework of M-estimation is in Sections 3.1 to 3.4 of Chapter 5. This chapter is about the new problems that nonlinearity creates: the coefficient of an interaction term is not the interaction effect; omitting a variable that is independent of the regressors attenuates the coefficients but leaves the average partial effect unchanged; endogenous regressors call for a control function; individual effects in panels raise the incidental parameters problem; and multinomial choice comes with a strong restriction, IIA (independence of irrelevant alternatives).

The chapter proceeds as follows. Section 2 writes down Helix's adoption model and the estimands of the seven data blocks. Section 3, the most detailed in the chapter, takes up in turn the LPM and the index model, MLE, average partial effects and interaction terms, neglected heterogeneity, endogenous regressors, weaker assumptions, panel binary models, and ordered and multinomial response. Section 4 covers computation and software, Section 5 examines three anchoring papers, Section 6 works through the Helix data in full, and Section 7 discusses failure modes. The equivalence of the control function and 2SLS in linear models is Theorem 4 in Section 3.4 of Chapter 9, the Mundlak form of CRE is in Section 3.5 of Chapter 10, and Chapter 27 will extend this chapter's random utility model to structural discrete choice.

## 1 A Probability Printed as 121%

::: {.case}
Helix sells cloud development tools to a large base of client firms. Each firm $i$ has a standardized engagement measure $X_i$ (constructed from past API usage, with mean zero and standard deviation one), and the platform cares about how it relates to whether the firm turns on the Pro tier, $Y_i \in \{0, 1\}$. The data are a cross section of 4,000 firms, with an adoption rate of 42.75%. The most natural first step is to treat adoption as an ordinary outcome variable and run OLS, that is, the LPM: $Y_i = \gamma_0 + \gamma_1 X_i + \varepsilon_i$. The estimate is $\hat\gamma_1 = 0.215$ (standard error 0.007), with an intercept of 0.427. Read literally, a one-standard-deviation increase in engagement raises the adoption probability by about 21.5 percentage points.
:::

The trouble lies in the fitted values. The highest predicted adoption probability is 1.213 and the lowest is -0.478; 2.93% of firms have a predicted "probability" outside $[0, 1]$, 2.5% below zero and 0.4% above one. This is not sampling error: in the population, 2.91% of the fitted values of the same linear projection fall out of bounds. A straight line with no ceiling and no floor cannot hold a quantity that must stay in $[0, 1]$, the conditional mean $\Pr(Y_i = 1 \mid X_i)$. It also locks the marginal effect of engagement at a constant, yet a firm already set on adopting, or on not adopting, will not be moved by a little more engagement; the firms that actually move are the ones wavering near the threshold.

![Helix's 4,000 firms. The gray dots are the adoption shares in 40 equal-sized groups formed by engagement; the orange solid line is the LPM's fitted line, the orange dots at its two ends mark the smallest and largest fitted values, -0.478 and 1.213, and the light orange bands are the regions outside $[0, 1]$; the blue solid line is the fitted probability of the logit MLE, and the gray dashed line is the true $\Pr(Y = 1 \mid X)$.](assets/fig/fig_12_lpm_logit.svg)

Yet the LPM slope itself is not bad. Estimated with a correctly specified logit (the data were in fact generated by a logit), the average partial effect of engagement on the probability is 0.2171, against a true value of 0.216, and the LPM's 0.2153 differs from it by only 0.0018. The coefficients of the three models, however, differ fivefold: the logit slope is 1.097, the probit slope 0.659, and the LPM slope 0.215. So the first question has two halves: why the LPM prints out-of-bounds probabilities and yet estimates the average effect almost exactly, and why three coefficients that look wildly different say the same thing. Both answers are in Section 3: the scale of the coefficients is set by the variance of the latent error, and it is the average partial effect that is comparable; the LPM estimates the best linear approximation to the conditional mean, and in a design like Helix's, where the regressor is normally distributed, the slope of this approximation exactly equals the average partial effect.

## 2 The Economic Model and the Estimand

### 2.1 Latent Variables and Random Utility

We can view the adoption decision through a random utility model. The net benefit firm $i$ obtains from adopting Pro is a latent variable, made up of an observable part and a part the platform cannot see:

$$Y_i^* = X_i'\beta + \varepsilon_i, \qquad Y_i = \mathbf{1}\{Y_i^* > 0\}.$$

In general $X_i$ is a $k \times 1$ vector of regressors that includes a constant; Helix has a single regressor, so $X_i'\beta$ is $\beta_0 + \beta_1 X_i$, and the scalar $X_i$ is engagement itself. The term $\varepsilon_i$ collects what the platform cannot observe: the procurement process, how tight the budget is, the decision maker's preferences, and so on. If $\varepsilon_i$ is independent of $X_i$, with a distribution symmetric about zero and distribution function $F$, then

$$\Pr(Y_i = 1 \mid X_i) = \Pr(\varepsilon_i > -X_i'\beta \mid X_i) = F(X_i'\beta).$$

This is the index model: the linear index $X_i'\beta$ is the only channel through which the probability depends on $X_i$, and $F$ squeezes the whole real line into $(0, 1)$. Taking $F$ to be the standard normal distribution function $\Phi$ gives probit, and taking the logistic distribution function $\Lambda(z) = e^z/(1 + e^z)$ gives logit. Wooldridge writes the link as $G$; in this book $G$ already denotes cohorts and cluster counts, so we write $F$ and its density $F'$.

::: {.warning}
The size of $\beta$ cannot be separated from the scale of $\varepsilon_i$. Double the utility, the threshold, and the shock all at once, and the adoption probability does not change at all, so the data cannot tell them apart; only $\beta$ divided by the standard deviation of the error is identified. The convention is to fix the scale of the error: the probit variance is set to 1, and the logit variance is $\pi^2/3 \approx 3.29$, a standard deviation of 1.81. On the same data, logit coefficients are therefore roughly 1.6 to 1.8 times probit coefficients: 1.664 times in the Helix sample; fitting a probit to probabilities generated by a logit gives the population pseudo-true value $(-0.2400, 0.6589)$, a slope ratio of 1.67. A latent variable usually has no natural unit of measurement, and Wooldridge (2010) p. 566 says that, except in special cases, the magnitude of the coefficients is not especially meaningful. Comparing coefficients across models or across groups compares their scales.
:::

::: {.intuition}
Think of the coefficients as measuring the same distance with different rulers. On the probit ruler one standard deviation of the error is 1 mark; on the logit ruler it is 1.81 marks, the standard deviation of the logistic error; the LPM measures the probability directly. The same distance, "engagement pushes up the adoption probability," measures 1.097 marks for logit, 0.659 marks for probit, and 0.215 for the LPM. Comparing the numbers of marks is meaningless; to compare, we have to convert the distance back into a common unit, the unit of probability, and that is the APE (average partial effect).
:::

### 2.2 Two Layers of Estimands and Seven Data Blocks

The first layer of estimands is the model parameter $\beta$, whose true values in Helix are $\beta_0 = -0.40$ and $\beta_1 = 1.10$. The second layer is the partial effect read off from the parameters. Write $P(x) = \Pr(Y_i = 1 \mid X_i = x)$; the derivative of the probability with respect to a continuous regressor $x_j$ is

$$\frac{\partial P(x)}{\partial x_j} = F'(x'\beta)\,\beta_j,$$

which varies with $x$, is largest near a probability of 0.5, and tends to zero at both ends. There are two conventions for collapsing it into one number: the APE averages each firm's partial effect, $\mathrm{APE}_j = \mathbb{E}[F'(X_i'\beta)]\,\beta_j$; the PEA (partial effect at the average) computes it only once, at the mean of the regressors, $F'(\mathbb{E}[X_i]'\beta)\,\beta_j$. In Helix the true APE is 0.216 and the PEA is 0.264: the average of a nonlinear function is not the function evaluated at the average. The target of this chapter is the APE, and `ame` in the inherited code refers to the same quantity.

Each of Helix's other six data blocks corresponds to one theoretical point of Section 3. The table below lists their designs and estimands; the true values are computed by numerical integration or in closed form (Section 6.1).

| Block | Question | Design | Estimand and true value |
|---|---|---|---|
| B2 | Interaction terms | 4,000 firms in the Helix Deploy beta; two z-standardized usage measures, $x_1$ (API engagement) and $x_2$ (deployment pipeline), with correlation 0.3; adoption probability $\Phi(1 + x_1 + 0.7x_2 + 0.2x_1x_2)$, so the interaction coefficient is positive and the two features are complements in latent utility | Average cross-partial effect $\mathbb{E}[\partial^2 P(X_i)/\partial x_1\,\partial x_2]$, true value -0.0194 |
| B3 | Endogenous regressor | 4,000 firms in a 30-day trial; $X_i$ is the log of monthly cloud spending relative to USD 5,000, exogenous; $D_i$ is engagement during the trial, endogenous, because projects in a rush to launch raise both engagement and the value of Pro; the platform sends an onboarding reminder $Z_i$ to a random half of the firms, and it affects adoption only through $D_i$ | APE of $D_i$, true value 0.1475 |
| B4 | Neglected heterogeneity | Adoption also depends on budget flexibility $a_i$, which is independent of engagement and unseen by the platform; structural coefficient 0.9 | APE, true value 0.2127 |
| B5 | Panel | A panel of 1,000 firms over 8 months, $\Pr(Y_{it} = 1 \mid X_{it}, \alpha_i) = \Lambda(\alpha_i + 0.8X_{it})$, with the firm effect $\alpha_i$ positively correlated with the firm's average engagement | $\beta = 0.8$ and the APE over all firm-months, true value 0.1245 (the exact value from two-dimensional Gauss-Hermite quadrature, the same for the subpanel of the first 4 months; 0.1214 for the 4-month DGP) |
| B6 | Ordered response | Quarterly satisfaction ratings of 1 to 5 stars from 3,000 firms; firms with an account manager (managed accounts) have a latent satisfaction mean higher by 0.2 and a standard deviation of 1.5, against 1 for self-serve firms | APE of engagement on the probability of "giving 5 stars," true value 0.1008 |
| B7 | Multinomial response | 5,000 firms each choose one plan among Free, Pro, and Enterprise, with engagement and the standardized log number of seats as regressors; one version satisfies IIA, and the other puts Pro and Enterprise in the same nest | Coefficients of the MNL (multinomial logit) and APEs of engagement on the probability of each plan; the engagement coefficient for Enterprise is positive, yet its APE on the Enterprise probability is -0.0407 |

The estimands of the seven blocks have one thing in common: apart from the $\beta$ of B5, they are all quantities on the probability scale. Section 3 explains why the coefficients themselves are often not the quantities we want.

## 3 Identification and Large-Sample Theory

This section follows the format of Chapters 5 and 6. It first uses Helix's numbers or figures to show what is to be proved, then states numbered assumptions and theorems, followed by proof sketches of two to six steps; each step names the assumptions it uses, and regularity conditions are only cited to their source. After each theorem, an "In words" paragraph answers four questions: under which assumptions it estimates which quantity, what it estimates when the specification is wrong, who sets the weights, and what Helix's numbers say. The parts of Wooldridge (2010) on binary and multinomial response have no numbered theorems, so the theorem boxes cite its equation numbers and pages, and the theorem numbering is this chapter's own. Sections 3.1 to 3.4 cover exogenous regressors, Sections 3.5 and 3.6 endogenous regressors and weaker assumptions, Sections 3.7 and 3.8 panels, and Section 3.9 ordered and multinomial response.

### 3.1 The LPM and the Index Model

Start with the population. Helix's engagement is standard normal. Computed by numerical integration, the population slope of the linear projection of $Y_i$ on $(1, X_i)$ is 0.2156, and the population APE of engagement on the probability is also 0.2156; Theorem 1 (ii) below shows that the two are exactly equal in this design, so the agreement to four decimal places is no coincidence. Yet 2.91% of the fitted values on this projection line still fall out of bounds. The LPM is both right and wrong, and the theorem below says where it is right.

::: {.assumption}
**LP.1 (Wooldridge 2010 equation (15.4), p. 562)** $\Pr(Y_i = 1 \mid X_i) = X_i'\gamma$.

**LP.2** $\{(Y_i, X_i)\}$ is a random sample, and $\mathbb{E}[X_iX_i']$ has full rank.
:::

LP.1 requires the probability to be linear over the entire support of the regressors, and it cannot hold when a regressor is unbounded and its slope is nonzero: as $x_j$ keeps increasing, the linear "probability" eventually crosses 0 or 1 (Wooldridge 2010 p. 563). LP.2 is just the condition for OLS to be defined.

::: {.theorem}
**Theorem 1 (the LPM; Wooldridge 2010 §15.2, equations (15.4) to (15.6), pp. 562-563; the result for normal regressors is in §15.6, p. 579, restated on p. 604)** Under LP.2, OLS converges to the coefficient of the linear projection of $Y_i$ on $X_i$, $\gamma_* = (\mathbb{E}[X_iX_i'])^{-1}\mathbb{E}[X_iY_i]$. (i) If LP.1 also holds, $\gamma_* = \gamma$ and OLS is unbiased and consistent, but $\mathrm{Var}(Y_i \mid X_i) = X_i'\gamma\,(1 - X_i'\gamma)$, so heteroskedasticity-robust standard errors are needed. (ii) If the regressors other than the constant are multivariate normal (so there are no squares, interactions, or discrete regressors), then, whatever the shape of $P(x)$, as long as it is differentiable and the expectations of its partial derivatives are finite, each slope in $\gamma_*$ equals the APE of the corresponding regressor.
:::

The proof sketch has three steps. Step one (LP.2): the law of large numbers gives the probability limit $\gamma_*$ of OLS, and this step does not need the model to be correct. Step two (LP.1): when the conditional mean is linear, the projection is the conditional mean, and the error has mean zero given $X_i$, which gives unbiasedness; the variance of a Bernoulli variable is equation (15.6). Step three (normality): Stein's lemma says that when $X_i$ is normal and $P$ is differentiable with a derivative of finite expectation, $\mathrm{Cov}(X_i, P(X_i)) = \mathrm{Var}(X_i)\,\mathbb{E}[P'(X_i)]$ (in the multivariate case, replace the variance by the covariance matrix and the derivative by the gradient). By iterated expectations, $\mathrm{Cov}(X_i, Y_i) = \mathrm{Cov}(X_i, P(X_i))$, so

$$\frac{\mathrm{Cov}(X_i, Y_i)}{\mathrm{Var}(X_i)} = \mathbb{E}[P'(X_i)] = \mathrm{APE}.$$

In words: under LP.2, the LPM estimates the best linear approximation to the conditional mean, $\gamma_*$, not a parameter of any probability model; only when LP.1 holds is it "the effect of $X$ on the probability." When the specification is wrong, for example when the true probability is a logit, OLS still converges to the projection, and the out-of-bounds fitted values are simply the price of approximating an S-shaped curve by a straight line. Who sets the weights: with a single regressor, write $P(X_i)$ as $P(-\infty) + \int P'(x)\mathbf{1}\{X_i > x\}\,dx$, substitute it into $\mathrm{Cov}(X_i, P(X_i))$, and interchange the order of integration; this gives the projection slope $= \int P'(x)\,\omega(x)\,dx$, where the weight $\omega(x) = \mathbb{E}[(X_i - \mathbb{E}X_i)\mathbf{1}\{X_i > x\}]/\mathrm{Var}(X_i)$ is nonnegative, integrates to one, and is determined by the distribution of $X_i$. When $X_i$ is normal (with any mean and variance), $\mathbb{E}[(X_i - \mathbb{E}X_i)\mathbf{1}\{X_i > x\}]$ equals $\mathrm{Var}(X_i)$ times the density of $X_i$ at $x$, so $\omega$ is exactly this density and the slope is the APE. Helix has this design: in the sample the LPM slope is 0.2153 and the logit APE is 0.2171, a difference of 0.0018. When $X_i$ is skewed or takes discrete values, $\omega$ departs from the density, and the LPM slope is only an approximation to the APE; Wooldridge (2010) p. 563 notes that the approximation is often good but unreliable at extreme values of the regressors. Horrace and Oaxaca (2006) give sufficient conditions for OLS to be unbiased and consistent under the LPM and, on that basis, suggest a trimming estimator that may reduce the bias. In substance, the condition can be understood through the statement in Wooldridge (2010) §15.2: if equation (15.4) holds over the entire support of the regressors, that is, if the linear index $X_i'\gamma$ of every observation lies in $[0, 1]$, then OLS is unbiased and consistent (p. 562); once some observations fall out of bounds, LP.1 already fails, and OLS estimates only the projection $\gamma_*$ (p. 563).

The index model repairs the LPM's flaws with a monotone $F$: the probability always stays in $(0, 1)$, and the partial effect $F'(x'\beta)\beta_j$ varies with $x$. The price is that the coefficients are identified only up to scale and must be converted back into probabilities for reporting.

### 3.2 MLE: The Likelihood, the Score, and Two Sets of Standard Errors

On Helix, a hand-coded Newton iteration started at zero converges in 6 steps to $(-0.367, 1.097)$, against the truth $(-0.40, 1.10)$, and differs from `glm` on the order of $10^{-13}$. The two sets of standard errors for the slope are almost the same: 0.0442 from the information matrix and 0.0439 from the sandwich, a ratio of 0.993. On the Poisson PML of Chapter 5 the two differed by nearly a factor of two; why do they agree here? Section 3.3 of Chapter 6 gave the general answer, and binary models have a reason of their own.

::: {.assumption}
**BR.1 (Wooldridge 2010 equations (15.8) and (15.9), p. 565; p. 568)** $Y_i = \mathbf{1}\{X_i'\beta + \varepsilon_i > 0\}$, where $\varepsilon_i$ is independent of $X_i$, and its distribution function $F$ is known, symmetric about zero, and twice continuously differentiable, and for all $z$, $0 < F(z) < 1$ and $F'(z) > 0$.

**BR.2 (Wooldridge 2010 §15.4, p. 568, via Theorems 13.1 and 13.2)** $\{(Y_i, X_i)\}$ is a random sample, $\mathbb{E}[X_iX_i']$ has full rank, and ML.2 and ML.3 of Chapter 6 hold for $q = -\ell_i$, with the log-likelihood given below.
:::

BR.1 specifies the entire conditional distribution: $Y_i$ given $X_i$ is Bernoulli with success probability $F(X_i'\beta)$, and the log-likelihood of a single observation is

$$\ell_i(\beta) = Y_i \ln F(X_i'\beta) + (1 - Y_i)\ln\{1 - F(X_i'\beta)\}.$$

The full-rank condition of BR.2 rules out perfect collinearity; a strictly monotone $F$ guarantees that different values of $\beta$ give different probabilities, so the KL argument of Section 3.2 of Chapter 5 delivers identification.

::: {.theorem}
**Theorem 2 (MLE for binary response; Wooldridge 2010 §15.4, equations (15.16) to (15.20), pp. 568-569; the general theory is Theorems 1 and 3 of Chapter 6, Wooldridge 2010 Theorems 13.1 and 13.2)** Under BR.1 and BR.2, the maximizer of $\sum_i \ell_i(\beta)$, denoted $\hat\beta$, is consistent, and $\sqrt n(\hat\beta - \beta) \xrightarrow{d} \mathcal{N}(0, A_0^{-1})$, where

$$A_0 = \mathbb{E}\left[\frac{F'(X_i'\beta)^2\,X_iX_i'}{F(X_i'\beta)\{1 - F(X_i'\beta)\}}\right],$$

and the information matrix equality holds, so the inverse of the information matrix and the sandwich are asymptotically the same. For logit, the score of the log-likelihood is $X_i\{Y_i - \Lambda(X_i'\beta)\}$, the Hessian $-\sum_i \Lambda(X_i'\beta)\{1 - \Lambda(X_i'\beta)\}X_iX_i'$ does not contain $Y_i$, and the log-likelihood is globally concave.
:::

The proof sketch has three steps. Step one (BR.1 and BR.2): the Bernoulli density of equation (15.16) is correctly specified and identified, and Theorem 1 of Chapter 6 gives consistency. Step two (BR.1 and BR.2): the score of the log-likelihood is $F'(X_i'\beta)X_i\{Y_i - F(X_i'\beta)\}/[F\{1 - F\}]$ (equation (15.18); it differs from the $s_i$ of Chapters 5 and 6 only by a sign), and because $\mathbb{E}[Y_i \mid X_i] = F(X_i'\beta)$, its mean given $X_i$ is zero; Theorem 3 of Chapter 6 gives asymptotic normality and the information matrix equality, and the expected Hessian is equation (15.19). Step three (BR.1, logit): $\Lambda'(z) = \Lambda(z)\{1 - \Lambda(z)\}$ cancels against the denominator, so the score becomes the residual times the regressors; the Hessian is a negative weighted sum of $X_iX_i'$ and is negative definite under full rank.

In words: under BR.1 and BR.2, MLE estimates the true parameter $\beta$ of the index model, but the scale of this parameter is set by the normalization of the error, and it is not itself an effect on the probability. When the specification is wrong, for example when a probit is fitted to data generated by a logit, MLE estimates the QMLE pseudo-true value of Section 3.6 of Chapter 6: in Helix's population the pseudo-true probit slope is 0.6589, but the APE it implies is 0.2151, very close to the true 0.2156. Who sets the weights: the first-order condition multiplies each residual $Y_i - F$ by $F'/[F(1 - F)]$; under logit this weight is identically one, and the first-order condition is just $\sum_i X_i\{Y_i - \Lambda(X_i'\hat\beta)\} = 0$, which are the "residual times regressor" moments of Section 3.3 of Chapter 7, so the logit MLE is the exactly identified GMM estimator with these moment conditions. Why the two sets of standard errors agree: once $\Pr(Y_i = 1 \mid X_i)$ is specified for a binary outcome, the whole conditional distribution is specified, and there is no room for "the mean right but the variance wrong" (Wooldridge 2010 p. 569). So, unlike NLS, in a binary model the gap between the information-matrix standard error and the sandwich is itself a signal of misspecification; Helix's 0.0442 and 0.0439 differ only by sampling error, the same agreement that Section 3.3 of Chapter 6 anticipated. Tests follow the trinity of Theorem 5 of Chapter 6 (Wooldridge 2010 §15.5), and Section 3.3 uses the LR test for the interaction term. §15.5 also gives LM tests that require estimation only under the null: for example, testing $\delta = 0$ in the heteroskedastic probit of Section 3.6 needs only an auxiliary regression after the ordinary probit (equations (15.24) to (15.28), pp. 571-573).

### 3.3 APE, PEA, and Interaction Terms

In the probit of block B2, the coefficient of the interaction term is 0.221 (standard error 0.042), $z = 5.24$, and the LR statistic for dropping the interaction is 27.56: highly significant. It seems natural to read this as "the two features are complements, and firms that use both gain more in adoption probability." Yet the average cross-partial effect is -0.0161, with a 95% interval of $[-0.0272, -0.0051]$ and a true value of -0.0194: on average, the two features are substitutes on the probability scale; the cross-partial effect is positive for 33.7% of firms and negative for the rest.

::: {.theorem}
**Theorem 3 (APEs and cross-partial effects; Wooldridge 2010 §15.6, equations (15.29) to (15.33), pp. 575-578; for the cross-partial effect see Ai and Norton 2003, with the formula and its four corollaries following Norton, Wang and Ai 2004, pp. 155-156)** Under BR.1 and BR.2: (i) $\widehat{\mathrm{APE}}_j = n^{-1}\sum_i F'(X_i'\hat\beta)\hat\beta_j \xrightarrow{p} \mathrm{APE}_j$, and it is asymptotically normal, with a variance that is the sum of the part transmitted from $\hat\beta$ through the Delta method and the part due to the sampling fluctuation of $X_i$; (ii) if the index contains an interaction term $\beta_{jj'}x_jx_{j'}$, then

$$\frac{\partial^2 P(x)}{\partial x_j\,\partial x_{j'}} = \beta_{jj'}F'(x'\beta) + (\beta_j + \beta_{jj'}x_{j'})(\beta_{j'} + \beta_{jj'}x_j)\,F''(x'\beta),$$

which is generally nonzero when $\beta_{jj'} = 0$; its sign can change with $x$ and can be opposite to that of $\beta_{jj'}$; and its significance cannot be judged by testing $\beta_{jj'}$ with a $z$ statistic.
:::

The proof sketch has three steps. Step one (Theorem 2): $\hat\beta$ is consistent and $F'$ is continuous, so by the uniform law of large numbers (Theorem 1 of Chapter 5) the sample average converges to $\mathbb{E}[F'(X_i'\beta)]\beta_j$ (equation (15.32)). Step two (BR.1 and BR.2): split $\widehat{\mathrm{APE}}_j - \mathrm{APE}_j$ into two terms: one is the error that, given the sample $X_i$, comes from $\hat\beta$ and is handled by the Delta method of Section 3.2 of Chapter 8 (equation (15.30)), and the other is the fluctuation of $n^{-1}\sum_i F'(X_i'\beta)\beta_j$ around its population mean; the score has mean zero given $X_i$, so the two terms are asymptotically uncorrelated and their variances add. Step three (BR.1: $F$ twice differentiable): differentiating $(\beta_j + \beta_{jj'}x_{j'})F'(x'\beta)$ once more with respect to $x_{j'}$, the product rule gives (ii); when $\beta_{jj'} = 0$, what remains is $\beta_j\beta_{j'}F''(x'\beta)$.

In words: under BR.1 and BR.2, the APE estimates the population average effect on the probability scale; it does not depend on the scale normalization and can therefore be compared across models (Wooldridge 2010 p. 578): Helix's logit, probit, and LPM coefficients are 1.097, 0.659, and 0.215, and their APEs are 0.2171, 0.2169, and 0.2153. When $F$ is misspecified, the APE is far more robust than the coefficients; the APE implied by the probit pseudo-true value in the previous section is an example. Who sets the weights: the APE gives equal weight to every firm in the distribution of $X_i$, while the PEA looks only at one "average firm." In the Helix sample the PEA is 0.265, more than a fifth higher than the APE, because the average firm happens to sit near where $F'$ is largest; when the regressors include interactions or squares, the PEA is not even clearly defined (pp. 575-576). The interaction term is the extreme case of this logic: $\beta_{jj'}$ acts only on latent utility, while the cross-partial effect also depends on $F''$. The adoption rate in B2 is 0.727, so most firms have a positive index, where the normal density is decreasing; the second term is negative and outweighs the first, so the positive 0.221 turns into -0.0161. A common mistake is to average only the first term, $\beta_{jj'}\phi(X_i'\beta)$; on B2 this gives 0.0456, which gets even the sign wrong. The figure below shows the mechanism.

::: {.intuition}
The interaction coefficient speaks of complementarity in latent utility, the cross-partial effect of complementarity in probability, and between the two stands a ceiling. For a firm already nine-tenths sure to adopt, using one more feature can raise the probability only a little; when it uses both more, the first has already used up most of the room the second could push into. The closer the probability is to 1, the stronger this "room used up" effect, and that is what $F'' < 0$ means. In B2 more than seven in ten firms are already close to the ceiling; the little complementarity in latent utility cannot withstand the squeeze of the ceiling, and on the probability scale the features become substitutes. For firms with an adoption probability below one half the situation is reversed: they are near the floor and far from the ceiling, and the cross effect is positive, which is why almost all the points on the negative-index side of the figure lie above zero.
:::

![Block B2: the gray dots are the estimated cross-partial effects $\partial^2 P/\partial x_1\,\partial x_2$ for one firm in five (800 firms in all), plotted against the estimated probit index $X_i'\hat\beta$ on the horizontal axis; the blue solid line is the average cross-partial effect over the 4,000 firms, and the orange dashed line is the first term alone, $\hat\beta_{jj'}\,\phi(X_i'\hat\beta)$.](assets/fig/fig_12_cross_partial.svg)

### 3.4 Neglected Heterogeneity

In block B4, budget flexibility $a_i$ is independent of engagement and unseen by the platform. A probit on engagement alone gives a coefficient of 0.638, whose interval $[0.590, 0.685]$ lies far from the structural coefficient 0.9; adding $a_i$ as well (possible only in a simulation) gives 0.908. It looks as if omitting a variable unrelated to the regressors can throw a probit off too, contrary to the intuition from linear regression. Yet the APEs of the two probits are 0.2135 and 0.2144, against a true value of 0.2127.

::: {.intuition}
Picture each firm as a steep S-shaped probability curve, with budget flexibility $a_i$ shifting each firm's curve left or right. The platform cannot see $a_i$; it sees only the average of all the curves stacked together, and the averaged curve is flatter than any single one. The short probit fits this flatter population curve, so its coefficient naturally shrinks; but "the slope of the population-average curve at a point" is the same as "the average of the firms' slopes at that point," so the APE is unaffected. The attenuated coefficient is not a wrong estimate; its object has changed, from the curve of an individual firm to the curve of the population.
:::

::: {.assumption}
**NH.1 (Wooldridge 2010 equation (15.34), pp. 582-583)** $\Pr(Y_i = 1 \mid X_i, a_i) = \Phi(X_i'\beta + a_i)$, $a_i$ is independent of $X_i$, and $a_i \sim \mathcal{N}(0, \sigma_a^2)$.
:::

Independence is the key assumption, and it is much stronger than zero correlation or $\mathbb{E}[a_i \mid X_i] = 0$ (p. 583).

::: {.theorem}
**Theorem 4 (neglected heterogeneity; Wooldridge 2010 §15.7.1, equations (15.35) to (15.38), pp. 583-584)** Under NH.1: (i) $\Pr(Y_i = 1 \mid X_i) = \Phi\big(X_i'\beta/\sqrt{1 + \sigma_a^2}\big)$, and the probit of $Y_i$ on $X_i$ consistently estimates $\beta/\sqrt{1 + \sigma_a^2}$; (ii) the structural partial effect $\beta_j\phi(x'\beta + a)$, averaged over the distribution of $a_i$, equals $\big(\beta_j/\sqrt{1 + \sigma_a^2}\big)\,\phi\big(x'\beta/\sqrt{1 + \sigma_a^2}\big)$, so this probit consistently estimates the APE.
:::

The proof sketch has three steps. Step one (NH.1): $a_i + \varepsilon_i$ is independent of $X_i$ and distributed as $\mathcal{N}(0, 1 + \sigma_a^2)$, so the probability is $\Phi\big(X_i'\beta/\sqrt{1 + \sigma_a^2}\big)$ (equation (15.35)), a correctly specified probit. Step two (NH.1): by iterated expectations, $\mathbb{E}_a[\Phi(x'\beta + a_i)] = \Phi\big(x'\beta/\sqrt{1 + \sigma_a^2}\big)$. Step three (NH.1, $\phi$ bounded): differentiate with respect to $x_j$; the derivative and the expectation over $a_i$ can be interchanged, which is equation (15.38).

In words: under NH.1, the short probit estimates $\beta/\sqrt{1 + \sigma_a^2}$. In B4, $\sqrt{1 + \sigma_a^2} = 1.414$, the probability limit is 0.636, and the sample estimate is 0.638; the coefficients are attenuated, but their signs and the ratio of any two of them are unchanged, and the APE is consistently estimated too: 0.2135 against the true 0.2127. When the specification is wrong, in the sense that $a_i$ is correlated with $X_i$, for example because $\mathbb{E}[a_i \mid X_i]$ varies with $X_i$, the probability limit absorbs this correlation and even the APE comes out wrong (p. 585); the pooled estimators of Section 3.7 face this problem. Who sets the weights: the APE averages the partial effect over the distributions of both $X_i$ and $a_i$; the larger $\sigma_a$, the more the probability curve is smeared flat, and the further the structural coefficient lies from "the population-average probability curve." This gives a practical rule: when two groups of firms have different heterogeneity variances, comparing their probit coefficients compares two different values of $\sqrt{1 + \sigma_a^2}$, and only the APEs are comparable. In B4 engagement is normal, and the LPM slope, 0.2131, again lands on the APE, another instance of Theorem 1 (ii).

### 3.5 Endogenous Regressors: The Control Function and IV Probit

In block B3, a probit run directly on $(X_i, D_i)$ gives an APE of engagement of 0.241, with interval $[0.2325, 0.2496]$, 1.63 times the true 0.1475: firms in a rush to launch are both more active and keener on Pro. The onboarding reminder is randomized, and with a first-stage F statistic of 642.2 the instrument is strong.

::: {.assumption}
**CF.1 (Wooldridge 2010 equations (15.39) to (15.41), p. 585)** $Y_i = \mathbf{1}\{X_i'\beta + \beta_D D_i + \varepsilon_i > 0\}$ and $D_i = X_i'\pi_1 + Z_i'\pi_2 + v_i$; $(\varepsilon_i, v_i)$ is bivariate normal with mean zero and independent of $(X_i, Z_i)$, with $\mathrm{Var}(\varepsilon_i) = 1$ and $\mathrm{Var}(v_i) = \sigma_v^2$.

**CF.2** $\pi_2 \neq 0$, and the second-moment matrix of $(X_i, Z_i)$ has full rank.
:::

CF.1 requires $D_i$ to be normal given the exogenous variables, so $D_i$ should be a continuous variable (p. 585); $\mathrm{Var}(\varepsilon_i) = 1$ is the normalization that makes the parameters correspond directly to APEs (p. 586). CF.2 is the relevance condition that comes on top of the exclusion restriction, namely 2SLS.2 of Section 3.2 of Chapter 9, on which the consistency of 2SLS in Theorem 2 of Chapter 9 also rests.

::: {.intuition}
The first-stage residual $\hat v_i$ is a shadow of "how urgent the project is." The onboarding reminder is random, so the part of engagement that the reminder and cloud spending can explain is unrelated to urgency, and the residual that is left over carries the urgency. Put this shadow into the probit as an extra regressor, and the coefficient of engagement is left with only the differences among "firms that are equally urgent but differ in engagement"; that is the whole idea of the control function. The coefficient of $\hat v_i$ also answers, in passing, whether engagement is endogenous: if the shadow has no explanatory power, there is nothing to control for.
:::

::: {.theorem}
**Theorem 5 (the control function and IV probit; Wooldridge 2010 §15.7.2, equations (15.42) to (15.50), pp. 586-592; Rivers and Vuong 1988)** Under CF.1 and CF.2, let $\rho = \mathrm{Cov}(v_i, \varepsilon_i)/\sigma_v^2$ and $\sigma_e^2 = 1 - \rho^2\sigma_v^2$: (i) $\varepsilon_i = \rho v_i + e_i$, where $e_i \sim \mathcal{N}(0, \sigma_e^2)$ is independent of $(X_i, Z_i, v_i)$, so $\Pr(Y_i = 1 \mid X_i, Z_i, v_i) = \Phi\big((X_i'\beta + \beta_D D_i + \rho v_i)/\sigma_e\big)$; running a first-stage OLS of $D_i$ to obtain the residual $\hat v_i$, and then a probit on $(X_i, D_i, \hat v_i)$, consistently estimates $(\beta, \beta_D, \rho)/\sigma_e$; (ii) for $\hat v_i$ in the second-step probit, the usual $t$ statistic is a valid test of $H_0: \rho = 0$; (iii) the average structural function $\mathrm{ASF}(x, d) = \mathbb{E}_v\big[\Phi\big((x'\beta + \beta_D d + \rho v_i)/\sigma_e\big)\big] = \Phi(x'\beta + \beta_D d)$ is consistently estimated from the second-step coefficients and the sample average over $\hat v_i$, and the APE of $D_i$ is the average of its derivative with respect to $d$; (iv) the full MLE of the same model (IV probit) is given by the joint likelihood of equation (15.50), and under exact identification it coincides with the two-step estimator.
:::

The proof sketch has four steps. Step one (CF.1): under bivariate normality the conditional mean of $\varepsilon_i$ given $v_i$ is linear, and the remainder is normal and independent of $v_i$ (equations (15.42) to (15.44)). Step two (CF.1 and CF.2): replacing $v_i$ by the first-stage residual is the two-step M-estimation of Section 3.6 of Chapter 5; when $\rho \neq 0$, the first-step error enters the variance of the second step (the $\Gamma_0 \neq 0$ of Theorem 6 of Chapter 5), and in practice the two steps are bootstrapped together (p. 588). Step three ($H_0$): when $\rho = 0$, the coefficient of $\hat v_i$ is zero, the first-step error does not affect the limiting distribution of the second-step score, and the usual $t$ statistic is valid; nor does this step need $v_i$ to be normal (p. 587). Step four (CF.1): integrating over $e_i$ and then $v_i$ gives (iii) (equations (15.46) and (15.47)); writing the joint density as the density of $Y_i$ given $(D_i, X_i, Z_i)$ times the density of $D_i$ given $(X_i, Z_i)$ gives (iv) (equations (15.48) to (15.50)).

In words: under CF.1 and CF.2, the two-step method estimates the coefficients scaled up by a factor of $1/\sigma_e$, and the ASF converts the scale back, so the resulting APE is the average change in the adoption probability when "engagement is raised a little exogenously while the distribution of unobserved project urgency stays the same." On B3 the control-function APE is 0.1659, with a bootstrap standard error of 0.0148 and an interval $[0.1370, 0.1949]$ that covers 0.1475; the scaled coefficient of $D_i$ is 0.647, whose true value is $\beta_D/\sigma_e = 0.577$; for $\hat v_i$ the $t$ statistic is 7.11, a clear rejection of exogeneity. The IV probit MLE estimates $\beta_D$ directly on the structural scale and gives 0.583, with an interval $[0.438, 0.728]$ that covers 0.5, and an APE of 0.1660. Under exact identification, the two-step coefficients converted back to the structural scale are the same as the MLE (Wooldridge 2010 p. 590); the difference between the two APEs in the fourth decimal place comes only from the computation: the control function averages over the sample $\hat v_i$, while the MLE uses the closed form of the ASF under normality, $\Phi(x'\beta + \beta_D d)$. When the specification is wrong: if $\rho \neq 0$, the normality of $v_i$ is essential; if $D_i$ is binary, CF.1 cannot hold, as the next section discusses. Who sets the weights: the ASF averages over the marginal distribution of $v_i$, whereas 2SLS weights by the variation that the instruments generate; 2SLS on the LPM gives 0.1559, with an interval $[0.1259, 0.1860]$ that also covers the truth, the comparison that Wooldridge (2010) p. 594 recommends. In linear models the control function and 2SLS are the same estimator (Theorem 4 (i) in Section 3.4 of Chapter 9); in a probit, copying the 2SLS recipe and plugging in the first-stage fitted values yields coefficients that carry an unknown positive scale factor (which can be any positive number, whereas the control function's factor must exceed one); this route offers no simple test of exogeneity, and its coefficients cannot be compared directly with those of an ordinary probit (p. 589).

### 3.6 Binary Endogenous Regressors, Heteroskedasticity, and the Special Regressor

This section has no numerical block; it only marks the boundaries of the theory.

When the endogenous regressor is binary (for example, whether a firm attended an onboarding training session), the equation $D_i = X_i'\pi_1 + Z_i'\pi_2 + v_i$ of CF.1 cannot have a $v_i$ that is normal and independent of the instruments. The approach of Wooldridge (2010) §15.7.3 is the bivariate probit: two probit equations with bivariate normal errors, estimated by full MLE (equations (15.51) to (15.56)). Two tempting shortcuts both fail: plugging the fitted probabilities of a first-stage probit into the second step is what Wooldridge (2010) §9.5 calls the "forbidden regression," and the parameter estimates are inconsistent; adding the probit residuals to the second step gives only a valid test of exogeneity, not consistent estimates (p. 597).

When the latent error is heteroskedastic, for example $\mathrm{Var}(\varepsilon_i \mid X_i) = \exp(2X_i'\delta)$ (with the constant component of $\delta$ set to zero), the probability becomes $\Phi\big(X_i'\beta\,e^{-X_i'\delta}\big)$ (equation (15.24)), and the sign of the partial effect is determined by $\beta_j - \delta_j\,X_i'\beta$, so it can be opposite to that of $\beta_j$ (equation (15.57)). But if heteroskedasticity is understood as "the distribution of $\varepsilon_i$ changes with $X_i$" and one averages over $\varepsilon_i$, the partial effects of the ASF are proportional to $\beta_j$ (equations (15.58) to (15.60)). The two interpretations are observationally equivalent, and Wooldridge (2010) pp. 603-604 admits that there is no convincing way to choose between them.

Go one step further and leave the distribution of $\varepsilon_i$ unspecified: what can still be identified? The special regressor of Lewbel (2000) gives one answer. The numbering below follows its Boston College Working Paper No. 454 (revised August 1999).

::: {.assumption}
**SR.1 (Lewbel 2000 working paper, Assumptions A.1 to A.4, p. 6)** $Y_i = \mathbf{1}\{X^{\mathrm{sr}}_i + X_i'\beta + \varepsilon_i > 0\}$, where the superscript sr marks the special regressor, whose coefficient in the latent variable is normalized to 1; $X^{\mathrm{sr}}_i$ is continuous given $X_i$, with conditional density $f_{\mathrm{sr}}(\cdot \mid X_i)$; the distribution of $\varepsilon_i$ given $(X^{\mathrm{sr}}_i, X_i)$ depends only on $X_i$; the support of $X^{\mathrm{sr}}_i$ given $X_i$ is an interval $[L, K]$ containing 0 ($L < 0 < K$; if necessary, first subtract a constant from $X^{\mathrm{sr}}_i$), and it covers the support of $-X_i'\beta - \varepsilon_i$; $\mathbb{E}[X_i\varepsilon_i] = 0$, and $\mathbb{E}[X_iX_i']$ has full rank.
:::

::: {.theorem}
**Theorem 6 (the special regressor; Lewbel 2000 working paper, Theorem 1, equations (2.1) and (2.2), p. 7)** Under SR.1, let

$$Y^{\mathrm{sr}}_i = \frac{Y_i - \mathbf{1}\{X^{\mathrm{sr}}_i > 0\}}{f_{\mathrm{sr}}(X^{\mathrm{sr}}_i \mid X_i)},$$

then $\mathbb{E}[Y^{\mathrm{sr}}_i \mid X_i] = X_i'\beta + \mathbb{E}[\varepsilon_i \mid X_i]$, and hence $\beta = (\mathbb{E}[X_iX_i'])^{-1}\mathbb{E}[X_iY^{\mathrm{sr}}_i]$.
:::

The proof sketch has two steps. Step one (SR.1): given $(X_i, \varepsilon_i)$, take the expectation over $X^{\mathrm{sr}}_i$; dividing by the density turns the expectation into an integral, $\int [\mathbf{1}\{\tilde x > -X_i'\beta - \varepsilon_i\} - \mathbf{1}\{\tilde x > 0\}]\,d\tilde x = X_i'\beta + \varepsilon_i$, and a support that contains 0 and covers the support of $-X_i'\beta - \varepsilon_i$ guarantees that the range of integration is not truncated. Step two (SR.1): then take the expectation over $\varepsilon_i$ and use $\mathbb{E}[X_i\varepsilon_i] = 0$; what remains is an ordinary OLS.

In words: under SR.1, the OLS of the transformed outcome on $X_i$ estimates, in units of the coefficient of $X^{\mathrm{sr}}_i$, the vector $\beta$, allowing $\varepsilon_i$ to be arbitrarily heteroskedastic in $X_i$, and Theorem 1' of the working paper extends it to instrumental variables. The price is equally clear. When the specification is wrong, most commonly when the distribution of $\varepsilon_i$ also depends on $X^{\mathrm{sr}}_i$, the result fails; Wooldridge (2010) pp. 606-607 points out that this requires a variable known to shift the mean of the latent variable while leaving its variance unchanged, and randomizing that variable does not make this any more credible. Wooldridge also says that the estimator identifies only relative coefficients and cannot deliver the APE. The first half of this claim holds, but the second goes too far: under the large support in SR.1, §5.2 of the working paper also consistently estimates the conditional distribution of $\varepsilon_i$ given $X_i$, from which the ASF and the APE follow, although this step requires a nonparametric regression and relies more heavily on the support condition. Who sets the weights: dividing by $f_{\mathrm{sr}}$ gives huge weights to observations where $X^{\mathrm{sr}}_i$ takes rare values, which makes the estimator very unstable in finite samples. For a problem like Helix's, whose target is the APE, the parametric route of Sections 3.4 to 3.5 actually delivers more.

### 3.7 Panel Binary Models: Pooled, RE, and CRE Probit

On the 8-month panel of block B5, a pooled logit that ignores the firm effects gives an engagement coefficient of 1.0006, whose interval $[0.942, 1.059]$ covers its own probability limit, 0.977, but lies far from $\beta = 0.8$; its APE is 0.1888, against a true value of 0.1245. The CRE probit gives an APE of 0.1235, with interval $[0.1140, 0.1329]$. The difference lies in how the correlation between $\alpha_i$ and engagement is treated.

::: {.assumption}
**PB.1 (strict exogeneity; Wooldridge 2010 equation (15.65), p. 610)** The distribution of $Y_{it}$ given $(X_{i1}, \ldots, X_{iT}, \alpha_i)$ depends only on $(X_{it}, \alpha_i)$.

**PB.2 (equations (15.66) and (15.77), pp. 611 and 619)** $\Pr(Y_{it} = 1 \mid X_{it}, \alpha_i) = F(X_{it}'\beta + \alpha_i)$, where $F(\cdot)$ is $\Phi(\cdot)$ or $\Lambda(\cdot)$.

**PB.3 (conditional independence; equation (15.67), p. 611)** Given $(X_{i1}, \ldots, X_{iT}, \alpha_i)$, $Y_{i1}, \ldots, Y_{iT}$ are mutually independent.

**PB.4 (Chamberlain's CRE; equation (15.73), p. 615)** Given the regressors of all periods, $\alpha_i$ is distributed as $\mathcal{N}(\xi_0 + \bar X_i'\xi, \sigma_a^2)$.
:::

PB.1 rules out lagged dependent variables (the dynamic panel case of Section 3.1 of Chapter 11) and regressors that respond to past shocks. PB.4 is the Mundlak projection of CRE.1 in Chapter 10 plus conditional normality; it allows $\alpha_i$ to be correlated with the regressors, at the price of having to write down the form of the correlation. The traditional RE probit is the special case $\xi = 0$ (equation (15.69)), which requires $\alpha_i$ to be independent of the regressors. When only the APE matters, PB.3 can be dropped.

::: {.theorem}
**Theorem 7 (pooled probit and CRE probit; Wooldridge 2010 §15.8.1 and §15.8.2, equations (15.62) to (15.64), (15.71), and (15.73) to (15.76), pp. 609-617)** (i) If $\Pr(Y_{it} = 1 \mid X_{it}) = F(X_{it}'\beta_{\mathrm{pool}})$ is correctly specified, the pooled MLE is the partial MLE of Section 3.8 of Chapter 6, consistent and asymptotically normal; with individual effects, dynamic completeness (equation (15.63)) fails, and the variance must be clustered by unit. (ii) Under PB.1, PB.2 (probit), and PB.4,

$$\Pr(Y_{it} = 1 \mid X_{i1}, \ldots, X_{iT}) = \Phi\big(\xi_{a0} + X_{it}'\beta_a + \bar X_i'\xi_a\big),$$

where the subscript $a$ means that a parameter is divided by $\sqrt{1 + \sigma_a^2}$, and a pooled probit on $(1, X_{it}, \bar X_i)$ estimates these parameters consistently. The APE over all firm-months, $\mathbb{E}\big[\beta_j\,\phi(X_{it}'\beta + \alpha_i)\big]$, equals $\mathbb{E}\big[\beta_{aj}\,\phi(\xi_{a0} + X_{it}'\beta_a + \bar X_i'\xi_a)\big]$, where the expectation is taken over the joint distribution of $(X_{it}, \bar X_i)$, and the sample average in which each firm-month uses its own $X_{it}$ and $\bar X_i$ estimates it consistently (by the argument of equations (15.75) and (15.76)). (iii) With PB.3 added, putting $\bar X_i$ into an RE probit estimated by MLE gives separate estimates of $\beta$ and $\sigma_a$.
:::

The proof sketch has three steps. Step one (Theorem 8 of Chapter 6): each period of the pooled likelihood is correctly specified, so the partial MLE is consistent; $\alpha_i$ makes the scores of different periods correlated, PL.2 fails, and the filling of the sandwich must be summed by unit. Step two (PB.1, PB.2, and PB.4): $\alpha_i = \xi_0 + \bar X_i'\xi + a_i$, with $a_i$ normal and independent of the regressors, and the argument of Section 3.4 integrates it out, giving the scaled probit (equation (15.74)). Step three (step two): $a_i$ is independent of $(X_{it}, \bar X_i)$, so, conditional on $(X_{it}, \bar X_i)$, taking the expectation over $a_i$ gives $\mathbb{E}[\beta_j\phi(X_{it}'\beta + \alpha_i) \mid X_{it}, \bar X_i] = \beta_{aj}\,\phi(\xi_{a0} + X_{it}'\beta_a + \bar X_i'\xi_a)$; iterated expectations then make the sample average a consistent estimate. Equations (15.75) and (15.76) hold $x$ fixed and average over $\bar X_i$ alone, which gives the APE with $\alpha_i$ drawn from its marginal distribution; averaging this further over $x = X_{it}$ yields another quantity, which in general differs from the APE above when $\alpha_i$ is correlated with $X_{it}$, as in B5.

In words: under PB.1, PB.2, and PB.4, the CRE probit estimates the coefficients divided by $\sqrt{1 + \sigma_a^2}$ and the APE over all firm-months, and it allows arbitrary serial correlation. When the specification is wrong: the pooled logit leaves $\alpha_i$ in the error, and $\alpha_i$ rises with average engagement, so its probability limit absorbs this correlation, 0.977 in the 8-month design and 1.045 in the 4-month design, and its APE, 0.1888, is one and a half times the truth. The DGP of B5 is in fact a logit with normal $a_i$, so the CRE probit has the wrong link, yet its APE, 0.1235, still lands very close to 0.1245, and the scaled coefficients of $X_{it}$ and $\bar X_i$, 0.410 and 0.506, are also close to the large-sample probability limits computed from 200,000 firms, 0.40 and 0.51; the random effects MLE of the CRE model gives an APE of 0.1238. The APE is insensitive to the shape of the link, the lesson of Section 3.3 once again. Who sets the weights: the APE averages over the joint distribution of $(X_{it}, \bar X_i)$ and gives every firm-month equal weight; clustering treats each firm as a single independent unit. When Section 3.8 of Chapter 6 said that Chapter 12 would meet cases in which PL.2 fails, this is what it had in mind: the cluster standard error of the pooled logit slope is 0.0299, against 0.0265 if the 8,000 firm-months are treated as independent observations.

### 3.8 FE Logit, Incidental Parameters, and Bias Correction

CRE has to specify the relation between $\alpha_i$ and the regressors. Can we, as FE does in Chapter 10, make no assumption at all about this relation? For logit we can, at the price of giving up the APE. The figure below shows the results of 1,000 simulated panels in the B5 design: treating each firm's $\alpha_i$ as a parameter and running a logit MLE with firm dummies overestimates by 0.311 (38.9%) on average with 4 months and by 0.130 with 8 months; the conditional MLE is almost unbiased at both lengths.

![Monte Carlo for block B5: 1,000 panels for each panel length, each panel with 1,000 firms. Left: 4 months (the 4-month DGP); right: 8 months. The orange solid line is the sampling density of the FE logit MLE with firm dummies, the green dashed line is its analytical bias correction (Fernández-Val and Weidner 2016), and the blue solid line is the conditional MLE; the vertical solid line is the true value 0.8, and the vertical dotted line is the probability limit of the pooled logit.](assets/fig/fig_12_fe_bias.svg)

::: {.assumption}
**FL.1 (Wooldridge 2010 equation (15.77), p. 619; p. 620)** PB.1 and PB.3 hold, $\Pr(Y_{it} = 1 \mid X_{it}, \alpha_i) = \Lambda(X_{it}'\beta + \alpha_i)$, and no assumption is made about the distribution of $\alpha_i$ given the regressors; every component of $X_{it}$ varies over time.
:::

::: {.theorem}
**Theorem 8 (FE logit as conditional MLE; Chamberlain 1980; Wooldridge 2010 §15.8.3, equations (15.81) to (15.83), pp. 620-622)** Under FL.1, the conditional distribution of $(Y_{i1}, \ldots, Y_{iT})$ given the regressors, $\alpha_i$, and $\sum_t Y_{it}$ does not contain $\alpha_i$. When $T = 2$, $\Pr(Y_{i2} = 1 \mid Y_{i1} + Y_{i2} = 1, X_{i1}, X_{i2}, \alpha_i) = \Lambda\big((X_{i2} - X_{i1})'\beta\big)$; for general $T$, the conditional log-likelihood of unit $i$ is

$$\ell_i(\beta) = \sum_{t=1}^T Y_{it}X_{it}'\beta - \ln \sum_{\mathbf{y} \in \mathcal{Y}_i} \exp\Big(\sum_{t=1}^T y_tX_{it}'\beta\Big),$$

where $\mathcal{Y}_i$ is the set of all 0/1 sequences with the same sum as $(Y_{i1}, \ldots, Y_{iT})$. Its maximizer $\hat\beta_{CMLE}$ is, with $T$ fixed and $N \to \infty$, consistent and $\sqrt N$-asymptotically normal; firms whose $\sum_t Y_{it}$ equals 0 or $T$ contribute no information; the APE is not identified.
:::

The proof sketch has three steps. Step one (FL.1 and PB.3): the probability of a sequence is $\exp\{\sum_t Y_{it}(X_{it}'\beta + \alpha_i)\}$ divided by a product that depends only on the regressors and $\alpha_i$; taking ratios among sequences with a fixed sum cancels both $\exp(\alpha_i\sum_t Y_{it})$ and the denominator (equation (15.82)), so $\sum_t Y_{it}$ is a sufficient statistic for $\alpha_i$. Step two (FL.1; Theorems 1 and 3 of Chapter 6): what remains is a conditional likelihood that contains only $\beta$ (equation (15.83)); when the sum is 0 or $T$, the set $\mathcal{Y}_i$ has a single element and the likelihood is identically one. Step three (FL.1: no assumption on the distribution of $\alpha_i$): the APE is $\mathbb{E}[\beta\Lambda'(X_{it}'\beta + \alpha_i)]$, which needs the distribution of $\alpha_i$, and that is what the conditional likelihood cancels (p. 622).

In words: under FL.1, $\hat\beta_{CMLE}$ estimates the coefficient on the log odds ratio scale, with no assumption about the relation between $\alpha_i$ and the regressors: the 8-month B5 sample gives 0.8006, with interval $[0.733, 0.868]$; the first 4 months give 0.7722, with interval $[0.669, 0.876]$; in the Monte Carlo, the coverage rates of the conditional MLE in the 4-month and 8-month designs are 0.951 and 0.952. It cannot deliver the APE, and this is the trade-off the chapter wants to stress. When the specification is wrong: if PB.3 fails, for example with serially correlated errors, the CMLE is generally no longer consistent (p. 622), whereas the CRE probit allows such correlation. Who sets the weights: only firms whose $Y_{it}$ changes enter the likelihood, and identification comes entirely from variation within firms; in the B5 sample these firms account for 0.802 of the firm-months over 8 months and 0.628 over the first 4 months.

Why does the logit MLE that treats $\alpha_i$ as parameters fail? Each $\hat\alpha_i$ is estimated from only $T$ observations, so its noise does not vanish however large $N$ is, and through the nonlinearity this noise turns into an $O(1/T)$ bias in $\hat\beta$: this is the incidental parameters problem posed by Neyman and Scott (1948) and surveyed by Lancaster (2000); the Nickell bias of Theorem 2 in Section 3.2 of Chapter 11 is its version in linear dynamic panels. The extreme case is $T = 2$, where this MLE converges to $2\beta$ (Wooldridge 2010 p. 621).

::: {.intuition}
Imagine a firm that uses Pro in 3 of its 4 months. Estimate its $\alpha_i$ as a parameter, with only 4 data points, and maximum likelihood will "read" it as a firm that likes Pro, pushing $\hat\alpha_i$ high; likewise, a firm that used Pro in only 1 month is read as one that shuns it. This noise in $\hat\alpha_i$ does not average out. Roughly speaking, the $\hat\alpha_i$, pushed toward the extremes, bring every firm's fitted probabilities closer to 0 or 1, and to explain the within-firm swings that go with engagement on so flat a stretch of the curve, the likelihood has no choice but to inflate $\beta$ as well. The more months, the more precise $\hat\alpha_i$ and the smaller the inflation, roughly in proportion to $1/T$. The conditional MLE simply does not estimate $\alpha_i$ at all: given how many months a firm used Pro, it asks only "which months," and the answer to that question does not depend on $\alpha_i$.
:::

::: {.theorem}
**Theorem 9 (incidental parameters and bias correction; Hahn and Newey 2004, Theorems 1 and 2 of cemmap working paper CWP17/03; Fernández-Val and Weidner 2016, Theorems 4.1, 4.3, and 4.4 of arXiv version v5)** Under FL.1, assume further that the regressors are independent and identically distributed over time (regularity conditions as in Hahn and Newey), and let $\hat\beta_{FE}$ be the MLE that treats $\alpha_i$ as parameters. With $T$ fixed it is inconsistent; as $N, T \to \infty$ with $N/T$ tending to a constant, $\sqrt{NT}(\hat\beta_{FE} - \beta)$ is asymptotically normal but with a nonzero mean, the bias being of order $1/T$; subtracting a consistent estimate of the bias gives $\hat\beta_{BC}$, whose asymptotic mean is zero. The analogous estimator of the APE also has a bias of order $1/T$ (Fernández-Val 2009, see Wooldridge 2010 p. 618), which can likewise be corrected analytically.
:::

The proof sketch has three steps. Step one (FL.1, $T \to \infty$): the error in $\hat\alpha_i(\beta)$ is $O_p(T^{-1/2})$, and once it is substituted into the score for $\beta$, the second-order term has a nonzero expectation, which gives a bias of $B/T$. Step two ($N/T$ tending to a constant): $\sqrt{NT}$ times $B/T$ converges to a nonzero constant, the asymptotic bias of Theorem 1. Step three (the sample counterpart of the bias is consistent): evaluate the sample counterpart at $\hat\beta$ to estimate $B$, subtract it, and the asymptotic mean becomes zero (Theorem 2).

In words: Theorem 9 rests on a large-$T$ approximation; what it estimates is $\beta$, and what it repairs is a bias that vanishes like $1/T$. When the specification is wrong, for example when PB.3 fails or the regressors are only predetermined (Section 3.7 of Chapter 11), the expression for the bias changes accordingly, and a correction written under strict exogeneity and conditional independence no longer hits the target: the correction of Fernández-Val and Weidner (2016) for predetermined regressors requires the additional choice of a trimming parameter for estimating spectral expectations (arXiv version v5, p. 22), and Wooldridge (2010) p. 618 also cautions that these results all maintain conditional independence (equation (15.67)) and regressors that are independent and identically distributed over time. The B5 Monte Carlo has the shape the theorem describes: without correction the overestimate is 0.311 with 4 months and 0.130 with 8 months, a ratio of 2.40, of the same order as the 2 predicted by a bias of order $1/T$ (the DGPs of the two lengths also differ slightly, with $\alpha_i$ built on the 4-month and 8-month means respectively); the coverage rates of the nominal 95% intervals are only 0.008 and 0.078. The correction of Fernández-Val and Weidner (2016) (`alpaca::biasCorr`; with individual effects only, it follows the same idea as the analytical correction of Hahn and Newey 2004) shrinks the bias to -0.015 and -0.0015, with coverage rates of 0.968 and 0.962.

The APE is a different matter and has to be read with care. Without correction, the all-firm APE is almost unbiased to begin with: lower by 0.70% with 4 months and lower by 0.22% with 8 months. This is a regularity, not a theorem: in the simulations of Hahn and Newey (2004) the bias of the probit APE is also small (as reported by Wooldridge 2010 p. 618), and Fernández-Val (2009) and Fernández-Val and Weidner (2016) §5.1 discuss the phenomenon further; in B5 it comes from two large errors that offset each other: for the firms that stay in the sample, the inflated coefficient makes their partial effects too large, contributing +0.0330 with 4 months; for the firms that are dropped, the estimated effect is zero, contributing -0.0339. The bias-corrected APE is instead too high by 2.70%: the analytical correction fixes the coefficient but does not improve the APE, and this +2.70% is not a robust regularity either. Who sets the weights: the all-firm APE gives every firm-month equal weight, and firms whose outcome never changes are counted with a zero effect; with 4 months the coverage rates of the APE intervals are 0.916 and 0.912, slightly below the nominal level, a short-panel phenomenon, and with 8 months they return to 0.947 and 0.943.

Dynamic models get only one sentence here: if $Y_{i,t-1}$ enters the equation, PB.1 no longer holds and $\alpha_i$ is correlated with the initial value $Y_{i0}$; the approach of Wooldridge (2005b) is to condition on $Y_{i0}$ and the exogenous variables of all periods and to specify a model like PB.4 for the distribution of $\alpha_i$ given them; telling state dependence apart from persistent preference heterogeneity is the subject of the anchoring paper of Section 5.3.

### 3.9 Ordered and Multinomial Response

In the single sample of block B6, the APE of engagement on the probability of "giving 5 stars" is 0.1118, and its interval $[0.1023, 0.1214]$ does not cover the true value 0.1008; the engagement coefficient, 0.562, exceeds the true 0.5 by 2.63 standard errors. In a Monte Carlo of 1,000 replications, the estimates average 0.1006 and the interval coverage rate is 0.953, and only 1.2% of the replications give an estimate of at least 0.1118: this sample is simply an unlucky draw, and there is nothing wrong with the estimator. B6 also hides a deeper problem: managed accounts give 5 stars at a rate of 0.279 and self-serve firms at 0.148, with average ratings of 3.343 against 3.194. Which group is more satisfied?

::: {.assumption}
**OR.1 (Wooldridge 2010 equations (16.21) and (16.22), p. 655; for heteroskedasticity, equation (16.29), p. 659)** $Y_i^* = X_i'\beta + \varepsilon_i$, where $X_i$ contains no constant and $\varepsilon_i \mid X_i \sim \mathcal{N}(0, \exp(2X_i'\delta))$; $\delta = 0$ gives the standard ordered probit (McKelvey and Zavoina 1975); the categories are $m = 0, 1, \ldots, M$, with $Y_i = m$ if and only if $\mu_m < Y_i^* \le \mu_{m+1}$, $\mu_0 = -\infty$ and $\mu_{M+1} = +\infty$, and the cut points $\mu_1 < \cdots < \mu_M$ are to be estimated. Helix's star rating is $m + 1$, with $M = 4$.
:::

::: {.theorem}
**Theorem 10 (ordered probit and the latent scale; Wooldridge 2010 §16.3.1 and §16.3.2, equations (16.23), (16.24), and (16.29), pp. 655-659; Bond and Lang 2019, Theorem 1 of NBER Working Paper 19950, p. 8)** Under OR.1, write $\sigma(X_i) = \exp(X_i'\delta)$: (i) $\Pr(Y_i = m \mid X_i) = \Phi\big((\mu_{m+1} - X_i'\beta)/\sigma(X_i)\big) - \Phi\big((\mu_m - X_i'\beta)/\sigma(X_i)\big)$, and the MLE is consistent and asymptotically normal; (ii) for a continuous regressor $x_j$ that enters only the mean, the APE on the highest category is $\mathbb{E}\big[(\beta_j/\sigma(X_i))\,\phi\big((\mu_M - X_i'\beta)/\sigma(X_i)\big)\big]$, with the same sign as $\beta_j$, while the signs for the middle categories are ambiguous; (iii) when the latent distributions of two groups have different scales, there is a monotone increasing transformation of the latent scale that reverses the ranking of the two group means, so ordered data by themselves cannot rank latent means.
:::

The proof sketch has three steps. Step one (OR.1): after standardization, $\varepsilon_i/\sigma(X_i)$ is standard normal, the category probabilities are differences of two distribution functions (equation (16.23)), the log-likelihood is equation (16.24), and Theorems 1 and 3 of Chapter 6 apply. Step two (OR.1, $x_j$ entering only the mean): differentiating the probability of the highest category leaves a single term; for a middle category it is the difference of two densities, whose sign depends on which cut point $x'\beta$ is closer to (p. 656). Step three (Bond and Lang): two normal distribution functions with different scales must cross, and on each side of the crossing point one group dominates; choose a monotone transformation that widens the gap on one side and compresses it on the other, and the ranking of the means flips.

In words: under OR.1, the ordered probit estimates the cut points and coefficients (whose scale is again fixed by normalization), and with them the category probabilities and APEs. In B6 the variance of the managed accounts is correctly specified, and the single sample estimates their latent mean to be higher by 0.210, with a standard deviation of 1.561 (true values 0.2 and 1.5); the LR statistic for the scale is 140.8. When the specification is wrong: if the variance equation is dropped, the engagement coefficient falls from 0.562 to 0.480, while its APE on the probability of a 5-star rating moves only from 0.1118 to 0.1148, so the coefficient moves more than the APE; the mean shift of the managed accounts also shrinks from 0.210 to 0.136. Who sets the weights, the question that matters most here: take $g(Y^*) = -\exp(-Y^*)$ as "true satisfaction" ($g$ is monotone increasing); when $Y^*$ is normal, $\mathbb{E}[g(Y^*)] = -\exp\{-\mathbb{E}[Y^*] + \mathrm{Var}(Y^*)/2\}$, and at the same engagement and bill the log ratio for the two groups is $-0.2 + (1.5^2 - 1)/2 = 0.425$; exponentiating, the $\mathbb{E}[g(Y^*)]$ of the managed group is 1.53 times that of the self-serve group, and since $g$ takes negative values, the managed group is in fact less satisfied; the sample estimate gives 0.508. Without conditioning there is no such ratio: the bill is lognormal, and $\mathbb{E}[g(Y^*)]$ diverges to minus infinity in both groups. Star ratings cannot deliver a cardinal comparison of latent satisfaction, and that is the warning of Bond and Lang (2019).

Multinomial responses have no order. The MNL of block B7 takes Free as the base category:

$$P_m(x) = \frac{\exp(x'\beta_{(m)})}{1 + \sum_{h=1}^M \exp(x'\beta_{(h)})}, \qquad m = 1, \ldots, M,$$

with $P_0(x) = 1/\{1 + \sum_h \exp(x'\beta_{(h)})\}$; in Helix, $m = 1, 2$ are Pro and Enterprise, with coefficients written $\beta_{\mathrm{Pro}}$ and $\beta_{\mathrm{Ent}}$.

::: {.theorem}
**Theorem 11 (MNL, IIA, and the Hausman-McFadden test; Wooldridge 2010 §16.2.1 and §16.2.2, equations (16.1) to (16.3) and (16.9) to (16.11), pp. 644-650; Hausman and McFadden 1984)** Suppose the MNL is correctly specified, $\{(Y_i, X_i)\}$ is a random sample, and $\mathbb{E}[X_iX_i']$ has full rank. (i) $\ln\{P_m(x)/P_h(x)\} = x'(\beta_{(m)} - \beta_{(h)})$ depends only on these two categories, and this is IIA; given that the choice falls in $m$ or $h$, choosing $m$ is a logit with parameter $\beta_{(m)} - \beta_{(h)}$. (ii) The log-likelihood of the MLE is globally concave, and the MLE is consistent and asymptotically normal. (iii) The partial effect of a continuous regressor $x_j$ on $P_m$ is $P_m(x)\{\beta_{(m)j} - \sum_h P_h(x)\beta_{(h)j}\}$ ($\beta_{(0)} = 0$), and its sign is not determined by $\beta_{(m)j}$. (iv) Dropping one category and re-estimating gives $\hat\beta_{\mathrm{sub}}$, which under IIA is consistent but less efficient than the full-sample $\hat\beta_{\mathrm{all}}$. Hausman and McFadden (1984) compare the two on this basis; written as the Hausman quadratic form of Theorem 7 of Chapter 10,

$$H = (\hat\beta_{\mathrm{sub}} - \hat\beta_{\mathrm{all}})'\big[\widehat{\mathrm{Avar}}(\hat\beta_{\mathrm{sub}}) - \widehat{\mathrm{Avar}}(\hat\beta_{\mathrm{all}})\big]^{-1}(\hat\beta_{\mathrm{sub}} - \hat\beta_{\mathrm{all}}) \xrightarrow{d} \chi^2_{\dim(\hat\beta_{\mathrm{sub}})}.$$
:::

The proof sketch has three steps. Step one (the MNL correctly specified, equation (16.1)): (i) and (iii) are its algebra and derivatives (equations (16.2) and (16.3)). Step two (the MNL correctly specified): the log-likelihood is $\sum_m \mathbf{1}\{Y_i = m\}\ln P_m$, and the Hessian is negative semidefinite (McFadden 1974, Wooldridge 2010 p. 645); Theorems 1 and 3 of Chapter 6 give (ii). Step three ((i) and Theorem 4 of Chapter 6): under IIA the conditional likelihood on the subset is still correctly specified, so $\hat\beta_{\mathrm{sub}}$ is consistent; the full-sample MLE uses all the information and is asymptotically more efficient, and the Hausman argument of Theorem 7 of Chapter 10 shows that the asymptotic variance of their difference equals the difference of their asymptotic variances, which is (iv).

In words: when the MNL is correctly specified (the conditions of Theorem 11), it estimates log-odds coefficients relative to the base category, together with the category probabilities and APEs derived from them. In B7 the engagement coefficient for Enterprise is 0.367 (true value 0.4), yet the APE on the Enterprise probability is -0.0443 (true value -0.0407): engagement pushes firms more strongly toward Pro, and the share of Enterprise is squeezed, which is what (iii) says. When the specification is wrong, that is, when IIA fails, the MNL estimates a pseudo-true value: in the nested version Pro and Enterprise are close to each other (the nest parameter is $\lambda = 0.4$; the smaller $\lambda \in (0, 1]$ is, the closer the alternatives within a nest are, and $\lambda = 1$ recovers the MNL), so among paying firms the choice is a logit with all coefficients divided by $\lambda$, and the engagement coefficient of Enterprise relative to Pro should be $(0.4 - 1.0)/0.4 = -1.5$; a logit on paying firms alone estimates -1.518, while the difference implied by the full-sample MNL is -1.329. Who sets the weights: the MLE weights by the observations in each category, the subsample without Free uses only paying firms, and IIA here means "the presence of Free does not affect the odds between Pro and Enterprise." On the nested sample the Hausman-McFadden test gives $H = 28.72$, with a $p$ value whose common logarithm is -5.59, a clear rejection of IIA; on the sample where IIA holds, $H = 3.15$ and $p = 0.369$. The finite-sample behavior of the classical statistic raises a further problem, taken up in Section 7. Chapter 27 will put this reduced-form MNL back into the structural random utility model, adding attributes of the alternatives, nested logit, and mixed logit.

## 4 Estimation: Computation, Estimators, and Software

### 4.1 Computation

The log-likelihoods of logit and probit are both globally concave, and Newton iterations starting from zero usually converge in a few steps, 6 for Helix; the MNL is globally concave as well (McFadden 1974). Concavity guarantees only that a maximum, if it exists, is unique, not that it exists: when some regressor or combination of regressors separates $Y_i = 1$ from $Y_i = 0$ perfectly (perfect separation), the likelihood keeps rising in that direction and has no interior maximum; the coefficients grow with every iteration, the standard errors become ever more absurd, and the software sometimes issues only a warning. Only a human check catches it: abnormally large coefficients and fitted probabilities that are all 0 or 1 most likely signal it; the remedy is to examine that regressor, not to switch algorithms. The full MLE of IV probit is no longer guaranteed to be concave, the estimated correlation coefficient can drift toward $\pm 1$, and the iterations may struggle to converge (Wooldridge 2010 p. 591), so run the two control-function steps first and use them as the starting point. The bootstrap can also run into trouble when one outcome is rare, with almost no $Y_i = 1$ in some resamples (p. 587). The conditional likelihood of FE logit sums over all sequences with the same sum, at most 70 terms when $T = 8$, so exact computation is easy; when $T$ is large the number of sequences grows rapidly, and a recursive algorithm is needed instead.

### 4.2 The Estimators at a Glance

| Estimator | What it estimates | Key assumptions | Helix's numbers |
|---|---|---|---|
| LPM (OLS) | Linear projection; equals the APE with normal regressors | LP.2 | 0.2153 (true APE 0.2156) |
| Logit and probit MLE | Normalized $\beta$; APE | BR.1, BR.2 | 1.097, 0.659; APE 0.2171, 0.2169 |
| Average cross-partial effect | Interaction effect on the probability scale | BR.1, BR.2 | -0.0161 (interaction coefficient 0.221) |
| Probit omitting $a_i$ | $\beta/\sqrt{1 + \sigma_a^2}$; APE | NH.1 | 0.638; APE 0.2135 |
| Control function | Scaled coefficients; APE via the ASF | CF.1, CF.2 | APE 0.1659 |
| IV probit MLE | $\beta$, $\beta_D$; APE | CF.1, CF.2 | APE 0.1660 |
| 2SLS on the LPM | Linear approximation | Instruments exogenous and relevant | 0.1559 |
| Pooled logit | Pseudo-true value contaminated by the correlation of $\alpha_i$ | None | 1.0006 (probability limit 0.977) |
| CRE probit | Scaled coefficients; APE | PB.1, PB.2, PB.4 | APE 0.1235 |
| FE logit (CMLE) | $\beta$ | FL.1 | 0.8006 |
| FE MLE with dummies | $\beta$ plus a bias of order $1/T$ | FL.1 | 0.931 |
| Bias-corrected FE | $\beta$ | FL.1, large $T$ | 0.7995 |
| Ordered probit | Cut points, coefficients; category APEs | OR.1 | APE on 5 stars 0.1118 |
| MNL | Log-odds coefficients; APE | IIA | APE of engagement on Enterprise -0.0443 |

The numbers in the panel rows all come from the 8-month sample. To choose among them, first answer "which number do you want": if it is the average effect on the probability scale, the APE is the target, and the LPM, probit, and logit often give similar answers; if you want to make no assumption about $\alpha_i$, FE logit gives $\beta$ but not the APE; if you want both, you have to write down the relation between $\alpha_i$ and the regressors, as CRE does.

### 4.3 Software

The computations in this chapter use the following R packages, and a few interface details affect the numbers. `glm(family = binomial("probit"))` and `binomial()` estimate probit and logit; their built-in variance is the inverse of the information matrix, `sandwich::sandwich` gives the sandwich, and for panels `sandwich::vcovCL` clusters by unit. `marginaleffects::avg_slopes` computes APEs, and its standard error contains only the Delta-method term, that is, the conditional variance given the sample regressors; when the target is the population APE, the sampling fluctuation of the regressors must be added, which on B2 changes the standard error of the APE of $x_1$ from 0.0052 to 0.0056. `survival::clogit` computes the CMLE of FE logit exactly, but refuses to give a cluster-robust variance under the exact method, so this chapter writes its own sandwich for the CMLE; on B5 it is 0.0353, against a model standard error of 0.0344. `alpaca::feglm` estimates FE logit with dummies, `alpaca::biasCorr` performs the bias correction of Fernández-Val and Weidner (2016), and `alpaca::getAPEs` gives the all-firm APE, whose variance is likewise the conditional variance given the regressors, to which this chapter adds the sampling variance of the partial effects; `marginaleffects` does not support `feglm` objects. `bife` is another implementation of FE logit and probit (Stammann, Heiss and McFadden 2016). `lme4::glmer` estimates the RE probit, with 8 quadrature points in this chapter. `ordinal::clm` estimates the ordered probit, `scale = ~ G` adds the variance equation, and without it the fit is the same as `MASS::polr`. `nnet::multinom` estimates the MNL; its default convergence tolerance shifts the value of the Hausman-McFadden statistic, so it must be tightened to `reltol = 1e-12`; `mlogit::hmftest` computes the classical Hausman-McFadden statistic and agrees with this chapter's hand-written version.

## 5 Anchoring Papers

Each of the three anchoring papers takes up one thread of the chapter: the coefficient of an interaction term is not the interaction effect; ordered categories cannot rank latent means; and state dependence has to be told apart from persistent heterogeneity in panel choice data. The first two are methods papers, each using examples to show why a common reading is wrong; the third is an application in marketing that uses a quasi-exogenous shock to test the conclusions of a standard discrete choice model. Each paper is laid out in five parts: paper, method, data, results, and limitations. The anchoring paper for multinomial response and IIA, McFadden (1974), appears in Section 5.1 of Chapter 27.

### 5.1 Ai and Norton (2003)

::: {.case}
Paper: "Interaction Terms in Logit and Probit Models," Economics Letters. It points out that in nonlinear models the coefficient of an interaction term is not the interaction effect, and it gives the cross-partial derivative together with its Delta-method standard error. The data and results below are taken from the companion paper by the same research group, Norton, Wang and Ai (2004, Stata Journal), which implements Ai and Norton's formulas and gives two examples.

Method: this is Theorem 3 (ii) of this chapter; when two dummy variables interact, the cross-partial derivative is replaced by a double difference, and when a continuous variable interacts with a dummy it is a difference of first derivatives; the standard errors come from the Delta method and are computed observation by observation.

Data: 12,365 adults aged 21 to 64 with complete variables in the 2000 Medical Expenditure Panel Survey (out of an original sample of 25,096); the outcome is whether the person visited an outpatient physician during the year, which about two thirds did.

Results: as the companion paper reports, Ai and Norton examined the 72 articles published from 1980 to 2000 in 13 economics journals on JSTOR that used interaction terms in nonlinear models, and none of them interpreted the coefficient of the interaction term correctly. The first example is a logit with an interaction of age and years of education: the interaction coefficient is -0.00137 ($z = -2.66$) and the average interaction effect is -0.00033, but person by person the interaction effect ranges from -0.00058 to 0.00016, some negative and some positive. The second example is a probit with an interaction of male and uninsured: the interaction coefficient is 0.054 ($z = 0.89$, insignificant), yet the average interaction effect is -0.0093, and the individual $z$ statistics range from -5.56 to 5.53, significant for many people.

Limitations: the cross-partial effect is itself a function of $x$, and reporting an average likewise hides changes in sign, so it is best plotted, as in the figure of Section 3.3; if the research question is about latent utility in the first place, for example whether two attributes are complements in a structural model, the interaction coefficient is then the quantity wanted, and the cross-partial effect on the probability scale is not the only correct convention. For further discussion of tests on interaction terms and of nonlinear DiD, see Greene (2010) and Puhani (2012).
:::

This short paper changed interpretation, not estimation: the coefficients are estimated as before, but the analysis cannot stop there. The 0.221 and -0.0161 of B2 are a simulated version of the same mechanism.

### 5.2 Bond and Lang (2019)

::: {.case}
Paper: "The Sad Truth about Happiness Scales," Journal of Political Economy. The numbering and numbers below follow its NBER Working Paper No. 19950 (March 2014 version). Its argument is that when a survey offers only a few ordered categories, two groups of people cannot be ranked by average happiness without strong auxiliary assumptions.

Method: this is Theorem 10 (iii) of this chapter. The working paper assumes that utility is unbounded, that the distribution functions are continuous with positive densities, and that the distributions belong to a location-scale family (Assumptions 1 to 3, p. 7). It first proves that the distribution functions of two groups must cross unless their scales are equal (Proposition 1, p. 8), and then that there always exists a monotone transformation of utility that fits the categorical data equally well yet reverses the ranking of the means (Theorem 1, p. 8); in large samples the estimated scales of two groups are almost never exactly equal, so such a reversal is almost always possible (Conclusion 1). Empirically, the authors fix two cut points at 0 and 1, estimate the mean and variance of each group's normal distribution with an ordered probit, and then apply exponential transformations to the latent variable to see whether the ranking flips; the $g(Y^*) = -\exp(-Y^*)$ of this chapter's B6 is one of the left-skewed ones.

Data: the authors reexamine three conclusions of the happiness literature, all based on questions with three categories, "very happy, pretty happy, not too happy." The first is the long-term follow-up of the Moving to Opportunity (MTO) housing voucher experiment: families in public housing in high-poverty neighborhoods were randomly divided into three groups, and the authors compare the control group with the experimental group, whose vouchers could be used only in low-poverty neighborhoods, using the published response distributions of the two groups. The second is the General Social Survey (GSS) file compiled by Stevenson and Wolfers, using the five surveys from 1972 to 1976 and the five from 1998 to 2006, to compare changes in the happiness gap between men and women. The third uses the same GSS and the annual series of U.S. real disposable income per capita to test the Easterlin paradox. The authors also use the World Values Survey to compare country rankings.

Results: in MTO, with normal distributions and the cut points fixed at 0 and 1, the control group does have a lower mean (0.44 against 0.60) but a higher variance (0.79 against 0.63), and the distribution functions of the two groups cross at the 83rd percentile, that is, at 1.20; raising the utility values above 1.20 is enough to reverse the ranking of the means. With lognormal distributions instead, the experimental group is still happier (2.22 against 2.14), but raising the upper cut point to 1.33 makes the two means equal. In the GSS, the ordered probit shows that women's happiness relative to men's is lower by 0.09 standard deviations in the later period than in the earlier one, a decline that shrinks to 0.07 once the variance is allowed to vary by sex and period; with a further left-skewed lognormal transformation, once the transformation parameter reaches 3.9, women instead become happier relative to men. The Easterlin paradox is no different: under normality, in a regression of average annual happiness on real disposable income per capita, each additional USD 10,000 of income goes with average happiness lower by 0.02 units ($p$ value 0.11), while the standard deviation of happiness is significantly lower by 0.04; when the parameter of the left-skewed transformation exceeds 2, the positive correlation between income and average happiness becomes significant. The rankings of countries are just as unstable: Ghana is the happiest country under a right-skewed distribution and ranks third from the bottom under a left-skewed one.

Limitations: this is an impossibility result and offers no alternative estimator by itself; the ways out are either to treat the survey categories themselves as the policy-relevant scale or to impose restrictions on the distribution of happiness that everyone can accept (Bond and Lang 2019, working paper version, §5). The transformation needed for a reversal differs from example to example: making the MTO effect negative requires a strongly right-skewed distribution, whereas dissolving the Easterlin paradox requires a left-skewed one, and the authors also note that there may be no distributional assumption that reverses both at once. The empirical part fits two-parameter distributions without covariates, whereas this chapter's B6 compares the groups given the covariates, as in Section 3.9.
:::

This paper matters for the chapter because it brings Theorem 10 (iii) to bear on real research conclusions: ordered categories tell us only between which two cut points the latent variable falls, the comparison of group means depends on the shape of the latent distributions, and the data say almost nothing about that shape. In B6 the managed accounts give more 5-star ratings and have a higher average rating, yet under a monotone transformation they are less satisfied, a simulated version of the same mechanism.

### 5.3 Levine and Seiler (2023)

::: {.case}
Paper: "Identifying State Dependence in Brand Choice: Evidence from Hurricanes," Marketing Science. The numbers below follow its CESifo Working Paper No. 9889 (version of September 1, 2022).

Method: a household that buys the same brand again and again may do so because its past purchases themselves changed its choices (state dependence), or simply because it has always preferred the brand (persistent preference heterogeneity). The standard approach estimates a random-coefficients logit with an indicator for "bought this brand last time," identified from variation in prices over time, and it has to model the initial conditions, which is the problem raised in the last sentence of Section 3.8. The authors take a different route: panic buying before a hurricane causes stockouts that force households to switch brands; without state dependence, the persistence of choices should return to its original path as soon as the hurricane has passed.

Data: bottled water purchases in the Nielsen household scanner panel, 14 hurricanes from 2006 to 2017, 15,838 affected household-hurricane pairs, and 31,676 control pairs; the main analysis keeps 2,201 treated households and 3,866 control households.

Results: in the week of the hurricane, the persistence of brand choice drops by 0.071; the effect in weeks 2 to 5 after the hurricane is 0.010, with 95% interval $(-0.010, 0.033)$, so "no state dependence" cannot be rejected. When a random-coefficients logit is estimated on the same category by the standard approach, the mean of the state dependence coefficient is 1.233, with a 95% posterior interval from 0.943 to 1.540, significantly positive; with these estimates, a simulation of the same stockout shock predicts an effect of -0.046 in weeks 2 to 5, outside the interval above.

Limitations: the test targets brand switching caused by stockouts, and extending it to switching caused by discounts requires assuming that the two kinds of variation act the same way; there is only one product category; and the comparison of treated and control households uses methods in the synthetic control family, which this book develops when Part III takes up causal designs.
:::

This paper matters for the chapter because it sets the identification problem of Sections 3.7 to 3.8 in a real marketing question: with the same panel choice data, a model that relies on functional form and on modeling the initial conditions says there is state dependence, while a quasi-exogenous shock says there is not. The credibility of the dynamic conclusions of a discrete choice model depends on how heterogeneity is specified.

## 6 A Full Walkthrough on Helix Data

The code below runs under R 4.6.1. It lists only the statements that generate the data, run the estimations, and summarize the Monte Carlo, with several statements merged onto one line; every number in the text comes from the output of an actual run, and the result lines in the code that start with `#` copy what the run printed. The code keeps the inherited variable names: adoption is written `D`, engagement `x`, the true values `a0` and `a1`, `ame` is the APE, and `mem` is the PEA; the endogenous regressor of B3 is also called `D` in the code and has nothing to do with the adoption of B1. `ci95` returns the estimate plus or minus the normal quantile times the standard error; `binary_mle` estimates a probit or logit with `glm.fit` and returns the coefficients, the inverse of the information matrix, and the observation-by-observation scores; `rivers_vuong` calls `rv_fit`, listed in Section 6.3, on the full sample, redoes the two steps on 499 pairs-bootstrap samples, and returns, for $\hat v_i$, the usual probit $t$ statistic; `ivprobit_mle` starts from the two-step estimates of `rv_fit` and maximizes the likelihood of equation (15.50); `fe_logit_ape` takes the all-firm APE of `alpaca::getAPEs` and adds the sampling variance of the partial effects to its conditional variance; `hm_test` fits the full model and the model with one category dropped and computes the $H$ of Section 3.9, and `hm_test_stacked` is the stacked-score version of Section 7; `mnl_ape` computes the APEs of the MNL.

### 6.1 The DGP

Block B1 is inherited code: seed 512 generates the 4,000 firms, and seed 9001 generates a reference population of five million firms, used to compute the true values. The inherited generator also draws a usage variable that this chapter does not use; it is omitted here, and since it is drawn after adoption, it does not affect `x` and `D`.

```r
set.seed(512)
a0 <- -0.40; a1 <- 1.10          # logit adoption:  P(D=1|x) = plogis(a0 + a1 x)
n_main <- 4000                   # main analysis sample
gen_helix <- function(n) {
  x  <- rnorm(n)                                   # standardized engagement
  ## (A) adoption -- logit
  pD <- plogis(a0 + a1 * x)
  D  <- rbinom(n, 1, pD)
  ## (the legacy generator then draws a usage count that this chapter does not use)
  data.frame(x = x, D = D, pD = pD)
}
dat <- gen_helix(n_main)
set.seed(9001)
big <- gen_helix(5e6)
ame_logit <- mean(a1 * dlogis(a0 + a1 * big$x))
mem_logit <- a1 * dlogis(a0 + a1 * 0)
base_rate <- mean(big$D)
cat(sprintf("adoption rate     = %.4f (oracle base rate %.4f)\n", mean(dat$D), base_rate))
cat(sprintf("oracle AME (logit)= %.4f\n", ame_logit))
cat(sprintf("oracle MEM (logit)= %.4f\n", mem_logit))
# adoption rate     = 0.4275 (oracle base rate 0.4206)
# oracle AME (logit)= 0.2156
# oracle MEM (logit)= 0.2643
```

Each of the other six blocks has its own seed (1202 to 1207 for B2 to B7 in turn, and 1217 for the nested version of B7), with the parameters below. The generator of B3 and the panel generator of B5 are listed unchanged, and the other four generate their blocks as described in Section 2.2; in B3, `r` is the correlation coefficient of $(\varepsilon_i, v_i)$, which, when $\sigma_v = 1$, equals the $\rho$ of the text; in B6, `G` indicates managed accounts, `delta_G` is their latent mean shift, unrelated to the $\delta$ of the variance function in the text, and `sigma` holds the standard deviations of the latent variable in the two groups.

```r
P2 <- list(n = 4000, rho_x = 0.3, beta = c(1.0, 1.0, 0.7, 0.2))                       # B2: x1, x2, x1 * x2
P3 <- list(n = 4000, spend_median_kusd = 5, sd_logspend = 0.8, p_reminder = 0.5, pi = c(-0.4, 0.3, 0.8),
           sigma_v = 1, beta = c(-0.6, 0.3), beta_D = 0.5, r = 0.5)                   # B3: Corr(eps, v) = r
P4 <- list(n = 4000, beta = c(-0.2, 0.9), sigma_a = 1)                               # B4: omitted a_i
P5 <- list(N = 1000, T = 8, T_short = 4, sigma_c = 0.7, sigma_a = 1, xi0 = -0.3, xi = 1.0, beta = 0.8,
           R_mc = 1000, N_large = 2e5)                                                # B5: monthly panel
P6 <- list(n = 3000, beta = c(0.5, -0.3), delta_G = 0.2, sigma = c(1, 1.5), cuts = c(-1.9, -1.1, -0.2, 0.8),
           bill_median_kusd = 1.2, sd_logbill = 0.5, p_managed = 0.3, R_mc = 1000)    # B6: 1 to 5 rating
P7 <- list(n = 5000, levels = c("Free", "Pro", "Enterprise"),
           coef = rbind(Pro = c(0.7, 1.0, 0.3), Enterprise = c(0.0, 0.4, 1.1)),
           lambda_nested = 0.4, drop = "Free", R_mc = 4000)                           # B7: plan choice
gen_B3 <- function(p, n = p$n) {
  X <- p$sd_logspend * rnorm(n)
  spend_kusd <- p$spend_median_kusd * exp(X)
  Z <- rbinom(n, 1, p$p_reminder)
  v <- p$sigma_v * rnorm(n)
  e <- rnorm(n)
  D <- p$pi[1] + p$pi[2] * X + p$pi[3] * Z + v
  eps <- p$r * v / p$sigma_v + sqrt(1 - p$r^2) * e
  Y <- as.integer(p$beta[1] + p$beta[2] * X + p$beta_D * D + eps > 0)
  data.frame(spend_kusd = spend_kusd, X = X, Z = Z, D = D, Y = Y, v_true = v, eps_true = eps)
}
gen_B5 <- function(p, N = p$N, T = p$T) {
  c_i <- p$sigma_c * rnorm(N)
  X <- rep(c_i, each = T) + rnorm(N * T)
  id <- rep(seq_len(N), each = T)
  Xbar <- ave(X, id)
  a <- p$sigma_a * rnorm(N)
  alpha <- p$xi0 + p$xi * Xbar + rep(a, each = T)
  Y <- rbinom(N * T, 1, plogis(alpha + p$beta * X))
  data.frame(id = id, t = rep(seq_len(T), times = N), X = X, Y = Y, alpha_true = alpha)
}
# gen_B2, gen_B4, gen_B6 and gen_B7 draw blocks B2, B4, B6 and B7 from P2, P4, P6 and P7 as described in the text
```

The true values of the blocks are computed by numerical integration: normal grids or closed forms for B2, B3, B4, B6, and B7, and two-dimensional Gauss-Hermite quadrature for B5, with the results stored in `truth_B2` to `truth_B7`; each block also draws a large sample as a check. B5 has two sets of true values: the 8-month design (probability limit of the pooled logit 0.97690, APE 0.124527, the same for the subpanel of the first 4 months) and the 4-month DGP (1.04498 and 0.121431); the 4-month Monte Carlo uses the latter.

### 6.2 Block B1: LPM, Logit, and Probit

The LPM and logit code is inherited; only references such as `oracle$` have been replaced by the `a0`, `a1`, and `ame_logit` defined above, and file reading and writing, along with printing that this chapter does not use, are omitted.

```r
lpm <- lm(D ~ x, data = dat)
b   <- coef(lpm)
phat <- fitted(lpm)
cat(sprintf("LPM intercept     = %.4f\n", b[1]))
cat(sprintf("LPM slope (x)     = %.4f (SE %.4f)\n",
            b[2], summary(lpm)$coefficients[2, 2]))
cat(sprintf("min fitted phat   = %.4f\n", min(phat)))
cat(sprintf("max fitted phat   = %.4f\n", max(phat)))
cat(sprintf("share phat<0 or >1= %.4f\n", mean(phat < 0 | phat > 1)))
cat(sprintf("  share phat > 1   = %.4f\n", mean(phat > 1)))
cat(sprintf("  share phat < 0   = %.4f\n", mean(phat < 0)))
# LPM intercept     = 0.4271
# LPM slope (x)     = 0.2153 (SE 0.0071)
# min fitted phat   = -0.4778
# max fitted phat   = 1.2133
# share phat<0 or >1= 0.0293
#   share phat > 1   = 0.0040
#   share phat < 0   = 0.0253
```

```r
X <- cbind(1, dat$x); y <- dat$D
newton_logit <- function(X, y, tol = 1e-10, maxit = 100) {
  b <- rep(0, ncol(X))
  for (it in 1:maxit) {
    Xb <- as.vector(X %*% b); F <- plogis(Xb)
    score <- crossprod(X, y - F)                 # k x 1
    W <- F * (1 - F)
    H <- -crossprod(X, X * W)                     # -sum F(1-F) X X'
    step <- solve(H, score)
    b <- b - step
    if (max(abs(step)) < tol) break
  }
  list(coef = b, iter = it, F = plogis(as.vector(X %*% b)), W = W)
}
fit <- newton_logit(X, y)
g <- glm(D ~ x, data = dat, family = binomial())
Fh   <- fit$F; Wh <- Fh * (1 - Fh)
I_hat <- crossprod(X, X * Wh)                      # observed = expected info (logit)
V_fisher <- solve(I_hat)
se_fisher <- sqrt(diag(V_fisher))
B <- crossprod(X, X * (y - Fh)^2)
V_sand <- V_fisher %*% B %*% V_fisher
se_sand <- sqrt(diag(V_sand))
cat(sprintf("iterations to converge = %d\n", fit$iter))
cat(sprintf("beta0 = %.4f   beta1 = %.4f   (truth: %.2f, %.2f)\n", fit$coef[1], fit$coef[2], a0, a1))
cat(sprintf("glm    beta0 = %.4f   beta1 = %.4f  (max abs diff %.2e)\n",
            coef(g)[1], coef(g)[2], max(abs(fit$coef - coef(g)))))
cat(sprintf("beta1 Fisher SE   = %.4f\n", se_fisher[2]))
cat(sprintf("beta1 sandwich SE = %.4f\n", se_sand[2]))
ame_hat <- mean(dlogis(as.vector(X %*% fit$coef)) * fit$coef[2])
cat(sprintf("AME (logit MLE)   = %.4f   (oracle %.4f)\n", ame_hat, ame_logit))
# iterations to converge = 6
# beta0 = -0.3672   beta1 = 1.0969   (truth: -0.40, 1.10)
# glm    beta0 = -0.3672   beta1 = 1.0969  (max abs diff 1.23e-13)
# beta1 Fisher SE   = 0.0442
# beta1 sandwich SE = 0.0439
# AME (logit MLE)   = 0.2171   (oracle 0.2156)
```

The standard errors of APEs are given by the function below: the Delta-method term conditions on the sample regressors, the other term is the sampling fluctuation of the individual partial effects, and their squares are added, as in Theorem 3 (i); for panels the latter term is clustered by unit.

```r
avg_effect_inference <- function(g, theta, V, cluster = NULL) {    # g(theta): the n per-observation effects
  gi <- g(theta); est <- mean(gi)
  G <- numDeriv::grad(function(t) mean(g(t)), theta)
  se_delta <- sqrt(drop(t(G) %*% V %*% G))                          # delta method, given the X_i
  se_sample <- if (is.null(cluster)) sd(gi) / sqrt(length(gi)) else sqrt(sum(rowsum(gi - est, cluster)^2)) / length(gi)
  se <- sqrt(se_delta^2 + se_sample^2)                              # standard error for the population APE
  list(est = est, se = se, ci = ci95(est, se))
}
ape_l <- avg_effect_inference(function(t) dlogis(drop(X %*% t)) * t[2], drop(fit$coef), V_fisher)
gp <- glm(D ~ x, data = dat, family = binomial("probit"), control = glm.control(epsilon = 1e-12, maxit = 100))
bp <- coef(gp)
ape_p <- avg_effect_inference(function(t) dnorm(drop(X %*% t)) * t[2], bp, vcov(gp))
cat(sprintf("LPM slope HC1 SE %.4f; logit APE %.4f (SE %.4f) [%.4f, %.4f], PEA %.4f\n",
            sqrt(sandwich::vcovHC(lpm, type = "HC1")[2, 2]), ape_l$est, ape_l$se, ape_l$ci[1], ape_l$ci[2],
            fit$coef[2] * dlogis(fit$coef[1] + fit$coef[2] * mean(dat$x))))
cat(sprintf("probit beta1 %.4f (SE %.4f), APE %.4f (SE %.4f), PEA %.4f; logit / probit coefficient ratio %.3f\n",
            bp[2], sqrt(vcov(gp)[2, 2]), ape_p$est, ape_p$se, bp[2] * dnorm(bp[1] + bp[2] * mean(dat$x)), fit$coef[2] / bp[2]))
# LPM slope HC1 SE 0.0059; logit APE 0.2171 (SE 0.0060) [0.2054, 0.2289], PEA 0.2653
# probit beta1 0.6593 (SE 0.0250), APE 0.2169 (SE 0.0059), PEA 0.2567; logit / probit coefficient ratio 1.664
```

The largest fitted value, 1.2133, is the 121% of Section 1. The conventional LPM standard error, 0.0071, is too large, and the heteroskedasticity-robust HC1 is 0.0059: LPM errors are heteroskedastic by nature (the variance formula of Theorem 1 (i); when LP.1 fails, the projection errors are heteroskedastic too), and here the conventional formula errs on the conservative side. The logit APE interval $[0.2054, 0.2289]$ covers the true value 0.2156.

The Monte Carlo for consistency is also inherited, with the same changes: the DGP is fixed, and each sample size is replicated 400 times.

```r
set.seed(778)
gen <- function(n){ x<-rnorm(n); D<-rbinom(n,1,plogis(a0+a1*x)); data.frame(x,D) }
ns  <- c(100, 400, 1600, 6400, 25600, 102400)
R   <- 400
cons <- sapply(ns, function(n){
  est <- replicate(R, coef(glm(D~x, gen(n), family=binomial()))[2])
  c(mean = mean(est), sd = sd(est), rmse = sqrt(mean((est-a1)^2)))
})
colnames(cons) <- ns
print(round(cons, 4))
#         100    400   1600   6400  25600 102400
# mean 1.1355 1.1091 1.0972 1.1007 1.1011 1.0997
# sd   0.2791 0.1473 0.0708 0.0353 0.0170 0.0085
# rmse 0.2810 0.1474 0.0708 0.0352 0.0171 0.0085
```

At $n = 100$ the average slope is 1.1355, clearly above 1.10, which is the finite-sample bias of MLE; each time the sample size is multiplied by 4, the standard deviation roughly halves, by a factor of 1.89 from 100 to 400 and of 2.01 from 25,600 to 102,400, the signature of $\sqrt n$ convergence. The average at a fixed sample size says nothing about consistency by itself; consistency means that the distribution shrinks toward 1.10, and the standard deviations in the table shrink in just this way.

### 6.3 Blocks B2 to B4: Interaction Terms, the Control Function, and Neglected Heterogeneity

```r
cross_partial_i <- function(beta, X, j, k, link = c("probit", "logit")) {   # d2 P / dx_j dx_k at each row of X
  L <- link_fns(link)                                    # F, F' and F'' of the link
  jk <- interaction_col(X, j, k); b_jk <- if (is.na(jk)) 0 else unname(beta[jk])   # coefficient of x_j * x_k
  u <- drop(X %*% beta)
  b_jk * L$Fp(u) + (beta[j] + b_jk * X[, k]) * (beta[k] + b_jk * X[, j]) * L$Fpp(u)
}
set.seed(1202); d2 <- gen_B2(P2)
X2 <- cbind("(Intercept)" = 1, x1 = d2$x1, x2 = d2$x2, "x1:x2" = d2$x1 * d2$x2)
m2 <- binary_mle(d2$y, X2, "probit"); b2 <- m2$coef; V2 <- m2$vcov_model
g2 <- glm(y ~ x1 * x2, data = d2, family = binomial("probit"), control = glm.control(epsilon = 1e-12, maxit = 100))
g2r <- glm(y ~ x1 + x2, data = d2, family = binomial("probit"), control = glm.control(epsilon = 1e-12, maxit = 100))
lr_int <- 2 * (as.numeric(logLik(g2)) - as.numeric(logLik(g2r)))
cross <- avg_effect_inference(function(t) cross_partial_i(t, X2, 2, 3, "probit"), b2, V2)
ci_i <- cross_partial_i(b2, X2, 2, 3, "probit")
naive <- mean(dnorm(drop(X2 %*% b2))) * b2[4]                            # beta_jj' phi(u) only: not the cross-partial
cat(sprintf("B2: beta_jj' %.4f (SE %.4f, z %.2f; LR %.2f); average cross-partial %.5f (SE %.5f) [%.4f, %.4f]; share positive %.3f; beta_jj' mean phi %.4f\n",
            b2[4], sqrt(V2[4, 4]), b2[4] / sqrt(V2[4, 4]), lr_int, cross$est, cross$se, cross$ci[1], cross$ci[2], mean(ci_i > 0), naive))
# B2: beta_jj' 0.2213 (SE 0.0422, z 5.24; LR 27.56); average cross-partial -0.01614 (SE 0.00563) [-0.0272, -0.0051]; share positive 0.337; beta_jj' mean phi 0.0456
```

```r
rv_fit <- function(y, X_exog, x_endog, Z, idx = seq_along(y)) {           # the two steps on the rows idx
  Xe <- X_exog[idx, , drop = FALSE]; d <- x_endog[idx]; W <- cbind(Xe, Z[idx, , drop = FALSE])
  pi_hat <- qr.coef(qr(W), d); vhat <- d - drop(W %*% pi_hat)             # step 1: OLS first stage, residual
  R <- cbind(Xe, d, vhat)
  b <- glm.fit(R, y[idx], family = binomial("probit"),
               control = glm.control(epsilon = 1e-12, maxit = 100))$coefficients   # step 2: probit with v-hat
  k <- ncol(Xe); a <- drop(Xe %*% b[1:k]) + b[k + 1] * d
  ape <- b[k + 1] * mean(dnorm(outer(a, b[k + 2] * vhat, "+")))           # ASF: average over i and over v-hat
  list(coef = unname(b), pi = unname(pi_hat), ape = unname(ape), vhat = vhat, R = R)
}
set.seed(1203); d3 <- gen_B3(P3)
Xe <- cbind("(Intercept)" = 1, X = d3$X); Xn <- cbind(Xe, D = d3$D)
mn <- binary_mle(d3$Y, Xn, "probit")                                      # naive probit of Y on (1, X, D)
naive3 <- avg_effect_inference(function(t) dnorm(drop(Xn %*% t)) * t[3], mn$coef, mn$vcov_model)
fs <- lm(D ~ X + Z, data = d3)
fs_F <- unname(coef(fs)["Z"] / sqrt(sandwich::vcovHC(fs, type = "HC0")["Z", "Z"]))^2
cf <- rivers_vuong(d3$Y, Xe, d3$D, cbind(Z = d3$Z), B = 499, seed = 1233)  # rv_fit on 499 pairs-bootstrap samples
iv <- ivprobit_mle(d3$Y, Xe, d3$D, cbind(Z = d3$Z))                       # maximises the likelihood (15.50)
iv_ape <- avg_effect_inference(function(t) t[3] * dnorm(drop(Xe %*% t[1:2]) + t[3] * d3$D), iv$coef, iv$vcov)
Xh <- cbind(Xe, D = fitted(fs)); b_iv <- solve(crossprod(Xh, Xn), crossprod(Xh, d3$Y))   # 2SLS of the LPM
e_iv <- drop(d3$Y - Xn %*% b_iv); A_iv <- solve(crossprod(Xh, Xn))
se_iv <- sqrt((A_iv %*% crossprod(Xh * e_iv) %*% t(A_iv))[3, 3])
cat(sprintf("B3: naive probit APE %.4f [%.4f, %.4f]; CF APE %.4f (bootstrap SE %.4f) [%.4f, %.4f]; t on v-hat %.2f; first-stage F %.1f\n",
            naive3$est, naive3$ci[1], naive3$ci[2], cf$ape, cf$ape_se, cf$ape_ci[1], cf$ape_ci[2], cf$rho_t, fs_F))
cat(sprintf("    CF scaled coefficients: X %.4f, D %.4f, v-hat %.4f\n", cf$beta[2], cf$beta[3], cf$beta[4]))
cat(sprintf("    IV-probit MLE: beta_D %.4f [%.4f, %.4f], Corr %.3f, APE %.4f [%.4f, %.4f]; 2SLS %.4f [%.4f, %.4f]\n",
            iv$beta[3], ci95(iv$beta[3], iv$se_beta[3])[1], ci95(iv$beta[3], iv$se_beta[3])[2], iv$r,
            iv_ape$est, iv_ape$ci[1], iv_ape$ci[2], b_iv[3], ci95(b_iv[3], se_iv)[1], ci95(b_iv[3], se_iv)[2]))
# B3: naive probit APE 0.2410 [0.2325, 0.2496]; CF APE 0.1659 (bootstrap SE 0.0148) [0.1370, 0.1949]; t on v-hat 7.11; first-stage F 642.2
#     CF scaled coefficients: X 0.3548, D 0.6474, v-hat 0.4799
#     IV-probit MLE: beta_D 0.5831 [0.4384, 0.7278], Corr 0.434, APE 0.1660 [0.1357, 0.1962]; 2SLS 0.1559 [0.1259, 0.1860]
```

The scaled true values are $(\beta, \beta_D, \rho)/\sigma_e$ with $\sigma_e = 0.866$: 0.346 for $X_i$, and 0.577 for both $D_i$ and $\hat v_i$, and the intervals of all three estimates cover them. The IV probit MLE estimates the error correlation at 0.434 (true value 0.5). In B3, $\Gamma_0$ is not zero, so the interval for the control-function APE comes from a bootstrap of the two steps together, not from the standard errors that the second-step probit reports.

```r
set.seed(1204); d4 <- gen_B4(P4)
Xs <- cbind(1, X = d4$X); ms <- binary_mle(d4$Y, Xs, "probit")            # short probit: a_i omitted
short <- avg_effect_inference(function(t) dnorm(drop(Xs %*% t)) * t[2], ms$coef, ms$vcov_model)
bci <- ci95(ms$coef[2], sqrt(ms$vcov_model[2, 2]))
Xl <- cbind(Xs, a = d4$a_true); ml <- binary_mle(d4$Y, Xl, "probit")     # infeasible long probit: a_i is simulation truth
cat(sprintf("B4: short probit beta %.4f [%.4f, %.4f]; APE %.4f [%.4f, %.4f]; long probit beta %.4f, APE %.4f; LPM slope %.4f\n",
            ms$coef[2], bci[1], bci[2], short$est, short$ci[1], short$ci[2], ml$coef[2],
            mean(dnorm(drop(Xl %*% ml$coef))) * ml$coef[2], coef(lm(Y ~ X, data = d4))[2]))
# B4: short probit beta 0.6376 [0.5901, 0.6852]; APE 0.2135 [0.2021, 0.2250]; long probit beta 0.9077, APE 0.2144; LPM slope 0.2131
```

The short-probit coefficient is 0.702 times the long-probit one, against a theoretical value of $1/\sqrt 2 = 0.707$. The long-probit APE, 0.2144, is averaged over the sample $a_i$ and estimates the same quantity as the short probit's 0.2135.

### 6.4 Block B5: The Panel

```r
library(survival)                                                            # clogit() needs strata() and coxph() attached
fit_panel <- function(d) {
  d$fid <- factor(d$id)
  fe <- alpaca::feglm(Y ~ X | fid, data = d, family = binomial("logit"))    # alpha_i as parameters
  bc <- alpaca::biasCorr(fe)                                                # analytical bias correction
  ap <- fe_logit_ape(fe, bc, d$id, d$Y)                                     # all-firm APEs and their SEs
  cl <- survival::clogit(Y ~ X + strata(id), data = d)                      # conditional MLE
  pl <- glm(Y ~ X, data = d, family = binomial("logit"), control = glm.control(epsilon = 1e-12, maxit = 100))
  Vpl <- sandwich::vcovCL(pl, cluster = d$id, type = "HC0", cadjust = FALSE)   # firm-clustered
  pa <- avg_effect_inference(function(t) dlogis(drop(cbind(1, d$X) %*% t)) * t[2], unname(coef(pl)), Vpl, cluster = d$id)
  mov <- !(ave(d$Y, d$id) %in% c(0, 1))                                     # firm-months of firms whose Y changes
  me_true <- P5$beta * dlogis(d$alpha_true + P5$beta * d$X)                 # true partial effects (simulation truth)
  list(pooled = coef(pl)[[2]], pooled_ci = ci95(coef(pl)[[2]], sqrt(Vpl[2, 2])), pooled_ape = pa$est, pooled_ape_ci = pa$ci,
       clogit = coef(cl)[[1]], clogit_se = sqrt(vcov(cl)[1, 1]), clogit_ci = ci95(coef(cl)[[1]], sqrt(vcov(cl)[1, 1])),
       fe = coef(fe)[[1]], fe_se = sqrt(vcov(fe)[1, 1]), bc = coef(bc)[[1]], bc_se = sqrt(vcov(bc)[1, 1]),
       ape_fe = ap$ape_fe, ape_fe_se = ap$se_fe, ape_fe_ci = ap$ci_fe, ape_bc = ap$ape_bc, ape_bc_se = ap$se_bc, ape_bc_ci = ap$ci_bc,
       share_kept = ap$share_kept, kept_part = sum(ap$g - me_true * mov) / nrow(d), drop_part = -sum(me_true[!mov]) / nrow(d))
}
set.seed(1205); d5 <- gen_B5(P5)
s8 <- fit_panel(d5); s4 <- fit_panel(d5[d5$t <= 4, ])                      # months 1 to 4 of the same panel
for (s in list(list("T = 8", s8), list("months 1-4", s4))) { o <- s[[2]]
  cat(sprintf("B5 (%s): pooled %.4f [%.4f, %.4f]; clogit %.4f [%.4f, %.4f]; FE MLE %.4f; BC %.4f\n", s[[1]],
              o$pooled, o$pooled_ci[1], o$pooled_ci[2], o$clogit, o$clogit_ci[1], o$clogit_ci[2], o$fe, o$bc))
  cat(sprintf("    APE over all firm-months: pooled %.4f [%.4f, %.4f]; FE %.4f [%.4f, %.4f]; BC %.4f [%.4f, %.4f]; share of firm-months kept %.3f\n",
              o$pooled_ape, o$pooled_ape_ci[1], o$pooled_ape_ci[2], o$ape_fe, o$ape_fe_ci[1], o$ape_fe_ci[2],
              o$ape_bc, o$ape_bc_ci[1], o$ape_bc_ci[2], o$share_kept))
}
cre_probit <- function(y, X, id) {                                          # pooled probit on (1, X_it, Xbar_i)
  X <- as.matrix(X); R <- cbind(1, X, apply(X, 2, function(x) ave(x, id)))
  m <- binary_mle(y, R, "probit")
  V <- m$vcov_model %*% crossprod(rowsum(m$scores, id)) %*% m$vcov_model   # firm-clustered sandwich
  a <- avg_effect_inference(function(t) dnorm(drop(R %*% t)) * t[2], m$coef, V, cluster = id)
  list(beta = m$coef, ape = a$est, ape_ci = a$ci)
}
cre <- cre_probit(d5$Y, cbind(X = d5$X), d5$id)
d5$Xbar <- ave(d5$X, d5$id)
gm <- lme4::glmer(Y ~ X + Xbar + (1 | id), data = d5, family = binomial("probit"), nAGQ = 8)   # CRE RE probit
sa <- as.numeric(attr(lme4::VarCorr(gm)$id, "stddev")); bg <- lme4::fixef(gm) / sqrt(1 + sa^2)
cat(sprintf("    CRE probit (T = 8): X %.4f, Xbar %.4f; APE %.4f [%.4f, %.4f]; CRE RE probit APE %.4f (sigma_a %.3f)\n",
            cre$beta[2], cre$beta[3], cre$ape, cre$ape_ci[1], cre$ape_ci[2],
            mean(dnorm(drop(cbind(1, d5$X, d5$Xbar) %*% bg))) * bg[2], sa))
# B5 (T = 8): pooled 1.0006 [0.9420, 1.0593]; clogit 0.8006 [0.7332, 0.8681]; FE MLE 0.9310; BC 0.7995
#     APE over all firm-months: pooled 0.1888 [0.1817, 0.1960]; FE 0.1236 [0.1143, 0.1328]; BC 0.1247 [0.1155, 0.1339]; share of firm-months kept 0.802
# B5 (months 1-4): pooled 0.9867 [0.9071, 1.0663]; clogit 0.7722 [0.6689, 0.8756]; FE MLE 1.0692; BC 0.7574
#     APE over all firm-months: pooled 0.1869 [0.1773, 0.1966]; FE 0.1208 [0.1080, 0.1337]; BC 0.1241 [0.1115, 0.1367]; share of firm-months kept 0.628
#     CRE probit (T = 8): X 0.4096, Xbar 0.5062; APE 0.1235 [0.1140, 0.1329]; CRE RE probit APE 0.1238 (sigma_a 0.547)
```

The subpanel of the first 4 months has the same true values as the 8-month panel: the pooled coefficient, 0.9867, aims at 0.977, and the interval for its APE, 0.1869, does not cover 0.1245; the APE intervals of FE and of the corrected FE both cover it. On a single panel the FE MLE with dummies is 1.0692 with 4 months and 0.9310 with 8 months, in line with the pattern the Monte Carlo below shows.

```r
library(survival)                                                            # clogit() inside fit_panel() needs strata() and coxph()
mc_one <- function(r, T) {                    # one simulated panel; seeds 1254000 + r (T = 4) and 1258000 + r (T = 8)
  set.seed((if (T == 4) 1254000 else 1258000) + r)
  f <- fit_panel(gen_B5(P5, T = T))
  c(fe = f$fe, fe_se = f$fe_se, bc = f$bc, bc_se = f$bc_se, cl = f$clogit, cl_se = f$clogit_se, ape_fe = f$ape_fe,
    ape_fe_se = f$ape_fe_se, ape_bc = f$ape_bc, ape_bc_se = f$ape_bc_se, kept = f$kept_part, drop = f$drop_part)
}
cover <- function(est, se, v) mean(abs(est - v) <= qnorm(0.975) * se)
for (T in c(4, 8)) {
  M <- do.call(rbind, parallel::mclapply(1:1000, mc_one, T = T, mc.cores = 8))
  ape_true <- if (T == 8) truth_B5$ape else truth_B5$ape_T4                 # quadrature: 0.124527 and 0.121431
  cat(sprintf("  T%d: bias of beta-hat FE %+.4f (%+.1f%%), BC %+.4f, CMLE %+.4f; coverage FE %.3f, BC %.3f, CMLE %.3f\n", T,
              mean(M[, "fe"]) - 0.8, 100 * (mean(M[, "fe"]) - 0.8) / 0.8, mean(M[, "bc"]) - 0.8, mean(M[, "cl"]) - 0.8,
              cover(M[, "fe"], M[, "fe_se"], 0.8), cover(M[, "bc"], M[, "bc_se"], 0.8), cover(M[, "cl"], M[, "cl_se"], 0.8)))
  cat(sprintf("      APE (truth %.4f): FE %+.2f%%, BC %+.2f%%; FE APE error = kept firms %+.4f + dropped firms %+.4f; coverage FE %.3f, BC %.3f\n",
              ape_true, 100 * (mean(M[, "ape_fe"]) / ape_true - 1), 100 * (mean(M[, "ape_bc"]) / ape_true - 1),
              mean(M[, "kept"]), mean(M[, "drop"]), cover(M[, "ape_fe"], M[, "ape_fe_se"], ape_true), cover(M[, "ape_bc"], M[, "ape_bc_se"], ape_true)))
}
#   T4: bias of beta-hat FE +0.3110 (+38.9%), BC -0.0153, CMLE +0.0021; coverage FE 0.008, BC 0.968, CMLE 0.951
#       APE (truth 0.1214): FE -0.70%, BC +2.70%; FE APE error = kept firms +0.0330 + dropped firms -0.0339; coverage FE 0.916, BC 0.912
#   T8: bias of beta-hat FE +0.1296 (+16.2%), BC -0.0015, CMLE -0.0004; coverage FE 0.078, BC 0.962, CMLE 0.952
#       APE (truth 0.1245): FE -0.22%, BC +0.57%; FE APE error = kept firms +0.0125 + dropped firms -0.0128; coverage FE 0.947, BC 0.943
```

The sum of the two terms of the error decomposition differs from the bias of the FE APE only by the sampling fluctuation of each panel's true APE around the population truth; with 4 months this term is below $10^{-4}$.

### 6.5 Blocks B6 and B7: Ordered and Multinomial Response

```r
set.seed(1206); d6 <- gen_B6(P6)
f6 <- ordinal::clm(rating ~ x1 + bill_kusd + G, scale = ~ G, data = d6, link = "probit")   # managed accounts get their own scale
th <- coef(f6); Xo <- cbind(d6$x1, d6$bill_kusd, d6$G)            # th: four cut points, three slopes, log scale of G
top_i <- function(t, s) t[5] / s * dnorm((t[4] - drop(Xo %*% t[5:7])) / s)        # effect of x1 on P(5 stars)
top <- avg_effect_inference(function(t) top_i(t, exp(t[8] * d6$G)), unname(th), vcov(f6))
f6h <- ordinal::clm(rating ~ x1 + bill_kusd + G, data = d6, link = "probit")       # no scale part (= MASS::polr)
toph <- avg_effect_inference(function(t) top_i(t, 1), unname(coef(f6h)), vcov(f6h))
lr_scale <- 2 * (as.numeric(logLik(f6)) - as.numeric(logLik(f6h)))
cat(sprintf("B6: APE of engagement on P(5) %.4f [%.4f, %.4f]; homoskedastic %.4f; delta_G %.3f, sigma_managed %.3f, LR (scale) %.1f\n",
            top$est, top$ci[1], top$ci[2], toph$est, th[7], exp(th[8]), lr_scale))
fit_B6 <- function(d) {                                           # the same estimator on a fresh sample
  f <- ordinal::clm(rating ~ x1 + bill_kusd + G, scale = ~ G, data = d, link = "probit"); X <- cbind(d$x1, d$bill_kusd, d$G)
  a <- avg_effect_inference(function(t) t[5] / exp(t[8] * d$G) * dnorm((t[4] - drop(X %*% t[5:7])) / exp(t[8] * d$G)),
                            unname(coef(f)), vcov(f))
  c(est = a$est, se = a$se, lo = a$ci[1], hi = a$ci[2])
}
M6 <- do.call(rbind, lapply(1:1000, function(r) { set.seed(1260000 + r); fit_B6(gen_B6(P6)) }))
cat(sprintf("    Monte Carlo (R = %d): mean APE %.5f, sd %.5f, mean SE %.5f, coverage %.3f\n", nrow(M6), mean(M6[, "est"]),
            sd(M6[, "est"]), mean(M6[, "se"]), mean(M6[, "lo"] <= truth_B6$ape_top & truth_B6$ape_top <= M6[, "hi"])))
cat(sprintf("    log of E[g | managed] / E[g | self-serve] at fixed covariates %.3f (truth %.3f)\n",
            -th[7] + (exp(th[8])^2 - 1) / 2, -P6$delta_G + (P6$sigma[2]^2 - P6$sigma[1]^2) / 2))
# B6: APE of engagement on P(5) 0.1118 [0.1023, 0.1214]; homoskedastic 0.1148; delta_G 0.210, sigma_managed 1.561, LR (scale) 140.8
#     Monte Carlo (R = 1000): mean APE 0.10063, sd 0.00491, mean SE 0.00487, coverage 0.953
#     log of E[g | managed] / E[g | self-serve] at fixed covariates 0.508 (truth 0.425)
```

In the Monte Carlo the standard errors average 0.00487, almost the same as the standard deviation of the estimates, 0.00491; the failure of this sample's interval to cover the truth is the kind of unlucky draw described in Section 3.9.

```r
set.seed(1207); d7 <- gen_B7(P7)
set.seed(1217); d7n <- gen_B7(P7, lambda = P7$lambda_nested)                 # Pro and Enterprise in one nest
mnl <- nnet::multinom(choice ~ x + s, data = d7, trace = FALSE, reltol = 1e-12, maxit = 1000, Hess = TRUE)
print(round(coef(mnl), 4))
cat(sprintf("B7 APE of engagement on P(Free, Pro, Ent): %s\n",
            paste(sprintf("%.4f", mnl_ape(coef(mnl), cbind(1, d7$x, d7$s), 2)), collapse = ", ")))
hm_iia <- hm_test(choice ~ x + s, d7, "Free"); hs_iia <- hm_test_stacked(hm_iia)       # drop Free: Pro vs Enterprise
mc_one <- function(r) {                                          # a fresh IIA sample of 5,000 firms
  set.seed(1270000 + r); h <- hm_test(choice ~ x + s, gen_B7(P7), "Free"); hs <- hm_test_stacked(h)
  c(stat = h$stat, p = h$p, min_eigen = h$min_eigen, s_p = hs$p)
}
M7 <- do.call(rbind, parallel::mclapply(1:4000, mc_one, mc.cores = 8))
cat(sprintf("B7 IIA sample: HM %.3f (p %.3f); Monte Carlo (R = %d) rejection at 1/5/10%%: %.4f %.4f %.4f; negative %.4f; V_r - V_f not PSD %.4f\n",
            hm_iia$stat, hm_iia$p, nrow(M7), mean(M7[, "p"] < 0.01), mean(M7[, "p"] < 0.05), mean(M7[, "p"] < 0.10),
            mean(M7[, "stat"] < 0), mean(M7[, "min_eigen"] < 0)))
cat(sprintf("    stacked-score: %.3f (p %.3f); rejection at 1/5/10%%: %.4f %.4f %.4f\n",
            hs_iia$stat, hs_iia$p, mean(M7[, "s_p"] < 0.01), mean(M7[, "s_p"] < 0.05), mean(M7[, "s_p"] < 0.10)))
log10p <- function(h) pchisq(h$stat, h$df, lower.tail = FALSE, log.p = TRUE) / log(10)
hm_nest <- hm_test(choice ~ x + s, d7n, "Free"); hs_nest <- hm_test_stacked(hm_nest)
mnl_n <- nnet::multinom(choice ~ x + s, data = d7n, trace = FALSE, reltol = 1e-12, maxit = 1000)
gpaid <- glm(I(choice == "Enterprise") ~ x + s, data = d7n[d7n$choice != "Free", ], family = binomial(),
             control = glm.control(epsilon = 1e-12))                                      # paying firms only
cat(sprintf("B7 nested sample: HM %.2f (log10 p %.2f); stacked-score %.2f (log10 p %.2f); Ent - Pro on x: full MNL %.4f, paying-firm logit %.4f\n",
            hm_nest$stat, log10p(hm_nest), hs_nest$stat, log10p(hs_nest),
            coef(mnl_n)["Enterprise", "x"] - coef(mnl_n)["Pro", "x"], coef(gpaid)[["x"]]))
#            (Intercept)      x      s
# Pro             0.6579 0.9972 0.2130
# Enterprise      0.0186 0.3666 1.0123
# B7 APE of engagement on P(Free, Pro, Ent): -0.1297, 0.1740, -0.0443
# B7 IIA sample: HM 3.151 (p 0.369); Monte Carlo (R = 4000) rejection at 1/5/10%: 0.0393 0.0700 0.1010; negative 0.1125; V_r - V_f not PSD 0.4128
#     stacked-score: 3.240 (p 0.356); rejection at 1/5/10%: 0.0083 0.0423 0.0855
# B7 nested sample: HM 28.72 (log10 p -5.59); stacked-score 30.97 (log10 p -6.06); Ent - Pro on x: full MNL -1.3289, paying-firm logit -1.5177
```

Of the six MNL coefficients, Pro's seats coefficient, 0.213, has an interval, $[0.135, 0.291]$, that misses the true value 0.3, while the other five all cover theirs. Nominal 95% intervals are bound to miss now and then, and one miss in six is not unusual; B7 has no Monte Carlo for the coefficients, and the table below lists this row anyway. The APE of engagement on the Pro probability is 0.1740, and on Free it is -0.1297.

### 6.6 The Grand Reconciliation

| Block | Estimator | Estimate | 95% interval | Target and true value | Covered |
|---|---|---|---|---|---|
| B1 | LPM slope | 0.2153 | [0.2037, 0.2269] | APE 0.2156 | Yes |
| B1 | Logit slope | 1.0969 | [1.0103, 1.1835] | 1.10 | Yes |
| B1 | Logit APE | 0.2171 | [0.2054, 0.2289] | 0.2156 | Yes |
| B2 | Average cross-partial effect | -0.0161 | [-0.0272, -0.0051] | -0.0194 | Yes |
| B3 | Naive probit APE | 0.2410 | [0.2325, 0.2496] | 0.1475 | No |
| B3 | Control-function APE | 0.1659 | [0.1370, 0.1949] | 0.1475 | Yes |
| B3 | IV probit $\beta_D$ | 0.583 | [0.438, 0.728] | 0.5 | Yes |
| B4 | Short-probit coefficient | 0.638 | [0.590, 0.685] | $\beta/\sqrt{1 + \sigma_a^2}$ 0.636 | Yes |
| B4 | Short-probit APE | 0.2135 | [0.2021, 0.2250] | 0.2127 | Yes |
| B5 | Pooled logit (8 months) | 1.0006 | [0.942, 1.059] | $\beta$ 0.8; probability limit 0.977 | No; yes |
| B5 | FE logit CMLE (8 months) | 0.8006 | [0.733, 0.868] | 0.8 | Yes |
| B5 | FE logit CMLE (first 4 months) | 0.7722 | [0.669, 0.876] | 0.8 | Yes |
| B5 | CRE probit APE (8 months) | 0.1235 | [0.1140, 0.1329] | 0.1245 | Yes |
| B5 | Pooled logit APE (8 months) | 0.1888 | [0.1817, 0.1960] | 0.1245 | No |
| B6 | APE on 5 stars | 0.1118 | [0.1023, 0.1214] | 0.1008 | No |
| B7 | Pro's seats coefficient | 0.213 | [0.135, 0.291] | 0.3 | No |
| B7 | Enterprise's engagement coefficient | 0.367 | [0.278, 0.455] | 0.4 | Yes |

The five rows that fail to cover fall into two kinds. B3's naive APE and B5's pooled estimates miss because of bias, not bad luck: the former ignores endogeneity, the latter the correlation between $\alpha_i$ and the regressors, and the pooled logit interval does cover its own probability limit. The rows of B6 and B7 are single-sample misses; B6's Monte Carlo coverage rate is 0.953, and conclusions about coverage are drawn only from the Monte Carlo. B1's LPM interval uses the HC1 standard error and aims at the slope of the population linear projection, which with normal regressors is the APE.

The reconciliation of this section comes down to one sentence: every estimator aims at its own population quantity, and where things go wrong, either a coefficient has been read as an effect or something omitted has been allowed to be correlated with the regressors.

## 7 Failure Modes and Robustness

The first kind of failure is numerical, most commonly the perfect separation of Section 4.1, followed by a nonconcave likelihood that stops at a false peak. The prescription is plain: check whether the coefficients and fitted probabilities are obviously abnormal, try several starting values, and use the two-step estimates as the starting point for the full MLE; when bootstrapping, count how many resamples fail to converge.

The second is mistaking coefficients for effects. The scale of the coefficients is set by the variance of the latent error: when two groups of firms have different heterogeneity variances (Theorem 4) or the latent error is heteroskedastic (Section 3.6), comparing the two groups' probit coefficients compares scales; the sign of an interaction coefficient is not the sign of the cross-partial effect, 0.221 against -0.0161 on B2; nor do MNL coefficients determine the direction of the partial effects, 0.367 against -0.0443 on B7. The remedy is to report APEs and, where necessary, to plot the observation-by-observation partial effects. Ordered response adds one more layer: Bond and Lang (2019) show that without the assumption that "the latent distributions of the two groups have the same scale," neither the mean nor the distribution of star ratings can rank latent satisfaction; B6's managed accounts give more 5-star ratings and have a higher average rating, yet under a reasonable monotone transformation they are less satisfied.

The third is endogeneity. The control function relies on the normality of the first-stage error and on $D_i$ being continuous; when $D_i$ is binary, a bivariate probit is needed, and the shortcut of plugging fitted probabilities into a probit is inconsistent. With weak instruments, the problems discussed in Section 3.5 of Chapter 9 (Theorem 5) arise in nonlinear models just the same, and the strength of the first stage has to be reported as usual; B3's first-stage F is 642.2, so this is no issue here. A very cheap cross-check is 2SLS on the LPM, which does not rely on normality; on B3 it gives 0.1559, close to the control function's 0.1659, and the intervals of both cover the truth.

The fourth is panels. A pooled estimator that ignores $\alpha_i$ absorbs the correlation between the heterogeneity and the regressors, and on B5 its APE is too high by half; even when only the APE is wanted, the variance must be clustered by unit. FE logit identifies $\beta$ but cannot deliver the APE; the FE MLE with dummies has a bias of order $1/T$, and the bias correction fixes the coefficient (leaving -0.015 with 4 months) but does not improve the APE, whose +2.70% is not robust. The all-firm FE APE has a bias of less than 1% on B5, but it owes this to two large errors offsetting each other and does not support the general claim that "the FE APE is always accurate"; with 4 months the coverage rate of the APE interval is only 0.916, so short panels call for caution. CRE needs the form of the distribution of $\alpha_i$ given the regressors, and B5 happens to satisfy Mundlak's linear form; when that form does not hold, the APE that CRE estimates is biased as well.

The fifth is IIA in multinomial response and its test. IIA is not credible when alternatives are similar. Wooldridge (2010) p. 648 attributes the red bus and blue bus example to McFadden (1974): commuters originally split evenly between the car and the red bus; when a blue bus that differs only in color is added, IIA forces the car's share to fall from one half to one third, whereas common sense says that bus riders simply split between the two buses. In the nested version of B7, the full-sample MNL estimates the engagement coefficient among paying firms at -1.329, against a true value of -1.5. The Hausman-McFadden test can detect this, but the classical statistic behaves poorly in finite samples: when IIA holds, across 4,000 replications it rejects 7.0% of the time at the 5% level, 5.8 Monte Carlo standard errors above the nominal level (standard error 0.0034), and 3.9% of the time at the 1% level. The cause lies in $\widehat{\mathrm{Avar}}(\hat\beta_{\mathrm{sub}}) - \widehat{\mathrm{Avar}}(\hat\beta_{\mathrm{all}})$: dropping one category loses little information, so the two variances are almost equal and their difference is close to singular; the median ratio of the smallest to the largest eigenvalue is only 0.004, the difference is not even positive semidefinite in 41.3% of the replications, and 11.25% of the statistics are negative. Stacking the scores of the two estimators to compute the variance of the difference (similar to Stata's suest) yields a variance matrix that is positive semidefinite by construction, with rejection rates of 0.8%, 4.2%, and 8.6% at the 1%, 5%, and 10% levels, close to nominal, and it rejects just as clearly on the nested sample (30.97, with a $p$ value whose common logarithm is -6.06). In practice, prefer this version, or estimate a nested logit directly and test $\lambda = 1$.

Taken together, these failure modes show that the robustness of this chapter's methods rests on two things: first, being clear about whether you want a parameter of latent utility or an effect on the probability scale, since the former depends on the scale normalization and only the latter can be compared across models, groups, and specifications; second, thinking through which unobserved factors are correlated with the regressors, since independent heterogeneity only attenuates the coefficients, while correlated heterogeneity and endogeneity make the effects themselves go wrong.

## 8 Further Reading

::: {.readings}
Required reading, in suggested reading order:

- Wooldridge (2010, *Econometric Analysis of Cross Section and Panel Data*, 2nd ed.), §15.2 to §15.8 and §16.2 to §16.3. The equation numbers and pages of Theorems 1 to 5, 7, 8, 10, and 11 of this chapter all come from here, and Table 15.4 (p. 625) lays out the estimators of panel binary models by their assumptions in a single table.
- McFadden (1974, *Frontiers in Econometrics*). The source of random utility and conditional logit, the template for translating the economic model of discrete choice into a likelihood, and the origin of the adoption model of Section 2 of this chapter.
- Ai and Norton (2003, Economics Letters) and Norton, Wang and Ai (2004, Stata Journal). The coefficient of an interaction term is not the interaction effect; the latter gives formulas, standard errors, and examples.
- Fernández-Val and Weidner (2016, Journal of Econometrics). Bias correction and APEs in nonlinear panels, and the source of the methods in `alpaca`.

Further reading:

- Rivers and Vuong (1988, Journal of Econometrics). The two-step control function and the exogeneity test for a continuous endogenous regressor.
- Chamberlain (1980, Review of Economic Studies). The source of the conditional likelihood of FE logit and of CRE.
- Neyman and Scott (1948, Econometrica) and Lancaster (2000, Journal of Econometrics). The original source of the incidental parameters problem, and a survey of it.
- Hahn and Newey (2004, Econometrica) and Fernández-Val (2009, Journal of Econometrics). Analytical and jackknife bias corrections for fixed effects MLE, and the APE in panel probit.
- Wooldridge (2005b, Journal of Applied Econometrics). The initial conditions problem in dynamic nonlinear panels.
- Lewbel (2000, Journal of Econometrics). The special regressor.
- Horrace and Oaxaca (2006, Economics Letters). The bias of OLS under the LPM and the trimming estimator.
- McKelvey and Zavoina (1975, Journal of Mathematical Sociology) and Bond and Lang (2019, Journal of Political Economy). The latent-variable model of the ordered probit, and the warning that comparisons on a latent scale are not identified.
- Hausman and McFadden (1984, Econometrica). A Hausman-type test of IIA.
- Greene (2010, Economics Letters) and Puhani (2012, Economics Letters). Further discussion of tests on interaction terms and of treatment effects in nonlinear DiD.
- Forman (2005, Management Science). Internet adoption decisions of a large sample of establishments in finance and services in 1998, an early example of binary adoption research in the IS literature.
- Levine and Seiler (2023, Marketing Science). Using stockouts caused by hurricanes to test for state dependence in brand choice.
- Stammann, Heiss and McFadden (2016). The source of the method in `bife`, and computation of FE logit in large panels.
:::

::: {.apa-refs}
- Ai, C., & Norton, E. C. (2003). Interaction terms in logit and probit models. *Economics Letters, 80*(1), 123-129. https://doi.org/10.1016/S0165-1765(03)00032-6
- Bond, T. N., & Lang, K. (2019). The sad truth about happiness scales. *Journal of Political Economy, 127*(4), 1629-1640. https://doi.org/10.1086/701679
- Chamberlain, G. (1980). Analysis of covariance with qualitative data. *The Review of Economic Studies, 47*(1), 225-238. https://doi.org/10.2307/2297110
- Fernández-Val, I. (2009). Fixed effects estimation of structural parameters and marginal effects in panel probit models. *Journal of Econometrics, 150*(1), 71-85. https://doi.org/10.1016/j.jeconom.2009.02.007
- Fernández-Val, I., & Weidner, M. (2016). Individual and time effects in nonlinear panel models with large N, T. *Journal of Econometrics, 192*(1), 291-312. https://doi.org/10.1016/j.jeconom.2015.12.014
- Forman, C. (2005). The corporate digital divide: Determinants of Internet adoption. *Management Science, 51*(4), 641-654. https://doi.org/10.1287/mnsc.1040.0343
- Greene, W. (2010). Testing hypotheses about interaction terms in nonlinear models. *Economics Letters, 107*(2), 291-296. https://doi.org/10.1016/j.econlet.2010.02.014
- Hahn, J., & Newey, W. (2004). Jackknife and analytical bias reduction for nonlinear panel models. *Econometrica, 72*(4), 1295-1319. https://doi.org/10.1111/j.1468-0262.2004.00533.x
- Hausman, J., & McFadden, D. (1984). Specification tests for the multinomial logit model. *Econometrica, 52*(5), 1219-1240. https://doi.org/10.2307/1910997
- Horrace, W. C., & Oaxaca, R. L. (2006). Results on the bias and inconsistency of ordinary least squares for the linear probability model. *Economics Letters, 90*(3), 321-327. https://doi.org/10.1016/j.econlet.2005.08.024
- Lancaster, T. (2000). The incidental parameter problem since 1948. *Journal of Econometrics, 95*(2), 391-413. https://doi.org/10.1016/S0304-4076(99)00044-5
- Levine, J., & Seiler, S. (2023). Identifying state dependence in brand choice: Evidence from hurricanes. *Marketing Science, 42*(5), 934-957. https://doi.org/10.1287/mksc.2022.1415
- Lewbel, A. (2000). Semiparametric qualitative response model estimation with unknown heteroscedasticity or instrumental variables. *Journal of Econometrics, 97*(1), 145-177. https://doi.org/10.1016/S0304-4076(00)00015-4
- McFadden, D. (1974). Conditional logit analysis of qualitative choice behavior. In P. Zarembka (Ed.), *Frontiers in econometrics* (pp. 105-142). Academic Press.
- McKelvey, R. D., & Zavoina, W. (1975). A statistical model for the analysis of ordinal level dependent variables. *The Journal of Mathematical Sociology, 4*(1), 103-120. https://doi.org/10.1080/0022250X.1975.9989847
- Neyman, J., & Scott, E. L. (1948). Consistent estimates based on partially consistent observations. *Econometrica, 16*(1), 1-32. https://doi.org/10.2307/1914288
- Norton, E. C., Wang, H., & Ai, C. (2004). Computing interaction effects and standard errors in logit and probit models. *The Stata Journal, 4*(2), 154-167. https://doi.org/10.1177/1536867X0400400206
- Puhani, P. A. (2012). The treatment effect, the cross difference, and the interaction term in nonlinear "difference-in-differences" models. *Economics Letters, 115*(1), 85-87. https://doi.org/10.1016/j.econlet.2011.11.025
- Rivers, D., & Vuong, Q. H. (1988). Limited information estimators and exogeneity tests for simultaneous probit models. *Journal of Econometrics, 39*(3), 347-366. https://doi.org/10.1016/0304-4076(88)90063-2
- Stammann, A., Heiss, F., & McFadden, D. (2016). *Estimating fixed effects logit models with large panel data* [Conference paper]. Beiträge zur Jahrestagung des Vereins für Socialpolitik 2016: Demographischer Wandel, Session: Microeconometrics, G01-V3. ZBW. http://hdl.handle.net/10419/145837
- Wooldridge, J. M. (2005b). Simple solutions to the initial conditions problem in dynamic, nonlinear panel data models with unobserved heterogeneity. *Journal of Applied Econometrics, 20*(1), 39-54. https://doi.org/10.1002/jae.770
- Wooldridge, J. M. (2010). *Econometric analysis of cross section and panel data* (2nd ed.). MIT Press.
:::
