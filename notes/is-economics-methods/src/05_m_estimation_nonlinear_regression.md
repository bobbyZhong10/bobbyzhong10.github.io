---
title: "M-Estimation & Nonlinear Regression"
subtitle: "The Extremum Estimator and Its Sandwich"
seriesline: "Foundations of Information Systems Economics · Chapter 5"
series: "Foundations of Information Systems Economics: Empirical Methods for Information Systems Economics"
sidebartitle: "Foundations of Information Systems Economics"
sidebarsub: "Theory · Methods · Applications"
chapterlabel: "Chapter 5 · M-Estimation & Nonlinear Regression"
---

## Introduction

Larkspur is a fictional online marketplace where small merchants sell home decor, apparel, jewelry, and stationery. The platform picks half of its merchants at random and switches on Spotlight, a shop-promotion feature, for them. After 30 days it wants to answer a question that could hardly be more ordinary: by what percentage did Spotlight raise merchants' sales? The trouble is that without Spotlight, about six in ten merchants fail to sell a single item in a month. With so many zeros the outcome cannot simply be logged, so some analysts add 1 to sales before taking logs, while others look only at merchants with positive sales. Section 1 shows that the first approach gives an answer that depends on whether sales are recorded in dollars or in cents.

This is the first chapter of Part II, and it asks how far a single idea can go: an estimator is the minimizer of some sample objective function. OLS (ordinary least squares), NLS (nonlinear least squares), Poisson regression, maximum likelihood, and quantile regression are all special cases of M-estimation: choose a loss $q(w,\theta)$ for a single observation, average it over the sample, and find the $\theta$ that makes the average loss smallest. The payoff of unifying them is not tidy notation. It is that three questions need answering only once: what the population minimizer is (identification), whether the sample minimizer comes close to it (consistency), and how it fluctuates around the population minimizer (asymptotic normality and the sandwich variance). The answer to the first question also tells us something that is often overlooked: when the model is misspecified, the estimator still converges; it just converges to the best approximation that the loss function itself picks out.

The chapter proceeds as follows. Section 2 writes down Spotlight's causal quantity and the estimands of Larkspur's other three data blocks. Section 3, the most detailed in the chapter, first gives the identification, consistency, and asymptotic normality theorems of M-estimation, then applies them to three concrete estimators: NLS, two-step estimation, and Poisson PML (Poisson pseudo maximum likelihood, PPML from here on). It also explains why the coefficient of log(1+y) is not a percentage effect. Section 4 covers algorithms and the choice of variance estimator, Section 5 examines three anchoring papers, Section 6 runs the whole toolkit on Larkspur, and Section 7 discusses failure modes. The objective function $q$, the sandwich's bread $A_0$ and filling $B_0$, and Wooldridge (2010) Theorem 12.3, all set up in this chapter, return in Chapter 6 (maximum likelihood), Chapter 7 (GMM, the generalized method of moments), Chapter 8 (inference), and Chapters 12 and 13.

## 1 An Effect That Changes with the Unit of Measurement

::: {.case}
Larkspur has 20,000 merchants, and Spotlight is assigned to exactly half of them by complete randomization. The outcome is sales in the 30 days after assignment, which the platform books in whole cents. The share of merchants with zero sales is 0.606 in the control group and falls to 0.505 in the Spotlight group. Average sales in the two groups are USD 175.25 and USD 237.82. Merchants with positive sales sell USD 465.07 on average, with a median of USD 229.43, and the top merchant sells USD 19,839.37.

Because this is a simulation, we know the answer: Spotlight multiplies average sales by 1.381, a rise of 38.1%, or 0.323 in logs. The effect comes through two channels: more merchants make a sale (the probability of selling anything is multiplied by 1.25), and merchants who were already selling sell 10.5% more.
:::

As many papers do, the analyst adds 1 to sales, takes logs, and runs OLS on a Spotlight indicator with category and tenure controls. The coefficient is 0.593 (robust standard error 0.039), read as "Spotlight raised sales by about 59%," which is more than half again as large as the truth. A colleague who checks the work uses the finance system's records in cents, and the same regression gives 1.058; in thousands of dollars, it gives 0.041. Same randomized experiment, same merchants, and three answers that differ by a factor of more than twenty. This is not sampling error: in a simulated population of one million merchants, the population counterparts of these three regressions are 0.583, 1.038, and 0.041 (figure below).

![How Spotlight's "effect" changes with the unit in which sales are measured: the coefficient of log(1+y) regressed on Spotlight with its 95% interval and its population counterpart, against the unit-free PPML interval and the true log ratio. The horizontal axis is the accounting unit, from USD 1,000 down to a ten-thousandth of a dollar.](assets/fig/fig_05_log1p_units.svg)

Change the estimator and the unit problem disappears. PPML, with the same right-hand-side variables, gives a coefficient of 0.302 whether sales are booked in dollars or in cents; the two differ on the order of $10^{-16}$. The 95% interval [0.224, 0.379] covers the truth, 0.323, and as a proportional effect the estimate is a rise of 35.2%. This chapter answers three questions: what exactly PPML estimates, and why it remains valid on data with many zeros; why the coefficient of log(1+y) moves with the unit; and how the standard errors of both should be computed. Section 3.8 works out the difference between the log(1+y) coefficients in dollars and in cents: it almost exactly equals the extensive-margin effect times $\log 100$.

## 2 The Economic Model and the Estimand

### 2.1 Spotlight's Two Margins and the Proportional Effect

Let $y_i(1)$ and $y_i(0)$ be the 30-day sales (in dollars) of merchant $i$ with and without Spotlight, and let $D_i$ be the random indicator of whether Spotlight is switched on. What we observe is $y_i = D_i\,y_i(1) + (1 - D_i)\,y_i(0)$. This is the first time potential outcomes appear in the book. Every merchant has two, and we see only one of them. Random assignment makes the two groups comparable before Spotlight is switched on, so the difference in group means is the average causal effect. In Larkspur's DGP (data-generating process), Spotlight acts on two margins at once. One is the extensive margin: the probability of making a sale rises from 0.395 to 0.494, a difference of 0.099. The other is the intensive margin: for merchants who were selling anyway, log sales rise by 0.10. Neither number is what management wants to know. Their question is "by what percentage did sales go up?"

The target of this chapter is the proportional effect on average sales:

$$\log\frac{\mathbb{E}[y_i(1)]}{\mathbb{E}[y_i(0)]} = \log 1.25 + 0.10 = 0.323,$$

that is, average sales rise by

$$\frac{\mathbb{E}[y_i(1)]}{\mathbb{E}[y_i(0)]} - 1 = 38.1\%.$$

This target has three advantages. First, it does not depend on the unit: multiply all sales by 100 and the numerator and denominator scale up together, leaving the ratio unchanged. Second, zeros count: a merchant who goes from 0 to USD 300 and another who goes from USD 300 to USD 600 make the same contribution, USD 300, to the average effect $\mathbb{E}[y_i(1)] - \mathbb{E}[y_i(0)]$. Third, numerator minus denominator is the average effect in levels: Larkspur's design values rise from USD 171.18 to USD 236.48, so the proportional effect converts directly back into dollars. The price should also be stated up front: it is "the percentage change in the average," not "the average of each merchant's percentage change." The latter is undefined for merchants who go from zero to positive sales. Section 3.8 shows that when zeros are present, no parameter can at once be an average of individual effects that respond to the intensive margin, be free of the unit, and be identified from the data.

How do the two margins combine into this number? In the DGP, the selling probability is multiplied by 1.25 and the sales of selling merchants by $e^{0.10}$. Because whether a merchant sells is independent of the shock to sales, the ratio is exactly the product of the two: in logs, 0.223 comes from the extensive margin and 0.10 from the intensive margin. For each category $c$ and standardized tenure, the conditional mean has the exponential form

$$\mathbb{E}[y_i(d) \mid c, \mathrm{tenure}] = \exp(\alpha_c + 0.323\,d + 0.25\,\mathrm{tenure}),$$

where $\alpha_c$ is the intercept of category $c$. This is why PPML is correctly specified on M1: what it estimates is precisely the coefficient of $D$ in an exponential mean.

### 2.2 Conditional Means, Loss Functions, and Estimands

Each of Larkspur's other three data blocks has a conditional mean to estimate, and each illustrates one theoretical point of the chapter. Block M2 has 2,000 active merchants. The conditional mean of monthly sales $y$ (in thousands of dollars) is $m(x,\beta) = \exp(\beta_0 + \beta_1 x_1 + \beta_2 x_2)$, where $x_1$ is the log of the number of listings relative to the median of 40, $x_2$ indicates free shipping, and the true value is $\beta = (\log 1.2,\ 0.45,\ 0.20)$. The error is multiplicative, and $\mathrm{Var}(y \mid x)$ is proportional to the square of the mean. Block M2m changes the true mean to $\exp(0.2 + 0.55x - 0.10x^2)$, with diminishing returns to listings, while the analyst fits only $\exp(\beta_0 + \beta_1 x)$; now there is no "true $\beta_1$" to estimate. Block M3 targets the elasticity of sales with respect to shop traffic, 0.8, but the traffic index is unobserved and must be estimated in a first step. Block M4 holds the monthly order counts of 2,000 merchants, with conditional mean $\exp(1.4 + 0.45x)$ ($x$ is a standardized rating) and the negative binomial NB2 variance $\mathrm{Var}(y \mid x) = \mathbb{E}[y \mid x] + \eta^2\,\mathbb{E}[y \mid x]^2$, $\eta^2 = 1/3$. The ratio of the average conditional variance to the average mean, $\mathbb{E}[\mathrm{Var}(y \mid x)]/\mathbb{E}[y]$, is 2.83 (the pointwise ratio $1 + \eta^2\,\mathbb{E}[y \mid x]$ varies with the rating).

The estimands of the four data blocks share one form, and that form is the core idea of the chapter: an estimand is the minimizer of some population loss. The coefficients of M1 and M4 minimize the expected Poisson loss, those of M2 minimize the mean squared error, and those of M3 minimize the mean squared error given the first step. When the mean is correctly specified, all these minimizers are the true parameters of the conditional mean, and it does not matter which loss is used. When it is misspecified, different losses pick out different best approximations: on M2m, the NLS and PPML slopes aim at 0.415 and 0.473. So "what the estimator estimates" cannot be separated from "what it minimizes," and Section 3 begins by turning this sentence into theorems.

## 3 Identification and Large-Sample Theory

This section first uses Larkspur's numbers or figures to show what is to be proved, then states numbered assumptions and theorems, followed by proof sketches of two to six steps; each step names the assumptions it uses, and regularity conditions are only cited to their source. After each theorem, an "In words" paragraph answers four questions: under which assumptions it estimates which quantity, what it estimates when the specification is wrong, who sets the weights, and what Larkspur's numbers say. Assumptions and theorems are numbered within this chapter; the parentheses give the theorem numbers, equation numbers, and page numbers of Wooldridge (2010). Sections 3.1 to 3.4 are general theory, Sections 3.5 to 3.7 cover three concrete estimators, and Section 3.8 discusses log(1+y) and transformations like it.

### 3.1 M-Estimation: A Unified Framework

Start with the intuition. Imagine a rolling terrain whose height at each point is set by the parameter $\theta$; we are looking for its lowest point. The population terrain $Q_0(\theta)$ cannot be seen; all we have is a noisy contour map $Q_n(\theta)$ drawn from the sample. The estimator is the lowest point of the sample terrain, and the estimand is the lowest point of the population terrain. The larger the sample, the more the two maps resemble each other and the closer the two lowest points are. That is the whole intuition of consistency.

Now turn this into a definition. Let $w_i$ collect all the variables of observation $i$, let the parameter be $\theta \in \Theta \subset \mathbb{R}^k$, and let $q(w,\theta)$ be the loss of a single observation. The M-estimator is

$$\hat\theta = \arg\min_{\theta \in \Theta} Q_n(\theta), \qquad Q_n(\theta) = \frac{1}{n}\sum_{i=1}^n q(w_i,\theta),$$

and its population counterpart is the point that minimizes $Q_0(\theta) = \mathbb{E}[q(w,\theta)]$, written $\theta_0$. Part II always writes problems as minimizations; an objective to be maximized, such as the log-likelihood, enters with its sign flipped. The commonly used estimators differ only in $q$:

| Estimator | Loss of a single observation $q(w,\theta)$ | Where this book develops it |
|---|---|---|
| OLS | $(y - x'\theta)^2$ | The familiar special case |
| NLS | $\big(y - m(x,\theta)\big)^2$ | Section 3.5 |
| Poisson PML | $\exp(x'\theta) - y\,x'\theta$ | Section 3.7 |
| MLE (maximum likelihood) | $-\ln f(w;\theta)$ | Section 3.2 and Chapter 6 |
| LAD (least absolute deviations) | $\lvert y - x'\theta \rvert$ | Section 3.4 and Chapter 13 |

GMM also belongs to the family of extremum estimators, except that its objective is not an average of single-observation losses but a quadratic form in sample moments:

$$\hat\theta_{GMM} = \arg\min_\theta\ \bar g_n(\theta)'\, W\, \bar g_n(\theta), \qquad \bar g_n(\theta) = \frac{1}{n}\sum_{i=1}^n g(w_i, \theta).$$

Its consistency and asymptotic normality follow the argument of this section, recast in the language of moment conditions; Chapter 7 develops it.

The payoff of this viewpoint is that once you prove, a single time, under what conditions an M-estimator is consistent and what form its variance takes, every estimator in the table is proved at once. It also has a more important corollary: $\theta_0$ is defined as "the minimizer of the population loss," not as "the parameter of the true model." Whether the two coincide depends on additional assumptions, such as a correctly specified conditional mean for NLS. When they do not, the theorems below still hold; the estimator simply converges to the point that $q$ itself picks out. Under correct specification this minimizer is the true parameter and is still written $\theta_0$; under misspecification it is written $\theta_*$ and called the pseudo-true value.

### 3.2 Identification: The Population Objective Must Have a Unique Valley Floor

Start with M2m. The true mean $\exp(0.2 + 0.55x - 0.10x^2)$ is not in the working model $\exp(\beta_0 + \beta_1 x)$, yet the mean squared error $\mathbb{E}[(\mathbb{E}[y \mid x] - \exp(\beta_0 + \beta_1 x))^2]$, as a function of $(\beta_0, \beta_1)$, still has a single lowest point, which numerical integration puts at $(0.170, 0.415)$; fitting a million noise-free points gives $(0.170, 0.416)$. Identification asks only whether the population problem has a unique solution, not whether the model is right.

::: {.assumption}
**M.1 (identification; Wooldridge 2010 equation (12.11), p. 402)** For all $\theta \in \Theta$ with $\theta \neq \theta_0$, $\mathbb{E}[q(w,\theta_0)] < \mathbb{E}[q(w,\theta)]$.
:::

M.1 is a restriction on the population; the sample size plays no part in it. When it fails, the data cannot distinguish two parameters even at the population level, and no sample, however large, would help: this is a failure of identification, not a shortfall of precision. Wooldridge (2010) §12.2 gives a typical example: in the mean model $\theta_1 + \theta_2 x_2 + \theta_3 x_3^{\theta_4}$, if the true $\theta_3 = 0$, every value of $\theta_4$ gives the same mean squared error, and the valley floor is a flat trench.

Why does $\theta_0$ land exactly on the valley floor? That depends on $q$ and on what is assumed about the data. For NLS the reason is that the conditional mean minimizes the mean squared error: if $\mathbb{E}[y \mid x] = m(x,\theta_0)$, write $\varepsilon = y - m(x,\theta_0)$; the cross term vanishes because $\mathbb{E}[\varepsilon \mid x] = 0$, so

$$\mathbb{E}\big[(y - m(x,\theta))^2\big] = \mathbb{E}[\varepsilon^2] + \mathbb{E}\big[(m(x,\theta_0) - m(x,\theta))^2\big],$$

and the second term is nonnegative, so $\theta_0$ is a minimizer; that it is the unique one is exactly assumption NLS.2 of Section 3.5. For MLE the reason is Jensen's inequality. Take $q = -\ln f(w;\theta)$, where $f(\cdot;\theta_0)$ is the true density of the data. Because $-\ln$ is convex,

$$\begin{aligned}
\mathbb{E}[q(w,\theta)] - \mathbb{E}[q(w,\theta_0)] &= \mathbb{E}\left[ -\ln \frac{f(w;\theta)}{f(w;\theta_0)} \right] \\
&\geq -\ln \mathbb{E}\left[ \frac{f(w;\theta)}{f(w;\theta_0)} \right] \\
&= -\ln \int f(w;\theta)\, dw = 0,
\end{aligned}$$

with equality if and only if the two densities are equal almost everywhere; if different values of $\theta$ give different densities, $\theta_0$ is the unique valley floor, which is M.1. The left side is the Kullback-Leibler (KL) divergence, which is nonnegative and zero only when the two distributions are the same. When the model family does not contain the true density, the minimizer is the member closest to the true distribution in KL, and this is where Chapter 6 starts its discussion of QMLE (quasi-maximum likelihood). For Poisson PML the reason is the most surprising: provided the conditional mean is correctly specified, and whatever the distribution of $y$, $\theta_0$ sits on the valley floor (Section 3.7 gives the proof).

Identification also comes in degrees, strong and weak. Even when the valley floor is unique, it may be very flat: a broad swath of $\theta$ values all give nearly the same low population loss. The steepness of the valley is measured by the expected Hessian (the matrix of second derivatives), $A_0 = \mathbb{E}[\nabla^2_\theta q(w,\theta_0)]$, and Section 3.4 shows that this is exactly the "bread" on both sides of the asymptotic variance.

::: {.intuition}
Picture the strength of identification as the steepness of the valley. A deep gorge with steep walls, where you could grope your way to the bottom with your eyes closed, is strong identification: move slightly away from the truth and the average loss rises noticeably, so the data have strong opinions about which $\theta$ is better. A shallow hollow sinking gently, where the bottom is anyone's guess, is weak identification: countless values of $\theta$ look equally good. The model above behaves this way when $\theta_3$ is close to zero: $\theta_4$ can then be pinned down only by the faint signal in $\theta_3 x_3^{\theta_4}$. For given gradient noise $B_0$, the steeper the valley, the larger $A_0$ and the smaller $A_0^{-1}B_0A_0^{-1}$. What really determines precision is the size of the curvature relative to the gradient noise (multiply $q$ by 2 and $A_0$ doubles while $B_0$ quadruples, leaving the variance unchanged; see Section 3.4). In a maximization problem, flip the sign of the objective and the valley turns back into a peak; the logic is the same.
:::

### 3.3 Consistency: How the Sample Valley Floor Is Pulled Toward the Population One

Begin by asking where the estimator aims on M4. In a Monte Carlo simulation that redraws M4 2,000 times, the PPML estimate of the rating coefficient averages 0.4507 against a truth of 0.45, a bias of 0.0007. As the ratio of the average conditional variance to the average mean is turned up from Poisson's 1 all the way to 8.33, the average estimate stays between 0.4488 and 0.4507 throughout. So the center of the sampling distribution sits at the population valley floor and does not move with the conditional variance. Consistency says that the estimator concentrates around this point as the sample grows, and that step is what the theorems below handle. The point estimate does not care what the variance looks like, only where the valley floor of the population loss lies.

The pointwise law of large numbers guarantees that at each fixed $\theta$, $Q_n(\theta) \xrightarrow{p} Q_0(\theta)$ ($\xrightarrow{p}$ denotes convergence in probability, $o_p(1)$ a term that converges to zero in probability, and $\xrightarrow{d}$ convergence in distribution). That is not enough: the minimizer is a property of the whole curve, so the two curves must be close everywhere at once.

::: {.intuition}
Why is pointwise convergence not enough? Imagine a narrow, deep pit in the sample terrain whose location moves with the sample size. Fix any $\theta$, and sooner or later the pit moves away from it, so $Q_n(\theta)$ converges at every point. Yet the lowest point of the sample terrain stays in the pit, travels with it, and never approaches $\theta_0$. Uniform convergence requires the maximum distance between the sample curve and the population curve to go to zero, which rules such a pit out.
:::

::: {.assumption}
**M.2 (conditions of Wooldridge 2010 Theorem 12.1, p. 403)** $\{w_i\}$ is a random sample from the population, and (a) $\Theta$ is compact; (b) for each $\theta \in \Theta$, $q(\cdot,\theta)$ is Borel measurable; (c) for each $w$, $q(w,\cdot)$ is continuous on $\Theta$; (d) there is a nonnegative function $b(w)$ such that for all $\theta \in \Theta$, $\lvert q(w,\theta) \rvert \le b(w)$, and $\mathbb{E}[b(w)] < \infty$.
:::

\(a) is almost never a restriction in applications, since $\Theta$ can be taken to be a closed and bounded set large enough to be certain to contain $\theta_0$; (b) is a technical condition that every objective function actually used in econometrics satisfies; (c) rules out discontinuous objectives such as maximum score estimation; (d) is a dominance condition, requiring $q$ to be held down everywhere by a single function with finite mean, which is stronger than "finite expectation at each $\theta$." Upgrading pointwise convergence to uniform convergence takes all four together: compactness lets $\Theta$ be covered by finitely many small neighborhoods, continuity keeps $q$ from varying much within each of them, dominance keeps the expected variation small as well, and measurability guarantees that these expectations are defined, so pointwise convergence at finitely many points is enough to control the whole curve.

::: {.theorem}
**Theorem 1 (uniform weak law of large numbers; Wooldridge 2010 Theorem 12.1)** Under M.2, $\max_{\theta \in \Theta} \lvert Q_n(\theta) - Q_0(\theta) \rvert \xrightarrow{p} 0$.

**Theorem 2 (consistency of the M-estimator; Wooldridge 2010 Theorem 12.2, p. 404)** Under M.1 and M.2, $Q_n$ has a minimizer $\hat\theta$, and $\hat\theta \xrightarrow{p} \theta_0$.
:::

The proof sketch has three steps; the complete proof is in Newey and McFadden (1994), to which Wooldridge (2010) §12.2 also refers. Step one (M.2): Theorem 1 gives uniform convergence; under the same conditions $Q_0$ is continuous, and a continuous function on a compact set attains its minimum, so both $\hat\theta$ and $\theta_0$ exist. Step two (M.1 and M.2): take any open neighborhood of $\theta_0$ and call it $N$. The set $\Theta$ minus $N$ is still compact, and the minimum of $Q_0$ on it is strictly greater than $Q_0(\theta_0)$, so the gap between the two is a positive number. Step three (steps one and two): with probability approaching one, $Q_n$ and $Q_0$ differ everywhere by less than half this gap. On that event, outside $N$, every $\theta$ has $Q_n(\theta) > Q_n(\theta_0) \ge Q_n(\hat\theta)$, so $\hat\theta$ can only lie in $N$. Since $N$ can be arbitrarily small, consistency follows.

In words: consistency needs only two things, a unique solution to the population problem (M.1) and a sample objective that is uniformly close to the population objective (Theorem 1). Under M.1 and M.2, the M-estimator estimates the unique minimizer of the population loss. When the model is misspecified, it converges to the pseudo-true value; that M2m's NLS converges to 0.415 rather than to the local slope 0.55 rests on this same theorem. Who sets the weights: $Q_0(\theta) = \mathbb{E}[q(w,\theta)]$ is a weighted average. The researcher's choice of $q$ fixes how much each observation is penalized at each $\theta$, the population distribution of $w$ fixes the share of each kind of observation, and the two jointly determine the valley floor; neither the sample size nor the estimator's variance plays any part. Moreover, only the part of the distribution of $y$ that $q$ uses can move the valley floor. The Poisson loss $\exp(x'\theta) - y\,x'\theta$ is linear in $y$, so after taking expectations given $x$ only $\mathbb{E}[y \mid x]$ remains, and the conditional variance has no way in. That is why M4's point estimate aims at 0.45 under every degree of overdispersion.

### 3.4 Asymptotic Normality and the Sandwich

Consistency says only that $\hat\theta$ gets close to $\theta_0$, not how it fluctuates around $\theta_0$. Stay with M4: across the 2,000 replications, the standard deviation of the rating-coefficient estimates is 0.0192, the true size of the sampling fluctuation. The Poisson model's built-in standard error averages only 0.0106, while the sandwich standard error averages 0.0187. The former understates the fluctuation by nearly half; the latter is almost exactly right. The derivation below shows where the difference comes from.

Start with the intuition. In each draw, $\hat\theta$ lands near $\theta_0$ but not exactly on it, the outcome of a tug-of-war between two forces. One force is the jitter of the sample gradient at the truth: in the population the first-order condition holds at $\theta_0$, but in a finite sample it carries noise, and the larger this noise, the farther $\hat\theta$ is pushed. The other force is the curvature of the objective function: the steeper the valley, the shorter the distance that the same gradient noise can push the minimizer. The shape of the variance is thus already in view: noise in the middle, and the reciprocal of the curvature on both sides.

The engine of the derivation is a Taylor expansion of the first-order condition. First define the score $s(w,\theta) = \nabla_\theta q(w,\theta)'$, the transpose of the gradient of $q$ with respect to $\theta$; it is a $k \times 1$ column vector. The sample gradient itself is a row vector, and its transpose $\nabla Q_n(\theta)' = n^{-1}\sum_i s(w_i,\theta)$ is exactly the sample average of the score. Since $\hat\theta$ minimizes the sample objective, the sample gradient is zero there. Expand it to first order around $\theta_0$:

$$0 = \nabla Q_n(\hat\theta)' \approx \nabla Q_n(\theta_0)' + \nabla^2 Q_n(\theta_0)\,(\hat\theta - \theta_0),$$

and rearranging and multiplying by $\sqrt{n}$ gives

$$\sqrt{n}(\hat\theta - \theta_0) \approx -\big[\nabla^2 Q_n(\theta_0)\big]^{-1} \sqrt{n}\, \nabla Q_n(\theta_0)'.$$

Each of the two pieces on the right has a destination: the sample Hessian converges to $A_0$, and $\sqrt n$ times the sample gradient is a sum of independent mean-zero terms, which converges to a normal distribution by the central limit theorem. The argument needs the following conditions:

::: {.assumption}
**M.3 (conditions of Wooldridge 2010 Theorem 12.3, p. 407)** In addition to M.1 and M.2: (a) $\theta_0$ is in the interior of $\Theta$; (b) for each $w$, $s(w,\cdot)$ is continuously differentiable on the interior of $\Theta$; (c) every element of $\nabla^2_\theta q(w,\theta)$ is bounded in absolute value by some $b(w)$ with $\mathbb{E}[b(w)] < \infty$; (d) $A_0 = \mathbb{E}[\nabla^2_\theta q(w,\theta_0)]$ is positive definite; (e) $\mathbb{E}[s(w,\theta_0)] = 0$; (f) every element of $s(w,\theta_0)$ has a finite second moment.
:::

Of the six conditions in M.3, three are the ones that cause real trouble in applications. (a) rules out a truth on the boundary of the parameter space. (b) and (d) require a differentiable score and a positive definite expected Hessian: LAD, discussed below, satisfies neither. Under weak identification, (d) nearly fails and the Hessian is nearly singular; most of the numerical difficulties of Section 4.1 come from this. (c), (e), and (f) are technical conditions, and (e) can also be derived from the fact that $\theta_0$ is an interior minimizer, as the last sentence of the proof below shows.

::: {.theorem}
**Theorem 3 (asymptotic normality of the M-estimator; Wooldridge 2010 Theorem 12.3, equations (12.18) to (12.21))** Under M.1 through M.3,

$$\sqrt{n}(\hat\theta - \theta_0) \xrightarrow{d} \mathcal{N}\big(0,\ A_0^{-1} B_0 A_0^{-1}\big),$$

where $A_0 = \mathbb{E}[\nabla^2_\theta q(w,\theta_0)]$ and $B_0 = \mathbb{E}[s(w,\theta_0)\,s(w,\theta_0)']$. Moreover (as the proof shows), $\sqrt n(\hat\theta - \theta_0) = -A_0^{-1}\, n^{-1/2}\sum_{i=1}^n s(w_i,\theta_0) + o_p(1)$ (equation (12.16)).
:::

The proof sketch has four steps. Step one (M.3(a) and (b), Theorem 2): with probability approaching one $\hat\theta$ lies in the interior, so the first-order condition $\sum_i s(w_i,\hat\theta) = 0$ holds (equation (12.14)). Step two (M.3(b) and (c)): take a mean-value expansion of each of its components around $\theta_0$ (equation (12.15)). The mean values lie between $\hat\theta$ and $\theta_0$ and so also converge to $\theta_0$, and the dominance condition makes the average sample Hessian converge to $A_0$ (Wooldridge 2010 Lemma 12.1). Step three (M.3(d)): $A_0$ is invertible, so the expansion can be solved for the linear representation above. Step four (M.3(e) and (f)): the $s(w_i,\theta_0)$ are independent and identically distributed with mean zero and finite variance, so the central limit theorem gives $\mathcal{N}(0,B_0)$; multiplying by $-A_0^{-1}$ completes the proof. Condition (e) usually need not be assumed separately: $\theta_0$ is an interior minimizer of $Q_0$, and when derivative and expectation can be interchanged, $\mathbb{E}[s(w,\theta_0)] = \nabla_\theta Q_0(\theta_0)' = 0$ (equations (12.22) and (12.23)).

In words: under M.1 through M.3, the M-estimator estimates the unique minimizer $\theta_0$ of the population loss, and its fluctuation is described by the sandwich. Nowhere does the theorem use the assumption that the model is correct: under misspecification, as long as the pseudo-true value $\theta_*$ satisfies the same conditions, the conclusion holds unchanged, with $A$ and $B$ evaluated at $\theta_*$ instead. The weights are set by $q$: the loss function decides which deviations are expensive and which are cheap, and the distribution of $w$ then averages them; together they determine both where the valley floor is and what $A_0$ and $B_0$ are. Back to M4: the Poisson model's built-in variance is $A_0^{-1}$, which amounts to assuming $B_0 = A_0$. Overdispersion, however, makes $B_0$ much larger than $A_0$, and the nominal 95% interval covers the truth in only 72.9% of the replications, against a coverage rate of 94.4% for the sandwich.

The linear representation of equation (12.16) recurs throughout what follows: $-A_0^{-1}s(w_i,\theta_0)$ is the influence of observation $i$ on the estimator, and the estimator's sampling error is approximately the average of these influences. The two-step estimator of Section 3.6 adds a term to this representation, and the Delta method and the cluster variance of Chapter 8 also start from it.

Right after the theorem, Wooldridge points out an important example that the theorem does not cover: LAD. Its loss $\lvert y - m(x,\theta) \rvert$ is not differentiable at zero, and the Hessian of linear median regression is the zero matrix under any reasonable definition, violating condition (d). LAD is still $\sqrt n$-asymptotically normal, but this cannot be proved with Theorem 3; Chapter 13 takes another route when it covers quantile regression.

To estimate $A_0$ and $B_0$, use their sample counterparts: $\hat B = n^{-1}\sum_i \hat s_i \hat s_i'$ (equation (12.45)); $\hat A$ comes in two versions (equations (12.42) and (12.44)), one the average of the sample Hessians, the other an average taken after first replacing each Hessian by its conditional expectation over $y$ given $x$. Combining them gives $\widehat{\mathrm{Avar}}(\hat\theta) = \hat A^{-1}\hat B\hat A^{-1}/n$ (equation (12.47)), the sandwich in the sense of Huber (1967) and White (1980): the two slices of bread are the reciprocal of the curvature, and the filling is the noise of the gradient. Wooldridge calls the version with the sample Hessian fully robust and the version with the conditional expectation semirobust. The latter requires the relevant features of the conditional distribution to be correctly specified; otherwise the conditional expectation itself is computed wrongly. The sandwich is also unaffected by the scale of $q$: multiply $q$ by 2 and $A_0$ doubles while $B_0$ quadruples, leaving $A_0^{-1}B_0A_0^{-1}$ unchanged.

When does the sandwich collapse? When $B_0$ is proportional to $A_0$:

$$B_0 = \sigma^2 A_0 \quad \Longrightarrow \quad A_0^{-1}B_0A_0^{-1} = \sigma^2 A_0^{-1}.$$

Here $\sigma^2$ is the scale in the conditional variance $\mathrm{Var}(y \mid x) = \sigma^2 h(x)$: for NLS under homoskedasticity $h(x) \equiv 1$, and for Poisson-type models $h(x) = \mathbb{E}[y \mid x]$. Wooldridge calls the condition on the left the generalized information matrix equality (equation (12.53)); the collapsed variance is in Lemma 12.2. NLS satisfies it under homoskedasticity (Section 3.5, writing the loss as $q/2$), Poisson PML satisfies it when the variance is proportional to the mean (Section 3.7), and a correctly specified likelihood satisfies the version with $\sigma^2 = 1$, which is the information matrix equality of Chapter 6. The collapsed formula is right only when this equality holds, while the sandwich is right whether it holds or not; the sandwich is the default choice throughout Part II.

### 3.5 Nonlinear Least Squares

On M2, the Gauss-Newton iteration (Section 4.1), started at all zeros, converges in 11 iterations to $(0.117, 0.454, 0.272)$, within $1.2 \times 10^{-8}$ of the result of `stats::nls`. The point estimates are fine; the trouble is in the standard errors. For the listings elasticity, the homoskedastic formula gives a standard error of 0.0197 and the semirobust sandwich (equation (12.52)) gives 0.0331, 1.68 times as large. The reason is that M2's error is multiplicative, so the sales of larger merchants also fluctuate more.

::: {.assumption}
**NLS.1 (Wooldridge 2010 §12.1, p. 399)** For some $\theta_0 \in \Theta$, $\mathbb{E}[y \mid x] = m(x,\theta_0)$.

**NLS.2 (Wooldridge 2010 §12.2, p. 401)** For all $\theta \in \Theta$ with $\theta \neq \theta_0$, $\mathbb{E}\big[(m(x,\theta_0) - m(x,\theta))^2\big] > 0$.

**NLS.3 (Wooldridge 2010 §12.5.1, p. 417)** $\mathrm{Var}(y \mid x) = \sigma^2$.
:::

NLS.1 is the only economic assumption that matters: the functional form of the mean is correctly specified. It allows the error $\varepsilon = y - m(x,\theta_0)$ to be heteroskedastic and arbitrarily distributed, requiring only $\mathbb{E}[\varepsilon \mid x] = 0$. NLS.2 is the form M.1 takes for NLS (see the decomposition of Section 3.2). In the linear model it says that $\mathbb{E}[xx']$ has full rank, and common means such as the exponential and the logistic generally satisfy it as long as $x$ has no perfect collinearity. NLS.3 affects only the variance formula, not consistency, and for nonnegative outcomes it is almost always wrong: the variances of sales and of order counts grow with the mean. Below, $\nabla_\theta m(x,\theta)$ denotes the gradient of the mean with respect to $\theta$, a $1 \times k$ row vector. When computing $A_0$ and $B_0$ for NLS we follow Wooldridge and take $q = (y - m(x,\theta))^2/2$: multiplying by a constant changes neither $\hat\theta$ nor the sandwich (Section 3.4). With the $(y - m(x,\theta))^2$ of the table in Section 3.1, homoskedasticity would give $B_0 = 2\sigma^2 A_0$, and halving the loss is what makes the collapse condition read exactly $B_0 = \sigma^2 A_0$.

::: {.theorem}
**Theorem 4 (NLS; Wooldridge 2010 §12.2, §12.3, and §12.5.1, equations (12.25) to (12.30) and (12.53) to (12.58))** When computing $A_0$ and $B_0$, take $q = (y - m(x,\theta))^2/2$, and assume that $m(x,\cdot)$ is twice continuously differentiable and that the regularity conditions of M.2 and M.3 hold.

\(i) Under NLS.1 and NLS.2, $\hat\theta \xrightarrow{p} \theta_0$ and $\sqrt n(\hat\theta - \theta_0) \xrightarrow{d} \mathcal{N}(0, A_0^{-1}B_0A_0^{-1})$, where $A_0 = \mathbb{E}[\nabla_\theta m(x,\theta_0)'\nabla_\theta m(x,\theta_0)]$ and $B_0 = \mathbb{E}[\varepsilon^2\,\nabla_\theta m(x,\theta_0)'\nabla_\theta m(x,\theta_0)]$;

\(ii) if NLS.3 also holds, $B_0 = \sigma^2 A_0$, and the asymptotic variance collapses to $\sigma^2 A_0^{-1}$.
:::

The proof sketch has three steps. Step one (NLS.1 and NLS.2): the decomposition of Section 3.2 shows that $\theta_0$ is the unique minimizer of the mean squared error, which is M.1, and Theorem 2 gives consistency. Step two (NLS.1): the score is $s = -\nabla_\theta m'(y - m)$, whose conditional expectation given $x$ is zero (equation (12.26)), so M.3(e) holds. The Hessian is $\nabla_\theta m'\nabla_\theta m - \nabla^2_\theta m\,(y - m)$; at $\theta_0$ its second term has conditional expectation zero given $x$, so only the outer product of gradients remains in $A_0$ (equations (12.29) and (12.30)), and Theorem 3 gives (i). Step three (NLS.3): $\mathbb{E}[\varepsilon^2 \mid x] = \sigma^2$ means that in $B_0$ the factor $\varepsilon^2$ can be replaced by a constant (equation (12.56)), which is (ii).

In words: under NLS.1 and NLS.2, NLS estimates the true parameters of the conditional mean, and M2's listings-elasticity estimate of 0.454 aims at the truth, 0.45. NLS.3 decides only whether the variance can collapse, not what is estimated. When the mean is misspecified, NLS no longer estimates any true parameter; it estimates the pseudo-true value of Theorem 5 below. Who sets the weights: the first-order condition $\sum_i \nabla_\theta m_i'(y_i - m_i) = 0$ multiplies each residual by the gradient of the mean. For an exponential mean that gradient is $m\,x'$, so observations with large fitted means weigh heavily in the first-order condition. When the mean is correctly specified, this weighting does not move the probability limit, because $\theta_0$ makes the conditional residual zero at every $x$, and zero stays zero however it is weighted. What the weighting moves is the variance: $A_0$ weights by $\nabla_\theta m'\nabla_\theta m$, $B_0$ multiplies this by $\varepsilon^2$ as well, and when the error variance grows with the mean, large merchants push $B_0$ up faster than $A_0$. M2's multiplicative error is of this kind: NLS.3 fails, and the sandwich standard error of the listings elasticity, 0.0331, is 1.68 times the homoskedastic 0.0197. In the M2 sample, the intercept interval [0.0523, 0.1809] just misses the truth, 0.1823, while the listings and free-shipping intervals both cover theirs; Section 6.6 returns to this point.

What happens when the mean is misspecified? Look first at the pseudo-true values on M2m (figure below). The true log mean is concave, with a local slope $0.55 - 0.2x$ that declines in $x$ and equals 0.55 at $x = 0$. The NLS pseudo-true slope is 0.415; the sample estimate is 0.441 with a fully robust standard error of 0.0163, and the interval [0.409, 0.473] covers 0.415 but not 0.55. On the same data, PPML estimates 0.484, and its interval [0.458, 0.510] covers its own pseudo-true value 0.473 but not NLS's 0.415. Both estimators are "consistent"; they are just consistently estimating different things.

![M2m: the true conditional mean (log scale) and the exponential means at the NLS and PPML pseudo-true values; the points show half of the sample.](assets/fig/fig_05_misspec.svg)

Wooldridge (2010) §12.3 restates the result of White (1981):

::: {.theorem}
**Theorem 5 (NLS with a misspecified mean; White 1981 as restated in Wooldridge 2010 p. 409, variance estimation in §12.5.1)** Suppose $\mathbb{E}[y \mid x]$ is not in $\{m(x,\theta): \theta \in \Theta\}$, but

$$\theta_* = \arg\min_{\theta \in \Theta} \mathbb{E}\big[(\mathbb{E}[y \mid x] - m(x,\theta))^2\big]$$

is unique, and M.2 and M.3 hold at $\theta_*$. Then $\hat\theta \xrightarrow{p} \theta_*$ and $\sqrt n(\hat\theta - \theta_*) \xrightarrow{d} \mathcal{N}(0, A_*^{-1}B_*A_*^{-1})$, where $A_* = \mathbb{E}[\nabla_\theta m'\nabla_\theta m - \nabla^2_\theta m\,(y - m)]$ and $B_* = \mathbb{E}[(y - m)^2\,\nabla_\theta m'\nabla_\theta m]$, both evaluated at $\theta_*$.
:::

The proof sketch has two steps. Step one ($\theta_*$ unique): $\mathbb{E}[(y - m)^2] = \mathbb{E}[(y - \mathbb{E}[y \mid x])^2] + \mathbb{E}[(\mathbb{E}[y \mid x] - m)^2]$, and the first term does not depend on $\theta$, so the population valley floor of the mean squared error is $\theta_*$ and Theorem 2 applies. Step two (M.3 at $\theta_*$): Theorem 3 applies, but $\nabla^2_\theta m$ and $y - m(x,\theta_*)$ are no longer conditionally uncorrelated, so the second term of the Hessian cannot be dropped and $\hat A$ must use the sample Hessian (equations (12.42) and (12.48)). The semirobust formula that uses only the outer product of gradients (equation (12.52)) no longer holds (Wooldridge 2010 p. 416).

In words: the premise of Theorem 5 is precisely a misspecified mean. As long as the pseudo-true value is unique and M.2 and M.3 hold there, NLS estimates $\theta_*$, the best approximation to the conditional mean in the mean-squared-error sense, rather than any structural parameter. Inference is still possible, but the bread must use the sample Hessian, second-derivative term included; the semirobust formula that uses only the outer product of gradients no longer holds (on M2m the two give 0.0163 and 0.0177). Who sets the weights: in the first-order condition, the approximation error is weighted by the distribution of $x$ times the fitted mean itself, so observations with large means dominate. This weight comes from the loss function, and a different loss gives a different pseudo-true value. On M2m the pseudo-true NLS slope is 0.415 and the PPML one is 0.473, a gap of 0.058 that comes entirely from the choice of loss function. The difference between NLS and PPML is therefore itself a diagnostic of the mean specification: when the mean is correctly specified, the two have the same probability limit.

Why do the two pseudo-true values differ? Look at the first-order conditions each one satisfies. The NLS pseudo-true value solves $\mathbb{E}[\nabla_\theta m(x,\theta_*)'(\mathbb{E}[y \mid x] - m(x,\theta_*))] = 0$, and since the gradient of an exponential mean is $m\,x'$, this condition becomes

$$\mathbb{E}\big[m(x,\theta_*)\,x\,\big(\mathbb{E}[y \mid x] - m(x,\theta_*)\big)\big] = 0,$$

whereas, as Section 3.7 will prove, the PPML pseudo-true value satisfies

$$\mathbb{E}\big[x\,\big(\mathbb{E}[y \mid x] - \exp(x'\theta_*)\big)\big] = 0.$$

In the first-order condition, the NLS approximation error carries an extra factor of the fitted mean, so regions with large means get large weights. On M2m the slope is positive, so a large mean means a large $x$, where the true slope is flatter; the NLS slope is therefore smaller. PPML's first-order condition multiplies the error only by $x$, not by the fitted mean, so regions with large means get no extra weight, and its slope is larger than the NLS slope (pseudo-true values 0.4733 for PPML and 0.4154 for NLS).

This is of a piece with the classic result for OLS. Even when the conditional mean is not linear, OLS still gives its best linear approximation in the mean-squared-error sense; NLS swaps "linear" for the exponential or some other shape and gives the best approximation within that family of functions. A regression can almost always tell you what the conditional mean looks like, but what it tells you is "the closest one among the shapes you allowed, measured by your loss." If what comes out is not what you wanted, the problem usually lies in the shape and the loss chosen by whoever asked the question, not in the regression itself.

### 3.6 Two-Step M-Estimation and Generated Regressors

M3 targets the elasticity of monthly sales with respect to shop traffic. The traffic index $t$ is log visits measured as a deviation from the median, and the true conditional mean is $\mathbb{E}[y \mid z] = \exp(\beta_0 + \beta_1 t)$ with $\beta_1 = 0.8$, where $t = 0.30\,z_1 + 0.20\,z_2$, $z_1$ is a search-position score, and $z_2$ is log followers. The platform cannot see $t$; it sees only one week of log visits, $v = t + 0.60\,e$, which is very noisy. With $v$ used directly as the regressor, NLS gives 0.209, only 0.261 of the truth: attenuation caused by measurement error. The two-step approach first regresses $v$ on $(1, z_1, z_2)$ by OLS ($R^2$ of 0.253), calling the coefficients $\hat\gamma$, and then uses the fitted value $\hat t = z'\hat\gamma$ in place of $t$. A regressor estimated in an earlier step in this way is called a generated regressor. The second-step NLS gives 0.806, almost exactly on the truth. The problem is in the standard errors: treating $\hat t$ as data gives a standard error of 0.0302; accounting for the first step's estimation error gives 0.0431, 1.425 times as large. Across 2,000 Monte Carlo replications, the true standard deviation of the estimate is 0.0432, the naive standard error averages 0.0298, and the nominal 95% interval covers the truth in only 82.35% of the replications; the corrected standard error averages 0.0423, with coverage of 94.9%.

Why is the naive standard error too small? It treats $\hat t$ as a number as certain as $t$, which amounts to pretending that the platform knows every merchant's true traffic. But $\hat\gamma$ comes from a single regression on 2,000 merchants, and its error points the same way for all of them: if $\hat\gamma$ is a little too high, everyone's $\hat t$ is inflated together, and the second-step slope shifts systematically with it. An error shared by all observations does not cancel in the second step's average, yet the naive formula treats it as if it did not exist. Using $v$ directly is a different matter: each merchant's own measurement noise mixes variation unrelated to sales into the regressor and pushes the slope toward zero, which is a bias problem, not a variance problem.

Wooldridge (2010) §12.4 casts the generated regressor, and problems like it, as two-step M-estimation: the first step yields $\hat\gamma$, and the second step solves $\min_\theta n^{-1}\sum_i q(w_i,\theta,\hat\gamma)$. The first-step parameter $\gamma$ enters only as an extra argument, and the score and Hessian are derivatives with respect to $\theta$ only, written $s(w,\theta,\gamma)$. In M3, $\theta = (\beta_0, \beta_1)$, $\gamma$ is the vector of first-step OLS coefficients, and $q$ is the squared-error loss of the second step. The $\psi_i$ in TS.1 below is the influence function of the first step: the first step's estimation error can be written as an average of the contributions of individual observations.

::: {.assumption}
**TS.1 (first step; Wooldridge 2010 equation (12.38))** $\hat\gamma \xrightarrow{p} \gamma_*$, and $\sqrt n(\hat\gamma - \gamma_*) = n^{-1/2}\sum_{i=1}^n \psi_i(\gamma_*) + o_p(1)$, with $\mathbb{E}[\psi_i(\gamma_*)] = 0$ and finite variance.

**TS.2 (identification given $\gamma_*$; Wooldridge 2010 §12.4.1)** For all $\theta \neq \theta_0$, $\mathbb{E}[q(w,\theta_0,\gamma_*)] < \mathbb{E}[q(w,\theta,\gamma_*)]$; M.2 holds for $q(w,\theta,\gamma)$ on the product of the parameter spaces of $\theta$ and $\gamma$.

**TS.3 (smoothness; Wooldridge 2010 §12.4.2)** M.3 holds at $(\theta_0,\gamma_*)$; $s(w,\theta,\gamma)$ is continuously differentiable in $\gamma$, with its derivative dominated by a function with finite mean.
:::

For OLS, for example, $\psi_i = (\mathbb{E}[zz'])^{-1} z_i (v_i - z_i'\gamma_*)$; the MLE and GMM of Chapters 6 and 7 also have such representations. TS.1 writes $\gamma_*$ rather than $\gamma_0$ because the first step may converge to a probability limit with no economic meaning, such as the coefficients of a misspecified variance model. When the first step is correctly specified, $\gamma_* = \gamma_0$, the first step's true parameter. TS.2 and TS.3 merely carry M.1 through M.3 over to $(\theta,\gamma)$; the extra requirement is that the second-step score be smooth in $\gamma$, so that $\Gamma_0$ below is defined.

::: {.theorem}
**Theorem 6 (the two-step M-estimator; Wooldridge 2010 §12.4.2, equations (12.33) to (12.41))** Under TS.1 through TS.3, $\hat\theta \xrightarrow{p} \theta_0$, and $\sqrt n(\hat\theta - \theta_0)$ is asymptotically normal with mean zero and variance

$$A_0^{-1}\,\mathrm{Var}\big(s_i(\theta_0,\gamma_*) + \Gamma_0\,\psi_i(\gamma_*)\big)\,A_0^{-1},$$

where $\Gamma_0 = \mathbb{E}[\nabla_\gamma s(w,\theta_0,\gamma_*)]$.

\(i) If $\Gamma_0 = 0$ (equation (12.37)), the middle term is just $B_0$, the first step does not affect the asymptotic variance, and the sandwich computed by treating $\hat\gamma$ as known is correct; (ii) if $\Gamma_0 \neq 0$ and $\mathbb{E}[s_i\psi_i'] = 0$, the middle term equals $B_0 + \Gamma_0\,\mathbb{E}[\psi_i\psi_i']\,\Gamma_0'$, and ignoring the first step overstates the precision (Wooldridge 2010 p. 413).
:::

The proof sketch has three steps. Step one (TS.2 and TS.3): expanding as for Theorem 3 gives $\sqrt n(\hat\theta - \theta_0) = -A_0^{-1}n^{-1/2}\sum_i s_i(\theta_0,\hat\gamma) + o_p(1)$ (equation (12.33)). Step two (TS.1 and TS.3): expand $s_i(\theta_0,\hat\gamma)$ around $\gamma_*$ by the mean-value theorem in $\gamma$, which gives $n^{-1/2}\sum_i s_i(\theta_0,\hat\gamma) = n^{-1/2}\sum_i s_i(\theta_0,\gamma_*) + \Gamma_0\sqrt n(\hat\gamma - \gamma_*) + o_p(1)$ (equation (12.35)). Step three (TS.1): substituting the first step's influence-function representation turns the sum into $n^{-1/2}\sum_i(s_i + \Gamma_0\psi_i)$, which has mean zero, and the central limit theorem gives the middle term (equations (12.39) to (12.41)). When $\Gamma_0 = 0$ the second term vanishes, which is (i); when the scores of the two steps are uncorrelated the variances add, which is (ii).

In estimation, $\Gamma_0$ is replaced by $\hat\Gamma = n^{-1}\sum_i \nabla_\gamma s_i(\hat\theta,\hat\gamma)$ (equation (12.61)), $\psi_i$ by its sample version, and the filling in the middle becomes $n^{-1}\sum_i (\hat s_i + \hat\Gamma\hat\psi_i)(\hat s_i + \hat\Gamma\hat\psi_i)'$ (equation (12.62)); this is the algorithm of `two_step_vcov` in Section 6.4. Other routes are to write the two steps as one stacked system of moment conditions (Newey 1984) or to bootstrap the entire two-step procedure; the formula of Murphy and Topel (1985) and bootstrap practice are left to Section 4.4 of Chapter 8.

Why is M3's $\Gamma_0$ not zero? The second-step score is $s = -m\,(1, t)'\,(y - m)$ with $m = \exp(\beta_0 + \beta_1 z'\gamma)$, and $\gamma$ enters in two places, through $m$ and through the regressor. Differentiate with respect to $\gamma$ and take expectations: the terms containing $y - m$ have expectation zero, leaving $\Gamma_0 = \beta_1\,\mathbb{E}[m^2 (1, t)'z']$. This contains the factor $\beta_1$: when $\beta_1 = 0$ the regressor plays no role, and the first step's error has no channel to pass through. The Monte Carlo grid shows this pattern (the figure in Section 6.5): at $\beta_1 = 0$ the naive and corrected coverage rates are both 0.952; at $\beta_1 = 0.4$ the naive rate falls to 0.895, and at $\beta_1 = 1.2$ to 0.804, while the corrected rate stays between 0.947 and 0.963 throughout.

In words: under TS.1 through TS.3, the two-step estimator estimates the minimizer of the second-step loss at $\gamma_*$. M3's first step is consistent for $\gamma$ ($\mathbb{E}[v \mid z] = t$), so 0.806 aims at the elasticity 0.8 itself. This also relies on M3's traffic index being exactly a linear function of $z$. If $t = z'\gamma + r$ with $r$ unobserved, then even with a correctly specified first step, the conditional expectation of $\exp(\beta_1 t)$ given $z$ is not $\exp(\beta_1 z'\gamma)$, and plugging fitted values into a nonlinear second step is generally inconsistent (Jensen's inequality again; only in special cases such as $r$ independent of $z$ does the bias of an exponential mean fall entirely into the intercept). Such cases call for approaches like the control function of Chapter 12. Here the generated regressor is a pure inference problem, whereas using the noisy $v$ directly, which gives 0.209, is a bias problem caused by measurement error. If the first step is misspecified, $\gamma_* \neq \gamma_0$, the second step aims at something else, and $\Gamma_0$ is generally not zero either. Who sets the weights: each first-step error term $\psi_i$ enters the second-step variance only after being multiplied by $\Gamma_0$; $\Gamma_0$ is determined by the second-step model, and $\psi_i$ by the first-step design. In M3, $\Gamma_0 = \beta_1\,\mathbb{E}[m^2 (1, t)'z']$ contains the factor $\beta_1$ and is also weighted by the square of the fitted mean, so merchants with large sales amplify the first-step error the most. The size of $\psi_i$ depends on the noise in visits and on the explanatory power of the first step, and with an $R^2$ of only 0.253 the first step is not precise. M3 is exactly the case of Theorem 6 (ii): the multiplicative error of the second step is independent of the visit noise $e$ of the first step, so the scores of the two steps are uncorrelated and the corrected asymptotic variance cannot be smaller than the naive one. That 0.0431 exceeds 0.0302 is therefore the direction the theory requires here, not an accident. A useful special case comes from Wooldridge (2010) p. 412: for weighted NLS with a correctly specified mean, $\Gamma_0 = 0$ whether or not the weight function is right, so estimated weights require no correction of the variance.

### 3.7 Poisson PML

Return to Section 1. It is no coincidence that PPML on M1 gives exactly the same coefficient with dollar and with cent accounting. In the M4 sample it estimates the rating coefficient at 0.493, with a Poisson model standard error of 0.0104 and a sandwich standard error of 0.0206, 1.97 times as large. The Pearson dispersion is 2.49, so the data are clearly overdispersed, yet the point estimate is still usable.

::: {.assumption}
**PML.1 (conditional mean; Wooldridge 2010 equation (18.7), Gourieroux, Monfort and Trognon 1984b)** For some $\theta_0$, $\mathbb{E}[y \mid x] = \exp(x'\theta_0)$.

**PML.2 (identification; Wooldridge 2010 p. 728)** $\mathbb{E}[xx']$ is positive definite, that is, $x$ has no perfect collinearity.

**PML.3 (variance proportional to the mean; Wooldridge 2010 equation (18.3))** $\mathrm{Var}(y \mid x) = \sigma^2\,\mathbb{E}[y \mid x]$, $\sigma^2 > 0$.
:::

PML.1 restricts only the mean: $y$ can be a count or a continuous nonnegative variable such as sales, it can have many zeros, and its variance can have any shape. PML.3 is the variance assumption commonly used in the GLM (generalized linear model) literature and is needed only for the collapsed formula; $\sigma^2 > 1$ means overdispersion, and the Poisson distribution itself has $\sigma^2 = 1$. M4's NB2 variance $\mathbb{E}[y \mid x] + \eta^2\,\mathbb{E}[y \mid x]^2$ does not satisfy PML.3.

::: {.theorem}
**Theorem 7 (Poisson PML; Wooldridge 2010 §18.2.2 and §18.2.3, derived from Theorems 12.2 and 12.3; general theory in Gourieroux, Monfort and Trognon 1984b, the Poisson case in 1984a)** Take $q(w,\theta) = \exp(x'\theta) - y\,x'\theta$, and assume the regularity conditions of M.2 and M.3 hold.

\(i) Under PML.1 and PML.2, $\hat\theta \xrightarrow{p} \theta_0$, with no assumption at all about $\mathrm{Var}(y \mid x)$, and $y$ need not be a count;

\(ii) $\sqrt n(\hat\theta - \theta_0) \xrightarrow{d} \mathcal{N}(0, A_0^{-1}B_0A_0^{-1})$, with $A_0 = \mathbb{E}[\exp(x'\theta_0)\,xx']$ and $B_0 = \mathbb{E}[\mathrm{Var}(y \mid x)\,xx']$;

\(iii) if PML.3 also holds, $B_0 = \sigma^2 A_0$, and the asymptotic variance collapses to $\sigma^2 A_0^{-1}$; with $\sigma^2 = 1$ this is the Poisson model's built-in variance.
:::

The proof sketch has three steps. Step one (PML.1): given $x$, $\mathbb{E}[q(w,\theta) \mid x] = \exp(x'\theta) - \mathbb{E}[y \mid x]\,x'\theta$ is strictly convex as a function of $x'\theta$, and its derivative $\exp(x'\theta) - \mathbb{E}[y \mid x]$ is zero at $x'\theta = \log \mathbb{E}[y \mid x] = x'\theta_0$, so $\theta_0$ minimizes the conditional expected loss at every $x$ and hence minimizes $\mathbb{E}[q]$. This step uses only the mean. Step two (PML.2): for $\theta \neq \theta_0$, $x'\theta \neq x'\theta_0$ happens with positive probability, and strict convexity makes $\mathbb{E}[q]$ strictly larger. This is M.1, and Theorem 2 gives (i). Step three (PML.1 and PML.3): the score is $s = (\exp(x'\theta) - y)\,x$, and the Hessian $\exp(x'\theta)\,xx'$ does not contain $y$, so $A_0$ does not depend on the variance. Moreover, $B_0 = \mathbb{E}\big[\mathbb{E}[(y - \exp(x'\theta_0))^2 \mid x]\,xx'\big] = \mathbb{E}[\mathrm{Var}(y \mid x)\,xx']$; Theorem 3 gives (ii), and substituting PML.3 gives (iii). The fully robust variance estimator is equation (18.14), and the estimator of $\sigma^2$ is equation (18.15). When the mean is misspecified, PML.1 no longer holds. Provided the pseudo-true value $\theta_*$ is the unique interior minimizer of $\mathbb{E}[q]$, and M.2 and M.3 hold at that point, the conclusion of Section 3.4 holds at $\theta_*$ instead: the expected score is zero there (M.3(e)), that is, $\mathbb{E}[(\exp(x'\theta_*) - y)\,x] = 0$, and the law of iterated expectations turns this into $\mathbb{E}[x(\mathbb{E}[y \mid x] - \exp(x'\theta_*))] = 0$, which is the first-order condition cited in Section 3.5.

On M4, $B_0$ can be written out explicitly, which shows where the model variance goes wrong (in M4, $\theta_0 = (1.4, 0.45)$). Under NB2, $\mathrm{Var}(y \mid x) = \exp(x'\theta_0) + \eta^2\exp(2x'\theta_0)$, so

$$B_0 = A_0 + \eta^2\,\mathbb{E}\big[\exp(2x'\theta_0)\,xx'\big],$$

the second term is positive semidefinite, so the sandwich can never be smaller than the model variance. This term is also weighted by the square of the mean: merchants with more orders weigh more, and it is precisely their fluctuations that the model variance understates. As $\eta^2 \to 0$ the second term vanishes and we are back to Poisson.

Two corollaries matter most for Section 1. First, unit invariance: replace $y$ by $a\,y$ ($a > 0$, with $a = 100$ when switching from dollars to cents); then $\mathbb{E}[a\,y \mid x] = \exp(\log a + x'\theta_0)$, and only the intercept shifts, by $\log a$. The sample first-order condition $\sum_i x_i(y_i - \exp(x_i'\hat\theta)) = 0$ likewise moves only the intercept, so M1's dollar and cent estimates differ by $4.4 \times 10^{-16}$, mere floating-point error. Second, the coefficient of a binary treatment is a log ratio: PPML on $D$ alone is a saturated model, and the coefficient of $D$ is exactly the log of the ratio of the two groups' sample means: 0.305 on M1, which is $\log(237.82/175.25)$. When $D$ is randomly assigned, it estimates $\log(\mathbb{E}[y_i(1)]/\mathbb{E}[y_i(0)])$. With category and tenure added, PML.1 still holds because the conditional mean of Section 2.1 is exponential, and the estimate 0.302 aims at the same 0.323.

The coefficients are also easy to read (Wooldridge 2010 §18.2.1). The coefficient of a continuous variable is a semi-elasticity of the conditional mean: M4's 0.45 is the semi-elasticity of the average order count with respect to the rating, so a rating one standard deviation higher multiplies the average order count by $e^{0.45}$. A logged regressor gives an elasticity, and for a dummy variable, the exponential of the coefficient minus one is the proportional change: M1's 0.302 corresponds to 35.2%. There is also a shortcut to the average partial effect in levels: with an intercept, the first-order condition makes the residuals sum to zero, so the sample average partial effect of a continuous variable is exactly $\bar y$ times its coefficient.

In words: under PML.1 and PML.2, PPML estimates the parameters of an exponential mean, which for a randomly assigned binary treatment is the log ratio of average outcomes, with zeros counted in dollars (on M1 the truth is 0.323 and the estimate 0.302). No assumption about the variance or the shape of the distribution is needed. When the mean is misspecified, PPML converges to the pseudo-true value satisfying $\mathbb{E}[x(\mathbb{E}[y \mid x] - \exp(x'\theta_*))] = 0$, which is 0.473 on M2m. Who sets the weights: the first-order condition counts residuals in levels without multiplying them by the fitted mean, so PPML treats every dollar alike, and going from zero to USD 300 counts as much as going from USD 300 to USD 600. In the variance, $B_0$ weights by $\mathrm{Var}(y \mid x)$, and under NB2 the extra term weights by the square of the mean, so merchants with many orders weigh the most. The sandwich collapses only under PML.3. When the data really are Poisson, the coverage rates of the model, sandwich, and quasi-Poisson (Poisson variance times the estimated dispersion) standard errors are 0.945, 0.949, and 0.944, and the three are indistinguishable. When the ratio of the average conditional variance to the average mean is 2.83, they become 0.729, 0.944, and 0.913. Quasi-Poisson is right only under PML.3, and the NB2 variance is not proportional to the mean, so its correction is not enough either.

### 3.8 log(1+y) and Chen and Roth

The other half of Section 1's problem is log(1+y): the same regression gives 0.593 with dollar accounting and 1.058 with cent accounting. First, consider a classic problem that has nothing to do with zeros.

Even without zeros, taking logs has a cost. Suppose $y > 0$. An OLS regression of $\log y$ estimates $\mathbb{E}[\log y \mid x]$, and

$$\mathbb{E}[\log y \mid x] = \log\mathbb{E}[y \mid x] + \mathbb{E}\big[\log\big(y/\mathbb{E}[y \mid x]\big) \mid x\big].$$

By Jensen's inequality, the last term is at most zero. As long as it does not vary with $x$, it is a constant that enters only the intercept; the most common sufficient condition for this is that the distribution of the multiplicative error $y/\mathbb{E}[y \mid x]$ does not depend on $x$. When the variance of the multiplicative error varies with $x$, this term generally varies with $x$ as well, and so it leaks into the slope. M2's multiplicative error happens to be independent of $x$, so OLS on $\log y$ also recovers the elasticity, a piece of luck that is not common in real data. This is the argument of Santos Silva and Tenreyro (2006), discussed in Section 5.1.

The log(1+y) approach takes a different road: first transform the outcome, then estimate the average effect on the transformed outcome. Chen and Roth (2024) prove that when there is an extensive margin, the effect of such a transformation can be pushed around at will by the choice of unit. The proposition numbers below follow the published version. Let $\lambda(\cdot)$ be a log-like transformation, and let $\mathrm{ATE}_\lambda(a) = \mathbb{E}[\lambda(a\,y_i(1)) - \lambda(a\,y_i(0))]$ be the ATE (average treatment effect) after the outcome is multiplied by $a$.

::: {.intuition}
Two merchants are enough to see how the unit turns into a weight. Because of Spotlight, merchant A goes from zero to USD 300, and merchant B from USD 300 to USD 600. Under the proportional effect on the average, each contributes USD 300 to $\mathbb{E}[y_i(1)] - \mathbb{E}[y_i(0)]$. Under log(1+y), B contributes about $\log 2$, almost the same whether sales are recorded in dollars or in cents; A contributes $\log 301$, which in cents becomes $\log 30001$, an extra $\log 100$ or so out of nowhere. The more merchants like A who go from zero to positive sales, the more the unit affects the result. The log(1+y) regression has not miscalculated anything; it has merely decided for you, quietly, how many "doublings of sales" one "first sale" is worth, and that exchange rate is set by the unit you choose.
:::

::: {.assumption}
**CR.1 (condition 1 of Chen and Roth 2024 Proposition 1)** $\lambda: [0,\infty) \to \mathbb{R}$ is continuous and weakly increasing.

**CR.2 (condition 2)** As $y \to \infty$, $\lambda(y)/\log y \to 1$; both $\log(1 + y)$ and $\mathrm{asinh}(y)$ (the inverse hyperbolic sine) satisfy it.

**CR.3 (condition 3)** The treatment affects the extensive margin: $\Pr(y_i(1) = 0) \neq \Pr(y_i(0) = 0)$.

**CR.4 (condition 4)** For $d = 0, 1$, $\mathbb{E}\big[\lvert \log y_i(d) \rvert \mid y_i(d) > 0\big] < \infty$.
:::

CR.1 and CR.2 restrict only the transformation itself, and both log(1+y) and asinh satisfy them. CR.4 is a moment condition. Only CR.3 concerns the data, and in a randomized experiment it can be tested directly: M1's extensive-margin effect, 0.1011 (standard error 0.0069), is far from zero.

::: {.theorem}
**Theorem 8 (unit dependence of log-like transformations; Chen and Roth 2024 Proposition 1)** Under CR.1 through CR.4, for any given positive number there exists $a > 0$ such that $\lvert \mathrm{ATE}_\lambda(a) \rvert$ equals it exactly; as $a \to 0$, $\mathrm{ATE}_\lambda(a) \to 0$, and as $a \to \infty$, $\lvert \mathrm{ATE}_\lambda(a) \rvert \to \infty$. Moreover (Proposition 4 in the paper's appendix), for large $a$

$$\begin{aligned}
\mathrm{ATE}_\lambda(a) &= \big[\Pr(y_i(1) > 0) - \Pr(y_i(0) > 0)\big]\log a \\
&\quad + o(\log a).
\end{aligned}$$
:::

The proof sketch has three steps (Appendix A.1 of the paper). Step one (CR.1, CR.2, and CR.4): $\mathrm{ATE}_\lambda(a)$ is a continuous function of $a$ (dominated convergence). Step two (CR.1 and CR.4): as $a \to 0$, $\lambda(a\,y) \to \lambda(0)$, the transformations of the two potential outcomes tend to the same constant, and the effect tends to zero. Step three (CR.2 and CR.3): for large $a$ and positive $y$, $\lambda(a\,y) \approx \log a + \log y$. Splitting by margin, merchants who are positive in both states contribute $\log(y_i(1)/y_i(0))$, which does not depend on $a$. Merchants who go from zero to positive contribute $\lambda(a\,y_i(1)) - \lambda(0) \approx \log a + \log y_i(1)$, which grows without bound in $a$, and merchants who go from positive to zero contribute symmetric negative terms. The net effect is approximately the difference in selling probabilities times $\log a$, so $\lvert \mathrm{ATE}_\lambda(a) \rvert \to \infty$. Continuity and the intermediate value theorem then deliver every positive number.

Larkspur's numbers match this approximation almost exactly. The sample estimate of the extensive-margin effect is 0.1011, which times $\log 100$ gives 0.4654; the difference between the cent and dollar log(1+y) coefficients is 0.4650. The figure in Section 1 sweeps $a$ from $10^{-3}$ to $10^{4}$, and the coefficient rises all the way from 0.041 to 1.52, closely tracking the population ATE.

Can a "fixed" transformation be built? Proposition 2 of Chen and Roth (2024) says no, in the form of a trilemma. When the outcome can be zero, the following three properties cannot all hold: (a) the parameter is an average of some individual-level effect that is weakly increasing in the post-treatment outcome and is not constant when both potential outcomes are positive, that is, it responds to changes at the intensive margin; (b) it does not depend on the unit; (c) it is identified from the marginal distributions of the two potential outcomes. (This qualification comes from the corrigendum of Okat, Chen and Roth (2026). As the published version stated the trilemma, the extensive-margin effect $\mathbb{E}[1\{y_i(1)>0\} - 1\{y_i(0)>0\}]$, which is the LPM (linear probability model) row in the table of Section 4.3, satisfies all three. Since it does not respond at all to changes at the intensive margin, the corrected trilemma excludes it, and the authors state that their practical recommendations are unchanged.) Both log(1+y) and asinh satisfy (a) and (c) and give up (b); PPML's log ratio satisfies (b) and (c) and gives up (a): it is the proportional change in the average, not the average of proportional changes. Keeping only the merchants with positive sales (sellers, from here on) and taking logs is a third option: it estimates the intensive margin and drops the zeros. On Larkspur this gives 0.089 (standard error 0.024), and the interval [0.043, 0.136] covers 0.10. This is a special feature of the DGP: whether a merchant sells is independent of the shock to sales, and Spotlight multiplies the selling probability by the same 1.25 in every category, so the sellers in the two groups have the same distribution of category, tenure, and shocks, and there is no selection bias. In general, the marginal merchants whom Spotlight draws into selling differ from those who were selling anyway, and comparisons among sellers are no longer causal comparisons.

In words: under random assignment, the OLS coefficient of log(1+y) or asinh on the treatment estimates $\mathrm{ATE}_\lambda(a)$, an average effect on a transformed outcome; this step uses only random assignment. CR.1 through CR.4 govern the next step: they guarantee that this quantity changes arbitrarily with the unit of measurement $a$, so it is not a percentage effect. No mean model here has to be specified correctly; what is easy to get wrong is the interpretation, namely reading $\mathrm{ATE}_\lambda(a)$ as a proportional effect. If no merchant switched between zero and positive sales because of the treatment, then as $a \to \infty$ only the merchants positive in both states would remain, contributing $\log(y_i(1)/y_i(0))$, and only in this limit would the influence of the unit vanish. Spotlight, however, moves about 0.099 of the merchants from zero to positive (in the DGP the same random draw decides whether a merchant sells in both states, and no merchant stops selling because of Spotlight, so 0.099 is exactly the share of merchants who go from zero to positive). Who sets the weights: the unit sets the weight of the extensive margin, and each tenfold increase in $a$ adds another $\log 10$ times the extensive-margin effect; PPML's log ratio treats every dollar alike and does not depend on the unit. Under the three units, thousands of dollars, dollars, and cents, the coefficients 0.041, 0.593, and 1.058 each correctly estimate their own population ATE, yet none of them is the proportional effect of 38.1%.

## 4 Estimation: Algorithms, Variances, and the Choice of Estimator

### 4.1 Numerical Optimization

Neither NLS nor PPML has a closed-form solution, so $\hat\theta$ must be found by iteration. The workhorse is Newton-Raphson: approximate the objective function at the current point by a quadratic, jump directly to the lowest point of that quadratic, and repeat. With the value at step $j$ written $\theta^{(j)}$, the update is as follows (Wooldridge 2010 equation (12.81)):

$$\theta^{(j+1)} = \theta^{(j)} - \Big[\sum_{i=1}^n \nabla^2_\theta q(w_i,\theta^{(j)})\Big]^{-1} \sum_{i=1}^n s(w_i,\theta^{(j)}).$$

The gradient tells you which direction to go, and the Hessian tells you how far: small steps where the curvature is high, large steps where it is flat. PPML is the ideal case: the Hessian $\sum_i \exp(x_i'\theta)\,x_ix_i'$ is positive definite everywhere, the objective is globally convex, and Newton-Raphson converges fast; with the step halving described after equation (12.83), it converges from any starting point. The IRLS (iteratively reweighted least squares) that `glm` performs for the Poisson log link is exactly Newton-Raphson.

The NLS Hessian contains the term $\nabla^2_\theta m\,(y - m)$ and need not be positive definite far from the truth. Gauss-Newton (Wooldridge 2010 §12.7.3, equation (12.86)) replaces the Hessian with its conditional expectation given $x$, $\nabla_\theta m'\nabla_\theta m$, and in return every step is a single OLS: regress the current residuals on the current gradient, and the coefficients give the direction to move in:

$$\begin{aligned}
\theta^{(j+1)} = \theta^{(j)} &+ \Big[\sum_{i=1}^n \nabla_\theta m_i'\nabla_\theta m_i\Big]^{-1} \\
&\times \sum_{i=1}^n \nabla_\theta m_i'\big(y_i - m(x_i,\theta^{(j)})\big),
\end{aligned}$$

where the gradient $\nabla_\theta m_i$ is evaluated at $\theta^{(j)}$. Replacing the Hessian does not change the end point, and the reason lies in the update itself: as long as the gradient matrix has full column rank, the update is zero if and only if $\sum_i \nabla_\theta m_i'(y_i - m_i) = 0$, that is, the NLS sample score is zero, so the point where Gauss-Newton stops is a solution of the NLS first-order condition, whatever the Hessian was replaced with. The dropped term has conditional expectation zero at the truth and matters only for speed: when the mean is correctly specified, its sample average is small near the solution, the approximate Hessian is almost the same as the true one, and Gauss-Newton converges nearly as fast as Newton-Raphson; when the residuals are large or the mean is misspecified, it converges somewhat more slowly, but the point where it stops is still a solution of the same first-order condition.

In practice, each step should check that the objective actually falls and halve the step if it does not (the practice described after equation (12.83)). On M2, starting from all zeros, Gauss-Newton converges in 11 iterations and `stats::nls` takes 8; the two agree to within $1.2 \times 10^{-8}$. Stopping rules usually look at the relative change in the parameters or at the size of the gradient; one convenient criterion Wooldridge gives is that $n$ times the uncentered $R^2$ of that OLS regression be small enough. After convergence, also check that the Hessian at the final point is positive definite, and that several starting points lead to the same point.

BHHH (Berndt, Hall, Hall and Hausman 1974) is another variant: it replaces the Hessian with the sum of score outer products $\sum_i s_is_i'$ (equation (12.84)), needs only first derivatives, and is naturally positive semidefinite. Its justification in likelihood problems is the information matrix equality of Chapter 6. One caution applies to all of them: if the objective function is nonconvex, or identification is so weak that the Hessian is nearly singular, any algorithm may get lost in a flat or multi-valley terrain and halt at a false valley floor. When that happens, what you should suspect is usually not the algorithm but the model. NLS should be tried from several starting points.

### 4.2 Variance Estimators and Software

| Variance estimator | Bread $\hat A$ | Filling | When it is valid |
|---|---|---|---|
| Fully robust (equation (12.48)) | Average of the sample Hessians | $n^{-1}\sum_i \hat s_i\hat s_i'$ | M.1 through M.3, including a misspecified mean |
| Semirobust (equation (12.49); for NLS, equation (12.52)) | Expected Hessian given $x$ | Same as above | The mean, or the relevant features of the conditional distribution, correctly specified |
| Collapsed (equations (12.54) and (12.55)) | Either | $\hat\sigma^2\hat A$ | In addition $B_0 = \sigma^2 A_0$, as under NLS.3 or PML.3 |
| Two-step correction (equations (12.61) and (12.62)) | Second-step $\hat A$ | Outer product after adding $\hat\Gamma\hat\psi_i$ | Required when $\Gamma_0 \neq 0$ |

All numbers in this chapter use HC0, that is, $\hat A^{-1}\hat B\hat A^{-1}/n$ without a small-sample factor. In software, `sandwich::sandwich()` gives this quantity for `lm`, `glm`, and `nls` objects, and for `nls` it is the semirobust version of equation (12.52); `glm(family = quasipoisson)` has the same point estimates as Poisson, and its built-in variance is the collapsed version, the Poisson variance times the Pearson dispersion; `fixest::fepois(vcov = "hetero")` includes a small-sample factor, which on M1 multiplies the standard error by 1.00015; `numDeriv::hessian` can be used to check a hand-written Hessian. M2m is a reminder that the fully robust variance is not necessarily larger than the semirobust one: there the former gives 0.0163 and the latter 0.0177.

### 4.3 Outcomes with Zeros: The Estimators at a Glance

| Approach | What it estimates | What a change of unit does | Zeros | Estimate on M1 |
|---|---|---|---|---|
| OLS on log(1+y) | $\mathrm{ATE}_\lambda(a)$, $\lambda(y) = \log(1+y)$ | Changes it arbitrarily with $a$ | Kept, with weight set by $\log a$ | Dollars 0.593, cents 1.058 |
| OLS on asinh | Same, with $\lambda = \mathrm{asinh}$ | Same | Same | Dollars 0.663, cents 1.128 |
| LPM | Difference in selling probabilities | Nothing | Only whether zero or not | 0.1011 |
| log y among sellers | Log difference among sellers | Nothing | Dropped | 0.089 |
| PPML | Log ratio of the averages | Shifts only the intercept | Counted in dollars | 0.302 |

The choice depends first on which number you want, not on which regression runs. If you want a percentage effect, PPML gives the percentage change in the average; if what you really care about is the two margins, estimate them separately, and acknowledge that comparisons among sellers are generally not causal comparisons.

## 5 Anchoring Papers

Each of the three anchoring papers corresponds to one thread of the chapter: the bias of log-linear regression under heteroskedasticity and PPML, the unit dependence of log(1+y), and how widespread these problems are in the empirical literature. In the IS and marketing literatures we have not yet found applications to outcomes with zeros that explicitly use PPML or specifically discuss log(1+y), so the third paper comes from finance, where the problem is closest to the one in this chapter. Each paper is laid out in five parts: paper, method, data, results, and limitations.

### 5.1 Santos Silva and Tenreyro (2006)

::: {.case}
Paper: "The Log of Gravity," Review of Economics and Statistics. The question is how the gravity equation of international trade should be estimated.

Method: the gravity equation is a constant-elasticity model $\mathbb{E}[T_{ij} \mid x_{ij}] = \exp(x_{ij}'\beta)$, where $T_{ij}$ is the trade flow from country $i$ to country $j$. The common practice is to take logs and run OLS. The authors point out that, by Jensen's inequality, the difference between $\mathbb{E}[\log T_{ij} \mid x_{ij}]$ and $\log\mathbb{E}[T_{ij} \mid x_{ij}]$ generally varies with $x_{ij}$, so when the multiplicative error is heteroskedastic, the log-linear OLS slope is generally not the elasticity. Taking logs also forces the analyst to drop the observations with zero trade. They advocate estimating the model in levels by PPML: it is consistent as long as the conditional mean is correctly specified (Theorem 7 of this chapter), and it keeps the zeros as well.

Data: a cross-section of bilateral trade between country pairs. The specification with exporter and importer fixed effects has 18,360 country pairs, of which 9,613 have positive trade.

Results: in this specification (Table 5 of the paper), the log-linear OLS estimate of the distance elasticity is -1.347 (standard error 0.031) and the PPML estimate is -0.750 (0.041); as the authors put it, the former is almost twice the latter. The coefficient on colonial ties is 0.666 (0.070) with OLS, implying that trade almost doubles, but only 0.079 (0.134) with PPML, which is statistically insignificant. Keeping the zeros with $\log(1 + T_{ij})$ does not help: the distance elasticity is -1.332 (0.036), almost the same as with log OLS.

Limitations: PPML's consistency rests only on the mean, but its efficiency rests on the variance being proportional to the mean. The authors' test does not reject proportionality at the 5% level, but other data need not behave the same way. The constant-elasticity form of the mean is itself an assumption; the authors run specification tests of the RESET kind, and applications should do the same.
:::

This paper turned "taking logs" from a harmless habit into an identifying assumption that has to be defended: log-linear regression requires the conditional mean of the log of the multiplicative error not to vary with the regressors (a distribution of the multiplicative error that does not depend on the regressors is the most common guarantee), while PPML requires only a correctly specified mean.

### 5.2 Chen and Roth (2024)

::: {.case}
Paper: "Logs with Zeros? Some Problems and Solutions," Quarterly Journal of Economics. The proposition numbers in this chapter follow the published version; the numerical results cited are taken from its arXiv version v7 (November 2023) and the accepted manuscript the authors posted in June 2024; the two versions report the same numbers.

Method: the theoretical part consists of this chapter's Theorem 8 and the trilemma, the latter in the form corrected by the corrigendum of Okat, Chen and Roth (2026) (Section 3.8). The applied part reanalyzes the randomized experiment of Carranza, Garlick, Orkin and Rankin (2022) in South Africa, in which treated job seekers received skill certificates that they could show to employers.

Data: the job-seeker sample of the original experiment; the outcome is weekly hours worked, recorded as zero for those not employed; standard errors are clustered at the level of the assessment date, as in the original paper.

Results: the original paper's estimate with asinh(weekly hours) is 0.201 (0.052), read as a 20% increase in hours. With only the unit changed, asinh(annual hours) gives 0.417 (0.096) and asinh(full-time equivalents) gives 0.031 (0.012). Average weekly hours in the treatment and control groups are 9.84 and 8.85, a ratio of 1.11. The proportional effect implied by PPML is 0.112 (0.081); with covariates it is 0.150 (0.069), with a smaller standard error. The extensive-margin effect is 0.055, an employment rate 5.5 percentage points higher. The authors also replicate 10 papers in the American Economic Review that use asinh; in 5 of them, multiplying the outcome by 100 changes the estimate by more than 100% of the original estimate.

Limitations: the trilemma means there is no free lunch. PPML's proportional effect is not an average of individual proportional effects, and separating the two margins requires additional assumptions. The argument concerns the population quantity an estimator aims at and does not resolve identification itself: in designs that rely on functional form (parallel trends, for example), switching to a proportional effect often means switching to a different set of identifying assumptions, a point the paper discusses explicitly.
:::

The force of this paper is that it turned a seemingly technical detail into a theorem: in the coefficient of log(1+y), the weight on the extensive margin is set by the unit of measurement, and the unit is something the researcher picks offhand.

### 5.3 Cohn, Liu and Wardlaw (2022)

::: {.case}
Paper: "Count (and Count-Like) Data in Finance," Journal of Financial Economics. It takes up the question of which regression to use for outcomes with zeros and for count outcomes, and it replicates real studies. The numbers below are taken from its accepted manuscript of July 2022.

Method: the paper uses econometric derivations and simulations to compare regression models for count and count-like outcomes, then replicates six papers from top finance journals, choosing one main-table specification from each and estimating it both by OLS on log(1+y) and by Poisson regression. The authors show that the coefficient of log(1+y) has no natural economic meaning and that its bias can be large enough to reverse the sign of the coefficient in expectation; Poisson regression with separable fixed effects requires only a correct conditional mean.

Data: four of the papers study firms' annual patent counts and two study plants' toxic emissions. In the replicated data, 69% of the patent counts and 87.6% of the emissions are zero.

Results: in all six replications the coefficients of the two regressions differ greatly, and in three of them even the signs are opposite. In the five replications with control variables, switching the model moves the coefficient more, in every case, than dropping the most important control variable. The authors also count 44 papers on firms' annual patent counts in the three top finance journals from 2011 to 2020, of which 25 use log(1+y) and 23 use only that approach.

Limitations: the authors picked the replicated specifications according to whether the control variables were easy to obtain, and "dropping the most important control variable" is only a benchmark of magnitude, not a causal decomposition. The conclusions concern two kinds of outcomes in finance data, and carrying them over to other fields requires rechecking the share of zeros and the skewness.
:::

Taken together: Santos Silva and Tenreyro show that log-linear regression under heteroskedasticity does not estimate the elasticity, and zeros merely make the problem more visible; Chen and Roth show that finding a log-like transformation that is defined at zero does not fix the problem, because it turns the unit into the weight on the extensive margin; Cohn et al. show that this is not theoretical nitpicking, since switching estimators can flip the sign of real studies' conclusions. The prescription common to all three is exactly Theorem 7: if the target is the proportional change in the conditional mean, estimate the conditional mean directly.

## 6 A Full Walkthrough on Larkspur Data

The code below runs under R 4.6.1. It lists only the statements that generate the data, estimate, compute variances, and summarize the Monte Carlo; several statements are merged onto one line, and statements used only for checks are omitted. Every number in the text comes from the output of an actual run, and the result lines in the code that start with `#` copy what the run printed (with run times omitted). Variances are always the HC0 sandwich. In the code, the coefficient vectors of M2 and M3 keep the name `theta`; they correspond to $\beta$ in the text.

### 6.1 The DGP

The key settings of M1 are these: the category determines the selling probability and median sales; the same uniform random draw `U` decides whether a merchant sells; Spotlight multiplies the selling probability by 1.25 and adds 0.10 to sellers' log sales; and sales are first rounded to whole cents, with dollars equal to cents divided by 100.

```r
P1 <- list(n = 20000, categories = c("home_decor", "apparel", "jewelry", "stationery"),
           pi = c(0.35, 0.25, 0.25, 0.15), p0 = c(0.40, 0.48, 0.36, 0.30), k = 1.25,
           median0 = c(180, 260, 320, 120), sigma = 1.1, beta_tenure = 0.25, delta = 0.10)
gen_M1 <- function(par, n = par$n) {
  category <- sample.int(length(par$pi), n, replace = TRUE, prob = par$pi)
  tenure_std <- rnorm(n)
  D <- sample(rep(0:1, length.out = n))          # complete randomisation: exactly half treated
  U <- runif(n)                                  # the same draw decides selling with and without Spotlight
  eps <- rnorm(n)
  p0 <- par$p0[category]
  log_pos <- log(par$median0[category]) + par$beta_tenure * tenure_std + par$sigma * eps
  y0c <- round(100 * ifelse(U < p0, exp(log_pos), 0))
  y1c <- round(100 * ifelse(U < par$k * p0, exp(log_pos + par$delta), 0))
  cents <- ifelse(D == 1, y1c, y0c)
  data.frame(category = factor(par$categories[category], levels = par$categories),
             tenure_std = tenure_std, D = D, cents = cents, usd = cents / 100)
  # (the generator also returns the merchant id and both potential outcomes)
}
set.seed(5101)
M1 <- gen_M1(P1)
```

The parameters of the other three blocks are as follows. M2's number of listings is a lognormal variable with median 40 and log standard deviation 0.9, rounded to an integer, and the multiplicative error is lognormal with mean one and log standard deviation 0.7; M2m's $x$ follows $\mathcal{N}(0, 0.9^2)$, with an error log standard deviation of 0.5; M3's one-week visits are $v = t + 0.6e$, the second-step error has log standard deviation 0.35, and the median merchant's monthly sales are USD 1,500; M4 is NB2, with `size` equal to $1/\eta^2$.

```r
P2  <- list(n = 2000, listings_median = 40, sd_log_listings = 0.9, p_free_ship = 0.4,
            theta = c(log(1.2), 0.45, 0.20), sigma_eta = 0.7)
P2m <- list(n = 2000, sd_x = 0.9, a0 = 0.2, a1 = 0.55, a2 = -0.10, sigma_eta = 0.5)
P3  <- list(n = 2000, gamma = c(0, 0.30, 0.20), theta = c(log(1.5), 0.8), sigma_e = 0.6, sigma_eta = 0.35,
            R_mc = 2000, grid_theta1 = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2), R_grid = 1000)
P4  <- list(n = 2000, beta = c(1.4, 0.45), size = 3, R_mc = 2000,
            grid_size = c(Inf, 12, 6, 3, 1.5, 0.75), R_grid = 1000)
# gen_M2, gen_M2m, gen_M3 draw blocks M2, M2m, M3 from P2, P2m, P3 as described above, like gen_M4
gen_M4 <- function(par, n = par$n, size = par$size) {
  x <- rnorm(n)
  mu <- exp(par$beta[1] + par$beta[2] * x)
  y <- if (is.infinite(size)) rpois(n, mu) else rnbinom(n, size = size, mu = mu)   # Var = mu + mu^2 / size
  data.frame(x = x, y = as.numeric(y), mu_true = mu)
}
```

### 6.2 Block M1: log(1+y) and PPML

All regressions have the same right-hand-side variables: Spotlight, category dummies, and standardized tenure. `ci95` returns the estimate plus or minus the normal quantile times the standard error.

```r
ols_D <- function(resp, data = M1) {
  d <- data; d$resp <- resp
  f <- lm(resp ~ D + category + tenure_std, data = d)
  est <- unname(coef(f)[["D"]]); se <- sqrt(sandwich::sandwich(f)["D", "D"])
  list(est = est, se = se, ci = ci95(est, se))
}
ppml_D <- function(resp, rhs = "D + category + tenure_std") {
  d <- M1; d$resp <- resp
  f <- glm(as.formula(paste("resp ~", rhs)), family = quasipoisson(), data = d,
           control = glm.control(epsilon = 1e-12, maxit = 100))
  est <- unname(coef(f)[["D"]]); se <- sqrt(sandwich::sandwich(f)["D", "D"])
  list(est = est, se = se, ci = ci95(est, se), fit = f)
}
l_kusd <- ols_D(log1p(M1$usd / 1000)); l_usd <- ols_D(log1p(M1$usd)); l_cents <- ols_D(log1p(M1$cents))
a_usd <- ols_D(asinh(M1$usd)); a_cents <- ols_D(asinh(M1$cents))
p_usd <- ppml_D(M1$usd); p_cents <- ppml_D(M1$cents); p_raw <- ppml_D(M1$usd, rhs = "D")
truth <- log(P1$k) + P1$delta                    # log(E[y(1)] / E[y(0)]) = log 1.25 + 0.10
fe <- fixest::fepois(usd ~ D + tenure_std | category, data = M1, vcov = "hetero", glm.tol = 1e-10, fixef.tol = 1e-10)
fe_est <- unname(coef(fe)[["D"]]); fe_se <- unname(fixest::se(fe)[["D"]])   # category fixed effects instead of dummies
ext <- ols_D(as.numeric(M1$usd > 0))
pos <- M1[M1$usd > 0, ]; lp <- ols_D(log(pos$usd), data = pos)
# log(1 + y) coefficient on D: 0.0409 (USD thousands) 0.5931 (USD) 1.0581 (cents)
# asinh coefficient on D: 0.6627 (USD) 1.1281 (cents)
# PPML: 0.301561 (USD) 0.301561 (cents), HC0 SE 0.0398, 95% CI [0.2236, 0.3795]; truth 0.3231; fepois 0.301561 (SE 0.0398)
# Extensive margin 0.1011 (SE 0.0069); log y among sellers 0.0894 (SE 0.0235; delta = 0.10)
```

The population ATEs corresponding to the three log(1+y) coefficients 0.0409, 0.5931, and 1.0581 are 0.0412, 0.5825, and 1.0378, and each of the three intervals covers its own population counterpart. The estimators are not wrong; what is wrong is reading them as the same percentage. PPML's 0.3016 is 0.022 below the truth, within one standard error. PPML on $D$ alone gives 0.305, the log of the ratio of the two group means. If the standard error is computed under the Poisson assumption that the variance equals the mean, it is only 0.0010 with dollar accounting, about one fortieth of the sandwich's 0.0398. This number also changes with the unit, because whether the variance equals the mean itself depends on the unit; the sandwich is unaffected.

### 6.3 Blocks M2 and M2m: Gauss-Newton and Pseudo-True Values

`gauss_newton` runs OLS of the residuals on the gradient at each step and uses step halving to keep the residual sum of squares from rising; `m_vcov` is the $\hat A^{-1}\hat B\hat A^{-1}/n$ of equation (12.47).

```r
gauss_newton <- function(y, X, mfun, jac, theta0, tol = 1e-10, maxit = 200) {
  theta <- as.numeric(theta0)
  ssr <- function(th) sum((y - mfun(X, th))^2)
  s_cur <- ssr(theta); small <- sqrt(.Machine$double.eps); converged <- FALSE
  for (iter in seq_len(maxit)) {
    step <- as.numeric(qr.coef(qr(jac(X, theta)), y - mfun(X, theta)))   # OLS of residuals on the gradient
    if (anyNA(step)) break
    lambda <- 1
    repeat {                                                             # halve until the SSR does not rise
      cand <- theta + lambda * step; s_new <- ssr(cand)
      if (is.finite(s_new) && (s_new <= s_cur * (1 + 1e-12) || max(abs(step)) < small * (1 + max(abs(theta))))) break
      lambda <- lambda / 2
      if (lambda < 2^-30) break
    }
    if (!is.finite(s_new) || lambda < 2^-30) break
    theta <- cand; s_cur <- s_new
    if (max(abs(step)) < tol * (1 + max(abs(theta)))) { converged <- TRUE; break }
  }
  list(theta = theta, iter = iter, converged = converged)
}
m_vcov <- function(S, A) {            # S: n x k scores at theta_hat; A: k x k mean Hessian
  Ai <- solve(A)
  Ai %*% (crossprod(S) / nrow(S)) %*% t(Ai) / nrow(S)
}
exp_mean <- function(X, b) exp(drop(X %*% b))          # m = exp(X b)
exp_jac <- function(X, b) exp(drop(X %*% b)) * X       # its n x k gradient, row i = m_i x_i'
X2 <- cbind(1, M2$x1, M2$free_ship)
gn <- gauss_newton(M2$y, X2, exp_mean, exp_jac, c(0, 0, 0))
fit <- nls(y ~ exp(b0 + b1 * x1 + b2 * free_ship), data = M2, start = list(b0 = 0, b1 = 0, b2 = 0),
           control = nls.control(tol = 1e-8, maxiter = 200))
p_out <- nls_parts(M2$y, X2, gn$theta)          # NLS scores and A = mean of grad m' grad m
V_robust <- m_vcov(p_out$S, p_out$A)            # equals sandwich::sandwich(fit)
V_model <- vcov(fit)                            # sigma^2 (sum grad m' grad m)^-1, needs NLS.3
# M2: Gauss-Newton (11 iterations) 0.11658 0.45367 0.27166; nls 0.11658 0.45367 0.27166; max difference 1.20e-08
# M2: SE model 0.0303 0.0197 0.0360 | robust 0.0328 0.0331 0.0523
Xm <- cbind(1, M2m$x)
gm <- gauss_newton(M2m$y, Xm, exp_mean, exp_jac, c(0, 0.5))
pm_full <- nls_parts(M2m$y, Xm, gm$theta, hessian = "full")   # the mean may be wrong: full Hessian
pm_out <- nls_parts(M2m$y, Xm, gm$theta)                       # semirobust: A = mean of grad m' grad m
Vm_full <- m_vcov(pm_full$S, pm_full$A); Vm_semi <- m_vcov(pm_out$S, pm_out$A)
gpm <- glm(y ~ x, family = quasipoisson(), data = M2m, control = glm.control(epsilon = 1e-12, maxit = 100))
m0 <- function(x) exp(0.2 + 0.55 * x - 0.10 * x^2)              # the true mean of M2m
pseudo_nls <- pseudo_true_exp(m0, 0.9, "nls"); pseudo_ppml <- pseudo_true_exp(m0, 0.9, "ppml")   # quadrature
# M2m: b1 = 0.4406, fully robust SE 0.0163 (semirobust 0.0177), CI [0.4086, 0.4727]; pseudo-true 0.4154; PPML 0.4838 (pseudo-true 0.4733)
```

The sandwich standard errors of the three coefficients are 1.083, 1.679, and 1.451 times the homoskedastic ones; the more the multiplicative error "amplifies" large merchants, the bigger the gap. The pseudo-true values are obtained by numerical integration on a grid of 40,001 points: $(0.1702, 0.4154)$ for NLS and $(0.1396, 0.4733)$ for PPML. Running Gauss-Newton directly on a million noise-free points gives $(0.1700, 0.4159)$, in agreement with the integration.

### 6.4 Block M3: Correcting the Variance of a Two-Step Estimator

`twostep_fit` runs OLS in the first step and stores the influence function, then runs NLS on $\hat t$ in the second step; `F2` is $\hat\Gamma$ computed analytically, and row $i$ of `R1` is $\hat\psi_i'$. `two_step_vcov` simply replaces the filling with the outer product of $\hat s_i + \hat\Gamma\hat\psi_i$.

```r
twostep_fit <- function(d, theta_start = c(0, 0.5)) {
  Z <- cbind(1, d$z1, d$z2); n <- nrow(Z); ZZ <- crossprod(Z) / n
  gamma <- as.numeric(solve(ZZ, crossprod(Z, d$v) / n))       # step 1: OLS of v on z
  e <- d$v - drop(Z %*% gamma)
  R1 <- (Z * e) %*% solve(ZZ)                                # row i: influence function psi_i'
  X <- cbind(1, drop(Z %*% gamma))                           # generated regressor t_hat
  gn <- gauss_newton(d$y, X, exp_mean, exp_jac, theta_start)  # step 2: NLS
  theta <- gn$theta
  p <- nls_parts(d$y, X, theta); mu <- p$m * p$u
  F2 <- (-theta[2] * crossprod(X * mu, Z) - rbind(0, crossprod(mu, Z)) +
           theta[2] * crossprod(X * p$m^2, Z)) / n             # Gamma_hat: n^-1 sum of grad_gamma s_i
  list(gamma = gamma, theta = theta, V_naive = m_vcov(p$S, p$A),
       V_corrected = two_step_vcov(p$S, p$A, F2, R1))         # (other returned fields omitted)
}
two_step_vcov <- function(S2, A2, F2, R1) m_vcov(S2 + R1 %*% t(F2), A2)
nls_on <- function(xvar) {                                   # one NLS of y on exp(theta0 + theta1 xvar)
  X <- cbind(1, xvar); g <- gauss_newton(M3$y, X, exp_mean, exp_jac, c(0, 0.5))
  p <- nls_parts(M3$y, X, g$theta); se <- sqrt(m_vcov(p$S, p$A)[2, 2])
  list(est = g$theta[2], se = se, ci = ci95(g$theta[2], se))
}
f <- twostep_fit(M3); direct <- nls_on(M3$v); oracle <- nls_on(M3$t_true)
# M3: theta1_hat = 0.8064; naive SE 0.0302, corrected SE 0.0431 (ratio 1.425); direct NLS on v 0.2088; on true t 0.7750
one_rep <- function(par) {                                   # a fresh M3 sample, both steps refitted
  r <- twostep_fit(gen_M3(par))
  c(est = r$theta[2], se_naive = sqrt(r$V_naive[2, 2]), se_corrected = sqrt(r$V_corrected[2, 2]))
}
set.seed(5114); sim <- t(vapply(seq_len(P3$R_mc), function(r) one_rep(P3), numeric(3)))
cover <- function(se) mean(abs(sim[, "est"] - P3$theta[2]) <= qnorm(0.975) * sim[, se])
c(cover("se_naive"), cover("se_corrected"), sd(sim[, "est"]), colMeans(sim[, c("se_naive", "se_corrected")]))
# M3 Monte Carlo (2000 reps): coverage naive 0.8235, corrected 0.9490; sd of estimates 0.0432; mean SE 0.0298 / 0.0423
```

The first-step coefficient estimates are $(-0.0106, 0.294, 0.185)$, against true values $(0, 0.30, 0.20)$. If the true $t$ could be used directly, NLS would give 0.775, and its interval [0.718, 0.832] also covers 0.8. The corrected standard error of the two-step estimator, 0.0431, is larger than the naive 0.0302, and the excess reflects the first step's uncertainty. With the true $t$, the standard error is 0.0291, close to the naive 0.0302.

### 6.5 Block M4: Overdispersion and the Sandwich

```r
poisson_parts <- function(y, X, b) {             # q = mu - y log(mu), mu = exp(X b)
  mu <- exp_mean(X, b)
  list(S = (mu - y) * X, A = crossprod(X * mu, X) / length(y), mu = mu)
}
g <- glm(y ~ x, family = poisson(), data = M4, control = glm.control(epsilon = 1e-12, maxit = 100))
pp <- poisson_parts(M4$y, cbind(1, M4$x), coef(g))
se_model <- sqrt(vcov(g)[2, 2])                  # A^-1 / n: assumes Var(y | x) = E[y | x]
se_sand <- sqrt(m_vcov(pp$S, pp$A)[2, 2])        # A^-1 B A^-1 / n
gq <- glm(y ~ x, family = quasipoisson(), data = M4, control = glm.control(epsilon = 1e-12, maxit = 100))
se_quasi <- sqrt(vcov(gq)[2, 2])                 # model variance times the Pearson dispersion
phi <- sum((M4$y - pp$mu)^2 / pp$mu) / (nrow(M4) - 2)   # Pearson dispersion (1 under Poisson)
# M4 sample: beta1_hat 0.4926; SE model 0.0104, sandwich 0.0206, quasi 0.0165; Pearson dispersion 2.491
one_rep <- function(size) {                                  # a fresh M4 sample at NB2 size `size`
  d <- gen_M4(P4, size = size); Xr <- cbind(1, d$x)
  b <- glm.fit(Xr, d$y, family = poisson(), control = list(epsilon = 1e-10, maxit = 50))$coefficients
  p <- poisson_parts(d$y, Xr, b); v_model <- solve(p$A)[2, 2] / P4$n
  phi <- sum((d$y - p$mu)^2 / p$mu) / (P4$n - 2)              # Pearson dispersion
  c(est = b[[2]], se_model = sqrt(v_model), se_sandwich = sqrt(m_vcov(p$S, p$A)[2, 2]), se_quasi = sqrt(v_model * phi))
}
set.seed(5115); sim <- t(vapply(seq_len(P4$R_mc), function(r) one_rep(P4$size), numeric(4)))
cover <- function(se) mean(abs(sim[, "est"] - P4$beta[2]) <= qnorm(0.975) * sim[, se])
c(sapply(c("se_model", "se_sandwich", "se_quasi"), cover), sd(sim[, "est"]), colMeans(sim[, c("se_model", "se_sandwich")]))
# M4 Monte Carlo (2000 reps, size 3): coverage model 0.7290, sandwich 0.9440, quasi 0.9130; sd 0.0192; mean SE 0.0106 / 0.0187
```

In the sample, merchants with zero orders make up a share of 0.091, and the average is 4.57 orders per merchant. The Pearson dispersion of 2.49 aims at the population quantity $\mathbb{E}[\mathrm{Var}(y \mid x)/\mathbb{E}[y \mid x]] = 2.50$, which is not the same number as the ratio of the average conditional variance to the average mean (2.83) used on the horizontal axis. This sample's sandwich interval [0.452, 0.533] does not cover 0.45: the estimate is 2.07 sandwich standard errors from the truth. Coverage has to be read from the Monte Carlo (figure below): the true standard deviation, 0.0192, is 1.81 times the average model standard error. When the ratio of the average conditional variance to the average mean is pushed up to 8.33, the coverage of the model standard error drops to 0.452, while that of the sandwich is still 0.948. When this ratio is 1.46, the sandwich's 0.929 is slightly below the error band in the figure, about 3 Monte Carlo standard errors below the nominal level (the Monte Carlo standard error of a coverage rate is about 0.0069 with 1,000 replications and 0.0049 with 2,000). This is most likely Monte Carlo noise, but it may also reflect a tendency of HC0 to be slightly too small in finite samples: at 5 of the 6 grid points, the average sandwich standard error is below the Monte Carlo standard deviation of the estimates, and at this grid point it is 0.0131 against 0.0138.

![Monte Carlo coverage of nominal 95% intervals. One panel is M4's Poisson PML, with $\mathbb{E}[\mathrm{Var}(y \mid x)]/\mathbb{E}[y]$ on the horizontal axis; the other is M3's two-step NLS, with the true traffic elasticity on the horizontal axis. The solid line is the nominal level 0.95, and the two dotted lines above and below it are the Monte Carlo error band for 1,000 replications.](assets/fig/fig_05_coverage.svg)

### 6.6 The Grand Reconciliation

| Block | Estimator | Estimate | 95% interval | Target and true value | Covered |
|---|---|---|---|---|---|
| M1 | log(1+y), dollars | 0.593 | [0.516, 0.670] | Population ATE 0.583 | Yes |
| M1 | log(1+y), cents | 1.058 | [0.920, 1.196] | Population ATE 1.038 | Yes |
| M1 | PPML, dollars or cents | 0.302 | [0.224, 0.379] | Log ratio 0.323 | Yes |
| M1 | Extensive-margin effect (LPM) | 0.1011 | [0.0875, 0.1146] | 0.099 | Yes |
| M1 | log y among sellers | 0.089 | [0.043, 0.136] | Intensive-margin effect 0.10 | Yes |
| M2 | NLS intercept | 0.1166 | [0.0523, 0.1809] | 0.1823 | No |
| M2 | NLS listings elasticity | 0.454 | [0.389, 0.518] | 0.45 | Yes |
| M2 | NLS free-shipping coefficient | 0.272 | [0.169, 0.374] | 0.20 | Yes |
| M2m | NLS slope | 0.441 | [0.409, 0.473] | Pseudo-true value 0.415 | Yes |
| M2m | PPML slope | 0.484 | [0.458, 0.510] | Pseudo-true value 0.473 | Yes |
| M3 | Two-step NLS, naive standard error | 0.806 | [0.747, 0.866] | 0.8 | Yes |
| M3 | Two-step NLS, corrected standard error | 0.806 | [0.722, 0.891] | 0.8 | Yes |
| M3 | NLS using $v$ directly | 0.209 | [0.176, 0.242] | 0.8 | No |
| M4 | PPML, model standard error | 0.493 | [0.472, 0.513] | 0.45 | No |
| M4 | PPML, sandwich | 0.493 | [0.452, 0.533] | 0.45 | No |

Two single-sample intervals in the table should have covered the truth and did not: M2's intercept and M4's sandwich interval. A nominal 95% interval is supposed to miss about once in twenty, and a hit or miss in one sample says nothing about the quality of a variance formula. So these two rows stay in the table as they are, and conclusions about coverage are drawn only from the Monte Carlo. On M4 the sandwich covers 94.4% and the model standard error 72.9%, which in a single sample shows up as a standard-error ratio of 1.97; on M3 the corrected standard error covers 94.9% and the naive one 82.35%, a standard-error ratio of 1.425. The miss in the row that uses $v$ directly, by contrast, is bias and has nothing to do with luck.

The reconciliation of this section comes down to one sentence: every estimator aims at the population quantity defined by its loss function, and when things go wrong, either that quantity was misread as something else or the variance was computed too small.

## 7 Failure Modes and Robustness

Mean misspecification is the most fundamental threat to M-estimation, but its consequences are often misstated. Theorem 5 and Section 3.7 both show that the point estimate does not "diverge"; it aims at the best approximation that the loss function picks out, and the trouble is that this approximation may not be the number you want to report. There are three practical responses. The first is to make the mean flexible enough: Wooldridge (2010) §18.2 recommends putting squares, interactions, and logs of the regressors inside the exponential, which is often enough. The second is to run specification tests: tests of the RESET kind, which Santos Silva and Tenreyro (2006) use, can be run for both NLS and PPML. The third is to look at NLS and PPML together: when the mean is correctly specified the two have the same probability limit. On M2m their pseudo-true values differ by 0.058 and the sample estimates are 0.441 and 0.484; a gap like this is itself an alarm. For inference, when the mean may be wrong, use the fully robust variance with the sample Hessian.

Reporting the wrong standard error is the most common trap, and it is harder to spot than a wrong point estimate, because the point estimate looks perfectly normal. The core discipline is to think through which variance you are reporting: the model's default standard error is almost always the collapsed version, which is right only when $B_0 = \sigma^2 A_0$, and this requires homoskedasticity, or a variance proportional to the mean, or a fully correct likelihood, conditions that almost never hold for outcomes such as sales, order counts, or patent counts. On M1 the Poisson default standard error is about one fortieth of the sandwich, and it also changes with the unit. The default choice should be the sandwich. Note that the heteroskedasticity-robust sandwich corrects only one thing, "the observations have unequal variances," and still assumes that different observations are independent. When the data have a group or time structure, a cluster-robust variance is needed as well; that is the subject of Chapter 8. A recurring lesson: a beautiful, suspiciously narrow standard error often reflects not a precise estimate but the wrong variance formula.

The most common mistake in two-step estimation is to treat the first step as known. When $\Gamma_0 \neq 0$ the naive standard error is generally wrong, but the direction of the error depends on the case. In the setting of Theorem 6 (ii), that is, when the scores of the two steps are uncorrelated, it is necessarily too small, as in M3, where the coverage of the naive interval is only 82.35%. In general, the scores of the two steps are correlated, and the difference between the corrected and uncorrected variances can be indefinite (larger in some directions and smaller in others) or even negative semidefinite, in which case estimating the first step actually makes the second step more precise (Wooldridge 2010 p. 413). Whatever the direction, the remedy is to correct with equations (12.61) and (12.62), or to bootstrap the whole two-step procedure, re-estimating the first step in every resample. When the first step itself is misspecified, $\gamma_* \neq \gamma_0$, the second-step point estimate aims at something else, and correcting the variance cannot save it. Even a correctly specified first step may not be enough: if the estimated index is not an exact function of $z$, plugging fitted values into a nonlinear second step also introduces bias (Section 3.6). Using the noisy visits directly in M3 is a different matter: there is no first step at all, and it is measurement error that attenuates the elasticity to 0.209, which, as in Section 3.6, is a bias problem rather than a variance problem.

The risk with outcomes that have zeros lies in interpretation. Any regression that uses log(1+y), log(c+y), or asinh should report the unit of the outcome and read the coefficient as "the average effect on some transformation" rather than as a percentage. Chen and Roth (2024) point out that checking how sensitive the coefficient is to the unit or to the constant $c$ amounts, in essence, to an indirect measure of the size of the extensive margin; on M1, the change of 0.465 in the coefficient from dollars to cents is almost exactly the extensive-margin effect times $\log 100$. If you want a percentage, estimate the proportional change in the average; if both margins matter, report them separately, and say clearly under what assumptions a comparison among sellers is causal. In panel data, the Poisson fixed effects estimator remains consistent with a fixed number of periods when the regressors are strictly exogenous conditional on the multiplicative individual effect; it again relies only on the conditional mean (Wooldridge 1999), which is one reason PPML is popular in applications with many fixed effects.

Finally, there are identification and numerical problems. A nearly singular Hessian, an iteration that does not converge, and an answer that changes with the starting point are often signs of weak identification rather than faults of the software. Try NLS from several starting points, and check whether Gauss-Newton and a general-purpose optimizer land on the same point.

String these failure modes together, and the robustness of this chapter's methods rests, in the end, on two things: first, knowing which part of the model you got right (the mean, the variance, or the whole distribution) and judging from that what the point estimate still estimates and which set of standard errors to use; second, thinking through which observations your loss function lets dominate. Think these two things through and estimation is stable; fail to, and even the most automated software is merely handing you, quickly, a number that may be wrong in a very respectable way.

## 8 Further Reading

::: {.readings}
Required reading, in suggested reading order:

- Wooldridge (2010, *Econometric Analysis of Cross Section and Panel Data*, 2nd ed.), Chapters 12 to 14. The unified treatment of M-estimation, MLE, and GMM; Theorems 1 to 3 of this chapter are its Theorems 12.1 to 12.3, Theorems 4 to 6 are taken from §12.2 to §12.5, and the Poisson part is in §18.2.
- Santos Silva and Tenreyro (2006, Review of Economics and Statistics). PPML and the Jensen bias of log-linear models.
- Chen and Roth (2024, Quarterly Journal of Economics). The unit dependence of log(1+y) and asinh, the trilemma, and how several alternative parameters can be estimated; for the correction to the trilemma, see the corrigendum of Okat, Chen and Roth (2026).
- Newey and McFadden (1994, *Handbook of Econometrics* Vol. 4). The authoritative technical treatment of the consistency and asymptotic normality of extremum estimators, where the regularity conditions that Sections 3.3 and 3.4 of this chapter pass over are proved in full.

Further reading:

- Huber (1967) and White (1980, Econometrica). The two sources of the sandwich variance: the former treats maximum likelihood estimation under nonstandard conditions, the latter gives the heteroskedasticity-robust variance for linear regression.
- White (1981, Journal of the American Statistical Association). NLS with a misspecified mean and its pseudo-true value.
- Gourieroux, Monfort and Trognon (1984a, 1984b, Econometrica). PML in the linear exponential family. The 1984b paper gives the general theory and is the theoretical source of the result that Poisson PML needs only a correct mean; the 1984a paper applies it to Poisson models.
- Newey (1984, Economics Letters), Pagan (1984, International Economic Review), and Murphy and Topel (1985, Journal of Business & Economic Statistics). Three classic treatments of two-step estimation and generated regressors: stacked moment conditions, inference with generated regressors, and the two-step variance formula.
- Wooldridge (1999, Journal of Econometrics). Poisson fixed effects in panel data likewise rely only on the conditional mean.
- Mullahy and Norton (2024, Oxford Bulletin of Economics and Statistics) and Bellemare and Wichman (2020, Oxford Bulletin of Economics and Statistics). The former discusses the pitfalls of transformed outcomes with a mass at zero; the latter is a widely used reference for interpreting asinh elasticities and one of the papers Chen and Roth discuss.
- Cohn, Liu and Wardlaw (2022, Journal of Financial Economics). The difference between log(1+y) and Poisson regression in real studies.
- Amemiya (1985, *Advanced Econometrics*). The classic textbook source for the name and the unified framework of the extremum estimator.
- Cameron and Trivedi (2013, *Regression Analysis of Count Data*, 2nd ed.). The standard reference for count models, overdispersion, and NB2.
- Hansen (2014) §9.1. Another textbook treatment of the asymptotic distribution of NLS, with a proof sketch.
:::

::: {.apa-refs}
- Amemiya, T. (1985). *Advanced econometrics*. Harvard University Press.
- Bellemare, M. F., & Wichman, C. J. (2020). Elasticities and the inverse hyperbolic sine transformation. *Oxford Bulletin of Economics and Statistics, 82*(1), 50-61. https://doi.org/10.1111/obes.12325
- Berndt, E. R., Hall, B. H., Hall, R. E., & Hausman, J. A. (1974). Estimation and inference in nonlinear structural models. *Annals of Economic and Social Measurement, 3*(4), 653-665.
- Cameron, A. C., & Trivedi, P. K. (2013). *Regression analysis of count data* (2nd ed.). Cambridge University Press.
- Carranza, E., Garlick, R., Orkin, K., & Rankin, N. (2022). Job search and hiring with limited information about workseekers' skills. *American Economic Review, 112*(11), 3547-3583. https://doi.org/10.1257/aer.20200961
- Chen, J., & Roth, J. (2024). Logs with zeros? Some problems and solutions. *The Quarterly Journal of Economics, 139*(2), 891-936. https://doi.org/10.1093/qje/qjad054
- Cohn, J. B., Liu, Z., & Wardlaw, M. I. (2022). Count (and count-like) data in finance. *Journal of Financial Economics, 146*(2), 529-551. https://doi.org/10.1016/j.jfineco.2022.08.004
- Gourieroux, C., Monfort, A., & Trognon, A. (1984a). Pseudo maximum likelihood methods: Applications to Poisson models. *Econometrica, 52*(3), 701-720. https://doi.org/10.2307/1913472
- Gourieroux, C., Monfort, A., & Trognon, A. (1984b). Pseudo maximum likelihood methods: Theory. *Econometrica, 52*(3), 681-700. https://doi.org/10.2307/1913471
- Hansen, B. E. (2014). *Econometrics* [Unpublished manuscript]. Department of Economics, University of Wisconsin.
- Huber, P. J. (1967). The behavior of maximum likelihood estimates under nonstandard conditions. In L. M. Le Cam & J. Neyman (Eds.), *Proceedings of the Fifth Berkeley Symposium on Mathematical Statistics and Probability* (Vol. 1, pp. 221-233). University of California Press.
- Mullahy, J., & Norton, E. C. (2024). Why transform Y? The pitfalls of transformed regressions with a mass at zero. *Oxford Bulletin of Economics and Statistics, 86*(2), 417-447. https://doi.org/10.1111/obes.12583
- Murphy, K. M., & Topel, R. H. (1985). Estimation and inference in two-step econometric models. *Journal of Business & Economic Statistics, 3*(4), 370-379. https://doi.org/10.1080/07350015.1985.10509471
- Newey, W. K. (1984). A method of moments interpretation of sequential estimators. *Economics Letters, 14*(2-3), 201-206. https://doi.org/10.1016/0165-1765(84)90083-1
- Newey, W. K., & McFadden, D. (1994). Large sample estimation and hypothesis testing. In R. F. Engle & D. L. McFadden (Eds.), *Handbook of econometrics* (Vol. 4, pp. 2111-2245). Elsevier. https://doi.org/10.1016/S1573-4412(05)80005-4
- Okat, D., Chen, J., & Roth, J. (2026). *Corrigendum: Logs with zeros? Some problems and solutions* [Unpublished manuscript]. https://www.jonathandroth.com/assets/files/LogUniqueHOD0-corrigendum.pdf
- Pagan, A. (1984). Econometric issues in the analysis of regressions with generated regressors. *International Economic Review, 25*(1), 221-247. https://doi.org/10.2307/2648877
- Santos Silva, J. M. C., & Tenreyro, S. (2006). The log of gravity. *The Review of Economics and Statistics, 88*(4), 641-658. https://doi.org/10.1162/rest.88.4.641
- White, H. (1980). A heteroskedasticity-consistent covariance matrix estimator and a direct test for heteroskedasticity. *Econometrica, 48*(4), 817-838. https://doi.org/10.2307/1912934
- White, H. (1981). Consequences and detection of misspecified nonlinear regression models. *Journal of the American Statistical Association, 76*(374), 419-433. https://doi.org/10.1080/01621459.1981.10477663
- Wooldridge, J. M. (1999). Distribution-free estimation of some nonlinear panel data models. *Journal of Econometrics, 90*(1), 77-97. https://doi.org/10.1016/S0304-4076(98)00033-5
- Wooldridge, J. M. (2010). *Econometric analysis of cross section and panel data* (2nd ed.). MIT Press.
:::
