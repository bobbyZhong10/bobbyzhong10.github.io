---
title: "GMM & Minimum Distance"
subtitle: "Estimating from Moment Conditions"
seriesline: "Foundations of Information Systems Economics · Chapter 7"
series: "Foundations of Information Systems Economics: Empirical Methods for Information Systems Economics"
sidebartitle: "Foundations of Information Systems Economics"
sidebarsub: "Theory · Methods · Applications"
chapterlabel: "Chapter 7 · GMM & Minimum Distance"
---

## Introduction

Tamarind is a fictional US online tutoring marketplace: households buy lesson credits, one credit books one lesson, and base prices lie between USD 8.99 and USD 11.49. All the pricing team wants is one number: the price elasticity of demand for lessons. To get it, the platform has run three batches of experiments. The first batch ran a price A/B test in each of 16 cities (households are split at random into two groups, each paying its own price), and each city reports a price slope. The second ran ten randomized pricing experiments simultaneously in 200 regions, but regional managers also run ad hoc promotions of their own, so prices are not set by the experiments alone. The third batch of regions has only list-price perturbations of a few cents. Section 1 shows that a direct regression of lessons on price yields an elasticity only about a fifth of the truth.

This chapter covers GMM and minimum distance (MD). Section 3.1 of Chapter 5 listed GMM as a member of the extremum-estimator family, and Chapter 6 started from a complete distribution. This chapter takes the other path: the economic model supplies only a few moment conditions, each saying that some expectation is zero at the true parameter, and estimation is built on them, with no need to write down a complete distribution. The three batches of experiments correspond to the chapter's three threads. First, when there are more moment conditions than parameters, how are they combined into one estimate, and how do we test whether they are mutually compatible? The answers are GMM, the efficient weight, and the J test. Second, how does conventional inference fail when the moment conditions are weak, and what kind of confidence set is still honest? The answer is the S confidence set of Stock and Wright (2000). Third, when several markets each provide a reduced-form estimate, how are these combined into one structural parameter? The answer is MD, which is asymptotically equivalent to GMM.

The chapter proceeds as follows. Section 2 sets out Tamarind's demand model and the estimands and moment conditions of the four data blocks. Section 3, the most detailed in the chapter, gives in turn the consistency and asymptotic normality of GMM, the efficient weight, the J test, the finite-sample problems of two-step estimation, the S confidence set under weak moments, and MD and its equivalence with GMM. Section 4 covers implementation, the choice of moments, and software, Section 5 examines three anchoring papers, Section 6 runs a full walkthrough on Tamarind's data, and Section 7 discusses failure modes. GMM, the J test, and MD all return in Chapter 9 (instrumental variables and system estimation), Chapter 10 (Chamberlain's approach for panel data), and Chapter 11 (dynamic panels).

## 1 An Elasticity Only a Fifth of the Truth

::: {.case}
Block G2 is one week of Tamarind data from 200 regions. The outcome is weekly lessons booked per thousand households, and the price is the average price per lesson actually paid in the region; both are in logs. Prices have two sources. One is ten concurrent randomized pricing experiments, such as a booking fee, a deep discount, and a weekend surcharge: each experiment puts some regions in a high-price arm and the rest in a low-price arm, and the two arms differ in price by 1% to 6%. The other is the regional managers' ad hoc promotions: in a slow week, managers offer discounts on their own initiative.

Because this is a simulation, we know the answer: the price elasticity of demand for lessons is $-1.6$, so a price 1% higher means lessons booked drop by 1.6%.
:::

The analyst starts with the most direct approach: an OLS regression of log lessons on log price. The slope is $-0.318$ (robust standard error 0.290), read as "raise prices 10% and lessons fall by only about 3%": demand is almost inelastic, and on this basis the pricing team argues for a large price increase. This is not sampling error: in the DGP's population, the probability limit of this regression is $-0.382$. The culprit is the managers' promotions. They fall in the very weeks when the demand shock is negative, so low prices and low demand occur together, diluting the negative price effect.

Change the estimator and most of the problem disappears. The ten experiments were randomly assigned, are unrelated to the demand shock, and affect lessons only through price. Write "experiment assignment is uncorrelated with the demand shock" as ten moment conditions, add one for the intercept, and two-step GMM gives $-2.012$ (standard error 0.358), with a 95% interval $[-2.713, -1.311]$ that covers $-1.6$; 2SLS (two-stage least squares) gives $-1.843$ (0.414). New questions follow at once. With eleven moment conditions and two parameters, no elasticity in the sample can set all the moment conditions to zero at once, so how should they be weighted? Can the standard error of 0.358 be trusted? In another batch of 1,000 regions, where the experiments are only list-price perturbations of a few cents, the same two-step GMM gives $-0.229$, and its 95% interval $[-1.482, 1.025]$ excludes the truth. These are the three questions this chapter answers: how to weight and test the surplus moment conditions; why the standard errors of two-step GMM are too small in small samples; and what kind of confidence set is still honest when the moment conditions are weak.

## 2 The Economic Model and the Estimand

### 2.1 Demand, Prices, and Moment Conditions

For region $i$, let $y_i$ be log weekly lessons booked and $\log p_i$ the log price paid; the demand equation is

$$y_i = \beta_0 + \theta\,\log p_i + \varepsilon_i,$$

where $\theta$ is the price elasticity, with true value $\theta_0 = -1.6$, and $\varepsilon_i$ is the demand shock. Prices are set jointly by the experiments and the promotions: $\log p_i = \log 10.49 + \sum_{j=1}^{10} \gamma_j z_{ij} + v_i$. Here $z_{ij}$ is the assignment of region $i$ in experiment $j$: positive in the experimental arm (the arm with the larger fee, surcharge, discount, rebate, or threshold) and negative otherwise, centered by the assignment share so that its sample mean is exactly zero. $\gamma_j$ is how far one unit of the experiment (a dollar or a percentage point) pushes the price; for the discount and rebate experiments $\gamma_j$ is negative, and the arm with the positive value actually has the lower price. $v_i$ is the managers' promotion, positively correlated with the demand shock, which is why OLS is biased toward zero. Among regions with the same noise, the correlation between the two is 0.8; pooled over all of G2's regions, the correlation is lower because the noise differs across regions.

Write $z_i = (1, z_{i1}, \ldots, z_{i,10})'$. Random assignment guarantees that experiment assignment is unrelated to the demand shock; if we assume in addition that the experiments affect lessons only through price (the exclusion restriction), we obtain 11 moment conditions:

$$\mathbb{E}\big[z_i\,(y_i - \beta_0 - \theta\,\log p_i)\big] = 0.$$

The parameters are $(\beta_0, \theta)$, so $k = 2$. The exclusion restriction is an economic assumption that randomization cannot guarantee: if the premium-tutor upcharge experiment also changed parents' perception of lesson quality, it would bypass price and affect lessons directly. In G2's DGP the exclusion restriction holds. Each experiment on its own yields an elasticity, its effect on lessons divided by its effect on price:

$$\theta_0 = \frac{\mathrm{Cov}(z_{ij}, y_i)}{\mathrm{Cov}(z_{ij}, \log p_i)}, \qquad j = 1, \ldots, 10.$$

The numerator is the experiment's effect on lessons and the denominator its effect on price; a nonzero denominator is what "the experiment really moved the price" means. The ten experiments give ten such ratios, all equal to $-1.6$ in the population. GMM answers the question of how to combine these ten ratios, which differ from one another in the sample, into a single estimate.

We should also be clear about what this elasticity is. In the DGP demand has constant elasticity, $-1.6$ everywhere; real demand need not. The experiments move prices only a few percentage points around the base price, so they identify the elasticity over this stretch of prices, not the shape of the demand curve far from it. In the anchoring paper of Section 5.3, the same demand curve has an elasticity of $-0.33$ at USD 99 and of $-1.15$ at USD 399; extrapolating a local elasticity to a large price increase is where the constant-elasticity assumption is most likely to go wrong.

The estimand here is the parameter that makes all the moment conditions hold at once. In the general statements of Section 3, $\theta$ denotes all the parameters ($\theta \in \mathbb{R}^k$), and in G2 it is the elasticity; the moment function is written $g(w,\theta)$, which in G2 is $g(w_i,\theta) = z_i(y_i - \beta_0 - \theta\log p_i)$ and has $r = 11$ components, so $r - k = 9$ of the moment conditions are "surplus."

### 2.2 The Other Three Data Blocks

Block G1 consists of price A/B tests in 16 cities. In city $m$, 2,000 to 8,000 households take part; a random half pays the base price (USD 8.99 to USD 11.49) and the other half pays a price 3% to 12% higher, and the outcome is the log of the household's average weekly lessons booked over 8 weeks. Price is randomized within the test, so the OLS slope of city $m$'s outcome on log price, $\hat\pi_m$, estimates that city's elasticity, written $\pi_m$, with standard error $\hat\sigma_m$. In this chapter and in Chapters 9 and 10, $\pi$ and $\Pi$ are reduced-form parameters, not the profits or payoffs of Chapters 1 to 4. In G1 every city's elasticity is $-1.6$, so $\pi_0 = h(\theta_0) = \theta_0\,\boldsymbol\iota_{16}$ (where $\boldsymbol\iota_{16}$ is the 16-dimensional vector whose entries all equal 1; the subscript gives the dimension, as in $\boldsymbol\iota_r$ and $\boldsymbol\iota_T$ later), $r = 16$, and $k = 1$; here $h$ is the map from structural to reduced-form parameters and has nothing to do with the hazard $h(t \mid x)$ of Chapter 6. The precision of the 16 slopes varies widely: the true standard errors range from 0.155 (Seattle) to 1.272 (Columbus), depending on the number of participating households, the size of the price difference, and the noise in the outcome.

Block G1h is the same 16 tests, but the elasticity varies across cities: cities with higher incomes are less elastic, $\theta_m = -1.6 + 3(\log \mathrm{income}_m - \overline{\log \mathrm{income}})$, ranging from $-2.223$ (Houston) to $-0.672$ (Seattle). Now no $\theta$ can make $\pi = \theta\,\boldsymbol\iota_{16}$ hold, and "the elasticity" needs to be redefined. There are at least three candidates: the simple average over the 16 cities, $-1.6$; the average weighted by the platform's number of active households, $-1.619$; and the MD target that Section 3.7 arrives at, $-1.395$. Which to choose is a decision problem: if the platform changes prices uniformly nationwide, what matters is the average weighted by revenue or by households; if it wants to pick a few cities in which to raise prices first, what matters is each city's own elasticity. The estimator, however, will quietly choose one for you.

Block G3 has 1,000 regions, and the only experiments are four list-price perturbations of plus or minus 0.8, 1.2, 1.6, and 2.0 cents per lesson, each with a 50/50 split; demand and promotions are as in G2, and the demand shock is homoskedastic. Each of the four perturbations moves the price up or down by 0.08% to 0.19%, an amount almost drowned out by the promotions. The concentration parameter $\mu^2$ measures this: the sample size times the price variance explained by the experiments, divided by the variance of the promotions. In G2, $\mu^2$ is about 500; in G3 it is only 8.72, and the corresponding population first-stage F statistic (the F statistic from regressing price on the four experiments and testing that the four coefficients are jointly zero) is about 3.18. The experiment assignments have a sample mean of exactly zero, so the intercept drops out of G3's moment conditions of its own accord, leaving $r = 4$ and $k = 1$.

The estimands of the four data blocks can be written in one form: the parameter that makes the population moment conditions hold (G2, G3), or that makes the reduced form equal $h(\theta)$ (G1). When the moment conditions cannot all hold at once (G1h), the estimand becomes the minimizer of some weighted distance, with the weights determined by the estimator. Section 3 begins by turning this sentence into theorems.

## 3 Identification and Large-Sample Theory

This section follows the format of Chapters 5 and 6: it first uses Tamarind's numbers or figures to show what is to be proved, then states numbered assumptions and theorems, followed by proof sketches of two to six steps; each step names the assumptions it uses, and regularity conditions are only cited to their source. After each theorem, an "In words" paragraph answers four questions: under which assumptions it estimates which quantity, what it estimates when the specification is wrong, who sets the weights, and what Tamarind's numbers say. Assumptions and theorems are numbered within this chapter; the parentheses give the source in Wooldridge (2010) or in the original paper. Sections 3.1 to 3.4 are the general theory of GMM, Sections 3.5 and 3.6 discuss finite samples and weak moments, and Sections 3.7 and 3.8 cover MD.

### 3.1 GMM and Its Consistency

First, does the estimator get close to the truth? In G2's Monte Carlo, as the number of regions rises from 100 to 800, the standard deviation of the two-step GMM estimates falls from 0.401 to 0.161, and the coverage of the nominal 95% interval rises from 0.876 to 0.941. A Monte Carlo at a fixed sample size can show only what the sampling distribution looks like; consistency is a statement about the limit as the sample size goes to infinity, and the theorem below handles that step.

The starting point is a set of population moment conditions: at the true parameter, the expectation of some function is zero,

$$\mathbb{E}[g(w_i, \theta_0)] = 0,$$

where $g$ is an $r$-dimensional vector-valued function, and its sample analogue is $\bar g_n(\theta) = n^{-1}\sum_{i=1}^n g(w_i,\theta)$. When $r = k$, one can generally make $\bar g_n(\hat\theta) = 0$ hold exactly, the case of exact identification. When $r > k$, this is generally impossible, the case of overidentification, and the best one can do is to bring the moment conditions as close to zero as possible "overall." Closeness is measured with an $r \times r$ positive definite weight matrix $W$ ($W$ may be estimated from the data, written $\hat W$; as long as it converges in probability to a positive definite matrix, the limits below write it as $W$), and the estimator is the parameter that minimizes the weighted distance:

$$\hat\theta_{GMM} = \arg\min_\theta\ \bar g_n(\theta)'\, W\, \bar g_n(\theta).$$

::: {.assumption}
**GMM.1 (identification; conditions (e) and (f) of Wooldridge 2010 Theorem 14.1, p. 526; equation (14.1), p. 525)** $\mathbb{E}[g(w,\theta_0)] = 0$, and for all $\theta \in \Theta$ with $\theta \neq \theta_0$, $\mathbb{E}[g(w,\theta)] \neq 0$; $\hat W \xrightarrow{p} W$, with $W$ positive definite.

**GMM.2 (conditions (a) to (d) of Wooldridge 2010 Theorem 14.1, p. 526)** $\{w_i\}$ is a random sample and $\Theta$ is compact; for each $\theta$, $g(\cdot,\theta)$ is Borel measurable; for each $w$, $g(w,\cdot)$ is continuous on $\Theta$; there is a nonnegative function $b(w)$ with $\mathbb{E}[b(w)] < \infty$ such that the absolute value of every component of $g$ is, for all $\theta$, no greater than $b(w)$.
:::

::: {.intuition}
Think of GMM as weighting the testimony of ten witnesses. Each experiment gives an elasticity on its own: its effect on lessons divided by its effect on price. The ten pieces of testimony contradict one another, because each carries sampling noise. The identity weight listens to the witnesses according to how loud they are, so the experiments recorded in the largest units shout the loudest; the efficient weight listens according to credibility, so experiments with little noise and a large push on price count for more. If the witnesses are all describing the same thing, any way of listening hears the truth in the end, and the only difference is how long you must listen; if they are describing different things, the way you listen decides what you end up believing.
:::

GMM.1 has two halves. The first half is an economic assumption: in G2, random assignment plus the exclusion restriction. The second half is identification: only $\theta_0$ sets the population moments to zero. It requires $r \ge k$; for linear moment conditions like those of G2, it is equivalent to $\mathbb{E}[z_i(1, \log p_i)]$ having full column rank, that is, to the experiments really moving the price. In G3 this matrix is already close to rank deficient, a problem Section 3.6 takes up. GMM.2 is the moment-function counterpart of M.2 in Chapter 5.

::: {.theorem}
**Theorem 1 (consistency of GMM; Wooldridge 2010 Theorem 14.1, p. 526)** Under GMM.1 and GMM.2, $\bar g_n(\theta)'\hat W\bar g_n(\theta)$ has a minimizer $\hat\theta$, and $\hat\theta \xrightarrow{p} \theta_0$.
:::

The proof sketch has three steps. Step one (GMM.2): the uniform weak law of large numbers of Theorem 1 of Chapter 5, applied to each component of $g$, makes $\bar g_n(\theta)$ converge uniformly on $\Theta$ to $\mathbb{E}[g(w,\theta)]$; together with $\hat W \xrightarrow{p} W$, the sample criterion also converges uniformly to $Q_0(\theta) = \mathbb{E}[g(w,\theta)]'W\,\mathbb{E}[g(w,\theta)]$ (equations (14.3) and (14.4)). Step two (GMM.1): $W$ is positive definite, so $Q_0(\theta) \ge 0$, with equality only where $\mathbb{E}[g(w,\theta)] = 0$, that is, at $\theta = \theta_0$; $\theta_0$ is the unique valley floor, which is M.1 of Chapter 5. Step three (steps one and two): the proof of Theorem 2 of Chapter 5 uses only uniform convergence and a unique valley floor, never the fact that "the criterion is an average of single-observation losses," so it carries over unchanged.

In words: under GMM.1 and GMM.2, GMM estimates the parameter that makes all the moment conditions hold at once, which in G2 is the demand elasticity defined by random assignment and the exclusion restriction. When the specification is wrong, for example when some experiment violates the exclusion restriction, or when different experiments aim at different elasticities (the experiments move prices in different regions or price ranges, where the elasticities also differ), no $\theta$ can set all the moments to zero. If the elasticity merely varies at random across regions, unrelated both to experiment assignment and to promotions, the ten ratios all still equal the average elasticity, and the moment conditions do not conflict. When there is no solution, $\hat\theta$ converges to the minimizer of $Q_0$ (assumed unique), the pseudo-true value $\theta_*$; it depends on $W$, so each weight has its own pseudo-true value, and G1h in Section 3.7 puts numbers on this. Who sets the weights: the researcher chooses $W$. Under correct specification it does not affect where the estimator converges, because $\theta_0$ sets every moment to zero, and zero stays zero however it is weighted; what it affects is precision. In the same G2 sample, the identity weight, the 2SLS weight, and the two-step efficient weight give $-1.307$, $-1.843$, and $-2.012$; across 2,000 replications the three average $-1.541$, $-1.523$, and $-1.517$ (all slightly biased toward OLS, for reasons given in Section 3.5; where an estimator converges is a large-sample statement), yet their standard deviations are 0.610, 0.368, and 0.302.

### 3.2 Asymptotic Normality and the GMM Sandwich

The Monte Carlo standard deviation of the identity-weight estimate is 2.02 times that of the two-step efficient weight, although both use the same set of moment conditions. Where does the difference come from? Look at the first-order condition. The estimate $\hat\theta$ sets the gradient of the criterion to zero, $\nabla\bar g_n(\hat\theta)'\hat W\bar g_n(\hat\theta) = 0$; that is, $k$ linear combinations of the sample moments are pushed to zero. Expand $\bar g_n(\hat\theta)$ around $\theta_0$ as $\bar g_n(\theta_0) + D_0(\hat\theta - \theta_0)$ ($D_0$ is the Jacobian defined in the next paragraph) and solve for $\hat\theta - \theta_0$: the sampling error of the estimator is a fixed linear combination of the sampling errors of the $r$ sample moments, with coefficients determined jointly by $W$ and by how sensitive the moment conditions are to the parameters.

Write the Jacobian as $D_0 = \mathbb{E}[\partial g(w,\theta_0)/\partial\theta']$; it is an $r \times k$ matrix of derivatives, not a treatment variable. Write $S_0 = \mathbb{E}[g(w,\theta_0)\,g(w,\theta_0)']$; $S_0$ is the variance matrix of the moment conditions and has nothing to do with the survival function $S(t \mid x)$ of Chapter 6.

::: {.assumption}
**GMM.3 (conditions of Wooldridge 2010 Theorem 14.2, p. 527)** In addition to GMM.1 and GMM.2: (a) $\theta_0$ is in the interior of $\Theta$; (b) for each $w$, $g(w,\cdot)$ is continuously differentiable on the interior of $\Theta$; (c) every element of $g(w,\theta_0)$ has a finite second moment; (d) every element of $\nabla_\theta g(w,\theta)$ is dominated by some $b(w)$ with finite expectation; (e) $D_0$ has rank $k$.
:::

Of the five conditions in GMM.3, two are the ones to watch in applications. (a) rules out a truth on the boundary of the parameter space. (e) is the local version of identification: $D_0$ has full column rank, so near the truth the moment conditions respond in every direction of the parameter; in G3 $D_0$ is already close to rank deficient, so this condition holds in name but is effectively violated in finite samples, the problem of Section 3.6. (b) through (d) are technical conditions: linear moment conditions satisfy (b) automatically, and (d) requires only that the expected product of instruments and regressors be finite.

::: {.theorem}
**Theorem 2 (asymptotic normality of GMM; Wooldridge 2010 Theorem 14.2, equations (14.8) to (14.12), p. 527)** Under GMM.1 through GMM.3, let $\Lambda_0 = -(D_0'WD_0)^{-1}D_0'W$. Then

$$\begin{aligned}
\sqrt{n}(\hat\theta - \theta_0) &= \Lambda_0\,\sqrt{n}\,\bar g_n(\theta_0) + o_p(1) \\
&\xrightarrow{d} \mathcal{N}\big(0,\ \Lambda_0 S_0 \Lambda_0'\big),
\end{aligned}$$

where $\Lambda_0 S_0 \Lambda_0' = (D_0'WD_0)^{-1}D_0'WS_0WD_0(D_0'WD_0)^{-1}$.
:::

The proof sketch has four steps. Step one (GMM.3(a) and (b), Theorem 1): with probability approaching one $\hat\theta$ lies in the interior, so the first-order condition holds (equation (14.5)). Step two (GMM.3(b) and (d)): take a mean-value expansion of $\bar g_n(\hat\theta)$ around $\theta_0$; the sample Jacobian, both at $\hat\theta$ and at the mean values, converges in probability to $D_0$ (uniform convergence plus consistency). Step three (the positive definite $W$ of GMM.1, and GMM.3(e)): $D_0'WD_0$ is invertible, and solving gives the linear representation above (equations (14.8) to (14.10)). Step four (GMM.1 and GMM.3(c)): $\mathbb{E}[g(w,\theta_0)] = 0$ with finite variance, so the central limit theorem gives $\sqrt n\,\bar g_n(\theta_0) \xrightarrow{d} \mathcal{N}(0, S_0)$, and multiplying by $\Lambda_0$ completes the proof.

To estimate this variance, use sample counterparts: $\hat D = n^{-1}\sum_i \nabla_\theta g(w_i,\hat\theta)$, and $\hat S_n(\hat\theta) = n^{-1}\sum_i g(w_i,\hat\theta)\,g(w_i,\hat\theta)'$ in place of $S_0$; substituting them gives the sandwich variance estimator (Wooldridge 2010 equations (14.13) to (14.15)), which is how the standard errors of the identity weight and of 2SLS in Section 6 are computed. This is the sandwich of Section 3.4 of Chapter 5 recast in the language of moment conditions: the bread is $D_0'WD_0$ and the filling is $D_0'WS_0WD_0$ (Wooldridge 2010 equations (14.9) and (14.11)). $D_0$ plays the role of curvature: the more sensitive the moment conditions are to the parameters, the steeper the valley and the smaller the variance, which is the moment-condition version of the strength of identification in Section 3.2 of Chapter 5. Andrews, Gentzkow and Shapiro (2017) call $\Lambda_0$ the sensitivity, and Section 4.2 uses it to answer "which moments identify the parameter."

In words: under GMM.1 through GMM.3, GMM estimates $\theta_0$, and its sampling error is the linear combination $\Lambda_0\bar g_n(\theta_0)$ of the sampling errors of the individual sample moments. When the specification is wrong, $\hat\theta$ converges to the pseudo-true value (assuming the minimizer of $Q_0$ is unique), but there the mean of $\bar g_n(\theta_*)$ is not zero, the expansion acquires extra terms related to the estimation error of $\hat W$, and the variance formula above is no longer guaranteed to hold; conventional and robust standard errors alike are reliable only when the moment conditions hold. Who sets the weights: each column of $\Lambda_0$ tells you the coefficient with which the noise of one moment enters $\hat\theta$, and these coefficients are set jointly by the $W$ the researcher chooses and the $D_0$ the design determines (how hard each experiment pushes the price). G2's identity weight measures each moment in its recorded units: experiments measured in dollars with large values, such as the plus-or-minus USD 10 fee-waiver threshold, weigh most in the criterion, although they do not necessarily push the price hardest. So the standard deviation of the identity-weight estimate is 2.02 times that of the efficient weight, even though both aim at $-1.6$.

### 3.3 The Efficient Weight and Two-Step GMM

Why does the identity weight lose out? Compute the asymptotic standard errors in the population at $n = 200$: 0.636 for the identity weight, 0.374 for the 2SLS weight, and 0.329 for the efficient weight. The best weight exists, and it has a closed form.

::: {.theorem}
**Theorem 3 (the efficient weight and two-step GMM; Wooldridge 2010 equations (14.16) to (14.18), p. 528, Theorem 14.3, p. 540; Hansen 2014 Theorem 13.3.2)** Let GMM.1 through GMM.3 hold, let $S_0$ be positive definite, and suppose that $\theta_0$ has a neighborhood $\mathcal{U}$ on which $\mathbb{E}[\sup_{\theta \in \mathcal{U}}\lVert g(w,\theta)\rVert^2] < \infty$. (i) For any positive definite $W$, the variance of Theorem 2 minus $(D_0'S_0^{-1}D_0)^{-1}$ is positive semidefinite, with equality when $W = S_0^{-1}$. (ii) Two-step GMM: the first step uses any weight $W^{(1)}$ satisfying GMM.1 to obtain $\hat\theta^{(1)}$ and computes $\hat S = \hat S_n(\hat\theta^{(1)})$, where $\hat S_n(\theta) = n^{-1}\sum_{i=1}^n g(w_i,\theta)\,g(w_i,\theta)'$; the second step uses $\hat S^{-1}$ as the weight to obtain $\hat\theta^{(2)}$. Then

$$\mathrm{Avar}\big(\sqrt n(\hat\theta^{(2)} - \theta_0)\big) = (D_0'S_0^{-1}D_0)^{-1},$$

and the conventional variance estimator is $(\hat D'\hat S^{-1}\hat D)^{-1}/n$, where $\hat D$ is the sample Jacobian.
:::

The proof sketch has three steps. Step one (Theorem 2): substitute $W = S_0^{-1}$; the filling $D_0'S_0^{-1}S_0S_0^{-1}D_0$ equals the bread $D_0'S_0^{-1}D_0$ and the sandwich collapses, by the same mechanism through which the information matrix equality of Section 3.3 of Chapter 6 collapses the sandwich. Step two (Wooldridge 2010 Theorem 14.3): GMM with weight $W$ amounts to an estimator whose score is $D_0'Wg$ and whose "bread" is $D_0'WD_0$; it and the estimator with $W = S_0^{-1}$ satisfy $\mathbb{E}[(D_0'Wg)(D_0'S_0^{-1}g)'] = D_0'WD_0$ for all $W$, which is the condition of equation (14.51) in Theorem 14.3 (with the proportionality constant equal to 1); the argument is of a piece with the proof of the efficiency of MLE in Theorem 4 of Chapter 6, and (i) follows. Step three (Theorem 1, GMM.2, and the moment condition on $\mathcal{U}$ in the theorem, the standard regularity condition needed for $\hat S_n$ to be consistent): $\hat\theta^{(1)}$ is consistent, and uniform convergence gives $\hat S_n(\hat\theta^{(1)}) \xrightarrow{p} S_0$, so $\hat S^{-1}$ meets GMM.1's requirement on the weight, and Theorem 2 with $W = S_0^{-1}$ gives (ii).

In words: under GMM.1 through GMM.3, two-step GMM estimates the same $\theta_0$ as any other weight; it merely has the smallest variance among all GMM estimators that use these $r$ moment conditions. "Efficient" is relative to the given moment conditions: under the efficient weight, adding further valid moment conditions does not make the asymptotic variance larger (Wooldridge 2010 p. 530); MLE amounts to GMM with the score as the moment conditions (the logit of Chapter 12 is an example: its score is the set of moments "residual times regressor"), and when the distribution is correctly specified the score is the best moment condition (Theorem 4 of Chapter 6); when the distribution is misspecified, MLE becomes the QMLE of Section 3.6 of Chapter 6, whereas GMM, which relies only on the moment conditions, remains consistent as long as the moment conditions hold. When the specification is wrong, the efficient weight puts the largest weight on the most precise moments, and the most precise moment may be the very one that is wrong, so the pseudo-true value moves with $S_0^{-1}$. Who sets the weights: the noise structure of the data, $S_0$, not the researcher; noisy moments are turned down. G2's deep-discount experiment puts only 20% of regions in the arm with the deeper discount and lower price, and the demand noise in those regions is four times that elsewhere: this moment carries more noise per unit of instrument variance than the other moments, and the efficient weight turns it down; the 2SLS weight $(n^{-1}\sum_i z_iz_i')^{-1}$ presumes that the noise of each moment is proportional to its instrument variance, which holds exactly under homoskedasticity. The ratio of the two asymptotic standard errors is 0.878. If all ten experiments used balanced 50/50 assignment with the same noise rule, both asymptotic standard errors would be 0.4062. Then $S_0$ is still not proportional to $\mathbb{E}[z_iz_i']$ (the intercept moment is correlated with the deep-discount moment), but the noise per unit of instrument variance is the same for every experiment moment, and the two weights give the elasticity the same variance; so what matters is whether this ratio is the same across moments (it is in the balanced design but not in G2), not whether $S_0$ is proportional to $\mathbb{E}[z_iz_i']$. The balanced design has larger standard errors mainly because half the regions fall into the noisy arm.

The "efficiency" of Theorem 3 carries one more qualification: the comparison is confined to the given $r$ moments. What random assignment actually delivers is the stronger conditional moment restriction $\mathbb{E}[\varepsilon_i \mid z_i] = 0$, so any function of the experiment assignments is a valid instrument, and far more than eleven moment conditions could be added. Wooldridge (2010) §14.4.3 uses Theorem 14.3 to give the optimal choice within this class (equations (14.54) to (14.57), p. 542): the optimal instrument is the conditional expectation, given $z_i$, of the derivative of the residual with respect to the parameters, divided by the conditional variance of the residual, and the resulting variance is the efficiency bound under the conditional moment restriction (equation (14.60), pp. 543-544). In G2 the variance of the demand shock depends on the deep-discount assignment, so the optimal instrument contains interactions of each experiment with the deep-discount assignment, and two-step GMM with only the eleven linear moments does not reach this bound. Adding these interactions raises asymptotic efficiency at the cost of more moment conditions, which aggravates the finite-sample problems of Section 3.5.

### 3.4 The J Test

With eleven moment conditions and two parameters, can the nine surplus ones be used to test the model? In the G2 sample, the objective function at the two-step GMM estimate, times $n$, is 9.81. Whether this number is large or small depends on its distribution when "all the moment conditions hold."

::: {.theorem}
**Theorem 4 (the J test; Hansen 1982; Wooldridge 2010 equation (14.19), p. 528; Hansen 2014 Theorem 13.6.1, p. 283)** Under the conditions of Theorem 3, if $r > k$,

$$J = n\,\bar g_n(\hat\theta^{(2)})'\,\hat S^{-1}\,\bar g_n(\hat\theta^{(2)}) \xrightarrow{d} \chi^2_{r-k}.$$

If no $\theta$ gives $\mathbb{E}[g(w,\theta)] = 0$ and $\hat S$ converges to a positive definite matrix, then $J$ goes to infinity in probability (this sentence follows from step three of the proof; the sources cited state only the limiting distribution under the null hypothesis).
:::

The proof sketch has three steps. Step one (Theorems 2 and 3): expand $\bar g_n(\hat\theta^{(2)})$ around $\theta_0$ and substitute the linear representation, which gives $\sqrt n\,\bar g_n(\hat\theta^{(2)}) = [I_r - D_0(D_0'S_0^{-1}D_0)^{-1}D_0'S_0^{-1}]\sqrt n\,\bar g_n(\theta_0) + o_p(1)$. Step two (GMM.3(c)): premultiply by $S_0^{-1/2}$; the right side becomes $I_r - P$ times $S_0^{-1/2}\sqrt n\,\bar g_n(\theta_0)$, which converges in distribution to $\mathcal{N}(0, I_r)$, where $P = S_0^{-1/2}D_0(D_0'S_0^{-1}D_0)^{-1}D_0'S_0^{-1/2}$ is a projection matrix of rank $k$. Step three (the $\hat S \xrightarrow{p} S_0$ of Theorem 3): $J$ is asymptotically equal to the quadratic form of this standard normal vector with kernel $I_r - P$, and $I_r - P$ is symmetric and idempotent with rank $r - k$, so the limit is $\chi^2_{r-k}$. If the moment conditions cannot all hold at once, $\bar g_n(\hat\theta^{(2)})$ converges to a nonzero vector, and $J$ grows at the rate $n$.

In words: J estimates no parameter; it tests the claim that "there exists a $\theta$ that makes the $r$ moment conditions hold at once." Estimating $\theta$ uses up $k$ directions, and the remaining $r - k$ directions are used for the checkup. Under GMM.1 through GMM.3 it follows $\chi^2_{r-k}$; when the specification is wrong it will reject, but it will not tell you which moment is wrong. The cause may be that some experiment violates the exclusion restriction, or that different experiments aim at different elasticities, each right in its own sense; G1h in Section 3.7 is a case of the latter. Who sets the weights: $J$ weights by $\hat S^{-1}$, so the more precise a moment, the more its departure from zero matters; its $\chi^2$ distribution holds only under the efficient weight, and the objective function value under the identity weight does not have this distribution. In the G2 sample, $J = 9.81$ with 9 degrees of freedom, $p = 0.366$. Across 2,000 replications, its rejection rate at the 5% level is 0.0385, below the nominal level by more than two Monte Carlo standard errors (0.0049), which shows that at $n = 200$ the $\chi^2$ approximation is not yet accurate enough.

Restrictions on $\theta$ can also be tested with the criterion function. Keeping the same weight $\hat S^{-1}$, minimize with and without the restrictions; the difference between the two criterion values, times $n$, follows a $\chi^2$ when the restrictions hold, with degrees of freedom equal to the number of restrictions (Wooldridge 2010 equation (14.20), p. 529). It is the GMM version of the LR test and, like it, does not change with how the restrictions are written (Section 3.5 of Chapter 6), a property the Wald statistic lacks. The J of MD in Section 3.7 is just such a difference: the unrestricted reduced form fits perfectly, with a criterion value of zero, so the criterion value after imposing "the 16 elasticities are equal" is the test statistic itself.

### 3.5 Finite Samples: Two-Step Standard Errors, Iterated GMM, CUE, and Many-Moment Bias

Theorem 3 says that the asymptotic standard error of two-step GMM at $n = 200$ is 0.329, but that is not the number to compare with. Across 2,000 replications, the true standard deviation of the two-step estimate is 0.302, while the conventional standard error averages only 0.278, and the nominal 95% interval has coverage of only 0.9115; the standard error corrected as in Windmeijer (2005) averages 0.303, with coverage of 0.935. Push the number of regions from 100 to 800 (see the figure in this section): at $n = 100$ the conventional interval covers 0.876 and the corrected one 0.928; at $n = 800$ the two are 0.941 and 0.945, and the gap shrinks as the sample grows.

Why is the conventional standard error too small? It treats $\hat S$ as known, but $\hat S$, computed at the first-step estimate $\hat\theta^{(1)}$, is an $11 \times 11$ matrix estimated from 200 regions, and $\hat\theta^{(2)}$ depends on $\hat\theta^{(1)}$ through the weight. This dependence vanishes in first-order asymptotic theory, for the same reason that the first step does not affect the variance when $\Gamma_0 = 0$ in Section 3.6 of Chapter 5: changes in the weight multiply $\bar g_n$, which has mean zero at the truth. In finite samples it does not vanish.

::: {.assumption}
**LIN.1 (linear moment conditions; Windmeijer 2005)** $g(w,\theta) = z\,(y - x'\theta)$, where $z$ is an $r \times 1$ vector of instruments and $x$ a $k \times 1$ vector of regressors.
:::

::: {.theorem}
**Theorem 5 (finite-sample variance correction for two-step GMM; Windmeijer 2005)** Under LIN.1 and the conditions of Theorem 3, let $\hat\theta(\bar\theta)$ denote the GMM estimate with weight $\hat S_n(\bar\theta)^{-1}$, so that $\hat\theta^{(2)} = \hat\theta(\hat\theta^{(1)})$, and

$$\hat\theta^{(2)} - \theta_0 = \big[\hat\theta(\theta_0) - \theta_0\big] + D_W\,(\hat\theta^{(1)} - \theta_0) + o_p(n^{-1}),$$

where $D_W = \partial\hat\theta(\bar\theta)/\partial\bar\theta'$, evaluated at $\theta_0$, is $O_p(n^{-1/2})$, so the second term is $O_p(n^{-1})$. The corrected variance based on this is

$$\hat V_c = \hat V_2 + \hat D_W\hat V_2 + \hat V_2\hat D_W' + \hat D_W\hat V_1\hat D_W',$$

where $\hat V_2 = (\hat D'\hat S^{-1}\hat D)^{-1}/n$ is the conventional variance, $\hat V_1$ is the sandwich variance of the first step, and $\hat D_W$ has column $j$ equal to $(\hat D'\hat S^{-1}\hat D)^{-1}\hat D'\hat S^{-1}\big(\partial\hat S_n(\theta)/\partial\theta_j\big)\hat S^{-1}\bar g_n(\hat\theta^{(2)})$, with the derivative evaluated at $\hat\theta^{(1)}$.
:::

The proof sketch has three steps. Step one (LIN.1): with linear moment conditions, $\bar g_n(\theta) = \bar g_n(\theta_0) + \hat D(\theta - \theta_0)$, where $\hat D = -n^{-1}\sum_i z_ix_i'$ does not depend on $\theta$, so $\hat\theta(\bar\theta) = \theta_0 - (\hat D'\hat S_n(\bar\theta)^{-1}\hat D)^{-1}\hat D'\hat S_n(\bar\theta)^{-1}\bar g_n(\theta_0)$ holds exactly. Step two (LIN.1, Theorem 1): expand this in the weight's argument $\bar\theta$ around $\theta_0$; the derivative falls on $\hat S_n(\bar\theta)^{-1}$ and is multiplied by $\bar g_n$, so $D_W = O_p(n^{-1/2})$, and the second term is one order smaller than the first. Step three (Theorem 3, and the discussion of Lemma 14.1 in Wooldridge 2010 p. 539): take the variance of the sum of the two terms; the asymptotic covariance between the efficient estimator and any other consistent GMM estimator equals the efficient estimator's own variance, so the cross terms are approximated by $\hat V_2$, and evaluating each term at the sample estimates gives $\hat V_c$.

In words: under LIN.1 and the conditions of Theorem 3, the Windmeijer correction changes neither the quantity estimated nor the estimate; it only builds the dependence of $\hat\theta^{(2)}$ on the first step into the variance. When the moment conditions are nonlinear in the parameters there is no corresponding theorem: whether the correction helps also depends on other remainder terms of the same order, so it cannot simply be carried over; and when the moment conditions themselves are wrong, correcting the variance cannot rescue the point estimate. Who sets the weights: the size of $\hat D_W$ depends on how sensitive the weight matrix is to $\bar\theta$, and the stronger the heteroskedasticity, the more moment conditions there are, and the smaller the sample, the larger the correction. In G2 the noise of the deep-discount arm is four times that elsewhere, so the weight is very sensitive to the first step: in the sample the corrected standard error is 0.418, 1.167 times the conventional 0.358; in the Monte Carlo the average corrected standard error, 0.303, is almost the same as the true standard deviation, 0.302.

Two steps are not the only option. Iterated GMM recomputes the weight at each new estimate and re-estimates until convergence: $\hat\theta^{(j+1)} = \arg\min_\theta \bar g_n(\theta)'\hat S_n(\hat\theta^{(j)})^{-1}\bar g_n(\theta)$; continuously updated GMM (CUE from here on) instead lets the weight move together with $\theta$:

$$\hat\theta_{CUE} = \arg\min_\theta\ \bar g_n(\theta)'\,\hat S_n(\theta)^{-1}\,\bar g_n(\theta).$$

CUE was proposed by Hansen, Heaton and Yaron (1996), who also compared the finite-sample performance of the three estimators (two-step, iterated, and CUE); Hansen (2014) §13.4 and §13.5 give brief definitions. Under strong identification the three have the same asymptotic distribution as in Theorem 3; CUE does not depend on a first-step weight, and it is invariant to a change in the units of the moment conditions. In the G2 sample, iterated GMM converges after 13 iterations to $-2.043$, and CUE gives $-2.304$. Across 2,000 replications, the rejection rates of J at the 5% level for two-step GMM, iterated GMM, and CUE are 0.0385, 0.0295, and 0.025, all below the nominal level.

![Block G2's Monte Carlo: coverage of the elasticity by nominal 95% intervals as a function of the number of regions $n$ (log scale; 100, 200, 400, 800), one line for each of six estimators: 2SLS (robust standard errors), identity weight (robust standard errors), two-step GMM (conventional and Windmeijer standard errors), and iterated GMM and CUE (conventional standard errors). The dark gray horizontal solid line is 0.95, and the two dotted lines above and below it are the 95% Monte Carlo error band for 1,000 replications.](assets/fig/fig_07_gmm_coverage.svg)

When there are many moment conditions, two-step GMM is also biased. Newey and Windmeijer (2009) point out that the expectation of the first-order condition of two-step GMM at the truth is not zero: the Jacobian of each observation is correlated with its own moments (as happens when the regressors are correlated with the error), each moment contributes one such term, and the bias, which points toward OLS, grows with the number of moment conditions, the strength of endogeneity, and the weakness of identification. CUE's weight moves with $\theta$, and its first-order condition removes the part of the Jacobian that is correlated with the moments, so its bias is much smaller, but its conventional standard errors are too small; Newey and Windmeijer give a larger asymptotic variance and the corresponding sandwich estimator under an asymptotic framework of "many weak moments." The two corrections target different things: Windmeijer corrects for the fluctuation of the weight matrix, Newey and Windmeijer for the fluctuation of the derivatives of the moment conditions. G2 has only 11 moments, and the bias is already visible: across 2,000 replications the two-step GMM estimates average $-1.517$, a bias of 0.083 toward OLS's $-0.382$; iterated GMM and 2SLS give $-1.514$ and $-1.523$. CUE's finite-sample moments need not exist, and individual samples can land far from the truth (Hansen, Heaton and Yaron 1996), so bias is better compared through the median: CUE's median is $-1.578$ and two-step GMM's is $-1.505$, which puts CUE much closer to the truth. But CUE is more dispersed: its Monte Carlo standard deviation, 0.321, is larger than two-step's, while its conventional standard error averages only 0.280, and its coverage is 0.9165. In the dynamic panels of Chapter 11, the number of moment conditions grows with the square of the number of periods, and this bias becomes a central problem.

### 3.6 Weak Moments and the S Confidence Set

G3 pushes the problem of Section 3.5 to the extreme. The sample first-stage F statistic is 3.61; two-step GMM gives $-0.229$ (standard error 0.640), and the 95% interval $[-1.482, 1.025]$ does not cover $-1.6$; neither do the Windmeijer-corrected interval $[-1.485, 1.027]$ and the 2SLS interval $[-1.483, 1.049]$. This one sample alone says little: even with correct coverage, a 95% interval has about a one-in-twenty chance of missing the truth. What is informative is repeated sampling: across 2,000 replications, the coverage rates of the three intervals are only 0.774, 0.7745, and 0.772; the median of the two-step estimates is $-0.854$, pulled toward the probability limit of OLS, 1.044. The experiments' push on price is much smaller than the noise of the promotions, so the denominator of the ratio "effect on lessons divided by effect on price" is almost all noise; the estimator is no longer approximately normal, and the standard errors lose their meaning.

::: {.intuition}
The trouble with weak moments lies in the denominator of the ratio. Each experiment's elasticity is "the change in lessons divided by the change in price"; price moves by at most 0.19%, so the denominator is almost all noise, and whenever it comes close to zero the ratio lands far from the truth. Such a ratio has very fat tails, and its center is pulled toward OLS by the correlation that the promotions create, so the normal approximation fails fundamentally. The S statistic sidesteps the division: given a candidate elasticity, it asks only whether the demand shocks computed with this elasticity are uncorrelated with experiment assignment, a question that involves only addition and subtraction, and its $\chi^2$ approximation is reliable whatever the strength of identification.
:::

Stock and Wright (2000) describe this situation with an asymptotic framework in which the population mean of the moment conditions varies with $\theta$ by an amount of the same order as the sampling noise in the sample moments. Then $\hat\theta$ is inconsistent, and Theorem 2 does not apply. How weak is too weak is a question for Chapter 9, which covers rules of thumb for the first-stage F statistic; in G3's 2,000 replications, 99.9% of the samples have a first-stage F below 10. The way out is not to estimate $\theta$ but to test each candidate value of it.

::: {.assumption}
**SW.1 (conditions of Stock and Wright 2000 Theorem 2, p. 1063)** $\mathbb{E}[g(w,\theta_0)] = 0$; $\sqrt n\,\bar g_n(\theta_0) \xrightarrow{d} \mathcal{N}(0, S_0)$, with $S_0$ positive definite; $\hat S_n(\theta_0) \xrightarrow{p} S_0$. No assumption at all is made about $D_0$.
:::

::: {.theorem}
**Theorem 6 (the S statistic and the S confidence set; Stock and Wright 2000 Theorem 2, p. 1063; for the S confidence set, see p. 1064)** Define the S statistic

$$\mathcal{S}_n(\theta) = n\,\bar g_n(\theta)'\,\hat S_n(\theta)^{-1}\,\bar g_n(\theta).$$

Under SW.1, whatever the strength of identification, $\mathcal{S}_n(\theta_0) \xrightarrow{d} \chi^2_r$; hence the S confidence set $\{\theta : \mathcal{S}_n(\theta) \le \chi^2_{r,1-\alpha}\}$ covers $\theta_0$ with asymptotic probability $1 - \alpha$.
:::

The proof sketch has three steps. Step one (SW.1): at $\theta_0$ no parameter is estimated, and the central limit theorem directly gives the normal limit of $\sqrt n\,\bar g_n(\theta_0)$. Step two (SW.1): $\hat S_n(\theta_0) \xrightarrow{p} S_0$. Step three (steps one and two): a quadratic form in an asymptotically normal vector, weighted by the inverse of its own variance, has the limit $\chi^2_r$; $\theta_0$ lies in the set if and only if $\mathcal{S}_n(\theta_0)$ does not exceed the critical value, and the coverage follows. $D_0$ never appears in the argument, which is why the result is robust to weak identification.

In words: under SW.1, the S confidence set does not estimate a point; it collects every $\theta$ compatible with the moment conditions, and when identification is weak it widens, even to the point of being unbounded. The price is that it tests the $r$ moments jointly, so under overidentification it is wider than the Wald interval even when identification is strong: G2's S confidence set $[-4.09, -0.995]$ (the intercept moment exactly identifies $\beta_0$; solving that moment for $\beta_0$ and substituting back, that is, partialling out the intercept, leaves 10 experiment moments, so there are 10 degrees of freedom, not $r = 11$) is much wider than the two-step GMM Wald interval $[-2.713, -1.311]$. When the specification is wrong, no $\theta$ makes the moment conditions hold, the S statistic is too large everywhere, and the set shrinks and may even be empty. An empty set is a signal of misspecification: under correct specification, its asymptotic probability is at most $\alpha$. Weak identification can only widen the set, even to the point of covering the entire parameter space; it cannot make the set empty. What is truly hard to read is a set that is very small but not empty: it may come from correct specification with precise estimation, or from misspecification with too little evidence to reject outright (Stock and Wright 2000, pp. 1064-1065). Who sets the weights: at every candidate $\theta$, $\hat S_n(\theta)$ is recomputed; this is CUE's weight, and CUE is the minimizer of the S statistic. G2's S curve reaches its minimum near $-2.3$, close to but not equal to CUE's $-2.304$: the curve uses the 10 moments left after the intercept is partialled out, whereas CUE estimates $\beta_0$ and $\theta$ jointly from all 11 moments. In the G3 sample, $\mathcal{S}_n(-1.6) = 3.06$ and $p = 0.547$, and the S confidence set is $[-7.04, 1.29]$, which covers the truth; across 2,000 replications its coverage is 0.948. CUE's median is $-1.582$, very close to the truth, but the coverage of its Wald interval is only 0.890. The J test fails here too: its $\chi^2_{r-k}$ approximation requires a consistent $\hat\theta$; under weak moments the first-step estimate is far from the truth, pulled toward OLS, and the residual variance computed there is too small, so $\hat S$ underestimates $S_0$, J is inflated, and the nominal 5% test has a rejection rate of 0.12. The Anderson-Rubin test of Chapter 9 is the special case of the S statistic in the linear instrumental variables model.

![The S statistic $\mathcal{S}_n(\theta)$ of blocks G2 and G3 as a function of the elasticity $\theta$ (solid line); the dashed line is the 95% critical value of the $\chi^2$ (10 degrees of freedom for G2, whose intercept is partialled out, and 4 for G3), the light band is the S confidence set, the error bar below is the two-step GMM 95% Wald interval with its point estimate, and the vertical dotted line is the truth, $-1.6$. The left panel, G2, has ten experiments and a first-stage F of about 54; the right panel, G3, has only four perturbations of a few cents and a first-stage F of 3.6.](assets/fig/fig_07_s_set.svg)

Across 2,000 replications, 36.85% of the S confidence sets are unbounded. This is not a defect of the method but an honest answer. As $|\theta|$ goes to infinity, $\mathcal{S}_n(\theta)$ tends to the S statistic computed from the "experiments move the price" moments alone, that is, to a test of the hypothesis that the experiments do not affect the price at all. When this statistic is below the critical value, the data cannot rule out that the experiments have no effect on price whatsoever; then every elasticity large enough in absolute value is compatible with the data, and the set must be unbounded. Reporting a finite Wald interval in such a sample amounts to claiming knowledge that the data do not contain. In this G3 sample, the limiting value is 14.01, above the critical value of 9.49, so the set is bounded.

### 3.7 Classical Minimum Distance

Return to G1. The price slopes of the 16 cities range from $-2.599$ (Nashville) to $-1.096$ (Columbus), and every city's 95% interval covers $-1.6$. The simple average is $-1.528$ (standard error 0.117); the precision-weighted average, which is the MD estimate developed below, is $-1.535$ (standard error 0.0645), with interval $[-1.662, -1.409]$. Both cover the truth, and the standard error of the former is 1.81 times that of the latter. The platform needs only the slope and standard error that each city reports, not household-level data.

::: {.assumption}
**MD.1 (reduced-form estimates; Wooldridge 2010 equation (14.63), p. 545)** $\sqrt n(\hat\pi - \pi_0) \xrightarrow{d} \mathcal{N}(0, \Sigma_\pi)$, with $\Sigma_\pi$ positive definite; $\hat\Sigma_\pi \xrightarrow{p} \Sigma_\pi$.

**MD.2 (structural restrictions; Wooldridge 2010 §14.5, pp. 545-546)** $\pi_0 = h(\theta_0)$, where $h: \mathbb{R}^k \to \mathbb{R}^r$ is known and continuously differentiable; on the compact set $\Theta$, $\theta_0$ is the unique solution of $h(\theta) = \pi_0$ and lies in the interior; $D_0 = -\partial h(\theta_0)/\partial\theta'$ has rank $k$.
:::

This $D_0$ is the Jacobian of Section 3.2 when MD is viewed as GMM: the "sample moments" of MD are $\hat\pi - h(\theta)$, whose derivative with respect to $\theta'$ is $-\partial h/\partial\theta'$. In G1, $h(\theta) = \theta\,\boldsymbol\iota_{16}$ and $D_0 = -\boldsymbol\iota_{16}$, so the rank condition holds automatically.

::: {.theorem}
**Theorem 7 (optimal MD; Wooldridge 2010 §14.5, equations (14.64) to (14.69), pp. 545-546; for the version with linear restrictions, see Hansen 2014 Theorem 7.6.1 and Theorem 8.12.1)** Under MD.1 and MD.2, let $\hat\theta_{MD}$ minimize $(\hat\pi - h(\theta))'\hat\Sigma_\pi^{-1}(\hat\pi - h(\theta))$. (i) $\hat\theta_{MD} \xrightarrow{p} \theta_0$ and $\sqrt n(\hat\theta_{MD} - \theta_0) \xrightarrow{d} \mathcal{N}\big(0, (D_0'\Sigma_\pi^{-1}D_0)^{-1}\big)$, the smallest asymptotic variance over all weight matrices; (ii) the overidentification statistic

$$J = n\,\big(\hat\pi - h(\hat\theta_{MD})\big)'\,\hat\Sigma_\pi^{-1}\,\big(\hat\pi - h(\hat\theta_{MD})\big) \xrightarrow{d} \chi^2_{r-k}.$$
:::

The proof sketch has three steps. Step one (MD.1 and MD.2): treat $\hat\pi - h(\theta)$ as the sample moments; it converges uniformly to $\pi_0 - h(\theta)$ ($h$ contains no randomness), which is zero only at $\theta_0$, and Theorem 1 gives consistency. Step two (MD.2): the first-order condition (equation (14.65)) and a mean-value expansion (equation (14.66)) give $\sqrt n(\hat\theta_{MD} - \theta_0) = \Lambda_0\sqrt n(\hat\pi - \pi_0) + o_p(1)$, which is Theorem 2 with $S_0$ replaced by $\Sigma_\pi$; Theorem 3 says that $W = \Sigma_\pi^{-1}$ is optimal, which is equation (14.67). Step three (MD.1): the projection argument of Theorem 4 applies unchanged, which is equation (14.69).

In G1 the reduced-form estimates are mutually independent, $\hat\Sigma_\pi/n = \mathrm{diag}(\hat\sigma_1^2, \ldots, \hat\sigma_{16}^2)$, and $h(\theta) = \theta\,\boldsymbol\iota_{16}$, so Theorem 7 reduces to two familiar formulas:

$$\hat\theta_{MD} = \frac{\sum_{m=1}^{r} \hat\pi_m/\hat\sigma_m^2}{\sum_{m=1}^{r} 1/\hat\sigma_m^2}, \qquad J = \sum_{m=1}^{r} \frac{(\hat\pi_m - \hat\theta_{MD})^2}{\hat\sigma_m^2},$$

with standard error $\big(\sum_m 1/\hat\sigma_m^2\big)^{-1/2}$. The former is the inverse-variance weighted average, that is, the pooled estimate of a fixed-effect meta-analysis, and the latter is its heterogeneity test, Cochran's Q. The restriction $\pi_1 = \cdots = \pi_{16}$ is linear, so here J also equals the Wald statistic for "the 16 slopes are equal" (Hansen 2014 equation (8.17)).

In words: under MD.1 and MD.2, MD estimates the structural parameter that reproduces the reduced form, which in G1 is the common elasticity: the MD estimate is $-1.535$, $J = 14.64$ (15 degrees of freedom, $p = 0.478$), and across 2,000 replications the rejection rate of J at the 5% level is 0.046, while the MD interval has coverage of 0.9435. When the specification is wrong, no $\theta$ gives $h(\theta) = \pi_0$, and MD converges to the minimizer of $(\pi_0 - h(\theta))'\Sigma_\pi^{-1}(\pi_0 - h(\theta))$. In G1h this is $\theta_* = \sum_m \omega_m\theta_m$ with $\omega_m \propto 1/\sigma_m^2$, a precision-weighted average: $\theta_* = -1.395$, whereas the simple average is $-1.6$ and the household-weighted average is $-1.619$. The MD estimate is $-1.406$, and its interval $[-1.532, -1.280]$ covers $\theta_*$ and excludes the other two; the equal-weight estimate is $-1.538$, and its interval $[-1.767, -1.308]$ covers the simple average. Who sets the weights: the experimental design, not the importance of the cities. $\sigma_m^2$ depends on the number of participating households, the size of the price difference, and the noise: Chicago has the most active households and also the most test participants, but its price difference is only 6%, and its weight is 0.0714; Seattle's price difference is 12%, and its weight is 0.1735; Columbus has only 2,000 households and a price difference of 3%, and its weight is 0.0026, so the largest weight is about 67 times the smallest. High-income, less elastic cities such as Seattle and Boston weigh the most, so the combined elasticity is lower in absolute value than that of any kind of "average city." J is now 74.85, with $p = 6.0 \times 10^{-10}$; computed from the DGP, its probability of rejecting at the 1% level is 0.998. This rejection does not show that any test was faulty: each city's slope is right for its own city, and what is wrong is the restriction "there is only one elasticity."

::: {.intuition}
Inverse-variance weighting has a consequence that is easy to overlook. Suppose the platform dares to run only small 3% tests in big cities but boldly tests 12% in small cities; when the tests are combined by precision, the small cities' tests are more precise and so weigh more. Under homogeneity this is no problem at all: every test estimates the same number, and you listen most to whoever is most precise. Under heterogeneity, the pooled elasticity becomes the elasticity of "the cities whose tests were run best," and how well a test was run depends on the original experimental design, not on how important the city is to the platform. When J rejects homogeneity, this is what it is reminding you of: you must first decide whose elasticity you want, and then choose the weights.
:::

![Blocks G1 (left) and G1h (right): the reduced-form price slopes $\hat\pi_m$ of the 16 cities with 95% intervals; point size is the optimal MD weight, and the cities are ordered from top to bottom by their G1 weights, largest first. The blue vertical line and light blue band are the MD estimate (the inverse-variance weighted average) and its 95% interval; orange crosses are the truth (the common elasticity $-1.6$ for G1, each city's $\theta_m$ for G1h); the panel titles give the J statistic.](assets/fig/fig_07_md_markets.svg)

### 3.8 The Asymptotic Equivalence of MD and GMM, and Chamberlain's $\pi$ Matrix

G1 can also be estimated without going through the city slopes, by GMM directly on the 65,000 households: two moments per city (one for the intercept and one for the price variation), 32 in all; the parameters are 16 city intercepts plus one common elasticity, 17 in all, leaving 15 overidentifying degrees of freedom, the same as MD. The reduced form is each city's intercept and slope, 32 parameters in all, and the within-city OLS moments are linear in them; the intercepts are also free city by city in the structure and MD absorbs them exactly, so it is enough to compare the 16 slopes. Two-step GMM gives $-1.535155$ and MD $-1.535176$, a difference of $-2.0 \times 10^{-5}$, about 0.0003 times the standard error; the two J statistics are 14.627 and 14.636. If the GMM weight is computed from the residuals of the city-by-city regressions, the differences in both the estimate and J fall below $10^{-8}$, leaving only rounding error.

::: {.theorem}
**Theorem 8 (the asymptotic equivalence of MD and GMM; the general statement and proof are this chapter's own; for the Chamberlain model, Wooldridge 2010 §14.6.2, p. 551, notes that MD is asymptotically equivalent to the GMM of §11.1.2 and says that this can be derived from Angrist and Newey 1991)** Suppose the reduced-form parameters are estimated from exactly identifying moment conditions $\mathbb{E}[g(w,\pi_0)] = 0$ ($r$ moments, $r$ parameters), $\hat\pi$ is the corresponding GMM estimate, $D_\pi = \mathbb{E}[\partial g(w,\pi_0)/\partial\pi']$ is invertible, $S_0 = \mathbb{E}[g(w,\pi_0)g(w,\pi_0)']$, and $\hat\Sigma_\pi = \hat D_\pi^{-1}\hat S\hat D_\pi^{-1\prime}$; and suppose MD.2 holds. Then two-step GMM with moment conditions $g(w, h(\theta))$ and optimal MD with input $\hat\pi$ satisfy $\sqrt n(\hat\theta_{GMM} - \hat\theta_{MD}) \xrightarrow{p} 0$, and the difference between the two J statistics goes to zero in probability. If $g$ is linear in $\pi$ and the GMM weight is $\hat S_n(\hat\pi)^{-1}$, the two estimates and the two J statistics are numerically identical.
:::

The proof sketch has three steps. Step one (exact identification): $\bar g_n(\hat\pi) = 0$, and expanding $\bar g_n(h(\theta))$ around $\hat\pi$ gives $\bar g_n(h(\theta)) = -\hat D_\pi(\hat\pi - h(\theta))$ plus higher-order terms; when $g$ is linear in $\pi$, $\hat D_\pi$ is a constant and the equality holds exactly. Step two (MD.1): substituting into the GMM criterion gives $(\hat\pi - h(\theta))'\hat D_\pi'\hat S^{-1}\hat D_\pi(\hat\pi - h(\theta))$; under exact identification $\hat\Sigma_\pi = \hat D_\pi^{-1}\hat S\hat D_\pi^{-1\prime}$ (Theorem 2 with $r = k$), so $\hat D_\pi'\hat S^{-1}\hat D_\pi = \hat\Sigma_\pi^{-1}$, and the GMM criterion is the MD criterion. Step three (Theorems 2 and 4): in a neighborhood of $\theta_0$ of size $n^{-1/2}$ the two criteria differ only by $o_p(n^{-1})$, so the linear representations are the same and the estimators and the J statistics are asymptotically equal; the $\hat S$ of standard two-step GMM is computed at the restricted first step and differs from $\hat S_n(\hat\pi)$ only by $o_p(1)$, so the difference vanishes asymptotically but is not exactly zero.

In words: under the conditions above, MD and GMM estimate the same $\theta_0$ with the same efficiency, and which one to use depends on the data at hand: MD needs only the reduced-form estimates and their variance, which suits combining published experimental results or aggregate data; GMM needs the raw data. When the specification is wrong, two-step GMM computes its weight from the residuals of the restricted first step and so picks up specification error, while MD takes its weight from the unrestricted reduced form; the two pseudo-true values generally differ. Who sets the weights: in both cases, the precision of the reduced form, $\Sigma_\pi^{-1}$, that is, the design. On G1 the difference between two-step GMM and MD is $-2.0 \times 10^{-5}$, and with the weight taken from the unrestricted fit the difference is below $10^{-8}$, which is the linear case in the last sentence of the theorem.

The classic application of MD comes from Chamberlain (1982). This example uses panel data, with variables written in capitals as is customary for panels: units $i = 1, \ldots, N$ and periods $t = 1, \ldots, T$,

$$Y_{it} = X_{it}'\beta + \alpha_i + \varepsilon_{it}, \qquad \mathbb{E}[\varepsilon_{it} \mid X_{i1}, \ldots, X_{iT}, \alpha_i] = 0,$$

where the latter condition is called strict exogeneity: the error is uncorrelated with the regressors of every period. Here $X_{it}$ is a $k \times 1$ vector of regressors, so $k$ is the dimension of $\beta$; the dimension of the structural parameter $\theta$ is counted separately below. The individual effect $\alpha_i$ may be correlated with the regressors. Take its linear projection on the regressors of all periods, $\alpha_i = \xi_0 + \sum_{t=1}^T X_{it}'\xi_t + a_i$; the projection always exists, and $a_i$ is uncorrelated with $X$ in every period. After substitution, in the period-$t$ regression of $Y_{it}$ on $(1, X_{i1}, \ldots, X_{iT})$, only the coefficient on the period-$t$ regressor carries an extra $\beta$. Stacking the periods, the slope part is

$$\Pi = I_T \otimes \beta' + \boldsymbol\iota_T(\xi_1', \ldots, \xi_T'),$$

and the intercept of every period is $\xi_0$. The advantage of the projection is that it needs no additional assumption: $\alpha_i$ may be correlated with the regressors in any way, the correlated part is absorbed by $\xi_t$, and the remainder $a_i$ is uncorrelated with the regressors by definition. The price is a set of restrictions: the true structure has only one $\beta$, yet it must explain $T$ periods and $T \times T$ blocks of coefficients at once. Period-by-period OLS gives the unrestricted $\hat\pi$, made up of the intercepts of each period and $\mathrm{vec}(\Pi')$, with $T(1 + Tk)$ elements in all; the structural parameter $\theta = (\xi_0, \beta', \xi_1', \ldots, \xi_T')'$ has $1 + k + Tk$ elements, and $\pi$ is a linear function of $\theta$ (Wooldridge 2010 equations (14.78) and (14.79)). Theorem 7 delivers all the conclusions directly: the optimal weight is $\hat\Sigma_\pi^{-1}$, where $\hat\Sigma_\pi$ uses the robust variance of period-by-period OLS (equation (14.80)), valid under arbitrary heteroskedasticity and serial correlation; the MD estimator is asymptotically normal; and the overidentification statistic follows a $\chi^2$ with $\dim\pi - \dim\theta = (T - 1)(1 + Tk) - k$ degrees of freedom, and what it tests is strict exogeneity (Wooldridge 2010 p. 551). Section 3.5 of Chapter 10 will use this machinery on panel data.

## 4 Estimation and Practice

### 4.1 Implementing Two-Step GMM

With linear moment conditions every step has a closed-form solution, and the first step usually takes $W^{(1)} = (n^{-1}\sum_i z_iz_i')^{-1}$: in the linear instrumental variables model this step is 2SLS, which first regresses price on the experiments to obtain fitted values and then regresses lessons on the fitted values; Chapter 9 is devoted to it. The identity weight is suitable at most for the first step: it weights the moments by their recorded units, so changing the units changes the answer (Section 3.2); the efficient weight absorbs changes of units into $\hat S$, and CUE is completely insensitive to the units of the moments. A complete procedure with linear moment conditions runs as follows:

1. Write down the moment conditions and explain, one by one, why each holds at the truth, especially the exclusion restriction;
2. Run the first step with the 2SLS weight to obtain $\hat\theta^{(1)}$, and compute $\hat S$ there;
3. Run the second step with $\hat S^{-1}$, and report the estimate, the Windmeijer standard error, and J;
4. Report the strength of the first stage; when it is in doubt, report the S confidence set instead, with CUE or iterated GMM as a check.

There are two ways to write the moment variance. This chapter and the code of Section 6 use the uncentered $\hat S_n(\theta) = n^{-1}\sum_i g(w_i,\theta)g(w_i,\theta)'$; the other way first subtracts the sample mean:

$$\hat S = \frac{1}{n}\sum_{i=1}^n \big(g(w_i,\hat\theta^{(1)}) - \bar g_n\big)\big(g(w_i,\hat\theta^{(1)}) - \bar g_n\big)',$$

where $\bar g_n$ is evaluated at $\hat\theta^{(1)}$. When the moment conditions hold, the two are asymptotically equivalent; Hansen (2014) §13.4 recommends the centered version, on the grounds that when the moment conditions fail, the uncentered version mixes $\bar g_n\bar g_n'$ into the weight and weakens the power of the J test. The CUE minimizer is unaffected by this choice.

Iterated GMM converges in 13 iterations in the G2 sample and in 10.4 on average in the Monte Carlo. The CUE criterion is nonlinear in $\theta$ and requires numerical optimization; under weak moments the criterion may keep falling at large $|\theta|$, and 1,992 of G3's 2,000 replications converge. This is a symptom of weak identification, not a problem with the optimizer.

### 4.2 Which Moments to Use: Where Identification Comes From

The $\Lambda_0$ of Theorem 2 answers a practical question: if some moment is slightly wrong, how far will the estimate move? Andrews, Gentzkow and Shapiro (2017) prove that when the moment conditions are perturbed by an amount of order $n^{-1/2}$, the asymptotic bias of the estimator is $\Lambda_0$ times the perturbation of the moments, and they advocate reporting $\Lambda_0$ together with the estimate. For two-step GMM on G2, multiplying the elasticity row by the standard deviation of each moment shows how far the estimate moves when that moment is perturbed by one standard deviation. The standard deviation here is that of the moment function for a single observation; if the perturbation is instead one standard error of the sample mean, the movement must also be divided by $\sqrt n$. The largest movements come from the premium-tutor upcharge experiment (3.36), the loyalty-rebate experiment ($-2.82$), and the booking-fee experiment (2.81); the deep-discount experiment gives $-2.53$, and the fee-waiver threshold only 0.27. So the elasticity is identified mainly by the top few experiments: if you worry that the premium-tutor upcharge also changed parents' perception of quality, you should know that it is the moment on which the estimate depends most. In MD this is even more transparent: when $h(\theta) = \theta\,\boldsymbol\iota_r$, $\Lambda_0$ is the weight vector itself, and in G1 the weights of Seattle and Boston are 0.1735 and 0.1534 (sample estimates 0.171 and 0.152). $\Lambda_0$ is local: it describes small perturbations of the moment conditions, and large misspecification calls for a separate analysis. Choosing which moments to use determines which part of the variation in the data identifies the parameter.

### 4.3 The Estimators at a Glance

| Estimator | Weight | Standard error | When it is credible | G2 sample | G2 coverage ($n = 200$) |
|---|---|---|---|---|---|
| OLS | None | Robust | Exogenous price | $-0.318$ (0.290) | |
| 2SLS | $(n^{-1}\sum_i z_iz_i')^{-1}$ | Sandwich | GMM.1 through GMM.3 | $-1.843$ (0.414) | 0.9305 |
| Identity weight | $I_r$ | Sandwich | As above; inefficient and unit-dependent | $-1.307$ (0.692) | 0.9395 |
| Two-step GMM | $\hat S_n(\hat\theta^{(1)})^{-1}$ | Conventional | As above; large samples | $-2.012$ (0.358) | 0.9115 |
| Two-step GMM | As above | Windmeijer | As above; linear moment conditions | $-2.012$ (0.418) | 0.935 |
| Iterated GMM | $\hat S_n^{-1}$ at convergence | Conventional | As for two-step | $-2.043$ (0.368) | 0.906 |
| CUE | $\hat S_n(\theta)^{-1}$, moving with $\theta$ | Conventional | Small bias; standard errors too small | $-2.304$ (0.382) | 0.9165 |
| S confidence set | $\hat S_n(\theta)^{-1}$ | Test inversion | Needs only SW.1; holds under weak moments | $[-4.09, -0.995]$ | |

The identity weight's coverage is not bad; where it loses is interval width, since its Monte Carlo standard deviation is 2.02 times that of two-step GMM.

### 4.4 Software and Simulated Moments

In R, the `gmm` package's `gmm()` gives the three estimators with `type = "twoStep"`, `"iterative"`, or `"cue"`; `vcov = "MDS"` corresponds to independent observations, and `centeredVcov` chooses whether to demean; the `momentfit` package's `gmmFit()` is another interface. The two-step variance that `gmm()` reports is computed at $\hat\theta^{(2)}$ with a recomputed $\hat S$, whereas this chapter follows Windmeijer's convention and uses $\hat S_n(\hat\theta^{(1)})$; the two are asymptotically the same. Neither package offers the Windmeijer correction, though `plm::pgmm` for linear panel GMM does; there are no ready-made general-purpose functions for MD or the S confidence set, and Section 6 writes them by hand. Stata's `gmm` command has two-step GMM and iterated GMM, and the community-contributed command `ivreg2` has CUE.

When the expectations in the moment conditions have no closed form but can be simulated from the model, replacing the model predictions in $\bar g_n(\theta)$ with simulated averages gives the simulated method of moments (SMM; McFadden 1989; Pakes and Pollard 1989); simulation error makes the variance larger but does not change the framework of this chapter. Indirect inference (Gourieroux, Monfort and Renault 1993) goes one step further: first estimate an auxiliary model on the data to obtain $\hat\pi$, then find the $\theta$ for which the auxiliary model estimated on simulated data comes closest to it; this is an MD whose $h(\theta)$ is computed by simulation. The demand estimation of Chapter 28 computes market shares by simulated integration, which is exactly this idea; indirect inference is also common in the literature on dynamic structural models.

## 5 Anchoring Papers

Each of the three anchoring papers corresponds to one thread of the chapter: the large-sample theory of GMM, how MD pulls the reduced form of a panel back to structural parameters, and a platform application that learns a demand curve from randomized price experiments. In the IS and marketing literatures we have not yet found applications that use GMM or MD to combine several experiments into a structural parameter; the third paper is the closest to Tamarind in design, but its estimation method is not GMM. Each paper is laid out in five parts: paper, method, data, results, and limitations.

### 5.1 Hansen (1982)

::: {.case}
Paper: "Large Sample Properties of Generalized Method of Moments Estimators," Econometrica 50(4):1029-1054. It laid down the large-sample theory of GMM and is the original source of Theorems 1 to 4 of this chapter; this chapter states these results using the numbering of Wooldridge (2010) and Hansen (2014).

Method: estimation is built on a set of population moment conditions given by the economic model, without requiring a complete distributional assumption, and the data may be a stationary ergodic time series. The paper proves that the estimator that minimizes the weighted quadratic form of the moment conditions is consistent and asymptotically normal, gives the efficient weight and the minimum variance it attains, and proposes the overidentification test constructed from the objective function value, later called the Hansen J test, with degrees of freedom equal to the degree of overidentification.

Data: this is a theoretical paper with no empirical data; its classic application is the Euler equation for consumption and asset returns under rational expectations, models that naturally give moment conditions rather than a complete likelihood.

Results: GMM nests OLS, instrumental variables, and 2SLS as special cases, so researchers who are willing to maintain only a few moment conditions and dare not write down a complete distribution can still estimate efficiently and test the specification. It went on to become the workhorse tool of structural estimation, especially in industrial organization and dynamic models.

Limitations: the properties are all asymptotic. The weight matrix of the two-step estimate carries noise in small samples, and the standard errors are too small (Section 3.5); with many moment conditions the two-step estimate is biased; with weak moment conditions the normal approximation fails (Section 3.6). Nor does a rejection by J say which moment is wrong.
:::

Hansen's contribution pushes Chapter 5's unified framework as far as it goes: estimation need not start from a distribution, and you can get to work as long as you can write down credible moment conditions. Its line of defense is also clear: it concentrates credibility on a few economically grounded moment conditions rather than on a complete likelihood that may be misspecified, which is of a piece with the robust spirit of the QMLE of Section 3.6 of Chapter 6.

### 5.2 Chamberlain (1982)

::: {.case}
Paper: "Multivariate Regression Models for Panel Data," Journal of Econometrics 18(1):5-46. The question is what panel data can tell us when the individual effect is correlated with the regressors. The description of the method below follows the NBER technical paper of the same title, No. 8 (December 1980), and the account in Wooldridge (2010).

Method: instead of eliminating the individual effect, the paper takes its linear projection on the regressors of all periods, so that the multivariate regression of each period's outcome on the regressors of all periods, the $\pi$ matrix, must satisfy a set of restrictions determined by the structural parameters (Section 3.8). The true conditional expectation is generally not linear, so inference is carried out for the linear projection, based on the asymptotic distribution of sample moment functions; the restrictions are imposed by minimum distance, which the author points out is generally more efficient than conventional estimators. The overidentification test is a test of strict exogeneity.

Data: the method is aimed at panels with many units and few periods. Wooldridge (2010) p. 551 attributes the empirical examples that test "the individual effect is unrelated to the regressors" in the MD framework to the same author's handbook chapter, Chamberlain (1984).

Results: the paper provides a framework for panel estimation and testing that does not rely on homoskedastic or independent errors, which Wooldridge (2010) §11.1.2 and §14.6.2 present as standard textbook material; in the linear model with an additive individual effect, its overidentification test is a test of strict exogeneity.

Limitations: the dimension of the reduced form grows with the square of the number of periods, so when $T$ is even moderately large $\hat\Sigma_\pi$ is poorly estimated, and the optimal weight is no longer optimal in small samples, a problem with the same root as that of two-step GMM in Section 3.5; the approach handles strict exogeneity; when the regressors respond to past shocks, other methods are needed, and Chapter 11 develops them.
:::

This paper made MD a standard tool of panel econometrics: first estimate an unrestricted reduced form that is robust to heteroskedasticity and serial correlation, then pull it back onto the structural parameters, and the distance covered in pulling it back is itself a specification test.

### 5.3 Dubé and Misra (2023)

::: {.case}
Paper: "Personalized Pricing and Consumer Welfare," Journal of Political Economy 131(1):131-189. The question is how profits and consumer welfare change when a digital platform uses machine learning for personalized pricing. The numbers below follow the published journal version.

Method: randomized price experiments on the recruiting platform ZipRecruiter. The experimental data are used to estimate a logit demand model whose price coefficient varies with customer characteristics; the coefficients are estimated by LASSO (regression with an L1 penalty) under a weighted likelihood bootstrap, and the optimal uniform price and personalized prices are then computed from the estimated model. The authors point out that because prices are randomly assigned, there is none of the price endogeneity that is common when demand is estimated from observational data; they also ran a second experiment to test these pricing recommendations out of sample.

Data: in the first experiment, 7,867 prospective customers were randomly assigned to ten price groups, with the control group at the then-current monthly fee of USD 99; the second experiment randomly assigned newly arriving prospective customers to three groups: USD 99 (25%), USD 249 (25%), and personalized pricing (50%).

Results: doubling the price from USD 99 to USD 199 lowers conversion by only 25%. At USD 99, the posterior mean of the price elasticity of demand is only $-0.33$ (95% credible interval $-0.41$ to $-0.26$); at USD 249 it is $-0.82$, and at USD 399 it is $-1.15$; the optimal uniform price is USD 327, with expected profit 55% higher than at USD 99. Personalized pricing raises profit by a further 19% over the optimal uniform price, and by 86% over USD 99; consumer surplus falls by about 25% relative to optimal uniform pricing, but more than 60% of consumers get a lower price. In the second experiment, the uniform price of USD 249 raises profit per customer by 67.74% relative to the control group, and personalized pricing by 84.4%; the test that personalized pricing beats the USD 249 uniform price has $p = 0.069$.

Limitations: the estimation methods are the Bayesian bootstrap and LASSO, not GMM or MD; the experiments were run on one platform over one period of time, competitors' responses are not in the model, and the authors also state that the analysis does not cover personalized pricing in competitive markets.
:::

The paper makes the same point as Tamarind: random price variation is the cleanest source from which to learn a demand elasticity, and a platform often has a great deal of such variation. Dubé and Misra put the ten price groups into a single likelihood; when the experiments are scattered across markets and leave behind only their own slopes, or when prices also reflect managers' promotions, the MD and GMM of Section 3 are the tools for combining them into one structural parameter, and whose weights you combine them with determines whose elasticity you end up reporting.

## 6 A Full Walkthrough on Tamarind Data

The code below runs under R 4.6.1; several statements are merged onto one line, and statements used only for checks are omitted. Every number in the text comes from the output of an actual run, and the result lines in the code that start with `#` copy what the run printed. The code follows the notation of the text: `V_pi` is the variance of $\hat\pi$ itself, that is, $\hat\Sigma_\pi/n$; `W1` is $W^{(1)}$; and `D` is $\hat D$.

### 6.1 The DGP

G1 and G1h share one test design, but each has its own seed; the city elasticities of G1h vary linearly with log income.

```r
P1 <- list(markets = c("Atlanta", "Austin", "Boston", "Charlotte", "Chicago", "Columbus", "Dallas", "Denver", "Houston",
                       "Miami", "Minneapolis", "Nashville", "Philadelphia", "Phoenix", "Portland", "Seattle"),
           n = c(4000, 2400, 6000, 2000, 8000, 2000, 5200, 3600, 6400, 3000, 2800, 2200, 5600, 4400, 2600, 4800),
           delta = c(0.08, 0.05, 0.10, 0.12, 0.06, 0.03, 0.10, 0.08, 0.05, 0.12, 0.10, 0.06, 0.08, 0.04, 0.10, 0.12),
           p0 = c(9.99, 10.49, 11.49, 9.49, 10.99, 8.99, 9.99, 10.49, 9.49, 10.99, 10.49, 9.99, 11.49, 9.49, 10.99, 11.49),
           q0 = c(3.6, 3.8, 3.2, 3.9, 3.4, 4.1, 3.7, 3.5, 3.9, 3.3, 3.6, 3.8, 3.1, 4.0, 3.4, 3.2),
           sigma = c(0.62, 0.70, 0.58, 0.75, 0.60, 0.80, 0.66, 0.64, 0.72, 0.68, 0.62, 0.74, 0.60, 0.78, 0.64, 0.58),
           het = 1.1, theta = -1.6, R_mc = 2000,
           income = c(72, 86, 94, 68, 78, 64, 70, 88, 62, 66, 84, 74, 76, 66, 82, 104),             # USD thousands
           households = c(310, 180, 290, 150, 520, 120, 400, 210, 450, 260, 170, 140, 380, 230, 160, 270),
           kappa = 3)
gen_G1 <- function(par, theta_m = rep(par$theta, length(par$markets))) {
  d <- do.call(rbind, lapply(seq_along(par$markets), function(m) {
    nm <- par$n[m]; B <- sample(rep(0:1, length.out = nm))           # arm B pays p0 (1 + delta)
    dlp <- B * log(1 + par$delta[m]); s <- ifelse(B == 1, par$het, 1)
    data.frame(market = m, B = B, log_price = log(par$p0[m]) + dlp,
               y = log(par$q0[m]) + theta_m[m] * dlp + par$sigma[m] * s * rnorm(nm))   # log weekly lessons booked
  }))
  d$market <- factor(par$markets[d$market], levels = par$markets); d
}
set.seed(7101); G1 <- gen_G1(P1)
theta_m <- P1$theta + P1$kappa * (log(P1$income) - mean(log(P1$income)))   # G1h: less elastic where income is higher
set.seed(7102); G1h <- gen_G1(P1, theta_m)
avar <- with(P1, sigma^2 * (1 / (n / 2) + het^2 / (n / 2)) / log(1 + delta)^2)   # exact variance of each market slope
w_md <- (1 / avar) / sum(1 / avar)                                                # optimal MD weights
cat(sprintf("weights %.4f (%s) to %.4f (%s); theta_star %.4f, mean %.4f, household-weighted %.4f\n",
            max(w_md), P1$markets[which.max(w_md)], min(w_md), P1$markets[which.min(w_md)], sum(w_md * theta_m),
            mean(theta_m), sum(P1$households * theta_m) / sum(P1$households)))
# weights 0.1735 (Seattle) to 0.0026 (Columbus); theta_star -1.3953, mean -1.6000, household-weighted -1.6192
```

G2 and G3 share one generating function: each experiment puts regions in its experimental arm according to its assignment share, and the instruments are the centered assignments; G3 merely replaces the experiments with four perturbations of a few cents.

```r
P2 <- list(n = 200, step = c(3, 10, 0.8, 4, 1.5, 6, 20, 8, 2, 5), share = c(0.5, 0.2, rep(0.5, 8)),
           gamma = c(1 / 83.92, -0.006, 1 / 83.92, 0.4 / 83.92, 0.3 / 10.49, -0.5 / 83.92, 0.001, -0.004, 0.6 / 83.92, -0.003),
           sigma_v = 0.03, rho = 0.8, sigma_eps = 0.10, het = 3, het_on = 2,
           p_base = 10.49, q_base = 60, theta = -1.6, R_mc = 2000, grid_n = c(100, 200, 400, 800), R_grid = 1000)
P3 <- modifyList(P2, list(n = 1000, step = c(0.016, 0.024, 0.032, 0.040), share = rep(0.5, 4), gamma = rep(1 / 10.49, 4),
                          het = 0, het_on = 1, R_mc = 2000))                 # G3: four penny-sized list-price perturbations
gen_G2 <- function(par, n = par$n) {
  K <- length(par$step); n_up <- round(n * par$share)
  Tm <- vapply(seq_len(K), function(j) sample(rep(c(1, 0), c(n_up[j], n - n_up[j]))), numeric(n))   # upper arm of experiment j
  Zs <- sweep(sweep(Tm, 2L, par$share), 2L, par$step, "*")              # z_j = step_j (T_j - share_j): sample mean exactly 0
  e <- rnorm(n); xi <- rnorm(n)
  v <- par$sigma_v * (par$rho * e + sqrt(1 - par$rho^2) * xi)          # managers' promotions, correlated with the demand shock
  dlp <- drop(Zs %*% par$gamma) + v
  sig <- par$sigma_eps * (1 + par$het * Tm[, par$het_on])             # demand noise four times larger in deep-discount regions
  d <- data.frame(y = log(par$q_base) + par$theta * dlp + sig * e, log_price = log(par$p_base) + dlp)
  for (j in seq_len(K)) d[[paste0("z", j)]] <- Zs[, j]
  d
}
set.seed(7103); G2 <- gen_G2(P2)
set.seed(7104); G3 <- gen_G2(P3)
```

### 6.2 Blocks G1 and G1h: Minimum Distance

`market_slopes` runs OLS within each city, with HC0 standard errors; `md_fit` minimizes the MD criterion by Gauss-Newton, and when `V_pi` is diagonal it reduces to the inverse-variance weighted average.

```r
market_slopes <- function(d) do.call(rbind, lapply(split(d, d$market), function(s) {
  x <- s$log_price - mean(s$log_price); b <- sum(x * s$y) / sum(x^2); e <- (s$y - mean(s$y)) - b * x
  data.frame(pi_hat = b, se = sqrt(sum(x^2 * e^2)) / sum(x^2))          # OLS slope within the test, HC0 SE
}))
md_fit <- function(pi_hat, V_pi, h, theta0, W = solve(V_pi), tol = 1e-12) {   # V_pi: variance of pi_hat
  theta <- theta0
  repeat {                                                                    # Gauss-Newton on (pi_hat - h)' W (pi_hat - h)
    H <- numDeriv::jacobian(h, theta); e <- pi_hat - h(theta)
    step <- drop(solve(crossprod(H, W %*% H), crossprod(H, W %*% e))); theta <- theta + step
    if (max(abs(step)) < tol * (1 + max(abs(theta)))) break
  }
  H <- numDeriv::jacobian(h, theta); B <- solve(crossprod(H, W %*% H)); e <- pi_hat - h(theta)
  J <- drop(crossprod(e, W %*% e)); df <- length(pi_hat) - length(theta)
  list(theta = theta, vcov = B %*% crossprod(H, W %*% V_pi %*% W %*% H) %*% B, J = J, df = df,
       p = pchisq(J, df, lower.tail = FALSE))
}
h_hom <- function(th) rep(th, 16)                                            # h(theta) = theta iota_16
md_block <- function(d) {
  ms <- market_slopes(d); V <- diag(ms$se^2)
  f <- md_fit(ms$pi_hat, V, h_hom, 0); eq <- md_fit(ms$pi_hat, V, h_hom, 0, W = diag(16))   # optimal and equal weights
  list(ms = ms, theta = f$theta, se = sqrt(f$vcov[1, 1]), J = f$J, p = f$p,
       ivw = sum(ms$pi_hat / ms$se^2) / sum(1 / ms$se^2), eq = eq$theta, se_eq = sqrt(eq$vcov[1, 1]))
}
for (b in list(md_block(G1), md_block(G1h)))
  cat(sprintf("MD %.4f (SE %.4f) [%.4f, %.4f]; IVW mean %.4f; J %.2f (p %.3g); equal weights %.4f (SE %.4f)\n",
              b$theta, b$se, b$theta - qnorm(0.975) * b$se, b$theta + qnorm(0.975) * b$se, b$ivw, b$J, b$p, b$eq, b$se_eq))
# MD -1.5352 (SE 0.0645) [-1.6617, -1.4087]; IVW mean -1.5352; J 14.64 (p 0.478); equal weights -1.5282 (SE 0.1169)
# MD -1.4059 (SE 0.0642) [-1.5317, -1.2802]; IVW mean -1.4059; J 74.85 (p 6.03e-10); equal weights -1.5376 (SE 0.1173)
set.seed(7111)
sim <- t(vapply(1:2000, function(i) {
  ms <- market_slopes(gen_G1(P1)); f <- md_fit(ms$pi_hat, diag(ms$se^2), h_hom, 0)
  c(theta = f$theta, se = sqrt(f$vcov[1, 1]), p = f$p)
}, numeric(3)))
cat(sprintf("Monte Carlo: J rejects at 5%% in %.4f; MD interval covers in %.4f; sd %.4f, mean SE %.4f\n",
            mean(sim[, "p"] < 0.05), mean(abs(sim[, "theta"] + 1.6) <= qnorm(0.975) * sim[, "se"]),
            sd(sim[, "theta"]), mean(sim[, "se"])))
# Monte Carlo: J rejects at 5% in 0.0460; MD interval covers in 0.9435; sd 0.0640, mean SE 0.0647
```

For G1, the MD estimate and the inverse-variance weighted average agree exactly at printed precision; for G1h, J rejects homogeneity, and the MD estimate aims at $\theta_* = -1.395$. The equal-weight estimate is less precise in both blocks: its standard errors are 1.81 and 1.83 times those of optimal MD.

### 6.3 GMM Tools, and MD Against GMM

`gmm_min` minimizes the GMM criterion for a fixed weight; `gmm_fit` runs the first step, then the second, and then iterated GMM or CUE if requested. Below we use it for two-step GMM on G1's household data and compare the result with MD.

```r
gbar_jac <- function(gfun, theta, d) numDeriv::jacobian(function(t) colMeans(gfun(t, d)), theta)   # D-hat
moment_var <- function(G) crossprod(G) / nrow(G)                             # S_n(theta) = n^-1 sum g g' (uncentred)
gmm_min <- function(gfun, theta, d, W, tol = 1e-10) {                         # argmin gbar' W gbar by Gauss-Newton
  Q <- function(th) { gb <- colMeans(gfun(th, d)); drop(crossprod(gb, W %*% gb)) }
  q_cur <- Q(theta)
  for (iter in 1:100) {
    D <- gbar_jac(gfun, theta, d); gb <- colMeans(gfun(theta, d))
    step <- -drop(solve(crossprod(D, W %*% D), crossprod(D, W %*% gb))); lambda <- 1
    repeat {                                                                  # halve the step until the criterion falls
      cand <- theta + lambda * step; q_new <- Q(cand)
      if (q_new <= q_cur * (1 + 1e-12) || max(abs(step)) < sqrt(.Machine$double.eps) * (1 + max(abs(theta)))) break
      lambda <- lambda / 2
    }
    if (iter > 1 && q_cur - q_new <= 1e-12 * q_cur) break                     # no further fall beyond rounding
    theta <- cand; q_cur <- q_new
    if (max(abs(step)) < tol * (1 + max(abs(theta)))) break
  }
  list(theta = theta, obj = q_cur, D = D)
}
gmm_fit <- function(gfun, theta0, d, type = c("two_step", "iterated", "cue"), W1) {
  type <- match.arg(type); n <- nrow(gfun(theta0, d))
  s1 <- gmm_min(gfun, theta0, d, W1); theta1 <- s1$theta                     # first step, weight W^(1)
  W <- solve(moment_var(gfun(theta1, d))); s <- gmm_min(gfun, theta1, d, W)  # second step, weight S_n(theta^(1))^-1
  theta <- s$theta; D <- s$D; iter <- 2
  if (type == "iterated") repeat {                                           # re-estimate the weight until theta settles
    W <- solve(moment_var(gfun(theta, d))); s <- gmm_min(gfun, theta, d, W); iter <- iter + 1
    change <- max(abs(s$theta - theta)); theta <- s$theta; D <- s$D
    if (change < 1e-10 * (1 + max(abs(theta)))) break
  }
  if (type == "cue") {                                                       # weight moves with theta
    crit <- function(th) { G <- gfun(th, d); gb <- colMeans(G); drop(crossprod(gb, solve(moment_var(G), gb))) }
    op <- optim(theta, crit, gr = function(th) numDeriv::grad(crit, th), method = "BFGS",
                control = list(reltol = 1e-15, maxit = 5000))
    theta <- op$par; W <- solve(moment_var(gfun(theta, d))); D <- gbar_jac(gfun, theta, d); ok <- op$convergence == 0
  } else ok <- TRUE
  gb <- colMeans(gfun(theta, d))
  list(theta = theta, theta1 = theta1, vcov = solve(crossprod(D, W %*% D)) / n, J = n * drop(crossprod(gb, W %*% gb)),
       df = ncol(W) - length(theta), D = D, W = W, iter = iter, converged = ok)
}
```

```r
mk <- as.integer(G1$market); Dm <- outer(mk, 1:16, "==") * 1
pdat <- list(y = G1$y, lp = G1$log_price, mk = mk, Z = cbind(Dm, Dm * G1$log_price))   # 32 instruments
g_pool <- function(b, d) d$Z * (d$y - b[d$mk] - b[17] * d$lp)            # b = (a_1, ..., a_16, theta)
n1 <- length(pdat$y)
g2 <- gmm_fit(g_pool, rep(0, 17), pdat, "two_step", W1 = solve(crossprod(pdat$Z) / n1))   # first step: pooled OLS
md1 <- md_block(G1); ms <- md1$ms
a_hat <- as.numeric(tapply(G1$y, G1$market, mean) - ms$pi_hat * tapply(G1$log_price, G1$market, mean))
e_u <- pdat$y - a_hat[mk] - ms$pi_hat[mk] * pdat$lp                         # residuals of the market-by-market fits
S_u <- moment_var(pdat$Z * e_u)
gu <- gmm_min(g_pool, g2$theta, pdat, solve(S_u))
cat(sprintf("two-step GMM %.6f (SE %.4f, J %.3f, df %d); MD %.6f (J %.3f); MD - GMM %.1e\n",
            g2$theta[17], sqrt(g2$vcov[17, 17]), g2$J, g2$df, md1$theta, md1$J, md1$theta - g2$theta[17]))
cat(sprintf("unrestricted weight: estimate and J equal to MD within 1e-8: %s\n",
            max(abs(c(md1$theta - gu$theta[17], md1$J - n1 * gu$obj))) < 1e-8))
# two-step GMM -1.535155 (SE 0.0645, J 14.627, df 15); MD -1.535176 (J 14.636); MD - GMM -2.0e-05
# unrestricted weight: estimate and J equal to MD within 1e-8: TRUE
```

Household-level two-step GMM and MD differ by $-2.0 \times 10^{-5}$; the first-step weight comes from the restricted pooled regression, so the two are only asymptotically equal. Once the weight uses the residuals of the city-by-city regressions instead, the estimate and J agree to within $10^{-8}$, which is the linear case of Theorem 8.

### 6.4 Block G2: Six Estimators on One Sample

`windmeijer_vcov` applies only to linear moment conditions; `s_stat` is the S statistic, and the S confidence set is found by inversion on a grid.

```r
iv_data <- function(df) {
  Zs <- unname(as.matrix(df[, grep("^z[0-9]+$", names(df))]))
  list(y = df$y, X = cbind(1, df$log_price), Z = cbind(1, Zs), Zs = Zs, yc = df$y - mean(df$y), pc = df$log_price - mean(df$log_price))
}
g_iv <- function(theta, d) d$Z * drop(d$y - d$X %*% theta)                  # z_i (y_i - beta0 - theta log p_i)
g_conc <- function(theta, d) d$Zs * (d$yc - theta * d$pc)                    # intercept partialled out
sandwich_gmm <- function(D, W, S, n) { B <- solve(crossprod(D, W %*% D)); B %*% crossprod(D, W %*% S %*% W %*% D) %*% B / n }
windmeijer_vcov <- function(gfun, theta1, theta2, d, W1) {                   # linear moments only
  G1m <- gfun(theta1, d); n <- nrow(G1m); r <- ncol(G1m); k <- length(theta1)
  S1i <- solve(moment_var(G1m)); D <- gbar_jac(gfun, theta2, d)
  V2 <- solve(crossprod(D, S1i %*% D)) / n                                   # conventional two-step variance
  V1 <- sandwich_gmm(gbar_jac(gfun, theta1, d), W1, moment_var(G1m), n)     # robust variance of the first step
  dG <- numDeriv::jacobian(function(t) as.vector(gfun(t, d)), theta1)
  A <- solve(crossprod(D, S1i %*% D), crossprod(D, S1i)); gb2 <- colMeans(gfun(theta2, d))
  Dw <- sapply(1:k, function(j) { Gj <- matrix(dG[, j], n, r)                # d theta^(2) / d theta^(1)
    A %*% ((crossprod(Gj, G1m) + crossprod(G1m, Gj)) / n) %*% S1i %*% gb2 })
  Dw <- matrix(Dw, k, k)
  V2 + Dw %*% V2 + V2 %*% t(Dw) + Dw %*% V1 %*% t(Dw)
}
s_stat <- function(gfun, theta, d) { G <- gfun(theta, d); gb <- colMeans(G); nrow(G) * drop(crossprod(gb, solve(moment_var(G), gb))) }
fit_all <- function(dl) {
  n <- length(dl$y); r <- ncol(dl$Z); W1 <- solve(crossprod(dl$Z) / n)    # W^(1) = (Z'Z/n)^-1: the first step is 2SLS
  f2 <- gmm_fit(g_iv, c(0, 0), dl, "two_step", W1 = W1)
  V1 <- sandwich_gmm(gbar_jac(g_iv, f2$theta1, dl), W1, moment_var(g_iv(f2$theta1, dl)), n)
  Vw <- windmeijer_vcov(g_iv, f2$theta1, f2$theta, dl, W1)
  fi <- gmm_fit(g_iv, f2$theta1, dl, "iterated", W1 = W1); fc <- gmm_fit(g_iv, f2$theta1, dl, "cue", W1 = W1)
  fid <- gmm_min(g_iv, f2$theta1, dl, diag(r))                               # identity weight, recorded units
  V_id <- sandwich_gmm(fid$D, diag(r), moment_var(g_iv(fid$theta, dl)), n)
  rbind(tsls = c(f2$theta1[2], sqrt(V1[2, 2])), identity = c(fid$theta[2], sqrt(V_id[2, 2])),
        two_step = c(f2$theta[2], sqrt(f2$vcov[2, 2])), windmeijer = c(f2$theta[2], sqrt(Vw[2, 2])),
        iterated = c(fi$theta[2], sqrt(fi$vcov[2, 2])), cue = c(fc$theta[2], sqrt(fc$vcov[2, 2])))
}
dl <- iv_data(G2); est <- fit_all(dl)
ci <- cbind(est, est[, 1] - qnorm(0.975) * est[, 2], est[, 1] + qnorm(0.975) * est[, 2])
for (k in rownames(ci)) cat(sprintf("%-10s %.4f  SE %.4f  [%.4f, %.4f]\n", k, ci[k, 1], ci[k, 2], ci[k, 3], ci[k, 4]))
# tsls       -1.8431  SE 0.4139  [-2.6543, -1.0319]
# identity   -1.3069  SE 0.6921  [-2.6634, 0.0495]
# two_step   -2.0121  SE 0.3579  [-2.7135, -1.3107]
# windmeijer -2.0121  SE 0.4176  [-2.8306, -1.1936]
# iterated   -2.0429  SE 0.3675  [-2.7633, -1.3226]
# cue        -2.3044  SE 0.3820  [-3.0532, -1.5557]
ols <- lm(y ~ log_price, data = G2)
f2 <- gmm_fit(g_iv, c(0, 0), dl, "two_step", W1 = solve(crossprod(dl$Z) / 200))
grid <- seq(-8, 2, by = 0.005); S <- sapply(grid, function(th) s_stat(g_conc, th, dl))
in_set <- grid[S <= qchisq(0.95, 10)]                                         # 10 shock moments, intercept partialled out
cat(sprintf("OLS %.4f (SE %.4f); J %.3f (df %d, p %.3f); S set [%.3f, %.3f]\n", coef(ols)[2],
            sqrt(sandwich::vcovHC(ols, type = "HC0")[2, 2]), f2$J, f2$df, pchisq(f2$J, f2$df, lower.tail = FALSE),
            min(in_set), max(in_set)))
# OLS -0.3183 (SE 0.2903); J 9.813 (df 9, p 0.366); S set [-4.090, -0.995]
Lambda <- -solve(crossprod(f2$D, f2$W %*% f2$D), crossprod(f2$D, f2$W))   # sensitivity of (beta0, theta) to each moment
round(Lambda[2, ] * sqrt(diag(moment_var(g_iv(f2$theta, dl)))), 2)            # theta row, per standard deviation of each moment
#  [1]  0.00  2.81 -2.53  0.39  1.77  3.36 -2.82  0.27 -0.38 -0.38 -1.27
```

The last line is the sensitivity of Section 4.2, listing the intercept first and then the ten experiments; the premium-tutor upcharge experiment (in position 6) has the largest value, 3.36. The intercept entry is exactly zero in theory, and the printed 0.00 is only rounding error: the experiment assignments have a sample mean of exactly zero, $\beta_0$ appears only in the intercept moment, and any perturbation of that moment is absorbed by $\beta_0$, leaving the elasticity unaffected.

### 6.5 Block G2: Monte Carlo

```r
run_mc <- function(nn, R, seed) {
  set.seed(seed)
  sims <- lapply(1:R, function(i) fit_all(iv_data(gen_G2(P2, n = nn))))
  sapply(rownames(sims[[1]]), function(k) {
    est <- sapply(sims, function(s) s[k, 1]); se <- sapply(sims, function(s) s[k, 2])
    c(mean = mean(est), sd = sd(est), mean_se = mean(se), cover = mean(abs(est + 1.6) <= qnorm(0.975) * se))
  })
}
mc200 <- run_mc(200, 2000, 7113)
round(mc200, 4)
#            tsls identity two_step windmeijer iterated     cue
# mean    -1.5233  -1.5407  -1.5167    -1.5167  -1.5140 -1.5976
# sd       0.3676   0.6097   0.3019     0.3019   0.3005  0.3208
# mean_se  0.3612   0.6065   0.2780     0.3030   0.2769  0.2799
# cover    0.9305   0.9395   0.9115     0.9350   0.9060  0.9165
cover_n <- sapply(c(100, 400, 800), function(nn) run_mc(nn, 1000, 73100 + match(nn, P2$grid_n))["cover", ])
colnames(cover_n) <- c("n100", "n400", "n800"); round(t(cover_n), 4)
#       tsls identity two_step windmeijer iterated   cue
# n100 0.920    0.940    0.876      0.928    0.876 0.870
# n400 0.939    0.931    0.922      0.932    0.918 0.922
# n800 0.947    0.943    0.941      0.945    0.942 0.939
```

The `two_step` and `windmeijer` columns of the table are the numbers cited in Section 3.5. The Monte Carlo standard error of a coverage rate is about 0.0049 with 2,000 replications and 0.0069 with 1,000, and the gap between the conventional interval's 0.9115 and the nominal level far exceeds this noise.

### 6.6 Block G3: Weak Moments

```r
dl3 <- iv_data(G3); n3 <- 1000; W1 <- solve(crossprod(dl3$Zs) / n3)
f3 <- gmm_fit(g_conc, 0, dl3, "two_step", W1 = W1)
se3 <- sqrt(f3$vcov[1, 1])
first_stage_F <- function(d) {                                                 # robust Wald of the four experiments on price / 4
  ZZi <- solve(crossprod(d$Zs)); b <- ZZi %*% crossprod(d$Zs, d$pc); e <- drop(d$pc - d$Zs %*% b)
  drop(t(b) %*% solve(ZZi %*% crossprod(d$Zs * e) %*% ZZi, b)) / ncol(d$Zs)
}
s_inf <- function(d) { G <- d$Zs * d$pc; nrow(G) * drop(crossprod(colMeans(G), solve(moment_var(G), colMeans(G)))) }
grid3 <- seq(-12, 8, by = 0.01); S3 <- sapply(grid3, function(th) s_stat(g_conc, th, dl3))
cat(sprintf("F %.2f; two-step %.4f (SE %.4f) [%.3f, %.3f]; S(-1.6) %.3f; S set [%.2f, %.2f]; S at infinity %.3f (critical %.3f)\n",
            first_stage_F(dl3), f3$theta, se3, f3$theta - qnorm(0.975) * se3, f3$theta + qnorm(0.975) * se3,
            s_stat(g_conc, -1.6, dl3), min(grid3[S3 <= qchisq(0.95, 4)]), max(grid3[S3 <= qchisq(0.95, 4)]), s_inf(dl3), qchisq(0.95, 4)))
# F 3.61; two-step -0.2289 (SE 0.6396) [-1.482, 1.025]; S(-1.6) 3.064; S set [-7.04, 1.29]; S at infinity 14.014 (critical 9.488)
set.seed(7114)
sim3 <- t(vapply(1:2000, function(i) {
  d <- iv_data(gen_G2(P3)); W1 <- solve(crossprod(d$Zs) / n3); f <- gmm_fit(g_conc, 0, d, "two_step", W1 = W1)
  fc <- tryCatch(gmm_fit(g_conc, f$theta1, d, "cue", W1 = W1), error = function(e) NULL)
  ok <- !is.null(fc) && fc$converged                                          # CUE can drift off under weak moments
  c(est = f$theta, se = sqrt(f$vcov[1, 1]), cue = if (ok) fc$theta else NA, cue_se = if (ok) sqrt(fc$vcov[1, 1]) else NA,
    S0 = s_stat(g_conc, -1.6, d), S_inf = s_inf(d), J = f$J)
}, numeric(7)))
cat(sprintf("coverage: Wald %.4f, CUE %.4f (%d converged), S set %.4f; unbounded S sets %.4f; median two-step %.3f; J rejects %.3f\n",
            mean(abs(sim3[, "est"] + 1.6) <= qnorm(0.975) * sim3[, "se"]),
            mean(abs(sim3[, "cue"] + 1.6) <= qnorm(0.975) * sim3[, "cue_se"], na.rm = TRUE), sum(!is.na(sim3[, "cue"])),
            mean(sim3[, "S0"] <= qchisq(0.95, 4)),
            mean(sim3[, "S_inf"] <= qchisq(0.95, 4)), median(sim3[, "est"]), mean(sim3[, "J"] > qchisq(0.95, 3))))
# coverage: Wald 0.7740, CUE 0.8901 (1992 converged), S set 0.9480; unbounded S sets 0.3685; median two-step -0.854; J rejects 0.120
```

In this sample the two-step GMM interval misses the truth, while $\mathcal{S}_n(-1.6) = 3.06$ is far below the critical value of 9.49. The limit of the S statistic at infinity, 14.01, is above the critical value, so the S confidence set $[-7.04, 1.29]$ is bounded; the endpoints are found on a grid over $[-12, 8]$.

### 6.7 The Grand Reconciliation

| Block | Estimator | Estimate | 95% interval | Target and true value | Covered |
|---|---|---|---|---|---|
| G1 | MD (inverse-variance weighted) | $-1.535$ | $[-1.662, -1.409]$ | $-1.6$ | Yes |
| G1 | Equal weights | $-1.528$ | $[-1.757, -1.299]$ | $-1.6$ | Yes |
| G1 | Two-step GMM (household data) | $-1.535155$ | | Difference from MD $-2.0 \times 10^{-5}$ | |
| G1h | MD | $-1.406$ | $[-1.532, -1.280]$ | $\theta_*$ $-1.395$ | Yes |
| G1h | MD | $-1.406$ | $[-1.532, -1.280]$ | Simple average $-1.6$ | No |
| G1h | Equal weights | $-1.538$ | $[-1.767, -1.308]$ | Simple average $-1.6$ | Yes |
| G2 | OLS | $-0.318$ | | Probability limit $-0.382$ | |
| G2 | 2SLS | $-1.843$ | $[-2.654, -1.032]$ | $-1.6$ | Yes |
| G2 | Identity weight | $-1.307$ | $[-2.663, 0.050]$ | $-1.6$ | Yes |
| G2 | Two-step GMM, conventional | $-2.012$ | $[-2.713, -1.311]$ | $-1.6$ | Yes |
| G2 | Two-step GMM, Windmeijer | $-2.012$ | $[-2.831, -1.194]$ | $-1.6$ | Yes |
| G2 | Iterated GMM | $-2.043$ | $[-2.763, -1.323]$ | $-1.6$ | Yes |
| G2 | CUE | $-2.304$ | $[-3.053, -1.556]$ | $-1.6$ | Yes |
| G2 | S confidence set | | $[-4.09, -0.995]$ | $-1.6$ | Yes |
| G3 | Two-step GMM, conventional | $-0.229$ | $[-1.482, 1.025]$ | $-1.6$ | No |
| G3 | Two-step GMM, Windmeijer | $-0.229$ | $[-1.485, 1.027]$ | $-1.6$ | No |
| G3 | 2SLS | $-0.217$ | $[-1.483, 1.049]$ | $-1.6$ | No |
| G3 | CUE | $-0.362$ | $[-1.657, 0.932]$ | $-1.6$ | Yes |
| G3 | S confidence set | | $[-7.04, 1.29]$ | $-1.6$ | Yes |

The G1h row misses the simple average by design: MD aims at the precision-weighted $\theta_*$. G3's three Wald intervals do not cover the truth. A miss in a single sample is not evidence in itself; the systematic failure of conventional inference under weak moments shows in the Monte Carlo coverage of 0.774. In the same sample, the S confidence set and the CUE interval cover the truth.

The reconciliation of this section can be summed up as follows: when the moment conditions hold and identification is strong enough, every weight aims at the truth, the efficient weight is the most precise, and its standard errors need the Windmeijer correction; when the moment conditions are mutually incompatible, the estimator aims at a pseudo-true value defined by the weights; when the moment conditions are very weak, only the S confidence set gives an honest interval.

## 7 Failure Modes and Robustness

The exclusion restriction is the most fundamental assumption of GMM on experimental data, and randomization does not guarantee it. The premium-tutor upcharge may change parents' perception of quality, and the weekend surcharge may shift lessons from weekends to weekdays without changing the total; either would let an experiment bypass price and affect the outcome. The J test can catch some of these problems, but only those that show up as contradictions among the moment conditions: if each experiment's direct effect happens to bear the same proportion to its push on price, the ten ratios all shift to the same wrong value, and J has no power; and when J does reject, it does not tell you which moment is wrong either. There are three practical responses: first think through what each experiment changes besides price, and set the suspect experiments apart; report the sensitivity of Section 4.2, so that readers can see which moments the estimate depends on most; and drop the experiments one at a time and re-estimate, to see whether the estimate is stable.

When J rejects, first ask whether the cause is misspecification or heterogeneity. In G1h no test was faulty; the rejection comes from the restriction "there is only one elasticity." It is then fine to go on reporting an MD estimate, but you must make clear that it aims at the precision-weighted $\theta_*$, with weights determined by the number of households and the price difference of each city's test, and not at any kind of "average city." If the decision needs the simple average or the household-weighted average, use the corresponding weights and accept a larger standard error; if heterogeneity itself is of interest, report the city elasticities one by one or, as Meager (2019) does in combining seven microcredit experiments, use a hierarchical model to compromise between complete pooling and no pooling.

Small samples and many moments are a common trap for two-step GMM. When the number of moment conditions is not small relative to the sample size, conventional standard errors are too small and point estimates are biased toward OLS. The responses: with linear moment conditions, report Windmeijer standard errors by default; also report CUE or iterated GMM as a check; and do not add moment conditions without limit for the sake of "efficiency," since the dynamic panels of Chapter 11 will show that the J test also breaks down when there are too many moment conditions.

Weak moments are the most insidious failure. G3's Wald intervals look normal in width, yet they systematically miss the truth. Diagnosis starts with the first stage: do the experiments push the price hard enough? G3's first-stage F is only 3.61. For inference, switch to the S confidence set and accept that it may be very wide or even unbounded; as a remedy, the first instinct should not be to switch to a more complex estimator but to return to the design: enlarge the price differences of the experiments, or switch to a set of experiments that really move the price.

The moment conditions of this chapter all assume independence across observations. If regions are linked, for example because regions in the same city share advertising and tutors, or because the same households appear repeatedly over several weeks, $\hat S$ must be computed by cluster, and J and the S statistic must also use cluster-robust variances; Chapter 8 covers the practice of clustering. In MD, $\hat\Sigma_\pi$ must then account for the correlation among the reduced-form estimates.

MD has risks of its own. $\hat\Sigma_\pi$ must be the full variance of $\hat\pi$: when the reduced-form estimates come from the same sample, the covariances among them cannot be left out; when the reduced form is high-dimensional and the sample is small, $\hat\Sigma_\pi$ is poorly estimated, and the optimal weight turns out to be less robust than diagonal or equal weights, in which case one can compute the point estimate with a simple weight first and then compute standard errors with the sandwich.

Put these failure modes together, and the robustness of this chapter's methods comes down to three things: whether the moment conditions really hold, especially the exclusion restriction; whether the moment conditions are mutually compatible, and if not, whose average the estimator aims at; and whether the moment conditions are strong enough. The first rests on economic argument, the second on J and a clear understanding of the weights, and the third on first-stage diagnostics and the S confidence set.

## 8 Further Reading

::: {.readings}
Required reading, in suggested reading order:

- Wooldridge (2010, *Econometric Analysis of Cross Section and Panel Data*, 2nd ed.) §14.1 to §14.6. The consistency and asymptotic normality of GMM (Theorems 14.1 and 14.2), the efficient weight and the J statistic, the efficiency framework (Theorem 14.3), classical minimum distance, and the MD formulation of Chamberlain's approach.
- Hansen (1982, *Econometrica*). The founding work of GMM; focus on the efficient weight and the construction of the overidentification test.
- Stock and Wright (2000, *Econometrica*). The asymptotic theory of GMM under weak identification, and the S statistic.
- Chamberlain (1982, *Journal of Econometrics*). The $\pi$ matrix of panel data and minimum distance.

Further reading:

- Hansen, Heaton and Yaron (1996, *Journal of Business & Economic Statistics*). Proposes CUE and compares the finite-sample performance of three GMM estimators: two-step, iterated, and CUE.
- Windmeijer (2005, *Journal of Econometrics*). The finite-sample correction of the standard errors of linear two-step GMM.
- Newey and Windmeijer (2009, *Econometrica*). The bias of two-step GMM and the variance correction for CUE under many weak moments.
- Andrews, Gentzkow and Shapiro (2017, *Quarterly Journal of Economics*). The sensitivity $\Lambda_0$: which moments an estimate depends on.
- Angrist and Newey (1991, *Journal of Business & Economic Statistics*). Overidentification tests in fixed effects models of earnings equations; Wooldridge (2010) p. 551 holds that the asymptotic equivalence of MD and GMM in the Chamberlain model can be derived from it.
- Hansen (2014) §7.4 to §7.6, §8.12 to §8.13, and §13.2 to §13.6. Minimum distance under linear restrictions, MD tests, and another textbook treatment of GMM.
- Chamberlain (1984, *Handbook of Econometrics* Vol. 2). A survey of panel data methods, with empirical examples of MD.
- McFadden (1989), Pakes and Pollard (1989), and Gourieroux, Monfort and Renault (1993). SMM and indirect inference.
- Meager (2019, *American Economic Journal: Applied Economics*). Combines several randomized experiments with a Bayesian hierarchical model, a counterpoint to fixed-effect meta-analysis.
- Dubé and Misra (2023, *Journal of Political Economy*). Randomized price experiments on a platform, demand estimation, and personalized pricing.
:::

::: {.apa-refs}
- Andrews, I., Gentzkow, M., & Shapiro, J. M. (2017). Measuring the sensitivity of parameter estimates to estimation moments. *The Quarterly Journal of Economics, 132*(4), 1553-1592. https://doi.org/10.1093/qje/qjx023
- Angrist, J. D., & Newey, W. K. (1991). Over-identification tests in earnings functions with fixed effects. *Journal of Business & Economic Statistics, 9*(3), 317-323. https://doi.org/10.1080/07350015.1991.10509857
- Chamberlain, G. (1982). Multivariate regression models for panel data. *Journal of Econometrics, 18*(1), 5-46. https://doi.org/10.1016/0304-4076(82)90094-X
- Chamberlain, G. (1984). Panel data. In Z. Griliches & M. D. Intriligator (Eds.), *Handbook of econometrics* (Vol. 2, pp. 1247-1318). North-Holland. https://doi.org/10.1016/S1573-4412(84)02014-6
- Dubé, J.-P., & Misra, S. (2023). Personalized pricing and consumer welfare. *Journal of Political Economy, 131*(1), 131-189. https://doi.org/10.1086/720793
- Gourieroux, C., Monfort, A., & Renault, E. (1993). Indirect inference. *Journal of Applied Econometrics, 8*(S1), S85-S118. https://doi.org/10.1002/jae.3950080507
- Hansen, B. E. (2014). *Econometrics* [Unpublished manuscript]. Department of Economics, University of Wisconsin.
- Hansen, L. P. (1982). Large sample properties of generalized method of moments estimators. *Econometrica, 50*(4), 1029-1054. https://doi.org/10.2307/1912775
- Hansen, L. P., Heaton, J., & Yaron, A. (1996). Finite-sample properties of some alternative GMM estimators. *Journal of Business & Economic Statistics, 14*(3), 262-280. https://doi.org/10.1080/07350015.1996.10524656
- McFadden, D. (1989). A method of simulated moments for estimation of discrete response models without numerical integration. *Econometrica, 57*(5), 995-1026. https://doi.org/10.2307/1913621
- Meager, R. (2019). Understanding the average impact of microcredit expansions: A Bayesian hierarchical analysis of seven randomized experiments. *American Economic Journal: Applied Economics, 11*(1), 57-91. https://doi.org/10.1257/app.20170299
- Newey, W. K., & Windmeijer, F. (2009). Generalized method of moments with many weak moment conditions. *Econometrica, 77*(3), 687-719. https://doi.org/10.3982/ECTA6224
- Pakes, A., & Pollard, D. (1989). Simulation and the asymptotics of optimization estimators. *Econometrica, 57*(5), 1027-1057. https://doi.org/10.2307/1913622
- Stock, J. H., & Wright, J. H. (2000). GMM with weak identification. *Econometrica, 68*(5), 1055-1096. https://doi.org/10.1111/1468-0262.00151
- Windmeijer, F. (2005). A finite sample correction for the variance of linear efficient two-step GMM estimators. *Journal of Econometrics, 126*(1), 25-51. https://doi.org/10.1016/j.jeconom.2004.02.005
- Wooldridge, J. M. (2010). *Econometric analysis of cross section and panel data* (2nd ed.). MIT Press.
:::
