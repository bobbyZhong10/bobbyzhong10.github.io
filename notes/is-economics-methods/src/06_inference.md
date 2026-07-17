---
title: "Inference: Delta Method, Bootstrap, Clustering & Testing"
subtitle: "Honest Standard Errors"
seriesline: "Foundations of Information Systems Economics · Chapter 6"
series: "Foundations of Information Systems Economics: Empirical Methods for Information Systems Economics"
sidebartitle: "Foundations of Information Systems Economics"
sidebarsub: "Theory · Methods · Applications"
chapterlabel: "Chapter 6 · Inference: Delta Method, Bootstrap, Clustering & Testing"
---

## Introduction

Halcyon's city experiment delivered a result so clean it was almost reassuring: $t=13.76$. The point estimate is fine and the regression converged without complaint. The trouble is that once you treat orders that are correlated within a city as if they were independent observations, the sample looks like it carries far more information than it actually does. Cluster by city, and the $t$ value falls to 3.50. Worse still, in a placebo experiment where the true effect is exactly zero, that seemingly precise iid standard error produces a "significant" result about two thirds of the time.

Inference errors have a dangerous kind of politeness: they leave the point estimate untouched and quietly shrink the number in parentheses. Nonlinear quantities require the Delta method or the bootstrap to propagate parameter uncertainty; correlation within cities, firms, or users demands cluster-robust variances, and correlation over time demands the corresponding long-run variance estimators; and when there are only a few clusters, when estimation happens in two steps, or when a test is reparameterized, the familiar asymptotic approximations keep breaking down. The fact that Wald, likelihood ratio, and Lagrange multiplier tests all converge to the same thing in large samples is no guarantee that they will agree in the finite sample sitting in front of you.

This chapter is organized around one question: what does a standard error actually promise? Halcyon's willingness-to-pay ratio, its city experiment, and a time series will each show, in turn, that a variance is not decoration that the regression command hands you for free, but part of the research design and the dependence structure of the data. A credible empirical result does not need the largest possible $t$ value; it needs a $t$ value that survives interpretation even after you have fully acknowledged the uncertainty.

## 1 A t Value of 13.76

::: {.case}
Halcyon is a fictional food-delivery platform. We tell this case with simulated data, which buys us something real data can never offer as a teaching device: the data-generating process is fully known, so every standard error can be reconciled against the truth. The platform operates in 40 cities and wants to evaluate the causal effect of a city-level operating policy (say, a new delivery-pricing scheme) on how users place orders. The policy launches in a randomly chosen half of the cities, with the other half as controls, and the outcome is measured at the level of the user session, with about 200 sessions per city. Management's question is direct: how much did the policy raise the outcome metric, and can we believe the number?

The platform's analyst ran the most direct regression, regressing the session-level outcome $y$ on an indicator $D$ for whether the city launched the policy, $y_{ic} = \mu + \tau D_c + \varepsilon_{ic}$, where $i$ indexes the session and $c$ the city. The estimated treatment effect is $\hat\tau = 0.318$. She computed the standard error the regression reports by default: the SE is 0.0231, so the $t$ value comes out to $0.318 / 0.0231 = 13.76$, with a $p$ value of $4 \times 10^{-43}$. That number says that if the policy truly had no effect, the chance of seeing a $\hat\tau$ this large would be unimaginably small. She was ready to write the ironclad conclusion into her report.
:::

That 13.76 is the first bad number of this chapter, and it is bad in a way that runs deeper than plain absurdity: the arithmetic is not wrong, the regression coefficient is not biased, and $\hat\tau = 0.318$ is a reasonable estimate of the true value 0.25. What is wrong is that 0.0231 standard error, which rests on an assumption that simply does not hold here: 40 cities $\times$ 200 sessions $=$ 8000 sessions, each one a separate piece of information.

They are not separate. The 200 sessions inside a single city share everything about that city in that period: local demand highs and lows, the weather, what competitors are doing. Once you count those common shocks, sessions within a city are strongly positively correlated. Positive correlation means redundant information: the 200th session in a city does not bring a fresh piece of information the way the 1st one did, it mostly just repeats the city's common state. So the effective information contained in 8000 sessions is far less than 8000 independent observations, and much closer to 40 cities. The iid standard error treats each session as an independent vote, which amounts to conjuring thousands of extra votes out of thin air, and naturally computes the uncertainty as far too small.

Getting the correlation into the standard error correctly is the job of the cluster-robust variance, which Section 3 derives. On Halcyon it gives a standard error of 0.0908, nearly four times the iid 0.0231, and the $t$ value collapses accordingly from 13.76 to $0.318 / 0.0908 = 3.50$. The effect is still significant, but it has dropped from "unbelievably certain" to "reasonably solid," and those are entirely different things to report.

::: {.warning}
There is a misconception worth puncturing on the spot. Many people believe that swapping the standard error for a heteroskedasticity-robust one (the so-called robust or White standard error) fixes everything. On Halcyon, the heteroskedasticity-robust HC1 standard error is 0.0231, identical to the iid one, with no improvement whatsoever. The reason was planted in the previous chapter: heteroskedasticity-robust standard errors correct only for the single fact that observations have unequal variances, and they still assume zero correlation between different observations. The disease here is not unequal variance, it is correlation among sessions in the same city. Robust standard errors and cluster standard errors are two different things, one handles variance and the other handles correlation, and treating the former as a cure-all is one of the most common errors in applied work.
:::

The shrinkage of one estimate's standard error is only the tip of the iceberg. What is truly alarming is a health check on the reliability of this iid inference. Suppose the policy in fact does nothing ($\tau = 0$). We repeatedly generate placebo data of this kind and each time use the iid standard error to test the correct null hypothesis "the effect is zero" at the 5% level. An honest 5% test should reject it by mistake only 5% of the time. Yet the actual rejection rate of the iid inference is 64.5%: two thirds of the time it "discovers" a significant effect that does not exist, out of pure noise. With cluster-robust standard errors this false-positive rate returns to 7.0%, close to the 5% it should be. This is the lesson the chapter hammers on: when the standard error is wrong, what suffers is not the precision of any single estimate but the credibility of your entire significance verdict, and it will show you a whole sky of signal where there is none.

## 2 Economic Model and Estimand

Inference does not change what you estimate, it only cares about how the thing you estimated fluctuates around the truth. So this section first picks up the deliverables of the previous chapter, then makes clear which quantities this chapter will attach standard errors to.

The previous chapter delivered a batch of consistent point estimates together with their sandwich variances. For an extremum estimator $\hat\theta$, we already know that $\sqrt{n}(\hat\theta - \theta_0)$ is asymptotically normal with variance $A^{-1} B A^{-1}$, where $A$ is the curvature of the objective function and $B$ is the variance of the gradient. The standard error of one coefficient $\hat\theta_k$ is the square root of the corresponding diagonal element of this asymptotic variance matrix divided by $\sqrt n$, and the confidence interval is $\hat\theta_k \pm 1.96\, \mathrm{SE}$. This machinery is already enough for a single coefficient. What this chapter handles are the three kinds of estimand for which it is not enough.

The first kind is a nonlinear function of the coefficients. Besides the city experiment, Halcyon also ran an analysis of delivery pricing: a user's satisfaction $sat$ with an order falls with the delivery time $time$ and the delivery fee $fee$, $sat = \beta_0 + \beta_{time}\, time + \beta_{fee}\, fee + \varepsilon$. Management wants to know how users value time, that is, how much extra delivery fee a user is willing to pay to save one minute of waiting. This quantity is a ratio of two coefficients, $\mathrm{WTP} = \beta_{time} / \beta_{fee}$, in units of currency per minute. The regression gives the variance of $\beta_{time}$ and $\beta_{fee}$ each, but the variance of the ratio is not a simple combination of them, and this is exactly what the Delta method and the bootstrap are for.

The second kind is a treatment effect on correlated data. The $\tau$ in the city experiment of Section 1 has a fine point estimate, but its standard error must fold in the correlation among same-city sessions, or else you get a mirage like the 13.76.

The third kind is a hypothesis itself. Often we want not just an interval for a coefficient but a test of a specific proposition: does delivery time affect repurchase at all, are two coefficients jointly zero. There are several ways to construct a test, and whether they agree or disagree is itself a question at the level of the estimand.

Before working through these three kinds, we need to think clearly about what a standard error actually promises, otherwise all the corrections that follow will read like a pile of technical detail. That is the task of the opening of the next section. The key points of this section can be summarized as follows: inference does not change the point estimate, it only characterizes its sampling fluctuation; this chapter attaches honest uncertainty to three kinds of quantity, namely a nonlinear function of coefficients (the ratio WTP), a treatment effect on correlated data (the $\tau$ of the city experiment), and hypothesis testing itself.

## 3 The Theory of the Standard Error: Delta, Bootstrap, and Correlation

This is the most finely analyzed section of the chapter. It answers a question that looks simple but is in fact deep: what does a standard error actually say, and how do you get it right in the three common situations? We first make clear what a standard error promises (3.1), then cover the two routes for a nonlinear function, the Delta method (3.2) and the bootstrap (3.3), then how correlation shatters the iid standard error and how the cluster-robust variance fixes it (3.4), and finally the counterpart in time series, HAC (3.5).

### 3.1 What a Standard Error Promises

Let us first take apart this word we say every day. Imagine we could redo Halcyon's entire experiment a thousand times under exactly the same conditions, each time obtaining a slightly different $\hat\tau$ because of randomness. These thousand values of $\hat\tau$ form a distribution called the sampling distribution, which has a center (near the truth, if the estimator is consistent) and a width. The standard error is an estimate of that width, the standard deviation of the sampling distribution. It is not a statement about "whether this particular estimate is right," but about "how spread out the estimate would be if we redid it many times."

::: {.intuition}
The promise of a 95% confidence interval can be understood this way, and it is often misread, so it is worth stating precisely. It does not say "there is a 95% probability that the truth lies in this interval," because the truth is a fixed number, either inside or not. It speaks about the procedure that builds the interval: if you redid the experiment many times, and built an interval each time by the same rules, about 95% of these intervals would cover the truth. The interval is random, the truth is fixed, and what is being resampled is the interval, not the truth. This "95% of the intervals cover the truth" is called coverage, and it is the gold standard for whether a standard error is honest. The disease of the iid inference in Section 1, put in the language of coverage, is that its nominally 95% intervals actually cover the truth far less than 95% of the time, and Section 6 will measure several such coverage rates.
:::

The previous chapter already gave the source of the standard error, namely the sandwich variance. But inside that sandwich there is one inconspicuous assumption that decides everything, and this chapter watches it specifically: the variance of the gradient $B = \mathbb{E}[\nabla q\, \nabla q']$ is estimated by $\frac{1}{n}\sum_i \nabla q_i\, \nabla q_i'$ under the assumption that different observations are independent. That sum adds each observation's contribution as if it were an independent piece of information. When observations are correlated, this estimate is wrong, and the direction of the error is almost always downward, which is the subject of Section 3.4. And when the quantity we want is a nonlinear function of the estimate, even getting the variance of $g(\hat\theta)$ from the variance of $\hat\theta$ becomes a problem, which is the subject of the next two subsections.

### 3.2 The Delta Method: A Variance for a Nonlinear Function

We have $\hat\theta$ and its variance matrix $\hat V$, and we want the variance of $g(\hat\theta)$, where $g$ is a nonlinear function, for example the ratio $\beta_{time}/\beta_{fee}$. The most direct idea is to straighten $g$ into a line near the truth, because a line is the only function that lets the "variance of the input" translate simply into the "variance of the output." This is the Delta method.

Linearize $g$ at $\theta_0$ with a first-order Taylor expansion, $g(\hat\theta) \approx g(\theta_0) + \nabla g(\theta_0)'\,(\hat\theta - \theta_0)$. The second term on the right is a known vector $\nabla g$ times an asymptotically normal random vector $\hat\theta - \theta_0$, and the variance of a linear combination is standard, so:

::: {.theorem}
**(Delta method)** If $\sqrt{n}(\hat\theta - \theta_0) \xrightarrow{d} \mathcal{N}(0, \Sigma)$, and $g$ is continuously differentiable at $\theta_0$ with gradient $\nabla g(\theta_0) \neq 0$, then

$$\sqrt{n}\big(g(\hat\theta) - g(\theta_0)\big) \xrightarrow{d} \mathcal{N}\big(0,\ \nabla g(\theta_0)'\, \Sigma\, \nabla g(\theta_0)\big).$$
:::

The intuition of this formula lives entirely in the gradient $\nabla g$. It is an amplification factor: the steeper $g$ is at the truth (the larger $\nabla g$), the more the same amount of jitter in the input $\hat\theta$ gets amplified by $g$ into jitter in the output. The variance formula $\nabla g'\Sigma \nabla g$ just sandwiches the input variance matrix $\Sigma$ between this amplification factor front and back, another sandwich shape.

Applied to a ratio it is very concrete. $\mathrm{WTP} = g(\beta_{time}, \beta_{fee}) = \beta_{time}/\beta_{fee}$, with gradient $\nabla g = \big(1/\beta_{fee},\ -\beta_{time}/\beta_{fee}^2\big)'$. Substituting into the formula, the variance of the ratio is

$$\mathrm{Var}(\mathrm{WTP}) \approx \frac{1}{\beta_{fee}^2}\, V_{tt} - 2\,\frac{\beta_{time}}{\beta_{fee}^3}\, V_{tf} + \frac{\beta_{time}^2}{\beta_{fee}^4}\, V_{ff},$$

where $V_{tt}, V_{ff}, V_{tf}$ are the variances and covariance of the two coefficients. Section 6 will show that when the delivery fee varies enough in the data so that $\beta_{fee}$ is estimated precisely, the standard error from this formula agrees almost exactly with the true sampling fluctuation.

::: {.warning}
The Delta method has a fatal weakness hidden right inside that gradient. Look at the ratio-variance formula: the denominator carries $\beta_{fee}^2$ and $\beta_{fee}^4$ (the gradient itself carries $\beta_{fee}$ and $\beta_{fee}^2$). When $\beta_{fee}$ is close to zero, that is, when the delivery fee does not really affect satisfaction and the coefficient is estimated very imprecisely, the gradient explodes and the linearization is completely distorted. The more fundamental problem is that in this case the true sampling distribution of the ratio $\beta_{time}/\beta_{fee}$ is severely skewed and may not even have a finite variance (the denominator occasionally grazes zero, and the ratio flies off to infinity), while the Delta method mechanically spits out a symmetric normal interval that erases the skew entirely. This "near-zero denominator" disease is called the Fieller problem, and it is the same ghost as the weak identification we will meet later in this series when we cover instrumental variables, and as the weak-denominator variant in Section 6. Moreover, the Delta method also fails when $g$ is not differentiable at the expansion point (as with an absolute value) or when the gradient happens to be zero. In these cases you need a method that does not rely on linearization, and that is the bootstrap.
:::

### 3.3 The Bootstrap: Simulating the Sampling Distribution Directly

The Delta method relies on an analytic linear approximation. The bootstrap takes an entirely different road: it does not approximate, it rebuilds the sampling distribution itself on the computer. The idea looks like a magic trick at first and turns out to be startlingly plain on reflection.

The difficulty is that the sampling distribution normally requires "redoing the experiment many times" to be observed, yet we have only one dataset. The bootstrap's core move is to use the empirical distribution of the sample in hand to impersonate the population: since the real population is out of reach, assume the population looks exactly like the sample, with each observation carrying $1/n$ of the probability mass. To "resample" from this impersonated population is to draw $n$ observations with replacement from the original sample, forming one bootstrap sample. Drawing with replacement is the key, because it makes each bootstrap sample slightly different from the original (some observations are drawn more than once, some not at all), which reproduces exactly the randomness of "redoing the experiment once."

::: {.intuition}
Think of the bootstrap as a photocopier. It has no population as the original printing press, so it treats the one sample in hand as the master copy and, by rearranging the observations on that master with replacement, prints out thousands of "parallel samples." Each parallel sample recomputes the statistic you care about, and the spread of those thousands of statistics is a portrait of the sampling distribution. Its magic is that it needs to know nothing about the analytic form of the statistic: whether what you compute is a ratio, a quantile, or a complex object that first estimates game parameters and then solves for a Nash equilibrium, the bootstrap treats them all the same and simply recomputes. This is also why it is so popular in structural estimation, where many quantities have no closed-form variance at all.
:::

The algorithm is just three steps: draw a bootstrap sample of the same size with replacement from the original sample, recompute the statistic $\hat\theta^*$ on it, and repeat $B$ times ($B$ is often a few thousand). The spread of these $B$ values of $\hat\theta^*$ is the estimate of the sampling distribution. With it in hand, there are several ways to construct a confidence interval, increasing in refinement. The most naive is the percentile interval, which simply takes the 2.5% and 97.5% quantiles of these $B$ values. The second is the basic or pivotal interval, which reflects around $\hat\theta$ to correct for the shift of the bootstrap distribution relative to $\hat\theta$. The most careful is the percentile-t (studentized) interval: for each bootstrap sample it computes not just $\hat\theta^*$ but also its own standard error, constructs a studentized $t^* = (\hat\theta^* - \hat\theta)/\mathrm{SE}^*$, and then uses the quantiles of $t^*$ to build the interval.

::: {.theorem}
**(Higher-order refinement of the bootstrap)** For an asymptotically normal, studentized statistic, and for a one-sided confidence interval, the coverage error of the percentile-t bootstrap is $O(n^{-1})$, which approaches the nominal level faster than the $O(n^{-1/2})$ of the first-order normal approximation (and of the unstudentized percentile interval). This section constructs equal-tailed two-sided intervals, for which the $O(n^{-1/2})$ skewness errors at the two ends cancel each other, so the coverage error of the normal approximation is already $O(n^{-1})$, and the equal-tailed percentile-t is likewise $O(n^{-1})$ but with a better constant; if instead you use a symmetric two-sided interval (taking the quantiles of $|t^*|$), the $O(n^{-1})$ term of the percentile-t also cancels and the coverage error drops further to $O(n^{-2})$.
:::

This asymptotic refinement is the most charming property of the bootstrap, and it deserves a sentence of intuition. The first-order normal approximation uses only the first two moments of the sampling distribution (the mean and variance) and treats it as a symmetric normal. But the true sampling distribution often has skew (a nonzero third moment), and the studentized bootstrap, because it resamples nonparametrically, automatically captures the skew term (mathematically the $n^{-1/2}$ order in the Edgeworth expansion), so its interval can lean toward the direction it should lean and its coverage converges to 95% faster. The weak-denominator example in Section 6 will demonstrate this very plainly: when the ratio's distribution is severely right-skewed, both the symmetric Delta interval and the percentile interval fit worse than the studentized interval.

The bootstrap has several variants that must be kept apart. The default, the pairs bootstrap, resamples whole rows $(y_i, x_i)$, is the most robust, and is immune to heteroskedasticity. The residual bootstrap holds $x$ fixed and resamples only the estimated residuals, is slightly more efficient but relies on homoskedasticity. The wild bootstrap multiplies each residual by a random sign (for example Rademacher weights) and is designed for heteroskedasticity, and it has a clustered version, the wild cluster bootstrap, which is the workhorse for the few-cluster problem in Section 4.

::: {.warning}
The bootstrap is not omnipotent, and its weakness is precisely its foothold. First, it relies on "the empirical distribution of the sample approximating the population," so it requires the sample to be large enough and to not systematically miss some part of the population, and it fails in small samples. Second, it assumes by default that observations are independent and identically distributed (resampling with replacement is a simulation of iid), and once the data are correlated (time series, clusters) the naive bootstrap will, like the iid standard error, understate the uncertainty, so you must switch to the block bootstrap (which resamples in whole blocks to preserve within-block correlation) or the wild cluster bootstrap. Third, for statistics that are not smooth or not asymptotically normal (such as a parameter on the boundary of the parameter space, or an extremum statistic like a maximum), the bootstrap fails, and here its relative, subsampling (drawing smaller samples without replacement), may work instead. In one sentence, the bootstrap photocopies an iid sample, so if the sample is not iid, or the statistic is not smooth, the copy comes out distorted.
:::

### 3.4 Correlation Shatters the iid Standard Error: The Cluster-Robust Variance

Now back to the root of the 13.76 disease in Section 1. Where the iid standard error goes wrong, and how to fix it, this subsection derives cleanly.

The problem is in the middle filling $B$ of the sandwich variance, that is, the variance of $\sqrt n$ times the sample gradient (for OLS this is $\frac{1}{\sqrt n}\sum_i x_i \varepsilon_i$). Under the iid assumption, the variance of this sum equals the sum of the variances of its terms, $\frac{1}{n}\sum_i x_i x_i' \varepsilon_i^2$, because the variance of a sum of independent random variables has no cross-covariance terms. But when observations are correlated in clusters by city, the within-cluster cross-covariance terms $\mathbb{E}[x_i x_j' \varepsilon_i \varepsilon_j]$ (with $i, j$ in the same city) are nonzero, and dropping them is exactly the understatement. What the cluster-robust variance does is add these within-cluster cross terms back.

::: {.theorem}
**(Cluster-robust variance)** Let the observations be split into $G$ clusters,^[This chapter uses $g$ and $G$ to denote a generic cluster index and cluster count, valid only within this chapter, differing from the convention elsewhere in the book where $g$ denotes a cohort (a treatment-group cohort).] independent across clusters and arbitrarily correlated within a cluster. Stack the regressors and residuals of the $g$-th cluster into $X_g$ and $\hat u_g$. The cluster-robust variance of OLS is

$$\hat V^{CR} = (X'X)^{-1} \left( \sum_{g=1}^{G} X_g'\, \hat u_g\, \hat u_g'\, X_g \right) (X'X)^{-1} = (X'X)^{-1} \left( \sum_{g=1}^{G} \Big(\sum_{i \in g} x_i \hat u_i\Big)\Big(\sum_{i \in g} x_i \hat u_i\Big)' \right)(X'X)^{-1}.$$
:::

Comparing the two makes the difference clear. The heteroskedasticity-robust (HC) filling is $\sum_i x_i x_i' \hat u_i^2$, which adds up only each observation's own square, equivalent to assuming the within-cluster cross terms are zero. The cluster-robust filling is $\sum_g (\sum_{i\in g} x_i \hat u_i)(\sum_{i\in g} x_i \hat u_i)'$, which first sums $x_i \hat u_i$ within each cluster and then takes the outer product, and this step brings in all the $i \neq j$ cross-covariances within the cluster. When there is positive within-cluster correlation, these cross terms are positive, the filling grows, and the standard error grows with it. This is the mechanism behind Halcyon's standard error rising from 0.0231 to 0.0908. The reason the heteroskedasticity-robust HC1 stays stubbornly at 0.0231 is precisely that its filling contains none of these cross terms.

::: {.intuition}
How much information within-cluster correlation eats up has a beautiful rule of thumb called the Moulton factor. When the treatment variable varies at the cluster level, the within-cluster correlation coefficient is $\rho$, and each cluster has on average $\bar n$ observations, the inflation of the correct standard error relative to the iid standard error is about

$$\sqrt{1 + (\bar n - 1)\,\rho}.$$

The intuition is "effective sample size." When $\rho = 0$ (no within-cluster correlation) the factor is 1, and $n$ observations are $n$ pieces of information. When $\rho = 1$ (identical within the cluster) the factor is $\sqrt{\bar n}$, and a whole cluster of $\bar n$ observations is really worth only one piece of information, so the effective sample size shrinks to the number of clusters. Reality is somewhere in between. On Halcyon the within-cluster correlation is about 0.083 and each cluster has 200 sessions, so the Moulton factor predicts an inflation of $\sqrt{1 + 199 \times 0.083} = 4.18$, quite close to the 3.93 actually observed. Note that the $\bar n$ component in this factor is large: even if $\rho$ is a mere 0.083, multiplying by 199 is enough to quadruple the standard error. This explains why "my within-cluster correlation is tiny, so clustering probably does not matter" is a dangerous self-consolation, because as long as each cluster has many observations, even a small correlation gets amplified into a large bias.
:::

A practical issue with clustering is "at which level to cluster." The principle is simple and also hard: cluster at the level where you believe the judgment "independent across clusters, correlated within" holds. In the city experiment, treatment is assigned at the city level, so you should cluster on city; if you also worry about regional shocks across cities, you might cluster at the coarser regional level. Clustering too finely misses correlation (repeating the iid mistake), while clustering too coarsely leaves too few clusters (the trouble of the next section). A common rule of thumb is that if the model includes fixed effects at some level, you usually need to cluster at least at that level.

### 3.5 HAC: The Same Thing Along the Time Dimension

Clustering handles clustered correlation in the cross section. Switch to a time series, and the same disease appears in the guise of serial correlation, with a cure cut from the same cloth. Daily metrics like a platform's daily active users or gross merchandise volume are strongly correlated between today and yesterday, so the regression residuals are serially correlated, and the iid standard error once again understates by treating correlated observations as independent.

The fix is the heteroskedasticity and autocorrelation consistent (HAC) variance, most commonly the Newey-West estimator. On top of the heteroskedasticity term, its filling adds the covariance between residuals and their lags, giving smaller weight to more distant lags:

$$\hat\Omega^{HAC} = \sum_{t} \hat u_t^2\, x_t x_t' + \sum_{\ell=1}^{L} w_\ell \sum_{t=\ell+1}^{T} \hat u_t \hat u_{t-\ell}\big(x_t x_{t-\ell}' + x_{t-\ell} x_t'\big),\qquad w_\ell = 1 - \frac{\ell}{L+1}.$$

The first term is the heteroskedasticity-robust filling, and the second brings in the time correlation up to lag $L$, where $w_\ell$ is the triangular weight given by the Bartlett kernel, letting the covariance contribution decay linearly to zero with lag distance. The bandwidth $L$ must be chosen: too small misses long-range correlation, too large brings in too much noisy covariance. This is the same idea as clustering in two incarnations, with within-cluster correlation replaced by time correlation, and the "cluster" replaced by "a window of time." Section 6 will show that on a series with AR(1) errors, Newey-West raises the standard error from the iid 0.0484 to 0.0613 and lifts the confidence interval's coverage from 84.9% to 92.0%, the right direction, but also honestly exposing the old complaint that HAC often does not fix things all the way in finite samples and still slightly undercovers.

The key points of this section can be summarized as follows: a standard error is an estimate of the width of the sampling distribution, and whether it is honest is checked by coverage; a nonlinear function of coefficients is handled either by linearizing with the Delta method (the ratio variance has a closed form, but distorts through skew when the denominator is near zero) or by simulating the sampling distribution directly with the bootstrap (percentile-t has a higher-order refinement, but requires iid and smoothness); the mechanism by which correlated data shatter the iid standard error is that the cross-covariances in the sandwich filling are dropped, the cluster-robust variance adds the within-cluster cross terms back and its inflation is characterized by the Moulton factor, and HAC is its counterpart along the time dimension.

## 4 Testing and Practice

Once the standard error is right, the next step is to use it for testing, and to handle a few practical headaches you cannot avoid. This section covers the trinity of hypothesis testing (4.1), the few-cluster correction (4.2), randomization inference (4.3), and standard errors for two-step estimation (4.4).

### 4.1 The Trinity of Hypothesis Testing

To test a hypothesis, say "delivery time does not affect repurchase" (some coefficient is zero), there are three classic constructions that approach the same problem from three different angles, asymptotically equivalent yet each with its own temperament.

The Wald test stands at the unrestricted-estimation end: first estimate the parameter freely, then see how far it is from the null value, standardizing that distance by the estimated variance. For a single restriction $\theta = \theta_0$, the Wald statistic is $W = (\hat\theta - \theta_0)^2 / \widehat{\mathrm{Var}}(\hat\theta)$, asymptotically $\chi^2_1$. The likelihood ratio (LR) test stands in between: estimate both the unrestricted and the restricted model and compare the difference in their log-likelihoods, $LR = 2(\ell_u - \ell_r)$. The Lagrange multiplier (LM, also called the score) test stands at the restricted end: estimate only the restricted model and see how far its score (the gradient of the log-likelihood) is from zero in the unrestricted direction, $LM = s(\hat\theta_r)'\, I(\hat\theta_r)^{-1}\, s(\hat\theta_r)$. All three are asymptotically $\chi^2$ with degrees of freedom equal to the number of restrictions.

::: {.intuition}
The geometry of the three tests can be remembered by climbing a hill. Picture the log-likelihood as a hill whose summit is the unrestricted MLE, and the null hypothesis as a specified point on the slope. Wald asks about horizontal distance: from the summit to that point, how far did you walk sideways (converted by the steepness at the summit). LM asks about slope: standing at that specified point, how steep is the ground under your feet and how much does it want to slide toward the summit (the steeper the slope, the less reasonable it is that this point is the summit). LR asks about height difference: how much higher is the summit than that point. The same hill, the same point, measured from three angles of distance, slope, and height, asymptotically give the same answer. They just need different information: Wald needs only the summit, LM needs only that point on the slope, LR needs both. Use whichever is easiest to compute, and that is the most practical thing about the trinity.
:::

Asymptotic equivalence does not mean finite-sample agreement. The three give different numbers in small samples, and Section 6 will show Wald, LR, and LM at 35.58, 36.93, and 36.47 respectively, very close but not identical. More worth guarding against is a distinctive flaw of Wald: it is not invariant to reparameterization. The same hypothesis, written as $\theta = 0.5$ or as $\log\theta = \log 0.5$, are two logically equivalent statements, yet Wald gives a different statistic and a different $p$ value, because it depends on the linearization of $g$ at the estimate (the Delta-method machinery again), and reparameterization changes the gradient. Section 6 will demonstrate that two ways of writing the same hypothesis give $p = 0.0011$ and $p = 0.0092$, an order of magnitude apart. The LR test has no such flaw, because it looks only at the likelihood value and not at how the parameters are written, so when in doubt LR is often more reliable. There is one further layer: if the variance Wald uses comes from a misspecified model, it distorts along with the standard error, and in that case you should use the previous chapter's sandwich robust variance to construct a robust Wald.

### 4.2 Few Clusters: When G Is Small

The good properties of the cluster-robust variance are asymptotic, and that asymptotics is in the sense of the number of clusters $G \to \infty$, not the number of observations $n \to \infty$. This is deadly: even if you have millions of observations, as long as they fall into very few clusters (say an experiment in only 10 cities), the cluster-robust standard error is still untrustworthy, systematically too small, and the test over-rejects. In the simulations of Section 6, with $G = 10$ the cluster test using normal critical values has a false-positive rate of 10.3%, clearly above 5%.

There are several escalating corrections. The cheapest is to switch the critical value from the normal to a $t$ distribution with $G - 1$ degrees of freedom, acknowledging that with few clusters the tails should be fatter. In Section 6 this step alone pushes the $G = 10$ false-positive rate from 10.3% down to 6.7%. A more careful layer is to use cluster variances with small-sample bias corrections like CR2 or CR3 (analogous to HC2 and HC3 for heteroskedasticity), which work better for a few large clusters. The most robust layer is the wild cluster bootstrap: re-estimate the model under the null (that is, set the coefficient being tested to zero), multiply each cluster's residuals by a random sign as a whole cluster to generate bootstrap samples, recompute the cluster-robust $t$, and use its bootstrap distribution rather than the normal to set the critical value. In one $G = 10$ example in Section 6, for the same $t = 2.21$, the normal gives $p = 0.027$ (looks significant), $t(9)$ gives 0.054, and the wild cluster bootstrap gives 0.073 (actually not significant), and the gap among the three is exactly the divide, under few clusters, between "trusting the normal" and "trusting the bootstrap."

### 4.3 Randomization Inference: A Route That Does Not Rely on Asymptotics

All the methods above rely on asymptotic normality. When the sample is small, or the clusters are few, or the statistic has a strange distribution, the asymptotic approximation is untrustworthy to begin with. Randomization inference (RI) offers an entirely different route: it does not ask what the sampling distribution looks like, it uses only the fact that treatment in the experiment was assigned at random.

Its logic is beautiful and thorough. Test a sharp null, for example "the effect of the policy on every single unit is exactly zero." If this null is true, then the outcome observed for each unit, whether it was actually assigned to treatment or control, would not change if the treatment assignment changed. So you can do one thing: hold the data fixed, randomly permute the treatment labels across units, and each time recompute a treatment-effect estimate. Collect the estimates corresponding to all possible (or many random) permutations, and you get the distribution of the estimate that arises purely from random assignment when the null is true, called the permutation distribution. If the actually observed estimate lands in the extreme tail of this distribution, reject the null. The $p$ value is the fraction of permutations with an estimate more extreme than the actual one.

The virtue of RI is that it is exact and relies on no large-sample approximation, holding no matter how small the sample, because what it tests is the design fact of randomization. It is naturally suited to experiments with few clusters: in the city experiment, permuting the treatment labels across cities (rather than across sessions) gives a test robust to the small number of cities. The price is that the null must be sharp (holding for every unit), a stronger null than one that only constrains the average effect; and for non-experimental data you also need a credible "assignment mechanism" to define how to permute. Later in this series, when we cover the design of digital experiments, we will develop RI together with design-based inference, and here you need only remember that it is an exact alternative outside asymptotic inference.

### 4.4 Two-Step Estimation and Generated Regressors

The last practical headache is two-step estimation. Many empirical pipelines first estimate something and then feed it as data into a second step, for example first estimating a propensity score or a predicted value and then using it in a regression. If the second step uses ordinary standard errors directly, it misses the uncertainty introduced by the first-step estimate itself and computes the SE as too small. This class of problem is collectively called the generated regressors problem.

The analytic fix is called the Murphy-Topel variance, which adds a term to the second-step variance to absorb the uncertainty transmitted from the first-step estimation error through the second step. But the analytic correction has to be rederived for each new two-step structure, which is tedious and error-prone. In modern practice the more convenient approach is to bootstrap the whole two-step pipeline as a single object: on each bootstrap sample, rerun the entire process from the first step to the second, so that the first step's sampling fluctuation is naturally folded into the bootstrap distribution of the second-step estimate. This is the same reasoning as Section 3.3's point that the bootstrap can handle complex objects that "first estimate game parameters and then solve for an equilibrium": as long as you can write the whole pipeline as a single function, the bootstrap can give it an honest standard error.

The key points of this section can be summarized as follows: Wald, LR, and LM test the same hypothesis from the three angles of distance, slope, and height, asymptotically equivalent but differing in finite samples, and Wald is additionally not invariant to reparameterization, so when in doubt use LR or a robust Wald; under few clusters the cluster variance over-rejects, corrected layer by layer with $t(G-1)$, CR2, and the wild cluster bootstrap; randomization inference uses the design fact of random assignment to give an exact test that does not rely on asymptotics; and two-step estimation should use Murphy-Topel or bootstrap the whole pipeline directly to fold in the first step's uncertainty.

## 5 Anchoring Papers

Methods only stand up when they land in real research. Three anchoring papers: one turned the mechanics of the clustering problem into a lesson, one used a placebo experiment to expose how untrustworthy a DiD that ignores serial correlation can be, and one carried the Delta method into industrial-scale online experiments. Each is laid out along the five elements of paper, method, data, results, and limitations.

### 5.1 Moulton (1990)

::: {.case}
Paper: "An Illustration of a Pitfall in Estimating the Effects of Aggregate Variables on Micro Units," The Review of Economics and Statistics 72(2):334-338. It turned the clustering problem of Section 3.4 into a lesson everyone should take, and the Moulton factor is named after it.

Method: Moulton focuses on an extremely common setup, placing a variable at an aggregate level (such as state or industry) into an individual-level regression as an explanatory variable, for example explaining individual wages with a state-level policy. He points out that individual errors are positively correlated within an aggregate group, and the aggregate explanatory variable is identical within the group, so the two together make ordinary OLS standard errors severely understated. He gives the inflation factor $\sqrt{1 + (\bar n - 1)\rho}$, quantifying the degree of understatement as a function of the within-group correlation $\rho$ and the group size $\bar n$.

Data: individual-level survey data paired with aggregate-level explanatory variables, exactly the most typical data structure in labor economics and policy evaluation.

Results: the conclusion is a repeatedly cited warning: when the explanatory variable varies at the aggregate level, ignoring within-group correlation makes the standard error absurdly small, and even if the within-group correlation looks tiny, as long as the groups are large the inflation factor will be large. His demonstration is especially striking, throwing into the regression some obviously irrelevant state-level variables, and even a computer-generated random number at the state level, and under naive OLS standard errors they still appear statistically significant. This short paper made the whole applied field realize how easily the significance of aggregate variables can be an illusion.

Limitations: the Moulton factor is an approximation under several simplifying assumptions (homogeneous within-group correlation, similar group sizes), and in real data with uneven group sizes and complex correlation structure it is only an order-of-magnitude reference, with precise inference still relying on the full cluster-robust variance and even the few-cluster corrections.
:::

The exemplary significance of this paper is that it did not stop at the slogan "you should cluster," but explained why, and by how much the understatement, through a single computable factor. The Moulton factor of Section 3.4 is its core, and the 4.18-fold predicted inflation on Halcyon is exactly this formula made flesh.

### 5.2 Bertrand, Duflo and Mullainathan (2004)

::: {.case}
Paper: "How Much Should We Trust Differences-in-Differences Estimates?," The Quarterly Journal of Economics 119(1):249-275. It used a carefully designed placebo experiment to nail to the table how untrustworthy DiD inference that ignores serial correlation is.

Method: the authors did something cruel and persuasive. They took real wage data and randomly assigned some states and some years a "fake policy" that in fact does not exist, then used standard DiD with ordinary standard errors to test the effect of this fake policy. Because the policy was randomly made up and the true effect is zero, an honest 5% test should reject only 5% of the time. They did this repeatedly and tallied the actual rejection rate.

Data: female wage data from the U.S. Current Population Survey (CPS), spanning many states and many years, the typical data for DiD studies.

Results: the actual rejection rate of ordinary DiD standard errors is about 45%, meaning that nearly half of the "significant policy effects" are fabricated out of pure noise. The root of the disease is that wages are strongly serially correlated within a state across years, while ordinary standard errors treat different years as independent information. They tested several remedies: clustering at the state level, the block bootstrap, and collapsing each state's time series into two points, before and after treatment. This paper almost single-handedly made "clustering at the state level" the standard operation for DiD studies.

Limitations: the clustering remedy itself needs enough clusters (states) to be reliable, and with few clusters it returns to the few-cluster problem of Section 4.2; their diagnosis is also mainly aimed at the single structure of serial correlation, and other forms of correlation need their own responses.
:::

This paper is the real-world version of the placebo health check in Section 1, and its inspiration. Halcyon's 64.5% false-positive rate under iid inference and Bertrand et al.'s 45% figure tell the same story: standard errors that ignore correlation make you see a sky full of signal in a field of noise. It elevated "the standard error is wrong" from a technical detail into a question about the credibility of an entire literature.

### 5.3 Deng, Knoblich and Lu (2018)

::: {.case}
Paper: "Applying the Delta Method in Metric Analytics: A Practical Guide with Novel Ideas," Proceedings of the 24th ACM SIGKDD International Conference on Knowledge Discovery and Data Mining (KDD '18), 233-242. It carried the Delta method of Section 3.2 into the daily work of industrial-scale online experiments, with authors from Microsoft's experimentation platform team.

Method: in large-platform A/B experiments, many core metrics are ratios (such as clicks per user, conversions per session), nonlinear functions of variables at the randomization-unit level, whose variance cannot be computed directly from the sample variance. The authors systematically use the Delta method to derive variances and confidence intervals for these ratio metrics and other nonlinear metrics, and handle two practical difficulties of online experiments: first, the within-randomization-unit (user) correlation of observations (a user's multiple visits are not independent, which is essentially the clustering problem of this chapter), and second, computational feasibility under massive data.

Data: large-scale A/B test metrics on Microsoft's online experimentation platform, with sample sizes routinely in the tens of millions.

Results: the paper gives a deployable Delta-method recipe that makes inference for ratio metrics both correct (folding within-user correlation into the variance) and efficient at industrial scale, becoming a practical reference for online experiment analysis. It demonstrates how the two main threads of this chapter (the Delta method for nonlinear functions and variance corrections for correlated data) converge in real platform data science.

Limitations: the Delta method relies on linearization, and when a ratio's denominator is near zero or a metric is highly skewed it still runs into the distortion of Section 3.2, where the bootstrap or a dedicated skewness correction is more stable; the paper's recipe is also mainly aimed at randomized experiments, and identification problems in observational data are out of its scope.
:::

The significance of this paper for the chapter is that it brings the Delta method from a textbook formula back to a high-frequency real setting: the platform has to compute confidence intervals for thousands of ratio metrics every day, and the Delta method, because it has a closed form and is fast, becomes the first choice, while its handling of within-user correlation echoes exactly the clustering thread of this chapter. It shows that the chapter's two tools are not isolated tricks but the foundation of the modern data-science inference stack.

Put the three together, and the point of anchoring becomes clear: Moulton explained the mechanics of clustering understatement, Bertrand, Duflo, and Mullainathan measured the cost of ignoring correlation with a placebo experiment, and Deng, Knoblich, and Lu industrialized the Delta method. They share the chapter's core, that a standard error must honestly reflect the data's correlation structure and the quantity's nonlinearity, or else significance is only an illusion.

## 6 A Full Walkthrough on Halcyon Data

Now let us run the earlier tools once through Halcyon from the top. The code below uses R 4.5.3 with the random seed fixed at 626, and every number cited in the text comes from the actual output of this code.

### 6.1 The Data-Generating Process

Three blocks of data, one seed. Delivery-pricing data: satisfaction $sat = \beta_0 + \beta_{time}\, time + \beta_{fee}\, fee + \varepsilon$, with $\beta_{time} = -0.08$ and $\beta_{fee} = -0.40$, so the true WTP is $\beta_{time}/\beta_{fee} = 0.20$ per minute; the main sample has ample variation in the delivery fee, while the weak-denominator variant holds the delivery fee nearly constant. City-experiment data: $G = 40$ cities with 200 sessions each, $y_{ic} = 2.0 + \tau D_c + a_c + \varepsilon_{ic}$, where the city shock $a_c$ manufactures within-cluster correlation and the true treatment effect is $\tau = 0.25$. Time-series data: $T = 400$ days with AR(1) errors.

```r
set.seed(626)
## Delivery pricing (WTP ratio)
time <- runif(3000, 10, 45); fee <- runif(3000, 1, 7)
sat  <- 5.0 - 0.08*time - 0.40*fee + rnorm(3000, 0, 1.5)
## City experiment (clustering)
a  <- rnorm(40, 0, 0.30)                          # city shock -> within-cluster correlation
D  <- rep(0:1, length.out = 40)[sample(40)]       # half the cities launch the policy
city <- rep(1:40, each = 200)
y  <- 2.0 + 0.25*D[city] + a[city] + rnorm(8000, 0, 1.0)
```

In the city experiment the within-cluster correlation is about $0.30^2/(0.30^2 + 1.0^2) = 0.083$. These few true values are the targets against which all later standard errors are reconciled: WTP 0.20, $\tau$ 0.25, and the width each sampling distribution should have.

### 6.2 The WTP Ratio: Delta versus Bootstrap

First estimate the delivery-pricing regression, where WTP is a ratio of two coefficients. The Delta method computes the variance with the gradient from Section 3.2, and the bootstrap resamples whole rows of data with replacement and recomputes the ratio each time.

```r
m <- lm(sat ~ time + fee)
wtp <- coef(m)["time"] / coef(m)["fee"]           # point estimate
grad <- c(0, 1/coef(m)["fee"], -coef(m)["time"]/coef(m)["fee"]^2)
se_delta <- sqrt(t(grad) %*% vcov(m) %*% grad)    # Delta standard error
```

In the main sample the delivery fee varies enough, $\beta_{fee}$ is estimated precisely, and the two routes agree strongly: the point estimate WTP is 0.207 (true value 0.20), the Delta standard error is 0.0105, and the bootstrap standard error is 0.0106, nearly equal. The three intervals also nearly coincide: Delta's $[0.186, 0.227]$, percentile's $[0.188, 0.229]$, and percentile-t's $[0.188, 0.228]$. When the denominator is strong, the Delta method is fast and accurate, and there is no need to bring in the bootstrap.

The weak-denominator variant compresses the variation in the delivery fee to almost nothing, $\beta_{fee}$ is estimated only at $-0.246$ (standard error 0.105, close to insignificant), the ratio's denominator hovers near zero, and the two routes immediately part ways. The Delta standard error is 0.138, giving a symmetric interval $[0.05, 0.59]$. But the bootstrap distribution is severely right-skewed (skewness as high as 8.0), because in a few resamples $\beta_{fee}$ grazes zero and the ratio flies off to the sky: the percentile interval is $[0.17, 1.51]$, a long right tail, and the percentile-t interval is $[0.24, 0.82]$. Delta's symmetric interval erases this right tail entirely, a misleading pretense of honesty.

![Left: when the delivery fee varies enough (strong denominator), the bootstrap sampling distribution of the WTP ratio (gray bars) and the Delta method's normal approximation (navy curve) almost coincide, the two routes agreeing. Right: when the delivery fee barely varies (weak denominator), the bootstrap distribution is severely right-skewed, dragging a long tail, while Delta's symmetric normal cannot cover that tail at all, and here only the bootstrap tells the truth. The red dashed line is the true value 0.20.](assets/fig/fig_06_delta_bootstrap.svg)

### 6.3 Clustering: How a t Value of 13.76 Shrinks

The city-experiment regression, comparing four standard errors.

```r
library(sandwich)
m <- lm(y ~ D)
se_iid <- sqrt(vcov(m)["D","D"])                          # iid default
se_hc1 <- sqrt(vcovHC(m, type = "HC1")["D","D"])          # heteroskedasticity-robust
se_cr1 <- sqrt(vcovCL(m, cluster = ~city)["D","D"])       # city cluster
```

The treatment-effect estimate is $\hat\tau = 0.318$ (true value 0.25). The four standard errors tell a complete story: the iid default is 0.0231 ($t = 13.76$), the heteroskedasticity-robust HC1 is also 0.0231 ($t = 13.76$, unmoved), the city cluster CR1 is 0.0908 ($t = 3.50$), and the small-sample-corrected CR2 is 0.0920. Clustering inflates the standard error 3.93-fold relative to iid, quite close to the Moulton factor's predicted 4.18-fold. The reason HC1 offers no improvement is exactly what Section 3.4 says, that it only fixes unequal variance and not within-cluster correlation, and the disease here is precisely the latter. From 13.76 to 3.50, an order-of-magnitude difference in significance, all in this one step of the standard error.

The real health check is to look at the false-positive rate. Repeatedly generate placebo data with $\tau = 0$ and test the correct null "the effect is zero" at the 5% level:

| Test method | $G = 40$ cities | $G = 10$ cities |
|---|---|---|
| naive iid | 0.645 | 0.636 |
| cluster (normal critical value) | 0.070 | 0.103 |
| cluster $t(G-1)$ | 0.064 | 0.067 |

The iid inference has a false-positive rate of about two thirds in both cases, more than ten times the nominal 5%, confirming the alarming conclusion of Section 1. City clustering pulls it back to near 5%, but with $G = 10$ the normal critical value is still too high (10.3%), and only switching to $t(G-1)$ presses it back to 6.7%, which is the lead-in to the few-cluster problem of the next subsection.

![Left: when repeatedly testing a null that is in fact true (the policy has no effect), the actual rejection rate of each method at the nominal 5% level. Naive iid inference has a false-positive rate of about two thirds at both G=40 and G=10 (gold bars), far above the 5% target marked by the red dashed line; city clustering (blue, navy) pulls it back near 5%. Right: the treatment-effect estimate from the same city experiment, with confidence intervals from three standard errors. The iid and heteroskedasticity-robust HC1 intervals are almost identically narrow, and only the city-cluster interval is honestly wide.](assets/fig/fig_06_cluster.svg)

### 6.4 Few Clusters: Trust the Normal or Trust the Bootstrap

When $G$ is small, even the cluster standard error is not enough, because its asymptotic validity is in the number of clusters $G \to \infty$. Take a $G = 10$ city experiment with a true effect. For the same cluster $t$ statistic ($t = 2.21$), three ways of setting the $p$ value give different answers:

```r
t10 <- coef(m10)["D"] / sqrt(vcovCL(m10, ~city)["D","D"])
p_z   <- 2 * pnorm(-abs(t10))                     # normal critical value
p_t   <- 2 * pt(-abs(t10), df = 9)                # t(G-1)
# wild cluster bootstrap: under the null, multiply each whole cluster by a random sign, recompute the distribution of t
```

The normal critical value gives $p = 0.027$, looking significant; $t(9)$ gives $p = 0.054$, already on the edge; the wild cluster bootstrap gives $p = 0.073$, not significant. The three $p$ values for the same $t$ differ entirely because of "how fat the tails should be when there are few clusters." The wild cluster bootstrap, because it simulates the distribution of $t$ directly from the cluster structure of the data without assuming normality, is the most trustworthy under few clusters, and its "actually not significant" is the honest conclusion.

### 6.5 HAC: The Same Lesson in Time Series

A daily-series regression, comparing iid and Newey-West standard errors.

```r
library(sandwich)
m <- lm(yt ~ xt)
se_iid <- sqrt(vcov(m)["xt","xt"])
se_nw  <- sqrt(NeweyWest(m, lag = 5, prewhite = FALSE, adjust = TRUE)["xt","xt"])
```

The slope estimate is 0.189. The iid standard error is 0.0484 ($t = 3.91$), and the Newey-West standard error is 0.0613 ($t = 3.08$, with lag $L = 5$), an inflation of 1.27-fold. Looking at coverage under repeated sampling is clearer: for a nominal 95% interval, iid actually covers only 84.9% and Newey-West lifts it to 92.0%. The direction is right, serial correlation does make iid undercover and HAC pulls it back some; but 92.0% is still short of 95%, honestly exposing the old complaint that HAC often does not fix things all the way in finite samples, which is also why people still consider the block bootstrap under strong serial correlation.

### 6.6 The Trinity of Tests and Wald's Reparameterization Trap

Use a Halcyon repurchase logit to test "user tenure does not affect repurchase" (some coefficient is zero), with the three tests side by side.

```r
mu <- glm(rep_order ~ tenure, family = binomial())  # unrestricted
mr <- glm(rep_order ~ 1,      family = binomial())  # restricted
W  <- coef(mu)["tenure"]^2 / vcov(mu)["tenure","tenure"]      # Wald
LR <- 2 * (logLik(mu) - logLik(mr))                          # likelihood ratio
```

The three statistics agree strongly: Wald 35.58, LR 36.93, LM 36.47, all far above the $\chi^2_1$ critical value, with $p$ values nearly zero, unanimously rejecting. Their tiny differences are exactly the finite-sample divergence Section 4.1 mentioned.

Wald's reparameterization trap can be demonstrated directly. Testing the same hypothesis "the coefficient equals 0.5," written in terms of the coefficient itself the Wald statistic is 10.68 ($p = 0.0011$), and written in terms of the $\log$ of the coefficient it is 6.79 ($p = 0.0092$). Two equivalent ways of writing the same statement, and the $p$ values differ by nearly an order of magnitude, simply because Wald relies on linearization and reparameterization changes the gradient. The LR test is immune, because it looks only at the likelihood value.

### 6.7 Grand Reconciliation

Put the key numbers together:

| Quantity | Method | Result | Note |
|---|---|---|---|
| WTP (strong denominator) | Delta / bootstrap SE | 0.0105 / 0.0106 | true value 0.20, two routes agree |
| WTP (weak denominator) | Delta / bootstrap SE | 0.138 / skewness 8.0 | Delta's symmetric interval erases the right tail |
| $\tau$ (city experiment) | iid / HC1 / cluster SE | 0.0231 / 0.0231 / 0.0908 | $t$ from 13.76 to 3.50 |
| false-positive rate $G=40$ | iid / cluster | 0.645 / 0.070 | nominal 5% |
| few clusters $G=10$ | normal / $t(9)$ / wild boot $p$ | 0.027 / 0.054 / 0.073 | same $t=2.21$ |
| HAC coverage | iid / Newey-West | 0.849 / 0.920 | nominal 95% |
| trinity | Wald / LR / LM | 35.58 / 36.93 / 36.47 | asymptotically equivalent |

The reconciliation of this section can be summarized as follows: when the denominator is strong, Delta and the bootstrap give almost identical standard errors, and when the denominator is weak only the bootstrap captures the right skew of the ratio's distribution; in the city experiment the iid standard error props the $t$ up to 13.76 and the heteroskedasticity-robust standard error is no help at all, and only clustering honestly lowers it to 3.50, while the iid inference reports a false positive two thirds of the time in the placebo test; with few clusters even clustering over-rejects, and $t(G-1)$ and the wild cluster bootstrap correct it back layer by layer; under serial correlation Newey-West raises the coverage from 84.9% toward 92.0%; and Wald, LR, and LM asymptotically give the same answer, but Wald changes with how the hypothesis is written.

## 7 Failure Modes and Robustness

The correlation structures in the simulations were manufactured, and in real data they hide deeper. This section works through the most common failure modes and their responses.

The choice of the level of the standard error is the first and most consequential judgment. The default standard error is almost always too small, because it assumes independence and homoskedasticity. Heteroskedasticity-robust standard errors correct for unequal variance but are powerless against correlation, which is the lesson of the unmoved HC1 in Section 1. Cluster-robust standard errors correct within-cluster correlation, but you have to pick the right level to cluster on, and the basis for picking is not a statistical test but a substantive judgment about the design: at which level do you believe "independent across clusters." When treatment is assigned at some level, or the model includes fixed effects at some level, you usually should cluster at that level or coarser. Better to cluster a bit coarsely (conservative, standard errors a bit large) than to cluster finely and miss correlation. This judgment has no automatic transmission, and the researcher must shoulder it head on.

Few clusters is the trap that follows immediately. Cluster-robust asymptotics is in the number of clusters, and with few clusters (empirically fewer than thirty or forty) it systematically over-rejects. Here you cannot just report a cluster standard error and be done, you should use $t(G-1)$ critical values, the CR2 correction, or the wild cluster bootstrap, and with very few or very unbalanced clusters the last is nearly required. A common mistake is thinking that millions of observations mean a large sample and stable inference, when in fact what decides the reliability of inference is the number of clusters, not the number of observations.

The Delta method for nonlinear functions distorts when the denominator is near zero, as warned repeatedly above. The signal to watch is whether the coefficient in the denominator is itself estimated precisely: if it is close to insignificant, the ratio's distribution is probably severely skewed, and you should abandon the symmetric Delta interval for the bootstrap's percentile-t or a dedicated Fieller interval. Likewise, any quantity on the boundary of the parameter space, or involving a maximum or an absolute value, should raise the alarm that the Delta method and the naive bootstrap fail together, and subsampling or a dedicated method is reliable.

The bootstrap itself can also deceive, especially on correlated data. Doing a naive bootstrap with replacement directly on time series or clustered data pretends the data are iid and, like the iid standard error, understates. Correlated data need the block bootstrap or the wild cluster bootstrap, preserving the correlation structure in the resampling. Moreover the number of bootstrap replications must be enough (a few thousand to start), because too few will let the interval endpoints themselves carry appreciable simulation noise.

At the testing level there are two things to watch. First, Wald is not invariant to reparameterization and relies on a variance that may be misspecified, so when in doubt use LR or a robust Wald based on a robust variance. Second, when you run many tests (a platform watching dozens of metrics at once), the 5% single-test false-positive rate accumulates into a high overall false-positive rate, and you need multiple-testing corrections (such as Bonferroni or controlling the false discovery rate), or you will always "discover" some significance that is purely accidental. More broadly, repeatedly testing until you get significance (p-hacking) systematically destroys the meaning of inference, and pre-registering the analysis plan and honestly reporting all tests is the prerequisite for a $p$ value to earn its literal meaning.

Stringing these failure modes together, the robustness of this chapter ultimately rests on one thing: think clearly about what correlations are in your data and how nonlinear your quantity is, and choose the variance algorithm and the test construction accordingly. The credibility of a standard error never lies in the number the software spits out by default, but in whether you have honestly told it the dependence structure of the data. Get this right, and significance is truly significant; get it wrong, and even the smallest $p$ value is only a byproduct of a miscomputed variance.

## 8 Further Reading

::: {.readings}
Required reading, in suggested order:

- Angrist and Pischke (2009, *Mostly Harmless Econometrics*), Chapter 8. The most readable introduction to clustering, serial correlation, and robust standard errors, and the source of this chapter's tone that "ignoring correlation is the same as overstating significance."
- Bertrand, Duflo and Mullainathan (2004, *QJE*). A placebo experiment exposing DiD that ignores serial correlation, the real-world version of the false-positive health check in Section 1, required.
- Moulton (1990, *RESTAT*). The original source of the clustering trap in aggregate-variable regressions and of the Moulton factor, short and far-reaching.
- Cameron and Miller (2015, *Journal of Human Resources*). A comprehensive practical guide to cluster-robust inference, from the level of clustering to few-cluster corrections, a self-check list to run before doing empirical work.

Further reading:

- Efron (1979, *Annals of Statistics*). The founding work of the bootstrap, for understanding the source of the resampling idea.
- White (1980, *Econometrica*) and Newey and West (1987, *Econometrica*). The original sources of the heteroskedasticity-robust and HAC variances, the two cornerstones of this chapter's sandwich filling.
- Cameron, Gelbach and Miller (2008, *RESTAT*) and MacKinnon and Webb (2017, *Journal of Applied Econometrics*). The wild cluster bootstrap and its behavior under unbalanced clusters, the modern workhorse of few-cluster inference.
- Imbens and Kolesár (2016, *RESTAT*). Practical advice on robust standard errors in small samples, the source of the CR2 and $t(G-1)$ corrections, the technical backing for the few-cluster section.
- Abadie, Athey, Imbens and Wooldridge (2023, *QJE*). Re-answering "should you cluster at all, and at which level" from the design perspective, upgrading clustering from a convention into a defensible judgment.
- Horowitz (2001, *Handbook of Econometrics* Vol. 5). The authoritative survey of the bootstrap, with a complete technical treatment of higher-order refinement and failure cases.
- Deng, Knoblich and Lu (2018, *KDD*). The deployment of the Delta method in industrial-scale online experiments, the interface between this chapter's methods and platform data science.
:::

::: {.apa-refs}
- Abadie, A., Athey, S., Imbens, G. W., & Wooldridge, J. M. (2023). When should you adjust standard errors for clustering? *The Quarterly Journal of Economics, 138*(1), 1-35. https://doi.org/10.1093/qje/qjac038
- Angrist, J. D., & Pischke, J.-S. (2009). *Mostly harmless econometrics: An empiricist's companion*. Princeton University Press. https://doi.org/10.1515/9781400829828
- Bertrand, M., Duflo, E., & Mullainathan, S. (2004). How much should we trust differences-in-differences estimates? *The Quarterly Journal of Economics, 119*(1), 249-275. https://doi.org/10.1162/003355304772839588
- Cameron, A. C., Gelbach, J. B., & Miller, D. L. (2008). Bootstrap-based improvements for inference with clustered errors. *The Review of Economics and Statistics, 90*(3), 414-427. https://doi.org/10.1162/rest.90.3.414
- Cameron, A. C., & Miller, D. L. (2015). A practitioner's guide to cluster-robust inference. *Journal of Human Resources, 50*(2), 317-372. https://doi.org/10.3368/jhr.50.2.317
- Deng, A., Knoblich, U., & Lu, J. (2018). Applying the Delta method in metric analytics: A practical guide with novel ideas. In *Proceedings of the 24th ACM SIGKDD International Conference on Knowledge Discovery & Data Mining* (pp. 233-242). Association for Computing Machinery. https://doi.org/10.1145/3219819.3219919
- Efron, B. (1979). Bootstrap methods: Another look at the jackknife. *The Annals of Statistics, 7*(1), 1-26. https://doi.org/10.1214/aos/1176344552
- Horowitz, J. L. (2001). The bootstrap. In J. J. Heckman & E. Leamer (Eds.), *Handbook of econometrics* (Vol. 5, pp. 3159-3228). Elsevier. https://doi.org/10.1016/S1573-4412(01)05005-X
- Imbens, G. W., & Kolesár, M. (2016). Robust standard errors in small samples: Some practical advice. *The Review of Economics and Statistics, 98*(4), 701-712. https://doi.org/10.1162/REST_a_00552
- MacKinnon, J. G., & Webb, M. D. (2017). Wild bootstrap inference for wildly different cluster sizes. *Journal of Applied Econometrics, 32*(2), 233-254. https://doi.org/10.1002/jae.2508
- Moulton, B. R. (1990). An illustration of a pitfall in estimating the effects of aggregate variables on micro units. *The Review of Economics and Statistics, 72*(2), 334-338. https://doi.org/10.2307/2109724
- Newey, W. K., & West, K. D. (1987). A simple, positive semi-definite, heteroskedasticity and autocorrelation consistent covariance matrix. *Econometrica, 55*(3), 703-708. https://doi.org/10.2307/1913610
- White, H. (1980). A heteroskedasticity-consistent covariance matrix estimator and a direct test for heteroskedasticity. *Econometrica, 48*(4), 817-838. https://doi.org/10.2307/1912934
:::
