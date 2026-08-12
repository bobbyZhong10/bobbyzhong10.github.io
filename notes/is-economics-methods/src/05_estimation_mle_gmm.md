---
title: "Estimation: MLE, M-Estimation & GMM"
subtitle: "The Extremum Estimator and Its Sandwich"
seriesline: "Foundations of Information Systems Economics · Chapter 5"
series: "Foundations of Information Systems Economics: Empirical Methods for Information Systems Economics"
sidebartitle: "Foundations of Information Systems Economics"
sidebarsub: "Theory · Methods · Applications"
chapterlabel: "Chapter 5 · Estimation: MLE, M-Estimation & GMM"
---

## Introduction

Helix wants to predict whether a developer will adopt a new tool. An analyst runs a linear regression, and for some firms the predicted adoption probability comes out to 121%. The software throws no error, the coefficients even carry standard errors, only the "probability" has wandered across a boundary it was supposed to respect. This small accident reminds us of something: estimation is not finished the moment you hand the data to an optimizer. The objective function, the parameter space, and the economic model together decide what the number you compute actually is.

Maximum likelihood, nonlinear least squares, and the generalized method of moments look like three separate methods, but they are all searching for an extremum of a sample objective function. Once you place them inside the unified framework of the extremum estimator, three questions keep recurring: whether the population objective function has a unique extremum at the truth, whether the sample objective function will settle near the population objective function, and how the curvature and the noise around the extremum determine sampling error. Identification, consistency, and asymptotic normality then become one continuous argument rather than three disconnected lists of theorems.

The Helix case also displays a subtler kind of error: the distribution is misspecified, the point estimate can still aim at the right conditional mean, yet the standard errors the model reports leave a nominal 95% confidence interval covering the truth only 71.5% of the time. Why the sandwich variance can fix this, when the information matrix equality holds, and how GMM's moment conditions and weighting matrix change efficiency, form the spine of this chapter. An estimate usually occupies a single cell of a results table. What you really need to understand is why that cell should land near its target, and on what grounds you trust the parentheses beside it.

## 1 A Probability Printed as 121%

::: {.case}
Helix is a fictional B2B developer platform that sells a range of cloud tools to a large base of client firms. We use simulated data to tell this story, and the payoff is that the data-generating process is fully known, so every estimator can be reconciled against the truth, a teaching condition that real data can never give you. Each client firm $i$ has an observable, standardized engagement measure $x_i$ (past API activity, z-standardized), with mean zero and standard deviation one. The platform cares about two related things. The first is adoption: whether a firm turns on a paid Pro-tier tool, recorded as a binary variable $D_i \in \{0,1\}$. The second is usage: after turning it on, the number of monthly API calls to that tool (in thousands), recorded as a nonnegative integer $y_i$. This section looks at adoption first.

The platform's data scientist wants to know something straightforward: are more engaged firms more likely to adopt Pro, and by how much. She has a cross section of $n = 4000$ firms, with both $D_i$ and $x_i$ in hand. The most natural first step is to treat adoption as an ordinary outcome variable and run a linear regression on $x$, that is, a linear probability model (LPM): $D_i = \gamma_0 + \gamma_1 x_i + u_i$, estimated by OLS.

It comes back as $\hat\gamma_1 = 0.215$ (SE 0.007), with intercept $\hat\gamma_0 = 0.427$. Read literally: a one-standard-deviation increase in engagement raises the adoption probability by about 21.5 percentage points. The number is not outrageous, but when she prints out the model's fitted values as "predicted adoption probabilities" firm by firm, something breaks.
:::

Among the fitted values, the highest firm has a predicted adoption probability of 1.213 and the lowest has $-0.478$. Across the full sample, 2.93% of firms are predicted to have a "probability" outside $[0,1]$, of which 2.5% fall below zero and 0.4% rise above one. A probability printed as 121% or as $-48\%$ is not a precision problem, it is a model that simply does not know a probability belongs in $[0,1]$. This is the chapter's first bad number, and it is bad in a blunter way than the usual identification bias: nothing was estimated crookedly, the objective function was mismatched to the problem from the very start.

Where the trouble comes from can be said in a sentence. OLS minimizes the residual sum of squares, which by default lets the outcome variable range freely over the whole real line, and it uses a straight line to approximate the conditional mean $\mathbb{E}[D \mid x]$. But when the outcome is binary, that conditional mean is exactly the adoption probability $\Pr(D = 1 \mid x)$, which must stay in $[0,1]$, and a straight line with no ceiling and no floor cannot manage this. The subtler second flaw is that the LPM locks engagement's marginal effect to a constant 0.215, the same for every firm. Intuitively that is wrong: a firm already dead set on adopting will not be moved by a bit more engagement, and neither will a firm that will never adopt; the firms actually pushed around by engagement are the ones already teetering on the edge of adoption. The marginal effect ought to vary across firms, and the LPM assumes it does not.

::: {.intuition}
Picture the LPM's predicament as using a ruler with no upper mark to measure something that can only sit between 0 and 1. The ruler itself is fine, the mistake is measuring a bounded quantity with it. The fix is not to file the ruler down, but to swap in a different apparatus: let the linear index $x_i'\beta$ vary freely over the real line first, then use a monotone function $F(\cdot)$ that squeezes the whole real line into $[0,1]$ to translate it into a probability. When $x_i'\beta$ is large the probability hugs 1, when it is small the probability hugs 0, and it never crosses the boundary. Take $F$ to be the standard logistic and you have logit; take the standard normal CDF and you have probit. Once the apparatus changes, "minimize the residual sum of squares" is no longer the right objective, and in its place comes one that fits the probability structure better: the likelihood. What this chapter turns to next is exactly how to write down the correct objective function for a model, how to maximize it, and how to attach an honest standard error to the result.
:::

The lesson of Section 1 can be summarized as follows: the first step of estimation is to choose the right objective function for your model, and if you choose wrong, no matter how elegant the algorithm, it is only optimizing the wrong thing. A binary outcome does not want a residual sum of squares, it wants a likelihood, and this pushes us straight toward MLE and the framework behind it that unifies all extremum estimators.

## 2 The Economic Model and the Estimand

Before touching estimation, get two things clear: what model sits behind Helix's adoption behavior, and which quantity in that model we actually want to estimate.

### 2.1 A Random Utility Model of Adoption

Whether a firm adopts can be viewed through a random utility model (RUM). The net benefit firm $i$ obtains from adopting Pro is written as an observable part plus an unobservable part:

$$U_i^* = \beta_0 + \beta_1 x_i - \varepsilon_i.$$

The observable part $\beta_0 + \beta_1 x_i$ varies linearly with engagement, and $\beta_1 > 0$ means a more active firm gets a higher net benefit from adopting. The term $\varepsilon_i$ is everything the platform cannot observe: the firm's internal procurement process, the decision maker's preferences, how tight the budget is this period. The firm adopts when the net benefit is positive, $D_i = \mathbf{1}\{U_i^* > 0\} = \mathbf{1}\{\varepsilon_i < \beta_0 + \beta_1 x_i\}$. So the adoption probability is the probability that the unobservable falls below the threshold:

$$\Pr(D_i = 1 \mid x_i) = \Pr(\varepsilon_i < \beta_0 + \beta_1 x_i) = F(\beta_0 + \beta_1 x_i),$$

where $F$ is the cumulative distribution function of $\varepsilon_i$. This step turns the patch of Section 1, "use a monotone function to squeeze the real line into $[0,1]$", into a model with a behavioral interpretation: $F$ is not some squashing function picked off the shelf, it is the distribution of the unobservable benefit shock. Take $\varepsilon_i$ to be standard logistic, $F(z) = \Lambda(z) = 1 / (1 + e^{-z})$, and you get logit; take standard normal and you get probit. The two shapes are very close, logit has slightly thicker tails and is a bit cleaner to compute, and this chapter works mainly with it.

::: {.warning}
Here lies a constraint you cannot get around in discrete choice models, one that Chapter 19 will develop in full when it takes up discrete choice, but for now we state only the minimum necessary point: the absolute size of $\beta$ and the scale of $\varepsilon$'s distribution are tangled together, and only their ratio is identified. Double the utility, the threshold, and the shock all at once, and the adoption probability does not move a hair, so the data cannot tell them apart. The convention is to fix the scale of $\varepsilon$ (the logistic variance set to $\pi^2/3$, the normal to 1), so the estimated $\beta$ is a coefficient "relative to this scale." This is also why logit coefficients and probit coefficients cannot be compared directly in magnitude, each has been divided by a different scale. This chapter treats the scale as already normalized and focuses on the estimation machinery itself.
:::

### 2.2 What Exactly Are We Estimating

With the model set up, the estimand has two layers. The first layer is the model parameter itself, $\theta_0 = (\beta_0, \beta_1)$. It is the true parameter of the data-generating process, and later in this chapter we will give it a more abstract and more useful definition: the point at which a population objective function attains its extremum. The second layer is a quantity read off from the parameters that carries direct managerial meaning, called the marginal effect.

Why do we need the second layer? Because in logit $\beta_1$ itself is not "how many points the adoption probability rises for each one-unit increase in engagement." The derivative of the probability with respect to $x$ is

$$\frac{\partial \Pr(D_i = 1 \mid x_i)}{\partial x_i} = f(\beta_0 + \beta_1 x_i)\, \beta_1,$$

where $f = F'$ is the density. This marginal effect varies with $x$, which is exactly the point made in Section 1: engagement's influence is largest near an adoption probability of 0.5 (where the density $f$ is highest) and tends to zero at both ends. To collapse this into a single number to report, there are two common conventions. One is the average marginal effect (AME), averaging each firm's marginal effect, $\mathbb{E}[f(\beta_0 + \beta_1 x)\,\beta_1]$; the other is the marginal effect at the mean (MEM), computed once at $x = 0$, $f(\beta_0)\,\beta_1$. In Helix's truth the AME is 0.216 and the MEM is 0.264, and the two differ because the marginal effect is nonlinear, so averaging the regressor first and then computing, versus computing first and then averaging, naturally give different answers. What this section wants you to remember is: $\beta_1$ is the model's internal parameter, and the AME is the "how many percentage points" you can write into a report, with the bridge between them being the derivative formula above.

As for usage $y_i$, Section 2 sets it aside for now. Its model (a count model with a log mean) waits until Section 3, where we discuss QMLE, because it happens to be the best example of "half misspecified."

The main points of Section 2 can be summarized as follows: Helix's adoption is a random utility model, the adoption probability equals the probability that the unobservable benefit shock falls below the threshold, and taking the shock distribution to be logistic gives logit; the first layer we want to estimate is the model parameter $\theta_0$, and the second layer is the AME read off from the parameters, the latter being the managerially meaningful "marginal percentage points." The next section takes "estimating $\theta_0$" itself and places it inside the unified framework of the extremum estimator, making clear why it is consistent and what its variance looks like.

## 3 Identification and Large-Sample Theory

By now we have a model and a target parameter $\theta_0$, but we have not yet said clearly how to estimate it from the data, still less why what we estimate can be trusted. This section answers both questions, and answers them at a unified altitude, letting MLE and GMM share one argument for consistency and asymptotic normality. Logically it precedes any specific algorithm: whether an estimator is consistent or not, and what its variance looks like, depends on the shape of the objective function and on whether the specification is correct, and has nothing to do with whether you use Newton or a quasi-Newton method to find the extremum. This is the most finely detailed section of the chapter. We first set up the extremum estimator as a unified framework (3.1), then state the identification condition (3.2), then use Jensen's inequality to prove consistency (3.3), derive asymptotic normality and read off the sandwich variance (3.4), see how MLE uses the information matrix equality to collapse the sandwich into $I^{-1}$ (3.5), see how QMLE survives once that equality breaks (3.6), and finally mount estimation directly on moment conditions to get GMM (3.7).

### 3.1 The Extremum Estimator: A Unified Framework

Start with the intuition. Imagine you are searching for the highest point on a rolling terrain, where the height at each point is determined by the parameter $\theta$, and all you have is a noisy contour map drawn from the sample. Your estimate is the highest point on this sample terrain map. The true parameter $\theta_0$ is the highest point on the population terrain map you cannot see. The larger the sample, the closer the sample terrain is to the population terrain, and the closer the two summits are to each other. That is the whole intuition of consistency.

Write this sentence as a definition. An extremum estimator is the extremum point of some sample objective function $Q_n(\theta)$:

$$\hat\theta = \arg\max_{\theta \in \Theta} Q_n(\theta).$$

Different $Q_n$ give different estimators, and the ones most used in econometrics are all special cases of this one line (below, $z_i$ denotes one observation, which in Helix is $(D_i, x_i)$ or $(y_i, x_i)$):

::: {.theorem}
**(Several estimators are the same sentence.)** Different sample objective functions $Q_n(\theta)$ give different estimators. MLE maximizes the average log-likelihood:

$$\hat\theta_{MLE} = \arg\max_\theta \frac{1}{n} \sum_{i=1}^n \ln f(z_i; \theta).$$

Nonlinear least squares (NLS) minimizes the residual sum of squares, where $m(x_i, \beta)$ is the specified form of the conditional mean:

$$\hat\beta_{NLS} = \arg\min_\beta \frac{1}{n} \sum_{i=1}^n \big(y_i - m(x_i, \beta)\big)^2.$$

GMM minimizes a weighted quadratic form of the moment conditions:

$$\hat\theta_{GMM} = \arg\min_\theta\ \bar g_n(\theta)'\, W\, \bar g_n(\theta), \qquad \bar g_n(\theta) = \frac{1}{n}\sum_{i=1}^n g(w_i, \theta).$$

Maximizing with a sign flip is minimizing, so all three belong to the same family of extremum estimators.
:::

OLS is the special case of NLS with $m(x,\beta) = x'\beta$, and also the special case of GMM with the moment condition $g(w_i,\beta) = x_i(y_i - x_i'\beta)$ and weight $W = I$. This "everything reduces to one argmax" viewpoint is not for show; its practical payoff is that once you prove, a single time, under what conditions an extremum estimator is consistent and what its variance looks like, MLE, GMM, and NLS are all proved at once, with no need to go through each separately. The rest of this section walks that unified path.

### 3.2 Identification: The Objective Function Must Have a Unique, Separable Peak

The sample objective function $Q_n(\theta)$ converges pointwise, as the sample grows, to a population objective function $Q_0(\theta) = \mathbb{E}[q(z, \theta)]$ (for MLE this is $\mathbb{E}[\ln f(z;\theta)]$). For the estimator to pin down the truth, the first prerequisite is that this population objective function attains its extremum at $\theta_0$, and only at $\theta_0$.

::: {.assumption}
**(Identification.)** The population objective function $Q_0(\theta)$ attains a unique maximum at $\theta_0$ over the parameter space $\Theta$: for any $\theta \neq \theta_0$, $Q_0(\theta) < Q_0(\theta_0)$.
:::

The word "unique" carries the entire weight of identification. If there were another $\theta_1 \neq \theta_0$ with $Q_0(\theta_1) = Q_0(\theta_0)$, then at the population level the data cannot distinguish the two parameters, and no amount of sample would help, which is a failure of identification, not a shortfall of precision. The scale non-identification of Section 2.1 is one example: scale up $\beta$ and $\varepsilon$ proportionally and the likelihood is identical, the objective function develops a ridge of equal height, and the peak is not unique.

But identification actually comes in degrees, strong and weak, rather than the crisp binary of present or absent, and this point matters enormously in practice. Even when the peak is unique, it may not be "separable enough": the objective function is nearly flat around $\theta_0$, with a broad swath of $\theta$ values all giving nearly the same objective value. In that case $\theta_0$ is barely unique in the population, yet in a sample it is extremely hard to tell apart from nearby values, and the estimator's variance grows alarmingly large.

::: {.intuition}
Picture the strength of identification as the steepness of the peak. A lone peak rising sharply from the ground, you could grope your way to the top with your eyes closed, is strong identification: the objective function has high curvature near the truth, so a slight departure from the truth makes the objective value drop noticeably, and the data has strong opinions about which $\theta$ is better. A plateau swelling gently upward, where the top is anyone's guess, is weak identification: the objective function is nearly flat, countless $\theta$ values look equally good, and the data has almost no opinion. The steepness of the peak, as we will see, is exactly the curvature $A$ behind the bread $A^{-1}$ in the asymptotic variance; the steeper the peak, the larger $A$, the smaller $A^{-1}$, and the smaller the variance. This is why the geometry of identification (how steep the peak is) and the precision of inference (how large the variance is) are two sides of the same thing. Section 4.3 will place a lone peak and a plateau side by side in Helix and measure their difference.
:::

### 3.3 Consistency: Why Jensen's Inequality Pins the Summit to the Truth

The identification assumption requires the population objective function to attain its extremum at the truth. But for MLE, on what grounds does this hold? Why does "maximizing the expected log-likelihood" lead us exactly to the true parameter and not somewhere else? The answer is an argument that is surprisingly short and surprisingly beautiful, using only Jensen's inequality.

Start with the intuition. You hold a family of candidate probability densities $f(\cdot; \theta)$, one of which, $\theta_0$, is true. The data are generated by $f(\cdot; \theta_0)$. If you use a wrong $\theta$ to explain this data, on average you will feel the data "does not look much like" it came from $f(\cdot; \theta)$, that is, the expected log-likelihood is low; only with the right $\theta_0$ does the data most "look like" itself, and the expected log-likelihood is highest. This measure of "a wrong density always explains worse" is precisely the Kullback-Leibler divergence from information theory, which is nonnegative and equals zero only when the two distributions are equal. The algebra below just makes this sentence precise.

::: {.theorem}
**(The core of MLE consistency.)** If the model is correctly specified and $\theta_0$ is identified, then the expected log-likelihood is maximized at the truth: for any $\theta$, $\mathbb{E}[\ln f(z; \theta_0)] \geq \mathbb{E}[\ln f(z; \theta)]$.
:::

The proof is only a few lines, and every line is worth staring at. Look at the likelihood ratio $\frac{f(z;\theta)}{f(z;\theta_0)}$, take its log and negate it, and use the fact that $g(a) = -\ln a$ is convex together with Jensen's inequality $\mathbb{E}[g(Y)] \geq g(\mathbb{E}[Y])$:

$$\mathbb{E}_z\left[ -\ln \frac{f(z;\theta)}{f(z;\theta_0)} \right] \geq -\ln \mathbb{E}_z\left[ \frac{f(z;\theta)}{f(z;\theta_0)} \right].$$

The expectation inside the brackets on the right is taken over the true distribution $f(\cdot;\theta_0)$, the density cancels, and it equals exactly 1:

$$\mathbb{E}_z\left[ \frac{f(z;\theta)}{f(z;\theta_0)} \right] = \int \frac{f(z;\theta)}{f(z;\theta_0)}\, f(z;\theta_0)\, dz = \int f(z;\theta)\, dz = 1.$$

And $-\ln 1 = 0$, so the left side is $\geq 0$, and splitting apart the log on the left gives

$$\mathbb{E}[\ln f(z; \theta_0)] - \mathbb{E}[\ln f(z; \theta)] \geq 0.$$

This is exactly what we wanted to prove: the expected log-likelihood is maximized at $\theta_0$. Done. Equality holds if and only if $f(\cdot;\theta) = f(\cdot;\theta_0)$, and the identification assumption guarantees this happens only at $\theta = \theta_0$, so the peak is unique.

::: {.intuition}
The weight of this proof lies in how it upgrades "maximize the likelihood" from an operational recipe into a guaranteed identification strategy. The sample log-likelihood $\frac{1}{n}\sum \ln f(z_i;\theta)$ converges pointwise, by the law of large numbers, to $\mathbb{E}[\ln f(z;\theta)]$, whose summit is pinned by Jensen at $\theta_0$. So the summit of the sample objective, $\hat\theta$, converges to the summit of the population objective, $\theta_0$, and this is consistency. Convexity (the convexity of $-\ln$) is the protagonist here, guaranteeing the direction of that decisive inequality. The entire argument uses nothing about the specific form of $f$, so it treats logit, probit, Poisson, and normal alike, as long as the likelihood is written correctly.
:::

To upgrade this pointwise convergence rigorously to $\hat\theta \xrightarrow{p} \theta_0$, the technical apparatus also needs uniform convergence of the objective function in $\theta$, plus compactness of the parameter space and other regularity conditions, so that the convergence of the summit can be pushed from "every point converges" to "the extremum point converges." These conditions are standard and this chapter does not spell them out; Newey and McFadden (1994) give a complete treatment. The point is that the inequality fixing the direction comes from Jensen, not from any fine technical assumption.

### 3.4 Asymptotic Normality and the Sandwich

Consistency says $\hat\theta$ will converge to $\theta_0$, but it does not say at what speed or in what distribution it fluctuates around $\theta_0$. This subsection gives the answer, and gives it in the form that applies generally to extremum estimators, with MLE and GMM as special cases.

Start with the intuition, which will guide the algebra below all the way through. The reason $\hat\theta$ lands, in each draw, near $\theta_0$ but not exactly on it, is the outcome of a tug-of-war between two forces. One force is the jitter of the sample gradient at the truth: the first-order condition should be zero at $\theta_0$, but in a finite sample it carries noise, and the larger this noise, the farther $\hat\theta$ is pushed away from $\theta_0$. The other force is the curvature of the objective function: the steeper the peak, the shorter the distance the same gradient noise can push the extremum point, because a steep peak clamps $\hat\theta$ tightly. So the shape of the variance is already suggested: the numerator is the gradient's noise, the denominator is the curvature, with the noise in the middle and the curvature clamping it on both sides. This is where the sandwich $A^{-1} B A^{-1}$ comes from, and the Taylor expansion below just writes it precisely.

The engine of the derivation is a Taylor expansion of the first-order condition. $\hat\theta$ is the extremum point of the sample objective, so it drives the sample gradient to zero: $\nabla Q_n(\hat\theta) = 0$. Expand this gradient to first order around the truth $\theta_0$,

$$0 = \nabla Q_n(\hat\theta) \approx \nabla Q_n(\theta_0) + \nabla^2 Q_n(\theta_0)\,(\hat\theta - \theta_0),$$

rearrange, multiply by $\sqrt{n}$, and get

$$\sqrt{n}(\hat\theta - \theta_0) \approx -\big[\nabla^2 Q_n(\theta_0)\big]^{-1} \sqrt{n}\, \nabla Q_n(\theta_0).$$

Each of the two pieces on the right has a destination. The sample Hessian $\nabla^2 Q_n(\theta_0)$ converges by the law of large numbers to its expectation, which we write as $-A = \mathbb{E}[\nabla^2 q(z,\theta_0)]$, so $A$ is the negative curvature of the population objective at the truth, positive definite, measuring how steep the peak is. The sample gradient $\sqrt{n}\,\nabla Q_n(\theta_0) = \frac{1}{\sqrt n}\sum_i \nabla q(z_i, \theta_0)$ is a sum of mean-zero random terms (the gradient's expectation is zero at $\theta_0$) and converges by the central limit theorem to $\mathcal{N}(0, B)$, where $B = \mathbb{E}[\nabla q(z,\theta_0)\,\nabla q(z,\theta_0)']$ is the variance of the gradient. Put the two pieces together:

::: {.theorem}
**(Asymptotic normality of the extremum estimator.)** Under regularity conditions,

$$\sqrt{n}(\hat\theta - \theta_0) \xrightarrow{d} \mathcal{N}\big(0,\ A^{-1} B A^{-1}\big),$$

where $A = -\,\mathbb{E}[\nabla^2 q(z, \theta_0)]$ is the negative expected Hessian of the objective function at the truth (the curvature), and $B = \mathbb{E}[\nabla q(z, \theta_0)\, \nabla q(z, \theta_0)']$ is the variance of the gradient.
:::

This $A^{-1} B A^{-1}$ is the protagonist of the chapter, the sandwich variance. The two slices of bread $A^{-1}$ are the reciprocal of the curvature, and the filling $B$ is the gradient's noise. Its shape says two things clearly: the steeper the peak (the larger $A$), the smaller the variance, because a steep peak clamps $\hat\theta$ tighter; the larger the gradient noise (the larger $B$), the larger the variance, because the more the first-order condition jitters from draw to draw, the farther the extremum point is pushed. Section 3.2 said the strength of identification is the steepness of the peak, and that corresponds exactly to $A$ here; when identification is weak and the peak is flat, $A$ tends to zero, $A^{-1}$ explodes, and the variance explodes with it, which is what weak identification looks like inside the variance formula.

These two slices of bread and this filling are not abstract symbols. Section 6 will measure them out on Helix: for the same usage parameter, using only one slice of bread, treating the filling as equal to the bread (that is, misusing $I^{-1}$), gives a standard error of 0.0090, while honestly computing the full sandwich gives 0.0182, a full factor of two apart, and the truth stands on the side of the latter.

The sandwich is the default shape of the variance of any extremum estimator, whether or not the objective function is a likelihood. The next subsection will show that MLE, when correctly specified, enjoys a special grace that lets this sandwich collapse; but that grace is conditional, and once the condition fails, the sandwich is the only reliable shape.

### 3.5 MLE's Special Grace: The Information Matrix Equality and Efficiency

When the objective function happens to be a correctly specified log-likelihood, there is a nontrivial equality between $A$ and $B$ that lets the sandwich collapse into a single slice of bread. This equality is called the information matrix equality.

State it first. For MLE, the gradient is the score $s(z,\theta) = \nabla \ln f(z;\theta)$, and both the expectation of the negative Hessian and the variance of the score equal the Fisher information $I(\theta_0)$:

$$I(\theta_0) = -\,\mathbb{E}\!\left[\frac{\partial^2 \ln f}{\partial\theta\,\partial\theta'}(z;\theta_0)\right] = \mathbb{E}\!\left[\frac{\partial \ln f}{\partial\theta}(z;\theta_0)\,\frac{\partial \ln f}{\partial\theta'}(z;\theta_0)\right].$$

The left side is the curvature $A$, the right side is the gradient variance $B$, and the equality says $A = B = I(\theta_0)$. Substitute into the sandwich, $A^{-1} B A^{-1} = I^{-1} I\, I^{-1} = I^{-1}$, and the sandwich collapses:

$$\sqrt{n}(\hat\theta_{MLE} - \theta_0) \xrightarrow{d} \mathcal{N}\big(0,\ I(\theta_0)^{-1}\big).$$

At first glance this equality looks strange: curvature and gradient variance are two seemingly unrelated quantities, so why should they happen to be equal? The intuition is that they are really two windows onto the same thing, both measuring how firmly the data pins down $\theta$. One window is the curvature: the more sharply the log-likelihood bends at its summit, the faster the likelihood drops as you depart slightly from the truth, and the stronger the data's opinion about $\theta$. The other window is the variance of the score: a sharply bending likelihood is equivalent to a score (slope) that changes fast along $\theta$, so it also jitters more across different samples. A sharply bending likelihood must come with a fiercely jittering score, and the two windows see the same number, the Fisher information. The two differentiations below just squeeze this "same thing" out algebraically.

Why does this equality hold? The root is an identity: the density integrates to 1, $\int f(z;\theta)\, dz = 1$, and differentiating with respect to $\theta$, the right side is a constant with derivative zero:

$$0 = \frac{\partial}{\partial\theta} \int f(z;\theta)\, dz = \int \frac{\partial f}{\partial\theta}(z;\theta)\, dz = \int \frac{\partial \ln f}{\partial\theta}(z;\theta)\, f(z;\theta)\, dz = \mathbb{E}\!\left[\frac{\partial \ln f}{\partial\theta}(z;\theta_0)\right],$$

where the third equality uses $\frac{\partial \ln f}{\partial\theta} = \frac{1}{f}\frac{\partial f}{\partial\theta}$. This shows the score has expectation zero at the truth, which is exactly the incarnation, in the likelihood case, of Section 3.4's "the gradient has expectation zero." Differentiate this identity once more (differentiate $\int \frac{\partial \ln f}{\partial\theta} f\, dz = 0$ with respect to $\theta$, using the product rule) and you get $\mathbb{E}[\frac{\partial^2 \ln f}{\partial\theta\partial\theta'}] + \mathbb{E}[\frac{\partial \ln f}{\partial\theta}\frac{\partial \ln f}{\partial\theta'}] = 0$, and rearranging gives the information matrix equality. The point is: this equality depends entirely on $f$ being the true distribution of the data, and only two differentiations can wring out the information in "integrates to 1"; once $f$ is not the true distribution, the first step $\int f = 1$ still holds for the wrong $f$, but it is no longer the distribution of the data, so the step where the score has expectation zero need not match the data's true expectation, and the equality breaks with it. This foreshadows Section 3.6.

The information matrix equality also incidentally crowns MLE with efficiency. Here $I(\theta_0)$ is the Fisher information contributed by a single observation, and $n$ independent observations together supply information $n\, I(\theta_0)$. The Cramér-Rao lower bound says that any unbiased estimator built from $n$ observations has variance no smaller than the inverse of the total information $[n\, I(\theta_0)]^{-1}$, equivalently, the asymptotic variance of the standardized estimator is no smaller than the inverse of the single-observation information:

$$\mathrm{Avar}\big(\sqrt{n}(\hat\theta - \theta_0)\big) \geq I(\theta_0)^{-1}.$$

And MLE's asymptotic variance attains exactly this bound (matching the $\mathcal{N}(0, I(\theta_0)^{-1})$ above). In other words, when the likelihood is correctly specified, no (regular) estimator can be more precise than MLE, and this is the answer to the question "if MLE is most efficient, why use anything else": as long as you dare to guarantee the distribution is fully correct, you indeed need nothing else. The collapse of the sandwich shows up directly on Helix's adoption data in Section 6: the logit likelihood is correctly specified there, and the slope's Fisher-information standard error and sandwich standard error are both 0.044, identical to the digit, so the two slices of bread and the filling really are equal. The trouble is in that "dare to guarantee," and the next section pulls it apart.

### 3.6 When the Information Matrix Equality Breaks: QMLE and the Sandwich's Redemption

In practice we almost never dare to guarantee the distribution is fully correct. A very common situation is: you got the shape of the conditional mean right, but got the rest of the distribution (especially the variance) wrong. In that case the "likelihood" you wrote down is really a wrong likelihood, and maximizing it gives the quasi-maximum likelihood estimator (QMLE). The question arises: is QMLE still usable?

The answer splits in two, and this is exactly the most important distinction in the chapter. The point-estimate half is often still salvageable. Look at Section 3.5's condition that the score has expectation zero, $\mathbb{E}[s(z,\theta_0)] = 0$. If the score of your (wrong) likelihood happens to depend only on the part of the data you specified correctly (say the conditional mean), then even with the distribution wrong, this expectation is still zero at the true parameter, so Section 3.3's consistency argument still goes through, and QMLE is still consistent. This kind of likelihood, where "getting one part right is enough," belongs to the linear exponential family, which includes Poisson, normal, and binomial, whose scores all take the shape "residual times some function," so as long as the conditional mean is right, the residual has expectation zero.

The variance half cannot be salvaged, and if left unfixed it will deceive you. With the distribution misspecified, the premise "f is the true distribution" that the information matrix equality relies on is gone, $A = B$ no longer holds, and the sandwich no longer collapses. If you still stubbornly use $I^{-1}$ (that is, only one slice of bread $A^{-1}$) as the variance, you are computing a wrong variance. The correct variance returns to its original, general shape: the sandwich $A^{-1} B A^{-1}$, using the sample to estimate $A$ (the mean of the negative Hessian) and $B$ (the mean of the score outer product) separately, without assuming the two are equal. This robust variance is called the Huber-White sandwich in econometrics, and it is White's (1982) contribution.

::: {.intuition}
Remember this half-chapter in one sentence: get the conditional mean right and the point estimate lives; get the distribution right too and $I^{-1}$ lives. The two "rights" govern two different things, and conflating them is the beginner's most common mistake. Helix's usage data will play this out in gory detail: the real usage is overdispersed (variance far larger than the mean), yet the data scientist, taking the easy route, fits a Poisson likelihood. Poisson's conditional mean $\exp(x'\beta)$ is correctly specified, so her $\hat\beta$ is consistent and estimates hug the truth; but Poisson rigidly assumes the variance equals the mean, which is wrong for overdispersed data, and the information matrix equality breaks. If she trusts the default standard errors Poisson spits out (which are $I^{-1}$), she badly underestimates the uncertainty. Section 6 will measure this cost: a nominal 95% confidence interval actually covers the truth only 71.5% of the time, and swapping in the sandwich returns it to 95.5%. Same point estimate, two sets of standard errors, one deceives and one does not, and the difference lies entirely in whether you dare to guarantee the distribution.
:::

This half-chapter also rewrites the scope of Section 3.5's efficiency crown. MLE attains the Cramér-Rao lower bound on the premise that the likelihood is correct. Misspecify it and the crown comes off: QMLE is no longer efficient, and what it gains is robustness, trading a little efficiency loss for "even if the distribution is wrong, as long as the mean is right, both the point estimate and the (robust) inference still stand." When you have no full confidence in the distribution, this trade is almost always worth it, which is also why the default move of modern empirical work is to report sandwich standard errors rather than the model's default one.

### 3.7 GMM: Mounting Estimation Directly on Moment Conditions

MLE starts from a complete distributional assumption, which is both its strength (most efficient when correctly specified) and its soft spot (needing QMLE to patch it up when misspecified). Is there a path that bypasses the complete distribution and estimates using only the few moment conditions the economic model is confident about? Yes, and it is GMM. It is another pillar of this chapter's unified framework, and it also gathers all the earlier estimators into one larger pocket.

The starting point is a set of population moment conditions: the economic model tells you that at the true parameter, the expectation of some function is zero,

$$\mathbb{E}[g(w_i, \theta_0)] = 0,$$

where $g$ is a $q$-dimensional vector-valued function. The sample analogue is $\bar g_n(\theta) = \frac{1}{n}\sum_i g(w_i, \theta)$. If the number of moment conditions $q$ equals the number of parameters $k$ (just-identified), you can generally find a $\hat\theta$ that makes $\bar g_n(\hat\theta) = 0$ hold exactly. If there are more moment conditions than parameters (over-identified), you generally cannot find a $\theta$ satisfying all moment conditions at once, and you settle for the next best: find a $\theta$ that makes the moment conditions "overall closest to zero," using a $q \times q$ positive-definite weight matrix $W$ to weight their quadratic form:

$$\hat\theta_{GMM} = \arg\min_\theta\ \bar g_n(\theta)'\, W\, \bar g_n(\theta).$$

This is the GMM objective listed in Section 3.1. Its asymptotic normality can be plugged straight into Section 3.4's machine, only with $A$ and $B$ swapped into the language of moment conditions. Write the Jacobian $D = \mathbb{E}[\partial g(w_i,\theta_0)/\partial\theta']$ ($q \times k$, how sensitive the moment conditions are to the parameters) and the moment variance $S = \mathbb{E}[g(w_i,\theta_0)\,g(w_i,\theta_0)']$ ($q \times q$, how noisy the moment conditions themselves are). GMM's asymptotic variance is a larger sandwich:

$$\mathrm{Var} = (D'WD)^{-1}\,(D'W S W D)\,(D'WD)^{-1}.$$

The $D$ inside the bread $(D'WD)^{-1}$ is the moment-condition version of curvature, and the larger $D$ (the more sensitive the moment conditions are to $\theta$), the stronger the identification and the smaller the variance, which is exactly the moment-condition incarnation of Section 3.2's steep peak.

When over-identified, how to choose the weight $W$ has a beautiful optimal answer. Take $W = S^{-1}$, weighting by the inverse of the moment variance, so noisy moment conditions get small weight and steady ones get large weight. Substitute into the sandwich above, and $D'S^{-1}D$ lets the bread and filling cancel, so the variance collapses to

$$\mathrm{Var}_{\text{eff}} = (D'S^{-1}D)^{-1}.$$

This is the minimum variance efficient GMM can attain, formally just like MLE's $I^{-1}$, and the collapse mechanism is the same: choose the right weight so the two slices of bread of the sandwich equal the filling. $S$ is generally unknown, so in practice it is estimated in two steps: the first step uses a convenient weight (such as the identity $I$) to get a preliminary $\hat\theta$ and estimate $\hat S$ from it, and the second step re-estimates with $\hat W = \hat S^{-1}$, which is called two-step feasible GMM.

::: {.theorem}
**(Hansen's $J$ test.)** When over-identified ($q > k$), the value of the objective function estimated with the efficient weight, multiplied by $n$, converges to a chi-square distribution when all moment conditions are correct:

$$n\, \bar g_n(\hat\theta)'\, \hat S^{-1}\, \bar g_n(\hat\theta) \xrightarrow{d} \chi^2_{q-k}.$$

The degrees of freedom $q - k$ is the degree of over-identification, that is, the number of extra moment conditions.
:::

The intuition of the $J$ test is: when just-identified, you can drive the $k$ moment conditions exactly to zero, leaving no residual to test; when over-identified, the extra $q - k$ moment conditions cannot be pressed to zero even at the optimal $\hat\theta$, and if all moment conditions are true, this residual is just sampling noise and $J$ is small; if some moment condition is actually wrong (the model is misspecified), the residual systematically departs from zero, $J$ grows large, and it rejects. It is a checkup on the model's specification, and Section 6 will show it not rejecting when correctly specified and firmly rejecting when the mean is misspecified.

Finally, close the loop. MLE is really a special case of GMM: the score's expectation being zero, $\mathbb{E}[s(z,\theta_0)] = 0$, is itself a set of moment conditions, and moreover the "optimal" moment conditions (the variance estimated with them attains the Cramér-Rao lower bound). Section 2's logit score, once simplified, takes the shape $\sum_i (D_i - \Lambda(x_i'\beta))\, x_i$, which is exactly the moment condition "adoption residual times $x$," $\mathbb{E}[(D - \Lambda(x'\beta))\,x] = 0$. Section 6 will also show that Poisson's score takes the shape "usage residual times $x$," and adding one more moment on $x^2$ upgrades it from just-identified to over-identified GMM. MLE, QMLE, and GMM are thus gathered into one sentence: all are extremum estimators making some set of moment conditions (or their weighted quadratic form) hold as much as possible in the sample, and the only difference is where the moment conditions come from and whether that sandwich variance collapses.

The main points of this section can be summarized as follows: the extremum estimator is the extremum point of a sample objective function, and MLE, GMM, and NLS are all special cases; identification requires the population objective function to have a unique, separable peak, and the steepness of the peak is the curvature $A$ in the asymptotic variance; consistency is guaranteed by Jensen's inequality, independent of the specific form of $f$; the asymptotic variance generally takes the sandwich shape $A^{-1}BA^{-1}$, and MLE collapses it into $I^{-1}$ via the information matrix equality when correctly specified, attaining the Cramér-Rao lower bound; when half misspecified (mean right, distribution wrong), the QMLE point estimate is still consistent, but inference must use the sandwich rather than $I^{-1}$; GMM mounts estimation directly on moment conditions, the efficient weight $W = S^{-1}$ gives the minimum variance $(D'S^{-1}D)^{-1}$, over-identification allows a $J$ test to check the specification, and MLE is nothing but GMM with the optimal moment conditions.

## 4 Estimation: How to Actually Compute It

Identification and large-sample theory made clear why the estimator is reliable and what its variance looks like. This section covers the three practical matters of getting it done: how to find that extremum point numerically, how to do two-step feasible GMM, and how identification slides in practice from "strong" toward "weak."

### 4.1 Numerical Optimization

The logit log-likelihood has no closed-form solution, and $\hat\beta$ must be found by iteration. The workhorse algorithm is Newton-Raphson, whose idea is to approximate the objective function at the current point with a quadratic surface, jump directly to the vertex of that quadratic surface, and repeat. Written out with the gradient $\nabla \ell$ and Hessian $\nabla^2 \ell$, it is

$$\theta_{k+1} = \theta_k - \big[\nabla^2 \ell(\theta_k)\big]^{-1}\, \nabla \ell(\theta_k).$$

The intuition is: the gradient tells you which direction to go, the Hessian (curvature) tells you how far to go, taking small steps where the curvature is high and large steps where it is flat. When the objective function is concave (the logit log-likelihood is globally concave), Newton's method converges extremely fast, usually reaching the answer in a few steps. The hand-written Newton iteration in Section 6 converges on Helix in 6 steps to a solution matching mature software (differing at the order of $10^{-13}$).

Newton's method's soft spot is that each step must compute and invert the Hessian, which is costly when there are many parameters, and the Hessian may also be non-invertible (which is often a signal that identification is in trouble). Hence a series of variants. BHHH (Berndt-Hall-Hall-Hausman) is designed for MLE, replacing the negative Hessian with the mean of the score outer product $\frac{1}{n}\sum s_i s_i'$, which is exactly the substitution the information matrix equality permits (the two have equal expectation when correctly specified), and its advantage is needing only first derivatives and being naturally positive semidefinite. The quasi-Newton family (BFGS and others) gradually approximates the inverse of the Hessian during iteration using historical gradients, bypassing the explicit inversion at each step, and is the default choice for general nonlinear optimization. The shared caution is: if the objective function is nonconcave, or identification is so weak that the Hessian is nearly singular, any algorithm may get lost in a flat or multi-peaked terrain, failing to converge or halting at a false peak, and at that point what you should suspect is usually not the algorithm but the model.

### 4.2 Two-Step Feasible GMM and Efficiency

Section 3.7 said the efficient weight is $W = S^{-1}$, but $S$ is unknown, so getting it done means two steps. The first step takes a convenient initial weight $W_0$, commonly the identity $I$, to get a consistent but not necessarily efficient $\hat\theta^{(1)}$. Use it to construct an estimate of the moment variance, and here is a practical detail worth pointing out: use the demeaned form,

$$\hat S = \frac{1}{n}\sum_{i=1}^n \big(g(w_i,\hat\theta^{(1)}) - \bar g_n\big)\big(g(w_i,\hat\theta^{(1)}) - \bar g_n\big)',$$

rather than the un-demeaned $\frac{1}{n}\sum g g'$, because the latter mixes specification error into the weight when over-identified and is almost always a bad idea. The second step re-minimizes with $\hat W = \hat S^{-1}$ to get the efficient estimate $\hat\theta_{GMM}$, and then computes sandwich standard errors from it.

The precision gain that efficiency brings is real. In Section 6, on Helix's usage data, over-identified GMM with the identity weight gives a slope SE of 0.023, and switching to the efficient weight brings the SE down to 0.015, a reduction of a third. The reasoning is the same as in Section 3.7: the efficient weight turns noisy moment conditions down and steady ones up, which amounts to using the information more fully. Note that efficient GMM is asymptotically efficient, and in small samples the two-step $\hat S$ itself carries noise, so it can sometimes be less stable than a well-designed identity weight, which is also the motivation for improvements like continuously-updating GMM (Hansen, Heaton, and Yaron 1996).

### 4.3 Identification Weakens in Practice

Section 3.2 said identification comes in degrees, and that under weak identification the peak is flat and the variance explodes. This subsection pulls it from the abstract to the concrete, because it is one of the most hidden traps in estimation practice: everything in the report looks normal, the point estimate is there and the standard error is there too, yet that standard error is so large the estimate carries almost no information, and in more extreme cases even the usual normal approximation itself fails.

The root of weak identification is that the bread $A$ in the asymptotic variance (or $D$ in GMM) tends to zero. Return to Helix's usage model, where the slope $\beta_1$ relies on the covariation between engagement $x$ and usage, and the strength of identification is proportional to how large the variation in $x$ is. When the variation in $x$ is ample, the objective function is a steep, lone peak around $\beta_1$ and estimation is precise; when $x$ barely varies (say all firms' engagement is clustered together), the objective function spreads out into a plateau along the $\beta_1$ dimension, countless $\beta_1$ values give nearly the same objective value, and the estimator loses accuracy accordingly. Section 6 will place these two situations side by side: when the variation in $x$ is ample, the slope's sampling standard deviation is 0.017, and when the variation in $x$ is squeezed to a minimum it soars to 0.294, a full 17 times, while the sensitivity of the moment condition to $\beta_1$ behind it (the Jacobian $D$) differs by 567 times.

What needs honesty is the two faces of weak identification, neither exaggerated nor understated. In this chapter's example, the sample size is large enough ($n = 4000$), so weak identification mainly shows up as a collapse of precision: the estimator is still approximately normal, the robust standard error can still cover what it should, only the interval is so wide as to be useless (the slope's 95% interval spans from negative values to above 1, while the truth is a mere 0.5). This is "computable but useless," not "computed wrong." But in settings with small samples, or where the moment conditions are highly nonlinear in the parameters, weak identification will further pull the sampling distribution away from normal, and the usual $t$ interval cannot even hold its coverage, at which point you need specialized weak-identification-robust inference (Stock and Wright 2000), whose spirit is just like this series' chapter on instrumental variables handling weak instruments: when identification itself does not hold up, the right response is not to switch to a fancier estimator, but to admit that this data and this specification have little to say about this parameter.

The main points of Section 4 can be summarized as follows: MLE like logit, having no closed-form solution, is solved by iterating with Newton's method and its variants, converging extremely fast under a concave objective, and a nearly singular Hessian is often an identification alarm; efficient GMM is implemented as two steps, and using the demeaned $\hat S^{-1}$ as the weight delivers a real gain in precision; the strength of identification manifests as the curvature of the objective function, and weak identification shows up in large samples as a collapse of precision and in small samples as also destroying the normal approximation, and diagnosing it relies on thinking through where the parameter's identifying variation comes from, not on staring at the asterisks on some statistic.

## 5 Anchoring Papers

A method stands only when it lands in real research. Three anchoring papers: one turned discrete-choice MLE from a statistical trick into an estimation paradigm carrying economic structure, one laid down the entire asymptotic theory of GMM, and one is a representative work in information systems using discrete-choice MLE to study technology adoption. Each is organized by the five elements of paper, method, data, results, and limitations, with the emphasis on how they make the credibility of the specification clear.

### 5.1 McFadden (1974)

::: {.case}
Paper: "Conditional Logit Analysis of Qualitative Choice Behavior," in Zarembka (Ed.), Frontiers in Econometrics (Academic Press), pages 105 to 142. It is the founding work that welded the random utility model and conditional logit together and handed them to MLE for estimation, and the bloodline of Section 2's adoption model runs through here.

Method: McFadden starts from random utility, where an individual chooses the highest-utility option among several alternatives, and utility equals a linear combination of observable attributes plus an unobservable shock. He proves that when the shocks follow an independent, identically distributed extreme value (type I) distribution, the choice probability takes exactly the clean closed form of conditional logit, that is, the probability of choosing alternative $j$ is $P_j = \exp(x_j'\beta) / \sum_k \exp(x_k'\beta)$. The parameter $\beta$ is estimated by maximum likelihood, turning an economic model of individual choice into a practicable econometric estimator.

Data: mode-of-travel choices of San Francisco Bay Area commuters, using the travel time, cost, and other attributes of alternatives such as driving and public transit as explanatory variables, to predict which mode an individual chooses to commute by.

Results: the model can not only estimate the marginal effects of time and cost on choice and compute the substitution among modes, but also predict the mode share of a not-yet-built new mode (Bay Area Rapid Transit, BART), demonstrating the predictive power of structured discrete choice relative to purely descriptive regression. This method later became the cornerstone of demand estimation and discrete choice, and McFadden won the Nobel Prize in Economics in 2000 for it.

Limitations: conditional logit comes with the independence of irrelevant alternatives (IIA) property, arising from the assumption that shocks are independent and identically distributed, which constrains the pattern of substitution (the famous red-bus-blue-bus paradox) and is often violated in reality. Relaxing it requires generalizations such as nested logit and random-coefficients (mixed) logit, which is the main line this series develops later when discussing discrete choice and demand.
:::

The exemplary significance of this paper is that it demonstrates MLE is not just a statistical procedure, but a complete chain that translates a behaviorally interpretable economic model into a likelihood and then hands it to maximum likelihood for estimation. The credibility of the specification rests on the assumption about the shock distribution (extreme value gives logit), and the IIA constraint it brings is exactly the object that an entire subsequent literature set out to loosen.

### 5.2 Hansen (1982)

::: {.case}
Paper: "Large Sample Properties of Generalized Method of Moments Estimators," Econometrica 50(4):1029-1054. It laid down the entire asymptotic theory of GMM and is the original source of all the conclusions in this chapter's Section 3.7.

Method: Hansen builds estimation on a set of population moment conditions $\mathbb{E}[g(w,\theta_0)] = 0$ given by the economic model, without needing a complete distributional assumption. He proves the consistency and asymptotic normality of the estimator that minimizes the weighted quadratic form of the moment conditions under over-identification, gives the efficient weight $W = S^{-1}$ and the minimum variance $(D'S^{-1}D)^{-1}$ attained thereby, and proposes the over-identification test constructed from the objective function value (later called the Hansen $J$ test), with degrees of freedom equal to the degree of over-identification.

Data: the paper itself is methodological, and the classic application it drives is the Euler equation of asset pricing under rational expectations, using macro and financial time series (consumption growth, asset returns) to estimate utility parameters (such as the coefficient of relative risk aversion), models that naturally give moment conditions rather than a complete likelihood, which is where GMM comes into its own.

Results: GMM provides a unified estimation framework that gathers OLS, IV, and 2SLS as special cases, letting researchers still make efficient estimates and specification tests when they are only willing to maintain a few moment conditions and dare not write down a complete distribution. It thereafter became the workhorse tool of structural estimation (especially dynamic models and industrial organization).

Limitations: GMM's good properties are asymptotic, and the efficient weight $\hat S$ of the two-step estimate carries noise in small samples, which can distort inference; when the moment conditions are weak (the Jacobian $D$ is near zero) it falls into weak identification, and the $J$ test may also reject because of discrepancies among the moment conditions rather than a genuine model error. These are the practical vigilances this chapter's Section 4 and Section 7 stress repeatedly.
:::

Hansen's contribution pushes this chapter's unified framework to its extreme: estimation need not start from a distribution, and you can get to work as long as you can write down credible moment conditions. Its defense strategy is also clear, concentrating the credibility on a few economically grounded moment conditions rather than on a possibly misspecified complete likelihood, which shares one lineage with the robust spirit of Section 3.6's QMLE.

### 5.3 Forman (2005)

::: {.case}
Paper: "The Corporate Digital Divide: Determinants of Internet Adoption," Management Science 51(4):641-654. It is a representative work in information systems using discrete-choice MLE to study firm technology adoption, and it shares the same skeleton as Helix's adoption model: writing a binary adoption decision as the probability that a latent benefit exceeds a threshold, estimated by maximum likelihood.

Method: Forman models whether a firm adopts Internet technology as a binary discrete-choice problem, where adoption occurs if and only if the net benefit of adopting is positive, and uses discrete-choice maximum likelihood to estimate the dependence of the adoption probability on firm characteristics, existing IT investments, and organizational structure, which is exactly the binary-choice MLE of this chapter's Section 2. The core economic question is what kinds of firms adopt the Internet earlier, that is, what factors shape the so-called corporate digital divide.

Data: cross-sectional IT survey data at the establishment level for a large number of U.S. commercial establishments, covering firms of different sizes, geographic locations, and existing technology endowments across finance and services, observing whether they adopted Internet technology.

Results: the adoption probability rises with firm size and external competitive pressure, and larger firms facing fiercer competition adopt earlier. An interesting finding is that employees' geographic dispersion is complementary to Internet adoption, with more geographically dispersed firms more inclined to adopt, consistent with the logic that the Internet lowers internal coordination costs. Existing client/server network investments have two sides: on the one hand complementary to adoption, on the other hand proprietary, platform-locked legacy investments raise switching costs and instead slow adoption. These marginal effects are all read off from the coefficients of the binary-choice MLE, exactly the step from parameters to marginal effects discussed in this chapter's Section 2.

Limitations: the cross-sectional adoption probability model identifies correlation rather than causation, as adoption and firm characteristics may both be driven by unobserved organizational capability; binary choice only depicts whether adoption happens, not its depth or intensity; and the conclusions also depend on the specification of the shock distribution. To get closer to causation requires the designs of this series' causal chapters (panel, quasi-experiment, instrumental variables) to reinforce it.
:::

This paper's significance for the chapter is putting MLE back into a real IS research question: a binary adoption decision, a random-utility-style threshold model, a maximum likelihood estimate, plus an honest statement of scope (it estimates the determinants of and correlations with adoption, not the causal effect of an intervention). It demonstrates the most common use of discrete-choice MLE in information systems empirical work, and also reminds us of its boundaries.

The three together make the significance of anchoring clear: McFadden translates the economic model of discrete choice into a likelihood, Hansen frees estimation from distributional assumptions into moment conditions, and Forman applies this MLE machine to an adoption problem in information systems. They share this chapter's core: think through the model and the quantity to be estimated first, then choose a correct objective function, and finally attach an honest, robust measure of uncertainty to the result.

## 6 A Full Walkthrough on Helix Data

Now run the preceding tools through Helix from start to finish. The code below uses R 4.5.3 with the random seed fixed at 512, and every number cited in the text comes from the actual output of running this code.

### 6.1 The Data-Generating Process

The design parameters are as follows: $n = 4000$ firms, engagement $x_i \sim \mathcal{N}(0,1)$. Adoption is generated by logit, $\Pr(D_i = 1 \mid x_i) = \Lambda(a_0 + a_1 x_i)$, with $a_0 = -0.40$ and $a_1 = 1.10$. Usage is generated by a negative binomial, with a log-linear conditional mean $\mathbb{E}[y_i \mid x_i] = \exp(b_0 + b_1 x_i)$, with $b_0 = 1.00$ and $b_1 = 0.50$, while the size parameter (dispersion) is taken to be 2.0, making the variance $\mu + \mu^2/2$ far larger than the mean and manufacturing overdispersion.

```r
set.seed(512); n <- 4000
x  <- rnorm(n)                              # standardized engagement
D  <- rbinom(n, 1, plogis(-0.40 + 1.10 * x))# adoption: logit
mu <- exp(1.00 + 0.50 * x)                  # log-linear conditional mean of usage
y  <- rnbinom(n, size = 2.0, mu = mu)       # usage: negative binomial (overdispersed)
```

The key design intent is to make $y$ overdispersed but with a correct conditional mean: this way Poisson will get the mean right (QMLE consistent) yet get the variance assumption wrong (information matrix equality broken). In the sample, $y$'s variance divided by its mean is about 4.0, far larger than the 1.0 that Poisson imposes. The adoption rate is about 42.75%. In the truth, logit's AME is 0.216 and MEM is 0.264, and the true value of the usage slope is 0.50, and these numbers are the targets for all the estimates that follow.

### 6.2 The Bad Numbers of the LPM

```r
lpm  <- lm(D ~ x)
phat <- fitted(lpm)
```

The naive LPM gives a slope of 0.215 (SE 0.007), an intercept of 0.427. The real problem is in the fitted values: the maximum is 1.213 and the minimum is $-0.478$, and 2.93% of firms are predicted to have an adoption probability outside $[0,1]$ (2.5% below zero, 0.4% above one). The probability printed as 121% from Section 1 is confirmed here. The LPM's objective function (residual sum of squares) is mismatched to a binary outcome, and this pushes us toward the logit likelihood.

### 6.3 Logit MLE: Newton Iteration, Consistency, and Two Sets of Standard Errors

The logit log-likelihood is solved by Newton's method. The hand-written iteration agrees with the solution from mature software:

```r
g <- glm(D ~ x, family = binomial())        # logit MLE
```

$\hat\beta = (-0.367,\ 1.097)$, hugging the truth $(-0.40,\ 1.10)$; the hand-written Newton iteration converges in 6 steps, differing from `glm` at the order of $10^{-13}$. The AME read off from the parameters is 0.217, matching the true value 0.216. Because the logit likelihood is correctly specified here, the information matrix equality holds, and the two sets of standard errors should agree: the slope's Fisher-information standard error is 0.044, and the Huber-White sandwich standard error is also 0.044, agreeing as expected. This is Section 3.5 made manifest: with correct specification, the sandwich collapses into a single slice of bread, and the robust and non-robust standard errors are no different.

Consistency can be demonstrated directly. Fix the DGP and let the sample size grow from 100 to 102400, repeatedly sampling at each sample size to see the distribution of the slope estimate:

| $n$ | 100 | 400 | 1600 | 6400 | 25600 | 102400 |
|---|---|---|---|---|---|---|
| mean of $\hat\beta_1$ | 1.136 | 1.109 | 1.097 | 1.101 | 1.101 | 1.100 |
| SD of $\hat\beta_1$ | 0.279 | 0.147 | 0.071 | 0.035 | 0.017 | 0.009 |

The mean converges steadily to the truth 1.10, and the standard deviation roughly halves each time the sample size quadruples (times 4), which is exactly the $\sqrt{n}$ convergence rate: the variance shrinks as $1/n$ and the standard deviation shrinks as $1/\sqrt n$. The consistency that Jensen's inequality promised in Section 3.3 is delivered in the numbers.

![Left panel: the logit average log-likelihood as a function of $b_1$ is a smooth, lone peak attaining its maximum near the truth, with the summit being the MLE, which is the extremum in the extremum estimator. Right panel: slice the GMM objective function along $b_1$, and when engagement variation is ample (navy) it is a steep valley with strong identification, while when the variation is squeezed to a minimum (gold) it is nearly flattened out with weak identification, a flat objective function unable to tell which $\theta$ is better.](assets/fig/fig_05_identification.svg)

### 6.4 Poisson QMLE: Mean Right, Distribution Wrong, the Sandwich to the Rescue

Now it is usage's turn. The data scientist, taking the easy route, fits a Poisson likelihood to the overdispersed $y$:

```r
pm <- glm(y ~ x, family = poisson())        # Poisson QMLE
se_model <- sqrt(diag(vcov(pm)))            # Fisher-information standard errors
se_sand  <- sqrt(diag(sandwich::sandwich(pm)))# Huber-White sandwich
```

The point estimate is good: $\hat b = (1.024,\ 0.498)$, with the slope hugging the truth 0.50. This is the QMLE consistency of Section 3.6, where Poisson's conditional mean $\exp(x'\beta)$ is correctly specified, the score has expectation zero at the truth, and the consistency argument goes through unimpeded, even with the distributional assumption wrong.

The bad news is in the standard errors. Poisson rigidly assumes the variance equals the mean, while the data's Pearson dispersion is 2.63 (should be 1 under Poisson), so the information matrix equality breaks. The slope's Fisher-information standard error is 0.0090, and the sandwich standard error is 0.0182, the latter 2.0 times the former. Which is right? Use simulation as the referee: resample this DGP 2000 times and look at $\hat b_1$'s true sampling standard deviation, which is 0.0165. The sandwich's average of 0.0167 is almost dead on, while the Fisher-information's average of 0.0090 is only half the truth. The cost shows up directly in coverage: a nominal 95% confidence interval, using the Fisher-information standard error, actually covers the truth only 71.5% of the time, while using the sandwich it returns to 95.5%. Same point estimate, and trusting the wrong standard error turns an interval that should be 95% credible into one that is 71.5%.

![The Monte-Carlo sampling distribution of the Poisson-QMLE slope $\hat b_1$ (grey histogram, 2000 resamples), overlaid with the normal curves implied by each of the two sets of standard errors. The navy sandwich normal hugs the true sampling distribution, while the gold Fisher-information normal is far narrower, systematically underestimating the uncertainty by half, which is the geometric source of the nominal 95% interval actually covering only 71.5%.](assets/fig/fig_05_sandwich.svg)

### 6.5 The Same Thing, Done Again with GMM

Poisson's conditional mean gives a set of moment conditions $\mathbb{E}[(y - \exp(b_0 + b_1 x))\, h(x)] = 0$. Taking $h(x) = (1, x)$ is just-identified, and these two moment conditions are exactly Poisson's score first-order conditions, estimating $\hat b_1 = 0.498$, identical to the Poisson QMLE to the digit, confirming Section 3.7's "MLE is GMM with the optimal moment conditions."

Add one more moment, taking $h(x) = (1, x, x^2)$, and it upgrades from just-identified to over-identified (3 moments, 2 parameters). Now the weight matrix begins to matter. Estimating with the identity weight gives a slope of 0.509 (SE 0.023), and switching to the efficient weight $\hat W = \hat S^{-1}$ gives a slope of 0.484 (SE 0.015), the standard error reduced by a third, which is the precision dividend of efficient weighting. Over-identification also allows a Hansen $J$ test: here the conditional mean is correctly specified, $J = 1.88$ (1 degree of freedom), $p = 0.17$, does not reject, and the checkup passes.

To see whether the $J$ test really does raise an alarm, secretly add a quadratic term to usage's true mean, $\exp(b_0 + b_1 x + 0.25 x^2)$, while still fitting with the two-parameter model without $x^2$, and that $x^2$ moment condition should catch the flaw: $J$ jumps to 25.60, $p < 10^{-5}$, and firmly rejects. The $J$ test is quiet when correctly specified and flares up when the mean is misspecified, which is exactly its value as a specification checkup.

### 6.6 Weak Identification: A Flattened-Out Objective Function

Finally, demonstrate the weak identification of Section 4.3. Turn off the identifying variation of the usage slope: when engagement $x$ has ample variation (SD 1.0) it is strong identification, and when the variation is squeezed to a minimum (SD 0.05) it is weak identification, each resampled 1000 times:

| | strong identification | weak identification |
|---|---|---|
| sampling SD of $\hat b_1$ | 0.017 | 0.294 |
| moment-condition Jacobian $\lvert D_{b_1}\rvert$ | 3.85 | 0.007 |

Under weak identification the sampling standard deviation is 17 times that under strong identification, while the sensitivity of the moment condition to $b_1$ behind it (the Jacobian $D$, the signal strength of identification) differs by 567 times. The direction is exactly what the bread of Section 3.4's sandwich foretold: the weaker the identifying signal and the flatter the peak, the larger the variance. The two multiples are not in simple proportion, because the sandwich carries two slices of $A^{-1}$, so with fixed noise the variance roughly scales with the square of the reciprocal of the curvature, while squeezing $x$'s variation also shrinks the noise of the moment condition itself, and the two forces stack to land at a 17-fold gap in standard deviation. To be honest, this weak identification at $n = 4000$ is mainly a collapse of precision rather than a failure of inference, and the robust standard error still maintains about 94% coverage, only the slope's 95% interval is so wide it spans from negative values to above 1 while the truth is a mere 0.5, and the estimate is "computable but useless." This is the same pathology as the weak instruments in this series' chapter on instrumental variables, in two incarnations.

### 6.7 The Grand Reconciliation

Place all the estimates side by side, with each one's own truth as the target:

| estimator | estimate | SE | target / note |
|---|---|---|---|
| LPM slope | 0.215 | 0.007 | wrong objective function, predicted probabilities out of bounds |
| logit MLE slope | 1.097 | 0.044 | truth 1.10; Fisher = sandwich (correctly specified) |
| logit AME | 0.217 | none | truth 0.216 |
| Poisson QMLE slope | 0.498 | 0.018 | truth 0.50; sandwich, not the Fisher 0.009 |
| GMM just-identified slope | 0.498 | none | equal to Poisson QMLE |
| GMM efficient over-identified slope | 0.484 | 0.015 | efficient weighting, tighter than the identity-weight 0.023 |

This section's reconciliation can be summarized as follows: the LPM prints out-of-bounds probabilities because the objective function is chosen wrong; the logit MLE hugs the truth, and because it is correctly specified the two sets of standard errors agree, with consistency delivered as $\sqrt n$ with the sample size; the Poisson QMLE has a correctly specified mean so the point estimate is consistent, but the misspecified distribution makes the Fisher standard error report only half the true uncertainty, and only the sandwich fixes the coverage from 71.5% back to 95.5%; the same conditional mean rewritten as GMM has the just-identified case reproduce the QMLE and the over-identified case gain precision through efficient weighting and check the specification with a $J$ test; and draining the identifying variation flattens the objective function, so the variance amplifies a hundredfold with the collapse of the curvature.

## 7 Failure Modes and Robustness

The specifications in the simulation are manufactured, and in real research they can fail at any moment. This section reviews the most common ways they fail and the operational responses.

Misspecification is the most pervasive threat, and distinguishing "which part of the specification is wrong" is the whole key. Getting the conditional mean right and the distribution wrong is the luckiest tier: the point estimate is still consistent via QMLE, and all you need is to fix the standard errors with the sandwich, with Helix's Poisson example as the template. What is truly dangerous is getting even the conditional mean wrong, at which point consistency itself is gone, and no matter how robust the standard error, it only attaches an honest interval to a wrong point estimate, which is meaningless. Practice therefore has two lines of defense: first, report sandwich robust standard errors by default, unless you dare to guarantee the distribution is fully correct (a confidence rarely warranted), because robust standard errors are harmless when correctly specified (collapsing back to Fisher) and lifesaving when misspecified; second, use over-identification tests or specification tests to check the conditional mean itself, with the $J$ test, the Hausman test, and the information matrix test all being this kind of checkup, checking the mean layer rather than the distribution layer.

Reporting the wrong standard error is a closely adjacent second trap, and it is more hidden than a wrong point estimate, because the point estimate looks perfectly normal. The core discipline is to think through which slice the variance you report is: the model's default standard error is mostly $I^{-1}$ (one slice of bread), which is only right when the information matrix equality holds, and this equality requires the distribution to be fully correct plus homoskedasticity, no clustering, and no serial correlation, premises that often do not hold. The robust default should be the sandwich, and when the data has a group structure you should further cluster, which is the subject of the next chapter on inference. A recurring lesson is: a beautiful, suspiciously narrow standard error is often not because the estimate is precise, but because the wrong variance formula was used.

Weak identification is the third trap, and its destructive power has already been quantified above. Diagnosing it cannot rely only on the asterisks of the point estimate, but on where the identifying variation comes from and whether it is enough. In GMM, check whether the Jacobian $D$ is bizarrely close to singular; in MLE, check whether the Hessian is nearly singular, whether the log-likelihood is nearly flat, and whether the optimization struggles to converge, all of which are identification alarms. The first response to a remedy should not be to switch to a fancier estimator, but to return to the design: is the identifying variation for this parameter simply too thin. When the sample is small or the model is highly nonlinear, weak identification will also destroy the normal approximation, requiring specialized weak-identification-robust inference.

Numerical failure is equally worth watching, especially in discrete choice. A classic trap is perfect separation: if some covariate can perfectly separate $D = 1$ from $D = 0$, the logit maximum likelihood estimate diverges to infinity, with the coefficient estimated ever larger and the standard error ever more absurd, because the likelihood has no interior maximum along that direction. The software may throw no error, merely spitting out huge coefficients and standard errors, and it takes a person to see it. Similarly, when the objective function is nonconcave and multi-peaked, the optimization may halt at a false peak, and rerunning from different starting values and checking convergence diagnostics is basic hygiene. These are not theoretical problems, they are practical problems you run into every day.

The interpretation of tests also requires a measure of care. Asymptotically the three tests of Wald, likelihood ratio, and Lagrange multiplier (together the trinity of tests) are equivalent at the null hypothesis, but in finite samples or under suspect specification they part ways: the Wald test relies on the estimate of $I^{-1}$, which is distorted along with the standard error when misspecified; the likelihood ratio needs a correct likelihood and may not be usable under QMLE; and at that point a Wald test based on the robust variance or a specialized robust test is more trustworthy. When the $J$ test rejects you also need to think through what it is saying: it rejects "all moment conditions jointly hold," which could be that some moment is genuinely wrong, or that different moment conditions should give slightly different answers under heterogeneity, and not necessarily a disaster.

Stringing these failure modes together, this chapter's robustness rests, in the end, on two things: first, distinguishing which part of the model you got right (the mean? the distribution?), and judging from that whether the point estimate still lives and which set of standard errors to use; second, thinking through what variation your parameter is identified by and how strong that identification is. Think these two things through and estimation is stable; fail to think them through and even the most automated software is merely quickly handing you a number that may be wrong in a very respectable way.

## 8 Further Reading

::: {.readings}
Required reading, in suggested reading order:

- Wooldridge (2010, *Econometric Analysis of Cross Section and Panel Data*, 2nd ed.), Chapters 12 to 14. The clearest and most readable textbook version of a unified treatment of M-estimation, MLE, and GMM, and the source of this chapter's sandwich and QMLE framework.
- White (1982, *Econometrica*). The original source of QMLE and the sandwich, and reading it clarifies why the point estimate still lives when "the mean is right and the distribution is wrong," while $I^{-1}$ must be swapped for the sandwich.
- Hansen (1982, *Econometrica*). The founding work of GMM, focusing on how the efficient weight $W = S^{-1}$ collapses the variance into $(D'S^{-1}D)^{-1}$ and the construction of the over-identification $J$ test.
- McFadden (1974, *Frontiers in Econometrics*). The source of conditional logit and random utility, the template for translating the economic model of discrete choice into a likelihood, and the origin of this chapter's Section 2 adoption model.

Further reading:

- Newey and McFadden (1994, *Handbook of Econometrics* Vol. 4). The authoritative technical treatment of the consistency and asymptotic normality of extremum estimators, where the regularity conditions skipped in this chapter's 3.3 and 3.4 have complete proofs, the standard citation when writing methods papers.
- Amemiya (1985, *Advanced Econometrics*). The classic textbook source of the unified name and framework of the extremum estimator, gathering MLE, NLS, and GMM into one argmax viewpoint.
- Cameron and Trivedi (2013, *Regression Analysis of Count Data*, 2nd ed.). The standard reference for count models and Poisson QMLE, where the full discussion of overdispersion and the sandwich in this chapter's usage example resides.
- Hansen, Heaton and Yaron (1996, *JBES*). The finite-sample problems of two-step GMM and continuously-updating GMM, for understanding why "asymptotically efficient" need not be best in small samples.
- Stock and Wright (2000, *Econometrica*). Robust inference for GMM under weak identification, the response to the small-sample breakdown of weak identification in this chapter's Section 4.3, sharing one lineage with the weak instruments of the instrumental variables chapter.
:::

::: {.apa-refs}
- Amemiya, T. (1985). *Advanced econometrics*. Harvard University Press.
- Berndt, E. R., Hall, B. H., Hall, R. E., & Hausman, J. A. (1974). Estimation and inference in nonlinear structural models. *Annals of Economic and Social Measurement, 3*(4), 653-665.
- Cameron, A. C., & Trivedi, P. K. (2013). *Regression analysis of count data* (2nd ed.). Cambridge University Press.
- Forman, C. (2005). The corporate digital divide: Determinants of Internet adoption. *Management Science, 51*(4), 641-654. https://doi.org/10.1287/mnsc.1040.0343
- Hansen, L. P. (1982). Large sample properties of generalized method of moments estimators. *Econometrica, 50*(4), 1029-1054. https://doi.org/10.2307/1912775
- Hansen, L. P., Heaton, J., & Yaron, A. (1996). Finite-sample properties of some alternative GMM estimators. *Journal of Business & Economic Statistics, 14*(3), 262-280. https://doi.org/10.1080/07350015.1996.10524656
- McFadden, D. (1974). Conditional logit analysis of qualitative choice behavior. In P. Zarembka (Ed.), *Frontiers in econometrics* (pp. 105-142). Academic Press.
- Newey, W. K., & McFadden, D. (1994). Large sample estimation and hypothesis testing. In R. F. Engle & D. L. McFadden (Eds.), *Handbook of econometrics* (Vol. 4, pp. 2111-2245). Elsevier. https://doi.org/10.1016/S1573-4412(05)80005-4
- Stock, J. H., & Wright, J. H. (2000). GMM with weak identification. *Econometrica, 68*(5), 1055-1096. https://doi.org/10.1111/1468-0262.00151
- White, H. (1982). Maximum likelihood estimation of misspecified models. *Econometrica, 50*(1), 1-25. https://doi.org/10.2307/1912526
- Wooldridge, J. M. (2010). *Econometric analysis of cross section and panel data* (2nd ed.). MIT Press.
:::
