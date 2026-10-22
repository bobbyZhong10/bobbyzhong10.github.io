---
title: "Nonparametric & Semiparametric Methods"
subtitle: "Letting the Data Choose the Shape"
seriesline: "Foundations of Information Systems Economics · Chapter 14"
series: "Foundations of Information Systems Economics: Empirical Methods for Information Systems Economics"
sidebartitle: "Foundations of Information Systems Economics"
sidebarsub: "Theory · Methods · Applications"
chapterlabel: "Chapter 14 · Nonparametric & Semiparametric Methods"
---

## Introduction

Solstice's demand regression looks like a solid success: the linear fit has an $R^2$ as high as 0.859, and the constant-elasticity model estimates a price elasticity of $-1.584$. Taken at face value, that number says raising the price would lower revenue. Yet in the low-price region the firm actually cares about, the true local elasticity is only $-0.372$, which points to exactly the opposite recommendation; and at mid-range prices the elasticity jumps steeply to $-3.585$. The regression did not miscompute anything. It simply answered a local question with an average slope.

Some studies care about a single coefficient; others care precisely about how a curve bends. In the latter case, imposing a linear or constant-elasticity form up front can delete the answer before estimation even begins. Kernel regression, local linear fitting, series, and splines let the shape of the conditional mean be decided by the data, but flexibility has a price: a bandwidth that is too large smooths away the turns, one that is too small chases noise, and as the number of explanatory variables grows, even the nearest neighbors become far away. The bias-variance tradeoff, boundary correction, and the curse of dimensionality are all different line items on this bill.

This chapter also takes up a more practical compromise. If a researcher cares only about an advertising effect $\beta$ but does not want to impose an arbitrary form on the relation between demand and price, $g(x)$, the partially linear model can keep the former as a parameter, hand the latter to nonparametric methods, and recover the $\sqrt n$ rate through residualization. In the Solstice case, a straight line, a kernel estimate, a local linear fit, and a semiparametric estimate are reconciled on the same data. Letting the data choose the shape does not mean letting the data choose the question: the researcher must still state first whether the object of interest is the whole function, a derivative at some point, or a low-dimensional parameter net of flexible controls.

## 1 An Elasticity That Lies

::: {.case}
Solstice is a fictional e-commerce platform. We tell this case with simulated data, and the advantage is that the data-generating process is fully known, so the estimated curve can be reconciled point by point against the true curve, a teaching condition that real data can never provide. Across a large number of products, the platform observes the price $x$ and the log demand $y$ that follows. The true demand curve is nonlinear: demand barely moves with price in the low-price region (customers are price-insensitive), drops sharply near a mid-range price (the sensitive region), and flattens out again at high prices. Management asks a direct operational question: should the price go up now or not.

The analyst uses the most standard approach, regressing log demand on log price to estimate a constant-elasticity model, and obtains a single price elasticity $\hat\eta = -1.584$. This number is greater than 1 in absolute value, meaning demand is elastic and a price increase would cut sales by a larger proportion, which yields the recommendation: do not raise the price. The linear benchmark, meanwhile, has an $R^2$ as high as 0.859, which looks quite respectable, and the analyst is ready to report it as is.
:::

That $-1.584$ is the first bad number of this chapter, and it is bad because it compresses a curve rich in shape into a single figure. Because this is a simulation, we know the true elasticity is not constant at all; it varies dramatically with price. Near the current low price (price around 2), the true local elasticity is $-0.372$, far below 1, so demand is actually rather inelastic there, sales barely fall when the price rises, and this is a price that should go up. In the mid-range sensitive region (price around 5.5), the true elasticity jumps steeply to $-3.585$, and this is where the price truly cannot rise. At high prices (price around 8) it falls back to $-1.487$. A quantity that spans a tenfold range from $-0.37$ to $-3.59$ is reported by that constant-elasticity model as a single blanket $-1.584$. And the analyst's product happens to be priced in the low-price region, where the true elasticity $-0.37$ is plainly calling for a price increase, yet that $-1.584$ keeps her standing pat, giving exactly the opposite recommendation.

That respectable $R^2 = 0.859$ is also deceptive. Take the residuals of the linear fit and look at them: they are not randomly scattered but carry an obvious systematic curvature, and a RESET-type specification test returns $p = 2 \times 10^{-20}$, decisively rejecting "linearity is enough." The data are shouting through the structure in the residuals: the true curve is bent, and your straight line failed to keep up.

::: {.intuition}
The constant-elasticity model's error is not that it estimates imprecisely but that it asks the wrong question. It assumes by default that price sensitivity is a universal constant, and so it is forced to average the sluggishness of the low-price region, the violence of the sensitive region, and the moderation of the high-price region all into a single number. This is like using a "national average temperature" to decide what each city should wear today: the number itself is not wrong, yet it is bad advice for any specific city. To answer "should the price go up here," what you need is not an average elasticity but the whole curve of how elasticity varies with price, that is, the local slope of the demand curve at every point. To estimate that curve, you can no longer presuppose its shape; you have to let the data speak the shape themselves. This is exactly what nonparametric regression sets out to do.
:::

The lesson of Section 1 can be summarized as follows: when the shape itself is the answer to the question, presupposing a functional form (linear, constant-elasticity) smooths the answer away in advance and reports an average number that may mislead every specific decision; a respectable $R^2$ cannot conceal the systematic curvature in the residuals, and to recover the truth of how elasticity varies with price, you need nonparametric regression that does not presuppose the shape.

## 2 The Economic Model and the Estimand

Nonparametric methods change not the mechanics of estimation but its target. This section makes clear what that target is, and how semiparametrics compromise between flexibility and precision.

In the parametric world, the estimand is a finite set of numbers, such as the elasticity $\eta$ in the constant-elasticity model (estimated by $\hat\eta$). In the nonparametric world, the estimand is a whole curve: the conditional expectation function $m(x) = \mathbb{E}[y \mid x]$ itself, together with the local slope $m'(x)$ derived from it (in the demand example, this is the price sensitivity that varies with price; multiplied by price, it is the local elasticity). We do not write a form like $m(x) = \alpha + \beta x$; we assume only that $m$ is a sufficiently smooth curve and leave its actual shape to the data. This is the whole claim of nonparametrics, and it is the common source of its power and its cost: if you do not presuppose the shape, you cannot be biased by choosing the wrong one, but because you must learn the whole curve from the data, you need far more information than you would to estimate a few parameters.

Real empirical questions often do not need to go all the way to pure nonparametrics. Beyond the demand curve, Solstice also wants to know one more thing: the effect of advertising on demand. Suppose each product has, besides its price $x$, an advertising exposure $z$ (in this chapter, lowercase $z$ denotes advertising exposure, a continuous regressor; it is unrelated to the instruments written $z_i$, $\mathbf{z}_i$ or $Z_i$ in other chapters), and management cares about the coefficient $\beta$ on $z$ in the demand equation. The trouble is that advertising exposure and price are correlated (the platform advertises with different intensity at different price levels, and this relationship is itself nonlinear). If you simply regress demand on $z$ and ignore price, $\beta$ will be contaminated by confounding from price; if you instead control for a linear price term, the control is incomplete, because the true effect of price is nonlinear. When you do not want to force a possibly wrong functional form on price and care only about the single parameter $\beta$, the semiparametric partially linear model is exactly the right tool:

$$y_i = \beta\, z_i + g(x_i) + \varepsilon_i,$$

which lets the price effect $g(x)$ bend freely and keeps only the advertising effect $\beta$ as a parameter. Here the estimand is the single number $\beta$, while $g$ is a nuisance function that must be handled flexibly but need not be interpreted. Section 4.3 will show that this division of labor ("parameterize what you care about, nonparameterize the nuisance") lets $\beta$ be estimated at the fast $\sqrt n$ rate, escaping the slow convergence and the curse of dimensionality of pure nonparametrics.

The main points of Section 2 are summarized as follows: the nonparametric target is the whole curve $m(x) = \mathbb{E}[y \mid x]$ and its local slope, at the cost of learning the shape from the data and converging more slowly; the semiparametric partially linear model keeps the coefficient of interest $\beta$ as a parameter and hands the nuisance function $g$ to nonparametrics, so the estimand is again a single number; the aim is to have both flexibility and precision.

## 3 The Core of Nonparametric Estimation

This is the most finely analyzed section of the chapter. It starts from the most naive local average (3.1), draws out the bias-variance tradeoff, which is the central tension of nonparametrics (3.2), then explains how local linear fitting corrects the boundary bias of kernel regression (3.3), the parallel route of series and splines (3.4), and finally why the curse of dimensionality forces us toward semiparametrics (3.5).

### 3.1 Local Averaging: From kNN to Kernels

To estimate the conditional expectation $m(x_0) = \mathbb{E}[y \mid x = x_0]$ at some price $x_0$, the most naive idea is to look at the products whose price is near $x_0$ and average their demand. This idea of "using the average of the neighbors as the estimate" is the seed of all nonparametric regression.

The most direct implementation is k-nearest-neighbors (kNN): take the $k$ observations closest to $x_0$, average their $y$, $\hat m(x_0) = \frac{1}{k}\sum_{i \in N_k(x_0)} y_i$. Here $k$ is the smoothing knob: small $k$ looks only at very close neighbors, so the estimate hugs the local behavior but jitters badly; large $k$ takes in far neighbors, so the estimate is steady but averages in irrelevant information too. kNN has two flaws: the membership of the neighborhood changes by jumps as $x_0$ moves, so the estimated curve is rough and discontinuous; and it treats every neighbor alike, whether it sits right next to $x_0$ or barely reaches it.

The kernel method fixes both flaws at once. It no longer uses the hard boundary of "the closest $k$"; instead it gives each observation a weight that decays smoothly with distance, larger the closer it is. Using a kernel function $K$ (a symmetric, integrate-to-one bell-shaped function, such as the Gaussian density) and a bandwidth $h$ to set the weights, you get the Nadaraya-Watson kernel regression estimator:

$$\hat m(x_0) = \frac{\sum_{i=1}^n y_i\, K\!\big(\frac{x_0 - x_i}{h}\big)}{\sum_{i=1}^n K\!\big(\frac{x_0 - x_i}{h}\big)}.$$

It is just a weighted average: the numerator is the sum of weights times $y$, the denominator is the sum of the weights, and the weights are given by the kernel function according to distance from $x_0$. The bandwidth $h$ replaces kNN's $k$ as the smoothing knob: small $h$ means only very close points have appreciable weight (local, jittery), large $h$ means even far points have weight (smooth, dull). The specific shape of the kernel function (Gaussian or Epanechnikov) matters little; what really decides everything is the bandwidth $h$.

::: {.intuition}
Picture kernel regression as sweeping a focusable flashlight along the price axis. The flashlight lights up a ring of products around $x_0$, brighter the closer they are (larger weight), and you take the weighted average demand of this lit-up ring of products as the demand estimate at $x_0$. The bandwidth $h$ is the size of the aperture: too small an aperture lights only a handful of products, the average is dominated by the random noise of those few, and as you sweep along the price axis the estimate bounces up and down; too large an aperture lights up half the price axis, averaging together the sluggishness of the low-price region and the violence of the sensitive region, and the true bend is smoothed away. The whole craft of nonparametric regression is tuning this aperture to be neither too large nor too small, and that is the bias-variance tradeoff of the next section.
:::

### 3.2 The Bias-Variance Tradeoff: The Dilemma of Smoothing

How large the bandwidth should be has no once-and-for-all answer, because it faces a tradeoff that gives with one hand and takes with the other, the central tension of nonparametric estimation, and the home turf of the bias-variance decomposition here.

Split the mean squared error of the kernel estimator into squared bias plus variance, and the two pieces respond to the bandwidth $h$ in exactly opposite directions. The bias comes from "using the average of a ring of neighbors to represent the center point": as long as the true curve $m$ is bent over that ring of neighbors, the average will systematically depart from the true value at the center, and the more the curve bends (the larger $m''$) and the larger the aperture (the larger $h$), the larger this bias, of order $h^2$,

$$\text{Bias}\big(\hat m(x_0)\big) \approx \frac{h^2}{2}\, m''(x_0) \int u^2 K(u)\, du.$$

The variance comes from "there are only finitely many points in the ring, and their noise has not been averaged away cleanly": there are effectively about $nh$ points in the aperture taking part in the average, and the more points there are the flatter the noise is pressed down, so the variance falls inversely with $nh$, and the larger the aperture (the larger $h$) the smaller the variance,

$$\text{Var}\big(\hat m(x_0)\big) \approx \frac{1}{nh}\, \frac{\sigma^2}{f(x_0)} \int K(u)^2\, du,$$

where $\sigma^2 = \mathrm{Var}(y \mid x)$ is the conditional variance (taken here to be homoskedastic) and $f(x_0)$ is the density of price at $x_0$; the variance is large where the density is low (where data are sparse). Putting the two pieces together, the mean squared error is $\text{MSE} \approx A h^4 + B/(nh)$: one term rises with $h$ and the other falls with $h$, which gives a U-shaped curve.

::: {.warning}
This inverse relationship between bias and variance means there is no bandwidth that is "best of both worlds," only a "least bad" compromise. Minimizing the U-shaped MSE over $h$, the optimal bandwidth is $h^* \propto n^{-1/5}$, shrinking slowly with sample size, and substituting it back gives an optimal MSE convergence rate of $n^{-4/5}$. This $n^{-4/5}$ is the nonparametric fate, and it is slower than the parametric model's $1/n$: to reach the same estimation precision, nonparametrics needs far more data. This is not a matter of the method being inadequate but the cost of "not presupposing the shape" itself; if you want the data to speak the whole curve, you must pay with convergence far slower than estimating just a few parameters. Bandwidth selection is therefore the most crucial and least fudgeable step in nonparametric estimation, and Section 4 explains how to get it right with cross-validation.
:::

### 3.3 Local Linear: Correcting Boundary Bias

Nadaraya-Watson kernel regression has a hidden defect, especially pronounced at boundaries and steep places. It essentially fits the neighborhood with a horizontal line (local constant) at every point, but when the true curve has slope at that point, fitting a horizontal line to a ring of sloped points makes the average lean toward the side with more neighbors. In the interior of the sample this bias can still be roughly canceled by the neighbors on the two sides, but at the boundary of the price axis, one side has no neighbors at all, the bias has nothing to cancel it, and kernel regression systematically drags the estimate inward at the boundary.

The correction is surprisingly natural: at every point, instead of a horizontal line, fit locally with a sloped line. This is local linear regression, which at each $x_0$ does a weighted least squares, with the weights still given by the kernel function, but fitting a line with slope:

$$\min_{a, b}\ \sum_{i=1}^n K\!\Big(\frac{x_0 - x_i}{h}\Big)\,\big(y_i - a - b\,(x_i - x_0)\big)^2,$$

the estimated intercept $\hat a$ is $\hat m(x_0)$, and the slope $\hat b$ is, as a bonus, the derivative $\hat m'(x_0)$. This is Chapter 5's M-estimation framework applied locally, with a squared-error loss weighted by the kernel, which merely supplies distance-decaying weights. But the kernel weights change as $h \to 0$, so the consistency and $\sqrt n$ asymptotic normality of Theorems 2 and 3 in Chapter 5 do not carry over directly. At a fixed $h$, the estimator targets a smoothed version of the curve, and this is exactly where the bias comes from; the rate of $\hat m(x_0) = \hat a$ is $\sqrt{nh}$ (Section 4.2), and the slope $\hat b$ converges more slowly, at $\sqrt{nh^3}$. Giving the local fit a slope of its own lets it follow the tilt of the true curve at that point, and the boundary bias disappears. The cost is nearly negligible. Local linear keeps the bias of order $h^2$ at both the boundary and the interior, whereas kernel regression not only degrades to order $h$ (worse) at the boundary but also has, in its interior bias of order $h^2$, an extra term related to the slope of the data density, an $f'/f$ term (the clean $m''$ bias of Section 3.2 is in fact the local linear form); local linear removes this term as well. It also comes with a bonus: the slope $\hat b$ directly estimates exactly the derivative that the local elasticity needs. In practice, local linear almost always beats Nadaraya-Watson and is the default choice for kernel regression.

::: {.intuition}
Boundary bias can be seen this way. Suppose the true demand curve slopes downward at the far left of the price axis (the lowest price). When kernel regression estimates at the far left, all the neighbors in the aperture are to its right, at higher prices and lower demand, so the average demand of these neighbors is dragged down by a string of lower values on the right side, and kernel regression estimates the far-left demand too low. Local linear does not fall for this: it sees that this ring of neighbors lines up along a downhill sloped line, and it extrapolates back to the far left along that line rather than taking their downward-biased average. In a word, kernel regression at the boundary "only takes averages," while local linear "reads the trend," and the boundary is exactly where the trend is more trustworthy than the average. Section 6 uses a simulation to measure this boundary bias: kernel regression has an obvious systematic bias at both ends, local linear almost zero.
:::

### 3.4 Series and Splines: Another Route

Local averaging starts from "look at the neighbors of each point," but there is a completely different route that arrives at the same destination: instead of looking locally, use a combination of a family of simple functions to approximate the whole curve. This is series or sieve estimation.

The idea is to write the unknown curve $m(x)$ as a linear combination of a set of basis functions, $m(x) \approx \sum_{j=1}^{K_n} \theta_j\, p_j(x)$, where the $p_j$ are a pre-chosen basis (polynomials $1, x, x^2, \ldots$, or trigonometric functions, or the splines discussed below), and then run ordinary least squares on this model (which is linear in the coefficients $\theta_j$). The number of basis functions $K_n$ is the smoothing knob here, playing a role dual to the bandwidth $h$: the larger $K_n$, the richer the basis and the more curved the functions it can approximate (small bias), but also the more coefficients there are to estimate and the more easily the fit chases noise (large variance). The essence of the sieve is to let $K_n$ grow with the sample size $n$, but more slowly than $n$, so that bias and variance go to zero together.

Using high-degree polynomials directly as the basis has a well-known flaw: a global polynomial oscillates wildly elsewhere in order to accommodate a bend in one place (the Runge phenomenon), and this gets especially out of hand at high degree. Splines are the more stable choice. They cut the price axis into small segments with a number of knots, use a low-degree polynomial (usually cubic) on each segment, and require adjacent segments to join smoothly at the knots (the function value, first derivative, and second derivative are all continuous). This gives both the flexibility of the pieces (able to hug local bends) and global smoothness (no oscillation). The number and position of knots control the degree of flexibility, more knots meaning more flexible, and this is again that bias-variance knob. An elegant variant is the smoothing spline, which does not directly choose the number of knots but instead minimizes "fitting error plus a penalty on curvature," $\sum_i (y_i - m(x_i))^2 + \lambda \int m''(x)^2 dx$, using the penalty coefficient $\lambda$ to tune smoothness continuously, large $\lambda$ leaning toward smoothness (small curvature), small $\lambda$ leaning toward hugging the data.

Series and spline methods on the one hand and kernel methods on the other each have their strengths, but they share the same bias-variance structure: kernel methods tune smoothness with the bandwidth $h$, and series with the number of terms $K_n$ (or splines with the number of knots, and smoothing splines with the penalty $\lambda$). The knobs run in opposite directions (a large $h$ corresponds to a small $K_n$, and both mean more smoothing), yet they face the same U-shaped tradeoff. Whichever route you take, you cannot avoid the same question in the end: where to set this smoothing knob. That is the subject of cross-validation in Section 4.

### 3.5 The Curse of Dimensionality

So far there has been only one regressor (price). Once the regressors multiply, nonparametric methods hit a wall called the curse of dimensionality, and it is the key to understanding why semiparametrics is needed.

The problem is that the very concept of "local" breaks down in high dimensions. Kernel methods and kNN both rely on "taking a ring of neighbors near the target point," but in high-dimensional space points are on average far from one another, and there simply are not enough neighbors nearby. A simple calculation quantifies this: to enclose within a $d$-dimensional unit cube a small cube containing $10\%$ of the data, how long must each edge of this small cube be? The answer is $e_d(0.1) = 0.1^{1/d}$. In one dimension it is 0.1 (enclosing a tenth of the price interval suffices), in two it is 0.32, in five it is 0.63, in ten it is 0.79. That is, in ten dimensions you must let the neighborhood span nearly $80\%$ of the coordinate range on every single axis just to hold a mere $10\%$ of the data. Such a "neighborhood" is no longer local at all, and the foundation of nonparametrics collapses.

Reflected in the convergence rate, the optimal MSE rate of $d$-dimensional kernel regression is $n^{-4/(4+d)}$, growing rapidly slower with dimension: one dimension is $n^{-4/5}$, still acceptable; by five dimensions it drops to $n^{-4/9}$; in ten dimensions it barely moves. Put differently, to reach in high dimensions the same precision as in low dimensions, the required sample size explodes exponentially with dimension. Silverman has a famous table: to estimate the value of a standard normal density at the origin to a fixed precision, one dimension needs just 4 observations, five dimensions need 768, and ten dimensions need 840,000. This is the fearsome face of the curse of dimensionality: purely nonparametric methods are basically unusable beyond three or four regressors.

::: {.intuition}
The curse of dimensionality forces a pragmatic division of labor. If you have ten regressors but truly care about the effect of only one, with the other nine being nuisances that merely need controlling, then nonparameterizing the entire ten-dimensional surface is both infeasible and unnecessary. The smarter move is to loosen the structure just a little: assume those nine nuisances enter in some simple way (say additively, or only through a single linear index), bring the dimension back down to something manageable, and keep the nonparametric treatment only where flexibility is truly needed. This is exactly the idea of semiparametric methods, and it is what the partially linear model of the next section delivers: let the parameter you care about enjoy fast $\sqrt n$ convergence, and hand only a one- or two-dimensional nuisance function to nonparametrics, thereby getting around the curse.
:::

The main points of this section are summarized as follows: the seed of nonparametric regression is local averaging, from kNN to kernel regression (Nadaraya-Watson is a kernel-weighted average), and the bandwidth is the smoothing knob; smoothing faces a U-shaped tradeoff between bias ($h^2 m''$) and variance ($1/nh$), with optimal bandwidth $h^* \propto n^{-1/5}$ and convergence rate $n^{-4/5}$ slower than the parametric $1/n$; local linear replaces kernel regression's local constant with a local sloped line, fixing boundary bias and throwing in the derivative for free; series and splines are a parallel route, tuning the same bias-variance knob with the number of terms or knots; the curse of dimensionality makes pure nonparametrics fail with several regressors, forcing us toward semiparametrics.

## 4 Estimation: Tuning Smoothness, Doing Inference, Semiparametrics

Identification and theory made clear what nonparametrics estimates and what tradeoff it faces. This section covers three things that land in practice: how to tune the smoothing knob right (4.1), how to attach an honest measure of uncertainty to a nonparametric estimate (4.2), and how semiparametrics escapes the curse of dimensionality (4.3).

### 4.1 Cross-Validation: Let the Data Choose the Smoothness

The bandwidth or number of knots is the make-or-break of nonparametric estimation, but it cannot be guessed by eye; the data must choose it themselves. The main tool is cross-validation.

The idea is to estimate directly "how large this bandwidth's prediction error is on new data," then choose the bandwidth with the smallest prediction error. Leave-one-out cross-validation takes this to the extreme: for each observation $i$, estimate the curve using the other $n-1$ points, then see how accurately it predicts the held-out $i$-th point, sum the squared prediction errors over all held-out points as this bandwidth's score, $CV(h) = \sum_i (y_i - \hat m_{(-i)}(x_i))^2$, sweep across bandwidths, and take the lowest score. Leave-one-out is necessary because if you do not hold out, using the same point both to estimate and to check, a small bandwidth will always look perfect (the curve passes through every point), and cross-validation punctures this illusion of overfitting by "predicting points that did not take part in estimation." Section 6 will show that on Solstice, the bandwidth chosen by cross-validation, 0.372, is very close to the bandwidth that minimizes the true mean squared error, 0.40, the data having successfully tuned the aperture to be neither too large nor too small.

Two practical points complement bandwidth selection. One is a quick rule of thumb, such as Silverman's $h \approx 1.06\,\hat\sigma_x\, n^{-1/5}$ (where $\hat\sigma_x$ is the sample standard deviation of price $x$), which is suitable for a first exploratory look, though formal estimation should still use cross-validation. The other is undersmoothing for inference. The optimal bandwidth makes the bias and the standard deviation of the same order of magnitude, but in a confidence interval that bias pushes the interval off center; so for inference one often deliberately sets the bandwidth a little below the optimum, suppressing the bias in exchange for honest coverage. This leads into the next subsection.

### 4.2 Inference for Nonparametrics

Attaching a standard error to a nonparametric estimate is trickier than in the parametric case, because of the added specter of bias. A pointwise confidence interval looks roughly like $\hat m(x_0) \pm z_{1-\alpha/2}\cdot \widehat{\text{sd}}$ (where $z_{1-\alpha/2}$ is the standard normal quantile and $\alpha$ is the significance level). The standard deviation is of order $1/\sqrt{nh}$ (recall that the variance in Section 3.2 is $\propto 1/nh$): the effective sample size is the $nh$ points inside the aperture rather than all $n$, so a nonparametric interval is inherently wider than a parametric one.

The real trouble is bias. At the optimal bandwidth, the bias and the standard deviation of a nonparametric estimate are of the same order of magnitude. This means that $\hat m(x_0)$ deviates systematically from the true value by about half a standard deviation (at the MSE-optimal bandwidth, the squared bias is about a quarter of the variance), so a symmetric interval centered at $\hat m(x_0)$ is shifted as a whole and its coverage falls short of the nominal level. There are two remedies. One is undersmoothing: setting the bandwidth below the optimum makes the bias go to zero faster, so that it becomes negligible relative to the standard deviation and the interval is honest (at the cost of a slightly larger variance and a slightly wider interval). The other is to estimate the bias explicitly and correct for it. Moreover, the analytic variance of a nonparametric estimator is often hard to write down, so in practice standard errors often come from the bootstrap of Chapter 8 (resampling with the undersmoothed estimator). The point is that a nonparametric confidence interval cannot simply copy the parametric approach: it must confront the bias, which is of the same order as the standard deviation, or the interval will lie.

### 4.3 Semiparametrics: Partially Linear and Robinson

Now we deliver on the promise of Section 3.5: how to escape the curse of dimensionality when you care about only one coefficient and the rest are nuisances that need flexible control. The answer is the partially linear model $y_i = \beta z_i + g(x_i) + \varepsilon_i$, which keeps the $\beta$ you care about as a parameter and hands the nuisance $g$ to nonparametrics.

Robinson (1988) proposed an elegant two-step estimator built on a "double-residual" trick, which is essentially a nonparametric version of the "partial out" idea of the Frisch-Waugh-Lovell theorem. Taking expectations conditional on $x$ on both sides of the partially linear model gives $\mathbb{E}[y \mid x] = \beta\, \mathbb{E}[z \mid x] + g(x)$, and subtracting the two equations eliminates $g(x)$:

$$y_i - \mathbb{E}[y_i \mid x_i] = \beta\,\big(z_i - \mathbb{E}[z_i \mid x_i]\big) + \varepsilon_i.$$

The term $g(x)$ is gone. So once you estimate the two conditional expectations $\mathbb{E}[y \mid x]$ and $\mathbb{E}[z \mid x]$ nonparametrically (with kernels or splines), form the residuals $\hat e_{y} = y - \hat{\mathbb{E}}[y\mid x]$ and $\hat e_{z} = z - \hat{\mathbb{E}}[z\mid x]$, and regress $\hat e_y$ on $\hat e_z$ by ordinary least squares, the slope is $\hat\beta$. The intuition is that, once you strip out of both $y$ and $z$, nonparametrically, everything that price $x$ can explain, the confounding by price (however curved $g$ may be) has been washed out of the remaining residual variation, so the relationship between the two residuals is the pure advertising effect $\beta$.

::: {.theorem}
**(Robinson's $\sqrt n$ convergence for the partially linear model)** In the partially linear model, assume that $\mathbb{E}[\varepsilon \mid x, z] = 0$, that $z$ is not a function of $x$ (that is, $\mathbb{E}[e_z^2] > 0$ with $e_z = z - \mathbb{E}[z \mid x]$), and that the estimation errors of the first-step estimates of $\mathbb{E}[y \mid x]$ and $\mathbb{E}[z \mid x]$ vanish faster than $n^{-1/4}$ (Robinson guarantees this through smoothness, bandwidth and trimming conditions). Then, although the nuisance function $g$ can be estimated only at the slow nonparametric rate, the double-residual estimator $\hat\beta$ still converges at the parametric $\sqrt n$ rate and is asymptotically normal:

$$\sqrt{n}\,(\hat\beta - \beta) \xrightarrow{d} \mathcal{N}(0, V).$$

Here $V = \mathbb{E}[\varepsilon^2 e_z^2]/(\mathbb{E}[e_z^2])^2$, which reduces to $\sigma^2/\mathbb{E}[e_z^2]$ under homoskedasticity ($\mathbb{E}[\varepsilon^2 \mid x, z] = \sigma^2$).
:::

This conclusion is close to magical, and it is the core selling point of semiparametric methods: $g$ is estimated slowly ($n^{-2/5}$), yet $\beta$ is estimated fast ($n^{-1/2}$); the slow part does not drag down the fast part. The intuitive reason is that the errors of the two first-step nonparametric estimates cancel to higher order in the subtraction and averaging, so they do not affect the first-order asymptotics of $\beta$. This is the nonparametric version of Theorem 6 (i) in Section 3.6 of Chapter 5 ($\Gamma_0 = 0$): the derivative of the second-step score with respect to the two conditional expectations has expectation zero. The first step, however, is nonparametric and does not satisfy TS.1 of that chapter, so the result relies on the extension in Newey (1994) and on the condition in the theorem box that the first step converges faster than $n^{-1/4}$. So you have neither forced a possibly wrong functional form on the price effect ($g$ bends freely) nor lost parametric precision for the advertising coefficient, and the curse of dimensionality is sidestepped. Section 6 will show that under repeated sampling, the estimate that controls for price linearly averages 0.758 (systematically low), while the Robinson double-residual estimate averages 0.800, right on the true value.

A neighboring family of semiparametric ideas is the single-index model, which compresses several regressors (here $x$ is a vector) into a single linear index $x'\gamma$ that then enters the outcome through an unknown nonparametric link function, $m(x) = F(x'\gamma)$. It is likewise a way to relax the structure and bring a high dimension back down to a low one; logit and probit (Section 2.1 of Chapter 12, where the index coefficient is written $\beta$; in this chapter $\beta$ is already the advertising coefficient) are special cases for particular choices of $F$. Such models share the spirit of the partially linear model: keep the nonparametric treatment only for the one dimension that truly needs flexibility, and use parametric structure to hold down the dimension of the rest.

The main points of this section are summarized as follows: cross-validation sets the smoothing knob at the bottom of the bias-variance U by "predicting held-out points," and on Solstice the chosen bandwidth lies very close to the optimum; the width of a nonparametric confidence interval is set by the effective sample size $nh$, and the bias, which is of the same order as the standard deviation, must also be confronted directly (by undersmoothing or bias correction), with the bootstrap commonly used for standard errors; the semiparametric partially linear model uses Robinson's double residuals to partial out the nuisance function nonparametrically, so that the $\beta$ of interest is estimated at the $\sqrt n$ rate and the curse of dimensionality is avoided.

## 5 Anchoring Papers

Methods prove themselves only in real research. This section anchors the chapter in three papers: one that founded the $\sqrt n$ semiparametric estimation of the partially linear model, one that established the boundary and efficiency advantages of local linear regression, and one that applied nonparametric demand estimation to a real policy evaluation. Each is laid out along the five elements of paper, method, data, results, and limitations.

### 5.1 Robinson (1988)

::: {.case}
Paper: "Root-N-Consistent Semiparametric Regression," Econometrica 56(4):931-954. It founded the semiparametric estimation of the partially linear model, and both the double-residual estimator of Section 4.3 of this chapter and the "$g$ slow, $\beta$ fast" conclusion come from it.

Method (using this chapter's notation, denoting the nuisance dimension by $x$ and the object of interest by $z$): Robinson studies the partially linear model $y = z'\beta + g(x) + \varepsilon$, where $\beta$ is the finite-dimensional coefficient to be estimated and $g$ is an unknown nonparametric function of $x$. His key step is to take expectations conditional on $x$ and subtract, which eliminates $g(x)$ and yields $y - \mathbb{E}[y\mid x] = (z - \mathbb{E}[z\mid x])'\beta + \varepsilon$, an equation in which $g$ no longer appears. The feasible double-residual estimator estimates the two conditional expectations by nonparametric kernel regression, forms the residuals, and then regresses the $y$ residuals on the $z$ residuals by ordinary least squares.

Data: the paper is methodological, not tied to a specific dataset, but the estimator it gives was thereafter used in countless empirical studies that carry a single flexible control term.

Results: the core conclusion is that although $g$ can be estimated only at the slow nonparametric rate, $\beta$ still converges at the parametric $\sqrt n$ rate and is asymptotically normal, with a consistently estimable variance (a trimming device handles the problem of the random denominator in the kernel estimate). This result, that the slow part does not drag down the fast part, is the cornerstone of the semiparametric partially linear model and a forerunner of the residual-orthogonalization logic of later double/debiased machine learning (Chapter 24).

Limitations: the partially linear model assumes that the variable of interest $z$ enters linearly, and if this assumption is wrong, $\beta$ is still biased; each of the two first-step nonparametric estimates requires its own bandwidth, which adds a layer of tuning; and trimming the random denominator involves practical choices of its own.
:::

This paper is exemplary because it demonstrated a deep and practical division of labor: parameterize the question you want to answer, nonparameterize the nuisance that gets in the way, and then "partial out" to separate the two cleanly. The double residuals of Section 4.3 of this chapter are its core. On Solstice, Robinson's estimator pulls the advertising coefficient, biased under linear control, from 0.752 back to 0.789 (repeated-sampling mean 0.800, unbiased), which is exactly the payoff of this logic.

### 5.2 Fan (1992)

::: {.case}
Paper: "Design-Adaptive Nonparametric Regression," Journal of the American Statistical Association 87(420):998-1004. It analyzed local linear regression and proved its design adaptivity and boundary advantage; the boundary bias correction of Section 3.3 of this chapter draws on it.

Method: Fan analyzes nonparametric regression based on weighted local linear fitting, that is, fitting a line at each target point by kernel-weighted least squares and taking the intercept as the estimate, and compares it with the classical Nadaraya-Watson (local constant) kernel estimator.

Data: the paper is methodological, mainly theoretical analysis and simulation.

Results: the core property is design adaptation: the local linear estimator can automatically adapt to random and fixed designs, clustered and uniform designs, and above all can adapt simultaneously to interior points and boundary points. Nadaraya-Watson has large bias at boundaries and where the data are uneven, while local linear absorbs these boundary and design effects into its slope term and corrects them automatically, needing no special boundary treatment. Fan further proves it has good efficiency properties, close to minimax optimal over standard smoothness classes. This paper, together with Fan and Gijbels (1996), made local polynomial smoothing displace Nadaraya-Watson as the modern default choice.

Limitations: local linear has larger variance where data are sparse; it too must choose a bandwidth; going up to local quadratic or higher order can further reduce bias, but the cost in variance and stability rises accordingly.
:::

Fan's contribution elevated a seemingly minor change (local constant replaced by local linear) into default practice. The boundary bias table of Section 6.3 of this chapter is its embodiment: Nadaraya-Watson has an obvious systematic bias at both ends, local linear is nearly unbiased, and the difference comes precisely from the slope that local fitting carries of its own.

### 5.3 Hausman and Newey (1995)

::: {.case}
Paper: "Nonparametric Estimation of Exact Consumers Surplus and Deadweight Loss," Econometrica 63(6):1445-1476. It applied nonparametric demand estimation to a real welfare evaluation, and it shares one kernel with Solstice: do not presuppose a functional form for the demand curve, let the data speak the shape, because the shape directly determines the conclusion.

Method: Hausman and Newey use nonparametric methods (kernel estimation and series/spline estimation together) to estimate the gasoline demand curve, and from it compute the exact consumer surplus and deadweight loss of a gasoline tax, avoiding the earlier practice of driving welfare estimates by functional-form assumptions. They also develop asymptotic normality theory for kernel and series estimators and for functionals such as surplus and deadweight loss, and test theoretical restrictions such as the downward slope and symmetry (Slutsky) of compensated demand.

Data: US household-level cross-sectional data, from the Department of Energy's household surveys of 1979 to 1981, containing gasoline demand, price, and income.

Results: the estimated demand curve is nonlinear in price, and this nonlinearity substantively changes the welfare calculation; the average deadweight loss of a gasoline tax is a considerable fraction of tax revenue, and the nonparametric estimate differs markedly from the conclusion implied by a traditional parametric (such as log-linear) specification. This is a benchmark illustration of "flexible nonparametric demand estimation changes substantive economic conclusions."

Limitations: cross-sectional identification requires price variation to be exogenous; the nonparametric welfare functionals are data-hungry; and the curse of dimensionality limits the number of covariates that can be included, since adding one more continuous covariate significantly worsens the estimate.
:::

The significance of this paper for this chapter is that it turned "the shape is the answer" from a slogan into an empirical study with policy consequences: whether the demand curve bends directly determines how large the deadweight loss of a tax comes out. It demonstrated where nonparametric methods earn their keep in real economic problems, and honestly exposed their exacting demands on sample size and dimension.

Taking the three together, the point of anchoring is clear: Robinson founded the $\sqrt n$ semiparametric estimation of the partially linear model, Fan clarified the advantages of local linear regression and fixed the boundary, and Hausman and Newey applied nonparametric demand estimation to welfare evaluation. They share this chapter's core idea: when the shape itself is the answer, only by not presupposing a functional form, and letting the data choose the shape, do you avoid smoothing away the question in advance; semiparametrics, for its part, strikes a pragmatic balance between flexibility and precision.

## 6 A Full Walkthrough on Solstice Data

Now we run the earlier tools on Solstice from start to finish. The code below uses R 4.6.1. The random seed of the data-generating process is fixed at 848, and each Monte Carlo experiment sets its own seed (9, 31, 55 and 72, in scripts 03 to 06 of `code/14_nonparametrics/`); every number quoted in the text comes from the actual output of this code.

### 6.1 The Data-Generating Process

We design two blocks of data with one seed. In the demand block, log demand is $y = g(x) + \varepsilon$, price $x$ is uniform on $[1, 9]$, the true curve $g(x) = 4 - 0.18x - 1.2\,\Lambda(2(x - 5))$ is a gentle overall decline with a steep drop near price 5 superimposed ($\Lambda$ is the logistic function), and the noise standard deviation is 0.35. In the partially linear block, $y = \beta z + g(x) + \varepsilon$, advertising exposure $z$ is nonlinearly related to price, and $\beta = 0.8$.

```r
set.seed(848); n <- 600
g <- function(x) 4 - 0.18*x - 1.2*plogis(2*(x - 5))  # true demand curve
x <- runif(n, 1, 9)
y <- g(x) + rnorm(n, 0, 0.35)                        # demand block
z <- 0.15*x + 1.1*plogis(1.8*(x - 5)) + rnorm(n, 0, 0.7)  # advertising, nonlinearly related to price
yp <- 0.8*z + g(x) + rnorm(n, 0, 0.35)               # partially linear block
```

The local slopes of the true curve at prices 2 / 5.5 / 8 are $-0.186 / -0.652 / -0.186$, and the corresponding local elasticities $x\,g'(x)$ are $-0.372 / -3.585 / -1.487$, a nearly tenfold range. In the partially linear block, the true value is $\beta = 0.8$. These numbers are the targets of all the estimates that follow.

### 6.2 A Straight Line Cannot Serve as a Demand Curve

First look at that lying constant-elasticity fit.

```r
cel <- lm(y ~ log(x))                                # constant elasticity: log-demand on log-price
lin <- lm(y ~ x)                                     # linear
```

The constant-elasticity model reports a single elasticity $-1.584$, and the linear model reports a single slope $-0.404$. But the true elasticity at price 2 is $-0.372$ (inelastic, should raise the price), and at 5.5 is $-3.585$ (elastic, cannot raise the price), so a single number simply cannot cover it. The $R^2$ of the linear fit is 0.859, which looks respectable, but the residuals carry a systematic curvature, and a RESET-type test returns a $p$-value of $2 \times 10^{-20}$, rejecting "linearity is enough." The lying elasticity of Section 1 is confirmed here.

### 6.3 Kernel Regression and Local Linear

Using a Gaussian kernel to do Nadaraya-Watson and local linear, with the bandwidth first taken as 0.6.

```r
Kg <- function(u) dnorm(u)
nw <- function(x0, h) { w <- Kg((x0 - x)/h); sum(w*y)/sum(w) }        # local constant
ll <- function(x0, h) { w <- Kg((x0 - x)/h); X <- cbind(1, x - x0)   # local linear
  solve(t(X) %*% (X*w), t(X) %*% (w*y)) }                            # [1]=m_hat, [2]=m'_hat
```

Both can trace out the whole bent demand curve, and local linear also gives the slope on the side: the slopes it estimates at prices 2 / 5.5 / 8 are $-0.249 / -0.620 / -0.198$, against the true values $-0.186 / -0.652 / -0.186$, and it captures especially that steep $-0.62$ in the sensitive region. The most important difference between kernel regression and local linear is at the boundary, and boundary bias is a systematic (repeated-sampling) property, so a simulation is needed to see it clearly. Resample the DGP 400 times and look at the average bias of the two estimators at the left end, interior, and right end:

| Location | Left end $x = 1.2$ | Interior $x = 5.0$ | Right end $x = 8.8$ |
|---|---|---|---|
| Nadaraya-Watson bias | $-0.0657$ | $-0.0026$ | $+0.0637$ |
| Local linear bias | $+0.0020$ | $-0.0016$ | $-0.0027$ |

Kernel regression has an obvious systematic bias at both boundaries (low at the left end, high at the right end, both dragged inward), and is nearly unbiased in the interior; local linear is nearly zero-bias at all three places. This confirms Section 3.3: kernel regression at the boundary "only takes averages" and is dragged off by the inner-side neighbors, while local linear "reads the trend" and extrapolates along the slope. On a single dataset this bias is covered by the variance and cannot be seen clearly, but repeated sampling brings it out into the open, which also restates that bias is a concept about repeated sampling.

![Left: the true demand curve (red), the linear fit (gold dashed), and the nonparametric fit (navy, GCV smoothing spline). The straight line is too steep in both flat zones and too shallow through the steep drop in the middle, leaving S-shaped residuals; the nonparametric fit almost coincides with the truth, recovering how elasticity varies with price. Right: the bias-variance tradeoff. As the bandwidth $h$ increases, the squared bias (red) rises and the variance (blue) falls, and the mean squared error (navy) is U-shaped, with the cross-validation-chosen bandwidth 0.37 (dotted line) sitting right at the bottom of the U.](assets/fig/fig_14_curve.svg)

### 6.4 Bias-Variance and Cross-Validation

Sweep the bandwidth from very small to very large, and use repeated sampling to split the mean squared error into squared bias and variance:

| Bandwidth $h$ | 0.15 | 0.25 | 0.40 | 0.60 | 0.90 | 1.30 | 2.00 |
|---|---|---|---|---|---|---|---|
| Bias$^2$ | 0.00001 | 0.00005 | 0.00035 | 0.00133 | 0.00481 | 0.01104 | 0.01831 |
| Variance | 0.00303 | 0.00196 | 0.00122 | 0.00080 | 0.00059 | 0.00045 | 0.00048 |
| MSE | 0.00304 | 0.00202 | 0.00157 | 0.00213 | 0.00540 | 0.01148 | 0.01879 |

The squared bias rises monotonically with the bandwidth, the variance falls basically monotonically (rising slightly at very large bandwidth), and the mean squared error bottoms out at $h = 0.40$, a textbook U shape. The key is whether cross-validation can automatically find this bottom: the bandwidth chosen by least-squares cross-validation is 0.372, very close to the 0.40 that minimizes the true MSE. The data tuned the aperture to be neither too large nor too small by "predicting held-out points." The spline route tells the same story: the smoothing spline chosen by generalized cross-validation (GCV) has an effective degrees of freedom of 8.3, at which point the MSE against the truth is 0.00109, whereas deliberately taking 3 degrees of freedom (oversmoothing) worsens the MSE to 0.01859 and taking 30 degrees of freedom (overfitting) gives 0.00422, both ends worse, with the optimum right in the middle chosen by GCV.

### 6.5 Semiparametrics: Robinson Recovers the Advertising Effect

In the partially linear block, we estimate the coefficient $\beta$ on advertising exposure $z$ (true value 0.8); $z$ is nonlinearly related to price.

```r
b1 <- coef(lm(yp ~ z))["z"]                          # only z, dropping price
b2 <- coef(lm(yp ~ z + x))["z"]                      # linear control for price
# Robinson double residuals: partial price nonparametrically out of y and z
ey <- yp - predict(smooth.spline(x, yp), x)$y
ez <- z  - predict(smooth.spline(x, z),  x)$y
b_rob <- coef(lm(ey ~ ez - 1))["ez"]
```

Four approaches give four numbers: only $z$ dropping price gives 0.159 (heavily biased by the confounding of price), linear control for price gives 0.752 (still biased, because the price effect is nonlinear and a linear term does not control cleanly), cubic-polynomial control for price gives 0.780 (better), and Robinson double residuals give 0.789 (standard error 0.021, hugging the true value). Repeated sampling makes it clearer: the estimate with linear control for price averages 0.758 (systematically low), the Robinson estimate averages 0.800 (unbiased), and Robinson's sampling standard deviation is only 0.019, a $\sqrt n$-level precision. Strip price out nonparametrically, and the advertising effect is clean, and without sacrificing precision for it.

![Left: four estimates of the advertising-exposure coefficient $\beta$ in the partially linear model (true value 0.80). Using only z (dropping price) gives 0.16; a linear control for price gives 0.75, still biased; a cubic polynomial gives 0.78; and Robinson's double residuals give 0.79, close to the truth (repeated-sampling mean 0.800). Only by controlling for price flexibly can the advertising effect be estimated cleanly. Right: the curse of dimensionality. To enclose 10% of the data in $d$ dimensions, the fraction each edge of the neighborhood must span climbs from 0.10 in one dimension all the way to 0.79 in ten dimensions, and "local" exists in name only in high dimensions.](assets/fig/fig_14_semipar.svg)

### 6.6 The Curse of Dimensionality

Finally, we measure the curse of dimensionality. To enclose $10\%$ of the data in a $d$-dimensional unit cube, the fraction each edge of the neighborhood must span, $0.1^{1/d}$, is: 0.100 in one dimension, 0.316 in two, 0.631 in five, 0.794 in ten. In ten dimensions you must span nearly eight-tenths on every coordinate axis to hold a tenth of the data, and "local" is entirely gone. In estimation, with the sample size fixed at $n = 2000$ and kernel regression used to estimate an additive true function at the origin, the mean squared error rises from 0.00015 in one dimension to 0.00067 in five (about 4.6 times as large), worsening steadily with dimension. This is why, with several regressors, one turns to semiparametric methods, using parametric structure to hold down the dimensions that do not need flexibility.

### 6.7 The Full Reconciliation

Putting the key numbers side by side:

| Quantity | Method | Result | Target |
|---|---|---|---|
| Price elasticity | Constant elasticity (single number) | $-1.584$ | true value ranges with price from $-0.37$ to $-3.59$ |
| Local slope | Local linear @ 2/5.5/8 | $-0.249 / -0.620 / -0.198$ | true value $-0.186 / -0.652 / -0.186$ |
| Boundary bias | NW / local linear @ left end | $-0.066 / +0.002$ | local linear unbiased |
| Bandwidth | MSE-optimal / CV | $0.40 / 0.372$ | CV hits it |
| Advertising coefficient $\beta$ | Linear price control / Robinson | $0.752 / 0.789$ | true value 0.80 |
| Curse of dimensionality | edge length $e_d(0.1)$ @ d=1/5/10 | $0.10 / 0.63 / 0.79$ | local fails |

The reconciliation of this section can be summarized as follows: constant elasticity compresses a true elasticity that varies tenfold into a single $-1.584$, which misleads pricing; kernel regression and local linear trace out the whole curved demand, with local linear nearly unbiased at the boundary and also giving the local slope; the bias-variance tradeoff is U-shaped in the bandwidth, the cross-validation choice of 0.372 lies very close to the optimum of 0.40, and the spline's GCV behaves likewise; Robinson's double residuals strip out price nonparametrically, pulling the advertising coefficient, biased under linear control, from 0.752 back to 0.789 (repeated-sampling mean 0.800, unbiased) while keeping $\sqrt n$ precision; the edge length in the curse-of-dimensionality calculation rises from 0.10 in one dimension to 0.79 in ten, confirming that pure nonparametrics fails in high dimensions.

## 7 Failure Modes and Robustness

The curves in the simulation were built to order; in real data they hide deeper. This section lays out the most common ways things fail and how to respond.

The choice of smoothness is the first and also the most fatal judgment. A bandwidth or number of knots chosen too large will smooth away the true bends (oversmoothing, large bias), and you will get a deceptively smooth curve that mistakes the sensitive region of demand for the flat region; chosen too small, it will chase the noise (undersmoothing, large variance), and you will mistake random fluctuation for real structure, reading off a "price kink" that does not exist at all. Cross-validation is the default tool for choosing smoothness, but it is not omnipotent: it optimizes overall prediction error, which is not necessarily optimal on the stretch you care most about (say near the current price), and it is itself unstable where the data are sparse or the sample is small. The sound practice is to report sensitivity to smoothness, checking whether the conclusion holds up over a reasonable range of bandwidths.

Boundaries and sparse regions are where nonparametric estimation is least trustworthy. Kernel regression has systematic bias at boundaries (which local linear can ease), while in sparse price intervals (say the extreme prices at which few products are ever set), the variance of any nonparametric method explodes, and the curve there is basically guesswork. The discipline of judgment is: a nonparametric curve is trustworthy where the data are dense, doubtful at boundaries and in sparse regions, and when reporting you should draw the confidence band so the reader can see which stretch is solid and which is empty, rather than handing over a smooth curve that is equally trustworthy everywhere.

For inference, nonparametrics has a trap of its own: the bias, which is of the same order as the standard deviation. If you simply copy the parametric approach and use a symmetric interval centered at the point estimate, coverage will fall short of the nominal level, because the estimate itself is systematically off by about half a standard deviation. The response is undersmoothing (tuning the bandwidth smaller than optimal to suppress the bias) or explicit bias correction, with standard errors commonly obtained by the bootstrap. Treating the nonparametric point estimate as unbiased and directly applying a parametric-style confidence interval is a common and insidious error.

Dimension is the sword hanging over all nonparametric methods. Beyond three or four continuous regressors, pure nonparametrics is basically unusable, and forcing it will only yield a surface with enormous variance that cannot be trusted anywhere. At that point you should not cling to nonparametrics but return to semiparametrics: think through which dimensions truly need flexibility and which can be held down by parametric structure (additive, partially linear, single-index), and keep the nonparametric treatment only for the former. Semiparametrics has its own costs too: the partially linear model assumes that the variable of interest enters linearly, and if this assumption is wrong, $\beta$ is still biased; each of Robinson's two first-step nonparametric estimates also needs its own bandwidth, which adds a layer of tuning.

Stringing these failure modes together, the robustness of nonparametrics and semiparametrics ultimately rests on two judgments: one is whether the smoothness is tuned right and whether the conclusion holds up over a reasonable range, which decides whether the curve you trace is the true shape or the false image of oversmoothing or the noise of overfitting; the other is whether the dimension can be sustained, and which parts should be nonparametric and which held down by structure, which decides whether you are trapped to death by the curse or slip cleverly around it with semiparametrics. Think these two judgments through, and the freedom of "letting the data choose the shape" truly becomes yours to use; fail to think them through, and nonparametrics may either uncover for you the truth smoothed away by a parametric form, or hand you a curve with enormous variance, boundary guesswork, and lying intervals.

## 8 Further Reading

::: {.readings}
Required reading, in suggested reading order:

- Hastie, Tibshirani and Friedman (2009, *The Elements of Statistical Learning*, 2nd ed.), Chapters 5 to 6. The most readable modern textbook treatment of kernel methods, local regression, splines, and the bias-variance tradeoff, and much of the intuition for kNN, kernels, and splines in this chapter is drawn from it.
- Fan (1992, *JASA*). The key source for local linear regression and its boundary advantage, and the source of the boundary bias correction in Section 3.3 of this chapter.
- Robinson (1988, *Econometrica*). The founding work on $\sqrt n$ semiparametric estimation of the partially linear model, and the original source of the double-residual estimator of Section 4.3 and the "$g$ slow, $\beta$ fast" conclusion.

Further reading:

- Stone (1977, *Annals of Statistics*) and Stone (1982, *Annals of Statistics*). The theoretical bedrock of consistency and optimal convergence rates for nonparametric regression, and that dimension-dependent rate $n^{-4/(4+d)}$ comes from the latter.
- Silverman (1986, *Density Estimation for Statistics and Data Analysis*). The classic monograph on kernel density estimation and bandwidth selection, and that daunting dimension-sample-size table comes from it.
- Nadaraya (1964) and Watson (1964). The two independent original sources of the kernel regression estimator, and reading them shows clearly where the "weighted average" step comes from.
- Newey (1994, *Econometrica*). The general theory of the asymptotic variance of semiparametric estimators, useful for understanding why the first-step nonparametric error does not drag down the second-step parametric convergence.
- Li and Racine (2007, *Nonparametric Econometrics: Theory and Practice*). The authoritative monograph on nonparametric and semiparametric econometrics, with complete coverage from kernels to splines to partially linear models; a desk reference to consult before doing the work yourself.
:::

::: {.apa-refs}
- Fan, J. (1992). Design-adaptive nonparametric regression. *Journal of the American Statistical Association, 87*(420), 998-1004. https://doi.org/10.1080/01621459.1992.10476255
- Fan, J., & Gijbels, I. (1996). *Local polynomial modelling and its applications*. Chapman & Hall.
- Hastie, T., Tibshirani, R., & Friedman, J. (2009). *The elements of statistical learning: Data mining, inference, and prediction* (2nd ed.). Springer.
- Hausman, J. A., & Newey, W. K. (1995). Nonparametric estimation of exact consumers surplus and deadweight loss. *Econometrica, 63*(6), 1445-1476. https://doi.org/10.2307/2171777
- Li, Q., & Racine, J. S. (2007). *Nonparametric econometrics: Theory and practice*. Princeton University Press.
- Nadaraya, E. A. (1964). On estimating regression. *Theory of Probability and Its Applications, 9*(1), 141-142. https://doi.org/10.1137/1109020
- Newey, W. K. (1994). The asymptotic variance of semiparametric estimators. *Econometrica, 62*(6), 1349-1382. https://doi.org/10.2307/2951752
- Robinson, P. M. (1988). Root-N-consistent semiparametric regression. *Econometrica, 56*(4), 931-954. https://doi.org/10.2307/1912705
- Silverman, B. W. (1986). *Density estimation for statistics and data analysis*. Chapman & Hall.
- Stone, C. J. (1977). Consistent nonparametric regression. *The Annals of Statistics, 5*(4), 595-620. https://doi.org/10.1214/aos/1176343886
- Stone, C. J. (1982). Optimal global rates of convergence for nonparametric regression. *The Annals of Statistics, 10*(4), 1040-1053. https://doi.org/10.1214/aos/1176345969
- Watson, G. S. (1964). Smooth regression analysis. *Sankhyā: The Indian Journal of Statistics, Series A, 26*(4), 359-372.
:::
