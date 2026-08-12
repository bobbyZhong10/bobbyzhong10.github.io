---
title: "Demand Estimation: Berry Inversion & BLP"
subtitle: "From Market Shares to Elasticities and Substitution"
seriesline: "Foundations of Information Systems Economics · Chapter 20"
series: "Foundations of Information Systems Economics: Empirical Methods for Information Systems Economics"
sidebartitle: "Foundations of Information Systems Economics"
sidebarsub: "Theory · Methods · Applications"
chapterlabel: "Chapter 20 · Demand Estimation: Berry Inversion & BLP"
---

## Introduction

A spreadsheet from the car market will tell you how many units of each model sold in each city, at what price, with how much horsepower, but it usually will not tell you where a consumer who did not buy a given model turned instead. Yet that substitution is exactly what merger review, pricing, and welfare analysis care about most. Once the individual choice records vanish, demand does not vanish with them; the researcher simply has to reassemble the underlying consumer preferences out of a single table of market shares.

Berry's key insight is to treat market shares as an invertible map of mean utilities. Observe the shares and you can back out the $\delta_{jt}$ that makes the model generate exactly those shares; then use instruments to separate price from unobserved quality. In plain logit this inversion is just a one-line difference of log shares. Add random coefficients and the closed form is gone, but you can still find that vector of mean utilities with a contraction mapping. A problem that looked like a high-dimensional integral equation is thereby organized into two pieces of work: an inner inversion and an outer GMM.

Computational success does not automatically buy economic credibility. Plain logit inherits IIA; if your price instrument is merely correlated with price in the sample while also reflecting product quality, no amount of second-stage polish will save you; and random coefficients, without enough variation in shares and valid instruments, can quietly collapse back to simple logit. Selva's six meal-kit brands and three hundred markets run through the whole chapter. Because the truth in this simulated world is known, we can watch exactly how the wrong substitution structure propagates from the price coefficient through elasticities, markups, and marginal costs.

This chapter builds on the random utility, IIA, and mixed logit of Chapter 19, and it uses GMM and instrumental variables. The goal is not to treat BLP as a black-box procedure but to leave you able to write down the Berry inversion, explain which variation identifies the price coefficient versus the random coefficients, check whether the contraction and the outer optimization actually converged, and judge the economic content of an instrument. The next chapter picks up the estimated demand system and backs out costs and market power from firm pricing.

## 1 A number that doesn't add up

::: {.case}
Selva is a fictional delivery platform, and in this chapter we look at one of its meal-kit categories. The data are at the market level: the platform operates in 300 city-week markets, each market offers six brands, labeled 1 through 6, plus an outside option of not ordering. For each brand in each market we observe only three things: the market share, the price, and two observable characteristics ($x_1$ records quality, $x_2$ records menu breadth). We see no individual, who they are or what they bought, only the aggregated shares. What the platform wants to know is the same old question: how sensitive are users to price, and if a given brand raises its price, how much demand does it lose and to whom.
:::

Start with Berry's most direct step. Market shares are a function of utility, and in plain logit that function can be inverted in one line: the log share of product $j$ minus the log share of the outside option is exactly its mean utility $\delta_{jt} = x_{jt}\beta - \alpha p_{jt} + \xi_{jt}$. So take $\ln s_{jt} - \ln s_{0t}$ as the dependent variable, run an ordinary least squares regression on characteristics and price, and out pops the price coefficient $\alpha$.

The result is absurd. OLS estimates a price coefficient of $-0.19$, and note the sign: utility is written as $-\alpha p$, so $\alpha$ should be positive to make price a disutility, yet it comes out negative, meaning price enters utility positively. Translated into demand language, this demand curve slopes upward, with an own-price elasticity of $+1.2$: the more Selva's meal kits cost, the more users buy. No sensible demand curve has this shape. Take it seriously and you can derive even more ridiculous conclusions. Back out marginal cost from the pricing first-order condition, and an upward-sloping demand gives negative markups, implying a median marginal cost of $13.6$ dollars per meal while these brands sell at a median of only $7.6$ dollars. That is to say every brand, in every market, on every order, is selling at a loss of five or six dollars. This is of course impossible.

The regression is not miswritten and the inversion is not miswritten; Berry's line of algebra is an exact identity. What is wrong is treating price as exogenous. Expensive brands tend to be better along dimensions the researcher cannot see, so $\xi_{jt}$ raises both a brand's share and its price. OLS cannot tell "expensive" from "good," swallows the entire disutility of price and then some, and forces out an upward-sloping demand curve. This is exactly the price endogeneity of the previous chapter, only in market-level data its consequences can be extreme enough to flip the sign of the demand slope. Section 4 will show that once a cost instrument tames price, $\alpha$ returns to a positive $0.51$ and demand slopes down again; but that is not enough, because plain logit's IIA makes elasticities systematically too small (own-price elasticity is only $-3.3$, against a truth of about $-4.4$) and distorts the substitution pattern too. To get both elasticities and substitution right you need full BLP, at which point $\alpha$ lands at $0.75$, right up against the truth of $0.70$.

## 2 Economic model and estimand

Before doing any inversion, let us write down the market-level demand model and the target quantity clearly. It is a direct aggregation of the previous chapter's individual model, and the notation is almost carried over verbatim.

What demand estimation must deliver is still those two things: own-price elasticities and the substitution pattern (the diversion ratio). What is special about market-level data is that neither can be read directly off the shares; you must first estimate the utility parameters and then plug them into the elasticity formulas. So the estimand of this chapter is the utility parameters $(\beta, \alpha)$ and the random-coefficient distribution parameters $\theta_2$, along with the elasticity matrix and diversion ratios derived from them. Only once these are in hand can the supply side be brought in to recover marginal costs and evaluate mergers, which is the next chapter's business.

The notation follows the previous chapter. Index markets by $t$, products by $j \in \{1, \dots, J\}$, with $j = 0$ the outside option. The utility consumer $i$ gets from product $j$ is

$$u_{ijt} = \underbrace{x_{jt}\beta - \alpha p_{jt} + \xi_{jt}}_{\delta_{jt}} + \underbrace{\mu_{ijt}}_{\text{random coef.}} + \varepsilon_{ijt}, \qquad u_{i0t} = \varepsilon_{i0t}$$

$\delta_{jt}$ is the mean utility of product $j$ in market $t$, common to all consumers, holding the observable characteristics $x_{jt}$, price $p_{jt}$, and the quality $\xi_{jt}$ that the researcher cannot observe but the platform and consumers both know. $\mu_{ijt}$ is the consumer-varying part of preferences, the random coefficients from the previous chapter: $\mu_{ijt} = -\sigma_p \nu_{ip} p_{jt} + \sum_k \sigma_k \nu_{ik} x_{jt}^{(k)}$, where $\nu_i$ is the individual's taste draw and $\sigma_p$ and $\sigma_k$ describe how dispersed tastes are over price and over the $k$-th characteristic respectively (price counts as one characteristic carrying a random coefficient). $\varepsilon_{ijt}$ is Type-I extreme value. Turn $\mu$ off (all $\sigma_k = 0$) and you fall back to plain logit.

Given $\delta_t$ and the nonlinear parameters $\theta_2 = \sigma$ (the vector of random-coefficient standard deviations), the predicted share of product $j$ in market $t$ is an integral over the taste distribution. Given $\delta_t$, the share depends on $\sigma$ only through $\mu$, whereas $\alpha$ and $\beta$ are linear parameters that enter through $\delta$ and are solved by concentrated 2SLS (see Section 4.3):

$$s_{jt}(\delta_t, \theta_2) = \int \frac{\exp[\delta_{jt} + \mu_{ijt}]}{1 + \sum_k \exp[\delta_{kt} + \mu_{ikt}]}\, dF(\nu)$$

This equation is the pivot of the whole chapter. It says the predicted share is the average, over the taste distribution, of each consumer's internal logit probability, exactly the mixed logit share from the previous chapter, except that now we do not observe individuals and only observe the realized value of this integral (the market share). The first step of estimation is to invert this map from $\delta$ to shares: given the observed shares, find the vector of $\delta_{jt}$ that makes predicted shares exactly equal observed shares.

::: {.intuition}
Why insist on first inverting out $\delta_{jt}$ rather than running a nonlinear regression on shares directly? Because $\xi_{jt}$ hides inside $\delta_{jt}$, and it is correlated with price, the entire source of endogeneity. If you can solve $\delta_{jt}$ out on its own, then $\delta_{jt} = x_{jt}\beta - \alpha p_{jt} + \xi_{jt}$ is a linear equation with $\xi_{jt}$ as its residual, and the whole arsenal of instrumental-variable tools we have for endogenous residuals becomes immediately available. The point of the inversion is to reorganize a nonlinear integral equation, in which $\xi$ is hopelessly entangled, into a linear moment condition where $\xi$ appears cleanly in the residual position. All of BLP's ingenuity is in the service of completing this inversion and making it run numerically.
:::

To summarize: the market-level demand model is the aggregation of the previous chapter's individual model; $\delta_{jt}$ holds the endogenous $\xi_{jt}$; the predicted share is the integral of a random-coefficient logit; the estimand is the utility and random-coefficient parameters plus the elasticities and diversion ratios derived from them; and the inversion reorganizes endogeneity into an instrumental-variable problem in the residual, which is the starting point of the next section's identification.

## 3 Identification

Identification comes before estimation in logic. This section answers: starting from observables like market shares, prices, and characteristics, under what assumptions can the utility parameters and the substitution pattern be pinned down. This is the most finely worked section in the chapter, split into numbered subsections, each giving intuition before formalization. Inversion, endogeneity, instruments, and the identification of random coefficients are dissected one by one.

### 3.1 Berry inversion: recovering mean utility from shares

Start with the most transparent plain-logit inversion, the starting point for understanding everything, and it is just algebra. In plain logit there are no random coefficients, $\mu_{ijt} = 0$, and the share degenerates to

$$s_{jt} = \frac{\exp(\delta_{jt})}{1 + \sum_k \exp(\delta_{kt})}, \qquad s_{0t} = \frac{1}{1 + \sum_k \exp(\delta_{kt})}$$

Divide the two and take logs, and the whole denominator cancels completely:

$$\ln s_{jt} - \ln s_{0t} = \delta_{jt} = x_{jt}\beta - \alpha p_{jt} + \xi_{jt}$$

This is the Berry (1994) inversion. The left side is all data (the observed shares), and the right side is the utility we want to estimate plus the residual $\xi_{jt}$. It turns an entire discrete-choice model into a linear regression whose dependent variable is a difference of log shares. A few points to make clear. First, this inversion is one-to-one: each observed share corresponds to a unique $\delta_{jt}$, one $\xi_{jt}$ per product, no more and no less. Second, using observed shares versus true shares on the left agrees as the sample grows infinite, which is why the inversion is stable when market-level data are plentiful. Third, and this is the crux, $p_{jt}$ on the right is correlated with the residual $\xi_{jt}$, so running OLS directly is the disaster of Section 1, and you must use IV.

The Berry inversion can be written in a unified abstract form, $\delta_t = \sigma^{-1}(s_t; \theta_2)$, read as "invert the share map to recover mean utility." Plain logit is the special case where this inverse has a closed form; the nested-logit inverse carries one extra term, the within-nest share (Section 4.2); the random-coefficient BLP inverse has no closed form and must be solved by iteration (Section 3.4). Wherever a share appears inside the inverse, it corresponds to an endogenous variable and calls for an instrument, a rule that holds throughout.

### 3.2 Why price is endogenous: correlation from the supply side

Section 1 said price is correlated with $\xi_{jt}$; here we explain the mechanism fully, because it is both the root of endogeneity and the guide to finding instruments. $\xi_{jt}$ is product quality outside the researcher's view. The platform and merchants see it when they price, better things get priced higher, so price naturally carries a piece of $\xi$. More precisely, price is set by the supply side's optimal pricing, and optimal pricing is a function of demand, and demand contains $\xi$, so the equilibrium price is a function of $\xi$. This is not the incidental correlation of measurement error but a systematic correlation forced by an economic mechanism. In Selva's data the correlation between price and $\xi$ is $0.23$, which looks small but is enough to flip the sign of the OLS demand slope.

This mechanism also points to where the instrument should be found. Since price is endogenous because it moves with $\xi$ (the demand side), a variable that moves only with cost (the supply side) and does not enter utility can supply clean price variation. Such cost shifters are the most reliable instruments in BLP. In Selva it is the local input cost $w_{jt}$ in each market (ingredients, delivery, labor), which pushes price up but does not directly affect users' preference for a given meal, satisfying the exclusion restriction. Section 3.3 works through the available instruments systematically.

::: {.warning}
The direction of price endogeneity in market-level data can be extreme enough to flip the sign of the demand slope, which is unlike individual data. In the previous chapter's individual logit, endogeneity merely pushes $\alpha$ toward zero (elasticities too small); here OLS gives $\alpha < 0$ outright, an upward-sloping demand. The reason is that in market-level data price is almost entirely determined endogenously by the supply side, and the share of price coming from $\xi$ through pricing is heavier. Do not underestimate its destructive power in market-level data just because "in individual data endogeneity is only mild attenuation."
:::

### 3.3 Instruments and the exclusion restriction

After the inversion, $\delta_{jt} = x_{jt}\beta - \alpha p_{jt} + \xi_{jt}$ is a linear IV problem, $p_{jt}$ is endogenous, and it needs an instrument satisfying $\mathbb{E}[\xi_{jt} \mid z_{jt}] = 0$. There are several usable sources, each with its own logic and its own weakness.

Cost shifters are the cleanest class. Any variable that enters cost but not utility qualifies: local factor prices, exchange rates, raw-material costs, taxes. It moves price directly while being unrelated to $\xi$, the tightest logic there is. The weakness is that good cost shifters are not always found, especially when a digital product's marginal cost is nearly zero.

The BLP instrument is the original device of Berry, Levinsohn, and Pakes, using rivals' product characteristics. The intuition is that the more competitors a product has around it, and the more similar they are, the greater the substitution pressure it faces, the lower the equilibrium markup, and the lower the price, so rival characteristics move the product's own price through the markup channel while not entering the product's own utility. The common construction is the sum of rival characteristics in the same market. The weakness is the one Armstrong (2016) points out: when the number of products per market is large, the markup tends to a constant, and this class of instrument loses its power to move price, degenerating into a weak instrument, at which point you still cannot do without a cost shifter.

The Hausman instrument uses the same product's price in other markets. The logic is that prices in different markets share cost shocks but not the demand shocks specific to each market. The weakness is that this "demand shocks are not correlated across markets" assumption is strong; national advertising or national taste changes both break it.

Gandhi and Houde's (2019) differentiation instrument is designed specifically to identify the random coefficients. It uses the distance between a product and its rivals in characteristic space, for example the sum of squared differences between each rival's and the product's characteristics in the same market. What it captures is "how densely packed the products are around this product," and the substitution strength governed by the random coefficients depends precisely on this density, so it is especially powerful for identifying $\sigma$. In Selva's estimation, it is exactly by adding rival costs and the differentiation instrument that the random coefficients get identified out of the trap of collapsing to zero.

Matching these instruments to what they identify is the key to understanding BLP identification. The mean price coefficient $\alpha$ is identified by exogenous variation in price, with cost shifters as the workhorse. The random coefficient $\sigma$ is identified by "how the substitution pattern varies with product characteristics": raise a product's price, and watch whether demand flows mainly to products similar to it in characteristics or spreads out evenly; the more it concentrates on similar products, the larger is $\sigma$. The differentiation instrument and the BLP instrument supply exactly this "who sits next to whom in characteristic space" variation. Identifying $\sigma$ is naturally harder than identifying $\alpha$, a point Section 7 will return to from the angle of failure modes.

### 3.4 Inversion with random coefficients: the contraction mapping

The plain-logit inversion has a closed form because the share is a simple function of $\delta$. Once random coefficients enter, the share is that integral, $\delta$ is wrapped inside the integral sign, and there is no closed-form inverse. BLP's key technical contribution is to provide an iteration guaranteed to converge that solves $\delta$ out.

For a fixed $\theta_2$, the $\delta_t$ we seek satisfies the $J$ equations "predicted shares equal observed shares." BLP proves that the following iteration is a contraction mapping:

$$\delta_t^{(h+1)} = \delta_t^{(h)} + \ln S_t^{\text{obs}} - \ln s_t\big(\delta_t^{(h)}, \theta_2\big)$$

The intuition is plain: if the current $\delta$ predicts a product's share below the observed share, raise its $\delta$ a bit, and vice versa, with the adjustment equal to exactly the difference of the two log shares. BLP proves this map has a contraction modulus less than 1, so no matter where you start, repeated iteration converges to the unique $\delta_t(\theta_2)$. In practice the convergence tolerance is set very tight (on the order of $10^{-12}$), because errors in the inner inversion get amplified in the outer GMM, a point Section 4.4 makes formally. There are acceleration tricks like SQUAREM, and you can also rewrite the inversion as a convex optimization problem and solve it with Newton's method, but the contraction iteration is the most classic and most robust approach.

![Starting from any point, the mean utilities δ for each product converge monotonically along the red staircase to the unique fixed point δ* under repeated application of the Berry contraction map T; at δ* the model's predicted shares exactly equal the observed shares and the inversion is complete. Because the slope of T is everywhere below 1 (contraction modulus less than 1), convergence holds for any starting point.](assets/fig/fig_20_contraction.svg)

With the ability to solve $\delta(\theta_2)$ for every $\theta_2$, the random-coefficient inversion is complete, and $\delta_{jt}(\theta_2) = x_{jt}\beta - \alpha p_{jt} + \xi_{jt}$ becomes the familiar linear IV equation again, only now $\delta$ depends on the $\theta_2$ we want to estimate. How to nest the inner and outer loops together for estimation is the business of Section 4.

### 3.5 What random coefficients identify, and how

The random coefficient $\sigma$ is the entire increment of BLP over plain logit, and its identification deserves its own careful treatment, because it is the hardest and most prone to failure.

$\sigma$ measures how far the substitution pattern departs from IIA. As the previous chapter explained, plain logit's IIA forces substitution to be allocated in proportion to shares, regardless of "who is like whom." Random coefficients break this: if many consumers prefer high quality, then two high-quality products are favored by the same group of people, and when one raises its price, demand flows more to the other high-quality one rather than scattering by shares. The larger is $\sigma$, the stronger this characteristic-similarity-based concentrated substitution. So all the identifying information for $\sigma$ is in "how the substitution pattern varies with product characteristics": the data must let us see that substitution strength really does differ across markets with different characteristic distributions, or across products with different characteristics.

This brings two identification conditions to keep firmly in mind. First, product characteristics must have enough cross-market variation in the data, otherwise "how substitution varies with characteristics" cannot be observed and $\sigma$ is not identified. Second, you need instruments that capture density in characteristic space, that is the differentiation and BLP instruments above. The random coefficient on quality $\sigma_{x_1}$ can be identified with this variation, and Selva estimates it quite precisely, which is exactly what holds up the near-neighbor substitution we will see later. The random coefficient on price $\sigma_p$ is another matter, and is famously the hardest parameter to identify in BLP: price is itself endogenous, and the exogenous variation available to identify it is very limited. Section 6 will show that even with rival costs and the differentiation instrument added, $\sigma_p$ in Selva is barely pinned down, with a point estimate near zero and a standard error far larger than the point estimate. This is not a program bug but a faithful reflection of the difficulty of identifying the price random coefficient, and the honest thing is to report this uncertainty rather than pretend it was estimated precisely. Fortunately the substitution pattern is governed mainly by $\sigma_{x_1}$, so an imprecise $\sigma_p$ does not ruin the elasticity and substitution structure this chapter must deliver.

### 3.6 Nonparametric identification: Berry-Haile

All the identification above hangs on functional forms like "utility is linear, $\varepsilon$ is Type-I EV." A natural question is, setting these parametric assumptions aside, how much market-level data can identify on its own. Berry and Haile (2014) give the answer, worth knowing because it builds BLP's credibility on a firmer foundation than "we happened to pick the right distribution."

Their conclusion is that under a few nonparametric conditions, demand (and even marginal cost) can be nonparametrically identified from market-level data. There are three core conditions. Connected substitutes guarantees the demand map is invertible, so shares can be inverted to mean utility, which is the generalization of the Berry inversion from logit to arbitrary models. Exogeneity of the instruments, $\mathbb{E}[\xi_{jt} \mid z_t] = 0$, supplies the exogenous variation needed to identify price. Completeness is the nonparametric version of the rank condition, guaranteeing the information carried by the instruments is enough to pin down the demand function. The key count is that a market has $J$ endogenous shares and $J$ endogenous prices, and in general you need $2J$ instruments, $J$ of them from exogenous product characteristics (the BLP-instrument type) and another $J$ from variables that enter only price, like cost shifters or Hausman instruments. If you are willing to accept the single parametric assumption that utility is linear in price, then the price part needs only one instrument, greatly lightening the instrumental-variable burden. This explains why a good cost shifter is so critical in practice.

The main points of this section can be summarized as follows: the Berry inversion turns shares into mean utility, with a closed form for plain logit and a contraction-mapping iteration for random coefficients; price is endogenous because of the supply side's optimal pricing and needs an instrument that enters only cost; the mean price coefficient is identified by cost shifters and the random coefficient by how substitution varies with characteristics, which requires a differentiation instrument and is naturally harder; and Berry-Haile places all of this on the solid ground of nonparametric identification. None of this yet says how to compute it, which is the next section.

## 4 Estimation: from IV logit to BLP

This section proceeds by increasing complexity, each step giving birth to the next where the last one fails. Plain logit's IV regression fixes endogeneity but locks down substitution; the random-coefficient BLP fixes substitution; and NFXP and MPEC are two ways of computing it.

### 4.1 IV logit: Berry inversion plus linear instruments

With the inversion from Section 3.1, estimating plain logit is a linear IV regression: run two-stage least squares of $\ln s_{jt} - \ln s_{0t}$ on $x_{jt}$ and $p_{jt}$, instrumenting price with cost shifters and the BLP instrument. On Selva, the first-stage coefficient on the cost shifter $w$ is $0.696$, with a partial-F as high as $1520$, an extremely strong instrument. The second-stage price coefficient returns to $0.51$ (against OLS's $-0.19$), demand finally slopes down, and the median own-price elasticity is $-3.3$. The endogeneity layer is cured.

But IV logit is still plain logit, and IIA remains. Its elasticities are systematically too small: the true median is about $-4.4$ while IV logit gives only $-3.3$, because the homogeneous-logit elasticity formula does not match the true substitution with random coefficients. Worse is the substitution pattern: IIA forces diversion to be allocated by shares, and Section 6 will show IV logit clearly underestimates near-neighbor substitution. Markups are then overstated too, and the backed-out marginal costs come in too low. To get both elasticities and substitution right, you must relax IIA.

### 4.2 Nested logit: one step away

The most economical step toward relaxing IIA is nested logit, sorting products into nests and allowing stronger within-nest substitution. Its inversion has just one term more than plain logit:

$$\ln s_{jt} - \ln s_{0t} = x_{jt}\beta - \alpha p_{jt} + \rho \ln s_{j \mid g, t} + \xi_{jt}$$

The extra $\ln s_{j \mid g, t}$ is the product's share within its nest, and $\rho$ is the within-nest correlation parameter (corresponding to the previous chapter's $1 - \lambda$). This term brings a new trouble: the within-nest share $s_{j \mid g, t}$ is itself endogenous, being a function of the dependent variable, which amounts to using a part of $Y$ to explain $Y$. Forgetting to find it its own instrument will, through simultaneity bias, push $\rho$ toward 1, faking a strong within-nest correlation that does not exist. A valid instrument is usually the number of products in the nest or rival characteristics. Nested logit is the transition from plain logit to BLP: it is more flexible than plain logit, but the nest structure must be specified by the researcher in advance, and mis-drawing the nests undoes everything, and it also assumes price sensitivity is the same for everyone. What truly relaxes all of this together is random coefficients.

### 4.3 Full BLP: inner inversion nested in outer GMM

BLP brings in random coefficients, and estimation becomes a two-layer nesting, a structure called the nested fixed point (NFXP). The parameters split into two groups: the linear parameters $\theta_1 = (\beta, \alpha)$ and the nonlinear parameters $\theta_2 = \sigma$.

The inner loop, for a fixed $\theta_2$, does the contraction mapping of Section 3.4 and solves out $\delta(\theta_2)$ making predicted shares equal observed shares. With $\delta(\theta_2)$, the linear equation $\delta_{jt}(\theta_2) = x_{jt}\beta - \alpha p_{jt} + \xi_{jt}$ can be solved for $\theta_1$ and the residual $\xi(\theta_2)$ in one step by IV (concentrating out the linear parameters, so the outer loop only needs to search over $\theta_2$). The outer loop uses the moment condition $\mathbb{E}[\xi_{jt} \mid z_{jt}] = 0$ to construct a GMM objective:

$$\hat\theta_2 = \arg\min_{\theta_2}\ \big(Z' \xi(\theta_2)\big)' W \big(Z' \xi(\theta_2)\big)$$

Every time the optimizer proposes a $\theta_2$, the inner loop re-solves the inversion, recomputes $\xi$, and recomputes the objective, until the outer loop converges. That is the full BLP. On Selva it estimates the price coefficient at $0.751$ (standard error $0.075$), right up against the truth of $0.70$; the quality coefficient $\beta_{x_1} = 1.661$ (truth $1.6$); the median own-price elasticity $-4.7$, close to the truth $-4.4$; and the backed-out median marginal cost $5.94$ dollars, near the truth $5.81$, with none negative. The elasticities, markups, and costs, the quantities that really need to be delivered, are all estimated correctly by BLP.

One implementation detail is worth stressing. That integral share has no analytic form and must be approximated by simulation or numerical integration: draw $R$ sets of tastes $\nu$ and average the $R$ individual logit probabilities. Fixing the number of draws introduces simulation error, and in practice you make $R$ large enough, or use quasi-Monte Carlo (such as Halton sequences) to push down the simulation variance. The convergence tolerance must be set very tight, because small errors in the inner loop get amplified in the outer loop.

In production research, few people hand-write this nested loop, and instead use mature software, most mainstream being Python's pyblp, which packages up the contraction mapping, optimal instruments, the supply side, and micro moments, so a single configuration estimates it:

```python
import pyblp
problem = pyblp.Problem(product_formulations, product_data)
results = problem.solve(sigma=initial_sigma, optimization=pyblp.Optimization('bfgs'))
```

To explain every step fully, this chapter hand-writes the contraction mapping and GMM in R, and every number in the text comes from an actual run of this hand-written code; once you understand the mechanics, switching to pyblp in production will be more stable and faster.

### 4.4 NFXP and MPEC: one estimator, two algorithms

BLP's nested fixed point has a well-known pain point: if the inner contraction is not solved accurately enough, the error seeps into the outer GMM and biases $\hat\theta_2$, so the tolerance must be set extremely tight, at high computational cost. Dubé, Fox, and Su (2012) give another computational path, MPEC (mathematical programming with equilibrium constraints). It no longer embeds the inversion inside the objective function but treats "predicted shares equal observed shares" as a set of constraints, letting the optimizer search over $\theta$ and $\delta$ (and even $\xi$) simultaneously, with the constraints holding only at the optimum.

The point is that NFXP and MPEC deliver the same statistical estimator, only by different algorithms, converging to the same set of FOCs. The difference is entirely in computational convenience: NFXP concentrates $\delta$ out and searches only over the low-dimensional $\theta_2$, which is efficient when $\sigma$ is not high-dimensional but is sensitive to the inner tolerance; MPEC has many more parameters ($\delta$ all enter the optimization) and a large Hessian, but its constraints hold only at the optimum, so it is insensitive to the inner tolerance and its derivatives are easier to compute too. Which to choose is a computational trade-off and does not change the estimator itself. Separating "identification" from "algorithm" is a general habit of mind for understanding structural estimation: the same identification, the same estimator, can be computed by many different methods.

### 4.5 Zero shares and other practicalities

A common headache with market-level data is zero shares. If a product sells zero in a market, $\ln s_{jt}$ is undefined and the left side of the Berry inversion cannot be computed. Simply dropping these observations or stuffing in a small constant both introduce selection bias, because products with observed zero shares systematically correspond to lower $\xi$, so dropping them non-randomly deletes the sample and biases demand toward being inelastic. Gandhi, Lu, and Shi give a proper treatment: treat shares as multinomial sampling under a finite market size and use it to derive bounds on the true shares (moment inequalities), or use a Laplace-type constant-added transform to make the log defined and analytically correct the bias it introduces. Zero shares are especially common in the long-tail categories of digital platforms, so this problem cannot be ignored.

### 4.6 Micro moments: reconnecting individual information

If, besides market shares, you also have some individual-level information (a small consumer survey, second-choice data, the interaction of demographics with choice), you can bring them in as extra moments for BLP, that is micro-BLP. Its value is in identifying power: random coefficients, especially the parts interacted with demographics, are hard to identify precisely from market shares alone, whereas the covariation of "who bought what" in the individual information pins them down directly. Which moment identifies which parameter is a matter of care: observable heterogeneity ($\Pi$, the interaction of characteristics with demographics) is identified by the covariance of characteristics and demographics, and unobservable heterogeneity ($\Sigma$) is identified by the covariance of second choices. Petrin's (2002) minivan study and Conlon and Gortmaker's (2023) modern treatment are the representatives of this path. The control function discussed at the end of the previous chapter is another way to tame individual price endogeneity, complementary to the moment method here.

The route of this section is now complete: Berry inversion plus IV is the estimation of plain logit, curing endogeneity but locking down substitution; nested logit relaxes it halfway; full BLP uses contraction-mapping inversion plus GMM to fix substitution, with NFXP and MPEC as two ways to compute it; zero shares need special treatment; and with individual information, micro moments can greatly strengthen identification.

## 5 Anchor papers

A method stands only when it lands in real research. Three anchor papers: one exposes the inversion idea, one pushes it to the general case and founds BLP, and one gives an application that is still a model of practice today. Each is laid out by the five elements of paper, method, data, results, and limitations, with a focus on how the identifying assumptions are defended.

### 5.1 Berry (1994)

::: {.case}
Paper and place in the history of methods: "Estimating Discrete-Choice Models of Product Differentiation," RAND Journal of Economics. This paper exposes the core idea of inversion: market shares can be inverted to recover mean utility, thereby turning demand estimation with only aggregate data into a regression with instrumental variables, and it is the starting point of all market-level demand estimation.

Method: Berry inversion. For plain logit it gives the closed-form inversion $\ln s_j - \ln s_0 = x_j\beta - \alpha p_j + \xi_j$, and for nested logit an inversion with one extra within-nest-share term, and points out that both price and within-nest share are endogenous and both need instruments. The paper also sketches a framework estimating the supply side's optimal-pricing first-order condition jointly with demand, paving the way for recovering marginal costs later.

Data: the paper demonstrates on product-level data from the U.S. automobile market, with each model's share, price, and characteristics. Automobiles are the classic testing ground for differentiated-product demand, with rich product characteristics and dispersed market shares.

Results: inversion plus instruments can consistently estimate demand parameters from market shares, and gives own-price elasticities and markups far more reasonable than approaches ignoring endogeneity. The paper established the basic paradigm of market-level demand estimation for the next two decades.

Limitations: the closed-form inversion holds only for logit and nested logit, which carry IIA or a preset nest structure and have restricted substitution. Generalizing the inversion to arbitrary random-coefficient models, and thereby fully relaxing substitution, requires a numerical method, which is exactly the contribution of BLP the next year. Half of this paper's historical significance is in exposing the inversion, and half is in marking the ceiling of the logit inversion.
:::

### 5.2 Berry, Levinsohn and Pakes (1995)

::: {.case}
Paper: "Automobile Prices in Market Equilibrium," Econometrica, that is BLP. It generalizes the Berry inversion to a general model with random coefficients, gives a contraction mapping to solve $\delta$ numerically, and puts demand together with the supply side's optimal pricing into GMM, the founding work of structural demand estimation.

Method: random-coefficient logit demand, contraction-mapping inversion, GMM estimation, paired with the first-order condition of supply-side multiproduct Bertrand pricing. Random coefficients let the substitution pattern be governed by product characteristics rather than shares, escaping IIA; the BLP instrument (rival characteristics) supplies the exogenous variation needed for identification.

Data: a product-level panel of the U.S. automobile market from 1971 to 1990, with each model's share, price, horsepower, size, fuel economy, and other characteristics.

Results: the estimated own-price elasticities and substitution patterns are far more reasonable than logit's, with substitution between similar models significantly stronger than the share-proportional prediction; estimating demand and supply together also recovers each model's marginal cost and evaluates market power. This method has since become the standard language of empirical IO.

Limitations: random coefficients, especially the one on price, are hard to identify and sensitive to both the instruments and the numerical implementation; mishandling the contraction tolerance or the number of simulation draws can bias the estimates. These implementation-level pitfalls gave rise to later work: Dubé-Fox-Su's MPEC, Reynaert-Verboven's optimal instruments, and software like pyblp that hardens best practice. BLP's greatness is half in the method and half in the long string of follow-up questions it opened.
:::

### 5.3 Nevo (2001)

::: {.case}
Paper: "Measuring Market Power in the Ready-to-Eat Cereal Industry," Econometrica. It is a model of BLP in practice, and it also fixed a key practical trick: use product and market fixed effects to absorb most of $\xi$, making identification more stable.

Method: random-coefficient logit demand, with random coefficients that include both interactions with observable demographics (observable heterogeneity) and pure unobservable heterogeneity. The paper uses brand fixed effects to absorb the product-level mean quality that does not vary across markets, concentrating the endogeneity burden onto the part of $\xi$ that varies across markets, and then uses instruments against it. The paper's appendix later became the textbook from which countless people learned BLP implementation.

Data: brand-level scanner data on ready-to-eat cereal across multiple cities and quarters, including shares, prices, characteristics, and city demographics. Cereal is the classic battlefield for differentiated-consumer-good demand, with many brands and a rich substitution structure.

Results: the estimated demand is quite elastic, and the backed-out markups show cereal makers do have considerable market power, but this power comes mainly from product differentiation and brand positioning rather than collusion. The paper uses structural demand to take apart the antitrust core question of "does a high price equal collusion."

Limitations: fixed effects, while absorbing a great deal of endogeneity, also eat up part of the variation needed to identify the price coefficient, raising the strength requirement on the remaining instruments; identifying the random coefficients is still tight. The paper accounts for these trade-offs clearly. It demonstrates all the craft BLP needs to go from method to credible empirics: fixed effects, instruments, and an honest accounting of the sources of identification.
:::

Put the three together and the point of anchoring becomes clear: Berry exposes the inversion and makes demand estimation from aggregate data possible, BLP uses random coefficients and the contraction mapping to fully relax substitution and connect the supply side, and Nevo polishes it into a credible empirical tool and demonstrates how to answer antitrust questions. Progress in market-level demand estimation has always turned on the same three threads: inverting shares into utility, finding exogenous variation in price, and getting the substitution pattern right.

## 6 A full walkthrough on Selva data

Now run all the tools of Section 4 through Selva's market-level data end to end. The code below uses R 4.5.3, fixing set.seed(16) for reproducibility, and every number cited in the text comes from an actual run of this code.

### 6.1 DGP

The design parameters are as follows: 300 city-week markets, six meal-kit brands sorted into high, medium, and low quality tiers, plus an outside option; the baseline is six single-product firms (each brand a firm of its own, and the next chapter merges the firms of brands 1 and 2). True demand is random-coefficient logit, with a random coefficient on quality $x_1$ and one on price; true supply is multiproduct Bertrand-Nash, with equilibrium price set by marginal cost plus markup and therefore correlated with the unobserved quality $\xi$.

```r
set.seed(16)
alpha <- 0.7; beta_x1 <- 1.6; beta_x2 <- 0.8           # demand: mean-utility parameters
sigma_x1 <- 1.1; sigma_p <- 0.15                        # random-coefficient standard deviations
# cost mc = gamma0 + gamma_w*w + gamma_x1*x1 + omega ; w = local input cost (instrument)
dt[, mc := 2.2 + 0.8*w + 0.7*x1 + omega]
# equilibrium price: multiproduct Bertrand-Nash fixed point p = mc + (H .* Delta)^{-1} s
```

The equilibrium is solved by finding a pricing fixed point for each market: given costs and demand, iterate $p = mc + (\mathcal{H} \odot \Delta)^{-1} s$ until convergence, where $\Delta$ is the negative matrix of demand price derivatives ($\Delta_{jk} = -\partial s_j / \partial p_k$, with positive diagonal) and $\mathcal{H}$ is the ownership matrix (the identity for single-product firms). The full derivation and use of this supply-side first-order condition is the subject of the next chapter; here it only generates the data. In the generated data, the median price is $7.6$ dollars, the median marginal cost $5.8$, the median markup $1.65$, the median Lerner index (markup as a share of price) $0.23$, the median inside share $0.099$, and the median outside share $0.18$. The correlation of price with unobserved quality $\xi$ is $0.23$ (the seed of endogeneity), and its correlation with the cost shifter $w$ is $0.61$ (the strength of the instrument).

### 6.2 OLS logit: reproducing the opening disaster

```r
dt[, y := log(s) - log(s0)]                              # Berry inversion: difference of log shares
ols <- feols(y ~ x1 + x2 + p, data = dt)
```

The price coefficient is $-0.19$, the sign is wrong, demand slopes upward, and the own-price elasticity is $+1.2$. The median marginal cost backed out from the single-product pricing FOC is $13.6$ dollars, far above the $7.6$ selling price. This is the disaster of Section 1, and the root cause is price endogeneity.

### 6.3 IV logit: curing endogeneity

```r
iv <- feols(y ~ x1 + x2 | p ~ w + sumx1_riv + sumx2_riv, data = dt)
```

Price is instrumented by the cost shifter $w$ and the sum of rival characteristics. The first-stage coefficient on $w$ is $0.696$ with a partial-F of $1520$, a strong instrument. The price coefficient returns to $0.51$, demand slopes down, the median own-price elasticity is $-3.3$, and the backed-out median marginal cost rises to $5.3$ dollars. The sign is right, but IIA leaves the elasticities still too small and the substitution still distorted.

### 6.4 Full BLP: inner contraction plus outer GMM

```r
# inner: for fixed theta2, solve delta by contraction mapping
contract <- function(mu, s_obs, d0){
  d <- d0
  repeat { d_new <- d + log(s_obs) - log(pred_share(d, mu))
           if (max(abs(d_new - d)) < 1e-12) break; d <- d_new }
  d_new
}
# outer: minimize the GMM objective over (sigma_x1, sigma_p), concentrating out the linear parameters by 2SLS
gmm_obj <- function(theta2){
  delta <- delta_theta(theta2)                          # one contraction per market
  xi    <- delta - Xlin %*% solve_2sls(Xlin, delta, Z)  # recover the residual
  g <- t(Z) %*% xi; as.numeric(t(g) %*% Wmat %*% g)
}
```

Estimation results (truth in parentheses): the price coefficient $0.751$ ($0.70$, standard error $0.075$), the quality coefficient $\beta_{x_1} = 1.661$ ($1.6$, SE $0.115$), $\beta_{x_2} = 0.711$ ($0.8$, SE $0.075$), the quality random coefficient $\sigma_{x_1} = 1.21$ ($1.1$, SE $0.281$), and the price random coefficient $\sigma_p = 0.033$ ($0.15$, SE $0.601$). Note that the standard error of $\sigma_p$ is far larger than the point estimate, so it is barely identified, exactly the faithful reflection of the difficulty of identifying the price random coefficient discussed in Section 3.5. But the quantities that really must be delivered are estimated quite accurately: the median own-price elasticity is $-4.7$ (truth $-4.4$) and the backed-out median marginal cost is $5.94$ dollars (truth $5.81$), with none negative.

![Recovered price coefficients: OLS gives the wrong negative sign (upward-sloping demand), IV logit swings it back to positive but too small because of IIA, and full BLP hits almost exactly the truth of 0.70. The dashed line is the truth.](assets/fig/fig_20_recovery.svg)

### 6.5 Substitution patterns: BLP recovers the near-neighbor substitution logit throws away

The most important difference among the three estimators is not the value of the price coefficient but the substitution pattern. Put side by side the diversion ratios of demand flowing to each brand after brand 1 (the highest quality) raises its price:

| Flows to | Truth | BLP | logit (IIA) |
|---|---|---|---|
| $\to$ Brand 2 (second-highest quality, near neighbor) | 0.271 | 0.278 | 0.206 |
| $\to$ Brand 3 | 0.197 | 0.200 | 0.164 |
| $\to$ Brand 4 | 0.190 | 0.194 | 0.162 |
| $\to$ Brand 5 | 0.149 | 0.150 | 0.137 |
| $\to$ Brand 6 | 0.100 | 0.100 | 0.096 |

Brand 1 and brand 2 are closest in quality, near neighbors of each other, and the true substitution ratio is $0.271$. BLP recovers this strong near-neighbor substitution ($0.278$, almost exactly the truth), whereas plain logit, constrained by IIA, allocates only by shares and gives $0.206$, systematically too low. This difference is not academic fastidiousness: when the next chapter simulates merging the firms of brand 1 and brand 2, it is exactly this diversion ratio that determines whether the merger raises price significantly. Using logit's $0.206$ underestimates the merger's upward price pressure, while BLP's $0.278$ comes close to the truth.

![Diversion ratios after brand 1 raises its price. Brand 1 and brand 2 are closest in quality, with the highest true substitution ratio (navy); BLP (blue) recovers this strong near-neighbor substitution, while plain logit (red), constrained by IIA, systematically underestimates it.](assets/fig/fig_20_diversion.svg)

### 6.6 The full reconciliation

Put the key products of the three estimators side by side with the truth:

| Estimator | Price coefficient $\alpha$ | Own-price elasticity (median) | Backed-out marginal cost (median) |
|---|---|---|---|
| Truth (DGP) | 0.700 | -4.4 | \$5.81 |
| OLS logit | -0.19 | +1.2 | \$13.6 |
| IV logit | 0.51 | -3.3 | \$5.28 |
| Full BLP | 0.751 | -4.7 | \$5.94 |

Three rows tell the entire main thread of this chapter. OLS gets even the sign of the demand slope wrong and backs out a marginal cost above the selling price, purely the work of endogeneity. IV logit uses a cost instrument to cure endogeneity, and demand slopes down, but plain logit's IIA makes it underestimate elasticities and flatten near-neighbor substitution. Full BLP both cures endogeneity and relaxes substitution, with the price coefficient, elasticities, and marginal costs all close to the truth, and the substitution pattern also recovers the near-neighbor structure logit throws away. These differences all come from model specification and identification strategy, not from sample noise. This backed-out marginal cost and substitution structure is exactly the raw material for the next chapter's supply-side analysis.

## 7 Failure modes and robustness of identification

In the simulation the identifying assumptions are constructed, while in real research they may fail at any moment. This section reviews the most common ways they fail and the operable responses.

The validity of the price instrument is the foundation of the whole building. A cost shifter qualifies only when it truly enters cost but not utility; a variable that both raises cost and lifts the platform's overall popularity (a peak season, say) enters both sides of supply and demand, and the exclusion restriction breaks. The BLP instrument's weakness is Armstrong's (2016) large-market degeneracy: once the number of products per market is large, the markup tends to a constant, rival characteristics lose their power to move price, and it becomes a weak instrument. Diagnostically, look at first-stage strength; a weak instrument makes the price coefficient's variance explode and can bias it worse than OLS. When a digital product's marginal cost is nearly zero, a good cost shifter is especially hard to find, a real obstacle to applying BLP to platform data.

The identification of random coefficients is the second weak link and the most easily glossed over. $\sigma$ is identified by how substitution varies with characteristics, and if product characteristics have insufficient cross-market variation, or a differentiation instrument is lacking, $\sigma$ collapses toward zero and quietly pulls the model back to IIA. The random coefficient on price $\sigma_p$ is especially hard, and in Selva its standard error is far larger than its point estimate, barely identified. The honest thing is to report $\sigma$'s standard error and the strength of the instrument identifying it, rather than only a point estimate that looks nonzero. If $\sigma$ really cannot be identified, frankly admitting the data only support logit-level substitution conclusions is better than forcing an untrustworthy random coefficient.

Many-instrument and many-weak-instrument bias is an invisible trap in market-level BLP. To identify the random coefficients, people tend to stuff many functions of rival characteristics into the instrument set, but when the instruments are many and each individually weak, the second-stage estimate biases toward OLS, and OLS is exactly the one that gives the wrong sign. Selva's estimation showed this bias when the number of markets was small, with the price coefficient pulled down, only returning near the truth after the number of markets was increased. The response is to compress the dimension with principal components or optimal instruments, not to pile up instruments mindlessly.

The numerical implementation itself manufactures bias. If the contraction-mapping tolerance is set loose, the error of the inner inversion seeps into the outer GMM and biases $\theta_2$, exactly the problem MPEC aims to sidestep; if the number of simulation draws is too small, simulation error contaminates the estimates, and the smaller the share the more severe it is. The robust practice is to set the tolerance extremely tight, the number of draws large enough, and start from multiple initial values to confirm it is not a local solution. These are not optional tuning knobs but part of whether BLP is credible.

There are several more easily overlooked gates. Zero shares, if simply discarded, introduce selection bias, worst in long-tail categories, and need the Gandhi-Lu-Shi treatment. The definitions of the outside option and market size determine the level of all shares and elasticities, and defining them too large or too small biases systematically. On functional form, whether price enters utility as $p$ or as $p / y$ ($y$ being income) implies different substitution and welfare content, and needs a mechanistic reason. The last gate is left to the next chapter: when backing out marginal cost from demand, the pricing FOC assumes some competitive behavior (conduct), and if the true competition is not the assumed kind, the backed-out cost is wrong, which ties demand estimation to the behavioral assumptions of the supply side.

Stringing these failure modes together, BLP's credibility comes down to two things: whether exogenous variation in price can be found, and whether the way substitution varies with characteristics can be identified. Instrument strength and exclusion, the standard errors of the random coefficients, the rigor of the numerical implementation, and the handling of zero shares and market definition are all evidence provided around these two, and none of them can be replaced by "the GMM converged and the coefficients look good."

## 8 Further reading

::: {.readings}
Required reading, in suggested order:

- Berry (1994, RAND Journal of Economics). The origin of the inversion idea; read it first to understand why market shares can be inverted to utility and why price needs an instrument.
- Berry, Levinsohn and Pakes (1995, Econometrica). The founding work of BLP, with random coefficients, the contraction mapping, and the supply side as one, the source of structural demand estimation.
- Nevo (2001, Econometrica). BLP's model of practice, with the practicalities of fixed effects and instruments, and an appendix that is the best textbook for learning the implementation.
- Conlon and Gortmaker (2020, RAND Journal of Economics). The best-practice guide for pyblp, hardening twenty years of pitfalls into operable norms, required reading before doing BLP by hand.
- Gandhi and Houde (2019). The differentiation instrument, the modern way to identify random coefficients, the direct source of Section 3.3 of this chapter.

Further reading:

- Berry and Haile (2014, Econometrica). The foundation for nonparametric identification of demand from market-level data, with connected substitutes and completeness, for understanding why BLP is credible.
- Dubé, Fox and Su (2012, Econometrica). The contrast of MPEC and NFXP, the classic example of one estimator and two algorithms.
- Armstrong (2016, Econometrica). The critique of the BLP instrument degenerating in large markets, and why you cannot do without a cost shifter.
- Petrin (2002, Journal of Political Economy). The micro-moment example of bringing individual information into BLP.
- Gandhi, Lu and Shi (2023, Quantitative Economics). The problem of zeroes in market shares and its treatment.
:::

::: {.apa-refs}
- Armstrong, T. B. (2016). Large market asymptotics for differentiated product demand estimators with economic models of supply. *Econometrica, 84*(5), 1961-1980. https://doi.org/10.3982/ECTA10600
- Berry, S. T. (1994). Estimating discrete-choice models of product differentiation. *RAND Journal of Economics, 25*(2), 242-262. https://doi.org/10.2307/2555829
- Berry, S. T., & Haile, P. A. (2014). Identification in differentiated products markets using market level data. *Econometrica, 82*(5), 1749-1797. https://doi.org/10.3982/ECTA9027
- Berry, S., Levinsohn, J., & Pakes, A. (1995). Automobile prices in market equilibrium. *Econometrica, 63*(4), 841-890. https://doi.org/10.2307/2171802
- Conlon, C., & Gortmaker, J. (2020). Best practices for differentiated products demand estimation with PyBLP. *RAND Journal of Economics, 51*(4), 1108-1161. https://doi.org/10.1111/1756-2171.12352
- Conlon, C., & Gortmaker, J. (2023). Incorporating micro data into differentiated products demand estimation with PyBLP (SSRN Working Paper No. 4553609). https://doi.org/10.2139/ssrn.4553609
- Dubé, J.-P., Fox, J. T., & Su, C.-L. (2012). Improving the numerical performance of static and dynamic aggregate discrete choice random coefficients demand estimation. *Econometrica, 80*(5), 2231-2267. https://doi.org/10.3982/ECTA8585
- Gandhi, A., & Houde, J.-F. (2019). *Measuring substitution patterns in differentiated-products industries* (Working Paper No. 26375). National Bureau of Economic Research.
- Gandhi, A., Lu, Z., & Shi, X. (2023). Estimating demand for differentiated products with zeroes in market share data. *Quantitative Economics, 14*(2), 381-418. https://doi.org/10.3982/QE1593
- Nevo, A. (2001). Measuring market power in the ready-to-eat cereal industry. *Econometrica, 69*(2), 307-342. https://doi.org/10.1111/1468-0262.00194
- Petrin, A. (2002). Quantifying the benefits of new products: The case of the minivan. *Journal of Political Economy, 110*(4), 705-729. https://doi.org/10.1086/340779
:::
