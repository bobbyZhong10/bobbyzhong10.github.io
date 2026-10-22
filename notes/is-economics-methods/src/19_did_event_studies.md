---
title: "Difference-in-Differences & Event Studies"
subtitle: "From the 2×2 Table to Staggered Adoption"
seriesline: "Foundations of Information Systems Economics · Chapter 19"
series: "Foundations of Information Systems Economics: Empirical Methods for Information Systems Economics"
sidebartitle: "Foundations of Information Systems Economics"
sidebarsub: "Theory · Methods · Applications"
chapterlabel: "Chapter 19 · Difference-in-Differences & Event Studies"
---

## Introduction

Kestrel Marketplace rolls out a new policy city by city, in batches. Cities grow faster once the policy is live, and the event-study plot shows the familiar rise after treatment. The policy seems to have worked. But the cities did not go live on the same day. The early-treated cities end up serving as the control group for the late-treated ones, while the policy effect keeps accumulating with time since launch. Here a standard two-way fixed effects (TWFE) regression can blend several local comparisons, each with a different meaning, into one number whose sign is not even reliable.

The original appeal of DiD comes from a table with two rows and two columns: the before-after change for the treated group, minus the contemporaneous change for the control group. This difference of differences removes fixed differences between the groups and common time shocks, but it does not conjure a counterfactual out of thin air. It still has to ask: had the policy not happened, how would the two groups have moved relative to each other? "The pre-treatment trends look about the same" is supporting evidence, not the parallel trends assumption itself.

This chapter first makes that counterfactual explicit in the most transparent case, the 2x2, and then moves to staggered adoption (treatment in batches), which is far more common in platform research. Every identification result comes with numbered assumptions, a theorem, and a proof sketch, followed by a plain sentence on whose effect it estimates and under which assumptions. The Goodman-Bacon decomposition and the de Chaisemartin-D'Haultfœuille weight diagnostic reveal, from two sides, exactly which comparisons get mixed into a traditional TWFE coefficient. Callaway-Sant'Anna, Sun-Abraham, Borusyak-Jaravel-Spiess, and the $DID_\ell$ of de Chaisemartin-D'Haultfœuille then reorganize the control groups and event time so that every weight has a clear meaning. Finally, the Rambachan-Roth sensitivity analysis asks how wrong parallel trends must be before the conclusion flips. The chapter assumes you know the potential outcomes of Chapter 15 and the fixed-effects panel models of Chapter 10; inference uses the cluster-robust standard errors of Chapter 8. The question running through the chapter is not "which command should I run" but: at each point in time, who plays the never-observed counterfactual on behalf of the treated group?

## 1 A Number That Does Not Add Up

::: {.case}
Kestrel Marketplace is a fictional two-sided e-commerce platform. We tell its story with simulated data. The payoff is a fully known DGP: every estimator can be checked against the truth, a teaching luxury that real data never provide. Here is the setup. The platform operates in 60 cities over a 36-month observation window, and the outcome is the average log GMV of sellers in each city and month. The platform cuts the seller commission rate from 12% to 8% in three batches: 17 cities in month 13, another 17 in month 19, and 18 more in month 25, while the remaining 8 cities keep the 12% rate throughout the window. Management's question is direct: how large is the effect of the commission cut on seller GMV?

Because this is a simulation, we know the answer. The true effect is positive and grows with event time: about 0.02 (log points) in the month of the cut, climbing to around 0.14 a year later, and drifting to about 0.18 at the longest horizon ($e = 23$). Averaging the realized effect over all 930 post-treatment city-month cells gives a true simple ATT of 0.088.
:::

![The raw path of city-average log GMV by calendar month for each cohort, with dashed lines marking the timing of the three commission cuts (months 13, 19, and 25).](assets/fig/fig_19_raw_trends.svg)

Start with the raw data. Three features of the figure are worth noting. First, the cities that cut commissions early sit at a higher level. This is by design: the platform moves first in its large markets, and a level difference by itself does not threaten identification. Second, the four groups of cities share an obvious common growth trend and seasonal pattern. Third, after the cuts the groups gradually fan out, and each late-treated group also starts to accelerate once it passes its own treatment date. Intuitively, this looks like a policy whose effect is positive and accumulates over time.

So the researcher runs the textbook TWFE regression:

$$Y_{it} = \alpha_i + \lambda_t + \beta\, D_{it} + \varepsilon_{it}$$

where $D_{it}$ equals 1 once city $i$ has cut its commission by month $t$, and standard errors are clustered at the city level. The result is $\hat\beta_{TWFE} = 0.0227$ (SE 0.0087, $t = 2.61$, $p = 0.011$; we round it to 0.023 when quoting it later). Taken on its own, the result looks clean: positive and statistically significant, with a magnitude of 2.3%. In a report to management it would read "the commission cut works, but the effect is small and probably not worth rolling out."

The problem is that we know the truth is 0.088. Static TWFE recovers only 25.8% of the true effect; nearly three quarters of it has vanished. Nothing is wrong with how the regression is written: the fixed effects are all there, the standard errors are clustered at the right level, and parallel trends holds exactly in the data (the DGP guarantees it). So what is $\hat\beta_{TWFE}$ actually estimating? This is the main question of the chapter. Section 4 gives the answer: under staggered adoption combined with effects that grow with event time, static TWFE estimates a weighted average with strange weights, and part of the weight sits on comparisons that use already-treated units as controls. Section 6 takes this 0.023 apart precisely and shows where each piece of the bias comes from.

## 2 The Economic Model and the Estimand

Before estimating anything, we need to settle two things: which causal quantity we want, and what language writes it down precisely. The first is this section's estimand; the second is the notation and model that express it, namely potential outcomes and the additive structure on $Y(0)$ at the end of this section.

The core difficulty of causal inference is the counterfactual. For a city that has cut its commission, we can never observe what would have happened had it not. Potential outcomes are the notation designed to describe such parallel-world outcomes. This section sets them up first, because all the identification and estimation that follow rest on them.

Start with notation. Denote a city by $i$ and a month by $t$. The cohort variable $G_i$ is the period in which city $i$ is first treated; in Kestrel $G_i \in \{13, 19, 25\}$, and never-treated cities are coded $G_i = \infty$ (0 in the data). The treatment indicator is $D_{it} = \mathbf{1}\{t \geq G_i\}$. Staggered adoption means that treatment is irreversible, $D_{it} \geq D_{i,t-1}$: once the commission is cut, it does not go back up. Event time $e = t - g$ aligns the cohorts on their own treatment dates.

Potential outcomes are indexed by the period of first treatment: $Y_{it}(g)$ is the outcome that city $i$, if first treated in period $g$, would realize in period $t$, and $Y_{it}(0)$ is the never-treated path. The observed outcome is

$$Y_{it} = \mathbf{1}\{G_i > t\}\, Y_{it}(0) + \mathbf{1}\{G_i \leq t\}\, Y_{it}(G_i)$$

The basic estimand of this chapter is the group-time average treatment effect:

$$ATT(g,t) = \mathbb{E}\big[Y_{t}(g) - Y_{t}(0) \mid G = g\big]$$

that is, the ATT of cohort $g$ at calendar time $t$. Management's question, "how large is the effect of the commission cut on GMV," is not really one number but some average over the whole $ATT(g,t)$ matrix. You can average by cohort (which batch of cities benefits most), by event time (how fast the effect builds up and whether it fades), or all the way down to a single total. First decide which average you want and with what weights, then choose the estimator, never the other way around. The modern estimators of Section 4 share one structure: "estimate $ATT(g,t)$ first, then aggregate with explicit weights." The root problem of static TWFE is exactly that its aggregation weights are neither chosen by the researcher nor economically interpretable.

The **parallel trends** assumption is the identification core of DiD. It pays to understand its economic content first; the formal statement comes in Section 3. It says: had the treated group not been treated, its outcome would have changed over time in the same way as the control group's. What kind of outcome-generating process guarantees this? Write the untreated state as $Y_{it}(0) = \alpha_i + \lambda_t + u_{it}$, where $\alpha_i$ absorbs a city's intrinsic level and $\lambda_t$ is a common time shock. As long as the group mean paths of $u_{it}$ are the same, parallel trends holds. In other words, this is an additive separability assumption about the process generating $Y(0)$: unobserved factors enter additively, and time shocks hit the treated and control groups by the same amount. This is far stronger than "I controlled for fixed effects." Writing fixed effects into a regression is a mechanical operation. Parallel trends is an economic assertion about the counterfactual world: how the treated cities' GMV would have moved had they not cut commissions. Regression syntax cannot guarantee that assertion.

Additivity has a corollary: parallel trends depends on functional form, so "log or level" is a choice between two identification assumptions. Section 3.2 develops this point.

Platform data give us special reason to worry about parallel trends. A platform does not draw lots to decide which cities get the cut first. The first to be treated are often the fastest-growing, most competitive, or strategically most important markets, so treatment timing becomes entangled with growth trajectories. Two things must be kept apart here. Parallel trends allows the two groups to sit at different levels, and it allows the platform to move first in the cities with the largest expected effects. The only thing it rules out is that the two groups' growth paths in the untreated state differ systematically. A platform that picks cities by growth potential may be picking exactly the cities whose growth paths differ to begin with, and that is the real threat. Section 3 dissects the three identification assumptions of DiD one by one.

::: {.intuition}
DiD is easier to see when placed back on the menu of ways to construct a counterfactual (matching is covered in Chapter 15 and synthetic control in Chapter 22; here you only need to see whom each one borrows its counterfactual from). To evaluate a treatment, you always have to borrow a counterfactual from somewhere. A before-after comparison borrows your own past. Matching borrows other units that resemble you on observable characteristics. DiD borrows other units' change as the counterfactual for your own change. Synthetic control borrows the change of a weighted combination of other units. DiD and matching are mirror images. Matching makes almost no functional-form assumption but requires selection on observables, and unobserved selection kills it. DiD makes a strong additive functional-form assumption and in exchange tolerates arbitrary unobserved heterogeneity in levels. Conditional parallel trends (Section 3) then applies the matching logic to the growth of the outcome rather than its level. Choosing a design means choosing which assumption you are willing to pay for your counterfactual.
:::

To summarize: the implicit model of DiD is an additively separable structure for $Y(0)$, the estimand is the ATT family (generally not the ATE), and under staggered adoption the basic object is $ATT(g,t)$. This section has not touched "how to estimate" at all, and that is deliberate. Whether the identification assumptions are credible depends on the design and the institutional background, not on how the regression is written, although different estimators do require parallel trends of different strength (Section 4.5).

## 3 Identification: Which Bias Each of the Three Assumptions Blocks

Identification logically comes before estimation. Identification asks whether, under a set of assumptions about the counterfactual world, the causal quantity we want can be written as a function of the distribution of observable data. The answer depends only on the design and the institutional background. It has nothing to do with whether we later use OLS or GMM, or with whether the sample is large or small. This section first proves identification line by line in the 2x2 setup, then dissects the three assumptions one at a time, and then extends the argument to $ATT(g,t)$ under staggered adoption and to triple differences. The 2x2 case deserves this much space because every difficulty of staggered adoption is a rearrangement of this skeleton.

### 3.1 The 2x2 Setup and the Identification Theorem

Start with the most naive comparison any platform operator would make. To learn whether the commission cut worked, the most direct approach is to look at the change in GMV in the treated cities before and after the cut. This is called a single difference. Its flaw is obvious. Even without the cut, GMV would move with the platform's overall trajectory, the seasons, and the macro environment. Subtracting before from after folds all of these policy-irrelevant time changes into the effect: the estimate is too high when the market is rising and too low when it is falling. To strip out the time change, you need a set of control cities that did not cut commissions over the same period. Their before-after change works like a ruler, measuring how much the market itself rose over that stretch. The treated cities' before-after change minus the control cities' contemporaneous before-after change leaves the net effect of the policy. This one step, a difference of differences, is the whole intuition of difference-in-differences.

A minimal numerical example makes this concrete. Suppose the treated cities' monthly GMV rises from 120 before the cut to 150 after, a gain of 30, while over the same period a set of control cities rises from 100 to 110, a gain of 10. The two groups start from different levels, since the treated cities are the larger market, but that does not matter: DiD cares about how much each group rose, not whose starting point is higher. The control cities rose by 10 without any cut, which says the market itself contributed about 10 over this stretch. Subtracting this 10 from the treated cities' 30 gives 20, the part attributable to the commission cut. Looking only at the treated cities would overstate the effect at 30, precisely because the market was not netted out. All of DiD's credibility rests on one question: is the 10 that the control cities rose the same 10 that the treated cities would have risen without the cut?

In the potential outcomes notation of Section 2, the 20 in this example is the average treatment effect on the treated, and the hard part is the counterfactual $Y_{i2}(0)$. Now write the 2x2 case out formally. There are two periods $t \in \{1, 2\}$; units in the treated group $D_i = 1$ are treated at $t = 2$, and the control group $D_i = 0$ is never treated. In the 2x2 case, group and treatment coincide: treatment status is $D_{it} = D_i \cdot \mathrm{Post}_t$ with $\mathrm{Post}_t = \mathbf{1}\{t = 2\}$, and in this section $D_i$ denotes only the group. With two periods, the potential outcomes are abbreviated $Y_{it}(1)$ and $Y_{it}(0)$. The target is $ATT = \mathbb{E}[Y_{i2}(1) - Y_{i2}(0) \mid D_i = 1]$. DiD substitutes the control group's before-after change for the treated group's missing counterfactual change; whether the substitution is valid depends on the three assumptions below.

::: {.assumption}
**A1 (parallel trends)** The treated and control groups have the same expected growth in their untreated outcomes:

$$\mathbb{E}\big[Y_{i2}(0) - Y_{i1}(0) \mid D_i = 1\big] = \mathbb{E}\big[Y_{i2}(0) - Y_{i1}(0) \mid D_i = 0\big]$$
:::

::: {.assumption}
**A2 (no anticipation)** The pre-treatment outcome is unaffected by future treatment: for the treated group, $Y_{i1} = Y_{i1}(0)$.
:::

::: {.assumption}
**A3 (SUTVA / no spillovers)** Each unit's outcome depends only on its own treatment status, and treatment does not spill over; in particular, the control group's observed outcome is its $Y(0)$.
:::

Each of the three assumptions has one job: A1 governs the counterfactual trend, A2 keeps the base period clean, and A3 keeps the control group clean. Sections 3.2 to 3.4 take them one at a time: what each constrains, what it allows, and how it fails.

::: {.theorem}
**Theorem 1 (2x2 identification)** Under A1, A2, and A3,

$$\Delta_{DD} \equiv \big(\mathbb{E}[Y_{i2} \mid D_i = 1] - \mathbb{E}[Y_{i1} \mid D_i = 1]\big) - \big(\mathbb{E}[Y_{i2} \mid D_i = 0] - \mathbb{E}[Y_{i1} \mid D_i = 0]\big) = ATT$$
:::

The proof sketch has three steps, each using exactly one assumption. Write $\Delta Y_i = Y_{i2} - Y_{i1}$.

Step one (A2): for the treated group, A2 guarantees $Y_{i1} = Y_{i1}(0)$; then split $Y_{i2} = Y_{i2}(1)$ into effect plus counterfactual:

$$\mathbb{E}[\Delta Y_i \mid D_i = 1] = \underbrace{\mathbb{E}\big[Y_{i2}(1) - Y_{i2}(0) \mid D_i = 1\big]}_{ATT} + \underbrace{\mathbb{E}\big[Y_{i2}(0) - Y_{i1}(0) \mid D_i = 1\big]}_{\gamma(1)}$$

Step two (A3): for the control group, A3 guarantees that its observed outcome is the $Y(0)$ path throughout:

$$\mathbb{E}[\Delta Y_i \mid D_i = 0] = \underbrace{\mathbb{E}\big[Y_{i2}(0) - Y_{i1}(0) \mid D_i = 0\big]}_{\gamma(0)}$$

Step three (A1): subtracting the two expressions gives

$$\Delta_{DD} = ATT + \big(\gamma(1) - \gamma(0)\big)$$

A1 says $\gamma(1) = \gamma(0)$, so the bias term vanishes, and identification is complete. Replacing the four population means with sample means gives consistency by the law of large numbers when units are sampled independently; the regression form is in Section 4.1.

In words: under A1 through A3, the double difference of four means estimates the treated group's own ATT, with every treated unit weighted equally; the control group contributes only a ruler, not an effect. The closest thing to this 2x2 table in Kestrel is the Goodman-Bacon component that compares cohort 13 with the 8 never-treated cities, before versus after treatment. It gives 0.095, against a target of cohort 13's average effect over its 24 post-treatment months, whose true value is 0.106. Parallel trends holds exactly in the DGP, so this small gap is sampling noise from having only 8 control cities, not bias.

The decomposition $\Delta_{DD} = ATT + (\gamma(1) - \gamma(0))$ is also the skeleton of the whole chapter. The term $\gamma(1) - \gamma(0)$ is the difference between the two groups' counterfactual trends, and given A2 and A3, DiD identifies the ATT if and only if it is zero. It is most likely to be nonzero when the treated group would have grown faster anyway. Every failure and repair that follows is about when this term is nonzero and how to push it back to zero.

### 3.2 Anatomy of Parallel Trends

The name parallel trends is easiest to grasp as a picture. Draw the average GMV of the treated group and of the control group as two lines over time. DiD never requires the two lines to sit at the same height (the treated cities are large markets and naturally sit higher, which is fine). It requires that, in the hypothetical world where nobody is treated, the two lines move in parallel. The trouble is that in real data we see only the control group's line and a short stretch of the treated group's line before the cut. After treatment, the path the treated cities would have followed without the cut is a dashed line that cannot be drawn. Parallel trends assumes exactly that this invisible dashed line stays parallel to the control group's line.

![Geometric intuition for parallel trends: the control group (navy) and the treated group's pre segment (red) move in parallel; after treatment, the gap between the treated group's actual path (solid red) and its counterfactual path absent the cut (dashed red, parallel to the control group) is the effect ATT.](assets/fig/fig_19_parallel_intuition.svg)

Now turn this geometric intuition into math precise enough to be right or wrong. Write the untreated state in additive form

$$Y_{it}(0) = \alpha_i + \lambda_t + u_{it}$$

where $\alpha_i$ is a city's intrinsic level, $\lambda_t$ is a time shock that treats all units alike, and $u_{it}$ is the remaining idiosyncratic shock. Both $\alpha_i$ and $\lambda_t$ cancel in the differencing, so A1 ultimately falls only on the idiosyncratic term and is equivalent to $\mathbb{E}[u_{i2} - u_{i1} \mid D_i = 1] = \mathbb{E}[u_{i2} - u_{i1} \mid D_i = 0]$. This equality draws the line between what parallel trends allows and what it forbids. It allows $\alpha_i$ to differ systematically between the two groups: treated cities being intrinsically large markets ($\mathbb{E}[\alpha_i \mid D=1] \neq \mathbb{E}[\alpha_i \mid D=0]$) is no obstacle at all. This is exactly the point Wooldridge (2010) §10.5.6 makes about fixed effects in policy analysis: the policy variable may be arbitrarily correlated with the persistent component of the error, as long as it is unrelated to the part that deviates from the unit mean. It allows $\lambda_t$ to be arbitrarily large and arbitrarily nonlinear, as long as it hits both groups by the same amount. It even allows $u_{it}$ to have arbitrary individual heterogeneity and fluctuation over time, as long as the average changes of these fluctuations are equal across groups. The one thing it forbids is a systematic deviation of the treated group's average growth path in the untreated state from the control group's.

This point separates parallel trends from selection on observables (unconfoundedness). The latter requires the levels of the treated and untreated states to be comparable after conditioning on covariates, and so forbids $\alpha_i$ from being correlated with treatment. Parallel trends is far more permissive on the level dimension: it allows arbitrary fixed unobserved heterogeneity, even heterogeneity that is highly correlated with treatment. The price is that on the trend dimension it is a knife-edge assumption, and the slightest divergence in group trends breaks it. If the worry is that the treated group is simply a different kind of unit, DiD often fits better than matching; if the difference between groups looks more like a level gap driven by observable characteristics, the matching logic is more natural.

Parallel trends depends on functional form. This is not a technical detail but an identification choice to be faced head-on. If A1 holds for log GMV, it generally does not hold for the GMV level. Parallel logs mean that the two groups grow by the same percentage, and when their starting points differ, the same percentage corresponds to different absolute increments. A minimal example: one city's GMV rises from 100 to 110 and another's from 200 to 220. In levels the gains are 10 versus 20, yet in logs both grow by about 10%: levels are not parallel, logs are. Roth and Sant'Anna (2023, Econometrica) give a sharp characterization: parallel trends holds for every monotone transformation of the outcome if and only if a very strong condition is met. Roughly, the population must split into two parts: in one part treatment is as good as randomly assigned (the distribution of $Y(0)$ is the same in both groups), and in the rest the distribution of $Y(0)$ does not shift over time. Outside this knife-edge case, choosing log or level is a choice between two identification assumptions, and also between two estimands: log corresponds to a percentage effect, level to an effect in absolute terms. The choice should follow the mechanism. The commission rate acts proportionally on sellers' marginal incentives, so log has a mechanistic justification.

Finally, A1 does not require the treatment effect to be the same across groups. A platform that deliberately cuts commissions first in the cities with the largest expected effects is engaging in selection on gains, and this does not violate parallel trends at all, because parallel trends constrains the path of $Y(0)$, not the size of the effect. The only price is that what gets identified is the treated group's own ATT; extrapolating it to "what would happen if we also cut commissions in the control cities" has no basis. This gap between the ATT and the ATE recurs throughout the staggered setting.

### 3.3 No Anticipation and Shifting the Base Period

No anticipation is a plain idea: before the treatment actually lands, the units about to be treated should not change their behavior just because they sense it coming. If the commission cut is announced in advance, sellers may hold back orders and release them in the low-commission period. Base-period GMV is then artificially depressed, the first difference counts the rebound as effect, and the effect is overstated. Conversely, if sellers ramp up advertising after the announcement but before the cut takes effect, building volume for the low-commission period, the base period is inflated and the effect is understated. Either way, the base period is no longer a world undisturbed by the policy.

The remedy is to move the base period earlier and use as the baseline only earlier periods that are surely clean. To build this into an identification result, we need a version of the assumption that allows limited anticipation. Let the number of anticipation periods be $\kappa \geq 0$. Callaway and Sant'Anna (2021) write it as $\delta$, but in this chapter $\delta$ is reserved for the Sun-Abraham coefficients of Section 4.5, so we write $\kappa$ instead. It has nothing to do with Abadie's complier weighting $\kappa$ in Chapter 18; in this chapter it denotes only the number of anticipation periods. The notation follows the cohort indexing of Section 2, whose observation equation in fact presumes no anticipation. Once anticipation is allowed, the outcome observed for cohort $g$ in every period is $Y_t(g)$.

::: {.assumption}
**A2' (limited anticipation)** There is a known $\kappa \geq 0$ such that for every cohort $g$ and all $t < g - \kappa$: $\mathbb{E}[Y_t(g) \mid G = g] = \mathbb{E}[Y_t(0) \mid G = g]$. Setting $\kappa = 0$ gives the multi-period version of A2.
:::

A2' constrains the length of anticipation, not its direction: inside the window the anticipatory response can be positive or negative, and outside it must be zero. Its observable footprint is an abnormal climb or dip in the few pre coefficients right next to the treatment date on the event-study plot. How $\kappa$ enters the identification formula has to wait until Section 3.6 writes out the multi-period setup: Theorem 2 puts the base period at $g - \kappa - 1$, and Corollary 1 answers what is estimated when the base period is chosen wrongly.

### 3.4 SUTVA, Spillovers, and Composition

On a platform, SUTVA has a plain meaning: a city's outcome should be determined only by whether it is treated, not by whether other cities are. A3 splits into two conditions. One is no interference: a unit's outcome is unaffected by whether others are treated. The other is a single version of treatment: there is no mixture of treatments of different intensities all going under the name "commission cut." For DiD, the troublesome condition is no interference. Once it fails, the control group's observed outcome is no longer a clean $Y(0)$, and the $\gamma(0)$ that serves as the counterfactual trend in the proof of Theorem 1 has treatment spillover mixed into it.

On platforms the channels of spillover are very concrete, and each calls for a different diagnostic. Sellers can move their registration to a city that cut commissions, so sellers are pulled straight out of the control cities. Buyers can order across cities, so low prices in treated cities siphon demand from control cities. The platform's traffic-allocation algorithm may pull impressions away from control cities because conversion rates rise in treated cities, a spillover mediated purely by the algorithm.

Behind spillovers lies a more fundamental distinction: DiD estimates a partial-equilibrium contrast. It takes for granted that the control group supplies the counterfactual for the treated group absent treatment, and this holds only if treatment does not change the equilibrium the control group lives in. If the commission cut changes the seller pricing equilibrium, the allocation of buyer traffic, or the competitive landscape of the whole platform, the control cities are themselves swept into the new equilibrium. The DiD framework then fails as a whole, and the question can be answered only with a structural model or a platform-level experiment. Which case you are in is decided by economic judgment about how the treatment is transmitted, not by any statistical test.

Adjacent to spillover but different in nature is the composition problem. Even if treatment does not spill over, as long as it changes the composition of the sample, the before-after comparison is no longer about the same set of units. If the cut attracts sellers from neighboring cities, the treated cities' seller population in the post-period is not the same as in the pre-period. Redefining the population as sellers who stay in the city throughout the window restores a clean comparison, but that is a different estimand: the effect on sellers who stay throughout.

### 3.5 Conditional Parallel Trends

Unconditional parallel trends is sometimes hard to believe, yet becomes more believable once we condition on pre-treatment covariates. The typical case is one where the growth path of the outcome depends systematically on certain city characteristics (population structure, logistics density, existing penetration), and these characteristics are imbalanced between the treated and control groups. Unconditional group trends then have no reason to be equal, but trends among cities of the same type should be.

::: {.assumption}
**A1' (conditional parallel trends)** Conditional on pre-treatment covariates $X$ (and with overlap: $P(D_i = 1 \mid X_i) < 1$),

$$\mathbb{E}\big[Y_{i2}(0) - Y_{i1}(0) \mid X_i, D_i\big] = \mathbb{E}\big[Y_{i2}(0) - Y_{i1}(0) \mid X_i\big]$$
:::

The right way to read A1' is as the matching assumption applied to a redefined outcome: the growth of the outcome rather than its level. This opens up matching's whole toolbox. You can give the time-invariant $X_i$ time-varying coefficients ($\beta_t' X_i$, which lets each type of city have its own trend), weight outcome growth by the propensity score, or combine the two into a doubly robust estimator (Section 4.6). Three rules must be kept. First, $X$ must take pre-treatment values; putting a variable that responds to treatment into the conditioning set is the bad control error discussed in Section 7. Second, conditioning can repair only trend differences driven by observable characteristics; if the treated and control groups differ in unobserved drivers of growth, no number of covariates can help. Third, simply adding covariates to a TWFE regression does not implement A1'. Caetano and Callaway (2026) point out that the transformation that removes the unit fixed effects also transforms the covariates, so the regression either quietly switches identification strategy or relies on a correct linear specification. They call the resulting bias hidden linearity bias.

### 3.6 From 2x2 to Staggered: Identifying $ATT(g,t)$

The 2x2 identification logic carries over directly to the staggered setting, and this is the key step to understanding all the modern estimators that follow. Index periods by $t = 1, \ldots, \mathcal{T}$, where $\mathcal{T}$ is the last period of the panel (in Kestrel $\mathcal{T} = 36$). For a cohort $g$ and any period $t$ after it is treated, replace the 2x2 treated group with cohort $g$ and the control group with valid units still untreated up to period $t$; the same double difference then identifies $ATT(g,t)$. Take Kestrel. To identify cohort 13's effect in month 20, the base period is month 12, and the control group can be the 8 never-treated cities, possibly joined by cohort 25, which has still not cut commissions in month 20. Cohort 19 was treated in month 19 and cannot be used. There are two ways to choose a valid control group, corresponding to two versions of parallel trends.

::: {.assumption}
**A4 (staggered adoption)** All units are untreated at $t = 1$, and once treatment switches on it never switches off: $D_{i1} = 0$, $D_{it} \geq D_{i,t-1}$. The cohort $G_i$ therefore fully describes unit $i$'s treatment path, and cohort $g$'s observed outcome in every period is $Y_t(g)$.
:::

::: {.assumption}
**A5 (parallel trends, based on never-treated)** For every cohort $g$ and all $t \geq g - \kappa$:

$$\mathbb{E}\big[Y_t(0) - Y_{t-1}(0) \mid G = g\big] = \mathbb{E}\big[Y_t(0) - Y_{t-1}(0) \mid C = 1\big]$$

where $C = \mathbf{1}\{G_i = \infty\}$ flags the never-treated.

**A5' (parallel trends, based on not-yet-treated)** For every cohort $g$ and all $t \geq g - \kappa$ and $s \geq t + \kappa$ (that is, every such pair $(s, t)$):

$$\mathbb{E}\big[Y_t(0) - Y_{t-1}(0) \mid G = g\big] = \mathbb{E}\big[Y_t(0) - Y_{t-1}(0) \mid D_s = 0,\, G \neq g\big]$$

The original paper also requires $s$ to precede the treatment period of the latest cohort, which guarantees that $\{D_s = 0, G \neq g\}$ is nonempty; this holds automatically when the population contains a never-treated group.
:::

Each version has a weak spot. A5 places no constraint on pre-treatment trends, but it requires the path of the 8 never-treated cities during the treatment period to represent the treated cities, and the cities a platform never touches may be marginal markets to begin with. A5' borrows cities that will eventually be selected, just scheduled later, and these are more comparable on the selection dimension. The price is that with $\kappa = 0$ it also requires the groups' pre-treatment trends to be parallel (Callaway and Sant'Anna 2021 §2.3), and the more the platform schedules by growth potential, the more suspect this condition becomes.

::: {.theorem}
**Theorem 2 (identification of $ATT(g,t)$; the covariate-free version of Callaway and Sant'Anna 2021 Theorem 1)** Suppose A2', A3, and A4 hold, the population contains a never-treated group, and cohort $g$ and the control group used both have positive population shares. Then for all $g \geq \kappa + 2$ (so that the base period $g - \kappa - 1$ lies inside the panel) and all $g - \kappa \leq t \leq \mathcal{T} - \kappa$:

(i) under A5, $ATT(g,t) = \mathbb{E}\big[Y_t - Y_{g-\kappa-1} \mid G = g\big] - \mathbb{E}\big[Y_t - Y_{g-\kappa-1} \mid C = 1\big]$;

(ii) under A5', $ATT(g,t) = \mathbb{E}\big[Y_t - Y_{g-\kappa-1} \mid G = g\big] - \mathbb{E}\big[Y_t - Y_{g-\kappa-1} \mid D_{t+\kappa} = 0,\, G \neq g\big]$.
:::

The proof sketch has four steps and the same structure as that of Theorem 1. Step one (A4): the observed outcome of cohort $g$ is $Y_t(g)$, so the treated group's change equals $ATT(g,t)$ plus $\mathbb{E}[Y_t(0) - Y_{g-\kappa-1}(g) \mid G = g]$. Step two (A2'): the base period $g - \kappa - 1$ comes before $g - \kappa$, so its expected outcome equals the expectation of $Y(0)$, and the expression becomes $ATT(g,t) + \mathbb{E}[Y_t(0) - Y_{g-\kappa-1}(0) \mid G = g]$. Step three (A5 or A5'): split this counterfactual increment into period-by-period increments from $g - \kappa$ to $t$, and replace each with the control group's increment. A5' holds for all $s \geq t + \kappa$, so the whole stretch can use the same control group $\{D_{t+\kappa} = 0, G \neq g\}$; this is why the assumption must hold for a whole string of values of $s$ at once. Step four (A3 and A2'): the observed outcome of the never-treated is $Y(0)$; the not-yet-treated group is still before its own anticipation window at period $t$, so by A2' its expected outcome equals the expectation of $Y(0)$. Positive shares guarantee that these conditional expectations are defined. Once population means are replaced by sample means, asymptotic normality under independent sampling of cities and the validity of the multiplier bootstrap are established in Callaway and Sant'Anna (2021) Theorem 2 and Theorem 3, which also give the regularity conditions for the case with covariates.

In words: $ATT(g,t)$ is the average effect of cohort $g$ in period $t$, with equal weight on every unit in the cohort. Under A5 the ruler is the never-treated; under A5' it is the units still untreated as of period $t + \kappa$. Neither ruler ever includes already-treated units. At this level there is no weighting problem yet. It appears when the whole block of $ATT(g,t)$ is collapsed into one number, and the weights for that step are chosen by the researcher (Section 4.5). In Kestrel, using the not-yet-treated ruler and aggregating by cohort gives 0.091, 0.086, and 0.059, against true values of 0.106, 0.085, and 0.058 (Section 6.5). Whichever ruler you use, one rule cannot be broken: never use already-treated units as the ruler. Static TWFE, dissected in Section 4.3, quietly breaks exactly this rule.

Theorem 2 says that anticipation does not block identification as long as it is bounded. What if $\kappa \geq 1$ but we still use $g - 1$ as the base period?

::: {.theorem}
**Corollary 1 (what a wrong base period estimates)** Under the conditions of Theorem 2 (i), if $\kappa \geq 1$ but $g - 1$ is still used as the base period, then for all $g \leq t \leq \mathcal{T} - \kappa$:

$$\mathbb{E}\big[Y_t - Y_{g-1} \mid G = g\big] - \mathbb{E}\big[Y_t - Y_{g-1} \mid C = 1\big] = ATT(g,t) - ATT(g, g-1)$$

where $ATT(g, g-1) = \mathbb{E}[Y_{g-1}(g) - Y_{g-1}(0) \mid G = g]$ is the anticipatory response in the base period itself.
:::

The proof sketch has three steps, matching those of Theorem 2. Step one (A4): as in Theorem 2. Step two (A2' does not reach the base period): the base period $g - 1$ lies inside the anticipation window, so $\mathbb{E}[Y_{g-1}(g) \mid G = g] = \mathbb{E}[Y_{g-1}(0) \mid G = g] + ATT(g, g-1)$, and there is no reason for this last term to be zero. Step three (A5 and A3): as in steps three and four of Theorem 2, the increment of $Y(0)$ from $g$ to $t$ is replaced by the never-treated increment and cancels, leaving only $ATT(g,t) - ATT(g, g-1)$.

In words: a wrong base period does not turn the result into noise; it systematically subtracts the anticipatory response of the base period. Holding back orders ($ATT(g, g-1) < 0$) leads to overstatement, and building volume ($ATT(g, g-1) > 0$) leads to understatement, consistent with the intuition of Section 3.3. The cost of the remedy is also clear: each time the base period moves one period earlier, A5 has to cover one more stretch of pre-treatment time, and Remark 1 of Callaway and Sant'Anna (2021) notes that the larger $\kappa$ is, the stronger the parallel trends assumption becomes. The Kestrel DGP has no anticipation ($\kappa = 0$), and in the CS event study the short period-by-period difference at $e = -1$ is -0.010 (SE 0.008), indistinguishable from zero. When this term is significantly different from zero in real data, the right response is to set $\kappa$ based on the mechanism and move the base period earlier, not to treat it as noise.

### 3.7 Triple Differences

Sometimes one DiD is not clean enough, but the difference between two DiDs is. Suppose Kestrel's commission cut applies only to small sellers whose monthly GMV is below some threshold. Comparing the change for small and large sellers within the treated cities removes the shocks common to the city (say, the local economy improving), but not a trend that small sellers share across the whole platform. The same comparison in the control cities measures exactly the "natural" trend gap between small and large sellers. Subtracting one from the other gives triple differences (DDD), whose regression form is in Wooldridge (2010) §6.5.2. The sample analog printed in that section has only three terms and omits the change for the ineligible group in the control cities; the complete four-term form is Theorem 3 below.

We follow Wooldridge's labels and write them in upright type: $\mathrm{B}$ and $\mathrm{A}$ denote treated and control cities, and $\mathrm{E}$ and $\mathrm{N}$ denote eligible and ineligible sellers. They serve only as cell names and subscripts; the upright $\mathrm{N}$ has nothing to do with the italic number of observations $N$ in Section 4.7. Only the $\mathrm{BE}$ cell is treated, and $\gamma_{xy} = \mathbb{E}[Y_{i2}(0) - Y_{i1}(0) \mid xy]$ is the counterfactual trend of cell $xy$, in the same sense as $\gamma(1)$ and $\gamma(0)$ in Theorem 1.

::: {.assumption}
**A6 (DDD parallel trends)** The counterfactual trend gap between the eligible and ineligible groups is the same in both types of city: $\gamma_{\mathrm{BE}} - \gamma_{\mathrm{BN}} = \gamma_{\mathrm{AE}} - \gamma_{\mathrm{AN}}$.
:::

A6 changes what has to be parallel. It no longer requires treated and control cities to be parallel, nor small and large sellers, only that the trend gap of "small sellers relative to large sellers" be the same in both types of city. It is most likely to fail when the treated cities simultaneously run another policy aimed only at small sellers, such as a concurrent traffic subsidy.

::: {.theorem}
**Theorem 3 (DDD identification)** Under A2, A3 (only the $\mathrm{BE}$ cell is treated, and the observed outcomes of the other three cells are $Y(0)$), and A6,

$$\Delta_{DDD} \equiv \big(\mathbb{E}[\Delta Y_i \mid \mathrm{BE}] - \mathbb{E}[\Delta Y_i \mid \mathrm{BN}]\big) - \big(\mathbb{E}[\Delta Y_i \mid \mathrm{AE}] - \mathbb{E}[\Delta Y_i \mid \mathrm{AN}]\big) = ATT_{\mathrm{BE}}$$
:::

The proof sketch has two steps. Step one (A2 and A3): as in steps one and two of Theorem 1, $\mathbb{E}[\Delta Y_i \mid \mathrm{BE}] = ATT_{\mathrm{BE}} + \gamma_{\mathrm{BE}}$, and for the other three cells $\mathbb{E}[\Delta Y_i \mid xy] = \gamma_{xy}$. Step two (A6): substituting gives $\Delta_{DDD} = ATT_{\mathrm{BE}} + (\gamma_{\mathrm{BE}} - \gamma_{\mathrm{BN}}) - (\gamma_{\mathrm{AE}} - \gamma_{\mathrm{AN}})$, and A6 makes the difference of the last two terms zero.

In words: DDD still estimates the treated cell's own ATT; the parallel trends requirement has merely moved from trends to differences in trends. Ortiz-Villavicencio and Sant'Anna (2025) raise two warnings. When conditioning on covariates is needed, neither "subtracting two DiDs" nor a three-way fixed effects regression generally identifies this ATT. In a staggered setting, pooling all not-yet-treated units into one control group introduces additional bias, which calls for a dedicated DDD estimator (the R package triplediff). The Kestrel simulation has no such structure, so no numbers are given here.

### 3.8 Testability and Pre-Trends

This section closes with testability, where misuse is most common. A1 is an assertion about the counterfactual, and the counterfactual is never observed, so parallel trends is fundamentally untestable within the DiD framework, just as exogeneity is untestable within the OLS framework. The usual pre-trends test checks something else: whether the two groups' observable outcomes were parallel before treatment actually happened. This is the observable shadow that parallel trends casts on the pre-period. It is a necessary condition only if you believe that parallelism spans both sides of the treatment date (A5, for instance, places no constraint on the pre-period). If the pre-period is already not parallel, that should of course make you nervous; but parallel pre-periods do not imply parallel counterfactuals after the treatment date.

This "do not imply" is not nitpicking. Imagine that Kestrel cuts commissions in a city because headquarters expects a strong competitor to enter that city in six months, and the cut is only one move in a package of responses that also includes subsidies and traffic tilting. Then, cut or no cut, this city and the control cities part ways from the moment the competitor enters, and none of this leaves any trace in the pre-period. Why did treatment happen here and now? Against this kind of worry the pre-trends test offers no defense: it can see only trend differences that have already happened, not those about to happen.

::: {.warning}
Flat pre-trends are supporting evidence, not proof of parallel trends, and low-power supporting evidence at that: a trend difference the test cannot detect may be large enough to flip the conclusion (Roth 2022; Section 7 has the details). Treating "the pre-trends test passed, so parallel trends holds" as a formal identification argument does not hold up.
:::

To sum up this section: Theorem 1 says that the 2x2 DiD identifies the treated group's ATT under A1 through A3, with each assumption blocking one route of bias. Theorem 2 extends the same double difference to $ATT(g,t)$ under staggered adoption; the ruler is the never-treated or the not-yet-treated, never already-treated units, and bounded anticipation is handled by moving the base period earlier. Corollary 1 says that a wrong base period estimates the effect minus the anticipatory response of the base period. Theorem 3 says that DDD replaces "parallel trends" with "parallel differences in trends." The conditional version A1' extends all of this to comparisons among cities of the same type. All of these assumptions concern the counterfactual, and none can be tested directly. What can be tested is only the shadow they cast on the pre-period, and between a parallel shadow and a parallel substance lies an inferential gap that data cannot fill.

## 4 Estimation: From 2x2 to Staggered

This section proceeds in order of failure: the method of each subsection works in the setting of the one before it, fails in a more complex setting, and so gives rise to the next method.

### 4.1 The 2x2 Case as a Regression

The sample analog estimator in the 2x2 case is the DiD of four means. It can be written equivalently as a regression:

$$Y_{it} = \alpha + \eta\, D_i + \lambda\, \mathrm{Post}_t + \beta\, (D_i \times \mathrm{Post}_t) + \varepsilon_{it}$$

The OLS estimate of the interaction coefficient, $\hat\beta$, exactly equals the sample DiD, $\hat\beta = \overline{\Delta Y}_T - \overline{\Delta Y}_C$ (Wooldridge 2010 §6.5.2). This is a numerical identity that does not depend on any distributional assumption, and its population version is Theorem 1. The regression form is convenient: adding covariates, computing clustered standard errors, and extending to many periods and groups all come easily. Card and Krueger (1994) implement exactly this equation plus restaurant characteristics. They also use a clever variant: the depth to which each restaurant is bound by the minimum wage law (the gap between its initial wage and the new minimum) defines a continuous treatment intensity, which turns the between-group comparison into a dose-response relationship (Section 5.1 has the details).

### 4.2 The Event Study Setup

Policy effects rarely arrive in one step. A commission cut first changes sellers' marginal profit, and only then passes through to stocking, advertising, and prices, so the effect accumulates over time. A two-period design cannot see any of this. With multi-period data the standard practice is to replace the static $\beta$ with a set of event-time coefficients, the event study (dynamic) specification:

$$Y_{it} = \alpha_i + \lambda_t + \sum_{\substack{e = -K \\ e \neq -1}}^{L} \theta_e\, \mathbf{1}\{t - G_i = e\} + \varepsilon_{it}$$

Here and below, $\theta$ follows the original notation of Callaway-Sant'Anna (Sun-Abraham write $\mu$) and denotes event-study coefficients and aggregated treatment effects; it is unrelated to the $\theta$ used elsewhere in the book for structural parameters. There are three syntactic details, each hiding a statistical reason. Base-period normalization: for units that are eventually treated, the full set of event-time dummies sums to 1 within the sample window and is perfectly collinear with the unit's fixed effect, so one dummy must be dropped. The convention is to drop $e = -1$ (the last pre-treatment period), so each $\theta_e$ reads as "the change relative to the period just before treatment," and $\theta_{-1} \equiv 0$ is a zero by definition, not an estimated zero.

Endpoint handling: leads and lags far from the treatment date are identified by fewer and fewer cohorts (only the earliest-treated cohort can contribute observations at $e = 20$). The standard practice is to bin the two ends ($\mathbf{1}\{t - G_i \leq -K\}$ and $\mathbf{1}\{t - G_i \geq L\}$ as one dummy each) or to trim them. Note that if you quietly omit some event times without binning them, those observations become implicit controls and the interpretation of every coefficient changes.

Without never-treated units there is a second collinearity among the event dummies, the time fixed effects, and the unit fixed effects, so one more coefficient must be normalized (a common choice is to also set the most distant lead to zero), and only relative dynamics are then identified: event time $e = t - G_i$ is a linear combination of calendar time and cohort, so a linear trend in event time is absorbed by the two sets of fixed effects.

::: {.warning}
The dynamic equation with $\sum_{e=-K}^{L}$ written out in full, omitting no base period, is perfectly collinear. The software will not throw an error. It silently drops a collinear column for you, which one depending on the column order, and your coefficients are then measured against a benchmark you do not even know. Always specify the reference (ref) period explicitly.
:::

When reading an event-study plot, check the following in order:

1. Whether the pre-period coefficients are close to zero, with no systematic climb or dip. This is supporting evidence, but its power is limited and it is not proof.
2. Whether the shape of the response at $e = 0$ matches the mechanism's timetable: the commission takes effect at the start of the month, so there should be a response in that same month; if the effect starts at $e = -2$, suspect anticipation.
3. The shape of the dynamic path: immediate, gradually accumulating, or overshooting and then falling back. The shape itself is information about the mechanism.
4. Whether the long run reaches a plateau or keeps drifting.
5. Whether the confidence intervals of the endpoint coefficients widen noticeably, and whether the cohorts identifying them have become too few to be representative.

### 4.3 How Static TWFE Fails Under Staggered Adoption

Now return to the number that did not add up at the start. With three batches of treatment dates, the static TWFE regression looks harmless, but it actually packages every 2x2 comparison the sample allows into one weighted average. There are three kinds of comparison: early-treated versus never-treated, late-treated versus never-treated, and early- and late-treated groups serving as each other's controls. The last kind hides the most dangerous comparison: the late-treated group's before-after change minus the contemporaneous change of the early-treated group, which is already treated. Here the early-treated group plays the "control group," yet it still carries a treatment effect. If the early group's effect has already settled at a constant, the differencing treats it as a level difference and cancels it, which is harmless. If the effect is still growing with event time (as it is in Kestrel), the growth is subtracted as "the control group's trend," and the late-treated group's estimate is pulled down. This kind of comparison is the so-called **forbidden comparison**: already-treated units mixed into the control group.

Goodman-Bacon (2021) turned this intuition into exact algebra.

::: {.theorem}
**Theorem 4 (Goodman-Bacon 2021 Theorem 1, the decomposition theorem)** Let timing group $k$ be treated first and $j$ later ($j > k$), and let $U$ be the never-treated group. Let $\mathrm{PRE}(k)$, $\mathrm{MID}(k,j)$, and $\mathrm{POST}(k)$ be the sample windows before $k$ is treated, between the two treatment dates, and after $k$ is treated, and let $\bar y$ be the group-by-window mean. Then the static TWFE estimator decomposes exactly as

$$\hat\beta_{TWFE} = \sum_{k \neq U} s_{kU}\, \hat\tau^{2\times2}_{kU} + \sum_{k \neq U} \sum_{j > k} \Big[ s^{k}_{kj}\, \hat\tau^{2\times2,k}_{kj} + s^{j}_{kj}\, \hat\tau^{2\times2,j}_{kj} \Big]$$

where the three types of 2x2 component are

$$\hat\tau^{2\times2}_{kU} = \big(\bar y_k^{\mathrm{POST}(k)} - \bar y_k^{\mathrm{PRE}(k)}\big) - \big(\bar y_U^{\mathrm{POST}(k)} - \bar y_U^{\mathrm{PRE}(k)}\big)$$

$$\hat\tau^{2\times2,k}_{kj} = \big(\bar y_k^{\mathrm{MID}(k,j)} - \bar y_k^{\mathrm{PRE}(k)}\big) - \big(\bar y_j^{\mathrm{MID}(k,j)} - \bar y_j^{\mathrm{PRE}(k)}\big)$$

$$\hat\tau^{2\times2,j}_{kj} = \big(\bar y_j^{\mathrm{POST}(j)} - \bar y_j^{\mathrm{MID}(k,j)}\big) - \big(\bar y_k^{\mathrm{POST}(j)} - \bar y_k^{\mathrm{MID}(k,j)}\big)$$

The weights satisfy $s \geq 0$ and $\sum s = 1$.
:::

The first type is the clean comparison of a treated group with the never-treated. The second uses the stretch in which the late-treated group is not yet treated as the control, and is also clean. The third type, $\hat\tau^{2\times2,j}_{kj}$, is the forbidden comparison: the post-treatment change of the late-treated group $j$ minus the contemporaneous change of the already-treated group $k$. The weights take the following form (write $n_k$ for a group's sample share, $n_{kj} = n_k / (n_k + n_j)$, $\bar D_k$ for the fraction of periods group $k$ spends in the treated state, and $\hat V^D$ for a normalizing constant):

$$s_{kU} = \frac{(n_k + n_U)^2\, n_{kU}(1 - n_{kU})\, \bar D_k (1 - \bar D_k)}{\hat V^{D}}$$

$$s^{k}_{kj} = \frac{\big((n_k + n_j)(1 - \bar D_j)\big)^2\, n_{kj}(1 - n_{kj})\, \frac{\bar D_k - \bar D_j}{1 - \bar D_j} \cdot \frac{1 - \bar D_k}{1 - \bar D_j}}{\hat V^{D}}$$

$$s^{j}_{kj} = \frac{\big((n_k + n_j)\bar D_k\big)^2\, n_{kj}(1 - n_{kj})\, \frac{\bar D_j}{\bar D_k} \cdot \frac{\bar D_k - \bar D_j}{\bar D_k}}{\hat V^{D}}$$

There is no need to memorize the formulas, but two of their implications matter. First, the weights are determined by group sizes and the variance of treatment status. The closer $\bar D_k$ is to 1/2 (the closer the treatment date is to the middle of the panel), the larger the weight: OLS allocates weight by "convenience for identification," not by "policy relevance." Second, the decomposition weights are all nonnegative, so the problem lies not in the weights but in what the forbidden comparisons estimate. Substituting potential outcomes into each 2x2 component and taking probability limits, one can rearrange the result into

$$\mathrm{plim}\; \hat\beta_{TWFE} = VWATT + VWCT - \Delta ATT$$

$VWATT$ is the variance-weighted average of the target effects of the 2x2 components; $VWCT$ is the variance-weighted violation of common trends, which is zero when parallel trends holds; and $\Delta ATT$ collects the within-group changes of the treatment effect over time, which enter with a negative sign through the forbidden comparisons. When the effect grows with event time, $\Delta ATT > 0$ and the static estimate is pulled down systematically, possibly even flipping its sign: all true effects can be positive while $\hat\beta_{TWFE} < 0$. In Kestrel it pulled 0.088 down to 0.023, and Section 6.3 redoes this accounting with numbers. The Goodman-Bacon decomposition is a diagnostic tool: it tells you what TWFE is actually estimating and where each piece of weight sits, but it does not provide a replacement estimator. That is the topic of Section 4.5.

::: {.warning}
Do not take "TWFE only mildly understates under staggered adoption" as a general rule. The direction and size of the bias are jointly determined by the effect dynamics and the weight structure. Increasing effects bias it downward (possibly to a sign flip), decreasing effects bias it upward, and instantaneous homogeneous effects leave it unbiased, because then every 2x2 component estimates the same constant without bias, and any nonnegative weights that sum to one recover the truth. Any intuition that "the TWFE result just needs a small correction" is wrong; you must do the decomposition or switch estimators outright.
:::

de Chaisemartin and D'Haultfœuille (2020) look at the same problem from the other side. Goodman-Bacon keeps the books by comparison: the weights are all nonnegative, and the problem hides in what the forbidden comparisons estimate. dCDH keep the books by cell, writing the coefficient directly as a weighted sum of the treated-cell effects $\tau_{it} = Y_{it}(G_i) - Y_{it}(0)$. Under parallel trends, $\beta_{TWFE} = \mathbb{E}\big[\sum_{(i,t):\, D_{it} = 1} w^{\mathrm{fe}}_{it}\, \tau_{it}\big]$ (their Theorem 1), where the weight $w^{\mathrm{fe}}_{it}$ is proportional to the residual from regressing $D_{it}$ on city and month fixed effects; the weights sum to one but can be negative. The residual is smallest in months when most cities are already treated and in the cities treated longest, so the negative weights fall on the late cells of the early-treated cohort (their Proposition 1).

In Kestrel, the cities of cohort 13 are treated for 24 of the 36 months, and from month 25 on, 52 of the 60 cities have cut commissions. The two sets of fixed effects then predict values of $D_{it}$ above 1 for these cells, so the residuals are negative. Of the 930 treated cells, 204 have negative weights (21.9%), and the negative weights sum to -0.1416. All of them fall in months 25 to 36 of cohort 13, that is, after the last batch of cities cut commissions, the stretch where cohort 13's effect is largest. The growth in cohort 13's late effects that the forbidden comparisons subtract in the Bacon decomposition corresponds here to those 204 negative weights: the two decompositions look at two sides of the same problem. Cell-level weight diagnostics take center stage again in the dose TWFE of Chapter 20.

### 4.4 The Contamination Problem of Dynamic TWFE

Since the static regression is flawed, many people's first reaction is "then surely a dynamic event study will do." Sun and Abraham (2021) prove that the dynamic version has a flaw of its own. The intuition: the dynamic regression assigns only one coefficient to each event time, while the two sets of fixed effects are fitted jointly on all cohorts. If one cohort's effect at period $e'$ differs from the others', the deviation is spread through the fixed effects into the coefficients of other event times, and not even the pre-period coefficients escape. Define the cohort-by-event-time effect

$$CATT_{g,e} = \mathbb{E}\big[Y_{i,g+e}(g) - Y_{i,g+e}(0) \mid G_i = g\big]$$

::: {.theorem}
**Theorem 5 (Sun and Abraham 2021 Propositions 1 to 3, the contamination result)** Consider the dynamic TWFE regression of Section 4.2 (with some relative periods excluded). Under parallel trends holding in all periods (stronger than A5) and no anticipation (A2' with $\kappa = 0$), the population regression coefficient on the dummy for event time $e$ satisfies

$$\theta_e = \sum_{e' \geq 0} \sum_{g} \omega_{e}^{g,e'}\, CATT_{g,e'}$$

The linear combination runs over all cohorts $g$ and all post-treatment relative periods $e' \geq 0$, including $e' \neq e$ and the excluded periods. The weights $\omega_e^{g,e'}$ are given in Proposition 1. They are determined solely by the design structure of treatment timing (and can be estimated from the data), and across all relative periods they satisfy the following: the own-period ($e' = e$) weights sum to 1 across cohorts; the weights of every other included period sum to 0 across cohorts; and the weights of all excluded periods sum to -1 in total. In the original Proposition 1, the never-treated cohort's weights are identically zero (its event-time dummies are all zero, so no assumption is needed), which is why the expression above omits it. Individual weights can be negative.
:::

This is a purely mechanical OLS result; it happens without any identification assumption failing. Three consequences are worth stating separately. If $CATT_{g,e}$ is homogeneous across cohorts and the excluded relative periods carry no effect themselves (for example, when only $e = -1$ is excluded), the cross terms cancel in the sum and $\theta_e$ is clean; otherwise the effects of the excluded periods still seep in, and homogeneity is no guarantee. Once effects are heterogeneous across cohorts, the coefficients for $e \geq 0$ mix in effects from other periods. More insidiously, the pre-period coefficients for $e < 0$ absorb other cohorts' post-period effects, so the pre coefficients can be significantly nonzero even when parallel trends and no anticipation both hold exactly. The pre-trends test then fails in both directions: it can raise an alarm when nothing is wrong (a false positive), or look flat because a contamination term offsets a real trend difference (a false negative). Because weights can be negative, $\theta_e$ can fall outside the convex hull of all the $CATT_{g,e}$, estimating an effect that no subgroup has.

In the Kestrel simulation this contamination is numerically mild (Section 6.2: the largest naive pre-period coefficient is only 0.009, and the post-period coefficients are only slightly below the truth), because the DGP's heterogeneity across cohorts is small and a never-treated group provides a backstop. This is a reminder to choose words carefully. Contamination is a theoretical risk determined by the design structure, and its severity varies with the data. One mild simulation does not show that it is always mild, and its theoretical existence does not license the claim that every TWFE event study is badly wrong.

### 4.5 Modern Estimators: CS, SA, BJS, and $DID_\ell$

The repairs share one idea: break the estimand down to the finest cells, $ATT(g,t)$ or $CATT_{g,e}$, use only clean comparisons within each cell, and then aggregate with visible weights that the researcher chooses. The four mainstream implementations differ in how the cells are estimated, which units serve as controls, and how strong the implied parallel trends assumption is.

Start with Callaway and Sant'Anna (2021). It takes the identification formula of Theorem 2 directly as its estimand: replace population means with sample means and estimate $\widehat{ATT}(g,t)$ cell by cell, using the never-treated as controls under A5 or the not-yet-treated under A5'. The anticipation argument of the did package is the $\kappa$ of A2' and Theorem 2; when it is set to a positive number, the base period moves automatically to $g - \kappa - 1$.

::: {.warning}
The CS control group can be the never-treated or the not-yet-treated (`control_group = "notyettreated"` in the did package); there is no restriction that "only the never-treated can be used." Accordingly, the right reading of the overlap condition is that each $(g,t)$ cell has a valid control group, not that every cohort must be paired with never-treated units. When the never-treated group is very small or suspected of selection (for example, the 8 cities the platform never touches may be marginal markets to begin with), the not-yet-treated version is often more credible.
:::

Once the $ATT(g,t)$ matrix is estimated, it has to be aggregated. Aggregating by cohort and then averaging overall:

$$\theta_{\mathrm{sel}}(g) = \frac{1}{\mathcal{T} - g + 1} \sum_{t = 2}^{\mathcal{T}} \mathbf{1}\{g \leq t\}\, ATT(g,t), \qquad \theta_{\mathrm{sel}}^{O} = \sum_{g = 2}^{\mathcal{T}} \theta_{\mathrm{sel}}(g)\, P(G = g \mid G \leq \mathcal{T})$$

Aggregating by event time into an event study:

$$\theta_{\mathrm{es}}(e) = \sum_{g = 2}^{\mathcal{T}} \mathbf{1}\{g + e \leq \mathcal{T}\}\, ATT(g, g + e)\, P(G = g \mid G + e \leq \mathcal{T})$$

A practical trap: in the did package, `aggte(type = "simple")` uses yet another weighting (by the number of post-treatment cells, so cohorts with longer exposure get more weight), and it is generally not equal to the $\theta_{\mathrm{sel}}^O$ that `type = "group"` reports. A paper must state clearly which one it reports. When there are only a few cohorts, the most honest approach is to report the $ATT(g,t)$ matrix itself and put the aggregate numbers after it.

The interaction-weighted (IW) estimator of Sun and Abraham (2021) implements the same philosophy as a regression, in three steps. Step one runs a saturated interaction regression that excludes $e = -1$ and uses the never-treated (or, if there are none, the latest-treated cohort, with its treated periods dropped from the sample) as the pure control:

$$Y_{it} = \alpha_i + \lambda_t + \sum_{g} \sum_{e \neq -1} \delta_{g,e}\, \mathbf{1}\{G_i = g\}\, \mathbf{1}\{t - g = e\} + \varepsilon_{it}$$

Because the cohort-by-event-time cells enter the regression fully saturated, the cross-cohort contamination of Section 4.4 cannot occur by construction, and each $\hat\delta_{g,e}$ estimates $CATT_{g,e}$ under parallel trends and no anticipation. Step two estimates the sample share of each cohort among those observed at event time $e$. Step three aggregates by these shares:

$$\hat\theta_e^{IW} = \sum_{g} \hat\delta_{g,e}\; \hat P\big(G = g \mid G \in \mathcal{G}_e\big), \qquad \mathcal{G}_e = \{g : g + e \in [1, \mathcal{T}]\}$$

The weights are cohort shares: naturally nonnegative, summing to one, and clearly interpretable, the three things the implicit weights of dynamic TWFE cannot deliver. With the never-treated as controls and $e = -1$ as the common base period, in a balanced panel without covariates, IW and the CS dynamic aggregation are numerically identical; only the way the standard errors are computed differs. fixest::sunab runs it in one line.

The imputation estimator of Borusyak, Jaravel and Spiess (2024) takes a different angle: rather than look for a control group to difference against, fill in the counterfactual directly. Assume for the untreated state $Y_{it}(0) = \alpha_i + \lambda_t + \varepsilon_{it}$, then proceed in three steps. Using only untreated observations (all periods of the never-treated cities plus the pre-periods of the treated cities), estimate $\hat\alpha_i, \hat\lambda_t$ by OLS. Impute the counterfactual $\hat Y_{it}(0) = \hat\alpha_i + \hat\lambda_t$ for each post-treatment observation, which gives the cell-level effect $\hat\tau_{it} = Y_{it} - \hat Y_{it}(0)$. Finally, average with weights chosen by the researcher, $\hat\tau_w = \sum_{(i,t):\, D_{it} = 1} w_{it}\, \hat\tau_{it}$; an event study, a cohort average, or an overall average is just a different choice of weights. Here the $w_{it}$ are nonnegative and chosen by the researcher in advance, in direct contrast with the $w^{\mathrm{fe}}_{it}$ implied by TWFE in Section 4.3.

Four of its properties must be stated plainly. The fixed effects are fitted on all pre-periods, so the parallel trends it implies covers the entire pre-period, not just the local comparison from base period to target period that CS uses. This is a strictly stronger assumption, and what it buys is efficiency: BJS prove that under homoskedastic spherical errors it is the most efficient linear unbiased estimator. Treated observations never enter the fitting regression, so the forbidden comparison does not exist by construction. The pre-trends test must be done separately (re-estimate on the untreated observations with lead dummies added and test that their coefficients are zero), because the estimation procedure mechanically flattens the pre-period fit, and looking directly at a residual plot would fool you. Cell-level effect heterogeneity is folded into the error term, so inference is conservative for the sample ATT. The R implementation is didimputation::did_imputation; the fect package of Liu, Wang and Xu (2024) extends the same idea to interactive fixed effects, and the other members of the family are in Section 4.8.

The fourth is the $DID_\ell$ of de Chaisemartin and D'Haultfœuille (2026) (R package DIDmultiplegtDYN). Here $\ell$ indexes the horizon of the effect in periods; it is not the log-likelihood $\ell(\theta)$ of Part II. For each unit first treated in period $g$, the estimator takes the change in the outcome from $g - 1$ to $g - 1 + \ell$, subtracts the average change over the same stretch among units whose treatment is still unchanged as of period $g - 1 + \ell$ (including the never-treated), and then averages with equal weights over all treated units that can be followed through their $\ell$-th period. So $DID_1$ is the effect in the month of treatment, and $\ell = e + 1$.

Its identifying assumptions are no anticipation plus period-by-period parallel trends in $Y(0)$ for all groups (Assumptions 1 and 2 of the original paper; throughout this chapter, numbering for that paper follows its arXiv working-paper version v14). This is slightly stronger than A5', because it requires parallel trends in all periods for every group, including the periods before a cohort's own treatment, and the accompanying placebo estimators exist precisely to test this. Under these two assumptions $DID_\ell$ is conditionally unbiased (Lemma 1 of the original paper) and asymptotically normal as the number of groups grows (Theorem 1 of the original paper).

In a design like this chapter's (binary, absorbing, every unit starting untreated), $DID_\ell$ is numerically identical to the CS estimate $\hat\theta_{\mathrm{es}}(\ell - 1)$ with not-yet-treated controls, an equivalence that de Chaisemartin and D'Haultfœuille (2026) §3.1 states explicitly. There are two reasons. The comparison sets are the same: for each cohort, period $g - 1$ against period $g - 1 + \ell$, minus the units still untreated at that time. The aggregation weights are the same: $DID_\ell$ weights every treated unit equally, which is exactly CS weighting by cohort size. In Kestrel the two point estimates differ only by floating-point error of about $4 \times 10^{-17}$; only the standard errors differ (analytic cluster standard errors versus the multiplier bootstrap). Since the numbers coincide, in this design the effect estimates of $DID_\ell$ in fact rely only on A5'; the paper's stronger all-period assumption serves its placebo estimators, which test the pre-treatment periods, and non-absorbing, non-binary designs. The distinctive strength of $DID_\ell$ lies elsewhere: treatment can be non-binary and can switch on and off.

The table below compares the four main estimators:

| Dimension | Callaway-Sant'Anna | Sun-Abraham IW | BJS imputation | dCDH $DID_\ell$ |
|---|---|---|---|---|
| Strength of parallel trends | Parallel trends required only for periods $t \geq g - \kappa$, the weakest of the four (never-treated version; the not-yet-treated version also needs parallel pre-treatment trends for later cohorts) | Parallel trends in all periods, slightly stronger than CS | Parallel trends in all periods, with every pre-period entering the fit | Parallel trends in all periods, slightly stronger than CS |
| Control group | Never-treated or not-yet-treated, the user's choice | Never-treated, or the latest-treated cohort | All observations untreated at the time (both types used) | Units whose treatment is unchanged through the outcome period (never-treated included) |
| Covariates | Outcome regression / IPW / doubly robust | Enter the fixed-effects structure linearly | Enter the imputation regression linearly | Linear residualization on first differences of the covariates |
| Anticipation | Built-in anticipation argument that moves the base period | Partly handled by moving the base period manually | Suspect periods must be removed explicitly from the fitting sample | Assumed away, checked with placebos |
| Efficiency | Robustness first, efficiency second | Robustness first | Most efficient under homoskedasticity | Robustness first |
| Under violations of parallel trends | Less sensitive (local comparisons only) | Less sensitive | Most sensitive (relies on extrapolation from a long pre-period) | Less sensitive (local comparisons only) |
| R implementation | did::att_gt plus aggte | fixest::sunab | didimputation::did_imputation | DIDmultiplegtDYN::did_multiplegt_dyn |

What does DiD identify when treatment can be exited, can take several values, or is actually taken up by only some units? That is the subject of Chapter 20.

### 4.6 Estimation with Covariates

Under A1', the three implementation routes differ in which nuisance model you trust. The outcome regression route models the control group's outcome growth, $m_{g,t}(X) = \mathbb{E}[Y_t - Y_{g-1} \mid X, C = 1]$ (here $\kappa = 0$), and uses it to predict the treated group's counterfactual growth. The IPW route models "being selected into cohort $g$," using the generalized propensity score $p_g(X) = P(G_g = 1 \mid X, G_g + C = 1)$ to reweight the control group to the treated group's $X$ distribution. The doubly robust route combines the two. In the version with never-treated controls ($G_g = \mathbf{1}\{G_i = g\}$), the doubly robust estimand of Callaway and Sant'Anna (2021), which extends Sant'Anna and Zhao (2020), is

$$ATT(g,t) = \mathbb{E}\left[\left(\frac{G_g}{\mathbb{E}[G_g]} - \frac{\dfrac{p_g(X)\, C}{1 - p_g(X)}}{\mathbb{E}\left[\dfrac{p_g(X)\, C}{1 - p_g(X)}\right]}\right)\big(Y_t - Y_{g-1} - m_{g,t}(X)\big)\right]$$

As long as one of the two models, $p_g$ or $m_{g,t}$, is correctly specified, the estimate is consistent. Setting $m_{g,t}$ to zero reduces it to pure IPW, and removing the weighting term reduces it to pure outcome regression. The did package's default, `est_method = "dr"`, applies this formula cell by cell to every $(g,t)$; the not-yet-treated version merely swaps the comparison set. Two rules: covariates may only take pre-treatment values; and conditioning changes the credibility of the parallel trends assumption, not the algebra of aggregation, so do not expect a few more control variables to repair a design-level selection problem.

### 4.7 Inference

That covers the point-estimation problem of DiD, but the standard-error problem can ruin a paper just as surely. The logic of clustering fits in one sentence: which observations do you believe are independent of each other? A cluster structure assumes independence across clusters and arbitrary correlation within them, and the asymptotics run on the number of clusters $N_{\text{cl}}$, not the number of observations $N$ (Chapter 8 writes the number of clusters as $G$; in this chapter $G$ is reserved for the cohort). The CR1 sandwich estimator is

$$\hat V^{CR1}_{\hat\beta} = \frac{N_{\text{cl}}}{N_{\text{cl}}-1}\,\frac{N-1}{N-\dim(\beta)}\,(X'X)^{-1}\left(\sum_{c=1}^{N_{\text{cl}}} X_c' \hat u_c \hat u_c' X_c\right)(X'X)^{-1}$$

where $\dim(\beta)$ is the number of regression coefficients. The leading scalar is a finite-sample correction; dropping it gives the uncorrected CR0 sandwich, and with few clusters the difference between the two is not negligible.

At which level should you cluster? At the level of treatment assignment. Kestrel's commission policy is switched on city by city, so the city is the independent unit of the design, and standard errors are clustered by city. Clustering at the seller or city-month level would severely understate the standard errors, because it pretends that shocks in the same city in different months are independent, which the AR(1) noise plainly contradicts.

Fixed effects can absorb only the time-invariant part of within-cluster correlation (city fixed effects absorb the part that is constant at the city level). Serial correlation that decays over time, such as AR(1), stays in the residuals, so you still cluster by city after including city fixed effects. The rule of thumb in Chapter 8 is that if the model includes fixed effects at some level, you usually need to cluster at least at that level. It is a conservative rule; the level is ultimately determined by the design of sampling and treatment assignment. Kestrel also includes month fixed effects, but they absorb the shocks common to all cities in a given month, leaving residuals correlated only across months within a city, so we cluster only at the city level. The wrong level costs more than a few percent: standard errors can be off by an order of magnitude.

The classic source of this discipline is Bertrand, Duflo and Mullainathan (2004). They observed that typical DiD data (state-by-year wage panels) are strongly serially correlated, that the treatment variable is itself serially correlated (once switched on, it stays on), and that papers at the time generally did not cluster. Their diagnostic was elegant and brutal. In CPS wage data they assigned placebo laws at random (random states, random effective years, a true effect of exactly zero) and recorded how often the unclustered test rejected the null. The answer: at the 5% nominal level, up to 45% of the placebos "found an effect." They also verified the remedies that work: cluster by state (the level of treatment assignment), use a block bootstrap by state, or collapse the time series into one pre mean and one post mean per state and run a 2x2. The last remedy also yields a nice intuition: the effective number of observations in DiD is the number of clusters, not the number of panel cells.

With few clusters (the rule-of-thumb threshold is roughly $N_{\text{cl}} < 30$ to 40, and things are worse when there are few treated clusters), CR1 with $t_{N_{\text{cl}}-1}$ critical values over-rejects. The remedies escalate with severity: CR2 (the Bell and McCaffrey 2002 small-sample correction plus Satterthwaite degrees of freedom); the wild cluster bootstrap-t (Cameron, Gelbach and Miller 2008, R package fwildclusterboot); and when the treated clusters number in the single digits, even the wild bootstrap is unreliable, and the honest fallback is randomization inference along the lines of Conley and Taber (2011) and MacKinnon and Webb (2020). Kestrel has 60 cities and three cohorts, so the baseline CR1 suffices; if your application has only a dozen or so cities, CR2 and the wild bootstrap belong in the main table, not the appendix.

### 4.8 The Estimator Family: Clean Comparisons plus Explicit Weights

Since 2021 the estimators have multiplied, but most are variations on one recipe: estimate cell-level or cohort-level effects using only clean comparisons, then aggregate with visible weights. The table below puts four other common implementations side by side; the last column gives their numerical relation, on Kestrel, to the estimators discussed earlier in this section (the numbers come from the actual runs in Section 6.6):

| Estimator | Where the clean comparisons come from | Weights | Numerical relation on Kestrel |
|---|---|---|---|
| Wooldridge ETWFE (etwfe) | All cohort-by-period interactions in TWFE, equivalent to a two-way Mundlak regression | Equal weight on each treated cell | Not-yet-treated controls: 0.078, within about $2 \times 10^{-8}$ of BJS; never-treated controls: 0.081, within about $2 \times 10^{-15}$ of Sun-Abraham |
| Gardner (2022) two-stage (did2s) | The first stage estimates the two sets of fixed effects on untreated observations only; the second regresses the outcome net of the fixed effects | Equal weight on each treated cell | 0.078, within about $7 \times 10^{-17}$ of BJS |
| Stacked DiD (Wing, Freedman and Hollingsworth 2024) | Each cohort is paired with clean controls inside its event window to form a sub-experiment; the stacked data are then regressed | With corrective sample weights, identifies a trimmed aggregate ATT | Not estimated |
| LP-DiD (Dube, Girardi, Jordà and Taylor 2025) | One long-difference regression per horizon, keeping only newly treated units and clean controls still untreated at that horizon | Variance-weighted and nonnegative by default; after reweighting, equal weights, the same as CS | Not estimated |

What the four share is that the estimand is first broken down to the level at which clean comparisons can identify it, and the weights are written out in the open. Wooldridge (2025) proves that TWFE and the two-way Mundlak regression are numerically identical: if the unit mean of each regressor $X_{it}$, $\bar X_i = \mathcal{T}^{-1}\sum_{t=1}^{\mathcal{T}} X_{it}$, and its cross-sectional mean for each period are added to pooled OLS, then in a balanced panel the coefficients equal those of TWFE (see Theorem 5 (iii) in Section 3.5 of Chapter 10 and Section 4.3 there; the one-way version is the Mundlak regression of Wooldridge 2010 §10.7.3). ETWFE, which lets the treatment effect vary freely by cohort and period, is in turn equivalent to the imputation estimator. So the flaw of static TWFE lies not in the fixed-effects regression itself but in forcing the effects of all treated cells into a single $\beta$. The two-stage method of Gardner (2022) and BJS are two ways of writing the same estimator. This chapter does not run stacked DiD or LP-DiD; the weighting properties in the table come from the original papers. Which one to choose depends mainly on software, covariates, and convenience of inference, not on which is "more correct."

The roadmap of this section is now complete: the 2x2 regression is equivalent to DiD; with many periods, the event study expresses the dynamics; under staggered adoption, static TWFE is exposed by Goodman-Bacon and dCDH from the comparison side and the cell side, and dynamic TWFE by Sun-Abraham; CS, SA-IW, BJS, and $DID_\ell$ repair the problem with clean comparisons plus transparent weights, each paying with parallel trends of a different strength, while ETWFE, two-stage, stacked, and LP-DiD are other ways of writing the same recipe; covariates are handled by doubly robust estimation; and inference clusters at the level of treatment assignment, with a separate set of remedies when clusters are few.

## 5 Anchoring Papers

Methods prove themselves only in real research. This section anchors the chapter in two papers: one that is the historical origin of the 2x2 method and one that is a staggered application in IS. Each is laid out along the five elements of paper, method, data, results, and limitations, with a focus on how the assumptions are defended.

### 5.1 Card and Krueger (1994, AER)

::: {.case}
The paper and its place in the history of method: "Minimum Wages and Employment: A Case Study of the Fast-Food Industry in New Jersey and Pennsylvania," American Economic Review. This paper turned DiD from a trick into the default design language of policy evaluation, and it is one of the emblematic texts of what was later called the credibility revolution: the identification argument shifted from "my model specification is right" to "where does my policy variation come from."

Method: a textbook 2x2. On April 1, 1992, New Jersey raised the state minimum wage from \$4.25 to \$5.05, while eastern Pennsylvania next door kept its minimum unchanged. The treated group is NJ fast-food restaurants, the control group is eastern PA fast-food restaurants, and employment is measured in one wave before and one after. Beyond the between-group comparison, the paper exploits variation in treatment intensity within NJ: restaurants whose starting wage was already above \$5.05 were effectively unconstrained, so the depth of the constraint (the wage gap) forms a continuous treatment variable, and the dose-response relationship provides independent corroboration of the between-group comparison.

Data: a telephone survey of about 410 fast-food chain restaurants (Burger King, KFC, Wendy's, Roy Rogers); the analysis sample with usable data in both waves is slightly smaller. The first wave ran in February and March 1992, before the increase, and the second in November and December, after it. The choice of the fast-food industry is itself a design argument: the minimum wage binds hardest here, and the jobs are homogeneous and easy to measure.

Results: contrary to the textbook prediction for a competitive labor market, employment in NJ did not fall relative to PA. The DiD point estimate is positive (full-time-equivalent employment per restaurant rose by 2.76 jobs in relative terms, about 13%, $t = 2.03$), and the continuous wage-gap version gives results in the same direction.

Limitations and aftermath: the telephone survey measured employment with a lot of noise. Neumark and Wascher (2000) re-estimated with restaurant payroll data and found the opposite sign, which set off a head-on exchange in the AER in 2000. Card and Krueger (2000) re-examined the question with administrative BLS payroll data and kept their conclusion; their working-paper version (Card and Krueger 1998) says the increase "probably had no effect on total employment ... possibly had a small positive effect." This kind of self-correction, narrowing a claim in the face of evidence, is worth learning from. Seen with today's eyes, the design-level limitations are clearer. There is essentially one treated state and one control state, so inference faces the extreme few-cluster problem of $N_{\text{cl}} = 2$ (none of the tools of Section 4.7 existed then). Parallel trends is untestable, and differences between the two states' macro cycles around 1992 are a real threat. Composition is another hidden risk: cross-border flows, such as NJ restaurants closing while PA restaurants open, would contaminate both groups at once.
:::

The paper's defense strategy is worth emulating: a neighboring state and the same industry hold down common shocks, a tight before-after window around the policy date holds down long-run trends, and continuous treatment intensity within the state provides an internal cross-check. Every step buys insurance for the credibility of the counterfactual, rather than piling control variables into the regression.

### 5.2 Greenwood and Wattal (2017, MISQ)

::: {.case}
The paper: "Show Me the Way to Go Home: An Empirical Investigation of Ride-Sharing and Alcohol Related Motor Vehicle Fatalities," MIS Quarterly. The question is whether the entry of ride-sharing reduces drunk-driving fatalities, and the paper is a representative staggered DiD study in IS. The data, magnitudes, and robustness details below follow the 2015 working-paper version of the same study (Greenwood and Wattal 2015).

Method: staggered DiD. Uber entered California cities in batches, the staggered entry dates provide the identifying variation, and the regression includes location and time fixed effects. Note that the paper was published in 2017, before the estimator revolution of 2021, and uses the TWFE-type specification that was standard at the time.

Data: a panel of California locations, with traffic fatalities from the statewide accident records system of the California Highway Patrol (CHP/SWITRS), matched to the entry timetable of each Uber product line.

Results: after the low-price product UberX enters, alcohol-related motor vehicle fatalities fall significantly, with the fatality rate declining by 3.6% to 5.6% per quarter; the high-price product Uber Black has no robust effect. The split along the price dimension is itself evidence on mechanism: if the decline were only an artifact of trends or selection, there would be no reason for it to appear only for the low-price product, so the product-line contrast doubles as a placebo.

Limitations: Uber's entry dates are not random. The decision to enter a city is related to city size, the regulatory environment, and demand growth. The authors control for the local characteristics that drive entry and check pre-trends, confronting the issue rather than dodging it. By today's standards, the staggered TWFE estimates should be redone with CS or SA, and effect heterogeneity (large cities are entered first and may have larger effects) is precisely a real-world version of the negative-weight risk of Section 4.3. Spillovers are also real: substitute travel in neighboring areas and drunk drivers crossing boundaries both erode the cleanliness of the control group.
:::

Taken together, the two papers make the point of anchoring clear. Card and Krueger show how a 2x2 design stands on institutional detail and internal cross-checks. Greenwood and Wattal show how staggered variation is used in platform research and how endogenous entry is defended, and also how methodological progress can make an old paper's estimation strategy obsolete while its design ideas stay current.

## 6 A Full Walkthrough on Kestrel Data

Now run the full toolkit of Section 4 on the Kestrel panel from start to finish. The code below uses R 4.6.1 with set.seed(11) fixed for reproducibility, and every number quoted in the text comes from the actual output of this code.

### 6.1 The DGP

The design parameters: 60 cities over 36 months, 2160 city-month observations in total; three cohorts cut commissions in months 13, 19, and 25, with 17, 17, and 18 cities, and 8 cities are never treated.

```r
set.seed(11)

T_months     <- 36
cohort_g     <- c(13, 19, 25)   # first treated month per cohort
cohort_size  <- c(17, 17, 18)   # cities per cohort
n_never      <- 8               # never-treated cities
n_cities     <- sum(cohort_size) + n_never   # 60

tau0         <- 0.015
tau5         <- 0.020
plateau      <- 0.140
slope_late   <- 0.004
cohort_boost <- c(0.000, 0.004, 0.008)
mult_sd      <- 0.10

tau_profile <- function(e) {
  ifelse(e <= 5, tau0 + (tau5 - tau0) * e / 5,
  ifelse(e <= 12, tau5 + (plateau - tau5) * (e - 5) / 7,
         plateau + slope_late * (e - 12)))
}
```

The effect is a piecewise-linear ramp in event time: a slow start (from 0.015 to 0.020 over $e = 0$ to 5), then rapid growth to a plateau of 0.140 at $e = 12$, then a further drift of 0.004 per month. The cohort boost gives later cohorts a slightly larger effect at the same event time, and a city-level multiplier creates heterogeneity across cities. Assembling the outcome and recording the true effect:

```r
dt[, event_time := ifelse(g > 0, month - g, NA_integer_)]
dt[, D := as.integer(g > 0 & month >= g)]
dt[, tau_true := 0]
dt[D == 1, tau_true := (tau_profile(event_time) +
                          cohort_boost[match(g, cohort_g)]) * tau_mult]
dt[, lgmv := 8 + city_fe + month_fe + tau_true + eps]
```

The city fixed effects carry level selection (early-treated cities sit higher), the month fixed effects contain linear growth plus a sinusoidal seasonal term, and the noise is AR(1) within city ($\rho = 0.35$). Parallel trends holds exactly by construction, with no anticipation and no spillover. Whether each estimator below succeeds or fails therefore depends only on aggregation algebra, not on identification, which is exactly the purpose of this simulation design.

The true values computed from realized effects: the simple ATT is 0.088; by cohort, 0.106 for g13, 0.085 for g19, and 0.058 for g25. Note the direction. Later cohorts have larger effects at the same event time (the boost is positive), but their exposure is short and misses the second half of the ramp, so their realized cohort ATT is smaller. The answer to "whose effect is larger" depends on how you ask, which is one more reason to fix the estimand before discussing estimation.

One design decision needs explaining. The never-treated group has only 8 cities, deliberately kept small. The weights of static TWFE depend only on the treatment matrix. With 20 never-treated cities, the clean treated-versus-untreated comparisons would take most of the weight, $\hat\beta_{TWFE}$ would be much closer to the truth, and the contrast would be less striking; still, as long as some weights are negative, a steep enough effect curve would pull it far down anyway. Shrinking the clean control group to 8 cities pushes about 61% of the weight onto the timing comparisons, and the teaching point becomes obvious at once. The algebraic fact that "the never-treated share drives the TWFE bias" deserves attention in itself: the same estimator and the same effect structure, with a different control-group size, give a very different bias.

### 6.2 Static TWFE and the Naive Event Study

```r
m_static <- feols(lgmv ~ D | city_id + month, data = dt, cluster = ~city_id)
```

The output is $\hat\beta_{TWFE} = 0.0227$ (SE 0.0087), only 25.8% of the truth of 0.088. This is the number from Section 1. The naive event study, the dynamic version, takes just one more step:

```r
dt[, ev_naive := fifelse(g > 0, pmin(pmax(month - g, -12L), 18L), -1L)]
m_naive <- feols(lgmv ~ i(ev_naive, ref = -1) | city_id + month,
                 data = dt, cluster = ~city_id)
```

The base period is $e = -1$, the endpoints are binned at $e \leq -12$ and $e \geq 18$, and never-treated cities are coded $e = -1$ so that all their dummies are zero and they fall into the reference group. The result is surprisingly respectable. The largest pre-period coefficient is only 0.009 ($e = -3$, SE 0.008), and none is significant. The post-period coefficients broadly track the truth, only slightly low across the board (the reconciliation table in Section 6.8 shows this systematic mild understatement).

This confirms the caution at the end of Section 4.4: the severity of Sun-Abraham contamination depends on how heterogeneous the cohorts are and on the control structure. In this DGP the cohort boost is small and a never-treated group provides a backstop, so the contamination is numerically mild. Do not conclude from this that dynamic TWFE is harmless, and do not, conversely, claim in your own paper that it must be badly wrong. The right move is to run SA or CS alongside it and see how much they differ.

That the static regression is badly wrong while the dynamic regression is broadly fine shows that the problem lies in the weight structure of the static coefficient. The next section takes it apart with the Bacon decomposition.

### 6.3 Goodman-Bacon Decomposition: Where This 0.023 Comes From

```r
df <- as.data.frame(dt[, .(city_id, month, lgmv, D)])
bd <- bacon(lgmv ~ D, data = df, id_var = "city_id", time_var = "month")
```

The decomposition produces nine 2x2 components whose weighted sum reproduces 0.023 exactly. Summarized by comparison type:

| Comparison type | Weight | Average 2x2 estimate | Contribution |
|---|---|---|---|
| Treated vs Untreated | 0.390 | 0.073 | 0.029 |
| Earlier vs Later Treated | 0.306 | 0.026 | 0.008 |
| Later vs Earlier Treated (forbidden) | 0.304 | -0.046 | -0.014 |

Each of the three rows tells a story. The clean treated-versus-never-treated comparisons average 0.073, and each of these components uses the full post window of its cohort. The average falls below 0.088 for two reasons: the Bacon weights give g25, the cohort with the smallest effect, more weight than the simple ATT does, and there is sampling noise. These comparisons receive only 39% of the weight. The earlier-versus-later comparisons (using the later cohort's pre-period as the control) are also clean; their average of 0.026 is low because they can see only the early part of the earlier cohort's ramp.

The real damage comes from the third row: 30.4% of the weight sits on forbidden comparisons, whose average estimate is negative, -0.046. The individual components make this more concrete: 25 versus 13 has weight 0.138 and estimate -0.054; 19 versus 13 has weight 0.098 and estimate -0.036; 25 versus 19 has weight 0.069 and estimate -0.043. The mechanism is exactly the algebra of Section 4.3. The earlier cohort's effect is still in the early part of its ramp during the MID window ($e = 0$ to 11 for 25 versus 13), and by the window after the later cohort is treated it has climbed a great deal. This growth is subtracted from the later cohort's change as "the control group's trend," producing a string of negative values. All true effects are strictly positive, yet the components carrying three tenths of the weight are estimated as negative, and that is how 0.088 is diluted to 0.023.

### 6.4 Sun-Abraham and the Event-Study Plot

```r
dt[, g_sunab := fifelse(g == 0, 10000L, as.integer(g))]
m_sa <- feols(lgmv ~ sunab(g_sunab, month) | city_id + month,
              data = dt, cluster = ~city_id)
sa_att <- summary(m_sa, agg = "att")
```

sunab requires the never-treated cohort value to fall outside the sample period, so we first recode 0 as 10000. The aggregated ATT is 0.081 (SE 0.008). fixest's built-in iplot(m_sa) produces a diagnostic plot in one line; the comparison plot below is drawn with ggplot2 so that the true-effect curve can be overlaid:

![Naive TWFE event study, Sun-Abraham, and the true dynamic effect compared, $e \in [-12, 18]$, with 95% confidence intervals.](assets/fig/fig_19_eventstudy.svg)

### 6.5 Callaway-Sant'Anna

```r
set.seed(11)
cs <- att_gt(yname = "lgmv", tname = "month", idname = "city_id",
             gname = "g", data = dt, control_group = "notyettreated",
             clustervars = "city_id", bstrap = TRUE, biters = 1000)

agg_simple  <- aggte(cs, type = "simple")
agg_dynamic <- aggte(cs, type = "dynamic")
agg_group   <- aggte(cs, type = "group")
```

The control group is the not-yet-treated, with doubly robust estimation by default. Standard errors come from the multiplier bootstrap, and fixing the seed before estimation is what makes the results reproduce digit for digit. The simple aggregation gives 0.082 (SE 0.007), with 95% confidence interval [0.068, 0.095], which covers the truth of 0.088. The overall estimate of the group aggregation is 0.078 (SE 0.007), and by cohort:

| cohort | CS group ATT (SE) | True value |
|---|---|---|
| g = 13 | 0.091 (0.010) | 0.106 |
| g = 19 | 0.086 (0.009) | 0.085 |
| g = 25 | 0.059 (0.012) | 0.058 |

g19 and g25 land almost exactly on the truth, while g13 is somewhat low (its long-horizon cells can be identified only by the tiny never-treated group, so they are noisy). The gap between the simple and group aggregations (0.082 versus 0.078) is also telling. Simple weights by the number of post-treatment cells, so g13, with the longest exposure, gets more weight. The two numbers answer different questions and both are correct, but a paper must state which one it uses.

### 6.6 BJS Imputation and Its Relatives

```r
imp <- did_imputation(data = dt, yname = "lgmv", gname = "g",
                      tname = "month", idname = "city_id")

imp_es <- did_imputation(data = dt, yname = "lgmv", gname = "g",
                         tname = "month", idname = "city_id",
                         horizon = TRUE, pretrends = -12:-2)
```

The overall ATT is 0.078 (SE 0.005), the smallest standard error in the grand reconciliation table of Section 6.8, consistent in direction with the efficiency result even though the AR(1) errors of this DGP are not spherical. The point estimate is about two standard errors below the realized truth, and the $e = 0$ coefficient (0.005) is also low. This is the noise this seed happened to draw, not a systematic bias. The lesson: when writing a paper, report estimates with their standard errors, and do not claim that any modern estimator "exactly equals" the truth.

The two relatives from Section 4.8 take only a few lines each on the same data:

```r
df_full <- as.data.frame(dt)
fit_d2s <- did2s(df_full, yname = "lgmv",
                 first_stage  = ~ 0 | city_id + month,
                 second_stage = ~ i(D, ref = 0),
                 treatment = "D", cluster_var = "city_id")

etw <- etwfe(fml = lgmv ~ 0, tvar = month, gvar = g, data = df_full,
             ivar = city_id, cgroup = "notyet", vcov = ~city_id)
emfx(etw, type = "simple")
```

The point estimates of did2s and etwfe are both 0.078, the same as imputation; switching etwfe's control group to `cgroup = "never"` gives 0.081 (SE 0.008), the same as Sun-Abraham. Yet the three cluster standard errors for the same estimand all differ (imputation 0.0047, etwfe 0.0048, did2s 0.0053). The differences come from how the variance estimators are constructed and have nothing to do with the point estimates.

### 6.7 dCDH: TWFE Weights and $DID_\ell$

```r
library(TwoWayFEWeights)
library(polars)            # must be attached before DIDmultiplegtDYN
library(DIDmultiplegtDYN)

tw <- twowayfeweights(df, Y = "lgmv", G = "city_id", T = "month", D = "D",
                      type = "feTR", summary_measures = TRUE)

dyn <- did_multiplegt_dyn(df, outcome = "lgmv", group = "city_id", time = "month",
                          treatment = "D", effects = 12, placebo = 6,
                          cluster = "city_id", graph_off = TRUE)
```

Here df is the four-column data of Section 6.3. The weight diagnostic reproduces the numbers of Section 4.3: 930 treated cells, 726 with positive weights and 204 with negative ones; the positive weights sum to 1.1416 and the negative weights to -0.1416, and all the negative weights fall in months 25 to 36 of cohort 13. The total weights of the three batches of cities are nearly equal (0.325, 0.339, 0.336), so the problem is not across cohorts but across event times: the total weight per period is 0.125 for $e = 0$ to 5, drops to 0.058 for $e = 6$ to 11, and turns negative from $e = 12$ on. TWFE puts the largest weights on the first six months, when the effect is smallest, and negative weights on the late periods, when the effect is largest. That is where the 0.023 comes from.

$DID_1$ through $DID_{12}$ match the CS event-study coefficients of Section 6.5 digit for digit (for example, $DID_7$, which is $e = 6$, equals 0.034), with slightly different standard errors. The six placebos correspond to $e = -2$ through $-7$; their estimates lie between 0.005 and 0.012, and their joint test gives p = 0.725, while the joint test of the effects gives p < 0.0001. The placebos are level differences relative to $e = -1$, equal to the pre-treatment CS coefficients under a universal base period (Section 6.9), not the short period-by-period differences of the default used in Section 6.5.

The software also reports an average total effect per unit of treatment, 0.052 (SE 0.007). In a binary, absorbing design, this is the average of the effects of all treated cells in the first 12 months after treatment, where each cell's effect is the cumulative result of having been treated for $\ell$ periods (de Chaisemartin and D'Haultfœuille 2026 §3.3). In Kestrel every $\ell$ has 52 switchers, so it is exactly the simple average of $DID_1$ through $DID_{12}$. It is neither the instantaneous effect of some month nor the effect after 6.5 months of exposure. The 6.5 is the average number of periods over which the effect of one treatment has accumulated, and dCDH note that only the total effect divided by it is the average effect of one period of exposure.

### 6.8 The Grand Reconciliation

| Estimator | Estimate | SE | Share of truth |
|---|---|---|---|
| True simple ATT (DGP) | 0.088 | | 100% |
| Static TWFE (feols) | 0.0227 | 0.0087 | 25.8% |
| CS simple ATT (att_gt, not-yet-treated) | 0.082 | 0.007 | 93.1% |
| SA average ATT (sunab) | 0.081 | 0.008 | 91.7% |
| BJS imputation (did_imputation) | 0.078 | 0.005 | 88.5% |
| Gardner two-stage (did2s) | 0.078 | 0.005 | 88.5% |
| ETWFE, not-yet-treated (etwfe) | 0.078 | 0.005 | 88.5% |
| ETWFE, never-treated (etwfe) | 0.081 | 0.008 | 91.7% |

Reconciliation of the dynamic coefficients (the truth is 0 at every $e < 0$):

| $e$ | True | naive TWFE (SE) | sunab (SE) | CS dynamic (SE) | imputation (SE) | $DID_\ell$ (SE) |
|---|---|---|---|---|---|---|
| -6 | 0.000 | 0.001 (0.009) | 0.009 (0.008) | -0.002 (0.006) | 0.006 (0.007) | 0.005 (0.008) |
| -3 | 0.000 | 0.009 (0.008) | 0.009 (0.009) | 0.001 (0.007) | 0.011 (0.008) | 0.012 (0.008) |
| 0 | 0.019 | 0.015 (0.007) | 0.014 (0.007) | 0.012 (0.006) | 0.005 (0.006) | 0.012 (0.006) |
| 3 | 0.022 | 0.020 (0.008) | 0.011 (0.010) | 0.015 (0.008) | 0.008 (0.007) | 0.015 (0.009) |
| 6 | 0.041 | 0.036 (0.009) | 0.032 (0.012) | 0.034 (0.011) | 0.027 (0.008) | 0.034 (0.011) |
| 9 | 0.093 | 0.089 (0.009) | 0.087 (0.011) | 0.089 (0.011) | 0.082 (0.008) | 0.089 (0.011) |

The $DID_\ell$ column reports, at $e \geq 0$, the effect for $\ell = e + 1$, and at $e < 0$, the placebos. The pre-treatment columns are defined differently: naive, sunab, and $DID_\ell$ give level differences relative to $e = -1$; CS gives short period-by-period differences; and imputation gives the coefficients from re-estimating on untreated observations with leads added. At longer horizons, the truth is 0.141 at $e = 12$, 0.164 at $e = 18$, and 0.184 at $e = 23$. Naive estimates 0.130 (0.010) at $e = 12$, and 0.162 (0.011) for the endpoint bin $e \geq 18$, where the true values within the bin range from 0.164 to 0.184. Both sunab and CS estimate 0.118 at $e = 12$ (SE 0.014 and 0.015 respectively). At the long end only one or two cohorts remain, so the noise grows and the intervals widen, which is exactly the endpoint problem of Section 4.2.

### 6.9 Rambachan-Roth: How Wrong Must Parallel Trends Be to Flip the Result?

Parallel trends in Kestrel holds by construction, but we can pretend not to know that and ask a question that only real research would ask: if the post-treatment deviation from parallel trends is not zero, how large must it be before the effects at $e = 6$ and $e = 12$ stop being significant? HonestDiD needs a set of event-study coefficients relative to a single base period, together with their covariance matrix, so the CS event study here uses `base_period = "universal"`. The pre-treatment coefficients then become level differences relative to $e = -1$; for $e = -6$ through $-2$ they are 0.005, 0.005, 0.011, 0.012, and 0.010, exactly the placebos of Section 6.7. The post-treatment coefficients are identical to those of Section 6.5.

Let $b_e$ be the counterfactual trend difference at event time $e$, that is, the difference between the paths of $Y(0)$ in the treated and control groups (Rambachan and Roth write it as $\delta$, but in this chapter $\delta$ is reserved for Sun-Abraham), normalized so that $b_{-1} = 0$. Under no anticipation, the pre-treatment $b_e$ are the population values of the pre-treatment coefficients. Smoothness (SD), for all $e$, and relative magnitudes (RM), for $e \geq -1$, require respectively

$$\big|(b_{e+1} - b_e) - (b_e - b_{e-1})\big| \leq M, \qquad \big|b_{e+1} - b_e\big| \leq \bar M \max_{s \leq -2} \big|b_{s+1} - b_s\big|$$

The value of $M$ or $\bar M$ at which the robust interval first includes zero is called the breakdown value.

```r
cs_u <- att_gt(yname = "lgmv", tname = "month", idname = "city_id", gname = "g",
               data = dt, control_group = "notyettreated", base_period = "universal",
               clustervars = "city_id", bstrap = FALSE, cband = FALSE)
es <- aggte(cs_u, type = "dynamic", min_e = -6, max_e = 12, bstrap = FALSE, cband = FALSE)

keep <- es$egt != -1                                  # drop the reference period
IF   <- es$inf.function$dynamic.inf.func.e[, keep]
V    <- crossprod(IF) / nrow(IF)^2                    # 5 pre + 13 post coefficients
beta <- es$att.egt[keep]
l_e6 <- matrix(as.numeric(0:12 == 6), ncol = 1)       # target: theta(6)

sd_e6 <- createSensitivityResults(
  betahat = beta, sigma = V, numPrePeriods = 5, numPostPeriods = 13, l_vec = l_e6,
  Mvec = c(0, 0.0001, 0.0002, 0.0005, 0.001, 0.0015, 0.002))
rm_e6 <- createSensitivityResults_relativeMagnitudes(
  betahat = beta, sigma = V, numPrePeriods = 5, numPostPeriods = 13, l_vec = l_e6,
  Mbarvec = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1),
  bound = "deviation from parallel trends", method = "C-LF",
  grid.lb = -0.6, grid.ub = 0.6, gridPoints = 6001)
```

For $e = 12$, simply replace the 6 in l_e6 with 12. The breakdown values come from a point-by-point scan on a finer grid:

| Target | Point estimate | Conventional 95% CI (analytic SE) | SD interval at $M = 0$ | SD breakdown $M$ | RM breakdown $\bar M$ |
|---|---|---|---|---|---|
| $e = 6$ | 0.034 | [0.013, 0.054] | [0.00, 0.07] | about 0.00016 | about 0.2 |
| $e = 12$ | 0.118 | [0.090, 0.146] | [0.07, 0.19] | about 0.00086 | 0.40 |

The scan step is 0.00001 for SD and 0.005 for RM. The two restrictions should be read together. SD constrains bends. At $M = 0$ it requires linear extrapolation of the pre-treatment slope, and the uncertainty of that extrapolation grows with the horizon. So the interval at $e = 6$ is already much wider than the conventional one, its lower bound nearly touches zero, and an $M$ of only about 0.00016 is enough to bring zero inside. This number looks alarmingly small, but it should be compared with the bends observed before treatment: across the five pre-treatment coefficients plus the zero at $e = -1$, the largest absolute second difference is about 0.008, fifty times as large. In other words, once you concede that post-treatment bends can be as large as pre-treatment ones, neither conclusion survives under SD.

RM uses a different ruler, with the largest pre-treatment period-by-period change (0.010) as its unit. The effect at $e = 6$ flips when $\bar M$ is about 0.2, and the effect at $e = 12$ only at 0.40; that is, the effect at $e = 12$ can tolerate post-treatment period-by-period changes up to forty percent of the largest pre-treatment change. Both are far below the natural benchmark $\bar M = 1$, "post-treatment is no worse than pre-treatment."

These benchmarks, however, are themselves noise here. Parallel trends holds exactly in the DGP, the standard errors of the pre-treatment coefficients are between 0.007 and 0.008, and both the 0.008 second difference and the 0.010 period-by-period change come from sampling error. Calibrating $M$ or $\bar M$ with noisy pre-treatment estimates overstates the plausible bends and deviations, so in real research, too, you should first ask how much of the pre-treatment variation is noise. One more qualification: $e = 12$ is identified only by cohorts 13 and 19 (cohort 25 is observed only up to $e = 11$), and its true value of 0.141 is the average over the 34 cells of these two batches of cities.

To sum up the reconciliation: on the same data, with every identification assumption holding exactly, static TWFE delivers a quarter of the truth and the modern estimators deliver around ninety percent; the CS and Sun-Abraham intervals cover the truth, while the BJS family comes in about two standard errors low. The differences come entirely from the aggregation weights, and the Goodman-Bacon decomposition and the dCDH weights explain, from the comparison side and the cell side, where the missing three quarters went. $DID_\ell$ and CS, did2s and BJS, and ETWFE and BJS or Sun-Abraham coincide numerically, which shows that they are different ways of writing the same recipe. HonestDiD adds a reminder: even on these data, where parallel trends holds exactly, how robust the conclusions are depends on how strongly you restrict post-treatment deviations from parallel trends.

## 7 Failure Modes of Identification and Robustness

In the simulation the identification assumptions hold by construction; in real research they can fail at any time. This section goes through the most common failure modes and the practical responses.

The problem with pre-trends tests runs one layer deeper than "it is not proof." Roth (2022) gives two quantitative results. In common event-study designs, the pre-trends test has low power against trend differences large enough to flip the conclusion, so a failure to reject often means only that the data are insufficient, not that the trends are parallel. More subtle is the problem of conditional inference. Proceeding with inference only after the test passes is itself a screen, and the screen distorts the distribution of the subsequent estimates. When there is a small real violation, the draws that pass the screen are exactly those in which pre-period noise offset the violation; and because the pre and post coefficients share a base period, their noise is positively correlated, so the post-period bias is amplified rather than reduced. The practical response is to replace a binary test with a quantitative statement: report the test's power against a violation of a given size, plot confidence intervals rather than just reporting stars, and do not drop parallel trends from the discussion of limitations just because the pre test passed.

Rambachan and Roth (2023) push this idea to its logical end: rather than test "whether the violation is zero," ask "how large the violation must be to flip the result." Their framework no longer assumes that parallel trends holds exactly. Instead it confines the possible post-treatment violations to a set and delivers robust confidence sets that are uniformly valid over that set. The two common sets correspond to two worries: SD worries about smoothly evolving long-run trends, RM about shocks of the same magnitude as those before treatment. The actual run on Kestrel in Section 6.9 shows how to read the two numbers. The effect at $e = 12$ flips under RM only at $\bar M = 0.40$, but under SD it flips at an $M$ of about 0.00086, while the bends observed before treatment are as large as about 0.008. Reporting only SD makes the conclusion look hopelessly fragile; reporting only RM hides the fragility of extrapolation. The honest approach is to report both and to say what scale $M$ and $\bar M$ are each measured against. When writing a paper, rather than saying "pre-trends are parallel," argue how implausible a trend deviation it would take to overturn the conclusion.

Adjacent to anticipation but different in mechanism is the Ashenfelter dip, which is at heart a failure of parallel trends rather than of no anticipation. In a training-program setting, people enter the program because their income has temporarily dropped ($Y_{i1}(0)$ falls below some threshold). This amounts to truncated sampling on the pre-period shock, so the treated group's untreated path contains a stretch of pure mean reversion that differs systematically from the control group's $Y(0)$ trend, and DiD counts the rebound as a training effect and overstates it. The line between the two is this. Anticipation means that units foresee future treatment and adjust their behavior in advance, a failure of no anticipation (A2). The dip means that the selection rule picks units that happen to be in a temporary trough, a failure of parallel trends (A1). Both contaminate the pre-period, and the responses are similar.

The platform version is not far-fetched at all: if Kestrel's operations team deliberately cuts commissions in cities that just had a bad quarter of GMV, hoping for a boost, the rebound will be credited to the commission cut. Diagnosis and response: plot the event study and check whether the last few periods with negative $e$ dip abnormally; move the base period from $e = -1$ to an earlier period (exactly what the CS anticipation argument is for); re-estimate excluding the dip periods; and use a longer pre window to see whether the trend returns to its course.

Endogenous treatment timing is the norm, not the exception, in platform research. A platform that cuts commissions first in its fastest-growing cities threatens both parallel trends and the interpretability of the weights. There are three lines of defense. First, use the not-yet-treated rather than the never-treated as controls: late-treated cities are also cities the platform will select, just scheduled later, and on the selection dimension they are usually far more comparable to the early-treated cities than cities that are never selected. Second, find quasi-exogenous variation in treatment timing to use as an IV or as a design argument: if, for example, the schedule of commission cuts is set by bottlenecks unrelated to demand growth, such as engineering resources or compliance approval, write that institutional detail out in full. Third, combine DiD with synthetic control or synthetic DiD (Chapter 22): build a weighted counterfactual for each treated city, relaxing the implicit constraint that "all control cities get equal weight."

Diagnosing spillover must be specific to the mechanism. Sellers re-registering across cities, buyers ordering across cities, and platform-wide reallocation of resources each leave their own observable footprint. A practical check is to stratify the control cities by geographic or economic proximity to the treated cities. If the estimated effect using "near controls" differs systematically from the one using "far controls," spillover is essentially confirmed. You can then use a donut design (dropping control cities in the border band), redefine market boundaries, or admit that SUTVA fails and make the spillover itself an estimand. If the suspicion is a platform-wide equilibrium effect (the commission structure changed the seller pricing equilibrium of the whole platform), the DiD framework does not apply at all, and a structural model or a platform-level experiment is needed.

Section 2 laid the groundwork for functional-form sensitivity: parallel trends in logs and parallel trends in levels are two different assumptions. The practical rule is to run the main result in both forms and put both in the appendix. If the conclusions conflict, do not pick the nicer one. Admit in the main text that identification is fragile to functional form, and explain which form your mechanism supports (a change in the commission rate acts proportionally on GMV, so log has a mechanistic basis).

Few clusters are dealt with in Section 4.7; here we stress only the point most easily overlooked. The most serious problem is not a small total number of clusters but a small number of treated clusters. If only 3 of 50 states passed some law, the effective sample size is close to 3, not 50, and CR1 will hand you a pretty but false star.

The rule on bad controls: do not control for any variable that responds to treatment. Controlling for "the platform's marketing spend in the city" seems to make the regression cleaner, but marketing spend is itself a downstream result of the commission policy. Controlling for it shuts off one channel through which the effect operates, and may also introduce collider bias. The same kind of error is stuffing a lagged dependent variable into a DiD regression to "absorb the trend." A lag together with fixed effects brings Nickell bias, the estimand changes, and it is no diagnostic for parallel trends. To deal with trend differences, the right approach is the conditioning of Section 4.6 or the Rambachan-Roth sensitivity analysis. If parallel trends holds only after controlling for a covariate that is itself affected by treatment, Caetano, Callaway, Payne and Sant'Anna (2026) point out that simply dropping it is often wrong. You can argue that its pre-treatment value alone suffices (which is the CS estimator with covariates), or, under an unconfoundedness condition on the covariate, switch to imputation and DML estimators (the R package badcontrols).

Taken together, these failure modes show that the credibility of DiD lies not in the ingenuity of the regression specification but in the design itself: why treatment happens to these units at these times, and whether that mechanism is unrelated to the counterfactual trend of the outcome. The pre-trends plot, the Rambachan-Roth sensitivity analysis, not-yet-treated controls, and functional-form robustness checks are all evidence on this one question, and none of them can replace a substantive judgment about the treatment assignment mechanism.

## 8 Further Reading

::: {.readings}
Required reading, in suggested order:

- Baker, Callaway, Cunningham, Goodman-Bacon and Sant'Anna (2026, Journal of Economic Literature). The practitioner's guide to read first: it organizes covariates, weighting, multiple periods, and staggered treatment in one unified framework, and the original papers come much more easily after it.
- Roth, Sant'Anna, Bilinski and Poe (2023, Journal of Econometrics). Organizes the econometric DiD literature of 2018 to 2023 into a single map; it complements the previous entry.
- Goodman-Bacon (2021, Journal of Econometrics). The decomposition theorem for static TWFE, the full version of Section 4.3 of this chapter; focus on the comparative statics of the weights.
- Callaway and Sant'Anna (2021, Journal of Econometrics). The original paper on the group-time ATT framework; focus on Theorem 1, the difference between the two control-group assumptions, and the aggregation section.
- Bertrand, Duflo and Mullainathan (2004, QJE). The origin of the inference problem; the design of its placebo experiment is itself worth borrowing.
- Baker, Larcker and Wang (2022, Journal of Financial Economics). Replays the failure of staggered TWFE with finance data; useful as a self-check before writing a paper.

Further reading:

- Wooldridge (2010) §6.5, §10.5.6, and §10.7.3. Textbook treatments of the regression forms of DiD and DDD, the conditions for using fixed effects in policy analysis, and the Mundlak regression.
- Sun and Abraham (2021, Journal of Econometrics). The rigorous statement of the contamination theorem and the IW estimator; read Propositions 1 to 3.
- Borusyak, Jaravel and Spiess (2024, Review of Economic Studies). The imputation framework and the efficiency result; note its stronger requirement on parallel trends.
- de Chaisemartin and D'Haultfœuille (2020, American Economic Review). The theorem and diagnostic for negative TWFE weights.
- de Chaisemartin and D'Haultfœuille (2023c, Econometrics Journal). The authors' own survey, comparing estimators that are robust to heterogeneous effects in a single notation.
- de Chaisemartin and D'Haultfœuille (2026, Review of Economics and Statistics). The original paper on $DID_\ell$: event studies with non-binary and switching treatments.
- de Chaisemartin and D'Haultfœuille (2023a), the textbook draft. Princeton University Press is scheduled to publish the final version in December 2026; until it appears, read this SSRN draft.
- Wooldridge (2025, Empirical Economics). The two-way Mundlak regression and ETWFE, showing that TWFE with a full set of interactions is imputation.
- Wing, Freedman and Hollingsworth (2024, NBER Working Paper) and Dube, Girardi, Jordà and Taylor (2025, Journal of Applied Econometrics). The original papers on stacked DiD and LP-DiD.
- Caetano and Callaway (2026, working paper). The hidden linearity bias of TWFE under conditional parallel trends.
- Caetano, Callaway, Payne and Sant'Anna (2026, working paper). How bad controls affected by treatment can still serve as controls in DiD; the source of the bad-controls paragraph in Section 7.
- Ortiz-Villavicencio and Sant'Anna (2025, working paper). When DDD cannot be written as the difference of two DiDs.
- Roth (2022, AER: Insights). The power of pre-trends tests and pre-test distortions.
- Rambachan and Roth (2023, Review of Economic Studies). The honest inference framework; the original source of the SD and RM sensitivity analyses.
- Roth and Sant'Anna (2023, Econometrica). When parallel trends is invariant to functional form.
- Sant'Anna and Zhao (2020, Journal of Econometrics). The doubly robust DiD estimator, the theoretical basis of the did package's default method.
- Chiu, Lan, Liu and Xu (2026, American Political Science Review). A reanalysis of 49 political science studies that use TWFE: robust estimators broadly agree with the original conclusions but vary a great deal, many studies lack evidence supporting parallel trends, and power is insufficient once heterogeneous effects are allowed for.
:::

::: {.apa-refs}
- Baker, A., Callaway, B., Cunningham, S., Goodman-Bacon, A., & Sant'Anna, P. H. C. (2026). Difference-in-differences designs: A practitioner's guide. *Journal of Economic Literature, 64*(2), 498-557. https://doi.org/10.1257/jel.20251650
- Baker, A. C., Larcker, D. F., & Wang, C. C. Y. (2022). How much should we trust staggered difference-in-differences estimates? *Journal of Financial Economics, 144*(2), 370-395. https://doi.org/10.1016/j.jfineco.2022.01.004
- Bell, R. M., & McCaffrey, D. F. (2002). Bias reduction in standard errors for linear regression with multi-stage samples. *Survey Methodology, 28*(2), 169-181.
- Bertrand, M., Duflo, E., & Mullainathan, S. (2004). How much should we trust differences-in-differences estimates? *The Quarterly Journal of Economics, 119*(1), 249-275. https://doi.org/10.1162/003355304772839588
- Borusyak, K., Jaravel, X., & Spiess, J. (2024). Revisiting event-study designs: Robust and efficient estimation. *The Review of Economic Studies, 91*(6), 3253-3285. https://doi.org/10.1093/restud/rdae007
- Caetano, C., & Callaway, B. (2026). *Difference-in-differences when parallel trends holds conditional on covariates* (arXiv:2406.15288v3). arXiv. https://doi.org/10.48550/arXiv.2406.15288
- Caetano, C., Callaway, B., Payne, S., & Sant'Anna, H. (2026). *Difference-in-differences with "bad controls"* (arXiv:2608.03881v2). arXiv. https://doi.org/10.48550/arXiv.2608.03881
- Callaway, B., & Sant'Anna, P. H. C. (2021). Difference-in-differences with multiple time periods. *Journal of Econometrics, 225*(2), 200-230. https://doi.org/10.1016/j.jeconom.2020.12.001
- Cameron, A. C., Gelbach, J. B., & Miller, D. L. (2008). Bootstrap-based improvements for inference with clustered errors. *The Review of Economics and Statistics, 90*(3), 414-427. https://doi.org/10.1162/rest.90.3.414
- Card, D., & Krueger, A. B. (1994). Minimum wages and employment: A case study of the fast-food industry in New Jersey and Pennsylvania. *American Economic Review, 84*(4), 772-793.
- Card, D., & Krueger, A. B. (1998). *A reanalysis of the effect of the New Jersey minimum wage increase on the fast-food industry with representative payroll data* (NBER Working Paper No. 6386). National Bureau of Economic Research. https://doi.org/10.3386/w6386
- Card, D., & Krueger, A. B. (2000). Minimum wages and employment: A case study of the fast-food industry in New Jersey and Pennsylvania: Reply. *American Economic Review, 90*(5), 1397-1420. https://doi.org/10.1257/aer.90.5.1397
- Chiu, A., Lan, X., Liu, Z., & Xu, Y. (2026). Causal panel analysis under parallel trends: Lessons from a large reanalysis study. *American Political Science Review, 120*(1), 245-266. https://doi.org/10.1017/S0003055425000243
- Conley, T. G., & Taber, C. R. (2011). Inference with "difference in differences" with a small number of policy changes. *The Review of Economics and Statistics, 93*(1), 113-125. https://doi.org/10.1162/rest_a_00049
- de Chaisemartin, C., & D'Haultfœuille, X. (2020). Two-way fixed effects estimators with heterogeneous treatment effects. *American Economic Review, 110*(9), 2964-2996. https://doi.org/10.1257/aer.20181169
- de Chaisemartin, C., & D'Haultfœuille, X. (2023a). *Causal inference with differences-in-differences: Credible answers to hard questions* [Unpublished manuscript]. SSRN. https://doi.org/10.2139/ssrn.4487202
- de Chaisemartin, C., & D'Haultfœuille, X. (2023c). Two-way fixed effects and differences-in-differences with heterogeneous treatment effects: A survey. *The Econometrics Journal, 26*(3), C1-C30. https://doi.org/10.1093/ectj/utac017
- de Chaisemartin, C., & D'Haultfœuille, X. (2026). Difference-in-differences estimators of intertemporal treatment effects. *The Review of Economics and Statistics, 108*(4), 863-880. https://doi.org/10.1162/rest_a_01414
- Dube, A., Girardi, D., Jordà, Ò., & Taylor, A. M. (2025). A local projections approach to difference-in-differences. *Journal of Applied Econometrics, 40*(7), 741-758. https://doi.org/10.1002/jae.70000
- Gardner, J. (2022). *Two-stage differences in differences* (arXiv:2207.05943v1). arXiv. https://doi.org/10.48550/arXiv.2207.05943
- Goodman-Bacon, A. (2021). Difference-in-differences with variation in treatment timing. *Journal of Econometrics, 225*(2), 254-277. https://doi.org/10.1016/j.jeconom.2021.03.014
- Greenwood, B. N., & Wattal, S. (2015). *Show me the way to go home: An empirical investigation of ride sharing and alcohol related motor vehicle homicide* [Working paper]. SSRN. https://doi.org/10.2139/ssrn.2557612
- Greenwood, B. N., & Wattal, S. (2017). Show me the way to go home: An empirical investigation of ride-sharing and alcohol related motor vehicle fatalities. *MIS Quarterly, 41*(1), 163-187. https://doi.org/10.25300/MISQ/2017/41.1.08
- Liu, L., Wang, Y., & Xu, Y. (2024). A practical guide to counterfactual estimators for causal inference with time-series cross-sectional data. *American Journal of Political Science, 68*(1), 160-176. https://doi.org/10.1111/ajps.12723
- MacKinnon, J. G., & Webb, M. D. (2020). Randomization inference for difference-in-differences with few treated clusters. *Journal of Econometrics, 218*(2), 435-450. https://doi.org/10.1016/j.jeconom.2020.04.024
- Neumark, D., & Wascher, W. (2000). Minimum wages and employment: A case study of the fast-food industry in New Jersey and Pennsylvania: Comment. *American Economic Review, 90*(5), 1362-1396. https://doi.org/10.1257/aer.90.5.1362
- Ortiz-Villavicencio, M., & Sant'Anna, P. H. C. (2025). *Better understanding triple differences estimators* (arXiv:2505.09942v3). arXiv. https://doi.org/10.48550/arXiv.2505.09942
- Rambachan, A., & Roth, J. (2023). A more credible approach to parallel trends. *The Review of Economic Studies, 90*(5), 2555-2591. https://doi.org/10.1093/restud/rdad018
- Roth, J. (2022). Pretest with caution: Event-study estimates after testing for parallel trends. *American Economic Review: Insights, 4*(3), 305-322. https://doi.org/10.1257/aeri.20210236
- Roth, J., & Sant'Anna, P. H. C. (2023). When is parallel trends sensitive to functional form? *Econometrica, 91*(2), 737-747. https://doi.org/10.3982/ECTA19402
- Roth, J., Sant'Anna, P. H. C., Bilinski, A., & Poe, J. (2023). What's trending in difference-in-differences? A synthesis of the recent econometrics literature. *Journal of Econometrics, 235*(2), 2218-2244. https://doi.org/10.1016/j.jeconom.2023.03.008
- Sant'Anna, P. H. C., & Zhao, J. (2020). Doubly robust difference-in-differences estimators. *Journal of Econometrics, 219*(1), 101-122. https://doi.org/10.1016/j.jeconom.2020.06.003
- Sun, L., & Abraham, S. (2021). Estimating dynamic treatment effects in event studies with heterogeneous treatment effects. *Journal of Econometrics, 225*(2), 175-199. https://doi.org/10.1016/j.jeconom.2020.09.006
- Wing, C., Freedman, S., & Hollingsworth, A. (2024). *Stacked difference-in-differences* (NBER Working Paper No. 32054). National Bureau of Economic Research. https://doi.org/10.3386/w32054
- Wooldridge, J. M. (2010). *Econometric analysis of cross section and panel data* (2nd ed.). MIT Press.
- Wooldridge, J. M. (2025). Two-way fixed effects, the two-way Mundlak regression, and difference-in-differences estimators. *Empirical Economics, 69*(5), 2545-2587. https://doi.org/10.1007/s00181-025-02807-z
:::
