---
title: "Difference-in-Differences & Event Studies"
subtitle: "From the 2×2 Table to Staggered Adoption"
seriesline: "Foundations of Information Systems Economics · Chapter 13"
series: "Foundations of Information Systems Economics: Empirical Methods for Information Systems Economics"
sidebartitle: "Foundations of Information Systems Economics"
sidebarsub: "Theory · Methods · Applications"
chapterlabel: "Chapter 13 · Difference-in-Differences & Event Studies"
---

## Introduction

Kestrel Marketplace rolls out a new policy city by city, in batches. Cities grow faster after they go live, and the event-study plot draws the familiar post-treatment upswing. It looks as though the policy worked. But the cities did not all go live on the same day: the early-treated cities end up serving as a control group for the late-treated ones, while the policy effect itself accumulates over time since going live. In this situation a single standard two-way fixed effects regression can blend several individually reasonable local comparisons into one number whose sign is not even reliable.

The original appeal of DiD comes from a two-row, two-column table: the before-after change for the treated group, minus the contemporaneous change for the control group. This difference of differences wipes out fixed between-group differences and common time shocks, but it does not conjure a counterfactual out of thin air. It still requires an answer to one question: had the policy not happened, how would the two groups have moved relative to each other? "The pre-treatment trends look about the same" is supporting evidence, not the parallel trends assumption itself.

This chapter first makes this counterfactual explicit in the most transparent 2×2 case, then moves to the staggered adoption that is far more common in platform research. The Goodman-Bacon decomposition will reveal exactly which comparisons get mixed into a traditional TWFE coefficient; the Callaway-Sant'Anna, Sun-Abraham, and Borusyak-Jaravel-Spiess methods reorganize the control groups and the event time so that every weight carries a clear meaning. The question running through it all is not "which command should I run" but rather: at each point in time, who plays the never-observed counterfactual on behalf of the treated group?

## 1 A Number That Does Not Add Up

::: {.case}
Kestrel Marketplace is a fictional two-sided e-commerce platform. We tell this case with simulated data, and the payoff is that the DGP is fully known, so the behavior of every estimator can be reconciled against the truth, a teaching condition that real data can never provide. The setup: the platform operates in 60 cities over a 36-month observation window, and the outcome variable is the average log GMV of sellers in each city each month. The platform cuts the seller commission rate from 12% to 8% in three batches: 17 cities in month 13, another 17 in month 19, and a final 18 in month 25, while the remaining 8 cities keep the 12% rate throughout the window. Management's question is straightforward: how large is the effect of the commission cut on seller GMV.

Because this is a simulation, we know the answer. The true effect is positive and grows with event time: about 0.02 (log points) in the month of the cut, climbing to around 0.14 a year later, and drifting to about 0.18 at the longest horizon ($e = 23$). Averaging the realized effect over all 930 post-treatment city-month cells, the true simple ATT is 0.088.
:::

![The raw path of city-average log GMV by calendar month for each cohort, with dashed lines marking the timing of the three commission cuts (months 13, 19, and 25).](assets/fig/fig_13_raw_trends.svg)

Start with the raw data. The figure above has three features worth noting. The cities that cut commissions early sit at a higher level, which is by design, because the platform moves first in its large markets, and the level difference by itself does not threaten identification. The four groups of cities show an obvious common growth trend and seasonal fluctuation. After the commission cut, the groups gradually fan out, and the late-treated groups also start accelerating after their own treatment dates. Intuitively this looks like a policy whose effect is positive and accumulates over time.

The researcher therefore runs the textbook two-way fixed effects (TWFE) regression:

$$Y_{it} = \alpha_i + \lambda_t + \beta\, D_{it} + \varepsilon_{it}$$

where $D_{it}$ equals 1 once city $i$ has already cut its commission in month $t$, with standard errors clustered at the city level. The result is $\hat\beta_{TWFE} = 0.0227$ (SE 0.0087, $t = 2.61$, $p = 0.011$, rounded to 0.023 when we quote it later). Taken in isolation, this is a clean conclusion: statistically significant, positive, with a magnitude of 2.3%. Written into a report to management, the takeaway would be "the commission cut has an effect, but a small one, probably not worth rolling out."

The problem is that we know the truth is 0.088. Static TWFE recovers only 25.8% of the true effect, and nearly three quarters of the effect has vanished into thin air. The regression is not miswritten, the fixed effects are all there, the standard errors cluster at the right level, and parallel trends holds exactly in the data (the DGP guarantees it). So what is $\hat\beta_{TWFE}$ actually estimating? This is the main thread of the chapter. The answer arrives in Section 4: under staggered adoption (treatment in batches) combined with effects that grow with event time, static TWFE estimates a weighted average with bizarre weights, part of which sits on comparisons that use already-treated units as controls. In Section 6 we will take this 0.023 apart precisely and see where each piece of the bias comes from.

## 2 The Economic Model and the Estimand

Before we start estimating, we need to pin down two things: which causal quantity we want, and the language that writes it down precisely. The former is this section's estimand (the target quantity to estimate), the latter is the notation and model that carries it, namely potential outcomes, along with the additive structure on $Y(0)$ that appears at the end of this section.

Consider the target quantity first. The effect of the commission cut on GMV looks like one number but is really a whole family of quantities that varies across cities and over time, because different cities benefit differently and the effect also accumulates. How this family should be collapsed into a single estimable target is what the second half of this section answers formally. Now the language that describes it. The core difficulty of causality is the counterfactual: a city that has already cut its commission, had it not cut it, would have realized an outcome we can never observe. Potential outcomes are exactly the notation invented to describe such parallel-world outcomes, and this section sets it up first, since all the identification and estimation that follows is built on top of it.

Introduce the notation. Denote a city by $i$ and a month by $t$. The cohort variable $G_i$ is the period in which city $i$ is first treated, and in Kestrel $G_i \in \{13, 19, 25\}$, with never-treated cities coded $G_i = \infty$ (0 in the data). The treatment indicator is $D_{it} = \mathbf{1}\{t \geq G_i\}$. Staggered adoption means treatment is irreversible, $D_{it} \geq D_{i,t-1}$: once commission is cut it does not go back up. Event time $e = t - g$ aligns different cohorts to their own treatment dates.

Potential outcomes are indexed by the first period of treatment: $Y_{it}(g)$ is the outcome city $i$ would realize in period $t$ if it were first treated in period $g$, and $Y_{it}(0)$ is the never-treated path. The observed outcome is

$$Y_{it} = \mathbf{1}\{G_i > t\}\, Y_{it}(0) + \mathbf{1}\{G_i \leq t\}\, Y_{it}(G_i)$$

The basic estimand of this chapter is the group-time average treatment effect:

$$ATT(g,t) = \mathbb{E}\big[Y_{t}(g) - Y_{t}(0) \mid G = g\big]$$

that is, the ATT of cohort $g$ at calendar time $t$. Management's question, "how large is the effect of the commission cut on GMV," is really not one number but some average over the whole $ATT(g,t)$ matrix: you can average by cohort (which batch of cities benefits most), by event time (how quickly the effect comes up, whether it fades), or all the way down to a single grand total. First decide which average you want and with what weights, then choose the estimator, and this order cannot be reversed. The common structure of the modern estimators discussed in Section 4 is exactly "estimate $ATT(g,t)$ first, then aggregate with explicit weights," and the root disease of static TWFE is precisely that its aggregation weights are neither chosen by the researcher nor economically interpretable.

**Parallel trends** is the identification core of DiD, and it is worth understanding first through its economic content, with the formal statement deferred to Section 3. It says: had the treated group not been treated, the change in the time path of its outcome would have been the same as the control group's. What kind of outcome-generating process guarantees this? Write the untreated state as $Y_{it}(0) = \alpha_i + \lambda_t + u_{it}$, where $\alpha_i$ absorbs a city's intrinsic level, $\lambda_t$ is a common time shock, and the group means of the $u_{it}$ paths are equal, and parallel trends holds. In other words this is an additive separability assumption about the process generating $Y(0)$: unobserved factors enter additively, and the time shock acts on the treated and control groups by the same amount. This is far stronger than "I controlled for fixed effects." Fixed effects are merely an estimation device, and writing them into a regression is a purely mechanical operation; parallel trends is an economic assertion about the counterfactual world, about how the treated cities' GMV would have moved in those months where commissions were not cut. Regression syntax cannot guarantee that assertion.

Additivity also carries an often-overlooked corollary: parallel trends depends on functional form. Parallel trends that holds for log GMV generally does not hold for the GMV level, and vice versa, so in general at most one of the two holds. Roth and Sant'Anna (2023, Econometrica) prove that parallel trends holds for every monotone transformation of the outcome if and only if a very strong condition is met, roughly that treatment is approximately randomly assigned, or that the distribution of $Y(0)$ itself does not shift over time. This condition is enough to make parallel trends hold on both the log and the level scales at once; conversely, the mere fact that both log and level are parallel is not enough to imply it, because two particular transformations being parallel is far weaker than all monotone transformations being parallel. So "log or level" is not a matter of taste but a choice of which identification assumption you are willing to maintain, and of whether your estimand is a percentage effect or an absolute-quantity effect. Section 7 returns to the choice of functional form from the angle of robustness.

The reality of platform data makes parallel trends especially worth watching. A platform does not decide by lottery which cities get the commission cut first, and the first to be treated are often the fastest-growing, most competitive, or strategically most important markets, so the choice of treatment timing gets entangled with the growth trajectory. Two things need to be distinguished here: parallel trends allows the two groups to sit at different levels, and it allows the platform to deliberately move first in the cities with the largest expected effect. The only thing it forbids is that the two groups' growth paths in the untreated state differ systematically. And a platform that selects cities by growth potential may be selecting exactly the cities whose growth paths differ to begin with, which is the real threat. Section 3 dissects these three assumptions one by one.

::: {.intuition}
It becomes clearer to put DiD back on the menu of counterfactual constructions (of which matching and synthetic control are methods detailed in later modules, so here you only need to see whom each of them borrows a counterfactual from). To evaluate a treatment, you always have to borrow a counterfactual from somewhere: an event study borrows your own past as the counterfactual; matching borrows other units similar to you on observable characteristics; DiD borrows other units' change as the counterfactual for your change; synthetic control borrows the change of a weighted combination of other units. DiD and matching are a mirror image pair: matching makes almost no functional form assumption but requires selection on observables, and unobservable selection kills it; DiD makes a strong additive functional form assumption, and in exchange it tolerates arbitrary unobserved level heterogeneity. Conditional parallel trends (Section 3) then applies matching's logic to the growth of the outcome rather than its level. Choosing a design is choosing which kind of assumption you are willing to pay for a counterfactual.
:::

A summary: the implicit model of DiD is the additive separable structure of $Y(0)$, the estimand is the ATT family (generally not the ATE), and under a staggered setup the basic object is $ATT(g,t)$. Nothing in this section has yet touched "how to estimate," and this is deliberate: the credibility of the identification assumptions depends only on the design and institutional background, and has nothing to do with which estimator we use later.

## 3 Identification: Which Bias Each of the Three Assumptions Blocks

Identification logically precedes estimation. Identification asks whether, under a set of assumptions about the counterfactual world, we can write the causal quantity we want as a function of the observable data distribution; the answer to this question depends only on the design and the institutional background, and has nothing to do with whether we later use OLS or GMM, large samples or small. This section first proves DiD identification line by line in the simplest 2x2 setup, then dissects the three identification assumptions one by one. The reason it is worth spending so much space on the 2x2 case is that all the troubles of staggered adoption are permutations and combinations of this skeleton, and any problem that cannot be untangled in 2x2 only gets thornier once there are many periods and many groups.

### 3.1 The 2x2 Setup and the Identification Theorem

Start from the most naive comparison any platform operator would make. To learn whether the commission cut had an effect, the most direct approach is to look at the GMV change in the treated cities before and after the cut, which is called a single difference. Its flaw is obvious: even without the cut, GMV would have risen anyway with the platform's overall trajectory, the seasons, and the macro environment, and subtracting before from after folds all of these policy-irrelevant time changes into the effect. Whether it ends up overstating or understating depends on the direction of the overall trajectory: a single difference is too high when the market is trending up, and too low when the market was in fact trending down. To strip out the time change, you need a batch of control cities that did not cut commissions over the same period, and their before-after change serves as a ruler measuring how much the market itself rose over this stretch. Take the treated cities' before-after change, subtract the control cities' contemporaneous before-after change, and what remains is the net effect of the policy. This one step, a difference of differences, is the whole intuition of difference-in-differences.

A minimal numerical example makes this concrete. Suppose the treated cities' monthly GMV rises from 120 before the cut to 150 after, a gain of 30; over the same period a batch of control cities rises from 100 to 110, a gain of 10. The two groups start at different levels, and the treated cities are the larger market, but this does not matter, because DiD cares about how much each rose, not whose starting point is higher. The control cities did not cut commissions yet still rose by 10, which tells us the market contributed roughly 10 of growth over this stretch. So subtract that 10 from the treated cities' 30 to get 20, which is the part attributable to the commission cut. Looking only at the treated cities would overstate the effect as 30, precisely because the market was not netted out. DiD does exactly this subtraction, and all of its credibility rests on one question: is the 10 the control cities rose the same 10 the treated cities would have risen had they not cut commissions? The rest of this section is about making this question precise.

That 20 in the example above, in the potential-outcomes notation of Section 2, is the average treatment effect on the treated; what makes it hard is entirely the counterfactual $Y_{i2}(0)$, the outcome the treated cities would have had absent the cut. Below we write this 2x2 case out formally.

Two periods $t \in \{1, 2\}$, the treated group $D_i = 1$ is treated in $t = 2$, and the control group $D_i = 0$ is never treated. The target is the average treatment effect on the treated in $t = 2$, $ATT = \mathbb{E}[Y_{i2}(1) - Y_{i2}(0) \mid D_i = 1]$, and the difficulty lies in $Y_{i2}(0)$, the outcome the treated group would realize had it not been treated, which is never observed. The DiD idea is to substitute the control group's before-after change for this missing counterfactual change of the treated group, and whether the substitution is valid depends on the three assumptions below.

::: {.assumption}
**A1 (parallel trends)** The treated and control groups have the same expected growth in their untreated-state outcomes:

$$\mathbb{E}\big[Y_{i2}(0) - Y_{i1}(0) \mid D_i = 1\big] = \mathbb{E}\big[Y_{i2}(0) - Y_{i1}(0) \mid D_i = 0\big]$$
:::

::: {.assumption}
**A2 (no anticipation)** The pre-treatment outcome is unaffected by future treatment: for the treated group, $Y_{i1} = Y_{i1}(0)$.
:::

::: {.assumption}
**A3 (SUTVA / no spillovers)** Each unit's outcome depends only on its own treatment status, and treatment does not spill over; in particular, the control group's observed outcome is its $Y(0)$.
:::

::: {.theorem}
Under A1 through A3, DiD identifies the ATT:

$$\Delta_{DD} \equiv \big(\mathbb{E}[Y_{i2} \mid D_i = 1] - \mathbb{E}[Y_{i1} \mid D_i = 1]\big) - \big(\mathbb{E}[Y_{i2} \mid D_i = 0] - \mathbb{E}[Y_{i1} \mid D_i = 0]\big) = ATT$$

where $ATT = \mathbb{E}[Y_{i2}(1) - Y_{i2}(0) \mid D_i = 1]$.
:::

The derivation is only three steps, each using up exactly one assumption. Write $\Delta Y_i = Y_{i2} - Y_{i1}$. For the treated group, A2 guarantees $Y_{i1} = Y_{i1}(0)$, and splitting $Y_{i2} = Y_{i2}(1)$ into effect plus counterfactual:

$$\mathbb{E}[\Delta Y_i \mid D_i = 1] = \underbrace{\mathbb{E}\big[Y_{i2}(1) - Y_{i2}(0) \mid D_i = 1\big]}_{ATT} + \underbrace{\mathbb{E}\big[Y_{i2}(0) - Y_{i1}(0) \mid D_i = 1\big]}_{\gamma(1)}$$

For the control group, A3 guarantees its observed outcome is the $Y(0)$ path throughout:

$$\mathbb{E}[\Delta Y_i \mid D_i = 0] = \underbrace{\mathbb{E}\big[Y_{i2}(0) - Y_{i1}(0) \mid D_i = 0\big]}_{\gamma(0)}$$

Subtracting gives

$$\Delta_{DD} = ATT + \big(\gamma(1) - \gamma(0)\big)$$

A1 says $\gamma(1) = \gamma(0)$, the bias term vanishes, and $\Delta_{DD} = ATT$. QED. This decomposition is the skeleton of the whole chapter: $\gamma(1) - \gamma(0)$ is the difference between the treated and control groups' counterfactual trends, and DiD identifies the ATT if and only if this difference is zero; every failure and repair that follows is about when this term is nonzero and how to push it back to zero.

### 3.2 Anatomy of Parallel Trends

The name parallel trends is best understood as a picture. Draw the average GMV of the treated and control groups each as a line over time. What DiD requires is never that the two lines sit at the same height, since the treated cities are large markets and naturally sit higher, and this does not matter; it requires that in the hypothetical world where no one is treated, the two lines move in parallel, that is, they rise and fall by the same amount. The trouble is that in real data we only see the control group's line, untreated throughout, and a short stretch of the treated group's line before the cut; after treatment, how the treated cities would have moved absent the cut is a dashed line that cannot be drawn. The parallel trends assumption is precisely that this invisible dashed line stays parallel to the control group's line. The entire success or failure of DiD hangs on whether this dashed line is parallel.

![Geometric intuition for parallel trends: the control group (navy) and the treated group's pre segment (red) move in parallel; after treatment, the gap between the treated group's actual path (solid red) and its counterfactual path absent the cut (dashed red, parallel to the control group) is the effect ATT.](assets/fig/fig_13_parallel_intuition.svg)

Below we translate this geometric intuition into math that can tell right from wrong.

A1 is the strongest of the three assumptions, and the one most worth turning over repeatedly. First look at what it actually constrains. Write the untreated state in additive form

$$Y_{it}(0) = \alpha_i + \lambda_t + u_{it}$$

where $\alpha_i$ is a city's intrinsic level, $\lambda_t$ is a time shock that treats all units alike, and $u_{it}$ is the remaining idiosyncratic shock. The two lines moving in parallel from the previous paragraph, put into algebra, is that the two groups' expected period-by-period change is equal; and since $\alpha_i$ and $\lambda_t$ both cancel in the differencing, this requirement finally falls only on the idiosyncratic term $u$. So under this formulation A1 is equivalent to a condition on $u$, $\mathbb{E}[u_{i2} - u_{i1} \mid D_i = 1] = \mathbb{E}[u_{i2} - u_{i1} \mid D_i = 0]$. This statement draws the line between what parallel trends allows and forbids. It allows $\alpha_i$ to differ systematically between the two groups, so treated cities being intrinsically large markets ($\mathbb{E}[\alpha_i \mid D=1] \neq \mathbb{E}[\alpha_i \mid D=0]$) is no obstacle at all, because the fixed level difference is wiped out by the double differencing. It allows $\lambda_t$ to be arbitrarily large and arbitrarily nonlinear, as long as it acts on the two groups by the same amount. It even allows $u_{it}$ to have arbitrary individual heterogeneity and time fluctuation, as long as the average changes of these fluctuations are equal across groups. The one thing it forbids is an inconsistent drift in the group means of $u$, that is, the treated group's average growth path in the untreated state deviating systematically from the control group's.

This distinguishes parallel trends from another assumption that sounds similar but is different in nature and often confused with it. The selection-on-observables class of assumptions (unconfoundedness) requires that the levels of the treated and untreated states be comparable after conditioning on covariates, and it forbids $\alpha_i$ from correlating with treatment; parallel trends is far more permissive on the level dimension, allowing arbitrary, even highly treatment-correlated, fixed unobserved heterogeneity, at the price of being a knife-edge assumption on the trend dimension, where a hair of differential drift breaks it. Which assumption is more credible depends on the setting. When the worry is that the treated group is simply a different kind of unit to begin with, DiD trades permissiveness on the level dimension for strictness on the trend dimension, and often fits better than matching; when the between-group difference looks more like a level difference driven by observable characteristics, matching's logic is more natural instead.

Parallel trends depends on functional form, and this is not a technical detail but an identification choice to be faced head-on. Suppose A1 holds on log GMV, that is, the two groups' expected increments in the log outcome are equal, then this generally does not hold for the GMV level. The reason is that the between-group increment of $\mathbb{E}[Y]$ depends on the whole distribution and how it shifts over time, whereas $\mathbb{E}[\log Y]$ and $\log \mathbb{E}[Y]$ differ by a quantity that depends on the shape of the distribution, and since the two groups' distributional shapes and evolutions differ, log parallelism does not buy level parallelism. A minimal example to feel this: one city's GMV rises from 100 to 110, another from 200 to 220, and the two are not equal in level gains (10 versus 20), but in logs both grow by about 10%, so level is not parallel while log is. So whether things are parallel depends on which scale you view them on. Roth and Sant'Anna (2023, Econometrica) give a sharp characterization: parallel trends holds for all monotone transformations of the outcome simultaneously if and only if a very strong condition is met, roughly, either treatment is approximately randomly assigned (the distribution of $Y(0)$ is the same across the two groups, however they co-evolve over time), or the distribution of $Y(0)$ does not shift over time within groups at all, that is, strict stationarity. Outside this knife-edge, using log or level is a choice between two different identification assumptions, and at the same time a choice between two different estimands, log corresponding to a percentage effect, level to an absolute-quantity effect. The basis for the choice should be the mechanism and the research question: the commission rate acts proportionally on sellers' marginal incentives, so log has a mechanistic reason, a point Section 7 revisits from the angle of robustness.

Finally, note that A1 does not require the treatment effect to be the same across groups. The platform deliberately cutting commissions first in the cities with the largest expected effect, this selection on gains, does not violate parallel trends at all, because parallel trends constrains the path of $Y(0)$, not the size of the effect. The only price is that the identified quantity is the treated group's own ATT, and extrapolating it to "what if we also cut commissions for the control cities" loses its basis, because the control cities are precisely the ones the platform judged to have a smaller effect and so left for later or did not treat at all. This gap between the ATT and the ATE recurs throughout the staggered setup.

### 3.3 Anatomy of No Anticipation

No anticipation is about a very plain thing: before the treatment actually lands, the party about to be treated should not change its behavior in advance just because it senses the treatment coming. Once there is such a jump-the-gun response, the clean stretch of pre-treatment history is muddied.

A2 written in 2x2 is $Y_{i1} = Y_{i1}(0)$, and generalized to many periods it is: for the treated group in all pre-treatment periods $t < G_i$, $Y_{it} = Y_{it}(0)$. It guarantees that pre-period observations are clean untreated outcomes and thus qualify as a baseline. What breaks it is forward-looking behavior. If the commission cut is announced in advance, sellers may hold back orders to release them in the low-commission period, and pre-period GMV is artificially depressed, so the first difference counts this rebound as effect and overstates. Conversely, if sellers ramp up advertising after the announcement but before it takes effect, to stock up demand for the low-commission period, the pre period is inflated and the effect is understated. Which direction depends on the specific form of the behavior, but the common point is that the pre period is no longer a world undisturbed by the policy.

No anticipation is directly tied to the choice of base period. Identification uses the last pre-period as the baseline precisely because that period is not yet contaminated by treatment. If you suspect an anticipatory response $\delta$ periods ahead, the remedy is to move the base period earlier to $g - \delta - 1$, using only earlier periods you are confident are clean as the baseline, at the cost of sacrificing some sample and precision. The observable footprint it leaves is an abnormal climb or dip in the last few periods where $e$ is negative on the event-study plot, which is why, when drawing an event study, you must watch the few pre coefficients right next to the treatment date closely.

### 3.4 SUTVA, Spillovers, and Composition

SUTVA sounds abstract, but landed on a platform it is quite plain: a city's outcome should be determined only by whether it itself was treated, and should not be pulled around by whether other cities were treated. Once it fails, the control group that should have been clean is contaminated by the spillover of treatment.

A3 is usually written as SUTVA, and it unpacks into two conditions. One is no interference, that a unit's outcome depends only on its own treatment status and is unaffected by whether others are treated; the other is that treatment has a single version, that there is no mixing under the name "commission cut" of several treatments of different intensities. The most problematic for DiD is no interference, because once it fails, the control group's observed outcome is no longer a clean $Y(0)$, and the $\gamma(0)$ in the earlier derivation, which should have served as the counterfactual trend, has treatment spillover mixed into it.

In platform settings the channels of spillover are very concrete and worth naming one by one, because each corresponds to a different diagnostic. Sellers can move their registration to a commission-cut city, and treatment moves directly onto the control city; buyers can order across cities, and the low prices in treated cities siphon demand from control cities; the platform's on-site traffic allocation algorithm may pull impressions away from control cities because treated cities' conversion rates rise, a spillover mediated purely by the algorithm. If any one holds, the control group's change has the shadow of the treatment effect mixed in.

Behind spillovers stands a more fundamental distinction: DiD estimates a partial-equilibrium contrast. It assumes by default that the control group provides the counterfactual of what the treated group would have been absent treatment, and this default holds only when treatment does not change the equilibrium the control group sits in. If the commission cut is not a local shock confined to one city but changes the whole platform's seller pricing equilibrium, buyer traffic allocation, or competitive landscape, then the control cities are themselves swept into the new equilibrium, and they are no longer a clean counterfactual for any unit, so the DiD framework fails wholesale, and the question has to be answered with a structural model or a platform-level experiment. Judging which case you are in relies on economic judgment about the treatment's transmission mechanism, not on any statistical test.

Adjacent to but different in nature from spillover is the composition problem. Even if treatment does not spill over onto the control group's outcome, as long as it changes the composition of the sample, the before-after comparison is no longer of the same set of units. If the commission cut attracts sellers from neighboring cities to move in, the treated cities' seller population in the post period is not the same set as in the pre period, and the observed change in city-month average GMV has compositional change mixed in. Redefining the population as sellers who stay in the city throughout the entire window can restore a clean comparison, but that is already a different estimand, answering the effect on stably retained sellers rather than the effect on all sellers.

By this point the division of labor by which the three assumptions each block one route of bias is clear: A1 blocks the counterfactual trend difference, A2 blocks the pre period from being contaminated by anticipation, A3 blocks the control group from being contaminated by spillover, plus a composition checkpoint shoulder to shoulder with A3 that blocks change in sample composition. They fail in different ways, and diagnosis and remedy each go their own way, with Section 7 giving actionable responses one by one.

### 3.5 Conditional Parallel Trends

Unconditional parallel trends is sometimes not credible, but conditioning on pre-treatment covariates makes it more so. The typical case is where the growth path of the outcome depends systematically on certain city characteristics (population structure, logistics density, existing penetration), and these characteristics are imbalanced between the treated and control groups, so the unconditional between-group trends should not be equal to begin with, but should be equal between cities of the same type.

::: {.assumption}
**A1' (conditional parallel trends)** Conditional on pre-treatment covariates $X$,

$$\mathbb{E}\big[Y_{i2}(0) - Y_{i1}(0) \mid X_i, D_i\big] = \mathbb{E}\big[Y_{i2}(0) - Y_{i1}(0) \mid X_i\big]$$
:::

The correct reading of A1' is to apply matching's assumption to a redefined outcome variable, namely the growth of the outcome rather than its level. This step opens up matching's whole toolbox: you can let the time-invariant $X_i$ carry time-varying coefficients in the regression ($\beta_t' X_i$, equivalent to allowing each type of city its own trend), you can weight outcome growth by the propensity score, or you can combine the two into a doubly robust estimator, with Section 4.6 detailing the implementation. Two disciplines must be kept here. First, $X$ must take pre-treatment values, and putting a variable that will respond to treatment into the conditioning set is the bad control error of Section 7, which treats part of the effect itself as a difference to be controlled away. Second, and more fundamental, conditioning can only repair trend differences driven by observable characteristics; if the treated and control groups differ in the unobservable drivers of growth, no number of covariates can help, and this is already a design-level problem, not something a regression can remedy.

### 3.6 From 2x2 to Staggered: Identifying $ATT(g,t)$

The 2x2 identification logic carries over directly to the staggered setup, and this is the key step to understanding all the modern estimators that follow. For some cohort $g$ and any period $t \geq g$ after it is treated, replace the treated group in 2x2 with cohort $g$ and the control group with a valid comparison group still untreated up to period $t$, and the same double difference identifies $ATT(g,t)$: take $g - 1$, the last period before cohort $g$ is treated, as the baseline, and subtract the comparison group's contemporaneous change from cohort $g$'s change from $g-1$ to $t$. The assumptions supporting this step are the natural generalizations of the three in 2x2: parallel trends requires cohort $g$ and the comparison group to have equal expected growth in $Y(0)$, no anticipation qualifies $Y_{g-1}$ as a baseline, and SUTVA keeps the comparison group clean.

The crux is whom to choose as a valid comparison group. There are two valid choices: never-treated units that are never treated, and not-yet-treated units untreated up to period $t$ (they will eventually be treated, but not yet). The two choices correspond to two versions of the parallel trends assumption, whose credibility varies by setting, and Section 4.5 writes both versions out formally and discusses the tradeoff. The core idea to establish here is that identifying $ATT(g,t)$ in a staggered setup never needs, and should never, use a unit already treated in period $t$ as a comparison, because that unit carries a treatment effect and its change is not a clean $Y(0)$ trend. It is exactly this discipline that the static TWFE regression, dissected in the next section, quietly violates, which is the master root of its failure.

### 3.7 Testability and Pre-Trends

The last part of this section discusses testability, which is where things are most easily misused. A1 is an assertion about the counterfactual, about how the treated group's average outcome would move had it not been treated, and the counterfactual is never observed, so parallel trends is fundamentally untestable within the DiD framework, just as exogeneity is untestable within the OLS framework. The common pre-trends test tests something else: whether, before treatment actually happens, the two groups' observable outcomes are parallel. This is a necessary but not sufficient observable shadow of parallel trends: if the pre period is already non-parallel that should of course make you nervous, but pre-period parallelism does not imply counterfactual parallelism after the treatment date.

This failure to imply is not nitpicking. Imagine a platform version of the scenario: Kestrel cuts commissions in a city because headquarters judges that a strong competitor will enter that city in half a year, and the cut is only one move in a package of responses that also includes subsidies and traffic tilting. Then this city and the comparison city will part ways from the moment the competitor enters, whether or not commissions are cut, and none of this shows any sign in the pre period. Once one thing is destined to start diverging after treatment, chances are other things are diverging at the same time. The pre-trends test has no defense against this class of worries about why treatment happens exactly here and now, because it can only see trend differences that have already occurred, not those about to occur.

::: {.warning}
Flat pre-trends are supporting evidence, not proof of parallel trends, and they are low-power supporting evidence: a trend difference the test cannot detect may be enough to flip the conclusion (Roth 2022, detailed in Section 7). Treating "the pre-trends test passed, therefore parallel trends holds" as a formal argument for identification does not stand.
:::

The key points of this section can be summarized as follows: the identification of 2x2 DiD is a precise division of labor among three assumptions each blocking one route of bias, A1 governs the counterfactual trend, A2 governs a clean pre period, A3 governs a clean control group, the conditional version A1' extends the design to trend parallelism between cities of the same type, and the double difference of 2x2 in turn generalizes to identifying $ATT(g,t)$ under staggering. All these assumptions are about the counterfactual, none directly testable, and all that can be tested is their observable shadow cast on the pre period, and between whether the shadow is parallel and whether the substance is parallel lies a layer of inference the data cannot fill.

## 4 Estimation: From 2x2 to Staggered

This section proceeds in order of failure: the method in each subsection suffices in the previous subsection's scenario, and fails when moved to a more complex one, which then gives rise to the next method.

### 4.1 The Regression Implementation of 2x2

The sample analog estimator of 2x2 is the DiD of four means. It can be written equivalently as a regression:

$$Y_{it} = \alpha + \gamma\, D_i + \lambda\, \mathrm{Post}_t + \beta\, (D_i \times \mathrm{Post}_t) + \varepsilon_{it}$$

The OLS estimate of the interaction coefficient $\hat\beta$ exactly equals the sample DiD, $\hat\beta = \overline{\Delta Y}_T - \overline{\Delta Y}_C$. This is a numerical identity, not dependent on any distributional assumption, and you can verify it by plugging into the formula yourself. The advantage of the regression form is operational convenience: adding covariates, computing clustered standard errors, and generalizing to many periods and groups are all easy. The implementation of Card and Krueger (1994) is exactly this equation plus store characteristics, and they also did a clever variant, using each store's depth of exposure to the minimum wage law (the gap between the initial wage and the new minimum) to construct a continuous treatment intensity, turning the between-group comparison into a dose-response relationship, detailed in Section 5.1.

### 4.2 The Event Study Setup

Policy effects rarely arrive in one step. A commission cut first changes sellers' marginal profit, and only then transmits to stocking, advertising, and pricing, so the effect accumulates over time. A two-period design cannot see any of this, and the standard practice with multi-period data is to replace the static $\beta$ with a set of event-time coefficients, that is, the event study (dynamic) setup:

$$Y_{it} = \alpha_i + \lambda_t + \sum_{\substack{e = -K \\ e \neq -1}}^{L} \theta_e\, \mathbf{1}\{t - g_i = e\} + \varepsilon_{it}$$

Here and in what follows, in §4.5 and §6.5 through §6.8, the $\theta$ follows the original notation of Callaway-Sant'Anna and Sun-Abraham, denoting event-study coefficients and aggregated treatment effects, and is unrelated to the $\theta$ used as a structural parameter elsewhere in this book. There are three syntactic details, each hiding a statistical reason. Base-period normalization: for units that will eventually be treated, the full set of event-time dummies sums to 1 within the sample window and is perfectly collinear with that unit's fixed effect, so one must be dropped; the convention drops $e = -1$ (the last pre-treatment period), so each $\theta_e$ reads as "the change relative to one period before treatment," and $\theta_{-1} \equiv 0$ is not an estimated zero but a zero by definition. Endpoint handling: leads and lags far from the treatment date are identified by fewer and fewer cohorts (only the earliest-treated cohort can contribute an $e = 20$ observation), and the standard practice is to bin the two ends ($\mathbf{1}\{t - g_i \leq -K\}$ and $\mathbf{1}\{t - g_i \geq L\}$ each as a dummy) or to trim outright; note that if you quietly omit certain event times without binning, those observations become implicit controls, and the interpretation of all coefficients changes. When there are no never-treated units, there is a second collinearity among the event dummies, time fixed effects, and unit fixed effects, requiring one more coefficient to be normalized (a common practice is to also set the most distant lead to zero), and only relative dynamics are then identified.

::: {.warning}
The dynamic equation with $\sum_{e=-K}^{L}$ written out in full, omitting no base period, is perfectly collinear. The software will not error; it will silently drop a collinear column for you, and which one it drops depends on the column ordering, so your coefficients are henceforth relative to a benchmark you yourself do not know. Always specify the ref period explicitly.
:::

When reading an event-study plot, you can check the following in order:

1. Whether the pre-period coefficients are close to zero with no systematic climb or dip. This is supporting evidence, but note its power is limited and it is no proof.
2. Whether the shape of the response at $e = 0$ is consistent with the mechanism's timetable: if the commission takes effect at the start of the month, there should be a reaction that very month, and if the effect starts from $e = -2$, suspect anticipation.
3. The shape of the dynamic path: immediate, gradually accumulating, or overshooting then falling back, and the shape itself is information about the mechanism.
4. Whether the long run enters a plateau or keeps drifting.
5. Whether the confidence intervals of the endpoint coefficients widen noticeably, and whether the cohorts identifying them have become too few to be representative.

### 4.3 How Static TWFE Fails Under Staggering

Now return to the number that did not add up at the start. With three batches of treatment dates, the static TWFE regression looks harmless, but in fact it packages all the 2x2 comparisons the sample can make into one weighted average. There are three types of comparison: early-treated versus never-treated, late-treated versus never-treated, and early and late treated groups serving as controls for each other. The last type hides the most dangerous kind: the late-treated group's before-after change, minus the early-treated group's contemporaneous change (the early group already treated). Here the early-treated group plays the "control group," yet it still carries a treatment effect. If the early group's effect has already settled into a constant, the differencing treats it as a level difference and cancels it, harmless; if the effect is still growing with event time (as it is in Kestrel), the growing part is subtracted off as "the control group's trend," and the late-treated group's estimate is pulled down. This class of comparison is the so-called **forbidden comparison**: already-treated units mixed into the control group.

Goodman-Bacon (2021) turned this intuition into precise algebra.

::: {.theorem}
(Goodman-Bacon decomposition) Let timing group $k$ be treated earlier and $\ell$ later ($\ell > k$), and $U$ the never-treated group; let $\mathrm{PRE}(k)$, $\mathrm{MID}(k,\ell)$, $\mathrm{POST}(k)$ be the sample windows before $k$ is treated, between the two treatment dates, and after $k$ is treated, and $\bar y$ the group-by-window mean. Then the static TWFE estimator decomposes exactly into

$$\hat\beta_{TWFE} = \sum_{k \neq U} s_{kU}\, \hat\tau^{2\times2}_{kU} + \sum_{k \neq U} \sum_{\ell > k} \Big[ s^{k}_{k\ell}\, \hat\tau^{2\times2,k}_{k\ell} + s^{\ell}_{k\ell}\, \hat\tau^{2\times2,\ell}_{k\ell} \Big]$$

where the three types of 2x2 component are, respectively,

$$\hat\tau^{2\times2}_{kU} = \big(\bar y_k^{\mathrm{POST}(k)} - \bar y_k^{\mathrm{PRE}(k)}\big) - \big(\bar y_U^{\mathrm{POST}(k)} - \bar y_U^{\mathrm{PRE}(k)}\big)$$

$$\hat\tau^{2\times2,k}_{k\ell} = \big(\bar y_k^{\mathrm{MID}(k,\ell)} - \bar y_k^{\mathrm{PRE}(k)}\big) - \big(\bar y_\ell^{\mathrm{MID}(k,\ell)} - \bar y_\ell^{\mathrm{PRE}(k)}\big)$$

$$\hat\tau^{2\times2,\ell}_{k\ell} = \big(\bar y_\ell^{\mathrm{POST}(\ell)} - \bar y_\ell^{\mathrm{MID}(k,\ell)}\big) - \big(\bar y_k^{\mathrm{POST}(\ell)} - \bar y_k^{\mathrm{MID}(k,\ell)}\big)$$

The weights $s \geq 0$ and $\sum s = 1$.
:::

The first type is the clean comparison of a treated group against never-treated. The second uses the late-treated group's still-untreated stretch as the control, also clean. The third, $\hat\tau^{2\times2,\ell}_{k\ell}$, is the forbidden comparison: the late-treated group $\ell$'s post-treatment change, minus the already-treated group $k$'s contemporaneous change. The specific form of the weights (write $n_k$ for the group's sample share, $n_{k\ell} = n_k / (n_k + n_\ell)$, $\bar D_k$ for the fraction of periods group $k$ is in the treated state, and $\hat V^D$ a normalizing constant) is

$$s_{kU} = \frac{(n_k + n_U)^2\, n_{kU}(1 - n_{kU})\, \bar D_k (1 - \bar D_k)}{\hat V^{D}}$$

$$s^{k}_{k\ell} = \frac{\big((n_k + n_\ell)(1 - \bar D_\ell)\big)^2\, n_{k\ell}(1 - n_{k\ell})\, \frac{\bar D_k - \bar D_\ell}{1 - \bar D_\ell} \cdot \frac{1 - \bar D_k}{1 - \bar D_\ell}}{\hat V^{D}}$$

$$s^{\ell}_{k\ell} = \frac{\big((n_k + n_\ell)\bar D_k\big)^2\, n_{k\ell}(1 - n_{k\ell})\, \frac{\bar D_\ell}{\bar D_k} \cdot \frac{\bar D_k - \bar D_\ell}{\bar D_k}}{\hat V^{D}}$$

The formulas need not be memorized, but two of their meanings must be grasped. The weights are determined by group size and the variance of treatment status, and the closer $\bar D_k$ is to 1/2 (that is, the closer the treatment date is to the middle of the panel) the larger the weight, which is to say OLS allocates weight by "convenience of identification" rather than by "policy relevance," and there is no economic reason to endorse this weighting. The decomposition weights are all nonnegative, so the negative-weight problem is not in the 2x2 estimates but in the treatment effects: substituting potential outcomes into each 2x2 component and taking the probability limit, one can rearrange into

$$\mathrm{plim}\; \hat\beta_{TWFE} = VWATT + VWCT - \Delta ATT$$

$VWATT$ is the variance-weighted average of the target effects of each 2x2; $VWCT$ is the variance-weighted violation of common trends, zero when parallel trends holds; $\Delta ATT$ collects the within-group change of the treatment effect over time, entering with a negative sign via the forbidden comparisons. When the effect grows with event time $\Delta ATT > 0$, and the static estimate is systematically pulled down, down to far less than any true effect, or even to a sign flip: all true effects positive while $\hat\beta_{TWFE} < 0$ is entirely possible. In Kestrel it pulled 0.088 down to 0.023, and Section 6.4 reproduces this arithmetic with numbers.

The positioning of the Goodman-Bacon decomposition must be accurate: it is a diagnostic tool, telling you what TWFE is actually estimating and where each weight sits, but it does not provide a replacement estimator. After the diagnosis you still need to switch tools, which is the topic of Section 4.5.

::: {.warning}
Do not take "TWFE only mildly understates under staggering" as a general rule. The direction and size of the bias are jointly determined by the effect dynamics and the weight structure: increasing effects bias downward (possibly to a sign flip), decreasing effects bias upward, and instantaneous homogeneous effects are unbiased, because then each 2x2 component is unbiasedly estimating the same constant, and any nonnegative weights summing to one recover the truth. Any intuition that "the TWFE result just needs a small correction" is wrong; you must do the decomposition or switch estimators outright.
:::

### 4.4 The Contamination Problem of Dynamic TWFE

Since the static regression has a flaw, many people's first reaction is "then surely I can just run a dynamic event study." Sun and Abraham (2021) prove the dynamic version has a flaw of its own. Define the cohort-by-event-time effect

$$CATT_{g,e} = \mathbb{E}\big[Y_{i,g+e}(g) - Y_{i,g+e}(0) \mid G_i = g\big]$$

::: {.theorem}
(Sun and Abraham 2021, contamination result) Consider the dynamic TWFE regression of Section 4.2 (some relative periods excluded). The population regression coefficient of the dummy for event time $e$ satisfies

$$\theta_e = \sum_{e' \in \mathcal{E}} \sum_{g} \omega_{e}^{g,e'}\, CATT_{g,e'}$$

The linear combination runs over all cohorts $g$ and all relative periods $e'$, including $e' \neq e$ and the excluded periods. The weights $\omega_e^{g,e'}$ are determined solely by the design structure of the treatment dates (and can be estimated from the data), and satisfy: the own-period ($e' = e$) weights sum to 1 across cohorts; the weights of other included periods sum to 0 across cohorts; the weights of excluded periods generally do not sum to zero. The weights can be negative.
:::

This is a purely mechanical OLS result, and it happens without any identification assumption failing. Three consequences are worth stating separately. If $CATT_{g,e}$ is homogeneous across cohorts (heterogeneity only in the event-time dimension), the cross terms cancel in the summation and $\theta_e$ is clean; the premise of contamination is cross-cohort heterogeneity. Once there is cross-cohort heterogeneity, the coefficients for $e \geq 0$ have effects from other periods mixed in, and more insidiously the pre-period coefficients for $e < 0$ absorb the post-period effects of other cohorts, so even when parallel trends and no anticipation both hold exactly, the pre coefficients may be significantly nonzero. The pre-trends test thereby fails at both ends: it may raise an alarm when there is no problem (false positive), and it may look flat because the contamination term offsets a real trend difference (false negative). The weights being able to be negative means $\theta_e$ can fall outside the convex hull of all $CATT_{g,e}$, estimating an effect that no subgroup possesses.

In the Kestrel simulation this contamination is numerically mild (Section 6.3 will show the naive event study sticking surprisingly close to the truth), because the DGP's cross-cohort heterogeneity is small and there is a never-treated group as backstop. This reminds us to be precise in wording: contamination is a theoretical risk determined by the design structure, and its severity varies by data, so one mild simulation does not mean it is always mild, and its theoretical existence does not warrant declaring all TWFE event studies wildly wrong.

### 4.5 The Three Modern Estimators

The overall idea of the repair is consistent: break the estimand down to the finest $ATT(g,t)$ or $CATT_{g,e}$, use only clean comparisons for each cell, then aggregate with weights the researcher chooses and can see. The three mainstream implementations differ in how the cells are estimated, whom the control group is, and how strong the implicit parallel trends is.

Start with Callaway and Sant'Anna (2021). It takes the $ATT(g,t)$ of Section 2 directly as the estimand. There are two valid choices of control group, corresponding to two versions of the parallel trends assumption.

::: {.assumption}
**A4 (parallel trends, based on never-treated)** For all $t \geq g$:

$$\mathbb{E}\big[Y_t(0) - Y_{t-1}(0) \mid G = g\big] = \mathbb{E}\big[Y_t(0) - Y_{t-1}(0) \mid C = 1\big]$$

where $C = 1$ denotes never-treated.

**A4' (parallel trends, based on not-yet-treated)** For all $t \geq g$, $s \geq t$:

$$\mathbb{E}\big[Y_t(0) - Y_{t-1}(0) \mid G = g\big] = \mathbb{E}\big[Y_t(0) - Y_{t-1}(0) \mid D_s = 0,\, G \neq g\big]$$
:::

The two versions give, respectively, the identification results

$$ATT(g,t) = \mathbb{E}\big[Y_t - Y_{g-1} \mid G = g\big] - \mathbb{E}\big[Y_t - Y_{g-1} \mid C = 1\big]$$

$$ATT(g,t) = \mathbb{E}\big[Y_t - Y_{g-1} \mid G = g\big] - \mathbb{E}\big[Y_t - Y_{g-1} \mid D_t = 0,\, G \neq g\big]$$

The base period is $g - 1$, the last period before cohort $g$ is treated. No anticipation ($Y_{it} = Y_{it}(0)$ when $t < G_i$) is the assumption that qualifies $Y_{g-1}$ as a baseline, and the CS framework also allows limited anticipation: if you suspect a reaction $\delta$ periods ahead, move the base period to $g - \delta - 1$.

::: {.warning}
CS's control group can be never-treated or not-yet-treated (the did package's control_group = "notyettreated"), and there is no restriction that "only never-treated can be used." The correct meaning of the overlap condition is accordingly: each $(g,t)$ cell has a valid control group (never-treated or a not-yet-treated cohort), not that every cohort must be paired with never-treated units. When the never-treated group is very small or suspected of being selected (for example, the 8 cities the platform never touches may have been marginal markets to begin with), the not-yet-treated version is often more credible, because the comparison cities are also cities that will eventually be selected, just later.
:::

After estimating the $ATT(g,t)$ matrix comes aggregation. Aggregate by cohort then take a grand average:

$$\theta_{S}(g) = \frac{1}{\mathcal{T} - g + 1} \sum_{t = 2}^{\mathcal{T}} \mathbf{1}\{g \leq t\}\, ATT(g,t), \qquad \theta_{S}^{O} = \sum_{g = 2}^{\mathcal{T}} \theta_{S}(g)\, P(G = g \mid G \leq \mathcal{T})$$

Aggregate by event time into an event study:

$$\theta_{D}(e) = \sum_{g = 2}^{\mathcal{T}} \mathbf{1}\{g + e \leq \mathcal{T}\}\, ATT(g, g + e)\, P(G = g \mid G + e \leq \mathcal{T})$$

A practical trap: in the did package, aggte(type = "simple") is yet another weighting (weighted by the number of post-treatment observation cells, so cohorts with long exposure get more weight), and it is generally not equal to the $\theta_S^O$ given by type = "group," so the paper must state clearly which one is reported. When the cohorts are few, the most honest thing to do is actually to report the $ATT(g,t)$ matrix itself, letting readers see the effect for every cohort at every period, with the aggregate numbers in back.

The interaction-weighted (IW) estimator of Sun and Abraham (2021) is the regression-form implementation of the same philosophy, constructed in three steps. First, run a saturated interaction regression, excluding $e = -1$, using never-treated (or, absent it, the latest-treated cohort, with its already-treated periods removed from the sample) as a pure control:

$$Y_{it} = \alpha_i + \lambda_t + \sum_{g} \sum_{e \neq -1} \delta_{g,e}\, \mathbf{1}\{G_i = g\}\, \mathbf{1}\{t - g = e\} + \varepsilon_{it}$$

Because cohort-by-event-time is put in saturatedly, the cross-cohort contamination of Section 4.4 cannot occur by construction, and each $\hat\delta_{g,e}$ estimates $CATT_{g,e}$ under parallel trends and no anticipation. Second, estimate each cohort's observable sample share at event time $e$. Third, aggregate by share:

$$\hat\theta_e^{IW} = \sum_{g} \hat\delta_{g,e}\; \hat P\big(G = g \mid G \in \mathcal{G}_e\big), \qquad \mathcal{G}_e = \{g : g + e \in [1, \mathcal{T}]\}$$

The weights are cohort shares, naturally nonnegative, summing to one, and clearly interpretable, exactly the three things the implicit weights of dynamic TWFE cannot do. Using never-treated as the control and $e = -1$ as the common base period, IW and CS's dynamic aggregation are usually very close numerically, differing in the details of base-period handling and share estimation. fixest::sunab runs it in one line.

The imputation estimator of Borusyak, Jaravel and Spiess (2024) takes a different angle: rather than find a control group to difference against, impute the counterfactual directly. Assume for the untreated state $Y_{it}(0) = \alpha_i + \lambda_t + \varepsilon_{it}$, then proceed in three steps. Using only untreated observations (all periods of never-treated cities plus the pre periods of treated cities), estimate $\hat\alpha_i, \hat\lambda_t$ by OLS; fill in a counterfactual $\hat Y_{it}(0) = \hat\alpha_i + \hat\lambda_t$ for each post-treatment observation, giving the unit-level effect $\hat\tau_{it} = Y_{it} - \hat Y_{it}(0)$; finally average with weights the researcher chooses, $\hat\tau_w = \sum_{(i,t) \in \Omega_1} w_{it}\, \hat\tau_{it}$, and wanting an event study, a cohort average, or a grand average is just a change of weights.

Its properties have four points that must be stated honestly. The fixed effects are fit with all pre periods, so the parallel trends it implies covers the entire pre period, not just CS's local comparison from base period to target period, which is a strictly stronger assumption, and in exchange comes efficiency: BJS prove it is the most efficient linear unbiased estimator under homoskedastic spherical errors. Treated observations never enter the fitting regression, so the forbidden comparison does not exist by construction. The pre-trends test must be done separately (regress the residuals of untreated observations on lead dummies), because the estimation procedure mechanically flattens the pre-period fit, and looking directly at the residual plot would fool yourself. Unit-level effect heterogeneity is folded into the error term, so inference is conservative for the sample ATT. The R implementation is didimputation::did_imputation; the same idea also underlies Gardner's (2022) two-stage DiD and Liu, Wang and Xu's fect.

One more method is worth a brief mention. The $DID_M$ estimator of de Chaisemartin and D'Haultfœuille (2020) compares switchers (units that change treatment status) between adjacent periods against non-switchers, targeting the switchers' instantaneous ATT in the period of the switch. It is not the first choice for most of this chapter's scenarios, but it has one special ability: it allows treatment to exit. The CS framework rests on the irreversibility assumption $D_{it} \geq D_{i,t-1}$, and a platform setting where the commission rate is repeatedly adjusted up and down directly violates it, and in that case $DID_M$ (the R and Stata packages DIDmultiplegt / did_multiplegt) is one of the few choices that still stands.

The comparison of the three main estimators is summarized in the table below, which is the one place in this chapter where a table is well suited for a method comparison.

| Dimension | Callaway-Sant'Anna | Sun-Abraham IW | BJS imputation |
|---|---|---|---|
| parallel trends strength | local: base period $g-1$ against each post period, period by period | same tier as CS ($e = -1$ base period) | global: the entire pre period must be parallel, strongest |
| control group | never-treated or not-yet-treated, optional | never-treated, or latest-treated cohort | all currently-untreated observations (both types used) |
| covariates | outcome regression / IPW / doubly robust | linearly into the fixed effects structure | linearly into the imputation regression |
| anticipation | built-in anticipation parameter, base period can move | manually moving base period partly handles it | needs explicit removal of suspect periods from the fitting sample |
| efficiency | robustness first, efficiency second | robustness first | most efficient under homoskedasticity |
| under parallel trends violation | lower sensitivity (only local comparisons) | lower | most sensitive (relies on long pre-period extrapolation) |
| R implementation | did::att_gt plus aggte | fixest::sunab | didimputation::did_imputation |

### 4.6 Estimation with Covariates

Under the conditional assumption A1', the three implementation routes correspond to which nuisance model you trust. The outcome regression route models the control group's outcome growth, $m_{g,t}(X) = \mathbb{E}[Y_t - Y_{g-1} \mid X, C = 1]$, and uses it to predict the treated group's counterfactual growth. The IPW route models "being selected into cohort $g$," using a generalized propensity score $p_g(X) = P(G_g = 1 \mid X, G_g + C = 1)$ to reweight the control group to the treated group's $X$ distribution. Doubly robust combines the two, and the DR estimand of Sant'Anna and Zhao (2020) (the never-treated control version, $G_g = \mathbf{1}\{G_i = g\}$) is

$$ATT(g,t) = \mathbb{E}\left[\left(\frac{G_g}{\mathbb{E}[G_g]} - \frac{\dfrac{p_g(X)\, C}{1 - p_g(X)}}{\mathbb{E}\left[\dfrac{p_g(X)\, C}{1 - p_g(X)}\right]}\right)\big(Y_t - Y_{g-1} - m_{g,t}(X)\big)\right]$$

As long as one of $p_g$ and $m_{g,t}$ is correctly specified, the estimate is consistent; setting $m_{g,t}$ to zero degenerates to pure IPW, and removing the weighting term degenerates to pure outcome regression. The did package's default est_method = "dr" applies this formula cell by cell to each $(g,t)$, and the not-yet-treated version merely swaps the comparison set. Two disciplines: covariates can only take pre-treatment values; conditioning changes the credibility of the parallel trends assumption, not the algebra of the aggregation, so do not expect that stuffing in a few more control variables can repair a design-level selection problem.

### 4.7 Inference

The point-estimate problem of DiD is now covered, and the standard-error problem can ruin a paper just as well. The logic of clustering is one sentence: which observations do you believe are independent. The cluster structure assumes groups are independent of one another and within-group correlation is arbitrary, and the asymptotic theory runs on the number of groups $G$ rather than the number of observations $N$. The CR1 sandwich estimator is

$$\hat V^{CR1}_{\hat\beta} = \frac{G}{G-1}\,\frac{N-1}{N-K}\,(X'X)^{-1}\left(\sum_{c=1}^{G} X_c' \hat u_c \hat u_c' X_c\right)(X'X)^{-1}$$

where $K$ is the number of regression coefficients. The scalar in front is the finite-sample correction, and dropping it gives the uncorrected CR0 sandwich, and with few clusters the difference between the two is not negligible.

Cluster at which level? At the level of treatment assignment. Kestrel's commission policy switches on and off by city, and the city is the independent unit in the design, so cluster at the city; clustering at the seller or city-month level would severely understate the standard errors, because it pretends that shocks in the same city in different months are mutually independent, and the AR(1) noise plainly denies that. Whatever level you put fixed effects at, you usually need to cluster at least at that level. The cost of choosing the wrong level is not a few percent, and standard errors can differ by an order of magnitude.

The classic source of this discipline is Bertrand, Duflo and Mullainathan (2004). Their observation is that DiD's typical data (a state-by-year wage panel) is strongly serially correlated, the treatment variable is itself serially correlated (once switched on it stays on), and papers at the time commonly used unclustered standard errors. Their diagnostic method is elegant and brutal: assign placebo laws at random in the CPS wage data (random states, random effective years, true effect identically zero) and count the frequency at which the unclustered test rejects the null. The answer is that at a 5% nominal level about 45% of the placebos "found an effect." They also verified the effective remedies: cluster at the state (the treatment assignment level), do a block bootstrap by state, or compress the time series to one pre mean and one post mean per state and do a 2x2. The last move also gives a nice intuition as a bonus: the effective number of observations of DiD is the number of clusters, not the number of panel cells.

When the number of clusters is small (the rule-of-thumb threshold is roughly $G < 40$ to 50, worse when the treated clusters are very few), CR1 with $t_{G-1}$ critical values over-rejects. The treatment escalates by severity: CR2 (the Bell-McCaffrey small-sample correction plus Satterthwaite degrees of freedom); the wild cluster bootstrap-t (Cameron, Gelbach and Miller 2008, the R package fwildclusterboot); and when the treated clusters number in the single digits even the wild bootstrap is unreliable, and the honest fallback is the randomization inference logic of Conley and Taber (2011) and MacKinnon-Webb. Kestrel has 60 cities and three cohorts, so the baseline CR1 already suffices; if your application has only a dozen cities, CR2 and the wild bootstrap belong in the main table, not the appendix.

The roadmap of this section is now complete: the 2x2 regression is equivalent to DiD; multi-period dynamics are expressed with an event study; under staggering the weights of static TWFE are exposed by Goodman-Bacon, and dynamic TWFE by Sun-Abraham; CS, SA-IW, and BJS repair with clean comparisons plus transparent weights, each paying a different strength of parallel trends; covariates go doubly robust; inference clusters at the treatment assignment level, with an emergency kit for few clusters.

## 5 Anchoring Papers

A method only holds up when it lands in real research. Two anchoring papers, one the historical origin of the 2x2 method, one a staggered application in the IS field, each combed through by five elements (paper, method, data, results, limitations), with the focus on how the assumptions are defended.

### 5.1 Card and Krueger (1994, AER)

::: {.case}
The paper and its place in method history: "Minimum Wages and Employment: A Case Study of the Fast-Food Industry in New Jersey and Pennsylvania," American Economic Review. This paper turned DiD from a trick into the default design language for policy evaluation, and is also one of the emblematic texts of the later so-called credibility revolution: the identification argument shifted from "my model specification is right" to "where does my policy variation come from."

Method: textbook-grade 2x2. New Jersey raised the state minimum wage from $4.25 to $5.05 on April 1, 1992, while eastern Pennsylvania next door held constant. The treated group is NJ fast-food stores, the control group is eastern PA fast-food stores, and employment is measured once before and once after. Beyond the between-group comparison, the paper also exploits within-NJ variation in treatment intensity: stores whose initial wage was already above 5.05 are effectively unconstrained, and the depth of the bite (the wage gap) forms a continuous treatment variable, so the dose-response relationship gives the between-group comparison an independent corroboration.

Data: a telephone survey of about 410 fast-food chain stores (Burger King, KFC, Wendy's, Roy Rogers), and the analysis sample with usable data in both waves is slightly smaller than this, with the first wave in February to March 1992 before the raise and the second in November to December after. The reason for using the fast-food industry is itself a design argument: the minimum wage bites hardest here, and the jobs are homogeneous and easy to measure.

Results: contrary to the textbook prediction of a competitive labor market, employment in NJ did not fall relative to PA, and the DiD point estimate is positive (full-time-equivalent employment rose by about 2.8 jobs per store relative to PA, with limited precision), and the continuous wage-gap treatment version gives results in the same direction.

Limitations and aftermath: the telephone-survey measurement of employment is very noisy, and Neumark and Wascher re-estimated with store payroll data and got the opposite sign, sparking a head-on exchange in the AER in 2000. Card and Krueger (2000) re-examined using the BLS administrative payroll data and maintained the conclusion, wording it that this raise "probably had no effect on total employment ... possibly had a small positive effect," and this self-correction that narrows the claim in the face of evidence is worth learning. The design-level limitations are clearer with today's eyes: there is essentially only one treated state and one control state, so inference faces the extreme few-cluster problem of $G = 2$ (the tools of Section 4.7 did not exist then); parallel trends is untestable, and the macro-cycle difference between the two states around 1992 is a real threat; composition is also a hidden risk, since cross-border flows like NJ stores closing and PA stores opening would contaminate both groups at once.
:::

The defense strategy of this paper is worth emulating: the neighboring state in the same industry presses down the common shock, the tight before-after timing around the policy date presses down the long-run trend, and the within-state continuous treatment intensity does an internal cross-check, and every step buys insurance for the credibility of the counterfactual, rather than piling control variables into the regression.

### 5.2 Greenwood and Wattal (2017, MISQ)

::: {.case}
The paper: "Show Me the Way to Go Home: An Empirical Investigation of Ride-Sharing and Alcohol Related Motor Vehicle Fatalities," MIS Quarterly. The question is whether the entry of ride-sharing reduces drunk-driving fatalities, and it is a representative work of staggered DiD in the IS field.

Method: staggered DiD, with Uber entering California cities in batches, the staggering of entry dates providing the identifying variation, and the regression carrying location and time fixed effects. Note that the paper was published in 2017, before the 2021 estimator revolution, and uses the standard TWFE-class setup of the time.

Data: a panel of multiple California locations, with traffic fatalities from the California Highway Patrol's state accident records (CHP/SWITRS), matched to the entry timetable of Uber's product lines.

Results: after the low-price product UberX enters, alcohol-related motor vehicle fatalities fall significantly, with the main estimates in the 3.6% to 5.6% range; the high-price product Uber Black has no robust effect. The split along the price dimension is itself mechanism evidence: if the decline were merely an artifact of trend or selection, there would be no reason for it to appear only in the low-price product, so the product-line contrast also serves as a placebo.

Limitations: Uber's entry dates are not random, and the decision to enter a city correlates with city size, regulatory environment, and demand growth, and the authors respond head-on with an analysis of entry determinants and pre-trends evidence rather than dodging. By today's standards, staggered TWFE should be re-estimated with CS or SA, and effect heterogeneity (large cities enter first and may have larger effects) is precisely a real-world instance of the negative-weight risk of Section 4.3. Spillover is also real: substitute travel in neighboring areas and drunk drivers crossing borders both dilute the cleanliness of the control group.
:::

Taken together, the significance of the two anchors is clear: Card and Krueger show how a 2x2 design holds up on institutional detail and internal cross-checks; Greenwood and Wattal show how staggered variation is used in platform research and how endogenous entry is defended, and also show how methodological progress makes an old paper's estimation strategy obsolete while its design ideas do not expire.

## 6 A Full Walkthrough on Kestrel Data

Now run the full toolkit of Section 4 through the Kestrel panel end to end. The code below uses R 4.5.3, with set.seed(11) fixed for reproducibility, and every number quoted in the text comes from the actual run output of this code.

### 6.1 The DGP

The design parameters are as follows: 60 cities over 36 months, 2160 city-month observations in total; three cohorts cut commissions in months 13, 19, and 25, of size 17, 17, and 18 cities, with 8 cities never-treated.

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

The effect is a piecewise-linear ramp in event time: a slow start ($e = 0$ to 5 from 0.015 to 0.020), then rapid growth to the plateau of 0.140 at $e = 12$, then a further drift of 0.004 per month. The cohort boost makes late cohorts have a slightly larger effect at the same event time, and the city-level multiplier creates between-city heterogeneity. The assembly of the outcome variable and the recording of the true effect:

```r
dt[, event_time := ifelse(g > 0, month - g, NA_integer_)]
dt[, D := as.integer(g > 0 & month >= g)]
dt[, tau_true := 0]
dt[D == 1, tau_true := (tau_profile(event_time) +
                          cohort_boost[match(g, cohort_g)]) * tau_mult]
dt[, lgmv := 8 + city_fe + month_fe + tau_true + eps]
```

The city fixed effects carry level selection (early-treated cities sit at a higher level), the month fixed effects contain linear growth plus a sinusoidal seasonal term, and the noise is within-city AR(1) ($\rho = 0.35$). Parallel trends holds exactly by construction, with no anticipation and no spillover: the success or failure of all the estimators that follow is only about the aggregation algebra, not about identification, which is exactly the intent of the experimental design.

The true values computed from realized effects: the simple ATT is 0.088; by cohort, g13 is 0.106, g19 is 0.085, g25 is 0.058. Note the direction: late cohorts have a larger effect at the same event time (the boost is positive), but their exposure is short and does not catch the later half of the ramp, so the realized cohort ATT is smaller instead. This phenomenon of "whose effect is larger" depending on how you ask is yet another reason to fix the estimand before discussing estimation.

There is one design decision to explain. The never-treated group has only 8 cities, which is deliberately kept small. The weights of static TWFE are determined only by the treatment matrix, and if we kept 20 never-treated cities, the clean treated-versus-untreated comparison would eat up the vast majority of the weight, and no matter how steep the effect dynamics, $\hat\beta_{TWFE}$ could not drop below half the truth, and the contrast would not be striking enough. Shrinking the clean control to 8 cities pushes about 61% of the weight onto the timing comparisons, and the teaching effect appears at once. This algebraic fact that "the never-treated share backstops the TWFE bias" is itself worth noting: the same estimator, the same effect structure, a different control-group size, and the bias can differ several fold.

### 6.2 Static TWFE: Reproducing the Opening Number

```r
m_static <- feols(lgmv ~ D | city_id + month, data = dt, cluster = ~city_id)
```

The output is $\hat\beta_{TWFE} = 0.0227$ (SE 0.0087), leaving only 25.8% of the truth 0.088. This is the number of Section 1.

### 6.3 The Naive Event Study

```r
dt[, ev_naive := fifelse(g > 0, pmin(pmax(month - g, -12L), 18L), -1L)]
m_naive <- feols(lgmv ~ i(ev_naive, ref = -1) | city_id + month,
                 data = dt, cluster = ~city_id)
```

The base period is $e = -1$, the endpoints are binned at $e \leq -12$ and $e \geq 18$, and never-treated cities are coded $e = -1$ so that all their dummies are zero and they fall into the reference group. The result is surprisingly respectable: the pre-period coefficients are at most 0.009 ($e = -3$, SE 0.008), all insignificant; the post-period coefficients broadly track the truth, only generally a touch low (the reconciliation table in Section 6.8 shows this systematic mild understatement). This confirms the reminder at the end of Section 4.4: the severity of Sun-Abraham contamination depends on the size of the cross-cohort heterogeneity and the control structure, and this DGP's cohort boost is small and there is a never-treated backstop, so the contamination is numerically mild. One should not conclude from this that dynamic TWFE is harmless, nor conversely claim in your own paper that it must be wildly wrong, and the correct thing to do is to run an SA or CS comparison and see how much they differ.

At the same time note: the static regression is wildly wrong on the same data while the dynamic regression is broadly fine, and this itself is a teaching point. The static coefficient is not simply a weighted average of the dynamic coefficients, and its weight structure (the Bacon decomposition) is an entirely different matter.

### 6.4 Goodman-Bacon Decomposition: Where This 0.023 Comes From

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

The three rows each tell a story. The clean treated-versus-never-treated comparison averages 0.073, close to the truth (it only covers part of each cohort's post window, so it is slightly below 0.088), but gets only 39% of the weight. The earlier-versus-later comparison (using the late cohort's pre period as control) is also clean, and its average 0.026 is low because this class of comparison can only see the effect in the early ramp of the early cohorts. The real damage comes from the third row: 30.4% of the weight sits on forbidden comparisons, whose average estimate is negative 0.046. Looking at individual components is even more intuitive: the 25-versus-13 comparison has weight 0.138 and estimate -0.054; 19 versus 13 has weight 0.098 and estimate -0.036; 25 versus 19 has weight 0.069 and estimate -0.043. The mechanism is exactly the algebra of Section 4.3: g13's effect is still drifting up from the plateau within the window after g25 is treated, and this stretch of growth is subtracted from g25's change as "the control group's trend," subtracting out a string of negative values. All true effects are strictly positive, three tenths of the weight has components estimated as negative, and 0.088 is thereby diluted to 0.023.

### 6.5 Sun-Abraham and the Event-Study Plot

```r
dt[, g_sunab := fifelse(g == 0, 10000L, as.integer(g))]
m_sa <- feols(lgmv ~ sunab(g_sunab, month) | city_id + month,
              data = dt, cluster = ~city_id)
sa_att <- summary(m_sa, agg = "att")
```

sunab requires the never-treated cohort value to fall outside the sample period, so first change 0 to 10000. The aggregated ATT is 0.081 (SE 0.008). fixest's built-in iplot(m_sa) produces a diagnostic plot in one line; the comparison plot below is drawn with ggplot2 in order to overlay the true-effect curve as well:

![Comparison of the naive TWFE event study, Sun-Abraham, and the true dynamic effect, $e \in [-12, 18]$, with 95% confidence intervals.](assets/fig/fig_13_eventstudy.svg)

A technical footnote: sunab's full event-study coefficient covariance matrix triggered a fixest warning (VCOV not positive semidefinite, automatically corrected), so the SA standard errors in the plot are the corrected ones, and the ATT aggregation is unaffected. When encountering such numerical warnings, one should state them honestly in the results rather than silently ignore them.

### 6.6 Callaway-Sant'Anna

```r
set.seed(11)
cs <- att_gt(yname = "lgmv", tname = "month", idname = "city_id",
             gname = "g", data = dt, control_group = "notyettreated",
             clustervars = "city_id", bstrap = TRUE, biters = 1000)

agg_simple  <- aggte(cs, type = "simple")
agg_dynamic <- aggte(cs, type = "dynamic")
agg_group   <- aggte(cs, type = "group")
```

The control group is not-yet-treated, doubly robust by default, with bootstrap standard errors; the standard errors come from a multiplier bootstrap, and fixing the seed before estimation is what makes it reproduce bit for bit. The simple aggregation gives 0.082 (SE 0.007), 95% confidence interval [0.068, 0.096], covering the truth 0.088. The group aggregation's overall is 0.078 (SE 0.007), and by cohort:

| cohort | CS group ATT (SE) | True value |
|---|---|---|
| g = 13 | 0.091 (0.010) | 0.106 |
| g = 19 | 0.086 (0.009) | 0.085 |
| g = 25 | 0.059 (0.012) | 0.058 |

g19 and g25 land almost exactly on the truth, while g13 is a bit low (its long-horizon cells can only be identified by the tiny never-treated group, so they are noisy). The difference between the simple and group aggregations (0.082 versus 0.078) is also telling: simple weights by the number of post-treatment cells, so g13 with the longest exposure gets more weight, and the two numbers answer different questions, both correct, but the paper must state which one is used.

### 6.7 BJS Imputation

```r
imp <- did_imputation(data = dt, yname = "lgmv", gname = "g",
                      tname = "month", idname = "city_id")

imp_es <- did_imputation(data = dt, yname = "lgmv", gname = "g",
                         tname = "month", idname = "city_id",
                         horizon = TRUE, pretrends = -12:-2)
```

The overall ATT is 0.078 (SE 0.005), the smallest standard error of the four estimators, and the efficiency advantage is cashed in. The point estimate is about two standard errors below the realized truth, and the $e = 0$ coefficient (0.005) is also low, which is the luck of this seed's noise draw, not a systematic bias; the lesson is that when writing a paper you report the estimate with its standard error, and do not claim that any modern estimator "exactly equals" the truth.

### 6.8 The Grand Reconciliation

| Estimator | Estimate | SE | Fraction of truth |
|---|---|---|---|
| True simple ATT (DGP) | 0.088 | | 100% |
| Static TWFE (feols) | 0.0227 | 0.0087 | 25.8% |
| CS simple ATT (att_gt, not-yet-treated) | 0.082 | 0.007 | 93.1% |
| SA average ATT (sunab) | 0.081 | 0.008 | 91.7% |
| BJS imputation (did_imputation) | 0.078 | 0.005 | 88.4% |

The reconciliation of the dynamic coefficients (the truth is 0 at all $e < 0$):

| $e$ | True | naive TWFE (SE) | sunab (SE) | CS dynamic (SE) | imputation (SE) |
|---|---|---|---|---|---|
| -6 | 0.000 | 0.001 (0.009) | 0.009 (0.008) | -0.002 (0.006) | 0.006 (0.007) |
| -3 | 0.000 | 0.009 (0.008) | 0.009 (0.009) | 0.001 (0.007) | 0.011 (0.008) |
| 0 | 0.019 | 0.015 (0.007) | 0.014 (0.007) | 0.012 (0.006) | 0.005 (0.006) |
| 3 | 0.022 | 0.020 (0.008) | 0.011 (0.010) | 0.015 (0.008) | 0.008 (0.007) |
| 6 | 0.041 | 0.036 (0.009) | 0.032 (0.012) | 0.034 (0.011) | 0.027 (0.008) |
| 9 | 0.093 | 0.089 (0.009) | 0.087 (0.011) | 0.089 (0.011) | 0.082 (0.008) |

At longer horizons, the truth is 0.141 at $e = 12$, 0.164 at $e = 18$, and 0.184 at $e = 23$; naive estimates 0.130 (0.010) at $e = 12$ and 0.162 (0.011) at $e = 18$; sunab and CS both estimate 0.118 at $e = 12$ (SE 0.014 and 0.015 respectively), and the long end is conservative with widening intervals because the usable cohorts become fewer, which is exactly the endpoint problem of Section 4.2.

The reconciliation of this section can be summarized as follows: on the same data, with all identification assumptions holding exactly, static TWFE delivers a quarter of the truth, while the three modern estimators deliver around ninety percent and give intervals that cover the truth; the difference all comes from the aggregation weights, and the Goodman-Bacon decomposition recovers the lost three quarters piece by piece.

## 7 Failure Modes of Identification and Robustness

In the simulation the identification assumptions are constructed, and in real research they can fail at any time. This section combs through the most common failure modes and their actionable responses.

The problem with the pre-trends test runs one layer deeper than "it is not proof." Roth (2022) gives two quantitative conclusions. In common event-study designs, the pre-trends test has low power against a trend difference large enough to flip the conclusion, and the test not rejecting often just means the data are insufficient rather than that the trends are parallel. More subtle is the conditional inference problem: continuing to infer only after the test passes, this screening itself distorts the distribution of the subsequent estimate, and when there is a small real violation, the sample that passes the screening is exactly the one whose violation direction favors you, so the bias is amplified instead. The actionable response is to replace the binary test with a quantitative statement: report the test's power against a violation of a given size, plot confidence intervals rather than just report stars, and do not delete parallel trends from the limitations discussion just because the pre test passed.

Rambachan and Roth (2023) push this idea all the way: rather than test "whether the violation is zero," ask "how large the violation must be to flip." Their framework no longer assumes parallel trends holds exactly, but constrains the possible post-period violation within a set whose size is calibrated by the trend difference observed in the pre period. The most commonly used version uses a parameter $M$ to limit the period-by-period change in the slope of the counterfactual trend difference: $M = 0$ corresponds to linearly extrapolating the pre-period trend, and larger $M$ allows more nonlinear bending. The output is a robust confidence set that widens as $M$ grows, and a breakdown value where the conclusion flips. When writing a paper, reporting this breakdown $M$ and arguing "how unreasonable a bend it would take to kill the result" is far more forceful than a one-liner "pre-trends are parallel." The R implementation is HonestDiD, and this chapter's simulation does not run it, and the argument's logic does not depend on code.

Adjacent to anticipation but different in mechanism is a class of threat called the Ashenfelter dip, which is essentially a failure of parallel trends rather than anticipation. In a training-program setting, people enter the program because their income temporarily dropped ($Y_{i1}(0)$ below some threshold), which is equivalent to a truncated sampling on the pre-period shock, so the treated group's untreated path contains a stretch of pure mean reversion that differs systematically from the control group's $Y(0)$ trend, and DiD counts this rebound as the training effect and overstates. Its dividing line from anticipation is: anticipation is a subject foreseeing future treatment and adjusting behavior in advance, a failure of no anticipation (A2); the dip is the selection rule picking units currently in a temporary trough, a failure of parallel trends (A1). Both contaminate the pre period, and their responses are also shared. The platform version is not at all far-fetched: if Kestrel's operations team deliberately cuts commissions in cities that just had a bad quarter of GMV in hope of a boost, the rebound gets credited to the commission. Diagnosis and response: plot the event study to see whether the last few periods with negative $e$ dip abnormally; move the base period from $e = -1$ to earlier (this is exactly what CS's anticipation parameter is for); exclude the dip period and re-estimate; use a longer pre window to see whether the trend returns to track.

Endogenous treatment timing is the norm rather than the exception in platform research. The platform cuts commissions first in the fastest-growing cities, which threatens both parallel trends and the interpretability of the weights. There are three lines of defense. Choose not-yet-treated rather than never-treated as the control group: the late-treated cities are also cities the platform will select, just scheduled later, and their comparability to the early-treated cities on the selection dimension is usually far better than that of cities never selected. Find quasi-exogenous variation in the treatment timing for an IV or a design argument: for example, if the scheduling of the commission cut is determined by bottlenecks unrelated to demand growth, such as engineering resources or compliance approval, then this institutional detail should be written out fully. Combine with synthetic control or synthetic DiD: weight together a counterfactual for each treated city, relaxing the implicit constraint that "all control cities get equal weight."

Spillover diagnosis must be specific to the mechanism. Sellers re-registering across cities, buyers ordering across cities, and platform-wide resource reallocation each have their own observable footprint. A practical check is to stratify the control cities by geographic or economic proximity to the treated cities, and if the estimated effect of "near controls" differs systematically from that of "far controls," spillover is basically confirmed, and you can then do a donut design (drop the border-band control cities), redefine the market boundary, or admit SUTVA has failed and build the spillover itself into an estimand. If what you suspect is a platform-wide equilibrium effect (the commission structure changed the whole platform's seller pricing equilibrium), the entire DiD framework does not apply, and a structural model or a platform-level experiment is needed.

Functional-form sensitivity was foreshadowed in Section 2: parallel trends under log and parallel trends under level are two assumptions. The operational discipline is to run the main result under both forms and put it in the appendix; if the two conclusions conflict, do not cherry-pick the pretty one, and admit in the main text that identification is fragile to functional form, and explain which form your mechanism supports more (a change in the commission rate acts proportionally on GMV, so log has a mechanistic basis).

The treatment of few clusters is in Section 4.7, and here we only stress one point most easily overlooked: the most severe problem is not that the total number of clusters is small, but that the treated clusters are few. If only 3 of 50 states passed some law, the effective sample size is close to 3 rather than 50, and CR1 gives you a pretty but false star.

The discipline of bad controls: do not control for any variable that will respond to treatment. Controlling for "the platform's marketing spend in that city" looks like it makes the regression cleaner, but marketing spend is itself a downstream result of the commission policy, and controlling for it amounts to blocking off a transmission path of the effect, and may also introduce collider bias. The same class of error is stuffing a lagged dependent variable into a DiD regression to "absorb the trend": a lagged term and fixed effects appearing together brings Nickell bias, and the estimand also changes flavor, and it is not a diagnostic tool for parallel trends. To deal with trend differences, the right road is the conditioning of Section 4.6 or the sensitivity analysis of Rambachan-Roth.

Stringing these failure modes together, the credibility of DiD ultimately lies not in the ingenuity of the regression specification but in the design itself: why treatment happens to these units at these times, and whether this mechanism is unrelated to the counterfactual trend of the outcome. The pre-trends plot, the Rambachan-Roth sensitivity analysis, not-yet-treated controls, and functional-form robustness checks are all evidence provided around this one question, and none of them can substitute for a substantive judgment about the treatment assignment mechanism.

## 8 Further Reading

::: {.readings}
Required reading, in suggested reading order:

- Roth, Sant'Anna, Bilinski and Poe (2023, Journal of Econometrics). The overview to read first, organizing the 2018 to 2023 DiD literature into a single map, and reading it first makes tackling the original papers far more efficient.
- Goodman-Bacon (2021, Journal of Econometrics). The decomposition theorem for static TWFE, the full version of Section 4.3 of this chapter, with the focus on the comparative statics of the weights.
- Callaway and Sant'Anna (2021, Journal of Econometrics). The original paper on the group-time ATT framework, with the focus on the difference between the two control-group assumptions and the aggregation part.
- Bertrand, Duflo and Mullainathan (2004, QJE). The origin of the inference problem, and the research design of the placebo experiment is itself worth borrowing.
- Baker, Larcker and Wang (2022, Journal of Financial Economics). A restaging of the failure of staggered TWFE using financial data, with the most complete practical checklist, suitable for a self-check before writing a paper.

Further reading:

- Sun and Abraham (2021, Journal of Econometrics). The rigorous statement of the contamination theorem and the IW estimator, read Propositions 1 to 3.
- Borusyak, Jaravel and Spiess (2024, Review of Economic Studies). The imputation framework and the efficiency result, note its stronger requirement on parallel trends.
- de Chaisemartin and D'Haultfœuille (2020, American Economic Review). The solution when treatment is reversible.
- Roth (2022, AER: Insights). The power of the pre-trends test and pre-test distortion.
- Rambachan and Roth (2023, Review of Economic Studies). The honest inference framework, the original paper on the $M$-parameter sensitivity analysis.
- Roth and Sant'Anna (2023, Econometrica). When parallel trends is invariant to functional form.
- Sant'Anna and Zhao (2020, Journal of Econometrics). The doubly robust DiD estimator, the theoretical basis of the did package's default method.
:::

::: {.apa-refs}
- Baker, A. C., Larcker, D. F., & Wang, C. C. Y. (2022). How much should we trust staggered difference-in-differences estimates? *Journal of Financial Economics, 144*(2), 370-395. https://doi.org/10.1016/j.jfineco.2022.01.004
- Bell, R. M., & McCaffrey, D. F. (2002). Bias reduction in standard errors for linear regression with multi-stage samples. *Survey Methodology, 28*(2), 169-181.
- Bertrand, M., Duflo, E., & Mullainathan, S. (2004). How much should we trust differences-in-differences estimates? *The Quarterly Journal of Economics, 119*(1), 249-275. https://doi.org/10.1162/003355304772839588
- Borusyak, K., Jaravel, X., & Spiess, J. (2024). Revisiting event-study designs: Robust and efficient estimation. *The Review of Economic Studies, 91*(6), 3253-3285. https://doi.org/10.1093/restud/rdae007
- Callaway, B., & Sant'Anna, P. H. C. (2021). Difference-in-differences with multiple time periods. *Journal of Econometrics, 225*(2), 200-230. https://doi.org/10.1016/j.jeconom.2020.12.001
- Cameron, A. C., Gelbach, J. B., & Miller, D. L. (2008). Bootstrap-based improvements for inference with clustered errors. *The Review of Economics and Statistics, 90*(3), 414-427. https://doi.org/10.1162/rest.90.3.414
- Card, D., & Krueger, A. B. (1994). Minimum wages and employment: A case study of the fast-food industry in New Jersey and Pennsylvania. *American Economic Review, 84*(4), 772-793.
- Card, D., & Krueger, A. B. (2000). Minimum wages and employment: A case study of the fast-food industry in New Jersey and Pennsylvania: Reply. *American Economic Review, 90*(5), 1397-1420. https://doi.org/10.1257/aer.90.5.1397
- Conley, T. G., & Taber, C. R. (2011). Inference with "difference in differences" with a small number of policy changes. *The Review of Economics and Statistics, 93*(1), 113-125. https://doi.org/10.1162/rest_a_00049
- de Chaisemartin, C., & D'Haultfœuille, X. (2020). Two-way fixed effects estimators with heterogeneous treatment effects. *American Economic Review, 110*(9), 2964-2996. https://doi.org/10.1257/aer.20181169
- Gardner, J. (2022). *Two-stage differences in differences* (arXiv:2207.05943). arXiv. https://doi.org/10.48550/arXiv.2207.05943
- Goodman-Bacon, A. (2021). Difference-in-differences with variation in treatment timing. *Journal of Econometrics, 225*(2), 254-277. https://doi.org/10.1016/j.jeconom.2021.03.014
- Greenwood, B. N., & Wattal, S. (2017). Show me the way to go home: An empirical investigation of ride-sharing and alcohol related motor vehicle fatalities. *MIS Quarterly, 41*(1), 163-187. https://doi.org/10.25300/MISQ/2017/41.1.08
- Liu, L., Wang, Y., & Xu, Y. (2024). A practical guide to counterfactual estimators for causal inference with time-series cross-sectional data. *American Journal of Political Science, 68*(1), 160-176. https://doi.org/10.1111/ajps.12723
- MacKinnon, J. G., & Webb, M. D. (2020). Randomization inference for difference-in-differences with few treated clusters. *Journal of Econometrics, 218*(2), 435-450. https://doi.org/10.1016/j.jeconom.2020.04.024
- Neumark, D., & Wascher, W. (2000). Minimum wages and employment: A case study of the fast-food industry in New Jersey and Pennsylvania: Comment. *American Economic Review, 90*(5), 1362-1396. https://doi.org/10.1257/aer.90.5.1362
- Rambachan, A., & Roth, J. (2023). A more credible approach to parallel trends. *The Review of Economic Studies, 90*(5), 2555-2591. https://doi.org/10.1093/restud/rdad018
- Roth, J. (2022). Pretest with caution: Event-study estimates after testing for parallel trends. *American Economic Review: Insights, 4*(3), 305-322. https://doi.org/10.1257/aeri.20210236
- Roth, J., & Sant'Anna, P. H. C. (2023). When is parallel trends sensitive to functional form? *Econometrica, 91*(2), 737-747. https://doi.org/10.3982/ECTA19402
- Roth, J., Sant'Anna, P. H. C., Bilinski, A., & Poe, J. (2023). What's trending in difference-in-differences? A synthesis of the recent econometrics literature. *Journal of Econometrics, 235*(2), 2218-2244. https://doi.org/10.1016/j.jeconom.2023.03.008
- Sant'Anna, P. H. C., & Zhao, J. (2020). Doubly robust difference-in-differences estimators. *Journal of Econometrics, 219*(1), 101-122. https://doi.org/10.1016/j.jeconom.2020.06.003
- Sun, L., & Abraham, S. (2021). Estimating dynamic treatment effects in event studies with heterogeneous treatment effects. *Journal of Econometrics, 225*(2), 175-199. https://doi.org/10.1016/j.jeconom.2020.09.006
:::
