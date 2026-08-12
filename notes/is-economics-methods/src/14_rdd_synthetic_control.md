---
title: "RDD & Synthetic Control"
subtitle: "From the Cutoff to the Synthetic Counterfactual"
seriesline: "Foundations of Information Systems Economics · Chapter 14"
series: "Foundations of Information Systems Economics: Empirical Methods for Information Systems Economics"
sidebartitle: "Foundations of Information Systems Economics"
sidebarsub: "Theory · Methods · Applications"
chapterlabel: "Chapter 14 · RDD & Synthetic Control"
---

## Introduction

The Aster marketplace hands out a "Preferred" badge to sellers whose rating just clears a line. It is no surprise that badged sellers sell more online: highly rated sellers were probably already better to begin with. What is genuinely useful is a narrower picture: compare the sellers who fell just short of the badge, on the left of the threshold, with those who just cleared it, on the right. If nothing other than the badge jumps discontinuously near the threshold, the break in sales carries a causal meaning.

The same platform is separately piloting a fee reform in a single city. Here there is no threshold, and no natural "identical twin city" to serve as a control. Synthetic control takes another route: it lets several untreated cities, combined by weights, jointly play the role of the counterfactual city that does not exist. Whether the pre-treatment path can be stitched together closely enough decides whether the gap that opens up after treatment is worth believing.

RDD and synthetic control seem to share no common estimating equation, yet they share a research habit: identifying power comes from the concrete texture of institutions. The former leans on continuity near a threshold and answers only a local effect; the latter leans on pre-treatment fit and a donor pool, and answers only the counterfactual for the treated unit. Badges, tiers, and eligibility lines on a platform do not automatically make a discontinuity design, and a single-city pilot does not automatically come with a synthetic control. This chapter follows Aster's two storylines and presses on three questions: whether the threshold can be manipulated, what else jumps at the cutoff, and when a weighted set of cities can truly stand in for "Aster without the reform."

## 1 Two Numbers That Look Wrong

::: {.case}
Aster is a fictional two-sided online marketplace. We use simulated data to tell two design stories, and the payoff is that the data-generating process is fully known, so every estimator can be reconciled against the truth.

The first storyline is the Preferred badge. Aster gives a "Preferred" badge, displayed on the storefront page, to sellers whose quality score reaches 80. Management wants to know how large the causal effect of this badge is on a seller's monthly sales (log GMV). The most natural calculation is to subtract the average log GMV of unbadged sellers from that of badged ones, and this difference is **1.891**. Put into a growth deck, it says the badge lifts sales by nearly 1.9 log points, a number juicy enough for a team to claim credit.

The second storyline is the single-city fee reform. Aster rolls out a seller fee reform in one metro only, with the reform date recorded as $T_0$, while twenty other metros carry on as usual. Management wants the effect of this city-level policy on that metro's log GMV. The most natural calculation is a two-period, two-group DiD: the change in the treated city before and after the reform, minus the average change over the twenty donor cities over the same window, a number equal to **0.244**.
:::

Neither number is credible, and the direction of the trouble is something earlier chapters have already taught. The badge's 1.891 mixes in the selection bias of Chapter 9: high-scoring sellers were already selling more, and even with no badge at all, stores above 80 have far higher log GMV than stores below 80. Because the badge happens to go only to high-scoring sellers, the "badged vs. unbadged" comparison tangles the effect of the badge together with the difference that "high-scoring sellers were simply stronger." Because this is a simulation, we know the true causal jump of the badge at the threshold is only **0.200**, so that 1.891 is about 9.5 times the truth, of which nearly nine-tenths is selection bias.

The fee reform's 0.244 stumbles on the parallel trends of Chapter 13. In this DGP the cities differ in their growth trends to begin with: the treated city's pre-reform monthly log-GMV trend slope is **0.131**, while the donor cities average only **0.117**. The simple average of twenty cities is not the path the treated city would have taken had it not reformed, and using it as the counterfactual folds the pre-existing trend difference into the policy effect. The true post-reform average effect is only **0.095**, so that 0.244 is about 2.6 times the truth.

Each of the two designs offers its own way around the trap. RDD says: do not compare all high-scoring sellers against all low-scoring ones, look only at a very narrow band around the 80-point line, where sellers at 79 and 81 barely differ in ability, category, or operations, and who happens to clear the line versus who just misses it is nearly random. The small log-GMV difference between those just above and just below the line is the badge's clean effect.

![Aster's Preferred badge discontinuity: the horizontal axis is quality score, the vertical axis is log GMV. Sellers below 80 are binned by score and averaged (blue dots), and those at or above 80 (orange dots); each side gets a local-linear fit within the CCT-optimal bandwidth, and the dashed line marks the cutoff. The vertical jump at the threshold is the badge effect, labeled RDD jump = 0.200 (truth 0.200). Note how log GMV rises steeply with score overall, which is exactly why the naive comparison is absurdly high, while the jump at the discontinuity is far smaller than this steep slope.](assets/fig/fig_14_rdd.svg)

SC says: do not use the simple average of twenty cities as the counterfactual, but pick a convex combination from the donors that reproduces the treated city's log-GMV trajectory before the reform as closely as possible; the closer the pre-reform fit, the more the gap that opens after the reform between the treated city's actual path and this synthetic version looks like the net effect of the policy.

These two ideas are the spine of the chapter. On what grounds they hold, and where each of them fails, are the questions the next seven sections answer.

## 2 The Economic Model and the Estimand

Before estimating anything, let us be clear about two things: exactly which causal quantity each design is after, and what notation writes it down. Both storylines are built on the potential outcomes of Chapter 9, but the target estimands have very different shapes.

Start with RDD. Write the seller's quality score as the running variable (also called the forcing variable) $X_i$, with cutoff $c$, here $c = 80$. The treatment indicator $D_i$ records whether the seller got the badge. The defining feature of RDD is that the treatment probability jumps at $c$:

$$\lim_{x\downarrow c}\Pr(D_i = 1 \mid X_i = x) \ \neq\ \lim_{x\uparrow c}\Pr(D_i = 1 \mid X_i = x).$$

A **sharp RD** is the special case where this jump goes exactly from 0 to 1, that is, clearing the line guarantees a badge and missing it guarantees none, $D_i = \mathbf{1}[X_i \geq c]$. Each seller has two potential outcomes $Y_i(1), Y_i(0)$, with observed outcome $Y_i = D_i Y_i(1) + (1 - D_i) Y_i(0)$. The structural picture behind RDD is this: the untreated outcome is a smooth function of the running variable, and the only jump at $c$ comes from treatment,

$$Y_i = Y_i(0) + \tau_i \cdot \mathbf{1}[X_i \geq c], \qquad Y_i(0) = f(X_i) + \varepsilon_i, \ \ f \text{ smooth}.$$

The RDD estimand is not the average treatment effect over everyone, but the local average effect at that one point, the cutoff:

$$\tau_{\mathrm{RD}} = \mathbb{E}\big[\,Y_i(1) - Y_i(0) \mid X_i = c\,\big].$$

This has to be made plain from the outset: without extra assumptions, RDD can only learn the effect at the single point $X = c$. How much a seller on the 80-point line gains from getting the badge, it can answer; what would happen if you gave the badge to a 60-point seller, it cannot. This is RDD's built-in limitation, and the price it pays for strong internal validity, discussed formally in Section 3.3.

Now SC. Here there is a single treated unit (write it $i = 1$, the reform city Aster-Metro), treated for $t > T_0$; another $j = 2, \dots, J+1$ donor cities are never treated. The treatment indicator $D_{it} = 1$ if and only if $i = 1$ and $t > T_0$ (with $T_0$ the last untreated period, here month 24). In potential outcomes,

$$Y_{it} = Y_{it}(0) + \tau_{it} D_{it}.$$

The SC estimand is the ATT path of the treated unit, an effect curve that varies over time:

$$\tau_{1t} = Y_{1t}(1) - Y_{1t}(0) = \underbrace{Y_{1t}}_{\text{observed}} - \underbrace{Y_{1t}(0)}_{\text{counterfactual}}, \qquad t > T_0.$$

$Y_{1t}$ is observed, and all of the difficulty lands on $Y_{1t}(0)$, the log GMV the reform city would have realized had it not reformed. The whole machinery of SC does one thing: it constructs this unobserved counterfactual path.

To sum up: the RDD target is the local effect at the cutoff, $\tau_{\mathrm{RD}}$, a single number pinned to $X = c$; the SC target is the ATT path of a single treated unit, $\tau_{1t}$, a curve that opens over time. Both are local, one local at the cutoff, the other local at the treated unit. The next section makes clear the assumptions under which each of these quantities is identified.

## 3 Identification

Identification logically precedes estimation. Identification asks whether, under a set of assumptions about the counterfactual world, the desired causal quantity can be written as a function of the observable data distribution; this question depends only on the design and the institutional background, and has nothing to do with which estimator or how large a sample you use afterward. This is the most fine-grained analytical section in the chapter, and it deliberately builds from simple to hard: first it makes the continuity assumption fully clear in sharp RDD, then discusses the one testable shadow it leaves and its local limitation, then connects fuzzy RD back to the IV of Chapter 12, and finally turns to SC, running from the intuition of a factor model all the way through the core algebra of "matching the pre-period is matching the unobserved loadings," and laying out where SC parts ways from DiD and the conditions under which it is credible.

### 3.1 A Break on a Line: The Continuity Assumption of Sharp RDD

Start with an example anyone understands. Suppose the government rules that households with annual income below $25{,}000 can collect a subsidy. To estimate the subsidy's effect on household consumption, directly comparing households that collect it against those that do not will not do, because poor households simply consume less to begin with. But turn it around: a household earning $24{,}999 and one earning $25{,}001 barely differ in anything except that one just qualifies for the subsidy and the other just misses. Within a very narrow neighborhood of this threshold, who lands above the line and who below is nearly a matter of luck. So the consumption difference between those thin layers of households on either side of the threshold can be attributed cleanly to the subsidy. That is the whole intuition of RDD, in one sentence: near the cutoff, treatment is as good as randomly assigned.

Aster's badge is the same idea. Sellers at 79 and 81 are almost identical in category, inventory, and operations, and the only difference is that one just cleared the 80-point line and got the badge while the other just missed. Turning this intuition into a formal identifying assumption rests on continuity.

::: {.assumption}
**Continuity** The conditional expectation functions of both potential outcomes are continuous at the cutoff:

$$x \mapsto \mathbb{E}[Y_i(0) \mid X_i = x] \ \text{ and } \ x \mapsto \mathbb{E}[Y_i(1) \mid X_i = x] \ \text{ are continuous at } x = c.$$
:::

What this sentence constrains is worth reading word by word. It says that if you took the badge away, the seller's log GMV as a function of score should pass smoothly through the 80-point line, and should not sprout a step out of nowhere at exactly 80. In other words, nothing other than the badge jumps at this one point, 80. It allows $\mathbb{E}[Y(0) \mid X]$ to climb steeply with score, allows it to bend and be nonlinear, allows high-scoring and low-scoring sellers to be worlds apart; the one thing it forbids is a break unrelated to the badge at the cutoff itself. Precisely because it is continuous, the sellers who just miss the line on the left are a legitimate counterfactual for those who just clear it on the right: their $\mathbb{E}[Y(0)]$ joins up seamlessly at the threshold.

With continuity, identification is one step of taking a limit. The jump in the observed conditional mean at the cutoff equals exactly the local treatment effect at the cutoff:

$$\tau_{\mathrm{RD}} = \lim_{x\downarrow c}\mathbb{E}[Y_i \mid X_i = x] - \lim_{x\uparrow c}\mathbb{E}[Y_i \mid X_i = x] = \mathbb{E}[\tau_i \mid X_i = c].$$

The derivation is only three lines. Approaching from the right of the cutoff ($x \downarrow c$), these sellers all have the badge, $D_i = 1$, so $\mathbb{E}[Y_i \mid X_i = x] \to \mathbb{E}[Y_i(1) \mid X_i = c]$, using the continuity of $\mathbb{E}[Y(1) \mid X]$. Approaching from the left ($x \uparrow c$), these sellers all lack the badge, $D_i = 0$, so $\mathbb{E}[Y_i \mid X_i = x] \to \mathbb{E}[Y_i(0) \mid X_i = c]$, using the continuity of $\mathbb{E}[Y(0) \mid X]$. Subtracting the two, the pile of confounding terms on the right cancels because it joins up continuously, leaving only $\mathbb{E}[Y_i(1) - Y_i(0) \mid X_i = c]$, exactly the local treatment effect at the cutoff. This limit difference is all that RDD sets out to estimate.

::: {.intuition}
The continuity assumption has a reassuring property: it shifts the burden of identification from "I controlled for all the confounders correctly" to "nothing else jumps at this one point." You do not need to observe the seller's ability, category, or operational quality, nor do they need to be balanced across the two sides of the threshold. You only need to believe that, apart from the badge, no second thing switches on at exactly 80. This is far easier to check than the "I controlled for everything I should have" faith of selection on observables, because whether something else jumps at the threshold can often be seen from a plot and argued from the institutions. The credibility of RDD is, to a large extent, the credibility of "this cutoff triggers only the one thing, the badge."
:::

It bears emphasizing that continuity itself is not testable. It constrains the counterfactual $\mathbb{E}[Y(0) \mid X = c]$, and at that point we observe only the badged side, never what those units would look like without the badge. So the right attitude is not to "test continuity" but to plot and to argue. In a word, this is not a testable assumption, go plot it.

Continuity is not the only framework that makes RDD clear. A parallel line of thought is called local randomization: instead of talking about whether the conditional expectation is continuous at the cutoff, it directly assumes that within a very narrow window $W$ around the cutoff, who lands above the line and who below is nearly randomly assigned, as if that thin layer of units on either side of the threshold took part in a tiny randomized experiment. This assumption is stronger than continuity and more local, in that it concerns only within the window $W$ and does not constrain the limiting behavior of the whole curve. What it buys is more direct estimation and inference: within the window, the mean difference between treated and control units is the effect estimate, and inference can go the Fisher route of randomization inference, repeatedly permuting the treatment labels within the window to compute a p-value, yielding finite-sample exact conclusions without leaning on large-sample limits. When should you reach for it? When there are very few observations near the cutoff, or when the running variable is discrete with obvious integer heaping, the smooth approximation that the continuity framework relies on is shaky, and local randomization paired with randomization inference is often more robust. The main line of this chapter stays with the continuity framework, but keep in mind that this second route is available.

### 3.2 No-Manipulation: The One Assumption That Leaves a Shadow

Continuity has one concrete behavioral premise worth singling out, because it is one of the few things in RDD that can be indirectly tested. This premise is no-manipulation, or no-sorting: sellers cannot precisely control which side of the threshold they land on.

Why manipulation kills identification is clearest back in the subsidy example. If a household can precisely under-report its income to slip below $25{,}000 and collect the subsidy, then just below the threshold you get a crowd of households whose income was actually higher but who deliberately dialed down their filing, and they are no longer comparable in other respects to the un-manipulating households just above the threshold. Worse, this active sorting makes the density of the running variable break at the threshold, bunching just below it, and the "nearly random comparability across the two sides" that continuity relies on collapses along with it. The corresponding threat in Aster's setting is: if sellers can precisely push their quality score up to 80 by tactics like fake orders or category changes, then just above the 80-point line you get a batch of sellers "on their tiptoes to reach the line," no longer comparable to those just below.

The observable footprint that manipulation leaves is a jump in the density of the running variable at the cutoff, and that gives us a test.

::: {.assumption}
**McCrary density test** Estimate the density of the running variable separately on each side of the cutoff, and test whether they are continuous at $c$:

$$\theta = \ln\Big(\lim_{x\downarrow c}\hat f_X(x)\Big) - \ln\Big(\lim_{x\uparrow c}\hat f_X(x)\Big), \qquad H_0:\ \theta = 0.$$
:::

A sudden lift in the density above the threshold (bunching) is evidence of sorting, showing that someone succeeded in moving from below the line to above it. McCrary (2008) gave the original implementation of this test, and the modern default is Cattaneo, Jansson and Ma's `rddensity`, which estimates the density directly with a local polynomial, avoiding pre-binning, and delivers a bias-corrected robust test. The point is to place this test correctly: it is neither a necessary nor a sufficient condition for continuity, but essentially a testable implication of the behavioral no-manipulation assumption, that is, a falsification test. Density continuity cannot prove the potential outcomes are continuous (not sufficient); and the density may break for reasons unrelated to the outcome, such as integer heaping or administrative rules, so a density break by itself does not logically amount to identification failure (not necessary). But density bunching above the threshold is strong evidence of sorting, and the moment you see it you should be highly alert that continuity has been broken.

Besides density, the RDD diagnostic checklist has a companion check standing alongside it: predetermined covariates $W$ (features fixed before the score was obtained, such as the seller's tenure or historical category) should not jump at the cutoff. The logic is the same as above: if even these features, unrelated to the badge and fixed before treatment, show a step at 80, that says the sellers on either side of the threshold were not comparable to begin with, and continuity is suspect. Covariates smooth at the cutoff are placebo-style evidence in favor of identification.

Put these two together and you get one of RDD's standout advantages over other quasi-experimental designs: although its core assumption cannot be tested directly, it leaves two rather powerful observable shadows, no density jump and no jump in predetermined covariates. The exclusion restriction of Chapter 12 leaves no such shadow at all and can only be argued by storytelling. On this point RDD is much luckier.

### 3.3 RDD Is a Local Parameter

RDD's strong internal validity comes with a cost there is no getting around: what it identifies is extremely local. It learns the effect only at the cutoff, on that thin layer of units near the cutoff. That a seller barely qualifying on the 80-point line gains $\tau_{\mathrm{RD}}$ from the badge tells you nothing directly about what would happen if you gave the badge to a 90-point seller or a 60-point seller, because those sellers are too far from the threshold for RDD's "nearly random" to cover them at all.

To extrapolate $\tau_{\mathrm{RD}}$ to units far from the cutoff, you must additionally assume that the treatment effect is homogeneous across the whole score range, or assume a specific functional form that lets the effect vary smoothly with score, and RDD itself neither provides nor tests these assumptions. So the honest practice is to confine the wording of your conclusion to the neighborhood of the threshold, making clear this is "the effect for sellers on the 80-point line," not "the effect of the badge for all sellers." This is of a piece with the extrapolation anxiety of LATE in Chapter 12: there you learned the effect for compliers, here you learn the effect for units at the threshold, and in both the design illuminates only a small group. The strength of identification is bought precisely with the narrowness of extrapolation.

### 3.4 Fuzzy RD: IV at the Cutoff

So far we have assumed that clearing the line guarantees the badge. Real platform rules are often not so clean. Aster's actual situation is: among sellers who reach 80, some never activate the badge for one reason or another, while among sellers below 80 a few obtain it through appeals, special operational approval, and the like. So the treatment probability at the threshold does not jump from 0 to 1, but from a lower level to a higher level. This is fuzzy RD.

Since the treatment probability does not go from 0 to 1, the jump in the outcome at the threshold reflects only that fraction of sellers whose badge status was nudged to change by the threshold. To recover "how much log GMV each change in badge status buys," you divide the jump in the outcome by the jump in treatment, which is exactly the Wald ratio:

$$\tau_{\mathrm{FRD}}(c) = \frac{\displaystyle\lim_{x\downarrow c}\mathbb{E}[Y_i \mid X_i = x] - \lim_{x\uparrow c}\mathbb{E}[Y_i \mid X_i = x]}{\displaystyle\lim_{x\downarrow c}\Pr(D_i = 1 \mid X_i = x) - \lim_{x\uparrow c}\Pr(D_i = 1 \mid X_i = x)} = \frac{\text{jump in outcome}}{\text{jump in treatment}}.$$

This shape should look familiar. When Chapter 12 covered IV, the Wald estimator was exactly the reduced form divided by the first stage. Fuzzy RD just moves that IV machinery onto the cutoff: it takes the indicator for crossing the threshold, $Z_i = \mathbf{1}[X_i \geq c]$, as the instrument, uses it to instrument the actual badge status $D_i$, with the numerator the reduced-form effect of crossing the threshold on log GMV and the denominator the first-stage effect of crossing on badge status. The assumptions required for identification correspond to those in Chapter 12: crossing the threshold must genuinely change the badge probability (relevance, the denominator nonzero), crossing can affect log GMV only through the badge (a cutoff version of exclusion), and monotonicity at the threshold (crossing the line can only increase, never decrease, the probability of getting the badge, that is, no defiers).

In sharp RD, everyone at the threshold is a complier nudged by the threshold, so what you estimate is the local effect at the threshold. In fuzzy RD, near the threshold there are also always-takers (who have the badge even without clearing the line) and never-takers (who do not activate it even after clearing), and under monotonicity, $\tau_{\mathrm{FRD}}$ identifies the LATE for compliers at the cutoff, the effect on those sellers "who got the badge only because they crossed the 80-point line, and would not have otherwise." This is a doubly local quantity: local at the cutoff, and local among the compliers. The line Chapter 12 kept returning to holds again here: IV learns only about compliers, and fuzzy RD is no exception.

::: {.warning}
Fuzzy RD estimation has one detail you must hold onto, a close cousin of the forbidden regression of Chapter 12. When doing 2SLS you should include a slope for the running variable on each side of the cutoff (that is, use $X_i - c$ rather than $X_i$, and let the slopes differ across the two sides), and you must enter $X_i - c$ rather than $X_i$ into the regression, that is, recenter the running variable at the cutoff. The reason is that if you do not recenter and use $X_i$ directly, the coefficient on the instrument $Z_i$ is no longer the pure effect of "crossing the threshold" on the badge, it gets tangled up with the intercept and slope of the polynomial at $X_i = 0$, and what you estimate mixes the jump at the threshold with the level of the polynomial at the origin, rather than the jump at the threshold. Only after recentering can the coefficient on $Z_i$ be read cleanly as the jump in treatment probability at that one point, the cutoff.
:::

### 3.5 The Identification Logic of Synthetic Control: Matching Unobserved Factor Loadings with a Convex Combination

Now turn to SC. Its identification logic is entirely different from RDD's, but no less elegant. Build the intuition first. When Chapter 13 covered DiD, the counterfactual was borrowed from the change in the control group: however the treated group would have changed, you assumed it changes by the same amount as the control group. SC does this "borrowing" more finely: rather than using the simple average of twenty donor cities as the counterfactual, pick a weighted combination from those twenty cities so that this "synthetic city" reproduces, as closely as possible, each feature and each period's outcome of the treated city before the reform. How closely it fits before the reform is visible and measurable; once it fits, you have reason to believe it keeps fitting after the reform, so the gap between the treated city's actual path and this synthetic version is the effect.

To turn "if it fits, it is credible" into an identification result, you need a model of how the untreated outcome is generated. What SC assumes is interactive fixed effects, also called a factor model:

$$Y_{it}(0) = \delta_t + \boldsymbol\theta_t \mathbf{Z}_i + \boldsymbol\lambda_t \mathbf{F}_i + u_{it}.$$

Here $\delta_t$ is a common time effect, $\mathbf{Z}_i$ are observable covariates (whose coefficients $\boldsymbol\theta_t$ may vary over time; here $\boldsymbol\theta_t$ follows the Abadie-Diamond-Hainmueller notation for time-varying covariate coefficients, unrelated to the $\theta$ used elsewhere in the book as a structural parameter; likewise $\mathbf{Z}_i$ here refers specifically to observable covariates, differing in meaning from the $Z_i = \mathbf{1}[X_i \geq c]$ used as an instrument in Section 3.4 and throughout the book, which merely reuses the same letter), $\boldsymbol\lambda_t$ is a set of time-varying common factors, $\mathbf{F}_i$ is each city's loadings on those factors, and $u_{it}$ is idiosyncratic noise. This model is far more general than the additive two-way fixed effects $\gamma_i + \gamma_t$ of Chapter 13. Additive fixed effects are a special case of it: for those time-varying common factors, the cities share the same loadings, so the trends they induce are parallel to one another, and cross-city heterogeneity lives only in a time-invariant intercept (equivalently, a constant unchanging factor paired with city-specific loadings plays the role of $\gamma_i$). The genuinely general case is that loadings may differ, meaning different cities respond by different amounts to the same time-varying common shock, so their trends are naturally not parallel. Aster's DGP is built exactly this way, which is also the root of why the two cities' trend slopes differ in Section 1 and DiD fails.

The key algebra is a single step. Subtract the synthetic city (a convex combination of donors, with weights $w_j$) from the treated city, in the untreated state:

$$Y_{1t}(0) - \sum_j w_j Y_{jt}(0) = \boldsymbol\theta_t\underbrace{\Big(\mathbf{Z}_1 - \sum_j w_j \mathbf{Z}_j\Big)}_{\text{match directly}} + \boldsymbol\lambda_t\underbrace{\Big(\mathbf{F}_1 - \sum_j w_j \mathbf{F}_j\Big)}_{\text{match indirectly}} + \sum_j w_j (u_{1t} - u_{jt}).$$

The $\mathbf{Z}$ in the first term is observable, and you can require the weights to match it directly. The real trouble is the $\mathbf{F}_i$ in the second term: the factor loadings are unobservable and cannot be matched directly. The core insight of SC (Abadie, Diamond and Hainmueller 2010) is: if a set of weights can simultaneously reproduce the treated city's entire pre-reform outcome path $(Y_{11}, \dots, Y_{1T_0})$ and its covariates, then it approximately also matches those unseen factor loadings, $\mathbf{F}_1 \approx \sum_j w_j \mathbf{F}_j$.

::: {.theorem}
(Loading matching, an approximate result) Under the factor model, if the weights $\mathbf{w}$ make the synthetic city reproduce all of the treated city's pre-period outcomes and covariates, that is, for all $t \leq T_0$ we have $Y_{1t} = \sum_j w_j Y_{jt}$ and $\mathbf{Z}_1 = \sum_j w_j \mathbf{Z}_j$, then when the pre-period is long enough, the idiosyncratic noise $u_{it}$ is small enough relative to the factor signal, and the pre-period fit is indeed good, we have approximately $\mathbf{F}_1 \approx \sum_j w_j \mathbf{F}_j$, with the matching error shrinking as the number of pre-periods grows and the noise-to-signal ratio falls. So for $t > T_0$, both the observable term and the loading term on the right-hand side above nearly cancel, and the gap

$$\hat\tau_{1t} = Y_{1t} - \sum_{j} \hat w_j Y_{jt}$$

estimates the treatment effect $\tau_{1t}$ approximately without bias, with the bias tending to zero as the pre-period lengthens. This is precisely why Section 3.6 lists "a long enough pre-period" as a credibility condition.
:::

The intuition is this: the pre-reform outcome path is itself an information-rich fingerprint of the factor loadings. How exposed a city is to each common factor is imprinted entirely in the trajectory of how its log GMV grew and fluctuated seasonally over the years. If two cities track each other month by month over a long enough pre-period, their exposures to the common factors are nearly identical, even though we have never directly observed those exposures. Matching the visible outcome path indirectly matches the invisible loadings, and this is the source of all of SC's persuasive power, which also explains why the quality of the pre-period fit is SC's lifeline.

The weights are not picked at random; they are pinned down by a constrained minimum-distance problem. Stack the covariates and pre-period outcomes into a predictor vector, written $\mathbf{X}_i = (\mathbf{Z}_i, Y_{i1}, \dots, Y_{iT_0})$, where the boldface $\mathbf{X}_i$ refers specifically to SC's stacked predictors, a different thing from the scalar running variable $X_i$ in the RDD half, merely reusing the same letter,

$$\hat{\mathbf{w}} = \arg\min_{\mathbf{w}}\ \Big(\mathbf{X}_1 - \textstyle\sum_j w_j \mathbf{X}_j\Big)' \boldsymbol\Omega \Big(\mathbf{X}_1 - \textstyle\sum_j w_j \mathbf{X}_j\Big) \quad\text{s.t. } \sum_{j} w_j = 1,\ w_j \geq 0,$$

where $\boldsymbol\Omega$ is a positive-definite predictor weighting matrix. The two constraints are SC's signature: weights nonnegative and summing to one, together called the simplex constraints. They guarantee that the synthetic city is a convex combination of donors, lying inside the donors' convex hull, with no extrapolation. This is a deliberate choice of SC relative to ordinary regression: the implicit weights of ordinary regression can be negative or larger than one, which amounts to allowing extrapolation into regions with no donor support, whereas SC uses nonnegative-summing-to-one to lock itself within the support of the data, at the cost that sometimes the convex hull is too small to reconstruct the treated city.

### 3.6 Where SC Parts Ways from DiD, and When SC Is Credible

SC and the DiD of Chapter 13 deserve a head-on comparison, because they answer the same kind of question while differing in the strictness of their underlying identifying assumptions. DiD's two-way fixed effects model $y_{it} = \delta D_{it} + \gamma_i + \gamma_t + u_{it}$ assumes the fixed effects are additively separable: the baseline levels across cities can differ ($\gamma_i$), but the time evolution treats all cities alike ($\gamma_t$). The moment different cities respond differently to a common shock, that is, once a loading-varying term like $\boldsymbol\lambda_t \mathbf{F}_i$ is present, parallel trends breaks. SC relaxes exactly this layer: it allows cities to have different loadings on the common factors, so trends can be non-parallel, and then it matches those loadings away by reproducing the pre-period path. In a word, DiD tolerates only time-invariant confounding, while SC can tolerate time-varying confounding that enters through factor loadings.

There is also a very practical complaint about DiD: on what grounds do you pick its control group? In DiD, which cities to use as controls is often decided by the researcher's gut. SC turns this choice into data, replacing a hand-picked control group with a convex combination optimized from pre-period fit, and the weights are computed, transparent, and auditable. This is a major methodological advance of SC over DiD.

The price is that SC has its own full set of credibility conditions, which a paper should put on the table one by one:

1. **The pre-period fit must be good.** The synthetic city must hug the treated city's pre-period path (a small pre-period MSPE). If the fit is poor, the loading-matching argument does not hold, and the later gap cannot be interpreted. This is the first-order check.
2. **The treated unit must lie inside the donors' convex hull.** If the treated city's features are so extreme that no convex combination can reach them (for instance, it is the largest market in the country and no donor combination can match its volume), the convex weight constraint refuses to extrapolate, and SC can only leave obvious residual imbalance and a worse pre-period fit; this is the extrapolation / convex-hull problem. The California example in Section 5.3 later illustrates exactly its boundary: ADH (2010) could synthesize California precisely because they matched per-capita cigarette sales rather than total market size, sidestepping the convex-hull obstacle of an out-of-reach volume by using a per-capita metric. Distinguish this from interpolation bias: the latter refers to the bias introduced by using a pool of donors quite different from the treated unit and interpolating over a wide range of a (possibly nonlinear) response surface, which can happen even when the treated unit lies inside the convex hull, and is mitigated by restricting the donor pool to similar units.
3. **The donor pool must be clean.** Donors must be genuinely untreated and unaffected by spillovers from the intervention on the treated unit. If the reform induces sellers in neighboring cities to migrate across cities or diverts traffic, the donors are contaminated and no longer a clean counterfactual. This is interference, where treatment spreads from the treated unit to the control units.
4. **No anticipation.** The treated unit must not react early before $T_0$ out of foreseeing the reform, otherwise the pre-period is contaminated, $T_0$ is misdated, and the loadings that get fit are dirty too.
5. **The pre-period must be long enough.** The loadings are pinned down by the pre-period path, and if the pre-period is too short, the weights merely fit noise rather than signal, leading to overfitting.

To sum up: sharp RDD uses continuity to identify the jump at the cutoff as a local effect, and its one behavioral premise, no-manipulation, leaves two testable shadows in density and covariates; RDD is inherently local, identifying only the effect at the cutoff; fuzzy RD is IV at the cutoff, identifying the LATE for compliers, pointing straight back to Chapter 12; SC, under a factor model, indirectly matches unobserved factor loadings by reproducing the pre-period path, and so tolerates time-varying confounding that DiD cannot, but it requires a good pre-period fit, the treated unit inside the convex hull, clean donors, no anticipation, and a long enough pre-period. All of these assumptions concern the counterfactual, very few can be tested directly, and what can be tested is only the shadows they cast.

## 4 Estimation

Identification has clarified what each of the two designs is estimating; this section covers how to estimate it. The spine is still "where the previous method fails is what gives rise to the next."

### 4.1 RDD: From Naive OLS to Local Linear

The laziest RDD estimate is a regression with the running variable, plugging $D_i = \mathbf{1}[X_i \geq c]$ straight in:

$$Y_i = \beta_0 + \tau D_i + X_i \beta + \varepsilon_i.$$

It has two flaws. First, it imposes a single global straight line on the entire score range, when the relationship between log GMV and score is not a straight line at all. Second, even if it were a straight line, $\tau$ still only has a clean interpretation at $X = c$. The most dangerous consequence comes from the first flaw: if the true relationship between $Y(0)$ and $X$ is curved and you fit it with a straight line, then the line's systematic deviation near the cutoff will be misread as a discontinuity. The warning here is blunt: otherwise you will mistake nonlinearity for a discontinuity.

The fix is to replace the global straight line with a local, flexible fit. The semiparametric way to write it is

$$Y_i = f(X_i) + \tau D_i + \varepsilon_i,$$

with $f$ fit by a flexible form. The workhorse estimator of modern RDD is local linear regression: use only observations within a bandwidth $h$ of the cutoff, fit a straight line on each side, and take the difference of the two lines' intercepts at the cutoff. Formally it is two weighted least squares problems, with a kernel $K(\cdot)$ giving smaller weight to points farther from the cutoff:

$$\big(\hat\alpha_+, \hat\beta_+\big) = \arg\min_{\alpha,\beta}\sum_{i: X_i \geq c}\big(Y_i - \alpha - \beta(X_i - c)\big)^2 K\!\Big(\tfrac{X_i - c}{h}\Big),$$

$$\big(\hat\alpha_-, \hat\beta_-\big) = \arg\min_{\alpha,\beta}\sum_{i: X_i < c}\big(Y_i - \alpha - \beta(X_i - c)\big)^2 K\!\Big(\tfrac{X_i - c}{h}\Big),$$

and the point estimate is $\hat\tau_{\mathrm{RD}} = \hat\alpha_+ - \hat\alpha_-$. Why a local "line" rather than a local "average"? Because at the boundary of the support, a local average has severe bias. The cutoff is the boundary of the data on each side, and if $f$ has a slope there, an estimate that only takes averages will fold the level difference caused by the slope into the jump, the so-called kernel boundary problem. Fitting a local straight line makes the slope explicit and estimated, which removes this first-order boundary bias, and this is why local linear beats local constant. As for the kernel, the triangular kernel $K(u) = (1 - |u|)\mathbf{1}[|u| \leq 1]$ is MSE-optimal at the boundary and is the default in `rdrobust`; using the uniform kernel degenerates to ordinary OLS within the window.

### 4.2 Bandwidth Selection and Bias-Corrected Robust Confidence Intervals

The bandwidth $h$ is the most consequential tuning parameter in local linear. It embodies a classic bias-variance tradeoff: the larger the bandwidth, the more points included and the smaller the variance, but the farther the points are from the cutoff, the less representative of the cutoff they are and the larger the bias; a smaller bandwidth does the reverse. The leading bias of the boundary local-linear estimate is $O(h^2)$ and the variance is $O(1/(nh))$, so the asymptotic MSE is approximately

$$\mathrm{MSE}(h) \approx \underbrace{C_B\, h^4\, f''(c)^2}_{\text{bias}^2} + \underbrace{\frac{C_V\, \sigma^2}{n\, h}}_{\text{variance}},$$

minimized at $h^\ast \propto n^{-1/5}$. The approach of Imbens and Kalyanaraman (2012) is to estimate the unknown curvature $f''(c)$ and error variance in the formula from the data, substitute them back into the MSE-optimal formula, and obtain a data-driven plug-in bandwidth $h_{\mathrm{IK}}$.

The trouble is that doing ordinary inference at exactly the MSE-optimal bandwidth is invalid. The interval $\hat\tau_{\mathrm{RD}} \pm z_{1-\alpha/2}\, \widehat{\mathrm{se}}$ assumes by default that the point estimate is asymptotically unbiased, but at the MSE-optimal bandwidth the leading bias and the standard error are of the same order, the center of the $t$-statistic is off, and the coverage of the interval falls short of the nominal level. Calonico, Cattaneo and Titiunik (2014) provide two things to fix this. One is their own MSE-optimal bandwidth (`mserd`). The other is bias-corrected robust inference: first estimate the leading bias $\widehat{\mathrm{Bias}}$ using a higher-order polynomial on an auxiliary bandwidth and subtract it from the point estimate, then fold the extra uncertainty from estimating that bias itself into the variance,

$$\text{CI}_{\text{robust}} = \big(\hat\tau_{\mathrm{RD}} - \widehat{\mathrm{Bias}}\big) \pm z_{1-\alpha/2}\sqrt{\widehat V_{\text{bc}}}, \qquad \widehat V_{\text{bc}} = \widehat V(\hat\tau) + \widehat V(\widehat{\mathrm{Bias}}) - 2\widehat{\mathrm{Cov}}(\hat\tau, \widehat{\mathrm{Bias}}).$$

This is the Robust line in the `rdrobust` output, and the current default standard. Note one feature of it: the center of the robust interval is the bias-corrected estimate, not the ordinary point estimate, so the ordinary point estimate occasionally falls outside the robust interval; this is not a bug, it is by design, and Section 6 will run into an example of it.

### 4.3 Do Not Use High-Order Global Polynomials

Some people, to "control the curvature of $f$," will stuff a fourth-, sixth-, or even higher-order global polynomial into the regression. Gelman and Imbens (2019) systematically argue why you should not, giving three reasons. First, a high-order global polynomial assigns extreme and bizarre weights to observations far from the cutoff, so the estimate at the cutoff is dragged around by distant data that ought to be irrelevant. Second, the point estimate is highly sensitive to the polynomial order, and there is no principled way to choose the order. Third, the resulting confidence intervals have poor coverage and are often too narrow. Beyond their argument, a related fact from numerical analysis is worth mentioning: the Runge phenomenon, where high-degree polynomials produce spurious wiggles near the boundary of the support and can conjure a fake jump at the cutoff out of thin air, which is the extreme version of "mistaking nonlinearity for a discontinuity" (this one is not among Gelman-Imbens's three, but an independent supplement). Section 6 will use numbers to demonstrate this instability: on the same data, fourth-, sixth-, and eighth-order polynomials give wildly different jump estimates, one of them even negative. The right road is local linear or local quadratic with a bandwidth, not a global high-order polynomial.

### 4.4 2SLS for Fuzzy RD

Fuzzy RD, in estimation, is 2SLS in a window around the cutoff, using $Z_i = \mathbf{1}[X_i \geq c]$ to instrument the actual treatment $D_i$, with both stages carrying the recentered running variable (a slope on each side). Copy over the two-equation structure of Chapter 12:

$$\text{first stage: } D_i = \pi_0 + \pi_1 Z_i + h^D(X_i - c) + \nu_i, \qquad \text{structural: } Y_i = \beta_0 + \beta_1 D_i + h^Y(X_i - c) + \varepsilon_i,$$

where $h^D, h^Y$ are flexible polynomials set separately on each side of the cutoff (an independent slope per side). Here $\beta_1$ is exactly $\tau_{\mathrm{FRD}}$, and $\hat\beta_1$ is the fuzzy-RD LATE at the cutoff. The first-stage coefficient $\pi_1$ is the effect of crossing the threshold on the probability of getting the badge, the jump in treatment; the coefficient on $D_i$ in the structural equation is the clean effect after instrumenting. `rdrobust` runs this in one line with the `fuzzy=` argument, internally doing local linear, the CCT bandwidth, and bias-corrected robust inference.

### 4.5 SC: Weight Optimization and Fit

The estimation core of SC is the simplex-constrained minimum-distance problem of Section 3.5. In practice it goes in a few steps: first pick the predictors (the covariates plus the values of several pre-period outcomes, with a common practice being the pre-period mean plus a few specific-month lags), then optimize the weights $\hat{\mathbf{w}}$ over the pre-period so that the synthetic city is as close as possible to the treated city on these predictors, then use this set of weights to construct the synthetic path $\sum_j \hat w_j Y_{jt}$ over the full sample, with the treated city's actual path minus it being the period-by-period gap $\hat\tau_{1t}$. There is an interesting observation here: the weights that solve this constrained quadratic problem are often sparse, with only a few donors receiving nonzero weight. Note that this sparsity does not come from the $\ell_1$ shrinkage of LASSO: on the simplex, $\sum_j |w_j| = \sum_j w_j = 1$ is always constant, so the $\ell_1$ norm does not vary with the weights at all and provides no extra shrinkage. Its true source is the vertex geometry of the constraint polytope (the simplex), where the optimum of the quadratic objective tends to land on a low-dimensional face with only a few nonzero coordinates. This is good for interpretation: the synthetic city is stitched from a few nameable donors, more transparent than a twenty-city average. In R, `Synth` and `tidysynth` are common implementations, and `MSCMT` fixes some numerical issues in `Synth`'s optimization.

### 4.6 SC Inference: Permutation and Placebo

With only one treated unit, where does the standard error come from? The answer is that there is no analytic standard error, and SC inference rests entirely on permutation and placebo.

The core is the in-space placebo. The procedure is to take each donor city in turn and pretend it is the "treated unit," run SC on it, compute a gap curve, and obtain a whole family of "fake effects." If the treated city truly has an effect, its post-period gap should look abnormally large within this family of fake effects. To quantify this "abnormality," use the post/pre MSPE ratio:

$$r_j = \frac{\mathrm{MSPE}^{\text{post}}_j}{\mathrm{MSPE}^{\text{pre}}_j}, \qquad \mathrm{MSPE}^{\text{pre}}_j = \frac{1}{T_0}\sum_{t \leq T_0}\big(Y_{jt} - \hat Y_{jt}(0)\big)^2.$$

The denominator is the pre-period fit error and the numerator is the post-period gap. A unit that truly has an effect fits well in the pre-period (small denominator) and has a large post-period gap (large numerator), so the ratio is large. Rank the treated unit's $r_1$ among all units' $\{r_j\}$, and the permutation p-value is the fraction of units with a ratio no smaller than the treated unit's, $\hat p = \frac{1}{J+1}\sum_j \mathbf{1}[r_j \geq r_1]$. Placebos that fit poorly in the pre-period to begin with will have inflated ratios, and the standard practice is to drop them before comparing.

Another diagnostic is the in-time placebo: assume the treatment date is some period before the actual reform and re-estimate. A design that holds up should show no gap around this fake, not-yet-reformed date. It tests not the size of the effect but the timing: if a gap appears before the actual reform, that says the fit is capturing something else, not the policy.

### 4.7 Two Modern Extensions

SC has two modern extensions worth knowing, which this chapter treats only at the conceptual level.

Augmented SC (Ben-Michael, Feller and Rothstein 2021) targets the case where the pre-period fit is imperfect. When the treated unit lies outside the donors' convex hull and residual imbalance remains no matter how you stitch, plain SC is biased. ASC adds a bias-correction term: use an outcome model (commonly a ridge-penalized regression) to fit that residual imbalance and subtract it out of the estimate,

$$\hat Y_{1t}^{\text{aug}}(0) = \sum_j \hat w_j Y_{jt} + \underbrace{\Big(\hat m_t(\mathbf{X}_1) - \sum_j \hat w_j \hat m_t(\mathbf{X}_j)\Big)}_{\text{correction for residual pre-period imbalance}}.$$

The ridge version of ASC allows small negative weights and limited extrapolation, trading a little interpretability for smaller bias, the same idea as the bias-corrected matching of Chapter 9.

Synthetic DiD (Arkhangelsky, Athey, Hirshberg, Imbens and Wager 2021) blends SC and DiD. It estimates two sets of weights at once: SC-style unit weights $\hat\omega_j$ (making donors match the treated unit's pre-period trend) and time weights $\hat\lambda_t$ (giving larger weight to pre-periods that look more like the post-period), then runs a weighted two-way fixed effects regression:

$$(\hat\tau, \cdot) = \arg\min_{\tau, \mu, \alpha_i, \beta_t}\sum_{i,t}\big(Y_{it} - \mu - \alpha_i - \beta_t - \tau D_{it}\big)^2\, \hat\omega_i\, \hat\lambda_t.$$

Unlike pure SC, it keeps DiD-style unit fixed effects, so it does not require an exact match in levels, only parallel trends after weighting; unlike DiD, it does not give all donors and all periods equal weight. It is doubly robust and also supports formal jackknife/bootstrap inference. The three-way comparison of these two extensions is in the table below.

| Dimension | Classic SC | Augmented SC | Synthetic DiD |
|---|---|---|---|
| Weight constraint | Nonnegative, sum to one | Small negative weights, limited extrapolation | Unit weights plus time weights |
| Unit fixed effects | None (requires level matching) | None (with bias correction) | Yes (relaxes level matching) |
| When pre-period fit is poor | Biased | Corrected by outcome model | Buffered by DiD structure |
| Inference | Permutation/placebo | Permutation plus correction | Analytic jackknife/bootstrap |
| R implementation | Synth / tidysynth | augsynth | synthdid |

To sum up: RDD estimation moves from naive global regression to local linear, using local lines to remove boundary bias, choosing the bandwidth data-drivenly with IK/CCT, giving robust intervals with bias correction, and firmly avoiding high-order global polynomials; fuzzy RD is 2SLS at the cutoff, and you must recenter. SC estimation is a simplex-constrained minimum-distance weight optimization, its lifeline in the pre-period fit, and its inference can only lean on in-space and in-time placebos; augmented SC and synthetic DiD are two modern ways out when the fit is inadequate.

## 5 Anchor Papers

A method only stands up when it lands in real research. Three anchor papers: one is the methodological history and classic vehicle of RDD, one is a leading RDD work in a platform setting, and one is the founding application of SC. Each is laid out along five elements, paper, method, data, results, and limitations, with the focus on how the assumptions get defended.

### 5.1 Thistlethwaite and Campbell (1960) and Lee (2008)

::: {.case}
The starting point of the methodological history: Thistlethwaite and Campbell (1960) first proposed the discontinuity idea, using a scholarship award score cutoff to study the effect of the award on students' later performance, the origin of RDD. The classic vehicle that turned it into a modern workhorse of causal identification is Lee (2008).

Paper: Lee, D. S. (2008). Randomized experiments from non-random selection in U.S. House elections. Journal of Econometrics, 142(2), 675-697. The question is incumbency advantage: does a party winning a House seat this time raise the probability that it wins again next time.

Method: the difficulty is that "good candidates simply get more votes," so directly comparing winners and losers mixes in candidate quality. Lee's running variable is the Democratic vote share in the previous election, with cutoff 50%. Districts that narrowly won by just clearing 50% and those that narrowly lost by just missing it are nearly randomly comparable in voter composition and candidate quality, and who lands above the line and who below depends on the chance of the last few votes. Sharp RD: crossing 50% delivers the seat.

Data: district-level data from U.S. House elections, vote share and next-election outcomes.

Results: districts that narrowly won have a significantly higher win probability and vote share for that party in the next election, so incumbency advantage is a real causal effect, not merely a reflection of candidate quality. Lee also showed this design satisfies the RDD diagnostics: the density is smooth at 50% and predetermined covariates do not jump.

Limitations: the estimate is local, identifying only the effect in evenly matched districts with vote shares close to 50%, and it does not extrapolate to safe, one-party-dominant districts. This is the classic embodiment of RDD's strong internal validity and narrow external validity.
:::

This paper's defense strategy is the RDD template: not by piling control variables into a regression, but by arguing the thing itself, that "the threshold is nearly random," with the density test and covariate smoothness as corroborating evidence.

### 5.2 Luca (2016)

::: {.case}
Paper: Luca, M. (2016). Reviews, reputation, and revenue: The case of Yelp.com (Harvard Business School Working Paper 12-016). It is a leading RDD work in a platform setting, sharing a core with Aster's badge case: the platform rounds a continuous rating to a discrete tier, manufacturing an artificial discontinuity.

Method: Yelp rounds a restaurant's average star rating to the nearest half star, so a restaurant with a true average of 4.24 displays as 4.0 stars and one at 4.26 displays as 4.5 stars, two restaurants of almost identical underlying quality yet a half star apart in displayed rating. Luca uses this rounding threshold to run RDD, comparing the revenue of restaurants just rounded up to the higher half star with those just rounded down to the lower half star, separating the causal effect of "half a star more" from "quality was already higher."

Data: Yelp.com restaurant rating data matched to restaurant revenue records from the Washington State Department of Revenue.

Results: holding underlying quality fixed, an extra half star brings about a 4.5% rise in that quarter's revenue; more generally, each additional star raises revenue by about 5% to 9%. The effect concentrates in independent restaurants, while chains are almost unaffected by stars, because consumers already have stable expectations of chain brands and do not rely on reviews to judge quality.

Limitations: the estimate is local to restaurants near the rounding thresholds; whether the effect extrapolates to other parts of the rating range requires extra assumptions; and the external validity of a single platform and single category is the usual caveat.
:::

This paper demonstrates the huge opportunity for RDD in platform data: platform systems are full of rounding, thresholds, and tiering rules, and each one is a natural discontinuity waiting to be used for clean causal identification.

### 5.3 Abadie, Diamond and Hainmueller (2010)

::: {.case}
Paper: Abadie, A., Diamond, A., & Hainmueller, J. (2010). Synthetic control methods for comparative case studies: Estimating the effect of California's tobacco control program. Journal of the American Statistical Association, 105(490), 493-505. It is the founding application of the synthetic control method and the direct prototype of this chapter's SC line.

Method: in 1988 California passed Proposition 99, raising the cigarette excise tax by 25 cents per pack along with a tobacco control program. To estimate this policy's effect on per-capita cigarette sales, the difficulty is that there is only one treated unit, California, and no single state can serve as California's counterfactual on its own. SC's solution is to use a convex combination of other states to synthesize a "synthetic California," making it reproduce California's per-capita cigarette sales and several predictors as closely as possible before 1988, so that after the policy the gap that opens between synthetic California and real California is the policy effect.

Data: a state-level panel of U.S. states from 1970 to 2000, per-capita cigarette sales and a set of predictors. The donor pool takes the remaining U.S. states, and after excluding states that during the same period launched their own large tobacco control programs or sharply raised cigarette taxes (such as Massachusetts, Arizona, Oregon, Florida), yields 38 states.

Results: synthetic California is a weighted combination of a few states, and after the policy real California's per-capita cigarette sales stay persistently below synthetic California, about 26 packs lower (per person per year) by 2000. Using the in-space placebo, taking each donor state in turn as the treated unit, California's gap looks abnormally large within the whole placebo distribution, and permutation inference is given on that basis.

Limitations: with only one treated unit, inference rests entirely on the permutation logic, with no analytic standard error; and the credibility of the conclusion hinges on the quality of the pre-period fit and the cleanliness of the donor pool.
:::

This paper set up SC's entire operational workflow, convex-combination weights, pre-period fit, and in-space placebo inference, and nearly all later SC studies follow this template, including the Aster case in Section 6.

## 6 A Full Walkthrough on Aster Data

Now run the earlier tools end to end on Aster's two design lines. The code below uses R 4.5.3, and every number cited in the text comes from the actual run output of this code. RDD uses `rdrobust` and `rddensity`, and SC uses `tidysynth`.

### 6.1 The Data-Generating Process for RDD

The design parameters are as follows: fix `set.seed(101)`, with $N = 12000$ sellers. The running variable is a continuous quality score drawn from a smooth truncated normal, guaranteeing the density is continuous at 80 (that is, no manipulation). $Y(0)$ is a smooth nonlinear function of score, rising steeply with score (high-scoring sellers simply sell more), with a sinusoidal wiggle that is smooth at the cutoff yet forces a global polynomial to bend across the whole range. The true sharp jump is 0.20.

```r
set.seed(101); N <- 12000; c_cut <- 80
score <- pmin(pmax(rnorm(N, 74, 11), 45), 100)     # smooth truncation, density continuous at cutoff
s     <- (score - c_cut) / 10
mu0   <- 5.60 + 0.95 * s + 0.20 * s^2 + 0.50 * sin(4 * s)  # Y(0): rises with score, smooth wiggle
tau_true <- 0.20                                   # true jump at cutoff
D_sharp  <- as.integer(score >= c_cut)
y_sharp  <- mu0 + tau_true * D_sharp + rnorm(N, 0, 0.30)
```

The fuzzy variant lets the badge probability jump from 0.10 to 0.85 at the threshold (not from 0 to 1), with the structural badge effect still the same 0.20:

```r
p_badge <- ifelse(score >= c_cut, 0.85, 0.10)      # badge prob jumps at 80, not 0 to 1
D_fuzzy <- rbinom(N, 1, p_badge)
y_fuzzy <- mu0 + tau_true * D_fuzzy + rnorm(N, 0, 0.30)
```

The naive "badged vs. unbadged" difference is 1.891, about 9.5 times the truth of 0.200, confirming the absurdly high number from Section 1, of which nearly nine-tenths is selection bias.

### 6.2 RDD Estimation: Local Linear, Bandwidth, Polynomials, Density, Fuzzy

Sharp RD uses `rdrobust`'s local-linear with the CCT bandwidth:

```r
library(rdrobust); library(rddensity)
rd <- rdrobust(y = y_sharp, x = score, c = 80,
               kernel = "triangular", p = 1, bwselect = "mserd")
```

The conventional point estimate is 0.200, the conventional 95% interval is [0.133, 0.267]; the bias-corrected robust 95% interval is [0.111, 0.261], with p-value 1.26e-06. The CCT/MSE-optimal bandwidth is 1.876, meaning it uses only sellers between 78.124 and 81.876, with effective sample sizes of 695 and 689 on the two sides. Local linear lands cleanly on the truth of 0.200, while the naive comparison of 1.891 is nearly ten times too high.

One reminder: in this DGP the badge effect is the same 0.20 for all sellers ($\tau_i$ does not vary with score), so here the local effect at the cutoff happens to equal the ATE over everyone, and $\tau_{\mathrm{RD}}$ coincides numerically with the ATE, which is only a convenience of the simulation for clean reconciliation. In the real world the effect almost always varies with score, and then what local-linear can cleanly recover is still only the effect at that one point, the cutoff, not the average effect of giving the badge to a seller at any arbitrary score. The local nature emphasized in Section 3.3 shows up numerically only in that case of heterogeneous effects.

Bandwidth sensitivity is a mandatory robustness check. Re-estimate at half and double the bandwidth:

```r
rd_half <- rdrobust(y = y_sharp, x = score, c = 80, kernel = "triangular", p = 1, h = 1.876 / 2)
rd_dbl  <- rdrobust(y = y_sharp, x = score, c = 80, kernel = "triangular", p = 1, h = 1.876 * 2)
```

The half bandwidth gives 0.186 and the double bandwidth gives 0.237, both floating modestly around the truth, showing the estimate is insensitive to the bandwidth, a good sign.

Now demonstrate why you cannot use a high-order global polynomial. Rescale the running variable to $[-1, 1]$ and use raw powers (so that all powers are zero at the cutoff and the coefficient on $D$ is exactly the jump), fitting fourth-, sixth-, and eighth-order:

```r
sc_c <- (score - 80) / max(abs(score - 80))        # rescale to [-1,1]
poly_jump <- function(k) coef(lm(y_sharp ~ D_sharp * poly(sc_c, k, raw = TRUE)))["D_sharp"]
poly_jump(4); poly_jump(6); poly_jump(8)
```

The three estimates are 0.556, -0.215, 0.254, swinging violently, with the sixth-order giving a meaningless negative number. This is Gelman-Imbens's warning made flesh: change the order and the estimate at the cutoff flips upside down, and can conjure a reversed fake jump out of nothing.

The manipulation test uses `rddensity`:

```r
rddensity(X = score, c = 80)                        # manipulation test
```

The p-value of the density jump is 0.250 (T = 1.151), far from significant, so the density is continuous at 80, there is no evidence that sellers pushed their scores over the line, and no-manipulation holds up.

Last is fuzzy RD, instrumenting the badge as an imperfect treatment:

```r
rd_fz <- rdrobust(y = y_fuzzy, x = score, c = 80, fuzzy = D_fuzzy,
                  kernel = "triangular", p = 1, bwselect = "mserd")
```

The first-stage jump estimate is 0.749 (truth 0.750), so crossing the threshold lifts the probability of getting the badge by about three-quarters. The reduced-form jump is 0.170. The hand-computed Wald ratio is $0.170 / 0.749 = 0.227$. `rdrobust`'s fuzzy-RD LATE point estimate is 0.241, with robust 95% interval [0.071, 0.236]. Here there is a detail to be honest about: the fuzzy estimate (0.241) is a bit higher than the sharp one (0.200), because fuzzy's auto-selected bandwidth (3.333) is wider than sharp's (1.876), sweeping in more of the wiggly region and picking up a bit of curvature bias; also note that the conventional point estimate 0.241 falls just outside the robust interval's upper bound 0.236, which is exactly the behavior Section 4.2 described of the robust interval being centered on the bias-corrected estimate rather than the conventional point estimate, not an error. Overall, whether the hand-computed ratio 0.227 or `rdrobust`'s local-linear fuzzy estimate, both land near the true structural effect of 0.20.

::: {.intuition}
This bit of roughness in fuzzy RD is not an error but its normal state. Sharp RD only needs to estimate one jump in the outcome, while fuzzy RD must estimate two jumps, in the outcome and in treatment, and then divide them, and the auto-selected bandwidth is usually wider too, so the point estimate is more easily pushed off the truth by curvature and noise. As for the conventional point estimate falling outside the robust interval, that too is by design: the robust interval is centered on the bias-corrected estimate and should not be expected to hold the conventional point estimate neatly inside. When reading a fuzzy RD conclusion, focus on the robust interval centered on the bias-corrected estimate, not on the last few decimals of the conventional point estimate. Fuzzy RD on real data mostly looks just like this.
:::

### 6.3 The Data-Generating Process for SC

The design parameters are as follows: fix `set.seed(12021)`, with 20 donor metros plus 1 treated metro (Aster-Metro), 36 months, and the reform launching in month 25, the first month after the text's $T_0$, so the post-period is months 25 through 36, 12 months in all (the `T0` variable in the code below refers to this first treated month, 25). The data are generated by a factor model, two common factors plus city-specific loadings, with very small noise ($\sigma = 0.02$). The treated city's loadings are a convex combination of four donors' loadings (numbers 3, 7, 12, 18), with weights 0.35/0.30/0.20/0.15, which guarantees it lies inside the donors' convex hull and SC can reproduce its factor exposure.

```r
set.seed(12021)
n_donors <- 20; Tt <- 36; T0 <- 25; t_idx <- 1:Tt
f1 <- 2.0 * (1 - exp(-t_idx / 12)) + 0.05 * t_idx  # common factor 1: platform adoption trend (concave)
f2 <- sin(2 * pi * t_idx / 12) + 0.02 * t_idx       # common factor 2: seasonal term
tau_path <- ifelse(t_idx >= T0, 0.12 * (1 - exp(-(t_idx - T0 + 1) / 3)), 0)  # true effect: smooth ramp
```

Because the cities have different loadings and their trends are naturally non-parallel, DiD's parallel trends fails: the treated city's pre-period slope is 0.131 and the donor average is 0.117. The naive two-period, two-group DiD (treated city against the twenty-city simple average) gives 0.244, while the true post-period average effect is only 0.095, about 2.6 times too high. The wrong-looking 0.244 from Section 1 is confirmed here.

### 6.4 SC Estimation: Weights, Fit, and the Gap

Build the synthetic control with `tidysynth`, taking as predictors the pre-period mean plus the lags at months 6, 12, 18, and 24, with weights optimized over months 1 through 24:

```r
library(tidysynth)
sc <- panel %>%
  synthetic_control(outcome = log_gmv, unit = metro, time = month,
                    i_unit = "Aster-Metro", i_time = 25, generate_placebos = TRUE) %>%
  generate_predictor(time_window = 1:24, mean_pre = mean(log_gmv)) %>%
  generate_predictor(time_window = 6,  gmv_m06 = log_gmv) %>%
  generate_predictor(time_window = 12, gmv_m12 = log_gmv) %>%
  generate_predictor(time_window = 18, gmv_m18 = log_gmv) %>%
  generate_predictor(time_window = 24, gmv_m24 = log_gmv) %>%
  generate_weights(optimization_window = 1:24) %>%
  generate_control()
```

The pre-period RMSPE is 0.017, hugging the 0.02 noise floor, an excellent fit, so the loading-matching argument holds. The post-period average gap is 0.086, which is SC's effect estimate, already quite close to the truth of 0.095 and far better than the naive DiD's 0.244; the last month's gap is 0.103, against the truth of 0.118. The donors with the largest weights in the synthetic city are Donor-10 (0.240), Donor-02 (0.155), Donor-04 (0.108), Donor-08 (0.059), Donor-12 (0.039), and Donor-13 (0.037).

::: {.warning}
A teaching point that must be spelled out: the weights SC estimates do not correspond to the four anchors used to generate the data (numbers 3, 7, 12, 18). The reason is that with only two factors plus a level, many different convex combinations can reproduce the treated city's factor exposure, so the weights are not unique. This is not a bug. What SC delivers is a good pre-period fit and a credible counterfactual, not an exact recovery of the weights, so do not read the weights as "the true control group."
:::

### 6.5 SC Placebo Inference

Inference leans on placebos. The in-space placebo takes each donor in turn as the treated unit and computes the post/pre MSPE ratio, and the treated city's ratio reaches 32.777, ranking first among all 21 units, with permutation p-value 0.048. Note that 0.048 is the strongest result 21 units can give ($1/21 = 0.048$), meaning the treated city has the single largest MSPE ratio of all. For inference with a single treated unit, all the force rests on this ranking.

The in-time placebo assumes the reform date is 8 months before the actual reform and re-estimates:

```r
sc_time <- panel %>%
  synthetic_control(outcome = log_gmv, unit = metro, time = month,
                    i_unit = "Aster-Metro", i_time = 25 - 8, generate_placebos = FALSE) %>%
  generate_predictor(time_window = 1:16, mean_pre = mean(log_gmv)) %>%
  generate_predictor(time_window = 6,  gmv_m06 = log_gmv) %>%
  generate_predictor(time_window = 12, gmv_m12 = log_gmv) %>%
  generate_weights(optimization_window = 1:16) %>%
  generate_control()
```

In the window after this fake, not-yet-reformed date and before the actual reform, the average gap is only -0.019, about zero, showing that the gap really does open only at the actual reform date and is not an artifact of the design.

![Synthetic control for Aster's single-city fee reform. Top panel: the actual log-GMV trajectory of the treated city Aster-Metro (orange) versus the synthetic city (blue dashed), with a vertical line marking the reform date $T_0$; before the reform the two lines nearly overlap (pre-period RMSPE 0.017), and after the reform a gap opens. Bottom panel: gap plot, where the treated city's gap (orange) hugs zero before the reform and opens to about 0.10 after, with the gray background lines the in-space placebo gaps from taking each donor in turn as the treated unit (placebos whose pre-period fit is more than five times worse than the treated city's have been dropped per standard practice), and the treated city's gap looks abnormally large within the whole placebo family.](assets/fig/fig_14_synth.svg)

The reconciliation of this section can be summarized as follows: on the RDD line, the naive comparison of 1.891 is inflated by selection to nearly ten times the truth of 0.200, local-linear cleanly recovers 0.200 and gives a robust interval covering the truth, the bandwidth is robust and the density does not jump, and the fuzzy variant instruments the badge to land the LATE near 0.20; on the SC line, the naive DiD of 0.244 overstates by 2.6 times because trends are not parallel, and synthetic control uses the tight pre-period fit of 0.017 to pull the estimate back to 0.086, close to the truth of 0.095, with the in-space and in-time placebos together confirming this gap as the effect of the reform.

## 7 Failure Modes and Robustness

In a simulation the identifying assumptions are built in, but in real research they can fail at any time. This section lays out the most common ways each of the two designs fails and the actionable responses.

Start with RDD. The number-one threat is manipulation and sorting. Once units can precisely control which side of the cutoff they land on, the two sides of the threshold are no longer comparable and continuity fails along with it. Diagnose it with the McCrary/`rddensity` density test to see whether the running variable bunches at the cutoff, and at the same time check whether predetermined covariates pass smoothly through the threshold. If there is obvious integer heaping at the threshold (for instance, scores artificially rounded up to 80), a common remedy is the donut RD, dropping the observations in a band $[c - \delta, c + \delta]$ right at the cutoff and using only the slightly outer points to extrapolate, sidestepping the most suspect heaping zone.

The second class of threat is bandwidth sensitivity and polynomial overfitting. The bandwidth is the linchpin of RDD results, so you must report the curve of the estimate as the bandwidth varies, and the standard checklist explicitly requires robustness to different windows. On polynomials, keep in mind Gelman-Imbens's lesson: never use a global polynomial above fourth order, and instead use local linear or local quadratic with a bandwidth; the set of estimates in Section 6 that swing violently with order and even turn negative is a living counterexample.

The third class is other policies co-jumping at the cutoff, that is, compound treatment. If some other program also hangs on the same threshold, for instance if 80 not only awards the badge but simultaneously unlocks a lower fee or more traffic slots, then the jump at the threshold mixes these several things together, and RDD estimates their joint force rather than the badge's effect alone. This is not a statistical problem but an institutional one, and can only be ruled out by pinning down exactly which treatments the threshold triggers. Adjacent to it is covariates jumping at the cutoff: if predetermined covariates are not smooth, that is a signal that the two sides of the threshold were not comparable to begin with, and should be taken seriously as falsifying evidence. Remember too RDD's local nature: the estimate is valid only for units near the cutoff, and extrapolating it far away requires additional, usually untested, homogeneity assumptions. A last purely technical but common point is plotting: RDD plots use an equal-count binned binscatter, and drawing one casually with covariates will draw it wrong, so the right practice is to treat the binscatter as a semiparametric regression $y_i = \mu(x_i) + \mathbf{w}_i'\gamma + \varepsilon_i$ (`binsreg`), which is a plotting-diagnostic issue and does not affect the estimator itself.

Now SC. The first-order failure is a poor pre-period fit. If the synthetic city cannot hug the treated city before the reform, the loading-matching argument collapses and the post-period gap cannot be interpreted. Check the pre-period RMSPE, and if the fit is poor consider augmented SC to patch the residual imbalance with an outcome model. The second class is interpolation bias and convex-hull violation: the treated unit's features are too extreme and lie outside the donors' convex hull, no convex combination can reconstruct it, and under the simplex constraint SC can only give a suboptimal fit clinging to the edge of the hull. The third class is donor pool contamination: if a donor is actually treated too, or splashed by spillovers from the intervention on the treated unit, or hit by the same shock in the same period, it is no longer clean. In a platform setting this is especially real, as the reform city's lower fees may draw in sellers from neighboring cities and prompt buyers to order across cities, sweeping donors into the same change, which is interference, treatment spreading from the treated unit to the control units. The response is to exclude the obviously entangled donors from the pool, or redefine the market boundary.

SC has two more subtle problems. One is too few or too many donors: too few and you cannot assemble a good convex combination, too many and you tend to overfit noise in the pre-period, stitching a synthetic city that resembles the treated one only in the pre-period and falls apart in the post-period. The other is overfitting the pre-period path itself, where a long and flexible fit may be matching idiosyncratic noise rather than the factor structure, and the in-time placebo is exactly the diagnostic against this: if a gap appears before the fake date, it is most likely overfitting. Finally, inference with a single treated unit can inherently only lean on permutation, with no analytic standard error, and the resolution of the permutation p-value is capped by the number of units, so in Section 6, 21 units can at best reach 0.048. When donors are too few, even this permutation inference runs out of power, and this is an inborn limitation of SC, to be admitted honestly in a paper rather than papered over with a spuriously precise-looking asterisk.

Stringing together the failure modes of the two designs, the common lesson is: their credibility lies not in the mechanics of the estimator but in the design and the institutions. RDD's success or failure hinges on "does this cutoff trigger only one thing, and can units manipulate it," while SC's hinges on "is the pre-period fit really good, and are the donors really clean." The density test, bandwidth curve, covariate smoothness, pre-period RMSPE, and in-space and in-time placebos are all evidence offered around these two questions, and none of them can substitute for a substantive judgment about the institutional background.

## 8 Further Reading

::: {.readings}
Required reading, in suggested order:

- Lee and Lemieux (2010, Journal of Economic Literature). The most readable RDD survey, covering continuity, bandwidth, and the diagnostic checklist in one sweep, the main model for this chapter's RDD half, and the first read for beginners.
- Imbens and Lemieux (2008, Journal of Econometrics). Another cornerstone of the RDD practical guide, complementary to Lee-Lemieux, with a focus on the discussion of bandwidth and functional form.
- Calonico, Cattaneo and Titiunik (2014, Econometrica). Bias-corrected robust inference and the MSE-optimal bandwidth, the theory behind `rdrobust`, the full version of this chapter's Section 4.2.
- Cattaneo, Idrobo and Titiunik (2020, Cambridge booklet). The practical manual for `rdrobust` and `rddensity`, walking step by step from plotting to estimation to testing, the operating checklist to have before writing an RDD paper.
- Abadie, Diamond and Hainmueller (2010, JASA). The founding application of synthetic control, California Proposition 99, the direct prototype of this chapter's SC line, with a focus on the in-space placebo and weight optimization.
- Abadie (2021, Journal of Economic Literature). The authoritative survey of the SC method, mapping the convex hull, loading matching, inference, and extensions into a single chart, the first read for an overview of the SC half.

Further reading:

- McCrary (2008, Journal of Econometrics). The original manipulation test, to understand why a density break is evidence of sorting.
- Gelman and Imbens (2019, Journal of Business & Economic Statistics). Why you should not use high-order global polynomials, the full argument of this chapter's Section 4.3, short and forceful.
- Imbens and Kalyanaraman (2012, Review of Economic Studies). The original data-driven optimal bandwidth, the derivation of the MSE-optimal bandwidth.
- Thistlethwaite and Campbell (1960, Journal of Educational Psychology). The historical origin of the RDD idea, worth a read to know where it came from.
- Lee (2008, Journal of Econometrics) and Luca (2016, HBS working paper). The former is the classic political-science vehicle of RDD, the latter a leading work in the platform setting, and reading the two together conveys how the same design is used across different fields.
- Abadie and Gardeazabal (2003, American Economic Review). The methodological origin of synthetic control, the effect of terrorism on the economy of the Basque Country.
- Ben-Michael, Feller and Rothstein (2021, JASA). Augmented SC, bias correction when the pre-period fit is imperfect.
- Arkhangelsky, Athey, Hirshberg, Imbens and Wager (2021, American Economic Review). Synthetic DiD, blending SC and DiD and providing formal inference.
:::

::: {.apa-refs}
- Abadie, A. (2021). Using synthetic controls: Feasibility, data requirements, and methodological aspects. *Journal of Economic Literature, 59*(2), 391-425. https://doi.org/10.1257/jel.20191450
- Abadie, A., Diamond, A., & Hainmueller, J. (2010). Synthetic control methods for comparative case studies: Estimating the effect of California's tobacco control program. *Journal of the American Statistical Association, 105*(490), 493-505. https://doi.org/10.1198/jasa.2009.ap08746
- Abadie, A., & Gardeazabal, J. (2003). The economic costs of conflict: A case study of the Basque Country. *American Economic Review, 93*(1), 113-132. https://doi.org/10.1257/000282803321455188
- Arkhangelsky, D., Athey, S., Hirshberg, D. A., Imbens, G. W., & Wager, S. (2021). Synthetic difference-in-differences. *American Economic Review, 111*(12), 4088-4118. https://doi.org/10.1257/aer.20190159
- Ben-Michael, E., Feller, A., & Rothstein, J. (2021). The augmented synthetic control method. *Journal of the American Statistical Association, 116*(536), 1789-1803. https://doi.org/10.1080/01621459.2021.1929245
- Calonico, S., Cattaneo, M. D., & Titiunik, R. (2014). Robust nonparametric confidence intervals for regression-discontinuity designs. *Econometrica, 82*(6), 2295-2326. https://doi.org/10.3982/ECTA11757
- Cattaneo, M. D., Idrobo, N., & Titiunik, R. (2020). *A practical introduction to regression discontinuity designs: Foundations*. Cambridge University Press. https://doi.org/10.1017/9781108684606
- Cattaneo, M. D., Jansson, M., & Ma, X. (2018). Manipulation testing based on density discontinuity. *The Stata Journal, 18*(1), 234-261. https://doi.org/10.1177/1536867X1801800115
- Gelman, A., & Imbens, G. (2019). Why high-order polynomials should not be used in regression discontinuity designs. *Journal of Business & Economic Statistics, 37*(3), 447-456. https://doi.org/10.1080/07350015.2017.1366909
- Imbens, G. W., & Kalyanaraman, K. (2012). Optimal bandwidth choice for the regression discontinuity estimator. *The Review of Economic Studies, 79*(3), 933-959. https://doi.org/10.1093/restud/rdr043
- Imbens, G. W., & Lemieux, T. (2008). Regression discontinuity designs: A guide to practice. *Journal of Econometrics, 142*(2), 615-635. https://doi.org/10.1016/j.jeconom.2007.05.001
- Lee, D. S. (2008). Randomized experiments from non-random selection in U.S. House elections. *Journal of Econometrics, 142*(2), 675-697. https://doi.org/10.1016/j.jeconom.2007.05.004
- Lee, D. S., & Lemieux, T. (2010). Regression discontinuity designs in economics. *Journal of Economic Literature, 48*(2), 281-355. https://doi.org/10.1257/jel.48.2.281
- Luca, M. (2016). *Reviews, reputation, and revenue: The case of Yelp.com* (Harvard Business School Working Paper No. 12-016). Harvard Business School.
- McCrary, J. (2008). Manipulation of the running variable in the regression discontinuity design: A density test. *Journal of Econometrics, 142*(2), 698-714. https://doi.org/10.1016/j.jeconom.2007.05.005
- Thistlethwaite, D. L., & Campbell, D. T. (1960). Regression-discontinuity analysis: An alternative to the ex post facto experiment. *Journal of Educational Psychology, 51*(6), 309-317. https://doi.org/10.1037/h0044319
:::
