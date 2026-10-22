---
title: "Regression Discontinuity Designs"
subtitle: "Local Effects at the Cutoff"
seriesline: "Foundations of Information Systems Economics · Chapter 21"
series: "Foundations of Information Systems Economics: Empirical Methods for Information Systems Economics"
sidebartitle: "Foundations of Information Systems Economics"
sidebarsub: "Theory · Methods · Applications"
chapterlabel: "Chapter 21 · Regression Discontinuity Designs"
---

## Introduction

The Aster marketplace awards a "Preferred" badge to sellers whose rating just clears a line. That sellers above the line sell more is no surprise: highly rated sellers may well have been better to begin with. The useful comparison is much narrower: sellers just to the left of the threshold, who narrowly missed the badge, against sellers just to the right, who narrowly got it. If nothing but the badge jumps near the threshold, the break in sales has a causal meaning.

The identifying power of a regression discontinuity design (hereafter RDD) comes from institutional detail: it rests on continuity near the threshold and speaks only to a local effect. Badges, tiers, and eligibility lines on a platform do not automatically make a discontinuity design. Following Aster's badge storyline, this chapter asks whether the threshold can be manipulated and what else jumps at it. The same platform's single-city fee reform, which has no threshold, is left to the synthetic control of Chapter 22.

The chapter proceeds as follows. Section 1 opens with two wrong numbers. Section 2 states the estimand. Section 3 takes up, in turn, sharp RD, manipulation tests, local randomization, discrete running variables and honest inference (inference whose coverage is guaranteed uniformly over a whole class of conditional mean functions), fuzzy RD, regression kinks, multiple cutoffs and geographic RD, and extrapolation away from the cutoff. Section 4 covers the theory of local polynomial estimation, along with covariates, donuts, and power. Section 5 reads three anchoring papers, Section 6 runs the toolkit on Aster's seven data blocks, and Section 7 discusses failure modes. The chapter uses the randomization inference and clustered standard errors of Chapter 8, the kernels and bandwidths of Chapter 14, the potential outcomes of Chapter 15, the LATE theorem of Chapter 18, and the parallel trends of Chapter 19.

## 1 Two Numbers That Do Not Add Up

::: {.case}
Aster is a fictional two-sided online marketplace. We use simulated data to follow one design line. The payoff is that the data-generating process is fully known, so every estimator can be checked against the truth.

This design line is the Preferred badge. Aster gives a "Preferred" badge, displayed on the storefront page, to sellers whose quality score reaches 80. Management wants to know how large the causal effect of this badge is on a seller's monthly sales (log GMV). The most natural calculation subtracts the average log GMV of unbadged sellers from that of badged sellers, and the difference is 1.891. In a growth report, it would say that the badge lifts sales by nearly 1.9 log points, enough for the team to claim credit.
:::

This number is not credible, and earlier chapters have already explained what is wrong with it. The 1.891 mixes in the selection bias of Chapter 15. High-scoring sellers were already selling more: even without any badge, stores above 80 would have far higher log GMV than stores below 80. The "badged vs. unbadged" comparison therefore tangles the badge effect together with the fact that "high-scoring sellers were simply stronger." Because this is a simulation, we know that the true causal jump of the badge at the threshold is only 0.200; the 1.891 is about 9.5 times the truth, and nearly nine-tenths of it is selection bias.

The second number comes from another rule on the same platform. Aster's search page has a "4.5 stars and up" filter. A seller's displayed star rating is the average of its 10 most recent reviews rounded to the nearest half star, so an average of 4.25 is enough to enter the filter. A researcher follows the defaults of `rdrobust`, reads the Conventional line of the output, and estimates the filter's effect on log GMV at 0.092, with 95% interval [0.041, 0.143] (the Robust line of the same output, [0.023, 0.151], only just covers the truth). Noting that the average takes only 21 values, the researcher then clusters the standard errors by score value, and the interval narrows to [0.079, 0.106]. The truth is 0.15, and neither interval covers it. This is not the luck of one sample: in the repeated simulations of Section 3.4, the conventional interval covers the truth only 78.4% of the time, and the clustered interval only 36.9%.

RDD offers a way around the first trap: do not compare all high-scoring sellers with all low-scoring ones; look only at a very narrow band around the 80-point line. Sellers at 79 and 81 barely differ in ability, category, or operations, and who just clears the line and who just misses it is nearly random. The small log-GMV difference between those just above and those just below the line is the badge's clean effect. The second number shows that when the same idea meets a discrete score, inference has to start over.

![Aster's Preferred badge discontinuity: the horizontal axis is quality score and the vertical axis is log GMV. Sellers are binned by score and averaged, with blue dots below 80 and orange dots at or above 80; each side gets a local linear fit within the optimal bandwidth of Calonico, Cattaneo and Titiunik (hereafter CCT), and the dashed line marks the cutoff. The vertical jump at the threshold is the badge effect, labeled RDD jump = 0.200 (truth 0.200). Note that log GMV rises steeply with score overall, which is exactly why the naive comparison is absurdly high, while the jump at the discontinuity is far smaller than this steep slope.](assets/fig/fig_21_rdd.svg)

This idea is the spine of the chapter. Why it holds and where it fails are the questions the next seven sections answer.

## 2 The Economic Model and the Estimand

Before estimating anything, we state which causal quantity RDD is after and the notation in which to write it down. It builds on the potential outcomes of Chapter 15.

Write the seller's quality score as the running variable (also called the forcing variable) $X_i$, and the cutoff as $c$, here $c = 80$. The treatment indicator $D_i$ records whether the seller got the badge. The defining feature of RDD is that the treatment probability jumps at $c$:

$$\lim_{x\downarrow c}\Pr(D_i = 1 \mid X_i = x) \ \neq\ \lim_{x\uparrow c}\Pr(D_i = 1 \mid X_i = x).$$

A **sharp RD** is the special case in which this jump goes exactly from 0 to 1: clearing the line guarantees the badge and missing it guarantees none, $D_i = \mathbf{1}\{X_i \geq c\}$. Each seller has two potential outcomes $Y_i(1), Y_i(0)$, and the observed outcome is $Y_i = D_i Y_i(1) + (1 - D_i) Y_i(0)$. The structural picture behind RDD is that the untreated outcome is a smooth function of the running variable and the only abrupt change at $c$ comes from treatment,

$$Y_i = Y_i(0) + \tau_i \cdot \mathbf{1}\{X_i \geq c\}, \qquad Y_i(0) = m(X_i) + \varepsilon_i, \ \ m \text{ smooth}.$$

The RDD estimand is not the average treatment effect over everyone, but the local average effect at that one point, the cutoff:

$$\tau_{\mathrm{RD}} = \mathbb{E}\big[\,Y_i(1) - Y_i(0) \mid X_i = c\,\big].$$

This has to be made plain from the outset: without extra assumptions, RDD can learn only the effect at the single point $X = c$. How much a seller on the 80-point line gains from getting the badge, it can answer; what would happen if the badge were given to a 60-point seller, it cannot. This is RDD's built-in limitation and the price it pays for strong internal validity, discussed formally in Section 3.8.

To summarize: the target of RDD is the local effect at the cutoff, $\tau_{\mathrm{RD}}$, a single number pinned to the point $X = c$. The next section covers the assumptions under which it is identified.

## 3 Identification

Identification logically comes before estimation: it asks whether, under a set of assumptions about the counterfactual, the causal quantity we want can be written as a function of the distribution of observables, and the answer depends only on the design and the institutional background. This section first uses Aster's numbers or figures to show what is to be proved, then states numbered assumptions and theorems, followed by a proof sketch that names the assumption used up at each step, and finally answers whose effect the quantity estimates, who determines the weights, and where the assumptions are most likely to fail. The order follows the direction of generalization: Sections 3.1 to 3.3 cover continuity, no-manipulation, and local randomization at a single cutoff; Section 3.4 deals with a discrete running variable; Sections 3.5 and 3.6 generalize treatment from a switch from 0 to 1 to a probability and an intensity; and Sections 3.7 and 3.8 handle multiple cutoffs, two-dimensional boundaries, and extrapolation away from the cutoff.

### 3.1 A Break on a Line: The Continuity Assumption of Sharp RDD

Start with an example anyone understands. Suppose the government rules that households with annual income below \$25,000 can collect a certain subsidy. Directly comparing households that collect it with those that do not will of course not do, because poor households consume less to begin with. But a household earning \$24,999 and one earning \$25,001 barely differ in anything except that one just qualifies and the other just misses. Within a very narrow neighborhood of the threshold, who lands above the line and who below is nearly a matter of chance, so the consumption difference between the thin layers of households on either side can be attributed cleanly to the subsidy. The whole intuition of RDD fits in one sentence: near the cutoff, treatment is as good as randomly assigned. The same holds for Aster's sellers at 79 and 81. Turning this intuition into a formal identifying assumption rests on continuity.

::: {.assumption}
**A1 (continuity; the identification condition of Wooldridge 2010 §21.5.1)** The conditional expectation functions of both potential outcomes are continuous at the cutoff:

$$x \mapsto \mathbb{E}[Y_i(0) \mid X_i = x] \ \text{ and } \ x \mapsto \mathbb{E}[Y_i(1) \mid X_i = x] \ \text{ are continuous at } x = c.$$
:::

What this assumption constrains is worth reading closely. If the badge were taken away, the seller's log GMV as a function of score should pass smoothly through the 80-point line, not sprout a step out of nowhere at exactly 80. It allows $\mathbb{E}[Y(0) \mid X]$ to rise steeply, to bend, and to be nonlinear, and it allows high-scoring and low-scoring sellers to be worlds apart; the one thing it forbids is a break unrelated to the badge at the cutoff itself. Precisely because of continuity, the sellers just left of the cutoff who narrowly missed the line are a legitimate counterfactual for the sellers just right of it who narrowly cleared it.

::: {.theorem}
**Theorem 1 (sharp RD identification; the identification result of Hahn, Todd and Van der Klaauw 2001, with a checkable statement in Wooldridge 2010 equation (21.98))** Suppose treatment is determined by $D_i = \mathbf{1}\{X_i \geq c\}$, $X_i$ has positive density on both sides near $c$, and A1 holds. Then

$$\tau_{\mathrm{RD}} = \lim_{x\downarrow c}\mathbb{E}[Y_i \mid X_i = x] - \lim_{x\uparrow c}\mathbb{E}[Y_i \mid X_i = x] = \mathbb{E}[\tau_i \mid X_i = c].$$
:::

The proof sketch has three steps. Step one (the sharp design): $Y_i = \mathbf{1}\{X_i < c\}Y_i(0) + \mathbf{1}\{X_i \geq c\}Y_i(1)$, so for $x \geq c$, $\mathbb{E}[Y_i \mid X_i = x] = \mathbb{E}[Y_i(1) \mid X_i = x]$, and for $x < c$ it equals $\mathbb{E}[Y_i(0) \mid X_i = x]$ (Wooldridge equation (21.97)); positive density guarantees that the conditional means on both sides are defined near $c$ and can be estimated. Step two (A1): approaching from the right gives $\mathbb{E}[Y_i(1) \mid X_i = c]$, and approaching from the left gives $\mathbb{E}[Y_i(0) \mid X_i = c]$. Step three (subtraction, no new assumption): the difference of the two is $\mathbb{E}[Y_i(1) - Y_i(0) \mid X_i = c]$. In a sharp design overlap fails at every point, since the treatment probability is 0 below the cutoff and 1 above it; propensity score weighting is therefore of no use, and identification rests entirely on taking the limit (Wooldridge 2010 §21.5.1).

In words: sharp RD estimates the average effect for the units sitting exactly at the cutoff, with weights determined by the composition of sellers at the point $X_i = c$, which leaves the researcher no choice; away from that point it says nothing. In Aster this quantity is the badge effect for 80-point sellers, with true value 0.200; the local linear estimate of Section 6.2 is also 0.200, while the naive comparison's 1.891 counts the whole steep slope. A1 is most likely to fail in two situations: the 80-point line triggers some other treatment at the same time, or sellers can push their scores precisely over the line (Section 3.2). When it fails, the difference in limits estimates the badge effect plus whatever else jumps at 80, and the two cannot be separated.

::: {.intuition}
The continuity assumption has a reassuring property: it shifts the burden of identification from "I controlled for all the confounders correctly" to "nothing else jumps at this one point." You do not need to observe the seller's ability, category, or operational quality; you only need to believe that, apart from the badge, no second thing switches at exactly 80, which can often be checked in a plot and argued from the institutions. But continuity itself is not testable: it constrains the counterfactual $\mathbb{E}[Y(0) \mid X = c]$, and at that point we see only the badged side.
:::

### 3.2 No-Manipulation: An Assumption That Leaves Shadows

Behind continuity there is usually a concrete behavioral premise: sellers cannot precisely control which side of the threshold they land on. This is no-manipulation, also called no-sorting. Back in the subsidy example, if households could precisely under-report their income to just below \$25,000, the area just below the threshold would fill up with households whose income was actually higher and who deliberately lowered their filings; they would no longer be comparable to the households just above the threshold, and the density of the running variable would break at the threshold, bunching on one side. The corresponding threat at Aster is sellers pushing their scores precisely up to 80 through fake orders or category changes. Writing "no precise manipulation" as an assumption requires introducing the seller's unobservable type $U_i$.

::: {.assumption}
**A2 (no precise manipulation; Conditions 2b and 1c of Lee 2008)** For every type of seller $u$, the conditional density of the score $f_{X \mid U = u}(x)$ is continuous at $x = c$; given the type, potential outcomes vary smoothly with the score.
:::

A2 allows sellers to influence their own scores, through hard work or by soliciting good reviews, as long as no type of seller can place itself precisely on one side of 80. It is not necessary for identification (Card, Lee, Pei and Weber 2015, footnote 8), but it is easier to argue than A1, and it implies A1. This is the approach of Lee (2008), which footnote 9 of Card, Lee, Pei and Weber (2015) records as Lee's Proposition 2(a): when the density of every type is continuous, the conditional distribution of types is continuous at $c$. The next step is this book's own, written on that basis: by Bayes' rule, $\mathbb{E}[Y_i(d) \mid X_i = x] = \int \mathbb{E}[Y_i(d) \mid X_i = x, U_i = u]\, f_{X \mid U = u}(x)\, dF_U(u) / f_X(x)$, and since the outcome and the density of every type are continuous at $c$, the whole conditional mean is continuous. The assumption also has two testable implications: $f_X$ is continuous at $c$, and the conditional distribution of covariates fixed before treatment, $\mathbf{W}_i$, is continuous at $c$ (Corollaries 1 and 2 of Card, Lee, Pei and Weber 2015 give the kink version of the same argument). The original test of McCrary (2008) first bins the data, then smooths on each side, and tests the difference in the log densities on the two sides:

$$\theta = \ln\Big(\lim_{x\downarrow c} f_X(x)\Big) - \ln\Big(\lim_{x\uparrow c} f_X(x)\Big), \qquad H_0:\ \theta = 0.$$

The statistic $\hat\theta$ replaces the two one-sided limits with density estimates on the two sides. The modern default is the `rddensity` of Cattaneo, Jansson and Ma, which needs no pre-binning.

::: {.theorem}
**Theorem 2 (local polynomial density estimation and the manipulation test; Cattaneo, Jansson and Ma 2020 Theorem 1 and §4)** Suppose the subsamples to the left and right of $c$ each satisfy the following: in a one-sided neighborhood of $c$, the distribution function of the score on that side is $p + 1$ times continuously differentiable, with a positive one-sided density (the paper's Assumption 1 applied to each side's subsample; under the alternative the overall distribution function is not differentiable at $c$, so the condition can only be stated side by side), and the kernel is nonnegative, symmetric, and continuous on $[-1, 1]$ (Assumption 2); this theorem does not use A2 of this book. Let $n$ be the sample size and $h$ the bandwidth (which may differ between the two sides, written $h_-$ and $h_+$). On each side of $c$, run a local polynomial regression of order $p$ on that side's empirical distribution function; the first-order coefficients are the density estimates $\hat f_{X,+}(c)$ and $\hat f_{X,-}(c)$. If $nh^2 \to \infty$ and $nh^{2p+1} = O(1)$, the density estimates are asymptotically normal at both interior and boundary points, with leading bias of order $h^p$. If, further, $n\min\{h_-^2, h_+^2\} \to \infty$ and $n\max\{h_-^{1+2p}, h_+^{1+2p}\} \to 0$ (negligible bias), the statistic $T_p(h)$, the difference of the two estimates divided by the pooled standard error, satisfies, under $H_0: f_{X,+}(c) = f_{X,-}(c)$, $\lim \Pr\big(|T_p(h)| \geq z_{1-\alpha/2}\big) = \alpha$, and its rejection probability tends to 1 under the alternative.
:::

Here $\alpha$ is the significance level. The proof sketch has two steps (the paper's supplementary material). Step one (the paper's Assumptions 1 and 2): on one side of $c$, the empirical distribution function is an unbiased estimator of $F$; running a local polynomial regression on it, the first-order coefficient approximates the derivative of $F$, that is, the density. Because only data from one side are used, the estimator adapts to the boundary automatically, with no need to bin first. Step two (Assumption 1 holding side by side, with nonoverlapping samples on the two sides): the two density estimates are asymptotically independent, and their difference divided by the pooled standard error is asymptotically standard normal under the null; `rddensity` sets $p = 2$ by default and uses $T_3$ for the bias-corrected robust test.

In words: the density test targets an implication of A2, not A1 itself. It is neither necessary nor sufficient: a continuous density does not imply continuous potential outcomes, and the density may break for reasons unrelated to the outcome, such as integer heaping or administrative rules. But bunching on one side of the threshold is strong evidence of sorting. Call the badge data of Section 1 block R1 (the blocks are listed in the table of Section 6.1). There `rddensity` gives a p-value of 0.250 (T = 1.151), with no sign of manipulation. Block R1h is the same market a year later: the rating formula has been made public, scores are reported to one decimal place, and some sellers between 76 and 80 push their reported scores to exactly 80.0. As a result, 460 sellers crowd onto 80.0 (against only 23 at 79.9 and 30 at 80.1), and the p-value of the density test drops to $5.25 \times 10^{-5}$; the same test on the unobserved true scores gives a p-value of 0.460. Alongside the density stands covariate balance: in block R1c (a fresh sample from the same world, with two predetermined covariates), the jump at 80 in log GMV in the year before the badge launched is -0.029 (robust p-value 0.222), and the jump in tenure on Aster is 0.137 (p-value 0.321), both indistinguishable from zero. Although the core assumption cannot be tested directly, it leaves two shadows: no jump in the density and no jump in predetermined covariates. The exclusion restriction of Chapter 18 leaves no such shadow at all.

### 3.3 Local Randomization: Treating the Neighborhood of the Cutoff as a Small Experiment

The continuity framework requires only that the conditional means join up seamlessly at the cutoff, and estimation must extrapolate to $c$ from both sides. The other route is local randomization: assume directly that within a very narrow window around the cutoff, $\mathcal{W}_0 = [c - w_0, c + w_0]$, who lands above the line and who below is nearly random, as if the thin layer of units on either side of the threshold had taken part in a small randomized experiment. Aster's block R2 shows best what this route is good for. The average rating takes only 21 values between 3.0 and 5.0; 4.2 displays as 4.0 stars, while 4.3 displays as 4.5 stars and enters the filter. The two groups of sellers differ by only one step in the average, and comparing them is like comparing the two arms of an experiment.

::: {.assumption}
**A3 (local randomization; Cattaneo, Frandsen and Titiunik 2015 Assumption 1)** There is a window $\mathcal{W}_0 = [c - w_0, c + w_0]$ such that for all units with $X_i \in \mathcal{W}_0$: (a) the conditional distribution of the score within the window is the same for every unit; (b) potential outcomes depend on the score only through the treatment indicator within the window, not on the particular value of the score, and are not affected by the scores of units outside the window.

**A4 (local SUTVA and a locally constant effect; Assumptions 2 and 3 of the same paper)** The outcome of a unit in the window depends only on its own treatment status, and $Y_i(1) = Y_i(0) + \tau$.
:::

The $\tau$ in A4 is the constant treatment effect within the window; a bare $\tau$ in this chapter is always a treatment effect parameter, never a quantile level. A3 is much stronger than continuity: within the window the conditional mean is not merely continuous but does not vary with the score at all, and the paper describes it as an approximation to the continuity condition. It allows outcomes outside the window to vary arbitrarily with the score; the price is that the window is usually very narrow and the sample small, so inference must use finite-sample exact methods. For choosing the window, the criterion in §3 of Cattaneo, Frandsen and Titiunik (2015) is to run balance tests on predetermined covariates and take the largest window within which no nested window is rejected (Assumptions 4 and 5 of the same paper give the conditions under which this procedure is valid). The paper's algorithm starts from the full sample and narrows step by step; `rdwinselect` works the other way, widening step by step from the narrowest window and taking the window just before the first rejection. Both use the same criterion.

::: {.theorem}
**Proposition 1 (Fisher exact test and test inversion; written by this book from Cattaneo, Frandsen and Titiunik 2015 §2.2 and §2.3)** Let the window contain $n_{\mathcal{W}}$ units, of which $n_{\mathcal{W}}^{+}$ are above the cutoff. The randomization mechanism is fixed margins (the number of treated units in the window is fixed, and only which units are treated is random): the treatment vector $\mathbf{Z}$ is uniformly distributed over all 0/1 vectors with $n_{\mathcal{W}}^{+}$ entries equal to 1, its $i$th component $Z_i$ is the treatment indicator of unit $i$, $T(\mathbf{z})$ is the difference in means within the window, and $\mathbf{z}^{\mathrm{obs}}$ is the observed treatment vector.

(i) Under A3 and the sharp null $H_0$ (the potential outcomes of every unit in the window do not vary with the treatment vector $\mathbf{z}$), $p_{\mathrm{F}} = \Pr\big(|T(\mathbf{Z})| \geq |T(\mathbf{z}^{\mathrm{obs}})|\big)$ is the exact Fisher p-value, valid for any sample size in the window; when it is approximated with $R$ random permutations, $\hat p_{\mathrm{F}} = \big(1 + \#\{r : |T(\mathbf{z}^{(r)})| \geq |T(\mathbf{z}^{\mathrm{obs}})|\}\big)/(1 + R)$ is a valid p-value.

(ii) Under A3 and A4, substituting $Y_i - \tau^{\circ} Z_i$ into (i) tests $\tau = \tau^{\circ}$, and all $\tau^{\circ}$ that are not rejected form an exact confidence set for $\tau$.
:::

The proof sketch has three steps. Step one (the premise of the proposition): given the number of treated units in the window, every combination of which $n_{\mathcal{W}}^{+}$ units land above the line is equally likely; this is the randomization mechanism that §2.1 of Cattaneo, Frandsen and Titiunik (2015) posits separately. A3(a) makes it plausible but does not imply it; exchangeable assignment across units is also needed. Step two (A3(b) and the sharp null): outcomes do not depend on treatment, so the observed outcomes are the same under every combination, and the distribution of $T(\mathbf{Z})$ is fully determined by the uniform distribution of step one, with no large-sample approximation needed; the observed combination is itself a legitimate draw, and counting it in both numerator and denominator gives the valid form, which can never equal 0. Step three (A4): at the true value, $Y_i - \tau^{\circ} Z_i = Y_i(0)$ does not vary across combinations, so applying (i) to the adjusted outcomes gives the test, and inverting the test gives the confidence set.

In words: local randomization estimates the average effect for the finite population in the window, with equal weight on every unit in the window, and the window is delimited by covariate balance tests. Its point of failure is also the window: as soon as $Y(0)$ tilts with the score inside the window, the difference between the two groups picks up that stretch of slope.

At Aster, `rdwinselect`, with seller tenure as the balance covariate, selects the window [4.2, 4.3]. By the criterion above it is the largest window that is not rejected, and here it also happens to be the narrowest (balance test p-value 0.215; the next window, [4.1, 4.4], already drops to 0.016). The window holds 355 sellers on the left and 377 on the right. The difference in means is 0.1295, and the 95% interval from test inversion is [0.0875, 0.1725]; none of 1000 permutations reaches the observed value, so the valid p-value for the sharp null is 1/1001. The truth of 0.15 lies inside the interval, but the window's estimand is in fact the difference in means between the two points 4.3 and 4.2, which equals 0.1596 in the DGP, 0.0096 more than the jump at the cutoff. This 0.0096 has two sources, each accounting for about half: the effect itself has already risen by 0.0048 from 4.25 to 4.3, and $Y(0)$ also rises by 0.0048 from 4.2 to 4.3, the latter being exactly the curve that A3 requires to be flat within the window. In 1000 repeated simulations, the Neyman interval within the fixed window (the large-sample interval built in Chapter 15 from Neyman's conservative variance) covers 0.15 with frequency 0.923.

Widen the window and the slope creeps in: the difference in means is 0.175 on [4.1, 4.4] and 0.212 on [4.0, 4.5], and the latter is rejected under the sharp null $\tau = 0.15$ (p-value at most 1/1001). When there are few observations near the cutoff or the running variable is discrete, this route is often safer than the continuity framework, which relies on smooth extrapolation; the randomization inference of Chapter 8 has found a new stage here.

### 3.4 Discrete Running Variables: Why Clustering Fails and How to Construct Honest Intervals

Block R2 is also the source of the second wrong number of Section 1. The seller's average rating $X_i$ has a support $\mathrm{Supp}(X_i)$ of only $N_{\mathrm{sp}} = 21$ points, and the single point 5.0 holds 18.6% of the sellers. Filter results are sorted by average rating, and sellers above 4.4 begin to crowd onto the first page, so the effect keeps rising to the right of the cutoff. The continuity framework fits a line on each side and extrapolates to 4.25 (the local linear estimator and the robust interval of CCT are defined formally in Sections 4.1 and 4.2; here we use their conclusions in advance), but within the conventional bandwidth of 0.370 only 4 support points remain on each side, and the curvature sits precisely at the scale of the grid of support points. Assumption 1(a) of Calonico, Cattaneo and Titiunik (2014) requires the score to be continuously distributed, and their Remark 1 states explicitly that this rules out a discrete running variable. The consequence of discreteness is that the extrapolation bias does not vanish as the sample grows: however many sellers there are, they all fall on these 21 points, so the center of the conventional interval is systematically off the truth. Lee and Card (2008) recommend clustering the standard errors by score value, absorbing the specification error into the variance as a random effect; Kolesár and Rothe (2018) prove that this route does not work and offer an alternative.

::: {.assumption}
**A5 (bounded second derivative; the function class of Kolesár and Rothe 2018 §5.1, arXiv version v4)** On each side of the cutoff, the first derivative of the conditional mean function $m(x) = \mathbb{E}[Y_i \mid X_i = x]$ satisfies a Lipschitz condition (the change in the derivative is at most $M$ times the change in the argument): for any $a, b$ on the same side of the cutoff, $|m'(a) - m'(b)| \leq M|a - b|$, with the constant $M \geq 0$ fixed in advance by the researcher.
:::

Here $M$ is a bound on the second derivative of the conditional mean function; it is the same kind of bound as, but a different object from, the $M$ in the HonestDiD of Chapter 19, which bounds the second differences of the trend difference. Kolesár and Rothe (2018) write it as $K$, and the corresponding argument of `RDHonest` (Kolesár 2026) is called `M`. A5 constrains curvature, not levels: it allows the conditional mean to jump at the cutoff and to slope arbitrarily on either side, and requires only that the slope not change too fast. The paper offers a way to read it: under A5, over an interval of length one the function deviates from the straight line connecting the endpoints by at most $M/8$.

::: {.theorem}
**Theorem 3 (inference with a discrete running variable; Kolesár and Rothe 2018 Theorem 1 and Proposition 1, numbered as in arXiv version v4)** Let $(Y_i, X_i)$, $i = 1, \ldots, n$, be i.i.d., with $X_i$ taking finitely many values, and let the local linear estimator $\hat\tau$ with bandwidth $h$ use $N_h$ observations (the paper's Assumptions 1 and 2 give the moment and identifiability conditions).

(i) If the number of support points within the bandwidth does not change with the sample size and, in the framework of §4 of Kolesár and Rothe (2018), the difference between the conditional mean and the polynomial specification shrinks at rate $N_h^{-1/2}$, then the variance estimator clustered by score value (cluster-robust variance, CRV), $\hat\sigma^2_{\mathrm{CRV}}$, does not converge to a constant but to a nondegenerate random limit; under homoskedasticity, the relative difference between the expectation of this limit and the true variance equals a positive term that grows with $\sqrt{N_h}$ times the specification error, plus a negative term that depends only on the distribution of the score. If the specification error is fixed rather than shrinking, the asymptotic coverage of the CRV interval is 0 or 1 (the paper's footnote 10).

(ii) Under A5, let $\hat\sigma_{\mathrm{NN}}$ be the nearest-neighbor variance estimator, $r_{\sup}$ the ratio of the largest conditional bias of $\hat\tau$ over the A5 function class to its standard error, and $\mathrm{cv}_{1-\alpha}(r)$ the quantile of $|\mathcal{N}(r, 1)|$ at level $1 - \alpha$. Then the interval

$$\Big[\hat\tau \pm \mathrm{cv}_{1-\alpha}(r_{\sup})\,\hat\sigma_{\mathrm{NN}}/\sqrt{N_h}\Big]$$

is honest over the whole A5 function class: the infimum of the coverage probability over this class is asymptotically no less than $1 - \alpha$, that is, $\liminf_{n \to \infty}\inf_{m}\Pr_m(\tau \in \mathrm{CI}) \geq 1 - \alpha$, where the infimum runs over all conditional mean functions satisfying A5 (the paper's equation (5.1)).
:::

The proof sketch has three steps (§4 and §5.1 of Kolesár and Rothe 2018, and its Appendices A and B). Step one (i.i.d. sampling and the moment conditions of the paper's Assumption 2): split the $t$-statistic into two terms. Let $\bar\tau = \mathbb{E}[\hat\tau \mid X_1, \ldots, X_n]$ be the conditional expectation of $\hat\tau$ given the scores; $(\hat\tau - \bar\tau)/\mathrm{se}$ is asymptotically standard normal by the central limit theorem, and $(\bar\tau - \tau)/\mathrm{se}$ is the ratio of bias to standard error. Step two ((i)): the clustered variance treats the specification error at each support point as a random effect, and the number of clusters is the number of support points within the bandwidth; with a fixed number of clusters, as with the few clusters of Chapter 8, the variance estimator is itself random; when the specification error is not large relative to the noise, the negative term dominates, and clustered standard errors are on average even smaller than EHW (Eicker-Huber-White, heteroskedasticity-robust) standard errors. Step three ((ii), A5): the absolute value of the bias term is at most $r_{\sup}$; the worst-case function is $-(M/2)(x - c)^2$ to the right of the cutoff and $(M/2)(x - c)^2$ to the left, with second derivatives $-M$ and $M$, exactly on the boundary of A5 (the paper's Appendix B.1), and substituting the local linear weights gives a closed form for $r_{\sup}$. The quantiles of $|\mathcal{N}(r, 1)|$ increase in $r$, and $r_{\sup}$ is the worst case over the whole function class, so with its critical value, coverage holds uniformly for every function in the class. This "uniformly" is the key: for any one fixed smooth function, the robust interval of CCT also covers asymptotically when the score is continuous, but pointwise coverage does not imply uniform coverage, and honesty requires the latter. The coverage guarantee is conditional on the values of the score, the bandwidth may be fixed, and discrete and continuous scores are treated alike.

In words: the honest interval still estimates the jump at the cutoff, but its coverage guarantee is conditional on $M$. Armstrong and Kolesár (2020) point out in §3.3.1 of their arXiv version (`arXiv:1606.01200v6`) that to remain honest over the whole function class, $M$ must be chosen in advance and cannot be data driven. The abstract of Armstrong and Kolesár (2018) (arXiv version v4) gives a more general result: over centrosymmetric convex function classes, data-driven tuning cannot yield intervals that are substantially narrower while maintaining coverage over the whole class. And §5.1 of the arXiv version v4 of Kolesár and Rothe (2018) adds that the data can provide only a lower bound on $M$, never an upper bound.

Aster's 1000 repeated simulations turn this passage into numbers. With the true bound of the DGP, $M = 7.79$, the honest interval covers 0.15 with frequency 0.998, at the cost of an average length of 0.234, about 2.14 times that of the conventional interval; on the frozen sample the interval is [-0.031, 0.202], which covers the truth and also covers zero. If each replication instead sets the bound to the one-sided lower confidence limit, computed by `RDSmoothnessBound`, for the lower bound on $M$, coverage falls to 0.779; using the point estimate of this lower bound gives coverage of 0.932, still below the nominal level, because it estimates a lower bound and never once reaches the true bound (on the frozen sample the two are 0.134 and 1.069).

For comparison, the coverage is 0.784 for the conventional interval, 0.832 for the bias-corrected robust interval, and only 0.369 for the interval clustered by score. Clustering does worst because there are only 8 clusters within the bandwidth: the clustered standard error of 0.0070 is only about a quarter of the conventional standard error of 0.0259, exactly the direction (i) describes ((i) approximates the specification error here, which is fixed and not large, by an error that shrinks with the sample). With a discrete running variable, bias correction cannot rescue coverage, and clustering makes matters worse; what belongs in a paper is an honest interval that states the $M$ used and reports a sensitivity analysis over $M$.

![Discrete running variable: the average rating at the 4.5-star filter. The left panel is one sample of 8,000 sellers; dots are sample means at each support point, crosses are true means, and the local linear fits on the two sides underestimate the jump at 4.25. The right panel shows six 95% intervals on the same sample against the truth of 0.15, with the coverage rates from 1000 repeated simulations in the right column: conventional 78.4%, robust bias-corrected 83.2%, clustered by support point 36.9%, honest interval under the true bound of the DGP 99.8%, bias-aware interval using as $M$ the one-sided lower confidence limit of the data-based lower bound on $M$ 77.9% (since $M$ is not fixed in advance, it does not count as honest by this chapter's definition), and local randomization 92.3% (the repeated simulations use the Neyman interval within a fixed window).](assets/fig/fig_21_discrete.svg)

### 3.5 Fuzzy RD: IV at the Cutoff

So far we have assumed that clearing the line guarantees the badge. Aster's actual situation is this: some sellers who reach 80 do not activate the badge, while a few sellers below 80 obtain it through appeals or special operational approval. The treatment probability at the threshold jumps not from 0 to 1 but from a lower level to a higher one; this is fuzzy RD. The jump in the outcome at the threshold reflects only the sellers whose badge status the threshold actually changes, and to convert it into "how much log GMV each change in badge status buys," you divide the jump in the outcome by the jump in treatment, which is exactly the Wald ratio:

$$\tau_{\mathrm{FRD}}(c) = \frac{\displaystyle\lim_{x\downarrow c}\mathbb{E}[Y_i \mid X_i = x] - \lim_{x\uparrow c}\mathbb{E}[Y_i \mid X_i = x]}{\displaystyle\lim_{x\downarrow c}\Pr(D_i = 1 \mid X_i = x) - \lim_{x\uparrow c}\Pr(D_i = 1 \mid X_i = x)} = \frac{\text{jump in outcome}}{\text{jump in treatment}}.$$

This shape should look familiar: fuzzy RD moves the IV framework of Chapter 18 onto the cutoff. The crossing indicator $Z_i = \mathbf{1}\{X_i \geq c\}$ is the instrument, the numerator is the reduced form of crossing on log GMV, and the denominator is the first stage of crossing on badge status. Let $D_i(1)$ and $D_i(0)$ be the seller's badge status when crossing and when not crossing the line.

::: {.assumption}
**A6 (relevance; the key assumption of Wooldridge 2010 §21.5.2)** The treatment probability jumps at $c$, that is, the inequality of Section 2 holds.

**A6′ (independent gains; Wooldridge 2010 §21.5.2)** Given $X_i$, the individual gain $\tau_i = Y_i(1) - Y_i(0)$ is independent of $D_i$.

**A7 (monotonicity at the cutoff; Wooldridge 2010 §21.5.3 only states the conclusion and refers to Imbens and Lemieux 2008; this is a local version written by this book following the LATE theorem of Chapter 18)** $D_i(1) \geq D_i(0)$; the shares of the three types of sellers (monotonicity rules out defiers), and the conditional means of the potential outcomes within each type, are continuous at $x = c$.
:::

A6 simply says that the first stage is nonzero. A6′ allows whether a seller has the badge to be correlated with $Y_i(0)$, for example because sellers who were stronger to begin with activate it more often, but it rules out self-selection on gains. A7 allows effects to vary across individuals, and the type shares and within-type means may take any level (as long as they are continuous at the cutoff); it rules out only defiers, who would go without the badge precisely because they crossed the line.

::: {.theorem}
**Theorem 4 (fuzzy RD identification; the identification result of Hahn, Todd and Van der Klaauw 2001, with a checkable statement in Wooldridge 2010 equation (21.106) and §21.5.3)** (i) If A1, A6, and A6′ hold, then $\tau_{\mathrm{FRD}}(c) = \mathbb{E}[\tau_i \mid X_i = c]$. (ii) If A6 and A7 hold, then $\tau_{\mathrm{FRD}}(c) = \mathbb{E}[\tau_i \mid X_i = c,\ D_i(1) > D_i(0)]$, the LATE of the compliers at the cutoff.
:::

Proof sketch. Part (i) takes two steps (Wooldridge 2010, p. 958). Step one (A6′): $\mathbb{E}[Y_i \mid X_i = x] = \mathbb{E}[Y_i(0) \mid X_i = x] + \Pr(D_i = 1 \mid X_i = x)\,\mathbb{E}[\tau_i \mid X_i = x]$. Step two (A1 and A6): take limits from both sides; A1 makes both the $Y(0)$ term and the effect term continuous, so they cancel, and the jump in the outcome equals the jump in the treatment probability times $\mathbb{E}[\tau_i \mid X_i = c]$; A6 guarantees that we can divide. Part (ii) takes two steps, matching Section 3.4 of Chapter 18 step by step. Step one (the continuity part of A7): to the right of the cutoff, the observed outcome is a mixture of the $Y(1)$ of always-takers and compliers and the $Y(0)$ of never-takers; to the left, it is a mixture of the $Y(1)$ of always-takers and the $Y(0)$ of compliers and never-takers; the type shares and within-type means are continuous at $c$. Step two (monotonicity): on subtraction the contributions of always-takers and never-takers cancel, the jump in the outcome equals the compliers' share times their effect, and the jump in the treatment probability is exactly the compliers' share; dividing gives the result.

In words: (i) identifies the average effect for all sellers at the cutoff, not only for compliers; if the sellers who take the badge are precisely those with large gains, A6′ fails, and as long as A7 holds, the ratio falls back on the reading of (ii). In sharp RD everyone at the threshold is a complier; under (ii) the Wald ratio of fuzzy RD reflects only the sellers "who got the badge because they crossed the 80-point line and would not have gotten it otherwise." It is a doubly local quantity, local at the cutoff and local among compliers, with weights determined by the composition of compliers at the cutoff. In Aster's fuzzy variant the badge probability jumps from 0.10 to 0.85 at the threshold; the first stage is estimated at 0.749, the reduced form at 0.170, and the Wald ratio is 0.227. In the DGP the effect is 0.20 for every seller, and the compliers' LATE is also 0.20; the `rdrobust` estimate of Section 6.2 is 0.241. The two numbers differ because the hand-computed ratio uses a separate MSE-optimal bandwidth for its numerator and its denominator (2.04 for the reduced form, 4.78 for the first stage), whereas the fuzzy estimate of `rdrobust` uses the same bandwidth, 3.333, for both. A7 is most likely to fail when crossing the line changes the terms of the badge: if the badge above the line comes with stricter shipping obligations while the specially approved badge below the line does not, some sellers would rather not activate the badge after crossing the line, yet would appeal for special approval when below it; these sellers are defiers. Let the shares of compliers and defiers at the cutoff be $s_{\mathrm{c}}$ and $s_{\mathrm{d}}$ and their average effects $\tau_{\mathrm{c}}$ and $\tau_{\mathrm{d}}$; the ratio then becomes $(s_{\mathrm{c}}\tau_{\mathrm{c}} - s_{\mathrm{d}}\tau_{\mathrm{d}})/(s_{\mathrm{c}} - s_{\mathrm{d}})$, a weighted combination of the two types' effects in which compliers get a weight greater than 1 and defiers a negative weight. Even if everyone's effect is positive, the ratio can be negative, and it can no longer be read as a LATE.

::: {.warning}
When fuzzy RD is carried out as 2SLS, there is one detail you must hold onto, a close cousin of the forbidden regression of Chapter 18: the running variable must enter the regression as $X_i - c$ rather than $X_i$, with slopes allowed to differ across the two sides. Without recentering, the coefficient on the instrument $Z_i$ gets tangled up with the intercept of the polynomial at $X_i = 0$, and what you estimate is a mixture of the jump at the threshold and the level of the polynomial at the origin; only after recentering can the coefficient on $Z_i$ be read as the jump at the cutoff itself.
:::

### 3.6 Regression Kink: A Break in the Slope

Some platform rules make the slope, rather than the level, jump at the threshold. Aster gives large sellers a commission discount: sellers whose average monthly GMV in the previous quarter is at most \$1,500 pay a commission rate of 12%; above that, the rate falls by 0.6 percentage points for every additional \$100, bottoming out at 6% at \$2,500. As a function of GMV, the commission rate is continuous at \$1,500 but bends there. If commissions affect sales, the curve of log GMV against previous-quarter GMV should also show a kink at \$1,500, and the ratio of the two kinks is the marginal effect of the commission. This is the regression kink design (hereafter RKD). Let $D_i = b(X_i)$ be the commission rate a seller faces, where $b(\cdot)$ is the policy function, unrelated to the balancing score $b(X)$ of Chapter 15. The estimand is

$$\tau_{\mathrm{RK}} = \frac{\displaystyle\lim_{x\downarrow c}\frac{d}{dx}\mathbb{E}[Y_i \mid X_i = x] - \lim_{x\uparrow c}\frac{d}{dx}\mathbb{E}[Y_i \mid X_i = x]}{\displaystyle\lim_{x\downarrow c} b'(x) - \lim_{x\uparrow c} b'(x)}.$$

::: {.assumption}
**A8 (a valid sharp RKD; Card, Lee, Pei and Weber 2015 Assumptions 1 to 4)** The outcome is generated by $Y_i = y(D_i, X_i, U_i)$, where $U_i$ is the seller's unobservable type, of unrestricted dimension. (1) The support of $U_i$ is bounded, $y$ is continuous, and its partial derivative with respect to the first argument, $y_1$, is continuous; (2) the partial derivative with respect to the second argument, $y_2$, is continuous; (3) $b(\cdot)$ is known and continuous everywhere, continuously differentiable except at $c$, $\lim_{x\downarrow c} b'(x) \neq \lim_{x\uparrow c} b'(x)$, and the types with positive score density near $c$ have positive measure; (4) the conditional score density of each type, $f_{X \mid U = u}(x)$, and its derivative with respect to $x$ are continuous near $c$.
:::

Condition (2) is much weaker than an exclusion restriction: previous-quarter GMV may affect this quarter's sales directly, as long as the effect is smooth. Condition (3) requires $b(\cdot)$ to be continuous at $c$; if the commission rate also jumped in level, the design would become an RD. Condition (4) is the key to RKD: in RD, a continuous density for every type (A2) is sufficient for identification, and footnote 8 of the paper says it is not necessary; in RKD a continuous density alone is not enough, and the slope of the density must be continuous too, because the slope of the outcome mixes in a contribution from the slope of the density.

::: {.theorem}
**Theorem 5 (RKD identifies a weighted marginal effect; Card, Lee, Pei and Weber 2015 Proposition 1)** Under A8,

$$\tau_{\mathrm{RK}} = \mathbb{E}\big[y_1(b(c), c, U_i) \mid X_i = c\big] = \int y_1(b(c), c, u)\,\frac{f_{X \mid U = u}(c)}{f_X(c)}\,dF_U(u),$$

that is, the average of the marginal effects of the types of sellers at the cutoff, weighted by the relative density with which each type falls near $c$; the paper calls it the TT (treatment on the treated) parameter.
:::

The proof sketch has three steps (§2.2 of Card, Lee, Pei and Weber 2015). Step one (Bayes): $\mathbb{E}[Y_i \mid X_i = x] = \int y(b(x), x, u)\,f_{X \mid U = u}(x)/f_X(x)\,dF_U(u)$. Step two (A8 (1), with regularity conditions that allow differentiation under the integral): differentiating with respect to $x$ gives three terms: the first is $b'(x)$ times a weighted average of $y_1$, the second is a weighted average of $y_2$, and the third comes from the change in the weights themselves with $x$. Step three (A8 (2) and (4)): the last two terms are continuous at $c$ and cancel when the left and right limits are subtracted; only the $b'(x)$ in the first term jumps at $c$, so the jump in the slope of the outcome equals the jump in $b'$ times the weighted average of $y_1$, and A8 (3) guarantees a nonzero denominator.

In words: RKD does not estimate the commission elasticity of any single seller but a weighted average of the marginal effects of the types of sellers at the cutoff, with weights that equal the relative density with which each type falls near \$1,500 and are determined by the market. Aster's block R3 has 60,000 sellers, and `rdrobust(deriv = 1, p = 2)` estimates that each 1 percentage point of additional commission leaves log GMV lower by 0.052 (conventional standard error 0.0116), with 95% robust interval [-0.083, -0.018] and true value -0.050. Local linear ($p = 1$) gives -0.033, identical to the hand-computed ratio of the slope differences on the two sides; its interval [-0.072, -0.012] also covers the truth, but the point estimate is too small by more than a third. The bias of a slope at the boundary is one order larger than that of a level, which is why the kink result of Calonico, Cattaneo and Titiunik (2014) (their Theorem 2) estimates slopes by local quadratic regression, and in practice kinks default to $p = 2$. A8 is most likely to fail on (4): if sellers can push their quarterly GMV roughly, but not precisely, past \$1,500 to get the discount, each type's density will kink there. The term from the slope of the weights in step three of the proof then no longer cancels across the two sides, and the estimator returns the commission effect plus a density-slope term: the sum of the types' outcome levels weighted by the jumps in their relative density slopes at \$1,500, divided by the jump in the commission slope. If the pushing is precise, the density jumps outright at \$1,500, the weights differ across the two sides, and the outcome level may jump as well. Corollary 1 of Card, Lee, Pei and Weber (2015) gives a testable implication: under a valid sharp RKD the score density is continuously differentiable; whether the outcome has a level jump at \$1,500 is a test for precise sorting, and at Aster this jump is 0.0066, only 1.35 standard errors.

![Regression kink: Aster's commission discount. The top panel is the commission rate, 12% up to \$1,500 and then 0.6 percentage points lower for every \$100; the bottom panel shows mean log GMV in \$25 bins with local quadratic fits within the bandwidth on each side, and the dotted line extends the left-side fit past \$1,500. The two lines bend to different slopes at the threshold, and the ratio of the kinks is the kink estimate.](assets/fig/fig_21_kink.svg)

### 3.7 Multiple Cutoffs and Geographic Boundaries

Aster's badge threshold in fact varies by category: 70 for vintage, 80 for home, and 90 for electronics. The most common practice is to normalize the score as $\tilde X_i = X_i - c_i$ and estimate the three categories together as a single RD with cutoff 0. What does this pooled number estimate? Let the cutoff faced by unit $i$ be $c_i$, with set of values $\mathrm{Supp}(c_i) = \{70, 80, 90\}$, and let $\tau_{\mathrm{RD}}(c) = \mathbb{E}[Y_i(1) - Y_i(0) \mid X_i = c, c_i = c]$ be the effect, at their own cutoff, for the units facing cutoff $c$.

::: {.assumption}
**A9 (sharp RD with multiple cutoffs; Cattaneo, Keele, Titiunik and Vazquez-Bare 2016 Assumptions 1 to 3)** For every $c \in \mathrm{Supp}(c_i)$: treatment jumps from 0 to 1 at $X_i = c_i$; $\mathbb{E}[Y_i(d) \mid X_i = x, c_i = c]$ is continuous at $x = c$ ($d = 0, 1$); and the score density of the units facing $c$, $f_{X \mid c_i = c}(x)$, is positive and continuous at $x = c$.
:::

A9 allows the categories to differ in their effects, score distributions, and numbers of sellers; it requires only that each category satisfy the sharp RD conditions at its own threshold and that there be sellers near the threshold.

::: {.theorem}
**Theorem 6 (the weights of the normalized pooled estimator; Cattaneo, Keele, Titiunik and Vazquez-Bare 2016 Lemma 1)** Under A9, the RD estimand at 0 for the normalized score $\tilde X_i$ is

$$\tau_{\mathrm{pool}} = \sum_{c \in \mathrm{Supp}(c_i)} \varpi(c)\,\tau_{\mathrm{RD}}(c), \qquad \varpi(c) = \frac{f_{X \mid c_i = c}(c)\Pr(c_i = c)}{\sum_{c'} f_{X \mid c_i = c'}(c')\Pr(c_i = c')}.$$
:::

The proof sketch has three steps (the paper's supplementary material S3). Step one (iterated expectations): $\mathbb{E}[Y_i \mid \tilde X_i = x] = \sum_c \mathbb{E}[Y_i \mid X_i = c + x, c_i = c]\Pr(c_i = c \mid \tilde X_i = x)$. Step two (the density condition of A9 and Bayes): $\Pr(c_i = c \mid \tilde X_i = x)$ is proportional to $f_{X \mid c_i = c}(c + x)\Pr(c_i = c)$ and, as $x \to 0$, converges from both sides to the same $\varpi(c)$. Step three (the continuity and sharp conditions of A9): for each $c$, the difference between the two one-sided limits of the conditional mean is $\tau_{\mathrm{RD}}(c)$, by the same argument as for Theorem 1, and pulling the weights outside gives the result.

In words: the pooled RD estimates not a simple average of the separate effects but an average weighted by "the score density at the cutoff times the share of units facing that cutoff," with weights set by how many sellers of each category sit near their own threshold; the researcher does not choose them.

Aster's block R4 has 24,000 sellers. The true effects at the three cutoffs are 0.40, 0.25, and 0.10, estimated at 0.404, 0.263, and 0.093. The sample weights are the sample version of the Lemma 1 weights, namely each category's share among all sellers within 4.80 points of their own cutoff (the MSE-optimal bandwidth of the pooled RD): 0.567, 0.306, and 0.127. Weighting by these, the Lemma 1 weighted estimator is 0.322 with 95% interval [0.294, 0.361]; running the pooled RD directly on the normalized score gives 0.293 with interval [0.223, 0.369]. Both target the truth weighted by the sample weights, 0.316; with the population weights 0.571, 0.307, and 0.122, the population target of Theorem 6 is 0.317. The simple average of the three estimates is only 0.254 (computed from the unrounded estimates; the truth is 0.25). The vintage category, which has the largest effect, has many sellers packed near 70, so it takes more than half of the weight. If management asks "what is the badge worth," 0.29 and 0.25 answer two different questions.

At Aster, A9 is most likely to fail when manipulation appears in one category, for instance if some electronics sellers can push their scores past 90, so that this category's density is no longer equal on the two sides of the threshold. If the sharp and continuity conditions of A9 still hold and only the density condition fails, the weights in step two of the proof of Theorem 6 converge from the two sides to different limits $\varpi_+(c)$ and $\varpi_-(c)$. Writing $m_{0,c}(c) = \mathbb{E}[Y_i(0) \mid X_i = c, c_i = c]$, the pooled RD estimates

$$\sum_c \varpi_+(c)\,\tau_{\mathrm{RD}}(c) + \sum_c \big[\varpi_+(c) - \varpi_-(c)\big]\, m_{0,c}(c),$$

that is, the average of the effects under the right-side weights, plus a sum of each category's $Y(0)$ level weighted by the difference between the two sides' weights. The categories already differ in sales levels at the same score (vintage is 0.25 below home, and electronics 0.30 above home), and their thresholds sit at different scores, so the latter term will not be small. But sorting usually destroys continuity as well: those who push past 90 tend to be stronger sellers, the conditional mean of this category no longer joins up across the threshold, and the pooled RD estimates only $\sum_c \big[\varpi_+(c)\lim_{x\downarrow c}\mathbb{E}[Y_i \mid X_i = x, c_i = c] - \varpi_-(c)\lim_{x\uparrow c}\mathbb{E}[Y_i \mid X_i = x, c_i = c]\big]$, a difference in outcomes between different populations on the two sides, from which not even the categories' $\tau_{\mathrm{RD}}(c)$ can be separated.

![Multiple cutoffs: the left panel shows the three categories' effect estimates and robust intervals at their own cutoffs, with their weights; the gray line is the true effect at each score, and the triangle is the vintage effect extrapolated to 75 (Section 3.8, conventional interval with a fixed bandwidth). The right panel compares the simple average of the three estimates, the Lemma 1 weighted estimator, and the pooled RD on the normalized score; the dashed line is the truth weighted by the sample weights, 0.316.](assets/fig/fig_21_multi.svg)

Another multidimensional threshold is a geographic boundary. Aster offers same-day delivery within 25 km of a fulfillment center; neighboring sellers just inside and just outside the zone are similar in other respects, so the boundary of the zone is a two-dimensional cutoff. Write the seller's location as $(X_{1i}, X_{2i})$ (kilometers east and north), the same-day zone as $\mathcal{A}$, and its boundary as $\mathcal{B}$, with $Z_i = \mathbf{1}\{(X_{1i}, X_{2i}) \in \mathcal{A}\}$; a point on the boundary is written $\mathbf{x} = (x_1, x_2)'$, a point in the plane rather than a value of a score. Multi-score designs have the same structure: when each of two scores has a threshold and treatment requires clearing both, $\mathcal{A}$ is the quadrant enclosed by the two threshold lines, $\mathcal{B}$ is its boundary, and the results below apply just the same.

::: {.assumption}
**A10 (geographic RD; Keele and Titiunik 2015 Assumptions 1 and 2)** (1) For every point on the boundary, $\mathbf{x} \in \mathcal{B}$, the conditional means of both potential outcomes given location are continuous as the location approaches $\mathbf{x}$ from any direction; (2) compound treatment is irrelevant: other treatments that switch at the same boundary do not affect potential outcomes, which depend only on the target treatment.
:::

A10 (1) is the two-dimensional version of A1 and allows the effect to vary along the boundary; (2) rules out other treatments switching at the boundary, which can be argued only from the institutions.

::: {.theorem}
**Theorem 7 (the geographic treatment effect curve; Keele and Titiunik 2015 Proposition 1)** If every unit inside the zone is treated and every unit outside is untreated, the density of locations is positive near the boundary, and A10 (1) holds, then for every $\mathbf{x} \in \mathcal{B}$,

$$\tau_{\mathrm{RD}}(\mathbf{x}) = \mathbb{E}[Y_i(1) - Y_i(0) \mid (X_{1i}, X_{2i}) = \mathbf{x}]$$

equals the difference between the limit of the conditional mean of the outcome as the location approaches $\mathbf{x}$ from inside the zone and its limit as the location approaches $\mathbf{x}$ from outside; A10 (2) guarantees that this difference is the effect of the target treatment.
:::

The proof sketch has two steps (the paper's Appendix C, which notes that it follows the argument of Hahn, Todd and Van der Klaauw). Step one (sharp): inside the zone the observed conditional mean is the conditional mean of $Y(1)$, and outside it is that of $Y(0)$. Step two (A10 (1)): from whichever direction $\mathbf{x}$ is approached, the limits on the two sides are the conditional means at $\mathbf{x}$ of $Y(1)$ and $Y(0)$, respectively, and subtracting gives the result. This mirrors the three steps of Theorem 1, with "left and right" replaced by "outside and inside the zone."

In words: geographic RD identifies an effect curve along the boundary, one number per point, and the question of weights arises only when the curve is aggregated. In Aster's block R4g, buyers on the east side value same-day delivery more, and the true effect falls from 0.30 due east to 0.10 due west. `rdms` takes six points on the boundary, with estimates 0.310, 0.271, 0.177, 0.289, 0.194, and 0.234 against true values 0.30, 0.25, 0.15, 0.10, 0.15, and 0.25; the interval at the due-west point, [0.160, 0.462], does not cover 0.10, which is sampling noise (across fresh draws, coverage at this point is close to the nominal rate), while the other five cover their true values. Pooling with the signed distance to the boundary as a one-dimensional score gives 0.204 with interval [0.171, 0.230]; this weights the boundary points by the density of sellers near each, and since sellers in this DGP are nearly uniform in direction, the target is the boundary average under a uniform distribution of directions, 0.20. The naive difference in means between inside and outside the zone is 0.322, which mixes in an urban effect: demand is stronger closer to the center. A10 (2) matters especially in geographic designs: if the delivery zone happens to coincide with a city boundary, tax rates and school districts also switch at the boundary, and geographic RD estimates their combined force.

![Geographic RD: the left panel shows a sample of 600 sellers, the 25 km same-day zone, and the six numbered points on the boundary; the right panel shows the robust 95% interval at each point, with crosses marking the true values; the interval at point 4 (due west) does not cover the truth (sampling noise).](assets/fig/fig_21_geo.svg)

### 3.8 RDD Is a Local Parameter: What Extrapolation Requires

RDD's strong internal validity comes with a cost there is no getting around: it learns only the effect at the cutoff. That a seller on the 80-point line gains $\tau_{\mathrm{RD}}$ from the badge tells you nothing directly about what would happen if the badge went to a 90-point or a 60-point seller. To extrapolate, you must additionally assume that the effect is homogeneous over the range of scores, or that it varies smoothly with the score in some specified form, and RDD itself neither provides nor tests these assumptions. This is of a piece with the extrapolation anxiety of LATE in Chapter 18: there you learn the effect for compliers, here the effect for units at the threshold. The stronger the identification, the narrower the extrapolation. Multiple cutoffs offer a principled way out.

Write $\tau_{\mathrm{RD}}(x, c) = \mathbb{E}[Y_i(1) - Y_i(0) \mid X_i = x, c_i = c]$: the first argument is the value of the score and the second is the cutoff the unit faces, so $\tau_{\mathrm{RD}}(c) = \tau_{\mathrm{RD}}(c, c)$. All vintage sellers above 70 have the badge, so their effect at 75, $\tau_{\mathrm{RD}}(75, 70)$, has no visible control group; all home sellers below 80 lack the badge, so their $Y(0)$ at 75 is visible. Write $m_{d,c}(x) = \mathbb{E}[Y_i(d) \mid X_i = x, c_i = c]$, and let $m_c(x) = \mathbb{E}[Y_i \mid X_i = x, c_i = c]$ be the observable conditional mean.

::: {.assumption}
**A11 (continuity and constant bias; Cattaneo, Keele, Titiunik and Vazquez-Bare 2021 Assumptions 1 and 2)** For a low cutoff $\underline{c}$ and a high cutoff $\bar{c}$: (1) $m_{d,c}(x)$ is continuous on $[\underline{c}, \bar{c}]$ ($d = 0, 1$); (2) the difference between the two groups' $Y(0)$ does not vary with the score: $m_{0,\underline{c}}(x) - m_{0,\bar{c}}(x) = m_{0,\underline{c}}(\underline{c}) - m_{0,\bar{c}}(\underline{c})$ for all $x \in (\underline{c}, \bar{c})$.
:::

A11 (2) is the cross-sectional version of the parallel trends of Chapter 19: the two groups' $Y(0)$ curves may sit at different heights, as long as they are parallel between 70 and 80. It restricts the control state, not the effect, so the badge effects of vintage and home sellers may be completely different. The paper stresses that extrapolating the other way, from the high cutoff downward, would require parallel treated states, which amounts to assuming that the effect does not change between the cutoffs.

::: {.theorem}
**Theorem 8 (multi-cutoff extrapolation; Cattaneo, Keele, Titiunik and Vazquez-Bare 2021 Theorem 1)** Under A11, for any $\bar x \in (\underline{c}, \bar{c})$,

$$\tau_{\mathrm{RD}}(\bar x, \underline{c}) = m_{\underline{c}}(\bar x) - \Big[m_{\bar{c}}(\bar x) + \lim_{x\uparrow \underline{c}} m_{\underline{c}}(x) - m_{\bar{c}}(\underline{c})\Big].$$
:::

The proof sketch has two steps. Step one (A11 (1) and the sharp design): at $\bar x$ the low-cutoff group is entirely treated, $m_{\underline{c}}(\bar x) = m_{1,\underline{c}}(\bar x)$, and the high-cutoff group is entirely untreated, $m_{\bar{c}}(\bar x) = m_{0,\bar{c}}(\bar x)$; at $\underline{c}$ the $Y(0)$ of the low-cutoff group is identified by the left limit and the $Y(0)$ of the high-cutoff group is directly visible, so at $\underline{c}$ the difference between the two groups' $Y(0)$ can be estimated. Step two (A11 (2)): constant bias carries this difference over to $\bar x$, $m_{0,\underline{c}}(\bar x) = m_{\bar{c}}(\bar x) + \lim_{x\uparrow \underline{c}} m_{\underline{c}}(x) - m_{\bar{c}}(\underline{c})$, and subtracting it from $m_{1,\underline{c}}(\bar x)$ gives the result.

In words: what is extrapolated is the effect for the low-cutoff group of sellers in the interior of the interval, and the entire cost rests on a single untestable assumption, constant bias, which can be falsified only by checking that the two groups' $Y(0)$ are parallel below the cutoff. At Aster the true effect for vintage sellers at 75 is 0.325, and the estimator of Cattaneo, Keele, Titiunik and Vazquez-Bare (2021) gives 0.349 with 95% interval [0.314, 0.384]. It is built from four local linear fits with a triangular kernel and a fixed bandwidth $h = 4$, and it uses the conventional interval rather than the robust bias-corrected one; with bandwidths of 3 and 5 it gives 0.345 and 0.347. Without the constant-bias correction, subtracting home's untreated outcome directly from vintage's treated outcome gives only 0.084, because at 70 the estimated gap between the two groups' $Y(0)$ is already -0.265 (true value -0.25): the vintage category sells less to begin with. When constant bias fails, the estimator returns $\tau_{\mathrm{RD}}(\bar x, \underline{c})$ plus the change in the gap between the two groups' $Y(0)$ from $\underline{c}$ to $\bar x$, $[m_{0,\underline{c}}(\bar x) - m_{0,\bar c}(\bar x)] - [m_{0,\underline{c}}(\underline{c}) - m_{0,\bar c}(\underline{c})]$. The DGP of R4 satisfies constant bias by construction: the categories share the same $Y(0)$ curve up to a constant level, so 0.349 against 0.325 shows only that the estimator works when the assumption holds, not that the assumption holds. In a real market the most likely failure is that vintage sellers' sales rise with the score faster than home sellers', say because vintage buyers care more about ratings; then the -0.265 at 70 would be smaller by the time the score reaches 75, and the extrapolation would count that shrinkage as effect.

Angrist and Rokkanen (2015) take a different route: they assume that, given a set of covariates, the conditional mean no longer depends on the running variable, so conditional comparisons can be made far from the cutoff on either side. This needs no second cutoff; the price is an assumption of the selection-on-observables kind.

To sum up this section: Theorem 1 says that under continuity sharp RD identifies the average effect at the cutoff, and Theorem 2 and the balance tests provide the two shadows of no-manipulation; Proposition 1 treats the neighborhood of the cutoff as a small experiment, at the price that $Y(0)$ must be flat within the window; Theorem 3 says that with a discrete running variable clustered standard errors are unreliable, while the honest interval guarantees coverage conditional on $M$; Theorem 4 says that fuzzy RD identifies the average effect for everyone at the cutoff under independent gains and the LATE of compliers under monotonicity, the latter pointing back to Chapter 18; Theorem 5 says that a kink identifies a density-weighted marginal effect; and Theorems 6 to 8 say that a pooled multi-cutoff RD is an average weighted by densities and shares, that geographic RD yields a curve, and that extrapolation relies on constant bias. All of these assumptions concern the counterfactual; very few can be tested directly, and all that can be tested are the shadows they cast.

## 4 Estimation

Identification has made clear what RDD is estimating; this section covers how to estimate it. The spine is still "where the previous method fails is what gives rise to the next."

### 4.1 RDD: From Naive OLS to Local Linear

The laziest RDD estimate puts $D_i = \mathbf{1}\{X_i \geq c\}$ straight into a regression with the running variable:

$$Y_i = \beta_0 + \tau D_i + X_i \beta + \varepsilon_i.$$

In this chapter's regressions the coefficient on $D_i$ is written $\tau$ and is a treatment effect parameter; the chapter never uses $\tau$ for a quantile level. This regression imposes a single global straight line on the whole score range, when the relationship between log GMV and score is not a straight line at all; and even if it were, $\tau$ would have a clean interpretation only at $X = c$. The greatest danger is this: if the true relationship between $Y(0)$ and $X$ is curved, the line's systematic deviation near the cutoff will be misread as a discontinuity, that is, nonlinearity will be mistaken for a discontinuity. The fix is a local, flexible fit; the semiparametric way to write it is

$$Y_i = m(X_i) + \tau D_i + \varepsilon_i,$$

with $m$ fit by a flexible form. The workhorse of modern RDD is local linear regression: use only observations within a bandwidth $h$ of the cutoff and fit a weighted straight line on each side, with a kernel $K(\cdot)$ giving smaller weight to points farther from the cutoff (the kernels and bandwidths of Chapter 14):

$$\big(\hat\alpha_+, \hat\beta_+\big) = \arg\min_{\alpha,\beta}\sum_{i: X_i \geq c}\big(Y_i - \alpha - \beta(X_i - c)\big)^2 K\!\Big(\tfrac{X_i - c}{h}\Big),$$

$$\big(\hat\alpha_-, \hat\beta_-\big) = \arg\min_{\alpha,\beta}\sum_{i: X_i < c}\big(Y_i - \alpha - \beta(X_i - c)\big)^2 K\!\Big(\tfrac{X_i - c}{h}\Big),$$

and the point estimate is $\hat\tau_{\mathrm{RD}} = \hat\alpha_+ - \hat\alpha_-$. Why a local line rather than a local average? The cutoff is the boundary of the data on each side, and if $m$ has a slope there, an estimate that only takes averages folds the level difference caused by the slope into the jump; this is the kernel boundary problem (Section 3.3 of Chapter 14). Fitting a local straight line estimates the slope explicitly and removes this first-order boundary bias. The triangular kernel $K(u) = (1 - |u|)\mathbf{1}\{|u| \leq 1\}$ is MSE-optimal at the boundary and is the default in `rdrobust`; with the uniform kernel the estimator reduces to ordinary OLS within the window.

### 4.2 Bandwidth Selection and Bias-Corrected Robust Confidence Intervals

The bandwidth $h$ is the most consequential tuning parameter in local linear regression, and it embodies a classic bias-variance tradeoff: the larger the bandwidth, the more points are included and the smaller the variance, but the farther the points are from the cutoff, the less they represent it and the larger the bias. The leading bias of the boundary local linear estimate is $O(h^2)$ and the variance is $O(1/(nh))$, so the asymptotic MSE is approximately

$$\mathrm{MSE}(h) \approx \underbrace{C_{\mathrm{B}}\, h^4\, m''(c)^2}_{\text{bias}^2} + \underbrace{\frac{C_{\mathrm{V}}\, \sigma^2}{n\, h}}_{\text{variance}},$$

where $m''(c)^2$ is shorthand for a combination of the one-sided second derivatives $m_+''(c)$ and $m_-''(c)$ on the two sides of the cutoff, which determine the constants in Theorem 9 (i). The MSE is minimized at $h^\ast \propto n^{-1/5}$. Imbens and Kalyanaraman (2012) estimate the unknown curvature $m''(c)$ and the error variance from the data and substitute them back into the optimal formula, obtaining the plug-in bandwidth $h_{\mathrm{IK}}$; the `mserd` of Calonico, Cattaneo and Titiunik (2014) is an improved version of the same idea. In Aster's block R1 this bandwidth is 1.876, so only sellers between 78.124 and 81.876 are used. Whether the intervals built with it are correct depends on the theorem below.

::: {.assumption}
**A12 (Calonico, Cattaneo and Titiunik 2014 Assumptions 1 and 2)** In a neighborhood of the cutoff: $\mathbb{E}[Y_i^4 \mid X_i = x]$ is bounded, and the density of $X_i$ is continuous and bounded away from zero (this condition rules out a discrete running variable, the paper's Remark 1); $m_-(x) = \mathbb{E}[Y_i(0) \mid X_i = x]$ and $m_+(x) = \mathbb{E}[Y_i(1) \mid X_i = x]$ have $S$ continuous derivatives; the conditional variances on both sides are continuous and bounded away from zero; and the kernel is bounded, nonnegative, and positive and continuous on its support.
:::

A12 allows the conditional means to bend arbitrarily as long as they are smooth enough; it requires the score to be continuously distributed with positive density near the cutoff, which is exactly what a discrete score violates.

::: {.theorem}
**Theorem 9 (MSE expansion and robust bias-corrected inference; Calonico, Cattaneo and Titiunik 2014 Lemma 1 and Theorem 1)** Under A12:

(i) If $S \geq p + 1$, $h \to 0$, and $nh \to \infty$, the conditional MSE of the local polynomial RD estimator of order $p$ is a squared bias of order $h^{2(p+1)}$ plus a variance of order $1/(nh)$; when the leading bias constant is nonzero, the MSE-optimal bandwidth is $h_{\mathrm{MSE}} = C_{\mathrm{MSE}}\, n^{-1/(2p+3)}$, with a constant that depends on the higher-order derivatives on the two sides, the conditional variances, the density, and the kernel. For local linear regression the leading bias is $O(h^2)$ and $h_{\mathrm{MSE}} \propto n^{-1/5}$.

(ii) Use a local quadratic fit on an auxiliary bandwidth $\tilde h$ to estimate the leading bias $\widehat{\mathrm{Bias}}$. If $S \geq 3$, $n\min\{h^5, \tilde h^5\}\max\{h^2, \tilde h^2\} \to 0$, and $n\min\{h, \tilde h\} \to \infty$, then

$$\frac{\hat\tau_{\mathrm{RD}} - \widehat{\mathrm{Bias}} - \tau_{\mathrm{RD}}}{\sqrt{V_{\text{bc}}}} \to_d \mathcal{N}(0, 1),$$

where $V_{\text{bc}}$ is the conditional variance of $\hat\tau_{\mathrm{RD}} - \widehat{\mathrm{Bias}}$, including the variance of the bias estimate itself and its covariance with $\hat\tau_{\mathrm{RD}}$; the limit of $h/\tilde h$ can be any value in $[0, \infty]$, and the MSE-optimal bandwidth satisfies these rates (the paper's Remark 10).
:::

The proof sketch has three steps (the paper's Appendix A and supplementary material). Step one (the smoothness in A12): the conditional bias of the local linear estimator is $h^2$ times a combination of the second derivatives on the two sides and observable constants, and the variance is of order $1/(nh)$; balancing the two gives (i). Step two (the rates in (i)): at the MSE-optimal bandwidth $nh^5$ converges to a positive constant, so the leading bias is of the same order as the standard error, the conventional $t$-statistic is not centered correctly, and the conventional interval undercovers, which is exactly the ailment (ii) repairs. Step three ($S \geq 3$): after the bias estimated by the local quadratic fit is subtracted, the remaining bias is of higher order; the bias estimate is random in finite samples, and putting its variance and its covariance with the original estimator into the standardization removes the need for $h/\tilde h \to 0$, so the MSE-optimal bandwidth can be used directly for inference. This yields the Robust line in the `rdrobust` output:

$$\text{CI}_{\text{robust}} = \big(\hat\tau_{\mathrm{RD}} - \widehat{\mathrm{Bias}}\big) \pm z_{1-\alpha/2}\sqrt{\widehat V_{\text{bc}}}, \qquad \widehat V_{\text{bc}} = \widehat V(\hat\tau) + \widehat V(\widehat{\mathrm{Bias}}) - 2\widehat{\mathrm{Cov}}(\hat\tau, \widehat{\mathrm{Bias}}).$$

In words: the robust interval still covers the effect at the cutoff, but it is centered on the bias-corrected estimate, so the conventional point estimate occasionally falls outside the robust interval; this is not a bug but by design, and the fuzzy estimate of Section 6.2 is an example. In Aster's block R1 the conventional interval is [0.133, 0.267], the robust interval is [0.111, 0.261], and the robust p-value is $1.26 \times 10^{-6}$. Two related choices also come from this theory. One is the coverage-error-optimal bandwidth (Calonico, Cattaneo and Farrell 2020), which is smaller than the MSE-optimal bandwidth: on R1 it is 1.173, about 0.625 of the latter, with estimate 0.186 and robust interval [0.095, 0.266]. The other is to set $\tilde h = h$, in which case the bias-corrected estimate is numerically equal to the local quadratic estimate at the local linear MSE-optimal bandwidth (Remark 7 of Calonico, Cattaneo and Titiunik 2014), which provides grounds for reporting local quadratic estimates directly. On R1 this chapter runs a different variant: local quadratic at its own MSE-optimal bandwidth of 4.185, with estimate 0.170 and robust interval [0.108, 0.253]. A12 requires a continuous score, and block R2 of Section 3.4 is exactly where it fails; there the robust interval's coverage is only 0.832. Honest inference is another route: the abstract and §1 of the arXiv version v6 of Armstrong and Kolesár (2020) point out that over the function class with a bounded second derivative, if the local linear estimator uses the bandwidth that minimizes worst-case MSE, the 95% interval only needs its critical value changed from 1.96 to 2.18.

### 4.3 Do Not Use High-Order Global Polynomials

Some researchers, to "control the curvature of $m$," stuff a fourth-, sixth-, or even higher-order global polynomial into the regression. Gelman and Imbens (2019) give three reasons against it: a high-order global polynomial assigns extreme and erratic weights to observations far from the cutoff, so the estimate at the cutoff is dragged around by distant data; the point estimate is highly sensitive to the order, and there is no principled way to choose the order; and the confidence intervals have poor coverage and are often too narrow. The Runge phenomenon of numerical analysis is an independent supplement (not one of Gelman and Imbens's three): high-degree polynomials produce spurious wiggles near the boundary of the support and can conjure a fake jump at the cutoff out of nothing, the extreme version of "mistaking nonlinearity for a discontinuity." Section 6.2 demonstrates this with numbers: on the same data, fourth-, sixth-, and eighth-order polynomials give wildly different jump estimates, one of them even negative. The right road is local linear or local quadratic with a bandwidth.

### 4.4 2SLS for Fuzzy RD

In estimation, fuzzy RD is 2SLS in a window around the cutoff, using $Z_i$ as an instrument for $D_i$, and both stages include the recentered running variable, with slopes that may differ across the two sides; the two-equation structure of Chapter 18 carries over directly:

$$\text{first stage: } D_i = \pi_0 + \pi_1 Z_i + m^D(X_i - c) + v_i, \qquad \text{structural: } Y_i = \beta_0 + \beta_1 D_i + m^Y(X_i - c) + \varepsilon_i,$$

where $m^D, m^Y$ are flexible polynomials specified separately on each side of the cutoff. The coefficient $\pi_1$ is the jump in treatment, and $\hat\beta_1$ is the fuzzy-RD LATE at the cutoff; when both stages use the same polynomial of order $p$, the same bandwidth, and equal weights within the window (a uniform kernel), it is numerically identical to the ratio formed from the four local intercepts of order $p$, with $p = 1$ being local linear (Wooldridge 2010 equations (21.107) and (21.108), citing Hahn, Todd and Van der Klaauw 2001). `rdrobust` runs it in one line with the `fuzzy=` argument, internally handling the local linear fit, the CCT bandwidth, and bias-corrected robust inference.

### 4.5 Covariates: For Precision Only

Predetermined covariates $\mathbf{W}_i$ play the same role in RD as in a randomized experiment: not identification, only precision. In Aster's block R1c, log GMV in the year before the badge launched and tenure on Aster explain about two-fifths of the residual variance of the outcome near 80 (a partial $R^2$ of 0.401); if this noise can be removed, the interval should shrink. The additive adjustment estimator recommended by Calonico, Cattaneo, Farrell and Titiunik (2019) adds the covariates to the local linear regression with a common coefficient on both sides:

$$Y_i = m(X_i) + \tau D_i + \mathbf{W}_i'\boldsymbol\gamma + \varepsilon_i.$$

The bold $\boldsymbol\gamma$ is the coefficient vector on the covariates, unrelated to the counterfactual trend $\gamma(\cdot)$ of Chapter 19.

::: {.assumption}
**A13 (smoothness of the covariates; Calonico, Cattaneo, Farrell and Titiunik 2019 Assumption SRD)** A12 holds for the outcome with $S = 3$; on each side of the cutoff, the conditional mean of $\mathbf{W}_i$ is three times continuously differentiable, the conditional mean of $\mathbf{W}_i Y_i$ is continuously differentiable, the conditional variance of $(Y_i, \mathbf{W}_i')'$ is continuously differentiable and invertible, and the fourth moments are continuous.
:::

A13 requires only that the conditional means of the covariates be smooth on each side separately, not that they be continuous at the cutoff; whether they jump is handled by Theorem 10 (i).

::: {.theorem}
**Theorem 10 (sharp RD with covariates; Calonico, Cattaneo, Farrell and Titiunik 2019 Lemma 1 and Theorem 2)** Under A13, if $h \to 0$ and $nh \to \infty$:

(i) The additive adjustment estimator satisfies $\tilde\tau \to_p \tau_{\mathrm{RD}} - \boldsymbol\tau_{\mathbf{W}}'\boldsymbol\gamma$, where $\boldsymbol\tau_{\mathbf{W}}$ is the jump at the cutoff in the conditional mean of $\mathbf{W}_i$ and $\boldsymbol\gamma$ is the coefficient of the local linear projection of the outcome on the covariates at the cutoff; $\boldsymbol\tau_{\mathbf{W}} = \mathbf{0}$ is sufficient for $\tilde\tau$ to be consistent, and the interacted version with a separate coefficient on each side requires stronger conditions.

(ii) If $\boldsymbol\tau_{\mathbf{W}} = \mathbf{0}$, $nh^7 \to 0$, and $\lim h/\tilde h < \infty$, the bias-corrected $t$-statistic is asymptotically standard normal and the variance estimator is consistent; both the MSE-optimal bandwidth and the bias correction must be recomputed for the adjusted estimator.
:::

The proof sketch has two steps. Step one (the algebra of the within-bandwidth covariate projection): $\tilde\tau = \hat\tau_{\mathrm{RD}} - \hat{\boldsymbol\tau}_{\mathbf{W}}'\tilde{\boldsymbol\gamma}$, that is, the RD estimate without covariates minus the RD estimate for the covariates themselves times the projection coefficients; the two terms converge to $\tau_{\mathrm{RD}}$ and $\boldsymbol\tau_{\mathbf{W}}'\boldsymbol\gamma$, respectively, which gives (i). Step two ($\boldsymbol\tau_{\mathbf{W}} = \mathbf{0}$ and A13): the adjustment term does not change the estimand, but it changes the constants of both the leading bias and the variance (the MSE expansion in the paper's Theorem 1 differs from the one without covariates), so redoing the MSE expansion and the bias correction along the lines of Theorem 9 gives (ii).

In words: covariate adjustment still estimates the effect at the cutoff, provided the covariates themselves do not jump at the cutoff, which is one more use of the balance tests of Section 3.2; what it buys is shorter intervals. Aster's block R1c is a fresh sample from the same world, with the two predetermined covariates above. Without covariates the robust interval has length 0.155; with them, 0.100, a reduction of 35.8% and a length ratio of 0.642. Part of this comes from the wider bandwidth after recomputation (3.20 against 1.89); compared at the same bandwidth, the length ratio is 0.796. The point estimate moves from 0.235 to 0.242, a difference of only 0.0069. If the covariates jump at the cutoff, $\tilde\tau$ estimates $\tau_{\mathrm{RD}} - \boldsymbol\tau_{\mathbf{W}}'\boldsymbol\gamma$, with a bias jointly determined by the jump in the covariates and the projection coefficients.

### 4.6 Heaping and Donut RD

Heaping of the running variable at certain integer values is another kind of trouble. Barreca, Lindo and Waddell (2016) point out that units at the heaping points are often not comparable to neighboring units, and local estimates are biased by them. A common remedy is the donut RD: drop the observations in the band $[c - w_{\mathrm{d}}, c + w_{\mathrm{d}}]$ right next to the cutoff, and use only the slightly outer points to extrapolate to $c$. There are two costs: the remaining points are farther from the cutoff, so the extrapolation bias grows; and the sample is smaller, so the interval widens. The working paper of Noack and Rothe (2026) (arXiv version v2) analyzes exactly this bias and variance cost of donut designs.

In Aster's block R1h, the sellers who push their reported score to 80.0 are weaker sellers with true scores between 76 and 80 (average true score 78.65); they get the badge, but their sales follow their true scores, so a batch of weak sellers crowds in just to the right of the cutoff. Direct estimation gives -0.034 with 95% interval [-0.129, 0.040]; even the sign is wrong. A donut that drops the 768 sellers with $|X_i - 80| \leq 0.5$ ($w_{\mathrm{d}} = 0.5$) gives 0.182, with interval [-0.110, 0.384] covering the truth of 0.20 but nearly three times wider; and dropping only the 460 sellers at 80.0 gives 0.273. There is no data-driven rule for the donut radius; the transparent practice is to report results for several radii.

![Heaping: the left panel is a histogram of reported scores, with 460 sellers at 80.0 and the region dropped by the donut in orange; the right panel shows robust 95% intervals for three estimates (naive, donut, and dropping only 80.0), with the dashed line at the truth of 0.20.](assets/fig/fig_21_heaping.svg)

### 4.7 Power

RD uses only observations near the cutoff, so the effective sample is far smaller than the total sample, and power should be computed in advance. `rdpower` (Cattaneo, Titiunik and Vazquez-Bare 2019) uses the bandwidth and robust variance of `rdrobust` (in practice taken from pilot or earlier-period data) and reports the power of the robust test and the minimum detectable effect (MDE). In Aster's block R1 the effective sample sizes are 695 and 689; at a true effect of 0.2 the power of the robust test is 0.9995, and when the effect is halved to 0.1 it is only 0.743; the MDE at 80% power is 0.107.

To sum up this section: RDD estimation moves from naive global regression to local linear regression, using local lines to remove boundary bias, data-driven bandwidths to balance bias and variance, and bias correction plus the extra variance to deliver robust intervals, while firmly avoiding high-order global polynomials; fuzzy RD is 2SLS at the cutoff, and you must recenter; covariates serve only precision, provided they do not jump at the cutoff; and heaping is sidestepped with a donut, at the cost of extrapolating farther.

## 5 Anchoring Papers

The three anchoring papers correspond to three threads of this chapter: the methodological history and classic vehicle of RDD, discontinuities manufactured by the rounding of platform ratings, and inference with a discrete running variable. Each is examined under five headings (paper, method, data, results, limitations), with the focus on how the assumptions are defended.

### 5.1 Thistlethwaite and Campbell (1960) and Lee (2008)

::: {.case}
The starting point of the method's history: Thistlethwaite and Campbell (1960) first proposed the discontinuity idea, using a score cutoff for a scholarship award to study the effect of the award on students' later performance. The classic vehicle that turned it into a modern workhorse of causal identification is Lee (2008).

The paper: "Randomized Experiments from Non-random Selection in U.S. House Elections," Journal of Econometrics. The question is incumbency advantage: does a party that wins a House seat this time have a higher probability of winning it again next time?

Method: the difficulty is that "good candidates get more votes anyway," so directly comparing winners and losers mixes in candidate quality. Lee's running variable is the Democratic margin in the previous election, that is, the party's vote share minus that of its strongest opponent, with cutoff 0, which in a two-party race is a 50% vote share. Districts narrowly won by just clearing 50% and districts narrowly lost by just missing 50% are nearly randomly comparable in voter composition and candidate quality, and which side of the line a district lands on depends on the chance of the last few votes. Sharp RD: a margin above 0 delivers the seat.

Data: district-level data from U.S. House elections, vote shares and next-election outcomes.

Results: in districts the party narrowly won, its probability of winning and its vote share in the next election are significantly higher, so incumbency advantage is a real causal effect, not merely a reflection of candidate quality. Lee's diagnostics: predetermined covariates (such as results from the preceding election) do not jump at the threshold; a continuous conditional density of each district's vote is his assumption, supported by the element of chance in vote counting.

Limitations: the estimate is local, identifying only the effect in closely contested districts with vote shares near 50%, and it does not extrapolate to safe districts dominated by one party.
:::

This paper's defense strategy is the RDD template: not piling control variables into a regression, but arguing the claim itself that "assignment at the threshold is nearly random," with the balance of predetermined covariates as corroborating evidence; the density test is a further shadow, supplied separately by McCrary (2008).

### 5.2 Anderson and Magruder (2012)

::: {.case}
The paper: "Learning from the Crowd: Regression Discontinuity Estimates of the Effects of an Online Review Database," The Economic Journal. It shares a core with Aster's rating filter: the platform rounds a continuous rating to discrete tiers, manufacturing an artificial discontinuity.

Method: Yelp rounds a restaurant's average rating to the nearest half star. Take two restaurants whose true averages fall on either side of the midpoint between two half-star tiers: the one slightly below the midpoint displays the lower half star, and the one slightly above it displays the higher one. Their underlying quality is almost the same, yet their displayed ratings differ by half a star. The authors take the rounding threshold as the cutoff and compare restaurants just rounded up to the higher half star with those just rounded down to the lower one.

Data: all Yelp reviews of San Francisco restaurants, used to reconstruct each restaurant's true average rating, matched to real-time reservation availability recorded by a large online reservation website; the outcome is whether a restaurant's online reservations are sold out.

Results: an extra half star makes a restaurant 19 percentage points (49%) more likely to sell out its reservations; robustness checks find no evidence that ratings are manipulated.

Limitations: the sample consists only of San Francisco restaurants listed on the reservation website, and if listing depends on customer traffic, the estimate should be read as a lower bound; only Yelp users are affected; and the estimate is local to the rounding thresholds. Rounded ratings often make the score discrete (Section 3.4), but these restaurants have many reviews, so their average ratings are nearly continuous.
:::

Platform systems are full of rounding, thresholds, and tiering rules, and each one is a natural discontinuity. Luca (2011), a Harvard Business School working paper (No. 12-016), uses the same Yelp rounding rule to study restaurant revenue. Information systems and marketing journals contain similar designs: Sahoo, Dellarocas and Srinivasan (2018, Information Systems Research) use the discontinuity in displayed ratings to study product returns, Zhong (2022, Management Science) studies the jump in sales when Taobao sellers cross rating thresholds, and Narayanan and Kalyanam (2015, Marketing Science) use RD to estimate the causal effect of position in search advertising, where positions are determined by auction.

### 5.3 Kolesár and Rothe (2018) Revisit Oreopoulos (2006)

::: {.case}
The paper: "Inference in Regression Discontinuity Designs with a Discrete Running Variable," American Economic Review (hereafter KR). The numbers and numbering in this section are taken from §6.1 and Table 2 of its arXiv version v4. The application revisits the study by Oreopoulos (2006, American Economic Review) of the British minimum school-leaving age.

Method: in 1947 Great Britain raised the minimum school-leaving age from 14 to 15. Oreopoulos uses a sharp RD in which the running variable is the year a worker turned 14, the threshold is 1947, the outcome is log annual earnings, and the specification is a global fourth-order polynomial. KR compare EHW standard errors, CRV standard errors clustered by the running variable, and two kinds of honest intervals (BME, for bounded specification error, and BSD, for the bounded second derivative of Section 3.4).

Data: workers in Great Britain who turned 14 between 1935 and 1965, from the UK General Household Survey. KR use the data accompanying Oreopoulos's later online correction, so the estimates for the original specification differ from those in the 2006 paper. The running variable has only 31 support points, 19 above the threshold and 12 below.

Results: the estimate for the original specification (Table 2, column 1) is 0.055, and the EHW standard error of 0.030 is twice the CRV standard error of 0.015; the EHW 95% interval [-0.003, 0.113] contains zero, while the CRV interval [0.026, 0.084] does not. Except in columns 1 and 2 (column 2 is a linear specification on the full sample, where the EHW standard error is slightly smaller than the CRV standard error), the ratio of EHW to CRV standard errors in the other specifications ranges from 1.6 to 28.9, and it is larger the narrower the bandwidth and the higher the polynomial order. The BSD intervals take four values of the smoothness constant (the $M$ of Section 3.4), 0.004, 0.02, 0.04, and 0.2, of which the authors consider 0.02 and 0.04 reasonable; the data-based lower bound on the smoothness constant has a point estimate of 0.012 and a 95% confidence interval of [0.0013, ∞).

Limitations: the authors explain that this analysis does not imply that the EHW interval has correct coverage, but any interval with correct coverage must be at least as wide as the EHW interval, so the significance conclusions based on CRV do not stand. The honest intervals depend on the researcher's judgment about the smoothness constant, which the authors calibrate with the intuition of "how far the function can deviate from a straight line within one year," reporting a range of values.
:::

This paper puts the two conclusions of Section 3.4 to work on real data: with few support points, clustered standard errors are systematically too small, and significance conclusions tip to the wrong side along with them; honest intervals offer no free precision but put the judgment of "how curved the function is" on the table.

Taken together: Lee (2008) turns near-randomness at the threshold into a defensible identification argument, Anderson and Magruder (2012) show that a platform's rounding rules by themselves manufacture discontinuities, and Kolesár and Rothe (2018) warn that the scores behind such discontinuities are often discrete, so inference must switch to honest intervals.

## 6 A Full Walkthrough on Aster Data

Now run the earlier tools end to end on Aster's seven data blocks. The code below uses R 4.6.1, and every number quoted in the text comes from the actual output of this code. RDD uses `rdrobust` and `rddensity`, local randomization `rdlocrand` (Cattaneo, Titiunik and Vazquez-Bare 2016), honest inference `RDHonest`, multiple cutoffs and geographic RD `rdmulti` (Cattaneo, Titiunik and Vazquez-Bare 2020), and power `rdpower`. The listings below show the main calls; auxiliary computations of a few lines, such as the extrapolation and the Lemma 1 weighted estimator, are not listed line by line.

### 6.1 The Data-Generating Process for RDD

The design parameters are as follows: fix `set.seed(101)`, with $N = 12000$ sellers. The running variable is a continuous quality score drawn from a normal distribution, with scores below 45 and above 100 pushed onto the endpoints (this step is what the code comment calls smooth truncation; each end gets a small lump of point mass, far from 80), and the density is continuous at 80 (that is, there is no manipulation). $Y(0)$ is a smooth nonlinear function of the score that rises steeply with it (high-scoring sellers simply sell more), with a sinusoidal wiggle that is smooth at the cutoff yet forces a global polynomial to bend across the whole range. The true sharp jump is 0.20.

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

The naive "badged vs. unbadged" difference is 1.891, about 9.5 times the truth of 0.200, confirming the number of Section 1. The seven data blocks are as follows. R1 is the inherited block above, with its random number stream preserved as is; each of the six new blocks uses its own seed and consists of column vectors of equal length plus a `meta` list, where `x` is the score, `y` is log GMV, and `d` is treatment; the calls in Section 6.3 read the fields in the last column directly.

| Block | Story | Sellers | Seed | Truth | Fields |
|---|---|---|---|---|---|
| R1 | The badge world above, in sharp and fuzzy variants | 12,000 | 101 | Jump 0.20 | `score`, `y_sharp`, `D_fuzzy`, `y_fuzzy` |
| R1c | The R1 world plus two predetermined covariates: log GMV in the year before launch and tenure on Aster | 12,000 | 2101 | Jump 0.20 | `x`, `d`, `y`, `gmv_prev`, `tenure` |
| R1h | R1 a year later: scores reported to one decimal place, and some sellers push their reported score to 80.0 | 12,000 | 2102 | Jump 0.20 | `x` (reported score), `x_true`, `d`, `y` |
| R2 | Average of the 10 most recent reviews; from 4.25 on, a seller enters the 4.5-star filter | 8,000 | 2103 | Jump 0.15; $M = 7.79$ | `x`, `d`, `y`, `tenure`, `meta$M_bound` |
| R3 | Commission rate of 12%, 0.6 percentage points lower for every \$100 above \$1,500 | 60,000 | 2104 | -0.05 per percentage point | `x` (previous-quarter GMV), `b` (commission rate), `y` |
| R4 | Thresholds of 70, 80, and 90 for vintage, home, and electronics | 24,000 | 2105 | 0.40, 0.25, 0.10; 0.325 at 75 | `x`, `cutoff`, `d`, `y` |
| R4g | Same-day delivery within 25 km of a fulfillment center | 40,000 | 2106 | 0.10 to 0.30 along the boundary | `x1`, `x2` (km), `d`, `y`, `meta$boundary` |

### 6.2 RDD Estimation: Local Linear, Bandwidth, Polynomials, Density, Fuzzy

Sharp RD uses the local linear fit of `rdrobust` with the CCT bandwidth:

```r
library(rdrobust); library(rddensity)
rd <- rdrobust(y = y_sharp, x = score, c = 80,
               kernel = "triangular", p = 1, bwselect = "mserd")
```

The conventional point estimate is 0.200, with conventional 95% interval [0.133, 0.267]; the bias-corrected robust 95% interval is [0.111, 0.261], with p-value $1.26 \times 10^{-6}$. The CCT/MSE-optimal bandwidth is 1.876, meaning that only sellers between 78.124 and 81.876 are used, with effective sample sizes of 695 and 689 on the two sides. Local linear lands cleanly on the truth of 0.200. One reminder: in this DGP the effect is 0.20 for every seller, so the effect at the cutoff happens to equal the ATE for everyone, which is only a convenience for reconciliation; in the real world the effect varies with the score, and local linear still recovers only the effect at the cutoff. Block R4 is such a world.

Bandwidth sensitivity is a mandatory robustness check. Re-estimate at half and double the bandwidth:

```r
rd_half <- rdrobust(y = y_sharp, x = score, c = 80, kernel = "triangular", p = 1, h = 1.876 / 2)
rd_dbl  <- rdrobust(y = y_sharp, x = score, c = 80, kernel = "triangular", p = 1, h = 1.876 * 2)
```

The half bandwidth gives 0.186 and the double bandwidth 0.237, both floating modestly around the truth, so the estimate is insensitive to the bandwidth. Next, a demonstration of why high-order global polynomials should not be used. Rescale the running variable to $[-1, 1]$ and use raw powers (so that all powers are zero at the cutoff and the coefficient on $D$ is exactly the jump), fitting fourth-, sixth-, and eighth-order polynomials:

```r
sc_c <- (score - 80) / max(abs(score - 80))        # rescale to [-1,1]
poly_jump <- function(k) coef(lm(y_sharp ~ D_sharp * poly(sc_c, k, raw = TRUE)))["D_sharp"]
poly_jump(4); poly_jump(6); poly_jump(8)
```

The three estimates are 0.556, -0.215, and 0.254, swinging violently, with the sixth-order polynomial giving a meaningless negative number: the warning of Gelman and Imbens comes true here. The manipulation test uses `rddensity`:

```r
rddensity(X = score, c = 80)                        # manipulation test
```

The p-value for the density jump is 0.250 (T = 1.151), far from significant, so no-manipulation holds up. Last comes fuzzy RD, instrumenting the badge as an imperfect treatment:

```r
rd_fz <- rdrobust(y = y_fuzzy, x = score, c = 80, fuzzy = D_fuzzy,
                  kernel = "triangular", p = 1, bwselect = "mserd")
```

The first-stage jump is estimated at 0.749 (truth 0.750): crossing the threshold raises the probability of getting the badge by about three-quarters. The reduced-form jump is 0.170, and the hand-computed Wald ratio is $0.170 / 0.749 = 0.227$, with numerator and denominator each using its own bandwidth (Section 3.5). The fuzzy-RD LATE point estimate from `rdrobust` is 0.241, with robust 95% interval [0.071, 0.236]. The fuzzy estimate (0.241) is a bit higher than the sharp one (0.200), partly because of the wider bandwidth (3.333 against 1.876); the rest is sampling variation from the separate draw of noise for the fuzzy outcome. The conventional point estimate of 0.241 falls just outside the robust interval's upper bound of 0.236, which is exactly how Section 4.2 said a robust interval centered on the bias-corrected estimate behaves. When reading the conclusion, rely on the robust interval; both the hand-computed ratio of 0.227 and the 0.241 from `rdrobust` land near the true structural effect of 0.20.

### 6.3 Calls on the New Data Blocks

```r
library(RDHonest); library(rdlocrand); library(rdmulti); library(rdpower)
# R2: rating average with 21 support points, cutoff 4.25
rd <- rdrobust(R2$y, R2$x, c = 4.25, masspoints = "adjust")
cr <- rdrobust(R2$y, R2$x, c = 4.25, masspoints = "adjust", h = rd$bws["h", "left"],
               b = rd$bws["b", "left"], cluster = R2$x, vce = "hc1")   # clustered by support point
ho <- RDHonest(y ~ x, data = data.frame(y = R2$y, x = R2$x), cutoff = 4.25, M = R2$meta$M_bound)
RDSmoothnessBound(ho, s = 2)              # lower bound on M: point estimate and one-sided CI
ws <- rdwinselect(R2$x, R2$tenure, cutoff = 4.25, wmasspoints = TRUE, reps = 1000, seed = 2103)
rdrandinf(R2$y, R2$x, cutoff = 4.25, wl = ws$w_left, wr = ws$w_right,
          ci = c(0.05, seq(0, 0.3, by = 0.0025)), reps = 1000, seed = 2103)
```

```r
# R3: kink at USD 1,500; the schedule's slope changes by -0.006 percentage points per USD
rk <- rdrobust(R3$y, R3$x, c = 1500, deriv = 1, p = 2)
rk$coef["Conventional", 1] / -0.006       # kink effect per percentage point of commission
rdrobust(R3$y, R3$x, c = 1500)            # level jump at the kink: should be zero
# R4: three cutoffs; R4g: six points on the boundary of the 25 km zone
rdmc(Y = R4$y, X = R4$x, C = R4$cutoff)   # per-cutoff and normalized pooled estimates
bp <- R4g$meta$boundary; r <- sqrt(R4g$x1^2 + R4g$x2^2)
rdms(Y = R4g$y, X = R4g$x1, C = bp$b1, X2 = R4g$x2, zvar = R4g$d, C2 = bp$b2, xnorm = 25 - r)
```

```r
# R1c: covariates; R1h: heaping and a donut of radius 0.5; R1: power
rdrobust(R1c$y, R1c$x, c = 80, covs = cbind(R1c$gmv_prev, R1c$tenure))
rddensity(X = R1h$x, c = 80)
keep <- abs(R1h$x - 80) > 0.5
rdrobust(R1h$y[keep], R1h$x[keep], c = 80)
rdpower(data = cbind(y_sharp, score), cutoff = 80, tau = 0.2)
rdmde(data = cbind(y_sharp, score), cutoff = 80, init.cond = 0.1)
```

The output of these calls has already been interpreted, one by one, in Sections 3 and 4. Two implementation notes. First, R2's average rating takes only 21 values, so `masspoints = "check"` and `"off"` throw errors because there are too few support points; the `rdrobust` default `masspoints = "adjust"` runs, but it only lets bandwidth selection operate on repeated values and does not repair coverage: in the repeated simulations the conventional interval's coverage is still 0.784. Second, `rdmde` may fail to converge from its default starting value, so a verified starting value should be fixed. The Monte Carlo for the honest intervals uses 1000 samples redrawn from the same DGP, recomputing all intervals for each.

### 6.4 The Grand Reconciliation

| Block | Estimator | Estimate | 95% interval | Truth | Covered |
|---|---|---|---|---|---|
| R1 | Naive comparison | 1.891 | | 0.200 | |
| R1 | Local linear | 0.200 | [0.111, 0.261] | 0.200 | Yes |
| R1 | Fuzzy RD | 0.241 | [0.071, 0.236] | 0.20 | Yes |
| R1c | With covariates | 0.242 | [0.186, 0.286] | 0.20 | Yes |
| R1h | Naive (heaping) | -0.034 | [-0.129, 0.040] | 0.20 | No |
| R1h | Donut, $w_{\mathrm{d}} = 0.5$ | 0.182 | [-0.110, 0.384] | 0.20 | Yes |
| R2 | Conventional | 0.092 | [0.041, 0.143] | 0.15 | No |
| R2 | Clustered by support point | 0.092 | [0.079, 0.106] | 0.15 | No |
| R2 | Honest, $M = 7.79$ | 0.085 | [-0.031, 0.202] | 0.15 | Yes |
| R2 | Local randomization | 0.1295 | [0.0875, 0.1725] | 0.15 | Yes |
| R3 | Kink, $p = 2$ | -0.052 | [-0.083, -0.018] | -0.050 | Yes |
| R4 | Lemma 1 weighted | 0.322 | [0.294, 0.361] | 0.316 | Yes |
| R4 | Normalized pooled | 0.293 | [0.223, 0.369] | 0.316 | Yes |
| R4 | Extrapolated to 75 | 0.349 | [0.314, 0.384] | 0.325 | Yes |
| R4g | Boundary pooled | 0.204 | [0.171, 0.230] | 0.20 | Yes |

Interval types depend on the estimator: the first two R2 rows use conventional intervals, the honest interval comes from `RDHonest`, local randomization uses the test-inversion interval, the extrapolation uses a conventional interval with a fixed bandwidth, the interval for the Lemma 1 weighted estimator is assembled from the bias-corrected estimates and robust variances at each cutoff with the weights treated as fixed, and the rest are robust bias-corrected intervals. The truth of 0.316 in the first two R4 rows is the target weighted by the sample weights; the population target of Theorem 6 is 0.317. The Covered column describes only this one sample; for the coverage of the methods, see the repeated simulations of Section 3.4.

To sum up the reconciliation: the naive comparison of 1.891 is inflated by selection bias to about 9.5 times the truth of 0.200; local linear cleanly recovers 0.200 with a robust interval that covers the truth; the estimate is robust to the bandwidth, and the density shows no jump; and the fuzzy variant, instrumenting the badge, puts the LATE near 0.20. The new data blocks show that with a discrete score both the conventional and the clustered intervals miss the truth, that heaping can push the estimate to the wrong sign, and that the pooled number across multiple cutoffs is a weighted average, not a simple one.

## 7 Failure Modes and Robustness

In a simulation the identifying assumptions are built in, but in real research they can fail at any time. This section lays out the most common ways RDD fails and the actionable responses.

The number-one threat is manipulation and sorting. Once units can precisely control which side of the cutoff they land on, the two sides of the threshold are no longer comparable. Diagnose it by using `rddensity` to check whether the density bunches, and at the same time check whether predetermined covariates pass smoothly through the threshold; R1h shows that once scores can be rounded up, the density test sounds a loud alarm and the naive estimate can get even the sign wrong. When the heaping sits right at the cutoff, the donut RD is a common remedy; report a set of radii, not only the one that looks best.

The second class is bandwidth sensitivity and polynomial overfitting. Always report how the estimate changes with the bandwidth; do not use global polynomials of order three or higher, but local linear or local quadratic fits instead; and in kink designs use local quadratic by default.

The third class is a discrete running variable. When the score takes only a few values, clustering by score value is no remedy but makes things worse, and bias correction cannot rescue coverage either. Either report an honest interval with a sensitivity analysis over $M$, or take the local randomization route between the nearest support points on either side of the cutoff and report the results for wider windows as well. $M$ cannot be estimated from the data: `RDSmoothnessBound` gives a lower bound, and using it as $M$ makes the interval too narrow.

The fourth class is other policies jumping at the cutoff at the same time, that is, compound treatment. If 80 not only awards the badge but also unlocks a lower fee or more traffic slots, the jump at the threshold is their combined force. This is not a statistical problem but an institutional one, and it can be ruled out only by pinning down exactly which treatments the threshold triggers. It is especially common at geographic boundaries, since administrative boundaries often switch tax rates, school districts, and law enforcement at the same time.

The fifth class is aggregation and extrapolation. The pooled number across multiple cutoffs is an average weighted by densities and shares; state the weights when reporting it, and report cutoff by cutoff when needed. Extrapolating beyond the cutoff relies on untestable assumptions such as constant bias, which can only be falsified with the parallelism below the cutoff. A last, purely technical but common problem is plotting: an RDD plot is a plot of binned means (`rdplot` uses evenly spaced bins by default, and quantile-spaced bins are also available), and drawing one casually with covariates gets it wrong; the right approach is to treat the binned means as a semiparametric regression $Y_i = m(X_i) + \mathbf{W}_i'\boldsymbol\gamma + \varepsilon_i$ (`binsreg`; Cattaneo, Crump, Farrell and Feng 2025).

Stringing these failure modes together, the lesson is that RDD's credibility lies not in estimation technique but in the design and the institutions, and its success or failure hinges on "whether this cutoff triggers only one thing, whether units can manipulate it, and how many values the score takes." The density test, the bandwidth curve, covariate smoothness, honest intervals, and window sensitivity are all evidence offered around these questions, and none of them can substitute for a substantive judgment about the institutional background.

## 8 Further Reading

::: {.readings}
Required reading, in suggested order:

- Lee and Lemieux (2010, Journal of Economic Literature). The most readable RDD survey, covering everything from continuity and bandwidths to the diagnostic checklist; the first read for beginners.
- Imbens and Lemieux (2008, Journal of Econometrics). Another cornerstone of practical guidance; focus on the discussion of bandwidths and functional form. The RD subsection of Wooldridge (2010) §21.5 draws on it, and its equations (21.96) to (21.108) are checkable versions of Theorems 1 and 4 of this chapter.
- Calonico, Cattaneo and Titiunik (2014, Econometrica). Bias-corrected robust inference and the MSE-optimal bandwidth: the theory behind `rdrobust` and the source of Theorem 9 of this chapter.
- Cattaneo, Idrobo and Titiunik (2020, 2024, Cambridge booklets). The first volume is a practical manual for the continuity framework, and the second covers local randomization, fuzzy RD, discrete running variables, and multidimensional RD: the operating checklist to work through before writing an RDD paper.
- Kolesár and Rothe (2018, American Economic Review) and Armstrong and Kolesár (2020, Quantitative Economics). Discrete running variables and honest intervals. The numbering of Theorem 3 of this chapter is taken from the arXiv version of the former (`arXiv:1606.04086v4`), and the section numbers and numerical results cited from the latter come from its arXiv version (`arXiv:1606.01200v6`). Armstrong and Kolesár (2018, Econometrica) give the general theory of the optimality of such intervals.

Further reading:

- McCrary (2008, Journal of Econometrics) and Cattaneo, Jansson and Ma (2018, Stata Journal; 2020, Journal of the American Statistical Association). The original manipulation test, and the software documentation and theory of `rddensity`.
- Gelman and Imbens (2019, Journal of Business & Economic Statistics). Why high-order global polynomials should not be used; short and forceful.
- Imbens and Kalyanaraman (2012, Review of Economic Studies) and Calonico, Cattaneo and Farrell (2020, Econometrics Journal). MSE-optimal and coverage-error-optimal bandwidths.
- Calonico, Cattaneo, Farrell and Titiunik (2019, Review of Economics and Statistics). Covariates in RD.
- Cattaneo, Frandsen and Titiunik (2015, Journal of Causal Inference). The original paper on local randomization and randomization inference.
- Card, Lee, Pei and Weber (2015, Econometrica). Identification and inference in regression kink designs.
- Cattaneo, Keele, Titiunik and Vazquez-Bare (2016, Journal of Politics; 2021, Journal of the American Statistical Association) and Keele and Titiunik (2015, Political Analysis). Weights and extrapolation with multiple cutoffs, and geographic RD.
- Angrist and Rokkanen (2015, Journal of the American Statistical Association). Extrapolation away from the cutoff using covariates.
- Barreca, Lindo and Waddell (2016, Economic Inquiry) and Noack and Rothe (2026, working paper, arXiv version v2). Heaping and donut designs.
- Hahn, Todd and Van der Klaauw (2001, Econometrica). The original paper on sharp and fuzzy RD identification.
- Thistlethwaite and Campbell (1960, Journal of Educational Psychology). The historical origin of the RDD idea, worth a read to know where it came from.
- Lee (2008, Journal of Econometrics) and Anderson and Magruder (2012, Economic Journal). Two of the anchoring papers of Section 5, one in political science and one on platform ratings.
:::

::: {.apa-refs}
- Anderson, M., & Magruder, J. (2012). Learning from the crowd: Regression discontinuity estimates of the effects of an online review database. *The Economic Journal, 122*(563), 957-989. https://doi.org/10.1111/j.1468-0297.2012.02512.x
- Angrist, J. D., & Rokkanen, M. (2015). Wanna get away? Regression discontinuity estimation of exam school effects away from the cutoff. *Journal of the American Statistical Association, 110*(512), 1331-1344. https://doi.org/10.1080/01621459.2015.1012259
- Armstrong, T. B., & Kolesár, M. (2018). Optimal inference in a class of regression models. *Econometrica, 86*(2), 655-683. https://doi.org/10.3982/ECTA14434
- Armstrong, T. B., & Kolesár, M. (2020). Simple and honest confidence intervals in nonparametric regression. *Quantitative Economics, 11*(1), 1-39. https://doi.org/10.3982/QE1199
- Barreca, A. I., Lindo, J. M., & Waddell, G. R. (2016). Heaping-induced bias in regression-discontinuity designs. *Economic Inquiry, 54*(1), 268-293. https://doi.org/10.1111/ecin.12225
- Calonico, S., Cattaneo, M. D., & Farrell, M. H. (2020). Optimal bandwidth choice for robust bias-corrected inference in regression discontinuity designs. *The Econometrics Journal, 23*(2), 192-210. https://doi.org/10.1093/ectj/utz022
- Calonico, S., Cattaneo, M. D., Farrell, M. H., & Titiunik, R. (2019). Regression discontinuity designs using covariates. *The Review of Economics and Statistics, 101*(3), 442-451. https://doi.org/10.1162/rest_a_00760
- Calonico, S., Cattaneo, M. D., & Titiunik, R. (2014). Robust nonparametric confidence intervals for regression-discontinuity designs. *Econometrica, 82*(6), 2295-2326. https://doi.org/10.3982/ECTA11757
- Card, D., Lee, D. S., Pei, Z., & Weber, A. (2015). Inference on causal effects in a generalized regression kink design. *Econometrica, 83*(6), 2453-2483. https://doi.org/10.3982/ECTA11224
- Cattaneo, M. D., Crump, R. K., Farrell, M. H., & Feng, Y. (2025). Binscatter regressions. *The Stata Journal, 25*(1), 3-50. https://doi.org/10.1177/1536867X251322960
- Cattaneo, M. D., Frandsen, B. R., & Titiunik, R. (2015). Randomization inference in the regression discontinuity design: An application to party advantages in the U.S. Senate. *Journal of Causal Inference, 3*(1), 1-24. https://doi.org/10.1515/jci-2013-0010
- Cattaneo, M. D., Idrobo, N., & Titiunik, R. (2020). *A practical introduction to regression discontinuity designs: Foundations*. Cambridge University Press. https://doi.org/10.1017/9781108684606
- Cattaneo, M. D., Idrobo, N., & Titiunik, R. (2024). *A practical introduction to regression discontinuity designs: Extensions*. Cambridge University Press. https://doi.org/10.1017/9781009441896
- Cattaneo, M. D., Jansson, M., & Ma, X. (2018). Manipulation testing based on density discontinuity. *The Stata Journal, 18*(1), 234-261. https://doi.org/10.1177/1536867X1801800115
- Cattaneo, M. D., Jansson, M., & Ma, X. (2020). Simple local polynomial density estimators. *Journal of the American Statistical Association, 115*(531), 1449-1455. https://doi.org/10.1080/01621459.2019.1635480
- Cattaneo, M. D., Keele, L., Titiunik, R., & Vazquez-Bare, G. (2016). Interpreting regression discontinuity designs with multiple cutoffs. *The Journal of Politics, 78*(4), 1229-1248. https://doi.org/10.1086/686802
- Cattaneo, M. D., Keele, L., Titiunik, R., & Vazquez-Bare, G. (2021). Extrapolating treatment effects in multi-cutoff regression discontinuity designs. *Journal of the American Statistical Association, 116*(536), 1941-1952. https://doi.org/10.1080/01621459.2020.1751646
- Cattaneo, M. D., Titiunik, R., & Vazquez-Bare, G. (2016). Inference in regression discontinuity designs under local randomization. *The Stata Journal, 16*(2), 331-367. https://doi.org/10.1177/1536867X1601600205
- Cattaneo, M. D., Titiunik, R., & Vazquez-Bare, G. (2019). Power calculations for regression-discontinuity designs. *The Stata Journal, 19*(1), 210-245. https://doi.org/10.1177/1536867X19830919
- Cattaneo, M. D., Titiunik, R., & Vazquez-Bare, G. (2020). Analysis of regression-discontinuity designs with multiple cutoffs or multiple scores. *The Stata Journal, 20*(4), 866-891. https://doi.org/10.1177/1536867X20976320
- Gelman, A., & Imbens, G. (2019). Why high-order polynomials should not be used in regression discontinuity designs. *Journal of Business & Economic Statistics, 37*(3), 447-456. https://doi.org/10.1080/07350015.2017.1366909
- Hahn, J., Todd, P., & Van der Klaauw, W. (2001). Identification and estimation of treatment effects with a regression-discontinuity design. *Econometrica, 69*(1), 201-209. https://doi.org/10.1111/1468-0262.00183
- Imbens, G., & Kalyanaraman, K. (2012). Optimal bandwidth choice for the regression discontinuity estimator. *The Review of Economic Studies, 79*(3), 933-959. https://doi.org/10.1093/restud/rdr043
- Imbens, G. W., & Lemieux, T. (2008). Regression discontinuity designs: A guide to practice. *Journal of Econometrics, 142*(2), 615-635. https://doi.org/10.1016/j.jeconom.2007.05.001
- Keele, L. J., & Titiunik, R. (2015). Geographic boundaries as regression discontinuities. *Political Analysis, 23*(1), 127-155. https://doi.org/10.1093/pan/mpu014
- Kolesár, M. (2026). *RDHonest: Honest inference in regression discontinuity designs* (R package version 1.0.2) [Computer software]. https://doi.org/10.32614/CRAN.package.RDHonest
- Kolesár, M., & Rothe, C. (2018). Inference in regression discontinuity designs with a discrete running variable. *American Economic Review, 108*(8), 2277-2304. https://doi.org/10.1257/aer.20160945
- Lee, D. S. (2008). Randomized experiments from non-random selection in U.S. House elections. *Journal of Econometrics, 142*(2), 675-697. https://doi.org/10.1016/j.jeconom.2007.05.004
- Lee, D. S., & Card, D. (2008). Regression discontinuity inference with specification error. *Journal of Econometrics, 142*(2), 655-674. https://doi.org/10.1016/j.jeconom.2007.05.003
- Lee, D. S., & Lemieux, T. (2010). Regression discontinuity designs in economics. *Journal of Economic Literature, 48*(2), 281-355. https://doi.org/10.1257/jel.48.2.281
- Luca, M. (2011). *Reviews, reputation, and revenue: The case of Yelp.com* (Harvard Business School Working Paper No. 12-016). SSRN. https://doi.org/10.2139/ssrn.1928601
- McCrary, J. (2008). Manipulation of the running variable in the regression discontinuity design: A density test. *Journal of Econometrics, 142*(2), 698-714. https://doi.org/10.1016/j.jeconom.2007.05.005
- Narayanan, S., & Kalyanam, K. (2015). Position effects in search advertising and their moderators: A regression discontinuity approach. *Marketing Science, 34*(3), 388-407. https://doi.org/10.1287/mksc.2014.0893
- Noack, C., & Rothe, C. (2026). *Donut regression discontinuity designs* (arXiv:2308.14464v2). arXiv. https://doi.org/10.48550/arXiv.2308.14464
- Oreopoulos, P. (2006). Estimating average and local average treatment effects of education when compulsory schooling laws really matter. *American Economic Review, 96*(1), 152-175. https://doi.org/10.1257/000282806776157641
- Sahoo, N., Dellarocas, C., & Srinivasan, S. (2018). The impact of online product reviews on product returns. *Information Systems Research, 29*(3), 723-738. https://doi.org/10.1287/isre.2017.0736
- Thistlethwaite, D. L., & Campbell, D. T. (1960). Regression-discontinuity analysis: An alternative to the ex post facto experiment. *Journal of Educational Psychology, 51*(6), 309-317. https://doi.org/10.1037/h0044319
- Wooldridge, J. M. (2010). *Econometric analysis of cross section and panel data* (2nd ed.). MIT Press.
- Zhong, Z. (2022). Chasing diamonds and crowns: Consumer limited attention and seller response. *Management Science, 68*(6), 4380-4397. https://doi.org/10.1287/mnsc.2021.4060
:::
