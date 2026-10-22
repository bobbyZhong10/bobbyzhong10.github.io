---
title: "DiD II: Fuzzy, Continuous & Switching Treatments"
subtitle: "Beyond the Clean 0-to-1 Switch"
seriesline: "Foundations of Information Systems Economics · Chapter 20"
series: "Foundations of Information Systems Economics: Empirical Methods for Information Systems Economics"
sidebartitle: "Foundations of Information Systems Economics"
sidebarsub: "Theory · Methods · Applications"
chapterlabel: "Chapter 20 · DiD II: Fuzzy, Continuous & Switching Treatments"
---

## Introduction

Tessera is a fictional food-delivery platform. After it launches a courier incentive program, management wants to know how many more orders the program got couriers to deliver. Kestrel in Chapter 19 had a tidy structure: a city either cut its commission or did not, a cut was never reversed, and treatment jumped from 0 to 1 once and for all at the city level. Tessera's treatment is not so tidy. The program opens in one batch of cities and couriers enroll voluntarily: in the launch cities some couriers enroll and others do not, and in the cities without the launch some couriers joined an older program long ago. In other cities the per-order bonus differs from city to city. In still others the bonus tier is adjusted week by week, switched on and paused, paused and switched back on. Treatment has changed from "whether" to "how much, for whom, and when."

These three forms are more common in empirical work than the tidy switch from 0 to 1. de Chaisemartin and D'Haultfœuille (hereafter dCDH) found that 10.1% of the papers published in the American Economic Review from 2010 to 2012 use a fuzzy DiD design, estimating a simple Wald-DID or a weighted average of several Wald-DIDs (dCDH 2018). They later surveyed the 100 most cited papers in the same journal from 2015 to 2019: 26 of them run TWFE regressions, and in only 4 is the treatment binary and never switched off once it is on (dCDH 2026). The core question of this chapter is therefore direct: when treatment does not move a whole group of units from 0 to 1 at once, but instead involves partial take-up, varying intensity, or switching on and off, what does DiD identify?

This chapter follows the approach of Chapter 19. Section 1 opens with a wrong number from Tessera. Section 2 writes down the causal quantity that each of the three treatment forms targets. Section 3 is the most detailed section of the chapter and gives the identification results one by one: the three Wald ratios of fuzzy DiD and partial identification when the control group's treatment rate moves; $\mathrm{ATT}(d \mid d)$, ACRT, and strong parallel trends under a continuous dose; two designs, with and without stayers; $DID_M$ and $DID_\ell$ under switching treatments; and finally a bridge to the IV of Chapter 18. Section 4 covers estimation, and the negative weights of TWFE dose regressions are gathered there too. Section 5 examines three anchoring papers, Section 6 runs the whole toolkit on Tessera, and Section 7 discusses failure modes. The chapter uses the potential outcomes of Chapter 15, the LATE theorem of Chapter 18, and the parallel trends and $DID_\ell$ of Chapter 19. The question running through the chapter is still the one from Chapter 19: who plays the counterfactual on behalf of the treated group? Now there is a second question: which part of the treatment changed, and for whom?

## 1 A Number That Does Not Add Up

::: {.case}
Tessera's courier incentive program pays enrolled couriers a bonus on every order. Start with the first data block (block A). The platform operates in 160 cities, observed for two months: $t = 0$ is the month before the program opens in new cities, $t = 1$ the month after, and the outcome is a courier's weekly completed orders (a four-week average). The program opens in 80 cities at $t = 1$. Before the launch, a group of invited beta-testing couriers is already in the program; after it, more couriers enroll voluntarily, and the enrollment rate rises from 0.082 to 0.536. The other 80 cities opened an older program two years ago, about three in ten of their couriers enrolled long ago, and their enrollment rate does not budge over the two months: it is 0.301 in both. The two groups contain 19,424 couriers in total, and enrollment is self-selected in both: couriers with larger expected gains, and couriers who were more hardworking to begin with, are more likely to enroll. One more thing: at $t = 1$ the platform raises the per-order bonus for every enrolled courier in every city, so the effect for incumbent enrollees rises too.

Because this is a simulation, we know the answer. For the couriers in the launch cities who enroll because the program opens (the switchers), the average treatment effect at $t = 1$ is 15.857 orders per week.
:::

The researcher does the textbook thing: use "launch city times post-launch" as an instrument for enrollment and run 2SLS with city and month fixed effects. With two groups and two periods, this is the DiD of the outcome divided by the DiD of the enrollment rate, that is, the Wald-DID. The result is 12.002 orders per week. Its 95% city bootstrap interval, [11.290, 12.774], excludes the truth; the estimate falls 3.855 orders short of it, about 9.9 standard errors.

The troubling thing about this number is that everything looks clean. In the DGP, parallel trends for the untreated outcome holds exactly, the control group's enrollment rate is exactly constant, there are no spillovers, and the sample is not small. Nothing is wrong with how the regression is written. The problem lies in what Wald-DID quietly requires beyond parallel trends: the effect for incumbent enrollees must not change over time, and when the control group's enrollment rate moves, the switchers in the two groups must also have the same effect. Tessera raises its bonus at $t = 1$, so the effect for incumbent enrollees rises, and the control group happens to have far more incumbent enrollees than the launch group. Section 3.2 works through the accounting. On the same data, the two alternative estimators of Sections 3.3 and 3.4, Wald-TC and Wald-CIC, give 15.857 and 15.947, both close to the truth.

## 2 The Economic Model and the Estimand

The three treatment forms correspond to three estimands. We write them down first and turn to identification afterward.

### 2.1 Partial Take-Up: The LATE of Switchers

The notation is as follows. $S_i \in \{0,1\}$ denotes the group, with $S_i = 1$ the group whose treatment rate rises more; $T_i \in \{0,1\}$ denotes the period of the observation, with 0 the base period, which corresponds to $t = 1$ in the 2x2 of Chapter 19 and is unrelated to the number of periods observed for unit $i$, which the panel-data chapters write as $T_i$; $D_i$ denotes the treatment actually received. This section does not use $G_i$, to avoid confusion with the cohort of Chapter 19; the 2x2 of Chapter 19 let $D_i$ double as the group, whereas here group and treatment are kept apart. Conditional quantities for cell $(s,t)$ are written $Y_{st}$ and $D_{st}$, with the group as the first subscript and the period as the second; for example, $\mathbb{E}(Y_{11}) = \mathbb{E}(Y_i \mid S_i = 1, T_i = 1)$. In Tessera's block A, $S_i = 1$ marks the launch cities and $D_i$ is enrollment.

Fuzzy means $D_i \neq S_i T_i$: some units in the control group are treated too, and some units in the treated group remain untreated in the later period; the treated group's treatment rate simply rises more. Let $D_i(t)$ be the treatment status unit $i$ would have if observed in period $t$; the argument in parentheses is the period, not the instrument value of Chapter 18. By $\big(D_i(0), D_i(1)\big)$, the treated group contains three types: always-takers (enrolled in both periods; in Tessera, the beta testers), **switchers** (not enrolled in the base period, enrolled in the later period), and never-takers (enrolled in neither period). The target of this chapter is the LATE of the treated-group switchers in the later period:

$$\mathrm{LATE}_{sw} = \mathbb{E}\big[Y_i(1) - Y_i(0) \mid S_i = 1,\, T_i = 1,\, D_i(0) < D_i(1)\big]$$

Why this group? Because the program launch changes the treatment of these people only: without the program they would not have enrolled, and whether anyone else enrolls has nothing to do with the launch. This is the same logic as the compliers of Chapter 18, except that what pushes treatment here is not a random coin but the combination of "treated group" and "later period." The true value in Tessera's block A is 15.857.

The most widely used estimator is the Wald-DID. For any variable $R$, let $DID_R = \mathbb{E}(R_{11}) - \mathbb{E}(R_{10}) - \big[\mathbb{E}(R_{01}) - \mathbb{E}(R_{00})\big]$. Then

$$W_{DID} = \frac{DID_Y}{DID_D}$$

It is numerically identical to the IV regression that uses $Z_i = S_i T_i$ as the instrument, with group and period as exogenous controls (de Chaisemartin and D'Haultfœuille 2018 §2.2). The question of Section 3 is exactly this: under what conditions does $W_{DID}$ equal $\mathrm{LATE}_{sw}$?

### 2.2 Dose: $\mathrm{ATT}(d \mid d)$ and ACRT

The second data block (block B) has 500 cities and two periods ($t = 1, 2$). In period 2, 370 cities receive a per-order bonus, with amounts ranging from \$0.20 to \$3.00, and the remaining 130 cities receive no bonus (the stayers). Let $D_i \geq 0$ denote the dose a city receives, with lowercase $d$ a value it takes, and let $Y_t(d)$ be the potential outcome at dose $d$ in period $t$. Section 5.1 of Chapter 19 used the term "dose-response" when discussing Card and Krueger (1994), and it referred to exactly the $d \mapsto Y_t(d)$ here; this chapter keeps that term.

Callaway, Goodman-Bacon and Sant'Anna (2025; hereafter CGBS; the paper is forthcoming in the American Economic Review, and this chapter's numbering follows its December 2025 arXiv working-paper version v8) separate two kinds of building blocks:

$$\mathrm{ATT}(d \mid d) = \mathbb{E}\big[Y_2(d) - Y_2(0) \mid D = d\big], \qquad \mathrm{ACRT}(d \mid d) = \frac{\partial\, \mathbb{E}[Y_2(d') \mid D = d]}{\partial d'}\Big|_{d' = d}$$

Both sides of the vertical bar are doses, not the cohort and period of $ATT(g,t)$ in Chapter 19. In general, write $\mathrm{ATT}(d \mid d') = \mathbb{E}[Y_2(d) - Y_2(0) \mid D = d']$: the left of the bar is a hypothetical dose, and the right marks the group of cities that actually received $d'$. $\mathrm{ATT}(d \mid d)$ is the average effect, relative to zero dose, for the cities that received dose $d$: the height of the curve. $\mathrm{ACRT}(d \mid d)$ (average causal response on the treated) is these cities' average marginal effect of a little more dose at $d$: the slope of the curve. Dropping the conditioning on $D = d$ and averaging over all treated cities instead gives $\mathrm{ATT}(d)$ and $\mathrm{ACRT}(d)$. Aggregating over the dose distribution of the treated cities gives $\mathrm{ATT}^{loc} = \mathbb{E}[\mathrm{ATT}(D \mid D) \mid D > 0]$ and $\mathrm{ACRT}^{loc} = \mathbb{E}[\mathrm{ACRT}(D \mid D) \mid D > 0]$; the same average of $\mathrm{ACRT}(D)$ gives $\mathrm{ACRT}^{glob}$. The true values in Tessera's block B are $\mathrm{ATT}^{loc} = 4.858$ and $\mathrm{ACRT}^{loc} = 1.854$. The two are in different units: the first is the total effect of a bonus averaging about \$1.60, the second the effect of \$1 more on top of the existing dose. The slope is far smaller than "total effect divided by dose" because the dose response is concave: the initial bonus does the most, and each further increment adds less.

Dividing the height by the dose gives two other common summary measures (de Chaisemartin, D'Haultfœuille, Pasquier, Sow and Vazquez-Bare, 2026; hereafter dCDH-PSV). When everyone has zero dose in period 1, the average of slopes and the weighted average of slopes are, respectively,

$$\mathrm{AS} = \mathbb{E}\Big[\frac{Y_2(D) - Y_2(0)}{D} \,\Big|\, D > 0\Big], \qquad \mathrm{WAS} = \frac{\mathbb{E}[Y_2(D) - Y_2(0) \mid D > 0]}{\mathbb{E}[D \mid D > 0]}$$

AS averages each city's "effect per USD of bonus" with equal weights. WAS weights by dose and equals $\mathrm{ATT}^{loc}$ divided by the average dose; it answers the question "on average, how many orders does each USD of per-order bonus buy?" When cost is proportional to dose, comparing it with the cost of one unit of dose is a cost-benefit analysis (dCDH-PSV 2026 §4.1). The true values in Tessera are AS 2.923 and WAS 3.038.

### 2.3 Switching Treatments: Actual-versus-Status-Quo Effects

The third data block (block C) has 200 cities over 30 weeks, with a weekly bonus tier $D_{it} \in \{0, 1, 2\}$ that can be switched on, raised, lowered, paused, and resumed. A city's outcome in week $t$ may then depend on the whole treatment path $\mathbf{d} = (d_1, \ldots, d_t)$, and the potential outcome is written $Y_{it}(\mathbf{d})$: the bonus takes a few weeks to build up its effect, and the effect lingers for a few weeks after a pause.

The approach of dCDH (2026) is to compare each city with its own "status quo" (the numbering and numbers this chapter cites from that paper follow its May 2026 arXiv working-paper version v14). Let $G_i$ be the period in which city $i$'s treatment first changes; for an absorbing treatment that starts from 0 this is the cohort of Chapter 19, and here it generalizes to the period of the first switch. If a city kept its period-1 tier throughout, its outcome in period $t$ would be the **status-quo outcome** $Y_{it}(D_{i1}, \ldots, D_{i1})$. The target is

$$\mathrm{AVSQ}_{i,\ell} = \mathbb{E}\big[Y_{i,G_i - 1 + \ell} - Y_{i,G_i - 1 + \ell}(D_{i1}, \ldots, D_{i1})\big]$$

that is, the difference between the actual outcome and the status-quo outcome in the $\ell$-th period after the first switch (actual-versus-status-quo). Here $\ell$ is the horizon of the effect, not the log-likelihood $\ell(\theta)$ of Part II; as in Chapter 19, $\ell = e + 1$, and $\ell = 1$ is the period of the first switch. Averaging over the cities observed at horizon $\ell$ (with a minus sign for cities whose treatment falls) gives $\mathrm{AVSQ}_\ell$, the average effect of "exposure to a dose no lower than the status quo for $\ell$ periods." In Tessera's block C every city starts from 0, so every switcher switches upward; the true values of $\mathrm{AVSQ}_1$ through $\mathrm{AVSQ}_8$ rise from 3.547 to 6.672 at horizon 4 and then decline slowly to 5.949.

Dividing $\mathrm{AVSQ}_{i,\ell}$ by the extra dose accumulated over these $\ell$ periods, $\sum_{k=0}^{\ell-1}(D_{i,G_i+k} - D_{i1})$, gives the normalized effect $\mathrm{AVSQ}^n_{i,\ell}$. In aggregate, dividing $\mathrm{AVSQ}_\ell$ by the switchers' average cumulative increment (in absolute value for those whose dose falls) gives the normalized effect at horizon $\ell$, $\mathrm{AVSQ}^n_\ell$, which is what $DID^n_\ell$ in Section 6.4 estimates. Summing $\mathrm{AVSQ}_{i,\ell}$ over all switchers and all horizons and dividing by the total extra dose gives the average total effect per unit of treatment; its true value in Tessera is 5.038. Another benchmark often set against TWFE is "the sum of effects divided by the sum of doses" over all treated cells, with a true value of 5.278.

To summarize: the estimand in each of the three forms is the average effect for some group of "units whose treatment was actually changed": the LATE of switchers, the height and slope of the dose curve, and the path effect relative to the status quo. None of them is the population ATE, and each must first state "relative to what": relative to no program, relative to zero dose, or relative to keeping the period-1 tier.

## 3 Identification: Which Part of the Treatment Changed, and for Whom

Identification asks whether, under a set of assumptions about the counterfactual, the causal quantities of Section 2 can be written as functions of the distribution of observables. This section follows the order of Chapter 19: first look at what is to be proved through numbers or a figure, then state numbered assumptions and theorems, followed by a proof sketch that names the assumption used up at each step, and finally answer whose effect the quantity estimates under which assumptions, who determines the weights, and where the assumptions are most likely to fail. In the fuzzy part, A1 through A8 have the same content as Assumptions 1 through 8 of the journal version of dCDH (2018) in a different order, and each is followed by its journal-version number.

### 3.1 The Fuzzy DiD Setup and What Wald-DID Identifies

First see how Wald-DID arrives at 12.002. Average weekly orders in the launch cities rise from 74.776 to 86.572, and in the control cities from 73.925 to 80.263, so the DiD of the outcome is 5.457. The DiD of the enrollment rate is 0.455, because the launch cities' rate rises by that much while the control cities' rate stays put. Dividing one by the other gives 12.002. The intuition is the same as for the Wald ratio of Chapter 18: the numerator is "the change in orders brought about by the launch," and the denominator is "the change in enrollment brought about by the launch." The difference is that the instrument of Chapter 18 was randomly assigned, whereas the "instrument" here is group times period: the two groups differ in levels, and time itself also changes the outcome, so assumptions of the parallel trends kind have to fill the gap.

::: {.assumption}
**A1 (treatment participation; dCDH 2018 Assumption 3)** $D_i = \mathbf{1}\{V_i \geq v_{S_i T_i}\}$, and $V_i \perp T_i \mid S_i$.

**A2 (first stage; Assumption 1)** $\mathbb{E}(D_{11}) > \mathbb{E}(D_{10})$, and $\mathbb{E}(D_{11}) - \mathbb{E}(D_{10}) > \mathbb{E}(D_{01}) - \mathbb{E}(D_{00})$.
:::

A1 is a latent index model: $V_i$ is the propensity to enroll, the threshold varies by group and period, and the distribution of propensities within a group does not change over periods. Accordingly, write $D_i(t) = \mathbf{1}\{V_i \geq v_{S_i t}\}$. Within a group the threshold moves in only one direction from the base period to the later period, so switches within the group can go in only one direction, and the classification into types is well defined. In a panel where the propensity does not change over time, the "distribution does not change over periods" half holds automatically, but requiring all switches within a group to go the same way is still a real restriction; the paper therefore says this assumption is better suited to repeated cross sections or cohort data (de Chaisemartin and D'Haultfœuille 2018 §3.4, paragraph "Panel data"). In Tessera's block A the only switches are from not enrolled to enrolled, which satisfies this. A2 merely defines which group is the treated group.

::: {.assumption}
**A3 (common trends, that is, the parallel trends of Chapter 19; Assumption 4)** $\mathbb{E}[Y_i(0) \mid S_i, T_i = 1] - \mathbb{E}[Y_i(0) \mid S_i, T_i = 0]$ does not depend on $S_i$.

**A4 (effect stable over time; Assumption 5)** $\mathbb{E}[Y_i(1) - Y_i(0) \mid S_i, T_i = 1, D_i(0) = 1] = \mathbb{E}[Y_i(1) - Y_i(0) \mid S_i, T_i = 0, D_i(0) = 1]$.
:::

::: {.assumption}
**A5 (homogeneous effect between groups; Assumption 6)** $\mathrm{LATE}_{sw} = \mathrm{LATE}'_{sw}$, where $\mathrm{LATE}'_{sw} = \mathbb{E}[Y_i(1) - Y_i(0) \mid S_i = 0, T_i = 1, D_i(0) \neq D_i(1)]$ is the effect for the control-group switchers.

**A6 (stable control-group treatment rate; Assumption 2)** $0 < \mathbb{E}(D_{01}) = \mathbb{E}(D_{00}) < 1$.
:::

A3 is a direct generalization of A1 of Chapter 19 and restricts only $Y(0)$. A4 and A5 are specific to the fuzzy setting and are dissected in the next subsection. Of A3 through A6, only A6 can be tested directly with data: the control group's treatment rate is the same in both periods (0.301 in Tessera's block A).

::: {.theorem}
**Theorem 1 (Wald-DID; dCDH 2018 Theorem 1)** Let $\alpha = [\mathbb{E}(D_{11}) - \mathbb{E}(D_{10})] / DID_D$.

(i) Under A1 through A4, $W_{DID} = \alpha\, \mathrm{LATE}_{sw} + (1 - \alpha)\, \mathrm{LATE}'_{sw}$;

(ii) if A5 or A6 also holds, $W_{DID} = \mathrm{LATE}_{sw}$.
:::

The proof sketch has four steps (see the paper's appendix). Step one (A1 and A2): the change in the treated group's enrollment rate, $\mathbb{E}(D_{11}) - \mathbb{E}(D_{10})$, equals exactly the share of switchers in the treated group, because the propensity distribution does not change over periods, so the extra enrollees are exactly those swept in by the falling threshold. Step two (A1 and A4): the change in the treated group's outcome splits into three pieces: the switchers' effect times their share, the trend in $Y(0)$, and the change in the always-takers' effect times their share; A4 makes the last piece zero. Step three (A1 and A4): splitting the control group the same way gives the control-group switchers' effect times the change in the control group's treatment rate, plus the trend in $Y(0)$. Step four (A3): the two groups' $Y(0)$ trends cancel on subtraction, $DID_Y = \mathrm{LATE}_{sw} \cdot \Pr(\text{switcher} \mid S = 1) - \mathrm{LATE}'_{sw} \cdot [\mathbb{E}(D_{01}) - \mathbb{E}(D_{00})]$, and dividing by $DID_D$ gives (i). A6 sets $\alpha = 1$ and A5 makes the two effects equal; either one gives (ii).

In words: even when parallel trends holds, Wald-DID isolates the treated group's switchers only if the control group's treatment rate stays put. When the control group's treatment rate is rising, $\alpha > 1$, and Wald-DID is a weighted difference of the two groups' switcher effects, which can turn negative even when everyone's effect is positive. And whether or not the control group moves, Wald-DID cannot do without A4.

### 3.2 Why the Two Extra Assumptions Often Fail

A4 requires the average effect for the incumbent treated (those already treated in the base period) to be the same in both periods; equivalently, their $Y(1)$ and $Y(0)$ evolve in step. A5 requires the switchers in the two groups to have the same effect. Both are needed because Wald-DID uses the control group's change in outcome as its ruler, and that change contains not only the trend in $Y(0)$ but also the change in the effect for the control group's incumbent treated and the effect for its newly treated. When A6 holds the latter term vanishes, but the former remains. Without A4, dCDH (2018) give the expression for the bias at the end of their proof:

::: {.theorem}
**Corollary 1 (bias when A6 holds but A4 does not)** Let $\mathrm{ATT}^{a}_{st} = \mathbb{E}[Y_i(1) - Y_i(0) \mid S_i = s, T_i = t, D_i(0) = 1]$ be the average effect for the incumbent treated in cell $(s,t)$. Under A1 through A3 and A6,

$$W_{DID} = \mathrm{LATE}_{sw} + \frac{\big(\mathrm{ATT}^{a}_{11} - \mathrm{ATT}^{a}_{10}\big)\, \mathbb{E}(D_{10}) - \big(\mathrm{ATT}^{a}_{01} - \mathrm{ATT}^{a}_{00}\big)\, \mathbb{E}(D_{00})}{DID_D}$$
:::

The proof sketch has two steps. Step one (A1 and A2): split the two groups' changes in outcome as in the first three steps of the proof of Theorem 1, but without using A4, so that the change in each group's incumbent effect, times that group's base-period treatment rate, stays in its change in outcome. Step two (A3 and A6): A6 removes the control-group switcher term, A3 cancels the two groups' $Y(0)$ trends, and subtracting and dividing by $DID_D$ gives the result.

In words: Wald-DID treats the change in the effect for the control group's incumbent treated as part of the trend. Tessera is exactly this case. At $t = 1$ the platform raises the bonus for all enrolled couriers, and the effect for incumbent enrollees in both groups rises by 8 orders. In the population, incumbent enrollees make up 0.08 of the launch cities and 0.30 of the control cities, and switchers make up 0.46, so the bias equals 8 times the difference between 0.08 and 0.30, divided by 0.46, that is, $-3.826$: the population Wald-DID is 12.042, while the population $\mathrm{LATE}_{sw}$ is 15.868 (the sample truth of 15.857 used elsewhere in this chapter differs from it only by sampling error). The actual bias in the sample, $-3.855$, almost coincides with this closed form.

The reason for the negative sign is plain to see: the control group has more incumbent enrollees, so the "rising incumbent effect" piece of its change in outcome is larger, and treating it as trend subtracts too much from the launch group's change. Conversely, if the launch group had more incumbent enrollees, the bias would be positive.

These two assumptions often fail in practice, for plain reasons. Effects have many reasons to drift over time: the platform adjusts prices and bonuses, couriers become more skilled with experience, competitors enter and exit, and the seasons change the mix of orders. In revisiting Duflo (2001), dCDH (2018) point out that returns to experience are diminishing. A4 requires the wage gap between high school graduates of the two birth cohorts to stay the same had they completed only primary school; but with only primary school they would have entered the labor market earlier and had more experience by 1995, when wages are measured, and with diminishing returns to experience the gap would then have been smaller. The wage gap between the two cohorts therefore differs between the "educated" and "uneducated" states, and A4 is hard to sustain (Section 5.1). A5 requires the switchers in the two groups to be the same kind of people, but if some control-group units are switching too, they have often been pushed in by a different mechanism.

Testability is also limited: if an earlier period is available and both groups' treatment rates were stable then, A3 and A4 together imply that the placebo DiD over the earlier periods is zero (de Chaisemartin and D'Haultfœuille 2018 §2.5), which can catch part of the problem; A5, by contrast, can hardly be tested with placebos.

### 3.3 Wald-TC: Parallel Trends by Baseline Treatment Status

The trouble with Wald-DID is that it treats the change in the incumbent effect as trend. A natural repair is to measure the trends of the two kinds of people separately. In the control group, those not enrolled in the base period tell us how $Y(0)$ moves, and those already enrolled tell us how $Y(1)$ moves. Adding each trend to the corresponding people in the treated group's base period gives the treated group's average outcome in the later period "had the switchers not switched."

::: {.assumption}
**A3' (conditional common trends; Assumption 4')** For $d \in \{0, 1\}$, $\mathbb{E}[Y_i(d) \mid S_i, T_i = 1, D_i(0) = d] - \mathbb{E}[Y_i(d) \mid S_i, T_i = 0, D_i(0) = d]$ does not depend on $S_i$.
:::

A3' requires $Y(0)$ to be parallel among those untreated in the base period, and $Y(1)$ to be parallel among those already treated in the base period. It does not require the effect to be stable over time: the incumbent effect may rise, as long as it rises by the same amount for the incumbents of both groups. Let $\delta_d = \mathbb{E}(Y_{01} \mid D = d) - \mathbb{E}(Y_{00} \mid D = d)$ be the change in outcome for control-group units with treatment status $d$. The time-corrected Wald ratio is

$$W_{TC} = \frac{\mathbb{E}(Y_{11}) - \mathbb{E}\big(Y_{10} + \delta_{D_{10}}\big)}{\mathbb{E}(D_{11}) - \mathbb{E}(D_{10})}$$

where $\delta_{D_{10}}$ means adding to each person in the treated group's base period the trend that matches that person's treatment status.

::: {.theorem}
**Theorem 2 (Wald-TC; dCDH 2018 Theorem 2)** Under A1, A2, A3', and A6, $W_{TC} = \mathrm{LATE}_{sw}$.
:::

The proof sketch has three steps (equations (A9) to (A11) in the paper's appendix). Step one (A1 and A2): the treated group's change in outcome equals the switchers' effect times their share, plus the $Y(1)$ trend of those already treated in the base period times their share, plus the $Y(0)$ trend of those untreated in the base period times their share. Step two (A1 and A6): a constant control-group treatment rate, together with a propensity distribution that does not change over periods, means that the $D = 1$ cell of the control group contains the same people in both periods, and so does the $D = 0$ cell; hence $\delta_1$ is exactly the $Y(1)$ trend of the control group's incumbent treated, and $\delta_0$ exactly the $Y(0)$ trend of its untreated. Step three (A3'): these two trends equal the trends of the corresponding populations in the treated group; substituting into step one leaves only the switchers' effect times their share in the numerator, and the denominator is exactly that share.

In words: Wald-TC still estimates the LATE of the treated group's switchers, with equal weight on every switcher. It does not need A4; the price is two new requirements, parallel trends by baseline treatment status and a constant control-group treatment rate. In Tessera's block A, $\delta_0 = 3.933$ measures the common seasonal rise in $Y(0)$; $\delta_1 = 11.911$ measures that rise plus the increase in the incumbent enrollees' effect; and $W_{TC} = 15.857$, only 0.0001 away from the truth. A3' is most likely to fail when couriers enroll because of a positive shock to $Y(1)$: if the treated group's enrollment rate rises in the later period because that group's $Y(1)$ improved in the later period, A3' and A4 fail together, and Wald-DID and Wald-TC are both inconsistent.

### 3.4 Wald-CIC, Changes-in-Changes, and Partial Identification

Wald-TC measures trends by additive shifts. What if a trend does not act on the distribution as a uniform shift? For example, a higher bonus may make couriers who already deliver a lot deliver even more, while low-volume couriers barely change; then a shift in the mean no longer applies to every individual. Changes-in-changes (Athey and Imbens 2006) takes a different approach: matching by rank. For the full identification of the CIC model (Theorem 10) and the LQTE of switchers (Theorem 11), see Section 3.9 of Chapter 13; A1, A2, A6, A7, and A8 of this chapter are exactly FZ.1 through FZ.5 there.

First, the bare minimum of quantile machinery (Chapter 13 is devoted to quantile regression); only two objects are needed here. $F_Y(y) = \Pr(Y \leq y)$ is the distribution function, which tells you where $y$ ranks in the distribution (as a proportion). The quantile function $F^{-1}_Y(u) = \inf\{y : F_Y(y) \geq u\}$, $u \in (0,1)$, works the other way and tells you which value sits at rank $u$. Suppose the median among the control group's non-enrolled couriers is 60 orders in the base period and 64 in the later period. Then for a non-enrolled courier in the treated group who delivers 60 orders in the base period, CIC predicts 64 orders in the later period if that courier still does not enroll: first find the same order count in the control group's base period and read off its rank (the median), then read off the value at the same rank in the control group's later period. This two-step mapping is

$$Q_d(y) = F^{-1}_{Y_{01} \mid D = d}\big(F_{Y_{00} \mid D = d}(y)\big)$$

Do this once for $d = 0$ and once for $d = 1$ (the mapping for $d = 0$ is written $Q_{d=0}$). Wald-CIC is

$$W_{CIC} = \frac{\mathbb{E}(Y_{11}) - \mathbb{E}\big[Q_{D_i}(Y_i) \mid S_i = 1, T_i = 0\big]}{\mathbb{E}(D_{11}) - \mathbb{E}(D_{10})}$$

::: {.assumption}
**A7 (monotonicity and time invariance of unobservables; Assumption 7)** $Y_i(d) = h_d(U_{id}, T_i)$, $h_d(u, t)$ is strictly increasing in $u$, and $U_{id} \perp T_i \mid S_i, D_i(0)$.

**A8 (data restrictions; Assumption 8)** The outcome has the same support, a closed interval, in each of the 8 treatment-by-group-by-period cells, and each cell's distribution function is continuous and strictly increasing on the support.
:::

A7 extends the CIC model of Athey and Imbens (2006), in which the outcome is a strictly increasing function of a scalar unobservable whose distribution within a group does not change over time, to the fuzzy setting: the distribution of the unobservable must not change over time within each subgroup defined by "group times baseline treatment status." It is invariant to any monotone transformation of the outcome, which is its advantage over A3'; the price is that it restricts the whole distribution, whereas A3' restricts only the mean.

::: {.theorem}
**Theorem 3 (Wald-CIC; dCDH 2018 Theorem 3)** Under A1, A2, A6, A7, and A8, $W_{CIC} = \mathrm{LATE}_{sw}$; the local quantile treatment effects (LQTE) of the switchers are also identified.
:::

The proof sketch has three steps. Step one (A1, A6, and A7): the control group's treatment rate is constant, so the $D = d$ cell contains the same people in both periods and the distribution of the unobservable does not change; hence $Q_d$ maps the outcome at rank $u$ in the control group's base period to the outcome at rank $u$ in its later period, that is, $h_d(\cdot, 1) \circ h_d^{-1}(\cdot, 0)$. Step two (A7 and A8): a person in the treated group with base-period treatment status $d$ and outcome $y = h_d(u, 0)$, and a person in the control group's base period with the same status $d$ and the same outcome $y$, have the same $u$; common support guarantees that this match exists, so the person's later-period outcome under unchanged treatment is $h_d(u, 1) = Q_d(y)$. Step three (A1 and A2): $\mathbb{E}[Q_{D_i}(Y_i) \mid S_i = 1, T_i = 0]$ is therefore the treated group's average later-period outcome "if the switchers did not switch," and subtracting and dividing by the switchers' share gives the result.

In words: Wald-CIC and Wald-TC estimate the same LATE and differ only in how they carry base-period outcomes over to the later period: Wald-TC by a shift in means, Wald-CIC by matching ranks. When the two groups' base-period outcome distributions are far apart, the scale of the outcome strongly affects Wald-TC, and Wald-CIC is preferable; when the distributions are similar, Wald-TC, which restricts only the first moment, is the safer choice (de Chaisemartin and D'Haultfœuille 2018 §2.5). Tessera's DGP is a pure shift model, so both assumptions hold; $W_{CIC} = 15.947$, 0.090 above the truth, and the truth remains inside its interval.

Both estimators depend on A6. When the control group's treatment rate moves, the $D = d$ cell no longer contains the same people in both periods, and $\delta_d$ measures not only the trend but also a change in composition.

::: {.theorem}
**Theorem 4 (partial identification; dCDH 2018 Theorem 4)** Suppose A1 and A2 hold and $0 < \mathbb{E}(D_{01}) \neq \mathbb{E}(D_{00}) < 1$. (i) If A3' holds and the outcome is bounded, $\mathrm{LATE}_{sw}$ lies within the TC bounds; (ii) if A7, A8, and a moment condition hold, $\mathrm{LATE}_{sw}$ lies within the CIC bounds, and under a further monotonicity condition these bounds are sharp. When the support of the outcome is unbounded, the TC bounds are always infinite.
:::

The proof sketch has two steps. Step one (A1): when the control group's treatment rate moves, the two treatment cells no longer contain the same people in both periods, and the difference is exactly the control-group switchers; since switches go in only one direction, their share in a cell can be computed from the ratio of the two periods' shares, though we do not know who they are. Step two (A3' and a bounded outcome, or A7, A8, and the moment condition): placing these people at either end of the outcome distribution and trimming for the best and worst cases gives the largest and smallest possible values of $\delta_d$ (or of the $Q_d$ mapping), and substituting into $W_{TC}$ (or $W_{CIC}$) gives the upper and lower bounds.

In words: once A6 fails, the target is still the LATE of the treated group's switchers, but the data can only box it into an interval, and the smaller the change in the control group's treatment rate, the narrower the interval. Whether the interval is finite also depends on the support of the outcome. The paper notes that with an unbounded outcome, the CIC bounds are still finite if the control group's untreated share rises; if that share falls, the distribution functions used to construct the bounds are missing probability mass.

Variant A′ turns this point into numbers. Same couriers, same launch cities, same true value of 15.857; the only change is that the control group runs a recruitment drive at $t = 1$ in which the non-enrolled couriers with the highest order levels enroll, raising the enrollment rate from 0.301 to 0.402. A3' still holds; A6 fails. Now $W_{TC} = 21.349$ and $W_{CIC} = 22.283$, 5.493 and 6.426 above the truth, and both intervals are far from it; in a large-sample simulation $W_{TC}$ is about 21.147, so this is not noise. The two rulers reveal the root of the problem: the high-volume couriers leave the control group's non-enrolled cell, so $\delta_0$ falls from 3.933 to 1.116 and the trend in $Y(0)$ is badly underestimated; they join the enrolled cell, and $\delta_1$ rises to 13.003. Decomposing the bias with large-sample cell means, "the change in composition of the non-enrolled cell" contributes 5.460 and "the change in composition of the enrolled cell" contributes $-0.156$. In A′ the outcome is unbounded and the control group's untreated share is falling, so the TC bounds are infinite, and the CIC bounds do not fall under the case in which the paper guarantees finiteness either.

::: {.warning}
A6 is not a technicality; it is the precondition for point identification by Wald-TC and Wald-CIC. Before using these two estimators, report the control group's treatment rates in both periods and test their difference. When the treatment rate moves, their point estimates can be far from the truth, as in A′, while the intervals look narrow; all that can be reported honestly then are the partial identification bounds.
:::

### 3.5 Continuous Treatments: $\mathrm{ATT}(d \mid d)$, ACRT, and Strong Parallel Trends

Block B poses the problem in a different form. Picture two cities, one with a bonus of \$1 per order and the other with \$2, where orders rise more in the second. Can we say that "raising the bonus from \$1 to \$2" has an effect this large? No: the platform gave the higher bonus to the city whose couriers are more responsive, so even at \$1 the couriers in the \$2 city would have increased their orders more. In block B the gain coefficient averages 0.82 in cities receiving \$1 and 1.12 in cities receiving \$2, so the same bonus yields an increment more than a third larger in the latter. Comparisons across doses are contaminated by selection.

::: {.assumption}
**A9 (parallel trends; CGBS Assumption PT)** For all $d > 0$, $\mathbb{E}[Y_2(0) - Y_1(0) \mid D = d] = \mathbb{E}[Y_2(0) - Y_1(0) \mid D = 0]$.

**A10 (strong parallel trends; CGBS Assumption SPT)** For all $d$, $\mathbb{E}[Y_2(d) - Y_1(0) \mid D > 0] = \mathbb{E}[Y_2(d) - Y_1(0) \mid D = d]$.
:::

A9 is a direct generalization of the parallel trends of Chapter 19 and restricts only $Y(0)$. A10 restricts $Y(d)$: how all treated cities would change on average if they all received dose $d$ must equal how the cities that actually received $d$ actually changed. Given A9, A10 is equivalent to the equality on the diagonal $\mathrm{ATT}(d \mid d) = \mathrm{ATT}(d)$: for the cities that received dose $d$, the average effect of $d$ equals the average effect of $d$ for all treated cities, which CGBS describe as ruling out selection into dose by the treated units. It does not require $\mathrm{ATT}(d \mid d')$ to be equal everywhere across dose groups $d'$. The name strong parallel trends is a reminder that in most applications it is much stronger than A9.

::: {.theorem}
**Theorem 5 (identification with a continuous treatment; CGBS Theorems 3.1 to 3.3 and Corollary 3.1)** Suppose Assumptions 1 to 3 of CGBS hold (random sampling; no unit treated in period 1 and a positive probability of zero dose; no anticipation), and let $\Delta Y = Y_2 - Y_1$.

(i) Under A9, $\mathrm{ATT}(d \mid d) = \mathbb{E}[\Delta Y \mid D = d] - \mathbb{E}[\Delta Y \mid D = 0]$ and $\mathrm{ATT}^{loc} = \mathbb{E}[\Delta Y \mid D > 0] - \mathbb{E}[\Delta Y \mid D = 0]$;

(ii) under A9 and a continuous dose,

$$\frac{\partial\, \mathbb{E}[\Delta Y \mid D = d]}{\partial d} = \mathrm{ACRT}(d \mid d) + \frac{\partial\, \mathrm{ATT}(d \mid d')}{\partial d'}\Big|_{d' = d}$$

(iii) under A10, the right-hand side of (i) identifies $\mathrm{ATT}(d)$, the left-hand side of (ii) identifies $\mathrm{ACRT}(d)$, and its average over the distribution of treated doses identifies $\mathrm{ACRT}^{glob}$.
:::

The proof sketch has four steps. Step one (no anticipation): $\mathbb{E}[\Delta Y \mid D = d] = \mathrm{ATT}(d \mid d) + \mathbb{E}[Y_2(0) - Y_1(0) \mid D = d]$. Step two (A9): replacing the last term with the change for zero-dose cities gives (i). Step three (differentiation): the total derivative of $\mathrm{ATT}(d \mid d)$ with respect to $d$ splits into two parts: the part that holds the population fixed and changes the dose is $\mathrm{ACRT}(d \mid d)$, and the part that holds the dose fixed and changes the population is local selection bias, which gives (ii). Step four (A10 and no anticipation): applying A10 at dose $d$ and at dose 0 and subtracting gives $\mathbb{E}[\Delta Y \mid D = d] - \mathbb{E}[\Delta Y \mid D = 0] = \mathbb{E}[Y_2(d) - Y_2(0) \mid D > 0] = \mathrm{ATT}(d)$; differentiating with respect to $d$ gives $\mathrm{ACRT}(d)$, and averaging over the distribution of treated doses gives $\mathrm{ACRT}^{glob}$, which is (iii).

Note that A10 does not make $\mathrm{ATT}(d \mid d')$ constant in $d'$, so local selection bias is generally not zero. What A10 does is change the target: under A10 the same derivative identifies $\mathrm{ACRT}(d)$, and when A9 and A10 both hold, it differs from $\mathrm{ACRT}(d \mid d)$ by exactly the local selection bias term in (ii).

In words: ordinary parallel trends suffices to identify the height of each dose group's own effect, but not the slope; differences across dose groups mix the causal response with selection on "who got this dose." Tessera's block B DGP makes A9 hold exactly while A10 fails: the platform gives larger bonuses to cities whose couriers are more responsive, and the conditional mean of the gain coefficient is $1 + 0.3(d - 1.6)$. At $d = 1.5$, the slope that parallel trends identifies is 3.529, of which the true causal response $\mathrm{ACRT}(d \mid d)$ is only 2.134, with local selection bias accounting for 1.395. Averaged over all treated cities, the identifiable slope is 3.186, whereas $\mathrm{ACRT}^{loc}$ is only 1.854. Section 6.3 shows that contdid's overall "ACRT" estimate is 3.526, with an interval that covers 3.186 but excludes 1.854.

::: {.warning}
A dose slope estimated under ordinary parallel trends is not a causal response: it equals $\mathrm{ACRT}(d \mid d)$ plus local selection bias. Only under strong parallel trends is it the $\mathrm{ACRT}(d)$ of all treated cities. Pre-trends tests involve only $Y(0)$ and cannot tell A9 from A10; whether A10 is credible can be argued only from institutional knowledge of how doses are assigned.
:::

### 3.6 Stayers and the Heterogeneous Adoption Design

Theorem 5 relies on a group of zero-dose cities. dCDH-PSV extend the idea to settings in which treatment is continuously distributed in every period: units are split into switchers, whose dose changed, and **stayers**, whose dose did not, and the stayers are the ruler. Let $D_1$ and $D_2$ be the doses in the two periods; in block B, $D_1 = 0$ and $D_2 = D$.

::: {.assumption}
**A11 (parallel trends given the baseline dose; dCDH-PSV Assumption 3)** For all $d_1$, $\mathbb{E}[Y_2(d_1) - Y_1(d_1) \mid D_1 = d_1, D_2] = \mathbb{E}[Y_2(d_1) - Y_1(d_1) \mid D_1 = d_1]$.

**A12 (least treated on parallel trends; de Chaisemartin, Ciccia, D'Haultfœuille and Knau 2026 Assumption 2)** Let $\underline{d} = \inf \mathrm{Supp}(D_2)$. Then $\lim_{d \downarrow \underline{d}} \mathbb{E}[Y_2(0) - Y_1(0) \mid D_2 \leq d] = \mathbb{E}[Y_2(0) - Y_1(0)]$.
:::

A11 requires parallel trends only among units with the same baseline dose, so it places no restriction on how effects change over time. Using stayers with a different baseline dose as the ruler would require the additional assumption that the average effect does not change over time (dCDH-PSV 2026 §2.5), the same lesson as A4 in Section 3.2. A12 serves the **heterogeneous adoption design** (HAD), in which no unit remains at zero dose (the policy applies universally but exposure differs); it requires only that the group with the lowest doses be parallel with the population as a whole. Identification and estimation for HAD come from de Chaisemartin, Ciccia, D'Haultfœuille and Knau (2026; hereafter the HAD paper).

::: {.theorem}
**Theorem 6 (stayers and HAD; dCDH-PSV Theorems 1 and 4, HAD paper Theorem 1)** (i) Under no anticipation, no dynamic effects, A11, boundedness, and overlap,

$$\mathrm{WAS} = \frac{\mathbb{E}\big[\mathrm{sgn}(D_2 - D_1)\big(\Delta Y - \mathbb{E}(\Delta Y \mid D_1, D_2 = D_1)\big)\big]}{\mathbb{E}|D_2 - D_1|}$$

If, in addition, every switcher's change in dose is at least some positive number (no quasi-stayers), AS is identified too, as $\mathbb{E}\big[(\Delta Y - \mathbb{E}(\Delta Y \mid D_1, D_2 = D_1))/(D_2 - D_1) \mid D_2 \neq D_1\big]$.

(ii) In the HAD, if $D_1 = 0$, $D_2 > 0$, a quasi-untreated group exists ($\underline{d} = 0$), A12 holds, and effects are bounded, then

$$\mathrm{WAS} = \frac{\mathbb{E}[\Delta Y] - \lim_{d \downarrow 0} \mathbb{E}[\Delta Y \mid D_2 \leq d]}{\mathbb{E}[D_2]}$$
:::

The proof sketch has three steps. Step one (no anticipation, no dynamic effects, and A11): for a switcher with baseline dose $d_1$, its change in outcome had its dose stayed the same equals the average change of the stayers with the same baseline dose; so the switcher's change minus this quantity is $Y_2(D_2) - Y_2(d_1)$, and overlap guarantees that every baseline dose has stayers. Step two (AS additionally uses no quasi-stayers): dividing by the change in dose gives slopes, and averaging them with equal weights gives AS; summing without dividing and then dividing by the total change in dose gives WAS. If switchers' changes in dose can be arbitrarily close to zero, the denominator of the first approach goes to zero, AS can be identified only through trimming and a limit, and as the smallest change in dose goes to zero, the ratio of the asymptotic variances of AS and WAS goes to infinity (dCDH-PSV Theorem 2 and Proposition 1); WAS has no such problem. Step three (HAD, A12, and bounded effects): for cities whose dose goes to zero, the effect is squeezed to zero by "bounded effect times dose," so their average change converges to the average change in $Y(0)$; A12 then replaces this with the population's change in $Y(0)$, the population's actual change minus it is the average effect, and dividing by the average dose gives (ii).

In words: with stayers, the stayers play the counterfactual; without stayers, cities with near-zero doses play the quasi-untreated. WAS is the quantity needed for cost-benefit analysis; AS is the average slope, and under a concave or convex dose response it can bound the effects of other dose changes (dCDH-PSV 2026 §3.1). The true values in Tessera's block B are AS 2.923 and WAS 3.038. The no-stayer variant turns the 130 zero-dose cities into cities with tiny doses between \$0.005 and \$0.25, with a true WAS of 3.059. HAD is most likely to fail in two places: the cities with the lowest doses are themselves marginal markets (A12 fails), or even a tiny dose already has a sizable effect (an extensive margin), so the assumption of bounded effects does not hold.

### 3.7 Non-Absorbing, Non-Binary Treatments: From $DID_M$ to $DID_\ell$

A city in block C switches on tier 1 in week 8, pauses in week 12, and resumes in week 14. What is its "effect" in week 15? Who can play its counterfactual?

The first answer is the $DID_M$ of dCDH (2020). For each pair of adjacent periods $(t-1, t)$, it compares "groups that switch treatment on in this period" (joiners) with "groups untreated in both periods," compares "groups treated in both periods" with "groups that switch treatment off in this period" (leavers), and then averages these DiDs with weights given by the number of switchers. Under conditions such as potential outcomes that depend only on current treatment, common trends for both $Y(0)$ and $Y(1)$, and the presence of the corresponding stable groups in every period, $\mathbb{E}[DID_M]$ equals the average effect of the switchers in the period they switch (dCDH 2020 Theorem 3). Its limitation lies in this very "depends only on current treatment": when effects take weeks to build and linger after a pause, a group treated in period 1 and untreated in periods 2 and 3 is no longer a clean control, because its past treatment still affects its change in outcome (de Chaisemartin and D'Haultfœuille 2026 §1).

The $DID_\ell$ of dCDH (2026) removes this limitation. For city $i$, let $\mathcal{C}_{i,\ell} = \{j : D_{j1} = D_{i1},\ G_j > G_i - 1 + \ell\}$ be the set of cities with the same period-1 tier whose tier has not yet changed by period $G_i - 1 + \ell$, and let $N^{i}_{\ell}$ be its size:

$$DID_{i,\ell} = \big(Y_{i,G_i - 1 + \ell} - Y_{i,G_i - 1}\big) - \frac{1}{N^{i}_{\ell}} \sum_{j \in \mathcal{C}_{i,\ell}} \big(Y_{j,G_i - 1 + \ell} - Y_{j,G_i - 1}\big)$$

$DID_\ell$ is the average over all cities for which horizon $\ell$ can be computed (with a minus sign for those whose treatment falls).

::: {.assumption}
**A13 (no anticipation; dCDH 2026 Assumption 1)** $Y_{it}(d_1, \ldots, d_{\mathcal{T}}) = Y_{it}(d_1, \ldots, d_t)$.

**A14 (parallel trends for the status-quo outcome; Assumption 2)** If $D_{i1} = D_{j1}$, then for all $t \geq 2$ the expected changes in the status-quo outcomes of $i$ and $j$ from $t - 1$ to $t$ are equal.
:::

In A13, $\mathcal{T}$ is the last period of the panel, as in Chapter 19. A14 restricts only one potential outcome path per city, the one that keeps the period-1 tier, so by itself it places no restriction on the heterogeneity or dynamics of effects. It compares only cities with the same period-1 tier. If cities with different period-1 tiers were compared, parallel trends would quietly require that "the effect of being treated for $t$ periods equals the effect of being treated for $t-1$ periods," which rules out lagged and time-varying effects (de Chaisemartin and D'Haultfœuille 2026 §3.1). This is the same trap as A4 for Wald-DID in Section 3.2.

::: {.theorem}
**Theorem 7 ($DID_\ell$; dCDH 2026 Lemma 1, Lemma 2, and Theorem 1)** Under A13 and A14: (i) conditional on the treatment design, and when $\mathcal{C}_{i,\ell}$ is nonempty, $\mathbb{E}[DID_{i,\ell}] = \mathrm{AVSQ}_{i,\ell}$, so $DID_\ell$ is conditionally unbiased for $\mathrm{AVSQ}_\ell$; (ii) if every city's treatment is always at least (or always at most) its period-1 tier, $\mathrm{AVSQ}^n_{i,\ell}$ is a weighted average of the effects of current treatment and of its first $\ell - 1$ lags, with nonnegative weights that sum to one and are proportional to each period's dose increment relative to the status quo; (iii) under independence across cities and regularity conditions, $DID_\ell$ is asymptotically normal as the number of cities goes to infinity, and the conventional confidence intervals are asymptotically conservative.
:::

The proof sketch has three steps. Step one (A13): in period $G_i - 1$ city $i$ has not yet switched, and no anticipation guarantees that its outcome in that period is the status-quo outcome; the cities in the control set have not switched by period $G_i - 1 + \ell$, so their outcomes in both periods are status-quo outcomes too. Step two (A14): the average change in the control set is an unbiased substitute for the change city $i$ would have experienced under the status quo, and subtracting gives $\mathrm{AVSQ}_{i,\ell}$. Step three (telescoping): write $\mathrm{AVSQ}_{i,\ell}$ as the sum of the increments produced by changing each period's treatment one at a time; each term divided by the corresponding dose increment is a slope, so dividing by the cumulative increment gives a weighted average of slopes, and "treatment never crosses the status quo" guarantees that the weights are nonnegative.

In words: $DID_\ell$ estimates the switchers' effect relative to the status quo in the $\ell$-th period after the first switch, with equal weight on every switcher; the control group uses only cities with the same period-1 tier that have not yet switched, never cities that have already switched. In the kind of design studied in Chapter 19 (binary, absorbing, every unit starting from 0), it is the Callaway-Sant'Anna event-study coefficient with not-yet-treated controls. In Tessera's block C, $DID_1$ is 3.628 with interval [3.120, 4.137], covering the true value of 3.547, and all eight horizons are covered (Section 6.4). Its cost lies in interpretation: at the same $\ell$, different cities have traveled different paths, and $DID_2$ may average "switched on, then off" together with "switched on, then raised," so the paper recommends also reporting the distribution of paths.

### 3.8 A Bridge to the IV/LATE of Chapter 18: Instrumented DiD

Return to block A. With $Z_i = S_i T_i$ as the instrument, Wald-DID is the DiD version of the Wald ratio in the LATE theorem of Chapter 18 (Imbens and Angrist 1994), with one DiD in the numerator and one in the denominator, except that the random assignment of that setting fails here: the launch cities were not drawn by lot. Miyaji (2026) calls this design instrumented DiD (DID-IV) and replaces random assignment with two parallel trends assumptions.

::: {.assumption}
**A15 (DID-IV parallel trends; Miyaji Assumptions 5 and 6)** Let $D^{u}_{it}$ be the treatment status of unit $i$ in period $t$ if it were never exposed to the instrument (its city never launched the program). Absent the launch, the two groups' average treatment rates change by the same amount, $\mathbb{E}[D^{u}_{i1} - D^{u}_{i0} \mid S_i = 1] = \mathbb{E}[D^{u}_{i1} - D^{u}_{i0} \mid S_i = 0]$, and the two groups' average "unexposed outcomes" $Y_{it}(D^{u}_{it})$ also change by the same amount.
:::

::: {.theorem}
**Theorem 8 (Wald-DID and the LATET; Miyaji Theorem 1)** Under no carryover, exclusion, monotonicity, no anticipation in the first stage, relevance (Miyaji Assumptions 1 to 4 and 7), and A15, $W_{DID} = \mathbb{E}[Y_{i1}(1) - Y_{i1}(0) \mid S_i = 1, \text{complier}]$, the average effect in the later period for the compliers in the exposed group, whom the instrument pushes into treatment (the LATET).
:::

The proof sketch has three steps, parallel to the LATE theorem of Chapter 18. Step one (exclusion, no carryover, and monotonicity): the exposed group's later-period outcome equals the unexposed outcome plus the compliers' effect times their share. Step two (monotonicity and no anticipation in the first stage): the change in the exposed group's treatment rate equals the change in the unexposed treatment rate plus the compliers' share. Step three (A15): the changes in the two groups' unexposed parts cancel in the DiD, and taking the ratio gives the result.

In words: fuzzy DiD is a kind of IV that swaps the exogeneity of the instrument for parallel trends. The key is the second part of A15: the "unexposed outcome" includes the incumbent treated, so its parallel trends is not just parallel trends in $Y(0)$ but also implicitly requires parallel changes in the incumbent effect.

In Tessera the compliers in the launch cities are exactly the switchers, and the LATET equals 15.857. Absent the launch, however, the bonus increase would raise the control cities' unexposed outcome by an extra 0.30 times 8 orders, but the launch cities' by only 0.08 times 8 orders; A15 fails and Wald-DID is biased downward, the same account as Corollary 1 told another way (Miyaji also notes this in a footnote: when standard parallel trends in $Y(0)$ holds, this assumption in effect restricts effect heterogeneity).

In biostatistics, the instrumented DiD of Ye, Ertefaie, Flory, Hennessy and Small (2023) takes the same route: their identification conditions explicitly include a treatment effect that does not change over time (their Assumption 2(d)), and with the further condition of no unmeasured common effect modifier, the Wald ratio identifies the average treatment effect (their Proposition 1); replacing the latter condition with monotonicity, it identifies the effect for compliers.

::: {.intuition}
The threads of this section share one lesson: borrow a ruler across treatment states, and parallel trends quietly turns into an assumption that effects do not change over time. Wald-DID uses the whole control group as its ruler, and the ruler is contaminated by the rise in the incumbent effect (A4). Using stayers with a different baseline dose as the ruler requires assuming that the average effect does not change over time. If $DID_\ell$ used cities with a different period-1 tier as its ruler, it would rule out lagged and time-varying effects. The repair is also the same: look for the ruler only among units "with the same baseline treatment status and whose own treatment has not changed."
:::

To sum up this section: beyond parallel trends, Wald-DID requires effects to be stable over time (A4), and when the control group's treatment rate moves, it also requires equal switcher effects in the two groups (A5). Wald-TC and Wald-CIC replace A4 with parallel trends by baseline treatment status or with the CIC model, but they need a control group with a constant treatment rate; otherwise only partial identification is possible. With a continuous dose, ordinary parallel trends identifies $\mathrm{ATT}(d \mid d)$, and comparing doses requires strong parallel trends; stayers and quasi-untreated cities serve as the ruler with and without zero-dose units, respectively. With switching treatments, $DID_M$ applies only when there are no dynamic effects, while $DID_\ell$ uses as its ruler the units with the same period-1 tier that have not yet switched and identifies path effects relative to the status quo.

## 4 Estimation

### 4.1 Sample Analogues and Inference for Fuzzy DiD

The estimators of the three Wald ratios are all sample analogues: cell means, and empirical distribution functions and empirical quantiles within cells (de Chaisemartin and D'Haultfœuille 2018 §4). Under i.i.d. sampling and A1, A2, and A6, plus their respective identifying assumptions, all three are $\sqrt{n}$-consistent and asymptotically normal, and the bootstrap is consistent (dCDH 2018 Theorem 7): Wald-DID needs A3, A4, and $\mathbb{E}(Y^2) < \infty$; Wald-TC needs A3' and $\mathbb{E}(Y^2) < \infty$; and Wald-CIC needs A7, A8, and regularity conditions under which the outcome has bounded support and each cell's distribution function is continuously differentiable with a positive derivative. The influence functions are complicated, so in practice one uses the bootstrap. When treatment is assigned at the city level and outcomes are correlated within cities, run a cluster bootstrap by city, drawing cities stratified by group so that all four cells are nonempty.

With more than two groups, the common practice is a regression or 2SLS with group and period fixed effects, which estimates a weighted average of several Wald-DIDs (de Chaisemartin and D'Haultfœuille 2018 §1). Theorem S1 in the Web Appendix of dCDH (2020) shows that a TWFE regression in a fuzzy design is a weighted sum of the effects for the treated in each group and period, and some of the weights can be negative. The paper recommends splitting the groups into three "supergroups" whose treatment rate rises, stays constant, or falls, computing Wald ratios separately with the constant-rate supergroup as the control, and then weighting them (its Theorem 5). Stata has the authors' fuzzydid command; R has a third-party port, Rfuzzydid. When the DiD of the first stage is close to zero, the ratio blows up, so $DID_D$ and its interval should be reported together with the ratio.

### 4.2 Estimation with Doses and Stayers

Estimation for Theorem 5 (i) is direct: for the treated cities, run a flexible regression of $\Delta Y$ on the dose (contdid uses a B-spline, which with cubic degree and no interior knots is a cubic polynomial in the dose), and subtract the average $\Delta Y$ of the zero-dose cities to obtain the $\widehat{\mathrm{ATT}}(d \mid d)$ curve. Differentiating with respect to the dose gives the slope curve, and averaging over the treated cities gives the overall ATT and the overall "ACRT." The overall ATT equals the binary DiD that collapses the dose to 0/1. When the number of knots grows with the sample, this curve is the series estimator of Chapter 14, and the specification in this chapter is the special case of fixed order. Confidence bands use the multiplier bootstrap and can be pointwise or uniform.

The estimators of AS and WAS are doubly robust (dCDH-PSV Theorems 1 and 4) and are implemented in the R package DIDmultiplegtSTAT. The HAD estimator is implemented in DIDHAD: it estimates $\lim_{d \downarrow 0} \mathbb{E}[\Delta Y \mid D_2 \leq d]$ by local linear regression at zero dose, with an MSE-optimal bandwidth and the bias-corrected robust confidence interval of Calonico, Cattaneo and Farrell (2018); the rate of convergence is $n^{-2/5}$, with $n$ the number of cities, slower than the parametric rate. Whether quasi-untreated cities really exist can be tested: the statistic is the smallest dose divided by the gap between the smallest and the second smallest doses, and the test rejects when the statistic exceeds $1/\alpha - 1$, which gives a test of nominal level $\alpha$ (HAD paper Theorem 2; here $\alpha$ is the level of the test, not the weight in Theorem 1).

### 4.3 Negative Weights in TWFE Dose Regressions

In practice the most common approach is still a TWFE regression of the outcome on the dose. Section 4.3 of Chapter 19 showed that under binary staggered adoption the TWFE coefficient is a weighted sum of treated-cell effects; this is where that result does its main work, because doses and switching make the weighting problem worse.

Theorem 1 of dCDH (2020) says that under potential outcomes that depend only on current treatment, common trends, and strong exogeneity, $\beta_{TWFE} = \mathbb{E}\big[\sum_{(i,t): D_{it} = 1} w^{\mathrm{fe}}_{it}\, \tau_{it}\big]$, where $\tau_{it} = Y_{it}(1) - Y_{it}(0)$ is the cell effect and the weight $w^{\mathrm{fe}}_{it}$ is proportional to the residual from regressing $D_{it}$ on unit and period fixed effects, $\varepsilon_{it}$; the weights sum to one and can be negative. The residual is small or even negative in periods when treatment is widespread and for units treated for a long time (its Proposition 1). When treatment is ordered and multivalued, Theorem S3 of its Web Appendix gives the same structure: the TWFE coefficient is a weighted sum of each treated cell's "effect per unit of dose," with weights proportional to $\varepsilon_{it} D_{it}$. The natural benchmark is dose weighting, that is, the sum of effects divided by the sum of doses, which is what the true value of 5.278 in block C measures.

Corollary 1 of the same paper also gives a robustness measure: the smallest standard deviation of the treatment effects under which the TWFE coefficient could coexist with an ATT of zero, which equals the absolute value of the coefficient divided by the standard deviation of the weights; the smaller this number, the more fragile the conclusion is to heterogeneous effects.

In a two-period dose design this structure has a clean closed form.

::: {.theorem}
**Theorem 9 (TWFE in a two-period dose design; CGBS Theorem 3.4(c), HAD paper equation (8))** Suppose everyone has zero dose in period 1 and dose $D \geq 0$ in period 2, no anticipation holds, the untreated trends of all dose groups are the same (A9), and $\mathrm{Var}(D) > 0$. With the convention that the term inside the expectation is zero when $D = 0$, the probability limit of the TWFE dose coefficient is

$$\beta_{TWFE} = \mathbb{E}\Big[\frac{(D - \mathbb{E}[D])\, D}{\mathbb{E}[(D - \mathbb{E}[D])\, D]} \cdot \frac{\mathrm{ATT}(D \mid D)}{D}\Big]$$

Zero-dose units get zero weight, the weights sum to one, and treated units with doses below the full-sample average dose $\mathbb{E}[D]$ get negative weights.
:::

The proof sketch has two steps. Step one (algebra): with two periods and unit fixed effects, TWFE is the OLS regression of $\Delta Y$ on $D$, $\beta_{TWFE} = \mathrm{Cov}(\Delta Y, D) / \mathrm{Var}(D) = \mathbb{E}[(D - \mathbb{E}[D])\, \Delta Y] / \mathbb{E}[(D - \mathbb{E}[D])\, D]$. Step two (no anticipation and A9): $\mathbb{E}[\Delta Y \mid D] = \mu_0 + \mathrm{ATT}(D \mid D)$, where the constant $\mu_0$ is the average change of the zero-dose cities, and the expectation of its product with $D - \mathbb{E}[D]$ is zero; substituting gives the result.

In words: TWFE estimates a weighted sum of the per-unit-dose effects $\mathrm{ATT}(d \mid d)/d$, with weights determined mechanically by the dose distribution rather than chosen by the researcher. It treats units above the average dose as the effective treated group and units below the average dose (some of which are treated) as the effective control group. When the per-unit-dose effect does not vary across doses, any such weighting recovers the same number and TWFE equals WAS; the negative weights cause real bias only when the per-unit-dose effect varies with the dose, or when the platform picks doses by expected gains.

The HAD paper draws a testable implication from this: if $\mathbb{E}[\Delta Y \mid D_2]$ is not a linear function of the dose, TWFE should not be used (its Theorem 3). The weights also change with the dose distribution itself: in their Medicare application, CGBS find that dropping the zero-dose group shrinks the TWFE coefficient by 78%. In Tessera's block B, 125 treated cities have negative weights, all of them low-dose cities (Section 6.3), yet TWFE exceeds WAS by only 0.114, because the per-unit-dose effect rises only mildly.

When treatment can switch and there are lagged effects, another layer of problems appears: the decomposition of dCDH (2020) assumes that potential outcomes depend only on current treatment, so cells that are paused but still carry residual effects enter the regression as "controls" and pull the coefficient down. The event-study version does not escape either: a TWFE event study with dose-by-period dummies puts negative weights on groups below the average dose (dCDH 2026 Proposition 1), and a local-projection panel regression under Design 4 (zero dose for everyone in period 1) can be biased even with homogeneous effects, with weights that can sum to far less than one (its Proposition 2).

Putting several treatment variables in the regression at once adds yet another layer: de Chaisemartin and D'Haultfœuille (2023b) prove that under parallel trends each treatment's coefficient equals a weighted sum of that treatment's effects (with possibly negative weights) plus a weighted sum of the other treatments' effects, so one treatment's coefficient is contaminated by the others, and dropping a treatment from the regression may even reduce the bias. The R package TwoWayFEWeights computes these weights for a single treatment.

Section 6 will show that the TWFE in block B has 125 negative weights yet is close to WAS, while the TWFE in block C has a smaller share of negative weights yet is 13.4 standard errors from the truth. The number of negative weights is not a measure of bias; what matters is how the weights correlate with the effects.

### 4.4 Estimation of $DID_\ell$, Placebos, and Inference

$DID_\ell$ and the normalized $DID^n_\ell$ are computed by did_multiplegt_dyn in the R package DIDmultiplegtDYN; the average total effect is reported in the output as the average total effect per unit of treatment, and standard errors are clustered at the city level. The placebo estimators compare the changes of switchers and their control sets over the periods before the switch, testing an observable implication of A14; report both the individual placebos and the joint test. For cities whose treatment first rises and then falls, crossing the status quo, the effect carries no sign guarantee; the paper's approach is to drop these cells (the drop_larger_lower option), which yields an unbalanced panel that satisfies its design restriction.

::: {.warning}
The average total effect is a cumulative quantity, not a per-period effect: a dose's effect is accumulated over several periods on average, the software also reports this average number of periods, and only the total effect divided by it is the average effect of one period of exposure (de Chaisemartin and D'Haultfœuille 2026 §3.3). Reading it as "how many more orders per tier per week" counts several weeks of cumulative effect as one week.
:::

### 4.5 The Estimators at a Glance

| Estimator | Treatment type | Ruler (control) | Key assumptions | Target | R implementation |
|---|---|---|---|---|---|
| Wald-DID | Binary, fuzzy | The whole control group | A3, A4, and A5 or A6 | LATE of treated-group switchers | By hand or 2SLS |
| Wald-TC | Binary, fuzzy | Control-group units with the same baseline treatment status | A3', A6 | Same as above | By hand, Rfuzzydid |
| Wald-CIC | Binary, fuzzy | Same as above, matched by rank | A7, A8, A6 | Same as above, plus LQTE | By hand, Rfuzzydid |
| contdid | Continuous dose, with zero dose | Zero-dose units | A9; slopes also need A10 | ATT(d\|d); ACRT(d) under A10 | contdid |
| AS / WAS | Continuous dose, with stayers | Stayers with the same baseline dose | A11; AS also needs no quasi-stayers | AS, WAS | DIDmultiplegtSTAT |
| HAD | Continuous dose, no zero dose | Quasi-untreated units | A12, bounded effects | WAS | DIDHAD |
| $DID_\ell$ | Non-binary, switching | Not-yet-switchers with the same period-1 tier | A13, A14 | AVSQ and normalized effects | DIDmultiplegtDYN |
| TWFE | Any | Implicit, may include treated cells | Homogeneous effects, no dynamics | A weighted sum with no clear interpretation | fixest; weights via TwoWayFEWeights |

To sum up the route of this section: the three fuzzy ratios are sample analogues of cell means or quantiles, with inference by cluster bootstrap. In dose designs, first estimate the $\mathrm{ATT}(d \mid d)$ curve; slopes additionally rely on strong parallel trends, and AS, WAS, and HAD each have dedicated estimators. TWFE weights are proportional to the treatment residual times the dose; in a two-period dose design the cells below the average dose get negative weights, and switching and lagged effects also let paused cells slip into the control group. Switching treatments call for $DID_\ell$, together with placebos and an interpretation of the normalized effects.

## 5 Anchoring Papers

The three anchoring papers correspond to the three threads of this chapter: a classic re-estimation with fuzzy DiD, a methods application with switching non-binary treatments, and an IS study that describes itself as following a dCDH design. Each is examined under five headings (paper, method, data, results, limitations).

### 5.1 de Chaisemartin and D'Haultfœuille (2018) Revisit Duflo (2001)

::: {.case}
The paper: "Fuzzy Differences-in-Differences," Review of Economic Studies. Its application revisits the study by Duflo (2001, American Economic Review) of the INPRES school construction program in Indonesia.

Method: Duflo takes people born from 1968 to 1972, who entered school after the program began, as cohort 1 and those born from 1957 to 1962 as cohort 0, splits districts into treated and control groups by the residual from a regression of school construction intensity, and estimates the returns to education with Wald-DID and with 2SLS including cohort and district fixed effects. dCDH first report Wald-DID, 2SLS, and the TC and CIC bounds on the original grouping, and then regroup along the supergroup idea of the paper's Theorem 5: for each district they test whether the education distributions of the two cohorts are the same, and districts with p-values above 0.5 go into the control group. Treatment here is years of schooling, not 0/1: dCDH extend the three ratios to ordered treatments (de Chaisemartin and D'Haultfœuille 2018 §3.3), and the target changes from the LATE of switchers to the average causal response to one more year of schooling (ACR; not the ACRT of Section 2.2).

Data: the 1995 Indonesian intercensal survey used by Duflo. dCDH could not fully reconstruct her 31,061 observations and use a sample of 30,828, so their Wald-DID and 2SLS differ slightly from the original; standard errors are clustered by district.

Results: on the original grouping, Wald-DID is 0.195 with interval [-0.102, 0.491], and 2SLS is 0.073. But schooling also rose by 0.36 years in the control group, against 0.47 years in the treated group, so the first stage is only 0.11; Wald-DID then equals a weighted difference of the two groups' ACRs and can be far from either group's return when the returns differ, and 2SLS likewise requires the returns to be the same in every district. The TC bounds [-3.70, 2.18] and CIC bounds [-5.60, 3.36] are too wide to be informative. After regrouping, schooling rises by 0.99 years in 123 districts and falls by 0.74 years in 97, while it barely changes (-0.05 years) in the 64 control districts. Wald-DID is now 0.140 (0.015), Wald-TC 0.101 (0.017), and Wald-CIC 0.099 (0.017); the differences between the latter two and Wald-DID are statistically significant, with t-statistics of -4.27 and -4.61. The authors prefer Wald-TC and Wald-CIC for exactly the reason of A4 in Section 3.2: diminishing returns to experience mean that the effect of education is not stable across the two cohorts.

Limitations: the control group is estimated from the data, and the sampling error of the first-step grouping does not fully enter the second-step inference; the authors acknowledge that this part is tentative and run robustness checks with different p-value thresholds and sample splits. With cohorts in place of time, A1 and A3' must be reinterpreted as "comparability between two birth cohorts in the same district."
:::

The value of this paper is that it takes apart a widely cited number: the change in conclusions comes from the choice of control group and the stance taken on A4, not from computational tricks.

### 5.2 de Chaisemartin and D'Haultfœuille (2026) Revisit Favara and Imbs (2015)

::: {.case}
The paper: "Difference-in-Differences Estimators of Intertemporal Treatment Effects," Review of Economics and Statistics. Its application revisits the study by Favara and Imbs (2015, American Economic Review) of banking deregulation and the housing market.

Method: after the interstate banking law of 1994, states gradually lifted four types of restrictions. Treatment is the number of restrictions a state has lifted, taking values from 0 to 4; it is non-binary, and some states deregulated more than once. The original paper estimates dynamic effects with local-projection panel regressions; dCDH re-estimate them using current treatment and without control variables, decompose the weights of these regressions with their Proposition 2, and finally re-estimate the effects with $DID_\ell$ and $DID^n_\ell$. In 1994 every state's treatment is 0, so the event-study parameters do not rely on an assumption that the effects of different policies are homogeneous.

Data: a county-level panel for 1994 to 2005; the outcomes are the log growth rates of mortgage originations and of house prices. Of the states, 8 never deregulate, 33 deregulate once, 8 twice, and 1 three times; the effects at $\ell = 1$ and $\ell = 8$ cover 905 and 773 counties, respectively.

Results: the first local-projection coefficient is a weighted sum of 7,626 effects, of which 4,670 have positive weights and 2,956 negative weights; the positive and negative weights sum to 1.067 and -0.125, respectively. By the fourth horizon the weights sum to -0.018, so even with homogeneous effects the coefficient could have the opposite sign from the true effect. $DID_\ell$ shows the effect on loan growth building up gradually: $DID_1$ is 0.043 (0.035), $DID_3$ is 0.081 (0.049), and $DID_5$ is 0.148 (0.064); the joint placebo test has a p-value of 0.400. A test of whether $DID^n_\ell$ is constant in $\ell$ has a p-value of 0.328 and does not reject, contrary to the original paper's conclusion that the effect is only temporary. For house prices, $DID_1$ is 0.003 (0.004) and $DID_5$ is 0.026 (0.010), with a joint placebo p-value of 0.703.

Limitations: the placebos can test only three years, while unbiasedness of $DID_8$ requires parallel trends to hold for eight years; citing Roth (2022), the authors warn that insignificant placebos do not rule out trend differences large enough to cause bias. Four types of restrictions are collapsed into one count, so the effect of "how many were lifted" is still, in interpretation, a mixture of different policies.
:::

This paper shows how severe the weighting problem of Section 4.3 is in real data, and how $DID_\ell$ rewrites a story of "short-lived effects" into one of "persistent effects."

### 5.3 Tunc, Cavusoglu and Zheng (2025, ISR)

::: {.case}
The paper: "Resale Royalty in Non-Fungible Token Marketplaces: Blessing or Burden for Creators and Platforms?" Information Systems Research. The question is how resale royalties in NFT marketplaces affect creators' pricing and market liquidity.

Method: on October 18, 2021, the NFT marketplace Rarible introduced lazy minting, under which creators no longer pay on-chain fees up front when minting an NFT. The authors take creators who moved from traditional minting to lazy minting after the policy as the treated group and creators who kept using traditional minting as the control group, and state that their design follows the $DID_M$ of dCDH (2020) because it allows heterogeneous treatment effects as well as staggered adoption. Robustness checks include a two-period doubly robust DiD (Sant'Anna and Zhao 2020), DiD after coarsened exact matching, and event-study-style parallel trends tests.

Data: minting records on Rarible from May 2020 to May 2023. To avoid changes in composition, the DiD analysis keeps only creators who were already minting before the policy: 2,446 in the treated group and 953 in the control group.

Results: lazy minting lowered the royalty rates that creators set. The full-sample fixed-effects regression coefficient is -0.042 (0.008, log royalty rate); the treated-versus-control DiD is -0.16 (0.04), and -0.248 (0.116) after adding collection fixed effects; the two-period doubly robust DiD is -0.103 (0.019), and the DiD after matching is -0.093 (0.045).

Limitations: the policy took effect for all creators at the same time, and the treated and control groups are defined by whether creators themselves adopted. In the language of this chapter, there is no control group "whose exposure did not change"; the comparison is between self-selected adopters and non-adopters, and identification rests entirely on parallel trends between these two types of people. The authors acknowledge that self-selection cannot be avoided. Point estimates range from -0.248 to -0.093 across DiD specifications, a sizable difference in magnitude, which suggests that the conclusion is sensitive to how the comparison set is defined.
:::

Taken together: dCDH (2018) show how the stability of the control group's treatment rate determines whether point identification is possible; dCDH (2026) show how bad the weights of TWFE-type regressions can get under non-binary, switching treatments; and the study by Tunc and coauthors shows that IS research already invokes dCDH designs, and also that when a policy applies universally and individuals select into treatment, the fuzzy structure is still there, merely hidden by "grouping by adoption."

## 6 A Full Walkthrough on Tessera Data

The code below uses R 4.6.1. Every number quoted in the text comes from actual output. The listings below show the main calls; auxiliary computations such as the city bootstrap checks and the TWFE weight decompositions are not listed line by line.

### 6.1 The DGP

The key settings of block A are as follows (in data.table syntax; inside the function PA is called P). Enrollment is determined by a latent index with both a gain component and a level component; in the launch cities, thresholds separate the beta testers from the switchers, and in the control cities the incumbent enrollees are the same in both periods. At $t = 1$ every enrollee's effect rises by lambda.

```r
set.seed(2001)
PA <- list(
  n_city     = 80,     # cities per group
  mu0        = 70,     # mean weekly orders without the program
  gamma1     = 4,      # common seasonal rise from t = 0 to t = 1
  mu_beta    = 6,      # courier gain beta_i ~ N(mu_beta, sd_beta^2) at the t = 0 bonus
  sd_beta    = 5,
  lambda     = 8,      # rise in every enrolled courier's effect at t = 1
  pA1        = 0.08,   # g = 1: always-takers (beta testers)
  pS1        = 0.46,   # g = 1: switchers (enroll when the program opens)
  pA0        = 0.30,   # g = 0: incumbent enrollees, the same couriers at t = 0 and t = 1
  rho_beta   = 0.7,    # latent index V = rho_beta z_beta + kappa_a z_a + s xi + v_city
  kappa_a    = 0.3,    #   (z_beta, z_a standardized gain and level; s = sqrt(1 - rho^2 - kappa^2))
  pS0_rising = 0.10    # variant A': control-group switchers (share of all control couriers)
  # ... city and courier noise parameters omitted here
)
cr[, V := P$rho_beta * z_beta + P$kappa_a * z_a + s_xi * rnorm(N) + v_c]
cr[, type_true := ifelse(g == 1,
                         ifelse(V >= thr_A1, "always", ifelse(V >= thr_S1, "switcher", "never")),
                         ifelse(V >= thr_A0, "always", "never"))]
pan[, y0 := P$mu0 + c_lvl + a + P$gamma1 * t + ifelse(t == 1, eta1, eta0) + rnorm(.N, 0, P$sd_eps)]
pan[, tau_true := beta + P$lambda * t]
x[, d := as.integer(type_true == "always" | (type_true == "switcher" & t == 1))]
x[, y := y0 + d * tau_true]
```

In variant A′, the non-enrolled couriers in the control cities with the highest order-level index enroll at $t = 1$; they make up 0.10 of the control couriers. Block B's doses and gains are generated as follows, with the gain coefficient rising in the dose, so that A9 holds and A10 fails; the dose_H line generates the no-stayer variant (hereafter B′):

```r
set.seed(2002)
uB[, stayer := rank(s) <= round(PB$share_zero * PB$N)]           # 130 lowest s_i get no bonus
uB[, dose_B := ifelse(stayer, 0, runif(.N, PB$d_lo, PB$d_hi))]    # USD per order, U(0.2, 3.0)
uB[, b := ifelse(stayer, 1, Eb_B(dose_B)) * exp(PB$sd_nu * rnorm(.N) - PB$sd_nu^2 / 2)]
uB[, dose_H := ifelse(stayer, runif(.N, PB$had_lo, PB$had_hi), dose_B)]   # no-stayer variant
x[, tau_true := post * b * m_B(dose)]                             # m(d) = 4 d - 0.6 d^2
x[, y := alpha + PB$gamma2 * post + ifelse(t == 1, eps1, eps2) + tau_true]
```

Block C's tiers follow a weekly Markov chain; effects build up over four weeks and linger after a pause, and cities that switch earlier gain more:

```r
PC <- list(
  N          = 200,     # cities
  nT         = 30,      # weeks
  # ... share of never-switchers, first-switch window and level parameters omitted here
  trans      = rbind(c(0.70, 0.25, 0.05),    # from 0 (paused): stay, resume at 1, resume at 2
                     c(0.10, 0.80, 0.10),    # from 1: pause, stay, raise
                     c(0.08, 0.12, 0.80)),   # from 2: pause, lower, stay
  f          = c(0, 1, 1.6),   # f(0), f(1), f(2)
  omega      = c(1.0, 0.6, 0.4, 0.2)
)
set.seed(2003)
for (k in seq_along(PC$omega)) {
  lag <- k - 1
  cols <- (lag + 1):PC$nT
  tau_m[, cols] <- tau_m[, cols] + PC$omega[k] * fD[, cols - lag, drop = FALSE]
}
tau_m <- tau_m * uC$b
```

In block C, 50 cities never pay a bonus, and the rest switch on for the first time between weeks 6 and 20, with a probability of 0.4 of going straight to tier 2 at the first switch-on; after that the tier can be paused, resumed, raised, or lowered, and 71.5% of the cities switch at least twice.

### 6.2 Block A: The Three Wald Ratios

All three estimators are hand-written functions that follow the formulas of dCDH (2018); their input is courier-level data with group, period, and treatment columns, and inference uses a city bootstrap stratified by group, with 999 replications.

```r
est3 <- function(X) c(did = wald_did(X, "y", "d", "g", "t")$estimate,
                      tc  = wald_tc (X, "y", "d", "g", "t")$estimate,
                      cic = wald_cic(X, "y", "d", "g", "t")$estimate)
set.seed(20201)
for (b in seq_len(999)) {
  drawn <- unlist(lapply(cities_by_g, function(cc) sample(cc, length(cc), replace = TRUE)), use.names = FALSE)
  rows  <- unlist(rows_by_city[drawn], use.names = FALSE)
  outA[b, ]  <- est3(A[rows])
  outAR[b, ] <- est3(AR[rows])
}
```

The core of wald_tc is the formula of Section 3.3: first compute $\delta_0$ and $\delta_1$ by treatment status in the control group, then add the corresponding trend to each person in the launch group's base period. The function wald_cic replaces this step with the rank mapping $Q_d$, computing quantiles exactly with the generalized inverse defined in the paper.

| Estimator | Block A (constant control enrollment rate) | 95% interval | Variant A′ (rising control enrollment rate) | 95% interval |
|---|---|---|---|---|
| True $\mathrm{LATE}_{sw}$ | 15.857 | | 15.857 | |
| Wald-DID | 12.002 | [11.290, 12.774] | 12.124 | [11.225, 13.149] |
| Wald-TC | 15.857 | [15.168, 16.578] | 21.349 | [20.537, 22.242] |
| Wald-CIC | 15.947 | [15.250, 16.661] | 22.283 | [21.357, 23.391] |

![The three Wald estimators and their 95% city bootstrap intervals in block A and variant A′; the dashed line is the true LATE of the treated-group switchers.](assets/fig/fig_20_fuzzy.svg)

In block A the intervals of Wald-TC and Wald-CIC both cover the truth, while the Wald-DID interval lies far below it, exactly as Theorems 1 through 3 and Corollary 1 predict. In variant A′ all three estimators fail: Wald-DID is still biased downward because A4 fails, and Wald-TC and Wald-CIC are biased sharply upward because A6 fails. Note that in A′ Wald-DID ends up closer to the truth, which does not make it more reliable: its two biases go in opposite directions and partly cancel, one from the failure of A4 and one from the difference between the two groups' switcher effects when $\alpha > 1$ (Theorem 1 (i)).

The bootstrap standard errors of the three are similar (0.389, 0.364, and 0.364 in block A), so Wald-TC pays no precision cost for relaxing A4. The third-party R port Rfuzzydid gives the same Wald-DID and Wald-TC as the hand-written functions; its Wald-CIC differs slightly from the paper's generalized inverse in its rounding convention for quantiles, with differences in the fourth decimal place, and the text relies on the hand-written implementation.

### 6.3 Block B: The Dose Curve, AS, WAS, and HAD

```r
df_B <- as.data.frame(B); dfH <- as.data.frame(H)
d_pos <- B[post == 1 & dose > 0, dose]
dv <- sort(unique(c(min(d_pos), seq(0.25, 2.95, by = 0.05), max(d_pos))))   # spans the treated doses exactly
set.seed(20301)
rc <- cont_did(yname = "y", dname = "dose", gname = "g", tname = "t", idname = "unit", data = df_B,
               target_parameter = "slope", aggregation = "dose", treatment_type = "continuous",
               dvals = dv, degree = 3, num_knots = 0,
               control_group = "nevertreated", biters = 9999, cband = FALSE)

B[, treat_post := as.integer(post == 1 & dose > 0)]
m_bin  <- feols(y ~ treat_post | unit + t, data = B, cluster = ~unit)
m_twfe <- feols(y ~ d | unit + t, data = B, cluster = ~unit)
tw <- twowayfeweights(df_B[, c("unit", "t", "y", "d")], Y = "y", G = "unit", T = "t", D = "d",
                      type = "feTR", summary_measures = TRUE)
st <- did_multiplegt_stat(df_B[, c("unit", "t", "y", "d")], Y = "y", ID = "unit", Time = "t", D = "d",
                          estimator = c("aoss", "waoss"), aoss_vs_waoss = TRUE)
hd <- did_had(dfH, outcome = "y", group = "unit", time = "t", treatment = "d",
              effects = 1, placebo = 0, yatchew = TRUE, graph_off = TRUE)
```

::: {.warning}
Two practical notes on contdid 0.1.1. First, the evaluation grid for the curve must run exactly from the smallest to the largest treated dose; otherwise the curve is built on a different set of basis functions and comes out wrong (the overall ATT and the overall slope are unaffected), which is why the grid above is constructed by hand. Second, the package's standard error for the overall "ACRT" is too narrow (0.228, against 0.295 from the city bootstrap), so for the overall "ACRT" we report the city bootstrap interval below; the package's standard error for the overall ATT agrees with the city bootstrap (0.261 versus 0.257), and we use the package's interval as is.
:::

The overall ATT is 4.690, with package interval [4.178, 5.202], covering the true value of 4.858; the binary DiD gives the same point estimate. The curve's pointwise coverage over the 57 dose points is 0.95; at $d = 1.5$ the estimate is 4.508 with interval [3.972, 5.043], against a true value of 4.51. The slope is another matter: the overall "ACRT" is 3.526 with city bootstrap interval [2.931, 4.070], which covers the slope identifiable under parallel trends, 3.186, but is far from the true $\mathrm{ACRT}^{loc}$ of 1.854. The pointwise coverage of the slope curve relative to the true $\mathrm{ACRT}(d \mid d)$ is only 0.11.

![Left: the ATT(d | d) curve estimated by contdid with 95% pointwise intervals, against the truth. Right: the slope of the fitted curve, against the slope identifiable under parallel trends (solid line) and the true causal response ACRT(d | d) (dashed line).](assets/fig/fig_20_dose.svg)

The right panel illustrates Theorem 5 (ii): the estimated slope tracks the gray solid line, and the distance between the two true curves is the selection bias. The stayer estimators give AOSS (that is, AS) of 2.586 with interval [2.132, 3.040] and WAOSS (that is, WAS) of 2.932 with interval [2.659, 3.206], covering the true values of 2.923 and 3.038, respectively; their difference is $-0.346$ (0.128), with a p-value of 0.0068, which says that the effect per USD really does rise with the dose. In block B everyone has zero dose in period 1, so WAOSS reduces exactly to the overall ATT divided by the treated cities' average dose (dividing the full-precision values of 4.690 and 1.599 gives exactly 2.932), and the two packages corroborate each other on the same quantity.

On the no-stayer variant, HAD estimates WAS at 2.643 with interval [1.922, 3.404], covering the true value of 3.059; the quasi-untreated test statistic is 3.100 with a p-value of 0.244, so the test does not reject that quasi-untreated cities exist; the local linear bandwidth is 0.574, with 181 cities inside it; and the interval is much wider than those of the stayer estimators, the price of a nonparametric rate.

The TWFE dose regression gives 3.152 with interval [2.919, 3.385]. TwoWayFEWeights shows that of the 370 treated cities, 245 have positive weights and 125 negative ones, with positive and negative weights summing to 1.066 and $-0.066$. The negative weights all fall on low-dose cities: the highest dose among cities with negative weights is \$1.174, and the lowest among cities with positive weights is \$1.202, so the dividing line is the full-sample average dose, consistent with the closed form of Section 4.3. Even so, TWFE is not far from WAS, only 0.114 higher, and its interval covers both WAS and AS: in this DGP the per-unit-dose effect rises only mildly with the dose, and the negative weights fall on low-dose cities whose effect per USD is slightly lower, so the coefficient is slightly above WAS and the bias is limited.

On the B′ data the Yatchew linearity test has a p-value of 0.189 and does not reject linearity. Yet the B′ DGP is not linear, since the per-unit-dose effect varies with the dose, and when quasi-untreated cities exist, linearity is exactly equivalent to homogeneous linear effects (the second part of Theorem 3 in the HAD paper). Failing to reject only says that the test lacks power with 500 cities; it is not evidence that linearity holds. What TWFE really misses is $\mathrm{ACRT}^{loc}$: reading 3.152 as "the marginal effect per USD of bonus" overstates the truth by 1.298.

### 6.4 Block C: TWFE and $DID_\ell$

```r
df <- as.data.frame(C[, .(unit, t, y, dose)])
m_twfe <- feols(y ~ dose | unit + t, data = C, cluster = ~unit)
tw <- twowayfeweights(df, Y = "y", G = "unit", T = "t", D = "dose", type = "feTR", summary_measures = TRUE)

library(polars)            # must be attached before DIDmultiplegtDYN
library(DIDmultiplegtDYN)
dyn  <- did_multiplegt_dyn(df, outcome = "y", group = "unit", time = "t", treatment = "dose",
                           effects = 8, placebo = 4, cluster = "unit", normalized = FALSE, graph_off = TRUE)
dynn <- did_multiplegt_dyn(df, outcome = "y", group = "unit", time = "t", treatment = "dose",
                           effects = 8, placebo = 4, cluster = "unit", normalized = TRUE, graph_off = TRUE)
```

Static TWFE gives 3.847 per tier (SE 0.107) with interval [3.636, 4.057], only 72.9% of the dose-weighted true effect of 5.278 and 13.4 standard errors away from it. Rerunning TWFE with the outcome replaced by the noise-free true effect gives 3.869, so the bias comes from the structure of the weights, not from sampling. This 3.869 can be decomposed exactly: the sum of "weight times per-tier effect" over treated cells is 4.388, and cells that are paused but still carry residual effects contribute $-0.519$ in their role as controls.

As for the weights themselves: of the 2,073 treated cells, 1,721 have positive and 352 negative weights; the share of negative weights is 0.170 by count but only 0.026 by absolute weight, and the negatively weighted cells contribute $-0.212$ in total. The bigger problem is how the weights are allocated: tier 2 cells take 0.868 of the total weight, yet their average effect per tier is only 4.408, against 6.214 for tier 1 cells; concave tier effects combined with weights tilted toward tier 2 drag the coefficient down.

| $\ell$ | $DID_\ell$ | 95% interval | True | $DID^n_\ell$ | 95% interval | True |
|---|---|---|---|---|---|---|
| 1 | 3.628 | [3.120, 4.137] | 3.547 | 2.604 | [2.240, 2.969] | 2.546 |
| 2 | 5.426 | [4.797, 6.054] | 5.366 | 2.030 | [1.794, 2.265] | 2.007 |
| 3 | 6.327 | [5.597, 7.057] | 6.390 | 1.628 | [1.440, 1.816] | 1.644 |
| 4 | 6.869 | [6.136, 7.601] | 6.672 | 1.376 | [1.229, 1.522] | 1.336 |
| 5 | 6.419 | [5.644, 7.195] | 6.369 | 1.055 | [0.927, 1.182] | 1.046 |
| 6 | 6.315 | [5.555, 7.075] | 6.243 | 0.878 | [0.772, 0.984] | 0.868 |
| 7 | 6.254 | [5.444, 7.064] | 6.147 | 0.757 | [0.658, 0.855] | 0.744 |
| 8 | 6.141 | [5.370, 6.912] | 5.949 | 0.663 | [0.579, 0.746] | 0.642 |

All sixteen intervals cover the truth, and every $\ell$ is identified from all 150 switching cities. The shape of $DID_\ell$ tells the DGP's story: the effect builds up over the first four weeks and then dips slightly because of pauses and reductions. The normalized effects decline monotonically, which is exactly the weighting implication of Theorem 7 (ii): $DID^n_\ell$ is a weighted average of the effects of current treatment and its first $\ell - 1$ lags, and in the DGP the lagged effects decline with the lag (the coefficients on successive lags are 1, 0.6, 0.4, and 0.2) and vanish beyond three periods, so the larger $\ell$ is, the more weight falls on distant lags with small or even zero effects.

The four placebos are 0.261, 0.245, 0.547, and 0.215, with a joint test p-value of 0.284. To be candid, the third placebo on its own has interval [0.058, 1.036], which excludes zero, although the DGP has no anticipation effect at all: this is sampling fluctuation in which the four placebos, which share a base period, drift upward together. The joint test should therefore take precedence; conversely, an insignificant joint test should not be taken as proof of parallel trends either. The joint test of the effects has p < 0.0001.

![$DID_\ell$ and placebos from did_multiplegt_dyn with 95% intervals, against the true effects; event time is $e = \ell - 1$, and the base period $e = -1$ is normalized to zero.](assets/fig/fig_20_switching.svg)

The average total effect is 5.113 (SE 0.262) with interval [4.599, 5.626], covering the true value of 5.038. It is a cumulative quantity: a dose's effect accumulates over 4.722 periods on average, and only the total effect divided by this number of periods is the average effect of one period of exposure. The figure below sets static TWFE beside it; the two vertical lines are two different true values, and TWFE is far from both.

![Per-tier effects in block C: static TWFE, static TWFE on the noise-free outcome, and the average total effect from did_multiplegt_dyn, against the true average total effect (dashed line) and the dose-weighted true per-tier effect (dotted line).](assets/fig/fig_20_twfe.svg)

### 6.5 The Grand Reconciliation

| Block | Estimator | Estimate | 95% interval | Target and true value | Covered |
|---|---|---|---|---|---|
| A | Wald-DID | 12.002 | [11.290, 12.774] | LATE 15.857 | No |
| A | Wald-TC | 15.857 | [15.168, 16.578] | LATE 15.857 | Yes |
| A | Wald-CIC | 15.947 | [15.250, 16.661] | LATE 15.857 | Yes |
| A′ | Wald-TC | 21.349 | [20.537, 22.242] | LATE 15.857 | No |
| B | contdid overall ATT | 4.690 | [4.178, 5.202] | ATT 4.858 | Yes |
| B | contdid overall "ACRT" | 3.526 | [2.931, 4.070] | ACRT 1.854 | No |
| B | WAOSS | 2.932 | [2.659, 3.206] | WAS 3.038 | Yes |
| B | TWFE dose regression | 3.152 | [2.919, 3.385] | WAS 3.038 | Yes |
| B′ | HAD | 2.643 | [1.922, 3.404] | WAS 3.059 | Yes |
| C | Static TWFE | 3.847 | [3.636, 4.057] | Per-tier effect 5.278 | No |
| C | $DID_1$ | 3.628 | [3.120, 4.137] | 3.547 | Yes |
| C | Average total effect | 5.113 | [4.599, 5.626] | 5.038 | Yes |

Interval types depend on the estimator: the three Wald ratios and the contdid overall "ACRT" use city bootstrap percentile intervals, HAD uses the bias-corrected robust interval, and the rest are the intervals reported by the respective packages or fixest.

To sum up the reconciliation: in block A, with a stable control-group enrollment rate, Wald-TC and Wald-CIC hit the truth while Wald-DID is too low, and once the control group's enrollment rate moves, all three fail. In block B, ordinary parallel trends reliably identifies the height of the ATT, but the slope follows the selection bias, and in this nearly linear world TWFE is close to WAS and far from the ACRT. In block C, TWFE is off by 13.4 standard errors because of mismatched weights and residual effects during pauses, while $DID_\ell$ and the normalized effects cover the truth throughout.

## 7 Failure Modes of Identification and Robustness

A moving control-group treatment rate is the primary threat to fuzzy DiD. A6 can be tested directly: report the treatment rate in each cell and test the difference between the control group's treatment rates in the two periods. If no naturally stable control group can be found, dCDH (2018) pick out the groups whose treatment rate does not change and use them as controls, erring on the side of declaring that a group's treatment rate moves, because misclassifying a group with a moving treatment rate as a control (a type II error) is costlier than the reverse; the sampling error of the first-step group selection must also be acknowledged. Covariates can help too: if the change in the treatment rate is driven entirely by changes in covariate composition, A6 conditional on covariates may hold. Failing all else, report the TC or CIC bounds and explain why they are wide.

The gap between Wald-DID on the one hand and Wald-TC and Wald-CIC on the other is itself a diagnostic. In data where A6 holds, the three should be close; when Wald-DID departs systematically from the other two, the most likely cause is a changing effect for the incumbent treated. If an earlier stretch with stable treatment rates is available, placebo versions of the three ratios can be computed (de Chaisemartin and D'Haultfœuille 2018 §2.5).

Conversely, Wald-TC is no cure-all: when couriers enroll because their $Y(1)$ improves in the later period, A3' and A4 fail together. Wald-TC depends on the scale of the outcome and, when the two groups' base-period distributions differ a lot, is sensitive to whether the outcome is logged; Wald-CIC is invariant to monotone transformations, but it restricts the whole distribution. In a panel the switchers are visible, and some researchers run a DiD directly between the switchers and the never-takers of their own group. That relies on "the two self-selected types being parallel in $Y(0)$," a different and often stronger assumption, and it is exactly the problem the design of Tunc and coauthors faces. When the DiD of the first stage is close to zero, all three ratios blow up, a problem with the same root as weak instruments in Chapter 18.

In dose designs, the greatest danger is reading the slope as a causal response. Platforms and governments rarely assign doses at random; they tend to give more resources to units that respond more, which is exactly the case in which A10 fails. Even when A10 holds, the slope identifies only the $\mathrm{ACRT}(d)$ of all treated cities, not the marginal effect for the cities that received $d$. Pre-trends tests look only at $Y(0)$ and have no power at all against violations of A10. What one can do is argue from institutional knowledge that dose assignment is unrelated to gains, or simply report only quantities such as $\mathrm{ATT}(d \mid d)$ and WAS that are identified under A9 alone.

Without zero-dose units, run the quasi-untreated test first; when it rejects, the HAD estimator cannot be used, and all one can do is test the linearity condition on which TWFE relies (the Stute or Yatchew test), while acknowledging that only one of its implications has been tested. HAD intervals are inherently wide, and in modest samples HAD may fail to detect even a sizable effect.

With switching treatments, the threat comes from dynamics. Paths on which treatment crosses the status quo carry no sign guarantee for their effects and should be dropped; when there are few kinds of paths, reporting $DID_\ell$ path by path is easier to interpret. Long-horizon $DID_\ell$ requires parallel trends to hold over a longer window, while the window that placebos can test is often shorter. The warning of Roth (2022) applies here as well: insignificant placebos may simply reflect a lack of power. In Tessera the third placebo is significant on its own while the joint test is not, which shows that reading the stars on placebos one at a time can produce both false alarms and misses. The average total effect and the normalized effects are both weighted quantities, and the weights should be spelled out when reporting them: the former weights by dose increments and accumulates over the window, while the weights of the latter tilt toward the lags as $\ell$ grows.

These failure modes come down to one thing: once treatment is no longer tidy, parallel trends must be stated for the right population and the right potential outcome. The ruler must be found among units "with the same baseline treatment status and whose own treatment has not changed." The control group's treatment rate, the mechanism that assigns doses, and the shape of the treatment paths all have to be put on the table before deciding which estimator to use.

## 8 Further Reading

::: {.readings}
Required reading, in suggested order:

- de Chaisemartin and D'Haultfœuille (2018, Review of Economic Studies). The original fuzzy DiD paper; focus on Theorems 1 to 3, on §2.5 for how to choose among the three ratios, and on §3.1 for partial identification.
- Callaway, Goodman-Bacon and Sant'Anna (2025, working paper arXiv v8, forthcoming in the American Economic Review). The benchmark reference for DiD with a continuous treatment; focus on the selection bias of Theorem 3.2 and the TWFE decomposition of Theorem 3.4.
- de Chaisemartin and D'Haultfœuille (2026, Review of Economics and Statistics). The original paper on $DID_\ell$ and normalized effects; the setup in its §2 and §3, together with its Lemmas 1 to 3 (arXiv v14 numbering), is the source of Section 3.7 of this chapter.
- de Chaisemartin and D'Haultfœuille (2020, American Economic Review). The negative-weights theorem for TWFE and $DID_M$; the Web Appendix has versions for non-binary treatments.
- Baker, Callaway, Cunningham, Goodman-Bacon and Sant'Anna (2026, Journal of Economic Literature). The practitioner's guide already recommended in Chapter 19, worth reading alongside this chapter.

Further reading:

- de Chaisemartin, D'Haultfœuille, Pasquier, Sow and Vazquez-Bare (2026, working paper). The stayers design, AS and WAS, and their IV extension.
- de Chaisemartin, Ciccia, D'Haultfœuille and Knau (2026, working paper). The heterogeneous adoption design, the quasi-untreated test, and the linearity test for TWFE.
- de Chaisemartin, D'Haultfœuille and Vazquez-Bare (2024, AEA Papers and Proceedings) and Callaway, Goodman-Bacon and Sant'Anna (2024, AEA Papers and Proceedings). Two short papers: continuous treatments without stayers, and event studies with a continuous treatment.
- Athey and Imbens (2006, Econometrica). The original changes-in-changes paper and the foundation of Wald-CIC.
- Miyaji (2026, working paper) and Ye, Ertefaie, Flory, Hennessy and Small (2023, Biometrics). Two formulations of instrumented DiD; the former includes an extension to staggered designs.
- de Chaisemartin, D'Haultfœuille and Guyonvarch (2019, Stata Journal). Documentation of the fuzzydid command.
- de Chaisemartin and D'Haultfœuille (2023b, Journal of Econometrics). Mutual contamination when several treatment variables enter a TWFE regression together, and alternative estimators robust to heterogeneous effects.
- de Chaisemartin and D'Haultfœuille (2023a), the textbook draft. Until the final Princeton University Press version appears, read this SSRN draft.
- Wooldridge (2010) §21.4.3. The textbook treatment of LATE and the Wald estimator, to read against Section 3.8 of this chapter.
:::

::: {.apa-refs}
- Athey, S., & Imbens, G. W. (2006). Identification and inference in nonlinear difference-in-differences models. *Econometrica, 74*(2), 431-497. https://doi.org/10.1111/j.1468-0262.2006.00668.x
- Baker, A., Callaway, B., Cunningham, S., Goodman-Bacon, A., & Sant'Anna, P. H. C. (2026). Difference-in-differences designs: A practitioner's guide. *Journal of Economic Literature, 64*(2), 498-557. https://doi.org/10.1257/jel.20251650
- Callaway, B., Goodman-Bacon, A., & Sant'Anna, P. H. C. (2024). Event studies with a continuous treatment. *AEA Papers and Proceedings, 114*, 601-605. https://doi.org/10.1257/pandp.20241047
- Callaway, B., Goodman-Bacon, A., & Sant'Anna, P. H. C. (2025). *Difference-in-differences with a continuous treatment* (arXiv:2107.02637v8). arXiv. https://doi.org/10.48550/arXiv.2107.02637
- Calonico, S., Cattaneo, M. D., & Farrell, M. H. (2018). On the effect of bias estimation on coverage accuracy in nonparametric inference. *Journal of the American Statistical Association, 113*(522), 767-779. https://doi.org/10.1080/01621459.2017.1285776
- Card, D., & Krueger, A. B. (1994). Minimum wages and employment: A case study of the fast-food industry in New Jersey and Pennsylvania. *American Economic Review, 84*(4), 772-793.
- de Chaisemartin, C., Ciccia, D., D'Haultfœuille, X., & Knau, F. (2026). *Difference-in-differences estimators when no unit remains untreated* (arXiv:2405.04465v8). arXiv. https://doi.org/10.48550/arXiv.2405.04465
- de Chaisemartin, C., & D'Haultfœuille, X. (2018). Fuzzy differences-in-differences. *The Review of Economic Studies, 85*(2), 999-1028. https://doi.org/10.1093/restud/rdx049
- de Chaisemartin, C., & D'Haultfœuille, X. (2020). Two-way fixed effects estimators with heterogeneous treatment effects. *American Economic Review, 110*(9), 2964-2996. https://doi.org/10.1257/aer.20181169
- de Chaisemartin, C., & D'Haultfœuille, X. (2023a). *Causal inference with differences-in-differences: Credible answers to hard questions* [Unpublished manuscript]. SSRN. https://doi.org/10.2139/ssrn.4487202
- de Chaisemartin, C., & D'Haultfœuille, X. (2023b). Two-way fixed effects and differences-in-differences estimators with several treatments. *Journal of Econometrics, 236*(2), Article 105480. https://doi.org/10.1016/j.jeconom.2023.105480
- de Chaisemartin, C., & D'Haultfœuille, X. (2026). Difference-in-differences estimators of intertemporal treatment effects. *The Review of Economics and Statistics, 108*(4), 863-880. https://doi.org/10.1162/rest_a_01414
- de Chaisemartin, C., D'Haultfœuille, X., & Guyonvarch, Y. (2019). Fuzzy differences-in-differences with Stata. *The Stata Journal, 19*(2), 435-458. https://doi.org/10.1177/1536867X19854019
- de Chaisemartin, C., D'Haultfœuille, X., Pasquier, F., Sow, D., & Vazquez-Bare, G. (2026). *Difference-in-differences for continuous treatments and instruments with stayers* (arXiv:2201.06898v7). arXiv. https://doi.org/10.48550/arXiv.2201.06898
- de Chaisemartin, C., D'Haultfœuille, X., & Vazquez-Bare, G. (2024). Difference-in-difference estimators with continuous treatments and no stayers. *AEA Papers and Proceedings, 114*, 610-613. https://doi.org/10.1257/pandp.20241049
- Duflo, E. (2001). Schooling and labor market consequences of school construction in Indonesia: Evidence from an unusual policy experiment. *American Economic Review, 91*(4), 795-813. https://doi.org/10.1257/aer.91.4.795
- Favara, G., & Imbs, J. (2015). Credit supply and the price of housing. *American Economic Review, 105*(3), 958-992. https://doi.org/10.1257/aer.20121416
- Imbens, G. W., & Angrist, J. D. (1994). Identification and estimation of local average treatment effects. *Econometrica, 62*(2), 467-475. https://doi.org/10.2307/2951620
- Miyaji, S. (2026). *Instrumented difference-in-differences with heterogeneous treatment effects* (arXiv:2405.12083v6). arXiv. https://doi.org/10.48550/arXiv.2405.12083
- Roth, J. (2022). Pretest with caution: Event-study estimates after testing for parallel trends. *American Economic Review: Insights, 4*(3), 305-322. https://doi.org/10.1257/aeri.20210236
- Sant'Anna, P. H. C., & Zhao, J. (2020). Doubly robust difference-in-differences estimators. *Journal of Econometrics, 219*(1), 101-122. https://doi.org/10.1016/j.jeconom.2020.06.003
- Tunc, M. M., Cavusoglu, H., & Zheng, Z. (2025). Resale royalty in non-fungible token marketplaces: Blessing or burden for creators and platforms? *Information Systems Research, 36*(3), 1543-1564. https://doi.org/10.1287/isre.2023.0035
- Wooldridge, J. M. (2010). *Econometric analysis of cross section and panel data* (2nd ed.). MIT Press.
- Ye, T., Ertefaie, A., Flory, J., Hennessy, S., & Small, D. S. (2023). Instrumented difference-in-differences. *Biometrics, 79*(2), 569-581. https://doi.org/10.1111/biom.13783
:::
