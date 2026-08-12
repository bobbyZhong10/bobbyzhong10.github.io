---
title: "Potential Outcomes, Experiments & Selection-on-Observables"
subtitle: "From Randomized Experiments to Selection on Observables"
seriesline: "Foundations of Information Systems Economics · Chapter 9"
series: "Foundations of Information Systems Economics: Empirical Methods for Information Systems Economics"
sidebartitle: "Foundations of Information Systems Economics"
sidebarsub: "Theory · Methods · Applications"
chapterlabel: "Chapter 9 · Potential Outcomes, Experiments & Selection-on-Observables"
---

## Introduction

Northwind's client firms can adopt an AI customer-service assistant on their own. Adopters run 0.575 log points more efficiently than non-adopters, about 78 percent; a randomized pilot of the very same product estimates only 0.206, about 23 percent. Write the first number into a marketing deck and it is spectacular; write it into a causal study and the first thing a referee will ask is whether the firms that were already more advanced are simply the ones more willing to adopt.

A causal effect asks us to compare the same firm's outcome when it adopts and when it does not, yet the data let us observe only one of the two. The potential-outcomes framework writes down that missing counterfactual explicitly, and it forces the researcher to distinguish up front among the population average effect, the average effect on adopters, and the average effect on non-adopters. Randomization makes the control group a credible counterfactual through the assignment mechanism itself; in observational data, matching, the propensity score, inverse probability weighting, and doubly robust estimation all try to rebuild that comparability once we condition on covariates. They differ in their arithmetic, but they lean on the same two substantive assumptions: that selection is no longer driven by unobserved factors, and that every kind of firm has positive probability of appearing on both the treated and the control side.

This chapter puts the randomized pilot and the self-selected sample face to face inside the Northwind data. With the right adjustment we can pull 0.575 back to something near the true effect on adopters, 0.254; yet the moment unobserved selection enters the adoption decision, the same polished machinery reports 0.613. That failure matters more than the success: the propensity score can balance the variables you have already observed, but it will not go looking for the variables you left out. Later causal chapters swap in different data structures and different identifying assumptions, but none of them escapes the principle we establish here: an estimator can compute a counterfactual, but only the identifying assumption decides whether that counterfactual is credible.

## 1 A Suspiciously Large Number

::: {.case}
Northwind is a fictional B2B SaaS vendor that markets an AI customer-service assistant to a large base of client firms. We tell the case with simulated data, and the payoff of doing so is that the data-generating process is fully known: every estimator can be reconciled against the truth, a teaching condition that real data can never give us. The outcome is a client firm's monthly ticket-resolution efficiency, measured as log(tickets resolved per agent). Each firm carries a set of observable characteristics fixed before treatment: size `size`, a pre-existing digitization level `dig`, industry `ind` (four categories, A through D), and baseline efficiency `eff`.

The vendor does two things. The first is a randomized pilot (Study A): draw 2000 firms and flip a coin to decide who gets the assistant, which gives a clean benchmark for the effect. The second is a larger observational sample (Study B): 6000 firms decide for themselves whether to adopt, adoption is correlated with the characteristics above, and firms that expect a larger payoff are more willing to adopt.

Start with the most natural computation in Study B: take the mean efficiency of firms that adopted the assistant and subtract the mean efficiency of firms that did not. That difference is 0.575. A log-point gap of 0.575 means adopters resolve about 78 percent more tickets per agent, a number good enough to print in marketing copy. But the randomized pilot, Study A, gives a difference in means of only 0.206, about a 23 percent lift. Same product, same population of clients, one says 78 percent and the other says 23 percent, off by more than a factor of three.
:::

![Two answers for one product: the naive difference in Study B gives 0.575, while the randomized pilot in Study A gives only 0.206. The dashed line marks the true ATT of 0.254, which the naive difference overstates by a wide margin.](assets/fig/fig_09_raw_contrast.svg)

At least one of the two numbers is lying. The 0.206 from the randomized pilot is credible: the coin splits firms into two groups with no systematic difference before the assistant is switched on, so any efficiency gap afterward can only be attributed to the assistant itself. One caveat to note up front is that 0.206 measures the average effect of pushing the assistant onto every firm (the ATE), whereas what Study B really cares about is the effect on the firms that chose to adopt (the ATT), and the two need not be equal in the first place; Section 2 makes that distinction precise. But whichever one you aim at, neither reaches anywhere near 0.575. The truly suspicious number is the 0.575 from Study B. It is inflated because, in the observational sample, the firms that adopt the assistant and the firms that do not were never the same kind of firm to begin with. The firms willing to be early adopters of an AI tool tend to be larger, better digitized, and simply better run, and firms like that would post higher ticket efficiency even without switching the assistant on. That pre-existing gap is folded into the 0.575 and counted, in its entirety, as the assistant's doing.

What this chapter sets out to do is take that 0.575 apart, see how much of it is a real effect and how much is an artifact of selection; spell out the conditions under which we can strip the artifact away and keep only the real effect; and give a full set of estimators that put it back. By Section 6 we will see that, under the right assumptions, covariate adjustment on Study B pulls the estimate from 0.575 back to somewhere between 0.22 and 0.25, hovering around the true ATT of 0.254. But we will also see that this repair rests on one very tightly stretched premise, and once that premise snaps, the same method returns a 0.613 that looks careful but is just as badly wrong as before.

## 2 The Economic Model and the Estimand

Before doing any estimation, we have to be clear about two things: which causal quantity we actually want, and what language we use to write it down precisely. The latter is the potential-outcomes notation this section will set up, and the former is the estimand (the target quantity) that rides on top of that notation.

### 2.1 Potential Outcomes and the Fundamental Problem

Consider a binary treatment $D_i \in \{0, 1\}$; in Northwind, $D_i = 1$ means firm $i$ adopted the AI assistant. Each unit carries two **potential outcomes**: $Y_i(1)$ is the efficiency it would realize if it adopted, and $Y_i(0)$ is the efficiency it would realize if it did not. Both quantities exist conceptually at the same time, whether or not the firm actually adopts. The individual treatment effect is their difference:

$$\tau_i \equiv Y_i(1) - Y_i(0).$$

The crucial point is that $\tau_i$ varies across units. Some firms take off once the assistant is switched on; others already run a smooth process and gain little on the margin; a few are even thrown off rhythm by the new tool. Treatment effects are heterogeneous, and ideally we want to describe the whole distribution $f(\tau_i)$, not just a single average.

The **fundamental problem of causal inference** is this: for any single unit, we can observe only one of $Y_i(1)$ and $Y_i(0)$; the other is the parallel-world path we never travel, called the **counterfactual**. For a firm that adopted the assistant, we see its post-adoption efficiency $Y_i(1)$ but not the $Y_i(0)$ it would have had without adopting; for a firm that did not adopt, it is exactly the reverse. A table makes this most vivid:

| Firm type | $Y_i(1)$ | $Y_i(0)$ | $D_i$ | observed $Y_i$ |
|---|---|---|---|---|
| Adopter | observed | ? | 1 | $Y_i(1)$ |
| Non-adopter | ? | observed | 0 | $Y_i(0)$ |

Every row is missing a cell, and the missing cell is exactly the counterfactual. The individual effect $\tau_i = Y_i(1) - Y_i(0)$ needs both cells in the same row, and we only ever hold one. Causal inference was never a missing-data problem that a simple imputation could patch; it is a problem in which every row is guaranteed to be missing half.

Observed outcomes and potential outcomes are tied together by the **switching equation**:

$$Y_i = D_i\, Y_i(1) + (1 - D_i)\, Y_i(0).$$

The treatment switch $D_i$ decides which potential outcome gets flipped onto the table as the observed value. The equation looks bland, but it is the hinge of everything: it is the only bridge between the data (the left side) and the causal object we actually care about (the right side).

### 2.2 Three Estimands and How They Differ

Since the individual effect $\tau_i$ is unavailable one by one, empirical work settles for identifying some moment of $f(\tau_i)$, that is, some kind of average effect. Three are common. The **average treatment effect** (ATE) averages over all units:

$$\text{ATE} = \mathbb{E}[\tau_i] = \mathbb{E}[Y_i(1) - Y_i(0)].$$

The **average treatment effect on the treated** (ATT) averages only over the units actually treated:

$$\text{ATT} = \mathbb{E}[\tau_i \mid D_i = 1].$$

The **average treatment effect on the untreated** (ATU) covers the other half: $\text{ATU} = \mathbb{E}[\tau_i \mid D_i = 0]$. Writing the adoption rate as $\pi \equiv \Pr(D_i = 1)$, the three are strung together by a mixing identity:

$$\text{ATE} = \pi \cdot \text{ATT} + (1 - \pi) \cdot \text{ATU}.$$

These three quantities answer different policy questions, and mixing them up leads to error. The ATT asks, "for the firms that already adopted the assistant, how much lift did it deliver?", and it evaluates an accomplished fact. The ATE asks, "if we made a randomly drawn firm adopt the assistant, how much lift would it deliver on average?", and it speaks to the hypothetical policy of rolling treatment out to everyone. The ATU asks, "for the firms that have not adopted so far, how much lift would they get if we persuaded them to?", and it is the most fitting target for the question of whether to promote adoption. In Northwind these three numbers are not equal: the true ATT in Study B is 0.254 while the true ATE is 0.200, so the ATT sits clearly above the ATE. The reason is spelled out in the next section, but in short, the firms that adopt are precisely the ones that gain the most from the assistant, so the average payoff among the treated is naturally higher than the population average. This gap is not a flaw in the estimation but information carried by the setting itself: a cancer drug is trialed specifically on cancer patients, and we should never expect its effect to extrapolate to a healthy population.

::: {.intuition}
Think of choosing the estimand as deciding "for whom you are writing the prescription." The ATT looks back at how the patients who already took the drug recovered; the ATE imagines the average effect of dosing the whole population; the ATU weighs whether to prescribe to those who have not taken it yet. When treatment effects are heterogeneous, the answers to the three questions can differ widely. Decide which one you want first, then pick the estimator, and never reverse the order; many empirical disputes look like "the method was wrong" on the surface, when underneath two people simply wanted different estimands.
:::

### 2.3 The Economic Structure of the Adoption Decision

In Northwind's observational sample, firms are not randomly assigned to treatment; they run the numbers and decide for themselves whether to adopt. A firm switches the assistant on roughly when it expects the net payoff of doing so to be positive. Written in potential outcomes, the adoption decision is approximately

$$D_i = \mathbf{1}\{\, Y_i(1) - Y_i(0) > c_i \,\},$$

where $c_i$ is the cost of, or threshold for, adoption. This is Roy-type self-selection: whether a unit enters treatment depends on its own payoff $\tau_i = Y_i(1) - Y_i(0)$. This structure has an immediate consequence, namely that adopters are a group of units with systematically higher payoffs, so that $\text{ATT} = \mathbb{E}[\tau_i \mid \tau_i > c_i]$ generally exceeds $\text{ATE} = \mathbb{E}[\tau_i]$ strictly. It is at the same time the economic engine of selection bias, and the decomposition in Section 3 will translate this intuition into precise algebra.

Finally, one premise runs through the whole chapter: **SUTVA** (the stable unit treatment value assumption). It has two layers. The first is no interference: each unit's outcome depends only on whether that unit itself is treated, not on whether others are treated. The second is that treatment comes in a single, well-defined version. Applied to Northwind, the first layer requires that whether one firm switches the assistant on does not, through any channel, change another firm's ticket efficiency. The concrete ways it could fail are specific: if Northwind's compute is shared, a wave of simultaneous adoptions would slow each firm's responses; if the assistant's effectiveness comes from knowledge accumulated across firms, early adopters would spill over and help late adopters. All the analysis in this chapter is built on SUTVA holding, and once it breaks, even the starting point that "there are only two potential outcomes" has to be rewritten.

To summarize: the potential-outcomes framework writes a causal effect as the difference between two parallel-world outcomes, and the fundamental problem is that each unit reveals only one; the switching equation connects observed outcomes to potential outcomes; the estimand comes in three forms, ATE, ATT, and ATU, which are unequal under heterogeneous effects, so we must decide which one we want first; and Northwind's adoption decision is Roy-type self-selection, which foreshadows both that the ATT exceeds the ATE and where selection bias comes from.

## 3 Identification: From Randomization to Selection-on-Observables

Identification logically precedes estimation. Identification asks whether, under a set of assumptions about the counterfactual world, the causal quantity we want can be written as a function of the observable data distribution; the answer depends only on the design and the institutional background, not on which estimator or how large a sample we use afterward. This is the most detailed analytical section of the chapter, and it works down layer by layer: first we drive the skeleton of identification all the way through in a randomized experiment, then we relax to observational data and see what extra assumptions it costs to re-identify the causal effect. The through-line of the whole section is one decomposition, and it will run like a red thread through everything that follows.

### 3.1 What the Naive Comparison Is Actually Estimating

Return to the 0.575 from Section 1. It is the observed mean of adopters minus the observed mean of non-adopters, $\mathbb{E}[Y_i \mid D_i = 1] - \mathbb{E}[Y_i \mid D_i = 0]$. This quantity looks like an effect, but it is not. To see what it is really estimating, all we need is the small trick of adding and subtracting the same term.

Warm up with a purely illustrative numerical example. Suppose that if adopters had not adopted, their average efficiency would have been 3.2; after adopting it is 3.5; and non-adopters average 3.0. Then the naive comparison gives $3.5 - 3.0 = 0.5$. But of that 0.5, only the adopters' own before-and-after difference, $3.5 - 3.2 = 0.3$, is genuinely due to treatment; the remaining $3.2 - 3.0 = 0.2$ is the amount by which adopters would have exceeded non-adopters even without adopting, purely because the two groups start from different places. The naive comparison counts that 0.2 as effect too.

Write the intuition as general algebra. Start from the observed difference, use the switching equation to turn observed means into potential-outcome means, then add and subtract the treated group's untreated mean $\mathbb{E}[Y_i(0) \mid D_i = 1]$:

$$
\underbrace{\mathbb{E}[Y_i(1) \mid D_i = 1] - \mathbb{E}[Y_i(0) \mid D_i = 0]}_{\text{observed difference}}
= \underbrace{\mathbb{E}[Y_i(1) \mid D_i = 1] - \mathbb{E}[Y_i(0) \mid D_i = 1]}_{\text{ATT}}
+ \underbrace{\mathbb{E}[Y_i(0) \mid D_i = 1] - \mathbb{E}[Y_i(0) \mid D_i = 0]}_{\text{selection bias}}.
$$

This is the skeleton decomposition of the whole chapter: the naive comparison equals the ATT plus a **selection bias** term. The definition of selection bias is worth studying again and again: it is the difference in means between adopters and non-adopters in the same untreated state $Y(0)$, that is, "how much the two groups of firms would have differed had neither adopted the assistant." In the numerical example above it is exactly $3.2 - 3.0 = 0.2$. Put into Northwind's real numbers, the naive comparison 0.575 minus the true ATT 0.254 leaves about 0.32, and that entire residual is selection bias. In that suspiciously large number from Section 1, more than half (0.322 out of 0.575) is held up by selection bias.

Selection bias has a counterintuitive but extremely important property: it can be nonzero even when treatment effects are perfectly homogeneous. Suppose every firm's effect is exactly the same constant $\tau_i = \tau$, and substitute $Y_i(0) = \alpha + \varepsilon_i$; then selection bias is $\mathbb{E}[\varepsilon_i \mid D_i = 1] - \mathbb{E}[\varepsilon_i \mid D_i = 0]$. As long as the baseline disturbance $\varepsilon_i$ is correlated with the adoption decision $D_i$, it is nonzero. This is precisely an **endogeneity** problem: adopters' baselines differ from non-adopters' to begin with. Private-school children score better not entirely because of the school, since the families that choose private schooling can often also afford tutors and have more education themselves, so the baseline is higher to begin with. That slice of baseline difference is something no naive "compare adopters with non-adopters" procedure can strip away.

### 3.2 The Second Bias from Heterogeneous Effects

The decomposition above split the naive comparison into ATT plus selection bias. If what we truly want is the ATT, the problem is clean: we only need to kill the selection bias. But if what we want is the ATE, there is one more account to settle, because the ATT itself need not equal the ATE. Using the mixing identity, we can swap the ATT for the ATE plus a difference term, $\text{ATT} = \text{ATE} + (1 - \pi)(\text{ATT} - \text{ATU})$, and substituting back gives the full decomposition aimed at the ATE:

$$
\mathbb{E}[Y_i \mid D_i = 1] - \mathbb{E}[Y_i \mid D_i = 0]
= \text{ATE}
+ \underbrace{(1 - \pi)\big(\text{ATT} - \text{ATU}\big)}_{\text{heterogeneity bias}}
+ \underbrace{\big(\mathbb{E}[Y_i(0) \mid D_i = 1] - \mathbb{E}[Y_i(0) \mid D_i = 0]\big)}_{\text{selection bias}}.
$$

Now there are two biases. Selection bias is as before, coming from the between-group difference in baseline levels. The newly appearing **heterogeneity bias** comes from the difference in average payoffs between the treated and untreated groups: $\text{ATT} - \text{ATU}$ measures exactly "how much more adopters gain from the assistant than non-adopters would." When effects are homogeneous, or when the two groups happen to have equal average payoffs, this term vanishes and the ATT equals the ATE; but once there is Roy-type payoff self-selection, so that adopters are precisely the highest-payoff group and $\text{ATT} > \text{ATU}$, this term is strictly positive. This is the mechanism behind the ATT exceeding the ATE noted in Section 2. In Northwind, the true ATT of 0.254 exceeds the true ATE of 0.200, and the excess is exactly this payoff self-selection. The Roy model in Section 3.8 will tell it as a complete economic story.

### 3.3 How Randomization Erases Selection Bias

Now introduce the chapter's first identification route: randomization. If $D_i$ is randomly assigned, independent of the potential outcomes,

$$D_i \perp \big(Y_i(0), Y_i(1)\big),$$

then whether a firm adopts has nothing to do with any of its potential outcomes. This one independence condition cuts straight through the decomposition in Section 3.1. Because $Y_i(0)$ is independent of $D_i$, the untreated means of the adopter and non-adopter groups must be equal, $\mathbb{E}[Y_i(0) \mid D_i = 1] = \mathbb{E}[Y_i(0) \mid D_i = 0] = \mathbb{E}[Y_i(0)]$, and selection bias drops to zero on the spot. What is more, both potential outcomes are independent of $D_i$, which means the effect estimated on the treated group can be generalized to everyone:

$$\mathbb{E}[Y_i(1) \mid D_i = 1] - \mathbb{E}[Y_i(0) \mid D_i = 1] = \mathbb{E}[Y_i(1)] - \mathbb{E}[Y_i(0)] = \text{ATE},$$

so heterogeneity bias vanishes too, and the ATT merges with the ATE. The power of randomization can be summarized in two sentences: first, randomization eliminates selection bias; second, randomization makes the ATT equal the ATE. This is the entire reason the 0.206 from the randomized pilot, Study A, is credible. It needs no model of firm characteristics whatsoever, since the act of flipping a coin by itself guarantees the two groups are comparable.

::: {.intuition}
The magic of randomization is that it takes selection out of the firms' hands and gives it to a coin. If firms can self-select, adopters come pre-loaded with characteristics correlated with the outcome, and the comparison is crooked from the start; a coin looks at none of the firms' characteristics, so the two groups it produces are, in expectation, two copies of the same firms, and any difference between the copies afterward can only come from treatment. This also explains why randomized experiments are held up as the gold standard: not because the sample is larger or the measurement more precise, but because at the design level they force the hardest bias, selection bias, to equal zero identically. Everything an observational study does is, in essence, an attempt to reconstruct by hand, without a coin, the property that randomization hands us for free.
:::

Clean as it is, randomization also has limits it cannot reach. Real experiments meet imperfect execution: attrition, subjects manipulating their group assignment, insufficient compliance. There is also a subtler class of problem, that treatment itself may directly change the disturbance term; for instance, merely "knowing you have been given a new feature" may make a firm operate more diligently, and this placebo-style behavioral response mixes into the effect, which is exactly what double-blind designs guard against. On top of this, the experimental sample need not represent the population, and external validity is always an open question. All of this reminds us that randomization solves the single problem of selection bias, not every problem.

### 3.4 Relaxing to Observational Data: The Conditional Independence Assumption

In reality a great many studies cannot obtain randomization. A platform does not flip a coin to decide which firms get the assistant, and firms do not submit to a lottery. Here the only hope is this: although adoption is not unconditionally random, once we control for enough observable characteristics, adoption is approximately random. This is the chapter's second identification route, selection-on-observables. Its core is the following assumption.

::: {.assumption}
**CIA (conditional independence assumption; also called unconfoundedness or ignorability):** conditional on observable covariates $X_i$, the potential outcomes are independent of treatment:

$$\{Y_i(1), Y_i(0)\} \perp D_i \mid X_i.$$
:::

The CIA says that once $X_i$ is held fixed, adoption is "as good as random." It does not require adopters and non-adopters to be comparable overall; it only requires that within every small cell where $X$ takes the same value, the potential-outcome distributions of the adopter and non-adopter groups are identical. Put differently, the only systematic difference the two groups are allowed is a different covariate composition $f(X_i \mid D_i)$: adopters may be larger on average and better digitized, but as long as within the cell of "same size, same digitization, same industry, same baseline efficiency" who adopts and who does not is unrelated to their potential efficiencies, the CIA holds. It is equivalent to saying that every factor that jointly influences adoption and outcome has already been packed into $X_i$.

To judge whether the CIA is credible, you must interrogate, one by one, whether any confounder has been omitted. In Northwind, if firms adopt mainly on the basis of size, digitization, industry, and baseline efficiency, and all four are in $X$, then the CIA has a chance of holding. But if there remains even one factor that drives both adoption and efficiency yet goes unobserved, say "the management team's operational drive," the CIA breaks: a driven management team is both more willing to try adopting the assistant and better at running the ticket process itself, and with that drive lying outside $X$, the adopter group's $Y(0)$ still systematically exceeds the non-adopter group's.

::: {.warning}
The CIA is untestable. It constrains the counterfactual arm we can never observe: how adopters would have fared had they not adopted. That information simply is not in the data, and no statistical test can touch it. This point cannot be stressed too much: if your identification rests on the propensity score, you owe a full paragraph arguing why the CIA is credible in your setting, spelling out which confounders you observed and why nothing important is omitted, rather than tossing out a balance table and calling it done. Balance tests, placebo tests, and pre-treatment outcome tests are only indirect circumstantial evidence; they check the observable shadow, not the object that casts it.
:::

### 3.5 Overlap: The Other Half of the Assumption

The CIA guarantees "if the two groups are comparable within the same $X$ cell," but we still need to guarantee "that both groups actually have members in each cell," or there is nothing to compare. That is the second assumption.

::: {.assumption}
**Overlap (also called common support or positivity):** for every possible covariate value $X_i$, both treated and untreated occur with positive probability:

$$0 < e(X_i) < 1, \qquad e(X) \equiv \Pr(D = 1 \mid X).$$
:::

Here $e(X)$ is the propensity score, which takes the stage formally in the next section. The meaning of overlap is plain: we can only learn something in the regions of $X$ that contain both adopters and non-adopters. If some class of firms (say the especially large and highly digitized ones) all adopted the assistant, so that $e(X) = 1$, then that cell has no non-adopting control at all, their counterfactual $Y(0)$ has nothing to draw on, and we are left to extrapolate from the model, and extrapolation is a guess with no data behind it. Conversely, if some class of firms never adopts, so that $e(X) = 0$, their $Y(1)$ is equally out of reach.

The CIA and overlap together are called **strong ignorability** (Rosenbaum and Rubin's term). The testability of these two assumptions is starkly different, and this is a practical distinction. The CIA is untestable, as the previous section stressed repeatedly. Overlap, by contrast, can be checked: plot the estimated propensity score as two distributions by adoption status and see whether their supports overlap and whether mass piles up near 0 or 1. Section 6 will draw this figure for Northwind, where the propensity estimates fall in $[0.018, 0.977]$ with no pile-up at either end, a sample with benign overlap. Section 7 then discusses what to do when overlap is inadequate, and what trouble "overcorrecting trimming" itself brings.

### 3.6 Identification Under the CIA

With these two assumptions in hand, we can write the causal quantity as a function of the observable data. Take the ATT, and write $F^1(x)$ for the covariate distribution of the adopter group. The ATT can be written as the difference of two conditional means over the adopter group, integrated against $F^1$:

$$\text{ATT} = \mathbb{E}_{F^1(x)}\big[\mathbb{E}(Y(1) \mid D = 1, X)\big] - \mathbb{E}_{F^1(x)}\big[\mathbb{E}(Y(0) \mid D = 1, X)\big].$$

The first term is the adopter group's own observed mean, directly visible. The second term is the counterfactual mean of adopters had they not adopted, and it is invisible. The CIA is exactly what does the work here: within each $X$ cell, the adopter group's $Y(0)$ distribution equals the non-adopter group's $Y(0)$ distribution, so that invisible counterfactual mean can be replaced by the non-adopter group's observed mean at the same $X$:

$$\mathbb{E}_{F^1(x)}\big[\mathbb{E}(Y(0) \mid D = 1, X)\big] = \mathbb{E}_{F^1(x)}\big[\mathbb{E}(Y(0) \mid D = 0, X)\big].$$

This substitution is the entire essence of selection-on-observables methods. It says that to build the adopters' counterfactual, go into the non-adopter group and find firms with the same $X$, use their outcomes as the estimate of what adopters would have had "if they had not adopted," then take a weighted average according to the adopter group's covariate distribution. Matching finds the most similar neighbor one at a time, the propensity score first squeezes $X$ into one dimension before matching, and IPW reweights the non-adopters to the adopter group's $X$ distribution: different roads to the same place, all of them different implementations of this one substitution. Overlap is necessary here: if some $X$ cell has no non-adopters, the conditional mean on the right simply does not exist.

::: {.warning}
The CIA holding on $X$ does not imply it still holds on $(X, Z)$ after adding one more variable $Z$, and vice versa. Worse still, one can construct situations in which within every $(X, Z)$ cell treatment is correlated with the outcome, yet once we average over $Z$ the correlation disappears from the aggregate. Piling covariates blindly into the matching variables can wreck an otherwise clean identification. This is the entry point to the bad-controls problem that Section 7 develops: control variables must be measured pre-treatment and must never be downstream variables that respond to treatment, or else controlling for one is like blocking a channel through which the effect travels, and can even introduce collider bias.
:::

### 3.7 The Propensity Score Theorem

Matching cell by cell on $X$ has a practical obstacle. Once $X$ is multidimensional, the number of cells swells exponentially with dimension, each cell holds fewer and fewer observations, and finding a control with exactly the same $X$ becomes all but impossible: this is the curse of dimensionality. Rosenbaum and Rubin's (1983) theorem offers a way out: we need not match on all of $X$; matching on a single scalar is enough.

First a definition. The **propensity score** is the conditional probability of being treated given the covariates, $e(X_i) \equiv \Pr(D_i = 1 \mid X_i)$. A function $b(X)$ satisfying $D \perp X \mid b(X)$ is called a **balancing score**: conditioning on it balances the distribution of the covariates across the two groups. The propensity score is the coarsest balancing score, and any balancing score is finer than it.

::: {.theorem}
**(Rosenbaum-Rubin 1983)** If treatment is strongly ignorable given $X$ (the CIA holds), then it is strongly ignorable given any balancing score, and in particular given the propensity score:

$$\{Y(1), Y(0)\} \perp D \mid X \;\Longrightarrow\; \{Y(1), Y(0)\} \perp D \mid e(X),$$

provided $0 < e(X) < 1$.
:::

The proof takes only two steps. First, show that conditional on $e(X)$ the treatment probability is $e(X)$ itself:

$$\Pr\big(D = 1 \mid Y(1), Y(0), e(X)\big) = \mathbb{E}\big[\Pr(D = 1 \mid Y(1), Y(0), X) \mid e(X)\big] = \mathbb{E}\big[\Pr(D = 1 \mid X) \mid e(X)\big] = e(X),$$

where the middle step uses the CIA (given $X$, $D$ is unrelated to the potential outcomes, so the potential outcomes in the conditioning set can be dropped) and the last step uses the definition of the propensity score. Since this conditional probability does not depend on the potential outcomes, $D$ is independent of the potential outcomes given $e(X)$. The value of the theorem is that it compresses matching from $K$ dimensions to one: as long as we compare units with equal $e(X)$, it is equivalent to having balanced the covariates on all of $X$. Identifying the counterfactual mean for the ATT can then be written as matching on the one-dimensional $e(X)$:

$$\mathbb{E}_{F^1(x)}\big[\mathbb{E}(Y(0) \mid D = 1, X)\big] = \mathbb{E}_{F^1(x)}\big[\mathbb{E}(Y(0) \mid D = 0, e(X))\big].$$

The price is that $e(X)$ itself must be estimated. The theorem assumes $e(X)$ is known; in practice it is unknown and is usually fit with a logit or probit. This optimism, "assume we have a can opener," deserves an honest label: the credibility of propensity-score methods rests both on the untestable CIA and on the propensity model being correctly specified. The doubly robust estimator of Section 4 exists precisely to ease the second fragility.

### 3.8 The Roy Selection Model: Why the ATT Is Not the ATE

The economic root of the heterogeneity bias in Section 3.2 is Roy-type payoff self-selection. Here we write it as a complete model. Let each potential outcome carry an idiosyncratic shock, $Y_i(1) = \mu_1 + U_{1i}$ and $Y_i(0) = \mu_0 + U_{0i}$, and let a unit choose treatment only when the net payoff exceeds a cost $c$:

$$D_i = \mathbf{1}\{\, Y_i(1) - Y_i(0) > c \,\}.$$

This is **selection on gains**: whether a unit enters treatment depends on its own payoff $\tau_i$. The adopter group is therefore the right tail of the payoff distribution cut out, $\text{ATT} = \mathbb{E}[\tau_i \mid \tau_i > c]$, strictly above the population average $\text{ATE} = \mathbb{E}[\tau_i]$, and so the heterogeneity bias term $(1 - \pi)(\text{ATT} - \text{ATU})$ from Section 3.2 is positive. Who goes to college are the people whose expected return is higher than most; who undergoes high-risk surgery are the ones sickest and most likely to benefit from it. The firms that adopt an AI assistant are, likewise, precisely those that predict they will gain the most from it.

There is a fine but important dividing line here. Selection on gains does not necessarily destroy the CIA. If the payoff $\tau_i$ can be predicted from observable $X_i$, then although selection falls on the gains, in essence it still falls on $X$; conditioning on $X$ leaves the CIA intact, and matching-type methods can still recover the ATT. What is truly fatal is selection falling on unobservable gains or unobservable levels: only when adoption also depends on some factor that both influences $Y(0)$ and lies outside $X$ does the CIA genuinely break. Northwind's simulation distinguishes these two cases precisely, and Sections 6 and 7 will show that when selection is on observable gains, adjustment can rescue the ATT, but once selection attaches to an unobservable factor, the same adjustment leaves a residual bias that cannot be wiped away. This crack is exactly the exit toward instrumental variables in Chapter 12.

::: {.intuition}
The ATT not equaling the ATE is not a bug that needs fixing; it is the setting speaking. "The higher a firm's payoff, the more it will adopt" already means that "the average payoff among adopters" exceeds "the average payoff of a random firm." If someone reports an ATT and treats it as an ATE when discussing rollout policy, it is like using top students' tutoring gains to predict the average lift from tutoring the whole class, and is bound to overstate. Heterogeneity is not a nuisance; it is valuable information about "who should most be treated," provided you do not confuse two different averages.
:::

### 3.9 The CIA and Parallel Trends Are Two Different Routes

Selection-on-observables is not the only route to rebuilding comparability. There is another common line of thinking, worth planting a marker for here and taking up in detail later in this series: rather than assuming the two groups' untreated levels are comparable after controlling for covariates, assume the two groups' untreated outcomes trend in the same way, which is parallel trends. The two borrow the counterfactual in different ways. The CIA works in the level dimension: it assumes that after controlling for observable characteristics, the untreated outcome levels $Y(0)$ of adopters and non-adopters are comparable, that is, it forbids any omitted, treatment-correlated difference in levels. Parallel trends works in the trend dimension: it allows the two groups to have arbitrary, even strongly treatment-correlated, fixed differences in the level of $Y(0)$, and requires only that the two groups' $Y(0)$ trend in the same way.

So the two are mirror images. Matching imposes almost no restriction on functional form, at the cost of requiring selection on observables, and any unobservable selection kills it; the trend-controlling route makes a strong additive-separability assumption, buying tolerance of arbitrary unobservable heterogeneity in levels, at the cost of a knife-edge requirement on the trend. Which is more credible depends on the setting: when you worry that adopters were a different kind of firm all along, that the difference shows up in levels, and you have panel data to line up before and after, controlling trends often fits better; when the between-group difference looks more like it is driven by observable characteristics and you have only a cross section, the matching line of thinking is more natural. The full mechanics of this route, and how its conditional version carries the matching-style conditional-independence idea over to the growth of the outcome, are left for the Difference-in-Differences discussion later in this series.

The key points of this section, summarized: the naive comparison equals the ATT plus selection bias, and for the ATE you must further subtract a heterogeneity bias; randomization erases both biases at once and is the gold standard; in observational data we rebuild comparability with the CIA plus overlap, where the CIA is untestable while overlap can be checked; the propensity score drops matching to one dimension; Roy-type payoff self-selection explains why the ATT exceeds the ATE and draws the boundary between where the CIA holds and where it fails; and the CIA and parallel trends are two different routes, controlling levels versus controlling trends.

## 4 Estimation

This section advances in the order of "where the previous method falls short, and so gives rise to the next." Every estimator serves the same substitution from Section 3: fill in the adopters' counterfactual with comparable control units.

### 4.1 Difference in Means and Randomization Inference

In a randomized experiment estimation is simplest, and the difference in means is the answer:

$$\hat\tau = \bar Y_{D=1} - \bar Y_{D=0}.$$

Section 3.3 already proved it consistently estimates the ATE (which also equals the ATT) under randomization. Its standard error can use Neyman's conservative formula $\sqrt{\operatorname{var}(Y_1)/n_1 + \operatorname{var}(Y_0)/n_0}$, paired with large-sample normal critical values for inference. This inference relies on asymptotic normality and need not be reliable when the sample is small.

An alternative that leans on no distributional assumption is **randomization inference** (also called permutation inference). It is built on Fisher's sharp null: $H_0: Y_i(1) = Y_i(0)$ for every unit, that is, treatment has no effect on anyone. This null is special: once it is true, both potential outcomes of every unit are known (both equal the observed value), so the distribution of any test statistic across all possible treatment assignments can be computed exactly. The procedure is to repeatedly permute the treatment vector $\mathbf{D}$ in the ways the design allows, recompute the statistic each time to obtain its exact distribution under the null, and read off how extreme the observed statistic is as an exact p-value:

$$p = \frac{1}{|\mathcal{D}|} \sum_{\mathbf{d} \in \mathcal{D}} \mathbf{1}\big\{ |T(\mathbf{d}, \mathbf{Y})| \ge |T(\mathbf{d}^{obs}, \mathbf{Y})| \big\}.$$

It is finite-sample exact, makes minimal assumptions, and does not lean on asymptotic normality. Section 6 will run 5000 permutations on Study A and obtain a p-value below 0.001.

Designing an experiment also requires answering "how large a sample do we need." Take the difference of two proportions: the Wald statistic is asymptotically standard normal, and given two-sided significance $\alpha$ and power $1 - \beta$, the sample size required per arm is

$$n = \frac{(z_{1-\alpha/2} + z_{1-\beta})^2 \big(p_1(1 - p_1) + p_0(1 - p_0)\big)}{(p_1 - p_0)^2}.$$

This formula holds a sobering number. If the baseline conversion rate is $p_0 = 2\%$, detecting a 10 percent relative lift (an absolute difference of only 0.2 percentage points) requires about 80,681 units per arm; if the relative lift is 50 percent, only about 3,823 per arm are needed. Power is extremely sensitive to the absolute difference, and at low baseline rates even a large relative lift can require an enormous sample because the absolute difference is so small. This is a bucket of cold water on "randomized experiments solve everything": before asking whether you can run an experiment, work out whether the sample size is feasible.

### 4.2 Regression Adjustment and Its Pitfalls

In observational data the first idea is usually regression. Run $Y = \alpha + \beta D + \gamma' X + v$ and use the coefficient $\hat\beta$ as the effect. Under conditional mean independence $\mathbb{E}(v \mid D, X) = \mathbb{E}(v \mid X)$, this $\hat\beta$ consistently estimates the effect, and OLS can be seen as a special case of matching (approximating the conditional mean with a linear model).

The pitfall is heterogeneous effects. The additive specification above implicitly makes the treatment effect a constant. When the effect really varies with $X$, the OLS $\hat\beta$ is no longer a clean ATE or ATT but a weighted average whose weights are the variances of the treatment probabilities, weights that pile onto units with propensity near 1/2, with no policy interpretation. In Northwind, the $D$ coefficient of this additive OLS is 0.184, which actually lands near the ATE or even lower; it is not the ATT.

To truly aim at the ATT, you have to let the effect vary with covariates, that is, run a regression fully interacting $D$ with $X$, then center all covariates (including the industry dummies) at their adopter-group means, so that the $D$ coefficient reads as the average effect on adopters. The "all" here brooks no discount: if you leave even one covariate uncentered (say leaving the industry dummies untouched), the $D$ coefficient reads only a conditional effect in the reference industry, at the means of the remaining covariates, not the ATT. In Northwind this fully centered interacted specification gives 0.239, higher than the additive OLS's 0.184, correctly landing back near the ATT of 0.254. The additive coefficient 0.184 and the interacted ATT 0.239 are two different quantities, and conflating them is a common mistake.

### 4.3 Nearest-Neighbor Matching and Its Variance Trap

Regression fills the counterfactual by model extrapolation; matching takes a different tack and directly finds, for each unit, the most similar neighbor from the opposite group. For each unit $i$, take its $M$ nearest neighbors $\mathcal{J}_M(i)$ in the other group and use the average of the neighbors' outcomes to fill in the missing potential outcome: an adopter's $\widehat{Y}_i(0)$ is the mean of its $M$ non-adopting neighbors, a non-adopter's $\widehat{Y}_i(1)$ is the mean of its $M$ adopting neighbors, and the other arm uses the observed value itself. With both arms in hand, one can compute various average effects. The "nearest" metric commonly uses the Mahalanobis distance $d(\mathbf{x}, \mathbf{y}) = (\mathbf{x} - \mathbf{y}) S^{-1} (\mathbf{x} - \mathbf{y})'$, which standardizes each dimension by the covariance and is a generalization of Euclidean distance.

Matching has two traps to call out. First, it suffers the curse of dimensionality just the same: as the dimension grows, the "nearest" neighbor is in fact far away, and matching degenerates into extrapolation, which is why one either uses the propensity score to drop the dimension to one, or applies a regression-style bias correction to the matches (correcting for the covariate difference between neighbor and unit using within-group regression fitted values). Second, and this is the most classic trap, matching standard errors cannot be computed as if they were ordinary regression standard errors. Under matching with replacement, the same control unit may be reused as the neighbor of many adopters, and this reuse makes the matching estimator a non-smooth weighted sum of the observed values.

::: {.warning}
Abadie and Imbens (2006) derived the large-sample distribution of nearest-neighbor matching with replacement and a consistent analytic variance estimator, which requires an additional within-group nearest-neighbor variance estimate; the ordinary regression standard-error formula is wrong here. Worse, Abadie and Imbens (2008) proved that the bootstrap is invalid for nearest-neighbor matching with a fixed number of matches: the matching functional is not smooth, the number of times a given unit is used as a neighbor does not converge to a smooth limit, and the bootstrap standard error is therefore inconsistent. The correct practice is to use the Abadie-Imbens analytic variance (directly available in `Matching` and `MatchIt`) and not to wrap a matching estimator in `boot`. Section 6 will show that in Northwind the Abadie-Imbens SE is 0.021 while running OLS directly on the matched sample gives an SE of 0.025; this time the two happen to be numerically close, but whether they are close is luck, and what matters is which one stands up methodologically, not how much they differ.
:::

### 4.4 Propensity Score Matching and IPW

With the propensity score, there is no need to match on multidimensional $X$; matching directly on the one-dimensional $e(X)$ suffices, which is propensity-score matching. It can be done as kernel matching, using non-adopters with nearby propensity to form a weighted-average counterfactual for each adopter:

$$\widehat{\text{ATT}} = \frac{1}{N_1} \sum_{i: D_i = 1} \left[ Y_i - \frac{\sum_{j: D_j = 0} Y_j\, K\big(e(X_i) - e(X_j)\big)}{\sum_{s: D_s = 0} K\big(e(X_i) - e(X_s)\big)} \right].$$

The bandwidth $h$ controls the bias-variance tradeoff, and common support is still the premise.

Another route that uses the propensity score is **inverse probability weighting** (IPW). It does not match; instead it weights each unit by the inverse of the probability of being observed, reweighting a biased sample into a balanced one. The most basic Horvitz-Thompson form targets the ATE:

$$\hat\tau_{HT} = \frac{1}{N} \sum_i \frac{D_i Y_i}{\hat e(X_i)} - \frac{1}{N} \sum_i \frac{(1 - D_i) Y_i}{1 - \hat e(X_i)}.$$

The reason it is unbiased is that $\mathbb{E}[D Y / e(X)] = \mathbb{E}[Y(1)]$ and $\mathbb{E}[(1-D)Y/(1-e(X))] = \mathbb{E}[Y(0)]$ (under the CIA plus overlap): give the rare adopters magnified weight to restore the share of the class of firms they represent. The trouble with Horvitz-Thompson is that the weights can run out of control; when some class of units has propensity very close to 0 or 1, the inverse blows up to astronomical size, and one or two extreme weights dominate the whole estimate, giving enormous variance. The remedy is the **Hajek** (self-normalized) form, dividing by the sum of the realized weights rather than by $N$:

$$\hat\tau_{\text{Hajek}} = \frac{\sum_i D_i Y_i / \hat e(X_i)}{\sum_i D_i / \hat e(X_i)} - \frac{\sum_i (1 - D_i) Y_i / (1 - \hat e(X_i))}{\sum_i (1 - D_i) / (1 - \hat e(X_i))}.$$

Hajek is invariant to an overall shift of the outcome, usually has smaller variance, and is the default choice in practice; but both blow up as the propensity approaches 0 or 1, and the overlap problem returns in the form of extreme weights. When targeting the ATT, the weights are adjusted accordingly: adopters get weight 1, non-adopters are weighted by $e/(1-e)$, reweighting the control group to the adopter group's covariate distribution. In Northwind, Horvitz-Thompson gives 0.217 and Hajek gives 0.230, the latter being steadier, and both land near the ATT.

### 4.5 AIPW and Double Robustness

Regression adjustment, matching, and IPW each bet on one model. Regression adjustment and matching bet that the outcome model $\mu_d(X) = \mathbb{E}[Y \mid D = d, X]$ is right; IPW bets that the propensity model $e(X)$ is right. Either model may be misspecified. **AIPW** (augmented IPW, that is, the doubly robust estimator) combines the two and buys insurance:

$$\hat\tau_{AIPW} = \frac{1}{N} \sum_{i=1}^N \left\{ \hat\mu_1(X_i) - \hat\mu_0(X_i) + \frac{D_i\,(Y_i - \hat\mu_1(X_i))}{\hat e(X_i)} - \frac{(1 - D_i)\,(Y_i - \hat\mu_0(X_i))}{1 - \hat e(X_i)} \right\}.$$

Read it as: first form a regression-adjustment estimate $\hat\mu_1 - \hat\mu_0$ from the outcome model, then use IPW to make one bias correction on the residuals $Y - \hat\mu_d$. Its key property is **double robustness**: as long as one of the outcome model and the propensity model is correctly specified, AIPW is consistent; both need not be right. If both are right, it also attains the semiparametric efficiency bound, and it is the cornerstone of modern double machine learning (paired with cross-fitting). It still needs overlap, and extreme propensities can still make it unstable.

The expression above is written for the ATE and uses both $\hat\mu_1$ and $\hat\mu_0$. As with IPW, the weights are adjusted when targeting the ATT: one only needs to build an outcome model $\hat\mu_0(X)$ for the control group, adopters enter with weight 1, and the non-adopters' augmentation term is weighted by $\hat e/(1 - \hat e)$,

$$\hat\tau_{AIPW}^{\text{ATT}} = \frac{1}{N_1} \sum_{i} \left\{ D_i\big(Y_i - \hat\mu_0(X_i)\big) - (1 - D_i)\,\frac{\hat e(X_i)}{1 - \hat e(X_i)}\big(Y_i - \hat\mu_0(X_i)\big) \right\},$$

where $N_1 = \sum_i D_i$ is the number of adopters. The 0.235 reported in Section 6 is exactly this ATT form: fit only $\hat\mu_0$, weight the controls by $\hat e/(1 - \hat e)$, then divide by the number of adopters.

Lining up this section's estimators along a few dimensions makes clear what each stands on and where each is soft:

| Estimator | Identifying assumption | Model needed | Robustness to misspecification | Sensitivity to weak overlap | Valid standard error |
|---|---|---|---|---|---|
| Difference in means | Randomization | None | Not applicable | Not applicable | Neyman / permutation |
| Regression adjustment | CIA | Outcome model | Low (biased if wrong) | Low | Conventional / robust |
| NN matching (incl. NN-on-PS) | CIA + overlap | Distance metric or propensity | Medium | Medium | Abadie-Imbens (no bootstrap) |
| Kernel matching | CIA + overlap | Propensity | Medium | High | Bootstrap usable |
| IPW Horvitz-Thompson | CIA + overlap | Propensity | Low | Very high | Bootstrap |
| IPW Hajek | CIA + overlap | Propensity | Low | High | Bootstrap |
| AIPW | CIA + overlap | Outcome or propensity (either suffices) | High (doubly robust) | High | Bootstrap / influence function |

The route of this section can be summarized as follows: under randomization, the difference in means plus permutation inference is enough; in observational data, regression adjustment is the least effort but is sensitive to heterogeneous effects and model specification; matching buys weak functional-form assumptions at the cost of the curse of dimensionality and a tricky variance; the propensity score drops the dimension to one; IPW replaces matching with weighting but is hostage to extreme weights, and Hajek is steadier than Horvitz-Thompson; AIPW uses double robustness to buy the insurance of "getting one of two models right is enough," and is the current default first choice. All of these methods share one vital weakness: they are all built on the untestable CIA.

## 5 Anchoring Papers

A method only stands firm once it lands in real research. Two anchoring papers, one a landmark methodological case on "whether observational methods can reproduce an experiment," the other a representative platform-domain work that treats experiments as a touchstone for observational methods, are each laid out along five elements (paper, method, data, results, limitations), with the focus on how the identifying assumption is defended.

### 5.1 LaLonde (1986) and Dehejia-Wahba (1999, 2002)

::: {.case}
Paper and its place in methodological history: LaLonde's "Evaluating the Econometric Evaluations of Training Programs with Experimental Data" (American Economic Review, 1986) is a touchstone for selection-on-observables methods. Its clever move is to take a setting that genuinely has an experimental benchmark and test whether non-experimental econometric methods can reproduce the experimental answer.

Method: LaLonde uses the National Supported Work (NSW) randomized training program, whose experiment itself provides a credible benchmark for the treatment effect. He then replaces the experimental control group with non-experimental comparison groups drawn from the CPS and PSID, estimates the same effect with the era's standard regression and selection-correction methods, and checks whether the experimental benchmark can be recovered. Dehejia and Wahba (Journal of the American Statistical Association 1999; Review of Economics and Statistics 2002) switch, on the same data, to propensity-score methods: first estimate the propensity to adopt, then match within the region of common support.

Data: the NSW experimental sample, plus large non-experimental comparison groups drawn from the CPS and PSID, with the outcome being trainees' subsequent annual earnings.

Results: LaLonde's finding was discouraging: standard non-experimental estimators generally failed to reproduce the experimental benchmark, with sign and magnitude drifting as the comparison group and specification changed. This was a heavy blow to the notion that "controlling for a few variables identifies the causal effect." Dehejia-Wahba brought a ray of hope: they argued that after filtering to a subsample with good common support, propensity-score matching could pull the estimate back near the experimental benchmark (the NSW experimental effect on earnings is around 1,800 dollars).

Limitations and aftermath: Smith and Todd (Journal of Econometrics 2005) pointed out that Dehejia-Wahba's success is quite sensitive to sample selection and specification, and swapping to a different subsample or propensity model may leave the reproduction not holding. This back-and-forth is exactly the theme of this chapter: whether selection on observables can rescue the truth depends on whether the CIA is credible and whether overlap is adequate, and in real data neither can often be guaranteed.
:::

The lesson of this case runs both ways: the propensity score is not a magic that turns observational data into gold, and its success or failure hangs on the untestable CIA and on whether the data truly contain comparable common support. LaLonde taught the whole field to put observational methods on trial using experiments, and this spirit of self-scrutiny is more worth inheriting than any single estimator.

### 5.2 Blake, Nosko and Tadelis (2015)

::: {.case}
Paper: "Consumer Heterogeneity and Paid Search Effectiveness: A Large-Scale Field Experiment" (Econometrica, 2015) is a representative work in platform empirics that uses an experiment to puncture an observational estimate.

Method: the authors, in cooperation with search engines at eBay, ran a large-scale field experiment, turning off paid search ads in some geographic markets while leaving the rest unchanged, using market-level randomization to measure the true incremental effect of advertising. This experimental benchmark is precisely the yardstick used to put on trial the practice of "estimating advertising ROI from observational data," in a line of thinking descended from LaLonde.

Data: eBay's paid-search placements across multiple US geographic markets and the resulting transaction data, with ads split into brand and non-brand keywords.

Results: on the observational reading, paid search looks richly profitable, because most of the people clicking the ad would have come anyway. The experiment peeled back this illusion: the incremental effect of brand-keyword ads was near zero, because those users would have found eBay through organic search even without the ad; the short-run return on non-brand keywords was also far below what the observational estimate implied. The gulf between the observational and experimental benchmarks is, in essence, the selection bias of Section 3 of this chapter, only here the selection is over "who will click the ad."

Limitations: the external validity of a field experiment is constrained by platform, category, and time period, and the conclusion of one experiment need not generalize to all advertisers. But precisely because the observational estimate is so systematically and so greatly wrong here, it powerfully demonstrates that in the selection-riddled platform environment, the default practice of treating correlation as causation can be far from the truth, and credible incremental measurement often demands an experiment.
:::

Put the two together and the point of anchoring is clear: LaLonde and Dehejia-Wahba show that observational methods can, in the best case, approach an experiment, and are unreliable in the general case; Blake-Nosko-Tadelis show that in platform settings selection bias can drag an observational estimate entirely off course, and that experiments are irreplaceable as a benchmark. This also loops back to Northwind's setup: the reason Study A's randomized pilot is the fixed reference point of the whole chapter is precisely that it is immune to the kind of selection both papers warn about.

## 6 A Full Walkthrough on the Northwind Data

Now run the tools of Section 4 on Northwind from start to finish. The code below uses R 4.5.3 with the random seed fixed at 20260711, and every number quoted in the text comes from the actual run output of this code.

### 6.1 The Data-Generating Process

The design is as follows. The two data lines share the same structural effects. The untreated outcome $Y_i(0)$ is composed of a covariate mean function $\mu(X)$, an unobserved firm-quality factor $U_i$ (loading $\phi = 0.30$), and noise; the individual effect $\tau_i$ is composed of a baseline payoff, an observable payoff $g_{obs}(X)$, and $U_i$ (loading $\delta = 0.15$); and $Y_i(1) = Y_i(0) + \tau_i$. The baseline $\tau_{base} = 0.20$ sets the true ATE at about 0.20.

```r
set.seed(20260711)

a0 <- 1.50; b_size <- 0.30; b_dig <- 0.25; b_eff <- 0.20
ind_shift <- c(A = 0.00, B = 0.10, C = 0.20, D = 0.30)
sd_eps <- 0.40; phi <- 0.30
tau_base <- 0.20; gg_dig <- 0.12; gg_size <- 0.10; gg_eff <- -0.05; delta <- 0.15

draw_firms <- function(n) {
  size <- rnorm(n); dig <- rnorm(n); eff <- rnorm(n)
  ind  <- factor(sample(c("A","B","C","D"), n, TRUE), levels = c("A","B","C","D"))
  U    <- rnorm(n); eps <- rnorm(n, 0, sd_eps)
  mu    <- a0 + b_size*size + b_dig*dig + b_eff*eff + ind_shift[as.character(ind)]
  g_obs <- gg_dig*dig + gg_size*size + gg_eff*eff
  tau   <- tau_base + g_obs + delta*U
  Y0 <- as.numeric(mu + phi*U + eps); Y1 <- Y0 + tau
  data.frame(size, dig, eff, ind, U, g_obs, tau, Y0, Y1)
}
observed_Y <- function(df, D) df$Y0 + D * (df$Y1 - df$Y0)
```

The two data lines differ in how they assign. Study A flips a coin, while Study B selects on an index containing $X$ and the observable payoff, and it produces two adoption vectors: one that lets selection fall only on $X$ and $g_{obs}$ ($U$ does not enter selection, so the CIA holds), and one that lets selection additionally attach to the unobserved $U$ (so the CIA fails).

```r
rct <- draw_firms(2000); rct$D <- rbinom(2000, 1, 0.5)
rct$Y <- observed_Y(rct, rct$D)

obs <- draw_firms(6000)
lin <- with(obs, -0.20 + 0.40*size + 0.50*dig + 0.40*eff + 1.50*g_obs)
obs$D_cia  <- rbinom(6000, 1, plogis(lin))              # CIA holds
obs$D_fail <- rbinom(6000, 1, plogis(lin + 1.60*obs$U)) # CIA fails
obs$Y_cia  <- observed_Y(obs, obs$D_cia)
obs$Y_fail <- observed_Y(obs, obs$D_fail)
```

The truth computed from the individual effects: Study B's population ATE is 0.200, the RCT sample's ATE is 0.205; the true ATT of the CIA arm is 0.254, and the true ATT of the CIA-failure arm is 0.311; the two arms' adoption rates are 0.467 and 0.478 respectively. The CIA arm's ATT exceeding its ATE is exactly the work of selection on gains; these truths are the targets against which every estimator below is reconciled. One deliberate design choice bears noting: the selection coefficients are tuned so that the propensity estimates stay in the benign range $[0.018, 0.977]$ with no mass piling up at either end, so that matching and IPW have clean common support to use and are not dragged off by extreme weights.

### 6.2 Study A: Difference in Means and Permutation Inference

```r
y1 <- rct$Y[rct$D == 1]; y0 <- rct$Y[rct$D == 0]
dim_est   <- mean(y1) - mean(y0)
se_neyman <- sqrt(var(y1)/length(y1) + var(y0)/length(y0))
ci <- dim_est + c(-1, 1) * qnorm(0.975) * se_neyman

set.seed(101)
perm <- replicate(5000, {
  Dp <- sample(rct$D)
  mean(rct$Y[Dp == 1]) - mean(rct$Y[Dp == 0])
})
p_perm <- mean(abs(perm) >= abs(dim_est))
```

The adopter group has 1019 firms and the control group 981. The difference in means gives 0.206, the Neyman SE is 0.034, and the 95 percent confidence interval [0.139, 0.274] covers the RCT sample's true value of 0.205 comfortably. Across 5000 permutations, not once was the statistic more extreme than the observed value, so the permutation p-value is below 0.001. This is the credible benchmark from Section 1, which used none of the firms' characteristics and rests entirely on the random assignment itself.

### 6.3 Study B: Naive Comparison and Regression Adjustment

```r
obs$D <- obs$D_cia; obs$Y <- obs$Y_cia
naive <- mean(obs$Y[obs$D==1]) - mean(obs$Y[obs$D==0])

m_add <- lm(Y ~ D + size + dig + eff + ind, data = obs)     # additive, D coefficient is not the ATT

# center all covariates (including industry dummies) at treated-group means so the D coefficient reads as the ATT
Xv <- c("size","dig","eff"); tm <- colMeans(obs[obs$D==1, Xv])
oc <- obs; for (v in Xv) oc[[paste0(v,"_c")]] <- oc[[v]] - tm[v]
ind_d <- model.matrix(~ ind, obs)[, -1]                     # indB, indC, indD
ind_c <- sweep(ind_d, 2, colMeans(ind_d[obs$D==1, ]), "-")  # center the industry dummies too
colnames(ind_c) <- paste0(colnames(ind_d), "_c")
oc <- cbind(oc, ind_c)
m_int <- lm(Y ~ D * (size_c + dig_c + eff_c + indB_c + indC_c + indD_c),
            data = oc)                                       # interaction + full centering -> ATT
```

The naive difference in means is 0.575 (SE 0.019), and subtracting the true ATT of 0.254 leaves a selection bias as large as 0.322, confirming the suspicious number from Section 1. The additive OLS $D$ coefficient is 0.184, a variance-weighted average of the heterogeneous effects that lands near the ATE, not the ATT. Once all covariates (including the industry dummies) are centered at their adopter-group means and $D$ is fully interacted with them, the $D$ coefficient reads an ATT of 0.239 (SE 0.016), correctly returning near the true ATT. The 0.184 and the 0.239 are two different quantities, and this is the confusion most easily tripped over in regression adjustment.

### 6.4 The Propensity Score and Overlap

```r
ps_fit <- glm(D ~ size + dig + eff + ind, data = obs, family = binomial)
obs$e  <- predict(ps_fit, type = "response")
```

The estimated propensity falls in $[0.018, 0.977]$. The figure below plots it as two densities by adoption status. The two distributions are clearly shifted apart, with adopters skewed to the right overall (higher adoption propensity), which is exactly selection at work; but the two supports overlap completely, with no mass piling up at 0 or 1, so overlap is benign and both matching and weighting have clean comparison objects. This overlap figure is a diagnostic you can and should draw, whereas the CIA it stands opposite to cannot draw any figure to vindicate itself, and this contrast is itself worth remembering.

![Propensity-score densities by adoption status, clearly shifted apart but with completely overlapping support.](assets/fig/fig_09_overlap.svg)

### 6.5 Nearest-Neighbor Matching and the Abadie-Imbens Variance

```r
mi   <- matchit(D ~ size + dig + eff + ind, data = obs, method = "nearest",
                distance = "glm", link = "logit", estimand = "ATT", replace = TRUE)
mdat <- match.data(mi)
match_pt <- coef(lm(Y ~ D, data = mdat, weights = weights))["D"]

lps  <- predict(ps_fit)   # linear predictor (logit PS)
mout <- Match(Y = obs$Y, Tr = obs$D, X = lps, estimand = "ATT",
              M = 1, replace = TRUE, ties = TRUE, BiasAdjust = FALSE)
```

Doing 1:1 matching with replacement on the propensity, MatchIt's point estimate is 0.248. Using `Matching::Match` on the same logit PS to take the Abadie-Imbens variance, the point estimate is 0.246 with AI-SE 0.021. Here we match on an estimated propensity, whereas the Abadie and Imbens (2006) variance treats the matching variable as known; Abadie and Imbens (2016) proved that matching on an estimated propensity changes the large-sample variance and usually makes it smaller, so the 0.021 should be read as a somewhat conservative estimate. For contrast, running OLS directly on the matched sample gives an SE of 0.025, which is methodologically wrong (and the bootstrap is outright invalid for NN matching). This time the two SEs happen to be numerically close, but the point is that the AI-SE is the one that stands up, not how much they differ.

### 6.6 IPW and AIPW

```r
ipw_att <- function(Y, D, e) {
  w <- e / (1 - e); m1 <- mean(Y[D == 1])
  c(ht    = m1 - sum((1-D)*w*Y) / sum(D),           # Horvitz-Thompson
    hajek = m1 - sum((1-D)*w*Y) / sum((1-D)*w))     # Hajek
}
ipw_pt <- ipw_att(obs$Y, obs$D, obs$e)

aipw_att <- function(dat, e) {
  m0 <- predict(lm(Y ~ size+dig+eff+ind, data = dat[dat$D==0, ]), newdata = dat)
  p  <- mean(dat$D); w <- e / (1 - e)
  mean((dat$D*(dat$Y-m0) - (1-dat$D)*w*(dat$Y-m0)) / p)
}
aipw_pt <- aipw_att(obs, obs$e)
```

The Horvitz-Thompson version of IPW gives 0.217 (bootstrap SE 0.036), and the Hajek version gives 0.230 (bootstrap SE 0.024). Horvitz-Thompson's standard error is clearly larger, because its weights are not normalized, making it the noisiest of these estimators; Hajek steadies things considerably after normalizing, which is exactly the lesson of Section 4.4, and even though overlap is benign here and extreme weights are not severe, the difference is still visible. AIPW gives 0.235 (bootstrap SE 0.019), the smallest standard error, and the efficiency advantage of double robustness is realized. The IPW and AIPW standard errors from the nonparametric bootstrap are valid, in direct contrast to NN matching where the bootstrap is prohibited.

### 6.7 The Full Reconciliation and the CIA-Failure Variant

Lining up all the estimators to reconcile against the target, the CIA arm's true ATT of 0.254 (with the true ATE of 0.200 as reference):

| Estimator | Estimate | SE | Target |
|---|---|---|---|
| Study A difference in means | 0.206 | 0.034 | ATE (truth 0.205) |
| Study B naive difference | 0.575 | 0.019 | None (biased by 0.322) |
| OLS interacted adjustment | 0.239 | 0.016 | ATT |
| NN matching (AI SE) | 0.246 | 0.021 | ATT |
| IPW Horvitz-Thompson | 0.217 | 0.036 | ATT |
| IPW Hajek | 0.230 | 0.024 | ATT |
| AIPW | 0.235 | 0.019 | ATT |

The conclusion is clean: the naive difference of 0.575 is absurdly high, and every adjusted estimator pulls it back to between 0.22 and 0.25, hovering around the true ATT of 0.254. To be honest about it, these adjusted values all sit slightly below the ATT, landing between the ATE of 0.200 and the ATT of 0.254, and this is finite-sample sampling variation, not bias; we cannot claim any estimator "hits the ATT exactly." In direction and magnitude, the teaching point that "the naive is far above the truth, and after adjustment it is roughly the ATT" holds cleanly.

![Point estimates and 95 percent confidence intervals for each estimator, against the true ATT (dashed) and ATE (dotted). The naive difference sits far to the right, while the randomized pilot and the adjusted estimators cluster near the truth.](assets/fig/fig_09_estimators.svg)

The real alarm bell rings in the last step. Change the selection so that it additionally attaches to the unobserved $U$ (the CIA-failure arm), and the true ATT becomes 0.311. Run the same procedure once more:

```r
dff <- obs; dff$D <- obs$D_fail; dff$Y <- obs$Y_fail   # cleanly swap in the failure-arm treatment and outcome
naive_f <- mean(dff$Y[dff$D==1]) - mean(dff$Y[dff$D==0])
ps <- glm(D ~ size+dig+eff+ind, data = dff, family = binomial)
e  <- predict(ps, type = "response"); w <- e/(1-e)
m0 <- predict(lm(Y ~ size+dig+eff+ind, data = dff[dff$D==0, ]), newdata = dff)
aipw_f <- mean((dff$D*(dff$Y-m0) - (1-dff$D)*w*(dff$Y-m0)) / mean(dff$D))
```

The naive difference is 0.835 (biased by 0.524). This time the same exquisite OLS adjustment reaches only 0.611 (residual bias 0.300), and the doubly robust AIPW reaches only 0.613 (residual bias 0.302). Adjustment shaves off only part of the bias, and the remaining residual of about 0.30 cannot be wiped away no matter what. The reason lies not in the method but in the assumption: $U$ both drives adoption and raises $Y(0)$ yet is not in $X$, so the CIA is broken at the root. AIPW's double robustness can rescue a misspecified model but not a wrong identifying assumption. This residual that cannot be wiped away is exactly the endogeneity that instrumental variables in Chapter 12 will handle.

## 7 Failure Modes and Robustness

In the simulation the identifying assumptions are manufactured, but in real research they can fail at any moment. This section lays out the most common ways they fail and the operational responses.

Inadequate overlap is the most common hard failure of selection-on-observables methods. When some class of units almost all adopts or almost all does not, $e(X)$ approaches 0 or 1, IPW's weights blow up, one or two extreme units dominate the whole estimate, and the variance is enormous and unstable. The diagnosis is direct: plot the grouped density of the propensity and see whether mass piles up at either end and whether the supports overlap. Northwind's overlap was made benign on purpose, and real data are often not so courteous. The response is usually trimming, dropping the units whose $e(X)$ falls in the extreme intervals (say outside $[0.1, 0.9]$) and estimating only in the region of good common support.

But trimming itself overcorrects, and this is a tradeoff to face head-on. Dropping units with extreme propensity quietly changes the estimand: what you estimate is no longer the ATT for all adopters but the ATT for "the portion of adopters in the region of good common support," and the units dropped are often exactly the most typical, most policy-relevant extreme units. The harder you trim, the further the estimand drifts, and while the variance shrinks, the question being answered has also changed. The honest practice is to report how the estimate changes before and after trimming, how many observations were dropped, and to state clearly who the remaining population is.

The matching standard-error trap is worth stressing once more. Under matching with replacement, both the ordinary regression SE and the bootstrap are wrong, the former ignoring the correlation induced by match reuse, the latter inconsistent because the matching functional is not smooth. What stands up is the Abadie-Imbens analytic variance; if matching is on an estimated propensity, it further needs the Abadie-Imbens (2016) correction, which usually gives a smaller variance. In Northwind the two SEs happen to be close, but this must not be taken to mean either one will do, since in real data they can differ by a good deal. Whenever you report the standard error of a matching estimate, first confirm it is the AI variance and not whatever the software spits out by default.

The CIA is untestable, which is the deepest weakness of the whole method, and it must be handled squarely in the argument rather than dodged. All you can do is provide indirect circumstantial evidence: balance tests to see whether the adjusted covariates are balanced across the two groups, placebo tests to see whether a fake treatment that should have no effect estimates zero, and pre-treatment outcome tests to see whether pre-treatment outcomes are already comparable. These are only necessary conditions, and passing them still does not prove the CIA. Northwind takes this to the extreme: the CIA-holding arm and the CIA-failure arm have nearly identical covariate distributions and nearly identical post-adjustment balance, yet one recovers the ATT while the other leaves a residual bias of 0.30. The balance test simply cannot tell the two cases apart, because the difference hides in the unobserved $U$.

Selection on unobservables is exactly the substance of CIA failure. As long as there remains one factor that influences both adoption and outcome yet is unobserved, any adjustment based on observable variables leaves a residual. This is not a problem that a better matching algorithm can solve; it is a failure at the level of identification. The way out is either to find quasi-exogenous variation in treatment at the design level, using instrumental variables to isolate the endogeneity on its own (Chapter 12), or to use panel structure to difference out the time-invariant unobserved heterogeneity (Difference-in-Differences, later in this series). The honest boundary of selection-on-observables methods is to lay out plainly whether the CIA is credible, rather than papering over it with a pile of estimator numbers.

Bad controls are a counterintuitive class of error: more control variables are not better. You must never control for any downstream variable that responds to treatment. In Northwind, if you were to stuff "the additional automation investment a firm makes after adopting the assistant" into the regression as a covariate, it would seem to make the comparison cleaner, but it in fact blocks a channel through which the effect travels, counting part of the real effect as a difference to be controlled away; worse still, if this downstream variable is a collider, controlling for it artificially manufactures a correlation between variables that were unrelated, introducing collider bias. The discipline is a single rule: covariates must be measured pre-treatment. To deal with trend differences or unobserved heterogeneity, the right road is design (IV, panel differencing), not stuffing more controls into the regression.

Stringing these failure modes together, the credibility of selection on observables comes down in the end not to the sophistication of the estimator but to one judgment: whether the factors influencing adoption have truly all been observed. The overlap figure, balance tests, placebo tests, the robustness of trimming, and the check for bad controls are all evidence marshaled around this one question, and none of them can substitute for a substantive judgment about the adoption mechanism. When that judgment does not hold, the right response is not to switch to a fancier estimator but to switch to a different identification route.

## 8 Further Reading

::: {.readings}
Required reading, in the suggested order:

- Imbens and Rubin (2015, *Causal Inference for Statistics, Social, and Biomedical Sciences*). The most authoritative textbook on the potential-outcomes framework, treating unconfoundedness, matching, and subclassification systematically and with restraint, the full backing for Sections 2 through 4 of this chapter.
- Rosenbaum and Rubin (1983, *Biometrika*). The original text of the propensity-score theorem; focus on the definitions of the balancing score and strong ignorability, the source of Section 3.7.
- LaLonde (1986, *American Economic Review*) and Dehejia and Wahba (1999, *JASA*; 2002, *ReStat*). The classic case of whether observational methods can reproduce an experiment, best read together with the critique of Smith and Todd (2005, *Journal of Econometrics*) to see clearly how fragile the CIA is.
- Imbens (2015, *Journal of Human Resources*). A practical guide to the propensity score, with the most down-to-earth treatment of overlap diagnostics, trimming, and the tradeoffs of weighting.

Further reading:

- Abadie and Imbens (2006, *Econometrica*). The large-sample distribution and analytic variance of matching estimators, the theoretical source of the SE trap in Section 4.3.
- Abadie and Imbens (2008, *Econometrica*). The proof that the bootstrap is invalid for nearest-neighbor matching with a fixed number of matches, the original reference for that classic trap.
- Abadie and Imbens (2016, *Econometrica*). The correct large-sample variance when matching on an estimated propensity score, pointing out that the 2006 variance is conservative here.
- Robins, Rotnitzky and Zhao (1994) and Bang and Robins (2005, *Biometrics*). The theoretical origins of AIPW and double robustness; read for the proof of the double-robustness property.
- Blake, Nosko and Tadelis (2015, *Econometrica*). The representative work of using an experiment to puncture an observational estimate in a platform setting, Section 5.2 of this chapter.
- Athey and Imbens (2017, *Journal of Economic Perspectives*) and Imbens (2024 review). A map of modern causal inference, placing matching, weighting, and machine-learning estimators on a single chart; reading it and then returning to this chapter makes each method's position clearer.
:::

::: {.apa-refs}
- Abadie, A., & Imbens, G. W. (2006). Large sample properties of matching estimators for average treatment effects. *Econometrica, 74*(1), 235-267. https://doi.org/10.1111/j.1468-0262.2006.00655.x
- Abadie, A., & Imbens, G. W. (2008). On the failure of the bootstrap for matching estimators. *Econometrica, 76*(6), 1537-1557. https://doi.org/10.3982/ECTA6474
- Abadie, A., & Imbens, G. W. (2016). Matching on the estimated propensity score. *Econometrica, 84*(2), 781-807. https://doi.org/10.3982/ECTA11293
- Athey, S., & Imbens, G. W. (2017). The state of applied econometrics: Causality and policy evaluation. *Journal of Economic Perspectives, 31*(2), 3-32. https://doi.org/10.1257/jep.31.2.3
- Bang, H., & Robins, J. M. (2005). Doubly robust estimation in missing data and causal inference models. *Biometrics, 61*(4), 962-973. https://doi.org/10.1111/j.1541-0420.2005.00377.x
- Blake, T., Nosko, C., & Tadelis, S. (2015). Consumer heterogeneity and paid search effectiveness: A large-scale field experiment. *Econometrica, 83*(1), 155-174. https://doi.org/10.3982/ECTA12423
- Dehejia, R. H., & Wahba, S. (1999). Causal effects in nonexperimental studies: Reevaluating the evaluation of training programs. *Journal of the American Statistical Association, 94*(448), 1053-1062. https://doi.org/10.1080/01621459.1999.10473858
- Dehejia, R. H., & Wahba, S. (2002). Propensity score-matching methods for nonexperimental causal studies. *Review of Economics and Statistics, 84*(1), 151-161. https://doi.org/10.1162/003465302317331982
- Imbens, G. W. (2015). Matching methods in practice: Three examples. *Journal of Human Resources, 50*(2), 373-419. https://doi.org/10.3368/jhr.50.2.373
- Imbens, G. W. (2024). Causal inference in the social sciences. *Annual Review of Statistics and Its Application, 11*(1), 123-152. https://doi.org/10.1146/annurev-statistics-033121-114601
- Imbens, G. W., & Rubin, D. B. (2015). *Causal inference for statistics, social, and biomedical sciences: An introduction*. Cambridge University Press. https://doi.org/10.1017/CBO9781139025751
- LaLonde, R. J. (1986). Evaluating the econometric evaluations of training programs with experimental data. *American Economic Review, 76*(4), 604-620. https://www.jstor.org/stable/1806062
- Robins, J. M., Rotnitzky, A., & Zhao, L. P. (1994). Estimation of regression coefficients when some regressors are not always observed. *Journal of the American Statistical Association, 89*(427), 846-866. https://doi.org/10.1080/01621459.1994.10476818
- Rosenbaum, P. R., & Rubin, D. B. (1983). The central role of the propensity score in observational studies for causal effects. *Biometrika, 70*(1), 41-55. https://doi.org/10.1093/biomet/70.1.41
- Smith, J. A., & Todd, P. E. (2005). Does matching overcome LaLonde's critique of nonexperimental estimators? *Journal of Econometrics, 125*(1-2), 305-353. https://doi.org/10.1016/j.jeconom.2004.04.011
:::
