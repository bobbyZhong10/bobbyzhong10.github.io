---
title: "Digital Experimentation"
subtitle: "Designing A/B Tests That Can Actually Detect What Matters"
seriesline: "Foundations of Information Systems Economics · Chapter 17"
series: "Foundations of Information Systems Economics: Empirical Methods for Information Systems Economics"
sidebartitle: "Foundations of Information Systems Economics"
sidebarsub: "Theory · Methods · Applications"
chapterlabel: "Chapter 17 · Digital Experimentation"
---

## Introduction

Lumen's new feature is assigned by a coin flip, and at the start of the experiment the treatment and control groups balance almost impeccably. Two weeks in, the estimate is positive but not significant. The product manager says "no effect," the analyst says "let it run a few more days," and the dashboard gets checked once a day, until one afternoon the $p$ value finally drops below 0.05 and the team declares victory. Randomization worked the whole time, yet the conclusion lost its promised error rate somewhere in the wait for significance.

Experiments spare you from selection bias, but not from the responsibility of design. When the sample is too small, "not significant" may mean only that the experiment cannot see an effect large enough to matter commercially. When the metric is too noisy, even a flood of traffic gets burned for nothing. And repeated peeking, cherry-picking metrics, and after-the-fact slicing of the population can reliably manufacture good news out of pure noise. Random assignment answers "can we compare without bias," while power, variance, and the stopping rule answer "how much information does this particular comparison carry."

This chapter follows Lumen's decision process. The MDE forces the team to state, before the experiment starts, how small a lift is worth detecting. CUPED uses pre-treatment information to cut noise without stealing any post-treatment variable. Sequential methods reconcile continuous monitoring with a pre-committed error rate. Finally we check a problem that randomization cannot automatically remove: whether treatment and control are fighting over the same pool of delivery capacity or transactions. If they are, the individual-level A/B test still answers some question without bias, only that question need not be "what happens when we roll out to everyone." Chapter 18 takes up that piece.

## 1 A Misread Null Result

::: {.case}
Lumen is a fictional delivery-marketplace platform. The product team has built a new checkout redesign and wants to know whether it lifts user spending. We tell the case with simulated data, whose advantage is that the true effect is fully known, so that the consequence of every design choice can be reconciled against the truth, a teaching condition that a real experiment can never grant. The outcome $Y$ is a user's revenue over the next two weeks. It is noisy: the mean is about 5 currency units but the standard deviation is 6, because a handful of heavy users account for most of the spending. This is exactly the situation Lewis and Reiley faced in advertising experiments, where the signal you are looking for is drowned in enormous user-level noise.

The true effect of the new checkout page is to lift revenue by 2%, that is, to add 0.10 on a baseline of 5. This is a genuine, worth-shipping improvement. The team runs the experiment at a size that sounds perfectly reasonable: 15000 users per arm, treatment on the new page, control on the old, for two weeks.
:::

The experiment finishes and the team runs the most standard analysis: subtract the control group's mean revenue from the treatment group's and run a $t$ test. The difference is 0.059 with a standard error of 0.069, so $t = 0.86$ and $p = 0.392$. That $p$ value is nowhere near the 0.05 line. The product manager takes one look, and the conclusion writes itself: the redesign has no effect, do not ship, move the team onto something else.

The trouble is that we know the truth. The real effect is 0.10, an honest 2% lift, and this experiment missed it completely. Where did things go wrong? Not in identification, since random assignment was perfect. Not in the analysis, since the $t$ test was done correctly. The fault is that the experiment was doomed from the start never to see what it was looking for. With 15000 per arm and noise this large, the minimum effect the experiment can detect with 80% power (the MDE) is 0.194, a 3.9% lift, nearly twice the true effect. The real 2% effect sits well below that threshold, and the experiment's power against it is only 0.303. In other words, even if the effect is entirely real, this experiment has only a three-in-ten chance of returning a significant result. Verifying by simulation, rerunning this design many times, indeed only 30.6% of the experiments reject the null, agreeing tightly with the theoretical 0.303.

This is the first, and the most common, error the chapter guards against: reading the null result of an underpowered experiment as evidence of "no effect." $p = 0.392$ does not say the redesign is useless, it says the experiment was too small to see clearly. To detect the true 2% effect at 80% power, each arm needs about 56500 users, nearly four times the actual investment. Section 3 will derive the formulas for MDE and sample size so that this kind of accident can be computed away before the experiment even starts. Section 6 will show that this same 15000-user dataset, once run through the CUPED of Section 3.2, has its noise cut in half, letting the missed effect reappear.

## 2 The Economic Model and the Estimand

The nice thing about a randomized experiment is that the estimand and the identification are so simple they barely need saying, but precisely because they are simple, it is easy to rush past the parts that should be said clearly. This section sets the language straight, with the emphasis on the variance of the experimental estimator, because all of the design effort in the rest of the chapter is aimed at that variance.

We keep the potential outcomes of Chapter 9. User $i$ has two potential revenues: $Y_i(1)$ is the revenue realized if assigned to the new page, $Y_i(0)$ is the revenue under the old page, and the individual effect is $\tau_i = Y_i(1) - Y_i(0)$. The experiment splits users at random into two groups, the treatment indicator $D_i \in \{0, 1\}$ is decided by a coin flip, and the observed revenue is $Y_i = D_i Y_i(1) + (1 - D_i) Y_i(0)$. The target estimand is the average treatment effect:

$$\tau = \mathbb{E}[Y_i(1) - Y_i(0)].$$

Identification is free, and this is exactly where the entire advantage of an experiment over an observational study lies. Because $D_i$ is randomly assigned, it is independent of the potential outcomes, $\big(Y_i(0), Y_i(1)\big) \perp D_i$, and so the difference in the two groups' mean observed revenue estimates $\tau$ without bias:

$$\hat\tau = \frac{1}{n_1}\sum_{i: D_i = 1} Y_i - \frac{1}{n_0}\sum_{i: D_i = 0} Y_i.$$

You need not control for any covariate, you need no assumption about confounding, randomization zeroes out selection bias in one stroke. This is the heart of design-based inference (Chapter 9): the causal conclusion rests on the random assignment that the experimenter controls, not on modeling assumptions about the world.

What really matters is the variance of this estimator, because it decides whether the experiment sees clearly or not. Under random assignment, the variance of $\hat\tau$ is

$$\text{Var}(\hat\tau) = \frac{\sigma_1^2}{n_1} + \frac{\sigma_0^2}{n_0},$$

where $\sigma_1^2$ and $\sigma_0^2$ are the variances of the potential outcomes in the two groups (Neyman's classic result; the version that estimates them by within-group sample variances is slightly conservative under heterogeneous effects, an honest direction to err in). When the two arms are equal in size and similar in variance, it simplifies to $2\sigma^2 / n$, where $n$ is the per-arm sample size. The formula is frighteningly plain, yet it is the axis of the chapter: the $\sigma^2$ in the numerator is the noise in the outcome, and the $n$ in the denominator is the sample size. Whether the experiment can detect an effect turns entirely on how large the effect $\tau$ is relative to the square root of this variance.

So there are only two ways to raise an experiment's sensitivity, and the design of the whole chapter travels along these two: either grow $n$ (expensive, and with linearly diminishing returns, since the variance falls only as $1/n$), or shrink $\sigma^2$ (this is where CUPED comes in, and it can halve the noise without adding a single user). The failed experiment of Section 1 was sick precisely because $2\sigma^2/n$ was too large and $\tau$ too small, so the signal-to-noise ratio $\tau / \sqrt{2\sigma^2/n}$ was too weak for the test to reject.

To summarize: identification in a randomized experiment is free, and the difference-in-means estimates the ATE without bias and without any modeling assumption. All of the technical tension in the chapter concentrates on the variance $2\sigma^2/n$ of the estimator, and lowering it runs along two roads, adding sample and cutting variance, which the later sections take up in turn.

## 3 Identification and Design: From Power to Interference

This is the most finely analyzed section of the chapter. Unlike the earlier chapters, it does not dissect identifying assumptions (in an experiment identification is free), it dissects design: how to compute, before the experiment starts, what it can see (power and MDE), how to cut noise so it sees more clearly (CUPED), how to hold the test level under continuous monitoring (sequential testing), and how to estimate the right thing even when units interfere with one another (cluster and switchback). The logic moves from shallow to deep, each subsection solving a new problem exposed by the one before.

### 3.1 Power and MDE: The Sums You Should Do Before the Experiment Starts

The accident of Section 1 could have been avoided before the experiment ever started, just by doing some arithmetic. That arithmetic is the power analysis, and its product is the MDE, the minimum effect the current design can reliably detect.

First let us get the logic straight. We want to test the null $H_0: \tau = 0$ against the alternative $\tau \neq 0$, at level $\alpha$ (0.05 by convention), rejecting when $|\hat\tau / \text{SE}| > z_{1 - \alpha/2}$. Power is the probability $1 - \beta$ (0.80 by convention) that the test rejects the null when the true effect is in fact some $\tau$. Since $\hat\tau$ is approximately normal, centered at $\tau$ with standard deviation SE, power is approximately

$$\text{power} \approx \Phi\!\left(\frac{\tau}{\text{SE}} - z_{1 - \alpha/2}\right).$$

The formula says it plainly: for power to be high, the effect must be large enough relative to the standard error. Inverting it, setting power exactly to $1 - \beta$, gives the smallest effect detectable at that confidence, the MDE:

$$\text{MDE} = (z_{1 - \alpha/2} + z_{1 - \beta}) \cdot \text{SE} = (z_{1 - \alpha/2} + z_{1 - \beta}) \cdot \sigma\sqrt{\frac{2}{n}}.$$

At $\alpha = 0.05$ and power $= 0.80$, we have $z_{1 - \alpha/2} + z_{1 - \beta} = 1.96 + 0.84 = 2.80$. In Section 1, $\sigma = 6$ and $n = 15000$, so SE $= 6\sqrt{2/15000} = 0.0693$ and MDE $= 2.80 \times 0.0693 = 0.194$, exactly the number behind "the smallest detectable effect is nearly twice the true one." Turning it around, replace the MDE with the effect $\tau$ you actually care about and solve for the required sample size:

$$n = \frac{2\sigma^2 (z_{1 - \alpha/2} + z_{1 - \beta})^2}{\tau^2}.$$

To detect Lumen's $\tau = 0.10$, plug $\sigma \approx 6$ into this, and each arm needs about 56500, nearly four times the actual investment. This formula should be run once at the design stage, because it turns "is this experiment worth running" into a question that can be answered before the experiment starts.

Three implications are worth spelling out. First, sample size grows as $1/\tau^2$, so halving the effect quadruples the required sample, which is why the cost of detecting small effects rises so steeply, and why Lewis and Reiley lamented that measuring an advertising effect can take hundreds of millions of observations. Second, MDE falls only slowly, as $1/\sqrt{n}$, so piling on sample to buy sensitivity has diminishing returns, which pushes the incentive toward cutting variance. Third, a power analysis must be done before the experiment starts, since taking a nonsignificant result after the fact and computing an "observed power" is circular reasoning and utterly meaningless; the right question is always "given my sample and my noise, how large an effect could I have seen."

### 3.2 Variance Is the Enemy: CUPED

The MDE formula points its finger at $\sigma^2$: cut the variance of the outcome and you lower the MDE without adding a single user. CUPED (controlled experiment using pre-experiment data, Deng, Xu, Kohavi, and Walker 2013) is the variance-reduction tool most used across platforms, and its idea is elegant to the point of feeling like a free lunch.

The key insight is that a large part of a user's revenue during the experiment can be predicted from before the experiment ever began: a user who spent heavily in the previous two weeks probably spends heavily in these two weeks as well. That predictable variation has nothing to do with treatment, it is pure noise, and if we can subtract it off, the remaining variance shrinks. Let $X_i$ be user $i$'s pre-experiment revenue (a covariate known before the experiment begins). CUPED uses it to build an adjusted outcome:

$$Y_i^{cv} = Y_i - \theta (X_i - \bar X), \qquad \theta = \frac{\text{Cov}(Y, X)}{\text{Var}(X)}.$$

Here $\theta$ is exactly the slope of the regression of $Y$ on $X$, and $Y^{cv}$ is the residual of $Y$ after subtracting the part linearly predictable from $X$ (with a constant $\theta \bar X$ added back to keep the mean unchanged). Two properties make it work. First, unbiasedness is untouched: because random assignment guarantees that $X$ has the same distribution in both arms, $\mathbb{E}[X - \bar X]$ is zero in each, so the difference-in-means computed on $Y^{cv}$ still estimates $\tau$ without bias. Second, a chunk of variance is cut away:

$$\text{Var}(Y^{cv}) = \text{Var}(Y)\,(1 - \rho^2), \qquad \rho = \text{corr}(Y, X).$$

The fraction of variance removed is exactly $\rho^2$, the square of the correlation between pre-experiment revenue and experiment-period revenue. In Lumen this correlation is about 0.70, so $\rho^2 \approx 0.49$, and nearly half the variance is cut, equivalently the variance drops to $1 - \rho^2 \approx 0.51$ times the original. Because sample size and variance enter the MDE formula equivalently, halving the variance is like doubling the sample for free: to reach the same power with CUPED you need only about 29000 per arm rather than about 56500, nearly a halving.

::: {.intuition}
CUPED is really something everyone has known for a long time, regression adjustment, which statisticians call ANCOVA. Regress $Y$ jointly on the treatment indicator $D$ and the pre-experiment covariate $X$, and the coefficient on $D$ is the CUPED estimate. Adding $X$ absorbs the predictable part of the outcome's variation, so the coefficient on $D$ is estimated more precisely. Put differently, CUPED is not new magic, it is the old truth that "controlling for a strong predictor improves precision," applied to a randomized experiment. It shares the spirit of using the pre-period as a baseline in the DiD of Chapter 13, both cushioning the experiment-period comparison with information from before the experiment.
:::

One discipline must be kept: the covariate $X$ must take a pre-experiment value. This is CUPED's one red line, and crossing it turns variance reduction into bias. If you adjust with a variable measured during or after the experiment (say a user's click count during the experiment), that variable may itself be affected by treatment, and subtracting it subtracts part of the treatment effect too, biasing the estimate toward zero. This is the same mistake as the bad control of Chapter 13. Safe covariates can come only from the world before the experiment started, the world the treatment has not yet touched. Section 6 will show that the same 15000-user data that gave $p = 0.392$, once adjusted by CUPED, has its $p$ fall to 0.032, so the missed effect reappears, all without recruiting a single extra user.

### 3.3 Peeking and Sequential Testing

Both power and CUPED assume the experiment runs to a preplanned sample size and the result is looked at once, at the end. Real platforms do not work that way: the experiment sits on a live dashboard, the analyst refreshes it daily or even hourly, and the moment it looks significant there is an urge to stop and ship. This continuous monitoring combined with "stop as soon as it is significant" is called peeking, and it sends the false positive rate completely out of control.

First see how bad the problem is. Suppose the new checkout page actually has no effect at all (an A/A test, true $\tau = 0$), the experiment runs two weeks, the result is checked once a day, and significance is declared as soon as any day shows $p < 0.05$. On any single day the type I error is indeed 5%, but over fourteen days, stopping the moment any one day crosses the line, the cumulative false positive rate is far above 5%. Computing this by simulation: one look is 4.9%, two looks rise to 8.2%, five looks 14.5%, ten looks 19.3%, and all fourteen looks reach 21.3%. In other words, a redesign with no effect at all, allowed to be peeked at daily, has better than a one-in-five chance of being "found" to have a significant effect. This is not analyst dishonesty, it is that a fixed-threshold test was never designed for repeated looking, and the more you look, the more likely you are to run into an accidental significance.

::: {.warning}
The harm of peeking goes beyond the inflated false positive rate, it also poisons everything that comes after. Once you stop on some day because $p < 0.05$, the stopping rule itself makes your sample a selected one, and the effect estimate reported at that moment will systematically overstate (the day you cross the line tends to be the day noise happened to favor you), and the confidence interval is distorted too. So "peeking just lets me find a real effect sooner" is self-deception, it also makes you more likely to find a false effect and to overstate a real one.
:::

The remedy is the idea of always-valid inference: design a test whose type I error does not exceed $\alpha$ at any stopping time, no matter how many times you peek. Johari, Koomen, Pekelis, and Walsh brought this to A/B testing, with the mixture sequential probability ratio test (mSPRT) as the core tool. It places a prior over the effect size (the mixing distribution), builds a likelihood-ratio statistic $\Lambda_n$, and rejects only when $\Lambda_n \geq 1/\alpha$. The key property is that under the null $\Lambda_n$ is a martingale, so by the optional stopping theorem, whenever you stop, $\Pr(\exists n: \Lambda_n \geq 1/\alpha) \leq \alpha$. The analyst can therefore look every day and stop at any time with a clear conscience, and the test level holds. In the simulation of Section 6, monitoring the same fourteen times continuously, mSPRT's type I error is only 0.007, comfortably below 0.05, while the naive approach is at 0.213.

There is no free lunch, and this look-anytime guarantee is bought with power. To be robust to any stopping time, mSPRT is more conservative than a fixed-sample test that looks once at the end, and against the same small effect it has lower power (in Section 6, mSPRT's power against the true effect is 0.312, while the fixed-sample test at the same terminal sample size has 0.654). This cost is especially sharp for small effects, because a small effect can only be seen once the full sample is collected, and the sequential method rarely gets to stop early. The real value of sequential design comes when the effect is large: when the effect is obvious $\Lambda_n$ crosses the line quickly, and the experiment can end substantially early, saving time and traffic. Besides mSPRT, the other mainstream route is group-sequential methods, which use the O'Brien-Fleming or Pocock alpha-spending boundaries to allocate the total type I error budget across a preset number of checkpoints, strict thresholds early and looser ones late, likewise holding the total level while allowing stopping partway. Both routes share one principle: to look anytime you must pay for the number of looks with test level or power, while naively reusing a fixed-threshold test simply gives the false positive rate away.

### 3.4 Multiple Testing

Peeking is repeated testing along the time dimension, and its twin problem is repeated testing along the metric dimension. An experiment often watches dozens of metrics at once (revenue, retention, clicks, time on site), or compares several variants at once, and declares success if any one is significant. The logic is exactly the same as peeking: the more tests you run, the higher the chance at least one is significant by accident, and the family-wise error rate (FWER, the overall type I error of a family of tests) balloons with the number of tests. Testing 20 metrics that are all truly null, the chance at least one is significant at the 0.05 level is close to $1 - 0.95^{20} \approx 0.64$.

The remedy comes in two kinds, depending on what you fear. If you must control the FWER (not a single false positive can be tolerated, as in a ship decision), the Bonferroni correction lowers each test's threshold to $\alpha / m$ ($m$ is the number of tests), guaranteeing the overall FWER stays below $\alpha$, at the cost of being conservative and low in power. If you can tolerate a few false positives and care more about the overall quality of discoveries (say, exploratory screening for which metrics moved), then controlling the false discovery rate (FDR, the expected fraction of false positives among the rejected hypotheses) fits better, and the Benjamini-Hochberg procedure, which sorts the $p$ values and applies increasing thresholds, is far more powerful than Bonferroni. The common practice on platforms is to layer the metrics, looking first at a pre-specified primary metric for the ship decision, and managing the rest as guardrail and exploratory metrics separately, avoiding by design the fishing of "report whichever is significant."

### 3.5 SUTVA and Interference

This subsection is responsible only for identifying the problem and building the interface to Chapter 18, and does not fully develop the interference estimands here. Chapter 18 will rewrite this intuition as exposure-dependent potential outcomes and distinguish the direct, indirect, total, and overall effects.

The variance and testing problems of the earlier sections all assume something that holds by default but is often ignored: a user's outcome is determined only by which arm that user is assigned to, and is unaffected by anyone else's assignment. This is the SUTVA (stable unit treatment value assumption) of Chapter 13, whose first clause is no interference. On many platforms, especially marketplace platforms, this clause systematically fails, and once it fails, the perfect identification of randomization cannot save the estimate.

Lumen's setting lays the problem out plainly. It is a delivery marketplace, and delivery capacity (couriers) is a limited shared resource. Suppose the new checkout page raises users' willingness to place orders. In a user-level A/B, treatment and control users are mixed in the same pool of couriers, the same market, competing for capacity. Treatment users order more aggressively, consume more capacity, and crowd out the couriers who should have served control, so control users' experience is dragged down and their spending is suppressed. The result is that treatment looks especially good because control was made to look worse, and the difference between the groups folds this see-saw (cannibalization) into the effect.

To see the direction and size of the bias, we first have to be clear about the estimand we actually want. Management cares about the global average treatment effect: if the new page were given to all users versus the old page to all users, how much does average spending differ. This is the quantity a ship decision needs. But the user-level A/B does not estimate this. It estimates how much higher treatment is than control in that particular world where half the population is treated, and in that contrast control is cannibalized by treatment while treatment enjoys the capacity control gave up, so both sides are biased. In the simulation of Section 6, Lumen's true global effect is 0.040, while the user-level A/B gives 0.129, overstating by a factor of 3.2. The direction is systematic: for demand-side interventions that worsen congestion, the user-level experiment usually overstates, because after a full rollout the congestion would be far worse than under a half rollout, and the experiment only experienced the mild congestion of a half rollout.

::: {.warning}
The insidiousness of interference is that it does not violate randomization, the random assignment is still perfect, the $p$ value is still pretty, yet SUTVA has quietly broken underneath. No statistical test can find it, because it is not a sampling problem but an estimand problem: you have precisely estimated a quantity you did not want. Judging whether you are contaminated by interference relies on economic judgment about the mechanism through which treatment propagates. Shared capacity, shared inventory, shared attention budget, social-network spillovers, general-equilibrium feedback through prices or auctions are all channels of interference, and they are everywhere on platforms. This is the same ghost as the SUTVA failure in the DiD of Chapter 13, only in an experiment it is more easily hidden behind the illusion that "I randomized, so all is well."
:::

### 3.6 Cluster and Switchback Designs

The cluster and switchback here serve only as a design preview, explaining why changing the granularity of randomization may come closer to the rollout effect. Assignment probabilities, two-stage saturation, carryover, exposure mapping, and formal inference are all consolidated in Chapter 18.

The root of interference is splitting into different treatment arms the units that should have moved together. The fix follows naturally: switch the granularity of randomization so that mutually interfering units land in the same treatment state, and each unit thereby experiences a globally consistent world. Two mainstream designs correspond to two granularities.

Cluster randomization assigns by whole markets: an entire city or region is assigned to treatment or control, instead of splitting the users inside it. This way, all users within a market use either the new page or the old, the cannibalization of capacity happens within a single treatment state, and no longer contaminates across arms. As long as markets do not interfere with one another, the cluster experiment estimates the global effect without bias. The cost is variance: the effective sample size collapses from the number of users to the number of markets, and an experiment with tens of millions of users may be left with only a few dozen independent markets, so the standard error is much larger. This is the fundamental interference-variance tradeoff, cluster buying unbiasedness with variance.

Switchback design randomizes by time blocks: it flips the whole market back and forth between the new and old page over a sequence of time blocks, comparing the average outcome in treated blocks against control blocks. Within each time block the entire platform is in a single treatment state, and the congestion in the block is the true congestion of that state, so comparing treated blocks against control blocks approaches the global effect. Switchback suits resources like delivery capacity that clear quickly on short time scales especially well, since one block's capacity state does not readily bleed into the next. Its risk is along the time dimension: if treatment has a lagged effect across blocks (carryover), or there is strong time autocorrelation, adjacent blocks are no longer clean, and you must leave buffers between blocks or build the time autocorrelation into the variance estimate. In the simulation of Section 6, the user-level A/B overstates by a factor of 3.2, while switchback and cluster both land almost exactly on the true global effect of 0.040, at the cost of sampling variability markedly larger than the wildly biased but very stable user-level estimate, the interference-variance tradeoff laid out in numbers.

The main points of this section are as follows: identification in an experiment is free, but design is not free. Power and MDE let you compute before the experiment starts how large an effect you can see, avoiding the misreading of underpowered as no effect. CUPED uses pre-experiment information to cut variance, equivalent to doubling the sample for free. Peeking inflates false positives along time and multiple testing along metrics, and the sequential always-valid methods and the FWER/FDR corrections remedy them respectively, at the cost of power. The interference of SUTVA makes the user-level A/B precisely estimate the wrong quantity, and cluster and switchback cage the interference with a different randomization granularity, trading variance for unbiasedness.

## 4 Estimation and Design Choices

The ideas are laid out, and this section is about putting them into practice: exactly how each design is computed, how it is implemented, and how they trade off against one another.

### 4.1 The Practice of Power and Sample Size

The pre-experiment power analysis, in practice, is just plugging your scenario into the formulas of Section 3.1. It takes three inputs: the minimum effect $\tau$ you care about (from business judgment, how small a lift is worth shipping), the standard deviation $\sigma$ of the outcome (from historical data), and the $\alpha$ and power targets. Plug into $n = 2\sigma^2(z_{1-\alpha/2} + z_{1-\beta})^2 / \tau^2$ for the per-arm sample size. If the outcome is binary (like a conversion rate $p$), plug in $\sigma^2 = p(1-p)$ directly. If the outcome is heavy-tailed (like revenue), be wary when estimating $\sigma$ from historical data that extreme values do not blow up the variance, which is precisely why the CUPED of the next section is worth doing. Conversely, if the sample size is constrained by traffic and cannot be freely enlarged, plug it into the MDE formula to get "the largest effect this experiment can reliably see," and use that to judge whether it is worth running, rather than running an experiment doomed to see nothing.

### 4.2 Implementing CUPED

Putting CUPED into practice takes only two steps. First, take from the pre-experiment data a covariate $X$ strongly correlated with the outcome (the pre-experiment value of the same metric is usually best, say using the prior two weeks' revenue to predict the experiment-period two weeks' revenue), and estimate $\theta = \text{Cov}(Y, X) / \text{Var}(X)$ on the pooled sample. Second, construct $Y^{cv} = Y - \theta(X - \bar X)$, then do the usual difference-in-means and $t$ test on $Y^{cv}$. Equivalently, regress $Y$ directly on $D$ and $X$ by OLS and read the coefficient on $D$, and this ANCOVA form also hands you the correct standard error. Two practical notes: $X$ can only take a pre-experiment value (the red line of Section 3.2), and you may use several pre-experiment covariates, or even use machine learning to combine them into a single predicted value $\hat Y$ and then treat that as $X$ (this is the generalization of CUPED, called CUPAC), and the stronger the correlation, the more variance is cut.

### 4.3 Implementing Sequential Testing

To look repeatedly during the experiment without breaking the test level, there are two mature routes. The mSPRT route: set a prior variance $\tau_0^2$ for the effect (prior $N(0, \tau_0^2)$, usually set to your expected effect magnitude, and it affects only power, not level), compute the mixture likelihood ratio $\Lambda_n$ at each checkpoint, stop and declare significance only when $\Lambda_n \geq 1/\alpha$, and the corresponding always-valid $p$ value is the running minimum of $1/\Lambda_n$. The group-sequential route: fix the number of checkpoints in advance, and use O'Brien-Fleming or Pocock boundaries to allocate the total $\alpha$ budget, where O'Brien-Fleming is extremely strict early and close to the usual threshold late, suiting cases that want to preserve late-stage power, while Pocock uses one (stricter) threshold at every checkpoint, suiting cases that want to stop as early as possible. Both have mature implementations (such as gsDesign), and in practice you choose by the checking frequency the business can tolerate and its preference for early stopping. The core is that the checkpoints and the $\alpha$ allocation must be fixed before the experiment starts, and cannot be added as you go.

### 4.4 Estimating Switchback and Cluster

The estimation of interference-robust designs turns on the variance, not the point estimate. The cluster design's point estimate is just the market-level difference-in-means, but the standard error must be clustered at the market level (the reasoning of Chapter 13), because the independent units are markets not users, and a user-level standard error would badly understate the uncertainty and produce spurious significance. The effective sample size is the number of markets, and when markets are few you also need the few-clusters first-aid kit of Chapter 13 (CR2, the wild cluster bootstrap). The switchback design's point estimate is the mean difference between treated and control blocks, and the variance estimate must confront time autocorrelation: adjacent blocks' outcomes are usually positively correlated, the independence assumption fails, and you need a block-level or autocorrelation-robust variance estimate. Bojinov, Simchi-Levi, and Zhao (2023) give the optimal design and the corresponding inference theory for switchbacks. Both designs trade sample size (the number of independent units) for robustness to interference, and when reporting you should honestly disclose this variance cost, rather than applying a user-level formula that understates the uncertainty.

The practical points of this section are as follows: a power analysis takes three inputs available before the experiment starts (the minimum effect you care about, the outcome standard deviation, the test level and target power), plug them into the formula for the sample size or MDE. CUPED needs only one pre-experiment covariate and a single residualization, equivalent to running ANCOVA on $D$ and $X$. Sequential testing uses mSPRT or group-sequential boundaries to hold the level while allowing stopping partway, and the checkpoints must be fixed before the experiment starts. The difficulty of interference-robust designs lies entirely in the variance, cluster must cluster to the market level and switchback must be robust to time autocorrelation, both trading the number of independent units for unbiasedness.

## 5 Anchoring Papers

Methods only stand up when grounded in real research. Three anchoring papers, one that pushes the pain of underpower to the extreme, one that founds CUPED, and one that exposes interference in marketplace experiments, each organized by the five elements of paper, method, data, results, and limitations, with the focus on how the design answers the difficulty.

### 5.1 Lewis and Reiley (2014)

::: {.case}
Paper: "Online Ads and Offline Sales: Measuring the Effect of Retail Advertising via a Controlled Experiment on Yahoo!", Quantitative Marketing and Economics. It is the textbook case of the predicament "identification is perfect, yet the effect is almost impossible to measure because the noise is too large," and the opening puzzle of Section 1 is its miniature.

Method: a clean large-scale randomized experiment. The authors worked at Yahoo! and matched Yahoo!'s user database to a large brick-and-mortar and online retailer, identifying about 1.5 million consumers, split 80/20 into treatment and control, with the treatment group shown the retailer's display ads while browsing Yahoo!, after which both groups' online plus offline purchases were tracked. Random assignment made ad exposure exogenous, and the identification itself is airtight.

Data: individual-level purchase records for about 1.5 million matched consumers, with per-person sales during the treatment period of about 1.89 (currency units) but a standard deviation as high as 19, and ad spend about 1% of the target consumers' sales.

Results: the core conclusion is less an effect value than an alarm bell. Individual sales are enormously noisy, and the advertising ROI is a tiny slice of that sales figure, so to detect a 10% ROI at significance the required standard error is about one ten-thousandth of the sales standard deviation, which translates into a sample of hundreds of millions of consumers. The sample of 1.5 million, huge as it sounds, is nowhere near enough for this signal-to-noise ratio, and the intent-to-treat (ITT) estimate is drowned in noise and not significant.

Limitations: the authors honestly attribute the nonsignificance to a power problem rather than to advertising being ineffective. They also show that a naive observational comparison of users who "saw the ad" against those who "did not" is badly distorted: in their data such a naive comparison gives an estimate opposite in sign and several times larger in magnitude, because whether a user saw the ad depends on browsing behavior and is by no means random. This is exactly why one must run an experiment, and why the experiment must also be large enough.
:::

The value of this paper is that it turns the $2\sigma^2/n$ of Section 2 into a tale of blood and tears: an experiment with perfect identification can still measure nothing because $\sigma$ is too large and $n$ too small, and a power analysis is not a formality, it is the first gate that decides whether the experiment lives or dies.

### 5.2 Deng, Xu, Kohavi and Walker (2013)

::: {.case}
Paper: "Improving the Sensitivity of Online Controlled Experiments by Utilizing Pre-Experiment Data", Proceedings of the Sixth ACM International Conference on Web Search and Data Mining (WSDM). It proposes CUPED, the standard tool across platforms for cutting variance and raising experimental sensitivity, and the source of the method in Section 3.2.

Method: the core is to build a control variable from pre-experiment data to cut the outcome variance. For each user take a pre-experiment covariate $X$ (the most effective is the pre-experiment value of the same metric), construct $Y^{cv} = Y - \theta(X - \bar X)$ with $\theta$ chosen to minimize variance. Because random assignment guarantees $X$ has the same distribution in both arms, the adjustment introduces no bias, yet the variance is cut by $1 - \rho^2$. The paper also discusses generalizations using multiple covariates and stratification.

Data: real A/B tests on Microsoft's Bing search engine's large-scale online experimentation platform, with metrics including various user engagement measures.

Results: in real experiments, using the same pre-experiment metric for CUPED usually cuts about half the variance, equivalent to halving the experiment time or sample size to reach the same sensitivity. This "double the sample for free" gain made CUPED quickly become a default component of the experimentation systems at major platforms.

Limitations: the gain depends entirely on how strongly the pre-experiment covariate correlates with the outcome, and a weak correlation gives limited reduction; for new users, or settings with no pre-experiment data, there is nothing to work with; and the covariate must come strictly from before the experiment, otherwise it introduces bias.
:::

This paper carries through the road of "lowering $\sigma^2$" from Section 2, and its beauty is that without recruiting a single extra user and without changing unbiasedness, it doubles the experiment's sensitivity just by using a pre-experiment record that was sitting in the database all along.

### 5.3 Interference in Marketplace Experiments

::: {.case}
Paper: Blake and Coey (2014), "Why Marketplace Experimentation Is Harder Than It Seems: The Role of Test-Control Interference", Proceedings of the Fifteenth ACM Conference on Economics and Computation (EC). It uses a platform experiment to expose the systematic bias of the user-level A/B on marketplace platforms, and it is the empirical support for Sections 3.5 and 3.6; the methodological counterpart is Bojinov, Simchi-Levi, and Zhao (2023) on the theory of switchback design.

Method: the core argument is that on a marketplace platform the treatment and control groups interfere with each other through the shared market mechanism, making the standard A/B estimate deviate from the global effect. Taking an e-commerce platform as the example, if treatment raises the bids or activity of treatment-group buyers, they crowd out control-group buyers in the auction or the matching, so control is suppressed, and the difference between the groups overstates the effect of a full rollout. The paper uses real platform experiments to show the magnitude of this interference, and discusses using market-level (cluster) randomization to mitigate it. Bojinov et al., for their part, give a formal framework, optimal design, and inference for switchback design.

Data: real experiment data on a large e-commerce/marketplace platform (Blake-Coey use a marketplace platform of the eBay kind), comparing estimates under user-level versus market-level randomization.

Results: the effect magnitudes given by the user-level A/B and the market-level design differ systematically, the former deviating from the global effect because of test-control interference; switching to a cluster or switchback design brings the estimate closer to the true impact of a full rollout, but the variance grows substantially. This "trade bias for variance" tradeoff is the inescapable core tension of marketplace experiments.

Limitations: the effective sample size of a cluster design collapses to the number of markets, so the variance is large and it is especially hard on platforms with few markets; switchback relies on the treatment effect clearing quickly in time with no long lag, otherwise cross-block contamination reappears; and judging whether a given intervention is subject to interference, and how much, still needs a substantive judgment about the market mechanism, with no once-and-for-all test.
:::

Taken together, the three make the meaning of the anchoring clear: Lewis and Reiley show that an experiment with perfect identification can still lose to noise, so power is the first gate; Deng et al. show how to use pre-experiment data to halve the noise and double the sensitivity; Blake and Coey, and Bojinov et al., reveal that choosing the wrong randomization granularity makes the experiment precisely estimate the wrong quantity, and interference must be treated with design rather than statistics. Their common lesson is that the success or failure of an experiment is mostly decided in the design before it ever starts.

## 6 A Full Walkthrough on the Lumen Data

Now we run the earlier tools on Lumen from start to finish. The code below uses R 4.5.3 with a fixed random seed for reproducibility, and every number quoted in the text comes from the actual output of running this code. The basic setup: user two-week revenue has mean 5 and standard deviation 6, the correlation between pre-experiment revenue and experiment-period revenue is 0.70, and the true effect of the new checkout page is 0.10 (a 2% lift).

### 6.1 The DGP and the Opening Null Result

The design parameters are as follows: each user has a persistent value $v$, and both the pre-experiment revenue $X$ and the experiment-period control revenue $Y(0)$ equal $5 + 5v$ plus independent noise, so the two correlate at about 0.70; treatment adds 0.10 to revenue.

```r
gen_users <- function(n) {
  v  <- rnorm(n)
  X  <- 5 + 5 * v + rnorm(n, 0, 3.317)   # pre-period revenue
  Y0 <- 5 + 5 * v + rnorm(n, 0, 3.317)   # control potential
  Y1 <- Y0 + 0.10                        # treated potential (2% lift)
  data.frame(X = X, Y0 = Y0, Y1 = Y1)
}
```

First reproduce the opening. With 15000 users per arm, random assignment, and a difference-in-means $t$ test: the difference is 0.059, the standard error 0.069, $t = 0.86$, $p = 0.392$, not significant. Yet the true effect is 0.10. This design's 80%-power MDE at $n = 15000$ is 0.194 (a 3.9% lift), and its power against the true 2% effect is only 0.303, and verifying by 2000 simulation repetitions, the measured rejection rate is 0.306, in agreement with the theory. The null result is not no effect, it is an experiment too small.

### 6.2 Power and Sample Size

```r
za <- qnorm(0.975); zb <- qnorm(0.80); sigma <- 6; delta <- 0.10
mde <- (za + zb) * sigma * sqrt(2 / 15000)              # 0.194
n80 <- 2 * sigma^2 * (za + zb)^2 / delta^2              # ~56500
```

To detect the true 0.10 effect with 80% confidence, each arm needs about 56500 users, nearly four times the actual investment. Power climbs with sample size as follows:

| Per-arm sample size | Power to detect the 2% effect |
|---|---|
| 5000 | 0.13 |
| 10000 | 0.22 |
| 15000 | 0.30 |
| 25000 | 0.46 |
| 40000 | 0.65 |
| 56500 | 0.80 |
| 80000 | 0.92 |
| 120000 | 0.98 |

Power climbs only slowly, as $1/\sqrt{n}$, and lifting it from 30% to 80% by piling on sample costs nearly triple the traffic, which pushes the incentive toward the variance reduction of the next section.

![Power to detect the true 2% effect as a function of per-arm sample size. The orange line is the naive difference-in-means, the blue line is after CUPED adjustment, and the dashed line is 80% power. The opening 15000-user experiment (orange dot) has only 30% power; CUPED shifts the whole curve left, nearly halving the sample needed for the same power.](assets/fig/fig_17_power.svg)

### 6.3 CUPED: Cutting the Noise in Half

```r
theta <- cov(Y, X) / var(X)                    # 0.692
Ycv   <- Y - theta * (X - mean(X))
est   <- mean(Ycv[D == 1]) - mean(Ycv[D == 0])
```

Adjust with the pre-experiment revenue $X$ via CUPED. Here $\theta = 0.692$, and the variance of the adjusted estimator drops to 0.526 times the original (the empirical variance-retention ratio), agreeing with the theoretical $1 - \rho^2 = 0.513$, so about 47% of the variance is cut. Back to the opening 15000-user data: before adjustment $p = 0.392$, after adjustment the estimate is 0.107 with standard error 0.050, $t = 2.14$, $p = 0.032$; the same data, the same experiment, and CUPED makes the missed effect significant. The adjustment introduces no bias, and in simulation the bias of the CUPED estimate is 0.0009, negligible. The gain in sensitivity is real: at $n = 15000$ the power rises from 0.305 to 0.502, and the sample needed for 80% power falls from about 56500 to about 29000, nearly a halving, without recruiting a single extra user.

### 6.4 Peeking and the Always-Valid Test

```r
# A/A test (no effect), look daily, stop if ever p < 0.05
for (nk in daily_looks) if (abs(diff_k / se_k) > 1.96) { reject <- TRUE; break }
```

First see the harm of peeking. Monitoring a truly null A/A test continuously, the naive fixed-threshold false positive rate balloons with the number of looks: 1 look 0.049, 2 looks 0.082, 5 looks 0.145, 10 looks 0.193, all 14 looks 0.213. A redesign with no effect, allowed to be peeked at daily, has better than a one-in-five chance of being "found." Switching to mSPRT, monitoring the same 14 times continuously, the type I error is only 0.007, comfortably below 0.05. The cost is power: mSPRT's power against the true 2% effect is 0.312, below the fixed-sample test's 0.654 at the same terminal size, and this conservatism is precisely the price of "stop anytime," especially sharp for a small effect like this one that can only be seen once the full sample is collected.

![The type I error rate under continuous monitoring as a function of the number of looks. The orange line is the naive fixed 0.05 threshold, ballooning from the nominal 0.05 all the way to 0.21 at 14 looks; the blue line is mSPRT, staying comfortably below 0.05 no matter how many looks (0.007).](assets/fig/fig_17_peeking.svg)

### 6.5 Interference and Switchback

```r
# marketplace: capacity C shared; treated users cannibalize control's supply
slot_outcome <- function(want_prob) {
  wants  <- rbinom(length(want_prob), 1, want_prob)
  serve  <- if (sum(wants) > C) C / sum(wants) else 1   # capacity rationing
  wants * rbinom(length(wants), 1, serve)
}
```

Lumen's delivery marketplace has a shared-capacity constraint: 80 markets, 12 blocks per market, 600 users per block, and capacity 150. The new page raises users' willingness to order from 0.20 to 0.35. The true global effect (all-new versus all-old) is 0.040, with the all-control world producing 0.207 per person and the all-treatment world 0.248, a difference of 0.040. The estimates from three designs:

| Design | Estimate | Bias | Sampling SD |
|---|---|---|---|
| User-level A/B | 0.129 | +0.089 | 0.001 |
| Switchback | 0.040 | 0.000 | 0.002 |
| Cluster | 0.040 | 0.000 | 0.002 |

The user-level A/B gives 0.129, overstating the true global effect by a factor of 3.2, because treatment users cannibalized control's capacity, and both sides are biased. Switchback (flipping whole blocks) and cluster (assigning whole markets) both land almost exactly on 0.040, because each block or market experiences a globally consistent treatment state, and the cannibalization happens within groups, not contaminating across them. The cost is variance: the sampling SD of the two interference-robust designs (0.002) is markedly larger than the wildly biased but very stable user-level estimate (0.001), the interference-variance tradeoff laid out in numbers.

![The effect estimates from the three designs against the true global effect (green dashed line). The user-level A/B overstates badly because of capacity cannibalization, while switchback and cluster land on the truth but with wider confidence intervals, trading variance for unbiasedness.](assets/fig/fig_17_interf.svg)

The walkthrough of this section can be summarized as follows: on one and the same platform, randomization is perfect throughout, and the opening 15000-user experiment, being underpowered, misreports the true 2% effect as null; a power analysis could compute this before the experiment starts, CUPED uses pre-experiment data to halve the noise and make the effect significant on the same data, the sequential mSPRT holds the test level while allowing continuous monitoring where naive peeking pushes the false positive rate to 0.21, and the user-level A/B overstates the global effect by a factor of 3.2 because of capacity interference while switchback and cluster trade variance back for unbiasedness. The success or failure of an experiment is mostly in the design, not in the analysis.

## 7 Failure Modes and Robustness

Everything in the simulation is constructed, but real experiments come with a long list of pitfalls, and this section lays out the most common failure modes and their actionable remedies.

The most basic failure is being underpowered yet treating the null result as settled, already detailed in Section 1. Its twin error is post-hoc power analysis: taking a nonsignificant result and computing an "observed power" for self-comfort, which is circular reasoning, since observed power is just a monotone transform of the $p$ value and provides no information at all. The right discipline is to fix the MDE and sample size before the experiment starts, and even if the result is nonsignificant, report the confidence interval and let the width of the interval speak. A wide interval covering both 0 and a business-important effect should honestly be read as "not measured precisely," not as "no effect."

The second kind is sample ratio mismatch (SRM). The experiment is designed 50/50, but the two groups' sample proportions actually received deviate significantly (say 50.5/49.5, which can be anomalous in a large sample), and this is almost always the signal of a system bug: a flaw in the splitting logic, lost logs for one group, a redirect that drops some users. Once SRM appears, the premise of randomization has broken, and any effect estimate is untrustworthy. The response is to make the SRM test (a chi-square test on the group counts) the first automatic guardrail of every experiment, and if it fails, check the pipeline first and do not look at results.

The third kind is time-related contamination of the effect. The novelty effect is when users, briefly curious about a new feature, show inflated short-term metrics that fall back after a few weeks; the primacy effect is the reverse, users unaccustomed to the new page show deflated short-term metrics that slowly recover. Both make early readings unrepresentative of the long-run truth, and the response is to run long enough to see whether the metrics stabilize, or to analyze new versus returning users specifically and slice by exposure duration. Related are seasonality and day-of-week effects, where too short an experiment window mistakes the fluctuation of some special period for an effect, which is why an experiment is usually required to span at least a full weekly cycle.

The fourth kind is hiding interference behind a pretty $p$ value, already detailed in Section 3.5. Here we add one diagnostic in practice: if you suspect interference from a shared resource like capacity or inventory, you can stratify the experiment by the treatment share (say using different treatment proportions in different cities) and see whether the effect estimate varies systematically with the treatment share; if it does, the interference is confirmed. The cure is still to change the design (cluster or switchback), not to do any after-the-fact correction on user-level data.

The fifth kind is the confluence of multiple comparisons and selective reporting. Beyond the metric multiplicity of Section 3.4, there is a more insidious selectivity: the experiment fails to show an overall effect, so you slice into subgroups ("the overall effect is not significant, but it is for female users"), and slicing enough always finds a significant subgroup, which is subgroup fishing, essentially uncorrected multiple testing. The response is to pre-register the subgroups you will examine, apply a multiplicity correction to subgroup analyses, and treat any subgroup effect found after the fact as a hypothesis to be verified rather than a conclusion. Related is the winner's curse: picking the largest, most significant effect from a batch of experiments to ship, and its effect estimate is systematically too high, because being selected is precisely because noise favored it, so the prudent approach is to expect shrinkage on the post-ship effect, or to re-verify on held-out data.

Finally, the wisdom of Twyman's law: any result that looks especially interesting or especially anomalous is probably wrong. An absurdly large effect, a metric moving in a direction that defies common sense, the first response should be to check instrumentation (is the logging correct), check SRM, check the data pipeline, not to rush it into a ship report. The credibility of platform experiments is half in the design and half in this engineering discipline of staying skeptical of anomalous results.

Stringing these failure modes together, the difficulty of digital experimentation has never been in identification, since randomization delivers identification as a free lunch; the difficulty is in designing the experiment so that it sees what it should see, is not fooled by noise or peeking, and does not mistake interference for effect. Power analysis, CUPED, sequential methods, the SRM guardrail, and interference-robust designs are all safeguards provided to make an experiment with perfect identification really answer the right question, and none of them can substitute for careful thought about the design before the experiment starts.

## 8 Further Reading

::: {.readings}
Required reading, in suggested reading order:

- Kohavi, Tang and Xu (2020, Trustworthy Online Controlled Experiments). The authoritative survey textbook on platform experimentation, the practical panorama of every section of this chapter, best read first as an overview.
- Lewis and Reiley (2014, Quantitative Marketing and Economics). The classic case of underpower, the real-world version of Section 1 of this chapter, and read its estimates of noise and required sample size closely.
- Deng, Xu, Kohavi and Walker (2013, WSDM). The original CUPED paper, and read the derivation of the variance reduction and the choice of pre-experiment covariate closely.
- Johari, Koomen, Pekelis and Walsh (2017, KDD) and (2022, Operations Research). The peeking problem and always-valid inference, the source of mSPRT, and read the type I error control under continuous monitoring closely.
- Blake and Coey (2014, ACM EC). The cautionary work on interference in marketplace experiments, the empirical support for Section 3.5.

Further reading:

- Kohavi, Deng, Frasca, Walker, Xu and Pohlmann (2013, KDD). Lessons learned from large-scale online experimentation, the source of guardrails like SRM and Twyman's law.
- Bojinov, Simchi-Levi and Zhao (2023, Management Science). The optimal design and inference theory for switchback design, the methodological basis for Sections 3.6 and 4.4.
- Athey, Eckles and Imbens (2018, Journal of the American Statistical Association). Exact tests under network interference, building interference into an object of estimation rather than only a threat.
- Xu, Chen, Fernandez, Sinno and Bhasin (2015, KDD). LinkedIn's experimentation platform in practice, how engineering and statistics combine.
- Deng, Lu and Litz (2017, WSDM). The pitfalls, challenges, and solutions in analyzing online A/B tests, a practical checklist for avoiding traps.
- Bakshy, Eckles and Bernstein (2014, WWW). Facebook's experimentation infrastructure PlanOut, the design and execution of large-scale experiments.
:::

::: {.apa-refs}
- Athey, S., Eckles, D., & Imbens, G. W. (2018). Exact p-values for network interference. *Journal of the American Statistical Association, 113*(521), 230-240. https://doi.org/10.1080/01621459.2016.1241178
- Bakshy, E., Eckles, D., & Bernstein, M. S. (2014). Designing and deploying online field experiments. In *Proceedings of the 23rd International Conference on World Wide Web* (pp. 283-292). Association for Computing Machinery. https://doi.org/10.1145/2566486.2567967
- Blake, T., & Coey, D. (2014). Why marketplace experimentation is harder than it seems: The role of test-control interference. In *Proceedings of the Fifteenth ACM Conference on Economics and Computation* (pp. 567-582). Association for Computing Machinery. https://doi.org/10.1145/2600057.2602837
- Bojinov, I., Simchi-Levi, D., & Zhao, J. (2023). Design and analysis of switchback experiments. *Management Science, 69*(7), 3759-3777. https://doi.org/10.1287/mnsc.2022.4583
- Deng, A., Lu, J., & Litz, J. (2017). Trustworthy analysis of online A/B tests: Pitfalls, challenges and solutions. In *Proceedings of the Tenth ACM International Conference on Web Search and Data Mining* (pp. 641-649). Association for Computing Machinery. https://doi.org/10.1145/3018661.3018677
- Deng, A., Xu, Y., Kohavi, R., & Walker, T. (2013). Improving the sensitivity of online controlled experiments by utilizing pre-experiment data. In *Proceedings of the Sixth ACM International Conference on Web Search and Data Mining* (pp. 123-132). Association for Computing Machinery. https://doi.org/10.1145/2433396.2433413
- Johari, R., Koomen, P., Pekelis, L., & Walsh, D. (2017). Peeking at A/B tests: Why it matters, and what to do about it. In *Proceedings of the 23rd ACM SIGKDD International Conference on Knowledge Discovery and Data Mining* (pp. 1517-1525). Association for Computing Machinery. https://doi.org/10.1145/3097983.3097992
- Johari, R., Koomen, P., Pekelis, L., & Walsh, D. (2022). Always valid inference: Continuous monitoring of A/B tests. *Operations Research, 70*(3), 1806-1821. https://doi.org/10.1287/opre.2021.2135
- Kohavi, R., Deng, A., Frasca, B., Walker, T., Xu, Y., & Pohlmann, N. (2013). Online controlled experiments at large scale. In *Proceedings of the 19th ACM SIGKDD International Conference on Knowledge Discovery and Data Mining* (pp. 1168-1176). Association for Computing Machinery. https://doi.org/10.1145/2487575.2488217
- Kohavi, R., Tang, D., & Xu, Y. (2020). *Trustworthy online controlled experiments: A practical guide to A/B testing*. Cambridge University Press.
- Lewis, R. A., & Reiley, D. H. (2014). Online ads and offline sales: Measuring the effect of retail advertising via a controlled experiment on Yahoo! *Quantitative Marketing and Economics, 12*(3), 235-266. https://doi.org/10.1007/s11129-014-9146-6
- Xu, Y., Chen, N., Fernandez, A., Sinno, O., & Bhasin, A. (2015). From infrastructure to culture: A/B testing challenges in large scale social networks. In *Proceedings of the 21th ACM SIGKDD International Conference on Knowledge Discovery and Data Mining* (pp. 2227-2236). Association for Computing Machinery. https://doi.org/10.1145/2783258.2788602
:::
