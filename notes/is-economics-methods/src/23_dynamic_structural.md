---
title: "Dynamic Structural Models"
subtitle: "Optimal Stopping, Rust, and CCP Estimation"
seriesline: "Foundations of Information Systems Economics · Chapter 23"
series: "Foundations of Information Systems Economics: Empirical Methods for Information Systems Economics"
sidebartitle: "Foundations of Information Systems Economics"
sidebarsub: "Theory · Methods · Applications"
chapterlabel: "Chapter 23 · Dynamic Structural Models"
---

## Introduction

A merchant who keeps enduring ever slower legacy software is not necessarily doing so because upgrading has no value. He may be waiting for the next release, or he may believe next month's subsidy will be larger; not upgrading today is itself a choice about the future. If a researcher treats each period as an unrelated static decision, procrastination gets misread as low payoff, and the fixed upgrade cost gets mis-estimated right along with it. Dynamic problems most easily disguise themselves as static data, because the data record actions, not the actor's calculations about tomorrow.

Forward-looking decisions are organized by the Bellman equation. Today's action changes tomorrow's state, and tomorrow's value in turn feeds back onto today, so the optimal policy is an intertemporal fixed point. This is both the burden of the dynamic structural model and the source of its value: estimation has to grapple with the value function, but once you succeed you can answer counterfactuals like "how many merchants would a subsidy pull into upgrading early" that a static model cannot. The point here is not merely to fit the upgrade rate, but to separate the value of waiting, the state transition, and the flow payoff.

This chapter compares two classic computational routes. Rust's nested fixed point solves the value function under each candidate parameter, and an outer likelihood searches for the parameters; Hotz-Miller's conditional choice probability method instead exploits the value information contained in observed choice probabilities, sidestepping the repeated solution of the inner fixed point. A case with one thousand merchants on the Nimbus platform making upgrade choices over thirty periods will put the static and dynamic models through the same test. Magnac-Thesmar's identification result then reminds us: the discount factor generally does not reveal itself just because the panel is long enough, it has to be fixed, calibrated, or pinned down by an extra exclusion restriction.

The computational routes can differ, but the review standard is the same: does the state contain the variables that genuinely affect the future, is the conditional independence assumption credible, where does the discount factor come from, and how do persistent heterogeneity and sparse states change the counterfactual. The choice between NFXP and CCP is ultimately not an algorithm race, but a tradeoff among computational cost, statistical error, and model transparency. The purpose of structural estimation is not to replicate historical choices, but to use a behavioral model that survives identification scrutiny in order to change, carefully, a future that has not yet happened.

## 1 A Number That Doesn't Add Up

::: {.case}
Nimbus is a fictional e-commerce platform. Merchants on the platform use a suite of store tools (software for ordering, inventory, and marketing), and the tools gradually go stale as the platform iterates. The degree of obsolescence is recorded as staleness $x$: the larger $x$ is, the higher the cost of stale tools dragging down sales and adding friction. Each period a merchant decides whether to keep using the old tools (bearing a friction cost that rises with $x$) or to upgrade to the latest version (paying a fixed upgrade cost and resetting $x$ to 0). We observe one thousand merchants over thirty periods, along with their staleness and upgrade decisions. The platform wants to know two structural parameters: the friction cost $\theta_1$ per unit of obsolescence, and the fixed cost $RC$ of a single upgrade, so it can answer policy questions like "how much would an upgrade subsidy boost upgrading."
:::

Let us first try the lazy approach: pretend merchants are myopic, looking only at the current period's gains and losses, planning nothing for the future. This amounts to setting the discount factor to zero, then using a static choice model to fit the upgrade decisions and estimating the friction cost $\theta_1$ and the upgrade cost $RC$. The data fit fine: the myopic model can roughly match the upgrade frequency at each level of staleness. But the estimated structural parameters are wildly wrong: the friction cost $\theta_1$ comes out at $0.42$, while the true value is only $0.15$, nearly triple the truth.

Where does it go wrong? It goes wrong because merchants are in fact forward-looking. They do not wait until the tools are unusable before upgrading; they act while obsolescence is still moderate, because they foresee that dragging on will let friction accumulate and future costs rise, so upgrading early is worthwhile. The myopic model sees merchants upgrading at not-very-high staleness and cannot make sense of this foresight, so it can only attribute it to "the immediate cost of obsolescence must be very high," and thus estimates $\theta_1$ too large. The truth is that the immediate friction cost is not high; it is the expectation of the future that drives the early upgrade. Put the discount factor back at $0.95$ and build the forward-looking behavior into the model, and Rust's dynamic estimation puts $\theta_1$ back at $0.154$ and $RC$ back at $4.08$, both hugging the truth.

Mis-estimating the structural parameters is not even the worst of it; the worst is that it will get the policy counterfactual wrong. Suppose the platform rolls out a subsidy that cuts the upgrade cost by $30\%$, how many merchants would upgrade early as a result? The true answer is that the per-period upgrade rate rises by $0.053$ (from $0.162$ to $0.215$). The dynamic model, using its estimated parameters to compute this counterfactual, also gets $0.053$, matching the truth. The myopic model gets $0.076$ instead, overstating the subsidy's effect by forty percent. A model that gets even the discount factor wrong will get both the parameter estimates and the counterfactual wrong, and the counterfactual is exactly what the platform needs most when making this subsidy decision.

What this chapter is about is precisely how to correctly estimate such a forward-looking model. Section 2 will write out the Bellman equation and the fixed point of the value function clearly, Section 3 solves it using a contraction mapping and also lays out a counterintuitive identification fact: the discount factor is generally not identified and has to lean on extra restrictions. Section 4 gives the two main routes to estimation, Rust's nested fixed point and Hotz-Miller's CCP method, the latter cleverly sidestepping the inner fixed point. Section 6 will run all of this on the Nimbus data, showing how the myopic errors get corrected one by one, and how the counterfactual gets computed credibly.

## 2 The Economic Model and the Estimand

Before rolling up our sleeves, let us write down the dynamic decision model and the target quantity clearly. At the core are the Bellman equation and the value function.

The target quantities this chapter delivers come in two kinds. The first is the structural parameters, the deep parameters that characterize merchant preferences and costs, such as the friction cost $\theta_1$ and the upgrade cost $RC$; they do not change with the environment and are the foundation for doing counterfactuals. The second is counterfactuals, that is, how behavior would change after altering some environmental variable (such as subsidizing upgrades or changing the obsolescence rate). The second kind is the ultimate goal of the structural method: estimating structural parameters is not for the sake of those few numbers themselves, but in order to predict, under the assumption that the parameters stay fixed, the world after a policy change. This is especially true for dynamic models, because only by building in forward-looking behavior can they credibly predict how a forward-looking decision-maker would respond to a new policy.

Let us first write the dynamic decision. Time is discrete, indexed by $t$. The state is the degree of obsolescence $x_t \in \{0, 1, \dots, K-1\}$. Each period a merchant chooses between two actions: keep using the old tools ($i = 0$) or upgrade ($i = 1$). A local notation for this chapter: we use $i \in \{0, 1\}$ for the action and $n$ for the merchant, replacing the book-wide convention of using $i$ for the individual. The flow utility (negative cost) is

$$u(x, i) + \varepsilon_i, \qquad u(x, 0) = -\theta_1\, x, \quad u(x, 1) = -RC$$

The cost of keeping the old tools rises linearly with obsolescence, and upgrading pays a fixed cost $RC$. $\varepsilon_i$ is an action-specific shock that the merchant sees but the researcher does not. The state transition is: keeping the old tools drifts $x$ upward with some distribution (tools get older with use), while upgrading resets $x$ to 0 (installing the latest version). The merchant discounts the future with discount factor $\beta$, maximizing the sum of expected discounted utility.

At the heart of forward-looking behavior is the value function. Define the value function $V(x)$ as the expected discounted utility obtainable from state $x$ under optimal decisions thereafter. It satisfies the Bellman equation

$$V(x) = \mathbb{E}_\varepsilon \max_{i \in \{0, 1\}}\Big\{ u(x, i) + \varepsilon_i + \beta\, \mathbb{E}\big[V(x') \mid x, i\big] \Big\}$$

The term in brackets, $u(x,i) + \beta\,\mathbb{E}[V(x') \mid x, i]$, is called the action-specific value $v(x, i)$; it is the flow utility of choosing action $i$ plus the discounted future value. The merchant chooses the action with the largest $v(x,i) + \varepsilon_i$. $V$ appears on both sides of its own equation (the future $V(x')$ depends in turn on the still more distant future), so it is a fixed point that has to be solved for before choice probabilities can be computed. This is exactly what makes the dynamic model harder than the static one, and it is the core around which the whole chapter's estimation techniques revolve.

::: {.intuition}
Feel the foresight through a two-period example. Imagine there are only two periods left. The last period has no future, so the merchant looks only at the current cost, which is a static problem. But the second-to-last period is different: whether he upgrades or keeps using the old tools has to account for how this choice affects his situation in the last period. If he keeps using this period and lets $x$ rise, he carries a higher friction cost into the last period; if he upgrades this period and zeroes out $x$, the last period is easy. So even if the immediate cost of upgrading is higher this period, as long as the future cost it saves is large enough, the forward-looking merchant will upgrade now. What the Bellman equation does is generalize this reasoning of "how this period's choice affects the future" from two periods to infinitely many, characterizing it all at once with a self-consistent value function. The myopic model sets $\beta$ to zero, which amounts to cutting out this entire chain of forward-looking reasoning, leaving only the last-period, purely static tradeoff, and that is the root of why it mis-estimates the cost.
:::

The summary is as follows: the model of dynamic decision-making is the optimal stopping problem characterized by the Bellman equation; the estimand is the structural parameters (friction cost, upgrade cost) and the counterfactuals; the value function is a fixed point that must be solved for, and the forward-looking behavior lives entirely inside it; myopia amounts to cutting out discounting and discarding foresight, which is the source of the error in Section 1 and the starting point for identification and estimation in the next section.

## 3 Identification

Identification logically precedes estimation. This section answers: under what assumptions can the dynamic discrete choice model be estimated, which parameters can be identified, and which are in principle not identifiable. This is the most detailed analysis in the chapter, split into numbered subsections, each of which gives the intuition first and then the formalization.

### 3.1 Conditional Independence and Type-I EV: Why There Is a Closed Form

The value function is a fixed point and generally has no closed form, but Rust uses two assumptions to give choice probabilities a closed form, making estimation feasible. These two assumptions are the technical bedrock of the whole chapter.

The first is conditional independence (CI). It says that, given the current state $x$ and the chosen action $i$, the transition to the next-period state $x'$ is independent of the current shock $\varepsilon$, and $\varepsilon$ is independent and identically distributed over time. The intuition is that the unobserved shock affects only the current choice, and does not quietly influence the state evolution through some other channel. This decouples a high-dimensional dynamic problem: how the state transitions depends only on the observable $(x, i)$, and the shock's influence is isolated to the current period.

The second is that $\varepsilon$ follows a Type-I extreme value distribution, just as in the logit of Chapter 19. With these two assumptions, once the action-specific value $v(x, i)$ is known, the choice probability is the familiar logit form

$$P(i \mid x) = \frac{\exp(v(x, i))}{\sum_j \exp(v(x, j))}$$

and the expected maximum has a closed-form log-sum expression, $\mathbb{E}_\varepsilon \max_j\{v(x,j) + \varepsilon_j\} = \log\sum_j \exp(v(x,j))$ (up to an Euler constant). This log-sum is exactly the incarnation in the dynamic setting of the inclusive value from Chapter 19; it computes that expected-maximum-over-$\varepsilon$ in the Bellman equation in closed form, letting the fixed point of the value function be solved numerically. CI plus Type-I EV is the key that makes dynamic discrete choice estimable.

### 3.2 Rust's Contraction Mapping: Solving the Value Function

With the closed-form log-sum, the fixed point of the value function can be solved numerically. Writing the Bellman equation out with the log-sum, define the operator $T_\theta$ that maps a value function to

$$(T_\theta V)(x) = \log\Big[ \exp\big(u(x,0) + \beta\, \mathbb{E}[V(x') \mid x, 0]\big) + \exp\big(u(x,1) + \beta\, \mathbb{E}[V(x') \mid x, 1]\big) \Big]$$

Rust's Theorem 1 proves that $T_\theta$ is a contraction mapping, with contraction modulus exactly the discount factor $\beta < 1$. This means that starting from any initial value function, repeatedly applying $T_\theta$ converges to the unique fixed point $V_\theta$. In practice this is just iterating the operator to convergence, or accelerating with a Newton-type method. The contraction mapping guarantees the value function exists and is unique, and it is the algorithmic cornerstone of the entire dynamic estimation. Solving this fixed point on Nimbus, the optimal policy we obtain is a curve where the upgrade probability rises with the degree of obsolescence: at $x = 0$ upgrading is almost never chosen (probability $0.018$), at $x = 10$ it exceeds half ($0.564$), and at $x = 19$ it reaches as high as $0.875$. This upgrade probability rising with the state is what an optimal stopping policy looks like.

### 3.3 Magnac-Thesmar: The Discount Factor Is Not Identified

Dynamic models have a counterintuitive but extremely important identification fact: the discount factor $\beta$ and the flow utility generally cannot be separately identified from data. This was clarified by Magnac and Thesmar (2002), and not understanding it will lead you to misread the estimation results.

The root of the problem lies in what the data can provide. The data can nonparametrically identify two things: the conditional choice probability $P(i \mid x)$ and the state transition $P(x' \mid x, i)$. But the dynamic model's flow utility has too many degrees of freedom for these two to pin it down. Do a count: with discrete states, the flow utility has $|X| \cdot |J|$ free values, while the data contain only $|X| \cdot (|J| - 1)$ linearly independent choice probabilities, so it is underidentified. More concretely, given any discount factor $\beta$, any utility for a reference action, and any shock distribution, you can back out a set of flow utilities that make the model perfectly reproduce the same set of choice probabilities. In other words, different combinations of "discount factor plus flow utility" are observationally equivalent, and the data cannot tell them apart.

To restore identification, restrictions must be added. The standard set (Magnac-Thesmar's sufficient conditions) is: the transition is identified, utility is additively separable, the shock distribution is known (such as Type-I EV), the utility of a reference action is normalized to a known value (such as setting some baseline for "keep using" to zero), and the discount factor $\beta$ is fixed at a known value. With these five in place, the flow utilities of the remaining actions are point identified. Two things deserve special emphasis. First, setting the reference action's utility to zero is not a harmless normalization in a dynamic model but a substantive restriction, because the utility of one state affects the incentives at all states through the value function, and getting it wrong will bias the structural parameters. Second, $\beta$ can usually only be fixed, not estimated. To truly identify $\beta$, an exclusion restriction is needed: a variable that affects only the future value and not the flow utility, which acts like an instrumental variable in the dynamic setting, using the variation of "moving the future without moving the present" to separate the discount factor from the flow utility. Like the vast majority of applications, this chapter fixes $\beta$ at $0.95$ and estimates the rest; this is not laziness but a limit imposed by identification.

The main points of this section can be summarized as follows: CI plus Type-I EV gives closed-form logit choice probabilities and a log-sum, letting the fixed point of the value function be solved numerically; Rust's contraction mapping guarantees the fixed point exists and is unique and provides the algorithm; but the discount factor and flow utility generally cannot be separately identified, so $\beta$ must be fixed and the reference utility normalized (a substantive restriction), and $\beta$ itself can only be estimated with an exclusion restriction. None of this yet says how to compute things, which is the next section.

## 4 Estimation

This section brings the identification of Section 3 down to operational estimation: Rust's nested fixed point, Su-Judd's MPEC, Hotz-Miller's CCP, finite dependence, and the counterfactual.

### 4.1 NFXP: Nested Fixed Point

Rust's nested fixed point (NFXP) is the most direct estimation method, with two nested layers. The outer layer does maximum likelihood over the structural parameters $(\theta_1, RC)$. The inner layer, for each set of parameters the outer layer proposes, uses the contraction mapping of Section 3.2 to solve the value-function fixed point $V_\theta$, computes the choice probability at each state accordingly, and then computes the log-likelihood

$$\ell(\theta) = \sum_{n, t} \Big[ i_{nt} \log P(1 \mid x_{nt}; \theta) + (1 - i_{nt}) \log P(0 \mid x_{nt}; \theta) \Big]$$

Each time the optimizer proposes a set of parameters, the inner layer re-solves the fixed point once and recomputes the likelihood once, until the outer layer converges. On Nimbus, fixing $\beta = 0.95$, NFXP estimates the friction cost as $\theta_1 = 0.154$ (standard error $0.004$) and the upgrade cost as $RC = 4.082$ (standard error $0.058$), both hugging the truths of $0.15$ and $4.0$. NFXP's advantage is that it is statistically efficient and direct; its drawback is that the inner fixed point must be re-solved at every outer iteration, which is computationally expensive when the state space is large, and if the inner layer is not solved accurately enough it will bias the outer layer, so the tolerance has to be set very tight.

### 4.2 MPEC: Turning the Fixed Point into a Constraint

NFXP's inner fixed point is the source of the computational burden. Su and Judd (2012) give another route, MPEC (mathematical programming with equilibrium constraints). Instead of nesting the fixed point inside the objective function, it treats the Bellman equation as a set of constraints, letting the optimizer search over the structural parameters and the value function simultaneously, with the constraints satisfied only at the optimum. The key point is completely parallel to Chapter 20's discussion of NFXP versus MPEC for BLP: the two are the same estimator, just with different algorithms, converging to the same solution, with the difference lying entirely in computational convenience. MPEC has many parameters (the whole value function enters the optimization) but is insensitive to the inner tolerance and has easy-to-compute derivatives; NFXP has few parameters but is sensitive to the inner precision. This again confirms that general mantra: identification and algorithm are two separate things, and the same estimator can have several ways of being computed.

### 4.3 Hotz-Miller: Using Choice Probabilities to Sidestep the Fixed Point

NFXP re-solves the fixed point at every outer iteration, and Hotz and Miller (1993) discovered that this step can often be skipped. Their insight is that the observed conditional choice probabilities themselves encode value information, which can be used to back out the value difference without solving the fixed point.

The key is an inversion formula. Under Type-I EV, the value difference between two actions corresponds one-to-one with the ratio of their choice probabilities:

$$v(x, 1) - v(x, 0) = \log P(1 \mid x) - \log P(0 \mid x)$$

This is the same class of idea as Chapter 20's Berry inversion backing out utility from shares, except that here we back out the dynamic value difference from choice probabilities. Going further, Hotz-Miller-Sanders-Smith prove that the value function can be expressed entirely in terms of choice probabilities: the shock expectation for the chosen action has the closed form $\mathbb{E}[\varepsilon \mid i, x] = \gamma - \log P(i \mid x)$, so the value function satisfies a linear equation with choice probabilities as known quantities, $V = (I - \beta F_0)^{-1} \psi_0(P)$, which can be solved with a single matrix inversion, with no need to iterate the fixed point.

CCP estimation is therefore two steps. First, directly estimate the conditional choice probabilities $\hat P(i \mid x)$ (frequencies) and the state transition from the data. Second, use the inversion above to express the value difference as a linear function of the structural parameters, then make it match the ratio of choice probabilities in the data to estimate the structural parameters. The whole process never solves a fixed point once. On Nimbus, CCP estimation puts $\theta_1$ at $0.149$ and $RC$ at $4.015$, both agreeing with NFXP's $0.154$ and $4.082$ and with the truths of $0.15$ and $4.0$, while requiring far less computation. The price is that the first step's choice probabilities have to be estimated accurately, and states visited rarely have noisy frequencies, requiring smoothing or more data; also, CCP's statistical efficiency is generally slightly lower than NFXP's.

### 4.4 Finite Dependence: Letting the Continuation Values Cancel

In the CCP method, the expression for the value function still involves the entire future path. Arcidiacono and Miller (2011) further simplify with finite dependence: if two different current choices lead to the same state distribution a few periods later, then their subsequent continuation values are the same and cancel when taking the value difference, so estimation need only look a finite number of periods ahead. The upgrade problem in this chapter naturally has the strongest form of finite dependence, namely the renewal property: no matter how high the current staleness, upgrading resets it to 0 and continuation begins from the same starting point. This means that "upgrade now" and "continue now, upgrade next period" quickly merge into the same state afterward, the continuation values cancel, and the value difference is determined only by the flow utilities and choice probabilities of the next few periods. The renewal property makes CCP estimation of optimal-stopping problems like upgrading especially clean, which is also why Rust's bus engines and this chapter's tool upgrade both revolve around renewal.

### 4.5 Counterfactuals: The Purpose of Structural Estimation

Estimating the structural parameters is only a means; the counterfactual is the end. With structural parameters that do not change with policy, you can change the environment, re-solve the model, and predict new behavior. The procedure is to plug in the changed environment (such as a reduced upgrade cost), re-solve the value-function fixed point, and compute the new choice probabilities and the steady-state behavior they produce. On Nimbus, simulating the platform subsidizing the upgrade cost by $30\%$: using the true parameters, the per-period upgrade rate rises from $0.162$ to $0.215$, an increase of $0.053$; using the dynamically estimated parameters, it rises by $0.053$, matching the truth; using the myopic model's parameters, it rises by $0.076$, overstating by forty percent. The credibility of the counterfactual rests entirely on whether the structural parameters are estimated correctly and on the assumption that the parameters really do not change after the policy change. The reason the dynamic model is worth that large computational cost is precisely that it builds in forward-looking behavior, so it can credibly predict how a forward-looking decision-maker would respond to a new policy, while the myopic model errs systematically at this step.

The route of this section is now complete: NFXP does outer maximum likelihood and inner fixed-point solving; MPEC turns the fixed point into a constraint, the same estimator with another algorithm; Hotz-Miller uses the choice-probability inversion to sidestep the fixed point, a two-step estimator; finite dependence further simplifies by letting continuation values cancel; and the counterfactual is the purpose of all this, its credibility tied to the structural parameters and the parameter-invariance assumption.

## 5 Anchor Papers

Methods only stand firm when they land in real research. Three anchor papers: one erects the whole method of dynamic discrete choice, one gives the CCP revolution that sidesteps the fixed point, and one pushes the dynamic to durable goods and technology adoption, highly relevant to platforms. Each is laid out along five elements, paper, method, data, results, and limitations, with attention to how the identification assumptions are defended.

### 5.1 Rust (1987)

::: {.case}
Paper: "Optimal Replacement of GMC Bus Engines: An Empirical Model of Harold Zurcher," Econometrica. It is the foundational work of dynamic discrete choice structural estimation, turning the optimal stopping problem into an estimable empirical model, and Harold Zurcher's bus engine replacement is the prototype for all who followed.

Method: model the bus maintenance superintendent's engine replacement decision as an optimal stopping problem, where the state is engine mileage and each period one chooses to keep using (maintenance cost rising with mileage) or to replace (paying a fixed replacement cost, mileage reset to zero). Under conditional independence and Type-I EV assumptions, the value function satisfies a contraction mapping, estimated with nested fixed point plus maximum likelihood. This chapter's tool upgrade problem is its direct analog.

Data: years of engine mileage and replacement records for the Madison, Wisconsin bus fleet, with Harold Zurcher as the fleet's maintenance superintendent.

Results: estimated the maintenance cost of mileage and the engine replacement cost, the model reproduces the observed replacement patterns, and it can be used for counterfactuals (such as how a change in the replacement cost affects replacement frequency). This paper demonstrated how to turn a dynamic decision that seems merely describable into a structural model that can estimate structural parameters and do counterfactuals.

Limitations: the discount factor is not identified in this framework, and Rust fixes it; the conditional independence assumption rules out unobserved heterogeneity operating through state evolution; the state space must be discretized, and computation explodes as the dimension grows. These limitations gave rise to Hotz-Miller's later CCP, Magnac-Thesmar's identification analysis, and various methods for handling unobserved heterogeneity.
:::

### 5.2 Hotz and Miller (1993)

::: {.case}
Paper: "Conditional Choice Probabilities and the Estimation of Dynamic Models," Review of Economic Studies. It provides the way to sidestep Rust's inner fixed point, a computational revolution in dynamic discrete choice estimation that made dynamic models with large state spaces feasible.

Method: prove that the value difference corresponds one-to-one with the conditional choice probabilities, so value can be backed out from observed choice probabilities without solving the value function's fixed point. On this basis propose a two-step estimator: the first step estimates choice probabilities from data, and the second step uses the inversion to estimate the structural parameters. At its core is the inversion idea that "reads" the dynamic value out.

Data: the paper uses the National Fertility Survey to estimate parents' dynamic contraception choices and a fertility model, demonstrating the CCP method with it.

Results: the CCP method gives estimates consistent with Rust's NFXP while the computation drops sharply, with the advantage especially clear when the state space is large and the fixed point is hard to solve. It pushed dynamic structural estimation from "can only handle small models" toward "can handle realistically sized models."

Limitations: the first step's choice probabilities have to be estimated accurately, and rarely visited states are noisy; the two-step estimator is generally slightly less statistically efficient than the one-step NFXP; the inversion depends on the same CI and distributional assumptions as Rust. The later nested pseudo-likelihood of Aguirregabiria and Mira (2002) connects CCP and NFXP into a spectrum, striking a middle ground between efficiency and computation.
:::

### 5.3 Gowrisankaran and Rysman (2012)

::: {.case}
Paper: "Dynamics of Consumer Demand for New Durable Goods," Journal of Political Economy. It applies dynamic discrete choice to the timing of durable goods purchases: consumers decide when to buy and which generation to buy, highly isomorphic to when merchants on a platform adopt new tools, a modern template for dynamic demand.

Method: dynamic random-coefficient logit demand, where consumers are forward-looking, weighing buying now against waiting for a better and cheaper next generation. The key computational trick is inclusive value sufficiency, which uses a scalar inclusive value to summarize the evolution of the future product environment, compressing a high-dimensional state space to a solvable size, and then uses Rust-like dynamic programming plus BLP-style demand estimation. It stitches together the BLP demand of Chapter 20 with the dynamic choice of this chapter.

Data: market-level sales, prices, and product characteristics panels for durable consumer electronics such as digital camcorders, whose product quality iterates rapidly and prices keep falling, an ideal battleground for dynamic purchase timing.

Results: consumers do have significant foresight and will delay purchase to wait for a better next generation; ignoring this dynamic mis-estimates both the price elasticity and the value of new products. The paper demonstrated how to estimate forward-looking behavior in a realistically sized dynamic demand problem.

Limitations: inclusive value sufficiency is an approximation, requiring the future environment to be adequately summarized by a scalar; the discount factor is likewise fixed; identification of dynamic demand places high demands on exogenous variation in price and quality. This paper is especially instructive for platform research: merchants adopting new tools and consumers adopting new technology are both this kind of "when to get on board" dynamic decision.
:::

Taken together, the three make the meaning of anchoring clear: Rust erects the whole method of dynamic discrete choice and fixed-point estimation, Hotz-Miller uses CCP inversion to sidestep the fixed point and make large models feasible, and Gowrisankaran-Rysman push the dynamic to durable goods and technology adoption, a domain highly relevant to platforms. Progress in dynamic structural methods has always revolved around one main thread: how to estimate forward-looking behavior at an acceptable computational cost, and to do credible counterfactuals on that basis.

## 6 A Full Walkthrough on the Nimbus Data

Now let us run the tools of Section 4 all the way through on Nimbus's merchant upgrade data. The code below uses R 4.5.3, fixing set.seed(19) to ensure reproducibility, and every number cited in the text comes from the actual run output of this code.

### 6.1 DGP and Model Solution

The design parameters are as follows: the degree of obsolescence $x \in \{0, \dots, 19\}$, the discount factor $\beta = 0.95$, the friction cost $\theta_1 = 0.15$, the upgrade cost $RC = 4.0$; when keeping the old tools, $x$ drifts up by $0$, $1$, or $2$ notches with probabilities $(0.30, 0.50, 0.20)$, and upgrading resets it to 0. The value function is solved by contraction-mapping iteration.

```r
solve_V <- function(theta1, RC, beta){          # Rust's EV contraction mapping
  u0 <- -theta1*xs; u1 <- rep(-RC, K); V <- rep(0, K)
  repeat {
    v0 <- u0 + beta*(Tkeep %*% V); v1 <- u1 + beta*(Tup %*% V)
    Vn <- log(exp(v0) + exp(v1))                 # log-sum (stabilized numerically)
    if (max(abs(Vn - V)) < 1e-12) break; V <- Vn
  }
  list(V = Vn, Pup = 1/(1 + exp(v0 - v1)))       # P(upgrade | x)
}
```

The optimal policy solved out is a curve where the upgrade probability rises with the degree of obsolescence: $0.018$ at $x = 0$, $0.231$ at $x = 5$, $0.564$ at $x = 10$, $0.781$ at $x = 15$, and $0.875$ at $x = 19$. Simulating one thousand merchants over thirty periods each, thirty thousand merchant-periods in all, gives an upgrade rate of $0.146$ and an average obsolescence of $3.22$.

![The optimal upgrade policy solved out: upgrade probability rises with the degree of obsolescence (the navy curve is the solved model policy), and the red dots are the actual upgrade frequencies at each obsolescence level in the data, the two agreeing. This upgrade probability rising with the state is what an optimal stopping (renewal) policy looks like.](assets/fig/fig_23_policy.svg)

### 6.2 The Myopic Model: Reproducing the Opening Error

```r
nll <- function(par, beta_use){                  # inner solves fixed point, outer maximizes likelihood
  s <- solve_V(par[1], par[2], beta_use)
  -sum(ifelse(dat$i==1, log(s$Pup[dat$x+1]), log(1-s$Pup[dat$x+1])))
}
fit_myopic <- optim(c(0.1,3), nll, beta_use = 0)  # beta = 0: cut out foresight
```

Setting the discount factor to zero, the myopic model estimates the friction cost as $\theta_1 = 0.418$ and the upgrade cost as $RC = 3.436$, with the friction cost nearly triple the truth of $0.15$. The reason was explained in Section 1: merchants are in fact forward-looking and upgrade at moderate obsolescence, and the myopic model cannot make sense of this foresight, so it can only attribute it to an excessively high immediate friction cost.

### 6.3 Rust NFXP: Estimating the Parameters Back

```r
fit_nfxp <- optim(c(0.1,3), nll, beta_use = 0.95)  # dynamic: beta = 0.95
```

Fixing $\beta = 0.95$, NFXP estimates $\theta_1$ as $0.154$ (standard error $0.004$) and $RC$ as $4.082$ (standard error $0.058$), both hugging the truth. Once the forward-looking behavior is built into the model, the structural parameters are estimated correctly.

### 6.4 Hotz-Miller CCP: Sidestepping the Fixed Point

```r
Phat <- tapply(dat$i, dat$x, mean)               # step 1: upgrade frequency at each state (CCP)
V    <- solve(diag(K) - beta*M) %*% psi0          # CCP expresses value, one matrix inversion
# v1 - v0 = log(P1/P0) is linear in the structural parameters, OLS estimates theta1, RC
```

Using the conditional choice probability inversion, never solving a fixed point once, CCP estimates $\theta_1$ as $0.149$ and $RC$ as $4.015$, agreeing with NFXP's $0.154$ and $4.082$ and with the truths of $0.15$ and $4.0$, while requiring far less computation. This confirms Hotz-Miller's core insight: the observed choice probabilities themselves measure value, and using them one can sidestep the expensive inner fixed point.

### 6.5 Counterfactual: The Upgrade Subsidy

```r
# cut the upgrade cost by 30%, re-solve the model, compare upgrade rates
cf <- function(th1, RC, bet)
  c(base = agg(solve_V(th1, RC, bet)), sub = agg(solve_V(th1, RC*0.7, bet)))
```

Simulate the platform subsidizing the upgrade cost by $30\%$. Using the true parameters, the per-period upgrade rate rises from $0.162$ to $0.215$, an increase of $0.053$ (the baseline rate $0.162$ is the model-implied steady-state upgrade rate, on a different footing from the finite-sample upgrade rate $0.146$ obtained in Section 6.1 by starting from $x=0$ and simulating only 30 periods). Using the dynamically estimated parameters, it rises by $0.053$, matching the truth. Using the myopic model's parameters, it rises by $0.076$, overstating the subsidy's effect by forty percent. This is the payoff of structural estimation: the correct dynamic model not only estimates the parameters right but also gives a credible prediction for the policy counterfactual, while the myopic model errs at both ends, parameters and counterfactual.

![The myopic model's twofold error. Left: it estimates the friction cost of obsolescence $\theta_1$ as $0.42$, while the truth is only $0.15$, and the dynamic model estimates it back at $0.15$. Right: it predicts a $30\%$ upgrade subsidy raises the upgrade rate by $0.076$, while the truth is $0.053$, and the dynamic model computes $0.053$. Cut out foresight, and both parameter estimation and the counterfactual are systematically biased.](assets/fig/fig_23_myopic.svg)

### 6.6 Full Reconciliation

The walkthrough of this section can be summarized as follows: the optimal policy solved out has the upgrade probability rising with obsolescence (from $0.018$ at $x=0$ to $0.875$ at $x=19$); the myopic model ($\beta=0$) estimates the friction cost as $0.418$, far above the truth of $0.15$; Rust NFXP ($\beta=0.95$) estimates it back at $0.154$ with $RC$ back at $4.082$, and Hotz-Miller CCP, without solving a fixed point, also gives $0.149$ and $4.015$, both close to the truth; and for the counterfactual of a $30\%$ upgrade subsidy, the dynamic model computes the true $0.053$, while the myopic model errs at $0.076$. Only when the forward-looking behavior is estimated correctly are the parameters and the counterfactual correct, and that is the meaning of the entire computational cost of a dynamic structural model.

## 7 Failure Modes and Robustness

In the simulation the model structure is constructed, but in real research the various assumptions can fail at any moment. This section lays out the most common ways of failing and their operational remedies.

The non-identification of the discount factor is the most easily misread point of dynamic models. As Section 3.3 explained, $\beta$ and the flow utility generally cannot be separated, and in practice $\beta$ is fixed (such as $0.95$). The danger is treating this fixing as harmless, when in fact the conclusions may be sensitive to $\beta$, especially the counterfactuals. The honest approach is to report how the conclusions change when $\beta$ takes several values, or to formally identify $\beta$ when there is an exclusion restriction (a variable that moves only the future value and not the flow utility). Fixing $\beta$ without a sensitivity analysis amounts to treating an unidentified parameter as known, a common form of overconfidence.

The conditional independence assumption rules out an important class of unobserved heterogeneity. CI requires that the unobserved shock not operate through the state evolution, but in reality merchants may differ persistently along dimensions the researcher cannot see (some are born early adopters, some conservative), and this persistent heterogeneity violates CI and biases the structural parameters. The remedy is to explicitly model the unobserved types (such as finite mixtures, the EM algorithm); Arcidiacono-Miller extend the CCP method to the case with unobserved heterogeneity. Ignoring persistent heterogeneity and forcing CI is a hidden and common error in dynamic models.

The validity of the counterfactual is the ultimate threat to dynamic structural methods, essentially the Lucas critique. The counterfactual assumes the structural parameters do not change after the policy change, but if the policy itself changes how merchants form expectations, or changes an unmodeled part of the environment, the parameters need not be stable and the counterfactual prediction fails accordingly. More subtly, even if the flow utility is only partially determined for want of identification, some counterfactuals may still be identified and others not; Kalouptsidi and others characterize which counterfactuals are immune to the unidentified utility. The remedy is to restrict the counterfactual to the range where there is reason for the parameters to stay fixed, and to make clear which assumptions, untested by the data, the counterfactual depends on.

There are a few more concrete gates. The state space must be discretized, and too few notches lose the dynamic while too many make computation explode; the curse of dimensionality is the perennial headache of dynamic models, and dimension-reducing tricks like inclusive value sufficiency exist precisely to deal with it. If the transition probabilities are estimated inaccurately, both the value function and the choice probabilities are biased, and the CCP method especially depends on the first-step transition and choice probability estimates. In finite samples, the choice probabilities at rarely visited states are noisy, and CCP's first step needs smoothing or more data. Finally, this chapter has covered only single-agent dynamics, while real platform competition is a multi-agent dynamic game where each firm's optimal policy depends on the rivals' policies, the equilibrium is a Markov perfect equilibrium, and estimation uses the two-step CCP methods of Bajari-Benkard-Levin or Aguirregabiria-Mira, extending this chapter's Hotz-Miller idea to games, at the frontier of dynamic structural methods.

Stringing these failure modes together, the credibility of a dynamic structural model comes down to two things: whether the discount factor and the identification restrictions are confronted head-on, and whether the parameter-invariance assumption on which the counterfactual depends holds up. The sensitivity to $\beta$, unobserved heterogeneity, the validity of the counterfactual, and the estimation of the state space and transitions are all evidence offered around these two, and none of them can be replaced by "the fixed point was solved and the parameters were estimated."

## 8 Further Reading

::: {.readings}
Required reading, in suggested order:

- Rust (1987, Econometrica). The foundational work of dynamic discrete choice, read it first to understand the full logic of optimal stopping, the CI assumption, the contraction mapping, and NFXP.
- Hotz and Miller (1993, Review of Economic Studies). The revolution of CCP inversion sidestepping the fixed point, understand why observed choice probabilities measure dynamic value.
- Magnac and Thesmar (2002, Econometrica). The identification analysis of dynamic discrete choice, read it to understand why the discount factor is not identified and why the normalization is a substantive restriction.
- Aguirregabiria and Mira (2002, Econometrica). Nested pseudo-likelihood (NPL), connecting Hotz-Miller one-step estimation and Rust NFXP into a spectrum, striking a middle ground between efficiency and computation.
- Aguirregabiria and Mira (2007, Econometrica). Extending two-step CCP to dynamic games.
- Gowrisankaran and Rysman (2012, Journal of Political Economy). A modern template for dynamic demand, how inclusive value sufficiency tames the curse of dimensionality, highly isomorphic to technology adoption on platforms.

Further reading:

- Su and Judd (2012, Econometrica). MPEC for dynamic models, another algorithm for the same estimator.
- Arcidiacono and Miller (2011, Econometrica). Finite dependence and CCP estimation with unobserved heterogeneity.
- Hendel and Nevo (2006, Econometrica). Dynamic demand for storable goods, consumer stockpiling and intertemporal price discrimination.
- Bajari, Benkard and Levin (2007, Econometrica). Two-step estimation of dynamic games, pushing the dynamic to multi-agent competition.
- Ericson and Pakes (1995, Review of Economic Studies). The Markov perfect equilibrium framework for dynamic oligopoly, the theoretical bedrock of dynamic games.
:::

::: {.apa-refs}
- Aguirregabiria, V., & Mira, P. (2002). Swapping the nested fixed point algorithm: A class of estimators for discrete Markov decision models. *Econometrica, 70*(4), 1519-1543. https://doi.org/10.1111/1468-0262.00340
- Aguirregabiria, V., & Mira, P. (2007). Sequential estimation of dynamic discrete games. *Econometrica, 75*(1), 1-53. https://doi.org/10.1111/j.1468-0262.2007.00731.x
- Arcidiacono, P., & Miller, R. A. (2011). Conditional choice probability estimation of dynamic discrete choice models with unobserved heterogeneity. *Econometrica, 79*(6), 1823-1867. https://doi.org/10.3982/ECTA7743
- Bajari, P., Benkard, C. L., & Levin, J. (2007). Estimating dynamic models of imperfect competition. *Econometrica, 75*(5), 1331-1370. https://doi.org/10.1111/j.1468-0262.2007.00796.x
- Ericson, R., & Pakes, A. (1995). Markov-perfect industry dynamics: A framework for empirical work. *The Review of Economic Studies, 62*(1), 53-82. https://doi.org/10.2307/2297841
- Gowrisankaran, G., & Rysman, M. (2012). Dynamics of consumer demand for new durable goods. *Journal of Political Economy, 120*(6), 1173-1219. https://doi.org/10.1086/669540
- Hendel, I., & Nevo, A. (2006). Measuring the implications of sales and consumer inventory behavior. *Econometrica, 74*(6), 1637-1673. https://doi.org/10.1111/j.1468-0262.2006.00721.x
- Hotz, V. J., & Miller, R. A. (1993). Conditional choice probabilities and the estimation of dynamic models. *The Review of Economic Studies, 60*(3), 497-529. https://doi.org/10.2307/2298122
- Magnac, T., & Thesmar, D. (2002). Identifying dynamic discrete decision processes. *Econometrica, 70*(2), 801-816. https://doi.org/10.1111/1468-0262.00306
- Rust, J. (1987). Optimal replacement of GMC bus engines: An empirical model of Harold Zurcher. *Econometrica, 55*(5), 999-1033. https://doi.org/10.2307/1911259
- Su, C.-L., & Judd, K. L. (2012). Constrained optimization approaches to estimation of structural models. *Econometrica, 80*(5), 2213-2230. https://doi.org/10.3982/ECTA7925
:::
