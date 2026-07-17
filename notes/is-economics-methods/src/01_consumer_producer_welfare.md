---
title: "Consumer, Producer & Welfare"
subtitle: "Duality, Demand, and the Measurement of Welfare"
seriesline: "Foundations of Information Systems Economics · Chapter 1"
series: "Foundations of Information Systems Economics: Empirical Methods for Information Systems Economics"
sidebartitle: "Foundations of Information Systems Economics"
sidebarsub: "Theory · Methods · Applications"
chapterlabel: "Chapter 1 · Consumer, Producer & Welfare"
---

## Introduction

A platform raises the price of a feature from 2 to 4. Users spend no more on that feature than they did before, so the first analyst declares the welfare loss to be zero. The second analyst integrates the demand curve and gets 27.73. Working more carefully still, the compensating variation and the equivalent variation come out to 31.95 and 24.21. All four numbers could go on a respectable slide. The trouble is that they are not answers to the same question.

Welfare is not a column that sits ready-made in the transaction log. To get from prices, quantities, and income to "how much did the consumer actually lose," you first have to say how a person chooses under a budget constraint, and then separate the substitution effect from the income effect inside a price change. That is exactly the job of consumer theory. Utility maximization and expenditure minimization describe one and the same choice from two directions, and duality is what ties Marshallian demand, Hicksian demand, indirect utility, and the expenditure function together. Roy's identity, Shephard's lemma, and the Slutsky equation are not a set of formulas to be memorized for an exam. They are the drivetrain that turns observable demand into welfare measures and empirical predictions. On the producer side, cost minimization and profit maximization follow the same structure, and choice under uncertainty extends the same language to adoption, investment, and risk decisions.

This chapter reconciles Lumen's price increase line by line. You will see why the three welfare numbers differ, when consumer surplus is only an approximation, and how to push theoretical comparative statics all the way into predictions you can take to data. Later chapters on demand estimation, market power, merger simulation, and the welfare analysis of digital goods borrow this machinery again and again. If we do not settle the accounting conventions here, the more precisely you compute later, the more you may just be reporting the wrong decimal places.

## 1 A Loss You Cannot Pin Down

::: {.case}
Lumen is a heavy user of a digital content platform. To shrink the problem down to something you can do by hand, we group her monthly spending into two categories. One is premium content (paid content time, call it good 1, unit price $p_1$, starting at 2 per hour), the other is an outside good that packages everything else she consumes (a composite outside good, call it good 2, unit price $p_2 = 1$, which you can read as "everything else you get for each dollar spent"). Her disposable budget each month is $w = 100$. Her preferences can be summarized by a Cobb-Douglas utility $u(x_1, x_2) = x_1^{2/5}\, x_2^{3/5}$. The virtue of this form is that every step can be checked with elementary algebra, and its drawback is discussed in Section 7.

Now the platform raises the price of premium content from 2 to 4. The operations team wants to know something that sounds simple: how much did this price increase cost Lumen?
:::

Start with the two answers two analysts might give. Neither is exactly wrong, yet they are wildly far apart.

The first analyst leafs through the transaction log. Before the increase, Lumen bought 20 hours of premium content a month, spending 40. After the increase she bought less, only 10 hours, but at double the price, still 40. Her spending on the feature did not move by a cent. So the conclusion is: her wallet is no lighter, the increase was revenue neutral for her, and there is no loss to speak of. The flaw in this reasoning is that it mistakes spending for welfare. Spending stayed flat precisely because she was forced to cut consumption sharply, offsetting the higher unit price by buying much less, and it is exactly that cut she has to endure.

The second analyst knows the trade better, and goes to integrate the demand curve. Lumen's demand for good 1 is $x_1 = (2/5)\,w / p_1 = 40 / p_1$. Integrating this curve as $p_1$ moves from 2 to 4 gives a consumer surplus loss of $\int_2^4 (40/p_1)\, dp_1 = 40\ln 2 \approx 27.73$. This is a proper welfare number, close in magnitude to three tenths of her monthly budget, and far more credible than the first analyst's zero.

The trouble is that there are a third and a fourth number, every bit as proper. If you ask "how much would we have to pay Lumen so that she is as well off after the increase as before," the answer is 31.95. If you ask "before the increase, how much could we take away from her to make her exactly as unhappy as the increase makes her," the answer is 24.21. The former is called compensating variation, the latter equivalent variation. Both are logically airtight welfare measures, yet they differ by nearly a third, with the 27.73 of consumer surplus sitting squarely between them. The phrase "how much did she lose" does not point to a unique number, unless we say which compensation we want and at which prices it is valued.

What this chapter builds is precisely the machine that makes clear what each of these numbers means exactly, how they rank against one another, and when they collapse into one. The machine is called duality, and its two engines are utility maximization and expenditure minimization. Section 6 recomputes every one of the numbers above from scratch, showing which part of the machine each falls out of.

## 2 The Economic Model and the Target Quantity

Before doing anything, pin down two things: which quantity we ultimately want, and what notation we use to write down the consumer's behavior.

Take the target quantity first. On the surface the operations team wants a welfare number, but welfare is not a primitive. It is derived from something more basic. The true primitives are preferences, and the choices those preferences generate under a budget constraint, that is, demand. Demand is observable (it is right there in the log), while welfare is not directly observable (no table records "how much the increase hurt her"). The logical spine of the whole chapter is therefore a bridge: from observable demand, back out unobservable welfare. Duality is the piers of that bridge. So the target quantity of this chapter really has two layers. The bottom layer is the demand function and the utility structure behind it, and the top layer is the welfare measures CV and EV built on top of it, which can only be computed once the bottom layer is in place.

Now the notation. A consumption bundle is $x = (x_1, x_2) \in \mathbb{R}^2_+$, prices are $p = (p_1, p_2) \gg 0$, and income is $w > 0$. The consumer's preferences are represented by a continuous, increasing, quasi-concave utility function $u(x)$. Increasing (more is better) guarantees she spends her whole budget, and quasi-concave (convex upper contour sets) guarantees that the first-order conditions characterize an optimum rather than a worst point. These two are the minimal regularity that makes all the differentiation below valid. We assume them throughout the chapter, and revisit what happens when they fail in Section 7.

Given these primitives, the consumer's problem is to pick the most preferred bundle among the affordable ones:

$$\max_{x \geq 0}\; u(x) \quad \text{s.t.} \quad p \cdot x \leq w$$

This is the utility maximization problem (UMP). Its solution is Walrasian demand (also called Marshallian demand) $x(p, w)$, and recording the utility value of that optimal bundle gives the indirect utility function $v(p, w) = u(x(p, w))$, which answers "at these prices and this income, what is the best she can do." The UMP is the engine on the demand side.

The welfare side needs another engine. To measure what a price change is worth in money, the natural question is "to hold some given level of welfare fixed, what is the cheapest way to pay for it." When prices change, that cheapest outlay changes, and the size of that change is the monetary value of the price change. This engine is the expenditure minimization problem (EMP):

$$\min_{x \geq 0}\; p \cdot x \quad \text{s.t.} \quad u(x) \geq \bar u$$

Its solution is Hicksian demand (also called compensated demand) $h(p, \bar u)$, and the minimal outlay itself is the expenditure function $e(p, \bar u) = p \cdot h(p, \bar u)$. Hicksian demand is called "compensated" because, as you move along it, utility is held fixed at $\bar u$, the real income effect of the price change is compensated away in lockstep, and what remains is pure substitution.

With the expenditure function in hand, the welfare measures have precise definitions. Let prices change from $p^0$ to $p^1$, with utility before the increase $u^0 = v(p^0, w)$ and after the increase $u^1 = v(p^1, w)$. The two monetized welfare measures are

$$\mathrm{CV} = e(p^1, u^0) - e(p^1, u^1) = e(p^1, u^0) - w, \qquad \mathrm{EV} = e(p^0, u^0) - e(p^0, u^1) = w - e(p^0, u^1)$$

CV is valued at the new prices $p^1$ and asks how much extra we must give her to compensate her back to the old utility $u^0$. EV is valued at the old prices $p^0$ and asks how much we could take away, absent the increase, to drop her to the new utility $u^1$. The two definitions differ only in which prices serve as the measuring stick, and it is precisely this difference in the ruler that makes them give different numbers in general. Section 4 sets them side by side with consumer surplus, and Section 6 plugs in Lumen's numbers.

Here is an idea that runs through the whole chapter, and through the whole series that follows: define the target quantity first, then talk about how to compute it. The operations team's "how much did she lose" must first be translated into "do you want CV or EV, valued at which prices," before computation even makes sense. Skip that step and integrate the demand curve directly, and whether the consumer surplus you get approximates CV or EV, and how well, is a muddle. This discipline of "estimand first" reappears with exactly the same face later in the series when we cover demand estimation and merger simulation, where the target quantity is the market-level change in consumer surplus, and the logic is identical to here.

To summarize: consumer theory has two engines, the UMP producing demand and indirect utility, the EMP producing compensated demand and the expenditure function. The welfare measures CV and EV are differences of two values of the expenditure function. They are the top-layer target quantities of this chapter, and to compute them you first need the bottom-layer demand and utility structure. The next section takes the two engines apart to see their internal structure and the drivetrain between them.

## 3 The Core Structure: UMP, EMP, and Duality

This is the most finely worked section of the chapter, and the load-bearing block in the theoretical foundation of the whole series. The goal is to work through the solutions and properties of the UMP and EMP, then build the duality relations between the two one by one, and finally land on three results that later chapters use over and over: the Slutsky equation, Roy's identity, and Shephard's lemma. The exposition goes from shallow to deep, giving each object a geometric or algebraic intuition first, then the formalization.

### 3.1 UMP: From Geometry to First-Order Conditions

Look at the UMP once on a two-dimensional picture. The budget constraint $p_1 x_1 + p_2 x_2 = w$ is a downward-sloping line with slope $-p_1/p_2$ and horizontal and vertical intercepts $w/p_1$ and $w/p_2$. The level curves of the utility function are a family of indifference curves, with utility rising toward the upper right. Picking the highest-utility point on the budget line is, geometrically, letting the budget line reach the highest indifference curve it can, which it does at a point of tangency. Tangency means the slope of the indifference curve equals the slope of the budget line, and the slope of the indifference curve is the marginal rate of substitution $\mathrm{MRS} = (\partial u/\partial x_1)/(\partial u/\partial x_2)$. So the characterization of an interior optimum is a very intuitive sentence: the subjective rate at which she is willing to trade good 2 for good 1 equals exactly the objective rate $p_1/p_2$ the market forces her to pay.

![The geometry of the UMP: the budget line is tangent to the highest attainable indifference curve at the optimal bundle $x(p,w)$, where the MRS equals the price ratio $p_1/p_2$.](assets/fig/fig_01_duality.svg)

Writing that sentence as algebra gives the Lagrangian first-order conditions. Form $\mathcal{L} = u(x) + \lambda(w - p \cdot x)$ and differentiate with respect to each $x_\ell$:

$$\frac{\partial u}{\partial x_\ell} = \lambda\, p_\ell, \qquad \ell = 1, 2$$

Dividing the two equations to cancel $\lambda$ returns $\mathrm{MRS} = p_1/p_2$. The multiplier $\lambda$ is not a pure mathematical accessory. It has a real economic meaning: $\lambda = \partial v(p,w)/\partial w$, the marginal utility of income, the extra utility one more dollar can buy. This identity is used when we derive Roy's identity in Section 3.4.

The solution of the first-order conditions, $x(p, w)$, is Walrasian demand, and it has three properties that matter, because they are used everywhere later. First, homogeneity of degree zero: $x(\alpha p, \alpha w) = x(p, w)$ for any $\alpha > 0$. Scaling prices and income together cannot move the real choice, because the budget set does not change at all, which is also why only relative prices matter in economics. Second, Walras' law: $p \cdot x(p, w) = w$. Under increasing preferences she spends every cent, and the budget constraint holds with equality. Third, when preferences are strictly quasi-concave and the solution is interior, $x(p,w)$ is single-valued and varies continuously with prices. Quasi-concavity alone guarantees only that the solution set is convex, and strict quasi-concavity (strictly convex upper contour sets) is needed to rule out the multivaluedness that a flat segment in an indifference curve would cause. Indirect utility $v(p,w) = u(x(p,w))$ is correspondingly increasing in $w$, decreasing in each $p_\ell$, homogeneous of degree zero, and quasi-convex in prices. The intuition for this last property is that the more extreme prices are, the worse off she is, and the harm is not averaged away by extreme price combinations.

### 3.2 EMP: The Expenditure Function and Compensated Demand

The EMP is the mirror image of the UMP: instead of spending a given amount of money to reach the highest utility, it holds a given utility fixed and finds the cheapest way to pay for it. Geometrically, fix an indifference curve $u(x) = \bar u$, and among all bundles that land on this curve (or above it), pick the one that makes the constant $k$ in an iso-expenditure line $p \cdot x = k$ smallest. The iso-expenditure line is also a line of slope $-p_1/p_2$, and pushing it outward from the origin until it just touches the indifference curve gives a tangency, whose point is Hicksian demand $h(p, \bar u)$. The tangency condition is again $\mathrm{MRS} = p_1/p_2$, exactly as in the UMP, and this is the root of duality: the two problems satisfy the same tangency condition at the optimum, only one treats the budget as the constraint and utility as the objective while the other reverses them.

The expenditure function $e(p, \bar u) = p \cdot h(p, \bar u)$ has a set of properties symmetric to those of $v$. It is increasing in each $p_\ell$ (when things get more expensive, holding utility fixed of course costs more), increasing in $\bar u$, homogeneous of degree one in prices ($e(\alpha p, \bar u) = \alpha\, e(p, \bar u)$, so doubling all prices doubles the cost of holding utility fixed), and concave in prices. Concavity is the least obvious property but the most useful. The intuition is that when a price rises, the consumer substitutes consumption toward goods that got less expensive, so the actual rise in expenditure is smaller than the linear extrapolation "with the bundle held fixed," and concavity is the mathematical expression of exactly this "substitution saves you money." This concavity translates directly into the sign of the substitution effect in Section 3.5.

Hicksian demand itself has its signature property, also stemming from concavity: it is nonincreasing in its own price, $\partial h_\ell/\partial p_\ell \leq 0$, which is the so-called compensated law of demand, and compensated demand slopes down without exception. Note that no "normal good" qualification is needed here, because utility is held fixed, the income effect is compensated away, and the pure substitution that remains must make one buy less of what got more expensive. Walrasian demand is different, since it can slope up because of a strong income effect (a Giffen good), and the whole difference is whether or not there is compensation.

### 3.3 The Four Duality Identities

Now connect the outputs of the two engines. Under continuous, increasing, quasi-concave preferences, there are four groups of identities among the solutions and value functions of the UMP and EMP that pin the four objects $x, v, h, e$ down to one another pairwise. Let $u^* = v(p, w)$. Then

$$e(p, v(p, w)) = w, \qquad v(p, e(p, \bar u)) = \bar u$$

$$x_\ell(p, w) = h_\ell(p, v(p, w)), \qquad h_\ell(p, \bar u) = x_\ell(p, e(p, \bar u))$$

The first row says $v$ and $e$ are inverses of each other (fix prices, convert utility into money and money back into utility, and return to where you started). The second row says Walrasian and Hicksian demand take the same value at the "matching income-utility point": if your income is exactly the cheapest outlay that holds utility $u^*$, then your uncompensated choice and your compensated choice are the same bundle. These four identities are not coincidences. They are the inevitable result of one tangency described in two ways. Their use is that the property of any one object can be carried to another object through these identities, and the three signature results in the next three subsections are all products of this carrying.

::: {.intuition}
Why is duality worth all the trouble? Because it swaps "quantities that are hard to find" for "quantities that are easy to find." Welfare wants a difference of $e$, but $e$ is defined over unobservable utility, whereas demand $x$ is observable. Of the four identities, $h_\ell = x_\ell(p, e(p,\bar u))$ is exactly the channel that lets us assemble compensated demand from observed Walrasian demand, and from there the expenditure function and welfare. Later, when we cover demand estimation, what is estimated is $x$ and what is wanted is welfare, and the procedure in between is the duality here. So this section is not theoretical fastidiousness, it is the construction blueprint for the whole of empirical welfare analysis.
:::

### 3.4 Roy's Identity and Shephard's Lemma

Both results are direct applications of the envelope theorem, so give the intuition for the envelope theorem first and then apply it. The envelope theorem says that once you already stand at an optimum, the effect of a small change in a parameter on the optimal value only requires the direct effect of the parameter, and you need not worry about how you would reoptimize, because at the optimum the first-order gain from reoptimizing is zero (guaranteed by the first-order conditions).

Shephard's lemma is the envelope theorem applied to the EMP. Differentiating the expenditure function $e(p, \bar u) = \min_x \{p \cdot x : u(x) \geq \bar u\}$ with respect to price $p_\ell$, the direct effect is the quantity held $x_\ell$, and the indirect effect of the consumer readjusting the bundle is zero at the optimum, so

$$\frac{\partial e(p, \bar u)}{\partial p_\ell} = h_\ell(p, \bar u)$$

In one sentence: the slope of the expenditure function with respect to a price is the compensated demand for that good. The result looks plain, and its power is that it turns the abstract object "compensated demand" into a single derivative of the expenditure function, which is how Section 6 computes Lumen's Hicksian demand. Its twin on the producer side (the derivative of the cost function with respect to a factor price gives conditional factor demand) is the star of Section 4.2.

Roy's identity is the envelope theorem applied to the UMP, but with one more step. Differentiating indirect utility $v(p, w)$ with respect to $p_\ell$, by the envelope theorem only the direct effect $-\lambda x_\ell$ remains (a price rose, so at the original bundle she spends $x_\ell\,dp_\ell$ more, at $\lambda$ utility lost per dollar), and differentiating with respect to $w$ gives $\partial v/\partial w = \lambda$. Dividing the two cancels the unobservable multiplier $\lambda$:

$$x_\ell(p, w) = -\frac{\partial v(p, w)/\partial p_\ell}{\partial v(p, w)/\partial w}$$

In one sentence: demand equals the gradient of indirect utility with respect to prices divided by the marginal utility of income, with the minus sign coming from the fact that a price increase is a bad thing. The practical value of Roy's identity is in running it backwards: if you have a functional form for $v$ (many demand systems specify $v$ first and derive demand), one differentiation gives you all the demand functions, without solving the UMP again.

### 3.5 The Slutsky Equation: Splitting the Slope of Demand in Two

Now to the centerpiece of the section. The response of Walrasian demand to price, $\partial x_\ell/\partial p_k$, is observable, and it is exactly what demand estimation spits out directly, but it mixes together two things of different natures: when a price rises, one thing is that relative prices change and prompt substitution (the substitution effect), and the other is that real purchasing power shrinks (the income effect). The Slutsky equation separates these two halves precisely.

The derivation starts from the duality identity $h_\ell(p, \bar u) = x_\ell(p, e(p, \bar u))$ and differentiates with respect to $p_k$ (chain rule, noting that $e$ on the right also contains $p_k$):

$$\frac{\partial h_\ell}{\partial p_k} = \frac{\partial x_\ell}{\partial p_k} + \frac{\partial x_\ell}{\partial w}\cdot\frac{\partial e}{\partial p_k}$$

Then use Shephard's lemma to replace $\partial e/\partial p_k$ with $h_k = x_k$ (the two are equal at the matching point), and rearranging gives the Slutsky equation:

$$\underbrace{\frac{\partial x_\ell(p, w)}{\partial p_k}}_{\text{total}} = \underbrace{\frac{\partial h_\ell(p, \bar u)}{\partial p_k}}_{\text{substitution}}\; \underbrace{-\; x_k(p, w)\,\frac{\partial x_\ell(p, w)}{\partial w}}_{\text{income}}$$

Read the terms as follows: on the left the total effect (observable), on the right the substitution effect first, then the income effect (together with the minus sign in front, written $-x_k\,\partial x_\ell/\partial w$).

This equation explains a lot. Take the own-price case ($\ell = k$): the substitution term $\partial h_\ell/\partial p_\ell \leq 0$ is always negative (the compensated law of demand from Section 3.2), so whether Walrasian demand slopes down at all comes entirely down to the income term $-x_\ell\,\partial x_\ell/\partial w$. For a normal good ($\partial x_\ell/\partial w > 0$) the income term is also negative, and demand is guaranteed to slope down. For an inferior good ($\partial x_\ell/\partial w < 0$) the income term turns positive, and if it is strong enough to overwhelm the substitution term, you get an upward-sloping Giffen good. This explains a question that often puzzles undergraduates: why does the law of demand have exceptions, and what exactly are the conditions for the exception, with the answer written entirely in this decomposition.

The substitution effect has two more structural properties, both stemming from the fact that it is really a second derivative of the expenditure function, $\partial h_\ell/\partial p_k = \partial^2 e/\partial p_\ell \partial p_k$. First, symmetry $\partial h_\ell/\partial p_k = \partial h_k/\partial p_\ell$ (mixed partials commute), which is completely invisible in Walrasian demand and is a gift from duality, and is also the core testable implication for judging, in Section 7, "whether a set of demand functions really comes from some rational preference." Second, the Slutsky matrix $\big[\partial h_\ell/\partial p_k\big]$ is negative semidefinite with nonpositive diagonal, which translates "substitution saves you money" concavity into the language of matrices.

### 3.6 Revealed Preference: Backing Preferences Out of Demand

Everything up to here has run forward: preferences first, then solve for demand. But the empirical worker gets the reverse end of things. What is observed is choices (the prices and quantities in the log), and what is wanted is preferences and the welfare hidden behind them. So an unavoidable question surfaces: on what grounds can a pile of observed choices be said to come from some rational preference? And if it can, how much can we back out from it? This is what revealed preference and integrability answer, and it is the identification foundation for all the demand estimation that follows, and the flip side that this section naturally turns up from forward duality.

Start with the finite-observation version. Suppose we see a consumer pick a bundle in each of several price-income situations. If at some point she faced prices $p^0$ and picked $x^0$ while another bundle $x^1$ was clearly affordable at the time ($p^0 \cdot x^1 \leq p^0 \cdot x^0$) yet was not picked, then rationality requires: whenever $x^0$ is affordable, $x^1$ should not be picked in its place. In terms of prices, as long as $x^1 \neq x^0$ and $p^0 \cdot x^1 \leq p^0 \cdot x^0$, we must have $p^1 \cdot x^0 > p^1 \cdot x^1$. This is the weak axiom of revealed preference (WARP). Its geometric meaning is plain: since you insisted on $x^0$ when $x^1$ was within reach, you value $x^0$ more, and from then on, unless $x^0$ really becomes unaffordable, you cannot turn around and pick $x^1$. Passing this "directly revealed preference" along a chain, and forbidding cycles ($x^0$ revealed preferred to $x^1$, $x^1$ to $x^2$, and $x^2$ back to $x^0$), gives the strong axiom of revealed preference (SARP). SARP is the rationalizability condition Houthakker (1950) gave for single-valued demand. To handle the more general case of finite observations (where demand may be multivalued and indifference curves may contain flat segments), the correct necessary and sufficient condition is the generalized axiom of revealed preference (GARP): Afriat (1967) and Varian (1982) proved that a set of finite observations can be rationalized by some well-behaved preference (locally nonsatiated, continuous, concave, monotone) if and only if it satisfies GARP. GARP and SARP differ only in whether the end of the chain takes a weak or a strict inequality, and when demand is single-valued the two are equivalent, which is where the running Cobb-Douglas case of this chapter falls. With two goods WARP already suffices, and it is only when there are more goods, where WARP cannot block cycles, that you need SARP or GARP.

The version for a continuously differentiable demand function is called integrability, and it translates the combinatorial condition above into a requirement on derivatives. The question is posed this way: given a demand function $x(p, w)$ that satisfies Walras' law and homogeneity of degree zero, when is it exactly the UMP solution of some rational preference? The answer is astonishingly clean: if and only if its Slutsky matrix $S(p,w) = \big[\partial x_\ell/\partial p_k + x_k\,\partial x_\ell/\partial w\big]$ is symmetric and negative semidefinite. Negative semidefiniteness we already know from Section 3.5 as a corollary of the substitution effect "sloping down," and what really brings the binding force is symmetry $\partial h_\ell/\partial p_k = \partial h_k/\partial p_\ell$, which is completely invisible in Walrasian demand and is a testable restriction that rationality imposes on the data. Once these two hold, you can integrate the expenditure function back: Shephard's lemma says $\partial e/\partial p_\ell = h_\ell$, which is a system of partial differential equations with $e$ as the unknown, and symmetry is exactly the consistency condition for it to have a solution (to be integrable), and solving out $e(p, \bar u)$ is equivalent to recovering preferences. This step is recoverability, and it is the fundamental reason welfare analysis can hold up: welfare wants $e$, and $e$ is defined over unobservable utility, and integrability guarantees that as long as demand is estimated well enough, $e$ can be recovered from it (up to at most an irrelevant utility scale). The "assemble the expenditure function from demand" done for Lumen in Section 6.2 is one concrete execution of this recovery.

::: {.intuition}
Revealed preference gives demand estimation two things. One is falsifiability: Slutsky symmetry and negative semidefiniteness are restrictions the data could have violated but, because of rationality, should not, which turns "the consumer is rational" from a slogan into a proposition testable against data. The other is the boundary of identification: all that can be recovered is ordinal preference (the shapes of the indifference curves), and the specific scale of utility is smoothed away by monotone transformations anyway, cannot be recovered, and need not be, because the welfare measures CV and EV depend only on $e$ and not on the scale. The legitimacy of the "estimate demand first, then recover welfare via duality" that recurs later in demand estimation is anchored right here.
:::

### 3.7 From Individual to Market: Aggregation and the Representative Consumer

Every step so far has been about a single consumer, but demand estimation and policy evaluation almost always play out at the market level, with market demand, the sum over everyone, in hand. So another unavoidable question: when you add up the demand of a crowd of people, does the resulting total demand still look like it came from a single rational consumer? If it does, we can legitimately model the whole market with one representative consumer, sparing ourselves the trouble of tracking each person. If it does not, doing so will compute welfare wrong.

The trouble is obvious at a glance. Market demand $X_\ell(p, w_1, \dots, w_I) = \sum_i x_{\ell i}(p, w_i)$ in general depends not only on total income $W = \sum_i w_i$ but also on how income is distributed across people: move money from one person's hands to another's, and even if $W$ is unchanged, as long as the two people have different marginal propensities to consume the good, total demand will move. To make total demand recognize only $W$ and not the distribution, so that it can be rationalized by a representative consumer, requires that every person spend each marginal dollar of income across goods in the same proportions, independent of the income level. In terms of Engel curves: each person's demand must be linear in income, with slopes that are the same across people and intercepts that may differ.

Writing this condition as a value function gives the Gorman polar form: each consumer's indirect utility can be written as

$$v_i(p, w_i) = a_i(p) + b(p)\, w_i$$

where $b(p)$ is the same for everyone and $a_i(p)$ may vary by person. Why exactly this form works is checked in one step with Roy's identity. Apply the Roy's identity of Section 3.4 to $v_i$:

$$x_{\ell i}(p, w_i) = -\frac{\partial v_i/\partial p_\ell}{\partial v_i/\partial w_i} = -\frac{\partial a_i/\partial p_\ell}{b(p)} - \frac{\partial b/\partial p_\ell}{b(p)}\, w_i$$

The key is the last term: the slope of $x_{\ell i}$ in $w_i$ is $-(\partial b/\partial p_\ell)/b(p)$, which depends only on the common $b(p)$ and is independent of $i$. So when we aggregate, that pile of person-specific intercepts each stands on its own and can be merged into a single constant, and the part that varies with income recognizes only the total:

$$X_\ell(p, W) = -\frac{\sum_i \partial a_i/\partial p_\ell}{b(p)} - \frac{\partial b/\partial p_\ell}{b(p)}\, W$$

Total demand indeed depends only on $(p, W)$, and it is exactly the Walrasian demand of the representative consumer whose indirect utility is $v(p, W) = \sum_i a_i(p) + b(p)\,W$. Several familiar faces fall inside the Gorman form: homothetic preferences ($a_i \equiv 0$, in which case even the income distribution is irrelevant), quasi-linear preferences (linear in some common numeraire), and the linear expenditure system (subsistence demands may differ across people, but the marginal budget share is the same). Section 6.6 uses this last class to work out by hand an example where "redistribution does not change total demand."

::: {.warning}
The crux of the Gorman condition is the "common slope $b(p)$," that is, everyone's marginal budget share must agree, not that price sensitivity is free to be heterogeneous. In reality consumers' responses to prices differ enormously, and that difference shows up precisely in the slope rather than the intercept, and once the slope varies by person, aggregate demand can no longer be rationalized by any single rational consumer, and forcing a representative consumer on it systematically distorts elasticities and welfare. This is precisely why the demand estimation later in this series would rather spend great effort modeling consumer heterogeneity (moving toward mixed logit and random coefficients) than take the easy road of a representative consumer.
:::

To summarize: this section started from two engines, used four duality identities to stitch demand, indirect utility, compensated demand, and the expenditure function into one, then used the envelope theorem to squeeze out Roy's identity and Shephard's lemma, and used one chain-rule differentiation to get the Slutsky equation, splitting the observable demand slope into the compensated substitution and the purchasing-power income halves. The two subsequent subsections used this machine in reverse and sideways: revealed preference and integrability made clear the conditions (symmetry plus negative semidefiniteness) for backing preferences out of observed demand, and thence recovering the expenditure function and welfare, and aggregation and the Gorman form made clear the legitimate boundary for summing individual demand into market demand and modeling it with a representative consumer. Every part of this machine has a corresponding use in later empirical chapters: demand estimation estimates the total effect, welfare computation wants the expenditure function behind the substitution effect, the necessity of modeling heterogeneity hinges on the failure of the Gorman condition, and duality is always the drivetrain between them. The next section drives this machine off to compute welfare, and copies the same logic onto the producer and onto uncertainty.

## 4 Welfare Measures, and Two Extensions of Duality

The last section built the machine, and this section puts it to work on three jobs: measuring consumer welfare, characterizing the producer, and handling uncertainty. The three jobs share the same optimization-plus-duality skeleton, and lining them up makes clear they are really one idea performing three times.

### 4.1 CV, EV, and Consumer Surplus

Return to the CV and EV defined in Section 2. Let prices change from $p^0$ to $p^1$ (in this chapter's example only $p_1$ moves). Shephard's lemma lets us write these two expenditure differences as areas under the compensated demand curve, and this step is the key to turning abstract definitions into computable integrals. Because $\partial e/\partial p_1 = h_1$, integrating along $p_1$ from $p_1^0$ to $p_1^1$:

$$\mathrm{CV} = \int_{p_1^0}^{p_1^1} h_1(p_1, p_{-1}, u^0)\, dp_1, \qquad \mathrm{EV} = \int_{p_1^0}^{p_1^1} h_1(p_1, p_{-1}, u^1)\, dp_1$$

Both are "the area under some Hicksian demand curve, between two prices," differing only in whether the utility pinned down is the old $u^0$ or the new $u^1$. And the change in consumer surplus everyone is familiar with, $\Delta \mathrm{CS} = \int_{p_1^0}^{p_1^1} x_1(p_1, p_{-1}, w)\, dp_1$, uses Walrasian demand. By the duality identity $h_\ell(p, u) = x_\ell(p, e(p, u))$, the Walrasian curve intersects some Hicksian curve only at the price where expenditure equals exactly $w$: it intersects the Hicksian pinned at $u^0$ at the initial price $p^0$, and the one pinned at $u^1$ at the final price $p^1$, and the three curves do not share a common point, and the two Hicksian curves never intersect at all, so the three areas are in general unequal.

![The three demand curves for good 1: the $u^0$ and $u^1$ Hicksian (compensated) demands, and the Walrasian demand in between; as $p_1$ rises from 2 to 4, the areas under the three curves give CV, EV, and $\Delta$CS respectively.](assets/fig/fig_01_cv_ev.svg)

For a single price increase, with good 1 a normal good, the three have a definite ranking $\mathrm{EV} \leq \Delta\mathrm{CS} \leq \mathrm{CV}$. The intuition is this: for a normal good the Walrasian demand is flatter than the Hicksian (the income effect amplifies the price response), and CV is pinned at the higher old utility $u^0$, whose Hicksian curve sits farthest out, while EV is pinned at the lower new utility $u^1$, sitting farthest in, with Walrasian squeezed in between. Consumer surplus is therefore an approximation lying between two exact welfare measures, and its convenience (only observable Walrasian demand is needed, with no knowledge of utility) has a cost: as long as the income effect is nonnegligible, it equals neither CV nor EV. Section 6 shows that in Lumen's example these three numbers are 24.21, 27.73, 31.95, ranked without the slightest deviation.

When do the three collapse into one? When good 1 has no income effect, the two Hicksian curves coincide with the Walrasian curve into a single curve, and $\mathrm{CV} = \Delta\mathrm{CS} = \mathrm{EV}$. To make good 1, the good whose price moves, free of income effects, preferences must be quasi-linear in the outside good, that is, linear in good 2, $u = \phi(x_1) + x_2$ (quasi-linear, in which case the first-order condition $\phi'(x_1) = p_1/p_2$ pins down $x_1$ alone, without $w$, and all income variation is absorbed by good 2). Note the direction of quasi-linearity: the linear term goes on the outside good whose price does not move, so that it is the good 1 whose price does move that is spared income effects. This explains why consumer surplus can be used with confidence in undergraduate textbooks: they default to quasi-linearity in the numeraire. In reality, when good 1 takes up only a small slice of the budget and the income effect is weak, this approximation is also good enough. The criterion for whether consumer surplus can be used is a single one: whether the income effect of good 1 can be neglected.

### 4.2 The Producer Side: The Duality of Cost and Production

Producer theory is almost a word-for-word translation of consumer theory: replace utility with output, income with cost, the consumer with the firm, and the whole of duality carries over unchanged. The firm uses inputs $z = (z_1, z_2)$ to make output via a production function $f(z)$, with factor prices $r = (r_1, r_2)$. Given the output $y$ to be produced, the cost minimization problem

$$c(r, y) = \min_{z \geq 0}\; r \cdot z \quad \text{s.t.} \quad f(z) \geq y$$

is the same mathematical problem as the EMP: replace $u$ with $f$, $\bar u$ with $y$, and $p$ with $r$, word for word. The solution $z(r, y)$ is conditional factor demand, corresponding to Hicksian demand, and the value function $c(r, y)$ is the cost function, corresponding to the expenditure function. So all the properties of the EMP hold unchanged: $c$ is homogeneous of degree one in $r$, concave, and increasing, and Shephard's lemma here reads

$$\frac{\partial c(r, y)}{\partial r_\ell} = z_\ell(r, y)$$

The slope of the cost function with respect to a factor price is the conditional demand for that factor. This is extremely important in structural IO: later in the series, when we recover marginal cost and estimate production functions, we repeatedly use this layer of duality between the cost function and factor demand.

The firm's problem one level up is profit maximization. A price-taking firm takes the output price $p$ and factor prices $r$ as given and chooses output to maximize profit, which, with the help of the cost function, can be written as a one-dimensional problem

$$\pi(p, r) = \max_{y \geq 0}\; \big[\, p\, y - c(r, y)\,\big]$$

The first-order condition for an interior solution is $p = \partial c(r, y)/\partial y$, that is, price equals marginal cost, the old refrain of a competitive firm's supply decision, only here derived in one step from the cost function. The second-order condition requires marginal cost to be increasing, that is, cost to be convex in output, so that the supply curve slopes up and the optimum is an optimum and not a worst point. This $\pi(p, r)$ is the profit function, playing the role that corresponds to indirect utility $v$ on the consumer side: both are value functions obtained by substituting the optimal decision back into the objective.

The envelope result for the profit function is Hotelling's lemma, the twin of Roy's identity and Shephard's lemma on the producer side, and even cleaner than Roy's identity, because here no multiplier need be canceled. Applying the envelope theorem to $\pi(p, r) = \max_{y, z}[\, p f(z) - r \cdot z\,]$, a small change in the price and factor prices leaves only the direct effect:

$$\frac{\partial \pi(p, r)}{\partial p} = y(p, r), \qquad \frac{\partial \pi(p, r)}{\partial r_\ell} = -\, z_\ell(p, r)$$

In one sentence: differentiating the profit function with respect to the output price gives supply, and with respect to a factor price gives negative factor demand (the true unconditional factor demand, without an output constraint, distinct from the conditional factor demand in cost minimization that holds $y$ fixed). The profit function also has a set of properties symmetric to, but opposite in sign from, the expenditure function: it is homogeneous of degree one in $(p, r)$ (scaling all prices scales profit in proportion), and convex in $(p, r)$. Convexity is the most useful, and its intuition mirrors the concavity of the expenditure function exactly: when the price environment changes, the firm reallocates output and inputs in a favorable direction, so profit rises by more than the linear extrapolation "with the decision held fixed." Taking second derivatives of convexity in $p$ and in $r_\ell$ separately immediately gives two comparative statics: $\partial y/\partial p \geq 0$, supply is nondecreasing in its own price, and $\partial z_\ell/\partial r_\ell \leq 0$, factor demand is nonincreasing in its own price. The latter is completely isomorphic to the compensated law of demand for Hicksian demand, both being direct consequences of the direction the value function bends.

The last piece is returns to scale, which determines what the cost function looks like and thus whether competitive supply is well defined, a point that bears directly on the later chapters on market structure. If the production function is homogeneous of degree one (constant returns to scale, CRS, doubling all inputs doubles output), then the cost function is linear in output, $c(r, y) = y\cdot c(r, 1)$, with marginal cost identically equal to average cost and independent of output. In that case the price-taking firm's profit is unbounded when $p$ is above marginal cost, zero when below (shut down and produce nothing), and exactly zero for any output when equal, so the scale of a single firm cannot be pinned down and the number and size of firms in the market are indeterminate. If returns are decreasing, cost is convex in output, marginal cost is increasing, the supply curve slopes up, and firm scale has an interior solution. If returns are increasing (increasing returns, common in digital products with large fixed costs and near-zero marginal cost), average cost falls with output, competitive equilibrium usually fails to exist, and the market naturally tips toward a few large firms, which is exactly the technological root of why platform and network industries later discussed are naturally concentrated. Returns to scale, you can see, are not a detail of the production function but the hinge connecting firm technology to market form. Producer theory has no new mathematics, only new names for old mathematics, yet this naming runs all the way to market structure, which is itself proof of the power of duality.

### 4.3 Choice Under Uncertainty

So far all choices have been made in a deterministic environment, but most decisions on a platform are not: whether a seller should switch to a new pricing scheme with uncertain payoffs, whether a user should prepay an annual fee for a feature they might not use, both carry randomness with them. The standard tool for characterizing such choices is expected utility. von Neumann and Morgenstern proved that as long as preferences over lotteries (probability distributions over outcomes) satisfy the four axioms of completeness, transitivity, continuity, and independence (the independence axiom), there exists a Bernoulli utility function $u(\cdot)$ such that preferences over any lotteries can be ranked by its expected utility $\mathbb{E}[u(\cdot)]$. The point is that the expectation is taken over $u$, not over money, and the shape of $u$ therefore encodes the attitude toward risk.

The concavity or convexity of $u$ is exactly the risk attitude. If $u$ is concave, then by Jensen's inequality $\mathbb{E}[u(w)] \leq u(\mathbb{E}[w])$, and the expected utility of a lottery is below the utility of its expected value, so this person would rather have the certain expected value than the lottery, that is, risk averse. How to quantify the degree of concavity? Applying an affine transformation to $u$ (multiply by a positive number, add a constant) does not change preferences, so one cannot use $u''$ directly, and must use a quantity immune to affine transformations, which is the Arrow-Pratt coefficient of absolute risk aversion

$$A(w) = -\frac{u''(w)}{u'(w)}$$

The $u'$ in the denominator exactly cancels the scaling from the affine transformation, so that $A(w)$ reflects only curvature and not the meaningless scale. Two quantities that bring risk attitude down to money are the certainty equivalent $\mathrm{CE}$, that is, "how much certain money she is willing to trade for the lottery," defined by $u(\mathrm{CE}) = \mathbb{E}[u(w)]$, and the risk premium $\pi = \mathbb{E}[w] - \mathrm{CE}$, that is, the expected return she is willing to give up to avoid the risk. For a small zero-mean risk $\widetilde\varepsilon$ added to wealth $w$, the risk premium and $A(w)$ have a lovely approximation $\pi \approx \tfrac{1}{2}A(w)\,\mathrm{Var}(\widetilde\varepsilon)$: the risk premium is roughly half the coefficient of absolute risk aversion times the variance of the risk, and the more afraid of risk ($A$ large) and the larger the risk (large variance), the more compensation is wanted. Section 6 uses a seller's example to compute CE, $\pi$, and this approximation and reconcile them.

Arrow-Pratt compares one person's attitude toward different risks, and there is a complementary question: can we, without specifying a particular $u$, conclude that everyone of a certain class will prefer one lottery to another? This is what stochastic dominance answers, giving ranking criteria that require no knowledge of the specific utility, layered by the shape of $u$. Let the payoff distribution functions of two lotteries be $F$ and $G$. The first layer is first-order stochastic dominance (FOSD): if $F(t) \leq G(t)$ for all $t$, then $F$ first-order dominates $G$, meaning $F$ presses more probability onto the high-payoff side. It can be shown that $F$ first-order dominates $G$ if and only if every expected-utility maximizer with more-is-better preferences ($u$ increasing) weakly prefers $F$. The second layer relaxes the condition for risk-averters: when the two lotteries have the same mean, if $\int_{-\infty}^{t}\!F(s)\,ds \leq \int_{-\infty}^{t}\!G(s)\,ds$ for all $t$, then $F$ second-order dominates $G$ (second-order stochastic dominance, SOSD), meaning $G$ is a mean-preserving spread of $F$ (probability spread out toward the two ends, more "dispersed"). Rothschild and Stiglitz proved that $F$ second-order dominates $G$ if and only if every expected-utility maximizer with concave $u$ (risk averse) weakly prefers $F$. The two layers of dominance give welfare comparison a handy tool: as long as a dominance relation can be verified, one can rank risky prospects without tangling with the specific utility form. Section 6.6 uses a pair of equal-mean lotteries to match up the SOSD criterion with the ordering of expected utility.

Finally, point out two extensions of the expected utility framework that matter in settings like insurance, health, and platform risk control. The first is state-dependent utility. Standard vNM assumes that what an outcome is worth is independent of which state produces it, but reality may be otherwise: a dollar in health and a dollar in illness need not have the same utility, and in that case utility should be written as varying with the state, $\sum_s \pi_s\, u_s(x_s)$, with the same wealth measured by a different $u_s$ in different states, which is key to understanding why people buy "seemingly unprofitable" insurance. The second is subjective expected utility. Earlier we treated the probability $\pi_s$ as an exogenously given objective number, but many business judgments have no ready probability at all, only personal belief. Savage proved that when preferences are defined over "acts" (maps from states to outcomes) and satisfy a set of axioms, one can simultaneously derive a set of subjective probabilities and a utility function from the preferences, so that choice behaves as expected utility taken over subjective probabilities. This extends choice under uncertainty from "probability known" to "probability also revealed by behavior," and is the theoretical starting point for later discussions of belief, learning, and the value of information.

To summarize: this section had the duality machine do three jobs in a row. Consumer welfare was written as areas under compensated demand curves, with CV, EV, and consumer surplus ordered by the size of the income effect and collapsing into one under quasi-linearity. Producer theory was recognized as a word-for-word translation of consumer theory, with cost-production duality isomorphic to utility-expenditure duality, and the convexity of the profit function, Hotelling's lemma, and returns to scale connecting firm technology all the way to market form. Choice under uncertainty was characterized by expected utility, with risk attitude measured by the curvature of the Bernoulli utility, that is, the Arrow-Pratt coefficient, stochastic dominance giving rankings that require no specification of utility, and state-dependence and subjective expected utility extending the framework to cases closer to reality. The common backdrop of all three jobs is "optimize under a constraint, then read the wanted quantity off the derivative of the value function," which is exactly the main thread of the whole chapter.

### 4.4 How Comparative Statics Become Empirical Predictions

The previous sections built the machine and put it to use, and this section explains how it finally lands on data. Theory is useful for empirical work not because it tells a story for some single coefficient, but because it simultaneously constrains the relationships among several observable quantities, and it is this joint restriction that data can bring to bear as a test. Take the Slutsky equation first. The Marshallian demand response to a price increase is, by construction, the sum of a substitution effect and an income effect:

$$\frac{\partial x_j}{\partial p_k}=\frac{\partial h_j}{\partial p_k}-x_k\frac{\partial x_j}{\partial w}.$$

For a normal good, the own income effect is positive, so when the price rises, its Marshallian own-price response is more negative than the compensated response. This comparison, derived purely from signs, gives a prediction you can take to data: for the same price increase, the uncompensated demand response should be larger among those for whom the good takes up a higher budget share and the income effect is stronger. Beware, though: even if you really do see this cross-population heterogeneity in the data, it is only consistent with the income-effect mechanism, and does not automatically rule out that liquidity constraints, attention, or measurement error differing across groups could produce the same pattern. A matching sign is a necessary condition, not proof of the mechanism.

The duality on the producer side likewise presses out a set of joint restrictions. By Shephard's lemma, the conditional factor demand when a factor price $r_k$ rises is $\partial c/\partial r_k$, and the cost function is concave in factor prices, so the Jacobian of compensated factor demand must be symmetric and negative semidefinite, both of which you can check against an estimated cross-substitution matrix. If the estimated matrix badly violates symmetry, that is most likely not something a phrase like "firms have frictions" can wave away, and what should be suspected instead is aggregation at the technology level, adjustment costs, endogenous procurement prices, or a misspecified functional form, each to be ruled out in turn.

Under uncertainty, the Arrow-Pratt absolute risk aversion $A(w)=-u''(w)/u'(w)$ gives the risk premium approximation for a small risk:

$$\pi(w)\approx \frac{1}{2}A(w)\operatorname{Var}(\widetilde\varepsilon).$$

This approximation also spits out testable comparative statics: the larger the variance of the risk, the higher the risk premium, and if absolute risk aversion falls with wealth (decreasing absolute risk aversion), then rising wealth lowers the premium demanded for the same absolute risk. These predictions can themselves be falsified by data, but when you actually test them with insurance purchases or conservative choices, you must separate risk preferences, beliefs, and constraints, otherwise a behavior that "looks more risk averse" may just be the result of different beliefs or borrowing limits.

| Theoretical change | Comparative static | Observable implication | Alternative explanations not identified by sign alone |
|---|---|---|---|
| Own price rises | Hicksian demand has only substitution; Marshallian also has an income effect | Groups with a higher budget share or a stronger income response react more | liquidity, search, salience |
| Input price rises | Conditional input demand moves along the cost gradient | Conditional demand for substitutable inputs adjusts in the opposite direction | technical change, quantity discounts, endogenous procurement price |
| Risk variance rises | Risk premium increases with $A(w)$ | Risk-averters demand higher compensation | beliefs, ambiguity, borrowing constraints |
| Digital good's price falls to zero | The gap among CV, EV, and $\Delta \mathrm{CS}$ depends on income effects | High users' welfare gains are larger, but usage time alone cannot value it | selection into use, zero-price privacy cost |

To gather this section into one sentence: empirical research should treat the predictions above as a set of joint restrictions, and make clear what exogenous variation is driving prices, income, risk, or factor costs. Merely seeing the sign go the right way in a cross section is at best evidence consistent with a mechanism, and to really identify the mechanism you need price experiments, cost shifters, randomized information interventions, or other variation that can separate the competing alternative explanations.

## 5 Anchoring Papers

Theory only stands up once it lands in real measurement. This section picks two studies that put this chapter's welfare machine to actual use, one the methodological origin of valuing new goods, and one that carries the same logic over to digital products, each laid out along five elements (paper, method, data, results, limitations), with a focus on how the welfare measure is computed out of demand estimation.

### 5.1 Hausman (1996)

::: {.case}
Paper and methodological position: "Valuation of New Goods under Perfect and Imperfect Competition," in The Economics of New Goods, edited by Bresnahan and Gordon (NBER/University of Chicago Press). It turns the CV computation of Section 4.1 of this chapter into the standard way to value the social worth of a new product, and points out the problem that official price indices systematically overstate inflation by ignoring new goods.

Method: a new good coming into existence from nothing is equivalent to its price falling from a "reservation price so high nobody buys" to the actual market price. The consumer welfare brought by the new good is the area under the demand curve over this price interval, that is, this chapter's CV. To compute that area, you must first estimate the entire demand curve (not just a single point at the current price), and Hausman estimates the demand for breakfast cereal with an almost ideal demand system, using prices in other cities as instruments to handle the endogeneity of price when identifying price elasticities.

Data: Nielsen supermarket scanner data, focused on a newly introduced cereal (Apple-Cinnamon Cheerios).

Results: the consumer surplus brought by this one new cereal is sizable, amounting to a few percent of category spending; more generally, ignoring new goods overstates the relevant subcomponent of the CPI by about a quarter. The methodological contribution is to pin the vague question "what is a new good worth" down to a CV area under a demand curve, making it estimable and comparable.

Limitations: the whole set of numbers hinges on the demand functional form and on the extrapolation to that reservation price, and the data are sparsest at the high-price end, exactly where the extrapolation is least grounded. The validity of the instruments also depends on the judgment that price shocks across cities are uncorrelated. These are exactly the ready-made illustration for Section 7's point that welfare measures are sensitive to functional form.
:::

### 5.2 Brynjolfsson, Collis and Eggers (2019)

::: {.case}
Paper: "Using massive online choice experiments to measure changes in well-being," PNAS. The question is direct: digital products with a zero price like search, social, and maps leave almost no trace in GDP, yet they clearly bring enormous welfare, so how should this welfare be measured.

Method: precisely this chapter's equivalent variation. The researchers run incentive-compatible choice experiments on a large sample of users, asking "how much money would you need to go without a product for the next month," and the minimum compensation someone will accept is the money that leaves the user indifferent between "having the product with no compensation" and "not having the product with compensation," which by definition is the EV of losing the product. Summing the individual compensations gives the total consumer welfare of the product.

Data: multiple rounds of online choice experiments covering thousands of users, involving digital products like Facebook, search engines, email, and maps, with some experiments cashing out choices with real payments to guarantee incentive compatibility.

Results: users value these zero-price products highly, with the median compensation demanded for search engines the highest, around seventeen thousand dollars a year, email and maps each on the order of a few thousand dollars a year, and Facebook's median around tens of dollars a month. The core message: measured by EV, the consumer welfare of free digital products is vast, and this welfare completely escapes a spending-based GDP.

Limitations: reported compensation demands carry hypothetical bias and the endowment effect, individual valuations are highly dispersed, and extrapolating the median to the population calls for care; these numbers measure private welfare, which need not equal social welfare once externalities and attention costs are counted in.
:::

Together the two papers show the empirical value of this chapter's machine: CV and EV are not symbols on a blackboard but the only rigorous way to answer the real question "what is a new product or a free service actually worth." They also foreshadow the theme of the subsequent demand estimation chapters, namely, estimate demand first, then use duality to translate demand into welfare.

## 6 Running Case: Settling Lumen's Accounts

Now drive the machine of Sections 3 and 4 onto Lumen's numbers, and the losses that could not be pinned down in Section 1 become clear one by one. Every number can be done by hand with elementary algebra, and the R code below only transcribes the hand computation faithfully and checks it, with every number cited in the text coming from its actual output.

### 6.1 UMP: Demand and Indirect Utility

The first-order conditions of Cobb-Douglas $u = x_1^{2/5} x_2^{3/5}$ give the famous constant-expenditure-share result: the fraction of the budget spent on each good is identically equal to its exponent. So Walrasian demand is

$$x_1(p, w) = \frac{2}{5}\frac{w}{p_1}, \qquad x_2(p, w) = \frac{3}{5}\frac{w}{p_2}$$

Plugging in the initial prices $p^0 = (2, 1)$ and $w = 100$: $x_1 = \tfrac{2}{5}\cdot\tfrac{100}{2} = 20$, $x_2 = \tfrac{3}{5}\cdot 100 = 60$. Indirect utility $v(p, w) = w\,(2/5\,/\,p_1)^{2/5}(3/5\,/\,p_2)^{3/5}$, which comes out to $u^0 = 20^{2/5}\,60^{3/5} = 20\cdot 3^{3/5} \approx 38.664$. After the increase, $p^1 = (4, 1)$: $x_1 = 10$, $x_2 = 60$, $u^1 \approx 29.302$. Note that $x_1$ drops from 20 to 10, but $p_1 x_1$ does not budge from 40 to 40. The unmoving 40 is exactly what the first analyst in Section 1 saw. The constant-share property of Cobb-Douglas makes spending completely unresponsive to price, and makes the illusion that "unchanged spending means no loss" especially stubborn.

```r
a  <- 2/5; w <- 100; p2 <- 1
p1_0 <- 2; p1_1 <- 4
x1 <- function(p1) a * w / p1
u  <- function(p1) x1(p1)^a * ((1 - a) * w / p2)^(1 - a)
c(x1_0 = x1(p1_0), x1_1 = x1(p1_1), u0 = u(p1_0), u1 = u(p1_1))
#>     x1_0     x1_1       u0       u1
#>  20.0000  10.0000  38.6636  29.3016
```

### 6.2 EMP: The Expenditure Function and Hicksian Demand

The expenditure function of Cobb-Douglas is $e(p, u) = u\,(p_1/a)^{a}(p_2/(1-a))^{1-a} = u\,(5p_1/2)^{2/5}(5p_2/3)^{3/5}$. Run one self-check: $e(p^0, u^0)$ should equal exactly $w = 100$ (spend 100 to hold the utility $u^0$ she could already reach), and plugging in does give 100; likewise $e(p^1, u^1) = 100$. These two self-checks are the numerical confirmation of the duality identity $e(p, v(p,w)) = w$. Shephard's lemma gives Hicksian demand $h_1 = \partial e/\partial p_1$, which at the reference point $(p^0, u^0)$ computes to $h_1 = 20$, equal to the Walrasian $x_1 = 20$ at that point, and this is the numerical confirmation of the duality identity $h_\ell(p, u^0) = x_\ell(p, e(p, u^0))$.

```r
e <- function(p1, u) u * (p1 / a)^a * (p2 / (1 - a))^(1 - a)
c(check_e0 = e(p1_0, u(p1_0)), check_e1 = e(p1_1, u(p1_1)))
#> check_e0 check_e1
#>      100      100
```

### 6.3 The Slutsky Decomposition

Do the own-price Slutsky decomposition at the initial point. The total effect $\partial x_1/\partial p_1 = -\tfrac{2}{5}\,w/p_1^2 = -(2/5)(100)/4 = -10$. The income effect $-x_1\,\partial x_1/\partial w = -20\cdot(a/p_1) = -20\cdot(1/5) = -4$. The substitution effect by subtraction is $-10 - (-4) = -6$, which can also be checked by differentiating $\partial h_1/\partial p_1$ directly, with the same result.

$$\underbrace{-10}_{\text{total}} = \underbrace{-6}_{\text{substitution}} + \underbrace{(-4)}_{\text{income}}$$

Both halves are negative, consistent with good 1 being a normal good: the increase both pushes her toward the outside good (substitution) and makes her really poorer and cut back across the board (income), two effects in the same direction, so demand slopes down for sure. The decomposition also tells us that six tenths of this price response is pure substitution and four tenths is shrinking purchasing power, and this ratio reappears when we measure welfare in the next section.

### 6.4 The Three Welfare Numbers

Now compute the loss left hanging in Section 1. The three areas are obtained by integrating the functions of 6.1 and 6.2 directly, or by taking expenditure differences:

$$\Delta\mathrm{CS} = \int_2^4 \frac{40}{p_1}\, dp_1 = 40\ln 2 \approx 27.726$$

$$\mathrm{EV} = e(p^0, u^0) - e(p^0, u^1) = 100 - e(p^0, u^1) \approx 24.214$$

$$\mathrm{CV} = e(p^1, u^0) - e(p^1, u^1) = e(p^1, u^0) - 100 \approx 31.951$$

```r
dCS <- 40 * log(2)
EV  <- w - e(p1_0, u(p1_1))
CV  <- e(p1_1, u(p1_0)) - w
c(EV = EV, dCS = dCS, CV = CV)
#>       EV      dCS       CV
#>  24.2142  27.7259  31.9508
```

The ranking $\mathrm{EV} < \Delta\mathrm{CS} < \mathrm{CV}$ fits the general result of Section 4.1 exactly. The three numbers each answer a different but equally legitimate question. CV says: the increase is a done deal, and to bring Lumen back to as satisfied as before, we must compensate her 31.95 at the new prices. EV says: if the increase had not happened, taking 24.21 from her would push her down to as badly off as after the increase. The consumer surplus of 27.73 is an approximation obtainable by integrating observable demand alone, with no knowledge of utility, and it honestly falls between the two exact answers. The disagreement between the two analysts in Section 1 is now resolved: the first's "zero loss" confused spending with welfare; the second's 27.73 is a good approximation, but neither CV nor EV, and which one to actually report depends on whether the operations team means "compensate to no loss" or "buy out in advance." The near-a-third disagreement is not anyone getting the arithmetic wrong, it is that the three questions were different all along.

### 6.5 Two Companion Sub-Accounts: Producer and Uncertainty

Quick computations of the same logic in two other places, to corroborate Sections 4.2 and 4.3. On the producer side, suppose the platform's content delivery uses two inputs to produce output via $y = z_1^{1/2} z_2^{1/2}$, with factor prices $r = (4, 1)$, and it must deliver $y = 10$ units. The cost function $c(r, y) = 2y\sqrt{r_1 r_2} = 2\cdot 10\cdot\sqrt{4} = 40$, Shephard's lemma gives conditional factor demand $z_1 = \partial c/\partial r_1 = y\sqrt{r_2/r_1} = 5$, $z_2 = y\sqrt{r_1/r_2} = 20$, and substituting back to check, $r_1 z_1 + r_2 z_2 = 4\cdot 5 + 1\cdot 20 = 40$, agrees with the value of the cost function.

On the uncertainty side, suppose a seller has Bernoulli utility $u(x) = \sqrt{x}$, where $x$ is the final wealth in hand. She considers switching to a new pricing scheme, whose outcome is a 50-50 lottery: 144 in the good state, 64 in the bad. Expected wealth $\mathbb{E}[x] = 104$, expected utility $\mathbb{E}[u] = \tfrac{1}{2}\sqrt{64} + \tfrac{1}{2}\sqrt{144} = \tfrac{1}{2}\cdot 8 + \tfrac{1}{2}\cdot 12 = 10$, certainty equivalent $\mathrm{CE} = 10^2 = 100$, risk premium $\pi = \mathbb{E}[x] - \mathrm{CE} = 104 - 100 = 4$. That is, this lottery with an expected value of 104 is, in her eyes, only as good as a certain wealth of 100, and she would rather give up 4 units of expected wealth to dodge this volatility. Check with the Arrow-Pratt approximation: $A(x) = 1/(2x)$, at $\mathbb{E}[x] = 104$ giving $A \approx 0.0048$, the lottery variance $\mathrm{Var}(x) = 1600$, and the approximate risk premium $\tfrac{1}{2}A\,\mathrm{Var} \approx \tfrac{1}{2}\cdot 0.0048\cdot 1600 \approx 3.846$, close to the exact 4, with the gap coming from the fact that this risk is not exactly "small," and the approximation formula was built for small risks.

### 6.6 Three Reverse Self-Checks: Recovery, Aggregation, and Stochastic Dominance

Three more small accounts, bringing the new objects of Sections 3.6, 3.7, and 4.3 down to numbers.

The first is integrability. Since Lumen's Cobb-Douglas demand comes from a rational preference, its Slutsky matrix must be symmetric and negative semidefinite, which is exactly the integrability condition of Section 3.6. Compute and check the whole $2\times 2$ Slutsky matrix at the initial point $(p^0, w) = ((2,1), 100)$:

```r
a <- 2/5; w <- 100; p1 <- 2; p2 <- 1
x1 <- a*w/p1; x2 <- (1-a)*w/p2
S <- matrix(c(-a*w/p1^2 + x1*(a/p1),   0 + x2*(a/p1),
              0 + x1*((1-a)/p2),      -(1-a)*w/p2^2 + x2*((1-a)/p2)),
            nrow = 2, byrow = TRUE)
S
#>      [,1] [,2]
#> [1,]   -6   12
#> [2,]   12  -24
c(sym = S[1,2] - S[2,1], Sp = (S %*% c(p1, p2))[1])   # symmetry & S p = 0
#> sym  Sp
#>   0   0
eigen(S)$values
#> [1]   0 -30
```

The diagonal element $s_{11} = -6$ is exactly the substitution effect computed by hand in Section 6.3, the off-diagonal elements $s_{12} = s_{21} = 12$ being equal means symmetry holds, the two eigenvalues $0$ and $-30$ both being nonpositive means negative semidefiniteness, and $S p = 0$ reflects that the Slutsky matrix always has the price vector as a zero direction (a corollary of compensated demand's homogeneity of degree zero). All three pass, showing this set of demands is indeed integrable and can have its expenditure function recovered, and Section 6.2 assembling $e$ from it was no accident.

The second is Gorman aggregation. Build two consumers, both with a linear expenditure system, with different subsistence demands (one prefers basic features $\gamma^A = (5, 10)$, one is a heavy user $\gamma^B = (15, 0)$), but the same marginal budget share for good 1, $b_1 = 2/5$. With total income $W = 300$, split first as $(100, 200)$, then changed to an even $(150, 150)$, and see whether total demand $X_1$ moves:

```r
b1 <- 2/5
x1i <- function(g1, g2, wi) g1 + b1*(wi - p1*g1 - p2*g2)/p1
c(X1 = x1i(5,10,100) + x1i(15,0,200), X1_redis = x1i(5,10,150) + x1i(15,0,150))
#>       X1 X1_redis
#>       70       70
x1i_b <- function(g1, g2, wi, b) g1 + b*(wi - p1*g1 - p2*g2)/p1     # make B's slope 0.7
c(orig = x1i_b(5,10,100,.4) + x1i_b(15,0,200,.7),
  redis = x1i_b(5,10,150,.4) + x1i_b(15,0,150,.7))
#>  orig redis
#>  95.5  88.0
```

With a shared slope, the two individual demands change from $(21, 49)$ to $(31, 39)$, one up and one down, yet total demand stays put at $70$, because it recognizes only $W$ and is independent of how the money is distributed, and this is the numerical face of the representative consumer's existence. But once B's marginal budget share is nudged to $0.7$, breaking the common slope, redistribution immediately pushes total demand from $95.5$ to $88.0$, and the representative consumer fails accordingly. The warning of Section 3.7 is seen here in full clarity.

The third is stochastic dominance. Using the seller $u(x) = \sqrt{x}$ of Section 4.3, compare two lotteries with the same mean of $104$: a tighter $A$ (50-50 on $\{100, 108\}$), and a mean-preserving spread $B$ that spreads it toward the two ends (50-50 on $\{64, 144\}$). SOSD asserts that any risk-averter should prefer the tighter $A$:

```r
u <- sqrt
c(mean_A = 0.5*100 + 0.5*108, mean_B = 0.5*64 + 0.5*144,
  EU_A = 0.5*u(100) + 0.5*u(108), EU_B = 0.5*u(64) + 0.5*u(144))
#>   mean_A   mean_B     EU_A     EU_B
#> 104.0000 104.0000  10.1962  10.0000
```

Both lotteries have a mean of $104$, but the expected utility $\mathrm{EU}_A = 10.196 > \mathrm{EU}_B = 10.000$, and the concave utility penalizes the more dispersed $B$, matching the SOSD judgment to the digit. It is worth noting that $B$ is exactly the $\{64, 144\}$ lottery of Section 6.5, whose certainty equivalent $100$ and risk premium $4$ were already computed, and here we swap in a tighter lottery of the same mean to compare, highlighting not "whether to take the risk" but "which of two risks is better," with the dominance criterion letting this comparison proceed without knowing the specific scale of $u$.

The accounts of this section are now all settled: one and the same Cobb-Douglas consumer, one and the same duality machine, with the UMP producing demand, the EMP producing expenditure, the envelope theorem producing Slutsky, Shephard, and Hotelling, and expenditure differences producing CV and EV, and the reverse self-checks verifying integrability, aggregation, and stochastic dominance in turn, every step doable by hand and self-checkable, while the two small producer and uncertainty accounts prove this logic can be reused unchanged. The confusion of Section 1 was not resolved by some cleverer formula, but by the discipline of "first ask clearly which welfare quantity you want."

## 7 Failure Modes of Measurement and Robustness

The blackboard duality machine assumes regular preferences, demand from a single rational individual, and a known functional form, and in real welfare analysis these premises loosen at any moment. This section lays out the most common failures and the workable responses.

The path dependence of consumer surplus is the first trap. When there are several simultaneous price changes (say the platform adjusts the prices of two features at once), writing consumer surplus as a line integral over prices generally makes the integral depend on the order in which you let prices change, and adjusting good 1 first or good 2 first gives different numbers. CV and EV do not have this affliction, because they are defined by the difference of two values of the expenditure function, independent of the path. Whether consumer surplus can escape path dependence depends on whether Slutsky symmetry $\partial h_\ell/\partial p_k = \partial h_k/\partial p_\ell$ holds, which is exactly the gift from duality of Section 3.5. The practical lesson: welfare analysis with multiple price changes should compute CV or EV directly from the expenditure function, and should not integrate Walrasian demand for convenience.

Sensitivity to functional form is the second, and also the soft spot of the paper in Section 5.1. Welfare measures want the entire demand curve, especially the high-price end far from the observed price (valuing a new good even requires extrapolating to the reservation price), while the data are often dense only near the current price. Different demand functional forms can fit equally well within the observed interval yet, extrapolated to where the data are sparse, give areas under curves that differ greatly, so the estimate of CV is rather fragile to the choice of functional form. The honest thing to do is to compute the main result under several reasonable functional forms, report a range rather than a point, and make clear which stretch of extrapolation the number is most sensitive to.

Aggregation is the third. This chapter has been about a single consumer throughout, but policy evaluation usually wants total market welfare, and summing individuals is done either by directly adding individual CV and EV, or by pretending a "representative consumer" exists. The latter is legitimate only under very strong conditions: Gorman proved that if and only if all consumers' expenditure functions have the Gorman polar form $e_i(p, u_i) = a_i(p) + b(p)\,u_i$ with a shared slope $b(p)$, does aggregate demand look as if it came from a single rational individual, and only then can welfare be aggregated unambiguously. The substance of this condition is that everyone has the same spending pattern for a marginal dollar of income (marginal propensities to consume independent of the individual), which is hard to satisfy in reality. This explains why the demand estimation later in this series invests great effort in modeling consumer heterogeneity and cannot adopt a representative consumer for convenience, since the legitimacy of welfare aggregation hinges precisely on this.

Quasi-linearity as a special case is worth singling out, because it is the common "exemption zone" for the preceding problems. When the good under study has no income effect, CV, EV, and consumer surplus collapse into one (Section 4.1), path dependence vanishes, and aggregation becomes clean too (quasi-linear preferences automatically satisfy the Gorman form). This is why a great many applied theory models simply assume quasi-linear utility: not because it is true, but because it makes welfare analysis unambiguous. Before using it, be clear-eyed that you are assuming the spending under study takes up only a small slice of the consumer's budget and the income effect is negligible, a premise that fails when studying large expenditures (housing, health care).

Last is the extra fragility that uncertainty brings. Expected utility relies on the independence axiom, and experiments like the Allais paradox repeatedly show that real individuals systematically violate it, and behavioral economics' prospect theory and other alternative models rewrote choice under risk on that basis. This chapter uses expected utility because it is still the default starting point for structural modeling and can consistently connect to the earlier duality framework, but the reader should know it is an assumption that can be rejected by data, not a logical necessity.

Threading these failures together, this chapter's welfare machine is precise and useful, but every one of its numbers carries premises: regular preferences, the right functional form, additive individuals, rationality under risk. Robust welfare analysis is not about how well you have memorized the formulas, but about putting these premises on the table and making clear how your conclusion would move as the premises loosen. This is entirely of a piece with the attitude toward identification assumptions in every later chapter of the series: the credibility of a method lies not in the ingenuity of the technique but in the transparency of the assumptions.

## 8 Further Reading

::: {.readings}
Required reading, in suggested reading order:

- Nolan Miller, Notes on Microeconomic Theory, Chapters 3 to 6. The intuition for the consumer and uncertainty parts of this chapter is mainly modeled on these notes, which cover the UMP, EMP, duality, and Slutsky more patiently than most textbooks, and are well suited for a first read to build geometric intuition.
- Mas-Colell, Whinston and Green, Microeconomic Theory, Chapter 3. The authoritative reference for duality and demand theory, with Section 3.G on integrability (backing preferences out of demand) the complete source for the relevant claims in Sections 3.6 and 7 of this chapter, worth close study.
- Mas-Colell, Whinston and Green, Chapter 5. Producer theory, read against the duality of Chapter 3 to see exactly where the "word-for-word translation" of Section 4.2 does its translating.
- Mas-Colell, Whinston and Green, Chapter 6. Choice under uncertainty, with rigorous statements of the expected utility representation theorem and the Arrow-Pratt measures.

Further reading:

- Hausman (1996). The methodological origin of valuing new goods, the classic application of the CV computation of Section 4.1 of this chapter, with attention to the handling and controversy of the reservation-price extrapolation step.
- Brynjolfsson, Collis and Eggers (2019). The representative work applying EV to zero-price digital products, especially to be read for the IS and digital economics direction.
- Deaton and Muellbauer (1980). The original text of the Almost Ideal Demand System, the source of the demand system used in Section 5.1, and an evergreen tool of empirical demand analysis.
- Diamond and McFadden (1974). The methodological foundation of using the expenditure function for welfare analysis, making clear why CV and EV are superior to consumer surplus.
- Samuelson (1938) and Houthakker (1950). The two founding works of revealed preference, the former proposing WARP and the latter completing the full rationalizability condition with SARP, and reading them together makes clear the lineage of the two axioms in Section 3.6.
- Afriat (1967) and Varian (1982). The classics that made revealed preference into an operational test, with Afriat giving a constructive criterion for the rationalizability of finite observations and Varian turning it into a nonparametric method for empirical demand analysis, the applied outlet for the integrability of Section 3.6.
- Gorman (1961). The original source of the existence condition for the representative consumer, from which the Gorman polar form of Section 3.7 comes, extremely short in length yet weighty in every word.
- Rothschild and Stiglitz (1970). Defining "more risky" precisely with the mean-preserving spread, the methodological basis for the SOSD criterion of Section 4.3.
:::

::: {.apa-refs}
- Afriat, S. N. (1967). The construction of utility functions from expenditure data. *International Economic Review, 8*(1), 67-77. https://doi.org/10.2307/2525382
- Brynjolfsson, E., Collis, A., & Eggers, F. (2019). Using massive online choice experiments to measure changes in well-being. *Proceedings of the National Academy of Sciences, 116*(15), 7250-7255. https://doi.org/10.1073/pnas.1815663116
- Deaton, A., & Muellbauer, J. (1980). An almost ideal demand system. *American Economic Review, 70*(3), 312-326.
- Diamond, P. A., & McFadden, D. L. (1974). Some uses of the expenditure function in public finance. *Journal of Public Economics, 3*(1), 3-21. https://doi.org/10.1016/0047-2727(74)90020-6
- Gorman, W. M. (1961). On a class of preference fields. *Metroeconomica, 13*(2), 53-56. https://doi.org/10.1111/j.1467-999X.1961.tb00819.x
- Hausman, J. A. (1996). Valuation of new goods under perfect and imperfect competition. In T. F. Bresnahan & R. J. Gordon (Eds.), *The economics of new goods* (NBER Studies in Income and Wealth, Vol. 58, pp. 207-248). University of Chicago Press.
- Houthakker, H. S. (1950). Revealed preference and the utility function. *Economica, 17*(66), 159-174. https://doi.org/10.2307/2549382
- Mas-Colell, A., Whinston, M. D., & Green, J. R. (1995). *Microeconomic theory*. Oxford University Press.
- Miller, N. H. (2006). *Notes on microeconomic theory* [Unpublished lecture notes]. Harvard Kennedy School.
- Rothschild, M., & Stiglitz, J. E. (1970). Increasing risk: I. A definition. *Journal of Economic Theory, 2*(3), 225-243. https://doi.org/10.1016/0022-0531(70)90038-4
- Samuelson, P. A. (1938). A note on the pure theory of consumer's behaviour. *Economica, 5*(17), 61-71. https://doi.org/10.2307/2548836
- Varian, H. R. (1982). The nonparametric approach to demand analysis. *Econometrica, 50*(4), 945-973. https://doi.org/10.2307/1912771
:::
