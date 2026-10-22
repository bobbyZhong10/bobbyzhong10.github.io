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

Plover is a fictional content-subscription platform. It raises the price of its premium content from USD 2 to USD 4 an hour and wants to know how much its subscribers lose as a result. The transaction log says zero; integrating the demand curve gives USD 27.73; with more care, the compensating variation (CV) and the equivalent variation (EV) come out at USD 31.95 and USD 24.21. All four numbers could go on a respectable slide. The trouble is that they do not answer the same question.

Welfare is not a column that sits ready-made in the transaction log. To get from prices, quantities, and budgets to "how much did subscribers actually lose," you first have to say how a person chooses under a budget constraint, and then split the response to a price change into a substitution effect and an income effect. That is the job of consumer theory. Utility maximization and expenditure minimization describe one and the same choice from two directions, and duality ties Marshallian demand, Hicksian demand, indirect utility, and the expenditure function together; Roy's identity, Shephard's lemma, and the Slutsky equation are the drivetrain that turns observable demand into welfare. On the producer side, cost minimization and profit maximization are the same mathematics, and choice under uncertainty extends the same language to insurance and risk decisions.

The chapter proceeds as follows. Section 1 sets out the four numbers. Section 2 writes down the two optimization problems and the chapter's target quantities. Section 3, the most detailed in the chapter, builds up preferences, duality, revealed preference, and aggregation piece by piece from numbered assumptions and theorems, each theorem followed by a proof sketch and an "In words" paragraph. Section 4 covers welfare measurement itself and two extensions, to producers and to uncertainty. Section 5 reads three anchoring papers. Section 6 settles the accounts on Plover's four data blocks: W0, a hand-computed block for a single subscriber; W1, a price experiment on ten thousand subscribers; W2, a revealed-preference panel; and W3, a risk block built around an add-on protection plan. Section 7 discusses failure modes. The demand estimation, market power, and merger simulation of Chapters 27 to 29 all borrow this machinery; if the accounting conventions are not settled here, computing more precisely later may only mean reporting more wrong decimal places.

## 1 A Loss You Cannot Pin Down

::: {.case}
To shrink the problem to a size that can be done by hand, take a single subscriber and group her monthly spending into two categories: hours of premium content watched (good 1, unit price $p_1$, initially USD 2 an hour), and a composite outside good that packages everything else (good 2, unit price $p_2 = 1$, which can be read as "what USD 1 buys of everything else"). The monthly budget she allots to the two is $w = 100$ dollars; this is the budget of a two-good problem, not household income. Her preferences can be summarized by the Cobb-Douglas utility $u(x_1, x_2) = x_1^{2/5}\, x_2^{3/5}$: every step can be done by hand, and Section 6.1 shows the price of that convenience.

Now Plover raises the price of premium content from USD 2 to USD 4. The operations team wants to know how much this price increase cost her.
:::

The first analyst leafs through the transaction log. Before the increase she watched 20 hours of premium content a month and spent USD 40; after it she watches only 10 hours, but at double the price she still spends USD 40. Her spending did not move by a cent, so the conclusion is that the increase cost her nothing. The flaw is mistaking spending for welfare: spending stayed flat precisely because she cut her usage in half, and that cut is exactly what she has to endure.

The second analyst integrates the demand curve. Her demand for premium content is $x_1 = (2/5)\,w / p_1 = 40 / p_1$; integrating it as $p_1$ moves from 2 to 4, the loss of consumer surplus is $\int_2^4 (40/p_1)\, dp_1 = 40\ln 2 \approx 27.73$ dollars, close to three tenths of her monthly budget.

There are two more numbers, every bit as proper. Ask "how much would we have to pay her after the increase to leave her as well off as before," and the answer is USD 31.95. Ask "before the increase, how much could we take from her to leave her exactly as badly off as the increase does," and the answer is USD 24.21. The former is the compensating variation, the latter the equivalent variation. Both are logically airtight welfare measures, yet they differ by USD 7.74, about 28% of the USD 27.73 of consumer surplus that sits between them. The question "how much did she lose" does not point to a unique number until we say which compensation we want and at which prices it is valued.

Plover's decision is not about one person, either. It has run a price experiment on ten thousand subscribers, each facing a randomly assigned price between USD 2 and USD 4. How should the average loss over ten thousand people be computed, how far off do we land by treating them as a single "representative subscriber," and can the results be trusted when extrapolated to USD 12, a price the experiment never tried? Sections 3 and 4 first build the machine, labeling every part with the assumptions it needs; Section 6 recomputes all these numbers from scratch.

## 2 The Economic Model and the Target Quantity

The operations team wants a welfare number, but welfare is not a primitive. The primitives are preferences and the choices they generate under a budget constraint, that is, demand. Demand is observable; it is right there in the log. Welfare is not directly observable; no table records "how much the increase hurt her." The chapter is therefore organized around a bridge from observable demand back to unobservable welfare, and duality supplies its piers. The target quantity has two layers: the bottom layer is the demand function and the utility structure behind it, and the top layer is CV and EV.

The notation is as follows: a consumption bundle is $x = (x_1, x_2) \in \mathbb{R}^2_+$, prices are $p = (p_1, p_2) \gg 0$, the budget is $w > 0$, and the utility function is $u(x)$. The consumer picks the most preferred bundle among those she can afford:

$$\max_{x \geq 0}\; u(x) \quad \text{s.t.} \quad p \cdot x \leq w$$

This is the utility maximization problem (UMP). Its solution is Walrasian demand (also called Marshallian demand) $x(p, w)$, and the utility of the optimal bundle is the indirect utility function $v(p, w) = u(x(p, w))$.

The welfare side needs a second engine: the least it costs to hold utility at a given level. When prices change, that least cost changes with them, and the size of the change is the money value of the price change. This is the expenditure minimization problem (EMP):

$$\min_{x \geq 0}\; p \cdot x \quad \text{s.t.} \quad u(x) \geq \bar u$$

Its solution is Hicksian demand (also called compensated demand) $h(p, \bar u)$, and the minimal outlay is the expenditure function $e(p, \bar u) = p \cdot h(p, \bar u)$. Along Hicksian demand, utility is pinned at $\bar u$, the change in real income caused by the price change is compensated away as it happens, and what remains is pure substitution.

With the expenditure function in hand, the welfare measures have precise definitions. Let prices change from $p^0$ to $p^1$, with old utility $u^0 = v(p^0, w)$ and new utility $u^1 = v(p^1, w)$:

$$\mathrm{CV} = e(p^1, u^0) - e(p^1, u^1) = e(p^1, u^0) - w, \qquad \mathrm{EV} = e(p^0, u^0) - e(p^0, u^1) = w - e(p^0, u^1)$$

CV is valued at the new prices and asks how much extra money it takes to compensate her back to the old utility; EV is valued at the old prices and asks how much could be taken from her, absent the increase, to drop her to the new utility. The two differ only in which prices serve as the measuring stick. This chapter signs CV, EV, and the loss of consumer surplus $\Delta\mathrm{CS}$ as losses, positive for a price increase; Mas-Colell, Whinston and Green (1995, hereafter MWG) 3.I, in equations (3.I.1) and (3.I.2), signs them as gains, negative for a price increase, so the two conventions differ by a sign. The papers of Section 5 measure welfare gains, and the text reports their absolute values.

The price experiment extends the target quantity to a population. Write the CV of subscriber $i$ as $\mathrm{CV}_i$. The platform cares about the average loss $\overline{\mathrm{CV}}$, the average of $\mathrm{CV}_i$ over all subscribers; a common shortcut treats aggregate demand as the demand of a single representative consumer and computes $\mathrm{CV}^{\mathrm{rep}}$. Section 3.7 takes up when the two are equal; EV and $\Delta\mathrm{CS}$ are written the same way.

Define the target quantity first, then discuss how to compute it: this discipline runs through the whole book. "How much did she lose" must first be translated into "CV or EV, individual or average, valued at which prices."

## 3 The Core Structure: Preferences, Duality, Revealed Preference, and Aggregation

This section first uses Plover's numbers to show what is to be proved, then states numbered assumptions and theorems, followed by proof sketches in which each step names the assumptions it uses; regularity conditions are only mentioned and cited to their source. After each theorem, an "In words" paragraph answers five questions: under which assumptions the result characterizes which quantity; what it becomes, and in which direction it errs, when the assumptions fail; whose choice or which primitive determines the conclusion; what the data can test and what variation that takes; and what Plover's numbers say. Assumptions carry the chapter's prefixes C. (consumer), P. (producer), and U. (uncertainty); theorems are numbered consecutively within the chapter, with MWG's section, proposition, and page numbers in parentheses, and results that are not in the book cite the original papers. MWG states these results for any number of goods; the chapter's examples use two. Sections 3.1 to 3.5 run forward, from preferences to demand and on to the three signature results of duality; Sections 3.6 and 3.7 run the machine backward (from choices back to preferences) and sideways (from individuals to the market).

### 3.1 Preferences, Utility, and the UMP

Summarizing the subscriber of Section 1 by a Cobb-Douglas utility is itself an assumption. At USD 2 she chooses 20 hours of premium content plus USD 60 of other consumption; at USD 4 she chooses 10 hours plus USD 60. For "she chose the bundle with the highest utility" to mean anything, her preferences must satisfy at least the following four conditions.

::: {.assumption}
**C.1 (rational preferences: completeness and transitivity; MWG Definition 3.B.1, p. 42)** The preference relation $\succsim$ on the consumption set is complete (for any two bundles $x$ and $y$, at least one of $x \succsim y$ and $y \succsim x$ holds) and transitive ($x \succsim y$ and $y \succsim z$ imply $x \succsim z$).

**C.2 (continuity; MWG Definition 3.C.1, p. 46)** Preferences are preserved under limits: if two sequences of bundles satisfy $\succsim$ term by term, so do their limits; equivalently, for each $x$, the upper contour set $\{y : y \succsim x\}$ and the lower contour set $\{y : x \succsim y\}$ are both closed (MWG p. 47).

**C.3 (local nonsatiation; MWG Definition 3.B.3, p. 42)** For any $x$ and any $\varepsilon > 0$, the ball around $x$ of radius $\varepsilon$ contains some $y$ with $y \succ x$.

**C.4 (convexity; MWG Definitions 3.B.4 and 3.B.5, p. 44)** (a) Convex: every upper contour set $\{y : y \succsim x\}$ is convex; (b) strictly convex: if $y \succsim x$, $z \succsim x$, and $y \neq z$, then every point in the interior of the segment joining $y$ and $z$ is strictly preferred to $x$.
:::

C.1 rules out cycles: if she prefers A to B, B to C, and C to A, "the most preferred bundle" is undefined, and the GARP test of Section 3.6 checks exactly this condition on finite data. C.2 rules out lexicographic preferences (compare premium hours first, and other consumption only in a tie); it asks only that indifference curves not break apart suddenly. C.3 says something slightly better is always nearby, and it makes the subscriber spend her whole budget; a user whose time is already saturated fails it for premium content taken alone, but as long as the outside good can absorb leftover money, it still holds overall. C.4(a) says averages beat extremes; a product with a fixed cost of use (use it heavily or not at all) creates a nonconvexity. C.4(b) further rules out straight segments on indifference curves, which makes demand single-valued. Cobb-Douglas satisfies all four on $\mathbb{R}^2_{++}$.

::: {.theorem}
**Theorem 1 (utility representation; MWG Proposition 3.C.1, p. 47)** Under C.1 and C.2, there is a continuous utility function $u$ that represents $\succsim$: $x \succsim y$ if and only if $u(x) \geq u(y)$.
:::

The proof sketch has two steps. MWG proves the case of monotone preferences (pp. 47-48); the proof of the general case is longer and is omitted here. Step one (the completeness part of C.1, C.2, and monotonicity): let $\mathbf{1}$ be the vector whose components all equal 1. For each $x$, continuity gives a point on the diagonal that is indifferent to $x$, and monotonicity ensures that there is only one such point, written $u(x)\,\mathbf{1}$. Step two (C.1, C.2, plus monotonicity): transitivity ensures that $x \succsim y$ if and only if $u(x) \geq u(y)$; C.2 ensures that this "diagonal scale" is a continuous function of $x$.

In words: under C.1 and C.2, "she chooses the bundle with the highest utility" and "she chooses the most preferred bundle" say the same thing, so calculus on $u$ is safe. Utility is only a label for preferences: any strictly increasing transformation represents the same preferences; only the ranking means anything, and welfare can therefore be measured only in money (CV and EV depend only on the expenditure function). When C.2 fails, preferences can be perfectly rational and yet have no real-valued utility representing them, and all the calculus tools that follow break down. Data cannot test "utility exists" itself, only its implications; the GARP test of Section 3.6 is one of them. Plover's numbers: the subscriber's utility at the two prices is 38.664 and 29.302. These two numbers cannot be read as "satisfaction fell by some fraction"; all that can be read is the ranking, and the CV and EV converted into dollars.

Picking the highest-utility point on the budget line is, geometrically, letting the budget line $p_1 x_1 + p_2 x_2 = w$ reach the highest indifference curve it can; at an interior optimum the two are tangent, and the marginal rate of substitution $\mathrm{MRS} = (\partial u/\partial x_1)/(\partial u/\partial x_2)$ equals the price ratio $p_1/p_2$. As Lagrangian first-order conditions, with $\mathcal{L} = u(x) + \lambda(w - p \cdot x)$:

$$\frac{\partial u}{\partial x_\ell} = \lambda\, p_\ell, \qquad \ell = 1, 2$$

The multiplier $\lambda = \partial v(p, w)/\partial w$ is the marginal utility of income, which Section 3.4 needs to derive Roy's identity. In the theorem below, $\alpha > 0$ is an arbitrary constant that scales prices and the budget in proportion.

::: {.theorem}
**Theorem 2 (solution of the UMP, Walrasian demand, and indirect utility; MWG Propositions 3.D.1, p. 50; 3.D.2, pp. 51-52; 3.D.3, p. 56)** Let $p \gg 0$ and $w > 0$, and let $u$ be continuous (C.1, C.2, and Theorem 1). Then (i) the UMP has a solution; (ii) $x(p, w)$ is homogeneous of degree zero: $x(\alpha p, \alpha w) = x(p, w)$; (iii) under C.3, Walras' law holds: every bundle in $x(p, w)$ satisfies $p \cdot x = w$; (iv) under C.4(a), $x(p, w)$ is a convex set, and under C.4(b) it contains a single element; (v) under C.3, $v(p, w)$ is homogeneous of degree zero, strictly increasing in $w$, nonincreasing in each $p_\ell$, quasiconvex, and continuous.
:::

The proof sketch has four steps. Step one (C.1, C.2, $p \gg 0$): the budget set is bounded and closed, hence compact, and a continuous function attains its maximum on a compact set, which is (i). Step two (no assumption needed): multiplying prices and the budget by $\alpha$ leaves the budget set unchanged, so the optimal bundle does not change either, which is (ii); the homogeneity of degree zero of $v$ follows in the same way. Step three (C.3): if some optimal bundle left part of the budget unspent, there would be an affordable and strictly better bundle nearby, a contradiction, which is (iii); this also implies that $v$ is strictly increasing in $w$. Step four (C.4): a convex combination of two optimal bundles is still affordable, and convexity makes it no worse than either, so the set of optima is convex; strict convexity makes it strictly better, contradicting optimality, so the optimal bundle is unique, which is (iv). For the quasiconvexity of $v$, see the proof that follows MWG Proposition 3.D.3; MWG does not prove continuity (p. 56).

In words: under C.1 to C.4(b), the subscriber's choice is a single-valued function $x(p, w)$ that depends only on relative prices and the real budget, and she spends the whole budget. Homogeneity of degree zero says that nominal units do not matter: keep the books for prices and the budget in cents, and the choice is the same. Walras' law leaves only one of the two goods' demands independent; that is how the USD 40 and the other USD 60 of Section 1 lock each other in place. When C.3 fails, part of the budget may go unspent; when C.4(b) fails, demand may be multivalued, a tiny change in price can make the choice jump among several bundles, and the later theorems that rely on differentiability have to be rewritten. The intuition for the quasiconvexity of indirect utility in prices is that averaging two price vectors cannot leave her better off than the better of the two, because any bundle affordable at the averaged prices is affordable at one or both of the original price vectors. These properties are determined by preferences; homogeneity of degree zero and Walras' law can be tested directly against data, and Section 3.6 takes them as given. Plover's numbers: at USD 2 she chooses $(20, 60)$, at USD 4 she chooses $(10, 60)$, and both times she spends exactly USD 100.

### 3.2 The EMP: The Expenditure Function and Compensated Demand

After the increase, what is the least it costs to bring her back to the satisfaction she had before it? To hold $u^0$ at a unit price of USD 4, the cheapest combination is 13.20 hours of premium content plus USD 79.17 of other consumption, USD 131.95 in all; the USD 31.95 beyond her USD 100 is exactly the CV of Section 1. This "least outlay," as a function of prices and the target utility, is the expenditure function. Geometrically, fix the indifference curve $u(x) = \bar u$ and shift the iso-expenditure line $p \cdot x = k$ outward from the origin until it just touches the curve; the point of tangency is $h(p, \bar u)$, and the tangency condition is again $\mathrm{MRS} = p_1/p_2$. The two problems satisfy the same tangency condition at the optimum; one treats the budget as the constraint and utility as the objective, the other reverses the roles, and this is the root of duality.

::: {.theorem}
**Theorem 3 (the expenditure function and Hicksian demand; MWG Propositions 3.E.2, p. 59; 3.E.3, p. 61)** Let $u$ be continuous, let C.3 hold, let $p \gg 0$, and let the target utility $\bar u$ be above $u(0)$ and attainable. Then (i) $e(p, \bar u)$ is homogeneous of degree one in $p$, strictly increasing in $\bar u$, nondecreasing in each $p_\ell$, concave in $p$, and continuous; (ii) $h(p, \bar u)$ is homogeneous of degree zero in $p$, and every bundle in it has utility exactly $\bar u$; under C.4(a), $h(p, \bar u)$ is a convex set, and under C.4(b) it contains a single element.
:::

The proof sketch has three steps. Step one (no assumption needed): the constraint set $\{x : u(x) \geq \bar u\}$ does not depend on prices, so multiplying prices by $\alpha$ leaves the optimal bundle unchanged and multiplies the minimal outlay by $\alpha$, which gives the homogeneity of degree one of $e$ and of degree zero of $h$. Step two ($u$ continuous, $\bar u > u(0)$): if the optimal bundle had utility strictly above $\bar u$, scaling it down slightly would, by the continuity of $u$, still satisfy the constraint and cost less, a contradiction, so there is no surplus utility; that $e$ is strictly increasing in $\bar u$ follows from the same argument. Step three (no assumption needed): take the optimal bundle at the average of two price vectors $p'$ and $p''$. It satisfies the constraint at $p'$ and at $p''$ as well, so the minimal outlay at each price vector is no more than the outlay on this bundle; averaging the two inequalities gives the midpoint concavity of $e(\cdot, \bar u)$, and any convex combination works the same way. Convexity and uniqueness are as in step four of Theorem 2.

In words: under C.3, the expenditure function converts "holding a given level of satisfaction" into money, and it is the only ingredient of CV and EV. Its concavity in prices means that when a price rises, the consumer substitutes toward goods whose prices rose less, so the actual rise in expenditure is smaller than the linear extrapolation "with the bundle held fixed"; concavity is "substitution saves you money." The compensated law of demand (MWG Proposition 3.E.4, p. 62) says that compensated demand is nonincreasing in its own price, without exception; Walrasian demand, by contrast, can slope up because of a strong income effect (a Giffen good), and Section 3.5 writes this sentence as an equation. When preferences are not convex, $e$ is still concave in prices, but it recovers only the convex hull of the preferences (MWG p. 78). The expenditure function is defined over unobservable utility, so data cannot test it directly; it can only be recovered from demand through duality (Section 3.6), and its shape is ultimately determined by the demand the analyst specifies. Plover's numbers: once compensated, she cuts premium content from 20 hours to 13.20 hours and raises other consumption from USD 60 to USD 79.17, and this move is pure substitution.

### 3.3 The Four Duality Identities

Plover's subscriber comes with two self-checks. The expenditure function evaluated at the utility she actually reaches returns exactly her budget: $e(p^0, u^0)$ and $e(p^1, u^1)$ are both USD 100. Compensated demand at the initial point is $h_1(p^0, u^0) = 20$, the same as Walrasian demand. This is no coincidence.

::: {.theorem}
**Theorem 4 (duality of the UMP and the EMP; MWG Proposition 3.E.1, p. 58; equations (3.E.1) and (3.E.4), pp. 60 and 62)** Let $u$ be continuous, let C.3 hold, and let $p \gg 0$. (i) If $x^\ast$ is an optimal bundle of the UMP with budget $w > 0$, it is also an optimal bundle of the EMP with target utility $u(x^\ast)$, and the minimal outlay is exactly $w$; (ii) if $x^\ast$ is an optimal bundle of the EMP with target utility $\bar u > u(0)$, it is also an optimal bundle of the UMP with budget $p \cdot x^\ast$, and the maximal utility is exactly $\bar u$. Hence

$$e(p, v(p, w)) = w, \qquad v(p, e(p, \bar u)) = \bar u$$

$$x_\ell(p, w) = h_\ell(p, v(p, w)), \qquad h_\ell(p, \bar u) = x_\ell(p, e(p, \bar u))$$

When demand is not single-valued, the second row is an equality between sets.
:::

The proof sketch has two steps; the first uses C.3 and the second the continuity of $u$. Step one (C.3): if the UMP's optimal bundle $x^\ast$ were not optimal in the EMP, there would be a bundle with utility no lower than $u(x^\ast)$ that costs less, call it $x'$; local nonsatiation puts near $x'$ a bundle that is still affordable and strictly better, contradicting the optimality of $x^\ast$ in the UMP. If the minimal outlay were below $w$, the same argument would apply once more. Step two ($u$ continuous, $\bar u > u(0)$): $\bar u > u(0)$ ensures that the minimal outlay is positive; if some bundle were better at the budget $p \cdot x^\ast$, scaling it down slightly would, by the continuity of $u$, leave it better yet cheaper, contradicting the optimality of $x^\ast$ in the EMP. The four identities restate (i) and (ii) in terms of value functions and optimal solutions.

In words: under C.3, the UMP and the EMP describe one and the same tangency. At fixed prices $v$ and $e$ are inverses of each other, and Walrasian and Hicksian demand coincide at the point where "the budget exactly equals the least outlay that holds that utility." When C.3 fails, an optimal bundle of the UMP need not be optimal in the EMP, the least outlay may fall below the budget, the chain of duality breaks, and inferring welfare from demand loses its footing. The identities are determined by preferences alone, and data have no way to test duality itself; its use is to carry the properties of one object over to another, and the signature results of the next two subsections are all products of this carrying. Plover's two self-checks are the numerical form of the first identity (at both price vectors) and of the last two identities (at the initial point).

::: {.intuition}
Why is duality worth all the trouble? Because it swaps quantities that are hard to find for quantities that are easy to find. Welfare wants a difference of $e$, but $e$ is defined over unobservable utility, whereas demand $x$ is observable. The identity $h_\ell = x_\ell(p, e(p, \bar u))$ is exactly the channel that lets us assemble compensated demand from observed Walrasian demand, and from there the expenditure function and welfare. Later, in demand estimation, what is estimated is $x$ and what is wanted is welfare, and the procedure in between is the duality here.
:::

### 3.4 Roy's Identity and Shephard's Lemma

As the unit price creeps up, how fast does the least outlay that holds $u^0$ rise? To a first approximation, the answer is the number of premium hours she buys now: the money saved by readjusting the bundle is a second-order quantity at the optimum. This is the envelope theorem: standing at an optimum, computing the effect of a small change in a parameter on the optimal value requires only the direct effect, with no need to track the reoptimization, because the first-order gain from reoptimizing is zero.

::: {.theorem}
**Theorem 5 (Shephard's lemma and Roy's identity; MWG Propositions 3.G.1, p. 68; 3.G.4, p. 73)** Let $u$ be continuous, and let C.3 and C.4(b) hold. (i) For all $p \gg 0$ and $\bar u$, Hicksian demand is the gradient of the expenditure function with respect to prices:

$$\frac{\partial e(p, \bar u)}{\partial p_\ell} = h_\ell(p, \bar u)$$

(ii) If $v$ is differentiable at $(p, w) \gg 0$ and $\partial v/\partial w \neq 0$, then

$$x_\ell(p, w) = -\frac{\partial v(p, w)/\partial p_\ell}{\partial v(p, w)/\partial w}$$

Both conclusions are local; they hold where the optimal solution is locally unique (MWG footnote 15).
:::

The proof sketch has three steps, following MWG's third proof (the envelope theorem, p. 69) and assuming that the optimal solution is differentiable in prices. Step one (C.4(b)): strict convexity makes $h$ and $x$ single-valued, so the derivative of the value function with respect to a parameter is meaningful. Step two (no assumption needed): prices appear only in the objective $p \cdot x$, not in the constraint $u(x) \geq \bar u$, so the derivative of the value function with respect to $p_\ell$ equals the direct derivative of the objective with respect to $p_\ell$, namely $x_\ell$, which at the optimum is $h_\ell$; this is (i). Step three (C.3, C.4(b)): apply the envelope theorem to the UMP, where prices appear in the constraint: $\partial v/\partial p_\ell = -\lambda\, x_\ell$, while $\partial v/\partial w = \lambda$; dividing one by the other cancels the unobservable $\lambda$, which is (ii).

In words: the slopes of the two value functions with respect to prices are demands. Shephard's lemma turns compensated demand into a single derivative of the expenditure function, and Roy's identity translates indirect utility into Walrasian demand in one step. Demand obtained by first positing an $e$ that satisfies Theorem 3 (i) (or a $v$ that satisfies Theorem 2 (v)) and then differentiating automatically satisfies the restrictions of rationality; if the posited $e$ is not concave, the substitution matrix is no longer negative semidefinite. When strict convexity fails (indifference curves with flat segments), $e$ has kinks at some prices and demand jumps with price. What sets the magnitude: the first-order effect of a price increase on outlay is valued at her current purchases, so "the extra money paid, computed on the old bundle" is only a first-order approximation to CV, and the second-order part is determined by substitution (MWG equation (3.I.9), p. 90). On the data side, the analyst can derive demand from a specified utility and then fit it, and Chapter 27, which derives logit choice probabilities from random utility, takes the same road. Plover's numbers: at $(p^0, u^0)$, $\partial e/\partial p_1 = h_1 = 20$, the same as Walrasian demand $x_1(p^0, w) = 20$.

### 3.5 The Slutsky Equation: Splitting the Slope of Demand in Two

The figure below splits the discrete increase from USD 2 to USD 4 into two legs. Top up her budget at the new prices to USD 131.95, and she returns to her original indifference curve, choosing $(13.20, 79.17)$: the move from $(20, 60)$ to this point is substitution, which cuts premium viewing by 6.80 hours. Take the compensation away, and she lands at $(10, 60)$; this leg is the income effect, which cuts premium viewing by another 3.20 hours. The Slutsky equation is the differential form of this decomposition.

![The decomposition for one Plover subscriber: as the unit price rises from USD 2 to USD 4, her choice moves from (20, 60) to (10, 60); topping up the budget to USD 131.95 at the new prices returns her to the original indifference curve, where the compensated choice is (13.20, 79.17). The first leg is substitution, the second the income effect.](assets/fig/fig_01_budget.svg)

::: {.theorem}
**Theorem 6 (the Slutsky equation and the substitution matrix; MWG Propositions 3.G.3, p. 71; 3.G.2, p. 69; the observable Slutsky matrix, p. 72)** Let $u$ be continuous, and let C.3 and C.4(b) hold. (i) For all $(p, w)$ and $\bar u = v(p, w)$,

$$\underbrace{\frac{\partial x_\ell(p, w)}{\partial p_k}}_{\text{total}} = \underbrace{\frac{\partial h_\ell(p, \bar u)}{\partial p_k}}_{\text{substitution}}\; \underbrace{-\; x_k(p, w)\,\frac{\partial x_\ell(p, w)}{\partial w}}_{\text{income}}$$

(ii) If $h(\cdot, \bar u)$ is continuously differentiable at $p$, the substitution matrix $[\partial h_\ell/\partial p_k]$ equals the Hessian (the matrix of second derivatives) of $e$ with respect to prices; it is symmetric and negative semidefinite, and its product with the price vector is zero. By (i), the Slutsky matrix computed from observable demand, $S(p, w) = [\partial x_\ell/\partial p_k + x_k\,\partial x_\ell/\partial w]$, equals it, and is therefore also symmetric and negative semidefinite, with $S(p, w)\,p = 0$.
:::

The proof sketch has three steps. Step one (Theorem 4, using the continuity of $u$): differentiate both sides of the identity $h_\ell(p, \bar u) = x_\ell(p, e(p, \bar u))$ with respect to $p_k$, noting that $e$ on the right also contains $p_k$:

$$\frac{\partial h_\ell}{\partial p_k} = \frac{\partial x_\ell}{\partial p_k} + \frac{\partial x_\ell}{\partial w}\cdot\frac{\partial e}{\partial p_k}$$

Step two (Theorem 5 (i), using C.4(b)): use Shephard's lemma to replace $\partial e/\partial p_k$ with $h_k$, which at the corresponding point equals $x_k$; rearranging gives (i). Step three (Theorems 3 and 5): the substitution matrix is the Hessian of $e$; the symmetry of mixed partial derivatives gives symmetry, and the concavity of $e$ in prices gives negative semidefiniteness. Since $e$ is homogeneous of degree one in prices, its gradient is homogeneous of degree zero, and by Euler's theorem the Hessian times $p$ is zero.

In words: demand estimation delivers the total effect on the left directly, and that effect mixes two things of different natures: substitution caused by the change in relative prices, and the income effect caused by shrinking purchasing power. Under C.3 and C.4(b) the substitution term is always nonpositive, so whether Walrasian demand slopes down depends entirely on the income term. For a normal good ($\partial x_\ell/\partial w > 0$) the two terms have the same sign and demand must slope down; for an inferior good the income term is positive, and when it is strong enough to overwhelm the substitution term, the good is a Giffen good. Symmetry cannot be seen in Walrasian demand; it is a gift of duality, and the strongest testable restriction that rationality places on data (Section 3.6). What sets the magnitude: the substitution term depends only on the curvature of the indifference curve, while the income term is the quantity purchased times the income slope, so the same price increase produces a larger income effect when it falls on a heavy user. Testing symmetry and negative semidefiniteness takes variation in two directions, price and income. Plover's experiment randomizes only prices, so the income slope must come from differences in budgets across subscribers (Section 6.4), and what is measured is average demand, which, when preferences are heterogeneous, may fail to be negative semidefinite even if every individual is rational (MWG 4.C). Plover's numbers: at $p^0$ the total effect $-10$ equals the substitution effect $-6$ plus the income effect $-4$, six tenths pure substitution and four tenths shrinking purchasing power; in the discrete increase the two legs are 6.80 and 3.20 hours.

### 3.6 Revealed Preference and Integrability: Backing Preferences Out of Choices

Everything up to here has run forward: preferences first, then solve for demand. The empirical researcher holds the other end: observed choices, from which she wants preferences and the welfare hidden behind them. Start with the version for continuous demand: the Slutsky matrix of W0 in Section 6.3 is symmetric with nonpositive eigenvalues; must this demand system come from some preference?

::: {.theorem}
**Theorem 7 (integrability; MWG 3.H text, pp. 75-80, no numbered proposition; for two goods, equation (3.H.1), p. 78)** Let demand $x(p, w)$ be single-valued, continuously differentiable, and homogeneous of degree zero, and let it satisfy Walras' law. Then there is a rational preference whose Walrasian demand is $x(p, w)$ if and only if the Slutsky matrix $S(p, w)$ is symmetric and negative semidefinite everywhere (for the technical conditions, see MWG footnote 19). With two goods, symmetry holds automatically, and only negative semidefiniteness needs checking.
:::

The proof sketch has three steps; MWG itself describes its treatment in this section as not fully rigorous. Necessity (Theorem 6, using C.3 and C.4(b)): if demand comes from such a preference, $S$ is the substitution matrix, symmetric and negative semidefinite. Sufficiency, step one (symmetry, plus homogeneity of degree zero and Walras' law): Shephard's lemma and duality give a system of partial differential equations with $e$ as the unknown function, $\partial e/\partial p_\ell = x_\ell(p, e)$, and symmetry is exactly the consistency condition for it to have a solution (with many goods by the Frobenius theorem; with two goods it is an ordinary differential equation, equation (3.H.1)). Sufficiency, step two (negative semidefiniteness): the solved $e$ is concave in prices, and MWG Proposition 3.H.1 (p. 77) builds from it the upper contour sets $\{x : p \cdot x \geq e(p, \bar u) \text{ for all } p \gg 0\}$; the preference these sets define has $e$ as its expenditure function.

In words: under the theorem's conditions, symmetry and negative semidefiniteness of the Slutsky matrix are the entire testable content of "consumer rationality" in demand data, and they are the basis of recoverability: if demand is estimated well enough, the expenditure function can be recovered, up to an irrelevant utility scale, and CV and EV depend only on $e$. When the conditions fail, no preference generates the demand: if $S$ is not symmetric, the computed "CV" depends on the path of integration; if $S$ has a positive eigenvalue, the solved $e$ is not concave. Neither result is anyone's welfare. The conclusion is determined by the shape of the demand function alone, and testing it requires estimating slopes in both directions, price and income; Hausman and Newey (1995), in Section 5.3 of Chapter 14, test restrictions of this kind on nonparametrically estimated gasoline demand. The Slutsky matrix of Plover's Cobb-Douglas demand is symmetric and negative semidefinite everywhere, with eigenvalues 0 and $-30$ at the initial prices (Section 6.3); the expenditure function written down in Section 6.1 is exactly what this demand recovers.

Real data are finite. In Plover's revealed-preference panel W2, each of 300 users faces 8 random budgets. Each budget is worth USD 60 in points, and the hourly prices of premium and standard content are each drawn at random between USD 1 and USD 3, so the budget lines cross often; each user has only 8 observations, and there is no derivative to take. Write the user's $t$th observation as prices $p^t$ and the chosen bundle $x^t$, $t = 1, \dots, \mathcal{T}$. If $p^t \cdot x^t \geq p^t \cdot x^s$ (when $x^t$ was chosen, $x^s$ was affordable too), $x^t$ is said to be directly revealed preferred to $x^s$, written $x^t \mathrel{R^{\mathrm{d}}} x^s$; when the inequality is strict, write $x^t \mathrel{P^{\mathrm{d}}} x^s$; $R$ is the transitive closure of $R^{\mathrm{d}}$. The generalized axiom of revealed preference (GARP) requires that whenever $x^t \mathrel{R} x^s$, $x^s \mathrel{P^{\mathrm{d}}} x^t$ does not hold. Its predecessors are the WARP (weak axiom of revealed preference) of Samuelson (1938) and the SARP (strong axiom of revealed preference) of Houthakker (1950); MWG Proposition 3.J.1 (p. 91) says that a single-valued demand function defined on all budgets and satisfying SARP can be rationalized by a rational preference, which is the version with infinitely many observations. The version for finitely many observations follows.

::: {.theorem}
**Theorem 8 (Afriat's theorem; Afriat 1967; Varian 1982)** For finitely many observations $(p^t, x^t)$ with $p^t \gg 0$, $t = 1, \dots, \mathcal{T}$, the following two statements are equivalent: (i) some locally nonsatiated utility function rationalizes the data, that is, for each $t$, $x^t$ is the bundle with the highest utility on the budget $\{x : p^t \cdot x \leq p^t \cdot x^t\}$; (ii) the data satisfy GARP. Moreover, when (ii) holds, the data can be rationalized by a utility that is continuous, concave, and monotone.
:::

The proof sketch has three steps; for the details of the constructive direction, see Varian (1982). Step one ((i) implies (ii), using C.3): when $x^t \mathrel{R^{\mathrm{d}}} x^s$, $x^s$ was affordable when $x^t$ was chosen, so $u(x^t) \geq u(x^s)$; when $x^s \mathrel{P^{\mathrm{d}}} x^t$, $x^t$ was strictly affordable, and local nonsatiation puts a better and still affordable bundle near it, so $u(x^s) > u(x^t)$. Passing along a chain of $R$, a cycle with a strict link would imply $u(x^t) > u(x^t)$. Step two ((ii) implies that the inequalities have a solution, the hardest step of the theorem): GARP guarantees numbers $u^t$ and positive numbers $\lambda^t$ satisfying the Afriat inequalities $u^s \leq u^t + \lambda^t\, p^t \cdot (x^s - x^t)$. Step three ($\lambda^t > 0$, $p^t \gg 0$): let $u(x) = \min_t \{u^t + \lambda^t\, p^t \cdot (x - x^t)\}$, the lower envelope of finitely many affine functions; it is continuous, concave, and monotone, and it rationalizes the data.

In words: on finite data, "rationality" has one and only one testable implication, namely that the choices contain no revealed-preference cycle with a strict link. When GARP passes, the data can be rationalized by some well-behaved utility, but the theorem says only that one exists, not which one: infinitely many utilities satisfy the Afriat inequalities, and they imply different indifference curves away from the observations, so welfare numbers still need a functional form or bounds (Section 4). Failing the test does not mean the preferences are absurd, either: measurement error and changes in taste both produce violations, so in practice one looks at how severe the violations are. The power of the test is set by the budget sets the researcher designs: when budget lines rarely cross, almost any choices pass. Plover's numbers: all 210 users who exactly optimize a CES (constant elasticity of substitution) utility pass, while only 67.8% of the 90 noisy users do.

The severity of violations is measured by the critical cost efficiency index (CCEI, Afriat's efficiency index; Afriat 1972; Varian 1990). Relax the direct relation to $\varsigma\, p^t \cdot x^t \geq p^t \cdot x^s$ (and the strict relation likewise), with $\varsigma \in [0, 1]$; the supremum of the $\varsigma$ at which the relaxed data still satisfy GARP is the CCEI. $\varsigma = 1$ is GARP itself, and $1 - \mathrm{CCEI}$ can be read as the fraction of the budget the user "wastes." The power of the test is gauged by the benchmark of Bronars (1987): let a group of uniformly random choosers choose on the same budgets, and see how many of them pass. On W2, 28.4% of 3,000 random choosers pass GARP, so the power of the test is 71.6%.

![The CCEI distribution in the W2 panel: all 210 CES optimizers have a CCEI of 1; the GARP pass rate of the 90 noisy users is 67.8%; the pass rate of 3,000 uniformly random choosers on the same budgets is 28.4%.](assets/fig/fig_01_ccei.svg)

### 3.7 From Individuals to the Market: The Gorman Form and the Representative Consumer

Demand estimation and policy evaluation almost always play out at the market level, with aggregated market demand $X_\ell(p, w_1, \dots, w_I) = \sum_i x_{\ell i}(p, w_i)$ in hand, where $I$ is the number of consumers. In general it depends not only on total wealth $W = \sum_i w_i$ but also on how the money is divided among people. Take two consumers with linear-expenditure-system preferences: their subsistence quantities differ, $(5, 10)$ for A and $(15, 0)$ for B, and their marginal budget share of premium content is the same, 2/5. At prices $p^0 = (2, 1)$, when a total budget of USD 300 is reallocated from $(100, 200)$ to $(150, 150)$, A's premium viewing rises from 21 hours to 31 and B's falls from 49 hours to 39, and total demand stays at 70 hours. Change B's share to 0.7, and the same redistribution moves total demand from 95.5 down to 88.0.

::: {.theorem}
**Theorem 9 (Gorman polar form; MWG Proposition 4.B.1, pp. 107-108; Example 4.D.2, p. 119)** Let each of $I$ consumers satisfy C.1 to C.4(b). (i) At every price vector, everyone's wealth expansion paths are straight lines parallel to one another (equivalently, for any distribution of wealth, aggregate demand depends only on total wealth $W$) if and only if each person's indirect utility can be written in the Gorman polar form

$$v_i(p, w_i) = a_i(p) + b(p)\, w_i$$

where $b(p)$ is the same for everyone. (ii) In this form, $\sum_i a_i(p) + b(p)\,W$ is the indirect utility of a normative representative consumer: it generates aggregate demand, and this representative consumer's preferences do not depend on which social welfare function is chosen.
:::

The proof sketch has two steps and proves only sufficiency: MWG leaves it as Exercise 4.B.1, does not prove necessity, and stresses that the result is local. Step one (Roy's identity of Theorem 5, using C.3 and C.4(b)): apply Roy's identity to $v_i$,

$$x_{\ell i}(p, w_i) = -\frac{\partial v_i/\partial p_\ell}{\partial v_i/\partial w_i} = -\frac{\partial a_i/\partial p_\ell}{b(p)} - \frac{\partial b/\partial p_\ell}{b(p)}\, w_i$$

The slope of individual demand in $w_i$ contains only the common $b(p)$. Step two (algebra): in the aggregate, the person-specific intercepts merge into a single term, and the part that varies with wealth recognizes only the total,

$$X_\ell(p, W) = -\frac{\sum_i \partial a_i/\partial p_\ell}{b(p)} - \frac{\partial b/\partial p_\ell}{b(p)}\, W$$

which is exactly the demand that Roy's identity gives for the indirect utility $\sum_i a_i(p) + b(p)\,W$. For the argument for (ii), see MWG Example 4.D.2.

In words: the Gorman form characterizes "aggregate demand recognizes only total wealth," that is, everyone spends an increment of income in the same way: the marginal budget shares are the same for everyone, while the intercepts may differ. Identical homothetic preferences, preferences quasilinear in a common numeraire, and linear expenditure systems with different subsistence quantities but the same marginal budget shares all belong to it. The existence of a representative consumer as such does not require it: as long as wealth is redistributed optimally according to some social welfare function, there is a positive representative consumer (MWG Proposition 4.D.1, pp. 116-117). His preferences, however, change with the welfare function, and using him for welfare comparisons tacitly assumes that wealth is always redistributed optimally according to that welfare function (MWG Definition 4.D.3, p. 118); the Gorman form makes him independent of the welfare function. The data needed to test it are discussed below. Plover's numbers: W1g, where everyone has the same shares, satisfies it, and $\mathrm{CV}^{\mathrm{rep}}$ falls only 0.017% short of $\overline{\mathrm{CV}}$; in W1, where the shares are heterogeneous, it falls 1.15% short.

::: {.intuition}
When Gorman fails, $\mathrm{CV}^{\mathrm{rep}}$ is in general not equal to $\overline{\mathrm{CV}}$. Expanding $\mathrm{CV}_i$ to second order in the price change, the difference between the two is approximately $\tfrac{1}{2}(\Delta p_1)^2\,\mathrm{Cov}\big(x_{1i}, \partial x_{1i}/\partial w\big)$: if heavy users also have steep income slopes, the representative consumer understates the average loss. The approximation serves only to judge the direction; when the price doubles, it overstates the size of the gap. The gap comes from the correlation between heterogeneous marginal budget shares and usage, and it does not require income to be correlated with the shares: in W1, the marginal budget share of subscriber $i$ on premium content, $\beta_i$ (Section 6.4; unrelated to the $\beta$ of later chapters), is independent of the budget, yet subscribers with a large $\beta_i$ also watch more premium content.
:::

Testing this requires the same person to face different prices or budgets. In a cross-sectional experiment with one price per person, a demand regression sees only average demand at each (price, budget) point and cannot reveal the heterogeneity in $\beta_i$. Hausman and Newey (2016) show that average exact consumer surplus is in general not point identified from cross-sectional data, but that it can be bounded using bounds on income effects; the general framework of bounds and identified sets is in Chapter 23.

## 4 Welfare Measures, and Two Extensions of Duality

Section 3 built the machine, and this section puts it to work on three jobs: measuring consumer welfare, characterizing the producer, and handling uncertainty.

### 4.1 CV, EV, and Consumer Surplus

The figure below is the source of the three numbers of Section 1: the areas to the left of three demand curves between USD 2 and USD 4. The CV curve is pinned to the old utility $u^0$ and lies farthest out; the EV curve is pinned to the new utility $u^1$ and lies farthest in; Walrasian demand crosses the former at USD 2 and the latter at USD 4, and lies in between.

![Three premium demand curves for one Plover subscriber: the two Hicksian demands at $u^0$ and $u^1$ bracket the Walrasian demand; as the unit price rises from USD 2 to USD 4, the areas to the left of the three curves are, in order, CV, $\Delta$CS, and EV, that is, USD 31.95, 27.73, and 24.21.](assets/fig/fig_01_welfare_areas.svg)

::: {.theorem}
**Theorem 10 (area representation of CV and EV, and the ordering of the three; MWG equations (3.I.3) and (3.I.4), p. 83; for the ordering, Figure 3.I.3 and the text on p. 89, no numbered proposition)** Let $u$ be continuous, let C.3 and C.4(b) hold, and let only $p_1$ change, from $p_1^0$ to $p_1^1$. (i) With losses signed positive, as in this chapter,

$$\mathrm{CV} = \int_{p_1^0}^{p_1^1} h_1(p_1, p_{-1}, u^0)\, dp_1, \qquad \mathrm{EV} = \int_{p_1^0}^{p_1^1} h_1(p_1, p_{-1}, u^1)\, dp_1$$

while $\Delta\mathrm{CS} = \int_{p_1^0}^{p_1^1} x_1(p_1, p_{-1}, w)\, dp_1$. (ii) If good 1 is a normal good, $\mathrm{EV} \leq \Delta\mathrm{CS} \leq \mathrm{CV}$, for price increases and decreases alike; for an inferior good the ordering is reversed. (iii) If good 1 has no income effect, for instance when preferences are quasilinear in good 2, $u = \phi(x_1) + x_2$, the three are equal.
:::

The proof sketch has three steps. Step one (Theorem 5 (i), using C.4(b)): $\partial e/\partial p_1 = h_1$; integrating over $p_1$ turns the CV and EV of Section 2 into areas under two Hicksian curves. Step two (Theorems 3 and 4, using the continuity of $u$ and C.3): by $h_1(p, u) = x_1(p, e(p, u))$, and taking a price increase as the example, $e(p_1, u^1) \leq w \leq e(p_1, u^0)$ between the two prices; for a normal good demand rises with income, so $h_1(\cdot, u^1) \leq x_1(\cdot, w) \leq h_1(\cdot, u^0)$, and integrating gives (ii). Step three (the hypothesis of (iii)): without an income effect, $x_1$ does not contain $w$, the first-order condition $\phi'(x_1) = p_1/p_2$ pins down $x_1$ by itself, and the three curves coincide, which is (iii).

In words: the three measures are areas under three demand curves over the same price interval, and they differ only in the utility at which each curve is pinned. Consumer surplus needs only observable Walrasian demand; the price of that convenience is that, when income effects are not negligible, it is neither CV nor EV: for a price increase on a normal good, it understates CV and overstates EV. Which one to report is the questioner's choice: to "compensate after the increase so that she loses nothing," report CV; to "buy her out in advance," report EV. Whether the three coincide is determined by preferences: the linear term of quasilinearity must sit on the outside good, whose price does not change, for good 1, whose price does, to be free of income effects, and this is exactly what undergraduate textbooks assume when they use consumer surplus with confidence. What data can test is the size of the income effect, which requires variation in budgets or income; when Section 4.4 of Chapter 29 computes consumer surplus with the log-sum formula, it implicitly assumes that the marginal utility of income (the price coefficient) does not vary with income. Plover's numbers: the three areas are USD 24.21, 27.73, and 31.95, in exactly the order the theorem predicts. Every one of the ten thousand subscribers satisfies $\mathrm{EV}_i \leq \Delta\mathrm{CS}_i \leq \mathrm{CV}_i$; the averages are USD 22.47, 25.55, and 29.23 in that order, and the gap between CV and $\Delta\mathrm{CS}$ is 12.6% of CV.

How far is consumer surplus from CV? The gap comes from the income effect, whose size is set by the income elasticity and by the size of the loss relative to the budget.

::: {.theorem}
**Theorem 11 (Willig approximation; Willig 1976, first-order derivation in this chapter)** Let the conditions of Theorem 10 hold, and let the income elasticity of good 1, $\mathcal{E}_{1w} = (\partial x_1/\partial w)(w/x_1)$, be approximately constant along the price path. Then, with losses signed positive and to first order in $\Delta\mathrm{CS}/w$,

$$\mathrm{CV} \approx \Delta\mathrm{CS} + \frac{\mathcal{E}_{1w}\,\Delta\mathrm{CS}^2}{2w}, \qquad \mathrm{EV} \approx \Delta\mathrm{CS} - \frac{\mathcal{E}_{1w}\,\Delta\mathrm{CS}^2}{2w}$$
:::

The proof sketch has three steps. Step one (Theorems 4 and 5, using the continuity of $u$, C.3, and C.4(b)): view CV as a function of the end price, $\mathrm{CV}(p_1) = e(p_1, u^0) - w$; Shephard's lemma and duality give $d\,\mathrm{CV}(p_1)/dp_1 = x_1\big(p_1, w + \mathrm{CV}(p_1)\big)$. Step two (constant elasticity): expand the right side to first order around the budget $w$, $x_1(p_1, w + \mathrm{CV}) \approx x_1(p_1, w)\,\big(1 + \mathcal{E}_{1w}\,\mathrm{CV}/w\big)$. Step three (to first order in $\Delta\mathrm{CS}/w$): replace $\mathrm{CV}(p_1)$ in the correction term by its zeroth-order approximation $\mathrm{CS}(p_1) = \int_{p_1^0}^{p_1} x_1(\tilde p, w)\, d\tilde p$ and integrate from $p_1^0$ to $p_1^1$ ($d\,\mathrm{CS} = x_1\,dp_1$); the correction term is $(\mathcal{E}_{1w}/w)\int \mathrm{CS}\; d\,\mathrm{CS} = \mathcal{E}_{1w}\,\Delta\mathrm{CS}^2/(2w)$. EV works the same way, starting from $u^1$, with the opposite sign. Willig (1976) turns this approximation into error bounds: when the income elasticity is bounded above and below along the path and $\Delta\mathrm{CS}/w$ is small, $(\mathrm{CV} - \Delta\mathrm{CS})/\Delta\mathrm{CS}$ lies between correction terms computed from the two elasticity bounds; for the exact form of the bounds and their range of validity, see the original paper. MWG does not give these bounds and says only, on p. 89, that the error is small when income effects are small.

In words: the relative gap between consumer surplus and CV is about half the income elasticity times the loss as a share of the budget. When the income effect is small, or when the spending takes up only a small slice of the budget, consumer surplus is a good approximation to CV and EV. The approximation is governed by two primitives, the income elasticity and the size of the loss relative to the budget, so the same price increase produces a larger error for someone whose budget share for the good is larger; with a positive elasticity, $\Delta\mathrm{CS}$ tends to understate CV. The Willig bounds limit the error of $\Delta\mathrm{CS}$ relative to CV, not deadweight loss: deadweight loss is the difference of two large numbers, and an equally small absolute error can be a large fraction of it (the warning from Hausman 1981 that MWG cites on p. 90). Using the correction requires an estimate of the income elasticity, that is, variation in budgets. Plover's numbers: the elasticity in W0 is 1, and the corrected CV is USD 31.57 against the exact 31.95, closing 91.0% of the consumer-surplus error; the mirror formula gives an EV of USD 23.88 against the exact 24.21; writing the coefficient 1/2 as 1 by mistake would give USD 35.41, an overcorrection. For the ten thousand subscribers of W1, correcting each with her own true elasticity closes 92.5% of the error in the average CV, and 83.8% even for the worst-served subscriber.

The Willig correction is only an approximation; integrability tells us that CV can in fact be computed exactly from Walrasian demand.

::: {.theorem}
**Theorem 12 (exact recovery of CV from Walrasian demand; MWG equation (3.H.1), p. 78; Hausman 1981; Vartia 1983)** Suppose there are two goods, $p_2$ is fixed, and demand is single-valued, homogeneous of degree zero, satisfies Walras' law, and is integrable (Theorem 7), with $x_1$ continuously differentiable in income. Then compensated income $e(p_1, u^0)$ is the unique solution of the ordinary differential equation

$$\frac{d\, e(p_1, u^0)}{d p_1} = x_1\big(p_1, e(p_1, u^0)\big), \qquad e(p_1^0, u^0) = w$$

and $\mathrm{CV} = e(p_1^1, u^0) - w$; starting from the initial value $w$ and integrating back from $p_1^1$ to $p_1^0$ gives $e(p_1^0, u^1) = w - \mathrm{EV}$.
:::

The proof sketch has two steps. Step one (Theorems 4 and 5): Shephard's lemma gives $d\,e(p_1, u^0)/dp_1 = h_1(p_1, u^0)$, duality gives $h_1(p_1, u^0) = x_1(p_1, e(p_1, u^0))$, and together they are this equation, with the initial value coming from $e(p^0, u^0) = w$. Step two ($x_1$ continuously differentiable in income): the right side is Lipschitz in the unknown function, so the solution is unique, and knowing Walrasian demand is enough to compute $e(\cdot, u^0)$ step by step along the price path. Hausman (1981) gives closed-form solutions for commonly used demand forms, Vartia (1983) gives a numerical algorithm that works for any demand system, and this chapter uses fixed-step fourth-order Runge-Kutta (Section 6.2).

In words: under the conditions of integrability, CV and EV need no approximation; they are determined entirely by Walrasian demand, and the only error comes from how accurately demand itself is estimated. The hard part of welfare measurement was therefore never that "consumer surplus is imprecise," but whether the demand function can be trusted over the stretch of prices whose area is being computed. When demand is not integrable, the equation can still be solved, but what comes out is nobody's expenditure function. The conclusion is set by the analyst's specification of demand: several demand forms give almost the same CV within the experimental price range, yet may part ways when extrapolated beyond it (Section 6.4). Estimation requires variation in two directions, price and budget (Section 5.3 of Chapter 14). Plover's numbers: solving this equation for W0's demand gives USD 31.9508, matching the CV computed directly from the expenditure function; solving it for the income-term form that the analyst fits to W1 gives an average CV of USD 28.90, against a true value of USD 29.23.

### 4.2 The Producer Side: Cost, Profit, and Markups

Producer theory is almost a word-for-word translation of consumer theory. Suppose Plover's content delivery turns two inputs into output via $y = z_1^{1/2} z_2^{1/2}$, factor prices are $r = (4, 1)$, and it must deliver $y = 10$ units: what is the least cost, and how much of each input does it use?

::: {.assumption}
**P.1 (technology; MWG 5.C, p. 135; the conditions of Propositions 5.C.1 and 5.C.2, pp. 138 and 141)** The production set is nonempty, closed, and satisfies free disposal; with a single output it is described by a production function $f(z)$, with inputs $z = (z_1, z_2)$ and output $y = f(z)$.

**P.2 (price taking; MWG 5.C, pp. 135 and 139)** (a) The firm takes factor prices $r$ as given and chooses inputs to minimize cost; (b) when profit maximization is discussed, it also takes the output price $p$ as given.

**P.3 (a flexible input; MWG equation (5.C.4), pp. 139-140)** At least one input $z_{\mathrm{v}}$ can be adjusted freely within the period; it takes an interior value in cost minimization, and $f$ is differentiable in it.
:::

P.1 is a regularity condition on technology. P.2(a) is the crucial one: input markets are competitive, and the firm's purchases do not affect input prices; if Plover has bargaining power over creators' pay, it fails. P.2(b) is needed only when discussing supply and the profit function; the markup identity below does not use it. P.3 says that at least one input (for example, freelance creator hours settled monthly) has no adjustment cost; when every input is locked in by long-term contracts, the first-order condition no longer holds period by period.

A price-taking firm can use the cost function to write profit maximization as a one-dimensional problem:

$$\pi(p, r) = \max_{y \geq 0}\; \big[\, p\, y - c(r, y)\,\big]$$

The first-order condition of an interior solution is $p = \partial c(r, y)/\partial y$: price equals marginal cost.

::: {.theorem}
**Theorem 13 (the profit function, the cost function, and two lemmas; MWG Propositions 5.C.1 (ii)(vi), p. 138; 5.C.2 (ii)(vi), p. 141)** Let P.1 hold. (i) Under P.2(b), the profit function $\pi(p, r)$ is a convex function of $(p, r)$; if the optimal output and inputs are unique, Hotelling's lemma holds:

$$\frac{\partial \pi(p, r)}{\partial p} = y(p, r), \qquad \frac{\partial \pi(p, r)}{\partial r_\ell} = -\, z_\ell(p, r)$$

(ii) Under P.2(a), the cost function

$$c(r, y) = \min_{z \geq 0}\; r \cdot z \quad \text{s.t.} \quad f(z) \geq y$$

is a concave function of $r$; if the conditional factor demand $z(r, y)$ is unique, Shephard's lemma holds:

$$\frac{\partial c(r, y)}{\partial r_\ell} = z_\ell(r, y)$$
:::

The proof sketch has three steps. Step one (P.2): the profit function is the supremum of the family of functions $p\,y - r \cdot z$, each linear in prices, over feasible $(y, z)$, and a supremum of linear functions is convex; the cost function is the infimum of a family of functions linear in $r$, and so is concave, by the same argument as in step three of Theorem 3. Step two (P.1, unique optimum): prices appear only in the objective, so the derivative of the value function with respect to a price equals the direct derivative of the objective, which gives the two lemmas; MWG also notes that they follow directly from the duality theorem (Proposition 3.F.1), a route that does not require the optimal solution to be differentiable. Step three (P.1, P.2(a)): cost minimization is the same mathematical problem as the EMP; replace $u$ with $f$, $\bar u$ with $y$, and $p$ with $r$, and Theorems 3 to 6 translate word for word, with conditional factor demand corresponding to Hicksian demand.

In words: the direction in which a value function bends delivers comparative statics directly. Convexity of the profit function implies that supply is nondecreasing in its own price and that factor demand is nonincreasing in its own price; concavity of the cost function plus Shephard's lemma implies that the Jacobian of conditional factor demand is symmetric and negative semidefinite, which can be tested against an estimated cross-substitution matrix. A price-taking firm chooses only quantities, and the conclusions are determined by technology and prices; when P.2 fails (the firm affects prices), the profit function is no longer a value function of prices, Hotelling's lemma does not hold, and Chapter 4 covers pricing under imperfect competition. Returns to scale connect technology to market structure: when the production function is homogeneous of degree one, cost is linear in output (MWG Proposition 5.C.2 (viii)), and the scale of a competitive firm is indeterminate; when $f$ is concave, cost is convex in output and marginal cost is nondecreasing (ix), so the supply curve slopes up; digital products with large fixed costs and near-zero marginal cost have increasing returns, under which a competitive equilibrium usually does not exist and the market tips toward a few large firms. Plover's numbers: the cost function $c(r, y) = 2y\sqrt{r_1 r_2}$ gives USD 40, Shephard's lemma gives the conditional factor demand $z = (5, 20)$, and substituting back, $4 \times 5 + 1 \times 20 = 40$, agrees with the cost function.

Cost minimization also hides an identity that is very useful for studying market power. Suppose Plover sells the delivered content at USD 6 a unit; marginal cost is $2\sqrt{r_1 r_2} = 4$ dollars, so price is 1.5 times marginal cost. Without looking at the demand curve, can this 1.5 be read off production data?

::: {.theorem}
**Theorem 14 (the markup identity from cost minimization; Hall 1988; De Loecker and Warzynski 2012; derived from MWG equation (5.C.4), pp. 139-140)** Let P.1, P.2(a), and P.3 hold. Write marginal cost as $mc = \partial c(r, y)/\partial y$, the markup ratio as $\mathcal{M} = p/mc$, the output elasticity of input $z_{\mathrm{v}}$ as $\mathcal{E}_{y\mathrm{v}} = \partial \ln f/\partial \ln z_{\mathrm{v}}$, and the share of its expenditure in revenue as $s_{\mathrm{v}} = r_{\mathrm{v}} z_{\mathrm{v}}/(p\,y)$. Then

$$\mathcal{M} = \frac{\mathcal{E}_{y\mathrm{v}}}{s_{\mathrm{v}}}$$

No assumption is made about the output market; the Lerner index is $(p - mc)/p = 1 - 1/\mathcal{M}$.
:::

The proof sketch has two steps. Step one (P.2(a), P.3): at an interior point the first-order condition of cost minimization holds with equality, $r_{\mathrm{v}} = mc \cdot \partial f/\partial z_{\mathrm{v}}$, and the multiplier is marginal cost (MWG equation (5.C.4) and the remarks after it). Step two (algebra): multiply both sides by $z_{\mathrm{v}}/(p\,y)$; the left side is $s_{\mathrm{v}}$ and the right side is $(mc/p)\,\mathcal{E}_{y\mathrm{v}}$, and rearranging gives the result.

In words: as long as one input can be adjusted flexibly and input markets are competitive, the markup ratio can be read off production data: an elasticity divided by a share, requiring neither a demand curve nor knowledge of how the firm competes in the output market. Its credibility rests on two primitives: the output elasticity comes from technology, the share from accounting data. When P.3 fails, the first-order condition does not hold period by period; when P.2(a) fails, $r_{\mathrm{v}}$ must be replaced by marginal expenditure, and the identity folds input-market power into $\mathcal{M}$, overstating it. The share must be input expenditure divided by revenue; dividing by cost instead gives the ratio of average to marginal cost, the scale elasticity, not $\mathcal{M}$. This book does not cover production function estimation; Section 3.1 of Chapter 29 takes another road, backing out marginal cost from the first-order conditions of demand and pricing, and the two roads can serve as checks on each other. Plover's numbers: $\mathcal{M} = 6/4 = 1.5$, and the Lerner index is 1/3; the output elasticity of $z_2$ is 1/2, and its expenditure of USD 20 is 1/3 of the USD 60 of revenue, so the ratio of the two is also 1.5.

### 4.3 Choice Under Uncertainty

When Plover sells a reading tablet, it offers a one-year protection plan with it: the probability of a claim within the year is 0.08, a claim costs USD 150, and the expected loss is USD 12. The W3 subscriber buys the plan when it sells for USD 20 and does not buy it at USD 24. Making sense of paying USD 20 to cover a risk whose expected cost is USD 12 calls for expected utility. This section writes lotteries (probability distributions over outcomes) as $F$, $F'$, $F''$, and $\alpha \in [0, 1]$ is the probability with which two lotteries are mixed, unrelated to the scaling constant of Section 3.1.

::: {.assumption}
**U.1 (rational and continuous preferences over lotteries; MWG 6.B, Definition 6.B.3, p. 171)** The set of outcomes is finite, and preferences $\succsim$ over lotteries are complete, transitive, and continuous: for any $F$, $F'$, $F''$, the condition $\alpha F + (1 - \alpha) F' \succsim F''$ holds on a closed set of $\alpha$, and the condition $F'' \succsim \alpha F + (1 - \alpha) F'$ also holds on a closed set of $\alpha$.

**U.2 (independence; MWG Definition 6.B.4, p. 171)** For any $F$, $F'$, $F''$ and $\alpha \in (0, 1)$, $F \succsim F'$ if and only if $\alpha F + (1 - \alpha) F'' \succsim \alpha F' + (1 - \alpha) F''$.
:::

The continuity in U.1 rules out breaks of the kind "the slightest risk of death makes any lottery unacceptable." U.2 is the core: mixing each of two lotteries with a third in the same proportion does not change their ranking; it allows any attitude toward risk but requires preferences to be linear in probabilities. Experiments such as the Allais paradox show again and again that real people violate it systematically (Section 7).

::: {.theorem}
**Theorem 15 (expected utility theorem; MWG Propositions 6.B.3, p. 176; 6.B.2, p. 173)** Under U.1 and U.2, there is a Bernoulli utility function $u$ assigning a value to each outcome such that, for any two lotteries, $F \succsim F'$ if and only if $\int u\, dF \geq \int u\, dF'$ (with finitely many outcomes, the integral is the probability-weighted sum). Another Bernoulli utility represents the same preferences in expected-utility form if and only if it is a positive affine transformation of $u$ (multiplication by a positive number plus a constant).
:::

The proof sketch has three steps and proves only existence, following MWG (pp. 176-178, which for simplicity assumes that a best and a worst lottery exist). Step one (U.1, U.2): for each lottery $F$ there is a unique $\alpha_F \in [0, 1]$ such that $F$ is indifferent to the mixture "the best lottery with probability $\alpha_F$, otherwise the worst lottery"; continuity gives existence, and independence gives uniqueness. Step two (U.2): in such mixtures, a higher probability on the best lottery is better, so $\alpha_F$ represents the preferences. Step three (U.2): independence ensures that $\alpha_F$ is linear in mixtures of lotteries, so the value of two lotteries mixed with weights $\alpha$ and $1 - \alpha$ equals the same weighting of their values; on finitely many outcomes such a linear function is the probability-weighted sum, and the weight attached to each outcome is $u$ (MWG Proposition 6.B.1, p. 173).

In words: under U.1 and U.2, preferences over lotteries are characterized by the expectation of a Bernoulli utility, the result of von Neumann and Morgenstern (1947); the expectation is of $u$, not of money, so the shape of $u$ encodes the attitude toward risk. This $u$ is unique up to positive affine transformations, so its curvature is meaningful, whereas the curvature of the ordinal utility of Section 3.1 is not. When independence fails, preferences are no longer linear in probabilities, and models such as prospect theory (Kahneman and Tversky 1979) weight probabilities and count gains and losses from a reference point. Testing independence requires observing the same person's choices among mixed lotteries. Plover's numbers: not buying the plan is a lottery that loses USD 150 with probability 0.08, buying it is paying the plan's price for sure, and expected utility reduces the comparison of the two to a single $u$; under CARA (constant absolute risk aversion), the threshold at a plan price of USD 20 sits at an Arrow-Pratt coefficient of absolute risk aversion (defined below) of $A = 0.00714$ (Section 6.6).

The concavity or convexity of $u$ is the attitude toward risk. When $u$ is concave, Jensen's inequality makes the expected utility of a lottery lower than the utility of its expected value, and the person is risk averse. The degree of concavity must be measured by a quantity immune to positive affine transformations, namely the Arrow-Pratt coefficient of absolute risk aversion

$$A(w) = -\frac{u''(w)}{u'(w)}$$

Two quantities translate risk attitude into money: the certainty equivalent $\mathrm{CE}$, defined by $u(\mathrm{CE}) = \mathbb{E}[u(w)]$, and the risk premium $\pi = \mathbb{E}[w] - \mathrm{CE}$, the expected wealth one is willing to give up to avoid the risk (unrelated to the profit function $\pi(p, r)$ of Section 4.2). For wealth $w$ with a small zero-mean risk $\widetilde\varepsilon$ added to it:

$$\pi(w)\approx \frac{1}{2}A(w)\operatorname{Var}(\widetilde\varepsilon).$$

Suppose a seller's Bernoulli utility is $u(x) = \sqrt{x}$ and she faces a 50-50 lottery in which her wealth is either 144 or 64: expected wealth is 104, expected utility 10, the certainty equivalent 100, and the risk premium 4. The approximation, evaluated at 104, gives $A \approx 0.0048$, a variance of 1600, and an approximate premium of 3.846; the gap arises because this risk is not really small. The theorem below states the equivalent meanings of "more risk averse"; write $A(x; u)$ for the coefficient of utility $u$ at wealth $x$.

::: {.theorem}
**Theorem 16 (equivalent characterizations of comparative risk aversion; MWG Proposition 6.C.2, p. 191; the Arrow-Pratt coefficient, Definition 6.C.3, p. 190)** Let $u_1$ and $u_2$ be two increasing, twice differentiable Bernoulli utilities. The following five statements are equivalent: (i) for all $x$, $A(x; u_2) \geq A(x; u_1)$; (ii) $u_2$ is an increasing concave transformation of $u_1$; (iii) for any distribution of payoffs, the certainty equivalent under $u_2$ is no higher than under $u_1$; (iv) for any wealth and any symmetric risk of any size, the probability premium required by $u_2$ is no lower than that required by $u_1$; (v) whenever $u_2$ finds a lottery no worse than some sure amount of wealth, so does $u_1$.
:::

The proof sketch has two steps and proves only the equivalence of (i) and (ii); for the rest, see MWG Exercises 6.C.6 and 6.C.7. Step one (increasing): the two utilities rank money the same way, so $u_2$ can always be written as some increasing function of $u_1$. Step two (twice differentiable): differentiating this composite function twice, the coefficient of $u_2$ equals the coefficient of $u_1$ minus "the ratio of the second to the first derivative of the transformation, times $u_1'$"; since $u_1' > 0$, (i) holds if and only if the second derivative of the transformation is nonpositive, which is (ii).

In words: within the expected-utility framework, the five ways of saying "more risk averse" are equivalent, and the measure is the Arrow-Pratt coefficient, whose division by $u'$ makes it immune to the scale of utility. The premium for a small risk is determined by curvature and variance; if absolute risk aversion decreases with wealth, the higher the wealth, the lower the premium demanded for the same absolute risk. The approximation holds only for small risks and is biased when the risk is not small; in the seller example above it understates the premium. Testing risk attitudes requires looking at the same person's choices over risks of different sizes, and separating risk preferences from beliefs and borrowing constraints. Plover's numbers: by (i) and (iii), the larger $A$, the more she is willing to pay for the plan; the W3 subscriber buys at USD 20 and not at USD 24, so $A$ lies in $[0.00714, 0.00967]$.

Two prospects can also be compared without specifying a particular $u$.

::: {.theorem}
**Theorem 17 (stochastic dominance and expected utility; MWG Propositions 6.D.1, p. 195; 6.D.2, p. 199)** Let the payoff distribution functions $F$ and $G$ of two lotteries both have support in $[0, x_{\max}]$. (i) Every expected-utility maximizer with a nondecreasing Bernoulli utility weakly prefers $F$ ($F$ first-order stochastically dominates $G$) if and only if, for all $x$, $F(x) \leq G(x)$. (ii) If $F$ and $G$ have the same mean, the following three statements are equivalent: every expected-utility maximizer with a nondecreasing and concave Bernoulli utility weakly prefers $F$ (second-order stochastic dominance); $G$ is a mean-preserving spread of $F$; for all $x$, $\int_0^x G(\tilde x)\, d\tilde x \geq \int_0^x F(\tilde x)\, d\tilde x$.
:::

The proof sketch has two steps. Step one ((i), using monotonicity only): integrating by parts, $\int u\, dF - \int u\, dG = -\int u'(x)\,\big[F(x) - G(x)\big]\, dx$, which is nonnegative when $u' \geq 0$ and $F \leq G$ everywhere; conversely, if $F > G$ somewhere, a step utility that jumps at that point gives a counterexample. Step two ((ii), adding concavity): integrating by parts once more, the difference becomes $\int u''(x)\,\big[\int_0^x (F - G)\big]\, dx$ plus a boundary term; equal means make the boundary term zero, and $u'' \leq 0$ together with the integral condition gives nonnegativity. MWG states (ii) without proof; for the full argument, see Rothschild and Stiglitz (1970).

In words: once a dominance relation can be verified, there is no need to know the particular utility to conclude how a whole class of people ranks two prospects: first-order dominance holds for everyone who prefers more to less, and second-order dominance for every risk-averse person. Second-order dominance must be stated for equal means; when the means differ, "the tighter distribution is better" no longer holds for every risk averter. The criterion is determined by the two distributions alone, which makes it a robust tool for welfare comparisons, for example comparing the risk in creators' income under two revenue-sharing schemes; testing it needs only samples from the two distributions, with no preference estimation. Plover's numbers: of two lotteries that both have mean 104, the tighter $\{100, 108\}$ has expected utility 10.196 and the more dispersed $\{64, 144\}$ has expected utility 10.000; the concave utility penalizes the more dispersed one.

Back to the protection plan. If the subscriber's Bernoulli utility has the CARA form $u(x) = -\exp(-A x)$, wealth plays no role: buying at USD 20 shows that $A \geq 0.00714$, and not buying at USD 24 shows that $A \leq 0.00967$. What does this interval imply for other bets? Write a 50-50 bet as "lose $L$ dollars, win $g$ dollars."

::: {.theorem}
**Theorem 18 (calibration under CARA; a direct computation in this chapter, with the general case in Rabin 2000)** Let the Bernoulli utility be $u(x) = -\exp(-A x)$, $A > 0$. The 50-50 bet "lose $L$, win $g$" is rejected (indifference included) if and only if $\tfrac{1}{2}\exp(A L) + \tfrac{1}{2}\exp(-A g) \geq 1$. Hence it is rejected for every $g$ if and only if $L \geq (\ln 2)/A$.
:::

The proof sketch has two steps. Step one (CARA): wealth $w$ enters $u(w - L)$, $u(w + g)$, and $u(w)$ through the same factor $\exp(-A w)$, which cancels to leave the rejection condition, so wealth plays no role. Step two ($A > 0$): the left side decreases in $g$ and, as $g$ goes to infinity, tends to $\tfrac{1}{2}\exp(A L)$; when this limit is at least 1, the bet is rejected for every $g$, and when it is below 1, the bet is accepted for $g$ large enough, so the boundary is exactly $L = (\ln 2)/A$.

In words: under CARA, rejecting a small favorable bet pins down a lower bound on $A$, and with it a cap on acceptable losses: a 50-50 bet whose loss exceeds $(\ln 2)/A$ is rejected however much it pays. This is a directly computable special case of the calibration theorem of Rabin (2000). Rabin's theorem requires only that $u$ be strictly increasing and weakly concave: if a person rejects some small favorable 50-50 bet at every wealth level within a range, she must reject some 50-50 bets with moderate losses and enormous gains; if she rejects it at all wealth levels, some bets are rejected whatever the gain. This book only states the theorem, because the proof of the general case accumulates the decline in marginal utility piece by piece along wealth and takes considerable space. The lesson is that aversion to small risks, if attributed to the curvature of utility over wealth, implies absurd behavior toward large bets, so it more likely comes from loss aversion, narrow framing, or misjudged probabilities. $A$ is revealed by choices at insurance prices, and testing the conclusion requires the same person's choices across different bets. Plover's numbers: rejecting "lose 10, win 11" requires $A \geq 0.00908$, whereas $(\ln 2)/100 = 0.00693$, so any CARA user who rejects this small bet rejects every 50-50 bet that loses USD 100; at $A = 0.00908$, a 50-50 bet that loses more than USD 76.30 is rejected however much it pays. The W3 subscriber's interval $[0.00714, 0.00967]$ contains 0.00908, so whether she rejects "lose 10, win 11" cannot be determined, but the lower bound already exceeds 0.00693 (Section 6.6).

### 4.4 How Comparative Statics Become Empirical Predictions

Theory is useful for empirical work because it constrains the relationships among several observable quantities at once. The table below writes several of the chapter's comparative statics as predictions that can be taken to data.

| Theoretical change | Comparative static | Observable implication | Alternative explanations not identified by sign alone |
|---|---|---|---|
| Own price rises | Hicksian demand has only substitution; Walrasian demand also has an income effect | Groups with a higher budget share or a stronger income response react more | liquidity, search, salience |
| Input price rises | Conditional factor demand moves along the cost gradient | Conditional demand for substitutable inputs adjusts in the opposite direction | technical change, quantity discounts, endogenous procurement prices |
| Risk variance rises | Risk premium rises by about $\tfrac{1}{2}A(w)\,\Delta\operatorname{Var}$ | Risk averters demand higher compensation | beliefs, ambiguity, borrowing constraints |
| Digital good's price falls to zero | The gaps among CV, EV, and $\Delta \mathrm{CS}$ depend on income effects | Heavy users' welfare gains are larger, but usage time alone cannot value them | selection into use, zero-price privacy cost |

These predictions should be treated as a set of joint restrictions, with a clear statement of what exogenous variation moves prices, income, risk, or factor costs. A sign that goes the right way in a cross section is only evidence consistent with a mechanism; identifying the mechanism takes variation that can separate competing explanations, such as price experiments, cost shifters, or randomized information interventions.

## 5 Anchoring Papers

Three studies put the chapter's welfare machine to real use: the valuation of a new good, of the variety of online book catalogs, and of zero-price digital products. Each is laid out by paper, method, data, results, and limitations, with the focus on how the welfare number is computed from demand or from choices.

### 5.1 Hausman (1996)

::: {.case}
Paper: "Valuation of New Goods under Perfect and Imperfect Competition," in *The Economics of New Goods*, edited by Bresnahan and Gordon. It writes the social value of a new good as this chapter's CV and, on that basis, argues that a price index that ignores new goods overstates the rise in the cost of living.

Method: a new good coming into existence is equivalent to its price falling from the virtual price (the reservation price at which demand is exactly zero) to the actual price, and the welfare gain is the area under the demand curve over this price interval (Theorems 10 to 12). Computing this area requires estimating the whole demand curve, not just the point at the current price: Hausman estimates the demand for ready-to-eat cereal with a multilevel demand system, using the prices of the same product in other cities as instruments to deal with price endogeneity (instrumental variables are covered in Chapter 9, and Section 3.3 of Chapter 28 calls instruments of this kind Hausman instruments).

Data: supermarket scanner data, focused on a newly introduced cereal, Apple-Cinnamon Cheerios.

Results: the consumer welfare from this new cereal is substantial, and ignoring new goods makes the category's price index systematically overstate the rise in the cost of living. The methodological contribution is to pin the vague question "what is a new good worth" down to a CV area under a demand curve, which makes it estimable and comparable.

Limitations: the numbers hinge on the functional form of demand and on the extrapolation to the virtual price, while the data are sparsest at the high-price end, exactly where the extrapolation is least grounded; the USD 12 extrapolation of W1 in Section 6.4 is a miniature of this problem. The validity of the instruments depends on demand shocks being uncorrelated across cities.
:::

### 5.2 Brynjolfsson, Hu and Smith (2003)

::: {.case}
Paper: "Consumer Surplus in the Digital Economy: Estimating the Value of Increased Product Variety at Online Booksellers," *Management Science*. The question is a classic of information systems research: online bookstores carry far more titles than physical ones, so how much surplus does this variety bring consumers?

Method: following Hausman's approach to valuing new goods, books that rank beyond 100,000 in Amazon's sales rankings, and that physical stores usually do not stock, are treated as new goods whose price falls from the virtual price to the online price; demand takes a log-linear form with income effects dropped, so CV has the closed form of Hausman (1981). The price elasticity is not estimated from demand data but inferred from publishers' gross margins via the Lerner index; sales of obscure titles are inferred from a log-linear regression of sales on sales rank.

Data: weekly sales of 321 titles from one publisher together with their Amazon sales ranks over the same period, plus industry statistics on gross margins and sales.

Results: in 2000, the consumer surplus from catalog variety was between USD 731 million and USD 1.03 billion, 7.3 to 10.0 times the surplus from the lower book prices that online competition brought; the bulk of the welfare comes from the availability of obscure titles, not from prices.

Limitations: log-linear demand puts no cap on valuations, and excluding consumers whose valuations exceed 5 times the current price cuts the surplus by 28.0% to 40.6%; the elasticity comes from gross margins, and obscure titles sit in the tail where the rank data are sparsest. Summing the surplus from obscure titles also amounts to assuming that these books really cannot be bought offline and that consumers can find them, so the roles of search costs and recommendations are all packed into this assumption.
:::

### 5.3 Brynjolfsson, Collis and Eggers (2019)

::: {.case}
Paper: "Using massive online choice experiments to measure changes in well-being," *PNAS*. Zero-price digital products such as search, social media, and maps leave almost no trace in GDP, yet they clearly bring enormous welfare; how should this welfare be measured?

Method: this chapter's EV. The researchers ask users "how much would we have to pay you to give up a product for a period of time," and the minimum compensation a user will accept leaves her indifferent between "the product without compensation" and "the compensation without the product"; by definition this is the EV of having the product (relative to not having it), that is, the CV of losing it. Summing individual compensations gives the product's total consumer welfare.

Data: several rounds of large-scale choice experiments covering digital products such as Facebook, search engines, email, and maps.

Results: users value these zero-price products highly; measured by EV, the consumer welfare of free digital products is vast, and it escapes spending-based GDP entirely.

Limitations: reported compensation demands carry hypothetical bias and the endowment effect, individual valuations are highly dispersed, and extrapolating the median to the population calls for care; these numbers measure private welfare, which need not equal social welfare once externalities and attention costs are counted in.
:::

Together the three papers show that CV and EV are the rigorous way to answer "what is a new product, a larger catalog, or a free service worth"; the first two compute the answer from demand curves, while the third elicits compensation directly in choice experiments, sidestepping functional form at the price of reporting bias.

## 6 Running Case: Settling Plover's Accounts

The code below runs under R 4.6.1. Only the key statements are shown, with several statements sometimes merged onto one line; every number in the text comes from the output of an actual run, and the result lines in the code that begin with `#` are copied verbatim from what the run printed.

### 6.1 W0: Demand, the Expenditure Function, and the Slutsky Decomposition

The first-order conditions of Cobb-Douglas $u = x_1^{2/5} x_2^{3/5}$ give constant expenditure shares: the fraction of the budget spent on each good always equals its exponent, so

$$x_1(p, w) = \frac{2}{5}\frac{w}{p_1}, \qquad x_2(p, w) = \frac{3}{5}\frac{w}{p_2}$$

At the initial prices $p^0 = (2, 1)$, $x = (20, 60)$ and utility is $u^0 \approx 38.664$; after the increase, at $p^1 = (4, 1)$, $x = (10, 60)$ and $u^1 \approx 29.302$. Spending on premium content is USD 40 both times, and this unmoving 40 is exactly what the first analyst saw: constant shares make spending completely unresponsive to price, and they make the illusion that "unchanged spending means no loss" especially stubborn. The expenditure function is $e(p, u) = u\,(p_1/a)^{a}\,(p_2/(1-a))^{1-a}$ with $a = 2/5$; the two self-checks $e(p^0, u^0) = e(p^1, u^1) = 100$ are the numerical form of the duality identity $e(p, v(p, w)) = w$. The Slutsky decomposition at the initial point is

$$\underbrace{-10}_{\text{total}} = \underbrace{-6}_{\text{substitution}} + \underbrace{(-4)}_{\text{income}}$$

```r
a <- 2 / 5; w <- 100; p2 <- 1; p1_0 <- 2; p1_1 <- 4
x1 <- function(p1, m = w) a * m / p1                                       # Walrasian demand for premium
v  <- function(p1, m = w) x1(p1, m)^a * ((1 - a) * m / p2)^(1 - a)          # indirect utility
e  <- function(p1, u) u * (p1 / a)^a * (p2 / (1 - a))^(1 - a)               # expenditure function
u0 <- v(p1_0); u1 <- v(p1_1)
round(c(x1_0 = x1(p1_0), x1_1 = x1(p1_1), u0 = u0, u1 = u1, e_check_0 = e(p1_0, u0), e_check_1 = e(p1_1, u1)), 4)
total <- -a * w / p1_0^2; income <- -x1(p1_0) * (a / p1_0)                  # Slutsky at the initial point
c(total = total, substitution = total - income, income = income)
#      x1_0      x1_1        u0        u1 e_check_0 e_check_1
#   20.0000   10.0000   38.6636   29.3016  100.0000  100.0000
#        total substitution       income
#          -10           -6           -4
```

### 6.2 W0: The Three Welfare Numbers, Willig, and Vartia

The three areas come from integrating the functions above directly or from differences in expenditure:

$$\Delta\mathrm{CS} = \int_2^4 \frac{40}{p_1}\, dp_1 = 40\ln 2 \approx 27.726$$

$$\mathrm{EV} = e(p^0, u^0) - e(p^0, u^1) = 100 - e(p^0, u^1) \approx 24.214$$

$$\mathrm{CV} = e(p^1, u^0) - e(p^1, u^1) = e(p^1, u^0) - 100 \approx 31.951$$

`vartia_cv` solves the ordinary differential equation of Theorem 12 by fixed-step fourth-order Runge-Kutta, and `willig_cv` is the first-order correction of Theorem 11.

```r
willig_cv <- function(dcs, eta, w) dcs + eta * dcs^2 / (2 * w)       # Willig (1976), first order, loss positive
vartia_cv <- function(x1_fun, p0, p1, w, steps = 2000L) {            # RK4 for dm/dp1 = x1(p1, m), m(p0) = w
  h <- (p1 - p0) / steps; m <- w
  for (k in seq_len(steps)) {
    p <- p0 + (k - 1L) * h
    k1 <- x1_fun(p, m); k2 <- x1_fun(p + h / 2, m + h / 2 * k1)
    k3 <- x1_fun(p + h / 2, m + h / 2 * k2); k4 <- x1_fun(p + h, m + h * k3)
    m <- m + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
  }
  m - w
}
EV <- w - e(p1_0, u1); dCS <- a * w * log(p1_1 / p1_0); CV <- e(p1_1, u0) - w
round(c(EV = EV, dCS = dCS, CV = CV, willig = willig_cv(dCS, 1, w), willig_ev = willig_cv(dCS, -1, w),
        vartia = vartia_cv(x1, p1_0, p1_1, w)), 4)
#        EV       dCS        CV    willig willig_ev    vartia
#   24.2142   27.7259   31.9508   31.5695   23.8823   31.9508
```

The ordering $\mathrm{EV} < \Delta\mathrm{CS} < \mathrm{CV}$ fits Theorem 10 exactly; the Willig correction pushes consumer surplus up to USD 31.57, and the Vartia equation recovers the exact CV from Walrasian demand. The disagreement between the two analysts of Section 1 is now resolved: the first confused spending with welfare; the second analyst's number is a good approximation, but it is neither CV nor EV, and which one to report depends on whether the operations team is asking to "compensate to no loss" or to "buy out in advance."

### 6.3 W0's Remaining Checks: Integrability, Aggregation, Production, and Risk

The remaining numbers of Sections 3 and 4 come from the lines below: the Slutsky matrix at $(p^0, w)$, which is symmetric, has two nonpositive eigenvalues, and gives zero when multiplied by the price vector; the two linear-expenditure-system consumers, for whom redistribution leaves total demand unmoved under a common slope and moves it when the slopes differ; the cost function, conditional factor demand, and the markup identity; and the lottery under $\sqrt{x}$ utility, together with the comparison of two lotteries that have the same mean.

```r
x1s <- a * w / p1_0; x2s <- (1 - a) * w / p2
S <- matrix(c(-a * w / p1_0^2 + x1s * (a / p1_0), 0 + x2s * (a / p1_0),
              0 + x1s * ((1 - a) / p2), -(1 - a) * w / p2^2 + x2s * ((1 - a) / p2)), nrow = 2, byrow = TRUE)
S
c(sym = S[1, 2] - S[2, 1], Sp = (S %*% c(p1_0, p2))[1], eigen = zapsmall(eigen(S, symmetric = TRUE)$values))
#      [,1] [,2]
# [1,]   -6   12
# [2,]   12  -24
#    sym     Sp eigen1 eigen2
#      0      0      0    -30
x1_les <- function(g1, g2, wi, b = 2 / 5) g1 + b * (wi - p1_0 * g1 - p2 * g2) / p1_0   # LES demand for premium
c(X1 = x1_les(5, 10, 100) + x1_les(15, 0, 200), X1_redis = x1_les(5, 10, 150) + x1_les(15, 0, 150),
  X1_b07 = x1_les(5, 10, 100) + x1_les(15, 0, 200, 0.7), X1_b07_redis = x1_les(5, 10, 150) + x1_les(15, 0, 150, 0.7))
#           X1     X1_redis       X1_b07 X1_b07_redis
#         70.0         70.0         95.5         88.0
r1 <- 4; r2 <- 1; y <- 10                                    # y = z1^(1/2) z2^(1/2)
c(cost = 2 * y * sqrt(r1 * r2), z1 = y * sqrt(r2 / r1), z2 = y * sqrt(r1 / r2))   # Shephard: z = dc/dr
# cost   z1   z2
#   40    5   20
p_y <- 6; mc <- 2 * sqrt(r1 * r2); z2 <- y * sqrt(r1 / r2)
c(markup = p_y / mc, elasticity_z2 = 1 / 2, revenue_share_z2 = r2 * z2 / (p_y * y), ratio = (1 / 2) / (r2 * z2 / (p_y * y)))
#           markup    elasticity_z2 revenue_share_z2            ratio
#        1.5000000        0.5000000        0.3333333        1.5000000
Ew <- 0.5 * 64 + 0.5 * 144; Eu <- 0.5 * sqrt(64) + 0.5 * sqrt(144)          # u(x) = sqrt(x), lottery {64, 144}
c(Ew = Ew, CE = Eu^2, premium = Ew - Eu^2, approx = 0.5 * (1 / (2 * Ew)) * (0.5 * (64 - Ew)^2 + 0.5 * (144 - Ew)^2))
#         Ew         CE    premium     approx
# 104.000000 100.000000   4.000000   3.846154
c(EU_F = 0.5 * sqrt(100) + 0.5 * sqrt(108), EU_G = Eu)        # F = {100, 108} against G = {64, 144}, both mean 104
#     EU_F     EU_G
# 10.19615 10.00000
```

The diagonal element $-6$ is exactly the substitution effect of Section 6.1, and the equal off-diagonal elements express the symmetry; $S p = 0$ is a corollary of the homogeneity of degree zero of compensated demand. All three checks pass everywhere ($S(p, w) = a(1-a)\,w \begin{pmatrix} -1/p_1^2 & 1/(p_1 p_2) \\ 1/(p_1 p_2) & -1/p_2^2 \end{pmatrix}$), so this demand system is integrable, and the expenditure function of Section 6.1 is exactly what it recovers.

### 6.4 W1: A Price Experiment on Ten Thousand Subscribers

W1 has 10,000 subscribers, each with Stone-Geary preferences (the utility behind the linear expenditure system) $u_i(x) = \beta_i \ln(x_1 - \gamma_{1i}) + (1 - \beta_i)\ln(x_2 - \gamma_{2i})$ and demand $x_{1i} = \gamma_{1i} + \beta_i\,(w_i - p \cdot \gamma_i)/p_1$; here $\beta_i$ is the marginal budget share on premium content of subscriber $i$, and $\gamma_i = (\gamma_{1i}, \gamma_{2i})$ is the vector of subsistence quantities. The monthly budget $w_i$ is lognormal with a median of USD 100; $\beta_i$ follows a Beta distribution with mean 0.35 and standard deviation 0.104, independent of the budget, and this heterogeneity breaks the Gorman condition. The premium price is assigned by complete randomization among USD 2, 2.5, 3, 3.5, and 4, with 2,000 subscribers at each price. The variant W1g has the same subscribers with every $\beta_i$ set to 0.35, so that the Gorman condition holds. The welfare question is the same as for W0, the average loss from a price increase from USD 2 to USD 4, plus one question outside the experiment: the loss from an increase to USD 12. In the code, `sg_demand` and `sg_welfare` are the Stone-Geary demand and the closed-form CV, EV, and $\Delta\mathrm{CS}$, `dcs_numeric` computes the area to the left of a demand curve by Simpson's rule, and `fit_demand_forms` runs the three OLS (ordinary least squares) fits described below and returns each form's demand function and income elasticity.

```r
gen_W1 <- function(P, common_beta = FALSE) {
  N <- P$N
  w <- exp(log(P$w_median) + P$w_sdlog * rnorm(N))                   # monthly budget, median 100 USD
  beta <- rbeta(N, P$beta_phi * P$beta_mean, P$beta_phi * (1 - P$beta_mean))
  if (common_beta) beta <- rep(P$beta_mean, N)                        # W1g: Gorman holds
  g1 <- runif(N, P$g1_range[1], P$g1_range[2]); g2 <- runif(N, P$g2_range[1], P$g2_range[2])
  p1 <- sample(rep(P$p_grid, each = N / length(P$p_grid)))            # 2,000 subscribers at each price
  data.frame(p1 = p1, w = w, x1 = sg_demand(p1, w, beta, g1, g2, P$p2),
             beta_true = beta, g1_true = g1, g2_true = g2)
}
P1 <- list(N = 10000, p_grid = c(2, 2.5, 3, 3.5, 4), p2 = 1, w_median = 100, w_sdlog = 0.3,
           beta_mean = 0.35, beta_phi = 20, g1_range = c(0.5, 1.5), g2_range = c(0, 10))
set.seed(1101); W1 <- gen_W1(P1)
set.seed(1101); W1g <- gen_W1(P1, common_beta = TRUE)                 # the same subscribers, common beta
s <- with(W1, sg_welfare(2, 4, w, beta_true, g1_true, g2_true))        # truth, closed form per subscriber
cat(sprintf("truth: mean EV %.4f  dCS %.4f  CV %.4f;  EV <= dCS <= CV for all: %s\n",
            mean(s$ev), mean(s$dcs), mean(s$cv), all(s$ev <= s$dcs & s$dcs <= s$cv)))
# truth: mean EV 22.4749  dCS 25.5530  CV 29.2318;  EV <= dCS <= CV for all: TRUE
```

The analyst sees only each subscriber's price, budget, and usage, and fits three demand forms by OLS. The linear form $x_1 = \beta_0 + \beta_{\mathrm{p}} p_1 + \beta_{\mathrm{w}} w + \beta_{\mathrm{pw}}\, p_1 w$ includes a price-by-budget interaction so that the income effect can vary with price; the log-linear form $\ln x_1 = \beta_0 + \beta_{\mathrm{p}} \ln p_1 + \beta_{\mathrm{w}} \ln w$ recovers levels with Duan's smearing (multiplying back the mean of the exponentiated residuals); the income-term form $x_1 = \beta_0 + (\beta_{\mathrm{p}} + \beta_{\mathrm{w}} w)/p_1$, in which the budget enters through an income term, is a linear expenditure system, and its family of functions contains the true average demand. All three forms include the budget, because computing CV requires knowing the income effect. For each form, $\Delta\mathrm{CS}$, the Willig correction, and the Vartia CV are computed at each subscriber's own budget and then averaged; $\mathrm{CV}^{\mathrm{rep}}$ treats the fitted demand at the average budget as that of a single person.

```r
fm <- fit_demand_forms(W1)     # OLS: x1 ~ p1 * w; log(x1) ~ log(p1) + log(w); x1 ~ I(1 / p1) + I(w / p1)
for (nm in names(fm)) {
  f <- fm[[nm]]; w_i <- W1$w
  dcs_i <- dcs_numeric(f$x1_fun, 2, 4, w_i)                           # area under the fitted demand
  cat(sprintf("%-9s dCS %.4f  Willig %.4f  CV %.4f  CV^rep %.4f  CV 2 -> 12 %.4f\n", nm, mean(dcs_i),
              mean(willig_cv(dcs_i, f$eta_fun(2, w_i), w_i)), mean(vartia_cv(f$x1_fun, 2, 4, w_i)),
              vartia_cv(f$x1_fun, 2, 4, mean(w_i)), mean(vartia_cv(f$x1_fun, 2, 12, w_i))))
}
sg <- with(W1g, sg_welfare(2, 4, w, beta_true, g1_true, g2_true))
cat(sprintf("W1g: true mean CV %.4f, income-form CV^rep %.4f\n", mean(sg$cv),
            vartia_cv(fit_demand_forms(W1g)$income$x1_fun, 2, 4, mean(W1g$w))))
# linear    dCS 26.0475  Willig 29.3028  CV 29.5797  CV^rep 29.5797  CV 2 -> 12 41.6533
# loglinear dCS 25.5345  Willig 28.6178  CV 28.8793  CV^rep 28.8910  CV 2 -> 12 96.7118
# income    dCS 25.5428  Willig 28.7152  CV 28.8966  CV^rep 28.8966  CV 2 -> 12 96.4233
# W1g: true mean CV 28.8795, income-form CV^rep 28.8746
```

Start inside the experimental range. The true average losses are EV USD 22.47, $\Delta\mathrm{CS}$ USD 25.55, and CV USD 29.23, and every subscriber satisfies the ordering of Theorem 10; their average income elasticity is 1.009, and $\Delta\mathrm{CS}$ averages 0.244 of the budget, exactly the kind of non-negligible income effect the Willig correction is meant to handle. With the income-term form, the analyst's $\Delta\mathrm{CS}$ is USD 25.54, which falls USD 3.69 short of the true average CV; the Willig correction raises it to USD 28.72, closing 86.0% of this error; Vartia gives USD 28.90, which still falls 1.15% short of the truth. The remaining gap is not the fault of Willig or Vartia: at the average budget the income-term form is a representative consumer, its Vartia CV equals $\mathrm{CV}^{\mathrm{rep}}$ exactly, and the heterogeneity of $\beta_i$ makes the representative consumer understate the average CV (Section 3.7). Run the same analysis on W1g, and the representative consumer falls only 0.017% short. The three forms' CVs within the range are USD 29.58, 28.88, and 28.90, and the gap between the largest and the smallest is only 2.4% of their average: within the range, functional form hardly matters.

![Average premium demand in W1: the sample means of the five price arms (points), the true average demand, and the three fitted forms, drawn out to USD 12; the shaded band marks the prices the experiment used. The linear form hits zero at about USD 6.06, while the other two forms follow the true demand.](assets/fig/fig_01_extrapolation.svg)

Now step outside the range. Extrapolating the price to USD 12, the true average loss is USD 98.77; the log-linear and income-term forms give USD 96.71 and USD 96.42, lower by 2.09% and lower by 2.38%; the linear form gives only USD 41.65, lower by 57.8%, and the spread across the three forms jumps to 70.4% of their average. The disagreement comes almost entirely from the linear form: at the average budget its demand hits zero at about USD 6.06, imposing a choke price that the data cannot see, so all of the loss from prices above USD 6 is cut off. The risk of extrapolation lies not in "leaving the data" but in the functional form making, beyond the data, an assumption on your behalf that you have never tested.

### 6.5 W2: The Revealed-Preference Panel

Of the 300 users, 210 drawn at random choose exactly according to their own CES preferences, and the other 90 add noise along the budget line: a normal disturbance with standard deviation 0.5 added to the log odds of the premium share. The Bronars benchmark places 10 uniformly random choosers on each user's own set of budgets, 3,000 in all. `gen_W2` and `gen_W2_bronars` generate the two sets of choices under these settings; `garp_e` builds the direct revealed-preference matrix, computes its transitive closure with Warshall's algorithm, and then looks for cycles with a strict link; `ccei_e` bisects on the efficiency level.

```r
garp_e <- function(P, X, e = 1) {                  # rows of P, X: observations t; E[t, s] = p^t . x^s
  E <- P %*% t(X); own <- diag(E)
  R <- e * own >= E; Pd <- e * own > E               # direct and strict direct revealed preference
  for (k in seq_len(nrow(R))) R <- R | outer(R[, k], R[k, ], "&")   # Warshall: transitive closure
  sum(R & t(Pd)) == 0                                # GARP: no x^t R x^s with x^s Pd x^t
}
ccei_e <- function(P, X, tol = 1e-8) {               # supremum of the efficiency levels that pass GARP
  if (garp_e(P, X)) return(1)
  lo <- 0; hi <- 1
  while (hi - lo > tol) { mid <- (lo + hi) / 2; if (garp_e(P, X, mid)) lo <- mid else hi <- mid }
  lo
}
P2 <- list(n_users = 300, n_budgets = 8, share_rational = 0.7, m = 60, price_range = c(1, 3),
           alpha_range = c(0.3, 0.7), rho_range = c(-1.5, 0.6), noise_sd = 0.5, n_bronars_per_user = 10)
set.seed(1102); W2 <- gen_W2(P2)                     # CES optimisers and noisy users, 8 budgets each
set.seed(1103); W2b <- gen_W2_bronars(W2, P2$n_bronars_per_user)    # uniform random choosers, same budgets
test <- function(d, id) t(sapply(split(d, d[[id]]), function(g) {
  P <- cbind(g$p1, g$p2); X <- cbind(g$x1, g$x2); c(pass = garp_e(P, X), ccei = ccei_e(P, X)) }))
u <- test(W2, "user"); b <- test(W2b, "chooser"); type <- W2$type_true[W2$t == 1]
cat(sprintf("pass rate: rational %.4f, noisy %.4f (min CCEI %.4f), Bronars %.4f (median CCEI %.4f)\n",
            mean(u[type == "rational", "pass"]), mean(u[type == "noisy", "pass"]), min(u[type == "noisy", "ccei"]),
            mean(b[, "pass"]), median(b[, "ccei"])))
# pass rate: rational 1.0000, noisy 0.6778 (min CCEI 0.8951), Bronars 0.2840 (median CCEI 0.9583)
```

All 210 optimizers pass, with a CCEI of 1 throughout, which is the (i)-implies-(ii) direction of Theorem 8. Of the noisy users, 67.8% pass; those who fail have 3.24 violating pairs on average, the lowest CCEI is 0.895, and 93.3% have an efficiency index above 95%. Of the fully random choosers, 28.4% pass, the median CCEI is 0.958, and 53.6% are above 95%: reporting only that "the CCEI is very close to 1" does not by itself show rationality.

### 6.6 W3: The Protection Plan and Risk Aversion

The plan is offered at four price tiers: USD 15, 20, 24, and 30. Under CARA, not buying at price $g$ amounts to accepting the bet "lose an extra $L = 150 - g$ with probability 0.08, otherwise save $g$," so the threshold coefficient between buying and not buying is computed by one and the same function: `cara_reject_threshold(loss, gain, prob)` finds the positive root of $\text{prob}\cdot\exp(A \cdot \text{loss}) + (1 - \text{prob})\exp(-A \cdot \text{gain}) = 1$; with the rejection condition of Theorem 18 generalized to probability $\text{prob}$, this is the $A$ at which it holds with equality.

```r
loss <- 150; q <- 0.08; prices <- c(15, 20, 24, 30)
A_star <- sapply(prices, function(pr) cara_reject_threshold(loss - pr, pr, prob = q))   # buy at pr iff A >= A_star
A_rabin <- cara_reject_threshold(10, 11)                     # just rejects the 50-50 bet "lose 10, win 11"
round(c(setNames(A_star, prices), A_rabin = A_rabin, ln2_over_100 = log(2) / 100), 6)
round(c(loss_never_at_A_rabin = log(2) / A_rabin, loss_never_at_A_20 = log(2) / A_star[2]), 4)
#           15           20           24           30      A_rabin ln2_over_100
#     0.003168     0.007143     0.009670     0.012865     0.009084     0.006931
# loss_never_at_A_rabin    loss_never_at_A_20
#               76.3039               97.0392
```

A price of USD 20 is 66.7% above the expected loss, and she still buys, so $A \geq 0.00714$; she does not buy at USD 24, so $A \leq 0.00967$. The lower bound is already absurd enough: at $A = 0.00714$, a 50-50 bet that loses more than USD 97.04 is rejected however much it pays. A more credible reading is that her aversion to small risks comes from loss aversion or from overestimating the probability of a claim; using the $A$ calibrated under CARA to predict demand for big-ticket products would go badly wrong.

## 7 Failure Modes and Robustness

The chapter's tools have five common traps. The first is the path dependence of consumer surplus. When the platform changes the prices of two features at once, consumer surplus is a line integral over prices, and its value generally depends on the order in which the prices change; CV and EV are defined by two values of the expenditure function and do not depend on the path. For consumer surplus to escape path dependence, the Walrasian cross-price effects must be symmetric. The symmetry of Theorem 6 governs only compensated demand; by the Slutsky equation, Walrasian symmetry also requires $x_k\,\partial x_\ell/\partial w = x_\ell\,\partial x_k/\partial w$, that is, equal income elasticities among the goods whose prices change. Welfare analysis of multiple price changes should therefore compute CV or EV directly from the expenditure function.

The second is functional form, the shared soft spot of the two papers in Sections 5.1 and 5.2. Welfare needs the whole demand curve, especially the stretch far from observed prices, while the data are dense only near the experimental or current prices. In W1 the log-linear form fits $\ln x_1$, so its $R^2$ cannot be compared directly with that of the other two; what the evidence supports is only that, within the experimental price range, the three forms' CVs differ by only 2.4%, while extrapolated to USD 12 they differ by 70.4%, with nearly all of the gap coming from the linear form. The honest course is to compute the result under several reasonable forms and report a range rather than a point; an experiment that pushes prices into the region to be extrapolated is more reliable than any functional form.

The third is aggregation. Treating aggregate demand as that of a representative consumer is harmless only under the Gorman condition. In W1 the representative consumer falls only 1.15% short of the average CV, because the heterogeneity in shares is modest and the price only doubles; with strong income effects and large differences in shares the bias can be much larger, and a cross section with one price per person cannot reveal it. This is one reason the demand estimation later in the book works so hard to model consumer heterogeneity (the random coefficients of Chapter 28). Quasilinearity is a safe harbor only from the first and third traps: without income effects, the three measures coincide, path dependence disappears, and aggregation becomes clean too. Applied models often assume it because it makes welfare unambiguous, but the assumption holds approximately only when the spending is a small slice of the budget.

The fourth is the power of tests. Passing GARP, or a CCEI close to 1, carries weight only in comparison with a random benchmark on the same set of budgets.

The last is uncertainty. Expected utility relies on the independence axiom, and experiments such as the Allais paradox show that real individuals violate it systematically; Theorem 18 shows from the other side that forcing aversion to small risks into the curvature of utility over wealth implies absurd behavior toward large risks. This chapter uses expected utility because it remains the default starting point of structural modeling; it is an assumption that data can reject, not a logical necessity.

The chapter's welfare machine is precise and useful, but every number it produces carries premises: well-behaved preferences, the right functional form, individuals who can be added up, and rationality under risk. Robust welfare analysis puts these premises on the table and says in which direction the conclusion moves as they loosen.

## 8 Further Reading

::: {.readings}
Required reading, in suggested reading order:

- Nolan Miller (2006, Notes on Microeconomic Theory), Chapters 3 to 6. The intuition for the consumer and uncertainty parts of this chapter is mainly modeled on these notes, which are well suited to a first read to build geometric intuition.
- Mas-Colell, Whinston and Green (1995, Microeconomic Theory), Chapter 3. The authoritative reference for duality and demand theory, and the source of most theorems in Sections 3.1 to 3.6 and 4.1 of this chapter; its Section 3.H on integrability and Section 3.I on welfare evaluation repay close study.
- Mas-Colell, Whinston and Green (1995), Chapter 5. Producer theory, read against the duality of Chapter 3 to see exactly where the "word-for-word translation" of Section 4.2 does its translating.
- Mas-Colell, Whinston and Green (1995), Chapter 6. Choice under uncertainty, with rigorous statements of the expected utility theorem, the Arrow-Pratt measures, and stochastic dominance.

Further reading:

- Willig (1976), Hausman (1981), and Vartia (1983). When consumer surplus is good enough, and how to recover CV exactly from demand; these correspond to Theorems 11 and 12.
- Hausman and Newey (2016). Average welfare with heterogeneous consumers: cross-sectional data do not point identify average exact consumer surplus, which can instead be bounded using bounds on income effects.
- Hausman (1996). The methodological origin of valuing new goods, worth reading for its handling of the virtual-price extrapolation and the controversy around it.
- Brynjolfsson, Hu and Smith (2003) and Brynjolfsson, Collis and Eggers (2019). Representative applications of this chapter's welfare measures to online book catalogs and to zero-price digital products, essential reading for information systems research.
- Samuelson (1938) and Houthakker (1950). The founding works of revealed preference, proposing WARP and SARP respectively.
- Afriat (1967) and Varian (1982). The revealed-preference test on finite data, that is, Theorem 8; the efficiency index of Afriat (1972) and Varian (1990) and the power benchmark of Bronars (1987) are two companion tools.
- Gorman (1953) and Gorman (1961). The condition under which aggregate demand depends only on aggregate wealth, and the explicit Gorman polar form.
- Deaton and Muellbauer (1980). The original paper on the almost ideal demand system.
- Diamond and McFadden (1974). The methodological foundation of welfare analysis with the expenditure function.
- Hall (1988) and De Loecker and Warzynski (2012). Reading markup ratios off cost minimization; the sources of Theorem 14.
- Rothschild and Stiglitz (1970). A precise definition of "more risky" through mean-preserving spreads; the source of Theorem 17 (ii).
- Rabin (2000) and Kahneman and Tversky (1979). The trouble expected utility has with small risks, and prospect theory; these correspond to Theorem 18 and Section 7.
- Savage (1954). Subjective expected utility: probabilities themselves are revealed by preferences, the starting point for discussing beliefs and the value of information.
:::

::: {.apa-refs}
- Afriat, S. N. (1967). The construction of utility functions from expenditure data. *International Economic Review, 8*(1), 67-77. https://doi.org/10.2307/2525382
- Afriat, S. N. (1972). Efficiency estimation of production functions. *International Economic Review, 13*(3), 568-598. https://doi.org/10.2307/2525845
- Bronars, S. G. (1987). The power of nonparametric tests of preference maximization. *Econometrica, 55*(3), 693-698. https://doi.org/10.2307/1913608
- Brynjolfsson, E., Collis, A., & Eggers, F. (2019). Using massive online choice experiments to measure changes in well-being. *Proceedings of the National Academy of Sciences, 116*(15), 7250-7255. https://doi.org/10.1073/pnas.1815663116
- Brynjolfsson, E., Hu, Y. J., & Smith, M. D. (2003). Consumer surplus in the digital economy: Estimating the value of increased product variety at online booksellers. *Management Science, 49*(11), 1580-1596. https://doi.org/10.1287/mnsc.49.11.1580.20580
- Deaton, A., & Muellbauer, J. (1980). An almost ideal demand system. *American Economic Review, 70*(3), 312-326.
- De Loecker, J., & Warzynski, F. (2012). Markups and firm-level export status. *American Economic Review, 102*(6), 2437-2471. https://doi.org/10.1257/aer.102.6.2437
- Diamond, P. A., & McFadden, D. L. (1974). Some uses of the expenditure function in public finance. *Journal of Public Economics, 3*(1), 3-21. https://doi.org/10.1016/0047-2727(74)90020-6
- Gorman, W. M. (1953). Community preference fields. *Econometrica, 21*(1), 63-80. https://doi.org/10.2307/1906943
- Gorman, W. M. (1961). On a class of preference fields. *Metroeconomica, 13*(2), 53-56. https://doi.org/10.1111/j.1467-999x.1961.tb00819.x
- Hall, R. E. (1988). The relation between price and marginal cost in U.S. industry. *Journal of Political Economy, 96*(5), 921-947. https://doi.org/10.1086/261570
- Hausman, J. A. (1981). Exact consumer's surplus and deadweight loss. *American Economic Review, 71*(4), 662-676.
- Hausman, J. A. (1996). Valuation of new goods under perfect and imperfect competition. In T. F. Bresnahan & R. J. Gordon (Eds.), *The economics of new goods* (pp. 207-248). University of Chicago Press.
- Hausman, J. A., & Newey, W. K. (1995). Nonparametric estimation of exact consumers surplus and deadweight loss. *Econometrica, 63*(6), 1445-1476. https://doi.org/10.2307/2171777
- Hausman, J. A., & Newey, W. K. (2016). Individual heterogeneity and average welfare. *Econometrica, 84*(3), 1225-1248. https://doi.org/10.3982/ecta11899
- Houthakker, H. S. (1950). Revealed preference and the utility function. *Economica, 17*(66), 159-174. https://doi.org/10.2307/2549382
- Kahneman, D., & Tversky, A. (1979). Prospect theory: An analysis of decision under risk. *Econometrica, 47*(2), 263-291. https://doi.org/10.2307/1914185
- Mas-Colell, A., Whinston, M. D., & Green, J. R. (1995). *Microeconomic theory*. Oxford University Press.
- Miller, N. H. (2006). *Notes on microeconomic theory* [Unpublished lecture notes]. Harvard University.
- Rabin, M. (2000). Risk aversion and expected-utility theory: A calibration theorem. *Econometrica, 68*(5), 1281-1292. https://doi.org/10.1111/1468-0262.00158
- Rothschild, M., & Stiglitz, J. E. (1970). Increasing risk: I. A definition. *Journal of Economic Theory, 2*(3), 225-243. https://doi.org/10.1016/0022-0531(70)90038-4
- Samuelson, P. A. (1938). A note on the pure theory of consumer's behaviour. *Economica, 5*(17), 61-71. https://doi.org/10.2307/2548836
- Savage, L. J. (1954). *The foundations of statistics*. Wiley.
- Varian, H. R. (1982). The nonparametric approach to demand analysis. *Econometrica, 50*(4), 945-973. https://doi.org/10.2307/1912771
- Varian, H. R. (1990). Goodness-of-fit in optimizing models. *Journal of Econometrics, 46*(1-2), 125-140. https://doi.org/10.1016/0304-4076(90)90051-t
- Vartia, Y. O. (1983). Efficient methods of measuring welfare change and compensated income in terms of ordinary demand functions. *Econometrica, 51*(1), 79-98. https://doi.org/10.2307/1912249
- von Neumann, J., & Morgenstern, O. (1947). *Theory of games and economic behavior* (2nd ed.). Princeton University Press.
- Willig, R. D. (1976). Consumer's surplus without apology. *American Economic Review, 66*(4), 589-597.
:::
