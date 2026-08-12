---
title: "Game Theory"
subtitle: "Strategic Interaction, Equilibrium, and Credible Threats"
seriesline: "Foundations of Information Systems Economics · Chapter 2"
series: "Foundations of Information Systems Economics: Empirical Methods for Information Systems Economics"
sidebartitle: "Foundations of Information Systems Economics"
sidebarsub: "Theory · Methods · Applications"
chapterlabel: "Chapter 2 · Game Theory"
---

## Introduction

Both platforms know that holding prices high earns each of them 10, while a price war leaves each with only 6. So the most natural forecast in the boardroom is: since management on both sides can do arithmetic, of course they will pick 10. Unfortunately, when each firm does its own calculation, it finds that no matter what the rival does, cutting price is better for itself. Two rational decision makers thus walk together toward an outcome neither of them likes. What is missing here is not the arithmetic, it is putting the rival's arithmetic into your own decision.

Once payoffs depend on how others respond, "which outcome is best" no longer coincides with "which outcome can be sustained." Game theory writes the players, actions, information, and payoffs as a set of primitives, then uses different equilibrium concepts to look for self-fulfilling predictions. This chapter follows that logic in order, working through dominance, Nash equilibrium, incomplete information, dynamic games, and credible threats, and then extends to correlated signals, bargaining, coordination, and strategic complementarity. Repeated interaction adds a twist: cooperation that cannot be sustained in a one-shot game may become an equilibrium because of future punishment, but the threat of "never cooperating again" must actually be worth carrying out.

Aster and Boreal will keep changing the stage on which they meet, across pricing, compatibility standards, order of moves, and private information. The same handful of small games will show why good intentions cannot substitute for incentives, why moving first is sometimes an advantage and sometimes a burden, and why, when equilibrium is not unique, empirical work needs to say something extra about the selection mechanism. The later models of pricing, entry, platform competition, and algorithmic coordination all come down to the same question: when everyone is responding to everyone else, which prediction will not be overturned by the players' own best responses?

## 1 A prediction taken for granted

::: {.case}
Two digital platforms, Aster and Boreal, compete in the same market. To shrink the strategic interaction down to something we can do by hand, we look only at their pricing stance, each choosing between two options: hold a premium (a high price, call it High) or launch an undercut (a low price, call it Low). The monthly profits under the four combinations (in millions, the first number for Aster, the second for Boreal) are in the table below. When both hold high prices they each earn 10, the best joint outcome they can reach; whoever cuts price unilaterally uses the low price to grab a large share of the other's users, so the cutter earns 14 and the firm that stays high is left with 2; when both cut, they fall into a price war and each earns 6.

| Aster \ Boreal | High | Low |
|---|---|---|
| High | 10, 10 | 2, 14 |
| Low | 14, 2 | 6, 6 |
:::

An analyst forecasts as follows: both firms are smart, both can see that when they each hold high prices they earn 10 apiece, 20 combined, which is best, and there is no reason to hand good profit over to a price war, so the market will settle at both firms high, with a predicted monthly profit of 10 each and an industry total of 20.

This prediction sounds perfectly reasonable, yet it is almost surely wrong. The trouble is that it assumes the two firms will hold back for the common good, but when each firm does its own calculation, holding back is not its best move. Take Aster's point of view. If Boreal is high, Aster can raise its payoff from 10 to 14 by cutting, so it should cut; if Boreal is low, Aster cutting too goes from 2 back to 6, so again it should cut. No matter what Boreal does, cutting is always better for Aster, cutting is its strictly dominant strategy. Boreal faces the symmetric tradeoff, and cutting dominates for it as well. So both cut, the market settles at the mutually damaging (Low, Low), each earns 6, and the industry total is 12. That forecast of 20 is fully two thirds above the actual 12, and it is wrong in direction: the high-price cooperation it predicts is exactly the outcome that cannot survive in this game.

The right question is not "which outcome is good for both sides," but "given that the rival will respond in the way that is best for itself, which outcome can sustain itself." Here the two questions point in opposite directions, and that is precisely why game theory exists. Section 3 will hammer the intuition of "self-sustaining" into the precise definition of a Nash equilibrium, and demonstrate it on this same batch of platform games; Section 4 will bring a turning point: if Aster and Boreal do not meet just once but month after month, that cooperative outcome of 20 can be brought back to life under a clear condition.

## 2 The primitives of a game and two representations

To analyze a strategic situation, we first have to write it down completely. A noncooperative game is pinned down by four sets of primitives: the set of players, the actions available to each player, the information each player holds when making a choice, and the payoff each player receives under each combination of actions. Payoffs are already numbers converted from final outcomes through von Neumann-Morgenstern utility (the expected utility of Chapter 1 comes in handy here), so under randomness players maximize expected payoff. The table in Section 1 is exactly such a game: the players are Aster and Boreal, the actions are each {High, Low}, the two firms decide simultaneously without knowing the other's current choice, and the payoffs are in the table.

Here we need to distinguish two easily confused concepts: an action (a single concrete choice) and a strategy. A strategy is a complete plan of action, specifying what a player does at every decision point where it might be called upon. In a simultaneous one-shot game, strategy equals action, the two are not separate; but in a dynamic game with order of moves and information updating, a strategy is much richer, it must spell out "what I would do in every situation I could reach," even situations that, under the equilibrium, never actually occur. This distinction is central when we discuss credibility in Section 3.4, since whether a threat is credible is exactly the question of whether the action that carries it out is the strategy's best move at the corresponding decision point.

The same game has two standard representations, each convenient for a class of analysis. The normal form (also called the strategic form) compresses the game into a payoff table: each player picks a strategy and the table gives the outcome directly, and the table in Section 1 is a normal form. It erases the details of timing and information, and suits the analysis of simultaneous games. The extensive form draws the game as a tree: nodes are decision points, the branches growing out of a node are actions, information is captured by information sets (grouping together nodes a player cannot tell apart), and payoffs hang on the leaves. The extensive form keeps the full structure of who moves first and who sees what, and is indispensable for analyzing dynamic games, as Section 3.4 will use. The two representations can be converted into each other, but the extensive form carries more information, and compressing it into the normal form loses the timing. This information loss is exactly why solution concepts based on the normal form miss the credibility problem in dynamic games.

The target quantity of this chapter is a prediction: given the primitives of a game, which strategy combination (or combinations) rational players will land on. This prediction is delivered by a solution concept, and there is more than one solution concept, corresponding to different strengths of the requirement of "rationality." The next section works its way up this ladder of strength, level by level. It is the most fine-grained part of the whole chapter's analysis, and the foundation for understanding every later structured game-theoretic model.

## 3 The core structure: a ladder of solution concepts

This section starts from the weakest requirement of rationality and strengthens it level by level, and at each level we see what additional wrong predictions it rules out in the platform games. Let us fix notation first: player $i$'s strategy is written $s_i$, the strategy combination of everyone else is $s_{-i}$, and the payoff function is $u_i(s_i, s_{-i})$.

### 3.1 Dominance, iterated deletion, and rationalizability

The minimal requirement of rationality is: do not choose a strategy that is worse no matter what others do. If for player $i$ there is a strategy $s_i'$ such that $u_i(s_i', s_{-i}) > u_i(s_i, s_{-i})$ for all $s_{-i}$, we say $s_i$ is strictly dominated by $s_i'$, and a rational player will never use a strictly dominated strategy. In the pricing game of Section 1, High is strictly dominated by Low for both firms, so a single step deletes down to (Low, Low), and dominance alone already delivers a unique prediction for this game, without even invoking Nash equilibrium.

Most games are not this clean, but we can delete iteratively. After deleting every player's dominated strategies, some strategies that were not dominated before may become dominated in the smaller game that remains, so we delete again, and so on. This is called iterated deletion of strictly dominated strategies. The assumptions behind it ratchet up step by step: the first round needs only that everyone is rational, the second needs "I know you are rational," the third needs "I know that you know that I am rational," and deleting all the way to the end uses common knowledge of rationality (all levels of "I know that you know ..." hold). The set of strategies left after deletion is the set of rationalizable strategies, the precise answer to "how much can be ruled out by common knowledge of rationality alone." The key point is that rationalizability usually does not delete cleanly, and a large set of strategies survives, which shows that common knowledge of rationality by itself is often not enough to pin down a unique prediction. We need something stronger, and that is Nash equilibrium.

### 3.2 Nash equilibrium: mutual best responses

Nash equilibrium writes the intuition of "self-sustaining" as a precise condition: a strategy combination $s^* = (s_i^*, s_{-i}^*)$ is a Nash equilibrium if each player's strategy is a best response to everyone else's,

$$u_i(s_i^*, s_{-i}^*) \geq u_i(s_i, s_{-i}^*) \quad \text{for all } s_i,\ \text{all } i$$

In plain words, standing at the equilibrium, no side can make itself better off by unilaterally changing strategy. That is exactly what "self-sustaining" means: an equilibrium is a strategy combination from which no one has an incentive to deviate. It is stronger than rationalizability because it additionally requires the players' beliefs to confirm one another, each person's expectation of others' behavior is exactly equal to others' actual equilibrium behavior, a version of rational expectations brought down to the level of strategies. In the pricing game, (Low, Low) is the unique Nash equilibrium, both firms earn 6, and either one that switches back to High unilaterally only falls to 2.

Nash equilibrium is not always unique, nor does it always exist in pure strategies, and the compatibility game exposes both points. Aster and Boreal each choose a technology standard, S1 or S2, and only when both use the same standard can they interconnect and each reap a network payoff, while going their separate ways leaves them unable to connect, with zero payoff. Aster prefers S1, Boreal prefers S2, but both would rather accommodate the other than fail to connect. The payoffs are in the table below.

| Aster \ Boreal | S1 | S2 |
|---|---|---|
| S1 | 3, 1 | 0, 0 |
| S2 | 0, 0 | 1, 3 |

This game has two pure-strategy Nash equilibria: (S1, S1) and (S2, S2), each a coordination onto the same standard, and no one wants to jump alone to the other standard and turn its payoff into zero. The trouble is that the two equilibria distribute differently, Aster prefers the former and Boreal the latter, and when they decide simultaneously neither knows which one to bet on. This is the heart of the coordination problem: the existence of equilibria does not guarantee that coordination is achieved.

Besides the two pure-strategy equilibria, there is also a mixed strategy Nash equilibrium, in which players randomize their choices with probabilities. Let Aster choose S1 with probability $p$ and Boreal choose S1 with probability $q$. A mixed strategy equilibrium requires each side's randomization to make the other exactly indifferent between its two pure strategies (otherwise the other would tilt toward the pure strategy with the higher payoff and stop mixing). Making Boreal indifferent: choosing S1 gives it $p\cdot 1 + (1-p)\cdot 0 = p$, choosing S2 gives it $p\cdot 0 + (1-p)\cdot 3 = 3(1-p)$, and setting the two equal gives $p = 3(1-p)$, that is $p = 3/4$. Symmetrically, making Aster indifferent: choosing S1 gives $3q$, choosing S2 gives $1-q$, and $3q = 1-q$ yields $q = 1/4$. So the mixed equilibrium has Aster choose its preferred S1 with probability $3/4$ and Boreal choose S1 with probability $1/4$ (that is, choose its preferred S2 with probability $3/4$).

The outcome of this mixed equilibrium is especially worth noticing. The probability that the two firms actually coordinate onto the same standard is $p q + (1-p)(1-q) = \tfrac{3}{4}\cdot\tfrac{1}{4} + \tfrac{1}{4}\cdot\tfrac{3}{4} = \tfrac{3}{8}$, meaning that with probability $5/8$ the two firms choose different standards and both get zero. Each side's expected payoff is only $3/4$, far below the 3 or 1 that the successful coordinator gets in either pure-strategy equilibrium. This is not a case of anyone making a mistake, but the inevitable product of rational randomization under simultaneous choice: each side is correctly guarding against the possibility that the other picks a different standard, and the result of that guarding is a high probability of coordination failure. This example also foreshadows why platform standards battles so often stall or split, and why "who moves first" matters so much, the latter being the subject of Section 3.4.

::: {.theorem}
(Nash 1950b existence) Every finite game (finitely many players, each with finitely many strategies) has at least one Nash equilibrium, if mixed strategies are allowed.
:::

The existence theorem guarantees that "no pure-strategy equilibrium can be found" never means the game has no equilibrium, a mixed strategy equilibrium is surely there, which is why mixed strategies are not a dispensable technical patch but a necessity for making the equilibrium concept complete. The proof uses a fixed-point theorem (Kakutani or Brouwer), stitching everyone's best-response map into a continuous correspondence from strategy combinations to strategy combinations, whose fixed point is the equilibrium. For details, see the appendix of MWG Chapter 8, which we do not develop here.

### 3.3 Incomplete information and Bayesian Nash equilibrium

So far both platforms fully understand each other's payoffs. In reality this is often not the case: a platform's cost, technological readiness, and strategic resolve are often private information that only it knows. Harsanyi's treatment is to introduce a type for each player, encoding its private information, and to assume that the prior distribution of types is common knowledge, on which basis players make Bayesian inferences about the rival's type. With this, a game of incomplete information is turned into a game of imperfect information over types, and the equilibrium concept is correspondingly generalized to Bayesian Nash equilibrium: each type of each player best-responds given the distribution of the rival's types (and the rival's type-contingent strategies).

Take a platform version of the example. Aster is uncertain about Boreal's cost structure: with probability $\beta$ Boreal is the Normal type, whose cost makes Low pricing dominant; with probability $1-\beta$ it is the Committed type (say it has just signed an exclusive for high-value content and monetizes through premium positioning), whose cost makes High dominant. So Boreal's equilibrium strategy is obvious: the Normal type prices Low, the Committed type prices High. The difficulty is on Aster, which must place its bet without knowing the other's type, best-responding to Boreal's type mixture. Suppose Aster's payoff structure is such that whether it wants to differentiate from or align with the rival's pricing stance depends on the specific numbers, and here we take: facing High, Aster gets 10 from High and 7 from Low; facing Low, it gets 3 from High and 6 from Low. Given that Boreal is Low with probability $\beta$ and High with probability $1-\beta$, Aster's expected payoffs are

$$\mathbb{E}[u_{\text{High}}] = \beta\cdot 3 + (1-\beta)\cdot 10 = 10 - 7\beta, \qquad \mathbb{E}[u_{\text{Low}}] = \beta\cdot 6 + (1-\beta)\cdot 7 = 7 - \beta$$

Aster chooses High if and only if $10 - 7\beta \geq 7 - \beta$, that is $\beta \leq 1/2$. The intuition is clear: the more likely Boreal is the high-pricing Committed type (the smaller $\beta$), the more Aster should price High to align; the more likely Boreal is the cutthroat Normal type (the larger $\beta$), the more Aster should retreat to Low. Taking $\beta = 2/5$ (Committed types more common), Aster gets $10 - 7\cdot\tfrac{2}{5} = 7.2$ from High and $7 - \tfrac{2}{5} = 6.6$ from Low, so the Bayesian Nash equilibrium is: Aster prices High, Boreal's Normal type prices Low and Committed type prices High. The key to a BNE is not the complexity of the algebra but a shift in thinking: Aster does not best-respond to some specific rival, it best-responds to a distribution of types, and its equilibrium strategy is optimal with respect to the expectation under uncertainty. This idea of "best-responding to a belief" is the direct forerunner of auctions and mechanism design later on.

### 3.4 Dynamic games, credibility, and subgame perfection

The simultaneous-choice framework misses a dimension that is extremely important in platform competition: order. Who announces a standard first, who cuts price first, who enters first, often determines the outcome, and the normal form flattens the timing so that this layer cannot be seen. Dynamic games are represented with the extensive-form tree, and their central new problem is credibility: a player can threaten "if you dare enter I will start a price war," but if, when it actually comes to that, a price war is a loss for it too, the threat is empty talk, and a rational rival will not be scared off. Nash equilibrium does not block empty threats, because it only requires strategies to be mutually best on the equilibrium path, without asking whether those threats are worth carrying out when they are actually triggered.

Turn the compatibility game into a sequential one and the problem shows up at once. Now let Aster choose a standard first and announce it, and Boreal choose after seeing it.

![Extensive form of the sequential compatibility game: Aster chooses S1 or S2 first, Boreal chooses after observing; backward induction (bold branches) gives the SPE: Aster chooses S1, Boreal follows, outcome (3, 1).](assets/fig/fig_02_gametree.svg)

The standard way to solve a dynamic game is backward induction: work back from the ends of the tree, at each decision point assuming the player there makes the currently optimal choice, replace the whole subtree with that choice, and roll back level by level to the root. Look first at Boreal's response. If Aster chose S1, Boreal gets 1 from following with S1 and 0 from jumping to S2, so following is better; if Aster chose S2, Boreal gets 3 from following with S2 and 0 from jumping to S1, so it also follows. So Boreal's best response is always to match Aster. Aster foresees this, so it knows that choosing S1 gets it 3 (Boreal will follow) while choosing S2 gets only 1, so it chooses S1. The backward-induction outcome is (S1, S1), payoffs (3, 1): Aster locks in its preferred standard by moving first. The headache of coordination failure and equilibrium multiplicity under simultaneous choice is dissolved in one stroke by "who moves first," and this is the precise source of the first-mover advantage.

The result of backward induction corresponds to a solution concept stronger than Nash equilibrium, subgame perfect equilibrium. It requires the strategy combination to be a Nash equilibrium not just in the whole game but also when restricted to every subgame (the complete continuation game starting from any single decision node). This requirement of "equilibrium everywhere" is exactly what rules out empty threats: if a threat is not optimal in the very subgame where it would have to be carried out, then the strategy combination it belongs to does not constitute a Nash equilibrium in that subgame, and so is thrown out by SPE. The sequential compatibility game in fact has other Nash equilibria, for example Boreal declaring "no matter what you choose I will choose S2," which can scare Aster into also choosing S2 and gives the Nash equilibrium (S2, S2), but this threat is not credible: in the subgame where Aster actually chose S1, Boreal gets 1 from following with S1, better than 0 from holding out with S2, so it will not actually carry it out. SPE eliminates all such Nash equilibria propped up by empty threats and leaves only the backward-induction (S1, S1).

::: {.intuition}
Credibility is the soul of dynamic games, and its test can be summed up in one sentence: looking only at what everyone does on the equilibrium path is not enough to judge whether a strategy is rational, we also have to look at what it has you do in those "should not happen" surprise situations, and that arrangement must be one you would actually be willing to carry out if you did get there. Nash equilibrium checks only the path, SPE checks every corner. Commitments, deterrence, and concessions in platforms sometimes work and sometimes are bluffs, and the dividing line is entirely here: whether you can make good on your words when the moment comes. The repeated games of Section 4 push this point to its limit, where cooperation is sustained precisely by a threat that "if you dare betray me I will cut price forever," and it is credible precisely because once the other side betrays, cutting price along with it is genuinely optimal for you too.
:::

### 3.5 Refinements under imperfect information

When a dynamic game contains information sets (a player, when it is her turn, does not fully know what happened before, for example she did not see the rival's private choice), there may be too few subgames for SPE to get any traction, because a legitimate subgame cannot be cut out from the middle of an information set. This calls for finer refinements, which fold "the beliefs a player holds at an information set" into the equilibrium as part of it. Sequential equilibrium (Kreps and Wilson 1982) requires strategies and beliefs to be consistent: beliefs are derived from equilibrium strategies via Bayes' rule, and strategies are in turn sequentially optimal given beliefs, the two supporting each other. Trembling-hand perfection (Selten 1975) refines from a different angle, requiring the equilibrium to be robust to a tiny "trembling hand" (choosing any strategy by mistake with vanishingly small probability), so that any strategy that is optimal only when the rival never errs is eliminated. The common motivation of these refinements is to give a reasonable constraint on off-path beliefs (what to believe after a deviation from the equilibrium path). We do not develop their technical details in this chapter, and only need to know that they are all extensions of the credibility idea under imperfect information, and later in this series, when we cover signaling games, we will make concrete use of off-path belief refinements.

### 3.6 Existence: how a fixed point guarantees equilibrium

The existence theorem in Section 3.2 gave only the conclusion, and here we fill in why it holds, because understanding the structure of the proof clarifies a common confusion: why mixed strategies are not a dispensable technical patch. View each player's best response as a map: when $i$ faces others' strategies $s_{-i}$, its best-response set is $\mathrm{BR}_i(s_{-i}) = \arg\max_{s_i} u_i(s_i, s_{-i})$. Stacking everyone's best responses together gives a correspondence from strategy combinations to strategy combinations, $\mathrm{BR}(s) = \big(\mathrm{BR}_1(s_{-1}), \dots, \mathrm{BR}_n(s_{-n})\big)$. The definition of Nash equilibrium now has an extremely compact equivalent statement: $s^*$ is a Nash equilibrium if and only if it is a fixed point of this best-response correspondence, $s^* \in \mathrm{BR}(s^*)$, that is, "everyone's best responses to each other lead right back to everyone themselves." Finding an equilibrium thus becomes finding a fixed point.

Fixed-point theorems are made precisely for this self-referential structure. The Brouwer fixed point theorem says: a continuous self-map on a compact convex set always has a fixed point. Its two-dimensional intuition is that if you crumple a map and lay it flat back on top of the original map, some point is always pressing on its original position. But best responses are usually not single-valued functions but set-valued correspondences (at indifference there are several best responses), which calls for the Kakutani fixed point theorem, a generalization of Brouwer to set-valued correspondences, requiring the correspondence to be nonempty, convex-valued, and closed-graphed (upper hemicontinuous). Matching these three conditions to the game explains exactly why mixed strategies are necessary: only by expanding the strategy space from the finite pure strategies to the simplex of mixed strategies does the domain become compact and convex, and expected payoff is linear in one's own mixing probabilities, and this linearity in turn guarantees that best responses are convex-valued. Without the step of mixing, the domain is a discrete set of points, neither convex nor continuous, and the fixed-point theorem has nowhere to land, so a pure-strategy equilibrium can indeed fail to exist (matching pennies has no pure-strategy equilibrium), but once mixing is allowed, Kakutani immediately guarantees a fixed point, that is, the existence of an equilibrium.

::: {.intuition}
Put the mixed equilibrium of the compatibility game into a coordinate system and this abstraction becomes concrete at once. Let the horizontal axis be the probability $p$ that Aster chooses S1, the vertical axis the probability $q$ that Boreal chooses S1, and draw a best-response curve for each: given the other's probability, with what probability should I choose S1. Both curves are step-shaped with a kink, and they must cross, and the crossing point $(p, q) = (3/4, 1/4)$ is exactly the mixed equilibrium solved in Section 3.2. The existence theorem is nothing but the general version of this geometric fact: best-response curves cannot dodge each other in a compact convex strategy space, they always have to cross somewhere, and the crossing point is the equilibrium. Finding an equilibrium and finding the intersection of curves are one and the same thing.
:::

### 3.7 Correlated equilibrium: correlated signals and algorithmic coordination

Nash equilibrium assumes by default that everyone randomizes independently, without reference to each other. But in reality players can often observe some common external signal, a traffic light, an industry-wide public price index, a recommendation algorithm everyone uses, and then coordinate their actions on it. Folding this "act on the signal" behavior into equilibrium gives a concept broader than Nash, correlated equilibrium (Aumann 1974). Imagine a mediator draws an action combination according to some joint distribution $\mu$, then privately tells each person only his own component as a "recommendation." This $\mu$ is a correlated equilibrium if and only if each person, after receiving the recommendation $s_i$ and using it to condition his inference about others' actions, finds following it optimal:

$$\sum_{s_{-i}} \mu(s_{-i}\mid s_i)\,\big[\,u_i(s_i, s_{-i}) - u_i(s_i', s_{-i})\,\big] \geq 0 \quad \text{for all } s_i',\ \text{all } i$$

That is, no one has an incentive to deviate from the recommendation. Every Nash equilibrium is a correlated equilibrium (just make the signals mutually independent), and public randomization over several Nash equilibria is also a correlated equilibrium, so the set of correlated equilibria contains at least the convex hull of Nash equilibria, and is often strictly larger.

Return to the compatibility game. Let a public fair coin serve as the coordinating device: heads recommends both firms choose S1, tails recommends both choose S2, each with probability one half, and both firms can see the toss. On receiving "S1," Aster knows Boreal is also recommended S1 and will follow, and following with S1 gets 3, better than 0 from jumping to S2, so following is optimal; the same for Boreal. So this is a correlated equilibrium, and each firm's expected payoff is $\tfrac{1}{2}\cdot 3 + \tfrac{1}{2}\cdot 1 = 2$, far above the pitiful $3/4$ of the mixed equilibrium, and fairly centered between the two pure equilibria the two firms prefer, resolving in one stroke the coordination failure of Section 3.2. If the public signal is replaced by a private correlated signal, a correlated equilibrium can even reach payoffs that no convex combination of Nash equilibria can attain (Aumann gave a famous example using the game of Chicken), which we do not develop here. Correlated equilibrium is especially apt for information systems: a pricing or recommendation algorithm used in common by many parties is itself a correlating device, and it can coordinate the market onto some outcome without any explicit collusion. Questions like "the algorithm is quietly coordinating" need exactly the language of correlated equilibrium to characterize, and later in this series, when we cover algorithms and markets, we will return to it.

### 3.8 Bargaining: Nash bargaining and Rubinstein alternating offers

The games above are mostly zero-sum confrontations, but a great deal of interaction in the platform economy is of another kind: two sides cooperating can create a surplus, and the only question is how to divide it. A platform negotiating a revenue split with a content supplier, two firms negotiating equity to form a joint venture, a deal makes everyone happy while a breakdown sends everyone back home, and how that surplus finally gets cut is the bargaining problem. It has two seemingly unrelated solutions that reach the same conclusion, and understanding this is key to the platform bargaining analysis on the supply side later on.

The first is Nash's (1950a) axiomatic solution, treating bargaining as a cooperative game. Given the set of feasible utility combinations $U$ (convex, compact) and the disagreement point $d = (d_1, d_2)$ that each side gets on breakdown, Nash proved that the solution satisfying four plain axioms (Pareto efficiency, symmetry, invariance to affine transformations of utility, and independence of irrelevant alternatives) exists and is unique, and it is the point that maximizes the Nash product:

$$\max_{(u_1, u_2)\in U}\; (u_1 - d_1)(u_2 - d_2)$$

When both sides are risk-neutral and dividing a sum of money of size $V$, this solution has a plain-language reading: each side first takes back the $d_i$ it would get on breakdown, then the remaining net surplus $V - d_1 - d_2$ is split evenly. Section 6.6 will compute the $50/50$ when the disagreement point is $(0,0)$, and the $60/40$ when one side has an outside option $d = (20, 0)$. Replacing the symmetry axiom with a weighted version gives generalized Nash bargaining, where one side gets $d_1 + \alpha(V - d_1 - d_2)$, and the weight $\alpha$ measures its relative bargaining power. This is exactly the origin of the bargaining-power parameter in the platform-supplier bargaining (Nash-in-Nash) later in this series.

The second is Rubinstein's (1982) strategic solution, treating bargaining as a dynamic game and solving it with the SPE machinery of Section 3.4. Two players take turns proposing a division of a surplus of size 1, one proposes and the other chooses to accept or reject, rejection moves to the next round where the other proposes, but each round of delay shrinks the surplus by a discount factor $\delta$. In infinite-horizon alternating offers, there seem to be infinitely many paths of delay and counteroffer, yet the SPE is unique. Solve it using stationarity: let the proposer's equilibrium share be $x^*$. To make the responder just willing to accept, the proposer must give the responder its continuation value after rejecting, and after rejecting the responder becomes the proposer next round and can get $x^*$ but discounted one period, so the continuation value is $\delta x^*$, and thus the proposer offers $\delta x^*$ and keeps $1 - \delta x^*$. Stationarity requires that what it keeps is exactly $x^*$:

$$x^* = 1 - \delta x^* \;\Longrightarrow\; x^* = \frac{1}{1+\delta}$$

So the unique SPE has the proposer take $1/(1+\delta)$ and the responder take $\delta/(1+\delta)$, and the first mover has the advantage, an advantage that comes from the responder taking the discount hit for delaying.

The relationship between the two solutions is one of the most beautiful conclusions in game theory, called the Nash program: when frictions vanish, that is $\delta \to 1$ (both sides extremely patient, delay nearly costless), Rubinstein's SPE share converges to $1/2$ against $1/2$, exactly the solution of symmetric Nash bargaining; while different discount factors for the two sides correspond to weighted generalized Nash bargaining. The strategic solution provides microfoundations for the axiomatic solution, and this is why, when we later use "Nash bargaining" as shorthand for the outcome of a negotiation, there is in fact a flesh-and-blood alternating-offer process standing behind it. Section 6.6 computes the numbers of both solutions side by side and checks that they match.

### 3.9 Global games and supermodular games: coordination, tipping, and equilibrium selection

Platform adoption and standards battles are a special class of games: the benefit of adopting a platform rises with the number of adopters, which is called strategic complementarities. The old ailment of this class of coordination games is equilibrium multiplicity, everyone adopting and no one adopting are often both equilibria, and theory cannot say whether a new platform will tip (be adopted by all in a cascade) or go ignored. Two tools treat this multiplicity.

The first is supermodular games (Milgrom and Roberts 1990). As long as a game has strategic complementarities, it has a largest Nash equilibrium and a smallest Nash equilibrium, both in pure strategies, and comparative statics are monotone: an exogenous change favorable to adoption pushes both the largest and the smallest equilibrium up at the same time. This gives coordination games a clean skeleton, but other equilibria may still be sandwiched between the largest and the smallest, so multiplicity is not eradicated.

The second goes to the root, and is global games (Carlsson and van Damme 1993; Morris and Shin 1998). It mixes a tiny bit of incomplete information into complete information: no one knows the fundamental $\theta$ exactly anymore, each seeing only a private signal with tiny noise. With this touch, common knowledge is broken, and the equilibrium surprisingly converges to a unique threshold strategy, that is, adopt if the signal is above a critical value and not otherwise, so that the fundamental critical value $\theta^*$ for tipping is uniquely pinned down. What lets it select a unique equilibrium is a property called the Laplacian belief: the marginal player standing at the critical signal has, about "what fraction of people will adopt," a posterior that is exactly uniform on $[0, 1]$. Take a setup we can compute by hand: the payoff to adopting is $\theta + \ell - 1$, where $\ell$ is the fraction of adopters, and the payoff to not adopting is $0$. Under complete information, as long as $\theta \in (0, 1)$, both all adopting (payoff $\theta > 0$ when $\ell = 1$) and none adopting (payoff $\theta - 1 < 0$ when $\ell = 0$) are equilibria, so multiplicity. In the global game, the marginal player expects $\mathbb{E}[\ell] = 1/2$, and the indifference condition $\theta^* + \tfrac{1}{2} - 1 = 0$ immediately gives the unique critical value $\theta^* = 1/2$: if the fundamental is above $1/2$ the platform is adopted by all in a cascade, and if below it sinks into silence, and the equilibrium selected is exactly the risk-dominant one. Section 6.8 reviews this critical value. Global games are the theoretical pillar for understanding coordination phenomena like platform tipping, bank runs, and currency attacks, and later in this series, when we cover platforms and network effects, we will use it head-on.

Here is a summary: solution concepts are a ladder arranged from weak to strong requirements of rationality. Dominance and rationalizability use only common knowledge of rationality, and often cannot delete down to a unique prediction; Nash equilibrium adds mutual confirmation of beliefs, giving a self-sustaining strategy combination, but it may be multiple and may exist only in mixed strategies, and its existence is guaranteed by the fixed point of the best-response correspondence, with mixing being exactly the step that makes the fixed-point theorem land; Bayesian Nash equilibrium generalizes it to private information, players best-responding to a distribution of types; subgame perfection enters dynamic games, using "equilibrium everywhere" to eliminate incredible empty threats; and finer refinements handle beliefs under imperfect information. Beyond the ladder there are a few companion tools: correlated equilibrium folds in a common signal to enable coordination and characterizes algorithms quietly coordinating, the two solutions of bargaining (Nash's axiomatic and Rubinstein's alternating offers) reach the same division of surplus under the Nash program, and global games use a bit of incomplete information to select a unique tipping critical value for coordination games. With each level climbed and each tool added, fewer predictions survive and they are more credible, and this map is the core takeaway of the chapter.

### 3.10 Topkis, organizational complementary systems, and empirical tests

Strategic complementarity does not only occur between different firms, it also occurs among the multiple organizational choices within the same firm. Suppose a firm chooses AI adoption $a$, employee skill investment $s$, and decentralization of decision rights $d$ simultaneously, with profit $\pi(a,s,d;\theta)$. If any two choices have increasing differences, for example

$$\frac{\partial^2\pi}{\partial a\,\partial s}\geq0,
\qquad
\frac{\partial^2\pi}{\partial a\,\partial d}\geq0,$$

then raising skill investment increases the marginal return to AI adoption, and AI adoption also increases the marginal return to decentralization. If the choice set is a lattice and the profit function is supermodular, the Topkis monotonicity theorem guarantees that the smallest and largest optimal choices rise monotonically in the parameter $\theta$ that raises marginal returns. This result requires neither concavity of the profit function nor uniqueness of the optimum, and is especially suited to studying organizational practices that come in bundles.

The theory also explains why "bought AI but productivity did not rise" cannot directly falsify the value of the technology: moving only $a$ without moving $s$ and $d$, the firm may still be stuck in a low-complementarity system. Empirically there are two common types of test. A correlation test checks whether complementary practices are adopted together, but common shocks also create positive correlation; a performance test checks whether joint adoption brings superadditive gains, but observational interaction is affected by selection. The most powerful design randomizes or quasi-randomizes two practices separately, identifying with a factorial contrast

$$\Delta_{as}=\{Y(1,1)-Y(0,1)\}-\{Y(1,0)-Y(0,0)\}.$$

$\Delta_{as}>0$ is the empirical counterpart of causal complementarity. If we can only observe firms' own choices, we should at least report the adoption bundles, adjust the timing, and state what selection restriction the identification of the interaction relies on. Chapter 26 will apply this theory to AI, workflow redesign, and human capital.

## 4 Repeated games and the possibility of cooperation

Section 1 left a jarring conclusion: the unique equilibrium of the pricing game is the mutually damaging (Low, Low), and the joint high-price 20 cannot survive. But that is only the outcome of a one-shot encounter. Aster and Boreal are long-lived coexisting platforms, meeting month after month, and whether this long-term relationship can bring cooperation back is the question of this section, and the answer is yes, but conditionally.

The key mechanism is: in repeated interaction, betrayal today invites retaliation starting tomorrow, and as long as the future matters enough, the immediate gain from a sneak attack cannot outweigh the loss from being retaliated against later. Let us formalize this intuition. Suppose both firms discount future profits by a discount factor $\delta \in (0,1)$, and the closer $\delta$ is to 1 the more they weight the long run. Consider the grim trigger strategy: cooperate with High at the start, keep playing High as long as the other cooperates, and once you detect that the other has cut price, switch to Low forever in retaliation. If both use this strategy, whether cooperation can be sustained depends on comparing the discounted total payoffs of two paths. Keeping the deal and cooperating forever gets $R = 10$ per period, with total value $R/(1-\delta)$. Launching a sneak cut in some period gets $T = 14$ that period, but triggers the other's permanent retaliation, after which every period falls to $P = 6$, with total value $T + \delta P/(1-\delta)$. Cooperation can be sustained if and only if keeping the deal is no worse than the sneak attack:

$$\frac{R}{1-\delta} \geq T + \frac{\delta P}{1-\delta}$$

Rearranging gives the threshold on the discount factor

$$\delta \geq \frac{T - R}{T - P}$$

Plugging in the numbers of the pricing game, the threshold is $\delta \geq (14-10)/(14-6) = 4/8 = 1/2$. That is, as long as the two platforms weight the future by more than half ($\delta \geq 0.5$), jointly holding high prices is a subgame perfect equilibrium, and the 20 that "could not survive" in Section 1 survives in a long-term relationship that is patient enough. The structure of the threshold is itself an insight: the larger the immediate temptation of the sneak attack $T - R$, or the lighter the punishment of the price war $T - P$, the more patience is needed to sustain cooperation. This explains many empirical regularities of platform competition, for example why markets with high entry barriers, few players, and close mutual monitoring more easily sustain high prices, while markets where one price cut can grab most of the market (a huge $T - R$) have especially fragile cooperation.

The general form of this conclusion is the folk theorem. It says: in an infinitely repeated game, as long as players are patient enough ($\delta$ close enough to 1), any feasible payoff combination that gives every player strictly more than its minmax (the lowest value the others can jointly hold it to) can be supported as the average payoff of some subgame perfect equilibrium.

![The folk theorem payoff set: the shaded region is the convex hull of feasible payoffs, the dashed lines are the two firms' respective minmax levels, and the darker region in the upper right (feasible payoffs above both minmax lines) can all be supported as SPE by a sufficiently patient repeated game.](assets/fig/fig_02_folk.svg)

The folk theorem is a double-edged sword. On the optimistic side, it proves that long-term relationships can nurture cooperation, flipping over the inefficient outcome that is doomed in a one-shot game. On the pessimistic side, the supportable payoffs it gives are an entire region rather than a single point, meaning that the equilibria of a repeated game are severely multiple, and the theory itself cannot select which one will be realized, and anything from mutual damage to full collusion can be an equilibrium. This multiplicity is not an oversight of the theory but a faithful reflection of reality: the course of long-term competition does indeed depend heavily on expectations, conventions, and communication, things outside the model. It also poses a difficulty for empirical work, since an observed high price could be collusion or could be some other equilibrium, and equilibrium concepts alone cannot tell them apart, and this is exactly the problem that conduct tests and collusion identification later in this series will have to confront head-on.

Here is a summary: cooperation that is sentenced to death by the logic of dominance in a one-shot game can be revived when the game is infinitely repeated and players are patient enough, sustained by the credible threat that "betrayal invites retaliation," with the threshold $\delta \geq (T-R)/(T-P)$ determined by the ratio of the sneak-attack temptation to the strength of retaliation; the folk theorem pushes this to the general case, at the cost of severe equilibrium multiplicity, and the theory can say what outcomes are supportable but cannot say which one will happen.

### 4.1 Comparative statics, equilibrium selection, and empirical implications

The difficulty of comparative statics in game theory is that the equilibrium itself changes. If actions are strategic complements, a parameter that raises the marginal return to one side's action pushes its best response up, and produces a multiplier effect through the rival's response. Supermodular games under appropriate conditions guarantee that the smallest and largest equilibria move monotonically with the parameter, but they do not guarantee which equilibrium reality selects.

Suppose a platform adopts action $a_i\in\{0,1\}$, with adoption payoff $\theta_i+\gamma a_j-c$. Here $\gamma>0$ represents strategic complementarity. A fall in cost weakens the appeal of not adopting, predicting a rise in the adoption rate; near the tipping region, a small change in cost can trigger a large coordination jump. Observing an adoption jump is consistent with complementarity, but it could also come from a common demand shock. Identifying $\gamma$ requires variation that moves one side's payoff without directly moving the other side's payoff, or exploits the structure of network connections and restrictions on timing.

Repeated games give another kind of comparative statics. The constraint for grim-trigger cooperation is

$$\frac{R}{1-\delta}\geq T+\frac{\delta P}{1-\delta},\qquad \delta\geq\frac{T-R}{T-P}.$$

When interaction is more frequent, the probability the relationship continues is higher, or punishment is harsher, cooperation is easier to sustain. Empirically, price wars are more likely to appear when demand drops sharply, monitoring worsens, or the relationship nears its end. But observed price parallelism does not equal collusion, since cost co-movement and common algorithms also produce the same pattern.

| Change in theoretical parameter | Equilibrium prediction | What can be observed empirically | Extra handle needed for identification |
|---|---|---|---|
| Complementarity $\gamma$ rises | Best responses steeper, multiple equilibria and tipping more likely | Adoption clusters among connected agents, shocks have network multipliers | Exogenous network edges, one-sided payoff shifter |
| Discount factor $\delta$ rises | Cooperation incentive strengthens | More stable prices, fewer punishment periods | Exogenous relationship duration or monitoring change |
| Private-noise variance falls | Global-game action closer to common threshold | Actions more concentrated in the fundamental critical region | Randomized information precision or disclosure rule |
| First-mover commitment strengthens | SPE may deviate from simultaneous Nash | After the timing changes, follower response varies systematically | Exogenous rollout order or institutional timing |

A test of theoretical consistency should specify equilibrium selection in advance. If a parameter change moves several equilibria at once, the direction of the observed outcome also depends on the selection rule. Calibrating some Nash equilibrium to the data does not mean that equilibrium is identified; at a minimum one should report whether other equilibria give the same observable predictions.

## 5 Anchoring papers

The abstract concepts of game theory show their power only when brought down to real strategic situations. This section picks two foundational works, one giving the most general answer to how cooperation is possible in repeated games, and one bringing the compatibility game into the reality of the information-technology market, each laid out by paper, method, result, significance, and limitations, with a focus on how theoretical assumptions correspond to real mechanisms.

### 5.1 Fudenberg and Maskin (1986)

::: {.case}
Paper and methodological place: "The Folk Theorem in Repeated Games with Discounting or with Incomplete Information," Econometrica. It nailed the long-circulating but never quite rigorous folk theorem into a theorem, answering the central question of Section 4: in a discounted infinitely repeated game, which payoffs exactly can be supported by a subgame perfect equilibrium.

Method: for an infinitely repeated game, they construct the strategies needed to support a target payoff (cooperation plus credible punishment for betrayal), and prove that when the discount factor is close enough to 1, any individually rational feasible payoff can be realized as the average payoff of an SPE. Compared with earlier folk theorems that spoke only of Nash equilibrium, their key advance is to make the punishment itself credible (there is a reward after the punishment phase, so the punisher is willing to carry it out), so that the conclusion holds for subgame perfection, not just Nash equilibrium. They also give the dimensionality condition (full dimensionality) needed to achieve this, ensuring there is enough room to impose differentiated punishments and rewards on different players.

Result: it established the folk theorem for discounted repeated games, where the set of supportable SPE payoffs fills the entire individually rational feasible region as patience tends to the limit.

Significance and limitations: this paper is both the strongest theoretical guarantee that repeated interaction can sustain cooperation and a plain statement of the problem of equilibrium multiplicity. It is a basic reference for collusion, relational contracts, and long-term competition among platforms, but precisely because the supportable payoffs are an entire region, it cannot by itself predict where reality will land, and this "theoretical surplus" is what subsequent empirical collusion analysis must narrow with extra structure.
:::

### 5.2 Katz and Shapiro (1985)

::: {.case}
Paper: "Network Externalities, Competition, and Compatibility," American Economic Review. It is the real-world prototype for the compatibility game of this chapter, pioneering the tradition of using game theory to analyze technology standards and network effects, and profoundly influencing information systems and platform economics.

Method: in a market where a consumer's valuation of a product rises with the number of users of the same technology (a network externality), it models whether firms make their products compatible as a strategic choice, and analyzes the market equilibrium under compatibility and incompatibility. Compatibility merges the firms' user bases into one large network and amplifies network payoffs, but it also weakens the competitive advantage a large firm gets from its existing user base, so firms of different sizes and different expectations differ systematically in their preferences over compatibility.

Result: firms' private incentives for compatibility are often inconsistent with the social optimum. Firms with a large user base or optimistic about their future share tend to refuse compatibility to preserve their competitive moat, even when compatibility would enlarge overall network payoffs; small firms or late entrants would rather have compatibility to ride the large network. The market may therefore settle in a less efficient incompatible configuration, or fall into the coordination failure of this chapter.

Significance and limitations: this paper brought the core phenomena of information-technology markets, "standards battles," "interconnection," and "platform moats," into rigorous game-theoretic analysis, and is the starting point for understanding platform compatibility strategies. Its limitation is that its static framework simplifies the formation of expectations and the dynamics of preemption, and the large body of later work on installed base, dynamic standards competition, and two-sided platforms started precisely from relaxing these simplifications, a line we will return to later in this series when we cover platforms and network effects.
:::

Taken together, the two papers mark the chapter's two landing points: Fudenberg-Maskin answers "how cooperation is possible" at the most abstract level, and Katz-Shapiro brings the most concrete platform compatibility decision into a game-theoretic framework. They also foreshadow the two major themes of the structured industrial organization to come, the sustainability and identification of collusion, and the strategic analysis of platform standards and network effects.

## 6 Running case: settling the two platforms' accounts

Now we run the solution concepts of Sections 3 and 4 completely through the game between Aster and Boreal, with every number doable by hand with elementary algebra. The R code below merely transcribes and checks the by-hand computation as is.

### 6.1 The pricing game: dominance and Nash equilibrium

The pricing game is solved in one step. Check whether High is strictly dominated by Low: for Aster, facing Boreal's High, Low gets 14, higher than High's 10; facing Boreal's Low, Low gets 6, higher than High's 2. In both cases Low is strictly better, so High is strictly dominated, and the same for both firms, so the only surviving strategy combination is (Low, Low), with payoffs (6, 6), which is also the unique Nash equilibrium. The cooperative (High, High) with payoffs (10, 10), though Pareto dominant, is not an equilibrium, because each side wants to cut price secretly.

```r
# Payoff matrix: dimensions [Aster action, Boreal action], values are (Aster, Boreal)
A <- matrix(c(10, 2, 14, 6), nrow = 2, byrow = TRUE,
            dimnames = list(c("High","Low"), c("High","Low")))  # Aster payoffs
B <- t(A)  # symmetric game, Boreal payoffs are the transpose
# Is High dominated by Low (for Aster, compare column by column)
dominated <- all(A["Low", ] > A["High", ])
c(High_dominated_for_Aster = dominated)   # TRUE -> NE = (Low, Low) = (6, 6)
#> High_dominated_for_Aster
#>                     TRUE
```

### 6.2 The compatibility game: pure and mixed strategy equilibria

The compatibility game has two pure-strategy Nash equilibria, (S1, S1) and (S2, S2), verified by checking cell by cell for any profitable unilateral deviation: at (S1, S1), Aster jumping to S2 drops its payoff from 3 to 0 and Boreal jumping to S2 drops from 1 to 0, so no one wants to leave and it is an equilibrium; (S2, S2) is symmetric. The mixed strategy equilibrium is solved from the indifference conditions as $p = 3/4$ (Aster's probability of choosing S1) and $q = 1/4$ (Boreal's probability of choosing S1), the probability the two firms coordinate onto the same standard is only $3/8$, and each expects a payoff of $3/4$.

```r
# Mixed strategy: Aster chooses S1 with prob p, Boreal chooses S1 with prob q
# Boreal indifferent: p*1 = (1-p)*3 -> p = 3/4 ; Aster indifferent: 3q = 1-q -> q = 1/4
p <- 3/4; q <- 1/4
coord_prob <- p*q + (1-p)*(1-q)      # probability of successful coordination
EU_A <- 3*q                          # Aster's expected payoff (choosing S1, equal to S2 in equilibrium)
c(p = p, q = q, coord_prob = coord_prob, EU_A = EU_A)
#>       p       q coord_prob   EU_A
#>  0.7500  0.2500     0.3750  0.7500
```

### 6.3 The sequential game: backward induction and SPE

Turn the compatibility game into one where Aster moves first. Backward induction first fixes Boreal's response: if Aster chose S1 then Boreal follows with S1 for 1, better than jumping to S2 for 0; if Aster chose S2 then Boreal follows with S2 for 3, better than jumping to S1 for 0; Boreal always matches. Aster foresees this, compares 3 from choosing S1 with 1 from choosing S2, and chooses S1. The subgame perfect equilibrium is (S1, S1), payoffs (3, 1), and the first mover locks in its preferred standard. The Nash equilibrium (S2, S2) propped up by Boreal declaring "I will surely choose S2" is eliminated by SPE, because in the subgame where Aster actually chose S1, Boreal holding out with S2 for 0 is worse than following for 1, so the threat is not credible.

```r
payoff <- list(S1 = c(A = 3, B = 1), S2 = c(A = 1, B = 3), mismatch = c(A = 0, B = 0))
# Boreal's best response to each Aster choice: match (either 1 or 3 > 0)
boreal_best <- c(if_A_S1 = "S1", if_A_S2 = "S2")
# Aster expects Boreal to match, compares its own payoffs
aster_choice <- if (payoff$S1["A"] > payoff$S2["A"]) "S1" else "S2"
c(aster_choice = aster_choice, SPE_payoff_A = 3, SPE_payoff_B = 1)
#> aster_choice SPE_payoff_A SPE_payoff_B
#>         "S1"          "3"          "1"
```

### 6.4 Incomplete information: Bayesian Nash equilibrium

Boreal is the Normal type (Low dominant) with probability $\beta = 2/5$ and the Committed type (High dominant) with probability $3/5$. Aster best-responds to the type distribution, comparing expected payoffs $10 - 7\beta$ (from High) and $7 - \beta$ (from Low): plugging in $\beta = 2/5$ gives 7.2 and 6.6, so Aster prices High. The BNE is (Aster: High; Boreal: Normal type Low, Committed type High). The threshold is $\beta = 1/2$: when Committed types are more common Aster aligns high, and when Normal types are more common it retreats low.

```r
beta <- 2/5
EU_high <- beta*3 + (1 - beta)*10    # = 10 - 7*beta
EU_low  <- beta*6 + (1 - beta)*7     # = 7 - beta
c(EU_high = EU_high, EU_low = EU_low,
  aster = ifelse(EU_high >= EU_low, "High", "Low"), threshold_beta = 0.5)
#>  EU_high   EU_low    aster threshold_beta
#>    "7.2"    "6.6"   "High"          "0.5"
```

### 6.5 The repeated game: the threshold for cooperation

Under grim trigger, the discount threshold for sustaining (High, High) is $\delta \geq (T-R)/(T-P) = (14-10)/(14-6) = 1/2$. At $\delta = 1/2$, the discounted total value of cooperating forever $R/(1-\delta) = 10/0.5 = 20$ exactly equals the discounted total value of a single sneak attack $T + \delta P/(1-\delta) = 14 + 0.5\cdot 6/0.5 = 14 + 6 = 20$, and the equality of the two is precisely the indifference threshold. When $\delta$ is above $1/2$ cooperation is strictly better and sustainable; below it the sneak attack is better and cooperation collapses.

```r
R <- 10; T <- 14; P <- 6
delta_star <- (T - R) / (T - P)
coop  <- function(d) R / (1 - d)
cheat <- function(d) T + d * P / (1 - d)
c(delta_star = delta_star,
  coop_at_star = coop(delta_star), cheat_at_star = cheat(delta_star))
#>    delta_star  coop_at_star cheat_at_star
#>           0.5          20.0          20.0
```

### 6.6 Bargaining: reconciling the Nash solution and the Rubinstein solution

Let Aster and Boreal stop confronting each other and instead build an interconnected technology alliance together, creating a surplus of size $V = 100$, and negotiate how to divide it. First use the axiomatic Nash bargaining solution. When the disagreement point is $(0, 0)$ (each gets zero on breakdown) they split evenly, 50 each; if Aster holds an outside option worth 20 while Boreal has none, the disagreement point becomes $(20, 0)$, Aster first takes back 20 then splits the remaining 80 evenly with the other, getting 60, and Boreal gets 40; replacing symmetry with a weighted version where Aster's bargaining power is $\alpha = 0.6$, Aster gets $20 + 0.6\times 80 = 68$:

```r
V <- 100
nb <- function(d1, d2, alpha = 0.5) d1 + alpha*(V - d1 - d2)   # Aster's share
c(sym = nb(0, 0), asym = nb(20, 0), weighted = nb(20, 0, 0.6))
#>      sym     asym weighted
#>       50       60       68
```

Then use Rubinstein's alternating-offer solution to verify the same surplus. The two firms take turns offering, with discount factor $\delta$, and in the unique SPE the proposer takes $1/(1+\delta)$:

```r
rub <- function(delta) c(mover = 1/(1+delta), responder = delta/(1+delta))
round(100 * rub(0.80), 4)     # delta = 0.8
#>     mover responder
#>   55.5556   44.4444
round(100 * rub(0.99), 4)     # delta -> 1, approaches the symmetric Nash solution
#>     mover responder
#>   50.2513   49.7487
```

At $\delta = 0.8$ the first mover takes 55.56 and the second mover 44.44, the first-mover advantage coming from the second mover taking the discount hit for delaying; pushing $\delta$ to 0.99, the shares converge to 50.25 against 49.75, almost exactly the $50/50$ of symmetric Nash bargaining. The two solutions match up under the Nash program: the patience limit in the strategic process is precisely the symmetry in the axiomatic solution.

### 6.7 Correlated equilibrium: the value of a public coin

Return to the compatibility game, using a public fair coin as the correlating device, heads recommending both firms go S1 and tails both S2. Following the recommendation is optimal (the other also follows, and following gets a positive payoff while deviating gets zero), so this is a correlated equilibrium, and each firm's expected payoff is the equal-probability average of the two pure-equilibrium payoffs:

```r
c(CE_A = 0.5*3 + 0.5*1, CE_B = 0.5*1 + 0.5*3, mixed_NE = 3/4)
#>     CE_A     CE_B mixed_NE
#>     2.00     2.00     0.75
```

A coin lifts each firm's expected payoff from the pitiful 0.75 of the mixed equilibrium to 2, while keeping fairness between the two firms' preferred standards. The value of coordination lies not in a cleverer strategy but in a common signal everyone can see.

### 6.8 The global game: the unique tipping critical value

For that platform adoption game with adoption payoff $\theta + \ell - 1$, where $\ell$ is the fraction of adopters, complete information has both all adopting and none adopting as equilibria for $\theta \in (0, 1)$. After mixing in tiny private noise, the marginal player's posterior about the fraction of adopters is uniform on $[0, 1]$, $\mathbb{E}[\ell] = 1/2$, and the indifference condition pins the critical value down uniquely:

```r
theta_star <- 1 - 1/2          # solve theta* + E[l] - 1 = 0, E[l] = 1/2
c(theta_star = theta_star)     # theta > 1/2 -> cascade adoption, otherwise silence
#> theta_star
#>        0.5
```

The multiple equilibria collapse into a single critical value $\theta^* = 1/2$, and whether a platform tips shifts from "depends on expectations" to "depends on whether the fundamental is above $1/2$," which is exactly the certainty that global games bring to the coordination problem.

The section's accounts are now settled: the same batch of platforms, the same set of primitives, and the logic of dominance rules that one-shot pricing is doomed to a price war, pure and mixed Nash equilibria reveal the multiplicity and fragility of standards coordination, backward induction shows the first mover's locking-in advantage, the Bayesian equilibrium handles betting under private information, the repeated game gives the exact threshold for the revival of cooperation, the two bargaining solutions cut the surplus by different routes to the same end, the correlated equilibrium resolves coordination failure with a single coin, and the global game collapses the multiple equilibria into a unique tipping critical value. Every step can be done by hand and self-checked, and the ladder of solution concepts and the tools around it are laid out in full on these few small games.

## 7 Failure modes of modeling and robustness

The machinery of game theory is precise, but there are a few places where it easily goes astray when applied to real strategic situations, which this section works through one by one.

Equilibrium multiplicity is the most common trouble. The compatibility game has three Nash equilibria, and the folk theorem gives an entire region of supportable payoffs, so the theory itself often cannot select which one will be realized. This is not all bad, since it faithfully reflects the roles of coordination, expectations, and historical conventions in reality, but it means that "the game has an equilibrium" is far from "I can predict the outcome." There are several ways to cope: use refinements like SPE to shrink the equilibrium set (Section 3.4), introduce out-of-model selection principles like a focal point or a historical installed base, or simply acknowledge the multiplicity and confront it empirically as an incompleteness of the model. This last route is exactly the central difficulty later in this series when we cover entry and market structure, where multiple equilibria make one parameter correspond to several possible observations, so we can only obtain partial identification (the parameter lies in a set rather than a point), a topic to which a dedicated chapter later in this series is devoted, and here we only need to know that multiplicity is not a flaw that can be dodged with a trick but an intrinsic property of strategic interaction.

It is worth stressing that the tools of Sections 3.7 and 3.9 are made precisely to address this multiplicity trouble, but in opposite directions, so we must keep them straight when using them. Global games take the "selection" route, mixing a bit of private noise into complete information to collapse an entire region of coordination equilibria into a unique critical value, at the cost that this uniqueness depends on the specific assumptions about the noise structure, and once the noise fails to satisfy the conditions (signals highly correlated, or a public signal present), uniqueness may collapse and multiplicity may come roaring back, so it is a conditional selection rather than a master key. Correlated equilibrium takes exactly the "expansion" route: once it folds in a common signal, the supportable outcomes are larger than the convex hull of Nash equilibria, which is an advantage when explaining "why algorithms can coordinate high prices," but if used for prediction it amounts to spreading the already multiple equilibrium set even wider, so that almost anything can be rationalized, and predictive power is even weaker. In a word, what these tools change is the size of the equilibrium set, and before using them think clearly about whether you want to narrow predictions or explain more phenomena, since the two cannot both be had.

The strength of the rationality assumption is the second place. Nash equilibrium assumes players are not only rational but also correctly anticipate each other's equilibrium strategies, with common knowledge holding layer upon layer. Real people may not reason this far, and behavioral game theory experiments repeatedly show that people have limited depth of reasoning, make systematic biases, and are swayed by how things are framed. We should be cautious about taking game theory as a literal description of real behavior, since it is better suited as a benchmark for "what fully rational players who fully understand the structure would do," and the direction in which reality deviates from the benchmark is itself often informative.

The interpretation of mixed strategies is a third subtlety. Would a platform really roll dice to decide a standard? Most of the time, no. A better reading of a mixed equilibrium is a population-level frequency (a fraction of firms choose S1 and a fraction choose S2, in proportions exactly equal to the equilibrium probabilities), or the rival's uncertainty about your choice (Harsanyi's purification theorem shows that a mixed equilibrium can be seen as the limit of pure-strategy equilibria in a game with tiny private perturbations to the payoffs). When using a mixed equilibrium, be clear in your mind about which one you mean, especially when matching it to data.

The real-world counterpart of credibility and commitment is a fourth place. SPE eliminates empty threats by relying on "carrying out the threat is not optimal when it actually comes to that." But in reality players can make an originally incredible threat credible through irreversible investment, public commitment, and reputation mechanisms, for example a platform sinking cost to build capacity turns "I will start a price war" from empty talk into credible deterrence. When analyzing dynamic competition, what to watch is often not the threat itself but whether there is a commitment device that welds it in place, and this is the main line later in this series when we cover entry deterrence and strategic commitment.

Stringing these together, game theory gives us an irreplaceable language for characterizing strategic interaction, but the strength of its predictive power depends heavily on the situation: it is razor-sharp in markets with clear structure, few and focused players, and repeated encounters, while in settings with multiple equilibria, bounded rationality, and variable commitment it calls for extra care, and we should treat the model's conclusions as a conditional benchmark rather than an iron law. This is consistent with the previous chapter's stance on welfare measurement and with this series' stance on every identification strategy later on: the credibility of a method lies not in the ingenuity of the trick but in putting the assumptions on the table and being clear about how the conclusion shifts when the assumptions loosen.

## 8 Further reading

::: {.readings}
Required reading, in the suggested reading order:

- Vincent Crawford, Lecture Slides on Game Theory / Industrial Organization (corresponding to MWG Chapters 7 to 9 and 12). The organization and intuition of this chapter's ladder of solution concepts mainly follow these slides, and reading them once first will straighten out the relationships among dominance, Nash, and SPE.
- Mas-Colell, Whinston and Green, Microeconomic Theory, Chapters 7 to 9. The authoritative reference for game theory, with Chapter 8's Nash equilibrium and BNE and Chapter 9's SPE and sequential equilibrium being the complete sources for Sections 3.2 to 3.5 of this chapter, and the proof of the existence theorem is here too.
- Gibbons, Game Theory for Applied Economists. An application-oriented introduction to game theory, dense with examples and low in barriers, good to read alongside MWG to fill in intuition.

Further reading:

- Fudenberg and Maskin (1986). The rigorous statement of the folk theorem for discounted repeated games, the complete version of Section 4 of this chapter, with a focus on how punishment is made credible.
- Katz and Shapiro (1985). The origin of the game-theoretic analysis of network externalities and compatibility, required reading for the information systems and platform direction.
- Kreps and Wilson (1982). The original paper on sequential equilibrium, the classic treatment of the consistency of beliefs and strategies under imperfect information.
- Harsanyi (1967). The founding of the type method for games of incomplete information and of BNE, the source of the ideas in Section 3.3.
- Nash (1950a) and Rubinstein (1982). The two cornerstones of bargaining, the axiomatic solution and the strategic solution, read with Osborne and Rubinstein (1990) to see clearly how the Nash program connects the two, the complete source of Section 3.8 and the origin of the primitives for the platform bargaining analysis later on.
- Aumann (1974, 1987). The proposal of correlated equilibrium and its reinterpretation as an "expression of Bayesian rationality," the basis of Section 3.7, unavoidable for understanding algorithmic coordination.
- Carlsson and van Damme (1993) and Morris and Shin (1998). The founding of global games and its classic application to self-fulfilling crises, the source of the uniqueness conclusion of Section 3.9.
- Milgrom and Roberts (1990). The general theory of supermodular games and strategic complementarity, the source of the largest and smallest equilibria and the monotone comparative statics of Section 3.9.
- Fudenberg and Tirole, Game Theory. The standard graduate-level reference, especially complete on dynamic games and refinements.
:::

::: {.apa-refs}
- Aumann, R. J. (1974). Subjectivity and correlation in randomized strategies. *Journal of Mathematical Economics, 1*(1), 67-96. https://doi.org/10.1016/0304-4068(74)90037-8
- Aumann, R. J. (1987). Correlated equilibrium as an expression of Bayesian rationality. *Econometrica, 55*(1), 1-18. https://doi.org/10.2307/1911154
- Carlsson, H., & van Damme, E. (1993). Global games and equilibrium selection. *Econometrica, 61*(5), 989-1018. https://doi.org/10.2307/2951491
- Crawford, V. P. (2012). *Game theory / industrial organization: First-year M.Phil. micro lecture slides* [Unpublished lecture slides]. University of Oxford.
- Fudenberg, D., & Maskin, E. (1986). The folk theorem in repeated games with discounting or with incomplete information. *Econometrica, 54*(3), 533-554. https://doi.org/10.2307/1911307
- Fudenberg, D., & Tirole, J. (1991). *Game theory*. MIT Press.
- Gibbons, R. (1992). *Game theory for applied economists*. Princeton University Press.
- Harsanyi, J. C. (1967). Games with incomplete information played by "Bayesian" players, I-III: Part I. The basic model. *Management Science, 14*(3), 159-182. https://doi.org/10.1287/mnsc.14.3.159
- Katz, M. L., & Shapiro, C. (1985). Network externalities, competition, and compatibility. *American Economic Review, 75*(3), 424-440.
- Kreps, D. M., & Wilson, R. (1982). Sequential equilibria. *Econometrica, 50*(4), 863-894. https://doi.org/10.2307/1912767
- Mas-Colell, A., Whinston, M. D., & Green, J. R. (1995). *Microeconomic theory*. Oxford University Press.
- Milgrom, P., & Roberts, J. (1990). Rationalizability, learning, and equilibrium in games with strategic complementarities. *Econometrica, 58*(6), 1255-1277. https://doi.org/10.2307/2938316
- Morris, S., & Shin, H. S. (1998). Unique equilibrium in a model of self-fulfilling currency attacks. *American Economic Review, 88*(3), 587-597.
- Nash, J. F. (1950a). The bargaining problem. *Econometrica, 18*(2), 155-162. https://doi.org/10.2307/1907266
- Nash, J. F. (1950b). Equilibrium points in n-person games. *Proceedings of the National Academy of Sciences, 36*(1), 48-49. https://doi.org/10.1073/pnas.36.1.48
- Osborne, M. J., & Rubinstein, A. (1990). *Bargaining and markets*. Academic Press.
- Rubinstein, A. (1982). Perfect equilibrium in a bargaining model. *Econometrica, 50*(1), 97-109. https://doi.org/10.2307/1912531
- Selten, R. (1975). Reexamination of the perfectness concept for equilibrium points in extensive games. *International Journal of Game Theory, 4*(1), 25-55. https://doi.org/10.1007/BF01766400
- Topkis, D. M. (1978). Minimizing a submodular function on a lattice. *Operations Research, 26*(2), 305-321. https://doi.org/10.1287/opre.26.2.305
:::
