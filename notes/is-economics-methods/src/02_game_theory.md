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

Both platforms know that holding prices high earns each of them 10, while a price war leaves each with only 6. The most natural forecast in the boardroom is that management on both sides can do arithmetic, so of course they will pick 10. When each firm does its own calculation, it finds that, whatever the rival does, cutting its own price is the better move. Two rational decision makers thus walk together toward an outcome neither of them likes. What is missing is not the arithmetic; it is putting the rival's arithmetic into your own decision.

Once payoffs depend on how others respond, "which outcome is best" no longer coincides with "which outcome can be sustained." Game theory writes the players, strategies, information, and payoffs as a set of primitives, then uses different equilibrium concepts to look for self-fulfilling predictions. Repeated interaction adds a twist: cooperation that cannot be sustained in a one-shot game may become an equilibrium because of future punishment, but the threat of "never cooperating again" must actually be worth carrying out; and when firms cannot see their rivals' actions clearly, that punishment keeps occurring even though no one has betrayed anyone.

The chapter proceeds as follows. Section 1 uses a pricing game between two fictional platforms, Osprey and Boreal, to set out a prediction that is taken for granted, and then presents two numbers that an empirical researcher will run into. Section 2 writes down the primitives of a game and the chapter's target quantity. Section 3 climbs the ladder of rationality requirements from weak to strong, using numbered assumptions and theorems: dominance, Nash equilibrium and its existence, Bayesian games, subgame perfection, weak perfect Bayesian equilibrium, correlated equilibrium and learning, quantal response equilibrium (QRE), bargaining, supermodular games, and global games; each theorem is followed by a proof sketch and an "In words" paragraph. Section 4 covers repeated games: the threshold for cooperation, the folk theorem, and price wars on the equilibrium path. Section 5 reads three anchoring papers. Section 6 settles the accounts on five data blocks: the hand-computed block G0, the QRE block G1 of 4,000 local markets, the price-war panel G2 of 200 markets, the learning block G3, and the global-game block G4. Section 7 discusses failure modes. The later models of pricing, entry, platform competition, and algorithmic coordination all come down to the same question: when everyone is responding to everyone else, which prediction will not be overturned by the players' own best responses?

## 1 A Prediction Taken for Granted

::: {.case}
Two digital platforms, Osprey and Boreal, compete in the same market. To shrink the strategic interaction to a size small enough to do by hand, look first only at their pricing stance, with each firm choosing one of two options: hold a premium (a high price, call it High) or launch an undercut (a low price, call it Low). The monthly profits under the four combinations (in million USD, the first number for Osprey, the second for Boreal) are in the table below. When both hold high prices, each earns 10, the best outcome they can reach together. Whoever cuts price unilaterally uses the low price to grab many of the other's users, so the cutter earns 14 and the firm that stays high is left with 2. When both cut, they fall into a price war and each earns 6.

| Osprey \ Boreal | High | Low |
|---|---|---|
| High | 10, 10 | 2, 14 |
| Low | 14, 2 | 6, 6 |
:::

An analyst forecasts as follows. Both firms are smart, and both can see that holding high prices together earns 10 each, 20 combined, which is the best they can do. There is no reason to hand profit over to a price war, so the market will settle with both firms high, for a predicted monthly profit of 10 each and an industry total of 20.

This prediction sounds perfectly reasonable, yet it is almost certainly wrong. The trouble is that it assumes the two firms will hold back for the common good, but when each firm does its own calculation, holding back is not its best choice. Take Osprey's side: if Boreal is high, Osprey can raise its payoff from 10 to 14 by cutting; if Boreal is low, matching the cut takes Osprey from 2 back to 6. Whatever Boreal does, cutting is always better, and by symmetry the same holds for Boreal. So both cut, the market settles at the mutually damaging (Low, Low), each earns 6, and the industry total is 12. The forecast of 20 is two thirds above the actual 12, and the mistake is one of direction: the high-price cooperation it predicts is exactly the outcome that cannot survive in this game.

The right question is not "which outcome is good for both sides," but "given that the rival will respond in the way that is best for itself, which outcome sustains itself." Section 3 writes "self-sustaining" as the precise definition of a Nash equilibrium; Section 4 brings a twist: if the two firms meet not once but month after month, that 20 can be brought back to life under a clear condition.

The empirical researcher's trouble does not stop there. The two platforms choose an interconnection standard in 4,000 local markets, and an analyst wants to use the choice frequencies to back out Osprey's payoff from coordinating on standard S1 (true value: USD 3 thousand a month). When markets that behave alike are grouped before inverting, the average estimate over 1,000 simulations is 3.01; when markets that play different equilibria are pooled in one group, the average estimate is only 1.09, and the middle 95% of the estimates lies far below the truth. Another data set holds prices for 200 markets over 300 months: in the tier with the largest demand shocks, 43.5% of market-months are in a price war, while no firm in the simulation ever betrayed the others. The first number says how data can mislead when the equilibrium is not unique; the second says that "seeing a price war" cannot be read as "collusion breaking down." Section 6 recomputes these numbers from scratch.

## 2 The Primitives of a Game and the Target Quantity

To analyze a strategic situation, we first have to write it down completely. A noncooperative game is pinned down by four sets of primitives: the players $i = 1, \dots, n$, the strategies available to each, the information each holds when choosing, and each player's payoff under each strategy profile. Payoffs are already measured in von Neumann-Morgenstern utility (the expected utility of Theorem 15 in Section 4.3 of Chapter 1), so when players face randomness they maximize expected payoff. The table in Section 1 is exactly such a game: the players are Osprey and Boreal, the strategies are each {High, Low}, the two firms decide simultaneously without knowing the other's current choice, and the payoffs are in the table.

We need to distinguish an action (a single concrete choice) from a strategy. A strategy is a complete plan of action: it specifies what the player does at every decision point at which it might have to move, even points that the equilibrium never reaches. In a simultaneous one-shot game the two coincide; in a game with an order of moves and information updating, whether a threat is credible is exactly the question of whether the action that carries it out is optimal at the corresponding decision point, and this is the heart of Section 3.5.

The same game has two common representations. The normal form (also called the strategic form) compresses the game into a payoff table, like the one in Section 1; it erases the details of timing and information. The extensive form draws the game as a tree: nodes are decision points, branches are actions, an information set groups together the nodes a player cannot tell apart, and payoffs hang on the leaves. An extensive form can be compressed into a normal form but not the other way around, and what is lost is exactly the timing, so solution concepts based on the normal form miss the credibility problem in dynamic games.

The notation is as follows: the pure strategy set of player $i$ is $\mathcal{S}_i$, a pure strategy is $s_i \in \mathcal{S}_i$, the strategy profile of the other players is $s_{-i}$, and the payoff function is $u_i(s_i, s_{-i})$. A mixed strategy $\sigma_i \in \Delta(\mathcal{S}_i)$ is a probability distribution over $\mathcal{S}_i$, and $u_i(s_i, \sigma_{-i})$ is the expected payoff when $i$ plays $s_i$ and the others randomize according to $\sigma_{-i}$. In a 2×2 game, $p$ and $q$ are the probabilities that Osprey and Boreal, respectively, choose their first strategy.

The chapter's target quantity is a prediction of which strategy profile or profiles rational players will land on, given the primitives. The prediction is delivered by a solution concept, and it is often a set rather than a point. Empirical work therefore needs two layers: the bottom layer is primitives such as payoffs, and the top layer is the equilibrium set, together with the selection rule that picks out the one actually played. Both numbers of Section 1 are stuck on the second layer: to back payoffs out of choices, you first have to know which equilibrium each market plays; to read price wars as evidence about collusion, you first have to know which collusive equilibrium the firms play. So say first whether you are predicting the equilibrium set, one equilibrium, or an outcome under a selection rule.

## 3 The Core Structure: A Ladder of Solution Concepts

This section first uses Osprey's and Boreal's numbers to show what needs to be proved, then states numbered assumptions and theorems, followed by proof sketches in which each step names the assumption it uses. After each theorem, an "In words" paragraph answers five questions: under which assumptions the result characterizes which quantity; what it becomes when the assumptions fail; whose choice or which primitive determines the conclusion; what the data can test; and what the case's numbers say. Assumptions carry the chapter's prefix G.; theorems are numbered consecutively within the chapter, with MWG's section, proposition, and page numbers in parentheses, and results that MWG does not contain cite the original papers. Sections 3.1 to 3.6 move from weak to strong requirements of rationality: dominance, Nash equilibrium and its existence, incomplete information, and credibility in dynamic games; Sections 3.7 to 3.11 are several companion tools: correlated equilibrium and learning, QRE, bargaining, supermodular games, and global games.

### 3.1 Dominance, Iterated Deletion, and Rationalizability

The minimal requirement of rationality is not to choose a strategy that is worse whatever the others do. If there is a mixed strategy $\sigma_i'$ such that for all $s_{-i}$ we have $u_i(\sigma_i', s_{-i}) > u_i(s_i, s_{-i})$, we say $s_i$ is strictly dominated (MWG Definition 8.B.4, p. 240); the dominating strategy may be mixed, and a pure strategy is sometimes dominated only by some mixture. In the pricing game of Section 1, High is strictly dominated by Low for both firms, and a single step deletes down to (Low, Low). Most games are not this clean, but we can delete iteratively: once dominated strategies are deleted, some strategies become dominated in the smaller game that remains, and we delete again. Another route starts from beliefs: repeatedly delete the strategies that are not a best response to any conjecture about the rivals, and what remains is the set of rationalizable strategies.

::: {.assumption}
**G.1 (finite game and common knowledge; MWG Definition 7.D.2, p. 230; text of 7.C, p. 226; text of 8.B, p. 239)** (a) The number of players is finite, each player's pure strategy set $\mathcal{S}_i$ is finite, and the payoff $u_i$ is a von Neumann-Morgenstern utility over strategy profiles; (b) the structure of the game, that is, the players, strategies, and payoffs, and under incomplete information the prior distribution of types, is common knowledge: everyone knows it, everyone knows that everyone knows it, and so on ad infinitum; (c) every player is rational, that is, maximizes expected payoff under its beliefs about the others' behavior, and this too is common knowledge.
:::

G.1(a) rules out continuous prices and quantities, which need G.2 of Section 3.3. (b) and (c) are knowledge assumptions stacked in layers: deleting one round of dominated strategies needs only that everyone is rational, the second round needs "I know you are rational," and round $k$ needs mutual knowledge of order $k - 1$ (MWG p. 239). Real firms often do not know their rivals' payoffs that well, so (b) is the condition that fails most often; when it fails, Section 3.4 writes the uncertainty about payoffs as types.

::: {.theorem}
**Theorem 1 (rationalizable strategies and iterated deletion; MWG Definitions 8.C.1 and 8.C.2, pp. 242-243; text on p. 245)** Under G.1, (i) every rationalizable strategy survives the iterated deletion of strictly dominated strategies; (ii) in two-player games the two sets are equal; with three or more players, when conjectures have the rivals randomize independently, the rationalizable set can be strictly smaller, and when conjectures allow the rivals' actions to be correlated, the two coincide.
:::

The proof sketch has two steps; MWG leaves the proof of (ii) as Exercise 8.C.3. Step one (G.1(c)): a strictly dominated strategy is not a best response to any conjecture, because the mixed strategy that dominates it does better against every profile of the rivals, and still does better after taking expectations under the conjecture; each round of deleting never-best-response strategies removes at least the strictly dominated ones, and deleting round after round under the common knowledge of G.1(b)(c) gives (i). Step two (G.1(a), two players): a strategy that is not strictly dominated is always a best response to some mixed conjecture; this is a separation argument for a convex set in finite dimensions, namely the rival's set of mixed strategies; with three or more players, independent conjectures are product distributions, the set they form is not convex, the separation argument fails, and so correlated conjectures must be allowed.

In words: with rationality and its common knowledge alone, what can be predicted is a set, the rationalizable strategy profiles. Deletion often leaves more than one profile standing, so G.1 by itself usually cannot lock in a unique prediction, and the next subsection brings in a stronger condition. When knowledge goes only to order $k$, only $k + 1$ rounds can be deleted, and the set is larger. The set is determined entirely by payoffs and involves no selection rule; testing it requires knowing the payoffs; laboratory experiments often use guessing games to measure how many rounds of reasoning people actually perform (Section 7). Osprey's numbers: in the pricing game a single round deletes down to (Low, Low), each getting 6; in the compatibility game (Section 3.2) no strategy is dominated, and the rationalizable set is all four cells.

### 3.2 Nash Equilibrium: Mutual Best Responses

Nash equilibrium writes "self-sustaining" as a precise condition: a strategy profile $s^*$ is a Nash equilibrium if each player's strategy is a best response to the others' strategies (MWG Definition 8.D.1, p. 246),

$$u_i(s_i^*, s_{-i}^*) \geq u_i(s_i, s_{-i}^*) \quad \text{for all } s_i,\ \text{all } i$$

At the equilibrium, no player can make itself better off by unilaterally changing strategy. It is stronger than rationalizability; what it adds is that beliefs confirm one another: each player's expectation of the others' behavior is exactly the others' actual equilibrium behavior. The definition for a mixed strategy profile $\sigma^*$ is the same (MWG Definition 8.D.2, p. 250). In the pricing game, (Low, Low) is the unique Nash equilibrium, and a firm that unilaterally switches back to High only falls to 2.

Nash equilibrium is not always unique, as the compatibility game shows. Osprey and Boreal each choose an interconnection standard, S1 or S2. Only with the same standard can they interconnect and each reap a network payoff; if they go separate ways, payoffs are zero. Osprey prefers S1 and Boreal prefers S2, but both would rather accommodate the other than fail to connect. Monthly payoffs (million USD) are as follows.

| Osprey \ Boreal | S1 | S2 |
|---|---|---|
| S1 | 3, 1 | 0, 0 |
| S2 | 0, 0 | 1, 3 |

(S1, S1) and (S2, S2) are both pure-strategy Nash equilibria, but they divide the payoffs differently, and when the firms decide simultaneously neither knows which one to bet on: the existence of equilibrium does not guarantee that coordination is achieved. There is also a mixed-strategy equilibrium. MWG Proposition 8.D.1 (p. 250) says that $\sigma$ is a Nash equilibrium if and only if every pure strategy in each player's support earns the same payoff, no lower than any pure strategy outside the support, so each side's randomization must leave the other exactly indifferent. Making Boreal indifferent: S1 gives $p$, S2 gives $3(1-p)$, and solving gives $p = 3/4$; making Osprey indifferent: $3q = 1 - q$, so $q = 1/4$. The probability that the two firms coordinate on the same standard is $pq + (1-p)(1-q) = \tfrac{3}{4}\cdot\tfrac{1}{4} + \tfrac{1}{4}\cdot\tfrac{3}{4} = \tfrac{3}{8}$; with probability $5/8$ they fail to connect, and each side's expected payoff is only $3/4$, below the 3 or 1 the two firms get in either pure equilibrium. This is not anyone's mistake but the inevitable result of rational randomization under simultaneous choice: each side is correctly guarding against the other choosing the other standard.

::: {.theorem}
**Theorem 2 (existence of Nash equilibrium in finite games; MWG Proposition 8.D.2, p. 252; Nash 1950b)** Under G.1(a), every finite game has at least one Nash equilibrium in mixed strategies.
:::

The proof sketch has four steps, along the route of MWG Chapter 8, Appendix A (Lemma 8.AA.1, p. 260). Let $\Sigma = \Delta(\mathcal{S}_1) \times \dots \times \Delta(\mathcal{S}_n)$ and let the best-response correspondence be $\mathrm{BR}(\sigma) = \big(\mathrm{BR}_1(\sigma_{-1}), \dots, \mathrm{BR}_n(\sigma_{-n})\big)$, with $\mathrm{BR}_i(\sigma_{-i}) = \arg\max_{\sigma_i} u_i(\sigma_i, \sigma_{-i})$; $\sigma^*$ is a Nash equilibrium if and only if it is a fixed point, $\sigma^* \in \mathrm{BR}(\sigma^*)$. Step one (G.1(a)): $\Sigma$ is a product of finitely many simplices, nonempty, compact, and convex; expected payoff is a polynomial in the probabilities, hence continuous, so best responses are nonempty (the Weierstrass theorem, MWG Theorem M.F.2). Step two (expected payoff is linear in one's own mixture): the best-response set is a face of the simplex, hence convex. Step three ($u_i$ continuous): if $\sigma^k \to \sigma$, $\tau^k \in \mathrm{BR}(\sigma^k)$, and $\tau^k \to \tau$, taking limits gives $\tau \in \mathrm{BR}(\sigma)$, that is, the graph is closed; MWG uses this limit argument directly, and it is Berge's maximum theorem (MWG Theorem M.K.6, p. 963) in the special case of a fixed constraint set. Step four (G.1(a)): a nonempty, convex-valued correspondence with a closed graph on a nonempty, compact, convex set has a fixed point (Kakutani's theorem, MWG Theorem M.I.2, p. 953).

In words: under G.1(a), the equilibrium set is nonempty. Mixed strategies are therefore not a technical patch: without mixing, the strategy space is a finite set of points, not convex, and matching pennies has no pure-strategy equilibrium. The theorem does not guarantee uniqueness, still less which equilibrium will be observed; in degenerate games (where some player is exactly indifferent among several of its strategies) there can even be a whole continuum of mixed equilibria. Which one is picked is decided by things outside the model: history, focal points, the order of moves (Section 3.5), or learning (Section 3.7). For data, multiple equilibria make one set of parameters correspond to several observable outcomes; Sections 3.1 and 3.4 of Chapter 30 explain, respectively, why this breaks point identification and how to turn to partial identification. Osprey's numbers: the compatibility game has three equilibria, (S1, S1) gives (3, 1), (S2, S2) gives (1, 3), and the mixed equilibrium $(p, q) = (3/4, 1/4)$ has a coordination probability of 0.375 and an expected payoff of 0.75 for each firm.

::: {.intuition}
Plot the mixed equilibrium in a coordinate system, with the probability $p$ that Osprey chooses S1 on the horizontal axis and the probability $q$ that Boreal chooses S1 on the vertical axis, and draw a best-response curve for each. Both are staircases with a kink: when $q > 1/4$ Osprey chooses only S1, when $q < 1/4$ only S2, and when $q = 1/4$ any mixture will do; Boreal's kink is at $p = 3/4$. The two staircases cross at three points: the two corners are the pure equilibria, and the interior point $(3/4, 1/4)$ is the mixed equilibrium. The existence theorem is just the general version of this geometric fact: best responses cannot dodge each other in a compact convex strategy space, and the crossing point is the equilibrium.
:::

### 3.3 Continuous Games: Existence of Pure-Strategy Equilibrium

Prices and quantities are continuous variables. In each local market of block G2 in Section 6, the two firms choose monthly output (thousand units), the price is USD 10 minus their total output, and unit cost is USD 2. Given the rival's output $q'$ (in this subsection $q$ is an output), a firm's best output is $(8 - q')/2$; the two best-response lines cross where each produces $8/3$, at a price of USD 4.67, each earning USD 7.11 thousand a month. There is no randomization here, and yet the crossing exists; the conditions below say why.

::: {.assumption}
**G.2 (regularity conditions for continuous games; conditions (i)(ii) of MWG Proposition 8.D.3, p. 253)** For each $i$: (a) $\mathcal{S}_i$ is a nonempty, convex, compact subset of some Euclidean space; (b) $u_i$ is continuous in the strategy profile and quasiconcave in the player's own strategy $s_i$.
:::

G.2(a) allows bounded continuous choices such as output and price. The quasiconcavity in (b) requires that a player's payoff not have several separate peaks as its own strategy varies. G.2 fails in several common settings: in homogeneous-goods Bertrand competition, payoffs are discontinuous where the two prices are equal (the slightly cheaper firm takes all demand); integer choices such as capacity or the number of stores are not convex; and pricing with network effects can give payoffs two peaks in the firm's own price.

::: {.theorem}
**Theorem 3 (existence of pure-strategy Nash equilibrium in continuous games; MWG Proposition 8.D.3, p. 253; Debreu 1952; Glicksberg 1952; Fan 1952)** Under G.2, the game has a pure-strategy Nash equilibrium.
:::

The proof sketch has three steps, along the same route as Theorem 2, with convexity coming from a different source (MWG Lemma 8.AA.1, p. 260). Step one (compactness in G.2(a), continuity in G.2(b)): $u_i(\cdot, s_{-i})$ is a continuous function on a compact set, so best responses are nonempty. Step two (quasiconcavity in G.2(b), convexity in G.2(a)): the set of maximizers of a quasiconcave function on a convex set is convex, so best responses are convex-valued; Theorem 2 relied on the linearity of expected payoff in the mixture, and here quasiconcavity does the job. Step three (continuity in G.2(b)): the same limit argument gives a closed graph, and Kakutani's theorem gives a fixed point on $\mathcal{S}_1 \times \dots \times \mathcal{S}_n$.

In words: under G.2, there is an equilibrium without randomization, and this is exactly what guarantees the existence of the Cournot equilibrium of Theorem 7 in Section 4.1 of Chapter 4. When G.2 fails, it does not follow that there is no equilibrium, only that the theorem no longer vouches for one: in homogeneous-goods Bertrand, price equal to marginal cost is still an equilibrium; and when quasiconcavity fails but payoffs are continuous, MWG (p. 253) notes that a mixed-strategy equilibrium still exists. The conditions are determined by technology and demand, and whether quasiconcavity holds depends on the shape of the demand curve. What the data can test is the equilibrium's first-order condition; Section 3.1 of Chapter 29 backs marginal cost out of the first-order condition for pricing. Osprey's numbers: output in G2 is restricted to $[0, 8]$ (producing more must lose money), and the payoff is strictly concave in the firm's own output, so G.2 holds; the two best-response lines both have slope $-1/2$ and cross only once, so producing $8/3$ each, at a price of USD 4.67 and a monthly profit of USD 7.11 thousand, is the unique pure equilibrium, and it is also the outcome to which the two firms revert during the price wars of Section 4.

### 3.4 Incomplete Information and Bayesian Nash Equilibrium

So far both firms fully know each other's payoffs. In reality, a platform's cost, technological readiness, and strategic resolve are often known only to itself. The approach of Harsanyi (1967) is to give each player a type $\theta_i$ that represents its private information, and to assume that the prior distribution of types is common knowledge (G.1(b)); in this subsection $\theta_i$ denotes only the type, and in Sections 3.10 and 3.11 the letter has other meanings. The equilibrium concept is generalized accordingly to Bayesian Nash equilibrium (BNE): each type of each player best-responds to the distribution of the rivals' types and to their type-dependent strategies.

Here is a platform version of the example: Osprey is uncertain about Boreal's cost structure. With probability $\beta$ Boreal is the Normal type, for which Low is a dominant strategy; with probability $1 - \beta$ it is the Committed type (it has just signed an exclusive deal for high-value content and monetizes through premium positioning), for which High is a dominant strategy. Boreal's equilibrium strategy is obvious; the difficulty is on Osprey's side (its payoff table differs from that of Section 1): if Boreal prices High, Osprey gets 10 from High and 7 from Low; if Boreal prices Low, it gets 3 from High and 6 from Low. The expected payoffs are

$$\mathbb{E}[u_{\text{High}}] = \beta\cdot 3 + (1-\beta)\cdot 10 = 10 - 7\beta, \qquad \mathbb{E}[u_{\text{Low}}] = \beta\cdot 6 + (1-\beta)\cdot 7 = 7 - \beta$$

Osprey chooses High if and only if $10 - 7\beta \geq 7 - \beta$, that is $\beta \leq 1/2$: the more likely Boreal is the high-pricing Committed type, the more Osprey should align with the high price.

::: {.theorem}
**Theorem 4 (characterization of BNE and existence with finitely many types; MWG Definition 8.E.1 and Proposition 8.E.1, pp. 255-256; existence follows from Theorem 2 and is not stated separately in MWG 8.E, see footnote 7, p. 255)** Let types be drawn from a common prior and seen only by their owner, and let a strategy be a decision rule $s_i(\theta_i)$. Under G.1, (i) when types are finite (for infinitely many types, see MWG footnote 8, p. 256), a profile of decision rules is a BNE, that is, a Nash equilibrium of the game among types, if and only if every type $\theta_i$ that occurs with positive probability chooses optimally given the conditional distribution of the others' types and the others' decision rules; (ii) when types and actions are both finite, a BNE in mixed decision rules exists.
:::

The proof sketch has two steps. Step one (the common prior of G.1(b)): ex ante expected payoff equals the sum, over types, of the product of each type's probability and its conditional expected payoff; each term depends only on that type's choice, so the sum is maximized if and only if each type with positive probability maximizes on its own, which is (i). Step two (G.1(a)): with finitely many types and actions there are only finitely many decision rules, the game among types is a finite game, and Theorem 2 gives a mixed equilibrium, which is (ii). A pure-strategy BNE need not exist; the purification of Harsanyi (1973) reads a mixed equilibrium as the limit of pure-strategy equilibria when payoffs carry tiny private perturbations (MWG p. 257).

In words: under G.1, BNE predicts the action of each type; a player best-responds not to a specific rival but to a distribution of rival types, and this idea is the direct forerunner of the auctions and mechanism design of Chapter 3. When the common prior fails (the two sides disagree about the type distribution), the very definition of equilibrium must be rewritten. The conclusion is determined by the prior: change the prior and the best response can flip entirely. For data, the information structure is a modeling choice: if each firm knows only its own cost shock, an entry game may go from multiple equilibria to a unique one, which is exactly what Section 3.3 of Chapter 30 uses. Osprey's numbers: at $\beta = 2/5$, Osprey gets 7.2 from High and 6.6 from Low, so in the BNE Osprey prices High, Boreal's Normal type prices Low, and its Committed type prices High; the threshold is $\beta = 1/2$.

### 3.5 Dynamic Games: Backward Induction, Subgame Perfection, and One-Shot Deviation

The simultaneous-choice framework misses one dimension: order. Turn the compatibility game into a sequential one: Osprey chooses a standard first and announces it, and Boreal chooses after seeing it.

![Extensive form of the sequential compatibility game: Osprey chooses S1 or S2 first, and Boreal chooses after observing it; backward induction gives the subgame perfect equilibrium: Osprey chooses S1 and Boreal follows, with payoffs (3, 1) (million USD a month).](assets/fig/fig_02_gametree.svg)

Backward induction works back from the ends of the tree: if Osprey chose S1, Boreal gets 1 from following with S1 and 0 from switching to S2; if Osprey chose S2, Boreal gets 3 from following with S2 and 0 from switching to S1, so Boreal always matches. Osprey foresees this: S1 would earn it 3 and S2 would earn it 1, so it chooses S1. The coordination failure under simultaneous choice is dissolved in one stroke by the order of moves, and this is exactly where the first-mover advantage comes from. The game has other Nash equilibria too: Boreal declares "whatever you choose, I will choose S2," Osprey is scared into choosing S2 as well, which gives the Nash equilibrium (S2, S2); but once Osprey has actually chosen S1, Boreal holding out with S2 gets 0, worse than the 1 from following, so this is an empty threat.

::: {.assumption}
**G.3 (observed actions and discounting; MWG Definition 7.C.1, p. 226; Fudenberg and Tirole 1991, §4.2)** (a) The actions of each stage are seen by everyone before the next stage begins (a game of perfect information is the special case in which only one player moves at each stage), so every history $h^t$, the sequence of actions in the first $t$ periods, starts a subgame; (b) the number of periods is finite, or it is infinite with payoffs discounted by a discount factor $\delta \in (0, 1)$ and bounded stage payoffs.
:::

G.3(a) rules out situations in which rivals' actions cannot be seen: oligopolists who see only the market price and not the rival's output do not satisfy it, and Section 4.2 switches to strategies that depend only on the public price. G.3(b) rules out infinite-horizon games without discounting, where Theorem 6 can fail.

::: {.theorem}
**Theorem 5 (Zermelo's theorem and subgame perfect equilibrium; MWG Proposition 9.B.1, p. 272; Proposition 9.B.2, p. 276)** Let the game be finite (G.1(a)) and of perfect information (a special case of G.3(a)). (i) There is a pure-strategy Nash equilibrium that can be obtained by backward induction; (ii) it is a subgame perfect equilibrium (SPE), that is, it is a Nash equilibrium in every subgame (MWG Definition 9.B.2, p. 275); if no player has the same payoff at two terminal nodes, the SPE is unique.
:::

The proof sketch has three steps. Step one (finiteness): a decision maker at the last level has finitely many choices, so an optimal choice exists; replace the whole subtree with that choice and roll back level by level to the root, so every step is well defined. Step two (perfect information): what the resulting strategy profile prescribes at every node is an optimal action given optimal continuation play; restricted to any subtree, it is that subtree's backward-induction outcome and hence a Nash equilibrium of that subgame; the case of the root node is (i). Step three (perfect information, no ties): every node starts a subgame, so any SPE is a backward-induction outcome; backward induction is unique at every step, so the SPE is unique.

In words: in a finite game of perfect information, SPE predicts not only the actions on the path but also the actions in every situation that "should not happen," and requires them to be worth carrying out if play does get there, so empty threats are eliminated. With ties in payoffs the SPE can be non-unique; under imperfect information there are too few subgames for SPE to get traction, and the refinements of Section 3.6 are needed. The conclusion is determined by the order of moves: with the same payoff table, the equilibrium tilts toward whoever moves first. Testing such predictions requires exogenous variation in the order of moves, for example launch dates set by regulation or by technical factors. Osprey's numbers: backward induction gives (S1, S1) with payoffs (3, 1); the (S2, S2) that rests on an empty threat is eliminated; Osprey gets 3, far above the 0.75 of the mixed equilibrium under simultaneous choice.

An infinite-horizon game has no "last level," and backward induction has nowhere to start. The principle below shrinks the check to a one-shot deviation.

::: {.theorem}
**Theorem 6 (one-shot deviation principle; Fudenberg and Tirole 1991, §4.2)** Under G.3, a strategy profile $s$ is an SPE if and only if no player has a profitable one-shot deviation: there is no $i$ and history $h^t$ such that $i$, changing its action only at $h^t$ and following $s_i$ thereafter, obtains a higher discounted payoff conditional on reaching $h^t$.
:::

The proof sketch has two steps. Necessity is obvious: an SPE allows no profitable deviation. Sufficiency, step one (the finite horizon of G.3(b)): if some deviating strategy is profitable, look at the history of its last deviation; only a one-shot deviation remains there, which by hypothesis is not profitable, so switching it back to $s_i$ loses nothing; removing the deviations backward one at a time, the deviating strategy is no better than $s_i$, a contradiction. Step two (discounting and boundedness in G.3(b)): take a deviation that gains $\varepsilon > 0$; after a late enough period $t_0$, switch it back to $s_i$; the loss is at most $\delta^{t_0}$ times the range of payoffs divided by $1 - \delta$, a bound that can be made smaller than $\varepsilon/2$, so we obtain a strategy that deviates only finitely often and still gains, and we are back in step one. The property that distant payoffs can be ignored is what Fudenberg and Tirole call continuity at infinity.

In words: under G.3, checking whether an infinite-horizon strategy is an equilibrium only requires checking, in every state, whether "deviating for just one step" pays. Without discounting the principle can fail, because the gains from deviating can be piled up at infinity. It turns the incentive constraints of a repeated game into finitely many inequalities, and both the grim trigger threshold of Section 4.1 and the price-war equilibrium of Section 4.2 rely on it. Osprey's numbers: for the grim trigger of Section 4.1 at $\delta = 1/2$, the discounted values of keeping the deal and of deviating once are both 20; the four trigger strategies chosen in block G2 for $\delta = 0.95$ all pass the one-shot deviation check, while the same strategies all fail at $\delta = 1/2$, where the most profitable one-shot deviations gain 0.125, 0.075, 0.141, and 0.125 thousand USD (after re-optimization, collusion is still possible when noise is small, Section 4.2).

### 3.6 Imperfect Information: Weak Perfect Bayesian Equilibrium and Sequential Equilibrium

Draw the simultaneous compatibility game as a tree: Osprey chooses first, Boreal cannot see its choice, and Boreal's two decision nodes fall into the same information set $\mathcal{I}$. Apart from the whole game, this tree has no other subgame, so SPE coincides with Nash equilibrium and all three equilibria pass. To constrain actions at information sets, the equilibrium has to include the players' beliefs at those sets: a system of beliefs gives each node $h \in \mathcal{I}$ a probability $\mu(h)$, with $\sum_{h \in \mathcal{I}} \mu(h) = 1$; sequential rationality requires the action at each information set to be optimal given these beliefs. In this subsection $\mu(h)$ is a belief over nodes; the $\mu$ of Section 3.7 is a distribution over action profiles.

::: {.theorem}
**Theorem 7 (weak perfect Bayesian equilibrium and sequential equilibrium; MWG Proposition 9.C.1, p. 285; Definition 9.C.3, p. 285; Definition 9.C.4 and Proposition 9.C.2, pp. 290-291; Kreps and Wilson 1982)** Let the extensive-form game be finite (G.1(a)). (i) $\sigma$ is a Nash equilibrium if and only if there are beliefs $\mu$ such that $\sigma$ is sequentially rational given $\mu$ at every information set reached with positive probability, and $\mu$ at these information sets is derived from $\sigma$ by Bayes' rule. (ii) If $\sigma$ is sequentially rational given $\mu$ at every information set and $\mu$ is derived from it by Bayes' rule wherever possible, $(\sigma, \mu)$ is called a weak perfect Bayesian equilibrium (weak PBE); it is a Nash equilibrium but not necessarily an SPE. (iii) If in addition $\mu$ is required to be the limit of the Bayes beliefs of a sequence of completely mixed strategies converging to $\sigma$, $(\sigma, \mu)$ is a sequential equilibrium; every sequential equilibrium is a weak PBE, and its $\sigma$ is an SPE.
:::

The proof sketch has two steps. Step one (G.1(a)): information sets reached with probability zero do not affect expected payoff, so how one acts there is irrelevant to Nash equilibrium; at information sets with positive probability, conditional optimality computed with Bayes beliefs is equivalent to ex ante optimality, which is (i); (ii) adds sequential rationality at off-path information sets on top of (i), so it is still a Nash equilibrium, while off-path beliefs are unconstrained, and MWG Example 9.C.5 (p. 289) gives a weak PBE that is not an SPE. Step two (G.1(a)): completely mixed strategies reach every information set, so consistent beliefs are, in every subgame, derived from the strategies themselves by Bayes' rule, and sequential rationality then gives a Nash equilibrium in every subgame, which is (iii) (MWG p. 291).

In words: these refinements share SPE's motive, requiring off-path actions to be genuinely worth carrying out under some belief; weak PBE leaves off-path beliefs completely free and is too weak, while sequential equilibrium ties beliefs to "a tiny trembling hand." The signaling games of Theorem 2 in Section 3.2 of Chapter 3 add conditions on top of weak PBE, and with several receivers they also require the receivers to hold common beliefs about the sender. Off-path beliefs are the modeler's choice, and data never reveal them, so the sensitivity of the conclusion to them should be reported separately. Osprey's numbers: in the simultaneous compatibility game, the probability that Boreal, at its information set, assigns to the node "Osprey chose S1" is 1 in (S1, S1), 0 in (S2, S2), and $p = 3/4$ in the mixed equilibrium, all three uniquely determined by Bayes' rule.

### 3.7 Correlated Equilibrium and Learning

Nash equilibrium assumes by default that everyone randomizes independently; in reality players can often see a common signal and act on it. Imagine a mediator who draws a strategy profile according to a joint distribution $\mu$ and privately tells each player only its own component as a recommendation; $\mu$ is a correlated equilibrium (Aumann 1974) if, for every player, after receiving the recommendation $s_i$ and forming its conditional belief about the others' actions, following it is optimal:

$$\sum_{s_{-i}} \mu(s_{-i}\mid s_i)\,\big[\,u_i(s_i, s_{-i}) - u_i(s_i', s_{-i})\,\big] \geq 0 \quad \text{for all } s_i',\ \text{all } i$$

In the compatibility game, let a fair coin that both firms can see serve as the coordinating device: heads recommends that both choose S1, tails that both choose S2. On receiving "S1," Osprey knows that Boreal will also choose S1, and following the recommendation gets it 3, better than the 0 from deviating; the same holds for Boreal. Each firm's expected payoff is 2, far above the 0.75 of the mixed equilibrium. Now take an entry game, Chicken: each firm decides whether to enter a new city; if both yield, each gets 5; a firm that enters alone gets 6 and the one that yields gets 2; if both enter, each gets 0. Its mixed equilibrium gives each firm 4, and in each of the three Nash equilibria the payoffs sum to at most 8; recommending "both yield," "only Osprey enters," and "only Boreal enters," each with probability 1/3, gives each firm an expected payoff of 4.3333 and a total of 8.6667.

::: {.theorem}
**Theorem 8 (linear characterization of correlated equilibrium and the convex hull of Nash equilibria; Aumann 1974; definition in the small-print parts of MWG 8.D and 8.E, p. 252, p. 257; direct derivation in this chapter)** Under G.1(a), (i) $\mu$ is a correlated equilibrium if and only if for all $i$, $s_i$, $s_i'$ we have $\sum_{s_{-i}} \mu(s_i, s_{-i})\,[u_i(s_i, s_{-i}) - u_i(s_i', s_{-i})] \geq 0$; these are finitely many linear inequalities in $\mu$, so the set of correlated equilibria is a convex polytope; (ii) the product distribution induced by every Nash equilibrium is a correlated equilibrium, so the set of correlated equilibria contains the convex hull of the Nash equilibrium distributions, and it can be strictly larger.
:::

The proof sketch has two steps. Step one (algebra): multiplying the conditional form of the constraint by $\mu(s_i) > 0$ gives the unconditional form in (i), and when $\mu(s_i) = 0$ both sides are zero. Step two (definitions): under a product distribution the others' actions are independent of one's own recommendation, so the conditional distribution is just the rivals' equilibrium mixture, and the optimality of Nash equilibrium is exactly the inequality in (i); the solution set of linear inequalities is convex, so convex combinations of Nash equilibrium distributions (first toss a public coin to choose which equilibrium) still satisfy it; the Chicken example gives a point outside the convex hull, which is "strictly larger."

In words: in a finite game, the set of correlated equilibria is characterized by a set of linear inequalities, so computing it is a linear program, much easier than finding Nash equilibria. It builds common signals into coordination, so the supportable outcomes can only grow: this is an advantage when explaining "why a shared algorithm can coordinate high prices," but used for prediction it stretches an already large equilibrium set even wider. Which correlated equilibrium occurs depends on who designs the signal; a pricing or recommendation algorithm shared by many parties is itself a correlating device. Osprey's numbers: under the public coin, both firms' gain from disobeying is 0, and each gets 2; in Chicken, the recommendations with probability 1/3 give each 4.3333, 8.6667 in total, more than 8, the most that the convex hull of Nash equilibria can offer.

Who plays the mediator? Hart and Mas-Colell (2000) prove that no mediator is needed and that no one even needs to know what the equilibrium is: as long as players adjust by a simple regret rule, they reach the set of correlated equilibria. Regret matching works as follows: a player who chose $j$ last period adds up, over all past periods in which it chose $j$, the extra payoff that switching to $k$ would have earned, divides by the total number of periods $t$, takes the positive part, divides again by a large enough constant (block G3 uses twice the largest absolute payoff plus 1), and switches to $k$ with that probability, staying at $j$ otherwise.

::: {.theorem}
**Theorem 9 (regret matching converges to the set of correlated equilibria; Hart and Mas-Colell 2000)** In a finite game satisfying G.1(a), if every player acts according to regret matching, the empirical distribution of joint actions converges almost surely to the set of correlated equilibria.
:::

For the proof, see Hart and Mas-Colell (2000). This book only states the result, because the proof uses Blackwell's approachability theorem, treating the regret vector as a stochastic process pushed toward the nonpositive orthant, which is beyond this chapter.

In words: in a finite game, regret matching predicts that the empirical distribution falls into the set of correlated equilibria, not that it converges to one point, and it says nothing about how the last period is played; the rule itself does not favor equilibria with higher payoffs. The conclusion is determined by the learning rule: switch to rules such as fictitious play and both the object of convergence and the conditions change. Section 2.1 of Chapter 34 lets Q-learning pricing algorithms interact in a repeated game and asks the same kind of question: where will learning algorithms stop? Osprey's numbers (block G3, 50,000 periods): over 100 independent runs of each of the compatibility game and Chicken, the largest gains from disobeying are 0.00048 and 0.00078 million USD, respectively; in every run more than 99% of the empirical distribution sits on one pure Nash equilibrium, 51% of the Chicken runs stop at (Yield, Enter) and 49% at (Enter, Yield), and not one reaches the correlated equilibrium with the 1/3 recommendations. A pure Nash equilibrium is itself a correlated equilibrium, so the theorem holds, but it does not guarantee that the better equilibrium is found.

![Block G3: under regret matching, the empirical distribution's largest gain from disobeying falls with the number of periods (log scale); mean and maximum over 100 independent runs each for the compatibility game and Chicken.](assets/fig/fig_02_regret_matching.svg)

### 3.8 QRE: Noisy Best Responses

In each local market of block G1, the two platforms play the compatibility game of Section 3.2 (with payoffs now in thousand USD a month), each making 50 interconnection decisions. QRE choice probabilities lie strictly between 0 and 1, but a single market's frequency is often exactly 0 or 1, which has no log odds. If each firm best-responds with logit noise of precision $\chi$, an equilibrium is a fixed point in probabilities: in a market without a payoff shift, at $\chi = 0.8$ the unique equilibrium has Osprey choose S1 with probability 0.6101 and Boreal with probability 0.3899; at $\chi = 2$ there are three equilibria, $(0.1249, 0.0067)$, $(0.6645, 0.3355)$, and $(0.9933, 0.8751)$. Logit QRE (McKelvey and Palfrey 1995) requires

$$\sigma_i(s_i) = \frac{\exp\{\chi\, u_i(s_i, \sigma_{-i})\}}{\sum_{s_i' \in \mathcal{S}_i} \exp\{\chi\, u_i(s_i', \sigma_{-i})\}}$$

Here $\chi \geq 0$ is the precision and $1/\chi$ the scale of the noise. In a 2×2 game this reads $p = \Lambda\big(\chi\,[u_1(1, q) - u_1(2, q)]\big)$, where $\Lambda(\cdot)$ is the logistic distribution function and $u_1(1, q)$ is Osprey's expected payoff when it chooses its first strategy and Boreal chooses its first strategy with probability $q$.

::: {.theorem}
**Theorem 10 (logit QRE: existence, two limits, and log odds; McKelvey and Palfrey 1995)** Under G.1(a), for every $\chi \geq 0$: (i) a logit QRE exists; (ii) in any logit QRE, $\ln\{\sigma_i(s)/\sigma_i(s')\} = \chi\,[u_i(s, \sigma_{-i}) - u_i(s', \sigma_{-i})]$; (iii) at $\chi = 0$ all strategies are equally likely; as $\chi \to \infty$, every limit point of a sequence of logit QRE is a Nash equilibrium.
:::

The proof sketch has three steps. Step one (G.1(a)): the logit response maps mixed strategy profiles continuously into mixed strategy profiles; $\Sigma$ is nonempty, compact, and convex; and Brouwer's fixed point theorem gives (i). Step two (algebra): dividing two logit probabilities cancels the denominator, which is (ii). Step three ((ii), G.1(a)): at $\chi = 0$ all exponents are zero; if at the limit some $s$ has an expected payoff strictly below that of another $s'$, then by (ii) $\sigma_i(s)/\sigma_i(s')$ is the exponential of $\chi$ times a negative number and tends to zero, so the limit puts probability only on best responses, which is (iii).

In words: QRE predicts choice probabilities, which are functions of payoffs and precision; every strategy occurs with positive probability, and the observed "mistakes" are part of the model. In a coordination game, once precision is high, QRE is no longer unique: in G1's game, the principal branch starting from the uniform point at $\chi = 0$ bifurcates at $\chi = 1.0699$, and the three branches tend to the three Nash equilibria (figure below); the bifurcation comes from the game's symmetry under swapping both the players and the strategies, and in a generic game the principal branch does not bifurcate. Because (ii) holds in every QRE, log odds are proportional to expected payoff differences, and this is exactly the bridge from choice frequencies back to payoffs: given the precision, the payoffs can be inverted by OLS after averaging over groups of markets; when the precision is unknown, it too can be identified from the coefficient on the payoff shifter. The bridge holds only when the equilibrium is unique or when it is known which equilibrium each market plays; Section 6.3 computes the cost of pooling them. Section 3.2 of Chapter 12 and Section 3.1 of Chapter 27 estimate logit choice models, and Chapter 30 handles estimation under multiple equilibria. Osprey's numbers: in the $\chi = 0.8$ block each shift value has only one QRE, and Osprey's probability of choosing S1 ranges from 0.171 to 0.943; in the $\chi = 2$ block every shift value has three.

![Logit QRE of block G1: without a payoff shift, the probabilities that Osprey and Boreal choose S1 as functions of the precision $\chi$; the principal branch bifurcates at $\chi = 1.070$, the dashed lines are the limiting Nash equilibria, and the dotted lines mark the values $\chi = 0.8$ and $\chi = 2$ used in G1.](assets/fig/fig_02_qre_branches.svg)

### 3.9 Bargaining: The Nash Solution and Rubinstein's Alternating Offers

Much of the interaction in the platform economy is not competition but dividing money: a platform negotiates a revenue split with a content supplier, two firms negotiate equity stakes to form a joint venture. Suppose Osprey and Boreal jointly form an interconnection technology alliance that creates a surplus of size $V = 100$. The first solution is the axiomatic one of Nash (1950a): given the set $U$ of feasible utility pairs and the disagreement point $d = (d_1, d_2)$ that each side gets if talks break down, the solution is the point that maximizes the Nash product:

$$\max_{(u_1, u_2)\in U}\; (u_1 - d_1)(u_2 - d_2)$$

::: {.theorem}
**Theorem 11 (Nash bargaining solution; MWG Proposition 22.E.1, p. 843; Nash 1950a)** Let $U$ be convex and closed and allow free disposal, let $U$ contain a point that makes both sides strictly better off than $d$, and let $\{u \in U : u \geq d\}$ be bounded. The solution that maximizes the Nash product is the only solution that satisfies all four of the following: invariance to changes in the origin and units of utility, Pareto efficiency, symmetry, and independence of irrelevant alternatives.
:::

The proof sketch has two steps (MWG pp. 843-844). Step one ($U$ convex and closed): the maximizer of the Nash product satisfies the four axioms: an affine transformation only multiplies the Nash product by a positive constant; the Nash product is increasing in each $u_i$; it is symmetric on a symmetric $U$; and on a subset that contains the original maximizer, the maximizer does not change. Step two ($U$ convex): conversely, use invariance to origin and units to move $d$ to the origin and scale the Nash solution to $(1, 1)$; $U$ lies entirely below the tangent line $u_1 + u_2 = 2$; the solution on the symmetric half-plane $\{u : u_1 + u_2 \leq 2\}$ is pinned at $(1, 1)$ by symmetry and Pareto efficiency, $U$ is a subset of it that contains $(1, 1)$, and independence of irrelevant alternatives implies that the solution on $U$ is also $(1, 1)$.

In words: under these axioms, the division is determined by the feasible set and the disagreement point alone; when both sides are risk-neutral and divide a sum of money, the reading is that each first takes back its $d_i$, and then the net surplus $V - d_1 - d_2$ is split in half. Dropping the symmetry axiom gives generalized Nash bargaining, in which one side gets $d_1 + \alpha(V - d_1 - d_2)$, the weight $\alpha$ being its bargaining power. The disagreement point and bargaining power are both primitives supplied by the modeler; what the data can show is how changes in outside options move the split. Osprey's numbers: with $d = (0, 0)$ the split is 50 to 50; when Osprey has an outside option worth 20, $d = (20, 0)$ and the split is 60 to 40; with bargaining power $\alpha = 0.6$, Osprey gets $20 + 0.6\times 80 = 68$.

The second solution is the strategic one of Rubinstein (1982): two players take turns proposing a division of a surplus of size 1; the other accepts or rejects, and a rejection moves play to the next period, when the other side proposes; the surplus shrinks by the discount factor with each period of delay. Guess the solution by stationarity: if the proposer takes $x^*$ in equilibrium, it must give the responder its continuation value from rejecting, $\delta x^*$, and keep $1 - \delta x^*$; stationarity requires that what it keeps be exactly $x^*$:

$$x^* = 1 - \delta x^* \;\Longrightarrow\; x^* = \frac{1}{1+\delta}$$

::: {.theorem}
**Theorem 12 (unique SPE of alternating offers and the limit as the interval between offers goes to zero; MWG Example 9.AA.2, pp. 298-299; Rubinstein 1982, including the case of different discount factors; for (ii), see Binmore, Rubinstein and Wolinsky 1986)** Under G.3, let the two discount factors be $\delta_1, \delta_2 < 1$, and let each side get 0 if they never agree. (i) The SPE is unique, agreement is immediate in the first period, and the first proposer gets $(1 - \delta_2)/(1 - \delta_1 \delta_2)$; when $\delta_1 = \delta_2 = \delta$ this is $1/(1+\delta)$. (ii) If offers are $\Delta$ years apart, $\delta_i = \exp(-r_i \Delta)$, with $r_i$ the yearly discount rate (not the factor price of Chapter 1), then as $\Delta \to 0$ player 1's share tends to $r_2/(r_1 + r_2)$, whoever proposes first: this is the generalized Nash solution with disagreement point $(0, 0)$ and weights inversely proportional to the discount rates.
:::

The proof sketch has three steps; for (i) the sketch follows MWG's sup-and-inf argument (Shaked and Sutton 1984) with $\delta_1 = \delta_2 = \delta$, and MWG leaves the unequal case as Exercise 9.B.12 (p. 303). Let $M$ and $m$ be the supremum and infimum of the first proposer's payoff over all SPE; the subgame from the second period on is isomorphic to the original game with the roles swapped, so the responder faces the same $M$ and $m$. Step one (G.3, subgame perfection): the responder must accept any offer of more than $\delta M$, so $m \geq 1 - \delta M$. Step two (G.3): the responder must reject any offer of less than $\delta m$; after a rejection, agreement comes at the earliest in the next period, the responder as proposer gets at least $m$, and the first proposer gets at most $\delta(1 - m) \leq 1 - \delta m$, so $M \leq 1 - \delta m$. Combining the two, $M \leq 1 - \delta(1 - \delta M)$, that is $M \leq 1/(1+\delta) \leq m$, and the supremum and infimum coincide. Step three ($\delta_i = \exp(-r_i \Delta)$): take the limit of the share in (i); expanding the numerator and denominator of $\{1 - \exp(-r_2 \Delta)\}/\{1 - \exp(-(r_1 + r_2)\Delta)\}$ to first order in $\Delta$, the ratio tends to $r_2/(r_1 + r_2)$, which is (ii).

In words: under the rules of alternating offers, patience is bargaining power, and the first-mover advantage comes from the responder paying the discount for delay; the more frequent the offers, the smaller this advantage, and the weight $\alpha$ of the generalized Nash solution thereby gets a microfoundation. The conclusion is sensitive to the bargaining protocol: MWG (p. 299) warns that with another order of offers or another specification of preferences, the shares can differ. The disagreement point in (ii) is $(0, 0)$, which is a different problem from the 60 to 40 split at $d = (20, 0)$ above; outside options need to be modeled separately. Osprey's numbers: at $\delta = 0.8$ the first proposer gets 55.56 and the other side 44.44; at $\delta = 99/100$ it is 50.25 to 49.75. Osprey's yearly discount rate is 0.10 and Boreal's is 0.05; with yearly offers Osprey gets 0.3501 when it proposes first and 0.3168 when it proposes second, and with daily offers the two are about 0.3334 and 0.3333, respectively, with limit 1/3, that is 33.33 of $V = 100$.

### 3.10 Topkis, Supermodular Games, and Organizational Complementary Systems

In each local market $m$ of block G1, a payoff shift $x_m$ goes to the S1 side (thousand USD a month; it can be read as the local installed base of S1), and Osprey gains $x_m + 4q - 1$ more from S1 than from S2: the more likely Boreal is on S1, the more S1 is worth, and the slope 4 is the sum of the two coordination payoffs. This pattern, "the more the rival does, the higher my marginal return from doing it," is called strategic complements. The same structure also appears among several organizational choices within one firm. Suppose a firm chooses AI adoption $a$, employee skill investment $s$, and decentralization of decision rights $d$ simultaneously, with profit $\pi(a,s,d;\theta)$, where $\theta$ is a parameter that raises marginal returns. In this subsection, $s$ and $d$ denote only these two choices and are unrelated to the strategy $s_i$ and the disagreement point $d$ used earlier in the chapter. If any two choices have increasing differences, for example

$$\frac{\partial^2\pi}{\partial a\,\partial s}\geq0,
\qquad
\frac{\partial^2\pi}{\partial a\,\partial d}\geq0,$$

then raising skill investment increases the marginal return to AI adoption, and AI adoption also increases the marginal return to decentralization.

::: {.theorem}
**Theorem 13 (Topkis's monotonicity theorem; Topkis 1978)** Let the choice set $X$ be a lattice (the componentwise maximum $x \vee x'$ and minimum $x \wedge x'$ of any two points are still in $X$), let the objective $f(x; \theta)$ be supermodular in $x$ ($f(x \vee x') + f(x \wedge x') \geq f(x) + f(x')$), and let it have increasing differences in $(x, \theta)$ (for $x' \geq x$, $f(x'; \theta) - f(x; \theta)$ is nondecreasing in $\theta$). Then for $\theta' \geq \theta$, if $x$ is optimal at $\theta$ and $x'$ is optimal at $\theta'$, then $x \vee x'$ is optimal at $\theta'$ and $x \wedge x'$ is optimal at $\theta$; hence the largest and smallest optimal choices (when they exist) are monotone nondecreasing in $\theta$.
:::

The proof sketch has two steps. Step one (increasing differences, then supermodularity): $f(x \vee x'; \theta') - f(x'; \theta') \geq f(x \vee x'; \theta) - f(x'; \theta) \geq f(x; \theta) - f(x \wedge x'; \theta) \geq 0$, the last inequality coming from the optimality of $x$ at $\theta$; so $x \vee x'$ is no worse than $x'$ and is optimal at $\theta'$; reading the same chain of inequalities backward, $x \wedge x'$ is optimal at $\theta$. Step two ($X$ a lattice): take $x$ and $x'$ to be the largest optimal choices at the two parameters; $x \vee x'$ is optimal at $\theta'$ and no smaller than $x'$, so it can only equal $x'$, hence $x' \geq x$; the smallest optimal choices are handled the same way.

In words: under the lattice and supermodular structure, a parameter that raises the marginal return to every practice at once pushes the optimal bundle of practices up as a whole; this requires no concavity of the profit function, no uniqueness of the optimum, and no differentiability, which makes it especially suited to organizational practices that come in bundles. As soon as one pair of choices is a pair of substitutes (a negative cross derivative), for example AI replacing low skills, monotonicity can fail. It also explains why "bought AI but productivity did not rise" cannot directly falsify the value of the technology: a firm that moves only $a$, without moving $s$ and $d$, may still be stuck in a low-complementarity system. Osprey's numbers: ordering S1 above S2, $x_m + 4q - 1$ is increasing in both $q$ and $x_m$, and Topkis's theorem implies that Osprey's best response is monotone in both: it chooses S1 when $x_m + 4q \geq 1$.

Empirically there are two common types of test. A correlation test checks whether complementary practices are adopted together, but common shocks also create positive correlation; a performance test checks whether joint adoption brings superadditive gains, but an interaction estimated from observational data is subject to selection. The most powerful design randomizes or quasi-randomizes the two practices separately. Writing $Y(a,s)$ for the firm's outcome (such as productivity) when the two practices are set at $a$ and $s$, such a design identifies the factorial contrast

$$\Delta_{as}=\{Y(1,1)-Y(0,1)\}-\{Y(1,0)-Y(0,0)\}.$$

$\Delta_{as}>0$ is the empirical counterpart of causal complementarity. If we can observe only the firms' own choices, we should at least report the bundles of practices they adopt and the timing of their adjustments, and state which selection restriction the identification of the interaction relies on. Chapter 34 applies this theory to AI, workflow redesign, and human capital.

Applying Topkis's theorem to each player's best response gives the skeleton of games with strategic complements.

::: {.theorem}
**Theorem 14 (smallest and largest equilibria of supermodular games; Milgrom and Roberts 1990; Vives 1990)** Let each $\mathcal{S}_i$ be a compact lattice, and let $u_i$ be supermodular in the player's own $s_i$ and have increasing differences in $(s_i, s_{-i})$ (for the semicontinuity conditions, see the original papers). Then (i) the set of pure-strategy Nash equilibria is nonempty, with a smallest equilibrium $s^{\min}$ and a largest equilibrium $s^{\max}$, and all other equilibria lie between the two; (ii) if $u_i$ also has increasing differences in $(s_i, \theta)$, both $s^{\min}$ and $s^{\max}$ are nondecreasing in $\theta$.
:::

The proof sketch has two steps. Step one (Theorem 13): each player's largest best response is increasing in the rivals' strategies, so together they form an increasing map on a lattice, and Tarski's fixed point theorem guarantees that its set of fixed points is nonempty and has a largest element (MWG Theorem M.I.3, p. 953, is only the special case of existence on $[0,1]^N$); the largest fixed point of the largest-best-response map is the largest equilibrium, and the smallest is handled the same way. Step two (Theorem 13): a rise in $\theta$ pushes the increasing map up pointwise, and its largest and smallest fixed points rise with it.

In words: in games with strategic complements, one can bound the whole equilibrium set from above and below, and know which way the two endpoints move, without knowing which equilibrium is actually played; equilibria in between may move the other way, so comparative statics speak only about the endpoints. Equilibrium selection remains open: an observed jump in adoption is consistent with strategic complements, but it could also come from a common demand shock. Osprey's numbers: in the compatibility game, (S2, S2) is the smallest equilibrium and (S1, S1) the largest, with the mixed equilibrium in between; adding G1's shift $x_m$, the two pure equilibria coexist when $-1 \leq x_m \leq 1$, and once $x_m$ exceeds USD 1 thousand, Osprey chooses S1 even when sure that Boreal is on S2, and the smallest equilibrium jumps to (S1, S1). The shifts in G1's high-precision block lie between $-0.2$ and $0.2$, so the two pure equilibria coexist everywhere.

### 3.11 Global Games: A Little Private Noise Selects a Unique Equilibrium

Merchants decide whether to join the new interconnection standard of Osprey and Boreal. The payoff from joining is $\theta + \ell - 1$, where $\ell$ is the share of merchants that join and $\theta$ is the standard's fundamental (in this subsection $\theta$ means only this), and not joining gives 0. Under complete information, as long as $\theta \in (0, 1)$, everyone joining (when $\ell = 1$ the payoff is $\theta > 0$) and no one joining (when $\ell = 0$ the payoff is $\theta - 1 < 0$) are both equilibria, and theory cannot say whether the new standard will tip (be adopted by everyone in a landslide) or be ignored. A global game (Carlsson and van Damme 1993; Morris and Shin 1998) mixes in a little incomplete information: merchant $i$ sees only a private signal $x_i = \theta + \sigma\varepsilon_i$, where $\varepsilon_i$ is independent standard normal and $\sigma$ is the scale of the private noise (not the mixed-strategy profile of Sections 3.2 and 3.6). The marginal merchant standing at the critical signal holds a posterior about "what share of merchants join" that is exactly uniform on $[0, 1]$; this is called the Laplacian belief, so $\mathbb{E}[\ell] = 1/2$, and the indifference condition $\theta^* + \tfrac{1}{2} - 1 = 0$ gives $\theta^* = 1/2$. Now add a public signal as well: the public prior of $\theta$ is normal, with mean $m_{\mathrm{pub}}$ and standard deviation $\sigma_{\mathrm{pub}}$.

::: {.theorem}
**Theorem 15 (unique threshold equilibrium of the global game; Carlsson and van Damme 1993; Morris and Shin 1998; the precision condition is a direct derivation for this chapter's setting)** In the setting above, let $m(x)$ be the posterior mean that the signal $x$ implies for $\theta$, and let $\Phi$ be the standard normal distribution function. (i) If $\sigma^2(\sigma^2 + \sigma_{\mathrm{pub}}^2) < 2\pi\,\sigma_{\mathrm{pub}}^4(\sigma^2 + 2\sigma_{\mathrm{pub}}^2)$, the threshold equilibrium is unique: a merchant joins when its signal is at least $x^*$, where $x^*$ is the unique root of $m(x^*) + \Phi\big((m(x^*) - x^*)/\sqrt{\sigma^2 + (\sigma_{\mathrm{pub}}^{-2} + \sigma^{-2})^{-1}}\big) - 1 = 0$, and when the fundamental is exactly $\theta^* = x^*$, half the merchants join; (ii) with $\sigma_{\mathrm{pub}}$ fixed and $\sigma \to 0$, $\theta^* \to 1/2$; (iii) without a public signal, for every $\sigma$, $\theta^* = 1/2$.
:::

The proof sketch has three steps; ruling out equilibria that are not of threshold form requires iterated deletion, and for the idea, see Carlsson and van Damme (1993) and Morris and Shin (1998). Step one (normal signals): when everyone else uses the threshold $x^*$, the share that joins at fundamental $\theta$ is $\Phi((\theta - x^*)/\sigma)$; the marginal merchant's posterior is normal with mean $m(x^*)$, and taking the expectation of the joining share gives the indifference equation; without a public signal $m(x) = x$, and the expected share is $\Phi(0) = 1/2$, which is (iii). Step two (the precision condition): the derivative of the left-hand side of the equation with respect to $x^*$ is no smaller than a quantity proportional to $\sigma^{-2} - \sigma_{\mathrm{pub}}^{-2}\big/\big(\sqrt{2\pi}\sqrt{\sigma^2 + (\sigma_{\mathrm{pub}}^{-2} + \sigma^{-2})^{-1}}\big)$; the condition in (i) makes it positive, so the left-hand side is strictly increasing and the root is unique. Step three ($\sigma_{\mathrm{pub}}$ fixed): $m(x^*) - x^*$ is a small quantity of order $\sigma^2$ and the denominator is of order $\sigma$, so the ratio tends to zero and the expected share tends to $1/2$, which is (ii).

In words: when private information is precise enough relative to public information, the global game gives a unique tipping prediction: whether the standard is adopted depends only on whether the fundamental is above the critical value, and the equilibrium selected is the risk-dominant one. Uniqueness has a precondition: when the public signal is too precise, multiple equilibria return, and the public signal moves the critical value away from 1/2. The conclusion is determined by the information structure, and a platform that discloses information is changing exactly this structure; testing it requires randomized information precision or disclosure rules. Osprey's numbers (block G4): with a public prior of mean 0.25 and standard deviation 0.5, the critical values for private noise $\sigma$ of 0.1, 0.05, and 0.01 are 0.5434, 0.5176, and 0.5030, approaching 1/2 from the pessimistic side; the sufficient condition for uniqueness holds when $\sigma < 0.7212$. If the public signal is as precise as the private one (both standard deviations 0.1, prior mean 0.5), three thresholds appear: $-0.49996$, 0.5, and 1.49996.

![Block G4: the share of 2,000 merchants that join at each fundamental. Lines are the equilibrium shares, points the simulated populations; the smaller the private noise, the more steeply the joining share jumps from 0 to 1 at the critical value, and the critical value approaches 1/2.](assets/fig/fig_02_global_game.svg)

## 4 Repeated Games: Cooperation, Punishment, and Price Wars

Section 1 left a jarring conclusion: the unique equilibrium of the pricing game is the mutually damaging (Low, Low), and high-price cooperation worth 20 in total cannot survive. But that is only the outcome of a one-shot encounter. Osprey and Boreal meet month after month, betrayal today invites retaliation from tomorrow on, and as long as the future matters enough, the immediate gain from a sneak attack cannot outweigh the later loss from retaliation. This section writes $R = 10$, $T = 14$, and $P = 6$ for the stage payoffs under cooperation, a unilateral sneak attack, and both pricing low (they are payoffs, not numbers of periods), and the two firms discount the future by a discount factor $\delta$.

### 4.1 Grim Trigger and the Folk Theorem

The grim trigger strategy works as follows: start with High, keep playing High as long as no one has ever cut price, and otherwise play Low forever. Keeping the deal earns $R$ each period, a total value of $R/(1-\delta)$; a sneak attack in some period earns $T$ in that period, and every period after that yields only $P$.

::: {.theorem}
**Theorem 16 (the cooperation threshold under grim trigger; computed directly from Theorem 6; for the Bertrand version, see MWG Proposition 12.D.1, p. 401)** Let the stage game satisfy $T > R > P$, with (Low, Low) a Nash equilibrium of it, and let the two players repeat it infinitely often (G.3). Grim trigger is an SPE if and only if

$$\frac{R}{1-\delta} \geq T + \frac{\delta P}{1-\delta}$$

Rearranging gives

$$\delta \geq \frac{T - R}{T - P}$$
:::

The proof sketch has two steps, using Theorem 6. Step one (the punishment state): both play the stage Nash equilibrium forever after, a one-shot deviation today changes only today's payoff, and Low is a best response to Low, so there is no profitable deviation. Step two (the cooperation state): the only one-shot deviation to check is cutting price today and being punished forever afterward according to the strategy; it is unprofitable if and only if the first inequality holds; multiplying both sides by $1 - \delta$ and rearranging gives the threshold.

In words: under G.3, whether permanent cooperation can be sustained depends only on whether patience exceeds the ratio of "the temptation to sneak attack, $T - R$," to "the strength of punishment, $T - P$." The punishment here already holds the rival to its minmax payoff of 6, and there is no harsher credible punishment, so when $\delta$ is below the threshold no SPE can have both firms price High in every period. The threshold is determined by market structure: $T - R$ is small when there are few players or when a sneak attack cannot grab most of the market, and $\delta$ is large when interaction is frequent, and either makes cooperation easier. Since $\delta$ is not seen in the data, what can be exploited is exogenous change in the duration of relationships, the frequency of interaction, or monitoring conditions. Osprey's numbers: the threshold is $(14 - 10)/(14 - 6) = 1/2$; at $\delta = 1/2$, the discounted value of keeping the deal is $R/(1-\delta) = 10/0.5 = 20$ and that of the sneak attack is $14 + 0.5\cdot 6/0.5 = 20$, the two being equal; at $\delta = 0.4$, 0.7, and 0.9 the slack in keeping the deal is $-1.3333$, 5.3333, and 32, respectively. In the Cournot version of block G2, the threshold for splitting the monopoly output under Cournot punishment is 0.5294.

The general form of this conclusion is the folk theorem. Write $u(s)$ for the stage payoff vector of a strategy profile $s$; the feasible payoff set $\mathcal{F} = \mathrm{conv}\{u(s) : s \in \mathcal{S}_1 \times \dots \times \mathcal{S}_n\}$ is the convex hull of these vectors, and player $i$'s minmax payoff $\underline u_i = \min_{\sigma_{-i}} \max_{s_i} u_i(s_i, \sigma_{-i})$ is the lowest value to which the rivals together can hold it.

::: {.theorem}
**Theorem 17 (folk theorem; Fudenberg and Maskin 1986; for the two-player case, see MWG Proposition 12.AA.5, p. 422, and for the dimension condition with more players, its footnote 34)** In an infinitely repeated game under G.1(a) and G.3, let $u^{\circ} \in \mathcal{F}$, and for every $i$ let $u^{\circ}_i > \underline u_i$. With two players, there is a threshold below 1 such that, for discount factors above it, $u^{\circ}$ is the average payoff of some SPE; for $n \geq 3$, if the dimension of $\mathcal{F}$ equals $n$, the same conclusion holds.
:::

For the proof, see Fudenberg and Maskin (1986). This book only states the result, because the case of more players requires constructing strategies under which the punishers themselves are willing to carry out the punishment: after the punishment phase, different players receive different rewards, the dimension condition guarantees that such rewards exist, and the construction is long. With two players the idea fits in one sentence: on the path, play the outcome that gives $u^{\circ}$, and if someone deviates, minmax each other for a number of periods and then return; a patient deviator finds the brief sneak attack not worth it, and the punishers are willing to sit through the punishment phase because the rewards afterward are higher.

In words: under G.3 (with three or more players, plus the dimension condition), sufficiently patient players can sustain as an SPE any payoff profile that is feasible and strictly above the minmax, so the inefficient outcome that is unavoidable in the one-shot game can be overturned. The cost is severe equilibrium multiplicity: any outcome from mutual damage to full collusion can be an equilibrium; theory can say what can be sustained but not which one will happen, and the outcome depends on expectations, conventions, and communication, things outside the model. For empirical work, a high price might be collusion or might be some other equilibrium, and equilibrium concepts alone cannot tell them apart; Section 3.4 of Chapter 29 discusses what variation a test of conduct requires. Osprey's numbers: the feasible set is the convex hull of (10, 10), (2, 14), (14, 2), and (6, 6), with area 48; both firms' minmax is 6, and the part that is feasible and on or above both minmax lines has vertices (6, 6), (12, 6), (10, 10), and (6, 12), with area 24.

![The folk theorem in the pricing game: the light region is the convex hull of feasible payoffs, and the dashed lines are the two firms' minmax levels of 6 each; the dark region holds the feasible payoffs above both minmax lines, all of which can be sustained by an SPE when players are patient enough.](assets/fig/fig_02_folk.svg)

### 4.2 Imperfect Public Monitoring: Price Wars on the Equilibrium Path

Grim trigger assumes that a betrayal is always seen. In reality firms often see only the market price, not the rival's output; a low price may mean that someone secretly produced more, or just that demand was weak. In each local market of block G2, the price in month $t$ is $p_t = 10 - Q_t + \sigma z_t$ dollars, where $Q_t$ is the two firms' total output (thousand units), the demand shock $z_t$ is an independent standard normal draw that the firms do not see, and a price war lasts 6 months; in this subsection $p_t$ is the price, $q$ is output, and $\sigma$ is the demand noise.

::: {.theorem}
**Theorem 18 (trigger-price equilibrium under imperfect public monitoring; Green and Porter 1984, equations (2) to (6) on pp. 91-93; the original has multiplicative shocks with a general distribution, written here with additive shocks)** Let the price of a homogeneous-goods oligopoly be $p_t = p(Q_t) + \sigma z_t$, where $p(\cdot)$ is inverse demand and $z_t$ is i.i.d. with a continuous density; firms see only the price and their own output and discount by $\delta$ (G.3(b)). Trigger-price strategy: in normal periods each produces $q_{\mathrm{col}}$; once the price falls below the trigger price $\underline{p}$, each produces the stage-game Cournot output for the next $N_{\mathrm{pun}}$ periods and then returns to normal. Let $\zeta(q)$ be the probability, when a firm produces $q$ and the rival $q_{\mathrm{col}}$, that the price is below $\underline{p}$; let $\zeta = \zeta(q_{\mathrm{col}})$, let $\pi(q, q_{\mathrm{col}})$ be the current-period profit, and let $V_{\mathrm{c}}$ and $V_{\mathrm{w}}$ be the discounted values in a normal period and at the start of a price war. (i) These strategies form an equilibrium if and only if in a normal period no output $q$ gives $\pi(q, q_{\mathrm{col}}) + \delta\,[(1 - \zeta(q)) V_{\mathrm{c}} + \zeta(q) V_{\mathrm{w}}] > V_{\mathrm{c}}$; the first-order condition at an interior point is that the marginal profit from producing one more unit exactly equals the increase in the probability of a price war it causes times $\delta(V_{\mathrm{c}} - V_{\mathrm{w}})$. (ii) On the equilibrium path no firm deviates, yet a price war breaks out in each normal period with probability $\zeta$, and the long-run fraction of time spent in price wars is $\zeta N_{\mathrm{pun}}/(1 + \zeta N_{\mathrm{pun}})$.
:::

The proof sketch has three steps. Step one (stationarity): the rival's actions depend only on the public price, so each firm faces a stationary dynamic program with only two public states, "normal" and "price war" (original p. 92), and $V_{\mathrm{c}}$ is given by the recursion of the original's equations (2) and (3). Step two (the same argument as for Theorem 6, applied to the public states): the continuation after a price-war period does not depend on today's output, so the stage Cournot output is optimal; in a normal period only a one-shot deviation needs checking, which is (i), corresponding to the original's equations (4) to (6). Step three ($z_t$ i.i.d.): on the path everyone produces $q_{\mathrm{col}}$, so a price below the trigger can only come from a demand shock, with probability $\zeta$ (original p. 93); the normal state and each period of a price war form a Markov chain, and its stationary distribution gives (ii).

In words: in the Green-Porter equilibrium, price wars are part of collusion, not its breakdown: a low price is triggered by bad demand, and the firms know full well that no one betrayed, yet they must punish as agreed; otherwise the temptation to produce more in normal periods could not be held in check. With perfect monitoring price wars never happen; with noise, how often they happen depends on which equilibrium is selected: block G2 takes the trigger price that maximizes the value of collusion, and the rise in price-war frequency with noise is a property of this choice; with another trigger price or punishment length, the frequency differs. Abreu, Pearce and Stacchetti (1990) characterize the equilibrium payoff set of this class of games. In block G2, the larger the noise, the more often price wars break out and the smaller the loss from a price war, $V_{\mathrm{c}} - V_{\mathrm{w}}$ (at the margin, the effect of a deviation on the probability of a price war in fact rises slightly), so punishment deters less and the colluders can only produce more in normal periods, which shrinks the gain from cheating. Footnote 9 of Green and Porter notes that the optimal output differs from the monopoly output. On p. 95 the original says that the trigger price and the punishment length can be estimated consistently from the price series, and Chapter 4 reads the estimates of Porter (1983) for a railroad freight cartel. Osprey's numbers: for demand noise $\sigma$ of USD 0, 0.1, 0.2, and 0.4, the optimal trigger prices are USD 6, 5.7203, 5.4854, and 5.0169, the probabilities that a price war breaks out in a normal month are 0, 0.0191, 0.0465, and 0.1298, and the long-run fractions of time in price wars are 0, 0.1026, 0.2181, and 0.4377; the values of collusion are 8, 7.9198, 7.8185, and 7.5632 thousand USD a month, all above the Cournot value of 7.1111.

![Block G2: the left panel shows the fraction of market-months in a price war as demand noise varies, the line being the stationary fraction in equilibrium and the points the simulation of 200 markets over 300 months (95% intervals); the right panel shows the trigger price, the expected price in collusive periods, and the Cournot price.](assets/fig/fig_02_price_wars.svg)

When the discount factor is lowered from 0.95 to 0.5, all four strategies chosen for 0.95 fail the one-shot deviation check; after re-optimizing the trigger price, collusion can still be sustained at noise 0, 0.1, and 0.2, with values of 7.9960, 7.8657, and 7.4860 thousand USD a month, between the Cournot value and the value under patience, while at noise 0.4 no collusive equilibrium exists within this family of trigger strategies, and the value falls back to the Cournot value of 7.1111. Within this family, when patience is low, collusion holds up only when monitoring noise is small.

### 4.3 Comparative Statics, Equilibrium Selection, and Empirical Implications

The difficulty with comparative statics in game theory is that the equilibrium itself changes, and there may be more than one equilibrium. The table below writes several of the chapter's comparative statics as predictions that can be taken to data; the last column gives the extra handle needed to identify each.

| Change in theoretical parameter | Equilibrium prediction | What can be observed empirically | Extra handle needed for identification |
|---|---|---|---|
| Complementarity rises | Best responses steeper, multiple equilibria and tipping more likely | Adoption clusters among connected agents, shocks have network multipliers | Exogenous network edges, a variable that moves only one side's payoff |
| Discount factor $\delta$ rises | Cooperation incentive strengthens (Theorem 16) | Fewer punishment periods where collusion is already sustainable, but price wars appear from none once the threshold is crossed (block G2) | Exogenous relationship duration or monitoring change |
| Demand noise rises | More frequent price wars under the value-maximizing trigger price (block G2) | Low-price periods follow one bad demand shock and last a fixed number of periods, and prices recover afterward | Observable demand-shock variables, monitoring information beyond prices |
| Private-noise variance falls | Global-game actions closer to the common threshold (Theorem 15) | Actions more concentrated in the critical region of the fundamental | Randomized information precision or disclosure rule |
| First-mover commitment strengthens | SPE may differ from the simultaneous-move Nash equilibrium (Theorem 5) | After the timing changes, follower responses vary systematically | Exogenous rollout order or institutional timing |

When testing whether theory and data agree, the equilibrium selection rule should be specified in advance. If a parameter change moves several equilibria at once, the direction of the observed outcome also depends on the selection rule; calibrating some Nash equilibrium to the data does not mean that this equilibrium is identified, and at a minimum one should report whether other equilibria give the same observable predictions. Observed price parallelism is not collusion either: co-moving costs and shared algorithms can produce the same pattern.

## 5 Anchoring Papers

Three studies mark the chapter's three landing points: one brings the compatibility game into information-technology markets, one reads price wars as part of collusion, and one fits a game of incomplete information to store data.

### 5.1 Katz and Shapiro (1985)

::: {.case}
Paper: "Network Externalities, Competition, and Compatibility," *American Economic Review*. It is the real-world prototype of the chapter's compatibility game, and it started the tradition of analyzing technology standards and network effects with game theory.

Method: consumers' valuation of a product rises with the number of users of the same technology (a network externality), and whether firms make their products compatible is a strategic choice; the paper compares market equilibria with and without compatibility. Compatibility merges the firms' user bases into one large network and amplifies network payoffs, but it also weakens the advantage a large firm gets from its existing user base.

Result: firms' private incentives for compatibility often diverge from the social optimum. Firms with a large user base, or ones optimistic about their future share, tend to refuse compatibility to protect their moat, while small firms and late entrants want compatibility more; the market may settle in a less efficient incompatible configuration, or fall into the kind of coordination failure described in Section 3.2.

Limitations: the static framework simplifies the formation of expectations and the dynamics of preemption; later work on the installed base, dynamic standards competition, and two-sided platforms started precisely from relaxing these simplifications.
:::

### 5.2 Green and Porter (1984)

::: {.case}
Paper: "Noncooperative Collusion under Imperfect Price Information," *Econometrica*. The question is how collusion enforces itself when demand is uncertain, and what price wars in the data mean.

Method: the model of Theorem 18. Firms see only a market price subject to multiplicative demand shocks, not the rivals' output; after the price falls below the trigger, they revert to Cournot for a number of periods and then resume collusion (pp. 91-92). The incentive constraint and first-order condition for normal periods are the original's equations (4) to (6); in equilibrium no firm deviates, and the frequency of price wars is given by the trigger probability (p. 93); in footnote 10 the authors point out that this equilibrium is sequentially rational.

Data: this is a theory paper. The authors list four industry characteristics under which the model applies: long-run stability, output as the only decision variable, public information about everything except one's own output, and monitoring information only imperfectly correlated with firms' behavior; they consider that the US railroad freight industry in the 1880s roughly fits them (p. 90, p. 94).

Result: it reverses the traditional reading of price wars. A sharp price drop that lasts for a while and then recovers used to be read as entry deterrence, or as a cartel breaking down and being rebuilt; in this model it is a necessary part of a stable collusive mechanism (p. 94). The price process is stationary and ergodic, and the trigger price and punishment length can be estimated consistently from the price series (p. 95); the authors cite the econometric work of Porter (1983) as support (p. 94).

Limitations: the authors themselves say that the assumptions on industry structure are restrictive, that what they study is only the simplest kind of trigger strategy (footnote 8), and that it is not the only collusive equilibrium either (p. 94); the i.i.d. shock specification is inappropriate for durable goods or for consumption that can be shifted across periods (footnote 7).
:::

### 5.3 Ellickson and Misra (2008)

::: {.case}
Paper: "Supermarket Pricing Strategies," *Marketing Science*. The question is whether supermarkets, in choosing their positioning between EDLP (everyday low pricing) and promotional pricing, respond to their competitors' choices.

Method: each store's choice of pricing strategy is written as a static discrete game of incomplete information, the BNE of Section 3.4: a store does not know its rivals' private shocks, forms probabilistic expectations of their choices, and chooses accordingly; the authors estimate it as a system of simultaneous discrete choice models. Section 3.3 of Chapter 30 shows that incomplete information can restore a unique symmetric equilibrium under a competition effect; Ellickson and Misra, however, estimate strategic complements, where uniqueness is not guaranteed.

Data: store-level data on supermarkets operating in the United States in 1998, with each store's pricing strategy reported by the firm itself.

Result: contrary to several theoretical predictions, stores tend to choose the same strategy as their rivals in the same market, so strategies cluster; demographics and store and chain characteristics also matter. This is exactly what a game model fitted to data looks like: strategic complementarity is estimated in the form of "how the expectation of rivals' choices moves one's own choice probabilities."

Limitations: the static game abstracts from the dynamic adjustment of pricing strategies; the strategies are self-reported categories.
:::

Taken together, the three papers show that for a game model to meet data, one must first say how equilibrium determines the observed outcomes: Katz and Shapiro give a source of multiple equilibria, Green and Porter tell us what price series a given equilibrium generates, and Ellickson and Misra use incomplete information to write choices as probabilities, so that they can be estimated.

## 6 Running Case: Settling the Two Platforms' Accounts

The code below runs under R 4.6.1. Only the key statements are shown, with several statements sometimes merged onto one line; every number in the text comes from the output of an actual run, and the result lines in the code that begin with `#` are copied verbatim from what the run printed.

### 6.1 G0: Pricing, Compatibility, Sequential Play, and BNE by Hand

`nash_2x2` checks pure-strategy equilibria cell by cell and solves for the interior mixed equilibrium from the two indifference conditions; in a degenerate game there can be a whole continuum of mixed equilibria, which it does not report, and none of the chapter's games is degenerate. The last line of output covers the sequential game and the BNE: Boreal always matches, and Osprey gets 3 from S1 and 1 from S2; in the BNE, Osprey gets 7.2 from High and 6.6 from Low.

```r
source("games_fns.R")                                   # nash_2x2, minmax_2x2, ce_violation, rubinstein_share, ...
A_pd <- matrix(c(10, 2, 14, 6), 2, byrow = TRUE)        # pricing game, [Osprey's action, Boreal's action]: High, Low
B_pd <- t(A_pd)                                         # symmetric game: Boreal's payoffs are the transpose
A_cp <- matrix(c(3, 0, 0, 1), 2, byrow = TRUE)          # compatibility game: S1, S2
B_cp <- matrix(c(1, 0, 0, 3), 2, byrow = TRUE)
ne_pd <- nash_2x2(A_pd, B_pd); ne_cp <- nash_2x2(A_cp, B_cp)
p <- ne_cp$mixed[1]; q <- ne_cp$mixed[2]; joint <- outer(c(p, 1 - p), c(q, 1 - q))
c(high_dominated = all(A_pd[2, ] > A_pd[1, ]), ne_pd = ne_pd$pure[[1]], minmax = minmax_2x2(A_pd),
  n_pure_cp = length(ne_cp$pure), p = p, q = q, coordinate = joint[1, 1] + joint[2, 2], eu = sum(joint * A_cp))
follow <- c(which.max(B_cp[1, ]), which.max(B_cp[2, ]))  # Boreal's reply to S1 and to S2: it matches
A_bne <- matrix(c(10, 3, 7, 6), 2, byrow = TRUE)          # Osprey's payoffs: rows High, Low; columns Boreal High, Low
beta <- 2 / 5; eu <- drop(A_bne %*% c(1 - beta, beta))    # Boreal is Normal (plays Low) with probability beta
c(follow = follow, osprey_S1 = A_cp[1, follow[1]], osprey_S2 = A_cp[2, follow[2]], eu_high = eu[1], eu_low = eu[2])
# high_dominated         ne_pd1         ne_pd2         minmax      n_pure_cp
#          1.000          2.000          2.000          6.000          2.000
#              p              q     coordinate             eu
#          0.750          0.250          0.375          0.750
#   follow1   follow2 osprey_S1 osprey_S2   eu_high    eu_low
#       1.0       2.0       3.0       1.0       7.2       6.6
```

### 6.2 G0: Repeated Games, Bargaining, and Correlated Equilibrium

```r
R <- 10; Tm <- 14; P <- 6                                # cooperation, temptation, punishment (million USD a month)
coop <- function(d) R / (1 - d); cheat <- function(d) Tm + d * P / (1 - d)
round(c(delta_star = (Tm - R) / (Tm - P), slack = coop(c(0.4, 0.5, 0.7, 0.9)) - cheat(c(0.4, 0.5, 0.7, 0.9))), 4)
V <- 100; nb <- function(d1, d2, alpha = 0.5) d1 + alpha * (V - d1 - d2)   # Osprey's Nash bargaining share
c(sym = nb(0, 0), asym = nb(20, 0), weighted = nb(20, 0, 0.6))
round(V * rubinstein_share(c(0.8, 0.99), c(0.8, 0.99)), 4)           # first proposer's share, equal patience
dt <- c(year = 1, month = 1 / 12, week = 1 / 52, day = 1 / 365)       # time between offers, in years
round(rbind(osprey_first = rubinstein_share(exp(-0.10 * dt), exp(-0.05 * dt)),
            boreal_first = 1 - rubinstein_share(exp(-0.05 * dt), exp(-0.10 * dt))), 5)   # Osprey's share
coin <- matrix(c(0.5, 0, 0, 0.5), 2)                                  # public coin: (S1, S1) or (S2, S2)
c(ce_osprey = sum(coin * A_cp), ce_boreal = sum(coin * B_cp), violation = ce_violation(coin, A_cp, B_cp),
  gg_flat_prior = gg_threshold(0.1))
# delta_star     slack1     slack2     slack3     slack4
#     0.5000    -1.3333     0.0000     5.3333    32.0000
#      sym     asym weighted
#       50       60       68
# [1] 55.5556 50.2513
#                 year   month    week     day
# osprey_first 0.35013 0.33472 0.33365 0.33338
# boreal_first 0.31681 0.33195 0.33301 0.33329
#     ce_osprey     ce_boreal     violation gg_flat_prior
#           2.0           2.0           0.0           0.5
```

The Rubinstein share converges quickly as the interval between offers shrinks: with monthly offers, Osprey gets 0.3347 when it proposes first, already close to 1/3.

### 6.3 G1: QRE Markets and the Log-Odds Inversion

G1 has 4,000 local markets. The $\chi = 0.8$ block has 2,000 markets, with shifts $x_m$ from $-1.2$ to 1.2 thousand USD, 80 markets at each shift value, each with a unique QRE; the $\chi = 2$ block has 2,000 markets, with shifts between $-0.2$ and 0.2 and three QRE at each shift value; each market plays the QRE that favors the standard that entered the local market first, and which standard entered first (the history) is drawn according to a rule that varies with $x_m$. In each market each of the two firms makes 50 decisions. The analyst knows $\chi$ and $x_m$, and by Theorem 10 (ii), the frequencies $p$ and $q$ with which the two firms choose S1 in a group of markets playing the same QRE satisfy

$$\frac{1}{\chi}\ln\frac{p}{1-p} - x_m = g_{\mathrm{S1}}\, q - g_{\mathrm{S2}}\,(1 - q)$$

Here $g_{\mathrm{S1}}$ is the extra payoff Osprey gets from S1 over S2 when Boreal is on S1 (true value 3), and $g_{\mathrm{S2}}$ is the extra payoff from S2 when Boreal is on S2 (true value 1), in thousand USD a month. OLS without an intercept, run on the groups' frequencies, gives both. There are three groupings: the unique-QRE block in 25 groups by $x_m$; the high-precision block in only 5 groups by $x_m$, each group mixing markets that play different QRE; and the high-precision block in 10 groups by $x_m$ and history, that is, grouped by the equilibrium played. Sampling uncertainty is reported through a simulation that regenerates the data 1,000 times.

```r
source("games_fns.R"); tr <- readRDS("truth02.rds"); P1 <- tr$params$G1; G1 <- readRDS("osprey_G1.rds")
round(rbind(qre_all(P1$A0, P1$B0, 0.8), qre_all(P1$A0, P1$B0, 2)), 4)   # every QRE at x_m = 0: chi 0.8, then chi 2
designs <- list(unique = list(block = "unique", by = "x"), pooled = list(block = "high", by = "x"),
                by_market = list(block = "high", by = c("x", "history")))   # by_market: grouped by the QRE played
groups_of <- function(d, k) {                            # group frequencies p, q of S1, with the known chi
  g <- g1_group_freq(d[d$block == designs[[k]]$block, ], designs[[k]]$by)
  g$chi <- P1$lambda[[designs[[k]]$block]]; g
}
sims <- setNames(lapply(designs, function(z) vector("list", 1000)), names(designs))
set.seed(2102)                                           # Monte Carlo: 1,000 fresh draws of block G1
for (r in 1:1000) {
  d <- gen_G1(P1, tr$G1$eq)
  for (k in names(designs)) { g <- groups_of(d, k); g$rep <- r; sims[[k]][[r]] <- g }
}
round(t(sapply(names(designs), function(k) {             # OLS of log(p/(1-p))/chi - x on q and -(1-q)
  mc <- qre_logodds_ols(do.call(rbind, sims[[k]]))
  c(draw = qre_logodds_ols(groups_of(G1, k))$coef[["u11"]], mc_mean = mc$coef[["u11"]], mc$ci, u22 = mc$coef[["u22"]])
})), 4)
#           p      q
# [1,] 0.6101 0.3899
# [2,] 0.1249 0.0067
# [3,] 0.6645 0.3355
# [4,] 0.9933 0.8751
#             draw mc_mean     lo     hi    u22
# unique    2.9946  3.0009 2.9399 3.0656 1.0008
# pooled    1.0648  1.0908 1.0257 1.1548 0.6262
# by_market 3.0207  3.0055 2.9358 3.0744 1.0005
```

In the code, `u11` and `u22` are $g_{\mathrm{S1}}$ and $g_{\mathrm{S2}}$. In the unique-QRE block, the simulation mean of $g_{\mathrm{S1}}$ is 3.0009, and the middle 95% of the 1,000 estimates, $[2.9399, 3.0656]$, covers the truth; grouped by the equilibrium played, the mean is 3.0055 and the range $[2.9358, 3.0744]$ covers it as well, so multiple equilibria by themselves do not break the inversion. When markets that play different equilibria are pooled into one group, the mean is only 1.0908, biased down 63.6%, the range $[1.0257, 1.1548]$ is far from 3, and $g_{\mathrm{S2}}$ also drops from 1 to 0.6262. The reason is that log odds are nonlinear in frequencies: each group mixes markets where S1 entered first, in which Osprey almost always chooses S1, with markets where S2 entered first, in which it rarely does, so the average frequencies no longer satisfy the logit equation; the probability that S1 entered first also rises with $x_m$, from about 17% to 83%, and the slope is flattened. The direction of the bias depends on where the two equilibria sit in this design; with another mixing proportion the bias need not be downward.

### 6.4 G2: The Price-War Panel

```r
source("games_fns.R"); tr <- readRDS("truth02.rds"); G2 <- readRDS("osprey_G2.rds")   # 200 markets x 300 months
round(t(sapply(seq_along(tr$G2$eq), function(l) {
  e <- tr$G2$eq[[l]]; v <- gp_values(e); d <- G2[G2$level == l, ]   # e: trigger and q_c chosen in 01_dgp.R
  c(sigma = e$sigma, trigger = e$trigger, q_c = e$q_c, beta = v$beta, war_share = v$beta * 6 / (1 + v$beta * 6),
    war_sim = mean(d$war), value = (1 - 0.95) * v$Vc, osd = osd_check(e))
})), 4)
low <- lapply(c(0, 0.1, 0.2, 0.4), function(s)               # re-optimise the trigger for impatient firms
  gp_optimal_trigger(list(a = 10, cost = 2, delta = 0.5, n_pun = 6L, sigma = s)))
round(sapply(low, function(e) (1 - 0.5) * gp_values(e)$Vc), 4)   # value a month at delta = 0.5; Cournot 7.1111
#      sigma trigger    q_c   beta war_share war_sim  value osd
# [1,]   0.0  6.0000 2.0000 0.0000    0.0000  0.0000 8.0000   1
# [2,]   0.1  5.7203 2.0362 0.0191    0.1026  0.1010 7.9198   1
# [3,]   0.2  5.4854 2.0893 0.0465    0.2181  0.2177 7.8185   1
# [4,]   0.4  5.0169 2.2661 0.1298    0.4377  0.4352 7.5632   1
# [1] 7.9960 7.8657 7.4860 7.1111
```

The `beta` column of the output is $\zeta$, and the `q_c` column is $q_{\mathrm{col}}$. At each noise level, `gp_optimal_trigger` picks, among the equilibria that satisfy the incentive constraint, the trigger price with the largest value of collusion, and `osd_check` uses Theorem 6 to verify that neither state has a profitable one-shot deviation. No firm deviates in the simulated panel, and the fractions of months in price wars, 0.101, 0.2177, and 0.4352, agree with the stationary fractions; the larger the noise, the more each colluder produces in normal periods, rising from the monopoly level of 2 thousand units to 2.0362, 2.0893, and 2.2661.

### 6.5 G3: Regret Matching and Fictitious Play

```r
source("games_fns.R"); P3 <- readRDS("truth02.rds")$params$G3; G3 <- readRDS("osprey_G3.rds")
mu <- regret_matching(P3$chicken_A, P3$chicken_B, T = 50000, seed = 2105)   # Chicken; actions Yield, Enter
attributes(mu) <- list(dim = c(2L, 2L)); mu                                  # empirical joint distribution
c(chicken = ce_violation(mu, P3$chicken_A, P3$chicken_B),                    # largest gain from disobeying
  compat = ce_violation(matrix(G3$rm_compat, 2), P3$compat_A, P3$compat_B))
path <- attr(G3$reps_chicken, "path"); keep <- P3$checkpoints %in% c(100, 1000, 10000, 50000)
round(rbind(mean = colMeans(path), max = apply(path, 2, max))[, keep], 6)    # 100 independent Chicken runs
fp <- fictitious_play(P3$mp_A, P3$mp_B, T = 50000, seed = 2106); fp$freq      # matching pennies
#         [,1] [,2]
# [1,] 0.00018    0
# [2,] 0.99982    0
# chicken  compat
# 0.00018 0.00000
#        [,1]    [,2]     [,3]    [,4]
# mean 0.0331 0.00331 0.000331 6.6e-05
# max  0.3900 0.03900 0.003900 7.8e-04
#       p       q
# 0.49720 0.50028
```

This Chicken run locks into (Enter, Yield) from period 10 on, and the compatibility-game run stays at (S1, S1) from period 1. The average gain from disobeying over the 100 runs is 0.0331, 0.00331, and 0.000331 at periods 100, 1,000, and 10,000, falling roughly like $1/t$. Fictitious play on matching pennies, which has no pure equilibrium, gives empirical frequencies of 0.4972 and 0.5003, close to the 1/2 of the mixed equilibrium.

### 6.6 G4: The Critical Value of the Global Game

```r
source("games_fns.R"); sig <- c(0.1, 0.05, 0.01)                        # private noise
x_star <- gg_threshold(sig, prior_mean = 0.25, prior_sd = 0.5)          # public signal theta ~ N(0.25, 0.5^2)
pop <- readRDS("osprey_G4.rds")$pop                                     # 2,000 merchants at each theta
cross <- sapply(sig, function(s) { d <- pop[pop$sigma == s, ]; sh <- d$n_adopt / d$n_pop
  i <- which(sh[-length(sh)] < 0.5 & sh[-1] >= 0.5)[1]                  # where the simulated share passes 1/2
  d$theta[i] + (0.5 - sh[i]) * (d$theta[i + 1] - d$theta[i]) / (sh[i + 1] - sh[i]) })
round(rbind(theta_star = x_star, simulated = cross, min_slope = gg_min_slope(sig, 0.5)), 4)
round(c(sigma_max = gg_sigma_max(0.5), roots = gg_roots(0.1, prior_mean = 0.5, prior_sd = 0.1)), 5)
#              [,1]   [,2]   [,3]
# theta_star 0.5434 0.5176 0.5030
# simulated  0.5427 0.5164 0.5026
# min_slope  0.8520 0.9341 0.9883
# sigma_max    roots1    roots2    roots3
#   0.72123  -0.49996   0.50000   1.49996
```

The joining share of the simulated populations crosses one half at 0.5427, 0.5164, and 0.5026, almost exactly at the solved critical values; the minimum slopes of the indifference equation are all positive, which is exactly the uniqueness condition of Theorem 15.

## 7 Failure Modes and Robustness

The first failure mode is equilibrium multiplicity. The compatibility game has three Nash equilibria, the folk theorem gives a whole region of sustainable payoffs, and the theory itself often cannot select which one will be realized. There are several ways to cope: shrink the equilibrium set with refinements such as SPE (Sections 3.5 and 3.6); bring in selection principles from outside the model, such as focal points or history (an installed base, for example); or acknowledge the multiplicity and treat it as an incompleteness of the model, which Chapter 30 handles through partial identification. The tools of Sections 3.7 and 3.11 work in opposite directions: the global game collapses a region of coordination equilibria into a unique critical value, and multiplicity returns as soon as the public signal is precise; correlated equilibrium spreads the equilibrium set wider, so explanatory power grows and predictive power shrinks.

The second is treating equilibrium selection as known. The inversion of Section 6.3 hits the truth when markets are grouped by the equilibrium played, and when they are pooled it falls 63.6% short; Bajari, Hong and Ryan (2010) estimate the equilibrium selection rule together with the payoff parameters, which is another route. On seeing a price war, one should likewise ask first which equilibrium is being played: in the Green-Porter equilibrium, price wars are part of collusion.

The third is the strength of the rationality assumption. Nash equilibrium assumes that players are not only rational but also correctly anticipate each other's equilibrium strategies. Experiments repeatedly show that people's depth of reasoning is limited; level-k and cognitive hierarchy models (Nagel 1995; Camerer, Ho and Chong 2004) write this into the model, while QRE writes mistakes as noise. Game theory is better used as a benchmark for "what fully rational players who fully understand the structure would do," and the direction in which reality deviates from the benchmark is itself often informative.

The fourth is the interpretation of mixed strategies. Platforms rarely roll dice to choose a standard: a better reading of a mixed equilibrium is a population-level frequency, or the rival's uncertainty about your choice, which is the purification of Harsanyi (1973) (Section 3.4). When taking a mixed equilibrium to data, say which reading you mean.

The fifth is credibility and commitment. SPE eliminates empty threats, but irreversible investment, public commitment, and reputation can all lock in a threat that would otherwise be empty; when analyzing dynamic competition, what matters is often whether such a commitment device exists.

The predictive power of game theory depends heavily on the situation: it is powerful in markets with a clear structure, few players, and repeated encounters, while in settings with multiple equilibria, bounded rationality, and variable commitment its conclusions should be treated as a conditional benchmark, and one should state which way they move when the assumptions loosen.

## 8 Further Reading

::: {.readings}
Required reading, in the suggested reading order:

- Mas-Colell, Whinston and Green (1995, Microeconomic Theory), Chapters 7 to 9. The authoritative reference for game theory: its Chapter 8 (Nash equilibrium and BNE) and Chapter 9 (SPE and sequential equilibrium) are the sources of most of the theorems in Sections 3.1 to 3.6 of this chapter.
- Gibbons (1992, Game Theory for Applied Economists). An application-oriented introduction, dense with examples, good to read alongside MWG.

Further reading:

- Fudenberg and Tirole (1991, Game Theory). The standard graduate-level reference, especially complete on dynamic games, refinements, and repeated games.
- Fudenberg and Maskin (1986). The rigorous statement of the folk theorem for discounted repeated games, Theorem 17.
- Green and Porter (1984) and Abreu, Pearce and Stacchetti (1990). Collusion under imperfect monitoring and its equilibrium payoff set.
- Katz and Shapiro (1985). The origin of the game-theoretic analysis of network externalities and compatibility, required reading for work on information systems and platforms.
- Ellickson and Misra (2008) and Bajari, Hong and Ryan (2010). Two ways of fitting discrete games to data: incomplete information, and estimation together with equilibrium selection.
- Harsanyi (1967, 1973), Kreps and Wilson (1982), and Selten (1975). The type method, purification, sequential equilibrium, and trembling-hand perfection.
- Debreu (1952), Glicksberg (1952), and Fan (1952). The original papers on the existence of pure-strategy equilibrium in continuous games, Theorem 3.
- Nash (1950a, 1950b), Rubinstein (1982), Binmore, Rubinstein and Wolinsky (1986), and Osborne and Rubinstein (1990). Existence, the axiomatic and strategic solutions of bargaining, and the Nash program that connects the two.
- Aumann (1974, 1987) and Hart and Mas-Colell (2000). Correlated equilibrium and a learning rule that leads to it.
- McKelvey and Palfrey (1995). The original QRE paper, the bridge from games to choice data.
- Topkis (1978), Milgrom and Roberts (1990), and Vives (1990). Monotone comparative statics and supermodular games.
- Carlsson and van Damme (1993) and Morris and Shin (1998). The founding paper on global games and its classic application to self-fulfilling crises.
- Nagel (1995) and Camerer, Ho and Chong (2004). Experimental evidence on, and models of, limited depth of reasoning.
:::

::: {.apa-refs}
- Abreu, D., Pearce, D., & Stacchetti, E. (1990). Toward a theory of discounted repeated games with imperfect monitoring. *Econometrica, 58*(5), 1041-1063. https://doi.org/10.2307/2938299
- Aumann, R. J. (1974). Subjectivity and correlation in randomized strategies. *Journal of Mathematical Economics, 1*(1), 67-96. https://doi.org/10.1016/0304-4068(74)90037-8
- Aumann, R. J. (1987). Correlated equilibrium as an expression of Bayesian rationality. *Econometrica, 55*(1), 1-18. https://doi.org/10.2307/1911154
- Bajari, P., Hong, H., & Ryan, S. P. (2010). Identification and estimation of a discrete game of complete information. *Econometrica, 78*(5), 1529-1568. https://doi.org/10.3982/ecta5434
- Binmore, K., Rubinstein, A., & Wolinsky, A. (1986). The Nash bargaining solution in economic modelling. *The RAND Journal of Economics, 17*(2), 176-188. https://doi.org/10.2307/2555382
- Camerer, C. F., Ho, T.-H., & Chong, J.-K. (2004). A cognitive hierarchy model of games. *The Quarterly Journal of Economics, 119*(3), 861-898. https://doi.org/10.1162/0033553041502225
- Carlsson, H., & van Damme, E. (1993). Global games and equilibrium selection. *Econometrica, 61*(5), 989-1018. https://doi.org/10.2307/2951491
- Debreu, G. (1952). A social equilibrium existence theorem. *Proceedings of the National Academy of Sciences, 38*(10), 886-893. https://doi.org/10.1073/pnas.38.10.886
- Ellickson, P. B., & Misra, S. (2008). Supermarket pricing strategies. *Marketing Science, 27*(5), 811-828. https://doi.org/10.1287/mksc.1080.0398
- Fan, K. (1952). Fixed-point and minimax theorems in locally convex topological linear spaces. *Proceedings of the National Academy of Sciences, 38*(2), 121-126. https://doi.org/10.1073/pnas.38.2.121
- Fudenberg, D., & Maskin, E. (1986). The folk theorem in repeated games with discounting or with incomplete information. *Econometrica, 54*(3), 533-554. https://doi.org/10.2307/1911307
- Fudenberg, D., & Tirole, J. (1991). *Game theory*. MIT Press.
- Gibbons, R. (1992). *Game theory for applied economists*. Princeton University Press.
- Glicksberg, I. L. (1952). A further generalization of the Kakutani fixed point theorem, with application to Nash equilibrium points. *Proceedings of the American Mathematical Society, 3*(1), 170-174. https://doi.org/10.1090/s0002-9939-1952-0046638-5
- Green, E. J., & Porter, R. H. (1984). Noncooperative collusion under imperfect price information. *Econometrica, 52*(1), 87-100. https://doi.org/10.2307/1911462
- Harsanyi, J. C. (1967). Games with incomplete information played by "Bayesian" players, I-III: Part I. The basic model. *Management Science, 14*(3), 159-182. https://doi.org/10.1287/mnsc.14.3.159
- Harsanyi, J. C. (1973). Games with randomly disturbed payoffs: A new rationale for mixed-strategy equilibrium points. *International Journal of Game Theory, 2*(1), 1-23. https://doi.org/10.1007/bf01737554
- Hart, S., & Mas-Colell, A. (2000). A simple adaptive procedure leading to correlated equilibrium. *Econometrica, 68*(5), 1127-1150. https://doi.org/10.1111/1468-0262.00153
- Katz, M. L., & Shapiro, C. (1985). Network externalities, competition, and compatibility. *American Economic Review, 75*(3), 424-440.
- Kreps, D. M., & Wilson, R. (1982). Sequential equilibria. *Econometrica, 50*(4), 863-894. https://doi.org/10.2307/1912767
- Mas-Colell, A., Whinston, M. D., & Green, J. R. (1995). *Microeconomic theory*. Oxford University Press.
- McKelvey, R. D., & Palfrey, T. R. (1995). Quantal response equilibria for normal form games. *Games and Economic Behavior, 10*(1), 6-38. https://doi.org/10.1006/game.1995.1023
- Milgrom, P., & Roberts, J. (1990). Rationalizability, learning, and equilibrium in games with strategic complementarities. *Econometrica, 58*(6), 1255-1277. https://doi.org/10.2307/2938316
- Morris, S., & Shin, H. S. (1998). Unique equilibrium in a model of self-fulfilling currency attacks. *American Economic Review, 88*(3), 587-597.
- Nagel, R. (1995). Unraveling in guessing games: An experimental study. *American Economic Review, 85*(5), 1313-1326.
- Nash, J. F. (1950a). The bargaining problem. *Econometrica, 18*(2), 155-162. https://doi.org/10.2307/1907266
- Nash, J. F. (1950b). Equilibrium points in n-person games. *Proceedings of the National Academy of Sciences, 36*(1), 48-49. https://doi.org/10.1073/pnas.36.1.48
- Osborne, M. J., & Rubinstein, A. (1990). *Bargaining and markets*. Academic Press.
- Porter, R. H. (1983). A study of cartel stability: The Joint Executive Committee, 1880-1886. *The Bell Journal of Economics, 14*(2), 301-314. https://doi.org/10.2307/3003634
- Rubinstein, A. (1982). Perfect equilibrium in a bargaining model. *Econometrica, 50*(1), 97-109. https://doi.org/10.2307/1912531
- Selten, R. (1975). Reexamination of the perfectness concept for equilibrium points in extensive games. *International Journal of Game Theory, 4*(1), 25-55. https://doi.org/10.1007/bf01766400
- Shaked, A., & Sutton, J. (1984). Involuntary unemployment as a perfect equilibrium in a bargaining model. *Econometrica, 52*(6), 1351-1364. https://doi.org/10.2307/1913509
- Topkis, D. M. (1978). Minimizing a submodular function on a lattice. *Operations Research, 26*(2), 305-321. https://doi.org/10.1287/opre.26.2.305
- Vives, X. (1990). Nash equilibrium with strategic complementarities. *Journal of Mathematical Economics, 19*(3), 305-321. https://doi.org/10.1016/0304-4068(90)90005-t
:::
