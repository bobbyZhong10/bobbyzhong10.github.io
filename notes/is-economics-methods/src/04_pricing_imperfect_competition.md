---
title: "Pricing, Monopoly & Imperfect Competition"
subtitle: "Market Power, Price Discrimination, and Oligopoly"
seriesline: "Foundations of Information Systems Economics · Chapter 4"
series: "Foundations of Information Systems Economics: Empirical Methods for Information Systems Economics"
sidebartitle: "Foundations of Information Systems Economics"
sidebarsub: "Theory · Methods · Applications"
chapterlabel: "Chapter 4 · Pricing, Monopoly & Imperfect Competition"
---

## Introduction

The same demand curve and the same marginal cost can give rise to wildly different profits. A firm that owns the market pushes its price up to wherever demand elasticity allows; two firms competing in quantities earn less, but still earn something positive; if what they sell is perfectly homogeneous and they compete directly in prices, profit can slide all the way to zero. The price schedule has not become any more "reasonable"; what has changed is the competitive rules behind it. So when you see a high price, the question to ask is not whether the firm is greedy, but what demand, what rivals, and what constraints it faces.

Market power turns price from a given number into something the firm chooses. Monopoly pricing needs the demand and welfare of Chapter 1, oligopoly needs the games of Chapter 2, and second-degree price discrimination needs the screening of Chapter 3. This chapter brings the first three chapters together in Cirrus's market: Cirrus is the dominant e-commerce platform in a particular niche, and Garnet is its rival.

The chapter proceeds as follows. Section 1 sets out a profit caught in the middle. Section 2 writes down the monopolist's problem and the Lerner inverse-elasticity theorem. Section 3 covers monopoly pricing and price discrimination, including the output test for third-degree price discrimination, the limits of price discrimination, nonlinear pricing and bundling, and durable goods and the Coase conjecture. Section 4 covers imperfect competition, running from Cournot, Bertrand, and product differentiation through collusion, vertical relations, and entry deterrence to pass-through (the response of price to a cost shock), which ties the theory to data, and price leadership. Sections 2 to 4 are written with numbered assumptions and theorems, each theorem followed by a proof sketch and an "In words" paragraph. Section 5 reads three anchoring papers. Section 6 settles the accounts on Cirrus's five data blocks: C0, a hand-computed block; C1, 500 logit duopoly markets; C2, personalized pricing for twenty thousand merchants; C3, bundling; and C4, price leadership. Section 7 discusses failure modes. The demand estimation, markup recovery, and merger simulation of Chapters 28 and 29 all stand on this chapter's first-order conditions: pick the wrong game, and even a precisely solved equilibrium is just the answer for some other market.

## 1 A Profit Caught in the Middle

::: {.case}
Cirrus is the dominant e-commerce platform in a particular niche. To shrink the problem to a size that can be done by hand, we abstract its core business into selling "operating services" to merchants: the more a merchant uses, the more the platform charges. We summarize this market with a single linear demand, $p = 12 - q$, where $p$ is the price per unit of operating activity (think of it as the price equivalent of a commission, in dollars) and $q$ is total merchant usage; the platform's marginal cost of providing the service is $c = 2$.
:::

An analyst runs the textbook monopoly-pricing calculation. Setting marginal revenue equal to marginal cost, $12 - 2q = 2$, gives quantity $q = 5$, price $p = 7$, and profit $(7 - 2)\times 5 = 25$. He writes down a conclusion: Cirrus uses its market power to extract a profit of 25, merchants are left with consumer surplus of 12.5, and market power creates a deadweight loss of 12.5, which he says is the unavoidable price of monopoly.

Every step of that 25 is computed correctly, yet it is neither a ceiling nor a floor. Start with the ceiling. Part of that "unavoidable" 12.5 is a choice of pricing instrument, not an iron law of market power: if Cirrus lays out a tiered commission menu and lets merchants of different sizes self-select, profit rises from the 25 of a single package to 26 (Section 6 works this out in a version with two types of merchants); if it could tailor a price to each merchant, it could pocket the 12.5 of consumer surplus together with the 12.5 that was deadweight loss, for a profit of 50. Now the floor. The 25 rests on the premise that Cirrus owns the market, and the moment its rival Garnet enters, the price is pushed down: if the two compete in quantities (Cournot), the price falls to about 5.33 and total industry profit shrinks to about 22.2; if the two offer an undifferentiated service and compete in prices (Bertrand), the price falls to marginal cost 2 and profit goes to zero.

So the answer to "how much does the platform earn" is not a number but a range: the upper bound is set by the pricing instrument (can it price discriminate, can it use a nonlinear tariff), and the lower bound by the competitive structure (are there rivals, and in what dimension do they compete). The isolated 25 is only the reading for one cell, "single price plus exclusive monopoly." Section 3 marks the upper bound, Section 4 the lower bound, and Section 6 works out each of these numbers on Cirrus's market and adds three things that do not fit in a hand calculation: how the mode of competition is written into the pass-through of cost shocks, who actually gets the surplus when prices are personalized by predicted willingness to pay, and how bundling and the timing of pricing can raise revenue or prices without collusion.

## 2 Market Power and How to Measure It

The price taker of Chapter 1 faces horizontal demand: ask for one cent more and it sells nothing, so it can only price at marginal cost. A firm with market power faces downward-sloping demand: raising its price loses some customers but not all of them, so it can trade off "earning more by raising the price" against "earning less by losing customers," and set its price above marginal cost. The monopolist chooses quantity (or equivalently price) to maximize profit $\pi(q) = p(q)\,q - c(q)$, where $p(q)$ is inverse demand; the first-order condition is that marginal revenue equals marginal cost:

$$\underbrace{p(q) + q\,\frac{dp}{dq}}_{MR} = \underbrace{c'(q)}_{MC}$$

Marginal revenue is price plus the negative term $q\,dp/dq$: to sell one more unit the firm must cut the price, the cut drags down every unit already sold, and this is exactly what pushes monopoly output below the competitive level. First we state the two assumptions this problem needs.

::: {.assumption}
**D.1 (regular demand; MWG 12.B, p. 385, footnotes 5 and 6)** Inverse demand $p(q)$ and cost $c(q)$ are continuous and twice differentiable for $q \geq 0$, with $p'(q) < 0$; $p(0) > c'(0)$, and there is a unique $q^0 > 0$ with $p(q^0) = c'(q^0)$; profit $p(q)\,q - c(q)$ is concave on $[0, q^0]$.

**D.2 (constant marginal cost; MWG Example 12.B.1, p. 386; the Bertrand setting of 12.C, p. 388)** The cost to firm $j$ of producing $q_j$ is $c_j\,q_j$, with constant marginal cost $c_j \geq 0$; in symmetric models all firms have the same marginal cost, written $c$.
:::

The $q^0$ of D.1 is the competitive output, where price equals marginal cost; beyond it price is below cost, so the optimal output need only be sought on the compact set $[0, q^0]$ (MWG footnote 5), and the continuity of profit guarantees that a solution exists. Concavity makes the first-order condition sufficient (footnote 6); under D.2 it amounts to the revenue function $p(q)\,q$ being concave. Strongly convex demand can break it, and then the first-order condition may find a local minimum or a suboptimal point. Log-concave demand is another common sufficient condition; MWG does not use it, and it reappears in Section 4.4 for the existence of differentiated Bertrand equilibrium. D.2 rules out capacity constraints and marginal cost that varies with output: digital services with near-zero reproduction cost have roughly constant marginal cost and satisfy it (the economies of scale that a fixed cost brings do not enter the first-order condition for pricing and matter only for entry), and it fails when capacity binds (Section 4.3). Cirrus's $p = 12 - q$ and $c = 2$ satisfy both: $q^0 = 10$, and profit $(10 - q)\,q$ is concave.

::: {.theorem}
**Theorem 1 (the monopoly first-order condition and the Lerner inverse-elasticity rule; MWG 12.B, equations (12.B.3) and (12.B.4), p. 385; Exercise 12.B.1, p. 428)** Let D.1 hold. Monopoly output $q^m$ is positive and satisfies $p'(q^m)\,q^m + p(q^m) = c'(q^m)$, so $p^m = p(q^m) > c'(q^m)$ and $q^m < q^0$. Writing the price elasticity of demand as $\mathcal{E} = (dq/dp)(p/q)$, at $p^m$

$$\frac{p - c'}{p} = \frac{1}{\lvert\mathcal{E}\rvert}$$

If marginal cost is everywhere positive, demand is elastic at the optimal price: $\lvert\mathcal{E}\rvert > 1$.
:::

The proof sketch has three steps. Step one (D.1): a solution exists on $[0, q^0]$, and the first-order condition reads $p'(q^m)\,q^m + p(q^m) \leq c'(q^m)$, with equality when $q^m > 0$; at $q = 0$ the left side is $p(0) > c'(0)$, so $q^m > 0$ and equality holds; concavity makes it a maximum. Step two ($p' < 0$ in D.1): $p'(q^m)\,q^m < 0$, so marginal revenue is below price and $p^m > c'(q^m)$; $q^m \leq q^0$, and at $q^0$ price exactly equals marginal cost, so $q^m < q^0$. Step three (algebra): dividing both sides of the first-order condition by $p^m$ gives $(p^m - c')/p^m = -p'(q^m)\,q^m/p^m$, and since $dq/dp = 1/p'(q)$, the right side is exactly $1/\lvert\mathcal{E}\rvert$; when $c' > 0$ the left side is less than 1, so $\lvert\mathcal{E}\rvert > 1$.

In words: under D.1, the Lerner index (the $(p - mc)/p$ of Section 4.2 of Chapter 1) is determined by the demand elasticity: the mode of competition is "there is only one firm," and the elasticity is evaluated at the price the firm itself chooses. In this book the markup is the difference $p - mc$, and its share of the price is the Lerner index, which equals $1 - 1/\mathcal{M}$ ($\mathcal{M} = p/mc$ is the markup ratio of Theorem 14 in Section 4.2 of Chapter 1). It is not a causal formula: low elasticity is not the cause of a high Lerner index; both are pinned down together by the demand curve at the optimum. Change the mode of competition and the formula changes: the Lerner index of $n$ Cournot firms is $s_j/\lvert\mathcal{E}\rvert$ ($s_j$ is firm $j$'s output share, Section 4.1), and the single-product logit markup is set by price sensitivity and the share (Section 4.4). Empirical work uses it in reverse: estimate the elasticity, assume a mode of competition, and back marginal cost out of the price, $c' = p\,(1 - 1/\lvert\mathcal{E}\rvert)$, which is exactly what Section 3.1 of Chapter 29 does. If the elasticity is misestimated or the mode of competition misspecified, the recovered marginal cost is wrong too (Sections 3.3 and 3.5 of Chapter 29), and estimating the elasticity takes exogenous price variation (Section 3.3 of Chapter 28). Cirrus's numbers: at $p = 7$ the Lerner index is $(7 - 2)/7 = 0.714$, the elasticity is $\lvert\mathcal{E}\rvert = 7/5 = 1.4$, and its inverse is $1/1.4 = 0.714$; with two Cournot firms each has a share of 1/2, the market elasticity is 0.8, and $0.5/0.8 = 0.625$ is exactly the Lerner index at the Cournot price 5.33.

The social cost of monopoly is stated in the welfare language of Chapter 1: consumer surplus (the area under the demand curve and above the price line; this chapter's $\mathrm{CS}$ is this level, not the loss $\Delta\mathrm{CS}$ of Chapter 1), profit, and the gap between their sum and the competitive level, that is, the deadweight loss (DWL). DWL comes from output being held below $q^0$: trades that should have happened, where the buyer's valuation exceeds marginal cost, do not happen. The welfare effects of price discrimination, competition, and collusion all come down in the end to whether this area shrinks or grows and how total surplus is split between buyers and sellers.

## 3 Core Structure: Monopoly Pricing and Price Discrimination

This section first uses Cirrus's numbers to show what is to be proved, then states numbered assumptions and theorems, followed by proof sketches in which each step names the assumptions it uses; regularity conditions are only mentioned and cited to their source. After each theorem, an "In words" paragraph answers five questions: under which assumptions the result characterizes which quantity; what it becomes, and in which direction it errs, when the assumptions fail; whose choice or which primitive determines the conclusion; what the data can test and what variation that takes; and what Cirrus's numbers say. Assumptions carry the chapter's prefixes D. (demand and cost) and O. (oligopoly); theorems are numbered consecutively within the chapter, with MWG's section, proposition, and page numbers in parentheses, and results that are not in the book cite the original papers. Section 3.1 is the single-price benchmark; Sections 3.2 and 3.3 take up price discrimination, organized by "how much the firm knows about buyers, and whether it can stop resale"; Sections 3.4 and 3.5 cover nonlinear pricing, bundling, and versioning; Section 3.6 turns to intertemporal pricing.

### 3.1 Single-Price Monopoly and Deadweight Loss

Let demand be linear, $p = a - bq$ ($b$ is the slope of demand, unrelated to the bids of Chapter 3), with marginal cost $c$ (MWG Example 12.B.1, p. 386). Then profit is a concave function of $q$, and the first-order condition of Theorem 1 gives $q^m = (a - c)/(2b)$ and $p^m = (a + c)/2$, with competitive output $q^0 = (a - c)/b$. Cirrus's $a = 12$, $b = 1$, and $c = 2$ give $q^m = 5$, $p^m = 7$, profit 25, and competitive output 10.

![Single-price monopoly: demand $p = 12 - q$, marginal revenue $MR = 12 - 2q$, marginal cost $MC = 2$. $MR = MC$ pins down $q_m = 5$, $p_m = 7$. Between demand and the price line is consumer surplus (12.5); the rectangle between price and MC up to $q_m$ is monopoly profit (25); the triangle from $q_m$ to the competitive quantity $q_c = 10$ (that is, $q^0$) is deadweight loss (12.5).](assets/fig/fig_04_monopoly.svg)

Monopoly produces 5 units fewer than competition; on each of those units the buyer's valuation exceeds marginal cost 2, so trades that should have happened did not, and the lost area $\frac{1}{2}\times(10 - 5)\times(7 - 2) = 12.5$ is the deadweight loss, a quarter of the total surplus of 50 under competition. Consumer surplus is $\frac{1}{2}\times(12 - 7)\times 5 = 12.5$. These are exactly the numbers the analyst reported; they assume the monopolist can charge only a single price, and a single price is exactly the self-imposed restriction a monopolist is most reluctant to accept.

### 3.2 Three Kinds of Price Discrimination

Price discrimination means charging different prices to different buyers, or for different units bought by the same buyer, with the goal of capturing the surplus that a single price leaves in buyers' hands. It needs two conditions: the firm must be able, directly or indirectly, to tell buyers' willingness to pay apart, and it must be able to stop low-price buyers from reselling to high-price buyers. Ranked by how much information the firm has, from most to least, price discrimination is traditionally divided into three degrees. We begin with the end where the firm has the most information:

::: {.theorem}
**Theorem 2 (first-degree price discrimination removes the quantity distortion; MWG 12.B, equation (12.B.5), p. 387)** Let buyer $i$'s utility be quasilinear, $u_i(q_i) + m_i$, with $u_i(0) = 0$ and $m_i$ the numeraire; the monopolist knows each $u_i$ and makes each buyer a take-it-or-leave-it offer $(q_i, T_i)$, and buyers cannot resell. Then the optimal offers satisfy $T_i = u_i(q_i)$, the chosen quantities maximize total surplus $\sum_i u_i(q_i) - c\big(\sum_i q_i\big)$, total output equals the socially optimal level, and every buyer's surplus is zero.
:::

The proof sketch has two steps. Step one (quasilinearity, take-it-or-leave-it offers): buyer $i$ accepts if and only if $u_i(q_i) - T_i \geq 0$, so the firm raises $T_i$ up to $u_i(q_i)$ and the buyer's surplus is zero. Step two (complete information, no resale): substituting $T_i = u_i(q_i)$, the firm's profit $\sum_i u_i(q_i) - c\big(\sum_i q_i\big)$ is exactly total surplus, and maximizing it means choosing the socially optimal quantities.

In words: under quasilinearity and complete information, first-degree price discrimination separates efficiency from distribution completely: total surplus is maximal and DWL is zero, and all the surplus goes to the firm. Quasilinearity is essential; it makes money transfers correspond one for one to surplus (Theorem 10 (iii) in Section 4.1 of Chapter 1). When information is incomplete or resale cannot be blocked, the firm falls back on second- or third-degree discrimination, and an efficiency loss usually follows; third-degree discrimination is not necessarily inefficient, since the consumer-optimal segmentation of Theorem 4 has no deadweight loss; how much information the firm sees, and who decides what it sees, is exactly the question of the limits of price discrimination in Section 3.3. Data do not reveal $u_i$; when willingness to pay is observable, what can be done is to compute this upper bound. Cirrus's numbers: with homogeneous merchants whose demand is $p = 12 - q$, first-degree discrimination (or the two-part tariff of Section 3.4) lifts profit from 25 to 50 and output from 5 to 10; on C2's evaluation half, the profit of first-degree discrimination is the efficient surplus, USD 32.35 per merchant per month.

Third-degree price discrimination is the most common kind: the firm cannot price person by person, but it can split buyers into several segments by observable characteristics (by region, student status, or business versus individual account) and charge each segment a single price. Theorem 1 applies within each segment, so the less elastic the segment, the higher its price. If Cirrus could split merchants into a sticky segment with demand $p = 12 - q_1$ and a price-sensitive segment with $p = 8 - q_2$, the two prices would be 7 and 5, the elasticities 1.4 and 1.67, and the Lerner indices 0.714 and 0.6.

Second-degree price discrimination addresses the case in which the firm cannot even see which segment a buyer belongs to and knows only the proportion of each type; it lays out a menu of quantities and total prices and lets buyers self-select, so that their choices reveal their types indirectly. This is the pricing version of the screening of Section 3.3 of Chapter 3: the menu must satisfy incentive compatibility and individual rationality, the high type at the top gets the efficient quantity, and the low type's quantity is pushed below the efficient level, to reduce the information rent that must be conceded to keep the high type from posing as the low type. The tariff $T(q)$ and the payment rule $t(\theta)$ of Chapter 3 are linked by the taxation principle $t(\theta) = T\big(q(\theta)\big)$, so the optimal nonlinear tariff with continuous types comes straight from the screening theorem, Theorem 3 in Section 3.3 of Chapter 3, and this chapter does not reprove it. Section 6 works it out for two types of merchants: the high type gets the efficient quantity 10, the low type is pushed down to 2 (its efficient quantity is 6), the high type's information rent is 8, the menu earns 26, a single package earns at most 25, and complete information would yield 34.

::: {.intuition}
The distinction between the three kinds of price discrimination falls into place once you grasp "how much the firm knows about buyers, and how it prevents resale." First-degree is omniscient: it prices person by person and unit by unit, taking all the surplus while producing the most. Third-degree sees a label and prices by type, charging the less elastic type more. Second-degree cannot even tell the type, so it lays out a menu and lets buyers self-select, trading distortion at the bottom for rent at the top. The less information, the less complete the extraction and the less visible the distortion. Platforms do both third-degree (pricing by region and account type) and second-degree (tiered subscriptions and stepped commissions); once you see which kind of information about buyers a platform currently holds, you know which degree it is using.
:::

### 3.3 The Welfare of Third-Degree Discrimination: The Output Test and the Limits of Price Discrimination

Third-degree price discrimination does not lower profit, because the firm can always revert to a uniform price; its effect on total welfare, however, has no universal sign. Put the two segments above together: at the uniform price 6 both segments are served, with outputs 6 and 2; at the discriminatory prices 7 and 5 the outputs are 5 and 3; total output is 8 in both cases, yet welfare falls from 52 to 51. If the price-sensitive segment is replaced by $p = 4 - q_2$, the uniform price 7 serves only the sticky segment, while under discrimination the sensitive segment is served at price 3; total output rises from 5 to 6, and welfare rises from 37.5 to 39. For welfare to rise, total output must rise first.

::: {.theorem}
**Theorem 3 (the output test for third-degree price discrimination; Schmalensee 1981; Varian 1985)** Let D.2 hold, let market $i$'s demand $X_i(p)$ slope downward, let the uniform price be $p^\ast$, and let the price in market $i$ under discrimination be $p_i^\ast$. (i) If third-degree price discrimination does not increase total output relative to the uniform price, total surplus does not rise, and it falls strictly as soon as output changes in some market; equivalently, for discrimination to raise welfare, total output must rise. (ii) If demand in every market is linear and all markets are served at the uniform price, total output under discrimination equals total output under the uniform price, so welfare falls strictly as long as the market prices are not all equal.
:::

The proof sketch has three steps. Step one (D.2, downward-sloping demand): when market $i$ switches from the uniform price to its discriminatory price, the units involved in the output change $\Delta x_i$ have valuations between $p^\ast$ and $p_i^\ast$, so the welfare change satisfies $(p^\ast - c)\,\Delta x_i \geq \Delta W_i \geq (p_i^\ast - c)\,\Delta x_i$, with strict inequalities when $\Delta x_i \neq 0$. Step two (adding up): $\Delta W \leq (p^\ast - c)\sum_i \Delta x_i$, so $\Delta W \leq 0$ when total output does not increase, which is (i). Step three (linear demand, all markets served): when $X_i(p) = A_i - B_i p$, the monopoly output in each market is $A_i/2 - B_i c/2$, and these add up to exactly the uniform-price output under aggregate demand $\sum_i A_i - \big(\sum_i B_i\big) p$; when the prices are not all equal, output changes in some market, so by (i) welfare falls strictly, which is (ii).

In words: under D.2, the welfare effect of third-degree price discrimination is bracketed by the change in output; when demand is linear and all markets are served, discrimination only moves the same total output from people who value it more to people who value it less, so welfare must fall. Under general demand only a necessary condition remains: for welfare to rise, total output must rise, and a rise in output does not guarantee a rise in welfare; Aguirre, Cowan and Vickers (2010) give conditions that use the demand curvature of each market to sign the changes in output and welfare, and this book only points to their existence. The most common rise in output comes from opening a new market, so the conclusion is determined by the curvature of each market's demand and by which markets the uniform price serves. The implication for platforms is direct: when rolling out personalized pricing, first check whether the number of people served has increased; the test needs only the quantities under the two pricing schemes, while judging the magnitude takes demand estimates market by market. Cirrus's numbers: the two pairs of segments above are exactly case (ii) and the case of "opening a new market"; in C2, fine segmentation lifts the share of merchants served from 0.550 to 0.686.

One step further: if any segmentation at all is allowed, how far can third-degree price discrimination push welfare? The next result answers this question completely.

::: {.theorem}
**Theorem 4 (the limits of price discrimination; Bergemann, Brooks and Morris 2015)** Consider a single good with unit demand and zero marginal cost, with valuations taking finitely many values $v_1 < v_2 < \dots < v_{N_{\mathrm{v}}}$ and a known distribution over the aggregate market; a market segmentation splits the aggregate market into segments, and the seller sets a profit-maximizing price in each segment (when several prices are equally optimal, the segmentation also specifies which one is taken). Write $\pi^u$ for the monopoly profit at the uniform price and $\mathbb{E}[v]$ for the efficient total surplus. Then a pair of consumer and producer surplus $(\mathrm{CS}, \mathrm{PS})$ can be achieved by some segmentation if and only if $\mathrm{CS} \geq 0$, $\mathrm{PS} \geq \pi^u$, and $\mathrm{CS} + \mathrm{PS} \leq \mathbb{E}[v]$.
:::

The proof sketch has three steps; for the full argument and the exact premises, see the original paper. Step one (necessity; unit demand, zero marginal cost): in every segment the seller can still charge the uniform price, so $\mathrm{PS} \geq \pi^u$; buyers purchase voluntarily, so $\mathrm{CS} \geq 0$; and with zero marginal cost, total surplus cannot exceed $\mathbb{E}[v]$. Step two (the consumer-optimal vertex; finitely many values, ties broken by the segmentation): on a support $S$, construct an "extremal market" with $\Pr(v \geq v_k) = \min(S)/v_k$, so that the seller is indifferent among all prices in $S$; split the aggregate market step by step into such segments, each of whose supports contains the uniform price, so that when each segment is priced at its lowest valuation everyone buys and total surplus is efficient, while the seller gets exactly $\pi^u$. Under the same segmentation, if the seller prices at the highest valuation instead, the result is the vertex $\mathrm{CS} = 0$, $\mathrm{PS} = \pi^u$; the segmentation that fully reveals valuations gives the first-degree vertex $\mathrm{PS} = \mathbb{E}[v]$. Step three (mixing; no further assumption): dividing the aggregate market in proportion among the three segmentations yields any convex combination of the three vertices, so the whole triangle is attainable.

In words: under these conditions, the welfare outcomes that market segmentation can achieve are exactly this triangle, and where in it the outcome falls is decided not by demand but by who designs the information: a seller choosing the segmentation for profit picks full revelation of valuations and ends up at the first-degree vertex; a consumer or regulator designing the segmentation can reach the other vertex, where trade is efficient and the whole gain goes to consumers. The welfare question of price discrimination thus becomes the information-design problem of Section 3.7 of Chapter 3. The theorem says only what is attainable, not what is common: personalized pricing by segments of predicted willingness to pay is a segmentation chosen by the seller, and it is only third-degree discrimination; second-degree discrimination by menus is outside what this triangle characterizes. The triangle itself cannot be tested; what data can do, when willingness to pay is observable, is locate where a given scheme falls. Cirrus's numbers: on C2's evaluation half the optimal uniform price is USD 23, with profit USD 13.25 per merchant per month and efficient surplus USD 32.35; the consumer-optimal segmentation, found by construction, uses 180 segments to push consumer surplus to USD 19.10 while leaving profit unchanged; coarse and fine segmentation both fall inside the triangle; the uniform price, by contrast, is set on the training half, at USD 24, and its profit of USD 13.21 lies USD 0.04 below the floor of the triangle.

![The welfare triangle on C2's evaluation half (ten thousand merchants): the horizontal axis is consumer surplus, the vertical axis Cirrus's profit, both in USD per merchant per month. The top vertex is first-degree discrimination (profit equals the efficient surplus 32.35); the right vertex is the consumer-optimal segmentation (consumer surplus 19.1, profit 13.25). Uniform price (K = 1): profit 13.21, consumer surplus 12.54; coarse segmentation (K = 4): 15.23 and 11.78; fine segmentation (K = 50): 15.61 and 11.63. Coarse and fine segmentation fall inside the triangle; the uniform price of USD 24 was set on the training half, and the evaluation half's best single price is USD 23, so the uniform-price point lies USD 0.04 below the floor. Gray dots are segmentations with K = 1 to 200.](assets/fig/fig_04_personalized.svg)

### 3.4 Nonlinear Pricing: Two-Part Tariffs and Bundling

In platform and SaaS pricing, the menu of Section 3.2 takes the form of two families of tools, both of them nonlinear pricing (total payment is not proportional to the quantity bought).

The first family is the two-part tariff (Oi 1971): $T(q) = F + p\,q$, a usage-independent access fee $F$ plus a per-unit price $p$. With homogeneous buyers the optimum is clean: set the per-unit price at marginal cost, let buyers buy the efficient quantity, then use $F$ to sweep up the entire consumer surplus in one go; the two instruments reproduce Theorem 2. With Cirrus's homogeneous merchant demand $p = 12 - q$, a per-unit price of 2 and an access fee of 50 give profit 50, twice the single-price 25, and output rises to the efficient 10. With heterogeneous buyers, a high $F$ extracts more from large buyers but shuts small buyers out, so the optimal two-part tariff trades off the markup per unit against coverage; with continuous types the optimum is a fully nonlinear tariff (Theorem 3 in Section 3.3 of Chapter 3), and the two-part tariff is only its affine approximation. Membership fees plus usage fees, listing fees plus transaction commissions, subscription floors plus overage charges: all are two-part tariffs.

The second family is bundling. Adams and Yellen (1976) show that even when each good is already priced optimally on its own, pure bundling can earn more when buyers' valuations of the two goods are negatively correlated, though not always; mixed bundling (selling the bundle and the individual goods at the same time) contains separate selling and pure bundling as special cases and does at least as well as either. The reason is that bundling smooths out differences in willingness to pay across buyers. Cirrus has two value-added services A and B with near-zero marginal cost; one type of merchant values A at 10 and B at 2, and the other type the reverse: sold separately, each service's optimal price is 10, sold only to the half that values it, for a combined revenue of 20; bundled into one package priced at 12, both types buy, for a revenue of 24. As the number of goods grows, the power of bundling comes from a different source: the law of large numbers.

::: {.theorem}
**Theorem 5 (revenue and efficiency of large bundles; Bakos and Brynjolfsson 1999; the lower bound written directly from Chebyshev's inequality)** Let there be $J$ goods with zero marginal cost, each buyer wanting at most one unit of each, and let each buyer's valuation of good $j$, $v_{ij} \geq 0$, be independent and identically distributed, with mean $\mathbb{E}[v_{ij}] > 0$ and finite variance $\operatorname{Var}(v_{ij})$. For any $\gamma \in (0, 1)$, bundling the $J$ goods into one package priced at $J\,\mathbb{E}[v_{ij}](1 - \gamma)$ yields an expected revenue per good of at least

$$\mathbb{E}[v_{ij}]\,(1 - \gamma)\left[1 - \frac{\operatorname{Var}(v_{ij})}{J\,\gamma^2\,\mathbb{E}[v_{ij}]^2}\right]$$

Hence, as $J \to \infty$, the expected revenue per good of the optimal pure bundle tends to $\mathbb{E}[v_{ij}]$, and DWL as a share of the efficient surplus tends to zero.
:::

The proof sketch has three steps. Step one (independent and identically distributed, finite variance): a buyer buys the whole bundle if and only if the average valuation $J^{-1}\sum_j v_{ij}$ is no lower than $\mathbb{E}[v_{ij}](1 - \gamma)$; by Chebyshev's inequality, the probability of not buying is at most $\operatorname{Var}(v_{ij})/\big(J\gamma^2\mathbb{E}[v_{ij}]^2\big)$. Step two (unit demand): revenue per good equals the price divided by $J$, times the probability of purchase, and substituting gives the lower bound; as $J \to \infty$ it tends to $\mathbb{E}[v_{ij}](1 - \gamma)$; $\gamma$ can be arbitrarily small, while revenue per good can never exceed the buyers' average valuation $\mathbb{E}[v_{ij}]$, so the revenue per good of the optimal pure bundle tends to $\mathbb{E}[v_{ij}]$. Step three (zero marginal cost): on each good, revenue, consumer surplus, and DWL add up to $\mathbb{E}[v_{ij}]$; revenue tends to it, so the other two tend to zero together.

In words: under these conditions, a large bundle makes each buyer's average valuation of the package predictable, so the seller can sell to almost everyone at a price slightly below the mean; revenue and efficiency rise together, and the remaining consumer surplus and DWL tend to zero together. The speed of convergence is set by the ratio of the standard deviation to the mean of the valuation distribution (the coefficient of variation): the larger it is, the more goods are needed. When valuations are positively correlated (a common taste factor), the average valuation no longer converges to a constant and the advantage of bundling shrinks; when marginal cost is positive, goods that go unused still cost something, so the conclusion of Bakos and Brynjolfsson (1999) applies mainly to information goods. Testing it needs exogenous changes in bundle size and buyer-level purchase data. Cirrus's numbers: in C3, fifty thousand consumers value 100 information goods independently, each valuation drawn from the uniform distribution on $[0, 1]$, with mean 1/2; sold separately, revenue per good is 0.249 and DWL is 0.254 of the efficient surplus; with 20 goods bundled the figures are 0.375 and 0.060, and with 100 goods 0.432 and 0.026, while the share buying the whole bundle rises from 0.497 to 0.970.

### 3.5 Versioning: Designing Quality Menus for Information Goods

In this subsection $q$ is quality, not quantity. The high fixed cost, near-zero reproduction cost, and designable quality of information goods make versioning a central form of second-degree price discrimination (Shapiro and Varian 1999): the firm first builds a high-quality master, then carves out versions by features, speed, timeliness, API limits, or service level, letting buyers with different willingness to pay self-select; the limitations of the low-end version need not come from technical cost, and may be a degradation added deliberately to maintain incentive compatibility. Let types be $\theta_H > \theta_L$, utility $\theta q - T$, and the unit cost of a product of quality $q$ be $c(q)$, strictly convex ($c' > 0$, $c'' > 0$), so that efficient quality satisfies $\theta = c'(q)$. When the menu $(q_L, T_L)$, $(q_H, T_H)$ satisfies both types' incentive compatibility and individual rationality, single crossing gives the same structure as in Section 3.2: no distortion at the top, and $q_L$ distorted downward to reduce the high type's information rent. The convex cost is essential: with zero cost of quality, profit is linear in $q_L$, and the optimal menu collapses to a single version or excludes the low type (Salant 1989); for a zero-cost information good to be versioned, valuations cannot take the form $\theta v(q)$, scaled up in proportion to type. Competition compresses the room for artificial degradation, because a rival can use a better low-end version to take the distorted users, so as competition intensifies we should expect the quality gap between versions to narrow. Observing only that premium users pay more does not establish that price discrimination exists; quality differences, selection into products, and genuine heterogeneity in willingness to pay must first be separated.

### 3.6 Durable Goods and Intertemporal Pricing: Coase, Stokey, and Hotelling's Rule

A monopolist selling a durable good is competing with its own future self: it sells to high-valuation buyers at a high price today; tomorrow, facing the remaining low-valuation buyers, it wants to cut the price; and rational high-valuation buyers foresee this and would rather wait. Coase (1972) conjectured that when the seller can reprice frequently and cannot commit not to cut prices, this competition with itself drives the price toward the competitive level.

::: {.theorem}
**Theorem 6 (the Coase conjecture; Gul, Sonnenschein and Wilson 1986, Theorem 1 and its corollary, pp. 161-162; Theorem 3 and footnote 14, p. 168)** Let marginal cost be zero (with a constant cost $c$, read valuations and prices as net of $c$), and let a continuum of buyers of unit measure each buy at most one unit of the durable good, with lowest valuation $v_{\min} \geq 0$; seller and buyers share the discount factor $\delta < 1$, the seller quotes one price each period and cannot commit to future prices, and each buyer who has not yet bought decides each period whether to buy or wait. (i) If $v_{\min} > 0$ (gap) and valuation, as a function of the buyer's rank, satisfies a Lipschitz condition at the lowest-valuation end, then for generic (almost all) valuation distributions the equilibrium path is unique, buyers' strategies are stationary, prices fall along the path, and the market is sold out in finitely many periods. (ii) For any $\gamma > 0$ there is $\delta_0 < 1$ such that, when $\delta > \delta_0$: if $v_{\min} = 0$ (no-gap), the seller's first price is below $\gamma$ in every equilibrium in which buyers' strategies are stationary; if $v_{\min} > 0$, the first price is below $v_{\min} + \gamma$ in every equilibrium.
:::

For the proof, see the appendix of Gul, Sonnenschein and Wilson (1986). This book only states the result, because the proof has to handle stationary strategies for a continuum of buyers, randomization off the path, and uniqueness of equilibrium, and would run far longer than this chapter allows. The intuition is the one the paper gives on p. 169: when buyers' strategies are stationary, the seller can at any time move tomorrow's price forward to today; the cost is that today's buyers pay a little less, and the gain is that all later sales arrive one period earlier. Equilibrium requires that such a move not pay, and as the period length goes to zero ($\delta \to 1$) that can hold only if the price itself is already close to the lowest valuation.

In words: with stationary buyer strategies, a durable-goods monopolist without commitment loses its pricing power over high-valuation buyers: the more frequently it reprices, the closer its opening price comes to the lowest valuation, and in the no-gap case profit tends to zero. The conclusion is determined by the lowest valuation relative to cost (gap or no-gap), the frequency of repricing, and equilibrium selection (stationarity). Section 6 of the paper itself points to two exits: a continuing inflow of new buyers (p. 173), and rising costs that let the seller commit to limiting output (pp. 172-173; with quadratic production costs profit is positive); other commonly cited ones are costs of changing prices, reputation, leasing instead of selling, and most-favored-customer clauses; the subscriptions and time-limited pricing of digital products are all ways of manufacturing commitment. Testing it means comparing the price paths of sellers with different commitment power (selling versus leasing). Cirrus's numbers: C2's Cirrus Pro (Section 6.3) is a monthly subscription, and each month is a new sale, so C2 can be treated as a static problem. If it were switched to a one-time purchase whose price could be changed daily, it would be the gap case, because the lowest willingness to pay in C2, USD 1, is above the zero marginal cost, and the opening price would be pushed down to around USD 1.

The value of commitment is clearest from the other side. Should a durable-goods monopolist that can commit to its prices skim (set a high price for high-valuation buyers first, then cut it step by step)? Stokey's (1979) answer is no: facing forward-looking rational buyers, when buyers and seller discount equally and cost does not fall over time, the optimum for a seller that can commit is a constant price, because an anticipated price cut only induces high-valuation buyers to wait. Skimming pays only under other conditions: buyers are myopic; the good is perishable and each period brings a fresh batch of buyers who cannot wait (airline seats, concert tickets); or cost or quality falls over time. The opposite, penetration pricing (a low price first, a higher price later), is often optimal under network effects, learning-by-doing, switching costs, or experience goods; a platform's freemium (a free basic version with paid add-ons) and low-price customer acquisition aim to push past the adoption threshold of Theorem 15 in Section 3.11 of Chapter 2, and the pricing of two-sided platforms is in Section 2.1 of Chapter 32.

The other pole of dynamic pricing is a fixed, exhaustible stock, and the classic result is Hotelling's Rule (Hotelling 1931). The owner can at any moment choose one of two options: sell a unit now, take the scarcity rent $p_0 - c$, and invest it at the interest rate; or keep it and sell it in period $\ell$ for $p_\ell - c$ (in this subsection $\ell$ is the period index). No arbitrage requires the two paths to be equivalent, so the scarcity rent must grow at the interest rate. With the discount factor $\delta$ (1 divided by 1 plus the interest rate), this reads

$$p_\ell - c = (p_0 - c)\,\delta^{-\ell}$$

If the rent grows more slowly than the interest rate, it is better to sell everything now and invest; if it grows faster, it is better to hold on and not sell. Here the owner is treated as a price taker; with pricing power, what grows at the interest rate is marginal revenue minus marginal cost. A limited stock of Cirrus's premium ad slots is such a stock: with placement cost 2, interest rate 5%, and initial rent 4 (initial price 6), the rent runs 4, 4.2, 4.41, the price 6, 6.2, 6.41, and the present value is always 4 (Section 6.1). Spectrum licenses, reserved server capacity, and API call quotas should all be released over time in the same way: dumping them all at once gives up the appreciation, and hoarding them unused gives up the interest.

## 4 Imperfect Competition: The Several Faces of Oligopoly

Monopoly is one firm dominating, perfect competition is countless small firms, and platform markets often fall in between: a few firms with market power competing with one another, which is oligopoly. The outcome depends heavily on the dimension in which firms compete and how many times the competition is repeated, and each of these forms is an application of the game theory of Chapter 2.

### 4.1 Cournot: Competing in Quantities

In Cournot competition firms choose quantities simultaneously, and total output determines the price through demand. Cirrus and Garnet face $p = 12 - Q$ ($Q = q_1 + q_2$), both with marginal cost 2: the first-order condition of Cirrus's profit $(12 - q_1 - q_2 - 2)\,q_1$ with respect to $q_1$ gives the best response $q_1 = (10 - q_2)/2$; the more the rival produces, the less a firm should produce itself. The two best-response lines meet at $q_1 = q_2 = 10/3$, with price $16/3 \approx 5.33$, profit per firm $100/9 \approx 11.11$, and $200/9 \approx 22.2$ in total.

![Cournot best responses: Cirrus's $q_1 = (10 - q_2)/2$ and Garnet's $q_2 = (10 - q_1)/2$ meet at the symmetric equilibrium $q_1 = q_2 = 10/3$, with price 5.33. For comparison: under collusion each produces 2.5 (half the joint monopoly output 5), at price 7; under competition each produces 5, at price 2.](assets/fig/fig_04_cournot.svg)

::: {.assumption}
**O.1 (conditions for a continuous game, that is, G.2 of Section 3.3 of Chapter 2; MWG Proposition 8.D.3, p. 253)** Each firm's strategy set is a nonempty, convex, compact subset of some Euclidean space, and its payoff is continuous in all firms' strategies and quasiconcave in its own strategy.
:::

The Cournot quantity set $[0, \infty)$ is not compact; once output is restricted to $[0, q^0]$ (beyond $q^0$ price is below cost, MWG footnote 5), the first two conditions hold. Quasiconcavity is not automatic: MWG Figure 12.C.3 gives an example in which profit is not quasiconcave, the best response jumps, and no pure-strategy equilibrium exists; there the two firms' costs differ. With linear demand, profit is concave in own output, and O.1 holds.

::: {.theorem}
**Theorem 7 (Cournot equilibrium: existence, uniqueness with linear demand, and convergence to competition as the number of firms grows; MWG 12.C, p. 393, citing Proposition 8.D.3; Example 12.C.1, pp. 391-392; equation (12.C.6) and Exercise 12.C.7, p. 393)** Let $n$ firms choose quantities $q_j \in [0, q^0]$ simultaneously, let inverse demand $p(Q)$ satisfy the conditions of D.1 on $p$, and let all marginal costs equal $c$ (D.2). (i) If each firm's profit is quasiconcave in its own output, a pure-strategy Nash equilibrium exists. (ii) If $p(Q) = a - bQ$ with $a > c \geq 0$ and $b > 0$, the equilibrium is unique and symmetric: $q_j = (a - c)/\big(b(n + 1)\big)$ and $p = (a + nc)/(n + 1)$. (iii) The symmetric equilibrium satisfies $p'(Q^\ast)\,Q^\ast/n + p(Q^\ast) = c$, and the price is above $c$; as $n \to \infty$, $p(Q^\ast) \to c$.
:::

The proof sketch has four steps. Step one (O.1): the strategy set $[0, q^0]$ is nonempty, convex, and compact, profit is continuous, and with quasiconcavity the existence theorem for continuous games, Theorem 3 in Section 3.3 of Chapter 2 (MWG Proposition 8.D.3), gives (i). Step two (D.2, linear demand): firm $j$'s first-order condition is $a - c - bQ - b\,q_j \leq 0$, with equality when $q_j > 0$; if some firm produced zero, then $a - c - bQ \leq 0$, so the other firms could only produce zero as well, giving $Q = 0$, which contradicts $a > c$; so every firm produces, all in the same amount, and adding up gives a unique $Q$, which is (ii); MWG Example 12.C.1 covers only two firms, so the uniqueness argument for $n$ firms is supplied here. Step three (D.1, D.2): under symmetry, combining the first-order condition $p'(Q)\,q_j + p(Q) = c$ with $q_j = Q/n$ gives (12.C.6), and $p' < 0$ puts the price above $c$. Step four (D.1: $p'$ is continuous on the compact set $[0, q^0]$, hence bounded): $Q^\ast < q^0$, so $p'(Q^\ast)\,Q^\ast/n \to 0$ and the price tends to $c$.

In words: under O.1 and D.2, Cournot is a mild form of competition: each firm holds back because producing more depresses the price on its own sales, the markup falls gradually with the number of firms, and each firm's Lerner index is its output share divided by the market elasticity, in contrast to the jump of Bertrand, where two firms already push the price down to cost. When quasiconcavity fails, an equilibrium may not exist; when costs differ, shares and markups are ordered by cost. The conclusion is determined by the strategic variable (quantity, or the capacity of Section 4.3), the number of firms, and demand curvature. In the data, the relation between price and the number of firms is not causal, because entry is itself endogenous: the Bresnahan-Reiss approach of Section 4.1 of Chapter 30 uses the relation between market size and the number of firms to back out the intensity of competition. Cirrus's numbers: with two firms the price is 5.33 and industry profit 22.2, 0.889 of monopoly profit; with 8 firms the price is 3.11 and DWL 0.62; with 100 firms the price is 2.099 and DWL is only 0.0049.

### 4.2 Bertrand: Competing in Prices and That Paradox

In Bertrand competition firms choose prices simultaneously, and buyers all turn to whoever quotes the lowest price. As long as the rival prices above cost, I can undercut slightly and take the whole market; the rival thinks the same, and the two keep cutting until neither dares go any lower.

::: {.theorem}
**Theorem 8 (the Bertrand paradox; MWG Proposition 12.C.1, p. 388)** Let two firms sell a homogeneous product, both with marginal cost $c > 0$ (D.2); let market demand $x(p)$ be continuous, strictly decreasing where $x(p) > 0$, and zero above some finite price, with $x(c)$ positive and finite; the two firms quote prices simultaneously, buyers all turn to the lower-priced firm, and at equal prices they split equally. Then the unique Nash equilibrium is $p_1 = p_2 = c$, with zero profit for both; it remains the unique equilibrium in mixed strategies (MWG footnote 9 and Exercise 12.C.2).
:::

The proof sketch has three steps. Step one (D.2, buyers all turn to the lower-priced firm): $(c, c)$ is an equilibrium, since raising the price sells nothing and cutting it loses money. Step two (D.2): when the lowest price is below $c$, the firm quoting it makes a loss, and raising it to $c$ does better. Step three (continuous demand): when the lowest price is above $c$, if the two prices are equal, either firm can cut slightly, take the whole market, and nearly double its profit; if they are unequal, the higher-priced firm earns zero and can take the whole market at a positive profit by quoting slightly below its rival; when the lowest price equals $c$ and the other firm's is above $c$, the firm at $c$ can raise its price slightly and earn a positive profit. So there is no other equilibrium.

In words: under the four conditions of a homogeneous product, equal costs, no capacity constraint, and one-shot play, two firms reproduce perfect competition. Each condition is an exit from the paradox: capacity (Section 4.3 pulls the outcome back to Cournot), differentiation (Section 4.4), repetition (Section 4.5), and unequal costs (the low-cost firm takes the whole market at its rival's cost). The conclusion is determined by whether the product is homogeneous and by price being the strategic variable; in the same duopoly, competing in prices gives 2 and competing in quantities 5.33, and where reality falls depends on whether firms behave more as if they set capacity first or as if they could change prices at any moment. A low markup observed in a homogeneous market does not by itself say which exit is at work. Cirrus's numbers: in C1's homogeneous-product variant the two firms' costs each vary between USD 1.5 and USD 2.5, the Bertrand price equals the higher of the two costs, averaging 2.16, with an average markup of 0.33; under the same costs the average Cournot price is 5.35 and the average markup 3.36.

### 4.3 Kreps-Scheinkman: Capacity First, Prices After

The first exit from the Bertrand paradox is capacity. Firms often build capacity $q_j^{\mathrm{cap}}$ first and set prices afterward; when capacity is limited, the logic that quoting the lowest price takes the whole market fails, and the rival can sell to residual demand and hold its price above cost. Residual demand depends on the rationing rule, that is, on which buyers get the low price when the low-price firm lacks the capacity to serve everyone. Kreps and Scheinkman (1983) prove that when firms choose capacities first and then compete in prices, the unique subgame perfect equilibrium of the two-stage game is the Cournot outcome, under several conditions. The one MWG names in the text of 12.C, p. 395, is efficient rationing, that is, the demand of high-valuation buyers is served first by the low-price firm; the other conditions include a concave demand curve and requirements on the shape of capacity costs. For the details, see the original paper. MWG states the result in words only, without proof, and this book does the same. It gives Cournot a more realistic reading: quantity competition stands for the long-run choice of capacity, and price competition takes place in the short run, given capacity. If the rationing rule is switched to random rationing or queuing, the conclusion may fail.

### 4.4 Product Differentiation: Hotelling and Logit

The second exit from the Bertrand paradox is differentiation: when products differ in buyers' eyes (brand, interface, ecosystem, location), some buyers are always willing to pay a bit more for what they prefer, and the logic of cutting a cent to grab all of demand fails.

::: {.theorem}
**Theorem 9 (price equilibrium in Hotelling's linear city; MWG Example 12.C.2, pp. 396-398; d'Aspremont, Gabszewicz and Thisse 1979)** Let a mass $M$ of buyers be distributed uniformly on a segment of length 1, with the two firms fixed at its two ends, both with marginal cost $c$ (D.2); each buyer buys at most one unit, with gross value $v$, and at distance $d$ from the firm she buys from she pays an additional transport cost $t\,d$ ($t > 0$); $v$ is large enough relative to $c$ and $t$ (MWG takes $v > c + 3t$) that everyone buys. Then the unique Nash equilibrium is $p_1 = p_2 = c + t$, each firm sells $M/2$, and each earns profit $tM/2$. If the transport cost is changed to $t\,d^2$ while the two firms stay fixed at the ends, the conclusion is unchanged (the same indifferent location; see step one).
:::

The proof sketch has three steps. Step one (everyone buys): the indifferent buyer is located at $z = (t + p_2 - p_1)/(2t)$ (MWG equation (12.C.7)); with quadratic cost, $p_1 + t z^2 = p_2 + t(1 - z)^2$ simplifies to the same $z$, so with the ends fixed the two cost forms give the same demand $Mz$. Step two (D.2): firm 1 maximizes $(p_1 - c)(t + p_2 - p_1)M/(2t)$, and the first-order condition gives the best response $p_1 = (t + p_2 + c)/2$ (the middle part of MWG equation (12.C.12)). Step three (symmetry): $p = (t + p + c)/2$, that is, $p = c + t$; MWG further shows that there is no asymmetric equilibrium.

In words: with the ends fixed and everyone buying, the markup exactly equals the degree of differentiation $t$: the more buyers care about mismatch, the less firms fear a rival's price cut; as $t \to 0$ we are back to Bertrand. When firms choose their own locations under linear transport cost, a price equilibrium fails to exist once the two firms locate close together, and d'Aspremont, Gabszewicz and Thisse (1979) restore existence by switching to quadratic transport cost; the degree of differentiation itself becomes a strategic choice. The conclusion is determined by how dispersed buyers' tastes are and by where the firms are located. In the data, $t$ corresponds to the substitution pattern of demand and has to be read from demand estimates (the elasticities and diversion ratios of Section 4.5 of Chapter 27; the diversion ratio is the share of a product's lost sales that goes to another product). Cirrus's numbers: with $c = 2$, $t = 1$, and a buyer mass of 1, the price 3 lies between the 2 of undifferentiated Bertrand and the 5.33 of Cournot, and profit per firm is 0.5.

A line segment has room for only two firms. Empirical work more often uses logit demand: each product has a quality $\xi_j$, buyers share the same price sensitivity $\alpha$ (in this chapter $\alpha$ is price sensitivity, unrelated to the bargaining power of the negotiation in Chapter 2), and there is also an option not to buy. In C1's baseline market Cirrus's quality is 4 and Garnet's 3, $\alpha = 1$, and both marginal costs are 2; the Bertrand-Nash prices are 3.746 and 3.318, with shares 0.427 and 0.241.

::: {.theorem}
**Theorem 10 (the single-product markup under logit demand; this chapter's own derivation; for existence, see Caplin and Nalebuff 1991)** Let each firm sell a single product, let buyer $i$'s utility from buying product $j$ be $u_{ij} = \xi_j - \alpha p_j + \varepsilon_{ij}$ and from not buying $u_{i0} = \varepsilon_{i0}$, with the $\varepsilon_{ij}$ independent and identically Gumbel distributed (type I extreme value), so that the share is $s_j = \exp(\xi_j - \alpha p_j)\big/\big(1 + \sum_k \exp(\xi_k - \alpha p_k)\big)$; marginal cost $c_j$ is constant (D.2), and firms set prices simultaneously. Then every interior Nash equilibrium satisfies $p_j - c_j = 1/\big(\alpha(1 - s_j)\big)$; the own-price elasticity is $\mathcal{E}_{jj} = -\alpha p_j(1 - s_j)$, and the Lerner index equals $1/\lvert\mathcal{E}_{jj}\rvert$.
:::

The proof sketch has two steps; the existence and uniqueness of equilibrium need additional conditions, which this book only mentions: Caplin and Nalebuff (1991) obtain existence from concavity-type conditions on the distribution of preferences (a generalization of log-concavity), and for the uniqueness conditions, see the original paper. Step one (logit): $\partial s_j/\partial p_j = -\alpha s_j(1 - s_j)$. Step two (D.2, single product): firm $j$ maximizes $(p_j - c_j)\,s_j$, with first-order condition $s_j + (p_j - c_j)\,\partial s_j/\partial p_j = 0$, and substituting gives the markup; the elasticity is $(\partial s_j/\partial p_j)(p_j/s_j) = -\alpha p_j(1 - s_j)$, and the markup divided by the price is exactly $1/\lvert\mathcal{E}_{jj}\rvert$.

In words: with a single product and logit demand, the markup is determined only by price sensitivity and the firm's own share: the smaller $\alpha$ and the larger the share, the higher the markup; this is the Lerner formula of Theorem 1 in the form it takes under logit elasticities. IIA (the independence of irrelevant alternatives, Section 3.3 of Chapter 27) makes substitution depend only on shares, so who the rival is does not affect the markup; with random-coefficients demand (Chapter 28) the markup depends on diversion. When one firm owns all the products, every product's markup becomes $1/(\alpha s_0)$ (this chapter's own derivation); multiproduct firms and mergers are left to Section 4.2 of Chapter 29. The conclusion is determined by $\alpha$, the shares, and ownership; in estimation, $\alpha$ is first estimated with instruments for price (Section 3.3 of Chapter 28), and marginal cost is then backed out of the first-order condition (Section 3.1 of Chapter 29). Cirrus's numbers: in the baseline market Cirrus's markup is $1.746 \approx 1/(1 - 0.427)$ and Garnet's $1.318 = 1/(1 - 0.241)$, and Cirrus's Lerner index 0.466 is exactly the inverse of its own-price elasticity 2.145. When the two firms price jointly, both prices are 4.163, with markup $2.163 \approx 1/0.462$; an analyst who reads this pair of joint prices as Bertrand recovers marginal costs of 2.515 and 2.994, while the true values are both 2.

### 4.5 Collusion: The Theory Is in Chapter 2, the Identification Is Here

The third exit from the Bertrand paradox is repetition. Theorem 16 in Section 4.1 of Chapter 2 showed that with infinite repetition and enough patience, cooperation can be sustained by the credible threat that "deviation triggers retaliation." Carried over to pricing: two homogeneous Bertrand firms that sustain the monopoly price with Nash reversion strategies (once anyone cuts the price, both price at cost forever after) are in a subgame perfect equilibrium if and only if the discount factor is at least 1/2 (MWG Proposition 12.D.1, p. 401). Tacit collusion is therefore first of all a product of repeated interaction, not a direct result of market structure in a one-shot game; whether it can be sustained still depends on structural factors such as the number of firms: with $n$ firms the threshold rises to $(n - 1)/n$ (MWG Exercise 12.D.2, p. 405). This poses a hard problem for empirical work: behind a high price there may be collusion, or differentiation, or Cournot under capacity constraints, and the price alone cannot tell them apart. Telling them apart takes additional structure: Porter (1983), read in Section 5.2, identified the behavioral pattern from the cartel's switches between collusion and price wars; Sections 3.3 and 3.4 of Chapter 29 cover what variation the identification of conduct (the mode of competitive behavior) requires; and Section 2.1 of Chapter 34 covers how learning algorithms learn to collude without any agreement.

### 4.6 Vertical Relations: Double Marginalization, Vertical Integration, and RPM

Platforms are often embedded in a vertical chain: buying and reselling (a reseller, corresponding to the wholesale model), or only matching and taking a commission while the seller sets the price (a marketplace, corresponding to the agency model). The most important phenomenon is double marginalization (Spengler 1950): an upstream monopolist sells at wholesale price $w$ (unrelated to the budget of Chapter 1) to a downstream monopolist retailer, which adds a markup of its own, so the retail price is above the integrated monopoly price, while combined upstream and downstream profit is actually lower. Cirrus's numbers: integrated, the retail price is 7 and profit 25; separated, upstream sets the wholesale price at 7, the retail price rises to 9.5, upstream earns 12.5 and downstream 6.25, 18.75 in total, and consumer surplus falls from 12.5 to 3.125 (Section 6.1). The root of the problem is that the downstream firm, when it marks up, ignores the drag on upstream sales. Vertical integration, a two-part wholesale contract ($w = c$ plus a fixed franchise fee, the vertical version of the two-part tariff of Section 3.4), or resale price maintenance (RPM: upstream sets the retail price directly at 7) can each pull the outcome back to the integrated one. The agency model takes a fixed share of sales revenue: a seller with zero marginal cost then maximizes only a fixed multiple of revenue, which, when the platform's own marginal cost is also zero, has the same solution as integration, and the second markup disappears (this chapter's own derivation); when the seller has a positive marginal cost, the commission amounts to scaling its marginal cost up, dividing it by "1 minus the commission rate," and the distortion does not disappear. The dispute between the wholesale and agency models in e-books and app stores is precisely over who removes this markup; integration may at the same time bring foreclosure, and the merging party's price cut alone cannot separate efficiency from exclusion.

### 4.7 Entry Deterrence: Dixit's Capacity Commitment and a Taxonomy

The high profit of monopoly attracts entrants, so the incumbent has an incentive to deter entry. A verbal threat, "if you come in I will start a price war," is not credible (the empty threat of Chapter 2); only a sunk and visible investment can turn the threat into a commitment.

::: {.theorem}
**Theorem 11 (capacity commitment and entry deterrence; MWG Appendix B, pp. 423-427, text and Figures 12.BB.1 to 12.BB.8, no numbered proposition; Dixit 1980)** A three-stage game: the incumbent first pays $c_{\mathrm{cap}}$ per unit to build capacity $q_I^{\mathrm{cap}}$; the potential entrant observes it and decides whether to pay the entry cost $F_{\mathrm{E}}$ (not the access fee $F$ of the two-part tariff) and enter; after entry the two choose quantities simultaneously, the entrant's cost per unit is $c + c_{\mathrm{cap}}$, and the incumbent cannot produce beyond its capacity and pays only $c$ per unit within it. Let both firms' profits be strictly concave in own output (concave inverse demand suffices), with downward-sloping best responses. Then (i) the incumbent's post-entry best response is the smaller of its best response at marginal cost $c$ and its capacity; it can credibly induce any outcome on the entrant's best-response line up to a certain point, and capacity beyond that point would not be used after entry. (ii) The subgame perfect equilibrium is one of three cases: entry is blockaded (the incumbent builds the monopoly capacity it would choose without any threat, and entry is still unprofitable); entry is accommodated (entry is unavoidable, or deterrence earns less than accommodation; the incumbent picks the Stackelberg point (the outcome when the incumbent sets its quantity first), where its own isoprofit curve is tangent to the entrant's best-response line, or an end point of the credible range); or entry is deterred (capacity is built up to the point where the entrant just cannot profit, provided this earns more than accommodation).
:::

The proof sketch has three steps, by backward induction. Stage three (strict concavity, downward-sloping responses): capacity is sunk, so within capacity the incumbent's marginal cost is only $c$, and its best response shifts out, but not beyond capacity; when capacity exceeds the output at "the intersection of the best response at marginal cost $c$ with the entrant's best-response line," the excess would not be produced after entry, so it is not a credible threat, which is (i). Stage two (the entrant observes capacity): the entrant anticipates the post-entry outcome and enters only if its profit is no lower than $F_{\mathrm{E}}$. Stage one (capacity sunk and visible; the credible range of (i)): the incumbent chooses capacity within the credible range: if the monopoly capacity without any threat already makes entry unprofitable, entry is blockaded; otherwise it compares the profit from "building capacity up to the point where the entrant just cannot profit" with that of the best accommodation point and takes the higher (MWG p. 427, Figure 12.BB.8), which is (ii).

In words: when profit is strictly concave and capacity is sunk and visible, a credible deterrent can come only from capacity that would actually be used after entry: the sunk cost lowers the incumbent's ex post marginal cost and thereby changes its ex post best response. If capacity can be shed or cannot be seen, or post-entry competition is in prices, the conclusion changes; the table below gives the corresponding taxonomy. The conclusion is determined by the entry cost, the capacity cost, and demand; even when deterrence succeeds and entry does not happen, the threat of entry raises output and welfare (MWG p. 427). Excess capacity observed in the data is not by itself proof of deterrence; it takes the counterfactual "what would happen if the rival entered," and Chapter 30 covers the estimation of entry games. Cirrus's numbers: if Garnet enters and the two play symmetric Cournot, each earns 11.11; when Garnet's entry cost is below this number and Cirrus makes no commitment, entry is profitable.

Dixit's quantity competition is only one cell. Fudenberg and Tirole (1984) cross the effect of investment (whether it makes the incumbent tougher or softer after entry) with the type of post-entry competition: to deter entry, the incumbent must always look tough; only when it accommodates entry does it matter whether competition is in strategic substitutes or strategic complements.

| When accommodating entry | Investment makes the incumbent tough | Investment makes the incumbent soft |
|---|---|---|
| Strategic complements (e.g., prices) | puppy dog: underinvest | fat cat: overinvest |
| Strategic substitutes (e.g., quantities) | top dog: overinvest | lean and hungry: underinvest |

When entry is accommodated, showing weakness is sometimes what is optimal: in price competition, looking tough provokes the rival to cut its price. The other three tools work through different mechanisms; what they share is that an action today changes the entrant's assessment of what entry will pay: limit pricing uses a low price, under asymmetric information about costs, to send the signal "my cost is low" (Milgrom and Roberts 1982), an application of the signaling games of Section 3.2 of Chapter 3; predatory pricing relies on cross-market reputation or a financing advantage; and product proliferation fills up the product space, leaving no gap for an entrant.

### 4.8 Pass-through, Conduct, and Empirical Predictions

The most useful empirical restrictions of price theory are often not directional judgments like "the price will rise," but the size of the derivative of price with respect to a cost shock. Under linear inverse demand $p = a - bQ$ the single-product monopoly price is $p = (a + c)/2$, so

$$\frac{\partial p}{\partial c}=\frac{1}{2}.$$

Switching to $n$ symmetric Cournot firms, it becomes $n/(n + 1)$. Each number is just one combination of demand shape and mode of competition.

::: {.theorem}
**Theorem 12 (three closed forms for pass-through; derived from MWG equation (12.B.4), p. 385; Example 12.B.1, p. 386; Example 12.C.1 and Exercise 12.C.7, pp. 391 and 393; (iii) follows from Theorem 10)** Write pass-through as $\mathrm{PT} = dp/dc$. (i) For a single-product monopoly, with D.1 and D.2 holding and the second-order condition $2p'(Q) + Q\,p''(Q) < 0$: $\mathrm{PT} = 1/\big(2 + Q\,p''(Q)/p'(Q)\big)$; linear demand gives $1/2$, and constant-elasticity demand ($\lvert\mathcal{E}\rvert > 1$) gives $\lvert\mathcal{E}\rvert/(\lvert\mathcal{E}\rvert - 1)$. (ii) For $n$ symmetric Cournot firms with linear demand, when all firms' costs move together: $\mathrm{PT} = n/(n + 1)$. (iii) For a logit monopoly with a single product (the setting of Theorem 10): $\mathrm{PT} = 1 - s$.
:::

The proof sketch has three steps. Step one (D.1, D.2, the second-order condition): differentiating the first-order condition $p(Q) + Q\,p'(Q) = c$ with respect to $c$ gives $(2p' + Q\,p'')\,dQ/dc = 1$, and $dp/dc = p'\,dQ/dc$, which is (i); linear demand has $p'' = 0$, constant-elasticity demand has $Q\,p''/p' = -(1 + 1/\lvert\mathcal{E}\rvert)$, and substituting gives the two special cases. Step two (D.2, linear demand): differentiating the $p = (a + nc)/(n + 1)$ of Theorem 7 (ii) with respect to $c$ gives (ii). Step three (logit): totally differentiating $p - c = 1/\big(\alpha(1 - s)\big)$ and using $ds = -\alpha s(1 - s)\,dp$ gives $dp\,\big(1 + s/(1 - s)\big) = dc$, which is (iii).

In words: pass-through is determined jointly by demand curvature and the mode of competition. The curvature term $Q\,p''/p'$ alone can push monopoly pass-through from 0.5 to 1 and above: with a constant elasticity of 2 monopoly pass-through is 2, and with exponential demand it is 1; the number of firms then moves it in different directions (figure below). So own pass-through alone cannot identify conduct: different demand shapes can give the same number (a monopoly with linear demand and a logit monopoly with a share of exactly one half both give 0.5), and the mode of competition barely changes it (the true values in C1 below: 0.592 under Bertrand and 0.613 under joint pricing). Weyl and Fabinger (2013) extend these results into general principles of tax incidence under imperfect competition (they write pass-through as $\rho$; this book writes PT), and this book states only the three closed forms above; Miklós-Thal and Shaffer (2021) point out that some of Weyl and Fabinger's formulas hold only under specific conditions, and their comment should be read before citing Weyl and Fabinger's conclusions on price discrimination and exogenous entry. To distinguish conduct, either pin down the shape of demand or use variables that rotate the demand curve (Bresnahan 1982), as Sections 3.3 and 3.4 of Chapter 29 develop. With several products, the difference lies in the cross terms. Write $\mathrm{PT}_{jk} = \partial p_k/\partial c_j$ (the pass-through of cost $c_j$ to price $p_k$); under logit, joint pricing makes every product's markup $1/(\alpha s_0)$, and total differentiation gives $\mathrm{PT}_{jj} = 1 - s_j$ and $\mathrm{PT}_{jk} = -s_j$ ($j \neq k$): cross pass-through is negative, whereas under Bertrand it is positive. In the data, all it takes is an exogenous, observable cost shifter and an OLS regression of price on it. Cirrus's numbers (block C1, 500 markets): under Bertrand, the OLS coefficient of Cirrus's price on its own cost is 0.609 (standard error 0.017, true value 0.592), and on Garnet's cost 0.144 (true value 0.137); under joint pricing, own is 0.638 (true value 0.613) and cross $-0.112$ (true value $-0.151$). Own pass-through barely moves, while cross pass-through changes sign; all eight coefficients, four for each firm's price (Garnet's four are not listed), are within 4 conventional standard errors of their true values. With Cirrus alone, its logit monopoly price is 4, its share 1/2, and its pass-through 0.5, indistinguishable from a monopoly with linear demand: this is a coincidence of demand shape and has nothing to do with conduct.

![Pass-through of a common cost shock as a function of the number $n$ of symmetric firms. Left: under Cournot, linear demand gives $n/(n + 1)$, rising from 0.5 to 0.909; exponential demand stays at 1; constant-elasticity demand with elasticity 2 gives $2n/(2n - 1)$, falling from 2 to 1.053. Right: under logit (quality 3.5, $\alpha = 1$, cost 2), Bertrand pass-through rises from 0.566 to 0.993, and with all products owned by one firm it falls from 0.566 to 0.325, exactly the no-purchase share $s_0$.](assets/fig/fig_04_passthrough.svg)

These comparative statics should be written as the joint response of several outcome variables, rather than fixating on a single price coefficient:

| Theory or institutional change | Comparative static | Joint empirical prediction | Main alternative explanation |
|---|---|---|---|
| Marginal cost rises | Price rises by the pass-through rate | Own price, sales, and rival price move together | Demand shock, menu cost |
| Products more differentiated | Bertrand competition softens | Markup rises, cross diversion falls | Quality changes marginal cost at the same time |
| Market becomes segmentable | Segment prices diverge by elasticity | Price, coverage, and composition all change | Selection, personalized quality |
| Vertical integration | Double markup internalized | Merging party's retail price falls, sales rise | Foreclosure, cost synergies |
| Repeated interaction more stable | Collusion more sustainable | Fewer price wars, more persistent punishment after shocks | Common costs, algorithmic synchronization |

### 4.9 Higher Prices without Collusion: Price Leadership and Pricing Algorithms

Block C4 returns to C1's baseline market and changes only one thing: who sets its price first. With simultaneous pricing the two prices are 3.746 and 3.318; when Cirrus sets its price first and Garnet sets its own after seeing it, both prices rise.

::: {.theorem}
**Theorem 13 (price leadership raises both firms' prices; this chapter's own derivation)** Let two single-product firms sell substitutes (firm 1's profit rises with its rival's price), and let firm 2's best response to firm 1's price, $p_2 = \mathrm{BR}_2(p_1)$, be unique, differentiable, and strictly increasing (prices are strategic complements); firm 1 sets its price first, and firm 2 observes it and best-responds; the simultaneous-pricing Bertrand-Nash equilibrium $(p_1^\ast, p_2^\ast)$ is interior, $\mathrm{BR}_2'(p_1^\ast) > 0$, and $\partial\pi_1/\partial p_2 > 0$ at that point; firm 1's profit along its rival's response, $\pi_1\big(p_1, \mathrm{BR}_2(p_1)\big)$, is differentiable and strictly quasiconcave in $p_1$ and has a maximum. Then under sequential pricing both firms' prices are above the simultaneous-pricing Bertrand-Nash prices $(p_1^\ast, p_2^\ast)$.
:::

The proof sketch has three steps. Step one (interior Bertrand-Nash equilibrium): at $p_1^\ast$, $\partial\pi_1/\partial p_1 = 0$, and $p_2^\ast = \mathrm{BR}_2(p_1^\ast)$. Step two (substitutes, strategic complements): the total derivative of the leader's profit at $p_1^\ast$ is $\partial\pi_1/\partial p_1 + (\partial\pi_1/\partial p_2)\,\mathrm{BR}_2'(p_1^\ast)$; the first term is zero and the second is the product of two positive numbers, so the total derivative is positive. Step three (strict quasiconcavity): a positive total derivative means the maximum lies to the right of $p_1^\ast$, so the leader's price is above $p_1^\ast$; since $\mathrm{BR}_2$ is strictly increasing, the follower's price is also above $p_2^\ast$. Under logit, $\mathrm{BR}_2' = s_1 s_2/(1 - s_2) > 0$, so strategic complementarity holds; the quasiconcavity of the leader's profit has to be checked case by case.

In words: under strategic complements, the first mover's commitment softens competition: the leader knows its rival will follow a price rise, so it raises its price first, and both prices end up above their simultaneous-pricing levels, without any agreement, punishment, or repeated game. Under quantity competition the first mover instead produces more (Stackelberg, that is, the accommodation point of Theorem 11), and in undifferentiated Bertrand moving first brings no advantage. The conclusion is determined by the timing of moves, the ability to commit, and the slope of the best response. Brown and MacKay (2023) carry this logic over to pricing algorithms: when firms reprice at different frequencies and algorithms respond to rivals automatically, they find, in Markov perfect equilibrium, that simple pricing algorithms can raise the price level (as the abstract puts it; for the exact propositions, see the original paper). The division of labor with the Q-learning collusion of Section 2.1 of Chapter 34 is that this section covers commitment and timing, while Chapter 34 covers learning; in the data, telling the two apart means looking at the timing of price changes and the response to shocks. Cirrus's numbers: when Cirrus sets its price first the two prices are 3.807 and 3.326, and when Garnet goes first they are 3.765 and 3.425, all below the joint-pricing 4.163; when Cirrus leads, Garnet's profit rises 2.58% and Cirrus's only 0.10%; that the follower gains more is only this example's number, not a proposition of Brown and MacKay.

## 5 Anchoring Papers

Three studies take this chapter's theory to real markets: price discrimination in Broadway theaters, collusion and price wars in a railroad cartel, and personalized pricing on a digital platform; each is laid out by paper, method, data, results, and limitations.

### 5.1 Leslie (2004)

::: {.case}
Paper: "Price Discrimination in Broadway Theater," *RAND Journal of Economics*. It measures the second- and third-degree price discrimination of Section 3.2 against real pricing data and is a model example of structural empirical work on price discrimination.

Method and data: Broadway theaters charge different prices for different seats and different purchase channels (full-price box office, discount coupons, same-day discount booth) for the same performance, which involves both third-degree discrimination by observable label and second-degree discrimination by self-selection. The author estimates a demand model with consumer heterogeneity and then simulates how the theater's profit and consumer welfare change under a uniform price or other alternatives.

Results: relative to a uniform price, the observed discrimination raises the theater's profit and has little effect on aggregate consumer welfare, and most of the gain in profit is not a transfer from consumers (see the paper for magnitudes). The welfare effect of price discrimination has to be estimated, not presumed.

Limitations: the conclusions depend on the functional form of demand and the specification of heterogeneity; Theorem 3 says welfare cannot rise unless output rises, and the change in output in the counterfactuals is determined precisely by these specifications.
:::

### 5.2 Porter (1983)

::: {.case}
Paper: "A Study of Cartel Stability: The Joint Executive Committee, 1880-1886," *Bell Journal of Economics*. It carries the game-theoretic characterization of collusion and price wars in Section 4.5 over to a real railroad-freight cartel, and is the starting point of empirical work on collusion and of conduct tests.

Method: the theory of Green and Porter (1984) (Theorem 18 in Section 4.2 of Chapter 2) holds that when cartel members cannot observe whether a rival secretly expands output and see only a market price mixed with demand noise, collusion must be sustained by "entering a stretch of punitive price war together whenever the price is abnormally low"; the price war is one link of the collusive equilibrium, not a breakdown of collusion. Porter uses a switching regression (a regression that switches between the two states of collusion and price war) to estimate pricing behavior in the two states.

Data and results: freight rates and volumes of a 19th-century American railroad freight cartel. Rates differ significantly between the two states and are markedly lower in price-war periods, supporting the counterintuitive conclusion that price wars are part of collusion, and grounding the "punishment phase" of the repeated games of Chapter 2 in data.

Limitations: identification relies on a specific model of how the cartel operates; distinguishing price wars under collusion from competition itself is still delicate work, and the identification problem of Section 4.5 is resolved only provisionally, by structural assumptions.
:::

### 5.3 Dubé and Misra (2023)

::: {.case}
Paper: "Personalized Pricing and Consumer Welfare," *Journal of Political Economy*. The question is the same as C2's: what personalized pricing implemented with machine learning (an extreme form of third-degree price discrimination) means for profit and consumer welfare.

Method and data: a randomized pricing field experiment on the paywall of the recruiting platform ZipRecruiter, in which 7,867 prospective customers were randomly assigned to 10 monthly price cells, with the original price of USD 99 as the control (p. 150); a LASSO logit (variable selection with an L1 penalty), run under a weighted likelihood bootstrap, estimates demand as a function of many customer characteristics, personalized prices are constructed from it, and a second experiment validates them out of sample (Section 5.3 of Chapter 7 uses the same experiment to discuss demand estimation).

Results: the optimal uniform price raises profit by 55%; by the model's predictions, personalization raises expected profit by a further 19% on top of that, and by 86% relative to the original price (pp. 131 and 135). Personalization lowers linearly aggregated consumer surplus by 25% relative to the optimal uniform price, and conversion and total surplus fall as well; but 63% of customers get a price below the uniform price (USD 327), and under an inequality-averse welfare function personalization may be better (p. 167). Conversion and total surplus falling together fits Theorem 3; C2's segmentation expands the number served, and total surplus usually rises, but rising output is only a necessary condition for rising welfare.

Limitations: a single firm and a single product; how good the personalized prices are depends on the regularization and functional form of the demand model; the experimental period is short, and longer-run retention is only a robustness check.
:::

Together, the three papers mark where this chapter's theory meets the data: the welfare of price discrimination has to be computed from demand estimates; collusion can be identified only with additional structure; and whether personalized pricing is good or bad for consumers starts with whether it expands the number of transactions.

## 6 Running Case: Pricing and Competition in the Cirrus Market

The code below runs under R 4.6.1. Only the key statements are shown, with several statements merged onto one line; every number in the text comes from the output of an actual run, and the result lines in the code that begin with `#` copy what the run printed. `cournot_n`, `linear_pd`, `bertrand_logit`, `joint_max_logit`, `passthrough_logit`, `price_leader`, `bundle_opt`, and `bbm_segment` implement the formulas of Sections 2 to 4 directly, `gen_C1` to `gen_C3` generate data under the settings described below, and `pd_run` carries out C2's fitting, segmentation, pricing, and evaluation.

### 6.1 C0: The Hand-Computed Cirrus Market

C0 has no random numbers, and every number in it can be computed by hand; it is used to check the formulas of Sections 2 to 4. The first excerpt covers monopoly, the output test for third-degree discrimination, the menu for two types of merchants, and the path of Cournot toward competition:

```r
a <- 12; b <- 1; cost <- 2                                    # demand p = 12 - q, marginal cost 2 (USD)
qm <- (a - cost) / (2 * b); pm <- a - b * qm; qc <- (a - cost) / b
round(c(pm = pm, profit_m = (pm - cost) * qm, cs_m = 0.5 * (a - pm) * qm, dwl_m = 0.5 * (qc - qm) * (pm - cost),
        lerner = (pm - cost) / pm, inv_elas = qm / pm), 4)
unlist(linear_pd(a = c(12, 8), b = 1, cost = 2)[c("p_uniform", "Q_uniform", "W_uniform", "Q_pd", "W_pd")])
unlist(linear_pd(a = c(12, 4), b = 1, cost = 2)[c("p_uniform", "Q_uniform", "W_uniform", "Q_pd", "W_pd")])
thH <- 10; thL <- 6; mu <- 0.5                                # two merchant types; platform marginal cost 0
qL <- thL - (mu / (1 - mu)) * (thH - thL); rentH <- (thH - thL) * qL
TL <- thL * qL - qL^2 / 2; TH <- (thH * thH - thH^2 / 2) - rentH
c(qL = qL, rentH = rentH, TL = TL, TH = TH, menu = mu * TH + (1 - mu) * TL,
  single = mu * (thH * thH - thH^2 / 2), first_best = mu * thH^2 / 2 + (1 - mu) * thL^2 / 2)
round(sapply(c(1, 2, 8, 100), function(n) { o <- cournot_n(a, b, cost, n)
  c(n = n, p = o$p, industry_profit = n * o$profit, dwl = 0.5 * (qc - n * o$q) * (o$p - cost)) }), 4)
#       pm profit_m     cs_m    dwl_m   lerner inv_elas
#   7.0000  25.0000  12.5000  12.5000   0.7143   0.7143
# p_uniform Q_uniform W_uniform      Q_pd      W_pd
#         6         8        52         8        51
# p_uniform Q_uniform W_uniform      Q_pd      W_pd
#       7.0       5.0      37.5       6.0      39.0
#         qL      rentH         TL         TH       menu     single first_best
#          2          8         10         42         26         25         34
#                 [,1]    [,2]   [,3]     [,4]
# n                1.0  2.0000 8.0000 100.0000
# p                7.0  5.3333 3.1111   2.0990
# industry_profit 25.0 22.2222 9.8765   0.9803
# dwl             12.5  5.5556 0.6173   0.0049
```

In the menu, the two types of merchants each make up half the market, the high type $\theta_H = 10$ and the low type $\theta_L = 6$, with gross value $\theta q - q^2/2$. The low type gets 2 and pays $T_L = 10$; the high type gets 10 and pays $T_H = 50 - 8 = 42$; expected profit is $\frac{1}{2}\times 42 + \frac{1}{2}\times 10 = 26$. The high type gets utility 8 by choosing its own tier and also 8 by posing as the low type, so incentive compatibility binds exactly; the low type gets 0, so individual rationality binds. Selling a single package, the best Cirrus can do is to offer only the high type a package with $q = 10$ at a price of 50, for an expected 25 (serving both types yields only 18, and the best per-unit price is 4, with profit 16); that it equals the monopoly profit of 25 in Section 3.1 is only a coincidence. The difference of 8 between the 34 under complete information and the menu's 26 is the cost of information: part flows to the high type as information rent, and part evaporates through the distortion at the bottom. The second excerpt covers Hotelling, the two-part tariff, bundling, Hotelling's Rule, and double marginalization:

```r
vA <- c(10, 2); vB <- c(2, 10)                                  # the two merchant types' values for A and B
revA <- max(2 * sum(vA >= 2), 10 * sum(vA >= 10)); revB <- max(2 * sum(vB >= 2), 10 * sum(vB >= 10))
c(hotelling_p = cost + 1, two_part_F = 0.5 * (a - cost)^2, separate = revA + revB, bundle = 12 * sum((vA + vB) >= 12))
rent <- 4 * 1.05^(0:2)                                          # interest rate 0.05, initial scarcity rent 4
rbind(rent = rent, price = 2 + rent, pv_rent = rent / 1.05^(0:2))
w <- (a + cost) / 2; qd <- (a - w) / 2; pd <- (a + w) / 2       # upstream maximises (w - 2)(12 - w) / 2
c(w = w, retail = pd, upstream = (w - cost) * qd, downstream = (pd - w) * qd, cs = 0.5 * (a - pd) * qd)
# hotelling_p  two_part_F    separate      bundle
#           3          50          20          24
#         [,1] [,2] [,3]
# rent       4  4.2 4.41
# price      6  6.2 6.41
# pv_rent    4  4.0 4.00
#          w     retail   upstream downstream         cs
#      7.000      9.500     12.500      6.250      3.125
```

### 6.2 C1: 500 Logit Duopoly Markets

In each local market Cirrus and Garnet sell differentiated operating services with qualities 4 and 3 and $\alpha = 1$; marginal cost is 2 plus a cost shifter drawn independently for the two firms from the uniform distribution on $[-0.5, 0.5]$; the shifter is the per-unit payment-processing fee, which the analyst can observe. A demand shock $\mathcal{N}(0, 0.25^2)$ is added to quality, independent of costs. In each market the equilibrium prices are solved under two kinds of conduct: firm $j$ maximizes $\pi_j + \kappa\sum_{k \neq j}\pi_k$, where $\kappa = 0$ is Bertrand-Nash and $\kappa = 1$ is joint profit maximization. The analyst runs OLS of Cirrus's price on the two firms' costs and reports conventional standard errors:

```r
P1 <- list(M = 500, delta = c(Cirrus = 4, Garnet = 3), alpha = 1, c_base = c(2, 2), cost_halfwidth = 0.5,
           xi_sd = 0.25, hom_a = 12, hom_b = 1, hom_eta_sd = 1)
set.seed(4101); C1 <- gen_C1(P1)        # cost shifters, demand shocks, Bertrand and joint prices in 500 markets
tr <- readRDS("truth04.rds")$C1           # true mean pass-through, PT[price k, cost j]
out <- NULL
for (cond in c("bertrand", "joint")) {
  fit <- lm(C1[[paste0("p1_", cond)]] ~ c1 + c2, data = C1)              # Cirrus's price on both costs
  cf <- summary(fit)$coefficients
  out <- rbind(out, data.frame(conduct = cond, cost = c("Cirrus", "Garnet"), estimate = cf[c("c1", "c2"), 1],
                               se = cf[c("c1", "c2"), 2], truth = tr[[paste0("pt_", cond, "_mean")]][1, ]))
}
print(out, digits = 3, row.names = FALSE)
pm1 <- bertrand_logit(4, 1, 2)            # Cirrus alone: a single-product logit monopoly
c(p = pm1, s = logit_share(pm1, 4, 1), pt = drop(passthrough_logit(4, 1, 2)))
#   conduct   cost estimate     se  truth
#  bertrand Cirrus    0.609 0.0167  0.592
#  bertrand Garnet    0.144 0.0162  0.137
#     joint Cirrus    0.638 0.0164  0.613
#     joint Garnet   -0.112 0.0159 -0.151
#   p   s  pt
# 4.0 0.5 0.5
```

The true values are the averages over the 500 markets of the equilibrium pass-through rates, computed market by market with the implicit function theorem. In the homogeneous-product variant under the same costs ($p = a_m - Q$, with the intercept $a_m$ fluctuating around 12), the true pass-through of each firm's cost to the Cournot price is 1/3, and OLS gives 0.410 and 0.261 (standard errors about 0.055); the deviation comes from the noise of the intercept shocks.

### 6.3 C2: Personalized Pricing of Cirrus Pro

Cirrus Pro is a monthly data-analytics and promotion tool for merchants; it is software, with zero marginal cost. Twenty thousand merchants each subscribe to at most one copy a month; the log of willingness to pay is 3.2 plus 0.45 times standardized sales, plus 0.15 times standardized store age, plus 0.2 for electronics, minus 0.15 for home goods, plus noise with standard deviation 0.5. The median for the baseline merchant (apparel, average sales and store age) is about USD 24.5, and the overall median about USD 25; the characteristics explain about half of the variance; willingness to pay is recorded in whole dollars. A random half is the training half: Cirrus once ran a survey eliciting willingness to pay, and these merchants' willingness to pay is on record; in reality a platform usually sees only whether a merchant buys and has to estimate demand first. The analyst fits log willingness to pay by OLS on the training half, splits merchants into $K$ segments by quantiles of the predicted values, prices each segment at the single price that is most profitable for that segment on the training half, and then evaluates on the other half: $K = 1$ is the uniform price, $K = 4$ coarse segmentation, and $K = 50$ fine segmentation, that is, "personalization."

```r
P2 <- list(N = 20000, beta = c(const = 3.2, sales = 0.45, tenure = 0.15, electronics = 0.2, home = -0.15),
           sigma = 0.5, cat_levels = c("apparel", "electronics", "home"), cat_probs = c(0.5, 0.3, 0.2))
set.seed(4102); C2 <- gen_C2(P2)           # whole-USD WTP; train marks the half whose WTP the survey recorded
run <- pd_run(C2, Ks = c(1, 4, 50))          # OLS of log WTP on the training half, K segments, evaluation half
print(run$ladder, digits = 4, row.names = FALSE)
wev <- C2$wtp[!C2$train]; vals <- sort(unique(wev))
bb <- bbm_segment(vals, as.numeric(table(factor(wev, levels = vals))) / length(wev))
round(c(uniform_profit = bb$uniform_profit, efficient = bb$efficient_total,
        cs_max = bb$efficient_total - bb$uniform_profit, segments = length(bb$segments$weight)), 3)
set.seed(4104)                               # Monte Carlo: 200 fresh populations and splits
draws <- replicate(200, {
  L <- pd_run(gen_C2(P2), c(1, 4, 50))$ladder
  c(profit_fine = L$profit[3] / L$profit[1] - 1, cs_fine = L$cs[3] / L$cs[1] - 1,
    ts_fine = L$ts[3] / L$ts[1] - 1, served_fine = L$served[3] - L$served[1])
})
round(t(apply(draws, 1, function(x) c(mean = mean(x), lo = quantile(x, 0.025, names = FALSE),
                                      hi = quantile(x, 0.975, names = FALSE), share_pos = mean(x > 0)))), 3)
#   K profit    cs    ts served mean_price_paid
#   1  13.21 12.54 25.75 0.5504           24.00
#   4  15.23 11.78 27.01 0.6555           23.23
#  50  15.61 11.63 27.23 0.6858           22.76
# uniform_profit      efficient         cs_max       segments
#         13.248         32.348         19.100        180.000
#               mean     lo     hi share_pos
# profit_fine  0.193  0.174  0.212     1.000
# cs_fine     -0.124 -0.197 -0.055     0.000
# ts_fine      0.033 -0.006  0.069     0.955
# served_fine  0.120  0.070  0.157     1.000
```

In this draw, fine segmentation's profit is 18.2% higher than the uniform price's, its consumer surplus lower by 7.3%, and the share served rises from 0.550 to 0.686; merchant by merchant, 60.6% are offered a price below the uniform price, 38.1% are better off, and 18.8% are newly served. Over 200 redraws of the population and the split, fine segmentation's profit gain averages 19.3% (the middle 95% of replications runs from 17.4% to 21.2%), consumer surplus falls 12.4% on average ($-19.7\%$ to $-5.5\%$), and total surplus usually rises: by 3.3% on average, positive in 95.5% of replications, with a range of $-0.6\%$ to 6.9%. Coarse and fine segmentation fall inside the triangle of Theorem 4, and the uniform-price point lies slightly below the floor because its price is set out of sample (the figure in Section 3.3).

### 6.4 C3: Bundling

Fifty thousand consumers value 100 information goods independently, each valuation drawn from the uniform distribution on $[0, 1]$, with zero marginal cost; for $J \in \{1, 2, 5, 20, 100\}$ we compute the optimal pure-bundle price:

```r
J <- c(1, 2, 5, 20, 100)
set.seed(4103); B <- gen_C3(list(N = 50000, J = J))       # bundle values: sums of independent U[0, 1] values
round(t(sapply(seq_along(J), function(k) {
  v <- B[, k]; b <- bundle_opt(v)                           # optimal pure-bundle price at zero marginal cost
  c(J = J[k], rev_per_good = b$rev / J[k], dwl_share = b$dwl / mean(v), share_buy = mean(v >= b$price))
})), 3)
#        J rev_per_good dwl_share share_buy
# [1,]   1        0.249     0.254     0.497
# [2,]   2        0.272     0.189     0.657
# [3,]   5        0.311     0.113     0.820
# [4,]  20        0.375     0.060     0.922
# [5,] 100        0.432     0.026     0.970
```

Revenue per good climbs from a quarter toward a half, and the DWL share falls from a quarter toward zero, just as Theorem 5 predicts; convergence is not fast, and with 100 goods bundled, revenue per good still reaches only 86.5% of the mean.

### 6.5 C4: Price Leadership

```r
d <- c(4, 3); a <- 1; cc <- c(2, 2)                         # the C1 baseline: Cirrus first, Garnet second
P <- rbind(simultaneous = bertrand_logit(d, a, cc), cirrus_leads = price_leader(d, a, cc),
           garnet_leads = rev(price_leader(rev(d), a, rev(cc))), joint = joint_max_logit(d, a, cc))
X <- cbind(P, t(apply(P, 1, logit_profit, delta = d, alpha = a, c = cc)))
colnames(X) <- c("p_Cirrus", "p_Garnet", "profit_Cirrus", "profit_Garnet"); round(X, 4)
#              p_Cirrus p_Garnet profit_Cirrus profit_Garnet
# simultaneous   3.7461   3.3179        0.7461        0.3179
# cirrus_leads   3.8067   3.3261        0.7468        0.3261
# garnet_leads   3.7649   3.4250        0.7649        0.3192
# joint          4.1626   4.1626        0.8499        0.3127
```

Under both orders of moves, both prices are above simultaneous pricing, which is Theorem 13; that they are all below joint pricing is this example's number, not a conclusion of the theorem; profits are per unit of market size.

### 6.6 Reconciling the Accounts

The 25 of Section 1, caught in the middle, now has its place. With two types of merchants a single package earns at most 25, the menu lifts it to 26, and complete information could reach 34: the upper bound is set by information. With homogeneous merchants, the two-part tariff sweeps up the entire surplus of 50: the upper bound is set by the pricing instrument. Once Garnet enters, the Cournot price is 5.33 and industry profit 22.2, Hotelling differentiation gives a price of 3, and undifferentiated Bertrand wipes out profit: the lower bound is set by competition. Split into a vertical chain, double marginalization pushes the retail price above the integrated monopoly price, and combined profit falls to 18.75. C1 to C4 add the three things a hand calculation cannot hold: the mode of competition is written in the sign of cross pass-through, not in own pass-through; over 200 redraws, personalized pricing raises profit by nearly a fifth on average and lowers consumer surplus by more than a tenth, and total surplus can rise only if the number served rises; and both bundling and the timing of pricing can raise revenue or prices without collusion.

## 7 Failure Modes and Robustness

The welfare of price discrimination has no universal sign. The only reliable direction is that profit does not fall, because the firm can always revert to a uniform price. Beyond that, if total output does not rise, total welfare does not rise (Theorem 3), and first-degree discrimination is the most efficient yet the least favorable to consumers. C2 and Dubé and Misra (2023) move total surplus in opposite directions, consistent with Theorem 3: one expanded transactions, the other lowered conversion. Judging personalized pricing means computing output, welfare, and distribution separately.

The Bertrand paradox reminds us that the model is extremely sensitive to its assumptions: zero profit rests entirely on homogeneity, no capacity constraint, and one-shot play; loosen any one of them and the price rises back substantially. When analyzing platform competition, first asking in which dimension the firms compete, whether there are capacity constraints, and whether they coexist over the long run matters far more than applying a ready-made model.

Conduct cannot be read directly from prices. A high price may come from collusion, differentiation, or capacity constraints; own pass-through cannot separate modes of competition either: in C1 the true values of Cirrus's own pass-through under Bertrand and joint pricing are 0.592 and 0.613; when marginal cost is backed out under the wrong conduct, C1's true cost of USD 2 is read as 2.515 and 2.994. Identification takes cross pass-through, demand rotation, or institutional changes (Sections 3.3 and 3.4 of Chapter 29).

Two-sided markets can turn single-sided intuition entirely on its head: a platform connects buyers and sellers, the size of one side affects the value of the other, and optimal pricing may charge one side below marginal cost; Section 2.1 of Chapter 32 gives the pricing formula under cross-side network effects. The problem of durable goods and commitment is especially real for digital products, and whether a credible commitment can be established often matters more for the outcome than the static monopoly formula does (Theorem 6). Finally, C2's training half observes willingness to pay directly, a convenience of this example; in reality, how good personalized prices are depends on how the demand model extrapolates where data are sparse.

This chapter's conclusions all hang on specific assumptions: the welfare of discrimination depends on output and demand curvature; Bertrand's zero profit depends on homogeneity and no capacity constraint; collusion and competition are hard to tell apart in the data; and single-sided pricing reverses in a two-sided world. Robust pricing analysis means recognizing your own market's dimension of competition, information structure, and commitment power, then choosing the model, and spelling out which way the conclusions move when the assumptions are relaxed. With this, the microeconomic-theory foundation of Part I is laid, and the identification and structural chapters to come will keep returning here to recover the behavioral primitives they set out to estimate.

## 8 Further Reading

::: {.readings}
Required reading, in suggested reading order:

- Mas-Colell, Whinston and Green (1995, Microeconomic Theory), Chapter 12 (Market Power), with second-degree discrimination as screening in Section 14.C. The authoritative reference on monopoly, price discrimination, and oligopoly, and the source of Theorems 1, 2, 7 to 9, and 11 of this chapter.
- Nolan Miller (2006, Notes on Microeconomic Theory), Chapter 9 (Monopoly). The intuition-building groundwork for monopoly pricing and price discrimination, more patient than MWG, good to read first.
- Tirole (1988, The Theory of Industrial Organization), Chapters 1 to 9. The standard graduate textbook of industrial organization; every section of this chapter can be taken further there.

Further reading:

- Bergemann, Brooks and Morris (2015). The limits of price discrimination, the original source of Theorem 4, to be read against the information design of Chapter 3.
- Schmalensee (1981), Varian (1985), and Aguirre, Cowan and Vickers (2010). The output test for third-degree price discrimination and the demand-curvature conditions: the source of Theorem 3 and its extension.
- Dubé and Misra (2023). A field experiment on personalized pricing and its welfare, to be read against block C2.
- Leslie (2004) and Porter (1983). Two empirical classics on price discrimination and collusion, on which Sections 5.1 and 5.2 are based.
- Gul, Sonnenschein and Wilson (1986) and Coase (1972). Durable-goods monopoly and the Coase conjecture, the original source of Theorem 6; Stokey (1979) proves that with equal discounting, costs that do not fall, and commitment, the optimum is a constant price.
- Oi (1971), Adams and Yellen (1976), and Bakos and Brynjolfsson (1999). Two-part tariffs, bundling, and large bundles, on which Section 3.4 is based.
- Shapiro and Varian (1999, Information Rules) and Salant (1989). Versioning of information goods, on which Section 3.5 is based.
- Kreps and Scheinkman (1983) and d'Aspremont, Gabszewicz and Thisse (1979). Price competition under capacity constraints, and the existence correction to the Hotelling model.
- Caplin and Nalebuff (1991) and Aguirregabiria (2021). Existence of equilibrium in differentiated price competition (the source of the existence claim in Theorem 10), and a derivation of the logit markup to check against.
- Weyl and Fabinger (2013) and Miklós-Thal and Shaffer (2021). General principles of pass-through, and a comment on their range of validity; Bresnahan (1982) is the origin of conduct identification.
- Brown and MacKay (2023). How algorithms that reprice at different frequencies raise prices, on which Section 4.9 is based.
- Spengler (1950), Hotelling (1931), and Hotelling (1929). The original arguments for double marginalization, the rent path of exhaustible resources, and the linear city.
- Dixit (1980), Milgrom and Roberts (1982), and Fudenberg and Tirole (1984). The three cornerstones of entry deterrence: capacity commitment, the signaling version of limit pricing, and the taxonomy of investment postures, on which Section 4.7 is based.
:::

::: {.apa-refs}
- Adams, W. J., & Yellen, J. L. (1976). Commodity bundling and the burden of monopoly. *The Quarterly Journal of Economics, 90*(3), 475-498. https://doi.org/10.2307/1886045
- Aguirre, I., Cowan, S., & Vickers, J. (2010). Monopoly price discrimination and demand curvature. *American Economic Review, 100*(4), 1601-1615. https://doi.org/10.1257/aer.100.4.1601
- Aguirregabiria, V. (2021). *Empirical industrial organization: Models, methods, and applications* [Unpublished manuscript]. University of Toronto.
- Bakos, Y., & Brynjolfsson, E. (1999). Bundling information goods: Pricing, profits, and efficiency. *Management Science, 45*(12), 1613-1630. https://doi.org/10.1287/mnsc.45.12.1613
- Bergemann, D., Brooks, B., & Morris, S. (2015). The limits of price discrimination. *American Economic Review, 105*(3), 921-957. https://doi.org/10.1257/aer.20130848
- Bresnahan, T. F. (1982). The oligopoly solution concept is identified. *Economics Letters, 10*(1-2), 87-92. https://doi.org/10.1016/0165-1765(82)90121-5
- Brown, Z. Y., & MacKay, A. (2023). Competition in pricing algorithms. *American Economic Journal: Microeconomics, 15*(2), 109-156. https://doi.org/10.1257/mic.20210158
- Caplin, A., & Nalebuff, B. (1991). Aggregation and imperfect competition: On the existence of equilibrium. *Econometrica, 59*(1), 25-59. https://doi.org/10.2307/2938239
- Coase, R. H. (1972). Durability and monopoly. *Journal of Law and Economics, 15*(1), 143-149. https://doi.org/10.1086/466731
- d'Aspremont, C., Gabszewicz, J. J., & Thisse, J.-F. (1979). On Hotelling's "Stability in competition". *Econometrica, 47*(5), 1145-1150. https://doi.org/10.2307/1911955
- Dixit, A. (1980). The role of investment in entry-deterrence. *The Economic Journal, 90*(357), 95-106. https://doi.org/10.2307/2231658
- Dubé, J.-P., & Misra, S. (2023). Personalized pricing and consumer welfare. *Journal of Political Economy, 131*(1), 131-189. https://doi.org/10.1086/720793
- Fudenberg, D., & Tirole, J. (1984). The fat-cat effect, the puppy-dog ploy, and the lean and hungry look. *American Economic Review, 74*(2), 361-366 (Papers and Proceedings).
- Green, E. J., & Porter, R. H. (1984). Noncooperative collusion under imperfect price information. *Econometrica, 52*(1), 87-100. https://doi.org/10.2307/1911462
- Gul, F., Sonnenschein, H., & Wilson, R. (1986). Foundations of dynamic monopoly and the Coase conjecture. *Journal of Economic Theory, 39*(1), 155-190. https://doi.org/10.1016/0022-0531(86)90024-4
- Hotelling, H. (1929). Stability in competition. *The Economic Journal, 39*(153), 41-57. https://doi.org/10.2307/2224214
- Hotelling, H. (1931). The economics of exhaustible resources. *Journal of Political Economy, 39*(2), 137-175. https://doi.org/10.1086/254195
- Kreps, D. M., & Scheinkman, J. A. (1983). Quantity precommitment and Bertrand competition yield Cournot outcomes. *The Bell Journal of Economics, 14*(2), 326-337. https://doi.org/10.2307/3003636
- Leslie, P. (2004). Price discrimination in Broadway theater. *The RAND Journal of Economics, 35*(3), 520-541. https://doi.org/10.2307/1593706
- Mas-Colell, A., Whinston, M. D., & Green, J. R. (1995). *Microeconomic theory*. Oxford University Press.
- Miklós-Thal, J., & Shaffer, G. (2021). Pass-through as an economic tool: On exogenous competition, social incidence, and price discrimination. *Journal of Political Economy, 129*(1), 323-335. https://doi.org/10.1086/711348
- Milgrom, P., & Roberts, J. (1982). Limit pricing and entry under incomplete information: An equilibrium analysis. *Econometrica, 50*(2), 443-459. https://doi.org/10.2307/1912637
- Miller, N. H. (2006). *Notes on microeconomic theory* [Unpublished lecture notes]. Harvard University.
- Oi, W. Y. (1971). A Disneyland dilemma: Two-part tariffs for a Mickey Mouse monopoly. *The Quarterly Journal of Economics, 85*(1), 77-96. https://doi.org/10.2307/1881841
- Porter, R. H. (1983). A study of cartel stability: The Joint Executive Committee, 1880-1886. *The Bell Journal of Economics, 14*(2), 301-314. https://doi.org/10.2307/3003634
- Salant, S. W. (1989). When is inducing self-selection suboptimal for a monopolist? *The Quarterly Journal of Economics, 104*(2), 391-397. https://doi.org/10.2307/2937854
- Schmalensee, R. (1981). Output and welfare implications of monopolistic third-degree price discrimination. *American Economic Review, 71*(1), 242-247.
- Shapiro, C., & Varian, H. R. (1999). *Information rules: A strategic guide to the network economy*. Harvard Business School Press.
- Spengler, J. J. (1950). Vertical integration and antitrust policy. *Journal of Political Economy, 58*(4), 347-352. https://doi.org/10.1086/256964
- Stokey, N. L. (1979). Intertemporal price discrimination. *The Quarterly Journal of Economics, 93*(3), 355-371. https://doi.org/10.2307/1883163
- Tirole, J. (1988). *The theory of industrial organization*. MIT Press.
- Varian, H. R. (1985). Price discrimination and social welfare. *American Economic Review, 75*(4), 870-875.
- Weyl, E. G., & Fabinger, M. (2013). Pass-through as an economic tool: Principles of incidence under imperfect competition. *Journal of Political Economy, 121*(3), 528-583. https://doi.org/10.1086/670401
:::
