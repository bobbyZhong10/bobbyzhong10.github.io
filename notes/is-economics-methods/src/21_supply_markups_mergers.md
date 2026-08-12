---
title: "Supply Side: Markups, Conduct, Mergers & Welfare"
subtitle: "Recovering Marginal Costs and Simulating Mergers"
seriesline: "Foundations of Information Systems Economics · Chapter 21"
series: "Foundations of Information Systems Economics: Empirical Methods for Information Systems Economics"
sidebartitle: "Foundations of Information Systems Economics"
sidebarsub: "Theory · Methods · Applications"
chapterlabel: "Chapter 21 · Supply Side: Markups, Conduct, Mergers & Welfare"
---

## Introduction

Two meal-kit brands are about to merge. They claim that consolidating delivery will save money, and the regulator worries that common ownership will push prices up. Both sides talk about "efficiency" and "competition," but the number that matters most, the marginal cost of each box, never shows up in the transaction data. A researcher can see prices and quantities, yet has to judge how much markup firms are earning and how they will reprice after the merger. This sounds like reading a cash-register receipt to infer what is going on in the kitchen, and the supply-side structural model is precisely what makes that inference testable.

The pivot of the whole argument is the firm's pricing first-order condition. Given the slope of demand, the pattern of product ownership, and the mode of competition, the optimal price locks markup and marginal cost together, so you can back out cost from price, then change ownership or cost and solve for the new equilibrium price. The previous chapter recovered demand; this chapter recovers cost. But there is one extra economic assumption here that no algebra can substitute away: do firms compete as multiproduct Bertrand players, or do they coordinate in some other way? Any conduct can be made to imply a set of costs, and the mere fact that the equations solve does not mean the story holds.

The Selva case puts this danger out in the open. If you wrongly assume collusion, the model hands back a marginal cost of negative eight dollars per box; if you simulate the merger of brand 1 and brand 2 with IIA-constrained logit demand, you underestimate the price increase by a quarter because you fail to recognize that they are close neighbors. We will back out costs in turn, compute UPP and GUPPI, solve for post-merger prices and welfare, and compare the conclusions delivered by different conduct assumptions and demand models. The efficiency defense will also be written as an explicit question: by how much must marginal cost fall to offset the unilateral upward pricing pressure of the merger?

A credible merger analysis, then, has to lay the entire chain of reasoning on the table: demand estimation determines substitution, substitution determines the incentive to raise prices, the first-order condition turns price into marginal cost, and all of this feeds jointly into the counterfactual and the welfare accounting. Whether the costs are reasonable, whether the mode of competition passes a test, and whether the conclusion changes with the demand model are themselves results, not technical footnotes tucked into an appendix. That a model can generate a counterfactual is nothing special; what is rare is that every layer of behavioral assumption the counterfactual rests on can withstand interrogation.

## 1 A Number That Does Not Add Up

::: {.case}
Selva's two meal-kit brands, brand 1 and brand 2, have manufacturers who want to merge. Each brand holds under 10% share in every market, and together they add up to about 20%, so a mechanical concentration screen might well wave them through. But to judge whether the merger actually harms consumers, you need to know each brand's marginal cost, and marginal cost is not in the data. All we can do is back it out from pricing behavior: firms are pricing optimally, price equals marginal cost plus a markup determined by demand elasticity, and plugging the demand estimated in the previous chapter into that relation gives us the cost.
:::

Backing out cost requires first assuming how firms compete. Try one assumption first: these six brands secretly collude, pricing jointly like a single monopolist to maximize total profit. Under this assumption, each brand recognizes when it sets its price that raising it pushes customers onto all the other brands (including its co-conspirators), so the internalized markup is enormous. Subtract this large markup from price, and the median backed-out marginal cost is negative $7.89$ dollars per box, with $99.8\%$ of products yielding negative cost. Negative marginal cost means that a firm not only spends nothing to produce a unit but earns money on it, earning more the more it produces, which is physically impossible.

This absurd number is not a computational error but a signal: the collusion conduct assumption is wrong. Switch to a different assumption, Bertrand competition where each firm goes it alone, and each brand considers only the effect of its own price increase on its own sales; the internalized markup is much smaller, and the median backed-out marginal cost is $5.94$ dollars per box, all positive and hugging the true value of $5.81$ (because this data was in fact generated under Bertrand competition). Same prices, same demand, two conduct assumptions delivering wildly different costs, one absurd and one reasonable.

This is where supply-side analysis differs most from the demand estimation of the previous chapter, and it is the through-line of this chapter. The first-order condition can always solve for marginal cost, since the number of equations exactly equals the number of unknowns, so whether or not the conduct is right, you can always back out a set of numbers. Whether those numbers are reasonable is the first gate for judging whether a conduct assumption is right. Section 3 will explain thoroughly the point that "conduct cannot be identified from the first-order condition itself and must rely on extra structural testing," and Section 6 will use a formal test to distinguish Bertrand from collusion. Once you have selected the correct conduct and backed out the costs, simulating mergers, evaluating welfare, and computing antitrust screening indices all follow naturally, and these make up the second half of the chapter.

## 2 The Economic Model and the Estimand

Before backing out cost, let us write down the supply-side model and the target quantity clearly. At its core is the optimal pricing first-order condition of a multiproduct firm.

There are three kinds of target quantities this chapter has to deliver. The first is marginal cost, which is not in the data and must be backed out from pricing behavior. The second is counterfactual prices and welfare, chiefly the post-merger equilibrium prices, the change in consumer surplus, and the change in producer surplus. The third is conduct itself, that is, how firms actually compete, which is both the premise for backing out cost and something that can in turn be tested. These three quantities cascade: conduct plus demand pins down cost, cost plus a new ownership structure pins down post-merger prices, and the price changes pin down welfare.

Start with the firm's pricing problem. Firm $f$ owns a set of products $\mathcal{J}_f$ and chooses the prices of these products to maximize joint profit:

$$\max_{\{p_j\}_{j \in \mathcal{J}_f}}\ \sum_{j \in \mathcal{J}_f} (p_j - mc_j)\, s_j(\mathbf{p})$$

Taking the first-order condition for each owned product $j$:

$$s_j(\mathbf{p}) + \sum_{k \in \mathcal{J}_f} (p_k - mc_k)\, \frac{\partial s_k}{\partial p_j}(\mathbf{p}) = 0$$

The intuition of this equation is that raising the price of $j$ a little earns a good marginal profit on that unit (the first term), but drives some customers away: those flowing to rivals are a pure loss, while those flowing to the firm's other products $k$ are recaptured, still earning the markup $p_k - mc_k$, so the net cost of raising the price shrinks (the same-firm cross terms in the summation). A single-product firm has only the $k = j$ term and no such recapture; a multiproduct firm dares to set a higher markup precisely because it internalizes the mutual cannibalization among its own products.

Stacking the first-order conditions of all products into matrix form gives the compact expression this chapter uses again and again. Define the negative matrix of demand price derivatives $\Delta$, with elements $\Delta_{jk} = -\partial s_j / \partial p_k$ (the own terms are on the diagonal and positive), then define the ownership matrix $\mathcal{H}$, with $\mathcal{H}_{jk} = 1$ when products $j$ and $k$ belong to the same firm and 0 otherwise. Then the first-order conditions of all firms together are

$$s(\mathbf{p}) = \big(\mathcal{H} \odot \Delta(\mathbf{p})\big)\,(\mathbf{p} - mc)$$

where $\odot$ is element-wise multiplication (the Hadamard product), and $\mathcal{H} \odot \Delta$ zeroes out the cross terms of products that do not belong to the same firm. Solving for the markup gives the formula for backing out marginal cost:

$$\boxed{\ mc = \mathbf{p} - \big(\mathcal{H} \odot \Delta(\mathbf{p})\big)^{-1} s(\mathbf{p})\ }$$

This is the pivot of the whole chapter. It says that marginal cost equals price minus markup, and the markup is determined entirely by the demand derivatives (the $\Delta$ estimated in the previous chapter) and the ownership structure $\mathcal{H}$. The $\Delta$ comes from demand estimation, the $\mathcal{H}$ comes from who owns what plus a conduct assumption, and putting the two together backs cost out. Note that cost is backed out, not estimated: given $\Delta$ and $\mathcal{H}$, this is $J$ equations solving for $J$ unknowns, and there is always a unique solution. This "always solvable" property is both a convenience and the trap that Section 3 warns about.

::: {.intuition}
The ownership matrix $\mathcal{H}$ is the switch that fits different market structures into a single formula. With single-product firms, $\mathcal{H}$ is the identity matrix, each firm minds only its own single product, and markups are lowest. With multiproduct firms, $\mathcal{H}$ is a block 0-1 matrix, cannibalization among same-firm products is internalized, and markups are higher. With full collusion, $\mathcal{H}$ is a matrix of all ones, all products price like a single monopolist, and markups are highest. A merger, in this framework, is just flipping the relevant zeros of $\mathcal{H}$ for the two firms into ones: before the merger they ignore each other, after the merger they internalize each other's cannibalization, and so both want to raise prices. The counterfactual analysis of the entire chapter is essentially changing this $\mathcal{H}$ or changing $mc$ and then re-solving the first-order condition.
:::

To summarize: the supply-side model is multiproduct Bertrand optimal pricing, the estimand is marginal cost, counterfactual prices and welfare, and conduct; the core formula $mc = \mathbf{p} - (\mathcal{H} \odot \Delta)^{-1} s$ combines demand derivatives and ownership structure into cost; cost is backed out rather than estimated, and the first-order condition always has a solution, which is the starting point for the identification discussion in the next section.

## 3 Identification

Identification logically precedes computation. This section answers: given prices, shares, and the demand estimated in the previous chapter, what exactly can the supply side identify, and what can it in principle not identify and must instead fill in with extra assumptions or tests. This is the most detailed section of the chapter's analysis, split into numbered subsections, each of which gives the intuition of a difficulty before formalizing it.

### 3.1 Markup Recapture: Cost Is Backed Out Exactly Given Conduct

Build the intuition from the simplest single-product case first. A firm selling only one product has a first-order condition that degenerates to $s_j + (p_j - mc_j)\, \partial s_j / \partial p_j = 0$, which rearranges to

$$p_j - mc_j = \frac{s_j}{-\partial s_j / \partial p_j} = \frac{p_j}{|\epsilon_{jj}|}$$

This is Lerner's inverse-elasticity pricing rule: the ratio of markup to price equals the reciprocal of the own-price elasticity. The less elastic demand is (the smaller $|\epsilon_{jj}|$), the higher the markup. In Selva the median Lerner index is $0.23$, corresponding to an own-price elasticity of about $-4.4$ and a median markup of about $1.65$ dollars. Subtracting the markup from price gives marginal cost, backed out one number at a time.

The multiproduct case just generalizes this logic. When a firm prices $j$, on top of the inverse-elasticity markup on that unit, it adds one more term: the share of customers driven away by the price increase who flow to the firm's other products $k$ is an internalized opportunity cost. Writing it with the diversion ratio $D_{j \to k}$ (the fraction of customers lost by $j$ that are captured by $k$), the multiproduct pricing rule is the single-product markup plus $\sum_{k \neq j, k \in \mathcal{J}_f} (p_k - mc_k) D_{j \to k}$. The more products a firm owns and the closer neighbors they are (the larger $D_{j \to k}$), the higher this internalized opportunity cost and the higher the markup. The matrix formula from Section 2, $mc = \mathbf{p} - (\mathcal{H} \odot \Delta)^{-1} s$, is exactly the compact form of this logic, with $\Delta$ holding all the price derivatives and $\mathcal{H}$ deciding which cannibalization is internalized.

The crux is that, given the demand derivatives $\Delta$ and the ownership-plus-conduct $\mathcal{H}$, this backing-out is exact and unique, $J$ equations solving for $J$ unknowns. This is a convenience, because cost is immediately available; it is also a trap, because whether or not $\mathcal{H}$ is set correctly, you can always back out a set of costs. Cost will not tell you whether the conduct is right, unless the backed-out cost is itself absurd (negative, say), a point Section 3.3 develops.

### 3.2 Diversion Ratio: The Demand Product the Supply Side Truly Relies On

The previous chapter kept stressing that the diversion ratio is the core product demand estimation has to deliver, and only on the supply side does it become clear why. The substitution rate from $j$'s price increase to $k$ is defined as

$$D_{j \to k} = \frac{\partial s_k / \partial p_j}{-\partial s_j / \partial p_j}$$

It is the fraction of each unit of sales $j$ loses that is captured by $k$, and also equals "the conditional probability of choosing $k$ among those who do not choose $j$," that is, the share of second choice. Multiproduct pricing, merger analysis, and UPP screening all rely on it and nothing else. The reason for using diversion rather than cross elasticity is that shares sum to one, so diversion is comparable across products and across markets while cross elasticity is not. The previous chapter offered a practical slogan: anywhere you are tempted to report a cross elasticity, report diversion instead.

Diversion is also exactly where the previous chapter's demand-model choice echoes on the supply side. The IIA of plain logit forces $D_{j \to k} = s_k / (1 - s_j)$, which depends only on shares and has nothing to do with "who resembles whom." In Selva, brand 1 and brand 2 are quality neighbors, the true $D_{1 \to 2}$ is $0.271$, BLP estimates $0.278$, while logit gives only $0.206$. Section 6 will show that it is precisely this diversion, underestimated by logit, that makes logit underestimate the merger's price increase by a quarter. Every supply-side conclusion, for better or worse, hangs on whether the previous chapter estimated diversion correctly.

### 3.3 Conduct Cannot Be Identified from the First-Order Condition Itself

Now develop the "always solvable" trap of Section 3.1 fully; it is the deepest identification problem on the supply side. Given the observed prices and shares, and given the demand derivatives $\Delta$, for any conduct assumption $\mathcal{H}(\kappa)$ you can back out a set of marginal costs $mc(\kappa) = \mathbf{p} - (\mathcal{H}(\kappa) \odot \Delta)^{-1} s$. Bertrand backs out one set, collusion backs out another, some partial internalization backs out a third, and each set makes the observed price exactly optimal under the corresponding conduct. So from prices and shares alone, you cannot judge how firms actually compete, and conduct is not identifiable from the first-order condition itself. This is the same kind of problem as in earlier chapters, where exogeneity cannot be tested from within OLS and parallel trends cannot be tested from within DiD: the model itself cannot verify the assumption that supports it.

Then why could Section 1 see at a glance that collusion was wrong? Because the cost it backed out was absurd, negative eight dollars, negative for ninety-nine percent. This is a kind of external plausibility check: the backed-out cost must satisfy economic common sense (nonnegative, a reasonable relationship with known cost shifters). But this check is soft, able only to rule out grossly outlandish conduct and unable to distinguish two assumptions that both yield positive costs. To truly identify conduct, you need more structure.

To distinguish conduct, the classic idea is Bresnahan (1982)'s demand rotation: different conducts respond differently to the same demand curve, and an exogenous variable that rotates the demand curve about some point (a demand rotator) pulls apart the optimal pricing under different conducts, thereby distinguishing them. The modern approach turns this into a formal test: given an instrument that affects only cost or only rotates demand, without entering the other side directly, see under which conduct assumption the backed-out cost can be better explained by a cost model. Section 6 will use one such test to cleanly distinguish Bertrand from collusion. The core idea to establish here is that conduct is an assumption that needs a positive defense or a formal test, not a setting you can adopt by default; defaulting to Bertrand without a test is as groundless as defaulting to collusion.

### 3.4 What Variation Testing Conduct Requires

Since conduct must be identified through external variation, we should be clear about what variation that is. The key is to find an instrument that makes markups under different conducts move differently while not entering the marginal cost equation directly. Such variation has two kinds of source. One is a demand rotator, a variable that changes the shape of the demand curve (not just shifts it), such as a demographic change that alters the composition of consumers; it rotates demand, and different conducts respond differently to the rotation. The other is a cost shifter, a variable that moves the cost of some products and thereby moves rivals' markups through competition, such as a rival's input cost; it transmits differently to a given product's markup under different conducts.

Modern conduct tests turn this idea into model selection. Rivers and Vuong (2002) give a comparison test for non-nested models, computing a supply-side goodness of fit (residual sum of squares or a GMM objective) for two candidate conducts and testing whether their difference is significant. Duarte, Magnolfi, Sølvsten and Sullivan (2024) systematize this framework, pointing out that using the wrong conduct amounts to omitting a variable in the supply equation (the difference between the two conducts' markups), that the optimal instrument is exactly the part of this markup difference that can be predicted beyond the cost shifters, and prove that in structural estimation the Rivers-Vuong type test performs far better than other schemes. Their common identification condition is that the instrument must, beyond the cost shifters, still predict the two conducts' markup difference (relevance) and not enter the cost equation (exogeneity); if the markup difference of the two conducts happens to be fully explained by the cost shifters, there is no way to distinguish them. The test in Section 6 is a minimal implementation of this logic.

### 3.5 How Demand Misspecification Transmits to the Supply Side

The last identification problem lies not inside the supply side but at its interface with the previous chapter. Everything on the supply side is built on the demand derivatives $\Delta$; if $\Delta$ is wrong, the backed-out cost, the simulated merger, and the computed welfare are all wrong, and wrong systematically. If the previous chapter estimated demand with logit, IIA underestimates the diversion between neighboring products, so the internalized opportunity-cost term in the multiproduct markup is underestimated, the backed-out cost is biased upward, and the merger simulation, having underestimated the mutual substitution between the merging parties, underestimates the price increase. This is not random error but a systematic bias that the single misspecification of IIA amplifies all the way down the first-order condition. Section 6 will place the logit and BLP merger simulations side by side to see how large this bias is. The supply side is therefore extremely sensitive to the choice of demand model, and getting demand right is not something the previous chapter finishes on its own; its payoff is only fully realized on the supply side.

The main points of this section can be summarized as follows: given demand and conduct, marginal cost is backed out exactly by the first-order condition, and markup is determined by diversion and ownership; conduct itself is not identifiable from the first-order condition, the plausibility of the backed-out cost is a soft check, and a formal distinction relies on exogenous variation like a demand rotator or cost shifter plus a Rivers-Vuong type test; and the entire supply side is extremely sensitive to the correctness of the previous chapter's demand (especially diversion), with misspecification transmitting systematically. None of this has yet said how to compute, which is the next section.

## 4 Estimation and Counterfactuals

This section brings the identification logic of Section 3 down to operational computation: backing out cost, simulating mergers, first-order screening, welfare measurement, and testing conduct.

### 4.1 Backing Out Marginal Cost

With the demand from the previous chapter in hand, backing out cost is just computing $mc = \mathbf{p} - (\mathcal{H}^0 \odot \Delta(\mathbf{p}))^{-1} s(\mathbf{p})$ once at the observed prices, where $\mathcal{H}^0$ is the pre-merger ownership structure. $\Delta(\mathbf{p})$ is obtained by differentiating the estimated demand parameters at the observed prices; for a random-coefficients model the derivative is the integral of individual derivatives over the taste distribution. In Selva, using BLP demand and assuming single-product Bertrand, the median backed-out marginal cost is $5.94$ dollars, all positive, hugging the true value of $5.81$. This set of costs is the raw material for all the counterfactuals that follow.

### 4.2 Simulating a Merger

Simulating a merger takes four steps. Step one, back out the pre-merger marginal cost (Section 4.1). Step two, change the ownership matrix from $\mathcal{H}^0$ to $\mathcal{H}^1$, flipping the merging parties' cross terms from 0 to 1, so the merged firm internalizes each other's cannibalization. Step three, holding marginal cost and demand fixed, re-solve the first-order condition for the new equilibrium prices $\mathbf{p}^{\text{post}}$. Step four, the price change is the merger's price effect, from which welfare is then computed.

Step three is the only technical difficulty, because price appears on both sides of the equation (both share and derivative depend on price). Naively iterating $\mathbf{p} \leftarrow mc + (\mathcal{H}^1 \odot \Delta(\mathbf{p}))^{-1} s(\mathbf{p})$ is not guaranteed to converge and may diverge. Morrow and Skerlos (2011) give a numerically far more robust fixed point. Decompose $\Delta$ into the diagonal $\Lambda$ (each product's own price sensitivity, summed) and the rest $\Gamma$, $\Delta = \Lambda - \Gamma$; the first-order condition satisfied by the markup $\eta = \mathbf{p} - mc$ can be rewritten as

$$\eta \leftarrow \Lambda(\mathbf{p})^{-1}\big[\, s(\mathbf{p}) + (\mathcal{H}^1 \odot \Gamma(\mathbf{p}))\, \eta \,\big], \qquad \mathbf{p} = mc + \eta$$

This so-called $\zeta$-map is a contraction mapping, converging from any starting value, and is the default solution method in mature software like pyblp. This chapter uses it to solve for post-merger prices.

Simulating the merger of brand 1 and brand 2's manufacturers on Selva, with BLP demand, brand 1's price rises $4.8\%$ after the merger, brand 2's rises $7.3\%$, and the four non-merging rival brands barely move (rising about $0.3\%$). Brand 2 rises more because it starts with a higher markup rate and absorbs more customers from brand 1. That the rivals barely move shows these two brands are not close neighbors of the other four, and the merger's price pressure concentrates within the merging pair. This is the typical picture of a differentiated-products merger: the harm concentrates on the merging parties and their nearest substitutes.

### 4.3 First-Order Screening: UPP and GUPPI

A full simulation has to re-solve the equilibrium, which is not cheap to compute. Antitrust practice often first uses a first-order index that does not re-solve the equilibrium as an initial screen. At its core is upward pricing pressure (UPP): the merger makes firm $j$ newly internalize, when pricing, the customers flowing to the merger partner's product, and this new opportunity cost is the pressure to raise prices,

$$UPP_j = \sum_{k \in \mathcal{J}_{f'}} (p_k - mc_k)\, D_{j \to k} - E_j$$

where $\mathcal{J}_{f'}$ is the merger partner's products and $E_j$ is the offset from cost efficiencies. Normalizing to price and dropping efficiencies gives GUPPI (the gross UPP index):

$$GUPPI_j = \frac{(p_k - mc_k)\, D_{j \to k}}{p_j}$$

It is dimensionless and is a threshold commonly used by antitrust agencies for initial screening. In Selva, brand 2's GUPPI is $0.074$ and brand 1's is $0.049$. Note that brand 2's GUPPI ($0.074$) is very close to the price increase for brand 2 from the full simulation ($7.3\%$), and this is no coincidence: GUPPI is exactly the first-order approximation of the upward pricing pressure, and when the pass-through rate is near 1 it approximately equals the equilibrium price increase. Mergers with high GUPPI are worth a full simulation, while low ones can be waved through in the initial screen.

### 4.4 Welfare

Once prices change, welfare can be computed. Consumer surplus is measured by the log-sum (inclusive value) of discrete choice, which is the result of Small and Rosen (1981): a consumer's expected maximum utility is $\log(1 + \sum_j e^{V_{ij}})$, divided by the marginal utility of price $\alpha_i$ to convert to dollars, then integrated over the taste distribution to give aggregate consumer surplus. The difference between the two log-sums before and after the merger is the change in consumer surplus:

$$\Delta CS = \int \frac{1}{\alpha_i}\Big[\log\big(1 + \textstyle\sum_j e^{V_{ij}^{\text{post}}}\big) - \log\big(1 + \textstyle\sum_j e^{V_{ij}^{\text{pre}}}\big)\Big]\, dF_i$$

Producer surplus is the sum of each product's markup times quantity $\sum_j (p_j - mc_j) s_j$, and the difference before and after the merger is the change in producer surplus. In Selva's merger, each consumer loses on average $0.22$ dollars of surplus (using true demand), while the merging parties' producer surplus rises, which is exactly the distributional consequence of a differentiated-products merger: consumers are hurt and the merging parties profit. The net change in total welfare depends on the two offsetting each other, and also on whether there are cost efficiencies.

### 4.5 The Efficiency Defense and Testing Conduct

Merging parties often defend with cost efficiencies: the merger raises prices but also lowers costs, so the net effect need not harm consumers. This can be quantified. Ask an inverse question: by how much must the merging parties' marginal cost fall to completely offset the price increase. In Selva, cutting the merging parties' marginal cost by $14.5\%$ returns brand 2's post-merger price to the pre-merger level. So whether the defense holds turns into a checkable question: can the merger really deliver a $14.5\%$ cost reduction? Turning the vague "efficiencies will offset the price increase" into a concrete threshold number is a major contribution of the structural model to antitrust analysis.

Last is testing conduct, the question Section 3 left open. The approach is to back out cost for each candidate conduct, then see which backed-out cost can be better explained by a cost model. Regress the backed-out cost on cost shifters (local input cost, product characteristics) and on rivals' cost shifters, and measure the fit with a GMM objective; the smaller the objective, the more the cost backed out under this conduct resembles the true cost structure. In Selva, Bertrand's GMM objective is four orders of magnitude smaller than collusion's (the cost backed out under collusion is negative for ninety-nine percent and simply cannot be explained by any reasonable cost model), and the test cleanly selects Bertrand. This is the minimal implementation of the Rivers-Vuong type conduct test: let the data pick the most defensible one out of several conduct assumptions, rather than defaulting to one.

The section's route is now complete: back out cost, change ownership and re-solve the first-order condition to simulate the merger, solve robustly with Morrow-Skerlos's $\zeta$-map, use UPP and GUPPI for first-order initial screening, measure welfare with the log-sum, quantify the efficiency threshold with an inverse problem, and test conduct with GMM fit.

## 5 Anchoring Papers

A method only stands firm once it lands in real research. Three anchoring papers, one bringing structural merger simulation into antitrust, one a modern exemplar of differentiated-products merger analysis, and one making the supply-side lifeline of diversion clear, each organized by the five elements of paper, method, data, results, and limitations, with the focus on how the identification assumptions are defended.

### 5.1 Nevo (2000)

::: {.case}
Paper: "Mergers with Differentiated Products: The Case of the Ready-to-Eat Cereal Industry," RAND Journal of Economics. It combines BLP-style structural demand with the supply-side first-order condition to demonstrate how to use estimated demand to simulate mergers of differentiated products, and is the landmark text bringing structural merger analysis into antitrust practice.

Method: first estimate cereal demand with random-coefficients logit (exactly the method of the previous chapter), then back out each brand's marginal cost with the multiproduct Bertrand first-order condition, then change the ownership matrix to simulate several actual or hypothetical mergers, computing post-merger prices and welfare. The method is of one lineage with this chapter.

Data: brand-level scanner data on ready-to-eat cereal across multiple cities and periods, containing shares, prices, and characteristics. The cereal industry, with highly differentiated brands and many multiproduct firms, is the classic battleground for merger analysis.

Results: the price effect of a merger depends heavily on whether the merging parties are close neighbors. Mergers of neighboring brands raise prices noticeably, while mergers of distant brands have little effect, consistent with the picture of this chapter's Selva. The paper shows how structural methods deliver far finer merger predictions than concentration indices.

Limitations: the conclusion is extremely sensitive to the correctness of demand estimation, and if diversion is estimated wrong the entire merger prediction is wrong; backing out cost relies on the Bertrand conduct assumption without a test; the simulation assumes the merger does not change the product portfolio or entry, ignoring long-run repositioning. These limitations are exactly what Sections 3 and 7 of this chapter stress repeatedly.
:::

### 5.2 Miller and Weinberg (2017)

::: {.case}
Paper: "Understanding the Price Effects of the MillerCoors Joint Venture," Econometrica. It studies the price effects after Miller and Coors formed a joint venture in the US beer market, and is the modern exemplar of differentiated-products merger analysis, especially for building a conduct test into the merger study.

Method: random-coefficients logit demand plus supply-side pricing, first estimate demand and back out cost, then test whether the post-merger firm's conduct changed from Bertrand to some coordinated pricing. The key innovation is not to presume conduct but to let the data choose between Bertrand and partial coordination, using the price change before and after the merger to identify the change in conduct.

Data: retail scanner data on the US beer market, covering before and after the joint venture, containing brand shares, prices, and regional variation.

Results: the actual post-merger price increase exceeds the prediction of a pure Bertrand simulation, and the data support the merger having brought about some coordinated pricing (the conduct parameter departs from static Bertrand), which is a warning to merger analysis that looks only at unilateral effects. The paper demonstrates that conduct can be tested rather than defaulted.

Limitations: identifying the change in conduct relies on the before-and-after comparison and suitable cost shifters, and extrapolating to settings without such a natural experiment calls for caution; the specific mechanism of coordinated pricing remains reduced-form. The value of this paper is in making the point that "conduct is an assumption that needs testing" into an operational empirical exercise.
:::

### 5.3 Conlon and Mortimer (2021)

::: {.case}
Paper: "Empirical Properties of Diversion Ratios," RAND Journal of Economics. It systematically makes clear the supply-side lifeline of the diversion ratio: how to define it, how to compute it from different demand models, how different algorithms differ, and its relationship to second-choice data.

Method: place the diversion ratio in a random-utility framework for unified treatment, give the aggregation formula for individual diversion, explain that different interventions like price changes, product removal, and quality changes correspond to different aggregation weights, and connect diversion to second-choice data and to UPP in merger analysis.

Data: the paper demonstrates with settings like snack foods (distributor data) and also relies heavily on analytical and simulation comparisons of random-coefficients models.

Results: the diversion given by different demand models varies widely, the IIA of plain logit systematically distorts it, and this distortion transmits directly to merger and UPP analysis. The paper therefore stresses that diversion should be directly reported and tested as a core product of demand estimation, rather than computed incidentally after the fact.

Limitations: precise diversion relies on the demand model being specified correctly, especially on random coefficients getting the substitution structure right; different interventions correspond to different aggregation weights, which must be stated clearly when reporting. This paper substantiates a sentence that recurs throughout this chapter: everything on the supply side, for better or worse, hangs on whether the previous chapter estimated diversion correctly.
:::

Taking the three together, the meaning of anchoring becomes clear: Nevo brings structural merger simulation into antitrust, Miller-Weinberg demonstrates that conduct can be tested rather than defaulted, and Conlon-Mortimer makes the supply-side lifeline of diversion clear. Progress in supply-side analysis has always turned on these three through-lines: correctly backing out cost from demand, treating conduct as an assumption to be tested, and getting diversion right.

## 6 A Full Walkthrough on the Selva Data

Now run all the tools of Section 4 fully on Selva. The following code uses R 4.5.3, fixing set.seed(17) for reproducibility, with demand taken directly from the previous chapter's BLP estimation, and every number cited in the text comes from the actual run output of this code.

### 6.1 Backing Out Marginal Cost

```r
# H0 = identity matrix (6 single-product firms pre-merger); Delta = negative matrix of demand price derivatives (BLP demand)
Delta_t <- function(sij, alpha, sigp){        # Delta_jk = -ds_j/dp_k, diagonal positive
  ai <- alpha - sigp*nu_p                      # individual price sensitivity
  ... # diagonal mean(ai*sij_j*(1-sij_j)); off-diagonal -mean(ai*sij_j*sij_k)
}
mc <- p - solve(H0 * Delta, s)                 # back out market by market
```

Using BLP demand and assuming single-product Bertrand, the median backed-out marginal cost is $5.94$ dollars, all positive, hugging the true value of $5.81$. The median markup is about $1.65$ dollars and the Lerner index is $0.23$. This set of costs is the raw material for all the counterfactuals that follow.

### 6.2 Conduct: Bertrand versus Collusion

```r
H_bert <- diag(J)                              # Bertrand: go it alone
H_coll <- matrix(1, J, J)                      # full collusion: joint monopoly
mc_bert <- p - solve(H_bert * Delta, s)
mc_coll <- p - solve(H_coll * Delta, s)
```

The median cost backed out under Bertrand is $5.94$, all positive. The median cost backed out under collusion is negative $7.89$, negative for $99.8\%$. This is the absurd number of Section 1: the wrong conduct assumption gives a physically impossible cost. Regressing the backed-out cost on cost shifters (local input cost $w$, quality $x_1$) and rivals' cost shifters and measuring fit with a GMM objective, Bertrand's objective is about four orders of magnitude smaller than collusion's (a ratio near $1.4 \times 10^4$), and the test cleanly selects Bertrand. Conduct is not defaulted, it is picked out by the data.

![Testing conduct: scatter the backed-out marginal cost against the cost shifter w. Bertrand (blue) backs out costs that rise steadily with w and are all positive, consistent with the true cost structure; full collusion (red) backs out costs that fall almost entirely below the zero line, physically impossible and falsified at a glance.](assets/fig/fig_21_conduct.svg)

### 6.3 Simulating the Merger: How the Demand Model Changes the Conclusion

```r
H1 <- diag(J); H1[1,2] <- 1; H1[2,1] <- 1      # brand 1 and 2's manufacturers merge
# Morrow-Skerlos zeta-map for post-merger equilibrium prices
solve_merger <- function(a, x1, H, mc, p0){
  p <- p0
  repeat {
    ... # Lam = diagonal price sensitivity; Gam = substitution terms
    eta_new <- (s + (H*Gam) %*% (p - mc)) / Lam
    p_new <- mc + eta_new
    if (max(abs(p_new - p)) < 1e-11) break; p <- 0.5*p + 0.5*p_new
  }
  p_new
}
```

Simulate this merger once each with true demand, BLP-estimated demand, and the previous chapter's IV-logit demand, and place the results side by side:

| Demand model | Brand 1 price rise | Brand 2 price rise | Rival price rise | Per-consumer $\Delta CS$ |
|---|---|---|---|---|
| True demand | +5.2% | +7.8% | +0.3% | -\$0.22 |
| BLP estimate | +4.8% | +7.3% | +0.3% | -\$0.21 |
| IV-logit | +3.8% | +5.8% | +0.2% | -\$0.16 |

This table is the confluence of the previous chapter and this one. BLP demand nearly reproduces the true merger effect, brand 2's price rise of $7.3\%$ against the true value of $7.8\%$, and the consumer surplus loss almost ties as well. IV-logit systematically underestimates: brand 2's price rise is only $5.8\%$, a quarter below the true value, and the consumer surplus loss is underestimated by about $27\%$. The reason is single and clear: logit's IIA underestimates the diversion between brand 1 and brand 2 as $0.206$ (true value $0.271$), so the mutual cannibalization internalized by the merger is underestimated, and the pressure to raise prices is underestimated along with it. The payoff of the previous chapter getting diversion right is only fully realized here: use the wrong demand model, and you misjudge as harmless a merger that ought to be scrutinized.

![Price effects of brand 1 and 2's manufacturers merging, under three demand models. The BLP estimate (blue) nearly reproduces the price increases of true demand (navy), while IV-logit (red) systematically underestimates the increases because it underestimates the neighboring diversion, and rival brands are almost unaffected.](assets/fig/fig_21_merger.svg)

### 6.4 GUPPI and the Efficiency Threshold

```r
guppi2 <- D_21 * (p1 - mc1) / p2               # brand 2's first-order upward pricing pressure
```

Brand 2's GUPPI is $0.074$ and brand 1's is $0.049$. Brand 2's GUPPI ($7.4\%$) is almost equal to the full simulation's price increase ($7.3\%$), confirming the accuracy of GUPPI here as a first-order approximation. As for the efficiency defense, solve an inverse problem: by how much must the merging parties' marginal cost fall to offset brand 2's price increase. The answer is $14.5\%$. So the defense that "the merger has efficiencies and will not harm consumers" turns into a checkable threshold: can this merger really deliver a $14.5\%$ cost reduction?

### 6.5 Full Reconciliation

The walkthrough of this section can be summarized as follows: the marginal cost backed out from the previous chapter's BLP demand (median $5.94$, hugging the true value of $5.81$) stands up under the correct Bertrand conduct and is absurd under the collusion assumption (negative $7.89$), on the basis of which simulating the brand 1 and 2 merger gives brand 2 a $7.3\%$ price rise and a per-consumer surplus loss of $0.21$ dollars, while using IV-logit demand underestimates the price rise to $5.8\%$ and underestimates the consumer loss by about $27\%$, GUPPI's $0.074$ precisely foreshadows the equilibrium price increase, and a $14.5\%$ cost efficiency exactly offsets it. All these numbers flow from the same demand estimation, and so the credibility of the supply side is entirely staked on the correctness of the demand estimation and the conduct assumption.

## 7 Failure Modes and Robustness

In the simulation both conduct and demand are constructed, but in real research they can fail at any moment. This section reviews the most common ways they fail and the operational responses.

Conduct misspecification is the most dangerous error on the supply side, because it is silent and consequential. The first-order condition always solves for cost, and a wrong conduct only shows its hand when the cost is absurd, while two conducts that both give positive costs (static Bertrand and mild coordination, say) cannot be told apart from cost alone yet give very different merger predictions. The response is to treat conduct as an assumption to be tested: report the plausibility of the backed-out cost (nonnegativity, the relationship with cost shifters), and run a formal Rivers-Vuong type test when a suitable instrument is available, rather than defaulting to Bertrand. Defaulting to an untested conduct is as groundless as defaulting to collusion in Section 1.

Demand misspecification amplified down the first-order condition is the second fatal link. Any error in the previous chapter that estimates diversion wrong is, on the supply side, amplified into wrong costs and wrong merger predictions. Estimating demand with logit rather than random coefficients, IIA underestimates neighboring diversion, so it underestimates multiproduct markups and underestimates the merger price increase, and the $5.8\%$ against $7.8\%$ of Section 6 is a live example. The supply side therefore cannot be severed from demand estimation, and the robustness checks of demand (instrument strength, random-coefficients identification, functional form) are at the same time robustness checks of the supply-side conclusions.

The validity of counterfactual extrapolation is the third threat, essentially the Lucas critique incarnate in merger analysis. Merger simulation assumes demand parameters, marginal costs, and the product portfolio are unchanged before and after the merger, with only the ownership structure changing. But real mergers often come with product repositioning, new entry, and adjustments to cost structure, all of which lie outside the model. If a merger would significantly change the market structure (large repositioning, inducing entry or exit), the static simulation's price prediction is unreliable and requires a dynamic model or one with entry. The simulation gives the short-run effect "with all else unchanged," and treating it as a long-run prediction calls for particular care.

There are a few more specific gates. A partial simulation solves only the merging parties' first-order conditions and fixes rivals' prices, saving computation over the full simulation but usually underestimating the price increase, and which one is used must be stated when reporting. The cost reduction in the efficiency defense must be a reduction in marginal cost to transmit to price; savings in fixed cost do not enter the pricing first-order condition and do not benefit consumers, and confusing the two is a common error. The pass-through rate determines how much of a cost change enters price, and it itself depends on the curvature of demand, with random coefficients and logit giving different pass-through, so the quantitative conclusion of the efficiency defense also depends on the demand specification. Finally, converting to dollars in the welfare measure via the log-sum requires the marginal utility of price, and if the model allows some consumers to have a positive price coefficient (the old problem of normal random coefficients), the welfare conversion breaks, which is one reason practice often sets the price coefficient as log-normal.

Stringing these failure modes together, the credibility of the supply side comes down in the end to two things: whether the conduct assumption is right, and whether the previous chapter's demand (especially diversion) is estimated accurately. The plausibility of the backed-out cost, the formal test of conduct, the robustness of demand, and the extrapolation limits of the counterfactual are all evidence provided around these two, and none of them can be replaced by "the first-order condition solved and the merger simulation ran."

## 8 Further Reading

::: {.readings}
Required reading, in suggested reading order:

- Nevo (2000, RAND Journal of Economics). The landmark text bringing structural merger simulation into antitrust; read it first to establish merger analysis that integrates demand and supply.
- Conlon and Mortimer (2021, RAND Journal of Economics). A systematic treatment of the diversion ratio, the supply-side lifeline, for understanding why getting substitution right is the premise of everything.
- Miller and Weinberg (2017, Econometrica). The modern exemplar of differentiated-products merger and conduct testing, demonstrating that conduct can be tested rather than defaulted.
- Farrell and Shapiro (2010, B.E. Journal of Theoretical Economics). The original formulation of UPP, the theoretical basis for first-order merger screening.
- Morrow and Skerlos (2011, Operations Research). The $\zeta$-map for numerically solving merger simulation, the computational method to read before doing merger simulation by hand.

Further reading:

- Duarte, Magnolfi, Sølvsten and Sullivan (2024, Quantitative Economics). The modern conduct-testing framework, a systematic treatment of Rivers-Vuong type tests in structural estimation.
- Bresnahan (1982, Economics Letters). The origin of conduct identification, how demand rotation distinguishes competition from collusion.
- Backus, Conlon and Sinkinson (2021, American Economic Journal: Microeconomics). The generalization of common ownership and the ownership matrix, how $\mathcal{H}$ encodes the internalization of interests at the investor level.
- Crawford, Lee, Whinston and Yurukoglu (2018, Econometrica). Nash-in-Nash bargaining, a framework for when the supply side is not simple pricing but upstream-downstream negotiation, an exemplar for platform-content negotiations.
- Small and Rosen (1981, Econometrica). The theoretical basis for discrete-choice welfare measurement, the bridge between the log-sum and consumer surplus.
:::

::: {.apa-refs}
- Backus, M., Conlon, C., & Sinkinson, M. (2021). Common ownership in America: 1980-2017. *American Economic Journal: Microeconomics, 13*(3), 273-308. https://doi.org/10.1257/mic.20190389
- Bresnahan, T. F. (1982). The oligopoly solution concept is identified. *Economics Letters, 10*(1-2), 87-92. https://doi.org/10.1016/0165-1765(82)90121-5
- Conlon, C., & Mortimer, J. H. (2021). Empirical properties of diversion ratios. *The RAND Journal of Economics, 52*(4), 693-726. https://doi.org/10.1111/1756-2171.12388
- Crawford, G. S., Lee, R. S., Whinston, M. D., & Yurukoglu, A. (2018). The welfare effects of vertical integration in multichannel television markets. *Econometrica, 86*(3), 891-954. https://doi.org/10.3982/ECTA14031
- Duarte, M., Magnolfi, L., Sølvsten, M., & Sullivan, C. (2024). Testing firm conduct. *Quantitative Economics, 15*(3), 571-606. https://doi.org/10.3982/QE2319
- Farrell, J., & Shapiro, C. (2010). Antitrust evaluation of horizontal mergers: An economic alternative to market definition. *The B.E. Journal of Theoretical Economics, 10*(1), Article 9. https://doi.org/10.2202/1935-1704.1563
- Miller, N. H., & Weinberg, M. C. (2017). Understanding the price effects of the MillerCoors joint venture. *Econometrica, 85*(6), 1763-1791. https://doi.org/10.3982/ECTA13333
- Morrow, W. R., & Skerlos, S. J. (2011). Fixed-point approaches to computing Bertrand-Nash equilibrium prices under mixed-logit demand. *Operations Research, 59*(2), 328-345. https://doi.org/10.1287/opre.1100.0894
- Nevo, A. (2000). Mergers with differentiated products: The case of the ready-to-eat cereal industry. *The RAND Journal of Economics, 31*(3), 395-421. https://doi.org/10.2307/2600994
- Rivers, D., & Vuong, Q. (2002). Model selection tests for nonlinear dynamic models. *The Econometrics Journal, 5*(1), 1-39. https://doi.org/10.1111/1368-423X.t01-1-00071
- Small, K. A., & Rosen, H. S. (1981). Applied welfare economics with discrete choice models. *Econometrica, 49*(1), 105-130. https://doi.org/10.2307/1911129
:::
