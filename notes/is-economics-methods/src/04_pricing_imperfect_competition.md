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

The same demand curve and the same marginal cost can give birth to wildly different profits. A firm that owns the market will push its price up to wherever demand elasticity allows; two firms competing in quantities earn less, but still earn something positive; if what they sell is perfectly homogeneous and they compete directly in prices, profit can slide all the way down to zero. The price schedule has not become any more "reasonable," only the competitive rules behind it have changed. So when you see a high price, the real question is not whether the firm is greedy, but what demand, what rivals, and what constraints it faces.

Market power turns price from a given number into something the firm chooses. Monopoly pricing needs consumer theory to supply demand and welfare, oligopoly needs game theory to describe how firms react to one another, and second-degree price discrimination writes a screening mechanism into a menu of plans. This chapter brings the first three chapters together in the Cirrus market: first we compute a single price, group pricing, and nonlinear tariffs, then we compare Cournot, Bertrand, differentiated competition, vertical relations, and dynamic pricing. Each model does more than hand you a price; it also delivers testable comparative statics: how a cost shock is passed through, how markups change as products become closer substitutes, and why a capacity constraint can soften a price war.

The point of a pricing model is not to find some magic price for the firm, but to explain where profit comes from, how much competition takes away, and what consumers pay for it. Later demand estimation will recover price elasticities, supply models will back out marginal costs, and merger analysis will change a firm's first-order condition; whether any of that work holds up depends on whether the competitive model chosen here looks like the market itself. Pick the wrong game, and even a precisely solved equilibrium is just the answer to some other market.

## 1 A Profit Caught in the Middle

::: {.case}
Cirrus is the dominant e-commerce platform in a particular niche. To shrink the problem down to something you can do by hand, we abstract its core business into selling "operating services" to merchants: the more a merchant uses, the more the platform charges. We summarize the market for this service with a single linear demand, $p = 12 - q$, where $p$ is the price per unit of operating activity (think of it as the price equivalent of a commission) and $q$ is total merchant usage, and the platform's marginal cost of providing the service is $c = 2$.
:::

An analyst runs the textbook monopoly-pricing calculation. As the monopolist for this service, the platform sets marginal revenue equal to marginal cost, $12 - 2q = 2$, giving the optimal quantity $q = 5$, the price $p = 7$, and profit $(7 - 2)\times 5 = 25$. On this basis the analyst writes down a conclusion: Cirrus uses its market power to extract profit of 25, merchants are left with consumer surplus of 12.5, and market power creates a deadweight loss of 12.5, a loss that is the unavoidable price of monopoly.

Every step of that 25 is computed correctly, yet it is neither a ceiling nor a floor; it is just a number sitting in the middle of a very wide range. Start with why it is not a ceiling. Part of that "unavoidable" deadweight loss of 12.5 is actually a choice of pricing instrument, not an iron law of market power. If Cirrus does not charge every merchant a single price, but instead lays out a tiered commission menu and lets merchants of different sizes self-select, it can earn more; Section 6 will show that this kind of second-degree price discrimination lifts profit from 25 to 26. And if it could price perfectly and individually to each merchant, it could pocket both the 12.5 in consumer surplus and the 12.5 that was previously deadweight loss (total profit of 50). A single price is far from the limit of extraction. Now for why it is not a floor. That 25 rests on the premise that Cirrus owns the market. The moment a rival platform, Nimbus, enters and competes, the price gets pushed down: if the two compete in quantities (Cournot), the price falls to about 5.33 and total industry profit shrinks to about 22.2; if the two offer an undifferentiated service and compete in prices (Bertrand), the price falls all the way to marginal cost 2 and profit goes to zero. In the same market, profit can slide from near-complete extraction all the way down to the competitive level of zero.

So the question "how much does the platform earn" has no single number as its answer, but a range: the upper bound is set by the pricing instrument (can it discriminate) and the lower bound by the competitive structure (are there rivals, and how do they compete). That isolated 25 is only the outcome under one specific case, "single price plus exclusive monopoly," and treating it as the last word on platform profit both understates the ceiling that discrimination can extract and ignores the floor that competition can push it toward. What this chapter does is mark out both bounds of that range and the ground in between: Section 3 covers monopoly pricing and price discrimination and marks the upper bound; Section 4 covers imperfect competition and marks the lower bound; Section 6 works every one of these numbers out on the Cirrus market.

## 2 Market Power and How to Measure It

Before rolling up our sleeves, let us be clear about what market power is and how to measure it. The price taker of Chapter 1 faces horizontal demand: ask for one cent more and you sell nothing, so you can only price at marginal cost. A firm with market power faces downward-sloping demand: it knows that raising its price loses some customers but not all of them, so it can trade off "earning more per unit" against "losing some customers," and set its price above marginal cost. The essence of market power is facing a downward-sloping (rather than horizontal) demand curve.

The monopolist's problem is therefore to choose quantity (or equivalently price) to maximize profit $\pi(q) = p(q)\, q - c(q)$, where $p(q)$ is inverse demand. The first-order condition is that marginal revenue equals marginal cost:

$$\underbrace{p(q) + q\,\frac{dp}{dq}}_{MR} = \underbrace{c'(q)}_{MC}$$

Marginal revenue is below price by an amount $q\,\frac{dp}{dq}$, and this piece is the loss from "to sell one more unit you must cut the price, and the cut drags down every unit already sold"; it is exactly what pushes monopoly quantity below the competitive level. Rearranging the first-order condition gives the most famous measure of market power, the Lerner index:

$$\frac{p - c'}{p} = \frac{1}{|\varepsilon|}$$

where $\varepsilon$ is the price elasticity of demand. The left side $(p - c')/p$ is the markup, the fraction by which price exceeds marginal cost; the right side says it equals exactly the inverse of demand elasticity. This formula quantifies market power: the more inelastic demand is ($|\varepsilon|$ smaller, customers more dependent on you), the higher the markup and the greater the market power; as elasticity goes to infinity (perfect competition) the markup goes to zero and we are back to price equal to marginal cost. The Lerner index is a quantity this series will keep recovering in later structural analysis: demand estimation gives the elasticity, and the supply-side markup is pinned down by it. Section 6 will verify this relationship on the Cirrus numbers.

The quantities we want to measure in this chapter follow the welfare language of Chapter 1: consumer surplus (the area under the demand curve and above the price line), producer surplus or profit, and the gap between their sum and the social optimum, that is, the deadweight loss (DWL). The social cost of market power is exactly this DWL, which comes from output being held below the level set by marginal-cost pricing: some trades that should have happened, where the buyer's valuation exceeds marginal cost, do not happen. The welfare effects of price discrimination, competition, and collusion all come down in the end to whether this area shrinks or grows.

To summarize: market power is facing downward-sloping demand; the monopolist prices by setting marginal revenue equal to marginal cost, pushing price above marginal cost, and the markup is characterized by the Lerner index equal to the inverse elasticity; the social cost of market power is the deadweight loss from underprovision. The next section starts from the single monopoly price and looks at how a firm uses price discrimination to extract more, and what that means for welfare.

## 3 Core Structure: Monopoly Pricing and Price Discrimination

This section is the most detailed part of the chapter's analysis. First we lay out single-price monopoly fully, then along the line of "how much information about buyers can be obtained, and how far arbitrage can be blocked" we unfold the three degrees of price discrimination one at a time, and finally we touch on a subtle constraint on monopoly: the Coase conjecture.

### 3.1 Single-Price Monopoly and Deadweight Loss

Start with the simplest case: the monopolist charges every buyer the same price. Use this chapter's linear demand $p = 12 - q$ and marginal cost $c = 2$ to set up the geometry. Inverse demand is the straight line from $(0, 12)$ to $(12, 0)$; the marginal revenue curve has twice its slope, starts from the same intercept, and falls more steeply, $MR = 12 - 2q$. Setting $MR = MC$, that is $12 - 2q = 2$, gives monopoly quantity $q_m = 5$, and substituting back into demand gives monopoly price $p_m = 7$ and profit $(7 - 2)\times 5 = 25$.

![Single-price monopoly: demand $p = 12 - q$, marginal revenue $MR$, marginal cost $MC = 2$. $MR = MC$ pins down $q_m = 5$, $p_m = 7$. Above the price line is consumer surplus (12.5); between price and MC up to $q_m$ is monopoly profit (25); the triangle from $q_m$ to the competitive quantity $q_c = 10$ is deadweight loss (12.5).](assets/fig/fig_04_monopoly.svg)

Compare this outcome with perfect competition and the cost of market power is plain to see. A competitive market would push price down to marginal cost $p = c = 2$ and raise quantity to $q_c = 10$. Monopoly produces 5 units fewer than competition, and every buyer of those 5 units values the good above marginal cost 2 (they lie on the segment of the demand curve from $q_m$ to $q_c$), so trades that should have happened did not, and the lost social surplus is that triangle, with area $\frac{1}{2}\times(10 - 5)\times(7 - 2) = 12.5$, which is the deadweight loss. Consumer surplus here is the triangle above the price line and below the demand curve, $\frac{1}{2}\times(12 - 7)\times 5 = 12.5$, which also happens to be 12.5. These are exactly the numbers the analyst reported. They are not wrong, but they assume the monopolist will only ever charge a single price, and that is precisely the self-imposed limit the monopolist is least inclined to accept. Below we see how it breaks through.

### 3.2 Three Kinds of Price Discrimination

Price discrimination means charging different prices to different buyers, or for different units bought by the same buyer, with the goal of extracting the surplus that a single price leaves in consumers' hands. Whether it can be implemented depends on two conditions: the firm must be able (directly or indirectly) to tell buyers' willingness to pay apart, and it must be able to stop low-price buyers from reselling to high-price buyers (no arbitrage). By how much information the firm has about buyers, price discrimination is traditionally divided into three degrees.

First-degree price discrimination (also called perfect price discrimination) is the ideal limit: the firm knows exactly each buyer's willingness to pay for each unit, and charges each unit the most that buyer is willing to pay. The result is that the firm takes for itself all the social surplus under the demand curve and above marginal cost, and consumer surplus goes to zero. Its efficiency implication is interesting: because the firm charges the last unit exactly the buyer's valuation, it is willing to sell as long as that valuation exceeds marginal cost, so output reaches the competitive level $q_c = 10$ and deadweight loss is zero. First-degree discrimination is the worst for consumers in terms of distribution (all their surplus is taken away), yet the best in terms of efficiency (output is largest, DWL is zero), and this contrast, that "the most thorough extraction is precisely the most efficient," is the most counterintuitive point in the welfare analysis of price discrimination. In reality perfect information is unattainable, so first-degree discrimination is an unreachable upper bound, but it marks the ceiling on extraction.

Third-degree price discrimination is the most common kind: the firm cannot price to each person, but it can split buyers into several segments by some observable characteristic and charge each segment a different single price, for instance pricing by region, student status, or whether the account is a business or an individual. The firm runs one monopoly-pricing problem within each segment, setting that segment's marginal revenue equal to the common marginal cost. The conclusion is a clean rule: the more inelastic the segment, the higher the price it is charged. Because each segment's markup is set by the Lerner index equal to that segment's inverse elasticity, the inelastic segment, the one that cannot do without the service, is charged more. Section 6 will work out this inverse elasticity rule using two segments.

Second-degree price discrimination is the most delicate kind, and the closest to platform business. It targets the case where the firm cannot even observe which segment a buyer belongs to, and only knows that several types of buyers exist in the market and their proportions. The firm therefore does not price by buyer identity, but designs a self-selecting menu, usually a set of quantity-and-total-price combinations (equivalent to nonlinear pricing or tiered plans), and lets buyers choose for themselves, so that their choices indirectly reveal their types. This is exactly the pricing version of the screening problem of Chapter 3: the menu must satisfy incentive compatibility (each type voluntarily chooses the tier designed for it) and individual rationality (voluntary purchase), and solving it forces out the same information rents and downward distortion at the bottom. The conclusion is also consistent with Chapter 3: the quantity given to the high type is undistorted (efficient at the top), while the quantity given to the low type is deliberately held below the efficient level in order to compress the information rent that must be conceded to the high type, and the high type enjoys a slice of information rent. This "distort the bottom, leave the rent at the top" structure is exactly the economic reason a platform deliberately limits its basic plan and splits its commission into several tiers, and Section 6 works it out fully on Cirrus's two types of merchants, which is also the core of this chapter's running case.

::: {.intuition}
The distinction between the three kinds of price discrimination falls into place once you grasp "how much the firm knows about buyers, and how it prevents arbitrage." First-degree is omniscient: it prices to each person and each unit, taking all the surplus while producing the most. Third-degree knows which type you belong to (a visible label) and prices by type, charging the inelastic type more. Second-degree cannot even tell your type, so it lays out a menu and lets you self-select, trading a crippled bottom for rent at the top. The less information, the less complete the extraction and the larger the distortion, which is why second-degree discrimination has the most subtle efficiency loss of the three. Platforms do both third-degree (pricing by region and account type) and second-degree (tiered subscriptions and stepped commissions), and once you see which kind of information it currently holds about you, you know which degree of discrimination it is using.
:::

### 3.3 The Coase Conjecture: The Durable-Goods Monopolist Trapped by Itself

Monopoly pricing has one more counterintuitive weak spot, which appears when selling a durable good and being unable to commit to future prices, and this is the Coase conjecture. A durable good, once bought, lasts a long time, so a monopolist selling a durable good is essentially competing with its own future self: today it sells to high-valuation buyers at a high monopoly price, and tomorrow, facing the remaining low-valuation buyers, it has an incentive to cut the price and sell another round. But rational high-valuation buyers foresee that tomorrow the price will fall, so they would rather wait. This expectation of holding out for a lower price in turn keeps the seller from daring to price too high today.

Coase's conjecture is this: when the seller can reprice very frequently (the interval between price changes goes to zero) and cannot commit not to cut prices, this competition with its future self becomes so fierce that it drives the price all the way toward marginal cost, and the monopolist loses almost all of its market power. In other words, a durable-goods monopolist without commitment power gets trapped by its own future incentive to cut prices and prices near the competitive level. The practical significance of this result is that it points out the value of commitment: if the monopolist can credibly commit not to cut prices (by leasing rather than selling, by most-favored-customer clauses, by planned scarcity), it can escape the Coase trap and preserve its market power. The recurring leasing, subscription, and time-limited pricing that platforms use when selling digital durables (software licenses, memberships) all reflect a concern with fighting the Coase conjecture and maintaining pricing power. This also confirms the weight of the commitment theme from Chapter 2 in pricing: whether you can turn "I will not cut prices" into a credible commitment directly determines how much market power a durable-goods monopolist can keep.

### 3.4 Nonlinear Pricing: Two-Part Tariffs and Bundling

The second-degree discrimination of Section 3.2 uses a menu of "quantity versus total price" to sort buyers, and in platform and SaaS pricing this idea lands as two of the most common families of tools, both belonging to nonlinear pricing (total payment is not proportional to quantity bought), worth unpacking separately.

The first family is the two-part tariff (Oi 1971): total charge $T(q) = F + p\, q$, a usage-independent access fee (or membership fee) $F$ plus a per-unit price $p$. Its beauty is in the optimum when buyers are homogeneous: set the per-unit price $p$ at marginal cost $c$, let buyers buy the efficient quantity (creating no usage distortion), then use the fixed fee $F$ to sweep up all the resulting consumer surplus in one go. The result is that the firm captures the entire surplus that should have gone to buyers, while output reaches the competitive level: with two pricing instruments it reproduces the effect of first-degree discrimination. Section 6.5 will substitute this chapter's demand and find $F = 50$ and profit 50, exactly twice the single-price 25, with output rising to the efficient 10. When buyers are heterogeneous there is a trade-off: a high $F$ extracts more from large buyers but shuts small buyers out, so the optimal two-part tariff must balance "the markup per unit" against "how many buyers to cover," which is exactly the extension of the discrete menu of Section 3.2 to continuous types. Membership fees plus usage fees, a platform's listing fee plus transaction commission, subscription floors plus overage charges, the shadow of the two-part tariff is everywhere.

The second family is bundling (Adams and Yellen 1976): selling several goods together as a package. Counterintuitively, even when marginal cost is zero and each good priced separately is already optimal, as long as buyers' valuations of the different goods are negatively correlated, pure bundling can still earn more. The reason is that bundling smooths out the differences in willingness to pay across buyers, letting the firm set the bundle price close to each buyer's total valuation and thereby capture the surplus that separate pricing misses. In the two-good two-buyer example of Section 6.5, selling separately yields 20 in total, while pure bundling grabs 24 at a stroke. Going one step further, mixed bundling (selling the bundle and the individual goods at the same time) is usually even better. This explains why software suites, content memberships, and tiered platform services always sell a pile of features bundled together: bundling extracts exactly the surplus that piece-by-piece pricing leaves on the table.

### 3.5 Versioning: Designing Quality Menus for Information Goods

The high fixed cost, near-zero reproduction cost, and designable quality of information goods make versioning a central form of second-degree price discrimination (Shapiro and Varian 1999). The firm first builds a high-quality master, then forms versions through features, speed, timeliness, API limits, or service levels, letting consumers with different willingness to pay self-select. The limitations of the low-end version do not necessarily come from technical cost, but may be a degradation added deliberately to maintain incentive compatibility.

Let types be $\theta_H>\theta_L$, quality be $s$, and utility be $\theta s-p$. The menu $(s_L,p_L),(s_H,p_H)$ must satisfy the two IC and IR constraints. The standard single-crossing result is still efficient at the top: $s_H$ is undistorted, $s_L$ is distorted downward, in order to reduce the information rent that must be left when the high type imitates the low-end version. So "why the free version is deliberately missing features" is not an accident of product design, but a testable prediction of the screening constraint.

Versioning also has boundaries. Too many versions increase menu complexity and maintenance cost; competition compresses the room for artificial degradation, because a rival can use a higher-quality low-end product to fight for the distorted users. Empirically, sharper external competition should predict a narrowing of the quality gap between versions and an improvement of the low-end version, which we can call unversioning pressure. A researcher cannot conclude that price discrimination works simply by observing that premium users pay more; they also need to distinguish quality differences, selection into products, and genuine willingness-to-pay heterogeneity.

### 3.6 Intertemporal and Dynamic Pricing: Skimming, Penetration, and Hotelling's Rule

The Coase conjecture of Section 3.3 is the durable-goods pricing dilemma under no commitment; widening the view to the intertemporal setting, dynamic pricing has a few more pictures worth spelling out.

First clear up a common misconception: given that a durable-goods monopolist can commit, should it price skim (set a high price for high-valuation buyers first, then gradually cut the price to harvest low-valuation buyers)? Stokey's (1979) answer is no. Facing forward-looking rational buyers, a durable-goods monopolist that can commit optimally commits to a constant price, and intertemporal price discrimination earns nothing extra, because any anticipated future price cut only induces high-valuation buyers to hold their cash and wait, cannibalizing the current period. Skimming is genuinely profitable only under other conditions: buyers are myopic, or the good is perishable and a fresh batch of buyers who cannot wait arrives each period (airline seats and concert tickets are exactly like this), or cost or quality falls over time. The opposite pattern to skimming is penetration pricing (low price first, higher price later), which is optimal only under network effects or learning-by-doing: use a low introductory price to build the installed base first, then raise the price once network value takes hold. A platform's freemium and low-price customer acquisition are essentially penetration pricing, aimed at breaking past that adoption threshold of Chapter 2, and this thread leads to the platform chapter later on.

The other pole of dynamic pricing is a fixed and exhaustible stock, and the classic result is Hotelling's Rule (Hotelling 1931). Consider an exhaustible resource with marginal extraction cost $c$ and interest rate $r$. The resource owner can at any moment choose one of two options: extract one unit now, take the net price (the scarcity rent) $p_0 - c$ and invest it at rate $r$, or leave it in the ground and sell it in period $t$ for $p_t - c$. No arbitrage requires the two paths to be equivalent, so the scarcity rent must grow at exactly the interest rate $r$:

$$p_t - c = (p_0 - c)(1 + r)^t$$

The intuition is a one-line arbitrage: if the rent grows more slowly than $r$, you might as well extract everything now and invest it; if it grows faster than $r$, you might as well keep it in the ground. In equilibrium the scarcity rent appreciates at $r$, its present value is constant, and Section 6.6 will trace out this rent path. The use of this rule goes far beyond mining: any fixed and scarce stock, spectrum licenses, reserved server capacity, the option value of a batch of one-time data, the quota cap on API calls, should be released over time so that its shadow rent appreciates at $r$. Dumping it all at once today gives up the appreciation; hoarding it forever unused gives up the interest. The dynamic pricing of scarce capacity, such as spot pricing and congestion pricing in cloud computing, echoes exactly Hotelling's logic: price scarcity itself, and let its rent chase the cost of capital.

To summarize: single-price monopoly holds output below the competitive level and creates deadweight loss; price discrimination is divided into three degrees by how much information the firm has, with first-degree extracting perfectly yet most efficient, third-degree sorting by observable label and charging the inelastic more, and second-degree using a self-selecting menu to sort indirectly, trading distortion at the bottom for information rent at the top, exactly the pricing incarnation of the screening of Chapter 3; nonlinear pricing lands second-degree discrimination as two practical families, the two-part tariff (per-unit price pushed to marginal cost, fixed fee sweeping up the surplus) and bundling (smoothing valuation differences to extract one more layer); intertemporally, a durable-goods monopolist that can commit optimally holds a constant price (Stokey), skimming requires myopia or perishability, penetration relies on network effects, while one that cannot commit is pushed by the Coase conjecture toward competitive pricing, and the scarcity rent of an exhaustible stock follows Hotelling's Rule, growing at the interest rate. The next section moves from a single firm to several firms and looks at how competition pushes prices down.

## 4 Imperfect Competition: The Several Faces of Oligopoly

Monopoly is one firm dominating, perfect competition is countless small firms, and real platform markets often fall in between: a few firms with market power competing with one another, which is oligopoly. The outcome of an oligopoly depends heavily on what dimension the firms compete in and how many times the competition is repeated, and this section lays out several main forms of competition one at a time; they are all direct applications of the game theory of Chapter 2.

### 4.1 Cournot: Competing in Quantities

In Cournot competition, firms choose quantities simultaneously and the market price is determined by total quantity through demand. Use this chapter's demand $p = 12 - Q$ ($Q = q_1 + q_2$), marginal cost $c = 2$, and two firms. Firm 1's profit is $(12 - q_1 - q_2 - 2)\,q_1$; differentiating with respect to $q_1$ and setting it to zero gives the first-order condition $12 - 2q_1 - q_2 - 2 = 0$, that is, the best response $q_1 = (10 - q_2)/2$. This is firm 1's best-response function, which falls as the rival's quantity rises: if the rival produces more, the market price is low, so I should produce less. Firm 2 is symmetric. The intersection of the two best-response curves is the Cournot-Nash equilibrium.

![Cournot best-response functions: two downward-sloping response curves $q_1 = (10 - q_2)/2$ and $q_2 = (10 - q_1)/2$ meet at the symmetric equilibrium $q_1 = q_2 = 10/3$. For comparison: collusion produces $5/2$ each (half the joint monopoly quantity), competition produces 5 each.](assets/fig/fig_04_cournot.svg)

Solving the symmetric equilibrium, set $q_1 = q_2 = q$, substitute to get $q = (10 - q)/2$, and solve to get $q = 10/3$, total quantity $Q = 20/3$, price $p = 12 - 20/3 = 16/3 \approx 5.33$, profit per firm $(16/3 - 2)\times(10/3) = 100/9 \approx 11.11$, and the two together $200/9 \approx 22.2$. Comparing this with monopoly is interesting: monopoly price 7, profit 25; two-firm Cournot competition pushes the price down to 5.33 and total profit down to 22.2. Competition did lower the price, raise output, and cut industry profit, but it did not push all the way to the bottom, and two Cournot firms each still keep considerable market power, because choosing quantities means each firm is wary that increasing its own output will depress the price and hurt itself. Quantity competition is a "mild" form of competition.

### 4.2 Bertrand: Competing in Prices and That Paradox

In Bertrand competition, firms choose prices rather than quantities simultaneously, and buyers all buy from whoever quotes the lowest price. The conclusion is startling: as long as the products are undifferentiated and marginal costs are equal, competition between two firms pushes the price all the way down to marginal cost $p = c = 2$ and profit goes to zero. The reason is an inescapable chain: as long as the rival prices above marginal cost, I can undercut slightly, grab the whole market, and still profit; the rival thinks the same, so both keep cutting until neither dares go below marginal cost and the price stops at $c$. This is the Bertrand paradox: price competition between just two firms reproduces the outcome of perfect competition and market power vanishes entirely.

Bertrand and Cournot give such disparate outcomes (price 2 versus 5.33) that it is often puzzling how the same duopoly can have two answers. The key is the dimension of competition: when competing in prices, grabbing customers is too easy (cut a penny and you take the whole market), which forces out a brutal zero profit; when competing in quantities, increasing output hurts your own price, so competition is mild. Which end reality falls at depends on whether firms actually behave more like they are choosing prices or choosing quantities. The three exits from the Bertrand paradox matter most: capacity constraints, product differentiation, and repeated play, which are the subjects of the next three subsections, and each can lift the price back up from marginal cost.

### 4.3 Kreps-Scheinkman: Capacity First, Prices After

The first exit from the Bertrand paradox is capacity. In reality firms often invest in capacity first and then compete in prices, and when capacity is limited the logic that "quoting the lowest price takes the whole market" breaks down: even at the lowest quote, capacity cannot absorb the whole market, so the rival can still sell to residual demand and hold a price above marginal cost. Kreps and Scheinkman (1983) turned this intuition into a beautiful two-stage result: firms choose capacity in the first stage and play Bertrand price competition under the capacity constraint in the second stage, and under a reasonable demand-rationing rule, the subgame-perfect equilibrium of this two-stage game reproduces exactly the Cournot quantities and price.

This result reconciles the apparent contradiction between Cournot and Bertrand: when competition contains the two stages of "build capacity first, then price," seemingly brutal price competition, by way of the prior commitment that is the capacity constraint, ends up at the mild Cournot outcome. It also again underscores the weight of commitment: the first-stage capacity investment is an irreversible commitment, and it is exactly what pulls the second-stage price competition back from the abyss of zero profit. So Cournot is not necessarily the product of the somewhat unrealistic assumption of "choosing quantities"; it can be the equilibrium outcome of the more realistic process of "invest in capacity first, compete in prices after."

### 4.4 Product Differentiation: Hotelling and the Softening of Competition

The second exit from the Bertrand paradox is product differentiation. The premise of the paradox is that products are completely undifferentiated and buyers look only at price, but once the two firms' products differ in buyers' eyes (brand, interface, ecosystem, location), the logic of cutting a penny to grab all demand fails, because there are always buyers who, out of preference, are willing to pay a bit more for you, and price competition softens accordingly. Hotelling's linear city model is the classic way to capture this.

Let buyers be uniformly distributed on a "street" of length 1 (this street can be real geography, or a product characteristic space, such as the spectrum of interface styles), with two firms at the two ends 0 and 1, and buyers, besides paying the price, also bear a mismatch cost that rises with "distance" (represented by the quadratic form $t\,d^2$, where $d$ is the distance from a buyer's ideal point to the product bought, and $t$ measures how sensitive buyers are to differences). A buyer's choice between the two firms depends on the trade-off between the price difference and the mismatch cost, which pins down an exactly indifferent boundary buyer and splits the market in two. Solving the symmetric price-competition equilibrium gives $p = c + t$: the equilibrium price equals marginal cost plus the degree of differentiation $t$. The larger the differentiation $t$, the further the price is above marginal cost and the more market power the firms have; as $t \to 0$ (products converge) the price falls back to marginal cost and the Bertrand paradox reappears. Section 6 substitutes $c = 2$ and $t = 1$ and computes the equilibrium price 3 and profit 0.5 per firm, providing a clean numerical contrast between undifferentiated Bertrand's zero profit and how differentiation rescues market power. Product differentiation is therefore an extremely crucial link in platform competition: creating differentiation is winning yourself a slice of market power that a price war cannot erase.

### 4.5 Collusion and Entry Deterrence

The third exit from the Bertrand paradox is repetition. The competition so far has been one-shot, but platforms coexist over the long run, and Chapter 2 showed that with infinite repetition and sufficiently patient firms, cooperation can be sustained by the credible threat that "deviation triggers retaliation." Carried over to pricing, this means several firms that should be competing can engage in tacit collusion: everyone holds a high price, and whoever cuts the price to grab customers triggers everyone dropping their prices in retaliation; as long as the discount factor exceeds the threshold given in Chapter 2, the collusive high price is a subgame-perfect equilibrium. Collusion is therefore not a product of market structure but of repeated interaction: the same set of firms, meeting once, falls into Bertrand zero profit, but coexisting over the long run may sustain a high price close to monopoly. This also poses a hard problem for empirical work: a high observed price could be collusion, or differentiation, or Cournot under capacity constraints, and looking at the price alone cannot tell them apart; additional structure and evidence are needed, which is the theme of the conduct tests later in this series.

Adjacent to competition is entry. The high profit of monopoly attracts entrants, so the incumbent has an incentive to engage in entry deterrence, for instance by expanding capacity in advance, occupying product space, or lowering price to signal that "there is no profit in entering." The key here is still commitment: a purely verbal threat that "if you dare enter I will start a price war" is not credible (the empty threat of Chapter 2), but if the incumbent sinks an irreversible cost in advance (excess capacity, R&D, locking in users), it turns the threat into a credible deterrent. Entry deterrence is therefore another battlefield for strategic commitment, and this series will take it up directly when it covers entry and market structure later on; here we need only know that market power depends not only on how many firms there are now, but also on how strongly potential entrants are deterred.

### 4.6 Vertical Relations: Double Marginalization, Vertical Integration, and RPM

The firms above are all horizontal rivals, but platforms are often embedded in a vertical chain: the platform buys from suppliers, or chooses between the marketplace mode (matching and taking a commission, with the seller setting the price) and the reseller mode (buying and reselling on its own account). The most important phenomenon in vertical relations is double marginalization (Spengler 1950). Suppose one upstream monopolist sells to one downstream monopolist retailer at a per-unit wholesale price $w$, and the retailer marks it up again and sells to consumers. Two monopoly markups stack on top of each other, and the result is that the final retail price is even higher than what a vertically integrated monopolist would set, while the combined upstream-plus-downstream profit is actually lower. In a word, two monopolists on one chain are worse than one: consumers pay a higher price and the firms earn less. Section 6.7 substitutes this chapter's demand and computes: when integrated, retail price 7 and profit 25; when upstream and downstream are separate, the wholesale price is pushed up to 7, the retail price shoots up to 9.5, and the two together earn only 18.75.

The root of the disease is that when the downstream firm marks up, it cares only about its own profit and ignores that this drags down the upstream firm's sales; this vertical externality has nothing to do with horizontal competition and is purely a product of the chain structure. Several remedies can cure it and pull the outcome back to the integrated monopoly optimum: direct vertical integration (the two merge and internalize the externality); a two-part contract between upstream and downstream (upstream pushes the wholesale price down to marginal cost $w = c$, eliminating the downstream markup, then uses a fixed franchise fee to sweep up the downstream profit, exactly the vertical version of the two-part tariff of Section 3.4); or resale price maintenance (RPM, with upstream directly stipulating the retail price equal to the integrated monopoly price 7). These moves look different in method, but all serve to remove for the downstream firm that redundant markup. This structure fits information systems especially well: whether a platform goes marketplace (the commission is roughly upstream, the seller's pricing roughly downstream) or reseller (integrated pricing on its own account) is essentially a choice of vertical structure, and the wholesale-versus-agency fight in e-books and app stores is asking exactly who resolves double marginalization. The Nash-in-Nash bargaining of the supply side later in this series is about how upstream and downstream haggle over that $w$.

### 4.7 A Taxonomy of Entry Deterrence

Section 4.5 noted that entry deterrence relies on credible commitment, and here we lay out the taxonomy from the Tirole line of work, because there is no one-size-fits-all answer to "what posture the incumbent should strike." There are four families of instruments. The first is capacity overinvestment (Dixit 1980): the incumbent sinks excess capacity in advance as a credible commitment that "if you enter I will fight hard," making entry unprofitable, with credibility premised on the capacity being irreversible and visible to the entrant. The second is limit pricing (Milgrom and Roberts 1982): keeping the price low before entry, so that under asymmetric cost information the low price credibly conveys the signal "my cost is low, you will not make money by entering," and it is exactly this signaling version that revives an otherwise hollow threat. The third is predatory pricing: pricing below cost to drive the entrant out, with credibility coming from cross-market reputation or the financing advantage of deep pockets. The fourth is product proliferation and preemption: filling up the product space so no profitable gap is left for the entrant.

The truly delicate part is the mnemonic taxonomy of Fudenberg and Tirole (1984): whether the incumbent should over-invest or under-invest depends on whether the post-entry actions are strategic substitutes (like quantities, the Cournot type) or strategic complements (like prices, the Bertrand type). Under strategic substitutes it is best to strike the top dog posture, growing yourself big and fierce (overinvestment) to intimidate; under strategic complements it is instead best to strike the puppy dog posture, making yourself small and meek (underinvestment) to avoid provoking a brutal price war, and the other two cells are lean and hungry and fat cat. The most counterintuitive lesson of this taxonomy is that deterring entry does not always mean "making yourself big and fierce"; sometimes showing weakness is optimal. The through-line across these four families is still the commitment of Chapter 2: only a sunk and observable commitment makes deterrence credible. So market power depends not only on how many firms there are right now, but also on how tightly potential entrants are shut out, and the full treatment is left to the entry chapter later in this series.

To summarize: the outcome of an oligopoly depends on the dimension of competition and the number of repetitions. Cournot quantity competition is mild and leaves market power; Bertrand price competition is brutal and reproduces the zero profit of perfect competition; but the three exits of capacity constraints (Kreps-Scheinkman pulls Bertrand back to Cournot), product differentiation (Hotelling's $p = c + t$), and repeated interaction (tacit collusion) can all lift the price back up from marginal cost. Vertical relations bring another problem, where double marginalization on the chain pushes the price even higher than integrated monopoly, resolved only by integration, a two-part tariff, or RPM; while entry deterrence relies on credible commitment to keep potential competition outside the door, and whether to act fierce or show weakness depends on whether competition is in strategic substitutes or strategic complements. The next section grounds these mechanisms in real research.

### 4.8 Pass-through, Markups, and Empirical Predictions

The most useful empirical restrictions of price theory are often not directional judgments like "the price will rise," but how large the derivative of price with respect to cost, demand, and rival shocks is. Take linear demand $Q=a-bp$: the single-product monopoly price is $p=(a+bc)/(2b)$, so the pass-through of marginal cost is

$$\frac{\partial p}{\partial c}=\frac{1}{2}.$$

Switching to symmetric $n$-firm Cournot competition with inverse demand $P=a-bQ$, the price becomes $p=(a+nc)/(n+1)$, and pass-through is $n/(n+1)$, approaching 1 as the number of firms grows. This difference gives a prediction consistent with conduct: the same exogenous input-cost shock should pass through to price by different amounts under different competitive structures. But beware a confound: demand curvature also affects pass-through, so seeing a pass-through of $0.5$ alone does not let you conclude the market is a monopoly; you also have to pin down the shape of demand.

In differentiated Bertrand competition, a rise in one product's cost not only raises its own price but also, through strategic complementarity (prices are strategic complements), pushes the rival's price up along with it. The sign and magnitude of this cross-price response depend on the demand substitution matrix and the ownership structure of the products. If a rival's price is observed to be perfectly still, it could genuinely be strategic independence, or the cost shock may simply not have been cleanly isolated, or price adjustment itself may be sticky, and these must be told apart.

The predictions of price discrimination must also put price, sales, and composition together. Once third-degree price discrimination is allowed, the prices of the segments each move by the inverse elasticity rule, but the direction of total output and welfare is generally indeterminate. So seeing that one group of consumers got a lower price does not let you declare a welfare improvement; at a minimum you must also estimate the demand response, whether the market served expanded, and whether other groups saw price increases.

In vertical relations, double marginalization predicts that vertical integration or an effective two-part tariff will lower the terminal price and raise sales; but integration may also bring foreclosure. If, while the merging party cuts price, rivals' access to inputs or their input prices are also worsening, then the merging party's price cut alone is not enough to separate efficiency from exclusion.

| Theory or institutional change | Comparative static | Joint empirical prediction | Main alternative explanation |
|---|---|---|---|
| Marginal cost rises | Price rises by the pass-through rate | Own price, sales, and rival price move together | Demand shock, menu cost |
| Products more differentiated | Bertrand competition softens | Markup rises, cross diversion falls | Quality changes marginal cost at the same time |
| Market becomes segmentable | Segment prices diverge by elasticity | Price, coverage, and composition all change | Selection, personalized quality |
| Vertical integration | Double markup internalized | Merging party's retail price falls, sales rise | Foreclosure, cost synergies |
| Repeated interaction more stable | Collusion more sustainable | Fewer price wars, more persistent punishment after shocks | Common costs, algorithmic synchronization |

The empirical predictions of this section should be written as a single joint response vector of several outcome variables, rather than fixating on one price coefficient alone. To truly identify conduct, you need exogenous cost shifters, demand rotators that can rotate the demand curve, ownership changes, or institutional experiments, forcing the candidate models to give distinct responses to the same shock.

## 5 Anchor Papers

The theory of pricing and competition only shows its weight when grounded in real markets. This section picks two classics, one that measured the theory of price discrimination against real data, and one that matched the game-theoretic characterization of collusion and price wars to a real cartel, each laid out by paper, method, results, significance, and limitations.

### 5.1 Leslie (2004)

::: {.case}
Paper and methodological position: "Price Discrimination in Broadway Theater," RAND Journal of Economics. It measures the second- and third-degree price discrimination theory of Section 3.2 of this chapter against real pricing data, and is a model of empirical price discrimination.

Method: Broadway theaters charge different prices for different seats and different purchase channels (full-price window, discount coupons, same-day discount booth) for the same performance, involving both third-degree discrimination by observable label and second-degree discrimination by self-selecting menu. The author builds and estimates a demand model, and on that basis simulates how box-office revenue and consumer welfare would change if pricing switched to a single price or to perfect price discrimination, thereby quantifying where the current discrimination scheme sits relative to these two extremes.

Results: Price discrimination significantly raised theater revenue relative to single pricing, on the order of a few percent, a sizable range; the effect on total social welfare, however, is small, and price discrimination mainly transfers surplus from consumers to the theater rather than growing the pie. This finding that "discrimination mainly transfers surplus and does not necessarily improve efficiency" is exactly the empirical confirmation of the welfare analysis of Section 3.2 of this chapter.

Significance and limitations: This paper demonstrates how to turn abstract price discrimination theory into an estimable structural model that supports counterfactuals, and is a forerunner of the structural pricing analysis later in this series. The limitation is that the conclusions depend on the functional form of demand and the specification of consumer heterogeneity, and the credibility of the counterfactuals (especially the unobservable extreme of perfect price discrimination) hinges on these specifications, which is the same class of problem as the sensitivity of welfare measurement to functional form discussed in Chapter 1.
:::

### 5.2 Porter (1983)

::: {.case}
Paper: "A Study of Cartel Stability: The Joint Executive Committee, 1880-1886," Bell Journal of Economics. It matches the game-theoretic characterization of collusion and price wars of Section 4.5 to a real railroad-freight cartel, and is a founding work of empirical collusion research and conduct tests.

Method: A 19th-century American railroad freight cartel (the Joint Executive Committee) at times maintained high freight rates and at times erupted into price wars. Green and Porter's theory holds that when the cartel cannot perfectly observe whether members are secretly cutting prices and can only observe prices contaminated by demand noise, rational collusion must be sustained by "entering a stretch of punitive price war whenever the price is abnormally low," so a price war is an equilibrium phenomenon that sustains collusion rather than a breakdown of it. Porter uses a switching regression (a regression that switches between the two states of collusion and price war) to estimate pricing behavior under these two states and tests whether the data fit this theory.

Results: The data support the pattern of alternating collusion and price wars, the estimated freight rates under the two states differ significantly, and rates in the price-war period are markedly lower. This provides empirical support for the counterintuitive game-theoretic conclusion that "a price war is part of collusion rather than its end," and it grounds the abstract "punishment phase" of the folk theorem of Chapter 2 in observable data.

Significance and limitations: This paper is the methodological starting point for identifying collusion and distinguishing competitive from collusive pricing, and leads directly to the conduct tests later in this series. The limitation is that it relies on a specific modeling of how the cartel operates, and distinguishing "a price war under collusion" from "competition itself" is always delicate in identification: observing a low price does not necessarily let you infer the behavioral pattern behind it.
:::

Together the two papers mark the empirical landing point of this chapter: Leslie measures price discrimination against real pricing and confirms that discrimination mainly transfers surplus, and Porter measures collusion and price wars against a real cartel and confirms that a price war can be part of collusion. Both foreshadow the core task of the structural industrial organization to come, inferring firms' pricing behavior and market power from observed prices and quantities.

## 6 Running Case: Pricing and Competition in the Cirrus Market

Now let us run the models of Sections 3 and 4 fully on the Cirrus market. Baseline demand $p = 12 - q$, marginal cost $c = 2$. Every number can be done by hand, and the R code below transcribes and checks the by-hand steps exactly.

### 6.1 Single-Price Monopoly and the Lerner Index

Set $MR = MC$: $12 - 2q = 2$, giving $q_m = 5$, $p_m = 7$, profit 25. Consumer surplus $\frac{1}{2}(12 - 7)\times 5 = 12.5$, deadweight loss $\frac{1}{2}(10 - 5)\times(7 - 2) = 12.5$. Check the Lerner index: markup $(7 - 2)/7 = 5/7 \approx 0.714$, demand elasticity at this point $|\varepsilon| = |dq/dp|\cdot p/q = 1\times 7/5 = 1.4$, and its inverse $1/1.4 \approx 0.714$ equals the markup, so the formula of Section 2 checks out numerically.

```r
a <- 12; b <- 1; cost <- 2
qm <- (a - cost)/(2*b); pm <- a - b*qm; profit_m <- (pm - cost)*qm
qc <- (a - cost)/b                         # competitive quantity (p=c)
CS_m <- 0.5*(a - pm)*qm; DWL <- 0.5*(qc - qm)*(pm - cost)
lerner <- (pm - cost)/pm; elas <- (pm/qm)  # |eps| = 1 * p/q
c(qm = qm, pm = pm, profit_m = profit_m, CS_m = CS_m, DWL = DWL,
  lerner = lerner, inv_elas = 1/elas)
#>       qm       pm profit_m     CS_m      DWL   lerner inv_elas
#>  5.00000  7.00000 25.00000 12.50000 12.50000  0.71429  0.71429
```

### 6.2 Third-Degree Price Discrimination

Suppose Cirrus can split merchants into two observable segments, with demand $p_1 = 12 - q_1$ (high stickiness, cannot do without the platform) and $p_2 = 8 - q_2$ (more price-sensitive) respectively, and marginal cost 2 in both. Setting marginal revenue equal to marginal cost in each segment: segment 1 gives $q_1 = 5$, $p_1 = 7$, and segment 2 gives $q_2 = 3$, $p_2 = 5$. The inelastic segment 1 is charged the higher price, exactly the inverse elasticity rule: segment 1's $|\varepsilon| = 7/5 = 1.4$ versus segment 2's $|\varepsilon| = 5/3 \approx 1.67$, markup $5/7$ versus $3/5$, the more inelastic the higher the markup.

```r
q1 <- (12 - cost)/2; p1 <- 12 - q1        # segment 1
q2 <- (8  - cost)/2; p2 <- 8  - q2        # segment 2
c(p1 = p1, e1 = p1/q1, markup1 = (p1 - cost)/p1,
  p2 = p2, e2 = p2/q2, markup2 = (p2 - cost)/p2)
#>       p1       e1  markup1       p2       e2  markup2
#>  7.00000  1.40000  0.71429  5.00000  1.66667  0.60000
```

### 6.3 Second-Degree Price Discrimination: A Tiered Commission Menu

Now to the core of this chapter's running case. Cirrus cannot tell which type a merchant belongs to; it only knows there are two types of merchants in the market, a high-activity type $\theta_H = 10$ and a low-activity type $\theta_L = 6$, each half of the market. Each merchant's gross value from usage $q$ is $\theta q - q^2/2$, paying total price $T$ (the platform's marginal cost is set to 0 to focus on the discrimination mechanism). Cirrus designs a two-tier menu $\{(q_L, T_L), (q_H, T_H)\}$ for merchants to self-select, which must satisfy the IC and IR of both types.

First set up the single-price benchmark for comparison, noting that it must be computed within the same two-type characterization to be comparable with the menu. If Cirrus can only offer a single plan rather than a menu, it earns at most 25: serving both types requires pushing the plan down to a level the low type will accept, dragged down by the low end (the optimal single plan serving both types earns only 18), so it is better to simply offer only the high type a plan of $q = 10$ priced at 50, with a one-half chance of a sale and an expectation of exactly 25. This 25 is the ceiling for a single price in this two-type market, and it happens to equal the monopoly profit under linear demand in Section 6.1. Whether the menu can break through it is the value of second-degree discrimination.

By the screening solution of Chapter 3, the high type gets the efficient quantity $q_H = \theta_H = 10$ (no distortion at the top), and the low type's quantity is held down to $q_L = \theta_L - \frac{\mu}{1-\mu}(\theta_H - \theta_L) = 6 - 1\times 4 = 2$ (distortion at the bottom, well below its efficient quantity 6). On pricing, the low type's IR binds, $T_L = \theta_L q_L - q_L^2/2 = 12 - 2 = 10$; the high type's IC binds, and it is left an information rent $(\theta_H - \theta_L)q_L = 4\times 2 = 8$, so its payment is $T_H = (\theta_H q_H - q_H^2/2) - 8 = 50 - 8 = 42$. Cirrus's expected profit $\frac{1}{2}\times 42 + \frac{1}{2}\times 10 = 26$ breaks through the single-price ceiling of 25.

```r
thH <- 10; thL <- 6; mu <- 0.5
qH <- thH                                    # no distortion at the top
qL <- thL - (mu/(1 - mu))*(thH - thL)        # distortion at the bottom
rentH <- (thH - thL)*qL                      # high type's information rent
TL <- thL*qL - qL^2/2                         # low type IR binds
TH <- (thH*qH - qH^2/2) - rentH               # high type IC binds
profit2 <- mu*TH + (1 - mu)*TL
# self-check: IC and IR
uH_own   <- thH*qH - qH^2/2 - TH
uH_mimic <- thH*qL - qL^2/2 - TL
uL_own   <- thL*qL - qL^2/2 - TL
c(qL = qL, qH = qH, rentH = rentH, TL = TL, TH = TH, profit2 = profit2,
  uH_own = uH_own, uH_mimic = uH_mimic, uL_own = uL_own)
#>      qL      qH   rentH      TL      TH profit2  uH_own uH_mimic  uL_own
#>       2      10       8      10      42      26       8        8       0
```

The self-check confirms the menu holds: the high type gets utility 8 by choosing its own tier, and only 8 by pretending to be the low type (IC binds exactly), and the low type gets 0 by choosing its own tier (IR binds). Compared with first-degree discrimination under full information, where Cirrus could price to each type, Cirrus could have earned $\frac{1}{2}\times\frac{\theta_H^2}{2} + \frac{1}{2}\times\frac{\theta_L^2}{2} = \frac{1}{2}\times 50 + \frac{1}{2}\times 18 = 34$, and information asymmetry costs it 8 in lost profit; this 8 is exactly the price of information, part flowing as information rent to the high type and part evaporating into deadweight loss from the distortion at the bottom. Second-degree discrimination lifts profit from the single-price 25 to 26, approaching but not reaching the full-discrimination 34, and this gap is the ceiling that private information sets for the platform.

### 6.4 Oligopoly: Cournot, Bertrand, and Differentiation

Now let the rival platform Nimbus enter. If the two play Cournot quantity competition, the symmetric equilibrium has each producing $q = 10/3$, total quantity $20/3$, price $16/3 \approx 5.33$, profit per firm $100/9 \approx 11.11$, industry total profit $200/9 \approx 22.2$, below the monopoly's 25. If the two offer an undifferentiated service and play Bertrand price competition, the price is pushed down to marginal cost 2 and profit goes to zero. If the two services are differentiated (Hotelling, differentiation sensitivity $t = 1$), the equilibrium price rises back to $c + t = 3$ and profit per firm is 0.5. The three outcomes mark out the spectrum of competition: an undifferentiated price war is the most brutal (price 2), quantity competition is mild (price 5.33), and differentiation rescues the price partway back from marginal cost (price 3).

```r
# Cournot: p = 12 - Q, two symmetric firms
q_c  <- (a - cost)/(3*b); Q <- 2*q_c
p_cournot <- a - b*Q; profit_each <- (p_cournot - cost)*q_c
# Bertrand (undifferentiated): p = MC
p_bertrand <- cost
# Hotelling (quadratic transport cost, t=1): p = c + t
t <- 1; p_hotelling <- cost + t; profit_hot <- (p_hotelling - cost)*0.5
c(cournot_p = p_cournot, cournot_profit_each = profit_each,
  bertrand_p = p_bertrand, hotelling_p = p_hotelling, hotelling_profit = profit_hot)
#>          cournot_p cournot_profit_each          bertrand_p         hotelling_p
#>           5.333333           11.111111            2.000000            3.000000
#>   hotelling_profit
#>           0.500000
```

### 6.5 Nonlinear Pricing: Two-Part Tariffs and Bundling

Switch pricing instrument and look at the ceiling on profit again. If Cirrus faces a group of homogeneous merchants, each with demand $p = 12 - q$, it need not charge a single price and can use a two-part tariff instead: push the per-unit price down to marginal cost $p = c = 2$, let each merchant buy the efficient quantity $q = 10$, then charge a membership fee equal to the entire consumer surplus $50$. Profit per merchant is therefore $50$, twice the single-price monopoly $25$ of Section 6.1, and output reaches the competitive level with no deadweight loss.

```r
c_ <- 2
q_star <- 12 - c_                 # efficient quantity at p = c
CS <- 0.5*(12 - c_)*q_star        # entire consumer surplus at this price
F  <- CS                          # fixed fee sweeps up all the surplus
c(q = q_star, F = F, profit_two_part = F + (c_ - c_)*q_star, profit_single = 25)
#>               q               F profit_two_part   profit_single
#>              10              50              50              25
```

Bundling is another family of tools. Suppose Cirrus has two value-added services A and B with near-zero marginal cost, and two types of merchants whose valuations of them are exactly negatively correlated: one type values A (valuation 10) and not B (valuation 2), and the other is the reverse. Sold separately, each service's optimal price is 10, sold only to the half that values it, for a combined revenue of 20; bundling A and B into one package priced at 12, both types buy, for a revenue of 24. Bundling smooths out the difference in willingness to pay between the two types and extracts one more layer.

```r
vA <- c(10, 2); vB <- c(2, 10)                          # two types' valuations of A and B
revA <- max(2*sum(vA >= 2), 10*sum(vA >= 10))           # price 2 selling two vs price 10 selling one
revB <- max(2*sum(vB >= 2), 10*sum(vB >= 10))
c(separate = revA + revB, pure_bundle = 12*sum((vA + vB) >= 12))
#>    separate pure_bundle
#>          20          24
```

### 6.6 Dynamic Pricing: The Rent Path of Hotelling's Rule

Treat one of Cirrus's scarce resources (say a limited stock of premium ad slots) as an exhaustible resource, with marginal placement cost $c = 2$, interest rate $r = 0.05$, and initial scarcity rent $4$ (initial price $6$). Hotelling's Rule says the rent must grow at $r$ and its present value stays constant:

```r
r <- 0.05; rent0 <- 4; cext <- 2; t <- 0:2
rent <- rent0*(1 + r)^t; price <- cext + rent; pv <- rent/(1 + r)^t
rbind(t = t, rent = round(rent, 4), price = round(price, 4), pv_rent = round(pv, 4))
#>         [,1] [,2] [,3]
#> t          0  1.0 2.00
#> rent       4  4.2 4.41
#> price      6  6.2 6.41
#> pv_rent    4  4.0 4.00
```

The rent rises from $4$ to $4.2$ and $4.41$, growing by exactly $5\%$ each period, while its present value stays pinned at $4$: this is exactly the numerical face of the no-arbitrage between "extract now and invest" and "keep it and sell later." Dynamic pricing of scarce capacity, letting its shadow rent chase the cost of capital, is the modern version of the same rule.

### 6.7 Vertical Relations: The Cost of Double Marginalization

Finally, let Cirrus no longer price as an integrated firm but split into an upstream service provider (marginal cost 2) plus a downstream retail end. Upstream supplies at wholesale price $w$, and downstream, facing marginal cost $w$, adds another monopoly markup. Upstream foresees downstream's reaction and chooses $w$ to maximize its own profit, solving to $w = 7$, so the retail price is pushed by two markups to $9.5$, far above the integrated $7$, and the combined upstream-plus-downstream profit is only $18.75$, short of the integrated $25$:

```r
qm <- (12 - 2)/2; pm <- 12 - qm               # integrated monopoly
w <- 7; qd <- (12 - w)/2; pd <- (12 + w)/2    # separated: downstream MR=w -> q=(12-w)/2
c(integrated_price = pm, integrated_profit = (pm - 2)*qm,
  wholesale = w, retail_price = pd, joint_profit = (w - 2)*qd + (pd - w)*qd)
#>  integrated_price integrated_profit         wholesale      retail_price
#>              7.00             25.00              7.00              9.50
#>      joint_profit
#>             18.75
```

The retail price is higher and the two firms earn less, with consumers and firms worse off at the same time: this is double marginalization. Any one of vertical integration, a two-part wholesale contract, or RPM can pull the outcome back to the integrated optimum of $p = 7$ and profit $25$.

The accounts of this section now balance, and they finally settle that 25 caught in the middle from Section 1. In the same Cirrus market, facing two types of merchants, single-price monopoly earns 25, the self-selecting menu of second-degree discrimination lifts it to 26, and first-degree discrimination under full information can in theory reach 34 (the upper bound set by information); switch to homogeneous merchants and the two-part tariff instrument, and it can sweep up all the surplus 50 at once (the upper bound set by the pricing instrument); once a rival enters, Cournot competition pushes the price down to 5.33 and profit to 22.2, differentiation gives a Hotelling price of 3, and an undifferentiated Bertrand price war zeroes out profit (the lower bound set by competition); and once split into a vertical chain, double marginalization pushes the price past integrated monopoly and total profit actually falls to 18.75. That 25 is neither a ceiling nor a floor; it is just the reading of the one cell "single price plus exclusive integration," and the pricing instrument, competitive structure, and vertical organization together determine where in this wide range the platform ends up.

## 7 Failure Modes and Robustness in Pricing Analysis

The pricing and competition models of this chapter are clean and handy, but there are several places where applying them to real platforms easily goes astray, and this section lays them out one at a time.

The welfare of price discrimination is ambiguous and cannot be taken for granted. Intuitively discrimination "harms consumers," but this judgment is both incomplete and often wrong. Third-degree discrimination relative to single pricing may raise or lower total social output and total welfare, with the direction depending on the curvature of the segments' demand, and there is no universal sign; its reliable effect is only to shift surplus from consumers to the firm (Leslie's empirical work is exactly this). More counterintuitive is that first-degree discrimination is the most efficient (largest output, zero DWL) yet the harshest on consumers. And one more layer is often overlooked: discrimination can sometimes let a low-valuation group that a single high price would shut out afford the good, thereby expanding the market and improving their welfare. So judging price discrimination cannot rest on "some people paid more"; you have to compute output, welfare, and distribution separately, consistent with the discipline of Chapter 1 that stresses defining the welfare measure first and passing judgment after.

The fragility of the Bertrand paradox reminds us that the model is extremely sensitive to its assumptions. The startling conclusion that price competition reproduces perfect competition rests entirely on the premises of completely undifferentiated products, no capacity constraints, and one-shot play, and Sections 4.3 through 4.5 show that relaxing any one (add capacity, add differentiation, add repetition) can lift the price back up substantially. This means treating Bertrand zero profit as the general prediction for a duopoly is dangerous, and real duopoly profit depends heavily on the specific form of competition. When analyzing platform competition, first getting clear on what dimension the firms actually compete in, whether there are capacity constraints, and whether they coexist over the long run matters far more than applying some single model.

Collusion and competition are hard to tell apart from observation, a perennial headache of empirical work. Section 4.5 said that behind a high price could be collusion, differentiation, or Cournot under capacity constraints, and looking at price data alone cannot tell them apart. Porter's study is a classic precisely because it achieves identification only with the help of the extra structure of the cartel switching between collusion and price war. This reminds us that inferring behavior from price (so-called conduct) is never something you can conclude from one glance at a markup; it needs extra assumptions or institutional evidence, which is the core identification problem the conduct tests later in this series will confront head-on.

Two-sided markets are an important complexity of platform pricing that this chapter has not developed. The pricing models of this chapter assume by default that the firm faces a single group of buyers, but platforms often connect two sides (buyers and sellers, users and advertisers), where the size of one side affects the value of the other, and this cross-side network effect can make optimal pricing deviate sharply from single-sided intuition, with the typical conclusion that a platform may charge one side below marginal cost, or even subsidize it, and earn it back on the other side. The single-sided framework of this chapter is the necessary foundation, but applying it directly to a two-sided platform yields the wrong pricing prescription. Two-sided pricing is the subject of the platform chapter later in this series, where the pricing formulas under cross-side effects will be given; here we need only know that single-sided optimal intuition can reverse entirely in a two-sided world.

The durable-goods and commitment problem (the Coase conjecture of Section 3.3) is especially real for digital products. Software, memberships, and digital content are often durable, and sellers can also reprice extremely frequently, so Coase's logic bites platforms particularly easily, forcing their pricing power to leak toward the competitive level. A platform's responses (subscriptions, time-limited pricing, version management) are all essentially about manufacturing commitment to fight its own future incentive to cut prices, and when analyzing digital-product pricing, whether a credible commitment can be established often determines the outcome more than the static monopoly-pricing formula.

Stringing these together: this chapter gives a clear set of models characterizing market power, price discrimination, and oligopoly, but their conclusions all hang on specific assumptions: the welfare sign of discrimination depends on demand curvature, Bertrand's zero profit depends on undifferentiated products and no capacity, collusion and competition are hard to distinguish in the data, and single-sided pricing reverses in a two-sided world. Robust pricing analysis is not about applying the right formula, but about recognizing your own market's dimension of competition, information structure, and commitment power, and then choosing the appropriate model. This is the same attitude that runs through this whole series: the credibility of a method lies not in the elegance of the technique, but in putting the assumptions on the table and spelling out how the conclusions move when those assumptions are relaxed. With this, the microeconomic-theory foundation of Part I is laid, from consumers and welfare, to games and information, to pricing and competition, and the identification and structural chapters to come will keep returning here to recover the behavioral primitives they need to estimate.

## 8 Further Reading

::: {.readings}
Required reading, in suggested reading order:

- Mas-Colell, Whinston and Green, Microeconomic Theory, Chapter 14. The authoritative reference on monopoly, price discrimination, and oligopoly, the full source for Sections 3 and 4 of this chapter, where the formal treatments of Cournot, Bertrand, and price discrimination all live.
- Nolan Miller, Notes on Microeconomic Theory, Chapter 9 (Monopoly). The intuition-building groundwork for monopoly pricing and price discrimination, more patient than MWG, good to read first to build the geometric picture.
- Tirole, The Theory of Industrial Organization, Chapters 1 to 5. The standard graduate textbook of industrial organization, with the fullest treatment of price discrimination, oligopoly, and dynamic pricing, where every section of this chapter can be deepened.

Further reading:

- Leslie (2004). The model of structural empirical price discrimination, essential reading for the information systems and platform pricing direction; focus on the construction of the two counterfactuals, single pricing and perfect discrimination.
- Porter (1983). The empirical classic on collusion and price wars, to be read against the folk theorem of Chapter 2 to understand why a price war can be part of collusion.
- Kreps and Scheinkman (1983). The classic on Bertrand competition under capacity constraints reproducing Cournot, the key to reconciling the two oligopoly models.
- d'Aspremont, Gabszewicz and Thisse (1979). The quadratic transport cost correction to the Hotelling model, the rigorous source of the $p = c + t$ result of Section 4.4 of this chapter.
- Coase (1972). The origin of the durable-goods monopoly and commitment problem, the source of the Coase conjecture.
- Oi (1971) and Adams and Yellen (1976). The two founding works on two-part tariffs and bundling, the source of the nonlinear pricing of Section 3.4, direct theory for SaaS and platform pricing.
- Shapiro and Varian (1999), Information Rules. The source of the framework for information-good versioning and quality menus, on which Section 3.5 is based.
- Hotelling (1931) and Stokey (1979). The former is the origin of intertemporal pricing of exhaustible resources (scarcity rent growing at $r$), the latter proves that a durable-goods monopolist that can commit optimally holds a constant price; read together they clarify the two poles of the dynamic pricing of Section 3.6.
- Spengler (1950). The original argument for double marginalization and vertical integration, on which Section 4.6 is based.
- Dixit (1980), Milgrom and Roberts (1982), and Fudenberg and Tirole (1984). The three cornerstones of entry deterrence: capacity commitment, the signaling version of limit pricing, and the strategic taxonomy of over/under-investment, the full source of Section 4.7.
:::

::: {.apa-refs}
- Adams, W. J., & Yellen, J. L. (1976). Commodity bundling and the burden of monopoly. *The Quarterly Journal of Economics, 90*(3), 475-498. https://doi.org/10.2307/1886045
- Coase, R. H. (1972). Durability and monopoly. *Journal of Law and Economics, 15*(1), 143-149. https://doi.org/10.1086/466731
- d'Aspremont, C., Gabszewicz, J. J., & Thisse, J.-F. (1979). On Hotelling's "stability in competition". *Econometrica, 47*(5), 1145-1150. https://doi.org/10.2307/1911955
- Dixit, A. (1980). The role of investment in entry-deterrence. *The Economic Journal, 90*(357), 95-106. https://doi.org/10.2307/2231658
- Fudenberg, D., & Tirole, J. (1984). The fat-cat effect, the puppy-dog ploy, and the lean and hungry look. *American Economic Review, 74*(2), 361-366.
- Hotelling, H. (1931). The economics of exhaustible resources. *Journal of Political Economy, 39*(2), 137-175. https://doi.org/10.1086/254195
- Kreps, D. M., & Scheinkman, J. A. (1983). Quantity precommitment and Bertrand competition yield Cournot outcomes. *The Bell Journal of Economics, 14*(2), 326-337. https://doi.org/10.2307/3003636
- Leslie, P. (2004). Price discrimination in Broadway theater. *RAND Journal of Economics, 35*(3), 520-541. https://doi.org/10.2307/1593706
- Mas-Colell, A., Whinston, M. D., & Green, J. R. (1995). *Microeconomic theory*. Oxford University Press.
- Milgrom, P., & Roberts, J. (1982). Limit pricing and entry under incomplete information: An equilibrium analysis. *Econometrica, 50*(2), 443-459. https://doi.org/10.2307/1912637
- Miller, N. H. (2006). *Notes on microeconomic theory* [Unpublished lecture notes]. Harvard Kennedy School.
- Oi, W. Y. (1971). A Disneyland dilemma: Two-part tariffs for a Mickey Mouse monopoly. *The Quarterly Journal of Economics, 85*(1), 77-96. https://doi.org/10.2307/1881841
- Porter, R. H. (1983). A study of cartel stability: The Joint Executive Committee, 1880-1886. *The Bell Journal of Economics, 14*(2), 301-314. https://doi.org/10.2307/3003634
- Spengler, J. J. (1950). Vertical integration and antitrust policy. *Journal of Political Economy, 58*(4), 347-352. https://doi.org/10.1086/256964
- Shapiro, C., & Varian, H. R. (1999). *Information rules: A strategic guide to the network economy*. Harvard Business School Press.
- Stokey, N. L. (1979). Intertemporal price discrimination. *The Quarterly Journal of Economics, 93*(3), 355-371. https://doi.org/10.2307/1883163
- Tirole, J. (1988). *The theory of industrial organization*. MIT Press.
:::
