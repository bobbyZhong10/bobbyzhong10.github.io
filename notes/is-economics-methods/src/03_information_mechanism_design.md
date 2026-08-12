---
title: "Information Economics, Contracts & Mechanism Design"
subtitle: "Private Information, Incentives, and the Design of Rules"
seriesline: "Foundations of Information Systems Economics · Chapter 3"
series: "Foundations of Information Systems Economics: Empirical Methods for Information Systems Economics"
sidebartitle: "Foundations of Information Systems Economics"
sidebarsub: "Theory · Methods · Applications"
chapterlabel: "Chapter 3 · Information Economics, Contracts & Mechanism Design"
---

## Introduction

A platform asks advertisers to report how much a click is worth to them, then hands the best ad slot to the highest reported value. The rule looks fair, except it skips one detail: the advertiser knows its own value and the platform does not, and when a single bid decides both the slot and the payment, telling the truth need not pay off. Any institution that asks participants to volunteer their private information, without first answering the question of whether lying is profitable, tends to collect strategy rather than information.

Asymmetric information wears two faces. What is hidden before contracting is a type: a low quality seller, a low risk borrower, and a high value advertiser can all try to disguise themselves. What is hidden after contracting is an action: effort, maintenance, and content moderation cannot be fully monitored. Signaling, screening, and contracts each address one of these difficulties, while mechanism design turns the problem around: instead of taking the rules as given and predicting an equilibrium, it starts from the outcome we want and searches for rules under which self interested participants cooperate of their own accord. Incentive compatibility and individual rationality are therefore not side conditions but the boundary within which an institution can function at all.

This chapter takes Vireo's ad slot auction as its through line. We start from a rule that cannot extract the truth, then compare VCG, the generalized second price auction, reserve prices, and revenue equivalence, and finally carry the same ideas into information design, cheap talk, and matching. The question running through the chapter is plain: the platform cannot order participants to hand over their true information, but can it arrange payments, menus, or signals so that telling the truth or self selecting turns out to be in their own interest? That question is also the first one worth asking when we later study platform governance, online advertising, rating systems, and digital contracts.

## 1 An Auction That Cannot Extract the Truth

::: {.case}
The search platform Vireo wants to sell the ad slots on a page of search results to advertisers. To shrink the problem to the smallest scale we can compute by hand, consider a concrete setup: there are two ad slots, the top slot delivers 100 clicks per month and the second slot 40 clicks (a higher position draws more clicks, and this gap in click volume is the position effect). Three advertisers have private per click values of 10, 6, and 4 yuan respectively (a larger number means the advertiser cares more about the exposure, and this value is known only to the advertiser itself). Vireo wants to give the slots to the advertisers who value them most, and to collect as much money as possible.
:::

Vireo's product manager proposes something straightforward: let the three advertisers each report a per click price they are willing to pay, give the top slot to the highest report and the second slot to the next highest, and charge each advertiser the price it reported. This is the pay-your-bid first price auction. He casually estimates the revenue: advertisers will report their true values 10, 6, and 4, so the top slot advertiser pays $10 \times 100 = 1000$, the second slot pays $6 \times 40 = 240$, and one auction brings in 1240 yuan.

Buried in that estimate is a fatal assumption: that advertisers will report truthfully. But a moment of thinking on the advertiser's behalf shows that when you pay what you report, nobody reports the truth. If the advertiser with value 10 truly reported 10, it would win the top slot but hand over all of its surplus, for a net payoff of zero. It will clearly shave its report down, low enough that it still just wins, while keeping as much surplus as possible for the position it lands. Exactly how far it shaves depends on how it guesses its rivals will report, and the rivals are running the same calculation, so bidding turns into a game of mutual second guessing in which nobody's report equals anybody's value. The product manager's 1240 rests entirely on the sand of "the report equals the value." Once we admit that advertisers will shade, that number can neither be computed accurately nor extract what he really wants to know: who values these slots the most.

The root of the problem is that the first price auction gives advertisers no reason to report truthfully. What Vireo really needs is a set of rules under which "reporting the true value" is itself each advertiser's best choice, no matter what the others report. Does such a rule exist? If it does, how much money does it collect, and is it the one that collects the most? These questions are the heart of mechanism design, and this chapter answers them one by one. Section 3.6 gives a VCG mechanism that makes truthful reporting a dominant strategy, and Section 6 applies it to Vireo's numerical example and finds that it should really collect only 680 yuan, far less than that inflated 1240, but a defensible number built on the advertisers' truth. Truth is worth more than an inflated figure on paper: that is the first idea this chapter aims to establish.

## 2 Mechanisms and Constraints Under Private Information

To analyze asymmetric information, we first write a mechanism down in full. Take adverse selection as the skeleton: each participant has a private type $\theta_i$, encoding information known only to itself (value, quality, cost), and the prior distribution of types is common knowledge (a reminder: throughout this chapter $\theta$ denotes a private type, unlike the convention elsewhere in the book where $\theta$ denotes a structural parameter, so read it from context). A mechanism specifies what messages participants may send (for instance, reporting a price), and the outcome under each combination of messages (including who gets what allocation, and who pays how much in transfers). Participants choose messages to maximize their own payoff given their types, so the outcome of the mechanism is the equilibrium of a game whose inputs are the types.

Two constraints run through all of mechanism design, and we must state them clearly up front. The first is incentive compatibility (IC): the mechanism must make each participant, out of self interest, make the choice the designer wants. In a direct mechanism (where participants report their types directly), IC means that reporting your type truthfully is no worse than misreporting as any other type. The second is individual rationality (IR, also called the participation constraint): participating voluntarily is no worse than not participating, or the participant will drop out. Almost every mechanism design problem optimizes some objective (efficiency or revenue) subject to the IC and IR constraints, and this "objective plus two classes of constraints" structure is the unifying template for all the analysis in this chapter, of a piece with the discipline of the previous chapter, which insisted on fixing the target quantity before discussing implementation.

The objective has two common choices, corresponding to two kinds of mechanism. A mechanism that pursues efficiency gives the allocation to the people who value it most, maximizing total surplus, and VCG is its representative. A mechanism that pursues revenue collects the most for the designer subject to IC and IR, and Myerson's optimal auction is its representative. The two generally do not coincide: the efficient mechanism need not collect the most revenue, and the revenue optimal mechanism often sacrifices a bit of efficiency (for instance, by setting a reserve price that keeps low value bidders out). The revenue equivalence and optimal auction results in Section 4 will make this tension precise.

The primitives of moral hazard differ slightly, but the template is the same. Here what is private is not a type but an action (usually effort), which affects the probability distribution of the outcome but is itself unobservable. The designer can only pay based on the observable outcome, so IC becomes "make the participant voluntarily choose the effort level the designer wants," and the constraint is written on the marginal return to effort. Section 3.4 develops this.

To summarize: a mechanism consists of messages and outcomes, participants best respond according to their private types or actions, and all mechanism design optimizes efficiency or revenue subject to IC (which aligns self interested choices with the designer's intent) and IR (which keeps participants in voluntarily). We now work through the two broad classes of asymmetric information one by one, in the most detailed section of the chapter.

## 3 Core Structure: From Adverse Selection to Mechanism Design

This section dissects the core models one by one, organized by the type of asymmetric information and by who moves first, and finally converges on the two grand theorems of mechanism design: the revelation principle and VCG.

### 3.1 Adverse Selection and Market Unraveling

The most basic consequence of asymmetric information is that it can make an entire market that should exist disappear. Akerlof's used car market is the prototype of this idea, and it holds equally well in platform language. Imagine a resale platform for used equipment, where the quality $\theta$ of a seller's equipment is, from the buyer's perspective, uniformly distributed on $[0, 1]$. The seller knows the exact $\theta$ and values it at exactly $\theta$, while a buyer who acquires equipment of quality $\theta$ is willing to pay $1.5\theta$ (the buyer puts it to better use, so values it more, and there should be room to trade).

The trouble is that the buyer cannot see $\theta$ and can only bid based on average quality. At price $p$, the only sellers willing to sell are those whose value does not exceed $p$, that is, the batch with $\theta \leq p$, whose average quality is $p/2$. A rational buyer anticipates this and is willing to pay only $1.5 \times (p/2) = 0.75p$, strictly below the asking price $p$. So at any positive price the buyer refuses to trade, and the market shrinks all the way down to only the worst equipment with $\theta = 0$, with the entire block of mutually beneficial trade eaten away by adverse selection. This is market unraveling: the presence of low quality sellers squeezes high quality sellers out of the market, because price cannot tell the two apart. The platform analogue is everywhere. Once low quality sellers, order faking merchants, and fake reviews cannot be distinguished from the good ones, they drive the good ones out, and this is the fundamental reason platforms pour so many resources into building rating, certification, and guarantee systems: at bottom, these are all information mechanisms for fighting adverse selection.

### 3.2 Signaling: The Informed Party Speaks Up

There are two directions for fighting adverse selection, depending on who breaks the asymmetry. The first is to let the party holding the private information do something to prove itself, and this is signaling, with Spence's education signaling model as the prototype.

The core mechanism is this: if an action is less costly for a high type than for a low type, then that action can credibly separate the two. An equipment seller offering a long warranty is the platform version of a signal. For a seller of good equipment, the warranty almost never triggers a payout and is cheap; for a seller of bad equipment, the warranty will most likely have to pay out and is expensive. So only good sellers are willing to offer a long warranty, and the warranty itself becomes a credible signal of quality, even though the buyer cannot see quality but can see the warranty. The key here is the single-crossing property: the cost of the signal varies systematically with type, so high and low types make different tradeoffs over "how strong a signal to send," which lets the signal separate them.

The equilibria of a signaling game come in two kinds. In a separating equilibrium different types send signals of different strength, and the buyer can perfectly infer type from the signal; in a pooling equilibrium all types send the same signal, and the signal carries no information. Which equilibrium arises depends on the cost structure and on the buyer's off-path beliefs (what to believe upon seeing a signal never observed), and the latter is exactly the object that, as noted in Section 3.5 of the previous chapter, needs refinements such as sequential equilibrium to pin down. Signaling has a welfare implication worth emphasizing: it is often pure waste. To separate itself from the low quality types, a high quality type must incur a heavy cost to send the signal (excessive warranties, excessive investment), and this cost creates no real value, spent purely to overcome the asymmetry, a social cost that adverse selection imposes on the market.

### 3.3 Screening: The Uninformed Party Designs a Menu

The second direction for fighting adverse selection runs the other way: let the party without the information design a self selection menu that induces each type to reveal itself, and this is screening. Its only difference from signaling is who moves first: in signaling the informed party speaks first, and in screening the uninformed party puts out a menu first.

Screening is everywhere in platform settings. A cloud service platform faces two classes of customers, high usage and low usage, but cannot tell who is who, so instead of setting a single price it lays out a tiered menu (a basic plan with a usage cap and a low price, a professional plan with no cap and a high price) and lets customers choose, and high usage customers naturally choose the professional plan. For this menu to work, it must satisfy IC and IR: each class of customer must find that choosing the tier designed for it is no worse than choosing another tier (IC) and no worse than choosing nothing (IR). Solving these constraints forces out a core concept that recurs throughout mechanism design: information rent.

The logic of information rent goes like this: to keep the high type customer from impersonating the low type and choosing the cheap tier, the platform must leave the high type a little extra benefit (otherwise it would go and grab the benefit of the low price tier), and this benefit the platform is forced to concede is the information rent, the rent the holder of private information captures by virtue of its informational advantage. The result is a general tradeoff running through contract theory, the rent-efficiency tradeoff: to hold down the information rent it must pay the high type, the designer deliberately distorts the allocation given to the low type (for instance, capping the basic plan's usage below the efficient level), trading an efficiency loss for a rent saving. This structure of "distort the bottom to squeeze the top" will be computed in detail in Section 6 and later in this series when we cover second degree price discrimination, and it is also the key to understanding why platforms so often make their low end products deliberately crippled.

![Single crossing and information rent in screening: the reservation utility curves of two customer classes ($t = \theta V(q)$, steeper for the high type). The low type receives the distorted, pushed down $q_L$ and sits on its own IR curve (no rent); the high type receives the efficient quantity $q_H$ and, to satisfy IC, is left a slice of information rent (its payment $t_H$ falls below the height of its own IR curve at $q_H$, the vertical gap in the figure).](assets/fig/fig_03_screening.svg)

::: {.intuition}
Signaling and screening are two sides of the same coin, and grasping their division of labor saves a lot of confusion: who has the information and who moves first determines which language to use. The informed party moving first and burning money to prove itself is signaling (a job seeker earning a degree, a seller offering a warranty); the uninformed party moving first and laying out a menu to induce self selection is screening (an employer setting a pay ladder, a platform tiering its service, a monopolist doing price discrimination). Both rely on single-crossing to make types separable, and both incur a cost: the cost of signaling is the good type's money burning, and the cost of screening is the distortion of the low end allocation plus the information rent at the high end. The platform is both the stage for signaling and the designer of screening, and seeing clearly which side of the coin it stands on at any moment gets you half the analysis for free.
:::

### 3.4 Moral Hazard: Incentivizing an Unseen Action

Everything above is hidden information before contracting. Moral hazard is a hidden action after contracting: once the contract is signed, the effort one party takes is unseen by the other, who observes only the outcome, which is shaped by effort but also laced with luck. In the classic principal-agent problem, a principal (say, a platform) hires an agent (say, an operations team) to do work, the agent's effort affects the probability distribution of performance, but effort itself is unobservable, so the platform can only pay based on observable performance.

Let us make the mechanism concrete with the smallest example we can compute by hand. The operations team can choose high effort or low effort, with high effort costing 5 and low effort costing 0. Under high effort the probability of high performance (worth 100) is 0.8, and under low effort only 0.2; when performance is low the value is 0. The platform pays by performance, $w_H$ when performance is high and $w_L$ when it is low, subject to limited liability ($w \geq 0$, so it cannot charge the team). If the platform wants to induce high effort, it must make high effort no worse than low effort for the team, which is the IC constraint:

$$0.8\, w_H + 0.2\, w_L - 5 \;\geq\; 0.2\, w_H + 0.8\, w_L$$

Rearranging gives $0.6(w_H - w_L) \geq 5$, that is $w_H - w_L \geq 25/3$. The pay gap between good and bad performance must be large enough for the team to have any motive to exert effort. Under limited liability the platform pushes $w_L$ down to 0, so $w_H = 25/3 \approx 8.33$. Substituting back, the team's expected payoff is $0.8 \times 25/3 - 5 = 5/3 \approx 1.67 > 0$, strictly positive, a rent the team collects for free.

This small example distills two general lessons of moral hazard. The first is the risk-incentive tradeoff: to incentivize effort you must make pay vary with performance (widening the gap between $w_H$ and $w_L$), but performance is laced with luck, so making pay track performance means loading risk onto a risk averse agent, and the stronger the incentive the more risk imposed, so the optimal contract compromises between incentives and insurance. The second is the limited-liability rent: when you cannot fine the agent, widening the pay gap enough to satisfy IC often forces you to leave the agent a rent it collects even if it does nothing extra, which is exactly the 5/3 above. Commission rates, ratings, bans, and deposit schemes on platforms are, at bottom, all incentive designs for solving moral hazard under limited liability.

Generalizing effort from two levels to a continuum is what makes clear how the shape of the optimal contract is determined by the performance signal. Let effort $e$ be a continuous variable with an increasing, convex cost $c(e)$, let performance $x$ have density $f(x \mid e)$, and let the agent's compensation contract be a curve $w(x)$ that varies with performance. Solving the inner optimization "which $e$ will the agent choose" directly is hard, and the standard move is the first-order approach: under suitable regularity conditions, replace the entire IC constraint with the first order condition of the agent's effort choice, $\int u(w(x))\,\partial_e f(x\mid e)\,dx = c'(e)$, collapsing the two level problem into a single level problem solvable by Lagrange methods. The optimal contract that comes out has a beautiful characterization, Holmström's (1979) informativeness principle: the optimal pay should hinge on the likelihood ratio $\partial_e f(x\mid e)/f(x\mid e)$, which measures the extent to which observing performance $x$ is evidence of high effort rather than luck. Its practical import is strong: any additional signal that reduces the noise in inferring effort, even if it creates no value directly, should be written into the contract (which is exactly why, when a platform evaluates a seller, it folds in side evidence like return rates, response times, and peer benchmarks); conversely, a metric unrelated to effort and purely noisy should not enter pay. The first-order approach has its fragility too: when the performance distribution fails conditions like the monotone likelihood ratio, replacing IC with the first order condition can break down and give a false optimum, a boundary Section 7 will flag.

### 3.5 The Revelation Principle: Why Direct Truthful Mechanisms Are Enough

The space of mechanisms is frighteningly large: bids, menus, multi round games, arbitrarily complex message rules, and the designer can choose any of them. The revelation principle is a powerful result that abruptly shrinks this space, and it says: any outcome a mechanism can implement in equilibrium can also be implemented by a mechanism that is direct (participants report only their own types) and truthful (reporting truthfully is an equilibrium).

The intuition is unexpectedly simple. Suppose some complex mechanism has an equilibrium in which participants follow some equilibrium strategies given their types. Now build a direct mechanism: let participants report only their types, and then have the designer execute, on their behalf, the equilibrium strategy corresponding to that type in the original mechanism. Since acting according to the true type was already an equilibrium in the original mechanism, reporting the type truthfully is an equilibrium in this new mechanism, and the outcome is exactly the same as the original. In other words, the designer can internalize "the strategy the participants would have played" into the mechanism and execute it itself, leaving participants with nothing to do but report their types truthfully.

The use of this principle is methodological: to argue whether some objective can be achieved and what the optimal mechanism looks like, we need not search across all the strange and exotic mechanisms, but only within the small class of "direct and truthfully incentive compatible" mechanisms, because anything any other mechanism can do, this small class can do too. Section 4's derivation of revenue equivalence and Myerson's derivation of the optimal auction both rely on this contraction. It also explains why the analysis of mechanism design always starts from "let participants report their types truthfully": that is not a simplifying assumption but a lossless reduction guaranteed by the revelation principle.

### 3.6 VCG: Making Truth a Dominant Strategy

Now we return to the question of Section 1: does there exist a mechanism that makes truthful reporting a dominant strategy while giving the allocation to those who value it most? There does, and it is the Vickrey-Clarke-Groves (VCG) mechanism, which this section presents and proves to be dominant strategy truthful.

The construction of VCG takes only two sentences. The allocation rule: choose the allocation that maximizes the sum of all participants' reported values, that is, the efficient allocation. The payment rule: each participant $i$ pays the externality it imposes on everyone else, namely "the maximum total value the others would have obtained without $i$" minus "the total value the others actually obtain under the chosen allocation with $i$." The second part is precisely the so-called Clarke pivot: you pay only when your presence actually changes others' allocation, and you pay exactly the value of the others you crowd out.

::: {.theorem}
(Dominant strategy truthfulness of VCG) Under the VCG mechanism, no matter how the other participants report, reporting one's own value truthfully is a (weakly) dominant strategy for every participant.
:::

The proof is short and pretty, and worth giving in full. Let participant $i$'s true value function be $v_i(\cdot)$, and let the mechanism choose allocation $a^*$ based on all reports. $i$'s payment is

$$t_i = \max_{a}\sum_{j \neq i} \hat v_j(a) \;-\; \sum_{j \neq i} \hat v_j(a^*)$$

where $\hat v_j$ is $j$'s report. $i$'s net utility (measured at its true value) is

$$u_i = v_i(a^*) - t_i = \Big[v_i(a^*) + \sum_{j \neq i}\hat v_j(a^*)\Big] - \max_{a}\sum_{j \neq i}\hat v_j(a)$$

The second term on the right, $\max_a \sum_{j\neq i}\hat v_j(a)$, contains none of $i$'s report and is a constant independent of what $i$ reports. So to maximize $u_i$, $i$ need only maximize the first term in brackets, $v_i(a^*) + \sum_{j\neq i}\hat v_j(a^*)$. But $a^*$ is the allocation the mechanism chooses to maximize $\hat v_i(a) + \sum_{j\neq i}\hat v_j(a)$ given the reports, and if $i$ reports truthfully with $\hat v_i = v_i$, the $a^*$ the mechanism selects maximizes exactly the sum $i$ truly cares about. Misreporting can only make the mechanism select an allocation that is worse, or at best equal, for $i$'s true utility. Hence truthful reporting is a dominant strategy. QED.

The essence of this proof is that "constant independent of one's own report": VCG designs each person's payment to cancel exactly the part of the allocation that concerns the others, so that each person's net utility differs from total social welfare only by a constant beyond their control, and individual self interest aligns perfectly with social efficiency, making truth telling a dominant choice. The second price auction (the Vickrey auction, in which the highest bidder wins and pays the second highest bid) is exactly the single object special case of VCG, and the second highest price the winner pays is precisely the externality it imposes on the bidder with the second highest value. Section 6 applies VCG to Vireo's two ad slots and shows how it yields the figure of 680 yuan, a revenue built on truth.

### 3.7 Bayesian Persuasion: Designing Information Rather Than Allocation

The mechanisms so far all design "who gets what and who pays how much," but the platform holds another, subtler lever: it often controls information, and can decide how to disclose information about quality, matching, and risk to the other side. Designing the information itself, so as to steer the other side's choice, is information design, and its foundational framework is Kamenica and Gentzkow's (2011) Bayesian persuasion. Its key difference from signaling and screening lies in where commitment sits: the sender (say, the platform) first publicly commits to a disclosure rule (an experiment mapping the true state to a signal), then nature draws the state and a signal is released according to the rule, and the receiver (say, the buyer) updates by Bayes' rule upon seeing the signal and then chooses an action in its own interest. The sender cannot lie (the signal distribution is locked by commitment), but can strategically choose "how much to reveal and how to reveal it."

A minimal example of a platform recommendation makes the point best. Suppose a product is in one of two states, good or bad, the buyer's prior is $\Pr(\text{good}) = 0.3$, the buyer buys only when the posterior $\Pr(\text{good}) \geq 0.5$, and the platform only wants to maximize the probability of purchase and does not care about good or bad. When no information is disclosed, the posterior is just the prior, $0.3 < 0.5$, the buyer does not buy, and the platform gets nothing. But if the platform designs a disclosure rule, the effect changes at once: when the product is good it always sends "recommend to buy," and when it is bad it also sends "recommend to buy" with probability $x$. The buyer's posterior after seeing "recommend to buy" is $0.3/(0.3 + 0.7x)$, which the platform pushes down to exactly the threshold $0.5$ (any higher wastes persuasive power), solving to $x = 3/7$. The probability of purchase is then $0.3 + 0.7\times\tfrac{3}{7} = 0.6$, lifting the transaction probability from $0$ to $0.6$ by pure information design. The point is that optimal disclosure is neither "tell everything" nor "tell nothing," but a controlled mixing in of just enough noise to leave the receiver half believing, exactly to the point of being willing to act.

::: {.intuition}
Bayesian persuasion has a geometric solution called concavification, and one picture is worth a thousand words. Draw the sender's payoff as a function of the posterior belief $\mu$, call it $v(\mu)$. In the example above $v(\mu) = 1$ when $\mu \geq 0.5$ (the buyer buys) and $0$ otherwise, a step that jumps up at $0.5$. The best expected payoff the sender can achieve equals the value at the prior of the concave closure of $v$ (the smallest concave function that lies above $v$). On the segment $[0, 0.5]$ this concave hull is the straight line $2\mu$ connecting $(0,0)$ to $(0.5,1)$, and plugging in the prior $0.3$ gives exactly $0.6$, not a hair off the purchase probability computed by hand above. The intuition of concavification is this: any information experiment is nothing but a way to split the prior into several posteriors whose weighted average is still the prior (Bayes plausibility), and what the sender can achieve is precisely the convex combinations of the payoffs at these posteriors, and taking the supremum over all splits gives the concave hull. When information has value is instantly visible: if and only if $v$ is non concave near the prior, so that concavification can raise the value; if $v$ is already concave, no disclosure beats disclosing nothing. This method is the theoretical backbone of ratings, recommendations, certification, and platform information policy, and it bears directly on research into how platforms shape user beliefs in the age of AI.
:::

### 3.8 Cheap Talk: Soft Information Without Commitment

The power of Bayesian persuasion comes from commitment, but a great deal of real world information transmission has no commitment at all: a seller's verbal description of a product, a user review, an analyst recommendation, are all soft information in the sense of cheap talk (Crawford and Sobel 1982), where the sender observes the state and offhandedly sends a costless, unverifiable message that binds it to nothing after the fact, and the receiver acts on it. With no commitment and a message that costs nothing, how much information can be transmitted depends on how aligned the two sides' interests are.

The conclusion is clean and deep. If the sender and receiver have perfectly aligned interests (bias $b = 0$), fully truthful communication is possible; once the sender has a systematic bias (say, a seller always wants you to think quality is higher), full communication collapses, because the sender always has a motive to exaggerate. What can be sustained then is only a partition equilibrium: the state space is cut into several intervals, the sender honestly reveals only which interval the state falls into rather than the exact value, and the receiver believes only down to the interval. The larger the bias $b$, the fewer intervals can be sustained and the coarser the communication, and past some threshold it degenerates to a babbling equilibrium, in which the message carries no information and the receiver simply ignores it. Take the standard setup of a uniformly distributed state and quadratic loss: the largest number of intervals $N$ that can be sustained is the largest integer for which $N(N-1) < 1/(2b)$ holds, and every increment in bias makes possible communication one notch coarser. Setting cheap talk and persuasion side by side is most instructive: both transmit information, but with commitment (persuasion) the sender can actively shape the other side's beliefs to its own advantage, while without commitment (cheap talk) it cannot even manage truthfulness, and its informativeness is pinned down hard by the conflict of interest. How much rating systems, seller descriptions, and influencer recommendations on a platform can be trusted, the answer hides in this divide between commitment and interest alignment.

To summarize: asymmetric information splits into adverse selection and moral hazard, the former addressed by signaling (the informed party burns money to speak up) or screening (the uninformed party lays out a menu to induce selection), at the cost of wasteful money burning or low end distortion plus information rent, and the latter addressed by pay that floats with performance to incentivize unseen effort, at the cost of the risk-incentive tradeoff plus a limited liability rent, with its optimal contract shaped by the likelihood ratio of the performance signal (the informativeness principle); the revelation principle guarantees that we need only study direct truthful mechanisms, and VCG gives a concrete mechanism that makes truth a dominant strategy while achieving efficiency; and when what the platform holds is information itself, Bayesian persuasion uses commitment plus concavification to design disclosure, while cheap talk without commitment is limited by interest alignment, and the two together map out the limits of information as a lever. The next section applies this mechanism language to auctions and matching, answering the other half of the target quantities, revenue and allocation.

## 4 Auctions and Revenue Equivalence

Auctions are the most mature application of mechanism design and the one closest to platform business. This section adds a revenue perspective to the previous section's efficiency perspective, answering two questions: whether different auction rules yield the same revenue, and what the revenue maximizing auction looks like.

### 4.1 Four Basic Auctions and Two Clues About Truth

Single object auctions come in four textbook forms: first-price sealed-bid (the highest bidder wins and pays its own bid), second-price sealed-bid (Vickrey, the highest bidder wins and pays the second highest bid), English (an ascending open outcry until only one bidder remains), and Dutch (a descending price from high until someone accepts). They are pairwise equivalent: Dutch and first-price are strategically equivalent (both require deciding, in one shot and without knowing whether others bid higher, how much to bid), and English and second-price are outcome equivalent (both make the winner effectively pay only up to the second highest value).

The second price auction is the single object special case of the previous section's VCG, and truthful reporting is a dominant strategy: report your true value, and if you win you pay the second highest price and lock in a nonnegative surplus, while misreporting can only lose a win you should have had or turn a nonloss into a loss. The first price is different: paying your own bid means reporting the true value gives zero surplus, so you must shade. How much to shade is a Bayesian game: in the standard setup of symmetric, independent private values, values uniform on $[0, 1]$, and $n$ bidders, the symmetric equilibrium bid function is $b(v) = \frac{n-1}{n} v$, with everyone shading their bid to $\frac{n-1}{n}$ of their true value. The more competitors (the larger $n$), the less they shade, because with more rivals a low bid is easily topped.

### 4.2 Revenue Equivalence

Two auctions, one that cannot extract the truth and forces shading, and one that can extract the truth, will their revenues differ wildly? The revenue equivalence theorem gives a startling answer: under fairly broad conditions, their expected revenues are exactly the same.

::: {.theorem}
(Revenue equivalence theorem) When bidders are risk neutral, private values are independent and identically distributed, and the mechanism awards the object to the bidder with the highest value while giving the lowest type an expected utility of zero, every such auction mechanism yields the seller the same expected revenue.
:::

Verify with the uniform example above. The second price auction's expected revenue is the expectation of the second highest value, and the second highest of $n$ uniform $[0,1]$ values has expectation $\frac{n-1}{n+1}$. The first price auction's expected revenue is the expectation of the highest bid, and writing the highest value as $v_{(1)}$, this is $\frac{n-1}{n}\,\mathbb{E}[v_{(1)}] = \frac{n-1}{n}\cdot\frac{n}{n+1} = \frac{n-1}{n+1}$, not a hair off the second price. When $n = 2$ both are $1/3$, and when $n = 3$ both are $1/2$. The first price relies on "everyone shades but occasionally wins with a high bid," and the second price on "report the true value but pay only the second highest," two different paths whose expected revenues meet at the same endpoint.

The reason behind revenue equivalence is deep. It shows that the seller's expected revenue does not depend on the auction's surface rules but only on the allocation rule (who wins) and the treatment of the lowest type. The proof of the theorem uses precisely the revelation principle: after reducing any auction to a direct truthful mechanism, each bidder's expected payment is uniquely pinned down by its winning probability curve (Myerson's envelope argument), so as long as the winning rule is the same, the expected revenue is the same. This conclusion is extremely useful for platforms designing auctions: since revenue does not come from gimmicks, the choice of auction form should turn on other dimensions, such as collusion resistance, simplicity, and resistance to strategic manipulation, rather than the fantasy that some exotic rule can conjure extra revenue out of thin air.

### 4.3 Optimal Auctions and Reserve Prices

Revenue equivalence says that auctions in the class that "awards the object to the highest value bidder" yield the same revenue, but it does not say this class earns the most revenue. Myerson's optimal auction proves that the revenue maximizing auction usually has to sacrifice a bit of efficiency, by setting a reserve price below which nothing is sold, even if the buyer's value is still positive.

The logic of the reserve price is monopoly pricing incarnate within an auction. The seller facing the highest value buyer is essentially a monopolist: raising the threshold (the reserve price) loses the trade with some probability, but collects more money when the trade happens, and the optimal reserve balances the two. For the uniform $[0, 1]$ distribution, the optimal reserve solves $r = \frac{1 - F(r)}{f(r)}$, that is $r = 1 - r$, giving $r^* = 1/2$. With this reserve set, when $n = 2$ the seller's expected revenue rises from $1/3$ without a reserve to $5/12$, at the cost of a positive probability that both buyers' values fall below $1/2$ and the object goes unsold, damaging efficiency. This is the precise embodiment in auctions of the tension between efficiency and revenue: a revenue seeking seller would rather risk a failed sale in order to use the reserve to squeeze out the money from low value trades. Myerson generalizes this logic into the virtual value, and the optimal auction is equivalent to awarding the object to the bidder with the highest, positive virtual value, with the reserve price being the point where the virtual value equals zero.

### 4.4 Matching and Market Design

Auctions use price to allocate objects, but the platform has a large class of allocation problems that use no monetary pricing at all and instead pair two sides: labor platforms match workers with employers, ride hailing platforms match drivers with riders, admissions systems match students with schools, and dating platforms match users with users. Each side has preferences over the other, and the question is how to match in a way that stays put. The core concept characterizing "matched in a way that stays put" is stability: a matching is stable if there is no blocking pair (a pair on the two sides who both prefer each other to their current partners). Stability is the primary objective because an unstable matching unravels on its own, since a blocking pair always has a motive to pair off privately, and in a decentralized market this coming apart is a real hazard.

Gale and Shapley's (1962) deferred acceptance (DA) algorithm gives a construction that always produces a stable matching. Let one side (say, workers) propose in order of preference, let the other side (employers) tentatively hold the best offer in hand and reject the rest, let the rejected turn to their next choice, and repeat until no new rejection occurs. The algorithm must terminate, and its endpoint must be a stable matching. It also has two important properties. The first is proposer-optimality: the proposing side gets the best partner it can get in any stable matching, and the side being proposed to gets the worst, so who moves first makes a big difference. The second is that strategy-proofness holds for the proposing side: reporting preferences truthfully is a dominant strategy for the proposer, though sadly no stable matching mechanism can be strategy-proof for both sides at once. Section 6.5 runs DA on a 3 by 3 numerical example and shows how proposer-optimality shows up in the numbers. This theory is no blackboard toy: hospital resident matching, public school choice, and the matchmaking of all two sided platforms rest on it at bottom, and Roth's line of market design took DA into these real markets, which is why it won a Nobel Prize.

### 4.5 Myerson-Satterthwaite: The Impossibility of Efficient Bilateral Trade

Mechanism design has good news like VCG, where "efficiency is attainable," and it has bad news that draws a boundary, worth every platform market maker's memory. Consider the simplest matchmaking: a seller privately knows its cost $c$, a buyer privately knows its value $v$, both drawn from overlapping distributions, and efficiency requires trade if and only if $v \geq c$. Myerson and Satterthwaite (1983) prove that no mechanism can simultaneously achieve efficiency (trade whenever it should), incentive compatibility (both sides willing to report truthfully), individual rationality (both sides participate voluntarily), and budget balance (no reliance on outside subsidy, the platform breaking even on its own). The four good properties cannot all be gathered at one table, and some efficient trades that should happen are always forced to be given up.

Take the classic setup with $v$ and $c$ each uniform on $[0, 1]$ to see it most clearly. The first-best should trade when $v \geq c$, with probability $1/2$, while the second best (the linear equilibrium of the buyer seller double auction) trades only when $v \geq c + 1/4$, with the trade probability falling to $9/32 \approx 0.281$, and all those efficient trades whose gains from trade fall between $0$ and $1/4$ are sacrificed. This wedge of $1/4$ is exactly the efficiency loss forced out by the combination of both sides hiding their private information, both retaining the right to walk away, and no outside subsidy allowed. Its significance is fundamental: a platform matchmaking a buyer and a seller, each holding a private value, cannot simultaneously consummate every mutually beneficial trade, keep participation voluntary, and not subsidize the deal out of its own pocket. This impossibility marks the limit of market making, and it explains why real platforms fall back on second best devices like posted prices, commissions, and thickening the market. It stands in exact contrast to VCG: VCG can achieve efficiency but cannot achieve budget balance (it needs an outside subsidy to break even), and Myerson-Satterthwaite says precisely that efficiency, budget balance, and voluntary participation cannot all be had at once.

To summarize: among the four basic auctions the second price extracts the truth and the first price forces shading, but the revenue equivalence theorem proves that as long as the winning rule is the same and bidders are symmetric and risk neutral, the expected revenue is equal, so the choice of auction form should turn on collusion resistance and simplicity rather than the fantasy of collecting more; pursuing maximum revenue requires sacrificing efficiency and setting a reserve, whose logic is the auction version of monopoly pricing; two sided allocation without money is handled by matching, and Gale-Shapley's deferred acceptance always produces a stable matching that is proposer-optimal and strategy-proof for the proposer; and Myerson-Satterthwaite draws a hard boundary of mechanism design, that efficiency, incentive compatibility, voluntary participation, and budget balance cannot all be had at once, so market making has its natural limits. These conclusions all rest on assumptions like independent private values and risk neutrality, and Section 7 discusses what happens when they fail.

### 4.6 From Incentive Constraints to Testable Predictions

The empirical content of information models comes from how the incentive constraints move as the primitives move. Take screening: single crossing makes the high type naturally more willing to choose a high quality contract; and to keep the high type from disguising itself as the low type to grab a bargain, the principal usually pushes down the allocation the low type receives, and the more expensive the information rent it pays, the more severe this downward distortion. In the standard model the high type end satisfies no distortion at the top and is not distorted. So when the type gap widens, or when the share of high types rises, the direction and magnitude of the low type's distortion move together with the objective function, and this can be read in the data through the shape of the contract menu, users' upgrade choices, and the pattern of bunching (different types clustering on the same contract).

For a separating equilibrium to hold in a signaling model, the premise is that sending the signal is less costly for the high type. Once education, certification, or a platform seller badge becomes affordable to everyone, more low types imitate along, the threshold needed for separation is pushed up, and it may even collapse entirely into pooling. Hidden here is a recurring identification trap: seeing that sellers holding a badge have higher quality only shows that selection exists, that good sellers are more willing and more able to obtain the badge, and it does not mean that the badge itself transmits value to buyers. To identify the true effect of signaling, one must exogenously change the visibility of the badge or the receiver's beliefs, rather than simply comparing holders with non holders.

On the moral hazard side, the informativeness principle says that as long as a performance signal carries information about effort even after controlling for other signals, the optimal contract should generally respond to it. The more precise the signal, the lower the risk cost of reaching the same incentive strength, so one can either strengthen the incentive or lower the compensation cost. In the data, after a monitoring technology goes live, the pay-performance sensitivity should change accordingly. But be careful: if this technology, while improving measurability, also directly raises productivity, then one must design a further layer of variation to separate the measurement channel from the productivity channel, or the two get mixed into a single coefficient and cannot be read apart.

Auction theory gives a sharper joint prediction. Under independent private values plus symmetry and risk neutrality, various standard auctions satisfy revenue equivalence and yield equal revenue; the bid shading in the first-price weakens as the number of competitors grows, while the second-price is always truthful bidding. And affiliation or risk aversion breaks this equivalence. So only when the revenue gap between different auction formats varies systematically with the number of bidders, the degree of information disclosure, or common-value exposure does the data truly have the power to distinguish these mechanisms, and looking at the revenue of a single auction is not enough.

| Mechanism | Comparative static | Observable prediction | Identification warning |
|---|---|---|---|
| Screening | Distortion of the low type rises as the information rent gets more expensive | Quality falls at the low end of the menu, bunching rises | Type distribution and cost are often endogenous together |
| Signaling | Signal cost falls for the low type | Pooling rises, signal-quality gradient flattens | Holder quality gap may be only selection |
| Moral hazard | Performance signal becomes more precise | Pay-performance sensitivity or contract cost changes | Monitoring may directly raise productivity |
| Auction competition | Bidder count rises | First-price shading shrinks, revenue rises | Bidder entry is endogenous, needs an entry shifter |
| Persuasion | Sender prior or receiver threshold changes | Disclosure partition and action rate jump | Requires observing the information policy and beliefs |

These predictions share the same discipline: IC, IR, or Bayes plausibility provide only necessary constraints, and the data satisfying these constraints does not mean any one mechanism has been uniquely singled out. To truly identify a mechanism, we rely on the kind of institutional change that makes competing models give different joint responses to the same shock.

## 5 Anchoring Papers

The abstract theorems of mechanism design gain their weight only when they land in real auctions and markets. This section selects two foundational works, one that laid the entire analytical framework for the optimal auction and revenue equivalence, and one that carried mechanism design into search advertising, the most important real world battlefield of platform economics, each organized by paper, method, result, significance, and limitation.

### 5.1 Myerson (1981)

::: {.case}
Paper and position of the method: "Optimal Auction Design," Mathematics of Operations Research. This paper rewrote auction theory in the language of mechanism design, giving both the revenue equivalence theorem and the complete characterization of the optimal auction, and it is the cornerstone of the entire field, with Myerson receiving the 2007 Nobel Prize in Economics for his contributions to mechanism design.

Method: using the revelation principle to reduce any auction to a direct truthful mechanism, using an envelope argument to prove that each bidder's expected payment is uniquely determined by its winning probability (which directly gives revenue equivalence), then introducing the virtual value to write the seller's expected revenue as an expectation over virtual values, thereby turning "designing the revenue maximizing auction" into a pointwise maximization problem, and reaching the conclusion that the optimal auction awards the object to the bidder with the highest, nonnegative virtual value, with the reserve price at the zero point of the virtual value.

Result: it established the revenue equivalence theorem (standard auctions yield the same revenue under broad conditions), and characterized the optimal auction (which generally needs a reserve price, trading efficiency for revenue). These two results, one showing that "gimmicks do not collect more money" and one showing that "collecting more requires a threshold," together lay the theoretical foundation of modern auction design.

Significance and limitation: Myerson's framework is the full source for Section 4 of this chapter, and it is also the theoretical premise for the structural auction econometrics later in this series, since empirical auction analysis recovers the distribution of private values precisely by inverting bidders' equilibrium strategies. The limitation lies in its core assumptions, independent private values, risk neutrality, and single dimensional types, and once buyers' values are correlated (common value), bidders are risk averse, or types are multidimensional, the clean conclusions of revenue equivalence and the optimal auction all loosen, as detailed in Section 7.
:::

### 5.2 Edelman, Ostrovsky and Schwarz (2007)

::: {.case}
Paper: "Internet Advertising and the Generalized Second-Price Auction: Selling Billions of Dollars Worth of Keywords," American Economic Review. It brought the search ad slot auction, a real mechanism worth tens of billions of dollars a year, into rigorous mechanism design analysis, and it is required reading for the information systems and platform economics direction, with Varian (2007) independently giving a closely related analysis.

Method: what search engines actually use is not VCG but the generalized second-price (GSP) auction: advertisers report a per click price, higher reports get the more forward slots with more clicks, and the winner of each slot pays the per click fee of the report ranked just below it. The paper points out that despite the "second price" in its name, GSP is not a generalization of VCG, and truthful reporting is not a dominant strategy in GSP. They characterize the equilibria of GSP using game theory, propose the refinement of locally envy-free equilibrium, and prove that the lowest such equilibrium reproduces exactly the allocation and payments of VCG.

Result: GSP has multiple Nash equilibria, and truthful reporting is generally not one of them; but there exists a locally envy-free equilibrium reproducing VCG's efficient allocation and payments, so at this equilibrium GSP is both efficient and equal to VCG in revenue. This explains why the real world search ad auction uses the simpler form GSP yet can still approximate VCG's theoretical properties.

Significance and limitation: this paper is the direct prototype of this chapter's running case, and it welds mechanism design to the real business of a platform. The limitation lies in its static, single keyword, exogenous click through rate setup, which simplifies reality, since real ad auctions have quality scores, budget constraints, cross keyword competition, and dynamic learning, and a large body of later work proceeds precisely by relaxing these, a thread we will return to later in this series when we cover platforms and advertising.
:::

Taken together the two papers mark out this chapter's theoretical depth and its real world landing point: Myerson gives the general theorems of auction design from the most abstract angle, and Edelman-Ostrovsky-Schwarz map these theorems onto search advertising, a core battlefield of platform economics. They also foreshadow the later theme of structural auction econometrics, inverting the distribution of private values from observed bids.

## 6 Running Case: Vireo's Ad Slot Auction

Now we run the mechanisms of Sections 3 and 4 in full on Vireo's two ad slots. Recall the setup: the top slot has 100 clicks and the second slot 40 clicks, and the three advertisers have per click values 10, 6, and 4. All the numbers can be computed by hand, and the R code below transcribes the hand computation verbatim and checks it.

### 6.1 The Efficient Allocation and VCG Payments

Efficiency requires giving the slots to those who value clicks most: the advertiser with value 10 gets the top slot (100 clicks), the advertiser with value 6 gets the second slot (40 clicks), and the advertiser with value 4 gets nothing.

![The efficient allocation and VCG payments in Vireo's position auction: three advertisers (per click values 10, 6, 4) compete for two slots (click rates 100, 40). The highest value wins the top slot and pays the externality of 520 it imposes on others (5.2 per click); the second highest gets the second slot and pays 160 (4 per click); total revenue 680.](assets/fig/fig_03_position.svg)

VCG payments are computed by externality. For advertiser 1 who gets the top slot, ask "how much would the others get without it": advertiser 2 would rise to the top slot and get $6 \times 100$, advertiser 3 would move up to the second slot and get $4 \times 40$, for a total of 760; and "how much the others actually get with it": advertiser 2 gets $6 \times 40 = 240$ in the second slot, advertiser 3 gets nothing, for a total of 240. So advertiser 1 pays $760 - 240 = 520$, or $520/100 = 5.2$ per click. For advertiser 2, without it advertiser 1 still sits in the top slot and advertiser 3 moves up to the second slot for $4 \times 40 = 160$, so the others total $1000 + 160 = 1160$; with it advertiser 1 gets 1000 in the top slot and advertiser 3 gets nothing, for a total of 1000; so advertiser 2 pays $1160 - 1000 = 160$, or $160/40 = 4$ per click. Advertiser 3 gets nothing and pays 0.

```r
v   <- c(10, 6, 4)            # per click values (sorted)
ctr <- c(100, 40, 0)         # slot click rates (3rd entry means no slot)
# efficient allocation: advertiser i -> slot i; VCG payment = others' max value in its absence - others' realized value
pay1 <- (v[2]*ctr[1] + v[3]*ctr[2]) - (v[2]*ctr[2] + v[3]*ctr[3])
pay2 <- (v[1]*ctr[1] + v[3]*ctr[2]) - (v[1]*ctr[1] + v[3]*ctr[3])
c(pay1 = pay1, pay2 = pay2, per_click_1 = pay1/ctr[1], per_click_2 = pay2/ctr[2],
  VCG_revenue = pay1 + pay2)
#>        pay1        pay2 per_click_1 per_click_2 VCG_revenue
#>         520         160         5.2         4.0         680
```

The total VCG revenue is $520 + 160 = 680$, and the advertisers' net utilities are $10\times100 - 520 = 480$ and $6\times40 - 160 = 80$, both nonnegative, so truthful reporting is a dominant strategy. Compared with the product manager's 1240 from Section 1, built on "the report equals the value, pay your bid," VCG collects less, but every cent of it stands on the advertisers' truth, and that is credible, predictable revenue.

### 6.2 GSP Versus VCG

The real world search engine uses GSP: the winner of each slot pays the per click fee of the report ranked just below it. If advertisers naively report truthfully 10, 6, and 4, the top slot pays the second report 6 and the second slot pays the third report 4, for revenue $6\times100 + 4\times40 = 760$, even higher than VCG. But this is not a stable outcome, because truthful reporting is not a dominant strategy in GSP, and advertisers strategically adjust their reports. Edelman-Ostrovsky-Schwarz prove that GSP has a locally envy-free equilibrium whose allocation and payments reproduce VCG exactly, that is, the per click prices return to 5.2 and 4 and the total revenue returns to 680. So the 760 is a paper figure under "naive truth," and the strategic equilibrium pulls it back to VCG's 680.

```r
# GSP naive truthful reporting: top slot pays the second report, second slot pays the third report
gsp_naive <- v[2]*ctr[1] + v[3]*ctr[2]
# GSP locally envy-free equilibrium reproduces VCG: revenue returns to 680
c(GSP_naive_truthful = gsp_naive, GSP_locally_envy_free = 680, VCG = 680)
#>    GSP_naive_truthful GSP_locally_envy_free                   VCG
#>                   760                   680                   680
```

### 6.3 Revenue Equivalence and the Reserve Price

Let us check Vireo's single slot case against the revenue equivalence framework. Suppose the per click values of two advertisers are independent and uniform on $[0, 1]$. The second price auction collects the second highest of the truthful reports, with expectation $\frac{n-1}{n+1} = 1/3$; the first price auction has each shade to $b(v) = \frac{1}{2}v$, and the expectation of the highest bid is also $1/3$, the two equal, exactly revenue equivalence. If Vireo pursues revenue rather than efficiency and sets the optimal reserve $r^* = 1/2$, the expected revenue rises to $5/12$, at the cost of a $1/4$ probability that both advertisers' values fall below $1/2$ and the slot goes unsold.

```r
n <- 2
spa <- (n - 1)/(n + 1)                 # second price expected revenue = expected second highest value
fpa <- (n - 1)/n * n/(n + 1)           # first price: shading (n-1)/n * expected highest value
rev_reserve <- 5/12                    # uniform[0,1], n=2, expected revenue with reserve 1/2
c(second_price = spa, first_price = fpa, with_reserve = rev_reserve)
#> second_price  first_price with_reserve
#>    0.3333333    0.3333333    0.4166667
```

### 6.4 Information Design: Persuasion at Vireo's Recommendation Slot

Vireo also runs a product recommendation slot, and this time it is not an auction but the design of information. Suppose a product is in one of two states, good or bad, the user's prior is $\Pr(\text{good}) = 0.3$, the user clicks to buy only when the posterior is at least $0.5$, and Vireo only wants to maximize purchases. When no information is disclosed the posterior stays at $0.3$ and the user does not buy. Vireo therefore commits to a disclosure rule: a good product always sends "recommend," and a bad product also sends "recommend" with probability $x$, pushing the posterior after seeing "recommend" down to exactly the threshold $0.5$:

```r
prior <- 0.3; thresh <- 0.5
x <- prior*(1 - thresh)/(thresh*(1 - prior))     # probability a bad product also sends "recommend"
c(false_positive_x = x, P_buy = prior + (1 - prior)*x,
  concavification = 2*prior, no_info_baseline = 0)
#> false_positive_x            P_buy  concavification no_info_baseline
#>        0.4285714        0.6000000        0.6000000        0.0000000
```

The probability that a bad product is mixed into "recommend" is $3/7 \approx 0.429$, and the purchase probability rises from $0$ under no disclosure to $0.6$. This $0.6$ is not a hair off the concave hull value of $2\times 0.3$ given by concavification, and the two paths (solving the signal structure by hand, taking the concave hull geometrically) meet at the same number. The power of information design lies in the controlled mixing in of noise, leaving the user half believing, exactly to the point of being willing to order.

### 6.5 Matching: Stable Pairing in Vireo's Talent Market

Vireo also opens a matchmaking service, pairing 3 workers with 3 companies, each side with its own preferences. Use deferred acceptance to find a stable matching and compare who proposes:

```r
worker_pref <- list(c(1,2,3), c(2,1,3), c(1,2,3))   # workers' preferences over companies
firm_pref   <- list(c(2,1,3), c(1,2,3), c(1,2,3))   # companies' preferences over workers
DA <- function(propP, recvP){                       # proposers propose, receivers hold/reject
  n <- length(propP); nextc <- rep(1, n); held <- rep(NA_integer_, n); free <- 1:n
  while(length(free)){
    p <- free[1]; r <- propP[[p]][nextc[p]]; nextc[p] <- nextc[p] + 1
    if(is.na(held[r])){ held[r] <- p; free <- free[-1] }
    else if(match(p, recvP[[r]]) < match(held[r], recvP[[r]])){
      free <- c(free[-1], held[r]); held[r] <- p }
  }
  held
}
h <- DA(worker_pref, firm_pref); worker_opt <- integer(3)
for(f in 1:3) worker_opt[h[f]] <- f                 # convert to the worker -> company view
firm_opt <- DA(firm_pref, worker_pref)              # companies propose, already worker -> company
rbind(worker_proposes = worker_opt, firm_proposes = firm_opt)
#>                 [,1] [,2] [,3]
#> worker_proposes    1    2    3
#> firm_proposes      2    1    3
```

The two proposing directions give two matchings that are both stable yet different. When workers propose, worker 1 and worker 2 both match their top choice (company 1 and company 2); when companies propose, these two workers match only their second choice (company 2 and company 1), and the two companies get their wish instead. Proposer-optimality is just this plain in the numbers: whoever moves first gets the best partner it can get in any stable matching.

The accounts of this section now balance out: on the same Vireo platform, in the ad auction VCG prices by externality to make truth a dominant strategy and collects 680, GSP is simpler in form and its locally envy-free equilibrium reproduces VCG's 680, and revenue equivalence proves the first and second price have equal expected revenue while pursuing revenue relies on the reserve to lift $1/3$ to $5/12$; putting on the hat of information design, Vireo uses a committed disclosure rule to persuade the purchase probability up to $0.6$, exactly the concave hull value; putting on the hat of matchmaking, deferred acceptance produces a stable matching and proposer-optimality is plain to see. Every step can be computed by hand and self checked, and the distance from abstract theorem to platform business is laid bare in these numerical examples.

## 7 Failure Modes and Robustness of Mechanisms

The theorems of mechanism design are beautiful, but carrying them to real platforms goes astray in a few places, which this section works through one by one.

Though VCG carries the theoretical halo of dominant strategy truthfulness, it is rarely adopted as is in reality, and the reasons are worth spelling out. The first is the non monotonicity and fragility of revenue: VCG revenue can fall rather than rise as competitors increase, can be extremely low under certain allocations, and is easily manipulated by shill bidding (a seller staging fake competitive bids) and bidder collusion. The second is complexity and opacity: VCG's payment rule requires solving for the counterfactual allocation in others' absence, advertisers struggle to understand why they pay a given price, and understandability is crucial to a real market's adoption. This is precisely the practical reason search engines choose the simpler form GSP: though GSP sacrifices dominant strategy truthfulness, it buys transparency and robustness in exchange, and its equilibrium can approximate VCG's properties. This seam between theoretical optimum and engineering usability is a theme that recurs whenever mechanism design lands in practice.

Common value and the winner's curse are a major source of failure for the revenue equivalence assumptions. Section 4's analysis assumes private value (each person's value concerns only itself). But in many platform auctions the value has a common value component: the true commercial value of a keyword is roughly the same for all advertisers, only each has a noisy estimate. Then the winner is often the one who overestimates the most, and winning is itself the bad news that "your estimate was high," so a rational bidder must shade in advance to account for it, and this is the winner's curse. It breaks revenue equivalence, and it also forces auction design to consider how information is revealed during bidding, with English auctions, because they progressively reveal exit information, often earning higher revenue than sealed bid auctions under common value.

Correlated types loosen revenue equivalence further. Revenue equivalence requires types to be independent and identically distributed, and once bidders' values are positively correlated, the seller can exploit the correlation to design a mechanism with higher revenue (Crémer-McLean's extreme result can even extract all the surplus), and standard auctions are no longer equivalent. In real ad auctions values are highly correlated (the value of the same keyword rises and falls with the broad market), so revenue equivalence should be treated as a benchmark rather than an iron law.

Multidimensional types and multidimensional screening are the soft spot of screening theory. The screening in Section 3.3 has a clean solution under a single dimensional type (one usage dimension), but real customers are often heterogeneous along multiple dimensions (usage, time of day, preference for reliability), and multidimensional screening problems generally have no clean characterization, with the optimal menu possibly extremely complex or even involving randomization. Many of the difficulties platforms face when designing pricing menus trace their roots to the multidimensionality of types.

The robustness of moral hazard also has its boundaries. The optimal contract in Section 3.4 is sensitive to the specification of the performance distribution and the agent's risk preferences, and in reality the performance measure can often be gamed, so a strong incentive may induce the distorted behavior of making the metric look good rather than actually doing the work well, and in a multitask setting a strong incentive on the measurable task also crowds out the unmeasurable but equally important tasks. When designing incentives, measurability and gameability are practical constraints that should be thought through before mathematical optimality.

Stringing these together, mechanism design gives platforms a powerful language for designing rules and steering behavior under private information, but its sharpest theorems all rest on assumptions of private values, independent types, single dimensionality, and commitability, and the collusion, correlated values, multidimensional heterogeneity, and metric gaming of real platforms erode these assumptions one by one. Robust mechanism design lies not in applying the optimal formula but in recognizing which assumption you deviate from and where that deviation pushes the conclusion, and often in preferring a second best but transparent, manipulation resistant mechanism. This is the same attitude that runs through this series: the credibility of a method lies not in the sophistication of the technique but in putting the assumptions on the table.

## 8 Further Reading

::: {.readings}
Required reading, in suggested order:

- Mas-Colell, Whinston and Green, Microeconomic Theory, Chapter 16 (adverse selection, signaling, screening), Chapter 17 (principal-agent and moral hazard), and Chapter 27 (mechanism design and auctions). The theoretical skeleton of this chapter comes mainly from these three chapters, and reading them through first sorts out the two broad classes of asymmetric information and the overall framework of mechanism design.
- Vickrey (1961). The founding work of the second price auction and auction theory, the origin of Section 4's clues about truth, short and classic.
- Krishna, Auction Theory. The standard textbook of auction theory, treating revenue equivalence and the optimal auction in more detail than MWG, suitable for deepening Section 4.

Further reading:

- Myerson (1981). The complete framework of the optimal auction and revenue equivalence, the full source for Section 4, focusing on the virtual value and the envelope argument.
- Edelman, Ostrovsky and Schwarz (2007). The mechanism analysis of the search ad GSP auction, required reading for the information systems and platform direction, directly corresponding to this chapter's running case.
- Varian (2007). An independent analysis of position auctions, complementary to Edelman-Ostrovsky-Schwarz, focusing on the equilibrium characterization.
- Spence (1973). The original text of the education signaling model, the foundation of signaling, for understanding how single-crossing makes a signal credible.
- Akerlof (1970). The origin of the lemons market and adverse selection, the source of the market unraveling logic.
- Holmström (1979) and Grossman and Hart (1983). Two cornerstones of moral hazard and optimal contracts, the former giving the informativeness principle and the latter giving a general analysis independent of the first-order approach, designed precisely to work around the possible failure of the first-order approach, the basis of Section 3.4; for the validity conditions of the first-order approach itself (MLRP plus CDFC), see further Rogerson (1985) and Jewitt (1988).
- Kamenica and Gentzkow (2011) and Bergemann and Morris (2019). The foundation of Bayesian persuasion and the unified survey of information design, the full source for Section 3.7, required reading for the AI and platform information policy direction.
- Crawford and Sobel (1982). The origin of cheap talk and the partition equilibrium, the source of Section 3.8's analysis of the credibility of soft information.
- Gale and Shapley (1962) and Roth and Sotomayor (1990). The foundation and synthesis of deferred acceptance and stable matching theory, to be read alongside Roth (1984) on the empirics of resident matching, the source of Section 4.4's market design, required reading for the two sided platform direction.
- Myerson and Satterthwaite (1983). The original text of the impossibility theorem for efficient bilateral trade, the source of Section 4.5's limit of market making.
:::

::: {.apa-refs}
- Akerlof, G. A. (1970). The market for "lemons": Quality uncertainty and the market mechanism. *The Quarterly Journal of Economics, 84*(3), 488-500. https://doi.org/10.2307/1879431
- Bergemann, D., & Morris, S. (2019). Information design: A unified perspective. *Journal of Economic Literature, 57*(1), 44-95. https://doi.org/10.1257/jel.20181489
- Crawford, V. P., & Sobel, J. (1982). Strategic information transmission. *Econometrica, 50*(6), 1431-1451. https://doi.org/10.2307/1913390
- Edelman, B., Ostrovsky, M., & Schwarz, M. (2007). Internet advertising and the generalized second-price auction: Selling billions of dollars worth of keywords. *American Economic Review, 97*(1), 242-259. https://doi.org/10.1257/aer.97.1.242
- Gale, D., & Shapley, L. S. (1962). College admissions and the stability of marriage. *American Mathematical Monthly, 69*(1), 9-15. https://doi.org/10.1080/00029890.1962.11989827
- Grossman, S. J., & Hart, O. D. (1983). An analysis of the principal-agent problem. *Econometrica, 51*(1), 7-45. https://doi.org/10.2307/1912246
- Holmström, B. (1979). Moral hazard and observability. *The Bell Journal of Economics, 10*(1), 74-91. https://doi.org/10.2307/3003320
- Jewitt, I. (1988). Justifying the first-order approach to principal-agent problems. *Econometrica, 56*(5), 1177-1190. https://doi.org/10.2307/1911363
- Kamenica, E., & Gentzkow, M. (2011). Bayesian persuasion. *American Economic Review, 101*(6), 2590-2615. https://doi.org/10.1257/aer.101.6.2590
- Krishna, V. (2009). *Auction theory* (2nd ed.). Academic Press.
- Mas-Colell, A., Whinston, M. D., & Green, J. R. (1995). *Microeconomic theory*. Oxford University Press.
- Myerson, R. B. (1981). Optimal auction design. *Mathematics of Operations Research, 6*(1), 58-73. https://doi.org/10.1287/moor.6.1.58
- Myerson, R. B., & Satterthwaite, M. A. (1983). Efficient mechanisms for bilateral trading. *Journal of Economic Theory, 29*(2), 265-281. https://doi.org/10.1016/0022-0531(83)90048-0
- Rogerson, W. P. (1985). The first-order approach to principal-agent problems. *Econometrica, 53*(6), 1357-1367. https://doi.org/10.2307/1913212
- Roth, A. E. (1984). The evolution of the labor market for medical interns and residents: A case study in game theory. *Journal of Political Economy, 92*(6), 991-1016. https://doi.org/10.1086/261272
- Roth, A. E., & Sotomayor, M. A. O. (1990). *Two-sided matching: A study in game-theoretic modeling and analysis* (Econometric Society Monographs No. 18). Cambridge University Press. https://doi.org/10.1017/CCOL052139015X
- Spence, M. (1973). Job market signaling. *The Quarterly Journal of Economics, 87*(3), 355-374. https://doi.org/10.2307/1882010
- Varian, H. R. (2007). Position auctions. *International Journal of Industrial Organization, 25*(6), 1163-1178. https://doi.org/10.1016/j.ijindorg.2006.10.002
- Vickrey, W. (1961). Counterspeculation, auctions, and competitive sealed tenders. *The Journal of Finance, 16*(1), 8-37. https://doi.org/10.1111/j.1540-6261.1961.tb02789.x
:::
