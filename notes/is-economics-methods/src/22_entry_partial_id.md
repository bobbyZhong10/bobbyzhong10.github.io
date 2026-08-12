---
title: "Entry, Market Structure & Partial Identification"
subtitle: "Entry Games, Multiplicity, and Partial Identification"
seriesline: "Foundations of Information Systems Economics · Chapter 22"
series: "Foundations of Information Systems Economics: Empirical Methods for Information Systems Economics"
sidebartitle: "Foundations of Information Systems Economics"
sidebarsub: "Theory · Methods · Applications"
chapterlabel: "Chapter 22 · Entry, Market Structure & Partial Identification"
---

## Introduction

A city may be able to support one ride-hailing platform but not necessarily two. Whether firm A enters depends on whether firm B will come, and B's own calculation depends on A. In the end the data record only one outcome: perhaps A is present and B is not. But the data do not tell the researcher whether this was the only possible equilibrium, or whether it was one of two symmetric equilibria that happened to be realized. How many firms operate in a market is not a background condition for demand analysis, it is an endogenous outcome of firms sizing one another up.

This is what makes entry models different from single-agent discrete choice. Complete-information games often admit multiple equilibria: one and the same set of market characteristics and profit parameters can support A entering and can also support B entering. If the theory says nothing about how an equilibrium gets selected, the model does not deliver a unique probability for the outcome, and a standard likelihood has nowhere to start. You can of course force in a convenient selection rule and obtain a point estimate, but the precision past the decimal point only dresses up an assumption the researcher never argued for.

This chapter looks at three more honest routes. The first ignores firm identities and models only the number of entrants: even when "who enters" is not unique, "how many enter" may still be unique, and this yields the Bresnahan-Reiss threshold model and Berry's simulation estimator. The second changes the information structure, letting firms hold private information and recovering a unique choice probability in a Bayesian equilibrium. The third accepts that the model predicts only a set of possible outcomes and uses moment inequalities to confine the parameters to an identified set. A simulated case with three Volt brands across six hundred local markets will show how point identification fails, and why an interval is sometimes more informative than an arbitrary point.

Which route you choose depends on how much the data and the institutions actually let the researcher know: sometimes the number of entrants is more robust than their identities, sometimes private information is enough to pin the choice probabilities down, sometimes the evidence can only rule out part of the parameter space. The profit and equilibrium conditions have to be written out clearly, and the gap left by multiple equilibria has to be reported faithfully too. Partial identification is not a fallback after computation fails, it is an accurate statement of the boundary of the evidence when the theory genuinely says no more.

## 1 A Number That Does Not Add Up

::: {.case}
Volt is a fictional ride-hailing platform with three brands, A, B, and C. Across six hundred local markets (small cities) the platform decides, one market at a time, whether each brand enters. A brand's profit from entering a given market depends on market size, on the brand's local operating conditions, and on how many rivals also come in, with more rivals meaning thinner profit. What we observe is, in each market, whether each of the three brands entered or not. What the platform wants to know is: exactly how much does the presence of one rival lower the appeal of entering? This competitive intensity governs which markets the platform should deploy in, and how a merger would reshape the entry pattern.
:::

Start with the laziest thing you could do: treat the three brands' entry decisions as mutually independent events, run a separate probit for each, and use market size and brand conditions to predict whether each brand enters. This amounts to pretending there is no strategic interaction between brands. Using this independent model to predict how many markets would simultaneously hold two or more brands, the answer is $24.4\%$. But in the data only $21.5\%$ of markets actually have two or more brands. The independent model systematically overstates how crowded markets are, because it never sees that brands are in fact avoiding one another: once a brand enters a market, it lowers a rival's willingness to come in. This overstatement is itself evidence that strategic entry deterrence exists, and the independent model cannot capture it.

So we know we have to build entry as a game. But the game brings a problem far deeper than "the independent model is biased." In about $10.3\%$ of these six hundred markets, the entry game has more than one equilibrium. Given exactly the same market conditions and firm conditions, "only A enters" and "only B enters" can both be self-consistent Nash equilibria in which no one wants to unilaterally change their choice, and the theory cannot tell us which one will be realized. This is not a matter of too little data or a poorly tuned model, it is a property of the model itself: what it delivers is a set of outcomes, not a unique outcome. Standard maximum likelihood needs a probability for each outcome, and in these markets the model refuses to provide one.

To insist on a point estimate, you have to supply the model with an assumption it does not contain, namely which of the multiple equilibria gets selected. Assume one selection rule (say, the most profitable brand always wins), and the estimated competitive intensity $\delta$ is one number; switch to an equally defensible rule (say, a fixed priority ordering across brands), and it comes out as another. On this chapter's data these two point estimates are $1.2$ and $1.3$. Both are compatible with the observed entry pattern, and the data cannot choose between them, because the difference arises only in those multiple-equilibrium markets, where the selection mechanism is simply never observed. The point estimate is therefore not pinned down by the data, it is pushed around by that arbitrary selection assumption.

What this chapter is about is precisely how to treat this situation honestly. Section 3 explains thoroughly why multiplicity breaks point identification, and also makes clear a life-saving fact: even when who enters is not unique, the number of entrants often is. Section 4 gives three ways out, of which the most principled, partial identification, does not fabricate a point but confines $\delta$ to an interval. On this chapter's data that interval is $[1.1, 1.4]$, it contains the true value $1.1$, and its width honestly reflects that data plus model can only bracket this parameter rather than pin it down. An interval estimate looks less tidy than a point estimate, but it is far more credible than a point that is being pushed around by an arbitrary assumption.

## 2 The Economic Model and the Estimand

Before doing anything, let us write out the entry game and the target quantity clearly.

The core target quantity this chapter delivers is the competitive intensity $\delta$, that is, how much one rival's entry lowers other firms' profit. It measures the source of market power: a large $\delta$ says that each additional rival sharply intensifies competition, a small $\delta$ says the market is close to contestable and even many entrants do not depress profit much. Around it there are also entry costs, market-size effects, and other parameters. Together these quantities are what let us answer which markets the platform should enter, and how a merger or a policy would reshape the entry pattern. Unlike in earlier chapters, the target here need not be point-identified: under partial identification it is a set rather than a point, and that fact is itself part of what we deliver.

Now the entry game. Index markets by $m$ and potential entrants by $f \in \{1, \dots, F\}$, with each firm choosing to enter, $a_{fm} = 1$, or stay out, $a_{fm} = 0$. If firm $f$ enters, its profit is

$$\Pi_{fm} = X_m \beta + Z_{fm} \alpha - \delta\, N_{-f,m} + \xi_m + \varepsilon_{fm}$$

Here $X_m$ is market characteristics (such as size), common to all firms; $Z_{fm}$ is a firm-market characteristic (such as the brand's local operating advantage) that varies by firm and enters firm $f$'s own calculation; $N_{-f,m}$ is the number of firms entering other than $f$, with $\delta > 0$ the competitive intensity, so more rivals means lower profit; $\xi_m$ is a market-level unobservable common to all firms in the market (invisible to the researcher); and $\varepsilon_{fm}$ is a firm-market-specific shock. A firm enters if and only if entry profit is nonnegative, $a_{fm} = \mathbf{1}\{\Pi_{fm} \geq 0\}$.

In this entry rule $N_{-f,m}$ depends on other firms' choices, and their choices in turn depend on $N$, so this is a game and not a single-agent decision. A pure-strategy Nash equilibrium is an entry profile $(a_{1m}, \dots, a_{Fm})$ satisfying two conditions: every entrant has nonnegative profit given the other entrants, and every non-entrant would earn negative profit if it entered. Under complete information all firms see all the $\varepsilon$ and $\xi$, only the researcher does not, and the equilibrium is a function of these shocks.

::: {.intuition}
The clearest way to see this is to draw the two-firm entry game in the $(\varepsilon_1, \varepsilon_2)$ plane. If both firms' shocks are very low, either one loses money by entering, and the unique equilibrium is that neither enters. If both are very high, either one profits by entering, and the unique equilibrium is that both enter. The real trouble is the middle band: both firms' shocks lie in the range where "if I were the only one it would pay to enter, but if we both enter we both lose." Here both "I enter, you stay out" and "you enter, I stay out" are equilibria, each is stable, and neither firm wants to change unilaterally. Game theory stops here. It tells you both outcomes are possible, but not which one will happen. This region is the source of multiplicity, and every difficulty and technique later in the chapter grows out of it.
:::

![The equilibrium partition of the two-firm entry game in the $(\varepsilon_1, \varepsilon_2)$ plane. When both shocks are low the unique equilibrium is that neither enters, when both are high it is that both enter; in the red region in the middle both "(1,0), only firm 1 enters" and "(0,1), only firm 2 enters" are Nash equilibria, the theory does not specify which occurs, and this is multiplicity.](assets/fig/fig_22_multiplicity.svg)

To summarize: entry is a game among firms, the estimand is the competitive intensity $\delta$ and related parameters; under complete information the pure-strategy Nash equilibrium may not be unique, and in the middle region who enters is undetermined; this multiplicity is the root of every identification difficulty in the chapter, and the starting point for the next section.

## 3 Identification

Identification comes logically before estimation. This section answers: what can the data from an entry game actually identify, how does multiplicity break point identification, and what structure can rescue identification. It is the most detailed analysis in the chapter, broken into numbered subsections, each of which gives the intuition before the formalism.

### 3.1 Why Multiplicity Breaks Point Identification

Standard maximum likelihood needs one thing: given the parameters, the model assigns a unique probability to each possible outcome, and the product of the probabilities of the observed outcomes is the likelihood. In multiple-equilibrium markets the entry game cannot supply this unique probability. Return to the middle region of Section 2, where both $(1,0)$ and $(0,1)$ are equilibria. Given parameters $\theta$, you can compute the probability that the shocks fall in this region, but whether that probability should be attributed to $(1,0)$ or to $(0,1)$ the model does not say. To write down the likelihood, you would need "the conditional probability that the outcome is $(1,0)$ when the shocks land in the multiple region," and this selection probability is not part of the model, it depends on an equilibrium-selection mechanism outside the model.

Formally, the entry game is an incomplete economic model: it restricts which outcomes can occur but does not map the shocks uniquely onto outcomes. Without that map the likelihood cannot be defined. This is a different kind of difficulty from the ones in earlier chapters. There the model always delivered a unique likelihood or moment condition, and the trouble was that the identifying assumptions (exogeneity, parallel trends) were untestable; here the trouble comes earlier, the model cannot even deliver a single prediction. Forcing in a selection mechanism (assuming some equilibrium is always chosen) makes the likelihood exist, but the selection mechanism itself is unobserved, get it wrong and the estimate is biased, and the data cannot tell you whether you got it right. This is exactly where the two point estimates $1.2$ and $1.3$ of Section 1 come from.

### 3.2 A Life-Saving Fact: The Number of Entrants Is Often Unique

Multiplicity seems to make entry models unestimable, but one key fact rescues half the situation: even when who enters is not unique, the number of entrants often is. Back in the middle region, $(1,0)$ and $(0,1)$ are two different equilibria, but the number of entrants in each is 1. In fact, as long as the competitive effect $\delta$ is the same for all firms (the amount by which one rival lowers profit does not depend on who that rival is), you can rank firms by profitability and let the most profitable enter first, until the last one that can still make money, and this construction yields a unique equilibrium number $N^*$, with multiplicity showing up only in which specific firms occupy those slots. In this chapter's data, uniqueness of the number of entrants holds in a full one hundred percent of the six hundred markets, and in the $10.3\%$ with multiple equilibria the difference between equilibria is entirely in identity and never in number.

This fact opens the first way out. If you model only the number of entrants $N_m$ and do not worry about who enters, multiplicity is no longer an obstacle, because the number has a unique prediction, the likelihood can be written down, and the parameters are point-identified. The price is discarding the identity information, so you cannot answer questions like "which brand is more competitive." Bresnahan-Reiss and Berry take exactly this route, detailed in Section 4.

### 3.3 Recovering a Unique Equilibrium with Incomplete Information

The second way out is to change the information structure. The earlier setup assumed complete information, with every firm seeing everyone's shocks, and multiplicity followed from that. If instead each firm knows only its own private shock $\varepsilon_f$ and cannot see its rivals', the entry game becomes an incomplete-information Bayesian game. Now each firm cannot be certain whether a rival will enter, it can only form a belief about the rival's entry probability and decide by expected profit. An equilibrium is a set of self-consistent beliefs: each firm, expecting others to enter with their equilibrium probabilities, has an optimal entry probability that exactly equals what others expect of it.

This dissolves the multiplicity problem. In a symmetric Bayesian Nash equilibrium, each firm's entry probability $P^*$ satisfies a scalar fixed-point equation

$$P^* = \Lambda\big(\pi(P^*, X)\big) = \frac{\exp(\pi(P^*, X))}{1 + \exp(\pi(P^*, X))}$$

where $\pi(P^*, X)$ is expected profit given that a rival enters with probability $P^*$. Under mild conditions this fixed point has a unique solution, so the equilibrium is unique and point identification is recovered. Seim (2006) generalizes it to entry with location choice across multiple positions, and is the representative of this route. Note that incomplete information is not an approximation to complete information, it is a different assumption, saying that firms genuinely do not know their rivals' private information when they decide. Which information structure fits reality better depends on the setting: retail chains with fixed locations and public decisions may be close to complete information, while fast-moving, mutually secretive platform expansion may look more like incomplete information. This is a modeling choice that must be argued, not a technical detail you can switch at will.

### 3.4 Partial Identification: Accept the Set, Bracket with Moment Inequalities

The third way out is the most principled. It neither dodges multiplicity nor changes the information structure, it accepts the fact that "the model gives only a set of outcomes" and, on that basis, confines the parameters to a set. This is partial identification.

The core idea is to translate the ambiguity of multiplicity into two computable bounds. For any outcome $y$ (say, "only A enters"), given parameters $\theta$ and market characteristics, you can compute two probabilities. One is the probability that $y$ is the unique equilibrium, call it $H_1$, which is a lower bound on the probability that $y$ occurs, because whenever the shocks make $y$ unique, $y$ must occur. The other is the probability that $y$ is some equilibrium (possibly one of several), call it $H_2$, which is an upper bound, because $y$ can occur only when it is at least an equilibrium. The true probability that $y$ occurs is squeezed between the two:

$$H_1(\theta, X) \leq P(y \mid X) \leq H_2(\theta, X)$$

For outcomes that are unique equilibria (such as both enter or neither enters), $H_1 = H_2$, the bounds collapse to a point, and these outcomes are point-identified; for outcomes involved in the multiple-equilibrium region (such as exactly one firm entering), $H_1 < H_2$, there is a genuine gap between the bounds, and it is exactly this gap that houses the unknown selection mechanism. The identified set is the set of parameters that make all these inequalities hold, $\Theta_I = \{\theta : \forall y,\ H_1(\theta,X) \leq P(y \mid X) \leq H_2(\theta,X)\}$. The data (the $P(y\mid X)$ on the left) plus the model (the two bounds) confine $\theta$ to this set, but in general they cannot confine it to a point.

::: {.intuition}
The philosophy of partial identification is worth pausing over. The default in traditional econometrics is point identification: data plus assumptions ought to pin the parameter down to a number. Partial identification concedes that sometimes data plus the (weaker) assumptions you are willing to maintain are only enough to bracket the parameter within a range, and insisting on a point requires adding stronger, often unfounded assumptions (such as an equilibrium-selection mechanism). Rather than sneak in an unbelievable assumption for the sake of a pretty point estimate, better to honestly report the range. The width of the identified set is itself informative: it tells you how much the data and the credible assumptions can jointly say, and the remaining ambiguity is honest ignorance, not noise that more data could remove. This attitude, "better an honest interval than a suspect point," is partial identification's lesson to empirical work.
:::

### 3.5 Inference Under Partial Identification

Under point identification, a confidence interval covers that true parameter point. Under partial identification the situation is subtle, and you have to distinguish two targets. One is to construct a confidence set that covers the entire identified set $\Theta_I$, the other is to cover the true parameter $\theta_0$ inside the identified set. The latter is generally smaller, because it need only cover one point rather than the whole set. This distinction was flagged by Imbens and Manski (2004), and Chernozhukov, Hong, and Tamer, as well as Andrews and Soares, among others, gave systematic inference methods. In practice you must be clear about which confidence set you report, because their interpretation and their size both differ. In addition, the identified set itself has to be estimated by simulation (the two bounds are not analytic functions of the parameters), so finite-sample and simulation error both enter the inference, and these are the technical burdens of putting partial identification into practice.

The key points of this section can be summed up as follows: multiplicity makes the entry game an incomplete model and breaks point identification; the number of entrants, however, is often unique, which opens the first way out of "modeling only the number"; incomplete information turns the game into a Bayesian game with a unique equilibrium, the second way out; partial identification accepts the outcome set and, with the moment inequalities $H_1 \leq P \leq H_2$, confines the parameters to an identified set, the most principled third way out; and its inference must distinguish covering the set from covering the truth. None of this has yet said how to compute anything, and that is the next section.

## 4 Estimation

This section brings the three ways out of Section 3 down to operational estimators: the number models of Bresnahan-Reiss and Berry, Seim's incomplete information, and Ciliberto-Tamer's partial identification.

### 4.1 Bresnahan-Reiss: An Ordered Probit on the Number

The first way out uses only the number of entrants. Bresnahan and Reiss's insight is that market size determines how many firms a market can support, so the observed number of firms is an ordered response to market size. Run an ordered probit of each market's firm count $N_m$ on market size: as size grows, the market successively crosses the thresholds to support 1, 2, 3 firms, and the thresholds correspond to a set of cutpoints on the size index. Estimate the cutpoints and you know how large a market each additional firm requires, and how competition changes with the number of firms.

The key structural quantity is each firm's entry threshold $s_N$, the per-capita market size needed to support $N$ firms. If $s_N$ rises with $N$, each additional firm intensifies competition and later entrants need a proportionally larger market to survive; if $s_N$ is roughly flat, entry no longer intensifies competition and the market is close to contestable. The ratio of adjacent thresholds $s_{N+1}/s_N$ is the measure of competitive intensity. On the Volt data, the ordered probit's market-size coefficient is $1.419$ (larger size means more firms), the estimated cutpoints are $-0.633$, $1.341$, $3.319$, monotonically increasing. The third firm's per-capita threshold relative to the second is $s_3/s_2 = 1.65$, clearly greater than 1, meaning each additional brand really does intensify competition, and the third entrant needs a per-capita market sixty-five percent larger than the second to hold up. This point-identifies the competitive intensity, at the cost of using none of the information about who enters.

### 4.2 Berry 1992: The Number Under Heterogeneous Firms

Bresnahan-Reiss assumes firms are homogeneous. Berry (1992) handles heterogeneous firms while still modeling only the number, relying precisely on the "number is unique" fact of Section 3.2. He studies airline entry on city-pair routes, where firms differ in entry cost, and writes entry value as

$$EV_{fm}(N) = X_m \beta - \delta \ln(N) + Z_{fm} \alpha - \varepsilon_{fm}$$

where $\varepsilon_{fm}$ contains a market-common component (creating within-market correlation) and a firm-specific component. Heterogeneity enters only the cost and not the continuation value that varies with the number of firms, and this restriction guarantees a unique equilibrium number. The boundary of the shock set that generates an $N$-firm equilibrium is complex and cannot be integrated analytically, so Berry uses simulation: draw several sets of shocks, and for each set and each parameter use the "enter in order of profitability" construction to compute the unique equilibrium number, then average over the draws to get $\Pr(N)$, and estimate by simulated moments or simulated likelihood. Berry's framework brings heterogeneous firms in while sidestepping multiplicity by modeling only the number, and is the culmination of this route.

### 4.3 Seim: A Fixed Point Under Incomplete Information

The second way out is Seim's (2006) incomplete-information model. Firms know only their own private shock, and a rival's entry is probabilistic. Estimation proceeds in two steps. First, given parameters, solve the belief fixed-point equation of Section 3.3 to get the equilibrium entry probability (a vector of entry probabilities, one per location, in the multi-location case), iterating by successive approximation to convergence. Second, fit the equilibrium-predicted entry probabilities to the observed entry frequencies and estimate the parameters by maximum likelihood. The great advantage of incomplete information is that the equilibrium is unique and point-identified, and the computation is far lighter than Berry's complete-information model, because the fixed point is low-dimensional. The price is maintaining the assumption that "firms genuinely do not know their rivals' private information when they decide," which suits some settings and not others.

### 4.4 Ciliberto-Tamer: Partial Identification and Moment Inequalities

The third way out is Ciliberto and Tamer's (2009) partial identification, which keeps complete information and firm heterogeneity and does not dodge multiplicity. The core of estimation is the two bounds of Section 3.4. For each outcome $y$ and each candidate parameter $\theta$, use simulation to compute $H_1$ (the probability that $y$ is the unique equilibrium) and $H_2$ (the probability that $y$ is some equilibrium): draw several sets of shocks, for each set enumerate all Nash equilibria, count $y$ toward $H_1$ if it is the unique equilibrium and toward $H_2$ if it is one of the equilibria, and average over draws. The identified set is the set of parameters that make the observed frequency fall inside $[H_1, H_2]$, characterized by a criterion function that penalizes violations:

$$Q(\theta) = \sum_y \Big[ \big(P(y) - H_1(y, \theta)\big)_-^2 + \big(P(y) - H_2(y, \theta)\big)_+^2 \Big]$$

where $(\cdot)_-$ penalizes an observed frequency below the lower bound and $(\cdot)_+$ penalizes one above the upper bound, and the parameters with $Q(\theta) = 0$ satisfy all the inequalities and constitute the identified set. On the Volt data, sweeping $\delta$ over a grid with the other parameters fixed at their true values, the criterion $Q(\delta)$ rises sharply as $\delta$ moves away from the truth, and hugs the floor (the inequalities are satisfiable) over an interval. That interval is $[1.1, 1.4]$, and it contains the true value $1.1$. This is the product of partial identification: not a point but an honest interval, whose width reflects the limit of what data plus model can say without assuming a selection mechanism. By contrast with the two point estimates that drifted with the selection rule in Section 1 ($1.2$ and $1.3$), the partial-identification interval neither needs that arbitrary assumption nor hides the range of uncertainty.

![The Ciliberto-Tamer criterion function $Q(\delta)$ as a function of the competitive intensity $\delta$. As $\delta$ moves away from the truth the criterion rises sharply, and over an interval (the green band, $[1.1, 1.4]$) it hugs the floor, where the moment inequalities are satisfiable, and this interval is the identified set, containing the true value $1.1$ (the red dashed line). The identified set is an interval rather than a point, which is exactly the signature of partial identification.](assets/fig/fig_22_identifiedset.svg)

### 4.5 Trade-offs Among the Three Routes

Each route has its price. Modeling only the number (Bresnahan-Reiss, Berry) is point-identified and light to compute, but discards the identity information and cannot answer "which brand is more competitive." Incomplete information (Seim) is point-identified and light to compute, but requires maintaining a specific information-structure assumption. Partial identification (Ciliberto-Tamer) makes the weakest and most honest assumptions, keeps complete information and heterogeneity, but its product is an interval rather than a point, is the heaviest to compute (enumerating equilibria, simulating the two bounds), and the interval may be wide enough to carry little information. Which route to pick depends on the question: if you care only about the total competitive intensity and the data are counts, use a number model; if you care about firms' heterogeneous competitive effects and are willing to accept an information assumption, use Seim; if you want the fewest assumptions and can tolerate an interval conclusion, use Ciliberto-Tamer. They do not replace one another, they are options placed at different points on the trade-off between strength of assumptions and precision of conclusions.

The route of this section is now complete: Bresnahan-Reiss runs an ordered probit on the number and measures competitive intensity by the threshold ratio; Berry uses simulation to bring heterogeneous firms into a number model; Seim uses an incomplete-information fixed point to recover a unique equilibrium and point identification; Ciliberto-Tamer uses moment inequalities to confine the parameters to an identified set; and the three each occupy a different slot between assumptions and precision.

## 5 Anchor Papers

Methods stand up only when they land in real research. Three anchor papers, one that erects the number-threshold model, one that brings heterogeneous firms into a number model, and one that founds partial identification for entry games, each laid out along the five elements of paper, method, data, results, and limitations, with attention to how the identifying assumptions are defended.

### 5.1 Bresnahan and Reiss (1991)

::: {.case}
Paper: "Entry and Competition in Concentrated Markets," Journal of Political Economy. It uses the number of entrants rather than price data to infer competitive intensity, founding an entire methodology for reading competitive behavior off market structure.

Method: model the firm count in isolated small markets as an ordered response to market size. Variable profit scales with market size while fixed cost does not, and the equilibrium number is determined by size crossing a series of thresholds. Estimate each firm's entry threshold (the per-capita market size needed to support $N$ firms) and use the ratio of adjacent thresholds to measure how competition changes with the number of firms, with no price or quantity data needed.

Data: a number of U.S. small towns isolated from the outside, covering five retail and professional service industries (doctors, dentists, druggists, plumbers, tire dealers), each industry estimated separately. Isolated markets are used to ensure each market is an independent game, free of spillovers from neighboring markets.

Results: competitive intensity varies widely across industries. In most industries the threshold ratio from 1 to 2 firms is clearly greater than 1 (about 1.98 for doctors, 1.78 for dentists, 1.81 for tire dealers), meaning the second firm sharply intensifies competition; but beyond that the threshold ratio quickly falls back toward 1 (for doctors, $s_3/s_2, s_4/s_3, s_5/s_4$ are 1.10, 1.00, 0.95 in turn), so from the third to the fifth firm competition no longer intensifies with entry. Plumbers are the other extreme, with the threshold ratio hugging 1 throughout (only about 1.06 from 1 to 2 firms), and prices barely falling with the number of firms. That the number of entrants alone can read out these competitive differences is the power of this paper.

Limitations: using only the number discards the information about who enters, and cannot tell whether the homogeneous-firm assumption holds; the threshold ratio as a competitive measure depends on the shape of the profit function, and without shape restrictions it is not nonparametrically identified. These limitations are exactly what Berry's introduction of heterogeneity and Ciliberto-Tamer's introduction of identity information later address.
:::

### 5.2 Berry (1992)

::: {.case}
Paper: "Estimation of a Model of Entry in the Airline Industry," Econometrica. It brings heterogeneous firms into an entry model while sidestepping multiplicity by modeling only the number, a pivotal step in the methodological history of entry models.

Method: a complete-information entry game with firms heterogeneous in entry cost. The key construction is that as long as the competitive effect is the same for all firms and heterogeneity enters only the entry cost, you can rank firms by profitability and let them enter in order, obtaining a unique equilibrium number even when who enters is not unique. The boundary of the shock region corresponding to a given equilibrium number is complex, so it is integrated by simulation and estimated by simulated moments.

Data: entry on city-pair routes after airline deregulation, with firm heterogeneity captured by each carrier's airport presence and the like. City-pair markets provide a large amount of cross-sectional variation.

Results: airport presence (a carrier's existing network at the endpoint airports) is a strong predictor of entry, and the heterogeneous-firm model fits the data better than the homogeneous one. The paper demonstrates how to sidestep multiplicity using the structural fact that "the number is unique" while retaining firm heterogeneity.

Limitations: it still uses only the number, and the information about who enters goes unused; uniqueness of the equilibrium number depends on the restriction that the competitive effect is homogeneous, and once you allow firm-specific competitive effects the number is no longer unique either, which then calls for Ciliberto-Tamer's partial identification. This paper is therefore both the peak of the number model and a marker of its boundary.
:::

### 5.3 Ciliberto and Tamer (2009)

::: {.case}
Paper: "Market Structure and Multiple Equilibria in Airline Markets," Econometrica. It founds partial identification for entry games, embracing multiplicity head-on and using moment inequalities to confine the parameters to an identified set, and is the direct source of Sections 3.4 and 4.4 of this chapter.

Method: a complete-information entry game that allows firm-specific competitive effects (the amount by which one rival lowers profit can depend on who it is), so that the equilibrium number is no longer unique either. Rather than assume any equilibrium-selection mechanism, compute two bounds for each outcome, the probability that it is the unique equilibrium (lower bound) and the probability that it is some equilibrium (upper bound), and require the observed frequency to fall between the two. The identified set is the set of parameters satisfying all these inequalities, estimated by a criterion function that penalizes violations together with simulation.

Data: U.S. domestic airline airport-pair markets (about twenty-seven hundred markets), the four major carriers plus mid-size and low-cost carriers, with airport presence and cost variables, and no price data.

Results: once firm-specific competitive effects are allowed, low-cost carriers (such as Southwest) exert markedly stronger competitive pressure on rivals than other carriers, and the role of airport presence is stronger too, things a homogeneous-competitive-effect model cannot see. The paper reports confidence intervals covering the identified set rather than point estimates, honestly reflecting the ambiguity that multiplicity brings.

Limitations: the identified set can be wide, and where there is too little variation the conclusions will be underpowered; enumerating equilibria and simulating the two bounds is computationally heavy; and inference under partial identification (cover the set or cover the truth) requires dedicated methods. The value of this paper lies not in producing more precise numbers but in demonstrating how, when the model can only bracket the parameter, to report that bracketing honestly and rigorously.
:::

Put together, the point of the three anchors is clear: Bresnahan-Reiss erects a method for inferring competition from market structure using the number, Berry brings heterogeneous firms in without losing uniqueness of the number, and Ciliberto-Tamer embraces multiplicity head-on and puts the identity information to use through partial identification. Progress in entry models has always turned on one main thread: how to identify competitive intensity honestly in the face of multiplicity.

## 6 A Full Walk-Through on the Volt Data

Now let us run the tools of Section 4 end to end on Volt's entry data. The code below uses R 4.5.3, fixes set.seed(18) for reproducibility, and every number cited in the text comes from the actual output of running this code.

### 6.1 The DGP

The design parameters are as follows: six hundred local markets, three brands A, B, and C. Brand $f$'s entry profit in market $m$ is $\Pi_{fm} = b_0 + b_S \cdot \text{size}_m + a_Z \cdot Z_{fm} - \delta \cdot N_{-f,m} + \xi_m + \varepsilon_{fm}$, and it enters if and only if profit is nonnegative. The true competitive intensity is $\delta = 1.1$. The equilibrium is solved by the rule "enter in order of standalone profitability," which is a pure-strategy Nash equilibrium.

```r
set.seed(18)
delta <- 1.1                                   # competitive intensity (core parameter to identify)
solve_entry <- function(size, Z, eps, xi){
  base <- b0 + bS*size + a_Z*Z + xi + eps       # profit if entering alone
  ord  <- order(base, decreasing = TRUE)        # most profitable enters first
  entered <- rep(FALSE, 3); nin <- 0
  for (f in ord) if (base[f] - delta*nin >= 0){ entered[f] <- TRUE; nin <- nin + 1 }
  entered
}
```

In the generated data, the distribution of the number of entering brands $N$ is 222 markets with 0, 249 with 1, 112 with 2, and 17 with 3; $63.0\%$ of markets have at least one brand entering; the three brands' entry rates are $0.283$, $0.318$, $0.272$; and the correlation between market size and the number of entrants is $0.75$.

### 6.2 The Independent Probit: Reproducing the Opening Overstatement

```r
pr <- glm(enter ~ size + Z, data = long, family = binomial("probit"))  # no strategic term
# under independence P(N>=2) = 1 - P(0) - P(1), computed from each brand's marginal entry probability
```

Treating the three brands as independent events, and using market size and brand conditions to predict entry, the independent model gives $P(N \geq 2) = 0.244$, while in the data it is actually $0.215$. The independent model overstates how crowded markets are, because it cannot see brands deterring one another. The raw correlation between A's and B's entry is only $0.047$, near zero, and this is exactly where deterrence hides: market size pushes the two up together, deterrence pulls them apart, and the two forces nearly cancel in a simple correlation, so only by building out the game can you see it clearly.

### 6.3 Bresnahan-Reiss: An Ordered Probit on the Number

```r
# ordered probit: N ~ market size, cutpoints on the size index separate 0/1/2/3 firms
fit_br <- optim(start, ord_probit_ll, y = Nentry, X = cbind(size))
```

The market-size coefficient is $1.419$, the cutpoints are $-0.633$, $1.341$, $3.319$, monotonically increasing. Converted into a per-firm per-capita entry threshold, the third firm's threshold relative to the second is $s_3/s_2 = 1.65$, clearly greater than 1, meaning each additional brand intensifies competition, and the third entrant needs a per-capita market about sixty-five percent larger than the second to hold up. This point-identifies the competitive intensity, but uses none of the information about who enters.

### 6.4 Multiplicity and Partial Identification

```r
# enumerate the 8 entry profiles, decide which are Nash equilibria; count the equilibria per market
nash_set <- function(base, delta) which(apply(subsets, 1, function(e) is_nash(e, base, delta)))
```

Counting the number of Nash equilibria per market: $538$ markets have a unique equilibrium, $46$ have two, and $16$ have three, that is, $10.3\%$ of markets have multiple equilibria. Crucially, in these multiple-equilibrium markets the number of entrants is exactly the same across equilibria, the difference being entirely in identity, and uniqueness of the number of entrants is a full one hundred percent, which is exactly why Bresnahan-Reiss and Berry can point-identify using the number alone.

Multiplicity makes the point estimate depend on the selection rule. Assuming "the most profitable enters first" (the correct rule), simulated maximum likelihood gives $\hat\delta = 1.2$; switching to "a fixed priority ordering across brands" (an equally arbitrary rule) gives $\hat\delta = 1.3$. Both are compatible with the data, and the data cannot choose. Partial identification instead assumes no selection rule and confines $\delta$ with Ciliberto-Tamer's moment inequalities:

```r
# for each delta, simulate H1 = P(y is the unique equilibrium), H2 = P(y is some equilibrium)
# identified set = the delta for which the observed frequency falls inside [H1, H2]
Q <- sapply(grid, ct_Q)                        # criterion penalizing violations of the two bounds
```

The criterion $Q(\delta)$ rises sharply as $\delta$ moves away from the truth, and hugs the floor over $[1.1, 1.4]$. This interval is the identified set, and it contains the true value $1.1$. It is not a point, and its width honestly reflects the limit of what data plus model can say without assuming a selection mechanism.

### 6.5 The Full Reconciliation

The walk-through of this section can be summarized as follows: the independent model overstates crowding ($P(N\geq2)$ estimated at $0.244$ against an actual $0.215$), exposing strategic deterrence and forcing a game model; but the game has multiple equilibria in $10.3\%$ of markets, with the number of entrants unique (a full one hundred percent) and identity not unique; Bresnahan-Reiss point-identifies competitive intensity from the number (threshold ratio $s_3/s_2 = 1.65$); and to use the identity information, multiplicity makes the point estimate drift with the selection rule ($1.2$ against $1.3$), while Ciliberto-Tamer's partial identification gives the identified set $[1.1, 1.4]$, an interval that contains the true value and whose width is honest, rather than a suspect point.

## 7 Failure Modes and Robustness

In the simulation the game structure is constructed by design, but in real research all sorts of assumptions can fail at any moment. This section lays out the most common failure modes and the operational responses.

The equilibrium-selection assumption is the most insidious one. Assuming some equilibrium is always chosen for the sake of point identification looks harmless, but in fact it stakes the conclusion on a mechanism that is never observed and cannot be tested. The drift between $1.2$ and $1.3$ in Section 6 is the warning: different selection rules give different competitive intensities, and the data cannot adjudicate. The response is either to admit frankly that the point estimate depends on this assumption and run a sensitivity analysis, or to switch to partial identification, which does not need it. Treating the selection rule as an inconsequential technical setting is as baseless as defaulting to some conduct in Chapter 21.

The information-structure assumption must be defended head-on too. Complete and incomplete information are not two technical options to be switched at will, they make substantively different claims about what firms know when they decide, and they yield different equilibrium structures (the former possibly multiple, the latter often unique) and different estimators. Choose the wrong information structure and the estimate is biased, and the data often cannot distinguish the two. The response is to argue from institutional detail which fits better: whether decisions are public, whether rivals' information is available, whether entry is reversible. This choice must be written out in the open and not hidden in a default setting.

Partial identification's identified set may be wide enough to carry no information, and this is its most practical soft spot. If there is too little variation in market characteristics, or too many potential entrants making equilibrium enumeration explode, the identified set may be wide enough to contain almost all reasonable parameters, and then "the honest interval" degenerates into "saying nothing." The response is to find stronger exogenous variation (such as an exclusion variable that affects only one firm's entry cost and not its profit) to tighten the bounds, or to combine with shape restrictions. A wide identified set is itself information: it says that data plus credible assumptions are not enough to answer the question precisely, and then what to do is admit the limitation or look for more variation, not retreat and sneak in a strong assumption in exchange for false precision.

There are several other specific pitfalls. The isolated-market assumption (each market is an independent game) may not hold in platform settings, where entry in neighboring cities can spill over and break the independence across markets, so you need to define market boundaries precisely or model the spillovers. Defining the set of potential entrants is also often overlooked, and who counts as a "potential entrant" directly affects the estimate, with too broad a definition diluting the competitive effect. The dynamic dimension of entry is left out by this chapter's static setup, and real entry and exit are intertemporal decisions involving sunk costs and expectations, so the static model gives an approximation of "in or out in the long run," and treating it as the timing of an entry-exit decision requires care, which is exactly what the next chapter's dynamic structural models handle.

Stringing these failure modes together, the credibility of an entry model comes down to two things in the end: whether, faced with multiplicity, you honestly partially identify or quietly assume a selection mechanism, and whether the modeling choices of information structure and market boundaries are defended head-on. The width of the identified set, the sensitivity to the selection rule, and the institutional grounds for the information structure are all evidence provided around these two points, and none of them can be replaced by "the model was estimated and there is a number."

## 8 Further Reading

::: {.readings}
Required reading, in suggested order:

- Bresnahan and Reiss (1991, Journal of Political Economy). The origin of the number-threshold model, read it first to understand how to infer competitive intensity from firm counts alone, without price data.
- Berry (1992, Econometrica). The heterogeneous-firm entry model, read it to understand how the "number is unique" fact sidesteps multiplicity.
- Ciliberto and Tamer (2009, Econometrica). The founding work of partial identification for entry games, the direct source of Sections 3.4 and 4.4 of this chapter, focus on the construction of the two bounds.
- Tamer (2003, Review of Economic Studies). The general treatment of incomplete models and multiplicity, the theoretical foundation of partial identification.
- Seim (2006, RAND Journal of Economics). The incomplete-information entry model, to understand the other route to recovering point identification.

Further reading:

- Berry and Reiss (2007, Handbook of Industrial Organization). A survey of entry models, placing the various methods in a unified framework.
- Pakes, Porter, Ho and Ishii (2015, Econometrica). The general framework of moment-inequality estimation, the workhorse tool of partial identification.
- Imbens and Manski (2004, Econometrica). Inference under partial identification, the distinction between covering the identified set and covering the truth.
- Andrews and Soares (2010, Econometrica). Construction of confidence sets for moment-inequality models.
- Aguirregabiria and Mira (2007, Econometrica). Estimation of dynamic games, pushing entry into the dynamic and connecting to the next chapter.
:::

::: {.apa-refs}
- Aguirregabiria, V., & Mira, P. (2007). Sequential estimation of dynamic discrete games. *Econometrica, 75*(1), 1-53. https://doi.org/10.1111/j.1468-0262.2007.00731.x
- Andrews, D. W. K., & Soares, G. (2010). Inference for parameters defined by moment inequalities using generalized moment selection. *Econometrica, 78*(1), 119-157. https://doi.org/10.3982/ECTA7502
- Berry, S. T. (1992). Estimation of a model of entry in the airline industry. *Econometrica, 60*(4), 889-917. https://doi.org/10.2307/2951571
- Berry, S. T., & Reiss, P. (2007). Empirical models of entry and market structure. In M. Armstrong & R. Porter (Eds.), *Handbook of industrial organization* (Vol. 3, pp. 1845-1886). Elsevier. https://doi.org/10.1016/S1573-448X(06)03029-9
- Bresnahan, T. F., & Reiss, P. C. (1991). Entry and competition in concentrated markets. *Journal of Political Economy, 99*(5), 977-1009. https://doi.org/10.1086/261786
- Chernozhukov, V., Hong, H., & Tamer, E. (2007). Estimation and confidence regions for parameter sets in econometric models. *Econometrica, 75*(5), 1243-1284. https://doi.org/10.1111/j.1468-0262.2007.00794.x
- Ciliberto, F., & Tamer, E. (2009). Market structure and multiple equilibria in airline markets. *Econometrica, 77*(6), 1791-1828. https://doi.org/10.3982/ECTA5368
- Imbens, G. W., & Manski, C. F. (2004). Confidence intervals for partially identified parameters. *Econometrica, 72*(6), 1845-1857. https://doi.org/10.1111/j.1468-0262.2004.00555.x
- Pakes, A., Porter, J., Ho, K., & Ishii, J. (2015). Moment inequalities and their application. *Econometrica, 83*(1), 315-334. https://doi.org/10.3982/ECTA6865
- Seim, K. (2006). An empirical model of firm entry with endogenous product-type choices. *The RAND Journal of Economics, 37*(3), 619-640. https://doi.org/10.1111/j.1756-2171.2006.tb00034.x
- Tamer, E. (2003). Incomplete simultaneous discrete response model with multiple equilibria. *The Review of Economic Studies, 70*(1), 147-165. https://doi.org/10.1111/1467-937X.00240
:::
