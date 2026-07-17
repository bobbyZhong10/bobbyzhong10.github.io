---
title: "Measurement, Digital Traces & Algorithm-Generated Data"
subtitle: "From Platform Logs and Machine Labels to Valid Inference"
seriesline: "Foundations of Information Systems Economics · Chapter 11"
series: "Foundations of Information Systems Economics: Empirical Methods for Information Systems Economics"
sidebartitle: "Foundations of Information Systems Economics"
sidebarsub: "Theory · Methods · Applications"
chapterlabel: "Chapter 11 · Measurement, Digital Traces & Algorithm-Generated Data"
---

## Introduction

Lyra wants to know whether its customer-service AI raises customer satisfaction. Human review is too expensive, so the platform lets an LLM score twenty thousand conversations and draws only a random sample of 1200 for human re-checking. The labeler looks rather good: the predicted scores are highly correlated with the human scores. The regression also comes out unusually "precise," with an estimated AI effect of $0.72$. Unfortunately, the effect obtained from the true human scores is only about $0.42$.

The problem is not the small sample; it is that the ruler changed. The AI makes the agents phrase things more completely, and the LLM then reads that writing style as extra satisfaction, so it is systematically more lenient toward the treatment group. Twenty thousand measurements do not dilute the bias; they just apply the same crooked ruler twenty thousand times. What the platform calls "naturally occurring data" is not natural either: whether an impression enters the log depends on whether the page loaded, whether the user allowed tracking, and whether the event pipeline dropped packets, and a column of labels may be the joint product of people, models, and a review regime.

This chapter follows the chain along which data are formed and asks four things: what is the concept the researcher cares about, what actually happened in reality, what did the system record, and what did the person or model encode the record into. If construct, event, log, and label are treated as synonyms, big data will only estimate the wrong object more precisely. Lyra's case will show how a random audit reveals this kind of differential measurement error, and how prediction-powered inference lets a large volume of cheap predictions contribute precision while a small number of reliable labels do the debiasing. The point is not to trust that the model is unbiased, but to know where it errs and to carry that error into the inference.

## 1 An LLM Label with Decent Accuracy and a Badly Wrong Regression

Lyra randomly gives some agents an AI copilot. If we could obtain a human satisfaction score for every conversation, regressing true satisfaction on treatment and baseline complexity gives a treatment coefficient of $0.419$ with an SE of $0.012$. This is close to the $0.400$ set in the DGP.

Human review is expensive, so the platform uses an LLM to generate a satisfaction score for every conversation. The scatter of model predictions against gold scores lines up roughly along the 45-degree line, which makes it look like a decent labeler. But if we put the LLM score directly on the left-hand side of the regression, the treatment coefficient becomes $0.724$, almost three-quarters larger.

The error comes from asymmetry across groups. The LLM's average label error is $0.038$ on control conversations and $0.251$ on treatment conversations. The model misreads the more complete phrasing that AI assistance brings as extra satisfaction, so the outcome measurement itself is affected by treatment. An ordinary test-set RMSE mixes the two groups together and need not expose this difference, which is fatal for the causal contrast.

::: {.warning}
High predictive accuracy does not guarantee valid downstream inference. Causal research cares about whether the error is correlated with treatment, covariates, and the score of the target parameter, not just about average prediction loss.
:::

## 2 A Four-Layer Measurement Model and the Estimand

### 2.1 Construct, Event, Log, and Label

Platform data pass through four steps between the concept in the researcher's mind and the final column of numbers, and treating them as synonyms is precisely where big data starts to estimate the wrong object precisely. At the top is the construct, the theoretical concept the research truly cares about, such as "customer satisfaction," which has no natural observable counterpart. One layer down are events, the things that actually happened in reality: a conversation, an escalated complaint, a refund. Below that are logs, the traces the system actually kept, the transcript and the event timestamps, because something happening does not mean it was recorded. At the very bottom are labels, the codes a person or model assigns to the log, a human rating scale or an LLM score. Every step between the four layers can go astray, and treating them as one thing is to assume by default that none of these distortions exist.

Each of these four layers corresponds to a question that must be asked separately. Construct validity asks whether the label at the bottom truly represents the theoretical concept at the top: a model that measures "polite tone" very accurately need not be measuring satisfaction. Logging completeness asks whether the step from events to logs treats all events alike, or whether some class of event is systematically dropped. Label validity asks how large the coding error is in the step from logs to labels, in which direction it leans, and whether it is correlated with treatment. Statistical correction can only reach the part of the error in the bottom two layers that the model writes out explicitly: it can correct a ruler that is off, but it cannot make a construct that is measuring the wrong object correct, and the latter is a design problem, not something a correction formula can fix.

### 2.2 Machine Labels as Outcome

The laziest use of a machine label is to put it on the left-hand side of the regression as the outcome, but whether it qualifies as a stand-in for the true outcome depends on where it errs, not on how much. Think of a machine label as one cheap measurement equal to the truth plus a reading error; let the true outcome be $Y_i$ and the machine label be

$$\widetilde Y_i = Y_i + e_i.$$

How large the error $e_i$ is is actually secondary; what matters is whether it is the same across the two groups. A ruler that is biased high everywhere by the same amount poses no threat to the causal contrast, because the treatment and control groups are lifted together, and that offset cancels cleanly when the two are differenced; what is truly fatal is a ruler that is off by more for one group and less for the other. To make that concrete: if the target regression is $Y_i = \alpha + \tau D_i + X_i'\beta + u_i$, the bias in the treatment coefficient from using $\widetilde Y$ in place of $Y$ is exactly the across-group gap in the error,

$$\mathbb{E}[e_i\mid D_i=1,X_i]-\mathbb{E}[e_i\mid D_i=0,X_i].$$

If this gap is zero, that is, if the error is classical mean-zero noise unrelated to treatment, it usually only raises the variance without changing the conditional-mean coefficient, and the regression still lands on the right place. Lyra's error is not like this: control averages an error of $0.038$ and treatment averages $0.251$, the across-group gap is large and positive, and this is nonclassical outcome error, which does not merely add noise to the estimate but directly displaces the empirical counterpart of the estimand. What truly drives the coefficient is the conditional gap after controlling for $X$: in Lyra the trickier conversations are more likely to use the AI, treatment is correlated with baseline complexity, and after controlling for it the differential error in the treatment group is about $0.30$, which is how the $0.419$ in Section 1 gets inflated to $0.724$ (a difference of exactly $0.305$), somewhat larger than the roughly $0.21$ raw group difference between $0.038$ and $0.251$. This also demonstrates in passing which number to watch: $0.038$ and $0.251$ each look like unimportant small measurement errors, and averaging them looks even more benign, but what is truly fatal is the differential error that stretches with treatment, and it is invisible in any overall metric that mixes the two groups together.

### 2.3 Machine Variables as Regressor

If a machine prediction is used not as the outcome but as a regressor on the right-hand side, the trouble takes on a different face. The most well-behaved case is classical errors-in-variables, where the error is independent of the truth and mean-zero, and it systematically pushes the slope toward zero, which is called attenuation. The intuition is simple: mixing independent noise into a regressor only raises its variance without touching its covariance with the outcome, and since the slope is roughly the covariance divided by the variance, the denominator is inflated and the slope is diluted. The trouble is that this clean conclusion of bias toward zero holds only when the error is classical; once the error is correlated with the truth, the outcome, or the treatment, the directionality of attenuation is gone, and the bias can be positive or negative and hard to size. Chapter 26 will continue this discussion of generated regressors from the angle of text-as-data and AI research; this chapter first establishes the general logic of measurement and validation. As a common IS example, using a sentiment score that an LLM extracts from reviews as an explanatory variable to predict sales makes that sentiment score a generated regressor: if its error is merely random noise, the slope of sales on sentiment will be systematically depressed, but once that score is consistently biased high for some class of product, attenuation is replaced by a bias of uncertain direction.

### 2.4 Target Population and Validation Population

Gold labels are the standard ruler used to gauge how much the cheap labels are biased, and a ruler gives a reading only where it has measured. If humans audit only the conversations the "model is least certain about," this ruler hugs only the trickiest corner of the data, and the bias it measures is that corner's bias, not the whole's; if humans audit only escalated complaints, it has never touched the error of unescalated conversations, and that region is naturally blank. So gold labels must correspond to the target population: you can use a stratified audit or active sampling to spread the ruler strategically, but the correction must carry each record's inclusion probabilities into the audit, or you will mistake one stratum's bias for the whole's. There is also a directional trap hiding here: active sampling deliberately picks the conversations the "model is least certain about" to audit, and that cluster is exactly where the model is most likely to err, so the measured error is systematically higher than the population average, and if it is not reweighted by the sampling probabilities, it will exaggerate the bias instead.

::: {.assumption data-label="Representative validation"}
Given a record's sampling strata, the probability of entering the gold-standard audit is known and positive, the human labels are reliable, and the audit is not selectively missing on account of the unobserved true label.
:::

### 2.5 Selective Labels

Some outcomes are by nature visible only after some prior decision has been taken, and this kind of missingness is entirely different from ordinary random loss. Take a small loan-approval scenario: the bank observes repayment only for applicants it approved, and the rejected batch never has a repayment label because they never received the money. So in the data used for training and evaluation, the column "will they default" exists only for the people the history approved, and whether they were approved was precisely a risk-based selection. The customer-service setting is the same: only conversations escalated to a human get a human satisfaction rating, and the decision to escalate deliberately picks the hard cases, so the large batch of conversations that the algorithm resolved smoothly and that never escalated have no human label.

The key is that it is the decision policy, not randomness, that determines which labels exist. If a model is trained and evaluated only in the region where the historical policy allowed the outcome to surface, the model's performance in other regions is in principle unidentified, because there is not even a single true value there. To reach those regions, you need exploration, a random audit, bounds, or extra structural assumptions to fill in the never-observed counterfactuals. This is exactly the opposite of what statistics usually calls missing at random: whether something is missing depends systematically on the outcome itself, which is missing not at random, and this is the root of why simple imputation cannot rescue it, because the block you want to impute is entangled with the mechanism that determines its missingness.

::: {.intuition}
The trap of selective labels is that it makes a model look very accurate in a world filtered by past decisions, and that world precisely omits the samples the decision judged "not worth approving." A bank's scorecard having a beautiful AUC among the people it approved does not mean it works among the people it rejected, because the rejected people never entered any test set. Mistaking this post-selection accuracy for the population accuracy is the most common self-deception of selective labels.
:::

The summary of this section is as follows: platform data undergo many rounds of selection from construct to label; machine variables cause different biases on the left and on the right; and the sampling mechanism and label quality of the validation sample are identification assumptions, not data-cleaning details.

## 3 Identification: When a Small Number of Gold Labels Is Enough

### 3.1 First Define the Gold Standard

The whole PPI machinery of "predictions for precision, labels for debiasing" rests on one premise: that the gold standard is truly gold. So before doing anything you have to define it clearly. A gold standard is not "finding another model to serve as reference," but an independent, reliable measurement protocol aligned with the target construct. None of these three words is decoration: independent means its errors cannot share a source with the machine labels being checked, or the two err together and the correction cannot see the problem; reliable means repeated measurement gives a stable result; aligned with the construct means it truly measures "satisfaction" itself, not a close relative like "polite tone." For subjective concepts, this protocol must also report annotator instructions, inter-rater agreement, adjudication, and the handling of uncertain labels. This step cannot be skimped: if the human labels are themselves systematically biased, the so-called correction just shuttles back and forth between two biased measurements without moving either toward the truth. The most common form of "gold that isn't really gold" is the shortcut of using a stronger LLM to score a weaker LLM as reference: the two are trained on similar corpora, their blind spots overlap heavily, and where the weak model errs the strong model often errs too, so this reference ruler is biased together with the ruler being checked, the correction cannot see the real problem, and independence is exactly what guards against this shared-source error.

### 3.2 Random Audit Identifies the Error Structure

A random audit can identify the error structure because it puts the truth and the prediction on the same batch of conversations for you to look at side by side. The randomly drawn gold sample simultaneously observes $(Y_i,\widetilde Y_i,D_i,X_i)$, so the error $e_i = \widetilde Y_i - Y_i$ goes from hidden in the dark to a column you can compute directly, and regressing it on treatment and covariates reads off the conditional mean of the label error. Lyra's 1200 audits are a simple random sample from the target population, and splitting this error column by group makes the gap of $0.038$ average in control and $0.251$ average in treatment appear at once, so the treatment group's positive label bias cannot hide. The point is: precisely because the sampling is random, these 1200 records represent all twenty thousand, and the bias measured in the audit can be extrapolated back to the population. The reason a little over a thousand records is enough is that what is being estimated here is the mean of the error, not the incidence of some rare event: the sampling precision of a mean narrows steadily with the square root of the sample size, and 1200 records is already enough to cleanly separate two group averages that differ several-fold.

### 3.3 Prediction-Powered Rectification

The basic structure of PPI is "a large volume of predictions for precision, a small number of labels for debiasing." For a parameter that can be written as an estimating equation, first estimate on the large sample using predictions, then add back the average score difference between truth and prediction on the gold sample. In the simplified expression for linear regression,

$$\widehat\beta_{PPI}=\widehat\beta_{pred,all}+\big(\widehat\beta_{gold,labeled}-\widehat\beta_{pred,labeled}\big).$$

This additive-coefficient expression is only a simplification to aid understanding; when the design moments of the full sample and the labeled sample differ, or the target parameter is nonlinear, you cannot add three off-the-shelf estimators term by term, and should instead construct the prediction score and the rectifier within a single estimating equation. If the LLM is perfect, the correction term is near zero, and the large-sample predictions provide high precision; if the LLM has a systematic bias, the gold sample estimates and corrects it. Validity comes from a representative gold sample, a known sampling mechanism, and estimating-equation correction, not from trust in the model's accuracy.

::: {.intuition}
Think of the LLM as a cheap ruler whose scale is off. Twenty thousand cheap measurements tell you the shape of the population, and a randomly drawn twelve hundred standard measurements tell you how far off the scale is. Looking only at the cheap ruler will be precisely wrong, looking only at the standard ruler will be correct but noisy, and only combining the two has a chance of getting both.
:::

This "combining the two" becomes very clear when it lands on the variance, and it is worth walking through with Lyra's numbers. The gold-only estimator uses only 1200 human labels, so its precision follows those 1200, giving an SE of 0.047. PPI takes a different route: it first uses all twenty thousand LLM predictions to compute a highly precise but biased estimate, then uses the gold sample only to estimate the correction term of "how far off it is." The key is that the correction term is an average of the difference between truth and prediction $Y_i - \widetilde Y_i$, and when the LLM broadly tracks the truth, the variance of this difference is far smaller than the variance of $Y_i$ itself. So with the same 1200 labels, PPI applies them to an object with smaller fluctuation, buying a tighter interval. In Lyra, PPI's SE is 0.037, about a fifth narrower than gold-only's 0.047, while pulling naive's 0.724 back to 0.418. This regularity is monotone: the more accurate the prediction, the smaller $Y_i - \widetilde Y_i$, the lower the variance of the correction term, and the closer PPI comes to the precision of "holding twenty thousand true labels"; conversely, when the prediction has no information, the extra term in basic PPI only adds variance and need not be better, and may even be slightly worse than gold-only; to guarantee that efficiency is never below gold-only, switch to the adaptive variant with power tuning (PPI++), which automatically pulls the weight back when the prediction is useless and falls back exactly to gold-only. This is the precise meaning of "predictions for precision, labels for debiasing" in terms of variance, and it is why it does not require believing in advance that the model is unbiased.

### 3.4 Model Training and Inferential Reuse

There is a lazy practice that is especially tempting and especially dangerous: using the same batch of gold labels both to train the model and to estimate the correction. The danger is double-dipping. The model has already memorized the answers of this batch of labels during training, and testing it again on the same data makes it look more accurate than on unfamiliar data, so the error measured in the audit is too small, PPI's debiasing amount is too small along with it, and the honesty of "not presupposing the model is unbiased" quietly leaks away. The safe practice is to let every observation's prediction come from a model that never used its own gold label: external training, sample splitting, or cross-fitting all work. At the same time, the model version, prompt, temperature, and post-processing must be fixed and recorded, or you cannot even say "which model was tested." The practice of cross-fitting is very concrete: split the data into K folds, and each fold's prediction comes from a model trained only on the remaining folds, so no one scores its own gold label, and only then is the error measured in the audit the model's true performance on unfamiliar data, rather than its score for memorizing the answers.

### 3.5 Treatment-Induced Measurement

When treatment changes how people speak, the length of the record, or the probability of being audited, the measurement process itself becomes a treatment channel, and this is the root of all of Lyra's troubles, worth slowing down to see clearly. Imagine two conversations with identical true satisfaction, one from an agent using an AI copilot with more complete and polished phrasing, the other from an agent not using it with plain phrasing. The LLM reads "complete and polished" as "more satisfied," so it gives the former a higher score. The problem is not that the model misjudges one particular conversation, but that it wears more lenient glasses on the treatment group as a whole: the same model, applied to the two groups, becomes two rulers with different scales. What Lyra's treatment-stratified audit measures is exactly this stretching: the control group averages an error of only $0.038$, but the treatment group errs by $0.251$, the same model's scale differs by nearly sevenfold across the two groups, and this difference falls exactly on the step that the causal contrast subtracts.

For this reason, the defense "I used the same model on both groups, so it's fair" does not hold. The same model, when the text style is changed by treatment, gives exactly a reading that stretches with the group, and this is how differential error arises. The only way to detect it is treatment-stratified audits, checking the machine labels against the truth separately within the two groups to see whether the error is the same size.

::: {.intuition}
Compare it to an ordinary blood-pressure-cuff calibration. If the cuff reads 5 high for everyone, that 5 is subtracted out in a group comparison, harmless. But if it reads 10 extra specifically for people who just finished running and have a high heart rate, and treatment happens to make one whole group go running, then the difference in readings between the two groups mixes in the instrument's reaction to treatment. That the measurement process becomes a treatment channel is exactly this.
:::

### 3.6 Logging and Observation Selection

Something happening does not mean it was recorded. Lyra's system is more inclined to keep conversations that were escalated, less satisfied, and from treatment, and dropping records all the way down, in the end only $24.2\%$ of conversations enter a certain class of operational log. This is like a funnel: what goes in is all conversations, and what comes out is a filtered handful. If the research reads only this log, the target population has unknowingly shrunk from "all conversations" to "recorded conversations," and the means of the two generally do not match. The remedy is either to redefine the estimand and honestly admit that you are estimating the effect on selected conversations, or to reweight each record inversely by the known logging probability, statistically filling back the batch that was dropped. But there is one situation that even weighting cannot save: if whether a conversation enters the log depends on its own unobserved satisfaction, then the population mean cannot be point-identified from observed logs alone, because what determines whether it stays or goes is exactly the quantity you want to estimate. The principle of weighting is not hard: give each recorded conversation the reciprocal of its probability of entering the log, letting it speak on behalf of the similar conversations that were not recorded, so the rarer the record, the larger the weight it receives. But for this reciprocal weight to be usable, the premise is that the logging probability is estimable and does not depend on the outcome itself, and this is exactly what the worst case above lacks.

### 3.7 The Testable and Untestable Parts

Here a line must be drawn to distinguish what the gold sample can and cannot prove. What it can prove is all inside the gold sample: since truth and prediction are paired there, you can check calibration, group-specific error, residual patterns, and model drift, all questions the label pairs can answer on their own. What it cannot prove is all in the places the gold sample cannot reach: those events that never entered the log have no shadow at all in the gold sample, so a question like "does the logging systematically drop a whole class of events" cannot be answered by it; likewise, whether the human labels represent the construct, and whether the policy quietly changed the observation process after the model was deployed, must rely on institutional evidence or extra audits, and expecting a proof to grow automatically out of the existing label pairs is unrealistic. As a contrast of what can and cannot be checked: whether the model is systematically high on the treatment group, a pair in the gold sample tells you, and this is checkable; but if there is a whole class of conversations so angry that they hang up directly and never leave a transcript, they do not even qualify to enter the gold sample, and this kind of dropped record is a dark spot the label pairs can never illuminate. Keeping this line firmly in mind saves a lot of overreaching inference of the "my calibration is very good so the measurement is fine" kind.

### 3.8 Annotator Disagreement Is Information, Not Just Noise

Annotator disagreement is often treated as noise to be averaged away, but it is mostly information. LLM labels are often trained or calibrated on a small number of human annotations, and humans need not have a unique answer for subjective concepts like ambiguity, toxicity, quality, or intent. Take a small toxicity example: for a slightly sarcastic comment, three of five annotators judge it toxic and two judge it harmless. The majority vote records "toxic" and casually throws away that three-to-two split. But this split is itself a fact: this comment really does sit in an ambiguous zone, and different but equally reasonable readers read out different meanings. Compressing it into a majority vote is to declare the construct clear, when the truth is that the construct itself carries heterogeneity. The same is true for quality judgments: for the same answer, a strict reviewer gives three stars and a lenient one gives five, and that two-star gap is not that someone misjudged, but that "what counts as good" varies from person to person.

So three situations must be kept apart: first, random annotation noise, where some annotator misses for a moment, and averaging is right; second, systematic subgroup disagreement, where annotators of some background systematically differ from another, which is real construct heterogeneity, and averaging hides it; third, genuinely set-valued labels, where there simply is no unique answer. A report therefore cannot stop at a single Cohen's kappa, and should at least give the annotation protocol, the number of annotators per sample, the adjudication rule, the group-specific confusion matrix, and whether the downstream estimate changes when the aggregation rule is switched. Especially when treatment changes the text style so that annotators disagree by different amounts on the two groups, the measurement process again becomes treatment-dependent, and this is the same thing as Section 3.5 discusses.

::: {.intuition}
Treating disagreement as noise is like forcibly reporting a bimodal distribution as its mean: the number is computed, but the most important fact, that "two camps hold opposite ends," is smoothed away by the mean. A majority vote is not measuring truth; it is forcibly manufacturing consensus for a question that had none.
:::

### 3.9 Data Provenance and a Reproducible Empirical Chain

High-standard AI measurement must withstand others rechecking it from scratch, and this requires a provenance table that pins down every link like a lab notebook: which system recorded the raw event and when, which filters it passed through, which model version, prompt, temperature, and post-processing were used, how the gold labels were drawn, and at what point in time the final variable was frozen. Any item that quietly changes may bury an invisible regime break in the data, and afterward the accounts can never be reconciled. A concrete picture is this: between two data extractions, the model behind the API was silently swapped to a new checkpoint, all scores after a certain date were lifted a notch, and because you did not record the version, you misread this purely measurement jump as a real rise in user satisfaction, and when you want to go back and verify, you cannot even answer "which model gave these scores."

What must especially be guarded against is that a reproducible model does not mean reproducible measurement. The model behind an API may update without your knowledge, the same prompt in a different context gives a different result, and one change to the privacy policy changes even the inputs fed in. When saving the entire model snapshot is impractical, the next best is to at least save the raw inputs, the full call parameters, hashes, sampled audits, and calibration outputs, so that later people can at least reconstruct the measurement error, rather than clutching only that bare column of scores in the final CSV without knowing how it came to be.

## 4 Estimation: From Naive Substitution to Valid Correction

### 4.1 Naive Plug-In

Directly using $\widetilde Y$ to replace $Y$ is the laziest move, and also the most dangerous. It quietly assumes that the machine error is mean-zero under the target estimating equation, that is, that the ruler has no systematic bias. Section 1 already played out how this assumption collapses: reporting a nice overall RMSE does not test this condition at all, because RMSE mixes the two groups' errors together and averages them, and a treatment group biased high and a control group biased low can both be covered by a moderate small RMSE, while what is fatal is exactly that averaged-away across-group gap. More precisely, for naive plug-in to hold, what is needed is not "small error" but that the error's weighted mean under the target estimating equation is zero; an error with a very small absolute value that happens to point in the same direction as treatment or the score can still bias the coefficient, whereas a very large but completely random error is harmless. Moving attention from "how large the error is" to "what the error is correlated with" is the first step in using machine labels well.

### 4.2 Gold-Only Estimation

Using only human labels is the most transparent practice, because it does not touch that biased machine ruler at all. Lyra's gold-only treatment coefficient is $0.422$ with an SE of $0.047$, essentially unbiased, landing right by the oracle's $0.419$. The cost is that it throws away all 18800 observations that have a prediction but no gold label, and the precision is discounted accordingly: the SE of $0.047$ is nearly four times the oracle's $0.012$. In other words, gold-only trades away precision for correctness, and the entire point of PPI's existence is to want to pick back up this large handful of cheap information that was thrown away, without paying with correctness. It is worth appreciating what this fourfold means: the SE shrinks with the square root of the sample size, and to catch gold-only up to the oracle's precision by piling on human labels alone, you would have to expand the 1200 records more than tenfold, an absurdly high cost, which is exactly the value of PPI taking the shortcut with cheap predictions.

### 4.3 PPI and Uncertainty Propagation

Lyra's PPI estimate is $0.418$ with a bootstrap SE of $0.037$, doing two things at once: pulling naive's $0.724$ back near the oracle, while being more precise than gold-only's $0.047$. Section 3.3 has already explained from the variance angle why it can be more precise, and here we add the place where practice most easily crashes, namely uncertainty propagation. Those twenty thousand LLM predictions are not constants that fell from the sky; they carry their own sampling error, and if you treat the trained labels as error-free fixed values and compute robust SEs only for the second-stage regression, the reported interval pretends to be precise. The correct practice is to bootstrap by resampling the large prediction sample and the small gold sample together, letting the jitter of both sources enter the final SE. Laying the three SEs side by side is most intuitive: the oracle holds all true labels, and $0.012$ is the ceiling of precision; gold-only has only 1200 true labels, $0.047$; PPI borrows twenty thousand predictions to fill in the information and compresses the SE to $0.037$, landing between gold-only and oracle, clearly tighter than relying on true labels alone, and not fooling anyone with false precision as naive does.

### 4.4 Stratified Audit and Weighting

If the platform does not sample uniformly but stratifies by treatment, language, merchant size, or model confidence, the correction follows the stratification too: estimate each stratum's own correction within each stratum, then aggregate them back by each stratum's weight in the target population. The benefit of doing this is that you can oversample rare failures, concentrating the audit firepower where problems are most likely and diagnostic value is highest. But one discipline cannot be broken: the analysis must carry each record's inclusion probability all the way through, or the oversampled rare samples will unknowingly swap your estimand for another. Stratified sampling turns "where to invest the audit" into a design problem, and the optimal solution to that problem is the subject of Section 4.7. As a concrete stratification: if Lyra worries about differential error in the treatment group, it should draw more audits within the treatment group and fewer in the control group whose error is already known to be small, then stitch the two strata's corrections back together by their respective population shares, which both saves the annotation budget and presses firepower onto the stratum that most sways the causal contrast.

### 4.5 Drift Monitoring

One validation is just a snapshot of the measurement at this moment. Models will upgrade, the user mix will change, and platform policy will adjust, so the bias measured accurately yesterday need not still be accurate today. The sound practice is to keep rolling random audits by version and time, define in advance which quantities the drift metrics should watch, and once the measurement rule changes, treat it as a regime change in the data-generating process rather than handing it to an ordinary time fixed effect in the hope that it will be absorbed automatically. A change in the measurement rule directly breaks the implicit premise that "the same number represents the same thing at different times," and a fixed effect cannot make it back.

### 4.6 Binary Labels, Misclassification, and the Causal Contrast

Many digital traces are not continuous scores but binary labels like churn, fraud, toxicity, or adoption. Let the true state be $Y\in\{0,1\}$ and the observed label be $\widetilde Y$. The false-positive rate and false-negative rate are, respectively,

$$\mathrm{FPR}_d=P(\widetilde Y=1\mid Y=0,D=d),\qquad
\mathrm{FNR}_d=P(\widetilde Y=0\mid Y=1,D=d).$$

If the misclassification rates do not change with treatment, the observed group difference is usually compressed; but as long as $\mathrm{FPR}_1\neq\mathrm{FPR}_0$ or $\mathrm{FNR}_1\neq\mathrm{FNR}_0$, the direction of the bias can be arbitrary.

Let a small binary example nail this down. Suppose some feature actually slightly harms satisfaction: the true satisfied share is $55\%$ in control and $50\%$ in treatment, so the true effect is $-5$ percentage points. The observed satisfied share is given by $P(\widetilde Y=1\mid D=d)=P(Y=1\mid D=d)(1-\mathrm{FNR}_d)+P(Y=0\mid D=d)\,\mathrm{FPR}_d$. If the classifier has both error types at $10\%$ in control, the observed share is $0.55\times0.9+0.45\times0.1=0.54$. But as in Lyra, the phrasing in the treatment group is more complete, and the classifier more easily judges the actually dissatisfied as satisfied, so $\mathrm{FPR}_1$ rises to $30\%$ (with $\mathrm{FNR}_1$ still $10\%$), and the observed share becomes $0.50\times0.9+0.50\times0.3=0.60$. So the observed group difference is $+6$ percentage points while the truth is $-5$: a feature that lowers satisfaction is reported as significantly raising it, and the sign is entirely flipped. What makes it deadly is that pooled accuracy or F1 sees none of this, because it mixes the two groups' errors together and averages them away. This is exactly why the confusion matrix must be reported by treatment group, rather than giving only an overall accuracy or F1 score.

When the validation sample can identify each group's misclassification matrix, you can do a matrix inversion on group prevalence; but if sensitivity and specificity are close to the indistinguishable boundary, the inversion is extremely unstable. Here the honest practice is to report the effect bounds induced by the confidence region of the misclassification rates, rather than giving a seemingly precise corrected point estimate.

### 4.7 Audit Allocation Is a Statistical Design Problem

With a human-annotation budget held fixed, a simple random audit need not be the most cost-effective way to spend it. The reasoning is plain: if some rare class, some particular language, or the treatment group's error has the largest effect on the target parameter, spreading the audit evenly across every stratum is a waste, and the money should be piled onto the stratum that matters most. As for how to pile it, the classic Neyman allocation gives the answer. Let stratum $h$ have population weight $W_h$, target-score standard deviation $S_h$, and per-audit cost $c_h$; the cost-adjusted optimal allocation satisfies

$$n_h\propto \frac{W_hS_h}{\sqrt{c_h}}.$$

This formula is actually a very understandable little account. The $W_h$ in the numerator says: the more heavily this stratum weighs in the population, the greater its bias's effect on the final number, and the more it should be audited. The $S_h$ in the numerator says: the more uneven this stratum's target score, the more a few audits fail to pin down its uncertainty, and the more it should be audited; conversely, a stratum highly consistent inside is understood after a few audits, and more investment is a waste. The $\sqrt{c_h}$ in the denominator says: the more expensive each audit in this stratum, the fewer samples the same money buys, and the less it should be audited, and since cost enters through a square root, the penalty on expensive strata is decreasing. Putting the three together, the conclusion is to invest the audit in the stratum that is "large in share, messy inside, and cheap." As a contrast: with a large and noisy treatment group and a small and tidy never-treated group before you, the budget clearly should be pressed onto the former. There is one more meaning to make explicit: if your goal is the treatment effect rather than overall prediction accuracy, the $S_h$ in the formula should not be computed by prediction precision but defined by the rectifier's contribution to the treatment score, and whoever contributes more to the variance of the causal contrast should get more audits. Finally, a reminder: active learning that samples by model uncertainty can feed the model itself better, but it does not automatically improve the precision of the causal estimand; and once it fails to record inclusion probabilities, it will also casually destroy the representativeness of the validation sample, bending the standard ruler of Section 2.4.

::: {.intuition}
Think of it as assigning quality inspectors. A high-output line needs more people because its defects directly drag down the overall pass rate; a line whose quality swings needs more people because a few spot checks cannot pin it down; and a line where each inspection is especially laborious needs fewer people because the same labor can check more elsewhere. Neyman allocation is nothing more than twisting these three intuitions into a single proportionality.
:::

### 4.8 Measurement Invariance, Drift, and Versioning

Any cross-period comparison quietly assumes one thing: that the same recorded value represents the same construct at different times, which is measurement invariance. Once this assumption breaks, a time trend in the metric may be purely the measurement drifting rather than the world changing. And platform metrics happen to be the least stable, with the UI, logging schema, model version, and user mix all moving. LLM labels in particular must be watched, and the model identifier, prompt, temperature, decoding rule, input cleaning, and call time should all be stored; keeping only the final labels makes the measurement process unreproducible.

Drift monitoring cannot look only at the prediction distribution. At a minimum it needs to check covariate drift, label drift, and error drift together. The first two can be seen from the large-sample logs; the third must rely on continuous gold audits. The distinction among the three is worth remembering: covariate drift is that the inputs fed in changed, label drift is that the distribution of scores spat out changed, and error drift is that the relationship between scores and truth changed, and only the last directly threatens the estimand, and it happens to be the only one that cannot be detected without seeing the truth, which is also why watching only the prediction distribution misses the most fatal kind. If the treatment rollout and the model upgrade happen at the same time, the version change is itself differential measurement and cannot be absorbed automatically by time fixed effects. Feasible designs include freezing the main-analysis model, running the old and new versions on an overlapping sample, building a bridge sample, and folding the calibration difference before and after the version switch into the final uncertainty.

## 5 Anchor Papers

### 5.1 Lakkaraju et al. (2017)

::: {.case}
Paper: "The Selective Labels Problem: Evaluating Algorithmic Predictions in the Presence of Unobservables," KDD.

Method: analyzes the problem of model evaluation when historical human decisions determine whether the outcome is observable, and exploits differences across decision-makers to expand the evaluable region.

Data: uses high-stakes settings such as judicial decisions to illustrate the selective-label problem, where outcomes are observed only for those permitted to act.

Result: in regions where the historical policy never took a certain action, prediction quality cannot be evaluated with ordinary held-out labels.

Limitation: exploiting decision-maker heterogeneity requires a corresponding assignment and exclusion argument, and not all judge variation can be automatically treated as exogenous.
:::

### 5.2 Wang, McCormick, and Leek (2020)

::: {.case}
Paper: "Methods for Correcting Inference Based on Outcomes Predicted by Machine Learning," Proceedings of the National Academy of Sciences.

Method: uses independent validation data to estimate the relationship between the predicted outcome and the true outcome, and corrects the bias and uncertainty of the downstream regression.

Data: the paper uses simulations and applications to show the error propagation when predicted outcomes enter statistical inference directly.

Result: treating predictions as the truth produces biased coefficients and overly narrow intervals; explicitly using validation information can restore valid inference.

Limitation: the correction depends on the comparability of the validation population and the target population, and revalidation is needed under distribution shift.
:::

### 5.3 Angelopoulos et al. (2023)

::: {.case}
Paper: "Prediction-Powered Inference."

Method: combines arbitrary black-box predictions with a small number of gold labels, and through rectification constructs valid intervals for means, quantiles, and regression parameters.

Data: covers many kinds of tasks including census, remote sensing, genomics, astronomy, and ecology.

Result: the more accurate the prediction, the narrower the interval usually is; even if the predictor is biased, as long as the gold sample is well designed, rectification still guarantees inferential validity.

Limitation: the sampling and quality of the gold labels are still the foundation. If the validation set enters selectively or comes from a different population, the PPI formula will not automatically fix transportability.
:::

## 6 A Full Walkthrough on the Lyra Data

### 6.1 Data Generation

Lyra has 20000 conversations, 1200 of which obtain gold labels via a random audit. The LLM's average error is $0.038$ for control and $0.251$ for treatment.

```r
true_satisfaction <- 3 + 0.40*treatment + 0.55*x + rnorm(N, sd=0.80)
llm_satisfaction <- true_satisfaction + 0.30*treatment - 0.20*x +
  rnorm(N, sd=0.55)
gold[sample.int(N, 1200)] <- 1L
```

![Human satisfaction against LLM labels. Both groups line up overall along the 45-degree line, but the treatment observations sit systematically above the line; this differential error enters the across-group causal contrast directly.](assets/fig/fig_11_label_error.svg)

### 6.2 Four Estimates

```r
b_naive <- coef_ols(dat$llm_satisfaction, X)
b_gold <- coef_ols(dat$true_satisfaction[g], X[g,])
b_pred_gold <- coef_ols(dat$llm_satisfaction[g], X[g,])
b_ppi <- b_naive + (b_gold - b_pred_gold)
```

| Estimator | Treatment coefficient | SE |
|---|---:|---:|
| Oracle true outcome | 0.419 | 0.012 |
| LLM label naive | 0.724 | 0.014 |
| Gold labels only | 0.422 | 0.047 |
| PPI rectified | 0.418 | 0.037 |

![Four estimates of the same treatment effect. Using LLM labels directly gives a precise but severely inflated coefficient; gold-only is correct but wider; PPI uses large-sample predictions plus a gold correction to return near the oracle.](assets/fig/fig_11_logging_funnel.svg)

### 6.3 The Overall Reconciliation

The real lesson is not that "machine labels cannot be used," but that predictions cannot be treated as observed truth. The small SE of the naive estimate $0.724$ only reflects that there are many twenty thousand machine labels, not that the labels have differential bias against treatment. The random audit discovers the error, gold-only obtains $0.422$, PPI obtains $0.418$, both close to the oracle's $0.419$.

## 7 Failure Modes and Robustness

First, construct mismatch. A model accurately predicting "polite tone" does not mean it measures satisfaction. You should write the construct definition before modeling and align the audit protocol with it.

Second, nonrepresentative gold sample. Auditing only the most extreme or most uncertain cases changes the target of the correction. If a stratified design is used, the sampling weights must be saved.

Third, label leakage. If the prompt or the model input contains treatment-specific metadata, the label mechanically encodes treatment and produces a differential error like Lyra's. You should do blinded labeling and input ablation.

Fourth, model drift. Model upgrades, prompt changes, and shifts in the platform's user mix all invalidate the old validation. The version number and the audit date should enter the reproduction record like the data date.

Fifth, selective logging. Events that never entered the log cannot be corrected from recorded data alone. Random instrumentation audits and system-level reconciliation are more fundamental than statistical imputation.

Sixth, error propagation being ignored. The first-stage uncertainty of generated outcomes, regressors, or treatments must enter the final SE; robust SEs for the second-stage regression alone are not enough.

Finally, prediction and causal identification are two different things. PPI can correct outcome measurement, but it will not make an endogenous treatment exogenous. The assignment mechanism is still determined by the research design of Chapter 10.

## 8 Further Reading

::: {.readings}
Required reading:

- Lakkaraju et al. (2017). Understand why selective labels is not an ordinary test-set missingness.
- Wang, McCormick, and Leek (2020). Understand the bias and correction after predicted outcomes enter inference.
- Angelopoulos et al. (2023). Understand the unified framework of predictions plus gold-label rectification.

Further reading:

- Fuller (1987). Classic measurement-error theory.
- Wei (2021). Selective-label learning in online decision environments.
:::

::: {.apa-refs}
- Angelopoulos, A. N., Bates, S., Fannjiang, C., Jordan, M. I., & Zrnic, T. (2023). Prediction-powered inference. *Science, 382*(6671), 669-674. https://doi.org/10.1126/science.adi6000
- Fuller, W. A. (1987). *Measurement error models*. Wiley.
- Lakkaraju, H., Kleinberg, J., Leskovec, J., Ludwig, J., & Mullainathan, S. (2017). The selective labels problem: Evaluating algorithmic predictions in the presence of unobservables. In *Proceedings of the 23rd ACM SIGKDD International Conference on Knowledge Discovery and Data Mining* (pp. 275-284). https://doi.org/10.1145/3097983.3098066
- Wang, S., McCormick, T. H., & Leek, J. T. (2020). Methods for correcting inference based on outcomes predicted by machine learning. *Proceedings of the National Academy of Sciences, 117*(48), 30266-30275. https://doi.org/10.1073/pnas.2001238117
- Wei, D. (2021). Decision-making under selective labels: Optimal finite-domain policies and beyond. *Proceedings of Machine Learning Research, 139*, 11035-11046.
:::
