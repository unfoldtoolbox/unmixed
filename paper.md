---
title: 'The unmixed EEG toolbox: Overlap-corrected Linear Mixed Models'
tags:
  - Matlab
  - EEG
  - overlap correction
  - item effects
  - linear mixed models
  - LMM
authors:
  - name: Benedikt V. Ehinger
    orcid: 0000-0002-6276-3332
    affiliation: 1
affiliations:
 - name: Donders Institute for Brain, Cognition and Behaviour, Radboud University, Nijmegen, Netherlands
   index: 1
date: 11 September 2019
bibliography: paper.bib
---

# Summary
Using the unmixed toolbox (an addon to the unfold toolbox, @Ehinger2019), researcher can analyse event-related potentals (ERPs) while taking into account temporal overlap of events and make use of the power of linear mixed models which allows to model random variability of subjects and items, allows to include subjects with missing data and removes the need to use two-stage statistics. Effectively this results in higher power and more reliable parameter estimates, but at the cost of vastly increased computational time.

In many EEG experiments, events follow each other in rapid succession. Often the resulting overlapping brain activity is unavoidable. Examples are multimodal experiments, eye-tracking coregistered studies or VR and mobile EEG studies. But also in classical studies  we can find overlap, e.g. buttonpress-ERPs have to be separated from stimulus ERPs.

Most experimental designs contain the repetition of a limited set of stimuli to a limited set of subjects. What most researchers mostly neglect is, that not only the subjects repeat, but also the stimuli. This additional source of dependence (the so called item-effect) is hardly corrected for in EEG research, and currently no tool exists to account for this, while at the same time account for overlapping ERPs.

The unmixed toolbox is an extension for the *unfold* toolbox [@Ehinger2019]. Based on the linear mixed model implementation of the matlab statistics toolbox. The toolbox has been described in more details in a recent conference publication [@Ehinger2019ccn].

# Usage
The following section is adapted from the conference paper [@Ehinger2019ccn]. What follows describes the toolbox in a tutorial style. Mathematical details on mixed models can be found in e.g. @gelman2007. The functions of the toolbox have unit-tests to verify their correctness. The toolbox is publicly available at https://github.com/unfoldtoolbox/unmixed and is continuously developed.

## Model set up 

[! Flowchart of unmixed](doc/flowchart.pdf) 
The function ``um_designmat`` is used to specify a separate linear mixed model for each event type (e.g. stimulus onset, button press, eye movement offset). This function constructs a design matrix X and prepares the random effect design matrices based on the Wilkinson formulas, random effect specifications and the EEG data (in eeglab file format). The data can be read channel-wise from the harddrive to reduce memory footprint.
In the case of mixed models, the specification of the design matrix is made by following the extended Wilkinson formulation (lme4, matlab, statsmodels). An example formula:

*rERP ~ A + cat(B) + spl(C,5) + (1+A|subject) + (1|item)*

In this specification we would model the following parameters. Fixed effects: A will be modeled as a continuous variable, cat(B) as a categorical effect (with automatic expansion to the number of levels, e.g. using reference coding), and C as a non-linear effect using five b-splines. Random Effects: We have two grouping variables, one indicating repeated measurements within subjects (multiple trials per subject) and one indicating repeated items (one item is shown multiple times throughout all experiments). In addition, we allow the intercept and the continuous variable to vary between subjects. That is, we do not assume each subject to be affected by A in the same way. Finally, we model the correlation between the random coefficients 1 (“intercept”) and A.

The function ``um_timeexpandDesignmat`` is used to time expand the design matrix $X$ to $X_dc$, the step that allows the overlap correction (see @Ehinger2019 for details). For this, the time-limits of the events (which define the resulting rERP-length) need to be specified. The toolbox will construct one such matrix for the fixed effects and one matrix for each of the random groupings. The random effects structure I use can be understood in the following lme4 syntax:
Y~ 1+A$_0$+...+A+(1+A$_0$ |sub$_0$)+...+(1+A$_\tau$ |sub$_\tau$)
where $1$ represents an intercept, $\tau$ is the number of timelags and the brackets indicate three random effects within one of the random grouping variables (intercept, slope and correlation). Therefore, we have n random effect groupings to be estimated by $\sum_{i=1}^{N_{sub}} N_{samples,i}$ datapoints. 
The model is fitted using the ``um_mmfit`` function. I use the MatLab statistics toolbox and some custom modifications to allow sparse fixed effect design matrices. Different optimizers to solve the mixed linear model are available. Simulation based testing revealed that the external bobyqa solver [@romer] is more efficient than the matlab default quasi-newton optimizer.


# Limitations
- While the current implementation allows to specify the covariance structure of the random effects, they need to be the same for all random effects.
- Currently nonlinear effects can only be included as fixed effects but not random effects
- The runtime for realistic datasets without resampling (e.g. to 50-100 Hz) is very long, possibly resulting in days of calculation time for a single channel
- Random effects could be used to estimate the flexibility of spline functions in non-linear effect models (which in unfold have to be a priori fixed), but is not yet implemented



# Acknowledgements

I am thankful for discussions with Phillip Alday and Olaf Dimigen.

# References
