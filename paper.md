---
title: 'The unmixed toolbox: Linear Mixed Models and overlap correction for EEG data'
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
Using the unmixed toolbox (an addon to the unfold toolbox, [@Ehinger2019], researcher can analyse ERPs while taking into account temporal overlap of events and random variance of subject and item effects. In addition, the use of linear mixed models allows for a complete analysis of all data, without the need to use two-stage statistics. Effectively this results in higher power and more reliable parameter estimates, but at the cost of vastly increased computational time.

In many EEG experiments, events follow each other in rapid succession. Often the resulting overlapping brain activity is unavoidable. Examples are multimodal experiments, eye-tracking coregistered studies or VR and mobile EEG studies. But also in classical studies  we can find overlap, e.g. buttonpress-ERPs have to be separated from stimulus ERPs.

Most experimental designs contain the repetition of a limited set of stimuli to a limited set of subjects. What most researchers mostly neglect is, that not only the subjects repeat, but also the stimuli. This additional source of dependence (the so called item-effect) is hardly corrected for in EEG research, and currently no tool exists to account for this, while at the same time account for overlapping ERPs.

The unmixed toolbox is an extension for the *unfold* toolbox ([@Ehinger2019]). Based on the linear mixed model implementation of the matlab statistics toolbox, data from all subjects are be analyzed concurrently. 

# Limitations
- While the current implementation allows to specify the covariance structure of the random effects, they need to be the same for all random effects.
- Currently nonlinear effects can only be included as fixed effects but not random effects
- The runtime for realistic datasets with high sampling rates is very long, possibly resulting in days of calculation time for a single channel
- While theoretically feasible, random effects could also be used to estimate the flexibility of spline functions in non-linear effect models (which in unfold have to be a priori fixed). This is not yet implemented.


# Acknowledgements

I am thankful for discussions with Phillip Alday and Olaf Dimigen.

# References
