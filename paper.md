---
title: 'Unmixed: Linear Mixed Models and Overlap correction. An extension to the Unfold-Toolbox'
tags:
  - Matlab
  - EEG
  - overlap correction
  - item effects
  - GAM
  - non-linear effects
authors:
  - name: Benedikt V. Ehinger
    orcid: 0000-0002-6276-3332
    affiliation: 1
affiliations:
 - name: Donders Institute for Brain, Cognition and Behaviour, Radboud University, Nijmegen, Netherlands
   index: 1
date: 10 May 2019
bibliography: paper.bib
---

# Summary
A summary describing the high-level functionality and purpose of the software for a diverse, non-specialist audience
A clear statement of need that illustrates the purpose of the software
A list of key references including a link to the software archive
Mentions (if applicable) of any ongoing research projects using the software or recent scholarly publications enabled by it

Using the unmixed toolbox, researcher can analyse ERPs while taking into account temporal overlap of events, non-linear effects and random variance of subject and item effects. In addition, the use of linear mixed models allows for a complete analysis of all data, without the need to use two-stage statistics. Effectively this results in higher power and more reliable parameter estimates, at the cost of vastly increased computational time.

Most experimental design contain the repetition of a limited set of stimuli to a limited set of subjects. We call this a repeated measures design. What most researchers neglect is, that not only the subjects repeat, but also the stimuli. This additional source of dependence (the so called item-effect) is hardly corrected for in EEG research, and no tool existed to do it on overlap corrected ERPs.


The unmixed toolbox is an extension for the *unfold* toolbox. It supports all model-features described specified in the unfold toolbox [@Ehinger2019]. Based on the linear mixed model implementation of the matlab statistics toolbox, data from all subjects can be analyzed concurrently. 

# Limitations
- While the current implementation allows to specify the covariance structure of the random effects, they need to be the same for all random effects.
- The runtime for realistic datasets with high sampling rates is very high, possibly resulting in days of calculation time for a single channel
- The new bobyqa solver is faster than the matlab's offered ones, but comes at the cost of hacky-overwriting of the matlab fmunic solver function using the matlab PATH. This can possibly lead to problems if the fitting process is stopped and the path is not resetted.
- While theoretically feasible, random effects could also be used to estimate the flexibility of spline functions in non-linear effect models (which in unfold have to be a priori fixed). This is not yet implemented


# Acknowledgements

I am thankful for discussions with Phillip Alday and Olaf Dimigen.

# References
