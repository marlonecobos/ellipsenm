ellipsenm: An R package for ecological niche charaterization using ellipsoids
================
Marlon E. Cobos, Luis Osorio-Olvera, Jorge Soberónn, Laura Jiménez, A. Townsend Peterson, Vijay Barve, and Narayani Barve

-   [Project description](#project-description)
    -   [Estatus of the project](#estatus-of-the-project)
-   [Package description](#package-description)
-   [Installing the package](#installing-the-package)

<style>
body {text-align: justify}
</style>
<br>

**This repository is for the project "Grinnellian ecological niches and ellipsoids in R" developed during the program GSoC 2019.**

<br>

Project description
-------------------

Student: *Marlon E. Cobos*

GSoC Mentors: *Luis Osorio-Olvera, Vijay Barve, and Narayani Barve*

Motivation:

Ecological niche modeling represents a set of methodological tools that allow researchers to characterize and analyze species ecological niches. Currently, these methods are used widely and their applications include disease risk mapping, climate change risk predictions, conservation biology, among others. Physiological data suggests that Grinnellian ecological niches are convex in nature and they may probably have an ellipsoidal form when multiple dimensions are considered. However, among the diverse available software to characterize ecological niches, algorithms to model these niches as ellipsoids in the environmental space are scarce. Given the need for more solutions, the **ellipsenm** package aimed for developing specialized tools to perform multiple analyses related to ecological niches of species.

### Estatus of the project

At the moment we have completed the first of three phases. We have made few modifications to the list of original products that have helped us to improve the package functionality. Next steps are to complete scheduled activities to obtain all products.

All commits made can be seen at the <a href="https://github.com/marlonecobos/ellipsenm/commits/master" target="_blank">complete list of commits</a>.

Following you can find a brief description of this R package, as well as some examples of its use (still in development).

<br>

Package description
-------------------

The **ellipsenm** R package implements multiple tools to help in using ellipsoid envelopes to model ecological niches of species. Handly options for calibrating and selecting models, producing models with replicates and projections, and assessing niche overlap are included as part of this package. Other functions implemented here are useful to perform a series of pre- and post-modeling analyses.

<br>

Installing the package
----------------------

**ellipsenm** is in a GitHub repository and can be installed and/or loaded using the following code (make sure to have Internet connection).

``` r
# Installing and loading packages
if(!require(devtools)){
    install.packages("devtools")
}
if(!require(ellipsenm)){
    devtools::install_github("marlonecobos/ellipsenm")
}
library(ellipsenm)
```

<br>
