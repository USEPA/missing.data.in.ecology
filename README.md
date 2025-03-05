# missing.data.in.ecology

### A Supplementary R Package to "Missing data in ecology: Some synthesis, clarifications, and recommendations"

##### Michael Dumelle<sup>1</sup>, Rob Trangucci<sup>2</sup>, Amanda M. Nahlik<sup>1</sup>, Anthony R. Olsen<sup>1</sup>, Kathryn M. Irvine<sup>3</sup>, Karen Blocksom<sup>1</sup>, Jay M. Ver Hoef<sup>4</sup>, Claudio Fuentes<sup>2</sup> 

##### <sup>1</sup>United States Environmental Protection Agency, Office of Research and Development, Corvallis, OR, USA.
##### <sup>2</sup>Department of Statistics, Oregon State University, Corvallis, OR, USA
##### <sup>3</sup>United States Geological Survey, Bozeman, MT, USA.
##### <sup>4</sup>Marine Mammal Laboratory, Alaska Fisheries Science Center, NOAA Fisheries, Seattle, WA, USA

##### For correspondence, please email Michael Dumelle at Dumelle.Michael@epa.gov

### Abstract

In science, missing data are common and easily mishandled. When mishandled, missing data obfuscate scientific understanding. We review and synthesize several approaches for handling missing ecological data.  Generally, missing data can be grouped into one of three categories: missing completely at random (MCAR); missing at random (MAR); or missing not at random (MNAR). We review each category and pay special attention to MAR, which is quite flexible and useful but often misunderstood. We compare the benefits and drawbacks of several modern missing data approaches, including complete case analysis, deterministic imputation, single imputation, multiple imputation, and data augmentation. Multiple imputation and data augmentation perform best but are more complex than the other approaches. We clarify the important distinction between imputation and prediction and argue that using predictive metrics to evaluate imputation methods is bad statistical practice and should be avoided. We introduce a novel framework called a contingency filter, which clarifies whether missing data have a basis for measurement, and highlight its utility several contexts. We illustrate concepts throughout using wetland data from the United States Environmental Protection Agency's 2016 National Wetland Condition Assessment (NWCA). We end by providing ten explicit recommendations for ecologists to consider while handling missing data. 

### Package Overview

This supplementary R package contains all files used to generate components of the manuscript.

### Installation

To install the supplementary R package, either install the tarball provided alongside the submission or run
```r
install.packages("remotes") # if you don't already have remotes installed
remotes::install_github("michaeldumelle/spmodel.manuscript", ref = "main", dependencies = TRUE)
```

The package must be installed to view any of the files we discuss throughout this `README` on your machine. This package does not have to be installed to view any of these files if you want to look at them interactively using this GitHub repository.

### Data Availability

The NWCA 2016 data are available in its entirety at this link. The subset used in this manuscript is available at
```r
system.file("data/nwca_2016.csv", package = "missing.data.in.ecology")
```

### Numeric Results

The files required to reproduce numeric results (and figures associated with them) are available at the file path found by running
```r
system.file("scripts", package = "missing.data.in.ecology")
```

The numeric results themselves are available at the file path found by running
```r
system.file("output", package = "missing.data.in.ecology")
```

### Numeric Results

The raw figures are available at the file path found by running
```r
system.file("figures", package = "spmodel.manuscript")
```

### Disclaimer

The views expressed in this manuscript are those of the authors and do not necessarily represent the views or policies of the U.S. Environmental Protection Agency. Any mention of trade names, products, or services does not imply an endorsement by the U.S. government or the U.S. Environmental Protection Agency. The U.S. Environmental Protection Agency does not endorse any commercial products, services, or enterprises.
