# Overview
## Sex-specific early survival drives adult sex ratio bias in snowy plovers and impacts mating system and population growth
#### Luke J. Eberhart-Phillips, Clemens Küpper, Tom E. X. Miller, Medardo Cruz-López, Kathryn Maher, Natalie dos Remedios, Martin A. Stoffel, Joseph I. Hoffman, Oliver Krüger, Tamás Székely

##### *Proceedings of the National Academy of Sciences of the United States of America* (2017) [doi: 10.1073/pnas.1620043114](http://doi.org/10.1073/pnas.1620043114)


In this repository you can find all the necessary files needed to reproduce the analyses presented in our paper.

**`R`** *(SI Dataset)*

  - **Eberhart_Phillips_et_al._SI_Dataset.pdf** and **Eberhart_Phillips_et_al._SI_Dataset.Rmd** contains the documented code for all analyses, which can be implemented after downloading the datasets provided in the **`data`** and **`output/bootstrap/`** folders.

**`data`**

  - **chick_mark-recapture_data.txt** contains the mark-recapture field data of chicks. Each row is a single uniquely marked chick identified by their *ring*. The daily encounter history of an individual is expressed in their *ch*, where a "1" indicates that an individual was encountered, "0" indicates it was not encountered, and "." indicates that no survey took place on that day. *year* indicates the year during which an individual was monitored and *day_of_season* indicates the number of days since the start of the breeding season that an individual hatched. *sex* describes the molecular sex-type of an individual with "M" for males and "F" for females. *brood_ID* is a unique brood identifier for the family from which a chick hatched.

  - **juvenile_adult_mark-recapture_data.txt** contains the mark-recapture field data of juveniles and adults. Each row is a single uniquely marked individual identified by their *ring*. The annual encounter history of an individual is expressed in their *ch*, where a "1" indicates that an individual was encountered and "0" indicates it was not encountered. *sex* describes the molecular sex-type of an individual with "M" for males and "F" for females. *stage* describes the stage at which an individual was initially captured, where "J" indicates it was first captured as a chick, and "A" indicates it was first captured as an adult.

  - **breeding_data.txt** contains the individual reproductive histories of all marked breeding adults in the population. Each row is a nesting attempt uniquely identified by the nest *ID*. *no_chicks* expresses the number of chicks that hatched from the nest. *clutch_size* indicates the number of eggs in the nest when it was initially discovered. *year* describes the year in which the nest was active. *male* and *female* indicates the unique identity of the father and mother, respectively, with "male_NA" and "female_NA" describing cases in which the other mate was not identified.
  
**`output/bootstrap/`**

  - **AIC_table_chick_boot_out.txt** contains the bootstrap output for model selection of chick survival based on the mark-recapture analysis run in Program MARK. Each row is a *model* fitted via maximum likelihood to the bootstrapped data sample of each iteration (*iter*). *Phi* describes the model structure for fitting daily survival. *p* describes the model structure for fitting daily encounter probability. *npar* reveals the number of parameters used in a given model. *AICc* is the Akaike Information Criteria statistic corrected for small sample size. *DeltaAICc* is the difference in AICc between a given model and the best fit model of a given iteration. *weight* describes the AIC weight of a given model. *Deviance* describes the deviance of a given model.

  - **AIC_table_juvenile_adult_boot_out.txt** contains the bootstrap output for model selection of juvenile and adult survival based on the mark-recapture analysis run in Program MARK. Each row is a *model* fitted via maximum likelihood to the bootstrapped data sample of each iteration (*iter*). *Phi* describes the model structure for fitting annual survival. *p* describes the model structure for fitting annual encounter probability. *npar* reveals the number of parameters used in a given model. *AICc* is the Akaike Information Criteria statistic corrected for small sample size. *DeltaAICc* is the difference in AICc between a given model and the best fit model of a given iteration. *weight* describes the AIC weight of a given model. *Deviance* describes the deviance of a given model.

  - **ASR_boot_out.txt** contains the adult sex ratio estimates (*ASR_boot*) of each iteration of the bootstrap procedure. Each row represents an iteration (*iter*).

  - **survival_rates_boot_out.txt** contains the sex- and stage-specific survival estimates (*estimate*) of each iteration (*iter*) in the bootstrap procedure. Each row represents a given sex and stage (*sex_stage*) in a given iteration.
  
  
Last updated: June 7, 2017
