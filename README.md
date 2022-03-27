# Code and data for "Drug addiction mutations unveil a methylation ceiling in EZH2-mutant lymphoma"


----

Code requires an up-to-date version of `R` with packages `tidyverse`, `ggrepel`, `purrr`, and `deSolve`. 

- `run_analysis.R`: This file runs the full analysis (run `source('run_analysis.R')` while in the `EZH2/` directory). It relies on the data stored in the `data/` directory and on the functions in `utils.R`. 
- `utils.R`: This file defines useful functions for running the analysis
- `AddictionProp_031022.csv`: this file contains a table of the raw guide proportions at each time (start=week 0, switch=week 5, end=week 8) for each drug condition.  