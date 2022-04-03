# Code and data for: <br> _Drug addiction mutations unveil a methylation ceiling in EZH2-mutant lymphoma_

__Hui Si Kwok(1), Allyson M. Freedy(1), Allison P. Siegenfeld, Julia M. Morris, Amanda L. Waterbury, Stephen M. Kissler, Brian B. Liau__

(1) These authors contributed equally

This repository was produced and is maintained by [__Stephen Kissler__](mailto:skissler@hsph.harvard.edu)

----

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

Code requires an up-to-date version of `R` with packages `tidyverse`, `ggrepel`, `purrr`, and `deSolve`. 

- `run_analysis.R`: This file runs the full analysis (run `source('run_analysis.R')` while in the `EZH2/` directory). It relies on the data stored in the `data/` directory and on the functions in `utils.R`. 
- `utils.R`: This file defines useful functions for running the analysis
- `AddictionProp_031022.csv`: this file contains a table of the raw guide proportions at each time (start=week 0, switch=week 5, end=week 8) for each drug condition.  