# FWCCA: Feature-Weighted Canonical Correlation Analysis for Spatial fMRI

This repository contains the MATLAB implementation and reproducibility code
for the simulation and fMRI analyses presented in our work on
**Feature-Weighted Canonical Correlation Analysis (FWCCA)**.

FWCCA extends classical canonical correlation analysis (CCA) by incorporating
feature-specific weighting and local dependence structures. The repository includes implementations of classical CCA and three FWCCA
variants:

- **CCA**: classical canonical correlation analysis;
- **G-FWCCA**: global feature-weighted CCA, which incorporates global
  feature weights;
- **LG-FWCCA**: local-global FWCCA, which first applies local dependence
  weighting and then global feature weighting; and
- **GL-FWCCA**: global-local FWCCA, which first applies global feature
  weighting and then local dependence weighting.

The reproducibility analyses are organized into three main components:

1. simulation studies;
2. task-based fMRI analysis using the OpenNeuro Emotional Music dataset
   (`ds000171`); and
3. resting-state fMRI analysis using the OpenNeuro dataset `ds005747`.

Each component has its own README with detailed instructions for reproducing
the corresponding analyses.


## 1. Repository Structure

The repository contains three independent reproducibility modules together
with the shared SPM dependency:

```text
FWCCA_JCGS_CODE/
│
├── fMRISimulation/
│   └── README_Simulation.md
│
├── openNeuro_ds000171/
│   └── README_ds000171.md
│
├── openNeuro_ds005747/
│   └── README_ds005747.md
│
├── spm_25/
│
└── README.md
```

The three analysis directories can be reproduced independently. Detailed
dataset information, analysis workflows, required inputs, generated outputs,
and execution instructions are provided in the README within each module.


## 2. Reproducibility Data on Zenodo

Large reproducibility files that are not distributed through GitHub because
of their size are available from the accompanying Zenodo archive:

https://doi.org/10.5281/zenodo.23074275

The Zenodo archive provides prepared inputs and reference outputs that allow
users to reproduce selected downstream analyses without rerunning the most
computationally intensive data-preparation or Monte Carlo stages.

The archive supports the three reproducibility modules as follows:

- `fMRISimulation/`: large Monte Carlo result files for the main simulation
  and the initial-label sensitivity analysis;
- `openNeuro_ds000171/`: prepared inputs and reference outputs for the
  Emotional Music analysis; and
- `openNeuro_ds005747/`: the prepared resting-state fMRI input file generated
  by Step 01.

The README within each analysis directory specifies the exact Zenodo file to
download, where it should be placed, and which analysis step can be started
from that file.


## 3. SPM25

The fMRI analyses require SPM25. SPM25 is not distributed as part of this
repository and should be downloaded separately from the official SPM website.

After downloading SPM25, place the `spm_25` directory at the same directory
level as the three analysis modules:

```text
FWCCA_JCGS_CODE/
│
├── fMRISimulation/
├── openNeuro_ds000171/
├── openNeuro_ds005747/
├── spm_25/
└── README.md
```

The analysis scripts use this relative directory structure to locate SPM25.


## 4. Running the Analyses

The three reproducibility modules can be run independently.

For the simulation studies, follow:

```text
fMRISimulation/README_Simulation.md
```

For the task-based Emotional Music fMRI analysis, follow:

```text
openNeuro_ds000171/README_ds000171.md
```

For the resting-state fMRI analysis, follow:

```text
openNeuro_ds005747/README_ds005747.md
```

Each module README provides the required data source, execution order,
analysis-specific dependencies, input/output locations, and available
reproducibility shortcuts.

