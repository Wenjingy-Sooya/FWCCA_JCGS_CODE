# FWCCA: Feature-Weighted Canonical Correlation Analysis for Spatial fMRI

This repository contains the MATLAB implementation and reproducibility code
for the simulation and fMRI analyses presented in our work on
**Feature-Weighted Canonical Correlation Analysis (FWCCA)**.

FWCCA extends classical canonical correlation analysis (CCA) by incorporating
feature-specific weighting and local dependence structures. The repository
includes implementations of four methods:

- CCA
- G-FWCCA
- LG-FWCCA
- GL-FWCCA

The reproducibility analyses are organized into three main components:

1. simulation studies;
2. task-based fMRI analysis using the OpenNeuro Emotional Music dataset
   (`ds000171`); and
3. resting-state fMRI analysis using the OpenNeuro dataset `ds005747`.

Each component has its own README with detailed instructions for reproducing
the corresponding analyses.


## 1. Repository Structure

The top-level directory is organized as follows:
```text
FWCCA_JCGS_CODE/
│
├── fMRISimulation/
│   ├── Functions/
│   ├── 01_MainSimulation/
│   ├── 02_z0SensitivityAnalysis/
│   ├── +z0Sensitivity/
│   └── README_Simulation.md
│
├── openNeuro_ds000171/
│   ├── Data/
│   ├── Functions/
│   ├── 01_Preprocessing/
│   ├── 02_FirstLevelGLM/
│   ├── 03_SecondViewSensitivityAnalysis/
│   ├── 04_SecondViewAnalysis_Testing/
│   ├── 05_Evaluation/
│   └── README_ds000171.md
│
├── openNeuro_ds005747/
│   ├── Data/
│   ├── Functions/
│   ├── 01_Prepare_rsfMRI_FWCCA_Inputs.m
│   ├── 02_LG_GL_FWCCA_CrossValidation.m
│   ├── 03_rsfMRI_Testing.m
│   ├── 04_ModelEvaluation/
│   └── README_ds005747.md
│
├── spm_25/
│
├── Atlas_FunctionMap_Multilabel.xlsx
│
└── README.md
```


The three main analysis directories are organized as independent
reproducibility modules. Detailed dataset information, directory structures,
analysis workflows, input/output descriptions, and execution instructions are
provided in the README within each module.


## 2. Simulation Studies

The simulation studies are located under:

```text
fMRISimulation/
```

and consist of two main analyses:

```text
fMRISimulation/
├── 01_MainSimulation/
└── 02_z0SensitivityAnalysis/
```

### Main simulation study

The main simulation study evaluates four CCA variants:

```text
CCA
G-FWCCA
LG-FWCCA
GL-FWCCA
```

across different signal-to-noise ratio (SNR) settings.

The simulation framework generates spatial activation patterns with
predefined temporal canonical signals and evaluates the recovered canonical
spatial maps using the performance measures reported in the manuscript.


### Initial-label sensitivity analysis

The `02_z0SensitivityAnalysis/` directory evaluates the sensitivity of
FWCCA to the initial label vector `z0`.

Four initialization strategies are considered:

```text
oracle
noisyOracle
dataDrivenCCA
randomMatched
```

These represent oracle, perturbed, data-driven, and random initializations,
respectively.

Detailed simulation settings, execution instructions, output descriptions,
and manuscript figure correspondence are provided in:

```text
fMRISimulation/README_Simulation.md
```


## 3. Task-Based fMRI: OpenNeuro ds000171

The task-based fMRI analysis is located under:

```text
openNeuro_ds000171/
```

This analysis uses the OpenNeuro Emotional Music dataset `ds000171` and
compares CCA, G-FWCCA, LG-FWCCA, and GL-FWCCA in participants with major
depressive disorder (MDD) and non-depressed controls (ND).

The analysis consists of five stages:

```text
Step 01: Preprocessing
   |
   v
Step 02: First-Level GLM Analysis
   |
   v
Step 03: Second-View Sensitivity Analysis
   |
   v
Step 04: Held-Out Testing
   |
   v
Step 05: Group-Level Evaluation
```

The `trainrun` is used for first-level GLM estimation, second-view
construction, model fitting, and hyperparameter selection, while the
independent `testrun` is reserved for held-out evaluation.

Two GLM-derived second-view definitions are considered:

```text
positive view
negative view
```

and sensitivity to the second-view definition is evaluated across six
percentile thresholds:

```text
p90
p85
p80
p75
p70
p65
```

Held-out canonical maps are summarized according to functional systems, and
the final group-level analysis compares the MDD and control groups for:

```text
K = 3, 4, 5
```

where `K` denotes the cumulative number of canonical components included in
the functional-system evaluation.

Detailed dataset information, directory structure, analysis workflow,
input/output descriptions, and manuscript figure correspondence are provided
in:

```text
openNeuro_ds000171/README_ds000171.md
```


## 4. Resting-State fMRI: OpenNeuro ds005747

The resting-state fMRI analysis is located under:

```text
openNeuro_ds005747/
```

This analysis evaluates CCA, G-FWCCA, LG-FWCCA, and GL-FWCCA using four
default-mode-network (DMN) seed regions:

```text
AnG
PCgG
PCu
SFG
```

The analysis consists of four stages:

```text
Step 01: Prepare FWCCA Inputs
   |
   v
Step 02: LG/GL-FWCCA Cross-Validation
   |
   v
Step 03: Held-Out Testing
   |
   v
Step 04: Model Evaluation
```

For each subject, the first half of the fMRI time series is used for model
fitting and hyperparameter selection, while the second half is reserved for
held-out testing.

The held-out evaluation includes:

```text
component-wise held-out correlations
representative-subject selection
representative spatial-overlap analysis
```

These analyses assess the correspondence between the FWCCA canonical spatial
maps and independently estimated held-out seed-correlation maps.

Detailed dataset information, directory structure, analysis workflow,
input/output descriptions, and manuscript figure correspondence are provided
in:

```text
openNeuro_ds005747/README_ds005747.md
```


## 5. Dependencies and Shared Resources

The analyses are implemented primarily in MATLAB and use SPM for fMRI
processing and atlas-based analyses.

SPM is included under:

```text
spm_25/
```

and provides the SPM functions required by the task-based and resting-state
fMRI pipelines.

The file:

```text
Atlas_FunctionMap_Multilabel.xlsx
```

is used primarily in the task-based Emotional Music analysis (`ds000171`).
It maps Neuromorphometrics atlas labels to the functional systems used in the
functional-system evaluation.

Analysis-specific functions and additional dependencies are described in the
README within each analysis directory.


## 6. Running the Analyses

The three analysis modules can be reproduced independently.

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

Each README provides the required execution order and describes the inputs
and outputs generated at each stage.

At a high level, the two real-data analyses follow the same general
train/test principle:

```text
Input preparation
       |
       v
Training-data model fitting
       |
       v
Hyperparameter selection
       |
       v
Held-out testing
       |
       v
Model evaluation and manuscript figures
```

The exact implementation differs between the task-based and resting-state
analyses and is described in the corresponding README.


## 7. Method Naming Convention

The following naming convention is used throughout the repository:

```text
CCA       Classical canonical correlation analysis

G-FWCCA   Global feature-weighted CCA

LG-FWCCA  Local-global FWCCA:
          local weighting followed by global feature weighting

GL-FWCCA  Global-local FWCCA:
          global feature weighting followed by local weighting
```

In the MATLAB implementation, the two local FWCCA variants correspond to:

```text
LG-FWCCA -> postscale
GL-FWCCA -> prescale
```

This naming convention is used consistently across the simulation,
task-based fMRI, and resting-state fMRI analyses.


# Reproducibility Note

The repository is organized so that the simulation, task-based fMRI, and
resting-state fMRI analyses can be reproduced independently.

For both real-data analyses, model construction and hyperparameter selection
are separated from held-out evaluation.

For the Emotional Music analysis (`ds000171`), the `trainrun` is used for
first-level GLM estimation, second-view construction, model fitting, and
hyperparameter selection. The independent `testrun` is reserved for held-out
evaluation.

For the resting-state analysis (`ds005747`), each subject's time series is
divided into training and testing segments. Model fitting and hyperparameter
selection are performed using the training segment, while the testing segment
is reserved for held-out evaluation.

The simulation studies use the experimental settings and random seeds
specified in the corresponding simulation scripts.

For exact parameter settings, execution order, generated outputs, and
manuscript figure correspondence, users should refer to the README within
the corresponding analysis directory.