# FWCCA Analysis for the OpenNeuro Emotional Music Dataset ds000171

This directory contains the MATLAB code for the task-based fMRI analysis of
Feature-Weighted Canonical Correlation Analysis (FWCCA) using the OpenNeuro
Emotional Music dataset (`ds000171`).



The analysis compares four CCA-based methods:

- **CCA**: classical Canonical Correlation Analysis;

- **G-FWCCA**: global Feature-Weighted Canonical Correlation Analysis;

- **LG-FWCCA**: local-global Feature-Weighted Canonical Correlation Analysis;

- **GL-FWCCA**: global-local Feature-Weighted Canonical Correlation Analysis.

The two local-global variants differ in the order in which the local kernel

and global feature weights are applied: LG-FWCCA applies local weighting

before global weighting, whereas GL-FWCCA applies global weighting before

local weighting.


## 1. Dataset

The task-based fMRI data were obtained from the OpenNeuro Emotional Music
dataset **ds000171, version 00001**:

https://openneuro.org/datasets/ds000171/versions/00001

The participants included in this analysis are organized into two groups:

- `MDD`: participants with major depressive disorder;
- `ND`: non-depressed control participants.

Each participant completed five functional runs, each consisting of 105
volumes with a total duration of approximately 5 minutes and 24 seconds.

For our analysis, we selected a subset of 19 participants (9 MDD and
10 ND) who completed two full music fMRI runs, each containing the
auditory conditions used in this study: neutral tones, negative music,
and positive music. The corresponding `events.tsv` files were used to
identify these two eligible runs for each participant.

Of the two eligible runs, the run with the smaller original OpenNeuro
run number was used for training, while the run with the larger run
number was reserved for held-out evaluation. Therefore, the original
run numbers assigned to training and testing are subject-specific and
are not necessarily the same across participants.

Throughout this repository, the two selected runs are referred to as
`trainrun` and `testrun`. These labels indicate their roles in the
analysis rather than fixed OpenNeuro run numbers.

For example, `sub-mdd02`, which is used as the Path-A reproducibility
example, has the following original OpenNeuro data structure:

```text
sub-mdd02/
├── anat/
│   └── ...
│
└── func/
    ├── sub-mdd02_task-music_run-1_bold.nii.gz
    ├── sub-mdd02_task-music_run-1_events.tsv
    ├── ...
    ├── sub-mdd02_task-music_run-3_bold.nii.gz
    ├── sub-mdd02_task-music_run-3_events.tsv
    ├── ...
    ├── sub-mdd02_task-music_run-5_bold.nii.gz
    └── sub-mdd02_task-music_run-5_events.tsv
```

For `sub-mdd02`, `run-3` and `run-5` are the two eligible runs.
According to the run-assignment rule described above, they are used as:

```text
run-3  -> trainrun -> preprocessing, first-level GLM, model fitting,
                      second-view construction, and hyperparameter selection

run-5  -> testrun  -> preprocessing followed by independent held-out
                      testing and evaluation
```

This training/testing assignment is maintained throughout the subsequent
analysis.


---

## 2. Repository Structure

The main analysis structure is:

```text
openNeuro_ds000171/
│
├── 01_Preprocessing/
│
├── 02_FirstLevelGLM/
│
├── 03_SecondViewSensitivityAnalysis/
│
├── 04_SecondViewAnalysis_Testing/
│
├── 05_Evaluation/
│
├── Functions/
│
├── Atlas_FunctionMap_Multilabel.xlsx
│
├── manuscript_group_results.zip
│
└── README_ds000171.md
```

The analysis is organized into five stages:

```text
01_Preprocessing
        |
        v
02_FirstLevelGLM
        |
        v
03_SecondViewSensitivityAnalysis
        |
        v
04_SecondViewAnalysis_Testing
        |
        v
05_Evaluation
```

The repository provides the MATLAB code required for all five stages.

Common functions used by the analysis are stored under:

```text
Functions/
```

The functional-system mapping used during held-out evaluation is provided
in:

```text
Atlas_FunctionMap_Multilabel.xlsx
```


---

## 3. Complete Analysis Workflow

### 3.1 Step 01 — Preprocessing

The preprocessing code is provided under:

```text
01_Preprocessing/
```

The main preprocessing script is:

```matlab
TwophasePreprocessing.m
```

Users wishing to reproduce the complete analysis from the original
OpenNeuro data should first download `ds000171, version 00001`.

Before preprocessing, the two eligible music runs for each included
participant should be identified from the corresponding `events.tsv`
files as described in Section 1.

The original OpenNeuro run numbers are subject-specific. Among the two
eligible runs, the run with the smaller original run number is assigned
to `trainrun`, whereas the run with the larger original run number is
assigned to `testrun`.

Both selected runs undergo preprocessing. The standardized `trainrun`
and `testrun` labels are then used throughout the downstream analysis.

For example, for `sub-mdd02`:

```text
run-3 --> trainrun --> preprocessing
run-5 --> testrun --> preprocessing
```

The resulting preprocessed images are subsequently used by the
training and held-out testing analyses.

For the Path-A `sub-mdd02` example, the relevant preprocessed functional
images are:

```text
swarsub-mdd02_task-music_trainrun_bold.nii
swarsub-mdd02_task-music_testrun_bold.nii
```

Intermediate preprocessing outputs from Step 01 are not distributed
for all subjects because these files are large and can be regenerated
from the publicly available OpenNeuro data using the supplied
preprocessing code.


### 3.2 Step 02 — First-Level GLM

The first-level GLM code is provided under:

```text
02_FirstLevelGLM/
```

Run:

```matlab
RunGLMwithThreeConditions.m
```

The first-level model contains three experimental conditions:

```text
tones
positive_music
negative_music
```

The first-level GLM used for constructing the FWCCA second views is
fitted using the `trainrun` only.

Therefore:
```text
trainrun -> first-level GLM -> training-derived statistical maps -> second-view construction
```

The independent `testrun` is not used for first-level GLM estimation
for second-view construction.

For the processed `sub-mdd02` example, the Step-03-ready first-level
files include:

```text
SPM.mat
mask.nii
spmT_0002.nii
spmT_0003.nii
```

Intermediate first-level GLM outputs from Step 02 are not distributed
for all subjects. These files can be regenerated from the public
OpenNeuro data by running Steps 01 and 02.

For users primarily interested in reproducing the FWCCA methodology,
the Path-A reproducibility route described in Sections 4–5 provides the
processed Step-03-ready inputs for `sub-mdd02`.


### 3.3 Step 03 — Second-View Sensitivity Analysis

The FWCCA training and cross-validation analyses are contained in:

```text
03_SecondViewSensitivityAnalysis/
```

Step 03 consists of two parts.


#### Step 03A — Second-view construction and CCA / G-FWCCA training

Run:

```matlab
SecondViewSensitivity_TrainingAnalysis_for_SingleSubject
```

For each subject, the script:

1. loads the preprocessed `trainrun`;
2. performs the temporal processing required by the FWCCA analysis;
3. constructs the analysis mask;
4. loads the training-derived positive- and negative-music GLM
   statistical maps;
5. constructs the positive and negative second views;
6. generates the global feature weights;
7. fits CCA and G-FWCCA.

The second-view sensitivity analysis considers six percentile
thresholds:

```text
p90
p85
p80
p75
p70
p65
```

corresponding to the top:

```text
10%
15%
20%
25%
30%
35%
```

of voxels ranked by the corresponding training-derived GLM statistic
within the analysis voxels.

The selected voxels are subsequently spatially dilated using a
3 × 3 × 3 structuring element when constructing the second view.

The main Step-03A outputs are:

```text
sub-*_SecondViewInputs.mat
sub-*_CCA_GFWCCA_TrainingResults.mat
```


#### Step 03B — LG-FWCCA / GL-FWCCA cross-validation

Run:

```matlab
SecondViewSensitivity_LG_GLFWCCA_CV_for_SingleSubject
```

Cross-validation is performed using the training data only.

The analysis considers cumulative canonical-component settings:

```text
K = 1, 2, 3, 4, 5
```

The local-kernel parameter search uses:

```matlab
h_list = 0.25:0.25:2;
r_list = 1;
nFold  = 10;
```

LG-FWCCA and GL-FWCCA differ in the ordering of the local and global
feature-weighting operations.

The resulting cross-validation output is:

```text
sub-*_LG_GLFWCCA_CVResults.mat
```

Therefore, after Step 03 the subject-level analysis produces:

```text
sub-*_SecondViewInputs.mat
sub-*_CCA_GFWCCA_TrainingResults.mat
sub-*_LG_GLFWCCA_CVResults.mat
```


### 3.4 Step 04 — Held-Out Testing

The held-out testing analysis is contained in:

```text
04_SecondViewAnalysis_Testing/
```

Run:

```matlab
SecondViewSensitivity_Testing_for_SingleSubject
```

The testing stage uses the Step-03 training results together with the
independent preprocessed `testrun`.

The `testrun` is not used to:

```text
construct the second views
estimate the global feature weights
select LG-FWCCA / GL-FWCCA hyperparameters
fit the canonical directions
```

It is used only for independent held-out evaluation.

For each subject, results are generated for:

```text
K1
K2
K3
K4
K5
```

and separately for the:

```text
positive
negative
```

second views.

The resulting directory structure is:

```text
04_SecondViewAnalysis_Testing/
└── Results/
    ├── MDD/
    │   └── sub-*/
    │       ├── K1/
    │       ├── K2/
    │       ├── K3/
    │       ├── K4/
    │       └── K5/
    │
    └── ND/
        └── sub-*/
            ├── K1/
            ├── K2/
            ├── K3/
            ├── K4/
            └── K5/
```

Each K-specific directory contains:

```text
positive_MethodComparison.xlsx
negative_MethodComparison.xlsx
```

Each workbook contains six worksheets:

```text
p90
p85
p80
p75
p70
p65
```

Each worksheet contains the functional-system results for:

```text
CCA
G_FWCCA
LG_FWCCA
GL_FWCCA
```

across the seven functional systems:

```text
EmotionControl
Visual
Auditory
Motor
LanguageCognition
Memory
Attention
```


### 3.5 Step 05 — Group-Level Evaluation

The group-level analysis is contained in:

```text
05_Evaluation/
```

Step 05 uses the subject-level held-out results generated in Step 04.


#### Threshold-averaged group summaries

Run:

```matlab
AverageAcrossThresholds_KSpecific_FunctionalSystems
```

The manuscript group analysis uses:

```text
K = 3, 4, 5
```

For each subject, functional-system counts are averaged across:

```text
p90
p85
p80
p75
p70
p65
```

separately for:

```text
positive / negative second views

CCA
G-FWCCA
LG-FWCCA
GL-FWCCA
```

Group means and standard deviations are then calculated separately for
the MDD and ND groups.

The resulting files are written under:

```text
05_Evaluation/
└── Results/
    ├── MDD/
    └── ND/
```

For each group, the script generates:

```text
<group>_SubjectLevel_ThresholdAveragedCounts.xlsx

<group>_Group_ThresholdAveraged_KSpecificSummary.xlsx

<group>_ThresholdAveraged_KSpecificSummary.mat
```


#### Manuscript figures

After generating the group summaries, run:

```matlab
Plot_ND_vs_MDD_KSpecific_FunctionalSystems
```

The resulting figures are saved under:

```text
05_Evaluation/Figures/
```

The correspondence with the manuscript is:

```text
K = 3
    Positive View -> Supplementary Figure 17(a)
    Negative View -> Supplementary Figure 17(b)

K = 4
    Positive View -> Supplementary Figure 18(a)
    Negative View -> Supplementary Figure 18(b)

K = 5
    Positive View -> Figure 11(a)
    Negative View -> Figure 11(b)
```

The figures are exported as:

```text
PDF (vector)
PNG (300 dpi)
```


---

## 4. Reproducibility Options

The repository provides code for the complete analysis from raw data.
In addition, two practical reproducibility paths are provided so that
users do not need to regenerate all large intermediate files.


### 4.1 Availability of Data and Intermediate Results

The original fMRI data are publicly available from OpenNeuro as:

```text
ds000171, version 00001
```

The code for Steps 01–05 is provided in this repository.

The complete intermediate preprocessing and first-level GLM outputs
from Steps 01 and 02 are not distributed for all subjects because
these files are large and can be regenerated from the public raw data
using the provided Step-01 and Step-02 scripts.

Two additional reproducibility routes are therefore provided:

- **Path A — Subject-level FWCCA reproducibility:** a processed
  Step-03-ready example for `sub-mdd02` and corresponding reference
  outputs are provided through Zenodo: https://doi.org/10.5281/zenodo.23074275

- **Path B — Manuscript group-level reproducibility:** the all-subject
  Step-04 functional-system results required by Step 05 are provided
  directly in this GitHub repository.


### 4.2 Reproducibility Overview

```text
PATH A — Subject-level FWCCA reproducibility

Zenodo
│
├── sub-mdd02_03-ready.zip
│
└── sub-mdd02_reference_outputs.zip
          |
          v
 extract 03-ready inputs
          |
          v
       Step 03
          |
          v
       Step 04
          |
          v
 generated sub-mdd02 results
          |
          v
 compare with reference outputs


PATH B — Manuscript group-level reproducibility

GitHub repository
│
└── manuscript_group_results.zip
          |
          v
        extract
          |
          v
04_SecondViewAnalysis_Testing/Results/
          |
          v
       Step 05
          |
          v
   group summaries
          |
          v
Figure 11
Supplementary Figures 17 and 18
```


---

## 5. Path A — Subject-Level FWCCA Reproducibility

Path A reproduces the core subject-level FWCCA workflow using the
processed example:

```text
Group   : MDD
Subject : sub-mdd02
```

For this subject, the original OpenNeuro runs used in the analysis are:

```text
OpenNeuro run-3 -> trainrun
OpenNeuro run-5 -> testrun
```

The Path-A package already contains the corresponding processed data.
Therefore, users following Path A do not need to rerun the run-selection,
preprocessing, or first-level GLM stages.

Path A starts directly from Step 03.


### 5.1 Download the Path-A Packages

The Path-A reproducibility files are provided through Zenodo:

```text
DOI: 10.5281/zenodo.23074275
```

Download:

```text
sub-mdd02_03-ready.zip
sub-mdd02_reference_outputs.zip
```

The two archives serve different purposes:

```text
sub-mdd02_03-ready.zip
    -> processed inputs required to start Step 03

sub-mdd02_reference_outputs.zip
    -> reference Step-03 and Step-04 outputs for comparison
```


### 5.2 Extract `sub-mdd02_03-ready.zip`

Extract:

```text
sub-mdd02_03-ready.zip
```

directly into the:

```text
openNeuro_ds000171/
```

project directory.

Do **not** extract it into:

```text
03_SecondViewSensitivityAnalysis/
```

For example, if the repository is located at:

```text
FWCCA_JCGS_CODE/
└── openNeuro_ds000171/
```

extract the archive into:

```text
FWCCA_JCGS_CODE/openNeuro_ds000171/
```

After extraction, the nine required processed input files should appear
under:

```text
openNeuro_ds000171/
└── Data/
    └── MDD/
        └── sub-mdd02/
            ├── anat/
            │   ├── c1sub-mdd02_T1w.nii
            │   ├── c2sub-mdd02_T1w.nii
            │   └── c3sub-mdd02_T1w.nii
            │
            ├── func/
            │   ├── swarsub-mdd02_task-music_trainrun_bold.nii
            │   └── swarsub-mdd02_task-music_testrun_bold.nii
            │
            └── first_level_analysis/
                ├── SPM.mat
                ├── mask.nii
                ├── spmT_0002.nii
                └── spmT_0003.nii
```

After confirming this structure, proceed directly to Step 03.


### 5.3 Run Step 03A

Navigate to:

```text
03_SecondViewSensitivityAnalysis/
```

and run:

```matlab
SecondViewSensitivity_TrainingAnalysis_for_SingleSubject
```

For the reproducibility example, the script should be configured for:

```text
MDD / sub-mdd02
```

The generated files include:

```text
03_SecondViewSensitivityAnalysis/
└── SubjectResults/
    └── MDD/
        └── sub-mdd02/
            ├── sub-mdd02_SecondViewInputs.mat
            └── sub-mdd02_CCA_GFWCCA_TrainingResults.mat
```


### 5.4 Run Step 03B

Next, run:

```matlab
SecondViewSensitivity_LG_GLFWCCA_CV_for_SingleSubject
```

The resulting cross-validation file is:

```text
03_SecondViewSensitivityAnalysis/
└── SubjectResults/
    └── MDD/
        └── sub-mdd02/
            └── sub-mdd02_LG_GLFWCCA_CVResults.mat
```

At the end of Step 03, the generated subject-level files are:

```text
sub-mdd02_SecondViewInputs.mat
sub-mdd02_CCA_GFWCCA_TrainingResults.mat
sub-mdd02_LG_GLFWCCA_CVResults.mat
```


### 5.5 Run Step 04

Navigate to:

```text
04_SecondViewAnalysis_Testing/
```

and run:

```matlab
SecondViewSensitivity_Testing_for_SingleSubject
```

The resulting held-out results should appear under:

```text
04_SecondViewAnalysis_Testing/
└── Results/
    └── MDD/
        └── sub-mdd02/
            ├── K1/
            │   ├── positive_MethodComparison.xlsx
            │   └── negative_MethodComparison.xlsx
            ├── K2/
            │   ├── positive_MethodComparison.xlsx
            │   └── negative_MethodComparison.xlsx
            ├── K3/
            │   ├── positive_MethodComparison.xlsx
            │   └── negative_MethodComparison.xlsx
            ├── K4/
            │   ├── positive_MethodComparison.xlsx
            │   └── negative_MethodComparison.xlsx
            └── K5/
                ├── positive_MethodComparison.xlsx
                └── negative_MethodComparison.xlsx
```

There should be:

```text
5 K values × 2 views = 10 Excel workbooks
```

with six threshold worksheets per workbook:

```text
10 workbooks × 6 thresholds = 60 worksheets
```


### 5.6 Compare with the Reference Outputs

For verification, use:

```text
sub-mdd02_reference_outputs.zip
```

This archive should be extracted into a **separate comparison
directory** rather than over the outputs generated by the user's own
Step-03 and Step-04 runs.

This prevents the generated results from being overwritten.

The reference archive contains:

```text
sub-mdd02_reference_outputs/
│
├── 03_SecondViewSensitivityAnalysis/
│   └── SubjectResults/
│       └── MDD/
│           └── sub-mdd02/
│               ├── sub-mdd02_CCA_GFWCCA_TrainingResults.mat
│               └── sub-mdd02_LG_GLFWCCA_CVResults.mat
│
└── 04_SecondViewAnalysis_Testing/
    └── Results/
        └── MDD/
            └── sub-mdd02/
                ├── K1/
                │   ├── positive_MethodComparison.xlsx
                │   └── negative_MethodComparison.xlsx
                ├── K2/
                │   ├── positive_MethodComparison.xlsx
                │   └── negative_MethodComparison.xlsx
                ├── K3/
                │   ├── positive_MethodComparison.xlsx
                │   └── negative_MethodComparison.xlsx
                ├── K4/
                │   ├── positive_MethodComparison.xlsx
                │   └── negative_MethodComparison.xlsx
                └── K5/
                    ├── positive_MethodComparison.xlsx
                    └── negative_MethodComparison.xlsx
```

The reference package contains:

```text
2 MAT reference files
10 Excel reference workbooks
60 Excel worksheets
```

The large intermediate file:

```text
sub-mdd02_SecondViewInputs.mat
```

is intentionally not duplicated in the reference package because it is
generated directly by Step 03A.

The supplied reference outputs can be compared with the results
generated by the Path-A workflow.


---

## 6. Path B — Manuscript Group-Level Reproducibility

Path B reproduces the manuscript group-level summaries and figures
without requiring Steps 03 and 04 to be rerun for all subjects.

The required all-subject Step-04 results are provided directly in this
GitHub repository as:

```text
manuscript_group_results.zip
```

The archive is approximately 603 KB.


### 6.1 Extract `manuscript_group_results.zip`

Extract:

```text
manuscript_group_results.zip
```

directly into:

```text
openNeuro_ds000171/
```

The archive is structured so that the supplied results are placed under:

```text
04_SecondViewAnalysis_Testing/Results/
```

After extraction, the relevant project structure should be:

```text
openNeuro_ds000171/
│
├── 04_SecondViewAnalysis_Testing/
│   │
│   ├── SecondViewSensitivity_Testing_for_SingleSubject.m
│   │
│   └── Results/
│       ├── MDD/
│       │   ├── sub-mdd02/
│       │   └── ...
│       │
│       └── ND/
│           ├── sub-control01/
│           └── ...
│
└── 05_Evaluation/
```

The extraction adds the supplied `Results/` contents under
`04_SecondViewAnalysis_Testing/`.

It does **not** replace or remove the MATLAB scripts already present in
`04_SecondViewAnalysis_Testing/`.

The supplied Path-B results contain the subject-level outputs required
by Step 05 for:

```text
K = 3, 4, 5
```

and both:

```text
positive
negative
```

second views.

Each supplied workbook contains the six threshold worksheets:

```text
p90
p85
p80
p75
p70
p65
```


### 6.2 Run Step 05

After extracting `manuscript_group_results.zip`, navigate to:

```text
05_Evaluation/
```

and run:

```matlab
AverageAcrossThresholds_KSpecific_FunctionalSystems
```

No subject-level Step-03 or Step-04 rerunning is required for Path B.

The script reads the supplied results from:

```text
04_SecondViewAnalysis_Testing/Results/
```

and generates the group summaries under:

```text
05_Evaluation/Results/
```


### 6.3 Reproduce the Manuscript Figures

After the group summaries have been generated, run:

```matlab
Plot_ND_vs_MDD_KSpecific_FunctionalSystems
```

The resulting figures are saved under:

```text
05_Evaluation/Figures/
```

The correspondence with the manuscript is:

```text
K = 3 -> Supplementary Figure 17
K = 4 -> Supplementary Figure 18
K = 5 -> Figure 11
```

with separate positive- and negative-view panels.


---

## 7. Reproducibility Routes at a Glance

### Full reproduction from the original OpenNeuro data

```text
OpenNeuro ds000171
        |
        v
identify subject-specific eligible runs
        |
        v
Step 01 — Preprocessing
        |
        v
Step 02 — First-level GLM
        |
        v
Step 03 — FWCCA training and cross-validation
        |
        v
Step 04 — Held-out testing
        |
        v
Step 05 — Group evaluation
```

This route reproduces the complete analysis from the original public
dataset.


### Path A — Core subject-level FWCCA workflow

```text
Zenodo:
sub-mdd02_03-ready.zip
        |
        v
Step 03
        |
        v
Step 04
        |
        v
compare against:
sub-mdd02_reference_outputs.zip
```

This route reproduces the core FWCCA analysis for the processed
`sub-mdd02` example.


### Path B — Manuscript group results

```text
GitHub:
manuscript_group_results.zip
        |
        v
extract into openNeuro_ds000171/
        |
        v
Step 05
        |
        v
Figure 11
Supplementary Figures 17 and 18
```

This route reproduces the manuscript group-level summaries and figures
without rerunning the computationally expensive subject-level analyses.


---

## 8. Reproducibility Notes

The training/testing separation is maintained throughout the FWCCA
analysis.

For each participant, the two eligible music runs are identified from
the original OpenNeuro data. The eligible run with the smaller original
run number is assigned to `trainrun`, while the eligible run with the
larger original run number is assigned to `testrun`.

The `trainrun` is used for:

```text
preprocessing
first-level GLM estimation
second-view construction
global feature-weight construction
CCA / G-FWCCA model fitting
LG-FWCCA / GL-FWCCA hyperparameter selection
LG-FWCCA / GL-FWCCA refitting
```

The independent `testrun` is used for:

```text
preprocessing
held-out projection
functional-system evaluation
```

The `testrun` is not used for model construction or hyperparameter
selection.

The six second-view thresholds are defined using the training-derived
GLM statistics and are kept fixed for the corresponding downstream
analysis.

For LG-FWCCA and GL-FWCCA, hyperparameter selection is performed using
the training data only. Because the selected parameters may depend on
the cumulative number of canonical components `K`, the corresponding
training model is refitted using the selected K-specific parameters
before projection onto the independent `testrun`.

The `sub-mdd02` Path-A package is provided as a processed
single-subject example to demonstrate the core FWCCA workflow from
Step 03 through independent held-out Step 04.

The Path-B package contains the all-subject Step-04 functional-system
results required for the manuscript group analysis, allowing Step 05
and the reported group figures to be reproduced without rerunning the
subject-level analyses for every participant.


---
