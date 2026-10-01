# FWCCA Analysis for the OpenNeuro Emotional Music Dataset

This directory contains the MATLAB code for the task-based fMRI analysis of
Feature-Weighted Canonical Correlation Analysis (FWCCA) using the OpenNeuro
Emotional Music dataset (ds000171).

Four methods are compared:

- CCA
- G-FWCCA
- LG-FWCCA
- GL-FWCCA

The analysis is organized into five stages:

1. preprocessing of the anatomical and functional MRI data;
2. first-level GLM analysis using the training run;
3. second-view sensitivity analysis and hyperparameter selection;
4. held-out testing using the independent testing run; and
5. group-level evaluation and visualization.


## 1. Dataset

The task-based fMRI data were obtained from the OpenNeuro Emotional Music
dataset **ds000171, version 00001**:

https://openneuro.org/datasets/ds000171/versions/00001

The data used in this analysis are organized into two groups:

- `MDD`: participants with major depressive disorder;
- `ND`: non-depressed control participants.

Each subject contains one anatomical image and two functional runs.

For example, the original data for `sub-mdd02` are organized as:

```text
Data/MDD/sub-mdd02/
├── anat/
│   └── sub-mdd02_T1w.nii
└── func/
    ├── sub-mdd02_task-music_trainrun_bold.nii
    └── sub-mdd02_task-music_testrun_bold.nii
```

The two functional runs are treated as independent datasets throughout the
analysis:

```text
trainrun -> preprocessing, first-level GLM, model fitting,
            second-view construction, and hyperparameter selection

testrun  -> independent held-out testing and evaluation
```

The held-out `testrun` is not used for second-view construction,
hyperparameter selection, or model fitting.


## 2. Directory Structure

The Emotional Music fMRI analysis directory is organized as follows:

```text
OpenNeuro_ds000171/
│
├── Data/
│   ├── MDD/
│   │   ├── sub-mdd02/
│   │   │   ├── anat/
│   │   │   ├── func/
│   │   │   └── first_level_analysis/
│   │   └── ...
│   │
│   └── ND/
│       ├── sub-control01/
│       │   ├── anat/
│       │   ├── func/
│       │   └── first_level_analysis/
│       └── ...
│
├── Functions/
│   ├── Run_FWCCA_Family.m
│   ├── makeLocalKernel.m
│   ├── CV_FWCCA_task.m
│   ├── BuildActiveCCAViews.m
│   ├── CanonicalMap_to_FunctionalSystem.m
│   ├── get_driftremoved_fMRI.m
│   ├── generate_clean_mask.m
│   └── reslice_to_mask.m
│
├── 01_Preprocessing/
│   └── TwophasePreprocessing.m
│
├── 02_FirstLevelGLM/
│   └── RunGLMwithThreeConditions.m
│
├── 03_SecondViewSensitivityAnalysis/
│   ├── SecondViewSensitivity_TrainingAnalysis_for_SingleSubject.m
│   ├── SecondViewSensitivity_LG_GLFWCCA_CV_for_SingleSubject.m
│   └── SubjectResults/
│       ├── MDD/
│       │   ├── sub-mdd02/
│       │   │   ├── sub-mdd02_SecondViewInputs.mat
│       │   │   ├── sub-mdd02_CCA_GFWCCA_TrainingResults.mat
│       │   │   └── sub-mdd02_LG_GLFWCCA_CVResults.mat
│       │   └── ...
│       └── ND/
│           ├── sub-control01/
│           │   ├── sub-control01_SecondViewInputs.mat
│           │   ├── sub-control01_CCA_GFWCCA_TrainingResults.mat
│           │   └── sub-control01_LG_GLFWCCA_CVResults.mat
│           └── ...
│
├── 04_SecondViewAnalysis_Testing/
│   ├── SecondViewSensitivity_Testing_for_SingleSubject.m
│   └── Results/
│       ├── MDD/
│       │   ├── sub-mdd02/
│       │   │   ├── K1/
│       │   │   ├── K2/
│       │   │   ├── K3/
│       │   │   ├── K4/
│       │   │   └── K5/
│       │   └── ...
│       └── ND/
│           ├── sub-control01/
│           │   ├── K1/
│           │   ├── K2/
│           │   ├── K3/
│           │   ├── K4/
│           │   └── K5/
│           └── ...
│
├── 05_Evaluation/
│   ├── AverageAcrossThresholds_KSpecific_FunctionalSystems.m
│   ├── Plot_ND_vs_MDD_KSpecific_FunctionalSystems.m
│   ├── Results/
│   │   ├── MDD/
│   │   └── ND/
│   └── Figures/
│
└── README.md
```

The `Functions/` directory contains the shared FWCCA fitting,
local-kernel construction, cross-validation, second-view construction,
mask-processing, and functional-system evaluation functions used throughout
the analysis.


## 3. Analysis Workflow

The task-based fMRI analysis is organized into five sequential stages:

1. preprocessing;
2. first-level GLM analysis;
3. second-view sensitivity analysis and LG/GL-FWCCA hyperparameter selection;
4. independent held-out testing; and
5. group-level evaluation and visualization.


### 3.1---Step 01: Preprocessing---

Run:

```matlab
TwophasePreprocessing
```

The preprocessing procedure is applied separately to the `trainrun` and
`testrun`.

The preprocessing pipeline consists of two phases:

```text
Phase 1:
Realign
   |
   v
Slice Timing
   |
   v
Coregister
   |
   v
Segment

Phase 2:
Normalize
   |
   v
Smooth
```

The resulting preprocessed functional images follow the naming convention:

```text
swar<subject>_task-music_trainrun_bold.nii
swar<subject>_task-music_testrun_bold.nii
```

and are stored in the corresponding subject-specific `func/` directory.

The preprocessed `trainrun` is subsequently used for first-level GLM
analysis and model training, whereas the independently preprocessed
`testrun` is reserved for held-out testing.


### 3.2---Step 02: First-Level GLM Analysis---

Run:

```matlab
RunGLMwithThreeConditions
```

A subject-level SPM first-level GLM is fitted using the preprocessed
`trainrun` data.

The model contains three experimental conditions:

```text
tones
positive_music
negative_music
```

The subject-specific first-level results are stored under:

```text
Data/<group>/<subject>/first_level_analysis/
```

including:

```text
SPM.mat
mask.nii
spmT_*.nii
```

The positive-music and negative-music GLM statistics obtained from the
`trainrun` are subsequently used to construct the positive and negative
second views in Step 03.

The held-out `testrun` is not used for first-level GLM estimation.


### 3.3---Step 03: Second-View Sensitivity Analysis---

Step 03 constructs the second CCA views, fits CCA and G-FWCCA, and performs
cross-validation for LG-FWCCA and GL-FWCCA.

Two scripts are run sequentially:

```text
SecondViewSensitivity_TrainingAnalysis_for_SingleSubject.m
                         |
                         v
SecondViewSensitivity_LG_GLFWCCA_CV_for_SingleSubject.m
```

The provided scripts use `sub-mdd02` as an illustrative example. The same
procedure is applied independently to all MDD and control subjects.


#### Second-view definitions

Two second-view definitions are considered:

```text
positive view -> positive-music GLM statistic
negative view -> negative-music GLM statistic
```

To assess sensitivity to the second-view definition, six percentile
thresholds are considered separately for the positive and negative views:

```text
Threshold    Selected proportion

p90          top 10%
p85          top 15%
p80          top 20%
p75          top 25%
p70          top 30%
p65          top 35%
```

For example, `p90` denotes the 90th-percentile threshold of the corresponding
training-derived GLM statistic and therefore retains the top 10% of voxels.

The selected voxels are subsequently used to construct the corresponding
second CCA view.

All second-view definitions are derived exclusively from the `trainrun`.


#### Step 03A: Construct second views and fit CCA/G-FWCCA

Run:

```matlab
SecondViewSensitivity_TrainingAnalysis_for_SingleSubject
```

For each positive/negative view and percentile threshold, this script:

1. loads the preprocessed `trainrun` and first-level GLM results;
2. constructs the anatomical analysis mask;
3. identifies voxels exceeding the specified GLM percentile threshold;
4. constructs the corresponding second CCA view;
5. constructs the global temporal weight matrix;
6. fits CCA; and
7. fits G-FWCCA.

The script generates two subject-level files:

```text
sub-*_SecondViewInputs.mat
sub-*_CCA_GFWCCA_TrainingResults.mat
```

stored under:

```text
03_SecondViewSensitivityAnalysis/
└── SubjectResults/
    └── <group>/
        └── <subject>/
```

`sub-*_SecondViewInputs.mat` contains the training inputs required for the
subsequent LG-FWCCA and GL-FWCCA cross-validation.

`sub-*_CCA_GFWCCA_TrainingResults.mat` contains the fitted CCA and G-FWCCA
results for all second-view thresholds.


#### Step 03B: LG/GL-FWCCA hyperparameter selection

Run:

```matlab
SecondViewSensitivity_LG_GLFWCCA_CV_for_SingleSubject
```

LG-FWCCA and GL-FWCCA require selection of the local-kernel bandwidth `h`
and neighborhood radius `r`.

The candidate values are:

```matlab
h_list = 0.25:0.25:2;
r_list = 1;
Kcomp_max = 5;
nFold = 10;
```

Cross-validation is performed separately for:

```text
positive / negative view
          x
six percentile thresholds
          x
K = 1, 2, 3, 4, 5
```

LG-FWCCA uses the `postscale` formulation:

```text
local kernel -> global weighting
```

whereas GL-FWCCA uses the `prescale` formulation:

```text
global weighting -> local kernel
```

The selected hyperparameters and cross-validation results are saved to:

```text
sub-*_LG_GLFWCCA_CVResults.mat
```

Thus, each subject has three Step 03 output files:

```text
sub-*_SecondViewInputs.mat
sub-*_CCA_GFWCCA_TrainingResults.mat
sub-*_LG_GLFWCCA_CVResults.mat
```

These three files are subsequently used for independent testing in Step 04.


### 3.4---Step 04: Held-Out Testing---

Run:

```matlab
SecondViewSensitivity_Testing_for_SingleSubject
```

This stage evaluates CCA, G-FWCCA, LG-FWCCA, and GL-FWCCA using the
independent `testrun`.

The script loads the three subject-level training files generated in Step 03:

```text
sub-*_SecondViewInputs.mat
sub-*_CCA_GFWCCA_TrainingResults.mat
sub-*_LG_GLFWCCA_CVResults.mat
```

together with the preprocessed held-out `testrun`.

For CCA and G-FWCCA, the canonical directions fitted using the `trainrun`
are directly applied to the held-out data.

For LG-FWCCA and GL-FWCCA, the CV-selected hyperparameters corresponding
to each cumulative number of canonical components are used to refit the
models on the complete `trainrun`. The resulting canonical directions are
then applied to the held-out `testrun`.

Testing is performed separately for:

```text
positive / negative view
          x
six percentile thresholds
          x
K = 1, 2, 3, 4, 5
```

For each method, the held-out canonical maps are mapped to functional
systems using the Neuromorphometrics atlas.

The functional systems considered are:

```text
Emotion Control
Visual
Auditory
Motor
Language/Cognition
Memory
Attention
```

Results are accumulated across canonical components. For example:

```text
K1 -> component 1

K2 -> components 1 + 2

K3 -> components 1 + 2 + 3

K4 -> components 1 + 2 + 3 + 4

K5 -> components 1 + 2 + 3 + 4 + 5
```

The subject-level held-out results are stored under:

```text
04_SecondViewAnalysis_Testing/Results/<group>/<subject>/
```

with the following structure:

```text
<subject>/
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

Each workbook contains one worksheet for each second-view threshold:

```text
p90
p85
p80
p75
p70
p65
```

and each worksheet contains:

```text
System | CCA | G_FWCCA | LG_FWCCA | GL_FWCCA
```

No averaging across thresholds or subjects is performed in Step 04.

The `K = 3, 4, 5` results are subsequently used for the group-level
evaluation in Step 05.


### 3.5---Step 05: Group-Level Evaluation---

Step 05 averages the held-out functional-system results across second-view
thresholds, summarizes the results across subjects, and compares the MDD and
control groups.

Two scripts are run sequentially:

```text
AverageAcrossThresholds_KSpecific_FunctionalSystems.m
                         |
                         v
Plot_ND_vs_MDD_KSpecific_FunctionalSystems.m
```


#### Threshold averaging and group summaries

Run:

```matlab
AverageAcrossThresholds_KSpecific_FunctionalSystems
```

The analysis is performed separately for:

```text
MDD
ND
```

and for:

```text
K = 3, 4, 5
```

For each subject, method, functional system, and positive/negative view,
the functional-system counts are averaged across the six second-view
thresholds:

```text
p90
p85
p80
p75
p70
p65
```

The resulting subject-level values are then summarized across subjects
within each group using the group mean and standard deviation.

For each group, three output files are generated under:

```text
05_Evaluation/Results/<group>/
```

For example:

```text
05_Evaluation/Results/MDD/
├── MDD_SubjectLevel_ThresholdAveragedCounts.xlsx
├── MDD_Group_ThresholdAveraged_KSpecificSummary.xlsx
└── MDD_ThresholdAveraged_KSpecificSummary.mat
```

with corresponding files generated for the ND group.

The group-level Excel workbook contains separate worksheets for:

```text
K3_positive
K3_negative
K4_positive
K4_negative
K5_positive
K5_negative
```


#### MDD vs Control visualization

Run:

```matlab
Plot_ND_vs_MDD_KSpecific_FunctionalSystems
```

The script compares the threshold-averaged functional-system counts between
the MDD and control groups for:

```text
K = 3, 4, 5
```

and separately for:

```text
positive view
negative view
```

Four methods are displayed:

```text
CCA
G-FWCCA
LG-FWCCA
GL-FWCCA
```

The figures display the following six functional systems:

```text
Emotion Control
Visual
Auditory
Language/Cognition
Memory
Attention
```

The Motor system is retained in the numerical results but excluded from the
figures because its representation is negligible.

A total of six group-comparison figures are generated:

```text
K3 Positive View
K3 Negative View
K4 Positive View
K4 Negative View
K5 Positive View
K5 Negative View
```

The figures are stored under:

```text
05_Evaluation/Figures/
```

and are saved in both vector PDF and 300-dpi PNG formats.


#### Manuscript figure correspondence

The six generated figures correspond to the main-text and Supplementary
figures as follows:

```text
K3 Positive View  -> Supplementary Figure 17(a)
K3 Negative View  -> Supplementary Figure 17(b)

K4 Positive View  -> Supplementary Figure 18(a)
K4 Negative View  -> Supplementary Figure 18(b)

K5 Positive View  -> Main-text Figure 11(a)
K5 Negative View  -> Main-text Figure 11(b)
```

The `K = 5` results are presented in the main manuscript as Figure 11,
whereas the corresponding `K = 3` and `K = 4` analyses are reported as
Supplementary Figures 17 and 18, respectively.


## 4. Running the Analysis

Run the analysis in the following order:

```text
01_Preprocessing
       |
       v
TwophasePreprocessing.m
       |
       v
Preprocessed trainrun / testrun
       |
       v
02_FirstLevelGLM
       |
       v
RunGLMwithThreeConditions.m
       |
       v
Training-derived first-level GLM
       |
       v
03_SecondViewSensitivityAnalysis
       |
       +-- SecondViewSensitivity_TrainingAnalysis_for_SingleSubject.m
       |              |
       |              +-- sub-*_SecondViewInputs.mat
       |              |
       |              +-- sub-*_CCA_GFWCCA_TrainingResults.mat
       |
       +-- SecondViewSensitivity_LG_GLFWCCA_CV_for_SingleSubject.m
                      |
                      +-- sub-*_LG_GLFWCCA_CVResults.mat
       |
       v
04_SecondViewAnalysis_Testing
       |
       v
SecondViewSensitivity_Testing_for_SingleSubject.m
       |
       v
Subject-level held-out functional-system results
       |
       v
05_Evaluation
       |
       +-- AverageAcrossThresholds_KSpecific_FunctionalSystems.m
       |
       +-- Plot_ND_vs_MDD_KSpecific_FunctionalSystems.m
       |
       v
Group summaries and manuscript figures
```


# Reproducibility Note

The separation between the training and held-out testing runs is maintained
throughout the analysis.

The `trainrun` is used for:

```text
preprocessing
first-level GLM estimation
second-view construction
CCA/G-FWCCA model fitting
LG/GL-FWCCA hyperparameter selection
LG/GL-FWCCA refitting
```

The independent `testrun` is used only for:

```text
held-out projection
functional-system evaluation
group-level comparison
```

In particular, the `testrun` is not used to construct the positive or
negative second views, select LG-FWCCA or GL-FWCCA hyperparameters, or fit
the canonical directions.

The six second-view percentile thresholds are defined from the
training-derived GLM statistics and are kept fixed throughout
cross-validation and held-out testing.

For LG-FWCCA and GL-FWCCA, hyperparameter selection is performed using the
training data only. The selected hyperparameters may depend on the
cumulative number of canonical components `K`; therefore, the corresponding
training model is refitted for each `K` before projection onto the held-out
`testrun`.

The provided Step 03 and Step 04 scripts use `sub-mdd02` as an illustrative
single-subject example. The same procedure is applied independently to all
MDD and control subjects by changing the group and subject identifiers.