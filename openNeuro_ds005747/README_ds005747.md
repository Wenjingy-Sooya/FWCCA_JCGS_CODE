# FWCCA Analysis for Resting-State fMRI

This directory contains the MATLAB code for the resting-state fMRI
analysis of Feature-Weighted Canonical Correlation Analysis (FWCCA).

Four methods are compared:

- CCA
- G-FWCCA
- LG-FWCCA
- GL-FWCCA

The analysis is organized into four stages:

1. preparation of subject- and ROI-specific analysis inputs;
2. cross-validation for LG-FWCCA and GL-FWCCA;
3. model fitting and held-out testing; and
4. evaluation and visualization of the held-out results.


## 1. Dataset

The resting-state fMRI data were obtained from OpenNeuro dataset
**ds005747, version 1.2.1**:

https://openneuro.org/datasets/ds005747/versions/1.2.1

The dataset contains 7T resting-state fMRI from 90 healthy adults.
The present analysis uses five subjects:

| Subject | Age | Gender |
|---------|-----|--------|
| sub-011 | 24  | M      | 
| sub-012 | 24  | M      | 
| sub-018 | 21  | M      | 
| sub-021 | 27  | F      | 
| sub-022 | 19  | F      | 

Only the **first 10-minute resting-state run** is used. The analysis
uses the preprocessed functional volumes released with the dataset.

For the volumetric analysis, these images are resliced to the
Neuromorphometrics atlas grid before ROI extraction.

The first run contains 256 time points and is split evenly into:

```text
TR 1--128   -> training data
TR 129--256 -> held-out test data
```
The held-out data are not used for hyperparameter selection.

### Software and Directory Structure

The analysis requires MATLAB and SPM25. SPM25 is expected to be located at the same directory level as openNeuro_ds005747:

## 2. Directory Structure

The resting-state fMRI analysis directory is organized as follows:

```text
openNeuro_ds005747/
│
├── Data/
│   ├── sub-011/
│   │   └── func/
│   │       ├── sub-011_task-rest_space-MNI305_preproc.nii
│   │       └── r_sub-011_task-rest_space-MNI305_preproc.nii
│   │
│   ├── sub-012/
│   │   └── func/
│   │       ├── sub-012_task-rest_space-MNI305_preproc.nii
│   │       └── r_sub-012_task-rest_space-MNI305_preproc.nii
│   │
│   ├── sub-018/
│   │   └── func/
│   │       ├── sub-018_task-rest_space-MNI305_preproc.nii
│   │       └── r_sub-018_task-rest_space-MNI305_preproc.nii
│   │
│   ├── sub-021/
│   │   └── func/
│   │       ├── sub-021_task-rest_space-MNI305_preproc.nii
│   │       └── r_sub-021_task-rest_space-MNI305_preproc.nii
│   │
│   └── sub-022/
│       └── func/
│           ├── sub-022_task-rest_space-MNI305_preproc.nii
│           └── r_sub-022_task-rest_space-MNI305_preproc.nii
│
├── Functions/
│   ├── Run_FWCCA_Family.m
│   ├── normalize_to_unit_interval.m
│   ├── makeLocalKernel.m
│   ├── CV_FWCCA_SeedCorrelation.m
│   └── BuildActiveCCAViews.m
│
├── 01_PrepareAnalysisInputs/
│   ├── Prepare_rsfMRI_FWCCA_Inputs.m
│   └── Results/
│       └── rsfMRI_FWCCA_Inputs.mat
│
├── 02_HyperparameterSelection/
│   ├── LG_GL_FWCCA_CrossValidation.m
│   └── Results/
│       └── rsfMRI_LG_GL_CV_Results.mat
│
├── 03_ModelTesting/
│   ├── rsfMRI_Testing.m
│   └── Results/
│       ├── rsfMRI_TestingResults.mat
│       └── rsfMRI_TestingGroupSummary.xlsx
│
└── 04_ModelEvaluation/
    ├── ComponentWiseHeldOutCorrelation.m
    ├── SelectRepresentativeSubjects.m
    ├── RepresentativeSpatialOverlap.m
    │
    └── Results/
        ├── ComponentWiseCorrResults.mat
        ├── RepresentativeSubjectResults.mat
        │
        ├── ComponentWiseCorrelationFigures/
        │   └── [component-wise correlation figures]
        │
        └── RepresentativeSpatialOverlap/
            ├── SpatialOverlapResults.mat
            ├── pdf/
            ├── png/
            └── fig/
```

The files prefixed by `r_` under `Data/` are generated in Step 01 by
reslicing the dataset-provided preprocessed functional images to the
Neuromorphometrics atlas grid.

The `Functions/` directory contains the shared FWCCA fitting,
normalization, local-kernel construction, cross-validation, and
second-view construction functions used throughout the analysis.

 For each subject, the required dataset-provided functional image is:

 ```
 Data/<subject>/func<subject>_task-rest_space-MNI305_preproc.nii
 ```
Subject-specific anatomical images are not required.

## 3. Analysis Workflow

The resting-state fMRI analysis is organized into four sequential stages:

1. preparation of the analysis inputs;
2. hyperparameter selection for LG-FWCCA and GL-FWCCA;
3. model fitting and held-out testing; and
4. evaluation and visualization of the held-out results.


### 3.1---Step 01: Prepare Analysis Inputs---

Run:

```matlab
Prepare_rsfMRI_FWCCA_Inputs
```

For each subject and ROI, this script:

For each subject and ROI, this script:

1. reslices the dataset-provided preprocessed functional image to the
   Neuromorphometrics atlas grid, producing

   ```text
   Data/<subject>/func/r_<subject>_task-rest_space-MNI305_preproc.nii
   ```

2. extracts the ROI voxel time series;
3. splits the first resting-state run into 128 training and 128 held-out
   time points;
4. computes voxel-wise seed correlations with the ROI mean time series;
5. selects the top 20% of training seed-correlated voxels and applies a
   `3 x 3 x 3` dilation;
6. constructs the second CCA view; and
7. constructs the global temporal weight matrix `Wg_0`.

The four bilateral DMN ROIs are:

| ROI  | Region                      | Neuromorphometrics labels |
|------|-----------------------------|----------------------------|
| AnG  | Angular Gyrus               | 106, 107                   |
| PCgG | Posterior Cingulate Gyrus   | 166, 167                   |
| PCu  | Precuneus                   | 168, 169                   |
| SFG  | Superior Frontal Gyrus      | 190, 191                   |

The prepared inputs are saved to:

```text
01_PrepareAnalysisInputs/Results/rsfMRI_FWCCA_Inputs.mat
```

under the MATLAB variable:

```text
rsfMRIInputs
```

The main stored quantities are:

| Variable | Description |
|----------|-------------|
| `X1seed_train` | First-view training data |
| `X2seed_train` | Second-view training data |
| `Wg_0` | Global temporal weight matrix |
| `Xseed_test` | Held-out ROI data |
| `rseedtrain_vals` | Training seed-correlation values |
| `rseedtest_vals` | Held-out seed-correlation values |
| `seed_labels` | Neuromorphometrics labels defining the ROI |
| `TR` | Repetition time |


### 3.2---Step 02: LG/GL-FWCCA Hyperparameter Selection---

Run:

```matlab
LG_GL_FWCCA_CrossValidation
```

LG-FWCCA and GL-FWCCA require selection of the local-kernel bandwidth
`h` and neighborhood radius `r`.

The candidate values are:

```matlab
h_list = 1:0.5:5;
r_list = 0:1:5;
Kmaxcomp = 4;
num_fold = 10;
```

The parameters are selected separately for each subject, ROI, and
number of canonical components using 10-fold voxel-wise
cross-validation on the training data.

LG-FWCCA uses the `postscale` formulation, whereas GL-FWCCA uses the
`prescale` formulation.

The cross-validation results are saved to:

```text
02_HyperparameterSelection/Results/rsfMRI_LG_GL_CV_Results.mat
```

under:

```text
rsfMRI_CV_Results
```

The main outputs are:

```text
LG_meanCumCorr
LG_BestParamTable
LG_CV_seedCorr

GL_meanCumCorr
GL_BestParamTable
GL_CV_seedCorr
```

The held-out data stored in Step 01 are not used during this
hyperparameter-selection procedure.


### 3.3---Step 03: Model Fitting and Held-Out Testing---

Run:

```matlab
rsfMRI_Testing
```

For each subject and ROI, this script:

1. loads the training and held-out inputs prepared in Step 01;
2. loads the LG-FWCCA and GL-FWCCA hyperparameters selected in Step 02;
3. fits CCA, G-FWCCA, LG-FWCCA, and GL-FWCCA using the complete training
   data;
4. projects the independent held-out data onto the fitted canonical
   directions;
5. computes the absolute Pearson correlation between each held-out
   canonical spatial map and the held-out seed-correlation map; and
6. computes cumulative absolute correlations for the first
   `K = 1, 2, 3, 4` canonical components.

For LG-FWCCA and GL-FWCCA, the cross-validation-selected `(h,r)` values
corresponding to each value of `K` are used for model refitting.

The principal outputs are:

- `AllTestMaps`: held-out seed-correlation and canonical spatial maps;
- `AllTestCorrResults`: subject-level cumulative correlations; and
- `GroupTestMeans`: group-average cumulative correlations across
  subjects.

These outputs are saved under the Step 03 `Results/` directory and are
used for the analyses in Step 04.

### 3.4---Step 04: Model Evaluation---

Step 04 evaluates and visualizes the held-out results generated in
Step 03. The evaluation consists of three stages: component-wise
correlation analysis, representative-subject selection, and spatial
overlap analysis.


#### Component-wise held-out correlations

Run:

```matlab
ComponentWiseHeldOutCorrelation
```

For each ROI and canonical component, the absolute Pearson correlation
between the held-out seed-correlation map and each canonical spatial map
is calculated separately for CCA, G-FWCCA, LG-FWCCA, and GL-FWCCA.

Group means and standard errors across the five subjects are then
computed and visualized for the first four canonical components. The
subject-level and group-level numerical results are saved in:

```text
Results/ComponentWiseCorrResults.mat
```

and the corresponding figures are saved under:

```text
Results/ComponentWiseCorrelationFigures/
```

These results correspond to **Figure 12 in the main manuscript**:

> **Figure 12:** Mean absolute correlations between the held-out
> seed-correlation map and the canonical maps for (a) AnG, (b) PCgG,
> (c) PCu, and (d) SFG. Error bars represent standard errors across
> the five subjects.


#### Representative-subject selection

Run:

```matlab
SelectRepresentativeSubjects
```

For each of the three selected ROI/component combinations, the
component-wise improvement relative to CCA is calculated separately for
LG-FWCCA and GL-FWCCA for each subject:

```text
LG gain = |corr_LG| - |corr_CCA|
GL gain = |corr_GL| - |corr_CCA|
```

The three ROI/component combinations used for the subsequent spatial
analysis are:

```text
AnG  : component 4
PCgG : component 2
PCu  : component 1
```

The representative subject is defined as the subject whose pair of
LG-FWCCA and GL-FWCCA gains is jointly closest, in Euclidean distance,
to the corresponding group-average gains.

The selected representative subjects are:

```text
AnG  / component 4 -> sub-022
PCgG / component 2 -> sub-011
PCu  / component 1 -> sub-012
```

The selected subjects and the associated component-wise results are
saved in:

```text
Results/RepresentativeSubjectResults.mat
```

These representative subjects are subsequently used in the spatial-overlap
analysis below.


#### Representative spatial-overlap analysis

Run:

```matlab
RepresentativeSpatialOverlap
```

For each representative subject/component combination, the held-out
seed-correlation map and the CCA, LG-FWCCA, and GL-FWCCA canonical maps
are first sign-aligned with the held-out seed map and normalized to
`[0,1]`.

Voxel-wise agreement between a canonical map and the held-out seed map
is measured by the perpendicular distance to the identity line:

```text
distance = |canonical map - seed map| / sqrt(2)
```

The voxel-wise improvement of each local method relative to CCA is then
defined as:

```text
LG improvement = distance_CCA - distance_LG
GL improvement = distance_CCA - distance_GL
```

Positive values indicate voxels for which the corresponding local FWCCA
map is closer to the held-out seed-correlation map than the CCA map.

For each local method, only voxels with positive improvement are
considered, and the strongest top 20% of these positive improvements
are retained. Voxels retained by both LG-FWCCA and GL-FWCCA are defined
as common-improvement voxels.

The spatial agreement between the two local methods is summarized by
the number of LG-FWCCA and GL-FWCCA selected voxels, the number of
common-improvement voxels, and the Jaccard index:

```text
Jaccard = |LG ∩ GL| / |LG ∪ GL|
```

The median spatial location of the common-improvement voxels is used to
define the representative slice location. The normalized held-out
seed-correlation map and common-improvement region are then visualized
on the SPM anatomical T1 template in sagittal, coronal, and axial views.

The numerical spatial-overlap results are saved in:

```text
Results/RepresentativeSpatialOverlap/SpatialOverlapResults.mat
```

and the corresponding figures are saved in PDF, PNG, and MATLAB FIG
formats under:

```text
Results/RepresentativeSpatialOverlap/
├── pdf/
├── png/
└── fig/
```


#### Manuscript figure correspondence

The representative spatial-overlap results correspond to one main-text
figure and two Supplementary figures:

```text
PCgG / sub-011 / component 2 -> Main-text Figure 13

AnG  / sub-022 / component 4 -> Supplementary Figure 19

PCu  / sub-012 / component 1 -> Supplementary Figure 20
```
The PCgG result is presented in the main manuscript as Figure 13. The corresponding AnG and PCu results are reported in the Supplementary Material as Figures 19 and 20, respectively.


## 4. Running the Analysis

Run the scripts in the following order:

```text
01_PrepareAnalysisInputs
        |
        v
Prepare_rsfMRI_FWCCA_Inputs.m
        |
        v
rsfMRI_FWCCA_Inputs.mat
        |
        v
02_HyperparameterSelection
        |
        v
LG_GL_FWCCA_CrossValidation.m
        |
        v
rsfMRI_LG_GL_CV_Results.mat
        |
        v
03_ModelTesting
        |
        v
rsfMRI_Testing.m
        |
        v
Held-out testing results
        |
        v
04_ModelEvaluation
        |
        +-- ComponentWiseHeldOutCorrelation.m
        |
        +-- Representative-subject selection
        |
        +-- Representative spatial-overlap analysis
```

# Reproducibility Note

The training/test separation is maintained throughout the analysis.

The first 128 time points are used to construct the training inputs,
fit the models, and select the LG-FWCCA and GL-FWCCA hyperparameters.
The remaining 128 time points are reserved for held-out evaluation.

In particular, `Xseed_test` and `rseedtest_vals` are stored during
Step 01 for use in Steps 03 and 04, but are not used during
hyperparameter selection in Step 02.