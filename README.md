# EPTScoringAlgorithms R Package

## EPT's Score Calculation

Calculate EPT scores using one of the 9 algorithms recommended in Segal-Gordon et al. (2024). Segal-Gordon et al. recommended Algorithms 1 and 7. Algorithm 1 is set as default.

## Algorithm Details

The table below summarizes the details of the nine algorithms recommended in Segal-Gordon et al. (2024):

| Algorithm | Max_error | Min_ms | Max_ms | Winsorize | Error_treatment | Log_transformation | Score_computation | Parcel_based |
|-----------|-----------|--------|--------|-----------|-----------------|---------------------|-------------------|--------------|
| 1         | 0.5       | none   | 2SD    | no        | 600ms penalty   | no                  | G score           | no           |
| 2         | 0.4       | 300    | 2SD    | no        | 600ms penalty   | no                  | G score           | yes          |
| 3         | 0.5       | 300    | 1000   | yes       | 600ms penalty   | no                  | G score           | no           |
| 4         | 0.4       | none   | 1000   | yes       | 600ms penalty   | no                  | G score           | yes          |
| 5         | 0.4       | 300    | 1000   | yes       | 600ms penalty   | no                  | G score           | yes          |
| 6         | 0.5       | none   | 2SD    | no        | 600ms penalty   | no                  | G score           | no           |
| 7         | 0.5       | 300    | 2SD    | no        | 600ms penalty   | no                  | G score           | no           |
| 8         | 0.4       | none   | 1000   | no        | 600ms penalty   | yes                 | G score           | no           |
| 9         | 0.4       | none   | 1000   | yes       | 600ms penalty   | no                  | G score           | yes          |

## How to Use the Package

### Download

In R studio enviorment:

```R
library(devtools)   # Make sure that the devtools library is loaded
remotes::install_github("NogaS1/EPTScoringAlgorithms”)
library(EPTScoringAlgorithms)
?calc_ept_score     # For "help"
```
### Prepare the Data
Afterwards, pre process the data to fit this shape:

- `algoNum` (default: 1): The chosen algorithm.
- `blockNum` (default: max(data$block)): The number of blocks.
- `data`: The processed data, with columns:
  - `sid` (Participant's unique ID)
  - `error` (Indicator - 1 is error, 0 is correct)
  - `rt` (Reaction time)
  - `prime_cat` (The primes category)
  - `target_cat` (pos / neg)
  - `block` (optional) (Block number. If column is missing, it is calculated automatically by the "blockNum" argument.)


### Run the Scores' Calculation

```R
scores1 <- calc_ept_score(myProcessedData);
scores2 <- calc_ept_score(myProcessedData, algoNum = 7);
scores3 <- calc_ept_score(myProcessedDataWithoutBlockCol, blockNum = 5, algoNum = 3);
```

### Results

The results include EPT scores for each prime_cat, for each sid, overall and by block (for computing internal consistency).
