# Matrix_ATLAS
This is a repo to test the approach used by ATLAS to extract the correlation matrix among various 1D differential XS. The approach basically relies on summing two covariance matrices:

$$
C = C_{\rm stat} + C_{\rm syst}
$$

The datacards used for the HIG-23-014 paper are needed, they can be downloaded from:
```
git clone https://gitlab.cern.ch/cms-analysis/hig/hig-23-014/datacards.git
```

The folder `fidXS` contains the XS in each bin, it is needed to extract the fiducial acceptance. It is copied from https://github.com/JaLuka98/Hgg-PartialRun3-3A-ETH-Analysis/tree/main/spectrum_plotter/fidXS [TBC they may not be the final numbers used in the results]

## Statistical covariance matrix $C_{\rm stat}$
This is computed starting from counting the events in sidebands in bootstrap replica $C_{\rm stat}^{\rm sb}$. The statistical covariance matrix is obtained through the formula:

$$
C_{\rm stat} = R^{-1}~C_{\rm stat}^{\rm sb}~(R^{-1})^T
$$

where $R$ is the reponse matrix.

The steps to extract this matrix are encoded in the following scripts:
* `bootstrap.py`: this script produces bootstrap replicas of the dataset and count events in the sidebands
* `fiducial_acceptance.py`: this script extracts the fiducial acceptance from the computed fidXS for HIG-23-014
* `extract_response_matrix.py`: this script extracts the $(\epsilon\times A)$ matrix from the workspaces
  * I am extracting $(\epsilon\times A)$, but we only need $\epsilon$. The previous script computed the fiducial acceptances $A$ and here we perform the division $(\epsilon\times A)/A$
* `correlation_matrix.ipynb`: this notebook computes $C_{\rm stat}^{sb}$ and $C_{\rm stat}$
* (optional) ```extract_hessian.py```: this script extract the hessian matrix from the combined fit (= the combined fit is obtained by combining the various variables as if they were independent and then computing the covariance matrix from the Hessian)


List of commands:

```
python3 bootstrap.py --n-replicas 10000 --output bootstrap/bootstrap_sideband_counts.txt
```

```
python3 fiducial_acceptance.py 
```

Make sure to run `cmsenv` before trying the following commands. 

```
python3 extract_response_matrix.py \
  --datacard datacards/Datacard_PTH_2022_2023_2024.root \
  --out-prefix response_matrix_combined/response_matrix_PTH \
  --order "PTH_0p0_5p0,PTH_5p0_10p0,PTH_10p0_15p0,PTH_15p0_20p0,PTH_20p0_25p0,PTH_25p0_30p0,PTH_30p0_35p0,PTH_35p0_45p0,PTH_45p0_60p0,PTH_60p0_80p0,PTH_80p0_100p0,PTH_100p0_120p0,PTH_120p0_140p0,PTH_140p0_170p0,PTH_170p0_200p0,PTH_200p0_250p0,PTH_250p0_350p0,PTH_350p0_450p0,PTH_450p0_10000p0" \
  --category 0

python3 extract_response_matrix.py \
  --datacard datacards/Datacard_PTH_2022_2023_2024.root \
  --out-prefix response_matrix_combined/response_matrix_PTH \
  --order "PTH_0p0_5p0,PTH_5p0_10p0,PTH_10p0_15p0,PTH_15p0_20p0,PTH_20p0_25p0,PTH_25p0_30p0,PTH_30p0_35p0,PTH_35p0_45p0,PTH_45p0_60p0,PTH_60p0_80p0,PTH_80p0_100p0,PTH_100p0_120p0,PTH_120p0_140p0,PTH_140p0_170p0,PTH_170p0_200p0,PTH_200p0_250p0,PTH_250p0_350p0,PTH_350p0_450p0,PTH_450p0_10000p0" \
  --category 1

python3 extract_response_matrix.py \
  --datacard datacards/Datacard_PTH_2022_2023_2024.root \
  --out-prefix response_matrix_combined/response_matrix_PTH \
  --order "PTH_0p0_5p0,PTH_5p0_10p0,PTH_10p0_15p0,PTH_15p0_20p0,PTH_20p0_25p0,PTH_25p0_30p0,PTH_30p0_35p0,PTH_35p0_45p0,PTH_45p0_60p0,PTH_60p0_80p0,PTH_80p0_100p0,PTH_100p0_120p0,PTH_120p0_140p0,PTH_140p0_170p0,PTH_170p0_200p0,PTH_200p0_250p0,PTH_250p0_350p0,PTH_350p0_450p0,PTH_450p0_10000p0" \
  --category 2



python3 extract_response_matrix.py \
  --datacard datacards/Datacard_NJ_2022_2023_2024.root \
  --out-prefix response_matrix_combined/response_matrix_NJ \
  --order "NJ_0p0_1p0,NJ_1p0_2p0,NJ_2p0_3p0,NJ_3p0_4p0,NJ_4p0_100p0" \
  --category 0

python3 extract_response_matrix.py \
  --datacard datacards/Datacard_NJ_2022_2023_2024.root \
  --out-prefix response_matrix_combined/response_matrix_NJ \
  --order "NJ_0p0_1p0,NJ_1p0_2p0,NJ_2p0_3p0,NJ_3p0_4p0,NJ_4p0_100p0" \
  --category 1

python3 extract_response_matrix.py \
  --datacard datacards/Datacard_NJ_2022_2023_2024.root \
  --out-prefix response_matrix_combined/response_matrix_NJ \
  --order "NJ_0p0_1p0,NJ_1p0_2p0,NJ_2p0_3p0,NJ_3p0_4p0,NJ_4p0_100p0" \
  --category 2



python3 extract_response_matrix.py \
  --datacard datacards/Datacard_rapidity_2022_2023_2024.root \
  --out-prefix response_matrix_combined/response_matrix_YH \
  --order "rapidity_0p0_0p15,rapidity_0p15_0p3,rapidity_0p3_0p45,rapidity_0p45_0p6,rapidity_0p6_0p75,rapidity_0p75_0p9,rapidity_0p9_1p2,rapidity_1p2_1p6,rapidity_1p6_2p0,rapidity_2p0_2p5" \
  --category 0

python3 extract_response_matrix.py \
  --datacard datacards/Datacard_rapidity_2022_2023_2024.root \
  --out-prefix response_matrix_combined/response_matrix_YH \
  --order "rapidity_0p0_0p15,rapidity_0p15_0p3,rapidity_0p3_0p45,rapidity_0p45_0p6,rapidity_0p6_0p75,rapidity_0p75_0p9,rapidity_0p9_1p2,rapidity_1p2_1p6,rapidity_1p6_2p0,rapidity_2p0_2p5" \
  --category 1

python3 extract_response_matrix.py \
  --datacard datacards/Datacard_rapidity_2022_2023_2024.root \
  --out-prefix response_matrix_combined/response_matrix_YH \
  --order "rapidity_0p0_0p15,rapidity_0p15_0p3,rapidity_0p3_0p45,rapidity_0p45_0p6,rapidity_0p6_0p75,rapidity_0p75_0p9,rapidity_0p9_1p2,rapidity_1p2_1p6,rapidity_1p6_2p0,rapidity_2p0_2p5" \
  --category 2



python3 extract_response_matrix.py \
  --datacard datacards/Datacard_PTJ0_2022_2023_2024.root \
  --out-prefix response_matrix_combined/response_matrix_PTJ0 \
  --order "PTJ0_m10000p0_30p0,PTJ0_30p0_40p0,PTJ0_40p0_55p0,PTJ0_55p0_75p0,PTJ0_75p0_95p0,PTJ0_95p0_120p0,PTJ0_120p0_150p0,PTJ0_150p0_200p0,PTJ0_200p0_10000p0" \
  --category 0

python3 extract_response_matrix.py \
  --datacard datacards/Datacard_PTJ0_2022_2023_2024.root \
  --out-prefix response_matrix_combined/response_matrix_PTJ0 \
  --order "PTJ0_m10000p0_30p0,PTJ0_30p0_40p0,PTJ0_40p0_55p0,PTJ0_55p0_75p0,PTJ0_75p0_95p0,PTJ0_95p0_120p0,PTJ0_120p0_150p0,PTJ0_150p0_200p0,PTJ0_200p0_10000p0" \
  --category 1

python3 extract_response_matrix.py \
  --datacard datacards/Datacard_PTJ0_2022_2023_2024.root \
  --out-prefix response_matrix_combined/response_matrix_PTJ0 \
  --order "PTJ0_m10000p0_30p0,PTJ0_30p0_40p0,PTJ0_40p0_55p0,PTJ0_55p0_75p0,PTJ0_75p0_95p0,PTJ0_95p0_120p0,PTJ0_120p0_150p0,PTJ0_150p0_200p0,PTJ0_200p0_10000p0" \
  --category 2

```

The response matrices can also be computed individually, although this is outdated and not used for the derivation of the EFT results. One need the signal models for it to work.

```
# 2022
python3 extract_response_matrix_individual.py \
  --signal-dir signal_models/2022/PTH \
  --out-prefix response_matrix_2022/response_matrix_PTH \
  --order "PTH_0p0_5p0,PTH_5p0_10p0,PTH_10p0_15p0,PTH_15p0_20p0,PTH_20p0_25p0,PTH_25p0_30p0,PTH_30p0_35p0,PTH_35p0_45p0,PTH_45p0_60p0,PTH_60p0_80p0,PTH_80p0_100p0,PTH_100p0_120p0,PTH_120p0_140p0,PTH_140p0_170p0,PTH_170p0_200p0,PTH_200p0_250p0,PTH_250p0_350p0,PTH_350p0_450p0" \
  --category 0

python3 extract_response_matrix_individual.py \
  --signal-dir signal_models/2022/PTH \
  --out-prefix response_matrix_2022/response_matrix_PTH \
  --order "PTH_0p0_5p0,PTH_5p0_10p0,PTH_10p0_15p0,PTH_15p0_20p0,PTH_20p0_25p0,PTH_25p0_30p0,PTH_30p0_35p0,PTH_35p0_45p0,PTH_45p0_60p0,PTH_60p0_80p0,PTH_80p0_100p0,PTH_100p0_120p0,PTH_120p0_140p0,PTH_140p0_170p0,PTH_170p0_200p0,PTH_200p0_250p0,PTH_250p0_350p0,PTH_350p0_450p0" \
  --category 1

python3 extract_response_matrix_individual.py \
  --signal-dir signal_models/2022/PTH \
  --out-prefix response_matrix_2022/response_matrix_PTH \
  --order "PTH_0p0_5p0,PTH_5p0_10p0,PTH_10p0_15p0,PTH_15p0_20p0,PTH_20p0_25p0,PTH_25p0_30p0,PTH_30p0_35p0,PTH_35p0_45p0,PTH_45p0_60p0,PTH_60p0_80p0,PTH_80p0_100p0,PTH_100p0_120p0,PTH_120p0_140p0,PTH_140p0_170p0,PTH_170p0_200p0,PTH_200p0_250p0,PTH_250p0_350p0,PTH_350p0_450p0,PTH_450p0_10000p0" \
  --category 2

python3 extract_response_matrix_individual.py \
  --signal-dir signal_models/2022/NJ \
  --out-prefix response_matrix_2022/response_matrix_NJ \
  --order "NJ_0p0_1p0,NJ_1p0_2p0,NJ_2p0_3p0,NJ_3p0_4p0,NJ_4p0_100p0" \
  --category 0

python3 extract_response_matrix_individual.py \
  --signal-dir signal_models/2022/NJ \
  --out-prefix response_matrix_2022/response_matrix_NJ \
  --order "NJ_0p0_1p0,NJ_1p0_2p0,NJ_2p0_3p0,NJ_3p0_4p0,NJ_4p0_100p0" \
  --category 1

python3 extract_response_matrix_individual.py \
  --signal-dir signal_models/2022/NJ \
  --out-prefix response_matrix_2022/response_matrix_NJ \
  --order "NJ_0p0_1p0,NJ_1p0_2p0,NJ_2p0_3p0,NJ_3p0_4p0,NJ_4p0_100p0" \
  --category 2

python3 extract_response_matrix_individual.py \
  --signal-dir signal_models/2022/YH \
  --out-prefix response_matrix_2022/response_matrix_YH \
  --order "rapidity_0p0_0p15,rapidity_0p15_0p3,rapidity_0p3_0p45,rapidity_0p45_0p6,rapidity_0p6_0p75,rapidity_0p75_0p9,rapidity_0p9_1p2,rapidity_1p2_1p6,rapidity_1p6_2p0" \
  --category 0

python3 extract_response_matrix_individual.py \
  --signal-dir signal_models/2022/YH \
  --out-prefix response_matrix_2022/response_matrix_YH \
  --order "rapidity_0p0_0p15,rapidity_0p15_0p3,rapidity_0p3_0p45,rapidity_0p45_0p6,rapidity_0p6_0p75,rapidity_0p75_0p9,rapidity_0p9_1p2,rapidity_1p2_1p6,rapidity_1p6_2p0" \
  --category 1

python3 extract_response_matrix_individual.py \
  --signal-dir signal_models/2022/YH \
  --out-prefix response_matrix_2022/response_matrix_YH \
  --order "rapidity_0p0_0p15,rapidity_0p15_0p3,rapidity_0p3_0p45,rapidity_0p45_0p6,rapidity_0p6_0p75,rapidity_0p75_0p9,rapidity_0p9_1p2,rapidity_1p2_1p6,rapidity_1p6_2p0,rapidity_2p0_2p5" \
  --category 2

python3 extract_response_matrix_individual.py \
  --signal-dir signal_models/2022/PTJ0 \
  --out-prefix response_matrix_2022/response_matrix_PTJ0 \
  --order "PTJ0_m10000p0_30p0,PTJ0_30p0_40p0,PTJ0_40p0_55p0,PTJ0_55p0_75p0,PTJ0_75p0_95p0,PTJ0_95p0_120p0,PTJ0_120p0_150p0,PTJ0_150p0_200p0,PTJ0_200p0_10000p0" \
  --category 0

python3 extract_response_matrix_individual.py \
  --signal-dir signal_models/2022/PTJ0 \
  --out-prefix response_matrix_2022/response_matrix_PTJ0 \
  --order "PTJ0_m10000p0_30p0,PTJ0_30p0_40p0,PTJ0_40p0_55p0,PTJ0_55p0_75p0,PTJ0_75p0_95p0,PTJ0_95p0_120p0,PTJ0_120p0_150p0,PTJ0_150p0_200p0,PTJ0_200p0_10000p0" \
  --category 1

python3 extract_response_matrix_individual.py \
  --signal-dir signal_models/2022/PTJ0 \
  --out-prefix response_matrix_2022/response_matrix_PTJ0 \
  --order "PTJ0_m10000p0_30p0,PTJ0_30p0_40p0,PTJ0_40p0_55p0,PTJ0_55p0_75p0,PTJ0_75p0_95p0,PTJ0_95p0_120p0,PTJ0_120p0_150p0,PTJ0_150p0_200p0,PTJ0_200p0_10000p0" \
  --category 2





# 2023
python3 extract_response_matrix_individual.py \
  --signal-dir signal_models/2023/PTH \
  --out-prefix response_matrix_2023/response_matrix_PTH \
  --order "PTH_0p0_5p0,PTH_5p0_10p0,PTH_10p0_15p0,PTH_15p0_20p0,PTH_20p0_25p0,PTH_25p0_30p0,PTH_30p0_35p0,PTH_35p0_45p0,PTH_45p0_60p0,PTH_60p0_80p0,PTH_80p0_100p0,PTH_100p0_120p0,PTH_120p0_140p0,PTH_140p0_170p0,PTH_170p0_200p0,PTH_200p0_250p0,PTH_250p0_350p0,PTH_350p0_450p0" \
  --category 0

python3 extract_response_matrix_individual.py \
  --signal-dir signal_models/2023/PTH \
  --out-prefix response_matrix_2023/response_matrix_PTH \
  --order "PTH_0p0_5p0,PTH_5p0_10p0,PTH_10p0_15p0,PTH_15p0_20p0,PTH_20p0_25p0,PTH_25p0_30p0,PTH_30p0_35p0,PTH_35p0_45p0,PTH_45p0_60p0,PTH_60p0_80p0,PTH_80p0_100p0,PTH_100p0_120p0,PTH_120p0_140p0,PTH_140p0_170p0,PTH_170p0_200p0,PTH_200p0_250p0,PTH_250p0_350p0,PTH_350p0_450p0" \
  --category 1

python3 extract_response_matrix_individual.py \
  --signal-dir signal_models/2023/PTH \
  --out-prefix response_matrix_2023/response_matrix_PTH \
  --order "PTH_0p0_5p0,PTH_5p0_10p0,PTH_10p0_15p0,PTH_15p0_20p0,PTH_20p0_25p0,PTH_25p0_30p0,PTH_30p0_35p0,PTH_35p0_45p0,PTH_45p0_60p0,PTH_60p0_80p0,PTH_80p0_100p0,PTH_100p0_120p0,PTH_120p0_140p0,PTH_140p0_170p0,PTH_170p0_200p0,PTH_200p0_250p0,PTH_250p0_350p0,PTH_350p0_450p0,PTH_450p0_10000p0" \
  --category 2

python3 extract_response_matrix_individual.py \
  --signal-dir signal_models/2023/NJ \
  --out-prefix response_matrix_2023/response_matrix_NJ \
  --order "NJ_0p0_1p0,NJ_1p0_2p0,NJ_2p0_3p0,NJ_3p0_4p0,NJ_4p0_100p0" \
  --category 0

python3 extract_response_matrix_individual.py \
  --signal-dir signal_models/2023/NJ \
  --out-prefix response_matrix_2023/response_matrix_NJ \
  --order "NJ_0p0_1p0,NJ_1p0_2p0,NJ_2p0_3p0,NJ_3p0_4p0,NJ_4p0_100p0" \
  --category 1

python3 extract_response_matrix_individual.py \
  --signal-dir signal_models/2023/NJ \
  --out-prefix response_matrix_2023/response_matrix_NJ \
  --order "NJ_0p0_1p0,NJ_1p0_2p0,NJ_2p0_3p0,NJ_3p0_4p0,NJ_4p0_100p0" \
  --category 2

python3 extract_response_matrix_individual.py \
  --signal-dir signal_models/2023/YH \
  --out-prefix response_matrix_2023/response_matrix_YH \
  --order "rapidity_0p0_0p15,rapidity_0p15_0p3,rapidity_0p3_0p45,rapidity_0p45_0p6,rapidity_0p6_0p75,rapidity_0p75_0p9,rapidity_0p9_1p2" \
  --category 0

python3 extract_response_matrix_individual.py \
  --signal-dir signal_models/2023/YH \
  --out-prefix response_matrix_2023/response_matrix_YH \
  --order "rapidity_0p0_0p15,rapidity_0p15_0p3,rapidity_0p3_0p45,rapidity_0p45_0p6,rapidity_0p6_0p75,rapidity_0p75_0p9,rapidity_0p9_1p2" \
  --category 1

python3 extract_response_matrix_individual.py \
  --signal-dir signal_models/2023/YH \
  --out-prefix response_matrix_2023/response_matrix_YH \
  --order "rapidity_0p0_0p15,rapidity_0p15_0p3,rapidity_0p3_0p45,rapidity_0p45_0p6,rapidity_0p6_0p75,rapidity_0p75_0p9,rapidity_0p9_1p2,rapidity_1p2_1p6,rapidity_1p6_2p0,rapidity_2p0_2p5" \
  --category 2

python3 extract_response_matrix_individual.py \
  --signal-dir signal_models/2023/PTJ0 \
  --out-prefix response_matrix_2023/response_matrix_PTJ0 \
  --order "PTJ0_m10000p0_30p0,PTJ0_30p0_40p0,PTJ0_40p0_55p0,PTJ0_55p0_75p0,PTJ0_75p0_95p0,PTJ0_95p0_120p0,PTJ0_120p0_150p0,PTJ0_150p0_200p0,PTJ0_200p0_10000p0" \
  --category 0

python3 extract_response_matrix_individual.py \
  --signal-dir signal_models/2023/PTJ0 \
  --out-prefix response_matrix_2023/response_matrix_PTJ0 \
  --order "PTJ0_m10000p0_30p0,PTJ0_30p0_40p0,PTJ0_40p0_55p0,PTJ0_55p0_75p0,PTJ0_75p0_95p0,PTJ0_95p0_120p0,PTJ0_120p0_150p0,PTJ0_150p0_200p0,PTJ0_200p0_10000p0" \
  --category 1

python3 extract_response_matrix_individual.py \
  --signal-dir signal_models/2023/PTJ0 \
  --out-prefix response_matrix_2023/response_matrix_PTJ0 \
  --order "PTJ0_m10000p0_30p0,PTJ0_30p0_40p0,PTJ0_40p0_55p0,PTJ0_55p0_75p0,PTJ0_75p0_95p0,PTJ0_95p0_120p0,PTJ0_120p0_150p0,PTJ0_150p0_200p0,PTJ0_200p0_10000p0" \
  --category 2




# 2024
python3 extract_response_matrix_individual.py \
  --signal-dir signal_models/2024/PTH \
  --out-prefix response_matrix_2024/response_matrix_PTH \
  --order "PTH_0p0_5p0,PTH_5p0_10p0,PTH_10p0_15p0,PTH_15p0_20p0,PTH_20p0_25p0,PTH_25p0_30p0,PTH_30p0_35p0,PTH_35p0_45p0,PTH_45p0_60p0,PTH_60p0_80p0,PTH_80p0_100p0,PTH_100p0_120p0,PTH_120p0_140p0,PTH_140p0_170p0,PTH_170p0_200p0,PTH_200p0_250p0,PTH_250p0_350p0,PTH_350p0_450p0,PTH_450p0_10000p0" \
  --category 0

python3 extract_response_matrix_individual.py \
  --signal-dir signal_models/2024/PTH \
  --out-prefix response_matrix_2024/response_matrix_PTH \
  --order "PTH_0p0_5p0,PTH_5p0_10p0,PTH_10p0_15p0,PTH_15p0_20p0,PTH_20p0_25p0,PTH_25p0_30p0,PTH_30p0_35p0,PTH_35p0_45p0,PTH_45p0_60p0,PTH_60p0_80p0,PTH_80p0_100p0,PTH_100p0_120p0,PTH_120p0_140p0,PTH_140p0_170p0,PTH_170p0_200p0,PTH_200p0_250p0,PTH_250p0_350p0,PTH_350p0_450p0,PTH_450p0_10000p0" \
  --category 1

python3 extract_response_matrix_individual.py \
  --signal-dir signal_models/2024/PTH \
  --out-prefix response_matrix_2024/response_matrix_PTH \
  --order "PTH_0p0_5p0,PTH_5p0_10p0,PTH_10p0_15p0,PTH_15p0_20p0,PTH_20p0_25p0,PTH_25p0_30p0,PTH_30p0_35p0,PTH_35p0_45p0,PTH_45p0_60p0,PTH_60p0_80p0,PTH_80p0_100p0,PTH_100p0_120p0,PTH_120p0_140p0,PTH_140p0_170p0,PTH_170p0_200p0,PTH_200p0_250p0,PTH_250p0_350p0,PTH_350p0_450p0,PTH_450p0_10000p0" \
  --category 2

python3 extract_response_matrix_individual.py \
  --signal-dir signal_models/2024/NJ \
  --out-prefix response_matrix_2024/response_matrix_NJ \
  --order "NJ_0p0_1p0,NJ_1p0_2p0,NJ_2p0_3p0,NJ_3p0_4p0,NJ_4p0_100p0" \
  --category 0

python3 extract_response_matrix_individual.py \
  --signal-dir signal_models/2024/NJ \
  --out-prefix response_matrix_2024/response_matrix_NJ \
  --order "NJ_0p0_1p0,NJ_1p0_2p0,NJ_2p0_3p0,NJ_3p0_4p0,NJ_4p0_100p0" \
  --category 1

python3 extract_response_matrix_individual.py \
  --signal-dir signal_models/2024/NJ \
  --out-prefix response_matrix_2024/response_matrix_NJ \
  --order "NJ_0p0_1p0,NJ_1p0_2p0,NJ_2p0_3p0,NJ_3p0_4p0,NJ_4p0_100p0" \
  --category 2

python3 extract_response_matrix_individual.py \
  --signal-dir signal_models/2024/YH \
  --out-prefix response_matrix_2024/response_matrix_YH \
  --order "rapidity_0p0_0p15,rapidity_0p15_0p3,rapidity_0p3_0p45,rapidity_0p45_0p6,rapidity_0p6_0p75,rapidity_0p75_0p9,rapidity_0p9_1p2" \
  --category 0

python3 extract_response_matrix_individual.py \
  --signal-dir signal_models/2024/YH \
  --out-prefix response_matrix_2024/response_matrix_YH \
  --order "rapidity_0p0_0p15,rapidity_0p15_0p3,rapidity_0p3_0p45,rapidity_0p45_0p6,rapidity_0p6_0p75,rapidity_0p75_0p9,rapidity_0p9_1p2" \
  --category 1

python3 extract_response_matrix_individual.py \
  --signal-dir signal_models/2024/YH \
  --out-prefix response_matrix_2024/response_matrix_YH \
  --order "rapidity_0p0_0p15,rapidity_0p15_0p3,rapidity_0p3_0p45,rapidity_0p45_0p6,rapidity_0p6_0p75,rapidity_0p75_0p9,rapidity_0p9_1p2,rapidity_1p2_1p6,rapidity_1p6_2p0,rapidity_2p0_2p5" \
  --category 2

python3 extract_response_matrix_individual.py \
  --signal-dir signal_models/2024/PTJ0 \
  --out-prefix response_matrix_2024/response_matrix_PTJ0 \
  --order "PTJ0_m10000p0_30p0,PTJ0_30p0_40p0,PTJ0_40p0_55p0,PTJ0_55p0_75p0,PTJ0_75p0_95p0,PTJ0_95p0_120p0,PTJ0_120p0_150p0,PTJ0_150p0_200p0,PTJ0_200p0_10000p0" \
  --category 0

python3 extract_response_matrix_individual.py \
  --signal-dir signal_models/2024/PTJ0 \
  --out-prefix response_matrix_2024/response_matrix_PTJ0 \
  --order "PTJ0_m10000p0_30p0,PTJ0_30p0_40p0,PTJ0_40p0_55p0,PTJ0_55p0_75p0,PTJ0_75p0_95p0,PTJ0_95p0_120p0,PTJ0_120p0_150p0,PTJ0_150p0_200p0,PTJ0_200p0_10000p0" \
  --category 1

python3 extract_response_matrix_individual.py \
  --signal-dir signal_models/2024/PTJ0 \
  --out-prefix response_matrix_2024/response_matrix_PTJ0 \
  --order "PTJ0_m10000p0_30p0,PTJ0_30p0_40p0,PTJ0_40p0_55p0,PTJ0_55p0_75p0,PTJ0_75p0_95p0,PTJ0_95p0_120p0,PTJ0_120p0_150p0,PTJ0_150p0_200p0,PTJ0_200p0_10000p0" \
  --category 2
```

