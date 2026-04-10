# Matrix_ATLAS
This is a repo to test the approach used by ATLAS to extract the correlation matrix among various 1D differential XS. The approach basically relies on summing two covariance matrices:

$$
C = C_{\rm stat} + C_{\rm syst}
$$

The datacards used for the HIG-23-014 paper are needed, they can be downloaded from:
```
git clone https://gitlab.cern.ch/cms-analysis/hig/hig-23-014/datacards.git
```

## Statistical covariance matrix $C_{\rm stat}$
This is computed starting from counting the events in sidebands in bootstrap replica $C_{\rm stat}^{\rm sb}$. The statistical covariance matrix is obtained through the formula:

$$
C_{\rm stat} = R^{-1}~C_{\rm stat}^{\rm sb}~(R^{-1})^T
$$

where $R$ is the reponse matrix.

The steps to extract this matrix are encoded in the following scripts:
* `bootstrap.py`: this script produces bootstrap replicas of the dataset and count events in the sidebands
* `extract_response_matrix.py`: this script extracts the $(\epsilon\times A)$ matrix from the workspaces
  * TODO: I am extracting now $(\epsilon\times A)$, I should only look at $\epsilon$. My plan is to take the fiducial acceptances $A$ and perform the division $(\epsilon\times A)/A$
* `correlation_matrix.ipynb`: this notebook computes $C_{\rm stat}^{sb}$ and $C_{\rm stat}$


List of commands:

```
python3 bootstrap.py --n-replicas 1000 --output bootstrap/bootstrap_sideband_counts.txt
```

```
python3 extract_response_matrix.py \
  --signal-dir datacards/differentials/Models_PTH/signal \
  --out-prefix response_matrix/response_matrix_pth \
  --order "PTH_0p0_15p0,PTH_15p0_30p0,PTH_30p0_45p0,PTH_45p0_80p0,PTH_80p0_120p0,PTH_120p0_200p0,PTH_200p0_350p0,PTH_350p0_10000p0" 

python3 extract_response_matrix.py \
  --signal-dir datacards/differentials/Models_Njets2p5/signal \
  --out-prefix response_matrix/response_matrix_Njets2p5 \
  --order "Njets2p5_0p0_1p0,Njets2p5_1p0_2p0,Njets2p5_2p0_3p0,Njets2p5_3p0_100p0" 

python3 extract_response_matrix.py \
  --signal-dir datacards/differentials/Models_ptJ0/signal \
  --out-prefix response_matrix/response_matrix_ptJ0 \
  --order "first_jet_pt_0p0_30p0,first_jet_pt_30p0_75p0,first_jet_pt_75p0_120p0,first_jet_pt_120p0_200p0,first_jet_pt_200p0_10000p0" 
```