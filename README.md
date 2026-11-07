# ZIO: Order-based Structure Learning for Zero-Inflated Count Data under Data Heterogeneity

This repository provides the R implementation of **ZIO** and code for the simulation, heterogeneity, and real-data experiments. ZIO estimates directed acyclic graphs (DAGs) from zero-inflated count data using an order-based search with stepwise parent selection. The repository also contains competing-method implementations and supplementary experiments.

## Quick start

**Run the commands below from the `ZIO/` repository root.** For example, after cloning the repository, run `cd ZIO`. 
### Requirements

- **R** with the packages: `pscl`, `igraph`, `MASS`, `foreach`, `doParallel`, `doRNG`, `dplyr`, `ggplot2`, `tidyr`, `tibble`, `glmnet`, `MXM`, `ggraph`, `tidygraph`, `viridis`, `ggbump`, `shadowtext`, and `readr`. The standard `parallel` and `grid` packages ship with R.
- The **Python benchmark notebooks** additionally require the packages in `requirements-python.txt`. Use `pip install -r requirements-python.txt` to in a compatible Python environment.

For example, install the R dependencies with:

```r
install.packages(c("pscl", "igraph", "MASS", "foreach", "doParallel", "doRNG",
                   "dplyr", "ggplot2", "tidyr", "tibble", "glmnet", "MXM",
                   "ggraph", "tidygraph", "viridis", "ggbump", "shadowtext", "readr"))
```

### 1. Main ZIO simulation — multiple datasets

The main implementation learns a shared ordering across `K` heterogeneous datasets and estimates a DAG for each dataset. Example:

```bash
Rscript code/simulation/4.1/Simulation_ZIO.R 25 1 1 1000 4
```

**Arguments (in order):**

- `p`: number of variables (here, `25`).
- `K`: number of datasets (here, `1`).
- `seed`: random seed (here, `1`).
- `N` *(optional)*: number of order-search iterations (here, `1000`).
- `workers` *(optional)*: number of parallel workers requested (here, `4`).
- `checkpoint_every`, `leave_free`, and `output_dir` are further optional arguments; see the script header.

Simulation results are written to `output/simulation/4.1/final_simulation_results/`. Each seed produces a result containing `ZIO_DAG`.

### 2. Zero-inflation proportion experiment

This experiment evaluates the method under different zero proportions:

```bash
Rscript code/simulation/4.2/Simulation_ZIO_zeros.R 1000 20 1 '50%' 1
```

**Arguments (in order):**

- `n`: sample size (here, `1000`).
- `p`: number of variables (here, `20`).
- `K`: number of datasets (here, `1`).
- `zeros`: target proportion of zeros (e.g., `25%`, `50%`, `75%`).
- `seed`: random seed (here, `1`).

Results are written under `output/simulation/zero_proportions/`.

### 3. Real-data analysis — joint order search

The real-data script reads `data/control_list_41.rds` and `data/suicide_list_41.rds` and learns across the selected subjects:

```bash
Rscript code/real_data/5/ZIO_realdata_joint.R
```

The script currently uses `p = 41`, `N = 3000` and excludes the subject at original combined index 14. Its output is `output/ZIO_DAG_P_41_K_<K>.RData`, where `<K>` is the number of subjects remaining after exclusion. The final DAGs are stored in `ZIO_DAG$E_estimate`.

For **individual-subject** analysis, pass the dataset index:

```bash
Rscript code/supplementary/C2/ZIO_realdata_individual.R 1
```

Run with the appropriate index for each subject. Individual outputs are saved as `output/ZIO_DAG_P_41_K_1_dataset_XX.RData`.

### 4. Operator comparison — R2R, RTS and ADJ

These scripts compare different order-search neighborhoods on fixed simulated data. Examples:

```bash
Rscript code/supplementary/B2/Simulation_ZIO_fixed_data.R 10 1 1 3000 4
Rscript code/supplementary/B2/Simulation_ZIO_fixed_data_RTS.R 10 1 1 3000 4
Rscript code/supplementary/B2/Simulation_ZIO_fixed_data_ADJ.R 10 1 1 3333 4
```

**Arguments:** `p`, `K`, `seed`, then optionally `N`, `workers`, `checkpoint_every`, `leave_free`, and `output_dir`. The three scripts write to separate subfolders of `output/simulation/operator_comparison/`.

### 5. Shared-order sensitivity — fixed Kendall's tau

This experiment varies ordering similarity across datasets:

```bash
Rscript code/supplementary/B3/Simulation_ZIO_fixed_tau.R 1000 10 5 1 0.5 500 4
```

**Arguments (in order):** `n`, `p`, `K`, `seed`, `tau_target`, then optionally `N`, `workers`, `checkpoint_every`, `leave_free`, and `output_dir`. Results are written under `output/simulation/fixed_tau/`.

### 6. Accelerated parent selection — shortlist

This experiment uses the conditional ZIP shortlist for stepwise parent selection:

```bash
Rscript code/supplementary/B4/Simulation_ZIO_shortlist.R 25 1 1 1250 4
```

**Arguments:** `p`, `K`, `seed`, then optionally `N`, `workers`, `checkpoint_every`, `leave_free`, and `output_dir`. Optional environment variables configure the shortlist: `ZIPBN_PARENT_SELECTION` (default `shortlist`), `ZIPBN_SHORTLIST_B` (default `10`), and `ZIPBN_DMAX` (default `5`). Results are written under `output/simulation/shortlist/`.

### 7. Model misspecification — ZIP versus ZINB fitting

To fit a negative-binomial zero-inflated model to an existing, saved negative-binomial-generated dataset:

```bash
Rscript code/supplementary/B3/Simulation_misspecification_parallal_seed.R 1000 25 1 negbin 1 negbin
```

**Arguments:** `n`, `p`, `K`, fitted distribution (`poisson` or `negbin`), `seed`, and optional data-generating distribution (`poisson` or `negbin`). The script **reuses the seed-specific data** saved in a corresponding Poisson-fit `.RData` file, so that fitted models are compared on the same observations. These source results must exist before rerunning fits. Outputs use `output/simulation/model_misspecification/Model_misspecification_<data_model>_data/`.

## Repository layout

```text
ZIO/
├── code/
│   ├── helper/                 # Shared R helpers and Poisson NOTEARS functions
│   ├── simulation/             # Main experiments and competing methods
│   ├── real_data/              # Joint real-data estimation and summaries
│   └── supplementary/          # Additional numerical experiments
├── data/                       # Input data, including simulation_data_{25,50,75}
├── output/                     # Saved results and newly generated outputs
├── requirements-python.txt
└── README.md
```

The supplementary code includes numerical order-separation experiments, benchmark evaluations, and figure/table scripts. **This README focuses on running the estimation and simulation programs**, rather than detailing every plot and table command.

## Notes

- R scripts should be invoked from `ZIO/` or its parent directory. The Python benchmark notebooks use the same repository-relative paths.
- The `code/simulation/4.1/notears-master/` folder contains third-party NOTEARS code and its original license. The `ZICO-main/` folder also contains third-party code; review its licensing and attribution before making the repository public.
- Data redistribution rights (particularly for real-data files) should be checked before public release.
- Parallel performance depends on the machine and the requested number of workers. Some full experiments are computationally intensive.
