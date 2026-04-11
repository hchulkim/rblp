# Package index

## Problem Setup

Define model formulations and create estimation problems.

- [`blp_formulation()`](https://hchulkim.github.io/rblp/reference/blp_formulation.md)
  : Create a BLP formulation
- [`BLPFormulation`](https://hchulkim.github.io/rblp/reference/BLPFormulation.md)
  : BLP Model Formulation
- [`blp_problem()`](https://hchulkim.github.io/rblp/reference/blp_problem.md)
  : Create a BLP Problem
- [`BLPProblem`](https://hchulkim.github.io/rblp/reference/BLPProblem.md)
  : BLP Problem

## Estimation Configuration

Configure integration, iteration, and optimization for estimation.

- [`blp_integration()`](https://hchulkim.github.io/rblp/reference/blp_integration.md)
  : Create integration configuration
- [`BLPIntegration`](https://hchulkim.github.io/rblp/reference/BLPIntegration.md)
  : Integration Configuration
- [`blp_iteration()`](https://hchulkim.github.io/rblp/reference/blp_iteration.md)
  : Create iteration configuration
- [`BLPIteration`](https://hchulkim.github.io/rblp/reference/BLPIteration.md)
  : Iteration Configuration
- [`blp_optimization()`](https://hchulkim.github.io/rblp/reference/blp_optimization.md)
  : Create optimization configuration
- [`BLPOptimization`](https://hchulkim.github.io/rblp/reference/BLPOptimization.md)
  : Optimization Configuration

## Results

Post-estimation analysis: elasticities, surplus, mergers, and more.

- [`BLPResults`](https://hchulkim.github.io/rblp/reference/BLPResults.md)
  : BLP Estimation Results

## Simulation

Simulate equilibrium data for testing and Monte Carlo studies.

- [`blp_simulation()`](https://hchulkim.github.io/rblp/reference/blp_simulation.md)
  : Create a BLP Simulation
- [`BLPSimulation`](https://hchulkim.github.io/rblp/reference/BLPSimulation.md)
  : BLP Simulation
- [`BLPSimulationResults`](https://hchulkim.github.io/rblp/reference/BLPSimulationResults.md)
  : BLP Simulation Results
- [`build_id_data()`](https://hchulkim.github.io/rblp/reference/build_id_data.md)
  : Build Balanced ID Data

## Micro Moments

Match micro-level data moments to model predictions.

- [`micro_dataset()`](https://hchulkim.github.io/rblp/reference/micro_dataset.md)
  : Create a Micro Dataset
- [`MicroDataset`](https://hchulkim.github.io/rblp/reference/MicroDataset.md)
  : Micro Dataset
- [`micro_part()`](https://hchulkim.github.io/rblp/reference/micro_part.md)
  : Create a Micro Part
- [`MicroPart`](https://hchulkim.github.io/rblp/reference/MicroPart.md)
  : Micro Part
- [`micro_moment()`](https://hchulkim.github.io/rblp/reference/micro_moment.md)
  : Create a Micro Moment
- [`MicroMoment`](https://hchulkim.github.io/rblp/reference/MicroMoment.md)
  : Micro Moment

## Instruments and Data Construction

Construct instruments and generate panel structures.

- [`build_blp_instruments()`](https://hchulkim.github.io/rblp/reference/build_blp_instruments.md)
  : Build BLP Instruments (Sums of Characteristics)
- [`build_differentiation_instruments()`](https://hchulkim.github.io/rblp/reference/build_differentiation_instruments.md)
  : Build Differentiation Instruments (Gandhi & Houde 2020)
- [`build_custom_ownership()`](https://hchulkim.github.io/rblp/reference/build_custom_ownership.md)
  : Build Custom Ownership Matrix

## Bundled Datasets

Bundled empirical datasets for replicating published results. The Nevo
(2000) cereal data and BLP (1995) automobile data are the standard
benchmarks in the differentiated products literature.

- [`load_nevo_products()`](https://hchulkim.github.io/rblp/reference/load_nevo_products.md)
  : Load Nevo (2000) cereal product data
- [`load_nevo_agents()`](https://hchulkim.github.io/rblp/reference/load_nevo_agents.md)
  : Load Nevo (2000) cereal agent data
- [`load_blp_products()`](https://hchulkim.github.io/rblp/reference/load_blp_products.md)
  : Load BLP (1995) automobile product data
- [`load_blp_agents()`](https://hchulkim.github.io/rblp/reference/load_blp_agents.md)
  : Load BLP (1995) automobile agent data
- [`nevo_products`](https://hchulkim.github.io/rblp/reference/nevo_products.md)
  : Nevo (2000) Cereal Product Data
- [`nevo_agents`](https://hchulkim.github.io/rblp/reference/nevo_agents.md)
  : Nevo (2000) Cereal Agent Data
- [`blp_products`](https://hchulkim.github.io/rblp/reference/blp_products.md)
  : BLP (1995) Automobile Product Data
- [`blp_agents`](https://hchulkim.github.io/rblp/reference/blp_agents.md)
  : BLP (1995) Automobile Agent Data
- [`load_mixtape_products()`](https://hchulkim.github.io/rblp/reference/load_mixtape_products.md)
  : Load Mixtape Sessions cereal product data
- [`load_mixtape_demographics()`](https://hchulkim.github.io/rblp/reference/load_mixtape_demographics.md)
  : Load Mixtape Sessions cereal demographic data
- [`prepare_mixtape_data()`](https://hchulkim.github.io/rblp/reference/prepare_mixtape_data.md)
  : Prepare Mixtape product data for rblp estimation
- [`mixtape_products`](https://hchulkim.github.io/rblp/reference/mixtape_products.md)
  : Mixtape Sessions Cereal Product Data
- [`mixtape_demographics`](https://hchulkim.github.io/rblp/reference/mixtape_demographics.md)
  : Mixtape Sessions Cereal Demographic Data
- [`load_otc_products()`](https://hchulkim.github.io/rblp/reference/load_otc_products.md)
  : OTC Drug Product Data (baby_BLP)

## Options

Package-wide configuration options.

- [`rblp_options()`](https://hchulkim.github.io/rblp/reference/rblp_options.md)
  : Get or Set rblp Package Options

## Package

- [`rblp-package`](https://hchulkim.github.io/rblp/reference/rblp-package.md)
  [`rblp`](https://hchulkim.github.io/rblp/reference/rblp-package.md) :
  rblp: BLP Demand Estimation for Differentiated Products
