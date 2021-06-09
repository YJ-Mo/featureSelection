# featureSelection

## Edited by Y.MO, last updates on 2021/06/09

## Main file descriptions:

1. normalize.R/normalize_server.R: Normalize count matrix by cross-check sample and gene name, sparsity filtering and RUVSeq normalize.

2. LR_discovery.R: Main script to 1) partition discovery count matrix, 2) cross validate with DE genes, 3) shuffle and iterate, 4) fit model to validation count matrix and 5) output inner AUC and external AUC plots.

3. Classifiers.R: Sealed functions apploed in main script.

## Ackownledgement

1. Thanks to my lab fellow @uaauaguga for normalization and LR ideas.

2. Thanks to my lab fellow @QingZhan98 for discussion.
