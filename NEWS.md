# scClustEval 1.0.0

## Initial Release

* First public release of scClustEval
* Core features:
  - `sc_assessment()` / `self_projection()` for clustering quality assessment
  - `sc_optimize()` / `sc_optimize_all()` for iterative clustering optimization
  - Support for multiple classifiers: LR, RF, SVM, NB, DT, XGBoost, Ranger
  - Seamless Seurat integration with `RunAssessment()` and `RunOptimization()`
  - Visualization functions: ROC curves, confusion heatmaps, Sankey diagrams
  - C++ acceleration for distance computation and matrix operations
  - Full compatibility with Seurat V4 and V5
  - Cross-platform support (Windows, macOS, Linux)
