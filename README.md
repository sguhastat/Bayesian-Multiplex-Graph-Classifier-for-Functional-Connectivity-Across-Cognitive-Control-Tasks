# Bayesian-Multiplex-Graph-Classifier-for-Functional-Connectivity-Across-Cognitive-Control-Tasks

A Bayesian Multiplex Graph Classifier of Functional Brain Connectivity Across Diverse Tasks of Cognitive Control (Jose Rodriguez-Acosta, Sharmistha Guha, Ivo D. Dinov)

Neuroinformatics, 2024 

(Special Issue: Data Science Methods and Neuroinformatics Applications)

We investigate the impact of aging on functional brain connectivity across cognitive control tasks, with a focus on identifying brain regions linked to early aging. Modeling functional connectivity as multiplex graphs, we address the challenge of predicting a binary outcome (aging vs. normal) using multiple graph predictors. Existing methods often struggle to fully leverage within- and across-layer information, particularly in small samples. To overcome this, we propose the Bayesian Multiplex Graph Classifier (BMGC), which models edge coefficients via bilinear interactions of node-specific latent effects and applies variable selection to identify influential nodes. BMGC offers computational efficiency and uncertainty quantification in node identification, coefficient estimation, and prediction. Simulations show superior performance over existing methods, and application to fMRI data reveals symmetric patterns in the sensory motor network and asymmetric aging-related effects in the default mode network.

Repository for the classification model outlined in: A Bayesian Multiplex Graph Classifier of Functional Brain Connectivity Across Diverse Tasks of Cognitive Control.
Refer to: https://github.com/jeroda7105/Classification-with-Multi-Layer-Graphs/tree/main


binary_multi_network_sim_data_gen.R contains a function for simulating the data as in the paper.
binary_multi_network_function_with_results.R contains a function for running the model on supplied data and obtaining results.
binary_multi_network_model_example.R contains an example of generating simulated data, running the classification model, and observing the results.

