<!-- ABOUT THE PROJECT -->
## About The Project
Framework to predict and to metabotype Omics time series data from several diet interventions described in depth in our article: Skantze et al. "Data-driven analysis and prediction of dynamic postprandial metabolic response to multiple dietary challenges using dynamic mode decomposition" 
https://www.frontiersin.org/articles/10.3389/fnut.2023.1304540/full


<!-- USAGE EXAMPLES -->
## Usage
Run metabotyping_using_pDMDc.m and prediction_using_pDMDc.m to see examples for metabotyping and prediction of metabolomics time series. prediction_using_pDMDc.m is showing the prediction performance when adding more diets in the training set. The function pDMDcPredictionValidation.m can be used for training and prediction on the validation set. The function pDMDcPredictionTest.m can be used for prediction of time series using the trained parameters A and B from pDMDcPredictionValidation.m. The data needs to be structured in tensor/multidimensional array format and organized as described in the example files. The simulated data available is generated from the human metabolic model: Hiroyuki Kurata, Virtual metabolic human dynamic model for pathological analysis and therapy design for diabetes, https://www.sciencedirect.com/science/article/pii/S2589004221000699
<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE.txt` for more information.

<p align="right">(<a href="#readme-top">back to top</a>)</p>
