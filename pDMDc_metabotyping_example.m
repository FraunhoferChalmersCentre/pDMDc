%% Plug and play script for identifying metabotypes using pDMDc
%
% Example: 
%    run pDMDc_prediction_example.m
%
% Other m-files required: pDMDcMetabotyping.m & ClusteringFunction.m
% Subfunctions: none
% MAT-files required: "100_diets_simulated_from_herring.mat"
%
% Author: Viktor Skantze
% Work address: 
% email: viktor.skantze@fcc.chalmers.se
%------------- BEGIN CODE --------------
%% Load data
clc; clear all;close all
load("100_diets_simulated_from_herring.mat")
data_metabolomics = data_simulated{1};
input_meal_all = std(data_metabolomics,0,'all')*data_simulated{2}(:,1:3)'; % scaling input to account for more variance in svd
size(data_metabolomics)                                                    % 100 diets, 79 metabolites, 8 time points, 17 individuals
size(input_meal_all)                                                       % 3 macronutrients(fat, protein, carbohydrates) for 100 diets
%% Choose number of diets
nr_diets = 100;                                                            % Chosen number of diets to include in analysis
data_metabolomics = data_metabolomics(1:nr_diets,:,:,:);                   % Subsetting data
input_meal_all = input_meal_all(:,1:nr_diets);                             % Subsetting data
%% Metabotyping using pDMDc
nr_of_components = 5;                                                      % Number of components in pDMDc
[A, B, prediction_train, X_states_tensor] = pDMDcMetabotyping(data_metabolomics, input_meal_all,nr_of_components); % Train and predict
%% Clustering state trajectories
metric = 1;                                                                % 1: Cosine similarity, 2: Sums of squares
state_of_interest = 1;                                                     % State to cluster
diet_of_interest = 1;                                                      % Diet to cluster
nr_clusters = 2;                                                           % Number of clusters
clustering_states_distinct = ClusteringFunction(X_states_tensor, metric,state_of_interest, diet_of_interest, nr_clusters);
%------------- END OF CODE --------------