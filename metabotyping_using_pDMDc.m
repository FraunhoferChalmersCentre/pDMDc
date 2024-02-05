%% Plug and play script for identifying metabotypes using pDMDc
%
% Example: 
%    run pDMDc_prediction_example.m
%
% Other m-files required: The Nway package by Rasmus Bro, pDMDcMetabotyping.m & ClusteringFunction.m
% Nway can be installed from https://se.mathworks.com/matlabcentral/fileexchange/1088-the-n-way-toolbox
% Subfunctions: none
% MAT-files required: "100_diets_simulated_from_herring.mat"
%
% Author: Viktor Skantze
% Work address: 
% email: viktor.skantze@fcc.chalmers.se
%------------- BEGIN CODE --------------
clc; clear all
load("Data\data_3_diets.mat")
data_metabolomics = data_collection{1};
input_meal_all =  data_collection{2};
cluster_ground_truth = data_collection{3};
metabolite_list = data_collection{4};
%% Choosing model complexity
data = permute(data_metabolomics, [4 3 2 1]);
nr_diets = size(data,4);
nr_individuals = size(data,1);
nr_time_points = size(data,2);

X = [];
for i = 1:nr_individuals
    for d = 1:nr_diets
        tmp_dataset = squeeze(data(i,:,:,d));
        X = [X, tmp_dataset'];
    end
end

[~,S,~] = svd(X,'econ');
s = diag(S);
plot(s(1:10))
xlabel('States')
ylabel('Singular values')
title('Scree plot')
%% Metabotyping using pDMDc
nr_of_components = 4;                                                      % Number of components in pDMDc
data_metabolomics_dmd = nprocess(data_metabolomics, [0 0 0 0], [0 1 0 0],[],[],1);
data_metabolomics_dmd = data_metabolomics_dmd(:,~isnan(std(data_metabolomics_dmd,[],[1 3 4])),:,:);
%data_metabolomics_dmd = data_metabolomics;
[A, B, prediction_train, X_states_tensor, U] = pDMDcMetabotyping(data_metabolomics_dmd, input_meal_all,nr_of_components); % Train and predict
%% PARAFAC/CANDECOMP (CP)
data_metabolomics_parafac = data_metabolomics + data_metabolomics.*0.05.*rand(size(data_metabolomics));
data_metabolomics_parafac = nprocess(data_metabolomics_parafac, [0 0 0 0], [0 1 0 0],[],[],1);
data_metabolomics_parafac = data_metabolomics_parafac(:,~isnan(std(data_metabolomics_parafac,[],[1 3 4])),:,:);
parafac_components = 3;
model_parafac = parafac(data_metabolomics_parafac, parafac_components);
%% Classifying metabolomics
[basis, classification_dynamics] = ClassificationDynamics(data_metabolomics_parafac, 4);
%% Plot clusters using DMD, CP, and ground truth data
%metabolite_names = {'Glygn_H', 'Gap_H', 'Ala_H', 'Accoa_HM','Bhb_H', 'Accoa_HC'};
metabolite_names = {'$Glucose_p$', '$Glycogen_h$', '$Glucose_l$', ...
    '$\beta-Hydroxybutyrate_h$', ...
    '$Alanine_h$','$NADP+_{GI}$'};
list_amino_acids = [2 147 11 141 145 136];%[147 148 141 142 145 146];
state_of_interest = 4;
diet_of_interest = 1;  
denom = 6;
PlotClustering(cluster_ground_truth, ...
    X_states_tensor, ...
    model_parafac, ...
    data_metabolomics_dmd, ...
    list_amino_acids, ...
    metabolite_names, ...
    state_of_interest, ...
    diet_of_interest, ...
    denom)
%% Plot metabolites 
% metabolite_list = struct2cell(load("Data\metabolite_names.mat"));
% metabolite_list = metabolite_list{1};
% metabolite_list = metabolite_list(~idx_remove);
classification_description_metabolites = {'Component 1','Component 2','Component 3','Component 4'};
PlotMetaboliteLoadings(metabolite_list, U,data_metabolomics_dmd, ...
    model_parafac{2},[],classification_description_metabolites, ...
    classification_description_metabolites)
%% Plot dynamic states of DMD, CP and the clustered data
list_indices = [4 1 3 2];
order_cp = [1 2 3];
classification_description_dynamics = {'Cluster 1', 'Cluster 2', 'Cluster 3', 'Cluster 4'};
PlotDynamicStates(nr_of_components, ...
    classification_dynamics, ...
    classification_description_dynamics, ...
    data_metabolomics_dmd, ...
    X_states_tensor, ...
    model_parafac{3}, ...
    list_indices, ...
    order_cp, ...
    denom)
%------------- END OF CODE --------------