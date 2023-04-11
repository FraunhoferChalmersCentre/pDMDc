%% Plug and play script for predicting response of a new diet using pDMDc
%
% Example: 
%    run pDMDc_prediction_example.m
%
% Other m-files required: pDMDcPrediction.m
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
%% Predict using pDMDc
nr_diets_training = 3;                                                     % Change number of diets in training to see how the predictions change 
data_training = squeeze(data_metabolomics(1:nr_diets_training,:,:,1));     % 4 diets in the training set for individual 1
data_validation = squeeze(data_metabolomics(nr_diets_training+1,:,:,1));   % 1 diet in the validation set set for individual 1
input_meal_training = input_meal_all(:,1:nr_diets_training);                               % meal information for training set
input_meal_validation = input_meal_all(:,nr_diets_training+1);                              % meal information for validation set

[A, B, prediction_val_best, prediction_train_best] = pDMDcPrediction(data_training, data_validation, input_meal_training, input_meal_validation); % Train and predict
%% Plot training reconstruction
random_metabolites = randperm(79,10);
colors = lines();
iter = 0;
figure()
for m = random_metabolites
    iter = iter + 1;
    subplot(2,5,iter)
    for j = 1:nr_diets_training
        plot(1:8, squeeze(prediction_train_best(j,m,:)), Color= colors(j,:), LineWidth=1); hold on
        plot(1:8, squeeze(data_training(j,m,:)),'.-', Color= colors(j,:));
    end
    title(["Metabolite #",num2str(m)])
end
p = [];
legend_names = {};
for j = 1:nr_diets_training
    p(j) = plot(nan,nan,  Color= colors(j,:));
    legend_names{j} = strcat("Diet ", num2str(j) );
end
p(end+1) = plot(nan,nan,  'k-');
p(end+1) = plot(nan,nan,  'k.-');
legend_names{end+1} = "Prediction";
legend_names{end+1} = "Data";
legend(p, legend_names)
%% Plot validation prediction
iter = 0;
figure()
for m = random_metabolites
    iter = iter + 1;
    subplot(2,5,iter)
    plot(1:8, squeeze(prediction_val_best(m,:)), Color= colors(1,:)); hold on
    plot(1:8, squeeze(data_validation(m,:)),'.-', Color= colors(1,:)); hold on
    title(["Metabolite #",num2str(m)])
end

legend("Prediction", "Data")
%% Cosine score
prediction_val_tmp = zeros(1,79,8);
prediction_val_tmp(1,:,:) = prediction_val_best;
prediction_val_tmp = permute(prediction_val_tmp,[1,3,2]);
data_validation_tmp = zeros(1,79,8);
data_validation_tmp(1,:,:) = data_validation;
data_validation_tmp = permute(data_validation_tmp,[1,3,2]);

data_training_tmp = permute(data_training,[1,3,2]);
prediction_train_best_tmp = permute(prediction_train_best,[1,3,2]);

cosine_training = CosineTensor(prediction_train_best_tmp, data_training_tmp);
disp([strcat("Cosine similarity score for reconstruction and data: ", num2str(cosine_training))])


cosine_validation = CosineTensor(prediction_val_tmp, data_validation_tmp);
disp([strcat("Cosine similarity score for prediction and validation data: ", num2str(cosine_validation))])


function cosine_sim_mean = CosineTensor(pred, data)
nr_diets = size(data,1);
nr_metabolites = size(data,3);
cosine_sim = zeros(nr_diets, nr_metabolites);
for d = 1:nr_diets
    for m = 1:nr_metabolites
        current_data = normalize( squeeze( data(d,:,m) )');
        current_pred = normalize( squeeze( pred(d,:,m) )');

        norm_pred = sqrt(current_pred'*current_pred);
        norm_data = sqrt(current_data'*current_data);
        cosine_sim(d,m) = current_data'*current_pred/(norm_data*norm_pred);

    end
end
cosine_sim_mean = mean(cosine_sim,"all");
end
%------------- END OF CODE --------------