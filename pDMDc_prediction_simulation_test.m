%% Simulation study of predicting the response of a new diet using pDMDc
%
% Example: 
%    run dmdc_prediction.m
%
% Other m-files required: none
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
data_included = data_simulated{1}();
input_meal_all = std(data_included,0,'all')*data_simulated{2}(:,1:3)'; % scaling input to account for more variance in svd
%% Set up scaling loop
scaling_factor = [0,0.1,0.2,0.3]; % each noise level per metabolite standard deviation
cosine_similairty_tot = zeros(4,100,50);
rank_best_save = zeros(4,100,50);
for iter_scaling = 1:4
    error_per_diet = zeros(100,50,79)*nan;
    cosine_similairty = zeros(100,50)*nan;
    for iter_nr_diets = 1:50 % training set is ranging from 1 to 50 diets
        for i_repititions = 1:100 % sampling new data 100 times
            %% Sampling data, random diet and individual
            individual_sampling = randi(17)*ones(100,1);
            diet_sampling = randi(100,100,1);

            sampled_data = zeros(100,79,8);
            for i = 1:100
                sampled_data(i,:,:) = squeeze( data_included(diet_sampling(i),:,:,individual_sampling(i)) );
            end

            nr_diets = iter_nr_diets;
            data_individual = sampled_data;
            data_individual = permute(data_individual,[1 3 2]);
            %% Add noise
            std_per_metabolite = scaling_factor(iter_scaling)*squeeze(std(data_individual,1, [ 1 2]));
            data_individual_noise = zeros(size(data_individual));
            for iter_noise = 1:100
                data_individual_noise(iter_noise, :,:) = squeeze( data_individual(iter_noise, :,:)) + (mvnrnd(zeros(8,1),eye(8),79).*std_per_metabolite)';
            end
            %% Set up data
            nr_diets_training = 1:nr_diets;
            diets_validation = nr_diets+1;
            data_simulated_training = data_individual_noise(nr_diets_training,:,:);
            data_simulated_training = permute(data_simulated_training, [1,3,2]);
            data_simulated_validation = data_individual_noise(diets_validation,:,:);
            data_simulated_validation = squeeze( permute(data_simulated_validation, [1,3,2]) );
            data_simulated_validation_noiseless = data_individual(diets_validation,:,:);
            data_simulated_validation_noiseless = squeeze(permute(data_simulated_validation_noiseless, [1,3,2]));

            input_meal_validation = input_meal_all(:,diets_validation);
            input_meal_training = input_meal_all(:,nr_diets_training);
            
            %% Set up the DMD data format
            [A, B, prediction_val_best, ~, rank_best, rank_max] = pDMDcPrediction(data_simulated_training, data_simulated_validation, input_meal_training, input_meal_validation);
            %% Calculate cosine similarity
            cosine_similairty(i_repititions,iter_nr_diets) = CosineTensor(data_simulated_validation, prediction_val_best);
            rank_best_save(iter_scaling,i_repititions, iter_nr_diets) = 100*(rank_best/rank_max);
        end
    end
    cosine_similairty_tot(iter_scaling,:,:) = cosine_similairty;
    disp(iter_scaling) % Counting up to 100 when its finished
end
%% Save result
save('cosine_similairty_tot.mat','cosine_similairty_tot')
save('rank_best_save.mat','rank_best_save')
%% Plot saved result
load('cosine_similairty_tot.mat')
load('rank_best_save.mat')
scaling_factor = [0,0.1,0.2,0.3];
for iter_scaling = 1:4
    subplot(2,2,iter_scaling)
    yyaxis left
    boxplot(squeeze( cosine_similairty_tot(iter_scaling,:,:) ),'PlotStyle','compact')
    ylim([min(cosine_similairty_tot(:)), 1.1])
    ylabel('Cosine similarity')

    yyaxis right
    tmp_rank = mean(squeeze(rank_best_save(iter_scaling,:,:)),1);
    plot(1:50, tmp_rank,'r-', LineWidth=1)
    ylabel('% of full rank')
    ylim([0, 100])
    if iter_scaling>2
        xticks(1:7:50)
        xticklabels(1:7:50)
        xlabel('# diets in training set')
    else
        xticks([])
        xticklabels([])
    end
    title(['Noise level: \sigma_m*',num2str(scaling_factor(iter_scaling))])
end
%% Function for cosine similarity
function cosine_sim_mean = CosineTensor(pred, data)

nr_metabolites = size(data,3);
cosine_sim = zeros(1, nr_metabolites);
data = data';
pred = pred';
for m = 1:nr_metabolites
    current_data = normalize( squeeze( data(:,m) )');
    current_pred = normalize( squeeze( pred(:,m) )');

    norm_pred = sqrt(current_pred*current_pred');
    norm_data = sqrt(current_data*current_data');
    cosine_sim(1,m) = current_data*current_pred'/(norm_data*norm_pred);

end
cosine_sim_mean = mean(cosine_sim,"all");
end
%------------- END OF CODE --------------