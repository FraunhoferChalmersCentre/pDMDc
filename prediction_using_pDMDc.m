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
%rng(123)
clear all;%close all

addpath("Functions")
load("Data\data_save.mat")
%% Preallocate

data_metabolomics = data_save{1};
data_metabolomics = permute(data_metabolomics,[4 2 3 1]);


idx_remove =  (std(data_metabolomics,[],[1 3 4 ],"omitnan") < 1e-3);
data_metabolomics = data_metabolomics(:,~idx_remove,:,:);
individuals_keep = sum( isnan(squeeze(data_metabolomics(:,1,1,:))), 1 ) == 0;
data_metabolomics = data_metabolomics(:,:,:,individuals_keep);


input_meal_all =  data_save{2};
size(data_metabolomics)                                                    % 3 diets, 79 metabolites, 8 time points, 17 individuals
size(input_meal_all)                                                       % 3 macronutrients(fat, protein, carbohydrates) for 100 diets
data_metabolomics = data_metabolomics(:,:,1:end,:);                        % Removing baseline
nr_time = size(data_metabolomics,3);
nr_metabolomics = size(data_metabolomics,2);
nr_individuals = size(data_metabolomics,4);
%% Test data
data_metabolomics_test = data_metabolomics(51:end,:,:,1:17);
nr_test_diets = size(data_metabolomics_test,1);
data_metabolomics_mix_test = zeros(nr_test_diets*17,nr_metabolomics, nr_time);
input_test= zeros( size(input_meal_all,1), nr_test_diets*17);
for inds = 1:17
    indx = ((inds-1)*(nr_test_diets-1))+(inds:inds+(nr_test_diets-1));
    data_metabolomics_mix_test(indx,:,:) = data_metabolomics_test(1:nr_test_diets,:,:,inds);
    input_test(:,indx) = input_meal_all(:,51:end);
end
data_test = data_metabolomics_mix_test;
for i = 1:size(data_test,1)
    data_test(i,:,:) = normalize( squeeze(data_test(i,:,:))' )';
end
%% Sampling
new_sampling = unique(round(logspace(0,log10(nr_time),10)));
cv_iters = 5;
data_test = data_test(:,:,new_sampling);
%% Training
r_squared = zeros(3,cv_iters);
iter_diets = 0;
diet_vector = [3 5 10];
for i_diets = diet_vector% [3 5 10 15 20]%[3 5 10]
    iter_diets = iter_diets + 1;
    nr_diets_chosen = i_diets;
    data_metabolomics_training = data_metabolomics(1:50,:,:,1:17);
    nr_time = size(data_metabolomics,3);
    data_metabolomics_mix = zeros(nr_diets_chosen*17,nr_metabolomics, nr_time);
    input_meal_all_mix = zeros( size(input_meal_all,1), nr_diets_chosen*17);
    for inds = 1:17
        indx = ((inds-1)*(nr_diets_chosen-1))+(inds:inds+(nr_diets_chosen-1));
        data_metabolomics_mix(indx,:,:) = data_metabolomics_training(1:nr_diets_chosen,:,:,inds);
        input_meal_all_mix(:,indx) = input_meal_all(:,1:nr_diets_chosen);
    end
    data_metabolomics_training = data_metabolomics_mix;
    input_meal_all = input_meal_all_mix;
    for i = 1:size(data_metabolomics_training,1)
        data_metabolomics_training(i,:,:) = normalize( squeeze(data_metabolomics_training(i,:,:))' )';
    end

    %% Sampling
    data_metabolomics_training = data_metabolomics_training(:,:,new_sampling);
    nr_time = size(data_metabolomics_training,3);
    %%
    %plot(1:7, normalize(squeeze(data_metabolomics(1,:,:))'))
    %% Randomize
    

    for i_cv = 1:cv_iters
        random_indices = randperm(size(data_metabolomics_training,1));
        data_metabolomics_random = data_metabolomics_training(random_indices,:,:);
        input_meal_all_random = input_meal_all(:,random_indices);
        %% Set the sets of data


        nr_diets_tot = size(data_metabolomics_random,1);
        ratio_train = 0.8;
        diets_training = 1:round(nr_diets_tot*ratio_train);                                                     % Change number of diets in training to see how the predictions change
        diets_validation = (round(nr_diets_tot*ratio_train)+1):nr_diets_tot;
        data_training = squeeze(data_metabolomics_random(diets_training,:,:));     % 4 diets in the training set for individual 1
        data_validation = squeeze(data_metabolomics_random(diets_validation,:,:));   % 1 diet in the validation set set for individual 1
        input_meal_training = input_meal_all(:,diets_training);                               % meal information for training set
        input_meal_validation = input_meal_all(:,diets_validation);                              % meal information for validation set
        %% pDMDc training
        stop_crit = 0.9;
        [A, B, prediction_val_best, prediction_train_best] = pDMDcPredictionValidation(data_training, ...
                                                                                       data_validation, ...
                                                                                       input_meal_training, ...
                                                                                       input_meal_validation, ...
                                                                                       stop_crit); % Train and predict
        %% pDMDc testing
        prediction_test = pDMDcPredictionTest(data_test, input_test, A, B); % Using estimated A and B from training to predict the test set from the input
        %% Plot training reconstruction
        % disp(['rank: ', num2str(rank(A))])
        % random_metabolites = randperm(nr_metabolomics,10);
        % colors = lines();
        % iter = 0;
        % figure()
        % for m = random_metabolites
        %     iter = iter + 1;
        %     subplot(2,5,iter)
        %     for j = 1:length(diets_training)
        %         plot(1:nr_time, squeeze(prediction_train_best(j,m,:)), Color= colors(j,:), LineWidth=1); hold on
        %         plot(1:nr_time, squeeze(data_training(j,m,:)),'.', Color= colors(j,:));
        %     end
        %     title(["Metabolite #",num2str(m)])
        % end
        % p = [];
        % legend_names = {};
        % for j = 1:length(diets_training)
        %     p(j) = plot(nan,nan,  Color= colors(j,:));
        %     legend_names{j} = strcat("Diet ", num2str(j) );
        % end
        % p(end+1) = plot(nan,nan,  'k-');
        % p(end+1) = plot(nan,nan,  'k.-');
        % legend_names{end+1} = "Prediction";
        % legend_names{end+1} = "Data";
        % legend(p, legend_names)
        %% Plot test prediction
        % random_metabolites = randperm(nr_metabolomics,10);
        % iter = 0;
        % colors = lines();
        % figure()
        % for m = random_metabolites
        %     iter = iter + 1;
        %     subplot(2,5,iter)
        %     for j = 1:size(data_test,1)
        %         plot(1:nr_time, squeeze(prediction_test(j,m,:)), Color= colors(j,:), LineWidth=1); hold on
        %         plot(1:nr_time, squeeze(data_test(j,m,:)),'o-', Color= colors(j,:));
        %         pause
        %     end
        %     title(["Metabolite #",num2str(m)])
        % end
        % p = [];
        % legend_names = {};
        % for j = 1:length(diets_test)
        %     p(j) = plot(nan,nan,  Color= colors(j,:));
        %     legend_names{j} = strcat("Diet ", num2str(j) );
        % end
        % p(end+1) = plot(nan,nan,  'k-');
        % p(end+1) = plot(nan,nan,  'k.-');
        % legend_names{end+1} = "Prediction";
        % legend_names{end+1} = "Data";
        % legend(p, legend_names)

        %% Prediction and data in scatter plot
        %figure()
        p_normalized = zeros(size(prediction_test));
        d_normalized = zeros(size(data_test));
        for i = 1:size(data_test,1)
            %     tmp = [squeeze(data_test(i,:,:))';
            %       squeeze(prediction_test(i,:,:))']
            d_normalized(i,:,:) = normalize(squeeze(data_test(i,:,:))')';
            p_normalized(i,:,:) = normalize(squeeze(prediction_test(i,:,:))')';
        end

        % plot(p_normalized(:), d_normalized(:),'.'); hold on
        % plot(p_normalized(:), p_normalized(:),'k')
        % xlabel('Prediction')
        % ylabel('Data')
        r_squared(iter_diets,i_cv) = 1 - sum([prediction_test(:)-data_test(:)].^2)./sum([mean(data_test(:)) - data_test(:)].^2);  %mdl.Rsquared.Adjusted;


    end
end
%%
figure
boxplot(r_squared')
xticklabels( diet_vector)
xlabel('Diets in training set','Interpreter','latex')
ylabel('$R^2$ on test set','Interpreter','latex')

%------------- END OF CODE --------------