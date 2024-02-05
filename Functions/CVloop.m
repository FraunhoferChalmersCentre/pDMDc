function output = CVloop(data_metabolomics, input_meal_all, ratio_train, ratio_test)



%% Randomize
cv_iters = 30;
r_squared = zeros(2,cv_iters);
for i = 1:size(data_metabolomics,1)
    data_metabolomics(i,:,:) = fillmissing(normalize(squeeze(data_metabolomics(i,:,:))')', ...
        'nearest');
end

for i_cv = 1:cv_iters
    random_indices = randperm(size(data_metabolomics,1));
    data_metabolomics_random = data_metabolomics(random_indices,:,:);
    input_meal_all_random = input_meal_all(:,random_indices);
    %% Set the sets of data
    %data_metabolomics = nprocess(data_metabolomics,[0 0 0],[0 1 0], [],[],1);
    nr_metabolomics = size(data_metabolomics,2);
    nr_time = size(data_metabolomics,3);

    nr_diets_tot = size(data_metabolomics_random,1);

    diets_training = 1:round(nr_diets_tot*ratio_train);                                                     % Change number of diets in training to see how the predictions change
    diets_validation = (diets_training(end)+1):round(nr_diets_tot*ratio_test);
    diets_test = (diets_validation(end)+1):nr_diets_tot;
    data_training = squeeze(data_metabolomics_random(diets_training,:,:));     % 4 diets in the training set for individual 1
    data_validation = squeeze(data_metabolomics_random(diets_validation,:,:));   % 1 diet in the validation set set for individual 1
    data_test = squeeze(data_metabolomics_random(diets_test,:,:));
    input_meal_training = input_meal_all_random(:,diets_training);                               % meal information for training set
    input_meal_validation = input_meal_all_random(:,diets_validation);                              % meal information for validation set
    input_meal_test = input_meal_all_random(:,diets_test);                              % meal information for validation set
    %% pDMDc training
    [A, B, ~, prediction_train_best] = pDMDcPredictionValidation(data_training, data_validation, input_meal_training, input_meal_validation); % Train and predict
    %% pDMDc testing
    prediction_test = pDMDcPredictionTest(data_test, input_meal_test, A, B);
    %% Plot training reconstruction
%     disp(['rank: ', num2str(rank(A))])
%     random_metabolites = randperm(nr_metabolomics,10);
%     colors = lines();
%     iter = 0;
%     figure()
%     for m = random_metabolites
%         iter = iter + 1;
%         subplot(2,5,iter)
%         for j = 1:length(diets_training)
%             plot(1:nr_time, squeeze(prediction_train_best(j,m,:)), Color= colors(j,:), LineWidth=1); hold on
%             plot(1:nr_time, squeeze(data_training(j,m,:)),'.', Color= colors(j,:));
%         end
%         title(["Metabolite #",num2str(m)])
%     end
%     p = [];
%     legend_names = {};
%     for j = 1:length(diets_training)
%         p(j) = plot(nan,nan,  Color= colors(j,:));
%         legend_names{j} = strcat("Diet ", num2str(j) );
%     end
%     p(end+1) = plot(nan,nan,  'k-');
%     p(end+1) = plot(nan,nan,  'k.-');
%     legend_names{end+1} = "Prediction";
%     legend_names{end+1} = "Data";
%     legend(p, legend_names)
    %% Plot test prediction
%     random_metabolites = randperm(nr_metabolomics,10);
%     colors = lines();
%     iter = 0;
%     figure()
%     for m = random_metabolites
%         iter = iter + 1;
%         subplot(2,5,iter)
%         for j = 1%:size(data_test,1)
%             plot(1:nr_time, squeeze(prediction_test(j,m,:)), Color= colors(j,:), LineWidth=1); hold on
%             plot(1:nr_time, squeeze(data_test(j,m,:)),'.', Color= colors(j,:));
%         end
%         title(["Metabolite #",num2str(m)])
%     end
%     p = [];
%     legend_names = {};
%     for j = 1:length(diets_test)
%         p(j) = plot(nan,nan,  Color= colors(j,:));
%         legend_names{j} = strcat("Diet ", num2str(j) );
%     end
%     p(end+1) = plot(nan,nan,  'k-');
%     p(end+1) = plot(nan,nan,  'k.-');
%     legend_names{end+1} = "Prediction";
%     legend_names{end+1} = "Data";
%     legend(p, legend_names)
%     pause()
%% Prediction and reconstruction
%     random_metabolites = randperm(nr_metabolomics,10);
%     colors = lines();
%     iter = 0;
%     figure()
%     for m = random_metabolites
%         iter = iter + 1;
%         subplot(2,5,iter)
%    
%         for j = 1:size(data_test,1)
%             plot(1:nr_time, squeeze(prediction_train_best(j,m,:)), Color= [0.8 0.8 .8], LineWidth=1); hold on
%             plot(1:nr_time, squeeze(data_training(j,m,:)),'.', Color= [0.8 0.8 .8]);
%         end
%         plot(1:nr_time, squeeze(prediction_test(j,m,:)), Color= colors(2,:), LineWidth=1); hold on
%         plot(1:nr_time, squeeze(data_test(j,m,:)),'.', Color= colors(2,:));
%         title(["Metabolite #",num2str(m)])
%     end
% 
%     p(1) = plot(nan,nan,  '-','color',[0.8 0.8 .8]);
%     p(2) = plot(nan,nan,  '-','color',colors(2,:));
%     p(end+1) = plot(nan,nan,  'k-');
%     p(end+1) = plot(nan,nan,  'k.-');
%     legend_names{1} = "Training";
%     legend_names{2} = "Prediction";
%     legend_names{end+1} = "Prediction";
%     legend_names{end+1} = "Data";
%     legend(p, legend_names)
%     pause()
    %% Prediction and data in scatter plot
    %figure()
%     p_normalized = zeros(size(prediction_test));
%     d_normalized = zeros(size(data_test));
%     for i = 1:size(data_test,1)
%         %     tmp = [squeeze(data_test(i,:,:))';
%         %       squeeze(prediction_test(i,:,:))']
%         d_normalized(i,:,:) = normalize(squeeze(data_test(i,:,:))')';
%         p_normalized(i,:,:) = normalize(squeeze(prediction_test(i,:,:))')';
%     end
    
%     figure()
%     plot(prediction_test(:), data_test(:),'.'); hold on
%     plot(data_test(:), data_test(:),'k')
%     xlabel('Prediction')
%     ylabel('Data')
    %pause()
    mdl = fitlm(prediction_test(:),data_test(:));
    r_squared(1,i_cv) = mdl.Rsquared.Adjusted;

    mdl = fitlm(prediction_train_best(:),data_training(:));
    r_squared(2,i_cv) = mdl.Rsquared.Adjusted;
end
%%
% figure()
% plot(1:cv_iters,r_squared)
output = mean(r_squared,2);
end