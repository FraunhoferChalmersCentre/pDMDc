%% pDMDc algorithm for prediction
%
% Example: 
%    run pDMDc_prediction_example.m
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: nan
%
% Author: Viktor Skantze
% Work address: 
% email: viktor.skantze@fcc.chalmers.se
%------------- BEGIN CODE --------------
function [A_save, B_save, prediction_val_best, prediction_train_best, rank_best, rank_max] = pDMDcPredictionValidation(data_training, data_validation, input_meal_training, input_meal_validation, stop_crit)      
%% Set up the DMD data format
            time = size(data_training,3); %nr time points in data
            nr_diets_training = size(data_training,1);
            nr_diets_validation = size(data_validation,1);
            X_1_train = [];
            X_2_train = [];
            nr_macros = size(input_meal_training,1);
            stacked_input = [];
            for i = 1:nr_diets_training
                tmp_training = squeeze( data_training(i,:,:) );
                X_1_train = [X_1_train, tmp_training(:,1:time-1)];
                X_2_train = [X_2_train, tmp_training(:,2:time)];
                stacked_input = [stacked_input [input_meal_training(:,i) zeros(nr_macros,time-2)]];
            end

            input_space = [X_1_train;stacked_input];


            iter = 0;
            rank_max = min(rank(input_space), rank(X_2_train));
            error = zeros(rank_max,1);
            reconstruction_error = zeros(rank_max,1);
            iter_2 = 0;
            A_saving = {};
            B_saving = {};
            prediction_train = {};
            prediction_val_best = {};
            idx_save = 0;
            %% DMDc with all ranks
            for r = 1:rank_max
                iter = iter + 1;
                [U_input, S_input, V_input] = svd(input_space); 
                U_input = U_input(:, 1:r); % truncate to rank-r
                U_input_X = U_input(1:size(X_1_train,1),:);
                U_input_u = U_input((size(X_1_train,1)+1):end,:);
                S_input = S_input(1:r, 1:r);
                V_input = V_input(:, 1:r);

                [U_2, S_2, V_2] = svd(X_2_train); % svd on the output space
                U_2 = U_2(:, 1:r); % truncate to rank-r
                S_2 = S_2(1:r, 1:r);
                V_2 = V_2(:, 1:r);

                A_train = X_2_train*V_input*inv(S_input)*U_input_X'; % dynamics
                B_train = X_2_train*V_input*inv(S_input)*U_input_u'; % input mapping

                %% Training set
                prediction_train = zeros(size(data_training));
                for jj = 1:nr_diets_training
                    x0_train = squeeze(data_training(jj,:,1))';
                    states_train = x0_train;
                    for j= 2:time
                        if j == 2
                            states_train(:,j) = A_train*states_train(:,j-1) + B_train*input_meal_training(:,jj);
                        else
                            states_train(:,j) = A_train*states_train(:,j-1);
                        end
                    end
                    prediction_train(jj,:,:) = states_train;
                end
                %% Validation set
                prediction_val = zeros(size(data_validation));
                for jj = 1:nr_diets_validation
                x0_val = squeeze(data_validation(jj,:,1))';
                states_val = x0_val;
                for j= 2:time
                    if j == 2
                        states_val(:,j) = A_train*states_val(:,j-1) + B_train*input_meal_validation(:,jj);
                    else
                        states_val(:,j) = A_train*states_val(:,j-1);
                    end
                end
                    prediction_val(jj,:,:) = states_val;
                end
                %% Error on validation set
                error(iter) = sqrt( sum( ( prediction_val(:) - data_validation(:) ).^2, 'all')/numel(prediction_val) );
                reconstruction_error(iter) = sqrt( sum( ( prediction_train(:) - data_training(:) ).^2, 'all')/numel(prediction_train) );
                if (reconstruction_error(iter)<error(iter)*stop_crit) %&& (iter>10) && (iter<80)
                    [~, idx] = min(error(1:iter));
                    %idx = r;
                    %idx = iter;
                    iter_2 = iter_2 +1;
                else 
                    idx = 0;
                end
                %[~, idx] = min(error(1:iter));
                %[~, idx] = min(error(1:iter));
                A_saving{r} = A_train;
                B_saving{r} = B_train;
                prediction_val_save{r} = prediction_val;
                prediction_train_save{r} = prediction_train;
                    if  (iter_2 == 1)% save if the best error
                        A_save = A_saving{idx};
                        B_save = B_saving{idx};
                        prediction_val_best = prediction_val;
                        prediction_train_best = prediction_train;
                        rank_best = rank(A_save);
                        idx_save = idx;
                        iter_2 = iter_2 + 1;
                    end
                 
            end
       

            if idx_save==0
                %disp('min')
                [~, idx_save] = min(error(1:iter));
                %idx_save = 80;
                A_save = A_saving{idx_save};
                B_save = B_saving{idx_save};
                prediction_val_best = prediction_val;
                prediction_train_best = prediction_train;
                rank_best = rank(idx_save);
            end
%             plot(1:iter, error,'g'); hold on
%             plot(1:iter, reconstruction_error,'k'); hold on
%             plot(idx_save, error(idx_save),'ro')
%             legend('Validation error', ...
%                     'Reconstruction error', ...
%                     'Chosen nr of components for test set')
%             xlabel('# pDMDc components')
%             ylabel('RMSE')
end
%------------- END OF CODE --------------