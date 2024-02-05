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
function [A_save, B_save, prediction_val_best, prediction_train_best, rank_best, rank_max] = pDMDcPrediction(data_training, data_validation, input_meal_training, input_meal_validation)      
%% Set up the DMD data format
            time = size(data_training,3); %nr time points in data
            nr_diets_training = size(data_training,1);
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
                x0_val = squeeze(data_validation(:,1));
                states_val = x0_val;
                for j= 2:time
                    if j == 2
                        states_val(:,j) = A_train*states_val(:,j-1) + B_train*input_meal_validation;
                    else
                        states_val(:,j) = A_train*states_val(:,j-1);
                    end
                end
                prediction_val = states_val;
                %% Error on validation set
                error(iter) = sqrt( sum( ( prediction_val(:) - data_validation(:) ).^2, 'all')/numel(prediction_val) );
                [~, idx] = min(error(1:iter));
                    if  idx == iter % save if the best error
                        A_save = A_train;
                        B_save = B_train;
                        prediction_val_best = prediction_val;
                        prediction_train_best = prediction_train;
                        rank_best = rank(A_save);
                    end
            end

end
%------------- END OF CODE --------------