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
function [prediction_test] = pDMDcPredictionTest(data_test, input_meal_test, A_train, B_train)      
%% Set up the DMD data format
            time = size(data_test,3); %nr time points in data
            nr_diets_test = size(data_test,1);

            %% Training set
            prediction_test = zeros(size(data_test));
            for jj = 1:nr_diets_test
                x0_test = squeeze(data_test(jj,:,1))';
                states_train = x0_test;
                for j= 2:time
                    if j == 2
                        states_train(:,j) = A_train*states_train(:,j-1) + B_train*input_meal_test(:,jj);
                    else
                        states_train(:,j) = A_train*states_train(:,j-1);
                    end
                end
                prediction_test(jj,:,:) = states_train;
            end


end
%------------- END OF CODE --------------