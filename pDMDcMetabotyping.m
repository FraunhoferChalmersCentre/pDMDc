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
function [A_train, B_train, prediction_train, X_states_tensor] = pDMDcMetabotyping(data_metabolomics, input_meal, r)
%% Normalization
data_metabolomics = permute(data_metabolomics, [4 3 2 1]);
data = zeros(size(data_metabolomics));
prediction_train = zeros(size(data_metabolomics));
nr_diets = size(data,4);
nr_individuals = size(data,1);
nr_time_points = size(data,2);
for i = 1:nr_diets
    for j = 1:79
        tmp = squeeze(data_metabolomics(:,:,j,i));
        data(:,:,j,i) = ( tmp - mean(tmp,'all') ) / std(tmp,0,'all');
    end
end
%% Set up the DMD data format
X1 = [];
X2 = [];
X = [];
for i = 1:nr_individuals
    for d = 1:nr_diets
        tmp_dataset = squeeze(data(i,:,:,d));
        X = [X, tmp_dataset'];
        X1 = [X1, tmp_dataset(1:end-1,:)'];
        X2 = [X2, tmp_dataset(2:end,:)'];
    end
end

stacked_input = [];
for d = 1:nr_diets
    stacked_input = [stacked_input, [input_meal(:,d) zeros(3,nr_time_points-2)]];
end
%% pDMDc

[U,~,~] = svd(X);

U = U(:,1:r);

X_states_tensor = zeros(17,r,nr_time_points,nr_diets);

length_individual_data = (nr_time_points-1)*nr_diets;
for i = 1:17
    X1_ind =  X1(:,(1:length_individual_data)+(i-1)*length_individual_data);
    X2_ind =  X2(:,(1:length_individual_data)+(i-1)*length_individual_data);
    X_ind_1_projected = U' * X1_ind;
    X_ind_2_projected = U' * X2_ind;

    input_space = [X_ind_1_projected;stacked_input];
    G = X_ind_2_projected*pinv(input_space);
    A_train = G(:,1:r);
    B_train = G(:,r+1:end);
%% State space model

       for jj = 1:nr_diets
            x0_train = U'*X1_ind(:,1+(jj-1)*(nr_time_points-1));
            states_train = x0_train;
            for j= 2:nr_time_points
                if j == 2
                    states_train(:,j) = A_train*states_train(:,j-1) + B_train*input_meal(:,jj);
                else
                    states_train(:,j) = A_train*states_train(:,j-1);
                end
            end
            prediction_train(i,:,:,jj) = (U*states_train)';

            X_states_tensor(i,:,:,jj) = states_train;
   
        end

if i == 1
    for ddd = 1:3
    plot(1:8, squeeze(data(i,:,:,ddd)),'k', 'linewidth',2); hold on
    plot(1:8, squeeze(prediction_train(i,:,:,ddd)))
    xlabel('Time points')
    ylabel("Metabolite concentration normalized")
    title(['Data (black) & Prediction (colors) of individual 1, diet #', num2str(ddd)])
    pause()
    clf
    end
end
end



end
%------------- END OF CODE --------------