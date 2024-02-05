function PlotDynamicStates(nr_components, ...
    classification_dynamics, ...
    classification_description_dynamics, ...
    data_metabolomics_dmd, ...
    X_states_tensor, ...
    parafac_loadings, ...
    list_indices, ...
    order_cp,...
    denominator)
%% Plot states and dynamic classification
figure()
gray = [0.8, 0.8, 0.8] ;
nr_clusters = 4;
nr_time = size(data_metabolomics_dmd,3);
basis = zeros(nr_time,nr_clusters);

iter = 0;
data = permute(data_metabolomics_dmd, [1 3 2 4]);
data_matrix = squeeze( mean( squeeze( mean(data,1) ), 3) );
data_matrix = normalize(data_matrix);

T = classification_dynamics;
for i = list_indices
    iter = iter + 1;
    subplot(3,4,iter)
    plot((1:nr_time)/denominator,(data_matrix(:,T==i)),'linewidth',1,'color',gray);hold on
    basis(:,i) = median((data_matrix(:,T==i)),2,"omitnan");
    plot((1:nr_time)/denominator,basis(:,i),'r','LineWidth',2)
    xticks((1:nr_time)/denominator)
%     xlim([1,nr_time/denominator])
    ylim([min((data_matrix),[],'all'), max((data_matrix),[],'all')])
    title(classification_description_dynamics{i}, 'Interpreter','latex')
    xlabel('time (h)','Interpreter','latex')
    ylabel('Scaled concentration','Interpreter','latex')
    if iter == 4
        colororder({'k','k'})
        yyaxis right
        ylabel('Average data','Color',[0.15 0.15 0.15], ...
            'Interpreter','latex', 'FontSize',12)
        yticks([])
    end
end


%% DMD states
diet = 1;
state_names = {'$\tilde{x}_{1,1:T,1:I,d=1}$','$\tilde{x}_{2,1:T,1:I,d=1}$','$\tilde{x}_{3,1:T,1:I,d=1}$','$\tilde{x}_{4,1:T,1:I,d=1}$'};
median_latent_state = zeros(nr_time, nr_components);
for i = 1:nr_components
    tmp_data = squeeze(X_states_tensor(:,i,:,diet));
    subplot(3,4,i+4)
    if i == 2
        plot((1:nr_time)/denominator,1*squeeze( X_states_tensor(:,i,:,diet)),'color',gray,'LineWidth',1 );hold on
    elseif  max(tmp_data,[],'all')<-1e6
        tmp_data_2 = normalize(tmp_data);
        plot((1:nr_time)/denominator,-1*tmp_data_2,'color',[0.5, 0.5, 0.5],'LineWidth',1 ); hold on
    else
        plot((1:nr_time)/denominator,-1*squeeze( X_states_tensor(:,i,:,diet)),'color',gray,'LineWidth',1 ); hold on
    end
    ylabel(state_names{i}, 'interpreter','latex','FontSize', 12)
    if i == 2
        tmp_basis = 1*squeeze( X_states_tensor(:,i,:,diet));
    elseif  max(tmp_data,[],'all')<-1e6
        tmp_basis = -1*tmp_data_2;
    else
        tmp_basis = -1*squeeze( X_states_tensor(:,i,:,diet));
    end
    plot((1:nr_time)/denominator, median(tmp_basis,1),'r','LineWidth',2)
    xticks((1:nr_time)/denominator)
    %xlim([1,nr_time/denominator])
    median_latent_state(:,i) = median(tmp_basis,1);

    if i == 2
        ylim([min(1*tmp_data,[],'all'), max(1*tmp_data,[],'all')])
    elseif max(tmp_data,[],'all')<-1e6
        ylim([min(tmp_data_2,[],'all')*10, max(-1*tmp_data_2,[],'all')*10])
    else
        ylim([min(-1*tmp_data,[],'all'), max(-1*tmp_data,[],'all')])
    end
    xlabel('time (h)','Interpreter','latex')
        
    if i == 4
        colororder({'k','k'})
        yyaxis right
        ylabel('pDMDc','Color',[0.15 0.15 0.15], ...
            'Interpreter','latex', 'FontSize',12)
        yticks([])
    end
end

%% Cp states
state_names = {'$\bf{q_1}$', '$\bf{q_2}$', '$\bf{q_3}$'};
iter_cp = 0;
for i = order_cp

    iter_cp = iter_cp + 1;
    if iter_cp == 3
        iter_cp = 4;
    end
%     if iter_cp == 2
%         iter_cp = 3;
%     end
    subplot(3,4,iter_cp+8)
    plot((1:nr_time)/denominator, parafac_loadings(:,i),'r','LineWidth',2)
    ylabel(state_names{i}, 'interpreter','latex')
    xticks((1:nr_time)/denominator)
    %xlim([1,nr_time/denominator])
    tmp_data = parafac_loadings;
    ylim([min(tmp_data,[],'all'), max(tmp_data,[],'all')])
    xlabel('time (h)','Interpreter','latex')
    if iter_cp == 4
        colororder({'k','k'})
        yyaxis right
        ylabel('CP','Color',[0.15 0.15 0.15], 'Interpreter','latex', 'FontSize',12)
        yticks([])
    end
end
%% Correlation 
basis = basis(:,list_indices);
disp(['Correlation between  x_1 and basis 1: ', ...
    num2str(corr( median_latent_state(:,1), basis(:,1)))])
disp(['Correlation between x_2 and basis 2: ', ...
    num2str(corr( median_latent_state(:,2), basis(:,2)))])
disp(['Correlation between x_3 and basis 3: ', ...
    num2str(corr( median_latent_state(:,3), basis(:,3)))])
disp(['Correlation between x_4 and basis 4: ', ...
    num2str(corr( median_latent_state(:,4), basis(:,4)))])
disp('---------------------------------------------------------------')

disp(['Correlation between  q_1 and basis 1: ', ...
    num2str(corr( parafac_loadings(:,1), basis(:,1)))])
disp(['Correlation between q_2 and basis 2: ', ...
    num2str(corr( parafac_loadings(:,2), basis(:,2)))])
disp(['Correlation between q_3 and basis 3: ', ...
    num2str(corr( parafac_loadings(:,3), basis(:,4)))])
disp('---------------------------------------------------------------')

disp(['Correlation between  x_1 and q_1: ', ...
    num2str(corr( median_latent_state(:,1), parafac_loadings(:,1)))])
disp(['Correlation between x_2 and q_2: ', ...
    num2str(corr( median_latent_state(:,2), parafac_loadings(:,2)))])
disp(['Correlation between x_4 and q_3: ', ...
    num2str(corr( median_latent_state(:,4), parafac_loadings(:,3)))])
disp(['Correlation between x_3 and q_3: ', ...
    num2str(corr( median_latent_state(:,3), parafac_loadings(:,3)))])


end