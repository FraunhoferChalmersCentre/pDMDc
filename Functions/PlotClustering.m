function PlotClustering(cluster_ground_truth, ...
                        X_states_tensor, model, ...
                        data, ...
                        list_amino_acids, ...
                        metabolite_names, ...
                        state_of_interest, ...
                        diet_of_interest, ...
                        denom)                                        
%% Preallocation
nr_time = size(X_states_tensor,3);

%% Clustering state trajectories
metric = 2;                                                                % 1: Cosine similarity, 2: Sums of squares
nr_clusters = 2;                                                           % Number of clusters
n_time = size(X_states_tensor,3);
clustering_states = ClusteringFunction(-X_states_tensor, metric,state_of_interest, diet_of_interest, nr_clusters);

res_dmd = max( 1- sum(clustering_states == grp2idx(categorical(cluster_ground_truth)))/length(clustering_states), ...
    sum(clustering_states == grp2idx(categorical(cluster_ground_truth)))/length(clustering_states) );
disp('-------------------------------------------------------------')
disp(['Clustering result of pDMDc: ', num2str(res_dmd*100),'% match'])
disp('-------------------------------------------------------------')
%% Plot clustered states
figure()
%subplot(121)
colors = lines(2);
line_markers = {'--','-'};
T = clustering_states;
for i = unique(T)'
    plot((1:nr_time)/denom, -1*real( squeeze( X_states_tensor(i==T,state_of_interest,:,diet_of_interest)) )', line_markers{i}, 'Color', colors(i,:), 'LineWidth', 2); hold on
end
xlabel('\textbf{Time (h)}', 'Interpreter','latex')
ylabel('\textbf{State value}', 'Interpreter','latex')
title(['Clustering in pDMDc state: ', num2str(state_of_interest)])
%% 
nr_comps = size(model{4},2);
save_list = {};
res_cp_best = 0;
for i = 1:nr_comps
    tmp_choose = nchoosek(1:nr_comps,i);    
    save_list = [save_list; mat2cell(tmp_choose, ones(size(tmp_choose, 1), 1), size(tmp_choose, 2))];
end

for i = 1:length(save_list)
    gmfit = fitgmdist(model{4}(:,save_list{i}),2);
    clusterX = cluster(gmfit,model{4}(:,save_list{i}));
    res_cp = max( 1- sum(clusterX == grp2idx(categorical(cluster_ground_truth)))/length(clusterX), ...
        sum(clusterX == grp2idx(categorical(cluster_ground_truth)))/length(clusterX) );
    if res_cp>res_cp_best
        res_cp_best = res_cp;
    end
end


disp(['Clustering result of CP: ', num2str(res_cp_best*100),'% match'])
%%
figure()
cluster_ground_truth_s = cellstr( cluster_ground_truth );
nr_comps = size(model{4},2);
label = {};
for i = 1:nr_comps
    label{i} = ['CP Component ', num2str(i)];
end

pairplot(model{4}, label, cluster_ground_truth_s, model{4})
sgtitle('Clustering in all PARAFAC components')
%%

iter = 0;
figure()
for i = list_amino_acids
    iter = iter + 1;
    subplot(3,2,iter)
    for j = unique(T)'
        plot((1:nr_time)/denom,squeeze(data(2,i,:,T == j)), line_markers{j}, 'Color', colors(j,:));hold on
        ylim([ min(squeeze(data(2,i,:,:)),[],'all') - std(squeeze(data(2,i,:,:)),[],'all'), max(squeeze(data(2,i,:,:)),[],'all') + std(squeeze(data(2,i,:,:)),[],'all')])
    end
    
    title(metabolite_names{iter},'Interpreter','latex')
    if iter>4
        xlabel('Time (h)','Interpreter','latex')
    end
    if any( iter == [1 3 5] )
        ylabel('Amplitude','Interpreter','latex')
    end
end
sgtitle('Clusters color coded in simulated response data')
end