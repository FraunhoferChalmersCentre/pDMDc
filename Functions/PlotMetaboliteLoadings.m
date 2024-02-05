function PlotMetaboliteLoadings(metabolite_list, ...
                                U, ...
                                data_metabolomics_dmd, ...
                                parafac_loadings, ...
                                classification_description_dynamics, ...
                                classification_description_metabolites, ...
                                classification_description_cp)

%% Classification dynamics in the sam
%figure()
[classification_list, classification_description] = ClassificationMetabolites(metabolite_list);
[classification_list_sorted, indexing_sort] = sort(classification_list);
nr_components = 4;
colors = distinguishable_colors(10);
[~, classification_dynamics] = ClassificationDynamics(data_metabolomics_dmd,4);
[classification_list_sorted_dynamics, indexing_sort_dynamics] = sort(classification_dynamics);
colors_dynamic = colors(7:end,:);
U_sorted_dynamics = U(indexing_sort_dynamics,:);
%U_sorted_dynamics = U_pooled_estimate(indexing_sort_dynamics,:);
[classification_dynamics_other_sorting] = classification_dynamics(indexing_sort);
U_sorted_dynamics_other_sorting = U(indexing_sort,:);

barwidth = 0.8;
list_indices = [1 2 4 3];
for i = list_indices
%     subplot(3,4,i)
%     for j = unique(classification_dynamics_other_sorting)'
%         %         if j== 2
%         %             barwidth = 0.25;
%         %         else
%         %             barwidth = 1;
%         %         end
%         %plot(find(classification_list_sorted_dynamics == j),U_sorted_dynamics(classification_list_sorted_dynamics == j,i),'o','MarkerFaceColor',colors_dynamic(j,:),'MarkerEdgeColor','k');hold on
%         if i == 2
%         bar(find(classification_list_sorted_dynamics == j), ...
%             U_sorted_dynamics(classification_list_sorted_dynamics == j,i), ...
%             'FaceColor',colors_dynamic(j,:), ...
%             'EdgeColor',colors_dynamic(j,:), ...
%             'BarWidth',barwidth);hold on
%         else
%                 bar(find(classification_list_sorted_dynamics == j), ...
%             -U_sorted_dynamics(classification_list_sorted_dynamics == j,i), ...
%             'FaceColor',colors_dynamic(j,:), ...
%             'EdgeColor',colors_dynamic(j,:), ...
%             'BarWidth',barwidth);hold on
%         end
%     end
%     %title(strcat('Column vector ', {' '},num2str(i),' in $ \hat{\bf{U}}_{tot} $'),'fontweight','bold', Interpreter="latex")
%     %title(classification_description_dynamics_title{i},Interpreter="latex")
%     xlabel('Dynamic index','interpreter','latex')
%     ylabel( strcat( [' $\bf{U}_{*', num2str(i),'}','^{(tot)}$']), 'interpreter',"latex")
% end
% legend(classification_description_dynamics, 'orientation','horizontal', Interpreter="latex")
%%

%colors = {'b','r','k','g','m','y','c'};

colors_metabolic = lines(6);%colors(1:6,:);

%U_sparse = spca(X', [], nr_components, inf, 40);
U_sorted = U(indexing_sort,:);
sum( U_sorted==0, 'all')
%title_list = {'Amino acid '}
for i = 1:nr_components
    subplot(2,4,i)
    for j = unique(classification_list_sorted)'
        %plot(find(classification_list_sorted == j),U_sorted(classification_list_sorted == j,i),'o','MarkerFaceColor',colors_metabolic(j,:),'MarkerEdgeColor','k');hold on
        if i == 2

            bar(find(classification_list_sorted == j), ...
                U_sorted(classification_list_sorted == j,i), ...
                'FaceColor',colors_metabolic(j,:), ...
                'EdgeColor',colors_metabolic(j,:));hold on
            xticks([]);
        else
            bar(find(classification_list_sorted == j), ...
                -U_sorted(classification_list_sorted == j,i), ...
                'FaceColor',colors_metabolic(j,:), ...
                'EdgeColor',colors_metabolic(j,:));hold on
            xticks([]);
        end

    end

    title(classification_description_metabolites{i},Interpreter="latex")
    xlabel('Metabolite index','interpreter','latex')
    ylabel( strcat( [' $\bf{U}_{*', num2str(i),'}','^{(tot)}$']), 'interpreter',"latex")

    if i == 4
        colororder({'k','k'})
        yyaxis right
        ylabel('pDMDc','Color',[0.15 0.15 0.15], ...
            'Interpreter','latex', 'FontSize',12)
        yticks([])
    end
end
legend(classification_description,'Orientation','horizontal')


%%
%colors = {'b','r','k','g','m','y','c'};


%U_sparse = spca(X', [], nr_components, inf, 40);
parafac_loadings_sorted = 1*parafac_loadings(indexing_sort,:);
%title_list = {'Amino acid '}
iter = 0;
for i = 1:nr_components-1
    iter = iter + 1;
    if i == 3
        i = 4;
        subplot(2,4,i+4)
        i = 3;
    else
        subplot(2,4,i+4)
    end
    for j = unique(classification_list_sorted)'

        bar(find(classification_list_sorted == j), ...
            parafac_loadings_sorted(classification_list_sorted == j,i), ...
            'FaceColor',colors_metabolic(j,:), ...
            'EdgeColor',colors_metabolic(j,:));hold on
        xticks([]);
    end

    title(classification_description_cp{i},Interpreter="latex")
    xlabel('Metabolite index','interpreter','latex')
    ylabel( strcat( [' $ \bf{r_', num2str(i),'}$']), 'interpreter',"latex")
    if iter == 3
        colororder({'k','k'})
        yyaxis right
        ylabel('CP','Color',[0.15 0.15 0.15], ...
            'Interpreter','latex', 'FontSize',12)
        yticks([])
    end
end
%legend(classification_description,'Orientation','horizontal')
%% Correlation 
disp(['Correlation between U_1 and r_1: ', ...
    num2str(corr( -U_sorted(:,1), parafac_loadings_sorted(:,1)))])
disp(['Correlation between U_2 and r_2: ', ...
    num2str(corr( U_sorted(:,2), parafac_loadings_sorted(:,2)))])
disp(['Correlation between U_4 and r_3: ', ...
    num2str(corr( -U_sorted(:,4), parafac_loadings_sorted(:,3)))])
disp(['Correlation between U_3 and r_3: ', ...
    num2str(corr( -U_sorted(:,3), parafac_loadings_sorted(:,3)))])
end