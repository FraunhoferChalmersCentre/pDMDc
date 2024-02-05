function [classification_list,classification_description] = ClassificationMetabolites(metabolite_list)
index_blood = find(contains(metabolite_list,'_B'));
index_liver = find(contains(metabolite_list,'_L'));
index_muscle = find(contains(metabolite_list,'_M'));
index_adipose_tissue = find(contains(metabolite_list,'_A'));
index_gut = find(contains(metabolite_list,'_G'));


classification_list = ones(length(metabolite_list),1);
classification_list(index_blood) = 2;
classification_list(index_liver) = 3;
classification_list(index_muscle) = 4;
classification_list(index_adipose_tissue) = 5;
classification_list(index_gut) = 6;


classification_description = {'Blood','Liver','Muslce','Adipose tissue','Gut'};

end


