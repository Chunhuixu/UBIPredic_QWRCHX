% for i=1:1600
%     P_Ubidata{i} = rmfield(P_Ubidata{i},'KW');
%     P_Ubidata{i} = rmfield(P_Ubidata{i},'PSSM');
%     P_Ubidata{i} = rmfield(P_Ubidata{i},'SQ');
% end
% for j=1:3
%     for i=1:1600
%         N_Ubidata{j}{i} = rmfield(N_Ubidata{j}{i},'KW');
%         N_Ubidata{j}{i} = rmfield(N_Ubidata{j}{i},'PSSM');
%         N_Ubidata{j}{i} = rmfield(N_Ubidata{j}{i},'SQ');
%     end
% end
%save P_Ubi_Predict 
Feature_All{1}(:,176:200)=[];
PTM_Ubi_01_model=classRF_train([Feature_All_PSSM{1},Feature_All{1}],labelset,100);
Feature_All{2}(:,176:200)=[];
PTM_Ubi_02_model=classRF_train([Feature_All_PSSM{2},Feature_All{2}],labelset,100);
Feature_All{3}(:,176:200)=[];
PTM_Ubi_03_model=classRF_train([Feature_All_PSSM{3},Feature_All{3}],labelset,100);
save P_Ubi_Predict PTM_Ubi_01_model PTM_Ubi_02_model PTM_Ubi_03_model '-append'
