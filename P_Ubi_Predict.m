%%Input your Protein sequence or ID of swissprot
%%%For example YourInput='P23396'
%%YourInput=YourInput= 'MAVQISKKRKFVADGIFKAELNEFLTRELAEDGYSGVE...
%%VRVTPTRTEIIILATRTQNVLGEKGRRIRELTAVVQKRFGFPEGSVELYAEKVATRGLCA...
%%IAQAESLRYKLLGGLAVRRACYGVLRFIMESGAKGCEVVVSGKLRGQRAKSMKFVDGLMI...
%%HSGDPVNYYVDTAVRHVLLRQGVLGIKVKIMLPWDPTGKIGPKKPLPDHVSIVEPKDEIL...
%%PTTPISEQKGGKPEPPAMPQPVPTA'

function Is_or_Not_UbiP = P_Ubi_Predict(YourInput)
%%%initializing the output
Is_or_Not_UbiP='Your input may not be ubiquitination protein';
%%%Get the Protein ID
if length(YourInput)==6
    ProteinID{1}=YourInput;
    db='G:/blast/blast2.2.28+/db/swissprot';
else
    db='G:/blast/blast2.2.28+/db/swissprot';
    [p,ProteinID]= blastpssm_seq_likeP(YourInput, db);
end
if length(ProteinID)<1
    disp('Wrong input');
else
    %%%Get the outline of input protein 
    [PROFea,ErroS]= Read_PRT_Infor_Outline_Ubi(ProteinID{1});
    ProteinInput.PRTID=PROFea.PRTID;
    ProteinInput.ID=PROFea.ID;
    ProteinInput.GO=PROFea.GO;
    ProteinInput.SQ=PROFea.SQ;
    ProteinInput.INTERPRO=PROFea.INTERPRO;
    ProteinInput.PFAM=PROFea.PFAM;
    ProteinInput.PRINTS=PROFea.PRINTS;  
    ProteinInput.PROSITE=PROFea.PROSITE; 
    ProteinInput.SMART=PROFea.SMART;
    ProteinInput.SUPFAM=PROFea.SUPFAM;
    PROFea=[];
    PROFea= Read_PRT_Infor_SubceL_Ubi(ProteinID{1});
    ProteinInput.CLS=PROFea.CLS;    
    temDatcls=ProteinInput.CLS;
    temsign=upper(ProteinInput.CLS);
    load P_Ubi_Predict
    %Refining the Subcellular localization Keywords
    for j=1:length(temDatcls)  
        for k=1:311
            if (strcmp(temsign{j},upper(NewLocList{k,1}))==1)&length(upper(NewLocList{k,2}))>0
                temsign{j}=upper(NewLocList{k,2});
            end
        end
    end
    ProteinInput.CLS_53=unique(temsign);
    %%% Presenting the protein with greyPSSM model
    temPSSM= blastpssm_seq(ProteinInput.SQ, db);
    ProteinInput.PSSM=temPSSM(:,1:20)';
    ProteinInput.greyPSSM=greyPsePssm_seq(ProteinInput.PSSM',2);
end
%%% Presenting the protein with functional domain annotation model
Domainset={'GO','PFAM','SMART','PROSITE','SUPFAM','INTERPRO','PRINTS','CLS_53'};
ProteinInput_feature=[];
for setith=1:3
    temfeat=[];
    for domith=1:8    
       temDisMatrix=[];
       Datatem=P_Ubidata;
       for i=1:1600
           Datatem{i+1600}=N_Ubidata{setith}{i};
       end
       Datatem=Datatem(resortc);
        for i=1:3200
            eval(['temDi=','Datatem{',num2str(i),'}.',Domainset{domith},';'])
            if length(temDi)==0
                temDisMatrix(i,:)=2;
            else     
                    eval(['temDj=','ProteinInput.',Domainset{domith},';'])
                    if length(temDj)==0 temDj=[];end
                    lengthU=length(union(temDi,temDj));
                    lengthI=length(intersect(temDi,temDj));
                    if lengthU>0
                       temDisMatrix(i,1)=1-lengthI/lengthU;
                    else
                        temDisMatrix(i,1)=1;
                    end                   
            end
        end
        isnull=find(temDisMatrix==2);isnull_not=find(temDisMatrix~=2);
        meanvalue=mean(temDisMatrix(isnull_not));
        temDisMatrix(isnull)=meanvalue;
        DisMatrix_feat=[];
        DisMatrix_feat(1,:)=CountKNNScore_Dis_Vector(temDisMatrix,labelset,[0.003:0.003:0.075]);
        temfeat=[temfeat,DisMatrix_feat];
    end
    ProteinInput_feature{setith}=temfeat;
end
temVotes=[];
[~,temVotes(1,:)]=classRF_predict([ProteinInput.greyPSSM,ProteinInput_feature{1}],PTM_Ubi_01_model);
[~,temVotes(2,:)]=classRF_predict([ProteinInput.greyPSSM,ProteinInput_feature{2}],PTM_Ubi_02_model);
[~,temVotes(3,:)]=classRF_predict([ProteinInput.greyPSSM,ProteinInput_feature{3}],PTM_Ubi_03_model);
votes=sum(temVotes);
if max(votes)==votes(2);
   Is_or_Not_UbiP='Your input may be ubiquitination protein!';
else
    Is_or_Not_UbiP='Your input may not be ubiquitination protein!';
end
