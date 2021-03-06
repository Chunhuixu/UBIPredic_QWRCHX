%Reading the FDA information one by one 
% ���� fileFA='P61981'
function [PROFea,ErroS]= Read_PRT_Infor_Outline_Ubi(fileFA)
PROFea=[];PROFea.PRTID=fileFA;ErroS=[];
%PROFea.S_F  get the desiring result from the web or not
%%%%��¼InterPro,Pfam,PRINTS,PROSITE,SMART,SUPFAM
DomainSet={'DR   InterPro;','DR   Pfam;','DR   PRINTS;','DR   PROSITE;','DR   SMART;','DR   SUPFAM'};
Domain_Sign={'INTERPRO','PFAM','PRINTS','PROSITE','SMART','SUPFAM'};
 temseq='';    seqflag=0;    GOth=1;   GOs=[]; 
%% Get the input information from www.uniprot.org/uniprot
fid0= fopen('tmp.txt', 'w');
[str,b] = urlread(['https://www.uniprot.org/uniprot/',fileFA,'.txt']);
if b==1
    PROFea.S_F=1;
    fprintf(fid0, '%s', str);
    fclose(fid0);
    fid1 = fopen('tmp.txt', 'r');
    count = 0;     DomainNum=num2cell(ones(1,6));Tem_line='';
    for idom=1:6
        eval([ 'Domain_S_',num2str(idom),'=[];'])
        eval(['PROFea.',Domain_Sign{idom},'=Domain_S_',num2str(idom),';'])
    end
     Note_end=0;Note_begin=0;Note_break=0;
    while ~feof(fid1)
        line = fgetl(fid1) ;       
        if (length(line)>2)&(length(strfind(line(1:3),'ID'))==1) 
            IDl=line(6:20);
            IDl(find(isspace(IDl))) = []; 
            PROFea.ID=IDl;      
        end
        %notes for the protein information
        Ids(1)=length(strfind(line,'FT   DISULFID'));
        Ids(2)=length(strfind(line,'FT   MOD_RES'));
        Ids(3)=length(strfind(line,'FT   LIPID'));
        Ids(4)=length(strfind(line,'FT   CARBOHYD'));
        Ids(5)=length(strfind(line,'FT   CROSSLNK'));        
       if sum(Ids)>0
          Tem_line=[Tem_line,line(6:end)];Note_end=0;Note_begin=1;Note_break=0;
       elseif (Note_begin==1)&(length(strfind(line,'FT       '))>0) & (Note_break==0)
          Tem_line=[Tem_line,line(6:end)];
          Note_end=0;Note_begin=1;
       else
           Note_end=1;Note_break=1;
       end
       if Note_end*Note_begin==1; 
            Ids_S1=strfind(Tem_line,'DISULFID');
            Ids_S2=strfind(Tem_line,'MOD_RES');
            Ids_S3=strfind(Tem_line,'LIPID');
            Ids_S4=strfind(Tem_line,'CARBOHYD');
            Ids_S5=strfind(Tem_line,'CROSSLNK'); 
            Setsite=unique([Ids_S1,Ids_S2,Ids_S3,Ids_S4,Ids_S5,length(Tem_line)+1]);
            for iids=1:length(Setsite)-1
                temsubline=Tem_line(Setsite(iids):Setsite(iids+1)-1);  
                ind5=strfind(temsubline, ';');
                ind6=strfind(temsubline, '.');
                ind7=length(temsubline);       
                if length(strfind(temsubline, ';'))>0
                    inds=ind5(1);
                elseif length(strfind(temsubline, '.'))>0
                    inds=ind6(1);
                else
                    inds=ind7;
                end       
                Labflag=temsubline(23:inds-1);          
                Labflag(find(isspace(Labflag))) = []; 
                ind8=strfind(Labflag, '(');
            end
       end
            for idom=1:6
                if length(strfind(line, DomainSet{idom}))>0
                    indidom = strfind(line, ';');
                    indend = strfind(line, '.');
                    if length(indidom)<2
                        indidom(2)=indend;
                    end
                    temstr=line(indidom(1)+1:indidom(2)-1);
                    temstr(find(isspace(temstr))) = []; 
                    eval([ 'Domain_S_',num2str(idom),'{',num2str(DomainNum{idom}),'}=temstr;'])
                    DomainNum{idom}=DomainNum{idom}+1;
                end
            end
        %%%%GO
        if length(strfind(line,'DR   GO; GO:'))>0
            indgo= strfind(line, 'GO:');
            GOs{GOth}=line(indgo(end)+3:indgo(end)+9); 
            GOth=GOth+1;
        end
        PROFea.GO=GOs;
             %Protein Sequence
        if strfind(line,'//')==1
            seqflag=0;
        end
        if seqflag==1
            temseq=[temseq,line];
        end   
        if strfind(line,'SQ')==1
            temseq='';
            seqflag=1;
        end
    end
    for idom=1:6
         eval(['PROFea.',Domain_Sign{idom},'=Domain_S_',num2str(idom),';'])
    end
    ablank=find(temseq==' ');
    temseq(ablank)=[];
    PROFea.SQ=temseq;
    fclose(fid1);
else
    PROFea.S_F=0;
    ErroS=fileFA;
    fclose(fid0);
end
