%Reading the PTM SUBCELLULAR LOCATION information one by one 
% АэИз fileFA='P61981'
function PROFea= Read_PRT_Infor_SubceL_Ubi(fileFA)
PROFea=[];PROFea.PRTID=fileFA;
fid0= fopen('tmp.txt', 'w');
[str,b] = urlread(['https://www.uniprot.org/uniprot/',fileFA,'.txt']);
if b==1
    PROFea.S_F=1;
    fprintf(fid0, '%s', str);
    fclose(fid0);
    fid1 = fopen('tmp.txt', 'r');
    %fid1 = fopen('tmptry1.txt', 'r');  
    while ~feof(fid1)       
        line = fgetl(fid1) ;

        if (length(line)>2)&(length(strfind(line(1:3),'ID'))==1) 
            IDl=line(6:20);
            IDl(find(isspace(IDl))) = []; CellLS=[];Seqset='';
            PROFea.ID=IDl;seqflag=0;
        end   
        %if (seqflag==1)&length(strfind(line,'CC   -!-'))>0
        if (seqflag==1)&(length(strfind(line,'CC   -!-'))+length(strfind(line,'CC   ----------------')))>0
            break
        end
        if length(strfind(line,'CC   -!- SUBCELLULAR LOCATION'))>0
           seqflag=1;
           Seqset=line;
        end
        if (seqflag==1)&length(line)>2&length(strfind(line(1:3),'CC'))>0&length(strfind(line,'CC   -!-'))==0
            Seqset=[Seqset,line(10:end)];
        end
    end
    fclose(fid1);
    sign0=strfind(Seqset,':'); 
    sign1=strfind(Seqset,'Note=');
    sign2=strfind(Seqset,'.');
    sign3=strfind(Seqset,';');
    sign4=strfind(Seqset,',');
    sign5_0=strfind(Seqset,'{');
    sign5_1=strfind(Seqset,'}');
    for isign5=min(length(sign5_0),length(sign5_1)):-1:1
        for isign4=length(sign4):-1:1
            if (sign4(isign4)<sign5_1(isign5))&(sign4(isign4)>sign5_0(isign5))
                sign4(isign4)=[];
            end           
        end
        Seqset(sign5_0(isign5):sign5_1(isign5))=char(32*ones(1,sign5_1(isign5)-sign5_0(isign5)+1));
    end
    sign23=union(union(sign2,sign3),sign4);
    if length(sign1)>0    
       endsign=find(sign23<sign1);
    else
       endsign=[1:length(sign23)];
    end
    if length(endsign)>0
        Cellnum=union(sign0(1),sign23(endsign));
    else
         Cellnum=[sign0,sign23];
    end
    for i=1:length(Cellnum)-1
        temcell=Seqset(Cellnum(i)+1:Cellnum(i+1)-1);
        temcell(find(isspace(temcell))) = []; 
        temcell(find(temcell=='0')) = []; 
        CellLS{i}=temcell;
    end
    PROFea.CLS=CellLS;
else
    PROFea.S_F=0;
    fclose(fid0);
    PROFea.CLS=[];
end
