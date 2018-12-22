% To find the similar proteins as input one
%makeblastdb -in uniprot_sprot.fasta -title swissprot -out swissprot -dbtype prot
%seq the sequence of protein
% db the dataset for blasting
function [p,ProteinID]= blastpssm_seq_likeP(seq, db)
heads = 'Lala';p=[];ProteinID=[];
seqs = seq;
% %%the path of swissprot data, for example:  G:/blast/blast2.2.28+/bin/psiblast -db 
cmd = ['G:/blast/blast2.2.28+/bin/psiblast -db ' db ' -query inputtmp.fasta -num_iterations 3 -evalue 0.001 -outfmt 6 -out pssmresult11  -out_ascii_pssm pssmresult'];
    fid_in = fopen('inputtmp.fasta','w');
    fprintf(fid_in,seqs);
    fclose(fid_in);
     system(cmd);
    fid_out = fopen('pssmresult','r');
    for i = 1 : 4
        tline=fgetl(fid_out);
    end
    i = 1; 
    while ischar(tline) && ~isempty(tline)
        A = sscanf(tline,'%d %s %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d');
        p(i,:) = A(3:42)';
        i=i+1;
        tline=fgetl(fid_out);
    end
     fclose(fid_out);     
    fid_out1 = fopen('pssmresult11','r');
    i=1;
    tline1=fgetl(fid_out1) ;    
    index1 = find(tline1=='|');
    if length(index1)>0
        temsign=ProteinID;
        ProteinID{i}=tline1(index1(1)+1:index1(1)+6);
         while (~feof(fid_out1))&(i<21)
                tline1=fgetl(fid_out1);     
                if length(tline1)>50
                    index1 = find(tline1=='|');
                    temsign=ProteinID;
                    ProteinID{i+1}=tline1(index1(1)+1:index1(1)+6);
                    temsign2=ProteinID;
                    if length(temsign)<length(unique(temsign2))
                       ProteinID=temsign2;  
                       i=i+1;
                    else
                       ProteinID=temsign;   
                    end
                end
         end
    end
             fclose(fid_out1);
             
             
             



