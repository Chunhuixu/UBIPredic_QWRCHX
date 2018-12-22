%makeblastdb -in uniprot_sprot.fasta -title swissprot -out swissprot -dbtype prot
function p = blastpssm_seq(seq, db)
heads = 'Lala';
seqs = seq;
%%
%%path for db  G:/blast/blast2.2.28+/bin/psiblast -db 
cmd = ['G:/blast/blast2.2.28+/bin/psiblast -db ' db ' -query inputtmp.fasta -num_iterations 3 -evalue 0.001 -out_ascii_pssm pssmresult'];
    p = zeros(length(seqs),40);
    fid_in = fopen('inputtmp.fasta','w');
    fprintf(fid_in,seqs);
    fclose(fid_in);
    system(cmd);
    %[~,~]=system(cmd);
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
