% forexample PercentS =[0.0025 0.005 0.01 0.02 0.04]
% Calculate the KNNScores features
function DisKNNScore_set=CountKNNScore_Dis_Vector(DisVector,Labelset,PercentS)
Len=length(Labelset);
DisKNNScore_set=zeros(1,length(PercentS));
[a,b]=sort(DisVector);
for i=1:length(PercentS);    
    NumNegb=ceil(PercentS(i)*Len) ; 
    DisKNNScore_set(1,i)=sum(Labelset(b(1:NumNegb))==1)/NumNegb;
end
