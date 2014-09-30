function [NumField,PeSizeAll,PeActAll]=FieldStat(Dat)



SSM=size(Dat);
NumField=zeros(SSM(2),1);
PeActAll=[];
PeSizeAll=[];

for U=1:SSM(2)
Map=reshape(Dat(:,U),sqrt(SSM(1)),sqrt(SSM(1)));

if(mean(Map)==0)
continue
end
SS=size(Map);
de=SS(1)/2;

[PeAct,PeSize,DBW,DT,CC]=FindPlFieldsData(Map,0.33,1);

NumField(U)=size(PeAct,1);
PeActAll=cat(1,PeActAll,PeAct);
PeSizeAll=cat(2,PeSizeAll,PeSize);



end
end