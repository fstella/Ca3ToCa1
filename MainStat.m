Ca1(0.001,0.4,0.5,3,1,1,1);
%%
Da=load('Ca1Gamma-3Sp8Conn5Pop3NEnv1S1SC1.mat','ChiTTTAft');
%Da=load('Ca3Pop3NEnv1S1SC1.mat','PhiTTT');
Dat=struct2array(Da);
%%
[Radius,BumpStd]=UnitCorr(Dat);
%%
[PopRadius,PopStd]=PopuCorr(Dat);

[NumField,PeSizeAll,PeActAll]=FieldStat(Dat);

save('Ca1Sp08Act','NumField','PeSizeAll','PeActAll')
save('Ca1Sp08Corr','BumpStd','Radius','PopBumpStd','PopRadius')