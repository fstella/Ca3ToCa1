n=1;
Phi=reshape(PhiTTT,400,2,1500);
Phi2=squeeze(Phi(:,n,:));
sum(Phi2'>0)
hist(mean(Phi2))
%%
n=2;
Phi=reshape(ChiTTTAft,400,2,2000);
Phi2=squeeze(Phi(:,n,:));
hist(mean(Phi2),100)
sum(Phi2'>0)


%%
NEnv=2;
MM=zeros(2000,1);
for n=1:1
Phi2=squeeze(Phi(:,n,:));
MM0=double(mean(Phi2)>0.05)';

end

for n=1:NEnv
Phi2=squeeze(Phi(:,n,:));
MM=MM+double(mean(Phi2)>0.05)';

end
hist(MM,0:10)
hold on
for n=0:10
Poi(n+1)=(10*1/3)^n*exp(-10*1/3)/factorial(n);

end
plot(0:10,Poi*1500);
hold off

