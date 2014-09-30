%%

Map1=PhiTTT;
Map2=PhiFFF;

corr2(Map1,Map2)
nanmean(diag(corr(Map1,Map2)))
nanmean(diag(corr(Map1',Map2')))
%%
Res=size(PhiTTT,1);

SpCorrCa3=zeros(1500,40);
PvCorrCa3=zeros(Res,40);
SpCorrCa1Bef=zeros(2000,40);
PvCorrCa1Bef=zeros(Res,40);
SpCorrCa1Aft=zeros(2000,40);
PvCorrCa1Aft=zeros(Res,40);

%% Correlazione vs. Time per CA3
for i=1:40
Map1=PhiTTT;
Map2=squeeze(PhiD(1,:,i,:));

%SpCorrCa3(:,i)=squeeze(diag(corr(Map1,Map2)));
PvCorrCa3(:,i)=squeeze(diag(corr(Map1',Map2')));

end

%% Correlazione vs. Time per CA1

for i=1:40
Map1=ChiTTTBef;
Map2=squeeze(ChiDBef(1,:,i,:));

%SpCorrCa1Bef(:,i)=(diag(corr(Map1,Map2)));
PvCorrCa1Bef(:,i)=(diag(corr(Map1',Map2')));
end

for i=1:40
Map1=ChiTTTAft;
Map2=squeeze(ChiDAft(1,:,i,:));

%SpCorrCa1Aft(:,i)=(diag(corr(Map1,Map2)));
PvCorrCa1Aft(:,i)=(diag(corr(Map1',Map2')));
end

%% Figura Distanza vs. Correlazione
L=10;
Pos=10;

Dist=squareform(pdist(ChiTTTAft,'correlation'));
figure(8)
[x,y]=meshgrid(L/Pos:L/Pos:L);
x=x(:);
y=y(:);

EX=pdist(x,'euclidean');
EX=squareform(min(EX,L-EX));

EY=pdist(y,'euclidean');
EY=squareform(min(EY,L-EY));


%EE=cat(2,x,y);
EDist=(EY.^2+EX.^2).^(1/2);
hold on
scatter(tril(EDist(:)),tril(Dist(:)),'red')
xlabel('Spatial Distance','FontWeight','bold','FontSize',20);
ylabel('Population Vector Distance','FontWeight','bold','FontSize',20);
set(gca,'FontSize',16,'FontWeight','bold');



Dist=squareform(pdist(PhiTTT,'correlation'));
scatter(tril(EDist(:)),tril(Dist(:)),'blue')


%% Distanza Posizioni Finali CA3 vs. CA1
figure(2)
DM3=pdist(PhiFFF,'correlation');
DM1=pdist(ChiFFFAft,'correlation');
[~,IX]=sort(DM3,'ascend');
plot(DM3(IX))
hold on
plot(DM1(IX),'red')
hold off




