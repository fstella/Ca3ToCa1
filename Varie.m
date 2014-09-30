imagesc(reshape(PhiTTT(:,245),20,20))


%%
figure(2)
DM3=pdist(PhiFFF,'correlation');
DM1=pdist(ChiFFFAft,'correlation');
[~,IX]=sort(DM3,'ascend');
plot(DM3(IX))
hold on
plot(DM1(IX),'red')
hold off

%%
DensCa1Aft=zeros(400,1);
DensCa1Bef=zeros(400,1);
DensCa3=zeros(400,1);

Cut=0.5;

PhiTTTCut=PhiTTT;
PhiTTTCut(PhiTTTCut<Cut)=0;

ChiTTTAftCut=ChiTTTAft;
ChiTTTAftCut(ChiTTTAftCut<Cut)=0;

ChiTTTBefCut=ChiTTTBef;
ChiTTTBefCut(ChiTTTBefCut<Cut)=0;

for i=1:400
DensCa1Aft(i)=nnz(ChiTTTAftCut(i,:));
DensCa1Bef(i)=nnz(ChiTTTBefCut(i,:));
DensCa3(i)=nnz(PhiTTTCut(i,:));
end






