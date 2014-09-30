Step=15;
Pos=400;

PvCorrCa3=zeros(Pos,Step);
PvCorrCa3Man=zeros(Pos,Step);
PvCorrCa1Bef=zeros(Pos,Step);
PvCorrCa1BefMan=zeros(Pos,Step);
PvCorrCa1Aft=zeros(Pos,Step);
PvCorrCa1AftMan=zeros(Pos,Step);

for i=1:Step
Map1=PhiTTT;
Map2=squeeze(PhiD(1,:,i,:));

PvCorrCa3Man(:,i)=squeeze(min(1-corr(Map1',Map2'),[],1));
PvCorrCa3(:,i)=squeeze(diag(1-corr(Map1',Map2')));

end

for i=1:Step
Map1=ChiTTTBef;
Map2=squeeze(ChiDBef(1,:,i,:));

PvCorrCa1BefMan(:,i)=squeeze(min(1-corr(Map1',Map2'),[],1));
PvCorrCa1Bef(:,i)=squeeze(diag(1-corr(Map1',Map2')));

end

for i=1:Step
Map1=ChiTTTAft;
Map2=squeeze(ChiDAft(1,:,i,:));

PvCorrCa1AftMan(:,i)=squeeze(min(1-corr(Map1',Map2'),[],1));
PvCorrCa1Aft(:,i)=squeeze(diag(1-corr(Map1',Map2')));

end

%%
figure(6)
hold on
P3=mean(PvCorrCa3);
P3M=mean(PvCorrCa3Man);
P1A=mean(PvCorrCa1Aft);
P1AM=mean(PvCorrCa1AftMan);
P1B=mean(PvCorrCa1Bef);
P1BM=mean(PvCorrCa1BefMan);

    
    scatter(P3(1:15),P3M(1:15),'blue','filled')
    scatter(P1A(1:15),P1AM(1:15),'red','filled')
    scatter(P1B(1:15),P1BM(1:15),'green','filled')



%%
% P3=mean(PvCorrCa3);
% P3M=mean(PvCorrCa3Man);
%%
P1B=mean(PvCorrCa1Bef);
P1BM=mean(PvCorrCa1BefMan);
%%

P1A=mean(PvCorrCa1Aft);
P1AM=mean(PvCorrCa1AftMan);
%%

figure(6)

%plot(P3,'blue')
hold on
%plot(P3M,'blue')
plot(P1B,'green')
plot(P1BM,'green')
plot(P1A,'red')
plot(P1AM,'red')
%%
figure(6)
hold on
plot(mean(PvCorrCa3),'blue')
plot(mean(PvCorrCa3Man),'blue')
%plot(mean(PvCorrCa1Bef),'green')
%plot(mean(PvCorrCa1BefMan),'green')
%plot(mean(PvCorrCa1Aft),'red')
%plot(mean(PvCorrCa1AftMan),'red')
