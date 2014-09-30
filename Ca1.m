function []=Ca1(Gamma,spratio,ConnProbCa1,Pop,NEnvCa3,Spread,SpreadConn)

load(['Ca3Pop' num2str(Pop) 'NEnv' num2str(NEnvCa3) 'S' num2str(Spread) 'SC' num2str(SpreadConn) '.mat']);

L=10;
sp=0.1*1/(Pop/3);


spCa1=sp*spratio;
NTot=1500;
NTotCa1=2000;


%NTot=500*Pop;
Pos=20;

NStep=40;

SessionLength=5000;
vel=0.25;
NLSess=2;

NEnv=NEnvCa3;

%NEnv=1;

ChiTBef=zeros(NEnv,Pos^2,NTotCa1);
ChiFBef=zeros(NEnv,Pos^2,NTotCa1);
ChiTAft=zeros(NEnv,Pos^2,NTotCa1);
ChiFAft=zeros(NEnv,Pos^2,NTotCa1);
ChiDBef=zeros(NEnv,Pos^2,NStep,NTotCa1);
ChiDAft=zeros(NEnv,Pos^2,NStep,NTotCa1);
DisT=zeros(Pos^2,NTotCa1);    
ChiExp=zeros(NLSess,SessionLength,NTotCa1);
ActFirst=zeros(NTotCa1,Pos^2,1);


%Connections
%ConnCa1=rand(NTotCa1,NTot)*0.1+0.1; %Random or Uniform Valued Connections
ConnCa1=ones(NTotCa1,NTot);

ConnExistCa1=rand(NTotCa1,NTot);

ConnExistCa1(ConnExistCa1>(1-ConnProbCa1))=1;
ConnExistCa1(ConnExistCa1<=(1-ConnProbCa1))=0;

ConnCa1=ConnCa1.*ConnExistCa1;
TotW=sum(ConnCa1,2);
TotW=repmat(TotW,1,NTot);
ConnCa1=ConnCa1./(TotW);


ConnCa1Bef=ConnCa1;





for t=1:NEnv

%Novel Environment

% [x,y]=meshgrid(1:1:L);
% x=x(:);
% y=y(:);

for z=1:Pos^2
Th=0;
Act=squeeze(PhiTTT(z+(t-1)*Pos^2,:));    
Act=ConnCa1*Act';



[Actt,Th]=Sparsity(Act,spCa1,Th);


ChiTBef(t,z,:)=Actt;


for g=1:NStep
    Act=squeeze(PhiD(t,z,g,:))';
    
    Act=ConnCa1*Act';
    
    [Actt,Th]=Sparsity(Act,spCa1,Th);
    ChiDBef(t,z,g,:)=Actt;

end


Act=squeeze(PhiFFF(z+(t-1)*Pos^2,:));    
%Act=squeeze(PhiD(1,z,20,:))'; 
Act=ConnCa1*Act';
[Actt,Th]=Sparsity(Act,spCa1,Th);
ChiFBef(t,z,:)=Actt;
z


end

%Save





%Learning
for q=1:NLSess

ConnVar=zeros(NTotCa1,NTot);

x=zeros(SessionLength,1);
y=zeros(SessionLength,1);
Dir=zeros(SessionLength,1);

x(1)=rand(1)*L;
y(1)=rand(1)*L;
Dir(1)=rand(1)*2*pi;
Th=0;
for z=1:SessionLength


ActP=squeeze(PhiExp(t,z+(q-1)*SessionLength,:));
%MM=mean(ActP);
Act=ConnCa1*ActP;
[Actt,Th]=Sparsity(Act,spCa1,Th);

ChiExp(t,z,:)=Actt;
if(z>20)
%Plasticity
Trace=sum(squeeze(PhiExp(t,z-15:z-1,:)),1)./15;
MM=mean(Actt);
%ConnVar=ConnVar+Gamma*(Actt)*(ActP-Trace')';
ConnVar=ConnVar+Gamma*(Actt-MM)*(ActP-Trace')';
%ConnVar=ConnVar+Gamma*(Actt-MM)*(ActP)';
%ConnVar=ConnVar+Gamma*(Actt)*(ActP)';
%ConnVar=ConnVar+Gamma*(Actt)*(ActP-MM)';
end
%New Position
if(z<SessionLength)
x(z+1)=x(z)+vel*cos(Dir(z));
y(z+1)=y(z)+vel*sin(Dir(z));
x(z+1)=abs(mod(x(z+1),L));
y(z+1)=abs(mod(y(z+1),L));
Dir(z+1)=Dir(z)+0.2*randn;

end
if(mod(z,100)==0)
z
end
end
%Apply Plasticity


ConnCa1=(ConnCa1+ConnVar).*ConnExistCa1;
ConnCa1(ConnCa1<0)=0;

TotW=sum(ConnCa1,2);
TotW=repmat(TotW,1,NTot);
ConnCa1=ConnCa1./(TotW);

end

%FamiliarEnv


for z=1:Pos^2
Th=0;
Act=squeeze(PhiTTT(z+(t-1)*Pos^2,:));    
Act=ConnCa1*Act';



[Actt,Th]=Sparsity(Act,spCa1,Th);


ChiTAft(t,z,:)=Actt;

for g=1:NStep
    Act=squeeze(PhiD(t,z,g,:))';
    
     Act=ConnCa1*Act';
    [Actt,Th]=Sparsity(Act,spCa1,Th);
    ChiDAft(t,z,g,:)=Actt;

end


Act=squeeze(PhiFFF(z+(t-1)*Pos^2,:));    
%Act=squeeze(PhiD(1,z,20,:))'; 
Act=ConnCa1*Act';
[Actt,Th]=Sparsity(Act,spCa1,Th);
ChiFAft(t,z,:)=Actt;
z


end
end
%Save

ChiTT=ones(1,NTotCa1);
for t=1:NEnv
ChiTT=cat(1,ChiTT,squeeze(ChiTBef(t,:,:)));
end

ChiTTTBef=ChiTT(2:NEnv*Pos^2+1,:);

ChiFF=ones(1,NTotCa1);
for t=1:NEnv
ChiFF=cat(1,ChiFF,squeeze(ChiFBef(t,:,:)));
end

ChiFFFBef=ChiFF(2:NEnv*Pos^2+1,:);


ChiTT=ones(1,NTotCa1);
for t=1:NEnv
ChiTT=cat(1,ChiTT,squeeze(ChiTAft(t,:,:)));
end

ChiTTTAft=ChiTT(2:NEnv*Pos^2+1,:);

ChiFF=ones(1,NTotCa1);
for t=1:NEnv
ChiFF=cat(1,ChiFF,squeeze(ChiFAft(t,:,:)));
end

ChiFFFAft=ChiFF(2:NEnv*Pos^2+1,:);





SNot=floor(log(abs(Gamma))/log(10));
save(['Ca1Gamma' num2str(SNot) 'Sp' num2str(spratio*10) 'Conn' num2str(ConnProbCa1*10) 'Pop' num2str(Pop) 'NEnv' num2str(NEnv) 'S' num2str(Spread) 'SC' num2str(SpreadConn) '.mat'],'ConnCa1Bef','ConnCa1','ChiTTTBef','ChiFFFBef','ChiTTTAft','ChiFFFAft','ChiExp','ChiDBef','ChiDAft');




end
function [NPopVec,Th]=Sparsity(PopVec,sp,Th)
k=0;
DD=size(PopVec,1);
MM=max(PopVec);
mm=min(PopVec);
if(Th==0)
Th=mm+(1-sp)*(MM-mm);
end
PopVecA=PopVec-Th;
PopVecA(PopVecA<0)=0;

Tot=sum(PopVecA(:))/DD;
Tot2=sum(PopVecA.^2)/DD;
ASp=Tot^2/Tot2;
Diff=ASp-sp;
while(abs(Diff)>sp/20)
if(Diff>0)
mm=Th;
Th=Th+(MM-Th)/2;
else  
MM=Th;
Th=mm+(Th-mm)/2;
end
PopVecA=PopVec-Th;
PopVecA(PopVecA<0)=0;

Tot=sum(PopVecA(:))/DD;
Tot2=sum(PopVecA.^2)/DD;
ASp=Tot^2/Tot2;
Diff=ASp-sp;

k=k+1;
if(k>1000)
break
end
end


NPopVec=PopVecA.*(0.1/Tot);


end




