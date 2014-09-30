function []=Ca3(Pop,NEnv,Spread,SpreadConn)


L=10;
N=23;
PopFrac=1/Pop;
NTot=1500;
%NTot=500*Pop;
Pos=20;
Step=L/N;
%Spread=1;
%SpreadConn=1;
%NEnv=10;
%NEnv=10;
CNoise=10;
SessionLength=10000;
vel=0.25;

NStep=40;
ConnProb=0.6;
DFrac=0;

sp=0.1*1/(Pop/3);

PhiT=zeros(NEnv,Pos^2,NTot);
PhiF=zeros(NEnv,Pos^2,NTot);
PhiD=zeros(NEnv,Pos^2,NStep,NTot);
PhiExp=zeros(NEnv,SessionLength,NTot);

DisT=zeros(Pos^2,NTot);    
Involved=zeros(NTot,NEnv);
DoubleField=zeros(NTot,NEnv);
XC=zeros(NTot,2,NEnv);
YC=zeros(NTot,2,NEnv);
Conn=zeros(NTot,NTot);
EnergyT=zeros(NEnv,Pos^2);
EnergyD=zeros(NEnv,Pos^2,NStep);

% [X1,Y1]=meshgrid(Step:Step:L);
% X1=X1(:);
% Y1=Y1(:);

% X1=L*rand(NTot,1)+CNoise*randn(NTot,1);
% Y1=L*rand(NTot,1)+CNoise*randn(NTot,1);
% 
% 
% X1=abs(mod(X1,L));
% Y1=abs(mod(Y1,L));

[x,y]=meshgrid(L/Pos:L/Pos:L);
x=x(:);
y=y(:);
NConn=ones(NTot,NTot);
%Peak=rand(N^2,1)*0+1;


    
Involv=cat(1,ones(1000,1),zeros(500,1));
Involved(:,1)=Involv;

Involv=cat(1,zeros(500,1),ones(1000,1));
Involved(:,2)=Involv;

for t=1:NEnv

% Involv=rand(NTot,1);
% Involv(Involv<PopFrac)=0;
% Involv(Involv>=PopFrac)=1;
% Involved(:,t)=Involv;    
    
    
DField=rand(NTot,1);
DField(DField<DFrac)=0;
DField(DField>=DFrac)=1;
DoubleField(:,t)=DField;


X2=L*rand(NTot,2)+CNoise*randn(NTot,2);
Y2=L*rand(NTot,2)+CNoise*randn(NTot,2);


X2=abs(mod(X2,L));
Y2=abs(mod(Y2,L));

XC(:,:,t)=X2;
YC(:,:,t)=Y2;

for q=1:NTot
    if(Involv(q)==0)
    for o=1:q-1
    if(Involv(o)==0)
    
    nq=2-DField(q);
    no=2-DField(o);
        
for k=1:nq    
for p=1:no   
    
DisX2=min(abs(X2(q,k)-X2(o,p)),L-abs(X2(q,k)-X2(o,p)));
DisY2=min(abs(Y2(q,k)-Y2(o,p)),L-abs(Y2(q,k)-Y2(o,p)));
Conn(q,o)=Conn(q,o)+exp(-sqrt(DisX2^2+DisY2^2)/(2*SpreadConn*SpreadConn));
Conn(o,q)=Conn(q,o);
    NConn(q,o)=NConn(q,o)+1;
    NConn(o,q)=NConn(o,q)+1;
    
end
end  
    
    end
    end
    end
    
    
    
    
end


end

Conn=Conn./NConn;

for q=1:NTot
    for o=1:q-1
    Dilu=rand(1);
    if(Dilu>ConnProb)
    Conn(q,o)=0;
        Conn(o,q)=0;
    end
        
    
    end
end



for t=1:NEnv
Involv=Involved(:,t);
    X2=XC(:,:,t);
    Y2=YC(:,:,t);
    %cp=zeros(N^2,1);
%p=randperm(N^2);

% for j=1:100
% I=find(p==j);
% cp(j)=I;
% end

% X2=X1(p);
% Y2=Y1(p);

[x,y]=meshgrid(L/Pos:L/Pos:L);
x=x(:);
y=y(:);



for z=1:Pos^2
Th=0;
    
for q=1:NTot
    if(Involv(q)==0)
    DisX1=min(abs(x(z)-X2(q,1)),L-abs(x(z)-X2(q,1)));
    DisY1=min(abs(y(z)-Y2(q,1)),L-abs(y(z)-Y2(q,1)));
    
    DisT(z,q)=sqrt(DisX1^2+DisY1^2);
    if(DisT(z,q)<2)
    PhiT(t,z,q)=1*exp(-(DisX1^2+DisY1^2)/(2*Spread*Spread));
    else
    PhiT(t,z,q)=0;    
    end
    if(DField(q)==0)
    DisX1=min(abs(x(z)-X2(q,2)),L-abs(x(z)-X2(q,2)));
    DisY1=min(abs(y(z)-Y2(q,2)),L-abs(y(z)-Y2(q,2)));
    
    DisT(z,q)=sqrt(DisX1^2+DisY1^2);
    if(DisT(z,q)<2)
    PhiT(t,z,q)=PhiT(t,z,q)+1*exp(-(DisX1^2+DisY1^2)/(2*Spread*Spread));
    else
     PhiT(t,z,q)=0;   
    end
    end
    
    else
    PhiT(t,z,q)=0;
    end

    
    
end

Act=squeeze(PhiT(t,z,:));

[Actt,Th,ga]=Sparsity(Act,sp,Th);


PhiT(t,z,:)=Actt;

%EnergyT(t,z)=-Actt'*Conn*Actt-sum(Actt)*Th+sum(Actt.^2)*1/(2*ga);


for g=1:NStep
    Actt=Conn*Actt;
    
    if(g<4)
    Actt=Actt+(1-g*0.3).*Act;
    
    end
    
    [Actt,Th,ga]=Sparsity(Actt,sp,Th);

    PhiD(t,z,g,:)=Actt;

%EnergyD(t,z,g)=-Actt'*Conn*Actt-sum(Actt)*Th+sum(Actt.^2)*1/(2*ga);
end
PhiF(t,z,:)=Actt;
z
end 


%Exploration

x=zeros(SessionLength,1);
y=zeros(SessionLength,1);
Dir=zeros(SessionLength,1);

x(1)=rand(1)*L;
y(1)=rand(1)*L;
Dir(1)=rand(1)*2*pi;
Th=0;
for z=1:SessionLength
for q=1:NTot
    if(Involv(q)==0)
    DisX1=min(abs(x(z)-X2(q)),L-abs(x(z)-X2(q)));
    DisY1=min(abs(y(z)-Y2(q)),L-abs(y(z)-Y2(q)));
    
    DisT(z,q)=sqrt(DisX1^2+DisY1^2);
    if(DisT(z,q)<2)
    PhiExp(t,z,q)=1*exp(-(DisX1^2+DisY1^2)/(2*Spread*Spread));
    else
     PhiExp(t,z,q)=0;   
    end
    else
    PhiExp(t,z,q)=0;
    end

end

Act=squeeze(PhiExp(t,z,:));

[Actt,Th]=Sparsity(Act,sp,Th);
PhiExp(t,z,:)=Actt;

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




end


PhiTT=ones(1,NTot);
for t=1:NEnv
PhiTT=cat(1,PhiTT,squeeze(PhiT(t,:,:)));
end

PhiTTT=PhiTT(2:NEnv*Pos^2+1,:);

PhiFF=ones(1,NTot);
for t=1:NEnv
PhiFF=cat(1,PhiFF,squeeze(PhiF(t,:,:)));
end

PhiFFF=PhiFF(2:NEnv*Pos^2+1,:);

%PairDT=squareform(pdist(PhiTTT,'euclidean'));
%PairDF=squareform(pdist(PhiFFF,'euclidean'));
%Tree=linkage(PairD,'median');
%[MultiT,str]=mdscale(PairDT,3,'Criterion','Sammon','Start','cmdscale');
%[MultiF,str]=mdscale(PairDF,3,'Criterion','Sammon','Start','cmdscale');

save(['Ca3Pop' num2str(1/PopFrac) 'NEnv' num2str(NEnv) 'S' num2str(Spread) 'SC' num2str(SpreadConn) '.mat'],'PhiTTT','PhiFFF','PhiD','PhiExp','XC','YC','Involved');
end


function [NPopVec,Th,ga]=Sparsity(PopVec,sp,Th)
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
while(abs(Diff)>sp/10)
if(Diff>0)
mm=Th;
Th=Th+(MM-Th)/2;
else  
MM=Th;
Th=mm+(Th-mm)/2;
end
if(Th<0)
Th=0;
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
ga=0.1/Tot;


end



