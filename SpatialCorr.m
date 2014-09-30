function[SpatCorrMean,SpatCorrNormMean]=SpatialCorr(Id,name,Co,S,nn)
N=2000;

if(strcmp(name,'Phi'))
N=1000;
end

%load ChiFiring1.mat;
%load ChiFiring2.mat;
N1=load(['Num' name '1' num2str(Id) 'C' Co 'S' S 'nn' nn '.dat']);
N2=load(['Num' name '2' num2str(Id) 'C' Co 'S' S 'nn' nn '.dat']);
load([name 'Firing1' num2str(Id) 'C' Co 'S' S 'nn' nn '.mat']);
FireNorm1=file;
load([name 'Firing2' num2str(Id) 'C' Co 'S' S 'nn' nn '.mat']);
FireNorm2=file;

l=0;
Mix=0;
Mean1=0;
Mean2=0;
Var1=0;
Var2=0;
for n=1:N
    if(N1(n)~=0 && N2(n)~=0)
    l=l+1;
    end
end
SpatCorr=zeros(l,1);
SpatCorrNorm=zeros(l,1);
l=0;
for n=1:N
    Mix=0;
    Mean1=0;
    Mean2=0;
    Var1=0;
    Var2=0;
    if(N1(n)~=0 && N2(n)~=0)
    l=l+1;
    for i=1:100
    for j=1:100
    Mix=Mix+FireNorm1(i,j,n)*FireNorm2(i,j,n);
    Mean1=Mean1+FireNorm1(i,j,n);
    Mean2=Mean2+FireNorm2(i,j,n);
    end
    end
    Mix=Mix/(100*100);
    Mean1=Mean1/(100*100);
    Mean2=Mean2/(100*100);
    SpatCorr(l,1)=Mix-Mean1*Mean2;
    for i=1:100
    for j=1:100
    Var1=Var1+(FireNorm1(i,j,n)-Mean1)*(FireNorm1(i,j,n)-Mean1);
    Var2=Var2+(FireNorm2(i,j,n)-Mean2)*(FireNorm2(i,j,n)-Mean2);
    end
    end
    Std1=sqrt(Var1/(100*100));
    Std2=sqrt(Var2/(100*100));
    SpatCorrNorm(l,1)=SpatCorr(l,1)/(Std1*Std2);
    if(SpatCorr(l,1)==0 && (Std1==0 || Std2==0))
    SpatCorrNorm(l,1)=0;
    end
    end
end
    SpatCorrMean=mean(SpatCorr);
    SpatCorrNormMean=mean(SpatCorrNorm);
