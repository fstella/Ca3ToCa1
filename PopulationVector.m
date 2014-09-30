function[PopVectMean,PopVectNormMean]=PopulationVector(Id,name,Co,S,nn)

N=2000;

if(strcmp(name,'Phi'))
N=1000;
end
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

PopVect=zeros(10000,1);
PopVectNorm=zeros(10000,1);
for i=1:100
for j=1:100
    
    Mix=0;
    Mean1=0;
    Mean2=0;
    Var1=0;
    Var2=0;
    
    for n=1:N
    Mix=Mix+FireNorm1(i,j,n)*FireNorm2(i,j,n);
    Mean1=Mean1+FireNorm1(i,j,n);
    Mean2=Mean2+FireNorm2(i,j,n);
    end

    Mix=Mix/N;
    Mean1=Mean1/N;
    Mean2=Mean2/N;
    if(Mean1~=0 && Mean2~=0)
    l=l+1;
    PopVect(l,1)=Mix-Mean1*Mean2;
    for n=1:N
    Var1=Var1+(FireNorm1(i,j,n)-Mean1)*(FireNorm1(i,j,n)-Mean1);
    Var2=Var2+(FireNorm2(i,j,n)-Mean2)*(FireNorm2(i,j,n)-Mean2);
    end
    Std1=sqrt(Var1/N);
    Std2=sqrt(Var2/N);
    PopVectNorm(l,1)=PopVect(l,1)/(Std1*Std2);
    if(PopVect(l,1)==0 && (Std1==0 || Std2==0))
    PopVectNorm(l,1)=0;
    end
    if(Std1==0 && Std2==0)
    PopVectNorm(l,1)=1;
    end
    end
end
end
    PopVectMean=mean(PopVect);
    PopVectNormMean=mean(PopVectNorm);