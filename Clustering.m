PhiM=cat(1,PhiTTT,squeeze(PhiFFF));
%PhiM=cat(1,ChiTTTAft,ChiFFFAft);
PairDM=squareform(pdist(PhiM,'correlation'));
Pair=PairDM(1:100,101:200);

[~,MI]=min(Pair,[],1);

HH=hist(MI,1:100);
nnz(HH)


for i=1:100
xf(i)=mod(MI(i),10);
yf(i)=floor(MI(i)/10)+1;

xi(i)=mod(i,10);
yi(i)=floor(i/10);

DrX(i)=min(abs(xf(i)-xi(i)),10-abs(xf(i)-xi(i)));
DrY(i)=min(abs(yf(i)-yi(i)),10-abs(yf(i)-yi(i)));

Dr(i)=sqrt(DrY(i)^2+DrX(i)^2);

end

mean(Dr)
Cl=zeros(100,1);
ClM=0;
Dc=zeros(100,100);
Di=zeros(100,100);
DReal=zeros(8,1);
DAtt=zeros(8,1);
DCount=zeros(8,1);


for i=1:100

    for j=1:100
    if(j~=i)
        DcX=min(abs(xf(i)-xf(j)),10-abs(xf(i)-xf(j)));
        DcY=min(abs(yf(i)-yf(j)),10-abs(yf(i)-yf(j)));
        
        DiX=min(abs(xi(i)-xi(j)),10-abs(xi(i)-xi(j)));
        DiY=min(abs(yi(i)-yi(j)),10-abs(yi(i)-yi(j)));
        
        Dc(i,j)=(DcY^2+DcX^2);
        Di(i,j)=(DiY^2+DiX^2);
        
        DR=floor(Di(i,j).^(1/2));
        DReal(DR)=DReal(DR)+Di(i,j).^(1/2);
        DAtt(DR)=DAtt(DR)+Dc(i,j).^(1/2);
        DCount(DR)=DCount(DR)+1;
        
        Cl(i)=Cl(i)+exp(-Dc(i,j));
    
    
    end
    end

ClM=ClM+Cl(i);
end
ClM=ClM/(100*99)
DReal=DReal./DCount;
DAtt=DAtt./DCount;





