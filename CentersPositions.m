XCC=[];
YCC=[];
fa=1;
n=fa*2+1;

[X,Y]=meshgrid(1:22);
[XT,YT]=meshgrid(22/120:22/120:22);
MACA3=zeros(20,20);
str=ones(n);
    str(fa+1,fa+1)=0;

Th= 0.2*max(ChiTTTAft(:));   
    
    
for i=1:1500

    A=reshape(ChiTTTAft(:,i),20,20);
    
%     Ma=(A>Th);
    
 AA=cat(2,A,A(:,1:fa));
 AA=cat(2,A(:,20-fa+1:20),AA);
 AA=cat(1,AA,AA(1:fa,:));
 AA=cat(1,AA(20-fa+1:20,:),AA);
 AA(AA<0.2)=0;
    
 bw = AA > imdilate(AA, str);
 
%  AA=reshape(AA,484,1);
%  
%  XX=X(:);
%  YY=Y(:);
%  f1=fit([XX,YY],AA,'loess','Span',0.2);
%  
%  z=f1(XT,YT);
%  z(z<0.2)=0;
 Ma=bw(1+fa:20+fa,1+fa:20+fa);

 MACA3=MACA3+double(Ma);
 
  [xc,yc]=find(Ma);
  XCC=cat(1,XCC,xc);
  YCC=cat(1,YCC,yc);

end

% kernel=ones(4)/16;
% MAav = conv2(MACA1Aft, kernel, 'same');
h = fspecial('average', 3);
MAavCA1Aft=filter2(h,MACA3);

figure
imagesc(MAavCA1Aft)







