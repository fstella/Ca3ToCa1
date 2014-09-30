h1 = figure(8);
h2 = get(h1,'CurrentAxes');
% Note: edit these numbers to change position
% and size of inset plot
h3 = axes('pos',[.500 .175 .35 .35]);








Dist=squareform(pdist(ChiTTTBef,'correlation'));

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
scatter(EDist(:),Dist(:),'green')

Dist=squareform(pdist(PhiTTT,'correlation'));

scatter(EDist(:),Dist(:),'blue')
