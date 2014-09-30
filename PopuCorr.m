function [PopRadius,PopBumpStd]=PopuCorr(Dat)
RadiusTh=0;
CorrRange=10;
SSM=size(Dat);
SSMs=sqrt(SSM(1));
PopCorr=1-squareform(pdist(Dat,'correlation'));
XMa=1:SSMs;
YMa=1:SSMs;

PopRadius=zeros(SSM(1),1);
PopBumpStd=zeros(SSM(1),1);

for j=1:SSM(1)
xj=mod(j-1,SSMs)+1;
yj=floor((j-1)/SSMs)+1;    

MM=reshape(PopCorr(j,:),SSMs,SSMs);

Map=circshift(MM,[SSMs/2+1-xj SSMs/2+1-yj]);



%% Fit: 'Gaussian Fit'.
[xData, yData, zData] = prepareSurfaceData( XMa, YMa, Map );

% Set up fittype and options.
ft = fittype( 'Ac+A1*exp(-((x-11)^2+(y-11)^2)/(2*(A2^2)))', 'independent', {'x', 'y'}, 'dependent', 'z' );
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf -Inf];
opts.StartPoint = [1 2.5 0];
opts.Upper = [Inf Inf Inf];

% Fit model to data.
[fitresult, gof] = fit( [xData, yData], zData, ft, opts );

AcE=fitresult.Ac;
A1E=fitresult.A1;
A2E=fitresult.A2;

syms x;
[Ra]=abs(solve(AcE+A1E*exp(-(x-11)^2/(2*(A2E^2)))==RadiusTh,x)-CorrRange-1);
PopRadius(j)=double(Ra(1));
PopBumpStd(j)=abs(A2E);


end
end