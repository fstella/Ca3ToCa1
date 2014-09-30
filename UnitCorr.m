function [Radius,BumpStd]=UnitCorr(Dat)
%%





RadiusTh=0;



SSM=size(Dat);

Radius=zeros(SSM(2),1);
BumpStd=zeros(SSM(2),1);

for U=1:SSM(2)

    
Map=reshape(Dat(:,U),sqrt(SSM(1)),sqrt(SSM(1)));
MapSh=Map;
if(mean(Map)==0)
continue
end
SS=size(Map);
de=SS(1)/2;

DExt=Map;
DExt=cat(2,DExt,Map(:,1:de));
DExt=cat(2,Map(:,de+1:SS(1)),DExt);
DExt=cat(1,DExt,DExt(1:de,:));
DExt=cat(1,DExt(de+1:SS(1),:),DExt);
Map=DExt;

CC=normxcorr2(Map,Map);
SSC=size(CC);
CorrRange=10;
UniCo=CC((SSC(1)+1)/2-CorrRange:(SSC(1)+1)/2+CorrRange,(SSC(1)+1)/2-CorrRange:(SSC(1)+1)/2+CorrRange);
% figure(1)
% subplot(1,2,1)
% imagesc(MapSh)
% subplot(1,2,2)
% imagesc(UniCo)
% pause(0.5)
% % Initialize arrays to store fits and goodness-of-fit.
% fitresult = cell( 2, 1 );
% gof = struct( 'sse', cell( 2, 1 ), ...
%     'rsquare', [], 'dfe', [], 'adjrsquare', [], 'rmse', [] );

%% Fit: 'Gaussian Fit'.
XCorr=1:CorrRange*2+1;
YCorr=1:CorrRange*2+1;
[xData, yData, zData] = prepareSurfaceData( XCorr, YCorr, UniCo );

% Set up fittype and options.
ft = fittype( 'Ac+A1*exp(-((x-11)^2+(y-11)^2)/(2*(A2^2)))', 'independent', {'x', 'y'}, 'dependent', 'z' );
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf -Inf];
opts.StartPoint = [1 2 0];
opts.Upper = [Inf Inf Inf];

% Fit model to data.
[fitresult, gof] = fit( [xData, yData], zData, ft, opts );

AcE=fitresult.Ac;
A1E=fitresult.A1;
A2E=fitresult.A2;

syms x;
[Ra]=abs(solve(AcE+A1E*exp(-(x-11)^2/(2*(A2E^2)))==RadiusTh,x)-CorrRange-1);
Radius(U)=double(Ra(1));
BumpStd(U)=A2E;
end
end




