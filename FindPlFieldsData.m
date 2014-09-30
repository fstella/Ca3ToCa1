function [PeAct,PeSize,DBW,DT,CC]=FindPlFieldsData(FiringRate,Thh,flag)

%Thh=0.66;

G = fspecial('gaussian',[4 4],1.5);

DT = imfilter(FiringRate,G,'same');

%DT=wiener2(FiringRate,[3 3]);
%DT=RM(:,:,48);

dd=size(DT);
de=dd(1)/2;

% flag=0; %% Set to 1 IF THE ENVIRONMENT IS A TORUS 


if(flag==1)
DExt=DT;
DExt=cat(2,DExt,DT(:,1:de));
DExt=cat(2,DT(:,de+1:dd(1)),DExt);
DExt=cat(1,DExt,DExt(1:de,:));
DExt=cat(1,DExt(de+1:dd(1),:),DExt);
D=DExt;
else
D=DT;    
end

DBW=zeros(size(D));
MM=mean(D(D>0));
%D(D<MM/5)=0;
D(D<MM*Thh)=0;
DBW(D~=0)=1;
CC=bwconncomp(DBW);
%% Number of Fields
NumFields=CC.NumObjects;
AverageAct=zeros(NumFields,1);
PeakAct=zeros(NumFields,1);
PeakActCorr=zeros(NumFields,1);
LisPeakAct=zeros(NumFields,1);
PosPeakAct=zeros(NumFields,1);

%% Average Activity in Fields, Value and Position of Peak Activity in the Field
for i=1:NumFields
AverageAct(i)=mean(D(CC.PixelIdxList{i}));
[PeakActCorr(i) LisPeakAct(i)]=max(D(CC.PixelIdxList{i}));
PosPeakAct(i)=CC.PixelIdxList{i}(LisPeakAct(i));
PeakAct(i)=D(PosPeakAct(i));
end
if(flag==1)
SS=size(D);
[P2PY,P2PX]=ind2sub(SS,PosPeakAct);
Upb=dd(1)+de;
Dob=dd(1)-de;
II=find(P2PX<=Upb & P2PX>Dob & P2PY<=Upb & P2PY>Dob);
NumFields=size(II,2);
AvAct=AverageAct(II);
PeActt=PeakAct(II);
PosPeakX=P2PX(II)-de;
PosPeakY=P2PY(II)-de;

SizeField = cellfun(@numel,CC.PixelIdxList);
Size=SizeField(II);

PeSize=Size(Size>2);
PeAct=PeActt(Size>2);

PosCol=sub2ind(SS,P2PY(II),P2PX(II));

%DBW(PosCol)=2;
%imagesc(DBW)

else
    SS=size(D);
    AvAct=AverageAct;

[PosPeakX,PosPeakY]=ind2sub(SS,PosPeakAct);
Size=cellfun(@numel,CC.PixelIdxList);
PeActCorr=PeakActCorr(Size>2);
PeAct=PeakAct(Size>2);
PePos=PosPeakAct(Size>2);
PeSize=Size(Size>2);
PosCol=sub2ind(SS,PosPeakX,PosPeakY);

DBW(PosCol)=2;
%figure(1)
%imagesc(DBW)

end
%Pos=zeros(size(PosPeakX,1),3);
%Pos(:,1)=PosPeakX;
%Pos(:,2)=PosPeakY;


%Th=find(Size>0);
% figure(2)
% scatter(Pos(Th,1),Pos(Th,2));
