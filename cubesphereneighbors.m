% load cubesphere305small.mat lon lat

sindlat=sind(lat);
theinds=1:length(lon);

box4=reshape(1:4*4,[4 4])';
box3=reshape(1:3*3,[3 3])';
box3row4=reshape(1:4*3,[4 3])';
box4row3=reshape(1:4*3,[3 4])';
box3row2=reshape(1:3*2,[2 3])';
box4row2=reshape(1:4*2,[2 4])';
box2row3=reshape(1:2*3,[3 2])';
box2=reshape(1:2*2,[2 2])';

% face 1. length, 8281, 91x91
% 
% f1inds=1:(91*(15*3*2+1));

f1mat=NaN(91,91);
icount=0;
for ibox=0:21
    f1mat(1:4,ibox*4+(1:4))=box4+icount;
    icount=icount+16; %sum(sum(~isnan(f1mat)));
end
    icount=sum(sum(~isnan(f1mat)));
    ibox=22;
    f1mat(1:4,ibox*4+(1:3))=box4row3+icount;

for irow=0:28
    icount=sum(sum(~isnan(f1mat)));
    f1mat(4+irow*3+(1:3),(1:4))=box3row4+icount;
    icount=icount+12; 
    for ibox=0:28
        f1mat(4+irow*3+(1:3),4+ibox*3+(1:3))=box3+icount;
        icount=icount+9; %sum(sum(~isnan(f1mat)));
    end
end
    


% f2: 90 by 91 

f2mat=NaN(91,90);
icount=91*91;
for ibox=0:29
    f2mat(1:4,ibox*3+(1:3))=box4row3+icount;
    icount=icount+12; %sum(sum(~isnan(f2mat)));
end
for irow=0:28
    for ibox=0:29
        f2mat(4+irow*3+(1:3),ibox*3+(1:3))=box3+icount;
        icount=icount+9; %sum(sum(~isnan(f2mat)));
    end
end

f3mat=f2mat+90*91;


f4mat=NaN(91,89);
icount=91*91+90*91*2;
for ibox=0:28
    f4mat(1:4,ibox*3+(1:3))=box4row3+icount;
    icount=icount+12; %sum(sum(~isnan(f2mat)));
end
ibox=29;
f4mat(1:4,ibox*3+(1:2))=box4row2+icount;
icount=icount+4*2;
for irow=0:28
    for ibox=0:28
        f4mat(4+irow*3+(1:3),ibox*3+(1:3))=box3+icount;
        icount=icount+9; %sum(sum(~isnan(f2mat)));
    end
    ibox=29;
    f4mat(4+irow*3+(1:3),ibox*3+(1:2))=box3row2+icount;
    icount=icount+3*2;
end


f5mat=NaN(89,89);
icount=91*91+90*91*2+91*89;
for irow=0:28
    for ibox=0:28
        f5mat(irow*3+(1:3),ibox*3+(1:3))=box3+icount;
        icount=icount+9; %sum(sum(~isnan(f2mat)));
    end
    ibox=29;
    f5mat(irow*3+(1:3),ibox*3+(1:2))=box3row2+icount;
    icount=icount+3*2;
end
irow=29;
for ibox=0:28
    f5mat(irow*3+(1:2),ibox*3+(1:3))=box2row3+icount;
    icount=icount+6; %sum(sum(~isnan(f2mat)));
end
ibox=29;
f5mat(irow*3+(1:2),ibox*3+(1:2))=box2+icount;

f6mat=f5mat+89*89;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
neighbors=NaN(length(lon),4); % [l r t b]

%%%%%%%%   91    90    90    89
fbodymat=[f1mat f2mat f3mat f4mat];
fbodymattbn =   [[NaN f5mat(end,:) NaN fliplr(f5mat(:,end)') NaN fliplr(f5mat(1,:)) NaN f5mat(:,1)'];...
    fbodymat;...
    [NaN f6mat(1,:) NaN f6mat(:,end)' NaN fliplr(f6mat(end,:)) NaN fliplr(f6mat(:,1)')];];
fbodymatn=[[NaN; f4mat(:,end); NaN;] fbodymattbn [NaN; f1mat(:,1); NaN;]]; 

nr=size(fbodymatn,1);
nc=size(fbodymatn,2);
for ii=2:nr-1
    for jj=2:nc-1
        i=fbodymatn(ii,jj);
        neighbors(i,4)=fbodymatn(ii-1,jj);  % bottom neighbor
        neighbors(i,3)=fbodymatn(ii+1,jj);  % top neighbor
        neighbors(i,1)=fbodymatn(ii,jj-1);  % left neighbor
        neighbors(i,2)=fbodymatn(ii,jj+1);  % right neighbor
    end
end

f6matn=[[NaN f1mat(end,2:end-1) NaN]; flipud(f4mat(end,:)') f6mat f2mat(end,1:end-1)'; [NaN fliplr(f3mat(end,1:end-1)) NaN]; ];
for i=min(f6mat(:)):max(f6mat(:))
    [ii,jj]=find(i==f6matn);
    neighbors(i,4)=f6matn(ii-1,jj);  % bottom neighbor
    neighbors(i,3)=f6matn(ii+1,jj);  % top neighbor
    neighbors(i,1)=f6matn(ii,jj-1);  % left neighbor
    neighbors(i,2)=f6matn(ii,jj+1);  % right neighbor
end

f5matn=[[NaN fliplr(f3mat(1,1:end-1)) NaN]; f4mat(1,:)' f5mat flipud(f2mat(1,1:end-1)');[NaN f1mat(1,2:end-1) NaN]; ];
for i=min(f5mat(:)):max(f5mat(:))
    [ii,jj]=find(i==f5matn);
    neighbors(i,4)=f5matn(ii-1,jj);  % bottom neighbor
    neighbors(i,3)=f5matn(ii+1,jj);  % top neighbor
    neighbors(i,1)=f5matn(ii,jj-1);  % left neighbor
    neighbors(i,2)=f5matn(ii,jj+1);  % right neighbor
end

% save newneighbors.mat neighbors
