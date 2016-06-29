re=6378.1; %km
ae=4/3*pi*re.^3;
sa=sum(area);


% Find connected regions
% Instead of iterating through each gridpoint above the threshold one by 
% one and if-ing through complicated situations and reassignments, start 
% with the first gridpoint, identify and mark *ALL* of its neighbors, and 
% those points' neighbors, iteratively until there are no more adjacent
% neighbors. Then move onto the next (unconnected) unmarked region. This
% should remove the conflict resolution steps. 

for day=1:size(pday,2)
    
    disp(day)

    thispday=pday(:,day);
        
    region=NaN(size(lon));
    region(thispday<=pthresh)=0;
    
    clear regionlist
    regionlist=[];
    
    ireg=0;
    rainind=find(isnan(region),1);
    while ~isempty(rainind)
        ireg=ireg+1;        
        activeinds=rainind;
        while ~isempty(activeinds)
            region(activeinds)=ireg;
            theneighbors=neighbors(activeinds,:);
            theneighbors=unique(theneighbors(:)); % In case two neighbors are also adjacent
            theneighbors(isnan(theneighbors))=[]; % For the 8 3-connected points
            adj=region(theneighbors); % Should only be 0 or NaN or ireg
            if sum(adj>0 & adj~=ireg)>0
                disp('ERROR: adjacent pixel already defined as something else. quitting.')
                return;break;break;break;
            end
            activeinds=theneighbors(isnan(adj));
        end 
        regionlist(ireg).inds=find(region==ireg);
        rainind=find(isnan(region),1);
    end

    for ireg=1:length(regionlist)
        inds=regionlist(ireg).inds;
        thisa=area(inds);
        regionlist(ireg).Area=sum(thisa)./sa*ae;
        regionlist(ireg).P=sum(thispday(inds).*thisa./sum(thisa));
        regionlist(ireg).meanlat=mean(lat(inds));
        theselons=lon(inds);
        if max(abs(diff(theselons)))>180
            theselons(theselons>180)=theselons(theselons>180)-360;
        end
        thismeanlon=mean(theselons);
        regionlist(ireg).meanlon=thismeanlon;
        if thismeanlon<0
            thismeanlon=thismeanlon+360;
        end
        regionlist(ireg).meanlon=thismeanlon;
    end
    
    regionsdaylist{day}=regionlist;
    regiondays(:,day)=region;
    
end
