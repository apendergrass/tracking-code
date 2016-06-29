%%% choose overlap parameter. I think the industry standard for 6-h data is
%%% 0.5 (eg Hodges? see literature list) but this is daily rainfall data... 
% overlapthreshold=0.25;

%%% set constants
re=6378.1; %km
ae=4/3*pi*re.^3;
sa=sum(area);
%%%

% INITIALIZE 
timeregiondays=NaN(size(regiondays));
timeregiondays(regiondays==0)=0;

clear timeregionlist
timeregionlist=[];

% Initialize
iday1=1;
iday=iday1;
%
todayregionlist=regionsdaylist{iday1};
if ~isempty(todayregionlist)
    todayregion=regiondays(:,iday);
    for iregt=1:length(todayregionlist)
        ireg=iregt;
        timeregionlist(iregt).when(1)=iday;
        timeregionlist(iregt).where{1}=[ireg];
        timeregionlist(iregt).inds{1}=todayregionlist(ireg).inds;
        timeregionlist(iregt).Area(1)=todayregionlist(ireg).Area;
        timeregionlist(iregt).P(1)=todayregionlist(ireg).P;
        timeregionlist(iregt).meanlat(1)=todayregionlist(ireg).meanlat;
        timeregionlist(iregt).meanlon(1)=todayregionlist(ireg).meanlon;
    end
    timeregiondays(:,iday)=todayregion;
end

% TRACK EACH DAY
for iday=iday1:(size(regiondays,2)-1)
    disp(iday)
    todayregionlist=regionsdaylist{iday};
    nextdayregionlist=regionsdaylist{iday+1};
    % tomorrow's regions.
    if isempty(todayregionlist)                    % if today AND tomorrow are empty, SKIP.
        if ~isempty(nextdayregionlist(:))          % if today IS empty but tomorrow is NOT, get ready to initialize
            nextdayregion=regiondays(:,iday+1);
            % give birth to new regions
            for i=1:length(nextdayregionlist)
                jreg=i;
                iregt=length(timeregionlist)+1;
                timeregionlist(iregt).when(1)=iday+1;
                timeregionlist(iregt).where{1}=[jreg];
                timeregionlist(iregt).inds{1}=nextdayregionlist(jreg).inds;
                timeregionlist(iregt).Area(1)=nextdayregionlist(jreg).Area;
                timeregionlist(iregt).P(1)=nextdayregionlist(jreg).P;
                timeregionlist(iregt).meanlat(1)=nextdayregionlist(jreg).meanlat;
                timeregionlist(iregt).meanlon(1)=nextdayregionlist(jreg).meanlon;
                timeregiondays(nextdayregionlist(jreg).inds,iday+1)=iregt;
            end
        end
    else
        todayregion=regiondays(:,iday);
        todayarea=[todayregionlist(:).Area].*sa./ae;
        nextdayregion=regiondays(:,iday+1);
        if isempty(nextdayregionlist(:))           % if today is NOT empty but tomorrow IS, kill all regions.
            % kill off all regions
            for i=1:length(todayregionlist)
                ireg=i;
                todayinds=todayregionlist(ireg).inds;
                iregt=timeregiondays(todayinds(1),iday);
                timeregionlist(iregt).Length=length(timeregionlist(iregt).when);
                timeregionlist(iregt).MeanA=mean(timeregionlist(iregt).Area);
                timeregionlist(iregt).MeanP=mean(timeregionlist(iregt).P);
            end
        else
            
            
            nextdayarea=[nextdayregionlist(:).Area].*sa./ae;
            
            [Na,Ta]=meshgrid(nextdayarea,todayarea);
            
            minarea=min(Ta,Na);
            
            % fill out the matrix of all the overlapping areas.
            overlapmatrix=zeros(size(minarea));
            timeregind=NaN(size(todayarea));
            for ireg=1:length(todayarea)
                todayinds=todayregionlist(ireg).inds;
                timeregind(ireg)=timeregiondays(todayinds(1),iday);
                tomorrowregions=nextdayregion(todayinds);
                tomorrowareas=area(todayinds);
                checkconnectregions=unique(tomorrowregions(tomorrowregions>0));
                if ~isempty(checkconnectregions)
                    for icon=1:length(checkconnectregions)
                        jreg=checkconnectregions(icon);
                        overlaparea=sum(tomorrowareas(tomorrowregions==jreg));
                        overlapmatrix(ireg,jreg)=overlaparea;
                    end
                end
            end
            
            % set insufficiently connected regions to zero.
            % check before and after: sum(overlapmatrix(:)>0)
            connected=(overlapmatrix./minarea)>overlapthreshold;
            overlapmatrix(connected==0)=0;
            
            %imagesc(overlapmatrix./minarea); colorbar
            
            % Next find and reconcile merges and splits.
            % merges: sum(overlapmatrix>0,1)>1
            % splits: sum(overlapmatrix>0,2)>1
            
            % Merges, accounting for convolved splits.
            merges=sum(overlapmatrix>0,1)>1;
            mergeinds=find(merges);
            while ~isempty(mergeinds)
                % disp([num2str(length(mergeinds)) ' merges remaining'])
                mergetotomorrowind=mergeinds(1);
                thesemerges=find(overlapmatrix(:,mergetotomorrowind)>0);
                thesesplits=[];
                for mi=1:length(thesemerges)
                    thismerge=thesemerges(mi);
                    thissplitset=find(overlapmatrix(thismerge,:)>0);
                    thesesplits(end+(1:length(thissplitset)))=thissplitset;
                    getmergesfrom=thissplitset(merges(thissplitset));
                    for icon=1:length(getmergesfrom) % there will be at least one with the col we're looking at
                        newmerges=find(overlapmatrix(:,getmergesfrom(icon))>0);
                        thesemerges(end+(1:length(newmerges)))=newmerges;
                    end
                end
                thesemerges=unique(thesemerges);
                thesesplits=unique(thesesplits);
                overlapsubset=overlapmatrix(thesemerges,thesesplits);
                % go through the subset and choose the winning ones iteratively
                mergesub=find(sum(overlapsubset>0,1)>1);
                while ~isempty(mergesub)
                    %	disp(['subsets remaining: ' num2str(length(mergesub))])
                    %	break; break; break
                    winningarea=max(max(overlapsubset));
                    [isubwin,jsubwin]=find(overlapsubset==winningarea);
                    if length(isubwin)>1
                        % break ties by choosing the longest-lived event
                        winlength=zeros(size(isubwin));
                        for icon=1:length(isubwin)
                            thistimeind=timeregind(thesemerges(isubwin(icon)));
                            winlength(icon)=length(timeregionlist(thistimeind).when);
                        end
                        % now there's just one winner.
                        isubwin=isubwin(find(winlength==max(winlength),1));
                        jsubwin=jsubwin(find(winlength==max(winlength),1));
                    end
                    overlapsubset(isubwin,:)=0;  % could keep track of these.
                    overlapsubset(:,jsubwin)=0;
                    overlapsubset(isubwin,jsubwin)=winningarea;
                    
                    if sum(jsubwin==mergesub)==0
                        % this is the case where the winning area was in some other
                        % column and removed one member of the merge but may not
                        % have solved it. furthermore, that one didn't remove
                        % itself from the game. this causes an endless loop.
                        %                 disp('Reached skip merge sub flag.')
                        %                 disp(length(mergeinds))
                        %                 return; break; break; break;
                        mergesub=[];
                    else
                        mergesub=find(sum(overlapsubset>0,1)>1);
                    end
                    %            mergesub=find(sum(overlapsubset>0,1)>1);
                end
                % all merges should be reconciled. put it back where it came from.
                overlapmatrix(thesemerges,thesesplits)=overlapsubset;
                merges=sum(overlapmatrix>0,1)>1;
                mergeinds=find(merges);
            end
            % Remaining splits, which should be easier because there are no more merges
            % left.
            splits=sum(overlapmatrix>0,2)>1;
            splitinds=find(splits);
            while ~isempty(splitinds)
                splitfromtoday=splitinds(1);
                thesesplits=find(overlapmatrix(splitfromtoday,:)>0);
                theseoverlaps=overlapmatrix(splitfromtoday,thesesplits);
                winningarea=max(theseoverlaps);
                thewinner=find(theseoverlaps==winningarea,1); % the one with the most areal overlap wins. break ties by taking the first.
                overlapmatrix(splitfromtoday,thesesplits)=0; % could keep track of these.
                overlapmatrix(splitfromtoday,thesesplits(thewinner))=winningarea;
                splits=sum(overlapmatrix>0,2)>1;
                splitinds=find(splits);
            end
            
            % Now we have a unique set of continuing region indices. three possibilities: continuation, death of today's, birth of tomorrow's.
            [icontinue,jcontinue]=find(overlapmatrix>0); % continue from icontinue to jcontinue
            ideath=find(sum(overlapmatrix>0,2)==0);
            jbirth=find(sum(overlapmatrix>0,1)==0);
            
            % kill off old regions
            for i=1:length(ideath)
                ireg=ideath(i);
                iregt=timeregind(ireg);
                timeregionlist(iregt).Length=length(timeregionlist(iregt).when);
                timeregionlist(iregt).MeanA=mean(timeregionlist(iregt).Area);
                timeregionlist(iregt).MeanP=mean(timeregionlist(iregt).P);
            end
            
            % continue the continuing regions
            for i=1:length(icontinue)
                ireg=icontinue(i);
                iregt=timeregind(ireg);
                jreg=jcontinue(i);
                timeregionlist(iregt).when(end+1)=iday+1;
                timeregionlist(iregt).where{end+1}=[jreg];
                timeregionlist(iregt).inds{end+1}=nextdayregionlist(jreg).inds;
                timeregionlist(iregt).Area(end+1)=nextdayregionlist(jreg).Area;
                timeregionlist(iregt).P(end+1)=nextdayregionlist(jreg).P;
                timeregionlist(iregt).meanlat(1)=nextdayregionlist(jreg).meanlat;
                timeregionlist(iregt).meanlon(1)=nextdayregionlist(jreg).meanlon;                
                timeregiondays(nextdayregionlist(jreg).inds,iday+1)=iregt;
            end
            
            % give birth to new regions
            for i=1:length(jbirth)
                jreg=jbirth(i);
                iregt=length(timeregionlist)+1;
                timeregionlist(iregt).when(1)=iday+1;
                timeregionlist(iregt).where{1}=[jreg];
                timeregionlist(iregt).inds{1}=nextdayregionlist(jreg).inds;
                timeregionlist(iregt).Area(1)=nextdayregionlist(jreg).Area;
                timeregionlist(iregt).P(1)=nextdayregionlist(jreg).P;
                timeregionlist(iregt).meanlat(1)=nextdayregionlist(jreg).meanlat;
                timeregionlist(iregt).meanlon(1)=nextdayregionlist(jreg).meanlon;                
                timeregiondays(nextdayregionlist(jreg).inds,iday+1)=iregt;
            end
        end
    end
end

% FINALIZE 
% get summary statistics for the regions that are undead at the end
for iregt=1:length(timeregionlist)
    if ~exist('timeregionlist(iregt).Length')
        timeregionlist(iregt).Length=length(timeregionlist(iregt).when);
        timeregionlist(iregt).MeanA=mean(timeregionlist(iregt).Area);
        timeregionlist(iregt).MeanP=mean(timeregionlist(iregt).P);
    end
end

% Object-oriented version notes
% split into three parts:
%  1. initialize 
%  2. finalize 
%  3. track one day (iday+1) 
%     Inputs : todayregionlist, todayregion,
%     nextdayregionlist, nextdayregion, timeregionlist,
%     timeregiondays(:,iday)
%     Outputs: timeregionlist (updated), timeregiondays(:,iday+1) 
%  The main code change would be moving timeregiondays around.


