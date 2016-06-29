%%% UPDATED 2/24/16, version 0.1
%%% Changes: 
%%% 1. spatialregionneighborsearch.m -> spatialregionneighborsearchcentroid.m : Added meanlat and meanlon to region list
%%% 2. trackcodematrix.m -> trackcodematrixempty.m : Restructured to handle times when there are no events (incidentally improved initialization). Also incorporated meanlat and meanlon
%%% 3. changed handling of overlapthreshold. Now it must already exist when trackcodematixempty.m is called. 

% Author: Angie Pendergrass. apgrass@uw.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % necessary preparation: cubesphereneighbors
% % % dependencies: spatialregionneighborsearch (identify regions at each time)
% % %               trackecodematrix (track regions through time) 


% % % 1. figure out which pixels are adjacent to which: make a neighor or
% % % connectivity matrix

load cubesphere305small.mat lon lat area pdayend % this has 5 days of global precip data from the end of a simulation
cubesphereneighbors 
% save newneighbors.mat neighbors


%%%%%% Calculate 95th percentile rain.
%%% Option 1. Get it from sorting. 

pday=pdayend; 
psort=sort(pday(:));
pthresh=psort(round(length(psort)*.95));


% %%% Option 2. Calculate the distribution of rain and then interpolate
% pday=pdayend;
% 
% %%% This is a quick and dirty script to calculate the 95th percentile of
% %%% rain.  It cheats about the area, and the code isn't well commented. 
% calc95thpercentile
% pthresh=prrates(13); %%% 95th percentile of rain. 15 or 30 mm/d for the 305 K simulation, depending on how you calculated it.


% % % now identify all the contiguous regions of rainfall above the threshold, separately at each time

tic
%spatialregionneighborsearch
spatialregionneighborsearchcentroid
toc

% clear pday 
% save regionsstatic.mat regiondays regionsdaylist
% load regionsstatic.mat


% % % Finally, link the regions that overlap from one time to the next into "events" to figure out their lifetime. 
% % % There is a crucial parameter, the "overlap," which determines how much overlap there needs to be for events to be counted as continuing. 
% % % There is also splitting and merging of events, and the algorithm for that was interesting to write... 

overlapthreshold = 0.25 ; % recommend: 0.25 for daily data (which was included with the distribution), 0.5 for 6 hourly data (historically more common basis for tracking). 
tic
%trackcodematrix
trackcodematrixempty
toc

% save regiontimes.mat timeregionlist timeregiondays

%regiondays=regiondays(:,1:11);
%regionsdaylist=regionsdaylist(1:11);
%save checkstuff.mat regiondays regionsdaylist 


%%% Uncomment the code below for a few plots to help you digest this. 

% 
% figure(1); 
% clf
% for day=1:size(pdayend,2)
% subplot(5,1,day)
% scatter(lon,sindlat,20,pdayend(:,day),'filled');colorbar; % 15 might be better
% xlim([0 360])
% ylim([-1 1])
% set(gca,'ytick',sind([-90:30:90]),'yticklabel',{'SP';'60 S';'30 S';'EQ';'30 N';'60 N';'NP';})
% title([num2str(length(regionsdaylist{day})) ' regions over 95th percentile of rain'])
% end

%  
% 
% figure(2); 
% clf
% 
% regiondays2=regiondays; 
% regiondays2(regiondays==0)=NaN;
% for day=1:size(regiondays,2)
% subplot(5,1,day)
% scatter(lon,sindlat,20,regiondays2(:,day),'filled');colorbar; % 15 might be better
% xlim([0 360])
% ylim([-1 1])
% set(gca,'ytick',sind([-90:30:90]),'yticklabel',{'SP';'60 S';'30 S';'EQ';'30 N';'60 N';'NP';})
% title([num2str(length(regionsdaylist{day})) ' regions over 95th percentile of rain'])
% end
% 
% % 
% timeregiondays2=timeregiondays; 
% timeregiondays2(timeregiondays==0)=NaN;
% figure(3);clf
% for day=1:size(regiondays,2)
% subplot(5,1,day)
% scatter(lon,sindlat,20,timeregiondays2(:,day),'filled');colorbar; % 15 might be better
% xlim([0 360])
% ylim([-1 1])
% caxis([1 length(timeregionlist)])
% set(gca,'ytick',sind([-90:30:90]),'yticklabel',{'SP';'60 S';'30 S';'EQ';'30 N';'60 N';'NP';})
% title([num2str(length(regionsdaylist{day})) ' regions over 95th percentile of rain'])
% end


figure(4);clf
subplot(3,1,1)
hist([timeregionlist(:).Length])
title(['Length, mean=' num2str(mean([timeregionlist(:).Length]))  ' days'])
subplot(3,1,2)
hist([timeregionlist(:).MeanP])
title(['Mean P, mean=' num2str(mean([timeregionlist(:).MeanP])) ' mm/d'])
subplot(3,1,3)
hist([timeregionlist(:).MeanA])
title(['Mean Area, mean=' num2str(mean([timeregionlist(:).MeanA])) ' km^2'])
length(timeregionlist)
%


