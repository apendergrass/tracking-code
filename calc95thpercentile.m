
L=2.5e6;
wm2tommd=1./L*3600*24;

nd=size(pday,2);


mp=1500;

nbins=100;
binrlog=linspace(0,log(mp),nbins);
dbinlog=diff(binrlog);
binllog=binrlog-dbinlog(1);
binr=exp(binrlog);
binl=exp(binllog);

dbin=dbinlog(1);
binrlogex=binrlog;
binrend=exp(binrlogex(end));


while length(binrlogex)<140
    binrlogex(end+1)=binrlogex(end)+dbin;
    binrend=exp(binrlogex(end));
end

binrlog=binrlogex;
binllog=binrlog-dbinlog(1);
binr=exp(binrlog);
binl=exp(binllog);



pdaywm2=pday*L*1000;
nd=size(pdaywm2,2);

[n]=histc(pdaywm2(:),[0 binl(1:(end-1)) inf]);
ppdf=n./nansum(n);



binr=exp(binrlog)./L*3600*24;
binl=exp(binllog)./L*3600*24;
bincrates=[0 (binl+binr)]/2*wm2tommd;


% calculate the rain rate as a function of percentile of distribution from a rain frequency distribution
bincrates=bincrates(:);

pbinm=(1-10.^(-.1:-.1:-4)); % percentile bins
cdf=cumsum(ppdf);
% find the monotonically increasing part of the cdf for interpolation
startind=2; % initially, start at bin 2; bin 1 is dry frequency and should be much larger than bin 2
lgz=find(diff(cdf)<=0); % initially, end where the cdf stops going up. there is no more rain in the higher bins.
if isempty(lgz); % in this case, there is rain even in the last bin
    lgz=length(cdf);
else   % otherwise, we need to find the smooth part by force.
    if ~isempty(find(lgz<20)); % in this case, something wonky is happening at light rain rates.
        startind=1+max(lgz(find(lgz<20))); % try again, skipping the light rain rates.
    end
    lgz=lgz(find(lgz>20,1)); % try again to find the smooth part, making sure it's not at the beginning.
    if isempty(lgz);
        lgz=length(cdf); % it may still be the case that all the bins have rain, despite something weird happening at light rain rates.
    end
end
mrrates=interp1(cdf(startind:lgz),startind:lgz,pbinm);% interpolate percentiles onto bin indices from the monotonically increasing part of the distribution
prrates=interp1(1:length(bincrates),bincrates,mrrates);% interpolate bin indices onto rain rates from the monotonically increasing part of the distribution
%         end

pthresh=prrates(13);  %% this is the 95th percentile of rain
