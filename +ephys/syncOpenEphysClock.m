function [OE2PTBfit, OE2PTB,PTB2OE, maxreconstructionerror ] = syncOpenEphysClock(PDS, filenameE)
%load events only
%filenameE='experiment1.kwe';

timestamps = hdf5read(filenameE, '/event_types/TTL/events/time_samples');
highlow = hdf5read(filenameE, '/event_types/TTL/events/user_data/eventID');
bitNumber = hdf5read(filenameE, '/event_types/TTL/events/user_data/event_channels');

% timestamps=[timestamps]

strobeSet=find(bitNumber==7 & highlow==1);
strobeUnset=find(bitNumber==7 & highlow==0);
strobeUnset=[1; strobeUnset];

value=nan(size(strobeSet));
for iStrobe=1:length(strobeSet)
     ts=timestamps <= timestamps(strobeSet(iStrobe)) & timestamps >= timestamps(strobeUnset(iStrobe)) & bitNumber~=7;
     
     value(iStrobe)=sum(2.^bitNumber(ts) .* highlow(ts));
    
end

% uts=unique(timestamps);
% for iTS=1:length(uts)
%    thesebitNumbers=bitNumber(timestamps==uts(iTS));
%    thesebitValues=highlow(timestamps==uts(iTS));
%    
%    value(iTS)=sum(2.^thesebitNumbers .* thesebitValues);
% end
% 
% value(value==0)=[];
% lagmatrix(value,6)
sixletsOE=fliplr(conv2(value,eye(6)));
sixletsOE=sixletsOE(6:end,:);

sixletsPTB=cellfun(@(X) X.unique_number, PDS.data,'uniformOutput',false);

[~, hind]=ismember(datenum(mod(vertcat(sixletsPTB{:}),2^7)),datenum(sixletsOE));

sixletsOEts=double([timestamps(strobeSet(hind)) timestamps(strobeSet(hind+1)) timestamps(strobeSet(hind+2)) timestamps(strobeSet(hind+3)) timestamps(strobeSet(hind+4)) timestamps(strobeSet(hind+5))]);

sixletsPTBts=cellfun(@(X) X.datapixx.unique_number_time(:,1), PDS.data,'uniformOutput',false);
sixletsPTBts=[sixletsPTBts{:}]';
sixletsDPts=cellfun(@(X) X.datapixx.unique_number_time(:,2), PDS.data,'uniformOutput',false);
sixletsDPts=[sixletsDPts{:}]';


OE2PTBfit=[sixletsOEts(:) ones(numel(sixletsOEts),1)]\sixletsPTBts(:);
OE2PTB=@(x) x*OE2PTBfit(1) + OE2PTBfit(2);
PTB2OE=@(x) (x - OE2PTBfit(2))/OE2PTBfit(1);

% OE2DPfit=[sixletsOEts(:) ones(numel(sixletsOEts),1)]\sixletsDPts(:);
% OE2DP=@(x) x*OE2DPfit(1) + OE2DPfit(2);
% DP2OE=@(x) (x - OE2DPfit(2))/OE2DPfit(1);


%maybe get a reconstruction estimate
% mean(abs(((sixletsDPts(:)-OE2DP(sixletsOEts(:))))))
maxreconstructionerror = max(((sixletsDPts(:)-OE2DP(sixletsOEts(:)))));
