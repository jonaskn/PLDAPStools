function [PL2PTBfit,PL2PTB,PTB2PL, maxreconstructionerror ] = syncPlexonClock(PDS, filenameE)
%load events only
plx=ephys.readPlx(filenameE, false);
switchbits=false;
invertedBits=false;
%%
                
b=plx.eventChannels.values;
b=double(typecast(int16(b),'uint16'));

if invertedBits 
    b=bitcmp(b,'uint16');
end

if ~switchbits
    flagBits = log(bitshift(b,-8))/log(2) + 1;
    flagData = mod(b,2^8);
else
    flagBits = log(mod(b,2^8));
    flagData = bitshift(b,-8);
end
flagPlexonTime = double(plx.eventChannels.events);

tstartinds=find(flagBits==PDS.initialParametersMerged.event.TRIALSTART);
tenopinds=find(flagBits==PDS.initialParametersMerged.event.TRIALEND);

flagData2=flagData(isinf(flagBits));
%                 sixlets((isinf(flagBits)),:)=[flagData2(1:end-5),flagData2(2:end-4),flagData2(3:end-3),flagData2(4:end-2),flagData2(5:end-1),flagData2(6:end); ones(5,6)];
%                 tmpdn=cellfun(@(X) X.unique_number, PDS.data,'uniformOutput',false);
%                 [has, hind]=ismember(datenum(mod(vertcat(tmpdn{:}),2^8)),datenum(sixlets));
fivelets=zeros(length(flagBits),5);
fivelets((isinf(flagBits)),:)=[flagData2(1:end-4),flagData2(2:end-3),flagData2(3:end-2),flagData2(4:end-1),flagData2(5:end); ones(4,5)];
tmpdn=cellfun(@(X) X.unique_number, PDS.data,'uniformOutput',false);
sixletsPlexon=[ones(size(fivelets,1),1)*mod(tmpdn{1}(1),2^8) fivelets];
sixletsPLDAPS=datenum(mod(vertcat(tmpdn{:}),2^8));
[has, sixletPlexonInd]=ismember(sixletsPLDAPS,datenum(sixletsPlexon));
matchedPLDAPSTrialNumbers=find(has);

flagDatapixxTrialStartTime=nan(max(matchedPLDAPSTrialNumbers),1);
flagComputerTrialStartTime=nan(max(matchedPLDAPSTrialNumbers),1);
flagPlexonTrialStartTime=nan(max(matchedPLDAPSTrialNumbers),1);
for imatchedTrial=matchedPLDAPSTrialNumbers'
    flagDatapixxTrialStartTime(imatchedTrial)=PDS.data{imatchedTrial}.timing.datapixxTRIALSTART(2);
    flagComputerTrialStartTime(imatchedTrial)=PDS.data{imatchedTrial}.timing.datapixxTRIALSTART(1);

    cadidatePlexonTrialStartInd=tstartinds(find(tstartinds>=min(sixletPlexonInd(imatchedTrial)),1,'first'));
    %check if trialNumbers match
    if mod(imatchedTrial(1),2^8)==flagData(cadidatePlexonTrialStartInd)
        flagPlexonTrialStartTime(imatchedTrial)=flagPlexonTime(cadidatePlexonTrialStartInd);
    else
        %display(imatchedTrial);
    end

end

validField=~isnan(flagPlexonTrialStartTime+flagDatapixxTrialStartTime+flagComputerTrialStartTime);

PL2PTBfit=[flagPlexonTrialStartTime(validField) flagPlexonTrialStartTime(validField)*0+1]\flagComputerTrialStartTime(validField);
PL2PTB=@(x) x*PL2PTBfit(1) + PL2PTBfit(2);
PTB2PL=@(x) (x - PL2PTBfit(2))/PL2PTBfit(1);

reconstructionerror=flagPlexonTrialStartTime(validField)-PTB2PL(flagComputerTrialStartTime(validField));
%                 plot(reconstructionerror)
maxreconstructionerror=max(abs(reconstructionerror))*1000; %ms

%                 %
%                 PL2DPfit=[flagPlexonTrialStartTime(validField) flagPlexonTrialStartTime(validField)*0+1]\flagDatapixxTrialStartTime(validField);
%                 PL2DP=@(x) x*PL2DPfit(1) + PL2DPfit(2);
%                 DP2PL=@(x) (x - PL2DPfit(2))/PL2DPfit(1);
%                 
%                 reconstructionerrorDP=flagPlexonTrialStartTime(validField)-DP2PL(flagDatapixxTrialStartTime(validField));
%                 plot(reconstructionerrorDP)
%                 max(abs(reconstructionerrorDP))*1000 %ms