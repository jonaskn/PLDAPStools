function plx = readPlx(filename, readSpikes)  
    if nargin<2
        readSpikes=true;
    end


                [plx.info.openedFileName, plx.info.version, plx.info.freq, plx.info.comment, plx.info.trodalness, plx.info.NPW, plx.info.preThresh, plx.info.spikePeakV, plx.info.spikeADResBits, plx.info.slowPeakV, plx.info.slowADResBits, plx.info.duration, plx.info.dateTime] = plx_information(filename);
                % get some counts
                [tscounts, wfcounts, evcounts, slowcounts] = plx_info(plx.info.openedFileName,1);
                % get some other info about the spike channels
                [~,plx.spikeChannels.filters] = plx_chan_filters(plx.info.openedFileName);
                [~,plx.spikeChannels.gains] = plx_chan_gains(plx.info.openedFileName);
                [~,plx.spikeChannels.threshsholds] = plx_chan_thresholds(plx.info.openedFileName);
                [~,plx.spikeChannels.names] = plx_chan_names(plx.info.openedFileName);
                % gives actual number of units (including unsorted) and actual number of
                % channels plus 1
                %nmaxunits, nmaxchannels] = size( tscounts );   
                [units, channels]=find(tscounts);
                units=units-1;
                channels=channels-1;
                nUnits=length(units);
                
                plx.spikeChannels.units=units;
                plx.spikeChannels.channels=channels;
                % we will read in the timestamps of all units,channels into a two-dim cell
                % array named allts, with each cell containing the timestamps for a unit,channel.
                % Note that allts second dim is indexed by the 1-based channel number.
                % preallocate for speed
%                 readSpikes=true;
                if all(readSpikes)
                    %spikes
                    if length(readSpikes)==1
                        plx.spikeChannels.startIndex=1;
                        plx.spikeChannels.maxIndex=plx.info.duration*plx.info.freq;
                    else
                        plx.spikeChannels.startIndex=readSpikes(1);
                        plx.spikeChannels.maxIndex=readSpikes(2);
                        
                    end
                    plx.spikeChannels.spikes = sparse(1,1,false,plx.spikeChannels.maxIndex-plx.spikeChannels.startIndex+1,nUnits,sum(sum(tscounts)));
                    
                    plx.spikeChannels.spikesTimes=@(x) (x+plx.spikeChannels.startIndex-2)/plx.info.freq;
                    plx.spikeChannels.spikesIndecies=@(x) x*plx.info.freq - plx.spikeChannels.startIndex +2;
                    for iunit = 1:length(units)   % starting with unit 0 (unsorted) 
                        [~,ts]=plx_ts(plx.info.openedFileName, channels(iunit) , units(iunit) );
                        if any((ts*plx.info.freq - round(ts*plx.info.freq))>1e-6)
                            error('check what the time actually is');
                        end
                        ts=round(ts*plx.info.freq);
                        ts(ts<plx.spikeChannels.startIndex | ts>plx.spikeChannels.maxIndex)=[];
                        ts=ts-plx.spikeChannels.startIndex+1;
                        plx.spikeChannels.spikes(ts,iunit)=true;
    %                     for ich = 1:nmaxchannels-1
    %                         if ( tscounts( iunit+1 , ich+1 ) > 0 )
    %                             % get the timestamps for this channel and unit 
    %                             [nts, allts{iunit+1,ich}] = plx_ts(plx.info.openedFileName, ich , iunit );
    %                          end
    %                     end
                    end
                end
                

                [u,evchans] = plx_event_chanmap(plx.info.openedFileName);
                [eventChannels]=evchans(evcounts>0);
                if ismember(257,eventChannels)
                    [~, plx.eventChannels.events, plx.eventChannels.values] = plx_event_ts(plx.info.openedFileName, 257); 
                else
                    %events (only strobed for now)
%                     [eventChannels]=find(evcounts);
                    nevchannels=length(eventChannels);

                    nevs=cell(1,16);
                    tsevs=cell(1,16);
                    if ( nevchannels > 0 ) 
                        % need the event chanmap to make any sense of these
%                         [u,evchans] = plx_event_chanmap(plx.info.openedFileName);
%%fix: not only reading the first 16 maybe
                        for iev = 1:16
%                         for iev = 1:nevchannels %wrong, but ok, today
                            if ( evcounts(iev) > 0 )
                                evch = evchans(iev);
                                if ( evch == 257 )
                                    [nevs{iev}, tsevs{iev}, svStrobed] = plx_event_ts(plx.info.openedFileName, evch); 
                                else
                                    [nevs{iev}, tsevs{iev}, svdummy] = plx_event_ts(plx.info.openedFileName, evch);
                                end
                            end
                        end
                    end

                    if ~all(cellfun(@isempty,nevs))
                        tm=vertcat(tsevs{:});
                        ech=cellfun(@(x,y) ones(x,1)*y, nevs, num2cell(1:length(nevs)),'UniformOutput',false);
                        ech=vertcat(ech{:});

    %%fix:some times time does not round with diff of up to -=3*10^-8
                        sev=sparse(round(tm*plx.info.freq),ech,ones(length(ech),1)); %save???
                        [i2,~]=find(sev);

                        ui2=unique(i2);
    %%fix: some events accoured twice at the exact same time.....why?...
                        ssev=~~full(sev(ui2,:));
                        strob=bin2dec(num2str(fliplr(ssev),'%i'));
                        strobt=ui2/plx.info.freq;
                        plx.eventChannels.events=strobt;
                        plx.eventChannels.values=strob;
                    else
                        plx.eventChannels.events=[];
                        plx.eventChannels.values=[];
                    end
                end