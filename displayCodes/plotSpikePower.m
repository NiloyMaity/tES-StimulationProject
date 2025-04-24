%% Please keep the Stimulation Specific Protocol List in MATLAB path.
% Area=V1 or V4, so AreaFlag= 1 or 2 corresponds to V1 & V4
% StimulationType={'tDCS','tACS'}
% condition={'Stim','Sham'};
% Polarity={'Cathodal','Anodal' or 'SG','FG','Alpha'};
% Session={'single','dual','dual60'};

%% This code plots Line plots: delta Firing rate plot and delta band power plot

function plotSpikePower(plotHandle,folderIn,session,normalizeFlag,transientFlag,FRFlag,PSDFlag,bandFlag,dates,protocols,sfVals,conVals,oriVals,colCode,titleString,stimBlockID)
baselineS  = [-0.5 0];          % Baseline period for computing change
stimPeriod = [0.25 0.75];
transientPeriod= [0 0.1];
binWidthMS=10;

for ses=1:length(stimBlockID)
    colCode{1,stimBlockID(ses)}=[];
    titleString{1,stimBlockID(ses)}=[];
end

Col=colCode(~cellfun('isempty',colCode));
colCode=reshape(Col,[1 length(colCode)-length(stimBlockID)]);
% Collecting Data
for c=1:2
    for q=1:length(dates{1,c})
        for r=1:length(stimBlockID)
            protocols{1,c}{q,stimBlockID(r)}=[];
        end
    end
    % clearing the stimulation block from protocol
    y=protocols{1,c}(~cellfun('isempty',protocols{1,c}));
    protocol=reshape(y,[length(dates{1,c}) length(protocols{1,c})-length(stimBlockID)]);

    % data collection loop
    for day=1:length(dates{1,c})
        for j=1:length(sfVals)
            for k=1:length(oriVals)
                for l=1:length(conVals)
                    DataBand = dir(fullfile(folderIn{1,c},append(num2str(sfVals(j)),'SF','_',num2str(oriVals(k)),'Ori','_',num2str(conVals(l)),'Con','_',dates{1,c}{1,day},'*BandData.mat')));
                    clear BandValues
                    BandValues = load(fullfile(folderIn{1,c},DataBand.name));
                    cSG=BandValues.CollectSG;
                    cFG=BandValues.CollectFG;
                    cBlSG=BandValues.CollectBlSG;
                    cStSG=BandValues.CollectStSG;
                    cBlFG=BandValues.CollectBlFG;
                    cStFG=BandValues.CollectStFG;
                    for r=1:length(stimBlockID)
                        cSG(stimBlockID(r),:)=0;
                        cFG(stimBlockID(r),:)=0;
                        cBlSG(stimBlockID(r),:)=0;
                        cStSG(stimBlockID(r),:)=0;
                        cBlFG(stimBlockID(r),:)=0;
                        cStFG(stimBlockID(r),:)=0;
                    end

                    cBlSG=reshape(nonzeros(cBlSG),[size(protocol,2) size(cBlSG,2)]);
                    cStSG=reshape(nonzeros(cStSG),[size(protocol,2) size(cStSG,2)]);
                    cBlFG=reshape(nonzeros(cBlFG),[size(protocol,2) size(cBlFG,2)]);
                    cStFG=reshape(nonzeros(cStFG),[size(protocol,2) size(cStFG,2)]);

                    for iProt =1:size(protocol,2)
                        dataSpike = dir(fullfile(folderIn{1,c},append(num2str(sfVals(j)),'SF','_',num2str(oriVals(k)),'Ori','_',num2str(conVals(l)),'Con','_',protocol{day,iProt},'_',dates{1,c}{1,day},'*PSTH.mat')));
                        % Loading the spike file
                        if ~isempty(dataSpike)
                            s1 = load(fullfile(folderIn{1,c},dataSpike.name)); % DataStructure(field).name
                            rawPSTH{iProt,j,k,l} = s1.PSTHGrid;
                            spktList=s1.xs; %Spike timevals
                        else
                            rawPSTH{iProt,j,k,l} =[];
                        end

                        if ~isempty(dataSpike)
                            blPos = find(spktList>=baselineS(1),1)+ (1:(diff(baselineS))/(binWidthMS/1000));
                            stPos = find(spktList>=stimPeriod(1),1)+ (1:(diff(stimPeriod))/(binWidthMS/1000));
                            transientPos=find(spktList>=transientPeriod(1),1)+ (1:(diff(transientPeriod))/(binWidthMS/1000));

                            rawFRBl{iProt,j,k,l}= mean(s1.PSTHGrid(:,blPos),2);
                            rawFRSt{iProt,j,k,l}= mean(s1.PSTHGrid(:,stPos),2);
                            rawFRTransient{iProt,j,k,l}= mean(s1.PSTHGrid(:,transientPos),2);
                        end

                        SGAcrossBl{iProt,j,k,l}=cBlSG(iProt,:)';
                        SGAcrossSt{iProt,j,k,l}=cStSG(iProt,:)';
                        FGAcrossBl{iProt,j,k,l}=cBlFG(iProt,:)';
                        FGAcrossSt{iProt,j,k,l}=cStFG(iProt,:)';
                    end
                end
            end
        end
        for iProt =1:size(protocol,2)
            if ~isempty(dataSpike)
                rawFRAcrossBl{day,iProt,:}=mean(cell2mat(rawFRBl(iProt,:)),2);
                rawFRAcrossSt{day,iProt,:}=mean(cell2mat(rawFRSt(iProt,:)),2);
                rawFRAcrossTransient{day,iProt,:}=mean(cell2mat(rawFRTransient(iProt,:)),2);
            end

            % Averaging across conditions
            bandSGAcrossBl{day,iProt,:}=mean(cell2mat(SGAcrossBl(iProt,:)),2);
            bandFGAcrossBl{day,iProt,:}=mean(cell2mat(FGAcrossBl(iProt,:)),2);
            bandSGAcrossSt{day,iProt,:}=mean(cell2mat(SGAcrossSt(iProt,:)),2);
            bandFGAcrossSt{day,iProt,:}=mean(cell2mat(FGAcrossSt(iProt,:)),2);
        end
    end

    %% Here pulling of all data from required dates has been completed.
    clear SGAcrossAll FGAcrossAll SGAcrossBlAll SGAcrossStAll FGAcrossBlAll FGAcrossStAll FRAcrossBlAll FRAcrossStAll FRAcrossTransientAll
    for b=1:size(protocol,2)
        if FRFlag==1
            FRAcrossBlAll{b,:}=cat(1,rawFRAcrossBl{:,b});
            FRAcrossStAll{b,:}=cat(1,rawFRAcrossSt{:,b});
            FRAcrossTransientAll{b,:}=cat(1,rawFRAcrossTransient{:,b});
        end
        SGAcrossBlAll{b,:}=cat(1,bandSGAcrossBl{:,b});
        FGAcrossBlAll{b,:}=cat(1,bandFGAcrossBl{:,b});
        SGAcrossStAll{b,:}=cat(1,bandSGAcrossSt{:,b});
        FGAcrossStAll{b,:}=cat(1,bandFGAcrossSt{:,b});
    end
    clear SGAcross FGAcross SGAcrossBl FGAcrossBl SGAcrossSt FGAcrossSt rawFRAcrossBl rawFRAcrossSt rawFRAcrossTransient BandSGAcross BandFGAcross bandSGAcrossBl bandFGAcrossBl bandSGAcrossSt bandFGAcrossSt
    % And This brings all electrode in one page

    %% Putting all electrodes' band specific data in a cell, so that we can later do subtraction between conditions (Stim-Sham)
    if FRFlag==1
        compareCellFRBl{c,:}=FRAcrossBlAll;
        compareCellFRSt{c,:}=FRAcrossStAll;
        compareCellFRTransient{c,:}=FRAcrossTransientAll;
    end
    separateCellSGBl{c,:}=SGAcrossBlAll;
    separateCellSGSt{c,:}=SGAcrossStAll;
    separateCellFGBl{c,:}=FGAcrossBlAll;
    separateCellFGSt{c,:}=FGAcrossStAll;
end
%% Collecting Stim and Sham data separately for later usage
if FRFlag==1
    stimDataFRBl(:,1)=compareCellFRBl{1,1}; shamDataFRBl(:,1)=compareCellFRBl{2,1};
    stimDataFRSt(:,1)=compareCellFRSt{1,1}; shamDataFRSt(:,1)=compareCellFRSt{2,1};
    stimDataFRTransient(:,1)=compareCellFRTransient{1,1}; shamDataFRTransient(:,1)=compareCellFRTransient{2,1};
    if normalizeFlag==1 % normalizing the firing rate data
        for iProt=1:size(protocol,2)
            normShamDataFrSt{iProt,1}=shamDataFRSt{iProt,1}./shamDataFRSt{1,1};
            normShamDataFrBl{iProt,1}=shamDataFRBl{iProt,1}./shamDataFRBl{1,1};
            normShamDataFrTransient{iProt,1}=shamDataFRTransient{iProt,1}./shamDataFRTransient{1,1};

            normStimDataFrSt{iProt,1}=stimDataFRSt{iProt,1}./stimDataFRSt{1,1};
            normStimDataFrBl{iProt,1}=stimDataFRBl{iProt,1}./stimDataFRBl{1,1};
            normStimDataFrTransient{iProt,1}=stimDataFRTransient{iProt,1}./stimDataFRTransient{1,1};
        end

        clear shamDataFRBl stimDataFRBl shamDataFRSt stimDataFRSt shamDataFRTransient stimDataFRTransient
        stimDataFRBl=normStimDataFrBl;
        shamDataFRBl=normShamDataFrBl;
        stimDataFRSt=normStimDataFrSt;
        shamDataFRSt=normShamDataFrSt;

        if transientFlag==1
            stimDataFRTransient=normStimDataFrTransient;
            shamDataFRTransient=normShamDataFrTransient;
            clear stimDataFRSt shamDataFRSt
            stimDataFRSt=stimDataFRTransient;
            shamDataFRSt=shamDataFRTransient;
        end
    end
end

stimDataSt(:,1)=separateCellSGSt{1,1}; shamDataSt(:,1)=separateCellSGSt{2,1};
stimDataBl(:,1)=separateCellSGBl{1,1}; shamDataBl(:,1)=separateCellSGBl{2,1};

stimDataSt(:,2)=separateCellFGSt{1,1}; shamDataSt(:,2)=separateCellFGSt{2,1};
stimDataBl(:,2)=separateCellFGBl{1,1}; shamDataBl(:,2)=separateCellFGBl{2,1};

%% Data assorting and taking care of unequal days of recording
if FRFlag==1
    if strcmp(session,'single')
        for protNum=1:size(protocol,2)
            if size(stimDataFRSt{protNum},1)~=size(shamDataFRSt{protNum},1)%If days of sham and stim does not match
                elecSize=size(cBlSG,2);
                dayStim=size(stimDataFRSt{protNum},1)./elecSize;
                daySham=size(shamDataFRSt{protNum},1)./elecSize;
                uFR{protNum}=(stimDataFRSt{protNum}-stimDataFRBl{protNum});
                vFR{protNum}=(shamDataFRSt{protNum}-shamDataFRBl{protNum});

                %%Averaging Stim days
                for d=1:dayStim
                    firstElecs(1,d)=(d*elecSize-elecSize+1);
                    U{1,d}=uFR{protNum}(firstElecs(1,d):firstElecs(1,d)+elecSize-1);
                end
                %%Averaging Sham days
                for d=1:daySham
                    firstElecs(1,d)=(d*elecSize-elecSize+1);
                    V{1,d}=vFR{protNum}(firstElecs(1,d):firstElecs(1,d)+elecSize-1);
                end
                dFR{protNum}=mean(cell2mat(U),2)-mean(cell2mat(V),2);
            elseif size(stimDataFRSt{protNum},1)==size(shamDataFRSt{protNum},1)
                dFR{protNum}=((stimDataFRSt{protNum})-(shamDataFRSt{protNum}));
            end

            avgdFR{protNum}=mean(dFR{protNum});
            deltaPlotDataFR{protNum}=avgdFR{protNum}-avgdFR{1};%Making the first point zero
            errorDeltaDataFR=dFR{protNum};
            errorDeltaGridFR{protNum}=std(errorDeltaDataFR)/sqrt(size(errorDeltaDataFR,1));
        end
    elseif strcmp(session,'dual')
        % We are doing bootstrapping here, as Stimulation day good spiking
        % elecs(n) =~ Sham day good spiking elecs(n)
        % get any index between 1 to 9 to chose from 20 units we have
        numOfStimElec=size(stimDataFRSt{1, 1},1);
        numOfShamElec=size(shamDataFRSt{1, 1},1);
        for iProt=1:size(protocol,2)
            for n=1:10000
                randomElecID = randperm(numOfStimElec,numOfShamElec);% Generate 9 random numbers between 1 and 20
                for i=1:numOfShamElec
                    deltaValue(i)=stimDataFRSt{iProt,1}(randomElecID(i))-shamDataFRSt{iProt,1}(i);% Here one of the 10,000 times, the subtraction is happening for each elec.
                end
                diffGrid(n,:)=deltaValue;
            end
            dFR{iProt}=mean(diffGrid,1)'; % Averaging across observations and getting the difference values across electrode, point to point
            avgdFR{iProt}=mean(dFR{iProt},1);% Averaging across elecrodes, getting a single value for plotting
            deltaPlotDataFR{iProt}=avgdFR{iProt}-avgdFR{1};
            errorDeltaDataFR=dFR{iProt};
            errorDeltaGridFR{iProt}=std(errorDeltaDataFR)/sqrt(size(errorDeltaDataFR,1));
            clear diffGrid deltaValue
        end
    end

    %% plotting the FR plot
    hold on
    point=1:size(protocol,2);
    hold(plotHandle(1,1),'on')
    e=errorbar(subplot(plotHandle(1,1)),point,[deltaPlotDataFR{:}],[errorDeltaGridFR{:}],'color','#ED4672','LineWidth',1.5);
    for f=1:size(protocol,2)
        w=[deltaPlotDataFR{:}];
        scatter(point(1,f),w(1,f),30,colCode{1, f}{1, 1},"filled")
    end
    Xax=gca().XAxis;
    Yax=gca().YAxis;
    set(Xax,'TickDirection','out');
    set(Xax,'TickLength',[0.02 0.025]);
    set(Yax,'TickDirection','out');
    set(Yax,'TickLength',[0.02 0.025]);
    % set(gca,'FontWeight','bold');
    set(Xax,'FontSize',10);
    set(Yax,'FontSize',10);

    if strcmp(session,'single')
        line([0 length(point)+1], [0.99 0.99],'color','k')
        line([length(point)+1 length(point)+1], [-1 1],'color','k')
    elseif strcmp(session,'dual')
        line([0 length(point)+2], [0.99 0.99],'color','k')
        line([length(point)+2 length(point)+2], [-1 1],'color','k')
    end

    e.Marker='.';
    e.CapSize=3;
    yline(deltaPlotDataFR{1,1},'--')
    text(0.4,0.88,append('n=',num2str(size(dFR{1,1},1))),'color','k','FontSize',9,'FontWeight','bold');
    text(4.225,0.88,append( 'p<0.05(*)'),'color','k','FontSize',9,'FontWeight','bold');
    set(plotHandle(1,:),'XTick',[]);

    %% The Stat part
    for protN=2:size(protocol,2)
        x1=dFR{1,1};% Post-Pre
        y1=dFR{protN};
        [p1,h1]=ranksum(y1,x1); %This checks if gamma power in a protocol is significantly different than the pre-stim/sham protocol
        SGdeltaStat(1,protN)=double(h1);
        SGdeltaStat(2,protN)=p1;
    end

    %% Plotting the significance stars for FR plot
    z=[deltaPlotDataFR{:}];
    asterixData=(z);
    xCentres=point;
    hold on
    stat1=SGdeltaStat(1,:);
    stat2=SGdeltaStat(2,:);

    delete(findall(subplot(plotHandle(1,1)), 'Tag', 'SigStar')); %deleting the old significance stars
    for d=1:size(protocol,2)
        if stat1(1,d)==1
            hold on
            subplot(plotHandle(1,1))
            ypoint=(asterixData(:,d));
            if ypoint<0
                if stat2(1,d)<0.0005
                    text(xCentres(:,d)-0.3,ypoint-(cell2mat(errorDeltaGridFR(:,d))+0.15),'\ast\ast\ast','fontWeight','bold','tag','SigStar',HandleVisibility='off'); hold on
                elseif stat2(1,d)<0.005
                    text(xCentres(:,d)-0.2,ypoint-(cell2mat(errorDeltaGridFR(:,d))+0.15),'\ast\ast','fontWeight','bold','tag','SigStar',HandleVisibility='off'); hold on
                elseif stat2(1,d)<0.05
                    text(xCentres(:,d)-0.1,ypoint-(cell2mat(errorDeltaGridFR(:,d))+0.15),'\ast','fontWeight','bold','tag','SigStar',HandleVisibility='off'); hold on
                end
            elseif ypoint>0
                % text(Xcentres(:,d)-0.2,ypoint+0.25,'\ast','fontWeight','bold',HandleVisibility='off'); hold on
                if stat2(1,d)<0.0005
                    text(xCentres(:,d)-0.3,ypoint+(cell2mat(errorDeltaGridFR(:,d))+0.2),'\ast\ast\ast','fontWeight','bold','tag','SigStar',HandleVisibility='off'); hold on
                elseif stat2(1,d)<0.005
                    text(xCentres(:,d)-0.2,ypoint+(cell2mat(errorDeltaGridFR(:,d))+0.2),'\ast\ast','fontWeight','bold','tag','SigStar',HandleVisibility='off'); hold on
                elseif stat2(1,d)<0.05
                    text(xCentres(:,d)-0.1,ypoint+(cell2mat(errorDeltaGridFR(:,d))+0.2),'\ast','fontWeight','bold','tag','SigStar',HandleVisibility='off'); hold on
                end
            end
        end
    end
end

if PSDFlag==1
    %% Data assorting and taking care of unequal days of recording
    for band=bandFlag
        for protNum=1:size(protocol,2)
            if size(stimDataSt{protNum,band},1)~=size(shamDataSt{protNum,band},1)%If days of sham and stim does not match
                elecSize=size(cSG,2);
                dayStim=size(stimDataSt{protNum,band},1)./elecSize;
                daySham=size(shamDataSt{protNum,band},1)./elecSize;
                u{protNum,band}=(stimDataSt{protNum,band}-stimDataBl{protNum,band});
                v{protNum,band}=(shamDataSt{protNum,band}-shamDataBl{protNum,band});

                %%Averaging Stim days
                for d=1:dayStim
                    firstElecs(1,d)=(d*elecSize-elecSize+1);
                    U{1,d}=u{protNum,band}(firstElecs(1,d):firstElecs(1,d)+elecSize-1);
                end
                %%Averaging Sham days
                for d=1:daySham
                    firstElecs(1,d)=(d*elecSize-elecSize+1);
                    V{1,d}=v{protNum,band}(firstElecs(1,d):firstElecs(1,d)+elecSize-1);
                end

                dPower{protNum,band}=mean(cell2mat(U),2)-mean(cell2mat(V),2);
            elseif size(stimDataSt{protNum,band},1)==size(shamDataSt{protNum,band},1)
                dPower{protNum,band}=((stimDataSt{protNum,band}-stimDataBl{protNum,band})-(shamDataSt{protNum,band}-shamDataBl{protNum,band}));
            end

            AvgdPower{protNum,band}=mean(dPower{protNum,band});
            deltaPlotDataF{protNum,band}=AvgdPower{protNum,band}-AvgdPower{1,band};%Making the first point zero
            errordeltadataF=10*dPower{protNum,band}; %Log power converted to dB
            errordeltagridF{protNum,band}=std(errordeltadataF)/sqrt(size(errordeltadataF,1));
        end
    end

    %% Plotting band plot
    hold on
    point=1:size(protocol,2);

    for band=bandFlag %%% For slow and fast gamma band
        hold(plotHandle(1,1),'on')
        axes(subplot(plotHandle(1,1)));
        e=errorbar(subplot(plotHandle(1,1)),point,10*[deltaPlotDataF{:,band}],[errordeltagridF{:,band}],'color','#ED4672','LineWidth',1.5);
        for f=1:size(protocol,2)
            w=10*[deltaPlotDataF{:,band}];
            scatter(point(1,f),w(1,f),30,colCode{1, f}{1, 1},"filled")
        end

        Xax=gca().XAxis;
        Yax=gca().YAxis;
        set(Xax,'TickDirection','out');
        set(Xax,'TickLength',[0.02 0.025]);
        set(Yax,'TickDirection','out');
        set(Yax,'TickLength',[0.02 0.025]);
        % set(gca,'FontWeight','bold');
        set(Xax,'FontSize',10);
        set(Yax,'FontSize',10);
        if strcmp(session,'single')
            line([0 length(point)+1], [1.99 1.99],'color','k')
            line([length(point)+1 length(point)+1], [-2 2],'color','k')
        elseif strcmp(session,'dual')
            line([0 length(point)+2], [1.99 1.99],'color','k')
            line([length(point)+2 length(point)+2], [-2 2],'color','k')
        end
        e.Marker='.';
        e.CapSize=3;
        yline(deltaPlotDataF{1,band},'--')
        text(4.225,1.8,append( 'p<0.05(*)'),'color','k','FontSize',9,'FontWeight','bold');
        set(plotHandle(1,:),'XTick',[]);
        text(0.4,1.8,append('n=',num2str(size(dPower{1,band},1))),'color','k','FontSize',9,'FontWeight','bold');
    end

    %% The stat part
    for band=bandFlag
        for protN=2:size(protocol,2)
            x1=dPower{1,band};% Post-Pre
            y1=dPower{protN,band};
            [p1,h1]=ranksum(y1,x1); %This checks if gamma power in a protocol is significantly different than the pre-stim/sham protocol

            if band==1
                SGdeltaStat(1,protN)=double(h1);
                SGdeltaStat(2,protN)=p1;
                clear x1 y1 h1 p1
            elseif band==2
                FGdeltaStat(1,protN)=double(h1);
                FGdeltaStat(2,protN)=p1;
                clear x1 y1 h1 p1
            end
        end
    end

    %% Plotting the significance star
    for freqstat=bandFlag %Band specific identity
        m=10*[deltaPlotDataF{:,freqstat}];
        asterixData=(m);
        xCentres=point;
        hold on
        if freqstat==1
            stat1=SGdeltaStat(1,:);
            stat2=SGdeltaStat(2,:);
        elseif freqstat==2
            stat1=FGdeltaStat(1,:);
            stat2=FGdeltaStat(2,:);
        end
        delete(findall(subplot(plotHandle(1,1)), 'Tag', 'SigStar'));%deleting older sigstars on plot
        for d=1:size(protocol,2)
            if stat1(1,d)==1
                hold on
                subplot(plotHandle(1,1))
                ypoint=(asterixData(:,d));
                if ypoint<0
                    if stat2(1,d)<0.0005
                        text(xCentres(:,d)-0.3,ypoint-(cell2mat(errordeltagridF(d,freqstat))+0.3),'\ast\ast\ast','fontWeight','bold','tag','SigStar',HandleVisibility='off'); hold on
                    elseif stat2(1,d)<0.005
                        text(xCentres(:,d)-0.2,ypoint-(cell2mat(errordeltagridF(d,freqstat))+0.3),'\ast\ast','fontWeight','bold','tag','SigStar',HandleVisibility='off'); hold on
                    elseif stat2(1,d)<0.05
                        text(xCentres(:,d)-0.1,ypoint-(cell2mat(errordeltagridF(d,freqstat))+0.3),'\ast','fontWeight','bold','tag','SigStar',HandleVisibility='off'); hold on
                    end
                elseif ypoint>0
                    if stat2(1,d)<0.0005
                        text(xCentres(:,d)-0.35,ypoint+(cell2mat(errordeltagridF(d,freqstat))+0.4),'\ast\ast\ast','fontWeight','bold','tag','SigStar',HandleVisibility='off'); hold on
                    elseif stat2(1,d)<0.005
                        text(xCentres(:,d)-0.25,ypoint+(cell2mat(errordeltagridF(d,freqstat))+0.4),'\ast\ast','fontWeight','bold','tag','SigStar',HandleVisibility='off'); hold on
                    elseif stat2(1,d)<0.05
                        text(xCentres(:,d)-0.15,ypoint+(cell2mat(errordeltagridF(d,freqstat))+0.4),'\ast','fontWeight','bold','tag','SigStar',HandleVisibility='off'); hold on
                    end
                end
            end
        end
    end
end
end

