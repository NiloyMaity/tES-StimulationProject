%% Please keep the Stimulation Specific Protocol List in MATLAB path.
% Area=V1 or V4, so AreaFlag= 1 or 2 corresponds to V1 & V4
% StimulationType={'tDCS','tACS'}
% condition={'Stim','Sham'};
% Polarity={'Cathodal','Anodal' or 'SG','FG','Alpha'};
% Session={'single','dual','dual60'};

%% This code plots Line plots: delta Firing rate plot and delta band power plot

function plotSpikePower(plothandle,folderin,session,normalizeFlag,transientFlag,FRFlag,PSDFlag,BandFlag,dates,protocols,SFVals,ConVals,OriVals,ColCode,titleString,StimblockID)
baselineS  = [-0.5 0];          % Baseline period for computing change
stimPeriod = [0.25 0.75];
transientPeriod= [0 0.1];
binWidthMS=10;

for ses=1:length(StimblockID)
    ColCode{1,StimblockID(ses)}=[];
    titleString{1,StimblockID(ses)}=[];
end

Col=ColCode(~cellfun('isempty',ColCode));
ColCode=reshape(Col,[1 length(ColCode)-length(StimblockID)]);
% Collecting Data
for c=1:2
    for q=1:length(dates{1,c})
        for r=1:length(StimblockID)
            protocols{1,c}{q,StimblockID(r)}=[];
        end
    end
    % clearing the stimulation block from protocol
    y=protocols{1,c}(~cellfun('isempty',protocols{1,c}));
    protocol=reshape(y,[length(dates{1,c}) length(protocols{1,c})-length(StimblockID)]);

    % data collection loop
    for day=1:length(dates{1,c})
        for j=1:length(SFVals)
            for k=1:length(OriVals)
                for l=1:length(ConVals)
                    DataBand = dir(fullfile(folderin{1,c},append(num2str(SFVals(j)),'SF','_',num2str(OriVals(k)),'Ori','_',num2str(ConVals(l)),'Con','_',dates{1,c}{1,day},'*BandData.mat')));
                    clear BandValues
                    BandValues = load(fullfile(folderin{1,c},DataBand.name));
                    cSG=BandValues.CollectSG;
                    cFG=BandValues.CollectFG;
                    cBlSG=BandValues.CollectBlSG;
                    cStSG=BandValues.CollectStSG;
                    cBlFG=BandValues.CollectBlFG;
                    cStFG=BandValues.CollectStFG;
                    for r=1:length(StimblockID)
                        cSG(StimblockID(r),:)=0;
                        cFG(StimblockID(r),:)=0;
                        cBlSG(StimblockID(r),:)=0;
                        cStSG(StimblockID(r),:)=0;
                        cBlFG(StimblockID(r),:)=0;
                        cStFG(StimblockID(r),:)=0;
                    end

                    cBlSG=reshape(nonzeros(cBlSG),[size(protocol,2) size(cBlSG,2)]);
                    cStSG=reshape(nonzeros(cStSG),[size(protocol,2) size(cStSG,2)]);
                    cBlFG=reshape(nonzeros(cBlFG),[size(protocol,2) size(cBlFG,2)]);
                    cStFG=reshape(nonzeros(cStFG),[size(protocol,2) size(cStFG,2)]);

                    for iProt =1:size(protocol,2)
                        DataSpike = dir(fullfile(folderin{1,c},append(num2str(SFVals(j)),'SF','_',num2str(OriVals(k)),'Ori','_',num2str(ConVals(l)),'Con','_',protocol{day,iProt},'_',dates{1,c}{1,day},'*PSTH.mat')));
                        % Loading the spike file
                        if ~isempty(DataSpike)
                            s1 = load(fullfile(folderin{1,c},DataSpike.name)); % DataStructure(field).name
                            RawPSTH{iProt,j,k,l} = s1.PSTHGrid;
                            SPtList=s1.xs; %Spike timevals
                        else
                            RawPSTH{iProt,j,k,l} =[];
                        end

                        if ~isempty(DataSpike)
                            blPos = find(SPtList>=baselineS(1),1)+ (1:(diff(baselineS))/(binWidthMS/1000));
                            stPos = find(SPtList>=stimPeriod(1),1)+ (1:(diff(stimPeriod))/(binWidthMS/1000));
                            transientPos=find(SPtList>=transientPeriod(1),1)+ (1:(diff(transientPeriod))/(binWidthMS/1000));

                            RawFRBl{iProt,j,k,l}= mean(s1.PSTHGrid(:,blPos),2);
                            RawFRSt{iProt,j,k,l}= mean(s1.PSTHGrid(:,stPos),2);
                            RawFRTransient{iProt,j,k,l}= mean(s1.PSTHGrid(:,transientPos),2);
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
            if ~isempty(DataSpike)
                RawFRAcrossBl{day,iProt,:}=mean(cell2mat(RawFRBl(iProt,:)),2);
                RawFRAcrossSt{day,iProt,:}=mean(cell2mat(RawFRSt(iProt,:)),2);
                RawFRAcrossTransient{day,iProt,:}=mean(cell2mat(RawFRTransient(iProt,:)),2);
            end

            % Averaging across conditions
            BandSGAcrossBl{day,iProt,:}=mean(cell2mat(SGAcrossBl(iProt,:)),2);
            BandFGAcrossBl{day,iProt,:}=mean(cell2mat(FGAcrossBl(iProt,:)),2);
            BandSGAcrossSt{day,iProt,:}=mean(cell2mat(SGAcrossSt(iProt,:)),2);
            BandFGAcrossSt{day,iProt,:}=mean(cell2mat(FGAcrossSt(iProt,:)),2);
        end
    end

    %% Here pulling of all data from required dates has been completed.
    clear SGAcrossAll FGAcrossAll SGAcrossBlAll SGAcrossStAll FGAcrossBlAll FGAcrossStAll FRAcrossBlAll FRAcrossStAll FRAcrossTransientAll
    for b=1:size(protocol,2)
        if FRFlag==1
            FRAcrossBlAll{b,:}=cat(1,RawFRAcrossBl{:,b});
            FRAcrossStAll{b,:}=cat(1,RawFRAcrossSt{:,b});
            FRAcrossTransientAll{b,:}=cat(1,RawFRAcrossTransient{:,b});
        end
        SGAcrossBlAll{b,:}=cat(1,BandSGAcrossBl{:,b});
        FGAcrossBlAll{b,:}=cat(1,BandFGAcrossBl{:,b});
        SGAcrossStAll{b,:}=cat(1,BandSGAcrossSt{:,b});
        FGAcrossStAll{b,:}=cat(1,BandFGAcrossSt{:,b});
    end
    clear SGAcross FGAcross SGAcrossBl FGAcrossBl SGAcrossSt FGAcrossSt RawFRAcrossBl RawFRAcrossSt RawFRAcrossTransient BandSGAcross BandFGAcross BandSGAcrossBl BandFGAcrossBl BandSGAcrossSt BandFGAcrossSt
    % And This brings all electrode in one page

    %% Putting all electrodes' band specific data in a cell, so that we can later do subtraction between conditions (Stim-Sham)
    if FRFlag==1
        CompareCellFRBl{c,:}=FRAcrossBlAll;
        CompareCellFRSt{c,:}=FRAcrossStAll;
        CompareCellFRTransient{c,:}=FRAcrossTransientAll;
    end
    SeparateCellSGBl{c,:}=SGAcrossBlAll;
    SeparateCellSGSt{c,:}=SGAcrossStAll;
    SeparateCellFGBl{c,:}=FGAcrossBlAll;
    SeparateCellFGSt{c,:}=FGAcrossStAll;
end
%% Collecting Stim and Sham data separately for later usage
if FRFlag==1
    StimDataFRBl(:,1)=CompareCellFRBl{1,1}; ShamDataFRBl(:,1)=CompareCellFRBl{2,1};
    StimDataFRSt(:,1)=CompareCellFRSt{1,1}; ShamDataFRSt(:,1)=CompareCellFRSt{2,1};
    StimDataFRTransient(:,1)=CompareCellFRTransient{1,1}; ShamDataFRTransient(:,1)=CompareCellFRTransient{2,1};
    if normalizeFlag==1 % normalizing the firing rate data
        for iProt=1:size(protocol,2)
            normShamDataFrSt{iProt,1}=ShamDataFRSt{iProt,1}./ShamDataFRSt{1,1};
            normShamDataFrBl{iProt,1}=ShamDataFRBl{iProt,1}./ShamDataFRBl{1,1};
            normShamDataFrTransient{iProt,1}=ShamDataFRTransient{iProt,1}./ShamDataFRTransient{1,1};

            normStimDataFrSt{iProt,1}=StimDataFRSt{iProt,1}./StimDataFRSt{1,1};
            normStimDataFrBl{iProt,1}=StimDataFRBl{iProt,1}./StimDataFRBl{1,1};
            normStimDataFrTransient{iProt,1}=StimDataFRTransient{iProt,1}./StimDataFRTransient{1,1};
        end

        clear ShamDataFRBl StimDataFRBl ShamDataFRSt StimDataFRSt ShamDataFRTransient StimDataFRTransient
        StimDataFRBl=normStimDataFrBl;
        ShamDataFRBl=normShamDataFrBl;
        StimDataFRSt=normStimDataFrSt;
        ShamDataFRSt=normShamDataFrSt;

        if transientFlag==1
            StimDataFRTransient=normStimDataFrTransient;
            ShamDataFRTransient=normShamDataFrTransient;
            clear StimDataFRSt ShamDataFRSt
            StimDataFRSt=StimDataFRTransient;
            ShamDataFRSt=ShamDataFRTransient;
        end
    end
end

StimDataSt(:,1)=SeparateCellSGSt{1,1}; ShamDataSt(:,1)=SeparateCellSGSt{2,1};
StimDataBl(:,1)=SeparateCellSGBl{1,1}; ShamDataBl(:,1)=SeparateCellSGBl{2,1};

StimDataSt(:,2)=SeparateCellFGSt{1,1}; ShamDataSt(:,2)=SeparateCellFGSt{2,1};
StimDataBl(:,2)=SeparateCellFGBl{1,1}; ShamDataBl(:,2)=SeparateCellFGBl{2,1};

%% Data assorting and taking care of unequal days of recording
if FRFlag==1
    if strcmp(session,'single')
        for Protnum=1:size(protocol,2)
            if size(StimDataFRSt{Protnum},1)~=size(ShamDataFRSt{Protnum},1)%If days of sham and stim does not match
                elecsize=size(cBlSG,2);
                dayStim=size(StimDataFRSt{Protnum},1)./elecsize;
                daySham=size(ShamDataFRSt{Protnum},1)./elecsize;
                uFR{Protnum}=(StimDataFRSt{Protnum}-StimDataFRBl{Protnum});
                vFR{Protnum}=(ShamDataFRSt{Protnum}-ShamDataFRBl{Protnum});

                %%Averaging Stim days
                for d=1:dayStim
                    FirstElecs(1,d)=(d*elecsize-elecsize+1);
                    U{1,d}=uFR{Protnum}(FirstElecs(1,d):FirstElecs(1,d)+elecsize-1);
                end
                %%Averaging Sham days
                for d=1:daySham
                    FirstElecs(1,d)=(d*elecsize-elecsize+1);
                    V{1,d}=vFR{Protnum}(FirstElecs(1,d):FirstElecs(1,d)+elecsize-1);
                end
                dFR{Protnum}=mean(cell2mat(U),2)-mean(cell2mat(V),2);
            elseif size(StimDataFRSt{Protnum},1)==size(ShamDataFRSt{Protnum},1)
                dFR{Protnum}=((StimDataFRSt{Protnum})-(ShamDataFRSt{Protnum}));
            end

            AvgdFR{Protnum}=mean(dFR{Protnum});
            deltaPlotDataFR{Protnum}=AvgdFR{Protnum}-AvgdFR{1};%Making the first point zero
            errordeltadataFR=dFR{Protnum};
            errordeltagridFR{Protnum}=std(errordeltadataFR)/sqrt(size(errordeltadataFR,1));
        end
    elseif strcmp(session,'dual')
        % We are doing bootstrapping here, as Stimulation day good spiking
        % elecs(n) =~ Sham day good spiking elecs(n)
        % get any index between 1 to 9 to chose from 20 units we have
        numofStimElec=size(StimDataFRSt{1, 1},1);
        numofShamElec=size(ShamDataFRSt{1, 1},1);
        for iProt=1:size(protocol,2)
            for n=1:10000
                random_ElecID = randperm(numofStimElec,numofShamElec);% Generate 9 random numbers between 1 and 20
                for i=1:numofShamElec
                    deltaValue(i)=StimDataFRSt{iProt,1}(random_ElecID(i))-ShamDataFRSt{iProt,1}(i);% Here one of the 10,000 times, the subtraction is happening for each elec.
                end
                diffGrid(n,:)=deltaValue;
            end
            dFR{iProt}=mean(diffGrid,1)'; % Averaging across observations and getting the difference values across electrode, point to point
            AvgdFR{iProt}=mean(dFR{iProt},1);% Averaging across elecrodes, getting a single value for plotting
            deltaPlotDataFR{iProt}=AvgdFR{iProt}-AvgdFR{1};
            errordeltadataFR=dFR{iProt};
            errordeltagridFR{iProt}=std(errordeltadataFR)/sqrt(size(errordeltadataFR,1));
            clear diffGrid deltaValue
        end
    end

    %% plotting the FR plot
    hold on
    point=1:size(protocol,2);
    hold(plothandle(1,1),'on')
    e=errorbar(subplot(plothandle(1,1)),point,[deltaPlotDataFR{:}],[errordeltagridFR{:}],'color','#ED4672','LineWidth',1.5);
    for f=1:size(protocol,2)
        w=[deltaPlotDataFR{:}];
        scatter(point(1,f),w(1,f),30,ColCode{1, f}{1, 1},"filled")
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
    set(plothandle(1,:),'XTick',[]);

    %% The Stat part
    for ProtN=2:size(protocol,2)
        x1=dFR{1,1};% Post-Pre
        y1=dFR{ProtN};
        [p1,h1]=ranksum(y1,x1); %This checks if gamma power in a protocol is significantly different than the pre-stim/sham protocol
        SGdeltaStat(1,ProtN)=double(h1);
        SGdeltaStat(2,ProtN)=p1;
    end

    %% Plotting the significance stars for FR plot
    z=[deltaPlotDataFR{:}];
    Asterixdata=(z);
    Xcentres=point;
    hold on
    stat1=SGdeltaStat(1,:);
    stat2=SGdeltaStat(2,:);

    delete(findall(subplot(plothandle(1,1)), 'Tag', 'SigStar')); %deleting the old significance stars
    for d=1:size(protocol,2)
        if stat1(1,d)==1
            hold on
            subplot(plothandle(1,1))
            ypoint=(Asterixdata(:,d));
            if ypoint<0
                if stat2(1,d)<0.0005
                    text(Xcentres(:,d)-0.3,ypoint-(cell2mat(errordeltagridFR(:,d))+0.15),'\ast\ast\ast','fontWeight','bold','tag','SigStar',HandleVisibility='off'); hold on
                elseif stat2(1,d)<0.005
                    text(Xcentres(:,d)-0.2,ypoint-(cell2mat(errordeltagridFR(:,d))+0.15),'\ast\ast','fontWeight','bold','tag','SigStar',HandleVisibility='off'); hold on
                elseif stat2(1,d)<0.05
                    text(Xcentres(:,d)-0.1,ypoint-(cell2mat(errordeltagridFR(:,d))+0.15),'\ast','fontWeight','bold','tag','SigStar',HandleVisibility='off'); hold on
                end
            elseif ypoint>0
                % text(Xcentres(:,d)-0.2,ypoint+0.25,'\ast','fontWeight','bold',HandleVisibility='off'); hold on
                if stat2(1,d)<0.0005
                    text(Xcentres(:,d)-0.3,ypoint+(cell2mat(errordeltagridFR(:,d))+0.2),'\ast\ast\ast','fontWeight','bold','tag','SigStar',HandleVisibility='off'); hold on
                elseif stat2(1,d)<0.005
                    text(Xcentres(:,d)-0.2,ypoint+(cell2mat(errordeltagridFR(:,d))+0.2),'\ast\ast','fontWeight','bold','tag','SigStar',HandleVisibility='off'); hold on
                elseif stat2(1,d)<0.05
                    text(Xcentres(:,d)-0.1,ypoint+(cell2mat(errordeltagridFR(:,d))+0.2),'\ast','fontWeight','bold','tag','SigStar',HandleVisibility='off'); hold on
                end
            end
        end
    end
end

if PSDFlag==1
    %% Data assorting and taking care of unequal days of recording
    for Band=BandFlag
        for Protnum=1:size(protocol,2)
            if size(StimDataSt{Protnum,Band},1)~=size(ShamDataSt{Protnum,Band},1)%If days of sham and stim does not match
                elecsize=size(cSG,2);
                dayStim=size(StimDataSt{Protnum,Band},1)./elecsize;
                daySham=size(ShamDataSt{Protnum,Band},1)./elecsize;
                u{Protnum,Band}=(StimDataSt{Protnum,Band}-StimDataBl{Protnum,Band});
                v{Protnum,Band}=(ShamDataSt{Protnum,Band}-ShamDataBl{Protnum,Band});

                %%Averaging Stim days
                for d=1:dayStim
                    FirstElecs(1,d)=(d*elecsize-elecsize+1);
                    U{1,d}=u{Protnum,Band}(FirstElecs(1,d):FirstElecs(1,d)+elecsize-1);
                end
                %%Averaging Sham days
                for d=1:daySham
                    FirstElecs(1,d)=(d*elecsize-elecsize+1);
                    V{1,d}=v{Protnum,Band}(FirstElecs(1,d):FirstElecs(1,d)+elecsize-1);
                end

                dPower{Protnum,Band}=mean(cell2mat(U),2)-mean(cell2mat(V),2);
            elseif size(StimDataSt{Protnum,Band},1)==size(ShamDataSt{Protnum,Band},1)
                dPower{Protnum,Band}=((StimDataSt{Protnum,Band}-StimDataBl{Protnum,Band})-(ShamDataSt{Protnum,Band}-ShamDataBl{Protnum,Band}));
            end

            AvgdPower{Protnum,Band}=mean(dPower{Protnum,Band});
            deltaPlotDataF{Protnum,Band}=AvgdPower{Protnum,Band}-AvgdPower{1,Band};%Making the first point zero
            errordeltadataF=10*dPower{Protnum,Band}; %Log power converted to dB
            errordeltagridF{Protnum,Band}=std(errordeltadataF)/sqrt(size(errordeltadataF,1));
        end
    end

    %% Plotting band plot
    hold on
    point=1:size(protocol,2);

    for Band=BandFlag %%% For slow and fast gamma band
        hold(plothandle(1,1),'on')
        axes(subplot(plothandle(1,1)));
        e=errorbar(subplot(plothandle(1,1)),point,10*[deltaPlotDataF{:,Band}],[errordeltagridF{:,Band}],'color','#ED4672','LineWidth',1.5);
        for f=1:size(protocol,2)
            w=10*[deltaPlotDataF{:,Band}];
            scatter(point(1,f),w(1,f),30,ColCode{1, f}{1, 1},"filled")
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
        yline(deltaPlotDataF{1,Band},'--')
        text(4.225,1.8,append( 'p<0.05(*)'),'color','k','FontSize',9,'FontWeight','bold');
        set(plothandle(1,:),'XTick',[]);
        text(0.4,1.8,append('n=',num2str(size(dPower{1,Band},1))),'color','k','FontSize',9,'FontWeight','bold');
    end

    %% The stat part
    for Band=BandFlag
        for ProtN=2:size(protocol,2)
            x1=dPower{1,Band};% Post-Pre
            y1=dPower{ProtN,Band};
            [p1,h1]=ranksum(y1,x1); %This checks if gamma power in a protocol is significantly different than the pre-stim/sham protocol

            if Band==1
                SGdeltaStat(1,ProtN)=double(h1);
                SGdeltaStat(2,ProtN)=p1;
                clear x1 y1 h1 p1
            elseif Band==2
                FGdeltaStat(1,ProtN)=double(h1);
                FGdeltaStat(2,ProtN)=p1;
                clear x1 y1 h1 p1
            end
        end
    end

    %% Plotting the significance star
    for freqstat=BandFlag %Band specific identity
        m=10*[deltaPlotDataF{:,freqstat}];
        Asterixdata=(m);
        Xcentres=point;
        hold on
        if freqstat==1
            stat1=SGdeltaStat(1,:);
            stat2=SGdeltaStat(2,:);
        elseif freqstat==2
            stat1=FGdeltaStat(1,:);
            stat2=FGdeltaStat(2,:);
        end
        delete(findall(subplot(plothandle(1,1)), 'Tag', 'SigStar'));%deleting older sigstars on plot
        for d=1:size(protocol,2)
            if stat1(1,d)==1
                hold on
                subplot(plothandle(1,1))
                ypoint=(Asterixdata(:,d));
                if ypoint<0
                    if stat2(1,d)<0.0005
                        text(Xcentres(:,d)-0.3,ypoint-(cell2mat(errordeltagridF(d,freqstat))+0.3),'\ast\ast\ast','fontWeight','bold','tag','SigStar',HandleVisibility='off'); hold on
                    elseif stat2(1,d)<0.005
                        text(Xcentres(:,d)-0.2,ypoint-(cell2mat(errordeltagridF(d,freqstat))+0.3),'\ast\ast','fontWeight','bold','tag','SigStar',HandleVisibility='off'); hold on
                    elseif stat2(1,d)<0.05
                        text(Xcentres(:,d)-0.1,ypoint-(cell2mat(errordeltagridF(d,freqstat))+0.3),'\ast','fontWeight','bold','tag','SigStar',HandleVisibility='off'); hold on
                    end
                elseif ypoint>0
                    if stat2(1,d)<0.0005
                        text(Xcentres(:,d)-0.35,ypoint+(cell2mat(errordeltagridF(d,freqstat))+0.4),'\ast\ast\ast','fontWeight','bold','tag','SigStar',HandleVisibility='off'); hold on
                    elseif stat2(1,d)<0.005
                        text(Xcentres(:,d)-0.25,ypoint+(cell2mat(errordeltagridF(d,freqstat))+0.4),'\ast\ast','fontWeight','bold','tag','SigStar',HandleVisibility='off'); hold on
                    elseif stat2(1,d)<0.05
                        text(Xcentres(:,d)-0.15,ypoint+(cell2mat(errordeltagridF(d,freqstat))+0.4),'\ast','fontWeight','bold','tag','SigStar',HandleVisibility='off'); hold on
                    end
                end
            end
        end
    end
end
end

