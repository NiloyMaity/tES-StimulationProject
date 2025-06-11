%% please keep the Stimulation Specific Protocol List in MATLAB path.
%area=V1 or V4, so AreaFlag= 1 or 2 corresponds to V1 & V4
%stimulationType={'tDCS','tACS'}
%condition={'Stim','Sham'};
%polarity={'Cathodal','Anodal' or 'SG','FG','Alpha'};
%session={'single','dual','dual60'};
%sessionID={0,1,2}

%option of selecting LFP and spike electrode is of importance, pending.

function displayStimulationAll(monkeyName,folderSource,stimulationString,condition,polarityString,session,bandString,sf,con,ori,badTrialNameStr)

if ~exist('folderSource','var');  folderSource='Z:\Projects\Niloy_tES-StimulationProject';        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display main options
% fonts
fontSizeSmall = 10; fontSizeMedium = 12; fontSizeLarge = 16;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make Panels
panelHeight = 0.15; panelStartHeight = 0.85;
stimulationPanelWidth = 0.11; stimulationStartPos = 0.025;
dynamicPanelWidth = 0.15; dynamicStartPos = 0.15;
rangesPanelWidth =0.17; rangesStartPos =0.35;
timingPanelWidth = 0.15; timingStartPos = 0.625;
plotOptionsPanelWidth = 0.15; plotOptionsStartPos = 0.79;
actualplotPanelWidth=0.03; actualplotStartPos=0.96;
backgroundColor = 'w';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Stimulation panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stimulationHeight = 0.07; stimulationGap=0.015; stimulationTextWidth = 0.6;
hstimulationPanel = uipanel('Title','StimType','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[stimulationStartPos panelStartHeight stimulationPanelWidth panelHeight]);

% Stimulation Type
uicontrol('Parent',hstimulationPanel,'Unit','Normalized', ...
    'Position',[0 1-3*(stimulationHeight+stimulationGap) stimulationTextWidth stimulationHeight*3], ...
    'Style','text','String','Stimulation Type','FontSize',fontSizeSmall);
hStimulationType = uicontrol('Parent',hstimulationPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [stimulationTextWidth 1-2.5*(stimulationHeight+stimulationGap) 1-stimulationTextWidth stimulationHeight*3], ...
    'Style','popup','String',stimulationString,'FontSize',fontSizeSmall-1);

% Polarity
uicontrol('Parent',hstimulationPanel,'Unit','Normalized', ...
    'Position',[0 1-6*(stimulationHeight+stimulationGap) stimulationTextWidth stimulationHeight*3], ...
    'Style','text','String','Polarity','FontSize',fontSizeSmall);
hPolarity = uicontrol('Parent',hstimulationPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [stimulationTextWidth 1-5.5*(stimulationHeight+stimulationGap) 1-stimulationTextWidth stimulationHeight*3], ...
    'Style','popup','String',polarityString,'FontSize',fontSizeSmall-1);

% Choice of Band
uicontrol('Parent',hstimulationPanel,'Unit','Normalized', ...
    'Position',[0 1-9*(stimulationHeight+stimulationGap) stimulationTextWidth stimulationHeight*3], ...
    'Style','text','String','Band','FontSize',fontSizeSmall);
hBand = uicontrol('Parent',hstimulationPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [stimulationTextWidth 1-8.5*(stimulationHeight+stimulationGap) 1-stimulationTextWidth stimulationHeight*3], ...
    'Style','popup','String',bandString,'FontSize',fontSizeSmall-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameter panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dynamicHeight = 0.06; dynamicGap=0.015; dynamicTextWidth = 0.6;
hDynamicPanel = uipanel('Title','Parameters','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[dynamicStartPos panelStartHeight dynamicPanelWidth panelHeight]);

% Spatial Frequency
spatialFreqString = getStringFromValues(sf,1);
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-3*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight*3], ...
    'Style','text','String','Spatial Freq (CPD)','FontSize',fontSizeSmall);
hSpatialFreq = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-2.5*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight*3], ...
    'Style','popup','String',spatialFreqString,'FontSize',fontSizeSmall-1);

% Orientation
orientationString = getStringFromValues(ori,1);
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-6*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight*3], ...
    'Style','text','String','Orientation (Deg)','FontSize',fontSizeSmall);
hOrientation = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-5.5*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight*3], ...
    'Style','popup','String',orientationString,'FontSize',fontSizeSmall-1);

% Contrast
contrastString = getStringFromValues(con,1);
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-9*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight*3], ...
    'Style','text','String','Contrast (%)','FontSize',fontSizeSmall);
hContrast = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-8.5*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight*3], ...
    'Style','popup','String',contrastString,'FontSize',fontSizeSmall-1);

% Analysis Type
analysisTypeString = 'Firing Rate|delta TF|PAC';
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-12*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight*3], ...
    'Style','text','String','Analysis Type','FontSize',fontSizeSmall);
hAnalysisType = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-11.5*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight*3], ...
    'Style','popup','String',analysisTypeString,'FontSize',fontSizeSmall-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Ranges panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rangesHeight = 0.06; rangesTextWidth = 0.01; rangesBoxWidth = 0.6;
hRangePanel = uipanel('Title','Ranges','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[rangesStartPos panelStartHeight rangesPanelWidth panelHeight]);
frequencyRange=[0 100];
timeRange=[-0.2 1];
phaseFrequencyRange=[4 68];
ampFrequencyRange=[7 157];

% ranges
uicontrol('Parent',hRangePanel,'Unit','Normalized', ...
    'Position',[0 1-3.5*rangesHeight 0.4-rangesTextWidth rangesHeight*3.5], ...
    'Style','text','String','Frequency Range','FontSize',fontSizeSmall);
hFreqMin=uicontrol('Parent',hRangePanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[0.4+rangesTextWidth 1-3.5*rangesHeight  0.8-rangesBoxWidth rangesHeight*3.5], ...
    'Style','edit','String',num2str(frequencyRange(1)),'FontSize',fontSizeSmall-1);
hFreqMax=uicontrol('Parent',hRangePanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[0.6+rangesTextWidth 1-3.5*rangesHeight  0.8-rangesBoxWidth rangesHeight*3.5], ...
    'Style','edit','String',num2str(frequencyRange(2)),'FontSize',fontSizeSmall-1);

uicontrol('Parent',hRangePanel,'Unit','Normalized', ...
    'Position',[0 1-7.5*rangesHeight 0.4-rangesTextWidth rangesHeight*3.5], ...
    'Style','text','String','Time Range','FontSize',fontSizeSmall);
hTimeMin=uicontrol('Parent',hRangePanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[0.4+rangesTextWidth 1-7.5*rangesHeight  0.8-rangesBoxWidth rangesHeight*3.5], ...
    'Style','edit','String',num2str(timeRange(1)),'FontSize',fontSizeSmall-1);
hTimeMax=uicontrol('Parent',hRangePanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[0.6+rangesTextWidth 1-7.5*rangesHeight  0.8-rangesBoxWidth rangesHeight*3.5], ...
    'Style','edit','String',num2str(timeRange(2)),'FontSize',fontSizeSmall-1);

% phase frequency range
uicontrol('Parent',hRangePanel,'Unit','Normalized', ...
    'Position',[0 1-11.5*rangesHeight 0.4-rangesTextWidth rangesHeight*3.5], ...
    'Style','text','String','Phase Freq Range','FontSize',fontSizeSmall);
hPhaseFreqMin=uicontrol('Parent',hRangePanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[0.4+rangesTextWidth 1-11.5*rangesHeight  0.8-rangesBoxWidth rangesHeight*3.5], ...
    'Style','edit','String',num2str(phaseFrequencyRange(1)),'FontSize',fontSizeSmall-1);
hPhaseFreqMax=uicontrol('Parent',hRangePanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[0.6+rangesTextWidth 1-11.5*rangesHeight  0.8-rangesBoxWidth rangesHeight*3.5], ...
    'Style','edit','String',num2str(phaseFrequencyRange(2)),'FontSize',fontSizeSmall-1);

% amplitude frequency range
uicontrol('Parent',hRangePanel,'Unit','Normalized', ...
    'Position',[0 1-15.5*rangesHeight 0.4-rangesTextWidth rangesHeight*3.5], ...
    'Style','text','String','Amp Freq Range','FontSize',fontSizeSmall);
hAmpFreqMin=uicontrol('Parent',hRangePanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[0.4+rangesTextWidth 1-15.5*rangesHeight  0.8-rangesBoxWidth rangesHeight*3.5], ...
    'Style','edit','String',num2str(ampFrequencyRange(1)),'FontSize',fontSizeSmall-1);
hAmpFreqMax=uicontrol('Parent',hRangePanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[0.6+rangesTextWidth 1-15.5*rangesHeight  0.8-rangesBoxWidth rangesHeight*3.5], ...
    'Style','edit','String',num2str(ampFrequencyRange(2)),'FontSize',fontSizeSmall-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Timing panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timingHeight = 0.05; timingTextWidth = 0.015; timingBoxWidth = 0.6;
hTimingPanel = uipanel('Title','Limits','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[timingStartPos panelStartHeight timingPanelWidth panelHeight]);

TFZRange = [-6 10];
PACZRange = [0 0.00015];

% parameters
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-3*timingHeight 0.4-timingTextWidth timingHeight*3.5], ...
    'Style','text','String','Parameter','FontSize',fontSizeSmall-1);
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0.4+timingTextWidth 1-3*timingHeight 0.8-timingBoxWidth timingHeight*3.5], ...
    'Style','text','String','Min','FontSize',fontSizeSmall-1);
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0.7+timingTextWidth 1-3*timingHeight 0.8-timingBoxWidth timingHeight*3.5], ...
    'Style','text','String','Max','FontSize',fontSizeSmall-1);

% firing rate
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-7*timingHeight 0.4-timingTextWidth timingHeight*3.5], ...
    'Style','text','String','Firing rate','FontSize',fontSizeSmall-1);
hYMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[0.4+timingTextWidth 1-7*timingHeight 0.8-timingBoxWidth timingHeight*3.5], ...
    'Style','edit','String','0','FontSize',fontSizeSmall-1);
hYMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[0.1+timingBoxWidth 1-7*timingHeight 0.8-timingBoxWidth timingHeight*3.5], ...
    'Style','edit','String','1','FontSize',fontSizeSmall-1);

% color limits for TF
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-11*timingHeight 0.4-timingTextWidth timingHeight*3.5], ...
    'Style','text','String','TF Color Lims','FontSize',fontSizeSmall-1);
hTFZMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[0.4+timingTextWidth 1-11*timingHeight 0.8-timingBoxWidth timingHeight*3.5], ...
    'Style','edit','String',num2str(TFZRange(1)),'FontSize',fontSizeSmall-1);
hTFZMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[0.1+timingBoxWidth 1-11*timingHeight 0.8-timingBoxWidth timingHeight*3.5], ...
    'Style','edit','String',num2str(TFZRange(2)),'FontSize',fontSizeSmall-1);

% color limits for PAC
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-15*timingHeight 0.4-timingTextWidth timingHeight*3.5], ...
    'Style','text','String','PAC Color Lims','FontSize',fontSizeSmall-1);
hPACZMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[0.4+timingTextWidth 1-15*timingHeight 0.8-timingBoxWidth timingHeight*3.5], ...
    'Style','edit','String',num2str(PACZRange(1)),'FontSize',fontSizeSmall-1);
hPACZMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[0.1+timingBoxWidth 1-15*timingHeight 0.8-timingBoxWidth timingHeight*3.5], ...
    'Style','edit','String',num2str(PACZRange(2)),'FontSize',fontSizeSmall-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotOptionsHeight = 0.1;
hPlotOptionsPanel = uipanel('Title','Plotting Options','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[plotOptionsStartPos panelStartHeight plotOptionsPanelWidth panelHeight]);
uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 6*plotOptionsHeight 1 plotOptionsHeight+0.25], ...
    'Style','pushbutton','String','cla','FontSize',fontSizeMedium, ...
    'Callback',{@cla_Callback});
uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 0.05*plotOptionsHeight 1 plotOptionsHeight+0.25], ...
    'Style','pushbutton','String','rescale All','FontSize',fontSizeMedium, ...
    'Callback',{@rescaleData_Callback});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% put plot in a separate small button
actualplotHeight = 1;
hactualplotPanel = uipanel('Title',[],'fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[actualplotStartPos panelStartHeight actualplotPanelWidth panelHeight]);

uicontrol('Parent',hactualplotPanel,'Unit','Normalized', ...
    'Position',[0 0 1 actualplotHeight], ...
    'Style','pushbutton','String','plot','FontSize',fontSizeMedium-1, ...
    'Callback',{@plotData_Callback});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get plot handles
if strcmp(session,'single')
    columns=6;
elseif strcmp(session,'dual')
    columns=8;
end

hFR=getPlotHandles(1,columns,[0.025,0.65,0.78 0.15],0.01,0.01,0);
hAllFR=getPlotHandles(1,1,[0.852,0.65,0.143 0.15],0.01,0.01,0);
hTF=getPlotHandles(2,columns,[0.025,0.36,0.78 0.27],0.01,0.01,0);
hAllTF= getPlotHandles(1,1,[0.852,0.36,0.143 0.27],0.01,0.01,0);
hPAC=getPlotHandles(2,columns,[0.025,0.05,0.78 0.27],0.01,0.01,0);
hAllPAC= getPlotHandles(1,1,[0.852,0.05,0.143 0.27],0.01,0.01,0);

%% Plot firing rate
    function plotData_Callback(~,~)
        f=get(hSpatialFreq,'val');
        o=get(hOrientation,'val');
        c=get(hContrast,'val');
        analysisType = get(hAnalysisType,'val');
        stimulationVal= get(hStimulationType,'val');
        polarityVal=get(hPolarity,'val');
        bandFlag=get(hBand,'val');

        %% construct Protocol List
        protocolID = strcat(monkeyName,stimulationString{stimulationVal}, '_', polarityString{polarityVal}, '_', condition, '_', session);

        %% Load Protocol Information
        expDates = cell(1, 2);
        protocolNames = cell(1, 2);
        protocols = cell(1, 2);
        dates= cell(1, 2);
        for i=1:2
            [expDates{i},protocolNames{i},~] = eval(['allProtocols' protocolID{i}]);
            dates{i}=unique(expDates{i},'stable');
            protN=size(expDates{i},2)./size(dates{i},2); %Gives us the number of Protocols, condition specific

            for day=1:length(dates{i})
                for n=1:protN
                    protocols{1,i}{day,n}=protocolNames{1,i}{1,n+(day-1)*protN}; %This loop rearranges all the protocols in a nice identifiable structure
                end
            end
        end

        folderIn= fullfile(folderSource,'programs','savedData',monkeyName,strcat(session, 'Stim'),stimulationString{stimulationVal}, badTrialNameStr, condition, polarityString{polarityVal});
        [colCode,titleString,stimBlockID]=pickIDs(session);
        % an electrode deciding loop here, LFP(HighRMS) or Spike(goodResp)

        if analysisType==1
            normalizeFlag=1;
            transientFlag=0;
            FRFlag=1;
            PSDFlag=0;
            plotHandle1=hFR;
            plotHandle2=hAllFR;
            for cnd=1:2
                plotSpikeRate(plotHandle1,folderIn(cnd),condition(cnd),session,normalizeFlag,dates(cnd),protocols(cnd),f,c,o)
            end
            text(-0.45,0.9,append('Session=',num2str(size(dates{1,1},2))),'color','k','FontSize',8,'FontWeight','bold');
            legend('','Stim','','Sham','Location','northeast')
            plotSpikePower(plotHandle2,folderIn,session,normalizeFlag,transientFlag,FRFlag,PSDFlag,bandFlag,dates,protocols,f,c,o,colCode,titleString,stimBlockID)
            title(plotHandle2,'\Delta Stim-Sham',FontSize=16)
        elseif analysisType==2
            normalizeFlag=1;
            transientFlag=0;
            FRFlag=0;
            PSDFlag=1;
            plotHandle1=hTF;
            plotHandle2=hAllTF;
            plotTF(plotHandle1,folderIn,dates,protocols,f,c,o,badTrialNameStr)
            plotSpikePower(plotHandle2,folderIn,session,normalizeFlag,transientFlag,FRFlag,PSDFlag,bandFlag,dates,protocols,f,c,o,colCode,titleString,stimBlockID)
        elseif analysisType==3
            plotHandle1=hPAC;
            plotHandle2=hAllPAC;
            plotPAC(plotHandle1,plotHandle2,folderIn,dates,protocols,session,f,c,o,badTrialNameStr,monkeyName,colCode,titleString,stimBlockID)
        end
        for n=1:length(plotHandle1)
            title(hFR(1,n),titleString{1,n},color=colCode{1,n}{1,1},FontSize=16)
        end
    end


    function cla_Callback(~,~)
        claGivenPlotHandle(hFR);
        claGivenPlotHandle(hTF);
        claGivenPlotHandle(hPAC);
        claGivenPlotHandle(hAllFR);
        claGivenPlotHandle(hAllTF);
        claGivenPlotHandle(hAllPAC);
        % cla(hRFMapPlot);cla(hcenterRFMapPlot);

        function claGivenPlotHandle(plotHandle)
            [numRows,numCols] = size(plotHandle);
            for i=1:numRows
                for j=1:numCols
                    cla(plotHandle(i,j));
                end
            end
        end
    end

    function rescaleData_Callback(~,~)

        analysisType = get(hAnalysisType,'val');
        if analysisType==1
            plotHandle=hFR;
            xMax = str2double(get(hTimeMax,'String'));
            xMin = str2double(get(hTimeMin,'String'));
            yMax = str2double(get(hYMax,'String'));
            yMin = str2double(get(hYMin,'String'));
            zMax = [];
            zMin = [];
        elseif analysisType==2
            plotHandle=hTF;
            xMax = str2double(get(hTimeMax,'String'));
            xMin = str2double(get(hTimeMin,'String'));
            yMax = str2double(get(hFreqMax,'String'));
            yMin = str2double(get(hFreqMin,'String'));
            zMax = str2double(get(hTFZMax,'String'));
            zMin = str2double(get(hTFZMin,'String'));
        elseif analysisType==3
            plotHandle=hPAC;
            xMax = str2double(get(hPhaseFreqMax,'String'));
            xMin = str2double(get(hPhaseFreqMin,'String'));
            yMax = str2double(get(hAmpFreqMax,'String'));
            yMin = str2double(get(hAmpFreqMin,'String'));
            zMax = str2double(get(hPACZMax,'String'));
            zMin = str2double(get(hPACZMin,'String'));
        end
        rescaleData(plotHandle,xMin,xMax,yMin,yMax,zMin,zMax,analysisType);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function rescaleData(plotHandle,xMin,xMax,yMin,yMax,zMin,zMax,analysisType)
        [numRows,numCols] = size(plotHandle);
        for i=1:numRows
            for j=1:numCols
                axis(plotHandle(i,j),[xMin xMax yMin yMax]);
                if analysisType==2||analysisType==3
                    clim(plotHandle(i,j),[zMin,zMax])              
                end
                if (i==numRows && rem(j,2)==1)
                    if j~=1
                        set(plotHandle(i,j),'YTickLabel',[]);
                    end
                elseif (rem(i,2)==0 && j==1)
                    set(plotHandle(i,j),'XTickLabel',[]);
                else
                    set(plotHandle(i,j),'XTickLabel',[],'YTickLabel',[]);
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outString = getStringFromValues(valsUnique,decimationFactor)

if length(valsUnique)==1
    outString = convertNumToStr(valsUnique(1),decimationFactor);
else
    outString='';
    for i=1:length(valsUnique)
        outString = cat(2,outString,[convertNumToStr(valsUnique(i),decimationFactor) '|']);
    end
    outString = [outString 'all'];
end

    function str = convertNumToStr(num,f)
        if num > 16384
            num=num-32768;
        end
        str = num2str(num/f);
    end
end


%% Spike plotting code%%
function plotSpikeRate(plotHandle,folderIn,condition,session,normalizeFlag,dates,protocols,sfVals,conVals,oriVals)

if ~exist('sfVals','var'); sfVals=5; end
if ~exist('conVals','var'); conVals=4; end
if ~exist('OriVals','var'); oriVals=5; end

folderIn=cell2mat(folderIn);
if strcmp(condition,'Stim')
    cnd=1;
elseif strcmp(condition,'Sham')
    cnd=2;
end

% Calling the color grid
colCode=pickIDs(session);
colLine={'#ff6f00','#5800FF'};
clear rawPSTH

%     Collecting the data
for iProt =1:size(protocols{1,1},2)
    for j=1:length(sfVals)
        for k=1:length(oriVals)
            for l=1:length(conVals)
                for day=1:length(dates{1,1})
                    dataSpike = dir(fullfile(folderIn,append(num2str(sfVals(j)),'sf','_',num2str(oriVals(k)),'Ori','_',num2str(conVals(l)),'con','_',protocols{1,1}{day,iProt},'_',num2str(dates{1,1}{1,day}),'*PSTH.mat')));

                    % Loading the spike file
                    if ~isempty(dataSpike)
                        s1 = load(fullfile(folderIn,dataSpike.name)); % DataStructure(field).name
                        rawPSTH{day,:,:} = s1.PSTHGrid;
                        spktList=s1.xs; %Spike timevals
                    else
                        rawPSTH{day,:,:} =[];
                    end
                end
                alldayPSTH{j,k,l}=cat(1,rawPSTH{:,1}); %concatenating across days
            end
        end
    end
    PSTHData{iProt,:}=mean(cell2mat(permute(alldayPSTH(:),[2 3 1])),3,'omitnan');%Averaging across all condition specific cells
end

%% Here pulling of all data from required dates has been completed.

for i=1:size(protocols{1,1},2)
    for n=1:size(PSTHData{1,1},1)
        maxFR=max(PSTHData{1,1}(n,:)); % takes the max firing rate of the first protocol, per electrode
        normFRData{i,1}(n,:)=PSTHData{i,1}(n,:)./maxFR;
    end
end

% Plot delta shade PSTH plot
for iProt =1:size(protocols{1,1},2)
    if normalizeFlag==1
        deltaShadeData=normFRData{iProt,:};
    else
        deltaShadeData=PSTHData{iProt,:};
    end
    options.handle = plotHandle(1,iProt);
    % options.color_area = ColCode{1,iProt}{1,1};
    if strcmp(condition,'Stim')
        options.color_area =[255, 111, 0]./255;
    elseif strcmp(condition,'Sham')
        options.color_area =[88, 0, 255]./255;
    end

    options.alpha      = 0.2;
    options.error      = 'sem';
    options.x_axis= spktList;
    plot_areaerrorbar(deltaShadeData, options);hold on;
    axis(plotHandle(1,iProt),[[-0.5 1] [-0.2 1]]);
end

for iProt=1:size(protocols{1,1},2)
    if ~isempty(PSTHData{iProt,:})
        if normalizeFlag==1
            plot(plotHandle(1,iProt),spktList,mean(normFRData{iProt,:},1,'omitnan'),LineWidth=1.6, color=colLine{1,cnd});hold (plotHandle(1,iProt),'on');
            if iProt==size(protocols{1,1},2)
                axes(plotHandle(1,iProt))
            end
        else
            plot(plotHandle(1,iProt),spktList,mean(PSTHData{iProt,:},1,'omitnan'),LineWidth=1.6, color=colCode{1,iProt}{1,2});hold (plotHandle(1,iProt),'on');
            if iProt==size(protocols{1,1},2)
                axes(plotHandle(1,iProt))
            end
        end
    end
end
set(plotHandle(1,2:end),'ytick',[])
set(plotHandle(1,1:end),'xtick',[])
end


%% TF plotting code %%
function plotTF(plotHandle,folderIn,dates,protocols,sfVals,conVals,oriVals,badTrialNameStr)

if ~exist('sfVals','var'); sfVals=5; end
if ~exist('conVals','var'); conVals=4; end
if ~exist('OriVals','var'); oriVals=5; end


if strcmp(badTrialNameStr,'V1') % We have to have a grid of different range, to comply with convalues
    SGRange={{12 28}, {16,28}, {28,36}, {16,28}};
    FGRange={{32,44}, {32,48}, {48,68}, {36,52}};
    cLimsDiff=[-6 10];
elseif strcmp(badTrialNameStr,'V4')
    SGRange={{16 28}, {16,28}, {28,40}, {20,40}};
    FGRange={{32,40}, {32,44}, {60,76}, {52,72}};
    cLimsDiff=[-8 8];
end

%% Getting Data
for c=1:2
    for iProt =1:size(protocols{1,c},2)
        for j=1:length(sfVals)
            for k=1:length(oriVals)
                for l=1:length(conVals)
                    for day=1:length(dates{1,c})
                        dataTF = dir(fullfile(folderIn{1,c},append(num2str(sfVals(j)),'SF','_',num2str(oriVals(k)),'Ori','_',num2str(conVals(l)),'Con','_',protocols{1,c}{day,iProt},'_',dates{1,c}{1,day},'*TF.mat')));
                        % Loading the TF file
                        d2= load(fullfile(folderIn{1,c},dataTF.name)); % DataStructure(field).name
                        rawTF(day,:,:) = (d2.TFDeltaPow);
                        tList=d2.tList; % Time frequency list
                        fList=d2.fList; % Frequency point list
                    end
                    allday{j,k,l}=squeeze(mean(rawTF,1,"omitnan"));%Averaging across days took place here
                end
            end
        end
        TFData(iProt,:,:)=mean(cell2mat(permute(allday(:),[2 3 1])),3);%Averaging across all condition specific cells
    end

    %% Plot TimeFrequency plot
    for iProt=1:size(protocols{1,c},2)
        pcolor(tList,fList,squeeze(TFData(iProt,:,:))','Parent',plotHandle(c,iProt));colormap('jet');
        shading(plotHandle(c,iProt),'interp');
        clim(plotHandle(c,iProt),cLimsDiff);
        axes(plotHandle(c,iProt))
        axis(plotHandle(1:2,iProt),[[-0.5 1] [0 100]]);
        Xax=gca().XAxis;
        Yax=gca().YAxis;
        set(gca,'FontWeight','bold');
        set(Xax,'FontSize',10);
        set(Yax,'FontSize',10);
        yline(cell2mat(SGRange{1,conVals}),"--");
        yline(cell2mat(FGRange{1,conVals}),"-");
    end
end
set(plotHandle(1,1:size(protocols{1,c},2)),'XTick',[]);
set(plotHandle(2,1:size(protocols{1,c},2)),'XTick',[0 0.5 1]);
set(plotHandle(1:end,2:end),'ytick',[])
ax=plotHandle(2,size(protocols{1,c},2));
cb=colorbar(ax);
cb.Position=[0.8064 0.359 0.0064 0.1302];
cb.TickDirection = 'none';
cb.Ticks = [-6 0 10];
end


%% This code will plot PAC values with Frequency(Phase) in X axis and Frequency(Amplitude) in Y axis
function plotPAC(plotHandle1,plotHandle2,folderIn,dates,protocols,session,sfVals,conVals,OriVals,badTrialNameStr,monkeyName,colCode,titleString,stimBlockID)
% For better visualization: splitting PAC plot into two parts: 2 to 150 Hz and 150Hz to 500Hz
partitionFreq = [7 157 157 487];

% V1=1:48;
% V4=49:96;
% BrainArea={V1,V4};
nWin=1;

for ses=1:length(stimBlockID)
    colCode{1,stimBlockID(ses)}=[];
    titleString{1,stimBlockID(ses)}=[];
end

Col=colCode(~cellfun('isempty',colCode));
titl=titleString(~cellfun('isempty',colCode));
colCode=reshape(Col,[1 length(colCode)-length(stimBlockID)]);
titleString=reshape(titl,[1 length(titleString)-length(stimBlockID)]);


% rfDataFileName = [monkeyName 'MicroelectrodeRFData.mat']; %This file is in DataMap/ReceptiveFieldData/{monkeyName} folder and should be in Matlab's path
% if exist(rfDataFileName,'file')
%     tmp = load(rfDataFileName);
%     electrodesToUse = tmp.highRMSElectrodes;
% end
%
% % Finding the brain area
% if strcmp(badTrialNameStr,'V1')
%     AreaFlag=1;
% elseif strcmp(badTrialNameStr,'V4')
%     AreaFlag=2;
% end

if strcmp(badTrialNameStr,'V1') % We have to have a grid of different range, to comply with con values
    % SGRange={{12 28}, {16,28}, {28,36}, {16,28}};
    % FGRange={{32,44}, {32,48}, {48,68}, {36,52}};
    SGRange={{16 28}, {18,32}, {24,40}, {16,28}};
    HGRange={{40,88},{50,98},{80,120},{80,120}};
elseif strcmp(badTrialNameStr,'V4')
    SGRange={{16 28}, {16,28}, {28,40}, {20,40}};
end

stimFreqWin = {[SGRange{1,conVals}  HGRange{1,conVals}] };

if mod(stimFreqWin{1,nWin}{1},4)==1
    xFreq1=stimFreqWin{1,nWin}{1}-1;
elseif mod(stimFreqWin{1,nWin}{1},4)==2
    xFreq1=stimFreqWin{1,nWin}{1}+2;
elseif mod(stimFreqWin{1,nWin}{1},4)==3
    xFreq1=stimFreqWin{1,nWin}{1}+1;
elseif mod(stimFreqWin{1,nWin}{1},4)==0
    xFreq1=stimFreqWin{1,nWin}{1};
end

if mod(stimFreqWin{1,nWin}{2},4)==1
    xFreq2=stimFreqWin{1,nWin}{2}-1;
elseif mod(stimFreqWin{1,nWin}{2},4)==2
    xFreq2=stimFreqWin{1,nWin}{2}+2;
elseif mod(stimFreqWin{1,nWin}{2},4)==3
    xFreq2=stimFreqWin{1,nWin}{2}+1;
elseif mod(stimFreqWin{1,nWin}{2},4)==0
    xFreq2=stimFreqWin{1,nWin}{2};
end

% goodElectrodes = intersect(electrodesToUse,BrainArea{1,AreaFlag});

% cLimsDiff = [0 1];   % colormap limits for change in power
colormap jet
climVals1=0;
climVals2=0.00015;

%Pre-defined determinants
% gridType='Microelectrode';
removeEvokedResponseFlag=1;
tapers=[1 1];
modality='LFP';
sVarName='sf';
pacMethod='klmi';
filterName='fir';
nSurrogates=0;
useMPFlag=0;

blGrid=[];
stGrid=[];
for c=1:2
    for day=1:length(dates{1,c})
        for iProt =1:size(protocols{1,c},2)
            % if iProt==1
            %     respFile = fullfile(folderin{1,c},'data',monkeyName,gridType,dates{1,c}{1,day},protocols{1,c}{day,iProt},'segmentedData', append('GoodUnits',badTrialNameStr,'.mat'));
            %     if isfile(respFile)
            %         respFileStruct=load(respFile);
            %         AllGoodUnit{1,iProt}=respFileStruct.goodSpikeElectrodes;
            %         AreaGoodUnit{1,iProt}=intersect(BrainArea{1,AreaFlag},AllGoodUnit{1,iProt});
            %     end
            % end

            % [~,elecId]=find(ismember(goodElectrodes,respFileStruct.goodSpikeElectrodes));

            %tmpData struc=(Elecnum,numPeriods,Amplitude bins, PhaseBins);
            tmpData = load(fullfile(folderIn{1,c},[monkeyName dates{1,c}{1,day} protocols{1,c}{day,iProt} '_removeMean' num2str(removeEvokedResponseFlag) '_Tapers' num2str(tapers(1)) '_' num2str(tapers(2)) '_' modality '_' sVarName num2str(sfVals) '_o' num2str(OriVals) '_c' num2str(conVals) '_' pacMethod '_' filterName '_nSur' num2str(nSurrogates) '_MP' num2str(useMPFlag) badTrialNameStr '.mat']), 'pac','normalisedPac','meanAmp','surrogatePac','tval','pval', 'centerAmpFreq', 'centerPhaseFreq');
            blGrid{day,iProt}=tmpData.pac(:,1,:,:);
            xFreqPos = intersect(find(tmpData.centerPhaseFreq>= xFreq1),find(tmpData.centerPhaseFreq<=  xFreq2));
            yFreqPos = intersect(find(tmpData.centerAmpFreq>= stimFreqWin{1,nWin}{3}),find(tmpData.centerAmpFreq<= stimFreqWin{1,nWin}{4}));
            stGrid{day,iProt}=tmpData.pac(:,2,:,:);
            pacSig{c,day,iProt} = squeeze(mean(tmpData.pac(:,2,yFreqPos,xFreqPos),[3,4]));
        end
    end

    % Initialize the result as a 1x6 cell array
    blResult = cell(1, size(protocols{1,c},2));
    stResult = cell(1, size(protocols{1,c},2));

    % Loop through each column to concatenate along the first dimension
    for col = 1:size(protocols{1,c},2)
        blResult{col} = cat(1, blGrid{:, col});
        stResult{col} = cat(1, stGrid{:, col});
    end

    avgblResult = cell(1, size(protocols{1,c},2));
    avgstResult = cell(1, size(protocols{1,c},2));

    % Loop through each cell and compute the mean along the first dimension
    for col = 1:size(protocols{1,c},2)
        avgblResult{col} = squeeze(mean(blResult{col}, 1)); % Compute the mean across the 1st dimension (size becomes 1x1x49x17)
        avgstResult{col} = squeeze(mean(stResult{col}, 1));
    end
    freqIdx = find(tmpData.centerAmpFreq == partitionFreq(2));
    for n=1:size(protocols{1,c},2) % n=Pre or Post
        pcolor(plotHandle1(c,n),tmpData.centerPhaseFreq,tmpData.centerAmpFreq(1:freqIdx),avgstResult{1,n}((1:freqIdx),:)-avgblResult{1,n}((1:freqIdx),:));
        shading(plotHandle1(c,n),'interp');
        clim(plotHandle1(c,n), [climVals1 climVals2]);
        set(plotHandle1(1:2,:),'Yscale','log');
        set(plotHandle1(1:2,:),'Xscale','log');
        if n==1
            set(plotHandle1(c,n),'yTick',[ 10 50 100 150],'FontSize',8,'FontWeight','bold'); set(plotHandle1(2,n),'yTickLabel',[ 10 50 100 150],'FontSize',8,'FontWeight','bold');
        else
            set(plotHandle1(1,n),'yTick',[150 250 350 450],'FontSize',8,'FontWeight','bold'); set(plotHandle1(1,n),'yTickLabel',[],'FontSize',8,'FontWeight','bold');
            set(plotHandle1(2,n),'yTick',[ 10 50 100 150],'FontSize',8,'FontWeight','bold'); set(plotHandle1(2,n),'yTickLabel',[],'FontSize',8,'FontWeight','bold');
        end
    end
    if c==1
        set(plotHandle1(1:2,:),'xTickLabel',[])
    elseif c==2
        set(plotHandle1(2,:),'xTick',[ 5 10 24 40],'FontSize',8,'FontWeight','bold'); set(plotHandle1(2,:),'xTickLabel',[ 5 10 24 40],'FontSize',8,'FontWeight','bold')
    end
    eleclength=size(stGrid{1,1},1); % for later use
    clear stGrid blGrid tmpData avgstResult avgblResult
    % clear plotHandle
end
ax=plotHandle1(2,size(protocols{1,c},2));
cb=colorbar(ax);
cb.Position=[0.8064 0.0499 0.0064 0.1302];
cb.Ruler.Exponent = -3;
cb.Ticks = [0 0.00015];
cb.TickDirection = 'none';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numProtocols = size(protocols{1}, 2);
if strcmp(session,'single')
    % Initialize the concatenated data array (2x128x6)
    concatenatedData = zeros(2, eleclength*length(dates{1,c}), numProtocols);
    data=pacSig;
    % Iterate over the first and third dimensions
    for i = 1:2
        for prot = 1:numProtocols
            % Extract and concatenate along the second dimension
            concatenated = cat(2, data{i, :, prot}); % concatenate across the day cells
            concatenated = reshape(concatenated, [], 1); % Reshape into 128x1
            concatenatedData(i, :, prot) = concatenated; % Store in 2x128x6 array
        end
    end

    for prot=1:numProtocols
        St=concatenatedData(1,:,prot);
        Sh=concatenatedData(2,:,prot);
        dMI{1,prot}=St-Sh;
        stError{1,prot}=std(dMI{1,prot})/sqrt(size(dMI{1,prot},2));
        avgMI{1,prot}=mean(dMI{1,prot});
    end
    point=1:numProtocols-1;
    avgMI(2)=[];
    stError(2)=[];
    dMI(2)=[];

elseif strcmp(session,'dual')
    data = pacSig;
    numStimDays = length(dates{1,1});
    numShamDays = length(dates{1,2});

    %concatenate Data
    for prot = 1:numProtocols
        stimConcat = reshape(cat(2, data{1, :, prot}), [], 1);
        shamConcat = reshape(cat(2, data{2, :, prot}), [], 1);

        concatenatedStimData(:, prot) = stimConcat;
        concatenatedShamData(:, prot) = shamConcat;
    end

    %Slice Each Day’s Data into Cells
    for prot = 1:numProtocols
        % Stim days
        for d = 1:numStimDays
            idx = (d-1)*eleclength + 1;
            U{prot, d} = concatenatedStimData(idx:idx+eleclength-1, prot);
        end
        % Sham days
        for d = 1:numShamDays
            idx = (d-1)*eleclength + 1;
            V{prot, d} = concatenatedShamData(idx:idx+eleclength-1, prot);
        end
    end

    % Dynamically Average Across Days
    for prot = 1:numProtocols
        % Stim: average across all U{prot, :}
        stimDays = [U{prot, 1:numStimDays}];
        avg_Stim{prot} = mean(stimDays, 2);

        % Sham: average across all V{prot, :}
        shamDays = [V{prot, 1:numShamDays}];
        avg_Sham{prot} = mean(shamDays, 2);
    end

    % Compute Differences, Std Error, Mean
    for prot = 1:numProtocols
        dMI{1, prot} = avg_Stim{prot} - avg_Sham{prot};
        stError{1, prot} = std(dMI{1, prot}) / sqrt(length(dMI{1, prot}));
        avgMI{1, prot} = mean(dMI{1, prot});
    end
    point=1:numProtocols-2;
    avgMI(2)=[];% making first stimblock(prot=2) zero
    avgMI(3)=[];% making second stimblock(prot=4) zero, as it is now 3rd protocol after deletion
    stError(2)=[];
    stError(3)=[];
    dMI(2)=[];
    dMI(3)=[];
end

%% Plotting data
if strcmp(session,'single')
    numProtocols=numProtocols-1;
elseif strcmp(session,'dual')
    numProtocols=numProtocols-2;
end

for prot=1:numProtocols
    errorbar(subplot(plotHandle2),point(prot),avgMI{1,prot}-avgMI{1,1},[stError{1,prot}],'color','#ED4672','LineWidth',1.5);
    hold on
end

xlim([0 numProtocols+1]);
ylim([-5*10^-5 5*10^-5])
plot(point,cell2mat(avgMI)-cell2mat(avgMI(1,1)),'color','#ED4672','LineWidth',1.5);
yline(0*10^-5,LineWidth=0.8,LineStyle="--")

for prot=1:numProtocols
    w{1,prot}=avgMI{1,prot}-avgMI{1,1};
    scatter(point(1,prot),w{1,prot},30,colCode{1, prot}{1, 1},"filled")
end


%% Collecting band specific stats
for Band=nWin
    for ProtN=2:numProtocols
        x1=dMI{1,1};% Post-Pre
        y1=dMI{1,ProtN};
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
for freqstat=nWin %Band specific identity
    m=[w{:}];
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

    delete(findall(subplot(plotHandle2(1,1)), 'Tag', 'SigStar'));
    for d=1:numProtocols
        if stat1(1,d)==1
            hold on
            subplot(plotHandle2(1,1))
            yPoint=(asterixData(:,d));
            if yPoint<0
                if stat2(1,d)<0.0005
                    text(xCentres(:,d)-0.3,yPoint-(cell2mat(stError(freqstat,d))+0.4*10^-5),'\ast\ast\ast','fontWeight','bold','tag','SigStar',HandleVisibility='off'); hold on
                elseif stat2(1,d)<0.005
                    text(xCentres(:,d)-0.2,yPoint-(cell2mat(stError(freqstat,d))+0.4*10^-5),'\ast\ast','fontWeight','bold','tag','SigStar',HandleVisibility='off'); hold on
                elseif stat2(1,d)<0.05
                    text(xCentres(:,d)-0.1,yPoint-(cell2mat(stError(freqstat,d))+0.4*10^-5),'\ast','fontWeight','bold','tag','SigStar',HandleVisibility='off'); hold on
                end
            elseif yPoint>0
                if stat2(1,d)<0.0005
                    text(xCentres(:,d)-0.35,yPoint+(cell2mat(stError(freqstat,d))+0.7*10^-5),'\ast\ast\ast','fontWeight','bold','tag','SigStar',HandleVisibility='off'); hold on
                elseif stat2(1,d)<0.005
                    text(xCentres(:,d)-0.25,yPoint+(cell2mat(stError(freqstat,d))+0.7*10^-5),'\ast\ast','fontWeight','bold','tag','SigStar',HandleVisibility='off'); hold on
                elseif stat2(1,d)<0.05
                    text(xCentres(:,d)-0.15,yPoint+(cell2mat(stError(freqstat,d))+0.7*10^-5),'\ast','fontWeight','bold','tag','SigStar',HandleVisibility='off'); hold on
                end
            end
        end
    end
end

set(plotHandle2(1,1),'XTick',1:length(point));
xticklabels(titleString);
line([0 length(point)+1], [5*10^-5 5*10^-5],'color','k')
line([length(point)+1 length(point)+1], [-5*10^-5 5*10^-5],'color','k')
end


