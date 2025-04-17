%% This function saves PSD,Band Powers,TF and PSTH named by SF,Con, Ori and polarity of stimulation
%%Please keep the Stimulation Specific Protocol List in MATLAB path.
%Area=V1 or V4, so AreaFlag= 1 or 2 corresponds to V1 & V4
%StimulationType={'tDCS','tACS'}
%condition={'Stim','Sham'};
%Polarity={'Cathodal','Anodal' or 'SG','FG','Alpha'};
%Session={'single','dual','dual60'};
%SessionID={0,1,2}

function [BandData, ShadeData, PSTHGrid] = processLFPSpikeData(MonkeyName, folderSource, gridType, AreaFlag,StimulationType, condition, Polarity, SessionID, PSDTFFlag, SFVals, ConVals, OriVals, badTrialNameStr)

if nargin < 10, SFVals = 1:5; end
if nargin < 11, ConVals = 1:4; end
if nargin < 12, OriVals = 1:5; end
if nargin < 13, badTrialNameStr = 'V1'; end

%% Define Brain Areas
V1 = 1:48;
V4 = 49:96;
BrainArea = {V1, V4};

%% Define Stimulation Sessions
SessionTypes = {'single', 'dual', 'dual60'};
if SessionID < 0 || SessionID > 2
    error('Invalid SessionID. Must be 0, 1, or 2.');
end
Session = SessionTypes{SessionID + 1};

%% Define Frequency Ranges Based on Condition
[SGRange, FGRange] = defineFrequencyRanges(badTrialNameStr);

%% Construct Protocol List
protocolID = strcat(StimulationType, '_', Polarity, '_', condition, '_', Session);
SaveFolder = fullfile(folderSource, 'Programs', 'Saved Data', MonkeyName, strcat(Session, '_Stim'));
createDirectory(SaveFolder);

%% Load Protocol Information
[expDates,protocolNames,~] = eval(['allProtocols' protocolID]);
dates=unique(expDates,'stable');
ProtN=size(expDates,2)./size(dates,2); %Gives us the number of Protocols, condition specific

clear protocols
for day=1:length(dates)
    for n=1:ProtN
        protocols{day,n}=protocolNames{1,n+(day-1)*ProtN}; %This loop rearranges all the protocols in a nice identifiable structure
    end
end

%% Load High RMS Electrodes (if available)
rfDataFile = fullfile(folderSource, 'RMS_Cutoff', [MonkeyName gridType 'RFData.mat']);
electrodesToUse = loadElectrodeData(rfDataFile);
electrodeList = intersect(electrodesToUse, BrainArea{AreaFlag});

%% Main Processing Loop
for SF = SFVals
    for Ori = OriVals
        for Con = ConVals
            %% Determine Gamma Ranges
            SGIdx = getFrequencyIndex(SGRange, Con);
            FGIdx = getFrequencyIndex(FGRange, Con);

            %% Prepare Data Source Folder
            dataout = fullfile(SaveFolder, StimulationType, badTrialNameStr, condition, Polarity);
            createDirectory(dataout);

            %% Process Each Experiment Date
            for dayIdx = 1:numel(dates)
                currentDate = dates{dayIdx};
                protocolList = protocols(dayIdx, :);

                %% Identify Good Electrodes
                AreaGoodUnit = getGoodUnits(folderSource, MonkeyName, gridType, currentDate, ...
                    BrainArea, AreaFlag, badTrialNameStr);
                %% Process Each Protocol
                for iProt = 1:numel(protocolList)
                    protocol = protocolList{iProt};
                    dataPath = fullfile(folderSource, 'data', MonkeyName, gridType, currentDate, protocol);

                    %% Load LFP & Stimulus Info
                    [timeVals, Fs,  HighImpElec] = loadLFPandImpedance(dataPath, BrainArea{AreaFlag});

                    %% Get Trial Data
                    [goodPos, electrodeList] = getValidTrials(dataPath, electrodeList, HighImpElec, badTrialNameStr, SF, Ori, Con);

                    %% Compute Power Spectral Density (PSD), Time-Frequency (TF)  and Spike Data
                    [ShadeData,f1,TFDeltaPow, TFStPower, fList, tList, BandData,CollectStSG, CollectBlSG, CollectStFG, CollectBlFG, CollectSG, CollectFG,PSTHGrid,xs,BandDataFileName,PSDFileName,TFFileName,PSTHFileName] = computePSDTFSpike(electrodeList, dataPath, dataout, timeVals, Fs, goodPos, SGIdx, FGIdx, ...
                        AreaGoodUnit, SF, Ori, Con, protocol, currentDate, PSDTFFlag);

                    %% Save Data
                    saveData(ShadeData,f1,TFDeltaPow, TFStPower, fList, tList, BandData,CollectStSG, CollectBlSG, CollectStFG, CollectBlFG, CollectSG, CollectFG,PSTHGrid,xs,BandDataFileName,PSDFileName,TFFileName,PSTHFileName);
                end
            end
        end
    end
end
end

%% Define Frequency range
function [SGRange, FGRange] = defineFrequencyRanges(region)
if strcmp(region, 'V1')
    SGRange = {[12, 28], [16, 28], [28, 36], [16, 28]};
    FGRange = {[32, 44], [32, 48], [48, 68], [36, 52]};
else
    SGRange = {[16, 28], [16, 28], [28, 40], [20, 40]};
    FGRange = {[32, 40], [32, 44], [60, 76], [52, 72]};
end
end

%% Get to electrodes
function electrodes = loadElectrodeData(filename)
if exist(filename, 'file')
    data = load(filename);
    electrodes = data.highRMSElectrodes;
else
    electrodes = [];
end
end

%% Getting frequency range
function idx = getFrequencyIndex(range, contrast)
idx = [(range{contrast}(1) / 2) + 1, (range{contrast}(2) / 2) + 1];
end
function createDirectory(folderPath)
if ~exist(folderPath, 'dir')
    mkdir(folderPath);
end
end

%% Getting timeVals, Impedance
function [timeVals, Fs, HighImpElec] = loadLFPandImpedance(folder, brainArea)
load(fullfile(folder, 'segmentedData', 'LFP','lfpInfo.mat'), 'timeVals');
Ts = timeVals(2) - timeVals(1);
Fs = round(1 / Ts);
parentfolder=fileparts(folder);
impedanceFile = fullfile(parentfolder, 'impedanceValues.mat');
if exist(impedanceFile, 'file')
    impedanceData = load(impedanceFile);
    HighImpElec = intersect(find(impedanceData.impedanceValues > 2500), brainArea);
else
    HighImpElec = [];
end
end

%% Get Good units
function  AreaGoodUnit = getGoodUnits(folderSource, MonkeyName, gridType, currentDate, ...
    BrainArea, AreaFlag, badTrialNameStr)
% for iProt =1:size(protocols,2)
respFile = fullfile(folderSource,'data',MonkeyName,gridType,currentDate,'GRF_001','segmentedData', append('GoodUnits',badTrialNameStr,'.mat'));
if isfile(respFile)
    respFileStruct=load(respFile);
    AllGoodUnit{1,1}=respFileStruct.goodSpikeElectrodes;
    AreaGoodUnit{1,1}=intersect(BrainArea{1,AreaFlag},AllGoodUnit{1,1});
else
    AreaGoodUnit=[];
    disp 'No unit found';
end
% end
end

%% Get good valid trials 
function [goodPos, electrodeList] = getValidTrials(dataPath, electrodeList, HighImpElec, badTrialNameStr, SF, Ori, Con)
% Load necessary files with checks
goodStimFile = fullfile(dataPath, 'extractedData', 'goodStimNums.mat');
paramCombFile = fullfile(dataPath, 'extractedData', 'parameterCombinations.mat');
badTrialsFile = fullfile(dataPath, 'segmentedData', append('badTrials', badTrialNameStr, '.mat'));

if exist(goodStimFile, 'file')
    load(goodStimFile, 'goodStimNums');
else
    error('File not found: %s', goodStimFile);
end

if exist(paramCombFile, 'file')
    load(paramCombFile, 'parameterCombinations');
else
    error('File not found: %s', paramCombFile);
end

% Load bad trials if file exists; otherwise, assume no bad trials
if exist(badTrialsFile, 'file')
    load(badTrialsFile, 'badTrials');
else
    badTrials = []; % No bad trials if file is missing
end

% Remove high-impedance electrodes
electrodeList = setdiff(electrodeList, HighImpElec);

% Get valid trials based on parameter combinations
goodPos = parameterCombinations{1,1,1,SF,Ori,Con,1};

% If `badTrials` and `goodStimNums` are of the same length, ignore badTrials
if length(goodStimNums) == length(badTrials)
    badTrials = [];
end

% Remove bad trials from goodPos
goodPos = setdiff(goodPos, badTrials);

end

%% Coumpute Spike, TF and PSD
function [ShadeData,f1,TFDeltaPow, TFStPower, fList, tList, BandData,CollectStSG, CollectBlSG, CollectStFG, CollectBlFG, CollectSG, CollectFG,PSTHGrid,xs,BandDataFileName,PSDFileName,TFFileName,PSTHFileName] = computePSDTFSpike(electrodeList, dataPath, dataout, timeVals, Fs, goodPos, SGIdx, FGIdx, ...
    AreaGoodUnit, SF, Ori, Con, protocol, currentDate, PSDTFFlag)
% Computes Power Spectral Density (PSD) and Time-Frequency (TF) analysis of spike and LFP data.
% Initialize output
ShadeData = [];
BandData = [];

% MTM Parameters
BWList = 1;                      % Bandwidth: 1 for STFT, 2 for MTM
baselineS = [-0.5, 0];           % Baseline period
stimPeriod = [0.25, 0.75];       % Stimulation period
params = struct(...
    'pad', -1, ...
    'Fs', Fs, ...
    'fpass', [0 100], ...
    'trialave', 1, ...
    'tapers', [BWList, 2 * BWList - 1] ...
    );
movingWinTF = [0.250 0.025];  % Time window for TF analysis

% Baseline & Stimulus Positions
blPos = sum(timeVals < baselineS(1)) + (1:round(diff(baselineS) * Fs));
stPos = sum(timeVals < stimPeriod(1)) + (1:round(diff(stimPeriod) * Fs));

%%% Spike Data Processing %%%
r = cell2mat(AreaGoodUnit);
if ~isempty(r)
    for chan = 1:length(r)
        x = load(fullfile(dataPath,'segmentedData','Spikes', sprintf('elec%d_SID0.mat', r(chan))));
        [psthVals, xs] = getPSTH(x.spikeData(goodPos), 10, [timeVals(1), timeVals(end)]);
        PSTHGrid(chan, :) = psthVals;
    end

    % Save PSTH Data
    PSTHFileName = fullfile(dataout, ...
        sprintf('%dSF_%dOri_%dCon_%s_%s_PSTH.mat', SF, Ori, Con, protocol, currentDate));
    % if exist("PSTHGrid", "var") && ~isempty(PSTHGrid)
    %     save(PSTHFileName, "PSTHGrid", "xs");
    % end
end

%%% Compute PSD and TF Spectrograms %%%
if PSDTFFlag == 1
    numElectrodes = length(electrodeList);
    SpectraPSDBl = zeros(numElectrodes, 51);
    SpectraPSDSt = zeros(numElectrodes, 51);

    % Preallocate TF storage
    SpectraTFBl = zeros(numElectrodes, 1,26);
    SpectraTFSt = zeros(numElectrodes, 72,26);

    for iElec = 1:numElectrodes
        x = load(fullfile(dataPath,'segmentedData','LFP', sprintf('elec%d.mat', electrodeList(iElec))));
        ElecData = x.analogData(goodPos, :);
        ElecData=ElecData-repmat(mean(ElecData,2),1,size(ElecData,2)); %DC shift correction

        PSDdataBL = ElecData(:, blPos);
        PSDdataST = ElecData(:, stPos);

        % Compute Time-Frequency Spectrum
        [TFSpec, tList, fList] = mtspecgramc(ElecData', movingWinTF, params);
        tList = tList + timeVals(1) - 1 / Fs;

        %Finds Baseline position
        BlPosList = tList >= baselineS(1) & tList < baselineS(2);
        SpectraTFBl(iElec, :, :) = mean(TFSpec(BlPosList, :), 1);
        SpectraTFSt(iElec, :, :) = TFSpec;

        % Compute PSD
        [S1, f1] = mtspectrumc(PSDdataBL', params);
        [S2, ~] = mtspectrumc(PSDdataST', params);

        SpectraPSDBl(iElec, :) = S1;
        SpectraPSDSt(iElec, :) = S2;
    end

    % Compute log power differences
    ShadeData = 10 * (log10(SpectraPSDSt) - log10(SpectraPSDBl));
    TFBlPower = squeeze(mean(log10(SpectraTFBl), 1));
    TFStPower = squeeze(mean(log10(SpectraTFSt), 1));
    TFDeltaPow= 10*(TFStPower-repmat(TFBlPower',length(tList),1));

    % Save PSD and TF data
    PSDFileName = fullfile(dataout, ...
        sprintf('%dSF_%dOri_%dCon_%s_%s_PSD.mat', SF, Ori, Con, protocol, currentDate));
    TFFileName = fullfile(dataout, ...
        sprintf('%dSF_%dOri_%dCon_%s_%s_TF.mat', SF, Ori, Con, protocol, currentDate));

    % save(PSDFileName, "ShadeData", "f1");
    % save(TFFileName, "TFDeltaPow", "TFStPower", "fList", "tList");

    % Compute and save band data
    CollectBlSG = log10(mean(SpectraPSDBl(:, SGIdx(1):SGIdx(2)), 2));
    CollectBlFG = log10(mean(SpectraPSDBl(:, FGIdx(1):FGIdx(2)), 2));
    CollectStSG = log10(mean(SpectraPSDSt(:, SGIdx(1):SGIdx(2)), 2));
    CollectStFG = log10(mean(SpectraPSDSt(:, FGIdx(1):FGIdx(2)), 2));

    CollectSG = 10 * (CollectStSG - CollectBlSG);
    CollectFG = 10 * (CollectStFG - CollectBlFG);
    AvgSG = mean(CollectSG, 2);
    AvgFG = mean(CollectFG, 2);
    BandData = [AvgSG, AvgFG];

    BandDataFileName = fullfile(dataout, ...
        sprintf('%dSF_%dOri_%dCon_%s_BandData.mat', SF, Ori, Con, currentDate));
    % save(BandDataFileName, 'BandData', 'CollectStSG', 'CollectBlSG', 'CollectStFG', 'CollectBlFG', 'CollectSG', 'CollectFG');
end
end

%% Save all data
function saveData(ShadeData,f1,TFDeltaPow, TFStPower, fList, tList, BandData,CollectStSG, CollectBlSG, CollectStFG, CollectBlFG, CollectSG, CollectFG,PSTHGrid,xs,BandDataFileName,PSDFileName,TFFileName,PSTHFileName)
save(BandDataFileName, 'BandData', 'CollectStSG', 'CollectBlSG', 'CollectStFG', 'CollectBlFG', 'CollectSG', 'CollectFG');
save(PSDFileName, "ShadeData", "f1");
save(TFFileName, "TFDeltaPow", "TFStPower", "fList", "tList");
save(PSTHFileName, "PSTHGrid", "xs");
end
