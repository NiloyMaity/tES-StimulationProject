%%Please keep the Stimulation Specific Protocol List in MATLAB path.
%Area=V1 or V4, so AreaFlag= 1 or 2 corresponds to V1 & V4
%StimulationType={'tDCS','tACS'}
%condition={'Stim','Sham'};
%Polarity={'Cathodal','Anodal' or 'SG','FG','Alpha'};
%Session={'single','dual','dual60'};
%SessionID={0,1,2}

%% Inputs
MonkeyName='dona';
session='single';
badTrialNameStr='V1';
folderSource='Z:\Students\Niloy\MonkeyData';
SF=[1,2,3,4];
Ori=[1,2,3,4];
Con=[1,2,3];

%% Presets
stimulationString={'tACS','tDCS'};
condition={'Stim','Sham'};
polarityString={'SG','FG','Alpha','Cathodal'};
bandString={'Slow Gamma','Fast Gamma'};

%% Define Brain Area(for now only V1)
V1 = 1:48;
V4 = 49:96;
BrainArea = {V1, V4};
AreaFlag=1; % for V1

displayStimulationAll(MonkeyName,folderSource,stimulationString, condition, polarityString,session,bandString,SF, Con, Ori,badTrialNameStr)
