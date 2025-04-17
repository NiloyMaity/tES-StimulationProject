%please keep the Stimulation Specific Protocol List in MATLAB path.
%area=V1 or V4, so areaFlag= 1 or 2 corresponds to V1 & V4
%stimulationType={'tDCS','tACS'}
%condition={'Stim','Sham'};
%polarity={'Cathodal','Anodal' or 'SG','FG','Alpha'};
%session={'single','dual','dual60'};
%sessionID={0,1,2}

%% inputs
monkeyName='dona';
session='single';
badTrialNameStr='V1';
folderSource='Z:\Projects\Niloy_tES-StimulationProject';
sf=[1,2,3,4];
ori=[1,2,3,4];
con=[1,2,3];

%% presets
stimulationString={'tACS','tDCS'};
condition={'Stim','Sham'};
polarityString={'SG','FG','Alpha','Cathodal'};
bandString={'Slow Gamma','Fast Gamma'};

%% define Brain Area(for now only V1)
V1 = 1:48;
V4 = 49:96;
brainArea = {V1, V4};
areaFlag=1; % for V1

displayStimulationAll(monkeyName,folderSource,stimulationString, condition, polarityString,session,bandString,sf, con, ori,badTrialNameStr)
