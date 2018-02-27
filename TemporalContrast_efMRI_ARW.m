% fMRI stimulus code rapid event related design presenting flickering
% windowed grating at a range of different contrasts and temporal
% frequencies. We call in runPlaid_withISI that generates our sine grating stimuli. 
% All params are set within this code. Hold q to quit during ISI.

% MMH: 23/8/2016
% ARW: 27/2/2018 : Cleaned a little and made this a demo script for sharing
% with collaborators.

close all;
clear all;

startTime=tic;
DUMMYRUN=0; % Do we want to run without displaying anything (useful for testing)
WhichScreen = 0;
Screen('Preference', 'VisualDebuglevel', 1)% disables welcome and warning screens
Screen('Preference', 'verbosity', 1)

WaitMRIPulse = 1;
TR=3;
DummyTR = 4; % How many frames to 'burn' at the start?
ScanWaitTime = TR * DummyTR;

% dpy contains info about display that is relevant
dpy.res = [1920 1080]; % Note:VpixxProjector resolution 1920 1080
dpy.size = [.41 .235]; % Meters
dpy.distance = [.57]; % Distance from eyes to screen
dpy.frameRate = 60; % Hz

sfList=[1]'; %grating cpd: Spatial frequency

%stim.paramaters
stim.spatial.internalRotation = 0; % Does the grating rotate within the envelope?
stim.rotateMode = []; % rotation of mask grating (1= horizontal, 2= vertical, etc?)
stim.spatial.angle = [0 0]; % angle of gratings on screen
stim.temporal.duration=TR; % Flicker for 1TR (3 seconds) 
nRepeats=1; %how many times to repeat each condition.

%% get session details and create associated save file for results
R = num2str(input('Enter participant R number\n'));
run = num2str(input('Enter run number\n'));
myfile = sprintf('fRF_TCRF_R%s_run%s_%s',R,run,datestr(now,'ddmmyy_HHMM'));
[~,mydir] = uiputfile(myfile,'Choose file directory');
out_file = [mydir,myfile,'.mat'];
HideCursor

if (~DUMMYRUN)
    Screen('Preference', 'SkipSyncTests', 1);
    WhichScreen = 0; %
    Screen('Preference', 'VisualDebuglevel', 1)% disables welcome and warning screens
    PsychDebugWindowConfiguration(1);
    % Open a gray window
    dpy.win = Screen('OpenWindow', WhichScreen, 128);

    % Make sure the GLSL shading language is supported:
    AssertGLSL;
    
    % Retrieve video redraw interval for later control of our animation timing:
    dpy.ifi = Screen('GetFlipInterval', dpy.win);
    
else
    disp('Dummy run - not initializing the screen');
    dpy.win=-1;
end

Screen('BlendFunction', dpy.win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

%% Here we input our parameter levels
myContrasts = [1,4,8,16,64]/100; %contrasts in range 0-100
myTempFreq = [1,5,10,20]; %temporal frequncy in Hz 

% Optional: Pseudo ISI (mean = 6s) or a randomized ISI
PseudoMyISI = [3,3,6,3,3,6,3,6,6,9,6,3,9,6,9,6,3,9,3,3,9];
%RandMyISI = Shuffle(PsuedoMyISI); % You had better do this if you plan to
%run more than a single rep ....

%Use a random stimulus order to avoid adaptation effects.
%This is a faff but let's keep it.
StimCondt = [1:21];
FinalStimOrder = Shuffle(StimCondt);

% Now we set up an array that holds all combinations of our two stimulus
% paramaters.
% We can access stimulusParams for analysis to work out what order stim was
% presented.

stimulusParams = [];

% adds all combinations
for i = 1:length(myContrasts);
    for j = 1:length(myTempFreq);
        stimulusParams(end+1,1:2) = [myContrasts(i), myTempFreq(j)];
    end
end

% Here we add the null (0% contrast at 1Hz) 
% It will always be presented wherever '21' falls within our RandStimulusOrder variable
stimulusParams(end+1,1:2) = [0, 1]; %0% contrast at 1Hz = blank

%% Set up the stimulus arrays used by frf_runPlaid
myStimulusList = {};
for i=1:length(stimulusParams)
    % Set up stim first
    myStimulusList{i}.stim.spatial.frequency=[sfList(1), sfList(1)]; % Same spatial freq for each component
    
    % Phase is the phase shift in degrees (0-360 etc.)applied to the sine grating:
    myStimulusList{i}.stim.spatial.phase=[rand(1)*360 rand(1)*360]; %;
    myStimulusList{i}.stim.spatial.pOffset=rand(2,1)*360;
    myStimulusList{i}.stim.spatial.internalRotation = 0; %always 0
    myStimulusList{i}.stim.spatial.angle = [0 0]; %always 0
    myStimulusList{i}.stim.temporal.duration = 3;
    myStimulusList{i}.stim.temporal.frequency = [stimulusParams(i,2) stimulusParams(i,2) ]; %extract from stimulusParams 2nd column
    myStimulusList{i}.stim.rotateMode = []; %N/A
    myStimulusList{i}.stim.cont = [stimulusParams(i,1);0]; %Contrast level is extracted from stimulusParams 1st column

    %now set up dpy struct -- We need to calculate increment for flicker and contrast
    myStimulusList{i}.dpy = dpy; 
    myStimulusList{i}.dpy.phaseincrement = [myStimulusList{i}.stim.temporal.frequency] * 360 * dpy.ifi;
    myStimulusList{i}.dpy.contrastincrement = [myStimulusList{i}.stim.cont] * dpy.ifi; 
end

%% Fixation Task! 
% Let's present x at fixation that randomly changes to a + at random
% frame intervals. This means we need to change this based on the frame
% rate (usually 60Hz or 120Hz)

% totalFrames = (63s of stim) + (126s ofISIs) * 120Hz = 22680 frames
% totalFrames = (63s of stim) + (126s ofISIs) * 60Hz = 11340 frames

eventvect = ones(50,1);
v = sprand(22700,1,.001); %120Hz % A sparse list of events. This is long because we are assuming 120Hz.
%v = sprand(11340,1,.004); %60Hz
vv = full(v);
eventList = conv(vv,eventvect,'same');
fixationList = (eventList > 0);
fixationList = +fixationList;



%% Right, it's set up. Experiment begins now!
mri_startTime=clock;

%Scanner warm up
Screen('FillRect',dpy.win, 128);
DrawFormattedText(dpy.win,'Scanner is warming up...','center','center',[0, 0, 0]);
Screen('Flip',dpy.win);
WaitSecs(ScanWaitTime);

%prepare for fixation task
countFrame = 0;
fixResults = [];

for i=1:length(myStimulusList)
    %Retreive index of the stimulus number we want
    thisStimulusIndex = FinalStimOrder(i);
    %Use Psuedo ISI
    thisISI = PseudoMyISI(i);
    stim = myStimulusList{thisStimulusIndex}.stim;
    dpy = myStimulusList{thisStimulusIndex}.dpy;
    
    % Call in the main display loop to present the stimuli
    if (~DUMMYRUN)
        
        %function to present a plaid that we flicker
        d=runPlaid_withISIARW(dpy,stim,thisISI,fixationList,countFrame,fixResults);
    
        countFrame = d.countFrame;
        fixResults = d.fixResults';
        
        %keyboard check to quit (q to quit)
        [keyIsDown, secs, keyCode, deltaSecs] = KbCheck;
        if keyIsDown % This is a vector of key codes
            keyDown=find(keyCode);% Which things are down? We only take the first one.
            if (KbName(keyDown(1))=='q')
                sca
                error('abort experiment')
            end
        end
    end
end

% Stimuli presentation is complete, we close the window
mri_endTime = clock;
elapsedTime = etime(mri_endTime,mri_startTime);

%Important! Save this and ensure that you label it for each subject for
%analysis and decoding the stimulus order
save(out_file)

Screen('CloseAll');
sca

%% If you want to check your subject was behaving themselves
%fixCompare = horzcat(fixationList(1:length(fixResults)),fixResults);
%imagesc(fixCompare)