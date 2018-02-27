function [dataOut]=runPlaid_withISI(dpy,stim,ISI,fixationList,countFrame,fixResults)
% function dataOut=runPlaid_withISI(dpy,stim,ISI,fixationList,countFrame,fixResults)

% Generates a 2-component plaid. We use 1 component as as our flickering sine grating. 
% The conversion from contrast to amplitudes in the code is computed by a
% separate function flytv_computeAlphaAmps. Flicker params are set in main
% script. In this version we have a sine wave modulation.

% History:
% MMH 6/11/16: I have taken flyTV code and converted it for fMRI

%% Let's present a grating
% We have to play with time to look at temporal frequency. Contrast is
% easy.
% Compute the alpha and amplitudes that we will use for contrast
[amps,alpha]=computeAlphaAmps(stim.cont);

%How big do we want our grating window to be (in pixels)? Do the maths.
%478.738 is 10 degrees visual angle in pixels based on scanner display res
gratingtex1 = CreateProceduralSineGrating(dpy.win, dpy.res(1), dpy.res(2),[.5,.5,.5, 1], 478.738); % Bottom grating
gratingtex2 = CreateProceduralSineGrating(dpy.win, dpy.res(1), dpy.res(2),[.5 .5 .5 alpha], 478.738); % Top grating blend 50%

%Vertical blank
vbl = Screen('Flip', dpy.win);
vblendtime = vbl + stim.temporal.duration;
i=0;

% Update some grating animation parameters:
phase=stim.spatial.phase;
degToRad=pi/180;
pixelsPerMeter=dpy.res(1)/dpy.size(1);
metersPerDegree=dpy.distance*tan(degToRad);
pixPerDegree=pixelsPerMeter*metersPerDegree;
stim.spatial.frequencyCPerPixel=stim.spatial.frequency/pixPerDegree;

%Fixation Task Keys to look for
scanListVals = [KbName('1!'), KbName('2@'), KbName('3#'), KbName('4$')];  %4 buttons.
scanList = zeros(1,256); scanList(scanListVals) = 1;

lTime=vbl;
while (vbl < vblendtime)
    % Increment phase by the appropriate amount for this time period:
    tSecs = vbl - lTime;
    
    %Fixation Task
    countFrame = countFrame + 1;
    currentFrame = fixationList(countFrame);
    
    if currentFrame == 0
        fixshape = '+';
    else
        fixshape = 'x';
    end
   
    phase = phase + dpy.phaseincrement;
   
    %Remove the rounding for drifting stimuli
    %pMod = 180*(round(phase/180));
    pMod = [0 0];
    %sine modulate our contrast in sine wave using sin(tSecs*2pi)*temporal frequency
    contrast = amps.*sin(tSecs*2*pi.*stim.temporal.frequency);
    
    %Present Screen
    Screen('DrawTexture', dpy.win, [gratingtex1], [], [], [stim.spatial.angle(1)], [], [0], [], [], [stim.rotateMode], [pMod(1),stim.spatial.frequencyCPerPixel(1),contrast(1),0]');
    Screen('DrawTexture', dpy.win, [gratingtex2], [], [], [stim.spatial.angle(2)], [], [0], [], [], [stim.rotateMode], [pMod(2),stim.spatial.frequencyCPerPixel(2),contrast(2),0]');
    Screen('TextSize', dpy.win, 30);
    DrawFormattedText(dpy.win,fixshape,'center','center',[0, 0, 0]);
    
    % Show it at next retrace:
    vbl = Screen('Flip', dpy.win, vbl + 0.5 * dpy.ifi);
    
    %KbCheck
    [keyIsDown, ~, keyCode] = KbCheck([],scanList);
    if keyIsDown
        %30, 31, 32, 33 are 1234 for Apple (dummy runs)
        %49, 50, 51, 52 are 1234 for FORP (in scan)
        
        if find(keyCode,1) == 49 %button 1
            fixResults(countFrame) = 1;
        elseif find(keyCode,1) == 50 %button 2
            fixResults(countFrame) = 1;
        elseif find(keyCode,1) == 51 %button 3
            fixResults(countFrame) = 1;
        elseif find(keyCode,1) == 52 %button 4
            fixResults(countFrame) = 1;
        else
        fixResults(countFrame) = 0;   
        end
    end
end

% Now draw the ISI (fixation only) for the required time
while (vbl < vblendtime+ISI)
    countFrame = countFrame + 1;
    currentFrame = fixationList(countFrame);
    if currentFrame == 0
        fixshape = '+';
    else
        fixshape = 'x';
    end
   
    Screen('TextSize', dpy.win, 30);
    DrawFormattedText(dpy.win,fixshape,'center','center',[0, 0, 0]);
    vbl = Screen('Flip', dpy.win, vbl + 0.5 * dpy.ifi);
    
    %KbCheck + Keycode during ISI
    [keyIsDown, ~, keyCode] = KbCheck([],scanList);
    if keyIsDown
        %30, 31, 32, 33 are 1234 for mac (numpad)
        %49, 50, 51, 101 are 1234 for scanner
        if find(keyCode,1) == 49 %button 1
            fixResults(countFrame) = 1;
        elseif find(keyCode,1) == 50 %button 2
            fixResults(countFrame) = 1;
        elseif find(keyCode,1) == 51 %button 3
            fixResults(countFrame) = 1;
        elseif find(keyCode,1) == 52 %button 4
            fixResults(countFrame) = 1;
        else
        fixResults(countFrame) = 0;   
        end
    end
end

% Leave the screen in place for the next run and spit out the results
dataOut.countFrame = countFrame;
dataOut.fixResults = fixResults;


