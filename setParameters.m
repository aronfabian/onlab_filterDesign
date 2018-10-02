function [] = setParameters(startFreq,endFreq,compLines,plotLines,aFilt,cFilt,noFilt,hpfButter,lpfButter,tolerance)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    START_FREQ = startFreq; % vizsgálati freki tartomány alja
    END_FREQ = endFreq; % vizsgálati freki tartomány teteje
    COMP_FLINES = compLines; % az átviteli fgv számításakor használt felbontás
    PLOT_FLINES = plotLines; % a megjelenítéshez használt felbontás

    A_FILT = aFilt; % compensation with A-filter
    C_FILT = cFilt; % compensation with C-filter
    NO_FILT = noFilt; % no compensation

    HPF_BUTTER = hpfButter; % use highpass butterworth filter 
    LPF_BUTTER = lpfButter; % use lowpass butterworth filter

    TOLERANCE = tolerance; % 1=Class 1, 2=Class 2, 3=Class 1 and Class 2

    PAUSE = 0; % 0=run without pausing, 1 = run with pausing

    fileNames = {
    %  'CutOf_iPhone_5S_MFRecorder_2017-04-27 17-08-07_SweepOnly_h_bSpectr'
    % 'CutOf_Samsung Galaxy Ace 20170427_173930_SweepOnly_h_bSpectr'
      'CutOf_Samsung Galaxy Alpha - 20170427_183001_SweepOnly_h_bSpectr'
    %'CutOf_Samsung Galaxy Mini 2 - 20170427_174623_SweepOnly_h_bSpectr'
    };

    %% ParamaetricEQEst

    EST_MODE = 1; % 1 = fc,gain,bw; 2 = random parameter select 

    % IF EST_MODE = 1
    OUTER_LOOP = 10;
    INNER_LOOP = 20;

    % IF EST_MODE = 2
    CYCLES = 2000;

    %% ToleranceOptim

    LIMIT = 0.1; % effective range limit
    MAX_CYCLE = 15; % maxumim number of cycles 
    EQUAL_TOLERANCE = 0.05; % stop if the difference between distance from lower and upper band is under this value

    %% ErrorCalc
    ERROR_SEL = 1; % sum of abs (ERROR_SEL = 1)
                  % sum of squares (ERROR_SEL = 2)
                  % sum of abs where H_mic is outside the tolerance band (ERROR_SEL = 3)
                  % sum of squares where H_mic is outside the tolerance band (ERROR_SEL = 4)


    %% SAVE variables

    save('config.mat')

end

