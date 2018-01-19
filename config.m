clear all

START_FREQ = 50; % vizsgálati freki tartomány alja
END_FREQ = 15600; % vizsgálati freki tartomány teteje
COMP_FLINES = 1000; % az átviteli fgv számításakor használt felbontás
PLOT_FLINES = 1000; % a megjelenítéshez használt felbontás

A_FILT = 1; % compensation with A-filter
C_FILT = 1; % compensation with C-filter
NO_FILT = 1; % no compensation

HPF_BUTTER = 0; % use highpass butterworth filter 
LPF_BUTTER = 0; % use lowpass butterworth filter

PAUSE = 1; % 0=run without pausing, 1 = run with pausing

fileNames = {
 'CutOf_iPhone_5S_MFRecorder_2017-04-27 17-08-07_SweepOnly_h_bSpectr'
'CutOf_Samsung Galaxy Ace 20170427_173930_SweepOnly_h_bSpectr'
 'CutOf_Samsung Galaxy Alpha - 20170427_183001_SweepOnly_h_bSpectr'
'CutOf_Samsung Galaxy Mini 2 - 20170427_174623_SweepOnly_h_bSpectr'
};

%% ParamaetricEQEst

EST_MODE = 2; % 1 = fc,gain,bw; 2 = random parameter select 

% IF EST_MODE = 1
OUTER_LOOP = 10;
INNER_LOOP = 20;

% IF EST_MODE = 2
CYCLES = 500;

%% ToleranceOptim

LIMIT = 0.1; % effective range limit
MAX_CYCLE = 15; % maxumim number of cycles 
EQUAL_TOLERANCE = 0.1; % stop if the difference between distance from lower and upper band is under this value

%% ErrorCalc
ERROR_SEL = 1; % sum of abs (ERROR_SEL = 1)
              % sum of squares (ERROR_SEL = 2)
              % sum of abs where H_mic is outside the tolerance band (ERROR_SEL = 3)
              % sum of squares where H_mic is outside the tolerance band (ERROR_SEL = 4)


%% SAVE variables

save('config.mat')