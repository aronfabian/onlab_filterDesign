

config();
f_interp = logspace(log10(START_FREQ),log10(END_FREQ),COMP_FLINES);
f_interp_plot = logspace(log10(START_FREQ),log10(END_FREQ),PLOT_FLINES);

if (A_FILT == 1)
    load('A_weight.mat');
    A_interp = interp1(A_weight(:,1), A_weight(:,2) , f_interp, 'pchip');
end
load('tolerance.mat');
    tol_interp = (interp1(tolerance(:,1), tolerance(:,2:3) , f_interp, 'pchip'))';
    tol_interp_plot = (interp1(tolerance(:,1), tolerance(:,2:3) , f_interp_plot, 'pchip'))';

load([fileNames{1}]);

X_interp = interp1(phone_bFreqs, phone_bX, f_interp, 'pchip');
X_interp_plot = interp1(phone_bFreqs, phone_bX, f_interp_plot, 'pchip');

% H_mic to dB scale
H_mic = 20*log10(X_interp);

% H_mic * A_weight
H_mic = H_mic - A_interp;

% find 1kHz
            trgt_f = f_interp(1:end-1)<1000 & f_interp(2:end)>1000;
            % H_mic offset with Magnitude[dB] at 1kHz 
            H_mic = H_mic - H_mic(trgt_f) ;


            % H_mic számolása megjelenítéshez
            H_mic_plot = interp1(f_interp, H_mic, f_interp_plot, 'pchip');

            H_trgt = zeros(1,length(f_interp));
            H_trgt_plot = zeros(1,length(f_interp_plot));
            figure()
            semilogx(f_interp_plot,H_mic_plot)
            title({fileNames{1}; '1'},'Interpreter', 'none')
            xlabel('Frequency (Hz)')
            ylabel('Magnitude (dB)')
            hold on
            semilogx(f_interp_plot,H_trgt_plot)
            semilogx(f_interp_plot,tol_interp_plot(1:2,:),'r--')
  H_mic_original = H_mic;
  fs=44100;
  gain = -11.22;
  fc = 57;
  bw = log2(1000/50);
  configFile='config.mat';
  hiba = 1.702569e+03;
  for n = 1:10
    [estGain,estFc,estBw,Ho,~] = parametricEQest(gain,fc,bw,fs,H_mic,50,1000,f_interp,H_trgt,tol_interp,configFile);
    [h,~] = errorCalc(H_mic+Ho,H_trgt,tol_interp,configFile);
    fprintf('A becslõ értékei: estFc: %0.0f Hz, estGain: %0.2f oct, estBw: %0.2f dB; eltérés a minimális hibától: %.2f  \n', estFc , estGain, estBw,abs(h-hiba))
  end             
 fs = 44100;
 fc = 50;
 gain = -10.36;
 bw_oct = 4.09;
 
  [Ho,~]=parametricEQ(gain,fc,bw_oct,fs,f_interp);
  H_mic = H_mic + Ho;
                    H_mic_plot = interp1(f_interp, H_mic, f_interp_plot, 'pchip');

                    % plot new transfer function
                    semilogx(f_interp_plot, H_mic_plot, 'DisplayName', string(1))
% fs = 44100;
%  fc = 3370;
%  gain = -8;
%  bw_oct = 1;
%  
%   [Ho,~]=parametricEQ(gain,fc,bw_oct,fs,f_interp);
%   H_mic = H_mic + Ho;
%                     H_mic_plot = interp1(f_interp, H_mic, f_interp_plot, 'pchip');
% 
%                     % plot new transfer function
%                     semilogx(f_interp_plot, H_mic_plot, 'DisplayName', string(1))
% 
%                     fs = 44100;
%  fc = 15700;
%  gain = -12;
%  bw_oct = 0.3;
%  
%   [Ho,~]=parametricEQ(gain,fc,bw_oct,fs,f_interp);
%   H_mic = H_mic + Ho;
%                     H_mic_plot = interp1(f_interp, H_mic, f_interp_plot, 'pchip');
% 
%                     % plot new transfer function
%                     semilogx(f_interp_plot, H_mic_plot, 'DisplayName', string(1))
%%
H_mic = H_mic_original;
% fc : 50-200
% gain: -5 - -12.5
% bw : 2 - 4.5
% cube = randi(2000,100,100,50);
% [X,Y,Z] = ndgrid(1:size(cube,1), 1:size(cube,2), 1:size(cube,3));
% pointsize = 30;
% scatter3(X(:), Y(:), Z(:), pointsize, cube(:));
errorCube = zeros(50,50,50);
fc = logspace(log10(50),log10(200),50);
gain = linspace(-5,-12.5,50);
bw = linspace(2,4.5,50);
fs = 44100;
configFile='config.mat';

for i = 1:50 %fc
    for j = 1:50 %gain 
        for k = 1:50 %bw
            [Ho,~] = parametricEQ(gain(j),fc(i),bw(k),fs,f_interp);
            [errorCube(i,j,k),~] = errorCalc(H_mic+Ho,H_trgt,tol_interp,configFile);
            
        end
    end
end
figure
[X,Y,Z] = ndgrid(1:size(errorCube,1), 1:size(errorCube,2), 1:size(errorCube,3));
pointsize = 30;
obj=scatter3(X(:), Y(:), Z(:), pointsize, errorCube(:));
xlabel('x (fc)')
ylabel('y (gain)')
zlabel('z (bw)')
c = colorbar;
c.Label.String = 'Hibaterület';
title('Hibaterületek a paramétertérben')

[minimum,minPlace]=min(errorCube(:));
[i,j,k] = ind2sub(size(errorCube),minPlace);
fprintf('Minimumhely f(x,y,z)=f(%d,%d,%d)=%d\n',i,j,k,minimum)
fprintf('Minimumhelyen fc=%.1f; gain=%.2f; bw=%.2f\n ',fc(i),gain(j),bw(k))
makedatatip(obj,minPlace)