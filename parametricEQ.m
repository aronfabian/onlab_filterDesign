function [H_dB,f] = parametricEQ(gain_dB,fc,bw_oct,fs,COMP_FREQS);
%BOOST - Design a digital boost filter at given gain g, 
%        center frequency fc in Hz,
%        bandwidth bw in Hz (default = fs/10), and
%        bandwidth bw in octaves (default = fs/10), and
%        sampling rate fs in Hz (default = 1).


% % Példa
% g=7;
% [Ho,f]=parametricEQ(10,1000,1/3,44100);
% [H1,f]=parametricEQ(g,1000/2^(1/9),1/9,44100);
% [H3,f]=parametricEQ(g,1000*2^(1/9),1/9,44100);
% [H2,f]=parametricEQ(g,1000,1/6,44100);
% semilogx(f,[Ho H1+H2+H3])

if ~exist('COMP_FREQS','var') COMP_FREQS = 65536; end

if nargin<4, fs = 1; end
if nargin<3, bw = fs/10; end

% g = 10^(gain_dB/20);
% bw = bw_oct * fs/g;

% c = cot(pi*fc/fs) % bilinear transform constant
% cs = c^2;  csp1 = cs+1; Bc=(bw/fs)*c; gBc=g*Bc;
% nrm = 1/(csp1 + Bc); % 1/(a0 before normalization)
% b0 =  (csp1 + gBc)*nrm;
% b1 =  2*(1 - cs)*nrm;
% b2 =  (csp1 - gBc)*nrm;
% a0 =  1;
% a1 =  b1;
% a2 =  (csp1 - Bc)*nrm;
% v2

% frequency wraping
fc = fs*2*tan(fc*2*pi/(2*fs))/(2*pi);

omega = 2 * pi * fc;
Omega0 = atan(omega/2/fs)*2;
gamma = sinh( log(2)/2*bw_oct*Omega0/sin(Omega0) ) * sin(Omega0);
K = 10^(gain_dB/20);

b0 = ( 1 + gamma * sqrt(K) ) / ( 1 + gamma / sqrt(K) );
b1 = (-2 * cos(Omega0) ) / ( 1 + gamma / sqrt(K) );
b2 = ( 1 - gamma * sqrt(K) ) / ( 1 + gamma / sqrt(K) ); 
a0 =  1;
a1 = (-2 * cos(Omega0) ) / ( 1 + gamma / sqrt(K) );
a2 = ( 1 - gamma / sqrt(K) ) / ( 1 + gamma / sqrt(K) ); 


A = [a0 a1 a2];
B = [b0 b1 b2];



% [H,f]=freqz(B,A,65536,fs);

[H,f]=freqz(B,A,COMP_FREQS,fs);

if nargout==0
  
  %[H,f]=freqz(B,A,65536,fs); % /l/mll/myfreqz.m
%   figure(1); semilogx(f,20*log10(abs(H)))
  hold on; semilogx(f,20*log10(abs(H)))
%   dstr=sprintf('boost(%0.2f,%0.2f,%0.2f,%0.2f)',g,fc,bw,fs);
%   subplot(2,1,1); title(['Boost Frequency Response: ',...
%                       dstr],'fontsize',24);
else 
    H_dB = 20*log10(abs(H));   

% N = 2;
% Wo = fc/(fs/2);
% Q = sqrt(2^(0.25))/(2^(0.25)-1);
% BW = Wo/Q;
% [B,A] = designParamEQ(N,gain_dB,Wo,BW);
% section1 = dsp.IIRFilter('Numerator',B(:,1)','Denominator',[1,A(:,1)']);
% [H,f] = freqz(section1,COMP_FREQS,fs);
% H_dB = 20*log10(abs(H));
end