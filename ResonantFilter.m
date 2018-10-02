function [H_dB,f] = ResonantFilter(fc,type,fs,Q,COMP_FREQS)
% [H_dB,f] = ResonantFilter(fc,type,fs,COMP_FREQS)
% Designs a digital resonant filter with given 
%        cut-off frequency fc in Hz,
%        type (high-pass/low-pass) valid values: 'HPF' or 'LPF',
%        sampling rate fs in Hz,
%        Q factor (default is 0.707) and
% and returns the transfer function at COMP_FREQS frequencies  in Hz
% (default are 65k points in the range of [0..fs/2[).
% Q and/or COMP_FREQS can be omitted.


% Ennek a f�ggv�nynek a magja egy anal�g rezon�ns sz�r�t tervez, a digit�lis
% megval�s�tott verzi�k ehhez k�pest a Superpowered-� a [0-v�gtelen[ frekvencia
% tartom�nyt belenyomja a [0 fs/2[ tartom�nyba. Ezt szimul�land� ennek a f�ggv�nynek
% a v�g�n a [0-v�gtelen[ tartom�nyt �ssze kell nyomni a [0 fs/2[ tartom�nyba.
% Annak �rdek�ben azonban, hogy a v�gs� �sszenyom�s k�vetkezt�ben a
% frekvenci�k majd a hely�k�n maradjanak, a f�ggv�ny elej�n sz�t kell h�zni
% a [0 fs/2[ tartom�nyt [0-v�gtelen[ -be.
% A v�g�n l�v� �sszeh�z�s az �gy t�rt�nik, hogy az el�torz�tott frekvencia
% pontokra kisz�m�tott �tvitel-�rt�keket az eredeti frekvenci�kra
% �rtelmezem.

if length(Q)>1, COMP_FREQS = Q; clear Q; end
if ~exist('Q','var'), Q = sqrt(2)/2; end
if ~exist('COMP_FREQS','var'), COMP_FREQS = [0:65535]/65536*fs/2; end

[b,a] = butter(2,fc/(fs/2), lower(type));
    [H_dB,~] = freqz(b,a,COMP_FREQS,fs);
    H_dB = 20*log10(abs(H_dB));
    
% Frequency pre-warping
fc = fs*2*tan(fc*2*pi/(2*fs))/(2*pi);
f = fs*2*tan(COMP_FREQS*2*pi/(2*fs))/(2*pi);

wn=2*pi*fc;
w=2*pi*f;

j=sqrt(-1);

switch upper(type)
    case 'HPF'
%         H_dB = 20*log10(abs( -(w/wn).^2 ./ ( 1 -(w/wn).^2 +j*w/(Q*wn)) ));
    case 'LPF'
%         H_dB = 20*log10(abs( 1 ./ ( 1 -(w/wn).^2 +j*w/(Q*wn)) ));
end

semilogx(COMP_FREQS,H_dB)

% disp(H_dB(find(abs(f-fc)==min(abs(f-fc)))))
axis([20 20000 -60 20])