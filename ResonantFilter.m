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


% Ennek a függvénynek a magja egy analóg rezonáns szûrõt tervez, a digitális
% megvalósított verziók ehhez képest a Superpowered-é a [0-végtelen[ frekvencia
% tartományt belenyomja a [0 fs/2[ tartományba. Ezt szimulálandó ennek a függvénynek
% a végén a [0-végtelen[ tartományt össze kell nyomni a [0 fs/2[ tartományba.
% Annak érdekében azonban, hogy a végsõ összenyomás következtében a
% frekvenciák majd a helyükön maradjanak, a függvény elején szét kell húzni
% a [0 fs/2[ tartományt [0-végtelen[ -be.
% A végén lévõ összehúzás az úgy történik, hogy az elõtorzított frekvencia
% pontokra kiszámított átvitel-értékeket az eredeti frekvenciákra
% értelmezem.

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